#pragma once
#pragma once
#include <string>
#include <vector>
#include <memory>
#include <map>
#include "Util.h"
#include "PolymerCollection.h"

namespace CLPLib
{
	/// <summary>
	/// Convenient container for information typically associated with a tree node.
	/// </summary>
	struct NodeInfo
	{
		/// <summary>
		/// The newick string for the node's child tree.
		/// </summary>
		std::string ChildrenNewick;
		/// <summary>
		/// A name assigned to the node.
		/// </summary>
		std::string Name;
		/// <summary>
		/// The length of the branch from the node to its parent.
		/// </summary>
		double BranchLength;
	};

	/// <summary>
	/// A struct for constructing newick instruction strings for building a tree.
	/// </summary>
	struct TreeBuildInstruction
	{
		/// <summary>
		/// The recognized symbols for parsing a newick file:  '(', ')', ',', ';' 
		/// </summary>
		char Token;
		/// <summary>
		/// A string for the name assigned to a node.
		/// </summary>
		std::string Name;
		/// <summary>
		/// The name given to the ASCII 0 character (for convenience).
		/// </summary>
		static const char NULL_CHAR = (char)0; // ASCII data

											   /// <summary>
											   /// Returns the next tree-building instruction in a partial newick string.
											   /// </summary>
											   /// <param name="descriptor">A newick string.</param>
											   /// <param name="cursorPosition">The position of the cursor upon entrance. The cursor position will be 
											   /// update upon exit.</param>
											   /// <returns>The next instruction extracted from the input string.</returns>
		static TreeBuildInstruction NextCommand(std::string descriptor, int cursorPosition)
		{
			TreeBuildInstruction newCommand;

			int positionOfNextToken = Utils::IndexOfAny(descriptor, tokens, cursorPosition);
			if (positionOfNextToken < 0)
			{
				newCommand.Token = (char)NULL;
				newCommand.Name = descriptor;
				return newCommand;
			}

			newCommand.Name = descriptor.substr(cursorPosition, positionOfNextToken + 1);
			newCommand.Token = descriptor.substr(cursorPosition, 1)[0];
			cursorPosition = positionOfNextToken + 1;
			return newCommand;
		}

		static std::vector<char> tokens;
	};


	template<class T> class Tree;
	/// <summary>
	/// Represents a tree-structured graph recursively.
	/// </summary>
	/// <typeparam name="T">The type of the contents associated with each node of the tree.</typeparam>
	template<class T>
	class Tree : public std::enable_shared_from_this<Tree<T>>
	{
	public:
		/// <summary>
		/// The children trees of the this tree.
		/// </summary>
		std::vector<std::shared_ptr<Tree<T>>> Children;
		/// <summary>
		/// The parental tree of this tree.
		/// </summary>
		std::weak_ptr<Tree<T>> Parent;
		/// <summary>
		/// The name assigned to the root node of this tree.
		/// </summary>
		std::string Name;

		/// <summary>
		/// The length of the branch from the root of this tree to its parent.
		/// </summary>
		double BranchLength;

		double Age()
		{
			double age = BranchLength;
			if (Parent)
				age += Parent->Age();

			return age;
		}

		/// <summary>
		/// The Newick string describing this tree.
		/// </summary>
		std::string Descriptor;

		/// <summary>
		/// The largest number of branches from the root of this tree to any of its descendants.
		/// </summary>
		int MaxLevel;

		/// <summary>
		/// The number of branches between this tree and the global root.
		/// </summary>
		int Level;

		/// <summary>
		/// The object associated with the root node of this tree.
		/// </summary>
		shared_ptr<T> Contents;

		Tree(std::string newick, std::shared_ptr<Tree<T>> parent = shared_ptr<Tree<T>>())
		{
			Descriptor = newick;

			// remove unnecessary semicolon if it is present
			if (newick.back() == ';')
			{
				newick = newick.substr(0, newick.length() - 1);
			}

			if (parent)
			{
				Parent = parent;
				Level = parent->Level + 1;
			}

			NodeInfo ni = ParseNewick(newick);
			Name = ni.Name;
			BranchLength = ni.BranchLength;

			if (ni.ChildrenNewick != "")
			{
				std::vector<std::string> childNewicks = ExtractChildNewicks(ni.ChildrenNewick);
				for (auto &n : childNewicks)
				{
					Children.push_back(std::make_shared<Tree<T>>(n, this->shared_from_this()));
				}
			}

		}

		void AddChild(std::shared_ptr<Tree<T>> child)
		{
			Children.Add(child);
			child->Parent = this->shared_from_this();
			ComputeMaxLevel(this);
		}
	private:
		friend class TestHelper;

		static int findMatchingRightParenthesis(std::string text, int startPosition)
		{
			if (text[startPosition] != '(')
			{
				throw std::runtime_error("ArgumentException(The character at startPosition is not a left parenthesis.)");
			}
			int count = 0;
			for (unsigned int i = startPosition; i < text.length(); i++)
			{
				if (text[i] == '(')
				{
					count++;
				}
				else if (text[i] == ')')
				{
					count--;
				}

				if (count == 0)
				{
					return i;
				}
			}
			return -1;
		}

		static NodeInfo parseName(std::string label)
		{
			NodeInfo ni;
			std::vector<std::string> tempString = Utils::splitString(label, ':');
			ni.Name = tempString[0];
			if (tempString.size() > 1)
			{
				if (!Utils::tryParse<double>(tempString[1], ni.BranchLength))
				{
					ni.BranchLength = std::numeric_limits<double>::quiet_NaN();;
				}
			}
			else
			{
				ni.BranchLength = std::numeric_limits<double>::quiet_NaN();;
			}
			return ni;
		}

		static std::vector<std::string> ExtractChildNewicks(std::string newick)
		{
			std::vector<std::string> childNewicks;
			int level = 0;
			int start = 0;
			for (unsigned int i = 0; i < newick.length(); i++)
			{
				if (newick[i] == ',')
				{
					if (level == 0)
					{
						childNewicks.push_back(newick.substr(start, i - start));
						start = i + 1;
					}
				}
				else if (newick[i] == '(')
				{
					level++;
				}
				else if (newick[i] == ')')
				{
					level--;
				}

			}
			childNewicks.push_back(newick.substr(start, newick.length() - start));
			return childNewicks;
		}

		static NodeInfo ParseNewick(std::string newick)
		{
			std::string nameandlength;
			int lPos = -1;
			int rPos = -1;
			if (newick[0] == '(')
			{
				lPos = 0;
				rPos = findMatchingRightParenthesis(newick, 0);
				if (rPos < 1)
				{
					throw std::runtime_error("FormatException(Newick string is malformed: " + newick);
				}
				if (rPos < (int)newick.length() - 1)
				{
					nameandlength = newick.substr(rPos + 1);
				}
				else
				{
					nameandlength = "";
				}
			}
			else
			{
				nameandlength = newick;
			}

			NodeInfo ni = parseName(nameandlength);

			if (lPos > -1)
			{
				ni.ChildrenNewick = newick.substr(1, rPos - 1);
			}
			else
			{
				ni.ChildrenNewick = "NULL";
			}

			return ni;
		}

	public:
		/// <summary>
		/// Performs rescaling of the lengths of all the branches in this tree.
		/// </summary>
		/// <param name="length">The length of the sequence contained as content in the root node.</param>
		/// <param name="alpha">Alpha in the equation NewBranchLength = (gamma * BranchLength * length + alpha) / (length + beta)</param>
		/// <param name="beta">Beta in the equation NewBranchLength = (gamma * BranchLength * length + alpha) / (length + beta)</param>
		/// <param name="gamma">Gamma in the equation NewBranchLength = (gamma * BranchLength * length + alpha) / (length + beta)</param>
		void RescaleBranchLength(int length, double alpha, double beta, double gamma)
		{
			for (auto& child : Children)
			{
				child->RescaleBranchLength(length, alpha, beta, gamma);
			}
			BranchLength = (gamma * BranchLength * length + alpha) / (length + beta);
		}

		/// <summary>
		/// Computes the Newick descriptor for the string. This descriptor is computed by traversing the tree, 
		/// and does not refer to the tree's Descriptor itself.
		/// </summary>
		/// <param name="tree">The tree to serialize.</param>
		/// <returns>A Newick string containing the tree's description.</returns>
		static std::string GetDescriptor(Tree<T>& tree)
		{
			std::string desc;
			if (tree.Children.size() > 0)
			{
				desc += "(";
				for (auto& child : tree.Children)
				{
					desc += GetDescriptor(*child) + ",";
				}
				desc.pop_back();
				desc.push_back(')');
			}

			desc += tree.Name;
			if (std::isnan(tree.BranchLength))
			{
				desc += ":" + std::to_string(tree.BranchLength);
			}
			else
			{
				desc += ";";
			}

			return desc;
		}

		std::string GetDescriptor()
		{
			return GetDescriptor(this);
		}

		/// <summary>
		/// Computes the total number of nodes in a tree.
		/// </summary>
		/// <param name="tree">The tree whose nodes are to be counted.</param>
		/// <returns>The total number of nodes in the tree.</returns>
		static int CountNodes(Tree<T>& tree)
		{
			int returnVal = 1;
			for (auto& child : tree.Children)
			{
				returnVal += CountNodes(*child);
			}
			return returnVal;
		}

		static int CountOccupiedNodes(Tree<T>& tree)
		{
			int numOccupied = tree.Contents ? 1 : 0;
			for (auto& child : tree.Children)
			{
				numOccupied += CountOccupiedNodes(*child);
			}
			return numOccupied;
		}

		static bool AllNodesOccupied(Tree<T>& tree)
		{
			bool occupied = (bool)tree.Contents;
			for (auto &child : tree.Children)
			{
				occupied &= AllNodesOccupied(*child);
			}
			return occupied;
		}

		/// <summary>
		/// Counts the total number of leaves, or nodes without children, in a tree.
		/// </summary>
		/// <param name="tree">The tree whose leaves are to be counted.</param>
		/// <returns>The total number of leaves, or nodes without children.</returns>
		static int CountLeaves(Tree<T> tree)
		{
			int returnVal = tree.Children.size() == 0 ? 1 : 0;
			for (auto& child : tree.Children)
			{
				returnVal += CountLeaves(*child);
			}
			return returnVal;
		}

		/// <summary>
		/// Attaches this tree to another tree as its parent.
		/// </summary>
		/// <param name="p">The tree to be treated as the parent to this tree.</param>
		void SetParent(std::shared_ptr<Tree<T>> p, double distance = 0)
		{
			Parent = p;
			Level = (*p).Level;
			BranchLength = distance;
		}

		/// <summary>
		/// Computes the maximum number of branch-crossings required to reach all a tree's children.
		/// </summary>
		/// <param name="t">The tree to be examined.</param>
		static void ComputeMaxLevel(Tree<T>& t)
		{
			t.MaxLevel = t.Level;
			for (auto& child : t.Children)
			{
				ComputeMaxLevel(*child);
				if (child->MaxLevel > t.MaxLevel)
				{
					t.MaxLevel = child->MaxLevel;
				}
			}
		}

		/// <summary>
		/// The number of nodes with non-null contents.
		/// </summary>
		static int NBaubles;

		/// <summary>
		/// Associates the nodes of a tree with items from a list of objects.
		/// </summary>
		/// <param name="t">The tree to be decorated.</param>
		/// <param name="baubles">The list of objects used to decorate the tree. The keys are 
		/// the names associated with nodes on the tree, and the values are the individual objects.</param>
		static void Decorate(Tree<T>& t, std::map<std::string, T>& baubles)
		{
			for (auto &child : t.Children)
			{
				Decorate(*child, baubles);
			}

			for (auto &kvp : baubles)
			{
				if (kvp.first == t.Name)
				{
					t.Contents = std::make_shared<T>(kvp.second);
					NBaubles++;
				}
			}
			return;
		}

		static void Decorate(Tree<Polymer<T>> t, PolymerCollection<T> pc)
		{
			for (auto& child : t.Children)
			{
				Decorate(*child, pc);
			}

			for (auto& p : pc)
			{
				if (p.second->Name == t.Name)
				{
					t.Contents = p.second;
					NBaubles++;
				}
			}
			return;
		}

		/// <summary>
		/// Finds all of the internal nodes of a tree and returns them in a list.
		/// </summary>
		/// <param name="t">The tree to be examined.</param>
		/// <param name="intNodes">A list comprising the internal nodes of the tree being examined.</param>
		static void GetInternalNodes(shared_ptr<Tree<T>> t, vector<shared_ptr<Tree<T>>>& intNodes)
		{
			if (t->Children.size() == 0)
			{
				return;
			}

			for (auto &child : t->Children)
			{
				GetInternalNodes(child, intNodes);
			}

			if ((t->Parent).expired())
			{
				return;
			}
			intNodes.push_back(t);
		}

		static void GetTips(std::shared_ptr<Tree<T>> tree, std::vector<shared_ptr<Tree<T>>>& tips)
		{
			for (auto & child : tree->Children)
			{
				GetTips(child, tips);
			}
			if (tree->Children.size() == 0)
			{
				tips.push_back(tree);
			}
		}

		static shared_ptr<std::map<std::string, T>> CollectAllContents(shared_ptr<Tree<T>> tree, shared_ptr<std::map<std::string, T>> contents = shared_ptr<std::map<std::string, T>>())
		{
			for (auto& child : tree->Children)
			{
				CollectAllContents(child, contents);
			}

			//contesnts is a shared_ptr of map<string, T>
			//so we are inserting into themap. byt T
			contents->insert(std::make_pair(tree->Name, *tree->Contents));
			return contents;
		}
	};

	template <class T>
	int Tree<T>::NBaubles = 0;
}
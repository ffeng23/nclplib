#pragma once

#include <string>
#include <vector>
#include <tuple>
#include "Util.h"
#include "EnumString.h"
#include "PolymerPair.h"
#include "Translator.h"
#include "Tree.h"

namespace CLPLib
{
	/// <summary>
	/// Describes the three types of nucleotide mutation.
	/// </summary>
	enum class MutationType
	{
		/// <summary>
		/// Simple nucleotide point-substitution.
		/// </summary>
		Substitution,
		/// <summary>
		/// Insertion
		/// </summary>
		Insertion,
		/// <summary>
		/// Deletion
		/// </summary>
		Deletion
	};

	Begin_Enum_String(MutationType)
	{
		RegisterEnumerator(MutationType::Substitution, "Substitution");
		RegisterEnumerator(MutationType::Insertion, "Insertion");
		RegisterEnumerator(MutationType::Deletion, "Deletion");
	}
	End_Enum_String;


	enum SubstitutionType
	{
		Identity,
		Transition,
		TransversionToComplement,
		TransversionToNonComplement
	};

	/// <summary>
	/// Container for information about a nucleotide mutation.
	/// </summary>
	class Mutation
	{
	public:
		/// <summary>
		/// The type of the mutation.
		/// </summary>
		MutationType Type;
		/// <summary>
		/// The position of the mutation within the nucleotide sequence.
		/// </summary>
		int Position;
		/// <summary>
		/// Further details about the mutation. The form of these details depend on the mutation type.
		/// </summary>
		std::string Detail;

		static std::string ToString(Mutation& mutation, char delimiter = '\t')
		{
			std::string tmpstr = EnumString<MutationType>::From(mutation.Type);
			return tmpstr + delimiter + std::to_string(mutation.Position) + delimiter + mutation.Detail;
		}

		static std::tuple<std::string, Mutation> FromString(std::string mutationString, char delimiter = '\t')
		{
			Mutation mutation;

			std::vector<std::string> parsed = Utils::splitString(mutationString, delimiter);
			EnumString<MutationType>::To(mutation.Type, parsed[1]);
			mutation.Position = atoi(parsed[2].c_str());

			std::string sb = parsed[3];
			for (unsigned int i = 4; i < parsed.size(); i++)
			{
				sb += (delimiter + parsed[i]);
			}
			mutation.Detail = sb;
			return std::make_tuple(parsed[0], mutation);
		}

		std::string ToString()
		{
			return Mutation::ToString(*this);
		}
	};



	/// <summary>
	/// Encapsulates the result of a simple mutation analysis.
	/// </summary>
	class MutationList
	{
	public:
		/// <summary>
		/// A list of mutations.
		/// </summary>
		std::vector<Mutation> Mutations;
		/// <summary>
		/// The total bases examined.
		/// </summary>
		int TotalBases;
		/// <summary>
		/// The total number of nucleotide substitutions.
		/// </summary>
		int TotalSubstitutions;
		/// <summary>
		/// The total number of nucleotides deleted.
		/// </summary>
		int NucleotidesDeleted;
		/// <summary>
		/// The total number of nucleotides inserted.
		/// </summary>
		int NucleotidesInserted;
		/// <summary>
		/// The header for saving the MutationList to file.
		/// </summary>
		static std::string Header;

		MutationList()
		{
		}

		static std::string ToString(MutationList mutationList, char delimiter = '\t')
		{
			std::string sb = MutationList::Header + "\n";
			for (auto &m : mutationList.Mutations)
			{
				sb += (Mutation::ToString(m, delimiter) + "\n");
			}
			return sb;
		}

		static std::string ToString(std::map<std::string, MutationList> mutationDict, char delimiter = '\t')
		{
			std::string sb = (MutationList::Header + "\n");
			for (auto &kvp : mutationDict)
			{
				for (auto& m : kvp.second.Mutations)
				{
					sb += (kvp.first + "\t" + Mutation::ToString(m, delimiter) + "\n");
				}
			}
			return sb;
		}

		static MutationList FromString(std::string mutationListString)
		{
			MutationList ml;
			std::vector<std::string> lines = Utils::splitString(mutationListString, '\n');
			for (auto& line : lines)
			{
				if (line == MutationList::Header || line.length() == 0)
					continue;
				std::tuple<std::string, Mutation> tuple = Mutation::FromString(Utils::TrimString(line));
				ml.Mutations.push_back(std::get<1>(tuple));
			}
			return ml;
		}

		static std::map<std::string, MutationList> FromStringToDictionary(std::string mutationDictionaryString)
		{
			std::map<std::string, MutationList> mutationDictionary;
			std::vector<std::string> lines = Utils::splitString(mutationDictionaryString, '\n');
			for (auto &line : lines)
			{
				if (line == MutationList::Header || line.length() == 0)
					continue;

				std::tuple<std::string, Mutation> tuple = Mutation::FromString(Utils::TrimString(line));
				std::string tkey = std::get<0>(tuple);
				if (mutationDictionary.find(tkey) == mutationDictionary.end())
				{
					mutationDictionary.insert(std::make_pair(tkey, MutationList()));
				}
				mutationDictionary[tkey].Mutations.push_back(std::get<1>(tuple));
			}
			return mutationDictionary;
		}
	};

	/// <summary>
	/// Encapsulates the methods required to analyze the mutations between the nucleotide sequences in an aligned PolymerPair.
	/// </summary>
	class MutationFinder
	{
	public:
		/// <summary>
		/// Returns a list of mutations identified in the specified PolymerPair.
		/// </summary>
		/// <param name="pair">An aligned Polymer pair to be examined for mutations.</param>
		/// <returns>A list of mutations.</returns>
		static MutationList GetMutations(PolymerPair<byte, byte>& pair, int readingFrame = 0, int lowerLimit = 0, int upperLimitSetBack = 0)
		{
			// Polymer0 is treated as the parental sequence
			MutationList ml;
			bool deletion = false;
			bool insertion = false;
			std::string ins;

			int examinedLength = (int)pair.Index.size() - upperLimitSetBack;
			ml.TotalBases = examinedLength;

			for (int i = lowerLimit; i < examinedLength; i++)
			{
				if (pair.Index[i][0] > -1 && (pair.Polymer0->Seq)[pair.Index[i][0]] != 16 && pair.Index[i][1] > -1 && (pair.Polymer1->Seq)[pair.Index[i][1]] != 16)
				{
					if (!MixedNucleotideByte::AreConsistent((pair.Polymer0->Seq)[pair.Index[i][0]], (pair.Polymer1->Seq)[pair.Index[i][1]]))
					{
						int aaPos = (int)floor((pair.Index[i][0] - readingFrame) / 3.0);
						int codonStart = (int)(3 * floor((pair.Index[i][0] - readingFrame) / 3.0) + readingFrame);

						std::string aa0;
						if ((unsigned)codonStart + 2 < pair.Polymer0->Seq.size() && codonStart > -1)
						{
							auto tseq = pair.Polymer0->Subsequence(codonStart, codonStart + 2);
							auto subseq = TypeConverter::ToChars(tseq);
							//auto subseq = TypeConverter::ToChars(pair.Polymer0->Subsequence(codonStart, codonStart + 2));
							std::string codon0(subseq.begin(), subseq.end());
							aa0 = Translator::Translate(codon0);
						}
						else
						{
							aa0 = "X";
						}

						std::string aa1;
						if (codonStart + 2 < (int)pair.Polymer1->Seq.size() && codonStart > -1)
						{
							auto subseq = TypeConverter::ToChars(pair.Polymer1->Subsequence(codonStart, codonStart + 2));
							std::string codon1(subseq.begin(), subseq.end());
							aa1 = Translator::Translate(codon1);
						}
						else
						{
							aa1 = "X";
						}

						Mutation m;
						m.Type = MutationType::Substitution;
						m.Position = pair.Index[i][1];
						m.Detail = std::string(1, TypeConverter::ToChar((pair.Polymer1->Seq)[pair.Index[i][1]]))
							+ "\t" + std::string(1, TypeConverter::ToChar((pair.Polymer0->Seq)[pair.Index[i][0]])) +
							"\t" + std::to_string(aaPos) +
							+"\t" + aa1
							+ "\t" + aa0;

						ml.Mutations.push_back(m);
					}
				}
				else
				{
					if ((pair.Index[i][0] < 0 || (pair.Polymer0->Seq)[pair.Index[i][0]] == 16) && (pair.Index[i][1] > -1 && (pair.Polymer1->Seq)[pair.Index[i][1]] != 16))
					{
						// enter the deletion state
						deletion = true;
					}
					else if ((pair.Index[i][0] > -1 && (pair.Polymer0->Seq)[pair.Index[i][0]] != 16) && (pair.Index[i][1] < 0 || (pair.Polymer1->Seq)[pair.Index[i][1]] == 16))
					{
						// enter the insertion state
						insertion = true;
						ins = std::string(1, TypeConverter::ToChar((pair.Polymer0->Seq)[pair.Index[i][0]]));
					}
					for (int k = 1; i + k < (int)pair.Index.size(); k++)
					{
						if (deletion)
						{
							if (pair.Index[i + k][0] > -1 && (pair.Polymer0->Seq)[pair.Index[i + k][0]] != 16)
							{
								// exit the deletion state
								deletion = false;
								Mutation m;
								m.Type = MutationType::Deletion;
								m.Position = i;
								m.Detail = std::to_string(k);
								ml.Mutations.push_back(m);
								i += k - 1;
								break;
							}
						}
						else if (insertion)
						{
							if ((pair.Index[i + k][1] < 0 || (pair.Polymer1->Seq)[pair.Index[i + k][1]] == 16) && (pair.Index[i][0] > -1 && (pair.Polymer0->Seq)[pair.Index[i][0]] != 16))
							{
								ins += TypeConverter::ToChar((pair.Polymer0->Seq)[pair.Index[i + k][0]]);
							}
							else
							{
								insertion = false;
								Mutation m;
								m.Type = MutationType::Insertion;
								m.Position = i;
								m.Detail = ins.length() + "\t" + ins;
								ml.Mutations.push_back(m);
								i += k - 1;
								break;
							}
						}
					}
				}
			}
			Summarize(ml);
			return ml;
		}

		static void GetMutations(Tree<Polymer<byte>>& tree, std::map<std::string, MutationList>& mutations)
		{
			for (auto &child : tree.Children)
			{
				GetMutations(*child, mutations);
			}
			if (tree.Parent.expired())
			{
				return;
			}

			PolymerPair<byte, byte> pair(tree.Contents, tree.Parent.lock()->Contents);
			if ((int)tree.Contents->Seq.size() != tree.Parent.lock()->Contents->Seq.size())
			{
				throw std::runtime_error("Sequences must be the same length");
			}
			int length = (int)tree.Contents->Seq.size();
			pair.Index.resize(length);
			for (int i = 0; i < length; i++)
			{
				pair.Index[i] = std::vector<int>{ i, i };
			}

			MutationList muList = GetMutations(pair);
			mutations.insert(std::make_pair(tree.Contents->Name, muList));
		}

		/// <summary> 
		/// Produces a summary of the mutations passed in a list.
		/// </summary>
		/// <param name="ml">A mutation list. The list itself contains the summary upon return.</param>
		static void Summarize(MutationList& ml)
		{
			ml.TotalSubstitutions = 0;
			ml.NucleotidesDeleted = 0;
			ml.NucleotidesInserted = 0;
			for (auto &m : ml.Mutations)
			{
				switch (m.Type)
				{
				case MutationType::Substitution:
					ml.TotalSubstitutions++;
					break;
				case MutationType::Insertion:
				{
					auto strings = Utils::splitString(m.Detail, '\t');
					ml.NucleotidesInserted += std::stoi(strings[0]);
				}
				break;
				case MutationType::Deletion:
					ml.NucleotidesDeleted += std::stoi(m.Detail);
					break;
				}
			}
		}

		/// <summary> 
		/// Produces a summary of the mutations passed in a list.
		/// this function does not seem to be complete - Axin
		/// </summary>
		/// <param name="ml">A mutation list. The list itself contains the summary upon return.</param>
		static void CountTransitionsAndTransversion(MutationList& ml)
		{
			for (auto &m : ml.Mutations)
			{
				if (m.Type == MutationType::Substitution)
				{
					auto detail = Utils::splitString(m.Detail, '\t');
					char alpha = detail[0][0];
					char beta = detail[1][0];
				}
			}
		}
	};

	//note this re declares MutationFinder as template
	//which does not seem to be valid in c++
	/*
	template<class T, class U>
	static class MutationFinder
	{
	public static double MutationFrequency(MonomerPair<T, U>& mPair, PolymerPair<T, U>& pPair, int lowerLimit = 0, int upperLimitSetBack = 0)
	{
	int nTotal = pPair.Index.Length - upperLimitSetBack;
	double nMutations = 0;
	for (int i = lowerLimit; i < nTotal; i++)
	{
	if (pPair.Index[i][0] > -1 && pPair.Index[i][1] > -1)
	{
	nMutations += mPair.Mutation(pPair.Polymer0.Seq[pPair.Index[i][0]], pPair.Polymer1.Seq[pPair.Index[i][1]]);
	}
	else if (pPair.Index[i][0] < 0)
	{
	nMutations += 1 - mPair.SecondIsGap(pPair.Polymer1.Seq[pPair.Index[i][1]]);
	}
	else
	{
	nMutations += 1 - mPair.FirstIsGap(pPair.Polymer0.Seq[pPair.Index[i][0]]);
	}
	}
	return nMutations / (nTotal - lowerLimit);
	}
	};
	*/
}
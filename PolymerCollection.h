#pragma once

#include <memory>
#include <map>
#include <vector>
#include <functional>

#include "Reader.h"
#include "Polymer.h"
#include "Annotations.h"

using namespace std;

namespace CLPLib
{

	/// <summary>
	/// Provides a container for Polymers.
	/// </summary>
	/// <typeparam name="T">The type used to encode nucleotides.</typeparam>
	template<class T>
	class PolymerCollection
	{
		shared_ptr<MonomerEncoding<T>> monomerEncoding;
		shared_ptr<Reader<T>> reader;

	public:
		/// <summary>
		/// The name used to identify the collection.
		/// </summary>
		string Name;
		static int NextId;
	protected:
		//keep track of this cause slow performance when removing reads after merge contig.
		//int maxLength;

	public:
		/// <summary>
		/// The length of the longest sequence in the collection.
		/// </summary>
		int get_MaxLength()
		{
			getMaxLength(*this);
		}

	protected:

		//it is more convinient to use Id as the key, instead of string
		map<int, shared_ptr<Polymer<T>>> members;

	public:



		/// <summary>
		/// Default constructor.
		/// </summary>
		PolymerCollection(shared_ptr<MonomerEncoding<T>> encoding)
			:monomerEncoding(encoding), reader(new Reader<T>(encoding))
		{
			Name = "NULL";
		}

		/// <summary>
		/// Constructor.
		/// </summary>
		/// <param name="_members">A list containing the polynucleotides of the appropriate type to be stored in this instance.</param>
		PolymerCollection(shared_ptr<MonomerEncoding<T>> encoding, vector<shared_ptr<Polymer<T>>>& _members)
			:monomerEncoding(encoding), reader(new Reader<T>(encoding))
		{
			Name = "NULL";
			for (auto &m : _members)
			{
				m->Id = NextId++;
				members.insert(std::make_pair(m->Name, m));
			}
		}

		/// <summary>
		/// Constructor.
		/// </summary>
		/// <param name="_members">A list containing the polynucleotides of the appropriate type to be stored in this instance.</param>
		PolymerCollection(shared_ptr<MonomerEncoding<T>> encoding, string dataString, SequenceFileFormat format = SequenceFileFormat::FASTA)
			:monomerEncoding(encoding)
		{
			this.reader = std::make_shared<Reader<T>>(encoding);
			Name = "NULL";
			switch (format)
			{
			case SequenceFileFormat::FASTA:
				FromFASTA(dataString);
				break;
			case SequenceFileFormat::FASTQ:
				FromFASTQ(dataString);
				break;
			}
		}

		/// <summary>
		/// Determines whether the collection contains a polynucleotide with a given name.
		/// </summary>
		/// <param name="name">The name to seek.</param>
		/// <returns>True if the name is the name of one of the polynucleotides. False if not.</returns>
		bool Contains(int id)
		{
			return (members.find(id) != members.end());
		}

		bool contains(string name)
		{
			for (auto it = members.begin(); it != members.end(); it++)
			{
				if (it->second->Name == name)return true;
			}
			return false;
		}

		/// <summary>
		/// Renames a polynucleotide in the collection.
		/// </summary>
		/// <param name="name">The original name.</param>
		/// <param name="newName">The new name to replace the original one. Throws an argument exception if the 
		/// original name is contained in the collection, or if the new name is already contained in the collection.</param>
		void Rename(string name, string newName)
		{
			if (newName == name)
				return;

			if (this->Contains(name))
			{
				if (!this->Contains(newName))
				{
					auto p = members[name];
					p.Name = newName;
					members.insert(std::make_pair(newName, p));
					members.erase(name);
				}
				else
				{
					throw std::runtime_error("The new name duplicates an existing name.");
				}
			}
			else
			{
				throw std::runtime_error("No item named " + name + " was found.");
			}
		}

		/// <summary>
		/// Renames a polynucleotide in the collection.
		/// </summary>
		/// <param name="name">The original name.</param>
		/// <param name="newName">The new name to replace the original one. Throws an argument exception if the 
		/// original name is contained in the collection, or if the new name is already contained in the collection.</param>
		void Rename(map<int, string>& renameKey)
		{
			for (auto& pair : renameKey)
			{
				if (!this->Contains(pair.first))
				{
					throw std::runtime_error(std::string("ArgumentException: the name ") + std::to_string(pair.first) + " does not exist.");
				}
				if (this->Contains(pair.second))
				{
					throw std::runtime_error(std::string("ArgumentException: The name ") + pair.second + " already exists.");
				}
			}
			for (auto& pair : renameKey)
			{
				renameKey.insert(std::make_pair(pair.first, pair.second));
			}
		}

		/// <summary>
		/// Adds an ID/Polynuc vpair to the collection. Updates the MaxLength properties as necessary.
		/// </summary>
		/// <param name="p">The polynucleotide to add.</param>
		void AddNewPolymer(string name, string& sequence, std::shared_ptr<vector<byte>> qualityScores = nullptr)
		{
			if (this->Contains(name))
			{
				throw std::runtime_error("ArgumentException: A polymer with the same name is already present.");
			}

			members->insert(name, make_shared<Polymer<T>>(new Polymer<T>(monomerEncoding, name, sequence, qualityScores)));
			SetMemberId(name);
		}

	private:
		void SetMemberId(shared_ptr<Polymer<T>> polymer)
		{
			if (polymer->Id == 0)
			{
				polymer->Id = NextId++;
			}
		}

	public:
		//return a unique id to use
		int GetId()
		{
			return NextId++;
		}
		/// <summary>
		/// Adds an ID/Polynuc vpair to the collection. Updates the MaxLength properties as necessary.
		/// </summary>
		/// <param name="p">The polynucleotide to add.</param>
		void AddMember(shared_ptr<Polymer<T>> polymer)
		{
			// note that we are not checking for consistency of encoding
			SetMemberId(polymer);
			if (this->Contains(polymer->Id))
			{
				throw std::runtime_error("ArgumentException: A polymer with the same Id is already present.");
			}
			members.insert(make_pair(polymer->Id, polymer));
		}

		string GetMemberName(int id)
		{
			if (members.find(id) == members.end())return nullptr;
			return members[id]->Name;
		}

		std::vector<shared_ptr<Polymer<T>>> GetAllMembers()
		{
			std::vector<shared_ptr<Polymer<T>>> vect;
			for (auto&item : members)
			{
				vect.push_back(item.second);
			}
			return vect;
		}

		//given read name get id
		int GetMemberId(string name)
		{
			for (auto& item : members)
			{
				if (item.second->Name == name)
				{
					return item.second->Id;
				}
			}
		}

		/// <summary>
		/// Removes the named polynucleotide from the collection. Updates the MaxLength property as necessary.
		/// </summary>
		/// <param name="name">The name of the polynucleotide to remove. Throws an exception if the name is not present in the collection.</param>
		/// <returns>True if the named polynucleotide was present; false otherwise.</returns>
		bool RemoveMember(int id)
		{
			if (members.find(id) == members.end())
			{
				return false;
			}
			members.erase(id);
			return true;
		}

		/// <summary>
		/// Removes several polynucleotides at once.
		/// </summary>
		/// <param name="toRemove">A list of names of polynucleotides to remove.</param>
		/// <returns>True if and only if all named polynucleotides were present.</returns>
		bool RemoveMembers(vector<int>& toRemove)
		{
			bool allpresent = true;
			for (auto &name : toRemove)
			{
				if (!this->Contains(name))
				{
					allpresent = false;
				}
				else
				{
					RemoveIdFromIdNameMap(members[name]->Id);
					members.Remove(name);
				}
			}
			return allpresent;
		}

		/// <summary>
		/// Gets the number of polynucs stored in the instance. Read only.
		/// </summary>
		int Count()
		{
			return (int)members.size();
		}

		/// <summary>
		/// Implements the GetEnumerator method for the IEnumerable interface.
		/// </summary>
		/// <returns></returns>
		//IEnumerator GetEnumerator()
		//	{
		//		return new PolymerEnum<T>(members);
		//	}

	private:
		static int getMaxLength(PolymerCollection<T> &pc)
		{
			int length = 0;
			for (auto &p : pc)
			{
				if (p.second->Seq.size() > (unsigned)length)
				{
					length = (int)p.second->Seq.size();
				}
			}
			return length;
		}

		static bool getHasUniformLength(PolymerCollection<T> pc)
		{
			if (pc.Count < 1) return true;
			int len = pc.members.First().Value.Seq.size();
			for (auto&p : pc)
			{
				if (p.Seq.size() != len)
				{
					return false;
				}
			}
			return true;
		}

	public:

		/// <summary>
		/// Determines whether all of the polynucleotides in are the same length.
		/// </summary>
		/// <returns>True if and only if all the polynucleotides in this instance have the same length.</returns>
		bool HasUniformLength()
		{
			return getHasUniformLength(this);
		}

		/// <summary>
		/// Gets a list of the names of all the polynucleotides contained in the collection.
		/// </summary>
		/// <returns>A list of names of the collection members.</returns>
		vector<string> GetNames()
		{

			vector<string> names;
			names.reserve(members.size());
			for (auto& item : members)
			{
				names.push_back(item.second->Name);
			}
			return names;
		}

		vector<int> GetIds()
		{
			vector<int> names;
			names.reserve(members.size());
			for (auto& item : members)
			{
				names.push_back(item.first);
			}
			return names;
		}

		void FromFASTA(string fastaString)
		{
			//we will convert string mapping to id mapping from here.
			auto data = reader->ReadFASTA(fastaString);
			for (auto &item : data)
			{
				AddMember(item.second);
			}
		}

		void FromFASTQ(std::string fastqString)
		{

			auto data = reader->ReadFASTQ(fastqString);

			//we will convert string mapping to id mapping from here.
			//map<string, shared_ptr<Polymer<T>>> data = reader->ReadFASTQ(fastqString);
			for (auto &item : data)
			{
				AddMember(item.second);
			}
		}

		void FromPMF(string pmfString)
		{
			//map<string, shared_ptr<Polymer<T>>> data = reader->ReadNPMF(pmfString);
			auto data = reader->ReadNPMF(pmfString);
			for (auto &item : data)
			{
				AddMember(item.second);
			}
		}


	public:
		typename map<int, shared_ptr<Polymer<T>>>::iterator begin()
		{
			return members.begin();
		}

		typename map<int, shared_ptr<Polymer<T>>>::iterator end()
		{
			return members.end();
		}

		typename map<int, shared_ptr<Polymer<T>>>::const_iterator begin() const
		{
			return members.begin();
		}

		typename map<int, shared_ptr<Polymer<T>>>::const_iterator end() const
		{
			return members.end();
		}

		const shared_ptr<Polymer<T>> & operator[](int key) const {
			return members.at(key);
		}

		shared_ptr<Polymer<T>>& operator[](int key) {
			return members[key];
		}

		//this is to help debug only, since [] not usable when debugging
		shared_ptr<Polymer<T>>& GetPolymer(int key)
		{
			return members[key];
		}

		//return a vector of references to members of the collection.
		std::vector<std::reference_wrapper<Polymer<T>>> getMemberReferences()
		{
			std::vector<std::reference_wrapper<Polymer<T>>> retvec;
			std::vector<Polymer<T>> tst;
			transform(members.begin(), members.end(), back_inserter(retvec),
				[](std::pair<int, std::shared_ptr<Polymer<T>>> const & p)
			{
				return std::ref(*p.second);
			});

			return retvec;
		}

		std::vector<int> getMemberIds()
		{
			std::vector<int> retvec;
			for (auto&p : members)
			{
				retvec.push_back(p.second->Id);
			}
			return retvec;
		}


	};

	template<class T>
	int PolymerCollection<T>::NextId = 1;

}

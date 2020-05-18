#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "Util.h"
#include "Monomer.h"

namespace CLPLib
{
	/// <summary>
	/// States for a method that reads a fastq file.
	/// </summary>
	enum FastqReaderState
	{
		/// <summary>
		/// The initial state 
		/// </summary>
		Start,
		/// <summary>
		/// The state indicating that the sequence is being read.
		/// </summary>
		ReadingSeq,
		/// <summary>
		/// The state indicating that the quality scores are being read.
		/// </summary>
		ReadingQScore,
		/// <summary>
		/// The state indicating that nothing is happening.
		/// </summary>
		Waiting
	};


	//template<class X> class PolymerCollection;

	template <class T>
	class Reader
	{
	protected:
		shared_ptr<MonomerEncoding<T>> encoding;

	public:
		string Message;

		Reader(shared_ptr<MonomerEncoding<T>> encoding)
			: encoding(encoding)
		{
		}

		map<string, shared_ptr<Polymer<T>>> ReadFASTA(string fastaString)
		{
			map<string, shared_ptr<Polymer<T>>> members;
			std::string seq;
			string name = "";

			ifstream reader(fastaString);
			if (!reader.is_open())
			{
				throw std::runtime_error("Error opening file: " + fastaString);
			}

			string line;
			while (getline(reader, line))
			{
				if (line == "")
				{
					continue; // skip all blank lines
				}
				if (line[0] == '>') // this line contains the polymer name
				{
					if (seq.length() > 0) // process the preceding polymer sequence
					{
						if (members.find(name) == members.end())
						{
							std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
							members.insert(std::make_pair(name, std::shared_ptr<Polymer<T>>(new Polymer<T>(encoding, name, seq))));
						}
					}
					seq.clear();

					//trim string
					//name = line.Substring(1).Trim();
					auto it = line.begin() + 1;
					while (it != line.end() && isspace(*it))it++;
					auto rit = line.rbegin();
					while (rit.base() != it && isspace(*rit))rit++;
					name = std::string(it, rit.base());
				}
				else
				{
					for (auto & c : line)
					{
						if (!isspace(c))
						{
							seq.push_back(c);
						}
					}
				}
			}
			if (seq.length() > 0)
			{
				if (members.find(name) == members.end())
				{
					std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
					members.insert(std::make_pair(name, std::shared_ptr<Polymer<T>>(new Polymer<T>(encoding, name, seq))));
				}
			}

			return members;
		}

		map<string, shared_ptr<Polymer<T>>> ReadFASTQ(string fastqString)
		{

			map<string, shared_ptr<Polymer<T>>> members;

			int offset = -33;

			FastqReaderState state = FastqReaderState::Start;
			std::vector<byte> qual;
			std::string seq;

			string name;
			unsigned id = 0;
			ifstream reader(fastqString);
			if (!reader.is_open())
			{
				throw std::runtime_error("Error opening file: " + fastqString);
			}
			string line;
			while (getline(reader, line))
			{
				Utils::TrimString(line);
				if (line != "")
				{
					switch (state)
					{
					case FastqReaderState::Start:
						if (line[0] != '@')
						{
							continue;
						}
						else
						{
							name = line.substr(1);
						}
						state = FastqReaderState::ReadingSeq;
						break;
					case FastqReaderState::ReadingSeq:
						seq = line;
						state = FastqReaderState::Waiting;
						break;
					case FastqReaderState::Waiting:
						if (line[0] != '+')
						{
							throw std::runtime_error("Unexpected input: " + line);
						}
						state = FastqReaderState::ReadingQScore;
						break;
					case FastqReaderState::ReadingQScore:
						qual.resize(line.length());
						for (unsigned int i = 0; i < line.length(); i++)
						{
							auto q = (byte)(line[i] + offset);
							qual[i] = q;
						}
						state = FastqReaderState::Waiting;
						if (members.find(name) == members.end())
						{
							if (seq.length() == qual.size())
							{
								members.insert(std::make_pair(name, std::shared_ptr<Polymer<T>>(new Polymer<T>(encoding, name, seq, vector<byte>(qual)))));
							}
							else
							{
								Message += name + "\n";
							}
						}
						id++;
						state = FastqReaderState::Start;
						break;
					}
				}
			}

			return members;
		}

		map<string, shared_ptr<Polymer<T>>> ReadNPMF(string npmfString)
		{
			map<string, shared_ptr<Polymer<T>>> members;

			ifstream reader(npmfString);
			if (!reader.is_open())
			{
				throw std::runtime_error("Error opening file: " + npmfString);
			}
			string line;
			while (getline(reader, line))
			{
				Utils::TrimString(line);
				while (line != "")
				{
					auto fields = Utils::splitString(line, '\t');
					if (fields.size() <= 1)continue;
					string name = fields[0];
					int seqlen = (int)fields.size() - 1;

					vector<vector<double>> sequence(seqlen, vector<double>(5));
					for (int iNuc = 0; iNuc < 5; iNuc++)
					{
						if (!getline(reader, line))break;
						fields = Utils::splitString(line, '\t');
						if (fields.size() - 1 != seqlen)
						{
							throw std::runtime_error("nmpf length");
						}
						for (int iPos = 0; iPos < sequence.size(); iPos++)
						{
							if (!Utils::tryParse<double>(fields[iPos + 1], sequence[iPos][iNuc]))
							{
								throw std::runtime_error("Error parsing nmpf value");
							}
						}
					}

					//skip the error line.
					if (!getline(reader, line))break;
					auto seqsp = std::make_shared<vector<vector<double>>>(sequence);
					members.insert(std::make_pair(name, std::shared_ptr<Polymer<T>>(new Polymer<T>(encoding, name, seqsp))));

					if (!getline(reader, line))break;
				}
			}

			return members;
		}
	};


}
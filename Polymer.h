#pragma once

#include <vector>
#include <string>
#include <memory>
#include <iterator>
#include "Complementer.h"
#include "Monomer.h"
#include "Writer.h"

using namespace std;
#define uid(x)  ((x) > 0 ? (x) : (-(x)))

namespace CLPLib
{

	template <class T>
	class Polymer
	{
	private:

	
	public:

		vector<T> Seq;

		//PureNucleotide encodeed.
		vector<byte> PureSeq;

		vector<byte> Qual;

		string UID;

		string Name;

		//fields not in oringial LPLIb - consider to remove
		int Id;
		//std::string charSeq;


		Complementer<T> complementer;

		Polymer(shared_ptr<MonomerEncoding<T>> encoding, string name = "", vector<T> sequence = vector<T>(), vector<byte> qual = vector<byte>())
			:Name(name), Seq(sequence), Qual(qual), Id(0)
		{
			//complementer.SetEncoding(encoding);
			if (Qual.size() > 0 && Qual.size() != Seq.size())
			{
				fprintf(stderr, "The quality and sequence arrarys are not the same length.\n");
				throw std::runtime_error("The quality and sequence arrarys are not the same length.");
			}

			//check if only acgtn
			if (checkPureNucleotide(Seq) == true)
			{
				PureSeq = to_pure(Seq);
			}
		}


		Polymer(shared_ptr<MonomerEncoding<T>> encoding, string& name, string& sequence, vector<byte> qual = vector<byte>())
			:Name(name), Qual(qual), Id(0)
		{

			//charSeq = sequence;
			if (sequence.length() > 0 && Qual.size() != sequence.length())
			{
				fprintf(stderr, "The quality and sequence arrarys are not the same length.\n");
				throw std::runtime_error("ArgumentException(The quality and sequence arrarys are not the same length.");
			}

			complementer.SetEncoding(encoding);

			if (sequence.length() != 0)
			{
				this->Seq = encoding->FromString(sequence);
			}
		}

		//check if the sequence is pureNucletode, i.e. acgtn
		bool isPureNucleotide()
		{
			return PureSeq.size() > 0;
		}

		//trim lowqual ends???
		//for now trim if the bases is N
		bool Trim()
		{
			//front
			int left_trim = 0;
			for (int i = 0; i < (int)this->Qual.size(); i++)
			{
				//qual <=2 very low quality, usually has base call of 'N'
				if ((Qual)[i] > 2)break;
				left_trim++;
			}
			if (left_trim > 0)
			{
				this->TrimLeft(left_trim);
			}
			//back.
			int right_trim = 0;
			for (int i = (int)this->Qual.size() - 1; i >= 0; i--)
			{
				if ((Qual)[i] > 2)break;
				right_trim++;
			}
			if (right_trim > 0)
			{
				this->TrimRight(right_trim);
			}
			return left_trim + right_trim > 0;
		}

		/// <summary>
		/// Produces a new nucleotide sequence by removing nucleotides from the left-hand (5') end of a given nucleotide sequence.
		/// </summary>
		/// <param name="seq">The nucleotide sequence to be trimmed.</param>
		/// <param name="numberToCut">The number of nucleotides to remove from the 5' end.</param>
		/// <returns>A new nucleotide sequence in the same encoding as the argument.</returns>
		static void TrimLeft(vector<T>& seq, int numberToCut)
		{
			if (numberToCut < 0 || numberToCut >(int)seq.size())
			{
				throw std::runtime_error("ArgumentOutOfRangeException");
			}
			seq.erase(seq.begin(), seq.begin() + numberToCut);
		}

		/// <summary>
		/// Produces a new nucleotide sequence by removing nucleotides from the right-hand (3') end of a given nucleotide sequence.
		/// </summary>
		/// <param name="seq">The nucleotide sequence to be trimmed.</param>
		/// <param name="numberToCut">The number of nucleotides to remove from the 3' end.</param>
		/// <returns>A new nucleotide sequence in the same encoding as the argument.</returns>
		static std::shared_ptr<vector<T>> TrimRight(vector<T>& seq, int numberToCut)
		{
			if (numberToCut > (int)seq.size() || numberToCut < 0)
			{
				throw std::runtime_error("ArgumentOutOfRangeException");
			}

			int len = (int)seq.size() - numberToCut;
			return std::make_shared<vector<T>>(seq.begin(), seq.begin() + len);
		}

		/// <summary>
		/// Produces a new nucleotide sequence as a contiguous subsequence of a give nucleotide sequence. 
		/// </summary>
		/// <param name="seq">The nucleotide sequence to be trimmed.</param>
		/// <param name="startIndex">The index of the first nucleotide to be included in the subsequence.</param>
		/// <param name="endIndex">The index of tha last nucleotide to be included in the subsequence.</param>
		/// <returns>A new nucleotide sequence in the same encoding as the argument.</returns>
		static vector<T> Subsequence(vector<T>& seq, int startIndex, int endIndex)
		{
			if (startIndex < 0 || endIndex >= (int)seq.size() || endIndex < startIndex - 1)
			{
				throw std::runtime_error("Indices are inappropriate for subsequence isolation: " + std::to_string(startIndex) + "," + std::to_string(endIndex) + "," + std::to_string(seq.size()));
			}
			int length = endIndex - startIndex + 1;
			if (length == 0)
			{
				return vector<T>(0);
			}
			return vector<T>(seq.begin() + startIndex, seq.begin() + endIndex + 1);
		}

		/// <summary>
		/// Trims this nucleotide sequence by removing nucleotides from the left-hand (5') end.
		/// </summary>
		/// <param name="numberToCut">The number of nucleotides to be removed from the 5' end.</param>
		void TrimLeft(int numberToCut)
		{
			this->Seq = TrimLeft(*this->Seq, numberToCut);
			/*
			if (!this->charSeq.empty())
			{
				this->charSeq = charSeq.substr(numberToCut);
			}
			*/
			//qual may or may not stored.
			if (this->Qual.size() > 0)
			{
				this->Qual = Polymer<byte>::TrimLeft(this->Qual, numberToCut);
			}
		}

		/// <summary>
		/// Trims this nucleotide sequence by removing nucleotides from the right-hand (3') end.
		/// </summary>
		/// <param name="numberToCut">The number of nucleotides to be removed from the 3' end.</param>
		void TrimRight(int numberToCut)
		{
			this->Seq = TrimRight(*this->Seq, numberToCut);
                        
                        //why this is here.
			if (this->Qual.size() > 0)
				this->Qual = Polymer<byte>::TrimRight(this->Qual, numberToCut);
		}

		/// <summary>
		/// Trims this nucleotide sequence by removing nucleotides from both ends.
		/// </summary>
		/// <param name="startIndex">The index of the first nucleotide to be retained.</param>
		/// <param name="endIndex">The index of the last nucleotide to be retained.</param>
		vector<T> Subsequence(int startIndex, int endIndex)
		{
			return Subsequence(this->Seq, startIndex, endIndex);
		}

		/// <summary>
		/// Produces the complement of a Polynucleotide.
		/// </summary>
		/// <param name="p">The polynucleotide to be complemented.</param>
		/// <returns>The polynucleotide whose sequence is the complement of the input polynucleotide. The name has ".C" appended.</returns>
		std::shared_ptr<Polymer<T>> Complement()
		{
			return complementer.Complement(*this);
		}

		static shared_ptr<Polymer<vector<double>>> Normalize(Polymer<vector<double>>& p)
		{
			vector<vector<double>> normalizedSequence = vector<vector<double>>(p.Seq.size(), vector<double>(5));
			for (int i = 0; i < (int)normalizedSequence.size(); i++)
			{
				double norm = 0;
				for (int n = 0; n < 5; n++)
				{
					norm += (p.Seq)[i][n];
				}

				if (norm == 0)
					continue;
				for (int n = 0; n < 5; n++)
				{
					normalizedSequence[i][n] = p.Seq[i][n] / norm;;
				}
			}
			return std::make_shared<Polymer<vector<double>>>(std::make_shared<NucleotidePMF>(), p.Name + "-normalized", normalizedSequence);
		}

		static double GetExpectedError(Polymer<vector<double>>& p)
		{
			double expectedError = (double)p.Seq.size();
			for (auto &pmf : p.Seq)
			{
				expectedError -= (NucleotidePMF::IDFlux(pmf, pmf) + pmf[4] * pmf[4]);
			}
			return expectedError;
		}

		int GetLength()
		{
			if (!Seq)return 0;
			return (int)Seq.size();
		}

		bool checkPureNucleotide(std::vector<vector<double>>& _seq)
		{
			return false;
		}

		bool checkPureNucleotide(std::vector<byte>& _seq)
		{
			for (byte i : _seq)
			{
				switch (i)
				{
				case 1:
				case 2:
				case 4:
				case 8:
				case 15:
					break;
				default:
					return false;
				}
			}
			return true;
		}

		std::vector<byte> to_pure(std::vector<byte>& _seq)
		{
			vector<byte> newseq;
			newseq.reserve(this->Seq.size());
			std::transform(Seq.begin(), Seq.end(), std::back_inserter(newseq),
				[](byte a) -> byte {
				switch (a)
				{
				case 1:return 0;
				case 2:return 1;
				case 4:return 2;
				case 8:return 3;
				default:return 4;
				}
			});
			return newseq;
		}

		std::vector<byte> to_pure(std::vector<vector<double>>& _seq)
		{
			throw "not implemented exception.\n";
		}

	};


	//template<class T>
	//Complementer<T> Polymer<T>::complementer;

}


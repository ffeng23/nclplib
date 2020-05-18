#pragma once

#include <string>
#include <memory>
#include <vector>

#include "TypeConverter.h"

namespace CLPLib
{

	//forward declaration
	template <class T> class MonomerEncoding;
	template <class T> class Polymer;
	template <class T> class PolymerCollection;
	template <class T, class U> class PolymerPair;
	class Trefoil;


	template <class T>
	class Writer
	{
	public:
		static int QScoreOffset;

		static std::string WriteFasta(std::shared_ptr<MonomerEncoding<T>> encoding, PolymerCollection<T>& pc);

		static std::string WriteFasta(std::shared_ptr<MonomerEncoding<T>> encoding, Polymer<T>& p, bool keepName =true);

		static std::string WriteFasta(std::shared_ptr<MonomerEncoding<T>>, PolymerPair<T, T>& pp);

		/*static std::string WriteFasta(PolymerPair<std::vector<double>, std::vector<double>>& pp);*/

		static std::string WriteFasta_xx(PolymerPair<T, T>& pp);


		static std::string WriteFasta(PolymerPair<std::vector<double>, std::vector<double>>& pp, int index);

		static std::string WriteFasta(Trefoil& trefoil);

		/// <summary>
		/// Reads a polynucleotide collection from a FASTQ file, returning it as a Dictionary of polynucs keyed 
		/// by auto-generated IDs.
		/// </summary>
		/// <param name="filename">The name of the FASTA file.</param>
		/// <returns>A Dictionary of polynucs keyed by auto-generated IDs.</returns>
		static std::string WriteFastq(std::shared_ptr<MonomerEncoding<T>> encoding, PolymerCollection<T>& pc);

		static std::string WriteNPMF(Polymer<std::vector<double>>& p);

		static std::string WriteNPMF(PolymerCollection<std::vector<double>>& pc);

		static std::string WriteNPMF(PolymerPair<std::vector<double>, std::vector<double>>& pair);


		static int get_number();

	};

	template<class T>
	int Writer<T>::QScoreOffset = -33;


	template <class T>
	std::string Writer<T>::WriteFasta_xx(PolymerPair<T, T>& pp)
	{
		std::vector<char> seq1 = TypeConverter::ToChars(pp.Polymer0->Seq);
		std::vector<char> seq2 = TypeConverter::ToChars(pp.Polymer1->Seq);

		int indexLength = (int)pp.Index.size();
		std::string s;

		s += ">" + pp.Polymer0->Name + "\n";
		for (int k = 0; k < indexLength; k++)
		{
			if (pp.Index[k][0] >= 0)
				s += seq1[pp.Index[k][0]];
			else
				s += '-';
		}
		s += "\n";

		s += ">" + pp.Polymer1->Name + "\n";;
		for (int k = 0; k < indexLength; k++)
		{
			if (pp.Index[k][1] >= 0)
				s += seq2[pp.Index[k][1]];
			else
				s += '-';
		}
		s += "\n";
		return s;
	}

}

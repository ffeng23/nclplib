
#include "stdafx.h"
#include "Writer.h"

#include "PolymerPair.h"
#include "Trefoil.h"
#include "PolymerCollection.h"

namespace CLPLib
{

	//testing
	template <class T>
	int Writer<T>::get_number()
	{
		return 42;
	}


	template <class T>
	std::string Writer<T>::WriteFasta(std::shared_ptr<MonomerEncoding<T>> encoding, Polymer<T>& p, bool keepName)
	{
		std::string s;
		if (keepName)
		{
			s += ">" + p.Name;
		}

		vector<T>& seqref = (p.Seq);
		for (auto &p1 : seqref)
		{
			s += encoding->ToChar(p1);
		}
		return s;
	}

	template <class T>
	std::string Writer<T>::WriteFasta(std::shared_ptr<MonomerEncoding<T>> encoding, PolymerPair<T, T>& pp)
	{

		int indexLength = (int)pp.Index.size();
		std::string s;
		s += ">" + pp.Polymer0->Name;
		for (int k = 0; k < indexLength; k++)
		{
			if (pp.Index[k][0] >= 0)
			{
				s += encoding->ToChar(((pp.Polymer0->Seq))[pp.Index[k][0]]);
			}
			else
				s += '-';
		}
		s += "\n";

		s += ">" + pp.Polymer1->Name + "\n";
		for (int k = 0; k < indexLength; k++)
		{
			if (pp.Index[k][1] >= 0)
				s += encoding->ToChar((pp.Polymer1->Seq)[pp.Index[k][1]]);
			else
				s += '-';
		}
		s += "\n";
		return s;
	}

	template <class T>
	std::string Writer<T>::WriteFasta(Trefoil& trefoil)
	{
		vector<vector<char>> seqs(4);
		seqs[0] = TypeConverter::ToChars(trefoil.Polymers[0]->Seq);
		seqs[1] = TypeConverter::ToChars(trefoil.Polymers[1]->Seq);
		if (!trefoil.MRCANormalized)
		{
			if (!trefoil.MRCALikelihoodSequence)
			{
				seqs[2] = vector<char>(trefoil.Index.size());
			}
			else
			{
				trefoil.MRCANormalized = Polymer<vector<double>>::Normalize(*trefoil.MRCALikelihoodSequence);
				seqs[2] = TypeConverter::ToChars(trefoil.MRCANormalized->Seq);
			}
		}
		else
		{
			seqs[2] = TypeConverter::ToChars(trefoil.MRCANormalized->Seq);
		}
		seqs[3] = TypeConverter::ToChars(trefoil.RootNormalized->Seq);

		string names[4];
		names[0] = trefoil.Polymers[0]->Name;
		names[1] = trefoil.Polymers[1]->Name;
		if (trefoil.MRCALikelihoodSequence)
		{
			names[2] = trefoil.MRCALikelihoodSequence->Name;
		}
		else
		{
			names[2] = "SKIP";
		}
		names[3] = trefoil.RootLikelihoodSequence->Name;

		std::string s;

		for (int i = 0; i < 4; i++)
		{
			if (names[i] == "SKIP")
				continue;

			s += (">" + names[i]);
			for (int k = 0; k < trefoil.AlignmentLength(); k++)
			{
				if (trefoil.Index[k][i] >= 0)
					s += (seqs[i][trefoil.Index[k][i]]);
				else
					s += ('-');
			}
			s += "\n";
		}
		return s;
	}


	//template <class T>
	//std::string Writer<T>::WriteFasta_xx(PolymerPair<T, T>& pp)
	//{
	//	vector<char> seq1 = TypeConverter::ToChars(*pp.Polymer0->Seq);
	//	vector<char> seq2 = TypeConverter::ToChars(*pp.Polymer1->Seq);

	//	int indexLength = (int)pp.Index.size();
	//	std::string s;

	//	s += ">" + pp.Polymer0->Name + "\n";
	//	for (int k = 0; k < indexLength; k++)
	//	{
	//		if (pp.Index[k][0] >= 0)
	//			s += seq1[pp.Index[k][0]];
	//		else
	//			s += '-';
	//	}
	//	s += "\n";

	//	s += ">" + pp.Polymer1->Name + "\n";;
	//	for (int k = 0; k < indexLength; k++)
	//	{
	//		if (pp.Index[k][1] >= 0)
	//			s += seq2[pp.Index[k][1]];
	//		else
	//			s += '-';
	//	}
	//	s += "\n";
	//	return s;
	//}

	/*
	template <class T>
	std::string Writer<T>::WriteFasta(PolymerPair<std::vector<double>, std::vector<double>>& pp)
	{
		vector<char> seq1 = TypeConverter::ToChars(*pp.Polymer0->Seq);
		vector<char> seq2 = TypeConverter::ToChars(*pp.Polymer1->Seq);

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
	*/

	template <class T>
	std::string Writer<T>::WriteFasta(PolymerPair<std::vector<double>, std::vector<double>>& pp, int index)
	{
		vector<char> seq1 = TypeConverter::ToChars(pp.Polymer0->Seq);
		vector<char> seq2 = TypeConverter::ToChars(pp.Polymer1->Seq);

		int indexLength = (int)pp.Index.size();
		std::string s;
		if (index == 0)
		{
			for (int k = 0; k < indexLength; k++)
			{
				if (pp.Index[k][0] >= 0)
					s += seq1[pp.Index[k][0]];
				else
					s += '-';
			}
		}
		else
		{
			for (int k = 0; k < indexLength; k++)
			{
				if (pp.Index[k][1] >= 0)
					s += seq2[pp.Index[k][1]];
				else
					s += '-';
			}
		}
		return s;
	}


	template <class T>
	std::string Writer<T>::WriteFastq(std::shared_ptr<MonomerEncoding<T>> encoding, PolymerCollection<T>& pc)
	{
		
		std::string sb;
		for (auto &pp : pc)
		{
			Polymer<T> &p = *pp.second;
			sb += "@" + p.Name + "\n";
			sb += encoding->ToString(p.Seq) + "\n";
			sb += "+" + p.Name;
			for (auto &q : p.Qual)
			{

				sb += (char)(q - QScoreOffset);
			}
			sb += "\n";;
		}
		return sb;
	}

	template <class T>
	std::string Writer<T>::WriteNPMF(Polymer<std::vector<double>>& p)
	{
		std::string s;
		s += p.Name;
		for (int iPos = 0; iPos < (int)p.Seq.size(); iPos++)
		{
			s += "\t" + to_string(iPos + 1);
		}
		s += "\n";
		for (int iNuc = 0; iNuc < 5; iNuc++)
		{
			s += SymbolSets::PureNucleotides[iNuc];
			for (int iPos = 0; iPos < (int)p.Seq.size(); iPos++)
			{
				s += "\t" + to_string((p.Seq)[iPos][iNuc]);
			}
			s += "\n";
		}
		s += "Error";
		for (int iPos = 0; iPos < (int)p.Seq.size(); iPos++)
		{
			double error = 0;
			for (int iNuc = 0; iNuc < 5; iNuc++)
			{
				error -= (p.Seq)[iPos][iNuc] * (p.Seq)[iPos][iNuc];
			}
			error += 1;
			s += ("\t" + to_string(error));
		}
		s += "\n";
		return s;
	}

	template <class T>
	std::string Writer<T>::WriteNPMF(PolymerCollection<std::vector<double>>& pc)
	{
		std::string s;
		for (auto &p : pc)
		{
			s += WriteNPMF(*(p.second));
		}
		return s;
	}

	template <class T>
	std::string Writer<T>::WriteNPMF(PolymerPair<std::vector<double>, std::vector<double>>& pair)
	{
		std::string s;
		s += pair.Polymer0->Name;
		for (int k = 0; k < (int)pair.Index.size(); k++)
		{
			s += ("\t" + to_string(k + 1));
		}
		s += "\n";
		for (int iNuc = 0; iNuc < 5; iNuc++)
		{
			s += SymbolSets::PureNucleotides[iNuc];
			for (int k = 0; k < (int)pair.Index.size(); k++)
			{
				if (pair.Index[k][0] > -1)
				{
					s += "\t" + to_string((pair.Polymer0->Seq)[pair.Index[k][0]][iNuc]);
					//s.Append("\t" + pair.Polymer0.Seq[pair.Index[k][0]][iNuc]);
				}
				else
				{
					s += "\t";
					s += (iNuc == 4 ? "1" : "0");
				}
			}
			s += "\n";
		}
		s += ("Error");
		for (int k = 0; k < (int)pair.Index.size(); k++)
		{
			if (pair.Index[k][0] > -1)
			{
				double error = 1;
				for (int iNuc = 0; iNuc < 5; iNuc++)
				{
					error -= pow((pair.Polymer0->Seq)[pair.Index[k][0]][iNuc], 2);
				}
				s += ("\t" + to_string(error));
			}
			else
			{
				s += ("\t0");
			}
		}
		s += "\n";
		s += pair.Polymer1->Name;
		for (int k = 0; k < (int)pair.Index.size(); k++)
		{
			s += ("\t" + to_string(k + 1));
		}
		s += "\n";
		for (int iNuc = 0; iNuc < 5; iNuc++)
		{
			s += SymbolSets::PureNucleotides[iNuc];
			for (int k = 0; k < (int)pair.Index.size(); k++)
			{
				if (pair.Index[k][1] > -1)
				{
					s += "\t" + to_string((pair.Polymer1->Seq)[pair.Index[k][1]][iNuc]);

					//s.Append("\t" + pair.Polymer0.Seq[pair.Index[k][1]][iNuc]);
				}
				else
				{
					s += "\t";
					s += (iNuc == 4 ? "1" : "0");
				}
			}
			s += "\n";
		}
		s += ("Error");
		for (int k = 0; k < (int)pair.Index.size(); k++)
		{
			if (pair.Index[k][1] > -1)
			{
				double error = 1;
				for (int iNuc = 0; iNuc < 5; iNuc++)
				{
					error -= pow((pair.Polymer0->Seq)[pair.Index[k][1]][iNuc], 2);
				}
				s += "\t" + to_string(error);
			}
			else
			{
				s += ("\t0");
			}
		}
		s += "\n";
		return s;
	}


	template class Writer < std::vector<double> >;
}
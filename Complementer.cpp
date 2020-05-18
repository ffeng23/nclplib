#include "stdafx.h"

#include "Complementer.h"
#include "Monomer.h"
#include "Polymer.h"

namespace CLPLib
{

	template <class T>
	shared_ptr<Polymer<T>> Complementer<T>::Complement(Polymer<T>& p)
	{
		string newName = p.Name + ".C";
		shared_ptr<Polymer<T>> p2;
		if (p.Qual.size() > 0)
		{
			//copy and reverse qual
			std::vector<byte> &qual = p.Qual;
			std::vector<byte> tmp(qual.size());
			std::reverse_copy(std::begin(qual), std::end(qual), std::begin(tmp));
			p2 = std::make_shared<Polymer<T>>(nucleotideEncoding, newName, Complement(p.Seq), tmp);
		}
		else
		{
			p2 = std::make_shared<Polymer<T>>(nucleotideEncoding, newName, Complement(p.Seq));
		}
		return p2;
	}

	template <class T>
	std::vector<T> Complementer<T>::Complement(std::vector<T>& seq)
	{
		std::vector<T> newSeq(seq.size());
		int seqlen = (int)seq.size();
		for (int i = 0; i < seqlen; i++)
		{
			//auto b1 = TypeConverter::ToChar(seq[seqlen - 1 - i]);
			newSeq[i] = nucleotideEncoding->Complement(seq[seqlen - 1 - i]);
			//auto b2 = TypeConverter::ToChar(newSeq[i]);
			//int check = 1;
		}
		return newSeq;
	}




	template class Complementer<std::vector<double>>;
}
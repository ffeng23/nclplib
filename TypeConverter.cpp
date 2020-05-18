
#include "stdafx.h"

#include <memory>
#include "TypeConverter.h"
#include "Polymer.h"
#include "PolymerCollection.h"
#include "ByteBits.h"



namespace CLPLib
{

	//partial implementation of TypeConverter class

	map<char, int> TypeConverter::charToInt = PureNucleotideInt::encoder;
	vector<char> TypeConverter::intToChar = SymbolSets::PureNucleotides;
	vector<char> TypeConverter::byteToChar;
	map<char, byte> TypeConverter::charToByte = MixedNucleotideByte::encoder;
	double TypeConverter::aLittleBit = 1.0E-4;
	double TypeConverter::almostOne = 1 - aLittleBit;

	TypeConverter::TypeConverter()
	{
		byteToChar.resize(charToByte.size());
		const vector<char> &v = SymbolSets::MixedNucleotides;
		for (auto &iterator : v)
		{
			byteToChar[charToByte[iterator]] = iterator;
		}
	}

	Polymer<byte> TypeConverter::ToPolynucleotideByte(Polymer<vector<double>>& ppmf)
	{
		auto xtz = std::make_shared<MixedNucleotideByte>();
		return Polymer<byte>(std::make_shared<MixedNucleotideByte>(), ppmf.Name, vector<byte>(ToBytes(ppmf.Seq)), ppmf.Qual);
	}


	PolymerCollection<byte> TypeConverter::ToPolynucleotideByteCollection(PolymerCollection<vector<double>>& ppmfc)
	{
		PolymerCollection<byte> pbc(std::make_shared<MixedNucleotideByte>());
		pbc.Name = ppmfc.Name;
		for (auto &pd : ppmfc)
		{
			Polymer<byte> p = ToPolynucleotideByte(*(pd.second));
			pbc.AddMember(make_shared<Polymer<byte>>(p));
		}
		return pbc;
	}

	PolymerCollection<vector<double>> TypeConverter::ToPolynucleotidePMFCollection(PolymerCollection<byte>& pc)
	{
		PolymerCollection<vector<double>> ppmf(std::make_shared<NucleotidePMF>());
		ppmf.Name = pc.Name;
		for (auto &pb : pc)
		{
			Polymer<vector<double>> p = ToPolynucleotidePMF(*(pb.second));
			ppmf.AddMember(make_shared<Polymer<vector<double>>>(p));
		}

		return ppmf;
	}


	vector<double> TypeConverter::ToPMF(int i)
	{
		int len = (int)SymbolSets::PureNucleotides.size();
		vector<double> pmf(len, 0);
		if (i < 0 || i >= len)
		{
			return pmf;
		}
		else
		{
			pmf[i] = 1;
		}
		return pmf;
	}


	vector<double> TypeConverter::ToPMF(byte b)
	{
		if (b > 16)
		{
			throw std::runtime_error(b + " is not a valid nucleotide byte code.");
		}

		int dim = NucleotidePMF::Dim;
		vector<double> p(dim, 0);
		ByteBits bb(b);
		for (byte i = 0; i < dim; i++)
		{
			p[i] = bb.get_bit(i) ? (1.0 / MixedNucleotideByte::Degeneracy(b)) : 0;
		}
		return p;
	}

	vector<double> TypeConverter::ToPMF(byte b, byte q)
	{
		if (b > 16)
		{
			throw std::runtime_error(b + " is not a valid nucleotide byte code.");
		}

		double prob = pow(10.0, -q / 10.0);
		int dim = NucleotidePMF::Dim;
		vector<double> p(dim, 0);
		int degen = MixedNucleotideByte::Degeneracy(b);
		if (degen == 4)
		{
			for (byte i = 0; i < 4; i++)
			{
				p[i] = 0.25;
			}
		}
		else
		{
			ByteBits bb(b);
			for (byte i = 0; i < 4; i++)
			{
				p[i] = bb.get_bit(i) ? (1 - prob) / degen : prob / (4 - degen);
			}
		}
		return p;
	}


	Polymer<vector<double>> TypeConverter::ToPolynucleotidePMF(Polymer<byte>& p)
	{

		if (p.Qual.size() == 0)
		{
			return Polymer<vector<double>>(make_shared<NucleotidePMF>(), p.Name, vector<vector<double>>(ToPMFs(p.Seq, p.Qual)), p.Qual);

		}
		else
		{
			return Polymer<vector<double>>(make_shared<NucleotidePMF>(), p.Name, vector<vector<double>>(ToPMFs(p.Seq)));
		}
	}


	Polymer<vector<double>> TypeConverter::ToPolynucleotidePMF(Polymer<int>& p)
	{
		return Polymer<vector<double>>(make_shared<NucleotidePMF>(), p.Name, vector<vector<double>>(ToPMFs(p.Seq)), p.Qual);
	}

	Polymer<vector<double>> TypeConverter::AsPolynucleotidePMF(Polymer<byte>& p)
	{
		return ToPolynucleotidePMF(p);
	}

	Polymer<vector<double>> TypeConverter::AsPolymerPMF(Polymer<int>& p)
	{
		return ToPolynucleotidePMF(p);
	}


	//this is to initalize the static field of byteToChar
	TypeConverter converter;





}
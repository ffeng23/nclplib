#include "stdafx.h"
#include "Monomer.h"

namespace CLPLib
{
	//staitc variable initialization
	const vector<char> SymbolSets::PureNucleotides = { 'A', 'G', 'T', 'C', '-' };
	//												     0    1    2    3    4    5    6    7    8    9   10    11   12   13  14    15   16                              
	const vector<char> SymbolSets::MixedNucleotides = { ' ', 'A', 'G', 'R', 'T', 'W', 'K', 'D', 'C', 'M', 'S', 'V', 'Y', 'H', 'B', 'N', '-' };
	const vector<char> SymbolSets::AminoAcids = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*', '-' };

	const vector<char> PureNucleotideInt::symbols = SymbolSets::PureNucleotides;
	int PureNucleotideInt::complement[] = { 2, 3, 0, 1, 4 };
	int PureNucleotideInt::transition[] = { 1, 0, 3, 2, 4 };
	int PureNucleotideInt::transversion2[] = { 3, 2, 1, 0, 4 };
	map<int, char> PureNucleotideInt::decoder;
	map<char, int> PureNucleotideInt::encoder;

	const vector<char> MixedNucleotideByte::symbols = SymbolSets::MixedNucleotides;
	map<byte, char> MixedNucleotideByte::decoder;
	map<char, byte> MixedNucleotideByte::encoder;
	vector<byte> MixedNucleotideByte::complement = { 0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 16 };
	vector<int> MixedNucleotideByte::degeneracy = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1 };
	double MixedNucleotideByte::innerProduct[NSymbols][NSymbols];

	vector<char> MixedNucleotideChar::symbolComplements = { ' ', 'T', 'C', 'Y', 'A', 'W', 'M', 'H', 'G', 'K', 'S', 'B', 'R', 'D', 'V', 'N', '-' };


	MixedNucleotideByte NucleotidePMF::mnbEncoding;



}

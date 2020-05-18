#pragma once

#include <string>
#include <map>
#include "Monomer.h"

namespace CLPLib
{
	/// <summary>
	/// Provides tools for the translation of polynucleotide sequences into polypeptide sequences and the inverse.
	/// </summary>
	class Translator
	{
	public:
		/// <summary>
		/// std::map mapping codons to amino acids.
		/// </summary>
		static std::map<std::string, char> CodonTable;

		/// <summary>
		/// std::map mapping amino acids to codons. Ambiguities are handled by using the IUPAC codes where possible, and by using the specialized codes
		/// 1 = first nucleotide in Leuceine codon
		/// 2 = last nucleotide in Leuceine codon
		/// 3 = first nucleotide in Arginine codon
		/// 4 = last nucleotide in Argeinine codon
		/// 5,6,7 = (ordered) nucleotides in serine codon
		/// </summary>
		static std::map<char, std::string> ReverseCodonTable;

		//simulating static constructor of c#
		friend class constructor;
		class constructor
		{
		public:
			constructor()
			{
				initializeTranslators();
			}
		};
		static constructor initializer;
		//simulation ends

	private:
		static std::map<char, int> encoder;

		Translator()
		{
		}

		static void initializeEncoder()
		{
			for (unsigned int i = 0; i < SymbolSets::AminoAcids.size(); i++)
			{
				encoder.insert(std::make_pair(SymbolSets::AminoAcids[i], i));
			}
		}

		static void initializeTranslators()
		{
			//CodonTable = new std::map<std::string, char>();

			CodonTable.insert(std::make_pair("AAA", 'K'));
			CodonTable.insert(std::make_pair("AAC", 'N'));
			CodonTable.insert(std::make_pair("AAG", 'K'));
			CodonTable.insert(std::make_pair("AAT", 'N'));
			CodonTable.insert(std::make_pair("ACA", 'T'));
			CodonTable.insert(std::make_pair("ACC", 'T'));
			CodonTable.insert(std::make_pair("ACG", 'T'));
			CodonTable.insert(std::make_pair("ACT", 'T'));
			CodonTable.insert(std::make_pair("AGA", 'R'));
			CodonTable.insert(std::make_pair("AGC", 'S'));
			CodonTable.insert(std::make_pair("AGG", 'R'));
			CodonTable.insert(std::make_pair("AGT", 'S'));
			CodonTable.insert(std::make_pair("ATA", 'I'));
			CodonTable.insert(std::make_pair("ATC", 'I'));
			CodonTable.insert(std::make_pair("ATG", 'M'));
			CodonTable.insert(std::make_pair("ATT", 'I'));
			CodonTable.insert(std::make_pair("CAA", 'Q'));
			CodonTable.insert(std::make_pair("CAC", 'H'));
			CodonTable.insert(std::make_pair("CAG", 'Q'));
			CodonTable.insert(std::make_pair("CAT", 'H'));
			CodonTable.insert(std::make_pair("CCA", 'P'));
			CodonTable.insert(std::make_pair("CCC", 'P'));
			CodonTable.insert(std::make_pair("CCG", 'P'));
			CodonTable.insert(std::make_pair("CCT", 'P'));
			CodonTable.insert(std::make_pair("CGA", 'R'));
			CodonTable.insert(std::make_pair("CGC", 'R'));
			CodonTable.insert(std::make_pair("CGG", 'R'));
			CodonTable.insert(std::make_pair("CGT", 'R'));
			CodonTable.insert(std::make_pair("CTA", 'L'));
			CodonTable.insert(std::make_pair("CTC", 'L'));
			CodonTable.insert(std::make_pair("CTG", 'L'));
			CodonTable.insert(std::make_pair("CTT", 'L'));
			CodonTable.insert(std::make_pair("GAA", 'E'));
			CodonTable.insert(std::make_pair("GAC", 'D'));
			CodonTable.insert(std::make_pair("GAG", 'E'));
			CodonTable.insert(std::make_pair("GAT", 'D'));
			CodonTable.insert(std::make_pair("GCA", 'A'));
			CodonTable.insert(std::make_pair("GCC", 'A'));
			CodonTable.insert(std::make_pair("GCG", 'A'));
			CodonTable.insert(std::make_pair("GCT", 'A'));
			CodonTable.insert(std::make_pair("GGA", 'G'));
			CodonTable.insert(std::make_pair("GGC", 'G'));
			CodonTable.insert(std::make_pair("GGG", 'G'));
			CodonTable.insert(std::make_pair("GGT", 'G'));
			CodonTable.insert(std::make_pair("GTA", 'V'));
			CodonTable.insert(std::make_pair("GTC", 'V'));
			CodonTable.insert(std::make_pair("GTG", 'V'));
			CodonTable.insert(std::make_pair("GTT", 'V'));
			CodonTable.insert(std::make_pair("TAA", '*'));
			CodonTable.insert(std::make_pair("TAC", 'Y'));
			CodonTable.insert(std::make_pair("TAG", '*'));
			CodonTable.insert(std::make_pair("TAT", 'Y'));
			CodonTable.insert(std::make_pair("TCA", 'S'));
			CodonTable.insert(std::make_pair("TCC", 'S'));
			CodonTable.insert(std::make_pair("TCG", 'S'));
			CodonTable.insert(std::make_pair("TCT", 'S'));
			CodonTable.insert(std::make_pair("TGA", '*'));
			CodonTable.insert(std::make_pair("TGC", 'C'));
			CodonTable.insert(std::make_pair("TGG", 'W'));
			CodonTable.insert(std::make_pair("TGT", 'C'));
			CodonTable.insert(std::make_pair("TTA", 'L'));
			CodonTable.insert(std::make_pair("TTC", 'F'));
			CodonTable.insert(std::make_pair("TTG", 'L'));
			CodonTable.insert(std::make_pair("TTT", 'F'));

			//ReverseCodonTable = new std::map<char, std::string>();
			// note the special codes:
			// 1 = first nucleotide in Leuceine codon
			// 2 = last nucleotide in Leuceine codon
			// 3 = first nucleotide in Arginine codon
			// 4 = last nucleotide in Argeinine codon
			// 5,6,7 = (ordered)) nucleotides in serine codon
			ReverseCodonTable.insert(std::make_pair('F', "TTY"));
			ReverseCodonTable.insert(std::make_pair('L', "1T2"));
			ReverseCodonTable.insert(std::make_pair('I', "ATH"));
			ReverseCodonTable.insert(std::make_pair('M', "ATG"));
			ReverseCodonTable.insert(std::make_pair('V', "GTN"));
			ReverseCodonTable.insert(std::make_pair('S', "567"));
			ReverseCodonTable.insert(std::make_pair('P', "CCN"));
			ReverseCodonTable.insert(std::make_pair('T', "CAN"));
			ReverseCodonTable.insert(std::make_pair('A', "GCN"));
			ReverseCodonTable.insert(std::make_pair('Y', "TAY"));
			ReverseCodonTable.insert(std::make_pair('H', "CAY"));
			ReverseCodonTable.insert(std::make_pair('Q', "CAR"));
			ReverseCodonTable.insert(std::make_pair('N', "AAY"));
			ReverseCodonTable.insert(std::make_pair('K', "AAR"));
			ReverseCodonTable.insert(std::make_pair('D', "GAY"));
			ReverseCodonTable.insert(std::make_pair('E', "GAR"));
			ReverseCodonTable.insert(std::make_pair('C', "TGY"));
			ReverseCodonTable.insert(std::make_pair('W', "TGG"));
			ReverseCodonTable.insert(std::make_pair('R', "3G4"));
			ReverseCodonTable.insert(std::make_pair('G', "GGN"));
		}

		/// <summary>
		/// Performs the translation of a three-letter codon into a single-letter amino acid.
		/// </summary>
		/// <param name="s">The codon to be translated.</param>
		/// <returns>If s represents a valid codon, the correct single-letter amino acid code is returned. If not, returns "X".</returns>
		static char translate(std::string s)
		{
			if (CodonTable.find(s) != CodonTable.end())
			{
				return CodonTable[s];
			}
			else return 'X';
		}


	public:
		/// <summary>
		/// Translates a polynucleotide std::string into a polypeptide std::string.
		/// </summary>
		/// <param name="s">The polynucleotide std::string.</param>
		/// <returns>The amino acid sequence that corresponds to the translation of s.</returns>
		static std::string Translate(std::string s)
		{
			std::string translation;
			for (unsigned int i = 0; i < s.length() - 2; i += 3)
			{
				translation += (translate(s.substr(i, 3)));
			}
			return translation;
		}

		static vector<double> Translate(std::vector<std::vector<double>>& codonpmf)
		{
			// codonpmf is an array comprising three independent nucleotide pmfs
			vector<double> apmf(SymbolSets::AminoAcids.size());
			for (int i0 = 0; i0 < 4; i0++)
			{
				for (int i1 = 0; i1 < 4; i1++)
				{
					for (int i2 = 0; i2 < 4; i2++)
					{
						std::string codon = std::string(1, SymbolSets::PureNucleotides[i0]) +
							std::string(1, SymbolSets::PureNucleotides[i1]) +
							std::string(1, SymbolSets::PureNucleotides[i2]);
						char aminoAcid = translate(codon);
						int iaminoacid = encoder[aminoAcid];
						apmf[iaminoacid] += codonpmf[0][i0] * codonpmf[1][i1] * codonpmf[2][i2];
					}
				}
			}
			return apmf;
		}
	};
}

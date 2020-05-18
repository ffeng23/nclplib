#pragma once

#include<map>
#include<vector>
#include<list>
#include<math.h>
#include<algorithm>
#include<string>
#include<exception>
#include<stdexcept>

#include "typedef.h"

using namespace std;

namespace CLPLib
{
	template <class T> class Polymer;
	template <class T> class PolymerCollection;

	class TypeConverter
	{
		static map<char, int> charToInt;		// = PureNucleotideInt.Encoder;
		static vector<char> intToChar;			// = SymbolSets.PureNucleotides;
		static vector<char> byteToChar;


		static map<char, byte> charToByte;		// = MixedNucleotideByte.Encoder;
		static double aLittleBit;				// = 1.0E-4;
		static double almostOne;				// = 1 - aLittleBit;

	public:
		TypeConverter();

		static byte ToByte(vector<double>& pmf)
		{
			double sum = 0;
			for (auto &d : pmf) sum += d;
			if (sum < almostOne)
				return 15; // the pmf is being used to hold the place for N nucleotides			
			vector<pair<byte, double>> cArray;
			for (byte i = 0; i < 4; i++)
			{
				cArray.push_back(pair<byte, double>((byte)pow(2, i), pmf[i]));
			}

			sort(cArray.begin(), cArray.end(), [](const pair<byte, double> &a, const pair<byte, double> &b)
			{
				return a.second > b.second;
			});

			byte nucByte = 0;

			if (pmf[4] > 0.5)
				return 16;

			double norm = 1 - pmf[4];

			if (cArray[0].second - cArray[1].second > 0.5 * norm)
			{
				nucByte = cArray[0].first;
			}
			else if (cArray[0].second + cArray[1].second - 2 * cArray[2].second > 0.5 * norm)
			{
				nucByte = (byte)(cArray[0].first + cArray[1].first);
			}
			else if (cArray[3].second > 0.125 * norm)
			{
				nucByte = 15;
			}
			else
			{
				nucByte = (byte)(cArray[0].first + cArray[1].first + cArray[2].first);
			}
			return nucByte;
		}

		static vector<byte> ToBytes(vector<vector<double> >& pmfs)
		{
			vector<byte> bytes;
			for (auto &pmf : pmfs)
			{
				bytes.push_back(ToByte(pmf));
			}
			return bytes;
		}

		static vector<byte> ToBytes(string sequence)
		{
			vector<byte> bytes;
			for (char &c : sequence)
			{
				bytes.push_back(charToByte[c]);
			}
			return bytes;
		}

		static Polymer<byte> ToPolynucleotideByte(Polymer<vector<double>>& ppmf);

		static PolymerCollection<byte> ToPolynucleotideByteCollection(PolymerCollection<std::vector<double>>& ppmfc);

		static PolymerCollection<vector<double>> ToPolynucleotidePMFCollection(PolymerCollection<byte>& pc);

		static vector<double> ToPMF(int i);

		static vector<vector<double>> ToPMFs(vector<int> &seq)
		{
			vector<vector<double>> pmfs;
			for (auto &s : seq)
			{
				pmfs.push_back(ToPMF(s));
			}
			return pmfs;
		}

		static vector<double> ToPMF(byte b);

		static vector<vector<double>> ToPMFs(vector<byte>& bytes)
		{
			vector< vector<double> > pmfs;
			for (auto& byte : bytes)
			{
				pmfs.push_back(ToPMF(byte));
			}
			return pmfs;
		}

		static vector<double> ToPMF(byte b, byte q);

		static vector<vector<double>> ToPMFs(vector<byte> &bytes, vector<byte> & qual)
		{
			if (bytes.size() != qual.size())
			{
				throw std::runtime_error("The sequence and quality score arrays are not the same length.");
			}

			vector<vector<double>> pmfs;
			for (size_t i = 0; i < bytes.size(); i++)
			{
				pmfs[i] = ToPMF(bytes[i], qual[i]);
			}
			return pmfs;
		}

		static vector<vector<double>> ToPMFs(string sequence)
		{
			auto intermediate = ToBytes(sequence);
			return ToPMFs(intermediate);
		}

		static Polymer<vector<double>> ToPolynucleotidePMF(Polymer<byte>& p);

		static Polymer<vector<double>> ToPolynucleotidePMF(Polymer<int>& p);

		static Polymer<vector<double>> AsPolynucleotidePMF(Polymer<byte>& p);

		static Polymer<vector<double>> AsPolymerPMF(Polymer<int>& p);

		static char ToChar(int i)
		{
			return intToChar[i];
		}

		static char ToChar(byte b)
		{
			return byteToChar[b];
		}

		static char ToChar(vector<double> p)
		{
			byte intermediate = ToByte(p);
			return ToChar(intermediate);
		}

		static vector<char> ToChars(const vector<byte>& seq)
		{
			vector<char> chars(seq.size());
			for (size_t i = 0; i < seq.size(); i++)
			{
				chars[i] = byteToChar[seq[i]];
			}
			return chars;
		}

		static vector<char> ToChars( vector<vector<double>>& seq)
		{
			auto intermediate = ToBytes(seq);
			return ToChars(intermediate);
		}
	};
}

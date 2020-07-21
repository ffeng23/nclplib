#pragma once

#include <stdlib.h>
#include <stdexcept>
#include <map>
#include <vector>
#include "typedef.h"
#include "TypeConverter.h"
#include "ByteBits.h"

using namespace std;

namespace CLPLib
{

	enum MonomerType { Nucleotide, AminoAcid };
	enum RepresentationType { Pure, Mixed, PMF };

	class SymbolSets
	{
	public:
		static const vector<char> PureNucleotides;
		static const vector<char> MixedNucleotides;
		static const vector<char> AminoAcids;
	};

	template <class T>
	class MonomerEncoding
	{
	protected:

		//System.Type does not seem have a C++ equivalent, and it does not seem to be used anywhere 
		//System.Type basicType;
		RepresentationType repType;
	public:
		virtual bool IsValid(T& monomer) = 0;
		virtual bool IsGap(T& monomer) = 0;
		virtual T Complement(T& monomer) = 0;
		virtual T Copy(T& monomer) = 0;
		virtual T FromChar(char c) = 0;
		virtual char ToChar(T& monomer) = 0;

		virtual string ToString(vector<T>& sequence) = 0;

		virtual vector<T> FromString(string& sequence) = 0;

		virtual T Gap() = 0;
	};

	template<class T>
	class NucleotideEncoding : public MonomerEncoding < T >
	{
	public:
		static const double LogPrior;	// = -1.38629436; = Log(4)
	};
	template<class T> const double NucleotideEncoding<T>::LogPrior = -1.38629436;


	class AminoAcidChar : public MonomerEncoding < char >
	{
	public:
		char Complement(char& monomer)override
		{
			throw std::logic_error("NotImplementedException");
		}

		char Copy(char& monomer)override
		{
			throw std::logic_error("NotImplementedException");
		}

		char FromChar(char c) override
		{
			return c;
		}

		bool IsGap(char& monomer) override
		{
			return monomer == '-';
		}

		bool IsValid(char& monomer) override
		{
			if (monomer == '*' || monomer == '-')return false;
			const vector<char> &v = SymbolSets::AminoAcids;
			return std::find(v.begin(), v.end(), monomer) != v.end();
		}

		char ToChar(char& monomer) override
		{
			return monomer;
		}

		string ToString(vector<char> &sequence) override
		{
			string str(sequence.begin(), sequence.end());
			return str;
		}

		vector<char> FromString(string& sequence) override
		{
			vector<char> retval(sequence.begin(), sequence.end());
			return retval;
		}

		char Gap()override
		{
			return '-';
		}
	};


	class AminoAcidPMF : public MonomerEncoding < vector<double> >
	{
	private:
		static std::map<char, int> encoder;


	public:
		unsigned int dim = SymbolSets::AminoAcids.size();

		//c++ does not have static constructor, need to constrcut an instance
		//to initialize the static member later.
		AminoAcidPMF()
		{
			const vector<char> &v = SymbolSets::AminoAcids;
			for (int i = 0; i < (int)v.size(); i++)
			{
				encoder.insert(std::pair<char, int>(v[i], i));
			}
		}

		vector<double> Complement(vector<double>& monomer)override
		{
			throw std::logic_error("NotImplementedException");
		}

		vector<double> Copy(vector<double>& monomer)override
		{
			throw std::logic_error("NotImplementedException");
		}

		vector<double> FromChar(char c) override
		{
			vector<double> pmf(dim, 0);
			if (encoder.count(c) > 0)
			{
				pmf[encoder[c]] = 1;
			}
			return pmf;
		}

		vector<vector<double>> FromString(string &sequence) override
		{
			vector<vector<double>> vector;
			for (size_t i = 0; i < sequence.length(); i++)
			{
				vector.push_back(FromChar(sequence[i]));
			}
			return vector;
		}

		bool IsGap(vector<double>& a) override
		{
			return a[21] == 1;
		}

		bool IsValid(vector<double>& a) override
		{
			//we cannot check length here
			//if (a.Length != Dim)return false;
			double sum = 0;
			for (unsigned int i = 0; i < dim; i++)
			{
				if (a[i] < 0) return false;
				sum += a[i];
			}
			//this test may need to be improved.
			if (sum <0.99999999 || sum > 1.00000001) //if (sum != 1)
				return false;

			return true;
		}

		///now we don't have the array length informaiton....
		string ToString(vector<vector<double>> &sequence) override
		{
			string retstr;
			for (size_t i = 0; i < sequence.size(); i++)
			{
				retstr.push_back(ToChar(sequence[i]));

			}
			return retstr;
		}

		char ToChar(vector<double>& a) override
		{

			if (a.size() != dim)
				throw std::logic_error("The length of the argument is incorrect.");
			double maxval = a[0];
			int maxi = 0;
			for (unsigned int i = 1; i < dim; i++)
			{
				if (a[i] > maxval)
				{
					maxval = a[i];
					maxi = i;
				}
			}
			return SymbolSets::AminoAcids[maxi];
		}

		vector<double> Gap() override
		{
			vector<double> gap(22, 0);
			gap[21] = 1;
			return gap;
		}
	};

	/// <summary>
	/// Encoding of the pure nucleotides A, G, T, C by the integers 0, 1, 2, 3.
	/// </summary>
	class PureNucleotideInt : public NucleotideEncoding < int >
	{
		// A : 0
		// G : 1
		// T : 2
		// C : 3
		// - : 4
		// gap included for convenience
		// n1, n2 related by transition iff n1/2 = n2/2 in integer arithmetic
		// n1, n2 complementary iff n1 % 2 = n2 % 2

	private:

		static int complement[];
		static int transition[];
		static int transversion2[];
		static std::map<int, char> decoder;


	public:
		static const vector<char> symbols;
		static std::map<char, int> encoder;




		//c++ don't have static constructor, create a static member somewhere before use
		PureNucleotideInt()
		{
			if (encoder.size() == 0)
			{
				const vector<char> &v = symbols;
				for (int i = 0; i < (int)v.size(); i++)
				{
					encoder.insert(std::pair<char, int>(v[i], i));
					decoder.insert(std::pair<int, char>(i, v[i]));
				}
				decoder.insert(std::pair<int, char>(-1, 'x'));
			}
		}
		/// <summary>
		/// Gets the Watson-Crick complement of the argument
		/// </summary>
		/// <param name="n">An integer representing a pure nucleotide.</param>
		/// <returns>The integer representing the Watson-Crick complement of the argument. Returns -1 if the argument is illegitimate.</returns>
		int Complement(int& n) override
		{
			if (n < 0 || n > 4)
				return -1;
			else
				return complement[n];
		}

		int Copy(int& n) override
		{
			return n;
		}

		/// <summary>
		/// Determines whether a given variable of the appropriate type represents a gap.
		/// </summary>
		/// <param name="n">An integer representing a pure nucleotide.</param>
		/// <returns>True if and only if n represets a gap.</returns>
		bool IsGap(int& n) override
		{
			return n == 4;
		}

		/// <summary>
		/// Determines whether a given integer represents a pure nucleotide.
		/// </summary>
		/// <param name="n">An integer representing a pure nucleotide.</param>
		/// <returns>True if and only if n represets a pure nucleotide.</returns>
		bool IsValid(int& n) override
		{
			return (n > -1 && n < 5);
		}

		/// <summary>
		/// Gets the nucleotide that is the transition partner of the argument.
		/// </summary>
		/// <param name="n">A nucleotide in integer representation.</param>
		/// <returns>The nucleotide that represents a transition mutation from the argument nucleotide, or -1 if the argument is illegitimate.</returns>
		static int Transition(int n)
		{
			if (n < 0 || n > 4)
				return -1;
			else
				return transition[n];
		}

		/// <summary>
		/// Gets the nucleotide that is the transversion2 (not-complement) partner of the argument.
		/// </summary>
		/// <param name="n">A nucleotide in integer representation.</param>
		/// <returns>The nucleotide that represents a transversion-not-complement mutation from the argument nucleotide, or -1 if the argument is illegitimate.</returns>
		static int Transversion2(int n)
		{
			if (n < 0 || n > 4)
				return -1;
			else
				return transversion2[n];
		}

		int FromChar(char symbol) override
		{
			for (size_t i = 0; i < symbols.size(); i++)
			{
				if (symbols[i] == symbol)return encoder[symbol];
			}
			return -1;
		}

		char ToChar(int& monomer) override
		{
			return TypeConverter::ToChar(monomer);
		}

		string ToString(vector<int>& sequence) override
		{
			string str;
			for (auto &n : sequence)
			{
				str.push_back(decoder[n]);
			}
			return str;
		}

		vector<int> FromString(string& sequence) override
		{
			vector<int> array(sequence.length(), 0);
			for (size_t i = 0; i < sequence.length(); i++)
			{
				if (encoder.count(sequence[i]) > 0)
					array[i] = encoder[sequence[i]];
				else
					array[i] = -1;
			}
			return array;
		}

		int Gap() override
		{
			return 4;
		}
	};


	/// <summary>
	/// Provides support for for nucleotides encoded by bytes.
	/// </summary>
	class MixedNucleotideByte : public NucleotideEncoding < byte >
	{
	private:
		static const vector<char> symbols;
		static map<byte, char> decoder;
		static vector<byte> complement;
		static vector<int> degeneracy;
	public:
		static map<char, byte> encoder;
		static const int NSymbols = 17;
	private:
		static double innerProduct[NSymbols][NSymbols];

	public:
		MixedNucleotideByte()
		{
			if (encoder.size() == 0)
			{
				for (byte b = 0; b < NSymbols; b++)
				{
					encoder.insert(std::pair<char, int>(symbols[b], b));
					decoder.insert(std::pair<int, char>(b, symbols[b]));
				}
				computeInnerProducts();
			}
		}

	private:
		void computeInnerProducts()
		{
			for (int i = 0; i < NSymbols; i++)
			{
				for (int j = 0; j < NSymbols; j++)
				{
					ByteBits b3((byte)(i & j));
					innerProduct[i][j] = b3.NSet() / ((double)degeneracy[i] * degeneracy[j]);
				}
			}
		}
	public:
		/// <summary>
		/// Reports whether a given byte legitimately encodes a nucleotide.
		/// </summary>
		/// <param name="n">The nucleotide to be tested.</param>
		/// <returns>True if n encodes a nucleotide. False otherwise.</returns>
		bool IsValid(byte& n)override
		{
			return (n < 17);
		}

		/// <summary>
		/// Computes the complement of a byte-encoded nucleotide.
		/// </summary>
		/// <param name="n">The byte-encoded nucleotide to be complemenented.</param>
		/// <returns>The Watson-Crick complement of n. If n is not valid, returns 0. </returns>
		byte Complement(byte& n) override
		{
			if (IsValid(n))
			{
				return complement[n];
			}
			else
			{
				return 0;
			}
		}

		byte Copy(byte& n) override
		{
			return n;
		}

		/// <summary>
		/// Computes the IUPAC degeneracy of a byte-encoded nucleotide symbol.
		/// </summary>
		/// <param name="n">A byte-encoded nucleotide symbol.</param>
		/// <returns>The number of pure nucleotides consistent with n.</returns>
		static int Degeneracy(byte n)
		{
			if (n < 17)
			{
				return degeneracy[n];
			}
			else
			{
				return 0;
			}
		}

		/// <summary>
		/// Determines whether a nucleotide symbol encodes a gap.
		/// </summary>
		/// <param name="n">A byte-encoded nucleotide symbol.</param>
		/// <returns>True if n encodes a gap. False otherwise, including if n is illegitimate..</returns>
		bool IsGap(byte& n) override
		{
			return n == 16;
		}

		/// <summary>
		/// Determines whether two byte-encoded nucleotide symbols are consistent.
		/// </summary>
		/// <param name="n1">A byte-encoded nucleotide symbol.</param>
		/// <param name="n2">A byte-encoded nucleotide symbol.</param>
		/// <returns>True if there is a nucleotide that is consistent with both n1 and n2. False otherwise. The result is meanngless if either n1 or n2 is illegitimate.</returns>
		static bool AreConsistent(byte n1, byte n2)
		{
			return (n1 & n2) > 0;
		}

		static double InnerProduct(byte n1, byte n2)
		{
			if (n1 > 17 || n2 > 17) return std::numeric_limits<double>::quiet_NaN();

			return innerProduct[n1][n2];
		}

		static bool RepresentsPureState(byte n)
		{
			return (Degeneracy(n) == 1);
		}

		static int ToPureState(byte n)
		{
			switch (n)
			{
			case 1:
				return 0;
			case 2:
				return 1;
			case 4:
				return 2;
			case 8:
				return 3;
			default:
				return -1;
			}
		}

		byte FromChar(char symbol) override
		{
			if (find(symbols.begin(), symbols.end(), symbol) != symbols.end())
			{
				return encoder[symbol];
			}
			return 0;
		}

		char ToChar(byte& monomer) override
		{
			return TypeConverter::ToChar(monomer);
		}

		vector<byte> FromString(string& sequence)override
		{
			vector<byte> array(sequence.length(), 0);
			for (size_t i = 0; i < sequence.length(); i++)
			{
				if (encoder.count(sequence[i]) > 0)
				{
					array[i] = encoder[sequence[i]];
				}
			}
			return array;
		}

		string ToString(vector<byte>& sequence) override
		{
			string s;
			for (auto &b : sequence)
			{
				s.push_back(decoder[b]);
			}
			return s;
		}

		byte Gap() override
		{
			return 16;
		}
	};


	/// <summary>
	/// Provides support for for nucleotides encoded by bytes.
	/// </summary>
	class MixedNucleotideChar : public NucleotideEncoding < char >
	{
	private:
		static vector<char> symbols;
		static vector<char> symbolComplements;
		static map<char, char> complement;

	public:
		static const int NSymbols = 17;

		MixedNucleotideChar()
		{
			//initialize only once.
			if (complement.size() == 0)
			{
				for (size_t i = 0; i < symbols.size(); i++)
				{
					complement.insert(std::make_pair(symbols[i], symbolComplements[i]));
				}
			}
		}

		/// <summary>
		/// Reports whether a given byte legitimately encodes a nucleotide.
		/// </summary>
		/// <param name="n">The nucleotide to be tested.</param>
		/// <returns>True if n encodes a nucleotide. False otherwise.</returns>
		bool IsValid(char& c) override
		{
			return (find(symbols.begin(), symbols.end(), c) != symbols.end());
		}

		/// <summary>
		/// Computes the complement of a byte-encoded nucleotide.
		/// </summary>
		/// <param name="n">The byte-encoded nucleotide to be complemenented.</param>
		/// <returns>The Watson-Crick complement of n. If n is not valid, returns 0. </returns>
		char Complement(char& c) override
		{
			if (IsValid(c))
			{
				return complement[c];
			}
			else
			{
				return ' ';
			}
		}

		char Copy(char& c) override
		{
			return c;
		}

		/// <summary>
		/// Determines whether a nucleotide symbol encodes a gap.
		/// </summary>
		/// <param name="n">A byte-encoded nucleotide symbol.</param>
		/// <returns>True if n encodes a gap. False otherwise, including if n is illegitimate..</returns>
		bool IsGap(char& c) override
		{
			return c == '-';
		}

		char FromChar(char symbol) override
		{
			return symbol;
		}

		char ToChar(char& monomer) override
		{
			return monomer;
		}

		string ToString(vector<char>& sequence) override
		{
			throw std::runtime_error("NotImplementedException");
		}

		vector<char> FromString(string& sequence) override
		{
			throw std::runtime_error("NotImplementedException");
		}

		char Gap() override
		{
			return '-';
		}
	};


	/// <summary>
	/// Provides support for nucleotides as probability mass functions (PMF), that is, double[].
	/// </summary>
	class NucleotidePMF : public NucleotideEncoding < vector<double> >
	{
	public:
		static const unsigned int Dim = 5;
		static MixedNucleotideByte mnbEncoding;
		/// <summary>
		/// Constructor. Doesn't do anything.
		/// </summary>
		NucleotidePMF()
		{
		}

		/// <summary>
		/// Determines whether a double[] encodes a legitimate nucleotide (PMF).
		/// </summary>
		/// <param name="n">A double array.</param>
		/// <returns>True if n encodes a legitimate nucleotidet PMF; false otherwise.</returns>
		bool IsValid(vector<double>& n) override
		{
			if (n.size() != Dim)
				return false;
			double sum = 0;
			for (unsigned int i = 0; i < Dim; i++)
			{
				if (n[i] < 0)
					return false;
				sum += n[i];
			}

			if (sum != 1)
				return false;

			return true;
		}

		/// <summary>
		/// The complement of a double[] variable treated as a nucleotide PMF.
		/// </summary>
		/// <param name="n">A double array.</param>
		/// <returns>The Watson-Crick complement of n. Return value is meaningless if n is illegitimate.</returns>
		vector<double> Complement(vector<double>& n) override
		{
			if (n.size() != 5)
			{
				throw std::runtime_error("array size exception");
			}
			vector<double> comp(Dim);
			comp[0] = n[2];
			comp[1] = n[3];
			comp[2] = n[0];
			comp[3] = n[1];
			comp[4] = n[4];
			return comp;
		}

		vector<double> Copy(vector<double>& n) override
		{
			if (n.size() != 5)
			{
				throw std::runtime_error("array size exception");
			}
			vector<double> comp(Dim);
			comp[0] = n[0];
			comp[1] = n[1];
			comp[2] = n[2];
			comp[3] = n[3];
			comp[4] = n[4];
			return comp;
		}
		/// <summary>
		/// Determines whether a double array represents a gap state.
		/// </summary>
		/// <param name="n">A double[] variable representing a nucleotide symbol.</param>
		/// <returns>True if the argument represents a pure gap state, False otherwise. Throws an exception of n.Length is less than 5</returns>
		bool IsGap(vector<double>& n)override
		{
			return n[4] == 1;
		}

		vector<double> FromChar(char c)override
		{
			byte b = mnbEncoding.FromChar(c);
			return TypeConverter::ToPMF(b);
		}

		vector<double> FromCharQual(char c, byte q)
		{
			byte b = mnbEncoding.FromChar(c);
			return TypeConverter::ToPMF(b, q);
		}

		char ToChar(vector<double>& monomer) override
		{
			return TypeConverter::ToChar(monomer);
		}

		vector<vector<double>> FromString(string& sequence)override
		{
			vector<vector<double>> seq;
			seq.reserve(sequence.length());
			for (auto& s : sequence)
			{
				auto tmp = FromChar(s);
				seq.push_back(tmp);
			}
			return seq;
		}

		vector<vector<double>> FromStringQual(string& sequence, vector<byte>& qual)
		{
			vector<vector<double>> seq;
			seq.reserve(sequence.length());
			for (int i = 0; i < (int)sequence.length(); i++)
			{
				seq.push_back(FromCharQual(sequence[i], qual[i]));
			}
			return seq;
		}

		string ToString(vector<vector<double>>& sequence)override
		{
			vector<char> seqvec = TypeConverter::ToChars(sequence);
			string seq = std::string(seqvec.begin(), seqvec.end());
			return seq;
		}

		static double IDFlux(const vector<double>& n1, const vector<double>& n2)
		{
			if (n1.size() != Dim || n2.size() != Dim)
				throw std::runtime_error("ArgumentException");

			return n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2] + n1[3] * n2[3] + n1[4] * n2[4];
		}

		static double TransitionFlux(vector<double>& n1, vector<double>& n2)
		{
			if (n1.size() != Dim || n2.size() != Dim)
				throw std::runtime_error("ArgumentException");

			return n1[0] * n2[1] + n1[1] * n2[0] + n1[2] * n2[3] + n1[3] * n2[2];
		}

		static double TransversionFlux(vector<double>& n1, vector<double>& n2)
		{
			if (n1.size() != Dim || n2.size() != Dim)
				throw std::runtime_error("ArgumentException");

			return n1[0] * (n2[2] + n2[3]) + n1[1] * (n2[2] + n2[3]) + n1[2] * (n2[0] + n2[1]) + n1[3] * (n2[0] + n2[1]);
		}


		static double GapFlux(vector<double>& n1, vector<double>& n2)
		{
			if (n1.size() != Dim || n2.size() != Dim)
				throw std::runtime_error("ArgumentException");
			return (1 - n1[4]) * n2[4] + (1 - n2[4]) * n1[4];
		}


		vector<double> Gap() override
		{
			return vector < double > { 0, 0, 0, 0, 1 };
		}
	};

}



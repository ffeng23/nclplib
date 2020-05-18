#pragma once

#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>

#include "typedef.h"
#include "Mutations.h"
#include "Monomer.h"


namespace CLPLib
{
	class GeneReplicator
	{
	private:
		//note: c# const is implicitly static
		static const int halfWidth = 2;
		static int width; // = 2 * halfWidth + 1;

		//static Troschuetz.Random.Generators.ALFGenerator mt = new Troschuetz.Random.Generators.ALFGenerator();
		//static Troschuetz.Random.Distributions.Discrete.BernoulliDistribution bern = new Troschuetz.Random.Distributions.Discrete.BernoulliDistribution(mt);

		static std::vector<std::vector<byte>> motifs;

		static std::unordered_map<std::vector<byte>, std::vector<double>>   spectrumProbabilities;
	public:
		static std::unordered_map<std::vector<byte>, std::unordered_map<SubstitutionType, double>> SpectrumDistribution;
		
		
		//static Dictionary<byte[], TBKMath.GeneralDiscreteDistribution<int>> SpectrumProbabilities;

		
		static std::unordered_map<std::vector<byte>, double> baseMutabilities;


		static bool matches(std::vector<byte>& m1, std::vector<byte>& m2)
		{
			if (m1.size() != m2.size()) return false;

			for (int iPos = 0; iPos < m1.size(); iPos++)
			{
				if (!MixedNucleotideByte::AreConsistent(m1[iPos], m2[iPos]))
				{
					return false;
				}
			}
			return true;
		}

	public:
		double MutationRate;

		static std::unordered_map<std::vector<byte>, double> MotifMutabilities;

		GeneReplicator(fstream& stream)
		{
			readMotifMutability(stream);
			InitializeMutabilities(1);
		}

		std::vector<byte> Motif(std::vector<byte>& subSeq)
		{
			for (auto& m : motifs)
			{
				if (matches(m, subSeq))
				{
					return m;
				}
			}
			return std::vector<byte>(subSeq.size());
		}

	private:
		void readMotifMutability(fstream& stream)
		{
			throw "not implemented exception";
		}
	public:

		static void InitializeMutabilities(double exponent)
		{

			throw "not implemented exception";

			//MotifMutabilities = new Dictionary<byte[], double>(new ByteArrayEqualityComparer());
			//foreach(KeyValuePair<byte[], double> kvp in baseMutabilities)
			//{
			//	if (kvp.Value > 0)
			//	{
			//		MotifMutabilities.Add(kvp.Key, Math.Pow(kvp.Value, exponent));
			//	}
			//	else
			//	{
			//		MotifMutabilities.Add(kvp.Key, 0);
			//	}
			//}
		}

		std::vector<byte> Mutate(std::vector<byte>& seq)
		{
			// note that the first few (number = half-width) and last few nucleotides
			// will not be mutated

			throw "not implemented exception";

			/*byte[] mutatedSeq = (byte[])seq.Clone();
			for (int pos = halfWidth; pos < seq.Length - halfWidth; pos++)
			{
				byte[] subSeq = new byte[width];
				Array.Copy(seq, pos - halfWidth, subSeq, 0, width);
				byte[] mo = Motif(subSeq);

				bern.Alpha = MutationRate * MotifMutabilities[mo];
				if (bern.Next() == 1)
				{
					mutatedSeq[pos] = newNucleotide[seq[pos]][SpectrumProbabilities[mo].Next()];
				}
			}
			return mutatedSeq;*/
		}

		static byte Mutate(byte nucleotide, SubstitutionType sType)
		{
			throw "not implemented exception";

			//byte mutated;
			//switch (sType)
			//{
			//case SubstitutionType.Transition:
			//	mutated = newNucleotide[nucleotide][0];
			//	break;
			//case SubstitutionType.TransversionToComplement:
			//	mutated = mutated = newNucleotide[nucleotide][1];
			//	break;
			//case SubstitutionType.TransversionToNonComplement:
			//	mutated = mutated = newNucleotide[nucleotide][2];
			//	break;
			//default:
			//	mutated = nucleotide;
			//	break;
			//}
			//return mutated;
		}

		void NormalizeMutationRate(double targetRate, std::vector<byte>& seq)
		{

			throw "not implemented exception";

			//double totalRate = 0;
			//int count = 0;
			//for (int pos = halfWidth; pos < seq.Length - halfWidth; pos++)
			//{
			//	byte[] subSeq = new byte[width];
			//	Array.Copy(seq, pos, subSeq, 0, width);
			//	byte[] mo = Motif(subSeq);
			//	totalRate += MotifMutabilities[mo];
			//	count++;
			//}
			//MutationRate = count * targetRate / totalRate;
		}

		double GetMutationRate(std::vector<byte>& motif)
		{
			throw "not implemented exception";

			//if (motif.Length != width)
			//	return -1;

			//byte[] mo = Motif(motif);
			//if (MotifMutabilities.ContainsKey(mo))
			//{
			//	return MotifMutabilities[mo];
			//}
			//else
			//{
			//	return -1;
			//}
		}

		std::vector<double> GetMutationRates(Polymer<byte>& p)
		{
			throw "not implemented exception";

			/*double[] mutationRates = new double[p.Seq.Length];
			for (int pos = 0; pos < halfWidth; pos++)
			{
				mutationRates[pos] = 1.0;
			}

			for (int pos = halfWidth; pos < p.Seq.Length - halfWidth; pos++)
			{
				byte[] subSeq = new byte[width];
				Array.Copy(p.Seq, pos - halfWidth, subSeq, 0, width);
				byte[] mo = Motif(subSeq);

				mutationRates[pos] = MotifMutabilities[mo];
			}

			for (int pos = p.Seq.Length - halfWidth; pos < p.Seq.Length; pos++)
			{
				mutationRates[pos] = 1;
			}

			return mutationRates;*/
		}
	private:
		static unordered_map<byte, std::vector<byte>> newNucleotide;
		//
		//= new Dictionary<byte, byte[]>()
		//{
		//	{1, new byte[]{ 2,4,8 }},
		//	{ 2,new byte[]{ 1,8,4 } },
		//	{ 4,new byte[]{ 8,1,2 } },
		//	{ 8,new byte[]{ 4,2,1 } }
		//};
	};
}

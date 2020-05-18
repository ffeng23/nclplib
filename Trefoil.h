#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include "Polymer.h"
#include "Annotations.h"
#include "PolymerPair.h"
#include "EvolutionModels.h"
#include "AlignmentKernel.h"
#include "PairwiseAligner.h"
#include "Writer.h"

namespace CLPLib
{

	class Trefoil
	{
	private:
		std::vector<std::vector<int>> index;
		//shared_ptr<Polymer<std::vector<double>>> MRCANormalized;
		//std::vector<Polymer<std::vector<double>>> polymers;

	public:
		shared_ptr<Polymer<std::vector<double>>> RootLikelihoodSequence;

		shared_ptr<Polymer<std::vector<double>>> MRCALikelihoodSequence;


		//private:
		//std::vector<double> distancesToMRCA;
		//double distanceToRoot;

		//public:
		int AlignmentLength()
		{
			return (int)index.size();
		}

		std::vector<std::vector<int>> SegmentStart;

		std::vector<int> SegmentLength;

		int NSegments;

		std::vector<std::vector<int>> Index;

	private:
		double nearlyOne;

	public:
		std::vector<std::shared_ptr<Polymer<std::vector<double>>>> Polymers;

		//todo: if we don't want to declare this as pointer, then we will
		//need to add a default constructor for polymer
		shared_ptr<Polymer<std::vector<double>>> RootNormalized;

		shared_ptr<Polymer<std::vector<double>>> MRCANormalized = MRCANormalized;

		std::vector<double> DistancesToMRCA;

		double DistanceToRoot;

		double IndelPartialLikelihood;

		Trefoil(PolymerPair<std::vector<double>, std::vector<double>>& p1, PolymerPair<std::vector<double>, std::vector<double>>& p2)
		{
			nearlyOne = 1 - 0.0001;
			// merge alignment indices
			if (p1.Polymer1->Seq.size() == p2.Polymer1->Seq.size())
			{
				index = MergeIndices(p1.Index, p2.Index);
			}
			else
			{
				index = AlignIndices(p1, p2);
			}

			// create the segment list
			AlignmentToSegments();

			//polymers = new List<Polymer<double[]>>();
			Polymers.push_back(p1.Polymer0);
			Polymers.push_back(p2.Polymer0);
		}

		static std::vector<std::vector<int>> MergeIndices(std::vector<std::vector<int>>& index1, std::vector<std::vector<int>>& index2)
		{
			std::vector<std::vector<int>> index;
			int k1 = 0;
			int k2 = 0;

			// make sure the common reference index points to the same position in both reference sequences
			while (index1[k1][1] != index2[k2][1])
			{
				if (k1 >= (int)index1.size() || k2 >= (int)index2.size())
				{
					//return null
					return std::vector<std::vector<int>>();
				}

				if (index1[k1][1] < index2[k2][1])
				{
					k1++;
				}
				else
				{

					k2++;
				}
			}

			int i1, i2, ia, ir; // entries for the index array
			int ca = index1[k1][1]; // coordinate for the ancestor. Initial value is the common coordinate for the two references
			while (k1 < (int)index1.size() && k2 < (int)index2.size())
			{
				if (index1[k1][1] > -1 && index2[k2][1] > -1)
				{
					// both references are ungapped here
					i1 = index1[k1][0];
					i2 = index2[k2][0];
					ir = index1[k1][1];
					k1++;
					k2++;
					if (i1 > -1 || i2 > -1)
					{
						// at least one of the observed sequences is ungapped here. Let the ancestor be ungapped, too.
						ia = ca;
						ca++;
					}
					else
					{
						ia = -1;
					}
				}
				else if (index1[k1][1] > -1) // the first reference sequence is not gapped, but the second is. Thus, the second observed sequence has an insertion
											 // but the other does not. Gap the mrca, the root, the first observed sequence.
				{
					i1 = -1;
					ir = -1;
					ia = -1;
					i2 = index2[k2][0];
					k2++;
				}
				else if (index2[k2][1] > -1) // vice-versa
				{
					i1 = index1[k1][0];
					i2 = -1;
					ia = -1;
					ir = -1;
					k1++;
				}
				else // both reference sequences have gaps, both observed sequences have insertions. 
				{
					i1 = index1[k1][0];
					i2 = index2[k2][0];
					ia = ca;
					ir = -1;
					ca++;
					k1++;
					k2++;
				}
				index.push_back(std::vector<int> { i1, i2, ia, ir });
			}
			return index;
		}

		static std::vector<std::vector<int>> AlignIndices(PolymerPair<std::vector<double>, std::vector<double>>& p1, PolymerPair<std::vector<double>, std::vector<double>>& p2)
		{
			double time = 0.2;
			unique_ptr<JCEvolutionModel> model(new JCEvolutionModel(time));
			Polymer<byte> p1b = TypeConverter::ToPolynucleotideByte(*p1.Polymer1);
			Polymer<byte> p2b = TypeConverter::ToPolynucleotideByte(*p2.Polymer1);

			//MNMNAlignmentKernelSymmetric kernel(std::move(model), std::make_shared<Polymer<byte>>(p1b), std::make_shared<Polymer<byte>>(p2b));
			std::unique_ptr<MNMNAlignmentKernelSymmetric> kernel(new MNMNAlignmentKernelSymmetric(std::move(model), std::make_shared<Polymer<byte>>(p1b), std::make_shared<Polymer<byte>>(p2b)));

			PairwiseAlignerSymmetric<byte> aligner(std::move(kernel));
			std::vector<int> start = aligner.FillScoreMatrix();
			PolymerPair<byte, byte> pair = aligner.TraceBack(start);

			// debug
			bool debug_flag = false;
			if (debug_flag)
			{
				string outputstr = Writer<byte>::WriteFasta(std::make_shared<MixedNucleotideByte>(), pair);
				std::ofstream fs("alignIndicesPair.fasta");
				if (!fs)
				{
					std::cerr << "Cannot open the output file: alignIndicesPair.fasta" << std::endl;
					exit(1);
				}
				fs << outputstr;
				fs.close();
			}

			std::vector<std::vector<int>> index;
			std::vector<std::vector<int>>& index0 = p1.Index;
			std::vector<std::vector<int>>& index1 = p2.Index;

			int k0 = 0;
			int k1 = 0;
			int kR = 0;
			int ka = 0;
			int k = 0;

			int i0 = 0;
			int i1 = 0;
			int iR = 0;
			int ia = 0;
			while (k0 < (int)index0.size() && k1 < (int)index1.size() && kR < (int)pair.Index.size())
			{
				k++;
				if (pair.Index[kR][1] > -1 && pair.Index[kR][0] > -1) // note that these indices cannot both be gaps
				{
					if (index0[k0][1] > -1)
					{
						if (index1[k1][1] > -1) // A: both references are ungapped
						{
							i0 = index0[k0][0];
							k0++;
							i1 = index1[k1][0];
							k1++;
							iR = kR;
							kR++;
							ia = ka;
							ka++;
						}
						else // B1: reference 0 is ungapped, reference 1 is gapped;
						{
							// only i1 advances
							i0 = -1;
							i1 = index1[k1][0];
							k1++;
							iR = -1;
							ia = -1;
						}
					}
					else // reference 0 is gapped
					{
						if (index1[k1][1] > -1) // B0: reference 1 is ungapped, reference 0 is gapped
						{
							i0 = index0[k0][0];
							k0++;
							i1 = -1;
							iR = -1;
							ia = -1;
						}
						else // C: both references are gapped
						{
							i0 = index0[k0][0];
							k0++;
							i1 = index1[k1][0];
							k1++;
							iR = -1;
							ia = ka;
							ka++;
						}
					}
				}
				else if (pair.Index[kR][0] < 0 && pair.Index[kR][1] > -1)/* one of the reference sequences has a gap in their mutual alignment
																		 One must choose the reference sequence appropriately in this case.  */
				{
					if (index0[k0][1] > -1)
					{
						if (index1[k1][1] > -1) // D0
						{
							i1 = index1[k1][0];
							k1++;
							i0 = -1;
							iR = kR;
							kR++;
							ia = ka;
							ka++;
						}
						else // E10
						{
							i0 = -1;
							i1 = -1;
							iR = kR;
							kR++;
							ia = -1;
						}
					}
					else // index0[k0][1] < 0
					{
						if (index1[k1][1] > -1) // E00
						{
							i0 = index0[k0][0];
							k0++;
							i1 = -1;
							iR = kR;
							kR++;
							ia = ka;
							ka++;
						}
						else // F0
						{
							i0 = index0[k0][0];
							k0++;
							i1 = index1[k1][0];
							k1++;
							iR = kR;
							kR++;
							ia = ka;
							ka++;
						}
					}
				}
				else // pair.Index[k][0] > -1 & pair.Index[k][1] < 0
				{
					if (index0[k0][1] > -1)
					{
						if (index1[k1][1] > -1) // D1
						{
							i0 = index0[k0][1];
							k0++;
							i1 = -1;
							iR = kR;
							kR++;
							ia = ka;
							ka++;
						}
						else // E11
						{
							i0 = -1;
							i1 = index1[k1][0];
							k1++;
							iR = kR;
							kR++;
							ia = ka;
							ka++;
						}
					}
					else // index0[k0][1] < 0
					{
						if (index1[k1][1] > -1) // E01
						{
							i0 = -1;
							i1 = -1;
							iR = kR;
							kR++;
							ia = -1;
						}
						else // F1
						{
							i0 = index0[k0][0];
							k0++;
							i1 = index1[k1][0];
							k1++;
							iR = kR;
							kR++;
							ia = ka;
							ka++;
						}
					}
				}
				index.push_back(std::vector<int>{ i0, i1, ia, iR });
			}
			return index;

		}

		void TransferNucleotides(NucleotideEvolutionModel& model, int cdr3Start = 0)
		{
			double sum = 0;
			for (int k = cdr3Start; k < (int)index.size(); k++)
			{
				if (index[k][3] > -1 && index[k][2] > -1)
				{
					auto tmp = ((RootNormalized->Seq))[index[k][3]];
					double sum = 0;
					for (auto &v : tmp)
					{
						sum += v;
					}
					if (sum < nearlyOne)
					{
						vector<double> &src = ((MRCANormalized->Seq))[index[k][2]];
						std::vector<double> pmf = model.DestinationPMF(src, (int)this->DistanceToRoot);
						double weight = 1 - sum;
						for (int j = 0; j < 5; j++)
						{
							((RootNormalized->Seq))[index[k][3]][j] += weight * pmf[j];
						}
					}
				}
			}
		}

		void AlignmentToSegments()
		{
			std::vector<std::vector<int>> segmentStartAsList;
			segmentStartAsList.push_back(std::vector<int>(4, 0));
			std::vector<int> segmentLengthAsList;
			segmentLengthAsList.push_back(0);

			int segment = 0;
			for (int i = 0; i < 4; i++)
			{
				if (Index[0][i] > -1)
				{
					// initial segment is present in polymer i
					segmentStartAsList[0][i] = 0;
				}
				else
				{
					segmentStartAsList[0][i] = -1;
				}
			}

			for (int k = 1; k < (int)Index.size() - 1; k++)
			{
				if (isBreakpoint(k))
				{
					segment++;
					segmentStartAsList.push_back(std::vector<int>(4, 0));
					segmentLengthAsList.push_back(0);

					for (int i = 0; i < 4; i++)
					{
						if (Index[k + 1][i] < 0) // the new segment is absent in sequence 0...
						{
							segmentStartAsList[segment][i] = -1;
						}
						else
						{
							segmentStartAsList[segment][i] = Index[k + 1][i];
						}
					}

					// compute length
					for (int i = 0; i < 4; i++)
					{
						if (segmentStartAsList[segment - 1][i] > -1)
						{
							segmentLengthAsList[segment - 1] = Index[k][i] + 1 - segmentStartAsList[segment - 1][i];
							break;
						}
					}
				}
			}

			// finish up
			NSegments = segment + 1;

			for (int i = 0; i < 4; i++)
			{
				if (segmentStartAsList[segment][i] > -1)
				{
					segmentLengthAsList[segment] = Index.back()[i] + 1 - segmentStartAsList[segment][i];
					break;
				}
			}

			//??
			SegmentLength = segmentLengthAsList;
			SegmentStart = segmentStartAsList;
		}

		bool isBreakpoint(int k)
		{
			for (int i = 0; i < 4; i++)
			{
				if ((Index[k][i] < 0 && Index[k + 1][i] > -1) || (Index[k][i] > -1 && Index[k + 1][i] < 0))
					return true;
			}
			return false;
		}
	};

}

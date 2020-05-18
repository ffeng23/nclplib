#pragma once

#include "Annotations.h"
#include "Polymer.h"

namespace CLPLib
{

	/// <summary>
	/// A pair of polynucleotides that may be of mixed type. This class provides the basis for pairwise alignment.
	/// </summary>
	/// <typeparam name="T">The nucleotide type of the first polynucleotide.</typeparam>
	/// <typeparam name="U">The nucleotide type of the second polynucleotide.</typeparam>
	template<class T, class U>
	class PolymerPair
	{
	private:
		const double logScoreConst = log(0.25);
	public:
		/// <summary>
		/// The first of two polynucleotides.
		/// </summary>
		shared_ptr<Polymer<T>> Polymer0;
		/// <summary>
		/// The second of two polynucleotides.
		/// </summary>
		shared_ptr<Polymer<U>> Polymer1;


		///// <summary>
		///// The total length of the pairwise alignment between the two polynucleotides.
		///// </summary>
		//int Length;
		/// <summary>
		/// The alignment score of the pairwise alignment.
		/// </summary>
		double AlignmentScore;

		//the rawalignment score is not the final score.
		//double AdjustScore()
		//{

		//	int alignLen = (int)this->Index.size();
		//	int sumLen = (int)Polymer0->Seq.size() + (int)Polymer1->Seq.size();

		//	auto sumscore = -sumLen * logScoreConst;
		//	sumscore += AlignmentScore;

		//	if (sumscore > 0)
		//	{
		//		//enforcing at least 15 bp match??
		//		if (sumLen - alignLen < 15)return 0;
		//	}
		//	return sumscore;
		//}


		/// <summary>
		/// The array of index values. The first index is the position in the alignment, 
		/// the second index determines the polynucleotide to which the index refers (0 or 1).
		/// </summary>
		vector<vector<int>> Index;
		/// <summary>
		/// The alignment score by position.
		/// </summary>
		vector<double> Score;
		

		int NSegments;

		vector<vector<int>> SegmentStart;

		vector<int> SegmentLength;

		/// <summary>
		/// Constructor. Sets AlignmentScore to NaN.
		/// </summary>
		PolymerPair()
		{
			AlignmentScore = std::numeric_limits<double>::quiet_NaN();
		}

		/// <summary>
		/// Constructor. Creates a pair with specified polynucleotides.
		/// </summary>
		/// <param name="p0">The first polynucleotide in the pair.</param>
		/// <param name="p1">The second polynucleotide in the pair.</param>
		PolymerPair(shared_ptr<Polymer<T>> p1, shared_ptr<Polymer<U>> p2)
			:Polymer0(p1), Polymer1(p2), AlignmentScore(std::numeric_limits<double>::quiet_NaN())
		{
		}

		/// <summary>
		/// Creates a PolynucleotidePair from two Polynucleotides under the assumption that they are already aligned.
		/// </summary>
		/// <param name="p0">The first polynucleotide.</param>
		/// <param name="p1">The second polynucleotide.</param>
		/// <returns>A PolynucleotidePair with the specified polynucleotides and a completed Index.</returns>
		static PolymerPair<T, U> PrealignedPair(shared_ptr<Polymer<T>> p1, shared_ptr<Polymer<U>> p2)
		{
			PolymerPair<T, U> pair(p1, p2);
			pair.Index = vector<int>(p1->Seq.length(), vector<int>(2));;
			for (int i = 0; i < pair.Index.size(); i++)
			{
				pair.Index[i] = vector<int>(2, i);
			}
			return pair;
		}

		static vector<vector<int>> MergeIndices(vector<vector<int>> &index1, vector<vector<int>>& index2)
		{

			vector<vector<int>> index;

			int k1 = 0;
			int k2 = 0;
			while (index1[k1][1] != index2[k2][1])
			{
				if (k1 >= index1.size() || k2 >= index2.size())
				{
					throw logic_error("index out of range");
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

			while (k1 < index1.size() && k2 < index2.size())
			{
				if (index1[k1][1] == index2[k2][1])
				{

					index.push_back(vector<int>{ index1[k1][0], index2[k2][0], index1[k1][1] });
					k1++;
					k2++;
				}
				else if (index1[k1][1] < 0)
				{
					index.push_back(vector<int>{index1[k1][0], -1, index1[k1][1]});
					k1++;
				}
				else
				{
					index.push_back(vector<int>{ -1, index2[k2][0], index2[k2][1] });
				}
			}
			return index;
		}

		static void TransferGaps(PolymerPair<T, U>& pair, MonomerEncoding<T>& encoding1, MonomerEncoding<U>& encoding2)
		{
			vector<T> newSeq1(pair.Index.size());
			vector<U> newSeq2(pair.Index.size());

			for (int k = 0; k < pair.Index.size(); k++)
			{
				if (pair.Index[k][0] > -1)
				{
					newSeq1[k] = (*pair.Polymer0->Seq)[pair.Index[k][0]];
				}
				else
				{
					newSeq1[k] = encoding1.Gap();
				}
				pair.Index[k][0] = k;

				if (pair.Index[k][1] > -1)
				{
					newSeq2[k] = (*pair.Polymer1->Seq)[pair.Index[k][1]];
				}
				else
				{
					newSeq2[k] = encoding2.Gap();
				}
				pair.Index[k][1] = k;
			}
			*pair.Polymer0->Seq = std::move(newSeq1);
			*pair.Polymer1->Seq = std::move(newSeq2);
		}

		void AlignmentToSegments()
		{
			vector<vector<int>> segmentStartAsList;
			segmentStartAsList.push_back(vector<int>(2, 0));
			//to do correct below
			std::vector<int> segmentLengthAsList;
			segmentLengthAsList.push_back(0);

			int segment = 0;
			if (Index[0][0] > -1)
			{
				// initial segment is present in polymer 0
				segmentStartAsList[0][0] = 0;
			}
			else
			{
				segmentStartAsList[0][0] = -1;
			}
			if (Index[0][1] > -1)
			{
				// initial segment is present in polymer 1
				segmentStartAsList[0][1] = 0;
			}
			else
			{
				segmentStartAsList[0][1] = -1;
			}

			for (int k = 1; k < Index.size() - 1; k++)
			{
				if ((Index[k][0] < 0 && Index[k + 1][0] > -1) || (Index[k][1] < 0 && Index[k + 1][1] > -1)
					|| (Index[k][0] > -1 && Index[k + 1][0] < 0) || (Index[k][1] > -1 && Index[k + 1][1] <0))
				{
					segment++;
					segmentStartAsList.push_back(vector<int>(2, 0));
					segmentLengthAsList.push_back(0);

					if (Index[k + 1][0] < 0) // the new segment is absent in sequence 0...
					{
						segmentStartAsList[segment][0] = -1;
					}
					else
					{
						segmentStartAsList[segment][0] = Index[k + 1][0];
					}

					if (Index[k + 1][1] < 0) // the new segment is absent in sequence 1
					{
						segmentStartAsList[segment][1] = -1;
					}
					else
					{
						segmentStartAsList[segment][1] = Index[k + 1][1];
					}

					if (segmentStartAsList[segment - 1][0] > -1)
					{
						segmentLengthAsList[segment - 1] = Index[k][0] + 1 - segmentStartAsList[segment - 1][0];
					}
					else
					{
						segmentLengthAsList[segment - 1] = Index[k][1] + 1 - segmentStartAsList[segment - 1][1];
					}
				}
			}
			// finish up
			NSegments = segment + 1;
			if (segmentStartAsList[segment][0] > -1)
			{
				segmentLengthAsList[segment] = Index[Index.size() - 1][0] + 1 - segmentStartAsList[segment][0];
			}
			else
			{
				segmentLengthAsList[segment] = Index[Index.size() - 1][1] + 1 - segmentStartAsList[segment][1];
			}
			SegmentLength = segmentLengthAsList;
			SegmentStart = segmentStartAsList;
		}

	};

}






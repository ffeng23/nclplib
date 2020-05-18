#pragma once

//#define thread_local _Thread_local

#include <vector>
#include <memory>
#include <thread>
#include <stack>
#include "AlignmentKernel.h"
#include "alignedAllocator.h"

#define MATRIX_INDEX(i, j, k) ((i) * dim2 * 3 + (j) * 3 + (k))


#define MATRIX_INDEX_4(i, j, k) ((i) * dim2 * 4 + (j) * 4 + (k))




using namespace std;

extern int AlignMinMatch;
extern int AlignMinScore;

extern int DEBUG_CODE;

namespace CLPLib
{

	typedef struct ScoreEntry
	{
		double insertScore;
		double deleteScore;
		double subScore;
		double maxScore;
	}ScoreEntryStr;


	template <class T>
	T FastMax(const T& left, const T& right)
	{
		return left > right ? left : right;
	}

	template <class T, class U>
	class PairwiseAlignerAsymmetric
	{
	public:
		/// <summary>
		/// The total score for the alignment as reported.
		/// </summary>
		double AlignmentScore;

	private:
		std::vector<double> scoreMatrix;
		vector<int> traceBackMatrix;

		vector<byte> TraceBackBytes;

		ScoreEntryStr *ScoreStrMatrix = nullptr;
		int ScoreStrMatrixSize = 0;

		vector<double> cumulativeScoreQuery;
		vector<double> cumulativeScore1;

		int dimQueryPhysical = -1;
		int dimReferencePhysical = -1;
		int dimQuery;
		int dimReference;


		//static _Thread_local unique_ptr<double> CompScoreStorage;
		static unique_ptr<double[]> CompScoreStorage;

		shared_ptr<AlignmentKernelAsymmetric<T, U>> kernel;

		bool use_index_version = false;

	public:

		int MaxRowIndex;
		int MaxColIndex;

		/// <summary>
		/// Sets the default scores for the first column.
		/// </summary>
		/// <param name="values">A vector of cumulative alignment scores taken as the alignment scores that occur 
		///  if the polynucleotide local alignment begins on a row other than the zeroth.</param>
		void SetFirstColumn(vector<double>& values)
		{
			cumulativeScoreQuery = values;
		}

		/// <summary>
		/// Sets the default scores for the first column.
		/// </summary>
		/// <param name="value">The score to be accumulated over rows as the alignment scores that occur if the polynucleotide local alignment begins on a row other than the zeroth.</param>
		void SetFirstColumn()
		{
			if (dimQuery < 1)
				return;

			cumulativeScoreQuery[0] = kernel->PriorQueryScore(0);
			for (int i = 1; i < dimQuery; i++)
			{
				cumulativeScoreQuery[i] = cumulativeScoreQuery[i - 1] + kernel->PriorQueryScore(i);
			}
		}

		/// <summary>
		/// Sets the default scores for the first column.
		/// </summary>
		/// <param name="value">The score to be accumulated over rows as the alignment scores that occur if the polynucleotide local alignment begins on a row other than the zeroth.</param>
		void SetFirstColumn(double defaultValue)
		{
			if (dimQuery < 1)
				return;

			cumulativeScoreQuery[0] = defaultValue;
			for (int i = 1; i < dimQuery; i++)
			{
				cumulativeScoreQuery[i] = cumulativeScoreQuery[i - 1] + defaultValue;
			}
		}

		/// <summary>
		/// Sets the default scores for the first row.
		/// </summary>
		/// <param name="values">A vector of cumulative alignment scores taken as the alignment scores that occur if the polynucleotide local alignment begins on a column other than the zeroth.</param>
		void SetFirstRow(vector<double>& values)
		{
			cumulativeScore1 = values;
		}

		/// <summary>
		/// Sets the default scores for the first row.
		/// </summary>
		/// <param name="value">The score to be accumulated over rows as the alignment scores that occur if the polynucleotide local alignment begins on a column other than the zeroth.</param>
		void SetFirstRow(double value = 0)
		{
			if (dimReference < 1)
				return;

			cumulativeScore1[0] = value;
			for (int i = 1; i < dimReference; i++)
			{
				cumulativeScore1[i] = cumulativeScore1[i - 1] + value;
			}
		}

		/// <summary>
		/// Constructor.
		/// </summary>
		/// <param name="_pair">The polynucleotide pair to be aligned.</param>
		/// <param name="_comparator">A comparator of the appropriate types for use in the alignment.</param>
		PairwiseAlignerAsymmetric(unique_ptr<AlignmentKernelAsymmetric<T, U>> _kernel)
			:kernel(std::move(_kernel))
		{
			if (!kernel->GetQuery())
			{
				dimQuery = dimQueryPhysical = 0;
			}
			else
			{
				dimQuery = dimQueryPhysical = kernel->GetQuery()->Seq.size();
			}

			if (!kernel->GetReference())
			{
				dimReference = dimReferencePhysical = 0;
			}
			else
			{
				dimReference = dimReferencePhysical = kernel->GetReference()->Seq.size();
			}

			// allocate matrices
			scoreMatrix.resize(dimQueryPhysical * dimReferencePhysical * 3);
			traceBackMatrix.resize(dimQueryPhysical * dimReferencePhysical * 3);
			cumulativeScoreQuery.resize(dimQueryPhysical, 0);
			cumulativeScore1.resize(dimReferencePhysical, 0);

			SetFirstColumn();
			SetFirstRow();
		}

		/// <summary>
		/// Constructor.
		/// </summary>
		/// <param name="_kernel">A comparator of the appropriate types for use in the alignment.</param>
		/// <param name="_dim0">The initial first dimension of the score and traceback matrices.</param>
		/// <param name="_dim1">The initial second dimension of the score and traceback matrices.</param>
		PairwiseAlignerAsymmetric(int _dim0, int _dim1, shared_ptr<AlignmentKernelAsymmetric<T, U>> _kernel)
			:kernel(_kernel)
		{
			dimQuery = dimQueryPhysical = _dim0;
			dimReference = dimReferencePhysical = _dim1;

			// allocate matrices
			scoreMatrix.resize(dimQueryPhysical * dimReferencePhysical * 3);
			traceBackMatrix.resize(dimQueryPhysical * dimReferencePhysical * 3);
			cumulativeScoreQuery.resize(dimQueryPhysical);
			cumulativeScore1.resize(dimReferencePhysical);
		}

		~PairwiseAlignerAsymmetric()
		{

			if (ScoreStrMatrix != nullptr)
			{
				free(ScoreStrMatrix);
			}
		}

	private:

		friend class TestHelper;

		void allocateMatrices()
		{
			scoreMatrix.resize(dimQueryPhysical * dimReferencePhysical * 3);
			traceBackMatrix.resize(dimQueryPhysical * dimReferencePhysical * 3, 0);

			cumulativeScoreQuery.resize(dimQueryPhysical);
			cumulativeScore1.resize(dimReferencePhysical);
		}


	public:

		void SetQuery(shared_ptr<Polymer<T>> p, vector<double> defaultScores = vector<double>())
		{
			kernel->SetQuery(p);
			dimQuery = (int)p->Seq.size() + 1;
			if (dimQuery > dimQueryPhysical)
			{
				dimQueryPhysical = dimQuery;
				//scoreMatrix.resize(dimQueryPhysical * dimReferencePhysical * 3);
				//traceBackMatrix.resize(dimQueryPhysical * dimReferencePhysical * 3);
				cumulativeScoreQuery.resize(dimQueryPhysical);
			}

			if (defaultScores.size() > 0)
			{
				SetFirstColumn(defaultScores);
			}
			else
			{
				SetFirstColumn();
			}
		}

		void SetReference(shared_ptr<Polymer<U>> p, vector<double> defaultScores = vector<double>())
		{
			kernel->SetReference(p);
			dimReference = (int)p->Seq.size();;
			if (dimReference > dimReferencePhysical)
			{
				dimReferencePhysical = dimReference;
				//scoreMatrix.resize(dimQueryPhysical * dimReferencePhysical * 3);
				//traceBackMatrix.resize(dimQueryPhysical * dimReferencePhysical * 3);
				cumulativeScore1.resize(dimReferencePhysical);
			}

			if (defaultScores.size() > 0)
			{
				SetFirstRow(defaultScores);
			}
			else
			{
				SetFirstRow();
			}
		}

		double ScoreAlignedPair(AlignmentKernelAsymmetric<T, U>& _kernel, PolymerPair<T, U>& _pair)
		{
			// This method assumes that a query-reference pair has already been aligned and that a second reference, sufficiently similar to the first that a complete 
			// realignment is not necessary, is being substituted for the first.
			//if (cumulativeScoreQuery != nullptr)
			//    SetFirstColumn();

			// for this method, we assume that _kernel->query = _pair.Polynucleotide1
			if (_pair.Index.size() == 0)
				return std::numeric_limits<double>::lowest();

			int offset = std::min(_pair.Index[0][0], _pair.Index[0][1]);

			// firstx is the first position (in seqx coordinates) in the aligned _kernel->Pair where both sequences have nucleotides
			int first0 = _pair.Index[0][0] - offset;
			int first1 = _pair.Index[0][1] - offset;
			vector<int>& lastIndexed = _pair.Index.back();

			// the overhang is computed because the new allele may align to the query beyond the point where the original allele aligned.
			int overhang = std::min((int)_pair.Polymer0->Seq.size() - lastIndexed[0], (int)_pair.Polymer1->Seq.size() - lastIndexed[1]) - 1;

			// analogously for lastx
			int last0 = lastIndexed[0] + overhang;
			int last1 = lastIndexed[1] + overhang;

			double unmatchedScore = 0;
			if (first0 > 0)
			{
				unmatchedScore = cumulativeScoreQuery[first0 - 1];
			}

			int newLength = (int)_pair.Index.size() + offset + overhang;

			vector<double> newScore(newLength, 0);
			vector<int> traceback(newLength, 0);
			vector<vector<int>> newIndex(newLength);

			newScore[0] = unmatchedScore + kernel->ComparisonScore(first0, first1);
			traceback[0] = 1;
			newIndex[0] = vector<int>{ first0, first1 };
			// compute scores up to and including the start of the existing index  
			for (int kNew = 1; kNew < offset + 1; kNew++)
			{
				int i0 = first0 + kNew;
				int i1 = first1 + kNew;
				newScore[kNew] = kernel->ComparisonScore(i0, i1);
				if (newScore[kNew - 1] > cumulativeScoreQuery[i0 - 1])
				{
					newScore[kNew] += newScore[kNew - 1];
					traceback[kNew] = 0;
				}
				else
				{
					newScore[kNew] += cumulativeScoreQuery[i0 - 1];
					traceback[kNew] = 1;
				}
				newIndex[kNew] = vector<int>{ i0, i1 };
			}

			// compute scores through the end of the existing index
			int endK = (int)_pair.Index.size() + std::min(0, overhang);
			for (int k = 1; k < endK; k++)
			{
				int i1 = _pair.Index[k][1];
				if (i1 == _pair.Polymer1->Seq.size()) // if there is a gap in the original allele between the end of the original index 
				{                                   // and the end of the new allele, this condition will get triggered
					newLength -= (endK - k);
					break;
				}

				int kNew = offset + k;
				int i0 = _pair.Index[k][0];
				if (i0 > -1 && i1 > -1) // both sequences have a nucleotide at this position
				{
					newScore[kNew] = kernel->ComparisonScore(i0, i1);
					if (newScore[kNew - 1] > cumulativeScoreQuery[i0 - 1])
					{
						newScore[kNew] += newScore[kNew - 1];
						traceback[kNew] = 0;
					}
					else
					{
						newScore[kNew] += cumulativeScoreQuery[i0 - 1];
						traceback[kNew] = 1;
					}
				}
				else if (_pair.Index[k][0] < 0) // Sequence0 has a gap at this position
				{
					// if this gap is aligned opposite a gap in sequence1, no penalty
					if (kernel->IsGap((*_pair.Polymer1->Seq)[i1]))
					{
						newScore[kNew] = 0;
					}

					// note that the first index value cannot be negative for either sequence
					if (_pair.Index[k - 1][0] < 0)
					{
						newScore[kNew] = kernel->ContinueDeletionScore(i0, i1);
					}
					else
					{
						newScore[kNew] = kernel->OpenDeletionScore(i0, i1);
					}
					newScore[kNew] += newScore[kNew - 1];
					traceback[kNew] = 0;
				}
				else // sequence1 has a gap at this position.
				{
					if (_pair.Index[k - 1][1] < 0)
					{
						newScore[kNew] = kernel->ContinueInsertionScore(i0, i1);
					}
					else
					{
						newScore[kNew] = kernel->OpenInsertionScore(i0, i1);
					}
					newScore[kNew] += newScore[kNew - 1];
					traceback[kNew] = 0;
				}
				newIndex[kNew] = _pair.Index[k];
			}

			// continue past the end of the old index
			for (int deltaKNew = 0; deltaKNew < overhang; deltaKNew++)
			{
				int kNew = (int)_pair.Index.size() + offset + deltaKNew;
				int i0 = lastIndexed[0] + deltaKNew + 1;
				int i1 = lastIndexed[1] + deltaKNew + 1;
				newScore[kNew] = newScore[kNew - 1] + kernel->ComparisonScore(i0, i1);
				traceback[kNew] = 0;
				newIndex[kNew] = vector<int>{ i0, i1 };
			}

			// do traceback
			int stop = newLength - 1;
			double cumScoreLast = cumulativeScoreQuery[_pair.Polymer0->Seq.size() - 1];
			double totalScore = newScore[stop] + cumScoreLast - cumulativeScoreQuery[last0];
			double maxTotalScore = totalScore;
			// we have to go back through the index here
			for (int kNew = newLength - 1; kNew > 0; kNew--)
			{
				int i0 = newIndex[kNew][0];
				if (i0 < 0)
					continue;

				totalScore = cumScoreLast - cumulativeScoreQuery[i0] + newScore[kNew];
				if (totalScore > maxTotalScore)
				{
					stop = kNew;
					maxTotalScore = totalScore;
				}
			}
			int start = stop;
			while (traceback[start] == 0)
			{
				start--;
			}

			// form new _kernel->Pair index and new _kernel->Pair score array
			_pair.Index = vector<vector<int>>(newIndex.begin(), newIndex.begin() + stop - start + 1);
			_pair.Score = vector<double>(newScore.begin(), newScore.begin() + stop - start + 1);
			_pair.AlignmentScore = maxTotalScore;
			return maxTotalScore;
		}


		//template <class T>
		//T FastMax(const T& left, const T& right)
		//{
		//	return left > right ? left : right;
		//}

		///supposed to be same as Tom's currurnet version. optimized
		double FillScoreMatrix_get_score()
		{
			auto query = kernel->GetQuery();
			auto reference = kernel->GetReference();

			bool isPure = query->isPureNucleotide() && reference->isPureNucleotide();

			auto& querySeq = isPure ? query->PureSeq : query->Seq;
			auto& refSeq = isPure ? reference->PureSeq : reference->Seq;

			dimQuery = (int)querySeq.size();
			dimReference = (int)refSeq.size();

			double scoreDiag, scoreDelete, scoreInsert;
			MNMNAlignmentKernelAsymmetric *mnmn_kernel = dynamic_cast<MNMNAlignmentKernelAsymmetric *>(kernel.get());
			auto& cmpMatrix = isPure ? mnmn_kernel->getPureComparisonMatrix() : mnmn_kernel->getComparisonMatrix();

			double OpenInsertionScore = kernel->OpenInsertionScore(0, 0);
			double ContinueInsertionScore = kernel->ContinueInsertionScore(0, 0);
			double OpenDeletionScore = kernel->OpenDeletionScore(0, 0);
			double ContinueDeletionScore = kernel->ContinueDeletionScore(0, 0);

			//first element save max(insertScore + gap_extension, diagScore + gap_open)
			//second save max of 3.
			std::vector<double, AlignmentAllocator<double, 32>> rowVec(dimReference * 2);
			auto rowPtr = rowVec.data();

			double DOUBLE_MIN = std::numeric_limits<double>::lowest();

			//compute (0, 0)
			scoreDiag = kernel->ComparisonScore(0, 0); //tmp
			*rowPtr++ = scoreDiag + OpenInsertionScore;
			*rowPtr++ = scoreDiag;

			int max_row = 0; //-1
			int max_col = 0; //-1
			double maxRowScore = scoreDiag - cumulativeScoreQuery[0];

			// compute zeroth row
			for (int i1 = 1; i1 < dimReference; i1++)
			{
				scoreDiag = cmpMatrix[querySeq[0]][refSeq[i1]];
				if (maxRowScore < scoreDiag - cumulativeScoreQuery[0])
				{
					maxRowScore = scoreDiag - cumulativeScoreQuery[0];
					max_row = 0;
					max_col = i1;
				}
				*rowPtr++ = scoreDiag + OpenInsertionScore; //the first slot will contain the next Insertion score.
				*rowPtr++ = scoreDiag;
			}
			rowPtr = rowVec.data(); //reset to the begining

			//rowPtr1[0] = rowPtr1[1] = DOUBLE_MIN;

			auto col_cmpscores = mnmn_kernel->ComparisonScoreReference(0);
			// compute remaining rows
			for (int i0 = 1; i0 < dimQuery; i0++)
			{
				auto& cmpVec = cmpMatrix[querySeq[i0]];
				double cumScore = cumulativeScoreQuery[i0 - 1]; //wrong: should be i0-1

				rowPtr++; //skip the fisrt, then point to max of 3.
				double last_max = *rowPtr; //keep the max of 3 of upper row.
				//setting max of 3 of current row.
				double left_max = *rowPtr++ = cumScore + col_cmpscores[i0]; //now rowPtr point to the begining of next slot.

				double max_score = left_max;
				int max_j = 0;

				double nextScoreDelete = left_max + OpenDeletionScore;

				auto refseq_ptr = refSeq.data() + 1;
				for (int i1 = 1; i1 < dimReference; ++i1)
				{
					//layer 0
					scoreInsert = *rowPtr;

					scoreDelete = nextScoreDelete;

					//layer 1: matches: diag
					scoreDiag = FastMax(last_max, cumScore) + cmpVec[*refseq_ptr++];
					if (scoreDiag > max_score)
					{
						max_score = scoreDiag;
						max_j = i1;
					}

					// layer 2: deletions from sequence 0 : left
					nextScoreDelete = FastMax(scoreDiag + OpenDeletionScore, scoreDelete + ContinueDeletionScore);

					//resert the first element.
					*rowPtr++ = FastMax(scoreInsert + ContinueInsertionScore, scoreDiag + OpenInsertionScore);
					last_max = *rowPtr;

					//all maximum
					if (scoreDiag > scoreInsert)
					{
						*rowPtr++ = scoreDiag > scoreDelete ? scoreDiag : scoreDelete;
					}
					else
					{
						*rowPtr++ = scoreInsert > scoreDelete ? scoreInsert : scoreDelete;

					}
				}

				if (max_score - cumulativeScoreQuery[i0] > maxRowScore)
				{
					maxRowScore = max_score - cumulativeScoreQuery[i0];
					max_row = i0;
					max_col = max_j;
				}
				//reset array to begining.
				rowPtr = rowVec.data();
			}

			AlignmentScore = maxRowScore + cumulativeScoreQuery[dimQuery - 1];
			this->MaxRowIndex = max_row;
			this->MaxColIndex = max_col;
			return AlignmentScore;
		}

		double FillScoreMatrix_tom_ref()
		{
			dimQuery = (int)kernel->GetQuery()->Seq.size();
			dimReference = (int)kernel->GetReference()->Seq.size();

			if (dimQuery * dimReference * 4 > scoreMatrix.size())
			{
				scoreMatrix.resize(dimQuery * dimReference * 4 * 2);
			}

			double scoreDiag, scoreDelete, scoreInsert;
			double *matrixPtr0 = scoreMatrix.data();
			MNMNAlignmentKernelAsymmetric *mnmn_kernel = dynamic_cast<MNMNAlignmentKernelAsymmetric *>(kernel.get());
			mnmn_kernel->computeComparisonScoreOfReference();

			// compute element [0,0]
			double DOUBLE_MIN = std::numeric_limits<double>::lowest();
			scoreMatrix[0] = DOUBLE_MIN;
			scoreMatrix[1] = kernel->ComparisonScore(0, 0);
			scoreMatrix[2] = DOUBLE_MIN;
			scoreMatrix[3] = scoreMatrix[1];

			// compute zeroth row
			int index = 4;
			int row_length = dimReference * 4;
			auto row_cmpscores = mnmn_kernel->ComparisonScoreQuery(0);
			for (int i1 = 1; i1 < dimReference; i1++)
			{
				scoreMatrix[index++] = DOUBLE_MIN;
				scoreMatrix[index++] = row_cmpscores[i1]; //kernel->ComparisonScore(0, i1);
				scoreMatrix[index++] = DOUBLE_MIN;
				scoreMatrix[index++] = row_cmpscores[i1];
			}

			// compute zeroth column
			index = row_length;
			auto col_cmpscores = mnmn_kernel->ComparisonScoreReference(0);
			for (int i0 = 1; i0 < dimQuery; i0++, index += row_length)
			{
				scoreMatrix[index] = DOUBLE_MIN;
				scoreMatrix[index + 1] = cumulativeScoreQuery[i0 - 1] + col_cmpscores[i0]; //kernel->ComparisonScore(i0, 0);
				scoreMatrix[index + 2] = DOUBLE_MIN;
				scoreMatrix[index + 3] = scoreMatrix[index + 1];
			}

			double OpenInsertionScore = kernel->OpenInsertionScore(0, 0);
			double ContinueInsertionScore = kernel->ContinueInsertionScore(0, 0);
			double OpenDeletionScore = kernel->OpenDeletionScore(0, 0);
			double ContinueDeletionScore = kernel->ContinueDeletionScore(0, 0);

			//std::cerr << "consts: " << OpenInsertionScore << " " << ContinueInsertionScore << " " << OpenDeletionScore << " " << ContinueDeletionScore << std::endl;
			int max_row = -1;
			int max_col = -1;
			double maxRowScore = DOUBLE_MIN;
			double *lastRowPtr = matrixPtr0;
			double *rowPtr = lastRowPtr + row_length;

			// compute remaining rows
			for (int i0 = 1; i0 < dimQuery; i0++)
			{
				double cumScore = cumulativeScoreQuery[i0 - 1]; //wrong: should be i0-1
				auto& ComparisonScore = mnmn_kernel->ComparisonScoreQuery(i0); //wrong should be i0-1
				double max_score = DOUBLE_MIN;
				int max_j = -1;


				double nextScoreDelete = FastMax(rowPtr[1] + OpenDeletionScore, rowPtr[2] + ContinueDeletionScore);

				rowPtr += 4; //skep first slot.
				lastRowPtr += 4;
				double *cmpScorePtr = ComparisonScore.data() + 1;
				for (int i1 = 1; i1 < dimReference; ++i1, lastRowPtr += 4)
				{
					// layer 0: insertions into sequence 0: UP
					//scoreInsert = std::max(*lastRowPtr + ContinueInsertionScore, lastRowPtr[1] + OpenInsertionScore);
					scoreInsert = FastMax(*lastRowPtr + ContinueInsertionScore, lastRowPtr[1] + OpenInsertionScore);

					*rowPtr++ = scoreInsert;

					//layer 1: matches: diag
					//scoreDiag = std::max(lastRowPtr[-1], cumScore) + *cmpScorePtr++; //ComparisonScore[i1];
					scoreDiag = FastMax(lastRowPtr[-1], cumScore) + *cmpScorePtr++; //ComparisonScore[i1];;
					*rowPtr++ = scoreDiag;
					if (scoreDiag > max_score)
					{
						max_score = scoreDiag;
						max_j = i1;
					}

					// layer 2: deletions from sequence 0 : left
					//scoreDelete = std::max(rowPtr[-5] + OpenDeletionScore, rowPtr[-4] + ContinueDeletionScore);
					*rowPtr++ = scoreDelete = nextScoreDelete;
					nextScoreDelete = FastMax(scoreDiag + OpenDeletionScore, scoreDelete + ContinueDeletionScore);

					//scoreDelete = FastMax(rowPtr[-5] + OpenDeletionScore, rowPtr[-4] + ContinueDeletionScore);
					//*rowPtr++ = scoreDelete;


					//all maximum
					if (scoreDiag > scoreInsert)
					{
						*rowPtr++ = FastMax(scoreDiag, scoreDelete);
					}
					else
					{
						//*rowPtr++ = scoreInsert > scoreDelete ? scoreInsert : scoreDelete; //std::max(scoreInsert, scoreDelete);
						*rowPtr++ = FastMax(scoreInsert, scoreDelete);;

					}
				}

				if (max_score - cumulativeScoreQuery[i0] > maxRowScore)
				{
					maxRowScore = max_score - cumulativeScoreQuery[i0];
					max_row = i0;
					max_col = max_j;
				}


			}

			AlignmentScore = maxRowScore + cumulativeScoreQuery[dimQuery - 1];
			this->MaxRowIndex = max_row;
			this->MaxColIndex = max_col;

			return AlignmentScore;
		}

		///supposed to be same as Tom's currurnet version. optimized
		double FillScoreMatrix()
		{
			//auto tid = std::this_thread::get_id();
			//std::cerr << "thread" << tid << " is working\n";
			if (use_index_version)
			{
				return FillScoreMatrix_index();
			}

			//auto x1 = FillScoreMatrix_tom_ref();

			auto score = FillScoreMatrix_tom();
			//auto x2 = FillScoreMatrix_get_score();
			//if (score != x2)
			//{
			//	std::cerr << "score diff\n";
			//}
			return score;
		}
	
		//this supposedly matching Tom's result
		double FillScoreMatrix_tom(int mx = -1, int my=-1)
		{
			auto query = kernel->GetQuery();
			auto reference = kernel->GetReference();

			bool isPure = query->isPureNucleotide() && reference->isPureNucleotide();

			auto& querySeq = isPure ? query->PureSeq : query->Seq;
			auto& refSeq = isPure ? reference->PureSeq : reference->Seq;

			dimQuery = (int)querySeq.size();
			dimReference = (int)refSeq.size();

			if (dimQuery * dimReference * 3 > scoreMatrix.size())
			{
				scoreMatrix.resize(dimQuery * dimReference * 3 * 2);
			}

			double scoreDiag, scoreDelete, scoreInsert;
			double *matrixPtr0 = scoreMatrix.data();
			MNMNAlignmentKernelAsymmetric *mnmn_kernel = dynamic_cast<MNMNAlignmentKernelAsymmetric *>(kernel.get());
			auto& cmpMatrix = isPure ? mnmn_kernel->getPureComparisonMatrix() : mnmn_kernel->getComparisonMatrix();

			// compute element [0,0]
			double DOUBLE_MIN = std::numeric_limits<double>::lowest();
			scoreMatrix[0] = DOUBLE_MIN;
			scoreMatrix[1] = kernel->ComparisonScore(0, 0);
			scoreMatrix[2] = DOUBLE_MIN;

			int max_row = 0;
			int max_col = 0;
			double maxRowScore = scoreMatrix[1]- cumulativeScoreQuery[0];

			// compute zeroth row
			int index = 3;
			int row_length = dimReference * 3;
			auto& row_cmpscores = cmpMatrix[querySeq[0]];;
			for (int i1 = 1; i1 < dimReference; i1++)
			{
				scoreMatrix[index++] = DOUBLE_MIN;
				double tmp_score = row_cmpscores[refSeq[i1]];
				scoreMatrix[index++] = tmp_score;
				scoreMatrix[index++] = DOUBLE_MIN;
				if (maxRowScore < (tmp_score = tmp_score - cumulativeScoreQuery[0]))
				{
					maxRowScore = tmp_score;
					max_row = 0;
					max_col = i1;
				}
			}

			// compute zeroth column
			index = row_length;
			auto col_cmpscores = mnmn_kernel->ComparisonScoreReference(0);
			for (int i0 = 1; i0 < dimQuery; i0++, index += row_length)
			{
				scoreMatrix[index] = DOUBLE_MIN;
				scoreMatrix[index + 1] = cumulativeScoreQuery[i0 - 1] + col_cmpscores[i0];
				scoreMatrix[index + 2] = DOUBLE_MIN;
			}

			double OpenInsertionScore = kernel->OpenInsertionScore(0, 0);
			double ContinueInsertionScore = kernel->ContinueInsertionScore(0, 0);
			double OpenDeletionScore = kernel->OpenDeletionScore(0, 0);
			double ContinueDeletionScore = kernel->ContinueDeletionScore(0, 0);

			//std::cerr << "consts: " << OpenInsertionScore << " " << ContinueInsertionScore << " " << OpenDeletionScore << " " << ContinueDeletionScore << std::endl;

			double *lastRowPtr = matrixPtr0;
			double *rowPtr = lastRowPtr + row_length;

			auto querySeqPtr = querySeq.data() + 1;
			// compute remaining rows

			int i0max = mx == -1 ? dimQuery : mx + 1;
			//maxy cannot be used now, because pointers related to row length won't work in this version.
			//int i1max = my == -1 ? dimReference : my + 1;

			//for (int i0 = 1; i0 < dimQuery; i0++)
			for (int i0 = 1; i0 < i0max; i0++)
			{
				auto& cmpVec = cmpMatrix[*querySeqPtr++];
				double cumScore = cumulativeScoreQuery[i0 - 1]; //wrong: should be i0-1

				double left_max = rowPtr[1];
				scoreDelete = left_max + OpenDeletionScore;
				rowPtr += 3;

				lastRowPtr++;
				double last_max = *lastRowPtr++; 
				lastRowPtr++; //pointing to first element of next slot. 

				//to find maximum...
				double max_score = left_max; //last_max;
				int max_j = 0;

				auto refseq_ptr = refSeq.data() + 1;
				//for (int i1 = 1; i1 < dimReference; ++i1)
				for (int i1 = 1; i1 < dimReference; ++i1)
				{
					//layer 1 - diagnal
					scoreDiag = FastMax(last_max, cumScore) + cmpVec[*refseq_ptr++];
					if (scoreDiag > max_score)
					{
						max_score = scoreDiag;
						max_j = i1;
					}

					double a = *lastRowPtr++; //layer 0 score
					double b = *lastRowPtr++; //layer 1 score
					double c = *lastRowPtr++; //layer 2 score
					scoreInsert = FastMax(a + ContinueInsertionScore, b + OpenInsertionScore);
					last_max = FastMax(a, FastMax(b, c));

					*rowPtr++ = scoreInsert;
					*rowPtr++ = scoreDiag;
					*rowPtr++ = scoreDelete;
					//this compute "next" scoreDelete
					scoreDelete = FastMax(scoreDiag + OpenDeletionScore, scoreDelete + ContinueDeletionScore);
				}

				if (max_score - cumulativeScoreQuery[i0] > maxRowScore)
				{
					maxRowScore = max_score - cumulativeScoreQuery[i0];
					max_row = i0;
					max_col = max_j;
				}
			}

			AlignmentScore = maxRowScore + cumulativeScoreQuery[dimQuery - 1];
			this->MaxRowIndex = max_row;
			this->MaxColIndex = max_col;

			return AlignmentScore;
		}

		//while try to possibly going through the matrix only once.
		//then we needed to extra store operation of the two fields.
		//it actually slow thigns down, not helping.
		//so we might just check if the _pure version does better.
		double FillScoreMatrix_pure2_obsolate()
		{
			if (use_index_version)
			{
				return FillScoreMatrix_index();
			}

			auto& query = kernel->GetQuery();
			auto& reference = kernel->GetReference();

			bool isPure = query->isPureNucleotide() && reference->isPureNucleotide();

			auto& querySeq = isPure ? query->ToPureNucleotide() : query->Seq;
			auto& refSeq = isPure ? reference->ToPureNucleotide() : reference->Seq;

			dimQuery = (int)querySeq.size();
			dimReference = (int)refSeq.size();

			if (dimQuery * dimReference * 4 > scoreMatrix.size())
			{
				scoreMatrix.resize(dimQuery * dimReference * 4 * 2);
			}

			double scoreDiag, scoreDelete, scoreInsert;
			double *matrixPtr0 = scoreMatrix.data();
			MNMNAlignmentKernelAsymmetric *mnmn_kernel = dynamic_cast<MNMNAlignmentKernelAsymmetric *>(kernel.get());
			auto& cmpMatrix = isPure ? mnmn_kernel->getPureComparisonMatrix() : mnmn_kernel->getComparisonMatrix();

			double OpenInsertionScore = kernel->OpenInsertionScore(0, 0);
			double ContinueInsertionScore = kernel->ContinueInsertionScore(0, 0);
			double OpenDeletionScore = kernel->OpenDeletionScore(0, 0);
			double ContinueDeletionScore = kernel->ContinueDeletionScore(0, 0);

			std::vector<double, AlignmentAllocator<double, 64>> tmpvec(2 * dimReference);
			auto tmp_ptr = tmpvec.data();
			//1st element will store nextInsertion
			//2nd with store the maximum of all three.

			// compute element [0,0]
			double DOUBLE_MIN = std::numeric_limits<double>::lowest();
			scoreMatrix[0] = DOUBLE_MIN;
			scoreMatrix[1] = kernel->ComparisonScore(0, 0);
			scoreMatrix[2] = DOUBLE_MIN;
			scoreMatrix[3] = scoreMatrix[1];
			*tmp_ptr++ = 0;
			*tmp_ptr++ = scoreMatrix[1];

			// compute zeroth row
			int index = 4;
			int row_length = dimReference * 4;
			auto& row_cmpscores = cmpMatrix[querySeq[0]];;
			for (int i1 = 1; i1 < dimReference; i1++)
			{
				double tval = row_cmpscores[refSeq[i1]];
				scoreMatrix[index++] = DOUBLE_MIN;
				scoreMatrix[index++] = tval;
				scoreMatrix[index++] = DOUBLE_MIN;
				scoreMatrix[index++] = tval;

				*tmp_ptr++ = tval + OpenInsertionScore;
				*tmp_ptr++ = tval;
			}
			tmp_ptr = tmpvec.data(); //reset to begining...


			// compute zeroth column
			index = row_length;
			auto col_cmpscores = mnmn_kernel->ComparisonScoreReference(0);
			for (int i0 = 1; i0 < dimQuery; i0++, index += row_length)
			{
				scoreMatrix[index] = DOUBLE_MIN;
				scoreMatrix[index + 1] = cumulativeScoreQuery[i0 - 1] + col_cmpscores[i0];
				scoreMatrix[index + 2] = DOUBLE_MIN;
				scoreMatrix[index + 3] = scoreMatrix[index + 1];
			}



			//std::cerr << "consts: " << OpenInsertionScore << " " << ContinueInsertionScore << " " << OpenDeletionScore << " " << ContinueDeletionScore << std::endl;
			int max_row = -1;
			int max_col = -1;
			double maxRowScore = DOUBLE_MIN;
			double *rowPtr = matrixPtr0 + row_length;

			auto querySeqPtr = querySeq.data() + 1;
			// compute remaining rows
			for (int i0 = 1; i0 < dimQuery; i0++)
			{
				auto& cmpVec = cmpMatrix[*querySeqPtr++];

				double cumScore = cumulativeScoreQuery[i0 - 1]; //wrong: should be i0-1
				double max_score = DOUBLE_MIN;
				int max_j = -1;

				rowPtr += 4; //skep first slot.
				
				//there seems some issue here...
				tmp_ptr++; //skip the first, which is not used anywhere
				double last_max = *tmp_ptr;	//maximum value of 3 - computed in last row.
				double left_max = *tmp_ptr++ = cumScore + col_cmpscores[i0]; //max of 3 set for this row.
				//now the tmp_ptr pointer point to the slot 1, which is the scorInsert.
				//*tmp_ptr = cumScore + col_cmpscores[i0]; //max of 3

				double nextScoreDelete = left_max + OpenDeletionScore;

				auto refseq_ptr = refSeq.data() + 1;
				for (int i1 = 1; i1 < dimReference; ++i1)
				{
					// layer 0: insertions into sequence 0: UP
					*rowPtr++ = scoreInsert = *tmp_ptr;

					//layer 1: matches: diag
					//scoreDiag = FastMax(lastRowPtr[-1], cumScore) + *cmpScorePtr++;
					scoreDiag = FastMax(last_max, cumScore) + cmpVec[*refseq_ptr++];
					*rowPtr++ = scoreDiag;

					if (scoreDiag > max_score)
					{
						max_score = scoreDiag;
						max_j = i1;
					}

					// layer 2: deletions from sequence 0 : left
					//scoreDelete = std::max(rowPtr[-5] + OpenDeletionScore, rowPtr[-4] + ContinueDeletionScore);
					*rowPtr++ = scoreDelete = nextScoreDelete;

					nextScoreDelete = FastMax(scoreDiag + OpenDeletionScore, scoreDelete + ContinueDeletionScore);

					//nextScoreInsert
					*tmp_ptr++ = FastMax(scoreInsert + ContinueInsertionScore, scoreDiag + OpenInsertionScore);

					last_max = *tmp_ptr;
					//all maximum
					if (scoreDiag > scoreInsert)
					{
						*tmp_ptr++ = FastMax(scoreDiag, scoreDelete);
						*rowPtr++ = FastMax(scoreDiag, scoreDelete);
					}
					else
					{
						*tmp_ptr++ = FastMax(scoreInsert, scoreDelete);
						*rowPtr ++ = FastMax(scoreDiag, scoreDelete);
					}
				}

				if (max_score - cumulativeScoreQuery[i0] > maxRowScore)
				{
					maxRowScore = max_score - cumulativeScoreQuery[i0];
					max_row = i0;
					max_col = max_j;
				}
				tmp_ptr = tmpvec.data();
			}

			AlignmentScore = maxRowScore + cumulativeScoreQuery[dimQuery - 1];
			this->MaxRowIndex = max_row;
			this->MaxColIndex = max_col;

			return AlignmentScore;
		}

		///supposed to be same as Tom's currurnet version.
		double FillScoreMatrix_op0(bool qFivePrimeEndFree = true, bool rFivePrimeEndFree = true, bool qThreePrimeEndFree = true, bool rThreePrimeEndFree = true)
		{

			dimQuery = (int)kernel->GetQuery()->Seq.size();
			dimReference = (int)kernel->GetReference()->Seq.size();

			//std::vector<byte>& s = kernel->GetReference()->Seq;
			//auto cv = TypeConverter::ToChars(s);
			//std::string seq2(cv.begin(), cv.end());
			//std::cerr << seq2 << std::endl;

			//std::cerr << "vector size: " << scoreMatrix.size() << std::endl;
			//std::cerr << "need size: " << (dimQuery * dimReference * 3) << std::endl;

			if (dimQuery * dimReference * 3 > scoreMatrix.size())
			{
				scoreMatrix.resize(dimQuery * dimReference * 3);
				traceBackMatrix.resize(dimQuery * dimReference * 3);
			}

			int dim0Minus1 = dimQuery - 1;
			int dim1Minus1 = dimReference - 1;
			int i0 = 0;
			int i1 = 0;

			double scoreDiag, scoreDiagMax;
			double scoreDelete, scoreDeleteMax;
			double scoreInsert, scoreInsertMax;
			double jumpScore;
			bool bothFivePrimeEndsFree = qFivePrimeEndFree && rFivePrimeEndFree;

			MNMNAlignmentKernelAsymmetric *mnmn_kernel = dynamic_cast<MNMNAlignmentKernelAsymmetric *>(kernel.get());
			mnmn_kernel->computeComparisonScoreOfReference();


			// compute element [0,0]
			double DOUBLE_MIN = std::numeric_limits<double>::lowest();
			//std::cerr << "doublemin=" << DOUBLE_MIN << std::endl;

			scoreMatrix[0] = DOUBLE_MIN;
			scoreMatrix[1] = kernel->ComparisonScore(0, 0);
			traceBackMatrix[1] = -1;
			scoreMatrix[2] = DOUBLE_MIN;

			// compute zeroth row
			int index = 3;
			auto row_cmpscores = mnmn_kernel->ComparisonScoreQuery(0);
			for (i1 = 1; i1 < dimReference; i1++)
			{
				scoreMatrix[index++] = DOUBLE_MIN;
				scoreMatrix[index] = row_cmpscores[i1]; //kernel->ComparisonScore(0, i1);
				traceBackMatrix[index++] = -1;
				scoreMatrix[index++] = DOUBLE_MIN;
			}

			// compute zeroth column
			index = dimReference * 3;
			auto col_cmpscores = mnmn_kernel->ComparisonScoreReference(0);
			for (i0 = 1; i0 < dimQuery; i0++)
			{
				scoreMatrix[index] = DOUBLE_MIN;
				scoreMatrix[index + 1] = cumulativeScoreQuery[i0 - 1] + col_cmpscores[i0]; //kernel->ComparisonScore(i0, 0);
				traceBackMatrix[index + 1] = -1;
				scoreMatrix[index + 2] = DOUBLE_MIN;
				index += dimReference * 3;
			}

			double OpenInsertionScore = kernel->OpenInsertionScore(0, 0);
			double ContinueInsertionScore = kernel->ContinueInsertionScore(0, 0);
			double OpenDeletionScore = kernel->OpenDeletionScore(0, 0);
			double ContinueDeletionScore = kernel->ContinueDeletionScore(0, 0);

			//std::cerr << "consts: " << OpenInsertionScore << " " << ContinueInsertionScore << " " << OpenDeletionScore << " " << ContinueDeletionScore << std::endl;

			int max_row = -1;
			int max_col = -1;
			double maxRowScore = DOUBLE_MIN;
			// compute remaining rows
			for (i0 = 1; i0 < dimQuery; i0++)
			{
				int i0Index = i0 * dimReference * 3 + 3;
				int i0IndexLast = i0Index - dimReference * 3;
				double cumScore = cumulativeScoreQuery[i0 - 1]; //wrong: should be i0-1

				auto& ComparisonScore = mnmn_kernel->ComparisonScoreQuery(i0); //wrong should be i0-1
				double max_score = DOUBLE_MIN;
				int max_j = -1;
				for (i1 = 1; i1 < dimReference; i1++, i0Index += 3, i0IndexLast += 3)
				{
					scoreInsertMax = scoreMatrix[i0IndexLast] + ContinueInsertionScore;
					scoreInsert = scoreMatrix[i0IndexLast + 1] + OpenInsertionScore;
					if (scoreInsertMax >= scoreInsert)
					{
						scoreMatrix[i0Index] = scoreInsertMax;
						traceBackMatrix[i0Index] = 0;
					}
					else
					{
						scoreMatrix[i0Index] = scoreInsert;
						traceBackMatrix[i0Index] = 1;
					}


					// layer 1: matches
					int tmpIndex = i0IndexLast - 3;
					scoreDiagMax = std::max(scoreMatrix[tmpIndex], scoreMatrix[tmpIndex + 1]);
					scoreDiag = std::max(scoreMatrix[tmpIndex + 2], cumScore);

					//if (i0 == 1 && i1 == 1)
					//{
					//	std::cerr << "==============================================\n";
					//	std::cerr << "matchcost= " << ComparisonScore[i1] << "\n";
					//	std::cerr << "score0= " << scoreMatrix[tmpIndex] << "\n";
					//	std::cerr << "score1= " << scoreMatrix[tmpIndex+1] << "\n";
					//	std::cerr << "score2= " << scoreMatrix[tmpIndex+2] << "\n";
					//	std::cerr << "jumpscore (before cost)= " << cumScore << "\n";
					//	std::cerr << "maxscore= " << std::max(scoreDiag, scoreDiagMax) + ComparisonScore[i1] << "\n";
					//	std::cerr << "==============================================\n";
					//}


					if (scoreDiagMax >= scoreDiag)
					{
						scoreMatrix[i0Index + 1] = scoreDiagMax + ComparisonScore[i1];
						traceBackMatrix[i0Index + 1] = scoreMatrix[tmpIndex] > scoreMatrix[tmpIndex + 1] ? 0 : 1;
					}
					else
					{
						scoreMatrix[i0Index + 1] = scoreDiag + ComparisonScore[i1]; //wrong should be i0-1 : we should compute it so that htere i sno need to-1
						traceBackMatrix[i0Index + 1] = scoreMatrix[tmpIndex + 2] > cumScore ? 2 : -1;
					}

					if (scoreMatrix[i0Index + 1] > max_score)
					{
						max_score = scoreMatrix[i0Index + 1];
						max_j = i1;
					}

					// layer 2: deletions from sequence 0
					scoreDeleteMax = scoreMatrix[i0Index - 3 + 1] + OpenDeletionScore;
					scoreDelete = scoreMatrix[i0Index - 3 + 2] + ContinueDeletionScore;
					if (scoreDeleteMax > scoreDelete)
					{
						scoreMatrix[i0Index + 2] = scoreDeleteMax;
						traceBackMatrix[i0Index + 2] = 1;
					}
					else
					{
						scoreMatrix[i0Index + 2] = scoreDelete;
						traceBackMatrix[i0Index + 2] = 2;
					}

				}


				if (max_score - cumulativeScoreQuery[i0] > maxRowScore)
					//if (max_score - cumScore > maxRowScore)
				{
					maxRowScore = max_score - cumulativeScoreQuery[i0];
					max_row = i0;
					max_col = max_j;
				}

			}


			AlignmentScore = maxRowScore + cumulativeScoreQuery[dimQuery - 1];
			this->MaxRowIndex = max_row;
			this->MaxColIndex = max_col;

			//debug
			//if (DEBUG_CODE != 0)
			{
				std::cerr << "=============native output===================\n";
				//write_matrix_to_file(dimQuery, dimReference);
				fprintf(stderr, "maxi=%d maxj=%d score=%f\n", max_row, max_col, AlignmentScore);
				std::cerr << "base score= " << maxRowScore << " addscore= " << cumulativeScoreQuery[dimQuery - 1] << "\n";
				std::cerr << "==============================================\n";
			}

			//traceBackStart = std::vector<int>{ max_row, max_col };
			return AlignmentScore;
		}


		//for debug only
		void write_matrix_to_file(int dim1, int dim2)
		{
			std::cerr << "writing matrix...\n";

			ofstream debug_stream;
			debug_stream.open("debug-2.txt");

			for (int i = 0; i < dim1; i++)
			{
				for (int j = 0; j < dim2; j++)
				{

					debug_stream << i << " " << j << ": " << scoreMatrix[MATRIX_INDEX_4(i, j, 0)] << "\t" << scoreMatrix[MATRIX_INDEX_4(i, j, 1)] << "\t" << scoreMatrix[MATRIX_INDEX_4(i, j, 2)] << "\n";
					//debug_stream << i0 << " " << i1 << ": " << traceBackMatrix[i0Index - 2] << "\t" << traceBackMatrix[i0Index - 1] << "\t" << traceBackMatrix[i0Index] << "\n";
				}
			}

			debug_stream.close();
		}

		//version "index modifed".
		//this is the version with indexes modified. notice- it has to match the traceback version.
		double FillScoreMatrix_index(bool qFivePrimeEndFree = true, bool rFivePrimeEndFree = true, bool qThreePrimeEndFree = true, bool rThreePrimeEndFree = true)
		{
			//std::cerr << "in FillScoreMatrix.\n";
			dimQuery = (int)kernel->GetQuery()->Seq.size() + 1;
			dimReference = (int)kernel->GetReference()->Seq.size() + 1;
			int dim2 = dimReference;

			//while (DEBUG_CODE == 1)
			//{
			//}

			if (dimQuery * dimReference * 3 > scoreMatrix.size())
			{
				scoreMatrix.resize(dimQuery * dimReference * 3);
				traceBackMatrix.resize(dimQuery * dimReference * 3);
			}

			//std::vector<byte>& s = kernel->GetReference()->Seq;
			//auto cv = TypeConverter::ToChars(s);
			//std::string seq2(cv.begin(), cv.end());
			//std::cerr << seq2 << std::endl;

			//std::cerr << "vector size: " << scoreMatrix.size() << std::endl;
			//std::cerr << "need size: " << (dimQuery * dimReference * 3) << std::endl;

			int dim0Minus1 = dimQuery - 1;
			int dim1Minus1 = dimReference - 1;
			int i0 = 0;
			int i1 = 0;

			double scoreDiag, scoreDiagMax;
			double scoreDelete, scoreDeleteMax;
			double scoreInsert, scoreInsertMax;
			double jumpScore;
			bool bothFivePrimeEndsFree = qFivePrimeEndFree && rFivePrimeEndFree;

			MNMNAlignmentKernelAsymmetric *mnmn_kernel = dynamic_cast<MNMNAlignmentKernelAsymmetric *>(kernel.get());
			mnmn_kernel->computeComparisonScoreOfReference();


			// compute element [0,0]
			double DOUBLE_MIN = std::numeric_limits<double>::lowest();
			//std::cerr << "doublemin=" << DOUBLE_MIN << std::endl;

			scoreMatrix[0] = DOUBLE_MIN;
			scoreMatrix[1] = 0; //kernel->ComparisonScore(0, 0);
			traceBackMatrix[1] = -1;
			scoreMatrix[2] = DOUBLE_MIN;

			// compute zeroth row
			int index = 3;
			auto& cmpscore = mnmn_kernel->ComparisonScoreQuery(0);
			for (i1 = 1; i1 < dimReference; ++i1)
			{
				scoreMatrix[index] = DOUBLE_MIN;
				traceBackMatrix[index++] = -1;
				scoreMatrix[index] = DOUBLE_MIN; //cmpscore[i1 - 1]; //kernel->ComparisonScore(0, i1);
				traceBackMatrix[index++] = -1;
				scoreMatrix[index] = 0;
				traceBackMatrix[index++] = rFivePrimeEndFree ? -1 : 2;
			}

			// compute zeroth column
			index = dimReference * 3;
			for (i0 = 1; i0 < dimQuery; ++i0)
			{
				scoreMatrix[index] = cumulativeScoreQuery[i0 - 1];
				traceBackMatrix[index] = qFivePrimeEndFree ? -1 : 0;
				scoreMatrix[index + 1] = DOUBLE_MIN;
				traceBackMatrix[index + 1] = -1;
				scoreMatrix[index + 2] = DOUBLE_MIN;
				traceBackMatrix[index + 1] = -1;
				index += dimReference * 3;
			}

			double OpenInsertionScore = kernel->OpenInsertionScore(0, 0);
			double ContinueInsertionScore = kernel->ContinueInsertionScore(0, 0);
			double OpenDeletionScore = kernel->OpenDeletionScore(0, 0);
			double ContinueDeletionScore = kernel->ContinueDeletionScore(0, 0);

			double insertionScoreDiff = ContinueInsertionScore - OpenInsertionScore;
			double deletionScoreDiff = ContinueDeletionScore - OpenDeletionScore;

			int max_row = -1;
			int max_col = -1;
			double maxRowScore = DOUBLE_MIN;

			// compute remaining rows
			for (i0 = 1; i0 < dimQuery; ++i0)
			{
				int i0Index = i0 * dimReference * 3 + 3;
				int i0IndexLast = i0Index - dimReference * 3;
				double cumScore = cumulativeScoreQuery[i0 - 1]; //wrong: should be i0-1

				auto& ComparisonScore = mnmn_kernel->ComparisonScoreQuery(i0 - 1); //wrong should be i0-1
				double max_score = DOUBLE_MIN;
				int max_j = -1;

				//const double *dataRowPtr = scoreMatrix.data();
				for (i1 = 1; i1 < dimReference; ++i1, ++i0Index, i0IndexLast += 3)
				{
					scoreInsertMax = scoreMatrix[i0IndexLast] + ContinueInsertionScore;
					scoreInsert = scoreMatrix[i0IndexLast + 1] + OpenInsertionScore;
					if (scoreInsertMax >= scoreInsert)
					{
						scoreMatrix[i0Index] = scoreInsertMax;
						traceBackMatrix[i0Index] = 0;
					}
					else
					{
						scoreMatrix[i0Index] = scoreInsert;
						traceBackMatrix[i0Index] = 1;
					}

					// layer 1: matches
					int tmpIndex = i0IndexLast - 3;
					i0Index++;
					scoreDiagMax = std::max(scoreMatrix[tmpIndex], scoreMatrix[tmpIndex + 1]);
					scoreDiag = std::max(scoreMatrix[tmpIndex + 2], cumScore- ComparisonScore[i1-1]); //<=== note the minus of cumScore and blow.

					if (scoreDiagMax >= scoreDiag)
					{
						scoreMatrix[i0Index] = scoreDiagMax + ComparisonScore[i1 - 1];;
						traceBackMatrix[i0Index] = scoreMatrix[tmpIndex] >= scoreMatrix[tmpIndex + 1] ? 0 : 1;
					}
					else
					{
						scoreMatrix[i0Index] = scoreDiag + ComparisonScore[i1 - 1]; //wrong should be i1-1 : we should compute it so that htere i sno need to-1
						traceBackMatrix[i0Index] = scoreMatrix[tmpIndex + 2] >= cumScore - ComparisonScore[i1 - 1] ? 2 : -1;
					}
					if (scoreMatrix[i0Index] > max_score)
					{
						max_score = scoreMatrix[i0Index];
						max_j = i1;
					}

					// layer 2: deletions from sequence 0
					i0Index++; //i0Index = I0Index + 2;
					scoreDeleteMax = scoreMatrix[i0Index - 5 + 1] + OpenDeletionScore;
					scoreDelete = scoreMatrix[i0Index - 3] + ContinueDeletionScore;
					if (scoreDeleteMax >= scoreDelete)
					{
						scoreMatrix[i0Index] = scoreDeleteMax;
						traceBackMatrix[i0Index] = 1;
					}
					else
					{
						scoreMatrix[i0Index] = scoreDelete;
						traceBackMatrix[i0Index] = 2;
					}
				}

				//when trace back, i0 th score is used, not i0-1????????
				if (max_score - cumScore > maxRowScore)
					//if (max_score - cumulativeScoreQuery[i0] > maxRowScore)
				{
					maxRowScore = max_score - cumScore; //cumulativeScoreQuery[i0]????;
					max_row = i0;
					max_col = max_j;
				}

			}

			//debug
			if (DEBUG_CODE != 0)
			{
				//write_matrix_to_file(dimQuery, dimReference);
			}

			AlignmentScore = maxRowScore + cumulativeScoreQuery[dimQuery - 2];
			MaxRowIndex = max_row;
			MaxColIndex = max_col;

			//this->kernel->print_debug_info();
			//fprintf(stderr, "dimQ=%d dimR=%d maxi=%d maxj=%d score=%f\n", dimQuery, dimReference, max_row, max_col, AlignmentScore);
			return AlignmentScore;
		}


		double FillScoreMatrix_original(bool qFivePrimeEndFree = true, bool rFivePrimeEndFree = true, bool qThreePrimeEndFree = true, bool rThreePrimeEndFree = true)
		{

			dimQuery = (int)kernel->GetQuery()->Seq.size();
			dimReference = (int)kernel->GetReference()->Seq.size();


			int dim0Minus1 = dimQuery - 1;
			int dim1Minus1 = dimReference - 1;
			int i0 = 0;
			int i1 = 0;

			double scoreDiag, scoreDiagMax;
			double scoreDelete, scoreDeleteMax;
			double scoreInsert, scoreInsertMax;
			double jumpScore;
			bool bothFivePrimeEndsFree = qFivePrimeEndFree && rFivePrimeEndFree;


			// compute element [0,0]
			scoreMatrix[0] = std::numeric_limits<double>::lowest();
			scoreMatrix[1] = kernel->ComparisonScore(0, 0);
			traceBackMatrix[1] = -1;
			scoreMatrix[2] = std::numeric_limits<double>::lowest();

			// compute zeroth row
			for (i1 = 1; i1 < dimReference; i1++)
			{
				scoreMatrix[i1 * 3] = std::numeric_limits<double>::lowest();
				scoreMatrix[i1 * 3 + 1] = kernel->ComparisonScore(0, i1);
				traceBackMatrix[i1 * 3 + 1] = rFivePrimeEndFree ? -1 : 1;;
				scoreMatrix[i1 * 3 + 2] = std::numeric_limits<double>::lowest();
			}

			// compute zeroth column
			for (i0 = 1; i0 < dimQuery; i0++)
			{
				scoreMatrix[i0 * dimReference * 3] = std::numeric_limits<double>::lowest();
				scoreMatrix[i0 * dimReference * 3 + 1] = cumulativeScoreQuery[i0 - 1] + kernel->ComparisonScore(i0, 0);
				traceBackMatrix[i0 * dimReference * 3 + 1] = qFivePrimeEndFree ? -1 : 1;;
				scoreMatrix[i0 * dimReference * 3 + 2] = std::numeric_limits<double>::lowest();
			}

			// compute remaining rows
			for (i0 = 1; i0 < dimQuery; i0++)
			{
				int i0Index = i0 * dimReference * 3;
				for (i1 = 1; i1 < dimReference; i1++)
				{

					//double tt1 = kernel->ContinueInsertionScore(i0, i1);
					//double tt2 = kernel->OpenInsertionScore(i0, i1);
					//double tt3 = kernel->ComparisonScore(i0, i1);
					//double tt4 = kernel->ContinueDeletionScore(i0, i1);
					//double tt5 = kernel->OpenDeletionScore(i0, i1);

					//if (i0 == 1 && i1 == 1)
					//{
					//	fprintf(stderr, "%d %d : %f %f %f %f %f\n", i0, i1, tt1, tt2, tt3, tt4, tt5);
					//}

					// layer 0: insertions into sequence 0
					scoreInsertMax = std::numeric_limits<double>::lowest();
					for (int iLayer = 0; iLayer < 2; iLayer++)
					{
						scoreInsert = scoreMatrix[(i0 - 1) * dimReference * 3 + i1 * 3 + iLayer] + (iLayer == 0 ? kernel->ContinueInsertionScore(i0, i1) : kernel->OpenInsertionScore(i0, i1));
						if (scoreInsert > scoreInsertMax)
						{
							scoreInsertMax = scoreInsert;
							traceBackMatrix[i0Index + i1 * 3] = iLayer;
						}
						scoreMatrix[i0Index + i1 * 3] = scoreInsertMax;
					}

					// layer 1: matches
					scoreDiagMax = std::numeric_limits<double>::lowest();
					double matchCost = 0;
					matchCost = kernel->ComparisonScore(i0, i1);
					for (int iLayer = 0; iLayer < 3; iLayer++)
					{
						scoreDiag = scoreMatrix[(i0 - 1) * dimReference * 3 + (i1 - 1) * 3 + iLayer] + matchCost;
						if (scoreDiag > scoreDiagMax)
						{
							scoreDiagMax = scoreDiag;
							traceBackMatrix[i0Index + i1 * 3 + 1] = iLayer;
						}
					}

					if (bothFivePrimeEndsFree)
					{
						jumpScore = cumulativeScoreQuery[i0 - 1] + matchCost;
						if (jumpScore > scoreDiagMax)
						{
							scoreDiagMax = jumpScore;
							traceBackMatrix[i0Index + i1 * 3 + 1] = -1;
						}
					}

					scoreMatrix[i0Index + i1 * 3 + 1] = scoreDiagMax;

					// layer 2: deletions from sequence 0
					scoreDeleteMax = std::numeric_limits<double>::lowest();
					for (int iLayer = 1; iLayer < 3; iLayer++)
					{
						scoreDelete = scoreMatrix[i0Index + (i1 - 1) * 3 + iLayer] + (iLayer == 2 ? kernel->ContinueDeletionScore(i0, i1) : kernel->OpenDeletionScore(i0, i1));
						if (scoreDelete > scoreDeleteMax)
						{
							scoreDeleteMax = scoreDelete;
							traceBackMatrix[i0Index + i1 * 3 + 2] = iLayer;
						}
						scoreMatrix[i0Index + i1 * 3 + 2] = scoreDeleteMax;
					}
				}
			}

			// find starting point and alignment score for traceback

			int i1Max = -1;
			int i0Max = -1;
			double currentScore = 0;
			double maxScore = 0;
			for (i0 = 0; i0 < dimQuery; i0++)
			{
				int i0Index = i0 * dimReference * 3;
				for (i1 = 0; i1 < dimReference; i1++)
				{
					currentScore = scoreMatrix[i0Index + i1 * 3 + 1] - cumulativeScoreQuery[i0];
					if (currentScore > maxScore)
					{
						maxScore = currentScore;
						i1Max = i1;
						i0Max = i0;
					}
				}
			}

			AlignmentScore = maxScore + cumulativeScoreQuery[dimQuery - 1];
			//traceBackStart = std::vector<int>{ i0Max, i1Max };
			this->MaxRowIndex = i0Max;
			this->MaxColIndex = i1Max;
			return AlignmentScore;
		}


		//save as oringail, not used in indexed version.
		std::vector<int> GetTracebackStart(bool qThreePrimeEndFree, bool rThreePrimeEndFree)
		{
			// find starting point and alignment score for traceback
			double maxScore = std::numeric_limits<double>::lowest();;
			int i0AtEnd = -1;
			int i1AtEnd = -1;
			int lastQIndex = dimQuery - 1;
			int lastRIndex = dimReference - 1;

			dimReference = (int)kernel->GetReference()->Seq.size();
			int dim2 = dimReference;

			if (qThreePrimeEndFree && rThreePrimeEndFree)
			{
				double currentScore = 0;
				maxScore = currentScore;
				for (int i0 = 0; i0 < dimQuery; i0++)
				{
					for (int i1 = 0; i1 < dimReference; i1++)
					{
						currentScore = scoreMatrix[MATRIX_INDEX(i0, i1, 1)] - cumulativeScoreQuery[i0];
						if (currentScore > maxScore)
						{
							maxScore = currentScore;
							i1AtEnd = i1;
							i0AtEnd = i0;
						}
					}
				}
				AlignmentScore = maxScore + cumulativeScoreQuery[lastQIndex];
			}
			else if (qThreePrimeEndFree) // the traceback must start at the end of the reference sequence
			{
				double currentScore = 0;
				i1AtEnd = lastRIndex;
				maxScore = currentScore;
				for (int i0 = 0; i0 < dimQuery; i0++)
				{
					currentScore = scoreMatrix[MATRIX_INDEX(i0, lastRIndex, 1)] - cumulativeScoreQuery[i0];
					if (currentScore > maxScore)
					{
						maxScore = currentScore;
						i0AtEnd = i0;
					}
				}
				AlignmentScore = maxScore + cumulativeScoreQuery[lastQIndex];
			}
			else if (rThreePrimeEndFree)// traceback starts from end of query
			{
				double currentScore = 0;
				maxScore = currentScore;
				i0AtEnd = lastQIndex;
				for (int i1 = 0; i1 < dimReference; i1++)
				{
					currentScore = scoreMatrix[MATRIX_INDEX(lastQIndex, i1, 1)];
					if (currentScore > maxScore)
					{
						maxScore = currentScore;
						i1AtEnd = i1;
					}
				}
				AlignmentScore = maxScore;
			}
			else
			{
				i0AtEnd = lastQIndex;
				i1AtEnd = lastRIndex;
				maxScore = scoreMatrix[MATRIX_INDEX(lastQIndex, lastRIndex, 1)];
			}
			return std::vector<int>{ i0AtEnd, i1AtEnd };
		}

		/// <summary>
		/// Performs traceback for a dynamic programming-based pairwise alignment.
		/// not used, 
		/// </summary>
		PolymerPair<T, U> TraceBack(bool extend = false)
		{

			auto start = std::vector<int>{ this->MaxRowIndex, this->MaxColIndex };
			// this method assumes that the Score matrix has been appropriately filled
			// in the initial stage of the process, the Index function returns a 
			// value that is displaced by one from the necessary value.
			// This condition is corrected in the last step

			PolymerPair<T, U> pair(kernel->GetQuery(), kernel->GetReference());
			pair.AlignmentScore = this->AlignmentScore;

			if (start[0] < 0 || start[1] < 0)
			{
				//vector is initialized with size of 0
				//pair.Index = new int[0][];
				//pair.Score = new double[0];
				return pair;
			}

			// start traceback from far corner
			int x = start[0];
			int y = start[1];
			int currentLayer = 1;
			if (start.size() > 2)
			{
				currentLayer = start[2];
			}

			// stop when the top row is passed
			stack<int> indexXStack;		//= new Stack<int>();
			stack<int> indexYStack;		// = new Stack<int>();
			stack<double> scoreStack; // = new Stack<double>();

			int dimReference = (int)kernel->GetReference()->Seq.size();
			int indexLength = 0;
			bool foundMatch = false;
			int nextLayer;
			while (currentLayer > -1)
			{
				int matrixIndex = x * dimReference * 3 + y * 3 + currentLayer;
				nextLayer = traceBackMatrix[matrixIndex];
				scoreStack.push(scoreMatrix[matrixIndex]);
				indexLength++;
				switch (currentLayer)
				{
				case -1: // termination
					indexXStack.push(x);
					indexYStack.push(y);
					x = 0;
					y = 0;
					break;
				case 0: // insertion
					indexXStack.push(x);
					indexYStack.push(-1);
					x -= 1;
					break;
				case 1: // diagonal move
					if (!foundMatch)
					{
						foundMatch = true;
					}
					indexXStack.push(x);
					indexYStack.push(y);
					x -= 1;

					// when the first column is reached, we must prevent y from becoming negative:
					y = max(0, y - 1);

					break;
				case 2: // deletion
					indexXStack.push(-1);
					indexYStack.push(y);
					y -= 1;
					break;
				}
				currentLayer = nextLayer;
			}

			if (extend)
			{
				// By default, the alignment ends when the beginning of either sequence is reached.
				// This block is entered if the user wants both sequences to be completed.
				// The score is not available for this portion of the alignment.
				int nextX = indexXStack.top() - 1;
				for (; nextX > -1; nextX--)
				{
					indexXStack.push(nextX);
					indexYStack.push(-1);
					scoreStack.push(std::numeric_limits<double>::quiet_NaN());
					indexLength++;
				}
				int nextY = indexYStack.top() - 1;
				for (; nextY > -1; nextY--)
				{
					indexXStack.push(-1);
					indexYStack.push(nextY);
					scoreStack.push(std::numeric_limits<double>::quiet_NaN());
					indexLength++;
				}
			}

			// declare the index array and fill it 
			pair.Index.resize(indexLength);
			pair.Score.resize(indexLength);

			for (int i = 0; i < indexLength; i++)
			{
				pair.Index[i].push_back(indexXStack.top());
				indexXStack.pop();
				pair.Index[i].push_back(indexYStack.top());
				indexYStack.pop();
				pair.Score.push_back(scoreStack.top());
				scoreStack.pop();
			}

			return pair;
		}

		///result matching Tom's vesion - this assumes
		//the scores are stored in the order of layers[0, 1, 2] 
		int TraceBack_tom(int *index1, int *index2, double *scores, int arrayLength, bool complete0 = false, bool complete1 = true)
		{

			if (MaxRowIndex < 0 || MaxColIndex < 0)
			{
				return 0;
			}

			// start traceback from far corner
			int x = MaxRowIndex;
			int y = MaxColIndex;
			int currentLayer = 1;

			// stop when the top row is passed
			stack<int> indexXStack;		//= new Stack<int>();
			stack<int> indexYStack;		// = new Stack<int>();
			stack<double> scoreStack; // = new Stack<double>();

			double OpenInsertionScore = kernel->OpenInsertionScore(0, 0);
			double ContinueInsertionScore = kernel->ContinueInsertionScore(0, 0);
			double OpenDeletionScore = kernel->OpenDeletionScore(0, 0);
			double ContinueDeletionScore = kernel->ContinueDeletionScore(0, 0);

			MNMNAlignmentKernelAsymmetric *mnmn_kernel = dynamic_cast<MNMNAlignmentKernelAsymmetric *>(kernel.get());
			if (mnmn_kernel == nullptr)
			{
				fprintf(stderr, "mnmn_kernel null!");
				exit(1);
			}
			int dimReference = (int)kernel->GetReference()->Seq.size();
			int indexLength = 0;
			bool foundMatch = false;
			int nextLayer;
			bool finishing0 = false;
			bool finishing1 = false;
			int dimRefer3 = dimReference * 3;
			double *matrixPtr = scoreMatrix.data();
			double *dptr;
			double cumu_score;
			int maxval_index;
			while (currentLayer > -1)
			{
				int matrixIndex = x * dimRefer3 + y * 3 + currentLayer;
				switch (currentLayer)
				{
				case 0:
					nextLayer = scoreMatrix[matrixIndex - dimRefer3] + ContinueInsertionScore > scoreMatrix[matrixIndex - dimRefer3 + 1] + OpenInsertionScore ? 0 : 1;
					break;
				case 1:
					//find maximum value of slot[i-1][j-1]
					dptr = matrixPtr + matrixIndex - dimRefer3 - 4;
					if (dptr >= matrixPtr)
					{
						maxval_index = 0;
						if (dptr[1] > dptr[maxval_index])
						{
							maxval_index = 1;
						}
						if (dptr[2] > dptr[maxval_index])
						{
							maxval_index = 2;
						}
						nextLayer = dptr[maxval_index] >= cumulativeScoreQuery[x - 1] ? maxval_index : -1;
					}
					else nextLayer = -1; //but as long as the current layer is not -1, we should be able to push scores????
					break;
				case 2:
					nextLayer = (scoreMatrix[matrixIndex - 3] + ContinueDeletionScore > scoreMatrix[matrixIndex-4] + OpenDeletionScore) ? 2 : 1;
					break;
				}
				scoreStack.push(scoreMatrix[matrixIndex]);
				indexLength++;
				if (indexLength > 2000)
				{
					std::cerr << "debug error indexLength?\n";
					exit(1);
				}
				switch (currentLayer)
				{
				case 0: // insertion
					indexXStack.push(x);
					indexYStack.push(-1);
					x--;
					break;
				case 1: // diagonal move
					indexXStack.push(finishing1 ? -1 : x);
					indexYStack.push(finishing0 ? -1 : y);
					if (x == 0) // reached top edge
					{
						finishing1 = true;
						if (!complete1) // end the alignment here
						{
							nextLayer = -1;
						}
						y--;
					}
					else if (y == 0) // reached left edge
					{
						finishing0 = true;
						if (!complete0) // end the alignment here
						{
							nextLayer = -1;
						}
						x--;
					}
					else
					{
						x--;
						y--;
					}
					break;
				case 2: // deletion
					indexXStack.push(-1);
					indexYStack.push(y);
					y--;
					break;
				}
				currentLayer = nextLayer;
			}

			if (indexLength > arrayLength)
			{
				fprintf(stderr, "Error native trace back, array too short.\n");
				return 0;
			}

			for (int i = 0; i < indexLength; i++)
			{
				index1[i] = indexXStack.top();;
				indexXStack.pop();
				index2[i] = indexYStack.top();;
				indexYStack.pop();
				scores[i] = scoreStack.top();
				scoreStack.pop();
			}
			return indexLength;
		}


		///result matching Tom's vesion
		int TraceBack(int *index1, int *index2, double *scores, int arrayLength, bool complete0 = false, bool complete1 = true)
		{
			if (use_index_version)
			{
				return TraceBack_index(index1, index2, scores, arrayLength);
			}

			return TraceBack_tom(index1, index2, scores, arrayLength);
		}


		//version: index modified.
		//this traceback is for the index modifed version. not same as Tom's current version
		int TraceBack_index(int *index1, int *index2, double *scores, int arrayLength, bool complete0 = false, bool complete1 = true)
		{
			if (MaxRowIndex < 0 || MaxColIndex < 0)
			{
				arrayLength = 0;
				return 0;
			}

			// start traceback from far corner
			int x = MaxRowIndex;
			int y = MaxColIndex;
			int currentLayer = 1;

			// stop when the top row is passed
			stack<int> indexXStack;		//= new Stack<int>();
			stack<int> indexYStack;		// = new Stack<int>();
			stack<double> scoreStack; // = new Stack<double>();

			int dimReference = (int)kernel->GetReference()->Seq.size() + 1;
			int indexLength = 0;
			bool foundMatch = false;
			int nextLayer;
			bool finishing0 = false;
			bool finishing1 = false;
			while (currentLayer > -1)
			{
				int matrixIndex = x * dimReference * 3 + y * 3 + currentLayer;
				//std::cerr << "check point 2\n" << matrixIndex << " " << traceBackMatrix.size() << "\n";
				nextLayer = traceBackMatrix[matrixIndex];
				if (nextLayer == -1) //??? && (x == 0 || y == 0))
				{
					break;
				}
				scoreStack.push(scoreMatrix[matrixIndex]);
				indexLength++;
				switch (currentLayer)
				{
				case 0: // insertion
					indexXStack.push(x);
					indexYStack.push(-1);
					x--;
					break;
				case 1: // diagonal move
					indexXStack.push(finishing1 ? -1 : x);
					indexYStack.push(finishing0 ? -1 : y);
					if (x == 1) // reached top edge
					{
						finishing1 = true;
						if (!complete1) // end the alignment here
						{
							nextLayer = -1;
						}
						x--;
						y--;
					}
					else if (y == 1) // reached left edge
					{
						finishing0 = true;
						if (!complete0) // end the alignment here
						{
							nextLayer = -1;
						}
						x--;
						y--;
					}
					else
					{
						x--;
						y--;
					}
					break;
				case 2: // deletion
					indexXStack.push(-1);
					indexYStack.push(y);
					y -= 1;
					break;
				}
				currentLayer = nextLayer;
			}

			if (indexLength > arrayLength)
			{
				fprintf(stderr, "Error native trace back, array too short.\n");
				exit(1);
			}

			for (int i = 0; i < indexLength; i++)
			{
				int pos1 = indexXStack.top();
				index1[i] = pos1 == -1 ? -1 : pos1 - 1;
				indexXStack.pop();
				int pos2 = indexYStack.top();
				index2[i] = pos2 == -1 ? -1 : pos2 - 1;
				indexYStack.pop();
				scores[i] = scoreStack.top();
				scoreStack.pop();
			}
			return indexLength;
		}

		void SetKernelDistance(double distance)
		{

			MNMNAlignmentKernelAsymmetric *mnmn_kernel = dynamic_cast<MNMNAlignmentKernelAsymmetric *>(kernel.get());
			if (mnmn_kernel != nullptr)
			{
				mnmn_kernel->set_Distance(distance);
				std::cerr << "distance set to " << distance << "\n";
			}
			else
			{
				throw "incorrect model cast exception.\n";
			}
		}
	};



	/// <summary>
	/// Provides the means to perform pairwise alignment.
	/// </summary>
	/// <typeparam name="T">The nucleotide type for the first polynucleotide.</typeparam>
	/// <typeparam name="U">The nucleotide type for the second polynucleotide.</typeparam>
	template <class T>
	class PairwiseAlignerSymmetric
	{
	public:
		/// <summary>
		/// The total score for the alignment as reported.
		/// </summary>
		double AlignmentScore;

		bool Verbose = false;
	private:

		std::vector<double> scoreMatrix;
		std::vector<int> traceBackMatrix;

		std::vector<std::vector<double>> priorScore;

		int dimPhysical[2]; //{ -1, -1 };
		int dim[2];

		unique_ptr<AlignmentKernelSymmetric<T>> kernel;

		const double logScoreConst = log(0.25);

	public:
		/// <summary>
		/// Sets the default scores for the first column.
		/// </summary>
		/// <param name="value">The score to be accumulated over rows as the alignment scores that occur if the polynucleotide local alignment begins on a row other than the zeroth.</param>
		void SetPriorScore(int sequenceNumber)
		{
			if (dim[sequenceNumber] < 1)
				return;
			priorScore[sequenceNumber][0] = 0;
			for (int i = 1; i < dim[sequenceNumber]; i++)
			{
				priorScore[sequenceNumber][i] = kernel->PriorScore(sequenceNumber, i - 1);
			}
		}

		/// <summary>
		/// Sets the default scores for the first column.
		/// </summary>
		/// <param name="value">The score to be accumulated over rows as the alignment scores that occur if the polynucleotide local alignment begins on a row other than the zeroth.</param>
		void SetPriorScore(int iSeq, double score)
		{
			if (dim[iSeq] < 1)
				return;
			priorScore[iSeq][0] = 0;
			for (int i = 1; i < dim[iSeq]; i++)
			{
				priorScore[iSeq][i] = score;
			}
		}

		/// <summary>
		/// Sets the default scores for the first row.
		/// </summary>
		/// <param name="scores">A vector of alignment scores taken as the alignment scores that occur if the polynucleotide local alignment begins on a column other than the zeroth.</param>
		void SetPriorScore(int iSeq, vector<double>& scores)
		{
			if (priorScore[iSeq].size() != scores.size())
			{
				throw std::runtime_error("Scores must have the same length as the existing prior score array.");
			}
			priorScore[iSeq] = scores;
		}


		//this constructor exists as a helper to create aligner.
		//without this, visual studio give error when trying to initalize an aligner
		PairwiseAlignerSymmetric()
			: PairwiseAlignerSymmetric(std::unique_ptr<PMFAlignmentKernelSymmetric>(new PMFAlignmentKernelSymmetric()))
		{
		}

		/// <summary>
		/// Constructor.
		/// </summary>
		/// <param name="_pair">The polynucleotide pair to be aligned.</param>
		/// <param name="_comparator">A comparator of the appropriate types for use in the alignment.</param>
		PairwiseAlignerSymmetric(unique_ptr<AlignmentKernelSymmetric<T>> _kernel)
			:kernel(std::move(_kernel))
		{
			dim[0] = dim[1] = -1;

			for (int iSeq = 0; iSeq < 2; iSeq++)
			{
				if (!kernel->GetSeq(iSeq))
				{
					dim[iSeq] = dimPhysical[iSeq] = 0;
				}
				else
				{
					dim[iSeq] = dimPhysical[iSeq] = (int)kernel->GetSeq(iSeq)->Seq.size() + 1;
				}
			}

			// allocate matrices
			scoreMatrix.resize(dimPhysical[0] * dimPhysical[1] * 3, 0);
			traceBackMatrix.resize(dimPhysical[0] * dimPhysical[1] * 3);
			priorScore.resize(2);
			for (int iSeq = 0; iSeq < 2; iSeq++)
			{
				priorScore[iSeq].resize(dimPhysical[iSeq], 0);
				SetPriorScore(iSeq);
			}
		}

		/// <summary>
		/// Constructor.
		/// </summary>
		/// <param name="kernel">A comparator of the appropriate types for use in the alignment.</param>
		/// <param name="_dim0">The initial first dimension of the score and traceback matrices.</param>
		/// <param name="_dim1">The initial second dimension of the score and traceback matrices.</param>
		PairwiseAlignerSymmetric(int dim0, int dim1, unique_ptr<AlignmentKernelSymmetric<T>> kernel)
			: kernel(kernel)
		{
			this.dim[0] = dimPhysical[0] = dim0 + 1;
			this.dim[1] = dimPhysical[1] = dim1 + 1;

			// allocate matrices
			scoreMatrix.resize(dimPhysical[0] * dimPhysical[1] * 3, 0);
			traceBackMatrix.resize(dimPhysical[0] * dimPhysical[1] * 3, 0);
			priorScore.resize(2, vector<double>(dimPhysical[0]));
		}

		void SetSequence(shared_ptr<Polymer<T>>& p, int sequenceNumber)
		{
			//the sequence might be same during repeated call...
			if (kernel->GetSeq(sequenceNumber) == p)
			{
				return;
			}
			kernel->SetSeq(p, sequenceNumber);
			dim[sequenceNumber] = (int)p->Seq.size() + 1;
			if (dim[sequenceNumber] > dimPhysical[sequenceNumber])
			{
				dimPhysical[sequenceNumber] = dim[sequenceNumber];

				scoreMatrix.resize(dimPhysical[0] * dimPhysical[1] * 3, 0);
				traceBackMatrix.resize(dimPhysical[0] * dimPhysical[1] * 3, 0);
				priorScore[sequenceNumber].resize(dimPhysical[sequenceNumber]);
			}
			SetPriorScore(sequenceNumber);
		}

		/// <summary>
		/// Fills the score matrix for pairwise alignment between the two Polynucleotides in the pair.
		/// </summary>

		/// <summary>
		/// Fills the score matrix for pairwise alignment between the two Polynucleotides in the pair.
		/// note: this is the translation of FillScoreMatrix2 not FillScoreMatrix
		/// </summary>
		vector<int> FillScoreMatrix()
		{
			dim[0] = (int)kernel->GetSeq(0)->Seq.size();
			dim[1] = (int)kernel->GetSeq(1)->Seq.size();
			int dim2 = dim[1] + 1; //for the extra column in the left

			int i0 = 0;
			int i1 = 0;

			double scoreDiag, scoreDiagMax;
			double scoreDelete, scoreDeleteMax;
			double scoreInsert, scoreInsertMax;

			// compute element [-1, -1]
			scoreMatrix[MATRIX_INDEX(0, 0, 0)] = 0;
			traceBackMatrix[MATRIX_INDEX(0, 0, 0)] = -1;
			scoreMatrix[MATRIX_INDEX(0, 0, 1)] = 0;
			traceBackMatrix[MATRIX_INDEX(0, 0, 1)] = -1;
			scoreMatrix[MATRIX_INDEX(0, 0, 2)] = 0;
			traceBackMatrix[MATRIX_INDEX(0, 0, 2)] = -1;

			double matchCost = std::numeric_limits<double>::lowest();

			// compute row -1
			for (i1 = 1; i1 <= dim[1]; i1++)
			{
				scoreMatrix[MATRIX_INDEX(0, i1, 0)] = scoreMatrix[MATRIX_INDEX(0, i1 - 1, 2)] + priorScore[1][i1]; //std::numeric_limits<double>::lowest();
				traceBackMatrix[MATRIX_INDEX(0, i1, 0)] = 2;
				scoreMatrix[MATRIX_INDEX(0, i1, 1)] = std::numeric_limits<double>::lowest();
				traceBackMatrix[MATRIX_INDEX(0, i1, 1)] = 2;
				scoreMatrix[MATRIX_INDEX(0, i1, 2)] = scoreMatrix[MATRIX_INDEX(0, i1 - 1, 2)] + priorScore[1][i1];
				traceBackMatrix[MATRIX_INDEX(0, i1, 2)] = 2;
			}

			// compute column -1
			for (i0 = 1; i0 <= dim[0]; i0++)
			{
				scoreMatrix[MATRIX_INDEX(i0, 0, 0)] = scoreMatrix[MATRIX_INDEX(i0 - 1, 0, 0)] + priorScore[0][i0];
				traceBackMatrix[MATRIX_INDEX(i0, 0, 0)] = 0;
				scoreMatrix[MATRIX_INDEX(i0, 0, 1)] = std::numeric_limits<double>::lowest();
				traceBackMatrix[MATRIX_INDEX(i0, 0, 1)] = 0;
				scoreMatrix[MATRIX_INDEX(i0, 0, 2)] = scoreMatrix[MATRIX_INDEX(i0 - 1, 0, 0)] + priorScore[0][i0]; // std::numeric_limits<double>::lowest();
				traceBackMatrix[MATRIX_INDEX(i0, 0, 2)] = 0;
			}

			//testing uniform comparison scores.
			double OpenGapScoreVal = kernel->OpenInsertionScore(0, 0);
			double ContinueGapScoreVal = kernel->ContinueGapScore(0, 0);
			double ComparisionScoreMin = OpenGapScoreVal;
			// compute entries up to last row and column
			for (i0 = 1; i0 < dim[0]; i0++)
			{
				for (i1 = 1; i1 < dim[1]; i1++)
				{
					// layer 0: insertions into sequence 0
					scoreInsertMax = std::numeric_limits<double>::lowest();;
					for (int iLayer = 0; iLayer < 2; iLayer++)
					{
						scoreInsert = scoreMatrix[MATRIX_INDEX(i0 - 1, i1, iLayer)] + (iLayer == 0 ? kernel->ContinueGapScore(i0 - 1, i1 - 1) : kernel->OpenInsertionScore(i0 - 1, i1 - 1));
						if (scoreInsert > scoreInsertMax)
						{
							scoreInsertMax = scoreInsert;
							traceBackMatrix[MATRIX_INDEX(i0, i1, 0)] = iLayer;
						}
						scoreMatrix[MATRIX_INDEX(i0, i1, 0)] = scoreInsertMax;
					}

					// layer 2: deletions from sequence 0
					scoreDeleteMax = std::numeric_limits<double>::lowest();;
					for (int iLayer = 1; iLayer < 3; iLayer++)
					{
						scoreDelete = scoreMatrix[MATRIX_INDEX(i0, i1 - 1, iLayer)] + (iLayer == 2 ? kernel->ContinueGapScore(i0 - 1, i1 - 1) : kernel->OpenDeletionScore(i0 - 1, i1 - 1));
						if (scoreDelete > scoreDeleteMax)
						{
							scoreDeleteMax = scoreDelete;
							traceBackMatrix[MATRIX_INDEX(i0, i1, 2)] = iLayer;
						}
						scoreMatrix[MATRIX_INDEX(i0, i1, 2)] = scoreDeleteMax;
					}


					// layer 1: matches
					scoreDiagMax = std::numeric_limits<double>::lowest();;
					matchCost = max(ComparisionScoreMin, kernel->ComparisonScore(i0 - 1, i1 - 1));
					for (int iLayer = 0; iLayer < 3; iLayer++)
					{
						scoreDiag = scoreMatrix[MATRIX_INDEX(i0 - 1, i1 - 1, iLayer)] + matchCost;
						if (scoreDiag > scoreDiagMax)
						{
							scoreDiagMax = scoreDiag;
							traceBackMatrix[MATRIX_INDEX(i0, i1, 1)] = iLayer;
						}
					}
					scoreMatrix[MATRIX_INDEX(i0, i1, 1)] = scoreDiagMax;
				}

			}

			// last row
			i0 = dim[0];
			for (i1 = 1; i1 < dim[1]; i1++)
			{
				// layer 0: insertions into sequence 0
				scoreInsertMax = std::numeric_limits<double>::lowest();;
				for (int iLayer = 0; iLayer < 2; iLayer++)
				{
					double a = scoreMatrix[MATRIX_INDEX(i0 - 1, i1, iLayer)];
					double b = (iLayer == 0 ? kernel->ContinueGapScore(i0 - 1, i1 - 1) : kernel->OpenInsertionScore(i0 - 1, i1 - 1));
					scoreInsert = scoreMatrix[MATRIX_INDEX(i0 - 1, i1, iLayer)] + (iLayer == 0 ? kernel->ContinueGapScore(i0 - 1, i1 - 1) : kernel->OpenInsertionScore(i0 - 1, i1 - 1));
					if (scoreInsert > scoreInsertMax)
					{
						scoreInsertMax = scoreInsert;
						traceBackMatrix[MATRIX_INDEX(i0, i1, 0)] = iLayer;
					}
					scoreMatrix[MATRIX_INDEX(i0, i1, 0)] = scoreInsertMax;
				}

				// layer 2: deletions from sequence 0
				scoreDeleteMax = std::numeric_limits<double>::lowest();;
				for (int iLayer = 1; iLayer < 3; iLayer++)
				{
					scoreDelete = scoreMatrix[MATRIX_INDEX(i0, i1 - 1, iLayer)] + priorScore[1][i1];
					if (scoreDelete > scoreDeleteMax)
					{
						scoreDeleteMax = scoreDelete;
						traceBackMatrix[MATRIX_INDEX(i0, i1, 2)] = iLayer;
					}
					scoreMatrix[MATRIX_INDEX(i0, i1, 2)] = scoreDeleteMax;
				}

				// layer 1: matches
				scoreDiagMax = std::numeric_limits<double>::lowest();;
				matchCost = max(ComparisionScoreMin, kernel->ComparisonScore(i0 - 1, i1 - 1));
				for (int iLayer = 0; iLayer < 3; iLayer++)
				{
					scoreDiag = scoreMatrix[MATRIX_INDEX(i0 - 1, i1 - 1, iLayer)] + matchCost;
					if (scoreDiag > scoreDiagMax)
					{
						scoreDiagMax = scoreDiag;
						traceBackMatrix[MATRIX_INDEX(i0, i1, 1)] = iLayer;
					}
				}
				scoreMatrix[MATRIX_INDEX(i0, i1, 1)] = scoreDiagMax;
			}

			// last column
			i1 = dim[1];
			for (i0 = 1; i0 < dim[0]; i0++)
			{
				// layer 0: insertions into sequence 0
				scoreInsertMax = std::numeric_limits<double>::lowest();;
				for (int iLayer = 0; iLayer < 2; iLayer++)
				{
					scoreInsert = scoreMatrix[MATRIX_INDEX(i0 - 1, i1, iLayer)] + priorScore[0][i0];
					if (scoreInsert > scoreInsertMax)
					{
						scoreInsertMax = scoreInsert;
						traceBackMatrix[MATRIX_INDEX(i0, i1, 0)] = iLayer;
					}
					scoreMatrix[MATRIX_INDEX(i0, i1, 0)] = scoreInsertMax;
				}

				// layer 2: deletions from sequence 0
				scoreDeleteMax = std::numeric_limits<double>::lowest();;
				for (int iLayer = 1; iLayer < 3; iLayer++)
				{
					scoreDelete = scoreMatrix[MATRIX_INDEX(i0, i1 - 1, iLayer)] + (iLayer == 2 ? kernel->ContinueGapScore(i0 - 1, i1 - 1) : kernel->OpenDeletionScore(i0 - 1, i1 - 1));
					if (scoreDelete > scoreDeleteMax)
					{
						scoreDeleteMax = scoreDelete;
						traceBackMatrix[MATRIX_INDEX(i0, i1, 2)] = iLayer;
					}
					scoreMatrix[MATRIX_INDEX(i0, i1, 2)] = scoreDeleteMax;
				}

				// layer 1: matches
				scoreDiagMax = std::numeric_limits<double>::lowest();;
				matchCost = max(ComparisionScoreMin, kernel->ComparisonScore(i0 - 1, i1 - 1));
				for (int iLayer = 0; iLayer < 3; iLayer++)
				{
					scoreDiag = scoreMatrix[MATRIX_INDEX(i0 - 1, i1 - 1, iLayer)] + matchCost;
					if (scoreDiag > scoreDiagMax)
					{
						scoreDiagMax = scoreDiag;
						traceBackMatrix[MATRIX_INDEX(i0, i1, 1)] = iLayer;
					}
				}
				scoreMatrix[MATRIX_INDEX(i0, i1, 1)] = scoreDiagMax;
			}

			// last entry
			i0 = dim[0];
			i1 = dim[1];
			// layer 0: insertions into sequence 0
			scoreInsertMax = std::numeric_limits<double>::lowest();;
			for (int iLayer = 0; iLayer < 2; iLayer++)
			{
				scoreInsert = scoreMatrix[MATRIX_INDEX(i0 - 1, i1, iLayer)] + priorScore[0][i0];
				if (scoreInsert > scoreInsertMax)
				{
					scoreInsertMax = scoreInsert;
					traceBackMatrix[MATRIX_INDEX(i0, i1, 0)] = iLayer;
				}
				scoreMatrix[MATRIX_INDEX(i0, i1, 0)] = scoreInsertMax;
			}

			// layer 2: deletions from sequence 0
			scoreDeleteMax = std::numeric_limits<double>::lowest();;
			for (int iLayer = 1; iLayer < 3; iLayer++)
			{
				scoreDelete = scoreMatrix[MATRIX_INDEX(i0, i1 - 1, iLayer)] + priorScore[1][i1];
				if (scoreDelete > scoreDeleteMax)
				{
					scoreDeleteMax = scoreDelete;
					traceBackMatrix[MATRIX_INDEX(i0, i1, 2)] = iLayer;
				}
				scoreMatrix[MATRIX_INDEX(i0, i1, 2)] = scoreDeleteMax;
			}

			// layer 1: matches
			scoreDiagMax = std::numeric_limits<double>::lowest();;
			matchCost = max(ComparisionScoreMin, kernel->ComparisonScore(i0 - 1, i1 - 1));
			for (int iLayer = 0; iLayer < 3; iLayer++)
			{
				scoreDiag = scoreMatrix[MATRIX_INDEX(i0 - 1, i1 - 1, iLayer)] + matchCost;
				if (scoreDiag > scoreDiagMax)
				{
					scoreDiagMax = scoreDiag;
					traceBackMatrix[MATRIX_INDEX(i0, i1, 1)] = iLayer;
				}
			}

			scoreMatrix[MATRIX_INDEX(i0, i1, 1)] = scoreDiagMax;

			return GetTracebackStart();
		}

		friend class TestHelper;

	private:
		vector<int> GetTracebackStart()
		{
			// find starting point and alignment score for traceback

			int lastQIndex = dim[0];
			int lastRIndex = dim[1];
			int dim2 = dim[1] + 1; //total number of columns
			int startLayer = 0;
			double maxScore = scoreMatrix[MATRIX_INDEX(lastQIndex, lastRIndex, 0)];
			for (int layer = 1; layer < 3; layer++)
			{
				if (scoreMatrix[MATRIX_INDEX(lastQIndex, lastRIndex, layer)] > maxScore)
				{
					startLayer = layer;
					maxScore = scoreMatrix[MATRIX_INDEX(lastQIndex, lastRIndex, layer)];
				}
			}
			AlignmentScore = maxScore;
			return vector<int>{ lastQIndex, lastRIndex, startLayer };
		}

	public:

		/// <summary>
		/// Performs traceback for a dynamic programming-based pairwise alignment.
		/// </summary>
		PolymerPair<T, T> TraceBack(vector<int>& start)
		{
			// traceback either to the beginning, indicated by -1 as the next layer.

			PolymerPair<T, T> pair(kernel->GetSeq(0), kernel->GetSeq(1));

			if (start[0] < 0 || start[1] < 0)
			{
				return pair;
			}

			int x = start[0];
			int y = start[1];
			int currentLayer = 1;
			if (start.size() > 2)
			{
				currentLayer = start[2];
			}

			stack<int> indexXstack;
			stack<int> indexYstack;
			stack<double> scorestack;

			int dim2 = this->dim[1] + 1;
			int indexLength = 0;
			int nextLayer;
			while (x > 0 || y > 0)
			{
				nextLayer = traceBackMatrix[MATRIX_INDEX(x, y, currentLayer)];
				scorestack.push(scoreMatrix[MATRIX_INDEX(x, y, currentLayer)]);
				indexLength++;
				switch (currentLayer)
				{
				case 0: // insertion
					indexXstack.push(x - 1);
					indexYstack.push(-1);
					x--;
					break;
				case 1: // diagonal move
					indexXstack.push(x - 1);
					indexYstack.push(y - 1);
					x--;
					y--;
					break;
				case 2: // deletion
					indexXstack.push(-1);
					indexYstack.push(y - 1);
					y--;
					break;
				}
				currentLayer = nextLayer;
				//fix the corner transition
				if (x == 0)currentLayer = 2;
				if (y == 0)currentLayer = 0;
			}

			// declare the index array and fill it 
			pair.Index.resize(indexLength, vector<int>(2));
			pair.Score.reserve(indexLength);
			for (int i = 0; i < indexLength; i++)
			{

				pair.Index[i][0] = indexXstack.top();
				indexXstack.pop();

				pair.Index[i][1] = indexYstack.top();
				indexYstack.pop();

				pair.Score.push_back(scorestack.top());
				scorestack.pop();
			}

			//adjust score.
			double score = this->AlignmentScore;
			double ComparisionScoreMin = kernel->OpenInsertionScore(0, 0);

			//the open gap score used during alingment.
			double OpenGapScoreVal = kernel->OpenInsertionScore(0, 0);

			//class PMFAlignmentKernelSymmetric : public AlignmentKernelSymmetric < vector<double> >
			auto pmfKenal = dynamic_cast<PMFAlignmentKernelSymmetric*>(kernel.get());


			int lastQIndex = dim[0] - 1;
			int lastRIndex = dim[1] - 1;
			int indexQ = -1;
			int indexR = -1;
			int gap_count = 0;
			for (int i = 0; i < pair.Index.size(); i++)
			{
				auto& indexi = pair.Index[i];

				if (indexi[0] != -1)indexQ++;
				if (indexi[1] != -1)indexR++;
				if (indexi[0] == -1 || indexi[1] == -1)
				{
					if ((indexi[0] == -1 && indexi[1] == -1) || i == 0)continue;

					//check if its gap in the  middle
					if (indexQ == -1 || indexQ == lastQIndex || indexR == -1 || indexR == lastRIndex)continue;
					gap_count++;
					if (indexi[0] == -1)
					{
						//only for gap open.
						if (pair.Index[i - 1][0] != -1)
						{
							double gapScore = pmfKenal->OpenInsertionScoreBounded(indexQ, indexR);
							if (gapScore < OpenGapScoreVal)
							{
								score -= (OpenGapScoreVal - gapScore);
							}
						}
					}
					else //indexi[1] == -1
					{
						//only for gap open.
						if (pair.Index[i - 1][1] != -1)
						{
							double gapScore = pmfKenal->OpenDeletionScoreBounded(indexQ, indexR);
							if (gapScore < OpenGapScoreVal)
							{
								score -= (OpenGapScoreVal - gapScore);
							}
						}
					}
					continue;
				}
				double comparisonScore = kernel->ComparisonScore(indexi[0], indexi[1]);
				if (comparisonScore >= ComparisionScoreMin)continue;
				score -= (ComparisionScoreMin - comparisonScore);
			}
			pair.AlignmentScore = score;

			//adjust score
			int alignLen = (int)pair.Index.size();
			int sumLen = (int)pair.Polymer0->Seq.size() + (int)pair.Polymer1->Seq.size();

			//enforce at least 20 bp match
			if (sumLen - alignLen < AlignMinMatch)
			{
				pair.AlignmentScore = 0; //considered no valid alignment
				return pair;
			}
			auto score_adjust = -sumLen * logScoreConst;
			pair.AlignmentScore += score_adjust;

			if (gap_count > 0)
			{
				if (pair.AlignmentScore > 0)
				{
					//std::cout << "accecpted gap count= " << gap_count << "\n";
				}
				else
				{
					//std::cout << "rejected gap count = " << gap_count << "\n";
				}
			}

			return pair;
		}
	};

}



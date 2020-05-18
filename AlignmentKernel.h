#pragma once

#include <math.h>
#include <memory>
#include <vector>
#include <limits>


#include "typedef.h"
#include "EvolutionModels.h"
#include "Polymer.h"
#include "PolymerPair.h"
#include "GeneReplicator.h"
#include "TypeConverter.h"
#include "alignedAllocator.h"

using namespace std;

namespace CLPLib
{

	class GapProbabilities;

	/// <summary>
	/// Data structure to contain gapscore components.
	/// </summary>
	class GapScores
	{
	public:
		/// <summary>
		/// The score for opening a new deletion.
		/// </summary>
		double OpenDeletionScore = -3.912; // log(0.02), 0.02 = GapProbabilities.OpenDeletionProbability default value.
										   /// <summary>
										   /// The score for continuing a deletion for one more base.
										   /// </summary>
		double ContinueDeletionScore = -0.1054; // log(0.9)
											   /// <summary>
											   /// The score for opening a new insertion.
											   /// </summary>
		double OpenInsertionScore = -3.912;
		/// <summary>
		/// The score for continuing an insertion for one more base.
		/// </summary>
		double ContinueInsertionScore = -0.1054;

		GapScores() // uses default values
		{

		}

		//method not in managed code, used for pass data from managed code
		GapScores(double ods, double cds, double ois, double cis)
		{
			this->OpenDeletionScore = ods;
			this->ContinueDeletionScore = cds;
			this->OpenInsertionScore = ois;
			this->ContinueInsertionScore = cis;
		}



		GapScores(GapProbabilities& gapProbabilities);
	};

	class GapProbabilities
	{
	public:
		/// <summary>
		/// The probability for opening a new deletion.
		/// </summary>
		double OpenDeletionProbability = 0.02;
		/// <summary>
		/// The probability for continuing a deletion for one more base.
		/// </summary>
		double ContinueDeletionProbability = 0.9;
		/// <summary>
		/// The probability for opening a new insertion.
		/// </summary>
		double OpenInsertionProbability = 0.02;
		/// <summary>
		/// The probability for continuing an insertion for one more base.
		/// </summary>
		double ContinueInsertionProbability = 0.9;
	};


	/// <summary>
	/// Abstract class that specifies what methods alignment kernel must provide. It may compare Polymers of different nucleotide encoding types.
	/// </summary>
	/// <typeparam name="T">The type of the first nucleotide encoding.</typeparam>
	/// <typeparam name="U">The type of the second nucleotide encoding.</typeparam>

	template <class T, class U>
	class AlignmentKernelAsymmetric
	{
	public:

		/// <summary>
		/// Reports the status of the symbol as a gap.
		/// </summary>
		/// <param name="n">The symbol to be tested.</param>
		/// <returns>True if the symbol represents a gap.</returns>
		virtual bool IsGap(U& n) = 0;
		/// <summary>
		/// Computes the comparison score for the alignment of the nucleotide at a given position in the query with the 
		/// nucleotide at a given position in the reference.
		/// </summary>
		/// <param name="posInQuery">The position in query.</param>
		/// <param name="posInReference">The position in the reference.</param>
		/// <returns></returns>
		virtual double ComparisonScore(int posInQuery, int posInReference) = 0;
		/// <summary>
		/// Computes the score for opening an insertion in the query at a given query position when that position 
		/// is aligned against a given position in the reference.
		/// </summary>
		/// <param name="posInQuery">The query position.</param>
		/// <param name="posInReference">The reference position.</param>
		/// <returns>The cost of opening an insertion at the given coordinates.</returns>
		virtual double OpenInsertionScore(int posInQuery, int posInReference) = 0;
		/// <summary>
		/// Computes the score for continuing an insertion in the query at a given query position when that position 
		/// is aligned against a given position in the reference.
		/// </summary>
		/// <param name="posInQuery">The query position.</param>
		/// <param name="posInReference">The reference position.</param>
		/// <returns>The cost of opening an insertion at the given coordinates.</returns>
		virtual double ContinueInsertionScore(int posInQuery, int posInReference) = 0;
		/// <summary>
		/// Computes the score for opening a deletion in the query at a given query position when that position 
		/// is aligned against a given position in the reference.
		/// </summary>
		/// <param name="posInQuery">The query position.</param>
		/// <param name="posInReference">The reference position.</param>
		/// <returns>The cost of opening an insertion at the given coordinates.</returns>
		virtual double OpenDeletionScore(int posInQuery, int posInReference) = 0;
		/// <summary>
		/// Computes the score for continuing a deletion in the query at a given query position when that position 
		/// is aligned against a given position in the reference.
		/// </summary>
		/// <param name="posInQuery">The query position.</param>
		/// <param name="posInReference">The reference position.</param>
		/// <returns>The cost of opening an insertion at the given coordinates.</returns>
		virtual double ContinueDeletionScore(int posInQuery, int posInReference) = 0;
		/// <summary>
		/// Sets the query Polymer.
		/// </summary>
		/// <param name="query">The Polymer to set as query.</param>
		virtual void SetQuery(shared_ptr<Polymer<T>> query) = 0;
		/// <summary>
		/// Gets the query.
		/// </summary>
		/// <returns>The query Polymer.</returns>
		virtual shared_ptr<Polymer<T>> GetQuery() = 0;
		/// <summary>
		/// Sets the reference Polymer.
		/// </summary>
		/// <returns>The Polymer to set as reference.</returns>
		virtual void SetReference(shared_ptr<Polymer<U>> reference) = 0;
		/// <summary>
		/// Gets the reference.
		/// </summary>
		/// <returns>The reference Polymer.</returns>
		virtual shared_ptr<Polymer<U>> GetReference() = 0;


		/// <summary>
		/// The score taken as the default for the query at the given position.
		/// </summary>
		/// <param name="posInQuery">The position in the query.</param>
		/// <returns>The default score.</returns>
		virtual double PriorQueryScore(int posInQuery) = 0;
		/// <summary>
		/// The score taken as default for the query.
		/// </summary>
		/// <returns>The default query score.</returns>
		virtual double PriorQueryScore() = 0;
		/// <summary>
		/// Gets the array of default query scores.
		/// </summary>
		/// <returns>An array of doubles giving the default query score by position.</returns>
		virtual vector<double> PriorQueryScoreArray() = 0;

		/// <summary>
		/// fill the compareScores array with indel score in the following order
		/// matchScore,
		/// OpenInsertionScore, ContinueInsertionScore
		/// OpenDeletionScore, ContinueDeletionScore
		/// for each [i, j] where i = 0..querySeq.Length and j = 0..reference.Length
		/// </summary>
		/// <param name="?"></param>

		//virtual void GetComparisonScore(vector<double>& compareScores) = 0;

		virtual void print_debug_info()
		{

		}
	};



	//note: this is the class in AlignmentKernel.cs of LPLib and intended to match the class of MNMNAlignmentKernelAsymmetric2
	//since that is the class that is used in Cloanalyst.
	class MNMNAlignmentKernelAsymmetric : public AlignmentKernelAsymmetric < byte, byte >
	{
	private:
		std::shared_ptr<NucleotideEvolutionModel> model;

		std::shared_ptr<Polymer<byte>> query;
		std::shared_ptr<Polymer<byte>> reference;

		double distance;

		std::vector<double> mutationRate;

	public:

		std::shared_ptr<GapScores> gapScores;

		bool local = false;

		double get_Distance()
		{
			return distance;
		}

		void set_Distance(double value)
		{
			distance = value;
			model->set_Time(distance);
			//recompute score array
			computeComparisonScoreArray();
			mutationCost = log(1 - exp(-distance));
		}

		double mutationCost;

		vector<vector<double, AlignmentAllocator<double, 64>>> comparisonScore;

		//smaller verson, using PureNucleotides, we are putting each element
		//in a 64 aligned memory so that it can fit one cache line.
		vector<vector<double, AlignmentAllocator<double, 64>>> pureComparisonScore;


		//obsolate, used to compute reference matching score.
		vector<vector<double>> comparisonScoreArray;


		MNMNAlignmentKernelAsymmetric(std::shared_ptr<NucleotideEvolutionModel> _model) : model(_model), gapScores(new GapScores())
		{
			this->distance = model->get_Time();
			mutationCost = log(1 - exp(-distance));
			computeComparisonScoreArray();
		}

		MNMNAlignmentKernelAsymmetric(std::shared_ptr<NucleotideEvolutionModel>& model, shared_ptr<Polymer<byte>> query, shared_ptr<Polymer<byte>> reference)
			:model(model), query(query), reference(reference), gapScores(new GapScores())
		{
			this->distance = model->get_Time();
			mutationCost = log(1 - exp(-distance));
			computeComparisonScoreArray();
		}

		MNMNAlignmentKernelAsymmetric(std::shared_ptr<NucleotideEvolutionModel>& _model, PolymerPair<byte, byte> &pair)
			:model(_model), query(pair.Polymer0), reference(pair.Polymer1), gapScores(new GapScores())
		{
			this->distance = model->get_Time();
			computeComparisonScoreArray();
		}

		void computeComparisonScoreArray()
		{
			comparisonScore.resize(MixedNucleotideByte::NSymbols, 
				vector<double, AlignmentAllocator<double, 64>>(MixedNucleotideByte::NSymbols, 0));
			for (byte bR = 0; bR < MixedNucleotideByte::NSymbols; bR++)
			{
				vector<double> pR = TypeConverter::ToPMF(bR);
				for (byte bQ = 0; bQ < MixedNucleotideByte::NSymbols; bQ++)
				{
					vector<double> pQ = TypeConverter::ToPMF(bQ);
					double x = 0;
					for (int r = 0; r < 4; r++) // pure nucleotide states in reference sequence
					{
						for (int q = 0; q < 4; q++) // pure nuclelotide states in query sequence
						{
							x += pR[r] * pQ[q] * model->TransitionProbability(r, q);
						}
					}
					comparisonScore[bQ][bR] = log(x);
				}
			}

			//for shorter versions - acgtn
			pureComparisonScore.resize(PureNucleotideInt::symbols.size(), 
				vector<double, AlignmentAllocator<double, 64>>(PureNucleotideInt::symbols.size(), 0));
			std::vector<int> tmap{ 1, 2, 4, 8, 15 }; //from pure to mixed map
			for (int i = 0; i < PureNucleotideInt::symbols.size(); i++)
			{
				int mi = tmap[i];
				for (int j = 0; j < PureNucleotideInt::symbols.size(); j++)
				{
					int mj = tmap[j];
					pureComparisonScore[i][j] = comparisonScore[mi][mj];
				}
			}

		}

		//may only need to do for acgt???
		void computeComparisonScoreOfReference()
		{
			comparisonScoreArray.resize(MixedNucleotideByte::NSymbols);
			if (comparisonScoreArray[0].size() < reference->Seq.size())
			{
				for (int i = 0; i < comparisonScoreArray.size(); i++)
				{
					comparisonScoreArray[i].resize(reference->Seq.size());
				}
			}

			for (byte i = 0; i < MixedNucleotideByte::NSymbols; i++)
			{
				auto& sai = comparisonScoreArray[i];
				auto& seq = reference->Seq;
				auto& csi = comparisonScore[i];
				for (int j = 0; j < seq.size(); j++)
				{
					sai[j] = csi[seq[j]];
				}
			}
		}

		bool IsGap(byte& n) override
		{
			return n == 16;
		}


		double ComparisonScore(int posInQuery, int posInReference) override
		{

			auto i = query->Seq[posInQuery];
			auto j = reference->Seq[posInReference];
			return comparisonScore[i][j];
		}

		//return precomputed comparescore for the given query position
		std::vector<double>& ComparisonScoreQuery(int posInQuery)
		{
			auto i = query->Seq[posInQuery];
			return comparisonScoreArray[i];
		}

		std::vector<double> ComparisonScoreReference(int posInReference)
		{
			std::vector<double> scorevec(query->Seq.size());
			auto r = reference->Seq[posInReference];
			for (int i = 0; i < query->Seq.size(); i++)
			{
				auto q = query->Seq[i];
				scorevec[i] = comparisonScore[q][r];
			}
			return scorevec;
		}


		double OpenInsertionScore(int posInQuery, int posInReference)override
		{
			return mutationCost + gapScores->OpenInsertionScore + PureNucleotideInt::LogPrior;
		}

		double ContinueInsertionScore(int posInQuery, int posInReference) override
		{
			return gapScores->ContinueInsertionScore + PureNucleotideInt::LogPrior;
		}

		double OpenDeletionScore(int posInQuery, int posInReference)override
		{
			return mutationCost + gapScores->OpenDeletionScore;
		}

		double ContinueDeletionScore(int posInQuery, int posInReference)override
		{
			return gapScores->ContinueDeletionScore;
		}

		void SetQuery(shared_ptr<Polymer<byte>> _query) override
		{
			query = _query;
		}

		shared_ptr<Polymer<byte>> GetQuery() override
		{
			return query;
		}

	public:

		void SetReference(shared_ptr<Polymer<byte>> reference)override
		{
			this->reference = reference;
		}

		shared_ptr<Polymer<byte>> GetReference()override
		{
			return reference;
		}

		double PriorQueryScore(int posInQuery)override
		{
			return PureNucleotideInt::LogPrior;
		}

		double PriorQueryScore() override
		{
			return PureNucleotideInt::LogPrior;
		}

		vector<double> PriorQueryScoreArray() override
		{
			vector<double> returnValue(query->Seq.size(), 0);

			returnValue[0] = PureNucleotideInt::LogPrior;
			for (unsigned int i = 1; i < query->Seq.size(); i++)
			{
				returnValue[i] = PureNucleotideInt::LogPrior + returnValue[i - 1];
			}
			return returnValue;
		}

		vector<vector<double, AlignmentAllocator<double, 64>>>& getComparisonMatrix()
		{
			return comparisonScore;
		}

		vector<vector<double, AlignmentAllocator<double, 64>>>& getPureComparisonMatrix()
		{
			return pureComparisonScore;
		}


		void print_debug_info() override
		{

			std::cerr << "distance: " << distance << "\n";
		}
	};

	class PMFAlignmentKernel : public AlignmentKernelAsymmetric < vector<double>, vector<double> >
	{
	private:
		shared_ptr<Polymer<vector<double>>> query;
		shared_ptr<Polymer<vector<double>>> reference;
		shared_ptr<GapScores> gapScores;
		double mutationCost = MixedNucleotideByte::LogPrior;

	public:
		PMFAlignmentKernel(shared_ptr<GapScores> gapScores) : gapScores(gapScores)
		{
		}

		PMFAlignmentKernel(shared_ptr<Polymer<vector<double>>> query, shared_ptr<Polymer<vector<double>>> reference)\
			: query(query), reference(reference), gapScores(new GapScores())
		{
		}

		PMFAlignmentKernel(PolymerPair<vector<double>, vector<double>> pair)
			:query(pair.Polymer0), reference(pair.Polymer1), gapScores(new GapScores())
		{
		}

		bool IsGap(vector<double>& n) override
		{
			return n[4] == 1;
		}

		double ComparisonScore(int posInQuery, int posInReference) override
		{
			return MixedNucleotideByte::LogPrior + log(NucleotidePMF::IDFlux((query->Seq)[posInQuery], (reference->Seq)[posInReference]));
		}


		double OpenInsertionScore(int posInQuery, int posInReference)override
		{
			return mutationCost + gapScores->OpenInsertionScore + 2 * PureNucleotideInt::LogPrior;
		}

		double ContinueInsertionScore(int posInQuery, int posInReference)override
		{
			return gapScores->ContinueInsertionScore + 2 * PureNucleotideInt::LogPrior;
		}

		double OpenDeletionScore(int posInQuery, int posInReference)override
		{
			return mutationCost + gapScores->OpenDeletionScore + PureNucleotideInt::LogPrior;
		}

		double ContinueDeletionScore(int posInQuery, int posInReference)override
		{
			return gapScores->ContinueDeletionScore + PureNucleotideInt::LogPrior;
		}

		void SetQuery(shared_ptr<Polymer<vector<double>>> query) override
		{
			this->query = query;
		}

		shared_ptr<Polymer<vector<double>>> GetQuery() override
		{
			return query;
		}

		void SetReference(shared_ptr<Polymer<vector<double>>> reference) override
		{
			this->reference = reference;
		}

		shared_ptr<Polymer<vector<double>>> GetReference() override
		{
			return reference;
		}

		double PriorQueryScore(int posInQuery)override
		{
			return PureNucleotideInt::LogPrior;
		}

		double PriorQueryScore() override
		{
			return PureNucleotideInt::LogPrior;
		}

		vector<double> PriorQueryScoreArray() override
		{
			vector<double> returnValue(query->Seq.size(), 0);
			returnValue[0] = PureNucleotideInt::LogPrior;
			for (unsigned int i = 1; i < query->Seq.size(); i++)
			{
				returnValue[i] = PureNucleotideInt::LogPrior + returnValue[i - 1];
			}
			return returnValue;
		}
	};

	class LikelihoodByteAlignmentKernel : public AlignmentKernelAsymmetric < vector<double>, byte >
	{
	private:
		shared_ptr<Polymer<vector<double>>> query;
		shared_ptr<Polymer<byte>> reference;
		shared_ptr<GapScores> gapScores;
		shared_ptr<PolymerPair<vector<double>, byte>> pair;

		vector<vector<double>> loglikelihood;
		double distance;

		double set_Distance(double value)
		{
			distance = value;
			mutationCost = log(1 - exp(-distance));
		}

		double mutationCost;

	public:
		void set_GapScores(shared_ptr<GapScores> gs)
		{
			gapScores = gs;
		}

		LikelihoodByteAlignmentKernel(double distance)
			: distance(distance), mutationCost(log(1 - exp(-distance))), gapScores(new GapScores()), pair(new PolymerPair<vector<double>, byte>())
		{
		}

		LikelihoodByteAlignmentKernel(shared_ptr<Polymer<vector<double>>> query, shared_ptr<Polymer<byte>> reference, double distance)
			:query(query), reference(reference), distance(distance), mutationCost(log(1 - exp(-distance))),
			gapScores(new GapScores())
		{
			makeLogLikelihood();
		}

		LikelihoodByteAlignmentKernel(shared_ptr<PolymerPair<vector<double>, byte>> pair, double distance)
			:query(pair->Polymer0), reference(pair->Polymer1), distance(distance), mutationCost(log(1 - exp(-distance))), gapScores(new GapScores())
		{
			makeLogLikelihood();
		}
	private:

		void makeLogLikelihood()
		{
			loglikelihood.resize(query->Seq.size(), vector<double>(17, 0));
			for (unsigned int pos = 0; pos < query->Seq.size(); pos++)
			{
				for (byte b = 0; b < 16; b++)
				{
					auto pb = TypeConverter::ToPMF(b);
					double averageloglikelihood = 0;
					for (int iNuc = 0; iNuc < 5; iNuc++)
					{
						if ((query->Seq)[pos][iNuc] > 0)
						{
							averageloglikelihood += pb[iNuc] * (query->Seq)[pos][iNuc];
						}
					}
					loglikelihood[pos][b] = averageloglikelihood;
				}
				loglikelihood[pos][16] = log((query->Seq)[pos][4]);
			}
		}
	public:

		bool IsGap(byte& n) override
		{
			return n == 16;
		}

		double ComparisonScore(int posInQuery, int posInReference) override
		{
			return loglikelihood[posInQuery][(reference->Seq)[posInReference]];
		}

		double OpenInsertionScore(int posInQuery, int posInReference)override
		{
			return loglikelihood[posInQuery][16];
			// compare: mutationCost + openInsertionScore + logPrior
		}

		double ContinueInsertionScore(int posInQuery, int posInReference) override
		{
			return loglikelihood[posInQuery][16] + gapScores->ContinueInsertionScore - gapScores->OpenInsertionScore;
			// compare: continueInsertionScore + logPrior
		}

		double OpenDeletionScore(int posInQuery, int posInReference) override
		{
			return mutationCost + gapScores->OpenDeletionScore;
			// compare: same. But make sure mutationCost is equivalent
		}

		double ContinueDeletionScore(int posInQuery, int posInReference) override
		{
			return gapScores->ContinueDeletionScore;
			// compare: same
		}

		void SetQuery(shared_ptr<Polymer<vector<double>>> query)override
		{
			this->query = query;
			makeLogLikelihood();
		}

		shared_ptr<Polymer<vector<double>>> GetQuery() override
		{
			return query;
		}

		void SetReference(shared_ptr<Polymer<byte>> reference) override
		{
			this->reference = reference;
		}

		shared_ptr<Polymer<byte>> GetReference() override
		{
			return reference;
		}

		double PriorQueryScore(int posInQuery) override
		{
			return loglikelihood[posInQuery][15] + MixedNucleotideByte::LogPrior;
		}

		double PriorQueryScore() override
		{
			return std::numeric_limits<double>::quiet_NaN();
		}

		vector<double> PriorQueryScoreArray() override
		{
			vector<double> returnValue(query->Seq.size(), 0);
			returnValue[0] = PriorQueryScore(0);
			for (unsigned int i = 1; i < query->Seq.size(); i++)
			{
				returnValue[i] = PriorQueryScore(i) + returnValue[i - 1];
			}
			return returnValue;
		}

		static vector<double> NScoreArray(Polymer<vector<double>>& queryLikelihood)
		{
			vector<double> loglikelihood(queryLikelihood.Seq.size(), 0);
			vector<double> pmfForN = { 0.25, 0.25, 0.25, 0.25, 0 };
			for (unsigned int pos = 0; pos < queryLikelihood.Seq.size(); pos++)
			{
				double likelihood = 0;
				likelihood = 0;
				for (int iNuc = 0; iNuc < 4; iNuc++)
				{
					likelihood += (queryLikelihood.Seq)[pos][iNuc];
				}
				loglikelihood[pos] = log(likelihood) + PureNucleotideInt::LogPrior;
			}

			// make cumulative
			for (unsigned int pos = 1; pos < queryLikelihood.Seq.size(); pos++)
			{
				loglikelihood[pos] += loglikelihood[pos - 1];
			}
			return loglikelihood;
		}
	};




	/// <summary>
	/// Abstract class that specifies what methods alignment kernel must provide. It may compare Polymers of different nucleotide encoding types.
	/// </summary>
	/// <typeparam name="T">The type of the first nucleotide encoding.</typeparam>
	/// <typeparam name="U">The type of the second nucleotide encoding.</typeparam>
	template<class T>
	class AlignmentKernelSymmetric
	{
	public:
		/// <summary>
		/// Reports the status of the symbol as a gap.
		/// </summary>
		/// <param name="n">The symbol to be tested.</param>
		/// <returns>True if the symbol represents a gap.</returns>
		virtual bool IsGap(T& n) = 0;
		/// <summary>
		/// Computes the comparison score for the alignment of the nucleotide at a given position in the query with the 
		/// nucleotide at a given position in the reference.
		/// </summary>
		/// <param name="pos0">The position in query.</param>
		/// <param name="posInReference">The position in the reference.</param>
		/// <returns></returns>
		virtual double ComparisonScore(int pos0, int pos1) = 0;
		/// <summary>
		/// Computes the score for continuing an insertion in the query at a given query position when that position 
		/// is aligned against a given position in the reference.
		/// </summary>
		/// <param name="posInQuery">The query position.</param>
		/// <param name="posInReference">The reference position.</param>
		/// <returns>The cost of opening an insertion at the given coordinates.</returns>
		virtual double OpenGapScore(int pos0, int pos1) = 0;


		virtual double OpenInsertionScore(int pos0, int pos1) = 0;

		virtual double OpenDeletionScore(int pos0, int pos1) = 0;



		/// <summary>
		/// Sets the query Polymer.
		/// </summary>
		/// <param name="query">The Polymer to set as query.</param>
		virtual double ContinueGapScore(int pos0, int pos1) = 0;
		/// <summary>
		/// Computes the score for opening a deletion in the query at a given query position when that position 
		/// is aligned against a given position in the reference.
		/// </summary>
		/// <param name="posInQuery">The query position.</param>
		/// <param name="posInReference">The reference position.</param>
		/// <returns>The cost of opening an insertion at the given coordinates.</returns>
		virtual void SetSeq(shared_ptr<Polymer<T>>& seq, int sequenceNumber) = 0;
		/// <summary>
		/// Gets the query.
		/// </summary>
		/// <returns>The query Polymer.</returns>
		virtual shared_ptr<Polymer<T>>& GetSeq(int sequenceNumber) = 0;

		virtual double PriorScore(int sequenceNumber, int position) = 0;

		virtual vector<double> GetPriorScoreArray(int sequenceNumber) = 0;
	};

	class PMFAlignmentKernelSymmetric : public AlignmentKernelSymmetric < vector<double> >
	{
	private:
		shared_ptr<Polymer<vector<double>>> sequences[2];

		double mutationCost = MixedNucleotideByte::LogPrior;
		vector<double> indelBaseProbability[2];
		double OpenGapScoreVal;
		double ContinueGapScoreVal;

	public:

		shared_ptr<GapScores> gapScores;

		PMFAlignmentKernelSymmetric() : gapScores(new GapScores())
		{
			//for symetric Kernel using same gapProbabilities.
			GapProbabilities gp;
			gp.ContinueDeletionProbability = gp.ContinueInsertionProbability = gp.OpenDeletionProbability = gp.OpenInsertionProbability = 0.01;
			gapScores = std::make_shared<GapScores>(gp);
			OpenGapScoreVal = mutationCost + gapScores->OpenInsertionScore + PureNucleotideInt::LogPrior;
			ContinueGapScoreVal = gapScores->ContinueInsertionScore + PureNucleotideInt::LogPrior;
		}

		PMFAlignmentKernelSymmetric(shared_ptr<Polymer<vector<double>>> sequence0, shared_ptr<Polymer<vector<double>>> sequence1)
			: gapScores(new GapScores())
		{
			GapProbabilities gp;
			gp.ContinueDeletionProbability = gp.ContinueInsertionProbability = gp.OpenDeletionProbability = gp.OpenInsertionProbability = 0.01;
			gapScores = std::make_shared<GapScores>(gp);

			sequences[0] = sequence0;
			sequences[1] = sequence1;
			computeIndelBaseProbabilities(0);
			computeIndelBaseProbabilities(1);
			OpenGapScoreVal = mutationCost + gapScores->OpenInsertionScore + PureNucleotideInt::LogPrior;
			ContinueGapScoreVal = gapScores->ContinueInsertionScore + PureNucleotideInt::LogPrior;
		}

		PMFAlignmentKernelSymmetric(PolymerPair<vector<double>, vector<double>>& pair)
			: gapScores(new GapScores())
		{
			GapProbabilities gp;
			gp.ContinueDeletionProbability = gp.ContinueInsertionProbability = gp.OpenDeletionProbability = gp.OpenInsertionProbability = 0.01;
			gapScores = std::make_shared<GapScores>(gp);
			sequences[0] = pair.Polymer0;
			sequences[1] = pair.Polymer1;
			computeIndelBaseProbabilities(0);
			computeIndelBaseProbabilities(1);
			OpenGapScoreVal = mutationCost + gapScores->OpenInsertionScore + PureNucleotideInt::LogPrior;
			ContinueGapScoreVal = gapScores->ContinueInsertionScore + PureNucleotideInt::LogPrior;
		}

	private:
		void computeIndelBaseProbabilities(int sequenceNumber)
		{
			indelBaseProbability[sequenceNumber].resize(sequences[sequenceNumber]->Seq.size());
			for (int i = 0; i < sequences[sequenceNumber]->Seq.size(); i++)
			{
				indelBaseProbability[sequenceNumber][i] = 1 - NucleotidePMF::IDFlux((sequences[sequenceNumber]->Seq)[i], (sequences[sequenceNumber]->Seq)[i]);
			}
		}

	public:
		bool IsGap(vector<double>& n) override
		{
			return n[4] == 1;
		}

		double ComparisonScore(int position0, int position1) override
		{
			auto& seq1 = sequences[0]->Seq;
			auto& seq2 = sequences[1]->Seq;
			return log(NucleotidePMF::IDFlux(seq1[position0], seq2[position1])) + PureNucleotideInt::LogPrior;
		}

		double OpenGapScore(int position0, int position1)override
		{
			return OpenGapScoreVal;
		}


		/// <summary>
		/// insertion into sequence 0
		/// </summary>
		/// <param name="position0"></param>
		/// <param name="position1"></param>
		/// <returns></returns>
		double OpenInsertionScore(int position0, int position1) override
		{
			return OpenGapScoreVal;
		}

		double OpenDeletionScore(int position0, int position1) override
		{
			return OpenGapScoreVal;
		}



		double OpenGapScoreBounded(int position0, int position1)
		{
			return log(indelBaseProbability[0][position0] + indelBaseProbability[1][position1]) + gapScores->OpenInsertionScore + PureNucleotideInt::LogPrior;
		}


		/// <summary>
		/// insertion into sequence 0
		/// </summary>
		/// <param name="position0"></param>
		/// <param name="position1"></param>
		/// <returns></returns>
		double OpenInsertionScoreBounded(int position0, int position1)
		{
			if (position0 + 1 < indelBaseProbability[0].size())
			{
				return log((indelBaseProbability[0][position0] + indelBaseProbability[0][position0 + 1]) / 2 + indelBaseProbability[1][position1]) +
					gapScores->OpenInsertionScore + PureNucleotideInt::LogPrior;
			}
			else return OpenGapScoreBounded(position0, position1);
		}

		double OpenDeletionScoreBounded(int position0, int position1)
		{
			if (position1 + 1 < indelBaseProbability[1].size())
			{
				double a = indelBaseProbability[0][position0];
				double b = indelBaseProbability[1][position1];
				double c = indelBaseProbability[1][position1 + 1];
				double d = a + (b + c) / 2;
				return log(indelBaseProbability[0][position0] + (indelBaseProbability[1][position1] + indelBaseProbability[1][position1 + 1]) / 2) +
					gapScores->OpenDeletionScore + PureNucleotideInt::LogPrior;
			}
			else return OpenGapScoreBounded(position0, position1);
		}



		double ContinueGapScore(int posInQuery, int posInReference)override
		{
			return ContinueGapScoreVal;
			/*return gapScores->ContinueInsertionScore + PureNucleotideInt::LogPrior;*/
		}

		void SetSeq(shared_ptr<Polymer<vector<double>>>& sequence, int sequenceNumber)override
		{
			sequences[sequenceNumber] = sequence;
			computeIndelBaseProbabilities(sequenceNumber);
		}

		shared_ptr<Polymer<vector<double>>>& GetSeq(int sequenceNumber) override
		{

			if (sequenceNumber >= 2)
			{
				throw "wrong argument value exception in GetSeq";
				//return std::shared_ptr<Polymer<vector<double>>>();
			}
			return sequences[sequenceNumber];
		}

		double PriorScore(int sequenceNumber, int position)override
		{
			return PureNucleotideInt::LogPrior;
		}

		vector<double> GetPriorScoreArray(int sequenceNumber)override
		{
			vector<double> returnValue(sequences[sequenceNumber]->Seq.size(), 0);
			returnValue[0] = PureNucleotideInt::LogPrior;
			for (unsigned int i = 1; i < sequences[sequenceNumber]->Seq.size(); i++)
			{
				returnValue[i] = PureNucleotideInt::LogPrior + returnValue[i - 1];
			}
			return returnValue;
		}
	};

	class MNMNAlignmentKernelSymmetric : public AlignmentKernelSymmetric < byte >
	{
	private:
		unique_ptr<NucleotideEvolutionModel> model;

		vector<shared_ptr<Polymer<byte>>> sequences;

		unique_ptr<GapScores> gapScores;

		double mutationCost;

	public:

		double distance;

		void set_Distance(double value)
		{
			distance = value;
			model->set_Time(distance);
			computeComparisonScoreArray();
			mutationCost = log(1 - exp(-distance));
		}

		vector<vector<double>> comparisonScore;

		MNMNAlignmentKernelSymmetric(unique_ptr<NucleotideEvolutionModel> model)
			:model(std::move(model)), gapScores(new GapScores())
		{
			distance = model->get_Time();
			mutationCost = log(1 - exp(-distance));
			computeComparisonScoreArray();
		}

		MNMNAlignmentKernelSymmetric(unique_ptr<NucleotideEvolutionModel> model, shared_ptr<Polymer<byte>> sequence0, shared_ptr<Polymer<byte>> sequence1)
			: model(std::move(model)), gapScores(new GapScores())
		{
			sequences.push_back(sequence0);
			sequences.push_back(sequence1);
			distance = model->get_Time();
			mutationCost = log(1 - exp(-distance));
			computeComparisonScoreArray();
		}

		MNMNAlignmentKernelSymmetric(unique_ptr<NucleotideEvolutionModel> model, PolymerPair<byte, byte> &pair)
			:model(std::move(model)), gapScores(new GapScores())
		{
			sequences.push_back(pair.Polymer0);
			sequences.push_back(pair.Polymer1);
			this->distance = model->get_Time();
		}

		void computeComparisonScoreArray()
		{
			comparisonScore.resize(MixedNucleotideByte::NSymbols, vector<double>(MixedNucleotideByte::NSymbols, 0));
			for (byte b1 = 0; b1 < MixedNucleotideByte::NSymbols; b1++)
			{
				vector<double> p1 = TypeConverter::ToPMF(b1);
				for (byte b2 = 0; b2 < MixedNucleotideByte::NSymbols; b2++)
				{
					vector<double> p2 = TypeConverter::ToPMF(b2);
					for (int i1 = 0; i1 < 4; i1++)
					{
						double x = 0;
						for (int i2 = 0; i2 < 4; i2++)
						{
							x += p2[i2] * model->TransitionProbability(i1, i2);
						}
						comparisonScore[b1][b2] += p1[i1] * log(x);
					}
				}
			}
		}

		bool IsGap(byte& n)override
		{
			return n == 16;
		}

		double ComparisonScore(int position0, int position1)override
		{
			auto i = (sequences[0]->Seq)[position0];
			auto j = (sequences[1]->Seq)[position1];
			return comparisonScore[i][j];

			return comparisonScore[i][j];
		}

		double ContinueGapScore(int position0, int position1)override
		{
			return gapScores->ContinueInsertionScore + PureNucleotideInt::LogPrior;
		}

		double OpenGapScore(int position0, int position1)override
		{
			return mutationCost + gapScores->OpenDeletionScore;
		}

		double OpenInsertionScore(int pos0, int pos1)
		{
			throw std::runtime_error("NotImplementedException");
		}

		double OpenDeletionScore(int pos0, int pos1) override
		{
			throw std::runtime_error("NotImplementedException");
		}

		shared_ptr<Polymer<byte>>& GetSeq(int sequenceNumber)override
		{
			return sequences[sequenceNumber];
		}
		void SetSeq(shared_ptr<Polymer<byte>>& sequence, int sequenceNumber)override
		{
			sequences[sequenceNumber] = sequence;
		}

		double PriorScore(int position0, int position1)override
		{
			return PureNucleotideInt::LogPrior;
		}

		vector<double> GetPriorScoreArray(int sequenceNumber)override
		{
			vector<double> returnValue(sequences[sequenceNumber]->Seq.size());
			returnValue[0] = PureNucleotideInt::LogPrior;
			for (unsigned int i = 1; i < sequences[sequenceNumber]->Seq.size(); i++)
			{
				returnValue[i] = PureNucleotideInt::LogPrior + returnValue[i - 1];
			}
			return returnValue;
		}
	};



}



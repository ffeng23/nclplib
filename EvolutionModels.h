#pragma once

#include <stdlib.h>
#include <math.h>
#include <vector>
#include "typedef.h"
#include "Monomer.h"
#include "Polymer.h"

namespace CLPLib
{
	/// <summary>
	/// Abstract class for Markov models of nucleotide evolution.
	/// </summary>
	class NucleotideEvolutionModel
	{
	public:

		NucleotideEvolutionModel(double time, const std::vector<double>& parameter)
		{
			this->time = time;
			this->Parameters = parameter;
		}

		/// <summary>
		/// The probability that a nucleotide inititially in state originNucleotide has state destinationNucleotide at the specified time.
		/// </summary>
		/// <param name="originNucleotide">The original nucleotide state.</param>
		/// <param name="destinationNucleotide">The state of the nucleotide at the observation time.</param>
		/// <returns>The probability that originNucleotide is destinationNucleotide at the (default) observation time.</returns>
		virtual double TransitionProbability(int originNucleotide, int destinationNucleotide) = 0;
		/// <summary>
		/// The probability that a nucleotide inititially in state originNucleotide has state destinationNucleotide at the specified time.
		/// </summary>
		/// <param name="originNucleotide">The original nucleotide state.</param>
		/// <param name="destinationNucleotide">The state of the nucleotide at the observation time.</param>
		/// <param name="_time">The branch length over which the observation occurs.</param>
		/// <returns>The probability that originNucleotide is destinationNucleotide at the observation time.</returns> 
		virtual double TransitionProbability(int originNucleotide, int destinationNucleotide, double _time) = 0;

		virtual double dTransitionProbabilityDT(int originNucleotide, int destinationNucleotide, double _time) = 0;

		virtual std::vector<double> TypeProbabilities(double _distance) = 0;

		virtual std::vector<double> DestinationPMF(std::vector<double>& originPMF, int len, double _time) = 0;

		virtual std::vector<double> DestinationPMF(std::vector<double>& originPMF, int len) = 0;

		virtual std::vector<double> DestinationPMF(std::vector<double>& originPMF) = 0;

		double get_Time()
		{
			return time;
		}

		virtual void set_Time(double value)
		{
			time = value;
		}

	protected:
		double time;
	public:
		std::vector<double> Parameters;
	};

	/// <summary>
	/// Implements the Jukes-Cantor evolution model.
	/// </summary>
	class JCEvolutionModel : public NucleotideEvolutionModel
	{
	private:
		double onDiagonal;
		double offDiagonal;
		double p;

	public:
		virtual void set_Time(double value) override
		{
			time = value;
			p = 0.75 * (1 - exp(-time));
			onDiagonal = 1 - p;
			offDiagonal = p / 3;
		}

		JCEvolutionModel(double _time) : NucleotideEvolutionModel(_time, std::vector < double > {1})
		{
			time = _time;
			p = 0.75 * (1 - exp(-_time));
			onDiagonal = 1 - p;
			offDiagonal = p / 3;
		}

		double TransitionProbability(int originNucleotide, int destinationNucleotide) override
		{
			// return zero if either is a gap
			if (originNucleotide > 3 || destinationNucleotide > 3)
				return 0;

			return originNucleotide == destinationNucleotide ? onDiagonal : offDiagonal;
		}

		double TransitionProbability(int originNucleotide, int destinationNucleotide, double _time) override
		{
			// return zero if either is a gap
			if (originNucleotide > 3 || destinationNucleotide > 3)
				return 0;

			double _p = 0.75 * (1 - exp(-_time));
			return originNucleotide == destinationNucleotide ? 1 - _p : _p / 3;
		}

		double dTransitionProbabilityDT(int originNucleotide, int destinationNucleotide, double _time) override
		{
			// return zero if either is a gap
			if (originNucleotide > 3 || destinationNucleotide > 3)
				return 0;

			if (originNucleotide == destinationNucleotide)
			{
				return -0.75 * exp(-_time);
			}
			else
			{
				return 0.25 * exp(-_time);
			}
		}

		static double ProbabilityToTime(double probability)
		{
			return -log(1 - 1.333333 * probability);
		}

		std::vector<double> TypeProbabilities(double _time) override
		{
			double _p = 0.75 * (1 - exp(-_time));
			return std::vector<double>{1 - _p, _p, 2 * _p };
		}


		/// <summary>
		/// Computes the PMF resulting from Juke-Cantor evolution over the default time acting on a given initial state.
		/// </summary>
		/// <param name="originNucleotide">The initial PMF.</param>
		/// <returns>The PMF resulting from evolution over time from the initial PMF.</returns>
		std::vector<double> DestinationPMF(std::vector<double>& originPMF) override
		{
			// TODO: test to see if originNucleotide is a gap

			std::vector<double> destinationPMF(originPMF.size());
			if (originPMF.size() > 4)
			{
				destinationPMF[4] = originPMF[4];
			}
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					if (i == j) // on-diagonal
					{
						destinationPMF[i] += onDiagonal * originPMF[j];
					}
					else
					{
						destinationPMF[i] += offDiagonal * originPMF[j];
					}
				}
			}
			return destinationPMF;
		}


		/// <summary>
		/// Computes the PMF resulting from Juke-Cantor evolution over the default time acting on a given initial state.
		/// </summary>
		/// <param name="originNucleotide">The initial PMF.</param>
		/// <returns>The PMF resulting from evolution over time from the initial PMF.</returns>
		std::vector<double> DestinationPMF(std::vector<double>& originPMF, int len) override
		{
			// TODO: test to see if originNucleotide is a gap
			std::vector<double> destinationPMF(len, 0);
			if (len > 4)
			{
				destinationPMF[4] = originPMF[4];
			}
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					if (i == j) // on-diagonal
					{
						destinationPMF[i] += onDiagonal * originPMF[j];
					}
					else
					{
						destinationPMF[i] += offDiagonal * originPMF[j];
					}
				}
			}
			return destinationPMF;
		}

		/// <summary>
		/// Computes the PMF resulting from Jukes-Cantor evolution over time acting on a given initial state.
		/// </summary>
		/// <param name="originNucleotide">The initial PMF.</param>
		/// <param name="_time">The evolutionary time.</param>
		/// <returns>The PMF resulting from evolution over time _time from the initial PMF.</returns>        
		std::vector<double> DestinationPMF(std::vector<double>& originPMF, int len, double _time) override
		{
			// TODO: test to see if originNucleotide is a gap

			double _p = 0.75 * (1 - exp(-_time));
			double _onDiagonal = 1 - _p;
			double _offDiagonal = _p / 3;

			std::vector<double> destinationPMF(len, 0);
			if (len > 4)
			{
				destinationPMF[4] = originPMF[4];
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					if (i == j)
					{
						destinationPMF[i] += _onDiagonal * originPMF[j];
					}
					else
					{
						destinationPMF[i] += _offDiagonal * originPMF[j];
					}
				}
			}
			return destinationPMF;
		}

	};

	/// <summary>
	/// Implements the Kimura evolution model.
	/// </summary>
	class KimuraEvolutionModel : public NucleotideEvolutionModel
	{
	private:
		double onDiagonal;
		double offDiagonalTransition;
		double offDiagonalTransversion;
		double tau;
	public:
		double Kappa;
		void set_Time(double value) override
		{

			time = value;
			tau = 2 * time / (Kappa + 1);
			offDiagonalTransversion = 0.25 * (1 - exp(-tau));
			offDiagonalTransition = 0.25 * (1 + exp(-tau) - 2 * exp(-time));
			onDiagonal = 1 - offDiagonalTransition - 2 * offDiagonalTransversion;

		}

		KimuraEvolutionModel(double _time, std::vector<double>& kappa) : NucleotideEvolutionModel(_time, kappa)
		{
			Kappa = kappa[0];
			time = _time;
			tau = 2 * time / (Kappa + 1);
			offDiagonalTransversion = 0.25 * (1 - exp(-tau));
			offDiagonalTransition = 0.25 * (1 + exp(-tau) - 2 * exp(-time));
			onDiagonal = 1 - offDiagonalTransition - 2 * offDiagonalTransversion;
		}



		std::vector<double> DestinationPMF(std::vector<double>& originPMF) override
		{
			throw "not implemented exception";
		}


		/// <summary>
		/// Computes the PMF resulting from Kimura evolution over the defaultBranchLength acting on a given initial pure state.
		/// </summary>
		/// <param name="originNucleotide">The initial PMF.</param>
		/// <param name="mu">The evolutionary distance mu.</param>
		/// <returns>The PMF resulting from evolution over distance mu from the initial PMF.</returns>
		std::vector<double> DestinationPMF(std::vector<double>& originPMF, int len) override
		{
			// TODO: test to see if originNucleotide is a gap

			std::vector<double> destinationPMF(len, 0);
			if (len > 4)
			{
				destinationPMF[4] = originPMF[4];
			}
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					if (i == j) // on-diagonal
					{
						destinationPMF[i] += onDiagonal * originPMF[j];
					}
					else if (i / 2 == j / 2) // transition
					{
						destinationPMF[i] += offDiagonalTransition * originPMF[j];
					}
					else
					{
						destinationPMF[i] += offDiagonalTransversion * originPMF[j];
					}
				}
			}
			return destinationPMF;
		}

		/// <summary>
		/// Computes the PMF resulting from Kimura evolution over distance branchLength acting on a given initial pure state.
		/// </summary>
		/// <param name="originNucleotide">The initial PMF.</param>
		/// <param name="mu">The evolutionary distance mu.</param>
		/// <returns>The PMF resulting from evolution over distance mu from the initial PMF.</returns>        
		std::vector<double> DestinationPMF(std::vector<double>& originPMF, int len, double _time) override
		{
			// TODO: test to see if originNucleotide is a gap

			double _tau = 2 * _time / (Kappa + 1);
			double x = exp(-_tau);
			double _offDiagonalTransversion = (1 - x) / 4;
			double _offDiagonalTransition = (1 + x - 2 * exp(-_time)) / 4;
			double _onDiagonal = 1 - _offDiagonalTransition - 2 * _offDiagonalTransversion;

			std::vector<double> destinationPMF(len, 0);
			if (len > 4)
			{
				destinationPMF[4] = originPMF[4];
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					if (i == j)
					{
						destinationPMF[i] += _onDiagonal * originPMF[j];
					}
					else if (i / 2 == j / 2)
					{
						destinationPMF[i] += _offDiagonalTransition * originPMF[j];
					}
					else
					{
						destinationPMF[i] += _offDiagonalTransversion * originPMF[j];
					}
				}
			}
			return destinationPMF;
		}

		double TransitionProbability(int originNucleotide, int destinationNucleotide) override
		{
			// return zero if either is a gap
			if (originNucleotide > 3 || destinationNucleotide > 3)
				return 0;

			if (originNucleotide == destinationNucleotide)
			{
				return onDiagonal;
			}
			else if (originNucleotide / 2 == destinationNucleotide / 2)
			{
				return offDiagonalTransition;
			}
			else
			{
				return offDiagonalTransversion;
			}
		}

		double TransitionProbability(int originNucleotide, int destinationNucleotide, double _time) override
		{

			// return zero if either is a gap
			if (originNucleotide > 3 || destinationNucleotide > 3)
				return 0;

			double x = exp(-2 * _time / (Kappa + 1));
			if (originNucleotide == destinationNucleotide) // identity
			{
				return (1 + x + 2 * exp(-_time)) / 4;
			}
			else if (originNucleotide / 2 == destinationNucleotide / 2) // transition
			{
				return (1 + x - 2 * exp(-_time)) / 4;
			}
			else // transversion
			{
				return (1 - x) / 4;
			}
		}

		//double TransitionProbability(byte originNucleotide, byte destinationNucleotide)
		//{
		//	// return zero if either is a gap
		//	if (originNucleotide > 3 || destinationNucleotide > 3)
		//		return 0;

		//	if (originNucleotide == destinationNucleotide)
		//	{
		//		return onDiagonal;
		//	}
		//	else if (originNucleotide / 2 == destinationNucleotide / 2)
		//	{
		//		return offDiagonalTransition;
		//	}
		//	else
		//	{
		//		return offDiagonalTransversion;
		//	}
		//}

		//double TransitionProbability(byte originNucleotide, byte destinationNucleotide, double _time)
		//{
		//	// return zero if either is a gap
		//	if (originNucleotide > 3 || destinationNucleotide > 3)
		//		return 0;

		//	if (originNucleotide == destinationNucleotide)
		//	{
		//		double _tau = 2 * _time / (Kappa + 1);
		//		double x = exp(-_tau);
		//		double _TwoTimesOffDiagonalTransversion = 0.5 * (1 - x);
		//		double _offDiagonalTransition = 0.25 * (1 + x - 2 * exp(-_time));
		//		double _onDiagonal = 1 - _offDiagonalTransition - _TwoTimesOffDiagonalTransversion;
		//		return _onDiagonal;
		//	}
		//	else if (originNucleotide / 2 == destinationNucleotide / 2)
		//	{
		//		double _tau = 2 * _time / (Kappa + 1);
		//		double _offDiagonalTransition = 0.25 * (1 + exp(-_tau) - 2 * exp(-_time));
		//		return _offDiagonalTransition;
		//	}
		//	else
		//	{
		//		double _tau = 2 * _time / (Kappa + 1);
		//		double _offDiagonalTransversion = 0.25 * (1 - exp(-_tau));
		//		return _offDiagonalTransversion;
		//	}
		//}

		double dTransitionProbabilityDT(int originNucleotide, int destinationNucleotide, double _time) override
		{
			// return zero if either is a gap
			if (originNucleotide > 3 || destinationNucleotide > 3)
				return 0;

			double x = exp(-2 * _time / (Kappa + 1));
			double y = exp(-_time);
			if (originNucleotide == destinationNucleotide) // identity
			{
				return -(y + x / (Kappa + 1)) / 2;
			}
			else if (originNucleotide / 2 == destinationNucleotide / 2) // transition
			{
				return (y - x / (Kappa + 1)) / 2;
			}
			else // transversion
			{
				return x / (2 * (Kappa + 1));
			}
		}

		std::vector<double> TypeProbabilities(double _time) override
		{
			double _tau = 2 * _time / (Kappa + 1);
			double x = exp(-_tau);
			double _TwoTimesOffDiagonalTransversion = 0.5 * (1 - x);
			double _offDiagonalTransition = 0.25 * (1 + x - 2 * exp(-_time));
			double _onDiagonal = 1 - _offDiagonalTransition - _TwoTimesOffDiagonalTransversion;
			return std::vector<double> { _onDiagonal, _offDiagonalTransition, _TwoTimesOffDiagonalTransversion};
		}
	};

	class Evolver
	{
	public:
		static std::shared_ptr<NucleotideEvolutionModel> Model;
		static std::shared_ptr<NucleotidePMF> encoder;

		static Polymer<std::vector<double>> Evolve(Polymer<std::vector<double>>& p)
		{

			Polymer<std::vector<double>> pe(encoder, p.UID);
			pe.Seq.reserve(p.Seq.size());			
			for (unsigned int i = 0; i < p.Seq.size(); i++)
			{
				pe.Seq[i] = Model->DestinationPMF(p.Seq[i]);
			}
			return pe;
		}

		static std::vector<std::vector<double>> Evolve(std::vector<std::vector<double>>& sequence)
		{
			std::vector<std::vector<double>> evolvedSeq(sequence.size(), std::vector<double>());
			for (unsigned int i = 0; i < sequence.size(); i++)
			{
				evolvedSeq[i] = Model->DestinationPMF(sequence[i]);
			}
			return evolvedSeq;
		}

	};



}


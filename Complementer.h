#pragma once

#include <memory>
#include <vector>

namespace CLPLib
{

	//forward declaration
	template <class T> class MonomerEncoding;
	template <class T> class Polymer;

	template <class T>
	class Complementer
	{
	private:

		std::shared_ptr<MonomerEncoding<T>> nucleotideEncoding;

	public:

		//std::make_shared<MixedNucleotideByte>() - set this as default
		Complementer() : nucleotideEncoding(std::shared_ptr<MonomerEncoding<T>>())
		{}

		/// <summary>
		/// Constructor.
		/// </summary>
		/// <param name="_nucleotide">An instance of the relevant NucleotideEncoding class.</param>
		//Complementer(std::shared_ptr<MonomerEncoding<T>> _nucleotide)
		//	:nucleotideEncoding(_nucleotide)
		//{
		//}

		//this exists so that we can avoid pointers
		//be careful when setting encoding like this.
		//if the Complementer is declared as static, and different threads are setting this
		//it cause cause a lot of problems here.
		void SetEncoding_ex(std::shared_ptr<MonomerEncoding<T>>& sp)
		{
			fprintf(stderr, "Warning: calling SetEncoding can cause access violation if this class is static member of Polymer and multithreading is used.\n");
			throw ("remove this if truely needed but notice warning above.");

			nucleotideEncoding = sp;
		}

		std::shared_ptr<MonomerEncoding<T>> GetEncoding()
		{
			return nucleotideEncoding;
		}

		/// <summary>
		/// Produces the complement of a Polynucleotide.
		/// </summary>
		/// <param name="p">The polynucleotide to be complemented.</param>
		/// <returns>The polynucleotide whose sequence is the complement of the input polynucleotide. The name has ".C" appended.</returns>
		std::shared_ptr<Polymer<T>> Complement(Polymer<T>& p);

		/// <summary>
		/// Produces the complement of a nucleotide array.
		/// </summary>
		/// <param name="seq">An array of nucleotides in the relevant encoding.</param>
		/// <returns>The complement of the input array.</returns>
		std::vector<T> Complement(std::vector<T>& seq);

		///// <summary>
		///// Complements all of the members of Polynucleotide collection.
		///// </summary>
		///// <param name="pc">A polynucleotide collection of the relevant type.</param>
		///// <returns>A polynucleotide collection whose members are all complements of those of the input array.
		///// The name of the return collection is that of the input collection with ".C" appended.</returns>
		//public PolynucleotideCollection<T> Complement(PolynucleotideCollection<T> pc)
		//{
		//    PolynucleotideCollection<T> pcc = new PolynucleotideCollection<T>();
		//    pcc.Name = pc.Name + ".C";
		//    foreach (Polynucleotide<T> p in pc)
		//    {
		//        pcc.AddMember(Complement(p));
		//    }
		//    return pcc;
		//}
	};

}


#pragma once

namespace CLPLib
{

	/// <summary>
	/// Specifies which end of a Polymer is being referred to.
	/// </summary>
	enum NucleicAcidEnd
	{
		/// <summary>
		/// The 5' end.
		/// </summary>
		FivePrimed,
		/// <summary>
		/// The 3' end.
		/// </summary>
		ThreePrimed
	};

	enum SequenceFileFormat { FASTA, FASTQ, Other };

}


#include "stdafx.h"
#include "AlignmentKernel.h"

namespace CLPLib
{


	GapScores::GapScores(GapProbabilities& gapProbabilities)
	{
		OpenDeletionScore = log(gapProbabilities.OpenDeletionProbability);
		OpenInsertionScore = log(gapProbabilities.OpenInsertionProbability);
		ContinueDeletionScore = log(gapProbabilities.ContinueDeletionProbability);
		ContinueInsertionScore = log(gapProbabilities.ContinueInsertionProbability);
	}




}


#include "definition.h"
#include <stdio.h>
#include <list>

#include "AlignmentKernel.h"
#include "PairwiseAligner.h"

using namespace CLPLib;


void print_line(const char * str)
{
	printf(str);
}

int DEBUG_CODE = 0;

std::shared_ptr<MixedNucleotideByte> MNB_encoding;


//global 
PolymerCollection<byte> VGeneList(MNB_encoding);

PolymerCollection<byte> JGeneList(MNB_encoding);

PolymerCollection<byte> DGeneList(MNB_encoding);

std::map <int, std::shared_ptr<Polymer<byte>>> RefSeqMap;

void set_debug_code(int code)
{
	DEBUG_CODE = code;
	fprintf(stderr, "Info: deubg code set to %d\n", code);
}


void *NewAligner_1(void *kernel_ptr, int type)
{
	if (type == 0)
	{
		std::shared_ptr<MNMNAlignmentKernelAsymmetric> _kernel((MNMNAlignmentKernelAsymmetric *)kernel_ptr);
		PairwiseAlignerAsymmetric<byte, byte> *aligner = new PairwiseAlignerAsymmetric<byte, byte>(0, 0, _kernel);

		if (!MNB_encoding)
		{
			MNB_encoding = std::make_shared<MixedNucleotideByte>();
			fprintf(stderr, "MixedNucleotideByte encoding created.\n");
		}
		return aligner;
	}

	std::cerr << "native aligner creation failed.\n";
	return nullptr;
}

void Aligner_1_SetQuery(void *aligner_ptr, unsigned char * barr, int length)
{
	PairwiseAlignerAsymmetric<byte, byte> *aligner = (PairwiseAlignerAsymmetric<byte, byte> *)aligner_ptr;
	std::vector<byte> seq(barr, barr + length);

	//auto encoding = std::make_shared<MixedNucleotideByte>();
	//auto polymer1 = std::make_shared<Polymer<byte>>(encoding, "query", seq);
	aligner->SetQuery(std::make_shared<Polymer<byte>>(MNB_encoding, "query", seq));
}

void Kernel_SetQuery(void *kernel_ptr, unsigned char * barr, int length)
{

	std::vector<byte> seq(barr, barr + length);
	//auto encoding = std::make_shared<MixedNucleotideByte>();
	//auto polymer1 = std::make_shared<Polymer<byte>>(encoding, "query", seq);

	int type = 0; //to do: pass in kernel type.
	if (type == 0)
	{
		auto kernel = (MNMNAlignmentKernelAsymmetric *)kernel_ptr;
		kernel->SetQuery(std::make_shared<Polymer<byte>>(MNB_encoding, "query", seq));
	}
	else
	{
		std::cerr << "kernel type not supporetd.\n";
	}
}


void Aligner_1_SetReference(void *aligner_ptr, unsigned char * barr, int length)
{
	PairwiseAlignerAsymmetric<byte, byte> *aligner = (PairwiseAlignerAsymmetric<byte, byte> *)aligner_ptr;
	std::vector<byte> seq(barr, barr + length);

	//auto encoding = std::make_shared<MixedNucleotideByte>();
	//auto polymer1 = std::make_shared<Polymer<byte>>(encoding, "query", seq);
	aligner->SetReference(std::make_shared<Polymer<byte>>(MNB_encoding, "query", seq));
}

bool Aligner_1_SetReferenceList(int segmentType, int refid, unsigned char * _seq, int length)
{
	auto& genelist = segmentType == 0 ? VGeneList : (segmentType == 1 ? DGeneList : JGeneList);
	std::vector<byte> seq(_seq, _seq + length);
	if (genelist.Contains(refid) == false)
	{
		auto pmptr = std::make_shared<Polymer<byte>>(MNB_encoding, "reference", seq);
		pmptr->Id = refid;
		genelist.AddMember(pmptr);
		RefSeqMap.insert(std::make_pair(refid, pmptr));
		return true;
	}
	return false;
}

void Kernel_SetReference(void *kernel_ptr, unsigned char * barr, int length)
{
	/*
	std::vector<byte> seq(barr, barr + length);

	auto encoding = std::make_shared<MixedNucleotideByte>();
	auto polymer1 = std::make_shared<Polymer<byte>>(encoding, "query", seq);
	*/

	auto polymer1 = std::make_shared<Polymer<byte>>(MNB_encoding, "query", 
		std::vector<byte>(barr, barr+length));

	int type = 0; //to do: pass in kernel type.
	if (type == 0)
	{
		auto kernel = (MNMNAlignmentKernelAsymmetric *)kernel_ptr;
		kernel->SetReference(polymer1);
	}
	else
	{
		std::cerr << "kernel type not supporetd.\n";
	}
}

void Aligner_1_SetFirstColumn(void *aligner_ptr, double * barr, int length)
{
	PairwiseAlignerAsymmetric<byte, byte> *aligner = (PairwiseAlignerAsymmetric<byte, byte> *)aligner_ptr;
	std::vector<double> defaultScores(barr, barr + length);

	aligner->SetFirstColumn(defaultScores);
	//fprintf(stderr, "firstColumn set: %d\n", length);
}


double Aligner_1_FillScoreMatrix(void *aligner_ptr, bool qFivePrimeEndFree, bool rFivePrimeEndFree, bool qThreePrimeEndFree, bool rThreePrimeEndFree)
{
	auto aligner = (PairwiseAlignerAsymmetric<byte, byte> *)aligner_ptr;
	//auto score = aligner->FillScoreMatrix(qFivePrimeEndFree, rFivePrimeEndFree, qThreePrimeEndFree, rThreePrimeEndFree);
	auto score = aligner->FillScoreMatrix();
	return score;
}

double Aligner_1_FillScoreMatrix2(void *aligner_ptr, const char *refname)
{

	throw("not used anymore exception");
	//auto aligner = (PairwiseAlignerAsymmetric<byte, byte> *)aligner_ptr;
	//aligner->SetReference(RefSeqMap[refname]);
	//auto score = aligner->FillScoreMatrix();
	//return score;
}

//segmentType - V(0), D(1), J(2)
//threshold - keep scores if > TopScore - threshold
//refIdList - list of Ids that align.
//scoreList - sumscores of each alignment
//indexLength - lengths of each index array. 
//index1, index2 - concatenated indexs; 
//scores - concatednated scores; - note, we put -2 between entries serving as debugging info.
//indexCapacity - total length of allocated array passed in.
int AlignToReference_0(void *aligner_ptr, int segmentType, double threshold, int *refIdList,
	double* scoreList, int* indexLength, int* index1, int* index2, double* scores, int indexCapacity)
{

	auto aligner = (PairwiseAlignerAsymmetric<byte, byte> *)aligner_ptr;
	double maxScore = std::numeric_limits<double>::lowest();

	int keepIndex = 0;
	int *index1_ptr = index1;
	int *index2_ptr = index2;
	double *score_ptr = scores;
	int capacity = indexCapacity;

	auto& genelist = segmentType == 0 ? VGeneList : (segmentType == 1 ? DGeneList : JGeneList);
	for (auto item : genelist)
	{
		auto& refseq = item.second;
		aligner->SetReference(refseq);
		auto alignScore = aligner->FillScoreMatrix();
		if (alignScore < maxScore - threshold)continue;
		if (alignScore > maxScore)
		{
			maxScore = alignScore;
		}
		refIdList[keepIndex] = refseq->Id;
		scoreList[keepIndex] = alignScore;

		int index_len = aligner->TraceBack(index1_ptr, index2_ptr, score_ptr, capacity);
		indexLength[keepIndex] = index_len;
		keepIndex++;

		index1_ptr += index_len;
		index2_ptr += index_len;
		score_ptr += index_len;
		capacity -= index_len;
		
		//separaters between entries.
		*index1_ptr++ = -2;
		*index2_ptr++ = -2;
		*score_ptr++ = -2;
		capacity -= 1;
	}
	//end of arrays.
	*index1_ptr++ = -1000;
	*index2_ptr++ = -1000;
	*score_ptr++ = -1000;
	return keepIndex;
}


int AlignToReference(void *aligner_ptr, int segmentType, double threshold, int *refIdList,
	double* scoreList, int* indexLength, int* index1, int* index2, double* scores, int indexCapacity)
{

	auto aligner = (PairwiseAlignerAsymmetric<byte, byte> *)aligner_ptr;
	double maxScore = std::numeric_limits<double>::lowest();

	int keepIndex = 0;
	int *index1_ptr = index1;
	int *index2_ptr = index2;
	double *score_ptr = scores;
	int capacity = indexCapacity;

	auto& genelist = segmentType == 0 ? VGeneList : (segmentType == 1 ? DGeneList : JGeneList);

	std::vector<std::pair<shared_ptr<Polymer<byte>>, double>> scorelist;

	//polymer id to maxscore position in scoreMatrix
	std::map<int, std::pair<int, int>> tracebackList;
	for (auto item : genelist)
	{
		auto& refseq = item.second;
		aligner->SetReference(refseq);
		auto alignScore = aligner->FillScoreMatrix_get_score();

		if (alignScore < maxScore - threshold)continue;
		if (alignScore > maxScore)maxScore = alignScore;

		scorelist.push_back(std::make_pair(refseq, alignScore));
		auto pair = std::make_pair(aligner->MaxRowIndex, aligner->MaxColIndex);

		tracebackList.insert(std::make_pair(refseq->Id, pair));
		//std::cerr << "ref_" << refseq->Id << "\tScore=\t" << alignScore << std::endl;
	}

	for (auto& kvp : scorelist)
	{
		if (kvp.second < maxScore - threshold)continue;
		auto& refseq = kvp.first;
		aligner->SetReference(refseq);

		auto maxpair = tracebackList[refseq->Id];
		int max_x = maxpair.first;
		int max_y = maxpair.second;
		//auto alignScore_x = aligner->FillScoreMatrix();
		
		//testing using max_x and max_y
		auto alignScore = aligner->FillScoreMatrix_tom(max_x, max_y);
		if (alignScore != kvp.second || aligner->MaxRowIndex != max_x || aligner->MaxColIndex != max_y)
		{
			std::cerr << "debug error in align to reference\n";
			auto tScore = aligner->FillScoreMatrix_get_score();
			std::cerr << "check point.";
		}

		if (alignScore < maxScore - threshold)continue;
		if (alignScore > maxScore)
		{
			maxScore = alignScore;
		}
		refIdList[keepIndex] = refseq->Id;
		scoreList[keepIndex] = alignScore;

		int index_len = aligner->TraceBack(index1_ptr, index2_ptr, score_ptr, capacity);
		indexLength[keepIndex] = index_len;
		keepIndex++;

		index1_ptr += index_len;
		index2_ptr += index_len;
		score_ptr += index_len;
		capacity -= index_len;

		//separaters between entries.
		*index1_ptr++ = -2;
		*index2_ptr++ = -2;
		*score_ptr++ = -2;
		capacity -= 1;
	}
	//end of arrays.
	*index1_ptr++ = -1000;
	*index2_ptr++ = -1000;
	*score_ptr++ = -1000;

	//std::cerr << "gene segment count = " << genelist.Count() << " alignCount= " << keepIndex << std::endl;

	return keepIndex;
}


int Aligner_1_TraceBack(void *aligner_ptr, int *index1, int *index2, double *scores, int arrayLength)
{
	auto aligner = (PairwiseAlignerAsymmetric<byte, byte> *)aligner_ptr;
	return aligner->TraceBack(index1, index2, scores, arrayLength);
}

//type is type of Kernel
void Aligner_1_SetKernelDistance(void *kernel_ptr, int type, double distance)
{
	if (type == 0)
	{
		auto kernel = (MNMNAlignmentKernelAsymmetric *)kernel_ptr;
		kernel->set_Distance(distance);
		//std::cerr << "kernel distance set at " << distance << "\n";
	}
	else
	{
		std::cerr << "kernel type not supporetd.\n";
	}
}

void Aligner_1_SetKernelGapScore(void *kernel_ptr, int type, double ods, double cds, double ois, double cis)
{
	if (type == 0)
	{
		auto kernel = (MNMNAlignmentKernelAsymmetric *)kernel_ptr;
		kernel->gapScores = std::make_shared<GapScores>(ods, cds, ois, cis);
	}
	else
	{
		std::cerr << "kernel type not supporetd.\n";
	}
}


void* Aligner_1_NewNativeKernel(int model_index, double _time, double _kappa, double ods, double cds, double ois, double cis)
{
	std::shared_ptr<NucleotideEvolutionModel> _model;
	std::vector<double> kv{ _kappa };

	switch (model_index)
	{
	case 0:		//JC 
		_model = std::make_shared<JCEvolutionModel>(_time);
		break;
	case 1:		//Kimura
		_model = std::make_shared<KimuraEvolutionModel>(_time, kv); //std::vector<double>{_kappa});
		break;
	default:
		fprintf(stderr, "unknown model: %d\n", model_index);
		return nullptr;
	}
	MNMNAlignmentKernelAsymmetric *_kernel = new MNMNAlignmentKernelAsymmetric(_model);
	_kernel->gapScores = std::make_shared<GapScores>(ods, cds, ois, cis);
	return _kernel;
}



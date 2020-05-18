#include <iostream>
#include <memory>
#include <vector>
#include "Monomer.h"
#include "PolymerCollection.h"
#include "EvolutionModels.h"
#include "AlignmentKernel.h"
#include "PairwiseAligner.h"
#include "Writer.h"

using namespace CLPLib;

void print_alignment(std::string seq1, std::string seq2, int index_len, int *index1, int *index2)
{
	std::string s1(index_len, ' ');
	std::string s2(index_len, 'X');
	std::string s3(index_len, ' ');

	for (int i = 0; i < index_len; i++)
	{
		s1[i] = index1[i] == -1 ? '-' : seq1[i];
		s3[i] = index2[i] == -1 ? '-' : seq2[i];
		s2[i] = s1[i] == s3[i] ? ' ' : 'X';
	}

	std::cout << "seq1\t" << s1 << std::endl;
	std::cout << "    \t" << s2 << std::endl;
	std::cout << "seq2\t" << s3 << std::endl;
}




int main(int argc, char *argv[])
{

	std::cerr << "starting program...\n";

	std::string seqFile = "C:/Projects/MMSeq/StatAssembler/Release/test2.fastq";

	//PolymerCollection<vector<double>> seqlist(std::make_shared<NucleotidePMF>());
	//seqlist.FromFASTQ(seqFile);

	PolymerCollection<byte> seqlist(std::make_shared<MixedNucleotideByte>());
	
	//seqlist.FromFASTQ(seqFile);
	string seq1 = "CACAGCCTCCCCTCGGTGCCCTGCTACCTCCTCAGGTCAGCCCTGGACATCCCGGGTTTCCCCAGGCCTGGCGGTAGGTTTGGGGTGAGGTCTGTGTCACTGTGGTATTACGATTTTTGGAGTGGTTATTATTCTTGCTACTATTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAGCC";

	string seq2 = "TAGTACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG";


	string seq3 = "ACCTTCAGTGACCACTACATGGACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTGGCCGTACTAGAAACAAAGCTAACAGCTACACCACAGAATACGCCGCGTCTGTGAAAGGCAGATTCACCATCTCAAGAGATGATTCAAAGAACTCACTGTATATAC";
	string seq4 = "GCTCCATCAGCAGTAGTAGTTACTACTGGGGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGAGTATCTATTATAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCCGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACAC";

	
	auto ps1 = std::make_shared<Polymer<byte>>(std::make_shared<MixedNucleotideByte>(), "seq1", TypeConverter::ToBytes(seq1));
	auto ps2 = std::make_shared<Polymer<byte>>(std::make_shared<MixedNucleotideByte>(), "seq2", TypeConverter::ToBytes(seq2));
	auto ps3 = std::make_shared<Polymer<byte>>(std::make_shared<MixedNucleotideByte>(), "seq3", TypeConverter::ToBytes(seq3));

	seqlist.AddMember(ps1);
	seqlist.AddMember(ps2);
	seqlist.AddMember(ps3);

	std::cerr << "sequences read: " << seqlist.Count() << "\n";

	std::vector<double> kappavec{2};
	std::unique_ptr<KimuraEvolutionModel> model(new KimuraEvolutionModel(0.4, kappavec));

	std::unique_ptr<MNMNAlignmentKernelAsymmetric> kernel(new MNMNAlignmentKernelAsymmetric(std::move(model)));

	//default for cloanalyst
	GapProbabilities gp;
	gp.OpenInsertionProbability = gp.OpenDeletionProbability = 0.001;
	gp.ContinueDeletionProbability = gp.ContinueInsertionProbability = 0.9;
	kernel->gapScores = std::make_shared<GapScores>(gp);

	PairwiseAlignerAsymmetric<byte, byte> aligner(std::move(kernel));

	//the seqlist id is 1 based, not 0 based.
	aligner.SetQuery(seqlist[1]);
	double sum_score = 0;

	int *index1 = (int *)malloc(2000 * sizeof(int));
	int *index2 = (int *)malloc(2000 * sizeof(int));
	double *scores = (double *)malloc(2000 * sizeof(double));
	for (int i = 2; i < 100000; i++)
	{
		int j = (i % 2) + 2;

	
		aligner.SetReference(seqlist[j]);

		auto score = aligner.FillScoreMatrix();
		//auto pair = aligner.TraceBack();
		//auto align_str = Writer<unsigned char>::WriteFasta_xx(pair); // WriteFasta(pair);

		int index_len = aligner.TraceBack(index1, index2, scores, 2000);

		if (j == 2)
		{
			print_alignment(seq1, seq2, index_len, index1, index2);
		}

		sum_score += score;
	}

	std::cout << sum_score;
	//auto& seq1 = seqlist[1]->charSeq;
	//auto& seq2 = seqlist[2]->charSeq;

	//std::cerr << seq1 << "\n" << seq2 << "\n";

	//std::cerr << "alignment: \n" << align_str << "\n";

	//getchar();

	return 0;

}

#pragma once

#if defined(_WIN32)
    //  Microsoft
    #define DllExport __declspec(dllexport)
    #define IMPORT __declspec(dllimport)
#elif defined(__GNUC__)
    //  GCC
    #define DllExport __attribute__((visibility("default")))
    #define IMPORT
#else
    //  do nothing and hope for the best?
    #define EXPORT
    #define IMPORT
    #pragma warning Unknown dynamic link import/export semantics.
#endif


#ifdef __cplusplus
extern "C" {
#endif

	DllExport void print_line(const char *str);

	DllExport void set_debug_code(int code);

	DllExport void* NewAligner_1(void *kernel, int type);

	DllExport void Aligner_1_SetQuery(void *inst, unsigned char * seq, int length);

	DllExport void Aligner_1_SetReference(void *inst, unsigned char * seq, int length);

	DllExport bool Aligner_1_SetReferenceList(int segmentType, int refid, unsigned char * seq, int length);


	DllExport void Kernel_SetQuery(void *inst, unsigned char * seq, int length);

	DllExport void Kernel_SetReference(void *inst, unsigned char * seq, int length);

	DllExport void Aligner_1_SetFirstColumn(void *inst, double* seq, int length);

	//return alignment score, as will as tracebackStart
	DllExport double Aligner_1_FillScoreMatrix(void *inst, bool qFivePrimeEndFree, bool rFivePrimeEndFree, bool qThreePrimeEndFree, bool rThreePrimeEndFree);

	DllExport double Aligner_1_FillScoreMatrix2(void *inst, const char *refname);


	DllExport int Aligner_1_TraceBack(void *inst, int *index1, int *index2, double *scores, int arrayLength);

	DllExport void * Aligner_1_NewNativeKernel(int model_index, double _time, double _kappa, double ods, double cds, double ois, double cis);

	DllExport void Aligner_1_SetKernelDistance(void *inst, int type, double distance);

	DllExport void Aligner_1_SetKernelGapScore(void *inst, int type, double ods, double cds, double ois, double cis);


	DllExport int AlignToReference(void *inst, int segtype, double threshold, int *refIdList,
		double* scoreList, int* indexLength, int* index1, int* index2, double* scores, int sumlen);


#ifdef __cplusplus
}
#endif

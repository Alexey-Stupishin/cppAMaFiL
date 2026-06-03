#pragma once

#include "agmScalarField.h"
#include "agmVectorField.h"
#include "agmMetrics.h"
#ifdef __GNUC__
#define __cdecl __attribute__((cdecl))
#define __declspec(dllexport) __attribute__((visibility ("default")))
#endif

void _proceedGlobals(bool bGet = true);

typedef uint32_t (*PROTO_mfoWiegelmannCallback) (double dstep, int d, int nChunks, int nTasks, 
                                                      int depth, double ddL, CagmMetrics *_metrics,
                                                      CagmVectorField *, int *stop);

extern "C" {
__declspec( dllexport ) uint32_t utilInitialize();
__declspec( dllexport ) int utilSetInt(char *, int);
__declspec( dllexport ) int utilGetInt(char *, int *);
__declspec( dllexport ) int utilSetDouble(char *, double);
__declspec( dllexport ) int utilGetDouble(char *, double *);
__declspec( dllexport ) int utilSetSetting(char *, double);

__declspec( dllexport ) int utilGetVersion(char *, int);

__declspec( dllexport ) uint32_t mfoWiegelmannProcedure(CagmVectorField *, CagmScalarField *,
                                                     CagmVectorField *, CagmVectorField *, CagmVectorField *, CagmVectorField *,
                                                     CagmScalarField *, CagmScalarField *, CagmScalarField *, CagmScalarField *,
                                                     CagmScalarField *, CagmScalarField *, CagmScalarField *, CagmScalarField *,
                                                     CagmScalarField *,
                                                     double *, PROTO_mfoWiegelmannCallback);

// Lines
__declspec(dllexport) uint32_t mfoGetLinesV(CagmVectorField *v,
    uint32_t _cond, double chromoLevel,
    double *_seeds, int _Nseeds, double relSeedsBound, 
    int nProc,
    double step, double tolerance, double absBoundAchieve, double relBoundAchieve,
    int _n_loop_control, double loop_abs_err,
    int *_nLines, int *_nPassed,
    int *_voxelStatus, double *_physLength, double *_avField,
    int *_linesLength, int *_codes,
    int *_startIdx, int *_endIdx, int *_apexIdx,
    uint64_t _maxCoordLength, uint64_t *_totalLength, double *_coords, uint64_t *_linesStart, int *_linesIndex, int *seedIdx, double *times);

__declspec(dllexport) uint32_t mfoGetLines(int *N,
    double *Bx, double *By, double *Bz,
    uint32_t _cond = 0x3, double chromoLevel = 0,
    double  *_seeds = nullptr, int _Nseeds = 0, double relSeedsBound = 0,
    int nProc = 0,
    double step = 1.0, double tolerance = 1e-3, double absBoundAchieve = 1e-3, double relBoundAchieve = 0,
    int _n_loop_control = 0, double loop_abs_err = 1.5e-1,
    int *_nLines = nullptr, int *_nPassed = nullptr,
    int *_voxelStatus = nullptr, double *_physLength = nullptr, double *_avField = nullptr,
    int *_linesLength = nullptr, int *_codes = nullptr,
    int *_startIdx = nullptr, int *_endIdx = nullptr, int *_apexIdx = nullptr,
    uint64_t _maxCoordLength = 0, uint64_t *_totalLength = nullptr, double *_coords = nullptr, uint64_t *_linesStart = nullptr, int *_linesIndex = nullptr, int *seedIdx = nullptr, double *times = nullptr);

__declspec(dllexport) uint32_t mfoGetLinesP(int *N,
    double *Bx, double *By, double *Bz,
    double *seeds, int Nseeds,
    int *nLines, int *nPassed,
    int *status, double *physLength, double *avField,
    int *linesLength, int *codes,
    int *startIdx, int *endIdx, int *apexIdx,
    uint64_t maxCoordLength, uint64_t *totalLength, double *coord, uint64_t *linesStart, int *linesIndex, int *seedIdx, double *times);

// IDL interface
__declspec(dllexport) int mfoNLFFF(int argc, void* argv[]);
__declspec(dllexport) int mfoLines(int argc, void* argv[]);
__declspec(dllexport) int mfoNLFFFVersion(int argc, void* argv[]);
__declspec(dllexport) double mfoNLFFFMemoryIDL(int /* argc */, void* argv[]);

// Core interface
__declspec(dllexport) uint32_t mfoNLFFFCore(int *N, double *Bx, double *By, double *Bz);
__declspec(dllexport) uint64_t mfoNLFFFMemory(int *N, int nProc);

// Disambiguation
__declspec(dllexport) int mfoDisambigStep(int *N, int *_state, double *_Fstate, double *_F, double *_kT, int *_prescribe
    , int *_natt, int *_nacc);
__declspec(dllexport) int mfoDisambig(int *N, double *absB, double *incl, double *azim, double *Brefx, double *Brefy, double *Brefz, double *step, int *prescribed
    , int *_state, int *_generation, int *_status = nullptr);
}
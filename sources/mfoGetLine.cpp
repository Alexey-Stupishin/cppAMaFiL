#include "stdDefinitions.h"

#include "mfoGlobals.h"

#include "MagFieldOps.h"
#include "LinesProcessor.h"
#include "LinesTask.h"
#include "agmRKF45.h"

#include "console_debug.h"
#include "DebugWrite.h"

static void internal_convert_indices(CubeXD *cube,
    int *status, double *physLength, double *avField,
    int *startIdx, int *endIdx, int *apexIdx, int *seedIdx,
    int *codes, double *times)
{
    if (status)
        cube->transposeXY(status);
    if (physLength)
        cube->transposeXY(physLength);
    if (avField)
        cube->transposeXY(avField);
    if (codes)
        cube->transposeXY(codes);
    if (times)
        cube->transposeXY(times);

    if (startIdx)
        cube->transposeIDXY(startIdx);
    if (endIdx)
        cube->transposeIDXY(endIdx);
    if (apexIdx)
        cube->transposeIDXY(apexIdx);
    if (seedIdx)
        cube->transposeIDXY(seedIdx);
}

__declspec(dllexport) uint32_t mfoGetLinesV(CagmVectorField *v,
    uint32_t _cond, double chromoLevel,
    double *_seeds, int _Nseeds, double relSeedsBound,
    int nProc,
    double step, double tolerance, double absBoundAchieve, double relBoundAchieve,
    int _n_loop_control, double loop_abs_cell,
    int *_nLines, int *_nPassed,
    int *_voxelStatus, double *_physLength, double *_avField, 
    int *_linesLength, int *_codes,
    int *_startIdx, int *_endIdx, int *_apexIdx,
    uint64_t _maxCoordLength, uint64_t *_totalLength, double *_coords, uint64_t *_linesStart, int *_linesIndex, int *seedIdx, double *_times)
{
    nProc = TaskQueueProcessor::getProcInfo(nProc);
    nProc = _Nseeds > 0 ? (_Nseeds < nProc ? _Nseeds : nProc) : nProc;
    TaskQueueProcessor proc(nProc);

    int maxResult = 50000;

    LQPTaskFactory factory;
    LQPSupervisor *supervisor = new LQPSupervisor(v, _cond, chromoLevel,
        _seeds, _Nseeds, relSeedsBound, lines_use_durstenfeld,
        _voxelStatus, _physLength, _avField,
        _linesLength, _codes, _times,
        _startIdx, _endIdx, _apexIdx,
        _maxCoordLength, _coords, _linesStart, _linesIndex, seedIdx, &factory, proc.get_sync());

    std::vector<ATQPProcessor *> processors;
    for (int i = 0; i < nProc; i++)
        processors.push_back(new LQPProcessor(supervisor, i, v, _cond, chromoLevel, 0, step, tolerance, 0
            , absBoundAchieve, relBoundAchieve, maxResult, _voxelStatus, proc.get_sync(), _n_loop_control, loop_abs_cell));

    proc.proceed(processors, supervisor, w_priority::low);

    console_debug("end of proceed")

    uint64_t cumLength;
    int nLines, nPassed, nNonStored;
    supervisor->getFinalState(&cumLength, &nLines, &nPassed, &nNonStored);
    if (_totalLength)
        *_totalLength = cumLength;
    if (_nLines)
        *_nLines = nLines;
    if (_nPassed)
        *_nPassed = nPassed;

    console_debug("end of assign")

    for (int i = 0; i < nProc; i++)
        delete processors[i];

    console_debug("processors deleted")

    DebugWriteData(v, "debug_lines_field");
    DebugWriteLines(v, "debug_lines", 
        _seeds, _Nseeds,
        nLines, nPassed,
        _voxelStatus, _physLength, _avField, 
        _linesLength, _codes,
        _startIdx, _endIdx, _apexIdx,
        cumLength, _coords, _linesStart, _linesIndex, seedIdx);

    delete supervisor;

    console_debug("supervisor deleted")

    return nNonStored;
}

__declspec(dllexport) uint32_t mfoGetLines(int *N,
    double *Bx, double *By, double *Bz,
    uint32_t conditions, double chromo_level,
    double *seeds, int Nseeds, double relSeedsBound,
    int nProc,
    double step, double tolerance, double absBoundAchieve, double relBoundAchieve,
    int _n_loop_control, double loop_abs_cell,
    int *nLines, int *nPassed,
    int *status, double *physLength, double *avField, 
    int *linesLength, int *codes,
    int *startIdx, int *endIdx, int *apexIdx,
    uint64_t maxCoordLength, uint64_t *totalLength, double *coord, uint64_t *linesStart, int *linesIndex, int *seedIdx, double *times)
{
    console_start();

    CagmVectorField *v;
    if (debug_input)
        v = new CagmVectorField(N, Bx, By, Bz);
    else
        v = new CagmVectorField(N, Bx, By, Bz, true);

    uint32_t rc = mfoGetLinesV(v,
        conditions, chromo_level,
        seeds, Nseeds, relSeedsBound, 
        nProc,
        step, tolerance, absBoundAchieve, relBoundAchieve,
        _n_loop_control, loop_abs_cell,
        nLines, nPassed,
        status, physLength, avField,
        linesLength, codes,
        startIdx, endIdx, apexIdx,
        maxCoordLength, totalLength, coord, linesStart, linesIndex, seedIdx, times);

    delete v;

    console_debug("vector box deleted")

    return rc;
}

__declspec(dllexport) uint32_t mfoGetLinesP(int *N,
    double *Bx, double *By, double *Bz,
    double *seeds, int Nseeds,
    int *nLines, int *nPassed,
    int *status, double *physLength, double *avField, 
    int *linesLength, int *codes,
    int *startIdx, int *endIdx, int *apexIdx,
    uint64_t maxCoordLength, uint64_t *totalLength, double *coord, uint64_t *linesStart, int *linesIndex, int *seedIdx, double *times)
{
    console_start();

    CagmVectorField *v;
    if (debug_input)
        v = new CagmVectorField(N, Bx, By, Bz);
    else
        v = new CagmVectorField(N, Bx, By, Bz, true);

    uint32_t rc = mfoGetLinesV(v,
        lines_conditions, lines_chromo_level,
        seeds, Nseeds, lines_rel_seeds_bound, 
        CommonThreadsN,
        lines_step, lines_tolerance, lines_abs_bound_achieve, lines_rel_bound_achieve,
        lines_n_loop_control, lines_loop_abs_cell,
        nLines, nPassed,
        status, physLength, avField,
        linesLength, codes,
        startIdx, endIdx, apexIdx,
        maxCoordLength, totalLength, coord, linesStart, linesIndex, seedIdx, times);

    if (!seeds && lines_internal_convert_indices)
        internal_convert_indices(v, 
                                 status, physLength, avField,
                                 startIdx, endIdx, apexIdx, seedIdx,
                                 codes, times);

    delete v;

    console_debug("vector box deleted")

    return rc;
}

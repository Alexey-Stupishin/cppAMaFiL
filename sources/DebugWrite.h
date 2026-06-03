#pragma once

#include "stdDefinitions.h"
#include "agmScalarField.h"
#include "agmVectorField.h"
#include "debug_data_trace_win.h"

void DebugWritePars(const char *fname, CagmVectorField * field, CagmScalarField * w);
void DebugWriteData(CubeXD *v, const char *fname, int depth = 0, int iter = 0);
void DebugWriteLines(CubeXD *v, const char *fname, 
    double *_seeds, int _Nseeds,
    int _nLines, int _nPassed,
    int *_voxelStatus, double *_physLength, double *_avField, 
    int *_linesLength, int *_codes,
    int *_startIdx, int *_endIdx, int *_apexIdx,
    uint64_t _totalLength, double *_coords, uint64_t *_linesStart, int *_linesIndex, int *seedIdx);

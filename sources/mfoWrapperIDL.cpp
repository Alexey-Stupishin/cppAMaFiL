#include <stdio.h>
#include "mfoGlobals.h"
#include "MagFieldOps.h"
#include "agpWiegelmannPar.h"
#include "availMem.h"

#include "console_debug.h"

#include "idl_export_ext.h"

int mfoNLFFFVersion(int /* argc */, void* argv[])
{
    IDL_STRING *idlstr = (IDL_STRING *)argv[0];
    char *buffer = new char[512];
    int L = utilGetVersion(buffer, 512);
    strcpy(idlstr->s, buffer);
    idlstr->slen = (IDL_STRING_SLEN_T)strlen(buffer);

    delete [] buffer;

    return L;
}

double mfoNLFFFMemoryIDL(int /* argc */, void* argv[])
{
    int c = 0;
    int *N = (int *)argv[c++];
    int n_proc = *(int *)argv[c++];
    uint64_t *p_est = (uint64_t *)argv[c++];
    uint64_t *p_avail = (uint64_t *)argv[c++];
    *p_est = CagpWiegelmann::estimateMemory(N, n_proc);
    *p_avail = availableMemory();

    return (double)*p_est / (double)*p_avail;
}

int mfoNLFFF(int argc, void* argv[])
{
    utilInitialize();

    CidlPassParameterMap *m = (CidlPassParameterMap *)argv[0];
    if (m->parse(&mapInt, &mapuint64_t, &mapDouble) != 0)
        return LIB_STATE_NO_PARAMETER;
    _proceedGlobals();

    int *N = (int *)argv[1];

    CagmVectorField *v = new CagmVectorField(N, (double *)argv[2], (double *)argv[3], (double *)argv[4]);
    // make potential!?

    // weight
    CagmScalarField *sW = new CagmScalarField(N);
    // WiegelmannWeightType == SWF_COS!
    int xBounds[2], yBounds[2], zBounds[2];
    sW->Weight(1, WiegelmannWeightBound, xBounds, yBounds, zBounds);

    // abs cond, weight (5,6)
    // los cond, weight (7,8, 9)
    // base cond, weight (10, 11), (12,13), (14,15)
    // abs2 cond, weight2 (16,17)
    // los2 cond, weight2 (18,19)
    // base2 cond, weight2 (20,21), (22,23), (24,25)
    double *fx = nullptr, *wx = nullptr, *fy = nullptr, *wy = nullptr, *fz = nullptr, *wz = nullptr, *fl = nullptr, *wl = nullptr, *fa = nullptr, *wa = nullptr;
    double *fx2 = nullptr, *wx2 = nullptr, *fy2 = nullptr, *wy2 = nullptr, *fz2 = nullptr, *wz2 = nullptr, *fl2 = nullptr, *wl2 = nullptr, *fa2 = nullptr, *wa2 = nullptr;
    double vc[] = {0, 0, 1};

    if (argc > 6 && argv[5] && argv[6])
    {
        fa = (double *)(argv[5]);
        wa = (double *)(argv[6]);
    }

    if (argc > 8 && argv[7] && argv[8])
    {
        fl = (double *)(argv[7]);
        wl = (double *)(argv[8]);
        if (argc > 9 && argv[9])
        {
            vc[0] = ((double *)argv[9])[0];
            vc[1] = ((double *)argv[9])[1];
            vc[2] = ((double *)argv[9])[2];
        }
    }

    if (argc > 11 && argv[10] && argv[11])
    {
        fx = (double *)(argv[10]);
        wx = (double *)(argv[11]);
    }

    if (argc > 13 && argv[12] && argv[13])
    {
        fy = (double *)(argv[12]);
        wy = (double *)(argv[13]);
    }

    if (argc > 15 && argv[14] && argv[15])
    {
        fz = (double *)(argv[14]);
        wz = (double *)(argv[15]);
    }

    if (argc > 17 && argv[16] && argv[17])
    {
        fa2 = (double *)(argv[16]);
        wa2 = (double *)(argv[17]);
    }

    if (argc > 19 && argv[18] && argv[19])
    {
        fl2 = (double *)(argv[18]);
        wl2 = (double *)(argv[19]);
    }

    if (argc > 21 && argv[20] && argv[21])
    {
        fx2 = (double *)(argv[20]);
        wx2 = (double *)(argv[21]);
    }

    if (argc > 23 && argv[22] && argv[23])
    {
        fy2 = (double *)(argv[22]);
        wy2 = (double *)(argv[23]);
    }

    if (argc > 25 && argv[24] && argv[25])
    {
        fz2 = (double *)(argv[24]);
        wz2 = (double *)(argv[25]);
    }

    // argv[26] reserved for user-defined weight

    CagmScalarField *absConditions = nullptr, *absWeight = nullptr;
    CagmScalarField *absConditions2 = nullptr, *absWeight2 = nullptr;
    if (fa)
    {
        absConditions = new CagmScalarField(N);
        absConditions->setField(fa);
        absWeight = new CagmScalarField(N);
        absWeight->setField(wa);
    }
    if (fa2)
    {
        absConditions2 = new CagmScalarField(N);
        absConditions2->setField(fa2);
        absWeight2 = new CagmScalarField(N);
        absWeight2->setField(wa2);
    }

    CagmScalarField *losConditions = nullptr, *losWeight = nullptr;
    CagmScalarField *losConditions2 = nullptr, *losWeight2 = nullptr;
    if (fl)
    {
        losConditions = new CagmScalarField(N);
        losConditions->setField(fl);
        losWeight = new CagmScalarField(N);
        losWeight->setField(wl);
    }
    if (fl2)
    {
        losConditions2 = new CagmScalarField(N);
        losConditions2->setField(fl2);
        losWeight2 = new CagmScalarField(N);
        losWeight2->setField(wl2);
    }

    CagmVectorField *baseConditions = nullptr, *baseWeight = nullptr;
    CagmVectorField *baseConditions2 = nullptr, *baseWeight2 = nullptr;
    if (fx || fy || fz)
    {
        baseConditions = new CagmVectorField(N);
        baseConditions->zero();
        baseWeight = new CagmVectorField(N);
        baseWeight->zero();
        if (fx)
        {
            baseConditions->CopyComp(fx, 0);
            baseWeight->CopyComp(wx, 0);
        }
        if (fy)
        {
            baseConditions->CopyComp(fy, 1);
            baseWeight->CopyComp(wy, 1);
        }
        if (fz)
        {
            baseConditions->CopyComp(fz, 2);
            baseWeight->CopyComp(wz, 2);
        }
    }
    if (fx2 || fy2 || fz2)
    {
        baseConditions2 = new CagmVectorField(N);
        baseConditions2->zero();
        baseWeight2 = new CagmVectorField(N);
        baseWeight2->zero();
        if (fx2)
        {
            baseConditions2->CopyComp(fx2, 0);
            baseWeight2->CopyComp(wx2, 0);
        }
        if (fy2)
        {
            baseConditions2->CopyComp(fy2, 1);
            baseWeight2->CopyComp(wy2, 1);
        }
        if (fz2)
        {
            baseConditions2->CopyComp(fz2, 2);
            baseWeight2->CopyComp(wz2, 2);
        }
    }

    // Callback, 27
    PROTO_mfoWiegelmannCallback callback = nullptr;
    if (argc > 27 && argv[27])
        callback = (PROTO_mfoWiegelmannCallback)(argv[27]);

    mfoWiegelmannProcedure(v, sW, baseConditions, baseWeight, baseConditions2, baseWeight2, 
        absConditions, absWeight, absConditions2, absWeight2, 
        losConditions, losWeight, losConditions2, losWeight2, nullptr, vc, callback);
    v->GetComp((double *)argv[2], 0);
    v->GetComp((double *)argv[3], 1);
    v->GetComp((double *)argv[4], 2);

    return 0;
}

int mfoLines(int /* argc */, void* argv[])
{
    utilInitialize();

    int c = 0;
    CidlPassParameterMap *m = (CidlPassParameterMap *)argv[c++];
    if (m->parse(&mapInt, &mapuint64_t, &mapDouble) != 0)
        return LIB_STATE_NO_PARAMETER;
    _proceedGlobals();

    
    int *N = (int *)argv[c++];
    double *Bx = (double *)argv[c++];
    double *By = (double *)argv[c++];
    double *Bz = (double *)argv[c++];

    int *status = (int *)argv[c++];
    double *physLength = (double *)argv[c++];
    double *avField = (double *)argv[c++];
    int *startIdx = (int *)argv[c++];
    int *endIdx = (int *)argv[c++];
    int *apexIdx = (int *)argv[c++];
    int *seedIdx = (int *)argv[c++];

    double *seeds = (double *)argv[c++];
    int *Nseeds_a = (int *)argv[c++];
    int Nseeds = 0;
    if (Nseeds_a)
        Nseeds = *Nseeds_a;
    uint64_t *maxLength_a = (uint64_t *)argv[c++];
    uint64_t maxLength = 0;
    if (maxLength_a)
        maxLength = *maxLength_a;

    uint64_t *totalLength = (uint64_t *)argv[c++];
    int *nLines = (int *)argv[c++];
    int *nPassed = (int *)argv[c++];
    double *coord = (double *)argv[c++];
    uint64_t *linesStart = (uint64_t *)argv[c++];
    int *linesLength = (int *)argv[c++];
    int *linesIndex = (int *)argv[c++];

    int *codes = (int *)argv[c++];
    double *times = (double *)argv[c++];

    return mfoGetLinesP(N, Bx, By, Bz,
        seeds, Nseeds,
        nLines, nPassed,
        status, physLength, avField,
        linesLength, codes,
        startIdx, endIdx, apexIdx,
        maxLength, totalLength, coord, linesStart, linesIndex, seedIdx, times);
}
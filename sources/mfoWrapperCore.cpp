#include "mfoGlobals.h"
#include "MagFieldOps.h"
#include "agpWiegelmannPar.h"

uint32_t mfoNLFFFCore(int *N, double *Bx, double *By, double *Bz)
{
    CagmVectorField *v = new CagmVectorField(N, Bx, By, Bz);
    // weight
    CagmScalarField *sW = new CagmScalarField(N); // WiegelmannWeightType == SWF_COS!
    int xBounds[2], yBounds[2], zBounds[2];
    sW->Weight(1, WiegelmannWeightBound, xBounds, yBounds, zBounds);
    double vc[] = {0, 0, 1};

    mfoWiegelmannProcedure(v, sW, nullptr, nullptr, nullptr, nullptr, 
            nullptr, nullptr, nullptr, nullptr, 
            nullptr, nullptr, nullptr, nullptr, nullptr, vc, nullptr);

    v->GetComp(Bx, 0);
    v->GetComp(By, 1);
    v->GetComp(Bz, 2);

    return 0;
}

uint64_t mfoNLFFFMemory(int *N, int nProc)
{
    return CagpWiegelmann::estimateMemory(N, nProc, true);
}

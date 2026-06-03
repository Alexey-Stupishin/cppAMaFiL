#include "stdDefinitions.h"
#include "MagFieldOps.h"
#include "mfoGlobals.h"
#include "agpWiegelmannPar.h"
#include "agmRotate3D.h"
#include "TimeTicToc.h"

#include "console_debug.h"
#include "debug_data_trace_win.h"
#include "DebugWrite.h"
#include "WndDebug.h"
#include "WProcCore.h"

uint32_t mfoWiegelmannProcedureCore(CagmVectorField *field, CagmScalarField *weight, 
    CagmVectorField *baseField, CagmVectorField *baseWeight, CagmVectorField *baseField2, CagmVectorField *baseWeight2,
    CagmScalarField *absField, CagmScalarField *absWeight, CagmScalarField *absField2, CagmScalarField *absWeight2,
    CagmScalarField *losField, CagmScalarField *losWeight, CagmScalarField *losField2, CagmScalarField *losWeight2,
    CagmScalarField *bottomWeight,
    double *vcos, int depth, PROTO_mfoWiegelmannCallback callback, double *maxStep)
{
    *maxStep = 0;

    CWProcCore core(field, weight, 
    baseField, baseWeight, baseField2, baseWeight2,
    absField, absWeight, absField2, absWeight2,
    losField, losWeight, losField2, losWeight2,
    bottomWeight,
    vcos, depth, callback);

    int reason  = core.proceed();

    return reason;
}

static int xfloor(double v)
{
    int c = (int)ceil(v);
    if (c - v < 0.5)
        return c;
    else
        return c-1;
}

__declspec( dllexport ) uint32_t mfoWiegelmannProcedure(CagmVectorField *field, CagmScalarField *weight, 
    CagmVectorField *baseField, CagmVectorField *baseWeight, CagmVectorField *baseField2, CagmVectorField *baseWeight2,
    CagmScalarField *absField, CagmScalarField *absWeight, CagmScalarField *absField2, CagmScalarField *absWeight2, 
    CagmScalarField *losField, CagmScalarField *losWeight, CagmScalarField *losField2, CagmScalarField *losWeight2,
    CagmScalarField *bottomWeight,
    double *vcos, PROTO_mfoWiegelmannCallback callback)
{
#ifdef _DEBUG_MEMORY
    debugtrace = new BWndDebug("C:\\temp\\W2Memory.log", "wb", TRUE, NULL, FALSE, TRUE);
    debugtrace->memory("Start");
#endif

    DebugWritePars("debug_input", field, weight);
    DebugWriteData(field, "debug_B");
    DebugWriteData(weight, "debug_W");

    console_start();

    if (!baseField || !baseWeight)
        WiegelmannProcCondBase = 0;
    if (!baseField2 || !baseWeight2)
        WiegelmannProcCondBase2 = 0;
    if (!absField || !absWeight)
        WiegelmannProcCondAbs = 0;
    if (!absField2 || !absWeight2)
        WiegelmannProcCondAbs2 = 0;
    if (!losField || !losWeight)
        WiegelmannProcCondLOS = 0;
    if (!losField2 || !losWeight2)
        WiegelmannProcCondLOS2 = 0;
    if (WiegelmannProcCondBase == 0 && WiegelmannProcCondBase2 == 0 &&
        WiegelmannProcCondAbs == 0 && WiegelmannProcCondAbs2 == 0 &&
        WiegelmannProcCondLOS == 0 && WiegelmannProcCondLOS2 == 0)
        WiegelmannProcCondType = 0;

    double maxStep;
    bool bMatr = WiegelmannMatryoshkaUse != 0;
    int depth = 1;
    if (bMatr)
    {
        depth = field->getMatryoshkaDepth(WiegelmannMatryoshkaDeepMinN, WiegelmannMatryoshkaFactor);
        if (depth < 2)
            bMatr = false;
    }

    uint32_t dwRC = 0;
    if (!bMatr)
    {
        WiegelmannProcStepIncr = WiegelmannProcStepIncrMain;
        WiegelmannProcStepDecr = WiegelmannProcStepDecrMain;
        WiegelmannProcStepLim = WiegelmannProcStepLimMain;
        WiegelmannProcdLStdVal = WiegelmannProcdLStdValMain;
        WiegelmannProcdLStdWin = WiegelmannProcdLStdWinMain;
        dwRC = mfoWiegelmannProcedureCore(field, weight, baseField, baseWeight, baseField2, baseWeight2,
            absField, absWeight, absField2, absWeight2,
            losField, losWeight, losField2, losWeight2, bottomWeight, vcos, 1, callback, &maxStep);
    }
    else
    {
#ifdef _DEBUG_MEMORY
    debugtrace->Say("*** depth %d", depth);
#endif
        double factor = WiegelmannMatryoshkaFactor;
        int *Ns = new int[3*depth];
        double *steps = new double[3*depth];
        field->dimensions(Ns);
        steps[0] = 1.0;
        steps[1] = 1.0;
        steps[2] = 1.0;
        for (int i = 1; i < depth; i++)
        {
            Ns[3*i  ] = xfloor( (Ns[3*(i-1)  ]-1.0)/factor + 1.0 ); 
            Ns[3*i+1] = xfloor( (Ns[3*(i-1)+1]-1.0)/factor + 1.0 ); 
            Ns[3*i+2] = xfloor( (Ns[3*(i-1)+2]-1.0)/factor + 1.0 ); 
        }
        double cf = factor;
        for (int i = 1; i < depth; i++)
        {
            steps[3*i  ] = (Ns[3*i  ]-1.0)/(Ns[0]-1.0);
            steps[3*i+1] = (Ns[3*i+1]-1.0)/(Ns[1]-1.0);
            steps[3*i+2] = (Ns[3*i+2]-1.0)/(Ns[2]-1.0);
            cf *= factor;
        }

        WiegelmannProcStepIncr = WiegelmannProcStepIncrInit;
        WiegelmannProcStepDecr = WiegelmannProcStepDecrInit;
        WiegelmannProcStepLim = WiegelmannProcStepLimInit;
        WiegelmannProcdLStdVal = WiegelmannProcdLStdValInit;
        WiegelmannProcdLStdWin = WiegelmannProcdLStdWinInit;

        int *pN0 = Ns+3*(depth-1);
        CagmVectorField *v0 = new CagmVectorField(pN0, steps+3*(depth-1));

        DebugWriteData(field, "field_init", depth);
        v0->stretch(field);
        DebugWriteData(v0, "field_stretch", depth);

        DebugWriteData(weight, "weight_init", depth);
        CagmScalarField *sW0 = new CagmScalarField(pN0, steps+3*(depth-1));
        sW0->stretch(weight);
        DebugWriteData(sW0, "weight_stretch", depth);

        mfoWiegelmannProcedureCore(v0, sW0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, vcos, depth, callback, &maxStep);
        delete sW0;
        CagmVectorField *v1;
        for (int i = depth-2; i > 0; i--)
        {
#ifdef _DEBUG_MEMORY
    debugtrace->Say("*** depth %d", i+1);
    debugtrace->memory("  pre core");
#endif
            WiegelmannProcStepIncr = WiegelmannProcStepIncrMatr;
            WiegelmannProcStepDecr = WiegelmannProcStepDecrMatr;
            WiegelmannProcStepLim = WiegelmannProcStepLimMatr;
            WiegelmannProcdLStdVal = WiegelmannProcdLStdValMatr;
            WiegelmannProcdLStdWin = WiegelmannProcdLStdWinMatr;
            pN0 = Ns+3*i;
            v1 = new CagmVectorField(pN0, steps+3*i);
            v1->stretch(v0);
            DebugWriteData(v1, "field_stretch", i+1);

            delete v0;
            sW0 = new CagmScalarField(pN0, steps+3*i);
            sW0->stretch(weight);
            DebugWriteData(sW0, "weight_stretch", i+1);
            mfoWiegelmannProcedureCore(v1, sW0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, vcos, i+1, callback, &maxStep);
            delete sW0;
            v0 = v1;
#ifdef _DEBUG_MEMORY
    debugtrace->memory("  end core");
#endif
        }

#ifdef _DEBUG_MEMORY
    debugtrace->Say("*** depth %d", 1);
#endif
        int NB[3];
        field->dimensions(NB);
        NB[0] = 2;
        CagmVectorField boundsx(NB);
        field->dimensions(NB);
        NB[1] = 2;
        CagmVectorField boundsy(NB);
        field->dimensions(NB);
        NB[2] = 2;
        CagmVectorField boundsz(NB);
        field->getBounds(&boundsx, &boundsy, &boundsz);
        field->stretch(v0);
        field->setBounds(&boundsx, &boundsy, &boundsz);
        DebugWriteData(field, "field_stretch", 1);

        delete v0;

        delete [] Ns;
        delete [] steps;

        WiegelmannProcStepIncr = WiegelmannProcStepIncrMain;
        WiegelmannProcStepDecr = WiegelmannProcStepDecrMain;
        WiegelmannProcStepLim = WiegelmannProcStepLimMain;

        WiegelmannProcdLStdVal = WiegelmannProcdLStdValMain;
        WiegelmannProcdLStdWin = WiegelmannProcdLStdWinMain;

        dwRC = mfoWiegelmannProcedureCore(field, weight, baseField, baseWeight, baseField2, baseWeight2, 
                                          absField, absWeight, absField2, absWeight2, 
                                          losField, losWeight, losField2, losWeight2, bottomWeight, vcos, 1, callback, &maxStep);
    }

#ifdef _DEBUG_MEMORY
    debugtrace->memory("End");
#endif

    return dwRC;
}

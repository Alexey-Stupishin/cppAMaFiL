#include "stdDefinitions.h"

//#include "MagFieldOps.h"
#include "mfoGlobals.h"
#include "agpWiegelmannPar.h"
#include "agmRotate3D.h"
#include "TimeTicToc.h"

//#include "console_debug.h"
//#include "debug_data_trace_win.h"
#include "DebugWrite.h"
//#include "WndDebug.h"

#include "WProcCore.h"

//--------------------------------------------------------------------
double CWProcCore::get_max_L_incr()
{
    return DBL_EPSILON * sqrt(N[0] * N[1] * N[2]);
}

//--------------------------------------------------------------------
int CWProcCore::proceedDL(double relL)
{
    if (stepN == 0)
        memoryAv[stepN] = 0.0;
    else if (stepN < WiegelmannProcdLStdWin)
        memoryAv[stepN] = relL;
    else
    {
        for (int k = 1; k < WiegelmannProcdLStdWin; k++)
            memoryAv[k-1] = memoryAv[k];
        memoryAv[WiegelmannProcdLStdWin-1] = relL;

        double s = 0, s2 = 0;
        for (int k = 0; k < WiegelmannProcdLStdWin; k++)
        {
            s += memoryAv[k];
            s2 += memoryAv[k]*memoryAv[k];
        }

        double w = WiegelmannProcdLStdWin;
        double var = sqrt(w/(w-1)*(w*s2/(s*s)-1));
        if (var < WiegelmannProcdLStdVal)
            return 7;
    }

     return 0;
}

//--------------------------------------------------------------------
double CWProcCore::init_step()
{
    FILE *fid = fopen("c:\\temp\\start.txt", "wb");
    fprintf(fid, "%f %f\r\n", step, L0);

    double xL = 0, xC = step, xR = step;
    double yL = L0, yC = L0, yR = L0;
    bool bLeft = false, bRight = false;
    double factor = 1.1;
    while (true)
    {
        //restore_field();
        next_field();

        L = proc->step(iterN);
        iterN++;

        fprintf(fid, "%f %f\r\n", step, L);

        if (L < yL)
        {
            if (bRight)
            { xC = step; yC = L; break; }
            if (!bLeft)
                bLeft = true;
            else
            { xL = xC; yL = yC; }
            xC = step; yC = L;
            step *= factor;
            store_field();
        }
        else
        {
            if (bLeft)
            { xR = step; yR = L; break; }
            if (!bRight)
                bRight = true;
            xR = step; yR = L;
            step /= factor;
            store_field();
        }
    }
    fclose(fid);

    double v1, v2, v3;
    v1 = yL*(xC-xR);
    v2 = yC*(xR-xL);
    v3 = yR*(xL-xC);

    double mba = (v1*(xC + xR) + v2*(xR + xL) + v3*(xL + xC)) / (v1+v2+v3);

    mba *= WiegelmannProcStepDecr;

    restore_field();

    return mba;
}

//--------------------------------------------------------------------
void CWProcCore::next_field()
{
    vF->mult(step);
    field->add(vF);
    field->setBounds(boundsx, boundsy, boundsz);
}

//--------------------------------------------------------------------
void CWProcCore::restore_field()
{
    field->Copy(prevField);
    vF->Copy(prevVF);
}

//--------------------------------------------------------------------
void CWProcCore::store_field()
{
    prevField->Copy(field);
    prevVF->Copy(vF);
}

//--------------------------------------------------------------------
CWProcCore::CWProcCore(CagmVectorField *_field, CagmScalarField *_weight, 
    CagmVectorField *_baseField, CagmVectorField *_baseWeight, CagmVectorField *_baseField2, CagmVectorField *_baseWeight2,
    CagmScalarField *_absField, CagmScalarField *_absWeight, CagmScalarField *_absField2, CagmScalarField *_absWeight2,
    CagmScalarField *_losField, CagmScalarField *_losWeight, CagmScalarField *_losField2, CagmScalarField *_losWeight2,
    CagmScalarField *_bottomWeight,
    double *vcos, int _depth, PROTO_mfoWiegelmannCallback _callback)
    : field(_field)
    , weight(_weight)
    , baseField(_baseField)
    , baseWeight(_baseWeight)
    , baseField2(_baseField2)
    , baseWeight2(_baseWeight2)
    , absField(_absField)
    , absWeight(_absWeight)
    , absField2(_absField2)
    , absWeight2(_absWeight2)
    , losField(_losField)
    , losWeight(_losWeight)
    , losField2(_losField2)
    , losWeight2(_losWeight2)
    , bottomWeight(_bottomWeight)
    , depth(_depth)
    , callback(_callback)
    , stepN(0)
    , stop(0)
    , iterN(0)
    , reason(0)
    , curr_incr(WiegelmannProcStepIncr)
    , curr_decr(WiegelmannProcStepDecr)
    , step0(WiegelmannProcStep0)
    , step(WiegelmannProcStep0)
{
    field->dimensions(N);

    memoryAv = new double[WiegelmannProcdLStdWin];
    memory_dL = new double[WiegelmannProcdLStdWin];
    n_appr = 0;

    vF = new CagmVectorField(field);

    int NB[3];
    memcpy(NB, N, 3*sizeof(int));
    NB[0] = 2;
    boundsx = new CagmVectorField (NB);
    memcpy(NB, N, 3*sizeof(int));
    NB[1] = 2;
    boundsy = new CagmVectorField (NB);
    memcpy(NB, N, 3*sizeof(int));
    NB[2] = 2;
    boundsz = new CagmVectorField (NB);
    field->getBounds(boundsx, boundsy, boundsz);

    rotator = new CagmRotate3D(vcos);

    proc = new CagpWiegelmann(N, CommonThreadsN, WiegelmannDerivStencil
        , field, weight
        , vcos
        , vF
        , baseField, baseWeight
        , absField, absWeight
        , losField, losWeight
        , bottomWeight
        , depth
        , WiegelmannThreadPriority
        );

    DebugWriteData(field, "debug_B", depth);

    proc->max_dL_incr = get_max_L_incr();

    proceedDL(0.0);

    prevField = new CagmVectorField(N);
    prevVF = new CagmVectorField(N);

    z_size = N[2];
    voxels = N[0]*N[1]*N[2];

    proc->step_failed[0] = 0;
}

//--------------------------------------------------------------------
CWProcCore::~CWProcCore()
{
    delete prevField;
    delete prevVF;
    delete proc;
    delete rotator;
    delete vF;
    delete[] memoryAv;
    delete[] memory_dL;
}

//--------------------------------------------------------------------
int CWProcCore::proceed()
{
    L0 = proc->step(iterN);
    Lprev = L0;
    iterN++;
    DebugWriteData(vF, "debug_vF", depth);

    store_field();
    step0 = init_step();
    step = step0;

    if (WiegelmannProcCondType == 2)
        field->relax(baseField, baseWeight);

    L0 = proc->step(iterN);
    Lprev = L0;
    store_field();

    tic.tic();
    while (true)
    {
        next_field();
        if (WiegelmannProcCondType == 2)
            field->relax(baseField, baseWeight);

        L = proc->step(iterN);
        iterN++;

        dL = Lprev - L;
        if (dL/Lprev <= - proc->max_dL_incr)
        {
            restore_field();
            if (step < step0*WiegelmannProcStepLim)
                reason = 1;
            else
            {
                step *= curr_decr;
                curr_incr = 1 + (curr_incr-1)*WiegelmannProcStepIncrFactor;
                curr_decr = 1 - (1-curr_decr)*WiegelmannProcStepDecrFactor;
            }
            (proc->step_failed[0])++;
        }
        else
        {
            stepN++;

            proc->setBase(L/L0, step, stepN, depth, iterN, tic.toc(), z_size, voxels);

            Lprev = L;
            store_field();

            step *= curr_incr;
            if (step > step0*WiegelmannProcStepMax)
                step = step0*WiegelmannProcStepMax;

            if (stepN > WiegelmannProcMaxSteps)
                reason = 3;
            else
                reason = proceedDL(L/L0);

            if (callback && stepN%WiegelmannProtocolStep == 0)
                callback(step/step0, 0, 0, 0, depth, dL/L0/step, proc, field, &stop);

            proc->step_failed[0] = 0;
            // theta etc. reasons from 8
        }

        if (reason != 0)
        {
            restore_field();
            break;
        }

        if (stop)
        {
            reason = 4;
            break;
        }
    }

    field->setBounds(boundsx, boundsy, boundsz);
    if (WiegelmannProcCondType == 2)
        field->relax(baseField, baseWeight);

    stop = reason;
    proc->setBase(L/L0, step, stepN, depth, iterN, tic.toc(), z_size, voxels);
    if (callback)
        callback(step/step0, 0, 0, 0, depth, dL/L0, proc, field, &stop);

    return reason;
}

#include "stdDefinitions.h"
#include <math.h>

#include "agpWiegelmannPar.h"
#include "agmScalarField.h"
#include "agmVectorField.h"

#include "WndDebug.h"
#include "DebugWrite.h"

//-----------------------------------------------------------------------
CagpWiegelmann::CagpWiegelmann(int *_N, int _n_threads, int _stencil
    , CagmVectorField *_sourceB, CagmScalarField *_sourceW
    , double *_vcos
    , CagmVectorFieldOps *_outF
    , CagmVectorField *_baseField, CagmVectorField *_baseWeight
    , CagmScalarField *_absField, CagmScalarField *_absWeight
    , CagmScalarField *_losField, CagmScalarField *_losWeight
    , CagmScalarField *_bottomWeight
    , int _depth
    , w_priority _priority
    )
    : CagmMetrics(1, _N[2], true)
    , vB(_sourceB)
    , sW(_sourceW)
    , vdircos(_vcos)
    , vF(_outF)
    , baseField(_baseField)
    , baseWeight(_baseWeight)
    , absField(_absField)
    , absWeight(_absWeight)
    , losField(_losField)
    , losWeight(_losWeight)
    , constrN(0)
    , bottomWeight(_bottomWeight)
    , Bmod(nullptr)
    , Jmod(nullptr)
    , Binv(nullptr)
    , Jinv(nullptr)
    , JxB(nullptr)
    , st(nullptr)
    , depth(_depth)
    , priority(_priority)
{
    vB->dimensions(N);
#ifdef _DEBUG_MEMORY
    debugtrace->memory("    start alloc");
    debugtrace->Say("    size: %I64d voxels", (uint64_t)N[0]*N[1]*N[2]);
#endif

    vgradW = new CagmVectorField(sW);
#ifdef _DEBUG_MEMORY
    debugtrace->memory("      gradient_W");
#endif

    DebugWriteData(sW, "debug_W", depth);

    vgradW->grad(sW);

    DebugWriteData(vgradW, "debug_gradW", depth);

    int Nb[3];
    Nb[0] = N[0];
    Nb[1] = N[1];
    Nb[2] = 1;
    vBottom = new CagmVectorField(Nb);

    rotB = new CagmVectorField(vB);
#ifdef _DEBUG_MEMORY
    debugtrace->memory("      (v)rotor");
#endif
    divB = new CagmScalarField(vB);
#ifdef _DEBUG_MEMORY
    debugtrace->memory("      (s)divergence");
#endif
    Wa = new CagmVectorField(vB);
#ifdef _DEBUG_MEMORY
    debugtrace->memory("      (v)Wa");
#endif
    Wb = new CagmVectorField(vB);
#ifdef _DEBUG_MEMORY
    debugtrace->memory("      (v)Wb");
#endif
    Wa2Wb2 = new CagmScalarField(vB);
#ifdef _DEBUG_MEMORY
    debugtrace->memory("      (s)Wa2Wb2");
#endif
    WaxB = new CagmVectorField(vB);
#ifdef _DEBUG_MEMORY
    debugtrace->memory("      (v)WaxB");
#endif
    WbxB = new CagmScalarField(vB);
#ifdef _DEBUG_MEMORY
    debugtrace->memory("      (s)WbxB");
#endif

    if (WiegelmannGetMetricsTheta)
    {
        st  =  new CagmScalarField(vB);
        Binv = new CagmScalarField(vB);
        Binv->setTolerances(WiegelmannInversionTolerance, WiegelmannInversionDenom);
        Jinv = new CagmScalarField(vB);
        Jinv->setTolerances(WiegelmannInversionTolerance, WiegelmannInversionDenom);
        Bmod = new CagmScalarField(vB);
        Jmod = new CagmScalarField(vB);
        JxB  = new CagmScalarField(vB);
    }

    Lfunc = new double[N[2]];
    F2max = new double[N[2]];
    B2sum = new double[N[2]];
#ifdef _DEBUG_MEMORY
    debugtrace->memory("      funcs");
#endif

    main_proc = new TaskQueueProcessor(_n_threads);
    nProc = main_proc->get_num_proc();

    factory = new PCOTaskFactory();
    supervisor = new PCOSupervisor(N, factory, main_proc->get_sync());

    for (int i = 0; i < nProc; i++)
        processors.push_back(new PCOProcessor(i, main_proc->get_sync(), this, N, _stencil));

#ifdef _DEBUG_MEMORY
    debugtrace->memory("    allocated");
#endif

    if (baseField)
    {
        int Nc[3];
        baseField->dimensions(Nc);
        constrN = Nc[2];
    }
}

//-----------------------------------------------------------------------
CagpWiegelmann::~CagpWiegelmann()
{
    for (int i = 0; i < nProc; i++)
        delete processors[i];

    delete supervisor;
    delete factory;
    delete main_proc;

    delete vgradW;

    delete vBottom;

    delete rotB;
    delete divB;
    delete Wa;
    delete Wb;
    delete WaxB;
    delete WbxB;

    delete Binv;
    delete Jinv;
    delete Bmod;
    delete Jmod;
    delete JxB;
    delete st;

    delete [] F2max;
    delete [] B2sum;

#ifdef _DEBUG_MEMORY
    debugtrace->memory("    deleted");
#endif
}

//-----------------------------------------------------------------------
bool PCOProcessor::setTask(ATQPTask * _task)
{
    bool finish = ATQPProcessor::setTask(_task);
    if (!finish)
        task = (PCOTask *)_task;

    return finish;
}

//-----------------------------------------------------------------------
double CagpWiegelmann::step(int _iterN)
{
    iterN = _iterN;
    counter++;

    //DebugWriteData(vB, "vB", depth, iterN);

    mode = 0;
    main_proc->proceed(processors, supervisor, priority);

    if (depth == 1 && iterN == 0)
    { 
        DebugWriteData(rotB, "rotB", depth, iterN);
        DebugWriteData(Wa, "Wa", depth, iterN);
        DebugWriteData(divB, "divB", depth, iterN);
        DebugWriteData(Wb, "Wb", depth, iterN);
        DebugWriteData(WaxB, "WaxB", depth, iterN);
        DebugWriteData(WbxB, "WbxB", depth, iterN);
    }

    double _s = 0;
    B2sumW[0] = 0;
    for (int kz = 0; kz < N[2]; kz++)
    {
        _s += Lfunc[kz];
        B2sumW[0] += B2sum[kz];
    }

    mode = 1;
    main_proc->proceed(processors, supervisor, priority);

    //DebugWriteData(vF, "vF", depth, iterN);

    return _s;
}

//-----------------------------------------------------------------------
bool PCOProcessor::proceed()
{
    if (host->mode == 0)
    {
        B2->abs2_plane_lev(host->vB, task->z, 0);
        host->B2sum[task->z] = B2->sum_plane(0, host->sW);

        s->inv_plane(B2, 0); // s1 == 1 / B^2

        // Wa
        host->rotB->rot_plane(host->vB, task->z, WiegelmannDerivStencil); // rotB ~ J
        host->Wa->cross_plane(host->rotB, host->vB, task->z); // rotB x B
        // Wb
        host->divB->div_plane(host->vB, task->z, WiegelmannDerivStencil);
        host->Wb->mult_plane(host->divB, host->vB, task->z); // divB * B
        if (WiegelmannGetMetricsTheta)
        {
            host->Bmod->sqrt_plane_lev(B2, 0, task->z);
            host->Binv->inv_plane(host->Bmod, task->z);
            host->Jmod->abs_plane(host->rotB, task->z);
            host->Jinv->inv_plane(host->Jmod, task->z);
            host->JxB->abs_plane(host->Wa, task->z);

            host->st->mult_plane(host->JxB, host->Binv, task->z);
            host->st->mult_plane(host->st, host->Jinv, task->z);
            host->sS[task->z] = host->st->sum_plane(task->z, host->sW);
            host->sSW[task->z] = host->sW->sum_plane(task->z);

            host->st->mult_plane(host->JxB, host->Binv, task->z);
            host->sSJ[task->z] = host->st->sum_plane(task->z, host->sW);
            host->sJ[task->z] = host->Jmod->sum_plane(task->z, host->sW);

            host->st->mult_plane(host->JxB, host->Jinv, task->z);
            host->sSB[task->z] = host->st->sum_plane(task->z, host->sW);
            host->sB[task->z] = host->Bmod->sum_plane(task->z, host->sW);
        }

        host->Wa->mult_plane_lev(s, host->Wa, 0, task->z, task->z); // rotB x B / B^2
        host->Wb->mult_plane_lev(s, host->Wb, 0, task->z, task->z); // divB * B / B^2

        // Wa^2 + Wb^2
        host->Wa2Wb2->abs2_plane(host->Wa, task->z); // here is Wa^2
        if (WiegelmannGetMetricsTheta)
            host->LF[task->z] = host->Wa2Wb2->sum_plane(task->z, host->sW);

        s->abs2_plane_lev(host->Wb, task->z, 0);
        if (WiegelmannGetMetricsTheta)
            host->LD[task->z] = s->sum_plane_lev(0, host->sW, task->z);

        if (WiegelmannWeightDivfree != 1.0)
            s->mult_plane(WiegelmannWeightDivfree, s, 0);
        host->Wa2Wb2->add_plane_lev(host->Wa2Wb2, s, task->z, 0);

        // functional
        s->mult_plane_lev(B2, host->Wa2Wb2, 0, task->z);

        //if (task->z == 0)
        //    host->Lfunc[task->z] = 0;
        //else
            host->Lfunc[task->z] = s->sum_plane_lev(0, host->sW, task->z);

        host->WaxB->cross_plane(host->Wa, host->vB, task->z);
        host->WbxB->dot_plane(host->Wb, host->vB, task->z);
    }
    else
    {
        host->vF->mult_plane(host->Wa2Wb2, host->vB, task->z);
        v->rot_plane_lev(host->WaxB, task->z, 0, WiegelmannDerivStencil);
        host->vF->add_plane_lev(v, task->z, 0);

        v->cross_plane_lev(host->Wa, host->rotB, task->z, 0);
        host->vF->sub_plane_lev(v, task->z, 0);

        v->grad_plane_lev(host->WbxB, task->z, 0, WiegelmannDerivStencil);
        if (WiegelmannWeightDivfree != 1.0)
            v->mult_plane(WiegelmannWeightDivfree, v, 0);
        host->vF->add_plane_lev(v, task->z, 0); // rot(Wa x B) - Wa x rotB + grad(Wb . B) + (Wa^2 + Wb^2)*B

        v->mult_plane_lev(host->divB, host->Wb, task->z, 0);
        if (WiegelmannWeightDivfree != 1.0)
            v->mult_plane(WiegelmannWeightDivfree, v, task->z);
        host->vF->sub_plane_lev(v, task->z, 0); //  rot(Wa x B) - Wa x rotB + grad(Wb . B) + (Wa^2 + Wb^2)*B - Wb*divB

        host->vF->mult_plane(host->sW, host->vF, task->z);

        // gradW
        v->cross_plane_lev(host->WaxB, host->vgradW, task->z, 0);
        host->vF->add_plane_lev(v, task->z, 0);

        v->mult_plane_lev(host->WbxB, host->vgradW, task->z, 0);
        if (WiegelmannWeightDivfree != 1.0)
            v->mult_plane(WiegelmannWeightDivfree, v, 0);
        host->vF->add_plane_lev(v, task->z, 0); // rot(Wa x B) x gradW + (Wb . B)*gradW
        host->F2max[task->z] = host->vF->max2_plane(task->z);

        if (WiegelmannProcCondType == 1 && host->baseField && host->baseWeight && task->z < host->constrN)
        {
            v->sub_plane_lev(host->vB, host->baseField, task->z, task->z, 0);
            B2->dot_plane_lev(v, v, 0, 0, host->baseWeight, task->z);
            host->Lfunc[task->z] += B2->sum_plane(0);

            v->mult_plane_lev(host->baseWeight, v, task->z, 0, 0);
            host->vF->sub_plane_lev(v, task->z, 0);
        }
    }

    return true;
}

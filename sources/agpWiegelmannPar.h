#pragma once

#include "stdDefinitions.h"
#include "mfoGlobals.h"
#include "agmMetrics.h"
#include "TaskQueueProcessor.h"
#include "MagFieldOps.h"

#ifdef _WINDOWS
#include <windows.h>
#endif

class CagmScalarFieldOps;
class CagmVectorFieldOps;
class CagmScalarField;
class CagmVectorField;

//---------------------------------------------------------------------------------------
class PCOTask : public ATQPTask
{
public:
    int z;

public:
    PCOTask(int _id) : ATQPTask(_id) {}
    void setData(int _z)
    {
        z = _z;
    }

    virtual ~PCOTask()
    {
    }
};

//---------------------------------------------------------------------------------------
class PCOTaskFactory : public ATQPTaskFactory
{
public:
    PCOTaskFactory() : ATQPTaskFactory() {}
    virtual ~PCOTaskFactory() {}

    virtual ATQPTask *create()
    {
        return new PCOTask(counter++);
    }
};

//---------------------------------------------------------------------------------------
class PCOSupervisor : public ATQPSupervisor
{
public:
    PCOSupervisor(int *N, ATQPTaskFactory* _factory, ATQPSynchonizer *_sync)
    : ATQPSupervisor(N[2], _factory, _sync)
    {
        int cnt = 0;
        for (int z = 0; z < N[2]; z++)
            ((PCOTask *)tasks[cnt++])->setData(z);
    }

    virtual ~PCOSupervisor() {}
};

class CagpWiegelmann;

//---------------------------------------------------------------------------------------
class PCOProcessor : public ATQPProcessor
{
protected:
    PCOTask *task;
    CagpWiegelmann* host;
    int stencil;

    CagmVectorField *v;
	CagmScalarField *B2, *s;

public:
    PCOProcessor(int id, ATQPSynchonizer *_sync, CagpWiegelmann* _host, int *N, int _stencil) 
        : ATQPProcessor(id, _sync)
        , host(_host)
        , stencil(_stencil)
    {

        int Nplane[3];
        Nplane[0] = N[0];
        Nplane[1] = N[1];
        Nplane[2] = 1;

        v = new CagmVectorField(Nplane);
        B2 = new CagmScalarField(Nplane);
        s = new CagmScalarField(Nplane);
        s->setTolerances(WiegelmannInversionTolerance, WiegelmannInversionDenom);
    }
    virtual ~PCOProcessor()
    {
        delete v;
        delete B2;
        delete s;
    }

    virtual bool setTask(ATQPTask * _task);
    virtual bool proceed();
};

//------------------------------------------------------------------
class CagpWiegelmann : public CagmMetrics // CagpWiegelmannBase
{
friend class PCOProcessor;
public:
    double *F2max, *B2sum;

protected:

    // calculate
    CagmVectorField *vgradW;

    // intermediate
  	CagmVectorField *rotB, *Wa, *Wb, *WaxB;
	CagmScalarField *divB, *Wa2Wb2, *WbxB;

    int N[3];
    double dircos[3];

    // sources input
    CagmVectorField *vB;
    CagmScalarField *sW;
    CagmScalarField *bottomWeight;
    double *vdircos;

    // sources output
    CagmVectorFieldOps *vF;
    CagmVectorFieldOps *vBottom;

    // sources constrains
    CagmVectorField *baseField;
    CagmVectorField *baseWeight;
    CagmScalarField *absField;
    CagmScalarField *absWeight;
    CagmScalarField *losField;
    CagmScalarField *losWeight;

    int constrN;

    // metrics
    CagmScalarField *Bmod, *Jmod, *Binv, *Jinv, *JxB, *st;
    double *Lfunc;

    PCOSupervisor *supervisor;
    std::vector<ATQPProcessor *> processors;
    TaskQueueProcessor *main_proc;
    PCOTaskFactory *factory;
    int nProc;

    int mode;

    int iterN;
    int depth;

protected:
    w_priority priority;

public:
    CagpWiegelmann(int *_N, int _n_threads, int _stencil
        , CagmVectorField *_sourceB, CagmScalarField *_sourceW
        , double *_vcos
        , CagmVectorFieldOps *_outF
        , CagmVectorField *_baseField, CagmVectorField *_baseWeight
        , CagmScalarField *_absField, CagmScalarField *_absWeight
        , CagmScalarField *_losField, CagmScalarField *_losWeight
        , CagmScalarField *_bottomWeight
        , int _depth
        , w_priority _priority = w_priority::low);
	virtual ~CagpWiegelmann();

    double step(int _iterN);

    static uint64_t estimateMemory(int *N, int nThreads, bool full = true)
    {
        int np = TaskQueueProcessor::getProcInfo(nThreads);
        uint64_t P = (N[0]+1) * (N[1]+1);
        uint64_t S = N[1] * N[2];
        uint64_t V = P * (N[2]+1);

        uint64_t M = (V*(8*3 + 4) + P*5*np)*sizeof(double) + (6*3 + 4)*S*sizeof(double *);

        if (full)
            M += V * (3+1) * sizeof(double);

        return M;
    }
};

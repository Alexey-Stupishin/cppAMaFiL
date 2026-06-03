#pragma once

#include "TaskQueueProcessor.h"
#include "TimeTicToc.h"

#ifdef _WINDOWS
#pragma warning (disable : 4996)

#define path_disamb_step "c:\\temp\\disabm_step.txt"
#define path_disamb "c:\\temp\\disamb\\disabm_"
#endif

//---------------------------------------------------------------------------------------
class DAMTask : public ATQPTask
{
public:
    int N;
    int *x, *y;

public:
    DAMTask(int _id = 0) 
        : ATQPTask(_id) 
        , N(0)
        , x(nullptr)
        , y(nullptr)
    {}
    virtual ~DAMTask() {}

    void setData(int _N, int *_x, int *_y)
    {
        N = _N;
        x = _x;
        y = _y;
    }
};

//---------------------------------------------------------------------------------------
class DAMTaskFactory : public ATQPTaskFactory
{
public:
    DAMTaskFactory() : ATQPTaskFactory() {}
    virtual ~DAMTaskFactory() {}

    virtual ATQPTask *create()
    {
        return new DAMTask(counter++);
    }
};

//---------------------------------------------------------------------------------------
class DAMSupervisor : public ATQPSupervisor
{
public:
    DAMSupervisor(int _n_tasks, ATQPTaskFactory* _factory, ATQPSynchonizer *_sync)
    : ATQPSupervisor(_n_tasks, _factory, _sync)
    {
    }

    virtual ~DAMSupervisor() {}

    void setData(int n_tasks, int chunk_lng, int n_poi, int *task_x, int *task_y);
};

class CmfoDisambig;

//---------------------------------------------------------------------------------------
class DAMProcessor : public ATQPProcessor
{
protected:
    DAMTask *task;
    CmfoDisambig* host;

    std::mutex mutex_store;

public:
    DAMProcessor(int id, ATQPSynchonizer *_sync, CmfoDisambig* _host) 
        : ATQPProcessor(id, _sync)
        , host(_host)
    {
    }
    virtual ~DAMProcessor()
    {
    }

    virtual bool setTask(ATQPTask * _task)
    {
        bool finish = ATQPProcessor::setTask(_task);
        if (!finish)
            task = (DAMTask *)_task;

        return finish;
    }

    virtual bool proceed();
};

//---------------------------------------------------------------------------------------
class CmfoDisambigRules
{
public:
    int limatt, limacc, min_generation;
    double ordpart;
    double temp_min, temp_decr;
    double acc2att;
    double Bmax_term, non_stable_term;
};

//---------------------------------------------------------------------------------------
class CmfoDisambig : public TaskQueueProcessor
{
friend class DAMProcessor;

public:
    enum Status { none = 0,  T_min = 1,  acc_2_att = 2, bmax_term = 3, non_stable_term = 4};

public:
    int nacc;
    int natt;
    Status status;

protected:
    std::mutex mutex_store;

    int sNx, sNy, fNx, fNy, tNx, tNy;
    int sN[2], fN[2], tN[2];

    int *prescribe;
    int *state, *generation;
    double *F, *kT;
    int *state_prev;
    double *Fstate;
    int non_prescr;
    int *nprescrKx, *nprescrKy;
    double *absB;
    int step_n;

    bool cleanup;
    CmfoDisambigRules *rules;
    bool parallel;
    int n_proc;

    int n_threads_init;

    aguTimeTicToc tic;

#ifdef _WINDOWS
#ifdef _DEBUG_DISAMB_STEP
    FILE *fid;
#endif
#endif

public:
    CmfoDisambig(int *N, int *_state, int *_generation, double *absB, double *_F, double *_kT, bool _parallel = true, int _n_threads_init = 0);
    virtual ~CmfoDisambig();
    
    int init(CmfoDisambigRules *rules, double *_Fstate, int *_prescribe = nullptr);
    int create(CmfoDisambigRules *rules, int *_prescribed = nullptr);
    int temp_step();
    int process();

protected:
    static void Durstenfeld(int n, int *x, int *y);
    int get_F_state();
    int step_proceed();
    int get_order();
    int temp_step_task_core(DAMTask *task);
};

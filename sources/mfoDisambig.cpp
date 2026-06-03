#include "stdDefinitions.h"
#include <random>
#include <chrono>

#include "mfoDisambig.h"

#ifdef _WINDOWS
#include "binUtilitiesW.h"
#endif

CmfoDisambig::CmfoDisambig(int *N, int *_state, int *_generation, double *_absB, double *_F, double *_kT, bool _parallel, int _n_threads_init)
        : TaskQueueProcessor(_n_threads_init)
        , parallel(_parallel)
        , cleanup(false)
        , status(CmfoDisambig::none)
        , absB(_absB)
        , state(_state)
        , generation(_generation)
        , F(_F)
        , kT(_kT)
{
    n_threads_init = _n_threads_init;
#ifdef _WINDOWS
#ifdef _DEBUG_DISAMB_STEP
    fid = fopen(path_disamb_step, "ab");
#endif
#endif
    tN[0] = N[0];
    tN[1] = N[1];
    tNx = tN[0];
    tNy = tN[1];
    fN[0] = tNx + 1;
    fN[1] = tNy + 1;
    fNx = fN[0];
    fNy = fN[1];
    sN[0] = fNx + 1;
    sN[1] = fNy + 1;
    sNx = sN[0];
    sNy = sN[1];
}

//-----------------------------------------------------------------------
CmfoDisambig::~CmfoDisambig()
{
    delete [] nprescrKx;
    delete [] nprescrKy;
    if (cleanup)
    {
        delete [] Fstate;
        delete [] state_prev;
        delete [] prescribe;
    }

#ifdef _WINDOWS
#ifdef _DEBUG_DISAMB_STEP
fclose(fid);
#endif
#endif
}

//-----------------------------------------------------------------------
int CmfoDisambig::init(CmfoDisambigRules *_rules, double *_Fstate, int *_prescribe)
{
    Fstate = _Fstate;
    prescribe = _prescribe;
    rules = _rules;

    nprescrKx = new int[sNx*sNy];
    nprescrKy = new int[sNx*sNy];

    return 0;
}

//-----------------------------------------------------------------------
int CmfoDisambig::create(CmfoDisambigRules *_rules, int *_prescribed)
{
    rules = _rules;

    cleanup = true;

    Fstate = new double[fNx*fNy];
    state_prev = new int[sNx*sNy];
    prescribe = new int[sNx*sNy];

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distrib(0.0, 2.0);

    memset(state, 0, sNx*sNy*sizeof(int));
    memset(prescribe, 0, sNx*sNy*sizeof(int));
    memset(generation, 0, sNx*sNy*sizeof(int));
    for (int ky = 1; ky < sNy-1; ky++)
        for (int kx = 1; kx < sNx - 1; kx++)
        {
            int p = _prescribed[(kx-1)+(ky-1)*tNx];
            state[kx + ky*sNx] = p > 0 ? 0 : (p < 0 ? 1 : (distrib(gen) < 1 ? 0 : 1));
            prescribe[kx + ky*sNx] = (p != 0);
        }

    nprescrKx = new int[sNx*sNy];
    nprescrKy = new int[sNx*sNy];

    return 0;
}

//-----------------------------------------------------------------------
bool DAMProcessor::proceed()
{
    int rc = host->temp_step_task_core(task);
    return rc == 0;
}

//-----------------------------------------------------------------------
void DAMSupervisor::setData(int n_tasks, int chunk_lng, int n_poi, int *task_x, int *task_y)
{
    int rest = n_poi;
    int pos = 0;
    for (int z = 0; z < n_tasks; z++)
    {   
        if (rest > 0)
        {
            int lng = chunk_lng;
            if (rest < chunk_lng)
                lng = rest;
            ((DAMTask *)tasks[z])->setData(lng, task_x + pos, task_y + pos);
            rest -= lng;
            pos += lng;
        }
    }
}

//-----------------------------------------------------------------------
void CmfoDisambig::Durstenfeld(int n, int *x, int *y)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int k = n-1; k > 0; k--)
    {
        std::uniform_int_distribution<int> distrib(0, k);
        int idx = distrib(gen);
        int t = x[k]; x[k] = x[idx]; x[idx] = t;
            t = y[k]; y[k] = y[idx]; y[idx] = t;
    }
}

//-----------------------------------------------------------------------
int CmfoDisambig::get_order()
{
    int nprescrL = 0;

    for (int ky = 1; ky < sNy-1; ky++)
        for (int kx = 1; kx < sNx-1; kx++)
        {
            if (prescribe[ky*sNx+kx] == 0)
            {
                nprescrKx[nprescrL] = kx;
                nprescrKy[nprescrL] = ky;
                nprescrL++;
            }
        }

    return nprescrL;
}

//-----------------------------------------------------------------------
int CmfoDisambig::temp_step()
{
#ifdef _WINDOWS
#ifdef _DEBUG_DISAMB_STEP
    static FILETIME FileTime0, FileTime1;
    ULARGE_INTEGER time0, time1;
    GetSystemTimeAsFileTime(&FileTime0);
    time0.HighPart = FileTime0.dwHighDateTime;
    time0.LowPart = FileTime0.dwLowDateTime;

    tic.tic();
#endif
#endif

    non_prescr = get_order();
    int n_poi = (int)floor(non_prescr * rules->ordpart);

    std::vector<ATQPProcessor *> processors;

    DAMSupervisor *supervisor = nullptr;
    DAMTaskFactory this_factory;
    int chunk_lng = n_poi;
    if (parallel)
    {
        n_proc = get_num_proc();
        chunk_lng = (int)ceil((double)n_poi / n_proc);
        supervisor = new DAMSupervisor(n_proc, &this_factory, get_sync());
        for (int i = 0; i < n_proc; i++)
            processors.push_back(new DAMProcessor(i, get_sync(), this));
    }

    nacc = 0;
    natt = 0;

    int try_cnt = 0;
    while (natt <= ceil(rules->limatt*n_poi) && nacc <= ceil(rules->limacc*n_poi))
    {
        Durstenfeld(non_prescr, nprescrKx, nprescrKy);

        if (parallel)
        {
            supervisor->setData(n_proc, chunk_lng, n_poi, nprescrKx, nprescrKy);
            proceed(processors, supervisor, w_priority::low);
        }
        else
        {
            DAMTask task;
            task.setData(n_poi, nprescrKx, nprescrKy);
            temp_step_task_core(&task);
        }

        try_cnt++;
    }

    delete supervisor;

#ifdef _WINDOWS
#ifdef _DEBUG_DISAMB_STEP
    GetSystemTimeAsFileTime(&FileTime1);
    time1.HighPart = FileTime1.dwHighDateTime;
    time1.LowPart = FileTime1.dwLowDateTime;
    ULONGLONG udt = time1.QuadPart-time0.QuadPart;
    double utt = (double)udt*1e-7;
    fprintf(fid, "%9.5e %9.5e\r\n", utt, tic.toc());
#endif
#endif

    return try_cnt;
}

//-----------------------------------------------------------------------
int CmfoDisambig::temp_step_task_core(DAMTask *task)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distrib(0.0, 1.0);

    int task_nacc = 0;
    int task_natt = 0;

    for (int k = 0; k < task->N; k++)
    {
        int kx = task->x[k] - 1;
        int ky = task->y[k] - 1;

        int sbase, tbase, fbase, base4s, base4f, s16, fNyfNx = fNy*fNx;
        double Fnew[4], Fold[4], FnewS, FoldS;
        bool accepted;

        tbase = ky*tNx + kx;
        fbase = (ky + 1)*fNx + kx + 1;
        sbase = (ky + 1)*sNx + kx + 1;

        state[sbase] = 1 - state[sbase];

        base4s = sbase;
        base4f = fbase - fNx - 1;
        s16 = state[base4s - sNx] | ((state[base4s - sNx - 1]) << 1) | ((state[base4s]) << 2) | ((state[base4s - 1]) << 3);
        Fold[0] = Fstate[base4f];
        Fnew[0] = F[s16*fNyfNx + base4f];

        base4s = sbase + sNx;
        base4f = fbase - 1;
        s16 = state[base4s - sNx] | ((state[base4s - sNx - 1]) << 1) | ((state[base4s]) << 2) | ((state[base4s - 1]) << 3);
        Fold[1] = Fstate[base4f];
        Fnew[1] = F[s16*fNyfNx + base4f];

        base4s = sbase + 1;
        base4f = fbase - fNx;
        s16 = state[base4s - sNx] | ((state[base4s - sNx - 1]) << 1) | ((state[base4s]) << 2) | ((state[base4s - 1]) << 3);
        Fold[2] = Fstate[base4f];
        Fnew[2] = F[s16*fNyfNx + base4f];

        base4s = sbase + sNx + 1;
        base4f = fbase;
        s16 = state[base4s - sNx] | ((state[base4s - sNx - 1]) << 1) | ((state[base4s]) << 2) | ((state[base4s - 1]) << 3);
        Fold[3] = Fstate[base4f];
        Fnew[3] = F[s16*fNyfNx + base4f];

        accepted = false;
        FnewS = Fnew[0] + Fnew[1] + Fnew[2] + Fnew[3];
        FoldS = Fold[0] + Fold[1] + Fold[2] + Fold[3];
        if (FnewS < FoldS || distrib(gen) < exp(-(FnewS - FoldS) / kT[tbase]))
            accepted = true;

        if (accepted)
        {
            task_nacc++;

            if (parallel)
                std::unique_lock<std::mutex> locker(mutex_store);

            Fstate[fbase - fNx - 1] = Fnew[0];
            Fstate[fbase - 1] = Fnew[1];
            Fstate[fbase - fNx] = Fnew[2];
            Fstate[fbase] = Fnew[3];
        }
        else
            state[sbase] = 1 - state[sbase];
        task_natt++;
    }
    {
        if (parallel)
            std::unique_lock<std::mutex> locker(mutex_store);

        nacc += task_nacc;
        natt += task_natt;
    }

    return 0;
}

//-----------------------------------------------------------------------
int CmfoDisambig::get_F_state()
{
    int sizeF = fNx*fNy;
    for (int ky = 0; ky < fNy; ky++)
        for (int kx = 0; kx < fNx; kx++)
        {
            int idx = state[kx + 1 + ky*sNx] + ((state[kx + ky*sNx]) << 1) + ((state[kx + 1 + (ky + 1)*sNx]) << 2) + ((state[kx + (ky + 1)*sNx]) << 3);
            Fstate[kx + ky*fNx] = F[kx + ky*fNx + sizeF*idx];
        }

    return 0;
}

//-----------------------------------------------------------------------
int CmfoDisambig::step_proceed()
{
    int n_maxgen = 0;
    for (int ky = 1; ky < sNy-1; ky++)
        for (int kx = 1; kx < sNx-1; kx++)
        {
            if (generation[kx + ky*sNx] == step_n)
                n_maxgen++;
        }

    for (int ky = 0; ky < sNy; ky++)
        for (int kx = 0; kx < sNx; kx++)
        {
            int p = kx + ky*sNx;
            if (state[p] == state_prev[p])
                generation[p] += 1;
            else
                generation[p] = 0;

            if (n_maxgen == 0 && generation[p] >= rules->min_generation)
                prescribe[p] = 1;
        }

    return 0;
}

//-----------------------------------------------------------------------
int CmfoDisambig::process()
{
    double Temp = 1;

    step_n = 0;
    while (true)
    {
#ifdef _WINDOWS
#ifdef _DEBUG_DISAMB
    char fn[256], buf[16];
    strcpy(fn, path_disamb);
    strcat(fn, itoa(step_n, buf, 10));
    strcat(fn, ".bin");
    FILE *fid_da = fopen(fn, "wb");
    CbinDataStruct::WriteHeader(fid_da);
    CbinDataStruct::Write(fid_da, &tNx, 1, "N0");
    CbinDataStruct::Write(fid_da, &tNy, 1, "N1");
    CbinDataStruct::Write(fid_da, state, sNx*sNy, "state");
    CbinDataStruct::Write(fid_da, prescribe, sNx*sNy, "prescribe");
    CbinDataStruct::Write(fid_da, generation, sNx*sNy, "generation");
    CbinDataStruct::WriteFooter(fid_da);
    fclose(fid_da);
#endif
#endif
        memcpy(state_prev, state, sNx*sNy*sizeof(int));
        get_F_state();
        temp_step();

        step_proceed();
        step_n++;

        if (Temp < rules->temp_min)
        {
            status = T_min;
            break;
        }

        if ((double)nacc / natt <= rules->acc2att)
        {
            status = acc_2_att;
            break;
        }

        double B_np_max = 0;
        for (int ky = 1; ky < sNy-1; ky++)
            for (int kx = 1; kx < sNx-1; kx++)
            {
                if (prescribe[kx + ky*sNx] == 0)
                {
                    double v = absB[kx-1 + (ky-1)*tNx];
                    B_np_max = v > B_np_max ? v : B_np_max;
                }
            }
        if (B_np_max <= rules->Bmax_term)
        {
            status = bmax_term;
            break;
        }

        if (non_prescr <= rules->non_stable_term*sNx*sNy)
        {
            status = non_stable_term;
            break;
        }

        for (int ky = 0; ky < tNy; ky++)
            for (int kx = 0; kx < tNx; kx++)
                kT[kx + ky*tNx] *= rules->temp_decr;

        Temp *= rules->temp_decr;
    }

    return step_n;
}

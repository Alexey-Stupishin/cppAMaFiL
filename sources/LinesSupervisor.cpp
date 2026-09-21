#include "stdDefinitions.h"

#include "LinesSupervisor.h"
#include "Durstenfeld.h"

LQPSupervisor::LQPSupervisor(CagmVectorField *v,
        uint32_t _cond, double _chromoLevel,
        double *_seeds, int _Nseeds, double relSeedsBound, int lines_use_durstenfeld,
        int *_voxelStatus, double *_physLength, double *_avField,
        int *_linesLength, int *_codes, double *_times,
        int *_startIdx, int *_endIdx, int *_apexIdx,
        uint64_t _maxCoordLength, double *_coords, uint64_t *_linesStart, int *_linesIndex, int *_seedIdx,
        LQPTaskFactory *factory, ATQPSynchonizer *_sync)
        : ATQPSupervisor(0, factory, _sync)
        , n(0)
        , field(v)
        , cond(_cond)
        , chromoLevel(_chromoLevel)
        , maxCoordLength(4*_maxCoordLength)
        , linesLength(_linesLength)
        , codes(_codes)
        , times(_times)
        , voxelStatus(_voxelStatus)
        , physLength(_physLength)
        , avField(_avField)
        , startIdx(_startIdx)
        , endIdx(_endIdx)
        , apexIdx(_apexIdx)
        , coords(_coords)
        , linesStart(_linesStart)
        , linesIndex(_linesIndex)
        , seedIdx(_seedIdx)
        , cumLength(0)
        , coordTol(1e-3)
        , closedTol(3e-3)
        , nPassed(0)
        , passed(nullptr)
        , globalID(nullptr)
        , nLines(0)
        , nNonStored(0)
        , seeds(_seeds)
        , params(nullptr)
{
    cube = (CubeXD *)field;

    field->dimensions(NF);
    autoParams = (_Nseeds <= 0 || !seeds);

    int from[3], to[3];
    if (!autoParams)
    {
        Nseeds = _Nseeds;
    }
    else
    {
        from[0] = 0;
        from[1] = 0;
        from[2] = 0;
        to[0] = NF[0]-1;
        to[1] = NF[1]-1;
        to[2] = NF[2]-1;
        if (relSeedsBound > 0)
        {
            from[0] = (int)ceil(relSeedsBound*NF[0]);
            to[0] = NF[0] - 1 - from[0];
            from[1] = (int)ceil(relSeedsBound*NF[1]);
            to[1] = NF[1] - 1 - from[1];
            to[2] = NF[2] - 1 - (int)(ceil(relSeedsBound*NF[2]));
            Nseeds = _Nseeds;
        }
        Nseeds = (to[0]-from[0]+1)*(to[1]-from[1]+1)*(to[2]-from[2]+1);
        if (cond != Conditions::NoCond)
            passed = new int[NF[0]*NF[1]*NF[2]];
    }

    InitQueue(Nseeds);

    globalID = new int[Nseeds];
    if (!autoParams)
    {
        params = new double[3 * Nseeds];
        for (int id = 0; id < Nseeds; id++)
        {
            InitOutput(id);
            params[3*id + 0] = seeds[3*id + 0];
            params[3*id + 1] = seeds[3*id + 1];
            params[3*id + 2] = seeds[3*id + 2];
            globalID[id] = cube->getGlobalID(seeds + 3*id);
        }
    }
    else
    {
        params = new double[3];

        int id = 0;
        for (int kz = from[2]; kz <= to[2]; kz++)
            for (int ky = from[1]; ky <= to[1]; ky++)
                for (int kx = from[0]; kx <= to[0]; kx++)
                {
                    globalID[id] = cube->getGlobalID(kx, ky, kz);
                    InitOutput(globalID[id]);
                    id++;
                }
        if (lines_use_durstenfeld)
            Durstenfeld(Nseeds, globalID);
    }
}

LQPSupervisor::~LQPSupervisor()
{
    delete[] params;
    delete[] globalID;
    delete[] passed;
}

//-------------------------------------------------------------------------------
bool LQPSupervisor::getTask(ATQPTask*& task)
{
    task = nullptr;
    long n_item;
    while (true)
    {
        n_item = NextQueueItem();
        if (n_item < 0)
            return false;
        //double *point = (double *)GetParams(n_item);
        GetParams(n_item);
        if (passed ? !passed[n_item] : true)
            break;
    }
    task = factory->create();
    ((LQPTask *)task)->setData((double *)GetParams(n_item));
    ((LQPTask *)task)->setID(n_item);
    return true;
}

//-------------------------------------------------------------------------------
long LQPSupervisor::NextQueueItem()
{
    if (n >= nQueue)
        return -1; //  ABaseTaskQueue::NoItems;

    long nret = autoParams ? globalID[n] : n;
    n++;
    return nret;
}

//-------------------------------------------------------------------------------
void * LQPSupervisor::GetParams(uint32_t _queueID)
{
    if (!autoParams)
        return params + _queueID * 3;
    else
    {
        int kx, ky, kz;
        cube->parseGlobalID(_queueID, &kx, &ky, &kz);
        params[0] = kx;
        params[1] = ky;
        params[2] = kz;
        return params;
    }
}

//-------------------------------------------------------------------------------
void LQPSupervisor::InitOutput(int id)
{
    if (physLength)
        physLength[id] = 0;
    if (avField)
        avField[id] = 0;
    if (startIdx)
        startIdx[id] = 0;
    if (endIdx)
        endIdx[id] = 0;
    if (voxelStatus)
        voxelStatus[id] = Status::None;
    if (passed)
        passed[id] = false;
}

//-------------------------------------------------------------------------------
uint32_t LQPSupervisor::InitQueue(int _nQueue)
{
    return nQueue = _nQueue;
}

//-------------------------------------------------------------------------------
uint32_t LQPSupervisor::SetResult(int /* proc_id */, uint32_t queueID, LQPLineResult *line, int _code4over, double time)
{
    std::unique_lock<std::mutex> locker(sync->mutex_query);
    return proceedLine(queueID, line, _code4over, time);
}

//-------------------------------------------------------------------------------
void LQPSupervisor::getDivPoint(double *p1, double *p2)
{
    double t = (chromoLevel - p1[2])/(p2[2] - p1[2]);
    p1[0] = p1[0] + t*(p2[0] - p1[0]);
    p1[1] = p1[1] + t*(p2[1] - p1[1]);
    p1[2] = chromoLevel;
}

//-------------------------------------------------------------------------------
uint32_t LQPSupervisor::proceedLine(uint32_t queueID, LQPLineResult *line, int _code4over, double time)
{
    if (voxelStatus)
        voxelStatus[queueID] = Status::Processed;
    if (apexIdx)
        apexIdx[queueID] = 0;
    if (seedIdx)
        seedIdx[queueID] = 0;
    if (startIdx)
        startIdx[queueID] = 0;
    if (endIdx)
        endIdx[queueID] = 0;
    if (physLength)
        physLength[queueID] = 0;
    if (avField)
        avField[queueID] = 0;
    if (codes)
        codes[queueID] = 0;
    if (times)
        times[queueID] = time;

    if (line->useful)
    {
        bool mark_passed = ((line->status & Status::Closed) && (cond & Conditions::PassClosed)) || (!(line->status & Status::Closed) && (cond & Conditions::PassOpen));
        for (int k = line->start; k <= line->end; k++)
            proceedVox(line->indices[k], line->status, line->phys_length, line->av_field, line->apex_idx, line->seed_idx, line->start_idx, line->end_idx, mark_passed);
        voxelStatus[queueID] |= Status::BaseVoxel;

        nPassed++;
        if (coords)
        {
            int length = line->end - line->start + 1;
            if (4 * (cumLength + length) > maxCoordLength)
            {
                line->code |= _code4over;
                nNonStored++;
            }
            else
            {
                for (int k = line->start; k <= line->end; k++)
                {
                    uint64_t gpos = 4 * (cumLength + k - line->start);
                    memcpy(coords + gpos, line->coords + 3*k, 3*sizeof(double));
                    coords[gpos + 3] = line->distance[k];
                }

                if (linesStart)
                    linesStart[nLines] = cumLength;
                if (linesLength)
                    linesLength[nLines] = length;
                if (linesIndex)
                    linesIndex[nLines] = queueID;
                nLines++;
                cumLength += length;
            }
        }

        if (codes)
            codes[queueID] = line->code;

    }

    return 0;
}

//-------------------------------------------------------------------------------
uint32_t LQPSupervisor::proceedThisVox(int queueID, int sClosed, double thisPhysLength, double thisAvField, int apexid, int seedid, int startid, int endid, bool mark_passed)
{
    if (voxelStatus)
        voxelStatus[queueID] |= Status::Processed | Status::Lined | sClosed;

    if (physLength)
        physLength[queueID] = thisPhysLength;
    if (avField)
        avField[queueID] = thisAvField;

    if (apexIdx)
        apexIdx[queueID] = apexid;
    if (seedIdx)
        seedIdx[queueID] = seedid;
    if (startIdx)
        startIdx[queueID] = startid;
    if (endIdx)
        endIdx[queueID] = endid;

    if (passed)
        passed[queueID] = mark_passed;

    return 0;
}

//-------------------------------------------------------------------------------
uint32_t LQPSupervisor::proceedVox(int queueID, int sClosed, double thisPhysLength, double thisAvField, int apexid, int seedid, int startid, int endid, bool mark_passed)
{
    if (autoParams)
        return proceedThisVox(queueID, sClosed, thisPhysLength, thisAvField, apexid, seedid, startid, endid, mark_passed);

    for (int k = 0; k < Nseeds; k++)
    {
        if (globalID[k] == queueID)
            proceedThisVox(k, sClosed, thisPhysLength, thisAvField, apexid, seedid, startid, endid, mark_passed);
    }

    return 0;
}

//-------------------------------------------------------------------------------
void LQPSupervisor::getFinalState(uint64_t *_cumLength, int *_nLines, int *_nPassed, int *_nNonStored)
{
    *_cumLength = cumLength;
    *_nLines = nLines;
    *_nPassed = nPassed;
    *_nNonStored = nNonStored;
}

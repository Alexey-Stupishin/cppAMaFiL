#pragma once

#include <thread>
#include <mutex>

#include "TaskQueueProcessor.h"
#include "NLFFFLinesTaskQueue.h"

//---------------------------------------------------------------------------------------
class LQPTask : public ATQPTask
{
protected:
    double data[3];
    unsigned long queueID;

public:
    LQPTask(int _id) : ATQPTask(_id) {}
    void setData(double *_data) { memcpy(data, _data, 3 * sizeof(double)); }
    double *getData() { return data; }
    void setID(unsigned long _id) { queueID = _id; }
    unsigned long getID() { return queueID; }
};

//---------------------------------------------------------------------------------------
class LQPTaskFactory : public ATQPTaskFactory
{
public:
    LQPTaskFactory() : ATQPTaskFactory() {}
    virtual ~LQPTaskFactory() {}

    virtual ATQPTask *create()
    {
        return new LQPTask(counter++);
    }
};

//---------------------------------------------------------------------------------------
class LQPSupervisor : public ATQPSupervisor
{
public:
    CNLFFFLinesTaskQueue *queue;

protected:
    std::mutex mutex_store;
    
public:
    LQPSupervisor(CagmVectorField *v,
        uint32_t _cond, double chromoLevel,
        double *_seeds, int _Nseeds, double relSeedsBound,
        int *_voxelStatus, double *_physLength, double *_avField,
        int *_linesLength, int *_codes, double *_times,
        int *_startIdx, int *_endIdx, int *_apexIdx,
        uint64_t _maxCoordLength, double *_coords, uint64_t *_linesStart, int *_linesIndex, int *seedIdx,
        LQPTaskFactory *factory, ATQPSynchonizer *_sync)
        : ATQPSupervisor(0, factory, _sync)
    {
        int maxResult = 50000;

        queue = new CNLFFFLinesTaskQueue(v, _seeds, _Nseeds, relSeedsBound, _cond, chromoLevel,
            _physLength, _avField,
            _voxelStatus, _codes, _times,
            _startIdx, _endIdx, _apexIdx,
            maxResult,
            4 * _maxCoordLength, _linesLength, _coords, _linesStart, _linesIndex, seedIdx);

        Reset();
    }

    void Reset()
    {
    }

    virtual ~LQPSupervisor()
    {
        delete queue;
    }

    virtual bool getTask(ATQPTask*& task)
    {
        task = nullptr;
        long n = queue->NextQueueItem();
        if (n < 0)
            return false;
        task = factory->create();
        ((LQPTask *)task)->setData((double *)queue->GetParams(n));
        ((LQPTask *)task)->setID(n);
        return true;
    }

    bool needProcessing(uint32_t queueID)
    {
        return queue->needProcessing(queueID);
    }

    uint32_t SetResult(uint32_t queueID, double *point, double *result, int resLength, int _code, int _code4over, double time)
    {
        std::unique_lock<std::mutex> locker(mutex_store);
        uint32_t res = queue->SetResult(queueID, point, result, resLength, _code, _code4over, time);

        return res;
    }
};

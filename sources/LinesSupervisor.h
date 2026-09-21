#pragma once

#include <thread>
#include <mutex>

#include "TaskQueueProcessor.h"
#include "LinesTask.h"

#include "agmVectorField.h"
#include "LinesSingleLine.h"

#define d_sdist(c1, c2) sqrt((0[c2]-0[c1])*(0[c2]-0[c1]) + (1[c2]-1[c1])*(1[c2]-1[c1]) + (2[c2]-2[c1])*(2[c2]-2[c1]))
#define d_snorm(c1) sqrt(0[c1]*0[c1] + 1[c1]*1[c1] + 2[c1]*2[c1])

#define FS_VOID 0
#define FS_IN   1
#define FS_OUT  2
#define FS_END  3

//---------------------------------------------------------------------------------------
class LQPSupervisor : public ATQPSupervisor
{
public:
    enum Status { None = 0, Processed = 1, Lined = 2, Closed = 4, BaseVoxel = 8, OnlyFootpoint = 16 };
    enum Conditions { NoCond = 0, PassClosed = 1, PassOpen = 2 };

private:
    int nQueue, n;

    CagmVectorFieldOps *field;
    CubeXD *cube;
    int NF[3];

    double chromoLevel;
    uint32_t cond;

    uint64_t maxCoordLength;

    double coordTol;
    double closedTol;

    // Caller allocated!
    // length = Nseeds, 1-D
    double *physLength, *avField;
    int *startIdx, *endIdx, *apexIdx;
    int *codes;
    double *times;
    int *voxelStatus;

    // length = maxCoordLength x 3, 2-D
    double *coords;
    // length = Nseeds, 1-D
    int *linesLength;
    uint64_t *linesStart;
    int *linesIndex;
    int *seedIdx;

    int *passed;

    int Nseeds;
    bool autoParams;
    int *globalID;

    uint64_t cumLength;
    int nLines, nPassed;
    int nNonStored;

    double *seeds;
    double *params;
    
public:
    LQPSupervisor(CagmVectorField *v,
        uint32_t _cond, double chromoLevel,
        double *_seeds, int _Nseeds, double relSeedsBound, int lines_use_durstenfeld,
        int *_voxelStatus, double *_physLength, double *_avField,
        int *_linesLength, int *_codes, double *_times,
        int *_startIdx, int *_endIdx, int *_apexIdx,
        uint64_t _maxCoordLength, double *_coords, uint64_t *_linesStart, int *_linesIndex, int *seedIdx,
        LQPTaskFactory *factory, ATQPSynchonizer *_sync);

    virtual ~LQPSupervisor();

    uint32_t SetResult(int proc_id, uint32_t queueID, LQPLineResult *line, int _code4over, double time);
    void getFinalState(uint64_t *_cumLength, int *_nLines, int *_nPassed, int *_nNonStored);

//-----------------------------------------------------------------
private:
    virtual bool getTask(ATQPTask*& task);
    void InitOutput(int id);

    virtual void * GetParams(uint32_t _queueID);
    virtual uint32_t InitQueue(int _nQueue);
    virtual long NextQueueItem();

    void getDivPoint(double *p1, double *p2);

    uint32_t LQPSupervisor::proceedLine(uint32_t queueID, LQPLineResult *line, int _code4over, double time);

    uint32_t proceedThisVox(int queueID, int sClosed, double thisPhysLength, double thisAvField, int apexid, int seedid, int startid, int endid, bool mark_passed);
    uint32_t proceedVox(int queueID, int sClosed, double thisPhysLength, double thisAvField, int apexid, int seedid, int startid, int endid, bool mark_passed);
};

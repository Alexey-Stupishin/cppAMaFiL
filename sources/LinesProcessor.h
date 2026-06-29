#pragma once

#include "agmRKF45.h"
#include "agmVectorFieldLineFuncs.h"
#include "LinesTask.h"
#include "LinesSupervisor.h"
#include "LinesSingleLine.h"
#include "TimeTicToc.h"

#define d_sdist(c1, c2) sqrt((0[c2]-0[c1])*(0[c2]-0[c1]) + (1[c2]-1[c1])*(1[c2]-1[c1]) + (2[c2]-2[c1])*(2[c2]-2[c1]))
#define d_snorm(c1) sqrt(0[c1]*0[c1] + 1[c1]*1[c1] + 2[c1]*2[c1])

#define FS_VOID 0
#define FS_IN   1
#define FS_OUT  2
#define FS_END  3

class LQPProcessor : public ATQPProcessor
{
public:
    enum Status { None = 0, Processed = 1, Lined = 2, Closed = 4, BaseVoxel = 8, OnlyFootpoint = 16 };

private:
    LQPTask *this_task;
    uint32_t queueID;
    CagmRKF45 *rkf45;
    CagmVectorFieldOps *v;
    CubeXD *cube;
    int NF[3];

    int dir;
    double step;
    double absBoundAchieve, relBoundAchieve;
    uint32_t cond;
    int maxLength;
    int *passed;
    double *coord;
    double *linesteps;
    int *indices;

    double point[3];

    LQPSupervisor *supervisor;

    double chromoLevel;
    double coordTol;
    double closedTol;

    double *distance;

    LQPLineResult line;

public:
    LQPProcessor(LQPSupervisor *_supervisor, int _id, CagmVectorFieldOps *_v, uint32_t _cond, double _chromoLevel, int _dir, double _step, double _relErr, double _absErr
        , double _absBoundAchieve, double _relBoundAchieve, int _maxLength, int *_passed, ATQPSynchonizer *_sync, int _n_loop_control, double loop_abs_cell);

    virtual ~LQPProcessor();

private:
    bool setTask(ATQPTask * _task);
    bool proceed();

    void getDivPoint(double *p1, double *p2);
    uint32_t proceedLine(int resLength, int _code);
};

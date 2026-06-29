#include "stdDefinitions.h"

#include "mfoGlobals.h"
#include "LinesProcessor.h"

bool fcond(void *p, const CagmRKF45Vect& v)
{
    return ((T_Lines *)p)->inBoundCube(v);
}

uint32_t fdata(void *p, const double /*t*/, const CagmRKF45Vect& v, CagmRKF45Vect& vp)
{
    return ((T_Lines *)p)->derivatives(v, vp);
}

//-------------------------------------------------------------
LQPProcessor::LQPProcessor(LQPSupervisor *_supervisor, int _id, CagmVectorFieldOps *_v, uint32_t _cond, double _chromoLevel, int _dir, double _step, double _relErr, double _absErr
          , double _absBoundAchieve, double _relBoundAchieve, int _maxLength, int *_passed, ATQPSynchonizer *_sync, int _n_loop_control, double loop_abs_cell)
    : ATQPProcessor(_id, _sync)
    ,supervisor(_supervisor)
    ,v(_v)
    ,dir(_dir)
    ,step(_step)
    ,absBoundAchieve(_absBoundAchieve) 
    ,relBoundAchieve(_relBoundAchieve)
    ,maxLength(_maxLength)
    ,passed(_passed)
    , chromoLevel(_chromoLevel)
    , coordTol(1e-3)
    , closedTol(3e-3)
    , cond(_cond)
{
    cube = (CubeXD *)v;
    v->dimensions(NF);
    rkf45 = new CagmRKF45(_absErr, _relErr, (RKF45_FUNCTION_VECTOR)fdata, 3, nullptr, (RKF45_FUNCTION_VECTOR_COND)fcond, absBoundAchieve, _n_loop_control, loop_abs_cell);
    coord = new double[3*maxLength];
    indices = new int[maxLength];
    distance = new double[maxLength];
}

//-------------------------------------------------------------
LQPProcessor::~LQPProcessor()
{
    delete rkf45;
    delete[] coord;
    delete[] distance;
}

//-------------------------------------------------------------
bool LQPProcessor::setTask(ATQPTask * _task)
{
    bool finish = ATQPProcessor::setTask(_task);
    if (!finish)
        this_task = (LQPTask *)_task;

    return finish;
}

//-------------------------------------------------------------
bool LQPProcessor::proceed()
{
    double *params = this_task->getData();
    memcpy(point, (double *)params, 3*sizeof(double)); 
    queueID = this_task->getID();

    int lineLength;
    int code;

    aguTimeTicToc tic;
    v->getOneFullLine(rkf45, point, dir, step, absBoundAchieve, relBoundAchieve, maxLength, &lineLength, coord, &code);

    proceedLine(lineLength, code);
    supervisor->SetResult(id, queueID, &line, CagmVectorFieldOps::Status::BufferOverload, tic.toc());

    delete this_task;

    return true;
}

//-------------------------------------------------------------------------------
void LQPProcessor::getDivPoint(double *p1, double *p2)
{
    double t = (chromoLevel - p1[2])/(p2[2] - p1[2]);
    p1[0] = p1[0] + t*(p2[0] - p1[0]);
    p1[1] = p1[1] + t*(p2[1] - p1[1]);
    p1[2] = chromoLevel;
}

//-------------------------------------------------------------------------------
uint32_t LQPProcessor::proceedLine(int resLength, int _code)
{
    line.init(coord, indices, distance);
    line.status = Status::Processed;

    if (resLength > 1) // && point[2] >= floor(chromoLevel))
    {
        line.av_field = 0;
        int state = FS_VOID;
        // find fragment
        bool found = false;
        int kPoint = -1;
        for (int k = 0; k < resLength; k++)
        {
            if (state == FS_END)
                break;

            if (fabs(coord[3*k] - point[0]) <= coordTol && fabs(coord[3*k+1] - point[1]) <= coordTol && fabs(coord[3*k+2] - point[2]) <= coordTol)
            {
                found = true;
                kPoint = k;
            }

            bool above = (coord[3*k+2] >= chromoLevel);
            switch (state)
            {
            case FS_VOID:
                if (above)
                {
                    line.start = k;
                    state = FS_IN;
                }
                else 
                    state = FS_OUT;
                break;
            case FS_IN:
                if (!above)
                {
                    line.end = k;
                    state = (found ? FS_END : FS_OUT);
                }
                break;
            case FS_OUT:
                if (above)
                {
                    line.start = k-1;
                    state = FS_IN;
                }
                break;
            }
        }
        if (state == FS_IN)
            line.end = resLength-1;
        else if (state == FS_OUT)
            line.start = line.end = 0;

        if (found && line.end-line.start+1 >= 2)
        {
            resLength = line.end-line.start+1;
            bool bStartClosed = false, bEndClosed = false;
            int sClosed;
            if (coord[3*line.start+2] <= chromoLevel && coord[3*(line.start+1)+2] >= chromoLevel)
                getDivPoint(coord + 3*line.start, coord + 3*(line.start+1));

            if (coord[3*line.start+2] <= chromoLevel + closedTol)
                bStartClosed = true;

            if (coord[3*(line.end-1)+2] >= chromoLevel && coord[3*line.end+2] <= chromoLevel && coord[3*line.end+2] != coord[3*(line.end-1)+2])
                getDivPoint(coord + 3*line.end, coord + 3*(line.end-1));

            if (resLength > 2 || coord[3*line.end + 2] != coord[3 * (line.end - 1) + 2])
            {
                if (coord[3*line.end + 2] <= chromoLevel + closedTol)
                    bEndClosed = true;

                sClosed = (bStartClosed && bEndClosed ? Status::Closed : Status::None) | (bStartClosed != bEndClosed ? Status::OnlyFootpoint : Status::None);

                double B0[3], B1[3], B, Bprev;

                v->getPoint(coord, B0);
                distance[line.start] = 0;

                line.pos_min = line.start;
                double Blng = 0, Bmin = d_snorm(B0);
                Bprev = d_snorm(B0);
                for (int k = line.start + 1; k <= line.end; k++)
                {
                    v->getPoint(coord + 3 * k, B1);
                    B = d_snorm(B1);
                    if (B < Bmin)
                    {
                        Bmin = B;
                        line.pos_min = k;
                    }
                    double dist = d_sdist(coord + 3 * k, coord + 3 * (k - 1));
                    distance[k] = distance[k - 1] + dist;
                    Blng += (Bprev + B)*0.5*dist;

                    Bprev = B;
                }

                line.start_idx = cube->getGlobalID(coord + 3 * line.start);
                line.end_idx = cube->getGlobalID(coord + 3 * line.end);
                line.apex_idx = cube->getGlobalID(coord + 3 * line.pos_min);
                line.seed_idx = cube->getGlobalID(point);
                line.phys_length = distance[line.end];
                line.av_field = Blng / line.phys_length;
                line.status |= Status::Lined | sClosed;
                line.closed = sClosed;
                line.code = _code;
                line.useful = true;
                double dmin = line.distance[line.pos_min];
                for (int k = line.start; k <= line.end; k++)
                {
                    line.indices[k] = cube->getGlobalID(coord + 3 * k);
                    line.distance[k] -= dmin;
                }
            }
        }
    }

    return 0;
}

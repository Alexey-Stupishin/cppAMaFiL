#include "stdDefinitions.h"
#include "LinesTaskQueue.h"
#include "LinesProcessor.h"
#include "agmVectorField.h"
#include "agmVectorFieldLineFuncs.h"
#include "agmRKF45Vect.h"

#define l_div_position(p, N, k1, tk) \
{ \
    if (p >= double((N)-1) || fabs(p-(double((N)-1))) < 1e-5) \
    { \
        k1 = (N)-2; \
        tk = 1; \
    } \
    else \
    { \
        k1 = (int)floor(p); \
        if (k1 < 0) \
        { \
            k1 = 0;  \
            tk = 0; \
        } \
        else \
            tk = (p) - k1; \
    } \
}

#define l_get_point(Ny, tfield, x1, y1, z1, tx, ty, tz) \
    ((1-(tz))* ((1-(ty))* ((1-(tx))*tfield[(y1)     +  (z1)   *Ny][x1] + (tx)*tfield[(y1)     +  (z1)   *Ny][x1+1]) + \
                    (ty)* ((1-(tx))*tfield[((y1)+1) +  (z1)   *Ny][x1] + (tx)*tfield[((y1)+1) +  (z1)   *Ny][x1+1]))  \
   +     (tz)* ((1-(ty))* ((1-(tx))*tfield[(y1)     + ((z1)+1)*Ny][x1] + (tx)*tfield[(y1)     + ((z1)+1)*Ny][x1+1]) + \
                    (ty)* ((1-(tx))*tfield[((y1)+1) + ((z1)+1)*Ny][x1] + (tx)*tfield[((y1)+1) + ((z1)+1)*Ny][x1+1]))  \
    )

bool fcond(void *p, const CagmRKF45Vect& v)
{
    return ((T_Lines *)p)->inBoundCube(v);
}

uint32_t fdata(void *p, const double /*t*/, const CagmRKF45Vect& v, CagmRKF45Vect& vp)
{
    return ((T_Lines *)p)->derivatives(v, vp);
}

//-----------------------------------------------------------------------------
CLinesProcessor::CLinesProcessor(LQPSupervisor *_supervisor, CagmVectorField *_v, int _dir, double _step, double _relErr, double _absErr
        , double _absBoundAchieve, double _relBoundAchieve, int _maxLength, int *_passed, int _n_loop_control, double loop_abs_cell)
    : supervisor(_supervisor)
      ,v(_v)
      ,dir(_dir)
      ,step(_step)
      ,absBoundAchieve(_absBoundAchieve) 
      ,relBoundAchieve(_relBoundAchieve)
      ,maxLength(_maxLength)
      ,passed(_passed)
{
    rkf45 = new CagmRKF45(_absErr, _relErr, fdata, 3, nullptr, fcond, absBoundAchieve, _n_loop_control, loop_abs_cell);
    coord = new double[3*maxLength];
}

//-----------------------------------------------------------------------------
CLinesProcessor::~CLinesProcessor()
{
    delete rkf45;
    delete[] coord;
}

//-------------------------------------------------------------------------------------
uint32_t CLinesProcessor::setTaskParams(void *params)
{
    memcpy(point, (double *)params, 3*sizeof(double)); 

    return 0;
}

//-----------------------------------------------------------------------------
uint32_t CLinesProcessor::ActionCore()
{
    int lineLength;
    int code;

    if ( !supervisor->needProcessing(queueID) )
            return 0;

    aguTimeTicToc tic;
    /* CagmVectorFieldOps::Status s = */
    v->getOneFullLine(rkf45, point, dir, step, absBoundAchieve, relBoundAchieve, maxLength, &lineLength, coord, &code);
    /* uint32_t rc = */
    supervisor->SetResult(queueID, point, coord, lineLength, (int)code, CagmVectorFieldOps::Status::BufferOverload, tic.toc());

    return 0;
}

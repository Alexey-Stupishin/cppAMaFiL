#pragma once

class CWProcCore
{
protected:
    CagmVectorField *field, *baseField, *baseWeight, *baseField2, *baseWeight2;
    CagmScalarField *weight, *absField, *absWeight, *absField2, *absWeight2,
        *losField, *losWeight, *losField2, *losWeight2, *bottomWeight;

    int N[3];

    int depth;
    PROTO_mfoWiegelmannCallback callback;

    CagmVectorField *vF, *boundsx, *boundsy, *boundsz;

    CagmRotate3D *rotator;
    CagpWiegelmann *proc;

    int stepN, stop, iterN, reason;
    double L0, step0, L, Lprev,  step, dL;

    int z_size, voxels;
    double curr_incr, curr_decr = WiegelmannProcStepDecr;

    CagmVectorField *prevField;
    CagmVectorField *prevVF;

    aguTimeTicToc tic;

    double *memoryAv;
    double *memory_dL;
    int n_appr;

public:
    CWProcCore(CagmVectorField *field, CagmScalarField *weight, 
    CagmVectorField *baseField, CagmVectorField *baseWeight, CagmVectorField *baseField2, CagmVectorField *baseWeight2,
    CagmScalarField *absField, CagmScalarField *absWeight, CagmScalarField *absField2, CagmScalarField *absWeight2,
    CagmScalarField *losField, CagmScalarField *losWeight, CagmScalarField *losField2, CagmScalarField *losWeight2,
    CagmScalarField *bottomWeight,
    double *vcos, int depth, PROTO_mfoWiegelmannCallback callback);

    virtual ~CWProcCore();
    int proceed();

protected:
    double get_max_L_incr();
    int proceedDL(double L);
    double init_step();
    void next_field();
    void store_field();
    void restore_field();
};
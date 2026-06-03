#include "stdDefinitions.h"
#include <random>
#include <chrono>

#include "mfoGlobals.h"
#include "MagFieldOps.h"
#include "mfoDisambig.h"

//-----------------------------------------------------------------------
__declspec( dllexport ) int mfoDisambigStep(int *N, int *_state, double *_Fstate, double *_F, double *_kT, int *_prescribe
                                          , int *_natt, int *_nacc)
{
    CmfoDisambig dam(N, _state, nullptr, nullptr, _F, _kT, Disambig_parallel != 0, CommonThreadsN);

    CmfoDisambigRules rules;
    rules.limatt = Disambig_limatt;
    rules.limacc = Disambig_limacc;
    rules.ordpart = Disambig_ordpart;
    dam.init(&rules, _Fstate, _prescribe);

    int try_cnt = dam.temp_step();

    if (_natt)
        *_natt = dam.natt;
    if (_nacc)
        *_nacc = dam.nacc;

    return try_cnt;
}

//-----------------------------------------------------------------------
__declspec(dllexport) int mfoDisambig(int *N, double *absB, double *incl, double *azim, double *Brefx, double *Brefy, double *Brefz, double *step, int *prescribed
                                    , int *_state, int *_generation, int *_status)
{
    int N3[3];
    N3[0] = N[0];
    N3[1] = N[1];
    N3[2] = 1;
    CagmVectorField FIA(N3, absB, incl, azim, false, true);
    CagmVectorField Bref(N3, Brefx, Brefy, Brefz);
    CagmVectorField B000(N3);
    B000.FIA2XYZ(&FIA);
    CagmVectorField B180(N3);
    B180.InvertAzimuth(&B000);

    double *F = FIA.disambigGetF(&B000, &B180, &Bref, step);
    double *kT = FIA.disambigGetT(Disambig_KTFactor_vp, Disambig_KTFactor_v0, Disambig_KTFactor_M, Disambig_KTFactor_p, Disambig_KTFactor_init);

#ifdef _WINDOWS
#ifdef _DEBUG_DISAMB
    char fn[256], buf[16];
    strcpy(fn, path_disamb);
    strcat(fn, "init.bin");
    FILE *fid_da = fopen(fn, "wb");
    CbinDataStruct::WriteHeader(fid_da);
    CbinDataStruct::Write(fid_da, N3, 1, "N0");
    CbinDataStruct::Write(fid_da, N3+1, 1, "N1");
    CbinDataStruct::Write(fid_da, F, (N3[0]+1)*(N3[1]+1)*16, "F");
    CbinDataStruct::Write(fid_da, kT, N3[0]*N3[1], "kT");
    CbinDataStruct::WriteFooter(fid_da);
    fclose(fid_da);
#endif
#endif

    CmfoDisambig dam(N, _state, _generation, absB, F, kT, Disambig_parallel != 0, CommonThreadsN);

    CmfoDisambigRules rules;
    rules.limatt = Disambig_limatt;
    rules.limacc = Disambig_limacc;
    rules.min_generation = Disambig_min_generation;
    rules.ordpart = Disambig_ordpart;
    rules.temp_min = Disambig_temp_min;
    rules.temp_decr = Disambig_temp_decr;
    rules.acc2att = Disambig_acc2att;
    rules.Bmax_term = Disambig_Bmax_term;
    rules.non_stable_term = Disambig_non_stable_term;

    dam.create(&rules, prescribed);

    int rc = dam.process();
    if (_status)
        *_status = (int)dam.status;

    return rc;
}

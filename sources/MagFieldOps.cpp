// MagFieldOps.cpp : Defines the initialization routines for the DLL.
//

#include "stdDefinitions.h"

#ifdef _WINDOWS
#include <windows.h>
#pragma warning (disable : 4996)
#endif

#include <ctime>
#include "string.h"

#include "MagFieldOps.h"

#include "mfoGlobals.h"
#include "VersionInfoResource.h"

#include "agmScalarField.h"
#include "agmVectorField.h"

//------------------------------------------------------------------
bool mapIntProceed(bool bGet, std::string name, int &value, int defaultValue = 0)
{
    auto search = mapInt.find(name);
    bool rc = search != mapInt.end();
    if (bGet) // return true, if exist
    {
        if (rc)
            value = search->second;
        else
            value = defaultValue;
    }
    else
    {
        mapInt[name] = defaultValue;
    }

    return rc;
}

//------------------------------------------------------------------
bool mapDoubleProceed(bool bGet, std::string name, double &value, double defaultValue = 0)
{
    auto search = mapDouble.find(name);
    bool rc = search != mapDouble.end();
    if (bGet) // return true, if exist
    {
        if (rc)
            value = search->second;
        else
            value = defaultValue;
    }
    else
    {
        mapDouble[name] = defaultValue;
    }

    return rc;
}

//------------------------------------------------------------------
bool mapSetProceed(std::string name, double &value, double defaultValue = 0)
{
    auto searchD = mapDouble.find(name);
    bool rc = searchD != mapDouble.end();
    if (rc) // exists in Double
        mapDoubleProceed(false, name, value, defaultValue);
    else
    {
        auto searchI = mapInt.find(name);
        rc = searchI != mapInt.end();
        if (rc) // exists in Int
        {
            int iv = (int)value;
            mapIntProceed(false, name, iv, (int)defaultValue);
        }
    }

    return rc;
}

//------------------------------------------------------------------
__declspec( dllexport ) uint32_t utilInitialize()
{
    _proceedGlobals(false);
    _proceedGlobals();
    return 0;
}

//------------------------------------------------------------------
__declspec( dllexport ) int utilSetInt(char *query, int value)
{
    bool exists = mapIntProceed(false, query, value, value);
    _proceedGlobals();
    return (exists ? 1 : 0);
}

//------------------------------------------------------------------
__declspec( dllexport ) int utilGetInt(char *query, int *result)
{ 
    return mapIntProceed(true, query, *result, 0);
}

//------------------------------------------------------------------
__declspec( dllexport ) int utilSetDouble(char *query, double value)
{
    bool exists = mapDoubleProceed(false, query, value, value);
    _proceedGlobals();
    return (exists ? 1 : 0);
}

//------------------------------------------------------------------
__declspec( dllexport ) int utilGetDouble(char *query, double *result)
{ 
    return mapDoubleProceed(true, query, *result, 0);
}

//------------------------------------------------------------------
__declspec( dllexport ) int utilSetSetting(char *query, double value)
{
    bool exists = mapSetProceed(query, value, value);
    _proceedGlobals();
    return (exists ? 1 : 0);
}

//------------------------------------------------------------------
void _proceedGlobals(bool bGet)
{
    mapIntProceed   (bGet, "n_processes", CommonThreadsN, 0);

    int wp;
    mapIntProceed   (bGet, "threads_priority", wp, (int)w_priority::low);
    WiegelmannThreadPriority = (w_priority)wp;

    mapIntProceed   (bGet, "bounds_correction", WiegelmannBoundsCorrection, 0);

    mapIntProceed   (bGet, "weight_type", WiegelmannWeightType, SWF_COS);
    mapDoubleProceed(bGet, "weight_bound_size", WiegelmannWeightBound, 0.1);
    mapDoubleProceed(bGet, "weight_relative_divfree", WiegelmannWeightDivfree, 1.0);

    mapIntProceed   (bGet, "derivative_stencil", WiegelmannDerivStencil, 3);
    mapDoubleProceed(bGet, "inversion_tolerance", WiegelmannInversionTolerance, 0);
    mapDoubleProceed(bGet, "inversion_denominator", WiegelmannInversionDenom, 0);
    
    mapDoubleProceed(bGet, "step_initial", WiegelmannProcStep0, 0.1); // units of Bav/F2max
    
    mapIntProceed   (bGet, "step_max", WiegelmannProcStepMax, 1000000000); // max. aval. step
    mapIntProceed   (bGet, "max_iterations", WiegelmannProcMaxSteps, 100000); // max. number of step

    mapDoubleProceed(bGet, "step_increment_init", WiegelmannProcStepIncrInit, 1.618);
    mapDoubleProceed(bGet, "step_increment", WiegelmannProcStepIncrMatr, 1.618);
    mapDoubleProceed(bGet, "step_increment_main", WiegelmannProcStepIncrMain, 1.618);
    mapDoubleProceed(bGet, "step_decrement_init", WiegelmannProcStepDecrInit, 0.382);
    mapDoubleProceed(bGet, "step_decrement", WiegelmannProcStepDecrMatr, 0.382);
    mapDoubleProceed(bGet, "step_decrement_main", WiegelmannProcStepDecrMain, 0.382);
    mapDoubleProceed(bGet, "step_terminate_init", WiegelmannProcStepLimInit, 0.01); // related to init. step
    mapDoubleProceed(bGet, "step_terminate", WiegelmannProcStepLimMatr, 0.01);
    mapDoubleProceed(bGet, "step_terminate_main", WiegelmannProcStepLimMain, 0.0001);

    mapDoubleProceed(bGet, "step_increment_factor", WiegelmannProcStepIncrFactor, 1.0);
    mapDoubleProceed(bGet, "step_decrement_factor", WiegelmannProcStepDecrFactor, 1.0);

    mapDoubleProceed(bGet, "d_functional_stdev_value_init", WiegelmannProcdLStdValInit, 5e-4); // min std(rel. dL) ...
    mapDoubleProceed(bGet, "d_functional_stdev_value", WiegelmannProcdLStdValMatr, 5e-4);
    mapDoubleProceed(bGet, "d_functional_stdev_value_main", WiegelmannProcdLStdValMain, 5e-4);
    mapIntProceed   (bGet, "d_functional_stdev_window_init", WiegelmannProcdLStdWinInit, 101); // ... std window size + 1
    mapIntProceed   (bGet, "d_functional_stdev_window_matr", WiegelmannProcdLStdWinMatr, 101);
    mapIntProceed   (bGet, "d_functional_stdev_window_main", WiegelmannProcdLStdWinMain, 101);

    mapIntProceed   (bGet, "dense_grid_use", WiegelmannMatryoshkaUse, 1);
    mapIntProceed   (bGet, "dense_grid_min_n", WiegelmannMatryoshkaDeepMinN, 25);
    mapDoubleProceed(bGet, "dense_grid_factor", WiegelmannMatryoshkaFactor, 2.0);

    mapIntProceed   (bGet, "add_conditions_mode", WiegelmannProcCondType,  1); // 0 - ignore, 1 - as functional, 2 - as fixed

    mapIntProceed   (bGet, "add_conditions_abs", WiegelmannProcCondAbs,   2);  // 0 - no limit, 1 - no less, 2 - exact
    mapIntProceed   (bGet, "add_conditions_abs_max", WiegelmannProcCondAbs2,  1); // 0 - no limit, 1 - no greater
    mapIntProceed   (bGet, "add_conditions_los", WiegelmannProcCondLOS,   2);
    mapIntProceed   (bGet, "add_conditions_los_max", WiegelmannProcCondLOS2,  1);
    mapIntProceed   (bGet, "add_conditions_xyz", WiegelmannProcCondBase,  2);
    mapIntProceed   (bGet, "add_conditions_xyz_max", WiegelmannProcCondBase2, 1);

    mapIntProceed   (bGet, "protocol_step", WiegelmannProtocolStep, 10);

    mapIntProceed   (bGet, "metrics_theta", WiegelmannGetMetricsTheta, 0);
    mapIntProceed   (bGet, "debug_input", debug_input, 0);

    mapIntProceed(bGet, "disambig_parallel", Disambig_parallel, 0);
    mapIntProceed(bGet, "disambig_limacc", Disambig_limacc, 3);
    mapIntProceed(bGet, "disambig_limatt", Disambig_limatt, 5);
    mapIntProceed(bGet, "disambig_min_generation", Disambig_min_generation, 100);

    mapDoubleProceed(bGet, "disambig_ktfactor_vp", Disambig_KTFactor_vp, 1.0);
    mapDoubleProceed(bGet, "disambig_ktfactor_v0", Disambig_KTFactor_v0, 10.0);
    mapDoubleProceed(bGet, "disambig_ktfactor_M", Disambig_KTFactor_M, 100.0);
    mapDoubleProceed(bGet, "disambig_ktfactor_p", Disambig_KTFactor_p, 0.25);
    mapDoubleProceed(bGet, "disambig_ktfactor_init", Disambig_KTFactor_init, 0.5);
    mapDoubleProceed(bGet, "disambig_ordpart", Disambig_ordpart, 1.0);
    mapDoubleProceed(bGet, "disambig_temp_decr", Disambig_temp_decr, 0.99);
    mapDoubleProceed(bGet, "disambig_temp_min", Disambig_temp_min, 0.005);
    mapDoubleProceed(bGet, "disambig_acc2att", Disambig_acc2att, 0.005);
    mapDoubleProceed(bGet, "disambig_Bmax_term", Disambig_Bmax_term, 200);
    mapDoubleProceed(bGet, "disambig_non_stable_term", Disambig_non_stable_term, 0.005);

    mapIntProceed(bGet, "lines_conditions", lines_conditions, 3);
    mapDoubleProceed(bGet, "lines_chromo_level", lines_chromo_level, 1);
    mapDoubleProceed(bGet, "lines_step", lines_step, 1);
    mapDoubleProceed(bGet, "lines_tolerance", lines_tolerance, 1e-4);
    mapDoubleProceed(bGet, "lines_abs_bound_achieve", lines_abs_bound_achieve, 1e-3);
    mapDoubleProceed(bGet, "lines_rel_bound_achieve", lines_rel_bound_achieve, 0);
    mapDoubleProceed(bGet, "lines_rel_seeds_bound", lines_rel_seeds_bound, 0);
    mapIntProceed(bGet, "lines_n_loop_control", lines_n_loop_control, 200);
    mapDoubleProceed(bGet, "lines_loop_abs_cell", lines_loop_abs_cell, 1.5);
}

//------------------------------------------------------------------
__declspec(dllexport) int utilGetVersion(char *fullvers, int buflength)
{
    std::string s;
    s = VIR_ProductName;
    s += " v.";
    s += VIR_QUOTE_SUBST(VIR_Ver1);
    s += ".";
    s += VIR_QUOTE_SUBST(VIR_Ver2);
    s += ".";
    s += VIR_QUOTE_SUBST(VIR_Ver3);
    s += ".";
    s += VIR_QUOTE_SUBST(VIR_Ver4);
    s += " (";
#ifdef _WINDOWS
    s += VIR_QUOTE_SUBST(VIR_T_REV);
#else
    s += VIR_T_REV;
#endif
    s += VIR_QUOTE_SUBST(VIR_Revision);
    s += "). ";
#ifdef _WINDOWS
    s += VIR_QUOTE_SUBST(VIR_COPYRIGHT);
#else
    s += VIR_COPYRIGHT;
#endif
//    s += ", ";
#ifdef _WINDOWS
    s += VIR_QUOTE_SUBST(VIR_FROM);
#else
    s += VIR_FROM;
#endif
    s += VIR_QUOTE_SUBST(VIR_Year);
    s += ", ";
    s += VIR_CompanyName;
    
    strncpy(fullvers, s.c_str(), buflength);

    return (int)s.length();
}

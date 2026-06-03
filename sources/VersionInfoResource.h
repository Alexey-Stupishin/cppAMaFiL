#pragma once

#define VIR_QUOTE(a) #a
#define VIR_QUOTE_SUBST(a) VIR_QUOTE(a)

//-------------------------------------------------------------------------
#define VIR_Ver1 4
#define VIR_Ver2 4
#define VIR_Ver3 26
#define VIR_Ver4 601
#define VIR_Revision 61
#define VIR_Year 2026

//-------------------------------------------------------------------------
#define VIR_CompanyName "Special Astrophysical Observatory, Russia"
#define VIR_InternalName "MagFieldOps.dll"
#define VIR_OriginalFilename "WWNLFFFReconstruction.dll"
#define VIR_ProductName "Advanced Magnetic Field Library"

#ifdef _WINDOWS
#define VIR_T_REV rev.
#define VIR_T_LIBNAME : Weighted Wiegelmann NLFFF Reconstruction Library
#define VIR_COPYRIGHT Copyright (C) Alexey G. Stupishin (agstup@yandex.ru)
#define VIR_FROM , 2017-
#else
#define VIR_T_REV "rev."
#define VIR_T_LIBNAME ": Advanced Magnetic Field Library"
#define VIR_COPYRIGHT "Copyright (C) Alexey G. Stupishin (agstup@yandex.ru)"
#define VIR_FROM ", 2017-"
#endif

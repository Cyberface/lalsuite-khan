/*
*  Copyright (C) 2007 Jolien Creighton, Josh Willis
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <config.h>

#include <complex.h>
#include <mkl_dfti.h>

#include <lal/SeqFactories.h>
#include <lal/ComplexFFT.h>

/** \cond DONT_DOXYGEN */


#define CHECKINTELFFTSTATUS( fstat )                    \
  if ( (fstat) != DFTI_NO_ERROR )                       \
  {                                                     \
    char *errmsg = DftiErrorMessage( fftStat );         \
    XLAL_ERROR( COMPLEXFFTH_EINTL, "%s", errmsg );      \
  }                                                     \
  else (void)(0)

#define CHECKINTELFFTSTATUS_NULL( fstat )               \
  if ( (fstat) != DFTI_NO_ERROR )                       \
  {                                                     \
    char *errmsg = DftiErrorMessage( fftStat );         \
    XLAL_ERROR_NULL( COMPLEXFFTH_EINTL, "%s", errmsg ); \
  }                                                     \
  else (void)(0)

#define CHECKINTELFFTSTATUS_VOID( fstat )               \
  if ( (fstat) != DFTI_NO_ERROR )                       \
  {                                                     \
    char *errmsg = DftiErrorMessage( fftStat );         \
    XLAL_ERROR_VOID( COMPLEXFFTH_EINTL, "%s", errmsg ); \
  }                                                     \
  else (void)(0)


/*
 * Plan to perform FFT of COMPLEX8 data.
 */
struct
tagCOMPLEX8FFTPlan
{
  INT4       sign; /* sign in transform exponential, -1 for forward, +1 for reverse */
  UINT4      size; /* length of the real data vector for this plan */
  DFTI_DESCRIPTOR *plan; /* the MKL plan */
};

/*
 * Plan to perform FFT of COMPLEX16 data.
 */
struct
tagCOMPLEX16FFTPlan
{
  INT4       sign; /* sign in transform exponential, -1 for forward, +1 for reverse */
  UINT4      size; /* length of the real data vector for this plan */
  DFTI_DESCRIPTOR *plan; /* the MKL plan */
};

/* single- and double-precision routines */

#define SINGLE_PRECISION
#include "IntelComplexFFT_source.c"
#undef SINGLE_PRECISION
#include "IntelComplexFFT_source.c"


/** \endcond */

/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Kipp Cannon, Patrick Brady, Tania Regimbau
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

#ifndef _TIMEFREQFFT_H
#define _TIMEFREQFFT_H

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \defgroup TimeFreqFFT_h Header TimeFreqFFT_h
 * \ingroup lal_fft
 *
 * \brief Performs real-to-complex, complex-to-real FFTs and average power spectrum estimation.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/TimeFreqFFT.h>
 * \endcode
 *
 * Perform time-to-frequency and frequency-to-time fast Fourier
 * transforms. Also provides a function to compute mean and median power
 * spectra with user specified windowning.
 *
 * ### Description ###
 *
 * The routines LALTimeFreqRealFFT() and LALTimeFreqComplexFFT()
 * transform time series \f$h_j\f$, \f$0\le j<n\f$, into a frequency series
 * \f$\tilde{h}_k\f$.  For LALTimeFreqRealFFT(),
 * \f[
 * \tilde{h}_k = \Delta t \times H_k \;
 * \mbox{for } 0\le k\le\lfloor n/2\rfloor
 * \f]
 * The packing covers the range from dc (inclusive) to Nyquist (inclusive if
 * \f$n\f$ is even).
 * For LALTimeFreqComplexFFT(),
 * \f[
 * \tilde{h}_k = \Delta t \left\{
 * \begin{array}{ll}
 * H_{k+\lfloor(n+1)/2\rfloor} &
 * \mbox{for } 0\le k<\lfloor n/2\rfloor, \\
 * H_{k-\lfloor n/2\rfloor} &
 * \mbox{for } \lfloor n/2\rfloor\le k<n. \\
 * \end{array}
 * \right.
 * \f]
 * The packing covers the range from negative Nyquist (inclusive if \f$n\f$ is
 * even) up to (but not including) positive Nyquist.
 * Here \f$H_k\f$ is the DFT of \f$h_j\f$:
 * \f[
 * H_k = \sum_{j=0}^{n-1} h_j e^{-2\pi ijk/n}.
 * \f]
 * The units of \f$\tilde{h}_k\f$ are equal to the units of \f$h_j\f$ times seconds.
 *
 * The routines LALFreqTimeRealFFT() and LALFreqTimeComplexFFT()
 * perform the inverse transforms from \f$\tilde{h}_k\f$ back to \f$h_j\f$.  This is
 * done by shuffling the data, performing the reverse DFT, and multiplying by
 * \f$\Delta f\f$.
 *
 * The routine LALREAL4AverageSpectrum() uses Welch's method to compute
 * the average power spectrum of the time series stored in the input structure
 * \c tSeries and return it in the output structure \c fSeries.  A
 * Welch PSD estimate is defined by an FFT length, overlap length, choice of
 * window function and averaging method. These are specified in the
 * parameter structure; the FFT length is obtained from the length of the
 * \c REAL4Window in the parameters.
 *
 * On entry the parameter structure \c params must contain a valid
 * \c REAL4Window, an integer that determines the overlap as described
 * below and a forward FFT plan for transforming data of the specified
 * window length into the time domain. The method used to compute the
 * average must also be set.
 *
 * If the length of the window is \f$N\f$, then the FFT length is defined to be
 * \f$N/2-1\f$. The input data of length \f$M\f$ is divided into \f$i\f$ segments which
 * overlap by \f$o\f$, where
 * \f{equation}{
 * i = \frac{M-o}{N-o}.
 * \f}
 *
 * The PSD of each segment is obtained. The Welch PSD estimate is the average
 * of these \f$i\f$ sub-estimates.  The average is computed using the mean or
 * median method, as specified in the parameter structure.
 *
 * Note: the return PSD estimate is a one-sided power spectral density
 * normalized as defined in the conventions document. When the averaging
 * method is choosen to be mean and the window type Hann, the result is the
 * same as returned by the LDAS datacondAPI <tt>psd()</tt> action for a real
 * sequence without detrending.
 *
 * ### Operating Instructions ###
 *
 * \code
 * const UINT4 n  = 65536;
 * const REAL4 dt = 1.0 / 16384.0;
 * static LALStatus status; compute average power spectrum
 * static REAL4TimeSeries         x;
 * static COMPLEX8FrequencySeries X;
 * static COMPLEX8TimeSeries      z;
 * static COMPLEX8FrequencySeries Z;
 * RealFFTPlan    *fwdRealPlan    = NULL;
 * RealFFTPlan    *revRealPlan    = NULL;
 * ComplexFFTPlan *fwdComplexPlan = NULL;
 * ComplexFFTPlan *revComplexPlan = NULL;
 *
 * LALSCreateVector( &status, &x.data, n );
 * LALCCreateVector( &status, &X.data, n / 2 + 1 );
 * LALCCreateVector( &status, &z.data, n );
 * LALCCreateVector( &status, &Z.data, n );
 * LALCreateForwardRealFFTPlan( &status, &fwdRealPlan, n, 0 );
 * LALCreateReverseRealFFTPlan( &status, &revRealPlan, n, 0 );
 * LALCreateForwardComplexFFTPlan( &status, &fwdComplexPlan, n, 0 );
 * LALCreateReverseComplexFFTPlan( &status, &revComplexPlan, n, 0 );
 *
 * x.f0 = 0;
 * x.deltaT = dt;
 * x.sampleUnits = lalMeterUnit;
 * strncpy( x.name, "x", sizeof( x.name ) );
 *
 * z.f0 = 0;
 * z.deltaT = dt;
 * z.sampleUnits = lalVoltUnit;
 * strncpy( z.name, "z", sizeof( z.name ) );
 *
 * <assign data>
 *
 * LALTimeFreqRealFFT( &status, &X, &x, fwdRealPlan );
 * LALFreqTimeRealFFT( &status, &x, &X, revRealPlan );
 * LALTimeFreqComplexFFT( &status, &Z, &z, fwdComplexPlan );
 * LALFreqTimeComplexFFT( &status, &z, &Z, revComplexPlan );
 *
 * LALDestroyRealFFTPlan( &status, &fwdRealPlan );
 * LALDestroyRealFFTPlan( &status, &revRealPlan );
 * LALDestroyComplexFFTPlan( &status, &fwdComplexPlan );
 * LALDestroyComplexFFTPlan( &status, &revComplexPlan );
 * LALCDestroyVector( &status, &Z.data );
 * LALCDestroyVector( &status, &z.data );
 * LALCDestroyVector( &status, &X.data );
 * LALSDestroyVector( &status, &x.data );
 * \endcode
 *
 * ### Notes ###
 *
 * <ol>
 * <li> The routines do not presently work properly with heterodyned data,
 * i.e., the original time series data should have \c f0 equal to zero.
 * </li></ol>
 *
 */
/*@{*/

/** UNDOCUMENTED */
typedef struct
tagLALPSDRegressor
{
  unsigned average_samples;
  unsigned median_samples;
  unsigned n_samples;
  REAL8Sequence **history;
  REAL8FrequencySeries *mean_square;
}
LALPSDRegressor;

/*
 *
 * XLAL Functions
 *
 */
int XLALREAL4TimeFreqFFT(
    COMPLEX8FrequencySeries *freq,
    const REAL4TimeSeries   *tser,
    const REAL4FFTPlan      *plan
    );

int XLALREAL4FreqTimeFFT(
    REAL4TimeSeries               *tser,
    const COMPLEX8FrequencySeries *freq,
    const REAL4FFTPlan            *plan
    );

int XLALREAL8TimeFreqFFT(
    COMPLEX16FrequencySeries *freq,
    const REAL8TimeSeries    *tser,
    const REAL8FFTPlan       *plan
    );

int XLALREAL8FreqTimeFFT(
    REAL8TimeSeries                *tser,
    const COMPLEX16FrequencySeries *freq,
    const REAL8FFTPlan             *plan
    );

int XLALCOMPLEX8TimeFreqFFT(
    COMPLEX8FrequencySeries  *freq,
    const COMPLEX8TimeSeries *tser,
    const COMPLEX8FFTPlan    *plan
    );

int XLALCOMPLEX8FreqTimeFFT(
    COMPLEX8TimeSeries            *tser,
    const COMPLEX8FrequencySeries *freq,
    const COMPLEX8FFTPlan         *plan
    );

int XLALCOMPLEX16TimeFreqFFT(
    COMPLEX16FrequencySeries  *freq,
    const COMPLEX16TimeSeries *tser,
    const COMPLEX16FFTPlan    *plan
    );

int XLALCOMPLEX16FreqTimeFFT(
    COMPLEX16TimeSeries            *tser,
    const COMPLEX16FrequencySeries *freq,
    const COMPLEX16FFTPlan         *plan
    );

int XLALREAL4ModifiedPeriodogram(
    REAL4FrequencySeries        *periodogram,
    const REAL4TimeSeries       *tseries,
    const REAL4Window           *window,
    const REAL4FFTPlan          *plan
    );

int XLALREAL8ModifiedPeriodogram(
    REAL8FrequencySeries        *periodogram,
    const REAL8TimeSeries       *tseries,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan
    );

int XLALREAL4AverageSpectrumWelch(
    REAL4FrequencySeries        *spectrum,
    const REAL4TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL4Window           *window,
    const REAL4FFTPlan          *plan
    );

int XLALREAL8AverageSpectrumWelch(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan
    );

REAL8 XLALMedianBias( UINT4 nn );

REAL8 XLALLogMedianBiasGeometric( UINT4 nn );

int XLALREAL4AverageSpectrumMedian(
    REAL4FrequencySeries        *spectrum,
    const REAL4TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL4Window           *window,
    const REAL4FFTPlan          *plan
    );

int XLALREAL8AverageSpectrumMedian(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan
    );

int XLALREAL4AverageSpectrumMedianMean(
    REAL4FrequencySeries        *spectrum,
    const REAL4TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL4Window           *window,
    const REAL4FFTPlan          *plan
    );

int XLALREAL8AverageSpectrumMedianMean(
    REAL8FrequencySeries        *spectrum,
    const REAL8TimeSeries       *tseries,
    UINT4                        seglen,
    UINT4                        stride,
    const REAL8Window           *window,
    const REAL8FFTPlan          *plan
    );

int XLALREAL4SpectrumInvertTruncate(
    REAL4FrequencySeries        *spectrum,
    REAL4                        lowfreq,
    UINT4                        seglen,
    UINT4                        trunclen,
    REAL4FFTPlan                *fwdplan,
    REAL4FFTPlan                *revplan
    );

int XLALREAL8SpectrumInvertTruncate(
    REAL8FrequencySeries        *spectrum,
    REAL8                        lowfreq,
    UINT4                        seglen,
    UINT4                        trunclen,
    REAL8FFTPlan                *fwdplan,
    REAL8FFTPlan                *revplan
    );

REAL4TimeSeries *XLALRespFilt(
    REAL4TimeSeries             *strain,
    COMPLEX8FrequencySeries     *transfer
    );

REAL4TimeSeries *XLALREAL4Convolution(
    REAL4TimeSeries             *strain,
    REAL4TimeSeries             *transfer
    );

COMPLEX8FrequencySeries *XLALWhitenCOMPLEX8FrequencySeries(
    COMPLEX8FrequencySeries     *fseries,
    const REAL4FrequencySeries  *psd
    );

COMPLEX16FrequencySeries *XLALWhitenCOMPLEX16FrequencySeries(
    COMPLEX16FrequencySeries    *fseries,
    const REAL8FrequencySeries  *psd
);

REAL8Sequence *XLALREAL8WindowTwoPointSpectralCorrelation(
    const REAL8Window *window,
    const REAL8FFTPlan *plan
);

LALPSDRegressor *
XLALPSDRegressorNew(
    unsigned average_samples,
    unsigned median_samples
);

void
XLALPSDRegressorFree(
    LALPSDRegressor *r
);

void
XLALPSDRegressorReset(
    LALPSDRegressor *r
);

int XLALPSDRegressorSetAverageSamples(
    LALPSDRegressor *r,
    unsigned average_samples
);

unsigned XLALPSDRegressorGetAverageSamples(
    const LALPSDRegressor *r
);

unsigned XLALPSDRegressorGetNSamples(
    const LALPSDRegressor *r
);

int XLALPSDRegressorSetMedianSamples(
    LALPSDRegressor *r,
    unsigned median_samples
);

unsigned XLALPSDRegressorGetMedianSamples(
    const LALPSDRegressor *r
);

int
XLALPSDRegressorAdd(
    LALPSDRegressor *r,
    const COMPLEX16FrequencySeries *sample
);

REAL8FrequencySeries *
XLALPSDRegressorGetPSD(
    const LALPSDRegressor *r
);

int
XLALPSDRegressorSetPSD(
    LALPSDRegressor *r,
    const REAL8FrequencySeries *psd,
    unsigned weight
);


/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _TIMEFREQFFT_H */

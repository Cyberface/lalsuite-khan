/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpFilterOutputVeto.c
 *
 * Author: A.C. Rodriguez
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpFilterOutputVetoCV">
Author: Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpFilterOutputVeto.c}}
\label{ss:FindChirpFilterOutputVeto.c}

Memory management functions for creating and destroying input data and 
workspace memory for findchirp.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpFilterOutputVetoCP}
\idx{LALFindChirpFilterOutputVeto()}

\subsubsection*{Description}

The function \texttt{LALFindChirpFilterOutputVeto()} implements a signal based   
veto, currently it is used primarily for testing. The function itself tests the 
consistency of the triggers that survive the bns and macho inspiral search pipeline 
by monitoring the behavior of the $r^2$ time series which is calculated for each 
segment of data (256 seconds). \\

  Two thresholds are currently employed in the bns and macho inspiral searches: a 
 signal to noise, $\rho$, threshold and a threshold on the consistency of the 
 template chirp waveform with the data, $r^2$. The search code compares 
 $\rho$(\emph{t$_j$}) and $r^2$(\emph{t$_j$}) for a given segment of data and looks 
 for $\rho(\emph{t$_j$})$ $>$ $\rho^*(\emph{t$_j$})$ and $r^2(\emph{t$_j$})$ $<$ 
 $r^{2*}(\emph{t$_j$})$, where $^*$ denotes the threshold used in the search. If 
 both these criteria are met at a given \emph{t$_j$}, an inspiral "trigger" is 
 recorded. 
\\
\\
$r^2$(\emph{t$_j$}) is defined as:

\begin{equation}
{r^2(\emph{t$_j$})} = \frac{\chi^2(\emph{t$_j$})}{p + \delta^2\rho^2(\emph{t$_j$})}  
\end{equation}
\\  
where:
\\
\\
p $\equiv$  number of  $\chi^2$ bins\\
$\delta$ = mismatch parameter between the template and data\\ 
$\rho$ = the signal to noise ratio. 

\subsubsection*{Algorithm}
  The algorithm inputs the the vector  \texttt{chisqVec} (which is actually $r^2$) 
  for the whole data segment and searches a time window prior to the inferred 
  coalescence time (\texttt{rsqvetoWindow}) of the trigger up to the trigger time 
  and counts the number of time samples above a given $\chi^2$ threshold, looser 
  than the search pipeline employs (\texttt{rsqvetoThresh}). The value is currently 
  stored under the \texttt{rsqveto\_duration} field in the \texttt{sngl\_inspiral} 
  xml table. Future implementation of this function will have it take the 
  calculated value and decide whether or not to store the trigger for future 
  analysis.
\subsubsection*{Uses}

\subsubsection*{Notes}
  The same test described here could also be employed for monitoring the behavior 
  of the signal to noise time series, $\rho(\emph{t$_j$})$, about a trigger, 
  therefore the inclusion of  \texttt{qVec} as input to the function for future work.

\vfill{\footnotesize\input{FindChirpFilterOutputVetoCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpFilterOutputVeto.h>

double rint(double x);

NRCSID (FINDCHIRPFILTEROUTPUTVETOC, "$Id$");

/* <lalVerbatim file="FindChirpFilterOutputVetoCP"> */
void 
LALFindChirpFilterOutputVeto( 
    LALStatus                          *status,
    SnglInspiralTable                 **eventList, 
    FindChirpSegment                   *segment,
    REAL4Vector                        *chisqVec,
    REAL8                               deltaT,
    COMPLEX8Vector                     *qVec,
    REAL4                               qNorm,
    FindChirpFilterOutputVetoParams    *params
    )
/* </lalVerbatim> */
{

  UINT4                 x; /* for loops */
  UINT4		        eventIndex; /* store the event as an index */
  UINT4		        windowIndex; /* r^2 veto window index */ 
  UINT4                 timeaboversqThresh; /* time spent above the r^2 veto threshold */
  REAL4                 rsqvetoWindow; /* the r^2 veto window */
  REAL4		        rsqvetoThresh; /* the r^2 veto threshold */	
  SnglInspiralTable    *event = NULL; 
  event = *eventList;

  INITSTATUS( status, "LALFindChirpFilterFilterOutputVeto", 
      FINDCHIRPFILTEROUTPUTVETOC );
  ATTATCHSTATUSPTR( status );


  /*
   *    
   *   check that the arguments are reasonable
   *          
   */
 	


  /* check that the filterOutputVeto parameter structure exist */	
  ASSERT( params->rsqvetoWindow, status, 
      FINDCHIRPFILTEROUTPUTVETOH_ENULL, FINDCHIRPFILTEROUTPUTVETOH_MSGENULL);
  ASSERT( params->rsqvetoThresh, status, 
      FINDCHIRPFILTEROUTPUTVETOH_ENULL, FINDCHIRPFILTEROUTPUTVETOH_MSGENULL);
	
  /* check that the workspace vectors exist */
  ASSERT( chisqVec->data, status, 
      FINDCHIRPFILTEROUTPUTVETOH_ENULL, FINDCHIRPFILTEROUTPUTVETOH_MSGENULL);
	
  /* check that deltaT is reasonable */	
  ASSERT( deltaT > 0, status,
      FINDCHIRPFILTEROUTPUTVETOH_EDTZO, FINDCHIRPFILTEROUTPUTVETOH_MSGEDTZO );
			
  /* make sure that the segment structure exists */
  ASSERT( segment, status, 
      FINDCHIRPFILTEROUTPUTVETOH_ENULL, FINDCHIRPFILTEROUTPUTVETOH_MSGENULL );

  /* make sure that the segment structure contains some input */
  ASSERT( segment->data->epoch.gpsSeconds, status, 
      FINDCHIRPFILTEROUTPUTVETOH_ENULL, FINDCHIRPFILTEROUTPUTVETOH_MSGENULL );
	
  /* If an event exists, begin the r^2 veto calculations */  
  while ( event )
  {
    /* convert the time of the event into an index */
    eventIndex = (UINT4) floor( ( (REAL8) ( event->end_time.gpsSeconds - 
            segment->data->epoch.gpsSeconds ) + 
          (REAL8) ( event->end_time.gpsNanoSeconds + 1 ) * 1.0e-9 ) / deltaT );

    /* convert the window duration into an index */
    windowIndex = (UINT4) floor( (REAL8) ( ( params->rsqvetoWindow ) / deltaT ) );
											
    /* Initialize output to zero */ 
    timeaboversqThresh = 0;
	     
    /* Count the number of time samples above the given r^2 veto threshold */
    for( x = eventIndex - windowIndex; x <= eventIndex; ++x )
    {
      if ( chisqVec->data[x] >= params->rsqvetoThresh )
          ++timeaboversqThresh;
    }
		
    /* Convert the output to seconds and record the computed values in the 
       sngl_inspiral event */
    event->rsqveto_duration = (REAL4) timeaboversqThresh * deltaT;
    event = event->next;
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

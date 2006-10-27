/*
 * Copyright (C) 2004, 2005 R. Prix, B. Machenschalk, A.M. Sintes
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

/** \file
 * \ingroup SFTfileIO
 * \author R. Prix, B. Machenschalk, A.M. Sintes, B. Krishnan
 * 
 * \brief IO library for reading/writing "Short Fourier transform" (SFT) data files.
 *
 * $Id$ 
 *
 */

/*---------- INCLUDES ----------*/
#include <sys/types.h>
#include <errno.h>
#include <string.h>
#include <strings.h>

#ifndef _MSC_VER
#include <dirent.h>
#else
#include <io.h>
#endif

#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include "lal/SFTfileIO.h"

NRCSID (SFTFILEIOC, "$Id$");
/*---------- DEFINES ----------*/

#define MIN_SFT_VERSION 1
#define MAX_SFT_VERSION 2

#define TRUE    1
#define FALSE   0

/** blocksize used in SFT-reading for the CRC-checksum computation (has to be multiple of 8 !!) */
#define BLOCKSIZE 8192 * 8

/*----- Macros ----- */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

/* rounding for positive numbers! */
#define MYROUND(x) ( floor( (x) + 0.5 ) )
      
/*---------- internal types ----------*/

/* NOTE: the locator is implemented as an OPQUE type in order to enforce encapsulation 
 * of the actual physical storage of SFTs and to ease future extensions of the interface.
 * DO NOT TRY TO USE THIS TYPE OUTSIDE OF THIS MODULE!!
 */
struct tagSFTLocator 
{
  CHAR *fname;		/* name of file containing this SFT */
  long offset;		/* SFT-offset with respect to a merged-SFT */
};

typedef struct
{
  REAL8 version;
  INT4 gps_sec;
  INT4 gps_nsec;
  REAL8 tbase;
  INT4 first_frequency_index;
  INT4 nsamples;
} _SFT_header_v1_t;

typedef struct
{
  REAL8 version;
  INT4 gps_sec;
  INT4 gps_nsec;
  REAL8 tbase;
  INT4 first_frequency_index;
  INT4 nsamples;
  UINT8 crc64;
  CHAR detector[2];
  CHAR padding[2];
  INT4 comment_length;
} _SFT_header_v2_t;


/*---------- Global variables ----------*/
/* empty struct initializers */
static LALStatus empty_status;
const SFTConstraints empty_SFTConstraints;
const SFTCatalog empty_SFTCatalog;

/*---------- internal prototypes ----------*/
static LALStringVector *find_files (const CHAR *fpattern);

static void DestroyStringVector (LALStringVector *strings);
static void SortStringVector (LALStringVector *strings);
static void endian_swap(CHAR * pdata, size_t dsize, size_t nelements);
static int amatch(char *str, char *p);	/* glob pattern-matcher (public domain)*/

static BOOLEAN is_valid_detector (const char *channel);
static BOOLEAN consistent_mSFT_header ( SFTtype header1, UINT4 version1, UINT4 nsamples1, SFTtype header2, UINT4 version2, UINT4 nsamples2 );
static BOOLEAN timestamp_in_list( LIGOTimeGPS timestamp, LIGOTimeGPSVector *list );
static long get_file_len ( FILE *fp );

static FILE * fopen_SFTLocator ( const struct tagSFTLocator *locator );
static BOOLEAN has_valid_v2_crc64 (FILE *fp );

static void read_one_sft_from_fp (  LALStatus *status, SFTtype **sft, REAL8 fMin, REAL8 fMax, FILE *fp );
static int read_sft_header_from_fp (FILE *fp, SFTtype  *header, UINT4 *version, UINT8 *crc64, BOOLEAN *swapEndian, CHAR **comment, UINT4 *numBins );
static int read_v2_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, UINT8 *header_crc64, UINT8 *ref_crc64, CHAR **comment, BOOLEAN swapEndian);
static int read_v1_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, BOOLEAN swapEndian);
static int compareSFTdesc(const void *ptr1, const void *ptr2);
static UINT8 calc_crc64(const CHAR *data, UINT4 length, UINT8 crc);
int read_SFTversion_from_fp ( UINT4 *version, BOOLEAN *need_swap, FILE *fp );

/*==================== FUNCTION DEFINITIONS ====================*/

/** Find the list of SFTs matching the \a file_pattern and satisfying the given \a constraints, 
 * return an \c SFTCatalog of the matching SFTs.
 * 
 * The optional \a constraints that can be specified are (type SFTConstraints)
 * - 'detector':      	which detector
 * - 'time-span':    	GPS start- and end-times
 * - 'timestamps':    	list of GPS start-times
 *
 * 
 * ==> The returned SFTCatalog can be used directly as input to LALLoadSFTs() 
 * to load a single-IFO SFTVector, or LALLoadMultiSFTs() to load a 
 * multi-IFO vector of SFTVectors 
 *
 * Except for the 'file_pattern' input, all the other constraints are optional
 * and can be passed as NULL (either globally constraings==NULL, or individually).
 *
 * Note that the constraints are combined by 'AND' and the resulting full constraint 
 * MUST be satisfied (in particular: if 'timestamps' is given, all timestamps within 
 * [startTime, endTime] MUST be found!.
 *
 * The returned SFTs in the catalogue are sorted by increasing GPS-epochs !
 * 
 */
void
LALSFTdataFind (LALStatus *status,
		SFTCatalog **catalog,		/**< [out] SFT-catalogue of matching SFTs */
		const CHAR *file_pattern,	/**< which SFT-files */
		SFTConstraints *constraints	/**< additional constraints for SFT-selection */
		)
{
  LALStringVector *fnames;
  UINT4 i, numFiles;
  SFTCatalog *ret = NULL;

  INITSTATUS (status, "LALSFTdataFind", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 

  /* ----- check input */
  ASSERT ( catalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( (*catalog) == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );

  ASSERT ( file_pattern, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );

  if ( constraints && constraints->detector )
    if ( strncmp(constraints->detector, "??", 2) && !is_valid_detector(constraints->detector) ) 
      {
	if ( lalDebugLevel ) LALPrintError( "\nInvalid detector-constraint '%s'\n\n", constraints->detector );
	ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );
      }

  /* prepare return-catalog */
  if ( (ret = LALCalloc ( 1, sizeof ( SFTCatalog ))) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  /* find matching filenames */
  if ( (fnames = find_files (file_pattern)) == NULL) {
    if ( lalDebugLevel ) LALPrintError ("\nFailed to get filelist for pattern '%s'.\n\n", file_pattern);
    ABORT (status, SFTFILEIO_EGLOB, SFTFILEIO_MSGEGLOB);
  }
  numFiles = fnames->length;


  /* ----- main loop: parse all matching files */
  for ( i = 0; i < numFiles; i ++ )
    {
      CHAR *fname = fnames->data[i];
      
      FILE *fp;      
      long file_len;

      /* merged SFTs need to satisfy stronger consistency-constraints (-> see spec) */
      BOOLEAN mfirst_block = TRUE;
      UINT4   mprev_version = 0;
      SFTtype mprev_header;
      REAL8   mprev_nsamples = 0;


      if ( ( fp = LALFopen( fname, "rb" ) ) == NULL )
	{
	  LALPrintError ( "\nFailed to open matched file '%s'\n\n", fname );
	  LALDestroyStringVector ( status->statusPtr, &fnames );
	  LALDestroySFTCatalog ( status->statusPtr, &ret );
	  ABORT( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	}
      if ( (file_len = get_file_len(fp)) == 0 )
	{
	  LALPrintError ( "\nGot file-len == 0 for '%s'\n\n", fname );
	  LALDestroyStringVector ( status->statusPtr, &fnames );
	  LALDestroySFTCatalog ( status->statusPtr, &ret );
	  fclose(fp);
	  ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	}
      /* go through SFT-blocks in fp */
      while ( ftell(fp) < file_len )
	{
	  SFTtype this_header;
	  UINT4 this_version;
	  UINT4 this_nsamples;
	  UINT8 this_crc;
	  CHAR *this_comment = NULL;
	  BOOLEAN endian;
	  BOOLEAN want_this_block = FALSE;

	  long this_filepos;
	  if ( (this_filepos = ftell(fp)) == -1 )
	    {
	      LALPrintError ("\nftell() failed for '%s'\n\n", fname );
	      LALDestroyStringVector ( status->statusPtr, &fnames );
	      LALDestroySFTCatalog ( status->statusPtr, &ret );
	      fclose (fp);
	      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	    }

	  if ( read_sft_header_from_fp (fp, &this_header, &this_version, &this_crc, 
					&endian, &this_comment, &this_nsamples ) != 0 )
	    {
	      if ( lalDebugLevel ) 
		LALPrintError ("\nERROR:File-block '%s:%ld' is not a valid SFT!\n\n", fname, ftell(fp) );
	      LALDestroyStringVector ( status->statusPtr, &fnames );
	      if ( this_comment ) LALFree ( this_comment );
	      LALDestroySFTCatalog ( status->statusPtr, &ret );
	      fclose(fp);
	      ABORT ( status, SFTFILEIO_EHEADER, SFTFILEIO_MSGEHEADER );
	    }
	  
	  /* if merged-SFT: check consistency constraints */
	  if ( !mfirst_block )
	    {
	      if ( ! consistent_mSFT_header ( mprev_header, mprev_version, mprev_nsamples, 
					      this_header, this_version, this_nsamples ) )
		{
		  if ( lalDebugLevel ) 
		    LALPrintError ("merged SFT-file '%s' contains inconsistent SFT-blocks!\n\n", fname);
		  if ( this_comment ) LALFree ( this_comment );
		  LALDestroyStringVector ( status->statusPtr, &fnames );
		  LALDestroySFTCatalog ( status->statusPtr, &ret );
		  fclose(fp);
		  ABORT ( status, SFTFILEIO_EMERGEDSFT, SFTFILEIO_MSGEMERGEDSFT );
		}
	    } /* if !mfirst_block */

	  mprev_header = this_header;
	  mprev_version = this_version;
	  mprev_nsamples = this_nsamples;

	  /* seek to end of SFT data-entries in file  */
	  if ( fseek ( fp, this_nsamples * 8 , SEEK_CUR ) == -1 )
	    {
	      if ( lalDebugLevel ) 
		LALPrintError ( "\nFailed to skip DATA field for SFT '%s': %s\n", fname, strerror(errno) );
	      if ( this_comment ) LALFree ( this_comment );
	      LALDestroyStringVector ( status->statusPtr, &fnames );
	      LALDestroySFTCatalog ( status->statusPtr, &ret );
	      fclose(fp);
	      ABORT ( status, SFTFILEIO_ESFTFORMAT, SFTFILEIO_MSGESFTFORMAT);
	    }

	  want_this_block = TRUE;	/* default */
	  /* but does this SFT-block satisfy the user-constraints ? */
	  if ( constraints )
	    {
	      if ( constraints->detector && strncmp(constraints->detector, "??", 2) ) 
		{
		  /* v1-SFTs have '??' as detector-name */
		  if ( ! strncmp (this_header.name, "??", 2 ) )
		    strncpy ( this_header.name, constraints->detector, 2 );	/* SET to constraint! */
		  else if ( strncmp( constraints->detector, this_header.name, 2) )
		    want_this_block = FALSE;
		}

	      if ( constraints->startTime && ( GPS2REAL8(this_header.epoch) < GPS2REAL8( *constraints->startTime)))
		want_this_block = FALSE;

	      if ( constraints->endTime && ( GPS2REAL8(this_header.epoch) > GPS2REAL8( *constraints->endTime ) ) )
		want_this_block = FALSE;

	      if ( constraints->timestamps && !timestamp_in_list(this_header.epoch, constraints->timestamps) )
		want_this_block = FALSE;

	    } /* if constraints */

	  if ( want_this_block )
	    {
	      UINT4 oldlen = ret->length;
	      SFTDescriptor *desc;
	      
	      ret->data = LALRealloc ( ret->data, (oldlen+1) * sizeof( *(ret->data) ) );

	      ret->length ++;
	      if ( ret->data == NULL )
		{
		  LALDestroyStringVector ( status->statusPtr, &fnames );
		  LALDestroySFTCatalog ( status->statusPtr, &ret );
		  if ( this_comment ) LALFree ( this_comment );
		  fclose(fp);
		  ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
		}
	      desc = &(ret->data[oldlen]);

	      desc->locator = LALCalloc ( 1, sizeof ( *(desc->locator) ) );
	      if ( desc->locator )
		desc->locator->fname = LALCalloc( 1, strlen(fname) + 1 );
	      if ( (desc->locator == NULL) || (desc->locator->fname == NULL ) )
		{
		  LALDestroyStringVector ( status->statusPtr, &fnames );
		  LALDestroySFTCatalog ( status->statusPtr, &ret );
		  fclose(fp);
		  ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
		}
	      strcpy ( desc->locator->fname, fname );
	      desc->locator->offset = this_filepos;

	      desc->header  = this_header;
	      desc->comment = this_comment;
	      desc->numBins = this_nsamples;
	      desc->version = this_version;
	      desc->crc64   = this_crc;

	    } /* if want_this_block */
	  else
	    {
	      if ( this_comment ) LALFree ( this_comment );
	    }

	  mfirst_block = FALSE;

	} /* while !feof */

      fclose(fp);

    } /* for i < numFiles */

  /* free matched filenames */
  LALDestroyStringVector ( status->statusPtr, &fnames );

  /* ----- final consistency-checks: ----- */

  /* did we find exactly the timestamps that lie within [startTime, endTime]? */
  if ( constraints && constraints->timestamps )
    {
      LIGOTimeGPSVector *ts = constraints->timestamps;
      UINT4 numRequested = 0;
      REAL8 t0, t1;
      if ( constraints->startTime )
	t0 = GPS2REAL8 ( (*constraints->startTime) );
      else
	t0 = -1;
      if ( constraints->endTime )
	t1 = GPS2REAL8 ( (*constraints->endTime) );
      else
	t1 = LAL_REAL4_MAX;	/* large enough */
      
      for (i=0; i < ts->length; i ++ )
	{
	  REAL8 ti = GPS2REAL8(ts->data[i]);
	  if ( (t0 <= ti) && ( ti <= t1 ) )
	    numRequested ++;
	}
      if ( numRequested != ret->length )
	{
	  if ( lalDebugLevel ) 
	    LALPrintError ("\nERROR: found %s SFTs (%d) than given timestamps within [%f, %f] (%d)\n\n",
			   (ret->length < numRequested )?"fewer":"more", 
			   ret->length, t0, t1, numRequested );
	  ABORT ( status, SFTFILEIO_ECONSTRAINTS, SFTFILEIO_MSGECONSTRAINTS );
	}
    } /* if constraints->timestamps */

  /* have all matched SFTs identical dFreq values ? */
  for ( i = 0; i < ret->length; i ++ )
    {
      SFTtype first_header = empty_SFTtype;	/* avoid warning */
      SFTtype this_header = ret->data[i].header;

      if ( i == 0 )
	first_header = this_header;

      /* dont give out v1-SFTs without detector-entry, except if constraint->detector="??" ! */
      if ( !constraints || !constraints->detector || strncmp(constraints->detector, "??", 2) )
	{
	  if ( !strncmp ( this_header.name, "??", 2 ) )
	    {
	      LALDestroySFTCatalog ( status->statusPtr, &ret );
	      if ( lalDebugLevel ) 
		LALPrintError ("\nERROR: '%s' matched v1-SFTs but no detector-constraint given!\n\n",
			       file_pattern);
	      ABORT ( status, SFTFILEIO_EDETECTOR, SFTFILEIO_MSGEDETECTOR );
	    }
	} /* if detector-constraint was not '??' */

      if ( this_header.deltaF != first_header.deltaF ) 
	{
	  LALDestroySFTCatalog ( status->statusPtr, &ret );
	  LALPrintError("\nERROR: file-pattern '%s' matched SFTs with inconsistent Tsft!\n\n",
			file_pattern);
	  ABORT ( status, SFTFILEIO_EDIFFTSFT, SFTFILEIO_MSGEDIFFTSFT );	  
	}

    } /* for i < numSFTs */


  /* sort catalog in order of increasing GPS-time */
  qsort( (void*)ret->data, ret->length, sizeof( ret->data[0] ), compareSFTdesc );


  /* return result catalog (=sft-vect and locator-vect) */
  (*catalog) = ret;
  
  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALSFTdataFind() */



/** Extract a timstamps-vector from the given SFTCatalog.
 */
void
LALSFTtimestampsFromCatalog (LALStatus *status,
			     LIGOTimeGPSVector **timestamps,	/**< [out] extracted timestamps */
			     const SFTCatalog *catalog )	/**< input SFT-catalogue */
{
  UINT4 numSFTs;
  UINT4 i;
  LIGOTimeGPSVector *ret = NULL;

  INITSTATUS (status, "LALSFTdataFind", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 

  ASSERT ( timestamps, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( catalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( catalog->length > 0, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( *timestamps == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );

  numSFTs = catalog->length;

  TRY ( LALCreateTimestampVector ( status->statusPtr, &ret, numSFTs ), status );

  for ( i=0; i < numSFTs; i ++ )
    ret->data[i] = catalog->data[i].header.epoch;


  /* done: return Ts-vector */
  (*timestamps) = ret;

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALTimestampsFromSFTCatalog() */



/** Load the given frequency-band <tt>[fMin, fMax]</tt> (inclusively) from the SFT-files listed in the 
 * SFT-'catalogue' ( returned by LALSFTdataFind() ).
 *
 * Note: \a fMin (or \a fMax) is allowed to be set to \c -1, which means to read in all
 * Frequency-bins from the lowest (or up to the highest) found in the SFT-file.
 * 
 * Note 2: The returned frequency-interval is guaranteed to contain <tt>[fMin, fMax]</tt>, 
 * but is allowed to be larger, as it must be an interval of discrete frequency-bins as found 
 * in the SFT-file. 
 */
void
LALLoadSFTs ( LALStatus *status,
	      SFTVector **sfts,			/**< [out] vector of read-in SFTs */
	      const SFTCatalog *catalog,	/**< The 'catalogue' of SFTs to load */
	      REAL8 fMin,		  /**< minumum requested frequency (-1 = read from lowest) */
	      REAL8 fMax		  /**< maximum requested frequency (-1 = read up to highest) */
	      )
{
  UINT4 i;
  UINT4 numSFTs;
  SFTVector *ret = NULL;
  /* CHAR prev_det[3] = {0,0,0}; */  /* equal-detector constraint dropped */

  INITSTATUS (status, "LALLoadSFTs", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 

  ASSERT ( sfts, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( *sfts == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );
  
  ASSERT ( catalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );

  numSFTs = catalog -> length;

  /* Allocate SFT-vector with zero length to append to */
  if ( ( ret = LALCalloc ( 1, sizeof( *ret ) )) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  for ( i=0; i < numSFTs; i++ )
    {
      FILE *fp;
      SFTtype *this_sft = NULL;


      /* equal-detector constraint dropped */
      if ( (fp = fopen_SFTLocator ( catalog->data[i].locator )) == NULL )
	{
	  LALDestroySFTVector (status->statusPtr, &ret );
	  ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
	}

      read_one_sft_from_fp ( status->statusPtr,  &this_sft, fMin, fMax, fp );
      BEGINFAIL ( status ) {
	LALDestroySFTVector (status->statusPtr, &ret );
	fclose(fp);
      } ENDFAIL(status);

      fclose(fp);

      /* NOTE: LALSFTdataFind() requires overriding v1 detector-name ('??')
       * by a proper detector-name using the constraint-entry.
       * This is field therefore updated here from the catalog */
      strncpy ( this_sft->name, catalog->data[i].header.name, 2 );

      LALAppendSFT2Vector ( status->statusPtr, ret, this_sft );
      BEGINFAIL ( status ) {
	LALDestroySFTtype ( status->statusPtr, &this_sft );
	LALDestroySFTVector (status->statusPtr, &ret );
	fclose(fp);
      } ENDFAIL(status);

      LALDestroySFTtype ( status->statusPtr, &this_sft );

    } /* for i < numSFTs */

  /* return final SFT-vector */
  (*sfts) = ret;
  
  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALLoadSFTs() */


/** Function to load a catalog of SFTs from possibly different detectors.
    This is similar to LALLoadSFTs except that the input SFT catalog is 
    allowed to contain multiple ifos.  The output is the structure 
    MultiSFTVector which is a vector of (pointers to) SFTVectors, one for 
    each ifo found in the catalog.   As in LALLoadSFTs, fMin and fMax can be 
    set to -1 to get the full SFT from the lowest to the highest frequency 
    bin found in the SFT.  
 */
void LALLoadMultiSFTs ( LALStatus *status,
			MultiSFTVector **out,             /**< [out] vector of read-in SFTs -- one sft vector for each ifo found in catalog*/
			const SFTCatalog *inputCatalog,   /**< The 'catalogue' of SFTs to load */
			REAL8 fMin,		          /**< minumum requested frequency (-1 = read from lowest) */
			REAL8 fMax		          /**< maximum requested frequency (-1 = read up to highest) */
			)
     
{
  UINT4 k, j, i, length;
  UINT4 numifo=0; /* number of ifos */
  UINT4 numifoMax, numifoMaxNew; /* for memory allocation purposes */
  CHAR  *name=NULL; 
  CHAR  **ifolist=NULL; /* list of ifo names */
  UINT4  *numsfts=NULL; /* number of sfts for each ifo */
  UINT4 **sftLocationInCatalog=NULL; /* location of sfts in catalog for each ifo */
  SFTCatalog **catalog=NULL;
  MultiSFTVector *multSFTVec=NULL;

  INITSTATUS (status, "LALLoadMultiSFTs", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 

  ASSERT ( out, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( *out == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );  
  ASSERT ( inputCatalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( inputCatalog->length, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );

  length = inputCatalog->length; 
  if ( (name = (CHAR *)LALCalloc(3, sizeof(CHAR))) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  /* the number of ifos can be at most equal to length */
  /* each ifo name is 2 characters + \0 */

  numifoMax = 3; /* should be sufficient -- realloc used later in case required */

  if ( (ifolist = (CHAR **)LALCalloc( 1, numifoMax * sizeof(CHAR *))) == NULL) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }
  if ( (sftLocationInCatalog = (UINT4 **)LALCalloc( 1, numifoMax * sizeof(UINT4 *))) == NULL) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }
  if ( (numsfts = (UINT4 *)LALCalloc( 1, numifoMax * sizeof(UINT4))) == NULL) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  for ( k = 0; k < numifoMax; k++) {
    if ( (ifolist[k] = (CHAR *)LALCalloc( 1, 3*sizeof(CHAR))) == NULL) {
      ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
    }
    if ( (sftLocationInCatalog[k] = (UINT4 *)LALCalloc( 1, length*sizeof(UINT4))) == NULL) {
      ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
    }
  }

  /* loop over sfts in catalog and look at ifo names and
     find number of different ifos and number of sfts for each ifo
     Also find location of sft in catalog  */
  for ( k = 0; k < length; k++)
    { 
      strncpy( name, inputCatalog->data[k].header.name, 3 );

      /* go through list of ifos till a match is found or list is exhausted */      
      for ( j = 0; ( j < numifo ) && strncmp( name, ifolist[j], 3); j++ ) 
	;

      if ( j < numifo ) 
	{
	  /* match found with jth existing ifo */
	  sftLocationInCatalog[j][ numsfts[j] ] = k;
	  numsfts[j]++;
	}      
      else
	{
	  /* add ifo to list of ifos */

	  /* first check if number of ifos is larger than numifomax */
	  /* and realloc if necessary */
	  if ( numifo >= numifoMax )
	    {
	      numifoMaxNew = numifoMax + 3;
	      if ( (ifolist = (CHAR **)LALRealloc( ifolist, numifoMaxNew * sizeof(CHAR *))) == NULL) {
		ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
	      }
	      if ( (sftLocationInCatalog = (UINT4 **)LALRealloc( sftLocationInCatalog, numifoMaxNew * sizeof(UINT4 *))) == NULL) {
		ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
	      }
	      if ( (numsfts = (UINT4 *)LALRealloc( numsfts, numifoMaxNew * sizeof(UINT4))) == NULL) {
		ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
	      }
	      
	      for ( i = numifoMax; i < numifoMaxNew; i++) {
		if ( (ifolist[i] = (CHAR *)LALCalloc( 1, 3*sizeof(CHAR))) == NULL) {
		  ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
		}
		if ( (sftLocationInCatalog[i] = (UINT4 *)LALCalloc( 1, length*sizeof(UINT4))) == NULL) {
		  ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
		}
	      } /* loop from numifoMax to numifoMaxNew */

	      /* reset numifoMax */	      
	      numifoMax = numifoMaxNew;
	    } /* if ( numifo >= numifoMax) -- end of realloc */
	  
	  strncpy( ifolist[numifo], name, 3);
	  sftLocationInCatalog[j][0] = k;
	  numsfts[numifo] = 1;
	  numifo++;

	} /* else part of if ( j < numifo ) */ 
    } /*  for ( k = 0; k < length; k++) */

  /* now we can create the catalogs */
  if ( (catalog = (SFTCatalog **)LALCalloc( numifo, sizeof(SFTCatalog *))) == NULL) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  for ( j = 0; j < numifo; j++)
    {
      if ( (catalog[j] = (SFTCatalog *)LALCalloc(1, sizeof(SFTCatalog))) == NULL) {
	ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
      }
      catalog[j]->length = numsfts[j];
      if ( (catalog[j]->data = (SFTDescriptor *)LALCalloc( numsfts[j], sizeof(SFTDescriptor))) == NULL) {
	ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
      }

      for ( k = 0; k < numsfts[j]; k++)
	{
	  UINT4 location = sftLocationInCatalog[j][k]; 
	  catalog[j]->data[k] = inputCatalog->data[location];
	}
    }

  /* create multi sft vector */
  if ( (multSFTVec = (MultiSFTVector *)LALCalloc(1, sizeof(MultiSFTVector))) == NULL){
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }
  multSFTVec->length = numifo;

  if ( (multSFTVec->data = (SFTVector **)LALCalloc(numifo, sizeof(SFTVector *))) == NULL) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }
  for ( j = 0; j < numifo; j++) {
    LALLoadSFTs ( status->statusPtr, multSFTVec->data + j, catalog[j], fMin, fMax );    
    BEGINFAIL ( status ) {
      /* free sft vectors created previously in loop */ 
      for ( i = 0; (INT4)i < (INT4)j-1; i++)
	LALDestroySFTVector ( status->statusPtr, multSFTVec->data + i);
      LALFree(multSFTVec->data);
      LALFree(multSFTVec);

      /* also free catalog and other memory allocated earlier */
      for ( i = 0; i < numifo; i++) {
	LALFree(catalog[i]->data);
	LALFree(catalog[i]);
      }
      LALFree( catalog);
      
      for ( i = 0; i < numifoMax; i++) {
	LALFree(ifolist[i]);
	LALFree(sftLocationInCatalog[i]);
      }
      LALFree(ifolist);
      LALFree(sftLocationInCatalog);
      
      LALFree(numsfts);
      LALFree(name);
      
    } ENDFAIL(status);
  }



  /* free memory and exit */
  for ( j = 0; j < numifo; j++) {
    LALFree(catalog[j]->data);
    LALFree(catalog[j]);
  }
  LALFree( catalog);

  for ( k = 0; k < numifoMax; k++) {
    LALFree(ifolist[k]);
    LALFree(sftLocationInCatalog[k]);
  }
  LALFree(ifolist);
  LALFree(sftLocationInCatalog);

  LALFree(numsfts);
  LALFree(name);

  *out = multSFTVec;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALLoadMultiSFTs() */


/** Function to check validity of SFTs listed in catalog. 
 * This function simply reads in those SFTs and checks their CRC64 checksum, which
 * is the only check that has not yet been done by the operations up to this point.
 * 
 * Returns the LAL-return code of a failure (or 0 on success) in 'check_result'.
 *
 * \note: because this function has to read the complete SFT-data into memory for the
 *  whole set of matching SFTs, it is potentially slow and memory-intensive.
 *
 * This function will NOT fail if one of the SFT-operations fails, instead it returns
 * the status-code of this failure in 'check_result'. The function will fail, however,
 * on invalid input.
 *
 */
void
LALCheckSFTs ( LALStatus *status, 
	       INT4 *check_result, 		/**< LAL-status of SFT-operations */
	       const CHAR *file_pattern,     /**< where to find the SFTs: normally a path+file-pattern */
	       SFTConstraints *constraints	/**< additional constraints for SFT-selection */
	       )
{
  LALStatus sft_status = empty_status;
  SFTCatalog *catalog = NULL;
  UINT4 i;

  INITSTATUS (status, "LALCheckSFTs", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 

  ASSERT ( check_result, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( file_pattern, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );

  if ( constraints && constraints->detector && ! is_valid_detector(constraints->detector) ) 
    {
      if ( lalDebugLevel ) LALPrintError( "\nInvalid detector-constraint '%s'\n\n", constraints->detector );
      ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );
    }

  /* Step 1: find the catalog of matching SFTs */
  LALSFTdataFind ( &sft_status, &catalog, file_pattern, constraints );
  if ( ((*check_result) = sft_status.statusCode) != 0 ) 
    goto sft_failed;

  /* Step 2: step through SFTs and check CRC64 */
  for ( i=0; i < catalog->length; i ++ )
    {
      FILE *fp;

      switch ( catalog->data[i].version  )
	{
	case 1:	/* version 1 had no CRC  */
	  continue;
	case 2:
	  if ( (fp = fopen_SFTLocator ( catalog->data[i].locator )) == NULL )
	    {
	      LALPrintError ( "Failed to open locator '%s'\n", 
			      XLALshowSFTLocator ( catalog->data[i].locator ) );
	      (*check_result) = SFTFILEIO_EFILE;
	      goto sft_failed;
	    }
	  if ( ! has_valid_v2_crc64 ( fp ) != 0 )
	    {
	      LALPrintError ( "CRC64 checksum failure for SFT '%s'\n", 
			      XLALshowSFTLocator ( catalog->data[i].locator ) );
	      (*check_result) = SFTFILEIO_ECRC64;
	      goto sft_failed;
	    }
	  break;

	default:
	  LALPrintError ( "Illegal SFT-version encountered : %d\n", catalog->data[i].version );
	  (*check_result) = SFTFILEIO_EVERSION;
	  goto sft_failed;
	  break;
	} /* switch (version ) */

    } /* for i < numSFTs */

 sft_failed:

  if ( catalog ) {
    TRY ( LALDestroySFTCatalog ( status->statusPtr, &catalog ), status );
  }

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* LALCheckSFTs() */


/** Read timestamps file and returns timestamps vector (alloc'ed in here!).
 */
void
LALReadTimestampsFile (LALStatus* status, LIGOTimeGPSVector **timestamps, const CHAR *fname)
{  
  FILE *fp;
  LIGOTimeGPSVector *ts = NULL;

  INITSTATUS( status, "LALReadTimestampsFile", SFTFILEIOC );
  ATTATCHSTATUSPTR (status);

  ASSERT (fname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT (timestamps, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (*timestamps == NULL, status, SFTFILEIO_ENONULL,  SFTFILEIO_MSGENONULL);

  if ( (fp = fopen( fname, "r")) == NULL) {
    LALPrintError("\nUnable to open timestampsname file %s\n\n", fname);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* initialize empty timestamps-vector*/
  if ( (ts = LALCalloc(1, sizeof(LIGOTimeGPSVector))) == NULL) {
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
  }

  while(1)
    {
      INT4 secs, ns;
      if (fscanf ( fp, "%d  %d\n", &secs, &ns ) != 2)
	break;

      if ( ( secs < 0 ) || ( ns < 0 ) ) {
	LALPrintError ("\nERROR: timestamps-file contained negative time-entry in line %d \n\n",
		       ts->length);
	ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
      }

      /* make space for the new entry */
      ts->length ++;
      if ( (ts->data = LALRealloc(ts->data, ts->length * sizeof(ts->data[0])) ) == NULL) {
	ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
      }
      
      ts->data[ts->length - 1].gpsSeconds = (UINT4) secs;
      ts->data[ts->length - 1].gpsNanoSeconds = (UINT4) ns;

    } /* while entries found */
  fclose(fp);

  /* hand over timestamps vector */
  (*timestamps) = ts;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadTimestampsFile() */


/** Write the given *v2-normalized* (i.e. dt x DFT) SFTtype to a v2-SFT file. 
 *  Add the comment to SFT if comment != NULL.
 *
 * NOTE: Currently this only supports writing v2-SFTs.
 * If you need to write a v1-SFT, you should use LALWrite_v2SFT_to_v1file() 
 */
void
LALWriteSFT2file (LALStatus *status,
		  const SFTtype *sft,		/**< SFT to write to disk */
		  const CHAR *fname,		/**< filename */
		  const CHAR *comment)		/**< optional comment (for v2 only) */
{
  FILE  *fp = NULL;
  UINT4 comment_len = 0;
  UINT4 pad_len = 0;
  CHAR pad[] = {0, 0, 0, 0, 0, 0, 0};	/* for comment-padding */
  _SFT_header_v2_t rawheader;

  INITSTATUS (status, "LALWriteSFTfile", SFTFILEIOC);
  ATTATCHSTATUSPTR (status);   
 
  /*   Make sure the arguments are not NULL and perform basic checks*/ 
  ASSERT (sft,   status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (sft->data,  status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (sft->deltaF > 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (sft->f0 >= 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT ( (sft->epoch.gpsSeconds >= 0) && (sft->epoch.gpsNanoSeconds >= 0), status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT ( sft->epoch.gpsNanoSeconds < 1000000000, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT ( sft->data->length > 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);

  ASSERT (fname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL); 

  if ( !is_valid_detector(sft->name) ) {
    if ( lalDebugLevel ) LALPrintError ("\nInvalid detector prefix '%c%c'\n\n", sft->name[0], sft->name[1] );
    ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );
  }

  /* open SFT-file for writing */
  if ( (fp = LALFopen ( fname, "wb" )) == NULL )
    {
      if ( lalDebugLevel ) LALPrintError ("\nFailed to open file '%s' for writing: %s\n\n", fname, strerror(errno));
      ABORT ( status, SFTFILEIO_ESFTWRITE, SFTFILEIO_MSGESFTWRITE );      
    }
  
  /* comment length including null terminator to string must be an
   * integer multiple of eight bytes. comment==NULL means 'no
   * comment'  
   */
  if (comment) 
    {
      comment_len = strlen(comment) + 1;
      pad_len = (8 - (comment_len % 8)) % 8;
    }
  
  /* ----- fill out header */
  rawheader.version        		= 2;
  rawheader.gps_sec        		= sft->epoch.gpsSeconds;
  rawheader.gps_nsec       		= sft->epoch.gpsNanoSeconds;
  rawheader.tbase          		= 1.0 / sft->deltaF;
  rawheader.first_frequency_index 	= MYROUND( sft->f0 / sft->deltaF );
  rawheader.nsamples       		= sft->data->length;
  rawheader.crc64          		= 0;	/* set to 0 for crc-calculation */
  rawheader.detector[0]    		= sft->name[0];
  rawheader.detector[1]    		= sft->name[1];
  rawheader.padding[0]     		= 0;
  rawheader.padding[1]     		= 0;
  rawheader.comment_length 		= comment_len + pad_len;
  
  /* ----- compute CRC */
  rawheader.crc64 = calc_crc64((const CHAR*)&rawheader, sizeof(rawheader), ~(0ULL));

  if ( comment )
    {
      rawheader.crc64 = calc_crc64((const CHAR*)comment, comment_len, rawheader.crc64);
      rawheader.crc64 = calc_crc64((const CHAR*)pad, pad_len, rawheader.crc64);
    }
  rawheader.crc64 = calc_crc64((const CHAR*) sft->data->data, sft->data->length * sizeof( *sft->data->data ), rawheader.crc64);
  
  /* ----- write the header to file */
  if (1 != fwrite( &rawheader, sizeof(rawheader), 1, fp) ) {
    ABORT ( status, SFTFILEIO_ESFTWRITE, SFTFILEIO_MSGESFTWRITE );
  }
  
  /* ----- write the comment to file */
  if ( comment )
    {
      if ( comment_len != fwrite( comment, 1, comment_len, fp) ) {
	ABORT ( status, SFTFILEIO_ESFTWRITE, SFTFILEIO_MSGESFTWRITE );
      }
      if (pad_len != fwrite( pad, 1, pad_len, fp) ) {
	ABORT ( status, SFTFILEIO_ESFTWRITE, SFTFILEIO_MSGESFTWRITE );
      }
    } /* if comment */

  /* write the data to the file.  Data must be packed REAL,IMAG,REAL,IMAG,... */
  if ( sft->data->length != fwrite( sft->data->data, sizeof(*sft->data->data), sft->data->length, fp) ) {
    ABORT ( status, SFTFILEIO_ESFTWRITE, SFTFILEIO_MSGESFTWRITE );
  }
  
  fclose(fp);

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* WriteSFTtoFile() */



/** Write the given *v2-normalized* (i.e. dt x DFT) SFTVector to a directory. 
 *  Add the comment to SFT if comment != NULL.
 *
 * NOTE: Currently this only supports writing v2-SFTs.
 * If you need to write a v1-SFT, you should use LALWriteSFTfile() 
 *
 * Output SFTs have naming convention following LIGO-T040164-01
 */
void
LALWriteSFTVector2Dir (LALStatus *status,
		       const SFTVector *sftVect,	/**< SFT vector to write to disk */
		       const CHAR *dirname,		/**< base filename (including directory path)*/
		       const CHAR *comment,		/**< optional comment (for v2 only) */
		       const CHAR *description)         /**< optional sft description to go in the filename */
{
  UINT4 length, k;
  CHAR *filename = NULL;
  CHAR filenumber[16];
  SFTtype *sft;
  UINT4 timeBase, duration;
  UINT4 filenamelen;
  LIGOTimeGPS time0;

  INITSTATUS (status, "LALWriteSFTVector2Dir", SFTFILEIOC);
  ATTATCHSTATUSPTR (status);   

  ASSERT (sftVect, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (sftVect->data, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (sftVect->length > 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (dirname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);

  length = sftVect->length;

  filenamelen = 128 + strlen(dirname);
  if ( description ) 
    filenamelen += strlen ( description );

  if ( (filename = (CHAR *)LALCalloc(1, filenamelen )) == NULL) {
    ABORT( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
  }

  /* will not be same as actual sft timebase if it is not
     an integer number of seconds */
  timeBase = ceil(1.0/sftVect->data[0].deltaF);

  for ( k = 0; k < length; k++) {

    sft = sftVect->data + k;
    if ( sft == NULL ) {
      ABORT( status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
    }

    if ( sft->name == NULL ) {
      ABORT( status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
    }


    time0 = sft->epoch;

    /* calculate sft duration -- may be different from timebase if nanosecond  
       of sft-epoch is non-zero or is timebase is not an integer number of seconds */ 
    XLALGPSAdd( &time0, timeBase);
    if ( time0.gpsNanoSeconds > 0) {
      duration = time0.gpsSeconds - sft->epoch.gpsSeconds + 1;
    }
    else {
      duration = time0.gpsSeconds - sft->epoch.gpsSeconds;
    }
    
    /* create the k^th filename following naming convention 
       -- try to simplify this*/
    strcpy( filename, dirname);
    strcat( filename, "/");
    strncat( filename, sft->name, 1);
    strcat( filename, "-1_"); /* single (not merged) sft */
    strncat( filename, sft->name, 2); /* full detector name */
    strcat( filename, "_");
    sprintf( filenumber, "%d", timeBase); /* sft timebase */
    strcat( filename, filenumber);
    strcat( filename, "SFT");
    if ( description ) {
      strcat( filename, "_");
      strcat( filename, description);
    }
    strcat( filename, "-");
    sprintf( filenumber, "%09d", sft->epoch.gpsSeconds);
    strncat( filename, filenumber, 9);
    strcat( filename, "-");
    sprintf( filenumber, "%d", duration);
    strcat( filename, filenumber);
    strcat( filename, ".sft");
    
    /* write the k^th sft */
    TRY ( LALWriteSFT2file ( status->statusPtr, sft, filename, comment), status);
  }

  LALFree(filename);

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* WriteSFTVector2Dir() */





/** For backwards-compatibility: write a *v2-normalized* (ie dt x DFT) SFTtype
 * to a v1-SFT file.
 *
 * NOTE: the only difference to WriteSFTfile() is that the data-normalization
 *      is changed back to v1-type 'DFT', by dividing the dt corresponding to the
 *      frequency-band contained in the SFTtype.
 */
void
LALWrite_v2SFT_to_v1file (LALStatus *status,
			  const SFTtype *sft,		/**< SFT to write to disk */
			  const CHAR *fname)		/**< filename */
{
  UINT4 i, numBins;
  REAL8 Band, dt;
  SFTtype v1SFT;

  INITSTATUS (status, "LALWriteSFTfile", SFTFILEIOC);
  ATTATCHSTATUSPTR (status);   
 
  /*   Make sure the arguments are not NULL and perform basic checks*/ 
  ASSERT (sft,   status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (sft->data,  status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (sft->deltaF > 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (sft->f0 >= 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT ( (sft->epoch.gpsSeconds >= 0) && (sft->epoch.gpsNanoSeconds >= 0), status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT ( sft->epoch.gpsNanoSeconds < 1000000000, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (sft->data->length > 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);

  ASSERT (fname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL); 

  if ( !is_valid_detector(sft->name) ) {
    ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );
  }

  numBins = sft->data->length;
  Band = sft->deltaF * numBins ;
  dt = 1.0 / (2.0 * Band);

  v1SFT.data = NULL;
  TRY ( LALCopySFT (status->statusPtr, &v1SFT, sft ), status );

  for ( i=0; i < numBins; i ++ )
    {
      v1SFT.data->data[i].re = (REAL4) ( (REAL8)v1SFT.data->data[i].re / dt );
      v1SFT.data->data[i].im = (REAL4) ( (REAL8)v1SFT.data->data[i].im / dt );
    }

  TRY ( LALWriteSFTfile (status->statusPtr, &v1SFT, fname ), status );

  XLALDestroyCOMPLEX8Vector ( v1SFT.data );

  DETATCHSTATUSPTR ( status );
  RETURN ( status );

} /* LALWrite_v2SFT_to_v1file() */




/** [OBSOLETE] Write a *v1-normalized* (i.e. raw DFT) SFTtype to a SFT-v1 file. 
 * 
 * \note:only SFT-spec v1.0 is supported, and the SFTtype must follow the 
 * *obsolete* v1-normalization. => Use LALWriteSFT2file() to write v2 SFTs !
 *
 */
void
LALWriteSFTfile (LALStatus  *status,
		 const SFTtype *sft,		/**< SFT to write to disk */
		 const CHAR *outfname)		/**< filename */
{
  FILE  *fp = NULL;
  COMPLEX8  *inData;
  INT4  i;
  UINT4 datalen;
  REAL4  *rawdata;
  CHAR *rawheader, *ptr;
  SFTHeader header;

  INITSTATUS (status, "LALWriteSFTfile", SFTFILEIOC);
  ATTATCHSTATUSPTR (status);   
 
  /*   Make sure the arguments are not NULL and perform basic checks*/ 
  ASSERT (sft,   status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (sft->data,  status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (sft->deltaF > 0, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);
  ASSERT (outfname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL); 

  /* fill in the header information */
  header.version = 1.0;
  header.gpsSeconds = sft->epoch.gpsSeconds;
  header.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  header.timeBase = 1.0 / sft->deltaF;
  header.fminBinIndex = (INT4) floor (sft->f0 / sft->deltaF + 0.5);	/* round to closest int! */
  header.length = sft->data->length; 

  /* build raw header for writing to disk */
  rawheader = LALCalloc (1, sizeof(_SFT_header_v1_t) );
  if (rawheader == NULL) {
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);    
  }
  ptr = rawheader;
  *(REAL8*) ptr = header.version;
  ptr += sizeof (REAL8);
  *(INT4*) ptr = header.gpsSeconds;
  ptr += sizeof (INT4);
  *(INT4*) ptr = header.gpsNanoSeconds;
  ptr += sizeof (INT4);
  *(REAL8*) ptr = header.timeBase;
  ptr += sizeof (REAL8);
  *(INT4*) ptr = header.fminBinIndex;
  ptr += sizeof (INT4);
  *(INT4*) ptr = header.length; 

  /* write data into a contiguous REAL4-array */
  datalen = 2 * header.length * sizeof(REAL4);	/* amount of bytes for SFT-data */

  rawdata = LALCalloc (1, datalen);
  if (rawdata == NULL) {
    LALFree (rawheader);
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);    
  }

  inData = sft->data->data;
  for ( i = 0; i < header.length; i++)
    {
      rawdata[2 * i]     = inData[i].re;
      rawdata[2 * i + 1] = inData[i].im;
    } /* for i < length */


  /* open the file for writing */
  fp = fopen(outfname, "wb");
  if (fp == NULL) {
    LALFree (rawheader);
    LALFree (rawdata);
    LALPrintError ("\nFailed to open file '%s' for writing!\n\n", outfname );
    ABORT (status, SFTFILEIO_EFILE,  SFTFILEIO_MSGEFILE);
  }

  /* write the header*/
  if( fwrite( rawheader, sizeof(_SFT_header_v1_t), 1, fp) != 1) {
    LALFree (rawheader);
    LALFree (rawdata);
    fclose (fp);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }
  
  /* write the data */
  if (fwrite( rawdata, datalen, 1, fp) != 1) {
    LALFree (rawheader);
    LALFree (rawdata);
    fclose (fp);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* done */
  fclose(fp);
  LALFree (rawheader);
  LALFree (rawdata);


  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* WriteSFTtoFile() */



/** Free an 'SFT-catalogue' */
void
LALDestroySFTCatalog ( LALStatus *status,
		       SFTCatalog **catalog )	/**< the 'catalogue' to free */
{
  INITSTATUS (status, "LALDestroySFTcatalog", SFTFILEIOC);

  ASSERT ( catalog, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );

  if ( *catalog )
    {
      if ( (*catalog) -> data )
	{
	  UINT4 i;
	  for ( i=0; i < (*catalog)->length; i ++ )
	    {
	      SFTDescriptor *ptr = &( (*catalog)->data[i] );
	      if ( ptr->locator ) 
		{
		  if ( ptr->locator->fname )
		    LALFree ( ptr->locator->fname );
		  LALFree ( ptr->locator );
		}
	      if ( ptr->comment ) 
		LALFree ( ptr->comment );
	      
	      /* this should not happen, but just in case: free data-entry in SFT-header */
	      if ( ptr->header.data )
		LALFree ( ptr->header.data );
	    } /* for i < length */
	
	  LALFree ( (*catalog)->data );

	} /* if *catalog->data */
      
      LALFree ( *catalog );

    } /* if *catalog */

  (*catalog) = NULL;

  RETURN ( status );

} /* LALDestroySFTCatalog() */



/** Mostly for *debugging* purposes: provide a user-API to allow inspecting the SFT-locator
 * [which is an OPAQUE entry in the SFTCatalog!]
 *
 * NOTE: this returns a STATIC string, so don't try to FREE it, and make a copy if you need
 * to keep it beyond one call of this function!
 * 
 */
const CHAR *
XLALshowSFTLocator ( const struct tagSFTLocator *locator )
{
  static CHAR ret[512];

  if ( !locator )
    return NULL;

  LALSnprintf ( ret, sizeof(ret), "%s : %ld", locator->fname, locator->offset );
  ret[ sizeof(ret) - 1 ] = 0;
  
  return ret;

} /* XLALshowSFTLocator() */


INT4 XLALCountIFOsInCatalog( const SFTCatalog *catalog)
{
  
  UINT4 k, j, numifo=0, length;
  CHAR  *name=NULL; 
  CHAR  **ifolist=NULL; /* list of ifo names */

  length = catalog->length;

  name = (CHAR *)LALCalloc(3, sizeof(CHAR));

  ifolist = (CHAR **)LALCalloc( length, sizeof(CHAR *));
  for ( k = 0; k < length; k++) 
    ifolist[k] = (CHAR *)LALCalloc( 3, sizeof(CHAR));

  /* go through catalog and look at each ifo name */  
  for ( k = 0; k < length; k++)
    { 
      strncpy( name, catalog->data[k].header.name, 3 );

      /* go through list of ifos till a match is found or list is exhausted */      
      for ( j = 0; ( j < numifo ) && strncmp( name, ifolist[j], 3); j++ ) 
	;

      if ( j >= numifo ) 
	{
	  /* add ifo to list of ifos */
	  strncpy( ifolist[numifo], name, 3);
	  numifo++;
	}
   
    }

  LALFree(name);
  for ( j = 0; j < catalog->length; j++) 
    LALFree(ifolist[j]);
  LALFree(ifolist);
  
  return numifo;

}

/*================================================================================
 * OBSOLETE and deprecated SFT-v1 API :
 *================================================================================*/

/** [DEPRECATED] Function to read and return a list of SFT-headers for given 
 * filepattern and start/end times.
 *
 * The \a startTime and \a endTime entries are allowed to be NULL, 
 * in which case they are ignored.
 *
 * \note We return the headers as an SFTVector, but with empty data-fields.
 */
void
LALGetSFTheaders (LALStatus *status,
		  SFTVector **headers,		/**< [out] Vector of SFT-headers */
		  const CHAR *fpattern,		/**< path/filepattern */
		  const LIGOTimeGPS *startTime,	/**< include only SFTs after this time (can be NULL) */
		  const LIGOTimeGPS *endTime)	/**< include only SFTs before this (can be NULL)*/
{
  UINT4 i, numSFTs;
  SFTVector *out = NULL;
  LALStringVector *fnames;
  REAL8 t0, t1;		/* start- and end-times as reals */
  UINT4 numHeaders;

  INITSTATUS (status, "LALGetSFTheaders", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 

  /* check input */
  ASSERT ( headers, status,  SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT ( *headers == NULL, status,  SFTFILEIO_ENONULL,  SFTFILEIO_MSGENONULL);
  ASSERT ( fpattern, status,  SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);

  if ( startTime )
    t0 = GPS2REAL8( *startTime );
  else
    t0 = 0;
  
  if ( endTime )
    t1 = GPS2REAL8 ( *endTime );
  else
    t1 = 0;

  /* get filelist of files matching fpattern */
  if ( (fnames = find_files (fpattern)) == NULL) {
    ABORT (status, SFTFILEIO_EGLOB, SFTFILEIO_MSGEGLOB);
  }

  /* prepare output-vector */
  if ( (out = LALCalloc (1, sizeof(*out) )) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  numSFTs = fnames->length;

  /* main loop: load all SFT-headers and put them into an SFTvector */
  numHeaders = 0;
  for (i=0; i < numSFTs; i++)
    {
      SFTHeader header;
      REAL8 epoch;
      SFTtype *sft;

      TRY ( LALReadSFTheader (status->statusPtr, &header, fnames->data[i]), status );
      epoch = 1.0 * header.gpsSeconds + 1.0e-9 * header.gpsNanoSeconds;

      if ( t0 && ( epoch < t0 ) )	/* epoch earlier than start-time ==> skip */
	continue;
      if ( t1 && ( epoch > t1 ) )	/* epoch later than end-time ==> skip */
	continue;

      /* found suitable SFT ==> add to SFTVector of headers */
      numHeaders ++;
      if ( ( out->data = LALRealloc(out->data, numHeaders * sizeof( *out->data) )) == NULL ) 
	{
	  LALFree (out);
	  ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
	}
      
      sft = &(out->data[numHeaders - 1]);

      /* fill in header-data */
      strncpy (sft->name, fnames->data[i], LALNameLength);
      sft->name[LALNameLength - 1 ] = '\0';	/* make sure it's 0-terminated */
      sft->deltaF  			= 1.0 / header.timeBase;
      sft->f0      			= header.fminBinIndex / header.timeBase;
      sft->epoch.gpsSeconds     	= header.gpsSeconds;
      sft->epoch.gpsNanoSeconds 	= header.gpsNanoSeconds;

      sft->data = NULL;	/* no SFT-data proper */

    } /* for i < numSFTs */

  out->length = numHeaders;

  DestroyStringVector (fnames);

  *headers = out;
  
  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALGetSFTheaders() */		  



/** [DEPRECATED] Basic SFT reading-function.
 *  Given a filename \a fname and frequency-limits [\a fMin, \a fMax], 
 *  returns an SFTtype \a sft containing the SFT-data.
 * 
 * \note 1) the actual returned frequency-band is 
 * <tt>[floor(Tsft * fMin), ceil(Tsft * fMax)] / Tsft</tt>,
 * i.e. we need to round to an integer frequency-bin within the SFT, but
 * the requested frequency-band is guaranteed to be contained in the output 
 * (if present in the SFT-file), but can be slightly larger.
 * 
 * 2) The special input <tt>fMin=fMax=0</tt> means to read and
 * return the <em>whole</em> frequency-band contained in the SFT-file.
 *
 * 3) Currently only SFTv1 are supported!!
 *
 */
void
LALReadSFTfile (LALStatus *status, 
		SFTtype **sft, 		/**< [out] output SFT */
		REAL8 fMin, 		/**< lower frequency-limit */
		REAL8 fMax,		/**< upper frequency-limit */
		const CHAR *fname)	/**< path+filename */
{ 
  SFTHeader  header;		/* SFT file-header version1 */
  REAL8 deltaF;
  UINT4 readlen;
  INT4 fminBinIndex, fmaxBinIndex;
  SFTtype *outputSFT = NULL;
  UINT4 i;
  REAL4 renorm;

  INITSTATUS (status, "LALReadSFTfile", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 
  
  ASSERT (sft, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (*sft == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL);
  ASSERT (fname,  status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (fMin <= fMax, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);

  /* read the header */
  TRY ( LALReadSFTheader (status->statusPtr, &header, fname), status);

  /* ----- figure out which data we want to read ----- */
  deltaF = 1.0 / header.timeBase;

  /* special case: fMin==fMax==0 means "read all" */
  if ( (fMin == 0) && (fMax == 0) )
    {
      fminBinIndex = header.fminBinIndex;
      fmaxBinIndex = fminBinIndex + header.length - 1;
    }
  else
    {
      /* find the right frequency-bin and number of bins
       * The rounding here is chosen such that the required 
       * frequency-interval is _guaranteed_ to lie within the 
       * returned range  */
      fminBinIndex = (INT4) floor (fMin * header.timeBase);  /* round this down */
      fmaxBinIndex = (INT4) ceil  (fMax * header.timeBase);  /* round up */
    }

  readlen = (UINT4)(fmaxBinIndex - fminBinIndex) + 1;	/* number of bins to read */

  /* allocate the final SFT to be returned */
  TRY ( LALCreateSFTtype (status->statusPtr, &outputSFT, readlen), status);


  /* and read it, using the lower-level function: */
  LALReadSFTdata (status->statusPtr, outputSFT, fname, fminBinIndex);
  BEGINFAIL (status) {
    LALDestroySFTtype (status->statusPtr, &outputSFT);
  } ENDFAIL (status);


  /*
   * NOTE: the following renormalization is necessary for v1-SFTs
   * as the data are not multiplied by dt, therefore extracting a 
   * sub-band B' of the total band B requires mulitplication of 
   * the data by B'/B.
   * SFTv2 will store FFT-data multiplied by dt and then this will
   * not be necessary any more.
   * 
   */
  renorm = 1.0 * readlen / header.length;

  /* let's re-normalize and fill data into output-vector */
  if (renorm != 1)
    for (i=0; i < readlen; i++)
      {
	outputSFT->data->data[i].re *= renorm;
	outputSFT->data->data[i].im *= renorm;
      }    
  
  /* that's it: return */
  *sft = outputSFT;

  DETATCHSTATUSPTR (status);
  RETURN(status);

} /* LALReadSFTfile() */



/** [DEPRECATED] Higher-level SFT-reading function to read a whole vector of SFT files 
 * and return an SFTvector \a sftvect. The handling of
 * [fMin, fMax] is identical to LALReadSFTfile().
 * 
 * \note 1) the file-pattern \a fpattern can use a wide range of glob-patterns.
 *
 * 2) currently the SFTs matching the pattern are required to have the same 
 * number of frequency bins, otherwise an error will be returned. 
 * (This might be relaxed in the future).
 *
 * 3) This function does <em>not</em> use <tt>glob()</tt> and should therefore
 * be safe even under condor.
 *
 */
void
LALReadSFTfiles (LALStatus *status,
		 SFTVector **sftvect,	/**< [out] output SFT vector */
		 REAL8 fMin,	       	/**< lower frequency-limit */
		 REAL8 fMax,		/**< upper frequency-limit */
		 UINT4 wingBins,	/**< number of frequency-bins to be added left and right. */
		 const CHAR *fpattern)	/**< path/filepattern */
{

  UINT4 i, numSFTs;
  SFTVector *out = NULL;
  SFTtype *oneSFT = NULL;
  SFTHeader header;
  REAL8 dFreq; 		/* frequency spacing in SFT */
  REAL8 fWing;		/* frequency band to be added as "wings" to the 'physical band' */
  LALStringVector *fnames;
  UINT4 firstlen = 0;

  INITSTATUS (status, "LALReadSFTfiles", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 
  
  ASSERT (sftvect, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (*sftvect == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL);
  ASSERT (fpattern,  status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (fMin <= fMax, status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL);

  /* make filelist 
   * NOTE: we don't use glob() as it was reported to fail under condor */
  if ( (fnames = find_files (fpattern)) == NULL) {
    ABORT (status, SFTFILEIO_EGLOB, SFTFILEIO_MSGEGLOB);
  }

  /* allocate head of output SFT-vector */
  if ( (out = LALCalloc ( 1, sizeof(SFTVector) )) == NULL ) {
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  /* read header of first sft to determine Tsft, and therefore dfreq */
  LALReadSFTheader(status->statusPtr, &header, fnames->data[0]);
  BEGINFAIL(status) {
    DestroyStringVector (fnames);    
  } ENDFAIL(status);

  numSFTs = fnames->length;
  dFreq = 1.0 / header.timeBase;
  fWing = wingBins * dFreq;

  /* main loop: load all SFTs and put them into the SFTvector */
  for (i=0; i < numSFTs; i++)
    {
      LALReadSFTfile (status->statusPtr, &oneSFT, fMin-fWing, fMax+fWing, fnames->data[i]);
      BEGINFAIL (status) {
	DestroyStringVector (fnames);    
	LALDestroySFTtype (status->statusPtr, &oneSFT);
	if (out) LALDestroySFTVector (status->statusPtr, &out);
      } ENDFAIL (status);

      if ( !firstlen ) 
	firstlen = oneSFT->data->length;
      /* make sure all SFTs have same length */
      if ( oneSFT->data->length != firstlen )
	{
	  DestroyStringVector (fnames);    
	  LALDestroySFTtype (status->statusPtr, &oneSFT);
	  LALDestroySFTVector (status->statusPtr, &out);
	  ABORT (status, SFTFILEIO_EDIFFLENGTH, SFTFILEIO_MSGEDIFFLENGTH);
	} /* if length(thisSFT) != common length */

      LALAppendSFT2Vector ( status->statusPtr, out, oneSFT );
      BEGINFAIL(status) {
	DestroyStringVector (fnames);    
	LALDestroySFTtype (status->statusPtr, &oneSFT);
	LALDestroySFTVector (status->statusPtr, &out);
      } ENDFAIL(status);

      LALDestroySFTtype (status->statusPtr, &oneSFT);
      oneSFT = NULL;	/* important for next call of LALReadSFTfile()! */

    } /* for i < numSFTs */

  DestroyStringVector (fnames);

  *sftvect = out;

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadSFTfiles () */



/*================================================================================
 * LOW-level internal SFT-handling functions, should *NOT* be used outside this file!
 *================================================================================*/

/** Open an "SFT" defined by the SFT-locator, return a FILE-pointer to the beginning of this SFT.
 *
 * \note The returned filepointer could point to an SFT-block within a merged SFT-file, 
 * so you should not assume that SEEK_SET takes you to the beginning of this block!
 * (instead you'd have to save the current position returned by this function, which \
 * points to the beginning of the block!)
 *
 * NOTE: Ideally this should be the *ONLY* function using the internal structure of the opaque
 * SFTLocator type
 *
 */
FILE *
fopen_SFTLocator ( const struct tagSFTLocator *locator )
{
  FILE *fp = NULL;
  CHAR *fname;

  if ( !locator )
    return NULL;

  fname = locator->fname;
  if ( (fp = LALFopen( fname, "rb" )) == NULL )
    {
      LALPrintError ("\nFailed to open SFT '%s' for reading: %s\n\n", fname, strerror(errno) );
      return NULL;
    }
      
  if ( fseek( fp, locator->offset, SEEK_SET ) == -1 )
    {
      LALPrintError ("\nFailed to set fp-offset to '%ld': %s\n\n", locator->offset, strerror(errno) );
      fclose(fp);
      return NULL;
    }

  return fp;

} /* fopen_SFTLocator() */


/** Check the v2 SFT-block starting at fp for valid crc64 checksum.
 *  Restores filepointer before leaving.
 */
BOOLEAN
has_valid_v2_crc64 ( FILE *fp )
{
  long save_filepos;
  UINT8 computed_crc, ref_crc;
  SFTtype header;
  UINT4 numBins;
  CHAR *comment = NULL;
  UINT4 data_len;
  char block[BLOCKSIZE];
  UINT4 version;
  BOOLEAN need_swap;

  /* input consistency */
  if ( !fp ) 
    {
      LALPrintError ("\nhas_valid_v2_crc64() was called with NULL filepointer!\n\n");
      return FALSE;
    }
  
  /* store fileposition for restoring in case of failure */
  if ( ( save_filepos = ftell(fp) ) == -1 )
    {
      LALPrintError ("\nERROR: ftell() failed: %s\n\n", strerror(errno) );
      return -1;
    }

  if ( read_SFTversion_from_fp ( &version, &need_swap, fp ) != 0 )
    return -1;

  if ( version != 2 ) 
    {
      LALPrintError ("\nhas_valid_v2_crc64() was called on non-v2 SFT.\n\n");
      return -1;
    }

  /* ----- compute CRC ----- */
  /* read the header, unswapped, only to obtain it's crc64 checksum */
  if ( read_v2_header_from_fp ( fp, &header, &numBins, &computed_crc, &ref_crc, &comment, need_swap ) != 0 )
    {
      if ( comment ) LALFree ( comment );
      return FALSE;
    }
  if ( comment ) LALFree ( comment );

  /* read data in blocks of BLOCKSIZE, computing CRC */
  data_len = numBins * 8 ;	/* COMPLEX8 data */
  while ( data_len > 0 )
    {
      /* read either BLOCKSIZE or amount remaining */
      int toread = (BLOCKSIZE < data_len) ? BLOCKSIZE : data_len;
      if (toread != (int)fread ( block, 1, toread, fp) ) 
	{
	  LALPrintError ("\nFailed to read all frequency-bins from SFT.\n\n");
	  return FALSE;
	}
      data_len -= toread;

      /* compute CRC64: don't endian-swap for that! */
      computed_crc = calc_crc64( (const CHAR*)block, toread, computed_crc );

    } /* while data */
     
  /* check that checksum is consistent */
  return ( computed_crc == ref_crc );
  
} /* has_valid_v2_crc64 */


/** [DEPRECATED]: Low-level function to read only the SFT-header of a given file. 
 * 
 * NOTE: don't use! This function is obsolete and SFT-v1-specific and only kept for 
 * backwards-compatibility with Hough-codes. 
 */
void 
LALReadSFTheader (LALStatus  *status,
		  SFTHeader   *header,	/**< [out] returned header */
		  const CHAR  *fname)	/**< path+filename */
{
  FILE *fp = NULL;
  SFTHeader  header1;
  CHAR *rawheader = NULL;
  CHAR *ptr = NULL;
  CHAR inVersion[8];
  REAL8 version;
  BOOLEAN swapEndian = 0;

  INITSTATUS (status, "LALReadSFTHeader", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (header, status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  ASSERT (fname,  status, SFTFILEIO_ENULL,  SFTFILEIO_MSGENULL);
  
  /* opening the SFT binary file */
  fp = LALOpenDataFile( fname );
  if (fp == NULL) {
    ABORT (status, SFTFILEIO_EFILE,  SFTFILEIO_MSGEFILE);
  }
  
  /* read version-number */
  if  (fread (inVersion, sizeof(inVersion), 1, fp) != 1) {
    fclose (fp);
    if (lalDebugLevel) LALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* try version 1.0 */
  version = 1.0;
  /* try swapping the version if it is not equal */
  if(memcmp(inVersion,&version,sizeof(version))){
    endian_swap (inVersion, sizeof(inVersion),1);
    swapEndian = 1;
    /* fail if still not true */
    if(memcmp(inVersion,&version,sizeof(version))){
      fclose (fp);
      if (lalDebugLevel) LALPrintError ("\nOnly v1-SFTs supported at the moment!: %s\n\n", fname);
      ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
    }
  }

  /* read the whole header */
  rawheader = LALCalloc (1, sizeof(_SFT_header_v1_t) );
  if (rawheader == NULL) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);    
  }
  
  rewind (fp);	/* go back to start */
  if (fread( rawheader, sizeof(_SFT_header_v1_t), 1, fp) != 1) {
    fclose (fp);
    LALFree (rawheader);
    if (lalDebugLevel) LALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  fclose(fp);

  /* now fill-in the header-struct with the appropriate fields */
  /* NOTE: we have to do it this way, because the struct can have 
   * padding in memory, so the fields are not guaranteed to lie 'close'
   * Endian-swapping ist done here if necessary
   */
  ptr = rawheader;
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(REAL8),1);
  header1.version	= *(REAL8*) ptr;
  ptr += sizeof(REAL8);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  header1.gpsSeconds 	= *(INT4*) ptr;
  ptr += sizeof(INT4);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  header1.gpsNanoSeconds= *(INT4*) ptr;
  ptr += sizeof(INT4);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(REAL8),1);
  header1.timeBase      = *(REAL8*) ptr;
  ptr += sizeof(REAL8);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  header1.fminBinIndex  = *(INT4*) ptr;
  ptr += sizeof(INT4);
  if (swapEndian) endian_swap((CHAR*)ptr,sizeof(INT4),1);
  header1.length     	= *(INT4*) ptr;

  LALFree (rawheader);

  /* ----- do some consistency-checks on the header-fields: ----- */

  /* gps_sec and gps_nsec >= 0 */
  if ( (header1.gpsSeconds < 0) || (header1.gpsNanoSeconds <0) ) {
    if (lalDebugLevel) LALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* tbase > 0 */
  if ( header1.timeBase <= 0 ) {
    if (lalDebugLevel) LALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* fminindex >= 0 */
  if (header1.fminBinIndex < 0) {
    if (lalDebugLevel) LALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }
  
  /* nsamples >= 0 */
  if (header1.length < 0) {
    if (lalDebugLevel) LALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* ok, the SFT-header seems consistent, so let's return it */
  *header = header1;
  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadSFTheader() */



/** [DEPRECATED] This is a function for low-level SFT data-reading:
 * the SFT-data is read starting from fminBinIndex and filled
 * into the pre-allocate vector sft of length N
 *
 * NOTE: !! NO re-normalization is done here!! this remains up
 *       to the caller of this function!!
 *
 */
void
LALReadSFTdata(LALStatus *status,
	       SFTtype    *sft,    /**< [out] output-SFT: assuming memory is allocated  */
	       const CHAR *fname,  /**< path+filename */
	       INT4 fminBinIndex)  /**< minimun frequency-index to read */
{ 
  FILE        *fp = NULL;
  SFTHeader  header;
  UINT4 offset, readlen;
  REAL4 *rawdata = NULL;
  CHAR inVersion[8];
  REAL8 version;
  BOOLEAN swapEndian = 0;
  UINT4 i;

  INITSTATUS (status, "LALReadSFTdata", SFTFILEIOC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (sft->data, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  ASSERT (fname, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL);
  
  /* Read header */
  TRY ( LALReadSFTheader (status->statusPtr, &header, fname), status);

  /* check that the required frequency-interval is part of the SFT */
  readlen = sft->data->length;
  if ( (fminBinIndex < header.fminBinIndex) 
       || (fminBinIndex + (INT4)readlen > header.fminBinIndex + header.length) ) {
    ABORT (status, SFTFILEIO_EFREQBAND, SFTFILEIO_MSGEFREQBAND);
  }

  /* how many frequency-bins to skip */
  offset = fminBinIndex - header.fminBinIndex;

  /* open file for reading */
  if ( (fp = LALOpenDataFile( fname )) == NULL) {
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* read version-number */
  if  (fread (inVersion, sizeof(inVersion), 1, fp) != 1) {
    fclose (fp);
    if (lalDebugLevel) LALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* set invalid version */
  version = -1.0;

  if (version < 0) {
    /* try version 1.0 */
    version = 1.0;
    /* try swapping the version if it is not equal */
    if(memcmp(inVersion,&version,sizeof(version))){
      endian_swap (inVersion, sizeof(inVersion),1);
      swapEndian = 1;
      /* set invalid version if still not true */
      if(memcmp(inVersion,&version,sizeof(version))){
	version = -1;
	endian_swap (inVersion, sizeof(inVersion),1);
	swapEndian = 0;
      }
    }
  }

  /* fail if the version is invalid */
  if (version < 0) {
    fclose (fp);
    if (lalDebugLevel) LALPrintError ("\nInvalid SFT-file: %s\n\n", fname);
    ABORT (status, SFTFILEIO_EHEADER,  SFTFILEIO_MSGEHEADER);
  }

  /* check compatibility of version with this function */
  if (version != 1) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EVERSION, SFTFILEIO_MSGEVERSION);
  }

  /* skip SFT-header in file */
  rewind (fp);
  if (fseek(fp, sizeof(_SFT_header_v1_t), SEEK_SET) != 0) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* skip offset data points to the correct frequency-bin */
  if (fseek(fp, offset * 2 * sizeof(REAL4), SEEK_CUR) != 0) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }

  /* ----- prepare memory for data-reading ----- */
  rawdata = LALCalloc (1, 2 * readlen *sizeof(REAL4) );
  if (rawdata == NULL) {
    fclose (fp);
    ABORT (status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM);
  }

  /* we don't rely on memory-packing, so we read into a REAL4 array first */
  if (fread( rawdata, 2 * readlen * sizeof(REAL4), 1, fp) != 1) {
    fclose(fp);
    LALFree (rawdata);
    ABORT (status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE);
  }
  fclose(fp);

  /* now fill data into output-vector */
  for (i=0; i < readlen; i++)
    {
      if (swapEndian)
	endian_swap((CHAR*)&rawdata[2*i], sizeof(REAL4), 2);
      sft->data->data[i].re = rawdata[2 * i];
      sft->data->data[i].im = rawdata[2 * i + 1];
    }    
  
  LALFree (rawdata);

  /* now fill in the header-info */
  strncpy (sft->name, fname, LALNameLength);
  sft->name[LALNameLength - 1 ] = '\0';	/* make sure it's 0-terminated */
  sft->deltaF  			= 1.0 / header.timeBase;
  sft->f0      			= fminBinIndex / header.timeBase;
  sft->epoch.gpsSeconds     	= header.gpsSeconds;
  sft->epoch.gpsNanoSeconds 	= header.gpsNanoSeconds;

  
  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadSFTdata() */



/***********************************************************************
 * internal helper functions
 ***********************************************************************/

static BOOLEAN
timestamp_in_list( LIGOTimeGPS timestamp, LIGOTimeGPSVector *list )
{
  UINT4 i;
  LIGOTimeGPS *el;

  if ( !list )
    return FALSE;

  el = &(list->data[0]);
  for ( i=0; i < list->length; i++, el++)
    {
      if ( (timestamp.gpsSeconds == el->gpsSeconds) && ( timestamp.gpsNanoSeconds == el->gpsNanoSeconds ) )
	return TRUE;
    } /* for i < length */

  return FALSE;

} /* timestamp_in_list() */


/* check consistency constraints for SFT-blocks within a merged SFT-file,
 * see SFT-v2 spec */
static BOOLEAN
consistent_mSFT_header ( SFTtype header1, UINT4 version1, UINT4 nsamples1, SFTtype header2, UINT4 version2, UINT4 nsamples2 )
{
  /* 1) identical detector */
  if ( (header1.name[0] != header2.name[0]) || (header1.name[1] != header2.name[1]) )
    {
      if ( lalDebugLevel ) LALPrintError ("\nInvalid merged SFT: non-identical detectors\n\n");
      return FALSE;
    }

  /* 2) identical version-number */
  if ( version1 != version2 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nInvalid merged SFT: non-identical version-numbers\n\n");
      return FALSE;
    }

  /* 3) increasing GPS-times */
  if ( GPS2REAL8 ( header1.epoch ) >= GPS2REAL8 ( header2.epoch ) )
    {
      if ( lalDebugLevel ) LALPrintError ("\nInvalid merged SFT: non-increasing GPS epochs \n\n" );
      return FALSE;
    }

  /* 4) identical tbase */
  if ( header1.deltaF != header2.deltaF )
    {
      if ( lalDebugLevel ) LALPrintError ("\nInvalid merged SFT: non-identical time baselines\n\n");
      return FALSE;
    }

  /* 5) identical start-frequency */
  if ( header1.f0 != header2.f0 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nInvalid merged SFT: non-identical start-frequencies\n\n");
      return FALSE;
    }

  /* 6) identical number of frequency-bins */
  if ( nsamples1 != nsamples2 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nInvalid merged SFT: non-identical number of frequency-bins\n\n" );
      return FALSE;
    }

  return TRUE;

} /* consistent_mSFT_header() */

/** Read one SFT, leave filepointer at the end of the read SFT if successful, 
 * leave fp at initial position if failure 
 *
 * If fMin == -1: read starting from first bin
 *    fMax == -1: read until last bin
 *
 *
 * This uses read_sft_header_from_fp() for reading the header, and then reads the data. 
 * 
 * NOTE: we DO NOT check the crc64 checksum in here, a separate "check-SFT" function
 * should be used for that.
 *
 * NOTE2:  The returned SFT is always normalized correctly according to SFT-v2 spec
 * (see LIGO-T040164-01-Z, and LIGO-T010095-00), irrespective of the input-files.
 *
 */
static void
read_one_sft_from_fp (  LALStatus *status, SFTtype **sft, REAL8 fMin, REAL8 fMax , FILE *fp )
{
  SFTtype *ret = NULL;
  UINT4 version;
  UINT8 crc64;
  BOOLEAN swapEndian;
  UINT4 firstBin2read, lastBin2read, numBins2read;
  UINT4 firstSFTbin, lastSFTbin, numSFTbins;
  INT4 offsetBins;
  long offsetBytes;
  volatile REAL8 tmp;	/* intermediate results: try to force IEEE-arithmetic */

  INITSTATUS (status, "read_one_sft_from_fp", SFTFILEIOC);
  ATTATCHSTATUSPTR ( status );

  ASSERT ( sft, status, SFTFILEIO_ENULL, SFTFILEIO_MSGENULL );
  ASSERT ( *sft == NULL, status, SFTFILEIO_ENONULL, SFTFILEIO_MSGENONULL );

  TRY ( LALCreateSFTtype ( status->statusPtr, &ret, 0 ), status );

  if ( read_sft_header_from_fp (fp, ret, &version, &crc64, &swapEndian, NULL, &numSFTbins ) != 0 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nFailed to read SFT-header!\n\n");
      LALDestroySFTtype ( status->statusPtr, &ret );
      ABORT ( status, SFTFILEIO_EHEADER, SFTFILEIO_MSGEHEADER ); 
    }

  tmp = ret->f0 / ret->deltaF;
  firstSFTbin = MYROUND ( tmp );
  lastSFTbin = firstSFTbin + numSFTbins - 1;

  /* figure out frequency-bounds to read in */
  if ( fMin > 0 )
    {
      tmp = fMin / ret->deltaF;
      firstBin2read = floor ( tmp );	/* round down! requested fMin is contained */
    }
  else
    {
      firstBin2read = firstSFTbin;
    }

  if ( fMax > 0 )
    {
      tmp = fMax / ret->deltaF;
      lastBin2read = ceil ( tmp );	/* round up! requested fMax is contained */
    }
  else
    {
      lastBin2read = lastSFTbin;
    }

  if ( firstBin2read > lastBin2read )
    {
      if ( lalDebugLevel ) 
	LALPrintError ("\nEmpty frequency-interval requested [%d, %d] bins\n\n", 
		       firstBin2read, lastBin2read );
      ABORT ( status, SFTFILEIO_EVAL, SFTFILEIO_MSGEVAL );
    }

  /* check that requested interval is found in SFT */
  if ( firstBin2read < firstSFTbin )
    {
      if ( lalDebugLevel ) 
	LALPrintError ( "\nRequested fMin=%f is not contained in SFT! (%d < %d)\n\n", 
			fMin, firstBin2read, firstSFTbin );
      LALDestroySFTtype ( status->statusPtr, &ret );
      ABORT ( status, SFTFILEIO_EFREQBAND, SFTFILEIO_MSGEFREQBAND );
    }
  if ( lastBin2read > lastSFTbin )
    {
      if ( lalDebugLevel ) 
	LALPrintError ( "\nRequested fMax=%f is not contained in SFT! (%d > %d)\n\n", 
			fMax, lastBin2read, lastSFTbin );
      LALDestroySFTtype ( status->statusPtr, &ret );
      ABORT ( status, SFTFILEIO_EFREQBAND, SFTFILEIO_MSGEFREQBAND );
    }
  
  offsetBins = firstBin2read - firstSFTbin;
  offsetBytes = offsetBins * 2 * sizeof( REAL4 );
  numBins2read = lastBin2read - firstBin2read + 1;

  if ( fseek ( fp, offsetBytes, SEEK_CUR ) != 0 )
    {
      if ( lalDebugLevel ) 
	LALPrintError ( "\nFailed to fseek() to first frequency-bin %d: %s\n\n", 
			firstBin2read, strerror(errno) );
      LALDestroySFTtype ( status->statusPtr, &ret );
      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
    }

  if ( (ret->data = XLALCreateCOMPLEX8Vector ( numBins2read )) == NULL ) {
    LALDestroySFTtype ( status->statusPtr, &ret );
    ABORT ( status, SFTFILEIO_EMEM, SFTFILEIO_MSGEMEM );
  }

  if ( numBins2read != fread ( ret->data->data, 2*sizeof( REAL4 ), numBins2read, fp ) )
    {
      if (lalDebugLevel) LALPrintError ("\nFailed to read %d bins from SFT!\n\n", numBins2read );
      LALDestroySFTtype ( status->statusPtr, &ret );
      ABORT ( status, SFTFILEIO_EFILE, SFTFILEIO_MSGEFILE );
    }

  /* update the start-frequency entry in the SFT-header to the new value */
  ret->f0 = 1.0 * firstBin2read * ret->deltaF;

  /* take care of normalization and endian-swapping */
  if ( version == 1 || swapEndian )
    {
      UINT4 i;
      REAL8 band = 1.0 * numSFTbins * ret->deltaF;/* need the TOTAL frequency-band in the SFT-file! */
      REAL8 fsamp = 2.0 * band;
      REAL8 dt = 1.0 / fsamp;

      for ( i=0; i < numBins2read; i ++ )
	{
	  REAL4 *rep, *imp;

	  rep = &(ret->data->data[i].re);
	  imp = &(ret->data->data[i].im);

	  if ( swapEndian )
	    {
	      endian_swap( (CHAR *) rep, sizeof ( *rep ), 1 );
	      endian_swap( (CHAR *) imp, sizeof ( *imp ), 1 );
	    }

	  /* if the SFT-file was in v1-Format: need to renormalize the data now by 'Delta t' 
	   * in order to follow the correct SFT-normalization 
	   * (see LIGO-T040164-01-Z, and LIGO-T010095-00)
	   */ 
	  if ( version == 1 )
	    {
	      (*rep) *= dt;
	      (*imp) *= dt;
	    }
	} /* for i < numBins2read */
    } /* if SFT-v1 */

  /* return resulting SFT */
  (*sft) = ret;

  DETATCHSTATUSPTR ( status );
  RETURN (status);
} /* read_one_sft_from_fp() */


/* Try to read an SFT-header (of ANY VALID SFT-VERSION) at the given FILE-pointer fp,
 * and return the SFT-header, SFT-version-number and number of frequency-samples in the SFT. 
 *
 * Sets the filepointer fp at the end of the header if successful, leaves it at 
 * initial position if not. 
 * 
 * RETURN 0 = OK, -1 on ERROR
 *
 * We do basic checking of compliance with the SFT-spec (<tt>LIGO-T04164-01-Z</tt>)
 * as far as a single header is concerned.
 *
 * NOTE: fatal errors will produce a LALPrintError() error-message, but
 * non-fatal 'SFT-format'-errors will only output error-messages if lalDebugLevel > 0.
 * --> this function can therefore be used to check if a given file actually contains SFTs
 * 
 *
 */
static int
read_sft_header_from_fp (FILE *fp, SFTtype  *header, UINT4 *version, UINT8 *crc64, BOOLEAN *swapEndian, CHAR **comment, UINT4 *numBins )
{
  SFTtype head = empty_SFTtype;
  UINT4 nsamples;
  CHAR *comm = NULL;
  UINT8 ref_crc = 0;
  UINT8 header_crc;

  UINT4 ver;
  BOOLEAN need_swap;
  long save_filepos;

  if ( !header || !version || !numBins || !fp  )
    {
      LALPrintError ("\nERROR read_sft_header_from_fp(): called with NULL input\n\n");
      return -1;
    }
  if ( comment && ((*comment) != NULL) )
    {
      LALPrintError ("\nERROR: Comment-string passed to read_sft_header_from_fp() is not empty!\n\n");
      return -1;
    }

  /* store fileposition for restoring in case of failure */
  if ( ( save_filepos = ftell(fp) ) == -1 )
    {
      LALPrintError ("\nftell() failed: %s\n\n", strerror(errno) );
      return -1;
    }

  if ( read_SFTversion_from_fp ( &ver, &need_swap, fp ) != 0 ) 
    return -1;


  /* read this SFT-header with version-specific code */
  head = empty_SFTtype;

  switch( ver )
    {
    case 1:
      if ( read_v1_header_from_fp ( fp, &head, &nsamples, need_swap ) != 0 )
	goto failed;
      break;
      
    case 2:
      if ( read_v2_header_from_fp ( fp, &head, &nsamples, &header_crc, &ref_crc, &comm, need_swap ) != 0 )
	goto failed;
      break;
      
    default:
      LALPrintError ("\nUnsupported SFT-version %d.\n\n", ver);
      goto failed;
      break;
    } /* switch(ver) */


  /* ----- some general SFT-header consistency-checks */
  if ( (head.epoch.gpsSeconds < 0 ) || ( head.epoch.gpsNanoSeconds < 0 ) || ( head.epoch.gpsNanoSeconds >= 1000000000)  )
    {
      if ( lalDebugLevel ) LALPrintError ("\nInvalid GPS-epoch in SFT : [%d, %d]!\n\n", 
					  head.epoch.gpsSeconds, head.epoch.gpsNanoSeconds );
      goto failed;
    }

  if ( head.deltaF <= 0 ) 
    {
      if ( lalDebugLevel ) LALPrintError ("\nNegative frequency-spacing in SFT!\n\n");
      goto failed;
    }
      
  if ( head.f0 < 0 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nNegative start-frequency in SFT!\n\n");
      goto failed;
    }
      
  /* ok */
  (*header) = head;
  (*version) = ver;

  if ( comment )	  /* return of comment is optional */
    (*comment) = comm;
  else
    if ( comm ) LALFree(comm);

  (*swapEndian) = need_swap;
  (*crc64) = ref_crc;
  (*numBins) = nsamples;
  return 0;

  /* ---------- */
 failed:
  /* restore filepointer initial position  */
  if ( fseek ( fp, save_filepos, SEEK_SET ) == -1 )
    LALPrintError ("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno) );
  
  /* free comment  if we allocated one */
  if ( comm )
    LALFree ( comm );
  
  return -1;

} /* read_sft_header_from_fp() */



/* ----- SFT v2 -specific header-reading function: 
 *
 * return general SFTtype header, place filepointer at the end of the header if it succeeds,
 * set fp to initial position if it fails.
 * RETURN: 0 = OK, -1 = ERROR 
 *
 * NOTE: fatal errors will produce a LALPrintError() error-message, but
 * non-fatal 'SFT-format'-errors will only output error-messages if lalDebugLevel > 0.
 * 
 */
static int
read_v2_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, UINT8 *header_crc64, UINT8 *ref_crc64, CHAR **comment, BOOLEAN swapEndian)
{
  _SFT_header_v2_t rawheader;
  long save_filepos;
  CHAR *comm = NULL;
  UINT8 crc;


  /* check input-consistency */
  if ( !fp || !header || !nsamples || !comment )
    {
      LALPrintError ( "\nERROR read_v2_header_from_fp(): called with NULL input!\n\n");
      return -1;
    }
  if ( comment && (*comment != NULL) )
    {
      LALPrintError ("\nERROR: Comment-string passed to read_v2_header_from_fp() is not NULL!\n\n");
      return -1;
    }
  
  /* store fileposition for restoring in case of failure */
  if ( ( save_filepos = ftell(fp) ) == -1 )
    {
      LALPrintError ("\nERROR: ftell() failed: %s\n\n", strerror(errno) );
      return -1;
    }

  /* read the whole header */
  if (fread( &rawheader, sizeof(rawheader), 1, fp) != 1) 
    {
      if (lalDebugLevel) LALPrintError ("\nCould not read v2-header. %s\n\n", strerror(errno) );
      goto failed;
    }

  /* ----- compute CRC for the header:
   * NOTE: the CRC checksum is computed on the *bytes*, not the numbers, 
   * so this must be computed before any endian-swapping.
   */
  {
    UINT8 save_crc = rawheader.crc64;
    rawheader.crc64 = 0;

    crc = calc_crc64((const CHAR*)&rawheader, sizeof(rawheader), ~(0ULL));

    rawheader.crc64 = save_crc;
    /* NOTE: we're not done with crc yet, because we also need to 
     * include the comment's CRC , see below
     */
  }/* compute crc64 checksum */

  /* ----- swap endian-ness if required ----- */
  if (swapEndian) 
    {
      endian_swap((CHAR*)(&rawheader.version), 			sizeof(rawheader.version) 		, 1);
      endian_swap((CHAR*)(&rawheader.gps_sec), 			sizeof(rawheader.gps_sec) 		, 1);
      endian_swap((CHAR*)(&rawheader.gps_nsec), 		sizeof(rawheader.gps_nsec) 		, 1);
      endian_swap((CHAR*)(&rawheader.tbase), 			sizeof(rawheader.tbase) 		, 1);
      endian_swap((CHAR*)(&rawheader.first_frequency_index), 	sizeof(rawheader.first_frequency_index) , 1);
      endian_swap((CHAR*)(&rawheader.nsamples), 		sizeof(rawheader.nsamples) 		, 1);

      /* v2-specific */
      endian_swap((CHAR*)(&rawheader.crc64),			sizeof(rawheader.crc64)			, 1);
      endian_swap((CHAR*)(&rawheader.comment_length),		sizeof(rawheader.comment_length)	, 1);
      /* ----- */

    } /* if endian_swap */

  /* double-check version-number */
  if ( rawheader.version != 2 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nWrong SFT-version %d in read_v2_header_from_fp()\n\n", rawheader.version );
      goto failed;
    }

  if ( rawheader.nsamples <= 0 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nNon-positive number of samples in SFT!\n\n");
      goto failed;
    }

  /* ----- v2-specific consistency-checks ----- */
  if ( rawheader.comment_length < 0 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nNegative comment-length in SFT!\n\n");
      goto failed;
    }

  if ( rawheader.comment_length % 8 != 0 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nComment-length must be multiple of 8 bytes!\n\n");
      goto failed;
    }

  if ( ! is_valid_detector ( rawheader.detector ) )
    {
      if ( lalDebugLevel ) LALPrintError ("\nIllegal detector-name in SFT: '%c%c'\n\n", 
					  rawheader.detector[0], rawheader.detector[1] );
      goto failed;
    }

  /* ----- Now read comment (if any) ----- */
  comm = NULL;
  if ( rawheader.comment_length )
    {
      CHAR *ptr;
      if ( (comm = LALCalloc(1, rawheader.comment_length )) == NULL )
	{
	  LALPrintError ("\nFATAL: out of memory ...\n\n");
	  goto failed;
	}
      if ( (size_t)rawheader.comment_length != fread ( comm, 1, rawheader.comment_length, fp ) )
	{
	  if ( lalDebugLevel ) LALPrintError ("\nCould not read %d-bytes comment\n\n");
	  goto failed;
	}

      /* check that comment is 0-terminated */
      if ( comm[ rawheader.comment_length - 1] != 0 )
	{
	  if ( lalDebugLevel ) LALPrintError ("\nComment is not properly 0-terminated!\n\n");
	  goto failed;
	}

      /* check that no NON-NULL bytes after first NULL in comment (->spec) */
      ptr = index ( comm, 0 );	/* guaranteed to find sth, after previous check */
      while ( ptr < (comm + rawheader.comment_length - 1) )
	if ( *ptr++ != 0 )
	  {
	    if ( lalDebugLevel ) LALPrintError ("\nNon-NULL bytes found after comment-end!\n\n");
	    goto failed;
	  }

      /* comment length including null terminator to string must be an
       * integer multiple of eight bytes. comment==NULL means 'no
       * comment'  
       */
      if (comment) 
	{
	  CHAR pad[] = {0, 0, 0, 0, 0, 0, 0};	/* for comment-padding */
	  UINT4 comment_len = strlen(comm) + 1;
	  UINT4 pad_len = (8 - (comment_len % 8)) % 8;
	  
	  crc = calc_crc64((const CHAR*)comm, comment_len, crc);
	  crc = calc_crc64((const CHAR*)pad, pad_len, crc);
	}
	      
    } /* if comment_length > 0 */

  /*  ok: */
  memset ( header, 0, sizeof( *header ) );

  header->name[0]		= rawheader.detector[0];
  header->name[1]		= rawheader.detector[1];
  header->name[2]		= 0;

  header->epoch.gpsSeconds 	= rawheader.gps_sec;
  header->epoch.gpsNanoSeconds 	= rawheader.gps_nsec;

  header->f0 			= rawheader.first_frequency_index / rawheader.tbase;
  header->deltaF 		= 1.0 / rawheader.tbase;

  (*nsamples) = rawheader.nsamples;
  (*ref_crc64) = rawheader.crc64;
  (*comment) = comm;
  (*header_crc64) = crc;


  return 0;

 failed:
  /* restore filepointer initial position  */
  if ( fseek ( fp, save_filepos, SEEK_SET ) == -1 )
    LALPrintError ("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno) );

  /* free comment  if we allocated one */
  if ( comm )
    LALFree (comm);

  return -1;

} /* read_v2_header_from_fp() */


/* ----- SFT v1 -specific header-reading function: 
 *
 * return general SFTtype header, place filepointer at the end of the header if it succeeds,
 * set fp to initial position if it fails.
 * RETURN: 0 = OK, -1 = ERROR 
 *
 * NOTE: fatal errors will produce a LALPrintError() error-message, but
 * non-fatal 'SFT-format'-errors will only output error-messages if lalDebugLevel > 0.
 *
 */
static int
read_v1_header_from_fp ( FILE *fp, SFTtype *header, UINT4 *nsamples, BOOLEAN swapEndian)
{
  _SFT_header_v1_t rawheader;
  long save_filepos;

  if ( !fp || !header || !nsamples)
    {
      LALPrintError ( "\nERROR read_v1_header_from_fp(): called with NULL input!\n\n");
      return -1;
    }
  
  /* store fileposition for restoring in case of failure */
  if ( ( save_filepos = ftell(fp) ) == -1 )
    {
      LALPrintError ("\nftell() failed: %s\n\n", strerror(errno) );
      return -1;
    }

  /* read the whole header */
  if (fread( &rawheader, sizeof(rawheader), 1, fp) != 1) 
    {
      if (lalDebugLevel) LALPrintError ("\nCould not read v1-header. %s\n\n", strerror(errno) );
      goto failed;
    }

  if (swapEndian) 
    {
      endian_swap((CHAR*)(&rawheader.version), 			sizeof(rawheader.version) 		, 1);
      endian_swap((CHAR*)(&rawheader.gps_sec), 			sizeof(rawheader.gps_sec) 		, 1);
      endian_swap((CHAR*)(&rawheader.gps_nsec), 		sizeof(rawheader.gps_nsec) 		, 1);
      endian_swap((CHAR*)(&rawheader.tbase), 			sizeof(rawheader.tbase) 		, 1);
      endian_swap((CHAR*)(&rawheader.first_frequency_index), 	sizeof(rawheader.first_frequency_index) , 1);
      endian_swap((CHAR*)(&rawheader.nsamples), 		sizeof(rawheader.nsamples) 		, 1);
    }

  /* double-check version-number */
  if ( rawheader.version != 1 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nWrong SFT-version %d in read_v1_header_from_fp()\n\n", rawheader.version );
      goto failed;
    }

  if ( rawheader.nsamples <= 0 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nNon-positive number of samples in SFT!\n\n");
      goto failed;
    }


  /* ok: */
  memset ( header, 0, sizeof( *header ) );

  /* NOTE: v1-SFTs don't contain a detector-name, in which case we set it to '??' */
  strcpy ( header->name, "??" );

  header->epoch.gpsSeconds 	= rawheader.gps_sec;
  header->epoch.gpsNanoSeconds 	= rawheader.gps_nsec;
  header->deltaF 		= 1.0 / rawheader.tbase;
  header->f0 			= rawheader.first_frequency_index / rawheader.tbase;
  
  (*nsamples) = rawheader.nsamples;

  return 0;

 failed:
  /* restore filepointer initial position  */
  if ( fseek ( fp, save_filepos, SEEK_SET ) == -1 )
    LALPrintError ("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno) );

  return -1;

} /* read_v1_header_from_fp() */

/* check that channel-prefix defines a 'known' detector.  The list of
 * known detectors implemented here for now follows the list in
 * Appendix D of LIGO-T970130-F-E:
 *
 * returns TRUE if valid, FALSE otherwise */
static BOOLEAN
is_valid_detector (const char *channel) 
{
  int i;
  const char *knownDetectors[] = 
    {
      "A1",       /* ALLEGRO */
      "B1",       /* NIOBE */
      "E1",       /* EXPLORER */
      "G1",       /* GEO_600 */
      "H1",       /* LHO_4k */
      "H2",       /* LHO_2k */
      "K1",       /* ACIGA */
      "L1",       /* LLO_4k */
      "N1",       /* Nautilus */
      "O1",       /* AURIGA */
      "P1",       /* CIT_40 */
      "T1",       /* TAMA_300 */
      "V1",       /* Virgo_CITF */
      "V2",       /* Virgo (3km) */
      "Z1",	  /* LISA effective IFO 1 */
      "Z2",	  /* LISA effective IFO 2 */
      "Z3",	  /* LISA effective IFO 3 */
      NULL
    };
  
  if ( !channel )
    return FALSE;

  if ( strlen(channel) < 2 )
    return FALSE;

  for ( i = 0; knownDetectors[i]; i ++ ) 
    {
      if ( ( knownDetectors[i][0] == channel[0] ) && ( knownDetectors[i][1] == channel[1] )  )
	return TRUE;
    }
  
  return FALSE;

} /* is_valid_detector() */


/* a little endian-swapper needed for SFT reading/writing */
static void 
endian_swap(CHAR * pdata, size_t dsize, size_t nelements)
{
  UINT4 i, j, indx;
  CHAR tempbyte;

  if (dsize <= 1) return;

  for (i=0; i<nelements; i++)
    {
      indx = dsize;
      for (j=0; j<dsize/2; j++)
	{
	  tempbyte = pdata[j];
	  indx = indx - 1;
	  pdata[j] = pdata[indx];
	  pdata[indx] = tempbyte;
	}
      
      pdata = pdata + dsize;
    }
  
  return;

} /* endian swap */

/*----------------------------------------------------------------------
 * glob() has been reported to fail under condor, so we use our own
 * function to get a filelist from a directory, using a glob-like pattern.
 *
 * looks pretty ugly with all the #ifdefs for the Microsoft C compiler
 *
 * NOTE: the list of filenames is returned SORTED ALPHABETICALLY !
 *
 *----------------------------------------------------------------------*/
static LALStringVector *
find_files (const CHAR *globdir)
{
#ifndef _MSC_VER
  DIR *dir;
  struct dirent *entry;
#else
  intptr_t dir;
  struct _finddata_t entry;
#endif
  CHAR *dname, *ptr1;
  CHAR *fpattern;
  size_t dirlen;
  CHAR **filelist = NULL; 
  UINT4 numFiles = 0;
  LALStringVector *ret = NULL;
  UINT4 j;
  UINT4 namelen;
  CHAR *thisFname;

  /* First we separate the globdir into directory-path and file-pattern */

#ifndef _MSC_VER
#define DIR_SEPARATOR '/'
#else
#define DIR_SEPARATOR '\\'
#endif

  /* any path specified or not ? */
  ptr1 = strrchr (globdir, DIR_SEPARATOR);
  if (ptr1)
    { /* yes, copy directory-path */
      dirlen = (size_t)(ptr1 - globdir) + 1;
      if ( (dname = LALCalloc (1, dirlen)) == NULL) 
	return (NULL);
      strncpy (dname, globdir, dirlen);
      dname[dirlen-1] = '\0';

      ptr1 ++; /* skip dir-separator */
      /* copy the rest as a glob-pattern for matching */
      if ( (fpattern = LALCalloc (1, strlen(ptr1) + 1)) == NULL )
	{
	  LALFree (dname);
	  return (NULL);
	}
      strcpy (fpattern, ptr1);   

    } /* if ptr1 */
  else /* no pathname given, assume "." */
    {
      if ( (dname = LALCalloc(1, 2)) == NULL) 
	return (NULL);
      strcpy (dname, ".");

      if ( (fpattern = LALCalloc(1, strlen(globdir)+1)) == NULL)
	{
	  LALFree (dname);
	  return (NULL);
	}
      strcpy (fpattern, globdir);	/* just file-pattern given */
    } /* if !ptr */
  

#ifndef _MSC_VER
  /* now go through the file-list in this directory */
  if ( (dir = opendir(dname)) == NULL) {
    LALPrintError ("Can't open data-directory `%s`\n", dname);
    LALFree (dname);
    return (NULL);
  }
#else
  if ((ptr1 = (CHAR*)LALMalloc(strlen(dname)+3)) == NULL)
    return(NULL);
  sprintf(ptr1,"%s\\*",dname);  
  dir = _findfirst(ptr1,&entry);
  LALFree(ptr1);
  if (dir == -1) {
    LALPrintError ("Can't find file for pattern `%s`\n", ptr1);
    LALFree (dname);
    return (NULL);
  }
#endif

#ifndef _MSC_VER
  while ( (entry = readdir (dir)) != NULL )
    {
      thisFname = entry->d_name;
#else
  do
    {
      thisFname = entry.name;
#endif
      /* now check if glob-pattern fpattern matches the current filename */
      if ( amatch(thisFname, fpattern) 
	   /* and check if we didnt' match some obvious garbage like "." or ".." : */
	   && strcmp( thisFname, ".") && strcmp( thisFname, "..") )
	{

	  numFiles ++;
	  if ( (filelist = LALRealloc (filelist, numFiles * sizeof(CHAR*))) == NULL) {
	    LALFree (dname);
	    LALFree (fpattern);
	    return (NULL);
	  }

	  namelen = strlen(thisFname) + strlen(dname) + 2 ;

	  if ( (filelist[ numFiles - 1 ] = LALCalloc (1, namelen)) == NULL) {
	    for (j=0; j < numFiles; j++)
	      LALFree (filelist[j]);
	    LALFree (filelist);
	    LALFree (dname);
	    LALFree (fpattern);
	    return (NULL);
	  }

	  sprintf(filelist[numFiles-1], "%s%c%s", dname, DIR_SEPARATOR, thisFname);

	} /* if filename matched pattern */

    } /* while more directory entries */
#ifdef _MSC_VER
  while ( _findnext (dir,&entry) == 0 );
#endif

#ifndef _MSC_VER
  closedir (dir);
#else
  _findclose(dir);
#endif

  LALFree (dname);
  LALFree (fpattern);

  /* ok, did we find anything? */
  if (numFiles == 0)
    return (NULL);

  if ( (ret = LALCalloc (1, sizeof (LALStringVector) )) == NULL) 
    {
      for (j=0; j<numFiles; j++)
	LALFree (filelist[j]);
      LALFree (filelist);
      return (NULL);
    }

  ret->length = numFiles;
  ret->data = filelist;

  /* sort this alphabetically (in-place) */
  SortStringVector (ret);

  return (ret);
} /* find_files() */


static void
DestroyStringVector (LALStringVector *strings)
{
  UINT4 i;

  for (i=0; i < strings->length; i++)
    LALFree (strings->data[i]);
  
  LALFree (strings->data);
  LALFree (strings);

} /* DestroyStringVector () */

/* comparison function for strings */
static int mycomp (const void *p1, const void *p2)
{
  const char *s1 = p1;
  const char *s2 = p2;
  return (strcmp ( s1, s2 ) );
}

/* sort string-vector alphabetically */
static void
SortStringVector (LALStringVector *strings)
{
  qsort ( (void*)(strings->data), (size_t)(strings->length), sizeof(CHAR*), mycomp );
} /* SortStringVector() */

/* portable file-len function */
static long get_file_len ( FILE *fp )
{
  long save_fp;
  long len;

  if ( (save_fp = ftell(fp)) == -1 )
    return 0;
  
  if ( fseek ( fp, 0, SEEK_END ) == -1 )
    return 0;

  len = ftell(fp);

  if ( fseek ( fp, save_fp, SEEK_SET ) == -1 )
    return 0;

  return len;

} /* get_file_len() */





/* ----- the following function crc64() was taken from SFTReferenceLibrary.c
 * and adapted to LAL .
 *
 *  The quantity below is: D800000000000000 (base-16) =
 *  1101100000000000000000000000000000000000000000000000000000000000
 *  (base-2).  The primitive polynomial is x^64 + x^4 + x^3 + x + 1.
 */
#define POLY64 0xd800000000000000ULL
#define TABLELEN 256

/* The crc64 checksum of M bytes of data at address data is returned
 * by crc64(data, M, ~(0ULL)). Call the function multiple times to
 * compute the checksum of data made in contiguous chunks, setting
 * final argument to the previously accumulated checksum value. */
static UINT8 
calc_crc64(const CHAR *data, UINT4 length, UINT8 crc) 
{
  UINT8 CRCTable[TABLELEN];
  UINT4 i;

  /* is there is no data, simply return previous checksum value */
  if (!length || !data )
    return crc;
  
  /* initialize the CRC table for fast computation.  We could keep
     this table in memory to make the computation faster, but that is
     not re-entrant for multi-threaded code.
  */
  for (i = 0; i < TABLELEN; i++) {
    UINT4 j;
    UINT8 part = i;
    for (j = 0; j < 8; j++) {
      if (part & 1)
        part = (part >> 1) ^ POLY64;
      else
        part >>= 1;
    }
    CRCTable[i] = part;
  }
  
  /* compute the CRC-64 code */
  for (i=0; i<length; i++) {
    UINT8 temp1 = crc >> 8;
    UINT8 temp2 = CRCTable[(crc ^ (UINT8) data[i]) & 0xff];
    crc = temp1 ^ temp2;
  }
  
  return crc;

} /* calc_crc64() */


/* compare two SFT-descriptors by their GPS-epoch */
static int 
compareSFTdesc(const void *ptr1, const void *ptr2)
{
  const SFTDescriptor *desc1 = ptr1;
  const SFTDescriptor *desc2 = ptr2;

  if ( GPS2REAL8( desc1->header.epoch ) < GPS2REAL8 ( desc2->header.epoch ) )
    return -1;
  else if ( GPS2REAL8( desc1->header.epoch ) > GPS2REAL8 ( desc2->header.epoch ) )
    return 1;
  else
    return 0;
  
} /* compareSFTdesc() */

/** Read valid SFT version-number at position fp, and determine if we need to 
 * endian-swap the data.
 * Restores filepointer to original position before returning.
 * 
 * RETURN: 0 = OK, -1 = ERROR
 */
int
read_SFTversion_from_fp ( UINT4 *version, BOOLEAN *need_swap, FILE *fp )
{
  long save_filepos;
  REAL8 ver;

  /* store fileposition for restoring in case of failure */
  if ( ( save_filepos = ftell(fp) ) == -1 )
    {
      LALPrintError ("\nftell() failed: %s\n\n", strerror(errno) );
      return -1;
    }

  /* read version-number */
  if  ( 1 != fread ( &ver, sizeof( ver ), 1, fp) )
    {
      if (lalDebugLevel) LALPrintError ("\nCould not read version-number from file\n\n");
      goto failed;
    }


  /* figure out endian-ness and check version-range */
  for ( *version = MAX_SFT_VERSION; *version >= MIN_SFT_VERSION; --(*version) )
    {
      REAL8 vertest = *version;
      if ( ! memcmp( &ver, &vertest, sizeof( ver ) ) ) {
	*need_swap = FALSE;
	break;
      }
      endian_swap( (char*)(&vertest), sizeof( vertest ), 1 );
      if ( ! memcmp( &ver, &vertest, sizeof( ver ) ) ) {
	*need_swap = TRUE;
	break;
      }
    }
  if ( *version < MIN_SFT_VERSION ) {
    if ( lalDebugLevel )
      LALPrintError( "\nERROR: illegal SFT-version not within [%.0f, %.0f]\n\n", MIN_SFT_VERSION, MAX_SFT_VERSION );
    goto failed;
  }
  
  /* restore initial filepointer position */
  if ( fseek ( fp, save_filepos, SEEK_SET ) == -1 )
    {
      LALPrintError ("\nfseek() failed to return to intial fileposition: %s\n\n", strerror(errno) );
      goto failed;
    }

  return 0;

 failed:
  fseek ( fp, save_filepos, SEEK_SET );
  return -1;

} /* read_SFTversion_from_fp() */


/*======================================================================*/
/*
 * robust glob pattern matcher
 * ozan s. yigit/dec 1994
 * public domain
 *
 * glob patterns:
 *	*	matches zero or more characters
 *	?	matches any single character
 *	[set]	matches any character in the set
 *	[^set]	matches any character NOT in the set
 *		where a set is a group of characters or ranges. a range
 *		is written as two characters seperated with a hyphen: a-z denotes
 *		all characters between a to z inclusive.
 *	[-set]	set matches a literal hypen and any character in the set
 *	[]set]	matches a literal close bracket and any character in the set
 *
 *	char	matches itself except where char is '*' or '?' or '['
 *	\char	matches char, including any pattern character
 *
 * examples:
 *	a*c		ac abc abbc ...
 *	a?c		acc abc aXc ...
 *	a[a-z]c		aac abc acc ...
 *	a[-a-z]c	a-c aac abc ...
 *
 */

#ifndef NEGATE
#define NEGATE	'^'			/* std cset negation char */
#endif

static int amatch(char *str, char *p);

static int
amatch(char *str, char *p)
{
	int negate;
	int match;
	int c;

	while (*p) {
		if (!*str && *p != '*')
			return FALSE;

		switch (c = *p++) {

		case '*':
			while (*p == '*')
				p++;

			if (!*p)
				return TRUE;

			if (*p != '?' && *p != '[' && *p != '\\')
				while (*str && *p != *str)
					str++;

			while (*str) {
				if (amatch(str, p))
					return TRUE;
				str++;
			}
			return FALSE;

		case '?':
			if (*str)
				break;
			return FALSE;
/*
 * set specification is inclusive, that is [a-z] is a, z and
 * everything in between. this means [z-a] may be interpreted
 * as a set that contains z, a and nothing in between.
 */
		case '[':
			if (*p != NEGATE)
				negate = FALSE;
			else {
				negate = TRUE;
				p++;
			}

			match = FALSE;

			while (!match && (c = *p++)) {
				if (!*p)
					return FALSE;
				if (*p == '-') {	/* c-c */
					if (!*++p)
						return FALSE;
					if (*p != ']') {
						if (*str == c || *str == *p ||
						    (*str > c && *str < *p))
							match = TRUE;
					}
					else {		/* c-] */
						if (*str >= c)
							match = TRUE;
						break;
					}
				}
				else {			/* cc or c] */
					if (c == *str)
						match = TRUE;
					if (*p != ']') {
						if (*p == *str)
							match = TRUE;
					}
					else
						break;
				}
			}

			if (negate == match)
				return FALSE;
/*
 * if there is a match, skip past the cset and continue on
 */
			while (*p && *p != ']')
				p++;
			if (!*p++)	/* oops! */
				return FALSE;
			break;

		case '\\':
			if (*p)
				c = *p++;
		default:
			if (c != *str)
				return FALSE;
			break;

		}
		str++;
	}

	return !*str;
}


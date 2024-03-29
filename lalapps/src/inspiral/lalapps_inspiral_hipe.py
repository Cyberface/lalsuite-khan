"""
inspiral_hipe.in - standalone inspiral pipeline driver script

This script generates the condor DAG necessary to analyze LIGO and GEO
data through the inspiral pipeline.  The DAG does datafind, tmpltbank,
inspiral, inca and sire/coire steps of the pipeline.  It analyzes the
single, double, triple and quadro ifo times accordingly.  It can also be
run with injections.
"""

## \file
#
# <dl>
# <dt>Name</dt><dd>
# \c lalapps_inspiral_hipe --- python script to generate Condor DAGs to
# run the inspiral hierarchical pipeline.</dd>
#
# <dt>Synopsis</dt><dd>
# \code
#   -h, --help            show this help message and exit
#   -v, --version         print version information and exit
#   -u  USERTAG, --user-tag= USERTAG
#                         tag the jobs with USERTAG (overrides value in ini
#                         file)
#   -g, --g1-data         analyze g1 data
#   -a, --h1-data         analyze h1 data
#   -b, --h2-data         analyze h2 data
#   -l, --l1-data         analyze l1 data
#   -S, --one-ifo         analyze single ifo data (not usable for GEO)
#   -D, --two-ifo         analyze two interferometer data
#   -T, --three-ifo       analyze three interferometer data
#   -Q, --four-ifo        analyze four intereferometer data
#   -A, --analyze-all     analyze all ifos and all data (over-rides above)
#   -d, --datafind        run LSCdataFind to create frame cache files
#   -t, --template-bank   run lalapps_tmpltbank to generate template banks
#   -i, --inspiral        run lalapps_inspiral to generate triggers
#   -c, --coincidence     run lalapps_thinca to test for coincidence
#   -B, --trigbank        run lalapps_trigbank for banks of coinc triggers
#   -V, --inspiral-veto   run lalapps_inspiral with vetos
#   -C, --second-coinc    run lalapps_thinca on the inspiral veto triggers
#   -j, --coherent-bank   run lalapps_coherentbank to make coherent bank
#   -k, --coherent-inspiral
#                         run lalapps_coherent_inspiral for coherent analysis
#   -s, --sire            do sires to sweep up triggers
#   -R, --read-cache      read cache file from ini-file (if LSCDataFind is
#                         broken)
#   -P  PRIO, --priority= PRIO
#                         run jobs with condor priority PRIO
#   -f  FILE, --config-file= FILE
#                         use configuration file FILE
#   -p  PATH, --log-path= PATH
#                         directory to write condor log file
#   -o, --output-segs     output the segment lists of analyzed data
#   -x, --dax             create a dax instead of a dag
# \endcode</dd>
#
# <dt>Description</dt><dd> \c lalapps_inspiral_hipe generates a Condor DAG to run
# the hierarchical inspiral analysis pipeline.  It currently works for
# the four LSC interferometers: G1, H1, H2, L1.
#
# The code reads in segment lists for the four instruments.  If one of the
# segment files is not specified or is empty, it is assumed that there is no
# data from that instrument.  From the segment files, the pipeline calculates
# four lists of single ifo segments, for G1, H1, H2 and L1; six lists of double
# ifo segments, for G1-H1, G1-H2, G1-L1, H1-H2, H1-L1 and H2-L1; four lists of
# three ifo data, for G1-H1-H2, G1-H1-L1, G1-H2-L1 and H1-H2-L1, and one list of
# four ifo segments for G1-H1-H2-L1.  The options <tt>--g1-data</tt>,
# <tt>--h1-data</tt>, <tt>--h2-data</tt> and <tt>--l1-data</tt> allow you to choose
# which of the interferometers' data to analyze.  Similarly, the
# <tt>--one-ifo</tt>, <tt>--two-ifo</tt>, <tt>--three-ifo</tt> and <tt>--four-ifo</tt>
# flags determine whether to analyze times during which one, two, three or four
# instruments respectively were operational.  Thus, by specifying
# <tt>--h1-data</tt>, <tt>--l1-data</tt> and <tt>--two-ifo</tt>, the pipeline will
# analyze only the H1-L1 double coincident times.  If the <tt>--analyze-all</tt>
# flag is set, the pipeline will analyze all data from all instruments.  If the
# <tt>--output-segments</tt> option is chosen, the pipeline will output segment
# lists for the non-empty data types.  The file names are
# "h1_l1_segs_analyzed.txt" etc, or if the analyis is restricted to
# playground, they are "h1_l1_play_segs_analyzed.txt".
#
# The pipeline uses a coincidence stage early on in order to cut down the number
# of triggers for which we have to perform the computationally costly chi
# squared and r-squared veto computations.  Thus, the pipeline performs a first
# inspiral stage (without vetoes) which is followed immediately by a coincidence
# stage.  The triggers which survive in coincidence are then passed back to the
# inspiral code where the chi-squared and r-squared signal based vetoes are
# computed.  The remaining triggers are then tested again for coincidence.  Any
# triggers surviving this second coincidence stage are passed into the coherent
# analysis.  Some of the steps described above require more than one code to be
# run.  For example, before running the inspiral code for the first time, we
# must first locate the data using datafind, then generate template-banks for
# the analysis.
#
# At present, the pipeline can perform the following steps of a hierarchical
# inspiral search: <tt>--datafind</tt>, <tt>--template-bank</tt>, <tt>--inspiral</tt>,
# <tt>--coincidence</tt>, <tt>--trigbank</tt>, <tt>--inspiral-veto</tt>,
# <tt>--second-coinc</tt>, <tt>--coherent-bank</tt> and <tt>--coherent-inspiral</tt>.
# Any or all of these options may be specified.  However, each step of the
# pipeline relies on results files produced by the previous step (and in the
# case of the \c inspiral, <tt>inspiral-veto</tt> and <tt>coherent-inspiral</tt>
# steps, the output from \c datafind is also required).
#
# The configuration file specifies the parameters needed to run the analysis jobs
# contained in the pipeline.  It is specified with the <tt>--config-file</tt>
# option.  A typical .ini file is the following:
#
# TODOinput{hipeinifile}
#
# The .ini file contains several sections.  The <tt>[condor]</tt> section contains
# the names of the executables which will run the various stages of the
# pipeline.  The <tt>[pipeline]</tt> section gives the CVS details of the
# pipeline, the usertag (which can be overwritten on the command line with the
# <tt>--user-tag</tt> option) and the <tt>playground-data-mask</tt> which must be
# set to one of \c playground_only, \c exclude_playground or
# \c all_data.  The \c input section contains the names of the segment
# files for the four interferometers, the channel names, frame types, injection
# file name and number of time slides.  If any of the segment files are left
# blank, it is assumed that there are no segments for that instrument.
# Similarly, a blank injection-file signifies that no injections are to be
# performed, and a blank num-slides signifies no time slides.
#
# The remaining sections set options for the various jobs to be run in the
# pipeline.  The options in the <tt>[datafind]</tt>, <tt>[tmpltbank]</tt>,
# <tt>[inspiral]</tt>, <tt>[inca]</tt>, <tt>[thinca]</tt>, <tt>[trigtotmplt]</tt>,
# <tt>[cohbank]</tt> and <tt>[chia]</tt> sections are added to every instance of the
# relevant executable.  Note that these options are set the same for all
# interferometers.  The options in the <tt>[data]</tt> section are added to all
# <tt>[inspiral]</tt> and <tt>[tmpltbank]</tt> jobs, while the <tt>[ligo-data]</tt>
# and <tt>[geo-data]</tt> commands are added to the LIGO/GEO jobs respectively.
# The <tt>[calibration]</tt> information is used to determine the calibration data
# for these jobs.  In the <tt>[ifo-thresholds]</tt> sections, the ifo specific
# snr, chi squared and r-squared veto thresholds for the various intererometers
# are set; while the <tt>[no-veto-inspiral]</tt> and <tt>veto-inspiral]</tt>
# sections contain arguments relevant for the first and second inspiral steps
# respectively.
#
# The science segments are read in from the segment files for each instrument.
# These science segments are split up into analysis chunks.  The analysis chunk
# size is determined from the number of data segments and their length and
# overlap specified in config file. Currently, we are using 256 second analysis
# segments.  We use 15 segments in a chunk, and overlap them by 128 seconds to
# give a chunk length of 2048 seconds.  The chunks are constructed for each of
# the interferometers independently.  Any science segment shorter than the
# length of a chunk is not analyzed.  Additionally, we cannot produce triggers
# for the first and last overlap/2 seconds of a science segment, due to the
# finite length of the inspiral templates.  Using this information, we construct
# segment lists of analyzable data for each of the fifteen types of data we may
# have (four single ifo, six two ifo, four three ifo and one four ifo).  If the
# playground only option is specified, the segments are restricted to playground
# times.  We decide which chunks should be analyzed by testing for overlap
# between the chunk and the data we need to analyze.  Note that if the pipeline
# is restricted to playground data, then only the playground times are analyzed
# for triggers in the inspiral code.  This is done by setting the
# <tt>trig-start-time</tt> and <tt>trig-end-time</tt> for the inspiral jobs
# appropriately.
#
# Once the DAG file has been created it should be submitted to the Condor pool
# with the \c condor_submit_dag command.</dd>
#
# <dt>Options</dt><dd>
# <ul>
#
# <li><tt>--help</tt>: Display a brief usage summary.</li>
#
# <li><tt>--version</tt>: Display the version information and exit.</li>
#
# <li><tt>--user-tag</tt> \c usertag: Set the user-tag to \c usertag.
# This overrides the user-tag which may have been set in the ini file.  The usertag
# will be added to the names of the output files of the pipeline.</li>
#
# <li><tt>--g1-data</tt>: Analyze the G1 data, the times of which are
# determined from the <tt>g1-segments</tt> file specified in the ini file.  If not
# set, then no data when G1 was operational will be analyzed.</li>
#
# <li><tt>--h1-data</tt>: Analyze the H1 data, the times of which are
# determined from the <tt>h1-segments</tt> file specified in the ini file.  If not
# set, then no data when H1 was operational will be analyzed.</li>
#
# <li><tt>--h2-data</tt>: Analyze the H2 data, the times of which are
# determined from the <tt>h2-segments</tt> file specified in the ini file.  If not
# set, then no data when H2 was operational will be analyzed.</li>
#
# <li><tt>--l1-data</tt>: Analyze the L1 data, the times of which are
# determined from the <tt>l1-segments</tt> file specified in the ini file.  If not
# set, then no data when H1 was operational will be analyzed.</li>
#
# <li><tt>--one-ifo</tt>: Analyze any times when one and only one instrument
# was operational.  Note that this option works together with the IFO options
# given above.  For example if <tt>--one-ifo</tt> and <tt>--h2-data</tt> were
# specified, then only the single IFO H2 times would be analyzed.</li>
#
# <li><tt>--two-ifo</tt>: Analyze any times when two instruments were
# operational.  Note that this option works together with the IFO options given
# above.  For example if <tt>--two-ifo</tt>, <tt>h1-data</tt> and <tt>--h2-data</tt>
# were specified, then the times when only H1 and H2 were operational would be
# analyzed.  However, if only <tt>--two-ifo</tt> and <tt>h1-data</tt> were
# specified, no data would be analyzed.</li>
#
# <li><tt>--three-ifo</tt>: Analyze any times when three instruments were
# operational.  Note that this option works together with the IFO options given
# above.  For example if <tt>--three-ifo</tt>, <tt>h1-data</tt>, <tt>--h2-data</tt>
# and <tt>--l1-data</tt> were specified, then the times when H1, H2 and L1
# were operational, but G1 was not, would be analyzed. </li>
#
# <li><tt>--four-ifo</tt>: Analyze any times when all four instruments were
# operational.  Note that this option works together with the IFO options given
# above.  For example if <tt>--four-ifo</tt>, <tt>g1-data</tt>, <tt>h1-data</tt>,
# <tt>--h2-data</tt> and <tt>--l1-data</tt> were specified, then the times when all
# of G1,H1, H2 and L1 were operational would be analyzed.</li>
#
# <li><tt>--analyze-all</tt>: Analyze all ifos and all data.  This is
# equivalent to setting all six of the options above.  Then, all the data is
# analyzed.</li>
#
# <li><tt>--datafind</tt>: Run the datafind step of the pipeline.</li>
#
# <li><tt>--template-bank</tt>: Run the template-bank step of the pipeline.
# Note that the template-bank jobs require the cache files created by datafind,
# so <tt>--datafind</tt> must either be run in the pipeline or have been run
# previously.</li>
#
# <li><tt>--inspiral</tt>: Run the inspiral step of the pipeline.  These jobs
# will take the arguments from the <tt>[no-veto-inspiral]</tt> section in the ini
# file.  Note that the inspiral jobs require the cache files created by datafind
# and template banks, so both <tt>--datafind</tt> and <tt>template-bank</tt> must
# either be run in the pipeline or have been run previously.</li>
#
# <li><tt>--coincidence</tt>: Run the coincidence step of the pipeline.  Times
# when two or more detectors were operational are analyzed by the thinca code.
# This determines all coincidences between two or more operational detectors.
# If <tt>num-slides</tt> is specified then time slides are also performed, and
# output in a separate file named THINCA_SLIDE.  Finally, for the times when
# only one instrument was on, inca is run (in single ifo mode) to sweep up the
# single ifo triggers.  These triggers are not used by later stages of the
# pipeline.  Note that the thinca and inca jobs require the inspiral triggers
# created by inspiral, so that <tt>--inspiral</tt> must either be run in the
# pipeline or have been run previously.  For each single IFO segment, the
# coincidence step simply creates a file containing all the triggers in the time
# interval.  For two/three/four IFO segments, the coincidence step performs
# coincidence and outputs the double, triple and quadruple coincidence triggers
# in one file, and time slide coincident triggers in a second file.</li>
#
# <li><tt>--trigbank</tt>: Run the triggered bank step of the pipeline.  This
# step takes in the coincident, and thinca slide, triggers.  It outputs those
# triggers which should be filtered by a specific instrument in the follow-up
# inspiral stage.  This is done by keeping those triggers from the relevant ifo,
# within the times of the chunk to be analyzed by the inspiral code.</li>
#
# <li><tt>--inspiral-veto</tt>: Run the second inspiral step of the pipeline.
# These jobs will take the arguments from the <tt>[veto-inspiral]</tt> section in
# the ini file, and are intended to do the computationally costly signal based
# vetoes of coincident triggers.  These jobs are differentiated from the first
# inspiral jobs by the inclusion of an ifo-tag, so an example job may be named
# <tt>H1-INSPIRAL_H1L1-GPSTIME-DURATION.xml</tt>.  Note that the inspiral jobs
# require the cache files created by datafind and trigbanks, so both
# <tt>--datafind</tt> and \c trigbank must either be run in the pipeline or
# have been run previously.</li>
#
# <li><tt>--second-coinc</tt>: Re-run the coincidence step of the pipeline.
# This runs the thinca code on all times when two or more instruments were
# operational.  As with the first coinc stage, this will perform both the zero
# lag and time slide analysis if requested.  The output from these jobs also
# have an ifo-tag added, so an example output file might be
# <tt>H1H2L1-THINCA_H1H2L1-GPSTIME-DURATION.xml</tt> for the zero lag and
# <tt>H1H2L1-THINCA_SLIDE_H1H2L1-GPSTIME-DURATION.xml</tt> for the time slides.</li>
#
# <li><tt>--coherent-bank</tt>: Run the coherent bank step of the pipeline.
# This generates template banks ready for the coherent analysis.  These banks
# are generated from the triggers output in the second coincidence stage of the
# pipeline.</li>
#
# <li><tt>--coherent-inspiral</tt>: Run the coherent inspiral step of the
# pipeline.  This runs the inspiral code, using the coherent-bank.  The inspiral
# code outputs frames containing time series of snr data around the time of a
# coincident event.  These frames are then read in by the coherent inspiral code
# which calculates the coherent signal to noise ratio of the event.</li>
#
# <li><tt>--sire</tt>: Run sire on the various stages of the pipeline.  This
# will collect together the triggers from jobs performed in various stages of
# the pipeline</li>
#
# <li><tt>--read-cache</tt>: Specify a static cache file to use during the
# analysis.  This should only be used if LSCdataFind is not available where the
# pipeline is being run.</li>
#
# <li><tt>--priority</tt> \c PRIO: Set the condor priority \c PRIO
# of the condor jobs to be run in the pipeline.</li>
#
# <li><tt>--config-file</tt> \c config_file: Set the name of the
# configuration file to be \c config_file.  This is the which is used to
# determine the parameters for the pipeline.  This is a required argument.</li>
#
# <li><tt>--log-path</tt>: The directory in which to write the condor log
# file.  This should generally be a local directory of the condor submit
# machine.  This is a required argument.</li>
#
# <li><tt>--output-segs</tt>: Output the segment lists of analyzed data.  Up
# to seven files will be output, one for each of the types of interferometer
# data (H1, H2, L1, H1_H2, H1_L1, H2_L1, H1_H2_L1).  Any segment lists
# which are non-empty will we written.</li>
#
# <li><tt>--dax</tt>: Output a dax rather than a dag. </li>
#
# </ul> </dd>
#
# <dt>Author</dt><dd>
# Steve Fairhurst, Darren Woods</dd>
# </dl>

from __future__ import print_function

__author__ = 'Stephen Fairhurst <sfairhur@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

##############################################################################
# import standard modules
import sys, os, copy, math, shutil
import socket, time
import re, string
from optparse import *
import tempfile
import urlparse

##############################################################################
# import the modules we need to build the pipeline
from glue import pipeline
from six.moves import configparser
from glue.pipeline import DeepCopyableConfigParser as dcConfigParser
from glue import lal
from lalapps import inspiral
from lalapps.inspiralutils import get_data_options

##############################################################################
# concrete dag setup functions
def write_local_tc(cp):
  """
  Function to add the executables that this DAX needs to the transformation 
  catalog so pegasus can stage them to the remote sites.
  """
  # open the transformation catalog file
  try:
    tc = open('tc.data','w')
  except:
    print("Cannot open transformation catalog for writing", file=sys.stderr)
    sys.exit(1)

  # write a line to the transformation catalog for each executable
  host = socket.getfqdn()
  programs = cp.items('condor')
  for prog in programs:
    if prog[0] != '/':
      # relative path, so append the pwd to the path
      prog = os.path.join(os.getcwd(), prog)
    tc.write("local ligo::%s:1.0 gsiftp://" + os.path.join(host,prog))
    tc.write(" STATIC_BINARY INTEL32::LINUX\n")
  tc.close()

def create_local_pool(cp):
  """
  Function to set up the local pool entries in to pool catalog.
  """
  pass

def create_static_pfn_cache(cp):
  """
  Function to create a static PFN cache for any files needed by not in RLS.
  """
  pass


##############################################################################
# some functions to make life easier later
class AnalyzedIFOData:
  """
  Contains the information for the data that needs to be filtered.
  """
  def __init__(self,chunk,node, inj=0):
    self.__analysis_chunk = chunk
    self.__dag_node = node
    self.__inj = inj

  def set_chunk(self,chunk):
    self.__analysis_chunk = chunk

  def get_chunk(self):
    return self.__analysis_chunk

  def set_dag_node(self,node):
    self.__dag_node = node

  def get_dag_node(self):
    return self.__dag_node

  def set_inj(self, inj):
    self.__inj = inj

  def get_inj(self):
    return self.__inj

##############################################################################
# function to generate the coincident segments from four input segment lists
def generate_segments(ifo1_data, ifo2_data, ifo3_data, ifo4_data):
  """
    compute the segments arising as the overlap of the four sets of single
    ifo segment lists.
    ifo1_data = data segments for ifo1
    ifo2_data = data segments for ifo2
    ifo3_data = data segments for ifo3
    ifo4_data = data segments for ifo4
  """ 

  segment_list = pipeline.ScienceData()
  segment_list = copy.deepcopy(ifo1_data)
  segment_list.intersect_4(ifo2_data, ifo3_data, ifo4_data)
   
  return segment_list

    
##############################################################################
# function to set up datafind, template bank and inspiral jobs for an ifo  
def analyze_ifo(ifo_name,ifo_data,ifo_to_do,tmplt_job,insp_job,df_job,\
  prev_df,dag, exttrigInjections, usertag=None, inspinjNode = None,
  insp_ckpt_job = None):
  """
  Analyze the data from a single IFO.  Since the way we treat all this data is
  the same, this function is the same for all interferometers. Returns the last
  LSCdataFind job that was executed and the chunks analyzed.
  
  ifo_name = the name of the IFO
  ifo_data = the master science segs 
  ifo_to_do = the science segments we need to analyze
  tmplt_job = if not GeoBank: template bank job we should use 
  insp_job = the condor job that we should use to analyze data
  df_job = the condor job to find the data
  prev_df = the previous LSCdataFind job that was executed
  dag = the DAG to attach the nodes to
  exttrigInjections = a two-element list specifying the range 
                      of injections for the external trigger case.
  usertag = the usertag to add to the job names
  inspinjNode = the inspinj node to be added as a parent to inspirals
  """

  # add the non veto inspiral options
  if cp.has_section('no-veto-inspiral'): 
    insp_job.add_ini_opts(cp,'no-veto-inspiral')
  
  # add the ifo specific options
  if cp.has_section(ifo_name.lower() + '-inspiral'): 
    insp_job.add_ini_opts(cp,ifo_name.lower() + '-inspiral')

  if cp.has_section(ifo_name.lower() + '-tmpltbank'):
    tmplt_job.add_ini_opts(cp,ifo_name.lower() + '-tmpltbank')

  GeoBank = None

  # to analyze GEO data we may use a fixed bank specified in ini file
  if ifo_name == 'G1':
    try: GeoBank = cp.get('input','geo-bank')
    except: GeoBank = None
    if GeoBank: print("For G1 we use bank ", GeoBank)

  data_opts, type, channel = get_data_options(cp,ifo_name)

  if cp.has_section('tmpltbank-1'):
    tmplt_job.add_ini_opts(cp, 'tmpltbank-1')
  if cp.has_section(data_opts):
    tmplt_job.add_ini_opts(cp,data_opts)
    insp_job.add_ini_opts(cp,data_opts)

  tmplt_job.set_channel(channel)
  insp_job.set_channel(channel)

  # see if we are using calibrated data
  if cp.has_section(data_opts) and cp.has_option(data_opts,'calibrated-data'):
    calibrated = True
    print("we use calibrated data for ", ifo_name) 
  else: calibrated = False

  # prepare the exttrig injection filename
  if ifo_data:
    exttrigStart = ifo_data[0].start()
    exttrigDuration = ifo_data[-1].end()-exttrigStart
    injectionFileTemplate = "HL-INJECTION_%%s-%d-%d.xml" % \
      (exttrigStart, exttrigDuration)

  chunks_analyzed = []
  # loop over the master science segments
  for seg in ifo_data:

    # loop over the master analysis chunks in the science segment
    for chunk in seg:
      done_this_chunk = False

      # now loop over all the data that we need to filter
      for seg_to_do in ifo_to_do:

        # if the current chunk is in one of the segments we need to filter
        if not done_this_chunk and inspiral.overlap_test(chunk,seg_to_do):

          # make sure we only filter the master chunk once
          done_this_chunk = True

          # make sure we have done one and only one datafind for the segment
          if not opts.read_cache:
            if not seg.get_df_node():
              df = pipeline.LSCDataFindNode(df_job)
              if not opts.disable_dag_categories:
                df.set_category('datafind')
              if not opts.disable_dag_priorities:
                df.set_priority(100)
              df.set_observatory(ifo_name[0])
              # add a padding time to the start of the datafind call (but don't change datafind output name)
              if ifo_name == 'G1':
                dfsect = 'geo-data'
              elif ifo_name == 'V1':
                dfsect = 'virgo-data'
              else:
                dfsect = 'ligo-data'
              if cp.has_option(dfsect,ifo_name.lower() + '-datafind-start-padding'):
                padding=cp.get(dfsect,ifo_name.lower()+'-datafind-start-padding')
              else:
                padding=0.
              df.set_start(seg.start(),padding)
              df.set_end(seg.end())
              seg.set_df_node(df)
              if type: df.set_type(type)
              if prev_df and opts.disable_dag_categories:
                df.add_parent(prev_df)
              if opts.datafind: dag.add_node(df)
              prev_df = df
          else:
            prev_df = None  

          # make a template bank job for the master chunk
          bank = inspiral.TmpltBankNode(tmplt_job)
          if not opts.disable_dag_categories:
            bank.set_category('tmpltbank')
          if not opts.disable_dag_priorities:
            bank.set_priority(1)
          bank.set_start(chunk.start())
          bank.set_end(chunk.end())
          bank.set_ifo(ifo_name)
          bank.set_vds_group(ifo_name[0] + str(chunk.start()))
          if not opts.read_cache: bank.set_cache(df.get_output())
          else: bank.set_cache(cp.get('datafind',ifo_name+"-cache"))
          if not calibrated: bank.calibration()
          if opts.datafind: bank.add_parent(df)
          if (opts.template_bank and not GeoBank): dag.add_node(bank)
                  
          # loop over the ext-trig injections
          for inj in range(exttrigInjections[0], exttrigInjections[1]+1):
             if inj>0: 
                # set the appropriate user-tag and the injection filename
                exttrigUserTag = usertag + "_" + str(inj)
                injectionFile = injectionFileTemplate % exttrigUserTag

             # make an inspiral job for the master chunk
             insp = inspiral.InspiralNode(insp_job)
             if not opts.disable_dag_categories:
               insp.set_category('inspiral1')
             if not opts.disable_dag_priorities:
               insp.set_priority(2)
             if inj>0:
               insp.set_injections(injectionFile)
               insp.set_user_tag(exttrigUserTag)
             insp.set_start(chunk.start())
             insp.set_end(chunk.end())
             insp.set_trig_start(chunk.trig_start())
             insp.set_trig_end(chunk.trig_end())
             insp.set_ifo(ifo_name)
             insp.set_ifo_tag("FIRST")
             insp.set_vds_group(ifo_name[0] + str(chunk.start()))
             if not opts.read_cache: insp.set_cache(df.get_output())
             else:  insp.set_cache(cp.get('datafind',ifo_name+"-cache"))
             if not calibrated: insp.calibration()
             if GeoBank:
                insp.set_bank(GeoBank)
             else:
                insp.set_bank(bank.get_output())
             if opts.datafind: insp.add_parent(df)
             if inspinjNode and opts.inspinj: insp.add_parent(inspinjNode) 
             if (opts.template_bank and not GeoBank): insp.add_parent(bank)
             if opts.inspiral: dag.add_node(insp)

             # make an inspiral checkpoint restore job
             if opts.data_checkpoint:  
               insp_job.set_universe("vanilla")
               insp.set_data_checkpoint()
               insp.set_post_script(cp.get('condor','checkpoint-post-script'))
               insp.add_post_script_arg(os.path.join(os.getcwd(),insp.get_checkpoint_image()))
               insp_ckpt = inspiral.InspiralCkptNode(insp_ckpt_job)
               insp_ckpt.set_output(insp.get_output())
               insp_ckpt.set_injections(insp.get_injections())
               insp_ckpt.set_checkpoint_image(insp.get_checkpoint_image())

               if cp.has_option('pipeline','remote-site'):
                 # additional requirements to launch jon on remote pool
                 insp_ckpt_job.set_universe("grid")
                 insp_ckpt.set_grid_start("pegasuslite")
                 insp_ckpt.add_pegasus_profile("condor","grid_resource","condor %s" % cp.get('pipeline','remote-site'))
                 insp_ckpt.add_pegasus_profile("condor","+remote_jobuniverse","5")
                 insp_ckpt.add_pegasus_profile("condor","+remote_requirements","True")
                 insp_ckpt.add_pegasus_profile("condor","+remote_ShouldTransferFiles","True")
                 insp_ckpt.add_pegasus_profile("condor","+remote_WhenToTransferOutput","ON_EXIT")
                 insp_ckpt.add_pegasus_profile("condor","+remote_TransferInputFiles",'"' + insp.get_checkpoint_image() + '"')
                 insp_ckpt.add_pegasus_profile("condor","+remote_PeriodicRelease",'( JobStatus == 5 && HoldReasonCode == 13 && NumSystemHolds < 3 )')
               else:
                 insp_ckpt_job.set_universe("vanilla")

               insp_ckpt.add_parent(insp)
               if opts.inspiral: dag.add_node(insp_ckpt)

               # ensure output is added to list of output files
               output = insp_ckpt.get_output()

               # store this chunk in the list of filtered data
               chunks_analyzed.append(AnalyzedIFOData(chunk,insp_ckpt,inj))
             else:
               # ensure output is added to list of output files
               output = insp.get_output()

               # store this chunk in the list of filtered data
               chunks_analyzed.append(AnalyzedIFOData(chunk,insp,inj))

  return tuple([prev_df,chunks_analyzed])


##############################################################################
# function to do inca on single IFO data  
def single_coinc(ifo_analyzed,ifo_name,ifo_single,single_coinc_job,dag,
  usertag=None):
  """
  Run inca on the single coincident times from each of the IFOs. Since the way 
  we treat all this data is the same, this function is the same for all cases.
  
  ifo_analyzed = the analyzed chunks for the IFO
  ifo_name = the name of the IFO
  ifo_single = the single IFO science segments 
  single_coinc_job = the condor job to do single IFO inca
  dag = the DAG to attach the nodes to
  """

  single_coinc_analyzed = []
  for seg in ifo_single:
    sinca = inspiral.IncaNode(single_coinc_job)
    if not opts.disable_dag_categories:
      sinca.set_category('inca')
    if not opts.disable_dag_priorities:
      sinca.set_priority(100)
    sinca.set_start(seg.start())
    sinca.set_end(seg.end())
    sinca.set_ifo_a(ifo_name)
    sinca.add_var_opt('single-ifo',' ')

    # add all master chunks that overlap with segment to input
    for ifo_done in ifo_analyzed:
      if inspiral.overlap_test(ifo_done.get_chunk(),seg):
        sinca.add_file_arg(ifo_done.get_dag_node().get_output())
        if opts.inspiral: sinca.add_parent(ifo_done.get_dag_node())
    
    if opts.coincidence: dag.add_node(sinca)
    single_coinc_analyzed.append(AnalyzedIFOData(seg,sinca))

  return single_coinc_analyzed
  
##############################################################################
# function to do thinca on coincident IFO data  
def thinca_coinc(ifos,single_data_analyzed,coinc_segs,coinc_job,
    dag,do_coinc,do_insp, exttrigInjections, usertag=None, coire_job=None,
    do_coire=None,ifotag=None,slides=None):
  """
  Run thinca on the coincident times from each of the sets of IFOs. 
  Since the way we treat all this data is the same, this function is the same 
  for all. 
  
  ifos = the ifos we are to analyze
  single_data_analyzed = dictionary of single ifo data analyzed
  coinc_segs = the double coincident IFO science segments 
  coinc_job = the condor job to do two IFO thinca
  dag = the DAG to attach the nodes to
  do_coinc = whether we should add the thinca jobs to the dag
  do_insp  = whether previous inspiral jobs are in the dag
  exttrigInjections = a two-element list specifying the range 
                      of injections for the external trigger case.
  usertag = the usertag to add to the output file name
  coire_job = coire node to join to the output of each thinca node
  do_coire  = whether we should add the coire jobs to the dag
  ifotag = the ifotag to add to the output file name
  slides = number of time slides req'd, used as a boolean
  """
  if ifotag:
    ifotag = "SECOND_" + ifotag
  else:
    ifotag = "FIRST"

  # add time slide arguments if doing a slide
  if slides and cp.has_section('thinca-slide') :
    coinc_job.add_ini_opts(cp,'thinca-slide')
    
  coinc_analyzed = []

  # loop over the exttrig injections
  for inj in range(exttrigInjections[0], exttrigInjections[1]+1):
    if inj>0:
      # set the appropriate user-tag and the injection filename
      exttrigUserTag = usertag + "_" + str(inj)

    for seg in coinc_segs:
      thinca = inspiral.ThincaNode(coinc_job)
      if not opts.disable_dag_categories: 
        thinca.set_category('thinca')
      if not opts.disable_dag_priorities:
        thinca.set_priority(3)
      if cp.has_option('pipeline', 'collapse-thinca'):
        thinca.set_dax_collapse(cp.get('pipeline','collapse-thinca'))
      if inj>0: thinca.set_user_tag(exttrigUserTag)
      thinca.set_start(seg.start())
      thinca.set_end(seg.end())
      thinca.set_ifo_tag(ifotag)
      if slides: thinca.set_num_slides(slides)

      # scroll through ifos, adding the appropriate ones
      for ifo in ifo_list:
        if ifo in ifos:
          thinca.set_ifo(ifo)
                
          # add all IFO A master chunks that overlap with segment to input
          for done in single_data_analyzed[ifo]:
            done_inj = done.get_inj()
            if (done_inj==0 or done_inj==inj) and \
                inspiral.overlap_test(done.get_chunk(),seg):
              # set special user-tag first
              thinca.add_file_arg(done.get_dag_node().get_output())
              if do_insp: thinca.add_parent(done.get_dag_node())

              # set the injection file for this coire
              inj_file = done.get_dag_node().get_injections()

      if do_coinc: dag.add_node(thinca)
      coinc_analyzed.append(AnalyzedIFOData(seg,thinca,inj))

      if do_coire and inj > 0: 
        coire_segments_individually([AnalyzedIFOData(seg,thinca,inj)],
            coire_job, dag, do_coire, do_coinc, ifos, inj_file=inj_file, 
            num_slides=slides, ifotag=ifotag, usertag=exttrigUserTag, 
            inspinjNode=inspinj)

  return coinc_analyzed



##############################################################################
# function to set up trigbank and inspiral jobs for an ifo  
def trig_inspiral(ifo_data,ifo_name,coinc_data,coinc_ifos,trig_job,\
                  insp_veto_job,dag, exttrigInjections, usertag=None,\
                  inspinjNode=None,insp_ckpt_job=None):

  """
  Perform the trigbank and second inspiral steps for the single ifo data.  
  Since the way we treat all this data is the same, this function is the same 
  for all interferometers. Returns the chunks analyzed.

  ifo_data = the master science segs for the IFO
  ifo_name = the name of the ifo
  coinc_data = the coincident science segments and thinca jobs
  coinc_ifos  = the list of ifos in the coincident times
  trig_job = the triggered bank job we should use 
  insp_veto_job = the condor job that we should use to analyze data
  dag = the DAG to attach the nodes to
  exttrigInjections = a two-element list specifying the range 
                      of injections for the external trigger case.
  usertag = the usertag to add to the job names
  """

  # add the veto option
  if cp.has_section('veto-inspiral'): 
    insp_veto_job.add_ini_opts(cp,'veto-inspiral')
  # add ifo specific options
  if cp.has_section(ifo_name.lower() + '-inspiral'):
    insp_veto_job.add_ini_opts(cp,ifo_name.lower() + '-inspiral')
  
  if ifo == 'G1':
    data_opts = 'geo-data'
    channel = cp.get('input','geo-channel')
  elif ifo == 'V1':
    data_opts = 'virgo-data'
    channel = cp.get('input','virgo-channel')
  else:
    data_opts = 'ligo-data'
    channel = cp.get('input','ligo-channel')
 
  if cp.has_section(data_opts): insp_veto_job.add_ini_opts(cp,data_opts)
  insp_veto_job.set_channel(channel)
  
     
  # see if we are using calibrated data
  if cp.has_section(data_opts) and cp.has_option(data_opts,'calibrated-data'):
    calibrated = True
  else: calibrated = False

  # prepare the exttrig injection filename
  if ifo_data:
    exttrigStart = ifo_data[0].start()
    exttrigDuration = ifo_data[-1].end()-exttrigStart
    injectionFileTemplate = "HL-INJECTION_%%s-%d-%d.xml" % \
      (exttrigStart, exttrigDuration)

  trig_chunks_analyzed = []

  # loop over the master science segments
  for seg in ifo_data:

    # loop over the master analysis chunks in the science segment
    for chunk in seg:

      # loop over the exttrig injections
      for inj in range(exttrigInjections[0], exttrigInjections[1]+1):
        done_this_chunk = 0
        
        if inj>0:
          exttrigUserTag=usertag + "_" + str(inj)
          injectionFile = injectionFileTemplate % exttrigUserTag

        # now loop over all the data that we need to filter
        for coinc_done in coinc_data:

          # if the current chunk is in one of the segments we need to filter
          if (inspiral.overlap_test(chunk,coinc_done.get_chunk()) and \
              (inj==0 or coinc_done.get_inj()==inj)):

            if not done_this_chunk:
              # make sure we only filter the master chunk once
              done_this_chunk = 1

              # make a trigbank job for the master chunk
              trigbank = inspiral.TrigbankNode(trig_job)
              if not opts.disable_dag_categories:
                trigbank.set_category('trigbank')
              if not opts.disable_dag_priorities: 
                trigbank.set_priority(4)
              if inj>0: trigbank.set_user_tag(exttrigUserTag)
              trigbank.set_ifo_tag("SECOND_" + coinc_ifos)              
              trigbank.set_output_ifo(ifo_name)
              trigbank.set_input_ifo(ifo_name)
              trigbank.set_start(chunk.start())
              trigbank.set_end(chunk.end())
              if opts.trigbank: dag.add_node(trigbank)
        
              # make an inspiral job for the master chunk
              insp = inspiral.InspiralNode(insp_veto_job)
              if not opts.disable_dag_categories:
                insp.set_category('inspiral2')
              if not opts.disable_dag_priorities:
                insp.set_priority(5)
              if inj>0:
                insp.set_injections(injectionFile)
                insp.set_user_tag(exttrigUserTag)
              insp.set_start(chunk.start())
              insp.set_end(chunk.end())
              insp.set_trig_start(chunk.trig_start())
              insp.set_trig_end(chunk.trig_end())
              insp.set_ifo(ifo_name)
              insp.set_ifo_tag("SECOND_" + coinc_ifos)
              insp.set_vds_group(ifo_name[0] + str(chunk.start()))
              if not opts.read_cache: 
                insp.set_cache(seg.get_df_node().get_output())
              else: insp.set_cache(cp.get('datafind',ifo_name+'-cache')) 
              if not calibrated: insp.calibration()
              insp.set_bank(trigbank.get_output())
              if opts.datafind: insp.add_parent(seg.get_df_node())
              if opts.trigbank: insp.add_parent(trigbank)
              if inspinjNode and opts.inspinj: insp.add_parent(inspinjNode)
              if opts.inspiral_veto: dag.add_node(insp)

              # make an inspiral checkpoint restore job
              if opts.data_checkpoint:  
                insp_veto_job.set_universe("vanilla")
                insp.set_data_checkpoint()
                insp.set_post_script(cp.get('condor','checkpoint-post-script'))
                insp.add_post_script_arg(os.path.join(os.getcwd(),insp.get_checkpoint_image()))
                insp_ckpt = inspiral.InspiralCkptNode(insp_ckpt_job)
                insp_ckpt.set_output(insp.get_output())
                insp_ckpt.set_injections(insp.get_injections())
                insp_ckpt.set_checkpoint_image(insp.get_checkpoint_image())

                if cp.has_option('pipeline','remote-site'):
                  # additional requirements to launch jon on remote pool
                  insp_ckpt_job.set_universe("grid")
                  insp_ckpt.set_grid_start("pegasuslite")
                  insp_ckpt.add_pegasus_profile("condor","grid_resource","condor %s" % cp.get('pipeline','remote-site'))
                  insp_ckpt.add_pegasus_profile("condor","transfer_executable","True")
                  insp_ckpt.add_pegasus_profile("condor","+remote_jobuniverse","5")
                  insp_ckpt.add_pegasus_profile("condor","+remote_requirements","True")
                  insp_ckpt.add_pegasus_profile("condor","+remote_ShouldTransferFiles","True")
                  insp_ckpt.add_pegasus_profile("condor","+remote_WhenToTransferOutput","ON_EXIT")
                  insp_ckpt.add_pegasus_profile("condor","+remote_TransferInputFiles",'"' + insp.get_checkpoint_image() + '"')
                  insp_ckpt.add_pegasus_profile("condor","+remote_PeriodicRelease",'( JobStatus == 5 && HoldReasonCode == 13 && NumSystemHolds < 3 )')
                else:
                  insp_ckpt_job.set_universe("vanilla")

                insp_ckpt.add_parent(insp)
                if opts.inspiral_veto: dag.add_node(insp_ckpt)

                # ensure output is added to list of output files
                output = insp_ckpt.get_output()

                # store this chunk in the list of filtered data (ckpt version)
                trig_chunks_analyzed.append(AnalyzedIFOData(chunk,insp_ckpt,inj))
              else:
                # ensure output is added to list of output files
                output = insp.get_output()

                # store this chunk in the list of filtered data (ckpt version)
                trig_chunks_analyzed.append(AnalyzedIFOData(chunk,insp,inj))

            # jobs already set up, just add the var_arg and parent to trigbank
            trigbank.add_file_arg(coinc_done.get_dag_node().get_output())
            if opts.coincidence: 
                trigbank.add_parent(coinc_done.get_dag_node())

  return trig_chunks_analyzed


##############################################################################
# function to TD follow-up tmpltbank and inspiral jobs for an ifo  
def follow_up(ifo_data,ifo_name,coinc_data,coinc_ifos,tmplt_follow_job,\
                  insp_veto_job,dag,approximant=None,\
                  usertag=None,banks=None):

  """
  Perform the tmpltbank and second inspiral steps for the single ifo data.  
  Since the way we treat all this data is the same, this function is the same 
  for all interferometers. Returns the chunks analyzed.

  ifo_data = the master science segs for the IFO
  ifo_name = the name of the ifo
  coinc_data = the coincident science segments and thinca jobs
  coinc_ifos  = the list of ifos in the coincident times
  trig_job = the triggered bank job we should use 
  insp_veto_job = the condor job that we should use to analyze data
  dag = the DAG to attach the nodes to
  usertag = the usertag to add to the job names
  """

  # Store the template bank jobs
  bankList = {} 

  # add the veto option
  if cp.has_section('veto-inspiral'):
    insp_veto_job.add_ini_opts(cp,'veto-inspiral')
  # add ifo specific options
  if cp.has_section(ifo_name.lower() + '-inspiral'):
    insp_veto_job.add_ini_opts(cp,ifo_name.lower() + '-inspiral')
  
  if ifo == 'G1':
    data_opts = 'geo-data'
    channel = cp.get('input','geo-channel')
  elif ifo == 'V1':
    data_opts = 'virgo-data'  
    channel = cp.get('input','virgo-channel')
  else:
    data_opts = 'ligo-data'
    channel = cp.get('input','ligo-channel')
 
  if cp.has_section(data_opts): 
    insp_veto_job.add_ini_opts(cp,data_opts)
    tmplt_follow_job.add_ini_opts(cp, data_opts)
  if cp.has_section('tmpltbank-2'):
    tmplt_follow_job.add_ini_opts(cp, 'tmpltbank-2')
  insp_veto_job.set_channel(channel)
  tmplt_follow_job.set_channel(channel)
  
     
  # see if we are using calibrated data
  if cp.has_section(data_opts) and cp.has_option(data_opts,'calibrated-data'):
    calibrated = True
  else: calibrated = False

  follow_chunks_analyzed = []
  # loop over the master science segments
  for seg in ifo_data:

    # loop over the master analysis chunks in the science segment
    for chunk in seg:
      done_this_chunk = 0

      # now loop over all the data that we need to filter
      for coinc_done in coinc_data:

        # if the current chunk is in one of the segments we need to filter
        if (inspiral.overlap_test(chunk,coinc_done.get_chunk()) ):

          if not done_this_chunk:
            # make sure we only filter the master chunk once
            done_this_chunk = 1

            bank = None
            # We only need to run the template bank for 1 approximant
            if not banks:
              # make a template bank job for the master chunk
              bank = inspiral.TmpltBankNode(tmplt_follow_job)
              if not opts.disable_dag_categories:
                bank.set_category('tmpltbank')
              if not opts.disable_dag_priorities:
                bank.set_priority(1)
              bank.set_start(chunk.start())
              bank.set_end(chunk.end())
              bank.set_ifo(ifo_name)
              bank.set_ifo_tag(coinc_ifos)
              bank.set_vds_group(ifo_name[0] + str(chunk.start()))
              if not opts.read_cache: 
                  bank.set_cache(seg.get_df_node().get_output())
              else: bank.set_cache(cp.get('datafind',ifo_name+"-cache"))
              if not calibrated: bank.calibration()
              if opts.datafind: bank.add_parent(seg.get_df_node())
              if opts.td_follow_bank: dag.add_node(bank)

            else:
              # Use the bank generated for a previous approximant
              bank = banks[chunk.start()]

            bankList[chunk.start()] = bank

            # make an inspiral job for the master chunk
            insp = inspiral.InspiralNode(insp_veto_job)
            if not opts.disable_dag_categories:
              insp.set_category('inspiral2')
            if not opts.disable_dag_priorities:
              insp.set_priority(5)
            insp.set_start(chunk.start())
            insp.set_end(chunk.end())
            insp.set_trig_start(chunk.trig_start())
            insp.set_trig_end(chunk.trig_end())
            insp.set_ifo(ifo_name)
            insp.set_ifo_tag(coinc_ifos+"_"+approximant)
            insp.set_vds_group(ifo_name[0] + str(chunk.start()))
            if not opts.read_cache: 
              insp.set_cache(seg.get_df_node().get_output())
            else: insp.set_cache(cp.get('datafind',ifo_name+'-cache')) 
            if not calibrated: insp.calibration()
            insp.set_bank(bank.get_output())
            if approximant: insp.add_var_opt('approximant', approximant)
            # XXX: ensure output is added to list of output files
            output = insp.get_output()
            if opts.datafind: insp.add_parent(seg.get_df_node())
            if opts.td_follow_bank: insp.add_parent(bank)
            if opts.td_follow_inspiral: dag.add_node(insp)

            # store this chunk in the list of filtered data -- old version
            follow_chunks_analyzed.append(AnalyzedIFOData(chunk,insp))

            # jobs already set up, just add the var_arg and parent to trigbank
            # and inspiral-veto 
          if 'macrotdfollowup' in insp.get_opts().keys():
            insp.add_var_opt('td-follow-up', 
                insp.get_opts()['macrotdfollowup'] +
                ' ' + coinc_done.get_dag_node().get_output())
          else:
            insp.add_var_opt('td-follow-up',
                coinc_done.get_dag_node().get_output())

          if not banks: 
            if 'macrotdfollowup' in bank.get_opts().keys():
              bank.add_var_opt('td-follow-up', 
                  bank.get_opts()['macrotdfollowup'] +
                  ' ' + coinc_done.get_dag_node().get_output())
            else:
              bank.add_var_opt('td-follow-up',
                  coinc_done.get_dag_node().get_output())
            if opts.coincidence: bank.add_parent(coinc_done.get_dag_node())

  return follow_chunks_analyzed, bankList



##############################################################################
# function to sweep up data using sire
def sire_segments(segments, sire_job, dag, do_sire, do_input, ifo, 
  sire_start, sire_end, inj_file=None, ifotag=None, usertag=None, 
  inspinjNode=None):
  """
  Do a sire to sweep up all the triggers of a specific kind
  
  segments  = list of segments to use
  sire_job  = sire job to use for the analysis
  dag       = name of the dag
  do_sire   = whether the sired jobs are to be run by the dag
  do_input  = whether the files to be sired were produced by this dag
  ifo       = which ifo triggers to keep
  sire_start= start time for sire job
  sire_end  = end time for sire job
  inj_file  = name of injection file
  ifotag    = ifotag used in naming the file
  usertag   = the usertag to add to the output file name
  inspinjNode = add an inspinj node as a parent (default false)
  """
  # only write the input file if we're doing the job
  sire = inspiral.SireNode(sire_job)
  if not opts.disable_dag_categories:
    sire.set_category('sire')
  if not opts.disable_dag_priorities:
    sire.set_priority(100)
  if cp.has_option('pipeline', 'collapse-sire'):
    sire.set_dax_collapse(cp.get('pipeline','collapse-sire'))

  # set the options:
  sire.set_ifo(ifo)
  if inj_file: 
    sire.set_inj_file(inj_file)
  if ifotag: sire.set_ifo_tag(ifotag)
  if usertag: sire.set_user_tag(usertag)

  # set the segment
  sire.set_start(sire_start)
  sire.set_end(sire_end)

  if inspinjNode and opts.inspinj: sire.add_parent(inspinjNode)

  for seg in segments:
    # set the input files for sire
    if do_input: sire.add_parent(seg.get_dag_node())

    if seg.get_dag_node().get_output():
      output = seg.get_dag_node().get_output()
    else: 
      output = seg.get_dag_node().get_output_a()

    sire.add_file_arg(output)

  if do_sire: dag.add_node(sire) 

  return sire


##############################################################################
# function to sweep up data using sire
def sire_segments_individually(segments, sire_job, dag, do_sire, do_input, ifo,
  inj_file=None, ifotag=None, usertag=None, inspinjNode=None):
  """
  Do a sire to sweep up all the triggers of a specific kind
  
  segments  = list of thinca segments to use
  sire_job  = sire job to use for the analysis
  dag       = name of the dag
  do_sire   = whether the sired jobs are to be run by the dag
  do_input  = whether the files to be sired were produced by this dag
  ifo       = which ifo triggers to keep
  inj_file  = name of injection file
  ifotag    = ifotag used in naming the file
  usertag   = the usertag to add to the output file name
  """
  sire_analyzed = []
  for seg in segments:
    sire_node = sire_segments([seg], sire_job, dag, do_sire, do_input, ifo, \
        seg.get_chunk().start(), seg.get_chunk().end(), inj_file, ifotag, 
        usertag,inspinjNode=inspinjNode)

    sire_analyzed.append(AnalyzedIFOData(seg.get_chunk(),sire_node))

  return sire_analyzed

##############################################################################
# function to sweep up data using coire  
def coire_segments(segments, coire_job, dag, do_coire, do_input, ifos,
  coire_start, coire_end, inj_file=None, num_slides=0, ifotag=None, 
  usertag=None, inspinjNode=None):
  """
  Do a coire to sweep up all the triggers of a specific kind
  
  segments  = list of thinca segments to use
  coire_job = coire job to use for the analysis
  dag       = name of the dag
  do_coire  = whether the coire node generated should be added to the dag
  do_input  = whether the files to be coired were produced by this dag
  ifos      = which ifos were analyzed
  inj_file  = name of injection file
  cluster   = whether we should add the cluster arguments
  num_slides = the number of time slides
  usertag   = the usertag to add to the output file name
  inspinjNode = add an inspinj node as a parent (default false)
  """
  # only write the input file if we're running the coires
  coire = inspiral.CoireNode(coire_job)
  if not opts.disable_dag_categories:
    coire.set_category('coire')
  if not opts.disable_dag_priorities:
    coire.set_priority(100)
  if cp.has_option('pipeline', 'collapse-coire'):
    coire.set_dax_collapse(cp.get('pipeline','collapse-coire'))

  # set the options:
  coire.set_ifos(ifos)
  if inj_file:
    coire.set_inj_file(inj_file)
  if ifotag is not None: coire.set_ifo_tag(ifotag)
  if num_slides: coire.set_slides(num_slides)
  if usertag: coire.set_user_tag(usertag)
 
  # set the segment
  coire.set_start(coire_start)
  coire.set_end(coire_end)

  if inspinjNode and opts.inspinj: coire.add_parent(inspinjNode)

  for seg in segments:
    # set the glob file for coire 
    if do_input: coire.add_parent(seg.get_dag_node())

    if seg.get_dag_node().get_output():
      output = seg.get_dag_node().get_output()
    else: 
      output = seg.get_dag_node().get_output_a()

    coire.add_file_arg(output)

  coire.finalize()
  if do_coire: dag.add_node(coire)

  return coire

##############################################################################
# function to sweep up data using coire  
def coire_segments_individually(segments, coire_job, dag, do_coire, do_input, 
  ifos, inj_file=None, num_slides=0, ifotag=None, usertag=None,inspinjNode=None):
  """
  Do a coire to sweep up all the triggers of a specific kind
  
  segments  = list of thinca segments to use
  coire_job = coire job to use for the analysis
  dag       = name of the dag
  do_coire  = whether the coire nodes should be added to the dag
  do_input  = whether the files to be coired were produced by this dag
  ifos      = which ifos were analyzed
  inj_file  = name of injection file
  cluster   = whether we should add the cluster arguments
  num_slides = the number of time slides
  usertag   = the usertag to add to the output file name
  """

  coire_analyzed = []
  for seg in segments:
    coire_node = coire_segments([seg], coire_job, dag, do_coire, do_input, ifos,
        seg.get_chunk().start(), seg.get_chunk().end(), inj_file = inj_file, 
        num_slides = num_slides, ifotag = ifotag, usertag = usertag,
        inspinjNode=inspinjNode)
    coire_analyzed.append(AnalyzedIFOData(seg.get_chunk(),coire_node))

  return coire_analyzed

##############################################################################
# function to set up the cohbank jobs for a network of ifos
def coherent_bank(ifos,cb_job,dag,coinc_nodes,exttrigInjections,
  num_slides,doSlide=None,usertag=None):
  """
  Create template banks from coincident trigger parameters for
  coherent analysis.

  ifos                = ordered set of ifo identifiers
  coinc_nodes         = the coincident nodes for this ifo combination
  cb_job              = coherent bank job
  dag                 = name of the dag
  """

  if doSlide:
    ifotag = 'SLIDE_COHERENT_' + ifos
  else:
    ifotag = 'COHERENT_' + ifos

  cohbank_analyzed = []
  inj = None

  for coinc_done in coinc_nodes:
    # set up cohbank
    cohbank = inspiral.CohBankNode( cb_job )
    if not opts.disable_dag_categories:
      cohbank.set_category('cohbank')
    if not opts.disable_dag_priorities:
      cohbank.set_priority(6)
    cohbank.set_start(coinc_done.get_chunk().start())
    cohbank.set_end(coinc_done.get_chunk().end())
    if num_slides: cohbank.set_num_slides(num_slides)
    if usertag is not None:
      ifousertag = ifotag + "_" + usertag
      if coinc_done.get_inj():
        inj = coinc_done.get_inj()
        ifousertag = ifousertag + "_" + str(coinc_done.get_inj())
      cohbank.set_user_tag(ifousertag)
    else:
      if coinc_done.get_inj():
        inj = coinc_done.get_inj()
        ifotag = ifotag + "_" + str(coinc_done.get_inj())
      cohbank.set_user_tag(ifotag)
    cohbank.set_ifos(ifos)
    cohbank.add_file_arg(coinc_done.get_dag_node().get_output())
    if opts.second_coinc: cohbank.add_parent(coinc_done.get_dag_node())
    if opts.coherent_bank: dag.add_node(cohbank)
    cohbank_analyzed.append(AnalyzedIFOData(coinc_done.get_chunk(),cohbank,inj))

  return cohbank_analyzed

##############################################################################
# function to set up trigbank jobs for coherent analysis with a network of ifos
def coherent_trig_inspiral(ifo_data,ifo_name,cohbank_data,ifos,trig_coh_jobs,
  insp_coh_jobs,dag,calibrated,exttrigInjections,doSlide=None,usertag=None,
  inspinjNode = None):
  """
  Perform the trigbank and third inspiral steps for the single ifo triggers
  appearing in the cohbank templates.

  ifo_data      = the master science segs for the IFO
  ifo_name      = the name of the ifo
  cohbank_data  = the coincident science segments and cohbank jobs
  ifos          = the list of ifos in the coincident times
  trig_coh_jobs = the triggered coherent bank job we should use
  insp_coh_jobs = inspiral jobs to produce c-data
  dag           = the DAG to attach the nodes to

  usertag = the usertag to add to the job names
  inspinjNode = the inspinj node to be added as a parent to inspirals
  """

  if doSlide:
    ifotag = 'SLIDE_COHERENT_' + ifos
  else:
    ifotag = 'COHERENT_' + ifos

  # prepare the exttrig injection filename
  if ifo_data:
    exttrigStart = ifo_data[0].start()
    exttrigDuration = ifo_data[-1].end()-exttrigStart
    injectionFileTemplate = "HL-INJECTION_%%s-%d-%d.xml" % \
      (exttrigStart, exttrigDuration)

  trig_coh_chunks_analyzed = []

  # loop over the master science segments
  for seg in ifo_data:

    # loop over the master analysis chunks in the science segment
    for chunk in seg:

      # loop over the exttrig injections
      for inj in range(exttrigInjections[0], exttrigInjections[1]+1):
        done_this_chunk = 0

        if inj>0:
          exttrigUserTag=usertag + "_" + str(inj)
          injectionFile = injectionFileTemplate % exttrigUserTag

        # now loop over all the data that we need to filter
        for coinc_done in cohbank_data:

          # if the current chunk is in one of the segments we need to filter
          if inspiral.overlap_test(chunk,coinc_done.get_chunk()) and \
              (inj==0 or coinc_done.get_inj()==inj):

            if not done_this_chunk:
              # make sure we only filter the master chunk once
              done_this_chunk = 1

              # make a trigbank job for the master chunk
              trigbank = inspiral.TrigbankNode(trig_coh_jobs)
              if not opts.disable_dag_categories:
                trigbank.set_category('trigbank2')
              if not opts.disable_dag_priorities:
                trigbank.set_priority(9)
              if inj>0: trigbank.set_user_tag(exttrigUserTag)
              trigbank.set_ifo_tag(ifotag)
              trigbank.set_output_ifo(ifo_name)
              trigbank.set_input_ifo(ifo_name)
              trigbank.set_start(chunk.start())
              trigbank.set_end(chunk.end())
              if opts.coherent_bank: dag.add_node(trigbank)

              # make an inspiral job for the master chunk
              insp = inspiral.InspiralNode(insp_coh_jobs)
              if not opts.disable_dag_categories:
                insp.set_category('inspiral3')
              if not opts.disable_dag_priorities:
                insp.set_priority(10)
              if inj>0:
                insp.set_injections(injectionFile)
                insp.set_user_tag(exttrigUserTag)
              insp.set_start(chunk.start())
              insp.set_end(chunk.end())
              insp.set_ifo(ifo_name)
              insp.set_ifo_tag(ifotag)
              insp.add_var_opt('write-cdata','')
              insp.set_vds_group(ifo_name[0] + str(chunk.start()))
              if not opts.read_cache: insp.set_cache(seg.get_df_node().get_output())
              else: insp.set_cache(cp.get('datafind',ifo_name+'-cache'))
              if not calibrated: insp.calibration()
              insp.set_bank(trigbank.get_output())
              # XXX: ensure output is added to list of output files
              output = insp.get_output()
              if opts.datafind: insp.add_parent(seg.get_df_node())
              if inspinjNode and opts.inspinj: insp.add_parent(inspinjNode)
              if opts.coherent_bank: insp.add_parent(trigbank)
              if opts.coherent_bank: dag.add_node(insp)

              # store this chunk in the list of filtered data
              trig_coh_chunks_analyzed.append(AnalyzedIFOData(chunk,insp,inj))

            # jobs already set up, just add the var_arg and parent to trigbank
            trigbank.add_file_arg(coinc_done.get_dag_node().get_output())
            if opts.coherent_bank:
              trigbank.add_parent(coinc_done.get_dag_node())

  return trig_coh_chunks_analyzed

##############################################################################
# function to set up the network trigger bank jobs for a set of ifos
def coherent_inspiral_bank(ifos,cohinspbank_job,dag,coinc_segs,single_data_analyzed,
    exttrigInjections,num_slides,doSlide=None,usertag=None):
  """
  Create template banks from inspiral-coherent trigger parameters for
  coherent analysis.

  ifos                = ordered set of ifo identifiers
  coinc_nodes         = the coincident nodes for this ifo combination
  cohinspbank_job     = network trigger bank job
  dag                 = name of the dag
  """

  if doSlide:
    ifotag = 'SLIDE_COHERENT_' + ifos
  else:
    ifotag = 'COHERENT_' + ifos

  cohinspbank_analyzed = []
  inj = None

  # loop over the exttrig injections
  for inj in range(exttrigInjections[0], exttrigInjections[1]+1):
    if inj>0:
      # set the appropriate user-tag and the injection filename
      if usertag is not None:
        exttrigUserTag  = ifotag + "_" + usertag + "_" + str(inj)
      else:
        exttrigUserTag  = ifotag + "_" + str(inj)
    else:
      if usertag is not None:
        exttrigUserTag  = ifotag + "_" + usertag
      else:
        exttrigUserTag  = ifotag

    for seg in coinc_segs:
      cohinspbank = inspiral.CohInspBankNode(cohinspbank_job)
      if not opts.disable_dag_categories:
        cohinspbank.set_category('cohinspbank')
      if not opts.disable_dag_priorities:
        cohinspbank.set_priority(6)
      cohinspbank.set_start(seg.start())
      cohinspbank.set_end(seg.end())
      cohinspbank.set_user_tag(exttrigUserTag)
      cohinspbank.set_ifos(ifos)
      if num_slides: cohinspbank.set_num_slides(num_slides)

      # scroll through ifos, adding the appropriate ones
      for ifo in ifo_list:
        if ifo in ifos:

          # add all IFO A master chunks that overlap with segment to input
          for done in single_data_analyzed[ifo]:
            done_inj = done.get_inj()
            if (done_inj==0 or done_inj==inj) and \
                inspiral.overlap_test(done.get_chunk(),seg):
              # set special user-tag first
              cohinspbank.add_file_arg(done.get_dag_node().get_output())
              if opts.coherent_bank:
                cohinspbank.add_parent(done.get_dag_node())

              # set the injection file for this coire
              inj_file = done.get_dag_node().get_injections()

      if opts.coherent_bank:
        dag.add_node(cohinspbank)
      cohinspbank_analyzed.append(AnalyzedIFOData(seg,cohinspbank,inj))

  return cohinspbank_analyzed

##############################################################################
# function to set up the cohinspiral jobs for a network of ifos
def coherent_analysis(ifos,single_ifo_segments,cohinsp_job,dag,cohbank_data,
  cohtrig_data,exttrigInjections,doSlide=None,usertag=None,inspinjNode = None):
  """
  Do a coherent analysis

  ifos                = ordered set of ifo identifiers
  single_ifo_segments = dictionary of data segments for the single ifos
  cohinsp_job         = the coherent_inspiral job
  dag                 = name of the dag
  cohbank_data        = the coicident science segments and cohbank jobs
  cohtrig_data        = dictionary of single ifo (inspiral3) data analyzed
  usertag = the usertag to add to the job names
  inspinjNode = the inspinj node to be added as a parent to inspirals
  """

  if doSlide:
    ifotag = 'SLIDE_COHERENT_' + ifos
    if cp.has_section('thinca-slide') :
      cohinsp_job.add_ini_opts(cp,'thinca-slide')
  else:
    ifotag = 'COHERENT_' + ifos

  cohinsp_analyzed = []

  inj = 0

  for cohbank_done in cohbank_data:

    # set up cohinsp
    cohinsp = inspiral.ChiaNode(cohinsp_job)
    if not opts.disable_dag_categories:
      cohinsp.set_category('cohinsp')
    if not opts.disable_dag_priorities:
      cohinsp.set_priority(7)
    cohinsp.set_bank(cohbank_done.get_dag_node().get_output())
    cohinsp.set_start(cohbank_done.get_chunk().start())
    cohinsp.set_end(cohbank_done.get_chunk().end())
    cohinsp.set_ifo_tag(ifos)
    if usertag is not None:
      ifousertag = ifotag + "_" + usertag
    else: ifousertag = ifotag
    if cohbank_done.get_inj():
      inj = cohbank_done.get_inj()
      cohinsptag = ifousertag + "_" + str(cohbank_done.get_inj())
    else: cohinsptag = ifousertag
    cohinsp.set_user_tag(cohinsptag)
    cohinsp.add_var_opt('ifo-tag',ifos)

    # add the coherent bank job as a parent
    if opts.coherent_bank:
      cohinsp.add_parent(cohbank_done.get_dag_node())

    # set up trigbank and inspiral
    for ifo in ifo_list:
      if ifo in ifos:
        # set up the frjoin job
        frjoin = inspiral.FrJoinNode(frjoin_job)
        if not opts.disable_dag_categories:
          frjoin.set_category('frjoin')
        if not opts.disable_dag_priorities:
          frjoin.set_priority(8)
        outputName = ifo + '-INSPIRAL_' + ifotag
        if usertag: outputName += '_' + usertag
        if cohbank_done.get_inj(): outputName += '_' +str(cohbank_done.get_inj())
        outputName += '-' + str(cohbank_done.get_chunk().start()) + \
            '-' + str(cohbank_done.get_chunk().end() - \
            cohbank_done.get_chunk().start()) + '.gwf'
        frjoin.set_output(outputName)

        for done in cohtrig_data[ifo]:
          if inspiral.overlap_test(cohbank_done.get_chunk(),done.get_chunk()) \
              and ( inj == 0 or (cohbank_done.get_inj() == done.get_inj())):
            frjoin.add_file_arg( done.get_dag_node().get_froutput() )
            if opts.coherent_bank: frjoin.add_parent(done.get_dag_node())

        cohinsp.add_file_opt( ifo.lower() + '-framefile', frjoin.get_output() )
        if opts.coherent_bank: dag.add_node(frjoin)
        if opts.coherent_inspiral and opts.coherent_bank:
           cohinsp.add_parent(frjoin)

    if inspinjNode and opts.inspinj: cohinsp.add_parent(inspinjNode)

    if opts.coherent_inspiral: dag.add_node(cohinsp)
    cohinsp_analyzed.append(AnalyzedIFOData(cohbank_done.get_chunk(),cohinsp))

  return cohinsp_analyzed


##############################################################################
# function to sweep up data using cohire
def cohire_segments(segments, cohire_job, dag, do_cohire, do_input, ifos,
  cohire_start, cohire_end, num_slides, inj_file=None, ifotag=None,
  usertag=None, inspinjNode=None):
  """
  Do a cohire to sweep up all the triggers of a specific kind

  segments  = list of thinca segments to use
  cohire_job = cohire job to use for the analysis
  dag       = name of the dag
  do_cohire  = whether the cohire node generated should be added to the dag
  do_input  = whether the files to be cohired were produced by this dag
  ifos      = which ifos were analyzed
  inj_file  = name of injection file
  cluster   = whether we should add the cluster arguments
  num_slides = the number of time slides
  usertag   = the usertag to add to the output file name
  """
  # only write the input file if we're running the cohires
  cohire = inspiral.CohireNode(cohire_job)
  if not opts.disable_dag_categories:
    cohire.set_category('cohire')
  if not opts.disable_dag_priorities:
    cohire.set_priority(100)

  # set the options:
  cohire.set_ifos(ifos)
  if inj_file:
    cohire.set_inj_file(inj_file)
  if ifotag is not None: cohire.set_ifo_tag(ifotag)
  if num_slides: cohire.set_slides(num_slides)
  if usertag: cohire.set_user_tag(usertag)

  # set the segment
  cohire.set_start(cohire_start)
  cohire.set_end(cohire_end)

  if inspinjNode and opts.inspinj: cohire.add_parent(inspinjNode)

  for seg in segments:
    # set the glob file for cohire
    if do_input: cohire.add_parent(seg.get_dag_node())

    if seg.get_dag_node().get_output():
      output = seg.get_dag_node().get_output()
    else:
      output = seg.get_dag_node().get_output_a()

    cohire.add_file_arg(output)

  cohire.finalize()
  if do_cohire: dag.add_node(cohire)

  return cohire


##############################################################################
# function to sweep up data using cohire
def cohire_segments_individually(segments, cohire_job, dag, do_cohire,
  do_input, ifos, num_slides, inj_file=None, ifotag=None, usertag=None,
  inspinjNode=None):
  """
  Do a cohire to sweep up all the triggers of a specific kind

  segments  = list of cohinsp segments to use
  cohire_job = cohire job to use for the analysis
  dag       = name of the dag
  do_cohire  = whether the cohire nodes should be added to the dag
  do_input  = whether the files to be cohired were produced by this dag
  ifos      = which ifos were analyzed
  inj_file  = name of injection file
  cluster   = whether we should add the cluster arguments
  num_slides = the number of time slides
  usertag   = the usertag to add to the output file name
  """

  cohire_analyzed = []
  for seg in segments:
    cohire_node = cohire_segments([seg], cohire_job, dag, do_cohire, do_input,
        ifos, seg.get_chunk().start(), seg.get_chunk().end(), num_slides,
        inj_file = inj_file, ifotag = ifotag, usertag = usertag,
        inspinjNode=inspinjNode)
    cohire_analyzed.append(AnalyzedIFOData(seg.get_chunk(),cohire_node))

  return cohire_analyzed


##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################
usage = """usage: %prog [options] 
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")
    
parser.add_option("-u", "--user-tag",action="store",type="string",\
    default=None,metavar=" USERTAG",\
    help="tag the jobs with USERTAG (overrides value in ini file)")

parser.add_option("-g", "--g1-data",action="store_true",default=False,\
    help="analyze g1 data")
parser.add_option("-a", "--h1-data",action="store_true",default=False,\
    help="analyze h1 data")
parser.add_option("-b", "--h2-data",action="store_true",default=False,\
    help="analyze h2 data")
parser.add_option("-l", "--l1-data",action="store_true",default=False,\
    help="analyze l1 data")
parser.add_option("-n", "--v1-data",action="store_true",default=False,\
    help="analyze v1 data")

parser.add_option("-S", "--one-ifo",action="store_true",default=False,\
    help="analyze single ifo data (not usable for GEO)")
parser.add_option("-D", "--two-ifo",action="store_true",default=False,\
    help="analyze two interferometer data")
parser.add_option("-T", "--three-ifo",action="store_true",default=False,\
    help="analyze three interferometer data")
parser.add_option("-Q", "--four-ifo",action="store_true",default=False,\
    help="analyze four intereferometer data")
  
parser.add_option("-A", "--analyze-all",action="store_true",default=False,\
    help="analyze all ifos and all data (over-rides above)" \
    "(option not available since there are 5 instruments and the code " \
    "only supports quadruple coincidence)")


parser.add_option("-d", "--datafind",action="store_true",default=False,\
    help="run LSCdataFind to create frame cache files")
parser.add_option("-I", "--inspinj",action="store_true",default=False,\
    help="run lalapps_inspinj to generate injection files")
parser.add_option("-t", "--template-bank",action="store_true",default=False,\
    help="run lalapps_tmpltbank to generate template banks")
parser.add_option("-i", "--inspiral" ,action="store_true",default=False,\
    help="run lalapps_inspiral (or lalapps_ring) to generate triggers")
parser.add_option("-s", "--sire-inspiral",action="store_true",default=False,\
    help="do sires to sweep up triggers")
parser.add_option("-r", "--ringdown" ,action="store_true",default=False,\
    help="run lalapps_ring insted of lalapps_inspiral")
parser.add_option("-c", "--coincidence",action="store_true",default=False,\
    help="run lalapps_thinca to test for coincidence")
parser.add_option("-e", "--coire-coincidence",action="store_true",
    default=False, help="do coires to sweep up triggers")
parser.add_option("-U", "--td-follow-bank", action="store_true",default=False,\
    help="Run tmpltbank for TD follow-up")
parser.add_option("-B", "--trigbank",action="store_true",default=False,\
    help="run lalapps_trigbank for banks of coinc triggers")
parser.add_option("-V", "--inspiral-veto",action="store_true",default=False,\
    help="run lalapps_inspiral with vetos")
parser.add_option("-m", "--sire-inspiral-veto",action="store_true",
    default=False, help="do sires to sweep up triggers")
parser.add_option("-W", "--td-follow-inspiral",action="store_true",\
    default=False,help="run lalapps_inspiral for TD follow-up")
parser.add_option("-C", "--second-coinc" ,action="store_true",default=False,\
    help="run lalapps_thinca on the inspiral veto triggers")
parser.add_option("-E", "--coire-second-coinc",action="store_true",
    default=False, help="do coires to sweep up triggers")
parser.add_option("-F", "--sire-second-coinc",action="store_true",
    default=False, help="do sires to sweep up triggers")
parser.add_option("-j", "--coherent-bank",action="store_true",default=False,\
    help="run lalapps_coherentbank to make coherent bank")
parser.add_option("-k", "--coherent-inspiral",action="store_true",
    default=False,help="run lalapps_coherent_inspiral for coherent analysis")
parser.add_option("-q", "--cohire",action="store_true",
    default=False, help="do cohires to sweep up coherent triggers")
parser.add_option("-J", "--summary-coherent-inspiral-triggers",action="store_true", \
    default=False, \
    help="Produce summary coherent triggers from run.  This will produce "
    "cohire files containing coherent triggers from the whole run. "
    "Use with care, as files can become very large.")
parser.add_option("-Y", "--summary-inspiral-triggers",action="store_true", \
    default=False, \
    help="Produce summary triggers from run.  This will produce "
    "sire (first inspiral) files containing triggers from the whole run. "
    "Use with care, as files can become very large.")
parser.add_option("-X", "--summary-coinc-triggers",action="store_true", \
    default=False, \
    help="Produce summary triggers from run.  This will produce "
    "coire (second coincidence) files containing triggers from the whole run. "
    "Use with care, as files can become very large.")
parser.add_option("-Z", "--summary-single-ifo-triggers",action="store_true", \
    default=False, \
    help="Produce summary triggers from run.  This will produce "
    "sire (second coincidence) files containing triggers from the whole run. "
    "Use with care, as files can become very large.")

parser.add_option("-R", "--read-cache",action="store_true",default=False,\
    help="read cache file from ini-file (if LSCDataFind is broken)")

parser.add_option("-P", "--priority",action="store",type="int",\
    metavar=" PRIO",help="run jobs with condor priority PRIO")

parser.add_option("", "--disable-dag-categories",action="store_true",
    default=False,help="disable the internal dag category maxjobs")
parser.add_option("", "--disable-dag-priorities",action="store_true",
    default=False,help="disable the depth first priorities")

parser.add_option("", "--noop-inspinj", action="store_true", default=False,
    help="create a DAG with fake (no-op) inspinj jobs")
 
parser.add_option("-f", "--config-file",action="store",type="string",\
    metavar=" FILE",help="use configuration file FILE")

parser.add_option("-p", "--log-path",action="store",type="string",\
    metavar=" PATH",help="directory to write condor log file")

parser.add_option("-o", "--output-segs",action="store_true",default=False,\
    help="output the segment lists of analyzed data")

parser.add_option("-x","--dax", action="store_true", default=False,\
    help="create a dax instead of a dag")

parser.add_option("-w", "--write-script", action="store_true", default=False,
      help="write the workflow to a locally executable script")

parser.add_option("--summary-first-coinc-triggers", action='store_true', default=False,\
    help="coire (first coincidence) files containing triggers from the whole run. "
    "Use with care, as files can become very large.")

parser.add_option("--data-checkpoint", action="store_true", default=False,\
    help="checkpoint the inspiral code")

parser.add_option("--use-gpus", action="store_true", default=False,\
    help="run inspiral jobs on GPU nodes")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()


#################################
# if --version flagged
if opts.version:
  print("$Id$")
  sys.exit(0)

#################################
# Sanity check of input arguments
if not opts.config_file:
  print("No configuration file specified.", file=sys.stderr)
  print("Use --config-file FILE to specify location.", file=sys.stderr)
  sys.exit(1)

if not opts.log_path:
  print("No log file path specified.", file=sys.stderr)
  print("Use --log-path PATH to specify a location.", file=sys.stderr)
  sys.exit(1)

if not opts.g1_data and not opts.h1_data and not opts.h2_data and \
    not opts.l1_data and not opts.v1_data and not opts.analyze_all:
  print("No ifos specified.  Please specify at least one of", file=sys.stderr)
  print("--g1-data, --h1-data, --h2-data, --l1-data, --v1-data", file=sys.stderr)
  print("or use --analyze-all to analyze all ifos all data", file=sys.stderr)
  sys.exit(1)
elif opts.analyze_all:
  print("The --analyze-all flag is currently not available.", file=sys.stderr)
  print("The code supports quadruple coincidence, so you can", file=sys.stderr)
  print("choose at most four instruments to analyze.", file=sys.stderr)
  sys.exit(1)

if opts.g1_data and opts.h1_data and opts.h2_data and opts.l1_data \
    and opts.v1_data:
  print("Too many IFOs specified. " \
      "Please choose up to four IFOs, but not five.", file=sys.stderr)
  sys.exit(1)

if not opts.one_ifo and not opts.two_ifo and not opts.three_ifo and \
    not opts.four_ifo and not opts.analyze_all:
  print("No number of ifos given. Please specify at least one of", file=sys.stderr)
  print("--one-ifo, --two-ifo, --three-ifo, --four-ifo", file=sys.stderr)
  print("or use --analyze-all to analyze all ifos all data", file=sys.stderr)
  sys.exit(1)
elif opts.analyze_all:
  print("The --analyze-all flag can not be used to specify the", file=sys.stderr)
  print("number of ifos to analyze. The code supports quadruple", file=sys.stderr)
  print("coincidence, so you can choose at most four instruments", file=sys.stderr)
  print("to analyze.", file=sys.stderr)
  sys.exit(1)

if opts.inspiral_veto and opts.td_follow_inspiral:
  print("Please specify only one of", file=sys.stderr)
  print("--inspiral-veto, --td-follow-inspiral.", file=sys.stderr)
 
if not opts.datafind and not opts.template_bank and not opts.write_script and \
     not opts.inspiral and not opts.sire_inspiral and \
     not opts.coincidence and not opts.coire_coincidence and \
     not opts.trigbank and not opts.inspiral_veto and \
     not opts.sire_inspiral_veto and not opts.td_follow_bank and \
     not opts.td_follow_inspiral and not opts.second_coinc and \
     not opts.coire_second_coinc and not opts.sire_second_coinc and \
     not opts.coherent_bank and not opts.coherent_inspiral and \
     not opts.summary_inspiral_triggers and not opts.summary_coinc_triggers and \
     not opts.cohire and not opts.summary_coherent_inspiral_triggers:
  print("""  No steps of the pipeline specified.
  Please specify at least one of
  --datafind, --template-bank, --inspiral, --sire-inspiral, --coincidence,
  --coire-coincidence, --trigbank, --inspiral-veto, --sire-inspiral-veto,
  --td-follow-bank, --td-follow-inspiral, --second-coinc,
  --coire-second-coinc, --sire-second-coinc, --coherent-bank,
  --coherent-inspiral, --summary-inspiral-triggers --summary-coinc-triggers,
  --cohire, --summary-coherent-inspiral-triggers""", file=sys.stderr)
  sys.exit(1)

ifo_list = ['H1','H2','L1','V1','G1']

#################################################################
# If using G1 data, rearrange ifo_list since it only uses
# the first four ifos named in ifo_list for quadruple coincidence

if opts.g1_data:
  if not opts.h1_data and ifo_list[4]=='G1':
    ifo_list[0]='G1'
    ifo_list[4]='H1'
  if not opts.h2_data and ifo_list[4]=='G1':
    ifo_list[1]='G1'
    ifo_list[4]='H2'
  if not opts.l1_data and ifo_list[4]=='G1':
    ifo_list[2]='G1'
    ifo_list[4]='L1'
  if not opts.v1_data and ifo_list[4]=='G1':
    ifo_list[3]='G1'
    ifo_list[4]='V1'

  ifo_list = ifo_list[:4]

ifotag = None
#################################
# store the values
do = {}
do['G1'] = opts.g1_data
do['H1'] = opts.h1_data
do['H2'] = opts.h2_data
do['L1'] = opts.l1_data
do['V1'] = opts.v1_data

ifo_analyze = []
for ifo in ifo_list:
  if do[ifo]: 
    ifo_analyze.append(ifo)
ifo_analyze.sort()

#################################
# analyze everything if --analyze-all set
if opts.analyze_all:
  for ifo in ifo_list:
    do[ifo] = True
  opts.one_ifo   = True
  opts.two_ifo   = True
  opts.three_ifo = True
  opts.four_ifo  = True

##############################################################################
# determine all possible coincident sets of ifos and those to analyze:
analyze = []
ifo_coincs = []

# one ifo
for ifo1 in ifo_list:
  if opts.one_ifo and do[ifo1]:
      analyze.append(ifo1)

# two ifo
for ifo1 in ifo_list:
  for ifo2 in ifo_list:
    if ifo1 < ifo2:
      ifo_coincs.append(ifo1 + ifo2)
      if opts.two_ifo and do[ifo1] and do[ifo2]:
        analyze.append(ifo1 + ifo2)

# three ifo
for ifo1 in ifo_list:
  for ifo2 in ifo_list:
    for ifo3 in ifo_list:
      if ifo1 < ifo2 and ifo2 < ifo3:
        ifo_coincs.append(ifo1 + ifo2 + ifo3)
        if opts.three_ifo and do[ifo1] and do[ifo2] and do[ifo3]:
          analyze.append(ifo1 + ifo2 + ifo3)

# four ifo
for ifo1 in ifo_list:
  for ifo2 in ifo_list:
    for ifo3 in ifo_list:
      for ifo4 in ifo_list:
        if ifo1 < ifo2 and ifo2 < ifo3 and ifo3 < ifo4:
          ifo_coincs.append(ifo1 + ifo2 + ifo3 + ifo4)
          if opts.four_ifo and do[ifo1] and do[ifo2] and do[ifo3] and do[ifo4]:
            analyze.append(ifo1 + ifo2 + ifo3 + ifo4)

ifo_combinations = copy.deepcopy(ifo_list)
ifo_combinations.extend(ifo_coincs)


##############################################################################
# try to make a directory to store the cache files and job logs
try: os.mkdir('cache')
except: pass
try: os.mkdir('logs')
except: pass

##############################################################################
# create the config parser object and read in the ini file
cp = dcConfigParser()
cp.read(opts.config_file)

##############################################################################
# if a usertag has been specified, override the config file
if opts.user_tag:
  usertag = opts.user_tag
  cp.set('pipeline','user-tag',usertag)
else:
  try:
    usertag = string.strip(cp.get('pipeline','user-tag'))
  except:
    usertag = None

##############################################################################
# if --use-gpus was specified, override the config file

if opts.use_gpus:
  cp.set('condor', 'use-gpus', '')

##############################################################################
# create a log file that the Condor jobs will write to
basename = re.sub(r'\.ini',r'',opts.config_file)
tempfile.tempdir = opts.log_path
if usertag:
  tempfile.template = basename + '.' + usertag + '.dag.log.'
else:
  tempfile.template = basename + '.dag.log.'
logfile = tempfile.mktemp()
fh = open( logfile, "w" )
fh.close()

##############################################################################
# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile, opts.dax)
if usertag:
  dag.set_dag_file(basename + '.' + usertag)
  dag.set_dax_file(basename + '.' + usertag)
else:
  dag.set_dag_file(basename)
  dag.set_dax_file(basename)

# set better submit file names than the default
if usertag:
  subsuffix = '.' + usertag + '.sub'
else:
  subsuffix = '.sub'

# set the worker package if neccesssary
try:
  dag.set_pegasus_worker(cp.get('pipeline','pegasus-worker'))
except:
  pass

##############################################################################
# create the Condor jobs that will be used in the DAG

# check number of exttrig injections
doExtTrig = cp.has_option('pipeline', 'exttrig-inj-start') and cp.has_option('pipeline', 'exttrig-inj-stop')
if doExtTrig:
  startExttrig=int(cp.get('pipeline','exttrig-inj-start'))
  stopExttrig=int(cp.get('pipeline','exttrig-inj-stop'))
  exttrigInjections=[startExttrig, stopExttrig]

  # check the values given
  if startExttrig < 1:
    print("exttrig-inj-start must be larger than 0.", file=sys.stderr)
    sys.exit(1)
  if startExttrig > stopExttrig:
    print("exttrig-inj-stop must be larger than "\
                         "exttrig-inj-start.", file=sys.stderr)
    sys.exit(1)
else:
  exttrigInjections=[0,0]


# datafind:
frame_types = []
try:
  lsync_file = cp.get('pipeline','lsync-cache-file')
  try: frame_types.append(cp.get('input','ligo-type'))
  except: pass
  try: frame_types.append(cp.get('input','virgo-type'))
  except: pass
  try: frame_types.append(cp.get('input','geo-type'))
  except: pass
  frame_types = [t for t in frame_types if t]
except:
  lsync_file = None
df_job = pipeline.LSCDataFindJob(
  'cache','logs',cp,opts.dax,lsync_file,'|'.join(frame_types))
df_job.set_sub_file( basename + '.datafind'+ subsuffix )

# tmpltbank:
tmplt_jobs = {}

for ifo in ifo_list:
  if opts.ringdown:
    tmplt_jobs[ifo] = inspiral.TmpltBankJob(cp,opts.dax)
    tmplt_jobs[ifo].set_exec_name('ringbank')
  else:
    tmplt_jobs[ifo] = inspiral.TmpltBankJob(cp,opts.dax)
  tmplt_jobs[ifo].set_sub_file( basename + '.tmpltbank_' + ifo + subsuffix )

# inspinj:
inspinj_job = inspiral.InspInjJob(cp) 
inspinj_job.set_sub_file( basename + '.inspinj' + subsuffix )

if opts.noop_inspinj:
  inspinj_job.add_condor_cmd("noop_job", "true")

# inspiral:
insp_jobs = {}

for ifo in ifo_list:
  if opts.ringdown:
    insp_jobs[ifo] = inspiral.InspiralJob(cp,opts.dax)
    insp_jobs[ifo].set_exec_name('ring')
  else:
    insp_jobs[ifo] = inspiral.InspiralJob(cp,opts.dax)
  insp_jobs[ifo].set_sub_file( basename + '.inspiral_' + ifo + subsuffix )

# create inspiral checkpoint job
insp_ckpt_job = inspiral.InspiralCkptJob(cp,opts.dax)
if cp.has_option('pipeline','remote-site'):
  insp_ckpt_job.set_executable_installed(False)

# single_inca:
single_inca_job = inspiral.IncaJob(cp)
single_inca_job.set_sub_file( basename + '.s_inca' + subsuffix )

# thinca:
thinca_jobs = {}

for ifos in ifo_coincs:
  if opts.ringdown:
    thinca_jobs[ifos] = inspiral.ThincaJob(cp,opts.dax)
    thinca_jobs[ifos].set_exec_name('rinca')
  else:
    thinca_jobs[ifos] = inspiral.ThincaJob(cp)
  thinca_jobs[ifos].set_sub_file(basename + '.thinca_' + ifos + 
      subsuffix )

# thinca_slide:
thinca_slide_jobs = {}

for ifos in ifo_coincs:
  if opts.ringdown:
    thinca_slide_jobs[ifos] = inspiral.ThincaJob(cp,opts.dax)
    thinca_slide_jobs[ifos].set_exec_name('rinca')
  else:
    thinca_slide_jobs[ifos] = inspiral.ThincaJob(cp)
  thinca_slide_jobs[ifos].set_sub_file(basename + '.thinca_slides_' + 
      ifos + subsuffix )

# td follow-up tmpltbank
tmplt_follow_jobs = {}

for ifo in ifo_list:
  tmplt_follow_jobs[ifo] = inspiral.TmpltBankJob(cp,opts.dax)
  tmplt_follow_jobs[ifo].set_sub_file( basename + '.followup_tmpltbank' + \
      subsuffix )

# trigbank:
trig_job = inspiral.TrigbankJob(cp)
trig_job.set_sub_file( basename + '.trigbank' + subsuffix )

# inspiral veto:
insp_veto_jobs = {}

for ifo in ifo_list:
  insp_veto_jobs[ifo] = inspiral.InspiralJob(cp,opts.dax)
  insp_veto_jobs[ifo].set_sub_file( basename + '.inspiral_veto_' + ifo +
    subsuffix )

# thinca:
thinca2_jobs = {}

for ifos in ifo_coincs:
  thinca2_jobs[ifos] = inspiral.ThincaJob(cp)
  thinca2_jobs[ifos].set_sub_file(basename + '.thinca2_' + ifos + 
      subsuffix )

# thinca_slide:
thinca2_slide_jobs = {}

for ifos in ifo_coincs:
  thinca2_slide_jobs[ifos] = inspiral.ThincaJob(cp)
  thinca2_slide_jobs[ifos].set_sub_file(basename + '.thinca2_slides_' + 
      ifos + subsuffix )

# coherent bank:
cb_job = inspiral.CohBankJob(cp)
cb_slide_job = inspiral.CohBankJob(cp)
if not opts.ringdown:
  cb_job.set_sub_file( basename + '.cohbank' + subsuffix )
  cb_slide_job.set_sub_file( basename + '.cohbank_slide' + subsuffix )

# trigbank coherent:
trig_coh_jobs = {}

if not opts.ringdown:
  for ifo in ifo_list:
    trig_coh_jobs[ifo] = inspiral.TrigbankJob(cp)
    trig_coh_jobs[ifo].set_sub_file( basename + '.trigbank_coherent_' + ifo +
        subsuffix )

# inspiral (coherent):
insp_coh_jobs = {}

if not opts.ringdown:
  for ifo in ifo_list:
    insp_coh_jobs[ifo] = inspiral.InspiralCoherentJob(cp,opts.dax)
    insp_coh_jobs[ifo].set_sub_file( basename + '.insp_coh_' + ifo +
      subsuffix )

# fr join:
frjoin_job = inspiral.FrJoinJob(cp)
frjoin_job.set_sub_file( basename + '.frjoin' + subsuffix )

# network trigger bank:
cohinspbank_job = inspiral.CohInspBankJob(cp)
cohinspbank_slide_job = inspiral.CohInspBankJob(cp)
if not opts.ringdown:
  cohinspbank_job.set_sub_file( basename + '.cohinspbank' + subsuffix )
  cohinspbank_slide_job.set_sub_file( basename + '.cohinspbank_slide' +
    subsuffix )

# coherent_inspiral:
cohinsp_jobs = {}
cohinsp_slide_jobs = {}

if not opts.ringdown:
  for ifos in ifo_coincs:
    cohinsp_jobs[ifos] = inspiral.ChiaJob(cp)
    cohinsp_jobs[ifos].set_sub_file(basename + '.chia_' + ifos +
        subsuffix )
    cohinsp_slide_jobs[ifos] = inspiral.ChiaJob(cp)
    cohinsp_slide_jobs[ifos].set_sub_file(basename + '.chia_slide_' + ifos +
        subsuffix )

# sire:
sire_job = inspiral.SireJob(cp)
sire_job.set_sub_file( basename + '.sire' + subsuffix )
sire_summary_job = inspiral.SireJob(cp)
sire_summary_job.set_sub_file( basename + '.sire_summary' + subsuffix )

# coire:
if opts.ringdown:
  coire_job = inspiral.CoireJob(cp)
  coire_job.set_exec_name('coincringread')
  coire_job.set_sub_file( basename + '.coincringread' + subsuffix )
  coire_slide_job = inspiral.CoireJob(cp)
  coire_slide_job.set_exec_name('coincringread')
  coire_slide_job.set_sub_file( basename + '.coincringread_slide' + subsuffix )
  coire_summary_job = inspiral.CoireJob(cp)
  coire_summary_job.set_exec_name('coincringread')
  coire_summary_job.set_sub_file( basename + '.coincringread_summary' + subsuffix )
  coire_slide_summary_job = inspiral.CoireJob(cp)
  coire_slide_summary_job.set_exec_name('coincringread')
  coire_slide_summary_job.set_sub_file( basename + '.coincringread_slide_summary' + \
    subsuffix )
else:
  coire_job = inspiral.CoireJob(cp)
  coire_job.set_sub_file( basename + '.coire' + subsuffix )
  coire_slide_job = inspiral.CoireJob(cp)
  coire_slide_job.set_sub_file( basename + '.coire_slide' + subsuffix )
  coire_summary_job = inspiral.CoireJob(cp)
  coire_summary_job.set_sub_file( basename + '.coire_summary' + subsuffix )
  coire_slide_summary_job = inspiral.CoireJob(cp)
  coire_slide_summary_job.set_sub_file( basename + '.coire_slide_summary' + \
      subsuffix )

# coire2:
coire2_job = inspiral.CoireJob(cp)
coire2_job.set_sub_file( basename + '.coire2' + subsuffix )
coire2_slide_job = inspiral.CoireJob(cp)
coire2_slide_job.set_sub_file( basename + '.coire2_slide' + subsuffix )
coire2_summary_job = inspiral.CoireJob(cp)
coire2_summary_job.set_sub_file( basename + '.coire2_summary' + subsuffix )
coire2_slide_summary_job = inspiral.CoireJob(cp)
coire2_slide_summary_job.set_sub_file( basename + '.coire2_slide_summary' + \
    subsuffix )

# cohire:
cohire_job = inspiral.CohireJob(cp)
cohire_slide_job = inspiral.CohireJob(cp)
cohire_job.set_sub_file( basename + '.cohire' + subsuffix )
cohire_slide_job.set_sub_file( basename + '.cohire_slide' + subsuffix )
cohire_summary_job = inspiral.CohireJob(cp)
cohire_slide_summary_job = inspiral.CohireJob(cp)
cohire_summary_job.set_sub_file( basename + '.cohire_summary' + subsuffix )
cohire_slide_summary_job.set_sub_file( basename + '.cohire_slide_summary' +
    subsuffix )

all_jobs = [inspinj_job, single_inca_job, trig_job, cb_job, sire_job, 
  sire_summary_job, coire_job, coire_slide_job, coire_summary_job,
  coire_slide_summary_job, coire2_job, coire2_slide_job, coire2_summary_job, 
  coire2_slide_summary_job, cohinspbank_job, cohire_job, cohire_summary_job ]
all_jobs.extend(tmplt_jobs.values())
all_jobs.extend(insp_jobs.values())
all_jobs.extend(thinca_jobs.values())
all_jobs.extend(thinca_slide_jobs.values())
all_jobs.extend(tmplt_follow_jobs.values())
all_jobs.extend(insp_veto_jobs.values())
all_jobs.extend(thinca2_jobs.values())
all_jobs.extend(thinca2_slide_jobs.values())
all_jobs.extend(trig_coh_jobs.values())
all_jobs.extend(insp_coh_jobs.values())
all_jobs.extend(cohinsp_jobs.values())

##############################################################################
# set the usertag in the jobs
if usertag:
  for job in all_jobs:
    if not job.get_opts().has_key('user-tag'):
      job.add_opt('user-tag',usertag)

all_jobs.append(df_job)
all_jobs.append(insp_ckpt_job)
#all_jobs.append(frjoin_job)

##############################################################################
# set the condor job priority
if opts.priority:
  for job in all_jobs:
    job.add_condor_cmd('priority',str(opts.priority))


##############################################################################
# read in the number of time-slides from the ini file 
try:
  num_slides = int(cp.get('input','num-slides'))
except:
  num_slides = None

if num_slides:
  coire_slide_job.add_opt('num-slides', str(num_slides))
  coire_slide_summary_job.add_opt('num-slides', str(num_slides))

##############################################################################
# read in the maximum length of a science segment for thinca from the ini file
try:
  max_thinca_segment = int(cp.get('input','max-thinca-segment'))
except:
  max_thinca_segment = None
     
##############################################################################
# read in the GPS start and end times from the ini file
# Only used to set the cache files and coire files.
# XXX Should we use this to cut segments??
try:
  gps_start_time = cp.getint('input','gps-start-time')
except:
  gps_start_time = None

try:
  gps_end_time = cp.getint('input','gps-end-time')
except:
  gps_end_time = None

#############################################################################
# read in playground data mask from ini file 
# set the playground_only option and add to inca and sire jobs
try:
  play_data_mask = string.strip(cp.get('pipeline','playground-data-mask'))
except:
  play_data_mask = None

play_jobs = [trig_job, sire_job, sire_summary_job, 
  coire_job, coire_slide_job, coire_summary_job, coire_slide_summary_job,
  coire2_job, coire2_slide_job, coire2_summary_job, coire2_slide_summary_job,
  cohire_job, cohire_summary_job]
play_jobs.extend(thinca_jobs.values())
play_jobs.extend(thinca_slide_jobs.values())
play_jobs.extend(thinca2_jobs.values())
play_jobs.extend(thinca2_slide_jobs.values())
play_jobs.extend(trig_coh_jobs.values())

if play_data_mask == 'playground_only':
  playground_only = 2
  
  single_inca_job.add_opt('playground-only','')
  for job in play_jobs:
    job.add_opt('data-type','playground_only')

elif play_data_mask == 'exclude_playground':
  playground_only = 0
  
  single_inca_job.add_opt('no-playground','')
  for job in play_jobs:
    job.add_opt('data-type','exclude_play')


elif play_data_mask == 'all_data':
  playground_only = 0

  single_inca_job.add_opt('all-data','')
  for job in play_jobs:
    job.add_opt('data-type','all_data')

else:
  print("Invalid playground data mask " + play_data_mask + " specified")
  sys.exit(1)

 
 
##############################################################################
# get the pad and chunk lengths from the values in the ini file
pad = int(cp.get('data', 'pad-data'))
if opts.ringdown:
  length = int(cp.get('data','block-duration'))
  overlap = int(cp.get('data','segment-duration'))/2
else:
  n = int(cp.get('data', 'segment-length'))
  s = int(cp.get('data', 'number-of-segments'))
  r = int(cp.get('data', 'sample-rate'))
  o = int(cp.get('inspiral', 'segment-overlap'))
  length = ( n * s - ( s - 1 ) * o ) / r
  overlap = o / r


##############################################################################
#  The meat of the DAG generation comes below
#
#
#  The various data sets we compute are:
# 
#  data[ifo] : the science segments and master chunks
#
#  data_out[ifo] : the analyzable data 
#
#  not_data_out[ifo] : non analyzable data
#
#  analyzed_data[ifos] : the 1,2,3,4 ifo coincident data 
#
#  data_to_do[ifo] : the data to analyze for each ifo
#       (depends upon which of single,double,triple, quadruple data we analyze) 
#
#  And the lists of jobs are:
#
#  chunks_analyzed[ifo] : list of chunks analyzed for each ifo
#
#  single_coinc_nodes[ifo] : the single coincident inca jobs
#
#  coinc_nodes[ifos] : the double, triple, quadruple coincident thinca nodes

#  coinc_slide_nodes[ifos] : the double, triple, quadruple coincident 
#                            thinca slide nodes
#
#  analyzed_coinc_data[ifos][ifo] : the single ifo analyzed data (post coinc)
#                                   for each coincident combination, we have
#                                   data from each of the active ifos
#
#  coinc2_nodes[ifos] : the double, triple, quadruple coincident thinca nodes
#                       from the second coincidence step
#
#
#  coinc2_slide_nodes[ifos] : the double, triple, quadruple coincident 
#                             thinca slide nodes from the second coincidence 
#                             step
##############################################################################



##############################################################################
#   Step 1: read science segs that are greater or equal to a chunk 
#   from the input file

print("reading in single ifo science segments and creating master chunks...", end=' ')
sys.stdout.flush()

segments = {}
data = {}

for ifo in ifo_list:
  try:
    segments[ifo] = cp.get('input', ifo +'-segments')
  except:
    segments[ifo] = None
  
  data[ifo] = pipeline.ScienceData() 
  if segments[ifo]:
    data[ifo].read(segments[ifo],length + 2 * pad) 
    data[ifo].make_chunks(length,overlap,playground_only,0,overlap/2,pad)
    data[ifo].make_chunks_from_unused(length,overlap/2,playground_only,
        0,0,overlap/2,pad)

print("done")

# work out the earliest and latest times that are being analyzed
if not gps_start_time:
  gps_start_time = 10000000000
  for ifo in ifo_list:
    if data[ifo] and (data[ifo][0].start() < gps_start_time):
      gps_start_time = data[ifo][0].start()
  print("GPS start time not specified, obtained from segment lists as " + \
    str(gps_start_time))


if not gps_end_time:
  gps_end_time = 0
  for ifo in ifo_list:
    if data[ifo] and (data[ifo][-1].end() > gps_end_time):
      gps_end_time = data[ifo][0].end()
  print("GPS end time not specified, obtained from segment lists as " + \
    str(gps_end_time))

##############################################################################
#   Step 2: determine analyzable times

data_out = {}
not_data_out = {}
for ifo in ifo_list:
  data_out[ifo] = copy.deepcopy(data[ifo])

  # remove start and end of science segments which aren't analyzed for triggers
  for seg in data_out[ifo]:
    seg.set_start(seg.start() + overlap/2 + pad)
    seg.set_end(seg.end() - overlap/2 - pad)

  if playground_only:
    data_out[ifo].play()

  not_data_out[ifo] = copy.deepcopy(data_out[ifo])
  not_data_out[ifo].coalesce()
  not_data_out[ifo].invert()

# determine the data we can analyze for various detector combinations
analyzed_data = {}

# determine the coincident data, if it is to be analyzed
for ifos in ifo_combinations:
  analyzed_data[ifos] = pipeline.ScienceData()
  if ifos in analyze:
    selected_data = []
    for ifo in ifo_list:
      if ifo in ifos:
        selected_data.append(data_out[ifo])
      else:
        selected_data.append(not_data_out[ifo])
    analyzed_data[ifos] = generate_segments(selected_data[0], 
        selected_data[1], selected_data[2], selected_data[3])

    if max_thinca_segment:
      analyzed_data[ifos].split(max_thinca_segment)
  
##############################################################################
# Step 3: Compute the Science Segments to analyze

data_to_do = {}

for ifo in ifo_list:
  data_to_do[ifo] = copy.deepcopy(analyzed_data[ifo])
  for ifos in analyze:
    if ifo in ifos:
      data_to_do[ifo].union(analyzed_data[ifos])
  data_to_do[ifo].coalesce() 

##############################################################################
# Step 3b: Set up the injection job

# if doing injections then set up the injection analysis
try: seed = cp.get("input","injection-seed")
except: seed = None

if cp.has_option("input","hardware-injection"):
  inj_file_loc = cp.get("input","hardware-inj-file")
  inspinj = inspiral.InspInjNode(inspinj_job)
  inspinj.set_start(gps_start_time)
  inspinj.set_end(gps_end_time)
  inspinj.set_seed(0)
  inj_file = inspinj.get_output()
  print(inj_file)
  shutil.copy( inj_file_loc, inj_file)
  inspinj = None
elif seed:
  inspinj = inspiral.InspInjNode(inspinj_job)
  inspinj.set_start(gps_start_time)
  inspinj.set_end(gps_end_time)
  inspinj.set_seed(seed)
  if opts.inspinj: dag.add_node(inspinj)
  inj_file = inspinj.get_output()
else: 
  inj_file = None
  inspinj = None


# add the injection details to inspiral/sire/coire jobs  
if inj_file:
  inj_jobs = [sire_job, sire_summary_job, 
    coire_job, coire_slide_job, coire_summary_job, coire_slide_summary_job, 
    coire2_job, coire2_slide_job, coire2_summary_job, coire2_slide_summary_job,
    cohire_job, cohire_summary_job]
  inj_jobs.extend(insp_jobs.values())
  inj_jobs.extend(insp_veto_jobs.values())
  inj_jobs.extend(insp_coh_jobs.values())

  for job in inj_jobs:  
     job.add_file_opt('injection-file',inj_file)

if inj_file or doExtTrig:
  # add the injection coincidence to the sire jobs
  try: 
    sire_job.add_ini_opts(cp,'sire-inj')
    sire_summary_job.add_ini_opts(cp,'sire-inj')
  except: pass

  try: 
    coire_job.add_ini_opts(cp,'coire-inj')
    coire_slide_job.add_ini_opts(cp,'coire-inj')
    coire_summary_job.add_ini_opts(cp,'coire-inj')
    coire_slide_summary_job.add_ini_opts(cp,'coire-inj')
    coire2_job.add_ini_opts(cp,'coire-inj')
    coire2_slide_job.add_ini_opts(cp,'coire-inj')
    coire2_summary_job.add_ini_opts(cp,'coire-inj')
    coire2_slide_summary_job.add_ini_opts(cp,'coire-inj')  
    cohire_job.add_ini_opts(cp,'cohire-inj')
    cohire_summary_job.add_ini_opts(cp,'cohire-inj')
  except: pass

##############################################################################
# Step 4: Determine which of the master chunks needs to be filtered
chunks_analyzed = {}

prev_df = None

for ifo in ifo_list:
  print("setting up jobs to filter " + ifo + " data...", end=' ')
  sys.stdout.flush()

  (prev_df,chunks_analyzed[ifo]) = analyze_ifo(ifo,data[ifo],data_to_do[ifo],  
      tmplt_jobs[ifo],insp_jobs[ifo],df_job,prev_df,dag,exttrigInjections,
      usertag,inspinj,insp_ckpt_job)

  # Step 4S: Run sire on the single ifo triggers
  sire_nodes = sire_segments_individually(chunks_analyzed[ifo], sire_job, 
      dag, opts.sire_inspiral, opts.inspiral, ifo, inj_file = inj_file, 
      ifotag="FIRST", usertag = usertag, inspinjNode = inspinj)

  if len(sire_nodes):
    sire_segments(sire_nodes, sire_summary_job, dag, 
        opts.summary_inspiral_triggers, opts.sire_inspiral, ifo, 
        gps_start_time, gps_end_time, inj_file = inj_file, 
        ifotag="SUMMARY_FIRST", usertag = usertag, 
        inspinjNode=inspinj)

  print("done") 


##############################################################################
# Step 5: Run inca in single ifo mode on the single ifo triggers.

if not opts.ringdown:
  print("setting up jobs to inca single IFO data...", end=' ')
  sys.stdout.flush()

  single_coinc_nodes = {}
  for ifo in ifo_list:
    single_coinc_nodes[ifo] = []
    single_coinc_nodes[ifo] = single_coinc(chunks_analyzed[ifo],ifo,
       analyzed_data[ifo],single_inca_job,dag,usertag)
   
    # Step 5S: Run sire on the single ifo triggers
    sire_segments_individually(single_coinc_nodes[ifo], sire_job, dag, 
        opts.coire_coincidence, opts.coincidence, ifo,
        inj_file=inj_file, usertag=usertag, ifotag="FIRST", inspinjNode = inspinj)

  print("done")
  
 
##############################################################################
# Step 6: Run thinca on each of the disjoint sets of coincident data

coinc_nodes = {}
coinc_slide_nodes = {}

for ifos in ifo_coincs:
  print("setting up thinca jobs on " + ifos + " data...", end=' ')
  sys.stdout.flush()

  if cp.has_section('thinca-1'): thinca_jobs[ifos].add_ini_opts(cp, 'thinca-1')
  if cp.has_section('coire-1'):
    coire_job.add_ini_opts(cp, 'coire-1')
    coire_slide_job.add_ini_opts(cp, 'coire-1')
    coire_summary_job.add_ini_opts(cp, 'coire-1')
    coire_slide_summary_job.add_ini_opts(cp, 'coire-1')

  coinc_nodes[ifos] = thinca_coinc(ifos,chunks_analyzed,analyzed_data[ifos],
      thinca_jobs[ifos],dag,opts.coincidence,opts.inspiral,exttrigInjections,
      usertag, coire_job=coire_job, do_coire=opts.coire_coincidence)

  # Step 6 coire: Run coire on coinc triggers
  if not doExtTrig:
    coire_segments_individually(coinc_nodes[ifos], coire_job, dag, 
        opts.coire_coincidence, opts.coincidence, ifos, inj_file=inj_file, 
        ifotag="FIRST", usertag=usertag, inspinjNode=inspinj)

    if opts.summary_first_coinc_triggers and len(coinc_nodes[ifos]):
      coire_segments(coinc_nodes[ifos], coire_summary_job, dag,
          opts.summary_first_coinc_triggers,
          opts.coire_coincidence, ifos, gps_start_time, gps_end_time,
          inj_file=inj_file, ifotag="SUMMARY_FIRST_" + ifos, usertag=usertag, 
          inspinjNode=inspinj)

  # Step 6 slide: Time slides 
  coinc_slide_nodes[ifos] = []
  if num_slides:
    if cp.has_section('thinca-1'): 
      thinca_slide_jobs[ifos].add_ini_opts(cp, 'thinca-1')
    if cp.has_section('coire-1'):
      coire_job.add_ini_opts(cp, 'coire-1')
      coire_slide_job.add_ini_opts(cp, 'coire-1')
      coire_summary_job.add_ini_opts(cp, 'coire-1')
      coire_slide_summary_job.add_ini_opts(cp, 'coire-1')

    coinc_slide_nodes[ifos] = thinca_coinc(ifos,chunks_analyzed,
        analyzed_data[ifos],thinca_slide_jobs[ifos],dag,opts.coincidence,
        opts.inspiral,exttrigInjections, usertag, coire_job=coire_slide_job, 
        do_coire=opts.coire_coincidence, ifotag=ifotag, slides=num_slides )

    if not doExtTrig:
      coire_segments_individually(coinc_slide_nodes[ifos], coire_slide_job, 
          dag, opts.coire_coincidence, opts.coincidence, ifos, 
          inj_file = inj_file, num_slides = num_slides, ifotag="FIRST", 
          usertag=usertag,inspinjNode=inspinj)

      if opts.summary_first_coinc_triggers and len(coinc_nodes[ifos]):
        coire_segments(coinc_slide_nodes[ifos], coire_slide_summary_job, dag,
            opts.summary_first_coinc_triggers,
            opts.coire_coincidence, ifos, gps_start_time, gps_end_time,
            inj_file=inj_file, ifotag="SLIDE_SUMMARY_FIRST_" + ifos, usertag=usertag,
            inspinjNode=inspinj)

  # Concatenate the zerolag and slide nodes
  coinc_nodes[ifos] = coinc_nodes[ifos] + coinc_slide_nodes[ifos]
  
  print("done")


##############################################################################
# Step 7: Run trigbank and inspiral on each instrument in two ifo times
#         or run the TD follow-up

# For following up several approximants, we get the list here...
if opts.td_follow_bank or opts.td_follow_inspiral:
  try:
    approximants = cp.get('veto-inspiral', 'approximant')
    cp.remove_option('veto-inspiral', 'approximant')
  except configparser.NoOptionError:
    print("No approximant in veto-inspiral - using the main one...")
    approximants = cp.get('inspiral', 'approximant')

  approximants = approximants.split(',')
else:
  approximants = None

# Initialise the banks
banks = {}

for ifos in ifo_coincs:
  banks[ifos] = {}

  for ifo in ifo_list:
    if ifo in ifos:
      banks[ifos][ifo] = None

analyzed_coinc_data = {}

for ifos in ifo_coincs:
  analyzed_coinc_data[ifos] = {}
 
  for ifo in ifo_list:
    if ifo in ifos:
      print("setting up jobs to filter " + ifo + \
        " data with coinc trigs from " + ifos + " times ...", end=' ')
      sys.stdout.flush()

      if opts.td_follow_bank or opts.td_follow_inspiral:
        for this_approx in approximants:
          this_approx = this_approx.strip()
          analyzed_coinc_data[ifos][ifo], banks[ifos][ifo] = \
            follow_up(data[ifo],ifo,coinc_nodes[ifos],ifos,
            tmplt_follow_jobs[ifo],insp_veto_jobs[ifo],dag,this_approx,usertag,
            banks[ifos][ifo])
      else:
        analyzed_coinc_data[ifos][ifo] = trig_inspiral(data[ifo],ifo,
            coinc_nodes[ifos],ifos,trig_job,insp_veto_jobs[ifo],dag,
            exttrigInjections,usertag,inspinjNode=inspinj,
            insp_ckpt_job=insp_ckpt_job)
      
        sire_segments_individually(analyzed_coinc_data[ifos][ifo], sire_job, 
              dag, opts.sire_inspiral_veto, opts.inspiral_veto, ifo, 
              inj_file = inj_file, usertag = usertag, 
              ifotag="SNGL_SECOND_" + ifo, inspinjNode = inspinj)

      print("done")


##############################################################################
# Step 8: Run thinca2 on each of the disjoint sets of coincident data

coinc2_nodes = {}
coinc2_slide_nodes = {}
coire2_nodes = {}
coire2_slide_nodes = {}

for ifos in ifo_coincs:
  print("setting up second thinca jobs on " + ifos + " data...", end=' ')
  sys.stdout.flush()

  if cp.has_section('thinca-2'):
    thinca2_jobs[ifos].add_ini_opts(cp, 'thinca-2')
  if cp.has_section('coire-2'):
    coire2_job.add_ini_opts(cp, 'coire-2')
    coire2_slide_job.add_ini_opts(cp, 'coire-2')
    coire2_summary_job.add_ini_opts(cp, 'coire-2')
    coire2_slide_summary_job.add_ini_opts(cp, 'coire-2')
  if opts.td_follow_inspiral:
    for this_approx in approximants:
      this_approx = this_approx.strip()
      coinc2_nodes[ifos] = thinca_coinc(ifos, analyzed_coinc_data[ifos],
          analyzed_data[ifos], thinca2_jobs[ifos], dag, opts.second_coinc,
          opts.td_follow_inspiral, exttrigInjections, usertag=usertag,
          coire_job=coire2_job, do_coire=opts.coire_second_coinc, 
          ifotag=ifos + "_" + this_approx)

  else:
    coinc2_nodes[ifos] = thinca_coinc(ifos, analyzed_coinc_data[ifos],
      analyzed_data[ifos], thinca2_jobs[ifos], dag, opts.second_coinc,
      opts.inspiral_veto, exttrigInjections, usertag=usertag,
      coire_job=coire2_job, do_coire=opts.coire_second_coinc, ifotag=ifos)

    if not doExtTrig:
      coire2_nodes[ifos] = coire_segments_individually(coinc2_nodes[ifos], 
          coire2_job, dag, opts.coire_second_coinc, opts.second_coinc, 
          ifos, inj_file=inj_file, ifotag="SECOND_" + ifos, usertag=usertag,
          inspinjNode=inspinj)
 
      if len(coire2_nodes[ifos]):
        coire_segments(coire2_nodes[ifos], coire2_summary_job, dag,
            opts.summary_coinc_triggers, 
            opts.coire_second_coinc, ifos, gps_start_time, gps_end_time, 
            inj_file=inj_file, ifotag="SUMMARY_SECOND_" + ifos, usertag=usertag,
            inspinjNode=inspinj)

      for ifo in ifo_list:
        if ifo in ifos: 

          sire_nodes = sire_segments_individually(coire2_nodes[ifos], sire_job, 
              dag, opts.sire_second_coinc, opts.coire_second_coinc, ifo, 
              inj_file=inj_file, ifotag="SECOND_" + ifos, usertag=usertag,
              inspinjNode = inspinj)

          if len(sire_nodes):
            sire_segments(sire_nodes, sire_summary_job, dag, 
                opts.summary_single_ifo_triggers, opts.sire_second_coinc, ifo, 
                gps_start_time, gps_end_time, inj_file = inj_file, 
                ifotag="SUMMARY_SECOND_" + ifos, usertag = usertag,
                inspinjNode=inspinj)


  # Step 8slide: Time slides on the data
  if num_slides:
    if cp.has_section('thinca-2'): 
      thinca2_slide_jobs[ifos].add_ini_opts(cp, 'thinca-2')
  
    if opts.td_follow_inspiral:
      for this_approx in approximants:
        this_approx = this_approx.strip()
        coinc2_slide_nodes[ifos] = thinca_coinc(ifos,
          analyzed_coinc_data[ifos], analyzed_data[ifos],
          thinca2_slide_jobs[ifos], dag, opts.second_coinc,
          opts.inspiral_veto, exttrigInjections, usertag=usertag,
          coire_job=coire2_slide_job, do_coire=opts.coire_second_coinc, 
          ifotag=ifos+'_'+this_approx, slides=num_slides)
    else: 
      coinc2_slide_nodes[ifos] = thinca_coinc(ifos,
        analyzed_coinc_data[ifos], analyzed_data[ifos],
        thinca2_slide_jobs[ifos], dag, opts.second_coinc, opts.inspiral_veto,
        exttrigInjections, usertag=usertag, coire_job=coire2_slide_job,
        do_coire=opts.coire_second_coinc, ifotag=ifos, slides=num_slides)

      if not doExtTrig:
        coire2_slide_nodes[ifos] = coire_segments_individually(
            coinc2_slide_nodes[ifos], coire2_slide_job, 
            dag, opts.coire_second_coinc, opts.second_coinc, ifos, 
            num_slides=num_slides, ifotag="SECOND_" + ifos,
            usertag=usertag, inspinjNode=inspinj)

        if len(coire2_slide_nodes[ifos]):
          coire_segments(coire2_slide_nodes[ifos], coire2_slide_summary_job, dag, 
              opts.summary_coinc_triggers, opts.coire_second_coinc, ifos, 
              gps_start_time, gps_end_time, num_slides=num_slides, 
              ifotag="SUMMARY_SECOND_" + ifos, usertag=usertag, 
              inspinjNode=inspinj)

        for ifo in ifo_list:
          if ifo in ifos:
            sire_nodes = sire_segments_individually(coire2_slide_nodes[ifos], 
                sire_job, dag, opts.sire_second_coinc, opts.coire_second_coinc,
                ifo, inj_file = inj_file, ifotag="SLIDE_SECOND_" + ifos, 
                usertag=usertag, inspinjNode = inspinj)

            if len(sire_nodes):
              sire_segments(sire_nodes, sire_summary_job, dag, 
                  opts.summary_single_ifo_triggers, opts.sire_second_coinc, 
                  ifo, gps_start_time, gps_end_time, inj_file = inj_file, 
                  ifotag="SUMMARY_SECOND_" + ifos, usertag = usertag)


  print("done")

##############################################################################
# Step 9: Create coherent (template) banks:

if not opts.ringdown:
  cohbank_nodes = {}
  if num_slides: cohbank_slide_nodes = {}

  for ifos in ifo_coincs:
    print("setting up coherent-bank jobs on " + ifos + " data...", end=' ')
    sys.stdout.flush()

    if cp.has_section('thinca-2') :
      cb_job.add_ini_opts(cp,'thinca-2')

    if num_slides:
      #CHECK:cb_job[ifos].add_opt('num-slides',cp.get('input','num-slides'))
      cohbank_nodes[ifos] = coherent_bank(ifos,cb_job,dag,coinc2_nodes[ifos],
          exttrigInjections,0,0,usertag=usertag)
      if cp.has_section('thinca-slide') :
        cb_slide_job.add_ini_opts(cp,'thinca-slide')
      cb_slide_job.add_opt('num-slides', str(num_slides))
      if cp.has_section('thinca-2') :
        cb_slide_job.add_ini_opts(cp,'thinca-2')
      #cb_job[ifos].add_var_opt('num-slides',num_slides)
      cohbank_slide_nodes[ifos] = coherent_bank(ifos,cb_slide_job,dag,
          coinc2_slide_nodes[ifos],exttrigInjections,num_slides,1,
          usertag=usertag)
    else:
      cohbank_nodes[ifos] = coherent_bank(ifos,cb_job,dag,coinc2_nodes[ifos],
          exttrigInjections,0,0,usertag=usertag)

    print("done")

##############################################################################
# Step 10: Construct single-ifo CData frame files for the coherent analyses:

if not opts.ringdown:
  analyzed_trig_coh_data = {}
  data_opts = {}
  if num_slides: analyzed_trig_coh_slide_data = {}

  for ifos in ifo_coincs:
    analyzed_trig_coh_data[ifos] = {}
    if num_slides: analyzed_trig_coh_slide_data[ifos] = {}

    for ifo in ifo_list:
      if ifo in ifos:
        data_opts[ifo] = {}

        print("setting up jobs to filter " + ifo + \
          " data with cohbank templates from " + ifos + " times ...", end=' ')
        sys.stdout.flush()

        # add ifo specific options
        if ifo == 'G1':
          data_opts[ifo] = 'geo-data'
          channel = cp.get('input','geo-channel')
        elif ifo == 'V1':
          data_opts[ifo] = 'virgo-data'
          channel = cp.get('input','virgo-channel')
        else:
          data_opts[ifo] = 'ligo-data'
          channel = cp.get('input','ligo-channel')

        # set up the trigbanks
        if cp.has_section('trigbank-coherent'):
          trig_coh_jobs[ifo].add_ini_opts(cp,'trigbank-coherent')
        # set up inspirals:
        if cp.has_option("input","output-path"):
          insp_coh_jobs[ifo].add_opt('output-path',cp.get("input","output-path"))
        if cp.has_option("input","username"):
          insp_coh_jobs[ifo].add_opt('username',cp.get("input","username"))
        if cp.has_section('inspiral-coherent'):
          insp_coh_jobs[ifo].add_ini_opts(cp,'inspiral-coherent')
        if cp.has_section('veto-inspiral'):
          insp_coh_jobs[ifo].add_opt('chisq-bins', \
                               cp.get('veto-inspiral', 'chisq-bins'))
        if cp.has_section(ifo.lower() + '-inspiral-coherent'):
          insp_coh_jobs[ifo].add_ini_opts(cp,ifo.lower() + '-inspiral-coherent')
        if cp.has_section(data_opts[ifo]):
          insp_coh_jobs[ifo].add_ini_opts(cp,data_opts[ifo])
        insp_coh_jobs[ifo].set_channel(channel)
        # see if we are using calibrated data
        if cp.has_section(data_opts[ifo]) and \
            cp.has_option(data_opts[ifo],'calibrated-data'):
          calibrated = True
        else: calibrated = False

        analyzed_trig_coh_data[ifos][ifo] = coherent_trig_inspiral(data[ifo],
          ifo,cohbank_nodes[ifos],ifos,trig_coh_jobs[ifo],insp_coh_jobs[ifo],
          dag,calibrated,exttrigInjections,0,usertag,inspinj)

	if num_slides:
	  analyzed_trig_coh_slide_data[ifos][ifo] = coherent_trig_inspiral(
	    data[ifo],ifo,cohbank_slide_nodes[ifos],ifos,trig_coh_jobs[ifo],
	    insp_coh_jobs[ifo],dag,calibrated,exttrigInjections,1,usertag,inspinj)

        print("done")

##############################################################################
# Step 11: Create network trigger banks with inspiral-coherent triggers:

if not opts.ringdown:
  cohinspbank_nodes = {}

  if num_slides:
    cohinspbank_slide_nodes = {}

  if cp.has_section('cohinspbank'):
    cohinspbank_job.add_ini_opts(cp, 'cohinspbank')

  for ifos in ifo_coincs:
    print("setting up network-trigger-bank jobs on " + ifos + " data...", end=' ')
    sys.stdout.flush()

    cohinspbank_nodes[ifos] = coherent_inspiral_bank(ifos,cohinspbank_job,dag,
        analyzed_data[ifos],analyzed_trig_coh_data[ifos],exttrigInjections,
        0,0,usertag=usertag)

    if num_slides:
      if cp.has_section('cohinspbank'):
        cohinspbank_slide_job.add_ini_opts(cp, 'cohinspbank')
      if cp.has_section('thinca-slide') :
        cohinspbank_slide_job.add_ini_opts(cp,'thinca-slide')

      cohinspbank_slide_job.add_opt('num-slides', str(num_slides))

      cohinspbank_slide_nodes[ifos] = coherent_inspiral_bank(ifos,
        cohinspbank_slide_job,dag,analyzed_data[ifos],
        analyzed_trig_coh_slide_data[ifos],exttrigInjections,
        num_slides,1,usertag=usertag)

    print("done")

##############################################################################
# Step 12: Do the coherent analyses:

if not opts.ringdown:
  cohinsp_nodes = {}
  cohire_nodes = {}

  if num_slides:
    cohinsp_slide_nodes = {}
    cohire_slide_nodes = {}

  if cp.has_section('cohire'):
    cohire_job.add_ini_opts(cp, 'cohire')
    cohire_summary_job.add_ini_opts(cp, 'cohire')

    if num_slides:
      cohire_slide_job.add_ini_opts(cp, 'cohire')
      cohire_slide_summary_job.add_ini_opts(cp, 'cohire')

  for ifos in ifo_coincs:
    print("setting up coherent-analysis jobs on " + ifos + " data...", end=' ')
    sys.stdout.flush()

    if cp.has_section('data') and not opts.ringdown:
      cohinsp_jobs[ifos].add_opt('segment-length',cp.get('data','segment-length'))
      cohinsp_jobs[ifos].add_opt('sample-rate',cp.get('data','sample-rate'))
    # FIXME: this exponent can be different for virgo;
    # discontinue use of this variable
    if cp.has_section('ligo-data') and not opts.ringdown:
      cohinsp_jobs[ifos].add_opt('dynamic-range-exponent', \
                                    cp.get('ligo-data','dynamic-range-exponent'))

    cohinsp_nodes[ifos] = coherent_analysis(ifos,data,cohinsp_jobs[ifos],dag,
        cohinspbank_nodes[ifos],analyzed_trig_coh_data[ifos],
        exttrigInjections,0,usertag=usertag,inspinjNode=inspinj)

    if num_slides:
      cohinsp_slide_nodes[ifos] = coherent_analysis(ifos,data,cohinsp_jobs[ifos],
          dag,cohinspbank_slide_nodes[ifos],analyzed_trig_coh_data[ifos],
          exttrigInjections,1,usertag=usertag,inspinjNode=inspinj)

    if not doExtTrig:
      cohire_nodes[ifos] = cohire_segments_individually(cohinsp_nodes[ifos],
          cohire_job, dag, opts.cohire, opts.coherent_inspiral,
          ifos, 0, inj_file=inj_file, ifotag="COHERENT_" + ifos,
          usertag=usertag,inspinjNode=inspinj)

      if num_slides:
        cohire_slide_nodes[ifos] = cohire_segments_individually(cohinsp_slide_nodes[ifos],
          cohire_slide_job, dag, opts.cohire, opts.coherent_inspiral,
          ifos, num_slides, inj_file=inj_file, ifotag="COHERENT_" + ifos,
          usertag=usertag,inspinjNode=inspinj)

    if len(cohire_nodes[ifos]):
      cohire_segments(cohire_nodes[ifos], cohire_summary_job, dag,
             opts.summary_coherent_inspiral_triggers,
             opts.cohire, ifos, gps_start_time, gps_end_time,
             0, inj_file=inj_file, ifotag="SUMMARY_COHERENT_" + ifos,
             usertag=usertag,inspinjNode=inspinj)

      if num_slides:
        cohire_segments(cohire_slide_nodes[ifos], cohire_slide_summary_job, dag,
             opts.summary_coherent_inspiral_triggers,
             opts.cohire, ifos, gps_start_time, gps_end_time,
             num_slides, inj_file=inj_file, ifotag="SUMMARY_COHERENT_" + ifos,
             usertag=usertag,inspinjNode=inspinj)

    print("done")

##############################################################################
# Step 13: Write out the LAL cache files for the various output data

if gps_start_time is not None and gps_end_time is not None:
  print("generating cache files for output data products...", end=' ')
  cache_fname = ''
  for ifo in ifo_analyze: 
    cache_fname += ifo
  cache_fname += '-INSPIRAL_HIPE' 
  if usertag: cache_fname += '_' + usertag
  cache_fname += '-' + str(gps_start_time) + '-' + \
    str(gps_end_time - gps_start_time) + '.cache'
  output_data_cache = lal.Cache()

  for node in dag.get_nodes():
    if opts.dax and isinstance(node,pipeline.LSCDataFindNode):
      # ignore datafind nodes, as their output is a cache file
      continue
    # add the data generated by the job to the output data cache
    output_file = node.get_output()

    output_data_cache.append(lal.Cache.from_urls([output_file])[0])
    if (isinstance(node,inspiral.CoireNode) or \
        isinstance(node,inspiral.SireNode)) and \
        node.get_missed():
      output_data_cache.append(lal.Cache.from_urls([node.get_missed()])[0])

  output_data_cache.tofile(open(cache_fname, "w"))
  print("done")
else:
  print("gps start and stop times not specified: cache files not generated")


##############################################################################
# Step 14: Setup the Maximum number of jobs for different categories
if not opts.disable_dag_categories:
  for cp_opt in cp.options('condor-max-jobs'):
    dag.add_maxjobs_category(cp_opt,cp.getint('condor-max-jobs',cp_opt))

# Add number of retries to jobs as specified by ihope.ini
if cp.has_option("pipeline", "retry-jobs"):
  num_retries = cp.getint("pipeline", "retry-jobs")
  for node in dag.get_nodes():
    node.set_retry(num_retries)

##############################################################################
# Step 15: Write out the DAG, help message and log file
dag.write_sub_files()
dag.write_dag()

if opts.write_script:
  dag.write_script()

##############################################################################  
# write a message telling the user that the DAG has been written
if opts.dax:
  
  print("\nCreated an abstract DAX file", dag.get_dag_file())
  print("which can be transformed into a concrete DAG with gencdag.")
  print("\nSee the documentation on http://www.lsc-group.phys.uwm.edu/lscdatagrid/griphynligo/pegasus_lsc.html")



else:
  print("\nCreated a DAG file which can be submitted by executing")
  print("\n   condor_submit_dag", dag.get_dag_file())
  print("""\nfrom a condor submit machine (e.g. hydra.phys.uwm.edu)\n
  If you are running LSCdataFind jobs, do not forget to initialize your grid 
  proxy certificate on the condor submit machine by running the commands
  
    unset X509_USER_PROXY
    grid-proxy-init -hours 72

  Enter your pass phrase when prompted. The proxy will be valid for 72 hours. 
  If you expect the LSCdataFind jobs to take longer to complete, increase the
  time specified in the -hours option to grid-proxy-init. You can check that 
  the grid proxy has been sucessfully created by executing the command:
  
    grid-cert-info -all -file /tmp/x509up_u`id -u`
  
  This will also give the expiry time of the proxy. You should also make sure
  that the environment variable LSC_DATAFIND_SERVER is set the hostname and
  optional port of server to query. For example on the UWM medusa cluster this
  you should use
  
    export LSC_DATAFIND_SERVER=dataserver.phys.uwm.edu
  
  Contact the administrator of your cluster to find the hostname and port of the
  LSCdataFind server.
  """)

##############################################################################
# write out a log file for this script
if usertag:
  log_fh = open(basename + '.pipeline.' + usertag + '.log', 'w')
else:
  log_fh = open(basename + '.pipeline.log', 'w')
  
# FIXME: the following code uses obsolete CVS ID tags.
# It should be modified to use git version information.
log_fh.write( "$Id$" + "\n" )
log_fh.write( "$Name$" + "\n\n" )
log_fh.write( "Invoked with arguments:" )
for arg in command_line:
  if arg[0] == '-':
    log_fh.write( "\n" )
  log_fh.write( arg + ' ')

log_fh.write( "\n" )
log_fh.write( "Config file has CVS strings:\n" )
log_fh.write( cp.get('pipeline','version') + "\n" )
log_fh.write( cp.get('pipeline','cvs-tag') + "\n\n" )

print("\n===========================================\n", file=log_fh)
print("Science Segments and master chunks:\n", file=log_fh)

for ifo in ifo_list:
  print("\n===========================================\n", file=log_fh)
  print(ifo + "Data\n", file=log_fh)
  for seg in data[ifo]:
    print(" ", seg, file=log_fh)
    for chunk in seg:
      print("   ", chunk, file=log_fh)


for ifo in ifo_list:
  print("\n===========================================\n", file=log_fh)
  log_fh.write( 
    "Filtering " + str(len(chunks_analyzed[ifo])) + " " + ifo + \
    " master chunks\n" )
  total_time = 0
  for ifo_done in chunks_analyzed[ifo]:
    print(ifo_done.get_chunk(), file=log_fh)
    total_time += len(ifo_done.get_chunk())
  print("\n total time", total_time, "seconds", file=log_fh)

for ifo in ifo_list:
  print("\n===========================================\n", file=log_fh)
  log_fh.write( "Writing " + str(len(analyzed_data[ifo])) + " " + ifo + \
    " single IFO science segments\n" )
  total_time = 0
  for seg in analyzed_data[ifo]:
    print(seg, file=log_fh)
    total_time += seg.dur()
  print("\n total time", total_time, "seconds", file=log_fh)

  if opts.output_segs and len(analyzed_data[ifo]):
    if playground_only:
      f = open(ifo + '_play_segs_analyzed.txt', 'w')
    else:  
      f = open(ifo + '_segs_analyzed.txt', 'w')
    for seg in analyzed_data[ifo]:
      f.write('%4d %10d %10d %6d\n' % (seg.id(), seg.start(), seg.end(), 
        seg.dur()))
    f.close()


for ifos in ifo_coincs:  
  print("\n===========================================\n", file=log_fh)
  log_fh.write( "Writing " + str(len(analyzed_data[ifos])) + " " + ifos + \
    " coincident segments\n" )
  total_time = 0
  for seg in analyzed_data[ifos]:
    print(seg, file=log_fh)
    total_time += seg.dur()
  print("\n total time", total_time, "seconds", file=log_fh)

  if opts.output_segs and len(analyzed_data[ifos]):
    if playground_only:
      f = open(ifos + '_play_segs_analyzed.txt', 'w')
    else:  
      f = open(ifos + '_segs_analyzed.txt', 'w')
    for seg in analyzed_data[ifos]:
      f.write('%4d %10d %10d %6d\n' % (seg.id(), seg.start(), seg.end(), 
        seg.dur()))
    f.close()


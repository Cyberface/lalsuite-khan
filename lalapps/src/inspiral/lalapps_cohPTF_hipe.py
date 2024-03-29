"""
inspiral_hipe.in - standalone inspiral pipeline driver script

This script generates the condor DAG necessary to analyze LIGO and GEO
data through the inspiral pipeline.  The DAG does datafind, tmpltbank,
inspiral, inca and sire/coire steps of the pipeline.  It analyzes the
single, double, triple and quadro ifo times accordingly.  It can also be
run with injections.
"""

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
from glue.pipeline import DeepCopyableConfigParser as dcConfigParser
import urlparse
import itertools

##############################################################################
# import the modules we need to build the pipeline
from glue import pipeline
from glue import lal
from lalapps import inspiral

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
def setup_coh_inspiral(ifo_name,ifo_char,insp_job,runSplitBank,calibrated,\
                       runSpinChecker,chunk,dag,bank,scSpinBank,scNoSpinBank,sbOutBanks,scNodes,sb_node,usertag):

  # make one single inspiral job for the master chunk
  insp = inspiral.PTFInspiralNode(insp_job)
  if not opts.disable_dag_categories:
    insp.set_category('inspiral1')
    insp.set_priority(2)
  insp.set_user_tag(usertag) 
  if runSpinChecker:
    insp.set_spin_bank(scSpinBank[bank])
    insp.set_no_spin_bank(scNoSpinBank[bank])
    insp.add_parent(scNodes[bank])
  else:
    insp.set_no_spin_bank(sbOutBanks[bank])
    if runSplitBank:
      insp.add_parent(sb_node)
  insp.set_seed(bank)
  insp.set_start(chunk.start())
  insp.set_end(chunk.end())
  insp.set_trig_start(chunk.trig_start())
  insp.set_trig_end(chunk.trig_end())
  insp.set_ifo(ifo_char)
  if runSplitBank:
    insp.set_ifo_tag("FIRST_" + str(bank))
  else:
    insp.set_ifo_tag("FIRST")
  insp.set_output()
  insp.set_vds_group(ifo_name[0] + str(chunk.start()))
  if not calibrated: insp.calibration()
  output = insp.get_output()
  if opts.inspiral: dag.add_node(insp)

  return insp
##############################################################################
#function to build the timeslides vector
def setup_timeslides(ifo_analyse):
  # Initialise output 
  slide_vector = []

  # Calculate the number of segments being used
  # Note that this is the number of segments that coh_PTF will use for filtering
  # This is not equal to the number of segments that tmpltbank might use for PSD
  # estimation
  block_duration = int(cp.get('coh_PTF_inspiral','block-duration'))
  segment_duration = int(cp.get('coh_PTF_inspiral','segment-duration'))
  number_segments = (block_duration * 2 / segment_duration) - 1

  # Create a dictionary to hold offsets between pairs
  offsetPairDict = {}
  for i,ifo1 in enumerate(ifo_analyse):
    for j,ifo2 in enumerate(ifo_analyse):
      if ifo1 != ifo2 and i < j:
        offsetPairDict[ifo1+ifo2] = [0]

  # Begin looping over the possible offsets
  start = [0 for i in range(len(ifo_analyse) - 1)]
  size = [number_segments for i in range(len(ifo_analyse) - 1)]
  # I don't really know how this works, but it gives [0,0],[0,1] .. [0,len],
  # [1,0],[1,1] ... [1,len],[2,0] ... [len,0] ... [len,len] for 3 ifos,
  # for 4 it would be [0,0,0],[0,0,1] ...  and so on
  for offsetList in itertools.product(*[xrange(i, i+j) \
                                        for i,j in zip(start, size)]):
    currOffsets = [0]
    currOffsets.extend(offsetList)
    acceptSlide = True
    for i,ifo1 in enumerate(ifo_analyse):
      for j,ifo2 in enumerate(ifo_analyse):
        if ifo1 != ifo2 and i < j:
          ifoOffset = currOffsets[i] - currOffsets[j]
          if ifoOffset < 0:
            ifoOffset += number_segments
          if ifoOffset in offsetPairDict[ifo1+ifo2]:
            acceptSlide = False
            break
      if not acceptSlide:
        break
    if acceptSlide:
      # Add slide to list
      slideDict = {}
      for i,ifo in enumerate(ifo_analyse):
        slideDict[ifo] = currOffsets[i]
      slide_vector.append(slideDict)
      # And update ifo-ifo delay lists
      for i,ifo1 in enumerate(ifo_analyse):
        for j,ifo2 in enumerate(ifo_analyse):
          if ifo1 != ifo2 and i < j:
            ifoOffset = currOffsets[i] - currOffsets[j]
            if ifoOffset < 0:
              ifoOffset += number_segments
            offsetPairDict[ifo1+ifo2].append(ifoOffset)
  
  return slide_vector

##############################################################################
# function to set up datafind, template bank and inspiral jobs for an ifo  
def analyze_coh(ifo_list,ifo_data,ifo_to_do,tmplt_job,insp_job,df_job,\
  prev_df,dag, exttrigInjections, usertag=None, inspinjNode = None,\
  runSplitBank=False,sb_job = None,sbBankFile=None,sbNumBanks = None,\
  runSpinChecker=False,sc_job=None):
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

  data_opts = {}
  type = {}
  channel = {}
  ifo_char = ''

  for ifo_name in ifo_list:
    ifo_char += ifo_name
    if ifo_name == 'G1':
      data_opts['G1'] = 'geo-data'
      try: type['G1'] = cp.get('input','geo-type')
      except: type['G1'] = None
      channel['G1'] = cp.get('input','geo-channel')
    elif ifo_name == 'V1':
      data_opts['V1'] = 'virgo-data'
      try: type['V1'] = cp.get('input','virgo-type')
      except: type['V1'] = None
      channel['V1'] = cp.get('input','virgo-channel')
    else:
      data_opts[ifo_name] = 'ligo-data'
      try: 
        type[ifo_name] = cp.get('input','ligo-type')
        if (type[ifo_name] == 'RDS_R_L4') or ('RDS_C' in type[ifo_name]) or \
            ('DMT_C' in type[ifo_name]) or ('LDAS_C' in type[ifo_name]):
          type[ifo_name] = ifo_name + '_' + type[ifo_name]
      except: type[ifo_name] = None
      channel[ifo_name] = cp.get('input','ligo-channel')

  # see if we are using calibrated data
  if cp.has_section(data_opts[ifo_name]) and \
      cp.has_option(data_opts[ifo_name],'calibrated-data'):
    calibrated = True
  else: calibrated = False

  if ifo_data:
    exttrigStart = ifo_data[ifo_name][0].start()
    exttrigDuration = ifo_data[ifo_name][-1].end()-exttrigStart
    injectionFileTemplate = "HL-INJECTION_%%s-%d-%d.xml" % \
      (exttrigStart, exttrigDuration)

  # loop over the master science segments
  for seg in ifo_data[ifo_name]:

    # loop over the master analysis chunks in the science segment
    for chunk in seg:
      done_this_chunk = False

      # now loop over all the data that we need to filter
      for seg_to_do in ifo_to_do[ifo_name]:

        # if the current chunk is in one of the segments we need to filter
        if not done_this_chunk and inspiral.overlap_test(chunk,seg_to_do):

          # make sure we only filter the master chunk once
          done_this_chunk = True

          if not sbBankFile:
            # Determine template bank file name
            ifo_name = cp.get('templatebank-meta','bank-ifo')
            tb_node = inspiral.TmpltBankNode(tmplt_job)
            tb_node.set_start(chunk.start())
            tb_node.set_end(chunk.end())
            tb_node.set_ifo(ifo_name)
            tb_node.set_vds_group(ifo_name[0] + str(chunk.start()))
            tb_node.set_user_tag((usertag.split('_')[0])+'_DATAFIND')
            os.symlink("../datafind/" + tb_node.get_output(),\
                       tb_node.get_output())
            sbBankFile=tb_node.get_output()           

          # Set up the bank splitting
          if runSplitBank:
            sb_node = inspiral.SplitBankNode(sb_job)            
            sb_node.set_bank(sbBankFile)
            sb_node.set_num_banks(sbNumBanks)
            sbOutBanks = sb_node.get_output()
            dag.add_node(sb_node)
          else:
            sbNumBanks = 1

          scSpinBank = []
          scNoSpinBank = []
          scNodes = []
          if runSpinChecker:
            for bank in range(sbNumBanks):
              sc_node = inspiral.PTFSpinCheckerNode(sc_job)
              scNodes.append(sc_node)
              sc_node.set_start(chunk.start())
              sc_node.set_end(chunk.end())
              sc_node.set_ifo(ifo_char)
              if runSplitBank:
                sc_node.set_bank(sbOutBanks[bank])
                sc_node.add_parent(sb_node)
                scSpinBank.append(sbOutBanks[bank].replace('.xml','_spin.xml'))
                scNoSpinBank.append(sbOutBanks[bank].replace('.xml',\
                                                             '_nospin.xml'))
                sc_node.set_ifo_tag("FIRST_" + str(bank))
              else:
                sc_node.set_ifo_tag("FIRST")
                scSpinBank.append(sc_node.get_output_base + '_spin.xml.gz')
                scNoSpinBank.append(sc_node.get_output_base + '_nospin.xml.gz')
              sc_node.set_spin_output(scSpinBank[bank])
              sc_node.set_nospin_output(scNoSpinBank[bank])
              dag.add_node(sc_node)


          if doExtTrig:
            for inj in range(exttrigInjections[0], exttrigInjections[1]+1): 
              #XXX: ensure output is added to list of output files
              exttrigUserTag = usertag + "_" + str(inj)
              injectionFile = injectionFileTemplate % exttrigUserTag
              for bank in range(sbNumBanks):
                insp = setup_coh_inspiral(ifo_name,ifo_char,insp_job,\
                    runSplitBank,calibrated,runSpinChecker,chunk,dag,bank,\
                    scSpinBank,scNoSpinBank,sbOutBanks,scNodes,sb_node,\
                    exttrigUserTag)
                insp.set_injections(injectionFile)
          elif doSlides:
            slide_vector = setup_timeslides(ifo_analyze)
            num_slides = len(slide_vector)
            for slide in range(int(num_slides)):
              vector = slide_vector[slide]
              slidesUserTag = usertag + "_" + "slide" + "_" + \
                  str('_'.join(map(str,[str(key) + "_" + str(vector[key])\
                                        for key in vector.keys()]))) 
              for bank in range(sbNumBanks):
                insp = setup_coh_inspiral(ifo_name,ifo_char,insp_job,\
                    runSplitBank,calibrated,runSpinChecker,chunk,dag,bank,\
                    scSpinBank,scNoSpinBank,sbOutBanks,scNodes,sb_node,\
                    slidesUserTag)
                for key in vector.keys():
                  insp.add_var_opt(key.lower()+'-slide', vector[key])
          else:
            for bank in range(sbNumBanks):
              insp = setup_coh_inspiral(ifo_name,ifo_char,insp_job,\
                  runSplitBank,calibrated,runSpinChecker,chunk,dag,bank,\
                  scSpinBank,scNoSpinBank,sbOutBanks,scNodes,sb_node,usertag)

          # store this chunk in the list of filtered data 
          for ifo_name in ifo_list:
            chunks_analyzed[ifo_name] = []
            chunks_analyzed[ifo_name].append(AnalyzedIFOData(chunk,insp,0))
  
  return chunks_analyzed

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
parser.add_option("-s", "--splitbank",action="store_true",default=False,\
    help="Run splitbank to split up template banks.")
parser.add_option("-C", "--spin-checker",action="store_true",default=False,\
    help="Run spin checker on the template banks.")
parser.add_option("-i", "--inspiral" ,action="store_true",default=False,\
    help="run lalapps_inspiral (or lalapps_ring) to generate triggers")

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

if not opts.datafind and not opts.template_bank and \
     not opts.inspiral and not opts.sire_inspiral and \
     not opts.coincidence and not opts.coire_coincidence and \
     not opts.trigbank and not opts.inspiral_veto and \
     not opts.sire_inspiral_veto and not opts.td_follow_bank and \
     not opts.td_follow_inspiral and not opts.second_coinc and \
     not opts.coire_second_coinc and not opts.sire_second_coinc and \
     not opts.coherent_bank and not opts.coherent_inspiral and \
     not opts.summary_inspiral_triggers and not opts.summary_coinc_triggers: 
  print("""  No steps of the pipeline specified.
  Please specify at least one of
  --datafind, --template-bank, --inspiral, --sire-inspiral, --coincidence, 
  --coire-coincidence, --trigbank, --inspiral-veto, --sire-inspiral-veto,
  --td-follow-bank, --td-follow-inspiral, --second-coinc, 
  --coire-second-coinc, --sire-second-coinc, --coherent-bank, 
  --coherent-inspiral, --summary-inspiral-triggers --summary-coinc-triggers""", file=sys.stderr)
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
  dag.set_dag_file(basename + '.' + usertag )
else:
  dag.set_dag_file(basename )

# set better submit file names than the default
if usertag:
  subsuffix = '.' + usertag + '.sub'
else:
  subsuffix = '.sub'

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

doSlides = cp.has_option('input','do-long-slides')

tmplt_job = inspiral.TmpltBankJob(cp,opts.dax)

# spinchecker
spincheck_job = inspiral.PTFSpinCheckerJob(cp,opts.dax)

# inspiral:
insp_jobs = inspiral.PTFInspiralJob(cp,opts.dax)

if doExtTrig:
  insp_jobs.add_opt('analyze-inj-segs-only','') 

for ifo1 in ifo_list:
  if do[ifo1]:
    ifo1 = ifo1.upper()
    spincheck_job.add_opt(ifo1.lower() + '-data','')
    insp_jobs.add_opt(ifo1.lower() + '-data','')
    duration = int(cp.get('input','gps-end-time')) - int(cp.get('input','gps-start-time'))
    gps_times = cp.get('input','gps-start-time') + '-' + str(duration)
    if ifo1 == 'H1' or ifo1 == 'H2' or ifo1 == 'L1':
      spincheck_job.add_opt(ifo1.lower()+'-channel-name',ifo1.upper()+':' + cp.get('input','ligo-channel'))
      spincheck_job.add_opt(ifo1.lower()+'-frame-cache','cache/' + ifo1[0] + '-' + ifo1 + '_' + cp.get('input','ligo-type') + '_CACHE-' + gps_times + '.lcf')
      insp_jobs.add_opt(ifo1.lower()+'-channel-name',ifo1.upper()+':' + cp.get('input','ligo-channel'))
      insp_jobs.add_opt(ifo1.lower()+'-frame-cache','cache/' + ifo1[0] + '-' + ifo1 + '_' + cp.get('input','ligo-type') + '_CACHE-' + gps_times + '.lcf')

    elif ifo1 == 'V1':
      spincheck_job.add_opt(ifo1.lower()+'-channel-name',ifo1.upper()+':' + cp.get('input','virgo-channel'))
      if cp.get('input','virgo-type').startswith('T1300121'): spincheck_job.add_opt(ifo1.lower()+'-frame-cache','cache/' + ifo1[0]  + '-' + ifo1 + '_' + cp.get('input','virgo-type') + '_CACHE-' + gps_times + '.lcf')
      else: spincheck_job.add_opt(ifo1.lower()+'-frame-cache','cache/' + ifo1[0]  + '-' + cp.get('input','virgo-type') + '_CACHE-' + gps_times + '.lcf')
      insp_jobs.add_opt(ifo1.lower()+'-channel-name',ifo1.upper()+':' + cp.get('input','virgo-channel'))
      if cp.get('input','virgo-type').startswith('T1300121'): insp_jobs.add_opt(ifo1.lower()+'-frame-cache','cache/' + ifo1[0] + '-' + ifo1 + '_' + cp.get('input','virgo-type') + '_CACHE-' + gps_times + '.lcf')
      else: insp_jobs.add_opt(ifo1.lower()+'-frame-cache','cache/' + ifo1[0] + '-' + cp.get('input','virgo-type') + '_CACHE-' + gps_times + '.lcf')

for ifo1 in ifo_list:
  if do[ifo1]:
    ifo1 = ifo1.upper()
    spincheck_job.add_opt(ifo1.lower() + '-data','')
    duration = int(cp.get('input','gps-end-time')) - int(cp.get('input','gps-start-time'))
    gps_times = cp.get('input','gps-start-time') + '-' + str(duration)
    if ifo1 == 'H1' or ifo1 == 'H2' or ifo1 == 'L1':
      spincheck_job.add_opt(ifo1.lower()+'-channel-name',ifo1.upper()+':' + cp.get('input','ligo-channel'))
      spincheck_job.add_opt(ifo1.lower()+'-frame-cache','cache/' + ifo1[0] + '-' + ifo1 + '_' + cp.get('input','ligo-type') + '_CACHE-' + gps_times + '.lcf')

    elif ifo1 == 'V1':
      spincheck_job.add_opt(ifo1.lower()+'-channel-name',ifo1.upper()+':' + cp.get('input','virgo-channel'))
      if cp.get('input','virgo-type').startswith('T1300121'): spincheck_job.add_opt(ifo1.lower()+'-frame-cache','cache/' + ifo1[0]  + '-' + ifo1 + '_'  + cp.get('input','virgo-type') + '_CACHE-' + gps_times + '.lcf')
      else: spincheck_job.add_opt(ifo1.lower()+'-frame-cache','cache/' + ifo1[0]  + '-' + cp.get('input','virgo-type') + '_CACHE-' + gps_times + '.lcf')

  #splitbank
  split_job = inspiral.SplitBankJob(cp,opts.dax)
  if opts.splitbank:
    if cp.has_option('splitbank-meta','bank-file'):
      splitBankFile = cp.get('splitbank-meta','bank-file')
    else:
      splitBankFile = None
    splitNumBanks = int(cp.get('splitbank-meta','num-banks'))
    if doExtTrig:
      try:
        splitNumBanks = int(cp.get('splitbank-injmeta','num-banks'))
      except:
        pass
  else:
    splitBankFile = None
    splitNumBanks = None

all_jobs = [insp_jobs, split_job]

##############################################################################
# set the usertag in the jobs
if usertag:
  for job in all_jobs:
    if not job.get_opts().has_key('user-tag'):
      job.add_opt('user-tag',usertag)


##############################################################################
# set the condor job priority
if opts.priority:
  for job in all_jobs:
    job.add_condor_cmd('priority',str(opts.priority))

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

if play_data_mask == 'playground_only':
  playground_only = 2
elif play_data_mask == 'exclude_playground':
  playground_only = 0
elif play_data_mask == 'all_data':
  playground_only = 0
else:
  print("Invalid playground data mask " + play_data_mask + " specified")
  sys.exit(1)

##############################################################################
# get the pad and chunk lengths from the values in the ini file
pad = int(cp.get('data', 'pad-data'))
#length = int(cp.get('data','block-duration'))
#overlap = int(cp.get('data','segment-duration'))/2
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

playground_only = 0

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
sys.stdout.flush()

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

# Select on-source vs off-source (for triggered searches)
if doExtTrig:
  exttrig_analyze = cp.get('input', 'exttrig-analyze')
else:
  exttrig_analyze = None

# Read segments for external trigger analysis
exttrig_segments = pipeline.ScienceData()
if doExtTrig:
  exttrig_analyze == 'all'
#  exttrig_segments.read(cp.get('input', 'exttrig-segments'), 0)
else:
  assert (exttrig_analyze == 'all') or (exttrig_analyze is None), \
    "To use exttrig-analyze, exttrig-segments must point to a valid "\
    "segment file."


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

  # set the exttrig segments
  if exttrig_analyze is None or exttrig_analyze == "all":
    pass  # Explicitly do nothing
  elif exttrig_analyze == "on_source":
    data_out[ifo].intersection(exttrig_segments)

  elif exttrig_analyze == "off_source":
    sys.stdout.flush()
    temp = copy.deepcopy(exttrig_segments)  # don't mutate exttrig_segments
    temp.invert()
    data_out[ifo].intersection(temp)
    sys.stdout.flush()
  else:
    raise ValueError("exttrig-analyze has been set to a bad value.  It "\
      "must be {on_source | off_source | all} or omitted.")

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
# Step 4: Determine which of the master chunks needs to be filtered
chunks_analyzed = {}

prev_df = None

print("setting up jobs to filter " + ifo + " data...", end=' ')
sys.stdout.flush()

chunks_analyzed = analyze_coh(ifo_analyze,data,data_to_do,  
    tmplt_job,insp_jobs,None,None,dag,exttrigInjections,
    usertag,None,runSplitBank=opts.splitbank,sb_job = split_job,
    sbBankFile=splitBankFile,sbNumBanks = splitNumBanks,
    runSpinChecker=opts.spin_checker,sc_job=spincheck_job)

print("done") 

##############################################################################
# Step 12: Write out the LAL cache files for the various output data

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
    if isinstance(node,pipeline.LSCDataFindNode):
      # ignore datafind nodes, as their output is a cache file
      continue

    # add the data generated by the job to the output data cache
    output_file = node.get_output()

    if output_file.__class__.__name__ == 'list':
      output_data_cache.append(lal.Cache.from_urls(output_file)[0])
    else:
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
# Step 13: Setup the Maximum number of jobs for different categories
if not opts.disable_dag_categories:
  for cp_opt in cp.options('condor-max-jobs'):
    dag.add_maxjobs_category(cp_opt,cp.getint('condor-max-jobs',cp_opt))

# Add number of retries to jobs as specified by ihope.ini
if cp.has_option("pipeline", "retry-jobs"):
  num_retries = cp.getint("pipeline", "retry-jobs")
  for node in dag.get_nodes():
    node.set_retry(num_retries)

##############################################################################
# Step 14: Write out the DAG, help message and log file
dag.write_sub_files()
dag.write_dag()
if opts.dax:
  pass
  #print "writing rls cache" 
  dag.write_pegasus_rls_cache(cp.get('ldgsubmitdax','gsiftp'),cp.get('ldgsubmitdax','pool'))

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


for ifo in ifo_analyze:
  print("\n===========================================\n", file=log_fh)
  log_fh.write( 
    "Filtering " + str(len(chunks_analyzed[ifo])) + " " + ifo + \
    " master chunks\n" )
  total_time = 0
  for ifo_done in chunks_analyzed[ifo]:
    print(ifo_done.get_chunk(), file=log_fh)
    total_time += len(ifo_done.get_chunk())
  print("\n total time", total_time, "seconds", file=log_fh)

for ifo in ifo_analyze:
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

sys.exit(0)


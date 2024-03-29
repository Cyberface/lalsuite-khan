"""

MakeSFTDAG.py - Creates DAGs to run jobs that generates SFTs; can act as a dag generator for use with onasys.
Loosely based on lalapps_strain_pipe


"""

from __future__ import print_function

__author__ = 'Greg Mendell<gmendell@ligo-wa.caltech.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

# REVISIONS:
# 12/02/05 gam; generate datafind.sub and MakeSFTs.sub as well as dag file in PWD, with log files based subLogPath and dag filename. 
# 12/28/05 gam; Add option --make-gps-dirs, -D <num>, to make directory based on this many GPS digits.
# 12/28/05 gam; Add option --misc-desc, -X <string> giving misc. part of the SFT description field in the filename.
# 12/28/05 gam; Add options --start-freq -F and --band -B options to enter these.
# 12/28/05 gam; Add in --window-type, -w options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window.
# 12/28/05 gam; Add option --overlap-fraction -P (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 0.0)
# 12/28/05 gam; Add --sft-version, -v option to select output SFT version (1 = default is version 1 SFTs; 2 = version 2 SFTs.
# 12/28/05 gam; Add --comment-field, -c option, for comment for version 2 SFTs.
# 12/28/05 gam; Remove sample rate option
# 01/09/06 gam; Add -Z option; write SFT to .*.tmp file, then move to final file name.
# 01/14/07 gam; Add -u option to specify frame struct and type; add -i option to specify IFO name.
# 07/24/07 gam; Add in -q option to read in list of nodes on which to output SFTs, -Q option to give node path, and -R option for number of jobs per node.
# 04/XX/13 eag; Add -y option to synchronize the start times of SFTs.
# 07/24/14 eag; Change default to version 2 SFTs

# import standard modules and append the lalapps prefix to the python path
import sys, os
import getopt, re, string
import tempfile
import math
#import ConfigParser
#sys.path.append('')

# import the modules we need to build the pipeline
#from glue import pipeline
#import strain

#
# USAGE FUNCTION
#
def usage():
  msg = """\

This script creates datafind.sub, MakeSFTs.sub, and a dag file that generates SFTs based on the options given. 

The script can be used to create dag files for stand-alone use with condor_submit_dag, or as a dag generator with onasys.
  
Usage: MakeSFTDAG [options]

  -h, --help                 display this message
  -s, --gps-start-time       (optional) GPS start time of data segment (unused except by onasys)
  -e, --gps-end-time         (optional) GPS end time of data segment (unused except by onasys)
  -a, --analysis-start-time  GPS start time of data from which to generate SFTs (optional and unused if a segment file is given)
  -b, --analysis-end-time    GPS end time of data from which to generate SFTs (optional and unused if a segment file is given)
  -f, --dag-file             filename for .dag file (should end in .dag)
  -t, --aux-path             (optional) path to auxiliary data files (unused except by onasys)
  -G, --tag-string           tag string used in names of various files unique to jobs that will run under the DAG
  -d, --input-data-type      input data type for use with the LSCdataFind --type option
  -x, --extra-datafind-time  (optional) extra time to subtract/add from/to start/end time arguments of LSCdataFind (default is 0)
  -M, --datafind-match       (optional) string to use with the LSCdataFind --match option
  -y, --synchronize-start    (optional) synchronize the start times of the SFTs so that the start times are synchronized when there are gaps in the data (0 or 1, default is 0)
  -k, --filter-knee-freq     high pass filter knee frequency used on time domain data before generating SFTs
  -T, --time-baseline        time baseline of SFTs  (e.g., 60 or 1800 seconds)
  -p, --output-sft-path      path to output SFTs
  -C, --cache-path           (optional) path to cache files that will be produced by LSCdataFind (default is $PWD/cache; this directory is created if it does not exist and must agree with that given in .sub files)
  -O, --log-path             (optional) path to some log, output, and error files (default is $PWD/logs; this directory is created if it does not exist and should agree with that given in .sub files)
  -o, --sub-log-path         (optional) path to log file to give in datafind.sub and MakeSFTs.sub (default is $PWD/logs; this directory must exist and usually should be under a local file system.)
  -N, --channel-name         name of input time-domain channel to read from frames
  -i, --ifo                  (optional) Name of IFO, i.e., H1, H2, L1, or G1; use if channel name begins with H0, L0, or G0; default: use first two characters from channel name.
  -v, --sft-version          (optional) sft version to output (1 or 2; default is 2)
  -c, --comment-field        (optional) comment for version 2 SFT header.
  -F, --start-freq           (optional) start frequency of the SFTs (default is 48 Hz).
  -B, --band                 (optional) frequency band of the SFTs (default is 2000 Hz).
  -D  --make-gps-dirs        (optional) make directories for output SFTs based on this many digits of the GPS time.
  -Z  --make-tmp-file        (optional) write SFT to .*.tmp file, then move to final filename.
  -X  --misc-desc            (optional) misc. part of the SFT description field in the filename (also used if -D option is > 0).
  -w, --window-type          (optional) type of windowing of time-domain to do before generating SFTs (1 = Tukey given by Matlab, 2 = Tukey given in lalapps/src/pulsar/make_sfts.c, 3 = Hann given by Matlab; default is 1)
  -P, --overlap-fraction     (optional) overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 0.0).
  -m, --max-num-per-node     maximum number of SFTs to generate on one node
  -L, --max-length-all-jobs  maximum total amount of data to process, in seconds (optional and unused if a segment file is given)
  -g, --segment-file         (optional) alternative file with segments to use, rather than the input times.
  -A, --accounting-group     (optional) accounting group tag to be added to the condor submit files.
  -U, --accounting-group-user  (optional) accounting group albert.einstein username to be added to the condor submit files.  
  -q, --list-of-nodes        (optional) file with list of nodes on which to output SFTs.
  -Q, --node-path            (optional) path to nodes to output SFTs; the node name is appended to this path, followed by path given by the -p option;
                                        for example, if -q point to file with the list node1 node2 ... and the -Q /data/ -p /frames/S5/sfts/LHO options
                                        are given, the first output file will go into /data/node1/frames/S5/sfts/LHO; the next node in the
                                        list is used in constructing the path when the number of jobs given by the -R option reached, and so on.
  -R, --output-jobs-per-node (optional) number of jobs to output per node in the list of nodes given with the -q option.
  -l, --min-seg-length       (optional) minimum length segments to process in seconds (used only if a segment file is given; default is 0)
  -S, --use-single           (optional) use single precision in MakeSFTs for windowing, plan and fft; filtering is always done in double precision. Use of double precision in MakeSFTs is the default.
  -H, --use-hot              (optional) input data is from h(t) calibrated frames (h of t = hot!) (0 or 1).
  -u  --frame-struct-type    (optional) string specifying the input frame structure and data type. Must begin with ADC_ or PROC_ followed by REAL4, REAL8, INT2, INT4, or INT8; default: ADC_REAL4; -H is the same as PROC_REAL8.
  -j, --datafind-path        (optional) string specifying a path to look for the gw_data_find executable; if not set, will use LSC_DATAFIND_PATH env variable or system default (in that order).
  -J, --makesfts-path        (optional) string specifying a path to look for the lalapps_MakeSFTs executable; if not set, will use MAKESFTS_PATH env variable or system default (in that order).
  -Y, --request-memory       (optional) memory allocation in MB to request from condor for lalapps_MakeSFTs step
"""
  print(msg)

#
# FUNCTION THAT WRITE ONE JOB TO DAG FILE
#
def writeToDag(dagFID, nodeCount, filterKneeFreq, timeBaseline, outputSFTPath, cachePath, startTimeThisNode, endTimeThisNode, channelName, site, inputDataType, extraDatafindTime, useSingle, useHoT, makeTmpFile, tagString, windowType, overlapFraction, sftVersion, makeGPSDirs, miscDesc, commentField, startFreq, freqBand, frameStructType, IFO):
  LSCdataFind = 'LSCdataFind_%i' % nodeCount
  MakeSFTs    = 'MakeSFTs_%i' % nodeCount
  startTimeDatafind = startTimeThisNode - extraDatafindTime
  endTimeDatafind = endTimeThisNode + extraDatafindTime
  tagStringOut = '%s_%i' % (tagString, nodeCount)
  cacheFile   = '%s/%s-%s-%s.cache' % (cachePath, site,startTimeDatafind,endTimeDatafind)
  dagFID.write('JOB %s datafind.sub\n' % LSCdataFind)
  dagFID.write('RETRY %s 10\n' % LSCdataFind)
  dagFID.write('VARS %s gpsstarttime="%s" gpsendtime="%s" observatory="%s" inputdatatype="%s" tagstring="%s"\n' % (LSCdataFind, startTimeDatafind, endTimeDatafind, site, inputDataType, tagStringOut))
  dagFID.write('JOB %s MakeSFTs.sub\n' % MakeSFTs)
  dagFID.write('RETRY %s 5\n' % MakeSFTs)
  argList = '-f %s -t %s -p %s -C %s -s %s -e %s -N %s' %(filterKneeFreq, timeBaseline, outputSFTPath, cacheFile, startTimeThisNode, endTimeThisNode, channelName)
  argList = argList + ' -v %s' % sftVersion
  if IFO != None: argList = argList + ' -i %s' % IFO
  if commentField != None: argList = argList + ' -c %s' % commentField
  if frameStructType != None: argList = argList + ' -u %s' % frameStructType
  if startFreq != 48.0: argList = argList + ' -F %s' % startFreq
  if freqBand != 2000.0: argList = argList + ' -B %s' % freqBand
  if makeGPSDirs != 0: argList = argList + ' -D %s' % makeGPSDirs
  if miscDesc != None: argList = argList + ' -X %s' % miscDesc
  if windowType != 1: argList = argList + ' -w %s' % windowType
  if overlapFraction != 0.0: argList = argList + ' -P %s' % overlapFraction
  if useSingle: argList = argList + ' -S'
  if useHoT: argList = argList + ' -H'
  if makeTmpFile: argList = argList + ' -Z'
  dagFID.write('VARS %s argList="%s" tagstring="%s"\n'%(MakeSFTs,argList, tagStringOut))
  dagFID.write('PARENT %s CHILD %s\n'%(LSCdataFind,MakeSFTs))
  
#
# MAIN CODE START HERE 
#

# parse the command line options
shortop = "s:e:a:b:f:t:G:d:x:M:y:k:T:p:C:O:o:N:i:w:P:u:v:c:F:B:D:X:m:L:g:q:Q:A:U:R:l:hSHZj:J:Y:"
longop = [
  "help",
  "gps-start-time=",
  "gps-end-time=",
  "analysis-start-time=",
  "analysis-end-time=",
  "dag-file=",
  "aux-path=",
  "tag-string=",
  "input-data-type=",
  "extra-datafind-time=",
  "datafind-match=",
  "synchronize-start=",
  "filter-knee-freq=",
  "time-baseline=",
  "output-sft-path=",
  "cache-path=",
  "log-path=",
  "sub-log-path=",
  "channel-name=",
  "ifo=",
  "window-type=",
  "overlap-fraction=",
  "frame-struct-type=",
  "sft-version=",
  "comment-field=",
  "start-freq=",
  "band=",
  "make-gps-dirs=",
  "misc-desc=",
  "max-num-per-node=",
  "max-length-all-jobs=",
  "segment-file=",
  "accounting-group=",
  "accounting-group-user=",
  "list-of-nodes=",
  "node-path=",
  "output-jobs-per-node=",
  "min-seg-length=",
  "use-single=",  
  "use-hot",
  "make-tmp-file",
  "datafind-path=",
  "makesfts-path=",
  "request-memory="
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

# initialize with default values or with None
gpsStartTime = None
gpsEndTime = None
analysisStartTime = None
analysisEndTime = None
dagFileName = None
auxPath = "."
tagString = None
inputDataType = None
extraDatafindTime = 0
datafindMatch = None
synchronizeStart = 0
filterKneeFreq = -1
timeBaseline = None
outputSFTPath = None
cachePath = "cache"
logPath = "logs"
subLogPath = "logs"
channelName = None
IFO = None
windowType = 1
overlapFraction = 0.0
sftVersion = 2
commentField = None
frameStructType = None
startFreq = 48.0
freqBand = 2000.0
makeGPSDirs = 0
miscDesc = None
maxNumPerNode = None
maxLengthAllJobs = None
segmentFile = None
nodeListFile = None
accountingGroup = None
accountingGroupUser = None
nodePath = None
outputJobsPerNode = 0
useNodeList = False
nodeListIndex = 0
savedOutputSFTPath = None
minSegLength = 0
useSingle = False
useHoT = False
makeTmpFile = False
datafindPath = None
makeSFTsPath = None
requestMemory = None

for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-s", "--gps-start-time"):
    gpsStartTime = int(a)
  elif o in ("-e", "--gps-end-time"):
    gpsEndTime = int(a)
  elif o in ("-a", "--analysis-start-time"):
    analysisStartTime = int(a)
  elif o in ("-b", "--analysis-end-time"):
    analysisEndTime = int(a)    
  elif o in ("-f", "--dag-file"):
    dagFileName = a   
  elif o in ("-t", "--aux-path"):
    auxPath = a
  elif o in ("-G", "--tag-string"):
    tagString = a
  elif o in ("-d", "--input-data-type"):
    inputDataType = a
  elif o in ("-x", "--extra-datafind-time"):
    extraDatafindTime = int(a)
  elif o in ("-M", "--datafind-match"):
    datafindMatch = a
  elif o in ("-y", "--synchronize-start"):
    synchronizeStart = int(a)
  elif o in ("-k", "--filter-knee-freq"):
    filterKneeFreq = int(a)
  elif o in ("-T", "--time-baseline"):
    timeBaseline = int(a)
  elif o in ("-p", "--output-sft-path"):
    outputSFTPath = a
  elif o in ("-C", "--cache-path"):
    cachePath = a
  elif o in ("-O", "--log-path"):
    logPath = a
  elif o in ("-o", "--sub-log-path"):
    subLogPath = a
  elif o in ("-N", "--channel-name"):
    channelName = a
  elif o in ("-i", "--ifo"):
    IFO = a
  elif o in ("-w", "--window-type"):
    windowType = int(a)
  elif o in ("-P", "--overlap-fraction"):
    overlapFraction = float(a)
  elif o in ("-u", "--frame-struct-type"):
    frameStructType = a
  elif o in ("-v", "--sft-version"):
    sftVersion = int(a)
  elif o in ("-c", "--comment-field"):
    commentField = a
  elif o in ("-F", "--start-freq"):
    startFreq = int(a)
  elif o in ("-B", "--band"):
    freqBand = int(a)
  elif o in ("-D", "--make-gps-dirs"):
    makeGPSDirs = int(a)
  elif o in ("-X", "--misc-desc"):
    miscDesc = a
  elif o in ("-m", "--max-num-per-node"):
    maxNumPerNode = int(a)
  elif o in ("-L", "--max-length-all-jobs"):
    maxLengthAllJobs = int(a)
  elif o in ("-g", "--segment-file"):
    segmentFile = a
  elif o in ("-q", "--list-of-nodes"):
    nodeListFile = a
  elif o in ("-A", "--accounting-group"):
    accountingGroup = a
  elif o in ("-U", "--accounting-group-user"):
    accountingGroupUser = a
  elif o in ("-Q", "--node-path"):
    nodePath = a
  elif o in ("-R", "--output-jobs-per-node"):
    outputJobsPerNode = int(a)
  elif o in ("-l", "--min-seg-length"):
    minSegLength = int(a)
  elif o in ("-S", "--use-single"):
    useSingle = True    
  elif o in ("-H", "--use-hot"):
    useHoT = True
  elif o in ("-Z", "--make-tmp-file"):
    makeTmpFile = True
  elif o in ("-j", "--datafind-path"):
    datafindPath = a
  elif o in ("-J", "--makesfts-path"):
    makeSFTsPath = a
  elif o in ("-Y", "--request-memory"):
    requestMemory = a
  else:
    print("Unknown option:", o, file=sys.stderr)
    usage()
    sys.exit(1)

if not dagFileName:
  print("No dag filename specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if not tagString:
  print("No tag string specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if not inputDataType:
  print("No input data type specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if extraDatafindTime < 0:
  print("Invalid extra datafind time specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if synchronizeStart<0 or synchronizeStart>1:
  print("Invalid use of synchronize-start.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)
  
if filterKneeFreq < 0:
  print("No filter knee frequency specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if not timeBaseline:
  print("No time baseline specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if not outputSFTPath:
  print("No output SFT path specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)
  
if not cachePath:
  print("No cache path specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if not logPath:
  print("No log path specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if not subLogPath:
  print("No sub log path specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if not channelName:
  print("No channel name specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if (windowType != 0) and (windowType != 1) and (windowType != 2) and (windowType != 3):
  print("Invalid window type specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if (overlapFraction < 0.0) or (overlapFraction >= 1.0):
  print("Invalid make overlap fraction specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if (sftVersion != 1) and (sftVersion != 2):
  print("Invalid SFT version specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if (startFreq < 0.0):
  print("Invalid start freq specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if (freqBand < 0.0):
  print("Invalid band specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if (makeGPSDirs < 0) or (makeGPSDirs > 10):
  print("Invalid make gps dirs specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

if not maxNumPerNode:
  print("No maximum number of SFTs per node specified.", file=sys.stderr)
  print("Use --help for usage details.", file=sys.stderr)
  sys.exit(1)

dataFindExe = 'gw_data_find'
if datafindPath:
    dataFindExe = os.path.join(datafindPath,dataFindExe)
elif 'LSC_DATAFIND_PATH' in os.environ:
    dataFindExe = os.path.join('$ENV(LSC_DATAFIND_PATH)',dataFindExe)

makeSFTsExe = 'lalapps_MakeSFTs'
if makeSFTsPath:
    makeSFTsExe = os.path.join(makeSFTsPath,makeSFTsExe)
elif 'MAKESFTS_PATH' in os.environ:
    makeSFTsExe = os.path.join('$ENV(MAKESFTS_PATH)',makeSFTsExe)

# try and make a directory to store the cache files and job logs
try: os.mkdir(logPath)
except: pass
try: os.mkdir(cachePath)
except: pass

# Check if list of nodes is given, on which to output SFTs.
nodeList = [];
if (nodeListFile != None):

    if not nodePath:
       print("Node file list given, but no node path specified.", file=sys.stderr)
       sys.exit(1)
    
    if (outputJobsPerNode < 1):
       print("Node file list given, but invalid output jobs per node specified.", file=sys.stderr)
       sys.exit(1)

    try:
         for line in open(nodeListFile):
             # split the line get rid of the \n at the end
             splitLine = line.split()
             nodeList.append(splitLine[0])
         # End for line in open(nodeListFile)
         if (len(nodeList) < 1):
             print("No nodes found in node list file: %s." % nodeListFile, file=sys.stderr)
             sys.exit(1)
    except:
         print("Error reading or parsing node list file: %s." % nodeListFile, file=sys.stderr)
         sys.exit(1)

    # Set flag to use list of nodes in constructing output files
    useNodeList = True
    savedOutputSFTPath = outputSFTPath
# END if (nodeListFile != None)

# Check if segment file was given, else set up one segment from the command line
segList = [];
adjustSegExtraTime = False
if (segmentFile != None):

    if minSegLength < 0:
      print("Invalid minimum segment length specified.", file=sys.stderr)
      print("Use --help for usage details.", file=sys.stderr)
      sys.exit(1)

    # the next flag causes extra time that cannot be processes to be trimmed from the start and end of a segment
    adjustSegExtraTime = True
    try:
         for line in open(segmentFile):
             try: 
                 splitLine = line.split();
                 try: 
                     oneSeg = [];
                     oneSeg.append(int(splitLine[0]));
                     oneSeg.append(int(splitLine[1]));
                     if ((oneSeg[1] - oneSeg[0]) >= minSegLength):
                         segList.append(oneSeg)
                     else:
                         pass
                 except:
                     pass
             except:
                 pass
         # End for line in open(segmentFile)
         if (len(segList) < 1):
             print("No segments found in segment file: %s." % segmentFile, file=sys.stderr)
             sys.exit(1)
    except:
         print("Error reading or parsing segment file: %s." % segmentFile, file=sys.stderr)
         sys.exit(1)
else:
    if not analysisStartTime:
      print("No GPS analysis start time specified.", file=sys.stderr)
      print("Use --help for usage details.", file=sys.stderr)
      sys.exit(1)

    if not analysisEndTime:
      print("No GPS analysis end time specified.", file=sys.stderr)
      print("Use --help for usage details.", file=sys.stderr)
      sys.exit(1)

    if not maxLengthAllJobs:
      print("No maximum length of all jobs specified.", file=sys.stderr)
      print("Use --help for usage details.", file=sys.stderr)
      sys.exit(1)

    # Make sure not to exceed maximum allow analysis
    if (analysisEndTime > analysisStartTime + maxLengthAllJobs):
       analysisEndTime = analysisStartTime + maxLengthAllJobs

    try:
        oneSeg = [];
        oneSeg.append(analysisStartTime);
        oneSeg.append(analysisEndTime);
        segList.append(oneSeg);
    except:
        print("There was a problem setting up the segment to run on: [%s, %s)." % (analysisStartTime,analysisEndTime), file=sys.stderr)
        sys.exit(1)
    
# END if (segmentFile != None)

# Get the IFO site, which is the first letter of the channel name.
site = channelName[0]

# initialize count of nodes
nodeCount         = 0

# create datafind.sub
datafindFID = file('datafind.sub','w')
datafindLogFile = subLogPath + '/' + 'datafind_' + dagFileName + '.log'
datafindFID.write('universe = vanilla\n')
datafindFID.write('executable = ' + dataFindExe + '\n')
if not datafindMatch:
   dataFindMatchString = ''
else:
   dataFindMatchString = '--match ' + datafindMatch
datafindFID.write('arguments = -r $ENV(LIGO_DATAFIND_SERVER) --observatory $(observatory) --url-type file --gps-start-time $(gpsstarttime) --gps-end-time $(gpsendtime) --lal-cache --type $(inputdatatype) %s\n' % dataFindMatchString)
datafindFID.write('getenv = True\n')
if (accountingGroup != None):
   datafindFID.write('accounting_group = %s\n' % accountingGroup)
if (accountingGroupUser != None):
   datafindFID.write('accounting_group_user = %s\n' % accountingGroupUser)
datafindFID.write('log = %s\n' % datafindLogFile)
datafindFID.write('error = %s/datafind_$(tagstring).err\n' % logPath)
datafindFID.write('output = %s/$(observatory)-$(gpsstarttime)-$(gpsendtime).cache\n' % cachePath)
datafindFID.write('notification = never\n')
datafindFID.write('queue 1\n')
datafindFID.close

# create MakeSFTs.sub
MakeSFTsFID = file('MakeSFTs.sub','w')
MakeSFTsLogFile = subLogPath + '/' + 'MakeSFTs_' + dagFileName + '.log'
MakeSFTsFID.write('universe = vanilla\n')
MakeSFTsFID.write('executable = '+ makeSFTsExe + '\n')
MakeSFTsFID.write('arguments = $(argList)\n')
MakeSFTsFID.write('getenv = True\n')
if (accountingGroup != None):
   MakeSFTsFID.write('accounting_group = %s\n' % accountingGroup)
if (accountingGroupUser != None):
   MakeSFTsFID.write('accounting_group_user = %s\n' % accountingGroupUser)
MakeSFTsFID.write('log = %s\n' % MakeSFTsLogFile)
MakeSFTsFID.write('error = %s/MakeSFTs_$(tagstring).err\n' % logPath)
MakeSFTsFID.write('output = %s/MakeSFTs_$(tagstring).out\n' % logPath)
MakeSFTsFID.write('notification = never\n')
if (requestMemory != None):
    MakeSFTsFID.write('RequestMemory = %s\n' % requestMemory)
MakeSFTsFID.write('RequestCpus = 1\n')
MakeSFTsFID.write('queue 1\n')
MakeSFTsFID.close

# create the DAG file with the jobs to run
#dagFileName = dagFileName + '.dag'
dagFID = file(dagFileName,'w')

startTimeAllNodes = None
firstSFTstartTime = 0
for seg in segList:
    # Each segment in the segList runs on one or more nodes; initialize the number SFTs produced by the current node:
    numThisNode = 0
    numThisSeg = 0
    if (adjustSegExtraTime and synchronizeStart==0):
       segStartTime = seg[0]
       segEndTime = seg[1]
       segExtraTime = (segEndTime - segStartTime) % timeBaseline
       if overlapFraction != 0.0:
          # handle overlap
          if (segEndTime - segStartTime) > timeBaseline:
             segExtraTime = (segEndTime - segStartTime - timeBaseline) % int((1.0 - overlapFraction)*timeBaseline)
          # else there will just one SFT this segment
       else:
          # default case, no overlap
          segExtraTime = (segEndTime - segStartTime) % timeBaseline
       segExtraStart =  int(segExtraTime / 2)
       segExtraEnd = segExtraTime - segExtraStart
       #print segStartTime,segEndTime, segExtraTime, segExtraStart, segExtraEnd
       analysisStartTime = segStartTime + segExtraStart
       if analysisStartTime > segEndTime: analysisStartTime = segEndTime
       analysisEndTime = segEndTime - segExtraEnd
       if analysisEndTime < segStartTime: analysisEndTime = segStartTime
    elif (synchronizeStart):
       segStartTime = seg[0]
       segEndTime = seg[1]
       if firstSFTstartTime == 0: firstSFTstartTime = segStartTime
       analysisStartTime = int(round(math.ceil((segStartTime - firstSFTstartTime)/((1.0 - overlapFraction)*timeBaseline))*(1.0 - overlapFraction)*timeBaseline)) + firstSFTstartTime
       if analysisStartTime > segEndTime: analysisStartTime = segEndTime
       analysisEndTime = int(round(math.floor((segEndTime - analysisStartTime - timeBaseline)/((1.0 - overlapFraction)*timeBaseline))*(1.0 - overlapFraction)*timeBaseline)) + timeBaseline + analysisStartTime
       if analysisEndTime < segStartTime: analysisEndTime = segStartTime
    else:
       analysisStartTime = seg[0]
       analysisEndTime = seg[1]     
    #print analysisStartTime, analysisEndTime
    # Loop through the analysis time; make sure no more than maxNumPerNode SFTs are produced by any one node
    startTimeThisNode = analysisStartTime
    endTimeThisNode   = analysisStartTime
    endTimeAllNodes   = analysisStartTime
    while (endTimeAllNodes < analysisEndTime):
         # increment endTimeAllNodes by the timeBaseline until we get past the analysisEndTime
         if overlapFraction != 0.0:
            # handle overlap
            if numThisSeg == 0:
               endTimeAllNodes = endTimeAllNodes + timeBaseline
            else:
               endTimeAllNodes = endTimeAllNodes + int((1.0 - overlapFraction)*timeBaseline)
         else:
            # default case, no overlap
            endTimeAllNodes = endTimeAllNodes + timeBaseline
         if (endTimeAllNodes <= analysisEndTime):
            # increment the number of SFTs output from this node, and update the end time this node.
            numThisNode = numThisNode + 1
            numThisSeg = numThisSeg + 1
            endTimeThisNode = endTimeAllNodes
            if (numThisNode < maxNumPerNode):
               continue
            else:
               # write jobs to dag for this node
               nodeCount = nodeCount + 1

               if (useNodeList):
                  outputSFTPath = nodePath + nodeList[nodeListIndex] + savedOutputSFTPath
                  if ((nodeCount % outputJobsPerNode) == 0):
                     nodeListIndex = nodeListIndex + 1
                  # END if ((nodeCount % outputJobsPerNode) == 0L)
               # END if (useNodeList)

               if (nodeCount == 1): startTimeAllNodes = startTimeThisNode
               writeToDag(dagFID,nodeCount, filterKneeFreq, timeBaseline, outputSFTPath, cachePath, startTimeThisNode, endTimeThisNode, channelName, site, inputDataType, extraDatafindTime, useSingle, useHoT, makeTmpFile, tagString, windowType, overlapFraction, sftVersion, makeGPSDirs, miscDesc, commentField, startFreq, freqBand, frameStructType, IFO)
               # Update for next node
               numThisNode       = 0
               if overlapFraction != 0.0:
                  # handle overlap
                  startTimeThisNode = endTimeThisNode - int((overlapFraction)*timeBaseline)
               else:
                  # default case, no overlap
                  startTimeThisNode = endTimeThisNode
    else:
         # we are at or past the analysisEndTime; output job for last node if needed.
         if (numThisNode > 0):
            # write jobs to dag for this node
            nodeCount = nodeCount + 1

            if (useNodeList):
               outputSFTPath = nodePath + nodeList[nodeListIndex] + savedOutputSFTPath
               if ((nodeCount % outputJobsPerNode) == 0):
                  nodeListIndex = nodeListIndex + 1
               # END if ((nodeCount % outputJobsPerNode) == 0L)
            # END if (useNodeList)

            if (nodeCount == 1): startTimeAllNodes = startTimeThisNode
            writeToDag(dagFID,nodeCount, filterKneeFreq, timeBaseline, outputSFTPath, cachePath, startTimeThisNode, endTimeThisNode, channelName, site, inputDataType, extraDatafindTime, useSingle, useHoT, makeTmpFile, tagString, windowType, overlapFraction, sftVersion, makeGPSDirs, miscDesc, commentField, startFreq, freqBand, frameStructType, IFO)
    # END while (endTimeAllNodes < analysisEndTime)
# END for seg in segList

# Close the DAG file
dagFID.close

# Update actual end time of the last job and print out the times all jobs will run on:
endTimeAllNodes = endTimeThisNode

if not startTimeAllNodes:
  print("The startTimeAllNodes == none; the DAG file contains no jobs!", file=sys.stderr)
  sys.exit(1)

if (endTimeAllNodes <= startTimeAllNodes):
  print("The endTimeAllNodes <= startTimeAllNodes; the DAG file contains no jobs!", file=sys.stderr)
  sys.exit(1)

print(startTimeAllNodes, endTimeAllNodes)

sys.exit(0)

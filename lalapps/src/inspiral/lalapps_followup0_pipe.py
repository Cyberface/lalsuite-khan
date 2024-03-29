"""
Something

This program creates cache files for the output of inspiral hipe
"""

from __future__ import print_function

__author__ = 'Chad Hanna <channa@ligo.caltech.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

##############################################################################
# import standard modules
import sys, os, copy, math, random
from subprocess import *
import socket, time
import re, string
from optparse import *
import tempfile
from six.moves import configparser
import urlparse
import urllib
from UserDict import UserDict

##############################################################################
# import the modules we need to build the pipeline
from ligo import segments
from ligo.segments import utils as segmentsUtils
from glue import pipeline
from glue import lal
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import CoincInspiralUtils
from pylal.fu_utils import *
from pylal.fu_writeXMLparams import * 
from pylal.fu_Condor import *
from pylal.webUtils import *
from pylal import Fr
from lalapps import inspiral
from lalapps import inspiralutils

##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################

######################## OPTION PARSING  #####################################
usage = """usage: %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")

parser.add_option("-l", "--log-path",action="store",type="string",\
    metavar=" PATH",help="directory to write condor log file")

parser.add_option("-f", "--config-file",action="store",type="string",\
    metavar=" FILE",help="ini file")

parser.add_option("-g", "--generate-fu-cache",action="store_true",\
    default=False, help="create hipe caches")

parser.add_option("", "--read-times", action="store_true",\
    default=False, help="If this option is activated, the pipeline \
    will determine the list of triggers to be analysed from text files \
    containing the list of GPS times per interferometers.\ The paths to these \
    text files are expected to be provided within the fields \"XXtimes\" \
    (where XX is an ifo name) in the section [triggers] of the .ini file.")

parser.add_option("-D", "--disable-followup",action="store_true",\
    default=False, help="disable candidate followup (for running qscans only)")

parser.add_option("-m", "--datafind",action="store_true",\
    default=False, help="use datafind to get qscan/trends data")

parser.add_option("-M", "--hoft-datafind",action="store_true",\
    default=False, help="use datafind to get hoft data (for qscan)")

parser.add_option("", "--single-qevent",action="store_true",default=False,\
    help="do single qevent")

parser.add_option("", "--H1H2-qevent",action="store_true",default=False,\
    help="do coherent qevent")

parser.add_option("-q", "--qscan",action="store_true",default=False,\
    help="do qscans")

parser.add_option("-G", "--generate-segments",action="store_true",\
    default=False, help="generate the science segments for background qscans")

parser.add_option("-Q", "--background-qscan",action="store_true",\
    default=False, help="do qscans over a list of times")

parser.add_option("-n", "--hoft-qscan",action="store_true",\
    default=False, help="do hoft qscans")

parser.add_option("-N", "--background-hoft-qscan",action="store_true",\
    default=False, help="do hoft qscans over a list of times")

parser.add_option("-s", "--seis-qscan",action="store_true",\
    default=False, help="do seismic qscans")

parser.add_option("-S", "--background-seis-qscan",action="store_true",\
    default=False, help="do seismic qscans over a list of times")

parser.add_option("", "--distrib-remote-q", action="store_true",\
    default=False, help="this argument needs to be called when you want to \
    analyse the V1 qscans and seismic qscans that have been run remotely at  \
    CC-Lyon. After the V1 qscans have been done, you have received a \
    \"V1_qscans_results.tar.gz\" file that contains all the qscan results. \
    Before launching the \"analyseQscan\" jobs, you first need to copy these \
    results at their expected location, ie the \"xxouput\" fields defined in \
    the qscan sections of your \"followup_pipe.ini\" file. This is taken care \
    of by the argument \"--distrib-remote-q\"")

parser.add_option("-a", "--analyse-qscan",action="store_true",\
    default=False, help="run the analyseQscan script to interprete the qscan")

parser.add_option("-b", "--analyse-seis-qscan",action="store_true",\
    default=False, help="run the analyseQscan script to interprete the seismic qscan")

parser.add_option("-e", "--analyse-hoft-qscan",action="store_true",\
    default=False, help="run the analyseQscan script to interprete the hoft qscan")

parser.add_option("-d", "--data-quality",action="store_true",default=False,\
    help="do data quality lookup - CURRENTLY BROKEN, DO NOT USE IT")

parser.add_option("-t", "--trig-bank",action="store_true",default=False,\
    help="generate a pseudo trigbank xml file for single triggers")

parser.add_option("", "--inspiral-datafind",action="store_true",default=False,\
    help="do datafind jobs to update the old cache files." + \
         "This option affects the way the inspiral, frame-check, and mcmc jobs" + \
         "accessing to the data.")

parser.add_option("-i", "--inspiral",action="store_true",default=False,\
    help="do inspirals for single triggers - output SNR/CHISQ/PSD")

parser.add_option("-F", "--frame-check",action="store_true",default=False,\
    help="do FrCheck jobs on the frames used by the inspiral code")

parser.add_option("-I", "--ifo-status-check",action="store_true",default=False,\
    help="download IFO summary plots")

parser.add_option("-p", "--plots",action="store_true",default=False,\
    help="plot SNR/CHISQ from inspiral stage")

parser.add_option("-C", "--mcmc",action="store_true",default=False,\
    help="Do MCMC on the followup trigs (experimental)")

parser.add_option("-P", "--plot-mcmc",action="store_true",default=False,\
    help="Plot MCMC results (experimental)")

parser.add_option("", "--spin-mcmc",action="store_true",default=False,\
    help="Do SPIN MCMC on the followup trigs (experimental)")

parser.add_option("-H", "--inspiral-head",action="store_true",default=False,\
    help="Run a job using inspiral from HEAD (to get bank and cont. chisq) (experimental)")

parser.add_option("", "--disable-ifarsorting",action="store_true",\
    default=False, help="this option disables the sorting of the candidates \
    according to their IFAR value. Instead the candidates will be sorted \
    according to their combined snr, if this option is used.")

parser.add_option("", "--convert-eventid",action="store_true",\
    default=False, help="For the analysis of the second year branch triggers \
    the event_id was defined as an \"int_8s\" instead of an \"ilwd:char\" \
    Therefore when analysing S5 second year data, one must activate this \
    option. The followup code will make a copy of the original xml files in \
    the local directory and use the \"ligolw_conv_inspid\" to convert the \
    event_id into the correct type. This option also enables event_id \
    remapping when reading trigger files for plots of triggers versus time")

parser.add_option("", "--create-localcopy",action="store_true",\
    default=False, help="If this option is called, the followup pipeline will \
    make a local copy of all the useful xml files. This option is ignored \
    when \"--convert-eventid\" is requested.")

parser.add_option("", "--add-trig-columns",action="store_true",\
    default=False, help="For the analysis of teh second year branch, \
    we want to be able to use the followup code of the TRUNK. This requires to \
    handle the discrepancies between the column definitions in the \
    sngl_inspiral table. This option can be called when the \"--trig-bank\" \
    argument is also specified. When the option \"--add-trig-columns\" is \
    specified the code check for the existence of extra columns in the \
    sngl_inspiral table definition of the installed version of glue. If extra \
    columns are detected, they are added to the followup template bank xml \
    files.")

parser.add_option("-w", "--write-to-iulgroup",action="store_true", \
    default=False, help="publish the page to the iulgroup")

parser.add_option("", "--sky-map",action="store_true",default=False,\
    help="generate sky map data for an event")

parser.add_option("", "--sky-map-plot",action="store_true",default=False,\
    help="plot sky map data for an event")

parser.add_option("", "--coh-inspiral",action="store_true",default=False,\
    help="run a lalapps_inspiral job much the same way as the coherent code")

parser.add_option("-O","--odds",action="store_true",default=False,\
    help="Run coherent model selection on followup trigs (experimental)")

parser.add_option("-c", "--plot-chia",action="store_true",default=False,\
    help="Plot Coherent Inspiral results (experimental)")

parser.add_option("", "--followup-triggers",action="store_true",default=False,\
    help="Make plots of triggers versus time")

parser.add_option("", "--make-checklist",action="store_true",default=False,\
    help="Make checklist html files")

parser.add_option("", "--disable-dag-categories",action="store_true",\
    default=False,help="disable the internal dag category maxjobs")

parser.add_option("","--do-remote-ifo",action="store",\
    type="string",default=False,metavar=" IFO",\
    help="specify ifo for which calculations should be done remotely")

parser.add_option("","--activate-cm-messages",action="store_true",\
    default=False,help="when this option is called with the option " \
    "--do-remote-ifo V1 the Virgo omega scans will be run remotely at " \
    "Cascina. The communication between the LIGO cluster and Cascina use " \
    "Cm messages.")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

if opts.version:
  print("$Id$")
  sys.exit(0)

####################### SANITY CHECKS #####################################

if not opts.config_file:
  print("No configuration file specified.", file=sys.stderr)
  print("Use --config-file FILE to specify location", file=sys.stderr) 
  sys.exit(1)

if not opts.log_path and not opts.write_to_iulgroup:
  print("No log file path specified", file=sys.stderr)
  print("Use --log-path PATH to specify a location", file=sys.stderr)
  sys.exit(1)

if not opts.write_to_iulgroup and not opts.generate_fu_cache and \
  not opts.datafind and not opts.qscan and not opts.background_qscan \
  and not opts.trig_bank and not opts.inspiral and not opts.plots and not \
  opts.hoft_qscan and not opts.seis_qscan and not opts.background_hoft_qscan \
  and not opts.background_seis_qscan and not opts.hoft_datafind and not \
  opts.generate_segments and not opts.analyse_qscan and not \
  opts.analyse_seis_qscan and not opts.analyse_hoft_qscan and not \
  opts.mcmc and not opts.frame_check and not opts.inspiral_head and not \
  opts.ifo_status_check and not opts.single_qevent and not opts.H1H2_qevent \
  and not opts.plot_mcmc and not opts.inspiral_datafind and not \
  opts.analyse_hoft_qscan and not opts.distrib_remote_q and not \
  opts.coh_inspiral and not opts.sky_map and not opts.sky_map_plot and not \
  opts.convert_eventid and not opts.odds and not opts.plot_chia \
  and not opts.make_checklist and not opts.create_localcopy \
  and not opts.followup_triggers and not opts.spin_mcmc:

  print("No steps of the pipeline specified.", file=sys.stderr)
  print("Please specify at least one of", file=sys.stderr)
  print("--generate-fu-cache, --trig-bank, --inspiral, --plots,", file=sys.stderr)
  print("--datafind, --qscan, --hoft-qscan, --seis-qscan,", file=sys.stderr)
  print("--background-qscan, --background-hoft-qscan,", file=sys.stderr)
  print("--background-seis-qscan, --hoft-datafind,", file=sys.stderr)
  print("--generate-segments, --frame-check, --inspiral-datafind,", file=sys.stderr)
  print("--analyse-qscan, --analyse-seis-qscan,", file=sys.stderr)
  print("--analyse-hoft-qscan, --distrib-remote-q,", file=sys.stderr)
  print("--mcmc, --plot-mcmc, --coh-inspiral, --plot-chia,", file=sys.stderr)
  print("--sky-map, --sky-map-plot, --followup-triggers,", file=sys.stderr)
  print("--convert-eventid, --create-localcopy,", file=sys.stderr)
  print("--ifo-status-check, --single-qevent, --H1H2-qevent", file=sys.stderr)
  print("--make-checklist or --spin-mcmc or --write-to-iulgroup", file=sys.stderr) 
  sys.exit(1)

if opts.disable_followup:
  print("Warning: the option disable-followup disables any followup jobs, only qscan datafind and background qscan jobs will be run...", file=sys.stderr)

if opts.read_times:
  print("Warning: the option read-times disables the standard behaviour of the pipeline. The \"hipe-output-cache\" or \"xml-glob\" files will be ignored. Instead the times to be analysed will be read within the text files specified by the fields \"XXtimes\" of the section [triggers] of the .ini file", file=sys.stderr)

#################### READ IN THE CONFIG (.ini) FILE ########################
cp = configparser.ConfigParser()
cp.read(opts.config_file)

## set the option to make remote calculations for some Virgo qscans
if opts.do_remote_ifo:
  depIfo = opts.do_remote_ifo
  cp.set('followup-q-datafind','remote-ifo',depIfo)
  cp.set('followup-foreground-qscan','remote-ifo',depIfo)
  cp.set('followup-background-qscan','remote-ifo',depIfo)
  cp.set('followup-foreground-seismic-qscan','remote-ifo',depIfo)
  cp.set('followup-background-seismic-qscan','remote-ifo',depIfo)

## get the directory where the code is run
currentPath = os.path.abspath('.')

############# TURN THE HIPE OUTPUT INTO LAL CACHE FILES #######################

if not opts.disable_followup and not opts.read_times:
  cache = getCache(opts,cp,currentPath)

##############################################################################
# create the DAG writing the log to the specified directory
dag = followUpDAG(opts.config_file, opts.log_path)
dag.setupDAGWeb('followup web page','index.html',cp,opts)

############# READ IN THE COIRE FILES #########################################

if not opts.disable_followup and not opts.read_times:

  # if the option convert-eventid is called, we need to prepare the directory 
  # where the .xml files will be copied
  if opts.convert_eventid or opts.create_localcopy:
    if not os.access('LOCAL_XML_COPY',os.F_OK):
      os.mkdir('LOCAL_XML_COPY')
    else: pass
  print("Reading files from cache...")
  numtrigs, found, coincs, search = cache.readTriggerFiles(cp,opts)
  missed = None
  if coincs: print("found " + str(len(coincs)) + " triggers...")
  else: print("WARNING: NO COINCS FOUND...")
  if opts.trig_bank: trigbank_test = 1
  else: trigbank_test = 0
  if opts.disable_ifarsorting: ifar = False
  else: ifar = True
  trigtype = cp.get("followup-triggers", "trigger-type").strip()
  followuptrigs = getfollowuptrigs(cp,numtrigs,trigtype,dag.page,coincs,missed,search,trigbank_test,ifar,opts.add_trig_columns)
  # write out the info on the found triggers
  trigInfo = open("trigger_info.txt","w")
  for trig in followuptrigs: trig.write_trigger_info(trigInfo)
  trigInfo.close()

  print("\n.......Found " + str(len(followuptrigs)) + " trigs to follow up") 

############ SET UP THE REQUESTED JOBS ########################################

if not opts.disable_followup:
  qscanFgJob      = qscanJob(opts,cp)
  anaQscanJob     = analyseQscanJob(opts,cp)
  if opts.activate_cm_messages:
    rQscanJob       = remoteQscanJob(opts,cp)

  if not opts.read_times:
    inspDataJob     = followupDataFindJob(cp,'inspiral')
    inspJob         = followUpInspJob(cp)
    inspJobNotrig   = followUpInspJob(cp,'notrig')
    plotJob         = plotSNRCHISQJob(opts,cp)
    frcheckJob      = FrCheckJob(opts,cp)
    statusJob       = IFOstatus_checkJob(opts,cp)
    MCMCJob 	    = followupmcmcJob(opts,cp)
    PLOTMCMCJob     = plotmcmcJob(opts,cp)
    SPINMCMCJob     = followupspinmcmcJob(opts,cp)
    #PLOTSPINMCMCJob = plotspinmcmcJob(opts,cp)
    headInspJob     = followUpInspJob(cp, 'head')
    cohInspJob      = followUpInspJob(cp, 'coh')
    cohInspJobNotrig= followUpInspJob(cp, 'notrig')
    h1h2QJob        = h1h2QeventJob(opts,cp)   
    skyMapJob       = lalapps_skyMapJob(opts,cp)
    skyPlotJob      = pylal_skyPlotJob(opts,cp)
    oddsJob         = followupoddsJob(opts,cp)
    oddsPostJob     = followupOddsPostJob(opts,cp)
    chiaCDataInspJob= followUpInspJob(cp, 'chia')
    chiaJob         = followUpChiaJob(opts,cp)
    cohireJob       = followUpCohireJob(opts,cp)
    chiaPlotJob     = plotChiaJob(opts,cp)
    fuTriggerJob    = followupTriggerJob(opts,cp)
    checkListJob    = makeCheckListJob(opts,cp)

dataJob         = followupDataFindJob(cp,'futrig')
qscanBgJob      = qscanJob(opts,cp,'QSCANLITE')
distribQJob     = distributeQscanJob(cp)

print("\n.......Setting up pipeline jobs")

dq_url_pattern = "http://ldas-cit.ligo.caltech.edu/segments/S5/%s/dq_segments.txt"

segFile = {}
ifos_list = ['H1','H2','L1','G1','V1','T1']

############ PREPARE FILES FOR REMOTE (DEPORTED) VIRGO QSCANS #################
###############################################################################

if opts.do_remote_ifo and not opts.activate_cm_messages:

  depIfoIniConfig = depIfo+'config-file'
  depIfoIniWeb = depIfo+'web'
  depIfoDir = depIfo+'_qscans_config'
  depQscanList = ['foreground-qscan', 'background-qscan', 'foreground-seismic-qscan', 'background-seismic-qscan']
  
  # Build directories
  os.system('mkdir -p '+depIfoDir)
  os.system('mkdir -p '+depIfoDir+'/CONFIG')
  os.system('mkdir -p '+depIfoDir+'/RESULTS')
  for depQscan in depQscanList:
    os.system('mkdir -p '+depIfoDir+'/RESULTS/results_'+depQscan)
  os.system('mkdir -p '+depIfoDir+'/SCRIPTS')
  os.system('mkdir -p '+depIfoDir+'/TIMES')

  # Copy the qscan configuration files
  for depQscan in depQscanList:
    if cp.has_option("followup-"+depQscan, depIfoIniConfig):
      qscanConfig = string.strip(cp.get("followup-"+depQscan, depIfoIniConfig))
      if qscanConfig!='':
        print('copy '+qscanConfig+' -----> '+depIfoDir+'/CONFIG/'+depQscan+'_config.txt')
        os.system('cp '+qscanConfig+' '+depIfoDir+'/CONFIG/'+depQscan+'_config.txt')
  
  # Copy the scripts used in the remote computing center
  #   first, get the path to the scripts
  scriptPath = Popen(["which", "analyseQscan.py"], stdout=PIPE).communicate()[0]

  scriptPath = scriptPath.strip('analyseQscan.py\n')
  os.system('cp '+scriptPath+'/qsub_wscan.sh '+depIfoDir+'/SCRIPTS/')
  os.system('cp '+scriptPath+'/qsub_wscanlite.sh '+depIfoDir+'/SCRIPTS/')
  os.system('cp '+scriptPath+'/wscan_in2p3.sh '+depIfoDir+'/SCRIPTS/')
  os.system('cp '+scriptPath+'/wscanlite_in2p3.sh '+depIfoDir+'/SCRIPTS/')
  os.system('cp '+scriptPath+'/prepare_sendback.py '+depIfoDir)
  os.system('cp '+scriptPath+'/virgo_qscan_in2p3.py '+depIfoDir)
  depIfoWebForeground = string.strip(cp.get('followup-foreground-qscan', depIfoIniWeb))
  if not depIfoWebForeground.startswith('http://virgo.in2p3.fr/followups/'):
    print('\nWARNING for foreground qscans:')
    print('   wrong web address : '+depIfoWebForeground)
    print('   The web address for Virgo qscans should start with \"http://virgo.in2p3.fr/followups/\"')
    print('   followed by the name of the submitter of the jobs')
  else:
    depIfoWebForeground = depIfoWebForeground.replace('http://virgo.in2p3.fr/followups/','') + "/qscan/"
    os.system('sed -e \'s|@foreground@|'+depIfoWebForeground+'|\' '+depIfoDir+'/virgo_qscan_in2p3.py > '+depIfoDir+'/virgo_qscan_in2p3.py.tmp')
    os.system('mv -f '+depIfoDir+'/virgo_qscan_in2p3.py.tmp '+depIfoDir+'/virgo_qscan_in2p3.py')
    os.system('chmod +x '+depIfoDir+'/virgo_qscan_in2p3.py')
  depIfoWebForegroundSeismic = string.strip(cp.get('followup-foreground-seismic-qscan', depIfoIniWeb))
  if not depIfoWebForegroundSeismic.startswith('http://virgo.in2p3.fr/followups/'):
    print('\nWARNING for foreground-seismic qscans:')
    print('   wrong web address : '+depIfoWebForegroundSeismic)
    print('   The web address for Virgo qscans should start with \"http://virgo.in2p3.fr/followups/\"')
    print('   followed by the name of the submitter of the jobs')
  else:
    depIfoWebForegroundSeismic = depIfoWebForegroundSeismic.replace('http://virgo.in2p3.fr/followups/','') + "/seismic_qscan/"
    os.system('sed -e \'s|@foreground-seismic@|'+depIfoWebForegroundSeismic+'|\' '+depIfoDir+'/virgo_qscan_in2p3.py > '+depIfoDir+'/virgo_qscan_in2p3.py.tmp')
    os.system('mv -f '+depIfoDir+'/virgo_qscan_in2p3.py.tmp '+depIfoDir+'/virgo_qscan_in2p3.py')
    os.system('chmod +x '+depIfoDir+'/virgo_qscan_in2p3.py')

############# LOOP OVER RANDOM BACKGROUND TIMES ###############################
###############################################################################
# Prepare the qscan background
if opts.background_qscan or opts.background_hoft_qscan or opts.background_seis_qscan:
  for ifo in ifos_list:
    times, timeListFile = getQscanBackgroundTimes(cp,opts,ifo,segFile)
    for time in times:
      # SETUP DATAFIND JOBS FOR BACKGROUND QSCANS (REGULAR DATA SET)
      dNode = followupDataFindNode(dataJob,'futrig','q-datafind',cp,time,ifo,opts,dag,'datafind')

      # SETUP DATAFIND JOBS FOR BACKGROUND QSCANS (HOFT)
      dHoftNode = followupDataFindNode(dataJob,'futrig','q-hoft-datafind',cp,time,ifo,opts,dag,'hoft_datafind')
    
      # SETUP BACKGROUND QSCAN JOBS
      qBgNode = qscanNode(qscanBgJob,time,cp,dHoftNode.outputFileName,ifo,'background-hoft-qscan',opts,dHoftNode,dag,'hoft_datafind','background_hoft_qscan')
      qBgNode = qscanNode(qscanBgJob,time,cp,dNode.outputFileName,ifo,'background-qscan',opts, dNode, dag,'datafind','background_qscan')
      qBgNode = qscanNode(qscanBgJob,time,cp,dNode.outputFileName,ifo,'background-seismic-qscan',opts, dNode, dag,'datafind','background_seis_qscan')
  
    # WRITE TIMES FOR REMOTE (DEPORTED) CALCULATIONS
    if opts.do_remote_ifo and not opts.activate_cm_messages:
      if ifo == depIfo and len(times)!=0:
        os.system('cp '+timeListFile+' '+depIfoDir+'/TIMES/background_qscan_times.txt')

############# Preparing deported qscans and distribution node #################
###############################################################################

if opts.do_remote_ifo and not opts.activate_cm_messages:
  # Initialize qscan times text file for deported qscans
  if opts.qscan or opts.seis_qscan:
    depQscanTFile = open(depIfoDir+'/TIMES/qscan_times.txt','w')

  # Distribute the results of the deported qscans to the expected output directories
  distribQNode = distributeQscanNode(distribQJob,'QSCAN/QSCAN.cache','QSCANLITE/QSCANLITE.cache',depIfo,depIfo+"_qscans_results.tar.gz",opts,dag)

########## LOOP OVER FOLLOWUP TIMES (if read-times is activated) ##############
###############################################################################

if opts.read_times:
  for ifo in ifos_list:
    dag.appendSection("Triggers in " + ifo)
    times,timeListFile = getForegroundTimes(cp,opts,ifo)
    for time in times:

      #Lets make a new section for this event
      event_title = ifo + ' @ GPS ' + repr(time)
      dag.appendSubSection(event_title)

      # Prepare deported qscans
      if opts.do_remote_ifo and ifo==depIfo:
         if opts.activate_cm_messages:
           rQFgNode = remoteQscanFgNode(rQscanJob,time,cp,ifo,'foreground-qscan',opts,dag,'qscan')
           rQFgNode = remoteQscanFgNode(rQscanJob,time,cp,ifo,'foreground-seismic-qscan',opts,dag,'seis_qscan')
         else:
           if opts.qscan or opts.seis_qscan:
             depQscanTFile.write(repr(time)+'\n')

      # SETUP DATAFIND JOBS FOR HOFT QCANS
      dHoftNode = followupDataFindNode(dataJob,'futrig','q-hoft-datafind',cp,time,ifo,opts,dag,'hoft_datafind')

      # SETUP DATAFIND JOBS FOR REGULAR QSCANS
      dNode = followupDataFindNode(dataJob,'futrig','q-datafind',cp,time,ifo,opts,dag,'datafind')

      # SETUP FOREGROUND QSCAN JOBS
      qFgNode = qscanNode(qscanFgJob,time,cp,dHoftNode.outputFileName,ifo,'foreground-hoft-qscan',opts,dHoftNode,dag,'hoft_datafind','hoft_qscan')
      qFgNode = qscanNode(qscanFgJob,time,cp,dNode.outputFileName,ifo,'foreground-qscan',opts,dNode,dag,'datafind','qscan')
      qFgNode = qscanNode(qscanFgJob,time,cp,dNode.outputFileName,ifo,'foreground-seismic-qscan',opts,dNode,dag,'datafind','seis_qscan')

    # Loop for analyseQscan jobs (which must be children of all qscan jobs)
    for time in times:

      anaQscanNode = analyseQscanNode(anaQscanJob,time,ifo,'foreground-qscan','QSCAN/QSCAN.cache','QSCANLITE/QSCANLITE.cache',cp,opts,dag,'analyse_qscan')
      anaQscanNode = analyseQscanNode(anaQscanJob,time,ifo,'foreground-seismic-qscan','QSCAN/QSCAN.cache','QSCANLITE/QSCANLITE.cache',cp,opts,dag,'analyse_seis_qscan')
      anaQscanNode = analyseQscanNode(anaQscanJob,time,ifo,'foreground-hoft-qscan','QSCAN/QSCAN.cache','QSCANLITE/QSCANLITE.cache',cp,opts,dag,'analyse_hoft_qscan')

################# LOOP OVER FOLLOWUP TRIGGERS #################################
###############################################################################

if not opts.disable_followup and not opts.read_times:
  
  # LOOP ON TRIGGERS (=CANDIDATES)
  for trig in followuptrigs:

    #Write table of trigger parameters
    writeParamTable(trig,opts)
 
    # Initialization of local variables
    dHoftNode = {}

    # We need a source localization and a chia node at the 
    # coincident trigger level. it takes many ifos as input. 
    skyNode = lalapps_skyMapNode(skyMapJob,trig,opts)
    chiaNode = followUpChiaNode(chiaJob,trig,opts,dag,cp)
    cohireNode = followUpCohireNode(cohireJob,chiaJob.outputPath,trig,chiaNode,dag,dag.page,opts)
    chiaPlotNode = plotChiaNode(chiaPlotJob,chiaJob.outputPath,trig,cohireNode,dag,dag.page,opts,cp)
 
    #Lets make a new section for this event
    event_title = 'Trigger ID = '+str(trig.eventID)+ ' Stat = '+str(trig.statValue)
    if not opts.disable_ifarsorting:
      event_title += ' FAR = ' + str(trig.far)
    dag.appendSection(event_title)

    # PLOT TRIGGERS VERSUS TIME
    fuTriggerNode = followupTriggerNode(fuTriggerJob,trig,cp,opts,dag)
 
    # TRY GETTING INSPIRAL PROCESS PARAMS... 
    # currently done even when inspiral jobs are not requested... fix?
    inspiral_process_params = cache.processFollowupCache(cp, opts, trig)
    bank_process_params, bank = cache.processFollowupCache(cp, opts, trig, 'TMPLTBANK_*') 

    # loop over ifos found in coincidence
    for ifo in trig.ifolist_in_coinc:
      try: trig.gpsTime[ifo]
      except: continue

      # Write time in file for deported qscans
      if opts.do_remote_ifo and ifo==depIfo:
         if opts.activate_cm_messages:
           rQFgNode = remoteQscanFgNode(rQscanJob,trig.gpsTime[ifo],cp,ifo,'foreground-qscan',opts,dag,'qscan')
           rQFgNode = remoteQscanFgNode(rQscanJob,trig.gpsTime[ifo],cp,ifo,'foreground-seismic-qscan',opts,dag,'seis_qscan')
         else:
           if opts.qscan or opts.seis_qscan:
             depQscanTFile.write(repr(trig.gpsTime[ifo])+'\n')
 
      #Lets append a new subsection to the last section
      dag.appendSubSection(str(ifo)+' @ GPS '+str(trig.gpsTime[ifo]))

      # SETUP DATAFIND JOBS FOR HOFT QCANS
      dHoftNode[ifo] = followupDataFindNode(dataJob,'futrig','q-hoft-datafind',cp,trig.gpsTime[ifo],ifo,opts,dag,'hoft_datafind')

      # SETUP DATAFIND JOBS FOR REGULAR QSCANS
      dNode = followupDataFindNode(dataJob,'futrig','q-datafind',cp,trig.gpsTime[ifo],ifo,opts,dag,'datafind')
    
      # SETUP FOREGROUND QSCAN JOBS
      qFgNode = qscanNode(qscanFgJob,trig.gpsTime[ifo],cp,dHoftNode[ifo].outputFileName,ifo,'foreground-hoft-qscan',opts,dHoftNode[ifo],dag,'hoft_datafind','hoft_qscan')
      qFgNode = qscanNode(qscanFgJob,trig.gpsTime[ifo],cp,dNode.outputFileName,ifo,'foreground-qscan',opts,dNode,dag,'datafind','qscan')
      qFgNode = qscanNode(qscanFgJob,trig.gpsTime[ifo],cp,dNode.outputFileName,ifo,'foreground-seismic-qscan',opts,dNode,dag,'datafind','seis_qscan')

      # SETUP DATAFIND JOBS FOR INSPIRAL, FRAME CHECKS, AND MCMC JOBS (REQUIRED 
      # ONLY IF THE INSPIRAL_HIPE CACHE FILES ARE OBSOLETE) 
      inspiralDataFindNode = followupDataFindNode(inspDataJob,'inspiral','inspiral-datafind',cp,trig.gpsTime[ifo],ifo,opts,dag,'inspiral_datafind',inspiral_process_params[ifo])

      # SETUP INSPIRAL JOBS
      inspiralNode = followUpInspNode(inspJob,inspiral_process_params[ifo], ifo, trig, cp, opts, dag, inspiralDataFindNode.outputFileName,inspiralDataFindNode,'inspiral_datafind')
      cohInspNode = followUpInspNode(cohInspJob,inspiral_process_params[ifo], ifo, trig, cp, opts, dag, inspiralDataFindNode.outputFileName,inspiralDataFindNode,'inspiral_datafind',type='coh')
      chiaCDataNode = followUpInspNode(chiaCDataInspJob,inspiral_process_params[ifo], ifo, trig, cp, opts, dag, inspiralDataFindNode.outputFileName,inspiralDataFindNode,'inspiral_datafind',type='chia')
      #skyNode.append_insp_node(inspiralNode,ifo)
      skyNode.append_insp_node(cohInspNode,ifo)
      chiaNode.append_insp_node(chiaCDataNode,ifo)     
      chiaPlotNode.append_insp_node(chiaCDataNode,ifo) 
      headInspiralNode = followUpInspNode(headInspJob,inspiral_process_params[ifo], ifo, trig, cp, opts, dag, inspiralDataFindNode.outputFileName,inspiralDataFindNode,'inspiral_datafind', 'head', bank)


      # SETUP PLOT JOBS
      plotNode = plotSNRCHISQNode(plotJob,ifo,inspiralNode.output_file_name,trig,dag.page,dag,inspiralNode,opts)

      # SETUP FRAME CHECK JOBS
      frcheckNode = FrCheckNode(frcheckJob,inspiral_process_params[ifo], ifo, trig, cp, opts, dag, inspiralDataFindNode.outputFileName, inspiralDataFindNode, 'inspiral_datafind')

      # MCMC JOBS
      # run multiple MCMC chains for each trigger (1 chain = 1 node)
      chainNumber = string.strip(cp.get('followup-mcmc','chain_nb'))
      mcmcIdList = []
      for k in range(int(chainNumber)):
        randomseed = str(trig.gpsTime[ifo] + k).split('.')[0][5:9]
        MCMCNode = followupmcmcNode(MCMCJob,inspiral_process_params, trig, randomseed, cp, opts, dag, ifo)
        mcmcIdList.append(MCMCNode.id)

      # PLOT MCMC JOBS
      PLOTMCMCNode = plotmcmcNode(PLOTMCMCJob,trig,mcmcIdList,cp,opts,dag,MCMCNode.ifoRef,MCMCNode.ifonames)

      # SETUP STATUS SUMMARY PLOT JOBS
      statusNode = IFOstatus_checkNode(statusJob, ifo, trig, cp, opts, dag)


    # SETUP H1H2 SPECIFIC JOBS
    # first make sure that the coincident trigger is found in at least one Hanford ifo and that both Hanford ifos are in the analysed times
    lho_flag = trig.ifolist_in_coinc.count('H1') + trig.ifolist_in_coinc.count('H2')
    if lho_flag == 2:
      #Lets append a new subsection to the last section
      dag.appendSubSection('H1H2 followup' +' @ GPS '+str(trig.gpsTime['H1']))
      h1h2Qevent = h1h2QeventNode(h1h2QJob,dHoftNode,trig.gpsTime,['H1','H2'],'qevent',cp,opts,dag,'H1H2_qevent')


    # LOOP OVER IFOS NOT FOUND IN COINCIDENCE (THOUGH IN THE ANALYSED TIMES)
    for j in range(0,len(trig.ifoTag)-1,2):
      ifo = trig.ifoTag[j:j+2]
      if not trig.ifolist_in_coinc.count(ifo):

        #Lets append a new subsection to the last section
        dag.appendSubSection(str(ifo)+" (non coincident ifo)")

        # Write time in file for deported qscans
        if opts.do_remote_ifo and ifo==depIfo and not opts.activate_cm_messages:
           if opts.qscan or opts.seis_qscan:
             depQscanTFile.write(repr(trig.gpsTime[trig.ifolist_in_coinc[0]])+'\n')

        # SETUP DATAFIND JOBS FOR HOFT QCANS
        dHoftNode = followupDataFindNode(dataJob,'futrig','q-hoft-datafind',cp,trig.gpsTime[trig.ifolist_in_coinc[0]],ifo,opts,dag,'hoft_datafind')

        # SETUP DATAFIND JOBS FOR REGULAR QCANS
        dNode = followupDataFindNode(dataJob,'futrig','q-datafind',cp,trig.gpsTime[trig.ifolist_in_coinc[0]],ifo,opts,dag,'datafind')

        # SETUP FOREGROUND QSCAN JOBS
        qFgNode = qscanNode(qscanFgJob,trig.gpsTime[trig.ifolist_in_coinc[0]],cp,dHoftNode.outputFileName,ifo,'foreground-hoft-qscan',opts,dHoftNode,dag,'hoft_datafind','hoft_qscan')
        qFgNode = qscanNode(qscanFgJob,trig.gpsTime[trig.ifolist_in_coinc[0]],cp,dNode.outputFileName,ifo,'foreground-qscan',opts,dNode,dag,'datafind','qscan')

        # SETUP DATAFIND JOBS FOR INSPIRAL JOBS 
        # (REQUIRED ONLY IF THE INSPIRAL_HIPE CACHE FILES ARE OBSOLETE)
        inspiralDataFindNode = followupDataFindNode(inspDataJob,'inspiral','inspiral-datafind',cp,trig.gpsTime[trig.ifolist_in_coinc[0]],ifo,opts,dag,'inspiral_datafind',inspiral_process_params[ifo])

        cohInspNodeNotrig = followUpInspNode(cohInspJob,inspiral_process_params[ifo], trig.ifolist_in_coinc[0], trig, cp, opts, dag, inspiralDataFindNode.outputFileName,inspiralDataFindNode,'inspiral_datafind',type='coh')
        chiaCDataNodeNotrig = followUpInspNode(chiaCDataInspJob,inspiral_process_params[ifo], trig.ifolist_in_coinc[0], trig, cp, opts, dag, inspiralDataFindNode.outputFileName,inspiralDataFindNode,'inspiral_datafind',type='chia')
        chiaNode.append_insp_node(chiaCDataNodeNotrig,ifo)
        chiaPlotNode.append_insp_node(chiaCDataNodeNotrig,ifo)
        skyNode.append_insp_node(cohInspNodeNotrig,ifo)

        # generate the snr/chisq time series for each triggered template
        for itf in trig.ifolist_in_coinc:
          # SETUP INSPIRAL JOBS
          inspiralNode = followUpInspNode(inspJobNotrig,inspiral_process_params[ifo], itf, trig, cp, opts, dag, inspiralDataFindNode.outputFileName,inspiralDataFindNode,'inspiral_datafind','notrig',None)
          # SETUP PLOT JOBS
          plotNode = plotSNRCHISQNode(plotJob,ifo,inspiralNode.output_file_name,trig,dag.page,dag,inspiralNode,opts,itf)

    # SOURCE LOCALIZATION.  This requires the complex frame data from the above
    # inspiral jobs for all ifos that are on.  This should be replaced with the
    # output of the coherent inspiral stages which use the same template
    # Set up a sky map node.  These take inputs from all the triggers and
    # must appear outside the loop over ifos.  But we need the output file 
    # names from several inspiral jobs since they are only run over one ifo
    # we could walk the graph and find this (maybe) but instead I'll just track
    # the names we need. 
    skyNode.add_node_to_dag(dag,opts,trig)
    skyPlotNode = pylal_skyPlotNode(skyPlotJob,trig,skyNode,dag,dag.page,opts)

    # MCMC JOBS
    # run multiple MCMC chains for each trigger (1 chain = 1 node)
    chainNumber = string.strip(cp.get('followup-mcmc','chain_nb'))
    mcmcIdList = []
    for k in range(int(chainNumber)):
      randomseed = str(trig.gpsTime[trig.ifolist_in_coinc[0]] + k).split('.')[0][5:9]
      MCMCNode = followupmcmcNode(MCMCJob,inspiral_process_params, trig, randomseed, cp, opts, dag)
      mcmcIdList.append(MCMCNode.id)

    # PLOT MCMC JOBS
    PLOTMCMCNode = plotmcmcNode(PLOTMCMCJob,trig,mcmcIdList,cp,opts,dag,MCMCNode.ifoRef,MCMCNode.ifonames)

    # SPIN MCMC JOBS
    chainNumber = string.strip(cp.get('followup-spin-mcmc','chain_nb'))
    spinmcmcIdList = []
    for k in range(int(chainNumber)):
      SPINMCMCNode = followupspinmcmcNode(SPINMCMCJob,inspiral_process_params, trig, cp, opts, dag, k)
      spinmcmcIdList.append(SPINMCMCNode.id)

    # PLOT SPIN MCMC JOBS
    # PLOTSPINMCMCNode = plotspinmcmcNode(PLOTSPINMCMCJob,trig,spinmcmcIdList,cp,opts,dag)

    #oddsNode = followupoddsNode(oddsJob,inspiral_process_params, trig, randomseed, cp, opts, dag)

    # ODDS RATIO. This calculates the odds ratio vs Gaussian noise and outputs a nested sampling chain
    # Will split across CPUs
    oddsIdList=[]
    Nlive_tot=int(string.strip(cp.get('followup-odds','live-points')))
    Nlive_base=int(string.strip(cp.get('followup-odds','min-live')))
    Nparallel=int(math.ceil(Nlive_tot/Nlive_base))
    for k in range(Nparallel):
        randomseed = [str(trig.gpsTime[trig.ifolist_in_coinc[0]] + k).split('.')[0][5:9],str(trig.gpsTime[trig.ifolist_in_coinc[0]]).split('.')[1]]
        oddsNode = followupoddsNode(oddsJob,inspiral_process_params, trig, randomseed, cp, opts, dag)
        oddsIdList.append(oddsNode.id)

    # plot odds output
    oddsPostNode = followupOddsPostNode(oddsPostJob,inspiral_process_params,trig,oddsIdList, cp, opts, dag)

    #SETUP CHIA JOBS
    #chiaPlotNode = plotChiaNode(chiaPlotJob,chiaJob.outputPath,trig,chiaNode,dag,dag.page,opts)
    #cohireNode = followUpCohireNode(cohireJob,chiaJob.outputPath,trig,chiaNode,dag,dag.page,opts)
    #chiaPlotNode = plotChiaNode(chiaPlotJob,chiaJob.outputPath,trig,cohireNode,dag,dag.page,opts,cp)

  # Loop for analyseQscan jobs (which must be children of all qscan jobs)
  for trig in followuptrigs:
    for ifo in trig.ifolist_in_coinc:
      anaQscanNode = analyseQscanNode(anaQscanJob,trig.gpsTime[ifo],ifo,'foreground-qscan','QSCAN/QSCAN.cache','QSCANLITE/QSCANLITE.cache',cp,opts,dag,'analyse_qscan')
      anaQscanNode = analyseQscanNode(anaQscanJob,trig.gpsTime[ifo],ifo,'foreground-seismic-qscan','QSCAN/QSCAN.cache','QSCANLITE/QSCANLITE.cache',cp,opts,dag,'analyse_seis_qscan')
      anaQscanNode = analyseQscanNode(anaQscanJob,trig.gpsTime[ifo],ifo,'foreground-hoft-qscan','QSCAN/QSCAN.cache','QSCANLITE/QSCANLITE.cache',cp,opts,dag,'analyse_hoft_qscan')


  # Prepare the makeCheckList jobs (must be children of all jobs)
  for trig in followuptrigs:
    checkListNode = makeCheckListNode(checkListJob,trig,cp,opts,dag) 

# Prepare for moving the deported qscan directory (tar-gzip)
if opts.do_remote_ifo and not opts.activate_cm_messages:
  if opts.qscan or opts.seis_qscan:
    depQscanTFile.close()
  os.system('tar zcvf '+depIfoDir+'.tar.gz '+depIfoDir)
    
##########################################################################
dag.writeAll()
if opts.write_to_iulgroup:
  dag.publishToHydra()

# write out a log file for this script
log_fh = open(opts.config_file.replace('.ini','.log'), 'w')

log_fh.write("lalapps_followup_pipe")
for arg in command_line:
  if arg[0] == '-':
    log_fh.write( "\n" )
  log_fh.write( arg + ' ')



sys.exit(0)

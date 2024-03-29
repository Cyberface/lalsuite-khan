"""
ihope.in - weekly automate pipeline driver script

This script generates the condor DAG necessary to analyze a given amount
of LIGO and GEO data.  It takes in the start and end times, generates the
segment lists, runs the zero-lag, time shift and injection analyses, generates
summary information and plots and follows up the loudest events.
"""

from __future__ import print_function

__author__ = 'Stephen Fairhurst <sfairhur@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

##############################################################################
# import standard modules
import os, sys, copy, shutil, re
import optparse
import tempfile
import urllib
import time
import calendar
import shutil
import M2Crypto

##############################################################################
# import the modules we need to build the pipeline
from ligo import segments
from ligo.segments import utils as segmentsUtils
from glue import pipeline
from glue.pipeline import DeepCopyableConfigParser as dcConfigParser
from lalapps import inspiralutils


##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################
usage = """usage: %prog [options] 

lalapps_ihope is designed to run the inspiral analysis end to end.  It
performs several distinct steps.  The default is for the full search to
be run, various steps can be skipped by specifying --skip options.  The
required arguments are:

--config-file 
--log path
--gps-start-time
--gps-end-time

The options to run the analysis are specified in the config-file
(normally ihope.ini) file.

1) Generate the segments for which the analysis should be done.  This
can be skipped using --skip-generate-segments.  The option
--use-available-data will restrict the search to those times for which
there is data on the cluster the analysis is being run.

2) Generate the veto segments.  Can be skipped with
--skip-generate-veto-segments.

3) Find the data on the cluster.  Can be skipped with --skip-datafind.

4) Generate the template banks.  Can be skipped with --skip-tmpltbank.

5) The search itself.  Can be skipped with --skip-search.

6) The application of the data quality vetoes.  Can be skipped with
--skip-data-quality.

7) Make plots of the output.  Can be skipped with --skip-plots.

8) Followup the loudest events (currently not working).  Can be skipped
with --skip-followup.

For options 5-8, you can run 3 types of analysis:
a) playground and playground time slides. Skipped with --skip-playground
b) full data and full data time slides.  Skipped with --skip-full-data
c) software injections. Skipped with --skip-injections
"""

parser = optparse.OptionParser( usage=usage, \
      version= "%prog CVS\n" +
      "$Id$\n" +
      "$Name$\n")

# arguments
parser.add_option("-f", "--config-file",action="store",type="string",\
    metavar=" FILE", help="use configuration file FILE")

parser.add_option("-p", "--log-path",action="store",type="string",\
    metavar=" PATH", \
    help="directory to write condor log file, should be a local directory")

parser.add_option("-l", "--node-local-dir",action="store",type="string",\
    metavar=" NODEPATH", \
    help='''directory to write the temporary sql databases. This must be a
directory that is local to all nodes.''')

parser.add_option("-o", "--node-access-path",action="store",type="string",\
    metavar=" OUTPATH", \
    help="node path for outputting the temporary coherent files to.")

parser.add_option("-u", "--username",action="store",type="string",\
    metavar=" USERNAME", help='''username on cluster for constructing the
full path-name for outputting the temporary coherent files.''')

parser.add_option("-s", "--gps-start-time",action="store",type="int",\
    metavar=" GPS_START", help="begin analysis at GPS_START")

parser.add_option("-e", "--gps-end-time",action="store",type="int",\
    metavar=" GPS_END", help="end analysis at GPS_END")

parser.add_option("-D", "--use-available-data",action="store_true",\
    default=False, help="analyze only the data available on the cluster")

parser.add_option("-R", "--reverse-analysis",action="store_true",\
    default=False, help="do the entrire analysis using reverse chirp method")

parser.add_option("-S", "--skip-generate-segments",action="store_false",\
    default=True, dest="generate_segments", \
    help="skip generation segments for analysis")

parser.add_option("-V", "--skip-generate-veto-segments",action="store_false",\
    default=True, dest="generate_veto_segments", \
    help="skip generation segments for analysis")

parser.add_option("-B", "--skip-datafind",action="store_false",\
    default=True, dest="run_datafind", help="skip the datafind step")

parser.add_option("-T", "--skip-tmpltbank",action="store_false",\
    default=True, dest="run_tmpltbank", help="skip the template bank generation")

parser.add_option("-F", "--skip-full-data",action="store_false",\
    default=True, dest="run_full_data", help="skip the full data + slides")

parser.add_option("-P", "--skip-playground",action="store_false",\
    default=True, dest="run_playground", help="skip the playground analysis")

parser.add_option("-I", "--skip-injections",action="store_false",\
    default=True, dest="run_injections", \
    help="skip the inspiral analysis with software injections")

parser.add_option("-A", "--skip-search",action="store_false",\
    default=True, dest="run_search",
    help="skip the search of the data (i.e. don't run inspiral hipe)")

parser.add_option("-Q", "--skip-data-quality",action="store_false",\
    default=True, dest="run_data_quality", \
    help="skip generation dq veto segments and use in analysis ")

parser.add_option("-Z", "--skip-plots",action="store_false",\
    default=True, dest="run_plots",  help="skip the plotting step")

parser.add_option("-U", "--skip-followup",action="store_false",\
    default=True, dest="run_followup",  help="skip the inspiral followup")

parser.add_option("-H", "--skip-pipedown",action="store_false",\
    default=True, dest="run_pipedown",  help="skip running pipedown")

parser.add_option("-W","--skip-hardware-injections",action="store_false",\
    default=True,dest="run_hardware_inj", help="Skip the hardware injection script")

parser.add_option("-w","--skip-omega-scans",action="store_false",\
    default=True,dest="do_omega_scans", help="Skip the omega scan setup.")

parser.add_option("-M","--run-mvsc",action="store_true",\
    default=False,dest="run_mvsc", help="Run the multivariate statisical classifier in pipedown")

parser.add_option("-x", "--dax",action="store_true",\
    default=False, help="Delete the ligo_data_find jobs and populate frame LFNs in the DAX")

parser.add_option("-r", "--reuse-data",action="store_true",\
    default=False, help="In DAX mode, reuse data from a previous run stored in the static PFN cache")

parser.add_option("--data-checkpoint", action="store_true",default=False,\
    help="checkpoint the inspiral code")

parser.add_option("--use-gpus", action="store_true", default=False,\
    help="run inspiral jobs on GPU nodes")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

##############################################################################
# Sanity check of input arguments
if not opts.config_file:
  print("No configuration file specified.", file=sys.stderr)
  print("Use --config-file FILE to specify location.", file=sys.stderr)
  sys.exit(1)

if not opts.config_file[-4:]==".ini":
  print("Configuration file name must end in '.ini'!", file=sys.stderr) 
  sys.exit(1)

if not opts.log_path:
  print("No log file path specified.", file=sys.stderr)
  print("Use --log-path PATH to specify a location.", file=sys.stderr)
  sys.exit(1)

if opts.run_pipedown:
  if not opts.node_local_dir:
    print("No local dir specified. If running with pipedown", file=sys.stderr)
    print("use --node-local-dir to specify a local directory", file=sys.stderr)
    print("Use --help for more information", file=sys.stderr)
    sys.exit(1)

if not opts.gps_start_time:
  print("No GPS start time specified for the analysis", file=sys.stderr)
  print("Use --gps-start-time GPS_START to specify a location.", file=sys.stderr)
  sys.exit(1)

if not opts.gps_end_time:
  print("No GPS end time specified for the analysis", file=sys.stderr)
  print("Use --gps-end-time GPS_END to specify a location.", file=sys.stderr)
  sys.exit(1)

if opts.gps_end_time < opts.gps_start_time:
  print("The GPS end time must be after the GPS start time", file=sys.stderr)
  sys.exit(1)

opts.complete_cache = (opts.run_data_quality or opts.run_plots or opts.run_pipedown or opts.run_search)

##############################################################################
# check the user has a valid grid proxy

def check_grid_proxy(path):
  try:
    proxy = M2Crypto.X509.load_cert(path)
  except Exception as e:
    msg = "Unable to load proxy from path %s : %s" % (path, e)
    raise RuntimeError(msg)

  try:
    proxy.get_ext("proxyCertInfo")
  except LookupError:
    subject = proxy.get_subject().as_text()
    if re.search(r'.+CN=proxy$', subject):
      msg = "Proxy %s is not RFC compliant" % path
      raise RuntimeError(msg)

  try:
    expireASN1 = proxy.get_not_after().__str__()
    expireGMT  = time.strptime(expireASN1, "%b %d %H:%M:%S %Y %Z")
    expireUTC  = calendar.timegm(expireGMT)
    now = int(time.time())
  except Exception as e:
    msg = "could not determine time left on proxy: %s" % e
    raise RuntimeError(msg)

  return expireUTC - now

time_needed = 3600 * 6
tmp_proxy_time_left = proxy_time_left = 0
my_uid = os.getuid()
proxy_path = "/tmp/x509up_u%d" % my_uid 

try:
  # check to see how much time is left on the environment's proxy
  if os.environ.has_key('X509_USER_PROXY'):
    tmp_proxy_path = os.environ['X509_USER_PROXY']
    if os.access(tmp_proxy_path, os.R_OK):
      try:
        tmp_proxy_time_left = check_grid_proxy(tmp_proxy_path)
      except:
        tmp_proxy_time_left = 0

  # check to see how much time there is left on the proxy in /tmp
  if os.access(proxy_path, os.R_OK):
    try:
      proxy_time_left = check_grid_proxy(proxy_path)
    except:
      proxy_time_left = 0

  # use the one with the most time left
  if tmp_proxy_time_left > 0 and tmp_proxy_time_left > proxy_time_left:
    shutil.copy(tmp_proxy_path,proxy_path)

  # check that the proxy is valid and that enough time remains
  time_left = check_grid_proxy(proxy_path)
  if time_left < 0:
    raise RuntimeError("Proxy has expired.")
  elif time_left < time_needed:
    msg = "Not enough time left on grid proxy (%d seconds)." % time_left
    raise RuntimeError(msg)
  else:
    os.environ['X509_USER_PROXY'] = proxy_path

except:
  print("""
Error: Could not find a valid grid proxy. Please run 

   ligo-proxy-init albert.einstein

with your LIGO.ORG username to generate a new grid proxy certificate.

You may need to log out and do this on your home machine if your .globus
directory is not present on the machine where you are running ihope.

You can check the status of your grid proxy by running

  grid-proxy-info

At least %d seconds must be left before your proxy expires to run ihope.
""" % time_needed, file=sys.stderr)
  raise

##############################################################################
# set up the analysis directory
analysisDirectory = str(opts.gps_start_time) + "-" + str(opts.gps_end_time)
inspiralutils.mkdir(analysisDirectory)

# copy the ini file into the directory
shutil.copy( opts.config_file, analysisDirectory )
opts.config_file = opts.config_file.split("/")[-1]

os.chdir(analysisDirectory)
inspiralutils.mkdir("logs")

# parse the ini file:
cp = dcConfigParser()
cp.read(opts.config_file)

# set gps start and end times in the cp object.
cp.set("input", "gps-start-time", str(opts.gps_start_time) )
cp.set("input", "gps-end-time", str(opts.gps_end_time) )
cp.set("input", "output-path", str(opts.node_access_path) )
cp.set("input", "username", str(opts.username) )

# Check for --disable-dag-categories
opts.do_dag_categories = not cp.has_option(
              "hipe-arguments","disable-dag-categories")

# Tell matplotlib to put its cache files in the the local dir
cp.set("pipeline", "matplotlibdir", opts.node_local_dir )

# Add local dir to pipeline
if opts.run_pipedown:
  cp.set("pipeline", "node-tmp-dir", opts.node_local_dir)

# Add GPU flag
if opts.use_gpus:
  cp.set("condor", "use-gpus", "")

##############################################################################
# create a directory called executables and copy them over
inspiralutils.mkdir("executables")
os.chdir("../")
for (job, executable) in cp.items("condor"):
  if job not in ("universe", "use-gpus") and executable != "":
    shutil.copy( executable, analysisDirectory + "/executables" )
    if job == "hardware_inj_page":
      executable = "executables/" + executable.split("/")[-1]
    else:
      executable = "../executables/" + executable.split("/")[-1]
    cp.set("condor", job, executable)

os.chdir(analysisDirectory)

##############################################################################
# create a log file that the Condor jobs will write to
# 
# we already checked that the last 4 characters of the config file are '.ini'
basename = opts.config_file[:-4]
logname = basename + '.dag.log.'
try:
  os.mkdir(opts.log_path)
except:
  pass
tempfile.tempdir = opts.log_path
logfile = tempfile.mktemp(prefix=logname)
fh = open( logfile, "w" )
fh.close()

##############################################################################
# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile,dax=opts.dax)
dag.set_dag_file(basename)
dag.set_dax_file(basename)

##############################################################################
# Set up the IFOs and get the appropriate segments
ifos = [] 

for option in ["g1-data","h1-data","h2-data","l1-data","v1-data"]:
  if cp.has_option("ifo-details",option): ifos.append(option[0:2].upper() )

if cp.has_option("ifo-details","analyze-all"): 
  print("The inspiral pipeline does not yet support coincidence between", file=sys.stderr)
  print("all five IFOs. Do not use the analyze-all option.", file=sys.stderr)
  sys.exit(1)

ifos.sort()

print("Setting up an analysis for " + str(ifos) + " from " + \
    str(opts.gps_start_time) + " to "  + str(opts.gps_end_time))
print()
sys.stdout.flush()


##############################################################################
# Create a static frame cache file for use by pegasus

ifo_frame_type_dict = {}
for ifo in ifos:
  ifo_frame_type_dict[ifo] = inspiralutils.get_data_options(cp,ifo)[1]

# determine which data server port to use
try:
  ldrserver = cp.get("datafind","server")
except:
  ldrserver = "internal_ldr_port"

pfnfile = inspiralutils.create_frame_pfn_file(ifo_frame_type_dict,opts.gps_start_time,opts.gps_end_time,server=ldrserver)
peg_frame_cache = os.path.join(os.getcwd(),inspiralutils.create_pegasus_cache_file(pfnfile))


##############################################################################
# create a local execute directory for the pipeline
tmp_exec_dir = tempfile.mkdtemp( prefix="%s-%s-%s-%d." % (''.join(ifos),
  basename, opts.gps_start_time, 
  int(opts.gps_end_time) - int(opts.gps_start_time)) )
os.chmod(tmp_exec_dir, 0o755)

##############################################################################
# Determine the segments to analyze

# first, copy the segment files and make sure that the ini file reflects the
# new path, relative to the DAG directories.

inspiralutils.mkdir("segments")
os.chdir("..")
for vetoes, infile in cp.items("segments"):
  inspiralutils.copyCategoryFiles(cp,vetoes,"segments",infile,\
                analysisDirectory)

os.chdir(analysisDirectory)
os.chdir("segments")

# Determine which veto categories to filter
veto_categories = [int(cat) for cat in \
    cp.get("segments", "veto-categories").split(",")]
if 1 not in veto_categories: 
  raise ValueError("The list of veto categories must include 1!")
veto_categories.sort()
# separate out non-CAT1 categories for some purposes
vetocats_without1 = [cat for cat in veto_categories if cat!=1]

# Download the veto definer xml file
vetoDefFile = inspiralutils.downloadVetoDefFile(cp, opts.generate_segments)

# Generate veto xml files for each ifo and for each category
inspiralutils.generate_veto_cat_files(cp, vetoDefFile, \
    opts.generate_veto_segments)

# Generate a combined veto file for each ifo for use by pipedown
gps_start = int(opts.gps_start_time)
gps_end = int(opts.gps_end_time)
duration = gps_end - gps_start

dqVetoes = {}
dqVetoes["combined-veto-file"] = {}
for category in veto_categories:
  ifos_string = ""
  veto_cat_files = []
  veto_cat_filename = str(gps_start) + "-" + str(duration) + ".xml"
  for ifo in ifos:
    veto_cat_files.append("-".join([ifo, "VETOTIME_CAT" + str(category), veto_cat_filename]))
    ifos_string += ifo
  dqVetoes["combined-veto-file"][category] = '-'.join([ifos_string, "VETOTIME_CAT_" + str(category), veto_cat_filename])
  concat_veto_call = " ".join([ cp.get("condor", "ligolw_add"),
    "--output", dqVetoes["combined-veto-file"][category] ] + veto_cat_files)
  inspiralutils.make_external_call(concat_veto_call)

  # FIXME: remove this block when the veto files from
  # ligolw_segments_from_cats have contain a segment table with start_time_ns
  # and end_time_ns, and when the process table doesn't need the domain, jobid,
  # and is_online columns added

  compat_veto_call = " ".join([cp.get("condor", "ligolw_segments_compat"), dqVetoes["combined-veto-file"][category]])
  inspiralutils.make_external_call(compat_veto_call)

segFile = {}

for ifo in ifos:
  segFile[ifo], dqVetoes[ifo] = inspiralutils.findSegmentsToAnalyze(
      cp, ifo, veto_categories,
      opts.generate_segments, opts.use_available_data,
      opts.generate_veto_segments)
  # set correct name as we're in the segments directory:
  segFile[ifo] = "../segments/" + segFile[ifo]
  cp.set("input", ifo.lower() + "-segments", segFile[ifo])
  for key in dqVetoes[ifo].keys():
    dqVetoes[ifo][key] = "../segments/" + dqVetoes[ifo][key]
for key in dqVetoes["combined-veto-file"].keys():
  dqVetoes["combined-veto-file"][key] = "../segments/" + dqVetoes["combined-veto-file"][key]

os.chdir("..")

hw_inj_dir = "hardware_injection_summary"
if opts.run_hardware_inj:
  inspiralutils.mkdir(hw_inj_dir)
  inspiralutils.get_hwinj_segments(cp,ifos,hw_inj_dir)

if opts.do_omega_scans:
  if not os.path.isdir("omega_setup"):
    os.mkdir("omega_setup")
  os.chdir("omega_setup")
  inspiralutils.omega_scan_setup(cp,ifos)
  os.chdir("..")

# create the list of files that will be concatenated into the ihope cache
cachelist = []
cachefilename = basename + ".cache"

##############################################################################
# Run lalapps_inspiral_hipe for datafind and template bank generation

if opts.run_datafind or opts.run_tmpltbank:
  print("Running inspiral hipe for datafind/template bank run")
  hipeDfNode = inspiralutils.hipe_setup("datafind", cp, ifos, \
      opts.log_path, dataFind = opts.run_datafind, tmpltBank = opts.run_tmpltbank,
      dax=opts.dax, local_exec_dir=tmp_exec_dir, static_pfn_cache=peg_frame_cache,
      reuse_data=opts.reuse_data)
  for cacheFile in hipeDfNode.get_output_files():
      cachelist.append("datafind/" + cacheFile)
  dag.add_node(hipeDfNode)

if opts.complete_cache:
  # We need to get the names of the template bank files that are generated in
  # the datafind step, so we can link them into all the search directories
  datafind_cache_filename = "datafind/" + inspiralutils.hipe_cache(ifos, \
      None, opts.gps_start_time, opts.gps_end_time)
  tmpltbank_cache = inspiralutils.tmpltbank_cache(datafind_cache_filename)

##############################################################################
# Set up the directories for each run and run lalapps_inspiral_hipe

if opts.complete_cache:
  # dictionaries for the each type of run; keys will be the veto categories
  hipePlayVetoNode = {}
  hipeAnalysisVetoNode = {}
  hipeInjVetoNode = {}

  # define directories for hipe to run
  playDir = "playground"
  fullDir = "full_data"
  injSection = "injections"

  if opts.reverse_analysis:
    print("Setting up dags for a reverse chirp analysis")
    cp.set("inspiral", "reverse-chirp-bank","")
    # Add reverse to directories to specify that a reverse chirp search is being done
    playDir += "_reverse"
    fullDir += "_reverse"
    injSection += "-reverse"

  # loop over the veto categories
  for category in veto_categories:
    if not opts.run_data_quality and category > 1:
      break

    if opts.run_search:
      print("Setting up the category " + str(category) + " veto dags")

    if opts.run_playground:
      if opts.run_search:
        print("Running inspiral hipe for playground with vetoes")
        hipePlayVetoNode[category] = inspiralutils.hipe_setup(
            playDir, cp, ifos, opts.log_path, \
            playOnly = True, vetoCat = category, vetoFiles = dqVetoes, \
            dax=opts.dax, tmpltbankCache = tmpltbank_cache, \
            local_exec_dir=tmp_exec_dir, static_pfn_cache=peg_frame_cache,
            reuse_data=opts.reuse_data)
        for cacheFile in hipePlayVetoNode[category].get_output_files():
          cachelist.append(playDir + "/" + cacheFile)
        dag.add_node(hipePlayVetoNode[category])
        if opts.run_datafind or opts.run_tmpltbank:
          hipePlayVetoNode[category].add_parent(hipeDfNode)
        if category > 1:
          hipePlayVetoNode[category].add_parent(hipePlayVetoNode[1])

      # make cache of necessary files if not running search
      else:
        if cp.get("pipeline","user-tag"):
          usertag = cp.get("pipeline", "user-tag") + "_" + playDir.upper()
        else:
          usertag = playDir.upper()
        if category > 1:
          usertag += "_CAT_" + str(category) + "_VETO"
        cacheFile = inspiralutils.hipe_cache( ifos,usertag, \
            cp.getint("input", "gps-start-time"), \
            cp.getint("input", "gps-end-time"))
        if os.path.isfile(playDir + "/" + cacheFile):
          cachelist.append(playDir + "/" + cacheFile)
        else:
          print("WARNING: Cache file " + playDir + "/" + cacheFile, file=sys.stderr)
          print("does not exist! This might cause later failures.", file=sys.stderr)


    if opts.run_full_data:
      if opts.run_search:
        print("Running inspiral hipe for full data with vetoes")
        hipeAnalysisVetoNode[category] = inspiralutils.hipe_setup(
            fullDir, cp, ifos, opts.log_path, \
            vetoCat = category, vetoFiles = dqVetoes, \
            dax = opts.dax, tmpltbankCache = tmpltbank_cache, \
            local_exec_dir=tmp_exec_dir, \
            data_checkpoint = opts.data_checkpoint,
            static_pfn_cache=peg_frame_cache,
            reuse_data=opts.reuse_data)
        for cacheFile in hipeAnalysisVetoNode[category].get_output_files():
          cachelist.append(fullDir + "/" + cacheFile)
        dag.add_node(hipeAnalysisVetoNode[category])
        if opts.run_datafind or opts.run_tmpltbank:
          hipeAnalysisVetoNode[category].add_parent(hipeDfNode)
        if category > 1:
          hipeAnalysisVetoNode[category].add_parent(hipeAnalysisVetoNode[1])

      # make cache of necessary files if not running search  
      else:
        if cp.get("pipeline","user-tag"):
          usertag = cp.get("pipeline", "user-tag") + "_" + fullDir.upper()
        else:
          usertag = fullDir.upper()
        if category > 1:
          usertag += "_CAT_" + str(category) + "_VETO"
        cacheFile = inspiralutils.hipe_cache( ifos,usertag, \
            cp.getint("input", "gps-start-time"), \
            cp.getint("input", "gps-end-time"))
        if os.path.isfile(fullDir + "/" + cacheFile):
          cachelist.append(fullDir + "/" + cacheFile)
        else:
          print("WARNING: Cache file " + fullDir + "/" + cacheFile, file=sys.stderr)
          print("does not exist! This might cause later failures.", file=sys.stderr)


    if opts.run_injections:
      # dictionary for individual injection runs
      hipeInjVetoNode[category] = {}

      for (injDir, injSeed) in cp.items( injSection ):
        if opts.run_search:
          print("Running inspiral hipe for " + injDir + " with vetoes")
          hipeInjVetoNode[category][injDir] = inspiralutils.hipe_setup(
              injDir, cp, ifos, opts.log_path, \
              injSeed = injSeed, vetoCat = category, vetoFiles = dqVetoes, \
              dax=opts.dax, tmpltbankCache = tmpltbank_cache, \
              local_exec_dir = tmp_exec_dir, \
              data_checkpoint = opts.data_checkpoint,
              static_pfn_cache=peg_frame_cache,
              reuse_data=opts.reuse_data)
          for cacheFile in hipeInjVetoNode[category][injDir].get_output_files():
            cachelist.append(injDir + "/" + cacheFile)
          dag.add_node(hipeInjVetoNode[category][injDir])
          if opts.run_datafind or opts.run_tmpltbank:
            hipeInjVetoNode[category][injDir].add_parent(hipeDfNode)
          if category > 1:
            hipeInjVetoNode[category][injDir].add_parent(hipeInjVetoNode[1][injDir])

        # make cache of necessary files if not running search
        else:
          if cp.get("pipeline","user-tag"):
            usertag = cp.get("pipeline", "user-tag") + "_" + injDir.upper()
          else:
            usertag = injDir.upper()
          if category > 1:
            usertag += "_CAT_" + str(category) + "_VETO"
          cacheFile = inspiralutils.hipe_cache( ifos,usertag, \
              cp.getint("input", "gps-start-time"), \
              cp.getint("input", "gps-end-time"))
          if os.path.isfile(injDir + "/" + cacheFile):
            cachelist.append(injDir + "/" + cacheFile)
          else:
            print("WARNING: Cache file " + injDir + "/" + cacheFile, file=sys.stderr)
            print("does not exist! This might cause later failures.", file=sys.stderr)


##############################################################################
# cat the cache files together
  
if len(cachelist):
  command = "cat " + " ".join(cachelist) + " > " + cachefilename
  os.popen(command)


##############################################################################
# Run lalapps_pipedown

if opts.run_pipedown:
  # Set up parents correctly
  parentNodes = []
  if opts.run_search:
    for category in veto_categories:
      if opts.run_playground:
        parentNodes.append(hipePlayVetoNode[category])
      if opts.run_full_data:
        parentNodes.append(hipeAnalysisVetoNode[category])
      if opts.run_injections:
        for (injDir, injSeed) in cp.items( injSection ):
          parentNodes.append(hipeInjVetoNode[category][injDir])

  if len(parentNodes) == 0:
    parentNodes = None

  playgroundOnly = False

  print("Running lalapps_pipedown")
  dag = inspiralutils.pipedownSetup(dag,cp,opts.log_path,"pipedown",\
                            "../" + cachefilename,parentNodes,playgroundOnly,opts.run_mvsc)


##############################################################################
# Set up the directories for plotting and run lalapps_plot_hipe

# Plotting dag setup still treats CAT1 differently from other CATs
# hence need to hack the parent nodes somewhat

if opts.run_playground and opts.run_data_quality:
  hipePlayVetoNodeWithout1 = dict( (cat, hipePlayVetoNode[cat]) \
      for cat in vetocats_without1 )
if opts.run_full_data and opts.run_data_quality:
  hipeAnalysisVetoNodeWithout1 = dict( (cat, hipeAnalysisVetoNode[cat]) \
    for cat in vetocats_without1 )

if opts.run_plots:

  # playground plots
  if opts.run_playground:

    # set up parents
    if opts.run_search:
      # define hipePlayNode as parent for playground CAT1 plotting nodes
      hipePlayNode = hipePlayVetoNode[1]
      parentNodes = [hipePlayNode]
    else: parentNodes = None

    if opts.run_data_quality:
      parentVetoNodes = [hipePlayVetoNodeWithout1]
    else: parentVetoNodes = None
    
    print("Making plots depending only on the playground")
    dag = inspiralutils.zeroSlidePlots(dag, "playground_summary_plots", cp, \
        opts.log_path, "PLAYGROUND", "PLAYGROUND", "../" + cachefilename, \
        opts.do_dag_categories, parentNodes, parentVetoNodes, \
        vetocats_without1, ifos)

  # full data plots
  if opts.run_full_data:

    # set up parents
    if opts.run_search:
      # define hipeAnalysisNode as parent for full_data CAT1 plotting nodes
      hipeAnalysisNode = hipeAnalysisVetoNode[1]
      parentNodes = [hipeAnalysisNode]
    else: parentNodes = None

    if opts.run_data_quality:
      parentVetoNodes = [hipeAnalysisVetoNodeWithout1]
    else: parentVetoNodes = None

    print("Making full data zero lag plots")
    dag = inspiralutils.zeroSlidePlots(dag, "full_data_summary_plots", cp, \
        opts.log_path, "FULL_DATA", "FULL_DATA", "../" + cachefilename, \
        opts.do_dag_categories, parentNodes, parentVetoNodes, \
        vetocats_without1, ifos)

  # playground/full-data-slide plots
  if opts.run_playground and opts.run_full_data:

    # set up parents
    if opts.run_search:
      parentNodes = [hipePlayNode, hipeAnalysisNode]
    else: parentNodes = None

    if opts.run_data_quality:
      parentVetoNodes = [hipePlayVetoNodeWithout1, hipeAnalysisVetoNodeWithout1]
    else: parentVetoNodes = None

    print("Making plots with full data slides, playground zero lag")
    dag = inspiralutils.zeroSlidePlots(dag, "full_data_slide_summary_plots", \
        cp, opts.log_path, "PLAYGROUND", "FULL_DATA", "../" + cachefilename, \
        opts.do_dag_categories, parentNodes, parentVetoNodes, \
        vetocats_without1, ifos)

  # software injection run plots
  if opts.run_injections:

    hipeInjNode = {}

    print("Making plots depending on the injection runs")

    if opts.run_full_data:
      slideSuffix = "FULL_DATA"
    elif opts.run_playground:
      slideSuffix = "PLAYGROUND"
    else:
      slideSuffix = "NONE_AVAILABLE"
    if opts.run_playground:
      zeroLagSuffix = "PLAYGROUND"
    else:
      zeroLagSuffix = "NONE_AVAILABLE"

    for (injDir, injSeed) in cp.items("injections"):

      # similar hack to separate out the CAT1 injection parent node
      hipeInjNode[injDir] = hipeInjVetoNode[1][injDir]

      # set up parents
      if opts.run_search:
        parentNodes = [hipeInjNode[injDir]]
        if opts.run_full_data:
          parentNodes.append(hipeAnalysisNode)
        if opts.run_playground:
          parentNodes.append(hipePlayNode)
      else: parentNodes = None

      # use the last veto category node as parent
      if opts.run_data_quality:
        parentVetoNodes = [hipeInjVetoNode[veto_categories[-1]][injDir]]
        if opts.run_full_data:
          parentVetoNodes.append(hipeAnalysisVetoNode[veto_categories[-1]])
        if opts.run_playground:
          parentVetoNodes.append(hipePlayVetoNode[veto_categories[-1]])
      else: parentVetoNodes = None

      dag = inspiralutils.injZeroSlidePlots(dag, injDir + "_summary_plots", \
          cp, opts.log_path, injDir.upper(), zeroLagSuffix, slideSuffix, \
          "../" + cachefilename, opts.do_dag_categories, parentNodes, 
          parentVetoNodes, vetocats_without1, ifos)

    print("Making plots of all injection runs together")

    # set up parents
    if opts.run_search:
      parentNodes = hipeInjNode.values()
      if opts.run_full_data:
        parentNodes.append(hipeAnalysisNode)
      if opts.run_playground:
        parentNodes.append(hipePlayNode)
    else: parentNodes = None

    if opts.run_data_quality:
      parentVetoNodes = hipeInjVetoNode[veto_categories[-1]].values()
      if opts.run_playground:
        parentVetoNodes.append(hipePlayVetoNode[veto_categories[-1]])
      if opts.run_full_data:
        parentVetoNodes.append(hipeAnalysisVetoNode[veto_categories[-1]])
    else: parentVetoNodes = None

    dag = inspiralutils.injZeroSlidePlots(dag, "allinj_summary_plots", cp, \
        opts.log_path, "*INJ", zeroLagSuffix, slideSuffix, \
        "../" + cachefilename, opts.do_dag_categories, parentNodes, 
        parentVetoNodes, vetocats_without1, ifos)

  if opts.do_dag_categories:
    dag.add_maxjobs_category('plotting',2)

##############################################################################
# Set up jobs for HW injection page

if opts.run_hardware_inj:

  print("Setting up hardware injection summary jobs")

  # set up parents
  parentNodes = []
  if opts.run_search:
    if opts.run_full_data:
      parentNodes.append(hipeAnalysisVetoNode[1])

    if opts.run_data_quality:
      for category in vetocats_without1:
        parentNodes.append(hipeAnalysisVetoNode[category])

  HWinjNodes = inspiralutils.hwinj_page_setup(cp,ifos,vetocats_without1,hw_inj_dir)
  for node in HWinjNodes:
    for parent in parentNodes:
      node.add_parent(parent)
    dag.add_node(node)

##############################################################################
# Set up the directories for followup and run lalapps_followup_pipe
if opts.run_followup:
  followupNode = inspiralutils.followup_setup("followup", cp, opts, "full_data")
  # FIXME: the cache files need to be appended to cachelist
  dag.add_node(followupNode)
  if opts.run_search:
    followupNode.add_parent(hipeAnalysisNode)

##############################################################################
# set the number of retries for each of the sub-dags run by ihope
try:
  num_retries = int(cp.get("pipeline","retry-subdag"))
  for node in dag.get_nodes():
    node.set_retry(num_retries)
except:
  # do not retry nodes
  pass

##############################################################################
# Write the dag and sub files

# see if any of the jobs use remote sites
exec_site_dict = {}
for sec in cp.sections():
  if cp.has_option(sec,'remote-sites'):
    exec_site_dict[cp.get(sec,'remote-sites')] = True

exec_site_list = exec_site_dict.keys()
if exec_site_list:
  exec_site = ','.join(['local'] + exec_site_list)
else:
  exec_site=None

# write the workflow control files
dag.prepare_dax(tmp_exec_dir=tmp_exec_dir,grid_site=exec_site)
dag.write_sub_files()
dag.write_dag()

print() 
print("Created a workflow file which can be submitted by executing")
print("\n    cd " + analysisDirectory)
print("""
and then

    ./pegasus_submit_dax

from a condor submit machine. 

From the analysis directory on the condor submit machine, you can run the
command

    pegasus-status --long -t -i `./pegasus_basedir`

to check the status of your workflow. Once the workflow has finished you 
can run the command

    pegasus-analyzer -t -i `./pegasus_basedir`

to debug any failed jobs.
""")

##############################################################################
# write out a log file for this script
log_fh = open(basename + '.pipeline.log', 'w')
  
# FIXME: the following code uses obsolete CVS ID tags.
# It should be modified to use git version information.
log_fh.write( "$Id$" \
    + "\n" )
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

##############################################################################

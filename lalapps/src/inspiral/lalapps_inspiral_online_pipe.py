"""
inspiral_online_pipe.py - online inspiral pipeline driver script

This script produced the necessary condor submit and dag files to run
a prototype online analysis in S6
"""

from __future__ import print_function

from pylal import git_version


# import standard modules
import sys, os
import popen2, time
import getopt, re, string
import socket
from six.moves import configparser
import urlparse

# import the modules we need to build the pipeline
from ligo import segments
from glue import pipeline
from lalapps import inspiral


class TrigscanJob(inspiral.InspiralAnalysisJob):
    """
    Stand-alone trigscan job
    """
    
    # lalapps_trigscan --ifo H1 --gps-start-time 939296936 --gps-end-time 939298856 --ts-cluster T0T3TcAS --ts-metric-scaling 0.12 H1-INSPIRAL_UNCLUSTERED-939296936-1920.xml

    def __init__(self,cp,dax=False):
      """
      cp = ConfigParser object from which options are read.
      """
      exec_name = 'trigscan'
      sections = ['trigscan']
      extension = 'xml'
      inspiral.InspiralAnalysisJob.__init__(self,cp,sections,exec_name,extension,dax)


class SiClusterJobCoarse(inspiral.InspiralAnalysisJob):
    # lalapps_trigscan --ifo H1 --gps-start-time 939296936 --gps-end-time 939298856 --ts-cluster T0T3TcAS --ts-metric-scaling 0.12 H1-INSPIRAL_UNCLUSTERED-939296936-1920.xml

    def __init__(self,cp,dax=False):
      """
      cp = ConfigParser object from which options are read.
      """
      exec_name = 'siclustercoarse'
      sections = ['siclustercoarse']
      extension = 'xml'
      inspiral.InspiralAnalysisJob.__init__(self,cp,sections,exec_name,extension,dax)

class SiClusterJobFine(inspiral.InspiralAnalysisJob):
    # lalapps_trigscan --ifo H1 --gps-start-time 939296936 --gps-end-time 939298856 --ts-cluster T0T3TcAS --ts-metric-scaling 0.12 H1-INSPIRAL_UNCLUSTERED-939296936-1920.xml

    def __init__(self,cp,dax=False):
      """
      cp = ConfigParser object from which options are read.
      """
      exec_name = 'siclusterfine'
      sections = ['siclusterfine']
      extension = 'xml'
      inspiral.InspiralAnalysisJob.__init__(self,cp,sections,exec_name,extension,dax)


class FilesystemJob(pipeline.CondorDAGJob):
    """
    Job for filesystem calls. Needed to move H1-TRIGSCAN... to H1-TS_CLUSTERED...
    and to copy TS_CLUSTERED to TSSI_CLUSTERED
    """
    
    def __init__(self, cmd):
        pipeline.CondorDAGJob.__init__(self,'vanilla','/bin/%s'  % cmd)
        self.set_stderr_file('logs/%s-$(cluster)-$(process).err' % cmd)
        self.set_stdout_file('logs/%s-$(cluster)-$(process).out' % cmd)
        self.set_sub_file('%s.sub' % cmd)



class MvNode(pipeline.CondorDAGNode):
    def __init__(self,job):
        pipeline.CondorDAGNode.__init__(self,job)



def get_all_files_in_range(dirname, starttime, endtime, pad=64):
    """Returns all files in dirname and all its subdirectories whose
    names indicate that they contain segments in the range starttime
    to endtime"""
    
    ret = []

    # Maybe the user just wants one file...
    if os.path.isfile(dirname):
        if re.match('.*-[0-9]*-[0-9]*\.gwf', dirname):
            return [dirname]
        else:
            return ret

    first_four_start = starttime / 100000
    first_four_end   = endtime   / 100000

    for filename in os.listdir(dirname):
        if re.match('.*-[0-9]{5}$', filename):
            dirtime = int(filename[-5:])
            if dirtime >= first_four_start and dirtime <= first_four_end:
                ret += get_all_files_in_range(os.path.join(dirname,filename), starttime, endtime, pad=pad)
        elif re.match('.*-[0-9]*-[0-9]*\.gwf', filename):
            file_time = int(filename.split('-')[-2])
            if file_time >= (starttime-pad) and file_time <= (endtime+pad):
                ret.append(os.path.join(dirname,filename))
        else:
            # Keep recursing, we may be looking at directories of
            # ifos, each of which has directories with times
            ret += get_all_files_in_range(os.path.join(dirname,filename), starttime, endtime, pad=pad)

    return ret


def get_valid_segments(segment_url, base_dir, ifo, science_flag, start_time, end_time):
    print("Finding valid analysis times for %s, please hold..." % ifo)

    cmd  = 'ligolw_segment_query --query-segments --segment-url %s --include-segments %s --gps-start-time %d --gps-end-time %d | ligolw_print -t segment -c start_time -c end_time' % (segment_url, science_flag, start_time, end_time)
    pipe = os.popen(cmd)

    print(cmd)

    results   = [x.strip().split(',') for x in pipe]
    science   = segments.segmentlist([segments.segment(int(x[0]), int(x[1])) for x in results])
    science.coalesce()

    print("Science: ")
    for s in science:
       print(s[0], s[1])

    framedir  = base_dir + '/' + ifo[0] + '1'
    chunks    = [f.split('.')[0].split('-') for f in get_all_files_in_range(framedir, start_time, end_time)]
    available = segments.segmentlist([ segments.segment( int(x[-2]), int(x[-2]) + int(x[-1]) ) for x in chunks if len(x) == 6 ])
    available.coalesce()

    print("Available:")
    for s in available:
       print(s[0], s[1])

    result = science & available

    result.coalesce()

    print("Result:")
    for s in result:
       print(s[0], s[1])

    print("done.")

    return result



def usage():
  msg = """\
Usage: lalapps_inspiral_online_pipe [options]

  -h, --help                     display this message
  -v, --version                  print version information and exit

  -s, --gps-start-time SEC       start time of segment to be analyzed
  -e, --gps-end-time SEC         end time of segment to be analyzed
  -a, --analysis-start-time SEC  only output triggers after GPS time
  -b, --analysis-end-time SEC    only output triggers before GPS time SEC

  -f, --dag-file-name NAME       output dag file name
  -t, --aux-data-path PATH       path to auxillary data

  -c, --config-file FILE         use configuration file FILE
  -l, --log-path PATH            directory to write condor log file
  -p, --coh-ptf                 Use coh_PTF_inspiral.c instead of inspiral.c
"""
  print(msg, file=sys.stderr)

# pasrse the command line options to figure out what we should do
shortop = "hvs:e:a:b:f:t:c:l:p:"
longop = [
  "help",
  "version",
  "gps-start-time=",
  "gps-end-time=",
  "analysis-start-time=",
  "analysis-end-time=",
  "dag-file-name=",
  "aux-data-path=",
  "config-file=",
  "log-path=",
  "coh-ptf"
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

gps_start_time = None
gps_end_time = None
analysis_start_time = None
analysis_end_time = None
dag_file_name = None
aux_data_path = None
config_file = None
log_path = None
doCohPTF = False

for o, a in opts:
  if o in ("-v", "--version"):
    print(git_version.verbose_msg)
    sys.exit(0)
  elif o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-s", "--gps-start-time"):
    gps_start_time = int(a)
  elif o in ("-e", "--gps-end-time"):
    gps_end_time = int(a)
  elif o in ("-a", "--analysis-start-time"):
    analysis_start_time = int(a)
  elif o in ("-b", "--analysis-end-time"):
    analysis_end_time = int(a)
  elif o in ("-f", "--dag-file-name"):
    dag_file_name = a
  elif o in ("-t", "--aux-data-path"):
    aux_data_path = a
  elif o in ("-c", "--config-file"):
    config_file = a
  elif o in ("-l", "--log-path"):
    log_path = a
  elif o in ("-p", "--coh-ptf"):
    doCohPTF = True
  else:
    print("Unknown option:", o, file=sys.stderr)
    usage()
    sys.exit(1)

if not gps_start_time:
  print("No GPS start time specified.", file=sys.stderr)
  print("Use --gps-start-time SEC to specify start time.", file=sys.stderr)
  sys.exit(1)

if not gps_end_time:
  print("No GPS end time specified.", file=sys.stderr)
  print("Use --gps-end-time SEC to specify end time.", file=sys.stderr)
  sys.exit(1)

if not dag_file_name:
  print("No DAG file name specified.", file=sys.stderr)
  print("Use --dag-file-name NAME to specify output DAG file.", file=sys.stderr)
  sys.exit(1)

if not aux_data_path:
  print("No auxiliary data path specified.", file=sys.stderr)
  print("Use --aux-data-path PATH to specify directory.", file=sys.stderr)
  sys.exit(1)

if not config_file:
  print("No configuration file specified.", file=sys.stderr)
  print("Use --config-file FILE to specify location.", file=sys.stderr)
  sys.exit(1)

if not log_path:
  print("No log file path specified.", file=sys.stderr)
  print("Use --log-path PATH to specify a location.", file=sys.stderr)
  sys.exit(1)

try: os.mkdir('cache')
except: pass
try: os.mkdir('logs')
except: pass

# create the config parser object and read in the ini file
cp = configparser.ConfigParser()
cp.read(config_file)

ifo_list = ['H1','L1','V1']

# get the usertag from the config file
try:
  usertag = string.strip(cp.get('pipeline','user-tag'))
  inspstr = 'INSPIRAL_' + usertag
except:
  usertag = None
  inspstr = 'INSPIRAL'

# create a log file that the Condor jobs will write to
basename = re.sub(r'\.dag',r'',dag_file_name)
logfile = os.path.join(log_path, basename + '.log')

# create a DAG generation log file
log_fh = open(basename + '.pipeline.log', 'w')

# create the DAG writing the log to the specified file
dag = pipeline.CondorDAG(logfile)
dag.set_dag_file(dag_file_name)

# Build a set of jobs for each ifo
for ifo in ifo_list:
    if not cp.has_option('segments',ifo.lower() + '-analyze'):
        continue
    
    # decide if we need to segment the data
    available_segments = get_valid_segments(cp.get('segfind','segment-url'), cp.get('framefind','base-dir'), ifo, cp.get('segments',ifo.lower() + '-analyze'), gps_start_time, gps_end_time)

    if not available_segments:
        print("No available segments for %s, skipping" % ifo)
        continue
    
    
    # create the Condor jobs that will be used in the DAG
    df_job    = pipeline.LSCDataFindJob('cache','logs',cp)
    tmplt_job = inspiral.TmpltBankJob(cp)

    # Based on S6A results ttrigscan clustering has
    # been replaced with 30-ms window clustering
    # ts_job = TrigscanJob(cp)

    si_job_coarse = SiClusterJobCoarse(cp)
    si_job_fine   = SiClusterJobFine(cp)
    cp_job        = FilesystemJob('cp')

    # Add ifo-specific template config
    if cp.has_section(ifo.lower() + '-tmpltbank'):
        tmplt_job.add_ini_opts(cp,ifo.lower() + '-tmpltbank')

    # Create a job to split the template into parallelization pieces
    split_job = inspiral.SplitBankJob(cp)

    # and one to run the inspirals
    if doCohPTF:
      insp_job  = inspiral.PTFInspiralJob(cp)
    else:
      insp_job  = inspiral.InspiralJob(cp)

    # template and inspiral jobs also need to know about the data
    if ifo == 'V1':
        channel_name = 'virgo-channel'
        data_name    = 'virgo-data'
        type_name    = 'virgo-type'
    else:
        channel_name = 'ligo-channel'
        data_name    = 'ligo-data'
        type_name    = 'ligo-type'

    channel = cp.get('input',channel_name)
    type    = cp.get('input',type_name)

    tmplt_job.set_channel(channel)
    tmplt_job.add_ini_opts(cp, data_name)

    # Add ifo-specific inspiral config
    if doCohPTF:
      insp_job.add_opt(ifo.lower() + '-data','')
      gps_times = str(gps_start_time) + '-' + str(gps_end_time - gps_start_time)
      if ifo == 'H1' or ifo == 'H2' or ifo == 'L1':
        insp_job.add_opt(ifo.lower()+'-channel-name',ifo.upper()+':' + cp.get('input','ligo-channel'))
        insp_job.add_opt(ifo.lower()+'-frame-cache','cache/' + ifo[0] + '-' + ifo + '_' + cp.get('input','ligo-type') + '_CACHE-' + gps_times + '.lcf')
      elif ifo == 'V1':
        insp_job.add_opt(ifo.lower()+'-channel-name',ifo.upper()+':' + cp.get('input','virgo-channel'))
        insp_job.add_opt(ifo.lower()+'-frame-cache','cache/' + ifo[0] + '-' + cp.get('input','virgo-type') + '_CACHE-' + gps_times + '.lcf')
    else:
      if cp.has_section(ifo.lower() + '-inspiral'):
          insp_job.add_ini_opts(cp,ifo.lower() + '-inspiral')

      insp_job.set_channel(channel)
      insp_job.add_ini_opts(cp, data_name)

      # First inspiral doesn't do vetos
      insp_job.add_ini_opts(cp, 'no-veto-inspiral')

    # add log configuration
    lwadd_job = pipeline.LigolwAddJob('logs',cp)

    # make the lwadd job a priority so results get reported quickly
    lwadd_job.add_condor_cmd('priority', '20')
    lwadd_job.add_opt('add-lfn-table','')

    # set the usertag in the template bank
    if usertag:
        tmplt_job.add_opt('user-tag',usertag)

    # set better submit file names than the default
    subsuffix = '.sub'
    df_job.set_sub_file(    basename + '_' + ifo + '.datafind'  + subsuffix )
    tmplt_job.set_sub_file( basename + '_' + ifo + '.tmpltbank' + subsuffix )
    split_job.set_sub_file( basename + '_' + ifo + '.splitbank' + subsuffix )
    insp_job.set_sub_file(  basename + '_' + ifo + '.inspiral'  + subsuffix )
    lwadd_job.set_sub_file( basename + '_' + ifo + '.ligolwadd' + subsuffix )

    split_job.add_opt('minimal-match', cp.get('tmpltbank','minimal-match'))

    # get the pad and chunk lengths from the values in the ini file
    pad = int(cp.get('data', 'pad-data'))
    n   = int(cp.get('data', 'segment-length'))
    s   = int(cp.get('data', 'number-of-segments'))
    r   = int(cp.get('data', 'sample-rate'))
    o   = int(cp.get('inspiral', 'segment-overlap'))
    length  = ( n * s - ( s - 1 ) * o ) / r
    overlap = o / r
    if doCohPTF:
      overlap = int(cp.get('coh_PTF_inspiral','segment-duration')) / 2
    job_analysis_time = length - overlap

    # find the data between the start time and the end time
    df = pipeline.LSCDataFindNode(df_job)
    df.set_start(gps_start_time)
    df.set_end(gps_end_time)
    df.set_observatory(ifo[0])
    if ifo == 'V1':
      df.set_type(type)
    else:
      df.set_type(ifo + '_' + type)
    dag.add_node(df)

    # modify the start and end time by pad seconds
    log_fh.write("gps_start_time = %d\n" % gps_start_time)
    log_fh.write("gps_end_time = %d\n" % gps_end_time)

    # Don't need to do these, since we'll pad each segment
    # gps_start_time += pad
    # gps_end_time   -= pad
    # log_fh.write("gps_start_time + pad = %d\n" % gps_start_time)
    # log_fh.write("gps_end_time - pad = %d\n" % gps_end_time)



    # Determine the times to analyze
    job_segs      = segments.segmentlist()
    analysis_segs = segments.segmentlist()

    for seg in available_segments:
        seg_start = seg[0] + pad
        seg_end   = seg[1] - pad

        if seg_end - seg_start < length:
            "Segment from %d to %d is too short, skipping" % (seg_start, seg_end)
            continue

        while seg_start + length < seg_end:
            job_segs.append( segments.segment(seg_start, seg_start + length) )
            analysis_segs.append( segments.segment(seg_start + overlap/2, seg_start + length - overlap/2) )

            seg_start += length - overlap

        # Get the last one that goes to the end of the segment
        if seg_start < seg_end:
            job_segs.append( segments.segment(seg_end - length, seg_end) )
            analysis_segs.append( segments.segment(seg_start + overlap/2, seg_end) )
  
    for seg, analysis_seg in zip(job_segs, analysis_segs):
        # create the template bank job
        bank = inspiral.TmpltBankNode(tmplt_job)
        bank.set_start(seg[0])
        bank.set_end(seg[1])
        bank.set_ifo(ifo)
        bank.set_cache(df.get_output())
        bank.add_parent(df)
        dag.add_node(bank)
    
        # split the template bank up into smaller banks
        split = inspiral.SplitBankNode(split_job)
        split.set_bank(bank.get_output())
        split.set_num_banks(cp.get('splitbank','number-of-banks'))
        split.add_parent(bank)
        dag.add_node(split)

        # create the inspiral jobs to do the analysis
        sub_insp = []

        for i,subbank in enumerate(split.get_output()):
            if doCohPTF:
              insp = inspiral.PTFInspiralNode(insp_job)
            else:
              insp = inspiral.InspiralNode(insp_job)
            insp.set_start(seg[0])
            insp.set_end(seg[1])
            insp.set_trig_start('trig-start-time')
            insp.set_trig_end('trig-end-time')

            insp.set_ifo(ifo)
            if doCohPTF:
              insp.set_no_spin_bank(subbank)
              insp.set_ifo_tag('FIRST')
            else:
              insp.set_cache(df.get_output())
              insp.set_bank(subbank)

            if usertag:
                insp.set_user_tag(usertag + "_%2.2d" % i)
            else:
                insp.set_user_tag("UNCLUSTERED_%2.2d" % i)

            if doCohPTF:
              insp.set_output()
            insp.add_parent(split)
            dag.add_node(insp)
            sub_insp.append(insp)


        # create a cache file name to hold the input to ligolw_add
        out_basename         = ifo + '-' + inspstr + '-' + str(analysis_seg[0]) + '-'
        out_basename        += str(abs(analysis_seg))
        insp_cache_file_name = out_basename + '.cache'
        insp_cache_file      = open(insp_cache_file_name, 'w')
        path                 = 'file://localhost' + os.getcwd()

        output_name          = '%s-INSPIRAL_UNCLUSTERED-%d-%d.xml.gz' % (ifo, analysis_seg[0], abs(analysis_seg))

        # use ligolw_add to join everything back together
        lwadd = pipeline.LigolwAddNode(lwadd_job)
        lwadd.add_var_opt('input-cache',insp_cache_file_name)
        lwadd.set_output(output_name)
        lwadd.add_var_opt('lfn-start-time',seg[0])
        lwadd.add_var_opt('lfn-end-time',seg[1])

        for insp in sub_insp:
            file = insp.get_output()
            tmp  = re.split("[.-]", file)

            if len(tmp) == 5:
                [instrument, type, start, duration, extension] = tmp
            else:
                [instrument, type, start, duration, extension, gz] = tmp
                extension = extension + "." + gz

            print(instrument[0], type, start, duration, \
                     os.path.join(path , file), file=insp_cache_file)

            lwadd.add_parent(insp)

        output_file = lwadd.get_output()
        #[instrument, type, start, duration, extension] = re.split("[.-]", output_file)
        #lwadd.add_var_opt('lfn-comment',type)
        #start_str = str(start)
        dag.add_node(lwadd)

        # TrigScan's out, 30-ms window clustering is in
        # tsnode = inspiral.InspiralAnalysisNode(ts_job)
        # tsnode.set_start(analysis_seg[0])
        # tsnode.set_end(analysis_seg[1])
        # tsnode.add_var_opt('ifo',ifo)
        # tsnode.add_file_arg(output_name)
        # 
        # tsnode.add_parent(lwadd)
        # dag.add_node(tsnode)

        # sicluster works in-place, so copy unclustered triggers to 
        # new files for 30 ms and 16 sec clustering
        clustered_30ms_name = output_name.replace('UNCLUSTERED','30MILLISEC_CLUSTERED')
        clustered_16s_name  = output_name.replace('UNCLUSTERED','16SEC_CLUSTERED')

        for cname in [clustered_30ms_name, clustered_16s_name]:
            cpnode = pipeline.CondorDAGNode(cp_job)
            cpnode.add_file_arg(output_name)
            cpnode.add_file_arg(cname)
            cpnode.add_parent(lwadd)
            dag.add_node(cpnode)

            if cname == clustered_16s_name:
                sinode = inspiral.InspiralAnalysisNode(si_job_coarse)
            else:
                sinode = inspiral.InspiralAnalysisNode(si_job_fine)

            sinode.add_file_arg(cname)
            sinode.add_parent(cpnode)
            dag.add_node(sinode)
        
# write the dag
dag.write_sub_files()
dag.write_script()
dag.write_dag()

log_fh.close()
  
# Fix up the types in the dag (how does real ihope do this?!)
in_f  = open('daily.dag')
out_f = open('daily.dag.fixed', 'w')

for a in in_f:
    if a.find('macroobservatory="H"') != -1:
        a = a.replace('macrotype="ER_C00_L1"','macrotype="H1_ER_C00_L1"')
        a = a.replace('macrotype="DMT_C00_L2"','macrotype="H1_DMT_C00_L2"')
    elif a.find('macroobservatory="L"') != -1:
        a = a.replace('macrotype="ER_C00_L1"','macrotype="L1_ER_C00_L1"')
        a = a.replace('macrotype="DMT_C00_L2"','macrotype="L1_DMT_C00_L2"')

    print(a, end=' ', file=out_f)

in_f.close()
out_f.close()
os.system('mv daily.dag.fixed daily.dag')

sys.exit(0)

from __future__ import print_function

from optparse import OptionParser

from matplotlib import use
use('Agg')

import sys
import os
import shutil

import pickle
import operator
import glob
import subprocess
import time
from random import uniform

import numpy.numarray as na

import pylab
from pylab import exp
from pylab import sqrt

from ligo.segments import segment, segmentlist

from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw import utils

from glue.segmentdb import segmentdb_utils

import sqlite3

params = {'axes.labelsize'  : 12,
          'font.size'       : 12,
          'text.fontsize'   : 12,
          'legend.fontsize' : 12,
          'xtick.labelsize' : 12,
          'ytick.labelsize' : 12,
          'text.usetex'     : True}

pylab.rcParams.update(params)

class DefaultContentHandler(ligolw.LIGOLWContentHandler):
    pass
lsctables.use_in(DefaultContentHandler)

###############
# Page template
###############
html_template = """
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">

<head>
  <link media="all" href="ihope_daily_style.css" type="text/css" rel="stylesheet" />
  <title>Low mass CBC - %(date)s</title>
  <script src="ihope_daily_toggle.js" type="text/javascript"></script>
</head>

<body>
  <div id="wrapper">
    <div id="menubar">
      <div id="menu">

  IFO:
  <form id="whichIFO">
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="H1" onclick="switchViewIframe()" checked>H1&nbsp;<input type="radio" name="box0" value="L1" onclick="switchViewIframe()">L1&nbsp;<input type="radio" name="box0" value="V1" onclick="switchViewIframe()">V1<br>
  </form>


  <p>Veto level:
  <form id="whichVetoLevel">
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="0" onclick="switchViewIframe()">Science<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="1" onclick="switchViewIframe()" checked>Cat 1<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="2" onclick="switchViewIframe()">Cats 1 and 2<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="4" onclick="switchViewIframe()">Cats 1,2 and 4<br>
  </form>

  <p>Cluster level:
  <form id="whichCluster">
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="UNCLUSTERED" onclick="switchViewIframe()" checked>Raw triggers<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="30MILLISEC_CLUSTERED" onclick="switchViewIframe()">30 msec window<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="16SEC_CLUSTERED" onclick="switchViewIframe()">16 sec window<br>
  </form>

  <p>Contents:
  <form id="whichContents">
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="Summary" onclick="switchViewIframe()" checked>Analysis time, veto <br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;usage<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="Glitch" onclick="switchViewIframe()">Loudest Triggers<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="Hwinj" onclick="switchViewIframe()">HW Injections<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="snr_hist" onclick="switchViewIframe()">SNR Histograms<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="glitchgram" onclick="switchViewIframe()">Glitchgram<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="rate_vs_time" onclick="switchViewIframe()">Rate vs. Time,<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Snr vs. Time,<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;high-rate times<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="template" onclick="switchViewIframe()">Breakdown by template<br>
    &nbsp;&nbsp;&nbsp;<input type="radio" name="box0" value="bankchisq" onclick="switchViewIframe()">&chi;<sup>2</sup><br>
  </form>



      </div>
    </div>

    <div id="ihope">
      <h2> i h o p e </h2>
      <img src="ihopeFrog.jpg">
    </div>

    <div id="header">
      <h1>Low mass CBC - %(date)s</h1>
      <h3>%(start_time)d &mdash; %(end_time)d</h3>
    </div>

    <iframe name="theframe" class="iframecontent" src="H1_CAT1_UNCLUSTERED_Summary.html" height="600">
    </iframe>

  </div>
</body>
</html>

"""


#############
# Constants
#############
MTSUN_SI = 4.92549095e-06

ifo_colored_circle = {'H1':'ro', 'L1':'go', 'V1':'mo'}
ifo_colors = {'H1':'red', 'L1':'green', 'V1':'magenta'}

##############################
# Dag-related classes and code
##############################

class Tree:
    id         = ''
    action     = ''
    config     = ''
    ifo        = ''
    level      = ''
    cluster    = ''
    start_time = ''
    end_time   = ''
    children   = []
    local      = False

    def __init__(self, action, config, ifo, level, cluster, start_time, end_time, local = False, flower = 40):
        self.id         = self.make_id()
        self.action     = action
        self.config     = config
        self.ifo        = ifo
        self.level      = level
        self.cluster    = cluster
        self.start_time = start_time
        self.end_time   = end_time
        self.children   = []
        self.local      = local
        self.flower     = flower

    def make_id(self):
        idchars = '0123456789ABCDEF'
        return ''.join([idchars[int(uniform(0,len(idchars)))] for i in range(32)])

    def add_child(self,child):
        self.children.append(child)
        return child

    def add_sub_job(self, action, config, ifo, level, cluster, start_time, end_time, local = False, flower=40):
        return self.add_child( Tree(action, config, ifo, level, cluster, start_time, end_time, local, flower) )


    def add_sub_jobs(self, *args):
        self.children += args
        return self

    def render_job(self, out, follow=True):
        if self.local:
            print('JOB %s daily_ihope_page_local.sub' % self.id, file=out)
        else:
            print('JOB %s daily_ihope_page.sub' % self.id, file=out)

        print('RETRY %s 3' % self.id, file=out)
        print('VARS %s macroaction="%s" macroconfig="%s" macroifo="%s" macrocategory="%d" macrocluster="%s" macrogpsstarttime="%s" macrogpsendtime="%s" macroflower="%f"' % (self.id, self.action, self.config['filename'], self.ifo, self.level, self.cluster, self.start_time, self.end_time, self.flower), file=out)

	if self.action in ['make_glitch_page','make_hwinj_page']:
	    print('CATEGORY %s database' % self.id, file=out)

        for c in self.children:
            c.render_job(out)

    def render_relationships(self, out):
        for c in self.children:
            print('PARENT %s CHILD %s' % (self.id, c.id), file=out)
            c.render_relationships(out)


def make_dag(config, ifo, category, cluster, start_time, end_time):
    try:
        os.makedirs(config['tmp_dir'])
    except:
        pass

    try:
        os.makedirs(config['tmp_dir'] + '/log')
    except:
        pass

    try:
        os.makedirs(config['tmp_dir'] + '/logs')
    except:
        pass

    try:
        os.makedirs(config['out_dir'])
    except:
        pass

    roots = []

    for ifo in ifos:
        veto_job = Tree('make_veto_files', config, ifo, 1, 'UNCLUSTERED', start_time, end_time, local = True)

        roots.append(veto_job)

        veto_job = veto_job.add_sub_job('make_veto_files', config, ifo, 2, 'UNCLUSTERED', start_time, end_time, local = True)
        veto_job = veto_job.add_sub_job('make_veto_files', config, ifo, 4, 'UNCLUSTERED', start_time, end_time, local = True)

        for cluster in ['UNCLUSTERED','30MILLISEC_CLUSTERED','16SEC_CLUSTERED']:
            root     = veto_job.add_sub_job('make_caches',   config, ifo, category, cluster, start_time, end_time)
            csv_jobs = [root.add_sub_job('make_csv',         config, ifo, 0, cluster, start_time, end_time)]
            csv_jobs += [csv_jobs[0].add_sub_job('make_csv', config, ifo, i, cluster, start_time, end_time) for i in [1,2,4]]

            for job in csv_jobs:
                job.add_sub_job('plot_rate_vs_time',        config, ifo, job.level, cluster, start_time, end_time)
                job.add_sub_job('plot_snr_hist',            config, ifo, job.level, cluster, start_time, end_time)
                job.add_sub_job('make_glitchy_times_table', config, ifo, job.level, cluster, start_time, end_time)
                job.add_sub_job('make_usage_table',         config, ifo, job.level, cluster, start_time, end_time)
                job.add_sub_job('make_summary_table',       config, ifo, job.level, cluster, start_time, end_time)

                job.add_sub_job('plot_snr_vs_time',         config, ifo, job.level, cluster, start_time, end_time)
                job.add_sub_job('plot_template_counts',     config, ifo, job.level, cluster, start_time, end_time)
                job.add_sub_job('plot_mass_hist',           config, ifo, job.level, cluster, start_time, end_time, flower = options.flower)
                job.add_sub_job('plot_glitchgram',          config, ifo, job.level, cluster, start_time, end_time)

                # Ugh.  This one is actually a child of all the csv jobs, needs some tweaking...
                # job.add_sub_job('plot_rate_vs_time_all',    config, ifo, job.level, cluster, start_time, end_time)

                job.add_sub_job('plot_hexmass',             config, ifo, job.level, cluster, start_time, end_time)
                job.add_sub_job('plot_bank_veto',           config, ifo, job.level, cluster, start_time, end_time)

                job.add_sub_job('make_hwinj_page',          config, ifo, job.level, cluster, start_time, end_time, local = True)

                if cluster == '16SEC_CLUSTERED':
                    dbjob = job.add_sub_job('make_database',    config, ifo, job.level, cluster, start_time, end_time)
                    dbjob.add_sub_job('make_glitch_page',       config, ifo, job.level, cluster, start_time, end_time, local = True)
                else:
                     job.add_sub_job('make_glitch_page',        config, ifo, job.level, cluster, start_time, end_time, local = True)

    out = open(config['tmp_dir'] + '/daily_page.dag','w')
    for root in roots:
        root.render_job(out)

    for root in roots:
        root.render_relationships(out)

    print("MAXJOBS database 4", file=out)
    out.close()


    out = open(config['tmp_dir'] + '/daily_ihope_page.sub','w')
    print("""universe = vanilla
executable = %s
arguments = " --config $(macroconfig) --action $(macroaction) --ifos $(macroifo) --veto-categories $(macrocategory) --cluster-categories $(macrocluster) --gps-start-time $(macrogpsstarttime) --gps-end-time $(macrogpsendtime) --flower $(macroflower)"
priority     = 20
getenv       = True
log          = %s/log/daily_plot.log
error        = logs/daily_plot-$(cluster)-$(process).err
output       = logs/daily_plot-$(cluster)-$(process).out
notification = never
queue 1""" % (config['ihope_daily_page'], config['tmp_dir']), file=out)
    out.close()

    out = open(config['tmp_dir'] + '/daily_ihope_page_local.sub','w')
    print("""universe = local
executable = %s
arguments = " --config $(macroconfig) --action $(macroaction) --ifos $(macroifo) --veto-categories $(macrocategory) --cluster-categories $(macrocluster) --gps-start-time $(macrogpsstarttime) --gps-end-time $(macrogpsendtime) "
priority = 20
getenv = True
log = %s/log/daily_plot.log
error = logs/daily_plot-$(cluster)-$(process).err
output = logs/daily_plot-$(cluster)-$(process).out
notification = never
queue 1""" % (config['ihope_daily_page'], config['tmp_dir']), file=out)
    out.close()



# Command line
def parse_command_line():
    """
    Parse the command line, return an options object
    """

    parser = OptionParser(
        version = "%prog CVS $Header$",
        usage   = "%prog --action action_name --gps-start-time start_time --ifos ifos --veto-categories vetoes --cluster-categories clusters [--gps-end-time end_time] [--config config_file]",
        description = "Performs the named action from start time to end time (default = start + 1 day)\n"
        "on the named ifos, at the named veto and cluster categories"
	)

    parser.add_option("-a", "--action",             metavar = "action",      help = "Action to perform")
    parser.add_option("-f", "--config",             metavar = "config_file", help = "Auxiliary config file")
    parser.add_option("-i", "--ifos",               metavar = "ifos",        help = "Comma-separated list of IFOs")
    parser.add_option("-v", "--veto-categories",    metavar = "veto_categories",      help = "Comma-separated list of veto levels (0=SCIENCE)")
    parser.add_option("-c", "--cluster-categories", metavar = "cluster_categories",   help = "Comma-separated list of clusters")

    parser.add_option("-s", "--gps-start-time",    metavar = "gps_start_time",    help="")
    parser.add_option("-e", "--gps-end-time",      metavar = "gps_end_time",      help="")
    parser.add_option("-F", "--flower", type=float, metavar = 'flower', default =40., help="The lower frequency cutoff used")

    options, filenames = parser.parse_args()

    if not options.action:
        print("Please specify --action", file=sys.stderr)
        sys.exit(-1)
    if not options.ifos:
        print("Please specify --ifos", file=sys.stderr)
        sys.exit(-1)

    return options


#########################
#
# Utilities
#
#########################
def parse_config(filename):
    if not os.path.exists(filename):
        print("Config file %s does not exist" % filename, file=sys.stderr)
        sys.exit(1)

    if not os.path.isfile(filename):
        print("Config file %s is not a file" % filename, file=sys.stderr)
        sys.exit(1)

    ret = {'filename':filename}
    f   = open(filename)

    for l in f:
        l = l.strip()

        if not l:
            continue

        if l[0] == '#':
            continue

        key, value = l.split('=')
        ret[key.strip()] = value.strip()

    return ret




def m1_m2_to_mtotal_eta(m1, m2):
    mtot = m1+m2
    eta  = (m1 * m2) / (mtot * mtot)

    return mtotal, eta

def make_timeless_title(ifo, veto_level, cluster, start_time, desc = None):
    veto_level_labels = ['Science', 'Cat 1', "Cat 1 & 2", "Cat 1,2 & 3", "Cat 1,2 & 4"]
    cluster_labels    = {'UNCLUSTERED':'Unclustered',
                         '30MILLISEC_CLUSTERED':'30 msec window',
                         '16SEC_CLUSTERED':'16 sec window'}


    if desc:
        return '%s (%s, %s, %s)' % (desc, ifo, veto_level_labels[veto_level], cluster_labels[cluster])
    else:
        return '%s, %s, %s' % (ifo, veto_level_labels[veto_level], cluster_labels[cluster])


def make_title(ifo, veto_level, cluster, start_time, desc = None):
    veto_level_labels = ['Science', 'Cat 1', "Cat 1 & 2", "Cat 1, 2 & 3", "Cat 1,2 & 4"]
    cluster_labels    = {'UNCLUSTERED':'Unclustered',
                         '30MILLISEC_CLUSTERED':'30 msec window',
                         '16SEC_CLUSTERED':'16 sec window'}


    tstring =  os.popen('lalapps_tconvert -f "%b %d, %Y" ' + str(start_time)).readlines()[0].strip()

    if desc:
        return '%s (%s: %s, %s, %s)' % (desc, tstring, ifo, veto_level_labels[veto_level], cluster_labels[cluster])
    else:
        return '%s: %s, %s, %s' % (tstring, ifo, veto_level_labels[veto_level], cluster_labels[cluster])


def new_snr(snr, chisq, chisq_dof, index=6.0):
    rchisq = chisq/(2 * chisq_dof - 2)
    nhigh  = 2.

    if rchisq > 1.:
        return snr / ((1. + rchisq**(index/nhigh))/2)**(1./index)
    else:
        return snr


def make_time_label(start_time):
    tstring =  os.popen('lalapps_tconvert -f "%b %d %Y" ' + str(start_time)).readlines()[0].strip()

    return 'Time [hours since %s 00:00:00 UTC]' % tstring

def chirplen(fstart, mtot, eta):
    c0 = 5*mtot*MTSUN_SI/(256*eta);
    c2 = 743.0/252.0 + eta*11.0/3.0;
    c3 = -32. * pylab.pi / 5.;
    c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0);
    x  = pow(pylab.pi * mtot * MTSUN_SI * fstart, 1.0/3.0);
    x2 = x*x;
    x3 = x*x2;
    x4 = x2*x2;
    x8 = x4*x4;

    chirpTime = c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8;

    return chirpTime


#########################
#
# Methods that set up data
#
#########################
def make_veto_files(engine, config, ifo, category, cluster, start_time, end_time):
    """
    Creates VETOCATS_(ifo)_(cat).p, which contains a dicionary mapping segment_definer
    names to segmentlists containing their active times.
    """

    def make_segdef(row, overall_interval):
        veto_interval = overall_interval & segmentlist([segment(row.start_time, row.end_time or end_time)])

        if not veto_interval:
            veto_interval = [[0,0]]

        return (row.ifo, row.name, row.version, veto_interval[0][0], veto_interval[0][1], row.start_pad, row.end_pad)


    overall       = segmentlist([segment(start_time, end_time)])
    veto_definers = []

    for f in config['veto_definer_file'].split(','):
        xmldoc  = utils.load_url(f, contenthandler=DefaultContentHandler)

        veto_definers += lsctables.VetoDefTable.get_table(xmldoc)

    cats = [make_segdef(row, overall) for row in veto_definers if row.category == category and row.ifo == ifo]

    vetoed_segments = segmentdb_utils.query_segments(engine, 'segment', cats)
    veto_arrays     = [ [(x[0], x[1]) for x in y] for y in vetoed_segments ]
    result          = map(lambda x,y: ('%s:%d' % (x[1], x[2]), y), cats, veto_arrays)

    f = open('%s/VETOCATS_%s_%d.p' % (config['tmp_dir'], ifo, category),'w')
    pickle.dump(result,f)
    f.close()


def make_caches(config, ifo, level, cluster, start, end):
    """
    Creates (ifo)-INSPIRAL_(cluster).cache and (ifo)-TMPLTBANK_(cluster).cache suitable for pasing to ligolw_print.

    This only needs the ifo and cluster, other arguments are provided so that all methods can be called uniformly
    """
    cmd   = 'lalapps_path2cache'

    # sngl_inspirals
    files = sorted(glob.glob('%s/%s-INSPIRAL_%s-*xml*' % (config['trigger_dir'], ifo, cluster)))
    proc  = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc.communicate("\n".join(files))

    sys.stdout.flush()

    f = open('%s/%s-INSPIRAL_%s.cache' % (config['tmp_dir'], ifo, cluster),'w')
    print(out.strip(), file=f)
    f.close()

    # templates
    files = sorted(glob.glob('%s/%s-TMPLTBANK-*xml*' % (config['trigger_dir'],ifo)))
    proc  = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = proc.communicate("\n".join(files))

    f = open('%s/%s-TMPLTBANK_%s.cache' % (config['tmp_dir'], ifo, cluster),'w')
    print(out.strip(), file=f)
    f.close()


def make_csv(config, ifo, category, cluster, start, end):
    """
    At category 0, reads (ifo)-INSPIRAL_(cluster).cache and (ifo)-TMPLTBANK_(cluster).cache
    and produces (ifo)-0-INSPIRAL_(cluster).csv,  (ifo)-TMPLTBANK.csv and (ifo)-0-SUMMARY_(cluster).csv

    At other categories, reads (ifo)-0-INSPIRAL_(cluster).csv and (ifo)-0-SUMMARY_(cluster).csv and
    VETOCATS_(ifo)_(category).p and produces (ifo)-(cat)-INSPIRAL.csv and (ifo)-(cat)-SUMMARY.csv
    """
    if category == 0:
        # ligolw_print doesn't like empty cache files
        cache_size = os.path.getsize('%s/%s-INSPIRAL_%s.cache' % (config['tmp_dir'], ifo, cluster))

        if cache_size < 10:
            out_f = open('%s/%s-0-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, cluster),'w')
            out_f.close()
            out_f = open('%s/%s-0-SUMMARY_%s.csv' % (config['tmp_dir'], ifo, cluster), 'w')
            out_f.close()
            out_f = open('%s/%s-TMPLTBANK_%s.csv' % (config['tmp_dir'], ifo, cluster), 'w')
            out_f.close()
        else:
            # Create unfiltered 
            os.system('ligolw_print -i %s/%s-INSPIRAL_%s.cache -t sngl_inspiral -c end_time -c end_time_ns -c ifo -c snr -c mass1 -c mass2 -c mtotal -c eta -c event_duration -c template_duration -c eff_distance -c chisq -c chisq_dof -c bank_chisq -c bank_chisq_dof  -c cont_chisq -c cont_chisq_dof > %s/%s-0-INSPIRAL_%s_tmp.csv' % (config['tmp_dir'],ifo, cluster, config['tmp_dir'], ifo, cluster))

            # Remove time outside the given range
            in_f  = open('%s/%s-0-INSPIRAL_%s_tmp.csv' % (config['tmp_dir'], ifo, cluster))
            out_f = open('%s/%s-0-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, cluster),'w')

            for l in in_f:
                end_time = int(l.split(',')[0])
                if end_time >= start and end_time < end:
                    print(l, end=' ', file=out_f)
                    print(l, end=' ')
            in_f.close()
            out_f.close()
            os.remove('%s/%s-0-INSPIRAL_%s_tmp.csv' % (config['tmp_dir'], ifo, cluster))

            # and the search_summary
            search_seg = segment(start, end)

            in_f  = open('%s/%s-INSPIRAL_%s.cache' % (config['tmp_dir'], ifo, cluster))
            out_f = open('%s/%s-0-SUMMARY_%s.csv' % (config['tmp_dir'], ifo, cluster),'w')

            for l in in_f:
                start_time, duration = l.split(' ')[2:4]
                start_time = int(start_time)
                end_time   = start_time + int(duration)
                seg        = segment(start_time, end_time)

                # if they don't overlap this will fail
                try:
                    overlap = search_seg & seg
                    print('%d,%d' % (overlap[0], overlap[1]), file=out_f)
                except:
                    pass

            in_f.close()
            out_f.close()

            os.system('ligolw_print -i %s/%s-TMPLTBANK_%s.cache -t sngl_inspiral -c ifo -c mass1 -c mass2 -c mtotal -c eta -c event_duration -c template_duration > %s/%s-TMPLTBANK_%s.csv' % (config['tmp_dir'], ifo, cluster, config['tmp_dir'], ifo, cluster))


        # Give files time to propogate over nfs
        time.sleep(120)

        return

    # Determine times that were vetoed 
    vetoed_times = segmentlist([])

    for c in range(1,category+1):
        if os.path.exists('%s/VETOCATS_%s_%d.p' % (config['tmp_dir'], ifo, c)):
            f = open('%s/VETOCATS_%s_%d.p' % (config['tmp_dir'], ifo, c))
            t = pickle.load(f)
            f.close()

            seglists        = [segmentlist([segment(x[0],x[1]) for x in y[1]]) for y in t]
            vetoed_segments = reduce(operator.or_, seglists).coalesce()
            vetoed_times   |= vetoed_segments

    # Filter triggers
    infile  = open('%s/%s-0-INSPIRAL_%s.csv'  % (config['tmp_dir'], ifo, cluster))
    outfile = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, category, cluster), 'w')

    for l in infile:
        end_time = int(l[0:9])   # assumes that end_time is the first field and has only 9 digits!
        if end_time not in vetoed_times:
            print(l, end=' ', file=outfile)

    infile.close()
    outfile.close()

    # filter search summary
    infile  = open('%s/%s-0-SUMMARY_%s.csv'  % (config['tmp_dir'], ifo, cluster))
    outfile = open('%s/%s-%d-SUMMARY_%s.csv' % (config['tmp_dir'], ifo, category, cluster), 'w')

    lines   = [l.strip().split(',') for l in infile.readlines()]
    summary = segmentlist([segment(int(l[0]), int(l[1])) for l in lines]).coalesce()

    new_summary = summary - vetoed_times
    for l in new_summary:
        print('%d,%d' % (l[0], l[1]), file=outfile)

    infile.close()
    outfile.close()

    # Give files time to propogate over nfs
    time.sleep(120)


def make_database(config, ifo, level, cluster, start, end):
    """
    Reads (ifo)-(level)-INSPIRAL_(cluster).csv and  (ifo)-(level)-TEMPLATE_(cluster).csv, 
    produces (ifo)-(level)-INSPIRAL_(cluster).sqlite with the data loaded and indicies created.
    """

    dbname = '%s/%s-%d-INSPIRAL_%s.sqlite' % (config['tmp_dir'], ifo, level, cluster)
    if os.path.exists(dbname):
        os.remove(dbname)

    # Create the database and sngl_inspiral table
    os.system("sqlite3 %s 'CREATE TABLE sngl_inspiral (end_time int, end_time_ns int, ifo char(2), snr float, mass1 float, mass2 float, mtotal float, eta float, event_duration float, template_duration float, eff_distance float, chisq float, chisq_dof int, bank_chisq float, bank_chisq_dof int, cont_chisq float, cont_chisq_dof int)'" % dbname)

    # populate it
    os.system('sqlite3 -separator , %s  ".import %s/%s-%d-INSPIRAL_%s.csv sngl_inspiral"' % (dbname, config['tmp_dir'], ifo, level, cluster))

    # make it fast
    os.system("sqlite3 %s 'CREATE INDEX sngl_inspiral_snr ON sngl_inspiral(snr)'" % dbname)
    os.system("sqlite3 %s 'CREATE INDEX sngl_inspiral_end_time ON sngl_inspiral(end_time)'" % dbname)
    os.system("sqlite3 %s 'CREATE INDEX sngl_inspiral_mtotal ON sngl_inspiral(mtotal)'" % dbname)


    # Create template table
    os.system("sqlite3 %s 'CREATE TABLE template (ifo char(2), mass1 float, mass2 float, mtotal float, eta float, event_duration float, template_duration float)'" % dbname)
    os.system('sqlite3 -separator , %s  ".import %s/%s-TMPLTBANK_%s.csv template"' %  (dbname, config['tmp_dir'], ifo, cluster))
    os.system("sqlite3 %s 'CREATE INDEX template_mtotal ON template(mtotal)'" % dbname)


    # Create search summary table
    os.system("sqlite3 %s 'CREATE TABLE search_summary(out_start_time int, out_end_time int)'" % dbname)
    os.system('sqlite3 -separator , %s  ".import %s/%s-%d-SUMMARY_%s.csv search_summary"' % (dbname, config['tmp_dir'], ifo, level, cluster))



def setup_db():
    from glue.segmentdb import query_engine

    db_location = os.environ['S6_SEGMENT_SERVER']
    connection  = segmentdb_utils.setup_database(db_location)
    engine      = query_engine.LdbdQueryEngine(connection)

    return engine

##########################
#
# Methods that extract data
#
##########################
def make_summary_table(config, ifo, veto_level, cluster, start_time, end_time):
    # count triggers
    num_triggers = 0
    f = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))
    for l in f:
        num_triggers += 1
    f.close()

    # Find active time
    f     = open('%s/%s-%d-SUMMARY_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))
    lines = [l.strip().split(',') for l in f]
    t     = segmentlist([segment(int(x[0]), int(x[1])) for x in lines]).coalesce()
    f.close()

    analyzed_time = abs(t)

    out = open('%s/%s-%d_%s_summary_table.html' % (config['out_dir'], ifo, veto_level, cluster), 'w')
    print("""<table>
<tr><th>Veto level</th><th>Analyzed Time (s)</th><th>Num. Triggers</th></tr>
<tr><td>%d</td><td>%d</td><td>%d</td></tr>
</table>""" % (veto_level, analyzed_time, num_triggers), file=out)
    out.close()


def make_glitchy_times_table(config, ifo, veto_level, cluster, start_time, end_time):
    """Find intervals with trigger rates > 200 Hz"""

    def tconvert(t):
        return os.popen('lalapps_tconvert -f %H:%M:%S ' + str(t)).readlines()[0].strip()


    out = open('%s/%s-%d_%s_glitchy_times_table.html' % (config['tmp_dir'], ifo, veto_level, cluster), 'w')

    seconds_in_day = range(end_time - start_time)
    counts = [0 for i in seconds_in_day]
    f      = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))

    for line in f:
        counts[int(line.split(',')[0]) - start_time] += 1

    seg_start = -1
    seg_end   = -1


    results = []

    for i in seconds_in_day:
        if counts[i] > 200:
            if seg_start == -1:
                seg_start = i
        else:
            if seg_start != -1:
                seg_end = i

                gps_start = seg_start + start_time
                gps_end   = seg_end   + start_time
                duration  = seg_end - seg_start
                total     = 0

                for i in range(seg_start, seg_end):
                    total += counts[i]

                avg = total / duration

                results.append( (gps_start, gps_end, avg, (gps_end - gps_start)) )

                seg_start = -1
                seg_end   = -1

    print("<h3>Times when trigger rate over 1 s interval exceeds 500 Hz</h3>", file=out)
    print('<table>', file=out)
    print('<tr><th>GPS Start</th><th>GPS End</th><th>Duration (sec)</th><th>Avg. rate (Hz)</th><th>UTC Start</th><th>UTC End</th></tr>', file=out)

    for r in sorted(results, cmp = lambda x,y: cmp(y[2], x[2])):
        if r[2] > 500:
            print('<tr><td>%d</td><td>%d</td><td>%d</td><td>%.0f</td><td>%s</td><td>%s</td></tr>' % (r[0], r[1], r[3], r[2], tconvert(r[0]), tconvert(r[1])), file=out) 

    print('</table>', file=out)


    print("<h3>Times when trigger rate exceedes 200 Hz for more than 10 seconds</h3>", file=out)
    print('<table>', file=out)
    print('<tr><th>GPS Start</th><th>GPS End</th><th>Duration (sec)</th><th>Avg. rate (Hz)</th><th>UTC Start</th><th>UTC End</th></tr>', file=out)

    long = sorted( [r for r in results if r[3] >= 10], cmp = lambda x, y: cmp(y[3], x[3]) )
    for r in long:
        print('<tr><td>%d</td><td>%d</td><td>%d</td><td>%.0f</td><td>%s</td><td>%s</td></tr>' % (r[0], r[1], (r[1] - r[0]), r[2], tconvert(r[0]), tconvert(r[1])), file=out) 

    print('</table>', file=out)

    print('<a href="%s-%d_%s_glitchy_times.txt">All times with rates greater than 200 Hz</a>' % (ifo, veto_level, cluster), file=out)
    out.close()

    out = open('%s/%s-%d_%s_glitchy_times.txt' % (config['out_dir'], ifo, veto_level, cluster), 'w')
    print("# start time, end time, avg rate, duration", file=out)
    for r in results:
        print("%d\t%d\t%.0f\t%d" % r, file=out)
    out.close()



def make_usage_table(config, ifo, veto_level, cluster, start_time, end_time):
    out = open('%s/%s-%d_%s_usage_table.html' % (config['out_dir'], ifo, veto_level, cluster), 'w')

    if veto_level == 0:
        out.close()
        return

    # Open the database at the previous level
    if veto_level == 4:
        connection = setup_memory_db(config, ifo, 2, cluster)
    else:
        connection = setup_memory_db(config, ifo, veto_level-1, cluster)

    cursor     = connection.cursor()

    # and the veto information at this level
    f      = open('%s/VETOCATS_%s_%d.p' % (config['tmp_dir'], ifo, veto_level))
    vetoes = pickle.load(f)
    f.close()

    total_triggers = cursor.execute('SELECT COUNT(*) FROM sngl_inspiral').fetchone()[0]
    analysis_time  = segmentlist([segment(row[0], row[1]) for row in cursor.execute('SELECT out_start_time, out_end_time FROM search_summary')]).coalesce()
    seglists       = [(y[0], segmentlist([segment(x[0],x[1]) for x in y[1]])) for y in vetoes]

    if veto_level == 1:
        print("<h3>Efficiency of category 1 vetoes</h3>", file=out)
    else:
        print("<h3>Efficiency of category %d vetoes on triggers that passed category %d</h3>" % (veto_level, veto_level - 1), file=out)

    print(file=out)

    if total_triggers == 0:
        print("Note: no triggers at this veto level<p>", file=out)

    atime = float(abs(analysis_time))

    if atime == 0.0:
        print("Note: no analysis time at this veto level<p>", file=out)

    print("<table>", file=out)
    print("<tr><th>Name:Version</th><th>Efficiency (%)</th><th>Deadtime (%)</th><th>Efficiency / Deadtime</th></tr>", file=out)

    for seg in seglists:
        if len(seg[1]) == 0:
            clause = '1 = 1'
        else:
            clause   = 'NOT (' + ' OR '.join(['(end_time BETWEEN %d AND %d)' % (s[0], s[1]) for s in seg[1]]) + ')'

        triggers = cursor.execute('SELECT COUNT(*) FROM sngl_inspiral WHERE ' + clause).fetchone()[0]
        active_time = analysis_time - seg[1] 

        if total_triggers == 0:
            efficiency = 'NA'
        else:
            eff = 100.0 * float(total_triggers - triggers) / float(total_triggers)
            efficiency = '%.2f' % eff

        if atime == 0.0:
            dead_time = 'NA'
        else:
            dtime  = 100.0 * (atime - float(abs(active_time))) / atime
            dead_time = '%.2f' % dtime

        if dead_time != 'NA' and efficiency != 'NA' and dtime > 0:
            ratio = '%.2f' % (eff / dtime)
        else:
            ratio = 'NA'

        print("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (seg[0], efficiency, dead_time, ratio), file=out)

    print("</table>", file=out)

    out.close()



def make_glitch_page(config, ifo, veto_level, cluster, gps_start_time, gps_end_time):
    """Makes a glitch table """

    out = open('%s/%s_CAT%d_%s_Glitch.html' % (config['out_dir'], ifo, veto_level, cluster), 'w')

    print('<html>', file=out)
    print('<head>', file=out)
    print('  <link media="all" href="../../auxfiles/ihope_daily_style.css" type="text/css" rel="stylesheet" />', file=out)
    print('</head>', file=out)

    print('<body>', file=out)

    if cluster != '16SEC_CLUSTERED':
        print("<p><p><H3>The glitch page is only available for triggers that have been clustered with a 16-second window</H3>", file=out)
    else:
        tmpfile = '%s/%s-%d_glitchout.html' % (config['tmp_dir'], ifo, veto_level)

        cmd  = '%s ' % config['ligolw_cbc_glitch_page']
        cmd += '--statistic newsnr --min-glitch-snr 7.0 '
        cmd += '--trigger-db %s-%d-INSPIRAL_16SEC_CLUSTERED.sqlite '  % (ifo, veto_level)
        cmd += '--segments %s '     % os.environ['S6_SEGMENT_SERVER']
        cmd += '--html-file %s ' % tmpfile
        cmd += '--ifo %s ' % ifo
        cmd += '--timestamp-file none --known-count 5 --unknown-count 5 '
        cmd += '--gps-start-time %d '  % gps_start_time
        cmd += '--gps-end-time %d '    % gps_end_time
        cmd += '--omega-conf %s '      % config['omega_conf']
        cmd += '--omega-frame-dir %s ' % config['omega_frame_dir']
        cmd += '--omega-out-dir %s/omega ' % config['out_dir']

        os.popen(cmd)

        glitch_text = "".join(open(tmpfile).readlines()[1:-2])

        os.remove(tmpfile)

        print(glitch_text, file=out)

    print('</body>', file=out)
    print('</html>', file=out)
    out.close()



def make_hwinj_page(config, ifo, veto_level, cluster, gps_start_time, gps_end_time):
    """Makes a hardware injection table """

    #out = open('%s_CAT%d_%s_Hwinj.html' % (ifo, veto_level, cluster), 'w')
    #print >>out, '<html>'
    #print >>out, '<head>'
    #print >>out, '  <link media="all" href="../../auxfiles/ihope_daily_style.css" type="text/css" rel="stylesheet" />'
    #print >>out, '</head>'

    #print >>out, '<body>'

    ifo_args = {'H1':'-i ', 'L1':'-l ', 'V1':'-v '}

    anydata = True

    if cluster != '16SEC_CLUSTERED':
        out = open('%s/%s_CAT%d_%s_Hwinj.html' % (config['out_dir'], ifo, veto_level, cluster), 'w')
        print('<html>', file=out)
        print('<head>', file=out)
        print('  <link media="all" href="ihope_daily_style.css" type="text/css" rel="stylesheet" />', file=out)
        print('</head>', file=out)

        print('<body>', file=out)
        print("<p><p><H3>The hardware injection page is only available for triggers that have been clustered with a 16-second window</H3>", file=out)
        print('</body>', file=out)
        print('</html>', file=out)
        out.close()
    else:
        cmd  = '%s ' % config['ligolw_cbc_hardware_inj_page']
        cmd += '-t %d -e %d ' % (gps_start_time, gps_end_time)
        cmd += ifo_args[ifo]
        cmd += '-x %s ' % config['hwinj_file']
        cmd += '-s %s ' % os.environ['S6_SEGMENT_SERVER']
        cmd += '-o %s/%s_CAT%d_%s_Hwinj.html ' % (config['out_dir'], ifo, veto_level, cluster)
        cmd += ' '.join(sorted(glob.glob('%s-INSPIRAL_16SEC_CLUSTERED-*xml' % ifo)))

        os.system(cmd)


def make_index_page(config, ifo, veto_level, cluster, gps_start_time, gps_end_time):
    map = {}

    summary_names = ['%s-%d_%s_summary_table.html','%s-%d_%s_usage_table.html']
    section_names = ['%s-%d_%s_glitch_page.html',  '%s-%d_%s_hwinj_page.html']
    image_names   = ['rate_vs_time','snr_hist', 'glitchgram']

    map['start_time'] = gps_start_time
    map['end_time']   = gps_end_time
    map['date']       = os.popen('lalapps_tconvert -f "%B %d, %Y" ' + str(gps_start_time)).readlines()[0].strip()
    template          = html_template % map

    out = open('%s/index.html' % config['out_dir'],'w')
    print(template, file=out)
    out.close()

    for ifo in ifos:
        for cat in [0,1,2,4]:
            for cluster in ['UNCLUSTERED','30MILLISEC_CLUSTERED','16SEC_CLUSTERED']:
                # Make the summary section

                out = open('%s/%s_CAT%d_%s_Summary.html' % (config['out_dir'], ifo, cat, cluster), 'w')
                print('<html>', file=out)
                print('<head>', file=out)
                print('  <link media="all" href="ihope_daily_style.css" type="text/css" rel="stylesheet" />', file=out)
                print('</head>', file=out)

                print('<body>', file=out)

                for name in summary_names:
                    in_f = open(config['out_dir'] + '/' + name % (ifo, cat, cluster))

                    for l in in_f:
                        print(l, file=out)

                    in_f.close()

                print('</body>', file=out)
                print('</html>', file=out)
                out.close()

                # Make the bank chisq page
                out = open('%s/%s_CAT%d_%s_bankchisq.html' % (config['out_dir'], ifo, cat, cluster), 'w')
                print('<html>', file=out)
                print('<head>', file=out)
                print('  <link media="all" href="ihope_daily_style.css" type="text/css" rel="stylesheet" />', file=out)
                print('</head>', file=out)
                print('<body>', file=out)

                #imgs = sorted(glob.glob('%s/%s_%d_%s_bank_veto_dof_*.png' % (config['out_dir'], ifo, cat, cluster)))
                #for i in imgs:
                #    print >>out, '<img src="%s"><p>' % i.split('/')[-1]
                print('<img src="%s_%d_%s_chisq.png"><p>' % (ifo, cat, cluster), file=out)
                print('<img src="%s_%d_%s_bank_veto.png"><p>' % (ifo, cat, cluster), file=out)
                print('<img src="%s_%d_%s_cont_veto.png"><p>' % (ifo, cat, cluster), file=out)

                print('</body>', file=out)
                print('</html>', file=out)
                out.close()

                # Add the images
                for name in image_names:
                    out = open('%s/%s_CAT%d_%s_%s.html' % (config['out_dir'], ifo, cat, cluster, name), 'w')

                    print('<html>', file=out)
                    print('<head>', file=out)
                    print('  <link media="all" href="ihope_daily_style.css" type="text/css" rel="stylesheet" />', file=out)
                    print('</head>', file=out)
                    print('<body>', file=out)

                    if name == 'glitchgram':
                        print('<p><img src="%s_%d_%s_new_glitchgram.png"><p>' % (ifo, cat, cluster), file=out)

                    if name == 'snr_hist':
                        print('<p><img src="%s_%d_%s_new_snr_hist.png"><p>' % (ifo, cat, cluster), file=out)
                        print('<p><img src="%s_%d_%s_%s_all.png">' % (ifo, cat, cluster, name), file=out)

                    print('<img src="%s_%d_%s_%s.png">' % (ifo, cat, cluster, name), file=out)

                    if name == 'rate_vs_time':
                        print('<img src="%s_%d_%s_newsnr_vs_time.png"><p>' % (ifo, cat, cluster), file=out)
                        print('<img src="%s_%d_%s_snr_vs_time.png">' % (ifo, cat, cluster), file=out)

                        in_f = open('%s/%s-%d_%s_glitchy_times_table.html' % (config['tmp_dir'], ifo, cat, cluster))
                        for l in in_f:
                            print(l, end=' ', file=out) 
                        in_f.close()

                    print('</body>', file=out)
                    print('</html>', file=out)
                    out.close()

                # And the template page
                out = open('%s/%s_CAT%d_%s_template.html' % (config['out_dir'], ifo, cat, cluster), 'w')
                print('<html>', file=out)
                print('<head>', file=out)
                print('  <link media="all" href="ihope_daily_style.css" type="text/css" rel="stylesheet" />', file=out)
                print('</head>', file=out)

                print('<body>', file=out)

                print('<img src="%s_%d_%s_mass_hist.png">' % (ifo, cat, cluster), file=out)
                print('<img src="%s_%d_%s_tmpl_hist.png">' % (ifo, cat, cluster), file=out)
                print('<img src="%s_%d_%s_mass_hist_norm.png">' % (ifo, cat, cluster), file=out)
                print('<img src="%s_%d_%s_template_counts.png">' % (ifo, cat, cluster), file=out)
                print('<img src="%s_%d_%s_hexmass.png">' % (ifo, cat, cluster), file=out)
                print('</body>', file=out)
                print('</html>', file=out)
                out.close()

    # Copy asset files
    shutil.copy('%s/ihopeFrog.jpg'         % config['asset_dir'], config['out_dir'])
    shutil.copy('%s/ihope_daily_toggle.js' % config['asset_dir'], config['out_dir'])
    shutil.copy('%s/ihope_daily_style.css' % config['asset_dir'], config['out_dir'])


#########################
#
# Methods that make plots
#
#########################
def plot_snr_vs_time(config, ifo, veto_level, cluster, start_time, end_time):
    """snr vs. time"""

    SECS_PER_HOUR = 60.0 * 60.0
    start_time2   = float(start_time)

    snr    = []
    newsnr = []
    gps    = []
    gpsns  = []
    time   = []
    infTimes = []

    f = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))
    for line in f:
        tmp = line.strip().split(',')
        gps.append(float(tmp[0]))
        gpsns.append(float(tmp[1]))
        if tmp[3] == 'inf':
            tmp[3] = 5000
            infTimes.append(float(tmp[0]))
        snr.append(float(tmp[3]))
        newsnr.append(new_snr(float(tmp[3]), float(tmp[11]), float(tmp[12])))
    time = (pylab.array(gps) - start_time) / SECS_PER_HOUR

    fig = pylab.figure(0, figsize=(8.1,4.0))
    axes = pylab.Axes(fig,[0.07, 0.1, 1.0 - 0.07 - 0.03, 0.8])
    fig.add_axes(axes)

    if len(snr):
        pylab.plot(time, snr, ifo_colored_circle[ifo], markersize = 2)
        pylab.grid()
        pylab.xlim(0, 24)

        if max(snr) > 11:
            pylab.yscale('log')
        pylab.ylim(5.5, max(snr)*1.1)
        maxindex = pylab.argmax(snr)
        pylab.text(time[maxindex], snr[maxindex]*1.04, 'Max: GPS %.3f' % (gps[maxindex]+1e-9*gpsns[maxindex]), size='medium')
        if len(infTimes):
            pylab.text(0.5, max(snr)*1.07, 'One or more SNRs were inf!! See .err file for a list of times')
            print('Times with SNR=inf:', infTimes, file=sys.stderr)

    pylab.title(make_timeless_title(ifo, veto_level, cluster, start_time))
    pylab.xlabel(make_time_label(start_time))
    pylab.xticks(range(0,24,2))
    pylab.ylabel('SNR')
    pylab.savefig('%s/%s_%d_%s_snr_vs_time.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()

    fig  = pylab.figure(1, figsize=(8.1,4.0))
    axes = pylab.Axes(fig,[0.07, 0.1, 1.0 - 0.07 - 0.03, 0.8])
    fig.add_axes(axes)

    if len(snr):
        pylab.plot(time, newsnr, ifo_colored_circle[ifo], markersize = 2)
        pylab.grid()
        pylab.xlim(0, 24)

        maxnew = max(newsnr)
        newsnr = pylab.array(newsnr)
        is_big = newsnr>0.95*maxnew
        pylab.plot(time[is_big], newsnr[is_big], 'kx', markersize=8)
        maxindex = pylab.argmax(newsnr)
        pylab.text(time[maxindex], newsnr[maxindex]*1.02, 'Max: GPS %.3f' % (gps[maxindex]+1e-9*gpsns[maxindex]), size='medium')
        pylab.ylim(6, maxnew*1.05)

    pylab.title(make_timeless_title(ifo, veto_level, cluster, start_time))
    pylab.xlabel(make_time_label(start_time))
    pylab.xticks(range(0,24,2))
    pylab.ylabel('NewSNR')
    pylab.savefig('%s/%s_%d_%s_newsnr_vs_time.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()

def wait_for(filename):
    # wait for it to exist
    count = 10

    while count > 0 and not os.path.exists(filename):
        time.sleep(30)
        count -= 1

    if not os.path.exists(filename):
        print("Needed file %s not found, aborting" % filename, file=sys.stderr)
        sys.exit(1)

    # wait for it to stabalize
    last_size = 0
    size      = os.path.getsize(filename)
    count     = 5 

    while count > 0 and (size != last_size or size == 0):
        time.sleep(30)
        count    -= 1
        last_size = size
        size      = os.path.getsize(filename)

    from random import uniform
    num = int(uniform(0,1000))
    os.system('cat %s > tmp_%d_%s' % (filename, num, filename))


def setup_memory_db(config, ifo, level, cluster):
    wait_for('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, level, cluster))
    conn = sqlite3.connect(':memory:')
    cursor = conn.cursor()

    cursor.execute('CREATE TABLE sngl_inspiral (end_time int, end_time_ns int, ifo char(2), snr float, mass1 float, mass2 float, mtotal float, eta float, event_duration float, template_duration float, eff_distance float, chisq float, chisq_dof int, bank_chisq float, bank_chisq_dof int, cont_chisq float, cont_chisq_dof int)')

    f = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, level, cluster))
    totypes = [int, int, str, float, float, float, float, float, float, float, float, float, int, float, int, float, int]

    for l in f:
        vals = l.strip().split(',')
        typed_vals = tuple( map(lambda x,y: x(y), totypes, vals) )

        cursor.execute('INSERT INTO sngl_inspiral VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', typed_vals)

    # make it fast
    cursor.execute('CREATE INDEX sngl_inspiral_snr ON sngl_inspiral(snr)')
    cursor.execute('CREATE INDEX sngl_inspiral_end_time ON sngl_inspiral(end_time)')
    cursor.execute('CREATE INDEX sngl_inspiral_mtotal ON sngl_inspiral(mtotal)')


    # Create template table
    cursor.execute('CREATE TABLE template (ifo char(2), mass1 float, mass2 float, mtotal float, eta float, event_duration float, template_duration float)')

    f = open('%s/%s-TMPLTBANK_%s.csv' % (config['tmp_dir'], ifo, cluster))
    totypes = [str, float, float, float, float, float, float]

    for l in f:
        vals = l.strip().split(',')
        typed_vals = tuple( map(lambda x,y: x(y), totypes, vals) )

        cursor.execute('INSERT INTO template VALUES(?,?,?,?,?,?,?)', typed_vals)

    cursor.execute('CREATE INDEX template_mtotal ON template(mtotal)')


    # Create search summary table
    cursor.execute('CREATE TABLE search_summary(out_start_time int, out_end_time int)')

    f = open('%s/%s-%d-SUMMARY_%s.csv' % (config['tmp_dir'], ifo, level, cluster))
    totypes = [int, int]

    for l in f:
        vals = l.strip().split(',')
        typed_vals = tuple( map(lambda x,y: x(y), totypes, vals) )

        cursor.execute('INSERT INTO search_summary VALUES(?,?)', typed_vals)

    return conn


def plot_rate_vs_time(config, ifo, veto_level, cluster, start_time, end_time):
    """Rate (per sec) vs. time"""

    num_min     = (end_time - start_time) / 60
    mins_in_day = range(num_min)
    counts = [1e-6 for i in mins_in_day]
    f      = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))

    for line in f:
        trig_end_time = int(line.split(',')[0]) 
        counts[ (trig_end_time - start_time) / 60 ] += 1

    fig  = pylab.figure(0, figsize=(8.1,4.0))
    axes = pylab.Axes(fig,[0.08, 0.1, 1.0 - 0.07 - 0.03, 0.8])
    fig.add_axes(axes)

    times  = [float(x) / 60.0 for x in mins_in_day]
    counts = [c / 60.0 for c in counts]

    pylab.semilogy(times, counts, ifo_colored_circle[ifo], markersize = 2)
    pylab.grid()

    pylab.xlim(0, 24)
    min_y = cluster.find('TSSI') > -1 and 1e-2 or 1
    pylab.ylim(min_y, max(counts) * 1.1)
    pylab.xticks(range(0,24,2))

    pylab.title(make_timeless_title(ifo, veto_level, cluster, start_time))
    pylab.xlabel(make_time_label(start_time))
    pylab.ylabel('Rate (Hz)')
    pylab.savefig('%s/%s_%d_%s_rate_vs_time.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()


def plot_bank_veto(config, ifo, veto_level, cluster, start_time, end_time):
    """Bank chi^2 vs snr"""

    f = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))

    snrs   = []
    bank_chisqs = []
    cont_chisqs = []
    chisqs      = []

    for l in f:
        data  = l.split(',')
        rho   = float(data[3])

        chisqs.append(float(data[11]))

        bank_chisq = float(data[13])
        bank_dof   = int(data[14])

        # Now using continuous chisq
        cont_chisq = float(data[15])
        cont_dof   = int(data[16])

        if cont_dof == 0:
            cont_dof = 1

        if bank_dof == 0:
            bank_dof = 1.5

        snrs.append(rho)
        bank_chisqs.append(bank_chisq / (2 * bank_dof - 2))
        cont_chisqs.append(cont_chisq / cont_dof)

    pylab.figure(0)
    pylab.yscale('log')

    if snrs and max(bank_chisqs)>0:
        pylab.plot(snrs,bank_chisqs,'rx')

        if max(snrs) > 11:
            pylab.xscale('log')

        pylab.xlim(5.0,max(snrs) * 1.1)
        pylab.ylim(min(bank_chisqs) * 0.9, max(bank_chisqs) * 1.1)

    pylab.grid()
    pylab.xlabel('snr')
    pylab.ylabel(r'bank\_chisq / (2p-2)')
    pylab.title(make_title(ifo, veto_level, cluster, start_time))

    pylab.savefig('%s/%s_%d_%s_bank_veto.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()

    pylab.figure(1)
    pylab.yscale('log')

    if snrs and max(cont_chisqs)>0:
        pylab.plot(snrs,cont_chisqs,'rx')

        if max(snrs) > 11:
            pylab.xscale('log')

        pylab.xlim(5.0,max(snrs) * 1.1)
        pylab.ylim(min(cont_chisqs) * 0.9, max(cont_chisqs) * 1.1)

    pylab.grid()
    pylab.xlabel('snr')
    pylab.ylabel('cont. chisq / dof')
    pylab.title(make_title(ifo, veto_level, cluster, start_time))

    pylab.savefig('%s/%s_%d_%s_cont_veto.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()

    pylab.figure(2)
    pylab.yscale('log')

    if snrs and max(chisqs)>0:
        pylab.plot(snrs,chisqs,'rx')

        if max(snrs) > 11:
            pylab.xscale('log')

        pylab.xlim(5.0,max(snrs) * 1.1)
        pylab.ylim(min(chisqs) * 0.9, max(chisqs) * 1.1)

    pylab.grid()
    pylab.xlabel('snr')
    pylab.ylabel('chisq')
    pylab.title(make_title(ifo, veto_level, cluster, start_time))

    pylab.savefig('%s/%s_%d_%s_chisq.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()


def plot_bank_veto_old(config, ifo, veto_level, cluster, start_time, end_time):
    """Bank chi^2 vs snr: Old version, made each DOF a separate plot"""

    f = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))

    snrs   = []
    chisqs = []

    for l in f:
        data = l.split(',')
        dof  = int(data[14])
        rho  = float(data[3])
        chisq = float(data[13])

        while len(snrs) < dof+1:
            snrs.append([])
            chisqs.append([])

        snrs[dof].append(rho)
        chisqs[dof].append(chisq)

    imgcount = 0

    for i in range(len(snrs)):
        if snrs[i]:
            pylab.figure(imgcount)
            imgcount += 1

            pylab.plot(snrs[i],chisqs[i],'rx')
            pylab.xscale('log')
            pylab.yscale('log')
            pylab.xlim(5.0,max(snrs[i]) * 1.1)
            pylab.ylim(min(chisqs[i]) * 0.9, max(chisqs[i]) * 1.1)
            pylab.grid()
            pylab.xlabel('snr')
            pylab.ylabel('chisq')
            pylab.title(make_title(ifo, veto_level, cluster, start_time, desc='dof = %d' % i))

            pylab.savefig('%s/%s_%d_%s_bank_veto_dof_%d.png' % (config['out_dir'], ifo, veto_level, cluster, i))
            pylab.close()

def plot_snr_hist(config, ifo, veto_level, cluster, start_time, end_time):
    """Histogram of SNRs"""
    f = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))

    remaining = 0
    huge_snrs = 0
    upper_lim = 1e7
    total_count = 0
    bins        = [1e-40 for i in range(int(upper_lim))]  # Start off non-empty to avoid problems with no triggers
    biggest     = 0

    for line in f:
        tmp    = line.strip().split(',')
        snrsq  = float(tmp[3])
        snrsq  = int(snrsq * snrsq)

        total_count += 1

        if snrsq > 200:
            remaining += 1

        if snrsq > upper_lim:
            huge_snrs += 1
        else:
            if snrsq > biggest:
                biggest = snrsq
            bins[snrsq] += 1

    if total_count == 0:
        return

    total_count = float(total_count)

    gaussian = None

    if len(bins) > 31 and cluster != '16SEC_CLUSTERED':
        slope    = -0.5  # 1.0 - pylab.sqrt(2)
        offset   = float(bins[30]) / pylab.exp(slope * 30.0)
        gaussian = [offset * pylab.exp(slope * x) for x in range(biggest)]
        gtotal   = 0

        for i in range(31,biggest):
            gtotal += bins[i] - gaussian[i]

        non_gaussianity = gtotal / total_count

        f_out = open('%s/%s-%d-%s_nongaussianity.txt' % (config['out_dir'], ifo, veto_level, cluster), 'w')
        print(non_gaussianity, file=f_out)
        f_out.close()

        # Plot up to 200, with a cumulative dot showing the remaining
        fig = pylab.figure(0, figsize=(8.1,4.0))
        axes = pylab.Axes(fig,[0.07, 0.1, 1.0 - 0.07 - 0.03, 0.8])
        fig.add_axes(axes)

        high_end = min(biggest, 200)

        pylab.semilogy(range(30,high_end), bins[30:high_end], color = ifo_colors[ifo], label='_nolegend_')
        pylab.plot([199.5], [remaining], 'ko', label='Count > 200')

        if gaussian:
            pylab.semilogy(range(30,high_end), gaussian[30:high_end], 'k--', label='_nolegend_')

        pylab.grid()
        pylab.legend()
        pylab.ylim(0.5,max(max(bins[30:high_end]), remaining) * 1.1)
        pylab.xlim(30, high_end)
        pylab.title(make_title(ifo, veto_level, cluster, start_time))
        pylab.xlabel('SNR^2')
        pylab.ylabel('Count')
        pylab.savefig('%s/%s_%d_%s_snr_hist.png' % (config['out_dir'], ifo, veto_level, cluster))
        pylab.close()

        # Plot everything!
        fig = pylab.figure(1, figsize=(8.1,4.3))
        axes = pylab.Axes(fig,[0.07, 0.12, 1.0 - 0.07 - 0.03, 0.8])
        fig.add_axes(axes)

        pylab.loglog(range(30,biggest), bins[30:biggest], color = ifo_colors[ifo])

        if gaussian:
            pylab.loglog(range(30,biggest), gaussian[30:biggest], 'k--', label='_nolegend_')

        pylab.grid()

        if huge_snrs > 0:
            pylab.loglog([upper_lim - 100], [huge_snrs], 'ko', label='Count > %.0e' % upper_lim)
            pylab.legend()

        pylab.ylim(0.5,max(max(bins), huge_snrs) * 1.1)
        pylab.xlim(30, len(bins))
        pylab.title(make_title(ifo, veto_level, cluster, start_time))
        pylab.xlabel('SNR^2')
        pylab.ylabel('Count')
        pylab.savefig('%s/%s_%d_%s_snr_hist_all.png' % (config['out_dir'], ifo, veto_level, cluster))
        pylab.close()
    #
    # Now plot the new snr histograms
    #
    f = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))
    bins = [1e-40 for i in range(0,5)]  # Start off non-empty to avoid problems with no triggers
    remaining = 0
    huge_snrs = 0
    upper_lim = 100
    total_count = 0

    for line in f:
        tmp    = line.strip().split(',')
        snr    = float(tmp[3])
        chisq  = float(tmp[11])
        dof    = float(tmp[12])
        newsnr = int(10 * new_snr(snr, chisq, dof))

        total_count += 1

        while len(bins) <= newsnr:
            try:
                bins.append(1e-40)
            except:
                print("Unable to allocate %f" % newsnr)

        bins[newsnr] += 1

    total_count = float(total_count)

    fig  = pylab.figure(2, figsize=(8.1,4.0))
    axes = pylab.Axes(fig,[0.07, 0.1, 1.0 - 0.07 - 0.03, 0.8])
    fig.add_axes(axes)

    xs = [x / 10.0 for x in range(50, len(bins))]
    pylab.semilogy(xs, bins[50:], color = ifo_colors[ifo])

    pylab.grid()

    if len(bins) > 5:
        pylab.ylim(0.5,max(bins[5:]) * 1.1)
        pylab.xlim(5,max(xs) * 1.1)
    pylab.title(make_title(ifo, veto_level, cluster, start_time))
    pylab.xlabel('New SNR')
    pylab.ylabel('Count')
    pylab.savefig('%s/%s_%d_%s_new_snr_hist.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()


def plot_template_counts(config, ifo, veto_level, cluster, start_time, end_time):
    """Template counts as a 2-d scatter plot"""

    connection = setup_memory_db(config, ifo, veto_level, cluster)

    res    = [row for row in connection.cursor().execute('SELECT mtotal, eta, COUNT(*) FROM sngl_inspiral GROUP BY mtotal, eta')]
    mtots  = [float(x[0]) for x in res]
    etas   = [float(x[1]) for x in res]

    #res    = [row for row in connection.cursor().execute('SELECT mass1, mass2, COUNT(*) FROM sngl_inspiral GROUP BY mass1, mass2')]
    #mtots  = [float(row[0]) + float(row[1]) for row in res]
    #etas   = [ (float(row[0]) * float(row[1])) / pow(float(row[0]) + float(row[1]), 2) for row in res]
    counts = [pylab.log10(float(x[2]) + 0.001) for x in res]

    pylab.figure(0)

    if len(mtots) > 0:
        pylab.scatter(mtots, etas, c=counts, cmap=pylab.cm.spectral, edgecolor='none', alpha=0.8)
    else:
        pylab.scatter([0], [0], c=[0], cmap=pylab.cm.hot)

    pylab.xlim(1,36)
    pylab.ylim(0.001, 0.26)

    pylab.colorbar()
    pylab.title('Log(count) of triggers per template')
    pylab.xlabel('Total Mass (M_sun)')
    pylab.ylabel('eta')
    pylab.savefig('%s/%s_%d_%s_template_counts.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()

    connection.close()


def plot_hexmass(config, ifo, veto_level, cluster, start_time, end_time):
    """Hexbined plot of total mass vs. time"""

    SECS_PER_HOUR = 60.0 * 60.0
    start_time2   = float(start_time)

    mtots = []
    times = []

    f = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))
    for line in f:
        tmp = line.strip().split(',')
        times.append( (float(tmp[0]) - start_time) / SECS_PER_HOUR )
        mtots.append( float(tmp[6]) )

    num_map = {'H1':0, 'L1':1, 'V1':2}
    pylab.figure(num_map[ifo] * 3 + veto_level, figsize=(8.1,4))

    if len(times) > 0:
        pylab.hexbin(times, mtots)
        pylab.colorbar()

    pylab.title(make_timeless_title(ifo, veto_level, cluster, start_time, 'Trigger count'))
    pylab.xlabel(make_time_label(start_time))
    pylab.xticks(range(0,24,2))
    pylab.ylabel('Mass (M_sun)')
    pylab.savefig('%s/%s_%d_%s_hexmass.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()


def plot_mass_hist(config, ifo, veto_level, cluster, start_time, end_time, flower):
    """Histogram of templates by mass"""

    def filename_to_bins(filename, s, e, d):
        f = open(filename)
        mtots_and_etas = [l.strip().split(',')[s:e] for l in f.readlines()]
        f.close()
        f = open(filename)
        durations      = [chirplen(flower, float(x[0]), float(x[1])) for x in mtots_and_etas]
        trig_durations = [float(l.strip().split(',')[d]) for l in f.readlines()]
        f.close()

        maxduration = max([max(trig_durations),max(durations)])
        rnge = range(int(maxduration) + 5)

        bins = [0 for x in rnge]
        chirpTbins = [0 for x in rnge]
        for d,d2 in zip(trig_durations,durations):
            if d2 < d or d < 0.0001:
                d = d2
            bins[int(d)] += 1
            chirpTbins[int(d2)] += 1

        return rnge,bins,chirpTbins

    rnge,bins,inspChirpBins = filename_to_bins('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster), 6, 8, 9)
    rnge_insp = rnge
    fig  = pylab.figure(0, figsize=(8.1,4.0))
    # axes = pylab.Axes(fig,[0.07, 0.1, 1.0 - 0.07 - 0.03, 0.8])
    # fig.add_axes(axes)

    #xlocations = na.array(range(len(bins)))+0.8
    #width      = 0.8
    #pylab.bar(xlocations, bins, width=width)
    #pylab.xticks(xlocations+ width/2, rnge)

    #pylab.gca().get_xaxis().tick_bottom()
    #pylab.gca().get_yaxis().tick_left()

    if max(bins) > 1000:
        pylab.semilogy(rnge, bins, color=ifo_colors[ifo])
    else:
        pylab.plot(rnge, bins, color=ifo_colors[ifo])

    pylab.grid()
    pylab.title(make_title(ifo, veto_level, cluster, start_time))
    pylab.legend(['Triggers by duration'])
    pylab.xlabel('Duration (sec)')
    pylab.ylabel('Num. triggers')
    pylab.savefig('%s/%s_%d_%s_mass_hist.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()

    rnge, weights, tmpltChirpWeights = filename_to_bins('%s/%s-TMPLTBANK_%s.csv' % (config['tmp_dir'], ifo, cluster), 3, 5, 6)

    pylab.figure(2, figsize=(8.1,4.0))

    if max(weights) > 0:
        pylab.semilogy(rnge, weights, color=ifo_colors[ifo])
        pylab.grid()
    else:
        pylab.plot([-2,-1],[-2,-1], color=ifo_colors[ifo])
        pylab.xlim(0,45)
        pylab.ylim(0,1)
        pylab.grid()

    pylab.title(make_title(ifo, veto_level, cluster, start_time))
    pylab.legend(['Templates by chirplen'])
    pylab.xlabel('Chirp length from %fHz (sec)' %(flower))
    pylab.ylabel('Num. Templates')
    pylab.savefig('%s/%s_%d_%s_tmpl_hist.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()

    tmpltChirpWeights = [w or 1 for w in tmpltChirpWeights]
    if len(rnge_insp) < len(rnge):
      rnge = rnge_insp
    wbins   = [inspChirpBins[i] / tmpltChirpWeights[i] for i in rnge]

    pylab.figure(1, figsize=(8.1,4))

    pylab.plot(rnge, wbins, color=ifo_colors[ifo])
    pylab.grid()
    pylab.title(make_title(ifo, veto_level, cluster, start_time))
    pylab.legend(['Avg. triggers per template'])
    pylab.xlabel('Chirp length from %fHz (sec)' %(flower))
    pylab.ylabel('Num. triggers / Num. templates')

    pylab.savefig('%s/%s_%d_%s_mass_hist_norm.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()

def plot_glitchgram(config, ifo, veto_level, cluster, start_time, end_time):
    f      = open('%s/%s-%d-INSPIRAL_%s.csv' % (config['tmp_dir'], ifo, veto_level, cluster))
    gps    = []
    gpsns  = []
    mtot   = []
    eta    = []
    snr    = []
    newsnr = []
    duration = []
    for line in f:
      d = line.split(',')
      gps.append(float(d[0]))
      gpsns.append(float(d[1]))
      mtot.append(float(d[6]))
      eta.append(float(d[7]))
      duration.append(float(d[9]))
      snr.append(float(d[3]))
      newsnr.append( new_snr(float(d[3]), float(d[11]), float(d[12])) )

    gps    = pylab.array(gps)
    gpsns  = pylab.array(gpsns)
    mtot   = pylab.array(mtot)
    eta    = pylab.array(eta)
    snr    = pylab.array(snr)
    newsnr = pylab.array(newsnr)
    seconds_in_hour = 60.0 * 60.0
    time = (gps-start_time) / seconds_in_hour
    duration = pylab.array(duration)

    if not len(gps):
        # If there's no data generate an empty plot
        fig  = pylab.figure(1, figsize=(8.1,4.0))
        axes = pylab.Axes(fig,[0.08, 0.1, 1.0 - 0.07 - 0.03, 0.8])
        fig.add_axes(axes)
        pylab.xlim(0, 24)
        pylab.xticks(range(0,24,2))
        pylab.ylabel('Template duration (sec)')
        pylab.title('No triggers')
        pylab.savefig('%s/%s_%d_%s_glitchgram.png' % (config['out_dir'], ifo, veto_level, cluster))
        pylab.close()

        fig  = pylab.figure(2, figsize=(8.1,4.0))
        axes = pylab.Axes(fig,[0.08, 0.1, 1.0 - 0.07 - 0.03, 0.8])
        fig.add_axes(axes)
        pylab.xlim(0, 24)
        pylab.xticks(range(0,24,2))
        pylab.ylabel('Template duration (sec)')
        pylab.title('No triggers')
        pylab.savefig('%s/%s_%d_%s_new_glitchgram.png' % (config['out_dir'], ifo, veto_level, cluster))
        pylab.close()

        return

    ranges = [[ 0, 8,'b', 3],
              [ 8,16,'g',15],
              [16,float('inf'),'r',15]]

    fig  = pylab.figure(1, figsize=(8.1,4.0))
    axes = pylab.Axes(fig,[0.08, 0.1, 1.0 - 0.07 - 0.03, 0.8])
    fig.add_axes(axes)

    for r in ranges:
        subset = pylab.logical_and(snr>=r[0], snr<r[1])
        print('Plotting trigs for snr between', r[0], r[1])
        if len(time[subset]):
            pylab.scatter(time[subset], duration[subset], edgecolor='none', c=r[2], s=r[3])
    snrmax_index = pylab.argmax(snr)
    pylab.scatter(time[snrmax_index], duration[snrmax_index], marker=(5,1,0), c='y', s=100)

    pylab.yscale('log')
    pylab.ylim(min(duration) * 0.9, max(duration) * 1.1)
    pylab.xlabel(make_time_label(start_time))
    pylab.xlim(0, 24)
    pylab.xticks(range(0,24,2))
    pylab.ylabel('Template duration (sec)')
    pylab.title('Loudest: GPS=%.3f SNR=%.1f' % (gps[snrmax_index]+1e-9*gpsns[snrmax_index], snr[snrmax_index]))
    pylab.savefig('%s/%s_%d_%s_glitchgram.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()

    # Make the new SNR glitchgram
    ranges = [[ 6, 7,'b', 3],
              [ 7, 8,'g',15],
              [ 8,float('inf'),'r',15]]

    fig  = pylab.figure(2, figsize=(8.1,4.0))
    axes = pylab.Axes(fig,[0.08, 0.1, 1.0 - 0.07 - 0.03, 0.8])
    fig.add_axes(axes)

    for r in ranges:
        subset = pylab.logical_and(newsnr>=r[0], newsnr<r[1])
        print('Plotting trigs for newsnr between', r[0], r[1])
        if len(time[subset]):
            pylab.scatter(time[subset], duration[subset], edgecolor='none', c=r[2], s=r[3])
    newsnrmax_index = pylab.argmax(newsnr)
    pylab.scatter(time[newsnrmax_index], duration[newsnrmax_index], marker=(5,1,0), c='y', s=100)

    pylab.yscale('log')
    pylab.ylim(min(duration) * 0.9, max(duration) * 1.1)
    pylab.xlabel(make_time_label(start_time))
    pylab.xlim(0, 24)
    pylab.xticks(range(0,24,2))
    pylab.ylabel('Template duration (sec)')
    pylab.title('Loudest: GPS=%.3f  newSNR=%.1f' % (gps[newsnrmax_index]+1e-9*gpsns[newsnrmax_index], newsnr[newsnrmax_index]))
    pylab.savefig('%s/%s_%d_%s_new_glitchgram.png' % (config['out_dir'], ifo, veto_level, cluster))
    pylab.close()


def plot_rate_vs_time_all(config, ifo, veto_level, cluster, start_time, end_time):
    """Rate (per sec) vs. time"""

    num_min     = (end_time - start_time) / 60
    mins_in_day = range(num_min)

    fig  = pylab.figure(0, figsize=(8.1,4.0))
    axes = pylab.Axes(fig,[0.08, 0.1, 1.0 - 0.07 - 0.03, 0.8])
    fig.add_axes(axes)

    times  = [float(x) / 60.0 for x in mins_in_day]
    maxcount = 0

    for ifo in ['H1','L1','V1']:
        counts = [1e-6 for i in mins_in_day]
        f      = open('%s-%d-INSPIRAL_%s.csv' % (ifo, veto_level, cluster))

        for line in f:
            trig_end_time = int(line.split(',')[0]) 
            counts[ (trig_end_time - start_time) / 60 ] += 1

        counts   = [c / 60.0 for c in counts]
        maxcount = max(maxcount, max(counts))

        foo = {'H1':'red', 'L1':'green', 'V1':'magenta'}
        pylab.scatter(times, counts, color=foo[ifo], edgecolor='none', alpha=0.8)

    pylab.xlim(0, 24)

    if maxcount > 10:
        pylab.yscale('log')

    min_y = cluster.find('16SEC') > -1 and 1e-2 or 1
    pylab.ylim(min_y, maxcount * 1.1)
    pylab.xticks(range(0,24,2))


    pylab.grid()
    pylab.title('Rate vs. Time, 30 msec cluster, Cat 1')
    pylab.xlabel(make_time_label(start_time))
    pylab.ylabel('Rate (Hz)')
    pylab.savefig('%s_%d_%s_rate_vs_time_all.png' % (ifo, veto_level, cluster))
    pylab.close()

if __name__ == '__main__':
    actions = {'make_dag'                : make_dag,
               'make_veto_files'         : make_veto_files,
               'make_caches'             : make_caches,
               'make_csv'                : make_csv,
               'make_database'           : make_database,
               'make_glitchy_times_table': make_glitchy_times_table,
               'make_summary_table'      : make_summary_table,
               'make_usage_table'        : make_usage_table,
               'make_glitch_page'        : make_glitch_page,
               'make_hwinj_page'         : make_hwinj_page,
               'make_index_page'         : make_index_page,
               'plot_snr_vs_time'        : plot_snr_vs_time,
               'plot_rate_vs_time'       : plot_rate_vs_time,
               'plot_snr_hist'           : plot_snr_hist,
               'plot_template_counts'    : plot_template_counts,
               'plot_hexmass'            : plot_hexmass,
               'plot_bank_veto'          : plot_bank_veto,
               'plot_mass_hist'          : plot_mass_hist,
               'plot_glitchgram'         : plot_glitchgram,
               'plot_rate_vs_time_all'   : plot_rate_vs_time_all}

    options = parse_command_line()

    action     = options.action
    start_time = int(options.gps_start_time)
    end_time   = options.gps_end_time and int(options.gps_end_time) or start_time + (24*60*60)
    ifos       = options.ifos and options.ifos.split(',') or ['none']
    cats       = options.veto_categories and options.veto_categories.split(',') or [-1]
    clusters   = options.cluster_categories and options.cluster_categories.split(',') or ['none']
    config     = options.config and parse_config(options.config) or {}

    for k in ['out_dir','trigger_dir','tmp_dir']:
        if k not in config:
            config[k] = '.'


    # Don't connect to the database unless we need to, since doing so
    # brings in all the pyglobus stuff that may not be available.
    if action == 'make_veto_files':
        engine = setup_db()
        actions['make_veto_files'] = lambda config, ifo, cat, cluster, start, end: make_veto_files(engine, config, ifo, cat, cluster, start, end)

    # Run over all the ifos and categories provided
    for ifo in ifos:
        for cat in map(int, cats):
            for cluster in clusters:
                if action == 'plot_mass_hist':
                    actions[action](config, ifo, cat, cluster, start_time, end_time, options.flower)
                else:
                    actions[action](config, ifo, cat, cluster, start_time, end_time)


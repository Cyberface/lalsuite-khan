include $(top_srcdir)/gnuscripts/lalapps.am

MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

SUBDIRS = posterior

TMPLTBANKSRC = tmpltbank.c inspiralutils.c inspiral.h
INSPINJSRC = inspinj.c inspiralutils.c inspiral.h
INJCUTSRC = injcut.c
MINJSRC = minj.c
BBHINJSRC = bbhinj.c
BLINDINJSRC = blindinj.c
SPLITBANKSRC = splitbank.c
TRIGSCANSRC = trigscan.c
INSPFRINJSRC = inspfrinj.c
GWF2XMLSRC = gwf2xml.c
SPININJSRC = spininj.c
RANDOMBANKSRC = randombank.c
NINJASRC = ninja.c
FRNINJASRC = fr_ninja.c
CALCEXPSNRSRC = calcexpsnr.c
COINJSRC = coinj.c
STOCHBANKSRC = stochasticbank.c

if LALMETAIO
LALMETAIO_PROGS = \
	lalapps_tmpltbank \
	lalapps_inspinj \
	lalapps_injcut \
	lalapps_splitbank \
	lalapps_minj \
	lalapps_bbhinj \
	lalapps_blindinj \
	lalapps_inspfrinj \
	lalapps_gwf2xml \
	lalapps_trigscan \
	lalapps_spininj \
	lalapps_randombank \
	lalapps_ninja \
	lalapps_fr_ninja \
	lalapps_calcexpsnr \
	lalapps_coinj \
	lalapps_cbc_stochasticbank
lalapps_tmpltbank_SOURCES = $(TMPLTBANKSRC)
lalapps_inspinj_SOURCES = $(INSPINJSRC)
lalapps_injcut_SOURCES = $(INJCUTSRC)
lalapps_minj_SOURCES = $(MINJSRC)
lalapps_bbhinj_SOURCES = $(BBHINJSRC)
lalapps_blindinj_SOURCES = $(BLINDINJSRC)
lalapps_splitbank_SOURCES = $(SPLITBANKSRC)
lalapps_trigscan_SOURCES = $(TRIGSCANSRC)
lalapps_spininj_SOURCES = $(SPININJSRC)
lalapps_inspfrinj_SOURCES = $(INSPFRINJSRC)
lalapps_gwf2xml_SOURCES = $(GWF2XMLSRC)
lalapps_randombank_SOURCES = $(RANDOMBANKSRC)
lalapps_ninja_SOURCES = $(NINJASRC)
lalapps_fr_ninja_SOURCES = $(FRNINJASRC)
lalapps_calcexpsnr_SOURCES = $(CALCEXPSNRSRC)
lalapps_coinj_SOURCES = $(COINJSRC)
lalapps_cbc_stochasticbank_SOURCES = $(STOCHBANKSRC)
endif

bin_PROGRAMS = \
	$(LALMETAIO_PROGS)

if HAVE_PYTHON
pybin_scripts = \
	lalapps_inspinjfind \
	lalapps_inspiral_pipe \
	lalapps_inspiral_hipe \
	lalapps_inspiral_ssipe \
	lalapps_multi_hipe \
	lalapps_trigger_hipe \
	lalapps_thinca \
	lalapps_inspiral_online_pipe \
	lalapps_inspiral_online_pipe_inject \
	lalapps_ihope_daily_page \
	lalapps_flag_triggers \
	lalapps_check_flag \
	lalapps_check_hipe_times \
	lalapps_ihope \
	lalapps_followup0_pipe \
	lalapps_plot_hipe \
	lalapps_ihope_status \
	lalapps_write_ihope_page \
	lalapps_compare_ini \
	lalapps_link_old_ihope \
	lalapps_generate_upper_limits \
	lalapps_followup0_page \
	lalapps_cbc_pipedown \
	lalapps_cbc_pipedown_ssipe \
	lalapps_glitch_probe \
	lalapps_cohPTF_hipe \
	lalapps_coh_PTF_post_processing \
	lalapps_cbc_sbank \
	lalapps_cbc_sbank_choose_mchirp_boundaries \
	lalapps_cbc_sbank_hdf5_choose_mchirp_boundaries \
	lalapps_cbc_sbank_hdf5_bankcombiner \
	lalapps_cbc_sbank_pipe \
	lalapps_cbc_sbank_plot_sim \
	lalapps_cbc_sbank_sim \
	lalapps_cbc_sbank_sim_pipe \
	lalapps_make_nr_hdf_catalog \
	lalapps_run_sqlite \
	lalapps_bank_plotter

pkgpython_PYTHON = \
	inspiral.py \
	inspiralutils.py
endif

dist_bin_SCRIPTS = \
	lalapps_makeblindinj \
	lalapps_makeblindinj_himass \
	lalapps_cbc_code_compare \
	lalapps_inspiral_ckpttest \
	followup/followup_InspiralDataMover.sh

pkgdata_DATA = \
	BNSMasses.dat \
	inspsrcs100Mpc.errors \
	ihope_ssipe.ini \
	inspiral_pipe.ini \
	inspiral_hipe.ini \
	online.ini \
	online_inject_bbh.ini \
	cbc_pipedown.ini \
	followup0_pipe.ini \
	plot_hipe.ini \
	write_ihope_page.ini \
	write_ihope_style.css \
	s5_background_H0H1-RDS_R_L1-cbc.txt \
	s5_background_H0H1-RDS_R_L1-seismic-cbc.txt \
	s5_background_H0H2-RDS_R_L1-cbc.txt \
	s5_background_H0H2-RDS_R_L1-seismic-cbc.txt \
	s5_background_H1_hoft_cbc.txt \
	s5_background_H2_hoft_cbc.txt \
	s5_background_L0L1-RDS_R_L1-cbc.txt \
	s5_background_L0L1-RDS_R_L1-seismic-cbc.txt \
	s5_background_L1_hoft_cbc.txt \
	s5_background_V1-raw-cbc.txt \
	s5_background_V1-raw-seismic-cbc.txt \
	s5_background_V1_hoft_cbc.txt \
	s5_foreground_H0H1-RDS_R_L1-cbc.txt \
	s5_foreground_H0H1-RDS_R_L1-seismic-cbc.txt \
	s5_foreground_H0H2-RDS_R_L1-cbc.txt \
	s5_foreground_H0H2-RDS_R_L1-seismic-cbc.txt \
	s5_foreground_H1_hoft_cbc.txt \
	s5_foreground_H2_hoft_cbc.txt \
	s5_foreground_L0L1-RDS_R_L1-cbc.txt \
	s5_foreground_L0L1-RDS_R_L1-seismic-cbc.txt \
	s5_foreground_L1_hoft_cbc.txt \
	s5_foreground_V1-raw-cbc.txt \
	s5_foreground_V1-raw-seismic-cbc.txt \
	s5_foreground_V1_hoft_cbc.txt \
	s6_background_H0H1-RDS_R_L1-cbc.txt \
	s6_background_H0H1-RDS_R_L1-seismic-cbc.txt \
	s6_background_H1_hoft_cbc.txt \
	s6_background_L0L1-RDS_R_L1-cbc.txt \
	s6_background_L0L1-RDS_R_L1-seismic-cbc.txt \
	s6_background_L1_hoft_cbc.txt \
	s6_background_V1-raw-cbc.txt \
	s6_background_V1-raw-seismic-cbc.txt \
	s6_background_V1_hoft_cbc.txt \
	s6_foreground_H0H1-RDS_R_L1-cbc.txt \
	s6_foreground_H0H1-RDS_R_L1-seismic-cbc.txt \
	s6_foreground_H1_hoft_cbc.txt \
	s6_foreground_L0L1-RDS_R_L1-cbc.txt \
	s6_foreground_L0L1-RDS_R_L1-seismic-cbc.txt \
	s6_foreground_L1_hoft_cbc.txt \
	s6_foreground_V1-raw-cbc.txt \
	s6_foreground_V1-raw-seismic-cbc.txt \
	s6_foreground_V1_hoft_cbc.txt \
	ratioTest.pickle \
	remove_h1h2.sql \
	cluster.sql \
	simplify.sql \
	mcmcsummary.R

if HAVE_GLUE_PREFIX
gluebindir=$(GLUE_PREFIX)/bin
else
gluebindir=/usr/bin
endif

if HAVE_PYLAL_PREFIX
pylalbindir=$(PYLAL_PREFIX)/bin
else
pylalbindir=/usr/bin
endif

ihope_ssipe.ini: $(srcdir)/ihope_ssipe.ini.in
	$(AM_V_GEN)$(SED) -e 's|@LALAPPSDIR@|$(bindir)|' -e 's|@GLUEDIR@|$(gluebindir)|' -e 's|@PYLALDIR@|$(pylalbindir)|' $< > $@

CLEANFILES = \
	ihope_ssipe.ini

MOSTLYCLEANFILES += \
	*INSPIRAL*.xml \
	*INSPIRAL*.gwf \
	*TMPLTBANK*.gwf \
	*TMPLTBANK*.xml

EXTRA_DIST += \
	$(DATA) \
	$(INSPINJSRC) \
	$(MINJSRC) \
	$(SPLITBANKSRC) \
	$(TMPLTBANKSRC) \
	effective_distance_comparison.png \
	ihope_ssipe.ini.in \
	inspiral_programs.dox \
	$(END_OF_LIST)

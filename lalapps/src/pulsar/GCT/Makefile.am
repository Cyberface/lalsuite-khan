MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

AM_CPPFLAGS += -I$(top_srcdir)/src/pulsar/Fstatistic -I$(top_srcdir)/src/pulsar/HoughFstat

bin_PROGRAMS = lalapps_HierarchSearchGCT
EXTRA_PROGRAMS = lalapps_HierarchSearchGCT_SSE2 lalapps_HierarchSearchGCT_SSE2_NONC

lalapps_HierarchSearchGCT_SOURCES = \
	GCTtoplist.c \
	GCTtoplist.h \
	HeapToplist.c \
	HierarchSearchGCT.c \
	HierarchSearchGCT.h \
	RecalcToplistStats.c \
	RecalcToplistStats.h \
	$(END_OF_LIST)

lalapps_HierarchSearchGCT_SSE2_SOURCES = $(HSGCTSources)
lalapps_HierarchSearchGCT_SSE2_CPPFLAGS = $(AM_CPPFLAGS) -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT
lalapps_HierarchSearchGCT_SSE2_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse

lalapps_HierarchSearchGCT_SSE2_NONC_SOURCES = $(HSGCTSources)
lalapps_HierarchSearchGCT_SSE2_NONC_CPPFLAGS = $(AM_CPPFLAGS) -DHS_OPTIMIZATION -DHIERARCHSEARCHGCT -DGC_SSE2_OPT -DEXP_NO_NUM_COUNT
lalapps_HierarchSearchGCT_SSE2_NONC_CFLAGS = $(AM_CFLAGS) -msse -msse2 -mfpmath=sse

EXTRA_DIST = \
	testGCT.sh \
	testGCT_inject.sh \
	$(END_OF_LIST)

TESTS = \
	testGCT.sh \
	testGCT_inject.sh \
	$(END_OF_LIST)

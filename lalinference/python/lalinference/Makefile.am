BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	bayestar \
	imrtgr \
	io \
	plot \
	rapid_pe \
	tiger \
	util \
	$(END_OF_LIST)

vcs_info_sources = git_version.py

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	distance.py \
	healpix_tree.py \
	lalinference_pipe_utils.py \
	moc.py \
	nest2pos.py \
	wrapper.py \
	bayespputils.py \
	git_version.py \
	$(END_OF_LIST)

if HAVE_CHEALPIX
if SWIG_BUILD_PYTHON
pymodule_LTLIBRARIES = _distance.la _bayespputils.la _moc.la

_distance_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_distance_la_CFLAGS = $(AM_CFLAGS) $(SWIG_PYTHON_CFLAGS) -Wno-error
_distance_la_LDFLAGS = $(top_builddir)/src/liblalinference.la -shared -module -avoid-version

_moc_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_moc_la_CFLAGS = $(AM_CFLAGS) $(SWIG_PYTHON_CFLAGS) -Wno-error
_moc_la_LDFLAGS = $(top_builddir)/src/liblalinference.la -shared -module -avoid-version

_bayespputils_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_bayespputils_la_CFLAGS = $(AM_CFLAGS) $(SWIG_PYTHON_CFLAGS) -Wno-error
_bayespputils_la_LDFLAGS = $(top_builddir)/src/liblalinference.la -shared -module -avoid-version

endif
endif

all-local: _distance.so _bayespputils.so _moc.so

_distance.so:
	rm -f $@ && $(LN_S) .libs/$@

_bayespputils.so:
	rm -f $@ && $(LN_S) .libs/$@

_moc.so:
	rm -f $@ && $(LN_S) .libs/$@

CLEANFILES = _distance.so _bayespputils.so _moc.so

endif

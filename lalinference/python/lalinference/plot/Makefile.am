BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/plot

pymodule_PYTHON = \
	__init__.py \
	allsky.py \
	angle.py \
	cmap.py \
	healpix.py \
	spindisk.py \
	$(END_OF_LIST)

dist_pymodule_DATA = \
	cylon.csv \
	$(END_OF_LIST)

all-local: data-all-local
data-all-local: data-mostlyclean-local
	$(AM_V_at)if test "x$(builddir)" != "x$(srcdir)"; then \
		for file in $(dist_pymodule_DATA); do \
			if test -r "$(srcdir)/$${file}"; then \
				$(LN_S) "$(srcdir)/$${file}" "$(builddir)/$${file}" || exit 1; \
			fi; \
		done; \
	fi

mostlyclean-local: data-mostlyclean-local
data-mostlyclean-local:
	$(AM_V_at)if test "x$(builddir)" != "x$(srcdir)"; then \
		for file in $(dist_pymodule_DATA); do \
			if test -r "$(srcdir)/$${file}"; then \
				rm -f "$(builddir)/$${file}" || exit 1; \
			fi; \
		done; \
	fi

endif

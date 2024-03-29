#
# Copyright (C) 2006  Kipp Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


from __future__ import print_function


import math
import matplotlib.cm
import numpy
from optparse import OptionParser
import os
import sqlite3
import sys


from ligo.lw import dbtables
from ligo.lw import utils as ligolw_utils
import lal
from lal import rate
from lalburst import git_version
from lalburst import SnglBurstUtils
from lalburst import SimBurstUtils
from ligo import segments


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#


def parse_command_line():
	parser = OptionParser(
		version = "Name: %%prog\n%s" % git_version.verbose_msg
	)
	parser.add_option("-b", "--base", metavar = "base", default = "plotburca_", help = "Set the prefix for output filenames (default = \"plotburca_\").")
	parser.add_option("-f", "--format", metavar = "format", default = "png", help = "Set the output image format (default = \"png\").")
	parser.add_option("-l", "--live-time-program", metavar = "program", default = "lalapps_power", help = "Set the name, as it appears in the process table, of the program whose search summary entries define the search live time (default = \"lalapps_power\").")
	parser.add_option("--plot", metavar = "number", action = "append", default = None, help = "Generate the given plot number.")
	parser.add_option("-s", "--skip", metavar = "number", default = 0, help = "Skip this many files at the start of the list (default = process all).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.plot:
		options.plot = map(int, options.plot)
	else:
		options.plot = range(11)

	options.skip = int(options.skip)

	return options, (filenames or [None])


#
# =============================================================================
#
#                                Rate Contours
#
# =============================================================================
#


class RateContours(object):
	def __init__(self, x_instrument, y_instrument):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("%s Offset (s)" % x_instrument, "%s Offset (s)" % y_instrument)
		self.fig.set_size_inches(6,6)
		self.x_instrument = x_instrument
		self.y_instrument = y_instrument
		self.tisi_rows = None
		self.seglists = segments.segmentlistdict()
		self.counts = None

	def add_contents(self, contents):
		if self.tisi_rows is None:
			# get a list of time slide dictionaries
			self.tisi_rows = contents.time_slide_table.as_dict().values()

			# find the largest and smallest offsets
			min_offset = min(offset for vector in self.tisi_rows for offset in vector.values())
			max_offset = max(offset for vector in self.tisi_rows for offset in vector.values())

			# a guess at the time slide spacing:  works if the
			# time slides are distributed as a square grid over
			# the plot area.  (max - min)^2 gives the area of
			# the time slide square in square seconds; dividing
			# by the length of the time slide list gives the
			# average area per time slide;  taking the square
			# root of that gives the average distance between
			# adjacent time slides in seconds
			time_slide_spacing = ((max_offset - min_offset)**2 / len(self.tisi_rows))**0.5

			# use an average of 3 bins per time slide in each
			# direction, but round to an odd integer
			nbins = int(math.ceil((max_offset - min_offset) / time_slide_spacing * 3))

			# construct the binning
			self.counts = rate.BinnedArray(rate.NDBins((rate.LinearBins(min_offset, max_offset, nbins), rate.LinearBins(min_offset, max_offset, nbins))))

		self.seglists |= contents.seglists

		for offsets in contents.connection.cursor().execute("""
SELECT tx.offset, ty.offset FROM
	coinc_event
	JOIN time_slide AS tx ON (
		tx.time_slide_id == coinc_event.time_slide_id
	)
	JOIN time_slide AS ty ON (
		ty.time_slide_id == coinc_event.time_slide_id
	)
WHERE
	coinc_event.coinc_def_id == ?
	AND tx.instrument == ?
	AND ty.instrument == ?
		""", (contents.bb_definer_id, self.x_instrument, self.y_instrument)):
			try:
				self.counts[offsets] += 1
			except IndexError:
				# beyond plot boundaries
				pass

	def finish(self):
		livetime = rate.BinnedArray(self.counts.bins)
		for offsets in self.tisi_rows:
			self.seglists.offsets.update(offsets)
			livetime[offsets[self.x_instrument], offsets[self.y_instrument]] += float(abs(self.seglists.intersection(self.seglists.keys())))
		zvals = self.counts.at_centres() / livetime.at_centres()
		rate.filter_array(zvals, rate.gaussian_window(3, 3))
		xcoords, ycoords = self.counts.centres()
		self.axes.contour(xcoords, ycoords, numpy.transpose(numpy.log(zvals)))
		for offsets in self.tisi_rows:
			if any(offsets):
				# time slide vector is non-zero-lag
				self.axes.plot((offsets[self.x_instrument],), (offsets[self.y_instrument],), "k+")
			else:
				# time slide vector is zero-lag
				self.axes.plot((offsets[self.x_instrument],), (offsets[self.y_instrument],), "r+")

		self.axes.set_xlim([self.counts.bins().min[0], self.counts.bins().max[0]])
		self.axes.set_ylim([self.counts.bins().min[1], self.counts.bins().max[1]])
		self.axes.set_title(r"Coincident Event Rate vs. Instrument Time Offset (Logarithmic Rate Contours)")


#
# =============================================================================
#
#                              Confidence Scatter
#
# =============================================================================
#

def incomplete_injection_coincs(contents):
	# determine burst <--> burst coincidences for which at least one
	# burst, *but not all*, was identified as an injection;  these are
	# places in the data where an injection was done, a coincident
	# event was seen, but where, later, the injection was not found to
	# match all events in the coincidence;  these perhaps indicate
	# power leaking from the injection into nearby tiles, or accidental
	# coincidence of a marginal injection with near-by noise, etc, and
	# so although they aren't "bang-on" reconstructions of injections
	# they are nevertheless injections that are found and survive a
	# coincidence cut.
	for values in contents.connection.cursor().execute("""
SELECT
	bb_coinc_event.*
FROM
	coinc_event AS bb_coinc_event
WHERE
	bb_coinc_event.coinc_def_id == ?
	AND EXISTS (
		SELECT
			*
		FROM
			coinc_event AS sb_coinc_event
			JOIN coinc_event_map AS a ON (
				a.coinc_event_id == sb_coinc_event.coinc_event_id
				AND a.table_name == 'sngl_burst'
			)
			JOIN coinc_event_map AS b ON (
				b.table_name == 'sngl_burst'
				AND b.event_id == a.event_id
			)
		WHERE
			sb_coinc_event.coinc_def_id == ?
			AND b.coinc_event_id == bb_coinc_event.coinc_event_id
			AND sb_coinc_event.time_slide_id == bb_coinc_event.time_slide_id
	)
EXCEPT
	SELECT
		c.event_id
	FROM
		coinc_event_map AS c
		JOIN coinc_event AS sc_coinc_event ON (
			c.coinc_event_id == sc_coinc_event.coinc_event_id
			AND c.table_name == 'coinc_event'
		)
	WHERE
		sc_coinc_event.coinc_def_id == ?
		AND sc_coinc_event.time_slide_id == bb_coinc_event.time_slide_id
	""", (contents.bb_definer_id, contents.sb_definer_id, contents.scn_definer_id)):
		yield contents.coinc_table.row_from_cols(values)


def magnitude_a(burst):
	return burst.ms_hrss


def magnitude_b(burst):
	gmst = lal.GreenwichMeanSiderealTime(burst.peak)
	fplus, fcross = lal.ComputeDetAMResponse(lal.cached_detector_by_prefix[burst.ifo].response, SimBurstUtils.MW_CENTER_J2000_RA_RAD, SimBurstUtils.MW_CENTER_J2000_DEC_RAD, 0.0, gmst)
	return (burst.ms_hrss**2.0 / (fplus**2.0 + fcross**2.0))**0.5


class ConfidenceContours(object):
	def __init__(self, x_instrument, y_instrument, magnitude, desc, min_magnitude, max_magnitude):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("%s %s" % (x_instrument, desc), "%s %s" % (y_instrument, desc))
		self.fig.set_size_inches(6,6)
		self.axes.loglog()
		self.x_instrument = x_instrument
		self.y_instrument = y_instrument
		self.magnitude = magnitude
		self.foreground_x = []
		self.foreground_y = []
		self.n_foreground = 0
		self.n_background = 0
		self.n_injections = 0
		self.foreground_bins = rate.BinnedArray(rate.NDBins((rate.LogarithmicBins(min_magnitude, max_magnitude, 1024), rate.LogarithmicBins(min_magnitude, max_magnitude, 1024))))
		self.background_bins = rate.BinnedArray(rate.NDBins((rate.LogarithmicBins(min_magnitude, max_magnitude, 1024), rate.LogarithmicBins(min_magnitude, max_magnitude, 1024))))
		self.coinc_injection_bins = rate.BinnedArray(rate.NDBins((rate.LogarithmicBins(min_magnitude, max_magnitude, 1024), rate.LogarithmicBins(min_magnitude, max_magnitude, 1024))))
		self.incomplete_coinc_injection_bins = rate.BinnedArray(rate.NDBins((rate.LogarithmicBins(min_magnitude, max_magnitude, 1024), rate.LogarithmicBins(min_magnitude, max_magnitude, 1024))))

	def coords(self, contents, coinc_event_id):
		x = y = None
		for burst in SnglBurstUtils.coinc_sngl_bursts(contents, coinc_event_id):
			if burst.ifo == self.x_instrument:
				x = self.magnitude(burst)
			elif burst.ifo == self.y_instrument:
				y = self.magnitude(burst)
		return x, y

	def add_contents(self, contents):
		if self.x_instrument not in contents.instruments or self.y_instrument not in contents.instruments:
			# fast path: these are not the instruments we're
			# looking for
			return
		for values in contents.connection.cursor().execute("""
SELECT
	coinc_event.*,
	EXISTS (
		SELECT
			*
		FROM
			time_slide
		WHERE
			time_slide.time_slide_id == coinc_event.time_slide_id
			AND time_slide.offset != 0
	)
FROM
	coinc_event
WHERE
	coinc_event.coinc_def_id == ?
		""", (contents.bb_definer_id,)):
			coinc = contents.coinc_table.row_from_cols(values)
			not_zero_lag = values[-1]
			if not_zero_lag:
				self.n_foreground += 1
				x, y = self.coords(contents, coinc.coinc_event_id)
				try:
					self.foreground_bins[x, y] += 1
					#self.foreground_x.append(x)
					#self.foreground_y.append(y)
				except IndexError:
					# not on plot axes
					pass
			else:
				self.n_background += 1
				try:
					self.background_bins[self.coords(contents, coinc.coinc_event_id)] += 1
				except IndexError:
					# not on plot axes
					pass
		if contents.sce_definer_id is not None:
			# this query assumes each injection can match at
			# most 1 coinc
			for values in contents.connection.cursor().execute("""
SELECT
	coinc_event.*
FROM
	coinc_event AS sim_coinc
	JOIN coinc_event_map ON (
		coinc_event_map.coinc_event_id == sim_coinc.coinc_event_id
	)
	JOIN coinc_event ON (
		coinc_event_map.table_name == 'coinc_event'
		AND coinc_event_map.event_id == coinc_event.coinc_event_id
	)
WHERE
	sim_coinc.coinc_def_id == ?
			""", (contents.sce_definer_id,)):
				coinc = contents.coinc_event_table.row_from_cols(values)
				self.n_injections += 1
				try:
					self.coinc_injection_bins[self.coords(contents, coinc.coinc_event_id)] += 1
				except IndexError:
					# not on plot axes
					pass
		for coinc in incomplete_injection_coincs(contents):
			try:
				self.incomplete_coinc_injection_bins[self.coords(contents, coinc.coinc_event_id)] += 1
			except IndexError:
				# not on plot axes
				pass

	def finish(self):
		self.axes.set_title(r"\begin{center}Distribution of Coincident Events (%d Foreground, %d Background Events, %d Injections Found in Coincidence, Logarithmic Density Contours)\end{center}" % (self.n_foreground, self.n_background, self.n_injections))
		xcoords, ycoords = self.background_bins.centres()

		# prepare the data
		rate.filter_array(self.foreground_bins.array, rate.gaussian_window(8, 8))
		rate.filter_array(self.background_bins.array, rate.gaussian_window(8, 8))
		rate.filter_array(self.coinc_injection_bins.array, rate.gaussian_window(8, 8))
		rate.filter_array(self.incomplete_coinc_injection_bins.array, rate.gaussian_window(8, 8))
		self.foreground_bins.logregularize()
		self.background_bins.logregularize()
		self.coinc_injection_bins.logregularize()
		self.incomplete_coinc_injection_bins.logregularize()

		# plot background contours
		max_density = math.log(self.background_bins.array.max())
		self.axes.contour(xcoords, ycoords, numpy.transpose(numpy.log(self.background_bins.array)), sorted(max_density - n for n in xrange(0, 10, 1)), cmap = matplotlib.cm.Greys)

		# plot foreground (zero-lag) contours
		max_density = math.log(self.foreground_bins.array.max())
		self.axes.contour(xcoords, ycoords, numpy.transpose(numpy.log(self.foreground_bins.array)), sorted(max_density - n for n in xrange(0, 10, 1)), cmap = matplotlib.cm.Reds)
		#self.axes.plot(self.foreground_x, self.foreground_y, "r+")

		# plot coincident injection contours
		max_density = math.log(self.coinc_injection_bins.array.max())
		self.axes.contour(xcoords, ycoords, numpy.transpose(numpy.log(self.coinc_injection_bins.array)), sorted(max_density - n for n in xrange(0, 10, 1)), cmap = matplotlib.cm.Blues)

		# plot incomplete coincident injection contours
		max_density = math.log(self.incomplete_coinc_injection_bins.array.max())
		self.axes.contour(xcoords, ycoords, numpy.transpose(numpy.log(self.incomplete_coinc_injection_bins.array)), sorted(max_density - n for n in xrange(0, 10, 1)), cmap = matplotlib.cm.Greens)

		# add diagonal line
		lower = max(binning.min for binning in self.background_bins.bins)
		upper = min(binning.max for binning in self.background_bins.bins)
		self.axes.plot([lower, upper], [lower, upper], "k:")

		# fix axes limits
		self.axes.set_xlim([self.background_bins.bins[0].min, self.background_bins.bins[0].max])
		self.axes.set_ylim([self.background_bins.bins[1].min, self.background_bins.bins[1].max])


class ConfidenceContourProjection(object):
	def __init__(self, x, y, magnitude, max_magnitude):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("X", "Y")
		self.fig.set_size_inches(6,6)
		self.x = x
		self.y = y
		self.magnitude = magnitude
		self.n_foreground = 0
		self.n_background = 0
		self.n_injections = 0
		max_magnitude = math.log10(max_magnitude)
		self.foreground_bins = rate.BinnedArray(rate.NDBins((rate.LinearBins(-max_magnitude, max_magnitude, 1024), rate.LinearBins(-max_magnitude, max_magnitude, 1024))))
		self.background_bins = rate.BinnedArray(rate.NDBins((rate.LinearBins(-max_magnitude, max_magnitude, 1024), rate.LinearBins(-max_magnitude, max_magnitude, 1024))))
		self.coinc_injection_bins = rate.BinnedArray(rate.NDBins((rate.LinearBins(-max_magnitude, max_magnitude, 1024), rate.LinearBins(-max_magnitude, max_magnitude, 1024))))
		self.incomplete_coinc_injection_bins = rate.BinnedArray(rate.NDBins((rate.LinearBins(-max_magnitude, max_magnitude, 1024), rate.LinearBins(-max_magnitude, max_magnitude, 1024))))

	def coords(self, contents, coinc_event_id):
		mag = numpy.zeros(3, "Float64")
		for burst in SnglBurstUtils.coinc_sngl_bursts(contents, coinc_event_id):
			if burst.ifo == "H1":
				mag[0] = math.log10(self.magnitude(burst))
			elif burst.ifo == "H2":
				mag[1] = math.log10(self.magnitude(burst))
			elif burst.ifo == "L1":
				mag[2] = math.log10(self.magnitude(burst))
		return numpy.inner(mag, self.x), numpy.inner(mag, self.y)

	def add_contents(self, contents):
		for values in contents.connection.cursor().execute("""
SELECT
	coinc_event.*,
	EXISTS (
		SELECT
			*
		FROM
			time_slide
		WHERE
			time_slide.time_slide_id == coinc_event.time_slide_id
			AND time_slide.offset != 0
	)
FROM
	coinc_event
WHERE
	coinc_event.coinc_def_id == ?
		""", (contents.bb_definer_id,)):
			coinc = contents.coinc_table.row_from_cols(values)
			not_zero_lag = values[-1]
			x, y = self.coords(contents, coinc.coinc_event_id)
			if not_zero_lag:
				self.n_background += 1
				try:
					self.background_bins[x, y] += 1
				except IndexError:
					# not on plot axes
					pass
			else:
				self.n_foreground += 1
				try:
					self.foreground_bins[x, y] += 1
					#self.foreground_x.append(x)
					#self.foreground_y.append(y)
				except IndexError:
					# not on plot axes
					pass
		if contents.sce_definer_id is not None:
			# this query assumes each injection can match at
			# most 1 coinc
			for values in contents.connection.cursor().execute("""
SELECT
	coinc_event.*
FROM
	coinc_event AS sim_coinc
	JOIN coinc_event_map ON (
		coinc_event_map.coinc_event_id == sim_coinc.coinc_event_id
	)
	JOIN coinc_event ON (
		coinc_event_map.table_name == 'coinc_event'
		AND coinc_event_map.event_id == coinc_event.coinc_event_id
	)
WHERE
	sim_coinc.coinc_def_id == ?
			""", (contents.sce_definer_id,)):
				coinc = contents.coinc_event_table.row_from_cols(values)
				self.n_injections += 1
				try:
					self.coinc_injection_bins[self.coords(contents, coinc.coinc_event_id)] += 1
				except IndexError:
					# not on plot axes
					pass
		for coinc in incomplete_injection_coincs(contents):
			try:
				self.incomplete_coinc_injection_bins[self.coords(contents, coinc.coinc_event_id)] += 1
			except IndexError:
				# not on plot axes
				pass

	def finish(self):
		self.axes.set_title(r"\begin{center}Distribution of Coincident Events (%d Foreground, %d Background Events, %d Injections Found in Coincidence, Logarithmic Density Contours)\end{center}" % (self.n_foreground, self.n_background, self.n_injections))
		xcoords, ycoords = self.background_bins.centres()

		# prepare the data
		rate.filter_array(self.foreground_bins.array, rate.gaussian_window(8, 8))
		rate.filter_array(self.background_bins.array, rate.gaussian_window(8, 8))
		rate.filter_array(self.coinc_injection_bins.array, rate.gaussian_window(8, 8))
		rate.filter_array(self.incomplete_coinc_injection_bins.array, rate.gaussian_window(8, 8))
		self.foreground_bins.logregularize()
		self.background_bins.logregularize()
		self.coinc_injection_bins.logregularize()
		self.incomplete_coinc_injection_bins.logregularize()

		# plot background contours
		max_density = math.log(self.background_bins.array.max())
		self.axes.contour(xcoords, ycoords, numpy.transpose(numpy.log(self.background_bins.array)), [max_density - n for n in xrange(0, 10, 1)], cmap = matplotlib.cm.Greys)

		# plot foreground (zero-lag) contours
		max_density = math.log(self.foreground_bins.array.max())
		self.axes.contour(xcoords, ycoords, numpy.transpose(numpy.log(self.foreground_bins.array)), [max_density - n for n in xrange(0, 10, 1)], cmap = matplotlib.cm.Reds)
		#self.axes.plot(self.foreground_x, self.foreground_y, "r+")

		# plot coincident injection contours
		max_density = math.log(self.coinc_injection_bins.array.max())
		self.axes.contour(xcoords, ycoords, numpy.transpose(numpy.log(self.coinc_injection_bins.array)), [max_density - n for n in xrange(0, 10, 1)], cmap = matplotlib.cm.Blues)

		# plot incomplete coincident injection contours
		max_density = math.log(self.incomplete_coinc_injection_bins.array.max())
		self.axes.contour(xcoords, ycoords, numpy.transpose(numpy.log(self.incomplete_coinc_injection_bins.array)), [max_density - n for n in xrange(0, 10, 1)], cmap = matplotlib.cm.Greens)

		# fix axes limits
		self.axes.set_xlim([self.background_bins.bins.min[0], self.background_bins.bins.max[0]])
		self.axes.set_ylim([self.background_bins.bins.min[1], self.background_bins.bins.max[1]])


#
# =============================================================================
#
#                             Rate vs. Confidence
#
# =============================================================================
#


class RateVsConfidence(object):
	def __init__(self, instrument):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("%s Confidence" % instrument, "Coincident Event Rate (Hz)")
		self.instrument = instrument
		self.foreground = []
		self.background = []
		self.foreground_segs = segments.segmentlist()
		self.background_segs = segments.segmentlist()
		self.axes.loglog()

	def add_contents(self, contents):
		for time_slide_id, not_zero_lag in contents.connection.cursor().execute("SELECT DISTINCT time_slide_id, EXISTS ( SELECT * FROM time_slide AS a WHERE a.time_slide_id == time_slide.time_slide_id AND a.offset != 0 ) FROM time_slide"):
			if not_zero_lag:
				bins = self.background
				self.background_segs.append(contents.seglists[self.instrument])
			else:
				bins = self.foreground
				self.foreground_segs.append(contents.seglists[self.instrument])
			bins.extend(contents.connection.cursor().execute("SELECT sngl_burst.confidence FROM coinc_event JOIN coinc_event_map ON (coinc_event_map.coinc_event_id == coinc_event.coinc_event_id) JOIN sngl_burst ON (coinc_event_map.table_name == 'sngl_burst' AND coinc_event_map.event_id == sngl_burst.event_id) WHERE coinc_event.time_slide_id == ? AND coinc_event.coinc_def_id == ? AND sngl_burst.ifo == ?""", (time_slide_id, contents.bb_definer_id, self.instrument)))

	def finish(self):
		self.axes.set_title("Cummulative Coincident Event Rate vs. Confidence in %s" % self.instrument)
		self.background.sort()
		self.foreground.sort()
		background_y = numpy.arange(len(self.background), 0.0, -1.0, "Float64") / float(abs(self.background_segs))
		foreground_y = numpy.arange(len(self.foreground), 0.0, -1.0, "Float64") / float(abs(self.foreground_segs))
		self.axes.plot(self.background, background_y, "ko-")
		self.axes.plot(self.foreground, foreground_y, "ro-", markeredgecolor = "r")


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


def new_plots(plots):
	deltat_seg = segments.segment(-0.3, +0.3)
	deltat_width = 0.03125
	l = [
		RateContours("H2", "H1"),
		ConfidenceContours("H2", "H1", magnitude_a, "Confidence", 1, 10**10),
		ConfidenceContours("H2", "L1", magnitude_a, "Confidence", 1, 10**10),
		ConfidenceContours("L1", "H1", magnitude_a, "Confidence", 1, 10**10),
		ConfidenceContours("H2", "H1", magnitude_b, r"Power / D.o.F. / ($F_{+}^{2} + F_{\times}^{2}$)", 1, 10**10),
		ConfidenceContours("H2", "L1", magnitude_b, r"Power / D.o.F. / ($F_{+}^{2} + F_{\times}^{2}$)", 1, 10**10),
		ConfidenceContours("L1", "H1", magnitude_b, r"Power / D.o.F. / ($F_{+}^{2} + F_{\times}^{2}$)", 1, 10**10),
		ConfidenceContourProjection(numpy.array((-1/math.sqrt(2), +1/math.sqrt(2), 0), "Float64"), numpy.array((-1/math.sqrt(4), -1/math.sqrt(4), +1/math.sqrt(2)), "Float64"), magnitude_b, 10**5),
		RateVsConfidence("H1"),
		RateVsConfidence("H2"),
		RateVsConfidence("L1")
	]
	return [l[i] for i in plots]


options, filenames = parse_command_line()


plots = new_plots(options.plot)


for n, filename in enumerate(ligolw_utils.sort_files_by_size(filenames, options.verbose, reverse = True)[options.skip:]):
	if options.verbose:
		print("%d/%d: %s" % (n + 1, len(filenames) - options.skip, filename), file=sys.stderr)

	database = SnglBurstUtils.CoincDatabase(sqlite3.connect(filename), options.live_time_program)
	if options.verbose:
		SnglBurstUtils.summarize_coinc_database(database)

	for n, plot in zip(options.plot, plots):
		if options.verbose:
			print("adding to burca plot %d ..." % n, file=sys.stderr)
		plot.add_contents(database)

	database.connection.close()


# delete the plots as we go to save memory
n = 0
format = "%%s%%0%dd.%%s" % (int(math.log10(max(options.plot) or 1)) + 1)
while len(plots):
	filename = format % (options.base, options.plot[n], options.format)
	if options.verbose:
		print("finishing plot %d ..." % options.plot[n], file=sys.stderr)
	plots[0].finish()
	if options.verbose:
		print("writing %s ..." % filename, file=sys.stderr)
	plots[0].fig.savefig(filename)
	del plots[0]
	n += 1

if options.verbose:
	print("done.", file=sys.stderr)

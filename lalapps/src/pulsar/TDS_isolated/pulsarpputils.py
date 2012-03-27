# -*- coding: utf-8 -*-
#
#       pulsarpputils.py
#
#       Copyright 2012
#       Matthew Pitkin <matthew.pitkin@ligo.org>
#
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

# known pulsar analysis post-processing utilities

# Many functions in this a taken from, or derived from equivalents available in
# the PRESTO pulsar software package http://www.cv.nrao.edu/~sransom/presto/

import sys
import math
import os
import numpy as np
#import matplotlib

from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.mlab import specgram, find
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

from types import StringType, FloatType

# some common constants taken from psr_constants.py in PRESTO
ARCSECTORAD = float('4.8481368110953599358991410235794797595635330237270e-6')
RADTOARCSEC = float('206264.80624709635515647335733077861319665970087963')
SECTORAD    = float('7.2722052166430399038487115353692196393452995355905e-5')
RADTOSEC    = float('13750.987083139757010431557155385240879777313391975')
RADTODEG    = float('57.295779513082320876798154814105170332405472466564')
DEGTORAD    = float('1.7453292519943295769236907684886127134428718885417e-2')
RADTOHRS    = float('3.8197186342054880584532103209403446888270314977710')
HRSTORAD    = float('2.6179938779914943653855361527329190701643078328126e-1')
PI          = float('3.1415926535897932384626433832795028841971693993751')
TWOPI       = float('6.2831853071795864769252867665590057683943387987502')
PIBYTWO     = float('1.5707963267948966192313216916397514420985846996876')
SECPERDAY   = float('86400.0')
SECPERJULYR = float('31557600.0')
KMPERPC     = float('3.0856776e13')
KMPERKPC    = float('3.0856776e16')
Tsun        = float('4.925490947e-6') # sec
Msun        = float('1.9891e30')      # kg
Mjup        = float('1.8987e27')      # kg
Rsun        = float('6.9551e8')       # m
Rearth      = float('6.378e6')        # m
SOL         = float('299792458.0')    # m/s
MSUN        = float('1.989e+30')      # kg
G           = float('6.673e-11')      # m^3/s^2/kg 
C           = SOL
KPC         = float('3.0856776e19')   # kiloparsec in metres
I38         = float('1e38')           # moment of inertia kg m^2 

# some angle conversion functions taken from psr_utils.py in PRESTO
def rad_to_dms(rad):
  """
  rad_to_dms(rad):
     Convert radians to degrees, minutes, and seconds of arc.
  """
  if (rad < 0.0): sign = -1
  else: sign = 1
  arc = RADTODEG * np.fmod(np.fabs(rad), math.pi)
  d = int(arc)
  arc = (arc - d) * 60.0
  m = int(arc)
  s = (arc - m) * 60.0
  if sign==-1 and d==0:
    return (sign * d, sign * m, sign * s)
  else:
    return (sign * d, m, s)

def dms_to_rad(deg, min, sec):
  """
  dms_to_rad(deg, min, sec):
     Convert degrees, minutes, and seconds of arc to radians.
  """
  if (deg < 0.0):
    sign = -1
  elif (deg==0.0 and (min < 0.0 or sec < 0.0)):
    sign = -1
  else:
    sign = 1
  return sign * ARCSECTORAD * \
    (60.0 * (60.0 * np.fabs(deg) + np.fabs(min)) + np.fabs(sec))

def dms_to_deg(deg, min, sec):
  """
  dms_to_deg(deg, min, sec):
     Convert degrees, minutes, and seconds of arc to degrees.
  """
  return RADTODEG * dms_to_rad(deg, min, sec)

def rad_to_hms(rad):
  """
  rad_to_hms(rad):
     Convert radians to hours, minutes, and seconds of arc.
  """
  rad = np.fmod(rad, 2.*math.pi)
  if (rad < 0.0): rad = rad + 2.*math.pi
  arc = RADTOHRS * rad
  h = int(arc)
  arc = (arc - h) * 60.0
  m = int(arc)
  s = (arc - m) * 60.0
  return (h, m, s)
  
def hms_to_rad(hour, min, sec):
  """
  hms_to_rad(hour, min, sec):
     Convert hours, minutes, and seconds of arc to radians
  """
  if (hour < 0.0): sign = -1
  else: sign = 1
  return sign * SECTORAD * \
         (60.0 * (60.0 * np.fabs(hour) + np.fabs(min)) + np.fabs(sec))

def coord_to_string(h_or_d, m, s):
  """
  coord_to_string(h_or_d, m, s):
     Return a formatted string of RA or DEC values as
     'hh:mm:ss.ssss' if RA, or 'dd:mm:ss.ssss' if DEC.
  """
  retstr = ""
  if h_or_d < 0:
    retstr = "-"
  elif abs(h_or_d)==0:
    if (m < 0.0) or (s < 0.0):
      retstr = "-"
  
  h_or_d, m, s = abs(h_or_d), abs(m), abs(s)
  if (s >= 9.9995):
    return retstr+"%.2d:%.2d:%.4f" % (h_or_d, m, s)
  else:
    return retstr+"%.2d:%.2d:0%.4f" % (h_or_d, m, s)

def ra_to_rad(ra_string):
  """
  ra_to_rad(ar_string):
     Given a string containing RA information as
     'hh:mm:ss.ssss', return the equivalent decimal
     radians.
  """
  h, m, s = ra_string.split(":")
  return hms_to_rad(int(h), int(m), float(s))

def dec_to_rad(dec_string):
  """
  dec_to_rad(dec_string):
     Given a string containing DEC information as
     'dd:mm:ss.ssss', return the equivalent decimal
     radians.
  """
  d, m, s = dec_string.split(":")
  if "-" in d and int(d)==0:
    m, s = '-'+m, '-'+s
  
  return dms_to_rad(int(d), int(m), float(s))

def p_to_f(p, pd, pdd=None):
  """
  p_to_f(p, pd, pdd=None):
    Convert period, period derivative and period second
    derivative to the equivalent frequency counterparts.
    Will also convert from f to p.
  """
  f = 1.0 / p
  fd = -pd / (p * p)
  if (pdd==None):
    return [f, fd]
  else:
    if (pdd==0.0):
      fdd = 0.0
    else:
      fdd = 2.0 * pd * pd / (p**3.0) - pdd / (p * p)
     
    return [f, fd, fdd]

def pferrs(porf, porferr, pdorfd=None, pdorfderr=None):
  """
  pferrs(porf, porferr, pdorfd=None, pdorfderr=None):
     Calculate the period or frequency errors and
     the pdot or fdot errors from the opposite one.
  """
  if (pdorfd==None):
    return [1.0 / porf, porferr / porf**2.0]
  else:
    forperr = porferr / porf**2.0
    fdorpderr = np.sqrt((4.0 * pdorfd**2.0 * porferr**2.0) / porf**6.0 +
                          pdorfderr**2.0 / porf**4.0)
    [forp, fdorpd] = p_to_f(porf, pdorfd)
  
    return [forp, forperr, fdorpd, fdorpderr]

# class to read in a pulsar par file - this is heavily based on the function
# in parfile.py in PRESTO
float_keys = ["F", "F0", "F1", "F2", "F3", "F4", "F5", "F6",
              "P", "P0", "P1", "P2", "P3", "P4", "P5", "P6",
              "PEPOCH", "POSEPOCH", "DM", "START", "FINISH", "NTOA",
              "TRES", "TZRMJD", "TZRFRQ", "TZRSITE", "NITS",
              "A1", "XDOT", "E", "ECC", "EDOT", "T0", "PB", "PBDOT", "OM",
              "OMDOT", "EPS1", "EPS2", "EPS1DOT", "EPS2DOT", "TASC", "LAMBDA",
              "BETA", "RA_RAD", "DEC_RAD", "GAMMA", "SINI", "M2", "MTOT",
              "FB0", "FB1", "FB2", "ELAT", "ELONG", "PMRA", "PMDEC", "DIST",
              # GW PARAMETERS
              "H0", "COSIOTA", "PSI", "PHI0", "THETA", "I21", "I31"]
str_keys = ["FILE", "PSR", "PSRJ", "NAME", "RAJ", "DECJ", "RA", "DEC", "EPHEM",
            "CLK", "BINARY"]

class psr_par:
  def __init__(self, parfilenm):
    self.FILE = parfilenm
    pf = open(parfilenm)
    for line in pf.readlines():
      # Convert any 'D-' or 'D+' to 'E-' or 'E+'
      line = line.replace("D-", "E-")
      line = line.replace("D+", "E+")
      # also check for lower case
      line = line.replace("d-", "e-")
      line = line.replace("d+", "e+")
      splitline = line.split()
      
      # get all upper case version in case lower case in par file
      key = splitline[0].upper()
      
      if key in str_keys:
        setattr(self, key, splitline[1])
      elif key in float_keys:
        try:
          setattr(self, key, float(splitline[1]))
        except ValueError:
          pass
      
      if len(splitline)==3: # Some parfiles don't have flags, but do have errors
        if splitline[2] not in ['0', '1']:
          setattr(self, key+'_ERR', float(splitline[2]))
      
      if len(splitline)==4:
        setattr(self, key+'_ERR', float(splitline[3]))
     
    # sky position
    if hasattr(self, 'RAJ'):
      setattr(self, 'RA_RAD', ra_to_rad(self.RAJ))
      
      # set RA error in rads (rather than secs)
      if hasattr(self, 'RAJ_ERR'):
        setattr(self, 'RA_RAD_ERR', hms_to_rad(0, 0, self.RAJ_ERR))
    if hasattr(self, 'DECJ'):
      setattr(self, 'DEC_RAD', dec_to_rad(self.DECJ))
      
      # set DEC error in rads (rather than arcsecs)
      if hasattr(self, 'DECJ_ERR'):
        setattr(self, 'DEC_RAD_ERR', dms_to_rad(0, 0, self.DECJ_ERR))
      
    # periods and frequencies
    if hasattr(self, 'P'):
      setattr(self, 'P0', self.P)
    if hasattr(self, 'P0'):
      setattr(self, 'F0', 1.0/self.P0)
    if hasattr(self, 'F0'):
      setattr(self, 'P0', 1.0/self.F0)
    if hasattr(self, 'FB0'):
      setattr(self, 'PB', (1.0/self.FB0)/86400.0)
    if hasattr(self, 'P0_ERR'):
      if hasattr(self, 'P1_ERR'):
        f, ferr, fd, fderr = pferrs(self.P0, self.P0_ERR,
                                       self.P1, self.P1_ERR)
        setattr(self, 'F0_ERR', ferr) 
        setattr(self, 'F1', fd) 
        setattr(self, 'F1_ERR', fderr) 
      else:
        f, fd, = p_to_f(self.P0, self.P1)
        setattr(self, 'F0_ERR', self.P0_ERR/(self.P0*self.P0))
        setattr(self, 'F1', fd) 
    if hasattr(self, 'F0_ERR'):
      if hasattr(self, 'F1_ERR'):
        p, perr, pd, pderr = pferrs(self.F0, self.F0_ERR,
                                    self.F1, self.F1_ERR)
        setattr(self, 'P0_ERR', perr) 
        setattr(self, 'P1', pd) 
        setattr(self, 'P1_ERR', pderr) 
      else:
        p, pd, = p_to_f(self.F0, self.F1)
        setattr(self, 'P0_ERR', self.F0_ERR/(self.F0*self.F0))
        setattr(self, 'P1', pd) 
    
    # binary parameters
    if hasattr(self, 'EPS1') and hasattr(self, 'EPS2'):
      ecc = math.sqrt(self.EPS1 * self.EPS1 + self.EPS2 * self.EPS2)
      omega = math.atan2(self.EPS1, self.EPS2)
      setattr(self, 'E', ecc)
      setattr(self, 'OM', omega * RADTODEG)
      setattr(self, 'T0', self.TASC + self.PB * omega/TWOPI)
    if hasattr(self, 'PB') and hasattr(self, 'A1') and not \
       (hasattr(self, 'E') or hasattr(self, 'ECC')):
      setattr(self, 'E', 0.0)  
    if hasattr(self, 'T0') and not hasattr(self, 'TASC'):
      setattr(self, 'TASC', self.T0 - self.PB * self.OM/360.0)
        
    pf.close()
  
  def __getitem__(self, key):
    try:
      par = getattr(self, key)
    except:
      par = None
      
    return par
  
  def __str__(self):
    out = ""
    for k, v in self.__dict__.items():
      if k[:2]!="__":
        if type(self.__dict__[k]) is StringType:
          out += "%10s = '%s'\n" % (k, v)
        else:
          out += "%10s = %-20.15g\n" % (k, v)
      
    return out

# class to read in a nested sampling prior file   
class psr_prior:
  def __init__(self, priorfilenm):
    self.FILE = priorfilenm
    pf = open(priorfilenm)
    for line in pf.readlines():
      splitline = line.split()
      
      # get all upper case version in case lower case in par file
      key = splitline[0].upper()
      
      if key in str_keys:
        # everything in a prior files should be numeric
        setattr(self, key, [float(splitline[1]), float(splitline[2])])
      elif key in float_keys:
        setattr(self, key, [float(splitline[1]), float(splitline[2])])        
 
    # get sky positions in rads as strings 'dd/hh:mm:ss.s'  
    if hasattr(self, 'RA'):
      hl, ml, sl = rad_to_hms(self.RA[0])
      rastrl = coord_to_string(hl, ml, sl)
      hu, mu, su = rad_to_hms(self.RA[1])
      rastru = coord_to_string(hu, mu, su)
      setattr(self, 'RA_STR', [rastrl, rastru])
      
    if hasattr(self, 'DEC'):
      dl, ml, sl = rad_to_dms(self.DEC[0])
      decstrl = coord_to_string(dl, ml, sl)
      du, mu, su = rad_to_dms(self.DEC[1])
      decstru = coord_to_string(du, mu, su)
      setattr(self, 'DEC_STR', [decstrl, decstru])
    
    pf.close()
  
  def __getitem__(self, key):
    try:
      atr = getattr(self, key)
    except:
      atr = None
    
    return atr
  
  def __str__(self):
    out = ""
    for k, v in self.__dict__.items():
      if k[:2]!="__":
        if type(self.__dict__[k]) is StringType:
          out += "%10s = '%s'\n" % (k, v)
        else:
          out += "%10s = %-20.15g, %-20.15g\n" % (k, float(v[0]), float(v[1]))
      
    return out

# Function to return a pulsar's strain spin-down limit given its spin frequency
#(Hz), spin-down (Hz/s) and distance (kpc). The canonical value of moment of
# inertia of 1e38 kg m^2 is used
def spin_down_limit(freq, fdot, dist):
  hsd = math.sqrt((5./2.)*(G/C**3)*I38*math.abs(fdot)*freq)/(dist*KPC)
  
  return hsd
  
# Function to convert a pulsar stain into ellipticity assuming the canonical
# moment of inertia
def h0_to_ellipticity(h0, freq, dist):
  ell = h0*C**4*dist*KPC/(16.*math.pi**2*G*I38*freq**2)
  
  return ell

# function to convert the psi' and phi0' coordinates used in nested sampling
# into the standard psi and phi0 coordinates (using vectors of those parameters
def phipsiconvert(phipchain, psipchain):
  chainlen=len(phipchain)
  
  phichain = []
  psichain = []
  
  theta = math.atan2(1,2);
  ct = math.cos(theta);
  st = math.sin(theta);
  
  for i in range(0,chainlen):  
    phi0 = (1/(2*st))*phipchain[i] - (1/(2*st))*psipchain[i];
    psi = (1/(2*ct))*phipchain[i] + (1/(2*ct))*psipchain[i];
    
    # put psi between +/-pi/4
    if math.fabs(psi) > math.pi/4.:
      # shift phi0 by pi
      phi0 = phi0 + math.pi;
      
      # wrap around psi
      if psi > math.pi/4.:
        psi = -(math.pi/4.) + math.fmod(psi+(math.pi/4.), math.pi/2.);
      else:
        psi = (math.pi/4.) - math.fmod((math.pi/4.)-psi, math.pi/2.);
  
    # get phi0 into 0 -> 2pi range
    if phi0 > 2.*math.pi: 
      phi0 = math.fmod(phi0, 2.*math.pi);
    else:
      phi0 = 2.*math.pi - math.fmod(2.*math.pi-phi0, 2.*math.pi);
    
    phichain.append(phi0)
    psichain.append(psi)

  return phichain, psichain
 
# function to create histogram plot of the 1D posterior (potentially for
# multiple IFOs) for a parameter (param). If an upper limit is given then
# that will be output
def plot_posterior_hist(poslist, param, ifos,
                        parambounds=[float("-inf"), float("inf")], 
                        nbins=50, upperlimit=0):
  # create list of figures
  myfigs = []
  
  # create a list of upper limits
  ulvals = []
  
  # set some matplotlib defaults
  rc('text', usetex=True) # use LaTeX for all text
  rc('axes', linewidth=0.5) # set axes linewidths to 0.5
  rc('axes', grid=True) # add a grid
  rc('grid', linewidth=0.5)
  rc('font', family='serif')
  rc('font', size=12)
  
  # ifos line colour specs
  coldict = {'H1': 'b', 'H2': 'r', 'L1': 'g', 'V1': 'c', 'G1': 'm'}
  
  # some parameter names for special LaTeX treatment in figures
  paramdict = {'H0': '$h_0$', 'COSIOTA': '$\cos{\iota}$', 'PSI':
               '$\psi~{\rm(rads)}$', 'PHI0': '$\phi_0~{\rm(rads)}$', 'RA':
               '$\alpha$', 'DEC': '$\delta$'}
 
  # param name for axis label
  try:
    paraxis = paramdict[param.upper()]
  except:
    paraxis = param
  
  # loop over ifos
  for idx, ifo in enumerate(ifos):
    myfig = plt.figure(figsize=(4,3.5),dpi=200)
    
    pos = poslist[idx]
    
    pos_samps = pos[param].samples
 
    # get a normalised histogram for each
    n, bins = hist_norm_bounds( pos_samps, int(nbins), parambounds[0], \
                                parambounds[1] )
    
    # plot histogram
    plt.plot(bins, n, color=coldict[ifo])
    plt.xlabel(r''+paraxis, fontsize=14, fontweight=100)
    plt.ylabel(r'Probability Density', fontsize=14, fontweight=100)
    myfig.subplots_adjust(left=0.18, bottom=0.15) # adjust size
    
    myfigs.append(myfig)
    
    # if upper limit is needed then integrate posterior using trapezium rule
    if upperlimit != 0:
      ct = cumtrapz(n, bins)
      
      # prepend a zero to ct
      ct = np.insert(ct, 0, 0)
     
      #plt.plot(bincentres[1:len(bincentres)], ct)
      
      # use spline interpolation to find the value at 'upper limit'
      intf = interp1d(ct, bins, kind='cubic')
      ulvals.append(intf(float(upperlimit)))
  
  return myfigs, ulvals
  
# a function that creates and normalises a histograms of samples, with nbins
# between an upper and lower bound a upper and lower bound. The values at the
# bin points (with length nbins+2) will be returned as numpy arrays
def hist_norm_bounds(samples, nbins, low=float("-inf"), high=float("inf")):
  # get histogram
  n, binedges = np.histogram( samples, nbins )
  
  # get bin width
  binwidth = binedges[1] - binedges[0]
  
  # create bin centres
  bincentres = np.array([])
  for i in range(0, len(binedges)-1):
    bincentres = np.append(bincentres, binedges[i]+binwidth/2)
  
  # if histogram points are not close to boundaries (i.e. within a bin of the
  # boundaries) then add zeros to histrogram edges
  if bincentres[0] - binwidth > low:
    # prepend a zero to n
    n = np.insert(n, 0, 0);
    
    # prepend a new bin centre at bincentres[0] - binwidth
    bincentres = np.insert(bincentres, 0, bincentres[0] - binwidth)
  else:
    # we're  closer to the boundary edge than the bin width then, so set a new
    # bin on the boundary with a value linearly extrapolated from the
    # gradiant of the adjacent points
    dx = bincentres[0] - low;
    
    # prepend low value to bins
    bincentres = np.insert(bincentres, 0, low)
    
    dn = n[1]-n[0]
    
    nbound = n[0] - (dn/binwidth)*dx
    
    # prepend to n
    n = np.insert(n, 0, nbound)
  
  # now the other end!
  if bincentres[-1] + binwidth < high:
    # append a zero to n
    n = np.append(n, 0)
    
    # append a new bin centre at bincentres[end] + binwidth
    bincentres = np.append(bincentres, bincentres[-1] + binwidth)
  else:
    dx = high - bincentres[-1];
    
    # prepend low value to bins
    bincentres = np.append(bincentres, high)
    
    dn = n[-1]-n[-2]
    
    nbound = n[-1] + (dn/binwidth)*dx
    
    # prepend to n
    n = np.append(n, nbound)
  
  # now calculate area and normalise
  area = np.trapz(n, x=bincentres)
  
  ns = np.array([])
  for i in range(0, len(bincentres)):
    ns = np.append(ns, float(n[i])/area)
  
  return ns, bincentres

# create a Tukey window of length N
def tukey_window(N, alpha=0.5):
  # if alpha >= 1 just return a Hanning window
  if alpha >= 1:
    return np.hanning(N)
  
  # get x values at which to calculate window
  x = np.linspace(0, 1, N)
  
  # initial square window
  win = np.ones(x.shape)

  # get the left-hand side of the window  0 <= x < alpha/2  
  lhs = x<alpha/2
  win[lhs] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[lhs] - alpha/2) ))

  # get right hand side condition 1 - alpha / 2 <= x <= 1
  rhs = x>=(1 - alpha/2)
  win[rhs] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[rhs] - 1 + alpha/2)))
  
  return win

# create a function for plotting the absolute value of Bk data (read in from
# data files) and an averaged 1 day amplitude spectral density spectrogram for
# each IFO
def plot_Bks_ASDs( Bkdata, ifos ):
  # create list of figures
  Bkfigs = []
  psdfigs = []
  
  # ifos line colour specs
  coldict = {'H1': 'b', 'H2': 'r', 'L1': 'g', 'V1': 'c', 'G1': 'm'}
  
  # set some matplotlib defaults
  rc('text', usetex=True) # use LaTeX for all text
  rc('axes', linewidth=0.5) # set axes linewidths to 0.5
  rc('axes', grid=True) # add a grid
  rc('grid', linewidth=0.5)
  #rc('font', family='sans-serif')
  #rc('font', family='monospace')
  rc('font', family='serif')
  rc('font', size=12)
  
  # there should be data for each ifo
  for i, ifo in enumerate(ifos):
    # get data for given ifo
    try:
      dfile = open(Bkdata[i])
    except:
      print "Could not open file ", Bkdata[i]
      exit(-1)
      
    # should be three lines in file
    gpstime = []
    Bk = []
      
    Bkabs = [] # absolute Bk value
      
    # minimum time step between points (should generally be 60 seconds)
    mindt = float("inf")
      
    for line in dfile.readlines():
      sl = line.split()
        
      gpstime.append(float(sl[0]))
      Bk.append(complex(float(sl[1]), float(sl[2])))
        
      # get absolute value
      Bkabs.append(math.sqrt(Bk[-1].real**2 + Bk[-1].imag**2))
        
      # check time step
      if len(gpstime) > 1:
        dt = gpstime[-1] - gpstime[-2]
          
        if dt < mindt:
          mindt = dt
      
    dfile.close()
      
    # plot the time series of the data
    Bkfig = plt.figure(figsize=(10,4), dpi=200)
    Bkfig.subplots_adjust(bottom=0.12, left=0.06, right=0.94)
    
    tms = map(lambda x: x-gpstime[0], gpstime)
    
    plt.plot(tms, Bkabs, '.', color=coldict[ifo])
    plt.xlabel(r'GPS + ' + str(gpstime[0]), fontsize=14, fontweight=100)
    plt.ylabel(r'$|B_k|$', fontsize=14, fontweight=100)
    plt.title(r'$B_k$s for ' + ifo.upper(), fontsize=14)
    plt.xlim(tms[0], tms[-1])
    
    Bkfigs.append(Bkfig)
      
    # create PSD by splitting data into days, padding with zeros to give a
    # sample a second, getting the PSD for each day and combining them
    totlen = gpstime[-1] - gpstime[0] # total data length
      
    # check mindt is an integer and greater than 1
    if math.fmod(mindt, 1) != 0. or mindt < 1:
      print "Error time steps between data points must be integers"
      exit(-1)
        
    # loop over data, splitting it up
    mr = int(math.ceil(totlen/86400))
    
    count = 0
    npsds = 0
    totalpsd = np.zeros(86400) # add sum of PSDs
    for i in range(0, mr):
      datachunk = np.zeros(86400, dtype=complex)
       
      gpsstart = int(gpstime[count])
      
      prevcount = count
      
      for gt in gpstime[count:-1]:
        if gt >= gpsstart+86400:
          break
        else:
          datachunk[gt-gpsstart] = Bk[count]
          count += 1
      
      # only include the PSD if the chunk is more than 25% full of data
      pf = float(count-prevcount)*mindt/86400.
      
      if pf > 0.25:
        # get the PSD using a Tukey window with alpha = 0.25
        win = tukey_window(86400, alpha=0.25)
        
        Fs = 1 # sample rate in Hz
        
        psd, freqs, t = specgram(datachunk, NFFT=86400, Fs=Fs, window=win)
      
        # add psd onto total value
        totalpsd = map(lambda x, y: x+y, totalpsd, psd)
      
        # count number of psds
        npsds = npsds+1
      
    # average the PSD and convert to amplitude spectral density
    totalpsd = map(lambda x: math.sqrt(x/npsds), totalpsd)
    
    # plot PSD
    psdfig = plt.figure(figsize=(4,3.5), dpi=200)
    psdfig.subplots_adjust(left=0.18, bottom=0.15)
    
    # get the indices to plot in the actual frequency range 
    df = freqs[1]-freqs[0]
    minfbin = int((math.fabs(freqs[0])-1./(2.*mindt))/df)
    maxfbin = len(freqs) - minfbin
    
    plt.plot(freqs[minfbin:maxfbin], totalpsd[minfbin:maxfbin],
             color=coldict[ifo])
    plt.xlim(freqs[minfbin], freqs[maxfbin])
    plt.xlabel(r'Frequency (Hz)', fontsize=14, fontweight=100)
    plt.ylabel(r'$h/\sqrt{\rm Hz}$', fontsize=14, fontweight=100)
    plt.title(r'ASD for ' + ifo.upper(), fontsize=14)
    
    psdfigs.append(psdfig)
      
  return Bkfigs, psdfigs
  
import sys
import numpy as np
import matplotlib.pyplot as plt
from integrate_signal import *

BASE = '../../../Old/S-14/33-350/sims/single_switches/k1-5e4_multiseed/attempt'
ATTEMPT = 1

def plot_polar():
  attempt = sys.argv[1]
  t,mx,my,mz,hx,hy,hz = np.loadtxt('../../../Old/S-14/33-350/sims/experiment_matching/best_so_far/k3e4_a0-1_t2e-11_seed1.txt',unpack=True)
  #t,mx,my,mz,hx,hy,hz = np.loadtxt(BASE+attempt+'.txt',unpack=True)
  m_phi = np.arctan(mz/my)
  m_theta = np.arccos(mx)
  h_phi = np.arctan(hz/hy)
  h_r = np.sqrt(hx*hx+hy*hy+hz*hz)
  h_theta = np.arccos(hx/h_r)
  tau_r = (m_theta*h_r-m_phi*h_theta)
  tau_theta = (m_phi*h_r-h_phi)
  tau_phi = h_theta-m_theta*h_r
  plt.figure(1)
  plt.subplot(211);plt.hist(m_theta,100)#plt.plot(t,m_theta)
  plt.ylabel("$M_{\\theta}$")
  plt.subplot(212);plt.hist(m_phi,100)#plt.plot(t,m_phi)
  plt.ylabel("$M_{\phi}$")
  
  plt.figure(2)
  plt.subplot(311);plt.hist(tau_r,100)#plt.plot(t,integrate_signal(t,tau_r))
  plt.ylabel("$\int \\tau_r$")
  plt.subplot(312);plt.hist(tau_theta,100)#plt.plot(t,integrate_signal(t,tau_theta))
  plt.ylabel("$\int \\tau_{\\theta}$")
  plt.subplot(313);plt.hist(tau_phi,100)#plt.plot(t,integrate_signal(t,tau_phi))
  plt.ylabel("$\int \\tau_{\phi}$")
  plt.show()

def plot_integrated():
  attempt = sys.argv[1]
  #t,mx,my,mz,hx,hy,hz = np.loadtxt('../../../Old/S-14/33-350/sims/experiment_matching/best_so_far/k3e4_a0-1_t2e-11_seed1.txt',unpack=True)
  t,mx,my,mz,hx,hy,hz = np.loadtxt(BASE+attempt+'.txt',unpack=True)
  tau_x = []; tau_y = []; tau_z = []
  for i in xrange(len(t)):
    tau_x += [my[i]*hz[i]-mz[i]*hy[i]]
    tau_y += [mz[i]*hx[i]-mx[i]*hz[i]]
    tau_z += [mx[i]*hy[i]-my[i]*hx[i]]
  tau_phi = [np.arctan(tau_z[i]/tau_y[i]) for i in xrange(len(tau_z))]
  #tau_theta = [np.arccos(tau_x[i]) for i in xrange(len(tau_x))]
  plt.subplot(311)#;p1, = plt.plot(t,mx)#;p2, = plt.plot(t,my);p3, = plt.plot(t,mz);plt.xlabel("t (s)");plt.ylabel("arbitrary units")
  plt.hist(tau_x,50)
  plt.xlabel("t (s)")
  plt.ylabel("$M_x$")
  #plt.legend([p1,p2,p3],["$M_x$","$M_y$","M_z"])
  plt.subplot(312);plt.plot(t,integrate_signal(t,tau_x));plt.xlabel("t(s)");plt.ylabel("$\int \\tau_x$ (arbitrary units)")
  #plt.subplot(313);plt.plot(t,tau_x);plt.xlabel("t(s)");plt.ylabel("$\\tau_x$ (arbitrary units)")
  
  plt.subplot(313)
  #plt.plot(t,[np.sqrt(tau_y[i]**2+tau_z[i]**2) for i in xrange(len(tau_y))])
  #plt.plot(t,[np.sqrt(mz[i]*mz[i]+my[i]*my[i]) for i in xrange(len(t))])
  #plt.plot(t,integrate_signal(t,tau_phi))
  plt.hist(tau_phi,50)
  plt.xlabel("t(s)");plt.ylabel("$\int \sqrt{\\tau_y^2+\\tau_z^2}$ (arbitrary units)")
  #plt.subplot(413);plt.plot(t,integrate_signal(t,tau_y))
  #plt.subplot(414);plt.plot(t,integrate_signal(t,tau_z))
  plt.show()

def main():
  t,mx,my,mz,hx,hy,hz = np.loadtxt(sys.argv[1],unpack = True)
  region = slice(40000,60000)
  t = t[region]
  mx = mx[region]; my = my[region]; mz = mz[region]
  hx = hx[region]; hy = hy[region]; hz = hz[region]
  
  tau_x = []
  tau_y = []
  tau_z = []

  for i in xrange(len(t)):
    tau_x += [my[i]*hz[i]-mz[i]*hy[i]]
    tau_y += [mz[i]*hx[i]-mx[i]*hz[i]]
    tau_z += [mx[i]*hy[i]-my[i]*hx[i]]

  plt.subplot(211)
  plt.plot(t,mx)
  plt.subplot(212)
  plt.plot(t,[hx[i] for i in xrange(len(t))])
  #plt.subplot(413)
  #plt.plot(t[1:],[np.sqrt(hx[i]*hx[i]+hy[i]*hy[i]+hz[i]*hz[i])\
  #    -np.sqrt(hx[i-1]*hx[i-1]+hy[i-1]*hy[i-1]+hz[i-1]*hz[i-1])\
  #    for i in xrange(1,len(t))])
  #plt.subplot(412)
  #plt.plot(t,tau_x)
  #plt.subplot(413)
  #plt.plot(t,[tau_y[i] for i in xrange(len(tau_x))])
  #plt.subplot(414)
  #plt.plot(t,[tau_z[i] for i in xrange(len(tau_x))])

  plt.show()

plot_polar()
#plot_integrated()

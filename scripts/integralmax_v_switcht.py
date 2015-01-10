import numpy as np, matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from integrate_signal import *
from get_torque import *

BASE = '../../../Old/S-14/33-350/sims/single_switches/k1-5e4_multiseed/'

def line(x, *args):
  m,b = args
  return m*x+b

def main():
  max_ints = []
  max_ts = []
  for f in os.listdir(BASE):
    if f.endswith('.txt'):
      t = np.loadtxt(BASE+f,usecols=[0])
      tau_x,tau_y,tau_z = get_torque(BASE+f)
      int_tau = integrate_signal(t,tau_x)
      max_ints += [max(int_tau)]
      max_ts += [t[int_tau.index(max(int_tau))]]

  #fitted_params,_ = curve_fit(line, max_ts, max_ints, p0=[3.0,0.1])

  #plt.subplot(211);plt.plot(range(len(max_ints)),max_ints,linestyle='',marker='o')
  #plt.subplot(212);
  plt.plot(max_ts,max_ints,linestyle='',marker='o')
  plt.xlabel("t of max $\\tau_x$ (s)")
  plt.ylabel("max $\\tau_x$ (arbitrary units)")
  #plt.plot([max_ts[0],max_ts[-1]],[line(max_ts[0],*fitted_params),line(max_ts[-1],*fitted_params)])
  plt.show()

main()


      
      

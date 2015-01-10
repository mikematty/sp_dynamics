import math, os, sys, numpy as np, matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.misc import factorial
from scipy.stats import poisson
from scipy.special import gamma

def one_sim_distro(f_name):
  t,mx,my,mz,hx,hy,hz = np.loadtxt(f_name,unpack = True)
  run_len = 0.0
  tau_ups,tau_dns = [],[]
  
  dist_mx = [1 if m > 0 else -1 for m in mx]
  t_step = t[1]-t[0]

  for i in xrange(len(dist_mx)):
    if (i and (dist_mx[i] != dist_mx[i-1]) ):#and (run_len > 50)):
      if(dist_mx[i-1] > 0): tau_ups += [run_len*t_step]
      else: tau_dns += [run_len*t_step]
      run_len = 0.0
    run_len += 1.0

  return (tau_ups,tau_dns)
  

def main():
  tau_ups,tau_dns = [],[]
  if(len(sys.argv) > 1):
    for sim in os.listdir(sys.argv[1]):
      if(sim != ".DS_Store"):
        res = one_sim_distro(sys.argv[1]+sim)
        (tau_ups,tau_dns) = (tau_ups+res[0],tau_dns+res[1])

    outf_up = open("tau_ups.txt",'w')
    outf_dn = open("tau_dns.txt",'w')
    for t in tau_ups: outf_up.write(str(t)+'\n')
    for t in tau_dns: outf_dn.write(str(t)+'\n')
    outf_up.close()
    outf_dn.close()

  else:
    tau_ups = np.loadtxt("tau_ups.txt")
    tau_dns = np.loadtxt("tau_dns.txt")
  tau_ups *= 1000
  tau_dns *= 1000

  up_hist_bins,up_hist_edges = np.histogram(tau_ups,bins=50)
  dn_hist_bins,dn_hist_edges = np.histogram(tau_dns,bins=50)
  up_hist_x = [(up_hist_edges[i-1]+up_hist_edges[i])/2. for i in xrange(1,len(up_hist_edges))]
  dn_hist_x = [(dn_hist_edges[i-1]+dn_hist_edges[i])/2. for i in xrange(1,len(dn_hist_edges))]

  def poisson_fit(x, a, b, c):
    #a,l = args
    return a*np.exp(-b*x)+c
    #return a*np.power(l,x)*np.exp(-l)/(x*gamma(x))
  params_up, err_up = curve_fit(poisson_fit, np.array(up_hist_x),np.array(up_hist_bins), p0 = np.array([80.0,1.0,0.0]))
  params_dn, err_dn = curve_fit(poisson_fit, np.array(dn_hist_x),np.array(dn_hist_bins), p0 = np.array([100.0,1.0,0.0]))
  print "params_up", params_up
  print "err_up", err_up[0][0], err_up[1][1]
  print "chisq_up", sum([(x-poisson_fit(x,*params_up))**2/poisson_fit(x,*params_up) for x in up_hist_x[1:10]])
  print "params_dn", params_dn
  print "err_dn", err_dn[0][0], err_dn[1][1]

  plt.subplot(211)
  plt.hist(tau_ups,50)
  plt.plot(up_hist_x,[poisson_fit(x,*params_up) for x in up_hist_x],linewidth=3)
  plt.xlabel("$\\tau_{up}$ (ms)")
  plt.ylabel("Frequency")
  plt.title("Histogram of \"up\" Lifetimes")

  plt.subplot(212);plt.hist(tau_dns,50)
  plt.plot(dn_hist_x,[poisson_fit(x,*params_dn) for x in dn_hist_x],linewidth=3)
  plt.xlabel("$\\tau_{dn}$ (ms)")
  plt.ylabel("Frequency")
  plt.title("Histogram of \"down\" Lifetimes")
  
  plt.suptitle("All Events")
  plt.show()

main()

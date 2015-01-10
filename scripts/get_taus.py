import matplotlib.pyplot as plt
import numpy as np
import sys
import os

CUTOFF = 0.0

def taus(file_name):
  t, mx = np.loadtxt(file_name, usecols = [0,1], unpack = True)
  #print t[1]-t[0]

  run_len = 0.
  up_total = 0.
  up_runs = 0.
  dn_total = 0.
  dn_runs = 0.
  
  dist_mx = [1 if m > CUTOFF else -1 for m in mx]
  
  for i in xrange(len(dist_mx)):
    if(i and (dist_mx[i] != dist_mx[i-1])):
      if(dist_mx[i-1] > 0):
        up_total += run_len
        up_runs += 1
      else:
        dn_total += run_len
        dn_runs += 1
      run_len = 0
    run_len += 1

  try: t_up = (up_total/up_runs)*(t[1]-t[0])
  except: t_up = 0; t_dn = t[-1]-t[0]
  try:t_dn = (dn_total/dn_runs)*(t[1]-t[0])
  except: t_up = t[-1]-t[0]; t_dn = 0
  return (t_up,t_dn)

def parse_filename(file_name):
  file_name = file_name[file_name.rfind('/')+1:]
  k = float(file_name[file_name.index('k')+1:file_name.find('_')].replace('-','.'))
  a = float(file_name[file_name.index('a')+1:file_name.rfind('_')].replace('-','.'))
  t = file_name[file_name.index('t')+1:file_name.find('.')]
  (num,unit) = t.split('m')
  if (unit == 'us'): t = float(num)*1e-6
  if (unit == 's'): t = float(t)*1e-3
  return(k,a,t)

def parse_filename2(file_name):
  file_name = file_name[file_name.rfind('/')+1:]
  k = float(file_name[file_name.index('k')+1:file_name.find('_')].replace('-','.'))
  a = float(file_name[file_name.index('a')+1:file_name.rfind('_')].replace('-','.'))
  tstep = file_name[file_name.index('t')+1:file_name.find('.')]
  if(tstep.count('-') > 1): tstep = tstep.replace('-','.',1)
  tstep = float(tstep)
  return(k,a,tstep)

def directory_data2(path,out_name):
  out_file = open(out_name, 'w')
  out_file.write('K a step t_up t_dn\n')
  files = os.listdir(path)
  for f in files:
    (k,a,step) = parse_filename2(path+f)
    (t_up,t_dn) = taus(path+f)
    out_file.write(str(k)+' '+str(a)+' '+str(step)+' '+\
                   str(t_up)+' '+str(t_dn)+'\n')
  out_file.close()

def directory_data(path,out_name):
  out_file = open(out_name, 'w')
  out_file.write('K a len t_up t_dn\n')
  files = os.listdir(path)
  for f in files:
    (k,a,run_len) = parse_filename(path+f)
    if k < 16000:
      (t_up,t_dn) = taus(path+f)
      out_file.write(str(k)+' '+str(a)+' '+str(run_len)+' '+\
                   str(t_up)+' '+str(t_dn)+'\n')

  out_file.close()

def plot_k_scaling(out_name):
  k,a,run_len,t_up,t_dn = np.loadtxt(out_name,skiprows = 1, unpack = True)
  p1, = plt.plot(k,t_up*10**6, linestyle = '', marker = 'o', markersize = 7)
  p2, = plt.plot(k,t_dn*10**6, linestyle = '', marker = 'o', markersize = 7)
  plt.legend([p1,p2],["$\\tau_{+}$","$\\tau_{-}$"],loc=4)
  plt.xlabel("K $(J/m^3)$")
  plt.ylabel("$\\tau_i (\mu s)$")

  #plt.show()

def plot_a_scaling(out_name):
  k,a,run_len,t_up,t_dn = np.loadtxt(out_name, skiprows = 1, unpack = True)
  p1, = plt.plot(a,t_up*10**6, linestyle = '', marker = 'o', markersize = 7)
  p2, = plt.plot(a,t_dn*10**6, linestyle = '', marker = 'o', markersize = 7)
  #p3, = plt.plot(a,[1e-09 for i in a],linestyle='--')
  plt.legend([p1,p2],["$\\tau_+$","$\\tau_-$"])#,"$\Delta\\tau$"])
  plt.xlabel("$\\alpha$")
  plt.ylabel("$\\tau_i (\mu s)$")

  #plt.show()

def plot_t_scaling(out_name):
  k,a,tstep,t_up,t_dn = np.loadtxt(out_name, skiprows = 1, unpack = True)
  p1, = plt.plot(tstep, t_up*10**6, linestyle = '', marker = 'o', markersize = 7)
  p2, = plt.plot(tstep, t_dn*10**6, linestyle = '', marker = 'o', markersize = 7)
  plt.legend([p1,p2],["$\\tau_+$","$\\tau_-$"],loc=4)
  plt.xlabel("$\Delta t$")
  plt.ylabel("$\\tau_i (\mu s)$")

#directory_data('../sims/a_scaling/','../sims/a_scaling_out.txt')
#plot_k_scaling('../sims/k_scaling_out.txt')
#plt.figure(2)
#plot_a_scaling('../sims/a_scaling_out.txt')

#directory_data2('../sims/tstep_scaling/','../sims/tstep_scaling_out.txt')
#plot_t_scaling('../sims/scaling/tstep_scaling_out.txt')
#plt.show()
print taus(sys.argv[1])

#def main():
#  (t_up,t_down) = taus(sys.argv[1])
#  print "t_up = ",t_up, "t_down = ",t_down

#if __name__ == "__main__": main()

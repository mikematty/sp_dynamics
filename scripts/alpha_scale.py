import os

def main():
  for a in [0.01,.015,0.02,.025,0.03,.035,0.04,.045,0.05,.055,\
      0.06,0.07,0.08,0.09,0.1]:
    f_name = 'v3e-24_t1e-14_a'+str(a).replace('.','-')+'.txt'
    os.system('../simulation_src/run '+f_name+' '+str(a))

main()

import numpy as np, matplotlib.pyplot as plt

def get_torque(inf):
  t,mx,my,mz,hx,hy,hz = np.loadtxt(inf,unpack = True)
  tau_x,tau_y,tau_z = [],[],[]
  for i in xrange(len(t)):
    tau_x += [my[i]*hz[i]-mz[i]*hy[i]]
    tau_y += [mz[i]*hx[i]-mx[i]*hz[i]]
    tau_z += [mx[i]*hy[i]-my[i]*hy[i]]
  return tau_x,tau_y,tau_z

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import os,subprocess

chunk = slice(48000,56000)
t,mx,my,mz,hx,hy,hz = np.loadtxt("../../../../Old/S-14/33-350/sims/single_switches/k1-5e4_multiseed/attempt1.txt",unpack = True)
mx = mx[chunk]
my = my[chunk]
mz = mz[chunk]

class Arrow3D(FancyArrowPatch):
  def __init__(self, xs, ys, zs, *args, **kwargs):
    FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
    self._verts3d = xs, ys, zs

  def draw(self, renderer):
    xs3d, ys3d, zs3d = self._verts3d
    xs, ys, zs = proj3d.proj_transform(xs3d,ys3d,zs3d,renderer.M)
    self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
    FancyArrowPatch.draw(self, renderer)

def plot_frame(i,label):
  fig = plt.figure(figsize = (8,8))
  ax = Axes3D(fig)
  ax.view_init(elev = 10, azim = 250)
  ax.scatter(mx,my,mz,c = mx)
  vec = Arrow3D([0,mx[i]],[0,my[i]],[0,mz[i]], mutation_scale = 20, lw = 3)#, arrowstyle = "-|>")
  ax.add_artist(vec)
  plt.savefig('frame'+label)

def main():
  label = '000'
  for i in xrange(0,len(mx),100):
    plot_frame(i,label)
    if(int(label) < 9): label = '00'+str(int(label)+1)
    else: label = '0'+str(int(label)+1)
    #label = str(int(label)+1)
  #subprocess.call(['ffmpeg', '-framerate','15','-i', 'frame%d.png','output.avi'])
  #subprocess.call(['ffmpeg', '-i', 'output.avi', '-t', '20', '-r', '40', 'out.gif'])

main()


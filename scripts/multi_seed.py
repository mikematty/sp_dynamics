import os

def main():
  for i in xrange(1,21):
    os.system('./run k3e4_a0-1_t2e-11_seed'+str(i)+'.txt '+str(i))
  return

main()

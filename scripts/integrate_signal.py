def integrate_signal(t,s):
  for i in xrange(len(s)):
    if(i == 0): integral = [0]
    else: integral += [integral[i-1]+0.5*(t[1]-t[0])*(s[i]+s[i-1])]
  return integral

import numpy as np
import scipy.integrate as integrate
import re

f=open("./g2parameters.h","r")
fullFile=f.read()
f.close()

vars={}
for line in fullFile.split('\n'):
    if line.endswith('PyVAR'):
        temp=re.split("\s|;",line)
        vars[temp[2]]=float(temp[4])
nn=vars['N']
ll=vars['L']
op=vars['omega']
rs=vars['rstar']
a=vars['alpha']
dx=ll/nn
        
intLimit=integrate.quad(lambda r: (np.sqrt(1 + a)*op*np.sqrt(np.pi)*(-3*r**2 + np.sqrt(9*r**4 + (32*r*rs**3)/np.pi)))/(8.*r*rs**3), 0, np.inf)
intVal=integrate.quad(lambda r: (np.sqrt(1 + a)*op*np.sqrt(np.pi)*(-3*r**2 + np.sqrt(9*r**4 + (32*r*rs**3)/np.pi)))/(8.*r*rs**3), 0, (nn-nn/2.-1./np.sqrt(2.))*dx)
bkgFld=intVal[0]-intLimit[0]
old=[ v for v in fullFile.split("\n") if v.endswith('PyWRITE') ][0]
new=old.replace(old.split()[4],str(bkgFld)+';')
newFile=fullFile.replace(old,new)

f=open("./g2parameters.h","w")
f.write(newFile)
f.close()
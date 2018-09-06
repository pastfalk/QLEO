#Plot dispersion relation and magnetic energy spectrum at given time step
import numpy
import pylab

step=87

file_name="evolution/omega/"+numpy.str(step)+".dat"

t = open(file_name, "r")
num_lines = sum(1 for line in t)
t.seek(0)

print(num_lines)

t.readline()

k=numpy.zeros(num_lines)
Bksq=numpy.zeros(num_lines)
omega=numpy.zeros(num_lines)
gamma=numpy.zeros(num_lines)

i=0

for line in t:

	numbers_str = line.split()
		
	k[i]=numpy.float(numbers_str[0])
	Bksq[i]=numpy.float(numbers_str[1])
	omega[i]=numpy.float(numbers_str[2])
	gamma[i]=numpy.float(numbers_str[3])
	
	i=i+1

t.close()

pylab.subplots(3,1,figsize=(7,12))
pylab.subplot(3,1,1)
pylab.plot(k,omega,'-r',linewidth=4.0)
pylab.xlabel('$\\tilde{k}$')
pylab.ylabel('$\\tilde{\omega}$')
pylab.subplot(3,1,2)
pylab.plot(k,gamma,'-b',linewidth=4.0)
pylab.xlabel('$\\tilde{k}$')
pylab.ylabel('$\\tilde{\gamma}$')
pylab.subplot(3,1,3)
pylab.plot(k,Bksq,'-g',linewidth=4.0)
pylab.xlabel('$\\tilde{k}$')
pylab.ylabel('$\\tilde{B}_k ^2$')

pylab.show()



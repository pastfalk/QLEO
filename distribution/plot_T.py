#Plot solutions of wavenumber integration for given particle species and at given time step
import numpy
import pylab

species="1"

step=87

nvpa=128


#file_name="evolution/T1/"+numpy.str(step)+"_"+species+".dat"
#icol=1 # T11 =-0.25*dt*(q*mu)**2 * dB_k^2 * |omega|^2 /k^2 /(omega-k*vpa-sgn*mu*q)
#icol=2 # T12 =-0.25*dt*(q*mu)**2 * dB_k^2 * omega^* /k /(omega-k*vpa-sgn*mu*q)
#icol=3 # T13 =-0.25*dt*(q*mu)**2 * dB_k^2 * omega/k /(omega-k*vpa-sgn*mu*q)
#icol=4 # T14 =-0.25*dt*(q*mu)**2 * dB_k^2 /(omega-k*vpa-sgn*mu*q)

file_name="evolution/T2/"+numpy.str(step)+"_"+species+".dat"
#icol=1 # T21 = -0.25*dt*(q*mu)**2 * dB_k^2 * omega /(omega-k*vpa-sgn*mu*q)^2
icol=2 # T22 = -0.25*dt*(q*mu)**2 * dB_k^2 * k /(omega-k*vpa-sgn*mu*q)^2


t = open(file_name, "r")
num_lines = sum(1 for line in t)
t.seek(0)

print(num_lines)
print(nvpa)

dist=numpy.zeros(nvpa)

vpara=numpy.zeros(nvpa)

i=0

for line in t:

	numbers_str = line.split()
		
	numbers_float=numpy.zeros(2)

	numbers_float[0]=numpy.float(numbers_str[0])
	numbers_float[1]=numpy.float(numbers_str[icol])

	vpara[i]=numbers_float[0]
	dist[i]=numbers_float[1]
	
	i=i+1

t.close()


pylab.rcParams['contour.negative_linestyle'] = 'solid'

pylab.plot(vpara,dist,'+r')
pylab.xlabel('$\ttilde{v}_\parallel$')
pylab.ylabel('$T$')

pylab.show()



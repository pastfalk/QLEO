#Plot parallel or perpendicular velocity derivative of distribution for given particle species and at given time step. 
import numpy
import pylab

species="1"

step=87

nvpa=128
nvpe=64


file_name="evolution/delfdelpa/"+numpy.str(step)+"_"+species+".dat"
icol=1 # d f / d vpara
#icol=2 # d^2 f / d vpara^2
#icol=3 # d^2 f / d vpara d vperp

#file_name="evolution/delfdelpe/"+numpy.str(step)+"_"+species+".dat"
#icol=1 # d f / d vperp
#icol=2 # d^2 f / d vperp^2


t1 = open(file_name, "r")
num_lines1 = sum(1 for line in t1)
t1.seek(0)

print(num_lines1)
print(nvpe*nvpa)

dist=numpy.zeros((nvpe,nvpa))

vpara=numpy.zeros(nvpa)
vperp=numpy.zeros(nvpe)

i=0
ivpara=-1
ivperp=0

for line in t1:

	numbers_str = line.split()

	numbers_float=numpy.zeros(3)

	numbers_float[0]=numpy.float(numbers_str[0])
	numbers_float[1]=numpy.float(numbers_str[1])
	numbers_float[2]=numpy.float(numbers_str[icol+1])

	if(ivpara==-1):
		ivpara=ivpara+1
		vpara[ivpara] = numbers_float[0]

	if (numbers_float[0] == vpara[ivpara]):
		vperp[ivperp] = numbers_float[1]
		dist[ivperp,ivpara] = numbers_float[2]


		ivperp=ivperp+1

	else:
		ivpara=ivpara+1
		vpara[ivpara]=numbers_float[0]
		ivperp=0

		vperp[ivperp] = numbers_float[1]
		dist[ivperp,ivpara] = numbers_float[2]

		ivperp=ivperp+1
		
	i=i+1
		
	if(i > nvpa*nvpe):
		i=0

t1.close()



levels =numpy.linspace(-41,-2,num=40)
pylab.contourf(vpara, vperp, numpy.log10( abs(dist) ),levels,cmap=pylab.cm.get_cmap("gist_heat"))
pylab.colorbar()
pylab.rcParams['contour.negative_linestyle'] = 'solid'
pylab.contour( vpara, vperp, numpy.log10( abs(dist)), levels,colors='k')
pylab.xlabel('$\\tilde{v}_\parallel$',fontsize=20)
pylab.ylabel('$\\tilde{v}_\perp$',fontsize=20)
pylab.show()



#Plot velocity distribution for given particle species and at given time step, and overplot with bi-Maxwellian of given beta parameters
import numpy
import pylab

species="1"

step=87

nvpa=255
nvpe=64

beta_para=4.0
beta_perp=2.0


file_name="evolution/dist/"+numpy.str(step)+"_"+species+".dat"

t = open(file_name, "r")
num_lines = sum(1 for line in t)
t.seek(0)

print(num_lines)
print(nvpe*nvpa)

dist=numpy.zeros((nvpe,nvpa))

vpara=numpy.zeros(nvpa)
vperp=numpy.zeros(nvpe)

i=0
ivpara=-1
ivperp=0

for line in t:

	numbers_str = line.split()
		
	numbers_float=numpy.zeros(3)

	numbers_float[0]=numpy.float(numbers_str[0])
	numbers_float[1]=numpy.float(numbers_str[1])
	numbers_float[2]=numpy.float(numbers_str[2])

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

t.close()




maxwell=numpy.zeros((nvpe,nvpa))

for i in range(0,nvpa):
	for j in range(0,nvpe):
		
		fore_MW= numpy.pi**(3.0/2.0) *numpy.sqrt(beta_para)*beta_perp
		maxwell[j,i]=numpy.exp(-vpara[i]**2 /beta_para) * numpy.exp(-vperp[j]**2 /beta_perp)/fore_MW


levels1 =numpy.linspace(-41,-2,num=40)
pylab.contourf(vpara, vperp, numpy.log10( dist ),levels1,cmap=pylab.cm.get_cmap("gist_heat"))
pylab.colorbar()
pylab.contour( vpara, vperp, numpy.log10(maxwell), levels1,colors='k')
pylab.rcParams['contour.negative_linestyle'] = 'solid'
pylab.contour( vpara, vperp, numpy.log10(dist), levels1,colors='k')
pylab.xlabel('$ \\tilde{v}_\parallel $',fontsize=20)
pylab.ylabel('$\\tilde{v}_\perp$',fontsize=20)

pylab.show()



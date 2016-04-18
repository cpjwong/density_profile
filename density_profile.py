import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from scipy.optimize import curve_fit
import math

def den_prof(z,a,b):
	tt = 1-(0.5+0.5*np.tanh((z-a)/(2*b)))
	return(tt)
def first_d(p,z,em,r,k,T):
	f_d = k*T*(r+np.log(p)*r-1-np.log(1-p))+z*em*p-0.5*z*em
	return(f_d)
T = [2.5,3,3.5,4]
T_label = ['$T=4.0$','$T=4.5$','$T=5.0$','$T=5.5$']
color = ['-b','-r','-g','-m']
gas_s = np.zeros(4)
gas_y = np.zeros(4)
for ss in range(len(T)):
	n = np.arange(0.01,1,0.01)
	z = 14
	em = -1
	k = 1
	m = (0.5/2)**0.5
	r = 1/200
	density = np.zeros(len(n))
	work = np.zeros(len(n))
	second_derivative = np.zeros(len(n))
	for i in range(len(n)):
		density[i] = k*T[ss]*(n[i]*r*np.log(n[i])+(1-n[i])*np.log(1-n[i]))+0.5*z*em*n[i]*(n[i]-1)
		second_derivative[i] = first_d(n[i],z,em,r,k,T[ss])
		if abs(n[i] - 0.4) < 0.001:
			dw = density[i]
		else:
			pass
	lower = np.zeros(int(len(density)))
	for i in range(int(len(density))):
		lower[i] = second_derivative[i]
	count = 0
	px = []
	py = []
	tro = 0
	for j in range(len(lower)):
		if j!= 0 and abs(lower[j] - 0) < 3.5*10**(-1) and lower[j]-lower[j-1] > 0:
			px.append(float(n[j])-0.03)
			py.append(float(density[j]))
			gas_s[ss] = n[j]
			gas_y[ss] = density[j]
			count = count + 1
		
		else:
			pass
	px.append(0)
	py.append(0)
	print(lower)
	plt.figure
	plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
	plt.rc('text', usetex=True)
	plt.gcf().set_size_inches(5,3,forward=True)
	plt.subplots_adjust(left=0.18,bottom=0.15)
	plt.plot(n,density,'-m')
	plt.plot(px,py,'-b')
	jw = (0-gas_y[ss])*0.4/(0-gas_s[ss])
	plt.annotate('', xy=(0.4,jw), xycoords='data',xytext=(0.4,dw),arrowprops={'arrowstyle': '<->'})#textcoords='data'
	#plt.plot((0.4,0.4),(jw,jw), 'k', marker='v',) # lower arrowhead
	#plt.plot((0.4,0.4),(dw,dw), 'k', marker='^',) # upper arrowhead
	plt.text(0.45,0.2,r'$-W(\rho)$')
	plt.xlabel(r'$\rho$')
	plt.ylabel('$\psi$')
	if ss == 2:
		plt.savefig('free_energy_200_1.png',dpi=300)
	else:
		pass
	plt.show()
print(gas_s)


T = [2.5,3,3.5,4]
#gas_s = [0.999,0.957996,0.93566,0.871727]#0.916389,0.871727,0.73699
n = np.arange(0.1,1,0.01)
T_label = ['$T=2.5$','$T=3.0$','$T=3.5$','$T=4.0$']
color = ['-b','-r','-g','-m']
z2 = 14
L = 0.96
cor = 10**(-3)
plt.figure
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.gcf().set_size_inches(5,3,forward=True)
plt.subplots_adjust(left=0.18,bottom=0.15)
for ii in range(len(T)):
	n = np.arange(0.01,gas_s[ii]+0.01,0.01)
	z = np.zeros(len(n))
	for i in range(len(z)):
		h = (n[i]-0.01)/100
		n2 = np.arange(0.01,n[i]+h,h)
		su = 0
		for j in range(len(n2)):
			L = gas_s[ii]
			#g = ((k*T[ii]*(n2[j]*r*np.log(n2[j])+(1-n2[j])*np.log(1-n2[j]))+0.5*z2*em*n2[j]*(n2[j]-1))-(k*T[ii]*(L*r*np.log(L)+(1-L)*np.log(1-L))+0.5*z2*em*L*(L-1)))**(-0.5)*(m)**0.5
			g = ((k*T[ii]*(n2[j]*r*np.log(n2[j])+(1-n2[j])*np.log(1-n2[j]))+0.5*z2*em*n2[j]*(n2[j]-1))- (0-gas_y[ii])*n2[j]/(0-gas_s[ii]))**(-0.5)*(m)**0.5
			if j == 0 or j == (len(n2)-1):
				su = su + g
			elif j%2 == 0:
				su = su+2*g
			else:
				su = su+4*g
		z[i] = -(h/3)*su
	
	z22 = []
	nn = []
	for i in range(len(z)):
		if math.isnan(z[i]) == True or abs(z[i]) == float('Inf'):
			pass
		else:
			z22.append(float(z[i]))
			nn.append(float(n[i]))
	#print(z22)
	norm_n = np.zeros(len(nn))
	#xx2 = np.sort(nn)
	#xx2 = nn[int(len(nn)-1)]
	for i in range(len(nn)):
		norm_n[i] = (nn[i]-0)/(gas_s[ii]-0)
	popt, pcov = curve_fit(den_prof,z22,norm_n)
	z23 = np.arange(-2,1.0,0.01)
	density1 = np.zeros(len(z23))
	for i in range(len(density1)):
		density1[i] = 0.5-0.5*np.tanh((z23[i]-popt[0])/(2*popt[1]))
	#if ii == 1:
	#	plt.plot(z22,norm_n,'*r')#label=T_label[ii])
	#	plt.plot(z23,density1,color[ii])#,label=T_label[ii])
	#else:
	#
	plt.plot(z23,density1,color[ii],label=T_label[ii])
	#plt.plot(z22,norm_n,'*r')
	print(popt[1])
plt.xlim(-1.5,1)
plt.legend(loc=1,fontsize=8)
plt.xlabel('$z$')
plt.ylabel(r'$\frac{\rho-\rho_{gas}}{\rho_{liquid}-\rho_{gas}}$')
#plt.xlim(-1.5,1)
plt.savefig('density_profile_PS.png',dpi=300)
plt.show()

# density_profile

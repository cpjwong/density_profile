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
r = [1/10,1/50,1/100,1/200]
r = np.arange(1/200,1/10,0.001)
T = 2.5
T_label = ['$T=4.0$','$T=4.5$','$T=5.0$','$T=5.5$']
color = ['-b','-r','-g','-m']
gas_s = np.zeros(len(r))
gas_y = np.zeros(len(r))
for ss in range(len(r)):
	n = np.arange(0.01,1,0.01)
	z = 14
	em = -1
	k = 1
	m = (0.5/2)**0.5
	density = np.zeros(len(n))
	work = np.zeros(len(n))
	second_derivative = np.zeros(len(n))
	for i in range(len(n)):
		density[i] = k*T*(n[i]*r[ss]*np.log(n[i])+(1-n[i])*np.log(1-n[i]))+0.5*z*em*n[i]*(n[i]-1)
		second_derivative[i] = first_d(n[i],z,em,r[ss],k,T)
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
	jw = (0-gas_y[ss])*0.4/(0-gas_s[ss])

n = np.arange(0.1,1,0.01)
T_label = ['$T=2.5$','$T=3.0$','$T=3.5$','$T=4.0$']
color = ['-b','-r','-g','-m']
z2 = 14
cor = 10**(-3)
for ii in range(len(r)):
	n = np.arange(0.01,gas_s[ii]+0.01,0.01)
	z = np.zeros(len(n))
	for i in range(len(z)):
		h = (n[i]-0.01)/100
		n2 = np.arange(0.01,n[i]+h,h)
		su = 0
		for j in range(len(n2)):
			L = gas_s[ii]
			g = ((k*T*(n2[j]*r[ii]*np.log(n2[j])+(1-n2[j])*np.log(1-n2[j]))+0.5*z2*em*n2[j]*(n2[j]-1))- (0-gas_y[ii])*n2[j]/(0-gas_s[ii]))**(-0.5)*(m)**0.5
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
	norm_n = np.zeros(len(nn))
	for i in range(len(nn)):
		norm_n[i] = (nn[i]-0)/(gas_s[ii]-0)
	popt, pcov = curve_fit(den_prof,z22,norm_n)
	z23 = np.arange(-2,1.0,0.01)
	density1 = np.zeros(len(z23))
	for i in range(len(density1)):
		density1[i] = 0.5-0.5*np.tanh((z23[i]-popt[0])/(2*popt[1]))

z2 = 14
cor = 10**(-3)
gamma = np.zeros(len(r))
for ii in range(len(r)):
	h = (0.86-0.01)/100
	n2 = np.arange(0.01,0.86+h,h)
	su = 0
	for j in range(len(n2)):
		g = ((k*T*(n2[j]*r[ii]*np.log(n2[j])+(1-n2[j])*np.log(1-n2[j]))+0.5*z2*em*n2[j]*(n2[j]-1))- (0-gas_y[ii])*n2[j]/(0-gas_s[ii]))**(0.5)*(2*m)**0.5
		if j == 0 or j == (len(n2)-1):
			su = su + g
		elif j%2 == 0:
			su = su+2*g
		else:
			su = su+4*g
	gamma[ii] = (h/3)*su
plt.figure
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.gcf().set_size_inches(5,3,forward=True)
plt.subplots_adjust(left=0.18,bottom=0.15)
r = [1/ri for ri in r]
plt.plot(r,gamma,'-*g')
plt.xlabel('$1/r$')
plt.ylabel('$\gamma$')
plt.savefig('surf_ten_r_3.png',dpi=300)
plt.show()

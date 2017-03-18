import numpy as np
import matplotlib.pyplot as plt
import math

#x_puntos =np.array([3,5,12])
x_puntos =np.array([1.5,1.7,2.0])
a=1
b=20
lmin=0.01
lmax=5.0

def Evidencia(datos):
	return 0.02


def P_verosim(datos, lamb,zet):
	prod =1.0
	for i in datos:
		prod=prod*math.exp(-i/lamb)/zet
	return prod

def Zeta_func(lamb,a,b):
	return -1.0*lamb*(math.exp(-b/lamb)-math.exp(-a/lamb))

def Prior(x,mini,maxi):
	return 1.0/(maxi-mini);
	
lambdas=np.linspace(lmin,lmax,200)

posterior=lambdas.copy()

for i in range(len(lambdas)):
	zeta=Zeta_func(lambdas[i],a,b)
	vero=P_verosim(x_puntos,lambdas[i],zeta)
	p=Prior(lambdas[i],lmin,lmax)
	Ev=Evidencia(x_puntos)
	posterior[i]=p*vero/Ev

def get_prob_lamb(lambdi):
	zeta=Zeta_func(lambdi,a,b)
	vero=P_verosim(x_puntos,lambdi,zeta)
	p=Prior(lambdi,lmin,lmax)
	Ev=Evidencia(x_puntos)
	return p*vero/Ev

def proba_dist(xs,distr,value):
	return sum(distr[xs<=value])

longitud=20000
sigma=0.1
histo_l=np.zeros([longitud])

paso=0.01
histo_l[0]=0.2
for i in range(1,len(histo_l)):
	next_p=np.random.normal(histo_l[i-1],sigma)
	alpha=get_prob_lamb(next_p)/get_prob_lamb(histo_l[i-1])
	if (alpha>1.0):
		histo_l[i]=next_p
	else:
		beta=np.random.random()
		if (beta<=alpha):
			histo_l[i]=next_p
		else:
			histo_l[i]=histo_l[i-1]

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(lambdas,posterior)
#ax.hist(histo_l, normed=1)
plt.show()

import numpy as np
import matplotlib.pyplot as plt

def Lax(U,dt,dx,tms):
	U_fin=np.zeros(shape(U)[0],shape(tms)[0])
	t=0
	for i in range(tms):
		while t<=tms[i]:
			U_fin[:,i]=U.copy()
			F=np.power(U,2)/2
			U_fin[1:-1,i]= (U[2:,i]-U[0:-2,i])/2-dt*(F[2:]-F[0:-2])/(2*dx)
			t=t+dt
			U=U_fin.copy()
	return U_fin

dx=0.05
dt=dx*0.5
x=linspace(0,4,dx)
U_ini=zeros(shape(x))
U_ini[x<=2]=1


import numpy as np
import matplotlib.pyplot as plt


def FTFS(U,a,tm,dx,dt):
	t=0
	U_new = np.zeros([61]);
	while t<=tm:
		for i in range(np.shape(U_ini)[0]-1):
			U_new[i]=U[i]-a*dt*(U[i+1]-U[i])/dx
		U=U_new.copy()
		t=t+dt
	return U

def FTCS(U,a,tm,dx,dt):
	t=0
	U_new = np.zeros([61]);
	while t<=tm:
		for i in range(1,np.shape(U_ini)[0]-1):
			U_new[i]=U[i]-a*dt*(U[i+1]-U[i-1])/(2*dx)
		U=U_new.copy()
		t=t+dt
	return U

a=-20
tmax=0.45
dx=5
dt=0.015

x = dx*np.array(range(0,12))
temp=100*np.sin(np.pi*x/60)
U_ini=np.zeros([61])
U_ini[50/dx:110/dx]=temp

x=dx*np.array(range(0,61))
U_final=FTCS(U_ini,a,tmax,dx,dt)
plt.plot(x,U_ini)
plt.plot(x,U_final)
plt.show()


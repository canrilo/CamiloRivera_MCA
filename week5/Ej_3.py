import numpy as np
import matplotlib.pyplot as plt

#Solucion de Mac Cormak para ecuaciones de Euler en zonas de alta y baja presion

def Mac_Cormack(gamma,U,F,dx,dt,tsteps):
	U_temp=U.copy()
	F_temp=F.copy()
	
	U_new=U.copy()
	F_new=F.copy()
	
	for i in range(tsteps):
		U_temp[:,1:-1]=U[:,1:-1]-dt*(F[:,2:]-F[:,1:-1])/dx
		F_temp=get_F(U_temp,gamma)
		U_new[:,1:-1]=(U[:,1:-1]+U_temp[:,1:-1]-dt*(F_temp[:,1:-1]-F_temp[:,0:-2])/dx)/2
		F_new= get_F(U_new,gamma)
		
		F=F_new.copy()
		U=U_new.copy()
	return U

def get_F(U,gamma):
	rho=U[0,:]
	u=U[1,:]/rho
	e=U[2,:]
	p=(gamma-1)*(e-rho*(u**2)/2.0)
	F=np.zeros(np.shape(U))
	F[0,:]=u*rho
	F[1,:]=rho*(u**2)+p
	F[2,:]=u*(e+p)
	return F

def get_e(gamma, rho,u,p):
	u2=u**2
	e=p/(gamma-1)+rho*u2/(2.0)
	return e

dx=0.1
dt=0.05
tmax=1.0
xsteps=int(1/dx)+1
tsteps=int(tmax/dt)+1
gamma=1.4
x=np.linspace(0,1,xsteps)
u=np.zeros(np.shape(x))

rho=u.copy()
rho[x<=0.5]=1.0
rho[x>0.5]=0.125

p=u.copy()
p[x<=0.5]=1.0
p[x>0.5]=0.1

e=get_e(gamma,rho,u,p)

U_ini=np.zeros([3,np.shape(x)[0]])
U_ini[0,:]=rho.copy()
U_ini[1,:]=rho*u
U_ini[2,:]=e.copy()

F_ini=np.zeros([3,np.shape(x)[0]])
F_ini[0,:]=u*rho
F_ini[1,:]=rho*(u**2)+p
F_ini[2,:]=u*(e+p)

U_final=Mac_Cormack(gamma,U_ini,F_ini,dx,dt,tsteps)

rho_final =U_final[0,:]
u_final = U_final[1,:]/rho_final
e_final = U_final[2,:]
p_final=(gamma-1)*(e_final-rho_final*(u_final**2)/2.0)

fig = plt.figure()

axu = fig.add_subplot(311)
axu.plot(x,u)
axu.plot(x,u_final)
axu.set_title('u')

axp = fig.add_subplot(312)
axp.plot(x,p)
axp.plot(x,p_final)
axp.set_title('p')

axr = fig.add_subplot(313)
axr.plot(x,rho)
axr.plot(x,rho_final)
axr.set_title('rho')

plt.show()

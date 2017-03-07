import numpy as np
import matplotlib.pyplot as plt

gamma=1.4

def get_F(U,der):
	rho=U[0,:,:,:]
	u=U[1,:,:,:]/rho
	v=U[2,:,:,:]/rho
	w=U[3,:,:,:]/rho
	velocidades=[u,v,w]
	E=U[4,:,:,:]/rho
	ei=E-(u**2+v**2+w**2)/2
	p=rho*ei*(gamma-1)

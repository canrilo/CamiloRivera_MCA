import numpy as np
import matplotlib.pyplot as plt

f = open ("Precisiones_float.txt")
a =[]
b = []
x = -12
for line in f:
	#print(line.rstrip())
	a.append(float(line.rstrip()))
	b.append(float(1*pow(10,x)))
	x=x+1
f.close()

fig=plt.figure()
grafica=fig.add_subplot(111)
grafica.plot(b,a,'.-')

grafica.set_yscale('log')
grafica.set_xscale('log')
grafica.set_title("Precision float")
grafica.set_xlabel("Float number order")
grafica.set_ylabel("Precision")
#plt.show()

fig.savefig("Precisiones_float.pdf")

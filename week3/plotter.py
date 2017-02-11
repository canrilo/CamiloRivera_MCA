import numpy as np
import matplotlib.pyplot as plt

infile = open('output.txt','r')
header_size = int(infile.readline().rstrip())
m = int(infile.readline().rstrip())
h = float(infile.readline().rstrip())
numproc = int(infile.readline().rstrip())
infile.close()

data = np.genfromtxt('output.txt',skip_header=header_size)
potential = data.reshape(m, m)

xi, yi = np.linspace(0, m, m), np.linspace(0, m, m)
xi, yi = np.meshgrid(xi, yi)

fig, ax = plt.subplots()

dx, dy = np.gradient(-potential)
ax.streamplot(xi, yi, dy, dx, color=potential, cmap=plt.cm.brg, density=[2,1])

plt.xlim(0, m)
plt.ylim(m, 0)

x = np.array([1,2,3,4,5])*m/5.0
b=['1','2','3','4','5']
plt.xticks(x,b)
plt.yticks(x,b)
#ax.xaxis.set_ticks_position('top')

plt.title('Capacitor ' + str(numproc) + ' Processors')
plt.savefig('Potential_' + str(numproc) + '.png')


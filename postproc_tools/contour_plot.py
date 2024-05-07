import read_data as rd

import numpy as np 
import matplotlib.pyplot as plt
import sys 

# It expects ORDERED DATA by columns in the form:  " x y z1 z2 z3 z4... "
# x and y should be in increasing order.

file_path = sys.argv[1] #str()?
data = rd.read_data_file(file_path)

x_v = data.iloc[:,0]
y_v = data.iloc[:,1]
z_v = data.iloc[:,2]
print(x_v, y_v, z_v)

Nx = len( set(x_v) ) - 1
Ny = len( set(y_v) ) - 1 


x = np.zeros( [ Nx+1 ], dtype=float )
y = np.zeros( [ Ny+1 ], dtype=float )
z = np.zeros( [ Ny+1, Nx+1 ], dtype=float )

for j in range(Ny+1):
	for i in range(Nx+1): 
		
		k = i + j * (Nx+1)
		x[i] = x_v[k]
		y[j] = y_v[k]
		z[j,i] = z_v[k]
		print(x[i],y[j], z[j,i])

levels = 30
fig, ax = plt.subplots(1,1)
figure = ax.contour( x, y, z, levels, colors = 'k' )
figure = ax.contourf( x, y, z, levels, alpha = 0.75 )
figure.set_cmap("jet") 
fig.colorbar(figure)
ax.set(xlabel=x_v.name, ylabel=y_v.name, title=file_path)
fig.savefig('contour_plot_' + file_path +'.png')


#print(list(data))
#print( data.iloc[:, [0,2]] ) 
#print(data)



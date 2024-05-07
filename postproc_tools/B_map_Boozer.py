import MONKES_magnetic as magnetic
import numpy as np 
import matplotlib.pyplot as plt
import sys 

# Read magnetic field modes
parameters, Magnetic_field = magnetic.read_monkes_Magnetic_configuration('monkes_Magnetic_configuration.dat')
 
# Fourier poloidal (m), toroidal (n) mode numbers and Fourier (cosine) modes B_mnc
Np = parameters["Numberofperiods"]   
iota = parameters["iota"]   
B_theta = parameters["B_theta"]   
B_zeta = parameters["B_zeta"]   
psi_p = parameters["psi_p"]   

m    = Magnetic_field["m"]
n    = Magnetic_field["n"] 
B_mn = 	Magnetic_field["B_mn"]
N_modes = len(m) - 1  

# Select resolutions in theta and zeta, by default 64 x 64
if len(sys.argv) == 2: # If only one input, both resolutions are the same
	N_theta = int(sys.argv[1])
	N_zeta = N_theta
if len(sys.argv) > 2:
	N_theta = int(sys.argv[1])
	N_zeta = int(sys.argv[2])
else:
    N_theta = 64  
    N_zeta = N_theta
	
# Equispaced grids in poloidal (theta) and toroidal (zeta) Boozer angles
theta = np.linspace( 0.0, 2*np.pi, num=N_theta + 1 )
zeta  = np.linspace( 0.0, 2*np.pi/Np, num=N_zeta + 1 )

B = np.zeros( [ N_theta+1, N_zeta+1 ], dtype=float )
temp = np.zeros( [ N_modes ], dtype=float )

with open( "B_map_Boozer.dat", "w" ) as output_file: 
	print( "theta",  "zeta", "B [T]", file=output_file )
	for i in range(N_theta+1):
		for j in range(N_zeta+1): 
			B[i,j] = np.dot( B_mn , np.cos( m * theta[i] + Np * n * zeta[j] ) )
			print( theta[i],  zeta[j], B[i,j], file=output_file ) 
levels = 30
fig, ax = plt.subplots(1,1)
figure = ax.contour( zeta, theta, B, levels, colors = 'k' )
figure = ax.contourf( zeta, theta, B, levels, alpha = 0.75 )
fig.colorbar(figure)
ax.set(xlabel=r"$\zeta$", ylabel=r"$\theta$", title=r"$|B|$ in Boozer coordinates")
fig.savefig('B_map_Boozer.png')

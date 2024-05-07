import read_data as rd

import numpy as np 
import matplotlib.pyplot as plt
import sys 
import os
 
#plt.rcParams['text.usetex'] = True

file_path = "monkes_Monoenergetic_theta_zeta.dat"  
data = rd.read_data_file(file_path)

nu_v = data.iloc[:,0]
Er_v = data.iloc[:,1]
N_theta_v = data.iloc[:,2]
N_zeta_v = data.iloc[:,3]
N_xi_v = data.iloc[:,4]

N_nu =    len( set(nu_v) ) - 1
N_Er =    len( set(Er_v) ) - 1
M_theta = len( set(N_theta_v) ) - 1
M_zeta  = len( set(N_zeta_v) ) - 1
M_xi  =   len( set(N_xi_v) ) - 1
print( N_nu, N_Er, M_theta, M_zeta, M_xi )

x_v = data.iloc[:,10]
y_v = data.iloc[:,11]
D11_v = data.iloc[:,12]
D31_v = data.iloc[:,13]
D13_v = data.iloc[:,14]
D33_v = data.iloc[:,15]  


nu =      list(set(nu_v)     )#np.zeros( [ N_nu+1 ], dtype=float )
E_r =     list(set(Er_v)     )#np.zeros( [ N_Er+1 ], dtype=float )
N_theta = list(set(N_theta_v))#np.zeros( [ M_theta+1 ], dtype=int )
N_zeta =  list(set(N_zeta_v) )#np.zeros( [ M_zeta+1 ], dtype=int )
N_xi =    list(set(N_xi_v)   )#np.zeros( [ M_xi+1 ], dtype=int )
 
if not os.path.exists('D_ij_fs_maps'):
	os.mkdir("D_ij_fs_maps")
for j in range(N_Er+1):
	for i in range(N_nu+1): 
		for kk in range(M_xi+1): 
			for ii in range(M_theta+1): 
				for jj in range(M_zeta+1):
					
					print( i, j, ii, jj, kk )
					print( nu[i], E_r[j], int(N_theta[ii]), int(N_zeta[jj]), int(N_xi[kk]) )
					
					
					theta = x_v[ N_theta_v == N_theta[ii] ]
					theta = theta[ N_zeta_v == N_zeta[jj] ]
					theta = theta[ N_xi_v == N_xi[kk] ]
					theta = theta[ nu_v == nu[i] ]
					theta = theta[ Er_v == E_r[j] ] 
					
					zeta = y_v[ N_theta_v == N_theta[ii] ]
					zeta = zeta[ N_zeta_v == N_zeta[jj] ]
					zeta = zeta[ N_xi_v == N_xi[kk] ]
					zeta = zeta[ nu_v == nu[i] ]
					zeta = zeta[ Er_v == E_r[j] ] 
					
					D11 = D11_v[ N_theta_v == N_theta[ii] ]
					D11 = D11[ N_zeta_v == N_zeta[jj] ]
					D11 = D11[ N_xi_v == N_xi[kk] ]
					D11 = D11[ nu_v == nu[i] ]
					D11 = D11[ Er_v == E_r[j] ] 
					
					D31 = D31_v[ N_theta_v == N_theta[ii] ]
					D31 = D31[ N_zeta_v == N_zeta[jj] ]
					D31 = D31[ N_xi_v == N_xi[kk] ]
					D31 = D31[ nu_v == nu[i] ]
					D31 = D31[ Er_v == E_r[j] ]
					
					D13 = D13_v[ N_theta_v == N_theta[ii] ]
					D13 = D13[ N_zeta_v == N_zeta[jj] ]
					D13 = D13[ N_xi_v == N_xi[kk] ]
					D13 = D13[ nu_v == nu[i] ]
					D13 = D13[ Er_v == E_r[j] ]
					
					D33 = D33_v[ N_theta_v == N_theta[ii] ]
					D33 = D33[ N_zeta_v == N_zeta[jj] ]
					D33 = D33[ N_xi_v == N_xi[kk] ]
					D33 = D33[ nu_v == nu[i] ]
					D33 = D33[ Er_v == E_r[j] ]
					
					Nx = len( set(theta) ) - 1
					Ny = len( set(zeta) ) - 1 
					x = np.zeros( [ Nx+1 ], dtype=float )		
					y = np.zeros( [ Ny+1 ], dtype=float )		
					z11 = np.zeros( [ Nx+1, Ny+1 ], dtype=float ) 
					z31 = np.zeros( [ Nx+1, Ny+1 ], dtype=float ) 
					z13 = np.zeros( [ Nx+1, Ny+1 ], dtype=float ) 
					z33 = np.zeros( [ Nx+1, Ny+1 ], dtype=float ) 
					
					kkk_offset = D11.index[0] 
					for jjj in range(Ny+1):
						for iii in range(Nx+1):
						    kkk = iii + jjj * (Nx+1)    
						    
						    x[iii] = theta[kkk+kkk_offset]
						    y[jjj] = zeta[kkk+kkk_offset]
						    z11[iii,jjj] = D11[kkk+kkk_offset]
						    z31[iii,jjj] = D31[kkk+kkk_offset]
						    z13[iii,jjj] = D13[kkk+kkk_offset]
						    z33[iii,jjj] = D33[kkk+kkk_offset]
						     
						    
					levels = 35
					fig, ax = plt.subplots(1,1)    
					figure = ax.contour( y, x, z11, levels, colors = 'k' )    
					figure = ax.contourf( y, x, z11, levels, alpha = 0.75 )    
					figure.set_cmap("jet")     
					fig.colorbar(figure)
					ax.set(xlabel=r"$\zeta$", ylabel=r'$\theta$', title=r"$d_{11}(\theta,\zeta)$"+r" for $(\nu,E_r)=($"+str(nu[i])+","+str(E_r[j])+r"$)$")
					fig.savefig('D_ij_fs_maps/D11_nu_' + str("{:.8f}".format(nu[i]) ) +"_Er_" + str( "{:.8f}".format(E_r[j]))  +"_Ntheta_"+ str(int(Nx)) + "_Nzeta_"+str(int(Ny)) + "_Nxi_"+ str(int(N_xi[kk])) + ".png" )
					plt.close()
						    
					levels = 35
					fig, ax = plt.subplots(1,1)    
					figure = ax.contour( y, x, z31, levels, colors = 'k' )    
					figure = ax.contourf( y, x, z31, levels, alpha = 0.75 )    
					figure.set_cmap("jet")     
					fig.colorbar(figure)
					ax.set(xlabel=r"$\zeta$", ylabel=r'$\theta$', title=r"$d_{31}(\theta,\zeta)$"+r" for $(\nu,E_r)=($"+str(nu[i])+","+str(E_r[j])+r"$)$")
					fig.savefig('D_ij_fs_maps/D31_nu_' + str("{:.8f}".format(nu[i])) +"_Er_" + str("{:.8f}".format(E_r[j])) +"_Ntheta_"+ str(int(Nx)) + "_Nzeta_"+str(int(Ny)) + "_Nxi_"+ str(int(N_xi[kk])) + ".png" )
					plt.close()
						    
					levels = 35
					fig, ax = plt.subplots(1,1)    
					figure = ax.contour( y, x, z13, levels, colors = 'k' )    
					figure = ax.contourf( y, x, z13, levels, alpha = 0.75 )    
					figure.set_cmap("jet")     
					fig.colorbar(figure)
					ax.set(xlabel=r"$\zeta$", ylabel=r'$\theta$', title=r"$d_{13}(\theta,\zeta)$"+r" for $(\nu,E_r)=($"+str(nu[i])+","+str(E_r[j])+r"$)$")
					fig.savefig('D_ij_fs_maps/D13_nu_' + str("{:.8f}".format(nu[i])) +"_Er_" + str("{:.8f}".format(E_r[j])) +"_Ntheta_"+ str(int(Nx)) + "_Nzeta_"+str(int(Ny)) + "_Nxi_"+ str(int(N_xi[kk])) + ".png" )
					plt.close()
						    
					levels = 35
					fig, ax = plt.subplots(1,1)    
					figure = ax.contour( y, x, z33, levels, colors = 'k' )    
					figure = ax.contourf( y, x, z33, levels, alpha = 0.75 )    
					figure.set_cmap("jet")     
					fig.colorbar(figure)
					ax.set(xlabel=r"$\zeta$", ylabel=r'$\theta$', title=r"$d_{33}(\theta,\zeta)$"+r" for $(\nu,E_r)=($"+str(nu[i])+","+str(E_r[j])+r"$)$")
					fig.savefig('D_ij_fs_maps/D33_nu_' + str("{:.8f}".format(nu[i])) +"_Er_" + str("{:.8f}".format(E_r[j])) +"_Ntheta_"+ str(int(Nx)) + "_Nzeta_"+str(int(Ny)) + "_Nxi_"+ str(int(N_xi[kk])) + ".png" )
					plt.close()

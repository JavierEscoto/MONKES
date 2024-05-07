import sys
import MONKES_input as input

# Read resolutions and collisionality 
if len(sys.argv) >=6:
	N_theta = int(sys.argv[1])
	N_zeta = int(sys.argv[2]) 
	N_xi = int(sys.argv[3]) 
	nu = float(sys.argv[4]) 
	E_r = float(sys.argv[5])

if len(sys.argv) == 8:
	N_xi_DF = int(sys.argv[6]) 
	N_lambda = int(sys.argv[7]) 
    # Create monkes input
	input.create_monkes_input_parameters_DF(N_theta, N_zeta, N_xi, nu, E_r, N_xi_DF, N_lambda)
elif len(sys.argv) == 6:
	input.create_monkes_input_parameters(N_theta, N_zeta, N_xi, nu, E_r)
else:
	input.create_monkes_input_parameters(40, 40, 140, 1e-4, 1e-4) 


def create_monkes_input_parameters(N_theta, N_zeta, N_xi, nu, E_r):	   	
    
    with open( "monkes_input.parameters", "w" ) as output_file: 
    	print( "&parameters", file=output_file )
    	print( "N_theta =",  N_theta, file=output_file )
    	print( "N_zeta =",  N_zeta, file=output_file )
    	print( "N_xi =",  N_xi, file=output_file )
    	print( "nu =",  "{:e}".format(nu), file=output_file )
    	print( "E_r =",  "{:e}".format(E_r), file=output_file )
    	print( "/", file=output_file ) 
      




def create_monkes_input_parameters_DF(N_theta, N_zeta, N_xi, nu, E_r, N_xi_DF, N_lambda):	   	
    
    with open( "monkes_input.parameters", "w" ) as output_file: 
    	print( "&parameters", file=output_file )
    	print( "N_theta =",  N_theta, file=output_file )
    	print( "N_zeta =",  N_zeta, file=output_file )
    	print( "N_xi =",  N_xi, file=output_file )
    	print( "nu =",  "{:e}".format(nu), file=output_file )
    	print( "E_r =",  "{:e}".format(E_r), file=output_file )
    	print( "/", file=output_file )     	
    	print( "&parameters_DF", file=output_file ) 
    	print( "N_xi_DF =",  N_xi_DF, file=output_file )
    	print( "N_lambda =",  N_lambda, file=output_file ) 
    	print( "/", file=output_file ) 
      
      

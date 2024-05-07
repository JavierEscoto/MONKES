import pandas as pd 

def read_ddkes2data(N_theta, N_zeta): 
    file=open("ddkes2.data")
    for x in file.readlines():
	    if 'borbi' in x:
	       borbi_line = x 
	       print(borbi_line)
  
    file.close()

##############################################
def read_monkes_Magnetic_configuration(file_path):	
	  
    # Read parameters
    parameters={} #Initialize dictionary
    with open( file_path ) as myfile: 
        for line in myfile.readlines(): 
            if "Fourier" in line: 
                break
            if "***" in line: 
                continue 
            temp = line.replace(' ', '')
            temp = temp.split('=')
            if len(temp) >= 2:
                parameters[temp[0]] = float( temp[1] ) 
     
    Fourier_modes_df=pd.read_table(file_path, skiprows=16, sep="\s+", engine='python') 
    
    return parameters, Fourier_modes_df 
    
	
	

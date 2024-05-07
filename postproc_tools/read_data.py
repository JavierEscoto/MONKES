import pandas as pd 

# Reads a file which has a header: x1 x2... xN
def read_data_file(file_path):	  
     
    data=pd.read_table(file_path, sep="\s+", engine='python') 
    
    return data 

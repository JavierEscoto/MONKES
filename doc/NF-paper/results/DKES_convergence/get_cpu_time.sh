rm cpu_time.plt
echo " N_xi    t_cpu [s] " > cpu_time.plt
for folder in */ ; do
        echo $folder
        cd $folder
        lanxi=$(echo "$folder" | cut -d'/' -f1)         
        
        nxi=$(echo "$lanxi" | cut -d'a' -f2) 
        echo $nxi
        
        line=$(grep 'time used' dkesout.stellopt)
        
        timesec=$(echo "$line" | cut -d'=' -f2)
        
        time=$(echo "$timesec" | cut -d's' -f1)
        echo $time
        
        cd ../ 
        echo  $nxi " " $time >> cpu_time.plt
done

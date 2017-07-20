


echo "Checking results..."

for i in q2D 2D 3D
do
    if ! paste CPU/$i.gdr GPU/$i.gdr | awk '{if($2-$5 > 1e-2){exit 1}}'
    then	
	echo "ERROR!: CPU and GPU versions in $i do not match!!"
	echo "Not Passed!!"
	exit 1
    fi
done    
    echo "CPU and GPU versions match!!"
    echo "Passed!"

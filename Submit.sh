#! /bin/bash
echo 'Starting Job' 

for i in $(seq 1 $2); 
do 
    name=$1\_$i
    echo $name; 
    eval ./submit $name --medium
done

#! /bin/bash

if [ "${1}"  == "--help" ]; then
    echo "Options for running this configuration file:"
    echo "--help                    Prints this message"
    echo "example:  ./Submit.sh  <Jobs output file prefix> <Numver of Jobs>"



else



    echo 'Starting Job' 

    for i in $(seq 1 $2); 
    do 
	name=$1\_$i
	echo $name; 
	eval ./submit $name --medium
    done
fi
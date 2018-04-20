#! /bin/bash
I_size=32
Mu=$1
T=$2

if (( $# != 2 ))
then
  echo "Need Mu and T set"
  exit 1
fi

for V in 0 8 16 24 32 40
do
	./ProteinEvolution -E -N 3 -P 10 -K 10 -B 20 -R 0 -F 1 -M $Mu -T $T -X .51 -I .25 -D 8 -V $V 
	python ~/Documents/PolyominoDev/scripts/interface_analysis.py $I_size $Mu $T 8
	XZ_OPT=-9 tar -Jcvf /rsratch/asl47/Bulk_Run/DataV$V.tar.xz /rscratch/asl47/Bulk_Run/Interfaces/*txt --remove-files
done







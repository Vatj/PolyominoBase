#! /bin/bash


I_size=64
Mu=$1
T=$2
STEP_SIZE=4
RUNS=40

let "END = $RUNS - $STEP_SIZE"

if (( $# < 2 ))
then
  echo "Need Mu and T set"
  exit 1
fi

argc=$#
argv=($@)

    


mkdir -p "/scratch/asl47/Data_Runs/Bulk_Data"
cd /scratch/asl47/Data_Runs/Bulk_Data

for (( j=1; j<argc; j++ )); do
	T=${argv[j]}
	echo "on temperature $T"
	for (( V=0; V<=$END; V+=$STEP_SIZE )); do

		echo "starting evolution set $V"
		nice -n 15 ~/Documents/PolyominoDev/bin/ProteinEvolution -E -N 3 -P 500 -K 1500 -B 20 -R 1 -F 1 -M $Mu -T $T -X .51 -I .25 -D $STEP_SIZE -V $V
	
		nice -n 15 python ~/Documents/PolyominoDev/scripts/interface_analysis.py $I_size $Mu $T $STEP_SIZE $V
	XZ_OPT=-9 tar -Jcf DataI${I_size}Mu${Mu}T${T}V${V}.tar.xz *txt --remove-files
	done
done







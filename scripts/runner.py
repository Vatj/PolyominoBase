import subprocess
from numpy import linspace,logspace,log10,append
import sys




## RUN PARAMETERS ##

executable='ProteinEvolution'
ex_flags='-E -N 4 -P 500 -K 10000 -B 20 -X 0.2 -S 1 -F 1 -M {0} -T {1}  -D {2} -V {3}'
Runs=12
command_sequence='./{0} {1}'.format(executable,ex_flags)
print command_sequence

#subprocess.call(command_sequence,shell=True)

def main():
    if len(sys.argv)<2:
        print "V value required"
        exit()
    print "Running executable ",executable
    print "with params ",ex_flags
    for M in logspace(log10(0.25),log10(4),15):
        for T in append(logspace(-2,log10(0.05),8),[0.75,1]):
            command_sequence='./../bin/{0} {1}'.format(executable,ex_flags.format(M,T,Runs,Runs*int(sys.argv[1])))
            print command_sequence
    #subprocess.call(command_sequence,shell=True)
    
    


if __name__ == "__main__":
    main()

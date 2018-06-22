import glob

def Compact_File():
    with open('/rscratch/asl47/Reg_Lengths_Full_T5.txt','w') as f:
        for i,fileN in enumerate(glob.glob('/scratch/asl47/*.txt')):
            if i%5000==0:
                print("on ",i)
            lengths=[line.rstrip('\n').split() for line in open(fileN)]
            for val in lengths[0]:
                f.write(val+" ")
            f.write('\n')

def main():
    Compact_File()

if __name__ == "__main__":
    main()

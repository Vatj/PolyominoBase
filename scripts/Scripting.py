import itertools
def ya():
    num=0
    num2=0
    for comb in itertools.combinations_with_replacement([0,1,2,3,4,5,6,7],4):
        num+=1
        #print comb
    print num
    for comb in itertools.product([0,1,2],repeat=3):
        print comb
        num2+=1
    print num2

    
    for i in xrange(3):
        for j in xrange(i,3):
            print "hi"


        
def w():
    print "w"


def pointless():
    import matplotlib.pyplot as plt
    from math import factorial
    upperL=31
    tiles=8
    plt.plot(xrange(1,upperL),[factorial(n+tiles-1)/(factorial(n-1)*factorial(tiles)*1.) for N in xrange(1,upperL)])
    plt.plot(xrange(1,upperL),[n**tiles for n in xrange(1,upperL)])
    plt.plot(xrange(1,upperL),[(n**tiles)/(factorial(n+tiles-1)/(factorial(n-1)*factorial(tiles)*1.)) for n in xrange(1,upperL)])
    plt.yscale('log')
    plt.show()

def z():
    from collections import deque
    x=set()
    N=0
    for i in xrange(0,2):
        jrange=(1-i%2)*(i+1) + (i%2)*(i+2)
        #print jrange
        for j in xrange(0,jrange+1):
            krange=(1-max(i,j)%2)*(max(i,j)+1) + (max(i,j)%2)*(max(i,j)+2)
            for k in xrange(0,krange+1):
                lrange=(1-max(k,max(i,j))%2)*(max(k,max(i,j))+1) + (max(k,max(i,j))%2)*(max(k,max(i,j))+2)
                for l in xrange(0,lrange+1):
                    qrange=(1-max([i,j,k,l])%2)*(max([i,j,k,l])+1) + (max([i,j,k,l])%2)*(max([i,j,k,l])+2)
                    for q in xrange(0,qrange+1):
                        wrange=(1-max([i,j,k,l,q])%2)*(max([i,j,k,l,q])+1) + (max([i,j,k,l,q])%2)*(max([i,j,k,l,q])+2)
                        for w in xrange(0,wrange+1):
                            erange=(1-max([i,j,k,l,q,w])%2)*(max([i,j,k,l,q,w])+1) + (max([i,j,k,l,q,w])%2)*(max([i,j,k,l,q,w])+2)
                            for e in xrange(0,erange+1):
                                rrange=(1-max([i,j,k,l,q,w,e])%2)*(max([i,j,k,l,q,w,e])+1) + (max([i,j,k,l,q,w,e])%2)*(max([i,j,k,l,q,w,e])+2)
                                #print max([i,j,k,l,q,w,e])
                                for r in xrange(0,rrange+1):
                                    d=deque([i,j,k,l,q,w,e,r])
                                    rot=0
                                    while rot>-8:
                                        if tuple(d) in x:
                                            break
                                        else:
                                            d.rotate(-1)
                                            rot-=1
                                            
                                    else:
                                        x.add((i,j,k,l,q,w,e,r))
                                        #print i,j,k,l
                                    N+=1
    print N
    print len(x)
    return x

from collections import deque
def zzz(depth,seq,inrange,maxdepth=12,maxcolor=13,N=0,y=set()):
    if depth>=maxdepth:
        
        #d=deque(seq)
        #rot=0

        #while rot>-maxdepth:
        #    if tuple(d) in y:
        #        break
        #    else:
        #        d.rotate(-1)
        #        rot-=1
        #                                    
        #else:
        #    y.add(tuple(seq))
        #    #print seq
        N+=1

        #print seq
        #N+=1
        return N

    for i in xrange(0,inrange+1):
        #print "on i:",i,"  max depth is ",depth
        if depth>=1:
            seq.append(i)
            jrange=min(maxcolor,(1-max(seq)%2)*(max(seq)+1) + (max(seq)%2)*(max(seq)+2))
 
            del seq[len(seq)-1]
        else:
            jrange=min(maxcolor,(1-i%2)*(i+1) + (i%2)*(i+2))

        le=len(seq)
        seq.append(i)
        N=zzz(depth+1,seq,jrange,maxdepth,maxcolor,N,y)
        del seq[le]

    #if depth==0:
    #    print x
    return N

    
def x():
    N=0
    for i in xrange(1,-1,-1):
        jrange=(1-i%2)*(i+1) + (i%2)*(i+2)
        for j in xrange(jrange,i-1,-1):
            krange=(1-j%2)*(j+1) + (j%2)*(j+2)
            
            for k in xrange(krange,j-1,-1):
                lrange=(1-k%2)*(k+1) + (k%2)*(k+2)
                for l in xrange(lrange,k-1,-1):
                    print i,j,k,l
                    N+=1

                    
    print N
    
def xxx(depth,seq,inrange,minV,maxdepth=12,maxcolor=13,N=0):
    
    if depth>=maxdepth:

        #print seq
        N+=1
        return N

    for i in xrange(inrange,minV-1,-1):
        jrange=min(maxcolor,(1-i%2)*(i+1) + (i%2)*(i+2))
        le=len(seq)
        seq.append(i)
        N=xxx(depth+1,seq,jrange,i,maxdepth,maxcolor,N)
        del seq[le]
        
    return N

        #for j in xrange(jrange,i-1,-1):
        #    krange=(1-j%2)*(j+1) + (j%2)*(j+2)
        #    for k in xrange(krange,j-1,-1):
        #        lrange=(1-k%2)*(k+1) + (k%2)*(k+2)
        #        for l in xrange(lrange,k-1,-1):
        #            print i,j,k,l
        #            N+=1

                    
#print N


import numpy as np
def mmm():

    #plt.contour(xrange(20),xrange(50),np.array([xxx(0,[],1,0,i,j,0) for i in xrange(40) for j in xrange(50)]).reshape(40,-1),20)
    for i in xrange(4,20,4):
        plt.plot(xrange(i*3),[xxx(0,[],1,0,i,j,0) for j in xrange(i*3)],label=i)
    plt.legend()
    plt.yscale('log')
    plt.show()



def ppp():
    #x=np.linspace(0.01,.1,10);
    x=np.logspace(-2,-1,15);
    print x
    #y=np.linspace(.1,1,10);
    y=np.logspace(-1,0,20);
    z=np.logspace(0,1,15);
    print y
    print z
    q=np.array([]);
    q=np.append(q,x);
    q=np.append(q,y);
    q=np.append(q,z);
    w=np.logspace(-2,1,30);
    plt.plot(q,[1]*q.shape[0],'ro')
    plt.plot(w,[1]*w.shape[0],'gx')
    plt.xscale('log')
    plt.show()



def fukm():
    lineF=[line.rstrip('\n') for line in open("RunConfigs/TileTest_2_16.txt")]
    lineI=[[int(i) for i in line.split(" ")] for line in lineF]
    print "Starting off with ",len(lineI)," possibilities"
    lineT=[[line[:4],line[4:]] for line in lineI]
    #return lineT[:10]
    uniques=[]
    dups=0
    for i,line in enumerate(lineT):
        #print "HERE",line
        if(i%1000==0):
            print "On ",i
            #print uniques
            
        for uniq in uniques:
            tuple1=line[0]
            tuple2=line[1]
            if tuple1 in uniq and tuple2 in uniq and tuple1!=tuple2:
                print tuple1,"and",tuple2," are in ",line
                #print line
                dups+=1
                break
        else:
            uniques.append([line[0],line[1]])
    #print x
    print "Reduced to ",len(uniques)," possibilities"
    print "Dups found is ",dups
    uniqs=[y[0]+y[1] for y in uniques]
    return uniqs
    with open("Output/Unique_2_8_Space.txt", "w") as text_file:
        for uniq in uniqs:
            for charac in uniq:
                text_file.write(str(charac)+" ")
                #text_file.write(" ")
            text_file.write('\n')
def tm(x):
    for ln in x:
        print ln
        print ln[:4]
        print ln[4:]
        if [0,0,0,1] in ln and [0,0,0,0] in ln:
            print "hi"


########################
#SEARCH SPACE SIZE CODE#
########################
def divisor_of(n):
    for i in xrange(1,n+1):
        if n%i==0:
            yield i
from fractions import gcd
def Euler_Totient(d):
    total=0
    for k in xrange(1,d+1):
        if gcd(d,k)==1:
            total+=1
    return total

def Necklace_Number(k,n):
    N_k=0.
    for d in divisor_of(n):
       N_k+=Euler_Totient(d)*k**(n*1./d) 
    N_k/=n
    return int(N_k)

def Effective_Numbering_Number(depth,seq,inrange,maxdepth=12,maxcolor=13,N=0):
    if depth>=maxdepth:
        N+=1
        return N

    for i in xrange(0,inrange+1):
        if depth>=1:
            seq.append(i)
            jrange=min(maxcolor,(1-max(seq)%2)*(max(seq)+1) + (max(seq)%2)*(max(seq)+2))
            del seq[len(seq)-1]
        else:
            jrange=min(maxcolor,(1-i%2)*(i+1) + (i%2)*(i+2))
        le=len(seq)
        seq.append(i)
        N=Effective_Numbering_Number(depth+1,seq,jrange,maxdepth,maxcolor,N)
        del seq[le]
    return N
def attempt(depth,max_D,C_max,C_MAX,inRange,N):
    if depth>=max_D:
        return N+1
    for i in xrange(inRange+1):
        newRange=min(C_MAX,max(C_max,(1-i%2)*(i+1) + (i%2)*(i+2)))

        N=attempt(depth+1,max_D,max(C_max,(1-i%2)*(i+1) + (i%2)*(i+2)),C_MAX,newRange,N)

    return N
def N_Necklace(colours,tiles):
    return Necklace_Number(colours,4)**tiles
def N_Brute_Force(colours,tiles):
    NUM_FACES=4
    return colours**(tiles*NUM_FACES)
def N_Effective_Numbering(colours,tiles):
    return 0


import matplotlib.pyplot as plt
def plot_Search_Spaces():
    f, (ax1,ax2) = plt.subplots(1, 2)
    Tile_Space=xrange(1,6)

    N_1_8=57
    N_1_16=57
    N_2_16=33955
    N_2_8=28060
    N_3_8=15984422
    N_3_16=57986025
    N_4_8=7076011589
    EN_1_8=208
    EN_2_8=361468
    EN_3_8=1264020328
    EN_1_16=208
    EN_2_16=411200
    EN_3_16=3055284360
    ax1.plot(Tile_Space,[N_Brute_Force(8,N_t) for N_t in Tile_Space],label='Enumerative',lw=2)
    ax1.plot(Tile_Space,[N_Necklace(8,N_t) for N_t in Tile_Space],label='Necklace',lw=2)
    #ax1.plot([1,2,3],[EN_1_8,EN_2_8,EN_3_8],label='Effective Numbering',lw=2)
    ax1.plot([1,2,3],[N_1_8,N_2_8,N_3_8],label='Neck+EN+Sym',lw=2)

    ax1.legend(loc=4)
    ax1.set_yscale('log')
    ax1.set_xticks([1,2,3,4,5,6])
    ax1.set_title('8 Colours (2N 6I)')
    ax1.set_xlabel('Number of Tiles')
    ax1.set_ylabel('Search Space Size')
    ax2.plot(Tile_Space,[N_Brute_Force(16,N_t) for N_t in Tile_Space],lw=2)
    ax2.plot(Tile_Space,[N_Necklace(16,N_t) for N_t in Tile_Space],lw=2)
    
    #ax2.plot([1,2,3],[EN_1_16,EN_2_16,EN_3_16],lw=2)
    ax2.plot([1,2,3],[N_1_16,N_2_16,N_3_16],lw=2)
    ax2.set_yscale('log')
    ax2.set_xticks([1,2,3,4,5,6])
    ax2.set_title('16 Colours (2N 14I)')
    ax2.set_xlabel('Number of Tiles')
    ax2.set_ylabel('Search Space Size')

    plt.show()
import seaborn as sns
def plot_Search_Spaces2():
    Tile_Space=xrange(1,6)
    sns.set_context("paper",font_scale=2.2)
    sns.set_style("white",rc={"xtick.major.size": 8, "ytick.major.size": 8,"xtick.minor.size":5, "ytick.minor.size": 5,"axes.linewidth": 2,"axes.edgecolor":"darkgray","font.size":8,"axes.titlesize":8,"axes.labelsize":5})
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    N_1_8=57
    N_1_16=57
    N_2_16=33955
    N_2_8=28060
    N_3_8=15984422
    N_3_16=57986025
    N_4_8=7076011589
    N_4_16=209294594809

    EN_1_8=208
    EN_2_8=361468
    EN_3_8=1264020328
    EN_1_16=208
    EN_2_16=411200
    EN_3_16=3055284360

    plt.figure()
    plt.plot(Tile_Space,[N_Brute_Force(8,N_t) for N_t in Tile_Space],lw=2,ls='--',color='firebrick',alpha=0.8)
    plt.plot(Tile_Space,[N_Brute_Force(8,N_t) for N_t in Tile_Space],ls='',markersize=7,marker='o',markeredgecolor='firebrick',markeredgewidth=1.25,markerfacecolor='white')
    plt.plot([0],[0],label='Enumerative',lw=2,ls='--',markersize=7,color='firebrick',marker='o',markeredgecolor='firebrick',markeredgewidth=1.25,markerfacecolor='white')


    
    plt.plot(Tile_Space,[N_Necklace(8,N_t) for N_t in Tile_Space],lw=2,ls='--',color='dodgerblue',alpha=0.8)
    plt.plot(Tile_Space,[N_Necklace(8,N_t) for N_t in Tile_Space],ls='',markersize=7,marker='s',markeredgecolor='dodgerblue',markeredgewidth=1.25,markerfacecolor='white',zorder=10)
    plt.plot(Tile_Space,[N_Necklace(8,N_t) for N_t in Tile_Space],label='Necklace',lw=2,ls='--',color='dodgerblue',markersize=7,marker='s',markeredgecolor='dodgerblue',markeredgewidth=1.25,markerfacecolor='white')

    
    #plt.plot([1,2,3],[EN_1_8,EN_2_8,EN_3_8],label='Effective Numbering',lw=2,color='g')
    plt.plot([1,2,3,4],[N_1_8,N_2_8,N_3_8,N_4_8],lw=2,ls='--',color='forestgreen')
    plt.plot([1,2,3,4],[N_1_8,N_2_8,N_3_8,N_4_8],marker='D',markersize=7,ls='',markeredgewidth=1.25,markerfacecolor='white',zorder=10,markeredgecolor='forestgreen')
    plt.plot([0],[0],marker='D',markersize=7,ls='--',label='Combined Approach',markeredgewidth=1.25,markerfacecolor='white',zorder=10,markeredgecolor='forestgreen',color='forestgreen',lw=2)
    plt.plot([1,2,3],[11,1145,3750206],label='Topologies',ls='',marker='o',color='goldenrod')
    plt.plot([3,3],[2000000,3750206],ls='-',color='goldenrod')
    plt.plot([2.975,3.025],[2000000,2000000],ls='-',color='goldenrod')
    #plt.plot([1,2,3,4],[n1*1./n2 for n1,n2 in zip([N_Brute_Force(8,N_t) for N_t in Tile_Space],[N_1_8,N_2_8,N_3_8,N_4_8])],lw=2,color='k',marker='o',ls='',alpha=0.8,label='Reduction Factor')
    #plt.plot([1,2,3,4],[n1*1./n2 for n1,n2 in zip([N_Brute_Force(16,N_t) for N_t in Tile_Space],[N_1_16,N_2_16,N_3_16,N_4_16])],lw=2,color='k',marker='o',ls='',alpha=0.8)
    #plt.legend(loc='upper left')
    plt.yscale('log')
    plt.xticks([1,2,3,4,5,6])
    #plt.title('Solid-8 Colours, Dashed-16 Colours ')
    plt.xlabel(r'$N_{T}$')
    plt.ylabel(r'$N_{\mathrm{g}}$')

    plt.text(4.25,7*10**15,"Enumerative",ha='right',fontsize=16)
    plt.text(4.25,3.5*10**12,"Necklace",ha='left',va='top',fontsize=16)
    plt.text(3.95,8.9*10**8,"Combined",ha='left',va='top',fontsize=16)
    plt.text(2.15,850,"Topologies",ha='left',va='top',fontsize=16)
    #plt.plot(Tile_Space,[N_Brute_Force(16,N_t) for N_t in Tile_Space],lw=2,color='cornflowerblue',ls='-')
    #plt.plot(Tile_Space,[N_Necklace(16,N_t) for N_t in Tile_Space],lw=2,color='cornflowerblue',ls='--')
    
    #plt.plot([1,2,3],[EN_1_16,EN_2_16,EN_3_16],lw=2,color='g',ls='--')
    #plt.plot([1,2,3,4],[N_1_16,N_2_16,N_3_16,N_4_16],lw=2,color='cornflowerblue',ls=':')
    #sns.despine()
    plt.show(block=False)





def Check_Necklace(N,j):
    if N%j==0:
        return True
def List2Int(nums):
    return int(''.join(map(str, nums)))

def Is_Necklaceble(baseTile,newTile):
    maxB=max(baseTile)
    Max_Possible=(1-maxB%2)*(maxB+1) + (maxB%2)*(maxB+2)
    for face in newTile:
        if face > Max_Possible:
            return False
        elif ((1-face%2)*(face+1) + (face%2)*(face+2))> Max_Possible:
            Max_Possible=((1-face%2)*(face+1) + (face%2)*(face+2))
    return True
            
            
    
def Check_Writing_Unique(currentTiles,checkingTile):
    if len(currentTiles)==4:
        if Is_Necklaceble([0,0,0,0],checkingTile):
            if List2Int(currentTiles[-4:]) > List2Int(checkingTile):
                return False
            else:
                return True
        else:
            return True
    else:
        possible=False
        earliest=len(currentTiles)
        for index,prevTile in enumerate([currentTiles[i*4:i*4+4] for i in xrange(len(currentTiles)/4 - 1)]):
            if Is_Necklaceble(prevTile,checkingTile):
                earliest=min(earliest,index)
                possible=True
        else:
            if not possible:
                return True
            for prevTile in [currentTiles[jj*4:jj*4+4] for jj in xrange(earliest+1,len(currentTiles)/4)]:
                if List2Int(prevTile) > List2Int(checkingTile):
                    return False
                
            else:
                return True



Numm=0
Master_L=[]
def FKM_CAT(N,K,K_prime,n,j,charray,shrink,shrinkMin,existingLine,fname):
    if n > N:
        if Check_Necklace(N,j):
            writeIt=True
            
            if type(existingLine)==list:
                if  Check_Writing_Unique(existingLine,charray[1:]):
                    writeIt=True
                    #existingStr=str()
                    #for exi in existingLine:
                    #    existingStr+=str(exi)+" "
                    #fname.write(existingStr)
                else:
                    #print charray[1:],existingLine[-8],existingLine[-6],existingLine[-4],existingLine[-2]
                    writeIt=False
            if writeIt:
                global Numm
                Numm+=1
                #print charray[1:]
                Master_L.append(charray[1:])
                #for byteCount in xrange(2):
                #    firstByteHalf=charray[1+byteCount*2]
                #    secondByteHalf=charray[2+byteCount*2]
                #    fullByte=(firstByteHalf<<4)+secondByteHalf
                #    fname.write(struct.pack('B',fullByte))
                #fname.write("{} {} {} {}\n".format(*charray[1:]))
            #print charray[1:]
    else:
        charray[n]=charray[n-j]
        maxE=max(charray[1:1+n])
        jump= (1-maxE%2)*(maxE+1) + (maxE%2)*(maxE+2)
        kp=-1
        if jump<2:
            kp=2
        else:
            kp=jump+1
        K_eff=min(kp,K)
        if not shrink:
            K_eff=min(max(shrinkMin,kp),K)

        
        FKM_CAT(N,K,K_eff,n+1,j,charray,shrink,shrinkMin,existingLine,fname)
        for i in xrange(charray[n-j]+1,K_eff):
            charray[n]=i
            FKM_CAT(N,K,K_eff,n+1,n,charray,shrink,shrinkMin,existingLine,fname)
def pr():
    print Numm

import struct
def Make_Primary_Necklaces(N,K):
    global Numm
    Numm=0
    #f=open("/scratch/asl47/Search_Spaces/{}.txt".format(fname),'w')
    squeezeFile=open("/scratch/asl47/Search_Spaces/SearchSpace_1_{}.txt".format(K),'wb')
    FKM_CAT(N,K,2,1,1,[0,0,0,0,0],True,0,-1,squeezeFile)
    print "did this many ", Numm
    #f.close()
    squeezeFile.close()
    
def Make_Further_Necklaces(N,K,fame,T):
    #lines = [[int(i) for i in line.rstrip('\n').split(" ")] for line in open("/scratch/asl47/Search_Spaces/{}.txt".format(fame))]
    #print "have ",len(lines)," to process"
    #f=open("RunConfigs/{}.txt".format("TileTest_2_16"),'w')
    f=open("/scratch/asl47/Search_Spaces/{}_{}_{}.txt".format("SearchSpace",T,K),'wb')
    fN=open("/scratch/asl47/Search_Spaces/{}_{}_{}.txt".format("Numbering",T,K),'w')
    #tileGenerator=Read_Packed_Binary(fame)
    tileGenerator=Combine_Three_Tiles(K)
    #tileGenerator=Combine_Two_Tiles(K)
    for line in tileGenerator:
        global Numm
        Numm=0
        maxE=max(line)
        shrunk= (1-maxE%2)*(maxE+1) + (maxE%2)*(maxE+2)        
        FKM_CAT(N,K,2,1,1,[0,0,0,0,0],False,shrunk+1,line,f)
        fN.write(str(Numm)+'\n')
    f.close()
    fN.close()
    print "!!!wooo finished!!!"
    print "did"
    pr()


def Read_Packed_Binary(fileName):
    with open('/scratch/asl47/Search_Spaces/{}.txt'.format(fileName), 'rb') as f:
        byte = f.read(1)
        tile=[]
        while byte != "":
            unpackedB= struct.unpack('B',byte)[0]
            unpB1=unpackedB>>4
            unpB2=unpackedB & (0b00001111)
            tile.append(unpB1)
            tile.append(unpB2)
            if len(tile)==4:
                yield tile
                tile=[]   
            byte = f.read(1)
            
def Combine_Tiles(K):
    tile1Gen=Read_Packed_Binary('SearchSpace_1_{}'.format(K))
    for tile1 in tile1Gen:
        yield tile1

def Combine_Two_Tiles(K):
    with open('/scratch/asl47/Search_Spaces/Numbering_2_{}.txt'.format(K),'r') as numbering:
        numbers=[int(line.rstrip('\n')) for line in numbering]
        #with open('/scratch/asl47/Search_Spaces/SearchSpace_1_8_Test.txt','rb') as Tile1:
            #with open('/scratch/asl47/Search_Spaces/SearchSpace_Test_8.txt','rb') as Tile2:
        tile1Gen=Read_Packed_Binary('SearchSpace_1_{}'.format(K))
        tile2Gen=Read_Packed_Binary('SearchSpace_2_{}'.format(K))
        with open('/scratch/asl47/Search_Spaces/Test13.txt','w') as f:
            for tile1 in tile1Gen:
                for t2 in xrange(numbers.pop(0)):
                    yield tile1+next(tile2Gen)
                    #fullT=tile1+next(tile2Gen)
                    #x=''
                    #for t in fullT:
                    #    x+=str(t)+" "
                    #f.write(x[:-1]+'\n')

def Combine_Three_Tiles(K):
    with open('/scratch/asl47/Search_Spaces/Numbering_3_{}.txt'.format(K),'r') as numbering2:
        numbers2=[int(line.rstrip('\n')) for line in numbering2]

        
        tile1Gen=Read_Packed_Binary('SearchSpace_1_{}'.format(K))
        tile2Gen=Read_Packed_Binary('SearchSpace_2_{}'.format(K)) 
        tile3Gen=Read_Packed_Binary('SearchSpace_3_{}'.format(K))
            #with open('/scratch/asl47/Search_Spaces/Test14.txt','w') as f:
        for tile1 in Combine_Two_Tiles(16):
            for t3 in xrange(numbers2.pop(0)):
                
                yield tile1+next(tile3Gen)
                    #fullT=tile1+tile2+next(tile3Gen)
                    #x=''
                    #for t in fullT:
                #    x+=str(t)+" "
                            #f.write(x[:-1]+'\n')
def Do_It():
    full_list=[]
    j=0
    with open('/scratch/asl47/Search_Spaces/Topology_3_14_RAW.txt','w') as f:
        for i in Combine_Three_Tiles(16):
            if 14 in i or 15 in i:
                continue
            else:
                x=''
                for t in i:
                    x+=str(t)+' '
                f.write(x[:-1]+'\n')


        
        
    return full_list
def main():
    print "Starting on Run"
    #print "doing 4_8"
    #Make_Further_Necklaces(4,8,"SearchSpace_3_8",4)
    #print Effective_Numbering_Number(0,[],2,12,7,0)
    #print "doing 4_16"
    #print Effective_Numbering_Number(0,[],2,12,15,0)
    print "Running on Local Scratch for 4 tiles 16 colours super mode"
    Make_Further_Necklaces(4,16,"SearchSpace993_998",4)
    #FKM_CAT(4,3,2,5,1,[0,0,0,0,0],False,1,3,2)
    #Make_Primary_Necklaces(4,16)
    
if __name__ == "__main__":
    #print "hi"
    main()
def ff(s,a):
    return 1./(a*(s-12)**2 +1)


def Find_Degeneracy(tile_Set,colours):
    degeneracy_Factor=1
    unused_Colours=colours-2
    for index,tile in enumerate(tile_Set):
        if tile==0 or tile==(colours-1):
            degeneracy_Factor*=2
        else:
            conjugate_Tile=(tile % 2)*(tile+1)+(1-tile % 2)*(tile-1)
            if conjugate_Tile in tile_Set[:index] or tile in tile_Set[:index]:
                continue
            else:
                degeneracy_Factor*=unused_Colours
                unused_Colours-=2
    return degeneracy_Factor*Symmetry_Factor(tile_Set)

from collections import Counter
def Symmetry_Factor2(tile_Set):
    if len(set(tile_Set))==1:
        return 1
    elif len(set(tile_Set))==2:
        counted=Counter(tile_Set)
        if counted.values()==[2,2]:
            if tile_Set[0]==tile_Set[2] or tile_Set[1]==tile_Set[3]:
                return 2
            else:
                return 4
        else:
            return 4
    else:
        return 4

def Symmetry_Factor(tile_Set):
    if 0 in tile_Set:
        if tile_Set.count(0)==4:
            return 1
        if tile_Set.count(0)==1 or tile_Set.count(0)==3:
            return 4
        elif tile_Set[0]==tile_Set[2] or tile_Set[1]==tile_Set[3]:
            return 2
        else:
            return 4
    else:
        counted=Counter(tile_Set)
        if max(counted)==3:
            return 1
        else:
            return 1
        
def IKNOW():
    les=[]
    for i in xrange(8):
        for j in xrange(8):
            for k in xrange(8):
                for l in xrange(8):
                    tile=[i,j,k,l]
                    les.append(tile)

    return les
from copy import deepcopy
def make_Degenerates(tile_SetR,colours,depth,swapsR,pairs):
    #print depth,tile_Set,swaps
    tile_Set=deepcopy(tile_SetR)
    swaps=deepcopy(swapsR)
    if depth==0:
        for i in xrange(1,colours,2):
            if i in tile_Set and (i+1) in tile_Set:
                pairs.append(i+1)
                
    if depth>=len(tile_Set):
        print tile_Set
        return

    tile=tile_Set[depth]
    conjugate_tile=(tile % 2)*(tile+1)+(1-tile % 2)*(tile-1)
    if tile==0:
        make_Degenerates(tile_Set,colours,depth+1,swaps,pairs)
        tile_Set[depth]=7
        make_Degenerates(tile_Set,colours,depth+1,swaps,pairs)
    elif tile==7:
        make_Degenerates(tile_Set,colours,depth+1,swaps,pairs)
        tile_Set[depth]=0
        make_Degenerates(tile_Set,colours,depth+1,swaps,pairs)
    elif tile in swaps:
        tile_Set[depth]=swaps[tile]
        make_Degenerates(tile_Set,colours,depth+1,swaps,pairs)
    elif conjugate_tile in swaps:
        tile_Set[depth]=swaps[conjugate_tile]
        make_Degenerates(tile_Set,colours,depth+1,swaps,pairs)
    else:
        #make_Degenerates(tile_Set,colours,depth+1,swaps,pairs)
        new_Colour=tile
        
        while(new_Colour<colours-1):
            new_conjugate=(new_Colour % 2)*(new_Colour+1)+(1-new_Colour % 2)*(new_Colour-1)
            if False and new_Colour in pairs:
                
                swaps[tile]=new_Colour
                swaps[conjugate_tile]=new_conjugate
                make_Degenerates(tile_Set,colours,depth+1,swaps,pairs)
                new_Colour+=1
                continue
            else:
                swaps[tile]=new_Colour
                swaps[conjugate_tile]=new_conjugate
                tile_Set[depth]=new_Colour
                make_Degenerates(tile_Set,colours,depth+1,swaps,pairs)
                new_Colour+=1
            
def Make_Colour_Space(N):
    return (N-2)*(N-4)*(8+12*(N-6)+6*(N-6)*(N-8)+(N-6)*(N-8)*(N-10))


def Brutal():
    lineF=[line.rstrip('\n') for line in open('/scratch/asl47/Search_Spaces/SearchSpace_2_8.txt')]
    lineFP=[[int(ite) for ite in ln.split()] for ln in lineF]
    lineFPP=[]
    for line in lineFP:
        if 7 not in line:
            lineFPP.append(line)

    with open('/scratch/asl47/Search_Spaces/SearchSpace_2_7.txt','w') as numbering:
            for line in lineFPP:
                for part in line[:-1]:
                    numbering.write(str(part)+" ")
                numbering.write(str(line[-1]))
                numbering.write('\n')
                
    return lineFPP

def B22():
    num=0
    lineF=[line.rstrip('\n') for line in open('/rscratch/asl47/SearchSpace_2_7_Clean.txt')]
    lineFP=[[int(ite) for ite in ln.split()] for ln in lineF]
    for i,line in enumerate(lineFP):
        for j,line2 in enumerate(lineFP[i+1:]):
            if line==line2:
                #print "dup at ",i,j+i+1
                lineFP[j+i+1]=-1
                
        #if line.count(0)==len(line):
        #    continue
        #minim=min(i for i in line if i>0)
        #if minim >1:
        #    num+=1
        #    print line
    #print num

    with open('/scratch/asl47/Search_Spaces/SearchSpace_2_7_CleanM.txt','w') as numbering:
        for line in lineFP:
            if line==-1:
                continue
            for part in line[:-1]:
                numbering.write(str(part)+" ")
            numbering.write(str(line[-1]))
            numbering.write('\n')
                
    return lineFP

def B33():
    lineF=[line.rstrip('\n') for line in open('/rscratch/asl47/SearchSpace_2_7_CleanM.txt')]
    lineFP=[[int(ite) for ite in ln.split()] for ln in lineF]
    for i,line in enumerate(lineFP):
        #print line
        if 5 in line:
            if 3 in line:
                if 1 in line:
                    continue
                else:
                    lineFP[i]=-1
            else:
                lineFP[i]=-1
        elif 3 in line:
            if 1 in line:
                continue
            else:
                lineFP[i]=-1

    for i,line in enumerate(lineFP):
        if line==-1:
            continue
        #print line
        if 5 in line:
            if line.count(6)>(line.count(5)+1):
                lineFP[i]=-1
        if 3 in line:
            if line.count(4)>(line.count(3)+1):
                lineFP[i]=-1
        if 1 in line:
            if line.count(2)>(line.count(1)+1):
                print line
                lineFP[i]=-1
    with open('/scratch/asl47/Search_Spaces/SearchSpace_2_7_CleanMR.txt','w') as numbering:
        for line in lineFP:
            if line==-1:
                continue
            for part in line[:-1]:
                numbering.write(str(part)+" ")
            numbering.write(str(line[-1]))
            numbering.write('\n')

def B44():
    lineF=[line.rstrip('\n') for line in open('/rscratch/asl47/SearchSpace_2_7_CleanMRCSS.txt')]
    lineFP=[[int(ite) for ite in ln.split()] for ln in lineF]

    for i,line in enumerate(lineFP):
        if line==-1:
            continue
        for j,line2 in enumerate(lineFP[i+1:]):
            line3=line[4:]+line[:4]
            if line3==line2:
                lineFP[j+i+1]=-1
    
    for i,line in enumerate(lineFP):
        if line==-1:
            continue
        lineSw=[7 if t==1 else t for t in line]
        lineSw2=[1 if t==2 else t for t in lineSw]
        lineSw3=[2 if t==7 else t for t in lineSw2]
        lineSw4=[7 if t==3 else t for t in lineSw3]
        lineSw5=[3 if t==4 else t for t in lineSw4]
        lineSw6=[4 if t==7 else t for t in lineSw5]
        lineSw7=[8 if t==5 else t for t in lineSw6]
        lineSw8=[5 if t==6 else t for t in lineSw7]
        lineSw9=[6 if t==8 else t for t in lineSw8]
        for j,line2 in enumerate(lineFP[i+1:]):
            if lineSw9==line2:
                lineFP[j+i+1]=-1
            line3=lineSw9[4:]+lineSw9[:4]
            if line3==line2:
                lineFP[j+i+1]=-1
    with open('/scratch/asl47/Search_Spaces/SearchSpace_2_7_CleanMRCC.txt','w') as numbering:
        for line in lineFP:
            if line==-1:
                continue
            for part in line[:-1]:
                numbering.write(str(part)+" ")
            numbering.write(str(line[-1]))
            numbering.write('\n')

def B0():
    lineF=[line.rstrip('\n') for line in open('/rscratch/asl47/SearchSpace_2_7_CleanMRCSS.txt')]
    lineFP=[[int(ite) for ite in ln.split()] for ln in lineF]
    ML=[]
    for line in lineFP:
        lp1=line[:4]
        lp2=line[4:]

        if lp1.count(0)==2 and lp2.count(0)==2:
            ML.append(line)
    return ML


def Tempester():
    lineRaw=[line.rstrip('\n') for line in open('/rscratch/asl47/SearchSpace_3_8.txt')]
    lineF=[[int(i) for i in line if i!=" "] for line in lineRaw]
    print lineF[0]
    Ls=[]
    for line in lineF:
        if line[0]==line[1] and line[4]==line[0] and line[5]==line[0] and line[8]==line[0] and line[9]==line[0] and line[0]==0 and line.count(0)==6:
            
            Ls.append(line)
    print len(Ls)
    return Ls

def write_list(lis):
    with open('/rscratch/asl47/Temp_Cl.txt','w') as numbering:
        for line in lis:
            for part in line[:-1]:
                numbering.write(str(part)+" ")
            numbering.write(str(line[-1]))
            numbering.write('\n')
def hamming2(s1, s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def roul():
    x=[0.92]
    y=[1]
    z=x*990+y*10
    return z
    

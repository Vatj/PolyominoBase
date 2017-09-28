from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from math import ceil
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import icy
icy.Use_Seaborn()
from collections import Counter


def LoadData(t=3,run='x'):
    details='{}_O{}'.format(t,run) if t==4 else '3'
    shape_data_raw=[line.rstrip('\n').split('*') for line in open('/rscratch/asl47/Bulk_Run/Modular/Shape_Data_T{}.txt'.format(details))]
    shape_data_proccessed=[]
    for shape in shape_data_raw:
        shape_data_proccessed.append([int(s) for s in shape[4].split()+shape[6].split()+shape[8].split()+shape[10].split()])
    return shape_data_proccessed

def LoadLengths(t,o):
    details='{}_O{}'.format(t,o) if t==4 else '3'
    return np.genfromtxt('/rscratch/asl47/Bulk_Run/Modular/Shape_Data_T{}_L.txt'.format(details))

def PlotCCLengths(data,label_in):
    cnt=Counter(data)
    
    plt.loglog(sorted(cnt.keys()),[cnt[k]*1./sum(cnt.values()) for k in sorted(cnt.keys())],lw=2,ls='-',nonposy='mask',label=label_in)
        
    plt.xlabel('CC Size')
    plt.ylabel('Normalised Frequency')
    plt.legend(ncol=2)
    
    plt.show(block=False)

def PlotShapeSizes(data,ret=False):
    sizes=[sum(shape[2:-2]) for shape in data]
    sym_sizes=[list() for i in xrange(4)]
    colours=['forestgreen','darkorange','mediumslateblue','k']
    for shape in data:
        if shape[-2]==0:
            sym_sizes[3].append(sum(shape[2:-2]))
        else:
            sym_sizes[shape[-1]/2].append(sum(shape[2:-2]))

    if ret:
        return sym_sizes
    plt.figure()
    (h,b,c)=plt.hist(sym_sizes,range=(min(sizes),max(sizes)),bins=max(sizes)-min(sizes)+1,normed=True,stacked=True,color=colours,label=[r'$C_0$',r'$C_2$',r'$C_4$','partial'])

    plt.xlabel('Shape Size')
    plt.ylabel('Count')
    plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.show(block=False)

def LoadShapeDistrs(Ts,Os):
    keys=['T{}{}'.format(t,'' if t==3 else '_O{}'.format(o)) for (t,o) in zip(Ts,Os)]
    data_dict={}
    for (key,t,o) in zip(keys,Ts,Os):
        data=PlotShapeSizes(LoadData(t,o),True)
        flat_sizes=[item for sublist in data[:3] for item in sublist]
        data_dict[key]=flat_sizes    
    return data_dict

def PlotShapeDistrs(data_dict):
    plt.figure()
    sorted_keys=sorted(data_dict.keys())
    sorted_keys[1:]=sorted(sorted_keys[1:],key=lambda x: int(x[4:]))
    for key in sorted_keys:
        #hist_vals,bin_edges=np.histogram(data,
        #plt.hist(data,range=(0,100),bins=100)
        #print key
        cnt=Counter(data_dict[key])
        
        plt.loglog(sorted(cnt.keys()),smooth([cnt[k]*1./sum(cnt.values()) for k in sorted(cnt.keys())],1),lw=2,ls='-',nonposy='mask',label=r'{}'.format('Dynamic ($\Omega$ {})'.format(int(key[4:])) if key[1]=='4' else 'Static'))
        
    plt.xlabel('Phenotype Size')
    plt.ylabel('Normalised Frequency')
    plt.legend(ncol=2)
    
    plt.show(block=False)
    
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
                    
def Visualise_Shape_From_Binary(Shape_String_X,count,currentAxis,X_Y_LIMS,tightBounded=False,num=""):
    d={-10:"Connected UND",-5:"Disjointed \n UND",-4:"Disjointed \n Divergence",-1:"Steric"}
    if type(Shape_String_X)==int:
        
        currentAxis.text(.5, .5,d[Shape_String],horizontalalignment='center',verticalalignment='top',transform = currentAxis.transAxes)

    else:
        Shape_String=Shape_String_X[:-2]
        valid=Shape_String_X[-2]
        DELTA_X=Shape_String[0]
        DELTA_Y=Shape_String[1]
        for i in xrange(DELTA_X):
            for j in xrange(DELTA_Y):
                if(Shape_String[2+i+j*DELTA_X]):
                    if valid:
                        currentAxis.add_patch(Rectangle((i-DELTA_X/2., DELTA_Y/2.-j-1), 1, 1, facecolor='lightsteelblue',edgecolor='slategrey',fill=True,hatch='////',lw=0))
                        currentAxis.add_patch(Rectangle((i-DELTA_X/2., DELTA_Y/2.-j-1), 1, 1,edgecolor='maroon',fill=False,lw=2.5))
                    else:
                        currentAxis.add_patch(Rectangle((i-DELTA_X/2., DELTA_Y/2.-j-1), 1, 1, facecolor='lightgray',edgecolor='slategrey',fill=True,hatch='////',lw=0))
                        currentAxis.add_patch(Rectangle((i-DELTA_X/2., DELTA_Y/2.-j-1), 1, 1,edgecolor='k',fill=False,lw=2.5))
                        

                    
    currentAxis.text(0, 1,"{})".format(num),horizontalalignment='left',verticalalignment='top',transform = currentAxis.transAxes)
    if type(count)==int:
        currentAxis.text(.5, 0,"{}".format(count),horizontalalignment='center',verticalalignment='top',transform = currentAxis.transAxes)
    else:
        currentAxis.text(.5, 0,"{:.3G}%".format(count),horizontalalignment='center',verticalalignment='top',transform = currentAxis.transAxes)
    if tightBounded and type(Shape_String)==list:
        maxB=max(DELTA_X,DELTA_Y)
        currentAxis.set_xlim([-maxB/2.-0.5,maxB/2.+0.5])
        currentAxis.set_ylim([-maxB/2.-0.5,maxB/2.+0.5])

        
        #currentAxis.set_xlim([-DELTA_X/2.-0.5,DELTA_X/2.+0.5])
        #currentAxis.set_ylim([-DELTA_Y/2.-0.5,DELTA_Y/2.+0.5])
    else:
        currentAxis.set_xlim([-X_Y_LIMS,X_Y_LIMS])
        currentAxis.set_ylim([-X_Y_LIMS,X_Y_LIMS])
    currentAxis.set_aspect('equal')
    currentAxis.set_axis_off()
    

    
def Make_Shape_Page(Shape_Strings,Shape_Counts,X_rows,Y_rows,X_Y_LIMS,saving,pdf,pgN):
    fig, axarr = plt.subplots(X_rows, Y_rows,figsize=(8.27,11.69))
    temInd=0
    for shape,count,axis in zip(Shape_Strings,Shape_Counts,axarr.flatten()):
        Visualise_Shape_From_Binary(shape,count,axis,X_Y_LIMS,True,pgN*X_rows*Y_rows+temInd)
        temInd+=1
    emptyAxis=X_rows*Y_rows-len(Shape_Strings)
    if emptyAxis>0:
        for axis in axarr.flatten()[len(Shape_Strings):]:
            axis.set_axis_off()
    plt.subplots_adjust(wspace=0,hspace=0,top=1,left=0,bottom=0.02,right=1)
    if(saving):
        pdf.savefig(fig)
        #plt.close('all')
    else:
        plt.show(block=False)
    

def VisualiseShapesNew(shapes,saving=False,fileSave='Test',valid_only=True):
    X_rows=4
    Y_rows=5
    fname="Output/{}.pdf".format(fileSave)
    pdf = PdfPages(fname)
    pgN=0
    if valid_only:
        print "cutting"
        shapes.sort(key=lambda x: x[-2])
        gener = (i for i,v in enumerate(shapes) if v[-2]==1)
        indexer=gener.next()
        print indexer
        shapes=shapes[:100]#indexer:]
    shapes.sort(key=lambda x: sum(x[2:-2]))
    #return shapes
    for pg in xrange(int(ceil(len(shapes)/(X_rows*Y_rows*1.)))):
        Make_Shape_Page(shapes[X_rows*Y_rows*pgN:X_rows*Y_rows*(pgN+1)],[sum(shape[2:-2]) for shape in shapes][X_rows*Y_rows*pgN:X_rows*Y_rows*(pgN+1)],X_rows,Y_rows,10,saving,pdf,pgN)
        pgN+=1
    pdf.close()

def Load_Shapes():
    shape_dict={int(key): [int(info) for info in deltas.split()+shape.split()] for (key,deltas,shape) in [line.rstrip('\n').split('*')[2::2] for line in open('/rscratch/asl47/Modular_Shapes.txt')]}
    print "done"
    print shape_dict
    return shape_dict

def TargetShapes(targets,labels):
    #defaults
    square=[4,4, -1,2,0,-1, 0,1,1,-1, -1,1,1,0, -1,0,-1,-1]
    cross=[6,5, -1,-1,-1,0,-1,-1, -1,0,0,1,0,-1, 0,1,1,1,1,2, -1,0,0,1,0,-1, -1,-1,-1,0,-1,-1]
    tetris=[5,4, -1,0,0,-1,-1,0,1,1,0,-1,-1,0,1,1,2,-1,-1,0,0,-1]

    targets=[square,cross,tetris]
    labels=['square','cross','tetris']
    
    fig, axarr = plt.subplots(len(targets)+1, 1,figsize=(8.27,11.69))
    currentAxis=axarr[0]
    
    currentAxis.add_patch(Rectangle((0,0), 1, 1, facecolor='red',edgecolor='slategrey',fill=True,hatch='////',lw=0))
    currentAxis.add_patch(Rectangle((0,2), 1, 1, facecolor='firebrick',edgecolor='slategrey',fill=True,lw=0,alpha=0.5))  
    currentAxis.add_patch(Rectangle((0,6), 1, 1, facecolor='green',edgecolor='slategrey',fill=True,lw=0,alpha=0.2))
    currentAxis.add_patch(Rectangle((0,4), 1, 1, facecolor='blue',edgecolor='slategrey',fill=True,lw=0))
    currentAxis.add_patch(Rectangle((0,0), 1, 1,edgecolor='darkgray',fill=False,lw=2.5))
    currentAxis.add_patch(Rectangle((0,2), 1, 1,edgecolor='darkgray',fill=False,lw=2.5))
    currentAxis.add_patch(Rectangle((0,6), 1, 1,edgecolor='darkgray',fill=False,lw=2.5))
    currentAxis.add_patch(Rectangle((0,4), 1, 1,edgecolor='darkgray',fill=False,lw=2.5))

    currentAxis.text(1.2, 0.5,"Occupied",horizontalalignment='left',verticalalignment='center')
    currentAxis.text(1.2, 2.5,"Binding",horizontalalignment='left',verticalalignment='center')
    currentAxis.text(1.2, 6.5,"Ignored",horizontalalignment='left',verticalalignment='center')
    currentAxis.text(1.2, 4.5,"Empty",horizontalalignment='left',verticalalignment='center')

    currentAxis.set_axis_off()
    currentAxis.set_xlim([-1,3])
    currentAxis.set_ylim([-1,8])
    
    for count,(target,label) in enumerate(zip(targets,labels)):
        currentAxis=axarr[count+1]
        currentAxis.set_axis_off()
        delta_x=target[0]
        delta_y=target[1]
        maxB=max(delta_x,delta_y)
        currentAxis.set_xlim([-delta_x/2.-0.5,delta_x/2.+0.5])
        currentAxis.set_ylim([-delta_y/2.-0.5,delta_y/2.+0.5])
        
        currentAxis.text(0, delta_y/2.+0.5,label,horizontalalignment='center',verticalalignment='center')
        for i in xrange(delta_y):
            for j in xrange(delta_x):
                if target[2+j+i*delta_x]==1:
                    currentAxis.add_patch(Rectangle((j-delta_x/2., delta_y/2.-i-1), 1, 1, facecolor='red',edgecolor='slategrey',fill=True,hatch='////',lw=0))
                elif target[2+j+i*delta_x]==2:
                    currentAxis.add_patch(Rectangle((j-delta_x/2., delta_y/2.-i-1), 1, 1, facecolor='firebrick',edgecolor='slategrey',fill=True,lw=0,alpha=0.5))  
                elif target[2+j+i*delta_x]==-1:
                    currentAxis.add_patch(Rectangle((j-delta_x/2., delta_y/2.-i-1), 1, 1, facecolor='green',edgecolor='slategrey',fill=True,lw=0,alpha=0.2))
                else:
                    currentAxis.add_patch(Rectangle((j-delta_x/2., delta_y/2.-i-1), 1, 1, facecolor='blue',edgecolor='slategrey',fill=True,lw=0))
                    
                currentAxis.add_patch(Rectangle((j-delta_x/2., delta_y/2.-i-1), 1, 1,edgecolor='darkgray',fill=False,lw=2.5))

    plt.show(block=False)
    
def sortShapes(Shapes,Counts):
    if not -10 in Shapes:
        Shapes.append(-10)
    if not -5 in Shapes:
        Shapes.append(-5)
    if not -4 in Shapes:
        Shapes.append(-4)
    if not -1 in Shapes:
        Shapes.append(-1)
    np_counts=np.array(Counts.values())
    inds=np_counts.argsort()
    np_counts=np_counts[inds]
    sortedShapes=[Shapes[i] for i in inds]
    
    return sortedShapes[::-1],np_counts[::-1]
    
def Visualise_Shapes(Shape_Strings,Shape_Counts,X_Y_LIMS,saving=False,fileSave="Test"):
    sortShape,np_count=sortShapes(Shape_Strings,Shape_Counts)
    X_rows=5
    Y_rows=4
    totalC=sum(Shape_Counts.values())*1.
    fname="Output/{}.pdf".format(fileSave)
    pdf = PdfPages(fname)
    pgN=0
    for pg in xrange(int(ceil(len(Shape_Strings)/(X_rows*Y_rows*1.)))):
        Make_Shape_Page(sortShape[pg*(X_rows*Y_rows):(pg+1)*(X_rows*Y_rows)],np_count[pg*(X_rows*Y_rows):(pg+1)*(X_rows*Y_rows)]/totalC*100,X_rows,Y_rows,X_Y_LIMS,saving,pdf,pgN)
        pgN+=1
    pdf.close()


    

def Clockwise_Rotation(Shape_String):
    rotated_String=Shape_String[1::-1]
    for column in xrange(Shape_String[0]):
        for row in xrange(Shape_String[1]-1,-1,-1):
            rotated_String.append(Shape_String[2+row*Shape_String[0]+column])
    return rotated_String



def Merge_Shape_Data_Files(runType,runNum,numParts):
    
    Data_Base_Fname="Output/S_{}_Run_R{}P{}_Data.txt".format(runType,runNum,1)
    Shape_Base_Fname="Output/S_{}_Run_R{}P{}_Shapes.txt".format(runType,runNum,1)
    Shape_Counts_Base={int(key):int(value) for (key,value) in [line.rstrip('\n').split("Code:")[-1].split(" N:") for line in open(Data_Base_Fname)]}
    

    Shape_Strings_Base=[[int(j) for j in u if j!=' '] for u in [q[0]+q[1] for q in [line.rstrip('\n').split("Deltas ")[-1].split("  Shape: ") for line in open(Shape_Base_Fname)]]]

    for i in xrange(2,numParts+1):
        print "Merging on part ",i
        Data_Fname="Output/S_{}_Run_R{}P{}_Data.txt".format(runType,runNum,i)
        Shape_Fname="Output/S_{}_Run_R{}P{}_Shapes.txt".format(runType,runNum,i)
        Shape_Strings=[[int(j) for j in u if j!=' '] for u in [q[0]+q[1] for q in [line.rstrip('\n').split("Deltas ")[-1].split("  Shape: ") for line in open(Shape_Fname)]]]
        Shape_Counts={int(key):int(value) for (key,value) in [line.rstrip('\n').split("Code:")[-1].split(" N:") for line in open(Data_Fname)]}


        for key,value in Shape_Counts.iteritems():
            if key>=0: 
                continue
            else:
                Shape_Counts_Base[key]+=value;
            
    
        for ind,Shape_String in enumerate(Shape_Strings):
            if Shape_String in Shape_Strings_Base:
                Shape_Counts_Base[Shape_Strings_Base.index(Shape_String)]+=Shape_Counts[ind]
            else:
                rotation=1
                while(rotation<4):
                    Shape_String=Clockwise_Rotation(Shape_String)
                    if Shape_String in Shape_Strings_Base:
                        Shape_Counts_Base[Shape_Strings_Base.index(Shape_String)]+=Shape_Counts[ind]
                        break
                    else:
                        rotation+=1
                else:
                    Shape_String=Clockwise_Rotation(Shape_String)
                    Shape_Counts_Base[len(Shape_Strings_Base)]=Shape_Counts[ind]
                    Shape_Strings_Base.append(Shape_String)
    
    return Shape_Counts_Base,Shape_Strings_Base

def find_Difference(stringSet1, stringSet2):
    for string in stringSet1:
        if string in stringSet2:
            continue
        else:
            rotation=1
            while(rotation<4):
                string=Clockwise_Rotation(string)
                if string in stringSet2:
                    break
                else:
                    rotation+=1
            else: #is different
                print string

                

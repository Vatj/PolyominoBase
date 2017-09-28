#!/usr/bin/env python

#import ctypes
from copy import copy
import sys

from collections import defaultdict
from random import choice
from random import uniform

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.animation import FuncAnimation,ImageMagickWriter

## UTILS ##

def cycleList(inputList,numCycles):
    return inputList[-numCycles:] + inputList[:-numCycles]

## WRAPPER SECTION ##

#Poly_Lib= ctypes.cdll.LoadLibrary('./Polyomino_Library.so')
#Poly_Lib.Graph_Assembly_Outcome.restype=ctypes.c_int
#Poly_Lib.Graph_Assembly_Outcome.argtype=[ctypes.c_int,ctypes.POINTER(ctypes.c_int)]

def GraphAssemblyOutcome(genotype):
    genotype_Pointer=(ctypes.c_int*len(genotype))(*genotype)
    return Poly_Lib.Graph_Assembly_Outcome(len(genotype),genotype_Pointer)

def TimeIt(tops):
    for top in tops:
        Graph_Assembly_Outcome(top)

## POLYOMINO BUILDER ##

def InteractionMatrix(input_face):
    return  (1-input_face%2)*(input_face-1)+(input_face%2)*(input_face+1) if input_face>0 else input_face

def PolyominoBuilder(genotype):
    SIZE_LIMIT=(len(genotype)/2)**2
    POLYOMINO_GRID=defaultdict(tuple)
    POSSIBLE_GRID=defaultdict(list)
    IMPOSSIBLE_GRID=set()
    TILE_TYPES=[genotype[i:i+4] for i in xrange(0, len(genotype), 4)]

    def placeTile(tType,position,orientation):
        POLYOMINO_GRID[position]=(tType,orientation)
        return position,(tType,orientation)

    def identifyValidNeighbours(position):
        centerType,centerOrientation=POLYOMINO_GRID[position]
        identifyValidNeighbour(position,centerType,centerOrientation,(position[0],position[1]+1),0)#Check Top
        identifyValidNeighbour(position,centerType,centerOrientation,(position[0]+1,position[1]),1)#Check Right
        identifyValidNeighbour(position,centerType,centerOrientation,(position[0],position[1]-1),2)#Check Bottom
        identifyValidNeighbour(position,centerType,centerOrientation,(position[0]-1,position[1]),3)#Check Left
                
    def identifyValidNeighbour(position,centerType,centerOrientation,checkPosition,increment):
        if checkPosition in POLYOMINO_GRID:
            return False
        bindingEdge=cycleList(TILE_TYPES[centerType],centerOrientation)[increment]
        oppositeBindingEdgeIndex=(increment+2)%4
        for i,tile in enumerate(TILE_TYPES):
            for cycNum in xrange(4):
                if bindingEdge!=0 and cycleList(tile,cycNum)[oppositeBindingEdgeIndex]==InteractionMatrix(bindingEdge):
                    POSSIBLE_GRID[checkPosition].append((i,cycNum))

    placement=placeTile(0,(0,0),0)
    identifyValidNeighbours((0,0))
    yield placement,copy(POSSIBLE_GRID)
    while len(POSSIBLE_GRID)>0:
        newPolyominoPosition,newPolyominoDetails=choice([(position, tileDetail) for position, tileDetails in POSSIBLE_GRID.iteritems() for tileDetail in tileDetails])
        POSSIBLE_GRID.pop(newPolyominoPosition)
        placement= placeTile(newPolyominoDetails[0],newPolyominoPosition,newPolyominoDetails[1])
        identifyValidNeighbours(newPolyominoPosition)
        if len(POLYOMINO_GRID)>SIZE_LIMIT:
            return
        else:
            yield placement,copy(POSSIBLE_GRID)
        
## ANIMATION SECTION ##

def GrowPoly(genotype,write_it=False):
    #if GraphAssemblyOutcome(genotype)<=0:
    #    print 'bad build, should reject...' 
    fig = plt.figure(figsize=(10,10))
    plt.axis('off')
    ax = plt.gca()
    ax.set_aspect(1)
    COLORS=['royalblue','darkgreen','firebrick','chocolate','orchid','goldenrod']
    HATCHES=['//','\\','+','O', '.']

    def init():
        pass
    def AnimateBuild(i,data,ft,tt):
        ## FIRST FRAME ##
        if i==0:
            ft.append(Rectangle((0,0),0.95,0.95,fill=True,alpha=1,facecolor=COLORS[0],edgecolor='k',lw=2,hatch='*'))
            ax.add_patch(ft[-1])
            ft.append(Rectangle((0,0),0.95,0.95,fill=False,alpha=1,edgecolor='r',lw=2))
            ax.add_patch(ft[-1])
            return ft+tt
        ## SECOND FRAME ##
        if i==1:
            ft[1].set_visible(False)
            return ft+tt
        ## THIRD FRAME ##
        if i==2:
            for key in data[0][1]:
                potential_tiles=[t_type[0] for t_type in data[0][1][key]]
                for j,pt in enumerate(set(potential_tiles)):
                    tt.append(Rectangle(key,0.95,0.95,fill=True,alpha=0.25,facecolor=COLORS[pt],edgecolor='k',lw=2,hatch=HATCHES[j%5]))
                    ax.add_patch(tt[-1])
            return ft+tt
        ## LAST FRAME ##
        if i==len(data)*2+1:
            ft[-1].set_visible(False)
            return ft+tt
        ## !!FURTHER FRAMES!! ##
        if i%2==0:
            for del_it in ft[1::2]:
                del_it.set_visible(False)
            for del_it in tt:
                del_it.set_visible(False)
            for key in data[i/2-1][1]:
                potential_tiles=[t_type[0] for t_type in data[i/2-1][1][key]]
                for j,pt in enumerate(set(potential_tiles)):
                    tt.append(Rectangle(key,0.95,0.95,fill=True,alpha=0.25,facecolor=COLORS[pt],edgecolor='k',lw=2,hatch=HATCHES[j%5]))
                    ax.add_patch(tt[-1])
            return ft+tt
        if i%2==1:
            current_tile=data[i/2][0]
            ft.append(Rectangle(current_tile[0],0.95,0.95,fill=True,alpha=1,facecolor=COLORS[current_tile[1][0]],edgecolor='k',lw=2))
            ax.add_patch(ft[-1])
            ft.append(Rectangle(current_tile[0],0.95,0.95,fill=False,alpha=1,edgecolor='r',lw=2))
            ax.add_patch(ft[-1])
            return ft+tt

    fixed_tiles=[]
    temporary_tiles=[]
    poly_generator=list(PolyominoBuilder(genotype))
    plt.axis([min([i[0][0][0] for i in poly_generator])-0.25,max([i[0][0][0] for i in poly_generator])+1.25,min([i[0][0][1] for i in poly_generator])-0.25,max([i[0][0][1] for i in poly_generator])+1.25])
    anim = FuncAnimation(fig, AnimateBuild,init_func=init,frames=len(poly_generator)*2+2, interval=800, blit=False,fargs=(poly_generator,fixed_tiles,temporary_tiles),repeat=False)
    plt.tight_layout()
    if write_it:
        writer = ImageMagickWriter(fps=1.25)
        anim.save('PolyominoAnimation.gif', writer=writer)
    else:
        plt.show(block=False)
        
## MAIN SECTION ##

def main():
    if len(sys.argv)==1 or sys.argv[1]=='-H' or sys.argv[1]=='-h':
        print "**Polyomino Animator**\nParameter is the genotype as a space separated sequence\nExample: 'python polyomino_animator.pyc 1 2 3 3 4 0 0 0'\nResulting output is saved as a gif in the local directory"
    else:
        genotype=[int(i) for i in sys.argv[1:]]
        print "genotype is: ",genotype
        assert(len(genotype)%4==0),  "Poorly formed genotype, see example by running 'python polyomino_animator.pyc -h'"
        GrowPoly(genotype,1)
        print "finished!"

if __name__ == "__main__":
    main()

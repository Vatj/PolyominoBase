#include "graph_analysis.hpp"
//#include <exception>

int Graph_Analysis_Fast(std::vector<int> genotype) {
  std::vector<int> Null_Check{0,0,0,0};
  if(genotype==Null_Check)
      return 1;
  if(Double_BP_Check(genotype))
    return -1;
  int SIF_result=SIF_Elimination(genotype);
  if(SIF_result!=0)
    return SIF_result;
  if(BP_Check(genotype))
    return -1;
  int Bound_Limit=0;
  int loop_result=loop_Analysis(genotype,Bound_Limit);
  return loop_result>0 ? 1 : -1;
}

int FastGraphAnalysis(std::vector<int> genotype) {
  std::vector<int> Null_Check{0,0,0,0};
  if(genotype==Null_Check)
      return 1;
  if(Double_BP_Check(genotype))
    return -1;
  int SIF_result=SIF_Elimination(genotype);
  if(SIF_result!=0)
    return SIF_result;  
  if(BP_Disjointed_Check(genotype))
    return -1;
  if(BP_Check(genotype))
    return -1;
  return FastLoopAnalysis(genotype); 
  
}



int Graph_Analysis(std::vector<int> genotype) {
  std::vector<int> Null_Check{0,0,0,0};
  if(genotype==Null_Check)
      return 1;
  
  if(Double_BP_Check(genotype))
    return -2;
  
  int SIF_result=SIF_Elimination(genotype);
  if(SIF_result!=0)
    return SIF_result;  

  if(BP_Disjointed_Check(genotype))
    return -4;
  
  int Bound_Limit=0;
  int loop_result=loop_Analysis(genotype,Bound_Limit);

  if(loop_result<0)
    return loop_result;
  if(loop_result>0 && !BP_Check(genotype))
    return 1;  

  Bound_Limit+=Bound_Limit >0 ? genotype.size()/4 : genotype.size()/2;
  return Bounded_Loop_Growth_Check(genotype,Bound_Limit) ? 1 : -13;  
}





bool Bounded_Loop_Growth_Check(std::vector<int>& genome,int Bound_Limit) { 
  std::vector<int> tilesInformation{0,0,0,0,-1};
  std::vector<int> processingQueue{0,0,0,0,-10};
  int currentX=0, currentY=0, currentTile=0, currentOrientation=0, comingFrom=-1;
  bool superBREAK=false;
  std::vector<bool> Tile_Use(genome.size()/4,false);
  Tile_Use[0]=true;
  while(!processingQueue.empty()) {
    comingFrom=processingQueue.back(); processingQueue.pop_back();
    currentOrientation=processingQueue.back(); processingQueue.pop_back();
    currentTile=processingQueue.back(); processingQueue.pop_back();
    currentY=processingQueue.back(); processingQueue.pop_back();
    currentX=processingQueue.back();  processingQueue.pop_back();

    for(int f=0;f<4;++f) {
      superBREAK=false;
      if((f+2)%4==comingFrom)
        continue;
      int leadingEdge=genome[currentTile*4+(f-currentOrientation+4)%4];
      if(!leadingEdge)
        continue;
      int bondingEdge=Interaction_Matrix(leadingEdge);
      int foundIndex=std::find(genome.begin(),genome.end(),bondingEdge)-genome.begin();
      int tileN=foundIndex/4;
      int faceN=foundIndex%4;
      int tileO=((f+2)%4-faceN +4)%4;
      int X_OFFSET,Y_OFFSET;
      switch(f) {
      case 0:
        X_OFFSET=0;
        Y_OFFSET=1;
        break;
      case 1:
        X_OFFSET=1;
        Y_OFFSET=0;
        break;
      case 2:
        X_OFFSET=0;
        Y_OFFSET=-1;
        break;
      case 3:
        X_OFFSET=-1;
        Y_OFFSET=0;
        break;
      }
      for(std::vector<int>::const_iterator it = tilesInformation.begin();it!=tilesInformation.end();it+=5) {
        if(*it==currentX+X_OFFSET && *(it+1)==currentY+Y_OFFSET) {
          if(*(it+2)==tileN) {
            if(*(it+3)==tileO) {
              superBREAK=true; 
              break;
            }
            else {
              for(int test_face=0; test_face<4;++test_face) {
                if(genome[*(it+2)*4+(*(it+3)+test_face)%4]!=genome[tileN*4+(tileO+test_face)%4])
                  return false;
              }
              superBREAK=true; 
              break;
            }
          }            
          else
            return false;
        }
      }
      if(superBREAK)
        continue;    
      if(!Tile_Use[tileN])
        Tile_Use[tileN]=true;
      tilesInformation.insert(tilesInformation.end(),{currentX+X_OFFSET,currentY+Y_OFFSET,tileN,tileO,f});
      processingQueue.insert(processingQueue.end(),{currentX+X_OFFSET,currentY+Y_OFFSET,tileN,tileO,f});
      if(tilesInformation.size()/5>(unsigned int)Bound_Limit) {
        return false;
      }
    } 
  }
  return std::all_of(Tile_Use.begin(),Tile_Use.end(),[](bool v) { return v; });
}


double Steric_Check(std::vector<int>& genome,int Shape_Based_Fitness) {
  std::vector<int> phenotype_Shape;
  int DELTA_X=0, DELTA_Y=0;
  std::vector<int> tilesInformation{0,0,0,0,-1};
  std::vector<int> processingQueue{0,0,0,0,-10};
  std::vector<int> occupied_spiral(genome.size()*genome.size()/4);occupied_spiral[0]=1;
  int currentX=0, currentY=0, currentTile=0, currentOrientation=0, comingFrom=-1;
  bool superBREAK=false;
  while(!processingQueue.empty()) {
    comingFrom=processingQueue.back();
    processingQueue.pop_back();
    currentOrientation=processingQueue.back();
    processingQueue.pop_back();
    currentTile=processingQueue.back();
    processingQueue.pop_back();
    currentY=processingQueue.back();
    processingQueue.pop_back();
    currentX=processingQueue.back();
    processingQueue.pop_back();

    for(int f=0;f<4;++f) {
      superBREAK=false;
      if((f+2)%4==comingFrom)
        continue;
      int leadingEdge=genome[currentTile*4+(f-currentOrientation+4)%4];
      if(!leadingEdge)
        continue;

      int bondingEdge=Interaction_Matrix(leadingEdge);
      int foundIndex=std::find(genome.begin(),genome.end(),bondingEdge)-genome.begin();
      int tileN=foundIndex/4;
      int faceN=foundIndex%4;
      int tileO=((f+2)%4-faceN +4)%4;

      int X_OFFSET,Y_OFFSET;
      switch(f) {
      case 0:
        X_OFFSET=0;
        Y_OFFSET=1;
        break;
      case 1:
        X_OFFSET=1;
        Y_OFFSET=0;
        break;
      case 2:
        X_OFFSET=0;
        Y_OFFSET=-1;
        break;
      case 3:
        X_OFFSET=-1;
        Y_OFFSET=0;
        break;
      }
      
      //occupied_spiral[SpiralCoordinate(currentX+X_OFFSET,currentY+Y_OFFSET)]==1) { //spot is already occupied
      for(std::vector<int>::const_iterator it = tilesInformation.begin();it!=tilesInformation.end();it+=5) {
        if(*it==currentX+X_OFFSET && *(it+1)==currentY+Y_OFFSET) {
          if(*(it+2)==tileN) {
            if(*(it+3)==tileO) {
              superBREAK=true; 
              break;
            }
            else {
              for(int test_face=0; test_face<4;++test_face) {
                if(genome[*(it+2)*4+(*(it+3)+test_face)%4]!=genome[tileN*4+(tileO+test_face)%4])
                  return -1;
              }
              superBREAK=true; 
              break;
            }
          }            
          else
            return -1;
        }
      }
      
      if(superBREAK)
        continue;
      tilesInformation.insert(tilesInformation.end(),{currentX+X_OFFSET,currentY+Y_OFFSET,tileN,tileO,f});
      //occupied_spiral[SpiralCoordinate(currentX+X_OFFSET,currentY+Y_OFFSET)]=1;
      processingQueue.insert(processingQueue.end(),{currentX+X_OFFSET,currentY+Y_OFFSET,tileN,tileO,f});
      if(tilesInformation.size()>1.25*genome.size()*genome.size()) {
#pragma omp critical(ster)
        {
          std::cout<<"HARD EXIT, STERIC OVERFLOW from"<<std::endl;
          for(auto& g: genome) {
            std::cout<<g<<" ";
          }
          std::cout<<std::endl;
        }
        return -1;        
      }
    }    
  }
  if(Shape_Based_Fitness<0)
    return tilesInformation.size()/5;
  else {
    std::vector<int> X_LOCs, Y_LOCs;               
    for(auto it = tilesInformation.begin();it!=tilesInformation.end();it+=5) {
      X_LOCs.emplace_back(*it);
      Y_LOCs.emplace_back(*(it+1));
    }
    std::vector<int>::iterator LEFT_X,RIGHT_X,TOP_Y,BOTTOM_Y;
  
    std::tie(LEFT_X,RIGHT_X)=std::minmax_element(X_LOCs.begin(),X_LOCs.end());
    std::tie(BOTTOM_Y,TOP_Y)=std::minmax_element(Y_LOCs.begin(),Y_LOCs.end());
  
    DELTA_X=*RIGHT_X-*LEFT_X+1;
    DELTA_Y=*TOP_Y-*BOTTOM_Y+1;

    std::vector<int> Spatial_Occupation(DELTA_X*DELTA_Y);
    phenotype_Shape.assign(((DELTA_X+2)*(DELTA_Y+2)),0);
    std::vector<int> Spatial_Grid(DELTA_X*DELTA_Y);
    int UNIQUE_TILES=0;
    for(unsigned int tileIndex=0;tileIndex<X_LOCs.size();++tileIndex) {
      Spatial_Occupation[(*TOP_Y-Y_LOCs[tileIndex])*DELTA_X + (X_LOCs[tileIndex]-*LEFT_X)]=1;
      ++Spatial_Grid[(*TOP_Y-Y_LOCs[tileIndex])*DELTA_X + (X_LOCs[tileIndex]-*LEFT_X)];
      if(Spatial_Grid[(*TOP_Y-Y_LOCs[tileIndex])*DELTA_X + (X_LOCs[tileIndex]-*LEFT_X)]==1) {
        ++UNIQUE_TILES;
      }
    }
    for(int x_shift=0;x_shift<DELTA_X;++x_shift) {
      for(int y_shift=0;y_shift<DELTA_Y;++y_shift) {
        phenotype_Shape[(DELTA_X+2)+1 + x_shift+y_shift*(DELTA_X+2)]=Spatial_Occupation[x_shift+y_shift*DELTA_X];
      }
    }
    DELTA_X+=2;
    DELTA_Y+=2;
    return Shape_Matching_Fitness_Function(phenotype_Shape,DELTA_X,DELTA_Y,Shape_Based_Fitness);
  } 
}

std::vector<int> Steric_Check_Occupation(std::vector<int>& genome) {
  std::vector<int> phenotype_Shape;
  int DELTA_X=0, DELTA_Y=0;
  std::vector<int> tilesInformation{0,0,0,0,-1};
  std::vector<int> processingQueue{0,0,0,0,-10};
  int currentX=0, currentY=0, currentTile=0, currentOrientation=0, comingFrom=-1;
  bool superBREAK=false;//, advanced_Steric_Check=false;
  while(!processingQueue.empty()) {
    comingFrom=processingQueue.back();
    processingQueue.pop_back();
    currentOrientation=processingQueue.back();
    processingQueue.pop_back();
    currentTile=processingQueue.back();
    processingQueue.pop_back();
    currentY=processingQueue.back();
    processingQueue.pop_back();
    currentX=processingQueue.back();
    processingQueue.pop_back();

    for(int f=0;f<4;++f) {
      superBREAK=false;
      if((f+2)%4==comingFrom)
        continue;
      int leadingEdge=genome[currentTile*4+(f-currentOrientation+4)%4];
      if(!leadingEdge)
        continue;

      int bondingEdge=Interaction_Matrix(leadingEdge);
      int foundIndex=std::find(genome.begin(),genome.end(),bondingEdge)-genome.begin();
      int tileN=foundIndex/4;
      int faceN=foundIndex%4;
      int tileO=((f+2)%4-faceN +4)%4;

      int X_OFFSET,Y_OFFSET;
      switch(f) {
      case 0:
        X_OFFSET=0;
        Y_OFFSET=1;
        break;
      case 1:
        X_OFFSET=1;
        Y_OFFSET=0;
        break;
      case 2:
        X_OFFSET=0;
        Y_OFFSET=-1;
        break;
      case 3:
        X_OFFSET=-1;
        Y_OFFSET=0;
        break;
      }
      for(std::vector<int>::const_iterator it = tilesInformation.begin();it!=tilesInformation.end();it+=5) {
        if(*it==currentX+X_OFFSET && *(it+1)==currentY+Y_OFFSET) {
          if(*(it+2)==tileN) {
            if(*(it+3)==tileO) {
              superBREAK=true; 
              break;
            }
            else {
              for(int test_face=0; test_face<4;++test_face) {
                if(genome[*(it+2)*4+(*(it+3)+test_face)%4]!=genome[tileN*4+(tileO+test_face)%4])
                  return {-1};
              }
              superBREAK=true; 
              break;
            }
          }            
          else
            return {-1};
        }
      }
      if(superBREAK)
        continue;
      tilesInformation.insert(tilesInformation.end(),{currentX+X_OFFSET,currentY+Y_OFFSET,tileN,tileO,f});
      processingQueue.insert(processingQueue.end(),{currentX+X_OFFSET,currentY+Y_OFFSET,tileN,tileO,f});
      if(tilesInformation.size()>1.25*genome.size()*genome.size()) {
#pragma omp critical(ster)
        {
          std::cout<<"HARD EXIT, STERIC OVERFLOW from"<<std::endl;
          for(auto& g: genome) {
            std::cout<<g<<" ";
          }
          std::cout<<std::endl;
        }
        return {-1};        
      }
    }    
  }

  std::vector<int> X_LOCs, Y_LOCs;               
  for(auto it = tilesInformation.begin();it!=tilesInformation.end();it+=5) {
    X_LOCs.emplace_back(*it);
    Y_LOCs.emplace_back(*(it+1));
  }
  std::vector<int>::iterator LEFT_X,RIGHT_X,TOP_Y,BOTTOM_Y;
  
  std::tie(LEFT_X,RIGHT_X)=std::minmax_element(X_LOCs.begin(),X_LOCs.end());
  std::tie(BOTTOM_Y,TOP_Y)=std::minmax_element(Y_LOCs.begin(),Y_LOCs.end());
  
  DELTA_X=*RIGHT_X-*LEFT_X+1;
  DELTA_Y=*TOP_Y-*BOTTOM_Y+1;

  std::vector<int> Spatial_Occupation(DELTA_X*DELTA_Y);
  phenotype_Shape.assign(((DELTA_X+2)*(DELTA_Y+2)),0);
  std::vector<int> Spatial_Grid(DELTA_X*DELTA_Y);
  int UNIQUE_TILES=0;
  for(unsigned int tileIndex=0;tileIndex<X_LOCs.size();++tileIndex) {
    Spatial_Occupation[(*TOP_Y-Y_LOCs[tileIndex])*DELTA_X + (X_LOCs[tileIndex]-*LEFT_X)]=1;
    ++Spatial_Grid[(*TOP_Y-Y_LOCs[tileIndex])*DELTA_X + (X_LOCs[tileIndex]-*LEFT_X)];
    if(Spatial_Grid[(*TOP_Y-Y_LOCs[tileIndex])*DELTA_X + (X_LOCs[tileIndex]-*LEFT_X)]==1) {
      ++UNIQUE_TILES;
    }
  }
  for(int x_shift=0;x_shift<DELTA_X;++x_shift) {
    for(int y_shift=0;y_shift<DELTA_Y;++y_shift) {
      phenotype_Shape[(DELTA_X+2)+1 + x_shift+y_shift*(DELTA_X+2)]=Spatial_Occupation[x_shift+y_shift*DELTA_X];
    }
  }
  DELTA_X+=2;
  DELTA_Y+=2;
  phenotype_Shape.emplace_back(DELTA_X);
  phenotype_Shape.emplace_back(DELTA_Y);
  return phenotype_Shape;
  
}

//STERIC METHOD
//RETURNS UNIQUE SHAPE ID
int Steric_Check_Table(std::vector<int>& genome,std::vector<int>& Known_Shapes,int& Num_Shapes) {

  std::vector<int> tilesInformation{0,0,0,0,-1};
  std::vector<int> processingQueue{0,0,0,0,-10};
  int currentX=0;
  int currentY=0;
  int currentTile=0;
  int currentOrientation=0;
  int comingFrom=-1;
  int ip=0;
  bool superBREAK=false;
  bool advanced_Steric_Check=false;
  while(!processingQueue.empty()) {
    ++ip;
    comingFrom=processingQueue.back();
    processingQueue.pop_back();
    currentOrientation=processingQueue.back();
    processingQueue.pop_back();
    currentTile=processingQueue.back();
    processingQueue.pop_back();
    currentY=processingQueue.back();
    processingQueue.pop_back();
    currentX=processingQueue.back();
    processingQueue.pop_back();

    for(int f=0;f<4;++f) {
      superBREAK=false;
      if((f+2)%4==comingFrom) {
        continue;
      }
      int leadingEdge=genome[currentTile*4+(f-currentOrientation+4)%4];

      if(!leadingEdge) {
        continue;
      }
      int bondingEdge=Interaction_Matrix(leadingEdge);//((1-leadingEdge%2)*(leadingEdge-1)+(leadingEdge%2)*(leadingEdge+1));
      int foundIndex=std::find(genome.begin(),genome.end(),bondingEdge)-genome.begin();
      int tileN=foundIndex/4;
      int faceN=foundIndex%4; //should be 3 since f=1
      int tileO=((f+2)%4-faceN +4)%4;

      int X_OFFSET,Y_OFFSET;
      switch(f) {
      case 0:
        X_OFFSET=0;
        Y_OFFSET=1;
        break;
      case 1:
        X_OFFSET=1;
        Y_OFFSET=0;
        break;
      case 2:
        X_OFFSET=0;
        Y_OFFSET=-1;
        break;

      case 3:
        X_OFFSET=-1;
        Y_OFFSET=0;
        break;
      }
      for(auto it = tilesInformation.begin();it!=tilesInformation.end();it+=5) {
        if(*it==currentX+X_OFFSET && *(it+1)==currentY+Y_OFFSET) {
          if(*(it+2)==tileN) {
            if(*(it+3)==tileO) {
              superBREAK=true; 
              break;
            }
            else {
              //check for something
              for(int test_face=0; test_face<4;++test_face) {
                if(genome[*(it+2)*4+(*(it+3)+test_face)%4]!=genome[tileN*4+(tileO+test_face)%4]) {
                  return -1;
                }
              }
              superBREAK=true; 
              break;
            }
          }
          else {
            return -1;
            if(std::count(genome.begin()+*(it+2)*4,genome.begin()+(1+*(it+2))*4,0)==std::count(genome.begin()+tileN*4,genome.begin()+(1+tileN)*4,0)) {
              advanced_Steric_Check=true;
              break;
            }
            else {
              return -1;
            }
          }
        }
      }
      if(superBREAK) {
        continue;
      }
      tilesInformation.insert(tilesInformation.end(),{currentX+X_OFFSET,currentY+Y_OFFSET,tileN,tileO,f});
      processingQueue.insert(processingQueue.end(),{currentX+X_OFFSET,currentY+Y_OFFSET,tileN,tileO,f});
    }    
  }

  
 
  std::vector<int> X_LOCs, Y_LOCs;
                          
  for(auto it = tilesInformation.begin();it!=tilesInformation.end();it+=5) {
    X_LOCs.emplace_back(*it);
    Y_LOCs.emplace_back(*(it+1));
    //std::cout<<"X "<<*it<<", Y "<<*(it+1)<<std::endl;
  }
  std::vector<int>::iterator LEFT_X,RIGHT_X,TOP_Y,BOTTOM_Y;
  
  std::tie(LEFT_X,RIGHT_X)=std::minmax_element(X_LOCs.begin(),X_LOCs.end());
  std::tie(BOTTOM_Y,TOP_Y)=std::minmax_element(Y_LOCs.begin(),Y_LOCs.end());
  
  int DELTA_X=*RIGHT_X-*LEFT_X+1;
  int DELTA_Y=*TOP_Y-*BOTTOM_Y+1;

  std::vector<int> Spatial_Occupation(DELTA_X*DELTA_Y);
  std::vector<int> Spatial_Grid(DELTA_X*DELTA_Y);
 
  for(unsigned int tileIndex=0;tileIndex<X_LOCs.size();++tileIndex) {
    Spatial_Occupation[(*TOP_Y-Y_LOCs[tileIndex])*DELTA_X + (X_LOCs[tileIndex]-*LEFT_X)]=1;
    ++Spatial_Grid[(*TOP_Y-Y_LOCs[tileIndex])*DELTA_X + (X_LOCs[tileIndex]-*LEFT_X)];
  }
  
  if(advanced_Steric_Check) {
    bool doubleBreak=false;
    for(int r=0;r<DELTA_Y;++r) {
      for(int c=0;c<DELTA_X;++c) {
        if(Spatial_Grid[r*DELTA_X+c]==1) {
          if(!Traverse_Numbered_Sterics(Spatial_Grid,DELTA_X,DELTA_Y,c,r,-1)) {
            return -1;
          }
          else {
            doubleBreak=true;
            break;
          }
        }
      }
      if(doubleBreak) {
        break;
      }
    }
  }
  
  if(Known_Shapes.empty()) {
    Known_Shapes.emplace_back(DELTA_X);
    Known_Shapes.emplace_back(DELTA_Y);   
    Known_Shapes.insert(Known_Shapes.end(),Spatial_Occupation.begin(),Spatial_Occupation.end());
    return 0;
  }

  
  
  int Phenotype_Num=0;
  for(auto it = Known_Shapes.begin();it!=Known_Shapes.end();) {
    if(*it==DELTA_X && *(it+1)==DELTA_Y && std::accumulate(it+2,it+2+*it* *(it+1),0)==std::accumulate(Spatial_Occupation.begin(),Spatial_Occupation.end(),0))  {
      if(std::equal(it+2,it+2+DELTA_X*DELTA_Y,Spatial_Occupation.begin())) {
        return Phenotype_Num;
      }
      Clockwise_Pi_Rotation(Spatial_Occupation,DELTA_X,DELTA_Y);
      if(std::equal(it+2,it+2+DELTA_X*DELTA_Y,Spatial_Occupation.begin())) {
        //std::cout<<"everywhere"<<std::endl;
        return Phenotype_Num;
      }

    }  
    if(*it==DELTA_Y && *(it+1)==DELTA_X && std::accumulate(it+2,it+2+*it* *(it+1),0)==std::accumulate(Spatial_Occupation.begin(),Spatial_Occupation.end(),0)) { //rotated version
      Clockwise_Rotation(Spatial_Occupation,DELTA_X,DELTA_Y);
      std::swap(DELTA_X,DELTA_Y);
      if(std::equal(it+2,it+2+DELTA_X*DELTA_Y,Spatial_Occupation.begin())) {
        //std::cout<<"here"<<std::endl;
        return Phenotype_Num;
      }
      Clockwise_Pi_Rotation(Spatial_Occupation,DELTA_X,DELTA_Y);
      if(std::equal(it+2,it+2+DELTA_X*DELTA_Y,Spatial_Occupation.begin())) {
        //std::cout<<"there"<<std::endl;
        return Phenotype_Num;
      }
    } 
    it+=*it* *(it+1)+2;
    ++Phenotype_Num;
  }
  //std::cout<<"For shape "<<Num_Shapes+1<<std::endl;
  //for(auto& g:genome) {
  //  std::cout<<g<<" ";
  //}
  //std::cout<<std::endl;
  Known_Shapes.emplace_back(DELTA_X);
  Known_Shapes.emplace_back(DELTA_Y);
  Known_Shapes.insert(Known_Shapes.end(),Spatial_Occupation.begin(),Spatial_Occupation.end());
  ++Num_Shapes;
  return Num_Shapes;
}


bool Traverse_Numbered_Sterics(std::vector<int>& Spatial_Grid,int DELTA_X,int DELTA_Y,int x,int y,int cameFrom) {
  for(int f=0;f<4;++f) {
    if(cameFrom!=-1 && f==(cameFrom+2)%4) {
      continue;
    }
    int nextX,nextY;
    switch(f) {
    case 0:
      nextX=x;
      nextY=y-1;
      break;
    case 1:
      nextX=x+1;
      nextY=y;
      break;
    case 2:
      nextX=x;
      nextY=y+1;
      break;
    case 3:
      nextX=x-1;
      nextY=y;
      break;
    }
    if((nextX<0||nextX>=DELTA_X)||(nextY<0||nextY>=DELTA_Y)) {
      continue;
    }
    if(Spatial_Grid[nextX+nextY*DELTA_X]>=2) {
      continue;
    }
    if(Spatial_Grid[nextX+nextY*DELTA_X]==1) {
      Spatial_Grid[x+y*DELTA_X]=0;
      Traverse_Numbered_Sterics(Spatial_Grid,DELTA_X,DELTA_Y,nextX,nextY,f);
    }  
  }
  if(cameFrom!=-1) {
    Spatial_Grid[x+y*DELTA_X]=0;
    return false;
  }
  if(cameFrom==-1) {
    if(std::count(Spatial_Grid.begin(),Spatial_Grid.end(),1)>=1) {
      return false;
    }
    else {
      return true;
    }
  }
  return false;
}
  
void Clockwise_Rotation(std::vector<int>& Spatial_Occupation,int DELTA_X,int DELTA_Y) {
  std::vector<int> swapper;
  for(int column=0;column<DELTA_X;++column) {
    for(int row=DELTA_Y-1;row>=0;--row) {
      swapper.emplace_back(Spatial_Occupation[row*DELTA_X+column]);
    }
  }
  Spatial_Occupation=swapper;
}

void Clockwise_Pi_Rotation(std::vector<int>& Spatial_Occupation,int DELTA_X,int DELTA_Y) {
  std::vector<int> swapper;
  for(int row=DELTA_Y-1;row>=0;--row) {
    for(int column=DELTA_X-1;column>=0;--column) {
      swapper.emplace_back(Spatial_Occupation[row*DELTA_X+column]);
    }
  }
  Spatial_Occupation=swapper;
}

bool GetMultiplePhenotypeFitness(std::vector<int> genome,std::vector<int> target_types,std::vector<double>& target_fitnesses,int active_targets) {
  std::vector<int> phenotype_information;
  Clean_Genome(genome,-1);
  if(Disjointed_Check(genome))
    genome.erase(std::find(genome.begin(),genome.end(),-1),genome.end());
  if(Graph_Analysis(genome)>0)
      phenotype_information=Steric_Check_Occupation(genome);
  else
    return false;
  if(phenotype_information.size()==1)
    return false;
  std::vector<int> phenotype_shape(phenotype_information.begin(),phenotype_information.end()-2);
  
  int dx=*(phenotype_information.end()-2);
  int dy=*(phenotype_information.end()-1);
  bool early_exit=false;
  for(int nth=0;nth<active_targets;++nth) {
    target_fitnesses[nth]=Shape_Matching_Fitness_Function(phenotype_shape,dx,dy,target_types[nth]);
    if(target_fitnesses[nth]<1.)
      early_exit=true;
  }

  for(unsigned int nth=active_targets;nth<target_types.size();++nth) {
    //if(early_exit) {
      //target_fitnesses[nth]=-1;
    //}
    //else {
      target_fitnesses[nth]=Shape_Matching_Fitness_Function(phenotype_shape,dx,dy,target_types[nth]);
      if(target_fitnesses[nth]<1.)
        early_exit=true;
      //} 
  }
  return true;
}


double Get_Phenotype_Fitness(std::vector<int> genome,int Shape_Based_Fitness, bool zeroth_seed) {
  Clean_Genome(genome,-1);
  if(Disjointed_Check(genome)) {
    std::vector<double> disjointed_results;
    while(std::find(genome.begin(),genome.end(),-1)!=genome.end()) {
      std::vector<int>::iterator foundAt=std::find(genome.begin(),genome.end(),-1);
      std::vector<int> partialGenome(genome.begin(),foundAt);
      genome.erase(genome.begin(),foundAt+1);
      if(Graph_Analysis(partialGenome)>0) {
        double steric_result=Steric_Check(partialGenome,Shape_Based_Fitness);
        if(steric_result>0) {
          if(zeroth_seed)
            return steric_result;
          else
            disjointed_results.emplace_back(steric_result);
        }
        else 
          return 0.0;
      }
      else
        return 0.0;
    }
    return *std::max_element(disjointed_results.begin(),disjointed_results.end());
  }
  else //Single connected component
    return Graph_Analysis(genome)>0 ? std::max(0.0,Steric_Check(genome,Shape_Based_Fitness)): 0.0;
}

    

double Shape_Matching_Fitness_Function(std::vector<int>& Spatial_Occupation,int& Delta_X,int& Delta_Y,int target_choice) {

  //HARDCODED PARAMETERS
  bool Size_Penalty_On=true;
  int asymmetry_factor=2;

  //TARGETS
  std::vector<int> Target_Occupation;
  int Target_Delta_X, Target_Delta_Y;
  int symmetric_rotations;
  /*Codes: 0 is empty, -1 is ignored, 1 is filled*/
  switch(target_choice) {
  case 1: //short T shape
    Target_Occupation={-1,0,-1,-1, 0,1,0,-1, 0,1,1,1, 0,1,0,-1, -1,0,-1,-1};
    Target_Delta_X=4;
    Target_Delta_Y=5;
    symmetric_rotations=4;
    break;
  case 2: //short L shape
    Target_Occupation={-1,0,-1,-1, 0,1,0,-1, 0,1,1,1, -1,0,0,-1};
    Target_Delta_X=4;
    Target_Delta_Y=4;
    symmetric_rotations=4;
    break;
  case 3: //long L shape
    Target_Occupation={-1,1,-1,-1, 0,1,0,-1, 0,1,0,-1, 0,1,1,0, -1,0,0,-1};
    Target_Delta_X=4;
    Target_Delta_Y=5;
    symmetric_rotations=4;
    break;
    
  case 4: //small square
    Target_Occupation={-1,1,0,-1, 0,1,1,-1, -1,1,1,0, -1,0,-1,-1};
    Target_Delta_X=4;
    Target_Delta_Y=4;
    symmetric_rotations=4;
    break;
  case 5: //long T shape
    Target_Occupation={-1,-1,-1,0,-1,-1, -1,0,0,1,0,-1, 0,1,1,1,1,1, -1,0,0,1,0,-1, -1,-1,-1,0,-1,-1};
    Target_Delta_X=6;
    Target_Delta_Y=5;
    symmetric_rotations=4;
    break;
  case 6: //tetris shape
    Target_Occupation={-1,0,0,-1,-1, 0,1,1,0,-1, -1,0,1,1,1, -1,-1,0,0,-1};
    Target_Delta_X=5;
    Target_Delta_Y=4;
    symmetric_rotations=4;
    break;
  case 7: //space lander shape
    Target_Occupation={-1,0,0,1,0,0,-1, 0,1,1,1,1,1,0, 0,1,0,0,0,1,0, -1,0,-1,-1,-1,0,-1};
    Target_Delta_X=7;
    Target_Delta_Y=4;
    symmetric_rotations=4;
    break;
  case 8: //stair shape
    Target_Occupation={-1,0,-1,-1,-1, 0,1,0,-1,-1, 0,1,1,0,-1, -1,0,1,1,0, -1,-1,0,1,-1};
    Target_Delta_X=5;
    Target_Delta_Y=5;
    symmetric_rotations=4;
    break;
  case 0: //small cross 
  default: 
    Target_Occupation={-1,-1,0,-1,-1, -1,0,1,0,-1, 0,1,1,1,1, -1,0,1,0,-1, -1,-1,0,-1,-1};
    Target_Delta_X=5;
    Target_Delta_Y=5;
    symmetric_rotations=1;
  }
  
  int max_possible_overlap=Target_Delta_X*Target_Delta_Y-std::count(Target_Occupation.begin(),Target_Occupation.end(),-1);
  std::vector<double> Fitness_Overlaps;


  for(int rotations=0;rotations<symmetric_rotations;++rotations) {
    if(Delta_X>=Target_Delta_X && Delta_Y>=Target_Delta_Y) {
      Fitness_Overlaps.reserve((Delta_X-Target_Delta_X+1)*(Delta_Y-Target_Delta_Y+1)*symmetric_rotations);
      for(int Y_slides=0;Y_slides<=Delta_Y-Target_Delta_Y; ++Y_slides) {
        for(int X_slides=0;X_slides<=Delta_X-Target_Delta_X; ++X_slides) {
          int Positive_Match=0;
          for(int Y_Check=0; Y_Check<Target_Delta_Y; ++Y_Check) {
            for(int X_Check=0; X_Check<Target_Delta_X; ++X_Check) {                       
              if(Spatial_Occupation[X_Check+Y_Check*Delta_X+X_slides+Y_slides*Delta_X]==Target_Occupation[X_Check+Y_Check*Target_Delta_X]) {
                ++Positive_Match;
              }
              else {
                if(Spatial_Occupation[X_Check+Y_Check*Delta_X+X_slides+Y_slides*Delta_X]==1 && Target_Occupation[X_Check+Y_Check*Target_Delta_X]==0) {
                  Positive_Match-=asymmetry_factor;
                }
              }
            }
          }
          if(Positive_Match==max_possible_overlap)
            return 1.;
          Fitness_Overlaps.emplace_back(static_cast<double>(Positive_Match)/max_possible_overlap);
        }
      }   
    }
    else {
      if(Delta_X<=Target_Delta_X && Delta_Y<=Target_Delta_Y) {
        Fitness_Overlaps.reserve((Target_Delta_X-Delta_X+1)*(Target_Delta_Y-Delta_Y+1)*symmetric_rotations);
        for(int Y_slides=0;Y_slides<=Target_Delta_Y-Delta_Y; ++Y_slides) {
          for(int X_slides=0;X_slides<=Target_Delta_X-Delta_X; ++X_slides) {
            int Positive_Match=0;
            for(int Y_Check=0; Y_Check<Target_Delta_Y; ++Y_Check) {
              for(int X_Check=0; X_Check<Target_Delta_X; ++X_Check) {
                if(X_Check>=(Delta_X+X_slides) || Y_Check>=(Delta_Y+Y_slides) || X_Check<X_slides || Y_Check<Y_slides) {
                  if(Target_Occupation[X_Check+Y_Check*Target_Delta_X]==0) {
                    ++Positive_Match;
                  }
                }
                else {
                  if(Spatial_Occupation[X_Check+Y_Check*Delta_X-X_slides-Y_slides*Delta_X]==Target_Occupation[X_Check+Y_Check*Target_Delta_X]) {
                    ++Positive_Match;
                  }
                  else {
                    if(Spatial_Occupation[X_Check+Y_Check*Delta_X+X_slides+Y_slides*Delta_X]==1 && Target_Occupation[X_Check+Y_Check*Target_Delta_X]==0) {
                      Positive_Match-=asymmetry_factor;
                    }
                  }
                }
              }
            }
            if(Positive_Match==max_possible_overlap)
              return 1.;
            Fitness_Overlaps.emplace_back(static_cast<double>(Positive_Match)/max_possible_overlap);
          }
        }
      }
      else { //Special non-overlapping cases
        if(Delta_X>=Target_Delta_X && Delta_Y<Target_Delta_Y) {
          Clockwise_Rotation(Spatial_Occupation,Delta_X,Delta_Y);
          std::swap(Delta_X,Delta_Y);
        }
        if(Delta_X<Target_Delta_X && Delta_Y>=Target_Delta_Y) {
          Fitness_Overlaps.reserve((Target_Delta_X-Delta_X+1)*(Delta_Y-Target_Delta_Y+1)*symmetric_rotations);
          for(int Y_slides=0;Y_slides<=Delta_Y-Target_Delta_Y; ++Y_slides) {
            for(int X_slides=0;X_slides<=Target_Delta_X-Delta_X; ++X_slides) {
              int Positive_Match=0;
              for(int Y_Check=0; Y_Check<Target_Delta_Y; ++Y_Check) {
                for(int X_Check=0; X_Check<Target_Delta_X; ++X_Check) {
                  if(X_Check>=(Delta_X+X_slides) || X_Check<X_slides) {
                    if(Target_Occupation[X_Check+Y_Check*Target_Delta_X]==0) {
                      ++Positive_Match;
                    }
                  }
                  else {
                    if(Spatial_Occupation[X_Check+Y_Check*Delta_X+Y_slides*Delta_X-X_slides]==Target_Occupation[X_Check+Y_Check*Target_Delta_X]) {
                      ++Positive_Match;
                    }
                    else {
                    if(Spatial_Occupation[X_Check+Y_Check*Delta_X+X_slides+Y_slides*Delta_X]==1 && Target_Occupation[X_Check+Y_Check*Target_Delta_X]==0) {
                      Positive_Match-=asymmetry_factor;
                    }
                  }
                    
                  }
                }
              }
              if(Positive_Match==max_possible_overlap)
                return 1.;
              Fitness_Overlaps.emplace_back(static_cast<double>(Positive_Match)/max_possible_overlap);
            }
          }
        }
      }
    }
    Clockwise_Rotation(Target_Occupation,Target_Delta_X,Target_Delta_Y);
    std::swap(Target_Delta_X,Target_Delta_Y);
  }
  if(Size_Penalty_On) {
    return std::max(0.,*std::max_element(Fitness_Overlaps.begin(), Fitness_Overlaps.end()))*std::min(1.,std::accumulate(Spatial_Occupation.begin(),Spatial_Occupation.end(),0.0)/std::count(Target_Occupation.begin(),Target_Occupation.end(),1));
  }
  else
    return *std::max_element(Fitness_Overlaps.begin(), Fitness_Overlaps.end());
  
}

/*
int SpiralCoordinate(int x, int y) {
  int diag_coord;
  if(abs(x) > abs(y))
    diag_coord=x;
  else {
    if(abs(y) > abs(x))
      diag_coord=y;
    else
      diag_coord= x>0? std::min(x,y) : std::max(x,y);
  }
  int diag_offset=4*diag_coord*diag_coord-2*diag_coord;
  if(x+y>0)
    return diag_offset+(y-x);
  else {
    if(x+y<0)
      return diag_offset-(y-x);
    else {
      return x>y? diag_offset+(x-y) : diag_offset+(y-x); 
    }      
  }    
}
double StericNew(std::vector<int>& genome) {
  std::vector<int> tilesInformation{0,0,0,0};
  std::vector<int> processing_queue{0,0,0};
  int new_x=0, new_y=0, new_interaction=0,
  //initial check
  int X_OFFSET,Y_OFFSET;
  for(int face=0; face<4;++face) {
    if(genome[face]) {
      switch(face) {
      case 0:
        X_OFFSET=0;Y_OFFSET=1;
        break;
      case 1:
        X_OFFSET=1;Y_OFFSET=0;
        break;
      case 2:
        X_OFFSET=0;Y_OFFSET=-1;
        break;
      case 3:
        X_OFFSET=-1;Y_OFFSET=0;
        break;
      }
      processing_queue.insert(processing_queue.end(),{X_OFFSET,Y_OFFSET,Interaction_Matrix(genome[face])})l
    }
    while(!proccessing_queue.empty()) {
      new_interaction=processing_queue.back();processing_queue.pop_back();
      new_y=processing_queue.back();processing_queue.pop_back();
      new_x=processing_queue.back();processing_queue.pop_back();
    }
  }
}
*/

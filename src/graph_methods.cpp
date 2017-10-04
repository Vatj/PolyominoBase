#include "graph_methods.hpp"


//Provides the conjugate face
//Positive interfaces: odd<-->even
//Negative interfaces: self-interaction 
int Interaction_Matrix(int input_face) {
  return input_face>0 ?  (1-input_face%2)*(input_face-1)+(input_face%2)*(input_face+1) : input_face;
}

int SIF_Elimination(std::vector<int>& genome,bool erase_ends) {
  if(genome.size()==(unsigned int)std::count(genome.begin(),genome.end(),0) || two_SIF_Check(genome))
    return 5;
  for(std::vector<int>::iterator genome_iter=genome.begin();genome_iter!=genome.end();genome_iter+=4) {
    if(std::count(genome_iter,genome_iter+4,0)==3) { //IS A SIF
      std::vector<int>::iterator non_zero_iter=std::find_if(genome_iter,genome_iter+4,[](int x) { return x != 0; });
      if(std::count(genome.begin(),genome.end(),*non_zero_iter)>=2)
        return -3; //can reject as dangling end
      int conjugate_face=Interaction_Matrix(*(non_zero_iter));
      std::vector<int>::iterator conjugate_iter=std::find(genome.begin(),genome.end(),conjugate_face);
      *(non_zero_iter)=0;
      int zeroed_conjugates=0;
      while(conjugate_iter!=genome.end()) {
        *(conjugate_iter)=0;
        ++zeroed_conjugates;
        conjugate_iter=std::find(conjugate_iter+1,genome.end(),conjugate_face);
      }
      if(erase_ends) 
        genome.erase(genome_iter,genome_iter+4);
      if(zeroed_conjugates>1 && erase_ends)
        return Disjointed_Check(genome) ? -3 : SIF_Elimination(genome,erase_ends);
      else
        return SIF_Elimination(genome,erase_ends);
    }
  }
  return 0;
}     

bool two_SIF_Check(std::vector<int>& genome) {
  if(genome.size()==8 && std::count(genome.begin(),genome.end(),0)==6)
    return *std::max_element(genome.begin(),genome.begin()+4)==Interaction_Matrix(*std::max_element(genome.begin()+4,genome.end()));
  else
    return false;
}

bool BP_Disjointed_Check(std::vector<int> genome) {
  if(genome.size()<=8)
    return false;
  for(int i=1;i<=*std::max_element(genome.begin(),genome.end());i+=2) {
    std::vector<int> Spare_genome(genome);
    std::set<int> BP_IDs;
    int lowerCount=std::count(Spare_genome.begin(),Spare_genome.end(),i);
    int BP_out=0,BP_in=0;
    if(lowerCount==1) {
      if(std::count(Spare_genome.begin(),Spare_genome.end(),i+1)>1) {
        BP_out=i;
        BP_in=i+1;      
      }
    }
    else {
      if(lowerCount>1) {
        BP_out=i+1;
        BP_in=i;
      }
    }
    if(BP_out==0)
      continue;
    std::replace(Spare_genome.begin(),Spare_genome.end(),BP_out,0);
    std::vector<int>::iterator face_iter=std::find(Spare_genome.begin(),Spare_genome.end(),BP_in);
    while(face_iter!=Spare_genome.end()) {
      BP_IDs.insert((face_iter-Spare_genome.begin())/4);
      *face_iter=0;
      face_iter=std::find(face_iter+1,Spare_genome.end(),BP_in);
    }
    
    if(BP_IDs.size()>1) {
      std::vector<int> BP_genome;
      for(int tile :BP_IDs)
        BP_genome.insert(BP_genome.end(),Spare_genome.begin()+tile*4,Spare_genome.begin()+tile*4+4);
      if(Disjointed_Check(BP_genome))
        return true;
    }
    else
      continue; //not disjointed, no point checking   
  }
  return false;
}


bool BP_Check(std::vector<int>& genome) {
  //Handle positive "asymmetric iteractions"
  for(int i=1;i<=*std::max_element(genome.begin(),genome.end());i+=2) {
    int lowerCount=std::count(genome.begin(),genome.end(),i);
    if(lowerCount>1)
      return true;
    if(lowerCount==1)
      if(std::count(genome.begin(),genome.end(),i+1)>1) 
        return true;   
  }
  return false;
}

bool Double_BP_Check(std::vector<int>& genome) {
  //Handle negative "self-interactions"
  for(int j=*std::min_element(genome.begin(),genome.end());j<0;++j) {
    if(std::count(genome.begin(),genome.end(),j)>=2)
      return true;
  }
  //Handle positive "asymmetric iteractions"
  for(int i=1;i<=*std::max_element(genome.begin(),genome.end());i+=2) {
    int lowerCount=std::count(genome.begin(),genome.end(),i);
    if(lowerCount>=2)
      if(std::count(genome.begin(),genome.end(),i+1)>=2)
        return true;
  }
  return false;
}

void Trim_Topology(std::vector<int> genome,std::vector<int>& looping_genome) {
  std::vector<int> replacements;
  for(std::vector<int>::iterator face_iter=genome.begin();face_iter!=genome.end();face_iter+=4) {
    int zero_count_offset=std::count(face_iter,face_iter+4,0)>0? 0:1;
    std::vector<int> tile(face_iter,face_iter+4);
    std::sort(face_iter,face_iter+4);
    if((std::unique(face_iter, face_iter+4)-face_iter)<(3-zero_count_offset)) {
      replacements.emplace_back(Interaction_Matrix(*std::find_if(face_iter,face_iter+4,[](int x) { return x != 0; })));
      looping_genome.insert(looping_genome.end(),{0,0,0,0});
    }
    else
      looping_genome.insert(looping_genome.end(),tile.begin(),tile.end());
  }
  for(std::vector<int>::iterator replace_iter=replacements.begin();replace_iter!=replacements.end();replace_iter++) {
    int dead_replace=*replace_iter;
    std::replace(looping_genome.begin(),looping_genome.end(),dead_replace,0);
  }
}

int loop_Analysis(std::vector<int>& genome, int& Bound_Limit) {
  std::vector<int> looping_genome;
  Trim_Topology(genome,looping_genome);
  int rank1=0, rank2=0, rank4=0;

  int NUM_TILES=looping_genome.size()/4;
  for(std::vector<int>::iterator tile_iter=looping_genome.begin();tile_iter!=looping_genome.end();tile_iter+=4) {
    for(int face=0;face<4;++face) {
      if(*(tile_iter+face)==0)
        continue;
      int conjugate_Face=Interaction_Matrix(*(tile_iter+face));
      if(conjugate_Face>0) {
        int found_Internal_Index=std::find(tile_iter+face+1,tile_iter+4,conjugate_Face)-tile_iter;
        while(found_Internal_Index!=4) {
          if(found_Internal_Index==face+2) //infite loop
            return -5;
          if((std::count(looping_genome.begin(),looping_genome.end(),*(tile_iter+face))+std::count(looping_genome.begin(),looping_genome.end(),*(tile_iter+found_Internal_Index)))>2) //internal BP
            return -6;
          else {
            *(tile_iter+face)=0;
            *(tile_iter+found_Internal_Index)=0;
          }
          ++rank4;
          Bound_Limit+=3;
          found_Internal_Index=std::find(tile_iter+found_Internal_Index+1,tile_iter+4,conjugate_Face)-tile_iter;
        }
      }
      else {
        ++rank2;
        Bound_Limit+=1;
        *(tile_iter+face)=0;
      }
    }
  }

  

  if((rank4+rank2)>=2) //Can reject as multiple internal H.O. loops
    return -7;

  if(NUM_TILES==1)
    return BP_Check(genome) ? -8 : 1;

  if(rank2+rank4>0)
    SIF_Elimination(looping_genome,false);

  bool found_harmful_loop=false;
  
  std::vector<int> loop_Path;
  std::vector<int> loop_Histories;

  for(int tile=0;tile<NUM_TILES;++tile) {
    for(int face=0;face<4;++face) {
      if(looping_genome[tile*4+face]==0) 
        continue;
      
      int temp_rank2=0,temp_rank4=0,temp_rankInf=0;
      std::vector<int> temp_loop_Histories;
      //std::cout<<"stepping for "<<tile<<","<<face<<"which is "<<looping_genome[tile*4+face]<<std::endl;
      std::tie(temp_rank2,temp_rank4,temp_rankInf)=Take_Loop_Step(looping_genome,tile,face,0,{},temp_loop_Histories,false);
      //std::cout<<"stepped for "<<tile<<","<<face<<std::endl;
      //std::cout<<"Found "<<temp_rank2<<","<<temp_rank4<<","<<temp_rankInf<<std::endl;
      if(temp_rank2>=999)
        found_harmful_loop=true;
      bool simplified_R1=false;
      if(temp_rankInf>0) {
        if(NUM_TILES<4) 
          return -9; //Can reject as found true rank infinity loop
        else {
          for(std::vector<int>::iterator it = temp_loop_Histories.begin();it!=temp_loop_Histories.end(); ) {
            std::vector<int>::iterator endIter=std::find(it,temp_loop_Histories.end()-1, -1);
            if(*it==0) { //INFINITE LOOP
              std::vector<int> partial_History(it+1,endIter);
              if(!Loop_Rank_One_Check(partial_History))
                return -9; //Reject as infinite loop
              else {
                //--rankInf;
                ++rank1;
                if(!Cut_Rank_One_Loop(genome,partial_History))
                  return -10;


                SIF_Elimination(genome,false);


                Cut_Rank_One_Loop(looping_genome,partial_History);
 
                SIF_Elimination(looping_genome,false);

                simplified_R1=true;
                //std::cout<<"breaking out"<<std::endl;
                break;
                  
              }                  
            }
            it=endIter+1;
          }
        }
      }

      if(temp_rankInf==0 && ((temp_rank2>1 && temp_rank4==0) || (temp_rank4>1 && temp_rank2==0))) {
        //std::cout<<"skip it "<<std::endl;
        continue;
      }

      //Zero the current step
      int replacing_main=looping_genome[tile*4+face];
      if(replacing_main==0)
        continue;
      std::replace(looping_genome.begin(),looping_genome.end(),replacing_main,0);
      std::replace(looping_genome.begin(),looping_genome.end(),Interaction_Matrix(replacing_main),0);

      if(!simplified_R1 || temp_rankInf==0) {
        rank2+=temp_rank2;
        rank4+=temp_rank4;
        loop_Histories.insert(loop_Histories.end(),temp_loop_Histories.begin(),temp_loop_Histories.end());
      }
      
    }
  } //END i iterations
  if(found_harmful_loop) {
    if((genome.size()-std::count(genome.begin(),genome.end(),0))<=2)
      return 1;
    else
      return -11;
  }
  
  
  //std::cout<<"RANKS "<<rank1<<","<<rank2<<","<<rank4<<","<<std::endl;
  if(rank1==0 && (rank2+rank4)<=1) {
    if((rank2+rank4)==0)
      return -12;
    if(!loop_Histories.empty()) 
      Bound_Limit+= (*loop_Histories.begin()-1)*(loop_Histories.size()-2)/4;
    return 1;
  }
  if(rank1>0 && (rank2+rank4)==0)
    return 1;




  for(std::vector<int>::iterator HO_iterator=loop_Histories.begin(); HO_iterator!=loop_Histories.end();) {

    bool Duplicate_Loop=false;
    std::vector<int>::iterator end_iter=std::find(HO_iterator,loop_Histories.end(), -1);
    if(end_iter==loop_Histories.end())
      break;
    for(std::vector<int>::iterator inside_iter=HO_iterator+1;inside_iter!=end_iter;inside_iter+=2) {
      if(genome[*(inside_iter)*4 + *(inside_iter+1)]==0) {
        Duplicate_Loop=true;
        break;
      }
    }
    if(!Duplicate_Loop) 
      Bound_Limit+=(*HO_iterator-1)*(end_iter-HO_iterator)/4;
    HO_iterator=end_iter+1;
  }
  Trim_Zero_Tiles(genome);

  //std::cout<<"final"<<std::endl;
  //for(int g: genome)
  //  std::cout<<g<<" ";
  //std::cout<<std::endl;
  //std::cout<<"B "<<Bound_Limit+(Bound_Limit >0 ? genome.size()/4 : genome.size()/2)<<std::endl;
  return 0;
}



bool Cut_Rank_One_Loop(std::vector<int>& genome, std::vector<int>& loopR1_Path) {
  for(std::vector<int>::iterator cutting_iter = loopR1_Path.begin(); cutting_iter!=loopR1_Path.end(); cutting_iter+=4 ) {
    if(genome[*(cutting_iter)*4 + *(cutting_iter+1)]==0)
      continue;
    else {
      int replacing_main=genome[*(cutting_iter)*4 + *(cutting_iter+1)];
      std::replace(genome.begin(),genome.end(),replacing_main,0);
      std::replace(genome.begin(),genome.end(),Interaction_Matrix(replacing_main),0);
      return true;
    }
  }
  return false;
}

void Trim_Zero_Tiles(std::vector<int>& genome) {
  for(std::vector<int>::iterator tile_iter=genome.begin(); tile_iter!=genome.end();) {
    if(std::count(tile_iter,tile_iter+4,0)==4)
      tile_iter=genome.erase(tile_iter,tile_iter+4);
    else
      tile_iter+=4;
  }
}


bool Check_If_Loop(int tile,int face, std::vector<int>& loop_Path) {
  for(auto it = loop_Path.begin();it!=loop_Path.end();it+=2) {
    if(*it==tile)
      return true;        
  }
  return false;
}

std::tuple<int,int,int> Take_Loop_Step(std::vector<int>& genome,int tile,int face,int deltaFace,std::vector<int> loop_Path,std::vector<int>& loop_Histories,bool Passed_BP) {
  std::vector<int>::iterator toIter=std::find(genome.begin(),genome.end(),Interaction_Matrix(genome[tile*4+face]));
  if(toIter==genome.end())
    return std::make_tuple(0,0,0);
  else {
    int rank2=0, rank4=0, rankInf=0;
    int toTile=(toIter-genome.begin())/4, toFace=(toIter-genome.begin())%4;
    
    if(!loop_Path.empty() && (toTile==*(loop_Path.end()-4) && toFace==*(loop_Path.end()-3)) && (tile==*(loop_Path.end()-2) && face==*(loop_Path.end()-1)))
      return std::make_tuple(0,0,0);
    if(genome[toTile*4+toFace]==0 || std::count(genome.begin()+toTile*4,genome.begin()+toTile*4+4,0)==3)
      return std::make_tuple(0,0,0);

    int toCount=std::count(genome.begin(),genome.end(),genome[toTile*4+toFace]);
    if(toCount>1) {
      if(toCount>2) {
        return std::make_tuple(999,0,0);
      }
      else {
        if(std::count(genome.begin()+toTile*4,genome.begin()+toTile*4+4,genome[toTile*4+toFace])==2 && genome[toTile*4+toFace]==genome[toTile*4+(toFace+2)%4] && genome[toTile*4+(toFace+1)%4]==genome[toTile*4+(toFace+3)%4])
          Passed_BP=true;
        else
          return std::make_tuple(999,0,0);
      }
    }

    int fromCount=std::count(genome.begin(),genome.end(),genome[tile*4+face]);
    if(fromCount>1) {
      if(fromCount>2) {
        return std::make_tuple(999,0,0);
      }
      else {
        if(std::count(genome.begin()+tile*4,genome.begin()+tile*4+4,genome[tile*4+face])==2 && genome[tile*4+face]==genome[tile*4+(face+2)%4] && genome[tile*4+(face+1)%4]==genome[tile*4+(face+3)%4])
          Passed_BP=true;
        else
          return std::make_tuple(999,0,0);
      }
    }

      
    if(Check_If_Loop(toTile,toFace,loop_Path)) {        
      if(toTile!=*loop_Path.begin() || toFace==*(loop_Path.begin()+1)) 
        return std::make_tuple(rank2,rank4,rankInf); //'discovered loop 'too early', return upstream for now
        
      //Clear to treat as a proper loop
      int rankType=-1;
      switch(std::abs((((face+deltaFace)%4)+2)%4-toFace)) {
      case 0 : ++rankInf;
        rankType=0;
        break;
      case 2 : ++rank2;
        rankType=2;
        break;
      case 1 :
      case 3 :  ++rank4;
        rankType=4;
        break;
      }
      if(Passed_BP && rankType==4)
          return std::make_tuple(999,0,0);
        
      loop_Path.insert(loop_Path.end(),{tile,face,toTile,toFace});
      loop_Histories.emplace_back(rankType);
      loop_Histories.insert(loop_Histories.end(),loop_Path.begin(),loop_Path.end());
      loop_Histories.emplace_back(-1);
      return std::make_tuple(rank2,rank4,rankInf);
    }
      

    loop_Path.insert(loop_Path.end(),{tile,face,toTile,toFace});
    deltaFace+=((face+2)%4-toFace+4)%4;
    deltaFace%=4;      

    //bool foundNextMove=false;
    for(int m=1;m<4;++m) {
      if(genome[toTile*4+(toFace+m)%4]>0) {
        int temp_rank2=0,temp_rank4=0,temp_rankInf=0;
        std::tie(temp_rank2,temp_rank4,temp_rankInf)=Take_Loop_Step(genome,toTile,(toFace+m)%4,deltaFace,loop_Path,loop_Histories,Passed_BP);
        if(temp_rank2>=999)
          return std::make_tuple(temp_rank2,temp_rank4,temp_rankInf);        
        rank2+=temp_rank2;
        rank4+=temp_rank4;
        rankInf+=temp_rankInf;
      }
    }
    return std::make_tuple(rank2,rank4,rankInf);
  }
} 



bool Loop_Rank_One_Check(std::vector<int>& loop_Path) {
  int delta_X=0,delta_Y=0,delta_Face=0;
  for(std::vector<int>::const_iterator step_iter=loop_Path.begin(); step_iter!=loop_Path.end(); step_iter+=4) {
    switch((*(step_iter+1)+delta_Face+4)%4) {
    case 0: ++delta_Y;
      break;
    case 1: ++delta_X;
      break;
    case 2: --delta_Y;
      break;
    case 3: --delta_X; 
      break;
    }
    delta_Face= ((*(step_iter+1)+delta_Face+2 +4)%4 -*(step_iter+3)+4)%4; 
  }
  return (delta_X==0 && delta_Y==0);
}


void Clean_Genome(std::vector<int>& genome,int secondNonInteracting,bool Remove_Duplicates) {
  //defaults set to -1 and true//
  //Removes any non-bonding interfaces and reduces additional non-interacting faces to 0s
  //Note, no need to check for negatives, as either present and able to self-interact, or absent and no need to act
  if(secondNonInteracting!=-1) {
    std::replace(genome.begin(),genome.end(),secondNonInteracting,0);
  }
  for(int t=1;t<=*std::max_element(genome.begin(),genome.end());t+=2) { 
    if(std::count(genome.begin(),genome.end(),t)==0) //genotype doens't contain this face
        std::replace(genome.begin(),genome.end(),t+1,0);
    else //genotype does contain this face
      if(std::find(genome.begin(),genome.end(),t+1)==genome.end()) //genotype doesn't contain conjugate
        std::replace(genome.begin(),genome.end(),t,0);
  }
  Minimize_Tile_Set(genome);
  if(Remove_Duplicates) {
    for(int check_index=genome.size()/4-1;check_index>0;--check_index) {
      std::vector<int> check_Genome(genome.begin()+check_index*4,genome.begin()+check_index*4+4);
      for(int compare_index=check_index-1;compare_index>=0;--compare_index) {
        std::vector<int> compare_Genome(genome.begin()+compare_index*4,genome.begin()+compare_index*4+4);
        if(check_Genome==compare_Genome) {
          genome.erase(genome.begin()+check_index*4,genome.begin()+check_index*4+4);
          break;
        }
      }
    }
  }  
}



void Minimize_Tile_Set(std::vector<int>& genome) {
  std::vector<int> Minimal_Genome;
  for(std::vector<int>::iterator genome_iter=genome.begin();genome_iter!=genome.end();genome_iter+=4) {
    Minimal_Genome.assign(genome_iter,genome_iter+4);
    if(std::count(Minimal_Genome.begin(),Minimal_Genome.end(),0)==4)
      continue;
    for(int rotation=1;rotation<=3;++rotation) {
      std::rotate(genome_iter,genome_iter+1,genome_iter+4);
      for(int index=0;index<4;++index) {
        if(Minimal_Genome[index]<*(genome_iter+index))
          break;
        if(Minimal_Genome[index]==*(genome_iter+index))
          continue;
        else {
          Minimal_Genome.assign(genome_iter,genome_iter+4);
          break;
        }
      }
    }
    std::vector<int>::iterator Minimal_beginning=Minimal_Genome.begin();
    for(int swapping=0;swapping<4;++swapping) {
      *(genome_iter+swapping)=*(Minimal_beginning+swapping);
    }
  }
}


bool Disjointed_Check(std::vector<int>& genome) {
  std::vector<int> Unvisited(genome.size()/4-1);
  std::map<int, std::vector<int> > Connected_Components;
  Connected_Components[0].insert(Connected_Components[0].end(),genome.begin(),genome.begin()+4);
  std::iota(Unvisited.begin(),Unvisited.end(),1);
  Search_Next_Tile(genome,Unvisited,Connected_Components[0],0);
  while(!Unvisited.empty()) {
    int CC_Size=Connected_Components.size();
    int New_Launch_Point=Unvisited[0];
    Connected_Components[CC_Size].insert(Connected_Components[CC_Size].end(),genome.begin()+ New_Launch_Point*4,genome.begin()+New_Launch_Point*4+4);
    Unvisited.erase(Unvisited.begin());
    Search_Next_Tile(genome,Unvisited,Connected_Components[CC_Size], New_Launch_Point);
  }
  
  if(Connected_Components.size()>1) {
    genome.clear();
    for(auto& KV_pair : Connected_Components) {
      genome.insert(genome.end(),KV_pair.second.begin(),KV_pair.second.end());
      genome.emplace_back(-1);
    }
    return true;
  }
  return false;  
}

void Search_Next_Tile(std::vector<int>& genome, std::vector<int>& Unvisited, std::vector<int>& Connected_Components,int tile) {
  for(int face=0;face<4;++face) {
    if(genome[tile*4+face]==0)
      continue;
    int conjugate_Face=Interaction_Matrix(genome[tile*4+face]);
    std::vector<int>::iterator conjugate_iter=std::find(genome.begin(),genome.end(),conjugate_Face);
    while(conjugate_iter!=genome.end()) {
      int corresponding_Tile=(conjugate_iter-genome.begin())/4;
      if(std::find(Unvisited.begin(),Unvisited.end(),corresponding_Tile)!=Unvisited.end()) { //new tile
        Connected_Components.insert(Connected_Components.end(),genome.begin()+corresponding_Tile*4,genome.begin()+corresponding_Tile*4+4);
        Unvisited.erase(std::find(Unvisited.begin(),Unvisited.end(),corresponding_Tile));
        Search_Next_Tile(genome,Unvisited,Connected_Components,corresponding_Tile);
      } 
      conjugate_iter=std::find(conjugate_iter+1,genome.end(),conjugate_Face);
    }
  }
}



int FastLoopAnalysis(std::vector<int>& genome) {
  std::vector<int> looping_genome;
  Trim_Topology(genome,looping_genome);
  int rank1=0, rank2=0, rank4=0, rankInf=0;

  int NUM_TILES=looping_genome.size()/4;
  for(std::vector<int>::iterator tile_iter=looping_genome.begin();tile_iter!=looping_genome.end();tile_iter+=4) {
    for(int face=0;face<4;++face) {
      if(*(tile_iter+face)==0)
        continue;
      int conjugate_Face=Interaction_Matrix(*(tile_iter+face));
      if(conjugate_Face>0) {
        int found_Internal_Index=std::find(tile_iter+face+1,tile_iter+4,conjugate_Face)-tile_iter;
        while(found_Internal_Index!=4) {
          if(found_Internal_Index==face+2) //infite loop
            return -11;
          if((std::count(looping_genome.begin(),looping_genome.end(),*(tile_iter+face))+std::count(looping_genome.begin(),looping_genome.end(),*(tile_iter+found_Internal_Index)))>2) //internal BP
            return -12;
          else {
            *(tile_iter+face)=0;
            *(tile_iter+found_Internal_Index)=0;
          }
          ++rank4;
          found_Internal_Index=std::find(tile_iter+found_Internal_Index+1,tile_iter+4,conjugate_Face)-tile_iter;
        }
      }
      else {
        ++rank2;
        *(tile_iter+face)=0;
      }
    }
  }

  if((rank4+rank2)>=2) //Can reject as multiple internal H.O. loops
    return -10;

  if(NUM_TILES==1)
    return BP_Check(genome) ? -13 : 1;

  std::vector<int> loop_Path;
  std::vector<int> loop_Histories;

  for(int tile=0;tile<NUM_TILES;++tile) {
    for(int face=0;face<4;++face) {
      if(looping_genome[tile*4+face]==0) 
        continue;
      
      int temp_rank2=0,temp_rank4=0,temp_rankInf=0;
      std::vector<int> temp_loop_Histories;
      std::tie(temp_rank2,temp_rank4,temp_rankInf)=Take_Loop_Step(looping_genome,tile,face,0,{},temp_loop_Histories,false);
      if(temp_rank2>=999)
        return -1;
      else {
        rank2+=temp_rank2;
        rank4+=temp_rank4;
        rankInf+=temp_rankInf;
        if(temp_rankInf>0) {
          if(NUM_TILES<4) 
            return -111;
          else {
            for(std::vector<int>::iterator it = temp_loop_Histories.begin();it!=temp_loop_Histories.end(); ) {
              std::vector<int>::iterator endIter=std::find(it,temp_loop_Histories.end()-1, -1);
              if(*it==0) { //INFINITE LOOP
                std::vector<int> partial_History(it+1,endIter);
                if(Loop_Rank_One_Check(partial_History)) {
                  --rankInf;
                  ++rank1;
                  if(!Cut_Rank_One_Loop(genome,partial_History))
                    return -20;
                  Cut_Rank_One_Loop(looping_genome,partial_History);
                  SIF_Elimination(looping_genome,false);
                  SIF_Elimination(genome,false);
                }
                else
                  return -111; //Reject as infinite loop
              }
              else
                if(*it!=-1)
                  loop_Histories.insert(loop_Histories.end(),it,endIter+1);
              it=endIter+1;
            }
          }
        }
        else
          loop_Histories.insert(loop_Histories.end(),temp_loop_Histories.begin(),temp_loop_Histories.end());
        int replacing_main=looping_genome[tile*4+face];
        if(replacing_main==0)
          continue;
        std::replace(looping_genome.begin(),looping_genome.end(),replacing_main,0);
        std::replace(looping_genome.begin(),looping_genome.end(),Interaction_Matrix(replacing_main),0);
      }
      
    }
  } //END i iterations
  
  
  if(rank1==0 && (rank2+rank4)<=1) {
    if((rank2+rank4)==0)
      return -101;
    return 1;
  }
  if(rank1>0 && (rank2+rank4)==0)
    return 1;
  else
    return -1;
}

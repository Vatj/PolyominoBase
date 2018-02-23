#include <stochastic_model.hpp>

namespace Stochastic
{
  std::random_device rd;
  std::mt19937 RNG_Generator(rd());

  int Analyse_Genotype_Outcome(std::vector<int> genome, int N_Repeated_Checks, PhenotypeTable* pt,int seed) {
    Clean_Genome(genome,-1); 
    const unsigned int THRESHOLD_SIZE=(genome.size()*genome.size())/4;
    std::vector<int> Placed_Tiles_Check=Stochastic_Polyomino_Builder(genome,THRESHOLD_SIZE,seed,0);
    
    if(Placed_Tiles_Check.empty() || Placed_Tiles_Check.size()/4 > THRESHOLD_SIZE) //STERIC
      return -1;
    
    int dx=-1,dy=-1,dx_prime=-1,dy_prime=-1;
    std::vector<int> rs_lattice=Generate_Spatial_Occupancy(Placed_Tiles_Check,dx,dy,3),Placed_Tiles_Compare, rs_lattice_prime;
    
    for(int nth_repeat=1;nth_repeat<N_Repeated_Checks;++nth_repeat) {
      Placed_Tiles_Compare=Stochastic_Polyomino_Builder(genome,THRESHOLD_SIZE,seed,0);
      if(Placed_Tiles_Compare.empty() || Placed_Tiles_Compare.size()/4 > THRESHOLD_SIZE)
        return -1;
      if(Placed_Tiles_Check.size()!=Placed_Tiles_Compare.size())
        return -2;
      
      rs_lattice_prime=Generate_Spatial_Occupancy(Placed_Tiles_Compare,dx_prime,dy_prime,3);
      if(!Compare_Two_Polyominoes_Shapes(rs_lattice,dx,dy,rs_lattice_prime,dx_prime,dy_prime))
        return -2;
      if(N_Repeated_Checks-nth_repeat>1) {
        Placed_Tiles_Check=Placed_Tiles_Compare;
        rs_lattice=rs_lattice_prime;
        dx=dx_prime;
        dy=dy_prime;
      }
    }
    return pt->GetPhenotypeID(Phenotype{dx,dy,rs_lattice}); //Placed_Tiles_Check.size()/4;
  }
  
  std::vector<int> Stochastic_Polyomino_Builder(std::vector<int> genome, unsigned int THRESHOLD_SIZE, int initial_Tile,int initial_Rotation) {
    
    std::vector<int> Placed_Tiles{0,0,initial_Tile,initial_Rotation},Interacting_Faces; //DEFINED AS (X,Y,Tile Type Number, Tile Rotation[in CW rotation])
    //DEFINED AS (X,Y,Target Face colour, Target Face Index) 
    for(int face=0; face<4; ++face) {
      if(genome[initial_Tile*4+ (face-initial_Rotation+4)%4]!=0)
        Stochastic_Interacting_Adjacency(Interacting_Faces,genome[initial_Tile*4+ (face-initial_Rotation+4)%4],face,0,0);
    }
    
    while(!Interacting_Faces.empty()) {  
      for(int number_possibilities=(Interacting_Faces.size()/4)-1; number_possibilities >=0; --number_possibilities) {
        bool collision_Detected=false;
        for(std::vector<int>::iterator occupied_iter = Placed_Tiles.begin(); occupied_iter!=Placed_Tiles.end(); occupied_iter+=4) {
          if(Interacting_Faces[number_possibilities*4]==*occupied_iter && Interacting_Faces[number_possibilities*4+1]==*(occupied_iter+1)) {
            if(Interacting_Faces[number_possibilities*4+2]!=genome[*(occupied_iter+2)*4+(Interacting_Faces[number_possibilities*4+3]-*(occupied_iter+3)+4)%4])
              return {}; //Steric mismatch, reject empty vector as sign
            Interacting_Faces.erase(Interacting_Faces.end()-4,Interacting_Faces.end()); 
            collision_Detected=true;
            break;
          }
        }
        if(collision_Detected)
          continue;
        else {
          int conjugate_count=std::count(genome.begin(),genome.end(),Interacting_Faces[number_possibilities*4+2]);
          std::uniform_int_distribution<int> Random_Count(0,conjugate_count-1);
          int nth_conjugate=Random_Count(RNG_Generator);
          std::vector<int>::iterator current_conjugate=std::find(genome.begin(),genome.end(),Interacting_Faces[number_possibilities*4+2]);
          for(int conj_cnt=1;conj_cnt<=nth_conjugate;++conj_cnt) {
            current_conjugate=std::find(current_conjugate+1,genome.end(),Interacting_Faces[number_possibilities*4+2]);
          }
          int new_Tile=(current_conjugate-genome.begin())/4;
          int new_Face=(current_conjugate-genome.begin())%4;
          int rotation= (Interacting_Faces[number_possibilities*4+3]-new_Face+4)%4;
          Placed_Tiles.insert(Placed_Tiles.end(),{Interacting_Faces[number_possibilities*4],Interacting_Faces[number_possibilities*4+1],new_Tile,rotation});
          if(Placed_Tiles.size()/4 >THRESHOLD_SIZE) 
            return Placed_Tiles;
          int placed_X=Interacting_Faces[number_possibilities*4];
          int placed_Y=Interacting_Faces[number_possibilities*4+1];
          Interacting_Faces.erase(Interacting_Faces.begin()+number_possibilities*4,Interacting_Faces.begin()+number_possibilities*4+4);
          for(int face=1;face<4;++face) {
            int temp_Face=genome[new_Tile*4+(new_Face+face)%4];
            if(temp_Face!=0)
              Stochastic_Interacting_Adjacency(Interacting_Faces,temp_Face,(new_Face+face+rotation)%4,placed_X,placed_Y);
          }
          break;
        }
      }
    }
    return Placed_Tiles;
  }

  void Stochastic_Interacting_Adjacency(std::vector<int>& Interacting_Faces, int interacting_Face, int face_index, int X, int Y) {
    int X_OFFSET=0,Y_OFFSET=0;
    switch(face_index) {
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
    int conjugate_Face=Interaction_Matrix(interacting_Face);
    std::uniform_int_distribution<int> Random_Insertion(0,Interacting_Faces.size()/4);
    Interacting_Faces.insert(Interacting_Faces.begin()+Random_Insertion(RNG_Generator)*4,{X+X_OFFSET,Y+Y_OFFSET,conjugate_Face,(face_index+2)%4});
  }

  
  std::vector<int> Generate_Spatial_Occupancy(std::vector<int>& Placed_Tiles_Check, int& DELTA_X_Check,int& DELTA_Y_Check, int generate_mode) {
    std::vector<int> X_Locs_Check, Y_Locs_Check, Tile_Type_Check,Tile_Orientation_Check;
    for(std::vector<int>::iterator check_iter = Placed_Tiles_Check.begin();check_iter!=Placed_Tiles_Check.end();check_iter+=4) {
      X_Locs_Check.emplace_back(*check_iter);
      Y_Locs_Check.emplace_back(*(check_iter+1));
      Tile_Type_Check.emplace_back(*(check_iter+2));
      Tile_Orientation_Check.emplace_back(*(check_iter+3));
    }
    std::vector<int>::iterator LEFT_X_Check,RIGHT_X_Check,TOP_Y_Check,BOTTOM_Y_Check;
    std::tie(LEFT_X_Check,RIGHT_X_Check)=std::minmax_element(X_Locs_Check.begin(),X_Locs_Check.end());
    std::tie(BOTTOM_Y_Check,TOP_Y_Check)=std::minmax_element(Y_Locs_Check.begin(),Y_Locs_Check.end());
    DELTA_X_Check=*RIGHT_X_Check-*LEFT_X_Check+1;
    DELTA_Y_Check=*TOP_Y_Check-*BOTTOM_Y_Check+1;
    std::vector<int> Spatial_Occupancy_Check(DELTA_X_Check*DELTA_Y_Check);

    for(unsigned int tileIndex=0;tileIndex<X_Locs_Check.size();++tileIndex) {
      if(generate_mode==0) {
        Spatial_Occupancy_Check[(*TOP_Y_Check-Y_Locs_Check[tileIndex])*DELTA_X_Check + (X_Locs_Check[tileIndex]-*LEFT_X_Check)]=1;
      }
      if(generate_mode==1) {
        Spatial_Occupancy_Check[(*TOP_Y_Check-Y_Locs_Check[tileIndex])*DELTA_X_Check + (X_Locs_Check[tileIndex]-*LEFT_X_Check)]=Tile_Type_Check[tileIndex];
      }
      if(generate_mode==2) {
        Spatial_Occupancy_Check[(*TOP_Y_Check-Y_Locs_Check[tileIndex])*DELTA_X_Check + (X_Locs_Check[tileIndex]-*LEFT_X_Check)]=Tile_Orientation_Check[tileIndex];
      }
      if(generate_mode==3) {
        Spatial_Occupancy_Check[(*TOP_Y_Check-Y_Locs_Check[tileIndex])*DELTA_X_Check + (X_Locs_Check[tileIndex]-*LEFT_X_Check)]=1+Tile_Type_Check[tileIndex]*4+Tile_Orientation_Check[tileIndex];
      } 
    }
    return Spatial_Occupancy_Check;
  }
 
  bool Compare_Two_Polyominoes_Shapes(std::vector<int>& rs_lattice,int dx,int dy, std::vector<int>& rs_lattice_prime,int dx_prime,int dy_prime) {
    if(std::count(rs_lattice.begin(),rs_lattice.end(),0)!=std::count(rs_lattice_prime.begin(),rs_lattice_prime.end(),0))
      return false;
    if(dx==dx_prime && dy==dy_prime && dx==dy) { //bounding boxes match, symmetric
      if(rs_lattice==rs_lattice_prime) 
        return true;
      
      else {
        for(int rotation=0;rotation<3;++rotation) {
          Clockwise_Rotation(rs_lattice,dx,dy);
          if(rs_lattice==rs_lattice_prime) 
            return true;
        }
        return false;
      }
    }
    if(dx==dx_prime && dy==dy_prime) { //bounding boxes match, asymmetric
      if(rs_lattice==rs_lattice_prime) 
        return true;
      else {
        Clockwise_Pi_Rotation(rs_lattice,dx,dy);
        return rs_lattice==rs_lattice_prime ? true : false; 
      }
    }
    if(dx==dy_prime && dy==dx_prime) { //bounding boxes pi/2 off, asymmetric
      Clockwise_Rotation(rs_lattice,dx,dy);
      std::swap(dx,dy);
      if(rs_lattice==rs_lattice_prime)
        return true;
      else {
        Clockwise_Pi_Rotation(rs_lattice,dx,dy);
        return rs_lattice==rs_lattice_prime ? true : false; 
      }
    }
    return false;
  }
}


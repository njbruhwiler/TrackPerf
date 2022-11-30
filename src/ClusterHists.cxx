#include "TrackPerf/ClusterHists.hxx"

#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <UTIL/LCTrackerConf.h>


using namespace TrackPerf;

ClusterHists::ClusterHists()
{
  h_size_theta    = new TH2F("cluster_size_vs_theta" , ";Cluster #theta; Cluster size" , 100, -3.14,  3.14,  20,  -0.5,  19.5  ); 
  h_cluster_pos   = new TH2F("cluster_position"      , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_0 = new TH2F("cluster_position_0"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_1 = new TH2F("cluster_position_1"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_2 = new TH2F("cluster_position_2"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_3 = new TH2F("cluster_position_3"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
}

void ClusterHists::fill(const EVENT::TrackerHit* trkhit)
{
  //Calculating theta
  float x = trkhit->getPosition()[0];
  float y = trkhit->getPosition()[1];
  float z = trkhit->getPosition()[2];
  float r = sqrt(pow(x,2)+pow(y,2));
  float incidentTheta = std::atan(r/z);
  if(incidentTheta<0)
    incidentTheta += M_PI;

  //Calculating cluster size
  const lcio::LCObjectVec &rawHits = trkhit->getRawHits();
  float max = -1000000;
  float min = 1000000;
  for (size_t j=0; j<rawHits.size(); ++j) {
    lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit*>( rawHits[j] );
    const double *localPos = hitConstituent->getPosition();
    float x_local = localPos[0];
    float y_local = localPos[1];
    if (y_local < min){
      min = y_local;
      }
    else if (y_local > max){
      max = y_local;          
      }
    }
  float cluster_size = (max - min)+1;

  //Get hit subdetector/layer 
  std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
  uint32_t systemID = decoder(trkhit)["system"];
  uint32_t layerID = decoder(trkhit)["layer"];
      
  h_size_theta->Fill(incidentTheta, cluster_size);
  h_cluster_pos->Fill(z,r);
  if(layerID==0 or layerID==1){
    h_cluster_pos_0->Fill(z,r);
    }
  if(layerID==2 or layerID==3){
    h_cluster_pos_1->Fill(z,r);
    }
  if(layerID==4 or layerID==5){
    h_cluster_pos_2->Fill(z,r);
    }
  if(layerID==6 or layerID==7){
    h_cluster_pos_3->Fill(z,r);
    }   
}
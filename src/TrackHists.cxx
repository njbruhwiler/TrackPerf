#include "TrackPerf/TrackHists.hxx"

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <set>
#include <UTIL/LCTrackerConf.h>

using namespace TrackPerf;

TrackHists::TrackHists()
{
  h_pt     = new TH1F("reco_pt"     , ";Track p_{T} [GeV];Tracks [/0.1 GeV]"     , 100,  0   , 10   );
  h_lambda = new TH1F("reco_lambda" , ";Track #lambda; Tracks"                   , 100, -3.14,  3.14);
  h_phi    = new TH1F("reco_phi"    , ";Track #phi; Tracks"                      , 100, -3.14,  3.14);
  h_d0     = new TH1F("reco_d0"     , ";Track d_{0} [mm]; Tracks [/0.2 mm]"      , 100,-10   , 10   );
  h_z0     = new TH1F("reco_z0"     , ";Track z_{0} [mm]; Tracks [/0.2 mm]"      , 100,-10   , 10   );
  h_nhit   = new TH1F("reco_nhit"   , ";Track Hits; Tracks [/hit]"               , 20 ,-0.5  , 19.5 );
  h_chi2_spatial      = new TH1F("spatial_chi2"        , ";Spatial chi squared values; Tracks"            , 100,  0, 100);
  h_chi2_temp_avg     = new TH1F("temporal_chi2_avg"   , ";Average temporal chi squared values; Tracks"   , 100,  0, 10 );
  h_chi2_temp_max     = new TH1F("temporal_chi2_max"   , ";Max temporal chi squared values; Tracks"       , 100,  0, 10 );
  h_chi2_temp_std     = new TH1F("temporal_chi2_std"   , ";Std dev of temporal chi squared values; Tracks", 100,  0, 10 );
  h_chi2_temp_vtx     = new TH1F("temporal_chi2_vertex", ";Average temporal chi2 vertex detector; Tracks" , 2000, -.5, .5 );
  h_chi2_temp_inner   = new TH1F("temporal_chi2_inner" , ";Average temporal chi2 inner tracker; Tracks"   , 2000, -.5, .5 );
  h_chi2_temp_outer   = new TH1F("temporal_chi2_outer" , ";Average temporal chi2 outer tracker; Tracks"   , 2000, -.5, .5 ); 
  h_cov_d0_d0         = new TH1F("cov_d0_d0"        , ";Cov(d0, d0);Tracks"        , 100,  -1  ,  1   );
  h_cov_d0_phi        = new TH1F("cov_d0_phi"       , ";Cov(d0, phi);Tracks"       , 100,  -1  ,  1   );
  h_cov_d0_omega      = new TH1F("cov_d0_omega"     , ";Cov(d0, omega);Tracks"     , 100,  -1  ,  1   );
  h_cov_d0_z0         = new TH1F("cov_d0_z0"        , ";Cov(d0, z0);Tracks"        , 100,  -1  ,  1   );
  h_cov_d0_lambda     = new TH1F("cov_d0_lambda"    , ";Cov(d0, lambda);Tracks"    , 100,  -1  ,  1   );
  h_cov_phi_phi       = new TH1F("cov_phi_phi"      , ";Cov(phi, phi);Tracks"      , 100,  -1  ,  1   );
  h_cov_phi_omega     = new TH1F("cov_phi_omega"    , ";Cov(phi, omega);Tracks"    , 100,  -1  ,  1   );
  h_cov_phi_z0        = new TH1F("cov_phi_z0"       , ";Cov(phi, z0);Tracks"       , 100,  -1  ,  1   );
  h_cov_phi_lambda    = new TH1F("cov_phi_lambda"   , ";Cov(phi, lambda);Tracks"   , 100,  -1  ,  1   );
  h_cov_omega_omega   = new TH1F("cov_omega_omega"  , ";Cov(omega, omega);Tracks"  , 100,  -1  ,  1   );
  h_cov_omega_z0      = new TH1F("cov_omega_z0"     , ";Cov(omega, z0);Tracks"     , 100,  -1  ,  1   );
  h_cov_omega_lambda  = new TH1F("cov_omega_lambda" , ";Cov(omega, lambda);Tracks" , 100,  -1  ,  1   );
  h_cov_z0_z0         = new TH1F("cov_z0_z0"        , ";Cov(z0, z0);Tracks"        , 100,  -1  ,  1   );
  h_cov_z0_lambda     = new TH1F("cov_z0_lambda"    , ";Cov(z0, lambda);Tracks"    , 100,  -1  ,  1   );
  h_cov_lambda_lambda = new TH1F("cov_lambda_lambda", ";Cov(lambda, lambda);Tracks", 100,  -1  ,  1   );
  h_lambda_nhit = new TH2F("lambda_vs_nhit" , ";Track #lambda; Track Hits"       , 100, -3.14,  3.14,  20,  -0.5,  19.5  );
  h_pt_nhit     = new TH2F("pt_vs_nhit" , ";Track p_{T} [GeV]; Track Hits"       , 100, 0    ,  10  ,  20,  -0.5,  19.5  );
  h_pt_lambda   = new TH2F("pt_vs_lambda" , ";Track p_{T} [GeV]; Track #lambda"  , 100, 0    ,  10  , 100, -3.14,  3.14  );
  h_tempchi2_nhit = new TH2F("temporalchi2_vs_nhit", ";Average temporal chi squared values; Hits on track" , 20,0,2,20,-.5,19.5 );
  h_nhit_vtx  = new TH1F("reco_nhit_vtx" , ";Vertex detector track hits; Tracks [/hit]"  , 20 ,-0.5  , 19.5 );
  h_nhit_inner  = new TH1F("reco_nhit_inner" , ";Inner tracker track hits; Tracks [/hit]"  , 20 ,-0.5  , 19.5 );
  h_nhit_outer  = new TH1F("reco_nhit_outer" , ";Outer tracker track hits; Tracks [/hit]"  , 20 ,-0.5  , 19.5 );

  h_low_pt_lambda = new TH1F("low_pt_lambda", ";Track #lambda; Tracks" , 100, -3.14,  3.14);
  
  h_nhit_24_2 = new TH1F("reco_nhit_24_2" , ";Outer tracker barrel inner layer hits;" , 20, -0.5, 19.5);
  h_nhit_24_4 = new TH1F("reco_nhit_24_4" , ";Outer tracker barrel middle layer hits;" , 20, -0.5, 19.5);
  h_nhit_24_6 = new TH1F("reco_nhit_24_6" , ";Outer tracker barrel outer layer hits;" , 20, -0.5, 19.5);
}


void TrackHists::fill(const EVENT::Track* track)
{  
  float pt=fabs(0.3*_Bz/track->getOmega()/1000);
  h_pt    ->Fill(pt);
 
  //Calculating temporal chi squared value
  float chi2_temp_sum = 0;
  std::vector<float> chi2_temp_values;
  float chi2_temp_sum_vertex = 0;
  float chi2_temp_sum_inner = 0;
  float chi2_temp_sum_outer = 0;

  int nhits_24_2 = 0;
  int nhits_24_4 = 0;
  int nhits_24_6 = 0;

  int nhits = track->getTrackerHits().size();
  int vertex_nhits = track->getSubdetectorHitNumbers()[1]+track->getSubdetectorHitNumbers()[2];
  int inner_nhits = track->getSubdetectorHitNumbers()[3]+track->getSubdetectorHitNumbers()[4];
  int outer_nhits = track->getSubdetectorHitNumbers()[5]+track->getSubdetectorHitNumbers()[6];

  for (int i=0; i<(nhits-1); ++i)
    {float time_0 = track->getTrackerHits()[i]->getTime();
     float x_pos_0 = track->getTrackerHits()[i]->getPosition()[0];
     float y_pos_0 = track->getTrackerHits()[i]->getPosition()[1];
     float z_pos_0 = track->getTrackerHits()[i]->getPosition()[2];

     float time_1 = track->getTrackerHits()[i+1]->getTime();
     float x_pos_1 = track->getTrackerHits()[i+1]->getPosition()[0];
     float y_pos_1 = track->getTrackerHits()[i+1]->getPosition()[1];
     float z_pos_1 = track->getTrackerHits()[i+1]->getPosition()[2];

     float delta_r = sqrt( pow((x_pos_1-x_pos_0),2) + pow((y_pos_1-y_pos_0),2) + pow((z_pos_1-z_pos_0),2) );  // units mm

     float c = 300;  // units mm/ns  (speed of light)

     float expected_time = delta_r/c;  //units nanoseconds
     float observed_time = time_1 - time_0;  // units nanoseconds
     float chi2_temp = observed_time - expected_time;
     chi2_temp_values.push_back(pow(chi2_temp,2));
     chi2_temp_sum += pow(chi2_temp,2);

     //Find what subdetector the hit is on and fill subdetector hists
     std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
     UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
     uint32_t systemID_0 = decoder(track->getTrackerHits()[i])["system"];
     uint32_t systemID_1 = decoder(track->getTrackerHits()[i+1])["system"];
     uint32_t layerID = decoder(track->getTrackerHits()[i])["layer"];
     
     if(systemID_0 == systemID_1){
       if(systemID_0 == 1 or systemID_0 == 2){h_chi2_temp_vtx    ->Fill(chi2_temp);}
       if(systemID_0 == 3 or systemID_0 == 4){h_chi2_temp_inner  ->Fill(chi2_temp);}
       if(systemID_0 == 5 or systemID_0 == 6){h_chi2_temp_outer  ->Fill(chi2_temp);
         if(layerID == 0){nhits_24_2 += 1;}
         if(layerID == 1){nhits_24_4 += 1;}
         if(layerID == 2){nhits_24_6 += 1;}
         if(fabs(chi2_temp) >= 0.5){
          std::cout << "\ntime 0: " << track->getTrackerHits()[i]->getTime();
          std::cout << "\nx pos 0: " << track->getTrackerHits()[i]->getPosition()[0];
          std::cout << "\ny pos 0: " << track->getTrackerHits()[i]->getPosition()[1];
          std::cout << "\nz pos 0: " << track->getTrackerHits()[i]->getPosition()[2];
          std::cout << "\ntime 1: " << track->getTrackerHits()[i+1]->getTime();
          std::cout << "\nx pos 1: " << track->getTrackerHits()[i+1]->getPosition()[0];
          std::cout << "\ny pos 1: " << track->getTrackerHits()[i+1]->getPosition()[1];
          std::cout << "\nz pos 1: " << track->getTrackerHits()[i+1]->getPosition()[2];
          }
        }}}
  
  float chi2_temp_avg = chi2_temp_sum / (nhits-1);
  float chi2_temp_std = sqrt(fabs(chi2_temp_sum/nhits - pow(chi2_temp_avg,2)));

  h_chi2_temp_avg   ->Fill(chi2_temp_avg);
  h_chi2_temp_max   ->Fill(*max_element(chi2_temp_values.begin(), chi2_temp_values.end()));
  h_chi2_temp_std   ->Fill(chi2_temp_std);
  
  h_chi2_temp_inner ->Fill(chi2_temp_sum_inner/(inner_nhits-1));
  h_chi2_temp_outer ->Fill(chi2_temp_sum_outer/(outer_nhits-1));
  
  float lambda=std::atan(track->getTanLambda());
  h_lambda ->Fill(lambda);
  h_phi    ->Fill(track->getPhi());
  h_d0     ->Fill(track->getD0());
  h_z0     ->Fill(track->getZ0());
  if (pt < 2){h_low_pt_lambda ->Fill(lambda);}
  h_chi2_spatial->Fill(track->getChi2());

  h_nhit   ->Fill(nhits);

  h_cov_d0_d0          ->Fill(track->getCovMatrix()[0,0]);
  h_cov_d0_phi         ->Fill(track->getCovMatrix()[0,1]);
  h_cov_d0_omega       ->Fill(track->getCovMatrix()[0,2]);
  h_cov_d0_z0          ->Fill(track->getCovMatrix()[0,3]);
  h_cov_d0_lambda      ->Fill(track->getCovMatrix()[0,4]);
  h_cov_phi_phi        ->Fill(track->getCovMatrix()[1,1]);
  h_cov_phi_omega      ->Fill(track->getCovMatrix()[1,2]);
  h_cov_phi_z0         ->Fill(track->getCovMatrix()[1,3]);
  h_cov_phi_lambda     ->Fill(track->getCovMatrix()[1,4]);
  h_cov_omega_omega    ->Fill(track->getCovMatrix()[2,2]);
  h_cov_omega_z0       ->Fill(track->getCovMatrix()[2,3]);
  h_cov_omega_lambda   ->Fill(track->getCovMatrix()[2,4]);
  h_cov_z0_z0          ->Fill(track->getCovMatrix()[3,3]);
  h_cov_z0_lambda      ->Fill(track->getCovMatrix()[3,4]);
  h_cov_lambda_lambda  ->Fill(track->getCovMatrix()[4,4]);

  h_lambda_nhit ->Fill(lambda, nhits);
  h_pt_nhit     ->Fill(pt, nhits);
  h_pt_lambda   ->Fill(pt, lambda);
  h_tempchi2_nhit ->Fill(chi2_temp_avg, track->getTrackerHits().size());

  h_nhit_vtx ->Fill(vertex_nhits);
  h_nhit_inner ->Fill(inner_nhits);
  h_nhit_outer ->Fill(outer_nhits);

  h_nhit_24_2 ->Fill(nhits_24_2);
  h_nhit_24_4 ->Fill(nhits_24_4);
  h_nhit_24_6 ->Fill(nhits_24_6);
}


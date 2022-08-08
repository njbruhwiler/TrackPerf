#include "TrackPerf/TrackHists.hxx"

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <set>

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
  h_chi2_temp_vtx     = new TH1F("temporal_chi2_vertex", ";Average temporal chi2 vertex detector; Tracks" , 100, -5, 5  );
  h_chi2_temp_inner   = new TH1F("temporal_chi2_inner" , ";Average temporal chi2 inner tracker; Tracks"   , 100, -5, 5  );
  h_chi2_temp_outer   = new TH1F("temporal_chi2_outer" , ";Average temporal chi2 outer tracker; Tracks"   , 100, -5, 5  ); 
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
  h_nhit1  = new TH1F("reco_nhit_vtx" , ";Vertex detector track hits; Tracks [/hit]"  , 20 ,-0.5  , 19.5 );
  h_nhit2  = new TH1F("reco_nhit_inner" , ";Inner tracker track hits; Tracks [/hit]"  , 20 ,-0.5  , 19.5 );
  h_nhit3  = new TH1F("reco_nhit_outer" , ";Outer tracker track hits; Tracks [/hit]"  , 20 ,-0.5  , 19.5 );

  h_low_pt_lambda = new TH1F("low_pt_lambda", ";Track #lambda; Tracks" , 100, -3.14,  3.14);

  h_tempchi2_nhit = new TH2F("temporalchi2_vs_nhit", ";Average temporal chi squared values; Hits on track" , 20,0,2,20,-.5,19.5 );
}

void TrackHists::fill(const EVENT::Track* track)
{
  float pt=fabs(0.3*_Bz/track->getOmega()/1000);
  h_pt    ->Fill(pt);
 
  //Calculating temporal chi squared value
  float chi2_temp_sum = 0;
  float chi2_temp_squared_sum = 0;
  std::vector<float> chi2_temp_values;
  float chi2_temp_sum_vertex = 0;
  float chi2_temp_sum_inner = 0;
  float chi2_temp_sum_outer = 0;

  int nhits = track->getTrackerHits().size();
  int vertex_nhits = track->getSubdetectorHitNumbers()[1]+track->getSubdetectorHitNumbers()[2];
  int inner_nhits = track->getSubdetectorHitNumbers()[3]+track->getSubdetectorHitNumbers()[4];
  int outer_nhits = track->getSubdetectorHitNumbers()[5]+track->getSubdetectorHitNumbers()[6];

  for (int i=0; i<(nhits-1); ++i)
    {float time_0 = track->getTrackerHits()[i]->getTime();
     float x_pos_0 = track->getTrackerHits()[i]->getPosition()[0];
     float y_pos_0 = track->getTrackerHits()[i]->getPosition()[1];
     float r_pos_0 = sqrt(pow(x_pos_0,2) + pow(y_pos_0,2));

     float time_1 = track->getTrackerHits()[i+1]->getTime();
     float x_pos_1 = track->getTrackerHits()[i+1]->getPosition()[0];
     float y_pos_1 = track->getTrackerHits()[i+1]->getPosition()[1];
     float r_pos_1 = sqrt(pow(x_pos_1,2) + pow(y_pos_1,2));

     float delta_r = r_pos_1 - r_pos_0;  // units mm

     float c = 300;  // units mm/ns  (speed of light)

     float expected_time = delta_r/c;  //units nanoseconds
     float observed_time = time_1 - time_0;  // units nanoseconds
     float chi2_temp = observed_time - expected_time;
     chi2_temp_values.push_back(fabs(chi2_temp));
     chi2_temp_sum += fabs(chi2_temp);
     chi2_temp_squared_sum += pow(chi2_temp,2);

     if(i<vertex_nhits){chi2_temp_sum_vertex += chi2_temp;}
     if(i>vertex_nhits and i<(vertex_nhits+inner_nhits)){chi2_temp_sum_inner += chi2_temp;}
     if(i>(vertex_nhits+inner_nhits)){chi2_temp_sum_outer += chi2_temp;}
     }
  
  float chi2_temp_avg = chi2_temp_sum / (nhits-1);
  float chi2_temp_std = sqrt(fabs(chi2_temp_squared_sum/nhits - pow(chi2_temp_avg,2)));

  h_chi2_temp_avg   ->Fill(chi2_temp_avg);
  h_chi2_temp_max   ->Fill(*max_element(chi2_temp_values.begin(), chi2_temp_values.end()));
  h_chi2_temp_std   ->Fill(chi2_temp_std);
  h_chi2_temp_vtx   ->Fill(chi2_temp_sum_vertex/(vertex_nhits-1));
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

  h_nhit1 ->Fill(vertex_nhits);
  h_nhit2 ->Fill(inner_nhits);
  h_nhit3 ->Fill(outer_nhits);

  h_tempchi2_nhit ->Fill(chi2_temp_avg, track->getTrackerHits().size());
}


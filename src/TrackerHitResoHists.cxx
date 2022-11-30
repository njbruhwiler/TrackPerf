#include "TrackPerf/TrackerHitResoHists.hxx"

#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>


using namespace TrackPerf;

TrackerHitResoHists::TrackerHitResoHists()
{
  h_x_pull = new TH1F("x_pull", ";(x_{reco} - x_{truth})/#sigma_{x}; Hits" , 50, -2, 2); 
  h_y_pull = new TH1F("y_pull", ";(y_{reco} - y_{truth})/#sigma_{y}; Hits" , 50, -2, 2); 
  h_cov_x  = new TH1F("cov_x" , ";X variance; Hits"                        , 50, 0, 0.01);
  h_cov_y  = new TH1F("cov_y" , ";Y variance; Hits"                        , 50, 0, 0.01);
}

void TrackerHitResoHists::fill(const EVENT::TrackerHit* trkhit, const EVENT::SimTrackerHit* simtrkhit, IMPL::TrackerHitPlaneImpl* trkhitplane)
{
  //Calculating theta
  float x = trkhit->getPosition()[0];
  float y = trkhit->getPosition()[1];
  float z = trkhit->getPosition()[2];
  float r = sqrt(pow(x,2)+pow(y,2));
  float incidentTheta = std::atan(r/z);
  if(incidentTheta<0)
    incidentTheta += M_PI;

  float x_reco = trkhit->getPosition()[0];
  float y_reco = trkhit->getPosition()[1];
        
  float x_truth = simtrkhit->getPosition()[0];
  float y_truth = simtrkhit->getPosition()[1];

  if(x_reco-x_truth==0){
    std::cout<< "x truth: " << x_truth << std::endl;
    std::cout<< "y truth: " << y_truth << std::endl;
    std::cout<< "yreco-ytruth: " << y_reco-y_truth << std::endl;
    std::cout<< "incident theta: " << incidentTheta << std::endl << std::endl;
  }

  float sigma_x = trkhitplane->getdU();
  float sigma_y = trkhitplane->getdV();

  h_x_pull->Fill((x_reco-x_truth)/sigma_x);
  h_y_pull->Fill((y_reco-y_truth)/sigma_y);
  
  h_cov_x->Fill(sigma_x);
  h_cov_y->Fill(sigma_y);

}


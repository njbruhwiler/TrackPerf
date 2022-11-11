#include "TrackPerf/TrackPerfHistProc.hxx"

#include "TrackPerf/TrackHists.hxx"
#include "TrackPerf/TruthHists.hxx"
#include "TrackPerf/ResoHists.hxx"

#include <set>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <UTIL/LCTrackerConf.h>

#include <marlin/AIDAProcessor.h>

#include <AIDA/ITree.h>

TrackPerfHistProc aTrackPerfHistProc ;

TrackPerfHistProc::TrackPerfHistProc()
  : Processor("TrackPerfHistProc")
{  
  // modify processor description
  _description = "TrackPerfHistProc creates a series of output histograms for tracking performance studies." ;

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("MatchProb",
                             "Minimum matching probabilty to be considered a good track-mc match.",
                             _matchProb,
                             _matchProb);

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the MCParticle collection"  ,
			   _mcpColName,
			   _mcpColName
			   );

  registerInputCollection( LCIO::TRACK,
			   "TrackCollection" ,
			   "Name of the Track collection" ,
			   _trkColName,
			   _trkColName
			   );

  registerInputCollection( LCIO::LCRELATION,
			   "MCTrackRelationCollection" , 
			   "Name of LCRelation collection with track to MC matching",
			   _trkMatchColName,
			   _trkMatchColName
			   );

  registerInputCollection( LCIO::TRACKERHIT,
			   "VBTrackerHitsCollection" , 
			   "Name of vertex barrel tracker hits collection",
			   _vbtrkhitColName,
			   _vbtrkhitColName
			   );     

  registerInputCollection( LCIO::TRACKERHIT,
			   "IBTrackerHitsCollection" , 
			   "Name of inner barrel tracker hits collection",
			   _ibtrkhitColName,
			   _ibtrkhitColName
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "OBTrackerHitsCollection" , 
			   "Name of outer barrel tracker hits collection",
			   _obtrkhitColName,
			   _obtrkhitColName
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "VETrackerHitsCollection" , 
			   "Name of vertex endcap tracker hits collection",
			   _vetrkhitColName,
			   _vetrkhitColName
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "IETrackerHitsCollection" , 
			   "Name of inner endcap tracker hits collection",
			   _ietrkhitColName,
			   _ietrkhitColName
			   );     

  registerInputCollection( LCIO::TRACKERHIT,
			   "OETrackerHitsCollection" , 
			   "Name of outer endcap tracker hits collection",
			   _oetrkhitColName,
			   _oetrkhitColName
			   );  
}

void TrackPerfHistProc::init()
{
  // Print the initial parameters
  printParameters() ;

  // Create histograms
  AIDA::ITree* tree=marlin::AIDAProcessor::tree(this);
  marlin::AIDAProcessor::histogramFactory(this);

  // Create ROOT histograms, with location setup by the above factory
  tree->mkdir("all" ); tree->cd("all" );
  _allTracks =std::make_shared<TrackPerf::TrackHists>();
  _allTruths =std::make_shared<TrackPerf::TruthHists>();
  h_number_of_tracks = new TH1F("number_of_tracks", ";Number of seeds/tracks;Events", 100, 0, 300000);
  h_trackerhit_timing = new TH1F("time_of_flight", ";Time of flight for tracker hits [ns];Tracker hits;", 1000, -1, 1);
  tree->mkdir("../real"); tree->cd("../real");
  _realTracks=std::make_shared<TrackPerf::TrackHists>();
  _realTruths=std::make_shared<TrackPerf::TruthHists>();
  _realReso    =std::make_shared<TrackPerf::ResoHists>();
  h_relation_weight_real = new TH1F("relation_weight", ";Track-truth relation weight;", 100, 0, 2);
  tree->mkdir("../fake"); tree->cd("../fake");
  _fakeTracks=std::make_shared<TrackPerf::TrackHists>();
  h_number_of_fakes = new TH1F("number_of_fakes", ";Number of fake tracks;Events", 100,  0, 300000);
  h_relation_weight_fake = new TH1F("relation_weight", ";Track-truth relation weight;Tracks;", 100, 0, 2);
  tree->mkdir("../unmt"); tree->cd("../unmt");
  _unmtTruths=std::make_shared<TrackPerf::TruthHists>(); 
  tree->mkdir("../clusters" ); tree->cd("../clusters" );
  h_size_theta = new TH2F("cluster_size_vs_theta" , ";Cluster #theta; Cluster size" , 100, -3.14,  3.14,  20,  -0.5,  19.5  ); 
  h_cluster_pos = new TH2F("cluster_position", ";z; r" , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_0 = new TH2F("cluster_position_0", ";z; r" , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_1 = new TH2F("cluster_position_1", ";z; r" , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_2 = new TH2F("cluster_position_2", ";z; r" , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_3 = new TH2F("cluster_position_3", ";z; r" , 100, -500, 500, 100, 0, 200);

}

void TrackPerfHistProc::processRunHeader( LCRunHeader* /*run*/)
{ } 

void TrackPerfHistProc::processEvent( LCEvent * evt )
{
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.

  // MCParticles

  LCCollection* mcpCol  =evt->getCollection(_mcpColName);

  if( mcpCol->getTypeName() != lcio::LCIO::MCPARTICLE )
    { throw EVENT::Exception( "Invalid collection type: " + mcpCol->getTypeName() ) ; }

  std::set<const EVENT::MCParticle*> mcpSet;
  for(uint32_t i=0;i<mcpCol->getNumberOfElements();i++)
    {
      const EVENT::MCParticle *mcp=static_cast<const EVENT::MCParticle*>(mcpCol->getElementAt(i));

      if(mcp->getGeneratorStatus()!=1)
	{ continue; }

      if(mcp->getCharge()==0)
	{ continue; }

      if(mcp->isDecayedInTracker())
	{ continue; }

      // Tracker acceptance
      const double* mom=mcp->getMomentum();
      double pt=std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
      double lambda=std::atan2(mom[2],pt);
      //if(fabs(lambda)>75./180*3.14)
      if(fabs(lambda)>0.8)
	{ continue; }

      mcpSet.insert(mcp);
      _allTruths->fill(mcp);
    }

  // Tracks

  LCCollection* trkCol  =evt->getCollection(_trkColName);

  if( trkCol->getTypeName() != lcio::LCIO::TRACK )
    { throw EVENT::Exception( "Invalid collection type: " + trkCol->getTypeName() ) ; }

  std::set<const EVENT::Track*> trkSet;
  for(uint32_t i=0;i<trkCol->getNumberOfElements();i++)
    {
      const EVENT::Track *trk=static_cast<const EVENT::Track*>(trkCol->getElementAt(i));
      
      float lambda=std::atan(trk->getTanLambda());
      if(fabs(lambda)>0.8)
      {continue;}

      trkSet.insert(trk);
      _allTracks->fill(trk);
    }
    h_number_of_tracks->Fill(trkSet.size());

  // Relations

  LCCollection* relCol  =evt->getCollection(_trkMatchColName);

  if( relCol->getTypeName() != lcio::LCIO::LCRELATION )
    { throw EVENT::Exception( "Invalid collection type: " + relCol->getTypeName() ) ; }

  std::set<const EVENT::LCRelation*> relSet;
  for(uint32_t i=0;i<relCol->getNumberOfElements();i++)
    {
      const EVENT::LCRelation *rel=static_cast<const EVENT::LCRelation*>(relCol->getElementAt(i));

      relSet.insert(rel);
    }

  //Tracker hits

  LCCollection* vbtrkhitCol  =evt->getCollection(_vbtrkhitColName);
  LCCollection* ibtrkhitCol  =evt->getCollection(_ibtrkhitColName);
  LCCollection* obtrkhitCol  =evt->getCollection(_obtrkhitColName);
  LCCollection* vetrkhitCol  =evt->getCollection(_vetrkhitColName);
  LCCollection* ietrkhitCol  =evt->getCollection(_ietrkhitColName);
  LCCollection* oetrkhitCol  =evt->getCollection(_oetrkhitColName);

  for(int i=0; i<vbtrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(vbtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());

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

  for(int i=0; i<ibtrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(ibtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());}      
  for(int i=0; i<obtrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(obtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());}
  for(int i=0; i<vetrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(vetrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());}
  for(int i=0; i<ietrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(ietrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());}
  for(int i=0; i<oetrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(oetrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());}


  //
  // Loop over track to MC associations to save matched objects
  LCCollection* tr2mcCol=evt->getCollection(_trkMatchColName);
  if( tr2mcCol->getTypeName() != lcio::LCIO::LCRELATION )
    { throw EVENT::Exception( "Invalid collection type: "+ tr2mcCol->getTypeName() ); }

  for(int i=0; i<tr2mcCol->getNumberOfElements(); ++i)
    {
      const EVENT::LCRelation *rel=static_cast<const EVENT::LCRelation*>(tr2mcCol->getElementAt(i));
      const EVENT::MCParticle *mcp=static_cast<const EVENT::MCParticle*>(rel->getFrom());
      const EVENT::Track      *trk=static_cast<const EVENT::Track     *>(rel->getTo  ());

      if(mcpSet.count(mcp)==0)
	{ continue; } // truth particle not selected

      if(rel->getWeight()>_matchProb)
	{  
    if(trkSet.find(trk) != trkSet.end())
    {
	    _realTracks->fill(trk);
	    _realTruths->fill(mcp);
      _realReso  ->fill(trk,mcp);
      h_relation_weight_real ->Fill(rel->getWeight());

	    mcpSet.erase(mcp);
	    trkSet.erase(trk);
      relSet.erase(rel);
    }
	}
    }  

  //
  // Save unmatched objects
  for(const EVENT::MCParticle* mcp : mcpSet)
    { _unmtTruths->fill(mcp); }
  for(const EVENT::Track* trk : trkSet)
    { _fakeTracks->fill(trk);} 
  for(const EVENT::LCRelation* rel : relSet)
    { h_relation_weight_fake->Fill(rel->getWeight());}
  h_number_of_fakes->Fill(trkSet.size());
} 


void TrackPerfHistProc::check( LCEvent * /*evt*/ )
{ }

void TrackPerfHistProc::end()
{ }

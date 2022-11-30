#pragma once

#include <TH2.h>

#include "../ACTSTracking/ACTSTracking/GeometryIdMappingTool.hxx"

namespace EVENT
{
  class TrackerHit;
}

namespace TrackPerf
{
  //! Histograms for cluster analysis
  class ClusterHists
  {
  public:
    ClusterHists(const ClusterHists &) = delete ;
    ClusterHists& operator =(const ClusterHists &) = delete ;

    //! Initialize empty histograms
    ClusterHists() ;

    // Fill histograms with a single track hit
    void fill(const EVENT::TrackerHit* trkhit);

  private:
    TH2* h_size_theta;
    TH2* h_cluster_pos;
    TH2* h_cluster_pos_0;
    TH2* h_cluster_pos_1;
    TH2* h_cluster_pos_2;
    TH2* h_cluster_pos_3;
  };
}

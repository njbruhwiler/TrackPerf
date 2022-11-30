#pragma once

#include <TH1.h>
#include <TH2.h>
#include <TVector2.h>

namespace EVENT
{
  class Track;
  class MCParticle;
}

namespace TrackPerf
{
  //! Histograms for reconstructed tracks
  class TrackResoHists
  {
  public:
    TrackResoHists(const TrackResoHists &) = delete ;
    TrackResoHists& operator =(const TrackResoHists &) = delete ;

    //! Initialize empty histograms
    TrackResoHists() ;

    // Fill histograms with a single track
    void fill(const EVENT::Track* track, const EVENT::MCParticle* particle);

  private:
    //! magnetic field to use for curvature -> pT conversion
    float _Bz=3.57;

    //! Reconstructed track pT
    TH2* h_track_truth_pt;
    TH1* h_reso_pt_rel;
    TH1* h_reso_lambda;
    TH1* h_reso_phi;
  };
}
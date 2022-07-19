///////////////////////////////////////////
//Author: Paola Mastrapasqua UCLouvain(BE)
///////////////////////////////////////////

#ifndef CSCsegAnalyzer_H
#define CSCsegAnalyzer_H

// Base Class Headers

#include <string>
#include <iostream>
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoMuon/TrackingTools/interface/MuonSegmentMatcher.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/DYTInfo.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>

using namespace std;
using namespace edm;

class CSCsegAnalyzer : public edm::EDAnalyzer {
public: 

  explicit CSCsegAnalyzer(const edm::ParameterSet&);
  ~CSCsegAnalyzer();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

   inline bool isTrackerHighPtMuonNoVtx(const reco::Muon& muon){
      bool muID =   muon.isTrackerMuon() && muon.track().isNonnull() && (muon.numberOfMatchedStations() > 1);
      if(!muID) return false;

      bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
                  muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

      bool momQuality = muon.tunePMuonBestTrack()->ptError()/muon.tunePMuonBestTrack()->pt() < 0.3;

     // bool ip = fabs(muon.innerTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.innerTrack()->dz(vtx.position())) < 0.5;
      bool ip = true;
      
      return muID && hits && momQuality && ip;
      }

   inline bool isHighPtMuonNoVtx(const reco::Muon& muon/*, const reco::Vertex& vtx*/){
      if(!muon.isGlobalMuon()) return false;

      bool muValHits = ( muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0 ||
                      muon.tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits()>0 );

      bool muMatchedSt = muon.numberOfMatchedStations()>1;
      if(!muMatchedSt) {
        if( muon.isTrackerMuon() && muon.numberOfMatchedStations()==1 ) {
          if( muon.expectedNnumberOfMatchedStations()<2 ||
              !(muon.stationMask()==1 || muon.stationMask()==16) ||
              muon.numberOfMatchedRPCLayers()>2
            )
              muMatchedSt = true;
         }
      }

      bool muID = muValHits && muMatchedSt;

      bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
      muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

      bool momQuality = muon.tunePMuonBestTrack()->ptError()/muon.tunePMuonBestTrack()->pt() < 0.3;

      //bool ip = fabs(muon.innerTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.innerTrack()->dz(vtx.position())) < 0.5;
      bool ip = true;
   
      return muID && hits && momQuality && ip;
  
   }

  //Tight ID modified to remove Vtx constraint (as this gives low effciency for high pt Muon)
  //taken from https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonReco/src/MuonSelectors.cc

   inline bool isTightMuonNoVtx(const reco::Muon& muon) {
      if (!muon.isPFMuon() || !muon.isGlobalMuon()) return false;

      bool muID = muon::isGoodMuon(muon, muon::GlobalMuonPromptTight) && (muon.numberOfMatchedStations() > 1);

      bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

      //bool ip = std::abs(muon.muonBestTrack()->dxy(vtx.position())) < 0.2 && std::abs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;
      bool ip = true;
      return muID && hits && ip;
} 
      
  edm::InputTag MuonTags_; 
  edm::InputTag VtxTags_;
  edm::InputTag simHitsTags_;
  edm::InputTag genPTags_; 

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<edm::PSimHitContainer> simHitsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genPToken_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> CSCgeomToken_;
  std::unique_ptr<MuonSegmentMatcher> theMatcher_;

  string out, open;
  TFile* hFile;
  TTree* t;
  
  // ntuple content
  int nHitsCSC;
  int nSeg;
  int isLoose, isTight;
  int isHighPtMuonID, isTrackerHighPtMuonID;  
  int isCollision;
  double p,pt,eta,phi,timeDT,timeCSC,timeMuon,pt_err;
  double p_gen,pt_gen;
  double deltaR_recogen;
  
  std::vector<int>     nHits_perSeg;
  std::vector<double>  recoSeg_x, recoSeg_y;
  std::vector<double>  recoSeg_Err_xx, recoSeg_Err_xy, recoSeg_Err_yy;
  std::vector<double>  recoSeg_dx, recoSeg_dy, recoSeg_dz;
  std::vector<double>  recoSeg_Err_dxdzdxdz, recoSeg_Err_dxdzdydz, recoSeg_Err_dydzdydz;
  std::vector<double>  pullSeg, resSeg;
  std::vector<std::vector<double>> recoHits_x, recoHits_y;
  std::vector<std::vector<double>> recoHits_xGlobal, recoHits_yGlobal;
  std::vector<std::vector<double>> recoHits_Err_xx, recoHits_Err_xy, recoHits_Err_yy;
  std::vector<std::vector<double>> simHits_x, simHits_y, simHits_z, z_layer;
  std::vector<std::vector<double>> simHits_dx, simHits_dy,simHits_dz;
  std::vector<std::vector<int>> ringID, stationID, chamberID, layerID, zendcapID;

  //check
  std::vector<double> MUrecoHits_x, MUrecoHits_y;
  std::vector<double> MUrecoHits_Err_xx, MUrecoHits_Err_xy, MUrecoHits_Err_yy;
  std::vector<int> MUringID, MUstationID, MUchamberID, MUlayerID, MUzendcapID;
 
};
#endif

//Script for taking recoSegments matched with recoMuons
//Then pick the recoHits associated to the segment and find the corrisponding simHits 
//Finally save segments and hits info (position and direction) in a Ntupla 

/////////////////////////////////////////////
//Author: Paola Mastrapasqua UCLouvain(BE)
////////////////////////////////////////////

#include "CSCsegAnalyzer.h"

//
// constructors and destructor
//
CSCsegAnalyzer::CSCsegAnalyzer(const edm::ParameterSet& iConfig) 
  :
  MuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("Muons")),
  VtxTags_(iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVertex")),
  simHitsTags_(iConfig.getUntrackedParameter<edm::InputTag>("CSCg4SimHits")),
  genPTags_(iConfig.getUntrackedParameter<edm::InputTag>("g4genParticles")),
  out(iConfig.getParameter<string>("out")),
  open(iConfig.getParameter<string>("open"))
{
  edm::ConsumesCollector collector(consumesCollector());

  edm::ParameterSet matchParameters = iConfig.getParameter<edm::ParameterSet>("MatchParameters");
  theMatcher_ = std::make_unique<MuonSegmentMatcher>(matchParameters, collector);

  muonToken_ = consumes<reco::MuonCollection>(MuonTags_);
  vertexToken_ = consumes<reco::VertexCollection>(VtxTags_);
  simHitsToken_ = consumes<edm::PSimHitContainer>(simHitsTags_);
  genPToken_ = consumes<reco::GenParticleCollection>(genPTags_);
  CSCgeomToken_ = esConsumes();
}


CSCsegAnalyzer::~CSCsegAnalyzer()
{
  if (hFile!=0) {
    hFile->Close();
    delete hFile;
  }
}


// ------------ method called to for each event  ------------
void
CSCsegAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //CSC geometry needed to get z-position oh hits
  const CSCGeometry* geom = &iSetup.getData(CSCgeomToken_);

  // We need to get the primary vertex as it is needed for the Tight Muon ID definition
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexToken_,recVtxs);

  // Get the muon collection from the Event
  edm::Handle<reco::MuonCollection> MuCollection;
  iEvent.getByToken(muonToken_,MuCollection);
 
  //Get the simHits collection from the Event 
  edm::Handle<edm::PSimHitContainer> simHitsCollection;
  iEvent.getByToken(simHitsToken_, simHitsCollection);

  //Get the gen particle collection from the Event
  edm::Handle<reco::GenParticleCollection> genP; 
  iEvent.getByToken(genPToken_, genP);

  // Loop over reconstructed muons
  const reco::MuonCollection muonC = *(MuCollection.product());
  for(reco::MuonCollection::const_iterator imuon = muonC.begin(); imuon != muonC.end(); ++imuon){

    // fill muon ID
    isLoose = muon::isLooseMuon(*imuon);
    //tight ID modified to remove vtx constraint
    isTight = isTightMuonNoVtx(*imuon);
    isHighPtMuonID = isHighPtMuonNoVtx(*imuon);
    isTrackerHighPtMuonID = isTrackerHighPtMuonNoVtx(*imuon);
    
    // we need the muon to have the standalone track
    if ((imuon->standAloneMuon().isNull()) ) continue;

    // fill muon kinematics
    pt = imuon->pt();
    eta = imuon->eta();
    phi = imuon->phi();
    p = imuon->p();

    pt_err = imuon->bestTrack()->ptError();
   
    pt_err = 0.0;

    cout<<"Muon"<<endl;



    ////////////////////////////just a check - remove after //////////////////////////////////////////////
    MUrecoHits_x.clear();
    MUrecoHits_y.clear();
    MUrecoHits_Err_xx.clear();
    MUrecoHits_Err_xy.clear();
    MUrecoHits_Err_yy.clear();
    MUringID.clear();
    MUstationID.clear();
    MUchamberID.clear();
    MUlayerID.clear();
    MUzendcapID.clear();

    //std::vector<const TrackingRecHit*> rec_hits_fromMu = (std::vector<TrackingRecHit*>)((imuon->bestTrack())->recHits());
    for(trackingRecHit_iterator h_i = (imuon->standAloneMuon())->recHitsBegin() ;
           h_i != (imuon->standAloneMuon())->recHitsEnd(); ++h_i) {

        //cout<<"start outermost cycle"<<endl;

        if ( !(*h_i)->isValid()){ cout<< "Hit not valid"<<endl;  continue;}
        if ((*h_i)->geographicalId().det() != DetId::Muon) {continue;}
        if ((*h_i)->geographicalId().subdetId() != MuonSubdetId::CSC) {continue;}
        if (!(*h_i)->isValid()){ continue;}
        if ((*h_i)->recHits().size() < 2) { continue;}      

        //cout<<"recHits size: "<<(*h_i)->recHits().size()<<endl;
        std::vector<TrackingRecHit*> a = (*h_i)->recHits();
        
        for (trackingRecHit_iterator i_r = a.begin();
           i_r!=a.end();++i_r){

        if ( !(*i_r)->isValid()){ cout<< "Hit not valid"<<endl;  continue;}
        

        LocalPoint hlocal = (*i_r)->localPosition();
        LocalError hlocalerr = (*i_r)->localPositionError();
        
        MUrecoHits_x.push_back(hlocal.x());
        MUrecoHits_y.push_back(hlocal.y());
        MUrecoHits_Err_xx.push_back(hlocalerr.xx());
        MUrecoHits_Err_xy.push_back(hlocalerr.xy());
        MUrecoHits_Err_yy.push_back(hlocalerr.yy());

        CSCDetId idr = (CSCDetId)(*h_i)->geographicalId();
        MUringID.push_back(idr.ring());
        MUstationID.push_back(idr.station());
        MUchamberID.push_back(idr.chamber());
        MUlayerID.push_back(idr.layer());
        MUzendcapID.push_back(idr.zendcap());
        }
    }  
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    //initial definition of "closest" is really bad
    float minDeltaR = 999.0;

    //loop over the generated particles
    TLorentzVector v1,v2;   
     
    const reco::GenParticleCollection p = *(genP.product());
    reco::GenParticleCollection::const_iterator closest_genP;

    //cout<<"-------------------------------------------"<<endl;
    //cout<<"minDeltaR prima del for: "<< minDeltaR <<endl; 
    
    for (reco::GenParticleCollection::const_iterator p_it = p.begin(); p_it != p.end(); p_it++) {
        
      if(abs(p_it->pdgId()) == 13 && p_it->status() == 1){

        //cout<<"Eta of muon: " << p_it->eta()<<endl;
        v1.SetPxPyPzE((p_it->p4()).px(),(p_it->p4()).py(), (p_it->p4()).pz(), (p_it->p4()).E());
        v2.SetPxPyPzE((imuon->p4()).px(),(imuon->p4()).py(), (imuon->p4()).pz(), (imuon->p4()).E()); 
    	
        float tmp = (v1).DeltaR(v2);
        //cout<<"DeltaR with current genP: "<<endl;
    	//if it's closer, overwrite the definition of "closest"
    	if (tmp < minDeltaR) {
      		minDeltaR = tmp;
      		closest_genP = p_it;
    		}
        }
      }
    //cout<<" deltaR min : " << minDeltaR <<endl;
    //cout<<"---------------------------------"<<endl;
    pt_gen = closest_genP->pt();
    p_gen = closest_genP->p();
    deltaR_recogen = minDeltaR;   
   
    // ignore low pT muons
    if (pt<5.0) continue;

    nHitsCSC = 0;
    int iSeg = 0;
    nSeg = 0;
    int nHit = 0;
    nHits_perSeg.clear(); recoSeg_x.clear(); recoSeg_y.clear(); recoSeg_Err_xx.clear(); recoSeg_Err_xy.clear(); recoSeg_Err_yy.clear();
    recoSeg_dx.clear(); recoSeg_dy.clear(); recoSeg_dz.clear(); recoSeg_Err_dxdzdxdz.clear(); recoSeg_Err_dxdzdydz.clear(); recoSeg_Err_dydzdydz.clear();
    pullSeg.clear(); resSeg.clear(); recoHits_xGlobal.clear(); recoHits_yGlobal.clear();
    recoHits_x.clear(); recoHits_y.clear(); recoHits_Err_xx.clear(); recoHits_Err_xy.clear(); recoHits_Err_yy.clear();
    simHits_x.clear(); simHits_y.clear(); simHits_dx.clear(); simHits_dy.clear(); simHits_dz.clear();
    ringID.clear(); stationID.clear(); chamberID.clear(); layerID.clear(); zendcapID.clear(); z_layer.clear();
       
    // get CSC segements attached to the muon track
    std::vector<const CSCSegment*> matchedSegmentsCSC = theMatcher_->matchCSC(*(imuon->standAloneMuon()),iEvent);
    
    nSeg = matchedSegmentsCSC.size();
    //cout<<"Number of segments matched to a reco muon"<< nSeg <<endl; 
    
    float x_prec = 0.0;
    for (std::vector<const CSCSegment*>::iterator hiti = matchedSegmentsCSC.begin(); hiti!=matchedSegmentsCSC.end();++hiti) {
       
      if ( !(*hiti)->isValid()) {cout<<"seg not valid! "<<endl; continue;} 
      if ((*hiti)->localPosition().x() == x_prec) {cout<<"double counting segment! Removed!"<<endl; continue;}
      x_prec = (*hiti)->localPosition().x();
      //cout<<"Iter Seg: "<<iSeg<<endl;
      nHit = (*hiti)->nRecHits();
      nHitsCSC+= nHit;
      //cout<<"Hits in a segment: "<< nHit <<endl;
      
      // fill array with number of hits per segment (array length is the number of segment associated to a muon)
      nHits_perSeg.push_back(nHit);
      
      //position of the segment
      recoSeg_x.push_back((*hiti)->localPosition().x());
      recoSeg_y.push_back((*hiti)->localPosition().y());
      recoSeg_Err_xx.push_back((*hiti)->localPositionError().xx());
      recoSeg_Err_xy.push_back((*hiti)->localPositionError().xy());
      recoSeg_Err_yy.push_back((*hiti)->localPositionError().yy());
      
    
      LocalVector segmentdir = (*hiti)->localDirection();
      LocalError  segmentdirErr = (*hiti)->localDirectionError();
  
      //direction of the segment
      recoSeg_dx.push_back((*hiti)->localDirection().x());
      recoSeg_dy.push_back((*hiti)->localDirection().y());
      recoSeg_dz.push_back((*hiti)->localDirection().z());
      recoSeg_Err_dxdzdxdz.push_back((*hiti)->localDirectionError().xx());
      recoSeg_Err_dxdzdydz.push_back((*hiti)->localDirectionError().xy());
      recoSeg_Err_dydzdydz.push_back((*hiti)->localDirectionError().yy());
      
      //for each segment take the recoHits
      //vector of hits
      std::vector<const TrackingRecHit*> rec_hits = (*hiti)->recHits();


      recoHits_x.push_back(std::vector<double>()); recoHits_y.push_back(std::vector<double>()); recoHits_xGlobal.push_back(std::vector<double>()); recoHits_yGlobal.push_back(std::vector<double>());
      recoHits_Err_xx.push_back(std::vector<double>()); recoHits_Err_xy.push_back(std::vector<double>()); recoHits_Err_yy.push_back(std::vector<double>());
      simHits_x.push_back(std::vector<double>()); simHits_y.push_back(std::vector<double>());
      simHits_dx.push_back(std::vector<double>()); simHits_dy.push_back(std::vector<double>()); simHits_dz.push_back(std::vector<double>());
      ringID.push_back(std::vector<int>()); stationID.push_back(std::vector<int>()); chamberID.push_back(std::vector<int>()); layerID.push_back(std::vector<int>());
      zendcapID.push_back(std::vector<int>()); z_layer.push_back(std::vector<double>());

      //for each hit save position and find the corresponding simHit
      for (std::vector<const TrackingRecHit*>::iterator rec_hits_i = rec_hits.begin();
           rec_hits_i!=rec_hits.end();++rec_hits_i) {
      
      	if ( !(*rec_hits_i)->isValid()){ cout<< "Hit not valid"<<endl;  continue;}
      	
  	     	 
	LocalPoint rhitlocal = (*rec_hits_i)->localPosition();
	LocalError rhitlocalerr = (*rec_hits_i)->localPositionError();
	
        recoHits_x[iSeg].push_back(rhitlocal.x());
        recoHits_y[iSeg].push_back(rhitlocal.y());
        recoHits_Err_xx[iSeg].push_back(rhitlocalerr.xx());
	recoHits_Err_xy[iSeg].push_back(rhitlocalerr.xy());
	recoHits_Err_yy[iSeg].push_back(rhitlocalerr.yy());
	
        //taken from https://github.com/cms-sw/cmssw/blob/master/RecoLocalMuon/CSCSegment/test/CSCSegmentVisualise.cc#L210
        
        float r_closest = 9999;

        CSCDetId idrec = (CSCDetId)(*rec_hits_i)->geographicalId();

        GlobalPoint rhitglobal = (geom->chamber(idrec)->layer(idrec.layer())->surface()).toGlobal(rhitlocal);

        recoHits_xGlobal[iSeg].push_back(rhitglobal.x());
        recoHits_yGlobal[iSeg].push_back(rhitglobal.y());
       
        ringID[iSeg].push_back(idrec.ring());
        stationID[iSeg].push_back(idrec.station()); 
        chamberID[iSeg].push_back(idrec.chamber());
        layerID[iSeg].push_back(idrec.layer());
        zendcapID[iSeg].push_back(idrec.zendcap());

        edm::PSimHitContainer::const_iterator closest_simHit;

        int counter_enter = 0; 
        //loop over simHits
        const edm::PSimHitContainer simHits = *(simHitsCollection.product());
        for (edm::PSimHitContainer::const_iterator sim_it = simHits.begin(); sim_it != simHits.end(); sim_it++) {
            CSCDetId idsim = (CSCDetId)(*sim_it).detUnitId();

            if (idrec.endcap() == idsim.endcap() && 
                idrec.station() == idsim.station() && 
                idrec.ring() == idsim.ring() &&
                idrec.chamber() == idsim.chamber() && 
                idrec.layer() == idsim.layer()) {
          
                LocalPoint shitlocal = (*sim_it).localPosition();
                counter_enter=1;

                float dx2 = (rhitlocal.x() - shitlocal.x()) * (rhitlocal.x() - shitlocal.x());
                float dy2 = (rhitlocal.y() - shitlocal.y()) * (rhitlocal.y() - shitlocal.y());
                float dr2 = dx2 + dy2;
                if (dr2 < r_closest) {
                   r_closest = dr2;
                   closest_simHit = sim_it;
                   }
                }
          } //end loop on simHits
        if (counter_enter ==0) {
            cout<<"no simhit!?!"<<endl;
            }
      
        LocalPoint shitlocal = (*closest_simHit).localPosition();
        LocalVector shitlocaldir = (*closest_simHit).momentumAtEntry().unit();

	simHits_x[iSeg].push_back(shitlocal.x());
	simHits_y[iSeg].push_back(shitlocal.y());

        simHits_dx[iSeg].push_back(shitlocaldir.x());
        simHits_dy[iSeg].push_back(shitlocaldir.y());
        simHits_dz[iSeg].push_back(shitlocaldir.z());


        CSCDetId simId = (CSCDetId)(*closest_simHit).detUnitId();

        //Get z of layer where the hit is 
        z_layer[iSeg].push_back(geom->chamber(simId)->layer(simId.layer())->surface().toGlobal(LocalPoint(0, 0, 0)).z());

    
	}//end loop on recHits associated to CSCsegment 
	
    //get pull and resolution of segments
    double res = segmentdir.x()/segmentdir.z()-simHits_dx[iSeg][0]/simHits_dz[iSeg][0];
    resSeg.push_back(res);
    pullSeg.push_back(res/sqrt(segmentdirErr.xx()));

    iSeg++; 
    }//end loop on matched segments
    

    iSeg = 0;
    // fill the ntuple
    t->Fill();
  }//end loop on muons

}


// ------------ method called once each job just before starting event loop  ------------
void 
CSCsegAnalyzer::beginJob()
{
   
   hFile = new TFile( out.c_str(), open.c_str() );
   hFile->cd();

   t = new TTree("Muons", "Muons");
   
   // Loose muon ID bit
   t->Branch("isLoose", &isLoose, "isLoose/I");            

   // Tight muon ID bit
   t->Branch("isTight", &isTight, "isTight/I");

   // High Pt Muon ID
   t->Branch("isHighPtMuonID",&isHighPtMuonID, "isHighPtMuonID/I");
   t->Branch("isTrackerHighPtMuonID",&isTrackerHighPtMuonID, "isTrackerHighPtMuonID/I");

   // basic kinematical quantities for the muon
   t->Branch("pt", &pt, "pt/D");
   t->Branch("pt_err", &pt_err, "pt_err/D");
   t->Branch("eta", &eta, "eta/D");
   t->Branch("phi", &phi, "phi/D");
   t->Branch("p", &p, "p/D");
  
   //gen pt
   t->Branch("pt_gen", &pt_gen, "pt_gen/D");
   t->Branch("p_gen", &p_gen, "p_gen/D");

   //deltaR(gen,reco)
   t->Branch("deltaR_recogen", &deltaR_recogen, "deltaR_recogen/D");

   // numbers of hits
   t->Branch("nHitsCSC", &nHitsCSC, "nHitsCSC/I");

   // number of segments associated to a muon
   t->Branch("nSeg", &nSeg, "nSeg/I");
    
   // number of hits per segment 
   t->Branch("nHits_perSeg", &nHits_perSeg);  

   //segment position
   t->Branch("recoSeg_x", &recoSeg_x);
   t->Branch("recoSeg_y", &recoSeg_y);
   t->Branch("recoSeg_Err_xx", &recoSeg_Err_xx);
   t->Branch("recoSeg_Err_xy", &recoSeg_Err_xy);
   t->Branch("recoSeg_Err_yy", &recoSeg_Err_yy);
   
   //segment direction
   t->Branch("recoSeg_dx", &recoSeg_dx);
   t->Branch("recoSeg_dy", &recoSeg_dy);
   t->Branch("recoSeg_dz", &recoSeg_dz);
   t->Branch("recoSeg_Err_dxdzdxdz", &recoSeg_Err_dxdzdxdz);
   t->Branch("recoSeg_Err_dxdzdydz", &recoSeg_Err_dxdzdydz);
   t->Branch("recoSeg_Err_dydzdydz", &recoSeg_Err_dydzdydz);

   t->Branch("resSeg", &resSeg);
   t->Branch("pullSeg", &pullSeg);
   
   //hit position
   t->Branch("recoHits_x", &recoHits_x);
   t->Branch("recoHits_y", &recoHits_y);
   t->Branch("recoHits_xGlobal", &recoHits_xGlobal);
   t->Branch("recoHits_yGlobal", &recoHits_yGlobal);
   t->Branch("recoHits_Err_xx", &recoHits_Err_xx);
   t->Branch("recoHits_Err_xy", &recoHits_Err_xy);
   t->Branch("recoHits_Err_yy", &recoHits_Err_yy);

   //hit location on det
   t->Branch("ringID", &ringID);
   t->Branch("stationID", &stationID);
   t->Branch("chamberID", &chamberID);
   t->Branch("layerID", &layerID);
   t->Branch("zendcapID", &zendcapID);
  
   //sim hit position
   t->Branch("simHits_x", &simHits_x);
   t->Branch("simHits_y", &simHits_y);
   t->Branch("simHits_dx", &simHits_dx);
   t->Branch("simHits_dy", &simHits_dy);
   t->Branch("simHits_dz", &simHits_dz);
   t->Branch("z_layer", &z_layer); 
  
   //check
   t->Branch("MUrecoHits_x",&MUrecoHits_x);
   t->Branch("MUrecoHits_y",&MUrecoHits_y);
   t->Branch("MUrecoHits_Err_xx",&MUrecoHits_Err_xx);
   t->Branch("MUrecoHits_Err_xy",&MUrecoHits_Err_xy);
   t->Branch("MUrecoHits_Err_yy",&MUrecoHits_Err_yy);
   t->Branch("MUringID",&MUringID);
   t->Branch("MUstationID",&MUstationID);
   t->Branch("MUchamberID",&MUchamberID);
   t->Branch("MUlayerID",&MUlayerID);
   t->Branch("MUzendcapID",&MUzendcapID);

}

// ------------- write ntuple to file  ------------
void CSCsegAnalyzer::endJob() {

  hFile->cd();
  t->Write();
  hFile->Write();
  delete t;  

}

//define this as a plug-in
DEFINE_FWK_MODULE(CSCsegAnalyzer);

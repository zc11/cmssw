// -*- C++ -*-
//
// Package:    L1Trigger/L1TCalorimeter
// Class:      Stage2InputPatternWriter
// 
/**\class Stage2InputPatternWriter Stage2InputPatternWriter.cc L1Trigger/L1TCalorimeter/plugins/Stage2InputPatternWriter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  James Brooke
//         Created:  Tue, 11 Mar 2014 14:55:45 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/L1TObjects/interface/CaloParams.h"
#include "CondFormats/DataRecord/interface/L1TCaloParamsRcd.h"

#include "DataFormats/L1TCalorimeter/interface/CaloTower.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include <fstream>
#include <iostream>
#include <iomanip>

//
// class declaration
//

class Stage2InputPatternWriter : public edm::EDAnalyzer {
public:
  explicit Stage2InputPatternWriter(const edm::ParameterSet&);
  ~Stage2InputPatternWriter();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  edm::EDGetToken m_towerToken;

  std::string filename_;

  // constants
  unsigned nChan_;  // number of channels per quad
  unsigned nQuad_;
  unsigned nLink_;
  unsigned nFramePerEvent_;
  unsigned nFrame_;
  
  // data arranged by link and frame
  std::vector< std::vector<int> > data_;

  // map of towers onto links/frames
  std::map< int, int > map_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Stage2InputPatternWriter::Stage2InputPatternWriter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

  // register what you consume and keep token for later access:
  m_towerToken   = consumes<l1t::CaloTowerBxCollection>  (iConfig.getParameter<edm::InputTag>("towerToken"));

  filename_ = iConfig.getUntrackedParameter<std::string>("filename", "pattern.txt");

  nChan_ = 4;
  nQuad_ = 18;

  nFramePerEvent_ = 32;

  nFrame_ = 0;

  nLink_ = nChan_ * nQuad_;
  data_.resize(nLink_);
  LogDebug("L1TDebug") << "Preparing for " << nLink_ << " links" << std::endl;

}


Stage2InputPatternWriter::~Stage2InputPatternWriter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Stage2InputPatternWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  // get towers
  Handle< BXVector<l1t::CaloTower> > towHandle;
  iEvent.getByToken(m_towerToken,towHandle);

  std::vector<l1t::CaloTower> towers;
  
  for(std::vector<l1t::CaloTower>::const_iterator tower = towHandle->begin(0);
      tower != towHandle->end(0);
      ++tower) {
    towers.push_back(*tower);
  }
  

  // loop over frames
  for ( unsigned iFrame=0; iFrame<nFramePerEvent_; ++iFrame ) {
    
    // loop over links
    for ( unsigned iQuad=0; iQuad<nQuad_; ++iQuad ) {
      for ( unsigned iChan=0; iChan<nChan_; ++iChan ) {

	int data=0;

        // get tower ieta, iphi for link
	unsigned iLink = (iQuad*nChan_)+iChan;
	int ietaSgn = (iLink % 2==0 ? +1 : -1);
	int ieta = ietaSgn * iFrame;
	int iphi = (int) iLink/2;

	// get tower 1 data
	l1t::CaloTower tower = l1t::CaloTools::getTower(towers, ieta, iphi);
	data |= tower.hwPt() & 0xff;
	data |= (tower.hwEtRatio() & 0x7)<<8;
	data |= (tower.hwQual() & 0xf)<<12;

	// get tower 2
	ieta = ietaSgn * iFrame;
	iphi = ((int) iLink/2) + 1;
	tower = l1t::CaloTools::getTower(towers, ieta, iphi);
	data |= (tower.hwPt() & 0xff)<<16;
	data |= (tower.hwEtRatio() & 0x7)<<24;
	data |= (tower.hwQual() & 0xf)<<28;
	
	// add data to output
	data_.at(iLink).push_back( data );

      }

    }

      nFrame_++;

  }
      

}


// ------------ method called once each job just before starting event loop  ------------
void 
Stage2InputPatternWriter::beginJob()
{


}

// ------------ method called once each job just after ending the event loop  ------------
void 
Stage2InputPatternWriter::endJob() 
{

  LogDebug("L1TDebug") << "Read " << nFrame_ << " frames" << std::endl;

  // write file
  std::ofstream file( filename_ );

  file << "Board MP7_TEST" << std::endl;

  // quad/chan numbers
  file << " Quad/Chan : ";
  for ( unsigned i=0; i<nQuad_; ++i ) {
    for ( unsigned j=0; j<nChan_; ++j ) {
      file << "   q" << i << "c" << j << "   ";
    }
  }
  file << std::endl;

  // link numbers
  file << "      Link : ";
  for ( unsigned i=0; i<nQuad_; ++i ) {
    for ( unsigned j=0; j<nChan_; ++j ) {
      file << "    " << (i*nChan_)+j << "       ";
    }
  }

  file << std::endl;

  // then the data
  for ( unsigned iFrame=0; iFrame<nFrame_; ++iFrame ) {
    file << "Frame " << std::setw(4) << std::setfill('0') << iFrame << " : ";
    for ( unsigned iQuad=0; iQuad<nQuad_; ++iQuad ) {
      for ( unsigned iChan=0; iChan<nChan_; ++iChan ) {
	unsigned iLink = (iQuad*nChan_)+iChan;
	if (iLink<data_.size() && iFrame<data_.at(iLink).size()) {
	  file << "1v" << std::hex << std::setw(8) << std::setfill('0') << data_.at(iLink).at(iFrame) << " ";
	}
	else {
	  std::cerr << "Out of range : " << iLink << ", " << iFrame << std::endl;
	}
      }
    }
    file << std::endl;
  }

  file.close();
  
}

// ------------ method called when starting to processes a run  ------------
/*
void 
Stage2InputPatternWriter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Stage2InputPatternWriter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Stage2InputPatternWriter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Stage2InputPatternWriter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Stage2InputPatternWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Stage2InputPatternWriter);

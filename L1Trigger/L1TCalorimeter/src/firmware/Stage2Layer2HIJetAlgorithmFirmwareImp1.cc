// 
/// \class l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1
///
/// \author: Adam Elwood and Matthew Citron
///
/// Description: Implementation of Jad's asymmetric map overlap algorithm with donut subtraction

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2HIJetAlgorithmFirmware.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"
#include "L1Trigger/L1TCalorimeter/interface/AccumulatingSort.h"
#include "L1Trigger/L1TCalorimeter/interface/BitonicSort.h"
#include "CondFormats/L1TObjects/interface/CaloParams.h"
#include "TMath.h"

#include <vector>
#include <map>
#include <algorithm>
#include <math.h>

// jet mask, needs to be configurable at some point
// just a square for now
// for 1 do greater than, for 2 do greater than equal to


namespace l1t {
  inline bool operator > ( const l1t::Jet& a, l1t::Jet& b ) {
    return  a.hwPt() > b.hwPt();
  }
}

int mask2_[9][9] = {
  { 1,2,2,2,2,2,2,2,2 },
  { 1,1,2,2,2,2,2,2,2 },
  { 1,1,1,2,2,2,2,2,2 },
  { 1,1,1,1,2,2,2,2,2 },
  { 1,1,1,1,0,2,2,2,2 },
  { 1,1,1,1,1,2,2,2,2 },
  { 1,1,1,1,1,1,2,2,2 },
  { 1,1,1,1,1,1,1,2,2 },
  { 1,1,1,1,1,1,1,1,2 },
};


std::vector<l1t::Jet>::iterator start2_, end2_;

l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1::Stage2Layer2HIJetAlgorithmFirmwareImp1(CaloParamsHelper* params) :
  params_(params){}


l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1::~Stage2Layer2HIJetAlgorithmFirmwareImp1() {}

void l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1::processEvent(const std::vector<l1t::CaloTower> & towers,
							     std::vector<l1t::Jet> & jets,
							     std::vector<l1t::Jet> & alljets) {
  
  // find jets
  create(towers, jets, alljets, params_->jetPUSType());

  // jets accumulated sort
  accuSort(jets);

  unsigned int pos = 0;
  while(pos < jets.size()){
    bool isGood = true;

    l1t::Jet jet = jets.at(pos);

    for(unsigned int i = pos+1; i < jets.size(); ++i){
      if(jet.hwPt() < jets.at(i).hwPt()){
	jets.at(pos) = jets.at(i);
	jets.at(i) = jet;

	isGood = false;
	break;
      }
    }
    
    if(isGood) ++pos;
  }

  //  std::cout << "Jets after Stage2Layer2HIJetAlgorithmFirmwareImp1..." << std::endl;
  for(unsigned int i = 0; i < jets.size(); ++i){
    //    std::cout << " " << i << ": " << jets.at(i).hwPt() << ", " << jets.at(i).hwEta() << ", " << jets.at(i).hwPhi() << std::endl;
  }

}


void l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1::create(const std::vector<l1t::CaloTower> & towers,
						       std::vector<l1t::Jet> & jets, 
						       std::vector<l1t::Jet> & alljets, 
						       std::string PUSubMethod) {
  
  std::map<int, int> etaPUEstimate, etaPUN;
  if(PUSubMethod == "PhiRingPP" || PUSubMethod == "PhiRingPPExclude" || PUSubMethod == "PhiRingPPTower" || PUSubMethod == "PhiRingPPTowerMask" || PUSubMethod == "PhiRingHITower" || PUSubMethod == "PhiRingHIRegion") phiRingPUEstimate(&etaPUEstimate, &etaPUN, towers, PUSubMethod);

  // etaSide=1 is positive eta, etaSide=-1 is negative eta
  for (int etaSide=1; etaSide>=-1; etaSide-=2) {
    
    // the 4 groups of rings
    std::vector<int> ringGroup1, ringGroup2, ringGroup3, ringGroup4;
    for (int i=1; i<=CaloTools::mpEta(CaloTools::kHFEnd); i++) {
      if      ( ! ((i-1)%4) ) ringGroup1.push_back( i * etaSide );
      else if ( ! ((i-2)%4) ) ringGroup2.push_back( i * etaSide );
      else if ( ! ((i-3)%4) ) ringGroup3.push_back( i * etaSide );
      else if ( ! ((i-4)%4) ) ringGroup4.push_back( i * etaSide );
    }
    std::vector< std::vector<int> > theRings = { ringGroup1, ringGroup2, ringGroup3, ringGroup4 };
    
    // the 24 jets in this eta side
    std::vector<l1t::Jet> jetsHalf;
       
    // loop over the 4 groups of rings
    for ( unsigned ringGroupIt=1; ringGroupIt<=theRings.size(); ringGroupIt++ ) {
      
      // the 6 accumulated jets
      std::vector<l1t::Jet> jetsAccu;
     
      // loop over the 10 rings in this group
      for ( unsigned ringIt=0; ringIt<theRings.at(ringGroupIt-1).size(); ringIt++ ) {
	
	int ieta = theRings.at(ringGroupIt-1).at(ringIt);
       
	if(ieta == 28 || ieta == -28) continue;
	if(PUSubMethod == "PhiRingHITower" || PUSubMethod == "PhiRingHIRegion" || PUSubMethod == "PhiRingPPTowerMask"){
	  if(ieta == 27 || ieta == -27) continue;
	  else if(ieta == 26 || ieta == -26) continue;
	  else if(ieta == 25 || ieta == -25) continue;	  
	}

	// the jets in this ring
	std::vector<l1t::Jet> jetsRing;
	
	// loop over phi in the ring
	for ( int iphi=1; iphi<=CaloTools::kHBHENrPhi; ++iphi ) {
	  
	  // no more than 18 jets per ring
	  if (jetsRing.size()==18) break;
	  
	  // seed tower
	  const CaloTower& tow = CaloTools::getTower(towers, CaloTools::caloEta(ieta), iphi); 
	  
	  int seedEt = tow.hwPt();
	  int iEt = seedEt;
	  if(etaPUN[ieta] > 0){
	    if(PUSubMethod == "PhiRingPPTower" || PUSubMethod == "PhiRingPPTowerMask"|| PUSubMethod == "PhiRingHITower" || PUSubMethod == "PhiRingHIRegion") iEt -= (int)(etaPUEstimate[ieta]/etaPUN[ieta]);
	  }

	  bool vetoCandidate = false;
	  
	  // check it passes the seed threshold
	  if(iEt < floor(params_->jetSeedThreshold()/params_->towerLsbSum())) continue;
	  
	  // loop over towers in this jet
	  for( int deta = -4; deta < 5; ++deta ) {
	    for( int dphi = -4; dphi < 5; ++dphi ) {
	      
	      int towEt = 0;
	      int ietaTest = ieta+deta;
	      int iphiTest = iphi+dphi;

	      if(ietaTest == 28 || ietaTest == -28) continue;
	      if(PUSubMethod == "PhiRingHITower" || PUSubMethod == "PhiRingHIRegion" || PUSubMethod == "PhiRingPPTowerMask"){
		if(ietaTest == 27 || ietaTest == -27) continue;
		else if(ietaTest == 26 || ietaTest == -26) continue;
		else if(ietaTest == 25 || ietaTest == -25) continue;	  
	      }
	      
	      // wrap around phi
	      while ( iphiTest > CaloTools::kHBHENrPhi ) iphiTest -= CaloTools::kHBHENrPhi;
	      while ( iphiTest < 1 ) iphiTest += CaloTools::kHBHENrPhi;
	      
	      // wrap over eta=0
	      if (ieta > 0 && ietaTest <=0) ietaTest -= 1;
	      if (ieta < 0 && ietaTest >=0) ietaTest += 1;
	   
	      // check jet mask and sum tower et
	      const CaloTower& towTest = CaloTools::getTower(towers, CaloTools::caloEta(ietaTest), iphiTest);
	      towEt = towTest.hwPt();

	      if(PUSubMethod == "PhiRingPPTower" || PUSubMethod == "PhiRingPPTowerMask" || PUSubMethod == "PhiRingHITower" || PUSubMethod == "PhiRingHIRegion"){
		if(etaPUN[ietaTest] > 0) towEt -= (int)(etaPUEstimate[ietaTest]/etaPUN[ietaTest]);
		if(towEt < 0) towEt = 0;
	      }
						    
              if      (mask2_[8-(dphi+4)][deta+4] == 0) continue;
	      else if (mask2_[8-(dphi+4)][deta+4] == 1) vetoCandidate = (seedEt < towEt);
	      else if (mask2_[8-(dphi+4)][deta+4] == 2) vetoCandidate = (seedEt <= towEt);
	      
	      if (vetoCandidate) break;
	      else iEt += towEt;
	   
	    }
	    if(vetoCandidate) break; 
	  }
	
	  // add the jet to the list
	  if (!vetoCandidate) {

	    int rawEt = iEt;
	    int puEt(0);
	
	    math::XYZTLorentzVector p4;
	    int caloEta = CaloTools::caloEta(ieta);
	    l1t::Jet jet( p4, -999, caloEta, iphi, 0);

	    if(!params_->jetBypassPUS()){
	      if (PUSubMethod == "Donut") {
		puEt = donutPUEstimate(ieta, iphi, 5, towers);	    
		iEt -= puEt;
	      }
	      
	      if (PUSubMethod == "ChunkyDonut"){
		puEt = chunkyDonutPUEstimate(jet, 5, towers);
		iEt -= puEt;
	      }

	      if(PUSubMethod == "PhiRingPP"){
		puEt = applyPhiRingPUEstimate(etaPUEstimate, etaPUN, jet, iEt);
		iEt -= puEt;
	      }

	      if(PUSubMethod == "PhiRingPPExclude"){
		puEt = applyPhiRingPUEstimateExclude(etaPUEstimate, etaPUN, jet, iEt);
		iEt -= puEt;
	      }
	    }
	    
	    if (iEt<=0) continue;

	    // if tower Et is saturated, saturate jet Et
	    if (seedEt == CaloTools::kSatHcal || seedEt == CaloTools::kSatEcal || seedEt == CaloTools::kSatTower) iEt = CaloTools::kSatJet;

	    jet.setHwPt(iEt);
	    jet.setRawEt( (short int) rawEt);
	    jet.setSeedEt((short int) seedEt);
	    jet.setTowerIEta((short int) caloEta);
	    jet.setTowerIPhi((short int) iphi);
	    jet.setPUEt((short int) puEt);
	    

	    jetsRing.push_back(jet);
	    alljets.push_back(jet);
	    
	  }
	  
	}

	// sort these jets and keep top 6
	start2_ = jetsRing.begin();  
	end2_   = jetsRing.end();
	BitonicSort<l1t::Jet>(down, start2_, end2_);
	if (jetsRing.size()>6) jetsRing.resize(6);
	
	// update jets
	jets.insert(jets.end(),jetsRing.begin(),jetsRing.end());
	  
      }
    }
  } 
}


//Accumulating sort
void l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1::accuSort(std::vector<l1t::Jet> & jets){

  math::PtEtaPhiMLorentzVector emptyP4;
  l1t::Jet tempJet (emptyP4, 0, 0, 0, 0);
  std::vector< std::vector<l1t::Jet> > jetEtaPos( 41 , std::vector<l1t::Jet>(18, tempJet));
  std::vector< std::vector<l1t::Jet> > jetEtaNeg( 41 , std::vector<l1t::Jet>(18, tempJet));
  
  for (unsigned int iJet = 0; iJet < jets.size(); iJet++)
    {
      if (jets.at(iJet).hwEta() > 0) jetEtaPos.at(jets.at(iJet).hwEta()-1).at((jets.at(iJet).hwPhi()-1)/4) = jets.at(iJet);
      else  jetEtaNeg.at(-(jets.at(iJet).hwEta()+1)).at((jets.at(iJet).hwPhi()-1)/4) = jets.at(iJet);
    }
  
  AccumulatingSort <l1t::Jet> etaPosSorter(7);
  AccumulatingSort <l1t::Jet> etaNegSorter(7);
  std::vector<l1t::Jet> accumEtaPos;
  std::vector<l1t::Jet> accumEtaNeg;
    
  for( int ieta = 0 ; ieta < 41 ; ++ieta)
    {
      // eta +
      std::vector<l1t::Jet>::iterator start2_, end2_;
      start2_ = jetEtaPos.at(ieta).begin();  
      end2_   = jetEtaPos.at(ieta).end();
      BitonicSort<l1t::Jet>(down, start2_, end2_);
      etaPosSorter.Merge( jetEtaPos.at(ieta) , accumEtaPos );
      
      // eta -
      start2_ = jetEtaNeg.at(ieta).begin();  
      end2_   = jetEtaNeg.at(ieta).end();
      BitonicSort<l1t::Jet>(down, start2_, end2_);
      etaNegSorter.Merge( jetEtaNeg.at(ieta) , accumEtaNeg );
      
    }

  //check for 6 & 7th jets with same et and eta. Keep jet with larger phi
  
  if(accumEtaPos.at(6).hwPt()==accumEtaPos.at(5).hwPt() && accumEtaPos.at(6).hwEta()==accumEtaPos.at(5).hwEta()
     && accumEtaPos.at(6).hwPhi() > accumEtaPos.at(5).hwPhi()){
    accumEtaPos.at(5)=accumEtaPos.at(6);
  }
  if(accumEtaNeg.at(6).hwPt()==accumEtaNeg.at(5).hwPt() && accumEtaNeg.at(6).hwEta()==accumEtaNeg.at(5).hwEta()
     && accumEtaNeg.at(6).hwPhi() > accumEtaNeg.at(5).hwPhi()){
    accumEtaNeg.at(5)=accumEtaNeg.at(6);
  }
  
  //truncate
  accumEtaPos.resize(6);
  accumEtaNeg.resize(6);
  //Rework for more jets
  //  accumEtaPos.resize(24);
  //  accumEtaNeg.resize(24);

  // put all 12 candidates in the original jet vector, removing zero energy ones
  jets.clear();
  for (l1t::Jet accjet : accumEtaPos)
    {
      if (accjet.hwPt() > 0) jets.push_back(accjet);
    }
  for (l1t::Jet accjet : accumEtaNeg)
    {
      if (accjet.hwPt() > 0) jets.push_back(accjet);
    }
  
   
}



//A function to return the value for donut subtraction around an ieta and iphi position for donut subtraction
//Also pass it a vector to store the individual values of the strip for later testing
//The size is the number of ieta/iphi units out the ring is (ie for 9x9 jets, we want the 11x11 for PUS therefore we want to go 5 out, so size is 5)
int l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1::donutPUEstimate(int jetEta, 
							       int jetPhi, 
							       int size, 
							       const std::vector<l1t::CaloTower> & towers){

  //ring is a vector with 4 ring strips, one for each side of the ring
  std::vector<int> ring(4,0);

  int iphiUp = jetPhi + size;
  while ( iphiUp > CaloTools::kHBHENrPhi ) iphiUp -= CaloTools::kHBHENrPhi;
  int iphiDown = jetPhi - size;
  while ( iphiDown < 1 ) iphiDown += CaloTools::kHBHENrPhi;

  int ietaUp = jetEta+size;   //(jetEta + size > CaloTools::mpEta(CaloTools::kHFEnd)) ? 999 : jetEta+size;
  int ietaDown = jetEta-size; //(abs(jetEta - size) > CaloTools::mpEta(CaloTools::kHFEnd)) ? 999 : jetEta-size;

  for (int ieta = jetEta - size+1; ieta < jetEta + size; ++ieta)   
  {
    
    if (abs(ieta) > CaloTools::mpEta(CaloTools::kHFEnd) || abs(ieta) < 1) continue;
    int towerEta;
    
    if (jetEta > 0 && ieta <=0){
      towerEta = ieta-1;
    } else if (jetEta < 0 && ieta >=0){
      towerEta = ieta+1;
    } else {
      towerEta=ieta;
    }
  
    const CaloTower& tow = CaloTools::getTower(towers, CaloTools::caloEta(towerEta), iphiUp);
    int towEt = tow.hwPt();
    ring[0]+=towEt;
    
    const CaloTower& tow2 = CaloTools::getTower(towers, CaloTools::caloEta(towerEta), iphiDown);
    towEt = tow2.hwPt();
    ring[1]+=towEt;
    
  } 
  
  for (int iphi = jetPhi - size+1; iphi < jetPhi + size; ++iphi)   
  {
      
    int towerPhi = iphi;
    while ( towerPhi > CaloTools::kHBHENrPhi ) towerPhi -= CaloTools::kHBHENrPhi;
    while ( towerPhi < 1 ) towerPhi += CaloTools::kHBHENrPhi;
    
    const CaloTower& tow = CaloTools::getTower(towers, CaloTools::caloEta(ietaUp), towerPhi);
    int towEt = tow.hwPt();
    ring[2]+=towEt;
    
    const CaloTower& tow2 = CaloTools::getTower(towers, CaloTools::caloEta(ietaDown), towerPhi);
    towEt = tow2.hwPt();
    ring[3]+=towEt;
  } 
  
  //for the Donut Subtraction we only use the middle 2 (in energy) ring strips
  std::sort(ring.begin(), ring.end(), std::greater<int>());
  
  return 4*( ring[1]+ring[2] ); // This should really be multiplied by 4.5 not 4.
}

int l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1::chunkyDonutPUEstimate(l1t::Jet & jet, int size, 
								     const std::vector<l1t::CaloTower> & towers){
 
  int jetPhi = jet.hwPhi();
  int jetEta = CaloTools::mpEta(jet.hwEta());

   // ring is a vector with 4 ring strips, one for each side of the ring
  // order is PhiUp, PhiDown, EtaUp, EtaDown
  std::vector<int> ring(4,0);
  
  // number of strips in donut - should make this configurable
  int nStrips = 3;

  // loop over strips
  for (int stripIt=0; stripIt<nStrips; stripIt++) {

    int iphiUp   = jetPhi + size + stripIt;
    int iphiDown = jetPhi - size - stripIt;
    while ( iphiUp > CaloTools::kHBHENrPhi )   iphiUp   -= CaloTools::kHBHENrPhi;
    while ( iphiDown < 1 ) iphiDown += CaloTools::kHBHENrPhi;

    int ietaUp   = jetEta + size + stripIt;
    int ietaDown = jetEta - size - stripIt;
    if ( jetEta<0 && ietaUp>=0 )   ietaUp   += 1;
    if ( jetEta>0 && ietaDown<=0 ) ietaDown -= 1;
    
    // do PhiUp and PhiDown
    for (int ieta=jetEta-size+1; ieta<jetEta+size; ++ieta) {
      
      if (abs(ieta) > CaloTools::mpEta(CaloTools::kHFEnd)) continue;
      
      int towEta = ieta;
      if (jetEta>0 && towEta<=0) towEta-=1;
      if (jetEta<0 && towEta>=0) towEta+=1;
            
      const CaloTower& towPhiUp = CaloTools::getTower(towers, CaloTools::caloEta(towEta), iphiUp);
      int towEt = towPhiUp.hwPt();
      ring[0] += towEt;
            
      const CaloTower& towPhiDown = CaloTools::getTower(towers, CaloTools::caloEta(towEta), iphiDown);
      towEt = towPhiDown.hwPt();
      ring[1] += towEt;

    } 
    
    // do EtaUp
    for (int iphi=jetPhi-size+1; iphi<jetPhi+size; ++iphi) {
      
      if (abs(ietaUp) <= CaloTools::mpEta(CaloTools::kHFEnd)) {    
        int towPhi = iphi;
        while ( towPhi > CaloTools::kHBHENrPhi ) towPhi -= CaloTools::kHBHENrPhi;
        while ( towPhi < 1 ) towPhi += CaloTools::kHBHENrPhi;

        const CaloTower& towEtaUp = CaloTools::getTower(towers, CaloTools::caloEta(ietaUp), towPhi);
        int towEt = towEtaUp.hwPt();
        ring[2] += towEt;
      }

    }

    // do EtaDown
    for (int iphi=jetPhi-size+1; iphi<jetPhi+size; ++iphi) {
      
      if (abs(ietaDown) <= CaloTools::mpEta(CaloTools::kHFEnd)) {
        int towPhi = iphi;
        while ( towPhi > CaloTools::kHBHENrPhi ) towPhi -= CaloTools::kHBHENrPhi;
        while ( towPhi < 1 ) towPhi += CaloTools::kHBHENrPhi;
	
        const CaloTower& towEtaDown = CaloTools::getTower(towers, CaloTools::caloEta(ietaDown), towPhi);
        int towEt = towEtaDown.hwPt();
        ring[3] += towEt;
      }
     
    }     
    
    
  }
    
  // for donut subtraction we only use the middle 2 (in energy) ring strips
  // std::sort(ring.begin(), ring.end(), std::greater<int>());
  // return ( ring[1]+ring[2] ); 

  // use lowest 3 strips as PU estimate
  std::sort( ring.begin(), ring.end() );
  
  for(unsigned int i=0; i<4; ++i) jet.setPUDonutEt(i, (short int) ring[i]);

  return ( ring[0] + ring[1] + ring[2] );
  
}


void l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1::phiRingPUEstimate(std::map<int, int>* etaSum, std::map<int, int>* etaN, const std::vector<l1t::CaloTower> & towers, std::string algoStr)
{
  for(int i = -41; i <= 41; ++i){
    (*etaSum)[i] = 0;
    (*etaN)[i] = 0;
  }

  for(unsigned int i = 0; i < towers.size(); ++i){
    (*etaSum)[towers.at(i).hwEta()] += towers.at(i).hwPt();
    (*etaN)[towers.at(i).hwEta()] += 1;
  }

  if(algoStr == "PhiRingHITower" || algoStr == "PhiRingHIRegion"){
    int reVal = (*etaSum)[-1] + (*etaSum)[1];
    reVal /= 2;
    (*etaSum)[-1] = reVal;
    (*etaSum)[1] = reVal;

    if(algoStr == "PhiRingHITower"){
      for(int i = -41; i < -1; ++i){
	if(i%2 == 0) continue;
	reVal = (*etaSum)[i] + (*etaSum)[i+1];
	reVal /= 2;
	(*etaSum)[i] = reVal;
	(*etaSum)[i+1] = reVal;       
      }
      
      for(int i = 2; i < 41; ++i){
	if(i%2 != 0) continue;
	reVal = (*etaSum)[i] + (*etaSum)[i+1];
	reVal /= 2;
	(*etaSum)[i] = reVal;
	(*etaSum)[i+1] = reVal;       
      }
    }
    else if(algoStr == "PhiRingHIRegion"){
      int pos = 1;
      reVal = 0;
      for(int i = -41; i < -1; ++i){
	reVal += (*etaSum)[i];

	if(pos%4 == 0){
	  reVal /= 4;
	  (*etaSum)[i] = reVal;
	  (*etaSum)[i-1] = reVal;
	  (*etaSum)[i-2] = reVal;
	  (*etaSum)[i-3] = reVal;

	  reVal = 0;
	}
	pos++;
      }

      pos = 1;
      reVal = 0;
      for(int i = 2; i <= 42; ++i){
	reVal += (*etaSum)[i];

	if(pos%4 == 0){
	  reVal /= 4;
	  (*etaSum)[i] = reVal;
	  (*etaSum)[i-1] = reVal;
	  (*etaSum)[i-2] = reVal;
	  (*etaSum)[i-3] = reVal;
	  reVal = 0;
	}
	pos++;
      }
    }
  }

  return;
}


int l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1::applyPhiRingPUEstimate(std::map<int, int> etaSum, std::map<int, int> etaN, l1t::Jet & jet, int iEt)
{
  int jetEta = jet.hwEta();
  int etaRingSum = 0;
  int etaRingN = 0;

  //  std::cout << "Jet eta, pt: " << jetEta << ", " << iEt << std::endl;

  for(int i = jetEta-4; i <= jetEta+4; ++i){
    if(i == 28 || i == -28) continue;

    etaRingSum += etaSum[i];
    etaRingN += etaN[i];
  }

  //  std::cout << " ringSum, N: " << etaRingSum << ", " << etaRingN << std::endl;

  etaRingSum *= 81;
  etaRingSum /= etaRingN;
  return etaRingSum;
}


int l1t::Stage2Layer2HIJetAlgorithmFirmwareImp1::applyPhiRingPUEstimateExclude(std::map<int, int> etaSum, std::map<int, int> etaN, l1t::Jet & jet, int iEt)
{
  int jetEta = jet.hwEta();
  int etaRingSum = 0;
  int etaRingN = 0;

  //  std::cout << "Jet eta, pt: " << jetEta << ", " << iEt << std::endl;

  for(int i = jetEta-4; i <= jetEta+4; ++i){
    if(i == 28 || i == -28) continue;

    etaRingSum += etaSum[i];
    etaRingN += etaN[i];
  }

  //  std::cout << " ringSum, N: " << etaRingSum << ", " << etaRingN << std::endl;

  etaRingSum -= iEt;

  etaRingSum *= 81;
  etaRingSum /= (etaRingN - 81);
  return etaRingSum;
}

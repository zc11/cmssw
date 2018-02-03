///
/// Description: Firmware headers
///
/// Implementation:
///    Concrete firmware implementations
///
/// \author: Jim Brooke - University of Bristol
/// Modified: Adam Elwood - ICL

//
//

#ifndef Stage2Layer2HIJetAlgorithmFirmware_H
#define Stage2Layer2HIJetAlgorithmFirmware_H

#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2HIJetAlgorithm.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"

#include <map>

namespace l1t {

  // Imp1 is for v1 and v2
  class Stage2Layer2HIJetAlgorithmFirmwareImp1 : public Stage2Layer2HIJetAlgorithm {
  public:
    Stage2Layer2HIJetAlgorithmFirmwareImp1(CaloParamsHelper* params);
    virtual ~Stage2Layer2HIJetAlgorithmFirmwareImp1();
    virtual void processEvent(const std::vector<CaloTower> & towers,
			      std::vector<Jet> & jets, std::vector<Jet> & alljets);

    void create(const std::vector<CaloTower> & towers,
	                      std::vector<Jet> & jets, std::vector<Jet> & alljets, std::string PUSubMethod);

    void accuSort(std::vector<Jet> & jets);

    int donutPUEstimate(int jetEta, int jetPhi, int size,
                        const std::vector<l1t::CaloTower> & towers);

    int chunkyDonutPUEstimate(Jet & jet, int pos,
                              const std::vector<l1t::CaloTower> & towers);
    void phiRingPUEstimate(std::map<int, int>* etaSum, std::map<int, int>* etaN, const std::vector<l1t::CaloTower> & towers);
    int applyPhiRingPUEstimate(const std::map<int, int> etaSum, const std::map<int, int> etaN, l1t::Jet & jet, int iEt);
    int applyPhiRingPUEstimateExclude(const std::map<int, int> etaSum, const std::map<int, int> etaN, l1t::Jet & jet, int iEt);

  private:

    CaloParamsHelper* const params_;

  };

}

#endif

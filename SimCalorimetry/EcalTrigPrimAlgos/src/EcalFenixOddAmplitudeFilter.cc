#include "CondFormats/EcalObjects/interface/EcalTPGGroups.h"
#include "CondFormats/EcalObjects/interface/EcalTPGWeightGroup.h"
#include "CondFormats/EcalObjects/interface/EcalTPGWeightIdMap.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixOddAmplitudeFilter.h>
#include <iostream>
#include <string>
#include <fstream>

EcalFenixOddAmplitudeFilter::EcalFenixOddAmplitudeFilter(bool TPinfoPrintout, std::string oddWeightsTxtFile) : inputsAlreadyIn_(0), stripid_{0}, shift_(6), TPinfoPrintout_(TPinfoPrintout), oddWeightsTxtFile_(oddWeightsTxtFile) {}


EcalFenixOddAmplitudeFilter::~EcalFenixOddAmplitudeFilter() {}

int EcalFenixOddAmplitudeFilter::setInput(int input, int fgvb) {
  if (input > 0X3FFFF) {
    std::cout << "ERROR IN INPUT OF AMPLITUDE FILTER" << std::endl;
    return -1;
  }
  if (inputsAlreadyIn_ < 5) {
    buffer_[inputsAlreadyIn_] = input;
    fgvbBuffer_[inputsAlreadyIn_] = fgvb;
    inputsAlreadyIn_++;
  } else {
    for (int i = 0; i < 4; i++) {
      buffer_[i] = buffer_[i + 1];
      fgvbBuffer_[i] = fgvbBuffer_[i + 1];
    }
    buffer_[4] = input;
    fgvbBuffer_[4] = fgvb;
  }
  return 1;
}

void EcalFenixOddAmplitudeFilter::process(std::vector<int> &addout,
                                       std::vector<int> &output,
                                       std::vector<int> &fgvbIn,
                                       std::vector<int> &fgvbOut) {
  // test
  inputsAlreadyIn_ = 0;
  for (unsigned int i = 0; i < 5; i++) {
    buffer_[i] = 0;  // FIXME: 5
    fgvbBuffer_[i] = 0;
  }

  // test end

  for (unsigned int i = 0; i < addout.size(); i++) {
    if (i>=4){
      if(TPinfoPrintout_) std::cout<<i<<std::dec;//by  RK // need to add the boolean
    } 
    setInput(addout[i], fgvbIn[i]);
    process();
    output[i] = processedOutput_;
    fgvbOut[i] = processedFgvbOutput_;
  }
  // shift the result by 1!
  for (unsigned int i = 0; i < (output.size()); i++) {
    if (i != output.size() - 1) {
      output[i] = output[i + 1];
      fgvbOut[i] = fgvbOut[i + 1];
    } else {
      output[i] = 0;
      fgvbOut[i] = 0;
    }
  }
  return;
}

void EcalFenixOddAmplitudeFilter::process() {
  // UB FIXME: 5
  processedOutput_ = 0;
  processedFgvbOutput_ = 0;
  if (inputsAlreadyIn_ < 5)
    return;
  int output = 0;
  int fgvbInt = 0;

  if(TPinfoPrintout_) std::cout<<" "<<stripid_;
  for (int i = 0; i < 5; i++) {
    output += (weights_[i] * buffer_[i]) >> shift_;
    // if(TPinfoPrintout_) std::cout<<" "<<output<<std::dec; // Removing this because the information can be deduced from the 5 digis and 5 weights 
    if ((fgvbBuffer_[i] == 1 && i == 3) || fgvbInt == 1) {
      fgvbInt = 1;
    }
  }

  // by RK 
  if(TPinfoPrintout_){
    for (int i = 0; i < 5; i++) {
      std::cout<<" "<<weights_[i]<<std::dec;}
    for (int i = 0; i < 5; i++) {
      std::cout<<" "<<weights_[i]/64.0<<std::dec;}
    for (int i = 0; i < 5; i++) {
      std::cout<<" "<<buffer_[i]<<std::dec;} // Digis 
    // for (int i = 0; i < 5; i++) {
      // std::cout<<" "<<(weights_[i] * buffer_[i])<<std::dec; // Removing this because the information can be deduced from the 5 digis and 5 weights    
    // }
    std::cout << " ODD";
    std::cout<<std::endl;
      // -- by RK 
  }
  
  if (output < 0)
    output = 0;
  if (output > 0X3FFFF)
    output = 0X3FFFF;
  processedOutput_ = output;
  processedFgvbOutput_ = fgvbInt;
  //std::cout<<" output after full processing: "<<output<<std::endl; // by RK 
}

void EcalFenixOddAmplitudeFilter::setParameters(uint32_t raw,
                                             const EcalTPGWeightIdMap *ecaltpgWeightMap,
                                             const EcalTPGWeightGroup *ecaltpgWeightGroup) {
  stripid_ = raw;    // by RK  

  // Want to set Odd weights here 
  // Can see from header files that even amplitude weights come from CondFormats --> conditions database?   
  // For initial testing will load odd weights from text file 

  uint32_t params_[5];
  const EcalTPGGroups::EcalTPGGroupsMap &groupmap = ecaltpgWeightGroup->getMap();
  EcalTPGGroups::EcalTPGGroupsMapItr it = groupmap.find(raw);
  if (it != groupmap.end()) {
    uint32_t weightid = (*it).second;
    const EcalTPGWeightIdMap::EcalTPGWeightMap &weightmap = ecaltpgWeightMap->getMap();
    EcalTPGWeightIdMap::EcalTPGWeightMapItr itw = weightmap.find(weightid);
    (*itw).second.getValues(params_[0], params_[1], params_[2], params_[3], params_[4]);

    // we have to transform negative coded in 7 bits into negative coded in 32
    // bits maybe this should go into the getValue method??
    // std::cout << "peak flag settings" << std::endl;

    // for (int i = 0; i < 5; ++i) {
    //   weights_[i] = (params_[i] & 0x40) ? (int)(params_[i] | 0xffffffc0) : (int)(params_[i]);

    //   // Construct the peakFlag for sFGVB processing
    //   // peakFlag_[i] = ((params_[i] & 0x80) > 0x0) ? 1 : 0;
    //   // std::cout << " " << params_[i] << std::endl;
    //   // std::cout << " " << peakFlag_[i] << std::endl;
    // }

    // Setting Odd Weights via input text file here 
    std::fstream oddWeightsLine(oddWeightsTxtFile_, std::ios_base::in);
    oddWeightsLine >> weights_[0] >> weights_[1] >> weights_[2] >> weights_[3] >> weights_[4];  

  } else
    edm::LogWarning("EcalTPG") << " could not find EcalTPGGroupsMap entry for " << raw;
}

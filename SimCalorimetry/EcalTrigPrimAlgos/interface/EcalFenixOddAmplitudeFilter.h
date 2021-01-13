#ifndef ECAL_FENIX_ODD_AMPLITUDE_FILTER_H
#define ECAL_FENIX_ODD_AMPLITUDE_FILTER_H

#include <cstdint>
#include <vector>
#include <string>

class EcalTPGWeightIdMap;
class EcalTPGWeightGroup;

/**
 \ class EcalFenixOddAmplitudeFilter
 *  The purpose of this class is to implement the second (odd) ECAL FENIX amplitude filter 
 *  Derived from SimCalorimetry/EcalTrigPrimAlgos/src/EcalFenixAmplitudeFilter.cc, interface/EcalFenixAmplitudeFilter.h
 *  input: 18 bits
 *  output: 18 bits
 *
 */
class EcalFenixOddAmplitudeFilter {
private:
  int peakFlag_[5];
  int inputsAlreadyIn_;
  uint32_t stripid_; // by RK 
  int buffer_[5];
  int fgvbBuffer_[5];
  int weights_[5];
  int shift_;
  bool debug_; 
  bool TPinfoPrintout_; 
  std::string oddWeightsTxtFile_; // When including odd weights via a text file 
  int setInput(int input, int fgvb);
  void process();

  int processedOutput_;
  int processedFgvbOutput_;

public:
  EcalFenixOddAmplitudeFilter();
  EcalFenixOddAmplitudeFilter(bool TPinfoPrintout, std::string oddWeightsTxtFile);
  virtual ~EcalFenixOddAmplitudeFilter();
  virtual void process(std::vector<int> &addout,
                       std::vector<int> &output,
                       std::vector<int> &fgvbIn,
                       std::vector<int> &fgvbOut);
  void setParameters(uint32_t raw,
                     const EcalTPGWeightIdMap *ecaltpgWeightMap,
                     const EcalTPGWeightGroup *ecaltpgWeightGroup);
};

#endif

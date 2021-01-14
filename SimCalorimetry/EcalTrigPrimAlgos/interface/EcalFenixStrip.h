#ifndef ECAL_FENIXSTRIP_H
#define ECAL_FENIXSTRIP_H

#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixAmplitudeFilter.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixOddAmplitudeFilter.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixEtStrip.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixLinearizer.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixPeakFinder.h>
#include <SimCalorimetry/EcalTrigPrimAlgos/interface/EcalFenixStripFgvbEE.h>

#include "DataFormats/EcalDetId/interface/EcalTriggerElectronicsId.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include <DataFormats/EcalDigi/interface/EBDataFrame.h>
#include <DataFormats/EcalDigi/interface/EEDataFrame.h>
#include <string>

class EBDataFrame;
class EcalTriggerPrimitiveSample;
class EcalTPGSlidingWindow;
class EcalTPGFineGrainStripEE;
class EcalFenixStripFgvbEE;
class EcalFenixStripFormatEB;
class EcalFenixStripFormatEE;
class EcalTPGStripStatus;

/**
    \class EcalFenixStrip
    \brief class representing the Fenix chip, format strip
*/
class EcalFenixStrip {
public:
  // constructor, destructor
  EcalFenixStrip(const edm::EventSetup &setup,
                 const EcalElectronicsMapping *theMapping,
                 bool debug,
                 bool famos,
                 int maxNrSamples,
                 int nbMaxXtals,
                 bool TPinfoPrintout,
                 uint TPmode);
  virtual ~EcalFenixStrip();

private:
  const EcalElectronicsMapping *theMapping_;
  bool debug_;
  bool famos_;
  int nbMaxXtals_; 
  bool TPinfoPrintout_;
  uint TPmode_;
  std::vector<EcalFenixLinearizer *> linearizer_;
  EcalFenixAmplitudeFilter *amplitude_filter_;
  EcalFenixOddAmplitudeFilter *oddAmplitude_filter_;

  EcalFenixPeakFinder *peak_finder_;

  EcalFenixStripFormatEB *fenixFormatterEB_;

  EcalFenixStripFormatEE *fenixFormatterEE_;

  EcalFenixEtStrip *adder_;

  EcalFenixStripFgvbEE *fgvbEE_;

  // data formats for each event
  std::vector<std::vector<int>> lin_out_;
  std::vector<int> add_out_;
  std::vector<int> filt_out_;
  std::vector<int> peak_out_;
  std::vector<int> format_out_;
  std::vector<int> fgvb_out_;
  std::vector<int> fgvb_out_temp_;

  // Data formats for odd filter, as data path is duplicated for odd filter 
  std::vector<std::vector<int>> odd_lin_out_; // I'm not sure if this should be duplicated for the odd filter, but I think it's necessary for separate sFGVB vectors 
  std::vector<int> odd_add_out_; // I'm not sure if this should be duplicated for the odd filter, but I think it's necessary for separate sFGVB vectors 
  std::vector<int> odd_filt_out_; 
  std::vector<int> odd_peak_out_;
  // std::vector<int> odd_format_out_;
  std::vector<int> odd_fgvb_out_;
  std::vector<int> odd_fgvb_out_temp_; 

  const EcalTPGPedestals *ecaltpPed_;
  const EcalTPGLinearizationConst *ecaltpLin_;
  const EcalTPGWeightIdMap *ecaltpgWeightMap_;
  const EcalTPGWeightGroup *ecaltpgWeightGroup_;
  const EcalTPGOddWeightIdMap *ecaltpgOddWeightMap_;
  const EcalTPGOddWeightGroup *ecaltpgOddWeightGroup_;
  const EcalTPGSlidingWindow *ecaltpgSlidW_;
  const EcalTPGFineGrainStripEE *ecaltpgFgStripEE_;
  const EcalTPGCrystalStatus *ecaltpgBadX_;
  const EcalTPGStripStatus *ecaltpgStripStatus_;

  bool identif_;
  bool mydebug_;   // DP ADDED
  int a;
  
public:
  void setPointers(const EcalTPGPedestals *ecaltpPed,
                   const EcalTPGLinearizationConst *ecaltpLin,
                   const EcalTPGWeightIdMap *ecaltpgWeightMap,
                   const EcalTPGWeightGroup *ecaltpgWeightGroup,
                   const EcalTPGOddWeightIdMap *ecaltpgOddWeightMap,
                   const EcalTPGOddWeightGroup *ecaltpgOddWeightGroup,
                   const EcalTPGSlidingWindow *ecaltpgSlidW,
                   const EcalTPGFineGrainStripEE *ecaltpgFgStripEE,
                   const EcalTPGCrystalStatus *ecaltpgBadX,
                   const EcalTPGStripStatus *ecaltpgStripStatus) {
    ecaltpPed_ = ecaltpPed;
    ecaltpLin_ = ecaltpLin;
    ecaltpgWeightMap_ = ecaltpgWeightMap;
    ecaltpgWeightGroup_ = ecaltpgWeightGroup;
    ecaltpgOddWeightMap_ = ecaltpgOddWeightMap;
    ecaltpgOddWeightGroup_ = ecaltpgOddWeightGroup;
    ecaltpgSlidW_ = ecaltpgSlidW;
    ecaltpgFgStripEE_ = ecaltpgFgStripEE;
    ecaltpgBadX_ = ecaltpgBadX;
    ecaltpgStripStatus_ = ecaltpgStripStatus;
  }

  // main methods
  // process method is splitted in 2 parts:
  //   the first one is templated, the same except input
  //   the second part is slightly different for barrel/endcap


  template <class T>
  void process(const edm::EventSetup &, std::vector<const T> &, int nrxtals, std::vector<int> &out);
  void process_part2_barrel(uint32_t stripid,
                            const EcalTPGSlidingWindow *ecaltpgSlidW,
                            const EcalTPGFineGrainStripEE *ecaltpgFgStripEE);

  void process_part2_endcap(uint32_t stripid,
                            const EcalTPGSlidingWindow *ecaltpgSlidW,
                            const EcalTPGFineGrainStripEE *ecaltpgFgStripEE,
                            const EcalTPGStripStatus *ecaltpgStripStatus);

  // getters for the algorithms  ;

  EcalFenixLinearizer *getLinearizer(int i) const { return linearizer_[i]; }
  EcalFenixEtStrip *getAdder() const { return adder_; }
  EcalFenixAmplitudeFilter *getFilter() const { return amplitude_filter_; }
  EcalFenixOddAmplitudeFilter *getOddFilter() const { return oddAmplitude_filter_; }
  EcalFenixPeakFinder *getPeakFinder() const { return peak_finder_; }

  EcalFenixStripFormatEB *getFormatterEB() const { return fenixFormatterEB_; }
  EcalFenixStripFormatEE *getFormatterEE() const { return fenixFormatterEE_; }

  EcalFenixStripFgvbEE *getFGVB() const { return fgvbEE_; }

  void setbadStripMissing(bool flag) { identif_ = flag; }
  bool getbadStripMissing() const { return identif_; }

  // ========================= implementations
  // ==============================================================
  // void process(const edm::EventSetup &setup, std::vector<EBDataFrame> &samples, int nrXtals, std::vector<int> &out, bool OddFilter) {
  void process(const edm::EventSetup &setup, std::vector<EBDataFrame> &samples, int nrXtals, std::vector<int> &out) {
    
    // now call processing
    if (samples.empty()) {
      std::cout << " Warning: 0 size vector found in EcalFenixStripProcess!!!!!" << std::endl;
      return;
    }
    const EcalTriggerElectronicsId elId = theMapping_->getTriggerElectronicsId(samples[0].id());
    uint32_t stripid = elId.rawId() & 0xfffffff8;  // from Pascal

    // by RK 
    if (false){
    const EBDetId & id = samples[0].id();
    const EcalTrigTowerDetId towid = id.tower();
    for (int cryst = 0; cryst < nrXtals; cryst++) {
      //" EcalFenixStrip.h::stripid nsamples, tower eta, phi, xtal eta phi, 10 samples:  "
      std::cout<<stripid<<" "<<samples.size()<<" "<<towid.ieta()<<" "<<towid.iphi()<<"  "<<id.iphi()<<" "<<id.ieta(); // add to TPinfo 
      for (int i = 0; i < samples[cryst].size(); i++) {
	      std::cout<< " " << std::dec << samples[cryst][i].adc();
      }
      std::cout<<std::endl;
    }
  }
  
    identif_ = getFGVB()->getMissedStripFlag();

    process_part1(identif_,
                  samples,
                  nrXtals,
                  stripid,
                  ecaltpPed_,
                  ecaltpLin_,
                  ecaltpgWeightMap_,
                  ecaltpgWeightGroup_,
                  ecaltpgBadX_);  // templated part
    process_part2_barrel(stripid, ecaltpgSlidW_,
                         ecaltpgFgStripEE_);  // part different for barrel/endcap
    out = format_out_; 
    // if(OddFilter) out = odd_format_out_;
    // else out = format_out_;  
  }

  // void process(const edm::EventSetup &setup, std::vector<EEDataFrame> &samples, int nrXtals, std::vector<int> &out, bool OddFilter) {
  void process(const edm::EventSetup &setup, std::vector<EEDataFrame> &samples, int nrXtals, std::vector<int> &out) { // input even and odd outs? 
    // now call processing
    if (samples.empty()) {
      std::cout << " Warning: 0 size vector found in EcalFenixStripProcess!!!!!" << std::endl;
      return;
    }
    const EcalTriggerElectronicsId elId = theMapping_->getTriggerElectronicsId(samples[0].id());
    uint32_t stripid = elId.rawId() & 0xfffffff8;  // from Pascal

    identif_ = getFGVB()->getMissedStripFlag();

    process_part1(identif_,
                  samples,
                  nrXtals,
                  stripid,
                  ecaltpPed_,
                  ecaltpLin_,
                  ecaltpgWeightMap_,
                  ecaltpgWeightGroup_,
                  ecaltpgBadX_);  // templated part
    process_part2_endcap(stripid, ecaltpgSlidW_, ecaltpgFgStripEE_, ecaltpgStripStatus_);
    out = format_out_; // FIXME: timing
    // if(OddFilter) out = odd_format_out_;
    // else out = format_out_; // FIXME: timing
    return;
  }

  template <class T>
  void process_part1(int identif,
                     std::vector<T> &df,
                     int nrXtals,
                     uint32_t stripid,
                     const EcalTPGPedestals *ecaltpPed,
                     const EcalTPGLinearizationConst *ecaltpLin,
                     const EcalTPGWeightIdMap *ecaltpgWeightMap,
                     const EcalTPGWeightGroup *ecaltpgWeightGroup,
                     const EcalTPGCrystalStatus *ecaltpBadX) {
    if (debug_)
      std::cout << "\n\nEcalFenixStrip input is a vector of size: " << nrXtals << std::endl;
    
    if(debug_) std::cout << "Trigger Primitive mode: " << TPmode_ << std::endl; 

    // loop over crystals
    for (int cryst = 0; cryst < nrXtals; cryst++) {
      
      if (debug_) {
        std::cout << std::endl;
        std::cout << "cryst= " << cryst << " EBDataFrame/EEDataFrame is: " << std::endl;
        for (int i = 0; i < df[cryst].size(); i++) {
          std::cout << " " << std::dec << df[cryst][i].adc();
        }
        std::cout << std::endl;
      }
      // call linearizer
      this->getLinearizer(cryst)->setParameters(df[cryst].id().rawId(), ecaltpPed, ecaltpLin, ecaltpBadX);
      this->getLinearizer(cryst)->process(df[cryst], lin_out_[cryst]);
      
    }
    
    if (debug_) {
      std::cout << "output of linearizer is a vector of size: " << std::dec << lin_out_.size() << " of which used "
                << nrXtals << std::endl;
      for (int ix = 0; ix < nrXtals; ix++) {
        std::cout << "cryst: " << ix << "  value : " << std::dec << std::endl;
        std::cout << " lin_out[ix].size()= " << std::dec << lin_out_[ix].size() << std::endl;
        for (unsigned int i = 0; i < lin_out_[ix].size(); i++) {
          std::cout << " " << std::dec << (lin_out_[ix])[i];
        }
        std::cout << std::endl;
      }

      std::cout << std::endl;
    }

    // Now call the sFGVB - this is common between EB and EE!
    getFGVB()->setParameters(identif, stripid, ecaltpgFgStripEE_);
    getFGVB()->process(lin_out_, fgvb_out_temp_);

    // for odd filter 
    

    if (debug_) {
      std::cout << "output of strip fgvb is a vector of size: " << std::dec << fgvb_out_temp_.size() << std::endl;
      for (unsigned int i = 0; i < fgvb_out_temp_.size(); i++) {
        std::cout << " " << std::dec << (fgvb_out_temp_[i]);
      }
      std::cout << std::endl;
    }

    // call adder
    this->getAdder()->process(lin_out_, nrXtals, add_out_);  // add_out is of size SIZEMAX=maxNrSamples

    if (debug_) {
      std::cout << "output of adder is a vector of size: " << std::dec << add_out_.size() << std::endl;
      for (unsigned int ix = 0; ix < add_out_.size(); ix++) {
        std::cout << "Clock: " << ix << "  value : " << std::dec << add_out_[ix] << std::endl;
      }
      std::cout << std::endl;
    }
    
    if (famos_) {
      filt_out_[0] = add_out_[0];
      peak_out_[0] = add_out_[0];
      return;
    } else {

      // This is where the amplitude filters are called 
      // the TPmode flag will determine which are called. Even, Odd, or both. 

      if(TPmode_ == "Run2"){ // when calling both, add 'or' here for Run2 or both 

        // Call even amplitude filter 
        this->getFilter()->setParameters(stripid, ecaltpgWeightMap, ecaltpgWeightGroup);
        this->getFilter()->process(add_out_, filt_out_, fgvb_out_temp_, fgvb_out_);

        // Print out even filter ET and sfgvb values 
        if(debug_){
          std::cout << "output of EVEN filter is a vector of size: " << std::dec << filt_out_.size() << std::endl;
          for (unsigned int ix = 0; ix < filt_out_.size(); ix++) {
            std::cout << "Clock: " << ix << "  value : " << std::dec << filt_out_[ix] << std::endl;
          }
          std::cout << std::endl;    
          std::cout << "output of EVEN sfgvb after filter is a vector of size: " << std::dec << fgvb_out_.size() << std::endl;
          for (unsigned int ix = 0; ix < fgvb_out_.size(); ix++) {
            std::cout << "Clock: " << ix << "  value : " << std::dec << fgvb_out_[ix] << std::endl;
          }
          std::cout << std::endl;            
        }

          // Call peak finder on even filter output 
          this->getPeakFinder()->process(filt_out_, peak_out_);   

          // Print out even filter peak finder values 
          if(debug_){
            std::cout << "output of EVEN peakfinder is a vector of size: " << peak_out_.size() << std::endl;
            for (unsigned int ix = 0; ix < peak_out_.size(); ix++) {
              std::cout << "Clock: " << ix << "  value : " << peak_out_[ix] << std::endl;
            }
            std::cout << std::endl;               
          }      
        
      }

      else if(TPmode_ == "Odd"){

        // loop over crystals
        for (int cryst = 0; cryst < nrXtals; cryst++) {
          
          // if (debug_) {
          //   std::cout << std::endl;
          //   std::cout << "cryst= " << cryst << " EBDataFrame/EEDataFrame is: " << std::endl;
          //   for (int i = 0; i < df[cryst].size(); i++) {
          //     std::cout << " " << std::dec << df[cryst][i].adc();
          //   }
          //   std::cout << std::endl;
          // }
          // call linearizer
          this->getLinearizer(cryst)->setParameters(df[cryst].id().rawId(), ecaltpPed, ecaltpLin, ecaltpBadX);
          this->getLinearizer(cryst)->process(df[cryst], odd_lin_out_[cryst]);
          
        }

        getFGVB()->process(odd_lin_out_, odd_fgvb_out_temp_);
        this->getAdder()->process(odd_lin_out_, nrXtals, odd_add_out_);  // for odd filter input 

        // Call odd amplitude filter 
        this->getOddFilter()->setParameters(stripid, ecaltpgWeightMap, ecaltpgWeightGroup);
        this->getOddFilter()->process(add_out_, odd_filt_out_, odd_fgvb_out_temp_, odd_fgvb_out_);      
        // this->getOddFilter()->process(odd_add_out_, odd_filt_out_, odd_fgvb_out_temp_, odd_fgvb_out_);    

        // Print out odd filter ET and sfgvb values 
        if(debug_){
          std::cout << "output of ODD filter is a vector of size: " << std::dec << odd_filt_out_.size() << std::endl;
          for (unsigned int ix = 0; ix < odd_filt_out_.size(); ix++) {
            std::cout << "Clock: " << ix << "  value : " << std::dec << odd_filt_out_[ix] << std::endl;
          }
          std::cout << std::endl; 

          std::cout << "output of ODD sfgvb after filter is a vector of size: " << std::dec << odd_fgvb_out_.size() << std::endl;
          for (unsigned int ix = 0; ix < odd_fgvb_out_.size(); ix++) {
            std::cout << "Clock: " << ix << "  value : " << std::dec << odd_fgvb_out_[ix] << std::endl;
          }
          std::cout << std::endl;
        }        

        // Call peak finder on odd filter output 
        this->getPeakFinder()->process(odd_filt_out_, odd_peak_out_);        
          
        // Print out odd filter peak finder values 
        if (debug_) {
          std::cout << "output of ODD peakfinder is a vector of size: " << odd_peak_out_.size() << std::endl;
          for (unsigned int ix = 0; ix < odd_peak_out_.size(); ix++) {
            std::cout << "Clock: " << ix << "  value : " << odd_peak_out_[ix] << std::endl;
          }
          std::cout << std::endl;        
        }

      }
      return;
    }
  }
};
#endif

#ifndef DD_TECOptoHybAlgo_h
#define DD_TECOptoHybAlgo_h

#include <map>
#include <string>
#include <vector>
#include "DetectorDescription/Base/interface/DDTypes.h"
#include "DetectorDescription/Algorithm/interface/DDAlgorithm.h"

class DDTECOptoHybAlgo : public DDAlgorithm {
 public:
  //Constructor and Destructor
  DDTECOptoHybAlgo(); 
  virtual ~DDTECOptoHybAlgo();
  
  void initialize(const DDNumericArguments & nArgs,
                  const DDVectorArguments & vArgs,
                  const DDMapArguments & mArgs,
                  const DDStringArguments & sArgs,
                  const DDStringVectorArguments & vsArgs);

  void execute();

private:

  std::string              idNameSpace;    //Namespace of this and ALL parts
  std::string              childName;      //Child name
  double                   rmin;           //ICC piece Minimum R
  double                   rmax;           //          Maximum R
  double                   zpos;           //Z position of the OptoHybrid
  int                      startCopyNo;    //Start copy number
  std::vector<double>      angles;         //Angular position of Hybrid
};

#endif

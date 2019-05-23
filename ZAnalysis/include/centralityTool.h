#ifndef CENTTOOL
#define CENTTOOL

class CentralityTool{

  public:

  int getCentBinLow(int i);
  int getCentBinHigh(int i);

  int findCentIndex(int centLow, int centHigh);

  private:
  const int nCentBins = 22;
  int centBinLow[22] = {0,5, 10,20,30,40,50,60,70,80,90, 0,  0, 10,30, 50, 70, 70, 0, 30, 0, 50};
  int centBinHigh[22] ={5,10,20,30,40,50,60,70,80,90,100,100,10,30,50, 70, 90, 80, 30,100,50,100};

};

int CentralityTool::findCentIndex(int centLow, int centHigh){
  for(int i = 0; i<nCentBins; i++){
    if(centLow == centBinLow[i] && centHigh == centBinHigh[i]) return i;
  }
  std::cout << "Warning, could not find index!" << std::endl;
  return -1;
}

int CentralityTool::getCentBinLow(int i){
  if(i>=nCentBins){
    std::cout << "Out of range!" << std::endl;
    return -1;
  }
  
  return centBinLow[i];
}

int CentralityTool::getCentBinHigh(int i){
  if(i>=nCentBins){
    std::cout << "Out of range!" << std::endl;
    return -1;
  }
  
  return centBinLow[i];
}


#endif

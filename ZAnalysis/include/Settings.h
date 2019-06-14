#ifndef SETTINGS
#define SETTINGS


class Settings{
public:

  int nZMassBins = 30;
  float zMassRange[2] = {60,120};

  int nMuEffBinsEta = 48;
  static const int nMuEffBinsPt = 11;
  double muEffBinsPt[nMuEffBinsPt+1] = {20,25,30,35,40,45,50,60,80,120,160,200};

private:

};

#endif

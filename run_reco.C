int run_reco(char* fInputName, char* fOutputName, int fRecoIt=0)
{
  gInterpreter->LoadMacro("LightReco.C+");
  
  char fRecoName[200];

  if(fRecoIt==0)
    LightReco(fInputName,fOutputName);

  if(fRecoIt==1)
    LightReco(fInputName,fOutputName,fRecoIt);

  if(fRecoIt==2)
    LightReco(fInputName,fOutputName,fRecoIt,"f.root");

  if(fRecoIt==12)
  {
    sprintf(fRecoName,"reco1_%s",fOutputName);
    LightReco(fInputName,fRecoName,1);
    LightReco(fInputName,fOutputName,2,fRecoName);
  }
}


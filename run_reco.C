// fRecoIt=0 central events only - this will reproduce JINST paper results
// fRecoIt=1 1st round of reco to get preliminary estimate on vertex position
//           (selecting early emitted photons in sphere segments)
// fRecoIt=2 2nd round reco for precise vtx - applies position-dependent time cut
//           (needs output of 1st round of reco)
// fRecoIt=12 does 1 and 2 in a single run

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


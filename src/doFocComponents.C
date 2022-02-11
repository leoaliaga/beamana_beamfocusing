
#include "ParseFocusing.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/dk2nu.h"

#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"

using namespace NeutrinoFocAna;

const int Nfoc   = 6;
const int Nhel   = 4;
const int Nebin  = 400;
const int Npzbin = 160;
const int Nptbin = 50;
const std::string foc[Nfoc] = {"Unfocused","H2only","H1only",
			       "Underfocused","Overfocused","Other"};
const std::string hel[Nhel] = {"numu","anumu","nue","anue"};
const double emin  = 0;
const double emax  = 40;
const double pzmin = 0;
const double pzmax = 80;
const double ptmin = 0;
const double ptmax = 1.;

int idx_hel(int pdg);
int idx_foc(std::string foc_comp);

void doFocComponents(const char* inputFile, const char* cdet, const char* outputFile){
  
  std::cout<<"starting()"<<std::endl;
  std::cout<<"inputFile : "<<inputFile <<std::endl;
  std::cout<<"cdet      : "<<cdet      <<std::endl;
  std::cout<<"outputFile: "<<outputFile<<std::endl;

  int idet = atoi(cdet);
  bool is_single_file = (std::string(inputFile).find(".root") < 10000);
  bool is_many_files  = (std::string(inputFile).find(".txt")  < 10000);
  if(!is_single_file && !is_many_files){
    std::cout<<"Ont 1 single root file or a text file with a list of root files are allowed"<<std::endl;
    exit (1);
  }
  
  ParseFocusing* parsefoc = ParseFocusing::getInstance();
  
  TFile* fOut = new TFile(outputFile,"recreate");
  std::cout<<"File name: "<<fOut->GetName()<<std::endl;
  TH1D* hflux_tot[Nhel];
  TH1D* hflux_pinupar[Nhel];
  TH1D* hflux_kanupar[Nhel];
  TH1D* hflux_foc_pinupar[Nhel][Nfoc];
  TH1D* hflux_foc_kanupar[Nhel][Nfoc];
  TH2D* hpzpt_foc_pinupar[Nhel][Nfoc];
  TH2D* hpzpt_foc_kanupar[Nhel][Nfoc];
  for(int i=0;i<Nhel;i++){
    hflux_tot[i]     = new TH1D(Form("hflux_tot_%s"    ,hel[i].c_str()),"",Nebin,emin,emax);
    hflux_pinupar[i] = new TH1D(Form("hflux_pinupar_%s",hel[i].c_str()),"",Nebin,emin,emax);
    hflux_kanupar[i] = new TH1D(Form("hflux_kanupar_%s",hel[i].c_str()),"",Nebin,emin,emax);
    for(int j=0;j<Nfoc;j++){
      hflux_foc_pinupar[i][j] = new TH1D(Form("hflux_foc_pinupar_%s_%s",hel[i].c_str(),foc[j].c_str()),"",Nebin,emin,emax);
      hflux_foc_kanupar[i][j] = new TH1D(Form("hflux_foc_kanupar_%s_%s",hel[i].c_str(),foc[j].c_str()),"",Nebin,emin,emax);
      hpzpt_foc_pinupar[i][j] = new TH2D(Form("hpzpt_foc_pinupar_%s_%s",hel[i].c_str(),foc[j].c_str()),"",Npzbin,pzmin,pzmax,Nptbin,ptmin,ptmax);
      hpzpt_foc_kanupar[i][j] = new TH2D(Form("hpzpt_foc_kanupar_%s_%s",hel[i].c_str(),foc[j].c_str()),"",Npzbin,pzmin,pzmax,Nptbin,ptmin,ptmax);
    }
  }
  
  //Loading ntuples:
  TChain* chain_evts   = new TChain("dk2nuTree");  
  TChain* chain_meta   = new TChain("dkmetaTree");  
  bsim::Dk2Nu*  dk2nu  = new bsim::Dk2Nu;  
  bsim::DkMeta* dkmeta = new bsim::DkMeta;

  if(is_single_file){
    std::cout<<" Adding ntuple at "<<inputFile<<std::endl;
    chain_evts->Add(inputFile);
    chain_meta->Add(inputFile);
  }
  if(is_many_files){
    std::ifstream ifs;
    ifs.open(inputFile);
    std::string line;
    int counter = 0;
    while (ifs.good()) {
       getline(ifs,line);
       if(line.find(".root")>10000)continue;
       chain_evts->Add(line.c_str());
       if(counter==0)chain_meta->Add(line.c_str());
       std::cout<<"Adding ntuple at : "<<line<<std::endl;
       counter++;
    }
    ifs.close();  
  }
  
  chain_evts->SetBranchAddress("dk2nu",&dk2nu);
  chain_meta->SetBranchAddress("dkmeta",&dkmeta);
  
  chain_meta->GetEntry(0); 
  std::cout<<"=> Detector: "<< (dkmeta->location)[idet].name <<std::endl;
  int nentries  = chain_evts->GetEntries();
  std::cout<<"N of entries: "<<nentries<<std::endl;

  for(int ii=0;ii<nentries;ii++){  
    if(ii%30000==0)std::cout<<ii/1000<<" k evts"<<std::endl;
     chain_evts->GetEntry(ii);
    double fluxWGT  = ( (dk2nu->nuray)[idet].wgt )*(dk2nu->decay.nimpwt)/3.1416;
    double nuenergy = (dk2nu->nuray)[idet].E; 
    
    parsefoc->CalculateComponents(dk2nu,dkmeta,"NuMI");
    int ihel = idx_hel(dk2nu->decay.ntype);
    int ifoc = idx_foc(parsefoc->GetCategory());
    
    if(ihel<0 || ifoc<0){
      std::cout<<"something wrong, exiting: (ihel,ifoc) = ("<<ihel<<", "<<ifoc<<")"<<std::endl;
      exit (1);
    }
    
    //parent:
    bool is_pinupar = ( abs(dk2nu->decay.ptype) == 211 );
    bool is_kanupar = ( abs(dk2nu->decay.ptype) == 321 || abs(dk2nu->decay.ptype) == 310);
    double nupar_pz = dk2nu->decay.pppz;
    double nupar_pt = (dk2nu->decay.pppz) * sqrt(pow(dk2nu->decay.ppdxdz,2)+pow(dk2nu->decay.ppdxdz,2));
    //
    hflux_tot[ihel]->Fill(nuenergy,fluxWGT);
    if(is_pinupar){
      hflux_pinupar[ihel]->Fill(nuenergy,fluxWGT);
      hflux_foc_pinupar[ihel][ifoc]->Fill(nuenergy,fluxWGT);
      hpzpt_foc_pinupar[ihel][ifoc]->Fill(nupar_pz,nupar_pt,fluxWGT);
    }
    if(is_kanupar){
      hflux_kanupar[ihel]->Fill(nuenergy,fluxWGT);
      hflux_foc_kanupar[ihel][ifoc]->Fill(nuenergy,fluxWGT);
      hpzpt_foc_kanupar[ihel][ifoc]->Fill(nupar_pz,nupar_pt,fluxWGT);
    }
  }
  
  fOut->cd();
  for(int i=0;i<Nhel;i++){
    fOut->mkdir(hel[i].c_str());
    fOut->cd(hel[i].c_str());
    hflux_tot[i]->Write();
    hflux_pinupar[i]->Write();
    hflux_kanupar[i]->Write();
    for(int j=0;j<Nfoc;j++){
      hflux_foc_pinupar[i][j]->Write();
      hflux_foc_kanupar[i][j]->Write();
      hpzpt_foc_pinupar[i][j]->Write();
      hpzpt_foc_kanupar[i][j]->Write();
    }
  }
  std::cout<<"ending()"<<std::endl;
  
}
int idx_hel(int pdg){
  int iout = -1;
  if(pdg == 14) iout = 0;
  if(pdg ==-14) iout = 1;
  if(pdg == 12) iout = 2;
  if(pdg ==-12) iout = 3;
  return iout;
}
int idx_foc(std::string foc_comp){
  int iout = -1;
  for(int ii=0;ii<Nfoc;ii++){
    if(foc_comp == foc[ii])iout = ii;
  }
  return iout;
}

////////////////////////////////
#ifndef __CINT__
int main(int argc, const char* argv[]){  
  doFocComponents(argv[1],argv[2],argv[3]);
  return 0;
}
#endif

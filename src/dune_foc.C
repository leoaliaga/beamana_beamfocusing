#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h> 

#include "trajana_tree.h"
//bins
const int Nebin    = 1200;
const int Npzbin   = 1200;
const int Nptbin   =  100;

const double emin  =   0.;
const double emax  = 120.;
const double pzmin =   0.;
const double pzmax = 120.;
const double ptmin =   0.;
const double ptmax =   1.;

const int Nhel      = 4;
const char* hel[Nhel] = {"numu","numubar","nue","nuebar"};

const int Nfoc        = 14;
const char* foc[Nfoc] = {"notgt",
			 "HABCunder","HABover-HBCunder","HABover-HBCover","HABunder-HBCover",
			 "NeckABC","NeckAB-HC","NeckA-HB-NeckC","NeckA-HBC",
			 "HA-NeckBC","HA-NeckB-HC","HA-HB-NeckC",
			 "NeckOther",
			 "no-classified"};

int idx_hel(int pdgdcode);
bool passByNeck(double neckpos,double neckR,double irr, double izz, double frr, double fzz);

class trajana_tree;

int get_category(trajana_tree& tevts);

void dune_foc(const char* fileIn, const char* cdet,const char* fileOut){
  
  //makes the focusing componets for DUNE
  std::cout<<"=> dune_foc()"         <<std::endl;
  std::cout<<"=> fileIn : " <<fileIn <<std::endl;
  std::cout<<"cdet      : " <<cdet   <<std::endl;
  std::cout<<"=> fileOut: " <<fileOut<<std::endl;
  
  int idet = atoi(cdet);
  if(idet<0 || idet>2)exit (1);
  
  TFile* fOut = new TFile(fileOut,"recreate");

  //Total:
  TH1D* hflux[Nhel];
  TH1D* hfluxpz[Nhel];
  TH1D* hfluxpt[Nhel];
  TH1D* hfluxptot[Nhel];
  TH2D* hpzpt[Nhel];
  TH2D* hpzenu[Nhel];
  TH1D* hflux_foc[Nhel][Nfoc];
  TH1D* hfluxpz_foc[Nhel][Nfoc];
  TH1D* hfluxpt_foc[Nhel][Nfoc];
  TH1D* hfluxptot_foc[Nhel][Nfoc];
  TH2D* hpzpt_foc[Nhel][Nfoc];
  TH2D* hpzenu_foc[Nhel][Nfoc];
  //Only for pip and pim neutrino parents
  TH1D* hflux_pi[Nhel];
  TH1D* hfluxpz_pi[Nhel];
  TH1D* hfluxpt_pi[Nhel];
  TH1D* hfluxptot_pi[Nhel];
  TH2D* hpzpt_pi[Nhel];
  TH2D* hpzenu_pi[Nhel];
  TH1D* hflux_foc_pi[Nhel][Nfoc];
  TH1D* hfluxpz_foc_pi[Nhel][Nfoc];
  TH1D* hfluxpt_foc_pi[Nhel][Nfoc];
  TH1D* hfluxptot_foc_pi[Nhel][Nfoc];
  TH2D* hpzpt_foc_pi[Nhel][Nfoc];
  TH2D* hpzenu_foc_pi[Nhel][Nfoc];
  //only for kap and kam neutrino parents
  TH1D* hflux_ka[Nhel];
  TH2D* hpzpt_ka[Nhel];
  TH1D* hflux_foc_ka[Nhel][Nfoc];
  TH2D* hpzpt_foc_ka[Nhel][Nfoc];
  //
  for(int i=0;i<Nhel;i++){
    hflux[i]    = new TH1D(Form("hflux_%s"   ,hel[i]),"",Nebin,emin,emax);
    hfluxpz[i]    = new TH1D(Form("hfluxpz_%s"   ,hel[i]),"",Npzbin,pzmin,pzmax);
    hfluxpt[i]    = new TH1D(Form("hfluxpt_%s"   ,hel[i]),"",Nptbin,ptmin,ptmax);
    hfluxptot[i]    = new TH1D(Form("hfluxptot_%s"   ,hel[i]),"",Nebin,emin,emax);
    hflux_pi[i] = new TH1D(Form("hflux_pi_%s",hel[i]),"",Nebin,emin,emax);
    hfluxpz_pi[i] = new TH1D(Form("hfluxpz_pi_%s",hel[i]),"",Nebin,emin,emax);
    hfluxpt_pi[i] = new TH1D(Form("hfluxpt_pi_%s",hel[i]),"",Nptbin,ptmin,ptmax);
    hfluxptot_pi[i] = new TH1D(Form("hfluxptot_pi_%s",hel[i]),"",Nebin,emin,emax);
    hflux_ka[i] = new TH1D(Form("hflux_ka_%s",hel[i]),"",Nebin,emin,emax);
    hpzpt[i]    = new TH2D(Form("hpzpt_%s"   ,hel[i]),"",Npzbin,pzmin,pzmax,Nptbin,ptmin,ptmax); 
    hpzenu[i]    = new TH2D(Form("hpzenu_%s"   ,hel[i]),"",Npzbin,pzmin,pzmax,Nebin,emin,emax); 
    hpzpt_pi[i] = new TH2D(Form("hpzpt_pi_%s",hel[i]),"",Npzbin,pzmin,pzmax,Nptbin,ptmin,ptmax);
    hpzenu_pi[i] = new TH2D(Form("hpzenu_pi_%s",hel[i]),"",Npzbin,pzmin,pzmax,Nebin,emin,emax);
    hpzpt_ka[i] = new TH2D(Form("hpzpt_ka_%s",hel[i]),"",Npzbin,pzmin,pzmax,Nptbin,ptmin,ptmax);
    //
    for(int j=0;j<Nfoc;j++){
      hflux_foc[i][j]    = new TH1D(Form("hflux_%s_%s"   ,hel[i],foc[j]),"",Nebin,emin,emax);
      hfluxpz_foc[i][j]    = new TH1D(Form("hfluxpz_%s_%s"   ,hel[i],foc[j]),"",Nebin,emin,emax);
      hfluxpt_foc[i][j]    = new TH1D(Form("hfluxpt_%s_%s"   ,hel[i],foc[j]),"",Nptbin,ptmin,ptmax);
      hfluxptot_foc[i][j]    = new TH1D(Form("hfluxptot_%s_%s"   ,hel[i],foc[j]),"",Nebin,emin,emax);
      hflux_foc_pi[i][j] = new TH1D(Form("hflux_pi_%s_%s",hel[i],foc[j]),"",Nebin,emin,emax);
      hfluxpz_foc_pi[i][j] = new TH1D(Form("hfluxpz_pi_%s_%s",hel[i],foc[j]),"",Nebin,emin,emax);
      hfluxpt_foc_pi[i][j] = new TH1D(Form("hfluxpt_pi_%s_%s",hel[i],foc[j]),"",Nptbin,ptmin,ptmax);
      hfluxptot_foc_pi[i][j] = new TH1D(Form("hfluxptot_pi_%s_%s",hel[i],foc[j]),"",Nebin,emin,emax);
      hflux_foc_ka[i][j] = new TH1D(Form("hflux_ka_%s_%s",hel[i],foc[j]),"",Nebin,emin,emax);
      //
      hpzpt_foc[i][j]    = new TH2D(Form("hpzpt_%s_%s"   ,hel[i],foc[j]),"",Npzbin,pzmin,pzmax,Nptbin,ptmin,ptmax); 
      hpzenu_foc[i][j]    = new TH2D(Form("hpzenu_%s_%s"   ,hel[i],foc[j]),"",Npzbin,pzmin,pzmax,Nebin,emin,emax);     
      hpzpt_foc_pi[i][j] = new TH2D(Form("hpzpt_pi_%s_%s",hel[i],foc[j]),"",Npzbin,pzmin,pzmax,Nptbin,ptmin,ptmax);
      hpzenu_foc_pi[i][j] = new TH2D(Form("hpzenu_pi_%s_%s",hel[i],foc[j]),"",Npzbin,pzmin,pzmax,Nebin,emin,emax);
      hpzpt_foc_ka[i][j] = new TH2D(Form("hpzpt_ka_%s_%s",hel[i],foc[j]),"",Npzbin,pzmin,pzmax,Nptbin,ptmin,ptmax);
    }
  }
  
  TChain* chain = new TChain("TrajAna");  
  std::ifstream ifs;
  ifs.open(fileIn);
  std::string line;
  int counter = 0;
  while (ifs.good()){
    getline(ifs,line);
    if(line.find(".root")>10000)continue;
    counter++; 
    //if(counter>1)continue;
    if(counter%50==0)std::cout<<counter<<" trees loaded"<<std::endl;
    chain->Add(line.c_str());
  }//end of list of files
  ifs.close();
  
  trajana_tree evts(chain);
  std::cout<<"=> ntrees  : "<< evts.ntrees   <<std::endl;
  std::cout<<"=> nentries: "<< evts.nentries <<std::endl;
  
  for(int i=0; i < evts.nentries; i++){
    //for(int i=0; i < 1000; i++){
    if(i%100000==0)std::cout<<i/1000<<" k entries"<<std::endl;
    
    evts.GetEntry(i);
    int ihel = idx_hel(evts.nupdg);
    int ifoc = get_category(evts);
    if(ihel<0 || ifoc<0){
      std::cout<<"something wrong, exiting: (ihel,ifoc) = ("<<ihel<<", "<<ifoc<<")"<<std::endl;
      exit (1);
    }
    double wgt = evts.impwgt * (evts.geowgt->at(idet));    
    double enu = evts.enu->at(idet);
    
    bool is_pi = ( abs(evts.parpdg) == 211 );
    bool is_ka = ( abs(evts.parpdg) == 321 );
    double pz  = evts.momz->at(0)/1000.;
    double pt  = sqrt(pow(evts.momx->at(0),2)+pow(evts.momy->at(0),2))/1000.;
    double ptot  = sqrt(pow(evts.momx->at(0),2)+pow(evts.momy->at(0),2)+pow(evts.momz->at(0),2))/1000.;
        
    hflux[ihel]->Fill(enu  ,wgt);
    hfluxpz[ihel]->Fill(pz  ,wgt);
    hfluxpt[ihel]->Fill(pt  ,wgt);
    hfluxptot[ihel]->Fill(ptot  ,wgt);
    hpzpt[ihel]->Fill(pz,pt,wgt);
    hpzenu[ihel]->Fill(pz,enu,wgt);
    if(is_pi)hflux_pi[ihel]->Fill(enu   ,wgt);
    if(is_pi)hfluxpz_pi[ihel]->Fill(pz   ,wgt);
    if(is_pi)hfluxpt_pi[ihel]->Fill(pt   ,wgt);
    if(is_pi)hfluxptot_pi[ihel]->Fill(ptot   ,wgt);
    if(is_pi)hpzpt_pi[ihel]->Fill(pz,pt ,wgt);
    if(is_pi)hpzenu_pi[ihel]->Fill(pz,enu ,wgt);
    if(is_ka)hflux_ka[ihel]->Fill(enu   ,wgt);
    if(is_ka)hpzpt_ka[ihel]->Fill(pz,pt ,wgt);
    //foc:
    if(ifoc>=0){
      hflux_foc[ihel][ifoc]->Fill(enu,wgt);
      hfluxpz_foc[ihel][ifoc]->Fill(pz,wgt);
      hfluxpt_foc[ihel][ifoc]->Fill(pt,wgt);
      hfluxptot_foc[ihel][ifoc]->Fill(ptot,wgt);
      hpzpt_foc[ihel][ifoc]->Fill(pz,pt,wgt);
      hpzenu_foc[ihel][ifoc]->Fill(pz,enu,wgt);
      if(is_pi)hflux_foc_pi[ihel][ifoc]->Fill(enu   ,wgt);
      if(is_pi)hfluxpz_foc_pi[ihel][ifoc]->Fill(pz   ,wgt);
      if(is_pi)hfluxpt_foc_pi[ihel][ifoc]->Fill(pt   ,wgt);
      if(is_pi)hfluxptot_foc_pi[ihel][ifoc]->Fill(ptot   ,wgt);
      if(is_pi)hpzpt_foc_pi[ihel][ifoc]->Fill(pz,pt ,wgt);
      if(is_pi)hpzenu_foc_pi[ihel][ifoc]->Fill(pz,enu ,wgt);
      if(is_ka)hflux_foc_ka[ihel][ifoc]->Fill(enu   ,wgt);
      if(is_ka)hpzpt_foc_ka[ihel][ifoc]->Fill(pz,pt ,wgt);
    }
  }

  //Storing
  fOut->cd();
  for(int i=0;i<Nhel;i++){
    fOut->mkdir(Form("%s",hel[i]));
    fOut->mkdir(Form("%s_pi",hel[i]));
    fOut->mkdir(Form("%s_ka",hel[i]));
  }    
  for(int i=0;i<Nhel;i++){
    fOut->cd(Form("%s",hel[i]));
    hflux[i]->Write();  
    hfluxpz[i]->Write(); 
    hfluxpt[i]->Write(); 
    hfluxptot[i]->Write();    
    hpzpt[i]->Write();
    hpzenu[i]->Write();    
    
    fOut->cd(Form("%s_pi",hel[i]));
    hflux_pi[i]->Write();
    hfluxpz_pi[i]->Write();
    hfluxpt_pi[i]->Write();
    hfluxptot_pi[i]->Write();     
    hpzpt_pi[i]->Write();
    hpzenu_pi[i]->Write(); 
    
    fOut->cd(Form("%s_ka",hel[i]));
    hflux_ka[i]->Write(); 
    hpzpt_ka[i]->Write(); 
    //
    for(int j=0;j<Nfoc;j++){
      fOut->cd(Form("%s",hel[i]));
      hflux_foc[i][j]->Write();
      hfluxpz_foc[i][j]->Write(); 
      hfluxpt_foc[i][j]->Write();
      hfluxptot_foc[i][j]->Write();   
      hpzpt_foc[i][j]->Write();
      hpzenu_foc[i][j]->Write();   
      
      fOut->cd(Form("%s_pi",hel[i]));
      hflux_foc_pi[i][j]->Write();
      hfluxpz_foc_pi[i][j]->Write();
      hfluxpt_foc_pi[i][j]->Write();
      hfluxptot_foc_pi[i][j]->Write();
      hpzpt_foc_pi[i][j]->Write();
      hpzenu_foc_pi[i][j]->Write();
      
      fOut->cd(Form("%s_ka",hel[i]));
      hflux_foc_ka[i][j]->Write();
      hpzpt_foc_ka[i][j]->Write();
    }
  }
  fOut->Close();
  std::cout<<"=> dune_foc()"<<std::endl;
}
//
int get_category(trajana_tree& tevts){
  int idx   = -1;
  int ntraj = tevts.prevol->size();

  //First check if the nu parent was born in the target
  bool born_in_tgt = ( (tevts.prevol->at(0)).find("TargetSimpleCylinder")==0 || 
		       (tevts.prevol->at(0)).find("TargetNoSplit")==0 );  
  if(!born_in_tgt){idx=0;return idx;}
  
  //i: enter to horn-i
  double defN = 99999999;
  //indices: [0, 1, 2] -> [A, B, C] 
  double enterHx[3] = {-1.*defN,-1.*defN,-1.*defN};
  double enterHy[3] = {-1.*defN,-1.*defN,-1.*defN};
  double enterHz[3] = {-1.*defN,-1.*defN,-1.*defN};
  bool   enterHH[3] = {false,false,false};
  bool   inHH[3]    = {false,false,false};
  bool   byNeck[3]  = {false,false,false};
  
  //Check if the nu parent entered or was inside the target:
  for(int ii=0;ii<ntraj;ii++){
    //
    bool enterHA = (tevts.prevol->at(ii)).find("TargetHallAndHorn1_P")==0 &&
      (tevts.postvol->at(ii)).find("Horn1PolyM1_P")==0;
    bool enterHB = (tevts.prevol->at(ii)).find("Tunnel_P")==0 &&
      (tevts.postvol->at(ii)).find("LBNFSimpleHorn2Container_P")==0;
    bool enterHC = (tevts.prevol->at(ii)).find("Tunnel_P")==0 &&
      (tevts.postvol->at(ii)).find("LBNFSimpleHorn3Container_P")==0;  
    if(enterHA){
      enterHx[0] = tevts.posx->at(ii);
      enterHy[0] = tevts.posy->at(ii);
      enterHz[0] = tevts.posz->at(ii);  
      enterHH[0] = true;
    }
    if(enterHB){
      enterHx[1] = tevts.posx->at(ii);
      enterHy[1] = tevts.posy->at(ii);
      enterHz[1] = tevts.posz->at(ii);   
      enterHH[1] = true;
    }
    if(enterHC){
      enterHx[2] = tevts.posx->at(ii);
      enterHy[2] = tevts.posy->at(ii);
      enterHz[2] = tevts.posz->at(ii);   
      enterHH[2] = true;
    }
    // 
    bool insideHA  = (tevts.prevol->at(ii)).find("Horn1PolyM1_P")==0 ||
      (tevts.prevol->at(ii)).find("LBNFConceptHornA")==0;
    bool insideHB  = (tevts.prevol->at(ii)).find("LBNFSimpleHorn2Container_P")==0 ||
      (tevts.prevol->at(ii)).find("LBNFConceptHornB")==0;
    bool insideHC  = (tevts.prevol->at(ii)).find("LBNFSimpleHorn3Container_P")==0 ||
      ( (tevts.prevol->at(ii)).find("LBNFConceptHornC")==0 && (tevts.prevol->at(ii)).find("LBNFConceptHornCStrpL")!=0);  
    if(insideHA) inHH[0] = true;
    if(insideHB) inHH[1] = true;
    if(insideHC) inHH[2] = true;  
  }
  
  //Per trajectory
  for(int ii=0;ii<(ntraj-1);ii++){  
    //Is still in target or has not born in the target:
    bool is_tgt       = (tevts.prevol->at(ii)).find("TargetNoSplit")==0 ||
      (tevts.prevol->at(ii)).find("TargetSimpleCylinder")==0;
    double posZ       = tevts.posz->at(ii);    
    if(is_tgt || (posZ < 0))continue;
    
    double posZp1  = tevts.posz->at(ii+1);    
    double disR    = sqrt( pow(tevts.posx->at(ii),2) + pow(tevts.posy->at(ii),2) ); 
    
    //Neck A
    if( !enterHH[0] && !inHH[0] && (enterHx[0]>(defN-1) || enterHx[0]<(-1.*defN+1)) ){
      if(posZ<2254 && disR<50 && posZp1>-52){
	enterHx[0] = defN;
	enterHy[0] = defN;
	enterHz[0] = defN;
	byNeck[0]  = true;      
      }	
    }
    //Neck B:
    if( !enterHH[1] && !inHH[1] && (enterHx[1]>(defN-1) || enterHx[1]<(-1.*defN+1)) ){
      if(posZ<7793 && posZp1>2835){
	enterHx[1] = defN;
	enterHy[1] = defN;
	enterHz[1] = defN;
	double disnextR = sqrt( pow(tevts.posx->at(ii+1),2) + pow(tevts.posy->at(ii+1),2) ); 
	double posnextZ = tevts.posz->at(ii+1);
	if(passByNeck(5100,85,disR,posZ,disnextR,posnextZ))byNeck[1] = true;
      }
    }
    
    //Neck C:
    if( !enterHH[2] && !inHH[2] && (enterHx[2]>(defN-1) || enterHx[2]<(-1.*defN+1)) ){
      if(posZ<19737 && posZp1>16234){
	enterHx[2] = defN;
	enterHy[2] = defN;
	enterHz[2] = defN;
	double disnextR = sqrt( pow(tevts.posx->at(ii+1),2) + pow(tevts.posy->at(ii+1),2) ); 
	double posnextZ = tevts.posz->at(ii+1);
	if(passByNeck(18200,140,disR,posZ,disnextR,posnextZ))byNeck[2] = true;
      }
    }
  }//end necks traj
  
  //Looking for under/over
  bool HABunder = false;
  bool HABover  = false;  
  bool HBCunder = false;
  bool HBCover  = false;  
  if(enterHH[0] && enterHH[1]){
    double mvalAB = -1.0*enterHx[0]/enterHy[0];    
    double auxHA  = enterHy[0] - mvalAB * enterHx[0];
    double auxHB  = enterHy[1] - mvalAB * enterHx[1];
    bool no_cross = (auxHA>0. && auxHB>0.) || (auxHA<0. && auxHB<0.);
    if(no_cross)HABunder = true;
    else	HABover  = true;
  }
  if(enterHH[1] && enterHH[2]){
    double mvalBC = -1.0*enterHx[1]/enterHy[1];    
    double auxHB  = enterHy[1] - mvalBC * enterHx[1];
    double auxHC  = enterHy[2] - mvalBC * enterHx[2];
    bool no_cross = (auxHB>0. && auxHC>0.) || (auxHB<0. && auxHC<0.);
    if(no_cross)HBCunder = true;
    else	HBCover  = true;
  }

  if(     HABunder   && HBCunder)                 idx = 1;
  else if(HABover    && HBCunder)                 idx = 2;
  else if(HABover    && HBCover )                 idx = 3;
  else if(HABunder   && HBCover )                 idx = 4;
  
  else if(byNeck[0]  && byNeck[1]  && byNeck[2])  idx = 5;
  else if(byNeck[0]  && byNeck[1]  && enterHH[2]) idx = 6;
  else if(byNeck[0]  && enterHH[1] && byNeck[2])  idx = 7;
  else if(byNeck[0]  && enterHH[1] && enterHH[2]) idx = 8;

  else if(enterHH[0] && byNeck[1]  && byNeck[2])  idx = 9;
  else if(enterHH[0] && byNeck[1]  && enterHH[2]) idx = 10;
  else if(enterHH[0] && enterHH[1]  && byNeck[2]) idx = 11;
    
  else if(byNeck[0] || byNeck[1] || byNeck[2])   idx = 12;
  
  else idx = 13;
  
  return idx;
}
//
bool passByNeck(double neckpos,double neckR,double irr, double izz, double frr, double fzz){
  bool pass = false;
  if(izz==fzz)return pass;
  double mm  = (frr-irr)/(fzz-izz);
  double val = frr + mm * (neckpos-fzz); 
  if(val>0 && val<neckR && izz<neckpos && neckpos<fzz)pass = true;
  return pass;
}
int idx_hel(int pdgcode){
  int idx = -1;
  if(pdgcode ==  14)idx = 0;
  if(pdgcode == -14)idx = 1;
  if(pdgcode ==  12)idx = 2;
  if(pdgcode == -12)idx = 3;
  return idx;
}
//
////////////////////////////////
#ifndef __CINT__
int main(int argc, const char* argv[]){  
  dune_foc(argv[1],argv[2],argv[3]);
  return 0;
}
#endif



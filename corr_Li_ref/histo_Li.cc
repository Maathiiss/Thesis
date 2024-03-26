#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TApplication.h"
#include "TMultiGraph.h"
#include "TFeldmanCousins.h"
#include "TGaxis.h"
#include "TLeaf.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TParameter.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TSystem.h>
#include "/home/granjon/Documents/stage_radon/Stage/Mathis/sndisplay/sndisplay.cc"


using namespace std;


string namer(int om){    /////// name for the ref om for legend
  string name;
  if (om == 712) {name = "Ref MW1";}
  if (om == 713) {name = "Ref MW2";}
  if (om == 714) {name = "Ref GV";}
  if (om == 715) {name = "Ref XW1";}
  if (om == 716) {name = "Ref XW2";}
  return name;
}

int intensity_chooser(int pic){  ///// to obtain intensity
  //associer intensite LED a un numero pic
  int intensity;
  if (pic == 1) {
    intensity = 70;
  }
  if (pic == 2) {
    intensity = 80;
  }
  if (pic == 3) {
    intensity = 90;
  }
  if (pic == 4) {
    intensity = 100;
  }
  if (pic == 5) {
    intensity = 110;
  }
  if (pic == 6) {
    intensity = 120;
  }
  return intensity;
}

std::vector<double> time_measurer(int run_number){
  TFile tree_file(Form("entree/root/snemo_run-%d_LI.root", run_number), "READ");
  int om_number;
  double time;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  gROOT->cd();
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);

  TH1D *time_spectre = new TH1D ("time_spectre", "", 1000, 0, 550);
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    time_spectre->Fill(time);
  }
  // tree->Project("time_time_spectre", "time");
  std::cout << Form("entree/root/snemo_run-%d_LI.root", run_number) << '\n';
  // time_spectre->Draw();

  int plus =0;
  int moins = 0;
  std::vector<double> time_measurer;
  int i;
  for (i = 0; i < 1000; i++) {
    if (time_spectre->GetBinContent(i) > 0) {
      plus = 1;
    }
    if (time_spectre->GetBinContent(i) == 0) {
      moins = 1;
    }
    if (plus - moins == 0) {
      time_measurer.push_back(time_spectre->GetBinCenter(i));
      plus = 0;
      moins = 0;
      // std::cout << "time = " << time_spectre->GetBinCenter(i) << '\n';
    }
  }
  std::vector<double> interval;
  interval.push_back(time_measurer.at(0));
  for (size_t j = 1; j < time_measurer.size()-2; j+=2) {
    interval.push_back((time_measurer.at(j)+time_measurer.at(j+1))/2);
  }
  interval.push_back(time_measurer.at(time_measurer.size()-1));
  delete time_spectre;
  return interval;
}



int pic_selector(TH1D* spectre) {//selection pic le plus proche d'1 MeV -> a refaire
  int selector = 0;
  if (spectre->GetMean() > 150 && spectre->GetMean() < 350)selector = 1;
  return selector;
}


//Ref_correcteor = ouvre fichier ref et recup dans tableau la valeur correction pour chaque OM (ici ref -> 5 seulement)
double* Ref_corrector(int run, string correction, double *gain_tab, double *gain_tab_error) {
  //corriger la variation du gain des OMs de reference
  //correction = nom fichier = ensemble run ref -> boucle dessus
  //int run = pour trouver bon run
  
  //TFile file1(Form("/home/granjon/Documents/these/Analyse/corr_OM_ref/sortie/root_file_final/Fit_Ref_%s.root", correction.c_str()), "READ");   //mini ou normal ?
  TFile file1(Form("/home/granjon/Documents/these/Analyse/corr_OM_ref/sortie/root_file_final/Fit_Ref_%s.root", correction.c_str()), "READ");
  std::cout<<"FILEEEE "<<Form("/home/granjon/Documents/these/Analyse/corr_OM_ref/sortie/root_file_final/Fit_Ref_%s.root", correction.c_str())<<std::endl;

  // std::cout << Form("calcul_gain_ref/root/Merged_Fit/Fit_Ref_%s.root", correction.c_str()) << '\n';

  int run_number;
  double gain, gain_error;
  TTree* tree = (TTree*)file1.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);
  tree->SetBranchStatus("gain_error_plus",1);
  tree->SetBranchAddress("gain_error_plus", &gain_error);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);

  int compteur = 0;
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);    
    //mettre abs ici ???
    if (run_number == run+1) {//pour associer run correction avec LI
      //std::cout<<" RUN FIT REF  "<<run_number<<" RUN LI  "<< run<<"\n";
      gain_tab[compteur] = gain;
      gain_tab_error[compteur] = gain_error;
      compteur++;
    }
  }
  file1.Close();
  return gain_tab;
}



int Ref_bundle_number(int om_number, string wall){ // associer un num ref a un bundle = partie du detecteur
  int bundle_number = 0;
  if (om_number == 712 && wall.compare("IT") == 0) {
    bundle_number = 1;
  }
  if (om_number == 713 && wall.compare("IT") == 0) {
    bundle_number = 3;
  }
  if (om_number == 714 && wall.compare("IT") == 0) {
    bundle_number = 4;
  }
  if (om_number == 715 && wall.compare("IT") == 0) {
    bundle_number = 2;
  }
  if (om_number == 716 && wall.compare("IT") == 0) {
    bundle_number = 5;
  }
  if (om_number == 712 && wall.compare("FR") == 0) {
    bundle_number = 6;
  }
  if (om_number == 713 && wall.compare("FR") == 0) {
    bundle_number = 7;
  }
  if (om_number == 714 && wall.compare("FR") == 0) {
    bundle_number = 9;
  }
  if (om_number == 715 && wall.compare("FR") == 0) {
    bundle_number = 8;
  }
  if (om_number == 716 && wall.compare("FR") == 0) {
    bundle_number = 10;
  }
  return bundle_number;
}


int bundle_number(int om_number){//associer un numero d'OM a un bundle 
  int bundle_number = 0;
  if ((om_number%13 > 7 && om_number/13 < 10 && om_number < 260) || (om_number > 663 && om_number < 672) || (om_number < 552 && om_number > 545) || (om_number < 536 && om_number > 529)) {
    bundle_number = 1;
  }
  else if ((om_number%13 < 8 && om_number/13 < 6 && om_number < 260) ||(om_number > 647 && om_number < 652) ||(om_number < 546 && om_number > 535) ||(om_number < 530 && om_number > 519)){
    bundle_number = 2;
  }
  else if ((om_number%13 > 7 && om_number/13 > 9 && om_number < 260) || (om_number > 671 && om_number < 680) || (om_number < 568 && om_number > 561) || (om_number < 584 && om_number > 577)) {
    bundle_number = 3;
  }
  else if ((om_number%13 < 8 && om_number/13 > 5 && om_number/13 < 14 && om_number < 260) || (om_number < 660 && om_number > 651)) {
    bundle_number = 4;
  }
  else if ((om_number%13 < 8 && om_number/13 > 13 && om_number < 260) || (om_number < 664 && om_number > 659) || (om_number < 562 && om_number > 551) ||(om_number < 578 && om_number > 567)) {
    bundle_number = 5;
  }
  else if ((om_number%13 > 7 && (om_number/13-20) < 10 && om_number < 520 && om_number > 259) ||(om_number > 695 && om_number < 704) ||(om_number < 600 && om_number > 593) ||(om_number < 616 && om_number > 609)){
    bundle_number = 6;
  }
  else if ((om_number%13 > 7 && (om_number/13-20) > 9 && om_number < 520 && om_number > 259) ||(om_number > 703 && om_number < 712) ||(om_number < 648 && om_number > 641) ||(om_number < 632 && om_number > 625)){
    bundle_number = 7;
  }
  else if ((om_number%13 < 8 && (om_number/13-20) < 6 && om_number < 520 && om_number > 259) ||(om_number > 679 && om_number < 684) ||(om_number < 594 && om_number > 583) ||(om_number < 610 && om_number > 599)){
    bundle_number = 8;
  }
  else if ((om_number%13 < 8 && (om_number/13-20) > 5 && (om_number/13-20) < 14 && om_number < 520 && om_number > 259) ||(om_number > 683 && om_number < 692)){
    bundle_number = 9;
  }
  else if ((om_number%13 < 8 && (om_number/13-20) > 13 && om_number < 520 && om_number > 259) ||(om_number > 691 && om_number < 696) ||(om_number < 642 && om_number > 631) ||(om_number < 626 && om_number > 615)){
    bundle_number = 10;
  }
  return bundle_number;
}



void Li_corrector(std::vector<int> run, int n_run, int start, int stop){
  double Amplitude, Amplitude_error, Amplitude_corr, Amplitude_corr_error, Khi2, time, gain, gain_error;
  int om_number, run_number, pic,bundle_ref;
  string wall;
  TFile *file_sortie = new TFile(Form("/home/granjon/Li/sortie/ref_Li/Fit_Ampl_Ref/final_gain_%d-%d.root",start,stop),"RECREATE");
  std::cout<<"file create "<<Form("/home/granjon/Li/sortie/ref_Li/Fit_Ampl_Ref/final_gain_%d-%d.root",start,stop)<<std::endl;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Amplitude", &Amplitude);
  Result_tree.Branch("Amplitude_error", &Amplitude_error);
  // Result_tree.Branch("Amplitude_corr", &Amplitude_corr);
  // Result_tree.Branch("Amplitude_corr_error", &Amplitude_corr_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("gain",&gain);
  Result_tree.Branch("gain_error",&gain_error);
  Result_tree.Branch("bundle_ref",&bundle_ref);
  
  double amp_scor[5][8][n_run];
  double amp_scor_error[5][8][n_run];
  double amp_cor[5][8][n_run];
  double amp_cor_error[5][8][n_run];
  int bundle_ref_vec[5][8][n_run];

  TFile *file = new TFile(Form("/home/granjon/Li/sortie/ref_Li/Fit_Ampl_Ref/Merged_Fit_%d-%d.root",start,stop), "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("pic",1);
  tree->SetBranchAddress("pic", &pic);
  tree->SetBranchStatus("Khi2",1);
  tree->SetBranchAddress("Khi2", &Khi2);
  tree->SetBranchStatus("Amplitude",1);
  tree->SetBranchAddress("Amplitude", &Amplitude);
  tree->SetBranchStatus("Amplitude_error",1);
  tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree->SetBranchStatus("Amplitude_corr",1);
  tree->SetBranchAddress("Amplitude_corr", &Amplitude_corr);
  tree->SetBranchStatus("Amplitude_corr_error",1);
  tree->SetBranchAddress("Amplitude_corr_error", &Amplitude_corr_error);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchStatus("bundle_ref",1);
  tree->SetBranchAddress("bundle_ref", &bundle_ref);

  int scompteur = 0;
  int srun = start;
   for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (run_number != srun) {
      scompteur++;
      srun = run_number;
    }
    if(pic<5){
      wall = "IT";
    }
    else if(pic>=5){
      wall = "FR";
    }
    amp_cor[om_number -712][pic-1][scompteur] = Amplitude_corr;
    amp_cor_error[om_number -712][pic-1][scompteur] = Amplitude_corr_error;
    amp_scor[om_number -712][pic-1][scompteur] = Amplitude;
    amp_scor_error[om_number -712][pic-1][scompteur] = Amplitude_error;
    bundle_ref_vec[om_number -712][pic-1][scompteur] = bundle_ref;

  }

   for (int l = 0; l < n_run; l++) {//run num
    for (int i = 0; i < 5; i++) {//om num    
      for (int j = 0; j < 8; j++) {//pic num      
          gain =  amp_cor[i][j][l]/amp_cor[i][j][0];	  
	  gain_error = gain * sqrt((amp_cor_error[i][j][l]/amp_cor[i][j][l])*(amp_cor_error[i][j][l]/amp_cor[i][j][l]) + (amp_cor_error[i][j][0]/amp_cor[i][j][0])*(amp_cor_error[i][j][0]/amp_cor[i][j][0]));
	  if(l==0){gain_error = 0;}
          om_number = i+712;
          pic = j+1;
	  if(j>3){pic =j-3;}
          run_number = run[l];
          Amplitude = amp_scor[i][j][l];
	  Amplitude_error = amp_scor_error[i][j][l];
          Amplitude_corr = amp_cor[i][j][l];
	  Amplitude_corr_error = amp_cor_error[i][j][l];
	  bundle_ref = bundle_ref_vec[i][j][l];
          Result_tree.Fill();        
      }
    }
  }
  file_sortie->cd();
  Result_tree.Write();
  file_sortie->Close();
}





void applied_Li_correction(std::vector<int> run, int n_run, int start, int stop){

  int run_number_ref, pic_ref, bundle_ref;
  double Amplitude, Amplitude_error, Amplitude_corr, Amplitude_corr_error, Khi2, time, gain, gain_error;
  int om_number, run_number, pic, bundle;

  TFile *file_sortie = new TFile(Form("/home/granjon/Li/sortie/ref_Li/Fit_Ampl_Ref/final_final_gain_%d-%d.root",start,stop),"RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Amplitude", &Amplitude);
  Result_tree.Branch("Amplitude_error", &Amplitude_error);
  Result_tree.Branch("Amplitude_corr", &Amplitude_corr);
  Result_tree.Branch("Amplitude_corr_error", &Amplitude_corr_error);  
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("gain",&gain);
  Result_tree.Branch("gain_error",&gain_error);
  Result_tree.Branch("bundle",&bundle);
  Result_tree.Branch("time",&time);


  
  //READ REF
 TFile *file = new TFile(Form("/home/granjon/Li/sortie/ref_Li/Fit_Ampl_Ref/final_gain_%d-%d.root",start,stop), "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("pic",1);
  tree->SetBranchAddress("pic", &pic_ref);
  tree->SetBranchStatus("bundle_ref",1);
  tree->SetBranchAddress("bundle_ref", &bundle_ref);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number_ref);
  tree->SetBranchStatus("gain_error",1);
  tree->SetBranchAddress("gain_error", &gain_error);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);


  TFile file_cor(Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/Merged_Fit_%d-%d.root",start,stop), "READ");
  TTree* tree_cor = (TTree*)file_cor.Get("Result_tree");
  tree_cor->SetBranchStatus("*",0);
  tree_cor->SetBranchStatus("Amplitude",1);
  tree_cor->SetBranchAddress("Amplitude", &Amplitude);
  tree_cor->SetBranchStatus("Amplitude_error",1);
  tree_cor->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree_cor->SetBranchStatus("Khi2",1);
  tree_cor->SetBranchAddress("Khi2", &Khi2);
  tree_cor->SetBranchStatus("run_number",1);
  tree_cor->SetBranchAddress("run_number", &run_number);
  tree_cor->SetBranchStatus("om_number",1);
  tree_cor->SetBranchAddress("om_number", &om_number);
  tree_cor->SetBranchStatus("pic",1);
  tree_cor->SetBranchAddress("pic", &pic);
  tree_cor->SetBranchStatus("time",1);
  tree_cor->SetBranchAddress("time", &time);

  tree_cor->SetBranchStatus("bundle",1);
  tree_cor->SetBranchAddress("bundle", &bundle);

  for (int j = 0; j < tree_cor->GetEntries(); j++) {
    tree_cor->GetEntry(j);    
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      //std::cout<<"i = "<< i << " j = "<<j<<std::endl;
      //std::cout<<"bundle " << bundle << " bundle ref "<<bundle_ref<<std::endl;
      if(pic_ref == pic && run_number == run_number_ref && bundle_ref == bundle){
	// std::cout<<"run FIND" << run_number << " run ref "<<run_number_ref<<std::endl;
	// std::cout<<"bundle FIND" << bundle << " bundle ref "<<bundle_ref<<std::endl;
	// std::cout<<"pic FIND" << pic << " pic ref "<<pic_ref<<std::endl;
	Amplitude_corr = Amplitude / gain;
	Amplitude_corr_error =  Amplitude_corr * sqrt(pow(Amplitude_error/Amplitude,2) + pow(gain_error/gain,2));
	Result_tree.Fill();
      }
    }
  }
  file_sortie->cd();
  Result_tree.Write();
  file_sortie->Close();
}




void fit_LI_amplitude(int run_number, double *correction_gain_table){ //pour tout le calo
  std::array<std::array<TH1D*,4>, 712> histograms; //712 OM et 4 pics
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  string wall=""; 
  TFile file(Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/Amplitude_Li_run_%d.root", run_number),"RECREATE");
  std::cout<<"file open "<<Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/Amplitude_Li_run_%d.root", run_number)<<std::endl;
  
  
  std::vector<double> interval;
  interval = time_measurer(run_number);
  for (size_t i = 0; i < interval.size(); i++) {
    std::cout << "interval[" << i << "] = " << interval.at(i) << '\n';
  }
  
  
  
  int om_number=-1,bundle;
  double constante;
  double mean, time, start_time;
  double mean_error, nevent;
  double sigma;
  int pic =0;
  double Khi2 = 0;
  int intensity = 0;
  int om;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("constante", &constante);
  Result_tree.Branch("Amplitude", &mean);
  Result_tree.Branch("Amplitude_error", &mean_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("intensity", &intensity);
  Result_tree.Branch("time", &start_time);
  Result_tree.Branch("nevent", &nevent);
  Result_tree.Branch("bundle", &bundle);

  double entries, saturation, under_threshold;

  TTree fail_tree("fail_tree","");
  fail_tree.Branch("entries", &entries);
  fail_tree.Branch("saturation", &saturation);
  fail_tree.Branch("under_threshold", &under_threshold);
  fail_tree.Branch("om_number", &om_number);
  fail_tree.Branch("pic", &pic);


  TCanvas* canvas = new TCanvas;
  TFile *tree_file = new TFile (Form("entree/root/snemo_run-%d_LI.root", run_number), "READ");

  double charge_tree;
  double amplitude_tree;
  TTree* tree = (TTree*)tree_file->Get("Result_tree");
  gROOT->cd();
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge_tree);
  tree->SetBranchStatus("amplitude_tree",1);
  tree->SetBranchAddress("amplitude_tree", &amplitude_tree);
  
  TParameter<double> *param = new TParameter<double>("start_time", time);
  param = (TParameter<double>*)(tree_file->Get("start_time"));
  start_time = param->GetVal();
  gROOT->cd();
  bool with_palette = true;
  //sncalo = new sndisplay::calorimeter ("sndiplay_dis",with_palette);
  for(int i=0;i<712;i++){
    //std::cout<<"interval size"<<interval.size()<<std::endl;
    for (size_t j = 0; j <interval.size()/2; j++){
      //std::cout<<"j"<<j<<std::endl;
      histograms[i][j] = new TH1D(Form("om_%d_pic_%lu",i,j+1),Form("om_%d_pic_%lu",i,j+1),1000, 0, 180000);
      //std::cout<<histograms[i][j]<<std::endl;
    }
  }
  for(int entry = 0; entry<tree->GetEntries(); entry++){
    int debut = 0;
    wall="IT";
    tree->GetEntry(entry);
    if(om_number<712 && om_number>=0){
      if ((om_number > 259 && om_number < 520) || (om_number > 583 && om_number < 647) || (om_number > 679 && om_number < 712) ) {
	debut = debut+interval.size()/2;
	wall="FR";
      }
      
      for (size_t j = debut; j < debut + (interval.size()/2) ; j++){
	if(time>interval.at(j) && time<interval.at(j+1) && amplitude_tree > 10){
	  if(wall == "FR"){	    
	    histograms[om_number][j-4]->Fill(charge_tree);
	  }
	  else{
	    histograms[om_number][j]->Fill(charge_tree);		      
	  }
	}
      }
    }
  }
  
  
  for(int om_number=0;om_number<712;om_number++){
    for (size_t j = 0; j <  interval.size()/2 ; j++){
      //file.cd();
      //histograms[om_number][j]->Write();
      //std::cout<<"histograms"<<om_number<<" j = "<<j<<std::endl;
      
      if (histograms[om_number][j]->GetEntries() < 300) {
	
	std::cout << "" << '\n';
	std::cout << "trop peu d'entries pour l'OM " << om_number << '\n';
	std::cout << "" << '\n';
	
	Khi2 = 0;
	mean = 0;
	mean_error = 0;
	entries = histograms[om_number][j]->GetEntries();
	saturation = 0;
	under_threshold = 0;
	pic = j+1;
	om=om_number;

	Result_tree.Fill();
	fail_tree.Fill();
	gSystem->mkdir(Form("sortie/SN_Li/fit_Li/run_%d",run_number));
	histograms[om_number][j]->Draw();
	//canvas->SaveAs(Form("sortie/SN_Li/fit_Li/run_%d/OM_%03d_pic_%d.png", run_number, om_number, pic));
	histograms[om_number][j]->Delete();
      }
      else if (histograms[om_number][j]->GetMean() > 180000) {
	std::cout << "" << '\n';
	std::cout << "the amplitude sature" << '\n';
	std::cout << "" << '\n';
	Khi2 = 0;
	mean = 0;
	mean_error = 0;
	entries = 0;
	saturation = histograms[om_number][j]->GetMean();
	under_threshold = 0;
	pic = j+1;
	om=om_number;
	bundle = bundle_number(om_number);
	Result_tree.Fill();
	fail_tree.Fill();
	gSystem->mkdir(Form("sortie/SN_Li/fit_Li/run_%d",run_number));
	histograms[om_number][j]->Draw();
	//canvas->SaveAs(Form("sortie/SN_Li/fit_Li/run_%d/OM_%03d_pic_%d.png", run_number, om_number,pic));
	histograms[om_number][j]->Delete();
      }
      else if (histograms[om_number][j]->GetMean() < -10){
	std::cout << "" << '\n';
	std::cout << "too few charge" << '\n';
	std::cout << "" << '\n';          
	Khi2 = 0;
	mean = 0;
	mean_error = 0;
	entries = 0;
	saturation = 0;
	under_threshold = histograms[om_number][j]->GetMean();
	pic = j+1;
	om=om_number;
	bundle = bundle_number(om_number);
	Result_tree.Fill();
	fail_tree.Fill(); //quelles sont les pics et le LI qui deconnent
	gSystem->mkdir(Form("sortie/SN_Li/fit_Li/run_%d",run_number));
	histograms[om_number][j]->Draw();
	//canvas->SaveAs(Form("sortie/SN_Li/fit_Li/run_%d/OM_%03d_pic_%d.png", run_number, om_number, pic));
	histograms[om_number][j]->Delete();
      }
      else{
	pic = j+1;
	nevent = histograms[om_number][j]->Integral();
	TF1 *f_Gaus = new TF1("f_Gaus", "gaus(0)", 0, 1200);
	f_Gaus->SetParNames("N_evt","mean_charge","Sigma");
	//std::cout<<"mean"<<histograms[om_number][j]->GetMean()<<"RMS"<<histograms[om_number][j]->GetRMS()<<std::endl;
	f_Gaus->SetParameters(25, histograms[om_number][j]->GetMean(), histograms[om_number][j]->GetRMS());
	f_Gaus->SetRange(histograms[om_number][j]->GetMean()-400, histograms[om_number][j]->GetMean()+400);
	f_Gaus->Draw("same");
	histograms[om_number][j]->Fit(f_Gaus, "RQ0");
	f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
	histograms[om_number][j]->Fit(f_Gaus, "RQ0");
	f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
	histograms[om_number][j]->Fit(f_Gaus, "RQ0");
	
	Khi2 = f_Gaus->GetChisquare()/f_Gaus->GetNDF();
	constante = (f_Gaus->GetParameter(0));
	mean = (f_Gaus->GetParameter(1));
	sigma = (f_Gaus->GetParameter(2));
	mean_error = f_Gaus->GetParError(1);
	
	
	bundle = bundle_number(om_number);
	//intensity = intensity_chooser(pic);
	histograms[om_number][j]->GetXaxis()->SetRangeUser(mean-10*sigma,mean+10*sigma);
	histograms[om_number][j]->Draw();
	f_Gaus->Draw("same");
	//gSystem->mkdir(Form("sortie/SN_Li/fit_Li/run_%d",run_number));
	//canvas->SaveAs(Form("sortie/SN_Li/fit_Li/run_%d/OM_%03d_pic_%d.png", run_number, om_number, pic));
	/*
	if(pic ==1){
	  sncalo->setcontent(om_number,mean);
	}
	*/
	om=om_number;
	Result_tree.Fill();
	histograms[om_number][j]->Delete();
	delete f_Gaus;
      }
    }
  }
  

  file.cd();
  Result_tree.Write();
  fail_tree.Write();
  file.Close();
  /*
  sncalo->draw();
  sncalo->draw();
  sncalo->canvas_it->SaveAs("sortie/SN_Li/graph_interessant/calo_pic_1_display_it.png");
  sncalo->canvas_fr->SaveAs("sortie/SN_Li/graph_interessant/calo_pic_1_display_fr.png");
  */
  return;
  
}






void fit_LI_amplitude_Ref(int run_number, double *ref_gain_table, double *ref_gain_table_error,string corr = ""){
  //Fit amplitude LI des OM de ref corrige de leurs gains ->fit en charge
  
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  

  TFile file(Form("/home/granjon/Li/sortie/ref_Li/Fit_Ampl_Ref/Amplitude_Li_run_%d.root", run_number),"RECREATE");

  
  std::vector<double> interval;
  interval = time_measurer(run_number); // sert a conserver que les pics qui t'interesse -> plus de 4 pics car utilisation de led secondaires
  // plutot que de faire 4 fits -> on coupe en temps en fonction de la duree du fichier -> te donne les bornes a partir de quand on est dans premier pic, celui d'apres ... 
  
  for (size_t i = 0; i < interval.size(); i++) {
    std::cout << "interval[" << i << "] = " << interval.at(i) << '\n';
  }
  // return;
  
  int om_number,bundle_ref;
  double constante;
  double mean, time, start_time, mean_corr;
  double mean_error, nevent, ref_error, mean_corr_error;
  double sigma;
  int pic =0;
  double Khi2 = 0;
  int intensity = 0;
  string wall = "IT";
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("constante", &constante);
  Result_tree.Branch("Amplitude", &mean);
  Result_tree.Branch("Amplitude_error", &mean_error);
  Result_tree.Branch("Ref_error", &ref_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("intensity", &intensity);
  Result_tree.Branch("time", &start_time);
  Result_tree.Branch("nevent", &nevent);
  Result_tree.Branch("wall", &wall);
  Result_tree.Branch("Amplitude_corr", &mean_corr);
  Result_tree.Branch("Amplitude_corr_error", &mean_corr_error);
  Result_tree.Branch("bundle_ref", &bundle_ref);
  
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
;
  TFile *tree_file = new TFile (Form("entree/root/snemo_run-%d_LI.root", run_number), "READ");
  double charge_tree;
  double amplitude_tree;
  TTree* tree = (TTree*)tree_file->Get("Result_tree");
  gROOT->cd();
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge_tree);
  tree->SetBranchStatus("amplitude_tree",1);
  tree->SetBranchAddress("amplitude_tree", &amplitude_tree);


  TParameter<double> *param = new TParameter<double>("start_time", time);
  param = (TParameter<double>*)(tree_file->Get("start_time"));
  start_time = param->GetVal();
  std::cout << "time = " << start_time << '\n';
  gROOT->cd();

  std::array<std::array<TH1D*,8>, 5> histograms;

  for(int i=712 ; i<717; i++){
    for (size_t j = 1; j <interval.size(); j++){
      wall = "IT";
      if (j > interval.size()/2) {
	wall = "FR";
      }
      pic = j;
      histograms[i-712][j-1] = new TH1D(Form("om_%d_pic_%lu_%s",i,j,wall.c_str()),Form("om_%d_pic_%lu_%s",i,j,wall.c_str()),1000, 0, 180000);
    }
  }

  for(int entry = 0; entry<tree->GetEntries(); entry++){
    tree->GetEntry(entry);
    if(om_number>711){            
      for (size_t j = 1; j < interval.size() ; j++){
	if(time>interval.at(j-1) && time<interval.at(j) && amplitude_tree > 10){
	    histograms[om_number-712][j-1]->Fill(charge_tree);
	}
      }
    }
  }  
  
  int number = 712;
  int indice = 0;
  for(int om = number; om < number + 5; om+=1)
    {
      for (size_t j = 1; j <  interval.size(); j++) 
	{//boucle pour selection side
	  //std::cout<<"SIDE = "<<j<<" size "<< interval.size()<<std::endl;
	  wall = "IT";
	  if (j > interval.size()/2) {
	    wall = "FR";
	  }	  
	  pic =j;
	  //serie de condition pour verifier que courbe LI est fittable
	  if (histograms[om-number][j-1]->GetEntries() < 300) {

	    std::cout << "" << '\n';
	    std::cout << "trop peu d'entries pour l'OM " << om << '\n';
	    std::cout << "" << '\n';

	    Khi2 = 0;
	    om_number = om;
	    mean = 0;
	    mean_error = 0;
	    bundle_ref = Ref_bundle_number(om, wall);
	    Result_tree.Fill();
	    delete histograms[om-number][j-1];
	  }
	  else if (histograms[om-number][j-1]->GetMean() > 180000) {//sature
	    std::cout << "" << '\n';
	    std::cout << "the amplitude sature" << '\n';
	    std::cout << "" << '\n';

	    Khi2 = 0;
	    om_number = om;
	    mean = 0;
	    mean_error = 0;
	    bundle_ref = Ref_bundle_number(om, wall);	    
	    Result_tree.Fill();
	    delete histograms[om-number][j-1];
	  }
	  else if (histograms[om-number][j-1]->GetMean() < 20){ //si pic trop bas -> charge<0ee

	    std::cout<<"mean ="<<histograms[om-number][j-1]->GetMean();
	    std::cout << "" << '\n';
	    std::cout << "too few charge" << '\n';
	    std::cout << "" << '\n';


	    Khi2 = 0;
	    om_number = om;
	    mean = 0;
	    mean_error = 0;
	    std::cout<<" pic "<<j<<"om_number "<<om_number<< " run num "<<run_number <<std::endl;      
	    bundle_ref = Ref_bundle_number(om, wall);
	    Result_tree.Fill();
	    delete histograms[om-number][j-1];
	  }
	  else{ // si trois conditions sont bonnes -> fait fit
	    nevent = histograms[om-number][j-1]->Integral();
	    TF1 *f_Gaus = new TF1("f_Gaus", "gaus(0)", 0, 1200);
	    f_Gaus->SetParNames("N_evt","mean_charge","Sigma");
	    // f_Gaus->SetParameters(25, spectre->GetMean(), 100);
	    f_Gaus->SetParameters(25, histograms[om-number][j-1]->GetMean(), histograms[om-number][j-1]->GetRMS());
	    f_Gaus->SetRange(histograms[om-number][j-1]->GetRMS()*5, histograms[om-number][j-1]->GetRMS()*5);
	    f_Gaus->Draw("same");
	    //fit puis range puis fit...
	    histograms[om-number][j-1]->Fit(f_Gaus, "RQ0");
	    f_Gaus->SetRange(f_Gaus->GetParameter(1)-5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
	    histograms[om-number][j-1]->Fit(f_Gaus, "RQ0");
	    f_Gaus->SetRange(f_Gaus->GetParameter(1)-5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
	    histograms[om-number][j-1]->Fit(f_Gaus, "RQ0");	    
	    Khi2 = f_Gaus->GetChisquare()/f_Gaus->GetNDF();
	    om_number = om;
	    constante = (f_Gaus->GetParameter(0));
	    mean = (f_Gaus->GetParameter(1));
	    mean_error = f_Gaus->GetParError(1);
	    mean_corr = (f_Gaus->GetParameter(1))/ref_gain_table[om-number];
	    //delta Q_ref/delta alpha_ref = delta_LED
	    mean_corr_error = sqrt(pow(f_Gaus->GetParError(1)*ref_gain_table[om-number],2) + pow(ref_gain_table_error[om-number]*f_Gaus->GetParameter(1),2) );
	    sigma = (f_Gaus->GetParameter(2));	   
	    intensity = intensity_chooser(pic);
	    canvas->cd();
	    histograms[om-number][j-1]->Draw();
	    f_Gaus->Draw("same");
	    //std::cout<<"om"<<om<<endl;
	    histograms[om-number][j-1]->GetXaxis()->SetRangeUser(mean-10*sigma,mean+10*sigma);
	    	    	    
	    if (om > 711 && run_number > 836) {	      
	      std::cout<<"run num = "<< run_number<<endl;
	      canvas->SaveAs(Form("sortie/ref_Li/fit_Li/om_%d/OM_%03d_pic_%d_run_%d%s%s.png", om, om, pic, run_number,wall.c_str(),corr.c_str()));
	      
	    }
	    bundle_ref = Ref_bundle_number(om, wall);
	    Result_tree.Fill();
	    
	    delete histograms[om-number][j-1];
	    delete f_Gaus;
	  }
	}
    }
  std::cout << "time = " << start_time << '\n';

  file.cd();
  Result_tree.Write();
  file.Close();
  return;
}





void file_merger_tot(std::vector<int> run){
 TFile *newfile = new TFile(Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/Merged_Fit_%d-%d.root", run[0],run[run.size()-1]), "RECREATE");
  
  double Amplitude, Amplitude_error, Amplitude_corr, Amplitude_corr_error, Khi2, time;
  int om_number, run_number, pic, bundle;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("Amplitude", &Amplitude);
  Result_tree.Branch("Amplitude_error", &Amplitude_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("time", &time);
  Result_tree.Branch("bundle", &bundle);
  Result_tree.Branch("Amplitude_corr", &Amplitude_corr);
  Result_tree.Branch("Amplitude_corr_error", &Amplitude_corr_error);

 
  for (size_t j = 0; j < run.size(); j++) {
    TFile *file1 = new TFile(Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/Amplitude_Li_run_%d.root",run[j]), "READ");    
    TTree* tree = (TTree*)file1->Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchStatus("pic",1);
    tree->SetBranchAddress("pic", &pic);
    tree->SetBranchStatus("Khi2",1);
    tree->SetBranchAddress("Khi2", &Khi2);
    tree->SetBranchStatus("Amplitude",1);
    tree->SetBranchAddress("Amplitude", &Amplitude);
    tree->SetBranchStatus("Amplitude_error",1);
    tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &run_number);
    tree->SetBranchStatus("time",1);
    tree->SetBranchAddress("time", &time);
    tree->SetBranchStatus("bundle",1);
    tree->SetBranchAddress("bundle", &bundle);
    
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      Result_tree.Fill();
    }

    file1->Close();
  }
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
}







void file_merger(std::vector<int> run, int ref = 0, string addfile = "",string corr="") {
  double Amplitude, Amplitude_error, Amplitude_corr, Amplitude_corr_error, Khi2, time;
  int om_number, run_number, pic, bundle_ref;
  string filename = "";
  string cut = "";
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("Amplitude", &Amplitude);
  Result_tree.Branch("Amplitude_error", &Amplitude_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("time", &time);
  Result_tree.Branch("Amplitude_corr", &Amplitude_corr);
  Result_tree.Branch("Amplitude_corr_error", &Amplitude_corr_error);
  Result_tree.Branch("bundle_ref", &bundle_ref);
    
  //  std::cout<< "size ="<< run.size();
  for (size_t j = 0; j < run.size(); j++) {
    TFile *file1 = new TFile(Form("sortie/ref_Li/Fit_Ampl_Ref/Amplitude_Li_%srun_%d%s.root", cut.c_str(), run[j],corr.c_str()), "READ");
    TTree* tree = (TTree*)file1->Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchStatus("pic",1);
    tree->SetBranchAddress("pic", &pic);
    tree->SetBranchStatus("Khi2",1);
    tree->SetBranchAddress("Khi2", &Khi2);
    tree->SetBranchStatus("Amplitude",1);
    tree->SetBranchAddress("Amplitude", &Amplitude);
    tree->SetBranchStatus("Amplitude_error",1);
    tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &run_number);
    tree->SetBranchStatus("time",1);
    tree->SetBranchAddress("time", &time);
    tree->SetBranchStatus("Amplitude_corr",1);
    tree->SetBranchAddress("Amplitude_corr", &Amplitude_corr);
    tree->SetBranchStatus("Amplitude_corr_error",1);
    tree->SetBranchAddress("Amplitude_corr_error", &Amplitude_corr_error);
    tree->SetBranchStatus("bundle_ref",1);
    tree->SetBranchAddress("bundle_ref", &bundle_ref);
    
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      Result_tree.Fill();
    }

    file1->Close();
  }
  pic = 0;
  Khi2 = 0;
  Amplitude = 0;
  Amplitude_error = 0;

  TFile *newfile = new TFile(Form("sortie/ref_Li/Fit_Ampl_Ref/Merged_Fit_%d-%d%s.root",run[0],run[run.size()-1],corr.c_str()), "RECREATE");
  std::cout << "file saved : " << Form("sortie/ref_Li/Fit_Ampl_Ref/Merged_Fit_%d-%d.root",run[0],run[run.size()-1])<< '\n';
  newfile->cd();
  Result_tree.Write();
  newfile->Close();

}









void Khi2selector(int file_run_number) { //selection du pic en fonction du chi2 -> virer les pics qui ont un mauvais chi2
  TFile *file = new TFile(Form("/home/granjon/Documents/these/Analyse/corr_Li_ref/sortie/ref_Li/Amplitude_Li/Amplitude_Li_run_%d.root", file_run_number), "READ");

  int om_number, pic, run_number, intensity;
  double Khi2, constante, Amplitude, Amplitude_error, sigma, time, nevent;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("Amplitude", &Amplitude);
  Result_tree.Branch("Amplitude_error", &Amplitude_error);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("intensity", &intensity);
  Result_tree.Branch("time", &time);
  Result_tree.Branch("nevent", &nevent);

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("pic",1);
  tree->SetBranchAddress("pic", &pic);
  tree->SetBranchStatus("Khi2",1);
  tree->SetBranchAddress("Khi2", &Khi2);
  tree->SetBranchStatus("Amplitude",1);
  tree->SetBranchAddress("Amplitude", &Amplitude);
  tree->SetBranchStatus("Amplitude_error",1);
  tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree->SetBranchStatus("sigma",1);
  tree->SetBranchAddress("sigma", &sigma);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);
  tree->SetBranchStatus("intensity",1);
  tree->SetBranchAddress("intensity", &intensity);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchStatus("nevent",1);
  tree->SetBranchAddress("nevent", &nevent);

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (Khi2 > 6) {
      Amplitude = 0;
      Amplitude_error = 0;
    }
    Result_tree.Fill();
  }
  TFile *newfile = new TFile(Form("/home/granjon/Documents/these/Analyse/corr_Li_ref/sortie/ref_Li/Amplitude_Li/Amplitude_Li_khicut_run_%d.root", run_number), "RECREATE");
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
}









void Evolution_Li_ref_graph(int n_run, int start, int stop,string wall,string corr=""){//comparaison graph corrige et pas corrige
  gSystem->mkdir(Form("sortie/ref_Li/variation_ref/run_%d_%d",start,stop));
  double amp_scor[5][8][n_run];
  double amp_scor_error[5][8][n_run];
  double amp_cor[5][8][n_run];
  double amp_cor_error[5][8][n_run];
  double Amplitude_error, Amplitude, time, Amplitude_corr, Amplitude_corr_error;
  int run_number, om_number, pic, Ref_error;
  double time_vec[n_run];

  int srun = start;//mettre le premier run ??
  int scompteur = 0;



  TFile file_cor(Form("sortie/ref_Li/Fit_Ampl_Ref/Merged_Fit_%d-%d.root",start,stop), "READ");
  TTree* tree_cor = (TTree*)file_cor.Get("Result_tree");
  tree_cor->SetBranchStatus("*",0);
  tree_cor->SetBranchStatus("Amplitude",1);
  tree_cor->SetBranchAddress("Amplitude", &Amplitude);
  tree_cor->SetBranchStatus("Amplitude_error",1);
  tree_cor->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree_cor->SetBranchStatus("Amplitude_corr",1);
  tree_cor->SetBranchAddress("Amplitude_corr", &Amplitude_corr);
  tree_cor->SetBranchStatus("Amplitude_corr_error",1);
  tree_cor->SetBranchAddress("Amplitude_corr_error", &Amplitude_corr_error);
  tree_cor->SetBranchStatus("run_number",1);
  tree_cor->SetBranchAddress("run_number", &run_number);
  tree_cor->SetBranchStatus("om_number",1);
  tree_cor->SetBranchAddress("om_number", &om_number);
  tree_cor->SetBranchStatus("pic",1);
  tree_cor->SetBranchAddress("pic", &pic);
  tree_cor->SetBranchStatus("time",1);
  tree_cor->SetBranchAddress("time", &time);

  srun = start;
  scompteur = 0;
    
  for (int i = 0; i < tree_cor->GetEntries(); i++) {
    tree_cor->GetEntry(i);
    if(i==0){
      time_vec[0]=time;	   
    }
    //while (run_number == 1060) { //to avoid a run number
    //i++;
    //tree_cor->GetEntry(i);
    //}
    if (run_number != srun) {
      scompteur++;
      srun = run_number;
      time_vec[scompteur]=time;
    }
    if(wall=="IT" && pic<5){
      amp_cor[om_number -712][pic-1][scompteur] = Amplitude_corr;
      amp_cor_error[om_number -712][pic-1][scompteur] = Amplitude_corr_error;
      amp_scor[om_number -712][pic-1][scompteur] = Amplitude;
      amp_scor_error[om_number -712][pic-1][scompteur] = Amplitude_error;
    }
    else if(wall=="FR" && pic>=5){
      //cout << "om " << om_number << " -> pic " << pic-5 << " run " << std::endl;
      amp_cor[om_number -712][pic-5][scompteur] = Amplitude_corr;
      amp_cor_error[om_number -712][pic-5][scompteur] = Amplitude_corr_error;      
      amp_scor[om_number -712][pic-5][scompteur] = Amplitude;
      amp_scor_error[om_number -712][pic-5][scompteur] = Amplitude_error;
    }
  }
  
   
  double norm_amp_cor[5][5][n_run];
  double norm_amp_cor_error[5][5][n_run];
  double norm_amp_scor[5][5][n_run];
  double norm_amp_scor_error[5][5][n_run];
  double color[5] = {kBlack,kBlue,kGreen+1,kOrange-3,kRed};
  for (int i = 0; i < 5; i++) {//om num
    for (int j = 0; j < 4; j++) {//pic num
      for (int l = 0; l < n_run; l++) {//run num
        if (amp_cor[i][j][0] > 0.1 && amp_scor[i][j][0] > 0.1) {
          norm_amp_cor[i][j][l] = amp_cor[i][j][l]/amp_cor[i][j][0];
          norm_amp_cor_error[i][j][l] = sqrt(pow(amp_cor_error[i][j][l]/amp_cor[i][j][0],2) + pow(amp_cor[i][j][l]*amp_cor_error[i][j][0]/pow(amp_cor[i][j][0],2),2));
	  norm_amp_scor[i][j][l] = amp_scor[i][j][l]/amp_scor[i][j][0];
          norm_amp_scor_error[i][j][l] = sqrt(pow(amp_scor_error[i][j][l]/amp_scor[i][j][0],2) + pow(amp_scor[i][j][l]*amp_scor_error[i][j][0]/pow(amp_scor[i][j][0],2),2));
        }
      }
    }
  }

  double xaxis[n_run];
  double xaxis_error[n_run];
  for(int indice = 0; indice<n_run; indice++){
    xaxis[indice]=time_vec[indice];
    xaxis_error[indice]=0;
  }  
  for (int i = 0; i < 5; i++) { //om num
  TMultiGraph *multiGraph = new TMultiGraph();
  TMultiGraph *multiGraph_amp = new TMultiGraph();
    string name = namer(712+i);   

  multiGraph->SetNameTitle(Form("fit_OM_ref_%d", 712+i), Form("LED variation for %s", name.c_str()));
  multiGraph->GetYaxis()->SetTitle("#Delta LED");
  multiGraph->GetXaxis()->SetTitle("date");
  multiGraph_amp->SetNameTitle(Form("fit_OM_ref_%d", 712+i), Form("LED variation for %s", name.c_str()));
  multiGraph_amp->GetYaxis()->SetTitle("#Delta LED");
  multiGraph_amp->GetXaxis()->SetTitle("date");
  double min_x = *std::min_element(xaxis,xaxis+n_run);
  double max_x = *std::max_element(xaxis,xaxis+n_run);
  multiGraph->GetYaxis()->SetRangeUser(0.9,1.1);
  multiGraph->GetXaxis()->SetTimeDisplay(1);    
  //multiGraph_amp->GetYaxis()->SetRangeUser(0.9,1.1);
  multiGraph_amp->GetXaxis()->SetTimeDisplay(1);    

  
  TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
  TLegend *legend2 = new TLegend(0.7,0.7,0.9,0.9);
    for (int j = 0; j < 4; j++) { // pic num
      TCanvas* c = new TCanvas;
      double yaxis[n_run];
      double yaxis_error[n_run];
      double yaxis_amp[n_run];
      double yaxis_error_amp[n_run];
      double yaxis_samp[n_run];
      double yaxis_error_samp[n_run];
      double syaxis[n_run];
      double syaxis_error[n_run];
      for (int l = 0; l < n_run; l++) {
        if (norm_amp_cor[i][j][l] >0.1 && norm_amp_scor[i][j][l] > 0.1) {
	  yaxis_amp[l] = amp_cor[i][j][l] - amp_cor[i][j][0];	  
	  yaxis_error_amp[l] = sqrt(amp_cor_error[i][j][l]*amp_cor_error[i][j][l] + amp_cor_error[i][j][0] * amp_cor_error[i][j][0]);
	  yaxis_samp[l] = amp_scor[i][j][l] - amp_scor[i][j][0];
	  yaxis_error_samp[l] = sqrt(amp_scor_error[i][j][l]*amp_scor_error[i][j][l] + amp_scor_error[i][j][0]*amp_scor_error[i][j][0]);
          yaxis[l] = norm_amp_cor[i][j][l];
          yaxis_error[l] = norm_amp_cor_error[i][j][l];	
          syaxis[l] = norm_amp_scor[i][j][l];
          syaxis_error[l] = norm_amp_scor_error[i][j][l];

        }
      }
      TGraphErrors *variation_scor = new TGraphErrors(n_run, xaxis, syaxis, xaxis_error, syaxis_error);
      TGraphErrors *variation_cor = new TGraphErrors(n_run, xaxis, yaxis, xaxis_error, yaxis_error);
      TGraphErrors *variation_amp = new TGraphErrors(n_run, xaxis, yaxis_amp, xaxis_error, yaxis_error_amp);
      TGraphErrors *variation_samp = new TGraphErrors(n_run, xaxis, yaxis_samp, xaxis_error, yaxis_error_samp);
    
      variation_amp->SetLineColor(color[j+1]);
      variation_amp->SetLineWidth(2);
      variation_samp->SetLineColor(color[j+1]);
      variation_samp->SetLineStyle(2);
      variation_samp->SetLineWidth(2);
      
      legend2->AddEntry(variation_amp, Form("Intensity %d", j+1), "l");
      variation_amp->Draw("APL"/*"same"*/);
      variation_samp->Draw("same");

      
      variation_scor->SetLineColor(color[j+1]);
      variation_scor->SetLineStyle(2);
      variation_scor->SetLineWidth(2);
      variation_cor->SetLineColor(color[j+1]);
      variation_cor->SetLineWidth(2);
      legend->AddEntry(variation_cor, Form("Intensity %d", j+1), "l");
      variation_cor->Draw("APL"/*"same"*/);
      variation_scor->Draw("same");
      //variation_cor->GetYaxis()->SetRangeUser(0.9,1.1);
      //variation_cor->GetXaxis()->SetTimeDisplay(1);
      //variation_cor->GetXaxis()->SetRangeUser(0,n_run);
  
  
      //c->SaveAs(Form("sortie/ref_Li/variation_ref/variation_om_%d_pic_%d.png",712+i,j+1));
      multiGraph->Add(variation_cor);
      multiGraph->Add(variation_scor);
      multiGraph_amp->Add(variation_amp);
      multiGraph_amp->Add(variation_samp);
      //variation_scor->SetTitle("test");
      //variation_scor->SetName("test");
      //variation_scor->SetNameTitle("test");
      //TFile newfile(Form("sortie/ref_Li/variation_ref/variation_om_%d_pic_%d.root", 712+i,j),"RECREATE");

      //newfile.cd();
      //variation_scor->Write();
      //variation_cor->Write();
      //newfile.Close();
      c->Close();
    }
    TCanvas *c = new TCanvas("c", "Graphique", 800, 600);
    c->cd();
    multiGraph_amp->GetYaxis()->SetRangeUser(-2500,2500);
    multiGraph_amp->Draw("APL");
    multiGraph_amp->GetXaxis()->SetLimits(min_x-10000,max_x+10000);
    //multiGraph_amp->GetYaxis()->UnZoom();
    c->Modified();

    
    TCanvas *canv = new TCanvas("canvas", "Graphique", 800, 600);
    canv->cd();
    multiGraph->GetXaxis()->UnZoom();
    multiGraph->Draw("APL");
    multiGraph->GetXaxis()->SetLimits(min_x-10000,max_x+10000);
    canv->Modified();
    legend->Draw("SAME");    
    if(wall == "IT"){
      canv->SaveAs(Form("sortie/ref_Li/variation_ref/run_%d_%d/variation_om_%d_wall_%s.png",start,stop,712+i,wall.c_str()));
      canv->SaveAs(Form("sortie/ref_Li/variation_ref/run_%d_%d/variation_om_%d_wall_%s.root",start,stop,712+i,wall.c_str()));
    }
    else if(wall=="FR"){
      canv->SaveAs(Form("sortie/ref_Li/variation_ref/run_%d_%d/variation_om_%d_wall_%s.png",start,stop,712+i,wall.c_str()));
      canv->SaveAs(Form("sortie/ref_Li/variation_ref/run_%d_%d/variation_om_%d_wall_%s.root",start,stop,712+i,wall.c_str()));
    }
    canv->Close();

    legend2->Draw("same");
    if(wall == "IT"){
      c->SaveAs(Form("sortie/ref_Li/variation_ref/run_%d_%d/variation_om_%d_wall_%s_amp.png",start,stop,712+i,wall.c_str()));
      c->SaveAs(Form("sortie/ref_Li/variation_ref/run_%d_%d/variation_om_%d_wall_%s_amp.root",start,stop,712+i,wall.c_str()));
    }
    else if(wall=="FR"){
      c->SaveAs(Form("sortie/ref_Li/variation_ref/run_%d_%d/variation_om_%d_wall_%s_amp.png",start,stop,712+i,wall.c_str()));
      c->SaveAs(Form("sortie/ref_Li/variation_ref/run_%d_%d/variation_om_%d_wall_%s_amp.root",start,stop,712+i,wall.c_str()));
    }
    c->Close();  
  }
}



void Evolution_Li_SN_graph(int n_run, int start, int stop){//comparaison graph corrige et pas corrige
  //plutot que multigraph tu peux dessiner les j pics de l'histo dans le meme canvas
  double amp_scor[712][4][n_run]; //om pic run
  double amp_scor_error[712][4][n_run];
  double amp_cor[712][4][n_run];
  double amp_cor_error[712][4][n_run];
  double Amplitude_error, Amplitude, time, Amplitude_corr, Amplitude_corr_error, Khi2;
  int run_number, om_number, pic, Ref_error;
  double time_vec[n_run];

  int srun = start;//mettre le premier run ??
  int scompteur = 0;
  TFile *file_sortie = new TFile("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/sortie.root","RECREATE");

  std::cout<<"file open "<<Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/Merged_Fit_%d-%d.root",start,stop)<<std::endl;
  TFile file_cor(Form("/home/granjon/Li/sortie/ref_Li/Fit_Ampl_Ref/final_final_gain_%d-%d.root",start,stop), "READ");
  TTree* tree_cor = (TTree*)file_cor.Get("Result_tree");
  tree_cor->SetBranchStatus("*",0);
  tree_cor->SetBranchStatus("Amplitude",1);
  tree_cor->SetBranchAddress("Amplitude", &Amplitude);
  tree_cor->SetBranchStatus("Amplitude_error",1);
  tree_cor->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree_cor->SetBranchStatus("Amplitude_corr",1);
  tree_cor->SetBranchAddress("Amplitude_corr", &Amplitude_corr);
  tree_cor->SetBranchStatus("Amplitude_corr_error",1);
  tree_cor->SetBranchAddress("Amplitude_corr_error", &Amplitude_corr_error);
  // tree_cor->SetBranchStatus("Khi2",1);
  // tree_cor->SetBranchAddress("Khi2", &Khi2);
  tree_cor->SetBranchStatus("run_number",1);
  tree_cor->SetBranchAddress("run_number", &run_number);
  tree_cor->SetBranchStatus("om_number",1);
  tree_cor->SetBranchAddress("om_number", &om_number);
  tree_cor->SetBranchStatus("pic",1);
  tree_cor->SetBranchAddress("pic", &pic);
  tree_cor->SetBranchStatus("time",1);
  tree_cor->SetBranchAddress("time", &time);
  
  srun = start;
  scompteur = 0;
    
  for (int i = 0; i < tree_cor->GetEntries(); i++) {
    tree_cor->GetEntry(i);
    if(i==0){
      time_vec[0]=time;	   
    }
    if (run_number != srun) {
      scompteur++;
      srun = run_number;
      time_vec[scompteur]=time;
    }   
      amp_cor[om_number][pic-1][scompteur] = Amplitude_corr;
      amp_cor_error[om_number][pic-1][scompteur] = Amplitude_corr_error;
      amp_scor[om_number][pic-1][scompteur] = Amplitude;
      amp_scor_error[om_number][pic-1][scompteur] = Amplitude_error;    
  }  
   
  double norm_amp_cor[712][4][n_run];
  double norm_amp_cor_error[712][4][n_run];
  double norm_amp_scor[712][4][n_run];
  double norm_amp_scor_error[712][4][n_run];
  double color[5] = {kBlack,kBlue,kGreen+1,kOrange-3,kRed};
  for (int i = 0; i < 712; i++) {//om num
    for (int j = 0; j < 4; j++) {//pic num
      for (int l = 0; l < n_run; l++) {//run num
        if (amp_scor[i][j][0] > 0.1 /*&& amp_cor[i][j][0] > 0.1*/) {
          norm_amp_cor[i][j][l] = amp_cor[i][j][l]/amp_cor[i][j][0];
          norm_amp_cor_error[i][j][l] = sqrt(pow(amp_cor_error[i][j][l]/amp_cor[i][j][0],2) + pow(amp_cor[i][j][l]*amp_cor_error[i][j][0]/pow(amp_cor[i][j][0],2),2));
	  norm_amp_scor[i][j][l] = amp_scor[i][j][l]/amp_scor[i][j][0];
          norm_amp_scor_error[i][j][l] = sqrt(pow(amp_scor_error[i][j][l]/amp_scor[i][j][0],2) + pow(amp_scor[i][j][l]*amp_scor_error[i][j][0]/pow(amp_scor[i][j][0],2),2));
        }
      }
    }
  }
  double xaxis[n_run];
  double xaxis_error[n_run];
  for(int indice = 0; indice<n_run; indice++){
    xaxis[indice]=time_vec[indice];
    xaxis_error[indice]=0;
    //std::cout<<" time_vec "<< time_vec[indice]<<std::endl;
  }
    double min_x = *std::min_element(xaxis,xaxis+n_run);
    double max_x = *std::max_element(xaxis,xaxis+n_run);

  for (int i = 0; i < 712; i++) { //om num
    TMultiGraph *multiGraph = new TMultiGraph(Form("om_%d",i),Form("om_%d",i));
    TMultiGraph *multiGraph_amp = new TMultiGraph(Form("om_amp_%d",i),Form("om_amp_%d",i));
    multiGraph->GetYaxis()->SetTitle("#Delta Q");
    multiGraph->GetXaxis()->SetTitle("date");
    multiGraph_amp->GetYaxis()->SetTitle("#Delta LED");
    multiGraph_amp->GetXaxis()->SetTitle("date");
    multiGraph->GetYaxis()->SetRangeUser(0.9,1.1);
    multiGraph->GetXaxis()->SetTimeDisplay(1);    
    //multiGraph_amp->GetYaxis()->SetRangeUser(0.9,1.1);
    multiGraph_amp->GetXaxis()->SetTimeDisplay(1);    
    //std::cout<<"min x = "<<min_x<<" max x "<<max_x<<std::endl;
    for (int j = 0; j < 4; j++) { // pic num
      TCanvas* c = new TCanvas;
      double yaxis[n_run];
      double yaxis_error[n_run];
      double yaxis_amp[n_run];
      double yaxis_error_amp[n_run];
      double yaxis_samp[n_run];
      double yaxis_error_samp[n_run];
      double syaxis[n_run];
      double syaxis_error[n_run];
      for (int l = 0; l < n_run; l++) {
        if (norm_amp_cor[i][j][l] >0.1 && norm_amp_scor[i][j][l] > 0.1) {
	  yaxis_amp[l] = amp_cor[i][j][l] - amp_cor[i][j][0];	  
	  yaxis_error_amp[l] = sqrt(amp_cor_error[i][j][l]*amp_cor_error[i][j][l] + amp_cor_error[i][j][0] * amp_cor_error[i][j][0]);
	  yaxis_samp[l] = amp_scor[i][j][l] - amp_scor[i][j][0];
	  yaxis_error_samp[l] = sqrt(amp_scor_error[i][j][l]*amp_scor_error[i][j][l] + amp_scor_error[i][j][0]*amp_scor_error[i][j][0]);
          yaxis[l] = norm_amp_cor[i][j][l];
          yaxis_error[l] = norm_amp_cor_error[i][j][l];	
          syaxis[l] = norm_amp_scor[i][j][l];
          syaxis_error[l] = norm_amp_scor_error[i][j][l];
	  //std::cout<<"var scor"<<norm_amp_scor[i][j][l]<<" var samp "<<amp_scor[i][j][l] - amp_scor[i][j][0]<<std::endl;
        }
      }

      TGraphErrors *variation_scor = new TGraphErrors(n_run, xaxis, syaxis, xaxis_error, syaxis_error);
      TGraphErrors *variation_cor = new TGraphErrors(n_run, xaxis, yaxis, xaxis_error, yaxis_error);
      TGraphErrors *variation_amp = new TGraphErrors(n_run, xaxis, yaxis_amp, xaxis_error, yaxis_error_amp);
      TGraphErrors *variation_samp = new TGraphErrors(n_run, xaxis, yaxis_samp, xaxis_error, yaxis_error_samp);
    
      variation_amp->SetLineColor(color[j+1]);
      variation_amp->SetLineWidth(2);
      variation_samp->SetLineColor(color[j+1]);
      variation_samp->SetLineStyle(2);
      variation_samp->SetLineWidth(2);
      variation_samp->SetTitle(Form("#Delta Q Intensity %d", j+1));
      variation_amp->SetTitle(Form("#Delta #alpha Intensity %d", j+1));
      variation_amp->Draw("APL"/*"same"*/);
      variation_samp->Draw("APL"/*"same"*/);
      variation_scor->SetLineColor(color[j+1]);
      variation_scor->SetLineStyle(2);
      variation_scor->SetLineWidth(2);
      variation_scor->SetTitle(Form("#Delta Q Intensity %d", j+1));
      variation_cor->SetLineColor(color[j+1]);
      variation_cor->SetLineWidth(2);
      variation_cor->SetTitle(Form("#Delta #alpha Intensity %d", j+1));
      

      variation_cor->Draw("APL"/*"same"*/);
      variation_scor->Draw("APL"/*"same"*/);
      multiGraph->Add(variation_cor);
      multiGraph->Add(variation_scor);
      multiGraph_amp->Add(variation_amp);
      multiGraph_amp->Add(variation_samp);
    
    }//boucle j
    file_sortie->cd();
    //std::cout<<"om number " <<i<<std::endl;
    multiGraph_amp->GetXaxis()->SetLimits(min_x-10000,max_x+10000);
    multiGraph_amp->GetYaxis()->SetRangeUser(-2500,2500);
    multiGraph->GetXaxis()->SetLimits(min_x-10000,max_x+10000);
    multiGraph->Write();
    multiGraph_amp->Write();
  }//boucle i
  file_sortie->Close();
}






int main(int argc, char const *argv[]){
  int n_run, run;
  std::vector<int> run_number, ref_run_number, ref_time, energy_run_number;
  int compteur = 0;
  string file, ref_correction, calo_correction;
  bool add = false;
  bool energy = false;

  for(int i = 0; i<argc; i++){
    if (std::string(argv[i]) == "-add" ) {
      file = argv[i+1];
      std::cout << file << '\n';
      add = true;
    }
    if (std::string(argv[i]) == "-energy" ) {
      energy = true;
    }
  }

  std::cout << "How many run do you want ?" << '\n';
  std::cin >> n_run;
  std::cout << "Write the runs you want" << '\n';
  while (compteur < n_run && cin >> run) {
    run_number.push_back(run);
    compteur++;
    if (compteur < n_run) {
      std::cout << "Write the runs you want" << '\n';
    }
  }
  compteur = 0;
  if (energy == true) {
    std::cout << "Write the energy runs you want" << '\n';
    while (compteur < n_run && cin >> run) {
      energy_run_number.push_back(run);
      compteur++;
      if (compteur < n_run) {
        std::cout << "Write the energy runs you want" << '\n';
      }
    }
  }

  std::cout << "Ref Correction file name ?" << '\n';
  std::cin >> ref_correction;
  // std::cout << "Calo Correction file name ?" << '\n';
  // std::cin >> calo_correction;

  
  double* ref_gain_tab_base = new double[5];
  double* ref_gain_tab = new double[5];
  double* ref_gain_tab_error = new double[5];
  double* ref_gain_tab_base_error = new double[5];


  
  //rempli ref_gain_tab_base
  Ref_corrector(run_number[0], ref_correction, ref_gain_tab_base, ref_gain_tab_base_error);
  //tu corriges par rapport a un run qui n'est pas celui de refence 
  //ref_gain_tab_base = tableau ref
  //Ref_correctrur = remplir tableau gain
  //le tableau est remplir pour le run 0 et pour les suivants
  //run 0 = ref
  for (size_t i = 0; i < 5; i++) {
   std::cout << "ref = " << ref_gain_tab_base[i] << '\n';
  }
  compteur = 0;
  n_run = 0;

  std::cout << "Code start running" << '\n';

  for (size_t i = 0; i < run_number.size(); i++) {
    //std::cout<<" ref corr" << ref_correction << " ref gain tab" << *ref_gain_tab << " ref gain tab error " << *ref_gain_tab_error<<endl;
    //rempli le ref_gain_tab
    Ref_corrector(run_number[i], ref_correction, ref_gain_tab, ref_gain_tab_error);    
    std::cout << "Ref_Corrector "<< run_number[i] << " is ok" << '\n';
    // if (i > 0) {
    
    // for (int j = 0; j < 5; j++) {
    //   std::cout<<"ancien gain = "<<ref_gain_tab[j]<<std::endl;
    //     ref_gain_tab[j] = ref_gain_tab[j]/ref_gain_tab_base[j];
    // 	//ref gain tab te donne la variation a corriger pour le run etudie (i)
    //     //cout << ref_gain_tab[j] << " - " << ref_gain_tab_base[j] << endl;
    //     ref_gain_tab_error[j] = sqrt(pow(ref_gain_tab_error[j]/ref_gain_tab_base[j],2) + pow(ref_gain_tab[j]*ref_gain_tab_base_error[j]/(pow(ref_gain_tab_base[j],2)),2));
    // }
    // }       
    fit_LI_amplitude_Ref(run_number[i], ref_gain_tab, ref_gain_tab_error);
    
    //fit_LI_amplitude_Ref(run_number[i], ref_gain_tab, ref_gain_tab_error,"_sans_corr");
	
    std::cout << "fit LI amplitud Ref "<< run_number[i] << " is ok" << '\n';
  }
  
  if (add == false) {
    file_merger(run_number, 1);
    //file_merger(run_number, 1, file,"_sans_corr");//add peut etre pas bon
  }
  else{
    file_merger(run_number, 1, file);//add peut etre pas bon
  }
  Evolution_Li_ref_graph(run_number.size(),run_number[0],run_number[run_number.size()-1],"IT");
  Evolution_Li_ref_graph(run_number.size(),run_number[0],run_number[run_number.size()-1],"FR");
  

  
  
  
  //Passage au calo 
   std::cout << "" << '\n';
   std::cout << "START OF THE CALORIMETER FIT" << '\n';
   std::cout << "" << '\n';
  
   double gain_tab[40];  //correction lumiere par OM ref
   for (size_t i = 0; i < run_number.size(); i++) {
     fit_LI_amplitude(run_number[i], gain_tab);
     std::cout<<"fit_LI_amplitude ok"<<std::endl;
   }
   //std::cout<<"ICI"<<run_number[0]<<" "<<run_number[run_number.size()-1]<<std::endl;
   file_merger_tot(run_number);
   std::cout<<"File merger tot ok"<<std::endl;
   Li_corrector(run_number,run_number.size(),run_number[0],run_number[run_number.size()-1]);
   std::cout<<"Li corrector ok"<<std::endl;
   applied_Li_correction(run_number,run_number.size(),run_number[0],run_number[run_number.size()-1]);
   std::cout<<"applied_Li_correction ok"<<std::endl;
   Evolution_Li_SN_graph(run_number.size(),run_number[0],run_number[run_number.size()-1]);
   std::cout<<"end"<<std::endl;
  return 0;
}








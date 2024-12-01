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
#include <TSpectrum.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TSystem.h>
#include "/home/granjon/Documents/stage_radon/Stage/Mathis/sndisplay/sndisplay.cc"


using namespace std;
//summer 2023 cara
// constexpr int n_run_bi = 17;                                           
// int vecteur_Bi[n_run_bi] = {1090,1097,1105,1112,1119,1126,1133,1149,1156,1163,1170,1175,1180,1190,1195,1200,1205};

//February 2022 cara
// constexpr int n_run_bi = 3;                                           
// int vecteur_Bi[n_run_bi] = {1048, 1055, 1059};

//March 2023 cara
// constexpr int n_run_bi = 22;
// int vecteur_Bi[n_run_bi] ={1274,1275,1276,1277,1281,1286,/*1287, pas chassis tracker centre*/1292,/*1297 only 1 min,*/1303,1313,1320,/*1326too short*/ /*1327 root is truncated*/1337,1351,1359,1367,/*1374 root has problem*/1381,1389,1396,1403,1412,1413,1414,1415/*1416 too short + source move*/};
// //int vecteur_Bi[3] = {1274,1275,1276};

//total 2023-2024 cara
constexpr int n_run_bi=42;
int vecteur_Bi[n_run_bi] = {1090,1097,1105,1112,1119,1126,1133,1149,1156,1163,1170,1175,1180,1190,1195,1200,1205,1048, 1055, 1059,1274,1275,1276,1277,1281,1286,/*1287, pas chassis tracker centre*/1292,/*1297 only 1 min,*/1303,1313,1320,/*1326too short*/ /*1327 root is truncated*/1337,1351,1359,1367,/*1374 root has problem*/1381,1389,1396,1403,1412,1413,1414,1415/*1416 too short + source move*/};
  


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
    intensity = 80;
  }
  if (pic == 2) {
    intensity = 90;
  }
  if (pic == 3) {
    intensity = 100;
  }
  if (pic == 4) {
    intensity = 110;
  }
  return intensity;
}

std::vector<double> time_measurer(int run_number){ // this function find the number of intensities -> when the time spectra goes to 0 -> we change our intensities
  TFile tree_file(Form("entree/root/snemo_run-%d_LI.root", run_number), "READ");
  int om_number;
  double time = 0.0;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  gROOT->cd();
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);

  TH1D *time_spectre = new TH1D ("time_spectre", "", 2000, 0, 1100);
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    time_spectre->Fill(time);
  }

  int plus =0;
  int moins = 0;
  std::vector<double> time_measurer;
  int i;
  for (i = 0; i < 2000; i++) {
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

  // std::cout << Form("calcul_gain_ref/root/Merged_Fit/Fit_Ref_%s.root", correction.c_str()) << '\n';

  int run_number;
  double gain= 0.0, gain_error_plus= 0.0, gain_error_moins= 0.0 ;
  TTree* tree = (TTree*)file1.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  //if pic only
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);  
  tree->SetBranchStatus("gain_error_plus",1);
  tree->SetBranchAddress("gain_error_plus", &gain_error_plus);
  tree->SetBranchStatus("gain_error_moins",1);
  tree->SetBranchAddress("gain_error_moins", &gain_error_moins);  
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);

  int compteur = 0;
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);    
    //mettre abs ici ???
    if (run_number == run+1 || run_number == run -1 ) {//pour associer run correction avec LI
      //std::cout<<" RUN FIT REF  "<<run_number<<" RUN LI  "<< run<<"\n";      
      gain_tab[compteur] = gain;
      gain_tab_error[compteur] = gain_error_plus+gain_error_moins;
      compteur++;
    }
  }

  file1.Close();
  //March things
  //string correction2 = "1231-1418";
  
  //Summer 2023 things
  //string correction2 = "1085-1209";

  //February things
  //  string correction2 = "999-1065";

  //Total 2023-2024 things
  string correction2 = "999-1418";
  
  TFile file2(Form("/home/granjon/Documents/these/Analyse/corr_OM_ref/sortie/root_file_final/Fit_Ref_%s.root", correction2.c_str()), "READ");
  int run_number2, om_number2;
  double gain2= 0.0, gain_error_plus2= 0.0, gain_error_moins2= 0.0;
  TTree* tree2 = (TTree*)file2.Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("gain",1);
  tree2->SetBranchAddress("gain", &gain2);  
  tree2->SetBranchStatus("gain_error_plus",1);
  tree2->SetBranchAddress("gain_error_plus", &gain_error_plus2);
  tree2->SetBranchStatus("gain_error_moins",1);
  tree2->SetBranchAddress("gain_error_moins", &gain_error_moins2);  
  tree2->SetBranchStatus("run_number",1);
  tree2->SetBranchAddress("run_number", &run_number2);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number2);
  int compteur2 = 0;
  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    if (run_number2 == run+1 || run_number2 == run -1 ) {
      if(om_number2==714){
      gain_tab[compteur2] = gain2;
      gain_tab_error[compteur2] = gain_error_plus2+gain_error_moins2;
      }
      compteur2++;
    }
  }
  file2.Close();
  return gain_tab;
}



int Ref_bundle_number(int om_number, string wall){ // associer un num ref a un bundle = partie du detecteur
  int bundle_number = 0;
  if (om_number == 712 && wall.compare("IT") == 0) {
    bundle_number = 1;
  }
  if (om_number == 713 && wall.compare("IT") == 0) {
    bundle_number = 2;
  }
  if (om_number == 714 && wall.compare("IT") == 0) {
    bundle_number = 4;
  }
  if (om_number == 715 && wall.compare("IT") == 0) {
    bundle_number = 3;
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





void suivi_gamma(int run_number){
  double mean_hist= 0.0;
  int om_number,nb_entries;
  TFile *file_create = new TFile(Form("~/Bi/sortie/gamma/result_run_gamma_%d.root",run_number), "RECREATE");
  TTree final_tree("final_tree","");
  final_tree.Branch("om_number",&om_number);
  final_tree.Branch("run_number",&run_number);
  final_tree.Branch("mean_hist",&mean_hist);
  final_tree.Branch("nb_entries",&nb_entries);

  TFile *file = new TFile(Form("~/Bi/entree/fichier_gamma_%d.root",run_number), "READ");
  TTree* Result_tree = (TTree*)file->Get("arbre_alpha");
  std::vector<int> *om_num_elec = new std::vector<int>;
  std::vector<double> *energy_elec = new std::vector<double>;
  
  gROOT->cd();
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("om_num_elec",1);
  Result_tree->SetBranchAddress("om_num_elec", &om_num_elec);
  Result_tree->SetBranchStatus("charsge_elec",1);
  Result_tree->SetBranchAddress("charge_elec", &energy_elec);

  std::array<TH1D*, 712> histograms;
  int nb_bins = 200;
  int max_charge = 80000;
  int min_charge = 0;
  for(int i=0; i<712; i++){
      histograms[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),nb_bins,min_charge,max_charge);
  }
  
  for(int entry=0; entry<Result_tree->GetEntries();entry++){
    Result_tree->GetEntry(entry);
    for(size_t i = 0; i < om_num_elec->size(); i++){
    histograms[om_num_elec->at(i)]->Fill(energy_elec->at(i)); 
    }
  }
    
  for(int i=0; i<712; i++){
    file_create->cd();
    histograms[i]->Write();  
    nb_entries = histograms[i]->GetEntries();
    mean_hist = histograms[i]->GetMean();
    om_number=i;
    run_number = run_number;
    final_tree.Fill();
  }
  final_tree.Write();
  file_create->Close();
}




void fit_Bi_energy(int run_number){
  double chi2= 0.0, mean_pic= 0.0, mean_pic_error= 0.0,time_bi= 0.0,ndf= 0.0;
  int om_number,nb_entries;
  TFile *file_create = new TFile(Form("~/Bi/sortie/result_run_Bi_%d.root",run_number), "RECREATE");
  TTree final_tree("final_tree","");
  final_tree.Branch("Chi2",&chi2);
  final_tree.Branch("om_number",&om_number);
  final_tree.Branch("run_number",&run_number);
  final_tree.Branch("mean_pic",&mean_pic);
  final_tree.Branch("mean_pic_error",&mean_pic_error);
  final_tree.Branch("time_bi",&time_bi);
  final_tree.Branch("nb_entries",&nb_entries);
  final_tree.Branch("ndf",&ndf);

  TFile *file = new TFile(Form("~/Bi/entree/fichier_Bi_%d.root",run_number), "READ");
  TTree* Result_tree = (TTree*)file->Get("arbre_alpha");
  int nb_calo_touch;
  std::vector<int> *om_num_elec = new std::vector<int>;
  std::vector<double> *energy_elec = new std::vector<double>;
  
  gROOT->cd();
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("nb_calo_touch",1);
  Result_tree->SetBranchAddress("nb_calo_touch", &nb_calo_touch);
  Result_tree->SetBranchStatus("om_num_elec",1);
  Result_tree->SetBranchAddress("om_num_elec", &om_num_elec);
  Result_tree->SetBranchStatus("charge_elec",1);
  Result_tree->SetBranchAddress("charge_elec", &energy_elec);
  Result_tree->SetBranchStatus("time_bi",1);
  Result_tree->SetBranchAddress("time_bi", &time_bi);


  std::array<TH1D*, 712> histograms;
  std::array<TH1D*, 712> histograms_peaks;
  int nb_bins = 200;
  int max_charge = 60000;
  int min_charge = 0;
  for(int i=0; i<712; i++){
    if(i<520){
    histograms[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),nb_bins,min_charge,max_charge);
    
    }
    else{
      histograms[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),nb_bins/2,min_charge,max_charge);
    }
  }
  
  for(int entry=0; entry<Result_tree->GetEntries();entry++){
    Result_tree->GetEntry(entry);
    if(nb_calo_touch==1){
      histograms[om_num_elec->front()]->Fill(energy_elec->front());
    }
  }
  
  int nb_down = 0, nb_pb_fit = 0, nb_good = 0;
  TF1* f_MultipleGaus = new TF1 ("f_MultipleGaus","[0]*(7.11*TMath::Gaus(x[0],[1],[2]) + 1.84*TMath::Gaus(x[0],[1]*(1047.8/975.7),[2]*sqrt((1047.8/975.7))) + 0.44*TMath::Gaus(x[0],([1]*(1059.8/975.7)),[2]*sqrt((1059.8/975.7))))", 0, 60000);
  file_create->cd();
  for(int i=0; i<712; i++){
    //fit                                                                                               
    //cout<<"hist"<<i<<endl;                                                                            
    //cout<<histograms[i][0]->GetEntries()<<endl;
    nb_entries = histograms[i]->GetEntries();
    if(histograms[i]->GetEntries()>500 && histograms[i]->GetMean()>0){
      TFitResultPtr fitResult;      
      for (int j = 0; j < 10; j++) {
        if(j!=0){
	  f_MultipleGaus->SetParLimits(2, 800, 8000);
	  //f_MultipleGaus->SetParLimits(1,histograms[i]->GetMaximumBin()*(max_charge/nb_bins)-histograms[i]->GetMaximumBin()*(max_charge/nb_bins)/5,histograms[i]->GetMaximumBin()*(max_charge/nb_bins)+histograms[i]->GetMaximumBin()*(max_charge/nb_bins)/5);
	  TSpectrum spectrum;
          if (spectrum.Search(histograms[i], 7, "", 0.11) == 2){
            Double_t *xPeaks = spectrum.GetPositionX();
            f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);
	  f_MultipleGaus->SetParLimits(1,xPeaks[0]-xPeaks[0]/5,xPeaks[0]+xPeaks[0]/5);
	  f_MultipleGaus->SetRange(f_MultipleGaus->GetParameter(1)-1.3*f_MultipleGaus->GetParameter(2), f_MultipleGaus->GetParameter(1)+3*f_MultipleGaus->GetParameter(2));
	  }
	  else{
	    Double_t *xPeaks = spectrum.GetPositionX();
            f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);
	    f_MultipleGaus->SetParLimits(1,xPeaks[0]-xPeaks[0]/5,xPeaks[0]+xPeaks[0]/5);
	    f_MultipleGaus->SetRange(f_MultipleGaus->GetParameter(1)-1.3*f_MultipleGaus->GetParameter(2), f_MultipleGaus->GetParameter(1)+3*f_MultipleGaus->GetParameter(2));
	  }
	}
	else{//if first time	  
	  TSpectrum spectrum;                                                                
	  if (spectrum.Search(histograms[i], 7, "", 0.11) == 2){
	    Double_t *xPeaks = spectrum.GetPositionX();
	    f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);		    
	  }
	  else{
	    Double_t *xPeaks = spectrum.GetPositionX();
            f_MultipleGaus->SetParameters(25,xPeaks[0],xPeaks[0]/10);
	    cout<<"histo "<<i<<"has " << spectrum.Search(histograms[i], 7, "", 0.11)<<" pics"<<endl;
	  }
	}
        fitResult = histograms[i]->Fit(f_MultipleGaus, "RQ0");
      }
      ndf = nb_bins*(f_MultipleGaus->GetParameter(1)+3*f_MultipleGaus->GetParameter(2)-(f_MultipleGaus->GetParameter(1)-1*f_MultipleGaus->GetParameter(2)))/max_charge; //500 bins de 0 a 4 MeV
     
      TCanvas* canvas = new TCanvas(Form("om_%d",i), Form("om_%d",i), 800, 600);
      if(fitResult == 1){
	//cout<<"om"<<i<<endl;	
        nb_pb_fit++;	
        canvas->cd();
        gStyle->SetOptFit(1111);
        histograms[i]->Draw();
	//canvas->SaveAs(Form("/home/granjon/Bi/histo/histo_%d.png",i));                     
	canvas->Write();
        canvas->Close();
	om_number=i;
        chi2 = 10000;
	mean_pic = 0;
	mean_pic_error = 0;
	run_number = run_number;
        time_bi = time_bi;
	ndf = 1;
        final_tree.Fill();
        //cout<<"hist "<<i<<"pb fit"<<endl;                                                             
      }
      else{
	//cout<<" om bon"<<i<<endl;
        //cout<<"histo"<<i<<"chi2 "<<f_MultipleGaus->GetChisquare()<<endl;                              
        nb_good++;       
	canvas->cd();
	gStyle->SetOptFit(1111);
	histograms[i]->Draw();
	f_MultipleGaus->Draw("same");
	//histograms[i]->ShowPeaks(5, "", 0.2);
	//gSystem->mkdir(Form("/home/granjon/Bi/histo/run_%d",run_number));
	//canvas->SaveAs(Form("/home/granjon/Bi/histo/run_%d/histo_%d.png",run_number,i));             
	canvas->Write();
	canvas->Close();
	om_number=i;
	chi2 = f_MultipleGaus->GetChisquare();
	mean_pic = f_MultipleGaus->GetParameter(1);
	mean_pic_error = f_MultipleGaus->GetParError(1);
	run_number = run_number;
	time_bi = time_bi;
	final_tree.Fill();
	//histograms[i][0]->Write();                                                           

        //cout<<"hist "<<i<<"good"<<endl;                                                               
      }
    }
    else{
      TCanvas* canvas = new TCanvas(Form("om_%d",i), Form("om_%d",i), 800, 600);
      canvas->cd();
      gStyle->SetOptFit(1111);
      histograms[i]->Draw();
      //canvas->SaveAs(Form("/home/granjon/Bi/histo/run_%d/histo_%d.png",run_number,i));
      canvas->Write();
      canvas->Close();
      f_MultipleGaus->SetParameters(0,0,0);
      om_number=i;
      chi2 = 2000;
      ndf = 1;
      run_number = run_number;
      time_bi = time_bi;
      mean_pic = 0;
      mean_pic_error = 0;
      final_tree.Fill();      
      //cout<<"hist without entry "<<i<<endl;                                                          
      nb_down++;
    }
    	histograms[i]->Delete();
  }
  final_tree.Write();
  file_create->Close();
  cout<<"nb good "<<nb_good<<" nb down "<<nb_down<<" nb pb fit "<<nb_pb_fit<<endl;
}






void fit_LI_amplitude(int run_number){ //pour tout le calo
  std::array<std::array<TH1D*,4>, 712> histograms; //712 OM et 4 pics
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  string wall=""; 
  TFile file(Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/Amplitude_Li_run_%d.root", run_number),"RECREATE");
  //std::cout<<"file open "<<Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/Amplitude_Li_run_%d.root", run_number)<<std::endl;
  
  
  std::vector<double> interval;
  interval = time_measurer(run_number);
  // for (size_t i = 0; i < interval.size(); i++) {
  //   std::cout << "interval[" << i << "] = " << interval.at(i) << '\n';
  // }

  
  
  int om_number=-1,bundle;
  double constante= 0.0;
  double mean= 0.0, time= 0.0, start_time= 0.0;
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
  fail_tree.Branch("om_number", &om);
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
  start_time = param->GetVal() /*- (25 * 365.25 * 86400)*/; //time shift because SetTimeDisplay from root has a different reference than the 1st february 1970
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
  //std::cout<<"1 ici "<<std::endl;
  for(int entry = 0; entry<tree->GetEntries(); entry++){
    int debut = 0;
    wall="IT";
    tree->GetEntry(entry);
    if(om_number<712 && om_number>=0){
      if ((om_number > 259 && om_number < 520) || (om_number > 583 && om_number < 647) || (om_number > 679 && om_number < 712) ) {
	debut = debut+interval.size()/2;
	wall="FR";
      }
      //std::cout<<"interval "<<interval.size()<<std::endl;
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
	/*	
	std::cout << "" << '\n';
	std::cout << "trop peu d'entries pour l'OM " << om_number << '\n';
	std::cout << "" << '\n';
	*/
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
	//histograms[om_number][j]->Draw();
	file.cd();
	TCanvas* canvas = new TCanvas(Form("om_%d_pic_%lu",om_number,j+1), Form("om_%d_pic_%lu",om_number,j+1), 800, 600);
	canvas->cd();
	histograms[om_number][j]->Draw();
	canvas->Write();
	canvas->Close();
	
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
	file.cd();
	TCanvas* canvas = new TCanvas(Form("om_%d_pic_%lu",om_number,j+1), Form("om_%d_pic_%lu",om_number,j+1), 800, 600);
	canvas->cd();
	histograms[om_number][j]->Draw();
	canvas->Write();
	canvas->Close();
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
	file.cd();
	//histograms[om_number][j]->Write();
	TCanvas* canvas = new TCanvas(Form("om_%d_pic_%lu",om_number,j+1), Form("om_%d_pic_%lu",om_number,j+1), 800, 600);
	canvas->cd();
	gStyle->SetOptFit(1111);
	histograms[om_number][j]->Draw();
	f_Gaus->Draw("same");
	canvas->Write();
	canvas->Close();
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
  start_time = param->GetVal()- (25 * 365.25 * 86400);
  gROOT->cd();

  std::array<std::array<TH1D*,8>, 5> histograms;

  for(int i=712 ; i<717; i++){
    for (size_t j = 1; j <interval.size(); j++){
      wall = "IT";
      if (j > interval.size()/2) {
	wall = "FR";
      }
      pic = j;
      histograms[i-712][j-1] = new TH1D(Form("om_%d_intensity_%lu_%s",i,j,wall.c_str()),Form("om_%d_intensity_%lu_%s",i,j,wall.c_str()),1000, 0, 180000);
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
  double color[5] = {kBlack,kBlue,kGreen+1,kOrange-3,kRed}; //change if more intensities
  for(int om = number; om < number + 5; om+=1)
    {
      //si plot des 4 en meme temps
      // TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
      // TH1D *his = new TH1D(Form("om_%d",om),Form("om_%d",om),1000, 0, 120000);
      // his->GetYaxis()->SetRangeUser(0,270);
      // c2->cd();
      // his->Draw();
      for (size_t j = 1; j <  interval.size(); j++){//boucle pour selection side
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
	    mean_corr = 0;			
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
	    mean_corr =	0;			
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
	    mean_corr = 0;
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
	    intensity = intensity_chooser(pic); //change if more intensities
	    canvas->cd();
	    histograms[om-number][j-1]->Draw();
	    f_Gaus->Draw("same");
	    //std::cout<<"om"<<om<<endl;
	    histograms[om-number][j-1]->GetXaxis()->SetRangeUser(mean-10*sigma,mean+10*sigma);
	    if (om > 711 && run_number > 836) {	      
	      std::cout<<"run num = "<< run_number<<endl;
	      canvas->SaveAs(Form("sortie/ref_Li/fit_Li/om_%d/OM_%03d_intensity_%d_run_%d%s%s.png", om, om, pic, run_number,wall.c_str(),corr.c_str()));	      
	    }
	    //si plot des 4 pics en meme temps
	    // if(j<5){
	    // c2->cd();
	    // gStyle->SetOptFit(0);
	    // gStyle->SetOptStat(1);
	    // histograms[om-number][j-1]->Draw("sames");
	    // TPaveStats *st = (TPaveStats*)histograms[om-number][j-1]->FindObject("stats");
	    // st->SetY1NDC(1-(j-1)/5.-0.1);
	    // st->SetY2NDC(1-(j-1)/5.+0.1);
	    // st->SetLineColor(color[j]);
	    // f_Gaus->SetLineColor(color[j]);
	    // f_Gaus->Draw("same");
	    // c2->Modified();
	    // c2->Update();
	    // c2->SaveAs(Form("test_dans_%d.root",j));
	    // }
	    bundle_ref = Ref_bundle_number(om, wall);
	    Result_tree.Fill();	    
	    delete histograms[om-number][j-1];
	    delete f_Gaus;
	  }
      }
    }

  file.cd();
  Result_tree.Write();
  file.Close();
  return;
}



void Li_corrector(std::vector<int> run, int n_run, int start, int stop){
  double Amplitude, Amplitude_error, Amplitude_corr, Amplitude_corr_error, Khi2, time, gain, gain_error;
  int om_number, run_number, pic,bundle_ref;
  string wall;
  TFile *file_sortie = new TFile(Form("/home/granjon/Li/sortie/ref_Li/Fit_Ampl_Ref/Li_corr_%d-%d.root",start,stop),"RECREATE"); 
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
	gain =  amp_cor[i][j][l]/amp_cor[i][j][0];//delta LED	  
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
  double Amplitude= 0.0, Amplitude_error= 0.0, Amplitude_corr= 0.0, Amplitude_corr_error= 0.0, Khi2= 0.0, time= 0.0, gain= 0.0, gain_error= 0.0;
  int om_number, run_number, pic, bundle;

  TFile *file_sortie = new TFile(Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/SN_gain_%d-%d.root",start,stop),"RECREATE");
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
 TFile *file = new TFile(Form("/home/granjon/Li/sortie/ref_Li/Fit_Ampl_Ref/Li_corr_%d-%d.root",start,stop), "READ");
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


void file_merger_Bi(std::vector<int> run){
 TFile *newfile = new TFile(Form("~/Bi/sortie/Merged_Bi_Fit_%d-%d.root", run[0],run[run.size()-1]), "RECREATE");
 double mean_pic= 0.0, mean_pic_error= 0.0, Chi2= 0.0, time_bi= 0.0, ndf= 0.0;
  int om_number, run_number, pic, bundle,nb_entries;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("mean_pic", &mean_pic);
  Result_tree.Branch("mean_pic_error", &mean_pic_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("time_bi", &time_bi);
  Result_tree.Branch("nb_entries",&nb_entries);
  Result_tree.Branch("ndf",&ndf);

  for (size_t j = 0; j < run.size(); j++) {
    TFile *file1 = new TFile(Form("~/Bi/sortie/result_run_Bi_%d.root",run[j]), "READ");    
    TTree* tree = (TTree*)file1->Get("final_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchStatus("Chi2",1);
    tree->SetBranchAddress("Chi2", &Chi2);
    tree->SetBranchStatus("mean_pic",1);
    tree->SetBranchAddress("mean_pic", &mean_pic);
    tree->SetBranchStatus("mean_pic_error",1);
    tree->SetBranchAddress("mean_pic_error", &mean_pic_error);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &run_number);
    tree->SetBranchStatus("time_bi",1);
    tree->SetBranchAddress("time_bi", &time_bi);
    tree->SetBranchStatus("nb_entries",1);
    tree->SetBranchAddress("nb_entries", &nb_entries);
    tree->SetBranchStatus("ndf",1);
    tree->SetBranchAddress("ndf", &ndf);
    
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




void file_merger_tot(std::vector<int> run){
 TFile *newfile = new TFile(Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/Merged_Fit_%d-%d.root", run[0],run[run.size()-1]), "RECREATE");
  
  double Amplitude= 0.0, Amplitude_error= 0.0, Amplitude_corr= 0.0, Amplitude_corr_error= 0.0, Khi2= 0.0, time= 0.0;
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
  string* wall;
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
  Result_tree.Branch("wall", &wall);
    
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
    tree->SetBranchStatus("wall",1);
    tree->SetBranchAddress("wall", &wall);
    
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
  newfile->cd();
  Result_tree.Write();
  newfile->Close();

}




void Evolution_Li_ref_graph(int n_run, int start, int stop,string wall,string corr=""){//comparaison graph corrige et pas corrige
  gSystem->mkdir(Form("sortie/ref_Li/variation_ref/run_%d_%d",start,stop));
  double amp_scor[5][8][n_run]; //always five because there are five ref OMs
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
    TCanvas *c = new TCanvas("c", "Graphique", 1600, 600);
    c->cd();
    multiGraph_amp->GetYaxis()->SetRangeUser(-2500,2500);
    multiGraph_amp->GetXaxis()->SetNdivisions(520);
    multiGraph_amp->Draw("APL");
    multiGraph_amp->GetXaxis()->SetLimits(min_x-10000,max_x+10000);
    //multiGraph_amp->GetYaxis()->UnZoom();
    c->Modified();

    
    TCanvas *canv = new TCanvas("canvas", "Graphique", 1600, 600);
    canv->cd();
    multiGraph->GetXaxis()->UnZoom();
    multiGraph->GetXaxis()->SetNdivisions(520);
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



void Evolution_Li_SN_graph(std::vector<int> run, int n_run, int start, int stop,std::vector<int> run_bi){//comparaison graph corrige et pas corrige
  //plutot que multigraph tu peux dessiner les j pics de l'histo dans le meme canvas
  //Bi part
  int n_run_bi = run_bi.size();
  double mean_pic_value[712][n_run_bi];
  double mean_pic_value_error[712][n_run_bi];
  double mean_pic_vec[712][n_run_bi];
  double Chi2_vec[712][n_run_bi];
  double ndf_vec[712][n_run_bi];
  double time_bi_vec[712][n_run_bi];
  int nb_entries_vec[712][n_run_bi];
  double mean_pic_error_vec[712][n_run_bi];
  double mean_pic_error, mean_pic, time_bi, Chi2_bi, bi_relative, bi_relative_error,mean_pic_bi, mean_pic_bi_error,Chi2, ndf_bi, ndf, time_bi_f;
  int run_number_bi, om_number_bi, start_bi, srun_bi,nb_entries,nb_entries_bi;
  double time_vec_bi[n_run_bi];
  
  TFile *bi_file = new TFile(Form("~/Bi/sortie/Bi_final_%d-%d.root",start,stop),"RECREATE");
  TTree bi_tree("bi_tree","");
  bi_tree.Branch("om_number_bi", &om_number_bi);
  bi_tree.Branch("mean_pic_bi", &mean_pic_bi);
  bi_tree.Branch("mean_pic_bi_error", &mean_pic_bi_error);
  bi_tree.Branch("run_number_bi", &run_number_bi);
  bi_tree.Branch("bi_relative",&bi_relative);
  bi_tree.Branch("bi_relative_error",&bi_relative_error);
  bi_tree.Branch("Chi2_bi",&Chi2_bi);
  bi_tree.Branch("nb_entries_bi",&nb_entries_bi);
  bi_tree.Branch("ndf_bi",&ndf_bi);
  bi_tree.Branch("time_bi_f",&time_bi_f);

  TFile file_bi(Form("~/Bi/sortie/Merged_Bi_Fit_%d-%d.root",run_bi[0],run_bi[run_bi.size()-1]), "READ");
  TTree* tree_bi = (TTree*)file_bi.Get("Result_tree");
  tree_bi->SetBranchStatus("*",0);
  tree_bi->SetBranchStatus("mean_pic",1);
  tree_bi->SetBranchAddress("mean_pic", &mean_pic);
  tree_bi->SetBranchStatus("mean_pic_error",1);
  tree_bi->SetBranchAddress("mean_pic_error", &mean_pic_error);
  tree_bi->SetBranchStatus("run_number",1);
  tree_bi->SetBranchAddress("run_number", &run_number_bi);
  tree_bi->SetBranchStatus("om_number",1);
  tree_bi->SetBranchAddress("om_number", &om_number_bi);
  tree_bi->SetBranchStatus("time_bi",1);
  tree_bi->SetBranchAddress("time_bi", &time_bi);
  tree_bi->SetBranchStatus("Chi2",1);
  tree_bi->SetBranchAddress("Chi2", &Chi2);
  tree_bi->SetBranchStatus("nb_entries",1);
  tree_bi->SetBranchAddress("nb_entries", &nb_entries);
  tree_bi->SetBranchStatus("ndf",1);
  tree_bi->SetBranchAddress("ndf", &ndf);

  int scompteur_bi=0;
  srun_bi = run_bi[0];
  for (int i = 0; i < tree_bi->GetEntries(); i++) {
    tree_bi->GetEntry(i);
    if(i==0){
      time_vec_bi[0]=time_bi;
    }
    if (run_number_bi != srun_bi) {
      srun_bi = run_number_bi;
      scompteur_bi++;
      time_vec_bi[scompteur_bi]=time_bi;
    } 
      mean_pic_vec[om_number_bi][scompteur_bi] = mean_pic;
      mean_pic_error_vec[om_number_bi][scompteur_bi] = mean_pic_error;
      Chi2_vec[om_number_bi][scompteur_bi] = Chi2;      
      nb_entries_vec[om_number_bi][scompteur_bi] = nb_entries;
      ndf_vec[om_number_bi][scompteur_bi] = ndf;
      time_bi_vec[om_number_bi][scompteur_bi] = time_bi;
  }
  
  double xaxis_bi[n_run_bi];
  double xaxis_error_bi[n_run_bi];
  for(int indice = 0; indice<n_run_bi; indice++){
    xaxis_bi[indice]=time_vec_bi[indice];
    xaxis_error_bi[indice]=0;
  }
    double min_x_bi = *std::min_element(xaxis_bi,xaxis_bi+n_run_bi);
    double max_x_bi = *std::max_element(xaxis_bi,xaxis_bi+n_run_bi);

    std::array<TGraphErrors*, 712> variation_bi_vec;

    for (int i = 0; i < 712; i++) { //om num       
      double yaxis_bi[n_run_bi];
      double yaxis_error_bi[n_run_bi];
      for (int l = 0; l < n_run_bi; l++) {
	if (mean_pic_vec[i][0] > 0.1){
	  yaxis_bi[l] = mean_pic_vec[i][l]/mean_pic_vec[i][0];
	  yaxis_error_bi[l] = sqrt(pow(mean_pic_error_vec[i][l]/mean_pic_vec[i][0],2) + pow(mean_pic_vec[i][l]*mean_pic_error_vec[i][0]/pow(mean_pic_vec[i][0],2),2));
	  mean_pic_value[i][l] = mean_pic_vec[i][l]/mean_pic_vec[i][0];
	  mean_pic_value_error[i][l] = sqrt(pow(mean_pic_error_vec[i][l]/mean_pic_vec[i][0],2) + pow(mean_pic_vec[i][l]*mean_pic_error_vec[i][0]/pow(mean_pic_vec[i][0],2),2));
	  om_number_bi = i;
	  mean_pic_bi = mean_pic_vec[i][l];
	  mean_pic_bi_error = mean_pic_error_vec[i][l];
	  run_number_bi = vecteur_Bi[l];
	  bi_relative = mean_pic_value[i][l];
	  bi_relative_error = mean_pic_value_error[i][l];
	  Chi2_bi = Chi2_vec[i][l];
	  ndf_bi = ndf_vec[i][l];
	  nb_entries_bi=nb_entries_vec[i][l];
	  time_bi_f = time_bi_vec[i][l];
	  //cout<<nb_entries_bi<<endl;
	  bi_tree.Fill();
	}
	else{
	  yaxis_bi[l]=0.;
	  yaxis_error_bi[l]=0.;
	  mean_pic_value[i][l] = 0.;
	  mean_pic_value_error[i][l] = 0.;
	  om_number_bi = i;
	  mean_pic_bi = mean_pic_vec[i][l];
	  mean_pic_bi_error = 0;
	  run_number_bi = vecteur_Bi[l];
	  bi_relative = 0.;
	  bi_relative_error = 0.;
	  Chi2_bi = Chi2_vec[i][l];	  
	  nb_entries_bi=nb_entries_vec[i][l];
	  ndf_bi = ndf_vec[i][l];
	  time_bi_f = time_bi_vec[i][l];
	  bi_tree.Fill();	  
	}
      }
      variation_bi_vec[i] = new TGraphErrors(n_run_bi, xaxis_bi, yaxis_bi, xaxis_error_bi, yaxis_error_bi);
      variation_bi_vec[i]->SetLineColor(kBlack);
      variation_bi_vec[i]->SetLineWidth(3);
      variation_bi_vec[i]->SetTitle("Bi 1 MeV pic");
    }
    bi_file->cd();
    bi_tree.Write();
    bi_file->Close();

  //Li part
  double amp_scor[712][4][n_run]; //om pic run
  double amp_scor_error[712][4][n_run];
  double amp_cor[712][4][n_run];
  double amp_cor_error[712][4][n_run];
  double Amplitude_error, Amplitude, time, Amplitude_corr, Amplitude_uncorr=0.0, Amplitude_corr_error, Khi2, gain, gain_non_corr, time_li;
  int run_number, om_number, pic, Ref_error, bundle;
  double time_vec[n_run];
  double time_li_vec[712][4][n_run];
  int bundle_ref_vec[712][4][n_run];
  int srun = start;//mettre le premier run ??
  int scompteur = 0;
  TFile *file_sortie = new TFile(Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/histo_graph_%d-%d.root",start,stop),"RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Amplitude", &Amplitude);
  Result_tree.Branch("Amplitude_uncorr", &Amplitude_uncorr);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("gain",&gain);
  Result_tree.Branch("gain_non_corr",&gain_non_corr);
  Result_tree.Branch("bundle",&bundle);
  Result_tree.Branch("time_li",&time_li);

  TFile file_cor(Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/SN_gain_%d-%d.root",start,stop), "READ");
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
  tree_cor->SetBranchStatus("bundle",1);
  tree_cor->SetBranchAddress("bundle", &bundle);  
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
      bundle_ref_vec[om_number][pic-1][scompteur] = bundle;
      time_li_vec[om_number][pic-1][scompteur] = time;
  }  

  double norm_amp_cor[712][4][n_run];
  double norm_amp_cor_error[712][4][n_run];
  double norm_amp_scor[712][4][n_run];
  double norm_amp_scor_error[712][4][n_run];
  double color[5] = {kBlack,kBlue,kGreen+1,kOrange-3,kRed};
  for (int i = 0; i < 712; i++) {//om num
    for (int j = 0; j < 4; j++) {//pic num
      for (int l = 0; l < n_run; l++) {//run num
        if (amp_scor[i][j][0] > 0.1 && amp_cor[i][j][0] > 0.1 && amp_scor[i][j][l]>0.1 && amp_cor[i][j][0]>0.1) {
          norm_amp_cor[i][j][l] = amp_cor[i][j][l]/amp_cor[i][j][0]; //divide the run l by the first one to normalise
          norm_amp_cor_error[i][j][l] = sqrt(pow(amp_cor_error[i][j][l]/amp_cor[i][j][0],2) + pow(amp_cor[i][j][l]*amp_cor_error[i][j][0]/pow(amp_cor[i][j][0],2),2));
	  norm_amp_scor[i][j][l] = amp_scor[i][j][l]/amp_scor[i][j][0]; //do the same without correction
          norm_amp_scor_error[i][j][l] = sqrt(pow(amp_scor_error[i][j][l]/amp_scor[i][j][0],2) + pow(amp_scor[i][j][l]*amp_scor_error[i][j][0]/pow(amp_scor[i][j][0],2),2));
	  gain = norm_amp_cor[i][j][l];
	  gain_non_corr = norm_amp_scor[i][j][l];
	  run_number = run[l];
	  Amplitude = amp_cor[i][j][l];
	  Amplitude_uncorr = amp_scor[i][j][l];
	  bundle = bundle_ref_vec[i][j][l];
	  time_li = time_li_vec[i][j][l];
	  om_number = i;
	  pic =j+1;
	  Result_tree.Fill();
        }
	else{
          gain = 0.;
          gain_non_corr = 0.;
          run_number = run[l];
          Amplitude = amp_cor[i][j][l];
	  Amplitude_uncorr = amp_scor[i][j][l];
          bundle = bundle_ref_vec[i][j][l];
          om_number = i;
          pic =j+1;
	  time_li = time_li_vec[i][j][l];		    
          Result_tree.Fill();
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
    multiGraph->GetYaxis()->SetTitle("Relative variation");
    multiGraph->GetXaxis()->SetTitle("date");
    multiGraph->GetYaxis()->SetRangeUser(0.9,1.1);
    multiGraph->GetXaxis()->SetTimeDisplay(1);    
    for (int j = 0; j < 4; j++) { // pic num
      TCanvas* c = new TCanvas("","",1200, 600);                        
      double yaxis[n_run];
      double yaxis_error[n_run];
      double syaxis[n_run];
      double syaxis_error[n_run];
      for (int l = 0; l < n_run; l++) {
        if (norm_amp_cor[i][j][l] >0.1 && norm_amp_scor[i][j][l] > 0.1) {
          yaxis[l] = norm_amp_cor[i][j][l];
          yaxis_error[l] = norm_amp_cor_error[i][j][l];	
          syaxis[l] = norm_amp_scor[i][j][l];
          syaxis_error[l] = norm_amp_scor_error[i][j][l];
        }
	else{
	  yaxis[l]=0;
	  yaxis_error[l]=0;
	  syaxis[l]=0;
	  syaxis_error[l]=0;
	}	
      }
      TGraphErrors *variation_scor = new TGraphErrors(n_run, xaxis, syaxis, xaxis_error, syaxis_error);
      TGraphErrors *variation_cor = new TGraphErrors(n_run, xaxis, yaxis, xaxis_error, yaxis_error);
      variation_scor->SetLineColor(color[j+1]);
      variation_scor->SetLineStyle(2);
      variation_scor->SetLineWidth(2);
      //variation_scor->SetTitle(Form("Intensity %d", j+1));
      variation_cor->SetLineColor(color[j+1]);
      variation_cor->SetLineWidth(2);
      variation_cor->SetTitle(Form("gain variation intensity %d", j+1));
      variation_cor->Draw("APL"/*"same"*/);
      variation_scor->SetTitle(Form("charge variation intensity %d", j+1));
      variation_scor->Draw("APL"/*"same"*/);
      multiGraph->Add(variation_cor);
      multiGraph->Add(variation_scor);
    }//boucle j
    multiGraph->Add(variation_bi_vec[i]);    
    file_sortie->cd();
    multiGraph->GetXaxis()->SetLimits(min_x-10000,max_x+10000);
    //cout<<"multi"<<multiGraph<<endl;
    multiGraph->Write();
    variation_bi_vec[i]->Delete();
  }//boucle i
  Result_tree.Write();
  file_sortie->Close();
}


void comparaison_Bi_Li(int start, int stop, int stop_bi){
  double gain, gain_non_corr, bi_relative, bi_relative_error, li_relative, li_relative_error, diff_li_bi,Chi2,Chi2_bi, li_relative_non_corr, diff_li_bi_non_corr, shift_Li_with_1, shift_Bi_with_1, ndf_bi, time_li, time_bi,time_bi_f, mean_pic_bi, li_amplitude, li_amplitude_uncorr=0.0,Amplitude, diff_bi_start_stop, diff_li_start_stop, Amplitude_uncorr=0.0;
  int bundle, run_number, pic, om_number, om_number_bi, om_number_li, pic_li,bundle_li, run_number_bi, run_number_li, nb_entries_bi;
  bool is_before;

  TFile *file_final = new TFile(Form("/home/granjon/Bi/sortie/Comparison_Bi_Li/final_comparison_%d-%d.root",start,stop),"RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number_bi", &om_number_bi);
  Result_tree.Branch("om_number_li", &om_number_li);
  Result_tree.Branch("pic_li", &pic_li);
  Result_tree.Branch("bundle_li", &bundle_li);
  Result_tree.Branch("run_number_li", &run_number_li);
  Result_tree.Branch("run_number_bi", &run_number_bi);
  Result_tree.Branch("bi_relative",&bi_relative);
  Result_tree.Branch("bi_relative_error",&bi_relative_error);
  Result_tree.Branch("li_relative",&li_relative);
  Result_tree.Branch("li_relative_error",&li_relative_error);
  Result_tree.Branch("li_relative_non_corr",&li_relative_non_corr);
  Result_tree.Branch("diff_li_bi",&diff_li_bi);
  Result_tree.Branch("is_before",&is_before);
  Result_tree.Branch("Chi2_bi",&Chi2_bi);
  Result_tree.Branch("nb_entries_bi",&nb_entries_bi);
  Result_tree.Branch("diff_li_bi_non_corr",&diff_li_bi_non_corr);
  Result_tree.Branch("shift_with_Li_1",&shift_Li_with_1);
  Result_tree.Branch("shift_with_Bi_1",&shift_Bi_with_1);
  Result_tree.Branch("ndf_bi",&ndf_bi);
  Result_tree.Branch("time_li",&time_li);
  Result_tree.Branch("time_bi",&time_bi);
  Result_tree.Branch("li_amplitude",&li_amplitude);
  Result_tree.Branch("li_amplitude_uncorr",&li_amplitude_uncorr);
  Result_tree.Branch("mean_pic_bi",&mean_pic_bi);
  Result_tree.Branch("diff_bi_start_stop",&diff_bi_start_stop);
  Result_tree.Branch("diff_li_start_stop",&diff_li_start_stop);  
  
  
  TFile file_cor(Form("/home/granjon/Li/sortie/SN_Li/Amplitude_Li/histo_graph_%d-%d.root",start,stop), "READ");    
  TTree* tree_cor = (TTree*)file_cor.Get("Result_tree");
  tree_cor->SetBranchStatus("*",0);
  tree_cor->SetBranchStatus("gain",1);
  tree_cor->SetBranchAddress("gain", &gain);
  tree_cor->SetBranchStatus("gain_non_corr",1);
  tree_cor->SetBranchAddress("gain_non_corr", &gain_non_corr);
  tree_cor->SetBranchStatus("bundle",1);
  tree_cor->SetBranchAddress("bundle", &bundle);
  tree_cor->SetBranchStatus("run_number",1);
  tree_cor->SetBranchAddress("run_number", &run_number);
  tree_cor->SetBranchStatus("om_number",1);
  tree_cor->SetBranchAddress("om_number", &om_number);
  tree_cor->SetBranchStatus("pic",1);
  tree_cor->SetBranchAddress("pic", &pic);
  tree_cor->SetBranchStatus("time_li",1);
  tree_cor->SetBranchAddress("time_li", &time_li);
  tree_cor->SetBranchStatus("Amplitude",1);
  tree_cor->SetBranchAddress("Amplitude", &Amplitude);
  tree_cor->SetBranchStatus("Amplitude_uncorr",1);
  tree_cor->SetBranchAddress("Amplitude_uncorr", &Amplitude_uncorr);
 
  TFile file_bi(Form("~/Bi/sortie/Bi_final_%d-%d.root",start,stop), "READ");
  TTree* tree = (TTree*)file_bi.Get("bi_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("bi_relative",1);
  tree->SetBranchAddress("bi_relative", &bi_relative);
  tree->SetBranchStatus("bi_relative_error",1);
  tree->SetBranchAddress("bi_relative_error", &bi_relative_error);
  tree->SetBranchStatus("run_number_bi",1);
  tree->SetBranchAddress("run_number_bi", &run_number_bi);
  tree->SetBranchStatus("om_number_bi",1);
  tree->SetBranchAddress("om_number_bi", &om_number_bi);
  tree->SetBranchStatus("Chi2_bi",1);
  tree->SetBranchAddress("Chi2_bi", &Chi2);
  tree->SetBranchStatus("nb_entries_bi",1);
  tree->SetBranchAddress("nb_entries_bi", &nb_entries_bi);
  tree->SetBranchStatus("ndf_bi",1);
  tree->SetBranchAddress("ndf_bi", &ndf_bi);
  tree->SetBranchStatus("time_bi_f",1);
  tree->SetBranchAddress("time_bi_f", &time_bi_f);
  tree->SetBranchStatus("mean_pic_bi",1);
  tree->SetBranchAddress("mean_pic_bi", &mean_pic_bi);

 for (int j = 0; j < tree_cor->GetEntries(); j++) {
   tree_cor->GetEntry(j);
   diff_li_start_stop = 0;
   diff_bi_start_stop = 0;
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      om_number_li = om_number;
      om_number_bi = om_number_bi;
      if(om_number_li==om_number_bi){
	pic_li = pic;
	bundle_li=bundle;
	run_number_li = run_number;
	run_number_bi = run_number_bi;
	bi_relative = bi_relative;
	li_relative = gain;
	//cout<<"run num bi stop "<<run_number_bi<<" "<<stop_bi<<"run num li "<<run_number_li<<" "<<stop<<endl;
	if(run_number_bi==stop_bi && run_number_li==stop){
	  diff_li_start_stop = li_relative;
	}
        if(run_number_bi==stop_bi && run_number_li==stop){
	diff_bi_start_stop = bi_relative;
	}
	shift_Li_with_1 = fabs(1-gain);
	shift_Bi_with_1 = fabs(1-bi_relative);
	li_relative_non_corr = gain_non_corr;
	Chi2_bi = Chi2;
	nb_entries_bi = nb_entries_bi;
	ndf_bi = ndf_bi;
	time_li = time_li;
	time_bi = time_bi_f;
	mean_pic_bi = mean_pic_bi;
	li_amplitude = Amplitude;
	li_amplitude_uncorr = Amplitude_uncorr;
	if(bi_relative!=0){
	diff_li_bi = 100*(li_relative-bi_relative)/bi_relative;	
	diff_li_bi_non_corr = 100*(li_relative_non_corr-bi_relative)/bi_relative;
	}
	else{
	  diff_li_bi=0;
	  diff_li_bi_non_corr = 0;
	}
	//March 2024 cara
	// if(((run_number-run_number_bi)<=4 && (run_number-run_number_bi)>0) || (run_number_bi == 1412 && run_number_li==1417)){
	//   //cout<<"DIFF RUN  NUMBER "<<run_number-run_number_bi<<endl;
	//   Result_tree.Fill(); 
      
      //summer cara
	// if(run_number_bi==run_number+2 || (run_number_bi == 1105 && run_number_li == 1102)){
	//   is_before = true;
	//   Result_tree.Fill();
	// }
	// else if(run_number_bi==run_number-1){
	//   is_before = false;
	//   Result_tree.Fill();
	// }    
	//February cara
	if(run_number_bi == run_number-1){
	  Result_tree.Fill();	 
	}
      }
    }
 }
 file_final->cd();
 Result_tree.Write();
 file_final->Close();
}
 



int main(int argc, char const *argv[]){
  int n_run, run;
  std::vector<int> ref_run_number, ref_time, energy_run_number/*, run_number*/;
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


  //Summer 2023 cara
  // n_run = 29;
  // int run_number_before[n_run] = {1088,1091,1095,1098,1102,1106,1110,1113,1117,1120,1124,1127,1131,1134,1141,1144,1147,1150,1154,1157,1161,1164,1168,1173,1178,1188,1193,1198,1203};
  // std::vector<int> run_number(run_number_before, run_number_before + n_run);
  // ref_correction = "1085-1209";


  //February 2023 cara
  // n_run = 2;
  // int run_number_before[n_run] = {1049,/*1053 switch with the ref one*/1056};
  // std::vector<int> run_number(run_number_before, run_number_before + n_run);
  // ref_correction = "999-1065";
 
  
  //March 2024 cara
   // n_run = 19;
   // int run_number_before[n_run] = {1278,1282,1288,1293,1304,1314,1321,1328,1338,1352,1360,1362,1368,1375,1382,1390,1397,1406,1417};              
   // //n_run=1;
   // //int run_number_before[n_run] = {1091};
   //  std::vector<int> run_number(run_number_before, run_number_before + n_run);
   //  //March correction
  
    //ref_correction = "1231-1418_pic";
  //ref_correction = "1231-1418_alpha";

  //total 2023-2024 cara
  n_run = 50;
  int run_number_before[n_run] = {1049,/*1053 switch with the ref one*/1056,1088,1091,1095,1098,1102,1106,1110,1113,1117,1120,1124,1127,1131,1134,1141,1144,1147,1150,1154,1157,1161,1164,1168,1173,1178,1188,1193,1198,1203, 1278,1282,1288,1293,1304,1314,1321,1328,1338,1352,1360,1362,1368,1375,1382,1390,1397,1406,1417};
  std::vector<int> run_number(run_number_before, run_number_before + n_run);
  ref_correction = "999-1418";
   
    
  double* ref_gain_tab_base = new double[5];
  double* ref_gain_tab = new double[5];
  double* ref_gain_tab_error = new double[5];
  double* ref_gain_tab_base_error = new double[5];


  cout<<run_number[0]<<endl;
  //rempli ref_gain_tab_base
  Ref_corrector(run_number[0], ref_correction, ref_gain_tab_base, ref_gain_tab_base_error);
  //tu corriges par rapport a un run qui n'est pas celui de refence 
  //ref_gain_tab_base = tableau ref
  //Ref_correctrur = remplir tableau gain
  //le tableau est remplir pour le run 0 et pour les suivants
  //run 0 = ref

  compteur = 0;
  n_run = 0;

  std::cout << "Code start running" << '\n';

  for (size_t i = 0; i < run_number.size(); i++) {
    //rempli le ref_gain_tab
    Ref_corrector(run_number[i], ref_correction, ref_gain_tab, ref_gain_tab_error);    
    std::cout << "Ref_Corrector "<< run_number[i] << " is ok" << '\n';
    // if (i > 0) {
    
    for (int j = 0; j < 5; j++) {
        ref_gain_tab[j] = ref_gain_tab[j]/ref_gain_tab_base[j];
    	//On se refixe par rapport a une ref differente !!
    	//ref gain tab te donne la variation a corriger pour le run etudie (i)
        ref_gain_tab_error[j] = sqrt(pow(ref_gain_tab_error[j]/ref_gain_tab_base[j],2) + pow(ref_gain_tab[j]*ref_gain_tab_base_error[j]/(pow(ref_gain_tab_base[j],2)),2));
    }
    
    fit_LI_amplitude_Ref(run_number[i], ref_gain_tab, ref_gain_tab_error);        	
    std::cout << "fit LI amplitud Ref "<< run_number[i] << " is ok" << '\n';
  }
  
  if (add == false) {
    file_merger(run_number, 1);    
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

     for(int run_value : vecteur_Bi){
     fit_Bi_energy(run_value);
     std::cout<<"fit_BI ok "<<run_value<<std::endl;     
     }

     //changer vecteur statique en dynamique pour la fonction file_merger
     std::vector<int> vecteur_Bi_vector(std::begin(vecteur_Bi), std::end(vecteur_Bi));
     file_merger_Bi(vecteur_Bi_vector);
     cout<<"merger Bi ok"<<endl;
     cout<<"file "<<vecteur_Bi_vector[0]<<" - "<<vecteur_Bi_vector[vecteur_Bi_vector.size()-1] <<" created "<<endl;
     
   for (size_t i = 0; i < run_number.size(); i++) {
     fit_LI_amplitude(run_number[i]);
     std::cout<<"fit_LI_amplitude ok "<<run_number[i]<<std::endl;
   }
   file_merger_tot(run_number);
   std::cout<<"File merger tot ok"<<std::endl;
   Li_corrector(run_number,run_number.size(),run_number[0],run_number[run_number.size()-1]);
   std::cout<<"Li corrector ok"<<std::endl;
   applied_Li_correction(run_number,run_number.size(),run_number[0],run_number[run_number.size()-1]);
   std::cout<<"applied_Li_correction ok"<<std::endl;
   Evolution_Li_SN_graph(run_number,run_number.size(),run_number[0],run_number[run_number.size()-1],vecteur_Bi_vector);
   std::cout<<"evolution Li ok"<<std::endl;
   comparaison_Bi_Li(run_number[0],run_number[run_number.size()-1], vecteur_Bi[n_run_bi-1]);
   cout<<"comparison ok"<<endl;
   std::cout<<"end"<<std::endl;
  return 0;
}








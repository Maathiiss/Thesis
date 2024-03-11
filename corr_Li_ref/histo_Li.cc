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
#include <boost/filesystem.hpp>
#include <sys/stat.h>
#include <sys/types.h>
#include <TSystem.h>

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
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
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
    

    if (run_number == run+1) {//pour associer run correction avec LI
      std::cout<<" RUN FIT REF  "<<run_number<<" RUN LI  "<< run<<"\n";
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
  if ((om_number%13 > 7 && om_number/13 > 9 && om_number < 260) || (om_number > 671 && om_number < 680) || (om_number < 568 && om_number > 561) || (om_number < 584 && om_number > 577)) {
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



void Li_corrector(int run, double *gain_tab, int run_ref) { //recup run LI de ref et rempli tableau gain avec -> permet de corriger les OMs du calo de la variation de lumiere

  // TFile file(Form("root/Merged_Ref/Fit_Ref_%s.root", correction.c_str()), "READ");
  TFile file(Form("sortie/ref_Li/Fit_Ampl_Ref/Amplitude_Li_run_%d.root", run), "READ");
  if (run == run_ref) {
    for (int i = 0; i < 712; i++) {
      gain_tab[i] = 1;
    }
  }
  else{
    int run_number, pic, om_number;
    double Amplitude;
    double ref_gain_tab[10];
    string* wall;
    TTree* tree = (TTree*)file.Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("Amplitude",1);
    tree->SetBranchAddress("Amplitude", &Amplitude);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &run_number);
    tree->SetBranchStatus("pic",1);
    tree->SetBranchAddress("pic", &pic);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchStatus("wall",1);
    tree->SetBranchAddress("wall", &wall);

    // for (int i = 0; i < tree->GetEntries(); i++) {
    //   tree->GetEntry(i);
    //   if (run_number == run_ref && pic == 3){
    //     ref_gain_tab[Ref_bundle_number(om_number, wall)] = Amplitude;
    //   }
    // }
    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      if (run_number == run && pic == 3) {
        gain_tab[Ref_bundle_number(om_number, *wall)] = Amplitude/ref_gain_tab[Ref_bundle_number(om_number, *wall)];
      }
    }
  }
  file.Close();
  return;
}






void fit_LI_amplitude(int run_number, double *correction_gain_table){ //pour tout le calo
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  // TFile file(Form("test.root", run_number),"RECREATE");

  TFile file(Form("/sortie/ref_Li/Amplitude_Li/Amplitude_Li_run_%d.root", run_number),"RECREATE");


  std::vector<double> interval;
  interval = time_measurer(run_number);

  int om_number;
  double constante;
  double mean, time, start_time;
  double mean_error, nevent;
  double sigma;
  int pic =0;
  double Khi2 = 0;
  int intensity = 0;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
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

  double entries, saturation, under_threshold;

  TTree fail_tree("fail_tree","");
  fail_tree.Branch("entries", &entries);
  fail_tree.Branch("saturation", &saturation);
  fail_tree.Branch("under_threshold", &under_threshold);
  fail_tree.Branch("om_number", &om_number);
  fail_tree.Branch("pic", &pic);

  int debut = 0;
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

  for(int om = 0; om < 520; om++)
  {
    if ((om > 259 && om < 520) || (om > 583 && om < 647) || (om > 679 && om < 712) ) {debut = debut + (interval.size()/2);} //si IT -> change interval de temps
    for (double j = debut; j < debut + (interval.size()/2) ; j++)
    {
      TH1D *spectre = new TH1D ("spectre_amplitude", "", 700, 0, 2300 );
      tree->Project("spectre_amplitude", Form("amplitude_tree*%f", correction_gain_table[bundle_number(om)]), Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j), interval.at(j+1))); //peut etre inutile un des deux
      //avec correction
      if (j < debut + interval.size()-1) tree->Project("spectre_amplitude", Form("amplitude_tree*%f", correction_gain_table[bundle_number(om)]), Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j), interval.at(j+1)));  // si t'es dans 3 premiers pics
      else tree->Project("spectre_amplitude", Form("amplitude_tree*%f", correction_gain_table[bundle_number(om)]), Form("om_number == %d && time > %f && amplitude_tree > 10", om, interval.at(j)));

      //sans correction
      // if (j < debut + interval.size()-1) tree->Project("spectre_amplitude", Form("amplitude_tree", 1), Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j), interval.at(j+1)));
      // else tree->Project("spectre_amplitude", Form("amplitude_tree",1), Form("om_number == %d && time > %f && amplitude_tree > 10", om, interval.at(j)));

      if ((om > 259 && om < 520) || (om > 583 && om < 647) || (om > 679 && om < 712)) pic = j - (interval.size()/2); //condition sur les OMs  -> 4 pics it puis 4 pics fr 
      else pic = j;

      std::cout << spectre->GetEntries() << '\n';

      spectre->Draw();
      if (pic_selector(spectre) == -98540){//peut etre inutile
        cout << "pic selector" << endl;
        Khi2 = -1;
        constante = -1;
        mean = -1;
        mean_error = -1;
        sigma = -1;
        delete spectre;
      }

      else {
        if (spectre->GetEntries() < 300) {
          std::cout << "" << '\n';
          std::cout << "trop peu d'entries pour l'OM " << om << '\n';
          std::cout << "" << '\n';

          Khi2 = 0;
          om_number = om;
          mean = 0;
          mean_error = 0;
          if (om > 711 && run_number < 837) {
            om_number = om-88;
          }
          entries = spectre->GetEntries();
          saturation = 0;
          under_threshold = 0;
          pic = j;
          Result_tree.Fill();
          fail_tree.Fill();
          delete spectre;
        }
        else if (spectre->GetMean() > 1900) {
          std::cout << "" << '\n';
          std::cout << "the amplitude sature" << '\n';
          std::cout << "" << '\n';
          pic = j;
          Khi2 = 0;
          om_number = om;
          mean = 0;
          mean_error = 0;
          if (om > 711 && run_number < 837) {
            om_number = om-88;
          }
          entries = 0;
          saturation = spectre->GetMean();
          under_threshold = 0;
          pic = j;
          Result_tree.Fill();
          fail_tree.Fill();
          delete spectre;
        }
        else if (spectre->GetMean() < -10){
          std::cout << "" << '\n';
          std::cout << "too few charge" << '\n';
          std::cout << "" << '\n';
          pic = j;
          Khi2 = 0;
          om_number = om;
          mean = 0;
          mean_error = 0;
          if (om > 711 && run_number < 837) {
            om_number = om-88;
          }
          entries = 0;
          saturation = 0;
          under_threshold = spectre->GetMean();
          pic = j;
          Result_tree.Fill();
          fail_tree.Fill(); //quelles sont les pics et le LI qui deconnent
          delete spectre;
        }
        else{
          nevent = spectre->Integral();
          TF1 *f_Gaus = new TF1("f_Gaus", "gaus(0)", 0, 1200);
          f_Gaus->SetParNames("N_evt","mean_charge","Sigma");
          // f_Gaus->SetParameters(25, spectre->GetMean(), 100);
          f_Gaus->SetParameters(25, spectre->GetMean(), spectre->GetRMS());
          f_Gaus->SetRange(spectre->GetMean()-400, spectre->GetMean()+400);
          f_Gaus->Draw("same");
          spectre->Fit(f_Gaus, "RQ0");
          f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
          spectre->Fit(f_Gaus, "RQ0");
          f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
          spectre->Fit(f_Gaus, "RQ0");

          pic = j;
          if (run_number > 873) pic = j+1;
          Khi2 = f_Gaus->GetChisquare()/f_Gaus->GetNDF();
          om_number = om;
          constante = (f_Gaus->GetParameter(0));
          mean = (f_Gaus->GetParameter(1));
          sigma = (f_Gaus->GetParameter(2));
          mean_error = f_Gaus->GetParError(1) ;
          intensity = intensity_chooser(pic);
          spectre->Draw();
          f_Gaus->Draw("same");
          canvas->SaveAs(Form("sortie/ref_Li/fit_Li/amplitude_fit/om_%d/OM_%03d_pic_%d_run_%d.png", om, om_number, pic, run_number));

          Result_tree.Fill();

          delete spectre;
          delete f_Gaus;
        }
      }
    }
  debut = 0;
  }

  file.cd();
  Result_tree.Write();
  fail_tree.Write();
  file.Close();
  return;
}





void fit_LI_energy(int run_number, double *correction_gain_table, int energy_run){ //->pas utiliser
  //calibrer en energie avant de faire le fit
// void fit_LI_amplitude(int run_number){
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  // TFile file(Form("test.root", run_number),"RECREATE");

  int om_number, LED_number;
  double constante;
  double mean, calo_time, start_time;
  double mean_error, nevent, energy;
  double sigma;
  int pic =0;
  double Khi2 = 0;
  int intensity = 0;
  double entries, saturation, under_threshold;

  TFile *file = new TFile(Form("sortie/Fit_Energy/Energy_Li_run_%d.root", run_number),"RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("LED_number", &LED_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("constante", &constante);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("mean_energy", &energy);
  Result_tree.Branch("sigma", &sigma);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("intensity", &intensity);
  Result_tree.Branch("time", &start_time);
  Result_tree.Branch("nevent", &nevent);

  TTree fail_tree("fail_tree","");
  fail_tree.Branch("entries", &entries);
  fail_tree.Branch("saturation", &saturation);
  fail_tree.Branch("under_threshold", &under_threshold);
  fail_tree.Branch("om_number", &om_number);
  fail_tree.Branch("pic", &pic);
  gROOT->cd();

  double mean_energy;
  float energy_convertor[712];
  memset (energy_convertor, 0, 712*sizeof(float));
  TFile *energy_file = new TFile (Form("../Bi_selection/Bi_fit/fitted_bi_%d.root", energy_run), "READ");
  TTree* energy_tree = (TTree*)energy_file->Get("Result_tree");

  gROOT->cd();
  energy_tree->SetBranchStatus("*",0);
  energy_tree->SetBranchStatus("om_number",1);
  energy_tree->SetBranchAddress("om_number", &om_number);
  energy_tree->SetBranchStatus("mean",1);
  energy_tree->SetBranchAddress("mean", &mean_energy);

  for (int i = 0; i < energy_tree->GetEntries(); i++) {
    energy_tree->GetEntry(i);
    energy_convertor[om_number] = mean_energy;
  }
  energy_file->Close();

  std::vector<double> interval = time_measurer(run_number);
  std::cout << "time measurer ok" << '\n';

  TCanvas* canvas = new TCanvas;
  TFile *tree_file = new TFile (Form("entree/root/snemo_run-%d_LI.root", run_number), "READ");

  double charge_tree;
  double amplitude_tree;
  TTree* tree = (TTree*)tree_file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &calo_time);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge_tree);
  tree->SetBranchStatus("amplitude_tree",1);
  tree->SetBranchAddress("amplitude_tree", &amplitude_tree);

  TParameter<double> *param = (TParameter<double>*)(tree_file->Get("start_time"));
  start_time = param->GetVal();
  gROOT->cd();

  TH2D *it1 = new TH2D("it1","it1",712,0,712, 1000,0,250000);
  TH2D *it2 = new TH2D("it2","it2",712,0,712, 1000,0,250000);
  TH2D *it3 = new TH2D("it3","it3",712,0,712, 1000,0,250000);
  TH2D *it4 = new TH2D("it4","it4",712,0,712, 1000,0,250000);
  TH2D *fr1 = new TH2D("fr1","fr1",712,0,712, 1000,0,250000);
  TH2D *fr2 = new TH2D("fr2","fr2",712,0,712, 1000,0,250000);
  TH2D *fr3 = new TH2D("fr3","fr3",712,0,712, 1000,0,250000);
  TH2D *fr4 = new TH2D("fr4","fr4",712,0,712, 1000,0,250000);
  // TH2D h[2][4]
  cout << "size interval " << interval.size() << endl;
  for (int z = 0; z < tree->GetEntries(); z++) {
    tree->GetEntry(z);
    LED_number = bundle_number(om_number);
    if (calo_time > interval[0]  && calo_time < interval[1]) {
      it1->Fill(om_number, charge_tree*correction_gain_table[bundle_number(om_number)]);
    }
    else  if (calo_time > interval[1]  && calo_time < interval[2]) {
      it2->Fill(om_number, charge_tree*correction_gain_table[bundle_number(om_number)]);
    }
    else if (calo_time > interval[2]  && calo_time < interval[3]) {
      it3->Fill(om_number, charge_tree*correction_gain_table[bundle_number(om_number)]);
    }
    else if (calo_time > interval[3]  && calo_time < interval[4]) {
      it4->Fill(om_number, charge_tree*correction_gain_table[bundle_number(om_number)]);
    }
    else if (calo_time > interval[4]  && calo_time < interval[5]) {
      fr1->Fill(om_number, charge_tree*correction_gain_table[bundle_number(om_number)]);
    }
    else if (calo_time > interval[5]  && calo_time < interval[6]) {
      fr2->Fill(om_number, charge_tree*correction_gain_table[bundle_number(om_number)]);
    }
    else if (calo_time > interval[6]  && calo_time < interval[7]) {
      fr3->Fill(om_number, charge_tree*correction_gain_table[bundle_number(om_number)]);
    }
    else if (calo_time > interval[7]  && calo_time < interval[8]) {
      fr4->Fill(om_number, charge_tree*correction_gain_table[bundle_number(om_number)]);
    }
  }
  tree_file->Close();
  it1->Draw();
  // return;
  TH1D *spectre = nullptr;

  for(int om = 0; om < 520; om++)
  {
    int debut = 0;
    if ((om > 259 && om < 520) || (om > 583 && om < 647) || (om > 679 && om < 712) ) {debut = debut + (interval.size()/2);}
    for (size_t j = debut; j < debut + (interval.size()/2) ; j++){
      if (energy_convertor[om] != 0){
        if (j == 0) {
          spectre = it1->ProjectionY("pic1",om+1,om+1);
        }
        else if (j == 1) {
          spectre = it2->ProjectionY("pic1",om+1,om+1);
        }
        else if (j == 2) {
          spectre = it3->ProjectionY("pic1",om+1,om+1);
        }
        else if (j == 3) {
          spectre = it4->ProjectionY("pic1",om+1,om+1);
        }
        else if (j == 4) {
          spectre = fr1->ProjectionY("pic1",om+1,om+1);
        }
        else if (j == 5) {
          spectre = fr2->ProjectionY("pic1",om+1,om+1);
        }
        else if (j == 6) {
          spectre = fr3->ProjectionY("pic1",om+1,om+1);
        }
        else if (j == 7) {
          spectre = fr4->ProjectionY("pic1",om+1,om+1);
        }
        else {
          std::cout << "problème grave !!!!" << '\n';
          continue;
        }

        if ((om > 259 && om < 520) || (om > 583 && om < 647) || (om > 679 && om < 712)) pic = j - (interval.size()/2);
        else pic = j;

        canvas->cd();
        spectre->Draw();
        // if (pic_selector(spectre) == -98540){
        //   cout << "pic selector" << endl;
        //   Khi2 = -1;
        //   constante = -1;
        //   mean = -1;
        //   mean_error = -1;
        //   sigma = -1;
        //   delete spectre;
        // }
        // else {
          if (spectre->GetEntries() < 300) {
            std::cout << "" << '\n';
            std::cout << "trop peu d'entries pour l'OM " << om << " : " << spectre->GetEntries() << '\n';
            std::cout << "" << '\n';

            Khi2 = 0;
            om_number = om;
            mean = 0;
            mean_error = 0;
            if (om > 711 && run_number < 837) {
              om_number = om-88;
            }
            entries = spectre->GetEntries();
            saturation = 0;
            under_threshold = 0;
            pic = j;
            Result_tree.Fill();
            fail_tree.Fill();
            delete spectre;
          }
          else if (spectre->GetMean() > 200000) {
            std::cout << "" << '\n';
            std::cout << "the energy sature" << '\n';
            std::cout << "" << '\n';
            pic = j;
            Khi2 = 0;
            om_number = om;
            mean = 0;
            mean_error = 0;
            if (om > 711 && run_number < 837) {
              om_number = om-88;
            }
            entries = 0;
            saturation = spectre->GetMean();
            under_threshold = 0;
            pic = j;
            Result_tree.Fill();
            fail_tree.Fill();
            delete spectre;
          }
          else if (spectre->GetMean() < 10){
            std::cout << "" << '\n';
            std::cout << "too few charge" << '\n';
            std::cout << "" << '\n';
            pic = j;
            Khi2 = 0;
            om_number = om;
            mean = 0;
            mean_error = 0;
            if (om > 711 && run_number < 837) {
              om_number = om-88;
            }
            entries = 0;
            saturation = 0;
            under_threshold = spectre->GetMean();
            pic = j;
            Result_tree.Fill();
            fail_tree.Fill();
            delete spectre;
          }
          else{
            nevent = spectre->Integral();
            TF1 *f_Gaus = new TF1("f_Gaus", "gaus(0)", 0, 250000);
            f_Gaus->SetParNames("N_evt","mean_charge","Sigma");
            // f_Gaus->SetParameters(25, spectre->GetMean(), 100);
            f_Gaus->SetParameters(spectre->GetMaximum(), spectre->GetMean(), spectre->GetRMS());
            f_Gaus->SetRange(spectre->GetMean()-4000, spectre->GetMean()+4000);
            f_Gaus->Draw("same");
            spectre->Fit(f_Gaus, "RQ0");
            f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+2.5*f_Gaus->GetParameter(2));
            spectre->Fit(f_Gaus, "RQ0");
            f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+2.5*f_Gaus->GetParameter(2));
            spectre->Fit(f_Gaus, "RQ0");

            pic = j;
            if (run_number > 873) pic = j+1;
            Khi2 = f_Gaus->GetChisquare()/f_Gaus->GetNDF();
            om_number = om;
            constante = (f_Gaus->GetParameter(0));
            mean = (f_Gaus->GetParameter(1));
            sigma = (f_Gaus->GetParameter(2));
            mean_error = f_Gaus->GetParError(1);
            energy = f_Gaus->GetParameter(1)/energy_convertor[om];
            intensity = intensity_chooser(pic);
            spectre->Draw();
            f_Gaus->Draw("same");

            canvas->SaveAs(Form("fit/fit_Li/energy_fit/om_%d/OM_%03d_pic_%d_run_%d.png", om, om_number, pic, run_number));

            Result_tree.Fill();

            delete spectre;
            delete f_Gaus;
          }
        // }
      }
    }
  }

  file->cd();
  Result_tree.Write();
  fail_tree.Write();
  file->Close();
  return;
}





void fit_LI_amplitude_Ref(int run_number, double *ref_gain_table, double *ref_gain_table_error){
  //Fit amplitude LI des OM de ref corrige de leurs gains ->fit en charge
  
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  //TFile file(Form("/home/granjon/Documents/these/Analyse/corr_Li_ref/sortie/Amplitude_Li/Amplitude_Li_run_%d_noncor.root", run_number),"RECREATE");

  TFile file(Form("/home/granjon/Li/sortie/ref_Li/Fit_Ampl_Ref/Amplitude_Li_run_%d.root", run_number),"RECREATE");
  std::vector<double> interval;
  interval = time_measurer(run_number); // sert a conserver que les pics qui t'interesse -> plus de 4 pics car utilisation de led secondaires
  // plutot que de faire 4 fits -> on coupe en temps en fonction de la duree du fichier -> te donne les bornes a partir de quand on est dans premier pic, celui d'apres ... 
  
  for (size_t i = 0; i < interval.size(); i++) {
    std::cout << "interval[" << i << "] = " << interval.at(i) << '\n';
  }
  // return;
  
  int om_number;
  double constante;
  double mean, time, start_time;
  double mean_error, nevent, ref_error;
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

  int number = 712;
  /*
  if (run_number > 873) { //saute le premier pic car avant on avait 6 pics maintenant on en a 4
    pic = 1;
  }
  */

  for(int om = number; om < number + 5; om+=1)
    {
      for (size_t j = 1; j <  interval.size(); j++) 
	{//boucle pour selection side
	  std::cout<<"SIDE = "<<j<<" size "<< interval.size()<<std::endl;
	  wall = "IT";
	  if (j > interval.size()/2) {
	    wall = "FR";
	  }
	  TH1D *spectre = new TH1D ("spectre_amplitude", "", 300, 0, 120000/*70000*/);//la ou tu vas stocker donnes LI
	  
	  // avec corrections
	  //if (j < interval.size()/2){//si mur FR
	    //attention ici a la division 
	    tree->Project("spectre_amplitude", Form("charge_tree/%f", ref_gain_table[om-number]), Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j-1), interval.at(j))); // divise par valeur de ref
	    std::cout<<Form("charge_tree/%f", ref_gain_table[om-number]) << " + "<< Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j-1), interval.at(j));
	    //} 

	  /*
	  else {//si mur IT
	    tree->Project("spectre_amplitude", Form("charge_tree/%f", ref_gain_table[om-number]), Form("om_number == %d && time > %f && time < 2000 && amplitude_tree > 10", om, interval.at(j)));
	  }
	  */
	  //sans corrections
	  //if (j < debut + interval.size()-1) tree->Project("spectre_amplitude", "charge_tree", Form("om_number == %d && time > %f && time < %f && amplitude_tree > 10", om, interval.at(j), interval.at(j+1)));
	  //project LI
	  //else tree->Project("spectre_amplitude", "charge_tree", Form("om_number == %d && time > %f && time < 2000 && amplitude_tree > 10", om, interval.at(j))); //
	  
	  
	  
	  //serie de condition pour verifier que courbe LI est fittable
	  
	  //std::cout << "ref_gain = " << ref_gain_table[om-number] << '\n';
	  //std::cout << "temps = " << interval.at(j) << " - " << interval.at(j+1) << '\n';
	  //std::cout<<"nb entries = "<<spectre->GetEntries()<<endl;
	  if (spectre->GetEntries() < 300) {

	    std::cout << "" << '\n';
	    std::cout << "trop peu d'entries pour l'OM " << om << '\n';
	    std::cout << "" << '\n';
	    pic = j;
	    Khi2 = 0;
	    om_number = om;
	    mean = 0;
	    mean_error = 0;
	    if (om > 711 && run_number < 837) {
	      om_number = om-88;
	    }
	    Result_tree.Fill();
	    delete spectre;
	  }
	  else if (spectre->GetMean() > 120000) {//sature
	    std::cout << "" << '\n';
	    std::cout << "the amplitude sature" << '\n';
	    std::cout << "" << '\n';
	    pic = j;
	    Khi2 = 0;
	    om_number = om;
	    mean = 0;
	    mean_error = 0;
	    if (om > 711 && run_number < 837) {
	      om_number = om-88;
	    }
	    Result_tree.Fill();
	    delete spectre;
	  }
	  else if (spectre->GetMean() < 20){ //si pic trop bas -> charge<0ee

	    std::cout<<"mean ="<<spectre->GetMean();
	    std::cout << "" << '\n';
	    std::cout << "too few charge" << '\n';
	    std::cout << "" << '\n';

	    pic=j;
	    Khi2 = 0;
	    om_number = om;
	    mean = 0;
	    mean_error = 0;
	    std::cout<<" pic "<<j<<"om_number "<<om_number<< " run num "<<run_number <<std::endl;
	    if (om > 711 && run_number < 837) {
	      om_number = om-88;
	    }
	    Result_tree.Fill();
	    delete spectre;
	  }
	  else{ // si trois conditions sont bonnes -> fait fit
	    nevent = spectre->Integral();
	    TF1 *f_Gaus = new TF1("f_Gaus", "gaus(0)", 0, 1200);
	    f_Gaus->SetParNames("N_evt","mean_charge","Sigma");
	    // f_Gaus->SetParameters(25, spectre->GetMean(), 100);
	    f_Gaus->SetParameters(25, spectre->GetMean(), spectre->GetRMS());
	    f_Gaus->SetRange(spectre->GetMean()-400, spectre->GetMean()+400);
	    f_Gaus->Draw("same");
	    //fit puis range puis fit...
	    spectre->Fit(f_Gaus, "RQ0");
	    f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
	    spectre->Fit(f_Gaus, "RQ0");
	    f_Gaus->SetRange(f_Gaus->GetParameter(1)-2.5*f_Gaus->GetParameter(2),f_Gaus->GetParameter(1)+5*f_Gaus->GetParameter(2));
	    spectre->Fit(f_Gaus, "RQ0");
	    pic = j;
	    Khi2 = f_Gaus->GetChisquare()/f_Gaus->GetNDF();
	    om_number = om;
	    constante = (f_Gaus->GetParameter(0));
	    mean = (f_Gaus->GetParameter(1));
	    sigma = (f_Gaus->GetParameter(2));
	    mean_error = sqrt(pow(f_Gaus->GetParError(1)*ref_gain_table[om-number],2) + pow(ref_gain_table_error[om-number]*f_Gaus->GetParameter(1),2) );
	    intensity = intensity_chooser(pic);
	    canvas->cd();
	    spectre->Draw();
	    f_Gaus->Draw("same");
	    //std::cout<<"om"<<om<<endl;
	    if (om < 712) {
	      canvas->SaveAs(Form("fit/fit_Li/amplitude_fit/om_%d/OM_%03d_pic_%d_run_%d.png", om, pic, om_number, run_number));
	    }
	    if (om > 711 && run_number < 837) {
	      canvas->SaveAs(Form("fit/fit_Li/charge_fit_refnocor/om_%d/OM_%03d_pic_%d_run_%d.png", om-88, om -88, pic, run_number));
	      om_number = om-88;
	    }
	    if (om > 711 && run_number > 836) {	      
	      //canvas->SaveAs("test.root");
	      std::cout<<"run num = "<< run_number<<endl;
	      canvas->SaveAs(Form("sortie/ref_Li/fit_Li/om_%d/OM_%03d_pic_%d_run_%d.png", om, om, pic, run_number));
	      
	    }
	    Result_tree.Fill();
	    
	    delete spectre;
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


void file_merger(std::vector<int> run, int ref = 0, string addfile = "") {
  double Amplitude, Amplitude_error, Khi2, time;
  int om_number, run_number, pic;
  string filename = "";
  string cut = "";
  if (ref == 1) filename = "_Ref";
  else cut = "_khi2cut";
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("pic", &pic);
  Result_tree.Branch("Khi2", &Khi2);
  Result_tree.Branch("Amplitude", &Amplitude);
  Result_tree.Branch("Amplitude_error", &Amplitude_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("time", &time);
  std::cout<< "size ="<< run.size();
  for (size_t j = 0; j < run.size(); j++) {
    TFile *file1 = new TFile(Form("sortie/ref_Li/Fit_Ampl%s/Amplitude_Li_%srun_%d.root", filename.c_str(), cut.c_str(), run[j]), "READ");
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

  if (addfile.compare("") != 0) {
    TFile *file1 = new TFile(Form("sortie/ref_Li/%s.root", addfile.c_str()), "READ");
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

    for (int i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      Result_tree.Fill();
    }
    file1->Close();
  }

  TFile *newfile = new TFile(Form("sortie/ref_Li/Fit_Ampl%s/Merged_Fit_%d-%d.root", filename.c_str(),run[0],run[run.size()-1]), "RECREATE");
  std::cout << "file saved : " << Form("sortie/ref_Li/Fit_Ampl%s/Merged_Fit_%d-%d.root", filename.c_str(),run[0],run[run.size()-1])<< '\n';
  newfile->cd();
  Result_tree.Write();
  newfile->Close();

}





void TGrapher(std::string file_name, int n_run) { // pas sur marcher bien

  int om_number, pic;
  double Amplitude, Amplitude_error, time;
  TFile file(Form("sortie/ref_Li/TGraph/TGraph_%s.root", file_name.c_str()), "RECREATE");

  TFile *tree_file = new TFile(Form("sortie/ref_Li/root_file_final/%s.root", file_name.c_str()), "READ");
  TTree* tree = (TTree*)tree_file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("pic",1);
  tree->SetBranchAddress("pic", &pic);
  tree->SetBranchStatus("Amplitude",1);
  tree->SetBranchAddress("Amplitude", &Amplitude);
  tree->SetBranchStatus("Amplitude_error",1);
  tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);

  std::cout << "nrun = " << n_run << '\n';

  TGraphErrors *gain_graph[5][6];
  double yaxis[n_run];
  double yaxis_error[n_run];
  double xaxis[n_run];
  double xaxis_error[n_run];

  int jmin = 0;
  int jmax =712;
  double Ampl_ref, Ampl_ref_error;
  file.cd();
  for (int j = jmin; j < jmax; j++){
    auto canvas = new TCanvas(Form("Allfit_om_%d", j),"",1600,800);
    canvas->Divide(3,2);
    for (int i = 0; i < 6; i++) {
      for (int k = 0; k < n_run; k++) {
        tree->GetEntry((jmax-jmin)*6*k + (j-jmin)*6 + i);// + (j-jmin)*i + i);
        if (k == 0) {
          Ampl_ref = Amplitude;
          Ampl_ref_error = Amplitude_error;
        }
        std::cout << "entry = " << (jmax-jmin)*6*k + i + (j-jmin)*6 << '\n';
        std::cout << "k = " << k << "; i = " << i << "; j = " << j << " and Amplitude = " <<  Amplitude << '\n';
        yaxis[k] = Amplitude/Ampl_ref;
        yaxis_error[k] = Amplitude_error/Ampl_ref + (Ampl_ref_error*Amplitude)/(Ampl_ref*Ampl_ref);
        if (k == 0) {
          yaxis_error[k] = Amplitude_error/Ampl_ref;
        }
        xaxis[k] = time;
        xaxis_error[k] = 0.01;
      }
      std::cout << "" << '\n';
      gain_graph[jmin-jmax][i] = new TGraphErrors(n_run, xaxis, yaxis, xaxis_error, yaxis_error);
      string name = to_string (j);
      if (j > 711) {
        name = namer(j);
      }
      gain_graph[jmin-jmax][i]->SetName(Form("fit_Tl_om_%s_pic_%d", name.c_str(), i));
      gain_graph[jmin-jmax][i]->SetNameTitle(Form("fit_Tl_om_%s_pic_%d", name.c_str(), i), Form("evolution du gain de l'OM %s, pic %d", name.c_str(), i));
      gain_graph[jmin-jmax][i]->GetXaxis()->SetTitle("Temps (h)");
      gain_graph[jmin-jmax][i]->GetYaxis()->SetTitle("Gain(t)/Gain(0)");
      gain_graph[jmin-jmax][i]->GetXaxis()->SetTimeDisplay(1);
      gain_graph[jmin-jmax][i]->SetMarkerColor(2);
      gain_graph[jmin-jmax][i]->SetMarkerStyle(3);
      gain_graph[jmin-jmax][i]->SetMarkerSize(2);
      gain_graph[jmin-jmax][i]->GetYaxis()->SetRangeUser(0.9, 1.1);

      canvas->cd(i+1);
      gain_graph[jmin-jmax][i]->Draw("AP");

      gain_graph[jmin-jmax][i]->Write();
      canvas->Update();

      gain_graph[jmin-jmax][i]->Write();

    }
    canvas->Update();
    canvas->Write();
    delete canvas;
  }



  file.Close();
  return;
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

void Li_eres() { // voir la resolution en energie avec le LI

  int om_number, run_number, events;
  double mean;
  double mean_error, mean_energy;
  double sigma, FWHM;
  int pic =0;
  double Khi2 = 0;

  TFile *file = new TFile("sortie/ERES/ERES_1067.root","RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_energy", &mean_energy);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("FWHM", &FWHM);
  Result_tree.Branch("events", &events);
  Result_tree.Branch("pic", &pic);
  gROOT->cd();

  TFile *energy_file = new TFile ("sortie/Fit_Energy/New_energy.root", "READ");
  TTree* energy_tree = (TTree*)energy_file->Get("Result_tree");

  energy_tree->SetBranchStatus("*",0);
  energy_tree->SetBranchStatus("run_number",1);
  energy_tree->SetBranchAddress("run_number", &run_number);
  energy_tree->SetBranchStatus("om_number",1);
  energy_tree->SetBranchAddress("om_number", &om_number);
  energy_tree->SetBranchStatus("Khi2",1);
  energy_tree->SetBranchAddress("Khi2", &Khi2);
  energy_tree->SetBranchStatus("om_number",1);
  energy_tree->SetBranchAddress("om_number", &om_number);
  energy_tree->SetBranchStatus("mean",1);
  energy_tree->SetBranchAddress("mean", &mean);
  energy_tree->SetBranchStatus("sigma",1);
  energy_tree->SetBranchAddress("sigma", &sigma);
  energy_tree->SetBranchStatus("mean_energy",1);
  energy_tree->SetBranchAddress("mean_energy", &mean_energy);
  energy_tree->SetBranchStatus("pic",1);
  energy_tree->SetBranchAddress("pic", &pic);
  gROOT->cd();

  TH1D* distrib = new TH1D("FWHM at 1 MeV", "FWHM at 1 MeV",100,5,25);

  double energy_selec[520];
  double dif[520];

  for (int i = 0; i < 520; i++) energy_selec[i] = 10;
  for (int i = 0; i < 520; i++) dif[i] = 10;
  for (int i = 0; i < energy_tree->GetEntries(); i++){
    energy_tree->GetEntry(i);
    if (run_number  == 1067) {
      cout << "om : " << om_number << " -> " << mean_energy << endl;
      if (abs(mean_energy-1) < dif[om_number] && mean_energy> 0) {energy_selec[om_number] = mean_energy;}
      }
      dif[om_number] = abs(mean_energy-1);
  }
  for (int i = 0; i < energy_tree->GetEntries(); i++) {
    energy_tree->GetEntry(i);
    if (energy_selec[om_number] == mean_energy && run_number == 1067 &&  235.482*sqrt(mean_energy)*sigma/mean >0 && om_number%13 !=0 && om_number%13 !=12) {
      FWHM = 235.482*sqrt(mean_energy)*sigma/mean;
      distrib->Fill(FWHM);
      Result_tree.Fill();
    }
  }
  distrib->Draw();

  file->cd();
  Result_tree.Write();
  file->Close();



}

void Li_eres_new() {//lui ou Li_eres

  int om_number, run_number;
  double mean;
  double mean_error, mean_energy;
  double sigma, FWHM;
  int pic =0;
  double Khi2 = 0;
  double nevent;
  
  TFile *file = new TFile("ERES/ERES_1068_test.root","RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_energy", &mean_energy);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("FWHM", &FWHM);
  Result_tree.Branch("nevent", &nevent);
  gROOT->cd();

  double gain_corrector[520];
  memset(gain_corrector,0,520*sizeof(double));
  int pic_corrector[520];
  memset(pic_corrector,0,520*sizeof(int));

  TFile *energy_correction_file = new TFile ("ERES/ERES_1067.root", "READ");
  TTree* correction_tree = (TTree*)energy_correction_file->Get("Result_tree");

  correction_tree->SetBranchStatus("*",0);
  correction_tree->SetBranchStatus("mean_energy",1);
  correction_tree->SetBranchAddress("mean_energy", &mean_energy);
  correction_tree->SetBranchStatus("om_number",1);
  correction_tree->SetBranchAddress("om_number", &om_number);
  correction_tree->SetBranchStatus("pic",1);
  correction_tree->SetBranchAddress("pic", &pic);

  for (int i = 0; i < correction_tree->GetEntries(); i++) {
    correction_tree->GetEntry(i);
    gain_corrector[om_number] = mean_energy;
    pic_corrector[om_number] = pic;
  }

  TFile *energy_file = new TFile ("sortie/Fit_Energy/New_energy.root", "READ");
  TTree* energy_tree = (TTree*)energy_file->Get("Result_tree");

  energy_tree->SetBranchStatus("*",0);
  energy_tree->SetBranchStatus("run_number",1);
  energy_tree->SetBranchAddress("run_number", &run_number);
  energy_tree->SetBranchStatus("om_number",1);
  energy_tree->SetBranchAddress("om_number", &om_number);
  energy_tree->SetBranchStatus("Khi2",1);
  energy_tree->SetBranchAddress("Khi2", &Khi2);
  energy_tree->SetBranchStatus("om_number",1);
  energy_tree->SetBranchAddress("om_number", &om_number);
  energy_tree->SetBranchStatus("mean",1);
  energy_tree->SetBranchAddress("mean", &mean);
  energy_tree->SetBranchStatus("sigma",1);
  energy_tree->SetBranchAddress("sigma", &sigma);
  energy_tree->SetBranchStatus("mean_energy",1);
  energy_tree->SetBranchAddress("mean_energy", &mean_energy);
  energy_tree->SetBranchStatus("pic",1);
  energy_tree->SetBranchAddress("pic", &pic);
  energy_tree->SetBranchStatus("nevent",1);
  energy_tree->SetBranchAddress("nevent", &nevent);
  gROOT->cd();

  TH1D* distrib = new TH1D("FWHM at 1 MeV", "FWHM at 1 MeV",100,5,25);

  double energy_selec[520];
  double dif[520];

  for (int i = 0; i < 520; i++) energy_selec[i] = 10;
  for (int i = 0; i < 520; i++) dif[i] = 10;
  for (int i = 0; i < energy_tree->GetEntries(); i++){
    energy_tree->GetEntry(i);
    if (run_number  == 1068) {
      if (abs(gain_corrector[om_number]-1) < dif[om_number] && gain_corrector[om_number]-1> 0) {
        std::cout << abs(gain_corrector[om_number]-1) << " < " <<  dif[om_number] << '\n';
        energy_selec[om_number] = gain_corrector[om_number];
      }

      }
      dif[om_number] = abs(gain_corrector[om_number]-1);
  }

  for (int i = 0; i < energy_tree->GetEntries(); i++) {
    energy_tree->GetEntry(i);
    // if (gain_corrector[om_number] > 0.5 && gain_corrector[om_number] < 1.5 && Khi2 < 2.1 && energy_selec[om_number] == gain_corrector[om_number] && run_number == 1068) {
    if (gain_corrector[om_number] > 0.5 && gain_corrector[om_number] < 2.5  && pic_corrector[om_number] == pic && run_number == 1068 && om_number != 259) {
      FWHM = 235.482*sqrt(gain_corrector[om_number])*sigma/mean;
      distrib->Fill(FWHM);
      Result_tree.Fill();
    }
    else{
      FWHM = 999;
      nevent = 0;
      Result_tree.Fill();
    }
  }
  distrib->Draw();

  file->cd();
  Result_tree.Write();
  file->Close();



}







void Evolution_Li_ref_graph(int n_run, int start, int stop,string wall){//comparaison graph corrige et pas corrige
  gSystem->mkdir(Form("sortie/ref_Li/variation_ref/run_%d_%d",start,stop));
  double amp_scor[5][8][n_run];
  double amp_scor_error[5][8][n_run];
  double amp_cor[5][8][n_run];
  double amp_cor_error[5][8][n_run];
  double Amplitude_error, Amplitude, time;
  int run_number, om_number, pic, Ref_error;
  double time_vec[n_run];

  /*
  TFile file("sortie/ref_Li/Fit_Ampl_Ref/Mergedfitnoncor.root", "READ");
  TTree* tree = (TTree*)file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Amplitude",1);
  tree->SetBranchAddress("Amplitude", &Amplitude);
  tree->SetBranchStatus("Amplitude_error",1);
  tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree->SetBranchStatus("Ref_error",1);
  tree->SetBranchAddress("Ref_error", &Ref_error);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);
  tree->SetBranchStatus("pic",1);
  tree->SetBranchAddress("pic", &pic);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);

  int scompteur = 0;
  int srun = 1049;

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    while (run_number == 1060) {
      i++;
      tree->GetEntry(i);
    }
    if (run_number != srun) {
      scompteur++;
      srun = run_number;
    }
    amp_scor[om_number -712][pic-1][scompteur] = Amplitude;
    amp_scor_error[om_number -712][pic-1][scompteur] = Amplitude_error;
  }
  */

  TFile file_cor(Form("sortie/ref_Li/Fit_Ampl_Ref/Merged_Fit_%d-%d.root",start,stop), "READ");
  TTree* tree_cor = (TTree*)file_cor.Get("Result_tree");
  tree_cor->SetBranchStatus("*",0);
  tree_cor->SetBranchStatus("Amplitude",1);
  tree_cor->SetBranchAddress("Amplitude", &Amplitude);
  tree_cor->SetBranchStatus("Amplitude_error",1);
  tree_cor->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree_cor->SetBranchStatus("run_number",1);
  tree_cor->SetBranchAddress("run_number", &run_number);
  tree_cor->SetBranchStatus("om_number",1);
  tree_cor->SetBranchAddress("om_number", &om_number);
  tree_cor->SetBranchStatus("pic",1);
  tree_cor->SetBranchAddress("pic", &pic);
  tree_cor->SetBranchStatus("time",1);
  tree_cor->SetBranchAddress("time", &time);


  int srun = start;//mettre le premier run ??
  int scompteur = 0;
  
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
    if(pic<=4 && wall=="IT"){
      amp_cor[om_number -712][pic-1][scompteur] = Amplitude;
      std::cout<<"Position "<<om_number -712<<" "<< pic-1 << " "<<scompteur<<std::endl;
      amp_cor_error[om_number -712][pic-1][scompteur] = Amplitude_error;
    }
    else if(pic>4 && wall=="FR"){
      amp_cor[om_number -712][pic-5][scompteur] = Amplitude;
      std::cout<<"Position "<<om_number -712<<" "<< pic-5 << " "<<scompteur<<std::endl;
      amp_cor_error[om_number -712][pic-5][scompteur] = Amplitude_error;
    }
  }
  
   
  double norm_amp_cor[5][5][n_run];
  double norm_amp_cor_error[5][5][n_run];
  double norm_amp_scor[5][5][n_run];
  double norm_amp_scor_error[5][5][n_run];

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 4; j++) {
      for (int l = 0; l < n_run; l++) {
        if (amp_cor[i][j][0] > 0.1 /*&& amp_scor[i][j][0] > 0.1*/) {
          norm_amp_cor[i][j][l] = amp_cor[i][j][l]/amp_cor[i][j][0];
	  //std::cout<<"gain final = "<<norm_amp_cor[i][j][l]<<std::endl;
          norm_amp_cor_error[i][j][l] = sqrt(pow(amp_cor_error[i][j][l]/amp_cor[i][j][0],2) + pow(amp_cor[i][j][l]*amp_cor_error[i][j][0]/pow(amp_cor[i][j][0],2),2));
	  //norm_amp_scor[i][j][l] = amp_scor[i][j][l]/amp_scor[i][j][0];
          //norm_amp_scor_error[i][j][l] = sqrt(pow(amp_scor_error[i][j][l]/amp_scor[i][j][0],2) + pow(amp_scor[i][j][l]*amp_scor_error[i][j][0]/pow(amp_scor[i][j][0],2),2));

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
  multiGraph->SetNameTitle(Form("fit_OM_ref_%d", 712+i), Form("Relative gain LED evolution from the ref OM %d measured with LI spectra", 712+i));
  multiGraph->GetYaxis()->SetTitle("gain evolution");
  multiGraph->GetXaxis()->SetTitle("time");
  double min_x = *std::min_element(xaxis,xaxis+n_run);
  double max_x = *std::max_element(xaxis,xaxis+n_run);

  multiGraph->GetYaxis()->SetRangeUser(0.9,1.1);
  multiGraph->GetXaxis()->SetTimeDisplay(1);    

  
  TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
    for (int j = 0; j < 4; j++) { // pic num
      TCanvas* c = new TCanvas;
      double yaxis[n_run];
      double yaxis_error[n_run];
      //double syaxis[n];
      //double syaxis_error[n];
      //memset (yaxis, 0, n_run*sizeof(double));  //remplissage memoire avec des 0
      //memset (yaxis_error, 0, n_run*sizeof(double));
      //memset (syaxis, 0, n*sizeof(double));
      //memset (syaxis_error, 0, n*sizeof(double));
      for (int l = 0; l < n_run; l++) {
        if (norm_amp_cor[i][j][l] >0.1 /*&& norm_amp_scor[i][j][l] > 0.1*/) {
          yaxis[l] = norm_amp_cor[i][j][l];
          yaxis_error[l] = norm_amp_cor_error[i][j][l];	
          //syaxis[l] = norm_amp_scor[i][j][l];
          //syaxis_error[l] = norm_amp_scor_error[i][j][l];
          cout << "om " << i+712 << " -> pic " << j+1 << " run " << l << " : " << norm_amp_cor[i][j][l] << " +- " << norm_amp_cor_error[i][j][l]  << endl;
        }
      }
      //TGraphErrors *variation_scor = new TGraphErrors(n, xaxis, syaxis, xaxis_error, syaxis_error);
      TGraphErrors *variation_cor = new TGraphErrors(n_run, xaxis, yaxis, xaxis_error, yaxis_error);
      //variation_scor->Draw();
      variation_cor->SetLineColor(j+1);
      legend->AddEntry(variation_cor, Form("pic number %d", j+1), "l");
      variation_cor->Draw("APL"/*"same"*/);
      //variation_cor->GetYaxis()->SetRangeUser(0.9,1.1);
      //variation_cor->GetXaxis()->SetTimeDisplay(1);
      //variation_cor->GetXaxis()->SetRangeUser(0,n_run);
  
  
      //c->SaveAs(Form("sortie/ref_Li/variation_ref/variation_om_%d_pic_%d.png",712+i,j+1));
      multiGraph->Add(variation_cor);
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
  }
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
    for (int j = 0; j < 5; j++) {
      std::cout<<"ancien gain = "<<ref_gain_tab[j]<<std::endl;
        ref_gain_tab[j] = ref_gain_tab[j]/ref_gain_tab_base[j];
	//ref gain tab te donne la variation a corriger pour le run etudie (i)
        cout << ref_gain_tab[j] << " - " << ref_gain_tab_base[j] << endl;
        ref_gain_tab_error[j] = sqrt(pow(ref_gain_tab_error[j]/ref_gain_tab_base[j],2) + pow(ref_gain_tab[j]*ref_gain_tab_base_error[j]/(pow(ref_gain_tab_base[j],2)),2));
    }
    // }
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
  return 0;

  //Passage au calo 
  // std::cout << "" << '\n';
  // std::cout << "START OF THE CALORIMETER FIT" << '\n';
  // std::cout << "" << '\n';
  //
  // double* gain_tab = new double[11];  //correction liumiere par OM ref
  // for (int i = 0; i < run_number.size(); i++) {
  //   Li_corrector(run_number[i], gain_tab, run_number[0]);  //rempli le tableau
  //   std::cout << energy << '\n';
  //   if (energy == false) {
  //     fit_LI_amplitude(run_number[i], gain_tab);
  //   }
  //   if (energy == true) {
  //     fit_LI_energy(run_number[i], gain_tab, energy_run_number[i]);
  //   }
  //   std::cout << "fit ok" << '\n';
  //   Khi2selector(run_number[i]);
  // }
  // if (add == false) {
  //   file_merger(run_number);
  // }
  // else{
  //   file_merger(run_number, 0, file);
  // }
  //
  // TGrapher("Merged_Fit", n_run);



  // amplitude_variation(run_number, time);

  return 0;
}







void Evolution_Li_graph(int n){

  double amp_scor[50][8][n];
  double amp_scor_error[50][8][n];
  double amp_cor[50][8][n];
  double amp_cor_error[50][8][n];
  double omnum[50];

  TFile file("variation_calo/cor/Amplitude_Li.root", "READ");
  double Amplitude_error, Amplitude;
  int run_number, om_number, pic;
  TTree* tree = (TTree*)file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Amplitude",1);
  tree->SetBranchAddress("Amplitude", &Amplitude);
  tree->SetBranchStatus("Amplitude_error",1);
  tree->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);
  tree->SetBranchStatus("pic",1);
  tree->SetBranchAddress("pic", &pic);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);

  int scompteur = 0;
  int srun = 1049;
  int buffer = 0;
  int omcompteur = 0;

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (buffer !=0 && buffer != om_number) omcompteur++;
    if (run_number != srun) {
      scompteur++;
      srun = run_number;
      buffer = 8;
      omcompteur = 0;
      std::cout << buffer << " : " << om_number << " -> " << omcompteur << '\n';
    }
    buffer = om_number;
    while (run_number == 1060) {
      i++;
      tree->GetEntry(i);
    }
    omnum[omcompteur] = om_number;
    amp_scor[omcompteur][pic-1][scompteur] = Amplitude;
    amp_scor_error[omcompteur][pic-1][scompteur] = Amplitude_error;
  }
  TFile file_cor("variation_calo/non_cor/Amplitude_Li.root", "READ");
  TTree* tree_cor = (TTree*)file_cor.Get("Result_tree");
  tree_cor->SetBranchStatus("*",0);
  tree_cor->SetBranchStatus("Amplitude",1);
  tree_cor->SetBranchAddress("Amplitude", &Amplitude);
  tree_cor->SetBranchStatus("Amplitude_error",1);
  tree_cor->SetBranchAddress("Amplitude_error", &Amplitude_error);
  tree_cor->SetBranchStatus("run_number",1);
  tree_cor->SetBranchAddress("run_number", &run_number);
  tree_cor->SetBranchStatus("om_number",1);
  tree_cor->SetBranchAddress("om_number", &om_number);
  tree_cor->SetBranchStatus("pic",1);
  tree_cor->SetBranchAddress("pic", &pic);

  srun = 1049;
  scompteur = 0;
  buffer = 0;
  omcompteur = 0;

  for (int i = 0; i < tree_cor->GetEntries(); i++) {
    tree_cor->GetEntry(i);
    if (buffer !=0 && buffer != om_number) omcompteur++;
    if (run_number != srun) {
      scompteur++;
      srun = run_number;
      buffer = 8;
      omcompteur = 0;
      std::cout << buffer << " : " << om_number << " -> " << omcompteur << '\n';
    }
    buffer = om_number;

    while (run_number == 1060) {
      i++;
      tree_cor->GetEntry(i);
    }
    amp_cor[omcompteur][pic-1][scompteur] = Amplitude;
    amp_cor_error[omcompteur][pic-1][scompteur] = Amplitude_error;
  }

  double norm_amp_cor[50][8][n];
  double norm_amp_cor_error[50][8][n];
  double norm_amp_scor[50][8][n];
  double norm_amp_scor_error[50][8][n];




  for (int i = 0; i < 50; i++) {
    for (int j = 0; j < 8; j++) {
      for (int l = 0; l < 4; l++) {
        if (amp_cor[i][j][0] > 0.1 && amp_scor[i][j][0] > 0.1) {
          norm_amp_cor[i][j][l] = amp_cor[i][j][l]/amp_cor[i][j][0];
          norm_amp_cor_error[i][j][l] = amp_cor_error[i][j][l]/amp_cor[i][j][0];
          norm_amp_scor[i][j][l] = amp_scor[i][j][l]/amp_scor[i][j][0];
          norm_amp_scor_error[i][j][l] = amp_scor_error[i][j][l]/amp_scor[i][j][0];
          // cout << l << " -> om : " << i << " ; pic : " << j << " ; compteur : " << amp_scor[i][j][l] << endl;
        }
      }
    }
  }
  // return;
  double xaxis[n];
  double xaxis_error[n];
  gROOT->cd();
  for (int i = 0; i < 50; i++) {
    for (int j = 0; j < n; j++) {
      TCanvas* c = new TCanvas;
      double yaxis[n];
      double yaxis_error[n];
      double syaxis[n];
      double syaxis_error[n];
      memset (yaxis, 0, n*sizeof(double));
      memset (yaxis_error, 0, n*sizeof(double));
      memset (syaxis, 0, n*sizeof(double));
      memset (syaxis_error, 0, n*sizeof(double));

      for (int l = 0; l < n; l++) {
        if (norm_amp_cor[i][j][l] >0.1 && norm_amp_scor[i][j][l] > 0.1) {
          /* code */
          yaxis[l] = norm_amp_cor[i][j][l];
          yaxis_error[l] = norm_amp_cor_error[i][j][l];
          syaxis[l] = norm_amp_scor[i][j][l];
          syaxis_error[l] = norm_amp_scor_error[i][j][l];
          cout << "om " << i << " -> pic " << j << " run " << l << " : " << norm_amp_cor[i][j][l] << " +- " << norm_amp_cor_error[i][j][l]  << endl;
        }
      }
    
      TGraphErrors *variation_scor = new TGraphErrors(n, xaxis, syaxis, xaxis_error, syaxis_error);
      TGraphErrors *variation_cor = new TGraphErrors(n, xaxis, yaxis, xaxis_error, yaxis_error);

      variation_scor->Draw();
      variation_scor->SetLineColor(kBlue);
      variation_cor->Draw("same");
      variation_scor->GetYaxis()->SetRangeUser(0.9,1.1);
      c->SaveAs(Form("variation_calo/png_graph/variation_om_%f_pic_%d.png",omnum[i],j));

      variation_scor->SetTitle("test");
      variation_scor->SetName("test");
      variation_scor->SetNameTitle("test");
      TFile newfile(Form("sortie/variation_om_%f_pic_%d.root", omnum[i],j),"RECREATE");

      newfile.cd();
      variation_scor->Write();
      variation_cor->Write();
      newfile.Close();
    }
  }
}





void Evolution_Bi_graph(){

  int n = 3;
  double amp_cor[520][3];
  double amp_cor_error[520][3];


  TFile file("../Bi_selection/Bi_fit/fused_bi.root", "READ");
  double mean_error, mean;
  int om_number;
  TTree* tree = (TTree*)file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("mean",1);
  tree->SetBranchAddress("mean", &mean);
  tree->SetBranchStatus("mean_error",1);
  tree->SetBranchAddress("mean_error", &mean_error);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);

  int scompteur = -1;

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (om_number == 0) scompteur++;
    amp_cor[om_number][scompteur] = mean;
    amp_cor_error[om_number][scompteur] = mean_error;
  }

  double norm_amp_cor[520][3];
  double norm_amp_cor_error[520][3];




  for (int i = 0; i < 520; i++) {
    for (int l = 0; l < 3; l++) {
      if (amp_cor[i][0] > 0.1) {
        norm_amp_cor[i][l] = amp_cor[i][l]/amp_cor[i][0];
        norm_amp_cor_error[i][l] = amp_cor_error[i][l]/amp_cor[i][0];
        // cout << l << " -> om : " << i << " ; pic : " << j << " ; compteur : " << amp_scor[i][j][l] << endl;
      }
    }

  }
  // return;
  double xaxis[3] = {0,2,4};
  double xaxis_error[3] = {0,0,0};
  gROOT->cd();
  for (int i = 0; i < 520; i++) {
    TCanvas* c = new TCanvas;
    double yaxis[3];
    double yaxis_error[3];
    memset (yaxis, 0, 3*sizeof(double));
    memset (yaxis_error, 0, 3*sizeof(double));

    for (int l = 0; l < 3; l++) {
      if (norm_amp_cor[i][l] >0.1) {
        yaxis[l] = norm_amp_cor[i][l];
        yaxis_error[l] = norm_amp_cor_error[i][l];
      }
    }

    TGraphErrors *variation_cor = new TGraphErrors(n, xaxis, yaxis, xaxis_error, yaxis_error);

    variation_cor->Draw();
    variation_cor->GetYaxis()->SetRangeUser(0.9,1.1);
    c->SaveAs(Form("variation_calo/bi_png/variation_om_%d.png", i));

    variation_cor->SetTitle("test");
    variation_cor->SetName("test");
    variation_cor->SetNameTitle("test");
    TFile newfile(Form("sortie/variation_om_%d.root", i),"RECREATE");

    newfile.cd();
    variation_cor->Write();
    newfile.Close();
  }
}

void energy_selec(int *energy_selec) {// ici selection est mauvaise -> on veut que tous les OMs qui ont 200 mV = 1MeV    -> pic le plus proche de 200 mV est celui de l'etalonnage -> a changer

  TFile file("variation_calo/cor/Amplitude_Li.root", "READ");
  double Amplitude;
  int om_number,run_number, pic;
  TTree* tree = (TTree*)file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Amplitude",1);
  tree->SetBranchAddress("Amplitude", &Amplitude);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);
  tree->SetBranchStatus("pic",1);
  tree->SetBranchAddress("pic", &pic);
  double dif[520];

  for (int i = 0; i < 520; i++) energy_selec[i] = 10;
  for (int i = 0; i < 520; i++) dif[i] = 1000;
  for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if (run_number  == 1049) {
      if (abs(Amplitude-200) < dif[om_number]) {
        // std::cout << abs(gain_corrector[om_number]-1) << " < " <<  dif[om_number] << '\n';
        energy_selec[om_number] = pic;
        dif[om_number] = Amplitude-200;
      }
    }
  }
}

void comp_bi_li() {// permet d'obtenir graph p101

  int energy_sele[520];
  energy_selec(energy_sele);
  TH1D* diff_li = new TH1D("diff_li","diff_li",100, 0,10);
  TH1D* diff_sli = new TH1D("sdiff_li","sdiff_li",100, 0,10);
  TH1D* diff_li_error = new TH1D("diff_li_error","diff_li_error",100, 0,10);
  TH1D* diff_sli_error = new TH1D("sdiff_li_error","sdiff_li_error",100, 0,10);
  for (int i = 8; i < 130; i++) {
    if (bundle_number(i) == 1 && i%13 != 12 && i!= 128) {

      TFile *file = new TFile(Form("variation_calo/bi_root/variation_om_%d.root", i), "READ");
      gROOT->cd();
      TGraphErrors *bi_graph = (TGraphErrors*)file->Get("test");

      TFile *file2 = new TFile(Form("variation_calo/root_graph/variation_om_%d.000000_pic_%d.root",i, energy_sele[i]), "READ");
      gROOT->cd();
      TGraphErrors *li_graph = (TGraphErrors*)file2->Get("Graph");
      TGraphErrors *sli_graph = (TGraphErrors*)file2->Get("test");

      // TCanvas* c = new TCanvas;
      // li_graph->Draw();
      // li_graph->SetLineColor(kBlue);
      // li_graph->SetLineWidth(2);
      // sli_graph->Draw("same");
      // sli_graph->SetLineStyle(2);
      // sli_graph->SetLineWidth(2);
      // bi_graph->Draw("same");
      // bi_graph->SetLineColor(kGreen);
      // bi_graph->SetLineWidth(2);
      // li_graph->GetYaxis()->SetRangeUser(0.9,1.1);
      // c->SaveAs(Form("variation_calo/comp_graph/comp_graph_om_%d.png",i));
      for (int j = 1; j < 4; j++) {
        if (li_graph->GetPointY(j) > 0 && bi_graph->GetPointY(j) > 0) {
          double diff = 100*sqrt(pow((li_graph->GetPointY(j)-bi_graph->GetPointY(j))/bi_graph->GetPointY(j),2));
          double sdiff = 100*sqrt(pow((sli_graph->GetPointY(j)-bi_graph->GetPointY(j))/bi_graph->GetPointY(j),2));
          double diff_error = 100*sqrt(pow(li_graph->GetErrorY(j)/bi_graph->GetPointY(j),2) + pow(li_graph->GetPointY(j)*bi_graph->GetErrorY(j)/pow(bi_graph->GetPointY(j),2),2));
          double sdiff_error = 100*sqrt(pow(sli_graph->GetErrorY(j)/bi_graph->GetPointY(j),2) + pow(sli_graph->GetPointY(j)*bi_graph->GetErrorY(j)/pow(bi_graph->GetPointY(j),2),2));
          diff_li->Fill(diff);
          diff_li_error->Fill(diff_error);
          cout << diff << endl;
          diff_sli->Fill(sdiff);
          diff_sli_error->Fill(sdiff_error);
        }
      }
      // return;
      // delete c;
    }
  }

  diff_li_error->Draw();
  diff_sli_error->Draw("sames");

}

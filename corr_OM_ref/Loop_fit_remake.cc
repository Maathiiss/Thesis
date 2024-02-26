#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooRealSumPdf.h"
#include "RooParamHistFunc.h"
#include "RooHistConstraint.h"
//#include "RooBinSamplingPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooRealSumPdf.h"
#include "RooParamHistFunc.h"
#include "RooHistConstraint.h"
//#include "RooBinSamplingPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooAddPdf.h"
#include "RooMinimizer.h"
//#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
using namespace RooFit;
#define _USE_MATH_DEFINES
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
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TParameter.h>
using namespace std;

RooRealVar charge_tree("charge_tree","charge_tree",0,2e5);

TH2F* charge_spectre = NULL;
TH2F* charge_spectre_template = NULL;

const int gain_n_bin = 10001; // Precision au 10 000 (me
const double gain_bin_min = 0.5;
const double gain_bin_max = 1.5;
const double gain_bin_width = (gain_bin_max-gain_bin_min)/(gain_n_bin-1);
double time_pdf = 0;

void Load_spectre(int run_number){
  TFile *file = new TFile(Form("entree/root/histo_ref_%d.root", run_number), "READ");
  gROOT->cd();
  charge_spectre = (TH2F*)file->Get("histo_pm_charge");
  return;
}

TH1D* spectre_charge(int om_number){
  TH1D* spectre_charge = charge_spectre->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  // spectre_charge->Rebin(4);
  return spectre_charge;

}

TH2F* spectre_charge_full_template(int om_number, int run_number_pdf){
  TFile *file = new TFile(Form("entree/Modele/Modele_OM_%d_%d.root", om_number, run_number_pdf), "READ");
  //a ajouter
  TParameter<double> *param1 = new TParameter<double>("start_time", time_pdf);
  param1 = (TParameter<double>*)(file->Get("start_time"));
  time_pdf = param1->GetVal();
  gROOT->cd();
  
  charge_spectre_template = (TH2F*)file->Get(Form("Modele_Ref_OM_%d",om_number));
  return charge_spectre_template;
}

double roofitter(TH1D* modele, TH1D* spectre_om, int om_number, double *rootab, int run_number, double gain, double compteur)
{
  using namespace RooFit;

  RooRealVar x("x", "x", 0, 200000);
  x.setBins(1024);
  RooDataHist Tl("Tl", "Tl", x, Import(*modele));

  RooHistPdf modele_pdf ("modele_pdf", "", x, Tl);//construire une fonction de distribution a partir d'un histo
  RooDataHist spectre_data("spectre_data", "spectre_data", x, Import(*spectre_om));

  RooChi2Var RooChi2("Chi2", "Chi2", modele_pdf, spectre_data, Range(compteur, 80000), DataError(RooAbsData::Poisson));   // create the variance
  RooMinimizer *miniChi = new RooMinimizer(RooChi2);

  double Chi2 = RooChi2.getVal();
  // double Chi2 = RooChi2.getVal()/(1024. - 1);
  std::cout << "Chi2 = " << RooChi2.getVal() << '\n';
  TCanvas* can = new TCanvas;
  can->cd();
  auto frame = x.frame(Title("Fit gain"));

  spectre_data.plotOn(frame,DataError(RooAbsData::SumW2), DrawOption("P"));

  modele_pdf.plotOn(frame);
  //spectre_data.plotOn(frame, DataError(RooAbsData::SumW2), DrawOption("P"));
  // Plot model components
  //std::cout<<"oui\n";
  //modele_pdf.plotOn(frame, LineColor(kRed), Name("sum_curve"), Range(3500, 200000));
  //std::cout<<"non\n";
  double Tl_int = x.getVal();
  double Tl_int_error = x.getError();
  frame->GetYaxis()->SetTitle("n events");
  frame->GetXaxis()->SetTitle("charge (u.a)");
  frame ->Draw();
  
  //frame->SaveAs("graph_test.root");
  can->SetLogy();
  frame->GetYaxis()->SetRangeUser(0.01, 5e6);
  //frame->GetYaxis()->SetRangeUser(0,7e4);
  //frame->GetXaxis()->SetRangeUser(0,30e3);
  //can->SaveAs("test.root");

  rootab[1] = Tl_int;
  rootab[2] = Tl_int_error;
  //std::cout<<" 1 :"<<rootab[1]<<" 2 : "<<rootab[2]<<"\n";
  // rootab[0] = RooChi2.getVal()/(1024. - 1);
  rootab[0] =Chi2;
  TLatex l;
  l.SetTextFont(40);
  l.DrawLatex(20000, 40000, Form("Khi2 = %.2f", Chi2));
  //l.DrawLatex(90000, 80, Form("Khi2 = %.2f", Chi2));
  // return *rootab;
  // if (Chi2 < 10000) {
  //gSystem->mkdir(Form("sortie/om_%d/run_%d", om_number, run_number));
  //can->SaveAs(Form("sortie/om_%d/run_%d/fit_run_%d_om_%d_gain_%f.png", om_number, run_number, run_number, om_number, gain));
  //gSystem->mkdir(Form("/home/granjon/Documents/these/Analyse/corr_OM_ref/sortie/graph interessant/Avec pic alpha/gif/om_%d", om_number));
  //can->SaveAs(Form("/home/granjon/Documents/these/Analyse/corr_OM_ref/sortie/graph interessant/Avec pic alpha/gif/om_%d/fit_run_%d_om_%d_gain_%f.png",om_number,run_number, om_number, gain));

  // } 


  delete miniChi;
  delete can;
  delete frame;
  // delete miniLog;
  return *rootab;

}

void Fit_Ref(int run_number, int run_number_pdf) {
  Load_spectre(run_number);
  TH1::SetDefaultSumw2();

  double Chi2, gain, gain_error, ndf, time;
  int om_number, compteur;
  double* rootab = new double[3];

  TH1D* modele = NULL;
  TH2F* TH2modele = NULL;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error", &gain_error);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("nsrun", &compteur);
  Result_tree.Branch("ndf", &ndf);
  Result_tree.Branch("time", &time);

  for (om_number = 712; om_number < 717; om_number++) {

    TH2modele = spectre_charge_full_template(om_number,run_number_pdf);
    //TH2modele->SaveAs("TH2modele.root");
    // for (compteur = 0; compteur < 30; compteur+=1) {
    TH1D* spectre_om = NULL;
    spectre_om = spectre_charge(om_number);
    TFile *file = new TFile(Form("entree/root/histo_ref_%d.root", run_number), "READ");
    // spectre_om = (TH1D*)file->Get("histo_pm_charge");

    time = 0;
    TParameter<double> *param = new TParameter<double>("start_time", time);
    param = (TParameter<double>*)(file->Get("start_time"));
    time = param->GetVal();



    gROOT->cd();

    string t = Form("sortie/om_%d", om_number, run_number);
    if (mkdir(t.c_str(), 0777) == -1){
      cout << "directory already exist" << endl;
    }
    else{
      gSystem->mkdir(Form("sortie/om_%d",om_number));      
      cout << "Directory created" << endl;
    }
    //std::cout<<"cicicicicic"<<endl;
    double gainmin = 0;
    double Chi2min = 1000000;
    for (int gain_count = 0; gain_count <gain_n_bin; gain_count+=100) {
      double compteur = 0;
      std::cout<<compteur<<endl;
      gain = (gain_bin_min + gain_bin_width*(gain_count-1));
      // if (gain < 1.01 && gain > 0.99) {
      modele = TH2modele->ProjectionY("modele", gain_count, gain_count);
      
      TCanvas* canv = new TCanvas;
      canv->cd();
      modele->Draw();
      //modele->SaveAs("test.root");
      spectre_om->Draw("same==========");
      //spectre_om->SaveAs("test_om.root");
      //canv->SaveAs("canv.png");
      //modele->SaveAs("test_2.root");
      
      for (int i = 0; i < 25; i++) {  //mettre les 0 au debut
        modele->SetBinContent(i,0);
        modele->SetBinError(i, 0);
        spectre_om->SetBinContent(i,0);
        spectre_om->SetBinError(i, 0);
      }
      
      /*
      for (int i = 100; i < 1000; i++) {
        modele->SetBinContent(i,0);
        modele->SetBinError(i, 0);
        spectre_om->SetBinContent(i,0);
        spectre_om->SetBinError(i, 0);
	}*/
      // if (om_number == 716) {
      /*
      int k = 10;
      //std::cout << " hauteur bin  :"<<k<<"  = " <<spectre_om->Integral(k, k) << '\n';
      while (spectre_om->Integral(k, k) > 11000) {//enlever effets saturation ? Si hauteur du bin k est >11000 -> mis a 0
        // for (int i = 0; i < 17; i++) {
	modele->SetBinContent(k,0);
	modele->SetBinError(k, 0);
	spectre_om->SetBinContent(k,0);
	spectre_om->SetBinError(k, 0);
	k++;
	compteur = k;
	}*/
      std::cout << "/* message */" << '\n';
      // }
      spectre_om->Draw();
      modele->Draw("same");
      // spectre_om->Scale(1./spectre_om->Integral());
      // modele->Scale(1./modele->Integral());
      // return;
    
      roofitter(modele, spectre_om, om_number, rootab, run_number, gain, compteur);
      //premiere recherche grossiere de chi2
      Chi2 = rootab[0];
      if (Chi2min > Chi2){
        gainmin = gain_count;
        Chi2min = Chi2;
      }
      modele->Reset();

    }
    // gainmin = 1000;



    
    //recherche plus fine de chi2
    for (int gain_count = gainmin - 400; gain_count < gainmin + 400; gain_count+=10) {
      gain = (gain_bin_min + gain_bin_width*(gain_count-1));
      // if (gain < 1.01 && gain > 0.99) {
      modele = TH2modele->ProjectionY("modele", gain_count, gain_count);
      // modele->Draw();
      // spectre_om->Draw();
      
      for (int i = 0; i < 25; i++) {
        modele->SetBinContent(i,0);
        modele->SetBinError(i, 0);
        spectre_om->SetBinContent(i,0);
        spectre_om->SetBinError(i, 0);
      }
      
      //modifier range
      /*
      for (int i = 100; i < 1000; i++) {
        modele->SetBinContent(i,0);
        modele->SetBinError(i, 0);
        spectre_om->SetBinContent(i,0);
        spectre_om->SetBinError(i, 0);
	}*/
  
      // if (om_number == 716) {
      /*
      int k = 10;
      while (spectre_om->Integral(k, k) > 11000) {
        // for (int i = 0; i < 17; i++) {
        modele->SetBinContent(k,0);
        modele->SetBinError(k, 0);
        spectre_om->SetBinContent(k,0);
        spectre_om->SetBinError(k, 0);
        k++;
	}*/
      // }

      int compte = 0;
      for (int i = 0; i < 1000; i++) {
	if (spectre_om->GetBinContent(i) > 0) {
	  compte++;
	}
      }
      ndf = compte - 1;
      int k =0;
      roofitter(modele, spectre_om, om_number, rootab, run_number, gain, k);

      Chi2 = rootab[0];
      modele->Reset();
      Result_tree.Fill();
      // }
      delete modele;
    }

    delete spectre_om;
    // }

  }
  TFile new_file(Form("sortie/root_file_final/Fit_Ref_%d.root", run_number), "RECREATE");
  // TFile new_file(Form("root/Complete_Fit/Fit_Ref_%d.root", run_number), "RECREATE");
  new_file.cd();
  Result_tree.Write();
  new_file.Close();
  std::cout << "file " << Form("Fit_Ref_%d.root", run_number) << " saved" << '\n';
  return;

}


void minerror_calculator(string file_name, int run_number) {
  TFile file(Form("sortie/root_file_final/%s.root", file_name.c_str()), "READ");
  gROOT->cd();
  double PolChi2, Chi2, gain, gainmin, chi2min, error_plus, error_moins, minChi2, mingain_tree, time, ndf;
  int om_number, nsrun;
  int i;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &i);
  Result_tree.Branch("PolChi2", &PolChi2);
  Result_tree.Branch("gainmin", &gainmin);
  Result_tree.Branch("chi2min", &chi2min);
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("error_moins", &error_moins);
  Result_tree.Branch("error_plus", &error_plus);
  Result_tree.Branch("mingain_tree", &mingain_tree);
  Result_tree.Branch("ndf", &ndf);
  Result_tree.Branch("time", &time);

  TTree* tree = (TTree*)file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("Chi2",1);
  tree->SetBranchAddress("Chi2", &Chi2);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);
  tree->SetBranchStatus("ndf",1);
  tree->SetBranchAddress("ndf", &ndf);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);


  TFile newfile(Form("sortie/root_file_final/%s_mini.root", file_name.c_str()), "RECREATE");
  gROOT->cd();

  for (i = 712; i < 717; i++) {
    minChi2 = 100000;
    for (int j = 0; j < tree->GetEntries(); j++) {
      tree->GetEntry(j);
      if (om_number == i) {
        if (minChi2 > Chi2) {
          minChi2 = Chi2;
          gainmin = gain;
        }
      }
    }

    // ndf = ndf_counter(i, run_number)-1;
    mingain_tree = gainmin;
    TF1* poly = new TF1 ("poly", "pol2", gainmin - 0.02, gainmin + 0.02);
    std::cout << "gainmin = " << gainmin << '\n';

    //TH2D *Chi2gain = new TH2D ("Chi2gain", "Chi2gain", 1000, gainmin - 0.1, gainmin + 0.1, 10000, 0, 100000);
    TGraph *Chi2gain = new TGraph();
    //TH2D *Chi2gain = new TH2D ("Chi2gain", "Chi2gain", 1000, gainmin - 0.3, gainmin + 0.3, 10000, 0, 100000);
    double min=0;
    double max = 0;
    for (int j = 0; j < tree->GetEntries(); j++) {
      tree->GetEntry(j);
      if (om_number == i) {
	//Chi2gain->Fill(gain,Chi2);
	Chi2gain->SetPoint(j,gain,Chi2);        
      }
    }
 
    TCanvas* can = new TCanvas;
    can->cd();
    Chi2gain->GetXaxis()->SetRangeUser(0.9,1.1);
    Chi2gain->GetXaxis()->SetTitle("gain");
    Chi2gain->GetYaxis()->SetTitle("Chi2");
    Chi2gain->SetTitle(Form("Fit gain in function of the Chi2 for the OM %d",i));
    Chi2gain->Draw("AP");
    poly->Draw("same");
    Chi2gain->Fit(poly, "RQ");
    PolChi2 = poly->GetChisquare();
    gainmin = poly->GetMinimumX();
    chi2min = poly->GetMinimum();
    error_plus = poly->GetX(poly->GetMinimum()+1, poly->GetMinimumX(), poly->GetMinimumX() + 0.05) - poly->GetMinimumX();
    error_moins = poly->GetMinimumX() - poly->GetX(poly->GetMinimum()+1, poly->GetMinimumX() - 0.05, poly->GetMinimumX());
    std::cout << "min Chi2 = " << poly->GetMinimum() << " min gain = " << poly->GetMinimumX() << " + " << error_plus << " - " << error_moins << '\n';
    std::cout << "Chi2 poly = " << PolChi2 << '\n';
    Chi2gain->GetYaxis()->SetRangeUser(0, 10000);
    // return;

    //can->SaveAs(Form("sortie/om_%d/run_%d/fit_Chi2_gain.root", i, run_number));
    can->SaveAs(Form("sortie/om_%d/fit_Chi2_gain_run_%d.png", i, run_number));
    std::cout<<Form("sortie/om_%d/run_%d/fit_Chi2_gain.png", i, run_number)<<"\n";
    
    Result_tree.Fill();
    delete Chi2gain;
    delete poly;
  }

  newfile.cd();
  Result_tree.Write();
  newfile.Close();

}


void file_merger(std::vector<int> run_number, int run_number_pdf, string previous_file_s = "") {
  
  TFile file(Form("sortie/root_file_final/Fit_Ref_%d-%d.root", run_number.at(0), run_number.at(run_number.size()-1)), "RECREATE");
  // TFile file(Form("root/Merged_Fit/Fit_Ref_%d-%d.root", 736, 836), "RECREATE");
  double Chi2, gain, gain_error_plus, gain_error_moins, ndf;
  int om_number, int_run;
  double time;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error_plus", &gain_error_plus);
  Result_tree.Branch("gain_error_moins", &gain_error_moins);
  Result_tree.Branch("run_number", &int_run);
  Result_tree.Branch("time", &time);
  Result_tree.Branch("ndf", &ndf);

  if (previous_file_s.compare("") != 0){
    TFile previous_file(Form("sortie/root_file_final/%s.root", previous_file_s.c_str()), "READ");
    TTree* previous_tree = (TTree*)previous_file.Get("Result_tree");
    previous_tree->SetBranchStatus("*",0);
    previous_tree->SetBranchStatus("om_number",1);
    previous_tree->SetBranchAddress("om_number", &om_number);
    previous_tree->SetBranchStatus("Chi2",1);
    previous_tree->SetBranchAddress("Chi2", &Chi2);
    previous_tree->SetBranchStatus("gain",1);
    previous_tree->SetBranchAddress("gain", &gain);
    previous_tree->SetBranchStatus("error_moins",1);
    previous_tree->SetBranchAddress("error_moins", &gain_error_moins);
    previous_tree->SetBranchStatus("error_plus",1);
    previous_tree->SetBranchAddress("error_plus", &gain_error_plus);
    previous_tree->SetBranchStatus("run_number",1);
    previous_tree->SetBranchAddress("run_number", &int_run);
    previous_tree->SetBranchStatus("ndf",1);
    previous_tree->SetBranchAddress("ndf", &ndf);
    previous_tree->SetBranchStatus("time",1);
    previous_tree->SetBranchAddress("time", &time);
    for (double i = 0; i < previous_tree->GetEntries(); i++) {
      previous_tree->GetEntry(i);
      Result_tree.Fill();
    }
  }
  else {
    for (int i = 712; i < 717; i++) {
      om_number = i;
      Chi2 = 0;
      gain = 1;
      gain_error_moins = 0/*0.0009755*/;
      gain_error_plus = 0/*0.0009755*/;
      int_run = run_number_pdf;
      time = time_pdf;
      Result_tree.Fill();
    }
  }

  for (size_t i = 0; i < run_number.size(); i++) {
    TFile tree_file(Form("sortie/root_file_final/Fit_Ref_%d_mini.root", run_number.at(i)), "READ");
    std::cout<<"ici "<<Form("sortie/root_file_final/Fit_Ref_%d_mini.root", run_number.at(i))<<"\n";
    std::cout << run_number.at(i) << '\n';
    TTree* tree = (TTree*)tree_file.Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchStatus("chi2min",1);
    tree->SetBranchAddress("chi2min", &Chi2);
    tree->SetBranchStatus("gainmin",1);
    tree->SetBranchAddress("gainmin", &gain);
    tree->SetBranchStatus("error_moins",1);
    tree->SetBranchAddress("error_moins", &gain_error_moins);
    tree->SetBranchStatus("error_plus",1);
    tree->SetBranchAddress("error_plus", &gain_error_plus);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &int_run);
    tree->SetBranchStatus("ndf",1);
    tree->SetBranchAddress("ndf", &ndf);
    tree->SetBranchStatus("time",1);
    tree->SetBranchAddress("time", &time);

    int_run = run_number.at(i);
    std::cout << "ok" << i+1 << '\n';
    for (int j = 0; j < 5; j++) {
      tree->GetEntry(j);
      Result_tree.Fill();
    }
  }
  file.cd();
  Result_tree.Write();
  file.Close();
}

void TGrapher(std::string file_name, int n_run) {
  std::cout << "ok" << '\n';
  TFile file(Form("sortie/root_file_final/TGraph_%s.root", file_name.c_str()), "RECREATE");

  TFile tree_file(Form("sortie/root_file_final/%s.root", file_name.c_str()), "READ");
  double time;
  int om_number, run_number;
  double gain;
  double gain_error_moins, gain_error_plus;
  TTree* tree = (TTree*)tree_file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);
  tree->SetBranchStatus("gain_error_moins",1);
  tree->SetBranchAddress("gain_error_moins", &gain_error_moins);
  tree->SetBranchStatus("gain_error_plus",1);
  tree->SetBranchAddress("gain_error_plus", &gain_error_plus);
  tree->SetBranchStatus("run_number",1);
  tree->SetBranchAddress("run_number", &run_number);
  tree->SetBranchStatus("time",1);
  tree->SetBranchAddress("time", &time);

  double yaxis[n_run];
  double yaxis_error_moins[n_run];
  double yaxis_error_plus[n_run];
  double xaxis[n_run];
  double xaxis_error_moins[n_run];
  double xaxis_error_plus[n_run];
  auto canvas = new TCanvas("Allfit","",1600,800);
  canvas->Divide(3,2);
  TGraphAsymmErrors *gain_graph[5]; //(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);

  file.cd();
  int compteur = 0;
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < n_run; j++){
      //tree->GetEntry(j);
      tree->GetEntry(i+j*5);
      compteur++;
       std::cout << "compteur " << i+j*5 << " and om = " << om_number << " and run = " << run_number << '\n';

      yaxis[j] = gain;
      std::cout << gain << '\n';
      yaxis_error_moins[j] = gain_error_moins;
      yaxis_error_plus[j] = gain_error_plus;
      xaxis[j] = time;
      xaxis_error_plus[j] = 0.00001;
      xaxis_error_moins[j] = 0.00001;
    }
     for (size_t k = 0; k < n_run; k++) {
       std::cout << yaxis[k] << '\n';
     }
    gain_graph[i] = new TGraphAsymmErrors(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);

    gain_graph[i]->SetName(Form("fit_OM_ref_%d", om_number));
    gain_graph[i]->SetNameTitle(Form("fit_OM_ref_%d", om_number), Form("Relative gain evolution of the ref OM %d measured with background spectra", om_number));
    gain_graph[i]->GetXaxis()->SetTimeDisplay(1);
    gain_graph[i]->GetXaxis()->SetTitle("Time");
    gain_graph[i]->GetYaxis()->SetTitle("Gain evolution");
    gain_graph[i]->GetYaxis()->SetRangeUser(0.9, 1.1);
    gain_graph[i]->GetYaxis()->SetTitleOffset(0.9);
    gain_graph[i]->SetMarkerColor(2);
    gain_graph[i]->SetMarkerStyle(5);
    gain_graph[i]->SetMarkerSize(2);
    canvas->cd(i+1);
    gain_graph[i]->Draw("AP");

    gain_graph[i]->Write();
    canvas->Update();
    // if (i == 4) {
    //   canvas->Write();
    // }
  }
  canvas->Update();
  canvas->Write();
  file.Close();

}



void comparator(int om){
  TFile file("/sortie/root_file_final/Fit_Ref_1050.root", "READ");
  // TFile file("variation_gain/variation_gain_716.root", "READ");
  gROOT->cd();
  double Chi2, gain, mingain, minChi2, gainmin, measured_gain, fitted_gain;
  double measured_min[29];
  int om_number, nsrun;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("measuredgain", &measured_gain);
  Result_tree.Branch("fittedgain", &fitted_gain);
  Result_tree.Branch("nsrun", &nsrun);

  TTree* tree = (TTree*)file.Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("Chi2",1);
  tree->SetBranchAddress("Chi2", &Chi2);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);
  tree->SetBranchStatus("nsrun",1);
  tree->SetBranchAddress("nsrun", &nsrun);


  TH2D *Measuredgain = new TH2D ("Measuredgain", "", 1000, 0, 30, 1000, 0.99, 1.01);
  for (int i = 0; i < 30; i+=1) {
    minChi2 = 10000;
    for (int j = 0; j < tree->GetEntries(); j++) {
      tree->GetEntry(j);
      if (nsrun == i && om_number == om) {
        if (minChi2 > Chi2) {
          minChi2 = Chi2;
          mingain = gain;
        }
      }
      measured_min[i] = mingain;
    }
    Measuredgain->Fill(i, measured_min[i]);
    std::cout << i << '\n';
  }
  // TFile file2("root/Fittedgain.root", "READ");
  TFile file2(Form("sortie/root_file_final/Fit_Ref_1050_%d_mini.root", om), "READ");
  gROOT->cd();
  TTree* tree2 = (TTree*)file2.Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number);
  tree2->SetBranchStatus("gainmin",1);
  tree2->SetBranchAddress("gainmin", &gainmin);
  tree2->SetBranchStatus("nsrun",1);
  tree2->SetBranchAddress("nsrun", &nsrun);
  TH2D *Tcomparator = new TH2D ("Chi2gain", "", 1000, 0, 30, 1000, -0.01, 0.01);
  TH2D *Fittedgain = new TH2D ("Fittedgain", "", 1000, 0, 30, 1000,  0.99, 1.01);
  for (int i = 0; i < tree2->GetEntries(); i++) {
    tree2->GetEntry(i);
    Tcomparator->Fill(nsrun, measured_min[i]-gainmin);
    Fittedgain->Fill(nsrun, gainmin);

    measured_gain = measured_min[i];
    fitted_gain = gainmin;
    Result_tree.Fill();

  }
  Tcomparator->Draw();
  Tcomparator->SetMarkerStyle(2);
  // TCanvas* c = new TCanvas;
  Measuredgain->Draw();
  Measuredgain->SetMarkerStyle(2);
  Measuredgain->SetMarkerColor(4);
  Fittedgain->Draw("same");
  Fittedgain->SetMarkerStyle(4);
  Fittedgain->SetMarkerColor(2);

  TFile newfile(Form("sortie/root_file_final/comparator_%d.root", om), "RECREATE");
  newfile.cd();
  Tcomparator->Write();
  Measuredgain->Write();
  Fittedgain->Write();
  Result_tree.Write();
  newfile.Close();
  return;
}


int main(int argc, char const *argv[]) {
  int n_run, run, run_number_pdf;
  std::vector<int> run_number, time, ref_run_number, ref_time;
  int compteur = 0;
  std::string file;
  bool add = false;

  for(int i = 0; i<argc; i++){
    if (std::string(argv[i]) == "-add" ) {
      // if (std::string(argv[i]).compare("-add") == 0 ) {
      // file = argv[i+1];
      // std::cout << file << '\n';
      add = true;
    }
  }

  if (add == true) {                  ///// Add run to the other runs
    string old_run;
    std::cout << "To what file do you want to add the new run(s)?" << '\n';
    std::cin >> old_run;
    std::cout << "What is your run of pdf reference?" << '\n';
    std::cin >> run_number_pdf;

    std::cout << "How many run do you want to add?" << '\n';
    std::cin >> n_run;
    std::cout << "Write the run(s) you want" << '\n';
    while (compteur < n_run && cin >> run) {
      run_number.push_back(run);
      compteur++;
      if (compteur < n_run) {
	std::cout << "Write the runs you want" << '\n';
      }
    }
    std::cout << "Code start running" << '\n';

    for (int i = 0; i < n_run; i++) {
      Fit_Ref(run_number.at(i),run_number_pdf);
      minerror_calculator(Form("Fit_Ref_%d", run_number.at(i)), run_number.at(i));
    }
    std::cout << "Fit_Ref ok and minerror ok" << '\n';

    file_merger(run_number, run_number_pdf, old_run);
    std::cout << "file_merger ok" << '\n';

    TGrapher(Form("Fit_ref_716-%d", run_number.at(run_number.size()-1)), n_run);

  }
  else {                            ///// Create new file
    std::cout << "What is your run of pdf reference?" << '\n';
    std::cin >> run_number_pdf;
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
    std::cout << "Code start running" << '\n';
    for (int i = 0; i < n_run; i++) {
      Fit_Ref(run_number.at(i),run_number_pdf);
      cout << Form("Fit_Ref_%d", run_number.at(i)) << endl;
      minerror_calculator(Form("Fit_Ref_%d", run_number.at(i)), run_number.at(i));
    }
    std::cout << "Fit_Ref and minerror ok" << '\n';
    //
    file_merger(run_number,run_number_pdf);
    std::cout << "file_merger ok" << '\n';
    //
    TGrapher(Form("Fit_Ref_%d-%d", run_number.at(0), run_number.at(run_number.size()-1)), n_run+1);
  }
  std::cout << "Finish !!!" << '\n';
  return 0;
}

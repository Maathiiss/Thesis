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
#include <TPaveStats.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TMarker.h>
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

void Load_spectre_alpha(int run_number){
  TFile *file = new TFile(Form("entree/root/histo_ref_%d.root", run_number), "READ");
  gROOT->cd();
  charge_spectre = (TH2F*)file->Get("histo_pm_charge_50");
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

  modele_pdf.plotOn(frame);   //plot on function can adjust automatically the two graphs


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
  //if(run_number == 999){
  // gSystem->mkdir(Form("sortie/om_%d/run_%d", om_number, run_number));
  // can->SaveAs(Form("sortie/om_%d/run_%d/fit_run_%d_om_%d_gain_%f.png", om_number, run_number, run_number, om_number, gain));
  //}

  // gSystem->mkdir(Form("/home/granjon/Documents/these/Analyse/corr_OM_ref/sortie/graph interessant/Avec pic alpha/gif/om_%d", om_number));
  // can->SaveAs(Form("/home/granjon/Documents/these/Analyse/corr_OM_ref/sortie/graph interessant/Avec pic alpha/gif/om_%d/fit_run_%d_om_%d_gain_%f.png",om_number,run_number, om_number, gain));
  // } 


  delete miniChi;
  delete can;
  delete frame;
  // delete miniLog;
  return *rootab;

}

void Fit_Ref(int run_number, int run_number_pdf, int bins_start, int bins_stop) {
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
      if(bins_start==25 && bins_stop == 55){//pic only
	for (int i = 0; i < bins_start; i++) {  //mettre les 0 au debut                        
          modele->SetBinContent(i,0);
          modele->SetBinError(i, 0);
	  spectre_om->SetBinContent(i,0);
          spectre_om->SetBinError(i, 0);
        }
	for (int i = bins_stop; i < 1024; i++) {  //mettre les 0 au debut                        
          modele->SetBinContent(i,0);
          modele->SetBinError(i, 0);
	  spectre_om->SetBinContent(i,0);
          spectre_om->SetBinError(i, 0);
        }
      }
      else{//depend if we start before or after pic
	for (int i = bins_start; i < bins_stop; i++) {  //mettre les 0 au debut
	  modele->SetBinContent(i,0);
	  modele->SetBinError(i, 0);
	  spectre_om->SetBinContent(i,0);
	  spectre_om->SetBinError(i, 0);
	}
      }
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
      if(bins_start==25 && bins_stop == 55){//if pic only
        for (int i = 0; i < bins_start; i++) {  //mettre les 0 au debut                                 
          modele->SetBinContent(i,0);
          modele->SetBinError(i, 0);
          spectre_om->SetBinContent(i,0);
          spectre_om->SetBinError(i, 0);
        }
        for (int i = bins_stop; i < 1024; i++) {  //mettre les 0 au debut                               
          modele->SetBinContent(i,0);
          modele->SetBinError(i, 0);
          spectre_om->SetBinContent(i,0);
          spectre_om->SetBinError(i, 0);
        }
      }
      else{//before or after pic
	for (int i = bins_start; i < bins_stop; i++) {
	  modele->SetBinContent(i,0);
	  modele->SetBinError(i, 0);
	  spectre_om->SetBinContent(i,0);
	  spectre_om->SetBinError(i, 0);
	}
      }
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
      cout<<Chi2<<endl;
      modele->Reset();
      Result_tree.Fill();
      // }
      delete modele;
    }

    delete spectre_om;
    // }
      file->Close();
  }
  string fileadd = "";
  if (bins_stop == 55 && bins_start== 0){
    fileadd = "_bckg";
  }
  if(bins_start==25 && bins_stop == 55){
    fileadd= "_pic";    
  }
  

  TFile new_file(Form("sortie/root_file_final/Fit_Ref_%d%s.root", run_number, fileadd.c_str()), "RECREATE");
  cout<<Form("sortie/root_file_final/Fit_Ref_%d%s.root", run_number, fileadd.c_str())<<endl;
  // TFile new_file(Form("root/Complete_Fit/Fit_Ref_%d.root", run_number), "RECREATE");
  new_file.cd();
  Result_tree.Write();
  new_file.Close();
  return;

}



double langau_func (double *x, double *par)   //pointeur vers tableau de double
//compute convolution between gaussian and Landau and add decay exponential
{
   const double LANDAU_MPV = par[1];
   const double LANDAU_WIDTH = par[0];
   const double TOT_INT = par[2];
   const double GAUS_WIDTH = par[3];
   const double EXP_CONST = par[4];
   const double EXP_LAMBDA = par[5];
   //Fit parameters:                                                                                    
   //par[0]=Width (scale) parameter of Landau density                                                   
   //par[1]=Most Probable (MP, location) parameter of Landau density                                    
   //par[2]=Total area (integral -inf to inf, normalization constant)                                   
   //par[3]=Width (sigma) of convoluted Gaussian function                                               
   //                                                                                                   
   //In the Landau distribution (represented by the CERNLIB approximation),                             
   //the maximum is located at x=-0.22278298 with the location parameter=0.                             
   //This shift is corrected within this function, so that the actual                                   
   //maximum is identical to the MP parameter.                                                          
      // Numeric constants                                                                              
      double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)                                             
      double mpshift  = -0.22278298;       // Landau maximum location                                   

      // Control constants                                                                              
      double np = 100.0;      // number of convolution steps                                            
      double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas                            

      // Variables                                                                                      
      double xx;
      double mpc;
      double fland;
      double sum = 0.0;
      double xlow,xupp;
      double step;
      double i;


      // MP shift correction                                                                            
      mpc = LANDAU_MPV - mpshift * LANDAU_WIDTH;

      // Range of convolution integral                                                                  
      xlow = x[0] - sc * GAUS_WIDTH;
      xupp = x[0] + sc * GAUS_WIDTH;
      
      step = (xupp-xlow) / np;
      
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
	xx = xlow + (i-.5) * step;
	fland = TMath::Landau(mpc - xx,0,LANDAU_WIDTH) / LANDAU_WIDTH ;
	sum += fland * TMath::Gaus(x[0],xx,GAUS_WIDTH);
	
	xx = xupp - (i-.5) * step;
	fland = TMath::Landau(mpc - xx,0,LANDAU_WIDTH) / LANDAU_WIDTH ;
	sum += fland * TMath::Gaus(x[0],xx,GAUS_WIDTH);
      }      
      return (TOT_INT * step * sum * invsq2pi / GAUS_WIDTH) + EXP_CONST*TMath::Exp(-x[0]/EXP_LAMBDA) ;
}



std::vector<double> mean_error_alpha_pdf(5);
std::vector<double> mean_alpha_pdf(5);


void Fit_alpha_pdf(int run_number){
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");
 
  double mean, mean_error, time, gain_error, gain;
  int om_num;

  for(int om_number=712; om_number<717; om_number++){
    if(om_number == 714){
      mean = 0;
      mean_error=0;
      gain = 0;
      gain_error=0;
    }
    TFile *file1 = new TFile(Form("/home/granjon/corr_om_ref/entree/Modele/Modele_OM_%d_%d.root",om_number, run_number), "READ");
    TH2D* TH2modele_alpha = (TH2D*)file1->Get(Form("Modele_Ref_OM_%d", om_number));
    TH1D* spectre_om_pdf = TH2modele_alpha->ProjectionY(Form("modele"), 5000, 5000);
    spectre_om_pdf->SetTitle(Form("Modele_Ref_OM_%d",om_number));
    spectre_om_pdf->GetXaxis()->SetRangeUser(0, 50000);
    
    TCanvas* canvas = new TCanvas;    
    if(om_number == 712 ){
      TF1 *f_langaus = new TF1 ("f_langaus", langau_func, 5000, 11000, 6);      
      f_langaus->SetParameters(240, 8700, 6e9, 600, 5e6, 6300);      
      spectre_om_pdf->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om_pdf->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om_pdf->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om_pdf->Draw("hist");
      f_langaus->Draw("lsame");
      TF1 *f_exp = new TF1 ("f_exp", "[0]*TMath::Exp(-x/[1])", 0, 20000);
      f_exp->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      f_exp->SetParameters(f_langaus->GetParameter(4), f_langaus->GetParameter(5));
      f_exp->Draw("same");
      f_exp->SetLineColor(kRed);
      f_exp->SetLineStyle(kDashed);
      mean_alpha_pdf[om_number-712] = f_langaus->GetParameter(1);
      mean_error_alpha_pdf[om_number-712] = f_langaus->GetParError(1);


    }
    else{
      TF1 *f_langaus = new TF1 ("f_langaus", langau_func, 6000, 15000, 6);
      f_langaus->SetParameters(300, 9000, 6e9, 600, 5e6, 5300);
      spectre_om_pdf->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om_pdf->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om_pdf->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om_pdf->Draw("hist");
      f_langaus->Draw("lsame");
      TF1 *f_exp = new TF1 ("f_exp", "[0]*TMath::Exp(-x/[1])", 0, 20000);
      f_exp->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      f_exp->SetParameters(f_langaus->GetParameter(4), f_langaus->GetParameter(5));
      f_exp->Draw("same");
      f_exp->SetLineColor(kRed);
      f_exp->SetLineStyle(kDashed);
      mean_alpha_pdf[om_number-712] = f_langaus->GetParameter(1);
      mean_error_alpha_pdf[om_number-712] = f_langaus->GetParError(1);

    }
    canvas->SaveAs(Form("sortie/root_file_final/alpha_fit/gif_fit/fit_run_%d_om_%d.png", run_number,om_number));
  }
}



void Fit_alpha(int run_number, int run_number_pdf){
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03,"x");
  gStyle->SetLabelSize(0.03,"y");
  Load_spectre_alpha(run_number);  
  double mean, mean_error, time, gain_error_moins, gain_error_plus, Chi2, ndf;
  int om_num, om_number;
  double gain;
  TFile file(Form("sortie/root_file_final/alpha_fit/Am_fit_%d.root", run_number), "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("om_number", &om_num);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error_plus", &gain_error_plus);
  Result_tree.Branch("gain_error_moins", &gain_error_moins);
  Result_tree.Branch("time", &time);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("ndf", &ndf);
  gROOT->cd();

  TFile *file1 = new TFile(Form("entree/root/histo_ref_%d.root", run_number), "READ");
  time = 0;
  TParameter<double> *param = new TParameter<double>("start_time", time);
  param = (TParameter<double>*)(file1->Get("start_time"));
  time = param->GetVal();
  file1->Close();


  
  for(int indice=0; indice<5; indice++){
    om_number = indice+712;
    if(om_number == 714){
      om_num = om_number;
      mean = 0;
      mean_error=0;
      gain = 0;
      gain_error_plus=0;
      gain_error_moins = 0;
      Chi2=10000;
      ndf = 1;
    }
    TH1D* spectre_om = NULL;
    spectre_om = spectre_charge(om_number);
    TCanvas* canvas = new TCanvas;    
    if(om_number == 712 ){
      TF1 *f_langaus = new TF1 ("f_langaus", langau_func, 4000, 10000, 6);
      f_langaus->SetParameters(240, 7200, 7e6, 600, 7000, 5000);
      spectre_om->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om->Draw();
      f_langaus->Draw("lsame");
      TF1 *f_exp = new TF1 ("f_exp", "[0]*TMath::Exp(-x/[1])", 0, 20000);
      f_exp->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      f_exp->SetParameters(f_langaus->GetParameter(4), f_langaus->GetParameter(5));
      f_exp->Draw("same");
      f_exp->SetLineColor(kRed);
      f_exp->SetLineStyle(3);
      mean = f_langaus->GetParameter(1);
      mean_error = f_langaus->GetParError(1);
      gain = mean/mean_alpha_pdf[indice];
      Chi2 = f_langaus->GetChisquare();
      ndf = f_langaus->GetNDF();
      gain_error_plus = (gain * sqrt((mean_error/mean)*(mean_error/mean)+(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice])*(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice])))/2;
      gain_error_moins = (gain * sqrt((mean_error/mean)*(mean_error/mean)+(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice])*(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice])))/2;
      om_num = om_number;
      Result_tree.Fill(); 
    }
    else{
      TF1 *f_langaus = new TF1 ("f_langaus", langau_func, 6000, 12000, 6);
      f_langaus->SetParameters(240, 8700, 1.5e7, 600, 10000, 6300);
      spectre_om->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      spectre_om->Fit("f_langaus", "R");
      f_langaus->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));      
      spectre_om->Draw();
      f_langaus->Draw("lsame");
      TF1 *f_exp = new TF1 ("f_exp", "[0]*TMath::Exp(-x/[1])", 0, 20000);
      f_exp->SetRange(f_langaus->GetParameter(1)-6*f_langaus->GetParameter(3), f_langaus->GetParameter(1)+4*f_langaus->GetParameter(3));
      f_exp->SetParameters(f_langaus->GetParameter(4), f_langaus->GetParameter(5));
      f_exp->Draw("same");
      f_exp->SetLineColor(kRed);
      f_exp->SetLineStyle(kDashed);
      Chi2 = f_langaus->GetChisquare();
      ndf = f_langaus->GetNDF();
      mean = f_langaus->GetParameter(1);
      mean_error = f_langaus->GetParError(1);
      gain = mean/mean_alpha_pdf[indice];
      gain_error_plus = (gain * sqrt((mean_error/mean)*(mean_error/mean)+(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice])*(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice])))/2;
      gain_error_moins = (gain * sqrt((mean_error/mean)*(mean_error/mean)+(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice])*(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice])))/2;
      om_num=om_number;
      Result_tree.Fill();
    }
    canvas->SaveAs(Form("sortie/root_file_final/alpha_fit/gif_fit/fit_run_%d_om_%d.png", run_number,om_number));
  }
  file.cd();
  Result_tree.Write();
  file.Close();
  for(double i : mean_alpha_pdf){
    std::cout<<"I = "<<i<<std::endl;
  }
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
      TFile *file_graph = new TFile("graph.root", "RECREATE");
      Chi2gain->SetName("nom");
      Chi2gain->Write();
      file_graph->Close();

    
    Result_tree.Fill();
    delete Chi2gain;
    delete poly;
  }

  newfile.cd();
  Result_tree.Write();
  newfile.Close();

}


void file_merger(std::vector<int> run_number, int run_number_pdf, string previous_file_s = "", string bckg = "") {
  TFile file(Form("sortie/root_file_final/Fit_Ref_%d-%d%s.root", run_number.at(0), run_number.at(run_number.size()-1),bckg.c_str()), "RECREATE");
  //TFile file(Form("sortie/root_file_final/%s.root", bckg.c_str()), "RECREATE");
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
    TFile tree_file(Form("sortie/root_file_final/Fit_Ref_%d%s_mini.root", run_number.at(i), bckg.c_str()), "READ");
    std::cout<<"ici "<<Form("sortie/root_file_final/Fit_Ref_%d%s_mini.root", run_number.at(i), bckg.c_str())<<"\n";
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




void file_merger_alpha(std::vector<int> run_number, int run_number_pdf, string previous_file_s = "") {
  TFile file(Form("sortie/root_file_final/alpha_fit/Fit_Ref_%d-%d_alpha.root", run_number.at(0), run_number.at(run_number.size()-1)), "RECREATE");
  double gain, gain_error_plus, gain_error_moins, mean, mean_error, Chi2, ndf;
  int om_number, int_run;
  double time;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error_plus", &gain_error_plus);
  Result_tree.Branch("gain_error_moins", &gain_error_moins);
  Result_tree.Branch("run_number", &int_run);
  Result_tree.Branch("time", &time);
  Result_tree.Branch("Chi2_alpha", &Chi2);
  Result_tree.Branch("ndf_alpha", &ndf);

  for (int i = 712; i < 717; i++) {
    om_number = i;
    gain = 1;
    mean = mean_alpha_pdf[i-712];
    mean_error = mean_error_alpha_pdf[i-712];
    gain_error_plus = mean_error/(2*mean);
    gain_error_moins = mean_error/(2*mean);	
    int_run = run_number_pdf;
    time = time_pdf;
    Chi2 = Chi2;
    ndf = ndf;
    Result_tree.Fill();
  }

  for (size_t i = 0; i < run_number.size(); i++) {
    TFile tree_file(Form("sortie/root_file_final/alpha_fit/Am_fit_%d.root", run_number.at(i)), "READ");
    TTree* tree = (TTree*)tree_file.Get("Result_tree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("om_number",1);
    tree->SetBranchAddress("om_number", &om_number);
    tree->SetBranchStatus("gain",1);
    tree->SetBranchAddress("gain", &gain);
    tree->SetBranchStatus("gain_error_plus",1);
    tree->SetBranchAddress("gain_error_plus", &gain_error_plus);
    tree->SetBranchStatus("gain_error_moins",1);
    tree->SetBranchAddress("gain_error_moins", &gain_error_moins);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &int_run);
    tree->SetBranchStatus("time",1);
    tree->SetBranchAddress("time", &time);
    tree->SetBranchStatus("mean",1);
    tree->SetBranchAddress("mean", &mean);
    tree->SetBranchStatus("mean_error",1);
    tree->SetBranchAddress("mean_error", &mean_error);
    tree->SetBranchStatus("Chi2",1);
    tree->SetBranchAddress("Chi2", &Chi2);
    tree->SetBranchStatus("ndf",1);
    tree->SetBranchAddress("ndf", &ndf);
    
    
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




/*
void tree_creation(std::string file_name, int n_run){
  int om_number, run_number;
  double time,gain;
  double gain_error_moins, gain_error_plus, alpha_gain_error_plus, alpha_gain_error_moins, gain_alpha, bckg_gain, bckg_gain_error_plus, bckg_gain_error_moins;
 TFile file(Form("/home/granjon/corr_om_ref/sortie/root_file_final/compare_tot/Compare_%s.root", file_name.c_str()), "RECREATE");
 TTree Result_tree("Result_tree","");
 Result_tree.Branch("om_number", &om_number);
 Result_tree.Branch("run_number", &run_number);
 Result_tree.Branch("time", &time); 
 Result_tree.Branch("gain", &gain);
 Result_tree.Branch("gain_error_plus", &gain_error_plus);
 Result_tree.Branch("gain_error_moins", &gain_error_moins);
 Result_tree.Branch("gain_alpha", &gain_alpha);
 Result_tree.Branch("alpha_gain_error_plus", &alpha_gain_error_plus);
 Result_tree.Branch("alpha_gain_error_moins", &alpha_gain_error_moins);
 Result_tree.Branch("bckg_gain", &bckg_gain);
 Result_tree.Branch("bckg_gain_error_plus", &bckg_gain_error_plus);
 Result_tree.Branch("bckg_gain_error_moins", &bckg_gain_error_moins);

 
  TFile tree_file(Form("sortie/root_file_final/%s.root", file_name.c_str()), "READ");  
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

  TFile tree_file2(Form("sortie/root_file_final/%s_bckg.root", file_name.c_str()), "READ");
  double time2;
  int om_number2, run_number2;
  double gain2;
  double gain_error_moins2, gain_error_plus2;
  vector<double> gain_store2, gain_error2;
  TTree* tree2 = (TTree*)tree_file2.Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("gain",1);
  tree2->SetBranchAddress("gain", &gain2);
  tree2->SetBranchStatus("gain_error_moins",1);
  tree2->SetBranchAddress("gain_error_moins", &gain_error_moins2);
  tree2->SetBranchStatus("gain_error_plus",1);
  tree2->SetBranchAddress("gain_error_plus", &gain_error_plus2);

 TFile tree_file1(Form("sortie/root_file_final/alpha_fit/%s_alpha.root", file_name.c_str()), "READ");
 double time1;
 int om_number1, run_number1;
 double gain1;
 double gain_error_plus1, gain_error_moins1;
 vector<double> gain_store1, gain_error1;
 TTree* tree1 = (TTree*)tree_file1.Get("Result_tree");
 tree1->SetBranchStatus("*",0);
 tree1->SetBranchStatus("gain",1);
 tree1->SetBranchAddress("gain", &gain1);
 tree1->SetBranchStatus("gain_error_plus",1);
 tree1->SetBranchAddress("gain_error_plus", &gain_error_plus1);
 tree1->SetBranchStatus("gain_error_moins",1);
 tree1->SetBranchAddress("gain_error_moins", &gain_error_moins1);

 for(int i=0; i<6; i++){
   string name = Form("OM_ref_%d",i);
   TMultiGraph *name.c_str() = new TMultiGraph();

   for(int j=0;j<6;j++){

   }
 }
 
 for(int entry = 0; entry<tree->GetEntries(); entry++){
   tree->GetEntry(entry);
   tree1->GetEntry(entry);
   gain_alpha = gain1;
   alpha_gain_error_plus = gain_error_plus1;
   alpha_gain_error_moins = gain_error_moins1;
   tree2->GetEntry(entry);
   bckg_gain = gain2;
   bckg_gain_error_plus = gain_error_plus2;
   bckg_gain_error_moins = gain_error_moins2;
   Result_tree.Fill();
 }
 file.cd();
 Result_tree.Write();
 file.Close();
}
*/  




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
      Fit_Ref(run_number.at(i),run_number_pdf, 0,25);
      minerror_calculator(Form("Fit_Ref_%d", run_number.at(i)), run_number.at(i));
    }
    std::cout << "Fit_Ref ok and minerror ok" << '\n';

    file_merger(run_number, run_number_pdf, old_run);
    std::cout << "file_merger ok" << '\n';

    //TGrapher(Form("Fit_ref_716-%d", run_number.at(run_number.size()-1)), n_run);

  }
  else {                            ///// Create new file
    //n_run = 1;
    n_run = 37;
    //int run_number_before[n_run] = {1305,1315,1322,1329,1353,1361,1369,1376,1383,1391,1398,1407,1418};
    int run_number_before[n_run] = {1231,1232,1233,1234,1235,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1252,1253,1254,1256,1279,1283,1289,1294,1305,1315,1322,1329,1353,1361,1369,1376,1383,1391,1398,1407,1418};
    //int run_number_before[n_run] = {1329};
    std::vector<int> run_number(run_number_before, run_number_before + n_run);
    run_number_pdf = 1230;
    
    // std::cout << "What is your run of pdf reference?" << '\n';
    // std::cin >> run_number_pdf;
    // std::cout << "How many run do you want ?" << '\n';
    // std::cin >> n_run;
    // std::cout << "Write the runs you want" << '\n';
    // while (compteur < n_run && cin >> run) {
    //   run_number.push_back(run);
    //   compteur++;
    //   if (compteur < n_run) {
    // 	std::cout << "Write the runs you want" << '\n';
    //   }
    // }
    std::cout << "Code start running" << '\n';
    for (int i = 0; i < n_run; i++) {
      //total shift
      // Fit_Ref(run_number.at(i),run_number_pdf, 0,25);
      // cout << Form("Fit_Ref_%d", run_number.at(i)) << endl;      
      // minerror_calculator(Form("Fit_Ref_%d", run_number.at(i)), run_number.at(i));

      // //without alpha pic shift
      // Fit_Ref(run_number.at(i),run_number_pdf, 0,55);
      // cout << Form("Fit_Ref_%d_bckg", run_number.at(i)) << endl;
      // minerror_calculator(Form("Fit_Ref_%d_bckg", run_number.at(i)), run_number.at(i));

      //only alpha pic shift
      Fit_Ref(run_number.at(i),run_number_pdf, 25,55);
      cout << Form("Fit_Ref_%d_pic", run_number.at(i)) << endl;
      minerror_calculator(Form("Fit_Ref_%d_pic", run_number.at(i)), run_number.at(i));

      //alpha pic fit
      //Fit_alpha_pdf(run_number_pdf); //to have the first value of the pic alpha
      //Fit_alpha(run_number.at(i),run_number_pdf);
    }

    
    std::cout << "Fit_Ref and minerror ok" << '\n';
    //
    // file_merger(run_number,run_number_pdf,"","");
    // std::cout << "file_merger ok" << '\n';

    
    // string bckg = "_bckg";
    // file_merger(run_number,run_number_pdf,"",bckg);
    // std::cout << "file_merger bckg ok" << '\n';
    
    string pic = "_pic";
    file_merger(run_number,run_number_pdf,"",pic);
    std::cout << "file_merger pic ok" << '\n';
    
    // file_merger_alpha(run_number,run_number_pdf);
    // std::cout << "file_merger_alpha ok" << '\n';
    //tree_creation(Form("Fit_Ref_%d-%d", run_number.at(0), run_number.at(run_number.size()-1)), n_run+1);
  }
  std::cout << "Finish !!!" << '\n';
  return 0;
}

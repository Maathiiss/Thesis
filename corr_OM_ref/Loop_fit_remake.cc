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
  //if(run_number == 999){
  //gSystem->mkdir(Form("sortie/om_%d/run_%d", om_number, run_number));
  //can->SaveAs(Form("sortie/om_%d/run_%d/fit_run_%d_om_%d_gain_%f.png", om_number, run_number, run_number, om_number, gain));
  //}

  //gSystem->mkdir(Form("/home/granjon/Documents/these/Analyse/corr_OM_ref/sortie/graph interessant/Avec pic alpha/gif/om_%d", om_number));
    //can->SaveAs(Form("/home/granjon/Documents/these/Analyse/corr_OM_ref/sortie/graph interessant/Avec pic alpha/gif/om_%d/fit_run_%d_om_%d_gain_%f.png",om_number,run_number, om_number, gain));
  // } 


  delete miniChi;
  delete can;
  delete frame;
  // delete miniLog;
  return *rootab;

}

void Fit_Ref(int run_number, int run_number_pdf, int bins) {
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
      
      for (int i = 0; i < bins; i++) {  //mettre les 0 au debut
        modele->SetBinContent(i,0);
        modele->SetBinError(i, 0);
        spectre_om->SetBinContent(i,0);
        spectre_om->SetBinError(i, 0);
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
    for (int gain_count = gainmin - 400; gain_count < gainmin + 400; gain_count+=1) {
      gain = (gain_bin_min + gain_bin_width*(gain_count-1));
      // if (gain < 1.01 && gain > 0.99) {
      modele = TH2modele->ProjectionY("modele", gain_count, gain_count);
      // modele->Draw();
      // spectre_om->Draw();
      
      for (int i = 0; i < bins; i++) {
        modele->SetBinContent(i,0);
        modele->SetBinError(i, 0);
        spectre_om->SetBinContent(i,0);
        spectre_om->SetBinError(i, 0);
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
      modele->Reset();
      Result_tree.Fill();
      // }
      delete modele;
    }

    delete spectre_om;
    // }

  }
  string fileadd = "";
  if (bins == 55){
    fileadd = "_bckg";
  }

  TFile new_file(Form("sortie/root_file_final/Fit_Ref_%d%s.root", run_number, fileadd.c_str()), "RECREATE");
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
      TF1 *f_langaus = new TF1 ("f_langaus", langau_func, 5000, 12000, 6);      
      f_langaus->SetParameters(240, 8700, 1.5e9, 600, 1e6, 6300);
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
  double mean, mean_error, time, gain_error;
  int om_num, om_number;
  double gain;
  TFile file(Form("sortie/root_file_final/alpha_fit/Am_fit_%d.root", run_number), "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("run_number", &run_number);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("om_number", &om_num);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error", &gain_error);
  Result_tree.Branch("time", &time);
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
      gain_error=0;
    }
    TH1D* spectre_om = NULL;
    spectre_om = spectre_charge(om_number);
    TCanvas* canvas = new TCanvas;    
    if(om_number == 712 ){
      TF1 *f_langaus = new TF1 ("f_langaus", langau_func, 5000, 12000, 6);
      f_langaus->SetParameters(240, 8700, 8e6, 600, 5000, 6300);
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
      mean = f_langaus->GetParameter(1);
      mean_error = f_langaus->GetParError(1);
      gain = mean/mean_alpha_pdf[indice];
      gain_error = gain * sqrt((mean_error/mean)*(mean_error/mean)+(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice])*(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice]));
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
      mean = f_langaus->GetParameter(1);
      mean_error = f_langaus->GetParError(1);
      gain = mean/mean_alpha_pdf[indice];
      gain_error = gain * sqrt((mean_error/mean)*(mean_error/mean)+(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice])*(mean_error_alpha_pdf[indice]/mean_alpha_pdf[indice]));
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
      time = time_pdf/4;
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
  double gain, gain_error, mean, mean_error;
  int om_number, int_run;
  double time;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("mean", &mean);
  Result_tree.Branch("mean_error", &mean_error);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("gain_error", &gain_error);
  Result_tree.Branch("run_number", &int_run);
  Result_tree.Branch("time", &time);

  for (int i = 712; i < 717; i++) {
    om_number = i;
    gain = 1;
    mean = mean_alpha_pdf[i-712];
    mean_error = mean_error_alpha_pdf[i-712];
    gain_error = mean_error/mean;  
    int_run = run_number_pdf;
    time = time_pdf;
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
    tree->SetBranchStatus("gain_error",1);
    tree->SetBranchAddress("gain_error", &gain_error);
    tree->SetBranchStatus("run_number",1);
    tree->SetBranchAddress("run_number", &int_run);
    tree->SetBranchStatus("time",1);
    tree->SetBranchAddress("time", &time);
    tree->SetBranchStatus("mean",1);
    tree->SetBranchAddress("mean", &mean);
    tree->SetBranchStatus("mean_error",1);
    tree->SetBranchAddress("mean_error", &mean_error);

    
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




void TGrapher(std::string file_name, int n_run){
  
  /************************************************************************************************************************************************************************************************************************************************** WITH SHIFT METHOD TOTAL *********************************************************************************************************************************************************************************************************************************************************/
  
  TFile file(Form("sortie/root_file_final/TGraph_%s.root", file_name.c_str()), "RECREATE");
  TFile tree_file(Form("sortie/root_file_final/%s.root", file_name.c_str()), "READ");
  double time;
  int om_number, run_number;
  double gain;
  double gain_error_moins, gain_error_plus;
  vector<double> gain_store;
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
  TMultiGraph *mg[5];
  file.cd();
  for (int i = 0; i < 5; i++) {
    mg[i] = new TMultiGraph();
    for (int j = 0; j < n_run; j++){
      tree->GetEntry(i+j*5);
      yaxis[j] = gain;
      gain_store.push_back(gain);
      yaxis_error_moins[j] = gain_error_moins;
      yaxis_error_plus[j] = gain_error_plus;
      xaxis[j] = time;
      xaxis_error_plus[j] = 0.00001;
      xaxis_error_moins[j] = 0.00001;
    }
    gain_graph[i] = new TGraphAsymmErrors(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);
    gain_graph[i]->SetMarkerColor(2);
    gain_graph[i]->SetMarkerStyle(5);
    gain_graph[i]->SetMarkerSize(2);
    canvas->cd(i+1);
    mg[i]->Add(gain_graph[i]);
    mg[i]->SetName(Form("fit_OM_ref_%d", om_number));
    mg[i]->SetNameTitle(Form("fit_OM_ref_%d", om_number), Form("Relative gain evolution of the ref OM %d measured with background spectra", om_number));
    mg[i]->GetXaxis()->SetTimeDisplay(1);
    mg[i]->GetXaxis()->SetTitle("Time");
    mg[i]->GetYaxis()->SetTitle("Gain evolution");
    mg[i]->GetYaxis()->SetRangeUser(0.9, 1.1);
    mg[i]->GetYaxis()->SetTitleOffset(0.9);
    if(i==2){      
      mg[i]->Draw("AP");
      
      TLine* line_ref_d = new TLine(xaxis[0],1,xaxis[n_run-1],1);
      line_ref_d->SetLineWidth(1);
      line_ref_d->SetLineStyle(2);
      line_ref_d->SetLineColor(2);
      TLine* line_ref_d1 = new TLine(xaxis[0],1.01,xaxis[n_run-1],1.01);
      line_ref_d1->SetLineWidth(1);
      line_ref_d1->SetLineStyle(2);
      line_ref_d1->SetLineColor(kBlack);
      TLine* line_ref_d2 = new TLine(xaxis[0],0.99,xaxis[n_run-1],0.99);
      line_ref_d2->SetLineWidth(1);
      line_ref_d2->SetLineStyle(2);
      line_ref_d2->SetLineColor(kBlack);
      line_ref_d->Draw("same");
      line_ref_d2->Draw("same");
      line_ref_d1->Draw("same");
            
    }
  }

  canvas->Update();


  /************************************************************************************************************************************************************************************************************************************************** WITH SHIFT METHOD BACKGROUND ONLY  *************************************************************************************************************************************************************************************************************************************************/

TFile tree_file2(Form("sortie/root_file_final/%s_bckg.root", file_name.c_str()), "READ");
  double time2;
  int om_number2, run_number2;
  double gain2;
  double gain_error_moins2, gain_error_plus2;
  vector<double> gain_store2;
  TTree* tree2 = (TTree*)tree_file2.Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_number2);
  tree2->SetBranchStatus("gain",1);
  tree2->SetBranchAddress("gain", &gain2);
  tree2->SetBranchStatus("gain_error_moins",1);
  tree2->SetBranchAddress("gain_error_moins", &gain_error_moins2);
  tree2->SetBranchStatus("gain_error_plus",1);
  tree2->SetBranchAddress("gain_error_plus", &gain_error_plus2);
  tree2->SetBranchStatus("run_number",1);
  tree2->SetBranchAddress("run_number", &run_number2);
  tree2->SetBranchStatus("time",1);
  tree2->SetBranchAddress("time", &time2);
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < n_run; j++){
        tree2->GetEntry(i+j*5);
        yaxis[j] = gain2;
        gain_store2.push_back(gain);
        yaxis_error_moins[j] = gain_error_moins2;
        yaxis_error_plus[j] = gain_error_plus2;
	xaxis[j] = time2;
        xaxis_error_plus[j] = 0.00001;
        xaxis_error_moins[j] = 0.00001;
    }
        gain_graph[i] = new TGraphAsymmErrors(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);
        gain_graph[i]->SetMarkerColor(kMagenta-4);
        gain_graph[i]->SetMarkerStyle(5);
        gain_graph[i]->SetMarkerSize(2);
        canvas->cd(i+1);
        mg[i]->Add(gain_graph[i]);
        mg[i]->Draw("AP");
    }
    canvas->Update();
 




/**************************************************************************************************************************************************************************************************************************************************WITH FIT ALPHA PIC METHOD***********************************************************************************************************************************************************************************************************************************************************/

  TFile tree_file1(Form("sortie/root_file_final/alpha_fit/%s_alpha.root", file_name.c_str()), "READ");
 double time1;
 int om_number1, run_number1;
 double gain1;
 double gain_error1;
 vector<double> gain_store1;
 TTree* tree1 = (TTree*)tree_file1.Get("Result_tree");
 tree1->SetBranchStatus("*",0);
 tree1->SetBranchStatus("om_number",1);
 tree1->SetBranchAddress("om_number", &om_number1);
 tree1->SetBranchStatus("gain",1);
 tree1->SetBranchAddress("gain", &gain1);
 tree1->SetBranchStatus("gain_error",1);
 tree1->SetBranchAddress("gain_error", &gain_error1);
 tree1->SetBranchStatus("run_number",1);
 tree1->SetBranchAddress("run_number", &run_number1);
 tree1->SetBranchStatus("time",1);
 tree1->SetBranchAddress("time", &time1);
 double xaxis_min = *std::min_element(xaxis, xaxis+n_run);
 double xaxis_max = *std::max_element(xaxis, xaxis+n_run);
  TLine* line_ref = new TLine(xaxis_min,1,xaxis_max,1);
  line_ref->SetLineWidth(1);
  line_ref->SetLineStyle(2); 
  line_ref->SetLineColor(2);
  TLine* line_ref_1 = new TLine(xaxis_min,1.01,xaxis_max,1.01);
  line_ref_1->SetLineWidth(1);
  line_ref_1->SetLineStyle(2);
  line_ref_1->SetLineColor(kBlack);
  TLine* line_ref_2 = new TLine(xaxis_min,0.99,xaxis_max,0.99);
  line_ref_2->SetLineWidth(1);
  line_ref_2->SetLineStyle(2);
  line_ref_2->SetLineColor(kBlack);

 for (int i = 0; i < 5; i++) {
   for (int j = 0; j < n_run; j++){
     tree1->GetEntry(i+j*5);
     gain_store1.push_back(gain1);
     yaxis[j] = gain1;
     yaxis_error_moins[j] = gain_error1/2;
     yaxis_error_plus[j] = gain_error1/2;
     xaxis[j] = time1;
     xaxis_error_plus[j] = 0.00001;
     xaxis_error_moins[j] = 0.00001;
   }
   if(i!=2){
   gain_graph[i] = new TGraphAsymmErrors(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);
   gain_graph[i]->SetMarkerColor(kBlue);
   gain_graph[i]->SetMarkerStyle(5);
   gain_graph[i]->SetMarkerSize(2);   
   canvas->cd(i+1);
   mg[i]->Add(gain_graph[i]);
   mg[i]->Draw("AP");
   }
   canvas->cd(i+1);
   line_ref->Draw("same");
   line_ref_2->Draw("same");
   line_ref_1->Draw("same");
 }
 canvas->cd(6);
 TLegend *legend = new TLegend(0.1,0.4,0.9,0.7);
 TMarker *marker = new TMarker();
 marker->SetMarkerColor(2);
 marker->SetMarkerStyle(5);
 marker->SetMarkerSize(2);
 TLine *line = new TLine(0, 0, 0, 0);
 line->SetLineColor(2);
 line->SetLineWidth(2);
 line->Draw(); 
 legend->AddEntry(marker, "with pdf shift total", "lp");
 TMarker *marker2 = new TMarker();
 marker2->SetMarkerColor(kMagenta-4);
 marker2->SetMarkerStyle(5);
 marker2->SetMarkerSize(2);
 TLine *line2 = new TLine(0, 0, 0, 0);
 line2->SetLineColor(kMagenta-4);
 line2->SetLineWidth(2);
 line2->Draw();
 legend->AddEntry(marker2, "with pdf shift background only", "lp");
 TMarker *marker1 = new TMarker();
 marker1->SetMarkerColor(kBlue);
 marker1->SetMarkerStyle(5);
 marker1->SetMarkerSize(2);

 TLine *line1 = new TLine(0, 0, 0, 0);
 line1->SetLineColor(2);
 line1->SetLineWidth(2);
 line1->Draw();
 legend->AddEntry(marker1, "with alpha fit pic", "lp");
 legend->Draw();
 TLine *line_ref_leg = new TLine(0, 0, 0, 0);
 line_ref_leg->SetLineColor(2);
 line_ref_leg->SetLineWidth(1);
 line_ref_leg->SetLineStyle(2);
 legend->AddEntry(line_ref_leg, "gain = 1","l");
 TLine *line_1_leg = new TLine(0, 0, 0, 0);
 line_1_leg->SetLineColor(kBlack);
 line_1_leg->SetLineWidth(1);
 line_1_leg->SetLineStyle(2);
 legend->AddEntry(line_1_leg, "#Delta gain = 1 %","l");

 file.cd();
 canvas->Update();
 canvas->Write();
 file.Close();






 
 TCanvas *c1 = new TCanvas("histo_compare","",1200,800);
 c1->Divide(3,2);
 std::vector<TH1D*> histograms(5); 
 for (int i = 0; i < 5; ++i) {
   std::string name = "h" + std::to_string(i+1);
   std::string title = "Difference between 2 methods for OM 71" + std::to_string(i+2) + "; Difference in %; Number of events";
   histograms[i] = new TH1D(name.c_str(), title.c_str(), 30, -1, 1);
   histograms[i]->SetTitleSize(0.035);
   //gStyle->SetStatFontSize(4);
   gStyle->SetStatW(0.35);
 }
 int compteur2 =0;
 for(int i=0; i<n_run*5;i++){
   //std::cout<<"ICI "<<gain_store.at(i)<<" "<<gain_store1.at(i)<<std::endl;
   //std::cout<<100*(gain_store.at(i)-gain_store1.at(i))/gain_store1.at(i)<<std::endl;
   if (gain_store.at(i)==gain_store1.at(i)){//changement OM lorsque gain est = 1                      
     compteur2++;
   }
   if(100*(gain_store.at(i)-gain_store1.at(i))/gain_store1.at(i)!=0){
     histograms[compteur2-1]->Fill(100*(gain_store.at(i)-gain_store1.at(i))/gain_store1.at(i));
   }
   c1->cd(compteur2);
   histograms[compteur2-1]->Draw("hist");   
 }

 c1->SaveAs(Form("/home/granjon/corr_om_ref/sortie/root_file_final/alpha_fit/compare_avec_sans/comparaison_%s.pdf",file_name.c_str()));

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
      Fit_Ref(run_number.at(i),run_number_pdf, 25);
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
      //total shift
      Fit_Ref(run_number.at(i),run_number_pdf, 25);
      cout << Form("Fit_Ref_%d", run_number.at(i)) << endl;      
      minerror_calculator(Form("Fit_Ref_%d", run_number.at(i)), run_number.at(i));

      //without alpha pic shift
      Fit_Ref(run_number.at(i),run_number_pdf, 55);
      cout << Form("Fit_Ref_%d", run_number.at(i)) << endl;
      minerror_calculator(Form("Fit_Ref_%d_bckg", run_number.at(i)), run_number.at(i));

      //alpha pic
      Fit_alpha_pdf(run_number_pdf); //to have the first value of the pic alpha
      Fit_alpha(run_number.at(i),run_number_pdf);
    }

    
    std::cout << "Fit_Ref and minerror ok" << '\n';
    //
    file_merger(run_number,run_number_pdf,"","");
    std::cout << "file_merger ok" << '\n';
    
    string bckg = "_bckg";
    file_merger(run_number,run_number_pdf,"",bckg);
    std::cout << "file_merger bckg ok" << '\n';
    
    file_merger_alpha(run_number,run_number_pdf);
    std::cout << "file_merger_alpha ok" << '\n';
    //
    TGrapher(Form("Fit_Ref_%d-%d", run_number.at(0), run_number.at(run_number.size()-1)), n_run+1);
  }
  std::cout << "Finish !!!" << '\n';
  return 0;
}

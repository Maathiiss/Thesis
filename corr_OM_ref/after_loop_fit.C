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



void after_loop_fit(){
  
  /************************************************************************************************************************************************************************************************************************************************** WITH SHIFT METHOD TOTAL *********************************************************************************************************************************************************************************************************************************************************/
  string file_name = "Fit_Ref_999-1418";
  int n_run = 50;
  TFile file(Form("sortie/root_file_final/TGraph_change_%s.root", file_name.c_str()), "RECREATE");
  TFile tree_file(Form("sortie/root_file_final/%s.root", file_name.c_str()), "READ");
  double time;
  int om_number, run_number;
  double gain;
  double gain_error_moins, gain_error_plus;
  vector<double> gain_store, time_diff, gain_error;
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
  TMultiGraph *mg_only[5];
  file.cd();
  for (int i = 0; i < 5; i++) {
    mg[i] = new TMultiGraph();
    mg_only[i] = new TMultiGraph();
    for (int j = 0; j < n_run; j++){
      tree->GetEntry(i+j*5);
      yaxis[j] = gain;
      gain_store.push_back(gain);
      gain_error.push_back(gain_error_moins);
      time_diff.push_back(time);      
      yaxis_error_moins[j] = gain_error_moins;
      yaxis_error_plus[j] = gain_error_plus;
      if(run_number==986){
	xaxis[j] = time/4;
      }
      else{
	xaxis[j] = time;
      }      
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
    std::string om_ref[5] = {"MW1","MW2","GV","XW1","XW2"};
    mg[i]->SetNameTitle(Form("fit_OM_ref_%d", om_number), Form("Relative gain evolution ref OM %s", om_ref[i].c_str()));
    mg[i]->GetXaxis()->SetTimeDisplay(1);
    mg[i]->GetXaxis()->SetTitle("Time");
    mg[i]->GetYaxis()->SetTitle("Gain evolution");
    mg[i]->GetYaxis()->SetTitleSize(0.05);
    mg[i]->GetYaxis()->SetLabelSize(0.03);
    mg[i]->GetYaxis()->SetTitleOffset(1.1);
    mg[i]->GetYaxis()->SetRangeUser(0.9, 1.1);
    mg[i]->GetXaxis()->SetTitleSize(0.05);
    mg[i]->GetXaxis()->SetLabelSize(0.03);
    mg[i]->GetXaxis()->SetRangeUser(0.9, 1.1);
    mg[i]->GetXaxis()->SetTitleOffset(0.9);

    mg_only[i]->Add(gain_graph[i]);
    mg_only[i]->SetName(Form("fit_OM_ref_%d", om_number));
    mg_only[i]->SetNameTitle(Form("fit_OM_ref_%d", om_number), Form("Relative gain evolution of the ref OM %d measured with background spectra", om_number));
    mg_only[i]->GetXaxis()->SetTimeDisplay(1);
    mg_only[i]->GetXaxis()->SetTitle("Time");
    mg_only[i]->GetYaxis()->SetTitle("Gain evolution");
    mg_only[i]->GetYaxis()->SetRangeUser(0.9, 1.1);
    mg_only[i]->GetYaxis()->SetTitleOffset(1.1);
    
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
  
// TFile tree_file2(Form("sortie/root_file_final/%s_bckg.root", file_name.c_str()), "READ");
//   double time2;
//   int om_number2, run_number2;
//   double gain2;
//   double gain_error_moins2, gain_error_plus2;
//   vector<double> gain_store2, gain_error2;
//   TTree* tree2 = (TTree*)tree_file2.Get("Result_tree");
//   tree2->SetBranchStatus("*",0);
//   tree2->SetBranchStatus("om_number",1);
//   tree2->SetBranchAddress("om_number", &om_number2);
//   tree2->SetBranchStatus("gain",1);
//   tree2->SetBranchAddress("gain", &gain2);
//   tree2->SetBranchStatus("gain_error_moins",1);
//   tree2->SetBranchAddress("gain_error_moins", &gain_error_moins2);
//   tree2->SetBranchStatus("gain_error_plus",1);
//   tree2->SetBranchAddress("gain_error_plus", &gain_error_plus2);
//   tree2->SetBranchStatus("run_number",1);
//   tree2->SetBranchAddress("run_number", &run_number2);
//   tree2->SetBranchStatus("time",1);
//   tree2->SetBranchAddress("time", &time2);
//     for (int i = 0; i < 5; i++) {
//       for (int j = 0; j < n_run; j++){
//         tree2->GetEntry(i+j*5);
//         yaxis[j] = gain2;
//         gain_store2.push_back(gain);
// 	gain_error2.push_back(gain_error_moins2);
//         yaxis_error_moins[j] = gain_error_moins2;
//         yaxis_error_plus[j] = gain_error_plus2;
// 	xaxis[j] = time2;
//         xaxis_error_plus[j] = 0.00001;
//         xaxis_error_moins[j] = 0.00001;
//     }
//         gain_graph[i] = new TGraphAsymmErrors(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);
//         gain_graph[i]->SetMarkerColor(kMagenta-4);
//         gain_graph[i]->SetMarkerStyle(5);
//         gain_graph[i]->SetMarkerSize(2);
//         canvas->cd(i+1);
//         mg[i]->Add(gain_graph[i]);
//         mg[i]->Draw("AP");
//     }
//     canvas->Update();






      /************************************************************************************************************************************************************************************************************************************************** WITH SHIFT METHOD PIC ONLY  *************************************************************************************************************************************************************************************************************************************************/
  
// TFile tree_file3(Form("sortie/root_file_final/%s_pic.root", file_name.c_str()), "READ");
//   double time3;
//   int om_number3, run_number3;
//   double gain3;
//   double gain_error_moins3, gain_error_plus3;
//   vector<double> gain_store3, gain_error3;
//   TTree* tree3 = (TTree*)tree_file3.Get("Result_tree");
//   tree3->SetBranchStatus("*",0);
//   tree3->SetBranchStatus("om_number",1);
//   tree3->SetBranchAddress("om_number", &om_number3);
//   tree3->SetBranchStatus("gainmin",1);
//   tree3->SetBranchAddress("gainmin", &gain3);
//   tree3->SetBranchStatus("error_moins",1);
//   tree3->SetBranchAddress("error_moins", &gain_error_moins3);
//   tree3->SetBranchStatus("error_plus",1);
//   tree3->SetBranchAddress("error_plus", &gain_error_plus3);
//   tree3->SetBranchStatus("run_number",1);
//   tree3->SetBranchAddress("run_number", &run_number3);
//   tree3->SetBranchStatus("time",1);
//   tree3->SetBranchAddress("time", &time3);
//     for (int i = 0; i < 5; i++) {
//       for (int j = 0; j < n_run; j++){
//         tree3->GetEntry(i+j*5);
//         yaxis[j] = gain3;
//         gain_store3.push_back(gain);
// 	gain_error3.push_back(gain_error_moins3);
//         yaxis_error_moins[j] = gain_error_moins3;
//         yaxis_error_plus[j] = gain_error_plus3;
// 	xaxis[j] = time3;
//         xaxis_error_plus[j] = 0.00001;
//         xaxis_error_moins[j] = 0.00001;
//     }
//          if(i!=2){	   
//         gain_graph[i] = new TGraphAsymmErrors(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);
//         gain_graph[i]->SetMarkerColor(kGreen+1);
//         gain_graph[i]->SetMarkerStyle(5);
//         gain_graph[i]->SetMarkerSize(2);
//         canvas->cd(i+1);
//         mg[i]->Add(gain_graph[i]);
//         mg[i]->Draw("AP");
// 	 }
//     }
//     canvas->Update();





/**************************************************************************************************************************************************************************************************************************************************WITH FIT ALPHA PIC METHOD***********************************************************************************************************************************************************************************************************************************************************/

  TFile tree_file1(Form("sortie/root_file_final/%s_alpha.root", file_name.c_str()), "READ");
 double time1;
 int om_number1, run_number1;
 double gain1;
 double gain_error_plus1, gain_error_moins1;
 vector<double> gain_store1, gain_error1;
 TTree* tree1 = (TTree*)tree_file1.Get("Result_tree");
 tree1->SetBranchStatus("*",0);
 tree1->SetBranchStatus("om_number",1);
 tree1->SetBranchAddress("om_number", &om_number1);
 tree1->SetBranchStatus("gain",1);
 tree1->SetBranchAddress("gain", &gain1);
 tree1->SetBranchStatus("gain_error_plus",1);
 tree1->SetBranchAddress("gain_error_plus", &gain_error_plus1);
 tree1->SetBranchStatus("gain_error_moins",1);
 tree1->SetBranchAddress("gain_error_moins", &gain_error_moins1);
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
     gain_error1.push_back(gain_error_plus1);
     yaxis[j] = gain1;
     yaxis_error_moins[j] = gain_error_moins1;
     yaxis_error_plus[j] = gain_error_plus1;
     if(j==0){
       yaxis_error_moins[j] = 0;
       yaxis_error_plus[j] = 0;
     }
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
TLine *line3 = new TLine(0, 0, 0, 0);
 TMarker *marker3 = new TMarker();
 marker3->SetMarkerColor(kGreen+1);
 marker3->SetMarkerStyle(5);
 marker3->SetMarkerSize(2);
 line3->SetLineColor(kGreen+1);
 line3->SetLineWidth(2);
 line3->Draw();
 legend->AddEntry(marker3, "with alpha shift pic", "lp");
 legend->Draw();

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
 canvas->Draw();
 file.Close();
}


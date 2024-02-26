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

void compare_methode(){
  //def du graph

  int n_run;
  std::cout<<"How many runs"<<std::endl;
  std::cin>>n_run;
  double yaxis[n_run];
  double yaxis_error_moins[n_run];
  double yaxis_error_plus[n_run];
  double xaxis[n_run];
  double xaxis_error_moins[n_run];
  double xaxis_error_plus[n_run];
  TGraphAsymmErrors *gain_graph[5];
  TMultiGraph *mg[5];
 
  string file_name="";
  std::cout<<"What is the file number"<<std::endl;
  std::cin>>file_name;  
  auto compare = new TCanvas("Comparison","",1200,800);
  compare->Divide(3,2);

  std::cout<<Form("sortie/root_file_final/Fit_Ref_%s.root", file_name.c_str())<<std::endl;
  TFile tree_file(Form("sortie/root_file_final/Fit_Ref_%s.root", file_name.c_str()), "READ");
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
  int compteur = 0;
    for (int i = 0; i < 5; i++) {
      mg[i] = new TMultiGraph();
      for (int j = 0; j < n_run; j++){
	tree->GetEntry(i+j*5);
	compteur++;
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
	compare->cd(i+1);
	mg[i]->Add(gain_graph[i]);
	mg[i]->SetName(Form("fit_OM_ref_%d", om_number));
	mg[i]->SetNameTitle(Form("fit_OM_ref_%d", om_number), Form("Relative gain evolution of the ref OM %d measured with background spectra", om_number));
	mg[i]->GetXaxis()->SetTimeDisplay(1);
	mg[i]->GetXaxis()->SetTitle("Time");
	mg[i]->GetYaxis()->SetTitle("Gain evolution");
	mg[i]->GetYaxis()->SetRangeUser(0.9, 1.1);
	mg[i]->GetYaxis()->SetTitleOffset(0.9);
        mg[i]->Draw("AP");
    }






  std::cout<<Form("sortie/root_file_final/Fit_Ref_%s_alpha.root", file_name.c_str())<<std::endl;
  TFile tree_file1(Form("sortie/root_file_final/Fit_Ref_%s_alpha.root", file_name.c_str()), "READ");
  double time1;
  int om_number1, run_number1;
  double gain1;
  double gain_error_moins1, gain_error_plus1;
  vector<double> gain_store1;
  TTree* tree1 = (TTree*)tree_file1.Get("Result_tree");
  tree1->SetBranchStatus("*",0);
  tree1->SetBranchStatus("om_number",1);
  tree1->SetBranchAddress("om_number", &om_number1);
  tree1->SetBranchStatus("gain",1);
  tree1->SetBranchAddress("gain", &gain1);
  tree1->SetBranchStatus("gain_error_moins",1);
  tree1->SetBranchAddress("gain_error_moins", &gain_error_moins1);
  tree1->SetBranchStatus("gain_error_plus",1);
  tree1->SetBranchAddress("gain_error_plus", &gain_error_plus1);
  tree1->SetBranchStatus("run_number",1);
  tree1->SetBranchAddress("run_number", &run_number1);
  tree1->SetBranchStatus("time",1);
  tree1->SetBranchAddress("time", &time1);
  int compteur1 = 0;
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < n_run; j++){
	tree1->GetEntry(i+j*5);
	compteur1++;
	yaxis[j] = gain1;
	yaxis_error_moins[j] = gain_error_moins1;
	yaxis_error_plus[j] = gain_error_plus1;
	xaxis[j] = time1;
	gain_store1.push_back(gain1);
	xaxis_error_plus[j] = 0.00001;
	xaxis_error_moins[j] = 0.00001;
    }
	gain_graph[i] = new TGraphAsymmErrors(n_run, xaxis, yaxis, xaxis_error_moins, xaxis_error_plus, yaxis_error_moins, yaxis_error_plus);
	gain_graph[i]->SetMarkerColor(kBlue);
	gain_graph[i]->SetMarkerStyle(5);
	gain_graph[i]->SetMarkerSize(2);
	compare->cd(i+1);
	mg[i]->Add(gain_graph[i]);
	
	mg[i]->Draw("AP");
    }
    compare->Update();



    compare->cd(6);
    
    TLegend *legend = new TLegend(0.1,0.4,0.9,0.7); 
    TMarker *marker = new TMarker();
    marker->SetMarkerColor(2); 
    marker->SetMarkerStyle(5); 
    marker->SetMarkerSize(2);  
    marker->DrawMarker(0, 0);
    
   
    TLine *line = new TLine(0, 0, 0, 0); 
    line->SetLineColor(2); 
    line->SetLineWidth(2); 
    line->Draw();
    
    legend->AddEntry(marker, "with alpha pic", "lp"); 
    
    TMarker *marker1 = new TMarker();
    marker1->SetMarkerColor(kBlue); 
    marker1->SetMarkerStyle(5); 
    marker1->SetMarkerSize(2);  
    marker1->DrawMarker(0, 0); 
    TLine *line1 = new TLine(0, 0, 0, 0); 
    line1->SetLineColor(2); 
    line1->SetLineWidth(2); 
    line1->Draw();
    legend->AddEntry(marker1, "without alpha pic", "lp");
    legend->Draw();



    for (int bins=40 ; bins <=200; bins +=1){
    TCanvas *c1 = new TCanvas("Histo compare","",1200,800);
    TH1D *h1 = new TH1D("h1", "Difference between with and without alpha peak; Difference in %; Number of events", bins,-1,1);      
    for(int i=0; i<n_run*5;i++){
      std::cout<<100*(gain_store.at(i)-gain_store1.at(i))/gain_store1.at(i)<<std::endl;
      if(100*(gain_store.at(i)-gain_store1.at(i))/gain_store1.at(i)!=0){
      h1->Fill(100*(gain_store.at(i)-gain_store1.at(i))/gain_store1.at(i));
      }
    }
    c1->cd();
    h1->Draw();
    //c1->SaveAs(Form("Comparaison_bins_%d.png",bins));
    }
}




#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <TPad.h>
#include <TPaveStats.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <limits>
#include <TH1D.h>
#include "/home/granjon/Documents/stage_radon/Stage/Mathis/sndisplay/sndisplay.cc"


void each_OM_mean(){
  double diff_li_bi = 0.0, Chi2_bi = 0.0, li_relative = 0.0, bi_relative = 0.0, ndf_bi = 0.0, Mean = 0.0, RMS = 0.0, variation_li = 0.0, Mean_li = 0.0, variation_bi = 0.0,Mean_bi = 0.0, diff_li_start_stop = 0.0, diff_bi_start_stop = 0.0, min_somme_li_bi_save = 0.0, min_somme_li_bi_save_test = 0.0, li_amplitude = 0.0, li_amplitude_uncorr=0.0, mean_pic_bi=0.0, best_charge_diff = 0.0, charge_diff = 0.0, bi_amplitude = 0.0, li_amplitude1 = 0.0, best_li_amplitude1 = 0.0, min_charge_diff=0.0;
  int pic_li, nb_entries_bi, om_number_li, pic_li_after, om_number_li_after, bundle_li, run_number_li, om_number_bi, pic_associated, bundle_li1;

  sncalo = new sndisplay::calorimeter ("sndiplay_dis",1);
  //TFile* file = TFile::Open("/home/granjon/Bi/sortie/Comparison_Bi_Li/final_comparison_1088-1203.root", "READ");
  TFile* file = TFile::Open("/home/granjon/Bi/sortie/Comparison_Bi_Li/final_comparison_1278-1417.root", "READ");
  TTree* Result_tree = (TTree*)file->Get("Result_tree");
  Result_tree->SetBranchStatus("*",0);
  Result_tree->SetBranchStatus("diff_li_bi",1);
  Result_tree->SetBranchAddress("diff_li_bi",&diff_li_bi);
  Result_tree->SetBranchStatus("pic_li",1);
  Result_tree->SetBranchAddress("pic_li",&pic_li);
  Result_tree->SetBranchStatus("Chi2_bi",1);
  Result_tree->SetBranchAddress("Chi2_bi",&Chi2_bi);
  Result_tree->SetBranchStatus("nb_entries_bi",1);
  Result_tree->SetBranchAddress("nb_entries_bi",&nb_entries_bi);
  Result_tree->SetBranchStatus("ndf_bi",1);
  Result_tree->SetBranchAddress("ndf_bi",&ndf_bi);
  Result_tree->SetBranchStatus("nb_entries_bi",1);
  Result_tree->SetBranchAddress("nb_entries_bi",&nb_entries_bi);
  Result_tree->SetBranchStatus("li_relative",1);
  Result_tree->SetBranchAddress("li_relative",&li_relative);
  Result_tree->SetBranchStatus("bi_relative",1);
  Result_tree->SetBranchAddress("bi_relative",&bi_relative);
  Result_tree->SetBranchStatus("om_number_li",1);
  Result_tree->SetBranchAddress("om_number_li",&om_number_li);
  Result_tree->SetBranchStatus("om_number_bi",1);
  Result_tree->SetBranchAddress("om_number_bi",&om_number_bi);
  Result_tree->SetBranchStatus("run_number_li",1);
  Result_tree->SetBranchAddress("run_number_li",&run_number_li);
  Result_tree->SetBranchStatus("bundle_li",1);
  Result_tree->SetBranchAddress("bundle_li",&bundle_li);
  Result_tree->SetBranchStatus("diff_li_start_stop",1);
  Result_tree->SetBranchAddress("diff_li_start_stop",&diff_li_start_stop);
  Result_tree->SetBranchStatus("diff_bi_start_stop",1);
  Result_tree->SetBranchAddress("diff_bi_start_stop",&diff_bi_start_stop);
  Result_tree->SetBranchStatus("li_amplitude",1);
  Result_tree->SetBranchAddress("li_amplitude",&li_amplitude);
  Result_tree->SetBranchStatus("li_amplitude_uncorr",1);
  Result_tree->SetBranchAddress("li_amplitude_uncorr",&li_amplitude_uncorr);
  Result_tree->SetBranchStatus("mean_pic_bi",1);
  Result_tree->SetBranchAddress("mean_pic_bi",&mean_pic_bi);

  
  TFile *file_create = new TFile("/home/granjon/Bi/sortie/Comparison_Bi_Li/each_OM/each_OM_save_1278-1417.root", "RECREATE");
  TTree final_tree("Result_tree","");
  final_tree.Branch("Mean",&Mean);
  final_tree.Branch("RMS",&RMS);
  final_tree.Branch("pic_li",&pic_li_after);
  final_tree.Branch("om_number_li",&om_number_li_after);
  final_tree.Branch("bundle_li",&bundle_li1);
  final_tree.Branch("variation_li",&variation_li);
  final_tree.Branch("Mean_li",&Mean_li);
  final_tree.Branch("diff_li_start_stop",&diff_li_start_stop);
  final_tree.Branch("variation_bi",&variation_bi);
  final_tree.Branch("Mean_bi",&Mean_bi);
  final_tree.Branch("diff_bi_start_stop",&diff_bi_start_stop);
  final_tree.Branch("pic_associated",&pic_associated); //which intensity is associated with 
  final_tree.Branch("min_somme_li_bi_save",&min_somme_li_bi_save);
  final_tree.Branch("min_somme_li_bi_save_test",&min_somme_li_bi_save_test);
  final_tree.Branch("best_charge_diff",&best_charge_diff); // charge diff with better result in relative part (li-bi compare to 1 diff)
  final_tree.Branch("charge_diff",&charge_diff);
  final_tree.Branch("min_charge_diff",&min_charge_diff); //min charge diff in terms of absolute comparison between Li and Bi charge valeur (ex : 20000-40000)
  final_tree.Branch("li_amplitude",&li_amplitude1);
  final_tree.Branch("li_amplitude_uncorr",&li_amplitude_uncorr);
  final_tree.Branch("best_li_amplitude",&best_li_amplitude1); //best in terms of relative results
  final_tree.Branch("bi_amplitude",&bi_amplitude); 
  
  std::array<std::array<TH1D*,4>, 712> histograms; //712 OM et 4 pics
  std::array<std::array<TH1D*,4>, 712> histograms_variation; //712 OM et 4 pics
  std::array<TH1D*,712> histograms_bi;
  double diff_li_start_stop_vec[712][4] = {};
  double somme_diff_li_bi_vec[712][4] = {};
  double somme_charge_li_bi_vec[712][4] = {};
  double diff_bi_start_stop_vec[712] = {};
  double bundle_li_vec[712] = {};
  double li_amplitude_vec[712][4] = {};
  double li_amplitude_uncorr_vec[712][4] = {};
  double bi_amplitude_vec[712] = {};
  double best_pic[712] = {};
  double best_charge_pic[712] = {};
  
  for(int i=0;i<712;i++){
    histograms_bi[i] = new TH1D(Form("om_%d",i),Form("om_%d",i),200,-100,100);
    for (size_t j = 0; j <4; j++){
      histograms[i][j] = new TH1D(Form("om_%d_pic_%lu",i,j+1),Form("om_%d_pic_%lu",i,j+1),200,-100,100);
      histograms_variation[i][j] = new TH1D(Form("om_var_%d_pic_%lu",i,j+1),Form("om_var_%d_pic_%lu",i,j+1),200,-100,100);      
    }
  }
  
  for (int j = 0; j < Result_tree->GetEntries(); j++) {
    Result_tree->GetEntry(j);
    if(diff_li_start_stop!=0){
      diff_li_start_stop_vec[om_number_li][pic_li-1] = diff_li_start_stop;
    }
    if(diff_bi_start_stop!=0){
      diff_bi_start_stop_vec[om_number_bi] = diff_bi_start_stop;
    }
    if(li_relative!=0 && li_relative!=1){
      if(!((bundle_li ==4 && run_number_li>1300 && run_number_li<1365)|| (bundle_li ==9 && run_number_li>1300 && run_number_li<1365))){
	histograms_variation[om_number_li][pic_li-1]->Fill(li_relative);
	li_amplitude_vec[om_number_li][pic_li-1] += li_amplitude;
	li_amplitude_uncorr_vec[om_number_li][pic_li-1] += li_amplitude_uncorr;
	bundle_li_vec[om_number_li] = bundle_li;
	if(Chi2_bi!=2000 && nb_entries_bi>50 && Chi2_bi/ndf_bi<3 && bi_relative!=1 && bi_relative>0.1){
	  histograms_bi[om_number_bi]->Fill(bi_relative);
	  histograms[om_number_li][pic_li-1]->Fill(diff_li_bi);
	  //histograms[om_number_li][pic_li-1]->Fill(abs((li_amplitude-mean_pic_bi)/mean_pic_bi));
	  somme_diff_li_bi_vec[om_number_li][pic_li-1]+=abs(diff_li_bi)/*diff_li_bi*diff_li_bi/(100*100)*/;
  	  somme_charge_li_bi_vec[om_number_li][pic_li-1]+=abs((li_amplitude-mean_pic_bi)/mean_pic_bi);	
	  bi_amplitude_vec[om_number_bi] += mean_pic_bi;
	  // if(om_number_li==23){
	  //   cout<<"intensity "<<pic_li << " li amplitude "<<li_amplitude<<" bi amplitude "<<mean_pic_bi<<endl;
	  // }
 	}
      }
    }
  }
  
  for(int i=0; i<712; i++){
    bundle_li1 = bundle_li_vec[i];
    if(histograms_bi[i]->GetEntries()>0){
    bi_amplitude = bi_amplitude_vec[i]/histograms_bi[i]->GetEntries();
    }
    else{
      bi_amplitude = 0;
    }
    charge_diff = 0;
    best_charge_diff = 0;
    om_number_li_after = i;
    variation_bi = histograms_bi[i]->GetRMS();
    Mean_bi = histograms_bi[i]->GetMean();
    diff_bi_start_stop = diff_bi_start_stop_vec[i];
    pic_associated = -1;
    double min_somme_li_bi = 1000;
    double min_charge_diff_temp = 1000;
    int j_store=0;
    for (size_t j = 0; j <4; j++){
      if(min_somme_li_bi>somme_diff_li_bi_vec[i][j]/histograms[i][j]->GetEntries() && somme_diff_li_bi_vec[i][j]>0.0001 && somme_diff_li_bi_vec[i][j]<10e3 && histograms[i][j]->GetEntries() > 0){
	min_somme_li_bi = somme_diff_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
	pic_associated = j+1;
	best_charge_diff = somme_charge_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
	best_li_amplitude1 = li_amplitude_vec[i][j]/histograms_variation[i][j]->GetEntries();
      }
      if(min_charge_diff_temp>somme_charge_li_bi_vec[i][j]/histograms[i][j]->GetEntries()){
	min_charge_diff_temp = somme_charge_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
	j_store = j+1;
      }      
    }
    for (size_t j = 0; j <4; j++){
      li_amplitude1 = li_amplitude_vec[i][j]/histograms_variation[i][j]->GetEntries();
      li_amplitude_uncorr = li_amplitude_uncorr_vec[i][j]/histograms_variation[i][j]->GetEntries();
            // if(i==23){
	    //   cout<<" li amplitude "<<li_amplitude1<<" bi amplitude "<<bi_amplitude<<endl;
	    // }
      pic_li_after = j+1;
      diff_li_start_stop = diff_li_start_stop_vec[i][j];
      Mean = histograms[i][j]->GetMean();
      RMS = histograms[i][j]->GetRMS();	
      variation_li = histograms_variation[i][j]->GetRMS();
      Mean_li = histograms_variation[i][j]->GetMean();
      // if(i==23){
      // 	cout<<"intensity "<<j<<" somme_charge_li_bi "<<somme_charge_li_bi_vec[i][j]/histograms[i][j]->GetEntries()<<endl;
      // 	cout<<"entries "<<histograms[i][j]->GetEntries()<<endl;
      // }
      if(somme_charge_li_bi_vec[i][j]>0.0001 && histograms[i][j]->GetEntries() > 0){
	charge_diff = somme_charge_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
	min_charge_diff = min_charge_diff_temp;
	sncalo->setcontent(i,j_store);
	best_charge_pic[i] = j_store;
      }
      else{
	charge_diff=-10.0;
	min_charge_diff=0;
      }
      if(j+1 == pic_associated){
	min_somme_li_bi_save = somme_diff_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
	best_pic[i] = pic_associated;
	best_charge_diff = best_charge_diff;
      }		  
      else{
	min_somme_li_bi_save = 0;
      }
      min_somme_li_bi_save_test = somme_diff_li_bi_vec[i][j]/histograms[i][j]->GetEntries();
    
      final_tree.Fill();
      file_create->cd();      
      histograms[i][j]->Write();    
    }
    sncalo->settext(i, Form("%d", i));
  }
  file_create->cd();
  final_tree.Write();
  

  




  double best_diff_li_bi=0.0, best_diff_li_bi_charge = 0.0, RMS_li=0.0, RMS_bi=0.0, best_li_amplitude = 0.0,best_li_amplitude_charge =0.0, bi_amplitude1 = 0.0; 
  int best_pic_li;
  TFile *file_result = new TFile("/home/granjon/Bi/sortie/Comparison_Bi_Li/each_OM/final_file_best_%.root", "RECREATE");
  TTree result_tree("Result_tree","");
  result_tree.Branch("li_relative",&li_relative);
  result_tree.Branch("bundle_li",&bundle_li);
  result_tree.Branch("run_number_li",&run_number_li);
  result_tree.Branch("Chi2_bi",&Chi2_bi);
  result_tree.Branch("nb_entries_bi",&nb_entries_bi);
  result_tree.Branch("ndf_bi",&ndf_bi);
  result_tree.Branch("bi_relative",&bi_relative);
  result_tree.Branch("om_number_li",&om_number_li);
  result_tree.Branch("pic_li",&pic_li);
  result_tree.Branch("diff_li_bi",&diff_li_bi);
  result_tree.Branch("best_diff_li_bi",&best_diff_li_bi);
  result_tree.Branch("best_diff_li_bi_charge",&best_diff_li_bi_charge);
  result_tree.Branch("RMS_li",&RMS_li);
  result_tree.Branch("RMS_bi",&RMS_bi);
  result_tree.Branch("mean_pic_bi",&mean_pic_bi);
  result_tree.Branch("li_amplitude",&li_amplitude);
  result_tree.Branch("best_li_amplitude",&best_li_amplitude);
  result_tree.Branch("best_li_amplitude_charge",&best_li_amplitude_charge);
  result_tree.Branch("best_pic_li",&best_pic_li);


  for (int j = 0; j < Result_tree->GetEntries(); j++) {
    Result_tree->GetEntry(j);
    if(li_relative!=0 /*&& li_relative!=1*/){
      if(!((bundle_li ==4 && run_number_li>1300 && run_number_li<1365)|| (bundle_li ==9 && run_number_li>1300 && run_number_li<1365))){
        if(Chi2_bi!=2000 && nb_entries_bi>50 && Chi2_bi/ndf_bi<3 && /*bi_relative!=1 &&*/ bi_relative!=0){
  	  if(best_pic[om_number_li]==pic_li && diff_li_bi!=0){
	    RMS_li = histograms_variation[om_number_li][pic_li-1]->GetRMS();	    
	    RMS_bi = histograms_bi[om_number_li]->GetRMS();	    
	    best_diff_li_bi = diff_li_bi;
	    best_li_amplitude = li_amplitude;
	    bi_amplitude1 = mean_pic_bi;
	    best_pic_li=pic_li;
  	  }
	  else{
	    best_diff_li_bi = 0.0;
	    //best_li_amplitude = 0.0;	    	    
	  }
	  if(best_charge_pic[om_number_li]==pic_li){
	    best_diff_li_bi_charge = diff_li_bi;
	    best_li_amplitude_charge = li_amplitude;
	  }
	  else{
	    best_diff_li_bi_charge = 0.0;
	  }
        }
      }
    }
    result_tree.Fill();
  }
  for(int i=0;i<712;i++){
    delete histograms_bi[i];
    for(int j=0; j<4; j++){
      delete histograms[i][j];
      delete histograms_variation[i][j];
    }
  }
  
  file_result->cd();
  result_tree.Write();
  file_result->Close();
  file_create->Close();
  delete file_create;
  
  delete file;
  delete file_result;
  
  sncalo->draw1();
  sncalo->setrange(0,4);
  for(int i =0; i<712; i++){sncalo->settext(i, Form("%d", i));}
  sncalo->draw1();
  
  sncalo->canvas->SaveAs("/home/granjon/Bi/sortie/Comparison_Bi_Li/each_OM/calo_display_li_bi.png");
}

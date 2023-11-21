#ifndef VARIATION_H
#define VARIATION_H

#include <vector>
#include "TString.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"

class Variation
{
 public:
  Variation(TString identifier, TString prefix, TString order_n_str);
  ~Variation();
  TString ID;

  TH1D *h_vn_pp;
  TH1D *h_vn_pm;
  TH1D *h_vn_kp;
  TH1D *h_vn_km;
  TH1D *h_vn_pr;

  TH1D *h_vn_yCM_HADES;
  
  TH1D *h_vn_yCM_00to10_pr;
  TH1D *h_vn_yCM_10to40_pr;
  TH1D *h_vn_yCM_40to60_pr;
  TH1D *h_vn_yCM_00to10_pr_symm;
  TH1D *h_vn_yCM_10to40_pr_symm;
  TH1D *h_vn_yCM_40to60_pr_symm;

  TH1D *h_vn_yCM_00to05_pr_symm;
  TH1D *h_vn_yCM_05to10_pr_symm;
  TH1D *h_vn_yCM_10to15_pr_symm;
  TH1D *h_vn_yCM_15to20_pr_symm;
  TH1D *h_vn_yCM_20to25_pr_symm;
  TH1D *h_vn_yCM_25to30_pr_symm;
  TH1D *h_vn_yCM_30to35_pr_symm;
  TH1D *h_vn_yCM_35to40_pr_symm;
  TH1D *h_vn_yCM_40to45_pr_symm;
  TH1D *h_vn_yCM_45to50_pr_symm;
  TH1D *h_vn_yCM_50to55_pr_symm;
  TH1D *h_vn_yCM_55to60_pr_symm;

  TH1D *h_vn_pT_00to10_pr;
  TH1D *h_vn_pT_10to40_pr;
  TH1D *h_vn_pT_40to60_pr;

  TProfile *p_vn_pp;
  TProfile *p_vn_pm;
  TProfile *p_vn_kp;
  TProfile *p_vn_km;
  TProfile *p_vn_pr;

  TProfile *p_vn_yCM_HADES;
  
  TProfile *p_vn_yCM_00to10_pr;
  TProfile *p_vn_yCM_10to40_pr;
  TProfile *p_vn_yCM_40to60_pr;
  TProfile *p_vn_yCM_00to10_pr_symm;
  TProfile *p_vn_yCM_10to40_pr_symm;
  TProfile *p_vn_yCM_40to60_pr_symm;

  TProfile *p_vn_yCM_00to05_pr_symm;
  TProfile *p_vn_yCM_05to10_pr_symm;
  TProfile *p_vn_yCM_10to15_pr_symm;
  TProfile *p_vn_yCM_15to20_pr_symm;
  TProfile *p_vn_yCM_20to25_pr_symm;
  TProfile *p_vn_yCM_25to30_pr_symm;
  TProfile *p_vn_yCM_30to35_pr_symm;
  TProfile *p_vn_yCM_35to40_pr_symm;
  TProfile *p_vn_yCM_40to45_pr_symm;
  TProfile *p_vn_yCM_45to50_pr_symm;
  TProfile *p_vn_yCM_50to55_pr_symm;
  TProfile *p_vn_yCM_55to60_pr_symm;

  TProfile *p_vn_pT_00to10_pr;
  TProfile *p_vn_pT_10to40_pr;
  TProfile *p_vn_pT_40to60_pr;

 private:
  TString fileName;
  TFile *file;
  void initialize(TString order_n_str);
  void fixAttributes(TString order_n_str);
};

#endif

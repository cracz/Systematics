#ifndef DETACHEDCOMPARISON_H
#define DETACHEDCOMPARISON_H

#include <vector>
#include "TH1D.h"
#include "Variation.h"

class DetachedComparison
{
 public:
  DetachedComparison(TString prefix, Variation* data1, Variation* data2);

  struct ComparisonPoint
  {
    Double_t var1Value = -999.0;
    Double_t var1Error = -999.0;
    Double_t var2Value = -999.0;
    Double_t var2Error = -999.0;
    Double_t deltaAbsVal;
    Double_t deltaAbsValSquared;
  };
  
  std::vector<ComparisonPoint> v_vn_pp;
  std::vector<ComparisonPoint> v_vn_pm;
  std::vector<ComparisonPoint> v_vn_kp;
  std::vector<ComparisonPoint> v_vn_km;
  std::vector<ComparisonPoint> v_vn_pr;

  std::vector<ComparisonPoint> v_vn_yCM_HADES;
  
  std::vector<ComparisonPoint> v_vn_yCM_00to10_pr;
  std::vector<ComparisonPoint> v_vn_yCM_10to40_pr;
  std::vector<ComparisonPoint> v_vn_yCM_40to60_pr;
  std::vector<ComparisonPoint> v_vn_yCM_00to10_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_10to40_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_40to60_pr_symm;

  std::vector<ComparisonPoint> v_vn_yCM_00to05_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_05to10_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_10to15_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_15to20_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_20to25_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_25to30_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_30to35_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_35to40_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_40to45_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_45to50_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_50to55_pr_symm;
  std::vector<ComparisonPoint> v_vn_yCM_55to60_pr_symm;

  std::vector<ComparisonPoint> v_vn_pT_00to10_pr;
  std::vector<ComparisonPoint> v_vn_pT_10to40_pr;
  std::vector<ComparisonPoint> v_vn_pT_40to60_pr;

  std::vector<ComparisonPoint> v_vn_pT_bin6_pr_symm;
  std::vector<ComparisonPoint> v_vn_pT_bin7_pr_symm;
  std::vector<ComparisonPoint> v_vn_pT_bin8_pr_symm;
  std::vector<ComparisonPoint> v_vn_pT_bin9_pr_symm;
  std::vector<ComparisonPoint> v_vn_pT_bin10_pr_symm;
  std::vector<ComparisonPoint> v_vn_pT_bin11_pr_symm;
  std::vector<ComparisonPoint> v_vn_pT_bin12_pr_symm;

  std::vector<std::vector<ComparisonPoint>> v2_vn_yCM_cent_pr_symmetry;

  TString ID;

 private:
  void initialize();
  void mergePoints(TH1D* histo1, TH1D* histo2, std::vector<ComparisonPoint>& vectorOfPoints);
  void combine(Variation* data1, Variation* data2);
  void addRawValuesToFile(TFile* file, TString histogramName, std::vector<ComparisonPoint> vectorOfPoints);
  void saveDetails();
};

#endif

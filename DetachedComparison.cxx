#include <iostream>
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "PlotUtils.h"
#include "DetachedComparison.h"

DetachedComparison::DetachedComparison(TString prefix, Variation* data1, Variation* data2)
{
  ID = prefix;
  //initialize();
  combine(data1, data2);
}

// Calculate necessary info for one histogram type with normal data and one variation
void DetachedComparison::mergePoints(TH1D* var1Histo, TH1D* var2Histo, std::vector<ComparisonPoint>& vectorOfPoints)
{
  for (int i = 1; i <= var1Histo->GetNbinsX(); i++)
  {
    ComparisonPoint point;
    point.var1Value = var1Histo->GetBinContent(i);
    point.var1Error = var1Histo->GetBinError(i);
    point.var2Value = var2Histo->GetBinContent(i);
    point.var2Error = var2Histo->GetBinError(i);

    if((point.var1Value == 0.0 && point.var1Error == 0.0) || (point.var2Value == 0.0 && point.var2Error == 0.0))
      {
	point.deltaAbsVal = 0.0;
	point.deltaAbsValSquared = 0.0;
      }
    else
      {
	point.deltaAbsVal = TMath::Abs(point.var2Value - point.var1Value);
	point.deltaAbsValSquared = TMath::Power(TMath::Abs(point.var2Value - point.var1Value), 2.0);
      }

    vectorOfPoints.push_back(point);
  }
}


// Combine one variation with the normal data to get the attributes for systematics
void DetachedComparison::combine(Variation* var1Data, Variation* var2Data)
{
  mergePoints(var1Data->h_vn_pp, var2Data->h_vn_pp, v_vn_pp);
  mergePoints(var1Data->h_vn_pm, var2Data->h_vn_pm, v_vn_pm);
  mergePoints(var1Data->h_vn_kp, var2Data->h_vn_kp, v_vn_kp);
  mergePoints(var1Data->h_vn_km, var2Data->h_vn_km, v_vn_km);
  mergePoints(var1Data->h_vn_pr, var2Data->h_vn_pr, v_vn_pr);
  
  mergePoints(var1Data->h_vn_yCM_00to10_pr, var2Data->h_vn_yCM_00to10_pr, v_vn_yCM_00to10_pr);
  mergePoints(var1Data->h_vn_yCM_10to40_pr, var2Data->h_vn_yCM_10to40_pr, v_vn_yCM_10to40_pr);
  mergePoints(var1Data->h_vn_yCM_40to60_pr, var2Data->h_vn_yCM_40to60_pr, v_vn_yCM_40to60_pr);

  mergePoints(var1Data->h_vn_yCM_HADES, var2Data->h_vn_yCM_HADES, v_vn_yCM_HADES);

  mergePoints(var1Data->h_vn_yCM_00to10_pr_symm, var2Data->h_vn_yCM_00to10_pr_symm, v_vn_yCM_00to10_pr_symm);
  mergePoints(var1Data->h_vn_yCM_10to40_pr_symm, var2Data->h_vn_yCM_10to40_pr_symm, v_vn_yCM_10to40_pr_symm);
  mergePoints(var1Data->h_vn_yCM_40to60_pr_symm, var2Data->h_vn_yCM_40to60_pr_symm, v_vn_yCM_40to60_pr_symm);

  mergePoints(var1Data->h_vn_yCM_00to05_pr_symm, var2Data->h_vn_yCM_00to05_pr_symm, v_vn_yCM_00to05_pr_symm);
  mergePoints(var1Data->h_vn_yCM_05to10_pr_symm, var2Data->h_vn_yCM_05to10_pr_symm, v_vn_yCM_05to10_pr_symm);
  mergePoints(var1Data->h_vn_yCM_10to15_pr_symm, var2Data->h_vn_yCM_10to15_pr_symm, v_vn_yCM_10to15_pr_symm);
  mergePoints(var1Data->h_vn_yCM_15to20_pr_symm, var2Data->h_vn_yCM_15to20_pr_symm, v_vn_yCM_15to20_pr_symm);
  mergePoints(var1Data->h_vn_yCM_20to25_pr_symm, var2Data->h_vn_yCM_20to25_pr_symm, v_vn_yCM_20to25_pr_symm);
  mergePoints(var1Data->h_vn_yCM_25to30_pr_symm, var2Data->h_vn_yCM_25to30_pr_symm, v_vn_yCM_25to30_pr_symm);
  mergePoints(var1Data->h_vn_yCM_30to35_pr_symm, var2Data->h_vn_yCM_30to35_pr_symm, v_vn_yCM_30to35_pr_symm);
  mergePoints(var1Data->h_vn_yCM_35to40_pr_symm, var2Data->h_vn_yCM_35to40_pr_symm, v_vn_yCM_35to40_pr_symm);
  mergePoints(var1Data->h_vn_yCM_40to45_pr_symm, var2Data->h_vn_yCM_40to45_pr_symm, v_vn_yCM_40to45_pr_symm);
  mergePoints(var1Data->h_vn_yCM_45to50_pr_symm, var2Data->h_vn_yCM_45to50_pr_symm, v_vn_yCM_45to50_pr_symm);
  mergePoints(var1Data->h_vn_yCM_50to55_pr_symm, var2Data->h_vn_yCM_50to55_pr_symm, v_vn_yCM_50to55_pr_symm);
  mergePoints(var1Data->h_vn_yCM_55to60_pr_symm, var2Data->h_vn_yCM_55to60_pr_symm, v_vn_yCM_55to60_pr_symm);

  mergePoints(var1Data->h_vn_pT_00to10_pr, var2Data->h_vn_pT_00to10_pr, v_vn_pT_00to10_pr);
  mergePoints(var1Data->h_vn_pT_10to40_pr, var2Data->h_vn_pT_10to40_pr, v_vn_pT_10to40_pr);
  mergePoints(var1Data->h_vn_pT_40to60_pr, var2Data->h_vn_pT_40to60_pr, v_vn_pT_40to60_pr);

  mergePoints(var1Data->h_vn_pT_bin6_10to40_pr_symm,  var2Data->h_vn_pT_bin6_10to40_pr_symm,  v_vn_pT_bin6_pr_symm);
  mergePoints(var1Data->h_vn_pT_bin7_10to40_pr_symm,  var2Data->h_vn_pT_bin7_10to40_pr_symm,  v_vn_pT_bin7_pr_symm);
  mergePoints(var1Data->h_vn_pT_bin8_10to40_pr_symm,  var2Data->h_vn_pT_bin8_10to40_pr_symm,  v_vn_pT_bin8_pr_symm);
  mergePoints(var1Data->h_vn_pT_bin9_10to40_pr_symm,  var2Data->h_vn_pT_bin9_10to40_pr_symm,  v_vn_pT_bin9_pr_symm);
  mergePoints(var1Data->h_vn_pT_bin10_10to40_pr_symm, var2Data->h_vn_pT_bin10_10to40_pr_symm, v_vn_pT_bin10_pr_symm);
  mergePoints(var1Data->h_vn_pT_bin11_10to40_pr_symm, var2Data->h_vn_pT_bin11_10to40_pr_symm, v_vn_pT_bin11_pr_symm);
  mergePoints(var1Data->h_vn_pT_bin12_10to40_pr_symm, var2Data->h_vn_pT_bin12_10to40_pr_symm, v_vn_pT_bin12_pr_symm);
}

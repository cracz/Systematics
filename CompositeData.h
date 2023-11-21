#ifndef COMPOSITEDATA_H
#define COMPOSITEDATA_H

#include <vector>
#include "TH1D.h"
#include "Variation.h"

class CompositeData
{
 public:
  CompositeData(TString prefix, Variation* normalData, Variation* var1Data);
  CompositeData(TString prefix, Variation* normalData, Variation* var1Data, Variation* var2Data);
  CompositeData(TString prefix, Variation* normalData, Variation* var1Data, Variation* var2Data, Variation* var3Data);
  CompositeData(TString prefix, Variation* normalData, Variation* var1Data, Variation* var2Data, Variation* var3Data, Variation* var4Data);
  ~CompositeData();

  struct DataPoint
  {
    Double_t normalValue = -999.0;
    Double_t normalError = -999.0;
    Double_t var1Value = -999.0;
    Double_t var1Error = -999.0;
    Double_t var2Value = -999.0;
    Double_t var2Error = -999.0;
    Double_t var3Value = -999.0;
    Double_t var3Error = -999.0;
    Double_t var4Value = -999.0;
    Double_t var4Error = -999.0;
    Double_t maxValue;
    Double_t maxError;
    Double_t minValue;
    Double_t minError;
    Double_t delta;
    Double_t deltaError;
    Double_t deltaByDeltaError;
    Double_t mean;
    Double_t stdDev;
    Double_t variance;
    Double_t stdDevPercentage;

    void getMax()
    {
      // Get max out of normal value and first 2 variations, 
      //which should always be present when using this function.
      maxValue = normalValue;
      maxError = normalError;
      
      if (var1Value > maxValue)
      {
        maxValue = var1Value;
        maxError = var1Error;
      }
      if (var2Value > maxValue)
      {
        maxValue = var2Value;
        maxError = var2Error;
      }
      ////

      // Add check for 3rd variation if present
      if (var3Value != -999.0 && var3Error != -999.0)
      {
        if (var3Value > maxValue)
        {
          maxValue = var3Value;
          maxError = var3Error;
        }

        // Add check for 4th variation if present
        if (var4Value != -999.0 && var4Error != -999.0)
        {
          if (var4Value > maxValue)
          {
            maxValue = var4Value;
            maxError = var4Error;
          }
        }
      }
    }

    void getMin()
    {
      // Get min out of normal value and first 2 variations, 
      //which should always be present when using this function.
      minValue = normalValue;
      minError = normalError;
      
      if (var1Value < minValue)
      {
        minValue = var1Value;
        minError = var1Error;
      }
      if (var2Value < minValue)
      {
        minValue = var2Value;
        minError = var2Error;
      }
      ////

      // Add check for 3rd variation if present
      if (var3Value != -999.0 && var3Error != -999.0)
      {
        if (var3Value < minValue)
        {
          minValue = var3Value;
          minError = var3Error;
        }

        // Add check for 4th variation if present
        if (var4Value != -999.0 && var4Error != -999.0)
        {
          if (var4Value < minValue)
          {
            minValue = var4Value;
            minError = var4Error;
          }
        }
      }
    }
  }; // End struct DataPoint

  std::vector<DataPoint> v_vn_pp;
  std::vector<DataPoint> v_vn_pm;
  std::vector<DataPoint> v_vn_kp;
  std::vector<DataPoint> v_vn_km;
  std::vector<DataPoint> v_vn_pr;

  std::vector<DataPoint> v_vn_yCM_HADES;
  
  std::vector<DataPoint> v_vn_yCM_00to10_pr;
  std::vector<DataPoint> v_vn_yCM_10to40_pr;
  std::vector<DataPoint> v_vn_yCM_40to60_pr;
  std::vector<DataPoint> v_vn_yCM_00to10_pr_symm;
  std::vector<DataPoint> v_vn_yCM_10to40_pr_symm;
  std::vector<DataPoint> v_vn_yCM_40to60_pr_symm;

  std::vector<DataPoint> v_vn_yCM_00to05_pr_symm;
  std::vector<DataPoint> v_vn_yCM_05to10_pr_symm;
  std::vector<DataPoint> v_vn_yCM_10to15_pr_symm;
  std::vector<DataPoint> v_vn_yCM_15to20_pr_symm;
  std::vector<DataPoint> v_vn_yCM_20to25_pr_symm;
  std::vector<DataPoint> v_vn_yCM_25to30_pr_symm;
  std::vector<DataPoint> v_vn_yCM_30to35_pr_symm;
  std::vector<DataPoint> v_vn_yCM_35to40_pr_symm;
  std::vector<DataPoint> v_vn_yCM_40to45_pr_symm;
  std::vector<DataPoint> v_vn_yCM_45to50_pr_symm;
  std::vector<DataPoint> v_vn_yCM_50to55_pr_symm;
  std::vector<DataPoint> v_vn_yCM_55to60_pr_symm;

  std::vector<DataPoint> v_vn_pT_00to10_pr;
  std::vector<DataPoint> v_vn_pT_10to40_pr;
  std::vector<DataPoint> v_vn_pT_40to60_pr;

  std::vector<std::vector<DataPoint>> v2_vn_yCM_cent_pr_symmetry;

  // Histograms that hold the values for the Barlow check related to this set of variations
  TH1D* barlow_vn_pp;
  TH1D* barlow_vn_pm;
  TH1D* barlow_vn_kp;
  TH1D* barlow_vn_km;
  TH1D* barlow_vn_pr;

  TH1D* barlow_vn_yCM_HADES;
  
  TH1D* barlow_vn_yCM_00to10_pr;
  TH1D* barlow_vn_yCM_10to40_pr;
  TH1D* barlow_vn_yCM_40to60_pr;

  TH1D* barlow_vn_yCM_00to10_pr_symm;
  TH1D* barlow_vn_yCM_10to40_pr_symm;
  TH1D* barlow_vn_yCM_40to60_pr_symm;

  TH1D* barlow_vn_yCM_00to05_pr_symm;
  TH1D* barlow_vn_yCM_05to10_pr_symm;
  TH1D* barlow_vn_yCM_10to15_pr_symm;
  TH1D* barlow_vn_yCM_15to20_pr_symm;
  TH1D* barlow_vn_yCM_20to25_pr_symm;
  TH1D* barlow_vn_yCM_25to30_pr_symm;
  TH1D* barlow_vn_yCM_30to35_pr_symm;
  TH1D* barlow_vn_yCM_35to40_pr_symm;
  TH1D* barlow_vn_yCM_40to45_pr_symm;
  TH1D* barlow_vn_yCM_45to50_pr_symm;
  TH1D* barlow_vn_yCM_50to55_pr_symm;
  TH1D* barlow_vn_yCM_55to60_pr_symm;

  TH1D* barlow_vn_pT_00to10_pr;
  TH1D* barlow_vn_pT_10to40_pr;
  TH1D* barlow_vn_pT_40to60_pr;

  TString ID;
  bool onlyOneVariation;
  int nVariations;
  

 private:
  void initialize();
  void mergePoints(TH1D* normalHisto, TH1D* var1Histo, std::vector<DataPoint>& vectorOfPoints);
  void mergePoints(TH1D* normalHisto, TH1D* var1Histo, TH1D* var2Histo, std::vector<DataPoint>& vectorOfPoints);
  void mergePoints(TH1D* normalHisto, TH1D* var1Histo, TH1D* var2Histo, TH1D* var3Histo, std::vector<DataPoint>& vectorOfPoints);
  void mergePoints(TH1D* normalHisto, TH1D* var1Histo, TH1D* var2Histo, TH1D* var3Histo, TH1D* var4Histo, std::vector<DataPoint>& vectorOfPoints);
  void addRawValuesToFile(TFile* file, TString histogramName, std::vector<DataPoint> vectorOfPoints);
  void addBarlowValuesToFile(TFile* file, TH1D* barlowHistogram, std::vector<DataPoint> vectorOfPoints);
  void saveDetails(Variation* normalData);
  void combine(Variation* normalData, Variation* var1Data);
  void combine(Variation* normalData, Variation* var1Data, Variation* var2Data);
  void combine(Variation* normalData, Variation* var1Data, Variation* var2Data, Variation* var3Data);
  void combine(Variation* normalData, Variation* var1Data, Variation* var2Data, Variation* var3Data, Variation* var4Data);
};

#endif

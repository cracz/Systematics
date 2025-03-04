#include <iostream>
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "PlotUtils.h"
#include "CompositeData.h"


// Constructor for one variation.
CompositeData::CompositeData(TString prefix, Variation* normalData, Variation* var1Data)
{
  onlyOneVariation = true;
  nVariations = 1;
  ID = prefix;
  initialize();
  combine(normalData, var1Data);
  saveDetails(normalData);
}


// Constructor for two related variations.
CompositeData::CompositeData(TString prefix, Variation* normalData, Variation* var1Data, Variation* var2Data)
{
  onlyOneVariation = false;
  nVariations = 2;
  ID = prefix;
  initialize();
  combine(normalData, var1Data, var2Data);
  saveDetails(normalData);
}

// Constructor for three related variations.
CompositeData::CompositeData(TString prefix, Variation* normalData, Variation* var1Data, Variation* var2Data, Variation* var3Data)
{
  onlyOneVariation = false;
  nVariations = 3;
  ID = prefix;
  initialize();
  combine(normalData, var1Data, var2Data, var3Data);
  saveDetails(normalData);
}

// Constructor for four related variations.
CompositeData::CompositeData(TString prefix, Variation* normalData, Variation* var1Data, Variation* var2Data, Variation* var3Data, Variation* var4Data)
{
  onlyOneVariation = false;
  nVariations = 4;
  ID = prefix;
  initialize();
  combine(normalData, var1Data, var2Data, var3Data, var4Data);
  saveDetails(normalData);
}


CompositeData::~CompositeData()
{
  delete barlow_vn_pp;
  delete barlow_vn_pm;
  delete barlow_vn_kp;
  delete barlow_vn_km;
  delete barlow_vn_pr;

  delete barlow_vn_yCM_HADES;
  
  delete barlow_vn_yCM_00to10_pr;
  delete barlow_vn_yCM_10to40_pr;
  delete barlow_vn_yCM_40to60_pr;
  delete barlow_vn_yCM_00to10_pr_symm;
  delete barlow_vn_yCM_10to40_pr_symm;
  delete barlow_vn_yCM_40to60_pr_symm;

  delete barlow_vn_yCM_00to05_pr_symm;
  delete barlow_vn_yCM_05to10_pr_symm;
  delete barlow_vn_yCM_10to15_pr_symm;
  delete barlow_vn_yCM_15to20_pr_symm;
  delete barlow_vn_yCM_20to25_pr_symm;
  delete barlow_vn_yCM_25to30_pr_symm;
  delete barlow_vn_yCM_30to35_pr_symm;
  delete barlow_vn_yCM_35to40_pr_symm;
  delete barlow_vn_yCM_40to45_pr_symm;
  delete barlow_vn_yCM_45to50_pr_symm;
  delete barlow_vn_yCM_50to55_pr_symm;
  delete barlow_vn_yCM_55to60_pr_symm;

  delete barlow_vn_pT_00to10_pr;
  delete barlow_vn_pT_10to40_pr;
  delete barlow_vn_pT_40to60_pr;

  delete barlow_vn_pT_bin6_pr_symm;
  delete barlow_vn_pT_bin7_pr_symm;
  delete barlow_vn_pT_bin8_pr_symm;
  delete barlow_vn_pT_bin9_pr_symm;
  delete barlow_vn_pT_bin10_pr_symm;
  delete barlow_vn_pT_bin11_pr_symm;
  delete barlow_vn_pT_bin12_pr_symm;
}


void CompositeData::initialize()
{
  // Centrality
  barlow_vn_pp = new TH1D("barlow_vn_pp_"+ID, "pp vs cent;Centrality;#Delta/#sigma_{#Delta}", 12, 0, 60);
  barlow_vn_pm = new TH1D("barlow_vn_pm_"+ID, "pm vs cent;Centrality;#Delta/#sigma_{#Delta}", 12, 0, 60);
  barlow_vn_kp = new TH1D("barlow_vn_kp_"+ID, "kp vs cent;Centrality;#Delta/#sigma_{#Delta}", 6, 0, 60);
  barlow_vn_km = new TH1D("barlow_vn_km_"+ID, "km vs cent;Centrality;#Delta/#sigma_{#Delta}", 6, 0, 60);
  barlow_vn_pr = new TH1D("barlow_vn_pr_"+ID, "pr vs cent;Centrality;#Delta/#sigma_{#Delta}", 12, 0, 60);

  // Proton y
  barlow_vn_yCM_00to10_pr = new TH1D("barlow_vn_yCM_00to10_pr_"+ID, "0-10% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_10to40_pr = new TH1D("barlow_vn_yCM_10to40_pr_"+ID, "10-40% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_40to60_pr = new TH1D("barlow_vn_yCM_40to60_pr_"+ID, "40-60% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);

  // Proton y symmetric
  barlow_vn_yCM_HADES = new TH1D("barlow_vn_yCM_HADES_"+ID, "20-30% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);

  barlow_vn_yCM_00to10_pr_symm = new TH1D("barlow_vn_yCM_00to10_pr_symm_"+ID, "0-10% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_10to40_pr_symm = new TH1D("barlow_vn_yCM_10to40_pr_symm_"+ID, "10-40% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_40to60_pr_symm = new TH1D("barlow_vn_yCM_40to60_pr_symm_"+ID, "40-60% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);

  barlow_vn_yCM_00to05_pr_symm = new TH1D("barlow_vn_yCM_00to05_pr_symm_"+ID, "0-5% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_05to10_pr_symm = new TH1D("barlow_vn_yCM_05to10_pr_symm_"+ID, "5-10% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_10to15_pr_symm = new TH1D("barlow_vn_yCM_10to15_pr_symm_"+ID, "10-15% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_15to20_pr_symm = new TH1D("barlow_vn_yCM_15to20_pr_symm_"+ID, "15-20% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_20to25_pr_symm = new TH1D("barlow_vn_yCM_20to25_pr_symm_"+ID, "20-25% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_25to30_pr_symm = new TH1D("barlow_vn_yCM_25to30_pr_symm_"+ID, "25-30% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_30to35_pr_symm = new TH1D("barlow_vn_yCM_30to35_pr_symm_"+ID, "30-35% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_35to40_pr_symm = new TH1D("barlow_vn_yCM_35to40_pr_symm_"+ID, "35-40% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_40to45_pr_symm = new TH1D("barlow_vn_yCM_40to45_pr_symm_"+ID, "40-45% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_45to50_pr_symm = new TH1D("barlow_vn_yCM_45to50_pr_symm_"+ID, "45-50% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_50to55_pr_symm = new TH1D("barlow_vn_yCM_50to55_pr_symm_"+ID, "50-55% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  barlow_vn_yCM_55to60_pr_symm = new TH1D("barlow_vn_yCM_55to60_pr_symm_"+ID, "55-60% pr vs yCM;y-y_{mid};#Delta/#sigma_{#Delta}", 20, -1, 1);
  
  // Proton pT
  barlow_vn_pT_00to10_pr = new TH1D("barlow_vn_pT_00to10_pr_"+ID, "0-10% pr vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);
  barlow_vn_pT_10to40_pr = new TH1D("barlow_vn_pT_10to40_pr_"+ID, "10-40% pr vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);
  barlow_vn_pT_40to60_pr = new TH1D("barlow_vn_pT_40to60_pr_"+ID, "40-60% pr vs p_{T};p_{T} (GeV);#Delta/#sigma_{#Delta}", 10, 0, 2);

  barlow_vn_pT_bin6_pr_symm  = new TH1D("barlow_vn_pT_bin6_pr_symm_"+ID, "bin 6;p_{T} (GeV);#Delta/#sigma_{#Delta}", 12, 0, 2.5);
  barlow_vn_pT_bin7_pr_symm  = new TH1D("barlow_vn_pT_bin7_pr_symm_"+ID, "bin 7;p_{T} (GeV);#Delta/#sigma_{#Delta}", 12, 0, 2.5);
  barlow_vn_pT_bin8_pr_symm  = new TH1D("barlow_vn_pT_bin8_pr_symm_"+ID, "bin 8;p_{T} (GeV);#Delta/#sigma_{#Delta}", 12, 0, 2.5);
  barlow_vn_pT_bin9_pr_symm  = new TH1D("barlow_vn_pT_bin9_pr_symm_"+ID, "bin 9;p_{T} (GeV);#Delta/#sigma_{#Delta}", 12, 0, 2.5);
  barlow_vn_pT_bin10_pr_symm = new TH1D("barlow_vn_pT_bin10_pr_symm_"+ID, "bin 10;p_{T} (GeV);#Delta/#sigma_{#Delta}", 12, 0, 2.5);
  barlow_vn_pT_bin11_pr_symm = new TH1D("barlow_vn_pT_bin11_pr_symm_"+ID, "bin 11;p_{T} (GeV);#Delta/#sigma_{#Delta}", 12, 0, 2.5);
  barlow_vn_pT_bin12_pr_symm = new TH1D("barlow_vn_pT_bin12_pr_symm_"+ID, "bin 12;p_{T} (GeV);#Delta/#sigma_{#Delta}", 12, 0, 2.5);
}


// Save normal and varied flow values for each point to the file supplied.
void CompositeData::addRawValuesToFile(TFile* file, TString histogramName, std::vector<DataPoint> vectorOfPoints)
{
  TH1D* temp;
  TString nameWithBinNo;
  file->cd();

  for (int i = 0; i < vectorOfPoints.size(); i++)
  {
    nameWithBinNo.Form(histogramName+"_bin%d", i+1);
    temp = new TH1D(nameWithBinNo, nameWithBinNo, 500, 1, 1);

    temp->Fill(vectorOfPoints.at(i).normalValue);
    temp->Fill(vectorOfPoints.at(i).var1Value);

    if (nVariations > 1) // at least 2 variations
    {
      temp->Fill(vectorOfPoints.at(i).var2Value);

      if (nVariations > 2) // at least 3 variations
      {
        temp->Fill(vectorOfPoints.at(i).var3Value);

        if (nVariations > 3) // 4 variations
        {
          temp->Fill(vectorOfPoints.at(i).var4Value);
        }
      }
    }
    
    temp->Write();
    delete temp;
  }
}


void CompositeData::addBarlowValuesToFile(TFile* file, TH1D* barlowHistogram, std::vector<DataPoint> vectorOfPoints)
{
  for (int i = 1; i <= vectorOfPoints.size(); i++)
    { barlowHistogram->SetBinContent(i, vectorOfPoints.at(i-1).deltaByDeltaError); }

  file->cd();
  barlowHistogram->Write();
}


// Save each flow value along with the same value from one variation to see the effect on every point.
//Also save the values used for the Barlow check.
//Normal Variation is only used to get names of histograms.
void CompositeData::saveDetails(Variation* normalData)
{
  TFile *detailsFile = new TFile("Details_"+ID+".root", "RECREATE");

  // Save the flow values 
  addRawValuesToFile(detailsFile, normalData->h_vn_pp->GetName(), v_vn_pp);
  addRawValuesToFile(detailsFile, normalData->h_vn_pm->GetName(), v_vn_pm);
  addRawValuesToFile(detailsFile, normalData->h_vn_kp->GetName(), v_vn_kp);
  addRawValuesToFile(detailsFile, normalData->h_vn_km->GetName(), v_vn_km);
  addRawValuesToFile(detailsFile, normalData->h_vn_pr->GetName(), v_vn_pr);

  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_00to10_pr->GetName(), v_vn_yCM_00to10_pr);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_10to40_pr->GetName(), v_vn_yCM_10to40_pr);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_40to60_pr->GetName(), v_vn_yCM_40to60_pr);

  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_HADES->GetName(), v_vn_yCM_HADES);

  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_00to10_pr_symm->GetName(), v_vn_yCM_00to10_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_10to40_pr_symm->GetName(), v_vn_yCM_10to40_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_40to60_pr_symm->GetName(), v_vn_yCM_40to60_pr_symm);

  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_00to05_pr_symm->GetName(), v_vn_yCM_00to05_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_05to10_pr_symm->GetName(), v_vn_yCM_05to10_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_10to15_pr_symm->GetName(), v_vn_yCM_10to15_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_15to20_pr_symm->GetName(), v_vn_yCM_15to20_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_20to25_pr_symm->GetName(), v_vn_yCM_20to25_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_25to30_pr_symm->GetName(), v_vn_yCM_25to30_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_30to35_pr_symm->GetName(), v_vn_yCM_30to35_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_35to40_pr_symm->GetName(), v_vn_yCM_35to40_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_40to45_pr_symm->GetName(), v_vn_yCM_40to45_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_45to50_pr_symm->GetName(), v_vn_yCM_45to50_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_50to55_pr_symm->GetName(), v_vn_yCM_50to55_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_yCM_55to60_pr_symm->GetName(), v_vn_yCM_55to60_pr_symm);

  addRawValuesToFile(detailsFile, normalData->h_vn_pT_00to10_pr->GetName(), v_vn_pT_00to10_pr);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_10to40_pr->GetName(), v_vn_pT_10to40_pr);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_40to60_pr->GetName(), v_vn_pT_40to60_pr);

  addRawValuesToFile(detailsFile, normalData->h_vn_pT_bin6_10to40_pr_symm->GetName(), v_vn_pT_bin6_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_bin7_10to40_pr_symm->GetName(), v_vn_pT_bin7_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_bin8_10to40_pr_symm->GetName(), v_vn_pT_bin8_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_bin9_10to40_pr_symm->GetName(), v_vn_pT_bin9_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_bin10_10to40_pr_symm->GetName(), v_vn_pT_bin10_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_bin11_10to40_pr_symm->GetName(), v_vn_pT_bin11_pr_symm);
  addRawValuesToFile(detailsFile, normalData->h_vn_pT_bin12_10to40_pr_symm->GetName(), v_vn_pT_bin12_pr_symm);
  ////

  // Save the significance values
  addBarlowValuesToFile(detailsFile, barlow_vn_pp, v_vn_pp);
  addBarlowValuesToFile(detailsFile, barlow_vn_pm, v_vn_pm);
  addBarlowValuesToFile(detailsFile, barlow_vn_kp, v_vn_kp);
  addBarlowValuesToFile(detailsFile, barlow_vn_km, v_vn_km);
  addBarlowValuesToFile(detailsFile, barlow_vn_pr, v_vn_pr);
  
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_00to10_pr, v_vn_yCM_00to10_pr);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_10to40_pr, v_vn_yCM_10to40_pr);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_40to60_pr, v_vn_yCM_40to60_pr);

  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_HADES, v_vn_yCM_HADES);

  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_00to10_pr_symm, v_vn_yCM_00to10_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_10to40_pr_symm, v_vn_yCM_10to40_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_40to60_pr_symm, v_vn_yCM_40to60_pr_symm);

  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_00to05_pr_symm, v_vn_yCM_00to05_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_05to10_pr_symm, v_vn_yCM_05to10_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_10to15_pr_symm, v_vn_yCM_10to15_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_15to20_pr_symm, v_vn_yCM_15to20_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_20to25_pr_symm, v_vn_yCM_20to25_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_25to30_pr_symm, v_vn_yCM_25to30_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_30to35_pr_symm, v_vn_yCM_30to35_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_35to40_pr_symm, v_vn_yCM_35to40_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_40to45_pr_symm, v_vn_yCM_40to45_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_45to50_pr_symm, v_vn_yCM_45to50_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_50to55_pr_symm, v_vn_yCM_50to55_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_yCM_55to60_pr_symm, v_vn_yCM_55to60_pr_symm);

  addBarlowValuesToFile(detailsFile, barlow_vn_pT_00to10_pr, v_vn_pT_00to10_pr);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_10to40_pr, v_vn_pT_10to40_pr);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_40to60_pr, v_vn_pT_40to60_pr);

  addBarlowValuesToFile(detailsFile, barlow_vn_pT_bin6_pr_symm, v_vn_pT_bin6_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_bin7_pr_symm, v_vn_pT_bin7_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_bin8_pr_symm, v_vn_pT_bin8_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_bin9_pr_symm, v_vn_pT_bin9_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_bin10_pr_symm, v_vn_pT_bin10_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_bin11_pr_symm, v_vn_pT_bin11_pr_symm);
  addBarlowValuesToFile(detailsFile, barlow_vn_pT_bin12_pr_symm, v_vn_pT_bin12_pr_symm);
  ////
  
  detailsFile->Close();
  delete detailsFile;
}


// Calculate necessary info for one histogram type with normal data and one variation
void CompositeData::mergePoints(TH1D* normalHisto, TH1D* var1Histo, std::vector<DataPoint>& vectorOfPoints)
{
  Double_t N = (Double_t)nVariations + 1.0; // Number of measurements (variations + 1 from the nominal results)
  
  for (int i = 1; i <= normalHisto->GetNbinsX(); i++)
  {
    DataPoint point;
    point.normalValue = normalHisto->GetBinContent(i);
    point.normalError = normalHisto->GetBinError(i);
    
    point.var1Value = var1Histo->GetBinContent(i);
    point.var1Error = var1Histo->GetBinError(i);

    if((point.normalValue == 0.0 && point.normalError == 0.0) ||
       (point.var1Value   == 0.0 && point.var1Error   == 0.0))
      {
	point.delta = 0.0;
	point.deltaError = 0.0;
	point.deltaByDeltaError = 0.0;
	point.mean = 0.0;
	point.variance = 0.0;
	point.stdDev = 0.0;
	point.stdDevPercentage = 0.0;
      }
    else
      {
	point.delta = TMath::Abs(point.var1Value - point.normalValue);
	point.deltaError = TMath::Sqrt(TMath::Abs(TMath::Power(point.var1Error, 2) - TMath::Power(point.normalError, 2)));
	point.deltaByDeltaError = (point.deltaError == 0.0)?0.0:point.delta/point.deltaError;
	point.mean = (point.normalValue + point.var1Value) / N;
	
	Double_t Nvariance = TMath::Power(point.normalValue - point.mean, 2.0) + 
	  TMath::Power(point.var1Value - point.mean, 2.0);

	point.variance = Nvariance / N;
	point.stdDev = TMath::Sqrt(point.variance);
	point.stdDevPercentage = 100.0 * point.stdDev / TMath::Abs(point.normalValue);
      }
    vectorOfPoints.push_back(point);
  }
}


// Calculate necessary info for one histogram type with normal data and two variations
void CompositeData::mergePoints(TH1D* normalHisto, TH1D* var1Histo, TH1D* var2Histo, std::vector<DataPoint>& vectorOfPoints)
{
  Double_t N = (Double_t)nVariations + 1.0; // Number of measurements (variations + 1 from the nominal results)
  
  for (int i = 1; i <= normalHisto->GetNbinsX(); i++)
  {
    DataPoint point;
    point.normalValue = normalHisto->GetBinContent(i);
    point.normalError = normalHisto->GetBinError(i);
    
    point.var1Value = var1Histo->GetBinContent(i);
    point.var1Error = var1Histo->GetBinError(i);

    point.var2Value = var2Histo->GetBinContent(i);
    point.var2Error = var2Histo->GetBinError(i);

    point.getMax();
    point.getMin();

    if((point.normalValue == 0.0 && point.normalError == 0.0) ||
       (point.var1Value   == 0.0 && point.var1Error   == 0.0) ||
       (point.var2Value   == 0.0 && point.var2Error   == 0.0))
      {
	point.delta = 0.0;
	point.deltaError = 0.0;
	point.deltaByDeltaError = 0.0;
	point.mean = 0.0;
	point.variance = 0.0;
	point.stdDev = 0.0;
	point.stdDevPercentage = 0.0;
      }
    else
      {
	point.delta = TMath::Abs(point.maxValue - point.minValue);
	point.deltaError = TMath::Sqrt(TMath::Abs(TMath::Power(point.maxError, 2) - TMath::Power(point.minError, 2)));
	point.deltaByDeltaError = (point.deltaError == 0.0)?0.0:point.delta/point.deltaError; 
	point.mean = (point.normalValue + point.var1Value + point.var2Value) / N;
    
	Double_t Nvariance = TMath::Power(point.normalValue - point.mean, 2.0) + 
	  TMath::Power(point.var1Value - point.mean, 2.0) + 
	  TMath::Power(point.var2Value - point.mean, 2.0);

	point.variance = Nvariance / N;
	point.stdDev = TMath::Sqrt(point.variance);
	point.stdDevPercentage = 100.0 * point.stdDev / TMath::Abs(point.normalValue);
      }
    vectorOfPoints.push_back(point);
  }
}


// Calculate necessary info for one histogram type with normal data and three variations
void CompositeData::mergePoints(TH1D* normalHisto, TH1D* var1Histo, TH1D* var2Histo, TH1D* var3Histo, std::vector<DataPoint>& vectorOfPoints)
{
  Double_t N = (Double_t)nVariations + 1.0; // Number of measurements (variations + 1 from the nominal results)
  
 for (int i = 1; i <= normalHisto->GetNbinsX(); i++)
  {
    DataPoint point;
    point.normalValue = normalHisto->GetBinContent(i);
    point.normalError = normalHisto->GetBinError(i);
    
    point.var1Value = var1Histo->GetBinContent(i);
    point.var1Error = var1Histo->GetBinError(i);

    point.var2Value = var2Histo->GetBinContent(i);
    point.var2Error = var2Histo->GetBinError(i);

    point.var3Value = var3Histo->GetBinContent(i);
    point.var3Error = var3Histo->GetBinError(i);

    point.getMax();
    point.getMin();

    if((point.normalValue == 0.0 && point.normalError == 0.0) ||
       (point.var1Value   == 0.0 && point.var1Error   == 0.0) ||
       (point.var2Value   == 0.0 && point.var2Error   == 0.0) ||
       (point.var3Value   == 0.0 && point.var3Error   == 0.0))
      {
	point.delta = 0.0;
	point.deltaError = 0.0;
	point.deltaByDeltaError = 0.0;
	point.mean = 0.0;
	point.variance = 0.0;
	point.stdDev = 0.0;
	point.stdDevPercentage = 0.0;
      }
    else
      {
	point.delta = TMath::Abs(point.maxValue - point.minValue);
	point.deltaError = TMath::Sqrt(TMath::Abs(TMath::Power(point.maxError, 2) - TMath::Power(point.minError, 2)));
	point.deltaByDeltaError = (point.deltaError == 0.0)?0.0:point.delta/point.deltaError;
	point.mean = (point.normalValue + point.var1Value + point.var2Value + point.var3Value) / N;
    
	Double_t Nvariance = TMath::Power(point.normalValue - point.mean, 2.0) + 
	  TMath::Power(point.var1Value - point.mean, 2.0) + 
	  TMath::Power(point.var2Value - point.mean, 2.0) + 
	  TMath::Power(point.var3Value - point.mean, 2.0);

	point.variance = Nvariance / N;
	point.stdDev = TMath::Sqrt(point.variance);
	point.stdDevPercentage = 100.0 * point.stdDev / TMath::Abs(point.normalValue);
      }
    vectorOfPoints.push_back(point);
  }
}

// Calculate necessary info for one histogram type with normal data and four variations
void CompositeData::mergePoints(TH1D* normalHisto, TH1D* var1Histo, TH1D* var2Histo, TH1D* var3Histo, TH1D* var4Histo, std::vector<DataPoint>& vectorOfPoints)
{
  Double_t N = (Double_t)nVariations + 1.0; // Number of measurements (variations + 1 from the nominal results)
  
  for (int i = 1; i <= normalHisto->GetNbinsX(); i++)
  {
    DataPoint point;
    point.normalValue = normalHisto->GetBinContent(i);
    point.normalError = normalHisto->GetBinError(i);
    
    point.var1Value = var1Histo->GetBinContent(i);
    point.var1Error = var1Histo->GetBinError(i);

    point.var2Value = var2Histo->GetBinContent(i);
    point.var2Error = var2Histo->GetBinError(i);

    point.var3Value = var3Histo->GetBinContent(i);
    point.var3Error = var3Histo->GetBinError(i);

    point.var4Value = var4Histo->GetBinContent(i);
    point.var4Error = var4Histo->GetBinError(i);

    point.getMax();
    point.getMin();

    if((point.normalValue == 0.0 && point.normalError == 0.0) ||
       (point.var1Value   == 0.0 && point.var1Error   == 0.0) ||
       (point.var2Value   == 0.0 && point.var2Error   == 0.0) ||
       (point.var3Value   == 0.0 && point.var3Error   == 0.0) ||
       (point.var4Value   == 0.0 && point.var4Error   == 0.0))
      {
	point.delta = 0.0;
	point.deltaError = 0.0;
	point.deltaByDeltaError = 0.0;
	point.mean = 0.0;
	point.variance = 0.0;
	point.stdDev = 0.0;
	point.stdDevPercentage = 0.0;
      }
    else
      {
	point.delta = TMath::Abs(point.maxValue - point.minValue);
	point.deltaError = TMath::Sqrt(TMath::Abs(TMath::Power(point.maxError, 2) - TMath::Power(point.minError, 2)));
	point.deltaByDeltaError = (point.deltaError == 0.0)?0.0:point.delta/point.deltaError;
	point.mean = (point.normalValue + point.var1Value + point.var2Value + point.var3Value + point.var4Value) / N;
    
	Double_t Nvariance = TMath::Power(point.normalValue - point.mean, 2.0) + 
	  TMath::Power(point.var1Value - point.mean, 2.0) + 
	  TMath::Power(point.var2Value - point.mean, 2.0) + 
	  TMath::Power(point.var3Value - point.mean, 2.0) +
	  TMath::Power(point.var4Value - point.mean, 2.0);

	point.variance = Nvariance / N;
	point.stdDev = TMath::Sqrt(point.variance);
	point.stdDevPercentage = 100.0 * point.stdDev / TMath::Abs(point.normalValue);
      }
    vectorOfPoints.push_back(point);
  }
}


// Combine one variation with the normal data to get the attributes for systematics
void CompositeData::combine(Variation* normalData, Variation* var1Data)
{
  mergePoints(normalData->h_vn_pp, var1Data->h_vn_pp, v_vn_pp);
  mergePoints(normalData->h_vn_pm, var1Data->h_vn_pm, v_vn_pm);
  mergePoints(normalData->h_vn_kp, var1Data->h_vn_kp, v_vn_kp);
  mergePoints(normalData->h_vn_km, var1Data->h_vn_km, v_vn_km);
  mergePoints(normalData->h_vn_pr, var1Data->h_vn_pr, v_vn_pr);
  
  mergePoints(normalData->h_vn_yCM_00to10_pr, var1Data->h_vn_yCM_00to10_pr, v_vn_yCM_00to10_pr);
  mergePoints(normalData->h_vn_yCM_10to40_pr, var1Data->h_vn_yCM_10to40_pr, v_vn_yCM_10to40_pr);
  mergePoints(normalData->h_vn_yCM_40to60_pr, var1Data->h_vn_yCM_40to60_pr, v_vn_yCM_40to60_pr);

  mergePoints(normalData->h_vn_yCM_HADES, var1Data->h_vn_yCM_HADES, v_vn_yCM_HADES);

  mergePoints(normalData->h_vn_yCM_00to10_pr_symm, var1Data->h_vn_yCM_00to10_pr_symm, v_vn_yCM_00to10_pr_symm);
  mergePoints(normalData->h_vn_yCM_10to40_pr_symm, var1Data->h_vn_yCM_10to40_pr_symm, v_vn_yCM_10to40_pr_symm);
  mergePoints(normalData->h_vn_yCM_40to60_pr_symm, var1Data->h_vn_yCM_40to60_pr_symm, v_vn_yCM_40to60_pr_symm);

  mergePoints(normalData->h_vn_yCM_00to05_pr_symm, var1Data->h_vn_yCM_00to05_pr_symm, v_vn_yCM_00to05_pr_symm);
  mergePoints(normalData->h_vn_yCM_05to10_pr_symm, var1Data->h_vn_yCM_05to10_pr_symm, v_vn_yCM_05to10_pr_symm);
  mergePoints(normalData->h_vn_yCM_10to15_pr_symm, var1Data->h_vn_yCM_10to15_pr_symm, v_vn_yCM_10to15_pr_symm);
  mergePoints(normalData->h_vn_yCM_15to20_pr_symm, var1Data->h_vn_yCM_15to20_pr_symm, v_vn_yCM_15to20_pr_symm);
  mergePoints(normalData->h_vn_yCM_20to25_pr_symm, var1Data->h_vn_yCM_20to25_pr_symm, v_vn_yCM_20to25_pr_symm);
  mergePoints(normalData->h_vn_yCM_25to30_pr_symm, var1Data->h_vn_yCM_25to30_pr_symm, v_vn_yCM_25to30_pr_symm);
  mergePoints(normalData->h_vn_yCM_30to35_pr_symm, var1Data->h_vn_yCM_30to35_pr_symm, v_vn_yCM_30to35_pr_symm);
  mergePoints(normalData->h_vn_yCM_35to40_pr_symm, var1Data->h_vn_yCM_35to40_pr_symm, v_vn_yCM_35to40_pr_symm);
  mergePoints(normalData->h_vn_yCM_40to45_pr_symm, var1Data->h_vn_yCM_40to45_pr_symm, v_vn_yCM_40to45_pr_symm);
  mergePoints(normalData->h_vn_yCM_45to50_pr_symm, var1Data->h_vn_yCM_45to50_pr_symm, v_vn_yCM_45to50_pr_symm);
  mergePoints(normalData->h_vn_yCM_50to55_pr_symm, var1Data->h_vn_yCM_50to55_pr_symm, v_vn_yCM_50to55_pr_symm);
  mergePoints(normalData->h_vn_yCM_55to60_pr_symm, var1Data->h_vn_yCM_55to60_pr_symm, v_vn_yCM_55to60_pr_symm);

  mergePoints(normalData->h_vn_pT_00to10_pr, var1Data->h_vn_pT_00to10_pr, v_vn_pT_00to10_pr);
  mergePoints(normalData->h_vn_pT_10to40_pr, var1Data->h_vn_pT_10to40_pr, v_vn_pT_10to40_pr);
  mergePoints(normalData->h_vn_pT_40to60_pr, var1Data->h_vn_pT_40to60_pr, v_vn_pT_40to60_pr);

  mergePoints(normalData->h_vn_pT_bin6_10to40_pr_symm,  var1Data->h_vn_pT_bin6_10to40_pr_symm,  v_vn_pT_bin6_pr_symm);
  mergePoints(normalData->h_vn_pT_bin7_10to40_pr_symm,  var1Data->h_vn_pT_bin7_10to40_pr_symm,  v_vn_pT_bin7_pr_symm);
  mergePoints(normalData->h_vn_pT_bin8_10to40_pr_symm,  var1Data->h_vn_pT_bin8_10to40_pr_symm,  v_vn_pT_bin8_pr_symm);
  mergePoints(normalData->h_vn_pT_bin9_10to40_pr_symm,  var1Data->h_vn_pT_bin9_10to40_pr_symm,  v_vn_pT_bin9_pr_symm);
  mergePoints(normalData->h_vn_pT_bin10_10to40_pr_symm, var1Data->h_vn_pT_bin10_10to40_pr_symm, v_vn_pT_bin10_pr_symm);
  mergePoints(normalData->h_vn_pT_bin11_10to40_pr_symm, var1Data->h_vn_pT_bin11_10to40_pr_symm, v_vn_pT_bin11_pr_symm);
  mergePoints(normalData->h_vn_pT_bin12_10to40_pr_symm, var1Data->h_vn_pT_bin12_10to40_pr_symm, v_vn_pT_bin12_pr_symm);
}



// Combine two variations with the normal data to get the attributes for systematics
void CompositeData::combine(Variation* normalData, Variation* var1Data, Variation* var2Data)
{
  mergePoints(normalData->h_vn_pp, var1Data->h_vn_pp, var2Data->h_vn_pp, v_vn_pp);
  mergePoints(normalData->h_vn_pm, var1Data->h_vn_pm, var2Data->h_vn_pm, v_vn_pm);
  mergePoints(normalData->h_vn_kp, var1Data->h_vn_kp, var2Data->h_vn_kp, v_vn_kp);
  mergePoints(normalData->h_vn_km, var1Data->h_vn_km, var2Data->h_vn_km, v_vn_km);
  mergePoints(normalData->h_vn_pr, var1Data->h_vn_pr, var2Data->h_vn_pr, v_vn_pr);
  
  mergePoints(normalData->h_vn_yCM_00to10_pr, var1Data->h_vn_yCM_00to10_pr, var2Data->h_vn_yCM_00to10_pr, v_vn_yCM_00to10_pr);
  mergePoints(normalData->h_vn_yCM_10to40_pr, var1Data->h_vn_yCM_10to40_pr, var2Data->h_vn_yCM_10to40_pr, v_vn_yCM_10to40_pr);
  mergePoints(normalData->h_vn_yCM_40to60_pr, var1Data->h_vn_yCM_40to60_pr, var2Data->h_vn_yCM_40to60_pr, v_vn_yCM_40to60_pr);

  mergePoints(normalData->h_vn_yCM_HADES, var1Data->h_vn_yCM_HADES, var2Data->h_vn_yCM_HADES, v_vn_yCM_HADES);

  mergePoints(normalData->h_vn_yCM_00to10_pr_symm, var1Data->h_vn_yCM_00to10_pr_symm, var2Data->h_vn_yCM_00to10_pr_symm, v_vn_yCM_00to10_pr_symm);
  mergePoints(normalData->h_vn_yCM_10to40_pr_symm, var1Data->h_vn_yCM_10to40_pr_symm, var2Data->h_vn_yCM_10to40_pr_symm, v_vn_yCM_10to40_pr_symm);
  mergePoints(normalData->h_vn_yCM_40to60_pr_symm, var1Data->h_vn_yCM_40to60_pr_symm, var2Data->h_vn_yCM_40to60_pr_symm, v_vn_yCM_40to60_pr_symm);

  mergePoints(normalData->h_vn_yCM_00to05_pr_symm, var1Data->h_vn_yCM_00to05_pr_symm, var2Data->h_vn_yCM_00to05_pr_symm, v_vn_yCM_00to05_pr_symm);
  mergePoints(normalData->h_vn_yCM_05to10_pr_symm, var1Data->h_vn_yCM_05to10_pr_symm, var2Data->h_vn_yCM_05to10_pr_symm, v_vn_yCM_05to10_pr_symm);
  mergePoints(normalData->h_vn_yCM_10to15_pr_symm, var1Data->h_vn_yCM_10to15_pr_symm, var2Data->h_vn_yCM_10to15_pr_symm, v_vn_yCM_10to15_pr_symm);
  mergePoints(normalData->h_vn_yCM_15to20_pr_symm, var1Data->h_vn_yCM_15to20_pr_symm, var2Data->h_vn_yCM_15to20_pr_symm, v_vn_yCM_15to20_pr_symm);
  mergePoints(normalData->h_vn_yCM_20to25_pr_symm, var1Data->h_vn_yCM_20to25_pr_symm, var2Data->h_vn_yCM_20to25_pr_symm, v_vn_yCM_20to25_pr_symm);
  mergePoints(normalData->h_vn_yCM_25to30_pr_symm, var1Data->h_vn_yCM_25to30_pr_symm, var2Data->h_vn_yCM_25to30_pr_symm, v_vn_yCM_25to30_pr_symm);
  mergePoints(normalData->h_vn_yCM_30to35_pr_symm, var1Data->h_vn_yCM_30to35_pr_symm, var2Data->h_vn_yCM_30to35_pr_symm, v_vn_yCM_30to35_pr_symm);
  mergePoints(normalData->h_vn_yCM_35to40_pr_symm, var1Data->h_vn_yCM_35to40_pr_symm, var2Data->h_vn_yCM_35to40_pr_symm, v_vn_yCM_35to40_pr_symm);
  mergePoints(normalData->h_vn_yCM_40to45_pr_symm, var1Data->h_vn_yCM_40to45_pr_symm, var2Data->h_vn_yCM_40to45_pr_symm, v_vn_yCM_40to45_pr_symm);
  mergePoints(normalData->h_vn_yCM_45to50_pr_symm, var1Data->h_vn_yCM_45to50_pr_symm, var2Data->h_vn_yCM_45to50_pr_symm, v_vn_yCM_45to50_pr_symm);
  mergePoints(normalData->h_vn_yCM_50to55_pr_symm, var1Data->h_vn_yCM_50to55_pr_symm, var2Data->h_vn_yCM_50to55_pr_symm, v_vn_yCM_50to55_pr_symm);
  mergePoints(normalData->h_vn_yCM_55to60_pr_symm, var1Data->h_vn_yCM_55to60_pr_symm, var2Data->h_vn_yCM_55to60_pr_symm, v_vn_yCM_55to60_pr_symm);

  mergePoints(normalData->h_vn_pT_00to10_pr, var1Data->h_vn_pT_00to10_pr, var2Data->h_vn_pT_00to10_pr, v_vn_pT_00to10_pr);
  mergePoints(normalData->h_vn_pT_10to40_pr, var1Data->h_vn_pT_10to40_pr, var2Data->h_vn_pT_10to40_pr, v_vn_pT_10to40_pr);
  mergePoints(normalData->h_vn_pT_40to60_pr, var1Data->h_vn_pT_40to60_pr, var2Data->h_vn_pT_40to60_pr, v_vn_pT_40to60_pr);

  mergePoints(normalData->h_vn_pT_bin6_10to40_pr_symm,  var1Data->h_vn_pT_bin6_10to40_pr_symm,  var2Data->h_vn_pT_bin6_10to40_pr_symm,  v_vn_pT_bin6_pr_symm);
  mergePoints(normalData->h_vn_pT_bin7_10to40_pr_symm,  var1Data->h_vn_pT_bin7_10to40_pr_symm,  var2Data->h_vn_pT_bin7_10to40_pr_symm,  v_vn_pT_bin7_pr_symm);
  mergePoints(normalData->h_vn_pT_bin8_10to40_pr_symm,  var1Data->h_vn_pT_bin8_10to40_pr_symm,  var2Data->h_vn_pT_bin8_10to40_pr_symm,  v_vn_pT_bin8_pr_symm);
  mergePoints(normalData->h_vn_pT_bin9_10to40_pr_symm,  var1Data->h_vn_pT_bin9_10to40_pr_symm,  var2Data->h_vn_pT_bin9_10to40_pr_symm,  v_vn_pT_bin9_pr_symm);
  mergePoints(normalData->h_vn_pT_bin10_10to40_pr_symm, var1Data->h_vn_pT_bin10_10to40_pr_symm, var2Data->h_vn_pT_bin10_10to40_pr_symm, v_vn_pT_bin10_pr_symm);
  mergePoints(normalData->h_vn_pT_bin11_10to40_pr_symm, var1Data->h_vn_pT_bin11_10to40_pr_symm, var2Data->h_vn_pT_bin11_10to40_pr_symm, v_vn_pT_bin11_pr_symm);
  mergePoints(normalData->h_vn_pT_bin12_10to40_pr_symm, var1Data->h_vn_pT_bin12_10to40_pr_symm, var2Data->h_vn_pT_bin12_10to40_pr_symm, v_vn_pT_bin12_pr_symm);
}


// Combine three variations with the normal data to get the attributes for systematics
void CompositeData::combine(Variation* normalData, Variation* var1Data, Variation* var2Data, Variation* var3Data)
{
  mergePoints(normalData->h_vn_pp, var1Data->h_vn_pp, var2Data->h_vn_pp, var3Data->h_vn_pp, v_vn_pp);
  mergePoints(normalData->h_vn_pm, var1Data->h_vn_pm, var2Data->h_vn_pm, var3Data->h_vn_pm, v_vn_pm);
  mergePoints(normalData->h_vn_kp, var1Data->h_vn_kp, var2Data->h_vn_kp, var3Data->h_vn_kp, v_vn_kp);
  mergePoints(normalData->h_vn_km, var1Data->h_vn_km, var2Data->h_vn_km, var3Data->h_vn_km, v_vn_km);
  mergePoints(normalData->h_vn_pr, var1Data->h_vn_pr, var2Data->h_vn_pr, var3Data->h_vn_pr, v_vn_pr);
  
  mergePoints(normalData->h_vn_yCM_00to10_pr, var1Data->h_vn_yCM_00to10_pr, var2Data->h_vn_yCM_00to10_pr, var3Data->h_vn_yCM_00to10_pr, v_vn_yCM_00to10_pr);
  mergePoints(normalData->h_vn_yCM_10to40_pr, var1Data->h_vn_yCM_10to40_pr, var2Data->h_vn_yCM_10to40_pr, var3Data->h_vn_yCM_10to40_pr, v_vn_yCM_10to40_pr);
  mergePoints(normalData->h_vn_yCM_40to60_pr, var1Data->h_vn_yCM_40to60_pr, var2Data->h_vn_yCM_40to60_pr, var3Data->h_vn_yCM_40to60_pr, v_vn_yCM_40to60_pr);

  mergePoints(normalData->h_vn_yCM_HADES, var1Data->h_vn_yCM_HADES, var2Data->h_vn_yCM_HADES, var3Data->h_vn_yCM_HADES, v_vn_yCM_HADES);

  mergePoints(normalData->h_vn_yCM_00to10_pr_symm, var1Data->h_vn_yCM_00to10_pr_symm, var2Data->h_vn_yCM_00to10_pr_symm, var3Data->h_vn_yCM_00to10_pr_symm, v_vn_yCM_00to10_pr_symm);
  mergePoints(normalData->h_vn_yCM_10to40_pr_symm, var1Data->h_vn_yCM_10to40_pr_symm, var2Data->h_vn_yCM_10to40_pr_symm, var3Data->h_vn_yCM_10to40_pr_symm, v_vn_yCM_10to40_pr_symm);
  mergePoints(normalData->h_vn_yCM_40to60_pr_symm, var1Data->h_vn_yCM_40to60_pr_symm, var2Data->h_vn_yCM_40to60_pr_symm, var3Data->h_vn_yCM_40to60_pr_symm, v_vn_yCM_40to60_pr_symm);

  mergePoints(normalData->h_vn_yCM_00to05_pr_symm, var1Data->h_vn_yCM_00to05_pr_symm, var2Data->h_vn_yCM_00to05_pr_symm, var3Data->h_vn_yCM_00to05_pr_symm, v_vn_yCM_00to05_pr_symm);
  mergePoints(normalData->h_vn_yCM_05to10_pr_symm, var1Data->h_vn_yCM_05to10_pr_symm, var2Data->h_vn_yCM_05to10_pr_symm, var3Data->h_vn_yCM_05to10_pr_symm, v_vn_yCM_05to10_pr_symm);
  mergePoints(normalData->h_vn_yCM_10to15_pr_symm, var1Data->h_vn_yCM_10to15_pr_symm, var2Data->h_vn_yCM_10to15_pr_symm, var3Data->h_vn_yCM_10to15_pr_symm, v_vn_yCM_10to15_pr_symm);
  mergePoints(normalData->h_vn_yCM_15to20_pr_symm, var1Data->h_vn_yCM_15to20_pr_symm, var2Data->h_vn_yCM_15to20_pr_symm, var3Data->h_vn_yCM_15to20_pr_symm, v_vn_yCM_15to20_pr_symm);
  mergePoints(normalData->h_vn_yCM_20to25_pr_symm, var1Data->h_vn_yCM_20to25_pr_symm, var2Data->h_vn_yCM_20to25_pr_symm, var3Data->h_vn_yCM_20to25_pr_symm, v_vn_yCM_20to25_pr_symm);
  mergePoints(normalData->h_vn_yCM_25to30_pr_symm, var1Data->h_vn_yCM_25to30_pr_symm, var2Data->h_vn_yCM_25to30_pr_symm, var3Data->h_vn_yCM_25to30_pr_symm, v_vn_yCM_25to30_pr_symm);
  mergePoints(normalData->h_vn_yCM_30to35_pr_symm, var1Data->h_vn_yCM_30to35_pr_symm, var2Data->h_vn_yCM_30to35_pr_symm, var3Data->h_vn_yCM_30to35_pr_symm, v_vn_yCM_30to35_pr_symm);
  mergePoints(normalData->h_vn_yCM_35to40_pr_symm, var1Data->h_vn_yCM_35to40_pr_symm, var2Data->h_vn_yCM_35to40_pr_symm, var3Data->h_vn_yCM_35to40_pr_symm, v_vn_yCM_35to40_pr_symm);
  mergePoints(normalData->h_vn_yCM_40to45_pr_symm, var1Data->h_vn_yCM_40to45_pr_symm, var2Data->h_vn_yCM_40to45_pr_symm, var3Data->h_vn_yCM_40to45_pr_symm, v_vn_yCM_40to45_pr_symm);
  mergePoints(normalData->h_vn_yCM_45to50_pr_symm, var1Data->h_vn_yCM_45to50_pr_symm, var2Data->h_vn_yCM_45to50_pr_symm, var3Data->h_vn_yCM_45to50_pr_symm, v_vn_yCM_45to50_pr_symm);
  mergePoints(normalData->h_vn_yCM_50to55_pr_symm, var1Data->h_vn_yCM_50to55_pr_symm, var2Data->h_vn_yCM_50to55_pr_symm, var3Data->h_vn_yCM_50to55_pr_symm, v_vn_yCM_50to55_pr_symm);
  mergePoints(normalData->h_vn_yCM_55to60_pr_symm, var1Data->h_vn_yCM_55to60_pr_symm, var2Data->h_vn_yCM_55to60_pr_symm, var3Data->h_vn_yCM_55to60_pr_symm, v_vn_yCM_55to60_pr_symm);

  mergePoints(normalData->h_vn_pT_00to10_pr, var1Data->h_vn_pT_00to10_pr, var2Data->h_vn_pT_00to10_pr, var3Data->h_vn_pT_00to10_pr, v_vn_pT_00to10_pr);
  mergePoints(normalData->h_vn_pT_10to40_pr, var1Data->h_vn_pT_10to40_pr, var2Data->h_vn_pT_10to40_pr, var3Data->h_vn_pT_10to40_pr, v_vn_pT_10to40_pr);
  mergePoints(normalData->h_vn_pT_40to60_pr, var1Data->h_vn_pT_40to60_pr, var2Data->h_vn_pT_40to60_pr, var3Data->h_vn_pT_40to60_pr, v_vn_pT_40to60_pr);

  mergePoints(normalData->h_vn_pT_bin6_10to40_pr_symm,  var1Data->h_vn_pT_bin6_10to40_pr_symm,  var2Data->h_vn_pT_bin6_10to40_pr_symm,   var3Data->h_vn_pT_bin6_10to40_pr_symm,  v_vn_pT_bin6_pr_symm);
  mergePoints(normalData->h_vn_pT_bin7_10to40_pr_symm,  var1Data->h_vn_pT_bin7_10to40_pr_symm,  var2Data->h_vn_pT_bin7_10to40_pr_symm,   var3Data->h_vn_pT_bin7_10to40_pr_symm,  v_vn_pT_bin7_pr_symm);
  mergePoints(normalData->h_vn_pT_bin8_10to40_pr_symm,  var1Data->h_vn_pT_bin8_10to40_pr_symm,  var2Data->h_vn_pT_bin8_10to40_pr_symm,   var3Data->h_vn_pT_bin8_10to40_pr_symm,  v_vn_pT_bin8_pr_symm);
  mergePoints(normalData->h_vn_pT_bin9_10to40_pr_symm,  var1Data->h_vn_pT_bin9_10to40_pr_symm,  var2Data->h_vn_pT_bin9_10to40_pr_symm,   var3Data->h_vn_pT_bin9_10to40_pr_symm,  v_vn_pT_bin9_pr_symm);
  mergePoints(normalData->h_vn_pT_bin10_10to40_pr_symm, var1Data->h_vn_pT_bin10_10to40_pr_symm, var2Data->h_vn_pT_bin10_10to40_pr_symm,  var3Data->h_vn_pT_bin10_10to40_pr_symm, v_vn_pT_bin10_pr_symm);
  mergePoints(normalData->h_vn_pT_bin11_10to40_pr_symm, var1Data->h_vn_pT_bin11_10to40_pr_symm, var2Data->h_vn_pT_bin11_10to40_pr_symm,  var3Data->h_vn_pT_bin11_10to40_pr_symm, v_vn_pT_bin11_pr_symm);
  mergePoints(normalData->h_vn_pT_bin12_10to40_pr_symm, var1Data->h_vn_pT_bin12_10to40_pr_symm, var2Data->h_vn_pT_bin12_10to40_pr_symm,  var3Data->h_vn_pT_bin12_10to40_pr_symm, v_vn_pT_bin12_pr_symm);
}



// Combine four variations with the normal data to get the attributes for systematics
void CompositeData::combine(Variation* normalData, Variation* var1Data, Variation* var2Data, Variation* var3Data, Variation* var4Data)
{
  mergePoints(normalData->h_vn_pp, var1Data->h_vn_pp, var2Data->h_vn_pp, var3Data->h_vn_pp, var4Data->h_vn_pp, v_vn_pp);
  mergePoints(normalData->h_vn_pm, var1Data->h_vn_pm, var2Data->h_vn_pm, var3Data->h_vn_pm, var4Data->h_vn_pm, v_vn_pm);
  mergePoints(normalData->h_vn_kp, var1Data->h_vn_kp, var2Data->h_vn_kp, var3Data->h_vn_kp, var4Data->h_vn_kp, v_vn_kp);
  mergePoints(normalData->h_vn_km, var1Data->h_vn_km, var2Data->h_vn_km, var3Data->h_vn_km, var4Data->h_vn_km, v_vn_km);
  mergePoints(normalData->h_vn_pr, var1Data->h_vn_pr, var2Data->h_vn_pr, var3Data->h_vn_pr, var4Data->h_vn_pr, v_vn_pr);

  mergePoints(normalData->h_vn_yCM_00to10_pr, var1Data->h_vn_yCM_00to10_pr, var2Data->h_vn_yCM_00to10_pr, var3Data->h_vn_yCM_00to10_pr, var4Data->h_vn_yCM_00to10_pr, v_vn_yCM_00to10_pr);
  mergePoints(normalData->h_vn_yCM_10to40_pr, var1Data->h_vn_yCM_10to40_pr, var2Data->h_vn_yCM_10to40_pr, var3Data->h_vn_yCM_10to40_pr, var4Data->h_vn_yCM_10to40_pr, v_vn_yCM_10to40_pr);
  mergePoints(normalData->h_vn_yCM_40to60_pr, var1Data->h_vn_yCM_40to60_pr, var2Data->h_vn_yCM_40to60_pr, var3Data->h_vn_yCM_40to60_pr, var4Data->h_vn_yCM_40to60_pr, v_vn_yCM_40to60_pr);

  mergePoints(normalData->h_vn_yCM_HADES, var1Data->h_vn_yCM_HADES, var2Data->h_vn_yCM_HADES, var3Data->h_vn_yCM_HADES, var4Data->h_vn_yCM_HADES, v_vn_yCM_HADES);

  mergePoints(normalData->h_vn_yCM_00to10_pr_symm, var1Data->h_vn_yCM_00to10_pr_symm, var2Data->h_vn_yCM_00to10_pr_symm, var3Data->h_vn_yCM_00to10_pr_symm, var4Data->h_vn_yCM_00to10_pr_symm, v_vn_yCM_00to10_pr_symm);
  mergePoints(normalData->h_vn_yCM_10to40_pr_symm, var1Data->h_vn_yCM_10to40_pr_symm, var2Data->h_vn_yCM_10to40_pr_symm, var3Data->h_vn_yCM_10to40_pr_symm, var4Data->h_vn_yCM_10to40_pr_symm, v_vn_yCM_10to40_pr_symm);
  mergePoints(normalData->h_vn_yCM_40to60_pr_symm, var1Data->h_vn_yCM_40to60_pr_symm, var2Data->h_vn_yCM_40to60_pr_symm, var3Data->h_vn_yCM_40to60_pr_symm, var4Data->h_vn_yCM_40to60_pr_symm, v_vn_yCM_40to60_pr_symm);

  mergePoints(normalData->h_vn_yCM_00to05_pr_symm, var1Data->h_vn_yCM_00to05_pr_symm, var2Data->h_vn_yCM_00to05_pr_symm, var3Data->h_vn_yCM_00to05_pr_symm, var4Data->h_vn_yCM_00to05_pr_symm, v_vn_yCM_00to05_pr_symm);
  mergePoints(normalData->h_vn_yCM_05to10_pr_symm, var1Data->h_vn_yCM_05to10_pr_symm, var2Data->h_vn_yCM_05to10_pr_symm, var3Data->h_vn_yCM_05to10_pr_symm, var4Data->h_vn_yCM_05to10_pr_symm, v_vn_yCM_05to10_pr_symm);
  mergePoints(normalData->h_vn_yCM_10to15_pr_symm, var1Data->h_vn_yCM_10to15_pr_symm, var2Data->h_vn_yCM_10to15_pr_symm, var3Data->h_vn_yCM_10to15_pr_symm, var4Data->h_vn_yCM_10to15_pr_symm, v_vn_yCM_10to15_pr_symm);
  mergePoints(normalData->h_vn_yCM_15to20_pr_symm, var1Data->h_vn_yCM_15to20_pr_symm, var2Data->h_vn_yCM_15to20_pr_symm, var3Data->h_vn_yCM_15to20_pr_symm, var4Data->h_vn_yCM_15to20_pr_symm, v_vn_yCM_15to20_pr_symm);
  mergePoints(normalData->h_vn_yCM_20to25_pr_symm, var1Data->h_vn_yCM_20to25_pr_symm, var2Data->h_vn_yCM_20to25_pr_symm, var3Data->h_vn_yCM_20to25_pr_symm, var4Data->h_vn_yCM_20to25_pr_symm, v_vn_yCM_20to25_pr_symm);
  mergePoints(normalData->h_vn_yCM_25to30_pr_symm, var1Data->h_vn_yCM_25to30_pr_symm, var2Data->h_vn_yCM_25to30_pr_symm, var3Data->h_vn_yCM_25to30_pr_symm, var4Data->h_vn_yCM_25to30_pr_symm, v_vn_yCM_25to30_pr_symm);
  mergePoints(normalData->h_vn_yCM_30to35_pr_symm, var1Data->h_vn_yCM_30to35_pr_symm, var2Data->h_vn_yCM_30to35_pr_symm, var3Data->h_vn_yCM_30to35_pr_symm, var4Data->h_vn_yCM_30to35_pr_symm, v_vn_yCM_30to35_pr_symm);
  mergePoints(normalData->h_vn_yCM_35to40_pr_symm, var1Data->h_vn_yCM_35to40_pr_symm, var2Data->h_vn_yCM_35to40_pr_symm, var3Data->h_vn_yCM_35to40_pr_symm, var4Data->h_vn_yCM_35to40_pr_symm, v_vn_yCM_35to40_pr_symm);
  mergePoints(normalData->h_vn_yCM_40to45_pr_symm, var1Data->h_vn_yCM_40to45_pr_symm, var2Data->h_vn_yCM_40to45_pr_symm, var3Data->h_vn_yCM_40to45_pr_symm, var4Data->h_vn_yCM_40to45_pr_symm, v_vn_yCM_40to45_pr_symm);
  mergePoints(normalData->h_vn_yCM_45to50_pr_symm, var1Data->h_vn_yCM_45to50_pr_symm, var2Data->h_vn_yCM_45to50_pr_symm, var3Data->h_vn_yCM_45to50_pr_symm, var4Data->h_vn_yCM_45to50_pr_symm, v_vn_yCM_45to50_pr_symm);
  mergePoints(normalData->h_vn_yCM_50to55_pr_symm, var1Data->h_vn_yCM_50to55_pr_symm, var2Data->h_vn_yCM_50to55_pr_symm, var3Data->h_vn_yCM_50to55_pr_symm, var4Data->h_vn_yCM_50to55_pr_symm, v_vn_yCM_50to55_pr_symm);
  mergePoints(normalData->h_vn_yCM_55to60_pr_symm, var1Data->h_vn_yCM_55to60_pr_symm, var2Data->h_vn_yCM_55to60_pr_symm, var3Data->h_vn_yCM_55to60_pr_symm, var4Data->h_vn_yCM_55to60_pr_symm, v_vn_yCM_55to60_pr_symm);
  
  mergePoints(normalData->h_vn_pT_00to10_pr, var1Data->h_vn_pT_00to10_pr, var2Data->h_vn_pT_00to10_pr, var3Data->h_vn_pT_00to10_pr, var4Data->h_vn_pT_00to10_pr, v_vn_pT_00to10_pr);
  mergePoints(normalData->h_vn_pT_10to40_pr, var1Data->h_vn_pT_10to40_pr, var2Data->h_vn_pT_10to40_pr, var3Data->h_vn_pT_10to40_pr, var4Data->h_vn_pT_10to40_pr, v_vn_pT_10to40_pr);
  mergePoints(normalData->h_vn_pT_40to60_pr, var1Data->h_vn_pT_40to60_pr, var2Data->h_vn_pT_40to60_pr, var3Data->h_vn_pT_40to60_pr, var4Data->h_vn_pT_40to60_pr, v_vn_pT_40to60_pr);

  mergePoints(normalData->h_vn_pT_bin6_10to40_pr_symm,  var1Data->h_vn_pT_bin6_10to40_pr_symm,  var2Data->h_vn_pT_bin6_10to40_pr_symm,   var3Data->h_vn_pT_bin6_10to40_pr_symm, var4Data->h_vn_pT_bin6_10to40_pr_symm,  v_vn_pT_bin6_pr_symm);
  mergePoints(normalData->h_vn_pT_bin7_10to40_pr_symm,  var1Data->h_vn_pT_bin7_10to40_pr_symm,  var2Data->h_vn_pT_bin7_10to40_pr_symm,   var3Data->h_vn_pT_bin7_10to40_pr_symm, var4Data->h_vn_pT_bin7_10to40_pr_symm,  v_vn_pT_bin7_pr_symm);
  mergePoints(normalData->h_vn_pT_bin8_10to40_pr_symm,  var1Data->h_vn_pT_bin8_10to40_pr_symm,  var2Data->h_vn_pT_bin8_10to40_pr_symm,   var3Data->h_vn_pT_bin8_10to40_pr_symm, var4Data->h_vn_pT_bin8_10to40_pr_symm,  v_vn_pT_bin8_pr_symm);
  mergePoints(normalData->h_vn_pT_bin9_10to40_pr_symm,  var1Data->h_vn_pT_bin9_10to40_pr_symm,  var2Data->h_vn_pT_bin9_10to40_pr_symm,   var3Data->h_vn_pT_bin9_10to40_pr_symm, var4Data->h_vn_pT_bin9_10to40_pr_symm,  v_vn_pT_bin9_pr_symm);
  mergePoints(normalData->h_vn_pT_bin10_10to40_pr_symm, var1Data->h_vn_pT_bin10_10to40_pr_symm, var2Data->h_vn_pT_bin10_10to40_pr_symm,  var3Data->h_vn_pT_bin10_10to40_pr_symm, var4Data->h_vn_pT_bin10_10to40_pr_symm, v_vn_pT_bin10_pr_symm);
  mergePoints(normalData->h_vn_pT_bin11_10to40_pr_symm, var1Data->h_vn_pT_bin11_10to40_pr_symm, var2Data->h_vn_pT_bin11_10to40_pr_symm,  var3Data->h_vn_pT_bin11_10to40_pr_symm, var4Data->h_vn_pT_bin11_10to40_pr_symm, v_vn_pT_bin11_pr_symm);
  mergePoints(normalData->h_vn_pT_bin12_10to40_pr_symm, var1Data->h_vn_pT_bin12_10to40_pr_symm, var2Data->h_vn_pT_bin12_10to40_pr_symm,  var3Data->h_vn_pT_bin12_10to40_pr_symm, var4Data->h_vn_pT_bin12_10to40_pr_symm, v_vn_pT_bin12_pr_symm);
}

#include "Variation.h"
#include "CompositeData.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "THStack.h"
#include <vector>
#include <iostream>


void printPointErrors(TGraphErrors *graph);

struct AverageContributionTracker
{
  Double_t epdPercentQuadSum = 0.0;
  Double_t nhitsPercentQuadSum = 0.0;
  Double_t nSigPiPercentQuadSum = 0.0;
  Double_t nSigKaPercentQuadSum = 0.0;
  Double_t nSigPrPercentQuadSum = 0.0;
  Double_t rvtxPercentQuadSum = 0.0;
  Double_t zvtxPercentQuadSum = 0.0;
  Double_t dcaPercentQuadSum = 0.0;
  Double_t nhitsdEdxPercentQuadSum = 0.0;
  Double_t nhitsratioPercentQuadSum = 0.0;
  Double_t m2PiPercentQuadSum = 0.0;
  Double_t m2KaPercentQuadSum = 0.0;

  Double_t trackQAPercentQuadSum = 0.0;
  Double_t eventQAPercentQuadSum = 0.0;

  Int_t epdNbins = 0;
  Int_t nhitsNbins = 0;
  Int_t nSigPiNbins = 0;
  Int_t nSigKaNbins = 0;
  Int_t nSigPrNbins = 0;
  Int_t rvtxNbins = 0;
  Int_t zvtxNbins = 0;
  Int_t dcaNbins = 0;
  Int_t nhitsdEdxNbins = 0;
  Int_t nhitsratioNbins = 0;
  Int_t m2PiNbins = 0;
  Int_t m2KaNbins = 0;

  Int_t trackQANbins = 0;
  Int_t eventQANbins = 0;

  void addContribution(TString ID, Double_t stdDevContributed)
  {
    if (ID == "nhits" ||
	ID == "dca" ||
	ID == "nhitsdEdx" ||
	ID == "nhitsratio")
      {
	trackQAPercentQuadSum += stdDevContributed;
	trackQANbins++;
      }
    else if (ID == "rvtx" ||
	     ID == "zvtx")
      {
	eventQAPercentQuadSum += stdDevContributed;
	eventQANbins++;
      }
    else if (ID == "epd")
      {
	epdPercentQuadSum += stdDevContributed;
	epdNbins++;
      }
    else if (ID == "nSigPi")
      {
	nSigPiPercentQuadSum += stdDevContributed;
	nSigPiNbins++;
      }
    else if (ID == "nSigKa")
      {
	nSigKaPercentQuadSum += stdDevContributed;
	nSigKaNbins++;
      }
    else if (ID == "nSigPr")
      {
	nSigPrPercentQuadSum += stdDevContributed;
	nSigPrNbins++;
      }
    /*
    else if (ID == "rvtx")
      {
	rvtxPercentQuadSum += stdDevContributed;
	rvtxNbins++;
      }
    else if (ID == "zvtx")
      {
	zvtxPercentQuadSum += stdDevContributed;
	zvtxNbins++;
      }
    else if (ID == "dca")
      {
	dcaPercentQuadSum += stdDevContributed;
	dcaNbins++;
      }
    else if (ID == "nhits")
      {
	nhitsPercentQuadSum += stdDevContributed;
	nhitsNbins++;
      }
    else if (ID == "nhitsdEdx")
      {
	nhitsdEdxPercentQuadSum += stdDevContributed;
	nhitsdEdxNbins++;
      }
    else if (ID == "nhitsratio")
      {
	nhitsratioPercentQuadSum += stdDevContributed;
	nhitsratioNbins++;
      }
    */
    else if (ID == "m2Pi")
      {
	m2PiPercentQuadSum += stdDevContributed;
	m2PiNbins++;
      }
    else if (ID == "m2Ka")
      {
	m2KaPercentQuadSum += stdDevContributed;
	m2KaNbins++;
      }
    else
      {
	std::cout << "AverageContributionTracker error, check spelling of ID!" << std::endl;
      }
  }

  void printContributions()
  {
    std::cout << std::endl 
	      << "Std dev of each contribution as a percent of the normal v3 measurement:" << std::endl
	      << std::endl
	      << "Track QA, " << trackQAPercentQuadSum / (Double_t)trackQANbins << std::endl
      	      << "Event QA, " << eventQAPercentQuadSum / (Double_t)eventQANbins << std::endl
            /*
	      << "nHits, " << nhitsPercentQuadSum / (Double_t)nhitsNbins << std::endl
	      << "nHits dEdx, " << nhitsdEdxPercentQuadSum / (Double_t)nhitsdEdxNbins << std::endl
	      << "nHits Ratio, " << nhitsratioPercentQuadSum / (Double_t)nhitsratioNbins << std::endl
	      << "DCA, " << dcaPercentQuadSum / (Double_t)dcaNbins << std::endl
	      << "r Vertex, " << rvtxPercentQuadSum / (Double_t)rvtxNbins << std::endl
	      << "z Vertex, " << zvtxPercentQuadSum / (Double_t)zvtxNbins << std::endl
	    */
	      << "nSigma Pi, " << nSigPiPercentQuadSum / (Double_t)nSigPiNbins << std::endl
	      << "nSigma Ka, " << nSigKaPercentQuadSum / (Double_t)nSigKaNbins << std::endl
	      << "nSigma Pr, " << nSigPrPercentQuadSum / (Double_t)nSigPrNbins << std::endl
	      << "m^2 Pi, " << m2PiPercentQuadSum / (Double_t)m2PiNbins << std::endl
	      << "m^2 Ka, " << m2KaPercentQuadSum / (Double_t)m2KaNbins << std::endl
	      << "EP Resolution, " << epdPercentQuadSum / (Double_t)epdNbins << std::endl
	      << std::endl;
  }
};






void calculateSystematics(TString order_n_str = "3")
{
  //TFile* newFile = new TFile("systematicErrors_3p0GeV_withSharang_NEW.root", "RECREATE");
  TFile* newFile = new TFile("systematicErrors_3p0GeV.root", "RECREATE");

  //TFile* sharangFile = TFile::Open("v3_proton_comparison.root");
  //TProfile* p_3p0GeV_1040 = (TProfile*)sharangFile->Get("3p0_p_1040");
  
  TString directory = "../../../Results_July2023/flowResults_3p0GeV/thirdDraft/";
  TString directory20Percent = "../../../Results_July2023/flowResults_3p0GeV/thirdDraft/20percentVariations/";
  TString directory30Percent = "../../../Results_July2023/flowResults_3p0GeV/thirdDraft/30percentVariations/";
  TString directoryEPR = "../../../Results_July2023/flowResults_3p0GeV/thirdDraft/eprVariations/";
  TString energyPrefix = "";
  
  Variation* Normal = new Variation("Normal",directory+energyPrefix+"Normal_averagedRes", order_n_str);
  //Variation* Normal = new Variation("Normal",directory+energyPrefix+"Normal", order_n_str);
  Variation* epd_high = new Variation("epd_high",directoryEPR+energyPrefix+"epd_high", order_n_str);
  Variation* epd_scaled = new Variation("epd_low",directoryEPR+energyPrefix+"epd_low", order_n_str);

  Variation* nSigPi_high_20 = new Variation("nSigPi_high",directory20Percent+energyPrefix+"nSigPi_high", order_n_str);
  Variation* nSigPi_low_20  = new Variation("nSigPi_low",directory20Percent+energyPrefix+"nSigPi_low", order_n_str);
  Variation* nSigKa_high_20 = new Variation("nSigKa_high",directory20Percent+energyPrefix+"nSigKa_high", order_n_str);
  Variation* nSigKa_low_20  = new Variation("nSigKa_low",directory20Percent+energyPrefix+"nSigKa_low", order_n_str);
  Variation* nSigPr_high_20 = new Variation("nSigPr_high",directory20Percent+energyPrefix+"nSigPr_high", order_n_str);
  Variation* nSigPr_low_20  = new Variation("nSigPr_low",directory20Percent+energyPrefix+"nSigPr_low", order_n_str);
  //Variation* rvtx_high_20 = new Variation("rvtx_high",directory20Percent+energyPrefix+"rvtx_high", order_n_str);
  //Variation* rvtx_low_20  = new Variation("rvtx_low",directory20Percent+energyPrefix+"rvtx_low", order_n_str);
  //Variation* zvtx_high_20 = new Variation("zvtx_high",directory20Percent+energyPrefix+"zvtx_high", order_n_str);
  //Variation* zvtx_low_20  = new Variation("zvtx_low",directory20Percent+energyPrefix+"zvtx_low", order_n_str);
  Variation* dca_high_20 = new Variation("dca_high",directory20Percent+energyPrefix+"dca_high", order_n_str);
  Variation* dca_low_20  = new Variation("dca_low",directory20Percent+energyPrefix+"dca_low", order_n_str);
  Variation* nhits_high_20 = new Variation("nhits_high",directory20Percent+energyPrefix+"nhits_high", order_n_str);
  Variation* nhits_low_20 = new Variation("nhits_low",directory20Percent+energyPrefix+"nhits_low", order_n_str);
  Variation* nhitsdEdx_high_20 = new Variation("nhitsdEdx_high",directory20Percent+energyPrefix+"nhitsdEdx_high", order_n_str);
  //Variation* nhitsdEdx_low_20  = new Variation("nhitsdEdx_low",directory20Percent+energyPrefix+"nhitsdEdx_low", order_n_str);
  Variation* nhitsratio_high_20 = new Variation("nhitsratio_high",directory20Percent+energyPrefix+"nhitsratio_high", order_n_str);
  Variation* nhitsratio_low_20  = new Variation("nhitsratio_low",directory20Percent+energyPrefix+"nhitsratio_low", order_n_str);
  Variation* m2Pi_high_20 = new Variation("m2Pi_high",directory20Percent+energyPrefix+"m2Pi_high", order_n_str);
  Variation* m2Pi_low_20  = new Variation("m2Pi_low",directory20Percent+energyPrefix+"m2Pi_low", order_n_str);
  Variation* m2Ka_high_20 = new Variation("m2Ka_high",directory20Percent+energyPrefix+"m2Ka_high", order_n_str);
  Variation* m2Ka_low_20  = new Variation("m2Pi_low",directory20Percent+energyPrefix+"m2Ka_low", order_n_str);

  Variation* nSigPi_high_30 = new Variation("nSigPi_high",directory30Percent+energyPrefix+"nSigPi_high", order_n_str);
  Variation* nSigPi_low_30  = new Variation("nSigPi_low",directory30Percent+energyPrefix+"nSigPi_low", order_n_str);
  Variation* nSigKa_high_30 = new Variation("nSigKa_high",directory30Percent+energyPrefix+"nSigKa_high", order_n_str);
  Variation* nSigKa_low_30  = new Variation("nSigKa_low",directory30Percent+energyPrefix+"nSigKa_low", order_n_str);
  Variation* nSigPr_high_30 = new Variation("nSigPr_high",directory30Percent+energyPrefix+"nSigPr_high", order_n_str);
  Variation* nSigPr_low_30  = new Variation("nSigPr_low",directory30Percent+energyPrefix+"nSigPr_low", order_n_str);
  //Variation* rvtx_high_30 = new Variation("rvtx_high",directory30Percent+energyPrefix+"rvtx_high", order_n_str);
  //Variation* rvtx_low_30  = new Variation("rvtx_low",directory30Percent+energyPrefix+"rvtx_low", order_n_str);
  //Variation* zvtx_high_30 = new Variation("zvtx_high",directory30Percent+energyPrefix+"zvtx_high", order_n_str);
  //Variation* zvtx_low_30  = new Variation("zvtx_low",directory30Percent+energyPrefix+"zvtx_low", order_n_str);
  Variation* dca_high_30 = new Variation("dca_high",directory30Percent+energyPrefix+"dca_high", order_n_str);
  Variation* dca_low_30  = new Variation("dca_low",directory30Percent+energyPrefix+"dca_low", order_n_str);
  Variation* nhits_high_30 = new Variation("nhits_high",directory30Percent+energyPrefix+"nhits_high", order_n_str);
  Variation* nhits_low_30 = new Variation("nhits_low",directory30Percent+energyPrefix+"nhits_low", order_n_str);
  Variation* nhitsdEdx_high_30 = new Variation("nhitsdEdx_high",directory30Percent+energyPrefix+"nhitsdEdx_high", order_n_str);
  //Variation* nhitsdEdx_low_30  = new Variation("nhitsdEdx_low",directory30Percent+energyPrefix+"nhitsdEdx_low", order_n_str);
  Variation* nhitsratio_high_30 = new Variation("nhitsratio_high",directory30Percent+energyPrefix+"nhitsratio_high", order_n_str);
  Variation* nhitsratio_low_30  = new Variation("nhitsratio_low",directory30Percent+energyPrefix+"nhitsratio_low", order_n_str);
  Variation* m2Pi_high_30 = new Variation("m2Pi_high",directory30Percent+energyPrefix+"m2Pi_high", order_n_str);
  Variation* m2Pi_low_30  = new Variation("m2Pi_low",directory30Percent+energyPrefix+"m2Pi_low", order_n_str);
  Variation* m2Ka_high_30 = new Variation("m2Ka_high",directory30Percent+energyPrefix+"m2Ka_high", order_n_str);
  Variation* m2Ka_low_30  = new Variation("m2Ka_low",directory30Percent+energyPrefix+"m2Ka_low", order_n_str);


  CompositeData* epd = new CompositeData("epd", Normal, epd_scaled, epd_high);
  CompositeData* nhits = new CompositeData("nhits", Normal, nhits_low_30, nhits_high_30, nhits_low_20, nhits_high_20);  
  CompositeData* nSigPi = new CompositeData("nSigPi", Normal, nSigPi_low_30, nSigPi_high_30, nSigPi_low_20, nSigPi_high_20);
  CompositeData* nSigKa = new CompositeData("nSigKa", Normal, nSigKa_low_30, nSigKa_high_30, nSigKa_low_20, nSigKa_high_20);
  CompositeData* nSigPr = new CompositeData("nSigPr", Normal, nSigPr_low_30, nSigPr_high_30, nSigPr_low_20, nSigPr_high_20);
  //CompositeData* rvtx = new CompositeData("rvtx", Normal, rvtx_low_30, rvtx_high_30, rvtx_low_20, rvtx_high_20);
  //CompositeData* zvtx = new CompositeData("zvtx", Normal, zvtx_low_30, zvtx_high_30, zvtx_low_20, zvtx_high_20);
  CompositeData* dca  = new CompositeData("dca", Normal, dca_low_30, dca_high_30, dca_low_20, dca_high_20);
  //CompositeData* nhitsdEdx = new CompositeData("nhitsdEdx", Normal, nhitsdEdx_high_30, nhitsdEdx_high_20, nhitsdEdx_low_20);
  CompositeData* nhitsdEdx = new CompositeData("nhitsdEdx", Normal, nhitsdEdx_high_30, nhitsdEdx_high_20);
  CompositeData* nhitsratio = new CompositeData("nhitsratio", Normal, nhitsratio_low_30, nhitsratio_high_30, nhitsratio_low_20, nhitsratio_high_20);
  CompositeData* m2Pi = new CompositeData("m2Pi", Normal, m2Pi_low_30, m2Pi_high_30, m2Pi_low_20, m2Pi_high_20);
  CompositeData* m2Ka = new CompositeData("m2Ka", Normal, m2Ka_low_30, m2Ka_high_30, m2Ka_low_20, m2Ka_high_20);

  
  // Any variations applied universally (like epd variation) should not be in this vector.
  std::vector<CompositeData*> composites;
  composites.push_back(nhits);
  composites.push_back(nSigPi);
  composites.push_back(nSigKa);
  composites.push_back(nSigPr);
  //composites.push_back(rvtx);
  //composites.push_back(zvtx);
  composites.push_back(dca);
  composites.push_back(nhitsdEdx);
  composites.push_back(nhitsratio);
  composites.push_back(m2Pi);
  composites.push_back(m2Ka);
  ////
  
  newFile->cd();

  //======= CALCULATION OF SYSTEMATIC ERRORS
  Double_t ithBinSysErr = 0.0;
  Double_t quadSum = 0.0;
  Int_t bins;

  AverageContributionTracker avgTracker_00to10;
  AverageContributionTracker avgTracker_10to40;
  AverageContributionTracker avgTracker_40to60;

  //=== pi+ vs centrality
  //std::cout << "PP VS CENTRALITY" << std::endl;
  std::vector<Double_t> v_sys_pp;
  bins = Normal->h_vn_pp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
    {
      quadSum = 0.0;
      quadSum += epd->v_vn_pp.at(ithBin).variance;

      if (ithBin < 2)
	avgTracker_00to10.addContribution(epd->ID, epd->v_vn_pp.at(ithBin).stdDevPercentage);
      else if (ithBin >= 2 && ithBin <= 7)
	avgTracker_10to40.addContribution(epd->ID, epd->v_vn_pp.at(ithBin).stdDevPercentage);
      else if (ithBin > 7)
	avgTracker_40to60.addContribution(epd->ID, epd->v_vn_pp.at(ithBin).stdDevPercentage);
 
      for (int jthCut = 0; jthCut < composites.size(); jthCut++)
	{
	  if (composites.at(jthCut)->v_vn_pp.at(ithBin).deltaByDeltaError > 1.0)
	    {
	      quadSum += composites.at(jthCut)->v_vn_pp.at(ithBin).variance;
	      if (ithBin < 2)
		avgTracker_00to10.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pp.at(ithBin).stdDevPercentage);
	      else if (ithBin >= 2 && ithBin <= 7)
		avgTracker_10to40.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pp.at(ithBin).stdDevPercentage);
	      else if (ithBin > 7)
		avgTracker_40to60.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pp.at(ithBin).stdDevPercentage);
	    }
	}
    
      ithBinSysErr = TMath::Sqrt(quadSum);
      v_sys_pp.push_back(ithBinSysErr);
    }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_pp loop


  
  //=== pi- vs centrality
  //std::cout << "PM VS CENTRALITY" << std::endl;
  std::vector<Double_t> v_sys_pm;
  bins = Normal->h_vn_pm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pm.at(ithBin).variance;

    if (ithBin < 2)
      avgTracker_00to10.addContribution(epd->ID, epd->v_vn_pm.at(ithBin).stdDevPercentage);
    if (ithBin >= 2 && ithBin <= 7)
      avgTracker_10to40.addContribution(epd->ID, epd->v_vn_pm.at(ithBin).stdDevPercentage);
    if (ithBin > 7)
      avgTracker_40to60.addContribution(epd->ID, epd->v_vn_pm.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_pm.at(ithBin).deltaByDeltaError > 1.0)
        {
          quadSum += composites.at(jthCut)->v_vn_pm.at(ithBin).variance;

	  if (ithBin < 2)
	    avgTracker_00to10.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pm.at(ithBin).stdDevPercentage);
	  if (ithBin >= 2 && ithBin <= 7)
	    avgTracker_10to40.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pm.at(ithBin).stdDevPercentage);
	  if (ithBin > 7)
	    avgTracker_40to60.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pm.at(ithBin).stdDevPercentage);

        }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_pm loop


  
  //=== K+ vs centrality
  //std::cout << "KP VS CENTRALITY" << std::endl;
  std::vector<Double_t> v_sys_kp;
  bins = Normal->h_vn_kp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_kp.at(ithBin).variance;

    if (ithBin < 2)
      avgTracker_00to10.addContribution(epd->ID, epd->v_vn_kp.at(ithBin).stdDevPercentage);
    if (ithBin >= 2 && ithBin <= 7)
      avgTracker_10to40.addContribution(epd->ID, epd->v_vn_kp.at(ithBin).stdDevPercentage);
    if (ithBin > 7)
      avgTracker_40to60.addContribution(epd->ID, epd->v_vn_kp.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_kp.at(ithBin).deltaByDeltaError > 1.0)
	{
	  quadSum += composites.at(jthCut)->v_vn_kp.at(ithBin).variance;

	  if (ithBin < 2)
	    avgTracker_00to10.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_kp.at(ithBin).stdDevPercentage);
	  if (ithBin >= 2 && ithBin <= 7)
	    avgTracker_10to40.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_kp.at(ithBin).stdDevPercentage);
	  if (ithBin > 7)
	    avgTracker_40to60.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_kp.at(ithBin).stdDevPercentage);

	}
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_kp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_kp loop

  //=== K- vs centrality
  //std::cout << "KM VS CENTRALITY" << std::endl;
  std::vector<Double_t> v_sys_km;
  bins = Normal->h_vn_km->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_km.at(ithBin).variance;

    if (ithBin < 2)
      avgTracker_00to10.addContribution(epd->ID, epd->v_vn_km.at(ithBin).stdDevPercentage);
    if (ithBin >= 2 && ithBin <= 7)
      avgTracker_10to40.addContribution(epd->ID, epd->v_vn_km.at(ithBin).stdDevPercentage);
    if (ithBin > 7)
      avgTracker_40to60.addContribution(epd->ID, epd->v_vn_km.at(ithBin).stdDevPercentage);

 
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_km.at(ithBin).deltaByDeltaError > 1.0)
	{
	  quadSum += composites.at(jthCut)->v_vn_km.at(ithBin).variance;

	  if (ithBin < 2)
	    avgTracker_00to10.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_km.at(ithBin).stdDevPercentage);
	  if (ithBin >= 2 && ithBin <= 7)
	    avgTracker_10to40.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_km.at(ithBin).stdDevPercentage);
	  if (ithBin > 7)
	    avgTracker_40to60.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_km.at(ithBin).stdDevPercentage);

	}
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_km.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_km loop


  //=== proton vs centrality
  //std::cout << "PROTON VS CENTRALITY" << std::endl;
  std::vector<Double_t> v_sys_pr;
  bins = Normal->h_vn_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pr.at(ithBin).variance;

    if (ithBin < 2)
      avgTracker_00to10.addContribution(epd->ID, epd->v_vn_pr.at(ithBin).stdDevPercentage);
    if (ithBin >= 2 && ithBin <= 7)
      avgTracker_10to40.addContribution(epd->ID, epd->v_vn_pr.at(ithBin).stdDevPercentage);
    if (ithBin > 7)
      avgTracker_40to60.addContribution(epd->ID, epd->v_vn_pr.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_pr.at(ithBin).deltaByDeltaError > 1.0)
        {
          quadSum += composites.at(jthCut)->v_vn_pr.at(ithBin).variance;

	  if (ithBin < 2)
	    avgTracker_00to10.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pr.at(ithBin).stdDevPercentage);
	  if (ithBin >= 2 && ithBin <= 7)
	    avgTracker_10to40.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pr.at(ithBin).stdDevPercentage);
	  if (ithBin > 7)
	    avgTracker_40to60.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pr.at(ithBin).stdDevPercentage);

        }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_pr loop


  //=== Proton vs rapidity 0 - 10%
  //std::cout << "PROTON VS RAPIDITY" << std::endl;
  std::vector<Double_t> v_sys_yCM_00to10_pr;
  bins = Normal->h_vn_yCM_00to10_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_00to10_pr.at(ithBin).variance;

    if (ithBin > 9)
      avgTracker_00to10.addContribution(epd->ID, epd->v_vn_yCM_00to10_pr.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_00to10_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthCut)->v_vn_yCM_00to10_pr.at(ithBin).variance;
        avgTracker_00to10.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_yCM_00to10_pr.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_00to10_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===



  //=== Proton vs rapidity 10 - 40%
  //std::cout << "PROTON VS RAPIDITY" << std::endl;
  std::vector<Double_t> v_sys_yCM_10to40_pr;
  bins = Normal->h_vn_yCM_10to40_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_10to40_pr.at(ithBin).variance;
    /*
    // Add difference with Sharang's points in systematics
    Double_t v3_diff = 0.0;
    if (ithBin > 9)
      {
        v3_diff = TMath::Abs( p_3p9GeV_1040->GetBinContent(ithBin+5) - Normal->h_vn_yCM_10to40_pr->GetBinContent(ithBin+1) ) / TMath::Sqrt(12);
        std::cout << v3_diff << std::endl;
        //quadSum += v3_diff*v3_diff;
      }
    */
    if (ithBin > 9)
      avgTracker_10to40.addContribution(epd->ID, epd->v_vn_yCM_10to40_pr.at(ithBin).stdDevPercentage);
 
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_10to40_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthCut)->v_vn_yCM_10to40_pr.at(ithBin).variance;
        avgTracker_10to40.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_yCM_10to40_pr.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    //ithBinSysErr += v3_diff;
    v_sys_yCM_10to40_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== Proton vs rapidity 40 - 60%
  //std::cout << "PROTON VS RAPIDITY" << std::endl;
  std::vector<Double_t> v_sys_yCM_40to60_pr;
  bins = Normal->h_vn_yCM_40to60_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_40to60_pr.at(ithBin).variance;
    
    if (ithBin > 9)
      avgTracker_40to60.addContribution(epd->ID, epd->v_vn_yCM_40to60_pr.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_40to60_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthCut)->v_vn_yCM_40to60_pr.at(ithBin).variance;
        avgTracker_40to60.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_yCM_40to60_pr.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_40to60_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== Proton vs rapidity 0 - 10% symmetric
  std::vector<Double_t> v_sys_yCM_00to10_pr_symm;
  bins = Normal->h_vn_yCM_00to10_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_00to10_pr_symm.at(ithBin).variance;
    
    if (ithBin > 4 && ithBin < 15)
      avgTracker_00to10.addContribution(epd->ID, epd->v_vn_yCM_00to10_pr_symm.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_00to10_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthCut)->v_vn_yCM_00to10_pr_symm.at(ithBin).variance;
        avgTracker_00to10.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_yCM_00to10_pr_symm.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_00to10_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== Proton vs rapidity 10 - 40% symmetric
  std::vector<Double_t> v_sys_yCM_10to40_pr_symm;
  bins = Normal->h_vn_yCM_10to40_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_10to40_pr_symm.at(ithBin).variance;
    
    if (ithBin > 4 && ithBin < 15)
      avgTracker_10to40.addContribution(epd->ID, epd->v_vn_yCM_10to40_pr_symm.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_10to40_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthCut)->v_vn_yCM_10to40_pr_symm.at(ithBin).variance;
        avgTracker_10to40.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_yCM_10to40_pr_symm.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_10to40_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== Proton vs rapidity 40 - 60% symmetric
  std::vector<Double_t> v_sys_yCM_40to60_pr_symm;
  bins = Normal->h_vn_yCM_40to60_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_40to60_pr_symm.at(ithBin).variance;
    
    if (ithBin > 4 && ithBin < 15)
      avgTracker_40to60.addContribution(epd->ID, epd->v_vn_yCM_40to60_pr_symm.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_40to60_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthCut)->v_vn_yCM_40to60_pr_symm.at(ithBin).variance;
        avgTracker_40to60.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_yCM_40to60_pr_symm.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_40to60_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===

  //=== Proton vs rapidity symmetric 5% bins
  std::vector<Double_t> v_sys_yCM_00to05_pr_symm;
  bins = Normal->h_vn_yCM_00to05_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_00to05_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_00to05_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_00to05_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_00to05_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_05to10_pr_symm;
  bins = Normal->h_vn_yCM_05to10_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_05to10_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_05to10_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_05to10_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_05to10_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_10to15_pr_symm;
  bins = Normal->h_vn_yCM_10to15_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_10to15_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_10to15_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_10to15_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_10to15_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_15to20_pr_symm;
  bins = Normal->h_vn_yCM_15to20_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_15to20_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_15to20_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_15to20_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_15to20_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_20to25_pr_symm;
  bins = Normal->h_vn_yCM_20to25_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_20to25_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_20to25_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_20to25_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_20to25_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_25to30_pr_symm;
  bins = Normal->h_vn_yCM_25to30_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_25to30_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_25to30_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_25to30_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_25to30_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_30to35_pr_symm;
  bins = Normal->h_vn_yCM_30to35_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_30to35_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_30to35_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_30to35_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_30to35_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_35to40_pr_symm;
  bins = Normal->h_vn_yCM_35to40_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_35to40_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_35to40_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_35to40_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_35to40_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_40to45_pr_symm;
  bins = Normal->h_vn_yCM_40to45_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_40to45_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_40to45_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_40to45_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_40to45_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_45to50_pr_symm;
  bins = Normal->h_vn_yCM_45to50_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_45to50_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_45to50_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_45to50_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_45to50_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_50to55_pr_symm;
  bins = Normal->h_vn_yCM_50to55_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_50to55_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_50to55_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_50to55_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_50to55_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_55to60_pr_symm;
  bins = Normal->h_vn_yCM_55to60_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_55to60_pr_symm.at(ithBin).variance;
    
    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_yCM_55to60_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      { quadSum += composites.at(jthCut)->v_vn_yCM_55to60_pr_symm.at(ithBin).variance; }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_55to60_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===

  //=== Proton vs pT 0 - 10%
  std::vector<Double_t> v_sys_pT_00to10_pr;
  bins = Normal->h_vn_pT_00to10_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_00to10_pr.at(ithBin).variance;
    
    if (ithBin > 1)
      avgTracker_00to10.addContribution(epd->ID, epd->v_vn_pT_00to10_pr.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_pT_00to10_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthCut)->v_vn_pT_00to10_pr.at(ithBin).variance;
        avgTracker_00to10.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pT_00to10_pr.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_00to10_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== Proton vs pT 10 - 40%
  std::vector<Double_t> v_sys_pT_10to40_pr;
  bins = Normal->h_vn_pT_10to40_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_10to40_pr.at(ithBin).variance;
    
    if (ithBin > 1)
      avgTracker_10to40.addContribution(epd->ID, epd->v_vn_pT_10to40_pr.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_pT_10to40_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthCut)->v_vn_pT_10to40_pr.at(ithBin).variance;
        avgTracker_10to40.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pT_10to40_pr.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_10to40_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== Proton vs pT 40 - 60%
  std::vector<Double_t> v_sys_pT_40to60_pr;
  bins = Normal->h_vn_pT_40to60_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_40to60_pr.at(ithBin).variance;
    
    if (ithBin > 1)
      avgTracker_40to60.addContribution(epd->ID, epd->v_vn_pT_40to60_pr.at(ithBin).stdDevPercentage);

    for (int jthCut = 0; jthCut < composites.size(); jthCut++)
    {
      if (composites.at(jthCut)->v_vn_pT_40to60_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthCut)->v_vn_pT_40to60_pr.at(ithBin).variance;
        avgTracker_40to60.addContribution(composites.at(jthCut)->ID, composites.at(jthCut)->v_vn_pT_40to60_pr.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_40to60_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===



  //std::cout << "0-10% Centrality" << std::endl;
  //avgTracker_00to10.printContributions();  
  
  std::cout << "10-40% Centrality" << std::endl;
  avgTracker_10to40.printContributions();

  //std::cout << "40-60% Centrality" << std::endl;
  //avgTracker_40to60.printContributions();

  newFile->cd();

  

  // SAVE V3 VALUES, SYSTEMATIC UNCERTAINTIES, AND # OF ENTRIES IN VECTOR FORM FOR ODD/EVEN PROPAGATION
  std::vector<Double_t> vn_yCM_10to40_pr_values;
  std::vector<Double_t> vn_yCM_10to40_pr_uncertainties;
  std::vector<Double_t> vn_yCM_10to40_pr_entries;

  std::vector<Double_t> vn_yCM_10to40_pr_symm_values;
  std::vector<Double_t> vn_yCM_10to40_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_10to40_pr_symm_entries;
  
  bins = Normal->h_vn_yCM_10to40_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
    {
      vn_yCM_10to40_pr_values.push_back(Normal->p_vn_yCM_10to40_pr->GetBinContent(ithBin+1));
      vn_yCM_10to40_pr_uncertainties.push_back(v_sys_yCM_10to40_pr.at(ithBin));
      vn_yCM_10to40_pr_entries.push_back(Normal->p_vn_yCM_10to40_pr->GetBinEffectiveEntries(ithBin+1));
    }

  bins = Normal->h_vn_yCM_10to40_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
    {
      vn_yCM_10to40_pr_symm_values.push_back(Normal->p_vn_yCM_10to40_pr_symm->GetBinContent(ithBin+1));
      vn_yCM_10to40_pr_symm_uncertainties.push_back(v_sys_yCM_10to40_pr_symm.at(ithBin));
      vn_yCM_10to40_pr_symm_entries.push_back(Normal->p_vn_yCM_10to40_pr_symm->GetBinEffectiveEntries(ithBin+1));
    }

  newFile->WriteObject(&vn_yCM_10to40_pr_values, "vn_yCM_10to40_pr_values");
  newFile->WriteObject(&vn_yCM_10to40_pr_uncertainties, "vn_yCM_10to40_pr_uncertainties");
  newFile->WriteObject(&vn_yCM_10to40_pr_entries, "vn_yCM_10to40_pr_entries");

  newFile->WriteObject(&vn_yCM_10to40_pr_symm_values, "vn_yCM_10to40_pr_symm_values");
  newFile->WriteObject(&vn_yCM_10to40_pr_symm_uncertainties, "vn_yCM_10to40_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_10to40_pr_symm_entries, "vn_yCM_10to40_pr_symm_entries");


  // 5% centrality bins here
  std::vector<Double_t> vn_yCM_00to05_pr_symm_values;
  std::vector<Double_t> vn_yCM_00to05_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_00to05_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_05to10_pr_symm_values;
  std::vector<Double_t> vn_yCM_05to10_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_05to10_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_10to15_pr_symm_values;
  std::vector<Double_t> vn_yCM_10to15_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_10to15_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_15to20_pr_symm_values;
  std::vector<Double_t> vn_yCM_15to20_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_15to20_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_20to25_pr_symm_values;
  std::vector<Double_t> vn_yCM_20to25_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_20to25_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_25to30_pr_symm_values;
  std::vector<Double_t> vn_yCM_25to30_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_25to30_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_30to35_pr_symm_values;
  std::vector<Double_t> vn_yCM_30to35_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_30to35_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_35to40_pr_symm_values;
  std::vector<Double_t> vn_yCM_35to40_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_35to40_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_40to45_pr_symm_values;
  std::vector<Double_t> vn_yCM_40to45_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_40to45_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_45to50_pr_symm_values;
  std::vector<Double_t> vn_yCM_45to50_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_45to50_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_50to55_pr_symm_values;
  std::vector<Double_t> vn_yCM_50to55_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_50to55_pr_symm_entries;
  
  std::vector<Double_t> vn_yCM_55to60_pr_symm_values;
  std::vector<Double_t> vn_yCM_55to60_pr_symm_uncertainties;
  std::vector<Double_t> vn_yCM_55to60_pr_symm_entries;
  
  bins = Normal->h_vn_yCM_00to05_pr_symm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    vn_yCM_00to05_pr_symm_values.push_back(Normal->p_vn_yCM_00to05_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_00to05_pr_symm_uncertainties.push_back(v_sys_yCM_00to05_pr_symm.at(ithBin));
    vn_yCM_00to05_pr_symm_entries.push_back(Normal->p_vn_yCM_00to05_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_05to10_pr_symm_values.push_back(Normal->p_vn_yCM_05to10_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_05to10_pr_symm_uncertainties.push_back(v_sys_yCM_05to10_pr_symm.at(ithBin));
    vn_yCM_05to10_pr_symm_entries.push_back(Normal->p_vn_yCM_05to10_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_10to15_pr_symm_values.push_back(Normal->p_vn_yCM_10to15_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_10to15_pr_symm_uncertainties.push_back(v_sys_yCM_10to15_pr_symm.at(ithBin));
    vn_yCM_10to15_pr_symm_entries.push_back(Normal->p_vn_yCM_10to15_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_15to20_pr_symm_values.push_back(Normal->p_vn_yCM_15to20_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_15to20_pr_symm_uncertainties.push_back(v_sys_yCM_15to20_pr_symm.at(ithBin));
    vn_yCM_15to20_pr_symm_entries.push_back(Normal->p_vn_yCM_15to20_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_20to25_pr_symm_values.push_back(Normal->p_vn_yCM_20to25_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_20to25_pr_symm_uncertainties.push_back(v_sys_yCM_20to25_pr_symm.at(ithBin));
    vn_yCM_20to25_pr_symm_entries.push_back(Normal->p_vn_yCM_20to25_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_25to30_pr_symm_values.push_back(Normal->p_vn_yCM_25to30_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_25to30_pr_symm_uncertainties.push_back(v_sys_yCM_25to30_pr_symm.at(ithBin));
    vn_yCM_25to30_pr_symm_entries.push_back(Normal->p_vn_yCM_25to30_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_30to35_pr_symm_values.push_back(Normal->p_vn_yCM_30to35_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_30to35_pr_symm_uncertainties.push_back(v_sys_yCM_30to35_pr_symm.at(ithBin));
    vn_yCM_30to35_pr_symm_entries.push_back(Normal->p_vn_yCM_30to35_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_35to40_pr_symm_values.push_back(Normal->p_vn_yCM_35to40_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_35to40_pr_symm_uncertainties.push_back(v_sys_yCM_35to40_pr_symm.at(ithBin));
    vn_yCM_35to40_pr_symm_entries.push_back(Normal->p_vn_yCM_35to40_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_40to45_pr_symm_values.push_back(Normal->p_vn_yCM_40to45_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_40to45_pr_symm_uncertainties.push_back(v_sys_yCM_40to45_pr_symm.at(ithBin));
    vn_yCM_40to45_pr_symm_entries.push_back(Normal->p_vn_yCM_40to45_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_45to50_pr_symm_values.push_back(Normal->p_vn_yCM_45to50_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_45to50_pr_symm_uncertainties.push_back(v_sys_yCM_45to50_pr_symm.at(ithBin));
    vn_yCM_45to50_pr_symm_entries.push_back(Normal->p_vn_yCM_45to50_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_50to55_pr_symm_values.push_back(Normal->p_vn_yCM_50to55_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_50to55_pr_symm_uncertainties.push_back(v_sys_yCM_50to55_pr_symm.at(ithBin));
    vn_yCM_50to55_pr_symm_entries.push_back(Normal->p_vn_yCM_50to55_pr_symm->GetBinEffectiveEntries(ithBin+1));

    vn_yCM_55to60_pr_symm_values.push_back(Normal->p_vn_yCM_55to60_pr_symm->GetBinContent(ithBin+1));
    vn_yCM_55to60_pr_symm_uncertainties.push_back(v_sys_yCM_55to60_pr_symm.at(ithBin));
    vn_yCM_55to60_pr_symm_entries.push_back(Normal->p_vn_yCM_55to60_pr_symm->GetBinEffectiveEntries(ithBin+1));
  }

  newFile->WriteObject(&vn_yCM_00to05_pr_symm_values, "vn_yCM_00to05_pr_symm_values");
  newFile->WriteObject(&vn_yCM_00to05_pr_symm_uncertainties, "vn_yCM_00to05_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_00to05_pr_symm_entries, "vn_yCM_00to05_pr_symm_entries");
  
  newFile->WriteObject(&vn_yCM_05to10_pr_symm_values, "vn_yCM_05to10_pr_symm_values");
  newFile->WriteObject(&vn_yCM_05to10_pr_symm_uncertainties, "vn_yCM_05to10_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_05to10_pr_symm_entries, "vn_yCM_05to10_pr_symm_entries");

  newFile->WriteObject(&vn_yCM_10to15_pr_symm_values, "vn_yCM_10to15_pr_symm_values");
  newFile->WriteObject(&vn_yCM_10to15_pr_symm_uncertainties, "vn_yCM_10to15_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_10to15_pr_symm_entries, "vn_yCM_10to15_pr_symm_entries");
  
  newFile->WriteObject(&vn_yCM_15to20_pr_symm_values, "vn_yCM_15to20_pr_symm_values");
  newFile->WriteObject(&vn_yCM_15to20_pr_symm_uncertainties, "vn_yCM_15to20_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_15to20_pr_symm_entries, "vn_yCM_15to20_pr_symm_entries");
  
  newFile->WriteObject(&vn_yCM_20to25_pr_symm_values, "vn_yCM_20to25_pr_symm_values");
  newFile->WriteObject(&vn_yCM_20to25_pr_symm_uncertainties, "vn_yCM_20to25_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_20to25_pr_symm_entries, "vn_yCM_20to25_pr_symm_entries");
  
  newFile->WriteObject(&vn_yCM_25to30_pr_symm_values, "vn_yCM_25to30_pr_symm_values");
  newFile->WriteObject(&vn_yCM_25to30_pr_symm_uncertainties, "vn_yCM_25to30_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_25to30_pr_symm_entries, "vn_yCM_25to30_pr_symm_entries");
  
  newFile->WriteObject(&vn_yCM_30to35_pr_symm_values, "vn_yCM_30to35_pr_symm_values");
  newFile->WriteObject(&vn_yCM_30to35_pr_symm_uncertainties, "vn_yCM_30to35_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_30to35_pr_symm_entries, "vn_yCM_30to35_pr_symm_entries");
  
  newFile->WriteObject(&vn_yCM_35to40_pr_symm_values, "vn_yCM_35to40_pr_symm_values");
  newFile->WriteObject(&vn_yCM_35to40_pr_symm_uncertainties, "vn_yCM_35to40_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_35to40_pr_symm_entries, "vn_yCM_35to40_pr_symm_entries");
  
  newFile->WriteObject(&vn_yCM_40to45_pr_symm_values, "vn_yCM_40to45_pr_symm_values");
  newFile->WriteObject(&vn_yCM_40to45_pr_symm_uncertainties, "vn_yCM_40to45_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_40to45_pr_symm_entries, "vn_yCM_40to45_pr_symm_entries");
  
  newFile->WriteObject(&vn_yCM_45to50_pr_symm_values, "vn_yCM_45to50_pr_symm_values");
  newFile->WriteObject(&vn_yCM_45to50_pr_symm_uncertainties, "vn_yCM_45to50_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_45to50_pr_symm_entries, "vn_yCM_45to50_pr_symm_entries");
  
  newFile->WriteObject(&vn_yCM_50to55_pr_symm_values, "vn_yCM_50to55_pr_symm_values");
  newFile->WriteObject(&vn_yCM_50to55_pr_symm_uncertainties, "vn_yCM_50to55_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_50to55_pr_symm_entries, "vn_yCM_50to55_pr_symm_entries");
  
  newFile->WriteObject(&vn_yCM_55to60_pr_symm_values, "vn_yCM_55to60_pr_symm_values");
  newFile->WriteObject(&vn_yCM_55to60_pr_symm_uncertainties, "vn_yCM_55to60_pr_symm_uncertainties");
  newFile->WriteObject(&vn_yCM_55to60_pr_symm_entries, "vn_yCM_55to60_pr_symm_entries");
  ////====



  // PLOTTING
  TCanvas *canvas = new TCanvas("canvas", "Canvas", 1200, 1200);
  //canvas->SetGrid();
  canvas->SetTicks();
  canvas->SetLeftMargin(0.15);
  canvas->cd();

  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(6);

  TLine *zeroLine = new TLine(0, 0, 60, 0);
  zeroLine->SetLineStyle(9);

  TLine *zeroLine_y = new TLine(0, 0, 1, 0);
  zeroLine_y->SetLineStyle(9);

  TLine *zeroLine_y_pr = new TLine(-1, 0, 1, 0);
  zeroLine_y_pr->SetLineStyle(9);

  TLine *zeroLine_pt = new TLine(0, 0, 2, 0);
  zeroLine_pt->SetLineStyle(9);


  TLegend *piLegend = new TLegend(0.67, 0.71, 0.805, 0.87);
  piLegend->AddEntry(Normal->h_vn_pp,"#pi^{+}");
  piLegend->AddEntry(Normal->h_vn_pm,"#pi^{-}");
  piLegend->SetFillColorAlpha(0,0);
  piLegend->SetLineColorAlpha(0,0);

  TLegend *kaLegend = new TLegend(0.7, 0.72, 0.8, 0.87);
  kaLegend->AddEntry(Normal->h_vn_kp,"K^{+}");
  kaLegend->AddEntry(Normal->h_vn_km,"K^{-}");
  kaLegend->SetFillColorAlpha(0,0);
  kaLegend->SetLineColorAlpha(0,0);      

  TLegend *prLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
  prLegend->AddEntry(Normal->h_vn_yCM_00to10_pr, "0 - 10%");
  prLegend->AddEntry(Normal->h_vn_yCM_10to40_pr, "10 - 40%");
  prLegend->AddEntry(Normal->h_vn_yCM_40to60_pr, "40 - 60%");
  prLegend->SetBorderSize(0);
  prLegend->SetFillColorAlpha(0,0);

      
  TLegend *prPtLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
  prPtLegend->AddEntry(Normal->h_vn_pT_00to10_pr, "0 - 10%");
  prPtLegend->AddEntry(Normal->h_vn_pT_10to40_pr, "10 - 40%");
  prPtLegend->AddEntry(Normal->h_vn_pT_40to60_pr, "40 - 60%");
  prPtLegend->SetBorderSize(0);
  prPtLegend->SetFillColorAlpha(0,0);


      
  TPaveText *piText = new TPaveText(5, 0.025, 38, 0.07, "NB");
  piText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  piText->AddText("0 < y_{CM} < 0.5 GeV");
  piText->AddText("0.18 < p_{T} < 1.6 GeV");
  piText->SetFillColorAlpha(0,0);
  piText->SetLineColorAlpha(0,0);

  TPaveText *kaText = new TPaveText(5, 0.025, 38, 0.07, "NB");
  kaText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  kaText->AddText("0 < y_{CM} < 0.5 GeV");
  kaText->AddText("0.4 < p_{T} < 1.6 GeV");
  kaText->SetFillColorAlpha(0,0);
  kaText->SetLineColorAlpha(0,0);
      
  TPaveText *prText = new TPaveText(5, 0.025, 38, 0.07, "NB");
  prText->AddText("Proton");
  prText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  prText->AddText("0 < y_{CM} < 0.5 GeV");
  prText->AddText("0.4 < p_{T} < 2.0 GeV");
  prText->SetFillColorAlpha(0,0);
  prText->SetLineColorAlpha(0,0);
      
  TPaveText *prText_y = new TPaveText(-0.2, 0.02, 0.9, 0.05, "NB");
  prText_y->AddText("Proton");
  prText_y->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  prText_y->AddText("0.4 < p_{T} < 2.0 GeV");
  prText_y->SetFillColorAlpha(0,0);
  prText_y->SetLineColorAlpha(0,0);
  prText_y->SetTextSize(.035);
      
  TPaveText *prText_y_symm = new TPaveText(-0.2, 0.02, 0.9, 0.05, "NB");
  prText_y_symm->AddText("Proton");
  prText_y_symm->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  prText_y_symm->AddText("1.0 < p_{T} < 2.5 GeV");
  prText_y_symm->SetFillColorAlpha(0,0);
  prText_y_symm->SetLineColorAlpha(0,0);
  prText_y_symm->SetTextSize(.035);
      
  TPaveText *prPtText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
  prPtText->AddText("Proton");
  prPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
  prPtText->AddText("0 < y_{CM} < 0.5");
  prPtText->SetFillColorAlpha(0,0);
  prPtText->SetLineColorAlpha(0,0);
  prPtText->SetTextSize(.04);

  
  Double_t centralityUpperBound = 0.08;
  Double_t centralityLowerBound = -0.08;
  Double_t rapidityUpperBound = 0.15;
  Double_t rapidityLowerBound = -0.1;
  Double_t rapidityUpperBound_pr = 0.06;
  Double_t rapidityLowerBound_pr = -0.1;
  Double_t ptUpperBound = 0.25;
  Double_t ptLowerBound = -0.25;
      
  TGraphErrors *copyWithNewErrors1;
  TGraphErrors *copyWithNewErrors2;
  TGraphErrors *copyWithNewErrors3;

  Double_t xAxisError = 0.0;
  
  //=== pi+- vs centrality
  copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pp->Clone());
  xAxisError = Normal->h_vn_pp->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_pp.size(); i++)
    { copyWithNewErrors1->SetPointError(i, xAxisError, v_sys_pp.at(i)); }

  copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pm->Clone());
  xAxisError = Normal->h_vn_pm->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_pm.size(); i++)
    { copyWithNewErrors2->SetPointError(i, xAxisError, v_sys_pm.at(i)); }

  copyWithNewErrors1->Write();
  copyWithNewErrors2->Write();

  THStack *piCentralityStack = new THStack("piCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  piCentralityStack->Add(Normal->h_vn_pp);
  piCentralityStack->Add(Normal->h_vn_pm);

  piCentralityStack->Draw();
  piCentralityStack->SetMinimum(centralityLowerBound);
  piCentralityStack->SetMaximum(centralityUpperBound);
  piCentralityStack->Draw("NOSTACK E1P");
  zeroLine->Draw("SAME");
  piLegend->Draw();
  piText->Draw();
  copyWithNewErrors1->Draw("[]");
  copyWithNewErrors2->Draw("[]");
  canvas->SaveAs("sys_h_vn_pi.png");
  delete copyWithNewErrors1;
  delete copyWithNewErrors2;
  canvas->Clear();
  //===

  //=== K+- vs centrality
  copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_kp->Clone());
  xAxisError = Normal->h_vn_kp->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_kp.size(); i++)
    { copyWithNewErrors1->SetPointError(i, xAxisError, v_sys_kp.at(i)); }

  copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_km->Clone());
  xAxisError = Normal->h_vn_km->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_km.size(); i++)
    { copyWithNewErrors2->SetPointError(i, xAxisError, v_sys_km.at(i)); }

  copyWithNewErrors1->Write();
  copyWithNewErrors2->Write();

  THStack *kaCentralityStack = new THStack("kaCentralityStack", ";Centrality (%);v_{"+order_n_str+"}");
  kaCentralityStack->Add(Normal->h_vn_kp);
  kaCentralityStack->Add(Normal->h_vn_km);

  kaCentralityStack->Draw();
  kaCentralityStack->SetMinimum(centralityLowerBound);
  kaCentralityStack->SetMaximum(centralityUpperBound);
  kaCentralityStack->Draw("NOSTACK E1P");
  zeroLine->Draw("SAME");
  kaLegend->Draw();
  kaText->Draw();
  copyWithNewErrors1->Draw("[]");
  copyWithNewErrors2->Draw("[]");
  canvas->SaveAs("sys_h_vn_ka.png");
  delete copyWithNewErrors1;
  delete copyWithNewErrors2;
  canvas->Clear();
  //===

  //=== Proton vs centrality
  copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pr->Clone());
  xAxisError = Normal->h_vn_pr->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_pr.size(); i++)
    { copyWithNewErrors1->SetPointError(i, xAxisError, v_sys_pr.at(i)); }

  copyWithNewErrors1->Write();
      
  Normal->h_vn_pr->Draw();
  Normal->h_vn_pr->SetMinimum(centralityLowerBound);
  Normal->h_vn_pr->SetMaximum(centralityUpperBound);
  Normal->h_vn_pr->Draw("E1");
  copyWithNewErrors1->Draw("[]");
  zeroLine->Draw("SAME");
  prText->Draw();
  canvas->SaveAs("sys_h_vn_pr.png");
  delete copyWithNewErrors1;
  canvas->Clear();
  //===




  //=== Proton vs rapidity
  THStack *prRapidityStack = new THStack("prRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
  prRapidityStack->Add(Normal->h_vn_yCM_00to10_pr);
  prRapidityStack->Add(Normal->h_vn_yCM_10to40_pr);
  prRapidityStack->Add(Normal->h_vn_yCM_40to60_pr);

  copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_pr->Clone());
  xAxisError = Normal->h_vn_yCM_00to10_pr->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_yCM_00to10_pr.size(); i++)
    { copyWithNewErrors1->SetPointError(i, xAxisError, v_sys_yCM_00to10_pr.at(i)); }

  copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_pr->Clone());
  xAxisError = Normal->h_vn_yCM_10to40_pr->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_yCM_10to40_pr.size(); i++)
    { copyWithNewErrors2->SetPointError(i, xAxisError, v_sys_yCM_10to40_pr.at(i)); }

  copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_pr->Clone());
  xAxisError = Normal->h_vn_yCM_40to60_pr->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_yCM_40to60_pr.size(); i++)
    { copyWithNewErrors3->SetPointError(i, xAxisError, v_sys_yCM_40to60_pr.at(i)); }

  copyWithNewErrors1->Write();
  copyWithNewErrors2->Write();
  copyWithNewErrors3->Write();
      
  prRapidityStack->Draw();
  prRapidityStack->GetYaxis()->SetTitleOffset(1.9);
  prRapidityStack->GetXaxis()->SetNdivisions(210);
  prRapidityStack->SetMaximum(rapidityUpperBound_pr);
  prRapidityStack->SetMinimum(rapidityLowerBound_pr);
  prRapidityStack->Draw("NOSTACK E1P");
  zeroLine_y_pr->Draw("SAME");
  copyWithNewErrors1->Draw("[]");
  copyWithNewErrors2->Draw("[]");
  copyWithNewErrors3->Draw("[]");
  prLegend->Draw();
  prText_y->Draw();
  canvas->SaveAs("sys_prRapidityStack.png");
  canvas->Clear();
  delete copyWithNewErrors1;
  delete copyWithNewErrors2;
  delete copyWithNewErrors3;
  //===
      

  //=== Proton vs rapidity symmetric across midrapidity
  THStack *prRapidityStack_symm = new THStack("prRapidityStack_symm", ";y-y_{mid};v_{"+order_n_str+"}");
  prRapidityStack_symm->Add(Normal->h_vn_yCM_00to10_pr_symm);
  prRapidityStack_symm->Add(Normal->h_vn_yCM_10to40_pr_symm);
  prRapidityStack_symm->Add(Normal->h_vn_yCM_40to60_pr_symm);

  copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_pr_symm->Clone());
  xAxisError = Normal->h_vn_yCM_00to10_pr_symm->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_yCM_00to10_pr_symm.size(); i++)
    { copyWithNewErrors1->SetPointError(i, xAxisError, v_sys_yCM_00to10_pr_symm.at(i)); }

  copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_pr_symm->Clone());
  xAxisError = Normal->h_vn_yCM_10to40_pr_symm->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_yCM_10to40_pr_symm.size(); i++)
    { copyWithNewErrors2->SetPointError(i, xAxisError, v_sys_yCM_10to40_pr_symm.at(i)); }

  copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_pr_symm->Clone());
  xAxisError = Normal->h_vn_yCM_40to60_pr_symm->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_yCM_40to60_pr_symm.size(); i++)
    { copyWithNewErrors3->SetPointError(i, xAxisError, v_sys_yCM_40to60_pr_symm.at(i)); }

  copyWithNewErrors1->Write();
  copyWithNewErrors2->Write();
  copyWithNewErrors3->Write();
      
  prRapidityStack_symm->Draw();
  prRapidityStack_symm->GetYaxis()->SetTitleOffset(1.9);
  prRapidityStack_symm->GetXaxis()->SetNdivisions(210);
  prRapidityStack_symm->SetMaximum(rapidityUpperBound_pr);
  prRapidityStack_symm->SetMinimum(rapidityLowerBound_pr);
  prRapidityStack_symm->Draw("NOSTACK E1P");
  zeroLine_y_pr->Draw("SAME");
  copyWithNewErrors1->Draw("[]");
  copyWithNewErrors2->Draw("[]");
  copyWithNewErrors3->Draw("[]");
  prLegend->Draw();
  prText_y_symm->Draw();
  canvas->SaveAs("sys_pr_symmRapidityStack.png");
  canvas->Clear();
  delete copyWithNewErrors1;
  delete copyWithNewErrors2;
  delete copyWithNewErrors3;
  //===






  //=== Proton vs pT
  THStack *prPtStack = new THStack("prPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
  prPtStack->Add(Normal->h_vn_pT_00to10_pr);
  prPtStack->Add(Normal->h_vn_pT_40to60_pr);
  prPtStack->Add(Normal->h_vn_pT_10to40_pr);

  copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_pr->Clone());
  xAxisError = Normal->h_vn_pT_00to10_pr->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_pT_00to10_pr.size(); i++)
    { copyWithNewErrors1->SetPointError(i, xAxisError, v_sys_pT_00to10_pr.at(i)); }

  copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_pr->Clone());
  xAxisError = Normal->h_vn_pT_10to40_pr->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_pT_10to40_pr.size(); i++)
    { copyWithNewErrors2->SetPointError(i, xAxisError, v_sys_pT_10to40_pr.at(i)); }

  copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_pr->Clone());
  xAxisError = Normal->h_vn_pT_40to60_pr->GetXaxis()->GetBinWidth(1)/2.0;
  for (int i = 0; i < v_sys_pT_40to60_pr.size(); i++)
    { copyWithNewErrors3->SetPointError(i, xAxisError, v_sys_pT_40to60_pr.at(i)); }

  copyWithNewErrors1->Write();
  copyWithNewErrors2->Write();
  copyWithNewErrors3->Write();
      
  prPtStack->Draw();
  prPtStack->GetYaxis()->SetTitleOffset(1.9);
  prPtStack->GetXaxis()->SetNdivisions(210);
  prPtStack->SetMaximum(ptUpperBound);
  prPtStack->SetMinimum(ptLowerBound);
  prPtStack->Draw("NOSTACK E1P");
  zeroLine_pt->Draw("SAME");
  copyWithNewErrors1->Draw("[]");
  copyWithNewErrors3->Draw("[]");
  copyWithNewErrors2->Draw("[]");
  prPtLegend->Draw();
  prPtText->Draw();
  canvas->SaveAs("sys_prPtStack.png");
  canvas->Clear();
  delete copyWithNewErrors1;
  delete copyWithNewErrors2;
  delete copyWithNewErrors3;
  //===

  delete canvas;



  delete epd;
  delete nhits;
  delete nSigPi;
  delete nSigKa;
  delete nSigPr;
  //delete rvtx;
  //delete zvtx;
  delete dca;
  delete nhitsdEdx;
  delete nhitsratio;
  delete m2Pi;
  delete m2Ka;

  delete Normal;
  //delete epd_high;
  //delete epd_low;
  //delete epd_scaled;
  
  delete nSigPi_high_30;
  delete nSigPi_low_30;
  delete nSigKa_high_30;
  delete nSigKa_low_30;
  delete nSigPr_high_30;
  delete nSigPr_low_30;
  //delete rvtx_high_30;
  //delete rvtx_low_30;
  //delete zvtx_high_30;
  //delete zvtx_low_30;
  delete dca_high_30;
  delete dca_low_30;
  delete nhits_high_30;
  delete nhits_low_30;
  delete nhitsdEdx_high_30;
  //delete nhitsdEdx_low_30;
  delete nhitsratio_high_30;
  delete nhitsratio_low_30;
  delete m2Pi_high_30;
  delete m2Pi_low_30;
  delete m2Ka_high_30;
  delete m2Ka_low_30;


  delete nSigPi_high_20;
  delete nSigPi_low_20;
  delete nSigKa_high_20;
  delete nSigKa_low_20;
  delete nSigPr_high_20;
  delete nSigPr_low_20;
  //delete rvtx_high_20;
  //delete rvtx_low_20;
  //delete zvtx_high_20;
  //delete zvtx_low_20;
  delete dca_high_20;
  delete dca_low_20;
  delete nhits_high_20;
  delete nhits_low_20;
  delete nhitsdEdx_high_20;
  //delete nhitsdEdx_low_20;
  delete nhitsratio_high_20;
  delete nhitsratio_low_20;
  delete m2Pi_high_20;
  delete m2Pi_low_20;
  delete m2Ka_high_20;
  delete m2Ka_low_20;

  newFile->Close();
  delete newFile;
}



void printPointErrors(TGraphErrors *graph)
{
  Int_t nPoints = graph->GetN();

  std::cout << graph->GetName() << ":" << std::endl;
  
  for (int i = 0; i < nPoints; i++)
    {
      std::cout << "Bin " << i+1 << ": " << graph->GetErrorY(i) << std::endl;
    }
}

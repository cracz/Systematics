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

  Int_t epdNbins = 0;
  Int_t nhitsNbins = 0;
  Int_t nSigPiNbins = 0;
  Int_t nSigKaNbins = 0;
  Int_t nSigPrNbins = 0.0;
  Int_t rvtxNbins = 0.0;
  Int_t zvtxNbins = 0.0;
  Int_t dcaNbins = 0.0;
  Int_t nhitsdEdxNbins = 0.0;
  Int_t nhitsratioNbins = 0.0;
  Int_t m2PiNbins = 0.0;
  Int_t m2KaNbins = 0.0;

  void addContribution(TString ID, Double_t stdDevContributed)
  {
    if (ID == "epd")
    {
      epdPercentQuadSum += stdDevContributed;
      epdNbins++;
    }
    else if (ID == "nhits")
    {
      nhitsPercentQuadSum += stdDevContributed;
      nhitsNbins++;
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
    << "EP Resolution, " << epdPercentQuadSum / (Double_t)epdNbins << std::endl
    << "nHits, " << nhitsPercentQuadSum / (Double_t)nhitsNbins << std::endl
    << "nHits dEdx, " << nhitsdEdxPercentQuadSum / (Double_t)nhitsdEdxNbins << std::endl
    << "nHits Ratio, " << nhitsratioPercentQuadSum / (Double_t)nhitsratioNbins << std::endl
    << "DCA, " << dcaPercentQuadSum / (Double_t)dcaNbins << std::endl
    << "r Vertex, " << rvtxPercentQuadSum / (Double_t)rvtxNbins << std::endl
    << "z Vertex, " << zvtxPercentQuadSum / (Double_t)zvtxNbins << std::endl
    << "nSigma Pi, " << nSigPiPercentQuadSum / (Double_t)nSigPiNbins << std::endl
    << "nSigma Ka, " << nSigKaPercentQuadSum / (Double_t)nSigKaNbins << std::endl
    << "nSigma Pr, " << nSigPrPercentQuadSum / (Double_t)nSigPrNbins << std::endl
    << "m^2 Pi, " << m2PiPercentQuadSum / (Double_t)m2PiNbins << std::endl
    << "m^2 Ka, " << m2KaPercentQuadSum / (Double_t)m2KaNbins << std::endl
    << std::endl;
  }
};

void calculateSystematics(TString order_n_str)
{
  TFile* newFile = new TFile("systematicErrors.root", "RECREATE");
  
  Variation* Normal = new Variation("Normal_30", order_n_str);

  Variation* epd_high_30 = new Variation("epd_high_30", order_n_str);
  //Variation* epd_low_30 = new Variation("E5FCCCABE8E37F613B8E2F76BE72E3AD", order_n_str);
  //Variation* epd_low_30 = new Variation("epd_low_30", order_n_str);
  Variation* epd_scaled = new Variation("epd_scaled", order_n_str);
  Variation* nSigPi_high_30 = new Variation("nSigPi_high_30", order_n_str);
  Variation* nSigPi_low_30  = new Variation("nSigPi_low_30", order_n_str);
  Variation* nSigKa_high_30 = new Variation("nSigKa_high_30", order_n_str);
  Variation* nSigKa_low_30  = new Variation("nSigKa_low_30", order_n_str);
  Variation* nSigPr_high_30 = new Variation("nSigPr_high_30", order_n_str);
  Variation* nSigPr_low_30  = new Variation("nSigPr_low_30", order_n_str);
  Variation* rvtx_high_30 = new Variation("rvtx_high_30", order_n_str);
  Variation* rvtx_low_30  = new Variation("rvtx_low_30", order_n_str);
  Variation* zvtx_high_30 = new Variation("zvtx_high_30", order_n_str);
  Variation* zvtx_low_30  = new Variation("zvtx_low_30", order_n_str);
  Variation* dca_high_30 = new Variation("dca_high_30", order_n_str);
  Variation* dca_low_30  = new Variation("dca_low_30", order_n_str);
  Variation* nhits_high_30 = new Variation("nhits_high_30", order_n_str);
  Variation* nhits_low_30 = new Variation("nhits_low_30", order_n_str);
  Variation* nhitsdEdx_high_30 = new Variation("nhitsdEdx_high_30", order_n_str);
  //Variation* nhitsdEdx_low_30  = new Variation("nhitsdEdx_low_30", order_n_str);
  Variation* nhitsratio_high_30 = new Variation("nhitsratio_high_30", order_n_str);
  Variation* nhitsratio_low_30  = new Variation("nhitsratio_low_30", order_n_str);
  Variation* m2Pi_high_30 = new Variation("m2Pi_high_30", order_n_str);
  Variation* m2Pi_low_30  = new Variation("m2Pi_low_30", order_n_str);
  Variation* m2Ka_high_30 = new Variation("m2Ka_high_30", order_n_str);
  Variation* m2Ka_low_30  = new Variation("m2Ka_low_30", order_n_str);


  //Variation* Normal_20 = new Variation("Normal_20", order_n_str);
  //Variation* epd_high_20 = new Variation("epd_high_20", order_n_str);
  //Variation* epd_low_20 = new Variation("epd_low_20", order_n_str);
  Variation* nSigPi_high_20 = new Variation("nSigPi_high_20", order_n_str);
  Variation* nSigPi_low_20  = new Variation("nSigPi_low_20", order_n_str);
  Variation* nSigKa_high_20 = new Variation("nSigKa_high_20", order_n_str);
  Variation* nSigKa_low_20  = new Variation("nSigKa_low_20", order_n_str);
  Variation* nSigPr_high_20 = new Variation("nSigPr_high_20", order_n_str);
  Variation* nSigPr_low_20  = new Variation("nSigPr_low_20", order_n_str);
  Variation* rvtx_high_20 = new Variation("rvtx_high_20", order_n_str);
  Variation* rvtx_low_20  = new Variation("rvtx_low_20", order_n_str);
  Variation* zvtx_high_20 = new Variation("zvtx_high_20", order_n_str);
  Variation* zvtx_low_20  = new Variation("zvtx_low_20", order_n_str);
  Variation* dca_high_20 = new Variation("dca_high_20", order_n_str);
  Variation* dca_low_20  = new Variation("dca_low_20", order_n_str);
  Variation* nhits_high_20 = new Variation("nhits_high_20", order_n_str);
  Variation* nhits_low_20 = new Variation("nhits_low_20", order_n_str);
  Variation* nhitsdEdx_high_20 = new Variation("nhitsdEdx_high_20", order_n_str);
  //Variation* nhitsdEdx_low_20  = new Variation("nhitsdEdx_low_20", order_n_str);
  Variation* nhitsratio_high_20 = new Variation("nhitsratio_high_20", order_n_str);
  Variation* nhitsratio_low_20  = new Variation("nhitsratio_low_20", order_n_str);
  Variation* m2Pi_high_20 = new Variation("m2Pi_high_20", order_n_str);
  Variation* m2Pi_low_20  = new Variation("m2Pi_low_20", order_n_str);
  Variation* m2Ka_high_20 = new Variation("m2Ka_high_20", order_n_str);
  Variation* m2Ka_low_20  = new Variation("m2Ka_low_20", order_n_str);


  //CompositeData* epd = new CompositeData("epd", Normal, epd_low_30, epd_high_30);
  CompositeData* epd = new CompositeData("epd", Normal, epd_scaled, epd_high_30);
  CompositeData* nhits = new CompositeData("nhits", Normal, nhits_low_30, nhits_high_30, nhits_low_20, nhits_high_20);  
  CompositeData* nSigPi = new CompositeData("nSigPi", Normal, nSigPi_low_30, nSigPi_high_30, nSigPi_low_20, nSigPi_high_20);
  CompositeData* nSigKa = new CompositeData("nSigKa", Normal, nSigKa_low_30, nSigKa_high_30, nSigKa_low_20, nSigKa_high_20);
  CompositeData* nSigPr = new CompositeData("nSigPr", Normal, nSigPr_low_30, nSigPr_high_30, nSigPr_low_20, nSigPr_high_20);
  CompositeData* rvtx = new CompositeData("rvtx", Normal, rvtx_low_30, rvtx_high_30, rvtx_low_20, rvtx_high_20);
  CompositeData* zvtx = new CompositeData("zvtx", Normal, zvtx_low_30, zvtx_high_30, zvtx_low_20, zvtx_high_20);
  CompositeData* dca  = new CompositeData("dca", Normal, dca_low_30, dca_high_30, dca_low_20, dca_high_20);
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
  composites.push_back(rvtx);
  composites.push_back(zvtx);
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

  AverageContributionTracker avgTracker;

  //=== pi+ vs centrality
  std::vector<Double_t> v_sys_pp;
  bins = Normal->h_vn_pp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    //std::cout << "Bin " << ithBin+1 << std::endl;
    quadSum = 0.0;
    quadSum += epd->v_vn_pp.at(ithBin).variance;

    avgTracker.addContribution(epd->ID, epd->v_vn_pp.at(ithBin).stdDevPercentage);
    //epdPercentQuadSum += epd->v_vn_pp.at(ithBin).stdDevPercentage;
    //epdNbins += 1;
    //std::cout << "epd: " << epd->v_vn_pp.at(ithBin).stdDevPercentage << std::endl;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pp.at(ithBin).deltaByDeltaError > 1.0)
        {
          quadSum += composites.at(jthVariation)->v_vn_pp.at(ithBin).variance;
          avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_pp.at(ithBin).stdDevPercentage);
        }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_pp loop

  //=== pi- vs centrality
  std::vector<Double_t> v_sys_pm;
  bins = Normal->h_vn_pm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    //std::cout << "Bin " << ithBin << std::endl;
    quadSum = 0.0;
    quadSum += epd->v_vn_pm.at(ithBin).variance;

    avgTracker.addContribution(epd->ID, epd->v_vn_pm.at(ithBin).stdDevPercentage);
    //epdPercentQuadSum += epd->v_vn_pm.at(ithBin).stdDevPercentage;
    //epdNbins += 1;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pm.at(ithBin).deltaByDeltaError > 1.0)
        {
          quadSum += composites.at(jthVariation)->v_vn_pm.at(ithBin).variance;
          avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_pm.at(ithBin).stdDevPercentage);
        }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pm.push_back(ithBinSysErr);
  } 
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_pm loop

  //=== K+ vs centrality
  std::vector<Double_t> v_sys_kp;
  bins = Normal->h_vn_kp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_kp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_kp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_kp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_kp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_kp loop

  //=== K- vs centrality
  std::vector<Double_t> v_sys_km;
  bins = Normal->h_vn_km->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_km.at(ithBin).variance;


    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_km.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_km.at(ithBin).variance;
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
    //std::cout << "Bin " << ithBin+1 << std::endl;
    quadSum = 0.0;
    quadSum += epd->v_vn_pr.at(ithBin).variance;
    //std::cout << 100 * epd->v_vn_pr.at(ithBin).stdDev / TMath::Abs(Normal->h_vn_pr->GetBinContent(ithBin+1)) << std::endl;
    //std::cout << "epd: " << epd->v_vn_pr.at(ithBin).stdDevPercentage << std::endl;
    
    avgTracker.addContribution(epd->ID, epd->v_vn_pr.at(ithBin).stdDevPercentage);
    //epdPercentQuadSum += epd->v_vn_pr.at(ithBin).stdDevPercentage;
    //epdNbins += 1;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pr.at(ithBin).deltaByDeltaError > 1.0)
        {
          quadSum += composites.at(jthVariation)->v_vn_pr.at(ithBin).variance;
          avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_pr.at(ithBin).stdDevPercentage);
        }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_pr loop

/*
  //=== proton alternate acceptance vs centrality
  std::cout << "PROTON VS CENTRALITY (ALTERNATE ACCEPTANCE)" << std::endl;
  std::vector<Double_t> v_sys_pr_alt;
  bins = Normal->h_vn_pr_alt->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    std::cout << "Bin " << ithBin << std::endl;
    quadSum = 0.0;
    quadSum += epd->v_vn_pr_alt.at(ithBin).variance;


    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pr_alt.at(ithBin).deltaByDeltaError > 1.0)
        {
          quadSum += composites.at(jthVariation)->v_vn_pr_alt.at(ithBin).variance;

          std::cout << composites.at(jthVariation)->ID << ": "
        << (composites.at(jthVariation)->v_vn_pr_alt.at(ithBin).stdDev/Normal->h_vn_pr_alt->GetBinContent(ithBin+1) ) * 100
        << std::endl;
        }
    }

    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pr_alt.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_pr_alt loop



  //=== Deuteron vs centrality
  std::cout << "DEUTERON VS CENTRALITY" << std::endl;
  std::vector<Double_t> v_sys_de;
  bins = Normal->h_vn_de->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    std::cout << "Bin " << ithBin << std::endl;
    quadSum = 0.0;
    quadSum += epd->v_vn_de.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_de.at(ithBin).deltaByDeltaError > 1.0)
        {
          quadSum += composites.at(jthVariation)->v_vn_de.at(ithBin).variance;

          std::cout << composites.at(jthVariation)->ID << ": "
        << (composites.at(jthVariation)->v_vn_de.at(ithBin).stdDev/Normal->h_vn_de->GetBinContent(ithBin+1) ) * 100
        << std::endl;
        }
    }

    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_de.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_de loop


  //=== Triton vs centrality
  std::cout << "TRITON VS CENTRALITY" << std::endl;
  std::vector<Double_t> v_sys_tr;
  bins = Normal->h_vn_tr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    std::cout << "Bin " << ithBin << std::endl;
    quadSum = 0.0;
    quadSum += epd->v_vn_tr.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_tr.at(ithBin).deltaByDeltaError > 1.0)
        {
          quadSum += composites.at(jthVariation)->v_vn_tr.at(ithBin).variance;

          std::cout << composites.at(jthVariation)->ID << ": "
        << (composites.at(jthVariation)->v_vn_tr.at(ithBin).stdDev/Normal->h_vn_tr->GetBinContent(ithBin+1) ) * 100
        << std::endl;
        }
    }

    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_tr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_tr loop
*/


  
  //=== pi+ vs rapidity 0 - 10%
  std::vector<Double_t> v_sys_yCM_00to10_pp;
  bins = Normal->h_vn_yCM_00to10_pp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_00to10_pp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_00to10_pp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_00to10_pp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_00to10_pp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_pp loop


  //=== pi+ vs rapidity 10 - 40%
  std::vector<Double_t> v_sys_yCM_10to40_pp;
  bins = Normal->h_vn_yCM_10to40_pp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_10to40_pp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_10to40_pp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_10to40_pp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_10to40_pp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_pp loop


  //=== pi+ vs rapidity 40 - 60%
  std::vector<Double_t> v_sys_yCM_40to60_pp;
  bins = Normal->h_vn_yCM_40to60_pp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_40to60_pp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_40to60_pp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_40to60_pp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_40to60_pp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //=== End h_vn_pp loop


  //=== pi- vs rapidity 0 - 10%
  std::vector<Double_t> v_sys_yCM_00to10_pm;
  bins = Normal->h_vn_yCM_00to10_pm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_00to10_pm.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_00to10_pm.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_00to10_pm.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_00to10_pm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===

  //=== pi- vs rapidity 10 - 40%
  std::vector<Double_t> v_sys_yCM_10to40_pm;
  bins = Normal->h_vn_yCM_10to40_pm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_10to40_pm.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_10to40_pm.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_10to40_pm.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_10to40_pm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===

  //=== pi- vs rapidity 40 - 60%
  std::vector<Double_t> v_sys_yCM_40to60_pm;
  bins = Normal->h_vn_yCM_40to60_pm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_40to60_pm.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_40to60_pm.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_40to60_pm.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_40to60_pm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K+ vs rapidity 0 - 10%
  std::vector<Double_t> v_sys_yCM_00to10_kp;
  bins = Normal->h_vn_yCM_00to10_kp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_00to10_kp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_00to10_kp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_00to10_kp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_00to10_kp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K+ vs rapidity 10 - 40%
  std::vector<Double_t> v_sys_yCM_10to40_kp;
  bins = Normal->h_vn_yCM_10to40_kp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_10to40_kp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_10to40_kp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_10to40_kp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_10to40_kp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K+ vs rapidity 40 - 60%
  std::vector<Double_t> v_sys_yCM_40to60_kp;
  bins = Normal->h_vn_yCM_40to60_kp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_40to60_kp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_40to60_kp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_40to60_kp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_40to60_kp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K- vs rapidity 0 - 10%
  std::vector<Double_t> v_sys_yCM_00to10_km;
  bins = Normal->h_vn_yCM_00to10_km->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_00to10_km.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_00to10_km.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_00to10_km.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_00to10_km.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===

  //=== K- vs rapidity 10 - 40%
  std::vector<Double_t> v_sys_yCM_10to40_km;
  bins = Normal->h_vn_yCM_10to40_km->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_10to40_km.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_10to40_km.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_10to40_km.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_10to40_km.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K- vs rapidity 40 - 60%
  std::vector<Double_t> v_sys_yCM_40to60_km;
  bins = Normal->h_vn_yCM_40to60_km->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_40to60_km.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_40to60_km.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_yCM_40to60_km.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_40to60_km.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== Proton vs rapidity 0 - 10%
  std::vector<Double_t> v_sys_yCM_00to10_pr;
  bins = Normal->h_vn_yCM_00to10_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_00to10_pr.at(ithBin).variance;
    if (ithBin > 9)
    {
      avgTracker.addContribution(epd->ID, epd->v_vn_yCM_00to10_pr.at(ithBin).stdDevPercentage);
      //epdPercentQuadSum += epd->v_vn_yCM_00to10_pr.at(ithBin).stdDevPercentage;
      //epdNbins += 1;
    }

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_00to10_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthVariation)->v_vn_yCM_00to10_pr.at(ithBin).variance;
        avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_yCM_00to10_pr.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_00to10_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== Proton vs rapidity 10 - 40%
  std::vector<Double_t> v_sys_yCM_10to40_pr;
  bins = Normal->h_vn_yCM_10to40_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_10to40_pr.at(ithBin).variance;
    
    if (ithBin > 9)
    {
      avgTracker.addContribution(epd->ID, epd->v_vn_yCM_10to40_pr.at(ithBin).stdDevPercentage);
      //epdPercentQuadSum += epd->v_vn_yCM_10to40_pr.at(ithBin).stdDevPercentage;
      //epdNbins += 1;
    }


    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_10to40_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthVariation)->v_vn_yCM_10to40_pr.at(ithBin).variance;
        avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_yCM_10to40_pr.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_10to40_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== Proton vs rapidity 40 - 60%
  std::vector<Double_t> v_sys_yCM_40to60_pr;
  bins = Normal->h_vn_yCM_40to60_pr->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_yCM_40to60_pr.at(ithBin).variance;
    
    if (ithBin > 9)
    {
      avgTracker.addContribution(epd->ID, epd->v_vn_yCM_40to60_pr.at(ithBin).stdDevPercentage);
      //epdPercentQuadSum += epd->v_vn_yCM_40to60_pr.at(ithBin).stdDevPercentage;
      //epdNbins += 1;
    }

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_40to60_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthVariation)->v_vn_yCM_40to60_pr.at(ithBin).variance;
        avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_yCM_40to60_pr.at(ithBin).stdDevPercentage);
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
    {
      avgTracker.addContribution(epd->ID, epd->v_vn_yCM_00to10_pr_symm.at(ithBin).stdDevPercentage);
      //epdPercentQuadSum += epd->v_vn_yCM_00to10_pr_symm.at(ithBin).stdDevPercentage;
      //epdNbins += 1;
    }

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_00to10_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthVariation)->v_vn_yCM_00to10_pr_symm.at(ithBin).variance;
        avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_yCM_00to10_pr_symm.at(ithBin).stdDevPercentage);
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
    {
      avgTracker.addContribution(epd->ID, epd->v_vn_yCM_10to40_pr_symm.at(ithBin).stdDevPercentage);
      //epdPercentQuadSum += epd->v_vn_yCM_10to40_pr_symm.at(ithBin).stdDevPercentage;
      //epdNbins += 1;
    }

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_10to40_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthVariation)->v_vn_yCM_10to40_pr_symm.at(ithBin).variance;
        avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_yCM_10to40_pr_symm.at(ithBin).stdDevPercentage);
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
    {
      avgTracker.addContribution(epd->ID, epd->v_vn_yCM_40to60_pr_symm.at(ithBin).stdDevPercentage);
      //epdPercentQuadSum += epd->v_vn_yCM_40to60_pr_symm.at(ithBin).stdDevPercentage;
      //epdNbins += 1;
    }

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_yCM_40to60_pr_symm.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthVariation)->v_vn_yCM_40to60_pr_symm.at(ithBin).variance;
        avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_yCM_40to60_pr_symm.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_yCM_40to60_pr_symm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== pi+ vs pT 0 - 10%
  std::vector<Double_t> v_sys_pT_00to10_pp;
  bins = Normal->h_vn_pT_00to10_pp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_00to10_pp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_00to10_pp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_00to10_pp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_00to10_pp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== pi+ vs pT 10 - 40%
  std::vector<Double_t> v_sys_pT_10to40_pp;
  bins = Normal->h_vn_pT_10to40_pp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_10to40_pp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_10to40_pp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_10to40_pp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_10to40_pp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== pi+ vs pT 40 - 60%
  std::vector<Double_t> v_sys_pT_40to60_pp;
  bins = Normal->h_vn_pT_40to60_pp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_40to60_pp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_40to60_pp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_40to60_pp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_40to60_pp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== pi- vs pT 0 - 10%
  std::vector<Double_t> v_sys_pT_00to10_pm;
  bins = Normal->h_vn_pT_00to10_pm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_00to10_pm.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_00to10_pm.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_00to10_pm.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_00to10_pm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== pi- vs pT 10 - 40%
  std::vector<Double_t> v_sys_pT_10to40_pm;
  bins = Normal->h_vn_pT_10to40_pm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_10to40_pm.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_10to40_pm.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_10to40_pm.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_10to40_pm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== pi- vs pT 40 - 60%
  std::vector<Double_t> v_sys_pT_40to60_pm;
  bins = Normal->h_vn_pT_40to60_pm->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_40to60_pm.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_40to60_pm.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_40to60_pm.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_40to60_pm.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K+ vs pT 0 - 10%
  std::vector<Double_t> v_sys_pT_00to10_kp;
  bins = Normal->h_vn_pT_00to10_kp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_00to10_kp.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_00to10_kp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_00to10_kp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_00to10_kp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K+ vs pT 10 - 40%
  std::vector<Double_t> v_sys_pT_10to40_kp;
  bins = Normal->h_vn_pT_10to40_kp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_10to40_kp.at(ithBin).variance;


    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_10to40_kp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_10to40_kp.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_10to40_kp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K+ vs pT 40 - 60%
  std::vector<Double_t> v_sys_pT_40to60_kp;
  bins = Normal->h_vn_pT_40to60_kp->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_40to60_kp.at(ithBin).variance;


    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_40to60_kp.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_40to60_kp.at(ithBin).variance;
    }

    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_40to60_kp.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K- vs pT 0 - 10%
  std::vector<Double_t> v_sys_pT_00to10_km;
  bins = Normal->h_vn_pT_00to10_km->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_00to10_km.at(ithBin).variance;


    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_00to10_km.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_00to10_km.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_00to10_km.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K- vs pT 10 - 40%
  std::vector<Double_t> v_sys_pT_10to40_km;
  bins = Normal->h_vn_pT_10to40_km->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_10to40_km.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_10to40_km.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_10to40_km.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_10to40_km.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===


  //=== K- vs pT 40 - 60%
  std::vector<Double_t> v_sys_pT_40to60_km;
  bins = Normal->h_vn_pT_40to60_km->GetNbinsX();
  for (int ithBin = 0; ithBin < bins; ithBin++)
  {
    quadSum = 0.0;
    quadSum += epd->v_vn_pT_40to60_km.at(ithBin).variance;

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_40to60_km.at(ithBin).deltaByDeltaError > 1.0)
        quadSum += composites.at(jthVariation)->v_vn_pT_40to60_km.at(ithBin).variance;
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_40to60_km.push_back(ithBinSysErr);
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
    {
      avgTracker.addContribution(epd->ID, epd->v_vn_pT_00to10_pr.at(ithBin).stdDevPercentage);
      //epdPercentQuadSum += epd->v_vn_pT_00to10_pr.at(ithBin).stdDevPercentage;
      //epdNbins += 1;
    }

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_00to10_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthVariation)->v_vn_pT_00to10_pr.at(ithBin).variance;
        avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_pT_00to10_pr.at(ithBin).stdDevPercentage);
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
    {
      avgTracker.addContribution(epd->ID, epd->v_vn_pT_10to40_pr.at(ithBin).stdDevPercentage);
      //epdPercentQuadSum += epd->v_vn_pT_10to40_pr.at(ithBin).stdDevPercentage;
      //epdNbins += 1;
    }

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_10to40_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthVariation)->v_vn_pT_10to40_pr.at(ithBin).variance;
        avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_pT_10to40_pr.at(ithBin).stdDevPercentage);
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
    {
      avgTracker.addContribution(epd->ID, epd->v_vn_pT_40to60_pr.at(ithBin).stdDevPercentage);
      //epdPercentQuadSum += epd->v_vn_pT_40to60_pr.at(ithBin).stdDevPercentage;
      //epdNbins += 1;
    }

    for (int jthVariation = 0; jthVariation < composites.size(); jthVariation++)
    {
      if (composites.at(jthVariation)->v_vn_pT_40to60_pr.at(ithBin).deltaByDeltaError > 1.0)
      {
        quadSum += composites.at(jthVariation)->v_vn_pT_40to60_pr.at(ithBin).variance;
        avgTracker.addContribution(composites.at(jthVariation)->ID, composites.at(jthVariation)->v_vn_pT_40to60_pr.at(ithBin).stdDevPercentage);
      }
    }
    
    ithBinSysErr = TMath::Sqrt(quadSum);
    v_sys_pT_40to60_pr.push_back(ithBinSysErr);
  }
  ithBinSysErr = 0;
  quadSum = 0.0;
  //===



  
  /*
  // Proton vs yCM scan through pT
  std::vector<Double_t> v_sys_yCM_pr_1;
  bins = Normal->h_vn_yCM_pr_1->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(epd_high->v_vn_yCM_pr_1.at(i).variance);
				 
      v_sys_yCM_pr_1.push_back(ithBinSysErr);
    }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_pr_2;
  bins = Normal->h_vn_yCM_pr_2->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(epd_high->v_vn_yCM_pr_2.at(i).variance);
				 
      v_sys_yCM_pr_2.push_back(ithBinSysErr);
    }
  ithBinSysErr = 0;
  quadSum = 0.0;

  std::vector<Double_t> v_sys_yCM_pr_3;
  bins = Normal->h_vn_yCM_pr_3->GetNbinsX();
  for (int i = 0; i < bins; i++)
    {
      ithBinSysErr = TMath::Sqrt(epd_high->v_vn_yCM_pr_3.at(i).variance);
				 
      v_sys_yCM_pr_3.push_back(ithBinSysErr);
    }
  ithBinSysErr = 0;
  quadSum = 0.0;
  */
  
  
  avgTracker.printContributions();


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


  if (order_n_str == "3")
    {      
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

      /*
      TLegend *ppExtLegend = new TLegend(0.4, 0.62, 0.7, 0.82);
      ppExtLegend->AddEntry(vn_pp,"#pi^{+}, 0 < y_{CM} < 0.5");
      ppExtLegend->AddEntry(vn_pp_ext,"#pi^{+}, 0.5 < y_{CM} < 1.0");
      ppExtLegend->SetFillColorAlpha(0,0);
      ppExtLegend->SetLineColorAlpha(0,0);

      TLegend *pmExtLegend = new TLegend(0.15, 0.67, 0.45, 0.9);
      pmExtLegend->AddEntry(vn_pm,"#pi^{-}, 0 < y_{CM} < 0.5");
      pmExtLegend->AddEntry(vn_pm_ext,"#pi^{-}, 0.5 < y_{CM} < 1.0");
      pmExtLegend->SetFillColorAlpha(0,0);
      pmExtLegend->SetLineColorAlpha(0,0);

      TLegend *kpExtLegend = new TLegend(0.55, 0.7, 0.85, 0.9);
      kpExtLegend->AddEntry(vn_kp,"K^{+}, 0 < y_{CM} < 0.5");
      kpExtLegend->AddEntry(vn_kp_ext,"K^{+}, 0.5 < y_{CM} < 1.0");
      kpExtLegend->SetFillColorAlpha(0,0);
      kpExtLegend->SetLineColorAlpha(0,0);

      TLegend *kmExtLegend = new TLegend(0.28, 0.68, 0.55, 0.85);
      kmExtLegend->AddEntry(vn_km,"K^{-}, 0 < y_{CM} < 0.5");
      kmExtLegend->AddEntry(vn_km_ext,"K^{-}, 0.5 < y_{CM} < 1.0");
      kmExtLegend->SetFillColorAlpha(0,0);
      kmExtLegend->SetLineColorAlpha(0,0);

      TLegend *prExtLegend = new TLegend(0.25, 0.16, 0.55, 0.3);
      prExtLegend->AddEntry(vn_pr_for,"Proton, -0.5 < y_{CM} < 0");
      prExtLegend->AddEntry(vn_pr,"Proton, 0 < y_{CM} < 0.5");
      prExtLegend->AddEntry(vn_pr_ext,"Proton, 0.5 < y_{CM} < 1.0");
      prExtLegend->SetFillColorAlpha(0,0);
      prExtLegend->SetLineColorAlpha(0,0);

      TLegend *etaLegend = new TLegend(0.65, 0.25, 0.9, 0.45);
      etaLegend->AddEntry(vn_EpdE, "EPD -5.6 < #eta < -3.3");
      etaLegend->AddEntry(vn_EpdF, "EPD -3.3 < #eta < -2.4");
      etaLegend->AddEntry(vn_TpcB, "TPC -1.0 < #eta < 0");
      */
      
      TLegend *ppLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      ppLegend->AddEntry(Normal->h_vn_yCM_00to10_pp, "0 - 10%");
      ppLegend->AddEntry(Normal->h_vn_yCM_10to40_pp, "10 - 40%");
      ppLegend->AddEntry(Normal->h_vn_yCM_40to60_pp, "40 - 60%");
      ppLegend->SetBorderSize(0);
      ppLegend->SetFillColorAlpha(0,0);

      TLegend *pmLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      pmLegend->AddEntry(Normal->h_vn_yCM_00to10_pm, "0 - 10%");
      pmLegend->AddEntry(Normal->h_vn_yCM_10to40_pm, "10 - 40%");
      pmLegend->AddEntry(Normal->h_vn_yCM_40to60_pm, "40 - 60%");
      pmLegend->SetBorderSize(0);
      pmLegend->SetFillColorAlpha(0,0);

      TLegend *kpLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      kpLegend->AddEntry(Normal->h_vn_yCM_00to10_kp, "0 - 10%");
      kpLegend->AddEntry(Normal->h_vn_yCM_10to40_kp, "10 - 40%");
      kpLegend->AddEntry(Normal->h_vn_yCM_40to60_kp, "40 - 60%");
      kpLegend->SetBorderSize(0);
      kpLegend->SetFillColorAlpha(0,0);

      TLegend *kmLegend = new TLegend(0.18, 0.77, 0.38, 0.87);
      //kmLegend->AddEntry(Normal->h_vn_yCM_00to10_km, "0 - 10%");
      kmLegend->AddEntry(Normal->h_vn_yCM_10to40_km, "10 - 40%");
      //kmLegend->AddEntry(Normal->h_vn_yCM_40to60_km, "40 - 60%");
      kmLegend->SetBorderSize(0);
      kmLegend->SetFillColorAlpha(0,0);

      TLegend *prLegend = new TLegend(0.19, 0.15, 0.39, 0.3);
      prLegend->AddEntry(Normal->h_vn_yCM_00to10_pr, "0 - 10%");
      prLegend->AddEntry(Normal->h_vn_yCM_10to40_pr, "10 - 40%");
      prLegend->AddEntry(Normal->h_vn_yCM_40to60_pr, "40 - 60%");
      prLegend->SetBorderSize(0);
      prLegend->SetFillColorAlpha(0,0);

      

      TLegend *ppPtLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      ppPtLegend->AddEntry(Normal->h_vn_pT_00to10_pp, "0 - 10%");
      ppPtLegend->AddEntry(Normal->h_vn_pT_10to40_pp, "10 - 40%");
      ppPtLegend->AddEntry(Normal->h_vn_pT_40to60_pp, "40 - 60%");
      ppPtLegend->SetBorderSize(0);
      ppPtLegend->SetFillColorAlpha(0,0);

      TLegend *pmPtLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      pmPtLegend->AddEntry(Normal->h_vn_pT_00to10_pm, "0 - 10%");
      pmPtLegend->AddEntry(Normal->h_vn_pT_10to40_pm, "10 - 40%");
      pmPtLegend->AddEntry(Normal->h_vn_pT_40to60_pm, "40 - 60%");
      pmPtLegend->SetBorderSize(0);
      pmPtLegend->SetFillColorAlpha(0,0);

      TLegend *kpPtLegend = new TLegend(0.18, 0.72, 0.38, 0.87);
      kpPtLegend->AddEntry(Normal->h_vn_pT_00to10_kp, "0 - 10%");
      kpPtLegend->AddEntry(Normal->h_vn_pT_10to40_kp, "10 - 40%");
      kpPtLegend->AddEntry(Normal->h_vn_pT_40to60_kp, "40 - 60%");
      kpPtLegend->SetBorderSize(0);
      kpPtLegend->SetFillColorAlpha(0,0);

      TLegend *kmPtLegend = new TLegend(0.19, 0.12, 0.39, 0.27);
      //kmPtLegend->AddEntry(Normal->h_vn_pT_00to10_km, "0 - 10%");
      kmPtLegend->AddEntry(Normal->h_vn_pT_10to40_km, "10 - 40%");
      //kmPtLegend->AddEntry(Normal->h_vn_pT_40to60_km, "40 - 60%");
      kmPtLegend->SetBorderSize(0);
      kmPtLegend->SetFillColorAlpha(0,0);

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


      TPaveText *ppText = new TPaveText(0.5, 0.07, 0.8, 0.14, "NB");
      ppText->AddText("#pi^{+}");
      ppText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      ppText->AddText("0.18 < p_{T} < 1.6 GeV");
      ppText->SetFillColorAlpha(0,0);
      ppText->SetLineColorAlpha(0,0);
      ppText->SetTextSize(.04);

      TPaveText *pmText = new TPaveText(0.5, 0.07, 0.8, 0.14, "NB");
      pmText->AddText("#pi^{-}");
      pmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pmText->AddText("0.18 < p_{T} < 1.6 GeV");
      pmText->SetFillColorAlpha(0,0);
      pmText->SetLineColorAlpha(0,0);
      pmText->SetTextSize(.04);
      
      TPaveText *kpText = new TPaveText(0.5, 0.07, 0.8, 0.14, "NB");
      kpText->AddText("K^{+}");
      kpText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kpText->AddText("0.4 < p_{T} < 1.6 GeV");
      kpText->SetFillColorAlpha(0,0);
      kpText->SetLineColorAlpha(0,0);
      kpText->SetTextSize(.04);
      
      TPaveText *kmText = new TPaveText(0.3, 0.05, 0.7, 0.12, "NB");
      kmText->AddText("K^{-}");
      kmText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kmText->AddText("0.4 < p_{T} < 1.6 GeV");
      kmText->SetFillColorAlpha(0,0);
      kmText->SetLineColorAlpha(0,0);
      kmText->SetTextSize(.04);
      
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
      


      TPaveText *ppPtText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      ppPtText->AddText("#pi^{+}");
      ppPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      ppPtText->AddText("0 < y_{CM} < 0.5");
      ppPtText->SetFillColorAlpha(0,0);
      ppPtText->SetLineColorAlpha(0,0);
      ppPtText->SetTextSize(.04);

      TPaveText *pmPtText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      pmPtText->AddText("#pi^{-}");
      pmPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      pmPtText->AddText("0 < y_{CM} < 0.5");
      pmPtText->SetFillColorAlpha(0,0);
      pmPtText->SetLineColorAlpha(0,0);
      pmPtText->SetTextSize(.04);

      TPaveText *kpPtText = new TPaveText(0.2, -0.22, 1.2, -0.1, "NB");
      kpPtText->AddText("K^{+}");
      kpPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kpPtText->AddText("0 < y_{CM} < 0.5");
      kpPtText->SetFillColorAlpha(0,0);
      kpPtText->SetLineColorAlpha(0,0);
      kpPtText->SetTextSize(.04);

      TPaveText *kmPtText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      kmPtText->AddText("K^{-}");
      kmPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      kmPtText->AddText("0 < y_{CM} < 0.5");
      kmPtText->SetFillColorAlpha(0,0);
      kmPtText->SetLineColorAlpha(0,0);
      kmPtText->SetTextSize(.04);

      TPaveText *prPtText = new TPaveText(0.2, 0.07, 1.2, 0.18, "NB");
      prPtText->AddText("Proton");
      prPtText->AddText("Au+Au #sqrt{s_{NN}} = 3.0 GeV FXT");
      prPtText->AddText("0 < y_{CM} < 0.5");
      prPtText->SetFillColorAlpha(0,0);
      prPtText->SetLineColorAlpha(0,0);
      prPtText->SetTextSize(.04);

      newFile->cd();
      
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
  
      //=== pi+- vs centrality
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pp->Clone());
      for (int i = 0; i < v_sys_pp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pp.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pm->Clone());
      for (int i = 0; i < v_sys_pm.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pm.at(i)); }

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
      for (int i = 0; i < v_sys_kp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_kp.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_km->Clone());
      for (int i = 0; i < v_sys_km.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_km.at(i)); }

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
      for (int i = 0; i < v_sys_pr.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pr.at(i)); }

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

/*
      //=== Proton Alternate Acceptance vs centrality
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pr_alt->Clone());
      for (int i = 0; i < v_sys_pr_alt.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pr_alt.at(i)); }

      copyWithNewErrors1->Write();
      
      Normal->h_vn_pr_alt->Draw();
      Normal->h_vn_pr_alt->SetMinimum(centralityLowerBound);
      Normal->h_vn_pr_alt->SetMaximum(centralityUpperBound);
      Normal->h_vn_pr_alt->Draw("E1");
      copyWithNewErrors1->Draw("[]");
      zeroLine->Draw("SAME");
      //prText->Draw();
      canvas->SaveAs("sys_h_vn_pr_alt.png");
      delete copyWithNewErrors1;
      canvas->Clear();
      //===


      //=== Deuteron vs centrality
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_de->Clone());
      for (int i = 0; i < v_sys_de.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_de.at(i)); }

      copyWithNewErrors1->Write();
      
      Normal->h_vn_de->Draw();
      Normal->h_vn_de->SetMinimum(centralityLowerBound);
      Normal->h_vn_de->SetMaximum(centralityUpperBound);
      Normal->h_vn_de->Draw("E1");
      copyWithNewErrors1->Draw("[]");
      zeroLine->Draw("SAME");
      //deText->Draw();
      canvas->SaveAs("sys_h_vn_de.png");
      delete copyWithNewErrors1;
      canvas->Clear();
      //===


      //=== Triton vs centrality
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_tr->Clone());
      for (int i = 0; i < v_sys_tr.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_tr.at(i)); }

      copyWithNewErrors1->Write();
      
      Normal->h_vn_tr->Draw();
      Normal->h_vn_tr->SetMinimum(centralityLowerBound);
      Normal->h_vn_tr->SetMaximum(centralityUpperBound);
      Normal->h_vn_tr->Draw("E1");
      copyWithNewErrors1->Draw("[]");
      zeroLine->Draw("SAME");
      //deText->Draw();
      canvas->SaveAs("sys_h_vn_tr.png");
      delete copyWithNewErrors1;
      canvas->Clear();
      //===
*/

      //=== Pi+ vs rapidity
      THStack *ppRapidityStack = new THStack("ppRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
      ppRapidityStack->Add(Normal->h_vn_yCM_00to10_pp);
      ppRapidityStack->Add(Normal->h_vn_yCM_10to40_pp);
      ppRapidityStack->Add(Normal->h_vn_yCM_40to60_pp);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_pp->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_pp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_pp.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_pp->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_pp.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_pp.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_pp->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_pp.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_pp.at(i)); }

      copyWithNewErrors1->Write();
      copyWithNewErrors2->Write();
      copyWithNewErrors3->Write();
      
      ppRapidityStack->Draw();
      ppRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      ppRapidityStack->GetXaxis()->SetNdivisions(210);
      ppRapidityStack->SetMaximum(rapidityUpperBound);
      ppRapidityStack->SetMinimum(rapidityLowerBound);
      ppRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      ppLegend->Draw();
      ppText->Draw();
      canvas->SaveAs("sys_ppRapidityStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===


      //=== Pi- vs rapidity
      THStack *pmRapidityStack = new THStack("pmRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
      pmRapidityStack->Add(Normal->h_vn_yCM_00to10_pm);
      pmRapidityStack->Add(Normal->h_vn_yCM_10to40_pm);
      pmRapidityStack->Add(Normal->h_vn_yCM_40to60_pm);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_pm->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_pm.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_pm.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_pm->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_pm.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_pm.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_pm->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_pm.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_pm.at(i)); }

      copyWithNewErrors1->Write();
      copyWithNewErrors2->Write();
      copyWithNewErrors3->Write();
      
      pmRapidityStack->Draw();
      pmRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      pmRapidityStack->GetXaxis()->SetNdivisions(210);
      pmRapidityStack->SetMaximum(rapidityUpperBound);
      pmRapidityStack->SetMinimum(rapidityLowerBound);
      pmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      pmLegend->Draw();
      pmText->Draw();
      canvas->SaveAs("sys_pmRapidityStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===

      //=== K+ vs rapidity
      THStack *kpRapidityStack = new THStack("kpRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
      kpRapidityStack->Add(Normal->h_vn_yCM_00to10_kp);
      kpRapidityStack->Add(Normal->h_vn_yCM_10to40_kp);
      kpRapidityStack->Add(Normal->h_vn_yCM_40to60_kp);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_kp->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_kp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_kp.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_kp->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_kp.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_kp.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_kp->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_kp.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_kp.at(i)); }

      copyWithNewErrors1->Write();
      copyWithNewErrors2->Write();
      copyWithNewErrors3->Write();
      
      kpRapidityStack->Draw();
      kpRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      kpRapidityStack->GetXaxis()->SetNdivisions(210);
      kpRapidityStack->SetMaximum(rapidityUpperBound);
      kpRapidityStack->SetMinimum(rapidityLowerBound);
      kpRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      kpLegend->Draw();
      kpText->Draw();
      canvas->SaveAs("sys_kpRapidityStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===



      //=== K- vs rapidity
      THStack *kmRapidityStack = new THStack("kmRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
      //kmRapidityStack->Add(Normal->h_vn_yCM_00to10_km);
      kmRapidityStack->Add(Normal->h_vn_yCM_10to40_km);
      //kmRapidityStack->Add(Normal->h_vn_yCM_40to60_km);

      /*
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_km->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_km.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_km.at(i)); }
      */

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_km->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_km.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_km.at(i)); }

      /*
      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_km->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_km.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_km.at(i)); }
      */

      //copyWithNewErrors1->Write();
      copyWithNewErrors2->Write();
      //copyWithNewErrors3->Write();

      kmRapidityStack->Draw();
      kmRapidityStack->GetYaxis()->SetTitleOffset(1.9);
      kmRapidityStack->GetXaxis()->SetNdivisions(210);
      kmRapidityStack->SetMaximum(rapidityUpperBound);
      kmRapidityStack->SetMinimum(rapidityLowerBound);
      kmRapidityStack->Draw("NOSTACK E1P");
      zeroLine_y->Draw("SAME");
      //copyWithNewErrors1->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      //copyWithNewErrors3->Draw("[]");
      kmLegend->Draw();
      kmText->Draw();
      canvas->SaveAs("sys_kmRapidityStack.png");
      canvas->Clear();
      //delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      //delete copyWithNewErrors3;
      //===


      //=== Proton vs rapidity
      THStack *prRapidityStack = new THStack("prRapidityStack", ";y-y_{mid};v_{"+order_n_str+"}");
      prRapidityStack->Add(Normal->h_vn_yCM_00to10_pr);
      prRapidityStack->Add(Normal->h_vn_yCM_10to40_pr);
      prRapidityStack->Add(Normal->h_vn_yCM_40to60_pr);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_00to10_pr->Clone());
      for (int i = 0; i < v_sys_yCM_00to10_pr.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_pr.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_pr->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_pr.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_pr.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_pr->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_pr.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_pr.at(i)); }

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
      for (int i = 0; i < v_sys_yCM_00to10_pr_symm.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_00to10_pr_symm.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_10to40_pr_symm->Clone());
      for (int i = 0; i < v_sys_yCM_10to40_pr_symm.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_10to40_pr_symm.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_40to60_pr_symm->Clone());
      for (int i = 0; i < v_sys_yCM_40to60_pr_symm.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_40to60_pr_symm.at(i)); }

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




      //=== Pi+ vs pT
      THStack *ppPtStack = new THStack("ppPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
      ppPtStack->Add(Normal->h_vn_pT_00to10_pp);
      ppPtStack->Add(Normal->h_vn_pT_40to60_pp);
      ppPtStack->Add(Normal->h_vn_pT_10to40_pp);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_pp->Clone());
      for (int i = 0; i < v_sys_pT_00to10_pp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pT_00to10_pp.at(i));}

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_pp->Clone());
      for (int i = 0; i < v_sys_pT_10to40_pp.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pT_10to40_pp.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_pp->Clone());
      for (int i = 0; i < v_sys_pT_40to60_pp.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_pT_40to60_pp.at(i)); }

      copyWithNewErrors1->Write();
      copyWithNewErrors2->Write();
      copyWithNewErrors3->Write();
      
      ppPtStack->Draw();
      ppPtStack->GetYaxis()->SetTitleOffset(1.9);
      ppPtStack->GetXaxis()->SetNdivisions(210);
      ppPtStack->SetMaximum(ptUpperBound);
      ppPtStack->SetMinimum(ptLowerBound);
      ppPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      ppPtLegend->Draw();
      ppPtText->Draw();
      canvas->SaveAs("sys_ppPtStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===



      //=== Pi- vs pT
      THStack *pmPtStack = new THStack("pmPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
      pmPtStack->Add(Normal->h_vn_pT_00to10_pm);
      pmPtStack->Add(Normal->h_vn_pT_40to60_pm);
      pmPtStack->Add(Normal->h_vn_pT_10to40_pm);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_pm->Clone());
      for (int i = 0; i < v_sys_pT_00to10_pm.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pT_00to10_pm.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_pm->Clone());
      for (int i = 0; i < v_sys_pT_10to40_pm.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pT_10to40_pm.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_pm->Clone());
      for (int i = 0; i < v_sys_pT_40to60_pm.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_pT_40to60_pm.at(i)); }

      copyWithNewErrors1->Write();
      copyWithNewErrors2->Write();
      copyWithNewErrors3->Write();

      pmPtStack->Draw();
      pmPtStack->GetYaxis()->SetTitleOffset(1.9);
      pmPtStack->GetXaxis()->SetNdivisions(210);
      pmPtStack->SetMaximum(ptUpperBound);
      pmPtStack->SetMinimum(ptLowerBound);
      pmPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      pmPtLegend->Draw();
      pmPtText->Draw();
      canvas->SaveAs("sys_pmPtStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===


      //=== K+ vs pT
      THStack *kpPtStack = new THStack("kpPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
      kpPtStack->Add(Normal->h_vn_pT_00to10_kp);
      kpPtStack->Add(Normal->h_vn_pT_40to60_kp);
      kpPtStack->Add(Normal->h_vn_pT_10to40_kp);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_kp->Clone());
      for (int i = 0; i < v_sys_pT_00to10_kp.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pT_00to10_kp.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_kp->Clone());
      for (int i = 0; i < v_sys_pT_10to40_kp.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pT_10to40_kp.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_kp->Clone());
      for (int i = 0; i < v_sys_pT_40to60_kp.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_pT_40to60_kp.at(i)); }

      copyWithNewErrors1->Write();
      copyWithNewErrors2->Write();
      copyWithNewErrors3->Write();

      kpPtStack->Draw();
      kpPtStack->GetYaxis()->SetTitleOffset(1.9);
      kpPtStack->GetXaxis()->SetNdivisions(210);
      kpPtStack->SetMaximum(ptUpperBound);
      kpPtStack->SetMinimum(ptLowerBound);
      kpPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      copyWithNewErrors1->Draw("[]");
      copyWithNewErrors3->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      kpPtLegend->Draw();
      kpPtText->Draw();
      canvas->SaveAs("sys_kpPtStack.png");
      canvas->Clear();
      delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      delete copyWithNewErrors3;
      //===


      //=== K- vs pT

      THStack *kmPtStack = new THStack("kmPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
      //kmPtStack->Add(Normal->h_vn_pT_00to10_km);
      //kmPtStack->Add(Normal->h_vn_pT_40to60_km);
      kmPtStack->Add(Normal->h_vn_pT_10to40_km);

      /*
      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_km->Clone());
      for (int i = 0; i < v_sys_pT_00to10_km.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pT_00to10_km.at(i)); }
      */

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_km->Clone());
      for (int i = 0; i < v_sys_pT_10to40_km.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pT_10to40_km.at(i)); }

      /*
      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_km->Clone());
      for (int i = 0; i < v_sys_pT_40to60_km.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_pT_40to60_km.at(i)); }
      */

      //copyWithNewErrors1->Write();
      copyWithNewErrors2->Write();
      //copyWithNewErrors3->Write();

      kmPtStack->Draw();
      kmPtStack->GetYaxis()->SetTitleOffset(1.9);
      kmPtStack->GetXaxis()->SetNdivisions(210);
      kmPtStack->SetMaximum(ptUpperBound);
      kmPtStack->SetMinimum(ptLowerBound);
      kmPtStack->Draw("NOSTACK E1P");
      zeroLine_pt->Draw("SAME");
      //copyWithNewErrors1->Draw("[]");
      //copyWithNewErrors3->Draw("[]");
      copyWithNewErrors2->Draw("[]");
      kmPtLegend->Draw();
      kmPtText->Draw();
      canvas->SaveAs("sys_kmPtStack.png");
      canvas->Clear();
      //delete copyWithNewErrors1;
      delete copyWithNewErrors2;
      //delete copyWithNewErrors3;
      //===


      //=== Proton vs pT
      THStack *prPtStack = new THStack("prPtStack", ";p_{T} (GeV);v_{"+order_n_str+"}");
      prPtStack->Add(Normal->h_vn_pT_00to10_pr);
      prPtStack->Add(Normal->h_vn_pT_40to60_pr);
      prPtStack->Add(Normal->h_vn_pT_10to40_pr);

      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_pT_00to10_pr->Clone());
      for (int i = 0; i < v_sys_pT_00to10_pr.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_pT_00to10_pr.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_pT_10to40_pr->Clone());
      for (int i = 0; i < v_sys_pT_10to40_pr.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_pT_10to40_pr.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_pT_40to60_pr->Clone());
      for (int i = 0; i < v_sys_pT_40to60_pr.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_pT_40to60_pr.at(i)); }

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

      /*
      // Proton vs yCM scan through pT


      copyWithNewErrors1 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_pr_1->Clone());
      for (int i = 0; i < v_sys_yCM_pr_1.size(); i++)
	{ copyWithNewErrors1->SetPointError(i, 0.0, v_sys_yCM_pr_1.at(i)); }

      copyWithNewErrors2 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_pr_2->Clone());
      for (int i = 0; i < v_sys_yCM_pr_2.size(); i++)
	{ copyWithNewErrors2->SetPointError(i, 0.0, v_sys_yCM_pr_2.at(i)); }

      copyWithNewErrors3 = new TGraphErrors((TH1D*)Normal->h_vn_yCM_pr_3->Clone());
      for (int i = 0; i < v_sys_yCM_pr_3.size(); i++)
	{ copyWithNewErrors3->SetPointError(i, 0.0, v_sys_yCM_pr_3.at(i)); }

      copyWithNewErrors1->Write();
      copyWithNewErrors2->Write();
      copyWithNewErrors3->Write();

      //===
      */
      
      delete canvas;
    }// End if order_n_str == 3


  delete epd;
  delete nhits;
  delete nSigPi;
  delete nSigKa;
  delete nSigPr;
  delete rvtx;
  delete zvtx;
  delete dca;
  delete nhitsdEdx;
  delete nhitsratio;
  delete m2Pi;
  delete m2Ka;

  delete Normal;
  delete epd_high_30;
  //delete epd_low_30;
  delete epd_scaled;
  delete nSigPi_high_30;
  delete nSigPi_low_30;
  delete nSigKa_high_30;
  delete nSigKa_low_30;
  delete nSigPr_high_30;
  delete nSigPr_low_30;
  delete rvtx_high_30;
  delete rvtx_low_30;
  delete zvtx_high_30;
  delete zvtx_low_30;
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

  //delete epd_high_20;
  //delete epd_low_20;
  delete nSigPi_high_20;
  delete nSigPi_low_20;
  delete nSigKa_high_20;
  delete nSigKa_low_20;
  delete nSigPr_high_20;
  delete nSigPr_low_20;
  delete rvtx_high_20;
  delete rvtx_low_20;
  delete zvtx_high_20;
  delete zvtx_low_20;
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

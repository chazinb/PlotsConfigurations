# include "TFile.h"
# include "TH1F.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>



//-----------------------------------------
// CheckYields
//-----------------------------------------


TH1F* CheckYields (TH1F* histo)
{     
      int nBins = histo -> GetNbinsX();
      std::cout<<"nbins "<< nBins <<std::endl;
      float histoIntegral = histo -> Integral();
      std::cout<<"Integral " << histoIntegral<<std::endl;
      for (int n = 1; n <= nBins; n++ )
      {
        float yValue = histo -> GetBinContent(n);
        if (yValue <= 0.) histo -> SetBinContent (n, 0.001);
      }
      float correctedIntegral = histo -> Integral();
      if (correctedIntegral != histoIntegral && histoIntegral!= 0) histo -> Scale (correctedIntegral/histoIntegral);
      std::cout << "corrected Integral" << correctedIntegral << std::endl; 
      return histo;

}


//-------------------------------------------------------
// Normalize the Histogram to Unity 
//-------------------------------------------------------


void MakeShape (TH1F *histoName){
    float integral;
    integral = histoName -> Integral(0, histoName -> GetNbinsX()+1);
    cout  << integral << endl; 
    histoName -> Scale(1/integral);
}


//---------------------------------------------------------------------
// Create_HistoTest
//--------------------------------------------------------------------- 


void Smearing_Shapes (TString inputDir = "Datacards/ValidationRegions/")
  
    /*
    *    This function does a "smearing up" of a set of .root histograms. 
    *    It returns a rootfile containing these smeared histos normalized 
    *    to unity.
    */ 

 {   std::cout << "init" << std::endl;
     // set the variables  
     const int nsample = 3;  TString sampleName [nsample] = {"06_WW", "05_ST", "04_TTTo2L2Nu"};
     const int ncut   = 2;  TString cutName [ncut]    = {"VR1_Tag", "VR1_NoTag"};
     enum { ee, mm, em, nchannel}; const TString schannel[nchannel] = {"ee", "mm","em"}; 
   
     // make tho output folder
     gSystem->mkdir("VRTesting/", kTRUE);

     std::cout << "start the loops" << std::endl;

     /* start the loop */
     
     // for each sample
     for (int i = 0; i < nsample; i ++)
     {
       std::cout << sampleName[i] << std::endl;
       // set the output rootfile name. Then recreate it  
       TString outFileName = ""; 
       if (sampleName[i] == "06_WW")          outFileName = "06_WW_smearing.root";
       if (sampleName[i] == "05_ST" )         outFileName = "05_ST_smearing.root";
       if (sampleName[i] == "04_TTTo2L2Nu")   outFileName = "04_TTTo2L2Nu_smearing.root";

       TFile* fileOut =  new TFile ("VRTesting/" + outFileName, "recreate"); 
       
       // for each cut
       for (int j = 0; j < ncut; j ++)
       {       
        std::cout << cutName[j] << std::endl;
        // for each channel
        for ( int k = 0; k < 3 ; k ++)
        {
         std::cout << schannel[k] << std::endl;
         // open the input file
         TString inputFile = inputDir + cutName[j] + "_" + schannel[k] + "/MT2ll/shapes/" + "histos_" +  cutName[j] + "_" + schannel[k] + ".root";
         TFile*  fileIn    = new TFile(inputFile, "read");
//         TH1F* histoOut = HistoTest ( fileIn, sampleName[i], cutName[j], schannel[k]);
         // take the histo
         TH1F*   histoIn   = (TH1F*) fileIn -> Get("histo_" + sampleName[i]);  
         
         // check the bad bins (Bin content < 0 )
         TH1F* histo = CheckYields (histoIn);

         // do the smearing up
         histo->Scale(1.);
         histo->SetBinContent(1, histo->GetBinContent(1)*1.000);
         histo->SetBinContent(2, histo->GetBinContent(2)*1.017);
         histo->SetBinContent(3, histo->GetBinContent(3)*1.032);
         histo->SetBinContent(4, histo->GetBinContent(4)*1.050);
         histo->SetBinContent(5, histo->GetBinContent(5)*1.067);
         histo->SetBinContent(6, histo->GetBinContent(6)*1.082);
         histo->SetBinContent(7, histo->GetBinContent(7)*1.100);
        
         // normalize histo to unity
         MakeShape(histo);
         
         // just clone to be safe
         TH1F* histoOut = (TH1F*)  histo -> Clone();
        
         // write the produced histo 
         fileOut->cd();
         histoOut -> SetName("h_" + sampleName[i] + "_" + schannel[k] + "_" + cutName[j]);
         histoOut -> Write();
         // close the input file
         fileIn -> cd();
         fileIn -> Close();
        }
       }      
       // close the input file
       fileOut->cd();
       fileOut -> Close();
       std::cout << "FIN" << std::endl;
      } 
}      
 

//--------------------------------------------------------------
// Compare_Shapes
//-------------------------------------------------------------- 

void Compare_Shapes()
{
 // open the input file1
 

 // open the input file2
 

 /* start the loop */

   // take histo1

   // take histo2

   // normalize histo2 to unity
 
} 
   


 

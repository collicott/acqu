// SVN Info: $Id: TCCalibEnergy.cxx 1038 2011-11-14 13:01:17Z werthm $

/*************************************************************************
 * Author: Irakli Keshelashvili, Dominik Werthmueller
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCCalibEnergy                                                        //
//                                                                      //
// Base energy calibration module class.                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TCCalibEnergy.h"

ClassImp(TCCalibEnergy)


//______________________________________________________________________________
TCCalibEnergy::TCCalibEnergy(const Char_t* name, const Char_t* title, const Char_t* data,
                             Int_t nElem)
    : TCCalib(name, title, data, nElem)
{
    // Constructor.
    
    // init members
    fPi0Pos = 0;
    fLine = 0;

}

//______________________________________________________________________________
TCCalibEnergy::~TCCalibEnergy()
{
    // Destructor. 
    
    if (fLine) delete fLine;
}

//______________________________________________________________________________
void TCCalibEnergy::Init()
{
    // Init the module.
    
    Char_t tmp[256];

    // init members
    fPi0Pos = 0;
    fLine = new TLine();
    
    // configure line
    fLine->SetLineColor(4);
    fLine->SetLineWidth(3);
 
    // get histogram name
    sprintf(tmp, "%s.Histo.Fit.Name", GetName());
    if (!TCReadConfig::GetReader()->GetConfig(tmp))
    {
        Error("Init", "Histogram name was not found in configuration!");
        return;
    }
    else fHistoName = *TCReadConfig::GetReader()->GetConfig(tmp);
    
    // read old parameters (only from first set)
    TCMySQLManager::GetManager()->ReadParameters(fData, fCalibration.Data(), fSet[0], fOldVal, fNelem);
    
    // copy to new parameters
    for (Int_t i = 0; i < fNelem; i++) fNewVal[i] = fOldVal[i];

    // sum up all files contained in this runset
    TCFileManager f(fData, fCalibration.Data(), fNset, fSet);
    
    // get the main calibration histogram
    fMainHisto = f.GetHistogram(fHistoName.Data());
    if (!fMainHisto)
    {
        Error("Init", "Main histogram does not exist!\n");
        return;
    }
    
    // create the overview histogram
    fOverviewHisto = new TH1F("Overview", ";Element;2#gamma inv. mass [MeV]", fNelem, 0, fNelem);
    fOverviewHisto->SetMarkerStyle(2);
    fOverviewHisto->SetMarkerColor(4);
    
    // draw main histogram
    fCanvasFit->Divide(1, 2, 0.001, 0.001);
    fCanvasFit->cd(1)->SetLogz();
    sprintf(tmp, "%s.Histo.Fit", GetName());
    TCUtils::FormatHistogram(fMainHisto, tmp);
    fMainHisto->Draw("colz");

    // draw the overview histogram
    fCanvasResult->cd();
    sprintf(tmp, "%s.Histo.Overview", GetName());
    TCUtils::FormatHistogram(fOverviewHisto, tmp);
    fOverviewHisto->Draw("P");
}

//______________________________________________________________________________
void TCCalibEnergy::Fit(Int_t elem)
{
    // Perform the fit of the element 'elem'.
    
    Char_t tmp[256];
    Char_t tmp2[256]; //Maria
    Char_t tmp3[256]; //Maria
    // create histogram projection for this element
    sprintf(tmp, "ProjHisto_%i", elem);
    TH2* h2 = (TH2*) fMainHisto;
    if (fFitHisto) delete fFitHisto;
    fFitHisto = (TH1D*) h2->ProjectionX(tmp, elem+1, elem+1, "e");
    
    // draw histogram
    fFitHisto->SetFillColor(35);
    fCanvasFit->cd(2);
    sprintf(tmp, "%s.Histo.Fit", GetName());
    TCUtils::FormatHistogram(fFitHisto, tmp);
    fFitHisto->Draw("hist");
     
    // check for sufficient statistics
    if (fFitHisto->GetEntries() > 500)
    {
        // delete old function
        /*if (fFitFunc) delete fFitFunc;
        sprintf(tmp, "fEnergy_%i", elem);
        //fFitFunc = new TF1(tmp, "gaus(0)+pol3(3)");
        fFitFunc = new TF1(tmp, "gaus(0)+[3]/(1+exp((x-[4])/[5]))");
        //fFitFunc = new TF1(tmp, "gaus(0)+expo(3)");  //Maria
        fFitFunc->SetLineColor(2);
        
        // estimate peak position
        fPi0Pos = fFitHisto->GetBinCenter(fFitHisto->GetMaximumBin());
        if (fPi0Pos < 100 || fPi0Pos > 160) fPi0Pos = 135;
		*/
        // configure fitting function
        if (this->InheritsFrom("TCCalibCBEnergy"))
        {
			// delete old function
			if (fFitFunc) delete fFitFunc;
			sprintf(tmp, "fEnergy_%i", elem);
			//fFitFunc = new TF1(tmp, "gaus(0)+pol3(3)");
			fFitFunc = new TF1(tmp, "gaus(0)+[3]/(1+exp((x-[4])/[5]))");
			//fFitFunc = new TF1(tmp, "gaus(0)+expo(3)");  //Maria
			fFitFunc->SetLineColor(2);
        
			// estimate peak position
			fPi0Pos = fFitHisto->GetBinCenter(fFitHisto->GetMaximumBin());
			if (fPi0Pos < 100 || fPi0Pos > 160) fPi0Pos = 135;
			
            Double_t Threshold;
            Double_t Low;
            Double_t High;
            Double_t Mampl;

            Mampl = fFitHisto->GetBinContent(95);
            Threshold = 0.1*fFitHisto->GetMaximum();
            Low = fFitHisto->GetBinCenter(fFitHisto->FindFirstBinAbove(Threshold,1));
            High = fFitHisto->GetBinCenter(fFitHisto->FindLastBinAbove(Threshold,1));

            fFitFunc->SetRange(95, 200);
            //fFitFunc->SetParameters(fFitHisto->GetMaximum(), fPi0Pos, 8, 1, 1, 1, 0.1);
            fFitFunc->SetParameters(fFitHisto->GetMaximum(), fPi0Pos, 8, 1, fPi0Pos,High-Low);  //Maria
            fFitFunc->SetParLimits(0, 0.1, fFitHisto->GetMaximum());
            fFitFunc->SetParLimits(1, 130, 140);  
            fFitFunc->SetParLimits(2, 3, 15); 
            fFitFunc->SetParLimits(3,1., Mampl);
            fFitFunc->SetParLimits(4,fPi0Pos-20,fPi0Pos+20);
            fFitFunc->SetParLimits(5,5,70);

            for (Int_t i = 0; i < 10; i++)
            if (!fFitHisto->Fit(fFitFunc, "RBQ0")) break;
        }
        else if (this->InheritsFrom("TCCalibTAPSEnergyLG"))
        {
			// delete old function
			if (fFitFunc) delete fFitFunc;
			sprintf(tmp, "fEnergy_%i", elem);
			//fFitFunc = new TF1(tmp, "gaus(0)+pol3(3)");
			fFitFunc = new TF1(tmp, "gaus(0)+exp([3]+x*[4])");  //Maria
			fFitFunc->SetLineColor(2);
        
			// estimate peak position
			fPi0Pos = fFitHisto->GetBinCenter(fFitHisto->GetMaximumBin());
			if (fPi0Pos < 100 || fPi0Pos > 160) fPi0Pos = 135;
        
			//fFitFunc->SetRange(80, 200);
			fFitFunc->SetRange(40, 240);
			
			fFitFunc->FixParameter(0,0);
			fFitFunc->FixParameter(1,0);
			fFitFunc->FixParameter(2,0);
			
			for (Int_t i = 0; i < 10; i++)
            if (!fFitHisto->Fit(fFitFunc, "RBQ0L")) break;
            
			fFitFunc->SetParLimits(0, 1, 2000);
			fFitFunc->SetParLimits(1, 120, 140);
			fFitFunc->SetParLimits(2, 3, 15);
			
			Double_t expo1;
			Double_t expo2;
			expo1 = fFitFunc->GetParameter(3);
			expo2 = fFitFunc->GetParameter(4);
			if(expo1>=0) fFitFunc->SetParLimits(3,0.8*expo1,1.2*expo1);
			else fFitFunc->SetParLimits(3,1.2*expo1,0.8*expo1);
			
			if(expo2>=0) fFitFunc->SetParLimits(4,0.8*expo2,1.2*expo2);
			else fFitFunc->SetParLimits(4,1.2*expo2,0.8*expo2);
			
			for (Int_t i = 0; i < 10; i++)
            if (!fFitHisto->Fit(fFitFunc, "RBQ0L")) break;
			
			
        }
/*
        // fit
        for (Int_t i = 0; i < 10; i++)
            if (!fFitHisto->Fit(fFitFunc, "RBQ0L")) break;
*/    
        // final results
        fPi0Pos = fFitFunc->GetParameter(1);         
        
        /*
        Double_t fPi0PosM = fFitFunc->GetParameter(1);
        Double_t sigmaM = fFitFunc->GetParameter(2);
        Double_t ampM = fFitFunc->GetParameter(0);

		fFitFunc->SetNpx(100000) ;
		
		//fPi0Pos = fFitFunc->GetMaximumX(fPi0PosM-10,fPi0PosM+10,1.E-10,100,false);
		fPi0Pos = fFitFunc->GetMaximumX(fPi0PosM-0.8*sigmaM,fPi0PosM+0.8*sigmaM,1.E-10,100,false);
		
		if((ampM/sigmaM)<2)fPi0Pos=fPi0PosM;
		*/
		
        // check if mass is in normal range
        if (fPi0Pos < 80 || fPi0Pos > 200) fPi0Pos = 135;
 
        // set indicator line
        fLine->SetY1(0);
        fLine->SetY2(fFitHisto->GetMaximum() + 20);
        fLine->SetX1(fPi0Pos);
        fLine->SetX2(fPi0Pos);
   
        // draw fitting function
        if (fFitFunc) fFitFunc->Draw("same");
        
        
        //-------------------------------------------------------  Maria  ---------------------------------------------------------------------------
        Double_t Gamp;
        Double_t Gmean;
        Double_t Gsigma;
        if (fgausFunc) delete fgausFunc;
        sprintf(tmp2, "Gaus_fEnergy_%i", elem);
        fgausFunc = new TF1(tmp2, "gaus(0)");
        //fgausFunc->SetRange(110, 180);
        fgausFunc->SetRange(95, 200);
        Gamp = fFitFunc->GetParameter(0);
        Gmean = fFitFunc->GetParameter(1);
        Gsigma = fFitFunc->GetParameter(2);
        fgausFunc->SetParameters(Gamp,Gmean,Gsigma);
		fgausFunc->SetLineColor(12);
		
		if (fbackFunc) delete fbackFunc;
		sprintf(tmp3, "Background_fEnergy_%i", elem);
		//fbackFunc = new TF1(tmp3, "pol3(0)");
		if (this->InheritsFrom("TCCalibCBEnergy"))
		{
			fbackFunc = new TF1(tmp3, "[0]/(1+exp((x-[1])/[2]))");
			//fbackFunc->SetRange(110, 180);
			fbackFunc->SetRange(95, 200);
			fbackFunc->SetParameters(fFitFunc->GetParameter(3),fFitFunc->GetParameter(4),fFitFunc->GetParameter(5));
			fbackFunc->SetLineColor(1);
		}
		else if (this->InheritsFrom("TCCalibTAPSEnergyLG"))
		{
			fbackFunc = new TF1(tmp3, "exp([0]+x*[1])");
			//fbackFunc->SetRange(110, 180);
			fbackFunc->SetRange(60, 200);
			fbackFunc->SetParameters(fFitFunc->GetParameter(3),fFitFunc->GetParameter(4));
			fbackFunc->SetLineColor(1);
		}
		if (fgausFunc) fgausFunc->Draw("same");
		if (fbackFunc) fbackFunc->Draw("same");
		//----------------------------------------------------- End Maria ---------------------------------------------------------------------------
       
       
        // draw indicator line
        fLine->Draw();
    }
    
    // update canvas
    fCanvasFit->Update();

    // update overview
    if (elem % 20 == 0)
    {
        fCanvasResult->cd();
        fOverviewHisto->Draw("E1");
        fCanvasResult->Update();
    }   
}

void TCCalibEnergy::Calculate(Int_t elem)
//______________________________________________________________________________
{
    // Calculate the new value of the element 'elem'.
    
    Bool_t unchanged = kFALSE;

    // check if fit was performed
    if (fFitHisto->GetEntries() > 500)
    {
        // check if line position was modified by hand
        if (fLine->GetX1() != fPi0Pos) fPi0Pos = fLine->GetX1();
        
        // calculate the new offset
        fNewVal[elem] = fOldVal[elem] * TCConfig::kPi0Mass / fPi0Pos;
        
        // if new value is negative take old
        if (fNewVal[elem] < 0) 
        {
            fNewVal[elem] = fOldVal[elem];
            unchanged = kTRUE;
        }

        // update overview histogram
        fOverviewHisto->SetBinContent(elem+1, fPi0Pos);
        fOverviewHisto->SetBinError(elem+1, 0.0000001);
    
        // update average calculation
        fAvr += fPi0Pos;
        fAvrDiff += TMath::Abs(fPi0Pos - TCConfig::kPi0Mass);
        fNcalc++;
    }
    else
    {   
        // do not change old value
        fNewVal[elem] = fOldVal[elem];
        unchanged = kTRUE;
    }

    // user information
    printf("Element: %03d    Pi0: %12.8f    "
           "old gain: %12.8f    new gain: %12.8f    diff: %6.2f %%",
           elem, fPi0Pos, fOldVal[elem], fNewVal[elem],
           TCUtils::GetDiffPercent(fOldVal[elem], fNewVal[elem]));
    if (unchanged) printf("    -> unchanged");
    if (this->InheritsFrom("TCCalibCBEnergy"))
    {
        if (TCUtils::IsCBHole(elem)) printf(" (hole)");
    }
    printf("\n");

    // show average
    if (elem == fNelem-1)
    {
        fAvr /= (Double_t)fNcalc;
        fAvrDiff /= (Double_t)fNcalc;
        printf("Average pi0 position           : %.3f MeV\n", fAvr);
        printf("Average difference to pi0 mass : %.3f MeV\n", fAvrDiff);
    }
}   


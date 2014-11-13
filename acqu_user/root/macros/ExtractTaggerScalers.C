void main()
{
ExtractTaggerScalers(242);
ExtractTaggerScalers(243);
ExtractTaggerScalers(305);
ExtractTaggerScalers(306);
ExtractTaggerScalers(354);
ExtractTaggerScalers(355);
ExtractTaggerScalers(427);
ExtractTaggerScalers(428);
ExtractTaggerScalers(429);
ExtractTaggerScalers(430);
ExtractTaggerScalers(449);
ExtractTaggerScalers(483);
ExtractTaggerScalers(484);
ExtractTaggerScalers(485);
ExtractTaggerScalers(547);
ExtractTaggerScalers(548);
ExtractTaggerScalers(549);
ExtractTaggerScalers(603);
ExtractTaggerScalers(604);
ExtractTaggerScalers(605);
ExtractTaggerScalers(659);
ExtractTaggerScalers(660);
ExtractTaggerScalers(661);
ExtractTaggerScalers(720);
ExtractTaggerScalers(721);
ExtractTaggerScalers(722);
ExtractTaggerScalers(769);
ExtractTaggerScalers(770);
ExtractTaggerScalers(771);
ExtractTaggerScalers(772);
ExtractTaggerScalers(783);
ExtractTaggerScalers(786);
}

void ExtractTaggerScalers(Int_t runNumber)
{
	Char_t* file_in  = Form("Hist_Compton_%d.root",runNumber);
	Char_t* file_out = Form("Compton_%d.scalDump", runNumber);

	TFile *file = new TFile(file_in);
	TH1D  *FPD_ScalerAcc = (TH1D*)file->Get("FPD_ScalerAcc");

	FILE  *fp=fopen(file_out,"w");

	if(FPD_ScalerAcc)
	{
	   for(int n=1;n<=352;n++)
	   {
	      fprintf(fp,"%d %d\n",n, FPD_ScalerAcc->GetBinContent(n));
	   }
	}
	fclose(fp);
} 

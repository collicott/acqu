//--Author	Cristina Collicott   June 2013   Compton and Pi0 analysis

#include "TA2OnlinePhys.h"

enum {};
static const Map_t kInputs[] = {
	{NULL,          -1}
};

ClassImp(TA2OnlinePhys)

//-----------------------------------------------------------------------------
TA2OnlinePhys::TA2OnlinePhys( const char* name, TA2Analysis* analysis )
	:TA2Physics( name, analysis ) 
{
// Initialise Detectors
	fTAGG		= NULL; // Tagger
	fLadder		= NULL; // Tagger ladder
	fCB			= NULL; // CB (PID, MWPC, NaI)
	fTAPS		= NULL; // TAPS

	// Particle arrays
	fTaggedPhoton	= NULL;

	fNTagg			= 0;
	fTaggerChannel	= NULL;
	fTaggerTime		= NULL;
	
	fNParticle		= 0;
	fParticleP4		= NULL;
	fParticleApp	= NULL;

	AddCmdList(kInputs);
}


//-----------------------------------------------------------------------------
TA2OnlinePhys::~TA2OnlinePhys()
{
}
	
//-----------------------------------------------------------------------------
void TA2OnlinePhys::SetConfig(Char_t* line, Int_t key)
{
	// Any special command-line input for Crystal Ball apparatus

	switch (key){
		default:
		// default main apparatus SetConfig()
		TA2Physics::SetConfig( line, key );
		break;
	}
}

//---------------------------------------------------------------------------
void TA2OnlinePhys::PostInit()
{

// Introduce Detectors

	// Tagger
	fTAGG = (TA2Tagger*)((TA2Analysis*)fParent)->GetChild("TAGG");
	if ( !fTAGG) PrintError("","<No Tagger class found>",EErrFatal);
	else {  printf("Tagger included in analysis\n");
		fTAGGParticles = fTAGG->GetParticles(); }

	// Ladder
	fLadder = (TA2Ladder*)((TA2Analysis*)fParent)->GetGrandChild( "FPD");
	if (!fLadder) PrintError( "", "<No Ladder class found>", EErrFatal);

	// Central Apparatus
	fCB = (TA2CrystalBall*)((TA2Analysis*)fParent)->GetChild("CB");	
	if (!fCB) PrintError( "", "<No CrystalBall class found>", EErrFatal);
	else {  printf("CB system included in analysis\n"); }

	// TAPS
	fTAPS = (TA2Calorimeter*)((TA2Analysis*)fParent)->GetChild("TAPS");
	if ( !fTAPS) printf("TAPS *NOT* included in analysis\n");
	else {  printf("TAPS included in analysis\n"); }

	printf("\n");

// Get max # of Particles from detectors, used for defining array sizes

	fCBMaxNParticle = fCB->GetMaxParticle();	
	if (fTAPS)  	fTAPSMaxNParticle = fTAPS->GetMaxParticle(); 
	else 			fTAPSMaxNParticle = 0;
	fMaxNParticle = fCBMaxNParticle + fTAPSMaxNParticle;  

	// Create arrays to hold Particles
	fParticleP4		= new TLorentzVector[fMaxNParticle];
	fParticleApp	= new Int_t[fMaxNParticle];

	// Make 2x the room (enough) for tagger for multi-hit

	fTaggedPhoton	= new TA2Particle*[352*2];
	fTaggerTime		= new Double_t[352*2];
	fTaggerChannel	= new Int_t[352*2];

	// Define hard coded histograms
	DefineHistograms();

	// Default physics initialisation
	TA2Physics::PostInit();
}

//-----------------------------------------------------------------------------
void TA2OnlinePhys::LoadVariable( )
{

// Input name - variable pointer associations for any subsequent cut/histogram setup

	TA2Physics::LoadVariable();

	return;
}

//-----------------------------------------------------------------------------
void TA2OnlinePhys::Reconstruct() 
{
	
	ZeroCounters();
	
	GetCBParticles();
	GetTAPSParticles();
	GetTagger();
	
	BasicPhysCheck();

}

void TA2OnlinePhys::DefineHistograms()
{
	// Hard coded histograms
	gROOT->cd();
	
	IM_2g 		= new TH1D("PHYS_IM_2g", 		"IM of 2 photon events", 1000, 0, 1000);
	IM_2g_CB 	= new TH1D("PHYS_IM_2g_CB", 	"IM of 2 photon events - CB only", 1000, 0, 1000);
	IM_2g_TAPS 	= new TH1D("PHYS_IM_2g_TAPS", 	"IM of 2 photon events - TAPS only", 1000, 0, 1000);
	IM_2g_mix 	= new TH1D("PHYS_IM_2g_mix", 	"IM of 2 photon events - CB & TAPS", 1000, 0, 1000);

	IM_3g 		= new TH1D("PHYS_IM_3g", 		"IM of 3 photon events", 1000, 0, 1000);
	IM_3g_CB 	= new TH1D("PHYS_IM_3g_CB", 	"IM of 3 photon events - CB only", 1000, 0, 1000);
	IM_3g_E300	= new TH1D("PHYS_IM_3g_E300", 	"IM of 3 photon events - all E > 300 MeV", 1000, 0, 1000);
	
	IM_6g 		= new TH1D("PHYS_IM_6g", 		"IM of 6 photon events", 1000, 0, 1000);
	IM_6g_CB 	= new TH1D("PHYS_IM_6g_CB", 	"IM of 6 photon events - CB only", 1000, 0, 1000);

		
}

void TA2OnlinePhys::BasicPhysCheck()
{
	// 2 Gamma invariant mass	
	if(fNParticle == 2)
	{
		TLorentzVector p4 = fParticleP4[0] + fParticleP4[1]; 
		
		IM_2g->Fill(p4.M());
		
		if((fParticleApp[0] == 1) && (fParticleApp[1] == 1)) 
		IM_2g_CB->Fill(p4.M());
		
		else if((fParticleApp[0] == 2) && (fParticleApp[1] == 2)) 
		IM_2g_TAPS->Fill(p4.M());
		
		else 
		IM_2g_mix->Fill(p4.M());
	}

	// 3g case
	if(fNParticle == 3)
	{
		TLorentzVector p4 = fParticleP4[0] 
						  + fParticleP4[1]  
						  +	fParticleP4[2];
		
		IM_3g->Fill(p4.M());
		
		if((fParticleApp[0] == 1) 
		&& (fParticleApp[1] == 1) 
		&& (fParticleApp[2] == 1) )
		IM_3g_CB->Fill(p4.M());
		
		if((fParticleP4[0].T() >= 300.0) 
		&& (fParticleP4[1].T() >= 300.0) 
		&& (fParticleP4[2].T() >= 300.0) )
		IM_3g_E300->Fill(p4.M());		
		
	}		
	
	// 6g case
	if(fNParticle == 6)
	{
		TLorentzVector p4 = fParticleP4[0] + fParticleP4[1] 
						  +	fParticleP4[2] + fParticleP4[3]
						  +	fParticleP4[4] + fParticleP4[5];
		
		IM_6g->Fill(p4.M());
		
		if((fParticleApp[0] == 1) && (fParticleApp[1] == 1) 
		&& (fParticleApp[2] == 1) && (fParticleApp[3] == 1)
		&& (fParticleApp[4] == 1) && (fParticleApp[5] == 1)	)
		IM_6g_CB->Fill(p4.M());
		
	}	
}

void TA2OnlinePhys::ZeroCounters()
{
	fNParticle = 0;
}

void TA2OnlinePhys::GetCBParticles()
{
	fCBNParticle = fCB->GetNParticle();
	for (UInt_t i = 0; i < fCBNParticle; i++ ) {
		
		 TLorentzVector* p4cb = fCB->GetP4();
		 
		fParticleP4[fNParticle+i]  = p4cb[i];
		fParticleApp[fNParticle+i] = 1;
	}
	
	// update # of particles
	fNParticle+=fCBNParticle;
}

void TA2OnlinePhys::GetTAPSParticles()
{
	fTAPSNParticle = fTAPS->GetNParticle();
	for (UInt_t i = 0; i < fTAPSNParticle; i++ ) {
		
		 TLorentzVector* p4taps = fTAPS->GetP4();
		 
		fParticleP4[fNParticle+i]  = p4taps[i];
		fParticleApp[fNParticle+i] = 2;
	}
	
	// update # of particles
	fNParticle+=fTAPSNParticle;	
}

void TA2OnlinePhys::GetTagger()
{
	// Tagger
	if(fTAGG && fLadder)
	{
		// Collect Tagger M0 Hits
		fNTagg	= fLadder->GetNhits();
		for(UInt_t i=0; i<fNTagg; i++)
		{
			fTaggerChannel[i]	= fLadder->GetHits(i);
			fTaggerTime[i]		= (fLadder->GetTimeOR())[i];
			fTaggedPhoton[i] 	= fTAGGParticles+i;
		}
	
		// Collect Tagger M+ Hits
		for(UInt_t m=1; m<fLadder->GetNMultihit(); m++)
		{
			for(UInt_t i=0; i<fLadder->GetNhitsM(m); i++)
			{
				fTaggerChannel[fNTagg+i] 	= (fLadder->GetHitsM(m))[i];
				fTaggerTime[fNTagg+i]	 	= (fLadder->GetTimeORM(m))[i];
			}
			fNTagg	+= fLadder->GetNhitsM(m);
		}
	}
	else fNTagg = 0;	
	
}

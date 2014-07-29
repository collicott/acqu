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

// Pi0Compton Variables

	// Particle Counters
	fNPhoton		= 0;
	fNProton		= 0;
	fNUnknown		= 0;	

	// Particle arrays
	fTaggedPhoton	= NULL;
	fPhoton 		= NULL;	
	fProton			= NULL;
	fUnknown		= NULL;

	fPhotonEnergy	= NULL;
	fPhotonTheta	= NULL;
	fPhotonPhi		= NULL;
	fPhotonTime		= NULL;
	fPhotonApp		= NULL;

	fProtonEnergy	= NULL;
	fProtonTheta	= NULL;
	fProtonPhi		= NULL;
	fProtonTime		= NULL;
	fProtonApp		= NULL;

	fNTagg			= 0;
	fTaggerChannel	= NULL;
	fTaggerTime		= NULL;
	
// Histogram variables


	
	

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
	fCB = (TA2CentralApparatus*)((TA2Analysis*)fParent)->GetChild("CB");	
	if (!fCB) PrintError( "", "<No Central Apparatus/CB class found>", EErrFatal);
	else {  printf("CB system included in analysis\n");
		fCBParticles  = fCB->GetParticles(); }

	// TAPS
	fTAPS = (TA2Taps*)((TA2Analysis*)fParent)->GetChild("TAPS");
	if ( !fTAPS) printf("TAPS *NOT* included in analysis\n");
	else {  printf("TAPS included in analysis\n");
		fTAPSParticles = fTAPS->GetParticles();	}

	printf("\n");

// Get max # of Particles from detectors, used for defining array sizes

	fCBMaxNParticle = fCB->GetMaxParticle();	
	if (fTAPS)  	fTAPSMaxNParticle = fTAPS->GetMaxParticle(); 
	else 			fTAPSMaxNParticle = 0;
	fMaxNParticle = fCBMaxNParticle + fTAPSMaxNParticle;  

// Create arrays to hold Particles

	fPhoton 		= new TA2Particle*[fMaxNParticle];	
	fProton			= new TA2Particle*[fMaxNParticle];
	fUnknown		= new TA2Particle*[fMaxNParticle];
	fTaggedPhoton	= new TA2Particle*[352*2];

	fPhotonEnergy	= new Double_t[fMaxNParticle];
	fPhotonTheta	= new Double_t[fMaxNParticle];
	fPhotonPhi		= new Double_t[fMaxNParticle];
	fPhotonTime		= new Double_t[fMaxNParticle];
	fPhotonApp		= new Int_t[fMaxNParticle];

	fProtonEnergy	= new Double_t[fMaxNParticle];
	fProtonTheta	= new Double_t[fMaxNParticle];
	fProtonPhi		= new Double_t[fMaxNParticle];
	fProtonTime		= new Double_t[fMaxNParticle];
	fProtonApp		= new Int_t[fMaxNParticle];
	

	// Make 2x the room (enough) for tagger for multi-hit
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

	TA2DataManager::LoadVariable("NPhoton", 		&fNPhoton,			EISingleX);
	TA2DataManager::LoadVariable("PhotonTheta", 	fPhotonTheta,		EDMultiX);
	TA2DataManager::LoadVariable("PhotonPhi", 		fPhotonPhi,			EDMultiX);
	TA2DataManager::LoadVariable("PhotonEnergy", 	fPhotonEnergy,		EDMultiX);
	TA2DataManager::LoadVariable("PhotonTime", 		fPhotonTime,		EDMultiX);

	TA2DataManager::LoadVariable("NProton", 		&fNProton,			EISingleX);
	TA2DataManager::LoadVariable("ProtonTheta", 	fProtonTheta,		EDMultiX);
	TA2DataManager::LoadVariable("ProtonPhi", 		fProtonPhi,			EDMultiX);
	TA2DataManager::LoadVariable("ProtonEnergy", 	fProtonEnergy,		EDMultiX);
	TA2DataManager::LoadVariable("ProtonTime", 		fProtonTime,		EDMultiX);

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
	FillEndBuffers();
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
	if(fNPhoton == 2)
	{
		TLorentzVector p4 = fPhoton[0]->GetP4() + fPhoton[1]->GetP4(); 
		
		IM_2g->Fill(p4.M());
		
		if((fPhotonApp[0] == 1) && (fPhotonApp[1] == 1)) 
		IM_2g_CB->Fill(p4.M());
		
		else if((fPhotonApp[0] == 2) && (fPhotonApp[1] == 2)) 
		IM_2g_TAPS->Fill(p4.M());
		
		else 
		IM_2g_mix->Fill(p4.M());
	}

	// 3g case
	if(fNPhoton == 3)
	{
		TLorentzVector p4 = fPhoton[0]->GetP4()
						  + fPhoton[1]->GetP4() 
						  +	fPhoton[2]->GetP4();
		
		IM_3g->Fill(p4.M());
		
		if((fPhotonApp[0] == 1) 
		&& (fPhotonApp[1] == 1) 
		&& (fPhotonApp[2] == 1) )
		IM_3g_CB->Fill(p4.M());
		
		if((fPhotonEnergy[0] >= 300.0) 
		&& (fPhotonEnergy[1] >= 300.0) 
		&& (fPhotonEnergy[2] >= 300.0) )
		IM_3g_E300->Fill(p4.M());		
		
	}		
	
	// 6g case
	if(fNPhoton == 6)
	{
		TLorentzVector p4 = fPhoton[0]->GetP4() + fPhoton[1]->GetP4() 
						  +	fPhoton[2]->GetP4() + fPhoton[3]->GetP4()
						  +	fPhoton[4]->GetP4() + fPhoton[5]->GetP4();
		
		IM_6g->Fill(p4.M());
		
		if((fPhotonApp[0] == 1) && (fPhotonApp[1] == 1) 
		&& (fPhotonApp[2] == 1) && (fPhotonApp[3] == 1)
		&& (fPhotonApp[4] == 1) && (fPhotonApp[5] == 1)	)
		IM_6g_CB->Fill(p4.M());
		
	}	
}

void TA2OnlinePhys::ZeroCounters()
{
	fNPhoton = 0;
	fNProton = 0;
	fNUnknown = 0;	
}

void TA2OnlinePhys::GetCBParticles()
{
	// CB
	fCBNParticle = fCB->GetNParticle();
	for ( i = 0; i < fCBNParticle; i++ ) {

		switch( (fCBParticles+i)->GetParticleID() ) { // Get PDG code

		case kGamma:                               	// Identified as a Photon
		fPhoton[fNPhoton] 	= fCBParticles+i;       // Add to Photon Array
		fPhotonEnergy[fNPhoton]	= fPhoton[fNPhoton]->GetT();
		fPhotonTheta[fNPhoton]	= fPhoton[fNPhoton]->GetThetaDg();
		fPhotonPhi[fNPhoton]	= fPhoton[fNPhoton]->GetPhiDg();
		fPhotonTime[fNPhoton]	= fPhoton[fNPhoton]->GetTime();	
		fPhotonApp[fNPhoton]	= 1;   

		fNPhoton++;					// Increase Photon counter
		break;

		case kProton:                               	// Identified as a Proton
		fProton[fNProton]	= fCBParticles+i;       // Add to Proton Array
		fProtonEnergy[fNPhoton]	= fProton[fNProton]->GetT();
		fProtonTheta[fNPhoton]	= fProton[fNProton]->GetThetaDg();
		fProtonPhi[fNPhoton]	= fProton[fNProton]->GetPhiDg();
		fProtonTime[fNPhoton]	= fProton[fNProton]->GetTime();	
		fProtonApp[fNPhoton]	= 1;   
		 
		fNProton++;					// Increase Proton counter
		break;

		default:
		fUnknown[fNUnknown]   	= fCBParticles+i;    	// include in "Unknown" list
		fNUnknown++; 					// Increase "Unknown" counter

		}						// End switch/case function
	}							// End # of Particle Loop
}

void TA2OnlinePhys::GetTAPSParticles()
{
	// TAPS
	if(fTAPS) 
	{
	fTAPSNParticle	= fTAPS->GetNparticle(); 
	for ( i = 0; i < fTAPSNParticle; i++ ) {

		switch( (fTAPSParticles+i)->GetParticleID() ) { // Get PDG code

		case kGamma:                               	// Identified as a Photon
	        fPhoton[fNPhoton] = fTAPSParticles+i;   // Add to Photon Array
	        
			fPhotonEnergy[fNPhoton]	= fPhoton[fNPhoton]->GetT();
			fPhotonTheta[fNPhoton]	= fPhoton[fNPhoton]->GetThetaDg();
			fPhotonPhi[fNPhoton]	= fPhoton[fNPhoton]->GetPhiDg();
			fPhotonTime[fNPhoton]	= fPhoton[fNPhoton]->GetTime();	     
			fPhotonApp[fNPhoton]	= 2;   

	        fNPhoton++;								// Increase Photon counter
	        break;

		case kProton:                               	// Identified as a Proton
	        fProton[fNProton] 	= fTAPSParticles+i;     // Add to Proton Array
			fProtonEnergy[fNPhoton]	= fProton[fNProton]->GetT();
			fProtonTheta[fNPhoton]	= fProton[fNProton]->GetThetaDg();
			fProtonPhi[fNPhoton]	= fProton[fNProton]->GetPhiDg();
			fProtonTime[fNPhoton]	= fProton[fNProton]->GetTime();	
			fProtonApp[fNPhoton]	= 2;   
			   	        
	        fNProton++;					// Increase Proton counter
	        break;

		default:
	        fUnknown[fNUnknown] 	= fTAPSParticles+i;    	// include in "Unknown" list
	        fNUnknown++; 									// Increase "Unknown" counter

	        }						// End switch/case function
	}								// End # of Particle Loop
	}
	else fTAPSNParticle 	= 0;	// End If(fTAPS)

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

void TA2OnlinePhys::FillEndBuffers()
{
	// Apply BufferEnd to the end of all arrays
	fPhotonEnergy[fNPhoton]		= EBufferEnd;
	fPhotonTheta[fNPhoton]		= EBufferEnd;
	fPhotonPhi[fNPhoton]		= EBufferEnd;
	fPhotonTime[fNPhoton]		= EBufferEnd;

	fProtonEnergy[fNProton]		= EBufferEnd;
	fProtonTheta[fNProton]		= EBufferEnd;
	fProtonPhi[fNProton]		= EBufferEnd;
	fProtonTime[fNProton]		= EBufferEnd;

	fTaggerTime[fNTagg]			= EBufferEnd;
	fTaggerChannel[fNTagg]		= EBufferEnd;
	
}

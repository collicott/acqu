//--Author	C Collicott   14th Nov 2004   Most Basic Form
//
// PDG codes of particles generlly observed MAMI-C
// kElectron 11,     kPositron -11
// kMuonMinus 13     kMuonPlus -13      kGamma 22
// kPi0 111          kPiPlus 211        kPiMinus -211       kEta 221
// kProton 2212      kNeutron 2112
// 

#ifndef __TA2OnlinePhys_h__
#define __TA2OnlinePhys_h__

#include "TAcquRoot.h"
#include "TAcquFile.h"
#include "TA2Physics.h"
#include "TA2Analysis.h"
#include "TA2Tagger.h"
#include "TA2CrystalBall.h"
#include "TA2CentralApparatus.h"
#include "TA2Taps.h"
#include "TA2Ladder.h"
//#include "TA2PhotoPhysics.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <iostream>

class TA2Tagger;

class TA2OnlinePhys : public TA2Physics {

	protected:

// Begin by initialising Detectors
	TA2Tagger*		fTAGG;	// Glasgow photon tagger
	TA2Ladder*		fLadder;	// Ladder
	TA2CentralApparatus*	fCB;	// CentApp (PID, MWPCs, NaI)
	TA2Taps*		fTAPS;  // TAPS

// Tree Files
	TFile* 		fFile;
	TTree* 		fTree;
	Int_t		fProduceTreeFile;
	Char_t   	fTreeFileName[256];

// Pi0Compton Class Variables
	UInt_t 		i,j;

	// Particle Information

	UInt_t 		fCBMaxNParticle; 	// Max# Particle from CentApp
	UInt_t		fTAPSMaxNParticle; 	// Max# Particle from TAPS 
	UInt_t 		fMaxNParticle;		// Max # Particle (CB + TAPS)

	TA2Particle* 	fCBParticles;	   	// Particles from CB system
	TA2Particle* 	fTAPSParticles;	   	// Particles from TAPS
	TA2Particle* 	fTAGGParticles;	   	// Particles from Tagger

	UInt_t 		fCBNParticle;    	// # Particle from CentApp
	UInt_t		fTAPSNParticle;    	// # Particle from TAPS
	UInt_t		fTAGGNParticle; 	// # Particle from Tagger
	UInt_t 		fNParticle;		// # of Particles (CB + TAPS)

	// Particle Arrays

	UInt_t 		fNPhotTemp;		// # of Photons before pi0 ident
	UInt_t 		fNPhoton;		// # of Photons
	UInt_t 		fNProton;		// # of Protons
	UInt_t 		fNPi0;			// # of Pi0
	UInt_t 		fNUnknown;		// # of Unknowns

	TA2Particle**   fTaggedPhoton;		// Array to hold Tagger photons
	TA2Particle**	fPhotTemp;		// Array to hold photons before pi0 ident
	TA2Particle** 	fPhoton; 		// Array to hold photons	
	TA2Particle** 	fProton;		// Array to hold protons
	TA2Particle** 	fPi0;			// Array to hold Pi0s
	TA2Particle** 	fUnknown;		// Array to hold unknowns

	Double_t*	fPhotonEnergy;		// Energy of all photons
	Double_t*	fPhotonTheta;		// Theta of all photons
	Double_t*	fPhotonPhi;		// Phi of all photons
	Double_t*	fPhotonTime;		// Time of all photons
	Int_t*		fPhotonApp;

	Double_t*	fProtonEnergy;		// Energy of all protons
	Double_t*	fProtonTheta;		// Theta of all protons
	Double_t*	fProtonPhi;		// Phi of all protons
	Double_t*	fProtonTime;		// Time of all protons
	Int_t*		fProtonApp;

/*	Double_t*	fPi0Energy;		// Energy of all Pi0
	Double_t*	fPi0Theta;		// Theta of all Pi0
	Double_t*	fPi0Phi;		// Phi of all Pi0
	Double_t*	fPi0Time;		// Time of all Pi0

	Int_t		fNPionPhoton;		// # of photons which reconstruted into pion
	Bool_t*		fIsPionPhoton;		// Array to mark pion photons
	Int_t 		fPi0InvMassCut1;	//
	Int_t		fPi0InvMassCut2;*/

	UInt_t		fNTagg;
	Int_t*		fTaggerChannel;
	Double_t*	fTaggerTime;
/*	Double_t*	fTaggerPhotonTime;
	Double_t*	fTaggerPi0Time;

	Int_t	 	fN2PhotonInvariantMass;
	Double_t* 	f2PhotonInvariantMass;

	Int_t		fPhotTimePL;
	Int_t		fPhotTimePR;
	Int_t		fPhotTimeRL1;
	Int_t		fPhotTimeRR1;
	Int_t		fPhotTimeRL2;
	Int_t		fPhotTimeRR2;

	Int_t		fNTaggNPhot;
	Int_t		fNPrompt;
	Int_t		fNRandom;

	Double_t	fPromptRandomRatio;
	Int_t*		fTaggerChannelPrompt;
	Int_t*		fTaggerChannelRandom;
	Double_t*	fMissingMassPrompt;
	Double_t*	fMissingMassRandom;
	Double_t*	fPhotonThetaPrompt;
	Double_t*	fPhotonThetaRandom;
	Double_t*	fPhotonPhiPrompt;
	Double_t*	fPhotonPhiRandom;

	Int_t		fPi0TimePL;
	Int_t		fPi0TimePR;
	Int_t		fPi0TimeRL1;
	Int_t		fPi0TimeRR1;
	Int_t		fPi0TimeRL2;
	Int_t		fPi0TimeRR2;

	Int_t		fNTaggNPi0;
	Int_t		fNPromptPi0;
	Int_t		fNRandomPi0;

	Double_t	fPromptRandomRatioPi0;
	Int_t*		fTaggerChannelPromptPi0;
	Int_t*		fTaggerChannelRandomPi0;
	Double_t*	fMissingMassPromptPi0;
	Double_t*	fMissingMassRandomPi0;
	Double_t*	fPi0ThetaPrompt;
	Double_t*	fPi0ThetaRandom;
	Double_t*	fPi0PhiPrompt;
	Double_t*	fPi0PhiRandom;

	Double_t	fCBESum;
	Int_t		fNaINCluster;
	Int_t		fBaF2NCluster;*/
	
// Histograms 
	TH1* 	IM_2g;
	TH1*	IM_2g_CB;
	TH1*	IM_2g_TAPS;	
	TH1*	IM_2g_mix;

	TH1*	IM_3g;
	TH1*	IM_3g_CB;
	TH1*	IM_3g_E300;
	
	TH1*	IM_6g;
	TH1*	IM_6g_CB;
	
	
	protected:
		void ZeroCounters();
		void GetCBParticles();
		void GetTAPSParticles();
		void GetTagger();
		void FillEndBuffers();
		void DefineHistograms();
		void BasicPhysCheck();

	public:

	TA2OnlinePhys( const char*, TA2Analysis* );
		virtual ~TA2OnlinePhys();
		virtual void LoadVariable();	// variables for display/cuts
		virtual void PostInit( );	// init after parameter input
		virtual void SetConfig(Char_t*, Int_t);
		virtual void Reconstruct();	// reconstruct detector info
		virtual TA2DataManager* CreateChild( const char*, Int_t ){ return NULL;}		

	ClassDef(TA2OnlinePhys,1)
};

#endif

// Physics Module for TsEmDNAPhysicsActivator

#include "TsEmDNAPhysicsActivator.hh"

#include "TsParameterManager.hh"

#include "G4EmConfigurator.hh"
#include "G4EmParameters.hh"
#include "G4LossTableManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4VEnergyLossProcess.hh"

#include "G4Alpha.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "G4DummyModel.hh"
#include "G4EmProcessSubType.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProcessManager.hh"

#include "G4BetheBlochModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BraggModel.hh"
#include "G4ICRU49NuclearStoppingModel.hh"
#include "G4IonCoulombScatteringModel.hh"
#include "G4IonFluctuations.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4LowECapture.hh"
#include "G4LowEWentzelVIModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UrbanMscModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4eeToTwoGammaModel.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"

// Processes and models for Geant4-DNA
#include "G4DNAGenericIonsManager.hh"

#include "G4DNAChampionElasticModel.hh"
#include "G4DNAELSEPAElasticModel.hh"
#include "G4DNAElastic.hh"
#include "G4DNAVacuumModel.hh"
//#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAIonElasticModel.hh"
#include "G4DNAModelInterface.hh"
#include "G4DNAPTBElasticModel.hh"
#include "G4DNAUeharaScreenedRutherfordElasticModel.hh"
#include "G4eMultipleScattering.hh"

#include "G4DNAAttachment.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"
#include "G4DNAExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAPTBExcitationModel.hh"
#include "G4DNAPTBIonisationModel.hh"
#include "G4DNAVibExcitation.hh"

#include "TsDNARuddIonisationExtendedModel.hh"
#include "TsDNARuddIonisationExtendedRITRACKSModel.hh"

#include "G4SystemOfUnits.hh"
#include <vector>

#include "G4DNAChemistryManager.hh"
#include "G4DNAElectronSolvation.hh"
#include "G4DNAEmfietzoglouExcitationModel.hh"
#include "G4DNAEmfietzoglouIonisationModel.hh"

#include "G4DNACPA100ExcitationModel.hh"
#include "G4DNACPA100IonisationModel.hh"
#include <G4DNACPA100ElasticModel.hh>

#include "G4Threading.hh"

#include "G4UAtomicDeexcitation.hh" // added by MJPietrzak

TsEmDNAPhysicsActivator::TsEmDNAPhysicsActivator(G4int verbose)
	: G4VPhysicsConstructor("TsEmDNAPhysicsActivator")
{
	SetVerboseLevel(verbose);

	fEmParams = G4EmParameters::Instance();
	fEmParams->ActivateDNA();
}

TsEmDNAPhysicsActivator::TsEmDNAPhysicsActivator(TsParameterManager* pM)
	: G4VPhysicsConstructor("TsEmDNAPhysicsActivator"), fPm(pM)
{
	SetVerboseLevel(2);

	fEmParams = G4EmParameters::Instance();
	fEmParams->ActivateDNA();
}

void TsEmDNAPhysicsActivator::ConstructParticle()
{
	// bosons
	G4Gamma::Gamma();

	// leptons
	G4Electron::Electron();
	G4Positron::Positron();

	// baryons
	G4Proton::Proton();

	G4GenericIon::GenericIonDefinition();
	G4Alpha::Alpha();

	G4DNAGenericIonsManager* genericIonsManager;
	genericIonsManager = G4DNAGenericIonsManager::Instance();
	genericIonsManager->GetIon("alpha+");
	genericIonsManager->GetIon("helium");
	genericIonsManager->GetIon("hydrogen");
	// genericIonsManager->GetIon("carbon");
	// genericIonsManager->GetIon("nitrogen");
	// genericIonsManager->GetIon("oxygen");
	// genericIonsManager->GetIon("iron");
}

void TsEmDNAPhysicsActivator::ConstructProcess()
{
	const std::vector<G4String>& regnamesDNA = fEmParams->RegionsDNA();
	G4int nreg = regnamesDNA.size();
	if (0 == nreg)
	{
		if (GetVerboseLevel() > 0)
		{
			G4cout << "TsEmDNAPhysicsActivator::ConstructProcess: No regions for " << GetPhysicsName()
				   << ", thus not constructing any processes" << G4endl;
		}
		return;
	}
	const std::vector<G4String>& typesDNA = fEmParams->TypesDNA();

	if (GetVerboseLevel() > 0)
	{
		G4cout << "### TsEmDNAPhysicsActivator::ConstructProcess for " << nreg
			   << " regions; DNA physics type " << typesDNA[0] << G4endl;
	}

	// processes are defined with dummy models for the world
	AddDummyModels();

	// loop over regions
	for (G4int i = 0; i < nreg; ++i)
	{
		G4String reg = regnamesDNA[i];
		if (GetVerboseLevel() > 0)
		{
			G4cout << "### DNA models type " << typesDNA[i]
				   << " are being activated for G4Region " << reg << G4endl;
		}

		if (typesDNA[i].find("G4EmDNAPhysics") != std::string::npos)
		{
			if (GetVerboseLevel() > 0)
			{
				G4cout << "Adding models for a 'standard' em-dna: " << typesDNA[i] << G4endl;
				G4cout << "The consistency of those models with the physics lists has not been verified." << G4endl;
			}
			AddModels_G4EmDNA(reg, typesDNA[i]);
		}
		else if (typesDNA[i].find("TsEmDNAPhysics") != std::string::npos)
		{
			if (GetVerboseLevel() > 0)
			{
				G4cout << "Adding models from TOPAS-nBio dna physics: " << typesDNA[i] << G4endl;
			}

			fEmParams->SetDefaults();
			fEmParams->SetFluo(true);
			fEmParams->SetAuger(true);
			fEmParams->SetAugerCascade(true);
			fEmParams->SetDeexcitationIgnoreCut(true);

			if (typesDNA[i].compare("TsEmDNAPhysics") == 0)
			{
				AddModels_TsEmDNA(reg, typesDNA[i]);
			}
			// else if (typesDNA[i].compare("TsEmDNAPhysics_opt1") == 0)
			// {}
			else
			{
				G4cerr << "Topas is exiting due to a serious error in modular physics setup." << G4endl;
				G4cerr << "It was requested to activate EM modules from " << typesDNA[i] << "," << G4endl;
				G4cerr << "Which were not found by the TsEmDNAPhysicsActivator." << G4endl;
				fPm->AbortSession(1);
			}
		}
	}

	auto ltman = G4LossTableManager::Instance();
	ltman->EmConfigurator()->AddModels();
}

void TsEmDNAPhysicsActivator::AddModels_TsEmDNA(const G4String& reg, const G4String& name)
{
	AddElectronModels_TsEmDNA(reg, name);
	AddProtonModels_TsEmDNA(reg, name);
	AddIonModels_TsEmDNA(reg, name);
	AddPositronModels_TsEmDNA(reg);
	AddPhotonModels_Livermore(reg);

	G4LossTableManager::Instance()->SetAtomDeexcitation(new G4UAtomicDeexcitation);
}

void TsEmDNAPhysicsActivator::AddModels_G4EmDNA(const G4String& reg, const G4String& type)
{
	// list of particles
	const G4ParticleDefinition* elec = G4Electron::Electron();
	const G4ParticleDefinition* prot = G4Proton::Proton();
	//	const G4ParticleDefinition* gion = G4GenericIon::GenericIon();

	G4DNAGenericIonsManager* genericIonsManager = G4DNAGenericIonsManager::Instance();
	const G4ParticleDefinition* alpha2 = G4Alpha::Alpha();
	const G4ParticleDefinition* alpha1 = genericIonsManager->GetIon("alpha+");
	//	const G4ParticleDefinition* alpha0 = genericIonsManager->GetIon("helium");
	//	const G4ParticleDefinition* h0 = genericIonsManager->GetIon("hydrogen");

	G4ProcessManager* eman = elec->GetProcessManager();
	G4ProcessManager* pman = prot->GetProcessManager();
	//	G4ProcessManager* iman = gion->GetProcessManager();
	G4ProcessManager* a2man = alpha2->GetProcessManager();
	G4ProcessManager* a1man = alpha1->GetProcessManager();
	//	G4ProcessManager* a0man = alpha0->GetProcessManager();
	//	G4ProcessManager* h0man = h0->GetProcessManager();

	// alpha1 standard processes
	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
	G4ParticleDefinition* alpha11 = const_cast<G4ParticleDefinition*>(alpha1);
	ph->RegisterProcess(new G4hMultipleScattering(), alpha11);
	ph->RegisterProcess(new G4hIonisation(), alpha11);

	// limits for DNA model applicability
	//  static const G4double elowest= 7.4 * eV;  // seems to be option dependent - MJPietrzak
	const G4double elimel = 1 * MeV;
	const G4double pminbb = 2 * MeV;
	//	const G4double pmin = 0.1 * keV;
	const G4double pmax = 100 * MeV;
	//	const G4double hemin = 1 * keV;
	//	const G4double ionmin = 0.5 * MeV;

	if (type.find("option0") != std::string::npos)
	{
		AddElectronModels_G4EmDNA_opt0(reg, HasMsc(eman), elimel); // should this be option1?
	}
	else if (type.find("option2") != std::string::npos)
	{
		AddElectronModels_G4EmDNA_opt2(reg, HasMsc(eman), elimel);
	}
	else if (type.find("option4") != std::string::npos)
	{
		AddElectronModels_G4EmDNA_opt4(reg, HasMsc(eman), elimel);
	}
	else if (type.find("option4a") != std::string::npos)
	{
		AddElectronModels_G4EmDNA_opt4a(reg, HasMsc(eman), elimel);
	}
	else if (type.find("option6") != std::string::npos)
	{
		AddElectronModels_G4EmDNA_opt6(reg, HasMsc(eman), elimel);
	}
	else if (type.find("option6a") != std::string::npos)
	{
		AddElectronModels_G4EmDNA_opt6a(reg, HasMsc(eman), elimel);
	}
	else if (type.find("option7") != std::string::npos)
	{
		AddElectronModels_G4EmDNA_opt7(reg, HasMsc(eman), elimel);
	}
	else
	{
		G4cerr << "Topas is exiting due to a serious error in modular physics setup." << G4endl;
		G4cerr << "It was requested to activate EM modules from " << type << "," << G4endl;
		G4cerr << "Which were not found by the TsEmDNAPhysicsActivator." << G4endl;
		fPm->AbortSession(1);
	}

	// models for for other particles seems to be  independent so I moved them here - MJPietrzak
	AddProtonModels_G4EmDNA(reg, HasMsc(pman), elimel, pminbb, pmax);
	AddHeliumModels_G4EmDNA(reg, HasMsc(a1man), HasMsc(a2man), elimel, pminbb, pmax);
	AddGenericIonModels_G4EmDNA(reg, pminbb);
	DeactivateNuclearStopping_G4EmDNA(pman, elimel);
	DeactivateNuclearStopping_G4EmDNA(a1man, elimel);
	DeactivateNuclearStopping_G4EmDNA(a2man, elimel);
}

void TsEmDNAPhysicsActivator::AddElectronModels_TsEmDNA(const G4String& reg, const G4String& name)
{
	/*
	 em_config->SetExtraEmModel(particle, process, model, region, model_lowlimit, model_highlimit);
	 At run time the inner min/max limits are applied, so for consistency and correct
	 display at initialization, set the limits to the model limits.
	 Activation limits must be >= / <= model limits.
	 */

	G4EmConfigurator* em_config = G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* model;

	// This will be the same as the scattering model's low energy limit
	G4double solvationHighEnergyLimit = 0;

	// Elastic Scattering
	G4String eScatteringModel = "champion";
	if (fPm->ParameterExists("Ph/" + name + "/Electron/SetElasticScatteringModel"))
		eScatteringModel = fPm->GetStringParameter("Ph/" + name + "/Electron/SetElasticScatteringModel");
	eScatteringModel.toLower();

	if (eScatteringModel == "champion")
	{
		model = new G4DNAChampionElasticModel();
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("e-", "e-_G4DNAElastic", model, reg,
								   model->LowEnergyActivationLimit(), model->HighEnergyActivationLimit());
		solvationHighEnergyLimit = model->LowEnergyLimit(); // 7.4 eV
	}
	else if (eScatteringModel == "elsepa")
	{
		model = new G4DNAELSEPAElasticModel();
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("e-", "e-_G4DNAElastic", model, reg,
								   model->LowEnergyActivationLimit(), model->HighEnergyActivationLimit());
		solvationHighEnergyLimit = model->LowEnergyLimit(); // 10 eV
	}
	else if (eScatteringModel == "screenedrutherford")
	{
		model = new G4DNAScreenedRutherfordElasticModel();
		model->SetLowEnergyLimit(9 * eV);
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("e-", "e-_G4DNAElastic", model, reg,
								   model->LowEnergyActivationLimit(), model->HighEnergyActivationLimit());
		solvationHighEnergyLimit = model->LowEnergyLimit(); // 9 eV
	}
	else if (eScatteringModel == "ueharascreenedrutherford")
	{
		model = new G4DNAUeharaScreenedRutherfordElasticModel();
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("e-", "e-_G4DNAElastic", model, reg,
								   model->LowEnergyActivationLimit(), model->HighEnergyActivationLimit());
		solvationHighEnergyLimit = model->LowEnergyLimit(); // 9 eV
	}
	else if (eScatteringModel == "cpa100")
	{
		model = new G4DNACPA100ElasticModel();
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(255.955 * keV);
		em_config->SetExtraEmModel("e-", "e-_G4DNAElastic", model, reg,
								   model->LowEnergyActivationLimit(), model->HighEnergyActivationLimit());
		solvationHighEnergyLimit = model->LowEnergyLimit(); // 11 eV

		model = new G4DNAChampionElasticModel();
		model->SetActivationLowEnergyLimit(255.955 * keV);
		model->SetActivationHighEnergyLimit(fEmMaxElectron);
		em_config->SetExtraEmModel("e-", "e-_G4DNAElastic", model, reg);
	}
	else if (eScatteringModel.find("dna") == std::string::npos)
	{
		if (eScatteringModel == "dnascreenedrutherford")
		{
			model = new G4DNAModelInterface("e-_elastic_interaction");
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAScreenedRutherfordElasticModel(), G4Electron::Electron());
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAPTBElasticModel("THF/TMP/PY", G4Electron::Electron()));
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAVacuumModel());
			em_config->SetExtraEmModel("e-", "e-_G4DNAElastic", model, reg, 0, fEmMaxElectron);
		}
		else if (eScatteringModel == "dnacpa100")
		{
			model = new G4DNAModelInterface("e-_elastic_interaction");
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNACPA100ElasticModel(), G4Electron::Electron());
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAPTBElasticModel("THF/TMP/PY", G4Electron::Electron()));
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAVacuumModel());
			em_config->SetExtraEmModel("e-", "e-_G4DNAElastic", model, reg, 0, fEmMaxElectron);
		}
		else if (eScatteringModel == "dnaueharascreenedrutherford")
		{
			model = new G4DNAModelInterface("e-_elastic_interaction");
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAUeharaScreenedRutherfordElasticModel(), G4Electron::Electron());
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAPTBElasticModel("THF/TMP/PY", G4Electron::Electron()));
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAVacuumModel());
			em_config->SetExtraEmModel("e-", "e-_G4DNAElastic", model, reg, 0, fEmMaxElectron);
		}
		else
		{
			G4cerr << "Topas is exiting due to a serious error in modular physics setup." << G4endl;
			G4cerr << "Could not find the elastic scattering model:" << eScatteringModel << G4endl;
			fPm->AbortSession(1);
		}
	}
	else
	{ // WentzelVI
		model = new G4LowEWentzelVIModel();
		em_config->SetExtraEmModel("e-", "msc", model, reg, 0, fEmMaxElectron);
	}

	// Excitation
	G4String eExcitationModel = "born";
	if (fPm->ParameterExists("Ph/" + name + "/Electron/SetExcitationModel"))
		eExcitationModel = fPm->GetStringParameter("Ph/" + name + "/Electron/SetExcitationModel");
	eExcitationModel.toLower();

	G4String eIonisationModel = "born";
	if (fPm->ParameterExists("Ph/" + name + "/Electron/SetIonisationModel"))
		eIonisationModel = fPm->GetStringParameter("Ph/" + name + "/Electron/SetIonisationModel");
	eIonisationModel.toLower();

	if (eExcitationModel == "emfietzoglou")
	{
		model = new G4DNAEmfietzoglouExcitationModel();
		model->SetLowEnergyLimit(8 * eV);
		model->SetHighEnergyLimit(10 * keV);
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation", model, reg, 8 * eV, 10 * keV);

		model = new G4DNABornExcitationModel();
		model->SetLowEnergyLimit(10 * keV);
		model->SetHighEnergyLimit(1 * MeV);
		model->SetActivationLowEnergyLimit(10 * keV);
		model->SetActivationHighEnergyLimit(1 * MeV);
		em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation", model, reg, 10 * eV, 1 * MeV);
	}
	else if (eExcitationModel == "cpa100")
	{
		model = new G4DNACPA100ExcitationModel();
		model->SetLowEnergyLimit(11 * eV);
		model->SetHighEnergyLimit(255.955 * keV);
		model->SetActivationLowEnergyLimit(11 * eV);
		model->SetActivationHighEnergyLimit(255.955 * keV);
		em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation", model, reg, 11 * eV, 255.955 * keV);

		model = new G4DNABornExcitationModel();
		model->SetLowEnergyLimit(255.955 * keV);
		model->SetHighEnergyLimit(1 * MeV);
		model->SetActivationLowEnergyLimit(255.955 * keV);
		model->SetActivationHighEnergyLimit(1 * MeV);
		em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation", model, reg, 255.955 * eV, 1 * MeV);
	}
	else if (eExcitationModel.contains("dna"))
	{
		if (eExcitationModel == "dnaemfietzoglou")
		{
			model = new G4DNAModelInterface("e-_excitation_interaction");
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAEmfietzoglouExcitationModel(), G4Electron::Electron());
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAPTBExcitationModel("THF/TMP/PY", G4Electron::Electron()));
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAVacuumModel());
			em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation", model, reg, 0, fEmMaxElectron);
		}
		else if (eExcitationModel == "dnaborn")
		{
			model = new G4DNAModelInterface("e-_excitation_interaction");
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNABornExcitationModel(), G4Electron::Electron());
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAPTBExcitationModel("THF/TMP/PY", G4Electron::Electron()));
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAVacuumModel());
			em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation", model, reg, 0, fEmMaxElectron);
		}
		else
		{
			G4cerr << "Topas is exiting due to a serious error in modular physics setup." << G4endl;
			G4cerr << "Could not find the excitation model:" << eExcitationModel << G4endl;
			fPm->AbortSession(1);
		}
	}
	else
	{ // Born
		model = new G4DNABornExcitationModel();
		em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation", model, reg, 0, fEmMaxElectron);
	}

	// Ionisation
	G4bool useVarianceReduction = false;
	if (fPm->ParameterExists("Vr/UseG4DNAVarianceReduction"))
		useVarianceReduction = fPm->GetBooleanParameter("Vr/UseG4DNAVarianceReduction");

	if (eIonisationModel == "emfietzoglou")
	{
		// G4DNAIonisation*
		model = new G4DNAEmfietzoglouIonisationModel();
		model->SetLowEnergyLimit(10 * eV);
		model->SetHighEnergyLimit(10 * keV);
		model->SetActivationLowEnergyLimit(10 * eV);
		model->SetActivationHighEnergyLimit(10 * keV);
		em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation", model, reg, 10 * eV, 10 * keV);

		model = new G4DNABornIonisationModel();
		model->SetLowEnergyLimit(10 * keV);
		model->SetHighEnergyLimit(1 * MeV);
		model->SetActivationLowEnergyLimit(10 * keV);
		model->SetActivationHighEnergyLimit(1 * MeV);
		em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation", model, reg, 10 * keV, 1 * MeV);

		if (!useVarianceReduction)
		{
		}
		else
		{
			G4cerr << "Topas is exiting due to a serious error in modular physics setup." << G4endl;
			G4cerr << "Variance reduction is not yet implemented in the TsEmDNAPhysicsActivator:" << eExcitationModel << G4endl;
			fPm->AbortSession(1);
			/*
			G4int numberOfSplit = fPm->GetIntegerParameter("Vr/ParticleSplit/NumberOfSplit");
			G4String splitRegion = fPm->GetStringParameter("Vr/ParticleSplit/SplitElectronsInRegionNamed");

			G4cout << "-- Secondary split for electrons created in ionisation process actived "
			<< "with split number " << numberOfSplit << "-- " << G4endl;

			//G4ProcessManager* eman = particle->GetProcessManager();
			TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
			splitProcess->RegisterProcess(theDNAIonisationProcess);
			//eman->AddDiscreteProcess(splitProcess);
			ph->RegisterProcess(splitProcess, particle);
			*/
		}
	}
	else if (eIonisationModel == "cpa100")
	{
		model = new G4DNACPA100IonisationModel();
		model->SetLowEnergyLimit(11 * eV);
		model->SetHighEnergyLimit(255.955 * keV);
		model->SetActivationLowEnergyLimit(11 * eV);
		model->SetActivationHighEnergyLimit(255.955 * keV);
		em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation", model, reg, 11 * eV, 255.955 * keV);

		model = new G4DNABornIonisationModel();
		model->SetLowEnergyLimit(255.955 * keV);
		model->SetHighEnergyLimit(1 * MeV);
		model->SetActivationLowEnergyLimit(255.955 * keV);
		model->SetActivationHighEnergyLimit(1 * MeV);
		em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation", model, reg, 255.955 * keV, 1 * MeV);

		if (!useVarianceReduction)
		{
		}
		else
		{
			G4cerr << "Topas is exiting due to a serious error in modular physics setup." << G4endl;
			G4cerr << "Variance reduction is not yet implemented in the TsEmDNAPhysicsActivator:" << eExcitationModel << G4endl;
			/*
			G4int numberOfSplit = fPm->GetIntegerParameter("Vr/ParticleSplit/NumberOfSplit");
			G4String splitRegion = fPm->GetStringParameter("Vr/ParticleSplit/SplitElectronsInRegionNamed");

			G4cout << "-- Secondary split for electrons created in ionisation process actived "
			<< "with split number " << numberOfSplit << "-- " << G4endl;

			G4ProcessManager* eman = particle->GetProcessManager();
			TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
			splitProcess->RegisterProcess(theDNAIonisationProcess);
			eman->AddDiscreteProcess(splitProcess);
			*/
		}
	}
	else if (eIonisationModel.contains("dna"))
	{
		if (eIonisationModel == "dnaemfietzoglou")
		{
			model = new G4DNAModelInterface("e-_ionisation_interaction");
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAEmfietzoglouIonisationModel(), G4Electron::Electron());
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAPTBIonisationModel("THF/TMP/PY", G4Electron::Electron()));
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAVacuumModel());
			em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation", model, reg, 0, fEmMaxElectron);

			if (!useVarianceReduction)
			{
			}
			else
			{
				G4cerr << "Topas is exiting due to a serious error in modular physics setup." << G4endl;
				G4cerr << "Variance reduction is not yet implemented in the TsEmDNAPhysicsActivator:" << eExcitationModel << G4endl;
				/*
				G4int numberOfSplit = fPm->GetIntegerParameter("Vr/ParticleSplit/NumberOfSplit");
				G4String splitRegion = fPm->GetStringParameter("Vr/ParticleSplit/SplitElectronsInRegionNamed");

				G4cout << "-- Secondary split for electrons created in ionisation process actived "
				<< "with split number " << numberOfSplit << "-- " << G4endl;

				G4ProcessManager* eman = particle->GetProcessManager();
				TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
				splitProcess->RegisterProcess(theDNAIonisationProcess);
				eman->AddDiscreteProcess(splitProcess);
				*/
			}
		}
		else if (eIonisationModel == "dnaborn")
		{
			model = new G4DNAModelInterface("e-_ionisation_interaction");
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNABornIonisationModel(), G4Electron::Electron());
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAPTBIonisationModel("THF/TMP/PY", G4Electron::Electron()));
			dynamic_cast<G4DNAModelInterface*>(model)->RegisterModel(new G4DNAVacuumModel());
			em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation", model, reg, 0, fEmMaxElectron);

			if (!useVarianceReduction)
			{
			}
			else
			{
				G4cerr << "Topas is exiting due to a serious error in modular physics setup." << G4endl;
				G4cerr << "Variance reduction is not yet implemented in the TsEmDNAPhysicsActivator:" << eExcitationModel << G4endl;
				/*
				G4int numberOfSplit = fPm->GetIntegerParameter("Vr/ParticleSplit/NumberOfSplit");
				G4String splitRegion = fPm->GetStringParameter("Vr/ParticleSplit/SplitElectronsInRegionNamed");

				G4cout << "-- Secondary split for electrons created in ionisation process actived "
				<< "with split number " << numberOfSplit << "-- " << G4endl;

				G4ProcessManager* eman = particle->GetProcessManager();
				TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
				splitProcess->RegisterProcess(theDNAIonisationProcess);
				eman->AddDiscreteProcess(splitProcess);
				*/
			}
		}
	}
	else
	{
		model = new G4DNABornIonisationModel();
		model->SetHighEnergyLimit(1 * MeV);
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		dynamic_cast<G4DNABornIonisationModel*>(model)->SelectFasterComputation(true);
		em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation", model, reg,
								   model->LowEnergyLimit(), model->HighEnergyLimit());

		if (!useVarianceReduction)
		{
		}
		else
		{
			G4cerr << "Topas is exiting due to a serious error in modular physics setup." << G4endl;
			G4cerr << "Variance reduction is not yet implemented in the TsEmDNAPhysicsActivator:" << eExcitationModel << G4endl;
			/*
			G4int numberOfSplit = fPm->GetIntegerParameter("Vr/ParticleSplit/NumberOfSplit");
			G4String splitRegion = fPm->GetStringParameter("Vr/ParticleSplit/SplitElectronsInRegionNamed");

			G4cout << "-- Secondary split for electrons created in ionisation process actived "
			<< "with split number " << numberOfSplit << "-- " << G4endl;

			//G4ProcessManager* eman = particle->GetProcessManager();
			TsSplitProcessG4DNA* splitProcess = new TsSplitProcessG4DNA(splitRegion, numberOfSplit);
			splitProcess->RegisterProcess(theDNAIonisationProcess);
			//eman->AddDiscreteProcess(splitProcess);
			ph->RegisterProcess(splitProcess, particle);
			*/
		}
	}

	// Solvation
	if (fPm->ParameterExists("Ph/" + name + "/Electron/SetHighEnergyLimitForSolvation"))
		solvationHighEnergyLimit = fPm->GetDoubleParameter("Ph/" + name + "/Electron/SetHighEnergyLimitForSolvation", "Energy");
	model = G4DNASolvationModelFactory::GetMacroDefinedModel();
	if (solvationHighEnergyLimit > 0)
		model->SetHighEnergyLimit(solvationHighEnergyLimit);
	model->SetActivationLowEnergyLimit(0);
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	// Vibration
	G4bool activeVibExcitation = true;
	if (fPm->ParameterExists("Ph/" + name + "/Electron/ActiveVibExcitation"))
		activeVibExcitation = fPm->GetBooleanParameter("Ph/" + name + "/Electron/ActiveVibExcitation");
	if (activeVibExcitation)
	{
		model = new G4DNASancheExcitationModel();
		em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation", model, reg, 0, fEmMaxElectron);
	}

	// Attachment
	G4bool activeAttachment = true;
	if (fPm->ParameterExists("Ph/" + name + "/Electron/ActiveAttachment"))
		activeAttachment = fPm->GetBooleanParameter("Ph/" + name + "/Electron/ActiveAttachment");
	if (activeAttachment)
	{
		model = new G4DNAMeltonAttachmentModel();
		model->SetLowEnergyLimit(4. * eV);
		model->SetHighEnergyLimit(13. * eV);
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment", model, reg,
								   model->LowEnergyLimit(), model->HighEnergyLimit());
	}
}

void TsEmDNAPhysicsActivator::AddProtonModels_TsEmDNA(const G4String& reg, const G4String& name)
{
	G4EmConfigurator* em_config = G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* model;

	// Elastic Scattering
	G4String pScatteringModel = "default";
	if (fPm->ParameterExists("Ph/" + name + "/Proton/SetElasticScatteringModel"))
		pScatteringModel = fPm->GetStringParameter("Ph/" + name + "/Proton/SetElasticScatteringModel");

	pScatteringModel.toLower();
	if (pScatteringModel == "wentzelvi")
	{
		model = new G4LowEWentzelVIModel();
		em_config->SetExtraEmModel("proton", "msc", model, reg, 0, fEmMaxProton);
	}
	else
	{
		model = new G4DNAIonElasticModel();
		model->SetLowEnergyLimit(0 * eV);
		model->SetHighEnergyLimit(1 * MeV);
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("proton", "proton_G4DNAElastic", model, reg, 0 * eV, 1 * MeV);
	}

	// Excitation
	model = new G4DNAMillerGreenExcitationModel();
	model->SetLowEnergyLimit(10 * eV);
	model->SetHighEnergyLimit(500 * keV);
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("proton", "proton_G4DNAExcitation", model, reg, 10 * eV, 500 * keV);

	model = new G4DNABornExcitationModel();
	model->SetLowEnergyLimit(500 * keV);
	model->SetHighEnergyLimit(100 * MeV);
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("proton", "proton_G4DNAExcitation", model, reg, 500 * keV, 100 * MeV);

	// Ionisation
	G4String pIonisationModel = "default";
	if (fPm->ParameterExists("Ph/" + name + "/Proton/SetIonisationModel"))
		pIonisationModel = fPm->GetStringParameter("Ph/" + name + "/Proton/SetIonisationModel");
	pIonisationModel.toLower();

	if (pIonisationModel == "ritracks")
	{
		model = new G4DNARuddIonisationExtendedModel();
		model->SetLowEnergyLimit(0 * eV);
		model->SetHighEnergyLimit(500 * keV);
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation", model, reg, 0 * eV, 500 * keV);

		model = new TsDNARuddIonisationExtendedRITRACKSModel();
		model->SetLowEnergyLimit(500 * keV);
		model->SetHighEnergyLimit(500 * MeV);
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation", model, reg, 500 * keV, 500 * MeV);
	}
	else
	{
		model = new G4DNARuddIonisationExtendedModel();
		model->SetLowEnergyLimit(0 * eV);
		model->SetHighEnergyLimit(500 * keV);
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation", model, reg, 0 * eV, 500 * keV);

		model = new G4DNABornIonisationModel();
		model->SetLowEnergyLimit(500 * keV);
		model->SetHighEnergyLimit(100 * MeV);
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		dynamic_cast<G4DNABornIonisationModel*>(model)->SelectFasterComputation(true);
		em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation", model, reg, 500 * keV, 100 * MeV);

		model = new TsDNARuddIonisationExtendedModel();
		model->SetLowEnergyLimit(100 * MeV);
		model->SetHighEnergyLimit(500 * MeV);
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation", model, reg, 100 * MeV, 500 * MeV);

		model = new G4DNADingfelderChargeDecreaseModel();
		model->SetLowEnergyLimit(100 * eV);
		model->SetHighEnergyLimit(100 * MeV);
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("proton", "proton_G4DNAChargeDecrease", model, reg, 100 * eV, 100 * MeV);
	}
}

void TsEmDNAPhysicsActivator::AddIonModels_TsEmDNA(const G4String& reg, const G4String& name)
{
	G4EmConfigurator* em_config = G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* model;

	//
	// hydrogen
	model = new G4DNAIonElasticModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAElastic", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNAMillerGreenExcitationModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAExcitation", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNARuddIonisationModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAIonisation", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNADingfelderChargeIncreaseModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAChargeIncrease", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	//
	// alpha
	G4String alphaEModel = "default";
	if (fPm->ParameterExists("Ph/" + name + "/Alpha/SetElasticScatteringModel"))
		alphaEModel = fPm->ParameterExists("Ph/" + name + "/Alpha/SetElasticScatteringModel");
	alphaEModel.toLower();

	if (alphaEModel == "default")
	{
		model = new G4DNAIonElasticModel();
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("alpha", "alpha_G4DNAElastic", model, reg,
								   model->LowEnergyLimit(), model->HighEnergyLimit());
	}
	else
	{
		model = new G4LowEWentzelVIModel();
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("alpha", "msc", model, reg,
								   model->LowEnergyLimit(), model->HighEnergyLimit());
	}

	model = new G4DNAMillerGreenExcitationModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("alpha", "alpha_G4DNAExcitation", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNARuddIonisationExtendedModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("alpha", "alpha_G4DNAIonisation", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNADingfelderChargeDecreaseModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("alpha", "alpha_G4DNAChargeDecrease", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	//
	// alpha1
	alphaEModel = "default";
	if (fPm->ParameterExists("Ph/" + name + "/AlphaPlus/SetElasticScatteringModel"))
		alphaEModel = fPm->ParameterExists("Ph/" + name + "/AlphaPlus/SetElasticScatteringModel");
	alphaEModel.toLower();

	if (alphaEModel == "default")
	{
		model = new G4DNAIonElasticModel();
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAElastic", model, reg,
								   model->LowEnergyLimit(), model->HighEnergyLimit());
	}
	else
	{
		model = new G4LowEWentzelVIModel();
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("alpha+", "msc", model, reg,
								   model->LowEnergyLimit(), model->HighEnergyLimit());
	}

	model = new G4DNAMillerGreenExcitationModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAExcitation", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNARuddIonisationExtendedModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAIonisation", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNADingfelderChargeDecreaseModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAChargeDecrease", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNADingfelderChargeIncreaseModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAChargeIncrease", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	//
	// alpha2 (helium)
	model = new G4DNAIonElasticModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("helium", "helium_G4DNAElastic", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNAMillerGreenExcitationModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("helium", "helium_G4DNAExcitation", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNARuddIonisationExtendedModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("helium", "helium_G4DNAIonisation", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4DNADingfelderChargeIncreaseModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("helium", "helium_G4DNAChargeIncrease", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	//
	// generic ion
	G4String ionEModel = "default";
	if (fPm->ParameterExists("Ph/" + name + "/GenericIon/SetElasticScatteringModel"))
		ionEModel = fPm->GetStringParameter("Ph/" + name + "/GenericIon/SetElasticScatteringModel");
	ionEModel.toLower();

	if (ionEModel == "wentzelvi")
	{
		model = new G4LowEWentzelVIModel();
		model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
		model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
		em_config->SetExtraEmModel("GenericIon", "msc", model, reg,
								   model->LowEnergyLimit(), model->HighEnergyLimit());
	}

	if (fPm->ParameterExists("Ph/" + name + "/GenericIon/SetInelasticScatteringModel"))
		ionEModel = fPm->GetStringParameter("Ph/" + name + "/GenericIon/SetInelasticScatteringModel");
	ionEModel.toLower();

	if (ionEModel == "ruddextended")
		model = new TsDNARuddIonisationExtendedModel();
	else
		model = new G4DNARuddIonisationExtendedModel();
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("GenericIon", "GenericIon_G4DNAIonisation", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());
}

void TsEmDNAPhysicsActivator::AddPhotonModels_Livermore(const G4String& reg)
{
	G4EmConfigurator* em_config = G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* model;

	model = new G4LivermorePhotoElectricModel();
	// model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	// model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("gamma", "phot", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4LivermoreComptonModel();
	// model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	// model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("gamma", "compt", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4LivermoreGammaConversionModel();
	em_config->SetExtraEmModel("gamma", "conv", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4LivermoreRayleighModel();
	// model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	// model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("gamma", "Rayl", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());
}

#include "G4eIonisation.hh"
void TsEmDNAPhysicsActivator::AddPositronModels_TsEmDNA(const G4String& reg)
{
	G4EmConfigurator* em_config = G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* model;

	model = new G4UrbanMscModel;
	dynamic_cast<G4UrbanMscModel*>(model)->SetStepLimitType(fUseDistanceToBoundary);
	em_config->SetExtraEmModel("e+", "msc", model, reg);

	model = new G4MollerBhabhaModel;
	model->SetLowEnergyLimit(fEmParams->MinKinEnergy());
	model->SetHighEnergyLimit(fEmParams->MaxKinEnergy());
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	em_config->SetExtraEmModel("e+", "eIoni", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit(),
							   new G4UniversalFluctuation);

	model = new G4SeltzerBergerModel();
	model->SetLowEnergyLimit(fEmParams->MinKinEnergy());
	auto energyLimit = std::min(model->HighEnergyLimit(), 1 * GeV);
	model->SetHighEnergyLimit(energyLimit);
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	dynamic_cast<G4SeltzerBergerModel*>(model)->SetSecondaryThreshold(fEmParams->BremsstrahlungTh());
	dynamic_cast<G4SeltzerBergerModel*>(model)->SetLPMFlag(false);
	em_config->SetExtraEmModel("e+", "eBrem", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4eBremsstrahlungRelModel();
	model->SetLowEnergyLimit(energyLimit);
	model->SetHighEnergyLimit(fEmParams->MaxKinEnergy());
	model->SetActivationLowEnergyLimit(model->LowEnergyLimit());
	model->SetActivationHighEnergyLimit(model->HighEnergyLimit());
	dynamic_cast<G4eBremsstrahlungRelModel*>(model)->SetSecondaryThreshold(fEmParams->BremsstrahlungTh());
	dynamic_cast<G4eBremsstrahlungRelModel*>(model)->SetLPMFlag(fEmParams->LPM());
	em_config->SetExtraEmModel("e+", "eBrem", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());

	model = new G4eeToTwoGammaModel();
	model->SetLowEnergyLimit(fEmParams->MinKinEnergy());
	model->SetHighEnergyLimit(fEmParams->MaxKinEnergy());
	em_config->SetExtraEmModel("e+", "ePairProd", model, reg,
							   model->LowEnergyLimit(), model->HighEnergyLimit());
}

void TsEmDNAPhysicsActivator::AddDummyModels()
{
	const G4ParticleDefinition* elec = G4Electron::Electron();
	const G4ParticleDefinition* prot = G4Proton::Proton();
	const G4ParticleDefinition* gion = G4GenericIon::GenericIon();

	G4DNAGenericIonsManager* genericIonsManager = G4DNAGenericIonsManager::Instance();
	const G4ParticleDefinition* alpha2 = G4Alpha::Alpha();
	const G4ParticleDefinition* alpha1 = genericIonsManager->GetIon("alpha+");
	const G4ParticleDefinition* alpha0 = genericIonsManager->GetIon("helium");
	const G4ParticleDefinition* h0 = genericIonsManager->GetIon("hydrogen");

	G4ProcessManager* eman = elec->GetProcessManager();
	G4ProcessManager* pman = prot->GetProcessManager();
	G4ProcessManager* iman = gion->GetProcessManager();
	G4ProcessManager* a2man = alpha2->GetProcessManager();
	G4ProcessManager* a1man = alpha1->GetProcessManager();
	G4ProcessManager* a0man = alpha0->GetProcessManager();
	G4ProcessManager* h0man = h0->GetProcessManager();

	// electron solvation
	// e-
	G4DNAElectronSolvation* theDNAeSolvationProcess = new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
	theDNAeSolvationProcess->SetEmModel(new G4DummyModel());
	eman->AddDiscreteProcess(theDNAeSolvationProcess);

	// elastic scattering
	// e-
	G4DNAElastic* theDNAeElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
	theDNAeElasticProcess->SetEmModel(new G4DummyModel());
	eman->AddDiscreteProcess(theDNAeElasticProcess);

	// proton
	G4DNAElastic* theDNApElasticProcess = new G4DNAElastic("proton_G4DNAElastic");
	theDNApElasticProcess->SetEmModel(new G4DummyModel());
	pman->AddDiscreteProcess(theDNApElasticProcess);

	// alpha
	G4DNAElastic* theDNAa2ElasticProcess = new G4DNAElastic("alpha_G4DNAElastic");
	theDNAa2ElasticProcess->SetEmModel(new G4DummyModel());
	a2man->AddDiscreteProcess(theDNAa2ElasticProcess);

	// alpha1
	G4DNAElastic* theDNAa1ElasticProcess = new G4DNAElastic("alpha+_G4DNAElastic");
	theDNAa1ElasticProcess->SetEmModel(new G4DummyModel());
	a1man->AddDiscreteProcess(theDNAa1ElasticProcess);

	// alpha0 (helium)
	G4DNAElastic* theDNAa0ElasticProcess = new G4DNAElastic("helium_G4DNAElastic");
	theDNAa0ElasticProcess->SetEmModel(new G4DummyModel());
	a0man->AddDiscreteProcess(theDNAa0ElasticProcess);

	// hydrogen
	G4DNAElastic* theDNAh0ElasticProcess = new G4DNAElastic("hydrogen_G4DNAElastic");
	theDNAh0ElasticProcess->SetEmModel(new G4DummyModel());
	h0man->AddDiscreteProcess(theDNAh0ElasticProcess);

	// excitation
	// e-
	G4DNAExcitation* theDNAeExcProcess = new G4DNAExcitation("e-_G4DNAExcitation");
	theDNAeExcProcess->SetEmModel(new G4DummyModel());
	eman->AddDiscreteProcess(theDNAeExcProcess);

	// proton
	G4DNAExcitation* theDNApExcProcess = new G4DNAExcitation("proton_G4DNAExcitation");
	theDNApExcProcess->SetEmModel(new G4DummyModel());
	pman->AddDiscreteProcess(theDNApExcProcess);

	// alpha2
	G4DNAExcitation* theDNAa2ExcProcess = new G4DNAExcitation("alpha_G4DNAExcitation");
	theDNAa2ExcProcess->SetEmModel(new G4DummyModel());
	a2man->AddDiscreteProcess(theDNAa2ExcProcess);

	// alpha1
	G4DNAExcitation* theDNAa1ExcProcess = new G4DNAExcitation("alpha+_G4DNAExcitation");
	theDNAa1ExcProcess->SetEmModel(new G4DummyModel());
	a1man->AddDiscreteProcess(theDNAa1ExcProcess);

	// alpha0 (helium)
	G4DNAExcitation* theDNAa0ExcProcess = new G4DNAExcitation("helium_G4DNAExcitation");
	theDNAa0ExcProcess->SetEmModel(new G4DummyModel());
	a0man->AddDiscreteProcess(theDNAa0ExcProcess);

	// hydrogen
	G4DNAExcitation* theDNAh0ExcProcess = new G4DNAExcitation("hydrogen_G4DNAExcitation");
	theDNAh0ExcProcess->SetEmModel(new G4DummyModel());
	h0man->AddDiscreteProcess(theDNAh0ExcProcess);

	// vibration excitation
	// e-
	G4DNAVibExcitation* theDNAeVibExcProcess = new G4DNAVibExcitation("e-_G4DNAVibExcitation");
	theDNAeVibExcProcess->SetEmModel(new G4DummyModel());
	eman->AddDiscreteProcess(theDNAeVibExcProcess);

	// ionisation
	// e-
	G4DNAIonisation* theDNAeIoniProcess = new G4DNAIonisation("e-_G4DNAIonisation");
	theDNAeIoniProcess->SetEmModel(new G4DummyModel());
	eman->AddDiscreteProcess(theDNAeIoniProcess);

	// proton
	G4DNAIonisation* theDNApIoniProcess = new G4DNAIonisation("proton_G4DNAIonisation");
	theDNApIoniProcess->SetEmModel(new G4DummyModel());
	pman->AddDiscreteProcess(theDNApIoniProcess);

	// alpha2
	G4DNAIonisation* theDNAa2IoniProcess = new G4DNAIonisation("alpha_G4DNAIonisation");
	theDNAa2IoniProcess->SetEmModel(new G4DummyModel());
	a2man->AddDiscreteProcess(theDNAa2IoniProcess);

	// alpha1
	G4DNAIonisation* theDNAa1IoniProcess = new G4DNAIonisation("alpha+_G4DNAIonisation");
	theDNAa1IoniProcess->SetEmModel(new G4DummyModel());
	a1man->AddDiscreteProcess(theDNAa1IoniProcess);

	// alpha0 (helium)
	G4DNAIonisation* theDNAa0IoniProcess = new G4DNAIonisation("helium_G4DNAIonisation");
	theDNAa0IoniProcess->SetEmModel(new G4DummyModel());
	a0man->AddDiscreteProcess(theDNAa0IoniProcess);

	// hydrogen
	G4DNAIonisation* theDNAh0IoniProcess = new G4DNAIonisation("hydrogen_G4DNAIonisation");
	theDNAh0IoniProcess->SetEmModel(new G4DummyModel());
	h0man->AddDiscreteProcess(theDNAh0IoniProcess);

	// general ion
	G4DNAIonisation* theDNAiIoniProcess = new G4DNAIonisation("GenericIon_G4DNAIonisation");
	theDNAiIoniProcess->SetEmModel(new G4DummyModel());
	iman->AddDiscreteProcess(theDNAiIoniProcess);

	// attachment
	// e-
	G4DNAAttachment* theDNAAttachProcess = new G4DNAAttachment("e-_G4DNAAttachment");
	theDNAAttachProcess->SetEmModel(new G4DummyModel());
	eman->AddDiscreteProcess(theDNAAttachProcess);

	// charge exchange
	// proton
	G4DNAChargeDecrease* theDNApChargeDecreaseProcess = new G4DNAChargeDecrease("proton_G4DNAChargeDecrease");
	theDNApChargeDecreaseProcess->SetEmModel(new G4DummyModel());
	pman->AddDiscreteProcess(theDNApChargeDecreaseProcess);

	// alpha2
	G4DNAChargeDecrease* theDNAa2ChargeDecreaseProcess = new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease");
	theDNAa2ChargeDecreaseProcess->SetEmModel(new G4DummyModel());
	a2man->AddDiscreteProcess(theDNAa2ChargeDecreaseProcess);

	// alpha1
	G4DNAChargeDecrease* theDNAa1ChargeDecreaseProcess = new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease");
	theDNAa1ChargeDecreaseProcess->SetEmModel(new G4DummyModel());
	a1man->AddDiscreteProcess(theDNAa1ChargeDecreaseProcess);

	// alpha1
	G4DNAChargeIncrease* theDNAa1ChargeIncreaseProcess = new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease");
	theDNAa1ChargeIncreaseProcess->SetEmModel(new G4DummyModel());
	a1man->AddDiscreteProcess(theDNAa1ChargeIncreaseProcess);

	// alpha0 (helium)
	G4DNAChargeIncrease* theDNAa0ChargeIncreaseProcess = new G4DNAChargeIncrease("helium_G4DNAChargeIncrease");
	theDNAa0ChargeIncreaseProcess->SetEmModel(new G4DummyModel());
	a0man->AddDiscreteProcess(theDNAa0ChargeIncreaseProcess);

	// hydrogen
	G4DNAChargeIncrease* theDNAh0ChargeIncreaseProcess = new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease");
	theDNAh0ChargeIncreaseProcess->SetEmModel(new G4DummyModel());
	h0man->AddDiscreteProcess(theDNAh0ChargeIncreaseProcess);
}

void TsEmDNAPhysicsActivator::AddElectronModels_G4EmDNA_opt0(const G4String& reg,
															 G4bool emsc,
															 G4double elimel)
{
	G4EmConfigurator* em_config = G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* mod;

	static const G4double elowest = 7.4 * eV; // seems to be option dependent, so I moved it here - MJPietrzak
	static const G4double elimin = 1 * MeV;
	static const G4double elimvb = 100 * eV;
	static const G4double elimat = 13 * eV;
	static const G4double elim1 = 10 * keV;

	// for e- 100 MeV is a limit between different msc models
	G4double emax = fEmParams->MaxKinEnergy();

	if (emsc)
	{
		G4UrbanMscModel* msc = new G4UrbanMscModel();
		msc->SetActivationLowEnergyLimit(elimel);
		G4double emaxmsc = std::min(100 * MeV, emax);
		em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, emaxmsc);
	}
	else
	{
		mod = new G4eCoulombScatteringModel();
		mod->SetActivationLowEnergyLimit(elimel);
		em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, emax);
	}

	// cuts and solvation
	mod = new G4DNAOneStepThermalizationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation",
							   mod, reg, 0., elowest);

	// elastic
	mod = new G4DNAChampionElasticModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
							   mod, reg, 0.0, elimel);
	// ionisation
	mod = new G4MollerBhabhaModel();
	mod->SetActivationLowEnergyLimit(elimin);
	em_config->SetExtraEmModel("e-", "eIoni",
							   mod, reg, 0.0, emax,
							   new G4UniversalFluctuation());

	mod = new G4DNABornIonisationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
							   mod, reg, elim1, elimin);

	mod = new G4DNAEmfietzoglouIonisationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
							   mod, reg, 0.0, elim1);

	// exc
	mod = new G4DNAEmfietzoglouExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
							   mod, reg, 0.0, elim1);

	mod = new G4DNABornExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
							   mod, reg, elim1, elimin);

	mod = new G4DNASancheExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation",
							   mod, reg, 0.0, elimvb);

	mod = new G4DNAMeltonAttachmentModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment",
							   mod, reg, 0.0, elimat);

	mod = G4DNASolvationModelFactory::GetMacroDefinedModel();
	if (elowest > 0)
		mod->SetHighEnergyLimit(elowest);
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation", mod, reg);
}

void TsEmDNAPhysicsActivator::AddElectronModels_G4EmDNA_opt2(const G4String& reg,
															 G4bool emsc,
															 G4double elimel)
{
	G4EmConfigurator* em_config =
		G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* mod;

	static const G4double elowest = 7.4 * eV; // seems to be option dependent, so I moved it here - MJPietrzak
	static const G4double elimin = 1 * MeV;
	static const G4double elimvb = 100 * eV;
	static const G4double elimat = 13 * eV;

	// for e- 100 MeV is a limit between different msc models
	G4double emax = fEmParams->MaxKinEnergy();

	if (emsc)
	{
		G4UrbanMscModel* msc = new G4UrbanMscModel();
		msc->SetActivationLowEnergyLimit(elimel);
		G4double emaxmsc = std::min(100 * MeV, emax);
		em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, emaxmsc);
	}
	else
	{
		mod = new G4eCoulombScatteringModel();
		mod->SetActivationLowEnergyLimit(elimel);
		em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, emax);
	}

	// cuts and solvation
	mod = new G4DNAOneStepThermalizationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation",
							   mod, reg, 0., elowest);

	// elastic
	mod = new G4DNAChampionElasticModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
							   mod, reg, 0.0, elimel);

	// ionisation
	mod = new G4MollerBhabhaModel();
	mod->SetActivationLowEnergyLimit(elimin);
	em_config->SetExtraEmModel("e-", "eIoni",
							   mod, reg, 0.0, emax,
							   new G4UniversalFluctuation());

	mod = new G4DNABornIonisationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
							   mod, reg, elowest, elimin);

	// excitation
	mod = new G4DNABornExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
							   mod, reg, 0., elimin);

	// vib excitation
	mod = new G4DNASancheExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation",
							   mod, reg, 0.0, elimvb);

	// attachment
	mod = new G4DNAMeltonAttachmentModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment",
							   mod, reg, 0.0, elimat);

	mod = G4DNASolvationModelFactory::GetMacroDefinedModel();
	if (elowest > 0)
		mod->SetHighEnergyLimit(elowest);
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation", mod, reg);
}

void TsEmDNAPhysicsActivator::AddElectronModels_G4EmDNA_opt4(const G4String& reg, G4bool emsc, G4double elimel)
{
	G4EmConfigurator* em_config =
		G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* mod;

	static const G4double elowest = 10 * eV; // seems to be option dependent, so I moved it here - MJPietrzak
	static const G4double elimin = 1 * MeV;
	//  static const G4double elimvb = 100 * eV;
	//  static const G4double elimat = 13 * eV;

	// for e- 100 MeV is a limit between different msc models
	G4double emax = fEmParams->MaxKinEnergy();

	if (emsc)
	{
		G4UrbanMscModel* msc = new G4UrbanMscModel();
		msc->SetActivationLowEnergyLimit(elimel);
		G4double emaxmsc = std::min(100 * MeV, emax);
		em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, emaxmsc);
	}
	else
	{
		mod = new G4eCoulombScatteringModel();
		mod->SetActivationLowEnergyLimit(elimel);
		em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, emax);
	}

	// cuts and solvation
	mod = new G4DNAOneStepThermalizationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation",
							   mod, reg, 0., elowest);

	// elastic
	mod = new G4DNAUeharaScreenedRutherfordElasticModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
							   mod, reg, 0.0, elimel);

	// ionisation
	mod = new G4MollerBhabhaModel();
	mod->SetActivationLowEnergyLimit(elimin);
	em_config->SetExtraEmModel("e-", "eIoni",
							   mod, reg, 0.0, emax,
							   new G4UniversalFluctuation());

	mod = new G4DNAEmfietzoglouIonisationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
							   mod, reg, elowest, elimin);

	// excitation
	mod = new G4DNAEmfietzoglouExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
							   mod, reg, 0., elimin);

	// todo - MJPietrzak
	//   I don't understand why vib excit. and attachment models are turned off in option 4 (and 6)
	//   therefore I have created option 4a (and 6a), which has these models turned on
	//
	//  // vib excitation
	//  mod = new G4DNASancheExcitationModel();
	//  em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation",
	//			     mod, reg, 0.0, elimvb);
	//
	//  // attachment
	//  mod = new G4DNAMeltonAttachmentModel();
	//  em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment",
	//			     mod, reg, 0.0, elimat);
	mod = G4DNASolvationModelFactory::GetMacroDefinedModel();
	if (elowest > 0)
		mod->SetHighEnergyLimit(elowest);
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation", mod, reg);
}

void TsEmDNAPhysicsActivator::AddElectronModels_G4EmDNA_opt4a(const G4String& reg,
															  G4bool emsc,
															  G4double elimel)
{
	G4EmConfigurator* em_config =
		G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* mod;

	static const G4double elowest = 10 * eV; // seems to be option dependent, so I moved it here - MJPietrzak
	static const G4double elimin = 1 * MeV;
	static const G4double elimvb = 100 * eV;
	static const G4double elimat = 13 * eV;

	// for e- 100 MeV is a limit between different msc models
	G4double emax = fEmParams->MaxKinEnergy();

	if (emsc)
	{
		G4UrbanMscModel* msc = new G4UrbanMscModel();
		msc->SetActivationLowEnergyLimit(elimel);
		G4double emaxmsc = std::min(100 * MeV, emax);
		em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, emaxmsc);
	}
	else
	{
		mod = new G4eCoulombScatteringModel();
		mod->SetActivationLowEnergyLimit(elimel);
		em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, emax);
	}

	// cuts and solvation
	mod = new G4DNAOneStepThermalizationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation",
							   mod, reg, 0., elowest);

	// elastic
	mod = new G4DNAUeharaScreenedRutherfordElasticModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
							   mod, reg, 0.0, elimel);

	// ionisation
	mod = new G4MollerBhabhaModel();
	mod->SetActivationLowEnergyLimit(elimin);
	em_config->SetExtraEmModel("e-", "eIoni",
							   mod, reg, 0.0, emax,
							   new G4UniversalFluctuation());

	mod = new G4DNAEmfietzoglouIonisationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
							   mod, reg, elowest, elimin);

	// excitation
	mod = new G4DNAEmfietzoglouExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
							   mod, reg, 0., elimin);

	// I don't understand why vib excit. and attachment models are turned off in option 4
	// therefore I have created option 4a, which has these models turned on
	// and here it is - MJPietrzak

	// vib excitation
	mod = new G4DNASancheExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation",
							   mod, reg, 0.0, elimvb);

	// attachment
	mod = new G4DNAMeltonAttachmentModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment",
							   mod, reg, 0.0, elimat);

	mod = G4DNASolvationModelFactory::GetMacroDefinedModel();
	if (elowest > 0)
		mod->SetHighEnergyLimit(elowest);
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation", mod, reg);
}

void TsEmDNAPhysicsActivator::AddElectronModels_G4EmDNA_opt6(const G4String& reg,
															 G4bool emsc,
															 G4double elimel)
{
	G4EmConfigurator* em_config =
		G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* mod;

	static const G4double elowest = 11 * eV; // seems to be option dependent, so I moved it here - MJPietrzak
	static const G4double elimin = 1 * MeV;
	//  static const G4double elimvb = 100 * eV;
	//  static const G4double elimat = 13 * eV;

	// for e- 100 MeV is a limit between different msc models
	G4double emax = fEmParams->MaxKinEnergy();

	if (emsc)
	{
		G4UrbanMscModel* msc = new G4UrbanMscModel();
		msc->SetActivationLowEnergyLimit(elimel);
		G4double emaxmsc = std::min(100 * MeV, emax);
		em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, emaxmsc);
	}
	else
	{
		mod = new G4eCoulombScatteringModel();
		mod->SetActivationLowEnergyLimit(elimel);
		em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, emax);
	}

	// cuts and solvation
	mod = new G4DNAOneStepThermalizationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation",
							   mod, reg, 0., elowest);

	// elastic
	mod = new G4DNACPA100ElasticModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
							   mod, reg, 0.0, elimel);

	// ionisation
	mod = new G4MollerBhabhaModel();
	mod->SetActivationLowEnergyLimit(elimin);
	em_config->SetExtraEmModel("e-", "eIoni",
							   mod, reg, 0.0, emax,
							   new G4UniversalFluctuation());

	mod = new G4DNACPA100IonisationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
							   mod, reg, elowest, elimin);

	// excitation
	mod = new G4DNACPA100ExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
							   mod, reg, 0., elimin);

	// I don't understand why vib excit. and attachment models are turned off in option 6 (and 4)
	// therefore I have created option 6a (and 4a), which has these models turned on
	// and here it is - MJPietrzak

	//  // vib excitation
	//  mod = new G4DNASancheExcitationModel();
	//  em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation",
	//                             mod, reg, 0.0, elimvb);
	//
	//  // attachment
	//  mod = new G4DNAMeltonAttachmentModel();
	//  em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment",
	//                             mod, reg, 0.0, elimat);

	mod = G4DNASolvationModelFactory::GetMacroDefinedModel();
	if (elowest > 0)
		mod->SetHighEnergyLimit(elowest);
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation", mod, reg);
}

void TsEmDNAPhysicsActivator::AddElectronModels_G4EmDNA_opt6a(const G4String& reg,
															  G4bool emsc,
															  G4double elimel)
{
	G4EmConfigurator* em_config =
		G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* mod;

	static const G4double elowest = 11 * eV; // seems to be option dependent, so I moved it here - MJPietrzak
	static const G4double elimin = 1 * MeV;
	static const G4double elimvb = 100 * eV;
	static const G4double elimat = 13 * eV;

	// for e- 100 MeV is a limit between different msc models
	G4double emax = fEmParams->MaxKinEnergy();

	if (emsc)
	{
		G4UrbanMscModel* msc = new G4UrbanMscModel();
		msc->SetActivationLowEnergyLimit(elimel);
		G4double emaxmsc = std::min(100 * MeV, emax);
		em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, emaxmsc);
	}
	else
	{
		mod = new G4eCoulombScatteringModel();
		mod->SetActivationLowEnergyLimit(elimel);
		em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, emax);
	}

	// cuts and solvation
	mod = new G4DNAOneStepThermalizationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation",
							   mod, reg, 0., elowest);

	// elastic
	mod = new G4DNACPA100ElasticModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
							   mod, reg, 0.0, elimel);

	// ionisation
	mod = new G4MollerBhabhaModel();
	mod->SetActivationLowEnergyLimit(elimin);
	em_config->SetExtraEmModel("e-", "eIoni",
							   mod, reg, 0.0, emax,
							   new G4UniversalFluctuation());

	mod = new G4DNACPA100IonisationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
							   mod, reg, elowest, elimin);

	// excitation
	mod = new G4DNACPA100ExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
							   mod, reg, 0., elimin);

	// I don't understand why vib excit. and attachment models are turned off in option 6 (and 4)
	// therefore I have created option 6a (and 4a), which has these models turned on
	// and here it is - MJPietrzak

	// vib excitation
	mod = new G4DNASancheExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation",
							   mod, reg, 0.0, elimvb);

	// attachment
	mod = new G4DNAMeltonAttachmentModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment",
							   mod, reg, 0.0, elimat);

	mod = G4DNASolvationModelFactory::GetMacroDefinedModel();
	if (elowest > 0)
		mod->SetHighEnergyLimit(elowest);
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation", mod, reg);
}

void TsEmDNAPhysicsActivator::AddElectronModels_G4EmDNA_opt7(const G4String& reg,
															 G4bool emsc,
															 G4double elimel)
{
	G4EmConfigurator* em_config =
		G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* mod;

	static const G4double elowest = 10 * eV; // seems to be option dependent, so I moved it here - MJPietrzak
	static const G4double elimin = 1 * MeV;
	static const G4double elimvb = 100 * eV;
	static const G4double elimat = 13 * eV;
	static const G4double elim1 = 10 * keV;

	// for e- 100 MeV is a limit between different msc models
	G4double emax = fEmParams->MaxKinEnergy();

	if (emsc)
	{
		G4UrbanMscModel* msc = new G4UrbanMscModel();
		msc->SetActivationLowEnergyLimit(elimel);
		G4double emaxmsc = std::min(100 * MeV, emax);
		em_config->SetExtraEmModel("e-", "msc", msc, reg, 0.0, emaxmsc);
	}
	else
	{
		mod = new G4eCoulombScatteringModel();
		mod->SetActivationLowEnergyLimit(elimel);
		em_config->SetExtraEmModel("e-", "CoulombScat", mod, reg, 0.0, emax);
	}

	// cuts and solvation
	mod = new G4DNAOneStepThermalizationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation",
							   mod, reg, 0., elowest);

	// elastic
	mod = new G4DNAUeharaScreenedRutherfordElasticModel(); // G4DNAChampionElasticModel();  in Opt_0 - the main difference between Opt0 and Opt7
	em_config->SetExtraEmModel("e-", "e-_G4DNAElastic",
							   mod, reg, 0.0, elimel);

	// ionisation
	mod = new G4MollerBhabhaModel();
	mod->SetActivationLowEnergyLimit(elimin);
	em_config->SetExtraEmModel("e-", "eIoni",
							   mod, reg, 0.0, emax,
							   new G4UniversalFluctuation());
	// todo -  MJPietrzak
	//   I don't understand why the MollerBhabhaModel is here, as there is no sign of this model in regular DNA lists.
	//   However, it seems that it is necessary for all the options.

	mod = new G4DNABornIonisationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
							   mod, reg, elim1, elimin);

	mod = new G4DNAEmfietzoglouIonisationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation",
							   mod, reg, elowest, elim1); // mod, reg, 0.0, elim1);

	// excitation
	mod = new G4DNAEmfietzoglouExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
							   mod, reg, 8 * eV, elim1);

	mod = new G4DNABornExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation",
							   mod, reg, elim1, elimin);

	// vib excitation
	mod = new G4DNASancheExcitationModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation",
							   mod, reg, 0.0, elimvb);

	// attachment
	mod = new G4DNAMeltonAttachmentModel();
	em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment",
							   mod, reg, 0.0, elimat);

	mod = G4DNASolvationModelFactory::GetMacroDefinedModel();
	if (elowest > 0)
		mod->SetHighEnergyLimit(elowest);
	em_config->SetExtraEmModel("e-", "e-_G4DNAElectronSolvation", mod, reg);
}

void TsEmDNAPhysicsActivator::AddProtonModels_G4EmDNA(const G4String& reg,
													  G4bool pmsc, G4double elimel,
													  G4double pminbb, G4double pmax)
{
	G4EmConfigurator* em_config = G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* mod;

	static const G4double gmmax = 500 * keV;

	G4double emax = fEmParams->MaxKinEnergy();

	// proton

	// if SS physics list msc process does not exist
	if (pmsc)
	{
		G4WentzelVIModel* msc = new G4WentzelVIModel();
		msc->SetActivationLowEnergyLimit(elimel);
		em_config->SetExtraEmModel("proton", "msc", msc, reg, 0.0, emax);
	}
	// single scattering always applied
	mod = new G4eCoulombScatteringModel();
	mod->SetActivationLowEnergyLimit(elimel);
	em_config->SetExtraEmModel("proton", "CoulombScat", mod, reg, 0.0, emax);

	mod = new G4BraggModel();
	mod->SetActivationLowEnergyLimit(std::min(pminbb, pmax));
	em_config->SetExtraEmModel("proton", "hIoni",
							   mod, reg, 0.0, pminbb,
							   new G4UniversalFluctuation());

	mod = new G4BetheBlochModel();
	mod->SetActivationLowEnergyLimit(pmax);
	em_config->SetExtraEmModel("proton", "hIoni",
							   mod, reg, pminbb, emax,
							   new G4UniversalFluctuation());

	mod = new G4DNARuddIonisationModel();
	em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation",
							   mod, reg, 0.0, gmmax);

	mod = new G4DNABornIonisationModel();
	em_config->SetExtraEmModel("proton", "proton_G4DNAIonisation",
							   mod, reg, gmmax, pmax);

	mod = new G4DNAMillerGreenExcitationModel();
	em_config->SetExtraEmModel("proton", "proton_G4DNAExcitation",
							   mod, reg, 0.0, gmmax);

	mod = new G4DNABornExcitationModel();
	em_config->SetExtraEmModel("proton", "proton_G4DNAExcitation",
							   mod, reg, gmmax, pmax);

	mod = new G4DNADingfelderChargeDecreaseModel();
	em_config->SetExtraEmModel("proton", "proton_G4DNAChargeDecrease",
							   mod, reg, 0.0, pmax);

	mod = new G4DNAIonElasticModel();
	em_config->SetExtraEmModel("proton", "proton_G4DNAElastic",
							   mod, reg, 0.0, elimel);

	// hydrogen
	mod = new G4DNARuddIonisationModel();
	em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAIonisation",
							   mod, reg, 0.0, pmax);

	mod = new G4DNAMillerGreenExcitationModel();
	em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAExcitation",
							   mod, reg, 0.0, gmmax);

	mod = new G4DNADingfelderChargeIncreaseModel();
	em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAChargeIncrease",
							   mod, reg, 0.0, pmax);

	mod = new G4DNAIonElasticModel();
	em_config->SetExtraEmModel("hydrogen", "hydrogen_G4DNAElastic",
							   mod, reg, 0.0, elimel);
}

void TsEmDNAPhysicsActivator::AddGenericIonModels_G4EmDNA(const G4String& reg,
														  G4double pminbb)
{
	G4EmConfigurator* em_config =
		G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* mod;

	G4double emax = fEmParams->MaxKinEnergy();
	G4double iemax = std::min(10 * MeV, emax);
	// G4double iemin = 100*eV;

	mod = new G4BraggIonModel();
	mod->SetActivationLowEnergyLimit(iemax);
	em_config->SetExtraEmModel("GenericIon", "ionIoni",
							   mod, reg, 0.0, pminbb,
							   new G4IonFluctuations());

	mod = new G4BetheBlochModel();
	mod->SetActivationLowEnergyLimit(iemax);
	em_config->SetExtraEmModel("GenericIon", "ionIoni",
							   mod, reg, pminbb, emax,
							   new G4IonFluctuations());

	mod = new G4DNARuddIonisationExtendedModel();
	em_config->SetExtraEmModel("GenericIon", "GenericIon_G4DNAIonisation",
							   mod, reg, 0.0, iemax);
}

void TsEmDNAPhysicsActivator::AddHeliumModels_G4EmDNA(const G4String& reg,
													  G4bool a1msc,
													  G4bool a2msc, G4double elimel,
													  G4double pminbb, G4double)
{
	G4EmConfigurator* em_config =
		G4LossTableManager::Instance()->EmConfigurator();
	G4VEmModel* mod;

	static const G4double hemax = 400 * MeV;
	static const G4double massRatio = G4Alpha::Alpha()->GetPDGMass() / CLHEP::proton_mass_c2;

	G4double emax = fEmParams->MaxKinEnergy();
	G4double pminbba = massRatio * pminbb;
	if (GetVerboseLevel() > 0)
	{
		G4cout << "AddHeliumModels_G4EmDNA for <" << reg << "> a1msc: " << a1msc << " a2msc: " << a2msc
			   << " elimel= " << elimel << " pminbba= " << pminbba << G4endl;
	}
	// alpha2
	if (elimel < emax)
	{
		if (a2msc)
		{
			G4UrbanMscModel* msc = new G4UrbanMscModel();
			msc->SetActivationLowEnergyLimit(elimel);
			em_config->SetExtraEmModel("alpha", "msc", msc, reg, 0.0, emax);
		}
		else
		{
			mod = new G4IonCoulombScatteringModel();
			mod->SetActivationLowEnergyLimit(elimel);
			em_config->SetExtraEmModel("alpha", "CoulombScat", mod, reg, 0.0, emax);
		}
	}

	mod = new G4BraggIonModel();
	mod->SetActivationLowEnergyLimit(hemax / massRatio);
	em_config->SetExtraEmModel("alpha", "ionIoni",
							   mod, reg, 0.0, pminbba,
							   new G4IonFluctuations());

	mod = new G4BetheBlochModel();
	mod->SetActivationLowEnergyLimit(hemax / massRatio);
	em_config->SetExtraEmModel("alpha", "ionIoni",
							   mod, reg, pminbba, emax,
							   new G4IonFluctuations());

	mod = new G4DNARuddIonisationModel();
	em_config->SetExtraEmModel("alpha", "alpha_G4DNAIonisation",
							   mod, reg, 0.0, hemax);

	mod = new G4DNAMillerGreenExcitationModel();
	em_config->SetExtraEmModel("alpha", "alpha_G4DNAExcitation",
							   mod, reg, 0.0, hemax);

	mod = new G4DNADingfelderChargeDecreaseModel();
	em_config->SetExtraEmModel("alpha", "alpha_G4DNAChargeDecrease",
							   mod, reg, 0.0, hemax);

	mod = new G4DNAIonElasticModel();
	em_config->SetExtraEmModel("alpha", "alpha_G4DNAElastic",
							   mod, reg, 0.0, elimel);

	// ---
	// alpha1
	if (elimel < emax)
	{
		if (a1msc)
		{
			G4UrbanMscModel* msc = new G4UrbanMscModel();
			msc->SetActivationLowEnergyLimit(elimel);
			em_config->SetExtraEmModel("alpha+", "msc", msc, reg, 0.0, emax);
		}
		else
		{
			mod = new G4IonCoulombScatteringModel();
			mod->SetActivationLowEnergyLimit(elimel);
			em_config->SetExtraEmModel("alpha+", "CoulombScat", mod, reg, 0.0, emax);
		}
	}

	mod = new G4BraggIonModel();
	mod->SetActivationLowEnergyLimit(hemax / massRatio);
	em_config->SetExtraEmModel("alpha+", "hIoni",
							   mod, reg, 0.0, pminbba,
							   new G4IonFluctuations());

	mod = new G4BetheBlochModel();
	mod->SetActivationLowEnergyLimit(hemax / massRatio);
	em_config->SetExtraEmModel("alpha+", "hIoni",
							   mod, reg, pminbba, emax,
							   new G4IonFluctuations());

	mod = new G4DNARuddIonisationModel();
	em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAIonisation",
							   mod, reg, 0.0, hemax);

	mod = new G4DNAMillerGreenExcitationModel();
	em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAExcitation",
							   mod, reg, 0.0, hemax);

	mod = new G4DNADingfelderChargeDecreaseModel();
	em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAChargeDecrease",
							   mod, reg, 0.0, hemax);

	mod = new G4DNADingfelderChargeIncreaseModel();
	em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAChargeIncrease",
							   mod, reg, 0.0, hemax);

	mod = new G4DNAIonElasticModel();
	em_config->SetExtraEmModel("alpha+", "alpha+_G4DNAElastic",
							   mod, reg, 0.0, elimel);

	// ---
	// helium
	mod = new G4DNARuddIonisationModel();
	em_config->SetExtraEmModel("helium", "helium_G4DNAIonisation",
							   mod, reg, 0.0, hemax);

	mod = new G4DNAMillerGreenExcitationModel();
	em_config->SetExtraEmModel("helium", "helium_G4DNAExcitation",
							   mod, reg, 0.0, hemax);

	mod = new G4DNADingfelderChargeIncreaseModel();
	em_config->SetExtraEmModel("helium", "helium_G4DNAChargeIncrease",
							   mod, reg, 0.0, hemax);

	mod = new G4DNAIonElasticModel();
	em_config->SetExtraEmModel("helium", "helium_G4DNAElastic",
							   mod, reg, 0.0, elimel);
}

void TsEmDNAPhysicsActivator::DeactivateNuclearStopping_G4EmDNA(G4ProcessManager* pman,
																G4double elimel)
{
	G4ProcessVector* pv = pman->GetProcessList();
	G4int nproc = pman->GetProcessListLength();
	for (G4int i = 0; i < nproc; ++i)
	{
		if (((*pv)[i])->GetProcessSubType() == fNuclearStopping)
		{
			G4VEmProcess* proc = static_cast<G4VEmProcess*>((*pv)[i]);
			if (proc)
			{
				G4VEmModel* mod = new G4ICRU49NuclearStoppingModel();
				mod->SetActivationLowEnergyLimit(elimel);
				proc->SetEmModel(mod);
			}
			break;
		}
	}
}

G4bool TsEmDNAPhysicsActivator::HasMsc(G4ProcessManager* pman) const
{
	G4bool res = false;
	G4ProcessVector* pv = pman->GetProcessList();
	G4int nproc = pman->GetProcessListLength();
	for (G4int i = 0; i < nproc; ++i)
	{
		if (((*pv)[i])->GetProcessSubType() == fMultipleScattering)
		{
			res = true;
			break;
		}
	}
	return res;
}

// void TsEmDNAPhysicsActivator::FindOrAddProcess(const G4ParticleDefinition* particle, const G4String& proc_name)
//{
//	G4ProcessManager* pm = particle->GetProcessManager();
//	G4ProcessVector* pv = pm->GetProcessList();
//	G4int nproc = pm->GetProcessListLength();
//	for (G4int i = 0; i < nproc; ++i)
//	{
//		if (((*pv)[i])->GetProcessName() == proc_name)
//		{
//			return;
//		}
//	}
//	if (proc_name == "CoulombScat")
//	{
//		G4CoulombScattering* cs = new G4CoulombScattering();
//		cs->SetEmModel(new G4DummyModel());
//		pm->AddDiscreteProcess(cs);
//	}
//	else if (proc_name == "Rayl")
//	{
//		G4RayleighScattering* rs = new G4RayleighScattering();
//		rs->SetEmModel(new G4DummyModel());
//		pm->AddDiscreteProcess(rs);
//	}
// }

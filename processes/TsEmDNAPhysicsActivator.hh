#ifndef TsEmDNAPhysicsActivator_hh
#define TsEmDNAPhysicsActivator_hh 1

#include "G4VPhysicsConstructor.hh"

class TsParameterManager;

class G4EmParameters;
class G4ProcessManager;

class TsEmDNAPhysicsActivator : public G4VPhysicsConstructor
{
public:
	explicit TsEmDNAPhysicsActivator(TsParameterManager*);
	explicit TsEmDNAPhysicsActivator(G4int ver = 2);

	~TsEmDNAPhysicsActivator() override = default;

	void ConstructParticle() override;
	void ConstructProcess() override;

protected:
	void FindOrAddFindOrAddProcess(const G4ParticleDefinition* part,
								   const G4String& name);

	void AddDummyModels();

	void AddModels_TsEmDNA(const G4String&, const G4String&);
	void AddModels_G4EmDNA(const G4String&, const G4String&);

private:
	TsParameterManager* fPm;

	G4EmParameters* fEmParams;

	void AddElectronModels_TsEmDNA(const G4String&, const G4String&);
	void AddProtonModels_TsEmDNA(const G4String&, const G4String&);
	void AddIonModels_TsEmDNA(const G4String&, const G4String&);

	void AddPhotonModels_Livermore(const G4String&);
	// void AddPhotonModels_Penelope(const G4String &);
	void AddPositronModels_TsEmDNA(const G4String&);

	// Unchanged from Geant4.10.07.p03
	void AddElectronModels_G4EmDNA_opt0(const G4String& region, G4bool emsc, G4double elimel);
	void AddElectronModels_G4EmDNA_opt2(const G4String& region, G4bool emsc, G4double elimel);
	void AddElectronModels_G4EmDNA_opt4(const G4String& region, G4bool emsc, G4double elimel);
	void AddElectronModels_G4EmDNA_opt4a(const G4String& region, G4bool emsc, G4double elimel);
	void AddElectronModels_G4EmDNA_opt6(const G4String& region, G4bool emsc, G4double elimel);
	void AddElectronModels_G4EmDNA_opt6a(const G4String& region, G4bool emsc, G4double elimel);
	void AddElectronModels_G4EmDNA_opt7(const G4String& region, G4bool emsc, G4double elimel);
	void AddProtonModels_G4EmDNA(const G4String& region, G4bool pmsc, G4double elimel, G4double pminbb, G4double pmax);
	void AddHeliumModels_G4EmDNA(const G4String& region, G4bool a1msc, G4bool a2msc, G4double elimel, G4double pminbb, G4double pmax);
	void AddGenericIonModels_G4EmDNA(const G4String& region, G4double pminbb);
	void DeactivateNuclearStopping_G4EmDNA(G4ProcessManager*, G4double elimel);

	G4bool HasMsc(G4ProcessManager*) const;

	const G4double fEmMaxElectron = 1 * CLHEP::MeV;
	const G4double fEmMaxProton = 100 * CLHEP::MeV;
};

#endif

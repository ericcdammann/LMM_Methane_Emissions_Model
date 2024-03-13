function methane = carbon_to_methane(carbon)

MMC = 12.011;
MMH = 1.008;
Avo = 6.02214076*10^23;
carbon_grams = carbon.*1000;
carbon_atoms = carbon_grams./MMC.*Avo;
hydrogen_atoms = 4*carbon_atoms;
hydrogen_moles = hydrogen_atoms./Avo;
hydrogen = hydrogen_moles.*MMH./1000;
methane = carbon+hydrogen;

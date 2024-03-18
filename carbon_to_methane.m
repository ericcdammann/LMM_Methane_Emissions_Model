function methane = carbon_to_methane(carbon)

MMC = 12.011;
MMH = 1.008;
scale = 1+(4*MMH/MMC);
methane = carbon.*scale;
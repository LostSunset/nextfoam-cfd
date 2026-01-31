// ------------------------------------------------------------------------- //

const word cloudName
(
    args.getOrDefault<word>("cloud", "cloud")
);

Info<< "Constructing " << cloudName << endl;

SLGThermo slgThermo(mesh, thermo);

basicThermoCloud parcels
(
    cloudName,
    rho,
    U,
    g,
    slgThermo
);

// ------------------------------------------------------------------------- //

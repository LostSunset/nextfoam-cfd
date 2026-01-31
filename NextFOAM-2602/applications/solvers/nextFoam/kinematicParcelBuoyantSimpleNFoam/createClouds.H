// ------------------------------------------------------------------------- //

const word cloudName
(
    args.getOrDefault<word>("cloud", "cloud")
);

Info<< "Constructing " << cloudName << endl;

basicKinematicCloud parcels
(
    cloudName,
    rho,
    U,
    thermo.mu(),
    g
);

// ------------------------------------------------------------------------- //

/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION
    PtrList<SLGThermo> slgThermoFluid(fluidRegions.size());

    PtrList<autoPtr<regionModels::surfaceFilmModel>> surfaceFilmFluid
    (
        fluidRegions.size()
    );
    List<bool> surfaceFilmEnabledFluid(fluidRegions.size(), false);

    forAll(fluidRegions, i)
    {
        slgThermoFluid.set(i, new SLGThermo(fluidRegions[i], thermoFluid[i]));

        Info<< "    Adding surface film model" << endl;
        surfaceFilmFluid.set
        (
            i,
            new autoPtr<regionModels::surfaceFilmModel>
            (
                regionModels::surfaceFilmModel::New(fluidRegions[i], g)
            )
        );

        #ifdef NO_POST_PROCESS
        surfaceFilmEnabledFluid[i] = surfaceFilmFluid[i]->surfaceFilmActive();
        #endif
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else
    Info<< "\nConstructing surface film model" << endl;

    SLGThermo slgThermo(mesh, thermo);

    autoPtr<regionModels::surfaceFilmModel> surfaceFilm
    (
        regionModels::surfaceFilmModel::New(mesh, g)
    );

    #ifdef NO_POST_PROCESS
    bool surfaceFilmEnabled(surfaceFilm->surfaceFilmActive());
    #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


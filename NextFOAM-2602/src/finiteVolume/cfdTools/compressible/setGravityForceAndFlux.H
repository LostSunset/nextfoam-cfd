bool gravityEnabled(mag(g.value()) > SMALL ? true : false);

/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<const dimensionedScalar> initialMassFluid(fluidRegions.size());

    PtrList<dimensionedScalar> totalVolumeFluid(fluidRegions.size());
    PtrList<dimensionedScalar> rho0Fluid(fluidRegions.size());
    List<bool> rho0SpecifiedFluid(fluidRegions.size(), false);

    PtrList<volScalarField> hydroStaticPressureFluid(fluidRegions.size());
    #if defined(PIMPLE_CONTROL)
    PtrList<volScalarField> hydroStaticDensityFluid(fluidRegions.size());
    #endif
    PtrList<surfaceScalarField> gravityFluxPotentialFluid(fluidRegions.size());
    PtrList<volVectorField> gravityForceFluid(fluidRegions.size());
    PtrList<surfaceScalarField> gravityFluxFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        uniformDimensionedScalarField initialMassIO
        (
            IOobject
            (
                "initialMass",
                runTime.timeName(),
                "uniform",
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::domainIntegrate(rhoFluid[i])
        );

        initialMassFluid.set
        (
            i,
            new const dimensionedScalar
            (
                "initialMass",
                dimMass,
                initialMassIO.value()
            )
        );

        hydroStaticPressureFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "hydroStaticPressure", 
                    runTime.timeName(), 
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar("hydroStaticPressure", dimPressure, Zero)
            )
        );

        #if defined(PIMPLE_CONTROL)
        hydroStaticDensityFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "hydroStaticDensity", 
                    runTime.timeName(), 
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar("hydroStaticDensity", dimDensity, Zero)
            )
        );
        #endif

        gravityFluxPotentialFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "gravityFluxPotential",
                    runTime.timeName(),
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar
                (
                    "gravityFluxPotential",
                    dimForce/dimLength,
                    Zero
                )
            )
        );

        gravityForceFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "gravityForce",
                    runTime.timeName(),
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedVector
                (
                    "gravityForce",
                    dimDensity*dimAcceleration,
                    Zero
                )
            )
        );

        gravityFluxFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "gravityFlux",
                    runTime.timeName(),
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar
                (
                    "gravityFlux",
                    dimDensity*dimVelocity*dimArea,
                    Zero
                )
            )
        );

        rho0Fluid.set
        (
            i,
            new dimensionedScalar("rho0", dimDensity, Zero)
        );

        rho0SpecifiedFluid[i] = 
            operatingConditionsFluid[i].readIfPresent("rho0", rho0Fluid[i]);

        totalVolumeFluid.set
        (
            i,
            new dimensionedScalar
            (
                "totalVolume", 
                dimVolume,
                (gravityEnabled && !rho0SpecifiedFluid[i])
              ? gSum(fluidRegions[i].V())
              : Zero
            )
        );

        if (gravityEnabled)
        {
            dimensionedScalar& totalVolume = totalVolumeFluid[i];
            dimensionedScalar& rho0 = rho0Fluid[i];
            volScalarField& rho = rhoFluid[i];
            volScalarField& gh = ghFluid[i];
            volScalarField& hydroStaticPressure = hydroStaticPressureFluid[i];
            surfaceScalarField& rhof = rhofFluid[i]; 
            surfaceScalarField& gravityFluxPotential = 
                gravityFluxPotentialFluid[i];
            volVectorField& gravityForce = gravityForceFluid[i];
            #if defined(PIMPLE_CONTROL)
            volScalarField& hydroStaticDensity = hydroStaticDensityFluid[i]; 
            const volScalarField& psi = thermoFluid[i].psi();
            const dimensionedScalar& operatingPressure = 
                operatingPressureFluid[i];
            #endif

            bool& rho0Specified = rho0SpecifiedFluid[i];
            const fvMesh& mesh = fluidRegions[i];

            #include "updateGravitationalForce.H"
        }
        else
        {
            hydroStaticPressureFluid[i].clear();
            gravityFluxPotentialFluid[i].clear();
            gravityForceFluid[i].clear();
            gravityFluxFluid[i].clear();

            #if defined(PIMPLE_CONTROL)
            hydroStaticDensityFluid[i].clear();
            #endif
        }
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    uniformDimensionedScalarField initialMassIO
    (
        IOobject
        (
            "initialMass",
            runTime.timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::domainIntegrate(rho)
    );

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassIO.value()
    );

    tmp<volScalarField> tHydroStaticPressure
    (
        new volScalarField
        (
            IOobject("hydroStaticPressure", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("hydroStaticPressure", dimPressure, Zero)
        )
    );

    volScalarField& hydroStaticPressure = tHydroStaticPressure.ref();

    #if defined(PIMPLE_CONTROL)
    tmp<volScalarField> tHydroStaticDensity
    (
        new volScalarField
        (
            IOobject("hydroStaticDensity", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("hydroStaticDensity", dimDensity, Zero)
        )
    );
    volScalarField& hydroStaticDensity = tHydroStaticDensity.ref();
    #endif

    tmp<surfaceScalarField> tGravityFluxPotential
    (
        new surfaceScalarField
        (
            IOobject("gravityFluxPotential", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("gravityFluxPotential", dimForce/dimLength, Zero)
        )
    );
    surfaceScalarField& gravityFluxPotential = tGravityFluxPotential.ref();

    tmp<volVectorField> tGravityForce
    (
        new volVectorField
        (
            IOobject("gravityForce", runTime.timeName(), mesh),
            mesh,
            dimensionedVector("gravityForce", dimDensity*dimAcceleration, Zero)
        )
    );
    volVectorField& gravityForce = tGravityForce.ref();

    tmp<surfaceScalarField> tGravityFlux
    (
        new surfaceScalarField
        (
            IOobject("gravityFlux", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar
            (
                "gravityFlux", 
                dimDensity*dimVelocity*dimArea, 
                Zero
            )
        )
    );
    #ifdef NO_POST_PROCESS
    surfaceScalarField& gravityFlux = tGravityFlux.ref();
    #endif

    dimensionedScalar rho0("rho0", dimDensity, Zero);
    bool rho0Specified(operatingConditions.readIfPresent("rho0", rho0));

    dimensionedScalar totalVolume
    (
        "totalVolume",
        dimVolume,
        (gravityEnabled && !rho0Specified) ? gSum(mesh.V()) : Zero
    );

    if (gravityEnabled)
    {
        #include "updateGravitationalForce.H"
    }
    else
    {
        tHydroStaticPressure.clear();
        tGravityFluxPotential.clear();
        tGravityForce.clear();
        tGravityFlux.clear();

        #if defined(PIMPLE_CONTROL)
        tHydroStaticDensity.clear();
        #endif
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    const label nFluid(fluidRegions.size());
    PtrList<rhoReactionThermo> thermoFluid(nFluid);
    PtrList<volScalarField> rhoFluid(nFluid);
    PtrList<dimensionedScalar> totalVolumeFluid(nFluid);
    PtrList<volVectorField> UFluid(nFluid);
    PtrList<volScalarField> UmagFluid(nFluid);
    PtrList<surfaceScalarField> phiFluid(nFluid);
    PtrList<surfaceScalarField> rhofFluid(nFluid);
    #if defined(PIMPLE_CONTROL)
    PtrList<pressureControl> pressureControlFluid(nFluid);
    PtrList<volScalarField> KFluid(nFluid);
    #endif
    PtrList<autoPtr<compressible::turbulenceModel>> turbulenceFluid(nFluid);
    PtrList<volScalarField> p_rghFluid(nFluid);

    List<label> pRefCellFluid(fluidRegions.size(), 0);
    List<scalar> pRefValueFluid(fluidRegions.size(), Zero);

    PtrList<uniformDimensionedScalarField> initialMeanTemperatureFluid(nFluid);

    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermoFluid\n" << endl;

        thermoFluid.set
        (
            i,
            rhoReactionThermo::New(fluidRegions[i]).ptr()
        );

        thermoFluid[i].validate(args.executable(), "h", "e");

        Info<< "    Adding to rhoFluid\n" << endl;
        rhoFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                thermoFluid[i].rho()
            )
        );

        totalVolumeFluid.set
        (       
            i,  
            new dimensionedScalar
            (   
                "totalVolume",
                dimVolume,
                gSum(fluidRegions[i].V())
            )
        );  

        initialMeanTemperatureFluid.set
        (
            i,
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "initialMeanTemperature",
                    runTime.timeName(),
                    "uniform",
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::domainIntegrate(thermoFluid[i].T())/totalVolumeFluid[i]
            )
        );

        Info<< "    Adding to UFluid\n" << endl;

        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to UmagFluid\n" << endl;
        UmagFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Umag",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mag(UFluid[i])
            )
        );

        #include "compressibleInitPhi.H"   

        Info<< "    Adding to rhofFluid\n" << endl;
        rhofFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "faceRho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(rhoFluid[i])
            )
        );

        #if defined(PIMPLE_CONTROL)
        pressureControlFluid.set
        (
            i,
            new pressureControl
            (
                thermoFluid[i].p(), 
                rhoFluid[i],
                fluidRegions[i].solutionDict().subDict("PIMPLE"),
                false
            )
        );
        #endif

        Info<< "    Adding to turbulenceFluid\n" << endl;
        turbulenceFluid.set
        (
            i,
            new autoPtr<compressible::turbulenceModel>
            (
                compressible::turbulenceModel::New
                (
                    rhoFluid[i],
                    UFluid[i],
                    phiFluid[i],
                    thermoFluid[i]
                )
            )
        );

        p_rghFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p_rgh",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

        #if defined(SIMPLE_CONTROL)
        setRefCell
        (
            thermoFluid[i].p(),
            p_rghFluid[i],
            //fluidRegions[i].solutionDict().subDict("SIMPLE"),
            pvCoupling.dict(fluidRegions.get(i)),
            pRefCellFluid[i],
            pRefValueFluid[i]
        );
        #else
        if (thermoFluid[i].incompressible())
        {
            setRefCell
            (
                thermoFluid[i].p(),
                p_rghFluid[i],
                //fluidRegions[i].solutionDict().subDict("PIMPLE"),
                pvCoupling.dict(fluidRegions.get(i)),
                pRefCellFluid[i],
                pRefValueFluid[i]
            );
        }

        // createK
        Info<< "    Adding to KFluid\n" << endl;
        KFluid.set
        (
            i,
            new volScalarField("K", 0.5*magSqr(UFluid[i]))
        );

        if (UFluid[i].nOldTimes())
        {
            volVectorField* Uold = &UFluid[i].oldTime();
            volScalarField* Kold = &KFluid[i].oldTime();
            *Kold == 0.5*magSqr(*Uold);

            while (Uold->nOldTimes())
            {
                Uold = &Uold->oldTime();
                Kold = &Kold->oldTime();
                *Kold == 0.5*magSqr(*Uold);
            }
        }
        #endif
    }

    const label nSolid(solidRegions.size());
    PtrList<coordinateSystem> coordinates(nSolid);
    PtrList<solidThermo> thermoSolid(nSolid);
    PtrList<volScalarField> betavSolid(nSolid);
    PtrList<volSymmTensorField> aniAlphas(nSolid);

    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermoSolid\n" << endl;
        thermoSolid.set(i, solidThermo::New(solidRegions[i]));

        if (!thermoSolid[i].isotropic())
        {
            Info<< "    Adding coordinateSystems\n" << endl;
            coordinates.set
            (
                i,
                coordinateSystem::New
                (
                    solidRegions[i],
                    thermoSolid[i],
                    coordinateSystem::typeName
                )
            );

            tmp<volVectorField> tkappaByCp =
                thermoSolid[i].Kappa()/thermoSolid[i].Cp();

            aniAlphas.set
            (
                i,
                new volSymmTensorField
                (
                    IOobject
                    (
                        "Anialpha",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedSymmTensor(tkappaByCp().dimensions(), Zero),
                    zeroGradientFvPatchSymmTensorField::typeName
                )
            );

            aniAlphas[i].primitiveFieldRef() =
                coordinates[i].transformPrincipal
                (
                    solidRegions[i].cellCentres(),
                    tkappaByCp()
                );
            aniAlphas[i].correctBoundaryConditions();
        }

        IOobject betavSolidIO
        (
            "betavSolid",
            runTime.timeName(),
            solidRegions[i],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (betavSolidIO.typeHeaderOk<volScalarField>(true))
        {
            betavSolid.set
            (
                i,
                new volScalarField(betavSolidIO, solidRegions[i])
            );
        }
        else
        {
            betavSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "betavSolid",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("1", dimless, scalar(1))
                )
            );
        }
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
    rhoReactionThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();

    const volScalarField& psi = thermo.psi();

    dimensionedScalar totalVolume("totalVolume", dimVolume, gSum(mesh.V()));

    uniformDimensionedScalarField initialMeanTemperature
    (
        IOobject
        (
            "initialMeanTemperature",
            runTime.timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::domainIntegrate(thermo.T())/totalVolume
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField Umag
    (
        IOobject
        (
            "Umag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(U)
    );

    #include "compressibleInitPhi.H"   

    surfaceScalarField rhof
    (
        IOobject
        (
            "faceRho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(rho)
    );

    #if defined(PIMPLE_CONTROL)
    pressureControl pressureControl(p, rho, pvCoupling.dict(), false);
    #endif

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    mesh.setFluxRequired(p_rgh.name());

    label pRefCell = 0;
    scalar pRefValue = 0.0;

    #if defined(SIMPLE_CONTROL)
    setRefCell
    (
        p,
        p_rgh,
        pvCoupling.dict(),
        pRefCell,
        pRefValue
    );

    scalar UUrf = mesh.equationRelaxationFactor(U.name());
    scalar pUrf = mesh.fieldRelaxationFactor(p_rgh.name());
    scalar hUrf = mesh.equationRelaxationFactor(thermo.he().name());
    scalar rhofUrf = sqr(UUrf)*pUrf*hUrf;
    #else
    if (thermo.incompressible())
    {
        setRefCell // only get pRefCell & pRefValue from dictionary
        (
            p,
            p_rgh,
            pvCoupling.dict(),
            pRefCell,
            pRefValue
        );
    }

    #include "createK.H"
    #endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

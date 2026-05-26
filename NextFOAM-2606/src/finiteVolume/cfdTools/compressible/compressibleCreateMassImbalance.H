#define MASSIMBALANCE

/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<volScalarField> massImbalanceFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        massImbalanceFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "massImbalance",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::surfaceIntegrate(phiFluid[i])
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    volScalarField massImbalance
    (
        IOobject
        (
            "massImbalance",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::surfaceIntegrate(phi)
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

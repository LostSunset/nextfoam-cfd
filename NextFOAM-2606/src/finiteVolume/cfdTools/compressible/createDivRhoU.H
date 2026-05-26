/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<autoPtr<volScalarField>> divrhoUFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        divrhoUFluid.set(i, new autoPtr<volScalarField>);

        if (correctPhiFluid[i])
        {
            divrhoUFluid[i].reset
            (
                new volScalarField("divrhoU", fvc::div(phiFluid[i]))
            );
        }
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    autoPtr<volScalarField> divrhoU;
    if (correctPhi)
    {
        divrhoU.reset
        (
            new volScalarField("divrhoU", fvc::div(phi))
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

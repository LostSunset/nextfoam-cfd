/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    List<scalar> pCorrLimitFluid(fluidRegions.size());
    PtrList<dimensionedScalar> pCorrMaxFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        pCorrLimitFluid[i] =
            pvCoupling.dict(fluidRegions.get(i)).getOrDefault<scalar>
            (
                "pCorrLimit", 
                GREAT
            );

        pCorrMaxFluid.set
        (
            i,
            new dimensionedScalar
            (
                pCorrLimitFluid[i]*max(thermoFluid[i].p().internalField())
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    scalar pCorrLimit
    (
        pvCoupling.dict().getOrDefault<scalar>("pCorrLimit", GREAT)
    );

    dimensionedScalar pCorrMax(pCorrLimit*max(p.internalField()));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

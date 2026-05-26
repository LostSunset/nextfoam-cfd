/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

{
    forAll(fluidRegions, i)
    {
        volScalarField& rDeltaT = trDeltaTFluid[i].ref();

        scalar rDeltaTSmoothingCoeff
        (
            runTime.controlDict().getOrDefault<scalar>
            (
                "rDeltaTSmoothingCoeff",
                0.01
            )
        );

        scalar rDeltaTDampingCoeff
        (
            runTime.controlDict().getOrDefault<scalar>
            (
                "rDeltaTDampingCoeff", 
                1.0
            )
        );

        volScalarField rDeltaT0("rDeltaT0", rDeltaT);

        // Set the reciprocal time-step from the local Courant number
        rDeltaT.ref() = max
        (
            1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
            fvc::surfaceSum(mag(phiFluid[i]))()()
            /((2*maxCo)*fluidRegions[i].V()*rhoFluid[i]())
        );

        // Update the boundary values of the reciprocal time-step
        rDeltaT.correctBoundaryConditions();

        Info<< "Flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

        if (rDeltaTSmoothingCoeff < 1.0)
        {
            fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
        }

        Info<< "Smoothed flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

        // Limit rate of change of time scale
        // - reduce as much as required
        // - only increase at a fraction of old time scale
        if
        (
            rDeltaTDampingCoeff < 1.0
         && runTime.timeIndex() > runTime.startTimeIndex() + 1
        )
        {
            rDeltaT =
                rDeltaT0
               *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

            Info<< "Damped flow time scale min/max = "
                << gMin(1/rDeltaT.primitiveField())
                << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

{
    volScalarField& rDeltaT = trDeltaT.ref();

    scalar rDeltaTSmoothingCoeff
    (
        runTime.controlDict().getOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.01
        )
    );

    scalar rDeltaTDampingCoeff
    (
        runTime.controlDict().getOrDefault<scalar>
        (   
            "rDeltaTDampingCoeff", 
            1.0
        )
    );

    // Maximum change in cell temperature per iteration
    // (relative to previous value)
    scalar alphaTemp(pvCoupling.dict().getOrDefault("alphaTemp", 0.05));

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(mag(phi))()()/((2*maxCo)*mesh.V()*rho())
    );

    // Reaction source time scale
    {
        volScalarField::Internal rDeltaTT
        (
            mag
            (
               parcels.hsTrans()/(mesh.V()*runTime.deltaT())
            )
           /(
               alphaTemp
              *rho()
              *thermo.Cp()()()
              *thermo.T()
           )
        );

        Info<< "    Temperature = "
            << gMin(1/(rDeltaTT.field() + VSMALL)) << ", "
            << gMax(1/(rDeltaTT.field() + VSMALL)) << endl;

        rDeltaT.ref() = max
        (
            rDeltaT(),
            rDeltaTT
        );
    }

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT =
            rDeltaT0
           *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

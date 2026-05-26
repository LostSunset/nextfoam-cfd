#if defined(SIMPLE_CONTROL)
    bool limitDeltaT
    (
        pvCoupling.dict().lookupOrDefault<bool>("limitDeltaT", true)
    );

    if (limitDeltaT)
    {
        scalar maxDeltaT
        (
            pvCoupling.dict().lookupOrDefault<scalar>("maxDeltaT", 15)
        );

        scalar localCorrection
        (
            pvCoupling.dict().lookupOrDefault<scalar>("localCorrection", 2e-2)
        );

        scalar Tmin
        (
            pvCoupling.dict().lookupOrDefault<scalar>("Tmin", 2.17)
        );

        thermo.limitedCorrect(maxDeltaT, localCorrection, Tmin);
    }
    else
    {
        thermo.correct();
    }
#else
    thermo.correct();
#endif

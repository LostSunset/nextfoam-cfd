{
    const dimensionedScalar rhoMax
    (
        "rhoMax",
        dimDensity,
        GREAT,
        pvCoupling.dict()
    );

    const dimensionedScalar rhoMin
    (
        "rhoMin",
        dimDensity,
        SMALL,
        pvCoupling.dict()
    );

    thermo.limitRho(rhoMin, rhoMax);

#if defined(PIMPLE_CONTROL)
    #include "compressibleContinuityErrs.H"

    rho = thermo.rho();

    rho.relax();

    rhof = fvc::interpolate(rho);
#else
    rho = thermo.rho();

    rho.relax();

    // A temporary, unjustified way to stabilize pressure-velocity coupling. 
    //
    // The interpolated face density is used as a lagged coupling field in the
    // pressure-density, momentum and energy coupling.  It is therefore not
    // updated as a purely bounded reconstruction of the current cell-centred
    // density.
    // Its automatic relaxation factor is estimated from the effective
    // thermo-mechanical relaxation of the coupled system.  The velocity
    // relaxation is applied quadratically to account for the U^2-type
    // mechanical-energy/dynamic-pressure coupling, while the pressure and
    // enthalpy relaxation factors represent the pressure-density and thermal
    // energy coupling scales.
    // Needs to be implemented as a method with a good rationale. Maybe...
    {
        rhof.storePrevIter();

        UUrf = mesh.equationRelaxationFactor(U.name());
        pUrf = mesh.fieldRelaxationFactor(p_rgh.name());
        hUrf = mesh.equationRelaxationFactor(thermo.he().name());
        rhofUrf = sqr(UUrf)*pUrf*hUrf;

        rhof = fvc::interpolate(rho);
        rhof.relax(rhofUrf);
    }
#endif

    if (gravityEnabled)
    {
        #include "updateGravitationalForce.H"
    }

    //----------------------------------------------------------------------//
    // Treatment to maintain initial mass when simulating closed domain
    // problems. ---> What if the mass source exists??
    bool maintainInitialMass
    (
        pvCoupling.dict().lookupOrDefault<bool>("maintainInitialMass", false)
    );

    if (p_rgh.needReference() && maintainInitialMass)
    {
        // Changing the density doesn't really affect the mass conservation
        //, but explicitly changing the pressure solution never causes 
        // the pressure equation to converge.

        dimensionedScalar compressibility(fvc::domainIntegrate(thermo.psi()));

        if (compressibility.value() > SMALL)
        {
            dimensionedScalar totalMass(fvc::domainIntegrate(rho));
            //rho += psi*(initialMass - totalMass)/compressibility;
            p += (initialMass - totalMass)/compressibility;
            p_rgh = p - operatingPressure - hydroStaticPressure;
            p_rgh.correctBoundaryConditions();

            Info << "Mass Imbalance = " 
                << (initialMass - totalMass).value()
                << nl << endl;
        }
    }
    //----------------------------------------------------------------------//
}

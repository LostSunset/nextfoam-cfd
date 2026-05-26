if (!rho0Specified)
{
    if (mesh.moving())
    {
        totalVolume.value() = gSum(mesh.V());
    }

    rho0 = fvc::domainIntegrate(rho)/totalVolume;
}

hydroStaticPressure = rho0*gh;
//gravityFluxPotential = fvc::interpolate(rho0 - rho)*g & mesh.Sf();
gravityFluxPotential = (rho0 - rhof)*g & mesh.Sf();
gravityForce = (rho0 - rho)*g;

//hydroStaticPressure = rho*gh;
//gravityFluxPotential = ghf*fvc::snGrad(rho)*mesh.magSf();
//gravityForce = fvc::reconstruct(gravityFluxPotential);
//gravityForce = gh*fvc::grad(rho);

#if defined(PIMPLE_CONTROL)
hydroStaticDensity = psi*(hydroStaticPressure + operatingPressure);
#endif

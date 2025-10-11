if (MRFEnabled)
{
    MRF.updateForce(rhof, phi, U);
}

volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

if (pvCoupling.corrPISO() == pvCoupling.nCorrPISO())
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    rhof*(momentumInterpolate(HbyA) & mesh.Sf())
  + fvc::relaxCorr(rhof, phi, U, p, U.select(pvCoupling.finalIter()))
  + fvc::ddtCorr(phi, rhof, rho.name(), U, rAU, p_rgh)
);

#include "compressibleSetConsistentFlux.H"

// Add Coriolis flux and make it relative to MRF
if (MRFEnabled)
{
    MRF.addCoriolisFlux(rhof, rhorAUf, phiHbyA);
}

//- phiHbyA should be relative to the mesh motion and should not contain phig 
//- before adjustment
//fvc::makeRelative(phiHbyA, rhof, U);
//adjustPhi(phiHbyA, U, p_rgh);

if (gravityEnabled)
{
    phiHbyA += constrainPhig(gravityFlux, U);
}

//- Update the pressure BCs to ensure flux consistency
//- phiHbyA should be absolute to the mesh motion before constrainPressure
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

fvc::makeRelative(phiHbyA, rhof, rho.name(), U);

// Non-orthogonal pressure corrector loop
while (pvCoupling.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::ddt(psi, p_rgh)
      - fvm::laplacian(rhorAUf, p_rgh) 
     == 
        fvOptions(psi, p_rgh, rho.name())
      - fvc::div(phiHbyA)
    );

    if (gravityEnabled)
    {
        p_rghEqn += fvc::ddt(hydroStaticDensity);
    }

    if (thermo.incompressible())
    {
        p_rghEqn.setReference(pRefCell, pRefValue);
    }

    p_rghEqn.solve(mesh.solver(p_rgh.select(pvCoupling.finalInnerIter())));

    if (pvCoupling.finalNonOrthogonalIter())
    {
        // Correct mass flow rate
        phi = phiHbyA + p_rghEqn.flux();
    }
}

// Momentum corrector
#include "compressibleCorrectVelocity.H"

// Explicitly relax pressure for momentum corrector
p_rgh.relax();
p_rgh.correctBoundaryConditions();

#include "updateAbsolutePressure.H"

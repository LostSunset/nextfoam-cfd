volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

tUEqn.clear();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    rhof*(momentumInterpolate(HbyA) & mesh.Sf())
  + fvc::relaxCorr(rhof, phi, U, p, U.name())
);

#include "compressibleSetConsistentFlux.H"

// Make flux relative to MRF
if (MRFEnabled)
{
    MRF.makeRelative(rhof, phiHbyA);
}

adjustPhi(phiHbyA, U, p_rgh);

if (gravityEnabled)
{
    phiHbyA += constrainPhig(gravityFlux, U);
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

// Reset massImbalance.
massImbalance = dimensionedScalar(massImbalance.dimensions(), Zero);

// Non-orthogonal pressure corrector loop
while (pvCoupling.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rhorAUf, p_rgh)
     == 
        fvc::div(phiHbyA)
      - fvOptions(psi, p_rgh, rho.name()) 
    );

    p_rghEqn.setReference(pRefCell, pRefValue);

    p_rghEqn.solve();

    if (pvCoupling.finalNonOrthogonalIter())
    {   
        // Correct mass flow rate
        phi = phiHbyA - p_rghEqn.flux();
    }
}

// Momentum corrector
#include "compressibleCorrectVelocity.H"

// Complete continuity residual including mass source
massImbalance += fvc::surfaceIntegrate(phi);

// Explicitly relax pressure for momentum corrector
p_rgh.relax();
p_rgh.correctBoundaryConditions();

#include "updateAbsolutePressure.H"

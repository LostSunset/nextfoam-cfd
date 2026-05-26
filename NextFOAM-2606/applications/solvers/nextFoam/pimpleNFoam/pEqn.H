volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pvCoupling.corrPISO() == pvCoupling.nCorrPISO())
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (momentumInterpolate(HbyA) & mesh.Sf())
  + fvc::relaxCorr(phi, U, p, U.select(pvCoupling.finalIter()))
  + fvc::ddtCorr(phi, U, rAU, p)
);

#include "setConsistentFlux.H"

// Make flux relative to MRF
if (MRFEnabled)
{
    MRF.makeRelative(phiHbyA);
}

// if (p.needReference())
// {
//     fvc::makeRelative(phiHbyA, U);
//     adjustPhi(phiHbyA, U, p);
//     fvc::makeAbsolute(phiHbyA, U);
// }

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAUf, MRF);

fvc::makeRelative(phiHbyA, U);

// Non-orthogonal pressure corrector loop
while (pvCoupling.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pvCoupling.finalInnerIter())));

    if (pvCoupling.finalNonOrthogonalIter())
    {
        // Correct mass flow rate
        phi = phiHbyA - pEqn.flux();
    }
}

// Momentum corrector
#include "correctVelocity.H"

// Explicitly relax pressure for momentum corrector
p.relax();
p.correctBoundaryConditions();

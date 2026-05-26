#include "compressibleUpdatePressureGradient.H"

volVectorField velocityCorrection(gradp_rgh);

if (gravityEnabled)
{
    velocityCorrection += gravityForce;
}

velocityCorrection *= rAtU();

U = HbyA - velocityCorrection;
U.correctBoundaryConditions();

fvOptions.correct(U);

Umag = postMag(U);


gradp = fvc::gaussGrad(p);

volVectorField velocityCorrection(gradp);

velocityCorrection *= rAtU();

U = HbyA - velocityCorrection;
U.correctBoundaryConditions();

fvOptions.correct(U);

Umag = postMag(U);

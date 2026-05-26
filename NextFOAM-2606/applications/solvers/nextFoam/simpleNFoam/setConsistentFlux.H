if (pvCoupling.consistent())
{
    volVectorField velocityCorrection(gradp);
    surfaceScalarField fluxCorrection(fvc::snGrad(p)*mesh.magSf());

    velocityCorrection *= drAU();
    fluxCorrection *= momentumInterpolate(drAU());

    HbyA = constrainHbyA(HbyA + velocityCorrection, U, p);
    phiHbyA += constrainPhi(fluxCorrection, U, p);
}

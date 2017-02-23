volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA("HbyA", U);
HbyA = rAU*UEqn.H();
surfaceScalarField phiHbyA
(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh.Sf())
  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
);

MRF.makeRelative(phiHbyA);

adjustPhi(phiHbyA, U, p);

// Non-orthogonal pressure corrector loop
while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phiHbyA) + fvc::div(rAU*f)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

    if (piso.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

//Update pressure gradient field
gradP = fvc::grad(p);

//Update velocity
U = HbyA - rAU*gradP + rAU*f;
U.correctBoundaryConditions();

//Update IB forcing term
HFDIB.updateIBForcing( U, Ui, UEqn, f , &gradP);

#include "continuityErrs.H"
fvOptions.correct(U);

// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

//if (piso.momentumPredictor())
//{

    //Update IB forcing term
    //HFDIB.updateIBForcing( U, Ui, UEqn, f , &gradP);

    solve(UEqn == -gradP + f );

    fvOptions.correct(U);
//}

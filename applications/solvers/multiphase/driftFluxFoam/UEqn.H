    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + fvc::div
        (
            (alpha1/(scalar(1.001) - alpha1))*((rho2*rho1)/rho)*Vdj*Vdj,
            "div(phiVdj,Vdj)"
        )
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev2(T(fvc::grad(U))))
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }

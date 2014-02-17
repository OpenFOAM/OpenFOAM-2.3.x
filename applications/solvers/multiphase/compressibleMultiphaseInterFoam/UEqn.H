    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(multiphaseProperties.rhoPhi(), U)
      + turbulence->divDevRhoReff(U)
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    multiphaseProperties.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        K = 0.5*magSqr(U);
    }

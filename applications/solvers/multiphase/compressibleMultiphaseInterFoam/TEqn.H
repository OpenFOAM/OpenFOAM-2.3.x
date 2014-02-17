{
    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(multiphaseProperties.rhoPhi(), T)
      - fvm::laplacian(multiphaseProperties.alphaEff(turbulence->mut()), T)
      + (
            fvc::div(fvc::absolute(phi, U), p)
          + fvc::ddt(rho, K) + fvc::div(multiphaseProperties.rhoPhi(), K)
        )*multiphaseProperties.rCv()
    );

    TEqn.relax();
    TEqn.solve();

    multiphaseProperties.correct();
}

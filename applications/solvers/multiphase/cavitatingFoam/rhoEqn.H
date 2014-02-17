{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvm::div(phiv, rho)
    );

    rhoEqn.solve();

    phi = rhoEqn.flux();

    Info<< "max-min rho: " << max(rho).value()
        << " " << min(rho).value() << endl;

    rho == max(rho, rhoMin);
}

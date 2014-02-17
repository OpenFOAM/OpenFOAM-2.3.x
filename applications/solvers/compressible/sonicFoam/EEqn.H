{
    solve
    (
        fvm::ddt(rho, e) + fvm::div(phi, e)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + fvc::div(fvc::absolute(phi/fvc::interpolate(rho), U), p, "div(phiv,p)")
      - fvm::laplacian(turbulence->alphaEff(), e)
    );

    thermo.correct();
}

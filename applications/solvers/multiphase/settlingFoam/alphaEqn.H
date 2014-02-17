{
    surfaceScalarField phiAlpha
    (
        IOobject
        (
            "phiAlpha",
            runTime.timeName(),
            mesh
        ),
        phi + rhoc*(mesh.Sf() & fvc::interpolate(Vdj))
    );

    fvScalarMatrix AlphaEqn
    (
          fvm::ddt(rho, Alpha)
        + fvm::div(phiAlpha, Alpha)
        - fvm::laplacian(mut, Alpha)
    );

    AlphaEqn.relax();
    AlphaEqn.solve();

    Info<< "Solid phase fraction = "
        << Alpha.weightedAverage(mesh.V()).value()
        << "  Min(Alpha) = " << min(Alpha).value()
        << "  Max(Alpha) = " << max(Alpha).value()
        << endl;

    Alpha.min(1.0);
    Alpha.max(0.0);

    rho == rhoc/(scalar(1) + (rhoc/rhod - 1.0)*Alpha);
    alpha == rho*Alpha/rhod;
}

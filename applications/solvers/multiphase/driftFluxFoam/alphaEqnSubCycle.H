{
    surfaceScalarField phiAlpha
    (
        IOobject
        (
            "phiAlpha",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", phi.dimensions(), 0)
    );

    surfaceScalarField phir
    (
        rhoc*(mesh.Sf() & fvc::interpolate(Vdj/rho))
    );

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField phiAlphaSum
        (
            IOobject
            (
                "phiAlphaSum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("0", phi.dimensions(), 0)
        );

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            phiAlphaSum += (runTime.deltaT()/totalDeltaT)*phiAlpha;
        }

        phiAlpha = phiAlphaSum;
    }
    else
    {
        #include "alphaEqn.H"
    }

    // Apply the diffusion term separately to allow implicit solution
    // and boundedness of the explicit advection
    {
        fvScalarMatrix alphaEqn
        (
            fvm::ddt(alpha) - fvc::ddt(alpha)
          - fvm::laplacian(mut/rho, alpha)
        );

        alphaEqn.solve();

        phiAlpha += alphaEqn.flux();
    }

    Info<< "Phase-1 volume fraction = "
        << alpha.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha) = " << min(alpha).value()
        << "  Max(alpha) = " << max(alpha).value()
        << endl;

    rhoPhi = phiAlpha*(rhod - rhoc) + phi*rhoc;
    rho == alpha*rhod + (scalar(1) - alpha)*rhoc;
}

{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        phiAlpha =
        (
            fvc::flux
            (
                phi,
                alpha,
                alphaScheme
            )
          + fvc::flux
            (
                phir,
                alpha,
                alpharScheme
            )
        );

        MULES::explicitSolve(alpha, phi, phiAlpha, 1, 0);
    }

    Info<< "Phase-1 volume fraction = "
        << alpha.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha) = " << min(alpha).value()
        << "  Max(alpha) = " << max(alpha).value()
        << endl;
}

if (turbulence)
{
    if (mesh.changing())
    {
        y.correct();
    }

    dimensionedScalar k0("k0", k.dimensions(), 0);
    dimensionedScalar kMin("kMin", k.dimensions(), SMALL);
    dimensionedScalar epsilon0("epsilon0", epsilon.dimensions(), 0);
    dimensionedScalar epsilonMin("epsilonMin", epsilon.dimensions(), SMALL);

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(mut*(tgradU() && dev(twoSymm(tgradU()))));
    tgradU.clear();

    volScalarField Gcoef
    (
        Cmu*k/sigmak*(g & fvc::grad(rho))/(epsilon + epsilonMin)
    );

    volScalarField mul(twoPhaseProperties.mu());

    #include "wallFunctions.H"

    // Dissipation equation
    fvScalarMatrix epsEqn
    (
        fvm::ddt(rho, epsilon)
      + fvm::div(rhoPhi, epsilon)
      - fvm::laplacian
        (
            mut/sigmaEps + mul, epsilon,
            "laplacian(DepsilonEff,epsilon)"
        )
     ==
        C1*G*epsilon/(k + kMin)
      - fvm::SuSp(C1*(1.0 - C3)*Gcoef, epsilon)
      - fvm::Sp(C2*rho*epsilon/(k + kMin), epsilon)
    );

    #include "wallDissipation.H"

    epsEqn.relax();
    epsEqn.solve();

    bound(epsilon, epsilon0);


    // Turbulent kinetic energy equation
    fvScalarMatrix kEqn
    (
        fvm::ddt(rho, k)
      + fvm::div(rhoPhi, k)
      - fvm::laplacian
        (
            mut/sigmak + mul, k,
            "laplacian(DkEff,k)"
        )
     ==
        G
      - fvm::SuSp(Gcoef, k)
      - fvm::Sp(rho*epsilon/(k + kMin), k)
    );

    kEqn.relax();
    kEqn.solve();

    bound(k, k0);


    //- Re-calculate viscosity
    mut = rho*Cmu*sqr(k)/(epsilon + epsilonMin);

    #include "wallViscosity.H"

    muEff = mut + mul;
}
else
{
    muEff = mut + twoPhaseProperties.mu();
}

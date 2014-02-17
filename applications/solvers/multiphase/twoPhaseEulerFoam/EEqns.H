{
    volScalarField& he1 = thermo1.he();
    volScalarField& he2 = thermo2.he();

    volScalarField Cpv1("Cpv1", thermo1.Cpv());
    volScalarField Cpv2("Cpv2", thermo2.Cpv());

    volScalarField heatTransferCoeff(fluid.heatTransferCoeff());

    fvScalarMatrix he1Eqn
    (
        fvm::ddt(alpha1, he1) + fvm::div(alphaPhi1, he1)
      + fvc::ddt(alpha1, K1) + fvc::div(alphaPhi1, K1)

        // Compressibity correction
      - fvm::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), he1)
      - (fvc::ddt(alpha1) + fvc::div(alphaPhi1))*K1

      + (
            he1.name() == thermo1.phasePropertyName("e")
          ? fvc::ddt(alpha1)*p + fvc::div(alphaPhi1, p)
          : -alpha1*dpdt
        )/rho1

      //***HGW- fvm::laplacian(alpha1*turbulence1->alphaEff(), he1)
      - fvm::laplacian(alpha1*phase1.turbulence().nuEff(), he1)
     ==
        heatTransferCoeff*(thermo2.T() - thermo1.T())/rho1
      + heatTransferCoeff*he1/Cpv1/rho1
      - fvm::Sp(heatTransferCoeff/Cpv1/rho1, he1)
    );

    fvScalarMatrix he2Eqn
    (
        fvm::ddt(alpha2, he2) + fvm::div(alphaPhi2, he2)
      + fvc::ddt(alpha2, K2) + fvc::div(alphaPhi2, K2)

        // Compressibity correction
      - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), he2)
      - (fvc::ddt(alpha2) + fvc::div(alphaPhi2))*K2

      + (
            he2.name() == thermo2.phasePropertyName("e")
          ? fvc::ddt(alpha2)*p + fvc::div(alphaPhi2, p)
          : -alpha2*dpdt
        )/rho2

      //***HGW- fvm::laplacian(alpha2*turbulence2->alphaEff(), he2)
      - fvm::laplacian(alpha2*phase2.turbulence().nuEff(), he2)
     ==
        heatTransferCoeff*(thermo1.T() - thermo2.T())/rho2
      + heatTransferCoeff*he2/Cpv2/rho2
      - fvm::Sp(heatTransferCoeff/Cpv2/rho2, he2)
    );

    he1Eqn.relax();
    he1Eqn.solve();

    he2Eqn.relax();
    he2Eqn.solve();

    thermo1.correct();
    thermo2.correct();
}

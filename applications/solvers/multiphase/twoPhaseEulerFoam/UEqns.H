mrfZones.correctBoundaryVelocity(U1);
mrfZones.correctBoundaryVelocity(U2);
mrfZones.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

volScalarField dragCoeff(fluid.dragCoeff());

{
    volScalarField virtualMassCoeff(fluid.virtualMassCoeff());

    volVectorField liftForce(fluid.liftForce());
    volVectorField wallLubricationForce(fluid.wallLubricationForce());
    volVectorField turbulentDispersionForce(fluid.turbulentDispersionForce());

    {
        U1Eqn =
        (
            fvm::ddt(alpha1, rho1, U1)
          + fvm::div(alphaRhoPhi1, U1)
          - fvm::Sp(fvc::ddt(alpha1, rho1) + fvc::div(alphaRhoPhi1), U1)
          + phase1.turbulence().divDevRhoReff(U1)
         ==
          - fvm::Sp(dragCoeff, U1)
          - (
                liftForce
              + wallLubricationForce
              + turbulentDispersionForce
            )
          - virtualMassCoeff
           *(
                fvm::ddt(U1)
              + fvm::div(phi1, U1)
              - fvm::Sp(fvc::div(phi1), U1)
              - DDtU2
            )
        );
        mrfZones.addCoriolis(alpha1*rho1 + virtualMassCoeff, U1Eqn);
        U1Eqn.relax();
    }

    {
        U2Eqn =
        (
            fvm::ddt(alpha2, rho2, U2)
          + fvm::div(alphaRhoPhi2, U2)
          - fvm::Sp(fvc::ddt(alpha2, rho2) + fvc::div(alphaRhoPhi2), U2)
          + phase2.turbulence().divDevRhoReff(U2)
         ==
          - fvm::Sp(dragCoeff, U2)
          + (
                liftForce
              + wallLubricationForce
              + turbulentDispersionForce
            )
          - virtualMassCoeff
           *(
                fvm::ddt(U2)
              + fvm::div(phi2, U2)
              - fvm::Sp(fvc::div(phi2), U2)
              - DDtU1
            )
        );
        mrfZones.addCoriolis(alpha2*rho2 + virtualMassCoeff, U2Eqn);
        U2Eqn.relax();
    }
}

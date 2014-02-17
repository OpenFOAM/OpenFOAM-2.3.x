{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

    surfaceScalarField phig
    (
        (
            multiphaseProperties.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryField(),
        (
            phiHbyA.boundaryField()
          - (mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );

    PtrList<fvScalarMatrix> p_rghEqnComps(multiphaseProperties.phases().size());

    label phasei = 0;
    forAllConstIter
    (
        PtrDictionary<phaseModel>,
        multiphaseProperties.phases(),
        phase
    )
    {
        const rhoThermo& thermo = phase().thermo();
        const volScalarField& rho = thermo.rho()();

        p_rghEqnComps.set
        (
            phasei,
            (
                fvc::ddt(rho) + thermo.psi()*correction(fvm::ddt(p_rgh))
              + fvc::div(phi, rho) - fvc::Sp(fvc::div(phi), rho)
            ).ptr()
        );

        phasei++;
    }

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        tmp<fvScalarMatrix> p_rghEqnComp;

        phasei = 0;
        forAllConstIter
        (
            PtrDictionary<phaseModel>,
            multiphaseProperties.phases(),
            phase
        )
        {
            tmp<fvScalarMatrix> hmm
            (
                (max(phase(), scalar(0))/phase().thermo().rho())
               *p_rghEqnComps[phasei]
            );

            if (phasei == 0)
            {
                p_rghEqnComp = hmm;
            }
            else
            {
                p_rghEqnComp() += hmm;
            }

            phasei++;
        }

        solve
        (
            p_rghEqnComp
          + p_rghEqnIncomp,
            mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            // p = max(p_rgh + multiphaseProperties.rho()*gh, pMin);
            // p_rgh = p - multiphaseProperties.rho()*gh;

            phasei = 0;
            forAllIter
            (
                PtrDictionary<phaseModel>,
                multiphaseProperties.phases(),
                phase
            )
            {
                phase().dgdt() =
                    pos(phase())
                  *(p_rghEqnComps[phasei] & p_rgh)/phase().thermo().rho();
            }

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }

    p = max(p_rgh + multiphaseProperties.rho()*gh, pMin);

    // Update densities from change in p_rgh
    multiphaseProperties.correctRho(p_rgh - p_rgh_0);
    rho = multiphaseProperties.rho();

    K = 0.5*magSqr(U);

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;
}

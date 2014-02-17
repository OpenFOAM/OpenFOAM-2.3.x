    fvMesh& mesh = solidRegions[i];
    solidThermo& thermo = thermos[i];
    const radiation::radiationModel& radiation = radiations[i];

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volSymmTensorField> tAnialpha;

    if (!thermo.isotropic())
    {
        tmp<volVectorField> tkappaByCp = thermo.Kappa()/cp;

        const coordinateSystem& coodSys = coordinates[i];
        tAnialpha =
            tmp<volSymmTensorField>
            (
                new volSymmTensorField
                (
                    IOobject
                    (
                        "Anialpha",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    tkappaByCp().dimensions(),
                    zeroGradientFvPatchVectorField::typeName
                )
            );

        volSymmTensorField& Anialpha = tAnialpha();
        Anialpha.internalField() = coodSys.R().transformVector(tkappaByCp());
        Anialpha.correctBoundaryConditions();
    }

    volScalarField& h = thermo.he();

    const volScalarField& betav = betavSolid[i];

    fv::IOoptionList& fvOptions = solidHeatSources[i];

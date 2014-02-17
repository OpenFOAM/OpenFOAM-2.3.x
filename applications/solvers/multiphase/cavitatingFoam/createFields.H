    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhiv.H"
    #include "compressibleCreatePhi.H"

    Info<< "Reading transportProperties\n" << endl;

    incompressibleTwoPhaseMixture twoPhaseProperties(U, phiv);

    volScalarField& alphav(twoPhaseProperties.alpha1());
    alphav.oldTime();

    volScalarField& alphal(twoPhaseProperties.alpha2());

    Info<< "Creating compressibilityModel\n" << endl;
    autoPtr<barotropicCompressibilityModel> psiModel =
        barotropicCompressibilityModel::New
        (
            thermodynamicProperties,
            alphav
        );

    const volScalarField& psi = psiModel->psi();

    rho == max
    (
        psi*p
      + alphal*rhol0
      + ((alphav*psiv + alphal*psil) - psi)*pSat,
        rhoMin
    );

    // Create incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phiv, twoPhaseProperties)
    );

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
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

    #include "createPhi.H"


    // Transport
    // ~~~~~~~~~

    Info<< "Reading transportProperties\n" << endl;
    incompressibleTwoPhaseInteractingMixture twoPhaseProperties(U, phi);

    volScalarField& alpha1(twoPhaseProperties.alpha1());
    volScalarField& alpha2(twoPhaseProperties.alpha2());

    const dimensionedScalar& rho1 = twoPhaseProperties.rhod();
    const dimensionedScalar& rho2 = twoPhaseProperties.rhoc();

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Mixture density
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1*rho1 + alpha2*rho2
    );
    rho.oldTime();

    // Mass flux
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );


    // Relative Velocity
    // ~~~~~~~~~~~~~~~~~

    autoPtr<relativeVelocityModel> UdmModelPtr
    (
        relativeVelocityModel::New
        (
            transportProperties,
            twoPhaseProperties
        )
    );

    relativeVelocityModel& UdmModel(UdmModelPtr());


    // Turbulence
    // ~~~~~~~~~~

    IOdictionary RASProperties
    (
        IOobject
        (
            "RASProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


    Switch turbulence(RASProperties.lookup("turbulence"));

    dictionary kEpsilonDict(RASProperties.subDictPtr("kEpsilonCoeffs"));

    dimensionedScalar Cmu
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cmu",
            kEpsilonDict,
            0.09
        )
    );

    dimensionedScalar C1
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C1",
            kEpsilonDict,
            1.44
        )
    );

    dimensionedScalar C2
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C2",
            kEpsilonDict,
            1.92
        )
    );

    dimensionedScalar C3
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C3",
            kEpsilonDict,
            0.85
        )
    );

    dimensionedScalar sigmak
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmak",
            kEpsilonDict,
            1.0
        )
    );

    dimensionedScalar sigmaEps
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEps",
            kEpsilonDict,
            1.3
        )
    );

    dictionary wallFunctionDict(RASProperties.subDictPtr("wallFunctionCoeffs"));

    dimensionedScalar kappa
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "kappa",
            wallFunctionDict,
            0.41
        )
    );

    dimensionedScalar E
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "E",
            wallFunctionDict,
            9.8
        )
    );

    if (RASProperties.lookupOrDefault("printCoeffs", false))
    {
        Info<< "kEpsilonCoeffs" << kEpsilonDict << nl
            << "wallFunctionCoeffs" << wallFunctionDict << endl;
    }

    nearWallDist y(mesh);

    Info<< "Reading field k\n" << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field epsilon\n" << endl;
    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Calculating field mut\n" << endl;
    volScalarField mut
    (
        IOobject
        (
            "mut",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Cmu*rho*sqr(k)/epsilon
    );

    Info<< "Calculating field muEff\n" << endl;
    volScalarField muEff
    (
        IOobject
        (
            "muEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mut + twoPhaseProperties.mu()
    );


    // Pressure
    // ~~~~~~~~

    Info<< "Calculating field (g.h)f\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("gh", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }


    // MULES Correction
    // ~~~~~~~~~~~~~~~~

    tmp<surfaceScalarField> tphiAlphaCorr0;

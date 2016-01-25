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

//New. Read density of the first liquid and put it into the rho1-file.
//The file is created by swak4foam
    Info<< "Reading field rho1\n" << endl;
    volScalarField rho1
    (
        IOobject
        (
            "rho1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//New. Read density of the second liquid and put it into the rho2-file.
//The file is created by swak4foam
    Info<< "Reading field rho2\n" << endl;
    volScalarField rho2
    (
        IOobject
        (
            "rho2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//New. Read concentration of Spice (colloidal particles) from the S-file.
//The file is created by swak4foam
    Info<< "Reading field S\n" << endl;
    volScalarField S
    (
        IOobject
        (
            "S",
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

    Info<< "Reading transportProperties\n" << endl;
//Nothing was here after this phrase "Reading transportProperties"
//Standard notation is used for the case file "case_directory/constant/transportProperties"

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

//The density of carrier liquid
    dimensionedScalar rhoker
    (
         transportProperties.lookup("rhoker")
    );

//The density of the magnetic core (magnetite)
    dimensionedScalar rhomag
    (
         transportProperties.lookup("rhomag")
    );

//Einstein diffusion coefficient "kd" of a single colloidal particle in kerosene
    dimensionedScalar kd
    (
         transportProperties.lookup("kd")
    );

//Constant of coupling of two colloidal particles
    dimensionedScalar lambda
    (
         transportProperties.lookup("lambda")
    );

//Gravitation acceleration vector
//    dimensionedVector grav
//    (
//         transportProperties.lookup("grav")
//    );

    Info<< "AFTER reading transportProperties\n" << endl;


    volScalarField keff
    (
        IOobject
        (
            "keff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    kd*pow((1.0-S),6.5)*(1.0+2.0*S*(4.0-S)/pow((1.0-S),4)-8.0*pow(lambda,2.0)*S/(3.0*pow((1.0+1.25*lambda*S),2)))
    );


//Previously 'mixture()' was an object of the immiscibleIncompressibeTwoPhaseMixture-class
//and transfered 2 parameters to its Constructor: U, phi.
//    immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

//Now mixture() changed its class and transferes densities rho1, rho2 also.
    Info<< "Before creating mixture object\n" << endl;

    immiscibleIncompressibleTwoPhaseColloidMixture mixture(U, phi, rho1, rho2, S);

    Info<< "After creating mixture object\n" << endl;


    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

//This 2 stings of code are a blast from the past. They are inactive.
//    const dimensionedScalar& rho1 = mixture.rho1();
//    const dimensionedScalar& rho2 = mixture.rho2();


//Calculate the result density.
    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + alpha2*rho2,
        alpha1.boundaryField().types()
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


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, mixture)
    );

    #include "readGravitationalAcceleration.H"

    /*
    dimensionedVector g0(g);

    // Read the data file and initialise the interpolation table
    interpolationTable<vector> timeSeriesAcceleration
    (
        runTime.path()/runTime.caseConstant()/"acceleration.dat"
    );
    */

//const volVectorField& centres = Ua.mesh.C().internalField();

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh( "gh", g & mesh.C() );

    Info<< "After calculating field g.h\n" << endl;

    surfaceScalarField ghf("ghf", g & mesh.Cf());

    Info<< "Before calc p\n" << endl;

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

    Info<< "After calc p\n" << endl;

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


    fv::IOoptionList fvOptions(mesh);


    // MULES Correction
    tmp<surfaceScalarField> tphiAlphaCorr0;

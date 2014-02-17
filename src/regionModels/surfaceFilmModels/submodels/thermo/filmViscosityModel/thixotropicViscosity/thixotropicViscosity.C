/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "thixotropicViscosity.H"
#include "kinematicSingleLayer.H"
#include "addToRunTimeSelectionTable.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thixotropicViscosity, 0);

addToRunTimeSelectionTable
(
    filmViscosityModel,
    thixotropicViscosity,
    dictionary
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void thixotropicViscosity::updateMu()
{
    const kinematicSingleLayer& film = filmType<kinematicSingleLayer>();

    // blend based on mass fraction of added- to existing film mass
    const dimensionedScalar m0("zero", dimMass, 0.0);
    const dimensionedScalar mSMALL("SMALL", dimMass, ROOTVSMALL);
    const volScalarField deltaMass("deltaMass", max(m0, film.deltaMass()));
    const volScalarField filmMass("filmMass", film.netMass() + mSMALL);

    // weighting field to blend new and existing mass contributions
    const volScalarField w
    (
        "w",
        max(scalar(0.0), min(scalar(1.0), deltaMass/filmMass))
    );

    // evaluate thixotropic viscosity
    volScalarField muThx("muThx", muInf_/(sqr(1.0 - K_*lambda_) + ROOTVSMALL));

    // set new viscosity based on weight field
    mu_ = w*muInf_ + (1.0 - w)*muThx;
    mu_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thixotropicViscosity::thixotropicViscosity
(
    surfaceFilmModel& owner,
    const dictionary& dict,
    volScalarField& mu
)
:
    filmViscosityModel(typeName, owner, dict, mu),
    a_(coeffDict_.lookup("a")),
    b_(coeffDict_.lookup("b")),
    c_(coeffDict_.lookup("c")),
    d_(coeffDict_.lookup("d")),
    mu0_(coeffDict_.lookup("mu0")),
    muInf_(coeffDict_.lookup("muInf")),
    K_(1.0 - Foam::sqrt(muInf_/mu0_)),
    lambda_
    (
        IOobject
        (
            typeName + ":lambda",
            owner.regionMesh().time().timeName(),
            owner.regionMesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        ),
        owner.regionMesh()
    )
{
    lambda_.min(1.0);
    lambda_.max(0.0);

    // initialise viscosity to inf value
    // - cannot call updateMu() since this calls film.netMass() which
    //   cannot be evaluated yet (still in construction)
    mu_ = muInf_;
    mu_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thixotropicViscosity::~thixotropicViscosity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thixotropicViscosity::correct
(
    const volScalarField& p,
    const volScalarField& T
)
{
    const kinematicSingleLayer& film = filmType<kinematicSingleLayer>();

    // references to film fields
    const volVectorField& U = film.U();
    const volVectorField& Uw = film.Uw();
    const volScalarField& delta = film.delta();
    const volScalarField& deltaRho = film.deltaRho();
    const surfaceScalarField& phi = film.phi();

    // gamma-dot (shear rate) raised to the power d
    volScalarField gDotPowD
    (
        "gDotPowD",
        pow(mag(U - Uw)/(delta + film.deltaSmall()), d_)
    );

    dimensionedScalar c0("SMALL", dimMass/sqr(dimLength)/dimTime, SMALL);
    volScalarField coeff(-deltaRho*c_*gDotPowD + c0);

    fvScalarMatrix lambdaEqn
    (
        fvm::ddt(deltaRho, lambda_)
      + fvm::div(phi, lambda_)
      - fvm::Sp(fvc::div(phi), lambda_)
      ==
        deltaRho*a_*pow((1.0 - lambda_), b_)
      + fvm::SuSp(coeff, lambda_)
    );

    lambdaEqn.relax();

    lambdaEqn.solve();

    lambda_.min(1.0);
    lambda_.max(0.0);

    updateMu();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //

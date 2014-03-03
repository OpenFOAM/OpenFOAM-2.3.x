/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "plastic.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "incompressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(plastic, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        plastic,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::plastic::calcNu() const
{
    const incompressibleTwoPhaseMixture& twoPhaseProperties =
        alpha_.mesh().lookupObject<incompressibleTwoPhaseMixture>
        (
            "transportProperties"
        );

    bool isThisIsPhase1(&twoPhaseProperties.nuModel1() == this);

    dimensionedScalar
        rhoc
        (
            isThisIsPhase1
          ? twoPhaseProperties.rho2()
          : twoPhaseProperties.rho1()
        );

    dimensionedScalar
        rhop
        (
            isThisIsPhase1
          ? twoPhaseProperties.rho1()
          : twoPhaseProperties.rho2()
        );

    volScalarField
        nuc
        (
            (
                isThisIsPhase1
              ? twoPhaseProperties.nuModel2()
              : twoPhaseProperties.nuModel1()
            ).nu()
        );

    volScalarField
        nup
        (
            correctionNu(rhoc, rhop, nuc)
        );

    return
        max
        (
            nuMin_,
            min
            (
                nuMax_,
                (
                    nup + (rhoc/rhop)*nuc*alpha_
                )
            )
        )
       /max(alpha_, SMALL);
}


Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::plastic::correctionNu
(
    const dimensionedScalar& rhoc,
    const dimensionedScalar& rhop,
    const volScalarField& nuc
) const
{
    return
        plasticViscosityCoeff_
       *(
            pow
            (
                scalar(10),
                plasticViscosityExponent_*alpha_
            ) - scalar(1)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::plastic::plastic
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word modelName
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    plasticCoeffs_(viscosityProperties.subDict(modelName + "Coeffs")),
    plasticViscosityCoeff_
    (
        plasticCoeffs_.lookup("plasticViscosityCoeff")
    ),
    plasticViscosityExponent_
    (
        plasticCoeffs_.lookup("plasticViscosityExponent")
    ),
    nuMin_(plasticCoeffs_.lookup("nuMin")),
    nuMax_(plasticCoeffs_.lookup("nuMax")),
    alpha_
    (
        U.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                viscosityProperties.lookupOrDefault<word>("alpha", "alpha"),
                viscosityProperties.dictName()
            )
        )
    ),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimViscosity, 0)
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::plastic::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    plasticCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    plasticCoeffs_.lookup("k") >> plasticViscosityCoeff_;
    plasticCoeffs_.lookup("n") >> plasticViscosityExponent_;
    plasticCoeffs_.lookup("nuMin") >> nuMin_;
    plasticCoeffs_.lookup("nuMax") >> nuMax_;

    return true;
}


// ************************************************************************* //

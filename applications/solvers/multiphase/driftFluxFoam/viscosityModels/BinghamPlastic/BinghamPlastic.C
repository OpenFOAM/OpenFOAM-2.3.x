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

#include "BinghamPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(BinghamPlastic, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        BinghamPlastic,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::BinghamPlastic::correctionNu
(
    const dimensionedScalar& rhoc,
    const dimensionedScalar& rhop,
    const volScalarField& nuc
) const
{
    volScalarField
        tauy
        (
            yieldStressCoeff_
           *(
                pow
                (
                    scalar(10),
                    yieldStressExponent_
                   *(max(alpha_, scalar(0)) + yieldStressOffset_)
                )
              - pow
                (
                    scalar(10),
                    yieldStressExponent_*yieldStressOffset_
                )
            )
        );

    volScalarField
        nup
        (
            plastic::correctionNu(rhoc, rhop, nuc)
        );

    dimensionedScalar tauySmall("tauySmall", tauy.dimensions(), SMALL);

    return
        tauy
       /(
            mag(fvc::grad(U_))
          + 1.0e-4*(tauy + tauySmall)/(nup + (rhoc/rhop)*nuc)
        )
      + nup;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::BinghamPlastic::BinghamPlastic
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    plastic(name, viscosityProperties, U, phi, typeName),
    yieldStressCoeff_(plasticCoeffs_.lookup("yieldStressCoeff")),
    yieldStressExponent_(plasticCoeffs_.lookup("yieldStressExponent")),
    yieldStressOffset_(plasticCoeffs_.lookup("yieldStressOffset")),
    U_(U)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::BinghamPlastic::read
(
    const dictionary& viscosityProperties
)
{
    plastic::read(viscosityProperties);

    plasticCoeffs_.lookup("yieldStressCoeff") >> yieldStressCoeff_;
    plasticCoeffs_.lookup("yieldStressExponent") >> yieldStressExponent_;
    plasticCoeffs_.lookup("yieldStressOffset") >> yieldStressOffset_;

    return true;
}


// ************************************************************************* //

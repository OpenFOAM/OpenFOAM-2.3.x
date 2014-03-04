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

Class
    Foam::viscosityModels::plastic

Description
     Viscosity correction model for a generic power-law plastic.

SourceFiles
    plastic.C

\*---------------------------------------------------------------------------*/

#ifndef plastic_H
#define plastic_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleTwoPhaseMixture;

namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class plastic Declaration
\*---------------------------------------------------------------------------*/

class plastic
:
    public viscosityModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary plasticCoeffs_;

        //- Plastic viscosity coefficient
        dimensionedScalar plasticViscosityCoeff_;

        //- Plastic viscosity exponent
        dimensionedScalar plasticViscosityExponent_;

        //- Minimum viscosity
        dimensionedScalar nuMin_;

        //- Maximum viscosity
        dimensionedScalar nuMax_;

        //- Plastic phase fraction
        const volScalarField& alpha_;

        //- Viscosity
        volScalarField nu_;


    // Protected Member Functions

        //- Calculate and return the laminar viscosity
        virtual tmp<volScalarField> calcNu() const;

        //- Calculate and return the laminar viscosity correction
        virtual tmp<volScalarField> correctionNu
        (
            const dimensionedScalar& rhoc,
            const dimensionedScalar& rhop,
            const volScalarField& nuc
        ) const;


public:

    //- Runtime type information
    TypeName("plastic");


    // Constructors

        //- Construct from components
        plastic
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word modelName=typeName
        );


    //- Destructor
    ~plastic()
    {}


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::viscosityModels::BinghamPlastic

Description
     Viscosity correction model for Bingham plastics.

SourceFiles
    BinghamPlastic.C

\*---------------------------------------------------------------------------*/

#ifndef BinghamPlastic_H
#define BinghamPlastic_H

#include "plastic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class BinghamPlastic Declaration
\*---------------------------------------------------------------------------*/

class BinghamPlastic
:
    public plastic
{
protected:

    // Protected data

        //- Yield stress coefficient
        dimensionedScalar yieldStressCoeff_;

        //- Yield stress exponent
        dimensionedScalar yieldStressExponent_;

        //- Yield stress offset
        dimensionedScalar yieldStressOffset_;

        //- Velocity
        const volVectorField& U_;


    // Protected Member Functions

        //- Calculate and return the laminar viscosity correction
        virtual tmp<volScalarField> correctionNu
        (
            const dimensionedScalar& rhoc,
            const dimensionedScalar& rhop,
            const volScalarField& nuc
        ) const;


public:

    //- Runtime type information
    TypeName("BinghamPlastic");


    // Constructors

        //- Construct from components
        BinghamPlastic
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~BinghamPlastic()
    {}


    // Member Functions

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

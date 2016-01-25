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
    Foam::immiscibleIncompressibleTwoPhaseColloidMixture

Description
    An immiscible incompressible two-phase mixture transport model with spice
    transport mechanism.

SourceFiles
    immiscibleIncompressibleTwoPhaseColloidMixture.C

\*---------------------------------------------------------------------------*/

#ifndef immiscibleIncompressibleTwoPhaseColloidMixture_H
#define immiscibleIncompressibleTwoPhaseColloidMixture_H

#include "incompressibleTwoPhaseColloidMixture.H"
#include "interfaceProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class immiscibleIncompressibleTwoPhaseColloidMixture Declaration
\*---------------------------------------------------------------------------*/

class immiscibleIncompressibleTwoPhaseColloidMixture
:
    public incompressibleTwoPhaseColloidMixture,
    public interfaceProperties
{

public:

    // Constructors

        //- Construct from components
        immiscibleIncompressibleTwoPhaseColloidMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& rho1,
            const volScalarField& rho2,
            const volScalarField& S
        );


    //- Destructor
    virtual ~immiscibleIncompressibleTwoPhaseColloidMixture()
    {}


    // Member Functions

        //- Correct the transport and interface properties
        virtual void correct()
        {
            incompressibleTwoPhaseColloidMixture::correct();
            interfaceProperties::correct();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

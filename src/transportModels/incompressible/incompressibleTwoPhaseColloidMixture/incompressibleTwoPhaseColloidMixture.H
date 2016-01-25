/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::incompressibleTwoPhaseColloidMixture

Description
    A two-phase incompressible transportModel with spice diffusion mechanism

SourceFiles
    incompressibleTwoPhaseColloidMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoPhaseColloidMixture_H
#define incompressibleTwoPhaseColloidMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class incompressibleTwoPhaseColloidMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhaseColloidMixture
:
    public IOdictionary,
    public transportModel,
    public twoPhaseMixture
{
protected:

    // Protected data

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;

//Now rho1, rho2 are functions, not constant values.
//        dimensionedScalar rho1_;
//        dimensionedScalar rho2_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;
        const volScalarField& rho1_;
        const volScalarField& rho2_;
        const volScalarField& s_;

        volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    TypeName("incompressibleTwoPhaseColloidMixture");


    // Constructors

        //- Construct from components
        incompressibleTwoPhaseColloidMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& rho1,
            const volScalarField& rho2,
            const volScalarField& S
        );


    //- Destructor
    virtual ~incompressibleTwoPhaseColloidMixture()
    {}


    // Member Functions

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase1 density
        const volScalarField& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const volScalarField& rho2() const
        {
            return rho2_;
        };

        //- Return colloid concentration
        const volScalarField& s() const
        {
            return s_;
        };

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcNu();
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

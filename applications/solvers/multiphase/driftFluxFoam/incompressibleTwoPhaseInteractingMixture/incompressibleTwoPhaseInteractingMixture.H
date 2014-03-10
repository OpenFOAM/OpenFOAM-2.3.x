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
    Foam::incompressibleTwoPhaseInteractingMixture

Description
    A two-phase incompressible transportModel for interacting phases
    requiring the direct evaluation of the mixture viscosity,
    e.g. activated sludge or slurry.

SourceFiles
    incompressibleTwoPhaseInteractingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoPhaseInteractingMixture_H
#define incompressibleTwoPhaseInteractingMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "mixtureViscosityModel.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class incompressibleTwoPhaseInteractingMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoPhaseInteractingMixture
:
    public IOdictionary,
    public transportModel,
    public twoPhaseMixture
{
protected:

    // Protected data

        autoPtr<mixtureViscosityModel> muModel_;
        autoPtr<viscosityModel> nucModel_;

        dimensionedScalar rhod_;
        dimensionedScalar rhoc_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        volScalarField mu_;


public:

    TypeName("incompressibleTwoPhaseInteractingMixture");


    // Constructors

        //- Construct from components
        incompressibleTwoPhaseInteractingMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~incompressibleTwoPhaseInteractingMixture()
    {}


    // Member Functions

        //- Return const-access to the mixture viscosityModel
        const mixtureViscosityModel& muModel() const
        {
            return muModel_();
        }

        //- Return const-access to the continuous-phase viscosityModel
        const viscosityModel& nucModel() const
        {
            return nucModel_();
        }

        //- Return const-access to the dispersed-phase density
        const dimensionedScalar& rhod() const
        {
            return rhod_;
        }

        //- Return const-access to continuous-phase density
        const dimensionedScalar& rhoc() const
        {
            return rhoc_;
        };

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic mixture viscosity
        tmp<volScalarField> mu() const
        {
            return mu_;
        }

        //- Return the mixture viscosity
        virtual tmp<volScalarField> nu() const
        {
            notImplemented("incompressibleTwoPhaseInteractingMixture::nu()");
            return volScalarField::null();
        }

        //- Return the mixture viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            notImplemented
            (
                "incompressibleTwoPhaseInteractingMixture::nu(const label)"
            );
            return scalarField::null();
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            mu_ = muModel_->mu(rhoc_*nucModel_->nu());
        }

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

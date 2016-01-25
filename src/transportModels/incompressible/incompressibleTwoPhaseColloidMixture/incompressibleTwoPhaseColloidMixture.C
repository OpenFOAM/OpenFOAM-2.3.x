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

\*---------------------------------------------------------------------------*/

#include "incompressibleTwoPhaseColloidMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseColloidMixture, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void Foam::incompressibleTwoPhaseColloidMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // Average kinematic viscosity calculated from dynamic viscosity
//this was before new viscosity calculation
//  nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);

    nu_ = mu()*(exp((2.5*s_)/(1-s_))+(4.67*pow(s_,2))/(1-0.605*4.67*pow(s_,2)))/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseColloidMixture::incompressibleTwoPhaseColloidMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& rho1,
    const volScalarField& rho2,
    const volScalarField& S
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),

//Now rho1, rho2 are not constants, - they are functions.
//    rho1_("rho", dimDensity, nuModel1_->viscosityProperties().lookup("rho")),
//    rho2_("rho", dimDensity, nuModel2_->viscosityProperties().lookup("rho")),

    U_(U),
    phi_(phi),
    rho1_(rho1),
    rho2_(rho2),
    s_(S),

    nu_
    (
        IOobject
        (
            "nu",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseColloidMixture::mu() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            (limitedAlpha1*rho1_*nuModel1_->nu()
          + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu())*(exp((2.5*s_)/(1-s_))+(4.67*pow(s_,2))/(1-0.605*4.67*pow(s_,2)))
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseColloidMixture::muf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "muf",
            (alpha1f*fvc::interpolate(rho1_)*fvc::interpolate(nuModel1_->nu())
          + (scalar(1) - alpha1f)*fvc::interpolate(rho2_)*fvc::interpolate(nuModel2_->nu()))*(exp((2.5*fvc::interpolate(s_))/(1-fvc::interpolate(s_)))
          + (4.67*pow(fvc::interpolate(s_),2))/(1-0.605*4.67*pow(fvc::interpolate(s_),2)))
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseColloidMixture::nuf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuf",
            (exp((2.5*fvc::interpolate(s_))/(1-fvc::interpolate(s_)))
              + (4.67*pow(fvc::interpolate(s_),2))/(1-0.605*4.67*pow(fvc::interpolate(s_),2)))*(
                alpha1f*fvc::interpolate(rho1_)*fvc::interpolate(nuModel1_->nu())
              + (scalar(1) - alpha1f)*fvc::interpolate(rho2_)*fvc::interpolate(nuModel2_->nu())
            )/(alpha1f*fvc::interpolate(rho1_) + (scalar(1) - alpha1f)*fvc::interpolate(rho2_))
        )
    );
}


bool Foam::incompressibleTwoPhaseColloidMixture::read()
{
    if (regIOobject::read())
    {
        if
        (
            nuModel1_().read
            (
                subDict(phase1Name_ == "1" ? "phase1": phase1Name_)
            )
         && nuModel2_().read
            (
                subDict(phase2Name_ == "2" ? "phase2": phase2Name_)
            )
        )
        {
//            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
//            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

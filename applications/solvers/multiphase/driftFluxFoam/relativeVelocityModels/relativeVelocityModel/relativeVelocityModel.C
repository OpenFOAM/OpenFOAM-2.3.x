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

#include "relativeVelocityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(relativeVelocityModel, 0);
    defineRunTimeSelectionTable(relativeVelocityModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::relativeVelocityModel
(
    const dictionary& dict,
    const incompressibleTwoPhaseMixture& mixture
)
:
    mixture_(mixture),

    continuousPhaseName_(dict.lookup("continuousPhase")),

    alphaC_
    (
        mixture.phase1Name() == continuousPhaseName_
      ? mixture.alpha1()
      : mixture.alpha2()
    ),

    alphaD_
    (
        mixture.phase1Name() == continuousPhaseName_
      ? mixture.alpha2()
      : mixture.alpha1()
    ),

    rhoC_
    (
        mixture.phase1Name() == continuousPhaseName_
      ? mixture.rho1()
      : mixture.rho2()
    ),

    rhoD_
    (
        mixture.phase1Name() == continuousPhaseName_
      ? mixture.rho2()
      : mixture.rho1()
    ),

    Udm_
    (
        IOobject
        (
            "Udm",
            alphaC_.time().timeName(),
            alphaC_.mesh()
        ),
        alphaC_.mesh(),
        dimensionedVector("Udm", dimVelocity, vector::zero),
        mixture.U().boundaryField().types()
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::relativeVelocityModel> Foam::relativeVelocityModel::New
(
    const dictionary& dict,
    const incompressibleTwoPhaseMixture& mixture
)
{
    word modelType(dict.lookup(typeName));

    Info<< "Selecting relative velocity model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "relativeVelocityModel::New"
            "("
                "const dictionary&"
            ")"
        )   << "Unknown time scale model type " << modelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid time scale model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return
        autoPtr<relativeVelocityModel>
        (
            cstrIter()
            (
                dict.subDict(modelType + "Coeffs"),
                mixture
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModel::~relativeVelocityModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> Foam::relativeVelocityModel::rho() const
{
    return alphaC_*rhoC_ + alphaD_*rhoD_;
}


tmp<volSymmTensorField> Foam::relativeVelocityModel::tauDm() const
{
    volScalarField betaC(alphaC_*rhoC_);
    volScalarField betaD(alphaD_*rhoD_);

    // Calculate the relative velocity of the continuous phase w.r.t the mean
    volVectorField Ucm(betaD*Udm_/betaC);

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            "tauDm",
            betaD*sqr(Udm_) + betaC*sqr(Ucm)
        )
    );
}


// ************************************************************************* //

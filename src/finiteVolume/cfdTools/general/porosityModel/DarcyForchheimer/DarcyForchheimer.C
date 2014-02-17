/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "DarcyForchheimer.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(DarcyForchheimer, 0);
        addToRunTimeSelectionTable(porosityModel, DarcyForchheimer, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::DarcyForchheimer::DarcyForchheimer
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    D_(cellZoneIDs_.size()),
    F_(cellZoneIDs_.size()),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho")),
    muName_(coeffs_.lookupOrDefault<word>("mu", "thermo:mu")),
    nuName_(coeffs_.lookupOrDefault<word>("nu", "nu"))
{

    dimensionedVector d(coeffs_.lookup("d"));
    dimensionedVector f(coeffs_.lookup("f"));

    adjustNegativeResistance(d);
    adjustNegativeResistance(f);

    if (coordSys_.R().uniform())
    {
        forAll (cellZoneIDs_, zoneI)
        {
            D_[zoneI].setSize(1, tensor::zero);
            F_[zoneI].setSize(1, tensor::zero);

            D_[zoneI][0].xx() = d.value().x();
            D_[zoneI][0].yy() = d.value().y();
            D_[zoneI][0].zz() = d.value().z();

            D_[zoneI][0] = coordSys_.R().transformTensor(D_[zoneI][0]);

            // leading 0.5 is from 1/2*rho
            F_[zoneI][0].xx() = 0.5*f.value().x();
            F_[zoneI][0].yy() = 0.5*f.value().y();
            F_[zoneI][0].zz() = 0.5*f.value().z();

            F_[zoneI][0] = coordSys_.R().transformTensor(F_[zoneI][0]);
        }

    }
    else
    {
        forAll(cellZoneIDs_, zoneI)
        {
            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            D_[zoneI].setSize(cells.size(), tensor::zero);
            F_[zoneI].setSize(cells.size(), tensor::zero);

            forAll(cells, i)
            {
                D_[zoneI][i].xx() = d.value().x();
                D_[zoneI][i].yy() = d.value().y();
                D_[zoneI][i].zz() = d.value().z();

                // leading 0.5 is from 1/2*rho
                F_[zoneI][i].xx() = 0.5*f.value().x();
                F_[zoneI][i].yy() = 0.5*f.value().y();
                F_[zoneI][i].zz() = 0.5*f.value().z();
            }

            D_[zoneI] = coordSys_.R().transformTensor(D_[zoneI], cells);
            F_[zoneI] = coordSys_.R().transformTensor(F_[zoneI], cells);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::DarcyForchheimer::~DarcyForchheimer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::DarcyForchheimer::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), vector::zero);
    const scalarField& V = mesh_.V();

    apply(Udiag, Usource, V, rho, mu, U);

    force = Udiag*U - Usource;
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);
        const volScalarField& mu =
            mesh_.lookupObject<volScalarField>(muName_);

        apply(Udiag, Usource, V, rho, mu, U);
    }
    else
    {
        const volScalarField& nu =
            mesh_.lookupObject<volScalarField>(nuName_);

        apply(Udiag, Usource, V, geometricOneField(), nu, U);
    }
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    apply(Udiag, Usource, V, rho, mu, U);
}


void Foam::porosityModels::DarcyForchheimer::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);
        const volScalarField& mu =
            mesh_.lookupObject<volScalarField>(muName_);

        apply(AU, rho, mu, U);
    }
    else
    {
        const volScalarField& nu =
            mesh_.lookupObject<volScalarField>(nuName_);

        apply(AU, geometricOneField(), nu, U);
    }
}


bool Foam::porosityModels::DarcyForchheimer::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //

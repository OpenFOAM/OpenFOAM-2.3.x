/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "explicitPorositySource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "porosityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(explicitPorositySource, 0);
    addToRunTimeSelectionTable
    (
        option,
        explicitPorositySource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::explicitPorositySource::initialise()
{
    if (selectionMode_ != smCellZone)
    {
        FatalErrorIn("void Foam::fv::explicitPorositySource::initialise()")
            << "The porosity region must be specified as a cellZone.  Current "
            << "selection mode is " << selectionModeTypeNames_[selectionMode_]
            << exit(FatalError);
    }

    porosityPtr_.reset
    (
        porosityModel::New
        (
            name_,
            mesh_,
            coeffs_,
            cellSetName_
        ).ptr()
    ),

    fieldNames_.setSize(1, UName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::explicitPorositySource::explicitPorositySource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    porosityPtr_(NULL),
    UName_(coeffs_.lookupOrDefault<word>("UName", "U")),
    rhoName_(coeffs_.lookupOrDefault<word>("rhoName", "rho")),
    muName_(coeffs_.lookupOrDefault<word>("muName", "thermo:mu"))
{
    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::explicitPorositySource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    fvMatrix<vector> porosityEqn(eqn.psi(), eqn.dimensions());

    if (eqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);
        const volScalarField& mu = mesh_.lookupObject<volScalarField>(muName_);

        porosityPtr_->addResistance(porosityEqn, rho, mu);
    }
    else
    {
        porosityPtr_->addResistance(porosityEqn);
    }

    eqn -= porosityEqn;
}


void Foam::fv::explicitPorositySource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::explicitPorositySource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.readIfPresent("UName", UName_);
        coeffs_.readIfPresent("rhoName", rhoName_);
        coeffs_.readIfPresent("muName", muName_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

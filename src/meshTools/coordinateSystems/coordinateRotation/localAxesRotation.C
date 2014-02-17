/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "localAxesRotation.H"
#include "axesRotation.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "tensorIOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(localAxesRotation, 0);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        localAxesRotation,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        localAxesRotation,
        objectRegistry
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::localAxesRotation::localAxesRotation
(
    const dictionary& dict,
    const objectRegistry& orb
)
:
    Rptr_(),
    origin_(point::zero),
    e3_(vector::zero)
{
    // If origin is specified in the coordinateSystem
    if (dict.parent().found("origin"))
    {
        dict.parent().lookup("origin") >> origin_;
    }

    // rotation axis
    dict.lookup("e3") >> e3_;

    const polyMesh& mesh = refCast<const polyMesh>(orb);

    Rptr_.reset
    (
        new tensorField(mesh.nCells())
    );
    init(dict, orb);
}


Foam::localAxesRotation::localAxesRotation
(
    const dictionary& dict
)
:
    Rptr_(),
    origin_(),
    e3_()
{
    FatalErrorIn("localAxesRotation(const dictionary&)")
        << " localAxesRotation can not be constructed from  dictionary "
        << " use the construtctor : "
           "("
           "    const dictionary&, const objectRegistry&"
           ")"
        << exit(FatalIOError);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::localAxesRotation::clear()
{
    if (!Rptr_.empty())
    {
        Rptr_.clear();
    }
}


Foam::vector Foam::localAxesRotation::transform(const vector& st) const
{
    notImplemented
    (
        "vector localAxesRotation::transform(const vector&) const"
    );
    return vector::zero;
}


Foam::vector Foam::localAxesRotation::invTransform(const vector& st) const
{
    notImplemented
    (
        "vector localAxesRotation::invTransform(const vector&) const"
    );
    return vector::zero;
}


Foam::tmp<Foam::vectorField> Foam::localAxesRotation::transform
(
    const vectorField& st
) const
{
    if (Rptr_->size() != st.size())
    {
        FatalErrorIn
        (
            "tmp<vectorField> localAxesRotation::transform(const vectorField&)"
        )
            << "vectorField st has different size to tensorField "
            << abort(FatalError);
    }

    return (Rptr_() & st);
}


Foam::tmp<Foam::vectorField> Foam::localAxesRotation::invTransform
(
    const vectorField& st
) const
{
    return (Rptr_().T() & st);
}


Foam::tmp<Foam::tensorField> Foam::localAxesRotation::transformTensor
(
    const tensorField& st
) const
{
    if (Rptr_->size() != st.size())
    {
        FatalErrorIn
        (
            "tmp<tensorField> localAxesRotation::transformTensor"
            "("
                "const tensorField&"
            ")"
        )
            << "tensorField st has different size to tensorField Tr"
            << abort(FatalError);
    }
    return (Rptr_() & st & Rptr_().T());
}


Foam::tensor Foam::localAxesRotation::transformTensor
(
    const tensor& st
) const
{
    notImplemented
    (
        "tensor localAxesRotation::transformTensor(const tensor&) const"
    );

    return tensor::zero;
}


Foam::tmp<Foam::tensorField> Foam::localAxesRotation::transformTensor
(
    const tensorField& st,
    const labelList& cellMap
) const
{
    if (cellMap.size() != st.size())
    {
        FatalErrorIn
        (
            "tmp<tensorField> localAxesRotation::transformTensor"
            "("
                "const tensorField&"
                "const labelList&"
            ")"
        )
            << "tensorField st has different size to tensorField Tr"
            << abort(FatalError);
    }

    const tensorField Rtr(Rptr_().T());
    tmp<tensorField> tt(new tensorField(cellMap.size()));
    tensorField& t = tt();
    forAll(cellMap, i)
    {
        const label cellI = cellMap[i];
        t[i] = Rptr_()[cellI] & st[i] & Rtr[cellI];
    }
    return tt;
}


Foam::tmp<Foam::symmTensorField> Foam::localAxesRotation::transformVector
(
    const vectorField& st
) const
{
    if (Rptr_->size() != st.size())
    {
        FatalErrorIn("localAxesRotation::transformVector(const vectorField&)")
            << "tensorField st has different size to tensorField Tr"
            << abort(FatalError);
    }

    tmp<symmTensorField> tfld(new symmTensorField(Rptr_->size()));
    symmTensorField& fld = tfld();

    forAll(fld, i)
    {
        fld[i] = transformPrincipal(Rptr_()[i], st[i]);
    }
    return tfld;
}


Foam::symmTensor Foam::localAxesRotation::transformVector
(
    const vector& st
) const
{
    notImplemented
    (
        "tensor localAxesRotation::transformVector(const vector&) const"
    );
    return symmTensor::zero;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::localAxesRotation::init
(
    const dictionary& dict,
    const objectRegistry& obr
)
{
    const polyMesh& mesh = refCast<const polyMesh>(obr);
    forAll(mesh.cellCentres(), cellI)
    {
        vector dir = mesh.cellCentres()[cellI] - origin_;
        dir /= mag(dir) + VSMALL;

        Rptr_()[cellI] = axesRotation(e3_, dir).R();
    }
}


void Foam::localAxesRotation::write(Ostream& os) const
{
     os.writeKeyword("e3") << e3() << token::END_STATEMENT << nl;
}


// ************************************************************************* //

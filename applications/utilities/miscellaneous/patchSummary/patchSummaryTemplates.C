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

#include "patchSummaryTemplates.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::addToFieldList
(
    PtrList<GeometricField<Type, fvPatchField, volMesh> >& fieldList,
    const IOobject& obj,
    const label fieldI,
    const fvMesh& mesh
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obj.headerClassName() == fieldType::typeName)
    {
        fieldList.set
        (
            fieldI,
            new fieldType(obj, mesh)
        );
        Info<< "    " << fieldType::typeName << tab << obj.name() << endl;
    }
}


template<class Type>
void Foam::outputFieldList
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& fieldList,
    const label patchI
)
{
    forAll(fieldList, fieldI)
    {
        if (fieldList.set(fieldI))
        {
            Info<< "    " << pTraits<Type>::typeName << tab << tab
                << fieldList[fieldI].name() << tab << tab
                << fieldList[fieldI].boundaryField()[patchI].type() << nl;
        }
    }
}


template<class Type>
void Foam::collectFieldList
(
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& fieldList,
    const label patchI,
    HashTable<word>& fieldToType
)
{
    forAll(fieldList, fieldI)
    {
        if (fieldList.set(fieldI))
        {
            fieldToType.insert
            (
                fieldList[fieldI].name(),
                fieldList[fieldI].boundaryField()[patchI].type()
            );
        }
    }
}


// ************************************************************************* //

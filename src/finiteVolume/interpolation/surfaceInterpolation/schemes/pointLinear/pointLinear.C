/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "pointLinear.H"
#include "fvMesh.H"
#include "volPointInterpolation.H"
#include "triangle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::pointLinear<Type>::
correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    GeometricField<Type, pointPatchField, pointMesh> pvf
    (
        volPointInterpolation::New(mesh).interpolate(vf)
    );

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr =
        linearInterpolate(vf);

    Field<Type>& sfCorr = tsfCorr().internalField();

    const pointField& points = mesh.points();
    const pointField& C = mesh.C().internalField();
    const faceList& faces = mesh.faces();
    const scalarField& w = mesh.weights().internalField();
    const labelList& owner = mesh.owner();
    const labelList& neighbour = mesh.neighbour();

    forAll(sfCorr, facei)
    {
        point pi =
            w[owner[facei]]*C[owner[facei]]
          + (1.0 - w[owner[facei]])*C[neighbour[facei]];

        scalar at = triangle<point, const point&>
        (
            pi,
            points[faces[facei][0]],
            points[faces[facei][faces[facei].size()-1]]
        ).mag();

        scalar sumAt = at;
        Type sumPsip = at*(1.0/3.0)*
        (
            sfCorr[facei]
          + pvf[faces[facei][0]]
          + pvf[faces[facei][faces[facei].size()-1]]
        );

        for (label pointi=1; pointi<faces[facei].size(); pointi++)
        {
            at = triangle<point, const point&>
            (
                pi,
                points[faces[facei][pointi]],
                points[faces[facei][pointi-1]]
            ).mag();

            sumAt += at;
            sumPsip += at*(1.0/3.0)*
            (
                sfCorr[facei]
              + pvf[faces[facei][pointi]]
              + pvf[faces[facei][pointi-1]]
            );

        }

        sfCorr[facei] = sumPsip/sumAt - sfCorr[facei];
    }

    tsfCorr().boundaryField() = pTraits<Type>::zero;

    return tsfCorr;
}


namespace Foam
{
    makeSurfaceInterpolationScheme(pointLinear);
}

// ************************************************************************* //

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

#include "solidBodyMotionFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "transformField.H"
#include "cellZoneMesh.H"
#include "boolList.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyMotionFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, solidBodyMotionFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFvMesh::solidBodyMotionFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    dynamicMeshCoeffs_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                io.time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    ),
    SBMFPtr_(solidBodyMotionFunction::New(dynamicMeshCoeffs_, io.time())),
    undisplacedPoints_
    (
        IOobject
        (
            "points",
            io.time().constant(),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    zoneID_(-1),
    pointIDs_()
{
    if (undisplacedPoints_.size() != nPoints())
    {
        FatalIOErrorIn
        (
            "solidBodyMotionFvMesh::solidBodyMotionFvMesh(const IOobject&)",
            dynamicMeshCoeffs_
        )   << "Read " << undisplacedPoints_.size()
            << " undisplaced points from " << undisplacedPoints_.objectPath()
            << " but the current mesh has " << nPoints()
            << exit(FatalIOError);
    }

    word cellZoneName =
        dynamicMeshCoeffs_.lookupOrDefault<word>("cellZone", "none");

    if (cellZoneName != "none")
    {
        Info<< "Applying solid body motion to cellZone " << cellZoneName
            << endl;

        zoneID_ = cellZones().findZoneID(cellZoneName);

        if (zoneID_ == -1)
        {
            FatalErrorIn
            (
                "solidBodyMotionFvMesh::solidBodyMotionFvMesh(const IOobject&)"
            )
                << "Unable to find cellZone " << cellZoneName
                << ".  Valid celLZones are:"
                << cellZones().names()
                << exit(FatalError);
        }

        const cellZone& cz = cellZones()[zoneID_];


        // collect point IDs of points in cell zone

        boolList movePts(nPoints(), false);

        forAll(cz, i)
        {
            label cellI = cz[i];
            const cell& c = cells()[cellI];
            forAll(c, j)
            {
                const face& f = faces()[c[j]];
                forAll(f, k)
                {
                    label pointI = f[k];
                    movePts[pointI] = true;
                }
            }
        }

        syncTools::syncPointList(*this, movePts, orEqOp<bool>(), false);

        DynamicList<label> ptIDs(nPoints());
        forAll(movePts, i)
        {
            if (movePts[i])
            {
                ptIDs.append(i);
            }
        }

        pointIDs_.transfer(ptIDs);
    }
    else
    {
        Info<< "Applying solid body motion to entire mesh" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFvMesh::~solidBodyMotionFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidBodyMotionFvMesh::update()
{
    static bool hasWarned = false;

    if (zoneID_ != -1)
    {
        pointField transformedPts(undisplacedPoints_);

        UIndirectList<point>(transformedPts, pointIDs_) =
            transform
            (
                SBMFPtr_().transformation(),
                pointField(transformedPts, pointIDs_)
            );

        fvMesh::movePoints(transformedPts);
    }
    else
    {
        fvMesh::movePoints
        (
            transform
            (
                SBMFPtr_().transformation(),
                undisplacedPoints_
            )
        );
    }


    if (foundObject<volVectorField>("U"))
    {
        const_cast<volVectorField&>(lookupObject<volVectorField>("U"))
            .correctBoundaryConditions();
    }
    else if (!hasWarned)
    {
        hasWarned = true;

        WarningIn("solidBodyMotionFvMesh::update()")
            << "Did not find volVectorField U."
            << " Not updating U boundary conditions." << endl;
    }

    return true;
}


// ************************************************************************* //

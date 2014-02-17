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

#include "trackedParticle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::trackedParticle::trackedParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const point& end,
    const label level,
    const label i,
    const label j
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    end_(end),
    level_(level),
    i_(i),
    j_(j)
{}


Foam::trackedParticle::trackedParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> end_;
            level_ = readLabel(is);
            i_ = readLabel(is);
            j_ = readLabel(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&end_),
                sizeof(end_) + sizeof(level_) + sizeof(i_) + sizeof(j_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "trackedParticle::trackedParticle"
        "(const Cloud<trackedParticle>&, Istream&, bool)"
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::trackedParticle::move
(
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar dtMax = tEnd;

    if (tEnd <= SMALL)
    {
        // Remove the particle
        td.keepParticle = false;
    }
    else
    {
        td.keepParticle = true;

        while (td.keepParticle && !td.switchProcessor && tEnd > SMALL)
        {
            // set the lagrangian time-step
            scalar dt = min(dtMax, tEnd);

            // mark visited cell with max level.
            td.maxLevel()[cell()] = max(td.maxLevel()[cell()], level_);

            dt *= trackToFace(end_, td);

            tEnd -= dt;
            stepFraction() = 1.0 - tEnd/trackTime;
        }
    }

    return td.keepParticle;
}


bool Foam::trackedParticle::hitPatch
(
    const polyPatch&,
    trackingData& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    return false;
}


void Foam::trackedParticle::hitWedgePatch
(
    const wedgePolyPatch&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitSymmetryPlanePatch
(
    const symmetryPlanePolyPatch&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitSymmetryPatch
(
    const symmetryPolyPatch&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitCyclicPatch
(
    const cyclicPolyPatch&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    // Move to different processor
    td.switchProcessor = true;
}


void Foam::trackedParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices&
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitPatch
(
    const polyPatch& wpp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const trackedParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.end_
            << token::SPACE << p.level_
            << token::SPACE << p.i_
            << token::SPACE << p.j_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.end_),
            sizeof(p.end_) + sizeof(p.level_) + sizeof(p.i_) + sizeof(p.j_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const trackedParticle&)");

    return os;
}


// ************************************************************************* //

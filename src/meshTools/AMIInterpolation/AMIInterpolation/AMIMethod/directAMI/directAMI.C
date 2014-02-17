/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "directAMI.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::directAMI<SourcePatch, TargetPatch>::appendToDirectSeeds
(
    boolList& mapFlag,
    labelList& srcTgtSeed,
    DynamicList<label>& srcSeeds,
    DynamicList<label>& nonOverlapFaces,
    label& srcFaceI,
    label& tgtFaceI
) const
{
    const labelList& srcNbr = this->srcPatch_.faceFaces()[srcFaceI];
    const labelList& tgtNbr = this->tgtPatch_.faceFaces()[tgtFaceI];

    const pointField& srcPoints = this->srcPatch_.points();
    const pointField& tgtPoints = this->tgtPatch_.points();

    const vectorField& srcCf = this->srcPatch_.faceCentres();

    forAll(srcNbr, i)
    {
        label srcI = srcNbr[i];

        if (mapFlag[srcI] && srcTgtSeed[srcI] == -1)
        {
            const face& srcF = this->srcPatch_[srcI];
            const vector srcN = srcF.normal(srcPoints);

            bool found = false;
            forAll(tgtNbr, j)
            {
                label tgtI = tgtNbr[j];
                const face& tgtF = this->tgtPatch_[tgtI];
                pointHit ray = tgtF.ray(srcCf[srcI], srcN, tgtPoints);

                if (ray.hit())
                {
                    // new match - append to lists
                    found = true;

                    srcTgtSeed[srcI] = tgtI;
                    srcSeeds.append(srcI);

                    break;
                }
            }

            if (!found)
            {
                // no match available for source face srcI
                mapFlag[srcI] = false;
                nonOverlapFaces.append(srcI);
            }
        }
    }

    if (srcSeeds.size())
    {
        srcFaceI = srcSeeds.remove();
        tgtFaceI = srcTgtSeed[srcFaceI];
    }
    else
    {
        srcFaceI = -1;
        tgtFaceI = -1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::directAMI<SourcePatch, TargetPatch>::directAMI
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const scalarField& srcMagSf,
    const scalarField& tgtMagSf,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool reverseTarget
)
:
    AMIMethod<SourcePatch, TargetPatch>
    (
        srcPatch,
        tgtPatch,
        srcMagSf,
        tgtMagSf,
        triMode,
        reverseTarget
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::directAMI<SourcePatch, TargetPatch>::~directAMI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::directAMI<SourcePatch, TargetPatch>::calculate
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    label srcFaceI,
    label tgtFaceI
)
{
    bool ok =
        this->initialise
        (
            srcAddress,
            srcWeights,
            tgtAddress,
            tgtWeights,
            srcFaceI,
            tgtFaceI
        );

    if (!ok)
    {
        return;
    }


    // temporary storage for addressing and weights
    List<DynamicList<label> > srcAddr(this->srcPatch_.size());
    List<DynamicList<label> > tgtAddr(this->tgtPatch_.size());


    // construct weights and addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // list of faces currently visited for srcFaceI to avoid multiple hits
    DynamicList<label> srcSeeds(10);

    // list to keep track of tgt faces used to seed src faces
    labelList srcTgtSeed(srcAddr.size(), -1);
    srcTgtSeed[srcFaceI] = tgtFaceI;

    // list to keep track of whether src face can be mapped
    boolList mapFlag(srcAddr.size(), true);

    DynamicList<label> nonOverlapFaces;
    do
    {
        srcAddr[srcFaceI].append(tgtFaceI);
        tgtAddr[tgtFaceI].append(srcFaceI);

        mapFlag[srcFaceI] = false;

        // Do advancing front starting from srcFaceI, tgtFaceI
        appendToDirectSeeds
        (
            mapFlag,
            srcTgtSeed,
            srcSeeds,
            nonOverlapFaces,
            srcFaceI,
            tgtFaceI
        );
    } while (srcFaceI >= 0);

    if (nonOverlapFaces.size() != 0)
    {
        Pout<< "    AMI: " << nonOverlapFaces.size()
            << " non-overlap faces identified"
            << endl;

        this->srcNonOverlap_.transfer(nonOverlapFaces);
    }

    // transfer data to persistent storage
    forAll(srcAddr, i)
    {
        scalar magSf = this->srcMagSf_[i];
        srcAddress[i].transfer(srcAddr[i]);
        srcWeights[i] = scalarList(1, magSf);
    }
    forAll(tgtAddr, i)
    {
        scalar magSf = this->tgtMagSf_[i];
        tgtAddress[i].transfer(tgtAddr[i]);
        tgtWeights[i] = scalarList(1, magSf);
    }
}


// ************************************************************************* //

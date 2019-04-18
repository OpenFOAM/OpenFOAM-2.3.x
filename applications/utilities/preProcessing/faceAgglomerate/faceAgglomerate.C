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

Application
    faceAgglomerate

Description

    Agglomerate boundary faces using the pairPatchAgglomerationMod algorithm.
    It writes a map from the fine to coarse grid.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "CompactListList.H"
#include "unitConversion.H"
#include "pairPatchAgglomeration.H"
#include "pairPatchAgglomerationMod.H"
#include "labelListIOList.H"
#include "syncTools.H"
#include "globalIndex.H"

using namespace Foam;

void checkAndAgglomerate //this function call the agglomeration class and choses the standard or improved version
(
    const polyPatch& pp,
    labelListIOList& finalAgglom,
    const IOdictionary& agglomDict,
    label& nCoarseFaces,
    const label& patchI,
    const bool improved
)
{
    if (!pp.coupled())
    {
        Info << "\nAgglomerating patch : " << pp.name() << endl;
        if (!improved)
        {
            Info << "Using default algorithm" << endl;
            pairPatchAgglomeration agglomObject
            (
                pp,
                agglomDict.subDict(pp.name())
            );
            agglomObject.agglomerate();
            finalAgglom[patchI] = agglomObject.restrictTopBottomAddressing();
        }
        else
        {
            Info << "Using improved algorithm" << endl;
            pairPatchAgglomerationMod agglomObject
            (
                pp,
//                agglomDict.subDict(pp.name()),
                agglomDict.subOrEmptyDict(pp.name()),
                agglomDict
            );
            agglomObject.agglomerate();
            finalAgglom[patchI] = agglomObject.restrictTopBottomAddressing();
        }

        if (finalAgglom[patchI].size())
        {
            nCoarseFaces += max(finalAgglom[patchI] + 1);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "addDictOption.H"
    Foam::argList::addBoolOption
    (
        "improved",
        "run the improved version of the utility"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName("viewFactorsDict");

    #include "setConstantMeshDictionaryIO.H"

    // Read control dictionary
    const IOdictionary agglomDict(dictIO);

    bool writeAgglom = readBool(agglomDict.lookup("writeFacesAgglomeration"));

    bool agglomerateAllPatches = agglomDict.lookupOrDefault<bool>("agglomerateAllPatches", false);

    const bool improved = args.optionFound("improved");

    const polyBoundaryMesh& boundary = mesh.boundaryMesh();

    labelListIOList finalAgglom
    (
        IOobject
        (
            "finalAgglom",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        boundary.size()
    );

    label nCoarseFaces = 0;

    if (!agglomerateAllPatches || !improved) //agglomerating user-specified list of patches
    {
        Info << "Agglomerating user-specified patches" << endl;
        forAllConstIter(dictionary, agglomDict, iter)
        {
            labelList patchIds = boundary.findIndices(iter().keyword());
            forAll(patchIds, i)
            {
                label patchI =  patchIds[i];
                const polyPatch& pp = boundary[patchI];
                checkAndAgglomerate(pp, finalAgglom, agglomDict, nCoarseFaces, patchI, improved);
            }
        }
    }
    else //agglomerating all patches without explicit user list
    {
        Info << "Agglomerating all non-coupled patches" << endl;
        forAll(boundary, patchI)
        {
            const polyPatch& pp = boundary[patchI];
            checkAndAgglomerate(pp, finalAgglom, agglomDict, nCoarseFaces, patchI, improved);
        }
    }

    // - All patches which are not agglomarated are identity for finalAgglom
    forAll(boundary, patchId)
    {
        if (finalAgglom[patchId].size() == 0)
        {
            finalAgglom[patchId] = identity(boundary[patchId].size());
        }
    }

    // Sync agglomeration across coupled patches
    labelList nbrAgglom(mesh.nFaces() - mesh.nInternalFaces(), -1);

    forAll(boundary, patchId)
    {
        const polyPatch& pp = boundary[patchId];
        if (pp.coupled())
        {
            finalAgglom[patchId] = identity(pp.size());
            forAll(pp, i)
            {
                nbrAgglom[pp.start() - mesh.nInternalFaces() + i] =
                    finalAgglom[patchId][i];
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh, nbrAgglom);
    forAll(boundary, patchId)
    {
        const polyPatch& pp = boundary[patchId];
        if (pp.coupled() && !refCast<const coupledPolyPatch>(pp).owner())
        {
            forAll(pp, i)
            {
                finalAgglom[patchId][i] =
                    nbrAgglom[pp.start() - mesh.nInternalFaces() + i];
            }
        }
    }

    finalAgglom.write();

    if (writeAgglom)
    {
        globalIndex index(nCoarseFaces);
        volScalarField facesAgglomeration
        (
            IOobject
            (
                "facesAgglomeration",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("facesAgglomeration", dimless, 0)
        );

        label coarsePatchIndex = 0;
        forAll(boundary, patchId)
        {
            const polyPatch& pp = boundary[patchId];
            if (pp.size() > 0)
            {
                fvPatchScalarField& bFacesAgglomeration =
                    facesAgglomeration.boundaryField()[patchId];

                forAll(bFacesAgglomeration, j)
                {
                    bFacesAgglomeration[j] =
                        index.toGlobal
                        (
                            Pstream::myProcNo(),
                            finalAgglom[patchId][j] + coarsePatchIndex
                        );
                }

                coarsePatchIndex += max(finalAgglom[patchId]) + 1;
            }
        }

        Info<< "\nWriting facesAgglomeration" << endl;
        facesAgglomeration.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //

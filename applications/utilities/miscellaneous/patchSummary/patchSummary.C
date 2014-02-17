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

Application
    patchSummary

Description
    Writes fields and boundary condition info for each patch at each requested
    time instance.

    Default action is to write a single entry for patches/patchGroups with the
    same boundary conditions. Use the -expand option to print every patch
    separately. In case of multiple groups matching it will print only the
    first one.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "IOobjectList.H"
#include "patchSummaryTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "addRegionOption.H"
    argList::addBoolOption
    (
        "expand",
        "Do not combine patches"
    );
#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool expand = args.optionFound("expand");


#   include "createNamedMesh.H"
    const polyBoundaryMesh& bm = mesh.boundaryMesh();


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Update the mesh if changed
        if (mesh.readUpdate() == polyMesh::TOPO_PATCH_CHANGE)
        {
            Info<< "Detected changed patches. Recreating patch group table."
                << endl;
        }


        const IOobjectList fieldObjs(mesh, runTime.timeName());
        const wordList objNames = fieldObjs.names();

        PtrList<volScalarField> vsf(objNames.size());
        PtrList<volVectorField> vvf(objNames.size());
        PtrList<volSphericalTensorField> vsptf(objNames.size());
        PtrList<volSymmTensorField> vsytf(objNames.size());
        PtrList<volTensorField> vtf(objNames.size());

        Info<< "Valid fields:" << endl;

        forAll(objNames, objI)
        {
            IOobject obj
            (
                objNames[objI],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );

            if (obj.headerOk())
            {
                addToFieldList<scalar>(vsf, obj, objI, mesh);
                addToFieldList<vector>(vvf, obj, objI, mesh);
                addToFieldList<sphericalTensor>(vsptf, obj, objI, mesh);
                addToFieldList<symmTensor>(vsytf, obj, objI, mesh);
                addToFieldList<tensor>(vtf, obj, objI, mesh);
            }
        }

        Info<< endl;


        if (expand)
        {
            // Print each patch separately

            forAll(bm, patchI)
            {
                Info<< bm[patchI].type() << "\t: " << bm[patchI].name() << nl;
                outputFieldList<scalar>(vsf, patchI);
                outputFieldList<vector>(vvf, patchI);
                outputFieldList<sphericalTensor>(vsptf, patchI);
                outputFieldList<symmTensor>(vsytf, patchI);
                outputFieldList<tensor>(vtf, patchI);
                Info<< endl;
            }
        }
        else
        {
            // Collect for each patch the bc type per field. Merge similar
            // patches.

            // Per 'group', the map from fieldname to patchfield type
            DynamicList<HashTable<word> > fieldToTypes(bm.size());
            // Per 'group' the patches
            DynamicList<DynamicList<label> > groupToPatches(bm.size());
            forAll(bm, patchI)
            {
                HashTable<word> fieldToType;
                collectFieldList<scalar>(vsf, patchI, fieldToType);
                collectFieldList<vector>(vvf, patchI, fieldToType);
                collectFieldList<sphericalTensor>(vsptf, patchI, fieldToType);
                collectFieldList<symmTensor>(vsytf, patchI, fieldToType);
                collectFieldList<tensor>(vtf, patchI, fieldToType);

                label groupI = findIndex(fieldToTypes, fieldToType);
                if (groupI == -1)
                {
                    DynamicList<label> group(1);
                    group.append(patchI);
                    groupToPatches.append(group);
                    fieldToTypes.append(fieldToType);
                }
                else
                {
                    groupToPatches[groupI].append(patchI);
                }
            }


            forAll(groupToPatches, groupI)
            {
                const DynamicList<label>& patchIDs = groupToPatches[groupI];

                if (patchIDs.size() > 1)
                {
                    // Check if part of a group
                    wordList groups;
                    labelHashSet nonGroupPatches;
                    bm.matchGroups(patchIDs, groups, nonGroupPatches);

                    const labelList sortedPatches(nonGroupPatches.sortedToc());
                    forAll(sortedPatches, i)
                    {
                        Info<< bm[sortedPatches[i]].type()
                            << "\t: " << bm[sortedPatches[i]].name() << nl;
                    }
                    if (groups.size())
                    {
                        forAll(groups, i)
                        {
                            Info<< "group\t: " << groups[i] << nl;
                        }
                    }
                    outputFieldList<scalar>(vsf, patchIDs[0]);
                    outputFieldList<vector>(vvf, patchIDs[0]);
                    outputFieldList<sphericalTensor>(vsptf, patchIDs[0]);
                    outputFieldList<symmTensor>(vsytf, patchIDs[0]);
                    outputFieldList<tensor>(vtf, patchIDs[0]);
                    Info<< endl;
                }
                else
                {
                    // No group.
                    forAll(patchIDs, i)
                    {
                        label patchI = patchIDs[i];
                        Info<< bm[patchI].type()
                            << "\t: " << bm[patchI].name() << nl;
                        outputFieldList<scalar>(vsf, patchI);
                        outputFieldList<vector>(vvf, patchI);
                        outputFieldList<sphericalTensor>(vsptf, patchI);
                        outputFieldList<symmTensor>(vsytf, patchI);
                        outputFieldList<tensor>(vtf, patchI);
                        Info<< endl;
                    }
                }
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

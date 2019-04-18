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

#include "pairPatchAgglomerationMod.H"
#include "meshTools.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pairPatchAgglomerationMod::compactLevels(const label nCreatedLevels)
{ //compactLevels resisizes list to the final number of corser level found
    if (debugMode_) Info << "\tcompacting levels" << endl;
    nFaces_.setSize(nCreatedLevels);
    restrictAddressing_.setSize(nCreatedLevels);
    patchLevels_.setSize(nCreatedLevels);
}


bool Foam::pairPatchAgglomerationMod::continueAgglomerating //OFv1806
(
    const label nLocal,
    const label nLocalOld
)
{
    // Keep agglomerating
    // - if global number of faces is still changing
    // - and if local number of faces still too large (on any processor)
    //       or if global number of faces still too large

    label nGlobal = returnReduce(nLocal, sumOp<label>());
    label nGlobalOld = returnReduce(nLocalOld, sumOp<label>());

    return
    (
        returnReduce(nLocal > nFacesInCoarsestLevel_, orOp<bool>()) //stop if coarsest level has less faces than user-specified
     || nGlobal > nGlobalFacesInCoarsestLevel_
    )
    && nGlobal != nGlobalOld //stop if previous agglomeration had same number of coarse faces
    && scalar(nGlobal) >= scalar(nFacesInFinestLevel_)/100.0*percentFacesInCoarsestLevel_; //stop based on ratio between faces in coarsest and finest levels
}


void Foam::pairPatchAgglomerationMod::setBasedEdgeWeights()
{//this function is run only once from the constructor, does first level agglomeration and feature identification
    const bPatch& coarsePatch = patchLevels_[0]; //take level 0 patch i.e. finest patch
    forAll(coarsePatch.edges(), i) //loop over patch edges
    {
        if (coarsePatch.isInternalEdge(i)) //check if this edge is not a patch boundary edge
        {
            scalar edgeLength =
                coarsePatch.edges()[i].mag(coarsePatch.localPoints()); //extract lenght of edge

            const labelList& eFaces = coarsePatch.edgeFaces()[i]; //extract faces that has this edge in commin

            if (eFaces.size() == 2) //2 faces for the same edge
            {
                scalar cosI = //compute cos of angle between normal vectors of the faces to be used for feature identification
                    coarsePatch.faceNormals()[eFaces[0]]
                  & coarsePatch.faceNormals()[eFaces[1]];

                const edge edgeCommon = edge(eFaces[0], eFaces[1]); //build a new edge as pair of faces ID

                if (facePairWeight_.found(edgeCommon)) //if edge already found, add lenght
                {
                    facePairWeight_[edgeCommon] += edgeLength;
                }
                else //if the edge was not found in the list, add it
                {
                    facePairWeight_.insert(edgeCommon, edgeLength);
                }

                // if the angle between normal vectors is higher than user-specified, this is a geometric feature and value of -1 is used to identify it when agglomerating faces and building next levels
//                if (mag(cosI) < Foam::cos(degToRad(featureAngle_)))
                if (cosI < Foam::cos(degToRad(featureAngle_))) // FIXME check if magnitude is necessary
                {
                    facePairWeight_[edgeCommon] = -1.0;
                }
            }
            else
            {
                forAll(eFaces, j)
                {
                    for (label k = j+1; k<eFaces.size(); k++)
                    {
                        facePairWeight_.insert
                        (
                            edge(eFaces[j], eFaces[k]),
                            -1.0
                        );
                    }
                }
            }
        }
    }
}


void Foam::pairPatchAgglomerationMod::setEdgeWeights
( //similar to setBasedEdgeWeights but this is called for every next agglomeration levels
    const label fineLevelIndex
)
{
    const bPatch& coarsePatch = patchLevels_[fineLevelIndex]; //extract coarse (this level) patch

    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex]; //extract addressing finer->coarse addressing
    const label nCoarseI =  max(fineToCoarse) + 1;
    labelListList coarseToFine(invertOneToMany(nCoarseI, fineToCoarse)); //invert addressing to get coarse->fine addrsssing

    HashSet<edge, Hash<edge> > fineFeaturedFaces(coarsePatch.nEdges()/10);//build hashset for feature identification of coarser levels

    // Map fine faces with featured edge into coarse faces
    forAllConstIter(EdgeMap<scalar>, facePairWeight_, iter) //loop over finer level edges
    {
        if (iter() == -1.0) //we found a feature edge inside finer level edges
        {
            const edge e = iter.key();//extract finer level edge
            const edge edgeFeatured //build coarser edge as pair of the coarser faces related to the two finer faces
            (
                fineToCoarse[e[0]],
                fineToCoarse[e[1]]
            );
            fineFeaturedFaces.insert(edgeFeatured); //add this feature edge to hashset
        }
    }

    // Clean old weitghs
    facePairWeight_.clear(); //after we identified coarser patch (this patch) feature edges we reset facePairWeight_
    facePairWeight_.resize(coarsePatch.nEdges());

    scalar featureAngleLevel = Foam::min(featureAngle_*Foam::pow(featureAngleLevelsFactor_, fineLevelIndex), 90.0); //compute feature edge angle of this level
    scalar cosFeatureAngleLevel = Foam::cos(degToRad(featureAngleLevel));

    if (debugMode_)
    {
        if (featureAngleOnEveryLevel_)
        {
            Info << "\t\t\tcoarse feature angle: " << featureAngleLevel << endl;
        }
    }

    forAll(coarsePatch.edges(), i) //loop over coarse patch edges
    {
        if (coarsePatch.isInternalEdge(i)) //if this is not a coarse patch boundary edge
        {
            scalar edgeLength = //extract lenght of edge
                coarsePatch.edges()[i].mag(coarsePatch.localPoints());

            const labelList& eFaces = coarsePatch.edgeFaces()[i]; //extract coarse faces related to this edge

            if (eFaces.size() == 2) //2 faces for this edge 
            {
                scalar cosI = //compute the cos of the anle between the normal vectors of the coarse faces of this edge
                    coarsePatch.faceNormals()[eFaces[0]]
                  & coarsePatch.faceNormals()[eFaces[1]];

                const edge edgeCommon = edge(eFaces[0], eFaces[1]); //build edge as pair of coarse faces 
                if (facePairWeight_.found(edgeCommon))//if edge was already stored, add lenght
                {
                    facePairWeight_[edgeCommon] += edgeLength;
                }
                else//if edge was not stores, store it
                {
                    facePairWeight_.insert(edgeCommon, edgeLength);
                }

                if (!featureAngleOnEveryLevel_) //if coarse patch angle-based feature identification is off
                {
                    // If the fine 'pair' faces was featured edge so it is
                    // the coarse 'pair'
                    if (fineFeaturedFaces.found(edgeCommon))
                    {
                        facePairWeight_[edgeCommon] = -1.0;
                    }
                }
                else
                {
                    if (fineFeaturedFaces.found(edgeCommon) ||
//                        mag(cosI) < cosFeatureAngleLevel)
                        cosI < cosFeatureAngleLevel) //FIXME check if magnitude is necessary
                        {
                            facePairWeight_[edgeCommon] = -1.0;
                        }
                }
            }
            else
            {
                // Set edge as barrier by setting weight to -1
                forAll(eFaces, j)
                {
                    for (label k = j+1; k<eFaces.size(); k++)
                    {
                        facePairWeight_.insert
                        (
                            edge(eFaces[j], eFaces[k]),
                            -1.0
                        );
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairPatchAgglomerationMod::pairPatchAgglomerationMod
(
    const polyPatch& patch,
    const dictionary& patchDict,
    const dictionary& viewFactorsDict,
    const bool additionalWeights
)
:
    mergeLevels_
    (
        viewFactorsDict.lookupOrDefault<label>("mergeLevels", 2)
    ),
    maxLevels_(50),
    nFacesInCoarsestLevel_
    (
        viewFactorsDict.lookupOrDefault<label>("nFacesInCoarsestLevel", 10) //myMod
    ),
    nFacesInFinestLevel_(patch.size()),
    nGlobalFacesInCoarsestLevel_(labelMax), //OFv1806, but anyway initialized at labelMax
    percentFacesInCoarsestLevel_
    (
        viewFactorsDict.lookupOrDefault<scalar>("percentFacesInCoarsestLevel", 0.0)
    ),
    featureAngle_
    (
        viewFactorsDict.lookupOrDefault<scalar>("featureAngle", 10)
    ),
    featureAngleLevelsFactor_
    (
        viewFactorsDict.lookupOrDefault<scalar>("featureAngleLevelsFactor", 1.1)
    ),
    featureAngleOnEveryLevel_
    (
        viewFactorsDict.lookupOrDefault<bool>("featureAngleOnEveryLevel", false)
    ),
    nFaces_(maxLevels_),
    restrictAddressing_(maxLevels_),
    restrictTopBottomAddressing_(identity(patch.size())),
    patchLevels_(maxLevels_),
    facePairWeight_(patch.size()),
    debugMode_
    (
        viewFactorsDict.lookupOrDefault<bool>("debugMode", false)
    )
{
    //overwrite default or viewFactorsDict global settings with patch-specific settings if found
    if (patchDict.found("mergeLevels"))
    {
        mergeLevels_ = readLabel(patchDict.lookup("mergeLevels"));
    }
    if (patchDict.found("nFacesInCoarsestLevel"))
    {
        nFacesInCoarsestLevel_ = readLabel(patchDict.lookup("nFacesInCoarsestLevel"));
    }
    if (patchDict.found("percentFacesInCoarsestLevel"))
    {
        percentFacesInCoarsestLevel_ = readScalar(patchDict.lookup("percentFacesInCoarsestLevel"));
    }
    if (patchDict.found("featureAngle"))
    {
        featureAngle_ = readScalar(patchDict.lookup("featureAngle"));
    }
    if (patchDict.found("featureAngleLevelsFactor"))
    {
        featureAngleLevelsFactor_ = readScalar(patchDict.lookup("featureAngleLevelsFactor"));
    }
    if (patchDict.found("featureAngleOnEveryLevel"))
    {
        featureAngleOnEveryLevel_ = readBool(patchDict.lookup("featureAngleOnEveryLevel"));
    }

    // Set base fine patch
    patchLevels_.set
    (
        0,
        new bPatch
        (
            patch.localFaces(),
            patch.localPoints()
        )
    );

    // Set number of faces for the base patch
    nFaces_[0] = patch.size();

    // Set edge weights for level 0
    setBasedEdgeWeights();

    if (debugMode_)
    {
        Info << "Using mergeLevels: " << mergeLevels_ 
             << ", nFacesInCoarsestLevel: " << nFacesInCoarsestLevel_
             << ", percentFacesInCoarsestLevel: " << percentFacesInCoarsestLevel_
             << ", featureAngle: " << featureAngle_
             << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pairPatchAgglomerationMod::~pairPatchAgglomerationMod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::bPatch& Foam::pairPatchAgglomerationMod::patchLevel
(
    const label i
) const
{//gives coarse patch of level i
    return patchLevels_[i];
}


void Foam::pairPatchAgglomerationMod::mapBaseToTopAgglom
(
    const label fineLevelIndex
)
{//this is called every time a new level is found in order to update the finest to coarsest addressing
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];
    forAll (restrictTopBottomAddressing_, i)
    {
        restrictTopBottomAddressing_[i] =
            fineToCoarse[restrictTopBottomAddressing_[i]];
    }
}


bool Foam::pairPatchAgglomerationMod::agglomeratePatch
(
    const bPatch& patch,
    const labelList& fineToCoarse,
    const label fineLevelIndex
)
{
    if (min(fineToCoarse) == -1)
    {
        FatalErrorIn
        (
            "pairPatchAgglomerationMod::agglomeratePatch"
            "("
                "const bPatch&, "
                "const labelList&, "
                "const label"
            ")"
        )   << "min(fineToCoarse) == -1" << exit(FatalError);
    }

    if (fineToCoarse.size() == 0)
    {
        return false; //something went wrong during agglomerateOneLevel if addressing list is empty
    }

    if (fineToCoarse.size() != patch.size())
    {
        FatalErrorIn
        (
            "pairPatchAgglomerationMod::agglomeratePatch"
            "("
                "const bPatch&, "
                "const labelList&, "
                "const label"
            ")"
        )   << "restrict map does not correspond to fine level. " << endl
            << " Sizes: restrictMap: " << fineToCoarse.size()
            << " nEqns: " << patch.size()
            << abort(FatalError);
    }

    const label nCoarseI =  max(fineToCoarse)+1;
    List<face> patchFaces(nCoarseI);


    // Patch faces per agglomeration
    labelListList coarseToFine(invertOneToMany(nCoarseI, fineToCoarse));

    for (label coarseI = 0; coarseI < nCoarseI; coarseI++) //loop over coarse faces
    {
        const labelList& fineFaces = coarseToFine[coarseI]; //extract list of fine faces that belongs to this coarse face

        // Construct single face
        indirectPrimitivePatch upp //build a coarse face from the fine faces
        (
            IndirectList<face>(patch, fineFaces),
            patch.points()
        );

        if (upp.edgeLoops().size() != 1) //agglomerated face has topological issues
        {
            if (debugMode_) Info << "\t\t\tupp.edgeLoops().size() != 1" << endl;
            if (fineFaces.size() == 2)
            {
                const edge e(fineFaces[0], fineFaces[1]);
                facePairWeight_[e] = -1.0;
            }
            else if (fineFaces.size() == 3)
            {
                const edge e(fineFaces[0], fineFaces[1]);
                const edge e1(fineFaces[0], fineFaces[2]);
                const edge e2(fineFaces[2], fineFaces[1]);
                facePairWeight_[e] = -1.0;
                facePairWeight_[e1] = -1.0;
                facePairWeight_[e2] = -1.0;
            }

            return false; //return patch agglomeration failed
        }

        patchFaces[coarseI] = face //if the coarse face is ok in term of edgeLoops() then we can add the face with its points/loop
        (
            renumber
            (
                upp.meshPoints(),
                upp.edgeLoops()[0]
            )
        );
    }

    patchLevels_.set //set the new found patch (a list of coarse faces) inside the list of coarse patches
    (
        fineLevelIndex,
        new bPatch
        (
            SubList<face>(patchFaces, nCoarseI, 0),
            patch.points()
        )
    );

    return true;
}


void Foam::pairPatchAgglomerationMod:: agglomerate()
{//this is the main agglomeration algorithm that calls the other functions of this class
    label nPairLevels = 0;
    label nCreatedLevels = 1; //0 level is the base patch
    label nCoarseFaces = 0;
    label nCoarseFacesOld = 0;

    while (nCreatedLevels < maxLevels_) //loop until max number of levels is reached
    {
        if (debugMode_)
        {
            Info << "\tnCreatedLevels: " << nCreatedLevels << endl;
            Info << "\t\tnPairLevels: " << nPairLevels << endl;
        }        

        const bPatch& patch = patchLevels_[nCreatedLevels - 1]; //extract finer patch
        tmp<labelField> finalAgglomPtr(new labelField(patch.size())); //initialize agglomeration for coarse (this) level

        bool createdLevel = false;
        while (!createdLevel) //loop until agglomeration is successfull
        {
            finalAgglomPtr = agglomerateOneLevel //first thing to do is faces agglomeration, the function returns the fine->coarse addressing
            (
                nCoarseFaces, //overwritten with number of faces found for coarse level
                patch
            );

            if (debugMode_) Info << "\t\tnCoarseFaces: " << nCoarseFaces << endl;

            if (nCoarseFaces == 0) //if not faces found for coarse level, agglomeration was not successful, break
            {
                if (debugMode_) Info << "\t\tCoarseFaces == 0" << endl;
                break;
            }
            else //if agglomeration of faces was ok then create the new coarse patch from coarse faces
            {
                if (debugMode_) Info << "\t\tagglomeratePatch()" << endl;
                createdLevel = agglomeratePatch
                (
                    patch,
                    finalAgglomPtr,
                    nCreatedLevels
                );
            }
        }

        if (debugMode_)
        {
            Info << "\t\tcreatedLevel: " << createdLevel << endl;
        }

        if (createdLevel) //if agglomeration was successful continue the work (here the algorithm was different in the original version leading to possible crashes)
        {
            restrictAddressing_.set(nCreatedLevels, finalAgglomPtr); //store the fine->coarse addressing
            mapBaseToTopAgglom(nCreatedLevels); //update finest->coarsest addressing
            setEdgeWeights(nCreatedLevels);//update feature edge identification for coarse patch

            if (nPairLevels % mergeLevels_)//can do multiple aggloemrations inside same level
            {
                if (debugMode_) Info << "\t\tcombining levels" << endl;
                combineLevels(nCreatedLevels);//update the lists depending on mergeLevels value
            }
            else
            {
                if (debugMode_) Info << "\t\tincreasing level ID" << endl;
                nCreatedLevels++; //increase level ID
            }

            nPairLevels++; //increase pair ID

            nFaces_[nCreatedLevels] = nCoarseFaces; //store total number of coarse faces for this level
        }

        if
        (!continueAgglomerating(nCoarseFaces, nCoarseFacesOld))//check if we need to stop agglomeration based on multiple criteria
        {
            if (debugMode_)
            {
                Info << "\tStop agglomeration, nCoarseFaces: " << nCoarseFaces << ", nCoarseFacesOld: " << nCoarseFacesOld << endl;
            }
            break;
        }

        nCoarseFacesOld = nCoarseFaces; //update previous agglomeration number of faces value
    }
    compactLevels(nCreatedLevels); //resize lists to the final numer of found levels
}


Foam::tmp<Foam::labelField> Foam::pairPatchAgglomerationMod::agglomerateOneLevel
(
    label& nCoarseFaces,
    const bPatch& patch
)
{//this function perform the faces agglomeration based on feature edges found in previous steps
    const label nFineFaces = patch.size();

    tmp<labelField> tcoarseCellMap(new labelField(nFineFaces, -1)); //initialize addressing fine->coarse
    labelField& coarseCellMap = tcoarseCellMap();

    const labelListList& faceFaces = patch.faceFaces();

    nCoarseFaces = 0;

    forAll(faceFaces, facei)//loop over fine faces
    {
        const labelList& fFaces = faceFaces[facei]; //extract neighbour faces of this face

        if (coarseCellMap[facei] < 0) //if this face was not already agglomerated
        {
            label matchFaceNo = -1;
            label matchFaceNeibNo = -1;
            scalar maxFaceWeight = -GREAT;

            // check faces to find ungrouped neighbour with largest face weight
            forAll(fFaces, i) //loop over neighbour faces of this face
            {
                label faceNeig = fFaces[i];//extract neighbour face id
                const edge edgeCommon = edge(facei, faceNeig); //exract edge between this face and neighbour face
                if
                (
                    facePairWeight_[edgeCommon] > maxFaceWeight //used for optimization, try to maximize edge lenght when mergind faces
                 && coarseCellMap[faceNeig] < 0 //if neighbour face was not already merged to another face
                 && facePairWeight_[edgeCommon] != -1.0 //we can merge only if the edge is not a feature edge
                )
                {
                    // Match found. Pick up all the necessary data
                    matchFaceNo = facei;
                    matchFaceNeibNo = faceNeig;
                    maxFaceWeight = facePairWeight_[edgeCommon];
                }
            }

            if (matchFaceNo >= 0) //if merging was ok then store data inside addressing, the new coarse face is composed by the two faces sharing the same non-feature edge, this means they have both same "nCoarseFaces" coarse face ID
            {
                // Make a new group
                coarseCellMap[matchFaceNo] = nCoarseFaces;
                coarseCellMap[matchFaceNeibNo] = nCoarseFaces;
                nCoarseFaces++;
            }
            else //if was now possible to merge this face try another strategy, try to merge it to a group of already merged faces
            {
                // No match. Find the best neighbouring cluster and
                // put the cell there
                label clusterMatchFaceNo = -1;
                scalar clusterMaxFaceCoeff = -GREAT;

                forAll(fFaces, i) //loop over neighbour faces
                {
                    label faceNeig = fFaces[i];
                    const edge edgeCommon = edge(facei, faceNeig);
                    if
                    (
                        facePairWeight_[edgeCommon] > clusterMaxFaceCoeff
                        && facePairWeight_[edgeCommon] != -1.0 //cannot merge if this is a feature edge
                        && coarseCellMap[faceNeig] > 0 //can merge if neighour face has already been merged
                    )
                    {
                        clusterMatchFaceNo = faceNeig;
                        clusterMaxFaceCoeff = facePairWeight_[edgeCommon];
                    }
                }

                if (clusterMatchFaceNo >= 0) //FIXME in OFv1812 they use > instead of >=, why?
                {
                    // Add the cell to the best cluster
                    coarseCellMap[facei] = coarseCellMap[clusterMatchFaceNo];
                }
                else
                {
                    // if not create single-cell "clusters" for each
                    coarseCellMap[facei] = nCoarseFaces;
                    nCoarseFaces ++;
                }
            }
        }
    }

    // Check that all faces are part of clusters,

    for (label facei=0; facei<nFineFaces; facei++)
    {
        if (coarseCellMap[facei] < 0)
        {
            FatalErrorIn
            (
                "pairPatchAgglomerationMod::agglomerateOneLevel "
                "(label&, const bPatch&) "
            ) << " face " << facei
            << " is not part of a cluster"
            << exit(FatalError);
        }
    }

    return tcoarseCellMap;
}


void Foam::pairPatchAgglomerationMod::combineLevels(const label curLevel)
{//this function is used when mergeLevels > 1 to combine agglomerations into a single level
    label prevLevel = curLevel - 1;

    // Set the previous level nCells to the current
    nFaces_[prevLevel] = nFaces_[curLevel];

    // Map the restrictAddressing from the coarser level into the previous
    // finer level

    const labelList& curResAddr = restrictAddressing_[curLevel];
    labelList& prevResAddr = restrictAddressing_[prevLevel];

    forAll(prevResAddr, i)
    {
        prevResAddr[i] = curResAddr[prevResAddr[i]];
    }

    // Delete the restrictAddressing for the coarser level
    restrictAddressing_.set(curLevel, NULL);

    patchLevels_.set(prevLevel, patchLevels_.set(curLevel, NULL));
}


// ************************************************************************* //

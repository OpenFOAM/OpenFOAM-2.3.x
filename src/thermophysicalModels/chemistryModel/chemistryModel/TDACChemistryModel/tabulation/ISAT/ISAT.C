/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) <year> <author name(s)>
   \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "ISAT.H"
#include "error.H"
#include "TDACChemistryModel.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include "SLList.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
template<class CompType, class ThermoType>
Foam::ISAT<CompType, ThermoType>::ISAT
(
    const dictionary& chemistryProperties,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
:
    tabulation<CompType,ThermoType>(chemistryProperties, chemistry),
    chemisTree_(this->coeffsDict_),
    scaleFactor_(chemistry.nEqns(),1.0),
    runTime_(&chemistry.time()),
    chPMaxLifeTime_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "chPMaxLifeTime",
            (runTime_->endTime().value()-runTime_->startTime().value())
           /runTime_->deltaT().value()
        )
    ),
    maxGrowth_(this->coeffsDict_.lookupOrDefault("maxGrowth", INT_MAX)),
    MRURetrieve_(this->coeffsDict_.lookupOrDefault("MRURetrieve", false)),
    MRUSize_(this->coeffsDict_.lookupOrDefault("MRUSize", 0)),
    lastSearch_(NULL),
    growPoints_(this->coeffsDict_.lookupOrDefault("growPoints", true))
{
    if (this->active_)
    {
        dictionary scaleDict(this->coeffsDict_.subDict("scaleFactor"));
        label Ysize = this->chemistry_->Y().size();
        scalar otherScaleFactor = readScalar(scaleDict.lookup("otherSpecies"));
        for (label i = 0; i<Ysize; i++)
        {
            if (!scaleDict.found(this->chemistry_->Y()[i].name()))
            {
                scaleFactor_[i] = otherScaleFactor;
            }
            else
            {
                scaleFactor_[i] =
                    readScalar(scaleDict.lookup(this->chemistry_->Y()[i].name()));
            }
        }
        scaleFactor_[Ysize] = readScalar(scaleDict.lookup("Temperature"));
        scaleFactor_[Ysize+1] = readScalar(scaleDict.lookup("Pressure"));
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::ISAT<CompType, ThermoType>::~ISAT()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::addToMRU
(
    chemPointISAT<CompType, ThermoType>* phi0
)
{
    if (MRUSize_ > 0)
    {
        //first search if the chemPoint is already in the list
        bool isInList = false;
        typename SLList <chemPointISAT<CompType, ThermoType>*>::iterator iter =
            MRUList_.begin();
        for ( ; iter != MRUList_.end(); ++iter)
        {
            if (iter()==phi0)
            {
                isInList = true;
                break;
            }
        }
        //if it is in the list, then move it to front
        if (isInList)
        {
            if (iter()!=MRUList_.last())
            {
                //iter hold the position of the element to move
                MRUList_.remove(iter);
                //insert the element in front of the list
                MRUList_.append(phi0);
            }
        }
        else //chemPoint not yet in the list
        {
            if (MRUList_.size()==MRUSize_)
            {
                MRUList_.removeHead();
                MRUList_.append(phi0);
            }
            else
            {
                MRUList_.append(phi0);
            }
        }
    }
}


template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::calcNewC
(
    chemPointISAT<CompType, ThermoType>* phi0,
    const scalarField& phiq,
    scalarField& Rphiq
)
{
    label nEqns = this->chemistry_->nEqns();//full set of species
    bool mechRedActive = this->chemistry_->mechRed()->active();
    Rphiq = phi0->Rphi();
    scalarField dphi=phiq-phi0->phi();
    const List<List<scalar> >& Avar = phi0->A();
    List<label>& completeToSimplified(phi0->completeToSimplifiedIndex());

    //Rphiq[i]=Rphi0[i]+A[i][j]dphi[j]
    for (label i=0; i<nEqns-2; i++)
    {
        if (mechRedActive)
        {
            label si=completeToSimplified[i];
            //the species is active
            if (si!=-1)
            {
                for (label j=0; j<nEqns-2; j++)
                {
                    label sj=completeToSimplified[j];
                    if (sj!=-1)
                    {
                        Rphiq[i] += Avar[si][sj]*dphi[j];
                    }
                }
                Rphiq[i] += Avar[si][phi0->NsDAC()]*dphi[nEqns-2];
                Rphiq[i] += Avar[si][phi0->NsDAC()+1]*dphi[nEqns-1];
                //As we use an approximation of A, Rphiq should be ckecked for
                //negative values
                Rphiq[i] = max(0.0,Rphiq[i]);
            }
            //the species is not active A[i][j] = I[i][j]
            else
            {
                Rphiq[i] += dphi[i];
                Rphiq[i] = max(0.0,Rphiq[i]);
            }
        }
        else //mechanism reduction is not active
        {
            for (label j=0; j<nEqns; j++)
            {
                Rphiq[i] += Avar[i][j]*dphi[j];
            }
            //As we use an approximation of A, Rphiq should be ckecked for
            //negative values
            Rphiq[i] = max(0.0,Rphiq[i]);
        }
    }
}


template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::grow
(
    chemPointISAT<CompType, ThermoType>* phi0,
    const scalarField& phiq,
    const scalarField& Rphiq
)
{
    if (!phi0)
    {
        return false;
    }

    if ((phi0->nGrowth() < maxGrowth_) && !phi0->toRemove())
    {
        return (phi0->checkSolution(phiq,Rphiq));
    }
    else if (!phi0->toRemove())
    {
        cleaningRequired_ = true;
        phi0->toRemove() = true;
        bool inList(false);
        forAll(toRemoveList_,tRi)
        {
            if (toRemoveList_[tRi]==phi0)
            {
                inList=true;
                break;
            }
        }
        if (!inList)
        {
            toRemoveList_.append(phi0);
        }
        return false;
    }
    else
    {
        return false;
    }
}


template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::cleanAndBalance()
{
    bool treeModified(false);
    //1- check if the tree should be cleaned (flag from nUsed or nGrowth)
    if (cleaningRequired_)
    {
        cleaningRequired_=false;
        MRUList_.clear();
        //2- remove the points that have raised a flag because of number of
        // growth or used (they are stored in the toRemoveList)
        forAll(toRemoveList_,trli)
        {
            chemisTree_.deleteLeaf(toRemoveList_[trli]);
        }
        //set size to 0, the pointers have been deleted in deleteLeaf function
        toRemoveList_.clear();
        treeModified=true;
    }

    chemPointISAT<CompType, ThermoType>* x = chemisTree_.treeMin();
    while(x!=NULL)
    {
        chemPointISAT<CompType, ThermoType>* xtmp =
            chemisTree_.treeSuccessor(x);
        scalar elapsedTime = runTime_->timeOutputValue() - x->timeTag();
        scalar maxElapsedTime =
            chPMaxLifeTime_
          * runTime_->timeToUserTime(runTime_->deltaTValue());

        if (elapsedTime > maxElapsedTime)
        {
            chemisTree_.deleteLeaf(x);
            treeModified=true;
        }
        x = xtmp;
    }

    if (treeModified)
    {
        MRUList_.clear();
    }

    //return a bool to specify if the tree structure has been modified
    return treeModified;
}


template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::computeA
(
 scalarSquareMatrix& A,
 const scalarField& Rphiq,
 const scalar rhoi
 )
{
    scalar dt = runTime_->deltaT().value();
    bool mechRedActive = this->chemistry_->mechRed()->active();
    label speciesNumber=this->chemistry_->nSpecie();

    scalarField Rcq(this->chemistry_->nEqns());
    forAll(Rcq,i)
    {
        label s2c = this->chemistry_->simplifiedToCompleteIndex()[i];
        Rcq[i] = rho*Rphiq[s2c]/this->chemistry_->specieThermo()[s2c].W();
    }

    this->chemistry_->jacobian(runTime_->value(), Rcq, A);

    //the jacobian is computed according to the molar concentration
    //the following conversion allow to use A with mass fraction
    for (label i=0; i<speciesNumber; i++)
    {
        label si=i;
        if (mechRedActive)
        {
            si = this->chemistry_->simplifiedToCompleteIndex()[i];
        }
        for (label j=0; j<speciesNumber; j++)
        {
            label sj=j;
            if (mechRedActive)
            {
                sj = this->chemistry_->simplifiedToCompleteIndex()[j];
            }
            A[i][j] *=
                -dt*this->chemistry_->specieThermo()[si].W()
              / this->chemistry_->specieThermo()[sj].W();
        }
        A[i][i] += 1;
        //columns for pressure and temperature
        A[i][speciesNumber] *=
            -dt*this->chemistry_->specieThermo()[si].W()/rhoi;
        A[i][speciesNumber+1] *=
            -dt*this->chemistry_->specieThermo()[si].W()/rhoi;
    }
    //Before inversion
    A[speciesNumber][speciesNumber] += 1;
    A[speciesNumber+1][speciesNumber+1] += 1;
    gaussj(A, speciesNumber+2);

    //After inversion the last two lines of A are set to 0
    // only A[this->nSpecie()][this->nSpecie()]
    //    and A[this->nSpecie()+1][this->nSpecie()+1] !=0
    A[speciesNumber][speciesNumber] = 0.0;
    A[speciesNumber+1][speciesNumber+1] = 0.0;
} //end computeA function


template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::gaussj
(
 List<List<scalar> >& A,
 List<List<scalar> >& B,
 label n
 )
{
    //icol and irow initialized to 0 to make compiler happy (see below)
    label i, icol(0), irow(0), j, k, l, ll;
    scalar big, dum, pivinv;
    Field<label> indxc(n), indxr(n), ipiv(n);
    for (j=0; j<n; j++)
    {
        ipiv[j]=0;
    }
    for (i=0; i<n; i++)
    {
        big=0.0;
        for (j=0; j<n; j++)
        {
            if (ipiv[j] != 1)
            {
                for (k=0; k<n; k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs(A[j][k]) >= big)
                        {
                            big=fabs(A[j][k]);
                            irow=j;
                            icol=k;
                        }
                    }
                }
            }
        }
        ++(ipiv[icol]);
        if (irow != icol)
        {
            for (l=0; l<n; l++)
            {
                Swap(A[irow][l],A[icol][l]);
            }
            for (l=0; l<n; l++)
            {
                Swap(B[irow][l],B[icol][l]);
            }
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (A[icol][icol] == 0.0)
        {
            Info << "singular" << endl;
        }
        pivinv = 1.0/A[icol][icol];
        A[icol][icol] = 1.0;
        for (l=0; l<n; l++)
        {
            A[icol][l] *= pivinv;
        }
        for (l=0; l<n; l++)
        {
            B[icol][l] *= pivinv;
        }
        for (ll=0; ll<n; ll++)
        {
            if (ll != icol)
            {
                dum = A[ll][icol];
                A[ll][icol] = 0.0;
                for (l=0; l<n; l++)
                {
                    A[ll][l] -= A[icol][l]*dum;
                }
                for (l=0; l<n; l++)
                {
                    B[ll][l] -= B[icol][l]*dum;
                }
            }
        }
    }
    for (l=n-1; l>=0; l--)
    {
        if (indxr[l] != indxc[l])
        {
            for (k=0; k<n; k++)
            {
                Swap (A[k][indxr[l]],A[k][indxc[l]]);
            }
        }
    }
}//end gaussj(A,B)


template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::gaussj
(
 List<List<scalar> >& A,
 label n
 )
{
    List<List<scalar> > B(n,List<scalar>(n, 0.0));
    gaussj(A,B, n);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::retrieve
(
    const Foam::scalarField& phiq,
    scalarField& Rphiq
)
{
    if (chemisTree_.size())//if the tree is not empty
    {
        chemPointISAT<CompType, ThermoType>* phi0 =
            chemisTree_.binaryTreeSearch(phiq, chemisTree_.root());
        if (phi0->inEOA(phiq))
        {
            scalar elapsedTime = runTime_->timeOutputValue() - phi0->timeTag();
            scalar maxElapsedTime =
                chPMaxLifeTime_
              * runTime_->timeToUserTime(runTime_->deltaTValue());

            if (elapsedTime > maxElapsedTime && !phi0->toRemove())
            {
                cleaningRequired_ = true;
                phi0->toRemove() = true;
                bool inList(false);
                forAll(toRemoveList_,tRi)
                {
                    if (toRemoveList_[tRi]==phi0)
                    {
                        inList=true;
                        break;
                    }
                }
                if (!inList)
                {
                    toRemoveList_.append(phi0);
                }
            }
            addToMRU(phi0);
            totRetrieve_++;
            calcNewC(phi0,phiq,Rphiq);
            lastSearch_ = phi0;
            nFound_++;
            return true;
        }
        else if (chemisTree_.secondaryBTSearch(phiq, phi0))
        {
            scalar elapsedTime = runTime_->timeOutputValue() - phi0->timeTag();
            scalar maxElapsedTime =
            chPMaxLifeTime_
            * runTime_->timeToUserTime(runTime_->deltaTValue());

            if (elapsedTime > maxElapsedTime && !phi0->toRemove())
            {
                cleaningRequired_ = true;
                phi0->toRemove() = true;
                bool inList(false);
                forAll(toRemoveList_,tRi)
                {
                    if (toRemoveList_[tRi]==phi0)
                    {
                        inList=true;
                        break;
                    }
                }
                if (!inList)
                {
                    toRemoveList_.append(phi0);
                }
            }
            closest = phi0;
            phi0->lastTimeUsed()=runTime_->timeOutputValue();
            addToMRU(phi0);
            nFailedFirst_++;
            totRetrieve_++;
            calcNewC(phi0,phiq,Rphiq);
            lastSearch_ = phi0;
            nFound_++;
            return true;
        }
        else if (MRURetrieve_)
        {
            typename SLList
            <
                chemPointISAT<CompType, ThermoType>*
            >::iterator iter = MRUList_.begin();

            for ( ; iter != MRUList_.end(); ++iter)
            {
                phi0=iter();
                if (phi0->inEOA(phiq))
                {
                    scalar elapsedTime =
                        runTime_->timeOutputValue() - phi0->timeTag();
                    scalar maxElapsedTime =
                        chPMaxLifeTime_
                      * runTime_->timeToUserTime(runTime_->deltaTValue());

                    if (elapsedTime > maxElapsedTime && !phi0->toRemove())
                    {
                        cleaningRequired_ = true;
                        phi0->toRemove() = true;
                        bool inList(false);
                        forAll(toRemoveList_,tRi)
                        {
                            if (toRemoveList_[tRi]==phi0)
                            {
                                inList=true;
                                break;
                            }
                        }
                        if (!inList)
                        {
                            toRemoveList_.append(phi0);
                        }
                    }
                    phi0->lastTimeUsed()=runTime_->timeOutputValue();
                    addToMRU(phi0);
                    nFailedFirst_++;
                    totRetrieve_++;
                    calcNewC(phi0,phiq,Rphiq);
                    lastSearch_ = phi0;
                    nFound_++;
                    return true;
                }
            }
        }
    }
    //this point is reached when every retrieve trials have failed
    //or if the tree is empty
    lastSearch_ = NULL;
    return false;
}


template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::add
(
    const scalarField& phiq,
    const scalarField& Rphiq,
    const scalar rho
)
{
    //If lastSearch_ holds a valid pointer to a chem point AND the growPoints_
    //option is on, the code first tries to grow the point hold by lastSearch_
    if (lastSearch_ && growPoints_)
    {
        if (grow(lastSearch_,phiq,Rphiq))
        {
            nGrowth_++;
            return false;
        }
    }

    bool treeCleanedOrCleared(false);
    //If the code reach this point, it is either because lastSearch_ is not
    //valid, OR because growPoints_ is not on, OR because the grow operation
    //has failed. In the three cases, a new point is added to the tree.
    if (chemisTree().isFull())
    {
        //If cleanAndBalance operation do not result in a reduction of the tree
        //size, the last possibility is to delete completely the tree.
        //It can be partially rebuild with the MRU list if this is used.
        if (!cleanAndBalance())
        {
            DynamicList<chemPointISAT<CompType, ThermoType>*> tempList;
            if (MRUSize_>0)
            {
                //create a copy of each chemPointISAT of the MRUList_
                typename SLList
                <
                    chemPointISAT<CompType, ThermoType>*
                >::iterator iter = MRUList_.begin();
                for ( ; iter != MRUList_.end(); ++iter)
                {
                    tempList.append
                    (
                        new chemPointISAT<CompType, ThermoType>(*iter())
                    );
                }
            }
            chemisTree().clear();
            toRemoveList_.clear();
            MRUList_.clear();

            chemPointISAT<CompType, ThermoType>* nulPhi=0;
            forAll(tempList,i)
            {
                chemisTree().insertNewLeaf
                (
                     tempList[i]->phi(),
                     tempList[i]->Rphi(),
                     tempList[i]->A(),
                     scaleFactor(),
                     tolerance(),
                     scaleFactor_.size(),
                     nulPhi
                );
                deleteDemandDrivenData(tempList[i]);
            }
        }

        lastSearch_ = NULL;
        treeCleanedOrCleared = true;
    }

    //Compute the A matrix needed to store the chem point.
    label ASize = this->chemistry_->nEqns(); //reduced when mechRed is active
    scalarSquareMatrix A(ASize, ASize,0.0);
    computeA(A, Rphiq, rho);

    chemisTree().insertNewLeaf
    (
        phiq,
        Rphiq,
        A,
        scaleFactor(),
        tolerance(),
        scaleFactor_.size(),
        lastSearch_ //lastSearch_ may be NULL (handled by binaryTree)
    );

    nAdd_++;
    return treeCleanedOrCleared;

}


template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::writePerformance(fileName path)
{
    OFstream nFound(path + "found_isat.out");
    nFound << runTime_->timeOutputValue() << "    " <<  nFound_ <<endl;
    nFound_ = 0;

    OFstream nGrowth(path + "growth_isat.out");
    nGrowth << runTime_->timeOutputValue() << "    " <<  nGrowth_ <<endl;
    nGrowth_ = 0;

    OFstream nAdd(path + "add_isat.out");
    nAdd << runTime_->timeOutputValue() << "    " <<  nAdd_ <<endl;
    nAdd_ = 0;
}


// ************************************************************************* //

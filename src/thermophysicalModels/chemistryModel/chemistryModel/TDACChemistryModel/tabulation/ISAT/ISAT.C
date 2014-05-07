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
    chemisTree_(chemistry,this->coeffsDict_),
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
    checkEntireTreeInterval_
    (
        this->coeffsDict_.lookupOrDefault("checkEntireTreeInterval", INT_MAX)
    ),
    maxDepthFactor_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "maxDepthFactor",
            (chemisTree_.maxNLeafs()-1)
                /(std::log(chemisTree_.maxNLeafs())/std::log(2.0))
        )
    ),
    minBalanceThreshold_
    (
        this->coeffsDict_.lookupOrDefault
        (
            "minBalanceThreshold",0.1*chemisTree_.maxNLeafs()
        )
    ),
    MRURetrieve_(this->coeffsDict_.lookupOrDefault("MRURetrieve", false)),
    maxMRUSize_(this->coeffsDict_.lookupOrDefault("maxMRUSize", 0)),
    lastSearch_(NULL),
    growPoints_(this->coeffsDict_.lookupOrDefault("growPoints", true)),
    nRetrieved_(0),
    nGrowth_(0),
    nAdd_(0),
    nRetrievedFile_(chemistryProperties.name().path() + "/../found_isat.out"),
    nGrowthFile_(chemistryProperties.name().path() + "/../growth_isat.out"),
    nAddFile_(chemistryProperties.name().path() + "/../add_isat.out"),
    sizeFile_(chemistryProperties.name().path() + "/../size_isat.out"),
    cleaningRequired_(false)

{
    if (this->active_)
    {
        dictionary scaleDict(this->coeffsDict_.subDict("scaleFactor"));
        label Ysize = this->chemistry_.Y().size();
        scalar otherScaleFactor = readScalar(scaleDict.lookup("otherSpecies"));
        for (label i = 0; i<Ysize; i++)
        {
            if (!scaleDict.found(this->chemistry_.Y()[i].name()))
            {
                scaleFactor_[i] = otherScaleFactor;
            }
            else
            {
                scaleFactor_[i] =
                    readScalar
                    (
                        scaleDict.lookup(this->chemistry_.Y()[i].name())
                    );
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
    if (maxMRUSize_ > 0 && MRURetrieve_)
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
            if (iter()!=MRUList_.first())
            {
                //iter hold the position of the element to move
                MRUList_.remove(iter);
                //insert the element in front of the list
                MRUList_.insert(phi0);
            }
        }
        else //chemPoint not yet in the list, iter is last
        {
            if (MRUList_.size()==maxMRUSize_)
            {
                if (iter() == MRUList_.last())
                {
                    MRUList_.remove(iter);
                    MRUList_.insert(phi0);
                }
                else
                {
                    FatalErrorIn
                    (
                        "ISAT::addToMRU"
                    )   << "wrong MRUList construction"
                        << exit(FatalError);
                }
            }
            else
            {
                MRUList_.insert(phi0);
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
    label nEqns = this->chemistry_.nEqns();//full set of species
    bool mechRedActive = this->chemistry_.mechRed()->active();
    Rphiq = phi0->Rphi();
    scalarField dphi=phiq-phi0->phi();
    const scalarRectangularMatrix& gradientsMatrix = phi0->A();
    List<label>& completeToSimplified(phi0->completeToSimplifiedIndex());

    //Rphiq[i]=Rphi0[i]+A[i][j]dphi[j]
    //where Aij is dRi/dphi_j
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
                        Rphiq[i] += gradientsMatrix[si][sj]*dphi[j];
                    }
                }
                Rphiq[i] +=
                    gradientsMatrix[si][phi0->nActiveSpecies()]*dphi[nEqns-2];
                Rphiq[i] +=
                    gradientsMatrix[si][phi0->nActiveSpecies()+1]*dphi[nEqns-1];
                //As we use an approximation of A, Rphiq should be checked for
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
                Rphiq[i] += gradientsMatrix[i][j]*dphi[j];
            }
            //As we use a first order gradient matrix, Rphiq should be checked
            //for negative values
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
    //if the pointer to the chemPoint is NULL, the function stops
    if (!phi0)
    {
        return false;
    }

    //raise a flag when the chemPoint used has been grown more than the
    //allowed number of time
    if (!phi0->toRemove() && phi0->nGrowth() > maxGrowth_)
    {
        cleaningRequired_ = true;
        phi0->toRemove() = true;
    }

    //if the solution RphiQ is still within the tolerance we try to grow it
    //in some cases this might result in a failure and the grow function of
    //the chemPoint returns false
    if (phi0->checkSolution(phiq,Rphiq))
    {
        return phi0->grow(phiq);
    }
    //the actual solution and the approximation given by ISAT are too different
    else
    {
        return false;
    }
}


template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::cleanAndBalance()
{
    bool treeModified(false);

    //check all chemPoints to see if we need to delete some of the chemPoints
    //according to the ellapsed time and number of growths
    chemPointISAT<CompType, ThermoType>* x = chemisTree_.treeMin();
    while(x!=NULL)
    {
        chemPointISAT<CompType, ThermoType>* xtmp =
            chemisTree_.treeSuccessor(x);
        //timeOutputValue returns timeToUserTime(value()), therefore, it should
        //be compare with timeToUserTime(deltaT)
        scalar elapsedTime = runTime_->timeOutputValue() - x->timeTag();
        scalar maxElapsedTime =
            chPMaxLifeTime_
          * runTime_->timeToUserTime(runTime_->deltaTValue());

        if ((elapsedTime > maxElapsedTime) || (x->nGrowth() > maxGrowth_))
        {
            chemisTree_.deleteLeaf(x);
            treeModified=true;
        }
        x = xtmp;
    }

    //check if the tree should be balanced according to criterion:
    //  -the depth of the tree bigger than a*log2(size), log2(size) being the
    //      ideal depth (e.g. 4 leafs can be stored in a tree of depth 2)
    if
    (
        chemisTree_.size() > minBalanceThreshold_
     && chemisTree_.depth() >
            maxDepthFactor_*std::log(chemisTree_.size())/std::log(2.0)
    )
    {
        chemisTree_.balance();
        MRUList_.clear();
        treeModified = true;
    }

    //return a bool to specify if the tree structure has been modified
    return treeModified;
}


template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::computeA
(
 scalarRectangularMatrix& A,
 const scalarField& Rphiq,
 const scalar rhoi
 )
{
    scalar dt = runTime_->deltaTValue();
    bool mechRedActive = this->chemistry_.mechRed()->active();
    label speciesNumber=this->chemistry_.nSpecie();
    scalarField Rcq(this->chemistry_.nEqns());
    for (label i=0; i<speciesNumber; i++)
    {
        label s2c = i;
        if (mechRedActive)
        {
            s2c = this->chemistry_.simplifiedToCompleteIndex()[i];
        }
        Rcq[i] = rhoi*Rphiq[s2c]/this->chemistry_.specieThermo()[s2c].W();
    }
    Rcq[speciesNumber] = Rphiq[Rphiq.size()-2];
    Rcq[speciesNumber+1] = Rphiq[Rphiq.size()-1];

    // Aaa is computed implicitely,
    // A is given by A=C(psi0,t0+dt), where C is obtained through solving
    // d/dt C(psi0,t) = J(psi(t))C(psi0,t)
    // If we solve it implicitely:
    // (C(psi0, t0+dt) - C(psi0,t0))/dt = J(psi(t0+dt))C(psi0,t0+dt)
    // The Jacobian is thus computed according to the mapping
    // C(psi0,t0+dt)*(I-dt*J(psi(t0+dt))) = C(psi0, t0)
    // A = C(psi0,t0)/(I-dt*J(psi(t0+dt)))
    // where C(psi0,t0)=I

    this->chemistry_.jacobian(runTime_->value(), Rcq, A);

    //the jacobian is computed according to the molar concentration
    //the following conversion allows the code to use A with mass fraction
    for (label i=0; i<speciesNumber; i++)
    {
        label si=i;
        if (mechRedActive)
        {
            si = this->chemistry_.simplifiedToCompleteIndex()[i];
        }
        for (label j=0; j<speciesNumber; j++)
        {
            label sj=j;
            if (mechRedActive)
            {
                sj = this->chemistry_.simplifiedToCompleteIndex()[j];
            }
            A[i][j] *=
                -dt*this->chemistry_.specieThermo()[si].W()
              / this->chemistry_.specieThermo()[sj].W();
        }
        A[i][i] += 1;
        //columns for pressure and temperature
        A[i][speciesNumber] *=
            -dt*this->chemistry_.specieThermo()[si].W()/rhoi;
        A[i][speciesNumber+1] *=
            -dt*this->chemistry_.specieThermo()[si].W()/rhoi;
    }
    //For temperature and pressure, only unity on the diagonal
    A[speciesNumber][speciesNumber] = 1;
    A[speciesNumber+1][speciesNumber+1] = 1;
    //inverse of (I-dt*J(psi(t0+dt)))
    gaussj(A, speciesNumber+2);
} //end computeA function


template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::gaussj
(
 scalarRectangularMatrix& A,
 scalarRectangularMatrix& B,
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
 scalarRectangularMatrix& A,
 label n
 )
{
    scalarRectangularMatrix B(n,n, 0.0);
    gaussj(A, B, n);
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::retrieve
(
    const Foam::scalarField& phiq,
    scalarField& Rphiq
)
{
    bool retrieved(false);
    chemPointISAT<CompType, ThermoType>* phi0;
    if (chemisTree_.size())//if the tree is not empty
    {
        chemisTree_.binaryTreeSearch(phiq, chemisTree_.root(), phi0);

        //lastSearch keeps track of the chemPoint we obtain by the regular
        //binary tree search
        lastSearch_ = phi0;
        if (phi0->inEOA(phiq))
        {
            retrieved = true;
        }
        //after a successful secondarySearch, phi0 store a pointer to the
        //found chemPoint
        else if (chemisTree_.secondaryBTSearch(phiq, phi0))
        {
            retrieved = true;
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
                    retrieved = true;
                    break;
                }
            }
        }
    }
    //the tree is empty, retrieved is still false
    else
    {
        //there is no chempoints that we can try to grow
        lastSearch_ = NULL;
    }

    if (retrieved)
    {
        scalar elapsedTime =
            runTime_->timeOutputValue() - phi0->timeTag();
        scalar maxElapsedTime =
            chPMaxLifeTime_
          * runTime_->timeToUserTime(runTime_->deltaTValue());

        //raise a flag when the chemPoint has been used more than the allowed
        //number of time steps
        if (elapsedTime > maxElapsedTime && !phi0->toRemove())
        {
            cleaningRequired_ = true;
            phi0->toRemove() = true;
        }
        lastSearch_->lastTimeUsed()=runTime_->timeOutputValue();
        addToMRU(phi0);
        calcNewC(phi0,phiq,Rphiq);
        nRetrieved_++;
        return true;
    }
    else
    {
        //this point is reached when every retrieve trials have failed
        //or if the tree is empty
        return false;
    }
}


template<class CompType, class ThermoType>
bool Foam::ISAT<CompType, ThermoType>::add
(
    const scalarField& phiq,
    const scalarField& Rphiq,
    const scalar rho
)
{
    //If lastSearch_ holds a valid pointer to a chemPoint AND the growPoints_
    //option is on, the code first tries to grow the point hold by lastSearch_
    if (lastSearch_ && growPoints_)
    {
        if (grow(lastSearch_,phiq,Rphiq))
        {
            nGrowth_++;
            //the structure of the tree is not modified, return false
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
            if (maxMRUSize_>0)
            {
                //create a copy of each chemPointISAT of the MRUList_ before
                //they are deleted
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
            //pointers to chemPoint are not valid anymore, clear the list
            MRUList_.clear();

            //construct the tree without giving a reference to attach to it
            //since the structure has been completely discarded
            chemPointISAT<CompType, ThermoType>* nulPhi=0;
            forAll(tempList,i)
            {
                chemisTree().insertNewLeaf
                (
                     tempList[i]->phi(),
                     tempList[i]->Rphi(),
                     tempList[i]->A(),
                     scaleFactor(),
                     this->tolerance(),
                     scaleFactor_.size(),
                     nulPhi
                );
                deleteDemandDrivenData(tempList[i]);
            }
        }
        //the structure has been changed, it will force the binary tree to
        //perform a new search and find the most appropriate point still stored
        lastSearch_ = NULL;
        //either cleanAndBalance has changed the tree or it has been cleared
        //in any case treeCleanedOrCleared should be set to true
        treeCleanedOrCleared = true;
    }

    //Compute the A matrix needed to store the chemPoint.
    label ASize = this->chemistry_.nEqns(); //reduced when mechRed is active
    scalarRectangularMatrix A(ASize, ASize,0.0);
    computeA(A, Rphiq, rho);

    chemisTree().insertNewLeaf
    (
        phiq,
        Rphiq,
        A,
        scaleFactor(),
        this->tolerance(),
        scaleFactor_.size(),
        lastSearch_ //lastSearch_ may be NULL (handled by binaryTree)
    );

    nAdd_++;
    return treeCleanedOrCleared;

}


template<class CompType, class ThermoType>
void Foam::ISAT<CompType, ThermoType>::writePerformance()
{

    nRetrievedFile_
        << runTime_->timeOutputValue() << "    " <<  nRetrieved_ <<endl;
    nRetrieved_ = 0;

    nGrowthFile_ << runTime_->timeOutputValue() << "    " <<  nGrowth_ <<endl;
    nGrowth_ = 0;

    nAddFile_ << runTime_->timeOutputValue() << "    " <<  nAdd_ <<endl;
    nAdd_ = 0;

    sizeFile_ << runTime_->timeOutputValue() << "    " <<  this->size() <<endl;
}


// ************************************************************************* //

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

#include "binaryTree.H"
#include "binaryNode.H"
#include "demandDrivenData.H"
#include "clockTime.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
Foam::binaryTree<CompType, ThermoType>::binaryTree
(
    TDACChemistryModel<CompType, ThermoType>& chemistry,
    dictionary coeffsDict
)
:
    chemistry_(chemistry),
    root_(NULL),
    maxNLeafs_(readLabel(coeffsDict.lookup("maxNLeafs"))),
    size_(0),
    n2ndSearch_(0),
    max2ndSearch_(coeffsDict.lookupOrDefault("max2ndSearch",0)),
    minBalanceThreshold_
    (
        coeffsDict.lookupOrDefault("minBalanceThreshold",0.1*maxNLeafs_)
    ),
    maxNbBalanceTest_
    (
        coeffsDict.lookupOrDefault("maxNbBalanceTest",0.01*chemistry_.nSpecie())
    ),
    balanceProp_(coeffsDict.lookupOrDefault("balanceProp",0.35)),
    coeffsDict_(coeffsDict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
void Foam::binaryTree<CompType, ThermoType>::insertNewLeaf
(
 const scalarField& phiq,
 const scalarField& Rphiq, 
 const scalarSquareMatrix& A,
 const scalarField& scaleFactor, 
 const scalar& epsTol,
 const label nCols,
 chP*& phi0
 )
{
    if (size_ == 0) //no points are stored
    {
        //create an empty binary node and root points to it
        root_ = new bn();
        //create the new chemPoint which holds the composition point
        //phiq and the data to initialize the EOA
        chP* newChemPoint =
            new chP
            (
                chemistry_,
                phiq,
                Rphiq,
                A,
                scaleFactor,
                epsTol,
                nCols,
                coeffsDict_,
                root_
            );
        root_->leafLeft()=newChemPoint;
    }
    else //at least one point stored
    {
        //no reference chemPoint, a BT search is required
        if (phi0 == NULL) 
        {
            chemPointBase* phi0Base;
            binaryTreeSearch(phiq, root_,phi0Base);
            phi0 = dynamic_cast<chP*>(phi0Base);
        }
        //access to the parent node of the chemPoint
        bn* parentNode = phi0->node();
        
        //create the new chemPoint which holds the composition point
        //phiq and the data to initialize the EOA
        chP* newChemPoint =
            new chP
            (
                chemistry_,
                phiq,
                Rphiq,
                A,
                scaleFactor,
                epsTol,
                nCols,
                coeffsDict_
            );
        //insert new node on the parent node in the position of the
        //previously stored leaf (phi0)
        //the new node contains phi0 on the left and phiq on the right
        //the hyper plane is computed in the binaryNode constructor
        bn* newNode;
        if (size_>1)
        {
            newNode = new bn(phi0, newChemPoint, parentNode);
            //make the parent of phi0 point to the newly created node
            insertNode(phi0, newNode);
        }
        else //size_ == 1 (because not equal to 0)
        {
            //when size is 1, the binaryNode is empty without hyperplane
            deleteDemandDrivenData(root_);
            newNode = new bn(phi0, newChemPoint, NULL);
            root_ = newNode;
        }   
        
        phi0->node()=newNode;
        newChemPoint->node()=newNode;
    }
    size_++;
}


template<class CompType, class ThermoType>
void Foam::binaryTree<CompType, ThermoType>::insertNode
(
    chP*& phi0,
    bn*& newNode
)
{
    if (phi0==phi0->node()->leafRight())//phi0 is on the right
    {
        phi0->node()->leafRight() = NULL;
        phi0->node()->nodeRight() = newNode;
    }
    else
    {
        phi0->node()->leafLeft() = NULL;
        phi0->node()->nodeLeft() = newNode;
        
    }
}


template<class CompType, class ThermoType>
void Foam::binaryTree<CompType, ThermoType>::binaryTreeSearch
(
    const scalarField& phiq,
    bn* node,
    chemPointBase*& nearest
)
{
    if (size_ > 1)
    {
        scalar vPhi=0.0;
        const scalarField& v = node->v();
        const scalar& a = node->a();
        //compute v*phi
        for (label i=0; i<phiq.size(); i++) vPhi += phiq[i]*v[i];
        
        
        if (vPhi > a) //on right side (side of the newly added point)
        {
            if (node->nodeRight()!=NULL)
            {
                binaryTreeSearch(phiq, node->nodeRight(), nearest);
            }
            else //the terminal node is reached, return element on right
            {
                nearest = node->leafRight();
            }
        }
        else //on left side (side of the previously stored point)
        {
            if (node->nodeLeft()!=NULL)
            {
                binaryTreeSearch(phiq, node->nodeLeft(), nearest);
            }
            else //the terminal node is reached, return element on right
            {
                nearest = node->leafLeft();
            }
        }
    }
    // only one point stored (left element of the root)
    else if (size_ == 1)
    {
        nearest = root_->leafLeft();
    }
    else //no point stored
    {
        nearest = NULL;
    }
}//end binaryTreeSearch


template<class CompType, class ThermoType>
Foam::label Foam::binaryTree<CompType, ThermoType>::depth(bn* subTreeRoot)
{
    if (subTreeRoot == NULL)
    {
        return 0;
    }
    else 
    {
        return 1+max(depth(subTreeRoot->nodeLeft()),depth(subTreeRoot->nodeRight()));
    }
}


template<class CompType, class ThermoType>
void Foam::binaryTree<CompType, ThermoType>::deleteLeaf(chP*& phi0)
{

    if (size_ == 1) //only one point is stored
    {
        deleteDemandDrivenData(phi0);
        deleteDemandDrivenData(root_);
    }
    else if (size_ > 1)
    {
        bn* z = phi0->node();
        bn* x;
        chP* siblingPhi0 = chemPSibling(phi0);

        if (siblingPhi0 != NULL)//the sibling of phi0 is a chemPoint
        {
            //z was root (only two chemPoints in the tree)
            if (z->parent() == NULL)
            {
                root_ = new bn();
                root_->leafLeft()=siblingPhi0;
                siblingPhi0->node()=root_;
            }
            else if (z==z->parent()->nodeLeft())
            {
                z->parent()->leafLeft() = siblingPhi0;
                z->parent()->nodeLeft() = NULL;
                siblingPhi0->node() = z->parent();
            }
            else
            {
                z->parent()->leafRight() = siblingPhi0;
                z->parent()->nodeRight() = NULL;
                siblingPhi0->node() = z->parent();
            }
        }
        else
        {
            x = nodeSibling(phi0);
            transplant(z,x);
        }
        
        deleteDemandDrivenData(phi0);
        deleteDemandDrivenData(z);
    }
    size_--;
}//end of deleteLeaf


template<class CompType, class ThermoType>
void Foam::binaryTree<CompType, ThermoType>::transplant(bn* u, bn* v)
{
    if (u->parent() == NULL)
    {
        root_ = v;
    }
    else if (u == u->parent()->nodeLeft())
    {
        u->parent()->nodeLeft() = v;
    }
    else
    {
        u->parent()->nodeRight()=v;
    }

    if (v != NULL)
    {
        v->parent() = u->parent();
    }
}


template<class CompType, class ThermoType>
void Foam::binaryTree<CompType, ThermoType>::clear()
{
    //recursively delete the element in the subTree
    deleteSubTree();
    //reset root node (should already be NULL)
    root_=NULL;
    //reset size_
    size_=0;
}//end cleanAll


template<class CompType, class ThermoType>
void Foam::binaryTree<CompType, ThermoType>::deleteSubTree(bn* subTreeRoot)
{
    //no need to check for nullness it delete everything in the tree even NULL
    //pointer (no effect)
    if (subTreeRoot != NULL)
    {
        deleteDemandDrivenData(subTreeRoot->leafLeft());
        deleteDemandDrivenData(subTreeRoot->leafRight());
        deleteSubTree(subTreeRoot->nodeLeft());
        deleteSubTree(subTreeRoot->nodeRight());
        deleteDemandDrivenData(subTreeRoot);
    }
}


template<class CompType, class ThermoType>
void Foam::binaryTree<CompType, ThermoType>::deleteAllNode(bn* subTreeRoot)
{
    if (subTreeRoot != NULL)
    {
        deleteAllNode(subTreeRoot->nodeLeft());
        deleteAllNode(subTreeRoot->nodeRight());
        deleteDemandDrivenData(subTreeRoot);
    }
}


//Check if the tree has reached the maximum number of elements
template<class CompType, class ThermoType>
bool Foam::binaryTree<CompType, ThermoType>::isFull()
{
    return size_ >= maxNLeafs_;   
}


template<class CompType, class ThermoType>
bool Foam::binaryTree<CompType, ThermoType>::secondaryBTSearch
(
    const scalarField& phiq, 
    chP*& x
)
{
    //initialize n2ndSearch_
    n2ndSearch_ = 0;
    if ((n2ndSearch_ < max2ndSearch_) && (size_ > 1))
    {
        chP* xS = chemPSibling(x);
        if (xS != NULL)
        {
            n2ndSearch_++;
            if (xS->inEOA(phiq))
            {
                x=xS;
                return true;
            }
        }
        else if (inSubTree(phiq,nodeSibling(x),x))
        {
            return true;
        }
        bn* y = x->node();
        while((y->parent()!= NULL) && (n2ndSearch_ < max2ndSearch_))
        {
            xS = chemPSibling(y);
            if (xS != NULL)
            {
                n2ndSearch_++;
                if (xS->inEOA(phiq))
                {
                    x=xS;
                    return true;
                }
            }
            else if (inSubTree(phiq,nodeSibling(y),x))
            {
                return true;
            }
            y=y->parent();
        }
        //if we reach this point it is either because 
        //we did not find another covering EOA or 
        //we reach the maximm number of secondary search
        return false;
    }
    else
    {
        return false;
    }
}


template<class CompType, class ThermoType>
bool Foam::binaryTree<CompType, ThermoType>::inSubTree
(
    const scalarField& phiq, 
    bn* y,
    chP*& x
)
{
    if ((n2ndSearch_ < max2ndSearch_) && (y!=NULL))
    {
        scalar vPhi=0.0;
        const scalarField& v = y->v();
        const scalar& a = y->a();
        //compute v*phi
        for (label i=0; i<phiq.size(); i++)
        {
            vPhi += phiq[i]*v[i];
        }
        if (vPhi<=a)//on the left side of the node
        {
            if (y->nodeLeft() == NULL)//left is a chemPoint
            {
                n2ndSearch_++;
                x=y->leafLeft();
                if (x->inEOA(phiq))
                {
                    return true;
                }
            }
            else//the left side is a node
            {
                if (inSubTree(phiq,y->nodeLeft(),x))
                {
                    return true;
                }
            }    
            
            if ((n2ndSearch_ < max2ndSearch_) && y->nodeRight() == NULL)
            {
                n2ndSearch_++;
                x=y->leafRight();
                return x->inEOA(phiq);
            }
            else//test for n2ndSearch is done in the call of inSubTree
            {
                return inSubTree(phiq,y->nodeRight(),x);
            }
        }//end of "on left side"
        else //on right side (symetric of above)
        {
            if (y->nodeRight() == NULL)
            {
                n2ndSearch_++;
                x=y->leafRight();
                if (x->inEOA(phiq))
                {
                    return true;
                }
            }
            else//the right side is a node
            {
                if (inSubTree(phiq,y->nodeRight(),x))
                {
                    return true;
                }
            }
            //if we reach this point, the retrieve has 
            //failed on the right side, explore the left side
            if ((n2ndSearch_ < max2ndSearch_) && y->nodeLeft() == NULL)
            {
                n2ndSearch_++;
                x=y->leafLeft();
                return x->inEOA(phiq);
            }
            else 
            {
                return inSubTree(phiq,y->nodeLeft(),x);
            }
        }
    }//end if ((n2ndSearch_ < max2ndSearch_) && (y!=NULL))
    else 
    {
        return false;
    }
}


template<class CompType, class ThermoType>
Foam::chemPointISAT<CompType, ThermoType>*
Foam::binaryTree<CompType, ThermoType>::chemPSibling(bn* y)
{
    if (y->parent()!=NULL)
    {
        if (y==y->parent()->nodeLeft())//y is on the left, return right side
        {
            return y->parent()->leafRight();
        }
        else //y is on the right, return left side
        {
            return y->parent()->leafLeft();
        }
    }
    else 
    {
        return NULL;
    }
}


template<class CompType, class ThermoType>
Foam::chemPointISAT<CompType, ThermoType>*
Foam::binaryTree<CompType, ThermoType>::chemPSibling(chP* x)
{
    if (size_>1)
    {
        bn* y = x->node();
        if (x==y->leafLeft())//x is on the left, return right side
        {
            return y->leafRight();
        }
        else//x is on the right, return left side
        {
            return y->leafLeft();
        }
    }
    else
    {
        return NULL;
    }
}


template<class CompType, class ThermoType>
Foam::binaryNode<CompType, ThermoType>*
Foam::binaryTree<CompType, ThermoType>::nodeSibling(bn* y)
{
    if (y->parent()!=NULL)
    {
        if (y==y->parent()->nodeLeft())//y is on the left, return right side
        {
            return y->parent()->nodeRight();
        }
        else //y is on the right, return left side
        {
            return y->parent()->nodeLeft();
        }
    }
    else 
    {
        return NULL;
    }
}


template<class CompType, class ThermoType>
Foam::binaryNode<CompType, ThermoType>*
Foam::binaryTree<CompType, ThermoType>::nodeSibling(chP* x)
{
    if (size_>1)
    {
        bn* y = x->node();
        if (x==y->leafLeft())//x is on the left, return right side
        {
            return y->nodeRight();
        }
        else//x is on the right, return left side
        {
            return y->nodeLeft();
        }
    }
    else 
    {
        return NULL;
    }
}


template<class CompType, class ThermoType>
bool Foam::binaryTree<CompType, ThermoType>::balance()
{
    if (size_ > minBalanceThreshold_)
    {
        scalarField mean(chemistry_.nEqns(),0.0);//use the size of the space
        SortableList<scalar> variance(chemistry_.nEqns(),0.0);
        
        //1) walk through the entire tree by starting with the tree's most left
        //chemPoint
        chP* x=treeMin();
        List<chP*> chemPoints(size_);
        label chPi=0;
        //2) compute the mean composition
        while(x!=NULL)
        {
            const scalarField& phij = x->phi();
            mean += phij;
            chemPoints[chPi++]=x;
            x=treeSuccessor(x);
        }    
        mean /= size_;

        //3) compute the variance for each space direction
        forAll(chemPoints,j)
        {
            const scalarField& phij = chemPoints[j]->phi();
            forAll(variance,vi)
            {
                variance[vi] += sqr(phij[vi]-mean[vi]);
            }
        }
        
        //4) analyze which direction of the cutting plane better separate the
        //  points (maximal variance)
        variance.sort();
        label maxDir(-1);
        //maxDir indicates the direction of maximum variance
        //the vector v[maxDir]=1; v[i!=maxDir]=0 indicates the perpendicular
        //direction and can be used to define the hyperplane (like in DOLFA)
        //instead we create the new root node by taking the two extreme points
        //in this direction if these extreme points were not deleted in the
        //cleaning that come before the equilibrate they are still important and
        //the tree should therefore take them into account
        label nbLeft(0);
        label nbTests(0);
        scalar bestBalance(size_);
        while
        (
            ((nbLeft<balanceProp_*size_) || (nbLeft>((1-balanceProp_)*size_)))
         && ((nbTests<maxNbBalanceTest_) && (nbTests<variance.size()-1))
        )
        {
            nbLeft=0;
            nbTests++;
            //variance.indices refers to the directions in the composition space
            //variance.indices[variance.size()-(nbTests)] refers to the
            //direction of the nbTests highest value of variance
            label curDir = variance.indices()[variance.size()-(nbTests)];
            forAll(chemPoints,j)
            {
                scalar phiMaxDir = chemPoints[j]->phi()[curDir];
                if (phiMaxDir < mean[curDir])
                {
                    nbLeft++;
                }
            }
            //with bestBalance starting at size
            //the first while loop will go in this if loop
            if (fabs(nbLeft-size_*0.5)<bestBalance)
            {
                bestBalance = fabs(nbLeft-size_*0.5);
                maxDir = curDir;
            }
        }

        scalar maxPhi(0);
        scalar minPhi(GREAT);
        label minId(-1);
        label maxId(-1);
        forAll(chemPoints,j)
        {
            scalar phiMaxDir = chemPoints[j]->phi()[maxDir];
            if (phiMaxDir>maxPhi)
            {
                maxId = j;
                maxPhi=phiMaxDir;
            }
            if (phiMaxDir<minPhi)
            {
                minId = j;
                minPhi = phiMaxDir;
            }
        }
        chP* minRef = chemPoints[minId];
        chP* maxRef = chemPoints[maxId];
        
        //delete reference to all node since the tree is reshaped
        deleteAllNode();
        root_=NULL;
        
        //add the node for minRef and maxRef
        
        bn* newNode = new bn(minRef, maxRef, NULL);
        root_ = newNode;
        minRef->node() = newNode;
        maxRef->node() = newNode;
        
        //construct the new tree by adding the chemPoints 
        //without using minRef and maxRef => test for maxId and minId
        //random access to the chemPoint to try to maintain depth O(logn)
        labelList chPIndex = identity(size_);//cellIndexTmp[i]=i
        Random randGenerator(unsigned(time(NULL)));
        label j;
        for (label i=0; i<size_; i++)
        {
            j=randGenerator.integer(i,size_-1);
            label tmp = chPIndex[i];
            chPIndex[i] = chPIndex[j];
            chPIndex[j] = tmp;
        }

        forAll(chemPoints,cpi)
        {
            if ((chPIndex[cpi]!=minId) && (chPIndex[cpi]!=maxId))
            {
                //search tree for position
                chemPointBase* phi0Base;
                binaryTreeSearch
                (
                    chemPoints[chPIndex[cpi]]->phi(),
                    root_,phi0Base
                );
                chP* phi0 = dynamic_cast<chP*>(phi0Base);
                //add the chemPoint
                bn* nodeToAdd =
                    new bn(phi0,chemPoints[chPIndex[cpi]], phi0->node());
                //make the parent of phi0 point to the newly created node
                insertNode(phi0, nodeToAdd);
                phi0->node()=nodeToAdd;
                chemPoints[chPIndex[cpi]]->node()=nodeToAdd;
            }
        }
        return true;
    }
    else
    { 
        return false;
    }
}


template<class CompType, class ThermoType>
Foam::chemPointISAT<CompType, ThermoType>*
Foam::binaryTree<CompType, ThermoType>::treeMin(bn* subTreeRoot)
{
    if (subTreeRoot!=NULL)
    {
        while(subTreeRoot->nodeLeft() != NULL)
        {
            subTreeRoot = subTreeRoot->nodeLeft();
        }
        return subTreeRoot->leafLeft();
    }
    else 
    {
        return NULL;
    }
}


template<class CompType, class ThermoType>
Foam::chemPointISAT<CompType, ThermoType>*
Foam::binaryTree<CompType, ThermoType>::treeSuccessor(chP* x)
{
    if (size_>1)
    {
        if (x==x->node()->leafLeft())
        {
            bn* parentNode = x->node();
            if (parentNode->nodeRight()==NULL)
            {
                return parentNode->leafRight();
            }
            else 
            {
                return treeMin(parentNode->nodeRight());
            }
        }
        else
        {
            bn* y = x->node();
            while((y->parent() !=NULL))
            {
                if (y==y->parent()->nodeLeft())
                {
                    if (y->parent()->nodeRight()==NULL)
                    {
                        return y->parent()->leafRight();
                    }
                    else 
                    {
                        return treeMin(y->parent()->nodeRight());
                    }
                }
                y=y->parent();
            }
            //when we reach this point, y points to the root and
            //never entered in the if loop (coming from the right)
            //so we are at the tree maximum and there is no successor
            return NULL; 
        }
    }
    else 
    {
        return NULL;
    }
}


// ************************************************************************* //
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOstream.H"
#include "dictionary.H"
#include "Switch.H"
#include "scalarField.H"
#include "chemPointISAT.H"
#include "binaryNode.H"
#include "TDACChemistryModel.H"
#include <limits>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
Foam::scalar Foam::chemPointISAT<CompType, ThermoType>::tolerance_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class CompType, class ThermoType>
Foam::chemPointISAT<CompType, ThermoType>::chemPointISAT
(
TDACChemistryModel<CompType, ThermoType>& chemistry,
const scalarField& phi,
const scalarField& Rphi,
const scalarSquareMatrix& A,
const scalarField& scaleFactor,
const scalar& tolerance,
const label& spaceSize,
const dictionary& coeffsDict,
binaryNode<CompType, ThermoType>* node
)
:
    chemistry_(&chemistry),
    phi_(phi),
    Rphi_(Rphi),
    A_(A),
    scaleFactor_(scaleFactor),
    node_(node),
    spaceSize_(spaceSize),
    nGrowth_(0),
    nActiveSpecies_(chemistry.mechRed()->NsSimp()),
    completeToSimplifiedIndex_(spaceSize-2),
    simplifiedToCompleteIndex_(nActiveSpecies_),
    inertSpecie_(-1),
    timeTag_(chemistry_->time().timeOutputValue()),
    lastTimeUsed_(chemistry_->time().timeOutputValue()),
    toRemove_(false),
    maxNumNewDim_(coeffsDict.lookupOrDefault("maxNumNewDim",0))
{
    tolerance_=tolerance;

    bool isMechRedActive = chemistry_->mechRed()->active();
    if (isMechRedActive)
    {
        for (label i=0; i<spaceSize-2; i++)
        {
            completeToSimplifiedIndex_[i] =
                chemistry.completeToSimplifiedIndex_[i];
        }
        for (label i=0; i<nActiveSpecies_; i++)
        {
            simplifiedToCompleteIndex_[i] =
                chemistry.simplifiedToCompleteIndex()[i];
        }
    }
    
    label dim = spaceSize;
    if (isMechRedActive)
    {
        dim = nActiveSpecies_+2;
    }
    
    LT_ = scalarSquareMatrix(dim, dim, 0.0);
    
    //SVD decomposition A= U*D*V^T 
    scalarSquareMatrix Atilde(A);
    scalarSquareMatrix B(dim,dim,0.0);
    scalarField diag(dim,0.0);
    svd(Atilde, dim-2, dim, diag, B);

    //replace the value of vector diag by max(diag, 1/2)
    for (label i=0; i<dim; i++)
    {
        diag[i]=max(diag[i], 0.5);
    }
    
    //reconstruct Atilde = U*D'*V (ellipsoid in with length d'[i] and principal
    //semi-axes in the direction of the column of )
    for (label i=0; i<dim-2; i++)
    {
        scalarField AtildeI(dim);
        for (label n=0; n<dim; n++)
        {
            AtildeI[n] = Atilde[i][n];
        }
        for (label j=0; j<dim; j++)
        {    
            Atilde[i][j]=0.0;
            for (label k=0; k<dim; k++)
            {
                Atilde[i][j] += AtildeI[k]*diag[k]*B[j][k];
            }
            //added to use qrDecompose on the reduced composition space
            label si=i;
            if (isMechRedActive)
            {
                si = simplifiedToCompleteIndex[i];
            }
            Atilde[i][j] /= (tolerance*scaleFactor[si]);
        }
    }

    qrDecompose(dim,Atilde);
    word inertSpecieName(chemistry.thermo().lookup("inertSpecie"));
    forAll(chemistry.Y(),Yi)
    {
        if (chemistry.Y()[Yi].name()==inertSpecieName)
        {
            inertSpecie_=Yi;
        }
    }
}


template<class CompType, class ThermoType>
Foam::chemPointISAT<CompType, ThermoType>::chemPointISAT
(
    Foam::chemPointISAT<CompType, ThermoType>& p
)
:
    chemPointBase(),
    phi_(p.phi()),
    Rphi_(p.Rphi()),
    LT_(p.LT()),
    A_(p.A()),
    scaleFactor_(p.scaleFactor()),
    node_(p.node()),
    spaceSize_(p.spaceSize()),
    nUsed_(p.nUsed()),
    nGrown_(p.nGrown()),
    nActiveSpecies_(p.nActiveSpecies()),
    completeToSimplifiedIndex_(p.completeToSimplifiedIndex()),
    simplifiedToCompleteIndex_(p.simplifiedToCompleteIndex()),
    inertSpecie_(p.inertSpecie()),
    timeTag_(p.timeTag()),
    lastTimeUsed_(p.lastTimeUsed()),
    toRemove_(p.toRemove()),
    maxNumNewDim_(p.maxNumNewDim())
{
   tolerance_ = p.tolerance();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
bool Foam::chemPointISAT<CompType, ThermoType>::inEOA(const scalarField& phiq)
{
    scalarField dphi=phiq-phi();
    bool isMechRedActive = chemistry_->mechRed()->active();
    label dim = (isMechRedActive) ? nActiveSpecies_ : spaceSize()-2;
    scalar epsTemp=0.0;
    
    if (dphi[spaceSize()-2]/phiq[spaceSize()-2]>0.05)
    {    
        return false;
    }

    for (label i=0; i<spaceSize()-2; i++)
    {
        //skip the inertSpecie
        if (i==inertSpecie_)
        {
            continue;
        }
        
        //When mechanism reduction is inactive OR on active species
        //multiply L by dphi to get the distance in the active species direction
        //else (for inactive species), just multiply the diagonal
        // element and dphi
        if
        (
            !(isMechRedActive)
          ||(isMechRedActive && completeToSimplifiedIndex_[i]!=-1)
        )
        {
            label si=(isMechRedActive) ? completeToSimplifiedIndex_[i] : i;
            for (label j=si; j<dim; j++)//LT is upper triangular
            {
                label sj=(isMechRedActive) ? simplifiedToCompleteIndex_[j] : j;
                epsTemp += LT_[si][j]*dphi[sj];
            }
            epsTemp += LT_[si][nActiveSpecies_]*dphi[spaceSize()-2];
            epsTemp += LT_[si][nActiveSpecies_+1]*dphi[spaceSize()-1];
        }
        else
        {
            epsTemp += dphi[i]/(tolerance_*scaleFactor_[i]);
        }
    }
    
    if (epsTemp > 1.0)
    {    
        return false;
    }
    else
    {   
        if (nUsed_ < INT_MAX)
        {
            nUsed_++;
        }
        return true;
    }
}


template<class CompType, class ThermoType>
bool Foam::chemPointISAT<CompType, ThermoType>::checkSolution
(
    const scalarField& phiq,
    const scalarField& Rphiq
)
{
    scalar eps2 = 0.0;
    scalarField dR = Rphiq - Rphi();
    scalarField dphi = phiq - phi();
    const scalarField& scaleFactorV = scaleFactor();
    const scalarSquareMatrix& Avar = A();
    bool isMechRedActive = chemistry_->mechRed()->active();
    scalar dRl = 0.0;
    label dim = spaceSize()-2;
    if (isMechRedActive)
    {
        dim = nActiveSpecies_;
    }
    
    for (label i=0; i<spaceSize()-2; i++)
    {
        if (i==inertSpecie_)
        {
            continue;
        }

        dRl = 0.0;
        if (isMechRedActive)
        {
            label si = completeToSimplifiedIndex_[i];
            if (si!=-1)
            {
                for (label j=0; j<dim; j++)
                {
                    label sj;
                    if (isMechRedActive)
                    {
                        sj =simplifiedToCompleteIndex_[j];
                    }
                    else
                    {
                        sj=j;
                    }
                    dRl += Avar[si][j]*dphi[sj];
                }
                dRl += Avar[si][nActiveSpecies_]*dphi[spaceSize()-2];
                dRl += Avar[si][nActiveSpecies_+1]*dphi[spaceSize()-1];
            }
            else
            {
                dRl = dphi[i];
            }
        }
        else
        {
            for (label j=0; j<spaceSize(); j++)
            {
                dRl += Avar[i][j]*dphi[j];
            }
        }
        eps2 += sqr((dR[i]-dRl)/scaleFactorV[i]);
    }
    eps2 = sqrt(eps2);
    if (eps2 > tolerance())
    {
        return false;
    }
    else
    {
        // if the solution is in the ellipsoid of accuracy
        // GROW operation performed
        if (grow(phiq)) //phiq is on the boundary of the EOA or grow has failed
        {
            nGrown_++;
            return true;
        }
        else
        {
            return false;
        }
    }
}


template<class CompType, class ThermoType>
bool Foam::chemPointISAT<CompType, ThermoType>::grow(const scalarField& phiq)
{
    scalarSquareMatrix& LTvar = LT();
    scalarSquareMatrix& Avar  = A();
    scalarField dphi = phiq - phi();
    label dim = spaceSize();
    label initNActiveSpecies(nActiveSpecies_);
    bool isMechRedActive = chemistry_->mechRed()->active();

    if (isMechRedActive)
    {
        label activeAdded(0);
        List<label> sAdded(spaceSize()-2);
        DynamicList<label> dimToAdd(0);

        //check if the difference of active species is lower than the maximum
        //number of new dimensions allowed
        for (label i=0; i<spaceSize()-2; i++)
        {
            //first test if the current chemPoint has an inactive species
            //corresponding to an active one in the query point
            if
            (
                completeToSimplifiedIndex_[i]==-1
             && chemistry_->completeToSimplifiedIndex_[i]!=-1
            )
            {
                activeAdded++;
                dimToAdd.append(i);
            }
            //then test if an active species in the current chemPoint
            //corresponds to an inactive on in the query
            if
            (
                completeToSimplifiedIndex_[i]!=-1
             && chemistry_->completeToSimplifiedIndex_[i]==-1
            )
            {
                activeAdded++;
                //we don't need to add a new dimension but we count it to have
                //control on the difference throuhg maxNumNewDim
            }
        }

        //if the number of added dimension is too large, growth fail
        if (activeAdded > maxNumNewDim_)
        {
            return false;
        }

        //should now reflect the number of added dimension to the current chemPoint
        activeAdded=0;
        forAll(dimToAdd,dimi)
        {
            label i(dimToAdd[dimi]);
            nActiveSpecies_++;
            //add the new active species
            simplifiedToCompleteIndex_.setSize(nActiveSpecies_,i);
            completeToSimplifiedIndex_[i]=nActiveSpecies_-1;
            sAdded[activeAdded++]=simplifiedToCompleteIndex_[nActiveSpecies_-1];
        }

        //update LT and A :
        //-add new column and line for the new active species
        //-transfer last two lines of the previous matrix (p and T) to the end
        //  (change the diagonal position)
        //-set all element of the new lines and columns to zero except diagonal
        //  (=1/(tolerance*scaleFactor))
        if (nActiveSpecies_ > initNActiveSpecies)
        {
            LTvar.setSize(nActiveSpecies_+2, List<scalar>(initNActiveSpecies+2,0.0));
            Avar.setSize(nActiveSpecies_+2, List<scalar>(initNActiveSpecies+2,0.0));
            for (label i=0; i<nActiveSpecies_+2; i++)
            {
                LTvar[i].setSize(nActiveSpecies_+2, 0.0);
                Avar[i].setSize(nActiveSpecies_+2, 0.0);
            }
            for (label i=0; i<nActiveSpecies_-activeAdded; i++)
            {
                //star with last column, otherwise problems when activeAdded=1
                for (label j=1; j>=0; j--)
                {
                    LTvar[i][nActiveSpecies_+j]=LTvar[i][nActiveSpecies_+j-activeAdded];
                    Avar[i][nActiveSpecies_+j]=Avar[i][nActiveSpecies_+j-activeAdded];
                    Avar[nActiveSpecies_+j][i]=Avar[nActiveSpecies_+j-activeAdded][i];
                    LTvar[i][nActiveSpecies_+j-activeAdded]=0.0;
                    Avar[i][nActiveSpecies_+j-activeAdded]=0.0;
                    Avar[nActiveSpecies_+j-activeAdded][i]=0.0;
                }
            }
            for (label i=nActiveSpecies_+1; i>=nActiveSpecies_; i--)
            {
                for (label j=nActiveSpecies_+1; j>=nActiveSpecies_; j--)
                {
                    LTvar[i][j]=LTvar[i-activeAdded][j-activeAdded];
                    Avar[i][j]=Avar[i-activeAdded][j-activeAdded];
                    LTvar[i-activeAdded][j-activeAdded]=0.0;
                    Avar[i-activeAdded][j-activeAdded]=0.0;
                }
            }
            for (label i=nActiveSpecies_-activeAdded; i<nActiveSpecies_;i++)
            {
                LTvar[i][i]=1.0/(tolerance_*scaleFactor_[simplifiedToCompleteIndex_[i]]);
                Avar[i][i]=1.0;
            }
        }//end if (nActiveSpecies_>initNActiveSpecies)
        dim = nActiveSpecies_+2;
    }//end if (isMechRedActive)
    //beginning of grow algorithm
    scalarField phiTilde(dim, 0.0);
    scalar normPhiTilde = 0.0;
    //p' = L^T.(p-phi)
    for (label i=0; i<dim; i++)
    {
        for (label j=i; j<dim-2; j++)//LT is upper triangular
        {
            label sj = j;
            if (isMechRedActive)
            {
                sj=simplifiedToCompleteIndex_[j];
            }
            phiTilde[i] += LTvar[i][j]*dphi[sj];
        }
        phiTilde[i] += LTvar[i][dim-2]*dphi[spaceSize()-2];
        phiTilde[i] += LTvar[i][dim-1]*dphi[spaceSize()-1];
        normPhiTilde += sqr(phiTilde[i]);
    }
    scalar invSqrNormPhiTilde = 1.0/normPhiTilde;
    normPhiTilde = sqrt(normPhiTilde);
    //gamma = (1/|p'| - 1)/|p'|^2
    scalar gamma = (1/normPhiTilde - 1)*invSqrNormPhiTilde;
    scalarField u(gamma*phiTilde);
    scalarField v(dim,0.0);
    for ( label i=0; i<dim; i++)
    {
        for (register label j=0; j<=i;j++)
        {
            v[i] += phiTilde[j]*LTvar[j][i];
        }
    }
    qrUpdate(dim, u, v);
    return true;
}


template<class CompType, class ThermoType>
void Foam::chemPointISAT<CompType, ThermoType>::qrDecompose
(
    const label nCols,
    scalarSquareMatrix& Q
)
{
    scalarField c(nCols);
    scalarField d(nCols);
    scalar scale, sigma, sum;
    
    for (label k=0; k<nCols-1; k++)
    {
        scale = 0.0;
        for (label i=k; i<nCols; i++)
        {
            scale=max(scale, fabs(Q[i][k]));
        }
        if (scale == 0.0)
        {
            c[k]=d[k]=0.0;
        }
        else
        {
            for (label i=k; i<nCols; i++)
            {
                Q[i][k] /= scale;
            }
            sum = 0.0;
            for (label i=k; i<nCols; i++)
            {
                sum += sqr(Q[i][k]);
            }
            sigma = sign(Q[k][k])*sqrt(sum);
            Q[k][k] += sigma;
            c[k]=sigma*Q[k][k];
            d[k]=-scale*sigma;
            for (label j=k+1; j<nCols; j++)
            {
                sum=0.0;
                for ( label i=k; i<nCols; i++)
                {
                    sum += Q[i][k]*Q[i][j];
                }
                scalar tau = sum/c[k];
                for ( label i=k; i<nCols; i++)
                {
                    Q[i][j] -= tau*Q[i][k];
                }
            }
        }
    }
    d[nCols-1] = Q[nCols-1][nCols-1];
    
    //form R
    scalarSquareMatrix& R(LT());
    for (label i=0; i<nCols; i++)
    {
        R[i][i] = d[i];
        for ( label j=0; j<i; j++)
        {
            R[i][j]=0.0;
        }
        for (label j=i+1; j<nCols; j++)
        {
            R[i][j]=Q[i][j];
        }
    }
}


template<class CompType, class ThermoType>
void Foam::chemPointISAT<CompType, ThermoType>::qrUpdate
(
    const label n,
    const Foam::scalarField &u,
    const Foam::scalarField &v
)
{
    label k,i;
    scalarField w(u);
    for (k=n-1;k>=0;k--)
    {
        if (w[k] != 0.0)
        {
            break;
        }
    }
    if (k < 0)
    {
        k=0;
    }
    for (i=k-1;i>=0;i--)
    {
        rotate(i,w[i],-w[i+1], n);
        if (w[i] == 0.0)
        {
            w[i]=fabs(w[i+1]);
        }
        else if (fabs(w[i]) > fabs(w[i+1]))
        {
            w[i]=fabs(w[i])*sqrt(1.0+sqr(w[i+1]/w[i]));
        }
        else
        {
            w[i]=fabs(w[i+1])*sqrt(1.0+sqr(w[i]/w[i+1]));
        }
    }
    scalarSquareMatrix&R(LT());
    for (i=0;i<n;i++)
    {
        R[0][i] += w[0]*v[i];
    }
    for (i=0;i<k;i++)
    {
        rotate(i,R[i][i],-R[i+1][i], n);
    }
}


//rotate function used by qrUpdate    
template<class CompType, class ThermoType>
void Foam::chemPointISAT<CompType, ThermoType>::rotate
(
    const label i,
    const scalar a,
    const scalar b,
    label n
)
{
    label j;
    scalar c, fact, s, w, y;
    if (a == 0.0)
    {
        c=0.0;
        s=(b >= 0.0 ? 1.0 : -1.0);
    }
    else if (fabs(a) > fabs(b))
    {
        fact = b/a;
        c=sign(a)/sqrt(1.0+(fact*fact));
        s=fact*c;
    }
    else
    {
        fact=a/b;
        s=sign(b)/sqrt(1.0+(fact*fact));
        c=fact*s;
    }
    scalarSquareMatrix& R(LT());
    for (j=i;j<n;j++)
    {
        y=R[i][j];
        w=R[i+1][j];
        R[i][j]=c*y-s*w;
        R[i+1][j]=s*y+c*w;
    }
}


template<class CompType, class ThermoType>
void Foam::chemPointISAT<CompType, ThermoType>::svd
(
    scalarSquareMatrix& A,
    label m,
    label n,
    scalarField& d,
    scalarSquareMatrix& V)
{
    //UPDATED VERSION NR3
    bool flag;
    label i,its,j,jj,k,l,nm;
    scalar anorm,c,f,g,h,s,scale,x,y,z;
    scalarField rv1(n);
    scalar eps = std::numeric_limits<scalar>::epsilon();
    g = scale = anorm = 0.0;

    //Householder reduction to bidiagonal form 
    for ( i = 0; i<n; i++)
    {
        l=i+2; //change from i+1 to i+2
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m)
        {
            for (k=i;k<m;k++)
            {
                scale += fabs(A[k][i]);
            }
            if (scale != 0.0)
            {
                for ( k=i;k<m;k++)
                {
                    A[k][i] /= scale;
                    s += A[k][i]*A[k][i];
                }
                f = A[i][i];
                g = -sign(f)*sqrt(s);
                h = f*g-s;
                A[i][i]=f-g;
                for (j=l-1;j<n;j++)
                {
                    for (s=0.0,k=i;k<m;k++)
                    {
                        s += A[k][i]*A[k][j];
                    }
                    f = s/h;
                    for (k=i; k<m;k++)
                    {
                        A[k][j] += f*A[k][i];
                    }
                }
                for (k=i; k<m;k++)
                {
                    A[k][i] *= scale;
                }
            }
        }
        d[i] = scale * g;
        g=s=scale=0.0;
        
        if (i+1 <= m && i+1 != n)
        {
            for (k=l-1; k<n; k++)
            {
                scale += fabs(A[i][k]);
            }
            if (scale != 0.0)
            {
                for (k=l-1; k<n; k++)
                {
                    A[i][k] /= scale;
                    s += A[i][k]*A[i][k];
                }
                f = A[i][l-1];
                g = -sign(f)*sqrt(s);
                h = f*g-s;
                A[i][l-1] = f-g;
                for (k=l-1; k<n; k++)
                {
                    rv1[k]=A[i][k]/h;
                }
                for (j=l-1; j<m; j++)
                {
                    for (s=0.0,k=l-1; k<n; k++)
                    {
                        s += A[j][k]*A[i][k];
                    }
                    for (k=l-1; k<n; k++)
                    {
                        A[j][k] += s*rv1[k];
                    }
                }
                for (k=l-1; k<n; k++)
                {
                    A[i][k] *= scale;
                }
            }
        }
        anorm = max(anorm, (fabs(d[i])+fabs(rv1[i])));
    }//end Householder reduction
    
    //Accumulation of right-hand transformations
    for (i=n-1; i>=0; i--)
    {
        if (i < n-1)
        {
            if (g != 0.0)
            {
                for (j=l; j<n; j++)
                {
                    V[j][i] = (A[i][j]/A[i][l])/g;
                }
                for (j=l; j<n; j++)
                {
                    for (s=0.0,k=l; k<n; k++)
                    {
                        s += A[i][k]*V[k][j];
                    }
                    for (k=l; k<n; k++)
                    {
                        V[k][j] += s*V[k][i];
                    }
                }
            }
            for (j=l; j<n; j++)
            {
                V[i][j]=V[j][i]=0.0;
            }
        }
        V[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    //Accumulation of left-hand transformations
    for (i = min(m,n)-1; i>=0; i--)
    {
        l=i+1;
        g=d[i];
        for (j=l; j<n; j++)
        {
            A[i][j] = 0.0;
        }
        if (g != 0.0)
        {
            g = 1.0/g;
            for (j=l; j<n; j++)
            {
                for (s=0.0, k=l; k<m; k++)
                {
                    s+= A[k][i]*A[k][j];
                }
                f = (s/A[i][i])*g;
                for (k=i; k<m; k++)
                {
                    A[k][j] += f*A[k][i];
                }
            }
            for (j=i; j<m; j++)
            {
                A[j][i] *= g;
            }
        }
        else
        {
            for (j=i; j<m; j++)
            {
                A[j][i]=0.0;
            }
        }
        ++A[i][i];
    }
    
    //Diagonalization of the bidiagonal form :
    //Loop over singular values, and over allowed iteration
    for (k=n-1; k>=0; k--)
    {
        for (its=0; its<30; its++)
        {
            flag=true;
            // Test for splitting (rv1[1] always zero)
            for (l=k; l>=0; l--)
            {
                nm = l-1;
                if (l == 0 || fabs(rv1[l]) <= eps*anorm)
                {
                    flag = false;
                    break;
                }
                if (fabs(d[nm]) <= eps*anorm)
                {
                    break;
                }
            }
            //Cancellation of rv1[l], if l>1
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i=l; i<k+1; i++)
                {
                    f = s*rv1[i];
                    rv1[i] = c*rv1[i];
                    if (fabs(f) <= eps*anorm)
                    {
                        break;
                    }
                    g = d[i];
                    h = pythag(f,g);
                    d[i] = h;
                    h = 1.0/h;
                    c = g*h;
                    s = -f*h;
                    for (j=0; j<m; j++)
                    {
                        y = A[j][nm];
                        z = A[j][i];
                        A[j][nm] = y*c + z*s;
                        A[j][i] = z*c - y*s;
                    }
                }
            }
            
            z = d[k];
            if (l == k) //Convergence
            {
                if (z < 0.0) //Singular value is made nonnegative
                {
                    d[k] = -z;
                    for (j=0; j<n; j++)
                    {
                        V[j][k] = -V[j][k];
                    }
                }
                break;
            }
            if (its == 29)
            {
                Info<< "No convergence in 30 iterations" << endl;
            }
            x = d[l];
            nm = k-1;
            y = d[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g = pythag(f,1.0);
            f = ((x-z)*(x+z)+h*((y/(f+sign(f)*g))-h))/x;
            c=s=1.0; 
            //Next QR transformation
            for (j=l; j<=nm; j++)
            {
                i = j+1;
                g = rv1[i];
                y = d[i];
                h = s*g;
                g = c*g;
                z = pythag(f,h);
                rv1[j] = z;
                c = f/z;
                s = h/z;
                f = x*c + g*s;
                g = g*c - x*s;
                h = y*s;
                y *= c;
                for (jj=0; jj<n; jj++)
                {
                    x = V[jj][j];
                    z = V[jj][i];
                    V[jj][j] = x*c + z*s;
                    V[jj][i] = z*c - x*s;
                }
                z = pythag(f,h);
                d[j] = z;
                if (z)
                {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }
                f = c*g + s*y;
                x = c*y - s*g;
                for (jj=0; jj<m; jj++)
                {
                    y = A[jj][j];
                    z = A[jj][i];
                    A[jj][j] = y*c + z*s;
                    A[jj][i] = z*c - y*s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            d[k] = x;
        }
    }
}//end svd


// pythag function used in svd
// compute (a^2+b^2)^1/2 without descrutive underflow or overflow
template<class CompType, class ThermoType>
Foam::scalar
Foam::chemPointISAT<CompType, ThermoType>::pythag(scalar a, scalar b)
{
    scalar absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb)
    {
        return absa*sqrt(1.0+sqr(absb/absa));
    }
    else
    {
        return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+sqr(absa/absb)));
    }
}

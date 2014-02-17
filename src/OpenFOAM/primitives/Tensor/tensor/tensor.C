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

#include "tensor.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* const tensor::typeName = "tensor";

    template<>
    const char* tensor::componentNames[] =
    {
        "xx", "xy", "xz",
        "yx", "yy", "yz",
        "zx", "zy", "zz"
    };

    template<>
    const tensor tensor::zero
    (
        0, 0, 0,
        0, 0, 0,
        0, 0, 0
    );

    template<>
    const tensor tensor::one
    (
        1, 1, 1,
        1, 1, 1,
        1, 1, 1
    );

    template<>
    const tensor tensor::max
    (
        VGREAT, VGREAT, VGREAT,
        VGREAT, VGREAT, VGREAT,
        VGREAT, VGREAT, VGREAT
    );

    template<>
    const tensor tensor::min
    (
        -VGREAT, -VGREAT, -VGREAT,
        -VGREAT, -VGREAT, -VGREAT,
        -VGREAT, -VGREAT, -VGREAT
    );

    template<>
    const tensor tensor::I
    (
        1, 0, 0,
        0, 1, 0,
        0, 0, 1
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::vector Foam::eigenValues(const tensor& t)
{
    scalar i = 0;
    scalar ii = 0;
    scalar iii = 0;

    if
    (
        (
            mag(t.xy()) + mag(t.xz()) + mag(t.yx())
          + mag(t.yz()) + mag(t.zx()) + mag(t.zy())
        )
      < SMALL
    )
    {
        // diagonal matrix
        i = t.xx();
        ii = t.yy();
        iii = t.zz();
    }
    else
    {
        scalar a = -t.xx() - t.yy() - t.zz();

        scalar b = t.xx()*t.yy() + t.xx()*t.zz() + t.yy()*t.zz()
            - t.xy()*t.yx() - t.xz()*t.zx() - t.yz()*t.zy();

        scalar c = - t.xx()*t.yy()*t.zz() - t.xy()*t.yz()*t.zx()
            - t.xz()*t.yx()*t.zy() + t.xz()*t.yy()*t.zx()
            + t.xy()*t.yx()*t.zz() + t.xx()*t.yz()*t.zy();

        // If there is a zero root
        if (mag(c) < ROOTVSMALL)
        {
            scalar disc = sqr(a) - 4*b;

            if (disc >= -SMALL)
            {
                scalar q = -0.5*sqrt(max(0.0, disc));

                i = 0;
                ii = -0.5*a + q;
                iii = -0.5*a - q;
            }
            else
            {
                FatalErrorIn("eigenValues(const tensor&)")
                    << "zero and complex eigenvalues in tensor: " << t
                    << abort(FatalError);
            }
        }
        else
        {
            scalar Q = (a*a - 3*b)/9;
            scalar R = (2*a*a*a - 9*a*b + 27*c)/54;

            scalar R2 = sqr(R);
            scalar Q3 = pow3(Q);

            // Three different real roots
            if (R2 < Q3)
            {
                scalar sqrtQ = sqrt(Q);
                scalar theta = acos(R/(Q*sqrtQ));

                scalar m2SqrtQ = -2*sqrtQ;
                scalar aBy3 = a/3;

                i = m2SqrtQ*cos(theta/3) - aBy3;
                ii = m2SqrtQ*cos((theta + twoPi)/3) - aBy3;
                iii = m2SqrtQ*cos((theta - twoPi)/3) - aBy3;
            }
            else
            {
                scalar A = cbrt(R + sqrt(R2 - Q3));

                // Three equal real roots
                if (A < SMALL)
                {
                    scalar root = -a/3;
                    return vector(root, root, root);
                }
                else
                {
                    // Complex roots
                    WarningIn("eigenValues(const tensor&)")
                        << "complex eigenvalues detected for tensor: " << t
                        << endl;

                    return vector::zero;
                }
            }
        }
    }


    // Sort the eigenvalues into ascending order
    if (i > ii)
    {
        Swap(i, ii);
    }

    if (ii > iii)
    {
        Swap(ii, iii);
    }

    if (i > ii)
    {
        Swap(i, ii);
    }

    return vector(i, ii, iii);
}


Foam::vector Foam::eigenVector(const tensor& t, const scalar lambda)
{
    if (mag(lambda) < SMALL)
    {
        return vector::zero;
    }

    // Construct the matrix for the eigenvector problem
    tensor A(t - lambda*I);

    // Calculate the sub-determinants of the 3 components
    scalar sd0 = A.yy()*A.zz() - A.yz()*A.zy();
    scalar sd1 = A.xx()*A.zz() - A.xz()*A.zx();
    scalar sd2 = A.xx()*A.yy() - A.xy()*A.yx();

    scalar magSd0 = mag(sd0);
    scalar magSd1 = mag(sd1);
    scalar magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        vector ev
        (
            1,
            (A.yz()*A.zx() - A.zz()*A.yx())/sd0,
            (A.zy()*A.yx() - A.yy()*A.zx())/sd0
        );
        ev /= mag(ev);

        return ev;
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        vector ev
        (
            (A.xz()*A.zy() - A.zz()*A.xy())/sd1,
            1,
            (A.zx()*A.xy() - A.xx()*A.zy())/sd1
        );
        ev /= mag(ev);

        return ev;
    }
    else if (magSd2 > SMALL)
    {
        vector ev
        (
            (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
            (A.yx()*A.xz() - A.xx()*A.yz())/sd2,
            1
        );
        ev /= mag(ev);

        return ev;
    }
    else
    {
        return vector::zero;
    }
}


Foam::tensor Foam::eigenVectors(const tensor& t)
{
    vector evals(eigenValues(t));

    tensor evs
    (
        eigenVector(t, evals.x()),
        eigenVector(t, evals.y()),
        eigenVector(t, evals.z())
    );

    return evs;
}


Foam::vector Foam::eigenValues(const symmTensor& t)
{
    scalar i = 0;
    scalar ii = 0;
    scalar iii = 0;

    if
    (
        (
            mag(t.xy()) + mag(t.xz()) + mag(t.xy())
          + mag(t.yz()) + mag(t.xz()) + mag(t.yz())
        )
      < SMALL
    )
    {
        // diagonal matrix
        i = t.xx();
        ii = t.yy();
        iii = t.zz();
    }
    else
    {
        scalar a = -t.xx() - t.yy() - t.zz();

        scalar b = t.xx()*t.yy() + t.xx()*t.zz() + t.yy()*t.zz()
            - t.xy()*t.xy() - t.xz()*t.xz() - t.yz()*t.yz();

        scalar c = - t.xx()*t.yy()*t.zz() - t.xy()*t.yz()*t.xz()
            - t.xz()*t.xy()*t.yz() + t.xz()*t.yy()*t.xz()
            + t.xy()*t.xy()*t.zz() + t.xx()*t.yz()*t.yz();

        // If there is a zero root
        if (mag(c) < ROOTVSMALL)
        {
            scalar disc = sqr(a) - 4*b;

            if (disc >= -SMALL)
            {
                scalar q = -0.5*sqrt(max(0.0, disc));

                i = 0;
                ii = -0.5*a + q;
                iii = -0.5*a - q;
            }
            else
            {
                FatalErrorIn("eigenValues(const tensor&)")
                    << "zero and complex eigenvalues in tensor: " << t
                    << abort(FatalError);
            }
        }
        else
        {
            scalar Q = (a*a - 3*b)/9;
            scalar R = (2*a*a*a - 9*a*b + 27*c)/54;

            scalar R2 = sqr(R);
            scalar Q3 = pow3(Q);

            // Three different real roots
            if (R2 < Q3)
            {
                scalar sqrtQ = sqrt(Q);
                scalar theta = acos(R/(Q*sqrtQ));

                scalar m2SqrtQ = -2*sqrtQ;
                scalar aBy3 = a/3;

                i = m2SqrtQ*cos(theta/3) - aBy3;
                ii = m2SqrtQ*cos((theta + twoPi)/3) - aBy3;
                iii = m2SqrtQ*cos((theta - twoPi)/3) - aBy3;
            }
            else
            {
                scalar A = cbrt(R + sqrt(R2 - Q3));

                // Three equal real roots
                if (A < SMALL)
                {
                    scalar root = -a/3;
                    return vector(root, root, root);
                }
                else
                {
                    // Complex roots
                    WarningIn("eigenValues(const symmTensor&)")
                        << "complex eigenvalues detected for symmTensor: " << t
                        << endl;

                    return vector::zero;
                }
            }
        }
    }


    // Sort the eigenvalues into ascending order
    if (i > ii)
    {
        Swap(i, ii);
    }

    if (ii > iii)
    {
        Swap(ii, iii);
    }

    if (i > ii)
    {
        Swap(i, ii);
    }

    return vector(i, ii, iii);
}


Foam::vector Foam::eigenVector(const symmTensor& t, const scalar lambda)
{
    if (mag(lambda) < SMALL)
    {
        return vector::zero;
    }

    // Construct the matrix for the eigenvector problem
    symmTensor A(t - lambda*I);

    // Calculate the sub-determinants of the 3 components
    scalar sd0 = A.yy()*A.zz() - A.yz()*A.yz();
    scalar sd1 = A.xx()*A.zz() - A.xz()*A.xz();
    scalar sd2 = A.xx()*A.yy() - A.xy()*A.xy();

    scalar magSd0 = mag(sd0);
    scalar magSd1 = mag(sd1);
    scalar magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        vector ev
        (
            1,
            (A.yz()*A.xz() - A.zz()*A.xy())/sd0,
            (A.yz()*A.xy() - A.yy()*A.xz())/sd0
        );
        ev /= mag(ev);

        return ev;
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        vector ev
        (
            (A.xz()*A.yz() - A.zz()*A.xy())/sd1,
            1,
            (A.xz()*A.xy() - A.xx()*A.yz())/sd1
        );
        ev /= mag(ev);

        return ev;
    }
    else if (magSd2 > SMALL)
    {
        vector ev
        (
            (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
            (A.xy()*A.xz() - A.xx()*A.yz())/sd2,
            1
        );
        ev /= mag(ev);

        return ev;
    }
    else
    {
        return vector::zero;
    }
}


Foam::tensor Foam::eigenVectors(const symmTensor& t)
{
    vector evals(eigenValues(t));

    tensor evs
    (
        eigenVector(t, evals.x()),
        eigenVector(t, evals.y()),
        eigenVector(t, evals.z())
    );

    return evs;
}


// ************************************************************************* //

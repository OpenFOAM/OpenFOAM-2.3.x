/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "SphereDragForce.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::SphereDragForce<CloudType>::CdRe(const scalar Re) const
{
    if (Re <= .01)
    {
    	return (9/2 + 24/Re) * Re;
    }

    else if (Re <= 20)
    {
    	return 24 * (1 + .1315 * pow(Re, .82 - .05 * log(Re)));
    }

    else if (Re <= 260)
    {
    	return 24 * (1 + .1935 * pow(Re, .6305));
    }

    else if (Re <= 1.5 * pow(10, 3))
    {
    	return Re * pow(1.6435 - 1.1242 * log(Re) + .1558 * pow(log(Re),2),10);
    }

    else if (Re <= 1.2 * pow(10, 4))
    {
    	return Re * pow(-2.4571 + 2.5558 * log(Re) - .9295 * pow(log(Re),2) 
    		            + .1049 * pow(log(Re),3),10);
    }

    else if (Re <= 4.4 * pow(10, 4))
    {
    	return Re * pow(-1.9181 + .6370 * log(Re) - .0636 * pow(log(Re),2),10);
    }

    else if (Re <= 3.38 * pow(10, 5))
    {
    	return Re * pow(-4.339 + 1.5809 * log(Re) - .1546 * pow(log(Re),2),10);
    }

    else if (Re <= 4 * pow(10, 5))
    {
    	return Re * (29.78 - 5.3 * log(Re));
    }

    else if (Re <= 1 * pow(10, 6))
    {
    	return Re * (.1 * log(Re) - .49);
    }

    else
    {
    	return Re * (.19 - 8 * pow(10, 4) / Re);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SphereDragForce<CloudType>::SphereDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
{}


template<class CloudType>
Foam::SphereDragForce<CloudType>::SphereDragForce
(
    const SphereDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SphereDragForce<CloudType>::~SphereDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::SphereDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(vector::zero, 0.0);

    value.Sp() = mass*0.75*muc*CdRe(Re)/(p.rho()*sqr(p.d()));

    return value;
}


// ************************************************************************* //

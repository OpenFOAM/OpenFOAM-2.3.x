/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2014 Francesco Contino
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

#include "mechanismReduction.H"
#include "Switch.H"
#include "error.H"
#include "TDACChemistryModel.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class CompType, class ThermoType>
Foam::mechanismReduction<CompType,ThermoType>::mechanismReduction
(
    const Foam::IOdictionary& dict,
    Foam::TDACChemistryModel<CompType,ThermoType>& chemistry
)
:
    dict_(dict),
    chemistry_(chemistry),
    activeSpecies_(chemistry.nSpecie(),false),
    NsSimp_(chemistry.nSpecie()),
    nSpecie_(chemistry.nSpecie()),
    coeffsDict_(dict.subDict("mechanismReduction")),
    tolerance_(readScalar(coeffsDict_.lookup("tolerance"))),
    active_(coeffsDict_.lookup("active"))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::mechanismReduction<CompType,ThermoType>::~mechanismReduction()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

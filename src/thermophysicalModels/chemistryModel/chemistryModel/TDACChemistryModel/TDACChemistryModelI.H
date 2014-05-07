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

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
inline const Foam::PtrList<Foam::volScalarField>&
Foam::TDACChemistryModel<CompType, ThermoType>::Y()
{
    return this->Y_;
}

template<class CompType, class ThermoType>
inline Foam::autoPtr<Foam::mechanismReduction<CompType, ThermoType> >&
Foam::TDACChemistryModel<CompType, ThermoType>::mechRed()
{
    return mechRed_;
}


template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::setActive(label i)
{
    this->Y_[i].writeOpt()=IOobject::AUTO_WRITE;
    activeSpecies_[i]=true;

}


template<class CompType, class ThermoType>
inline void
Foam::TDACChemistryModel<CompType, ThermoType>::setNsDAC(Foam::label newNsDAC)
{
    NsDAC_ = newNsDAC;
}


template<class CompType, class ThermoType>
inline void
Foam::TDACChemistryModel<CompType, ThermoType>::setNSpecie(Foam::label newNs)
{
    this->nSpecie_ = newNs;
}


template<class CompType, class ThermoType>
inline Foam::DynamicList<Foam::label>&
Foam::TDACChemistryModel<CompType, ThermoType>::simplifiedToCompleteIndex()
{
    return simplifiedToCompleteIndex_;
}


template<class CompType, class ThermoType>
inline Foam::Field<Foam::label>&
Foam::TDACChemistryModel<CompType, ThermoType>::completeToSimplifiedIndex()
{
    return completeToSimplifiedIndex_;
}


template<class CompType, class ThermoType>
inline const Foam::Field<Foam::label>&
Foam::TDACChemistryModel<CompType, ThermoType>::
completeToSimplifiedIndex() const
{
    return completeToSimplifiedIndex_;
}


template<class CompType, class ThermoType>
inline Foam::Field<bool>&
Foam::TDACChemistryModel<CompType, ThermoType>::reactionsDisabled()
{
    return reactionsDisabled_;
}


template<class CompType, class ThermoType>
inline Foam::scalarField&
Foam::TDACChemistryModel<CompType, ThermoType>::completeC()
{
    return completeC_;
}


template<class CompType, class ThermoType>
inline Foam::scalarField&
Foam::TDACChemistryModel<CompType, ThermoType>::simplifiedC()
{
    return simplifiedC_;
}


template<class CompType, class ThermoType>
bool Foam::TDACChemistryModel<CompType, ThermoType>::isActive(label i)
{
    return activeSpecies_[i];
}


template<class CompType, class ThermoType>
inline Foam::List<Foam::List<Foam::chemkinReader::specieElement> >&
Foam::TDACChemistryModel<CompType, ThermoType>::specieComp()
{
    return specieComp_;
}


// ************************************************************************* //
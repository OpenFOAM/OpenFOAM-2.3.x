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
inline Foam::autoPtr<Foam::mechanismReduction<CompType, ThermoType> >
Foam::TDACChemistryModel<CompType, ThermoType>::mechRed()
{
    return mechRed_;
}

/*
template<class CompType, class ThermoType>
inline const Foam::PtrList<Foam::Reaction<ThermoType> >&
Foam::TDACChemistryModel<CompType, ThermoType>::reactions() const
{
    return reactions_;
}


template<class CompType, class ThermoType>
inline const Foam::PtrList<ThermoType>&
Foam::TDACChemistryModel<CompType, ThermoType>::specieThermo() const
{
    return specieThermo_;
}


template<class CompType, class ThermoType>
inline const Foam::label
Foam::TDACChemistryModel<CompType, ThermoType>::nSpecie() const
{
    return nSpecie_;
}


template<class CompType, class ThermoType>
inline Foam::label
Foam::TDACChemistryModel<CompType, ThermoType>::nReaction() const
{
    return nReaction_;
}


template<class CompType, class ThermoType>
inline const Foam::chemistrySolverTDAC<CompType, ThermoType>&
Foam::TDACChemistryModel<CompType, ThermoType>::solver() const
{
    return solver_;
}


template<class CompType, class ThermoType>
inline Foam::scalar
Foam::TDACChemistryModel<CompType, ThermoType>::solveChemistryCpuTime()
{
    return solveChemistryCpuTime_;
}


template<class CompType, class ThermoType>
inline Foam::scalar
Foam::TDACChemistryModel<CompType, ThermoType>::reduceMechCpuTime()
{
    return reduceMechCpuTime_;
}


template<class CompType, class ThermoType>
inline Foam::scalar
Foam::TDACChemistryModel<CompType, ThermoType>::searchISATCpuTime()
{
    return searchISATCpuTime_;
}


template<class CompType, class ThermoType>
inline Foam::scalar
Foam::TDACChemistryModel<CompType, ThermoType>::addNewLeafCpuTime()
{
    return addNewLeafCpuTime_;
}


template<class CompType, class ThermoType>
inline Foam::label
Foam::TDACChemistryModel<CompType, ThermoType>::NsDAC() const
{
    return NsDAC_;
}


template<class CompType, class ThermoType>
inline Foam::label
Foam::TDACChemistryModel<CompType, ThermoType>::meanNsDAC() const
{
    return meanNsDAC_;
}


template<class CompType, class ThermoType>
inline Foam::Switch
Foam::TDACChemistryModel<CompType, ThermoType>::DAC() const
{
    return DAC_;
}

*/
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
/*
template<class CompType, class ThermoType>
inline Foam::tmp<Foam::volScalarField>
Foam::TDACChemistryModel<CompType, ThermoType>::RR
(
    const label i
) const
{
    tmp<volScalarField> tRR
    (
        new volScalarField
        (
            IOobject
            (
                "RR(" + this->Y_[i].name() + ')',
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        tRR().internalField() = RR_[i];
        tRR().correctBoundaryConditions();
    }
    return tRR;
}


template<class CompType, class ThermoType>
inline const Foam::label
Foam::TDACChemistryModel<CompType, ThermoType>::nFound() const
{
    return nFound_;
}


template<class CompType, class ThermoType>
inline const Foam::label
Foam::TDACChemistryModel<CompType, ThermoType>::nGrown() const
{
    return nGrown_;
}


template<class CompType, class ThermoType>
inline bool Foam::TDACChemistryModel<CompType, ThermoType>::exhaustiveSearch()
{
    return exhaustiveSearch_;
}



template<class CompType, class ThermoType>
inline Foam::DynamicList<Foam::scalar,1,1,1>&
Foam::TDACChemistryModel<CompType, ThermoType>::timeLabel()
{
    return timeLabel_;
}


template<class CompType, class ThermoType>
inline Foam::DynamicList<Foam::labelPair*,1,1,1>&
Foam::TDACChemistryModel<CompType, ThermoType>::edgesST()
{
    return edgesST_;
}


template<class CompType, class ThermoType>
inline Foam::DynamicList<Foam::DynamicList<Foam::label,1,1,1>*,1,1,1>&
Foam::TDACChemistryModel<CompType, ThermoType>::timeUsed()
{
    return timeUsed_;
}


template<class CompType, class ThermoType>
inline Foam::DynamicList<Foam::DynamicList<Foam::scalar,1,1,1>*,1,1,1>&
Foam::TDACChemistryModel<CompType, ThermoType>::edgeWeight()
{
    return edgeWeight_;
}


template<class CompType, class ThermoType>
inline Foam::List<Foam::scalar>&
Foam::TDACChemistryModel<CompType, ThermoType>::nodeWeight()
{
    return nodeWeight_;
}


template<class CompType, class ThermoType>
inline Foam::SquareMatrix<Foam::label>&
Foam::TDACChemistryModel<CompType, ThermoType>::speciesToEdges()
{
    return speciesToEdges_;
}


template<class CompType, class ThermoType>
inline bool Foam::TDACChemistryModel<CompType, ThermoType>::CADSwitch()
{
    return CADSwitch_;
}


template<class CompType, class ThermoType>
inline bool Foam::TDACChemistryModel<CompType, ThermoType>::cellSwitch()
{
    return cellSwitch_;
}


template<class CompType, class ThermoType>
inline Foam::List<bool>&
Foam::TDACChemistryModel<CompType, ThermoType>::nodeActivated()
{
    return nodeActivated_;
}


template<class CompType, class ThermoType>
inline Foam::List<Foam::label>&
Foam::TDACChemistryModel<CompType, ThermoType>::nodeTime()
{
    return nodeTime_;
}


template<class CompType, class ThermoType>
inline bool Foam::TDACChemistryModel<CompType, ThermoType>::isGraphWritten()
{
    return isGraphWritten_;
}
*/

// ************************************************************************* //
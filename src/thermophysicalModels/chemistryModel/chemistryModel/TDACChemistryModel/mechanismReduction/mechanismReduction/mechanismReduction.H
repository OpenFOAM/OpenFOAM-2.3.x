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
 
Class
    Foam::mechanismReduction

Description
    An abstract class for methods of chemical mechanism reduction

SourceFiles
    mechanismReduction.C

\*---------------------------------------------------------------------------*/

#ifndef mechanismReduction_H
#define mechanismReduction_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<class CompType, class ThermoType>
class TDACChemistryModel;

/*---------------------------------------------------------------------------*\
                           Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
class mechanismReduction
{


protected:

    const IOdictionary& dict_;
    TDACChemistryModel<CompType, ThermoType>& chemistry_;

    //List of active species (active = true)
    List<bool> activeSpecies_;
    //Number of active species
    label NsSimp_;
    //Number of species
    const label nSpecie_;

    //Dictionary that store the algorithm data
    const dictionary coeffsDict_;

private:

    //Tolerance for the mechanism reduction algorithm
    const scalar tolerance_;

    //Is mechanism reduction active?
    const Switch active_;

public:

    //- Runtime type information
    TypeName("mechanismReduction");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        mechanismReduction,
        dictionary,
        (
            const IOdictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        ),
        (dict, chemistry)
    );


    // Constructors

        //- Construct from components
        mechanismReduction
        (
            const IOdictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Selector

        static autoPtr<mechanismReduction<CompType, ThermoType> > New
        (
            const IOdictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Destructor

        virtual ~mechanismReduction();


    // Member Functions

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalarField &c,
            const scalar T,
            const scalar p
        )  = 0;
    
        //- Return the active species
        inline const List<bool>& activeSpecies() const;

        //- Return the number of active species
        inline label NsSimp();

        //- Return the initial number of species
        inline label nSpecie();

        //- Return the tolerance
        inline scalar tolerance() const;

        inline Switch active() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mechanismReductionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mechanismReduction.C"
#   include "mechanismReductionNew.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

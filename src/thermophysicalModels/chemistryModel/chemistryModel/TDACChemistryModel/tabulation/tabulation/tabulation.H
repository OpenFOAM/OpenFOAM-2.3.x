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


Class
    Foam::tabulation

Description
    An abstract class for tabulation

SourceFiles
    tabulation.C

\*---------------------------------------------------------------------------*/

#ifndef tabulation_H
#define tabulation_H

#include "IOdictionary.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "fileName.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
template<class CompType, class ThermoType>
class TDACChemistryModel;

/*---------------------------------------------------------------------------*\
                           Class tabulation Declaration
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
class tabulation
{


protected:

    const dictionary& dict_;
    TDACChemistryModel<CompType, ThermoType>& chemistry_;
    const dictionary coeffsDict_;
    const Switch active_;
    scalar tolerance_;

public:

    //- Runtime type information
    TypeName("tabulation");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        tabulation,
        dictionary,
        (
            const dictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        ),
        (dict, chemistry)
    );


    // Constructors

        //- Construct from components
        tabulation
        (
            const dictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Selectors

        static autoPtr<tabulation> New
        (
            const IOdictionary& dict,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Destructor

        virtual ~tabulation();


    // Member Functions

        inline Switch active()
        {
            return active_;
        }

      	inline scalar tolerance() const
        {
            return tolerance_;
        }

    // Virtual functions

        virtual label size() = 0;

        virtual void writePerformance() = 0;

        // Retrieve function: (only virtual here)
        // Try to retrieve a stored point close enough (according to tolerance)
        // to a stored point. If successful, it returns true and store the
        // results in RphiQ, i.e. the result of the integration of phiQ
        virtual bool retrieve
        (
             const scalarField& phiQ,
             scalarField& RphiQ
        ) = 0;

        // Add function: (only virtual here)
        // Add information to the tabulation algorithm. Give the reference for
        // future retrieve (phiQ) and the corresponding result (RphiQ).
        virtual bool add
        (
            const scalarField& phiQ,
            const scalarField& RphiQ,
            const scalar rho
        ) = 0;

        // Update function: (only virtual here)
        // The underlying structure of the tabulation is updated/cleaned
        // to increase the performance of the retrieve
        virtual bool update() = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tabulation.C"
#   include "tabulationNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::ISAT

Description
    Implementation of the ISAT (In-situ adaptive tabulation), for chemistry
    calculation. 
 
    For information:
    S.B. Pope (1997)
    ``Computationally efficient implementation of combustion chemistry using
    in situ adaptive tabulation,''
    Combustion Theory and Modelling,  1, 41--63.

\*---------------------------------------------------------------------------*/

#ifndef ISAT_H
#define ISAT_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tabulation.H"
#include "IOstream.H"
#include "dictionary.H"
#include "Switch.H"
#include "scalarField.H"
#include "binaryTree.H"
#include "Time.H"
#include "fileName.H"
#include "OFstream.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ISAT Declaration
\*---------------------------------------------------------------------------*/
template<class CompType, class ThermoType>
class ISAT
:
    public tabulation<CompType, ThermoType>
{

private:
    // Private data
        //- List of the stored 'points' organized in a binary tree
        binaryTree<CompType, ThermoType> chemisTree_;

        //- List of scale factors for species, temperature and pressure
        scalarField scaleFactor_;

        //- Lifetime (number of time steps) of a stored point
        const Time* runTime_;
        label chPMaxLifeTime_;

        //- Maximum number of growths before removing from the tree
        label maxGrowth_;

        //- Check the binary tree for leafs to remove every interval
        label checkEntireTreeInterval_;

        //- Factor that multiply the ideal depth of a binary tree to decide
        //wheter to try to balance of not
        scalar maxDepthFactor_;

        //- Minimal size before trying to balance the tree
        label minBalanceThreshold_;

        //- After a failed primary retrieve, look in the MRU list
        Switch MRURetrieve_;

        //- Most Recently Used (MRU) list of chemPoint
        SLList<chemPointISAT<CompType, ThermoType>*> MRUList_;

        //- Maximum size of the MRU list
        label maxMRUSize_;

        //- Store a pointer to the last chemPointISAT found
        chemPointISAT<CompType, ThermoType>* lastSearch_;

        //- Switch to allow growth (on by default)
        Switch growPoints_;

        //- Statistics on ISAT usage
        label nRetrieved_;
        label nGrowth_;
        label nAdd_;
        OFstream nRetrievedFile_;
        OFstream nGrowthFile_;
        OFstream nAddFile_;
        OFstream sizeFile_;

        bool cleaningRequired_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ISAT(const ISAT&);

        //- Add a chemPoint to the MRU list
        void addToMRU(chemPointISAT<CompType, ThermoType>* phi0);


        //- Compute and return the mapping of the composition phiq
        //  Input : phi0 the nearest chemPoint used in the linear interpolation
        //  phiq the composition of the query point for which we want to
        //  compute the mapping
        //  Rphiq the mapping of the new composition point (given as empty)
        //  Output: void (the mapping is stored in the Rphiq array)
        //  Rphiq = Rphi0 + A * (phiq-phi0)
        void calcNewC
        (
            chemPointISAT<CompType, ThermoType>* phi0,
            const scalarField& phiq,
            scalarField& Rphiq
        );


        //- Check if the composition of the query point phiq lies in the
        //- ellipsoid of accuracy approximating the region of accuracy of the
        //- stored chemPoint phi0
        //  Input : phi0 the nearest chemPoint used in the linear interpolation
        //  phiq the composition of the query point for which we want to
        //  compute the mapping
        //  Output: true if phiq is in the EOA, false if not
        bool grow
        (
            chemPointISAT<CompType, ThermoType>* phi0,
            const scalarField& phiq,
            const scalarField& Rphiq
        );

        //- Clean and balance the tree
        bool cleanAndBalance();

        //- Functions to construct the gradients matrix
        //  When mecanism reduction is active, the A matrix is given by
        //        Aaa Aad
        //  A = ( Ada Add ), where the sub gradient matrices are:
        //  (Aaa) active species according to active species, (Aad) active
        //  species according to disabled species, (Ada) disabled species
        //  according to active species, and (Add) disabled species according to
        //  disabled species.
        //  The current implementation computes Aaa with the Jacobian of the
        //  reduced set of species. Aad = 0, Ada = 0, and Add = I.
        //  To be implemented: add options to compute the A matrix for different
        //  strategies
        void computeA
        (
            scalarRectangularMatrix& A,
            const scalarField& Rphiq,
            const scalar rho
         );

        void gaussj
        (
            scalarRectangularMatrix& A,
            scalarRectangularMatrix& B, label n
        );
        void gaussj(scalarRectangularMatrix& A, label n);

public:

    //- Runtime type information
    TypeName("ISAT");

    // Constructors

        //- Construct from dictionary
        ISAT
        (
            const dictionary& chemistryProperties,
            TDACChemistryModel<CompType, ThermoType>& chemistry
        );


    // Destructor

        ~ISAT();


    // Member Functions

        // Access
            inline binaryTree<CompType, ThermoType>& chemisTree()
            {
                return chemisTree_;
            }

            inline const scalarField& scaleFactor() const
            {
                return scaleFactor_;
            }

            //- Find the closest stored leaf of phiQ and store the result in
            // RphiQ or return false.
            bool retrieve
            (
                const Foam::scalarField& phiq,
                scalarField& Rphiq
            );

        // Database

            //- Return the size of the binary tree
            inline label size()
            {
                return chemisTree_.size();
            }

            void writePerformance();

        // Edit


            //Add information to the tabulation.This function can grow an
            //existing point or add a new leaf to the binary tree
            //Input : phiq the new composition to store
            //Rphiq the mapping of the new composition point
            bool add
            (
                const scalarField& phiq,
                const scalarField& Rphiq,
                const scalar rho
            );

            inline bool update()
            {
                return cleanAndBalance();
            }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ISAT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

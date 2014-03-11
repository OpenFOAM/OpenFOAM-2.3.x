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
    Foam::TDACChemistryModel

Description
    Extends chemistryModel by adding the TDAC method.
    Ref:  
    F. Contino, H. Jeanmart, T. Lucchini, and G. D’Errico. Coupling of in situ 
    adaptive tabulation and dy- namic adaptive chemistry: An effective method 
    for solving combustion in engine simulations. Proceedings of the Combustion 
    Institute, 33(2):3057–3064, 2011

SourceFiles
    TDACChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef TDACChemistryModel_H
#define TDACChemistryModel_H

#include "chemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class TDACChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class TDACChemistryModel
:
    public chemistryModel<CompType, ThermoType>
{
    // Private Member Functions

        //- Disallow copy constructor
        TDACChemistryModel(const TDACChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const TDACChemistryModel&);

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


public:

    //- Runtime type information
    TypeName("TDACChemistryModel");


    // Constructors

        //- Construct from mesh
        TDACChemistryModel(const fvMesh& mesh);


    //- Destructor
    virtual ~TDACChemistryModel();


    // Member Functions

        //- dc/dt = omega, rate of change in concentration, for each species
        //  adapted to work with a variable number of species
        virtual tmp<scalarField> omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        //  adapted to work with a variable number of species
        virtual scalar omega
        (
            const Reaction<ThermoType>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;


        // Chemistry model functions (overriding functions in
        // chemistryModel.H to take into account the private solve function)

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);
    

        // ODE functions (overriding functions in chemistryModel.H to take into
        //  account the variable number of species)

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt,
                scalarSquareMatrix& dfdc
            ) const;

            virtual void solve
            (
                scalarField &c,
                scalar& T,
                scalar& p,
                scalar& deltaT,
                scalar& subDeltaT
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "TDACChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

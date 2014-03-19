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

#include "makeTabulationTypes.H"

#include "thermoPhysicsTypes.H"

#include "psiChemistryModel.H"
#include "rhoChemistryModel.H"

#include "ISAT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeTabulationTypes(psiChemistryModel, constGasHThermoPhysics);
    makeTabulationTypes(psiChemistryModel, gasHThermoPhysics);
    makeTabulationTypes
    (
     psiChemistryModel,
     constIncompressibleGasHThermoPhysics
     );
    makeTabulationTypes
    (
     psiChemistryModel,
     incompressibleGasHThermoPhysics)
    ;
    makeTabulationTypes(psiChemistryModel, icoPoly8HThermoPhysics);
    makeTabulationTypes(rhoChemistryModel, constGasHThermoPhysics);
    makeTabulationTypes(rhoChemistryModel, gasHThermoPhysics);
    makeTabulationTypes
    (
     rhoChemistryModel,
     constIncompressibleGasHThermoPhysics
     );
    makeTabulationTypes
    (
     rhoChemistryModel,
     incompressibleGasHThermoPhysics
     );
    makeTabulationTypes(rhoChemistryModel, icoPoly8HThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy
    makeTabulationTypes(psiChemistryModel, constGasEThermoPhysics);
    makeTabulationTypes(psiChemistryModel, gasEThermoPhysics);
    makeTabulationTypes
    (
     psiChemistryModel,
     constIncompressibleGasEThermoPhysics
     );
    makeTabulationTypes
    (
     psiChemistryModel,
     incompressibleGasEThermoPhysics
     );
    makeTabulationTypes(psiChemistryModel, icoPoly8EThermoPhysics);
    makeTabulationTypes(rhoChemistryModel, constGasEThermoPhysics);
    makeTabulationTypes(rhoChemistryModel, gasEThermoPhysics);
    makeTabulationTypes
    (
     rhoChemistryModel,
     constIncompressibleGasEThermoPhysics
     );
    makeTabulationTypes
    (
     rhoChemistryModel,
     incompressibleGasEThermoPhysics
     );
    makeTabulationTypes(rhoChemistryModel, icoPoly8EThermoPhysics);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

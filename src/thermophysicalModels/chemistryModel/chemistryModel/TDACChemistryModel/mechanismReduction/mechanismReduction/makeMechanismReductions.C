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

#include "makeMechanismReductionTypes.H"

#include "thermoPhysicsTypes.H"

#include "psiChemistryModel.H"
#include "rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeMechanismReductionTypes(psiChemistryModel, constGasHThermoPhysics);
    makeMechanismReductionTypes(psiChemistryModel, gasHThermoPhysics);
    makeMechanismReductionTypes
    (
     psiChemistryModel,
     constIncompressibleGasHThermoPhysics
     );
    makeMechanismReductionTypes
    (
     psiChemistryModel,
     incompressibleGasHThermoPhysics)
    ;
    makeMechanismReductionTypes(psiChemistryModel, icoPoly8HThermoPhysics);
    makeMechanismReductionTypes(rhoChemistryModel, constGasHThermoPhysics);
    makeMechanismReductionTypes(rhoChemistryModel, gasHThermoPhysics);
    makeMechanismReductionTypes
    (
     rhoChemistryModel,
     constIncompressibleGasHThermoPhysics
     );
    makeMechanismReductionTypes
    (
     rhoChemistryModel,
     incompressibleGasHThermoPhysics
     );
    makeMechanismReductionTypes(rhoChemistryModel, icoPoly8HThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy
    makeMechanismReductionTypes(psiChemistryModel, constGasEThermoPhysics);
    makeMechanismReductionTypes(psiChemistryModel, gasEThermoPhysics);
    makeMechanismReductionTypes
    (
     psiChemistryModel,
     constIncompressibleGasEThermoPhysics
     );
    makeMechanismReductionTypes
    (
     psiChemistryModel,
     incompressibleGasEThermoPhysics
     );
    makeMechanismReductionTypes(psiChemistryModel, icoPoly8EThermoPhysics);
    makeMechanismReductionTypes(rhoChemistryModel, constGasEThermoPhysics);
    makeMechanismReductionTypes(rhoChemistryModel, gasEThermoPhysics);
    makeMechanismReductionTypes
    (
     rhoChemistryModel,
     constIncompressibleGasEThermoPhysics
     );
    makeMechanismReductionTypes
    (
     rhoChemistryModel,
     incompressibleGasEThermoPhysics
     );
    makeMechanismReductionTypes(rhoChemistryModel, icoPoly8EThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

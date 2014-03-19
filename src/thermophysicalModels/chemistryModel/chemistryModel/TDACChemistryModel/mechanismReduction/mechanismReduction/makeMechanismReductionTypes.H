/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2014 Francesco Contino
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

#ifndef makeMechanismReductionTypes_H
#define makeMechanismReductionTypes_H

#include "mechanismReduction.H"


#include "DAC.H"
#include "DRG.H"
#include "DRGEP.H"
#include "EFA.H"
#include "PFA.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeMechanismReductionType(SS, Comp, Thermo)                          \
                                                                              \
    typedef SS<Comp, Thermo> SS##Comp##Thermo;                                \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        SS##Comp##Thermo,                                                     \
        (#SS"<" + word(Comp::typeName_())                                     \
      + "," + Thermo::typeName() + ">").c_str(),                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    mechanismReduction<Comp, Thermo>::                                        \
        adddictionaryConstructorToTable<SS<Comp, Thermo> >                    \
            add##SS##Comp##Thermo##ConstructorToTable_;


#define makeMechanismReductionTypes(CompChemModel,Thermo)                     \
                                                                              \
    typedef mechanismReduction<CompChemModel, Thermo>                         \
        mechanismReduction##CompChemModel##Thermo;                            \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        mechanismReduction##CompChemModel##Thermo,                            \
        "mechanismReduction<"#CompChemModel","#Thermo">",                     \
        0                                                                     \
    );                                                                        \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        mechanismReduction##CompChemModel##Thermo,                            \
        dictionary                                                            \
    );                                                                        \
                                                                              \
    makeMechanismReductionType                                                \
    (                                                                         \
        DAC,                                                                  \
        CompChemModel,                                                        \
        Thermo                                                                \
    );                                                                        \
                                                                              \
    makeMechanismReductionType                                                \
    (                                                                         \
        DRG,                                                                  \
        CompChemModel,                                                        \
        Thermo                                                                \
    );                                                                        \
                                                                              \
    makeMechanismReductionType                                                \
    (                                                                         \
        DRGEP,                                                                \
        CompChemModel,                                                        \
        Thermo                                                                \
    );                                                                        \
                                                                              \
    makeMechanismReductionType                                                \
    (                                                                         \
        EFA,                                                                  \
        CompChemModel,                                                        \
        Thermo                                                                \
    );                                                                        \
                                                                              \
    makeMechanismReductionType                                                \
    (                                                                         \
        PFA,                                                                  \
        CompChemModel,                                                        \
        Thermo                                                                \
    );                                                                        \

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

#ifndef makeTabulationTypes_H
#define makeTabulationTypes_H

#include "tabulation.H"

#include "ISAT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeTabulationType(SS, Comp, Thermo)                                  \
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
    tabulation<Comp, Thermo>::                                                \
        adddictionaryConstructorToTable<SS<Comp, Thermo> >                    \
            add##SS##Comp##Thermo##ConstructorToTable_;


#define makeTabulationTypes(CompChemModel,Thermo)                             \
                                                                              \
    typedef tabulation<CompChemModel, Thermo>                                 \
        tabulation##CompChemModel##Thermo;                                    \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        tabulation##CompChemModel##Thermo,                                    \
        "tabulation<"#CompChemModel","#Thermo">",                             \
        0                                                                     \
    );                                                                        \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        tabulation##CompChemModel##Thermo,                                    \
        dictionary                                                            \
    );                                                                        \
                                                                              \
    makeTabulationType                                                        \
    (                                                                         \
        ISAT,                                                                 \
        CompChemModel,                                                        \
        Thermo                                                                \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

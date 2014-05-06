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

#include "tabulation.H"


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
Foam::autoPtr<Foam::tabulation<CompType, ThermoType> >
Foam::tabulation<CompType, ThermoType>::New
(
    const IOdictionary& dict,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
{
    IOdictionary thermoDict
    (
         IOobject
         (
              "thermophysicalProperties",
              dict.path(),
              dict.db(),
              IOobject::MUST_READ_IF_MODIFIED,
              IOobject::NO_WRITE,
              false
         )
     );

    word thermoTypeName;

    if (thermoDict.isDict("thermoType"))
    {
        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));
        thermoTypeName =
            word(thermoTypeDict.lookup("transport")) + '<'
          + word(thermoTypeDict.lookup("thermo")) + '<'
          + word(thermoTypeDict.lookup("equationOfState")) + '<'
          + word(thermoTypeDict.lookup("specie")) + ">>,"
          + word(thermoTypeDict.lookup("energy")) + ">";
    }
    else
    {
        FatalIOErrorIn
        (
             (tabulation::typeName + "::New(const mesh&)").c_str(),
             thermoDict
        )<< "thermoType is in the old format and must be upgraded"
         << exit(FatalIOError);
    }

    dictionary tabdict(dict.subDict("tabulation"));

    word tabulationTypeName =
        word(tabdict.lookup("tabulationMethod")) + '<'
      + word(dict.subDict("chemistryType").lookup("chemistryThermo")) + ','
      + thermoTypeName + '>';
    
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(tabulationTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "tabulation::New(const dictionary&, const chemistryModel&)"
        )   << "Unknown tabulationType type " << tabulationTypeName
            << endl << endl
            << "Valid tabulationType types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }
    return autoPtr<tabulation<CompType, ThermoType> >
        (cstrIter()(dict, chemistry));
}


// ************************************************************************* //

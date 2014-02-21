/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Perfect gas equation of state.

\*---------------------------------------------------------------------------*/

#include "PengRobinsonGas.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
PengRobinsonGas<Specie>::PengRobinsonGas(Istream& is)
:
    Specie(is),
    Tc_(readScalar(is)),
    Pc_(readScalar(is)),
    omega_(readScalar(is))
{
    is.check("PengRobinsonGas<Specie>::PengRobinsonGas(Istream& is)");
}


template<class Specie>
Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const dictionary& dict
)
:
    Specie(dict),
    Tc_(readScalar(dict.subDict("equationOfState").lookup("Tc"))),
    Pc_(readScalar(dict.subDict("equationOfState").lookup("Pc"))),
    omega_(readScalar(dict.subDict("equationOfState").lookup("omega")))
{}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
PengRobinsonGas<Specie> operator+
(
    const PengRobinsonGas<Specie>& pg1,
    const PengRobinsonGas<Specie>& pg2
)
{
    scalar nMoles = pg1.nMoles() + pg2.nMoles();
    scalar molr1 = pg1.nMoles()/nMoles;
    scalar molr2 = pg2.nMoles()/nMoles;

    return PengRobinsonGas<Specie>
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2),
        molr1*pg1.Tc_ + molr2*pg2.Tc_,
        molr1*pg1.Pc_ + molr2*pg2.Pc_,
        molr1*pg1.omega_ + molr2*pg2.omega_
    );
}


template<class Specie>
PengRobinsonGas<Specie> operator-
(
    const PengRobinsonGas<Specie>& pg1,
    const PengRobinsonGas<Specie>& pg2
)
{
    scalar nMoles = pg1.nMoles() + pg2.nMoles();
    scalar molr1 = pg1.nMoles()/nMoles;
    scalar molr2 = pg2.nMoles()/nMoles;

    return PengRobinsonGas<Specie>
    (
        static_cast<const Specie&>(pg1)
      - static_cast<const Specie&>(pg2),
        molr1*pg1.Tc_ - molr2*pg2.Tc_,
        molr1*pg1.Pc_ - molr2*pg2.Pc_,
        molr1*pg1.omega_ - molr2*pg2.omega_
    );
}


template<class Specie>
PengRobinsonGas<Specie> operator*
(
    const scalar s,
    const PengRobinsonGas<Specie>& pg
)
{
    return PengRobinsonGas<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.Tc_,
        pg.Pc_,
        pg.omega_
    );
}


template<class Specie>
PengRobinsonGas<Specie> operator==
(
    const PengRobinsonGas<Specie>& pg1,
    const PengRobinsonGas<Specie>& pg2
)
{
    return pg2 - pg1;
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Ostream& operator<<
(
    Ostream& os,
    const PengRobinsonGas<Specie>& pg
)
{
    os  << static_cast<const Specie&>(pg)
        << token::SPACE << pg.Tc_
        << token::SPACE << pg.Pc_
        << token::SPACE << pg.omega_;

    os.check
    (
        "Ostream& operator<<(Ostream& os, const PengRobinsonGas<Specie>& st)"
    );
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

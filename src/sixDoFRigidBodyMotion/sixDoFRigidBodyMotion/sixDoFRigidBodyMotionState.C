/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "sixDoFRigidBodyMotionState.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionState::sixDoFRigidBodyMotionState()
:
    centreOfMass_(vector::zero),
    Q_(I),
    v_(vector::zero),
    a_(vector::zero),
    pi_(vector::zero),
    tau_(vector::zero)
{}


Foam::sixDoFRigidBodyMotionState::sixDoFRigidBodyMotionState
(
    const point& centreOfMass,
    const tensor& Q,
    const vector& v,
    const vector& a,
    const vector& pi,
    const vector& tau
)
:
    centreOfMass_(centreOfMass),
    Q_(Q),
    v_(v),
    a_(a),
    pi_(pi),
    tau_(tau)
{}


Foam::sixDoFRigidBodyMotionState::sixDoFRigidBodyMotionState
(
    const dictionary& dict
)
:
    centreOfMass_(dict.lookup("centreOfMass")),
    Q_(dict.lookupOrDefault("orientation", tensor::I)),
    v_(dict.lookupOrDefault("velocity", vector::zero)),
    a_(dict.lookupOrDefault("acceleration", vector::zero)),
    pi_(dict.lookupOrDefault("angularMomentum", vector::zero)),
    tau_(dict.lookupOrDefault("torque", vector::zero))
{}


Foam::sixDoFRigidBodyMotionState::sixDoFRigidBodyMotionState
(
    const sixDoFRigidBodyMotionState& sDoFRBMS
)
:
    centreOfMass_(sDoFRBMS.centreOfMass()),
    Q_(sDoFRBMS.Q()),
    v_(sDoFRBMS.v()),
    a_(sDoFRBMS.a()),
    pi_(sDoFRBMS.pi()),
    tau_(sDoFRBMS.tau())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionState::~sixDoFRigidBodyMotionState()
{}


// ************************************************************************* //

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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::tensor
Foam::sixDoFRigidBodyMotion::rotationTensorX(scalar phi) const
{
    return tensor
    (
        1, 0, 0,
        0, Foam::cos(phi), -Foam::sin(phi),
        0, Foam::sin(phi), Foam::cos(phi)
    );
}


inline Foam::tensor
Foam::sixDoFRigidBodyMotion::rotationTensorY(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), 0, Foam::sin(phi),
        0, 1, 0,
        -Foam::sin(phi), 0, Foam::cos(phi)
    );
}


inline Foam::tensor
Foam::sixDoFRigidBodyMotion::rotationTensorZ(scalar phi) const
{
    return tensor
    (
        Foam::cos(phi), -Foam::sin(phi), 0,
        Foam::sin(phi), Foam::cos(phi), 0,
        0, 0, 1
    );
}


inline Foam::Tuple2<Foam::tensor, Foam::vector>
Foam::sixDoFRigidBodyMotion::rotate
(
    const tensor& Q0,
    const vector& pi0,
    const scalar deltaT
) const
{
    Tuple2<tensor, vector> Qpi(Q0, pi0);
    tensor& Q = Qpi.first();
    vector& pi = Qpi.second();

    tensor R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorY(0.5*deltaT*pi.y()/momentOfInertia_.yy());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorZ(deltaT*pi.z()/momentOfInertia_.zz());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorY(0.5*deltaT*pi.y()/momentOfInertia_.yy());
    pi = pi & R;
    Q = Q & R;

    R = rotationTensorX(0.5*deltaT*pi.x()/momentOfInertia_.xx());
    pi = pi & R;
    Q = Q & R;

    return Qpi;
}


inline const Foam::sixDoFRigidBodyMotionState&
Foam::sixDoFRigidBodyMotion::state() const
{
    return motionState_;
}


inline const Foam::PtrList<Foam::sixDoFRigidBodyMotionRestraint>&
Foam::sixDoFRigidBodyMotion::restraints() const
{
    return restraints_;
}


inline const Foam::PtrList<Foam::sixDoFRigidBodyMotionConstraint>&
Foam::sixDoFRigidBodyMotion::constraints() const
{
    return constraints_;
}


inline const Foam::point&
Foam::sixDoFRigidBodyMotion::initialCentreOfMass() const
{
    return initialCentreOfMass_;
}


inline const Foam::tensor&
Foam::sixDoFRigidBodyMotion::initialQ() const
{
    return initialQ_;
}


inline const Foam::tensor& Foam::sixDoFRigidBodyMotion::Q() const
{
    return motionState_.Q();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::v() const
{
    return motionState_.v();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::a() const
{
    return motionState_.a();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::pi() const
{
    return motionState_.pi();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::tau() const
{
    return motionState_.tau();
}


inline const Foam::tensor& Foam::sixDoFRigidBodyMotion::Q0() const
{
    return motionState0_.Q();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::v0() const
{
    return motionState0_.v();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::a0() const
{
    return motionState0_.a();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::pi0() const
{
    return motionState0_.pi();
}


inline const Foam::vector& Foam::sixDoFRigidBodyMotion::tau0() const
{
    return motionState0_.tau();
}


inline Foam::point& Foam::sixDoFRigidBodyMotion::initialCentreOfMass()
{
    return initialCentreOfMass_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyMotion::initialQ()
{
    return initialQ_;
}


inline Foam::tensor& Foam::sixDoFRigidBodyMotion::Q()
{
    return motionState_.Q();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotion::v()
{
    return motionState_.v();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotion::a()
{
    return motionState_.a();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotion::pi()
{
    return motionState_.pi();
}


inline Foam::vector& Foam::sixDoFRigidBodyMotion::tau()
{
    return motionState_.tau();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::point Foam::sixDoFRigidBodyMotion::currentPosition
(
    const point& initialPoints
) const
{
    return
    (
        centreOfMass()
      + (Q() & initialQ_.T() & (initialPoints - initialCentreOfMass_))
    );
}


inline Foam::vector Foam::sixDoFRigidBodyMotion::currentOrientation
(
    const vector& vInitial
) const
{
    return (Q() & initialQ_.T() & vInitial);
}


inline const Foam::tensor&
Foam::sixDoFRigidBodyMotion::orientation() const
{
    return Q();
}


inline Foam::vector Foam::sixDoFRigidBodyMotion::omega() const
{
    return  Q() & (inv(momentOfInertia_) & pi());
}


inline Foam::point Foam::sixDoFRigidBodyMotion::currentVelocity
(
    const point& pt
) const
{
    return (omega() ^ (pt - centreOfMass())) + v();
}


inline const Foam::point& Foam::sixDoFRigidBodyMotion::centreOfMass() const
{
    return motionState_.centreOfMass();
}


inline const Foam::point& Foam::sixDoFRigidBodyMotion::centreOfMass0() const
{
    return motionState0_.centreOfMass();
}


inline const Foam::diagTensor&
Foam::sixDoFRigidBodyMotion::momentOfInertia() const
{
    return momentOfInertia_;
}


inline Foam::scalar Foam::sixDoFRigidBodyMotion::mass() const
{
    return mass_;
}


inline bool Foam::sixDoFRigidBodyMotion::report() const
{
    return report_;
}


inline void Foam::sixDoFRigidBodyMotion::newTime()
{
    motionState0_ = motionState_;
}


inline Foam::point& Foam::sixDoFRigidBodyMotion::centreOfMass()
{
    return motionState_.centreOfMass();
}


inline Foam::diagTensor& Foam::sixDoFRigidBodyMotion::momentOfInertia()
{
    return momentOfInertia_;
}


inline Foam::scalar& Foam::sixDoFRigidBodyMotion::mass()
{
    return mass_;
}

// ************************************************************************* //

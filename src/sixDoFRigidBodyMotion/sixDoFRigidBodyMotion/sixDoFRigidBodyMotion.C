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

#include "sixDoFRigidBodyMotion.H"
#include "septernion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::applyRestraints()
{
    if (restraints_.empty())
    {
        return;
    }

    if (Pstream::master())
    {
        forAll(restraints_, rI)
        {
            if (report_)
            {
                Info<< "Restraint " << restraints_[rI].name() << ": ";
            }

            // restraint position
            point rP = vector::zero;

            // restraint force
            vector rF = vector::zero;

            // restraint moment
            vector rM = vector::zero;

            restraints_[rI].restrain(*this, rP, rF, rM);

            a() += rF/mass_;

            // Moments are returned in global axes, transforming to
            // body local to add to torque.
            tau() += Q().T() & (rM + ((rP - centreOfMass()) ^ rF));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion()
:
    motionState_(),
    motionState0_(),
    restraints_(),
    constraints_(),
    tConstraints_(tensor::I),
    rConstraints_(tensor::I),
    initialCentreOfMass_(vector::zero),
    initialQ_(I),
    momentOfInertia_(diagTensor::one*VSMALL),
    mass_(VSMALL),
    aRelax_(1.0),
    aDamp_(1.0),
    report_(false)
{}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const point& centreOfMass,
    const tensor& Q,
    const vector& v,
    const vector& a,
    const vector& pi,
    const vector& tau,
    scalar mass,
    const point& initialCentreOfMass,
    const tensor& initialQ,
    const diagTensor& momentOfInertia,
    scalar aRelax,
    scalar aDamp,
    bool report
)
:
    motionState_
    (
        centreOfMass,
        Q,
        v,
        a,
        pi,
        tau
    ),
    motionState0_(motionState_),
    restraints_(),
    constraints_(),
    tConstraints_(tensor::I),
    rConstraints_(tensor::I),
    initialCentreOfMass_(initialCentreOfMass),
    initialQ_(initialQ),
    momentOfInertia_(momentOfInertia),
    mass_(mass),
    aRelax_(aRelax),
    aDamp_(aDamp),
    report_(report)
{}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const dictionary& dict,
    const dictionary& stateDict
)
:
    motionState_(stateDict),
    motionState0_(motionState_),
    restraints_(),
    constraints_(),
    tConstraints_(tensor::I),
    rConstraints_(tensor::I),
    initialCentreOfMass_
    (
        dict.lookupOrDefault
        (
            "initialCentreOfMass",
            vector(dict.lookup("centreOfMass"))
        )
    ),
    initialQ_
    (
        dict.lookupOrDefault
        (
            "initialOrientation",
            dict.lookupOrDefault("orientation", tensor::I)
        )
    ),
    momentOfInertia_(dict.lookup("momentOfInertia")),
    mass_(readScalar(dict.lookup("mass"))),
    aRelax_(dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0)),
    aDamp_(dict.lookupOrDefault<scalar>("accelerationDamping", 1.0)),
    report_(dict.lookupOrDefault<Switch>("report", false))
{
    addRestraints(dict);
    addConstraints(dict);
}


Foam::sixDoFRigidBodyMotion::sixDoFRigidBodyMotion
(
    const sixDoFRigidBodyMotion& sDoFRBM
)
:
    motionState_(sDoFRBM.motionState_),
    motionState0_(sDoFRBM.motionState0_),
    restraints_(sDoFRBM.restraints_),
    constraints_(sDoFRBM.constraints_),
    initialCentreOfMass_(sDoFRBM.initialCentreOfMass_),
    initialQ_(sDoFRBM.initialQ_),
    momentOfInertia_(sDoFRBM.momentOfInertia_),
    mass_(sDoFRBM.mass_),
    aRelax_(sDoFRBM.aRelax_),
    aDamp_(sDoFRBM.aDamp_),
    report_(sDoFRBM.report_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotion::~sixDoFRigidBodyMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotion::addRestraints
(
    const dictionary& dict
)
{
    if (dict.found("restraints"))
    {
        const dictionary& restraintDict = dict.subDict("restraints");

        label i = 0;

        restraints_.setSize(restraintDict.size());

        forAllConstIter(IDLList<entry>, restraintDict, iter)
        {
            if (iter().isDict())
            {
                restraints_.set
                (
                    i++,
                    sixDoFRigidBodyMotionRestraint::New
                    (
                        iter().keyword(),
                        iter().dict()
                    )
                );
            }
        }

        restraints_.setSize(i);
    }
}


void Foam::sixDoFRigidBodyMotion::addConstraints
(
    const dictionary& dict
)
{
    if (dict.found("constraints"))
    {
        const dictionary& constraintDict = dict.subDict("constraints");

        label i = 0;

        constraints_.setSize(constraintDict.size());

        pointConstraint pct;
        pointConstraint pcr;

        forAllConstIter(IDLList<entry>, constraintDict, iter)
        {
            if (iter().isDict())
            {
                constraints_.set
                (
                    i,
                    sixDoFRigidBodyMotionConstraint::New
                    (
                        iter().keyword(),
                        iter().dict()
                    )
                );

                constraints_[i].constrainTranslation(pct);
                constraints_[i].constrainRotation(pcr);

                i++;
            }
        }

        constraints_.setSize(i);

        tConstraints_ = pct.constraintTransformation();
        rConstraints_ = pcr.constraintTransformation();

        Info<< "Translational constraint tensor " << tConstraints_ << nl
            << "Rotational constraint tensor " << rConstraints_ << endl;
    }
}


void Foam::sixDoFRigidBodyMotion::updatePosition
(
    scalar deltaT,
    scalar deltaT0
)
{
    // First leapfrog velocity adjust and motion part, required before
    // force calculation

    if (Pstream::master())
    {
        v() = tConstraints_ & aDamp_*(v0() + 0.5*deltaT0*a());
        pi() = rConstraints_ & aDamp_*(pi0() + 0.5*deltaT0*tau());

        // Leapfrog move part
        centreOfMass() = centreOfMass0() + deltaT*v();

        // Leapfrog orientation adjustment
        Tuple2<tensor, vector> Qpi = rotate(Q0(), pi(), deltaT);
        Q() = Qpi.first();
        pi() = rConstraints_ & Qpi.second();
    }

    Pstream::scatter(motionState_);
}


void Foam::sixDoFRigidBodyMotion::updateAcceleration
(
    const vector& fGlobal,
    const vector& tauGlobal,
    scalar deltaT
)
{
    static bool first = false;

    // Second leapfrog velocity adjust part, required after motion and
    // acceleration calculation

    if (Pstream::master())
    {
        // Save the previous iteration accelerations for relaxation
        vector aPrevIter = a();
        vector tauPrevIter = tau();

        // Calculate new accelerations
        a() = fGlobal/mass_;
        tau() = (Q().T() & tauGlobal);
        applyRestraints();

        // Relax accelerations on all but first iteration
        if (!first)
        {
            a() = aRelax_*a() + (1 - aRelax_)*aPrevIter;
            tau() = aRelax_*tau() + (1 - aRelax_)*tauPrevIter;
        }
        first = false;

        // Correct velocities
        v() += tConstraints_ & aDamp_*0.5*deltaT*a();
        pi() += rConstraints_ & aDamp_*0.5*deltaT*tau();

        if (report_)
        {
            status();
        }
    }

    Pstream::scatter(motionState_);
}


void Foam::sixDoFRigidBodyMotion::updateVelocity(scalar deltaT)
{
    // Second leapfrog velocity adjust part, required after motion and
    // acceleration calculation

    if (Pstream::master())
    {
        v() += tConstraints_ & aDamp_*0.5*deltaT*a();
        pi() += rConstraints_ & aDamp_*0.5*deltaT*tau();

        if (report_)
        {
            status();
        }
    }

    Pstream::scatter(motionState_);
}


void Foam::sixDoFRigidBodyMotion::updateAcceleration
(
    const pointField& positions,
    const vectorField& forces,
    scalar deltaT
)
{
    vector fGlobal = vector::zero;

    vector tauGlobal = vector::zero;

    if (Pstream::master())
    {
        fGlobal = sum(forces);

        forAll(positions, i)
        {
            tauGlobal += (positions[i] - centreOfMass()) ^ forces[i];
        }
    }

    updateAcceleration(fGlobal, tauGlobal, deltaT);
}


Foam::point Foam::sixDoFRigidBodyMotion::predictedPosition
(
    const point& pInitial,
    const vector& deltaForce,
    const vector& deltaMoment,
    scalar deltaT
) const
{
    vector vTemp = v() + deltaT*(a() + deltaForce/mass_);
    vector piTemp = pi() + deltaT*(tau() + (Q().T() & deltaMoment));
    point centreOfMassTemp = centreOfMass0() + deltaT*vTemp;
    Tuple2<tensor, vector> QpiTemp = rotate(Q0(), piTemp, deltaT);

    return
    (
        centreOfMassTemp
      + (QpiTemp.first() & initialQ_.T() & (pInitial - initialCentreOfMass_))
    );
}


Foam::vector Foam::sixDoFRigidBodyMotion::predictedOrientation
(
    const vector& vInitial,
    const vector& deltaMoment,
    scalar deltaT
) const
{
    vector piTemp = pi() + deltaT*(tau() + (Q().T() & deltaMoment));
    Tuple2<tensor, vector> QpiTemp = rotate(Q0(), piTemp, deltaT);

    vector o(QpiTemp.first() & initialQ_.T() & vInitial);
    o /= mag(o);

    return o;
}


void Foam::sixDoFRigidBodyMotion::status() const
{
    Info<< "Centre of mass: " << centreOfMass() << nl
        << "Linear velocity: " << v() << nl
        << "Angular velocity: " << omega()
        << endl;
}


Foam::tmp<Foam::pointField> Foam::sixDoFRigidBodyMotion::currentPosition
(
    const pointField& initialPoints
) const
{
    return
    (
        centreOfMass()
      + (Q() & initialQ_.T() & (initialPoints - initialCentreOfMass_))
    );
}


Foam::tmp<Foam::pointField> Foam::sixDoFRigidBodyMotion::scaledPosition
(
    const pointField& initialPoints,
    const scalarField& scale
) const
{
    // Calculate the transformation septerion from the initial state
    septernion s
    (
        centreOfMass() - initialCentreOfMass(),
        quaternion(Q() & initialQ().T())
    );

    tmp<pointField> tpoints(new pointField(initialPoints));
    pointField& points = tpoints();

    forAll(points, pointi)
    {
        // Move non-stationary points
        if (scale[pointi] > SMALL)
        {
            // Use solid-body motion where scale = 1
            if (scale[pointi] > 1 - SMALL)
            {
                points[pointi] = currentPosition(initialPoints[pointi]);
            }
            // Slerp septernion interpolation
            else
            {
                septernion ss(slerp(septernion::I, s, scale[pointi]));

                points[pointi] =
                    initialCentreOfMass()
                  + ss.transform(initialPoints[pointi] - initialCentreOfMass());
            }
        }
    }

    return tpoints;
}


bool Foam::sixDoFRigidBodyMotion::read(const dictionary& dict)
{
    dict.lookup("momentOfInertia") >> momentOfInertia_;
    dict.lookup("mass") >> mass_;
    aRelax_ = dict.lookupOrDefault<scalar>("accelerationRelaxation", 1.0);
    aDamp_ = dict.lookupOrDefault<scalar>("accelerationDamping", 1.0);
    report_ = dict.lookupOrDefault<Switch>("report", false);

    restraints_.clear();
    addRestraints(dict);

    constraints_.clear();
    addConstraints(dict);

    return true;
}


// ************************************************************************* //

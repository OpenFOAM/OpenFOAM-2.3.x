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

Class
    Foam::sixDoFRigidBodyMotion

Description
    Six degree of freedom motion for a rigid body.  Angular momentum
    stored in body fixed reference frame.  Reference orientation of
    the body (where Q = I) must align with the cartesian axes such
    that the Inertia tensor is in principle component form.

    Symplectic motion as per:

    title = {Symplectic splitting methods for rigid body molecular dynamics},
    publisher = {AIP},
    year = {1997},
    journal = {The Journal of Chemical Physics},
    volume = {107},
    number = {15},
    pages = {5840-5851},
    url = {http://link.aip.org/link/?JCP/107/5840/1},
    doi = {10.1063/1.474310}

    Can add restraints (i.e. a spring) and constraints (i.e. motion
    may only be on a plane).

SourceFiles
    sixDoFRigidBodyMotionI.H
    sixDoFRigidBodyMotion.C
    sixDoFRigidBodyMotionIO.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyMotion_H
#define sixDoFRigidBodyMotion_H

#include "sixDoFRigidBodyMotionState.H"
#include "pointField.H"
#include "sixDoFRigidBodyMotionRestraint.H"
#include "sixDoFRigidBodyMotionConstraint.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class sixDoFRigidBodyMotion;
Istream& operator>>(Istream&, sixDoFRigidBodyMotion&);
Ostream& operator<<(Ostream&, const sixDoFRigidBodyMotion&);


/*---------------------------------------------------------------------------*\
                      Class sixDoFRigidBodyMotion Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyMotion
{
    // Private data

        //- Motion state data object
        sixDoFRigidBodyMotionState motionState_;

        //- Motion state data object for previous time-step
        sixDoFRigidBodyMotionState motionState0_;

        //- Motion restraints
        PtrList<sixDoFRigidBodyMotionRestraint> restraints_;

        //- Motion constaints
        PtrList<sixDoFRigidBodyMotionConstraint> constraints_;

        //- Translational constraint tensor
        tensor tConstraints_;

        //- Rotational constraint tensor
        tensor rConstraints_;

        //- Centre of mass of initial state
        point initialCentreOfMass_;

        //- Orientation of initial state
        tensor initialQ_;

        //- Moment of inertia of the body in reference configuration
        //  (Q = I)
        diagTensor momentOfInertia_;

        //- Mass of the body
        scalar mass_;

        //- Acceleration relaxation coefficient
        scalar aRelax_;

        //- Acceleration damping coefficient (for steady-state simulations)
        scalar aDamp_;

        //- Switch to turn reporting of motion data on and off
        Switch report_;


    // Private Member Functions

        //- Calculate the rotation tensor around the body reference
        //  frame x-axis by the given angle
        inline tensor rotationTensorX(scalar deltaT) const;

        //- Calculate the rotation tensor around the body reference
        //  frame y-axis by the given angle
        inline tensor rotationTensorY(scalar deltaT) const;

        //- Calculate the rotation tensor around the body reference
        //  frame z-axis by the given angle
        inline tensor rotationTensorZ(scalar deltaT) const;

        //- Apply rotation tensors to Q0 for the given torque (pi) and deltaT
        //  and return the rotated Q and pi as a tuple
        inline Tuple2<tensor, vector> rotate
        (
            const tensor& Q0,
            const vector& pi,
            const scalar deltaT
        ) const;

        //- Apply the restraints to the object
        void applyRestraints();


        // Access functions retained as private because of the risk of
        // confusion over what is a body local frame vector and what is global

        // Access

            //- Return access to the restraints
            inline const PtrList<sixDoFRigidBodyMotionRestraint>&
                restraints() const;

            //- Return access to the constraints
            inline const PtrList<sixDoFRigidBodyMotionConstraint>&
                constraints() const;

            //- Return access to the initial centre of mass
            inline const point& initialCentreOfMass() const;

            //- Return access to the initial orientation
            inline const tensor& initialQ() const;

            //- Return access to the orientation
            inline const tensor& Q() const;

            //- Return access to acceleration
            inline const vector& a() const;

            //- Return access to angular momentum
            inline const vector& pi() const;

            //- Return access to torque
            inline const vector& tau() const;

            //- Return access to the orientation at previous time-step
            inline const tensor& Q0() const;

            //- Return access to velocity at previous time-step
            inline const vector& v0() const;

            //- Return access to acceleration at previous time-step
            inline const vector& a0() const;

            //- Return access to angular momentum at previous time-step
            inline const vector& pi0() const;

            //- Return access to torque at previous time-step
            inline const vector& tau0() const;


        // Edit

            //- Return access to the centre of mass
            inline point& initialCentreOfMass();

            //- Return access to the centre of mass
            inline tensor& initialQ();

            //- Return non-const access to the orientation
            inline tensor& Q();

            //- Return non-const access to vector
            inline vector& v();

            //- Return non-const access to acceleration
            inline vector& a();

            //- Return non-const access to angular momentum
            inline vector& pi();

            //- Return non-const access to torque
            inline vector& tau();


public:

    // Constructors

        //- Construct null
        sixDoFRigidBodyMotion();

        //- Construct from components
        sixDoFRigidBodyMotion
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
            scalar aRelax = 1.0,
            scalar aDamp = 1.0,
            bool report = false
        );

        //- Construct from constant and state dictionaries
        sixDoFRigidBodyMotion
        (
            const dictionary& dict,
            const dictionary& stateDict
        );

        //- Construct as copy
        sixDoFRigidBodyMotion(const sixDoFRigidBodyMotion&);


    //- Destructor
    ~sixDoFRigidBodyMotion();


    // Member Functions

        //- Add restraints to the motion, public to allow external
        //  addition of restraints after construction
        void addRestraints(const dictionary& dict);

        //- Add restraints to the motion, public to allow external
        //  addition of restraints after construction
        void addConstraints(const dictionary& dict);

        //- First leapfrog velocity adjust and motion part, required
        //  before force calculation.  Takes old timestep for variable
        //  timestep cases.
        void updatePosition
        (
            scalar deltaT,
            scalar deltaT0
        );

        //- Second leapfrog velocity adjust part
        //  required after motion and force calculation
        void updateAcceleration
        (
            const vector& fGlobal,
            const vector& tauGlobal,
            scalar deltaT
        );

        //- Second leapfrog velocity adjust part
        //  required after motion and force calculation
        void updateVelocity(scalar deltaT);

        //- Global forces supplied at locations, calculating net force
        //  and moment
        void updateAcceleration
        (
            const pointField& positions,
            const vectorField& forces,
            scalar deltaT
        );

        //- Transform the given initial state pointField by the current
        //  motion state
        tmp<pointField> currentPosition(const pointField& initialPoints) const;

        //- Transform the given initial state pointField by the current
        //  motion state
        tmp<pointField> scaledPosition
        (
            const pointField& initialPoints,
            const scalarField& scale
        ) const;

        //- Transform the given initial state point by the current motion
        //  state
        inline point currentPosition(const point& initialPoints) const;

        //- Transform the given initial state direction by the current
        //  motion state
        inline vector currentOrientation(const vector& vInitial) const;

        //- Access the orientation tensor, Q.
        //  globalVector = Q & bodyLocalVector
        //  bodyLocalVector = Q.T() & globalVector
        inline const tensor& orientation() const;

        //- Predict the position of the supplied initial state point
        //  after deltaT given the current motion state and the
        //  additional supplied force and moment
        point predictedPosition
        (
            const point& initialPoints,
            const vector& deltaForce,
            const vector& deltaMoment,
            scalar deltaT
        ) const;

        //- Predict the orientation of the supplied initial state
        //  vector after deltaT given the current motion state and the
        //  additional supplied moment
        vector predictedOrientation
        (
            const vector& vInitial,
            const vector& deltaMoment,
            scalar deltaT
        ) const;

        //- Return the angular velocity in the global frame
        inline vector omega() const;

        //- Return the velocity of a position
        inline point currentVelocity(const point& pt) const;

        //- Report the status of the motion
        void status() const;


        // Access

            //- Return access to the motion state
            inline const sixDoFRigidBodyMotionState& state() const;

            //- Return const access to the centre of mass
            inline const point& centreOfMass() const;

            //- Return const access to the centre of mass at previous time-step
            inline const point& centreOfMass0() const;

            //- Return access to the inertia tensor
            inline const diagTensor& momentOfInertia() const;

            //- Return const access to the mass
            inline scalar mass() const;

            //- Return the report Switch
            inline bool report() const;

            //- Return access to velocity
            inline const vector& v() const;


        // Edit

            //- Store the motion state at the beginning of the time-step
            inline void newTime();

            //- Return non-const access to the centre of mass
            inline point& centreOfMass();

            //- Return non-const access to the inertia tensor
            inline diagTensor& momentOfInertia();

            //- Return non-const access to the mass
            inline scalar& mass();


        //- Write
        void write(Ostream&) const;


        //- Read coefficients dictionary and update system parameters,
        //  constraints and restraints but not the current state
        bool read(const dictionary& dict);


    // IOstream Operators

        friend Istream& operator>>(Istream&, sixDoFRigidBodyMotion&);
        friend Ostream& operator<<(Ostream&, const sixDoFRigidBodyMotion&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDoFRigidBodyMotionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

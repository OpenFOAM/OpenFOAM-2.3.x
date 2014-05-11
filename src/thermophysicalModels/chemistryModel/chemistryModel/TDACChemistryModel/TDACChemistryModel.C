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

#include "TDACChemistryModel.H"
#include "reactingMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::TDACChemistryModel<CompType, ThermoType>::TDACChemistryModel
(
    const fvMesh& mesh
)
:
    chemistryModel<CompType, ThermoType>(mesh),
    NsDAC_(this->nSpecie_),
    completeC_(this->nSpecie_,0.0),
    reactionsDisabled_(this->reactions_.size(), false),
    activeSpecies_(this->nSpecie_,false),
    completeToSimplifiedIndex_(this->nSpecie_,-1),
    simplifiedToCompleteIndex_(this->nSpecie_),
    specieComp_(this->nSpecie_)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
         )
     );

    // Store the species composition according to the species index
    speciesTable speciesTab = this->thermo().composition().species();
    chemkinReader tchemRead(thermoDict, speciesTab);
    const HashTable<List<chemkinReader::specieElement> >& specComp =
    tchemRead.specieComposition();
    forAll(specieComp_,i)
    {
        specieComp_[i] = specComp[this->Y()[i].name()];
    }

    mechRed_ =
        mechanismReduction<CompType, ThermoType>::New
        (
            *this,
            *this
        );

    tabulation_ =
        tabulation<CompType, ThermoType>::New
        (
            *this,
            *this
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::TDACChemistryModel<CompType, ThermoType>::~TDACChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
Foam::tmp<Foam::scalarField>
Foam::TDACChemistryModel<CompType, ThermoType>::omega
(
    const scalarField& c,
    const scalar T,
    const scalar p
) const
{
    scalar pf, cf, pr, cr;
    label lRef, rRef;

    //tom will have a reduced size when mechanism reduction is active
    //because nEqns() takes into account the reduced set of species
    tmp<scalarField> tom(new scalarField(this->nEqns(), 0.0));
    scalarField& om = tom();

    //However we need a vector of the complete set of species for the
    //third-body reactions
    scalarField c2(completeC_.size(), 0.0);
    if (mechRed_->active())
    {
        c2 = completeC_;
        //Update the concentration of the species in the simplified mechanism
        //the other species remain the same and are used only for third-body
        //efficiencies
        for(label i=0; i<NsDAC_; i++)
        {
            c2[simplifiedToCompleteIndex_[i]] = max(0.0, c[i]);
        }
    }
    else
    {
        for(label i=0; i<this->nSpecie(); i++)
        {
            c2[i] = max(0.0, c[i]);
        }
    }

    forAll(this->reactions_, i)
    {
        if (!reactionsDisabled_[i])
        {
            const Reaction<ThermoType>& R = this->reactions_[i];

            scalar omegai = omega
            (
                R, c2, T, p, pf, cf, lRef, pr, cr, rRef
            );

            forAll(R.lhs(), s)
            {
                label si = R.lhs()[s].index;
                if (mechRed_->active())
                {
                    si = completeToSimplifiedIndex_[si];
                }
                const scalar sl = R.lhs()[s].stoichCoeff;
                om[si] -= sl*omegai;
            }
            forAll(R.rhs(), s)
            {
                label si = R.rhs()[s].index;
                if (mechRed_->active())
                {
                    si = completeToSimplifiedIndex_[si];
                }
                const scalar sr = R.rhs()[s].stoichCoeff;
                om[si] += sr*omegai;
            }
        }
    }
    return tom;
}



template<class CompType, class ThermoType>
Foam::scalar Foam::TDACChemistryModel<CompType, ThermoType>::omega
(
    const Reaction<ThermoType>& R,
    const scalarField& c,//contains all species even when mechRed is active
    const scalar T,
    const scalar p,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    scalarField c2(completeC_.size(), 0.0);
    for (label i = 0; i < completeC_.size(); i++)
    {
        c2[i] = max(0.0, c[i]);
    }

    const scalar kf = R.kf(p, T, c2);
    const scalar kr = R.kr(kf, p, T, c2);

    pf = 1.0;
    pr = 1.0;

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s = 1; s < Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(0.0, c[lRef]), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(0.0, c[si]), exp);
        }
    }
    cf = max(0.0, c[lRef]);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1.0)
        {
            if (cf > SMALL)
            {
                pf *= pow(cf, exp - 1.0);
            }
            else
            {
                pf = 0.0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1.0);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // find the matrix element and element position for the rhs
    pr = kr;
    for (label s = 1; s < Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(0.0, c[rRef]), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(0.0, c[si]), exp);
        }
    }
    cr = max(0.0, c[rRef]);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1.0)
        {
            if (cr>SMALL)
            {
                pr *= pow(cr, exp - 1.0);
            }
            else
            {
                pr = 0.0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1.0);
        }
    }

    return pf*cf - pr*cr;
}


template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::derivatives
(
    const scalar time,
    const scalarField &c,
    scalarField& dcdt
) const
{
    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    tmp<scalarField> tom(omega(c, T, p));
    for (label i=0; i<this->nEqns(); i++)
    {
        dcdt[i] = tom()[i];
    }

    scalarField c2(completeC_.size(), 0.0);
    if (mechRed_->active())
    {
        //when using DAC, the ODE solver submit a reduced set of species
        //the complete set is used and only the species in the simplified
        //mechanism are updated
        c2 = completeC_;

        //update the concentration of the species in the simplified mechanism
        //the other species remain the same and are used only for third-body
        //efficiencies
        for(label i=0; i<NsDAC_; i++)
        {
            c2[simplifiedToCompleteIndex_[i]] = max(0.0, c[i]);
        }
    }
    else
    {
        for(label i=0; i<this->nSpecie(); i++)
        {
            c2[i] = max(0.0, c[i]);
        }
    }

    // constant pressure
    // dT/dt = ...
    scalar rho = 0.0;
    for (label i = 0; i < c2.size(); i++)
    {
        const scalar W = this->specieThermo_[i].W();
        rho += W*c2[i];
    }
    scalar cp = 0.0;
    for (label i=0; i<c2.size(); i++)
    {
        //cp function returns [J/(kmol K)]
        cp += c2[i]*this->specieThermo_[i].cp(p, T);
    }
    cp /= rho;
    scalar dT = 0.0;
    //when mechanism reduction is active
    //dT is computed on the reduced set since dcdt is null
    //for species not involved in the simplified mechanism
    for (label i = 0; i < this->nSpecie_; i++)
    {
        label si;
        if (mechRed_->active())
        {
            si = simplifiedToCompleteIndex_[i];
        }
        else
        {
            si = i;
        }
        //ha function returns [J/kmol]
        const scalar hi = this->specieThermo_[si].ha(p, T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;
    dcdt[this->nSpecie_] = -dT;

    // dp/dt = ...
    dcdt[this->nSpecie_ + 1] = 0.0;
}


template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    scalarField& dcdt,
    scalarSquareMatrix& dfdc
) const
{
    //if the mechanism reduction is active, the computed Jacobian
    //is compact (size of the reduced set of species)
    //but according to the informations of the complete set
    //(i.e. for the third-body efficiencies)
    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    scalarField c2(completeC_.size(), 0.0);
    if (mechRed_->active())
    {
        c2 = completeC_;
        for(label i=0; i<NsDAC_; i++)
        {
            c2[simplifiedToCompleteIndex_[i]] = max(0.0, c[i]);
        }
    }
    else
    {
        forAll(c2, i)
        {
            c2[i] = max(c[i], 0.0);
        }
    }

    for (label i=0; i<this->nEqns(); i++)
    {
        for (label j=0; j<this->nEqns(); j++)
        {
            dfdc[i][j] = 0.0;
        }
    }
    // length of the first argument must be nSpecie()
    // (reduced when mechanism reduction is active)
    tmp<scalarField> tom(omega(c, T, p));
    for (label i=0; i<this->nEqns(); i++)
    {
        dcdt[i] = tom()[i];
    }

    forAll(this->reactions_, ri)
    {
        if (!reactionsDisabled_[ri])
        {
            const Reaction<ThermoType>& R = this->reactions_[ri];
            const scalar kf0 = R.kf(p, T, c2);
            const scalar kr0 = R.kr(kf0, p, T, c2);
            forAll(R.lhs(), j)
            {
                label sj = R.lhs()[j].index;
                if (mechRed_->active())
                {
                    sj = completeToSimplifiedIndex_[sj];
                }
                scalar kf = kf0;
                forAll(R.lhs(), i)
                {
                    const label si = R.lhs()[i].index;
                    const scalar el = R.lhs()[i].exponent;
                    if (i == j)
                    {
                        if (el < 1.0)
                        {
                            if (c2[si] > SMALL)
                            {
                                kf *= el*pow(c2[si] + VSMALL, el - 1.0);
                            }
                            else
                            {
                                kf = 0.0;
                            }
                        }
                        else
                        {
                            kf *= el*pow(c2[si], el - 1.0);
                        }
                    }
                    else
                    {
                        kf *= pow(c2[si], el);
                    }
                }

                forAll(R.lhs(), i)
                {
                    label si = R.lhs()[i].index;
                    if (mechRed_->active())
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sl = R.lhs()[i].stoichCoeff;
                    dfdc[si][sj] -= sl*kf;
                }
                forAll(R.rhs(), i)
                {
                    label si = R.rhs()[i].index;
                    if (mechRed_->active())
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sr = R.rhs()[i].stoichCoeff;
                    dfdc[si][sj] += sr*kf;
                }
            }

            forAll(R.rhs(), j)
            {
                label sj = R.rhs()[j].index;
                if (mechRed_->active())
                {
                    sj = completeToSimplifiedIndex_[sj];
                }
                scalar kr = kr0;
                forAll(R.rhs(), i)
                {
                    const label si = R.rhs()[i].index;
                    const scalar er = R.rhs()[i].exponent;
                    if (i == j)
                    {
                        if (er < 1.0)
                        {
                            if (c2[si] > SMALL)
                            {
                                kr *= er*pow(c2[si] + VSMALL, er - 1.0);
                            }
                            else
                            {
                                kr = 0.0;
                            }
                        }
                        else
                        {
                            kr *= er*pow(c2[si], er - 1.0);
                        }
                    }
                    else
                    {
                        kr *= pow(c2[si], er);
                    }
                }

                forAll(R.lhs(), i)
                {
                    label si = R.lhs()[i].index;
                    if (mechRed_->active())
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sl = R.lhs()[i].stoichCoeff;
                    dfdc[si][sj] += sl*kr;
                }
                forAll(R.rhs(), i)
                {
                    label si = R.rhs()[i].index;
                    if (mechRed_->active())
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sr = R.rhs()[i].stoichCoeff;
                    dfdc[si][sj] -= sr*kr;
                }
            }
        }
    }

    // Calculate the dcdT elements numerically
    const scalar delta = 1.0e-3;
    const scalarField dcdT0(omega(c, T - delta, p));
    const scalarField dcdT1(omega(c, T + delta, p));

    for (label i = 0; i < this->nEqns(); i++)
    {
        dfdc[i][this->nSpecie()] = 0.5*(dcdT1[i] - dcdT0[i])/delta;
    }
}


template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::jacobian
(
 const scalar t,
 const scalarField& c,
 scalarRectangularMatrix& dfdc
 ) const
{
    //if the mechanism reduction is active, the computed Jacobian
    //is compact (size of the reduced set of species)
    //but according to the informations of the complete set
    //(i.e. for the third-body efficiencies)
    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    scalarField c2(completeC_.size(), 0.0);
    if (mechRed_->active())
    {
        c2 = completeC_;
        for(label i=0; i<NsDAC_; i++)
        {
            c2[simplifiedToCompleteIndex_[i]] = max(0.0, c[i]);
        }
    }
    else
    {
        forAll(c2, i)
        {
            c2[i] = max(c[i], 0.0);
        }
    }
    for (label i=0; i<this->nEqns(); i++)
    {
        for (label j=0; j<this->nEqns(); j++)
        {
            dfdc[i][j] = 0.0;
        }
    }

    forAll(this->reactions_, ri)
    {
        if (!reactionsDisabled_[ri])
        {
            const Reaction<ThermoType>& R = this->reactions_[ri];

            const scalar kf0 = R.kf(p, T, c2);
            const scalar kr0 = R.kr(kf0, p, T, c2);

            forAll(R.lhs(), j)
            {
                label sj = R.lhs()[j].index;
                if (mechRed_->active())
                {
                    sj = completeToSimplifiedIndex_[sj];
                }
                scalar kf = kf0;
                forAll(R.lhs(), i)
                {
                    const label si = R.lhs()[i].index;
                    const scalar el = R.lhs()[i].exponent;
                    if (i == j)
                    {
                        if (el < 1.0)
                        {
                            if (c2[si] > SMALL)
                            {
                                kf *= el*pow(c2[si] + VSMALL, el - 1.0);
                            }
                            else
                            {
                                kf = 0.0;
                            }
                        }
                        else
                        {
                            kf *= el*pow(c2[si], el - 1.0);
                        }
                    }
                    else
                    {
                        kf *= pow(c2[si], el);
                    }
                }

                forAll(R.lhs(), i)
                {
                    label si = R.lhs()[i].index;
                    if (mechRed_->active())
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sl = R.lhs()[i].stoichCoeff;
                    dfdc[si][sj] -= sl*kf;
                }
                forAll(R.rhs(), i)
                {
                    label si = R.rhs()[i].index;
                    if (mechRed_->active())
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sr = R.rhs()[i].stoichCoeff;
                    dfdc[si][sj] += sr*kf;
                }
            }

            forAll(R.rhs(), j)
            {
                label sj = R.rhs()[j].index;
                if (mechRed_->active())
                {
                    sj = completeToSimplifiedIndex_[sj];
                }
                scalar kr = kr0;
                forAll(R.rhs(), i)
                {
                    const label si = R.rhs()[i].index;
                    const scalar er = R.rhs()[i].exponent;
                    if (i == j)
                    {
                        if (er < 1.0)
                        {
                            if (c2[si] > SMALL)
                            {
                                kr *= er*pow(c2[si] + VSMALL, er - 1.0);
                            }
                            else
                            {
                                kr = 0.0;
                            }
                        }
                        else
                        {
                            kr *= er*pow(c2[si], er - 1.0);
                        }
                    }
                    else
                    {
                        kr *= pow(c2[si], er);
                    }
                }

                forAll(R.lhs(), i)
                {
                    label si = R.lhs()[i].index;
                    if (mechRed_->active())
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sl = R.lhs()[i].stoichCoeff;
                    dfdc[si][sj] += sl*kr;
                }
                forAll(R.rhs(), i)
                {
                    label si = R.rhs()[i].index;
                    if (mechRed_->active())
                    {
                        si = completeToSimplifiedIndex_[si];
                    }
                    const scalar sr = R.rhs()[i].stoichCoeff;
                    dfdc[si][sj] -= sr*kr;
                }
            }
        }
    }

    // Calculate the dcdT elements numerically
    const scalar delta = 1.0e-3;
    const scalarField dcdT0(omega(c, T - delta, p));
    const scalarField dcdT1(omega(c, T + delta, p));

    for (label i = 0; i < this->nEqns(); i++)
    {
        dfdc[i][this->nSpecie()] = 0.5*(dcdT1[i] - dcdT0[i])/delta;
    }
}


template<class CompType, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::TDACChemistryModel<CompType, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    CompType::correct();

    scalar deltaTMin = GREAT;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    const volScalarField rho
    (
        IOobject
        (
            "rho",
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->thermo().rho()
    );

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField c(this->nSpecie_);
    scalarField c0(this->nSpecie_);

    forAll(rho, celli)
    {
        const scalar rhoi = rho[celli];
        scalar pi = p[celli];
        scalar Ti = T[celli];

        scalarField phiq(this->nEqns());//composition vector (Yi, T, p)
        for (label i=0; i<this->nSpecie_; i++)
        {
            c[i] = rhoi*this->Y_[i][celli]/this->specieThermo_[i].W();
            c0[i] = c[i];
            phiq[i] = this->Y()[i][celli];
        }
        phiq[this->nSpecie()]=Ti;
        phiq[this->nSpecie()+1]=pi;

        // Initialise time progress
        scalar timeLeft = deltaT[celli];

        scalarField Rphiq(this->nEqns(),0.0);
        // When tabulation is active (short-circuit evaluation for retrieve)
        // It first tries to retrieve the solution of the system with the
        // information stored through the tabulation method
        if (tabulation_->active() && tabulation_->retrieve(phiq, Rphiq))
        {
            // Retrieved solution stored in Rphiq
            for (label i=0; i<this->nSpecie(); i++)
            {
                c[i] = rhoi*Rphiq[i]/this->specieThermo_[i].W();
            }
        }
        // This position is reached when tabulation is not used OR
        // if the solution is not retrieved.
        // In the latter case, it adds the information to the tabulation
        // (it will either expand the current data or add a new stored poin).
        else
        {
            if (mechRed_->active())
            {
                //reduce mechanism change the number of species (only active)
                mechRed_->reduceMechanism(c,Ti,pi);
            }
            // Calculate the chemical source terms
            while (timeLeft > SMALL)
            {
                scalar dt = timeLeft;
                if (mechRed_->active())
                {
                    //completeC_ used in the overridden ODE methods
                    //to update only the active species
                    completeC_ = c;
                    //solve the reduced set of ODE

                    this->solve
                    (
                        simplifiedC_, Ti, pi, dt, this->deltaTChem_[celli]
                    );
                    for (label i=0; i<NsDAC_; i++)
                    {
                        c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                    }
                }
                else
                {
                    this->solve(c, Ti, pi, dt, this->deltaTChem_[celli]);
                }
                timeLeft -= dt;
            }

            // If tabulation is used, we add the information computed here to
            // the stored points (either expand or add)
            if (tabulation_->active())
            {
                forAll(c,i)
                {
                    Rphiq[i] = c[i]/rhoi*this->specieThermo_[i].W();
                }
                Rphiq[Rphiq.size()-2] = Ti;
                Rphiq[Rphiq.size()-1] = pi;
                tabulation_->add(phiq, Rphiq, rhoi);
            }

            // When operations are done and if mechanism reduction is active,
            // the number of species (which also affects nEqns) is set back
            // to the total number of species (stored in the mechRed object)
            if (mechRed_->active())
            {
                this->nSpecie_ = mechRed_->nSpecie();
            }
            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);
        }

        // Set the RR vector (used in the solver)
        for (label i=0; i<this->nSpecie_; i++)
        {
            this->RR_[i][celli] =
                (c[i] - c0[i])*this->specieThermo_[i].W()/deltaT[celli];
        }
    }
    if (tabulation_->active())
    {
        //every time-step, look if the tabulation should be updated
        tabulation_->update();
        //write the performance of the tabulation
        tabulation_->writePerformance();
    }

    return deltaTMin;
}


template<class CompType, class ThermoType>
Foam::scalar Foam::TDACChemistryModel<CompType, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar> >(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class CompType, class ThermoType>
Foam::scalar Foam::TDACChemistryModel<CompType, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class CompType, class ThermoType>
void Foam::TDACChemistryModel<CompType, ThermoType>::solve
(
    scalarField &c,
    scalar& T,
    scalar& p,
    scalar& deltaT,
    scalar& subDeltaT
) const
{
    notImplemented
    (
        "TDACChemistryModel::solve"
        "("
            "scalarField&, "
            "scalar&, "
            "scalar&, "
            "scalar&, "
            "scalar&"
        ") const"
    );
}


// ************************************************************************* //

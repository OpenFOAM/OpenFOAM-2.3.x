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

#include "externalWallHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
const char*
NamedEnum
<externalWallHeatFluxTemperatureFvPatchScalarField::operationMode, 3>::names[]=
{
    "fixed_heat_flux",
    "fixed_heat_transfer_coefficient",
    "unknown"
};

const NamedEnum
<
    externalWallHeatFluxTemperatureFvPatchScalarField::operationMode, 3
>
externalWallHeatFluxTemperatureFvPatchScalarField::operationModeNames;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined-K"),
    mode_(unknown),
    q_(p.size(), 0.0),
    h_(p.size(), 0.0),
    Ta_(p.size(), 0.0),
    thicknessLayers_(),
    kappaLayers_()
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf.KMethod(), ptf.kappaName()),
    mode_(ptf.mode_),
    q_(ptf.q_, mapper),
    h_(ptf.h_, mapper),
    Ta_(ptf.Ta_, mapper),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_)
{}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mode_(unknown),
    q_(p.size(), 0.0),
    h_(p.size(), 0.0),
    Ta_(p.size(), 0.0),
    thicknessLayers_(),
    kappaLayers_()
{
    if (dict.found("q") && !dict.found("h") && !dict.found("Ta"))
    {
        mode_ = fixedHeatFlux;
        q_ = scalarField("q", dict, p.size());
    }
    else if (dict.found("h") && dict.found("Ta") && !dict.found("q"))
    {
        mode_ = fixedHeatTransferCoeff;
        h_ = scalarField("h", dict, p.size());
        Ta_ = scalarField("Ta", dict, p.size());
        if (dict.found("thicknessLayers"))
        {
            dict.lookup("thicknessLayers") >> thicknessLayers_;
            dict.lookup("kappaLayers") >> kappaLayers_;
        }
    }
    else
    {
        FatalErrorIn
        (
            "externalWallHeatFluxTemperatureFvPatchScalarField::"
            "externalWallHeatFluxTemperatureFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n patch type '" << p.type()
            << "' either q or h and Ta were not found '"
            << "\n for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    temperatureCoupledBase(tppsf),
    mode_(tppsf.mode_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_)
{}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    temperatureCoupledBase(patch(), tppsf.KMethod(), tppsf.kappaName()),
    mode_(tppsf.mode_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    q_.autoMap(m);
    h_.autoMap(m);
    Ta_.autoMap(m);
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const externalWallHeatFluxTemperatureFvPatchScalarField& tiptf =
        refCast<const externalWallHeatFluxTemperatureFvPatchScalarField>(ptf);

    q_.rmap(tiptf.q_, addr);
    h_.rmap(tiptf.h_, addr);
    Ta_.rmap(tiptf.Ta_, addr);
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField q(size(), 0.0);
    const scalarField Tc(patchInternalField());
    const scalarField Tp(*this);
    const scalarField KWall(kappa(Tp));
    const scalarField KDelta(KWall*patch().deltaCoeffs());

    switch (mode_)
    {
        case fixedHeatFlux:
        {
            q = q_;
            break;
        }
        case fixedHeatTransferCoeff:
        {
            scalar totalSolidRes = 0.0;
            if (thicknessLayers_.size() > 0)
            {
                forAll (thicknessLayers_, iLayer)
                {
                    const scalar l = thicknessLayers_[iLayer];
                    if (kappaLayers_[iLayer] > 0.0)
                    {
                        totalSolidRes += l/kappaLayers_[iLayer];
                    }
                }
            }
            q = (Ta_ - Tp)*(1.0/h_ + totalSolidRes);
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "externalWallHeatFluxTemperatureFvPatchScalarField"
                "::updateCoeffs()"
            )   << "Illegal heat flux mode " << operationModeNames[mode_]
                << exit(FatalError);
        }
    }

    forAll(*this, i)
    {
        if (q[i] > 0) //in
        {
            this->refGrad()[i] = q[i]/KWall[i];
            this->refValue()[i] = 0.0;
            this->valueFraction()[i] = 0.0;
        }
        else //out
        {
            this->refGrad()[i] = 0.0;
            this->refValue()[i] = q[i]/KDelta[i] + Tc[i];
            this->valueFraction()[i] = 1.0;
        }
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(KWall*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    temperatureCoupledBase::write(os);

    switch (mode_)
    {
        case fixedHeatFlux:
        {
            q_.writeEntry("q", os);
            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.writeEntry("h", os);
            Ta_.writeEntry("Ta", os);
            os.writeKeyword("thicknessLayers")<< thicknessLayers_
                << token::END_STATEMENT << nl;
            os.writeKeyword("kappaLayers")<< kappaLayers_
                << token::END_STATEMENT << nl;
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "void externalWallHeatFluxTemperatureFvPatchScalarField::write"
                "("
                    "Ostream& os"
                ") const"
            )   << "Illegal heat flux mode " << operationModeNames[mode_]
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalWallHeatFluxTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //

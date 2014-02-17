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

#include "temperatureCoupledBase.H"
#include "volFields.H"
#include "fluidThermo.H"
#include "solidThermo.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::temperatureCoupledBase::KMethodType,
        4
    >::names[] =
    {
        "fluidThermo",
        "solidThermo",
        "directionalSolidThermo",
        "lookup"
    };
}


const Foam::NamedEnum<Foam::temperatureCoupledBase::KMethodType, 4>
    Foam::temperatureCoupledBase::KMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const word& kappaName
)
:
    patch_(patch),
    method_(KMethodTypeNames_[calculationType]),
    kappaName_(kappaName)
{}


Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(KMethodTypeNames_.read(dict.lookup("kappa"))),
    kappaName_(dict.lookup("kappaName"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::temperatureCoupledBase::kappa
(
    const scalarField& Tp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    switch (method_)
    {
        case mtFluidThermo:
        {
            typedef compressible::turbulenceModel turbulenceModel;

            if (mesh.foundObject<turbulenceModel>("turbulenceModel"))
            {
                const turbulenceModel& turbModel =
                    mesh.lookupObject<turbulenceModel>("turbulenceModel");

                return turbModel.kappaEff(patch_.index());
            }
            else if (mesh.foundObject<fluidThermo>("thermophysicalProperties"))
            {
                const fluidThermo& thermo =
                    mesh.lookupObject<fluidThermo>("thermophysicalProperties");

                return thermo.kappa(patch_.index());
            }
            else
            {
                FatalErrorIn
                (
                    "temperatureCoupledBase::kappa(const scalarField&) const"
                )
                    << "Kappa defined to employ " << KMethodTypeNames_[method_]
                    << " method, but thermo package not available"
                    << exit(FatalError);
            }

            break;
        }

        case mtSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>("thermophysicalProperties");

            return thermo.kappa(patch_.index());
            break;
        }

        case mtDirectionalSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>("thermophysicalProperties");

            const vectorField kappa(thermo.Kappa(patch_.index()));

            tmp<scalarField> tmeanKappa(Tp);
            scalarField& meanKappa = tmeanKappa();
            forAll(meanKappa, i)
            {
                meanKappa[i] = (kappa[i].x() + kappa[i].y() + kappa[i].z())/3.0;
            }

            return meanKappa;
            break;
        }

        case mtLookup:
        {
            if (mesh.foundObject<volScalarField>(kappaName_))
            {
                return patch_.lookupPatchField<volScalarField, scalar>
                (
                    kappaName_
                );
            }
            else if (mesh.foundObject<volSymmTensorField>(kappaName_))
            {
                const symmTensorField& KWall =
                    patch_.lookupPatchField<volSymmTensorField, scalar>
                    (
                        kappaName_
                    );

                const vectorField n(patch_.nf());

                return n & KWall & n;
            }
            else
            {
                FatalErrorIn
                (
                    "temperatureCoupledBase::kappa(const scalarField&) const"
                )
                    << "Did not find field " << kappaName_
                    << " on mesh " << mesh.name() << " patch " << patch_.name()
                    << nl
                    << "Please set 'kappa' to one of "
                    << KMethodTypeNames_.toc()
                    << " and 'kappaName' to the name of the volScalar"
                    << " or volSymmTensor field (if kappa=lookup)"
                    << exit(FatalError);
            }

            break;
        }

        default:
        {
            FatalErrorIn
            (
                "temperatureCoupledBase::kappa(const scalarField&) const"
            )
                << "Unimplemented method " << KMethodTypeNames_[method_] << nl
                << "Please set 'kappa' to one of " << KMethodTypeNames_.toc()
                << " and 'kappaName' to the name of the volScalar"
                << " or volSymmTensor field (if kappa=lookup)"
                << exit(FatalError);
        }
    }

    return scalarField(0);
}


void Foam::temperatureCoupledBase::write(Ostream& os) const
{
    os.writeKeyword("kappa") << KMethodTypeNames_[method_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("kappaName") << kappaName_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //

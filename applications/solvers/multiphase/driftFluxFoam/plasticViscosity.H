volScalarField plasticViscosity
(
    const dimensionedScalar& plasticViscosityCoeff,
    const dimensionedScalar& plasticViscosityExponent,
    const volScalarField& Alpha
)
{
    tmp<volScalarField> tfld
    (
         plasticViscosityCoeff*
         (
             pow
             (
                 10.0,
                 plasticViscosityExponent*Alpha + SMALL
             ) - scalar(1)
         )
    );

    return tfld();
}

volScalarField plasticViscosity
(
    const dimensionedScalar& plasticViscosityCoeff,
    const dimensionedScalar& plasticViscosityExponent,
    const volScalarField& alpha
)
{
    tmp<volScalarField> tfld
    (
         plasticViscosityCoeff*
         (
             pow(10.0, plasticViscosityExponent*alpha + SMALL) - scalar(1)
         )
    );

    return tfld();
}

volScalarField yieldStress
(
    const dimensionedScalar& yieldStressCoeff,
    const dimensionedScalar& yieldStressExponent,
    const dimensionedScalar& yieldStressOffset,
    const volScalarField& alpha
)
{
    tmp<volScalarField> tfld
    (
        yieldStressCoeff*
        (
            pow(10.0, yieldStressExponent*(alpha + yieldStressOffset))
          - pow(10.0, yieldStressExponent*yieldStressOffset)
        )
    );

    return tfld();
}

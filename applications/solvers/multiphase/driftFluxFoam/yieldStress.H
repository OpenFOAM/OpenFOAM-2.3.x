volScalarField yieldStress
(
    const dimensionedScalar& yieldStressCoeff,
    const dimensionedScalar& yieldStressExponent,
    const dimensionedScalar& yieldStressOffset,
    const volScalarField& alpha1
)
{
    tmp<volScalarField> tfld
    (
        yieldStressCoeff*
        (
            pow
            (
                10.0,
                yieldStressExponent*(max(alpha1, scalar(0)) + yieldStressOffset)
            )
          - pow
            (
                10.0,
                yieldStressExponent*yieldStressOffset
            )
        )
    );

    return tfld();
}

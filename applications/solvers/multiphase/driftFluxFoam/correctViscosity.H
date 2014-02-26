{
    mul = muc +
        plasticViscosity
        (
            plasticViscosityCoeff,
            plasticViscosityExponent,
            alpha1
        );

    if (BinghamPlastic)
    {
        volScalarField tauy = yieldStress
        (
            yieldStressCoeff,
            yieldStressExponent,
            yieldStressOffset,
            alpha1
        );

        mul =
            tauy/
            (
                mag(fvc::grad(U))
              + 1.0e-4*
                (
                    tauy
                  + dimensionedScalar
                    (
                        "deltaTauy",
                        tauy.dimensions(),
                        1.0e-15
                    )
                )/mul
            )
          + mul;
    }

    mul = min(mul, muMax);
}

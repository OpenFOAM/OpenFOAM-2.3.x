{
    const scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    const scalar kappa_ = kappa.value();
    const scalar E_ = E.value();

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            scalarField& mutw = mut.boundaryField()[patchi];
            const scalarField& mulw = mul.boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                scalar yPlus =
                    Cmu25*y[patchi][facei]*::sqrt(k[faceCelli])
                   /(mulw[facei]/rho2.value());

                if (yPlus > 11.6)
                {
                    mutw[facei] =
                        mulw[facei]*(yPlus*kappa_/::log(E_*yPlus) - 1);
                }
                else
                {
                    mutw[facei] = 0.0;
                }
            }
        }
    }
}

{
    const scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    const scalar kappa_ = kappa.value();
    const scalar E_ = E.value();
    const scalar muc_ = muc.value();
    const scalar nuc_ = muc_/rhoc.value();

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            scalarField& mutw = mut.boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                scalar yPlus =
                    Cmu25*y[patchi][facei]*::sqrt(k[faceCelli])
                   /nuc_;

                if (yPlus > 11.6)
                {
                    mutw[facei] =
                        muc_*(yPlus*kappa_/::log(E_*yPlus) - 1);
                }
                else
                {
                    mutw[facei] = 0.0;
                }
            }
        }
    }
}

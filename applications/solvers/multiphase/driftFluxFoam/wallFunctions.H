{
    labelList cellBoundaryFaceCount(epsilon.size(), 0);

    const scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    const scalar Cmu75 = ::pow(Cmu.value(), 0.75);
    const scalar kappa_ = kappa.value();
    const scalar muc_ = muc.value();

    const fvPatchList& patches = mesh.boundary();

    //- Initialise the near-wall P field to zero
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon[faceCelli] = 0.0;
                G[faceCelli] = 0.0;
            }
        }
    }

    //- Accumulate the wall face contributions to epsilon and G
    //  Increment cellBoundaryFaceCount for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            const scalarField& mutw = mut.boundaryField()[patchi];

            scalarField magFaceGradU
            (
                mag(U.boundaryField()[patchi].snGrad())
            );

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                // For corner cells (with two boundary or more faces),
                // epsilon and G in the near-wall cell are calculated
                // as an average

                cellBoundaryFaceCount[faceCelli]++;

                epsilon[faceCelli] +=
                    Cmu75*::pow(k[faceCelli], 1.5)
                   /(kappa_*y[patchi][facei]);

                G[faceCelli] +=
                    (mutw[facei] + muc_)
                   *magFaceGradU[facei]
                   *Cmu25*::sqrt(k[faceCelli])
                   /(kappa_*y[patchi][facei]);
            }
        }
    }


    // perform the averaging

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                G[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                cellBoundaryFaceCount[faceCelli] = 1;
            }
        }
    }
}

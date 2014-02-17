{
    const scalarField& V = mesh.V();

    forAll(mesh.cellZones(), czi)
    {
        const labelList& cellLabels = mesh.cellZones()[czi];

        scalar phaseVolume = 0;

        forAll(cellLabels, cli)
        {
            label celli = cellLabels[cli];
            phaseVolume += alpha1[celli]*V[celli];
        }

        reduce(phaseVolume, sumOp<scalar>());

        Info<< "Phase volume in zone " << mesh.cellZones()[czi].name()
            << " = " << phaseVolume*1e6 << " ml " << endl;
    }
}

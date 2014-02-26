if (VdjModel == "general")
{
    Vdj = V0*
    (
        exp(-a*max(alpha1 - alphaMin, scalar(0)))
      - exp(-a1*max(alpha1 - alphaMin, scalar(0)))
    );
}
else if (VdjModel == "simple")
{
    Vdj = V0*pow(10.0, -a*max(alpha1, scalar(0)));
}
else
{
    FatalErrorIn(args.executable())
        << "Unknown VdjModel : " << VdjModel
        << abort(FatalError);
}

Vdj.correctBoundaryConditions();

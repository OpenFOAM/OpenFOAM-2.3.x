if (VdjModel == "general")
{
    Vdj = V0*
    (
        exp(-a*max(alpha - alphaMin, scalar(0)))
      - exp(-a1*max(alpha - alphaMin, scalar(0)))
    );
}
else if (VdjModel == "simple")
{
    Vdj = V0*pow(10.0, -a*alpha);
}
else
{
    FatalErrorIn(args.executable())
        << "Unknown VdjModel : " << VdjModel
        << abort(FatalError);
}

Vdj.correctBoundaryConditions();

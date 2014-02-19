const dictionary& alphaControls = mesh.solverDict(alpha.name());

label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));
label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));

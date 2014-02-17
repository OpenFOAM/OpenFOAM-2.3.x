fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
);

solve(UEqn == -fvc::grad(p));
K = 0.5*magSqr(U);

//The simpliest way is not to take strong dependence of the diffusion koefficient ks on the concentration S into account
//ks=kd;
//The better way is to count it

keff = kd*pow((1.0-S),6.5)*(1.0+2.0*S*(4.0-S)/pow((1.0-S),4)-8.0*pow(lambda,2)*S/(3.0*pow((1.0+1.25*lambda*S),2)));

fvScalarMatrix SEqn
    (
      fvm::ddt(S)
     +fvm::div(phi,S)
     -fvm::laplacian(keff,S)
    );
    SEqn.solve();

rho1 = rhoker + (rhomag - rhoker)*S;
rho2 = rhoker + (rhomag - rhoker)*S;
rho = alpha1*rho1 + alpha2*rho2;

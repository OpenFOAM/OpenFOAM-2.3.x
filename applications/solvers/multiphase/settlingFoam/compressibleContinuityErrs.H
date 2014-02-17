    scalar sumLocalContErr =
        runTime.deltaTValue()*
        mag
        (
            fvc::ddt(rho)
          + fvc::div(phi)
        )().weightedAverage(rho*mesh.V()).value();

    scalar globalContErr =
          runTime.deltaTValue()*
          (
              fvc::ddt(rho)
            + fvc::div(phi)
          )().weightedAverage(rho*mesh.V()).value();

    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr
        << endl;

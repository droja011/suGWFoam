    dictionary simulation = mesh.solutionDict().subDict("SIMULATION");

    //maximum nonlinear iteration steps
    int nNonLinear =
        simulation.lookupOrDefault<int>("nNonlinear", 100);


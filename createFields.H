    Info<< "Reading pressure head field h\n" << endl;
    volScalarField h
    (
        IOobject
        (
            "h",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //whether to do project at each step to help convergence
    Switch Projection
    (
        transportProperties.lookup("Projection")
    );

    //whether to switch between different forms: h-based and mixed
    Switch formSwitch
    (
        transportProperties.lookup("formSwitch")
    );

    //only when formSwitch is turned off: which form to use
    word equationForm
    (
        transportProperties.lookup("equationForm")
    );

    //only when formSwitch is turned on: the threshold head change to switch between 
    //h-based or mixed form
    dimensionedScalar hswitch
    (
        transportProperties.lookup("hswitch")
    );


    //Source/sink field, such as pumping and injection
    Info<< "Creating source/sink field h_ss\n" << endl;

    volScalarField h_ss
    (
        IOobject
        (
            "h_ss",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
            (
              "0.0", 
              dimensionSet(0,0,-1,0,0,0,0), 
              0.0 
            ),
        zeroGradientFvPatchScalarField::typeName
    );


    //moisture content theta
    //Initial value does not matter as it will be updated.
    volScalarField theta
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
             "0.0",
             dimensionSet(0,0,0,0,0,0,0),
             0.0
        )
    );

    //relative permeability
    volScalarField kr
    (
        IOobject
        (
            "kr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
             "0.0",
             dimensionSet(0,0,0,0,0,0,0),
             0.0
        )
    );


    //specific moisture capacity
    volScalarField Ch
    (
        IOobject
        (
            "Ch",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
             "0.0",
             dimensionSet(0,-1,0,0,0,0,0),
             0.0
        )
    );

    //cell center coordinate
    volScalarField z
    (
        IOobject
        (
            "z",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C().component(vector::Z)
    );

    //gradient of z
    volVectorField grad_z
    (
        IOobject
        (
            "grad_z",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(z)
    );


    //some aux fields

    //h_old is the value at the beginning of each time step
    volScalarField h_old
    (
           IOobject
           (
               "h_old",
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           h
    );

    //h_old_tmp is the value at the beginning of each time step
    volScalarField h_old_tmp
    (
           IOobject
           (
               "h_old_tmp",
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::NO_WRITE
           ),
           h
    );


    //theta_old is the value at the beginning of each time step
    volScalarField theta_old
    (
           IOobject
           (
               "theta_old",
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           theta
    );

    //h_n is the value at the previous Picard iteration
    volScalarField h_n
    (
           IOobject
           (
               "h_n",
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           h
    );
//        h_n.boundaryField().updateCoeffs();


    //h_diff is the change between two Picard iterations
    volScalarField h_diff
    (
           IOobject
           (
               "h_diff",
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           h-h_n
    );


    //create the soil model
    autoPtr<soilModel> soil
    (
       soilModel::New("soil", transportProperties, h, theta, kr, Ch)
    );

    //update some of the coefficients.
    soil->correct();

    //theta_n is the value at the previous Picard iteration
    volScalarField theta_n
    (
           IOobject
           (
               "theta_n",
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
           ),
           theta
    );

    //interstitial velocity 
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -(soil->Ks()&(fvc::grad(h+z)))*kr/soil->theta_s()
    );

#   include "createPhi.H"


    //log file
    mkDir(runTime.path()/"simuLog");
    OFstream simuLogFile(runTime.path()/"simuLog"/"simuLog.dat");

    scalar totalFlux_through_boundary = 0.0;
    scalar totalMassChange_in_domain = 0.0;
    scalar totalWaterCompress_in_domain = 0.0;

    //Decide which convergence criterion to use:
    //create the convergence test model
    autoPtr<convergenceCriterion> converge
    (
       convergenceCriterion::New("converge", transportProperties, h, h_n, theta, theta_n)
    );

    //temporary field for diagonostic purpose
    volScalarField thetatt
    (
        IOobject
        (
            "thetatt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
             "0.0",
             dimensionSet(0,0,0,0,0,0,0),
             0.0
        )
    );


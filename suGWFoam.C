/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    suGWFoam (i.e., saturated-unsaturated GroundWater Foam)

Description
    Solves the nonlinear Richards equation for groundwater flow. 

Reference: 
X. Liu (2012). "suGWFoam: An Open Source Saturated-Unsaturated GroundWater Flow
  Solver based on OpenFOAM", Civil and Environmental Engineering Studies
  Technical Report No. 01-01, University of Texas at San Antonio, San Antonio,
  Texas, U.S.A.

Author
    Xiaofeng Liu, Ph.D., P.E.
    Assistant Professor
    Department of Civil and Environmental Engineering
    University of Texas at San Antonio
    email: xiaofeng.liu@utsa.edu
    web:http://engineering.utsa.edu/~xiaofengliu

    November, 2012

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "soilModel.H"
#include "convergenceCriterion.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "readTimeParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nSimulating unsaturated groudwater flow\n" << endl;

    //nonlinear iteration counter
    int nonLinear=0; 
    
    //equation form indicator: 1: h-based, 2:mixed form
    int nForm=0;

    //set a fake initial nonLinear iteration number so it can start
//    nonLinear = (int) (m_lower+m_upper)/2;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readSimulationControls.H"
#       include "readTimeParameters.H"
#       include "adjustDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        Info<< "deltaT = " <<  runTime.deltaT().value() << endl;

        //The following is for the choice of equation form
        if(formSwitch)   //if formSwitch is turn on
        {
          //Switch forms: each time step only choose one form based on the change 
          //              of head between h_n and h_(n-1)
          if(max(mag((h.internalField()-h_old.internalField())))<hswitch.value())
          {
                Info << "Selecting h-based form, code: 1" << endl;
                nForm = 1;
          }
          else
          {
                Info << "Selecting mixed form, code: 2" << endl;
                nForm = 2;
          }
        }
        else           //if formSwitch is turned off
        {
          if(equationForm.match("h-based"))      //h-based form
          {
                Info << "Selecting h-based form, code: 1" << endl;
                nForm = 1;
          }
          else if(equationForm.match("mixed")) //mixed form
          {
                Info << "Selecting mixed form, code: 2" << endl;
                nForm = 2;
          }
          else                     //wrong choice
          {
                FatalErrorIn("equationForm options")
                  << "The option for equationForm in transportDict: "
                  << " is not valid." << nl
                  << "Valid options are: 1 and 2"
                  << abort(FatalError);
          }
        }
        //End of choice of equation form
        
        //Save a copy of h before the projection step
        h_old_tmp = h;
        h_old_tmp.boundaryField().updateCoeffs();


        //Whether do a linear projection at the beginning?
        //In some cases, it might help the convergence.
        if(Projection)
        {
             Info << "Doing projection at the beginning of time step!" << endl;
             h = 2.0*h - h_old; 
             h.boundaryField().updateCoeffs();
        }

        //restore the old time h value
        h_old = h_old_tmp;
        h_old.boundaryField().updateCoeffs();

        //Save a copy of theta
        theta_old = theta;

        //At the beginning of Picard iteration, use the current value of h
        //to start with.
        h_n = h;

        //Picard iteration
        for (nonLinear=0; nonLinear<=nNonLinear; nonLinear++)
        {
            Info << "nonlinear iteration # = " << nonLinear << endl;

            //updated the soil parameters
            soil->correct();
 
            volScalarField coeff_A1 = Ch+soil->Sw()*soil->Ss();
            volScalarField coeff_A2 = Ch;

            //prepare the dispersion coefficient K=Ks*kr
            volTensorField coeff_B = soil->Ks()*kr;            

            volVectorField gravity = coeff_B & grad_z;

            surfaceScalarField phiG
            ( 
               "phiG",
               fvc::interpolate(coeff_B & grad_z) & mesh.Sf()
            );

            //switch forms 
            if(nForm==1)  //h-based form
            {
               solve
               (
                   fvm::ddt(coeff_A1, h) 
                 - fvm::laplacian(coeff_B, h)
                 - fvc::div(gravity)
               );
            }
            else if(nForm==2) //mixed form
            {
               dimensionedScalar rDeltaT = 1.0/runTime.deltaT();
               h_ss = rDeltaT*(theta_old - theta_n)+
                      rDeltaT*soil->Ss()*soil->Sw()*(h_old-h_n);

               //Trick OF to use the value of h at previous iteration as 
               //the values at previous time step. At least it works
               //for Euler scheme. Need to test other schemes where
               //multiple previous steps are used.
               //In "mixed" form, the "old" time step value should be 
               //the value at the previous iteration. See notes.
               h.oldTime() = h_n;

               solve
               (
                   fvm::ddt(coeff_A2, h)
                ==
                   fvm::laplacian(coeff_B, h)
                 + fvc::div(gravity)
                 + h_ss
               );

               //restore h.oldTime()
               h.oldTime() = h_old;
            }
            else //wrong form choice
            {
                FatalErrorIn(args.executable())
                        << "Wrong choice of form!"
                        << abort(FatalError);
            }

            h_diff = h - h_n;
            h_diff.boundaryField().updateCoeffs();
 
            //update moisture content theta for convergent test
            soil->update_theta();

            //convergence test
            if(converge->convergent() && 
               nonLinear>=2)  //force at least 1 nonlinear iteration
            {
                     h_n = h;
                     h_n.boundaryField().updateCoeffs();

                     theta_n = theta;

                     break;
            }
            else if(nonLinear==nNonLinear) //nonlinear iteration reached maxium
            {
                      h_n = h;
                      h_n.boundaryField().updateCoeffs();

                      theta_n = theta;

                      //Either
//                    FatalErrorIn(args.executable())
//                        << "Nonlinear iteration didn't converge !" 
//                        << abort(FatalError);
                      //Or
                      Info
                        << "Nonlinear iteration didn't converge !"
                        << endl;

            }
            else   //not converged nor reached the maximum iteration yet, continue
            {
                      h_n = h;
                      h_n.boundaryField().updateCoeffs();

                      theta_n = theta;
            }
        }

        Info << "Total nonlinear iterations: " << nonLinear << endl;
        //End of Picard iteration

        //update velocity and flux
        U = -(soil->Ks()&(fvc::grad(h+z)))*kr/soil->theta_s();
        phi = linearInterpolate(U) & mesh.Sf();

        //calculate mass balance at each time step
#       include "calcMassBalance.H"

        //write result
#       include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

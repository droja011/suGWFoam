   //output mass balance infomation       
   //current step
   Info << "\nMass balance (current step): " << endl;
   Info << "       volume change         = " << volMassChange << endl;
   Info << "       water compress        = " << volWaterCompress << endl;
   Info << "       flux through boundary = " << fluxBoundary  << endl;
   Info << "       absolute error        = " << (volMassChange+volWaterCompress+fluxBoundary) << endl;

   //cumulative
   Info << "\nMass balance (cumulative): " << endl;
   Info << "       volume change         = " << totalMassChange_in_domain << endl;
   Info << "       water compress        = " << totalWaterCompress_in_domain << endl;
   Info << "       flux through boundary = " << totalFlux_through_boundary  << endl;
   Info << "       absolute error        = "
        << (totalMassChange_in_domain+totalWaterCompress_in_domain+totalFlux_through_boundary) << endl;
   Info << "       relative error (\%)   = "
        << 100*(1.0-mag(totalFlux_through_boundary/(totalMassChange_in_domain+totalWaterCompress_in_domain+SMALL)))
        << endl;
   Info << endl;

   if (runTime.outputTime())
   {
        //write log file
#       include "output_log.H"

        runTime.write();
   }

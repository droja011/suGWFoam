        //output some information to the log file at current time
        simuLogFile << runTime.value() << tab;
        simuLogFile << nForm << tab;
        simuLogFile << nonLinear+1 << tab;
        simuLogFile << volMassChange << tab
                    << volWaterCompress << tab
                    << fluxBoundary  << tab
                    << volMassChange+volWaterCompress+fluxBoundary << tab;
        simuLogFile << totalMassChange_in_domain << tab
                    << totalWaterCompress_in_domain << tab
                    << totalFlux_through_boundary << tab
                    << totalFlux_through_boundary+totalMassChange_in_domain << tab
                    << 100*(1.0-mag(totalFlux_through_boundary/
                       (totalMassChange_in_domain+totalWaterCompress_in_domain+SMALL)))
                    << tab;
        simuLogFile << endl;




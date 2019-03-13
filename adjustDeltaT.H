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

Global
    adjustDeltaT

Description
    Reset the timestep to maintain a moderate nonlinear iteration number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    Info << "maxDeltaT = " << maxDeltaT << endl;
    Info << "minDeltaT = " << minDeltaT << endl;
    Info << "deltaT = " << runTime.deltaT().value() << endl;
    Info << "nonLinear,m_lower,m_upper = " << nonLinear << " " 
                                           << m_lower << " " 
                                           << m_upper << endl;

    if(nonLinear<m_lower) //time step is increased (slowly)
    {
      Info << "nonLinear < m_lower " 
           <<  min
               (
                 runTime.deltaT().value()*increase_factor,
                 maxDeltaT
               )
           << endl;

      runTime.setDeltaT
      (
          min
          (
              runTime.deltaT().value()*increase_factor,
              maxDeltaT
          )
      );
    }
    else if(nonLinear>m_upper) //time step is decreased
    {
      Info << "nonLinear > m_upper " 
           << max
              (
                  runTime.deltaT().value()*decrease_factor,
                  minDeltaT
              )
           << endl;
         

      runTime.setDeltaT
      (
          max
          (
              runTime.deltaT().value()*decrease_factor,
              minDeltaT
          )
      );
    }
    
    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}

// ************************************************************************* //

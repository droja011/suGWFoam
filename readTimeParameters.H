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
    readTimeParameters

Description
    Read the control parameters used by ajustDeltaT

\*---------------------------------------------------------------------------*/

Switch adjustTimeStep
(
    runTime.controlDict().lookup("adjustTimeStep")
);

//lower iteration number limit
label m_lower
(
    readLabel(runTime.controlDict().lookup("m_lower"))
);

//upper iteration number limit
label m_upper
(
    readLabel(runTime.controlDict().lookup("m_upper"))
);

//check whether m_lower < m_upper?
if(m_lower >= m_upper)
{
   Info << "m_lower is larger than m_upper!" << endl;
   exit(1);
}

//acceleration factor (>=1, but not too large, say 1.1)
scalar increase_factor
(
    readScalar(runTime.controlDict().lookup("increase_factor"))
);

//deceleration factor (<=1, say 0.5)
scalar decrease_factor
(
    readScalar(runTime.controlDict().lookup("decrease_factor"))
);

scalar minDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", SMALL);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

// ************************************************************************* //

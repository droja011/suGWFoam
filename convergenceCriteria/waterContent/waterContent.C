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

\*---------------------------------------------------------------------------*/

#include "waterContent.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace convergenceCriteria
{
    defineTypeNameAndDebug(waterContent, 0);

    addToRunTimeSelectionTable
    (
        convergenceCriterion,
        waterContent,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::convergenceCriteria::waterContent::waterContent
(
    const word& name,
    const dictionary& convergenceProperties,
    volScalarField& h,
    volScalarField& h_n,
    volScalarField& theta,
    volScalarField& theta_n
)
:
    convergenceCriterion(name, convergenceProperties, h, h_n, theta, theta_n),
    waterContentCoeffs_(convergenceProperties.subDict(typeName + "Coeffs"))
{
    delta_theta_ = readScalar(waterContentCoeffs_.lookup("delta_theta"));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::convergenceCriteria::waterContent::read
(
    const dictionary& convergenceProperties
)
{
    convergenceCriterion::read(convergenceProperties);

    waterContentCoeffs_ = convergenceProperties.subDict(typeName + "Coeffs");

    delta_theta_ = readScalar(waterContentCoeffs_.lookup("delta_theta"));

    return true;
}

// check whether convergent 
bool Foam::convergenceCriteria::waterContent::convergent()
{
   scalar maxChange = max(mag(theta_.internalField()-theta_n_.internalField()));
  
   Info << "Max change in theta: " << maxChange << endl;

   if(maxChange < delta_theta_)
   {
        return true;
   }
   else
   {
        return false;
   }
}

//*********************************************************** //

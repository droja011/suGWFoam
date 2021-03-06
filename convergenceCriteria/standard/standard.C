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

#include "standard.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace convergenceCriteria
{
    defineTypeNameAndDebug(standard, 0);

    addToRunTimeSelectionTable
    (
        convergenceCriterion,
        standard,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::convergenceCriteria::standard::standard
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
    standardCoeffs_(convergenceProperties.subDict(typeName + "Coeffs"))
{
    delta_h_ = readScalar(standardCoeffs_.lookup("delta_h"));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::convergenceCriteria::standard::read
(
    const dictionary& convergenceProperties
)
{
    convergenceCriterion::read(convergenceProperties);

    standardCoeffs_ = convergenceProperties.subDict(typeName + "Coeffs");

    delta_h_ = readScalar(standardCoeffs_.lookup("delta_h"));

    return true;
}

// check whether convergent 
bool Foam::convergenceCriteria::standard::convergent()
{
//   scalar maxChange = max(mag(h_.internalField()-h_n_.internalField()));
   scalar maxChange = max(mag(h_-h_n_)).value();

   Info << "Max change in h: " << maxChange << endl;

   if(maxChange < delta_h_)
   {
        return true;
   }
   else
   {
        return false;
   }
}

//*********************************************************** //

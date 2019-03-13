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

#include "Haverkamp.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace soilModels
{
    defineTypeNameAndDebug(Haverkamp, 0);

    addToRunTimeSelectionTable
    (
        soilModel,
        Haverkamp,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soilModels::Haverkamp::Haverkamp
(
    const word& name,
    const dictionary& soilProperties,
    volScalarField& h,
    volScalarField& theta,
    volScalarField& kr,
    volScalarField& Ch
)
:
    soilModel(name, soilProperties, h, theta, kr, Ch),
    HaverkampCoeffs_(soilProperties.subDict(typeName + "Coeffs")),
    Ks_(HaverkampCoeffs_.lookup("Ks")),
    theta_s_(HaverkampCoeffs_.lookup("theta_s")),
    theta_r_(HaverkampCoeffs_.lookup("theta_r")),
    alpha_(HaverkampCoeffs_.lookup("alpha")),
    beta_(HaverkampCoeffs_.lookup("beta")),
    gamma_(HaverkampCoeffs_.lookup("gamma")),
    A_(HaverkampCoeffs_.lookup("A")),
    Ss_(HaverkampCoeffs_.lookup("Ss")),
    Sw_
    (
        IOobject
        (
            "Sw",
            theta_.time().timeName(),
            theta_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        theta_/theta_s_
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::soilModels::Haverkamp::read
(
    const dictionary& soilProperties
)
{
    soilModel::read(soilProperties);

    HaverkampCoeffs_ = soilProperties.subDict(typeName + "Coeffs");

    HaverkampCoeffs_.lookup("Ks") >> Ks_;
    HaverkampCoeffs_.lookup("theta_s") >> theta_s_;
    HaverkampCoeffs_.lookup("theta_r") >> theta_r_;
    HaverkampCoeffs_.lookup("alpha") >> alpha_;
    HaverkampCoeffs_.lookup("beta") >> beta_;
    HaverkampCoeffs_.lookup("gamma") >> gamma_;
    HaverkampCoeffs_.lookup("A") >> A_;
    HaverkampCoeffs_.lookup("Ss") >> Ss_;

    return true;
}

// update soil moisture
void Foam::soilModels::Haverkamp::update_theta()
{  
   //update moisture content theta
   theta_ == neg(h_)*(alpha_*(theta_s_-theta_r_)/(alpha_+pow(dimensionedScalar("one", dimless/dimLength, 1.0)*mag(h_),beta_))+theta_r_)
            +pos(h_)*theta_s_;
}

// update the soil properties
void Foam::soilModels::Haverkamp::correct()
{
   //update the soil moisture
   update_theta();

   //update saturation ratio Sw
   Sw_ == theta_/theta_s_;

   //update relative permeability kr
   kr_ == neg(h_)*A_/(A_+pow(dimensionedScalar("one", dimless/dimLength, 1.0)*mag(h_),gamma_))
         +pos(h_)*1.0;

   //update specific moisture capacity (including the specific storage term)
   dimensionedScalar unit_dthetadh
   (
       "0.0",
       dimensionSet(0,0,-1,0,0,0,0),
       1.0
   );

   Ch_ == alpha_*beta_*(theta_s_-theta_r_)*pow(-h_*dimensionedScalar("one", dimless/dimLength, 1.0),beta_-1)/pow(alpha_+pow(-h_*dimensionedScalar("one", dimless/dimLength, 1.0),beta_),2.0)*unit_dthetadh
         +Ss_*Sw_;
}

//*********************************************************** //

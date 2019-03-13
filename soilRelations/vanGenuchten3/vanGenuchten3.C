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

#include "vanGenuchten3.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace soilModels
{
    defineTypeNameAndDebug(vanGenuchten3, 0);

    addToRunTimeSelectionTable
    (
        soilModel,
        vanGenuchten3,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soilModels::vanGenuchten3::vanGenuchten3
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
    vanGenuchten3Coeffs_(soilProperties.subDict(typeName + "Coeffs")),
    Ks_(vanGenuchten3Coeffs_.lookup("Ks")),
    theta_s_(vanGenuchten3Coeffs_.lookup("theta_s")),
    theta_r_(vanGenuchten3Coeffs_.lookup("theta_r")),
    alpha_(vanGenuchten3Coeffs_.lookup("alpha")),
    beta_(vanGenuchten3Coeffs_.lookup("beta")),
    ha_(vanGenuchten3Coeffs_.lookup("ha")),
    Ss_(vanGenuchten3Coeffs_.lookup("Ss")),
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
    ),
    H_
    (
        IOobject
        (
            "H",
            theta_.time().timeName(),
            theta_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (theta_-theta_r_)/(theta_s_-theta_r_)
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::soilModels::vanGenuchten3::read
(
    const dictionary& soilProperties
)
{
    soilModel::read(soilProperties);

    vanGenuchten3Coeffs_ = soilProperties.subDict(typeName + "Coeffs");

    vanGenuchten3Coeffs_.lookup("Ks") >> Ks_;
    vanGenuchten3Coeffs_.lookup("theta_s") >> theta_s_;
    vanGenuchten3Coeffs_.lookup("theta_r") >> theta_r_;
    vanGenuchten3Coeffs_.lookup("alpha") >> alpha_;
    vanGenuchten3Coeffs_.lookup("beta") >> beta_;
    vanGenuchten3Coeffs_.lookup("ha") >> ha_;
    vanGenuchten3Coeffs_.lookup("Ss") >> Ss_;

    return true;
}

// update soil moisture
void Foam::soilModels::vanGenuchten3::update_theta()
{  
   //update moisture content theta
   theta_ == neg(h_-ha_)*(theta_r_+(theta_s_-theta_r_)/(1.0+pow(alpha_*mag(ha_-h_),beta_)))
            +pos(h_-ha_)*theta_s_;
}

// update the soil properties
void Foam::soilModels::vanGenuchten3::correct()
{
   //update the soil moisture
   update_theta();

   //update saturation ratio Sw
   Sw_ == theta_/theta_s_;

   //update dimensionless moisture content
   H_ == (theta_-theta_r_)/(theta_s_-theta_r_);

   //update relative permeability kr
   kr_ == pow(H_,2.0);

   //update specific moisture capacity (including the specific storage term)
   Ch_ ==  alpha_*beta_*(theta_s_-theta_r_)*pow(H_,2.0)*pow((1.0-H_)/(H_+VSMALL), (beta_-1.0)/beta_)
         + Ss_*Sw_;
}

//*********************************************************** //

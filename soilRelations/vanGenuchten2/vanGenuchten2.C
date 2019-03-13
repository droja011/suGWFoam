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

#include "vanGenuchten2.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace soilModels
{
    defineTypeNameAndDebug(vanGenuchten2, 0);

    addToRunTimeSelectionTable
    (
        soilModel,
        vanGenuchten2,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soilModels::vanGenuchten2::vanGenuchten2
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
    vanGenuchten2Coeffs_(soilProperties.subDict(typeName + "Coeffs")),
    Ks_(vanGenuchten2Coeffs_.lookup("Ks")),
    theta_s_(vanGenuchten2Coeffs_.lookup("theta_s")),
    theta_r_(vanGenuchten2Coeffs_.lookup("theta_r")),
    n_(vanGenuchten2Coeffs_.lookup("n")),
    m_(1.0-1.0/n_),
    alpha_(vanGenuchten2Coeffs_.lookup("alpha")),
    a_(vanGenuchten2Coeffs_.lookup("a")),
    b_(vanGenuchten2Coeffs_.lookup("b")),
    Ss_(vanGenuchten2Coeffs_.lookup("Ss")),
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

bool Foam::soilModels::vanGenuchten2::read
(
    const dictionary& soilProperties
)
{
    soilModel::read(soilProperties);

    vanGenuchten2Coeffs_ = soilProperties.subDict(typeName + "Coeffs");

    vanGenuchten2Coeffs_.lookup("Ks") >> Ks_;
    vanGenuchten2Coeffs_.lookup("theta_s") >> theta_s_;
    vanGenuchten2Coeffs_.lookup("theta_r") >> theta_r_;
    vanGenuchten2Coeffs_.lookup("n") >> n_;
 
    m_ = 1.0 - 1.0/n_;
 
    vanGenuchten2Coeffs_.lookup("alpha") >> alpha_;
    vanGenuchten2Coeffs_.lookup("a") >> a_;
    vanGenuchten2Coeffs_.lookup("b") >> b_;
    vanGenuchten2Coeffs_.lookup("Ss") >> Ss_;

    return true;
}

// update soil moisture
void Foam::soilModels::vanGenuchten2::update_theta()
{  
   //update moisture content theta
   theta_ == neg(h_)*(theta_r_+(theta_s_-theta_r_)*pow(1.0+pow(alpha_*mag(h_),n_),-m_))
            +pos(h_)*theta_s_;
}

// update the soil properties
void Foam::soilModels::vanGenuchten2::correct()
{
   //update the soil moisture
   update_theta();

   //update saturation ratio Sw
   Sw_ == theta_/theta_s_;

   //update dimensionless moisture content
   H_ == (theta_-theta_r_)/(theta_s_-theta_r_);

   //update relative permeability kr
   kr_ == a_*pow(theta_,b_)/Ks_[0];

   //update specific moisture capacity (including the specific storage term)
   Ch_ ==  alpha_*m_*(theta_s_-theta_r_)/(1.0-m_)*pow(H_,1.0/m_)*pow(1.0-pow(H_,1.0/m_),m_)
         + Ss_*Sw_;
}

//*********************************************************** //

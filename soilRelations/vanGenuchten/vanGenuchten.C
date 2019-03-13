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

#include "vanGenuchten.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace soilModels
{
    defineTypeNameAndDebug(vanGenuchten, 0);

    addToRunTimeSelectionTable
    (
        soilModel,
        vanGenuchten,
        dictionary
    );


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soilModels::vanGenuchten::vanGenuchten
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
    vanGenuchtenCoeffs_(soilProperties.subDict(typeName + "Coeffs")),
    Ks_(vanGenuchtenCoeffs_.lookup("Ks")),
    theta_s_(vanGenuchtenCoeffs_.lookup("theta_s")),
    theta_r_(vanGenuchtenCoeffs_.lookup("theta_r")),
    n_(vanGenuchtenCoeffs_.lookup("n")),
    m_(1.0-1.0/n_),
    alpha_(vanGenuchtenCoeffs_.lookup("alpha")),
    Ss_(vanGenuchtenCoeffs_.lookup("Ss")),
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

bool Foam::soilModels::vanGenuchten::read
(
    const dictionary& soilProperties
)
{
    soilModel::read(soilProperties);

    vanGenuchtenCoeffs_ = soilProperties.subDict(typeName + "Coeffs");

    vanGenuchtenCoeffs_.lookup("Ks") >> Ks_;
    vanGenuchtenCoeffs_.lookup("theta_s") >> theta_s_;
    vanGenuchtenCoeffs_.lookup("theta_r") >> theta_r_;
    vanGenuchtenCoeffs_.lookup("n") >> n_;
 
    m_ = 1.0 - 1.0/n_;
 
    vanGenuchtenCoeffs_.lookup("alpha") >> alpha_;
    vanGenuchtenCoeffs_.lookup("Ss") >> Ss_;

    return true;
}

// update soil moisture
void Foam::soilModels::vanGenuchten::update_theta()
{  
   //update moisture content theta
   theta_ == neg(h_)*(theta_r_+(theta_s_-theta_r_)*pow(1.0+pow(alpha_*mag(h_),n_),-m_))
            +pos(h_)*theta_s_;
}

// update the soil properties
void Foam::soilModels::vanGenuchten::correct()
{
   //update the soil moisture
   update_theta();

   //update saturation ratio Sw
   Sw_ == theta_/theta_s_;

   //update dimensionless moisture content
   H_ == (theta_-theta_r_)/(theta_s_-theta_r_);

   //update relative permeability kr
   kr_ == neg(h_)*pow(1.0+pow(alpha_*mag(h_),n_),-m_/2.0)
                 *pow(1.0-pow(alpha_*mag(h_),n_-1)*pow(1.0+pow(alpha_*mag(h_),n_),-m_),2.0)
         +pos(h_)*1.0;

   //update specific moisture capacity (including the specific storage term)
   Ch_ == alpha_*m_*(theta_s_-theta_r_)/(1.0-m_)*pow(H_,1.0/m_)*pow(1.0-pow(H_,1.0/m_),m_)
        +Ss_*Sw_;
}

} //End of namespace soilModels
} //End of namespace Foam

//*********************************************************** //

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

#include "PaniconiEtAl1991.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace soilModels
{
    defineTypeNameAndDebug(PaniconiEtAl1991, 0);

    addToRunTimeSelectionTable
    (
        soilModel,
        PaniconiEtAl1991,
        dictionary
    );


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soilModels::PaniconiEtAl1991::PaniconiEtAl1991
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
    PaniconiEtAl1991Coeffs_(soilProperties.subDict(typeName + "Coeffs")),
    Ks_(PaniconiEtAl1991Coeffs_.lookup("Ks")),
    theta_s_(PaniconiEtAl1991Coeffs_.lookup("theta_s")),
    theta_r_(PaniconiEtAl1991Coeffs_.lookup("theta_r")),
    h_s_(PaniconiEtAl1991Coeffs_.lookup("h_s")),
    h_0_(PaniconiEtAl1991Coeffs_.lookup("h_0")),
    n_(PaniconiEtAl1991Coeffs_.lookup("n")),
    m_(1.0-1.0/n_),
    beta_0_(pow((h_0_/h_s_).value(),n_.value())),
    Ss_(PaniconiEtAl1991Coeffs_.lookup("Ss")),
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
    beta_
    (
        IOobject
        (
            "beta",
            theta_.time().timeName(),
            theta_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        theta_.mesh(),
        dimensionedScalar("beta", dimless, 0.0) 
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::soilModels::PaniconiEtAl1991::read
(
    const dictionary& soilProperties
)
{
    soilModel::read(soilProperties);

    PaniconiEtAl1991Coeffs_ = soilProperties.subDict(typeName + "Coeffs");

    PaniconiEtAl1991Coeffs_.lookup("Ks") >> Ks_;
    PaniconiEtAl1991Coeffs_.lookup("theta_s") >> theta_s_;
    PaniconiEtAl1991Coeffs_.lookup("theta_r") >> theta_r_;
    PaniconiEtAl1991Coeffs_.lookup("h_s") >> h_s_;
    PaniconiEtAl1991Coeffs_.lookup("h_0") >> h_0_;
    PaniconiEtAl1991Coeffs_.lookup("n") >> n_;
 
    m_ = 1.0 - 1.0/n_;
    
    beta_0_ = pow(h_0_/h_s_,n_).value();
 
    PaniconiEtAl1991Coeffs_.lookup("Ss") >> Ss_;

    return true;
}

// update soil moisture
void Foam::soilModels::PaniconiEtAl1991::update_theta()
{  
   //update dimensionless moisture content
   beta_ == pow(mag(h_/h_s_),n_);

   //update moisture content theta
   theta_ == neg(h_-h_0_)*(theta_r_+(theta_s_-theta_r_)*pow(1.0+beta_,-m_))
            +pos(h_-h_0_)*(theta_r_+(theta_s_-theta_r_)*pow(1.0+beta_0_,-m_)+Ss_*(h_-h_0_));
}

// update the soil properties
void Foam::soilModels::PaniconiEtAl1991::correct()
{
   //update the soil moisture
   update_theta();

   //update saturation ratio Sw
   Sw_ == theta_/theta_s_;

   //update relative permeability kr
   kr_ == neg(h_)*pow(1.0+beta_, -5.0*m_/2.0)
                 *pow( (pow(1.0+beta_, m_)-pow(beta_,m_)), 2)
         +pos(h_)*1.0;

   //Info << "Kr = " << kr_ << endl;

   //update specific moisture capacity (including the specific storage term)
   Ch_ == neg(h_-h_0_)*(n_-1.0)*(theta_s_-theta_r_)*pow(mag(h_),(n_-1.0))/pow(mag(h_s_),n_)/
                       pow(1.0+beta_,(m_+1.0))
         +pos(h_-h_0_)*Ss_;
}

} //End of namespace soilModels
} //End of namespace Foam
//*********************************************************** //

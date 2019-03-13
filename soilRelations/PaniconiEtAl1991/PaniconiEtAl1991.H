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

Class
    Foam::soilModels::PaniconiEtAl1991

Description
    soil constitutive relationship as in Paniconi et al. (1991).

    The Paniconi et al. (1991) model is a modified version of van Genuchten and
    Nielsen (1985). The difference is that Paniconi et al. (1991) model permits
    a nonzero value of specific moisture capacity to occure in the satureated zone. 
    For Ss=0 and h0=0, it recovers the orgiginal van Genuchten and Nielsen (1985) model. 

    Reference:
    Paniconi, C., Aldama, A.A., and Wood, E.F. (1991). "Numerical evalutation of iterative
    and noniterative methods for the solution of the nonlinear Richards equation", WRR 27(6):
    1147-1163.

    van Genuchten, M.T. and Nielsen, D.R. (1985). "On describing and predicting the hydraulic
    properties of unsaturated soils. Ann. Geophys. 3(5):615-628.

    van Genuchten, M.T. (1980). "A closed-form equation for predicting the hydraulic
    conductivity of unsaturated soils.". Soil Science Society of America Journal,
    44(5): 892-898


SourceFiles
    PaniconiEtAl1991.C

\*---------------------------------------------------------------------------*/

#ifndef PaniconiEtAl1991_H
#define PaniconiEtAl1991_H

#include "soilModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace soilModels
{

/*---------------------------------------------------------------------------*\
                           Class PaniconiEtAl1991 Declaration
\*---------------------------------------------------------------------------*/

class PaniconiEtAl1991
:
    public soilModel
{
    // Private data

        dictionary PaniconiEtAl1991Coeffs_;

        dimensionedTensor Ks_;        //saturated hydraulic conductivity
        dimensionedScalar theta_s_;   //saturated moisture content
        dimensionedScalar theta_r_;   //residual moisture content
        dimensionedScalar h_s_;       //bubbling (or air entry) pressure head
        dimensionedScalar h_0_;       //parameter (can be calculated based on Ss
        dimensionedScalar n_;         //van Genuchten parameter
        dimensionedScalar m_;         //van Genuchten paramter m=1-1/n
        scalar beta_0_;               //beta_0 = pow(h_0/h_s,n)
        dimensionedScalar Ss_;        //specific storage
       
        volScalarField Sw_;           //water saturation ratio Sw=theta/porosity
        volScalarField beta_;         // =pow(abs(h/hs),n)

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("PaniconiEtAl1991");


    // Constructors

        //- Construct from components
        PaniconiEtAl1991
        (
            const word& name, 
            const dictionary& soilProperties,
            volScalarField& h,
            volScalarField& theta,
            volScalarField& kr,
            volScalarField& Ch
        );


    // Destructor

        ~PaniconiEtAl1991()
        {}


    // Member Functions
       
        //- update the soil moisture content
        void update_theta();

        //- Correct the soil properties 
        void correct();

        //- Read transportProperties dictionary
        bool read(const dictionary& soilProperties);

        dimensionedTensor Ks()
        {
            return Ks_;
        }

        dimensionedScalar theta_s()
        {
            return theta_s_;
        }

        dimensionedScalar Ss()
        {
            return Ss_;
        }

        volScalarField Sw()
        {
            return Sw_;
        }

        volScalarField beta()
        {
            return beta_;
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace soilModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

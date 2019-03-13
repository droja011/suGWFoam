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
    Foam::soilModels::vanGenuchten3

Description
    soil moisture constitutive relationship as a variant of van Genuchten (1980).

SourceFiles
    vanGenuchten3.C

\*---------------------------------------------------------------------------*/

#ifndef vanGenuchten3_H
#define vanGenuchten3_H

#include "soilModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace soilModels
{

/*---------------------------------------------------------------------------*\
                           Class vanGenuchten3 Declaration
\*---------------------------------------------------------------------------*/

class vanGenuchten3
:
    public soilModel
{
    // Private data

        dictionary vanGenuchten3Coeffs_;

        dimensionedTensor Ks_;
        dimensionedScalar theta_s_;
        dimensionedScalar theta_r_;
        dimensionedScalar alpha_;
        dimensionedScalar beta_;
        dimensionedScalar ha_;
        dimensionedScalar Ss_;

        volScalarField Sw_; //water saturation ration Sw=theta/porosity
        volScalarField H_; //dimensionless moisture content

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("vanGenuchten3");


    // Constructors

        //- Construct from components
        vanGenuchten3
        (
            const word& name, 
            const dictionary& soilProperties,
            volScalarField& h,
            volScalarField& theta,
            volScalarField& kr,
            volScalarField& Ch
        );


    // Destructor

        ~vanGenuchten3()
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

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace soilModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

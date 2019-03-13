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
    Foam::soilModels::Haverkamp

Description
    soil constitutive relationship as in Haverkamp et al. (1977).
   
    Reference:
    Haverkamp R, Vauclin M, Touma J, Wierenga P, Vachaud G. Comparision of numerical
    simulation models for one dimensional infiltration. Soil Sci. Soc. Am J. 1977,
    41:285-294

SourceFiles
    Haverkamp.C

\*---------------------------------------------------------------------------*/

#ifndef Haverkamp_H
#define Haverkamp_H

#include "soilModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace soilModels
{

/*---------------------------------------------------------------------------*\
                           Class Haverkamp Declaration
\*---------------------------------------------------------------------------*/

class Haverkamp
:
    public soilModel
{
    // Private data

        dictionary HaverkampCoeffs_;

        dimensionedTensor Ks_;
        dimensionedScalar theta_s_;
        dimensionedScalar theta_r_;
        dimensionedScalar alpha_;
        dimensionedScalar beta_;
        dimensionedScalar gamma_;
        dimensionedScalar A_;
        dimensionedScalar Ss_;

        volScalarField Sw_; //water saturation ration Sw=theta/porosity

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("Haverkamp");


    // Constructors

        //- Construct from components
        Haverkamp
        (
            const word& name, 
            const dictionary& soilProperties,
            volScalarField& h,
            volScalarField& theta,
            volScalarField& kr,
            volScalarField& Ch
        );


    // Destructor

        ~Haverkamp()
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

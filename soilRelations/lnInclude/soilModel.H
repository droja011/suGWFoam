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

Namespace
    Foam::soilModels

Description
    A namespace for various soilModel implementations.

Class
    Foam::soilModel

Description
    An abstract base class for soilModels

    Constitutive relations for soils: theta(h) and K(h)

SourceFiles
    soilModel.C
    newSoilModel.C

\*---------------------------------------------------------------------------*/

#ifndef soilModel_H
#define soilModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class soilModel Declaration
\*---------------------------------------------------------------------------*/

class soilModel
{

protected:

    // Protected data

        word name_;
        dictionary soilProperties_;

        volScalarField& h_;           //head
        volScalarField& theta_;       //moisture content
        volScalarField& kr_;          //relative permeability
        volScalarField& Ch_;          //specific moisture capacity


    // Private Member Functions

        //- Disallow copy construct
        soilModel(const soilModel&);

        //- Disallow default bitwise assignment
        void operator=(const soilModel&);


public:

    //- Runtime type information
    TypeName("soilModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            soilModel,
            dictionary,
            (
                const word& name,
                const dictionary& soilProperties,
                volScalarField& h,
                volScalarField& theta,
                volScalarField& kr,
                volScalarField& Ch
            ),
            (name, soilProperties, h, theta, kr, Ch)
        );


    // Selectors

        //- Return a reference to the selected soil model
        static autoPtr<soilModel> New
        (
            const word& name,
            const dictionary& soilProperties,
            volScalarField& h,
            volScalarField& theta,
            volScalarField& kr,
            volScalarField& Ch
        );


    // Constructors

        //- Construct from components
        soilModel
        (
            const word& name,
            const dictionary& soilProperties,
            volScalarField& h,
            volScalarField& theta,
            volScalarField& kr,
            volScalarField& Ch
        );


    // Destructor

        virtual ~soilModel()
        {}


    // Member Functions

        //- Return the soil properties dictionary
        const dictionary& soilProperties() const
        {
            return soilProperties_;
        }

        //- Correct the soil moisture content 
        virtual void update_theta() = 0;

        //- Correct the soil properties 
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& soilProperties) = 0;

        virtual dimensionedTensor Ks() = 0;
        virtual dimensionedScalar theta_s() = 0;
        virtual dimensionedScalar Ss() = 0;
        virtual volScalarField Sw() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

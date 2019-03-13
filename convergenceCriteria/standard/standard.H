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
    Foam::convergenceCriterion::standard

Description
    standard convergence criterion 

SourceFiles
    standard.C

\*---------------------------------------------------------------------------*/

#ifndef standard_H
#define standard_H

#include "convergenceCriterion.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace convergenceCriteria
{

/*---------------------------------------------------------------------------*\
                           Class standard Declaration
\*---------------------------------------------------------------------------*/

class standard
:
    public convergenceCriterion 
{
    // Private data

        dictionary standardCoeffs_;

        scalar delta_h_;


    // Private Member Functions

public:

    //- Runtime type information
    TypeName("standard");


    // Constructors

        //- Construct from components
        standard
        (
            const word& name, 
            const dictionary& convergenceProperties,
            volScalarField& h,
            volScalarField& h_n,
            volScalarField& theta,
            volScalarField& theta_n
        );


    // Destructor

        ~standard()
        {}


    // Member Functions
       
        bool convergent();

        //- Read transportProperties dictionary
        bool read(const dictionary& convergenceProperties);


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace convergenceCriteria
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

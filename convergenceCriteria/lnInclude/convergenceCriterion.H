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
    Foam::convergenceCriterion

Description
    A namespace for various convergenceCriterion implementations.

Class
    Foam::convergenceCriterion

Description
    An abstract base class for convergenceCriterion

    See the reference:
    Huang, K., Mohanty, B.P., and van Genuchten, M.T. (1996). "A new convergence
    criterion for the modified Picard iteration method to solve the variably 
    saturated flow equation." Journal of Hydrology, 178(1-4): 69-91

SourceFiles
    convergenceCriterion.C
    newConvergenceCriterion.C

\*---------------------------------------------------------------------------*/

#ifndef convergenceCriterion_H
#define convergenceCriterion_H

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
                Class convergenceCriterion Declaration
\*---------------------------------------------------------------------------*/

class convergenceCriterion
{

protected:

    // Protected data

        word name_;
        dictionary convergenceProperties_;

        volScalarField& h_;             //head
        volScalarField& h_n_;           //head
        volScalarField& theta_;         //moisture content
        volScalarField& theta_n_;       //moisture content


    // Private Member Functions

        //- Disallow copy construct
        convergenceCriterion(const convergenceCriterion&);

        //- Disallow default bitwise assignment
        void operator=(const convergenceCriterion&);


public:

    //- Runtime type information
    TypeName("convergenceCriterion");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            convergenceCriterion,
            dictionary,
            (
                const word& name,
                const dictionary& convergenceProperties,
                volScalarField& h,
                volScalarField& h_n,
                volScalarField& theta,
                volScalarField& theta_n
            ),
            (name, convergenceProperties, h, h_n, theta, theta_n)
        );


    // Selectors

        //- Return a reference to the selected soil model
        static autoPtr<convergenceCriterion> New
        (
            const word& name,
            const dictionary& convergenceProperties,
            volScalarField& h,
            volScalarField& h_n,
            volScalarField& theta,
            volScalarField& theta_n
        );


    // Constructors

        //- Construct from components
        convergenceCriterion
        (
            const word& name,
            const dictionary& convergenceProperties,
            volScalarField& h,
            volScalarField& h_n,
            volScalarField& theta,
            volScalarField& theta_n
        );


    // Destructor

        virtual ~convergenceCriterion()
        {}


    // Member Functions

        //- Return the convergence properties dictionary
        const dictionary& convergenceProperties() const
        {
            return convergenceProperties_;
        }

        //- check whether convergent 
        virtual bool convergent() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& convergenceProperties) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

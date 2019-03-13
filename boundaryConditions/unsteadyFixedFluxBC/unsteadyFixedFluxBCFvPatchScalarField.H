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
    Foam::unsteadyFixedFluxBCFvPatchScalarField.H

Description
    Set the pressure head gradient boundary condition appropriately for fixed flux.

    The flux changes with time. 

    This special B.C. is hard-coded for the Paniconi et al. (1991) cases. For other
    unsteady flux B.C., it should be easy to adapt. 

SourceFiles
    unsteadyFixedFluxBCFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef unsteadyFixedFluxBCFvPatchScalarField_H
#define unsteadyFixedFluxBCFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class unsteadyFixedFluxBCFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class unsteadyFixedFluxBCFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        scalarField Ks_normal_;

public:

    //- Runtime type information
    TypeName("unsteadyFixedFluxBC");


    // Constructors

        //- Construct from patch and internal field
        unsteadyFixedFluxBCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        unsteadyFixedFluxBCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  unsteadyFixedFluxBCFvPatchScalarField onto a new patch
        unsteadyFixedFluxBCFvPatchScalarField
        (
            const unsteadyFixedFluxBCFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        unsteadyFixedFluxBCFvPatchScalarField
        (
            const unsteadyFixedFluxBCFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new unsteadyFixedFluxBCFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        unsteadyFixedFluxBCFvPatchScalarField
        (
            const unsteadyFixedFluxBCFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new unsteadyFixedFluxBCFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


        //- access
        virtual const scalarField& Ks_normal() const
        {
                return Ks_normal_;
        }

        virtual  scalarField& Ks_normal()
        {
                return Ks_normal_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

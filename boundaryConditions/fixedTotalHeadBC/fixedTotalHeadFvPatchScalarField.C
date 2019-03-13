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

#include "fixedTotalHeadFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedTotalHeadFvPatchScalarField::fixedTotalHeadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    h0_(p.size(), 0.0)
{}


Foam::fixedTotalHeadFvPatchScalarField::fixedTotalHeadFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    h0_("h0", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(h0_);
    }
}


Foam::fixedTotalHeadFvPatchScalarField::fixedTotalHeadFvPatchScalarField
(
    const fixedTotalHeadFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    h0_(ptf.h0_, mapper)
{}


Foam::fixedTotalHeadFvPatchScalarField::fixedTotalHeadFvPatchScalarField
(
    const fixedTotalHeadFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    h0_(tppsf.h0_)
{}


Foam::fixedTotalHeadFvPatchScalarField::fixedTotalHeadFvPatchScalarField
(
    const fixedTotalHeadFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    h0_(tppsf.h0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedTotalHeadFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    h0_.autoMap(m);
}


void Foam::fixedTotalHeadFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const fixedTotalHeadFvPatchScalarField& tiptf =
        refCast<const fixedTotalHeadFvPatchScalarField>(ptf);

    h0_.rmap(tiptf.h0_, addr);
}


void Foam::fixedTotalHeadFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& z =
        patch().lookupPatchField<volScalarField, scalar>("z");

//  Info << "z " << z << endl;
//  Info << "h0_ " << h0_ << endl;

    operator== 
    (
        h0_-z
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}



void Foam::fixedTotalHeadFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    h0_.writeEntry("h0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedTotalHeadFvPatchScalarField
    );
}

// ************************************************************************* //

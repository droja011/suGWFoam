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

#include "fixedFluxBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedFluxBCFvPatchScalarField::fixedFluxBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    Ks_normal_(p.size(), 0.0),
    flux_(p.size(), 0.0)
{}


fixedFluxBCFvPatchScalarField::
fixedFluxBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    Ks_normal_("Ks_normal", dict, p.size()),
    flux_("flux", dict, p.size())
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


fixedFluxBCFvPatchScalarField::
fixedFluxBCFvPatchScalarField
(
    const fixedFluxBCFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    Ks_normal_(ptf.Ks_normal_, mapper),
    flux_(ptf.flux_, mapper)
{
}


fixedFluxBCFvPatchScalarField::
fixedFluxBCFvPatchScalarField
(
    const fixedFluxBCFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    Ks_normal_(ptf.Ks_normal_),
    flux_(ptf.flux_)
{
}


fixedFluxBCFvPatchScalarField::
fixedFluxBCFvPatchScalarField
(
    const fixedFluxBCFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    Ks_normal_(ptf.Ks_normal_),
    flux_(ptf.flux_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedFluxBCFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& gradz =
        patch().lookupPatchField<volVectorField, vector>("grad_z");

    const fvPatchField<scalar>& kr =
        patch().lookupPatchField<volScalarField, scalar>("kr");

 
//   Info << "gradz= " << gradz << endl;
//   Info << "nf= " << patch().nf() << endl;
//   Info << "gradz&nf= " << (gradz & patch().nf()) << endl;

     //assume the influx boundary always saturated 
//   gradient() = -flux_/Ks_normal_ - (gradz & patch().nf());
     gradient() = -flux_/Ks_normal_/kr - (gradz & patch().nf());
//   gradient() = -(-this->db().time().value())/64.0/Ks_normal_/kr - (gradz & patch().nf());

//   Info << "gradient()= " << gradient() << endl;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void fixedFluxBCFvPatchScalarField::write(Ostream& os) const
{
//  Info << "patch name = " << patch().name() << endl;
//  Info << "gradient()= " << gradient() << endl;

    fixedGradientFvPatchScalarField::write(os);
    Ks_normal_.writeEntry("Ks_normal",os);
    flux_.writeEntry("flux",os); 
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedFluxBCFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

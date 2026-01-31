/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pasquillAtmBoundaryLayerInletKFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pasquillAtmBoundaryLayerInletKFvPatchScalarField::
pasquillAtmBoundaryLayerInletKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    pasquillAtmBoundaryLayer(iF.time(), p.patch())
{}


pasquillAtmBoundaryLayerInletKFvPatchScalarField::
pasquillAtmBoundaryLayerInletKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    pasquillAtmBoundaryLayer(iF.time(), p.patch(), dict)
{
    phiName_ = dict.getOrDefault<word>("phi", "phi");

    refValue() = k(patch());
    refGrad() = 0;
    valueFraction() = 1;

    if (!initABL_)
    {
        this->readValueEntry(dict, IOobjectOption::MUST_READ);
    }
    else
    {
        scalarField::operator=(refValue());
        initABL_ = false;
    }
}


pasquillAtmBoundaryLayerInletKFvPatchScalarField::
pasquillAtmBoundaryLayerInletKFvPatchScalarField
(
    const pasquillAtmBoundaryLayerInletKFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(psf, p, iF, mapper),
    pasquillAtmBoundaryLayer(psf, p, mapper)
{}


pasquillAtmBoundaryLayerInletKFvPatchScalarField::
pasquillAtmBoundaryLayerInletKFvPatchScalarField
(
    const pasquillAtmBoundaryLayerInletKFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(psf, iF),
    pasquillAtmBoundaryLayer(psf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pasquillAtmBoundaryLayerInletKFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    refValue() = k(patch());

    inletOutletFvPatchScalarField::updateCoeffs();
}


void pasquillAtmBoundaryLayerInletKFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);
    pasquillAtmBoundaryLayer::autoMap(m);
}


void pasquillAtmBoundaryLayerInletKFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(psf, addr);

    const pasquillAtmBoundaryLayerInletKFvPatchScalarField& blpsf =
        refCast<const pasquillAtmBoundaryLayerInletKFvPatchScalarField>(psf);

    pasquillAtmBoundaryLayer::rmap(blpsf, addr);
}


void pasquillAtmBoundaryLayerInletKFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    pasquillAtmBoundaryLayer::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    pasquillAtmBoundaryLayerInletKFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

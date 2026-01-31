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

#include "pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField::
pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    pasquillAtmBoundaryLayer(iF.time(), p.patch())
{}


pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField::
pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    pasquillAtmBoundaryLayer(iF.time(), p.patch(), dict)
{
    fvPatchVectorField::operator==(U(patch()));
}


pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField::
pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField
(
    const pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(pvf, p, iF, mapper),
    pasquillAtmBoundaryLayer(pvf, p, mapper)
{}


pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField::
pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField
(
    const pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    pasquillAtmBoundaryLayer(pvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fvPatchVectorField::operator==(U(patch()));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    pasquillAtmBoundaryLayer::autoMap(m);
}


void pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& pvf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(pvf, addr);

    const pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField& blpvf =
        refCast<const pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField>(pvf);

    pasquillAtmBoundaryLayer::rmap(blpvf, addr);
}


void pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    pasquillAtmBoundaryLayer::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
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

Class
    Foam::pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a log-law type ground-normal inlet
    boundary condition for the specific dissipation rate (i.e. \c omega)
    for homogeneous, two-dimensional, dry-air, equilibrium and neutral
    atmospheric boundary layer modelling.

    The ground-normal \c omega profile expression:

        \f[
            \omega = \frac{u^*}{\kappa \sqrt{C_\mu}} \frac{1}{z - d + z_0}
        \f]

    where
    \vartable
      \omega | Ground-normal specific dissipation rate profile [m^2/s^3]
      u^*    | Friction velocity                               [m/s]
      \kappa | von Kármán constant                             [-]
      C_\mu  | Empirical model constant                        [-]
      z      | Ground-normal coordinate component              [m]
      d      | Ground-normal displacement height               [m]
      z_0    | Aerodynamic roughness length                    [m]
    \endvartable

    Required fields:
    \verbatim
        omega   | Specific dissipation rate                    [1/s]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        // Mandatory entries (unmodifiable)
        type            pasquillAtmBoundaryLayerInletOmega;

        // Mandatory/Optional (inherited) entries (unmodifiable)
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type | Reqd | Deflt
      type      | Type name: pasquillAtmBoundaryLayerInletOmega  | word | yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link pasquillAtmBoundaryLayer.H \endlink
     - \link inletOutletFvPatchField.H \endlink

See also
    - Foam::pasquillAtmBoundaryLayer
    - Foam::pasquillAtmBoundaryLayerInletVelocityFvPatchVelocityField
    - Foam::pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField
    - Foam::pasquillAtmBoundaryLayerInletKFvPatchScalarField

SourceFiles
    pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField_H
#define pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField_H

#include "fvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "pasquillAtmBoundaryLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField
:
    public inletOutletFvPatchScalarField,
    public pasquillAtmBoundaryLayer
{

public:

    //- Runtime type information
    TypeName("pasquillAtmBoundaryLayerInletOmega");


    // Constructors

        //- Construct from patch and internal field
        pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField onto a new patch
        pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField
        (
            const pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField
        (
            const pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField&,
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
                new pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

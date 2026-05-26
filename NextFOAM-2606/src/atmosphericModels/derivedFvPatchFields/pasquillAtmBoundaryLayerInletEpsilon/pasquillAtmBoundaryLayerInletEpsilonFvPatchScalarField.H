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

Class
    Foam::pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a log-law type ground-normal inlet
    boundary condition for the turbulent kinetic energy dissipation rate
    (i.e. \c epsilon) for homogeneous, two-dimensional, dry-air, equilibrium
    and neutral atmospheric boundary layer modelling.

    The ground-normal \c epsilon profile expression:

        \f[
            \epsilon = \frac{(u^*)^3}{\kappa (z - d + z_0)}
            \sqrt{C_1 \ln \left( \frac{z - d + z_0}{z_0} \right) + C_2}
        \f]

    where
    \vartable
      \epsilon | Ground-normal TKE dissipation rate profile     [m^2/s^3]
      u^*      | Friction velocity                              [m/s]
      \kappa   | von Kármán constant                            [-]
      z        | Ground-normal coordinate component             [m]
      z_0      | Aerodynamic roughness length                   [m]
      d        | Ground-normal displacement height              [m]
      C_1      | Curve-fitting coefficient for \c YGCJ profiles [-]
      C_2      | Curve-fitting coefficient for \c YGCJ profiles [-]
    \endvartable

    Required fields:
    \verbatim
        epsilon   | Turbulent kinetic energy dissipation rate   [m2/s3]
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        // Mandatory entries (unmodifiable)
        type            pasquillAtmBoundaryLayerInletEpsilon;

        // Mandatory/Optional (inherited) entries (unmodifiable)
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type | Reqd | Deflt
      type      | Type name: pasquillAtmBoundaryLayerInletEpsilon | word | yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link pasquillAtmBoundaryLayer.H \endlink
     - \link inletOutletFvPatchField.H \endlink

See also
    - Foam::pasquillAtmBoundaryLayer
    - Foam::pasquillAtmBoundaryLayerInletKFvPatchScalarField
    - Foam::pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField
    - Foam::pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField

SourceFiles
    pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField_H
#define pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField_H

#include "fvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "pasquillAtmBoundaryLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField
:
    public inletOutletFvPatchScalarField,
    public pasquillAtmBoundaryLayer
{

public:

    //- Runtime type information
    TypeName("pasquillAtmBoundaryLayerInletEpsilon");


    // Constructors

        //- Construct from patch and internal field
        pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField onto a new patch
        pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField
        (
            const pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField&,
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
                new pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField(*this, iF)
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

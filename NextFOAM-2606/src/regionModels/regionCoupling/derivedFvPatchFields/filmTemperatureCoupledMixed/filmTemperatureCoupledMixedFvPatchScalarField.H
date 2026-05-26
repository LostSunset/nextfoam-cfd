/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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
    Foam::filmTemperatureCoupledMixedFvPatchScalarField

Description
    Mixed temperature boundary condition for coupled walls when a film region
    model is applied on one side.

    Example usage:
    \verbatim
    // Note: this patch belongs to the film region
    patchNameOnFilmRegion
    {
        type            filmTemperatureCoupledMixed;
        filmRegion      filmRegionName;
        coupledPatch    patchNameOnCoupledRegion;
        kappaMethod     lookup;
        kappa           kappa;
        Tnbr            T;
        value           $internalField;
    }

    // Note: this patch does not belong to the film primary region
    patchNameOnCoupledRegion
    {
        type            filmTemperatureCoupledMixed;
        filmRegion      filmRegionName;
        coupledPatch    patchNameOnFilmRegion;
        kappaMethod     fluidThermo;
        Tnbr            Tf;
        value           $internalField;
    }

    // Conditional entries

        thicknessLayers <scalarList>;
        kappaLayers     <scalarList>;

    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.
    It calculates local field as:

    Notes:
    - \c kappaMethod and \c kappa are inherited from temperatureCoupledBase.


See also
    Foam::temperatureCoupledBase

SourceFiles
    filmTemperatureCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef filmTemperatureCoupledMixedFvPatchScalarField_H
#define filmTemperatureCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "thermoSingleLayer.H"
#include "mappedFieldFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class filmTemperatureCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class filmTemperatureCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    typedef Foam::regionModels::surfaceFilmModels::thermoSingleLayer
        filmModelType;

    typedef filmTemperatureCoupledMixedFvPatchScalarField myType;


private:

    // Private Data

        //- Name of film region
        const word filmRegionName_;

        //- Name of temperature coupled patch on the neighbour region
        //  (Note: this is not necessarily the same as the patch name of the
        //  samplePatch in polyMesh)
        const word coupledPatchName_;

        //- Name of field on the neighbour region
        const word TnbrName_;
    
        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;

        //- Pointer to AMI interpolator
        mutable autoPtr<AMIPatchToPatchInterpolation> AMIPtr_;


    // Private Member Functions

        //- Retrieve film model from the database, or FatalError
        const filmModelType& filmModel() const;

        //- Reset the AMI
        void resetAMI
        (
            const polyPatch& source,
            const polyPatch& target,
            const word& sourceMeshName
        ) const;


public:

    //- Runtime type information
    TypeName("filmTemperatureCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        filmTemperatureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        filmTemperatureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        filmTemperatureCoupledMixedFvPatchScalarField
        (
            const
            filmTemperatureCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        filmTemperatureCoupledMixedFvPatchScalarField
        (
            const filmTemperatureCoupledMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

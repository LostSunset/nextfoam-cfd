/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::functionObjects::wallRadHeatFluxRatio

Group
    grpForcesFunctionObjects

Description
    Computes the wall radiative heat flux ratio relative to total heat flux at 
    selected wall patches.

    \table
      Operand        | Type              | Location
      input          | -                 | -
      output file    | dat               <!--
                 --> | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<field\>
      output field   | volScalarField (only boundaryField) <!--
                 --> | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    wallRadHeatFluxRatio1
    {
        // Mandatory entries (unmodifiable)
        type        wallRadHeatFluxRatio;
        libs        (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        patches     (<patch1> ... <patchN>); // (wall1 "(wall2|wall3)");
        qr          qr;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: wallRadHeatFluxRatio               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      patches    | Names of operand patches   | wordList | no | all wall patches
      qr         | Name of radiative heat flux field       | word | no | qr
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        <solver> -postProcess -func wallRadHeatFluxRatio
    \endverbatim

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::wallRadHeatFluxRatio

SourceFiles
    wallRadHeatFluxRatio.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_wallRadHeatFluxRatio_H
#define Foam_functionObjects_wallRadHeatFluxRatio_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class wallRadHeatFluxRatio Declaration
\*---------------------------------------------------------------------------*/

class wallRadHeatFluxRatio
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Wall patches to process (optionally filtered by name)
        labelList patchIDs_;

        //- Name of radiative heat flux name
        word qrName_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the radiative heat flux ratio
        void calcFluxRatio
        (
            const volScalarField& alpha,
            const volScalarField& he,
            volScalarField& wallRadHeatFluxRatio
        );


public:

    //- Runtime type information
    TypeName("wallRadHeatFluxRatio");


    // Constructors

        //- Construct from Time and dictionary
        wallRadHeatFluxRatio
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        wallRadHeatFluxRatio(const wallRadHeatFluxRatio&) = delete;

        //- No copy assignment
        void operator=(const wallRadHeatFluxRatio&) = delete;


    //- Destructor
    virtual ~wallRadHeatFluxRatio() = default;


    // Member Functions

        //- Read the wallRadHeatFluxRatio data
        virtual bool read(const dictionary& dict);

        //- Calculate the wall radiative heat flux ratio
        virtual bool execute();

        //- Write the wall radiative heat flux ratio
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

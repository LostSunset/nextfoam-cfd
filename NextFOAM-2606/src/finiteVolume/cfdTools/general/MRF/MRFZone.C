/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "MRFZone.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "faceSet.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "fvc.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MRFZone, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::MRFZone::setMRFFaces()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Type per face:
    //  0:not in zone
    //  1:moving with frame
    //  2:other
    labelList faceType(mesh_.nFaces(), Zero);

    // Determine faces in cell zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (without constructing cells)

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Cells in zone
    boolList zoneCell(mesh_.nCells(), false);

    if (cellZoneID_ != -1)
    {
        const labelList& cellLabels = mesh_.cellZones()[cellZoneID_];
        forAll(cellLabels, i)
        {
            zoneCell[cellLabels[i]] = true;
        }
    }


    // label nZoneFaces = 0;

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (zoneCell[own[facei]] || zoneCell[nei[facei]])
        {
            faceType[facei] = 1;
            // ++nZoneFaces;
        }
    }


    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (rotatingPatchNames_.contains(pp.name()))
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;

                if (!zoneCell[own[facei]])
                {
                    faceType[facei] = 3;
                }
                else
                {
                    faceType[facei] = 1;
                }
            }
        }

        if 
        (
            pp.coupled()
         || excludedPatchLabels_.contains(patchi)
         || ( 
                mesh_.foundObject<volScalarField>("rUAp")
             && (
                    pp.type() == "patch"
                 || isA<cyclicAMIPolyPatch>(pp)
                 || isA<symmetryPolyPatch>(pp)
                 || isA<symmetryPlanePolyPatch>(pp)
                )
            )
        )
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;

                if (zoneCell[own[facei]])
                {
                    faceType[facei] = 2;
                    // ++nZoneFaces;
                }
            }
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;

                if (zoneCell[own[facei]])
                {
                    faceType[facei] = 1;
                    // ++nZoneFaces;
                }
            }
        }
    }

    // Synchronize the faceType across processor patches
    syncTools::syncFaceList(mesh_, faceType, maxEqOp<label>());

    // Now we have for faceType:
    //  0   : face not in cellZone
    //  1   : internal face or normal patch face
    //  2   : coupled patch face or excluded patch face
    //  3   : face not in cellZone but rotating with MRF zone (by Gill)

    // Sort into lists per patch.

    internalFaces_.setSize(mesh_.nFaces());
    label nInternal = 0;

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceType[facei] == 1)
        {
            internalFaces_[nInternal++] = facei;
        }
    }
    internalFaces_.setSize(nInternal);

    labelList nIncludedFaces(patches.size(), Zero);
    labelList nExcludedFaces(patches.size(), Zero);
    labelList nRotatingFaces(patches.size(), Zero);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, patchFacei)
        {
            label facei = pp.start() + patchFacei;

            if (faceType[facei] == 1)
            {
                nIncludedFaces[patchi]++;
            }
            else if (faceType[facei] == 2)
            {
                nExcludedFaces[patchi]++;
            }
            else if (faceType[facei] == 3)
            {
                nRotatingFaces[patchi]++;
            }
        }
    }

    includedFaces_.setSize(patches.size());
    excludedFaces_.setSize(patches.size());
    rotatingFaces_.setSize(patches.size());
    forAll(nIncludedFaces, patchi)
    {
        includedFaces_[patchi].setSize(nIncludedFaces[patchi]);
        excludedFaces_[patchi].setSize(nExcludedFaces[patchi]);
        rotatingFaces_[patchi].setSize(nRotatingFaces[patchi]);
    }
    nIncludedFaces = 0;
    nExcludedFaces = 0;
    nRotatingFaces = 0;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        forAll(pp, patchFacei)
        {
            label facei = pp.start() + patchFacei;

            if (faceType[facei] == 1)
            {
                includedFaces_[patchi][nIncludedFaces[patchi]++] = patchFacei;
            }
            else if (faceType[facei] == 2)
            {
                excludedFaces_[patchi][nExcludedFaces[patchi]++] = patchFacei;
            }
            else if (faceType[facei] == 3)
            {
                rotatingFaces_[patchi][nRotatingFaces[patchi]++] = patchFacei;
            }
        }
    }

    forAll(includedFaces_, patchi)
    {
        if (includedFaces_[patchi].size() > 0)
        {
            rotatingFaces_[patchi] = includedFaces_[patchi];
        }
    }

    if (debug)
    {
        faceSet internalFaces(mesh_, "internalFaces", internalFaces_);
        Pout<< "Writing " << internalFaces.size()
            << " internal faces in MRF zone to faceSet "
            << internalFaces.name() << endl;
        internalFaces.write();

        faceSet MRFFaces(mesh_, "includedFaces", 100);
        forAll(includedFaces_, patchi)
        {
            forAll(includedFaces_[patchi], i)
            {
                label patchFacei = includedFaces_[patchi][i];
                MRFFaces.insert(patches[patchi].start()+patchFacei);
            }
        }
        Pout<< "Writing " << MRFFaces.size()
            << " patch faces in MRF zone to faceSet "
            << MRFFaces.name() << endl;
        MRFFaces.write();

        faceSet excludedFaces(mesh_, "excludedFaces", 100);
        forAll(excludedFaces_, patchi)
        {
            forAll(excludedFaces_[patchi], i)
            {
                label patchFacei = excludedFaces_[patchi][i];
                excludedFaces.insert(patches[patchi].start()+patchFacei);
            }
        }
        Pout<< "Writing " << excludedFaces.size()
            << " faces in MRF zone with special handling to faceSet "
            << excludedFaces.name() << endl;
        excludedFaces.write();

        faceSet rotatingFaces(mesh_, "rotatingFaces", 100);
        forAll(rotatingFaces_, patchi)
        {
            forAll(rotatingFaces_[patchi], i)
            {
                label patchFacei = rotatingFaces_[patchi][i];
                rotatingFaces.insert(patches[patchi].start()+patchFacei);
            }
        }
        Pout<< "Writing " << rotatingFaces.size()
            << " rotating faces not in MRF zone to faceSet "
            << rotatingFaces.name() << endl;
        rotatingFaces.write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZone::MRFZone
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    mesh_(mesh),
    name_(name),
    coeffs_(dict),
    active_(true),
    initRelativeVelocity_(false),
    cellZoneName_(cellZoneName),
    cellZoneID_(-1),
    excludedPatchNames_(),
    rotatingPatchNames_(),
    origin_(Zero),
    axis_(Zero),
    omega_(nullptr),
    Ut_(nullptr)
{
    read(dict);

    setInitialVelocity(); // by Gill
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::MRFZone::Omega() const
{
    return omega_->value(mesh_.time().timeOutputValue())*axis_;
}


Foam::vector Foam::MRFZone::Ut() const
{
    return Ut_->value(mesh_.time().timeOutputValue());
}


void Foam::MRFZone::addCoriolis
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    vectorField& ddtUc = ddtU.primitiveFieldRef();
    const vectorField& Uc = U;

    const vector Omega = this->Omega();
    const vector Ut = this->Ut();

    forAll(cells, i)
    {
        label celli = cells[i];
        ddtUc[celli] += (Omega ^ (Uc[celli] - Ut));
    }
}


void Foam::MRFZone::addCoriolis(fvVectorMatrix& UEqn, const bool rhs) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    const vector Omega = this->Omega();
    const vector Ut = this->Ut();

    if (rhs)
    {
        forAll(cells, i)
        {
            label celli = cells[i];
            Usource[celli] += V[celli]*(Omega ^ (U[celli] - Ut));
        }
    }
    else
    {
        forAll(cells, i)
        {
            label celli = cells[i];
            Usource[celli] -= V[celli]*(Omega ^ (U[celli] - Ut));
        }
    }
}


void Foam::MRFZone::addCoriolis
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn,
    const bool rhs
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    const vector Omega = this->Omega();
    const vector Ut = this->Ut();

    if (rhs)
    {
        forAll(cells, i)
        {
            label celli = cells[i];
            Usource[celli] += V[celli]*rho[celli]*(Omega ^ (U[celli] - Ut));
        }
    }
    else
    {
        forAll(cells, i)
        {
            label celli = cells[i];
            Usource[celli] -= V[celli]*rho[celli]*(Omega ^ (U[celli] - Ut));
        }
    }
}


void Foam::MRFZone::makeRelative(volVectorField& U) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const volVectorField& C = mesh_.C();

    const vector Omega = this->Omega();
    const vector Ut = this->Ut();

    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll(cells, i)
    {
        label celli = cells[i];
        U[celli] -= (Ut + (Omega ^ (C[celli] - origin_)));
    }

    // Included patches

    volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];
            Ubf[patchi][patchFacei] = Zero;
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];

            const vector& Cb = C.boundaryField()[patchi][patchFacei];

            Ubf[patchi][patchFacei] -= (Ut + (Omega ^ (Cb - origin_)));
        }
    }
}


void Foam::MRFZone::makeRelative(surfaceScalarField& phi) const
{
    makeRelativeRhoFlux(geometricOneField(), phi);
}


void Foam::MRFZone::makeRelative(FieldField<fvsPatchField, scalar>& phi) const
{
    makeRelativeRhoFlux(oneFieldField(), phi);
}


void Foam::MRFZone::makeRelative(Field<scalar>& phi, const label patchi) const
{
    makeRelativeRhoFlux(oneField(), phi, patchi);
}


void Foam::MRFZone::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeRelativeRhoFlux(rho, phi);
}


void Foam::MRFZone::setInitialVelocity() const // by Gill
{
    if 
    (
        !active_ 
     || !initRelativeVelocity_ 
     || mesh_.time().timeIndex() != 0
     || cellZoneID_ == -1
    )
    {
        return;
    }

    IOobject UIO
    (
        "U",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (UIO.typeHeaderOk<volVectorField>(true))
    {
        const localIOdictionary dict(UIO, "volVectorField");

        const word initType(dict.lookup("internalField"));

        const labelList& cells = mesh_.cellZones()[cellZoneID_];

        volVectorField& U
        (
            const_cast<volVectorField&>
            (
                mesh_.lookupObject<volVectorField>("U")
            )
        );

        surfaceScalarField& phi 
        (
            const_cast<surfaceScalarField&>
            (
                mesh_.lookupObject<surfaceScalarField>("phi")
            )
        );

        if (initType == "uniform" && cells.size() > 0)
        {
            const volVectorField& C = mesh_.C();

            const vector Omega = this->Omega();
            const vector Ut = this->Ut();

            forAll(cells, i)
            {
                label celli = cells[i];
                U[celli] += (Ut + (Omega ^ (C[celli] - origin_)));
            }
        }

        U.correctBoundaryConditions();
        correctBoundaryVelocity(U);

        if (phi.dimensions() == dimVelocity*dimArea)
        {
            phi = fvc::flux(U);
        }
        else if (mesh_.foundObject<volScalarField>("rho"))
        {
            const volScalarField& rho
            (
                mesh_.lookupObject<volScalarField>("rho")
            );

            phi = linearInterpolate(rho*U) & mesh_.Sf();
        }
    }
}


void Foam::MRFZone::makeAbsolute(volVectorField& U) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const volVectorField& C = mesh_.C();

    const vector Omega = this->Omega();
    const vector Ut = this->Ut();

    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll(cells, i)
    {
        label celli = cells[i];
        U[celli] += (Ut + (Omega ^ (C[celli] - origin_)));
    }

    // Included patches
    volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            const vector& Cb = C.boundaryField()[patchi][patchFacei];

            Ubf[patchi][patchFacei] = Ut + (Omega ^ (Cb - origin_));
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];

            const vector& Cb = C.boundaryField()[patchi][patchFacei];

            Ubf[patchi][patchFacei] += (Ut + (Omega ^ (Cb - origin_)));
        }
    }
}


void Foam::MRFZone::makeAbsolute(surfaceScalarField& phi) const
{
    makeAbsoluteRhoFlux(geometricOneField(), phi);
}


void Foam::MRFZone::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    makeAbsoluteRhoFlux(rho, phi);
}


void Foam::MRFZone::correctBoundaryVelocity(volVectorField& U) const
{
    if (!active_)
    {
        return;
    }

    const vector Omega = this->Omega();
    const vector Ut = this->Ut();

    // Included patches
    volVectorField::Boundary& Ubf = U.boundaryFieldRef();

    forAll(rotatingFaces_, patchi)
    {
        const vectorField& patchC = mesh_.Cf().boundaryField()[patchi];

        vectorField pfld(Ubf[patchi]);

        forAll(rotatingFaces_[patchi], i)
        {
            label patchFacei = rotatingFaces_[patchi][i];

            pfld[patchFacei] = Ut + (Omega ^ (patchC[patchFacei] - origin_));
        }

        Ubf[patchi] == pfld;
    }
}


void Foam::MRFZone::writeData(Ostream& os) const
{
    os  << nl;
    os.beginBlock(name_);

    os.writeEntry("active", active_);
    os.writeEntry("initRelativeVelocity", initRelativeVelocity_); // by Gill
    os.writeEntry("cellZone", cellZoneName_);
    os.writeEntry("origin", origin_);
    os.writeEntry("axis", axis_);
    omega_->writeData(os);
    Ut_->writeData(os);

    if (excludedPatchNames_.size())
    {
        os.writeEntry("nonRotatingPatches", excludedPatchNames_);
    }

    if (rotatingPatchNames_.size())
    {
        os.writeEntry("rotatingPatches", rotatingPatchNames_);
    }

    os.endBlock();
}


bool Foam::MRFZone::read(const dictionary& dict)
{
    coeffs_ = dict;

    coeffs_.readIfPresent("active", active_);

    { // by Gill
        coeffs_.readIfPresent("initRelativeVelocity", initRelativeVelocity_);
    }

    if (!active_)
    {
        cellZoneID_ = -1;
        return true;
    }

    coeffs_.readIfPresent("nonRotatingPatches", excludedPatchNames_);
    coeffs_.readIfPresent("rotatingPatches", rotatingPatchNames_);

    { // by Gill -- Need to replace nonRotatingPatches
        wordRes stationaryPatches;
        coeffs_.readIfPresent("stationaryPatches", stationaryPatches);
        excludedPatchNames_.push_back(stationaryPatches);
    }

    origin_ = coeffs_.get<vector>("origin");
    axis_ = coeffs_.get<vector>("axis").normalise();
    omega_.reset(Function1<scalar>::New("omega", coeffs_, &mesh_));
    Ut_.reset(new Function1Types::Constant<vector>("default", Zero));

    if (coeffs_.found("translationalVelocity"))
    {
        Ut_.reset
        (
            Function1<vector>::New
            (
                "translationalVelocity",
                coeffs_,
                &mesh_
            )
        );
    } // TODO: Need to make it mandatory // by Gill

    const word oldCellZoneName = cellZoneName_;
    if (cellZoneName_.empty())
    {
        coeffs_.readEntry("cellZone", cellZoneName_);
    }
    else
    {
        coeffs_.readIfPresent("cellZone", cellZoneName_);
    }

    if (cellZoneID_ == -1 || oldCellZoneName != cellZoneName_)
    {
        cellZoneID_ = mesh_.cellZones().findZoneID(cellZoneName_);

        excludedPatchLabels_ =
            mesh_.boundaryMesh().indices(excludedPatchNames_);

        if (!returnReduceOr(cellZoneID_ != -1))
        {
            FatalErrorInFunction
                << "cannot find MRF cellZone " << cellZoneName_
                << exit(FatalError);
        }

        setMRFFaces();
    }

    return true;
}


void Foam::MRFZone::update()
{
    if (mesh_.topoChanging())
    {
        setMRFFaces();
    }
}


// ************************************************************************* //

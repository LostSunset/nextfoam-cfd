/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "filmTemperatureCoupledMixedFvPatchScalarField.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const filmTemperatureCoupledMixedFvPatchScalarField::filmModelType&
filmTemperatureCoupledMixedFvPatchScalarField::filmModel() const
{
    typedef filmModelType ModelType;

    const UPtrList<const ModelType> models
    (
        db().time().csorted<ModelType>()
    );

    for (const ModelType& mdl : models)
    {
        if (mdl.regionMesh().name() == filmRegionName_)
        {
            return mdl;
        }
    }

    DynamicList<word> modelNames(models.size());
    for (const ModelType& mdl : models)
    {
        modelNames.push_back(mdl.regionMesh().name());
    }

    FatalErrorInFunction
        << "Unable to locate film region " << filmRegionName_
        << ".  Available regions include: " << modelNames
        << abort(FatalError);

    return models.front();  // Failure
}


void filmTemperatureCoupledMixedFvPatchScalarField::resetAMI
(
    const polyPatch& source,
    const polyPatch& target,
    const word& sourceMeshName
) const
{
    if (!AMIPtr_)
    {
        AMIPtr_.reset
        (
            AMIInterpolation::New
            (
                faceAreaWeightAMI::typeName,
                true,
                sourceMeshName == filmRegionName_
            )
        );

        AMIPtr_->calculate(source, target);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmTemperatureCoupledMixedFvPatchScalarField::
filmTemperatureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch()),  // default method (fluidThermo)
    filmRegionName_("surfaceFilmProperties"),
    coupledPatchName_(),
    TnbrName_("undefined-Tnbr"),
    AMIPtr_(nullptr)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 1.0;
}


filmTemperatureCoupledMixedFvPatchScalarField::
filmTemperatureCoupledMixedFvPatchScalarField
(
    const filmTemperatureCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    filmRegionName_(psf.filmRegionName_),
    coupledPatchName_(psf.coupledPatchName_),
    TnbrName_(psf.TnbrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    AMIPtr_(psf.AMIPtr_ ? psf.AMIPtr_->clone() : nullptr)
{}


filmTemperatureCoupledMixedFvPatchScalarField::
filmTemperatureCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    filmRegionName_
    (
        dict.getOrDefault<word>("filmRegion", "surfaceFilmProperties")
    ),
    coupledPatchName_(dict.get<word>("coupledPatch")),
    TnbrName_(dict.lookup("Tnbr")),
    AMIPtr_(nullptr)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    // Read list of layers
    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        dict.readEntry("kappaLayers", kappaLayers_);
    }

    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


filmTemperatureCoupledMixedFvPatchScalarField::
filmTemperatureCoupledMixedFvPatchScalarField
(
    const filmTemperatureCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    filmRegionName_(psf.filmRegionName_),
    coupledPatchName_(psf.coupledPatchName_),
    TnbrName_(psf.TnbrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    AMIPtr_(psf.AMIPtr_ ? psf.AMIPtr_->clone() : nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void filmTemperatureCoupledMixedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);
}


void filmTemperatureCoupledMixedFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const auto& fiptf = refCast<const myType>(ptf);

    temperatureCoupledBase::rmap(fiptf, addr);
}


void filmTemperatureCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const filmModelType& film = filmModel();

    const polyMesh& currentMesh = patch().boundaryMesh().mesh();
    const fvMesh& filmPrimaryMesh = film.primaryMesh();

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    // Get the coupling information from the mappedPatchBase
    const auto& mpp = refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();

    scalarField TP;
    scalarField KDeltaP;

    scalarField TN;
    scalarField KDeltaN;

    scalarField Tf;
    scalarField hf;
//    scalarField alphaf;

    if (currentMesh.name() == filmRegionName_)
    {
        const label filmPatchi = patch().index();

//        const label ownPatchi = mpp.samplePolyPatch().index();

//        const fvPatch& ownPatch = filmPrimaryMesh.boundary()[ownPatchi];

//        const auto& filmPrimaryTurb =
//            filmPrimaryMesh.lookupObject<compressible::turbulenceModel>
//            (
//                turbulenceModel::propertiesName
//            );
//
//        const scalarField filmPrimaryKappa = 
//            filmPrimaryTurb.kappaEff()().boundaryField()[ownPatchi];
        
        const label nbrPatchi =
            nbrMesh.boundaryMesh().findPatchID(coupledPatchName_);

        const fvPatch& nbrPatch =
            refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchi];

        const auto& nbrField =
            refCast<const myType>
            (
                nbrPatch.lookupPatchField<volScalarField>(TnbrName_)
            );

        resetAMI(patch().patch(), nbrPatch.patch(), currentMesh.name());

        // TP
        TP = film.TPrimary().boundaryField()[filmPatchi].patchInternalField();
        mpp.distribute(TP);

        // KDeltaP
//        KDeltaP = filmPrimaryKappa*ownPatch.deltaCoeffs();
//        mpp.distribute(KDeltaP);

        // TN
        TN = AMIPtr_->interpolateToSource(nbrField.patchInternalField());

        // KDeltaN
        KDeltaN = 
            AMIPtr_->interpolateToSource
            (
                nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs()
            );

        // Tf
        Tf = patchInternalField();

        // hf
        hf = film.htcw().h()().boundaryField()[filmPatchi];

        // alphaf
//        alphaf = film.alpha().boundaryField()[filmPatchi];

    }
    else if (currentMesh.name() == filmPrimaryMesh.name())
    {
        const label filmPatchi = mpp.samplePolyPatch().index();
            //nbrMesh.boundaryMesh().findPatchID(coupledPatchName_);

        const fvPatch& filmPatch =
            refCast<const fvMesh>(nbrMesh).boundary()[filmPatchi];

        const auto& filmField =
            refCast<const myType>
            (
                filmPatch.lookupPatchField<volScalarField>(TnbrName_)
            );

        const label nbrPatchi = film.nbrCoupledPatchID(film, filmPatchi);
        const fvPatch& nbrPatch = filmPrimaryMesh.boundary()[nbrPatchi];

        resetAMI(patch().patch(), nbrPatch.patch(), currentMesh.name());

        // TP
        TP = patchInternalField();

        // KDeltaP
        KDeltaP = this->kappa(*this)*patch().deltaCoeffs();

        // TN
        TN = film.TPrimary().boundaryField()[filmPatchi].patchInternalField();
        mpp.distribute(TN);

//        const auto& filmPrimaryTurb =
//            filmPrimaryMesh.lookupObject<compressible::turbulenceModel>
//            (
//                turbulenceModel::propertiesName
//            );
//
//        const scalarField filmPrimaryKappa = 
//            filmPrimaryTurb.kappaEff()().boundaryField()[nbrPatchi];
//
//        // KDeltaN
//        KDeltaN = 
//            AMIPtr_->interpolateToSource
//            (
//                filmPrimaryKappa*nbrPatch.deltaCoeffs()
//            );

        // Tf
        Tf = filmField.patchInternalField();
        mpp.distribute(Tf);

        // hf
        hf = film.htcw().h()().boundaryField()[filmPatchi];
        mpp.distribute(hf);

        // alphaf
//        alphaf = film.alpha().boundaryField()[filmPatchi];
//        mpp.distribute(alphaf);
    }
    else
    {
        FatalErrorInFunction
            << type() << " condition is intended to be applied to either the "
            << "primary or film regions only"
            << exit(FatalError);
    }

    scalarField rKDeltaC(this->size(), Zero);
    if (thicknessLayers_.size())
    {
        forAll(thicknessLayers_, iLayer)
        {
            rKDeltaC += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
        }
    }

    if (currentMesh.name() == filmRegionName_)
    {
        scalarField F1(KDeltaN);
        scalarField F2((1 + rKDeltaC*KDeltaN)*hf);

        // valueFraction() must be set to 1.0 because thermoSingleLayer 
        // enforces the hs boundary type for this patch as fixedValue.
        valueFraction() = 1.0;
        refValue() = (F1*TN + F2*Tf)/(F1 + F2);
    }
    else
    {
        scalarField F1(hf);
        scalarField F2((1.0 + hf*rKDeltaC)*KDeltaP);

        valueFraction() = F1/(F1 + F2);
        refValue() = Tf; 
    }

// Uncomment these lines if a suitable method is available to evaluate kappaEff 
// in the primary region for the case where alphaf = 0 (i.e. fully dry film).
//    if (currentMesh.name() == filmRegionName_)
//    {
//        scalarField F1 = alphaf*KDeltaN*KDeltaP + (1 - alphaf)*hf*KDeltaN;
//        scalarField F2 = (1 + rKDeltaC*KDeltaN)*hf*KDeltaP;
//        scalarField F3 = F1 + F2;
//
//        // valueFraction() must be set to 1.0 because thermoSingleLayer 
//        // enforces the hs boundary type for this patch as fixedValue.
//        valueFraction() = 1.0;
//        refValue() = (F1*TN + alphaf*F2*Tf + (1 - alphaf)*F2*TP)/F3;
//    }
//    else
//    {
//        scalarField F1 = hf*KDeltaN;
//        scalarField F2 = 
//            (alphaf + hf*rKDeltaC)*KDeltaN*KDeltaP 
//          + (1 - alphaf)*hf*KDeltaP;
//
//        valueFraction() = F1/(F1 + F2);
//        refValue() = (1 - alphaf)*TN + alphaf*Tf; 
//    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        //scalar Qt = gSum(qConv*patch().magSf());

        //Info<< mesh.name() << ':'
        //    << patch().name() << ':'
        //    << this->internalField().name() << " <- "
        //    << nbrMesh.name() << ':'
         //   << nbrPatch.name() << ':'
         //   << this->internalField().name() << " :" << nl
         //   << "     total heat     [W] : " << Qt << nl
         //   << "     wall temperature "
         //   << " min:" << gMin(*this)
         //   << " max:" << gMax(*this)
         //   << " avg:" << gAverage(*this)
         //   << endl;
    }

    UPstream::msgType(oldTag);  // Restore tag
}


void filmTemperatureCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>
    (
        "filmRegion",
        "surfaceFilmProperties",
        filmRegionName_
    );
    os.writeEntry("coupledPatch", coupledPatchName_);
    os.writeEntry("Tnbr", TnbrName_);

    if (thicknessLayers_.size())
    {
        thicknessLayers_.writeEntry("thicknessLayers", os);
        kappaLayers_.writeEntry("kappaLayers", os);
    }

    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    filmTemperatureCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //

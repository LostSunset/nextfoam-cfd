/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "heRhoThermo.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::calculate
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& he,
    volScalarField& psi,
    volScalarField& rho,
    volScalarField& mu,
    volScalarField& alpha,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculate
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            psi.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            alpha.oldTime(),
            true
        );
    }

    const scalarField& hCells = he.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& TCells = T.primitiveFieldRef();
    scalarField& psiCells = psi.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();
    scalarField& alphaCells = alpha.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        if (this->updateT())
        {
            TCells[celli] = mixture_.THE
            (
                hCells[celli],
                pCells[celli],
                TCells[celli]
            );
        }

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    volScalarField::Boundary& heBf = he.boundaryFieldRef();
    volScalarField::Boundary& muBf = mu.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = alpha.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                if (this->updateT())
                {
                    pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);
                }

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::calculateTemperature
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& he,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculateTemperature
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            true
        );
    }

    const scalarField& hCells = he.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& TCells = T.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        if (this->updateT())
        {
            TCells[celli] = mixture_.THE
            (
                hCells[celli],
                pCells[celli],
                TCells[celli]
            );
        }
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& heBf = he.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                if (this->updateT())
                {
                    pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);
                }
            }
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::calculateProperties
(
    const volScalarField& p,
    const volScalarField& T,
    volScalarField& psi,
    volScalarField& rho,
    volScalarField& mu,
    volScalarField& alpha,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculateProperties
        (
            p.oldTime(),
            T.oldTime(),
            psi.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            alpha.oldTime(),
            true
        );
    }

    const scalarField& pCells = p.primitiveField();
    const scalarField& TCells = T.primitiveField();

    scalarField& psiCells = psi.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();
    scalarField& alphaCells = alpha.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);
        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    const volScalarField::Boundary& TBf = T.boundaryField();

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    volScalarField::Boundary& muBf = mu.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = alpha.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        const fvPatchScalarField& pT = TBf[patchi];

        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        forAll(pT, facei)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->patchFaceMixture(patchi, facei);

            ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
            prho[facei] = mixture_.rho(pp[facei], pT[facei]);
            pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
            palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::calculateRhoMu
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& rho,
    volScalarField& mu,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculateRhoMu
        (
            p.oldTime(),
            T.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            true
        );
    }

    const scalarField& pCells = p.primitiveField();

    scalarField& TCells = T.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);
        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    volScalarField::Boundary& muBf = mu.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];

        forAll(pT, facei)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->patchFaceMixture(patchi, facei);

            prho[facei] = mixture_.rho(pp[facei], pT[facei]);
            pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::calculateInternal
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& he,
    volScalarField& psi,
    volScalarField& rho,
    volScalarField& mu,
    volScalarField& alpha,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculateInternal
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            psi.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            alpha.oldTime(),
            true
        );
    }

    const scalarField& hCells = he.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& TCells = T.primitiveFieldRef();
    scalarField& psiCells = psi.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();
    scalarField& alphaCells = alpha.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        if (this->updateT())
        {
            TCells[celli] = mixture_.THE
            (
                hCells[celli],
                pCells[celli],
                TCells[celli]
            );
        }

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);

        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::calculateBoundary
(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& he,
    volScalarField& psi,
    volScalarField& rho,
    volScalarField& mu,
    volScalarField& alpha,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if T.oldTime() is
    // created from T, it starts from the unconverted T
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        calculateBoundary
        (
            p.oldTime(),
            T.oldTime(),
            he.oldTime(),
            psi.oldTime(),
            rho.oldTime(),
            mu.oldTime(),
            alpha.oldTime(),
            true
        );
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    volScalarField::Boundary& TBf = T.boundaryFieldRef();
    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
    volScalarField::Boundary& heBf = he.boundaryFieldRef();
    volScalarField::Boundary& muBf = mu.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = alpha.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                phe[facei] = mixture_.HE(pp[facei], pT[facei]);

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                if (this->updateT())
                {
                    pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);
                }

                ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
                prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::addHydrostaticPressure
(
    volScalarField& p,
    const volScalarField& rho
)
{
    if
    (   
        p.mesh().time().foundObject<uniformDimensionedVectorField>("g")
     && p.mesh().foundObject<IOdictionary>("operatingConditions")
    )
    {
        IOobject pIO
        (
            "p",
            p.mesh().time().timeName(),
            p.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!pIO.typeHeaderOk<volScalarField>(true))
        {
            const volScalarField& gh
            (
                p.mesh().lookupObject<volScalarField>("gh")
            );

            p += rho*gh;
        }
    }

} // by Gill

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoThermo<BasicPsiThermo, MixtureType>::heRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );

    //addHydrostaticPressure      // by Gill
    //(
    //    this->p_, 
    //    this->rho_
    //);
}


template<class BasicPsiThermo, class MixtureType>
Foam::heRhoThermo<BasicPsiThermo, MixtureType>::heRhoThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName, dictName)
{
    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );

    //addHydrostaticPressure      // by Gill
    //(
    //    this->p_, 
    //    this->rho_
    //);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoThermo<BasicPsiThermo, MixtureType>::~heRhoThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::correct()
{
    DebugInFunction << endl;

    calculate
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::limitedCorrect
(
    const scalar maxDeltaT,
    const scalar localCorrection,
    const scalar minTemp
)
{
    DebugInFunction << endl;

    // Use the current T as the reference state for limiting
    // the change implied by the newly solved he.
    const scalarField& TCells = this->T_.primitiveField();
    const scalarField& pCells = this->p_.primitiveField();

    scalarField& hCells = this->he_.primitiveFieldRef();

    label nLimited = 0;
    scalar maxMagTemperatureCorrection = 0;

    // Limit the enthalpy correction so that the implied temperature change
    // remains within a fraction of both the reference temperature scale and
    // the local temperature scale.
    forAll(hCells, celli)
    {
        const scalar T0 = TCells[celli];

        const scalar localMaxDeltaT =
            min(maxDeltaT, localCorrection*max(mag(T0), scalar(1)));
            //max(maxDeltaT, localCorrection*max(mag(T0), scalar(1)));

        const scalar TminAllowed = max(T0 - localMaxDeltaT, minTemp);
        const scalar TmaxAllowed = T0 + localMaxDeltaT;

        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        const scalar hMin = mixture_.HE(pCells[celli], TminAllowed);
        const scalar hMax = mixture_.HE(pCells[celli], TmaxAllowed);

        const scalar hRaw = hCells[celli];
        const scalar hLimited = min(max(hRaw, hMin), hMax);

        if (mag(hLimited - hRaw) > SMALL)
        {
            hCells[celli] = hLimited;

            ++nLimited;

            // Diagnostic only: estimate the implied temperature correction
            // using the corresponding allowed bound.
            const scalar Tlimited =
                (hRaw < hMin ? TminAllowed : TmaxAllowed);

            maxMagTemperatureCorrection =
                max(maxMagTemperatureCorrection, mag(Tlimited - T0));
        }
    }

    reduce(nLimited, sumOp<label>());
    reduce(maxMagTemperatureCorrection, maxOp<scalar>());

    if (nLimited)
    {
        // Let he boundary conditions decide patch values.
        this->he_.correctBoundaryConditions();

        Info<< "Temperature limited " << nLimited
            << " cells, max allowed temperature change = "
            << maxMagTemperatureCorrection << " K"
            << nl << endl;
    }

    // Now calculate T from the accepted he and update properties.
    this->correct();

    DebugInFunction << "Finished" << endl;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::correctTemperature()
{
    DebugInFunction << endl;

    calculateTemperature
    (
        this->p_,
        this->T_,
        this->he_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::correctProperties()
{
    DebugInFunction << endl;

    calculateProperties
    (
        this->p_,
        this->T_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::correctRhoMu()
{
    DebugInFunction << endl;

    calculateRhoMu
    (
        this->p_,
        this->T_,
        this->rho_,
        this->mu_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::correctInternal()
{
    DebugInFunction << endl;

    calculateInternal
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}


template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoThermo<BasicPsiThermo, MixtureType>::correctBoundary()
{
    DebugInFunction << endl;

    calculateBoundary
    (
        this->p_,
        this->T_,
        this->he_,
        this->psi_,
        this->rho_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}


// ************************************************************************* //

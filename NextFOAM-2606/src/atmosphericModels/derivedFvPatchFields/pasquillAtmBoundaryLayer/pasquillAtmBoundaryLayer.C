/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "pasquillAtmBoundaryLayer.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Enum<typename pasquillAtmBoundaryLayer::PasquillClass>
pasquillAtmBoundaryLayer::PasquillClassNames_
({
    { PasquillClass::A, "extremelyUnstable" },
    { PasquillClass::B, "moderatelyUnstable" },
    { PasquillClass::C, "slightlyUnstable" },
    { PasquillClass::D, "neutral" },
    { PasquillClass::E, "slightlyStable" },
    { PasquillClass::F, "stable" },
    { PasquillClass::CLASS_A, "A" },
    { PasquillClass::CLASS_B, "B" },
    { PasquillClass::CLASS_C, "C" },
    { PasquillClass::CLASS_D, "D" },
    { PasquillClass::CLASS_E, "E" },
    { PasquillClass::CLASS_F, "F" }
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pasquillAtmBoundaryLayer::pasquillAtmBoundaryLayer
(
    const Time& time,
    const polyPatch& pp
)
:
    initABL_(false),
    stability_(PasquillClass::D),
    kappa_(0.41),
    Cmu_(0.09),
    ppMin_((boundBox(pp.localPoints())).min()),
    time_(time),
    patch_(pp),
    flowDir_(nullptr),
    zDir_(nullptr),
    Uref_(nullptr),
    Zref_(nullptr),
    z0_(nullptr),
    d_(nullptr),
    latitude_(37),
    qdot_s_(0.0),
    Cp_(1006.0),
    rho_(1.225),
    Tref_(300.0)
{}


pasquillAtmBoundaryLayer::pasquillAtmBoundaryLayer
(
    const Time& time,
    const polyPatch& pp,
    const dictionary& dict
)
:
    initABL_(dict.getOrDefault<bool>("initABL", true)),
    stability_
    (
        PasquillClassNames_.getOrDefault
        (
            "stability", 
            dict, 
            PasquillClass::D
        )
    ),
    kappa_
    (
        dict.getCheckOrDefault<scalar>("kappa", 0.41, scalarMinMax::ge(SMALL))
    ),
    Cmu_(dict.getCheckOrDefault<scalar>("Cmu", 0.09, scalarMinMax::ge(SMALL))),
    ppMin_((boundBox(pp.localPoints())).min()),
    time_(time),
    patch_(pp),
    flowDir_(Function1<vector>::New("flowDir", dict, &time)),
    zDir_(Function1<vector>::New("zDir", dict, &time)),
    Uref_(Function1<scalar>::New("Uref", dict, &time)),
    Zref_(Function1<scalar>::New("Zref", dict, &time)),
    z0_(PatchFunction1<scalar>::New(pp, "z0", dict)),
    d_(PatchFunction1<scalar>::New(pp, "d", dict)),
    latitude_(dict.get<scalar>("latitude")),
    qdot_s_(dict.getOrDefault<scalar>("surfaceHeatFlux", 0.0)),
    Cp_(dict.getOrDefault<scalar>("Cp", 1006)),
    rho_(dict.getOrDefault<scalar>("rho", 1.225)),
    Tref_(dict.getOrDefault<scalar>("Tref", 300))
{}


pasquillAtmBoundaryLayer::pasquillAtmBoundaryLayer
(
    const pasquillAtmBoundaryLayer& abl,
    const fvPatch& patch,
    const fvPatchFieldMapper& mapper
)
:
    initABL_(abl.initABL_),
    stability_(abl.stability_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    ppMin_(abl.ppMin_),
    time_(abl.time_),
    patch_(patch.patch()),
    flowDir_(abl.flowDir_.clone()),
    zDir_(abl.zDir_.clone()),
    Uref_(abl.Uref_.clone()),
    Zref_(abl.Zref_.clone()),
    z0_(abl.z0_.clone(patch_)),
    d_(abl.d_.clone(patch_)),
    latitude_(abl.latitude_),
    qdot_s_(abl.qdot_s_),
    Cp_(abl.Cp_),
    rho_(abl.rho_),
    Tref_(abl.Tref_)
{}


pasquillAtmBoundaryLayer::pasquillAtmBoundaryLayer
(
    const pasquillAtmBoundaryLayer& abl
)
:
    initABL_(abl.initABL_),
    stability_(abl.stability_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    ppMin_(abl.ppMin_),
    time_(abl.time_),
    patch_(abl.patch_),
    flowDir_(abl.flowDir_.clone()),
    zDir_(abl.zDir_.clone()),
    Uref_(abl.Uref_.clone()),
    Zref_(abl.Zref_.clone()),
    z0_(abl.z0_.clone(patch_)),
    d_(abl.d_.clone(patch_)),
    latitude_(abl.latitude_),
    qdot_s_(abl.qdot_s_),
    Cp_(abl.Cp_),
    rho_(abl.rho_),
    Tref_(abl.Tref_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<vectorField> pasquillAtmBoundaryLayer::flowDir(const fvPatch& fvp) const
{
    const scalar t = time_.timeOutputValue();
    const vector dir(flowDir_->value(t));
    const scalar magDir = mag(dir);

    // If magnitude of flowDir is 0, face-normal with zero z-dir component
    // is used    
    if (magDir < SMALL)
    {
        tmp<vectorField> flowDir = -zDir() ^ fvp.nf() ^ zDir();

        return flowDir()/mag(flowDir());
    }

    return tmp<vectorField>::New(patch_.size(), dir/magDir);
}


vector pasquillAtmBoundaryLayer::zDir() const
{
    const scalar t = time_.timeOutputValue();
    const vector dir(zDir_->value(t));
    const scalar magDir = mag(dir);

    if (magDir < SMALL)
    {
        FatalErrorInFunction
            << "magnitude of " << zDir_->name() << " = " << magDir
            << " vector must be greater than zero"
            << abort(FatalError);
    }

    return dir/magDir;
}


tmp<scalarField> pasquillAtmBoundaryLayer::MOL(const scalarField& z0) const
{
    scalar Zs = 0;
    scalar Ls = 1;

    switch (stability_)
    {
        case PasquillClass::A:
        case PasquillClass::CLASS_A:
        {
            Ls = 33.162;
            Zs = 1117;

            break;
        }
        case PasquillClass::B:
        case PasquillClass::CLASS_B:
        {
            Ls = 33.258;
            Zs = 11.46;

            break;
        }
        case PasquillClass::C:
        case PasquillClass::CLASS_C:
        {
            Ls = 51.787;
            Zs = 1.324;

            break;
        }
        case PasquillClass::D:
        case PasquillClass::CLASS_D:
        {
            break;
        }
        case PasquillClass::E:
        case PasquillClass::CLASS_E:
        {
            Ls = -48.33;
            Zs = 1.262;

            break;
        }
        case PasquillClass::F:
        case PasquillClass::CLASS_F:
        {
            Ls = -31.323;
            Zs = 19.36;

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration "
                << PasquillClassNames_[stability_]
                << abort(FatalError);
        }
    }

    if (Zs > SMALL)
    {
        return Ls/log10(z0/Zs);
    }

    return tmp<scalarField>::New(patch_.size(), 1.0e6);
}


tmp<scalarField> pasquillAtmBoundaryLayer::Ustar
(
    const scalarField& z0,
    const scalarField& d,
    const scalarField& L
) const
{
    const scalar t = time_.timeOutputValue();
    const scalar Uref = Uref_->value(t);
    const scalar Zref = Zref_->value(t);

    if (Zref < 0)
    {
        FatalErrorInFunction
            << "Negative entry in " << Zref_->name() << " = " << Zref
            << abort(FatalError);
    }

    tmp<scalarField> zRefField(tmp<scalarField>::New(patch_.size(), Zref));

    return kappa_*Uref/(log((Zref - d + z0)/z0) - PhiU(zRefField, L));
}


tmp<scalarField> pasquillAtmBoundaryLayer::H
(
    const scalarField& uStar,
    const scalarField& L
) const
{
    tmp<scalarField> th(tmp<scalarField>::New(patch_.size(), Zero));
    scalarField& h = th.ref();

    forAll(L, i)
    {
        if (mag(L[i]) >= 1.0e5)
        {
            scalar fcoeff = 0.3/(2*7.2921*1e-5*sin(latitude_*pi/180.0));

            h[i] = min(500, fcoeff*uStar[i]);
        }
        else if (L[i] < 0) 
        { 
            if (L[i] <= -100)
            {
                h[i] = 1000;
            }
            else
            {
                h[i] = 1500;
            }
        }
        else
        {
            scalar fcoeff = L[i]/(2*7.2921*1e-5*sin(latitude_*pi/180.0));

            h[i] = 0.4*sqrt(fcoeff*uStar[i]);
        }
    }

    return th;
}


tmp<scalarField> pasquillAtmBoundaryLayer::PhiU
(
    const scalarField& z,
    const scalarField& L
) const
{
    tmp<scalarField> tPhi(tmp<scalarField>::New(patch_.size(), Zero));
    scalarField& Phi = tPhi.ref();


    forAll(L, i)
    {
        if (mag(L[i]) >= 1e5)
        {
            Phi[i] = 0.0;
        }
        else if (L[i] < 0)
        {
            const scalar X = pow(1.0 - 16*z[i]/L[i], 0.25);

            Phi[i] = 
                2*log(0.5*(1 + X)) 
              + log(0.5*(1 + sqr(X))) 
              - 2*atan(X) 
              + 0.5*pi;
        }
        else
        {
            Phi[i] = -17*(1.0 - exp(-0.29*z[i]/L[i]));
        }
    }

    return tPhi;
}


void pasquillAtmBoundaryLayer::autoMap(const fvPatchFieldMapper& mapper)
{}


void pasquillAtmBoundaryLayer::rmap
(
    const pasquillAtmBoundaryLayer& abl,
    const labelList& addr
)
{}


tmp<vectorField> pasquillAtmBoundaryLayer::U(const fvPatch& fvp) const
{
    const vectorField& pCf = fvp.Cf();
    const scalar t = time_.timeOutputValue();
    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));   
    const scalar groundMin = zDir() & ppMin_;

    tmp<scalarField> tL = MOL(z0);
    const scalarField& L = tL();

    tmp<scalarField> tuStar = Ustar(z0, d, L);
    const scalarField& uStar = tuStar();

    tmp<scalarField> tz = min((zDir() & pCf) - groundMin, H(uStar, L));
    const scalarField& z = tz();

    scalarField Un
    (
        uStar*(log((z - d + z0)/z0) - PhiU(z, L))/kappa_
    );

    return flowDir(fvp)*Un;
}


tmp<scalarField> pasquillAtmBoundaryLayer::k(const fvPatch& fvp) const
{
    const vectorField& pCf = fvp.Cf();
    const scalar t = time_.timeOutputValue();
    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));   
    const scalar groundMin = zDir() & ppMin_;

    tmp<scalarField> tL = MOL(z0);
    const scalarField& L = tL();

    tmp<scalarField> tuStar = Ustar(z0, d, L);
    const scalarField& uStar = tuStar();

    tmp<scalarField> tz = (zDir() & pCf) - groundMin;
    const scalarField& z = tz();

    tmp<scalarField> th = H(uStar, L);
    const scalarField& h = th();

    tmp<scalarField> tzbyh = min(z/h, 1.0);
    const scalarField& zbyh = tzbyh();

    tmp<scalarField> tk(tmp<scalarField>::New(patch_.size(), Zero));
    scalarField& k = tk.ref();

    forAll(L, i)
    {
        if ((L[i] > 0.0) || (L[i] < -1.0e5)) 
        {
            if (zbyh[i] <= 0.1)
            {
                k[i] = 6.0*sqr(uStar[i]);
            }
            else
            {
                k[i] = 6.0*sqr(uStar[i])*pow(1 - zbyh[i], 1.75);
            }
        }
        else
        {
            const scalar g = 9.81;
            const scalar wStar(pow(g*qdot_s_*h[i]/Tref_/rho_/Cp_, 1/3));

            if (zbyh[i] <= 0.1)
            {
                k[i] = 
                    0.36*sqr(wStar)
                  + 0.85*sqr(uStar[i])*pow(1.0 - 3*z[i]/L[i], 2/3);
            }
            else
            {
                k[i] =
                    (0.36 + 0.9*pow(zbyh[i], 2/3)*sqr(1 - 0.8*zbyh[i]))
                    *sqr(wStar);
            }
        }
    }

    return tk;
}


tmp<scalarField> pasquillAtmBoundaryLayer::epsilon(const fvPatch& fvp) const
{
    const vectorField& pCf = fvp.Cf();
    const scalar t = time_.timeOutputValue();
    const scalarField d(d_->value(t));
    const scalarField z0(max(z0_->value(t), ROOTVSMALL));   
    const scalar groundMin = zDir() & ppMin_;

    tmp<scalarField> tL = MOL(z0);
    const scalarField& L = tL();

    tmp<scalarField> tuStar = Ustar(z0, d, L);
    const scalarField& uStar = tuStar();

    tmp<scalarField> tz = (zDir() & pCf) - groundMin;
    const scalarField& z = tz();

    tmp<scalarField> tzbyh = min(z/H(uStar, L), 1.0);
    const scalarField& zbyh = tzbyh();

    const scalar g = 9.81;
    const scalarField cubeUstar(pow(uStar, 3.0));
    const scalar cubeWStarbyh(g*qdot_s_/Tref_/rho_/Cp_);

    tmp<scalarField> tEpsilon(tmp<scalarField>::New(patch_.size(), Zero));
    scalarField& epsilon = tEpsilon.ref();

    forAll(L, i)
    {
        if ((L[i] > 0.0) || (L[i] < -1.0e5)) 
        {
            if (zbyh[i] <= 0.1)
            {
                epsilon[i] = cubeUstar[i]*(1.24 + 4.3*z[i]/L[i])/kappa_/z[i];
            }
            else
            {
                epsilon[i] = 
                    cubeUstar[i]*(1.24 + 4.3*z[i]/L[i])
                    *pow(1 - 0.85*zbyh[i], 1.5)/kappa_/z[i];
            }
        }
        else 
        {
            if (zbyh[i] <= 0.1)
            {
                epsilon[i] = 
                    cubeUstar[i]*pow(1.0 + 0.5*pow(mag(z[i]/L[i]), 2/3), 1.5)
                    /kappa_/z[i];
            }
            else
            {
                epsilon[i] = cubeWStarbyh*(0.8 - 0.3*zbyh[i]);
            }
        }
    }

    return tEpsilon;
}


tmp<scalarField> pasquillAtmBoundaryLayer::omega(const fvPatch& fvp) const
{
    // Currently, "omega = epsilon/k/Cmu" is being used.
    // need to verify whether "omega = (1/sqrt(Cmu))*(du/dz)" can be used or
    // confirm if any other suitable formulation exists.

    return epsilon(fvp)/max(k(fvp), SMALL)/Cmu_;
}


void pasquillAtmBoundaryLayer::write(Ostream& os) const
{
    os.writeEntry("initABL", initABL_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("Cmu", Cmu_);
    if (flowDir_)
    {
        flowDir_->writeData(os);
    }
    if (zDir_)
    {
        zDir_->writeData(os);
    }
    if (Uref_)
    {
        Uref_->writeData(os);
    }
    if (Zref_)
    {
        Zref_->writeData(os);
    }
    if (z0_)
    {
        z0_->writeData(os) ;
    }
    if (d_)
    {
        d_->writeData(os);
    }
    os.writeEntry("stability", PasquillClassNames_[stability_]);
    os.writeEntry("latitude", latitude_);

    os.writeEntry("surfaceHeatFlux", qdot_s_);
    os.writeEntry("Cp", Cp_);
    os.writeEntry("rho", rho_);
    os.writeEntry("Tref", Tref_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

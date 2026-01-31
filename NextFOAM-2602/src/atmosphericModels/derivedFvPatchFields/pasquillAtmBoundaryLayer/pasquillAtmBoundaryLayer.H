/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2018 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::pasquillAtmBoundaryLayer

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    Base class to set log-law type ground-normal inlet boundary conditions for
    wind velocity and turbulence quantities for homogeneous, two-dimensional,
    dry-air, equilibrium and neutral atmospheric boundary layer (ABL) modelling.

    The ground-normal profile expressions are due to \c YGCJ
    (refer to references below) whereat \c RH expressions were generalised:

    \f[
        u = \frac{u^*}{\kappa} \ln \left( \frac{z - d + z_0}{z_0} - \phi_u(z) \right)
    \f]

    \f[
        v = w = 0
    \f]
    
    \f[
        u^* =
            \frac{u_{ref} \kappa}{\ln\left(\frac{z_{ref} - d + z_0}{z_0} - \phi_u(z_{ref}) \right)}
    \f]

    \f[
        \phi_u(z) = \left\{\begin{matrix} 
                    0 & Neutral \\
                    2 \ln \left(\ frac{1 + \xi}{2} \right) + \ln \left( \frac{1 + \xi ^ 2}{2} \right)
                     - 2 \arctan (\xi) + \frac{\pi}{2} & L < 0 \\
                    -17 \left( 1 - exp(-0.29 \frac{z}{L}) \right) & L > 0
                    \end{matrix}\right.
    \f] 

    \f[
        \xi = (1 - 16 z / L) ^ {1/4}
    \f]

    \f[
        k = \frac{(u^*)^2}{\sqrt{C_\mu}}
        \sqrt{C_1 \ln \left( \frac{z - d + z_0}{z_0} \right) + C_2}
    \f]

    \f[
        \frac{1}{L} = \frac{1}{L_s} \log \frac{z_0}{z_s}
    \f]

    for A, B, C

    \f[
        k(z) = \left\{\begin{matrix} 
            0.36w^{*2} + 0.85u^{*2} (1 - 3 \frac{z}{L})^{2/3}  & z \leq 0.1h \\
            \left( 0.36 + 0.9 (\frac{z}{h})^{2/3} (1 - 0.8 \frac{z}{h})^2 \right) w^{*2}  & z > 0.1h
            \end{matrix}\right.
    \f]

    \f[
        \epsilon(z) = \left\{\begin{matrix} 
            \frac{u^{*3}}{kz} \left( 0.5 \left | \frac{z}{L} \right | ^ {2/3} \right) ^ {3/2}  & z \leq 0.1h \\
            \frac{w^{*3}}{h} (0.8 - 0.3 \frac{z}{h})  & z > 0.1h
            \end{matrix}\right.
    \f]

    \f[
        w^* = \left( \frac{g \dot{q}_s h}{T_0 \rho c_p} \right) ^ {1/3}
    \f]

    for D, E, F

    \f[
        k(z) = \left\{\begin{matrix} 
            6 u^{*2}  & z \leq 0.1h \\
            6 u^{*2} ( 1 - \frac{z}{h} ) ^ {1.75}  & z > 0.1h
            \end{matrix}\right.
    \f]

    \f[
        \epsilon(z) = \left\{\begin{matrix} 
            \frac{u^{*3}}{kz} \left( 1.24 + 4.3 \frac{z}{L} \right) & z \leq 0.1h \\
            \frac{u^{*3}}{kz} \left( 1.24 + 4.3 \frac{z}{L} \right)(1 - 0.85 \frac{z}{h}) ^ {3/2}  & z > 0.1h
            \end{matrix}\right.
    \f]

    \f[
        \omega(z) = \frac{1}{C_{\mu}} \frac{\epsilon(z)}{k(z)}
    \f]

    where
    \vartable
      u        | Ground-normal streamwise flow speed profile          [m/s]
      v        | Spanwise flow speed                                  [m/s]
      w        | Ground-normal flow speed                             [m/s]
      k        | Ground-normal turbulent kinetic energy (TKE) profile [m^2/s^2]
      \epsilon | Ground-normal TKE dissipation rate profile           [m^2/s^3]
      \omega   | Ground-normal specific dissipation rate profile      [m^2/s^3]
      u^*      | Friction velocity                                    [m/s]
      \kappa   | von Kármán constant                                  [-]
      C_\mu    | Empirical model constant                             [-]
      z        | Ground-normal coordinate component                   [m]
      d        | Ground-normal displacement height                    [m]
      z_0      | Aerodynamic roughness length                         [m]
      u_{ref}  | Reference mean streamwise wind speed at \f$z_{ref}\f$ [m/s]
      z_{ref}  | Reference height being used in \f$u^*\f$ estimations [m]
      h        | Boundary layer height for Golder profiles by Han     [m]
      Ls       | Charicteristic length for Golder porifles by BW      [m]
      Zs       | Charicteristic height for Golder porifles by BW      [m]
      \phi_u   | 
      \xi      | 
      w^*      | 
    \endvartable

    Reference:
    \verbatim
        The ground-normal profile expressions (tag:RH):
            Richards, P. J., & Hoxey, R. P. (1993).
            Appropriate boundary conditions for computational wind
            engineering models using the k-ε turbulence model.
            In Computational Wind Engineering 1 (pp. 145-153).
            DOI:10.1016/B978-0-444-81688-7.50018-8

    \endverbatim

Usage
    Example of the entries provided for the inherited boundary conditions:
    \verbatim
    inlet
    {
        // Mandatory and other optional entries
        ...

        // Mandatory (inherited) entries (runtime modifiable)
        flowDir         (1 0 0);
        zDir            (0 0 1);
        Uref            10.0;
        Zref            0.0;
        z0              uniform 0.1;
        d               uniform 0.0;

        // Optional (inherited) entries (unmodifiable)
        kappa           0.41;
        Cmu             0.09;
        initABL         true;
        phi             phi;

        Zs              0.0;
        Ls              1.0;

        qdot            10.0;

        // Conditional mandatory (inherited) entries (runtime modifiable)
        value           uniform 0;    // when initABL=false
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description              | Type          | Reqd  | Deflt
      flowDir   | Flow direction           | Function1<vector> | yes | -
      zDir      | Ground-normal direction  | Function1<vector> | yes | -
      Uref      | Reference mean streamwise flow speed being used in <!--
             --> \f$u^*\f$ estimations [m/s] | Function1<scalar> | yes | -
      Zref      | Reference height being used in \f$u^*\f$ estimations [m]  <!--
                                       --> | Function1<scalar> | yes | -
      z0        | Surface roughness length [m] <!--
                                       --> | PatchFunction1<scalar> | yes | -
      d         | Displacement height [m] - see Notes <!--
                                       --> | PatchFunction1<scalar> | yes | -  

      kappa     | von Kármán constant      | scalar | no  | 0.41
      Cmu       | Empirical model constant | scalar | no  | 0.09
      initABL   | Flag to initialise profiles with the theoretical <!--
            --> ABL expressions, otherwise use "value" list        <!--
                                       --> | bool   | no  | true
      value     | ABL profile content when initABL=false <!--
                                       --> | scalarList | conditional | -
      phi       | Name of the flux field   | word         | no        | phi

      h         | Boundary layer height    | scalar | no  | min(500, 0.3 * u^* / f)
      Ls        | Characteristic length    | scalar | no  | 1.0
      Zs        | Characteristic height    | scalar | no  | 0.0
      qdot      | Mean surface heat flux   | scalar | no  | 0.0
    \endtable

Note
    - The \c RH expressions are special cases of those in \c YGCJ when \c C1=0
    and \c C2=1. Both \c C1 and \c C2 can be determined by nonlinear fitting
    of (\c YGCJ:Eqs. 19-20) with an experimental dataset for \c k. By default,
    \c pasquillAtmBoundaryLayerInlet boundary conditions compute \c RH expressions.
    - \c z is the ground-normal height relative to the global minimum height
    of the inlet patch; therefore, the minimum of \c z is always zero
    irrespective of the absolute z-coordinate of the computational patch.
    - The derived ABL expressions automatically satisfy the simplified transport
    equation for \c k. Yet the same expressions only satisfy the simplified
    transport equation for \c epsilon when the model constants \c sigmaEpsilon
    is 1.11 with \c kappa=0.4 (\c HW:p. 358).
    - \c pasquillAtmBoundaryLayerInlet boundary conditions inherit \c inletOutlet
    traits, so that a given inlet condition can be supplied from all sides of
    the domain, e.g. a ground-normal cylinder domain having a single
    inlet/outlet boundary where the changes between inlet and outlet depend
    on the wind direction and patch normals, so that any change in inflow
    orientation can be handled with the same mesh.
    - \c d is the displacement height, and "is relevant for flows over forests
    and cities" (E:p. 28). "The displacement height gives the vertical
    displacement of the entire flow regime over areas which are densely covered
    with obstacles such as trees or buildings" (E:p. 28).

See also
    - Foam::pasquillAtmBoundaryLayerInletVelocityFvPatchVectorField
    - Foam::pasquillAtmBoundaryLayerInletKFvPatchScalarField
    - Foam::pasquillAtmBoundaryLayerInletEpsilonFvPatchScalarField
    - Foam::pasquillAtmBoundaryLayerInletOmegaFvPatchScalarField

SourceFiles
    pasquillAtmBoundaryLayer.C

\*---------------------------------------------------------------------------*/

#ifndef pasquillAtmBoundaryLayer_H
#define pasquillAtmBoundaryLayer_H

#include "fvPatchFields.H"
#include "Function1.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class pasquillAtmBoundaryLayer Declaration
\*---------------------------------------------------------------------------*/

class pasquillAtmBoundaryLayer
{
protected:

    // Protected Data

        //- Flag to initialise profiles with the theoretical ABL expressions,
        //- otherwise initialises by using "value" entry content
        bool initABL_;

        //- Enumeration defining the Pasquill class
        enum class PasquillClass
        {
            A,           
            B,   
            C,      
            D,
            E,
            F,
            CLASS_A,    
            CLASS_B,    
            CLASS_C,    
            CLASS_D,    
            CLASS_E,    
            CLASS_F
        };

        //- Names for the Pasquill class 
        static const Enum<PasquillClass> PasquillClassNames_;


private:

    // Private Data

        //- Pasquill stability
        PasquillClass stability_;

        //- von Kármán constant
        const scalar kappa_;

        //- Empirical model constant
        const scalar Cmu_;

        //- Minimum coordinate vector of this patch
        const vector ppMin_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the patch
        const polyPatch& patch_;

        //- Streamwise flow direction
        autoPtr<Function1<vector>> flowDir_;

        //- Direction of the ground-normal coordinate
        autoPtr<Function1<vector>> zDir_;

        //- Reference mean streamwise flow speed being used in Ustar estimations
        autoPtr<Function1<scalar>> Uref_;

        //- Reference height being used in Ustar estimations
        autoPtr<Function1<scalar>> Zref_;

        //- Surface roughness length
        autoPtr<PatchFunction1<scalar>> z0_;

        //- Displacement height
        autoPtr<PatchFunction1<scalar>> d_;

        //- Latitude
        const scalar latitude_;

        //- Mean surface heat flux
        const scalar qdot_s_;

        const scalar Cp_;

        const scalar rho_;

        const scalar Tref_;


public:

    // Constructors

        //- Construct null from time
        pasquillAtmBoundaryLayer(const Time& time, const polyPatch& pp);

        //- Construct from the time database and dictionary
        pasquillAtmBoundaryLayer
        (
            const Time& time,
            const polyPatch& pp,
            const dictionary& dict
        );

        //- Construct by mapping given pasquillAtmBoundaryLayer onto a new patch
        pasquillAtmBoundaryLayer
        (
            const pasquillAtmBoundaryLayer& abl,
            const fvPatch& patch,
            const fvPatchFieldMapper& mapper
        );

        //- Construct as copy
        pasquillAtmBoundaryLayer(const pasquillAtmBoundaryLayer&);


    // Member Functions

        // Access

            //- Return flow direction
            tmp<vectorField> flowDir(const fvPatch& fvp) const;

            //- Return the ground-normal direction
            vector zDir() const;

            //- Calculate Monin-Obukhov length from Golder profiles
            tmp<scalarField> MOL(const scalarField& z0) const;

            //- Return friction velocity
            tmp<scalarField> Ustar
            (
                const scalarField& z0,
                const scalarField& d,
                const scalarField& L
            ) const;

            tmp<scalarField> H
            (
                const scalarField& uStar,
                const scalarField& L
            ) const;

            tmp<scalarField> PhiU
            (
                const scalarField& z,
                const scalarField& L
            ) const;


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            void rmap(const pasquillAtmBoundaryLayer&, const labelList&);


        // Evaluation

            //- Return the velocity distribution for the ATM
            tmp<vectorField> U(const fvPatch& fvp) const;

            //- Return the turbulent kinetic energy distribution for the ATM
            tmp<scalarField> k(const fvPatch& fvp) const;

            //- Return the turbulent dissipation rate distribution for the ATM
            tmp<scalarField> epsilon(const fvPatch& fvp) const;

            //- Return the specific dissipation rate distribution for the ATM
            tmp<scalarField> omega(const fvPatch& fvp) const;


        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

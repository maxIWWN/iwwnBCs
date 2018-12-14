/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField

Group
    grpOutletBoundaryConditions

Description
    This boundary condition provides a hydrostatic pressure distribution + dynamic pressure (for inflow) 
    with free surface (two phase flows) at outlet for p_rgh when using interFoam. 
    with zFS = free surface z-coordinate

    calculation of p_rgh in multiphase flows:
        \f[
            p_rgh = p0 + \rho*9.81*zFS - 0.5*\rho*(1-pos(phi))*|U|^2
        \f]
        where
        \vartable
            zFS             | z-coordinate of Free Surface [m]
            p0              | total pressure [Pa] - for interFoam p0 = 0
            \rho            | density [kg/m3]
            U               | velocity [m/s]
            (1-pos(phi))    | if phi (flux) positive -> flow out of domain  -> (1-pos(phi)) = 0
                            | if phi (flux) negative -> flow into of domain -> (1-pos(phi)) = 1
        \endvartable

    \heading Patch usage

    \table
        Property     | Description                     | Required    | Default value
        U            | velocity field name             | no          | U
        phi          | flux field name                 | no          | phi
        rho          | density field name              | no          | none
        zFS          | z-coord. of free surface        | yes         |
        p0           | total pressure                  | yes
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            totalPressureFreeSurfaceP_rghInlet;
        zFS             uniform 0.2;    // Free Surface at z = 20 cm
        p0              uniform 0;      // same as p_rghRefValue (fvSolution) usually 0 - alternativ: 101325 Pa
        U               U;
        phi             phi;
        rho             rho;
    }
    \endverbatim

Author
   Max Heß 11/2018
   Institut für Waserbau und Wasserwirtschaft Nürnberg (IWWN)
   Institute for Hydraulic Engineering and Water Resources Management 

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    totalPressureFreeSurfaceP_rghInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef totalPressureFreeSurfaceP_rghInletFvPatchScalarField_H
#define totalPressureFreeSurfaceP_rghInletFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class totalPressureFreeSurfaceP_rghInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class totalPressureFreeSurfaceP_rghInletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;

        //- free surface
        scalarField zFS_;

        //- density of phase1
        scalarField rhoPhase1_;

        //- density of phase2
        scalarField rhoPhase2_;

        //- total pressure
        scalarField p0_;


public:

    //- Runtime type information
    TypeName("totalPressureFreeSurfaceP_rghInlet");


    // Constructors

        //- Construct from patch and internal field
        totalPressureFreeSurfaceP_rghInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        totalPressureFreeSurfaceP_rghInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given totalPressureFreeSurfaceP_rghInletFvPatchScalarField
        //  onto a new patch
        totalPressureFreeSurfaceP_rghInletFvPatchScalarField
        (
            const totalPressureFreeSurfaceP_rghInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        totalPressureFreeSurfaceP_rghInletFvPatchScalarField
        (
            const totalPressureFreeSurfaceP_rghInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new totalPressureFreeSurfaceP_rghInletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        totalPressureFreeSurfaceP_rghInletFvPatchScalarField
        (
            const totalPressureFreeSurfaceP_rghInletFvPatchScalarField&,
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
                new totalPressureFreeSurfaceP_rghInletFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& UName()
            {
                return UName_;
            }

            //- Return the name of the flux field
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of the flux field
            //  to allow adjustment
            word& phiName()
            {
                return phiName_;
            }

            //- Return the name of the density field
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of the density field
            //  to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return the free surface
            const scalarField& zFS() const
            {
                return zFS_;
            }

            //- Return reference to the free surface to allow adjustment
            scalarField& zFS()
            {
                return zFS_;
            }

            //- Return density of phase1
            const scalarField& rhoPhase1() const
            {
                return rhoPhase1_;
            }

            //- Return reference to density of phase1 to allow adjustment
            scalarField& rhoPhase1()
            {
                return rhoPhase1_;
            }

            //- Return density of phase2
            const scalarField& rhoPhase2() const
            {
                return rhoPhase2_;
            }

            //- Return reference to density of phase2 to allow adjustment
            scalarField& rhoPhase2()
            {
                return rhoPhase2_;
            }

            //- Return the total pressure
            const scalarField& p0() const
            {
                return p0_;
            }

            //- Return reference to the total pressure to allow adjustment
            scalarField& p0()
            {
                return p0_;
            }


        // Mapping functions

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


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  using the given patch free surface, total pressure and velocity fields
            virtual void updateCoeffs
            (
                const scalarField& zFSp,
                const scalarField& rhoPhase1p,
                const scalarField& rhoPhase2p,
                const scalarField& p0p,
                const vectorField& Up
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

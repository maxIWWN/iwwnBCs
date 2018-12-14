/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Author
   Max Heß 11/2018
   Institut für Waserbau und Wasserwirtschaft Nürnberg (IWWN)
   Institute for Hydraulic Engineering and Water Resources Management 

\*---------------------------------------------------------------------------*/

#include "totalPressureFreeSurfaceP_rghInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField::totalPressureFreeSurfaceP_rghInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("none"),
    zFS_(p.size(), 0.0),
    rhoPhase1_(p.size(), 0.0),
    rhoPhase2_(p.size(), 0.0),
    p0_(p.size(), 0.0)
{}


Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField::totalPressureFreeSurfaceP_rghInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "none")),
    zFS_("zFS", dict, p.size()),
    rhoPhase1_("rhoPhase1", dict, p.size()),
    rhoPhase2_("rhoPhase2", dict, p.size()),
    p0_("p0", dict, p.size())
{
    fvPatchField<scalar>::operator=(zFS_);
    fvPatchField<scalar>::operator=(rhoPhase1_);
    fvPatchField<scalar>::operator=(rhoPhase2_);
    fvPatchField<scalar>::operator=(p0_);
}


Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField::totalPressureFreeSurfaceP_rghInletFvPatchScalarField
(
    const totalPressureFreeSurfaceP_rghInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    zFS_(ptf.zFS_, mapper),
    rhoPhase1_(ptf.rhoPhase1_, mapper),
    rhoPhase2_(ptf.rhoPhase2_, mapper),
    p0_(ptf.p0_, mapper)
{}


Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField::totalPressureFreeSurfaceP_rghInletFvPatchScalarField
(
    const totalPressureFreeSurfaceP_rghInletFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    zFS_(tppsf.zFS_),
    rhoPhase1_(tppsf.rhoPhase1_),
    rhoPhase2_(tppsf.rhoPhase2_),
    p0_(tppsf.p0_)
{}


Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField::totalPressureFreeSurfaceP_rghInletFvPatchScalarField
(
    const totalPressureFreeSurfaceP_rghInletFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    zFS_(tppsf.zFS_),
    rhoPhase1_(tppsf.rhoPhase1_),
    rhoPhase2_(tppsf.rhoPhase2_),
    p0_(tppsf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    zFS_.autoMap(m);
    rhoPhase1_.autoMap(m);
    rhoPhase2_.autoMap(m);
    p0_.autoMap(m);
}


void Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const totalPressureFreeSurfaceP_rghInletFvPatchScalarField& tiptf =
        refCast<const totalPressureFreeSurfaceP_rghInletFvPatchScalarField>(ptf);

    zFS_.rmap(tiptf.zFS_, addr);
    rhoPhase1_.rmap(tiptf.rhoPhase1_, addr);
    rhoPhase2_.rmap(tiptf.rhoPhase2_, addr);
    p0_.rmap(tiptf.p0_, addr);
}


void Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField::updateCoeffs
(
    const scalarField& zFSp,
    const scalarField& rhoPhase1p,
    const scalarField& rhoPhase2p,
    const scalarField& p0p,
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<scalar>& rho =
        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    if (rho > (rhoPhase1p+rhoPhase2p)/2)
    {
        operator==(p0p + rho*9.81*zFSp - 0.5*rho*(1.0 - pos(phip))*magSqr(Up));

        // info: pos(phi)       -> if      phi >= 0 -> pos(phi) = 1
        //                                 phi <  0 -> pos(phi) = 0
        // ->   if phi (flux) negativ - flow into domain - dynamic pressure is taken into account with 0.5*rho*|U|²
        //      if phi (flux) positiv - flow out of domain - no dynamic pressure
    }
    else
    {

        zeroGradientFvPatchScalarField::updateCoeffs();
    
    }
}

void Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField::updateCoeffs()
{
    updateCoeffs
    (
        zFS(),
        rhoPhase1(),
        rhoPhase2(),
        p0(),
        patch().lookupPatchField<volVectorField, vector>(UName())
    );
}


void Foam::totalPressureFreeSurfaceP_rghInletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    zFS_.writeEntry("zFS", os);
    rhoPhase1_.writeEntry("rhoPhase1", os);
    rhoPhase2_.writeEntry("rhoPhase2", os);
    p0_.writeEntry("p0", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalPressureFreeSurfaceP_rghInletFvPatchScalarField
    );
}

// ************************************************************************* //

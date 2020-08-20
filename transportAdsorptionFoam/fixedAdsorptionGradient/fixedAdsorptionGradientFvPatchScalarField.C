/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "fixedAdsorptionGradientFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "unitConversion.H"



 // * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //
 Foam::volScalarField&
 Foam::fixedAdsorptionGradientFvPatchScalarField::YadsField
 (
     const word& fieldName,
     const fvMesh& mesh
 )
 {
     volScalarField* ptr = mesh.getObjectPtr<volScalarField>(fieldName);
  
     if (!ptr)
     {
         ptr = new volScalarField
         (
             IOobject
             (
                 fieldName,
                 mesh.time().timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh,
             dimensionedScalar(dimLength, Zero)
         );
  
         ptr->store();
     }
  
     return *ptr;
 }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedAdsorptionGradientFvPatchScalarField::
fixedAdsorptionGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchField<scalar>(p, iF),
    Kads_(Zero),
    Kdes_(Zero),
    Gamma_(Zero),
    Yads_(patch().size(), Zero),
    curTimeIndex_(-1)
{}


Foam::fixedAdsorptionGradientFvPatchScalarField::
fixedAdsorptionGradientFvPatchScalarField
(
    const fixedAdsorptionGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<scalar>(ptf, p, iF, mapper),
    Kads_(ptf.Kads_),
    Kdes_(ptf.Kdes_),
    Gamma_(ptf.Gamma_),
    Yads_(ptf.Yads_, mapper),
    curTimeIndex_(-1)
{}


Foam::fixedAdsorptionGradientFvPatchScalarField::
fixedAdsorptionGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<scalar>(p, iF, dict),
    Kads_(readScalar(dict.lookup("Kads"))),
    Kdes_(readScalar(dict.lookup("Kdes"))),
    Gamma_(readScalar(dict.lookup("Gamma"))),
    Yads_("Yads", dict, p.size()),
    curTimeIndex_(-1)
{
    // const scalarField& magSf = patch().magSf();
    // Info<< "Initial mass on surface " 
    //     << gSum(Yads_*magSf)
    //     << endl;
}


Foam::fixedAdsorptionGradientFvPatchScalarField::
fixedAdsorptionGradientFvPatchScalarField
(
    const fixedAdsorptionGradientFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchField<scalar>(ptf),
    Kads_(ptf.Kads_),
    Kdes_(ptf.Kdes_),
    Gamma_(ptf.Gamma_),
    Yads_(ptf.Yads_),
    curTimeIndex_(-1)
{}


Foam::fixedAdsorptionGradientFvPatchScalarField::
fixedAdsorptionGradientFvPatchScalarField
(
    const fixedAdsorptionGradientFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchField<scalar>(ptf, iF),
    Kads_(ptf.Kads_),
    Kdes_(ptf.Kdes_),
    Gamma_(ptf.Gamma_),
    Yads_(ptf.Yads_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedAdsorptionGradientFvPatchScalarField::updateCoeffs()
{
    //if (updated())
    if (updated() || (curTimeIndex_ == this->db().time().timeIndex()))
    {
        return;
    }


    const polyMesh& mesh = patch().boundaryMesh().mesh();

    const dictionary& transportProperties =  this->db().objectRegistry::lookupObject<IOdictionary>("transportProperties");
    const dimensionedScalar DT_("DT", dimViscosity, transportProperties);

    // need to update the gradient before the evaluation
    // Info << "******************************************************* " << endl;
    // Info << "Adsorption model applied for " << this->patch().name() << endl;
	// const label patchid =  this->patch().index();
    // const volScalarField& T_ = this->db().objectRegistry::lookupObject<volScalarField>("T");


    scalarField& Tp(*this); //this give the present T i think on surface???
    const scalar deltaT(db().time().deltaTValue());

    scalarField Nads(patch().size(), Zero);
    scalarField Ndes(patch().size(), Zero);
    scalarField Theta(patch().size(), Zero);

    Theta = Yads_/Gamma_;
    Theta = min(Theta, 1.0); // limit the theta value


    Nads = Kads_*Tp*(1-Theta); // adsorption flux
    Ndes = Kdes_*Yads_; // desorption flux

    //  calculate surface concentration
    Yads_ += (Nads-Ndes)*deltaT;
    //Yads_ += 1; Info<< "Yads = " << Yads_ << endl;

    gradient() = (Nads-Ndes)/-DT_.value(); // return gradient for the bulk transport

    // Output the surface concentration data as volumeScalar 
    // for visualization (only)
    scalarField& YadsVis_ =
        YadsField
        (
            "YadsVis",
            refCast<const fvMesh>(mesh)
        ).boundaryFieldRef()[patch().index()];
    YadsVis_ = Yads_; // pass the value for visualization

    // return
    fixedGradientFvPatchField<scalar>::updateCoeffs();

    // if (debug)
    // {
        // forAll(Yads_, faceid)
        // {
        // Info<< "faceid = " << faceid 
        //     << " Tp = "<< Tp[faceid] 
        //     << " T_boundary = " 
        //     << T_.boundaryField()[patchid][faceid] 
        //     << " Yads = " 
        //     << Yads_[faceid] 
        //     << endl;
        // }

        // some debugs
        const scalarField& magSf = patch().magSf();
        const scalarField q(-DT_.value()*snGrad());
        const scalar massflowrate = gSum(q*magSf); // mass flowrate outof bulk
        const scalar massflowrate_surface = gSum((Nads-Ndes)*magSf); // mass into surface
        const scalar mass_on_surface = gSum(Yads_*magSf); // mass into surface

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " wall temperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << " massflowrate outof bulk [M/T]:" << massflowrate
            << " massflowrate to surface [M/T]:" << massflowrate_surface
            << " mass on surface [M/T]:" << mass_on_surface
            << endl;
    // }
    
    curTimeIndex_ = this->db().time().timeIndex();
}


void Foam::fixedAdsorptionGradientFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchField<scalar>::write(os);
    os.writeEntry("Kads", Kads_);
    os.writeEntry("Kdes", Kdes_);
    os.writeEntry("Gamma", Kdes_);
    Yads_.writeEntry("Yads", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       fixedAdsorptionGradientFvPatchScalarField
   );
}


// ************************************************************************* //

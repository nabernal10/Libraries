/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "specHumSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"
#include "fvm.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = e97fe6ced1b0ca40c9289d9da8b1add0a8482aa2
//
// unique function name that can be checked if the correct library version
// has been loaded

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(specHumSource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    specHumSource,
    dictionary
);

}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::specHumSource::specHumSource
(
    const Foam::word& name,
    const Foam::word& modelType,
    const Foam::dictionary& dict,
    const Foam::fvMesh& mesh

)
:
    fv::cellSetOption(name, modelType, dict, mesh),
	Cp0_("Cp0", dimSpecificHeatCapacity, 1003.5),  // Specific heat capacity of air [J/(kg K)]
	L_("L", dimLength, 0.1),  // Characteristic leaf size [m]
	C_("C", sqrt(dimTime)/dimLength, 130), // Proportionality factor [s^0.5/m]
	L_v_("L_v", dimensionSet(0,2,-2,0,0,0,0), 2.5e+6)  // Latent heat of vaporization [J/kg]

{
    read(dict);           // <<-- Add
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::specHumSource::~specHumSource()
{
    Info << "Destructor: specHumSource" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::fv::specHumSource::read(const dictionary& dict)
{
	Info << " 🟣 access to specHumSource::read()" << endl;
    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }
	
    coeffs_.readEntry("fields", fieldNames_);
    applied_.resize(fieldNames_.size(), false);
	
	Info << " 🟣 return = true from read" << endl;

    return true;
}


void Foam::fv::specHumSource::addSup
(
    const volScalarField& rho,
	fvMatrix<scalar>& eqn,
    const label fieldi
)
    {
        
        Info << " 🟣 specHumSource::addSup() aplicando fuente en el campo escalar W" << endl;
		
		const dimensionedScalar timeUnit
    (
        dimTime,
        1
    );


		//const volScalarField& specHum = mesh_.lookupObject<volScalarField>("specHum");
		const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
		const volScalarField& TLeaf = mesh_.lookupObject<volScalarField>("TLeaf");
		const volScalarField& rho_ = mesh_.lookupObject<volScalarField>("rho");
		const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
		const volScalarField& LAD = mesh_.lookupObject<volScalarField>("LAD");
		const volScalarField& GLeaf = mesh_.lookupObject<volScalarField>("GLeaf");
		
		// Aerodinamic resistance - r_a [s/m]
        const dimensionedScalar Umin(dimVelocity, 0.5);
        const dimensionedScalar Umax(dimVelocity, 1000);
        volScalarField Umag(mag(U));
        Umag.clip(Umin, Umax);
        volScalarField r_a = C_ * sqrt(L_ / Umag);		
		
		// Convective heat transfer coefficient - h_ch [W/(m^2 K)]
		volScalarField h_ch = (2.0 * rho_ * Cp0_ ) / r_a;
		
		// Sensible heat flux from the leaf - qPlantSen [W/m^2]
		volScalarField qPlantSen = h_ch * (TLeaf - T);
		
		// Latent heat flux from the leaf - qPlantLat [W/m^2]
		volScalarField qPlantLat = GLeaf - qPlantSen;
		
		// Vapour mass flux from the leaf (TRANSPIRATION MODEL) - g_vLeaf [kg/(s m^2)]       
		volScalarField g_vLeaf = qPlantLat / L_v_;
		
		// Humidity ratio source term - Sw [kg_waterVapourMass / (kg_dryAirMass s)]
		volScalarField Sw = LAD * g_vLeaf / rho_;
		
		eqn += Sw;
		
    }
	
void Foam::fv::specHumSource::addSup //(Para que fvOptions detecte specHumSource)
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    return this->addSup(volScalarField::null(), eqn, fieldi);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
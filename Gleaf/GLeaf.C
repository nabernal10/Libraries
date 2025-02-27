/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "GLeaf.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(GLeaf, 0);
    addToRunTimeSelectionTable(functionObject, GLeaf, dictionary);
}
}

// * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::functionObjects::GLeaf::GLeaf
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
	
    // Specific heat capacity of air - Cp0 [J/(kg·K)]
    Cp0_("Cp0", dimSpecificHeatCapacity, 1003.5),  

    // Leaf absorptance - α_leaf [-]
    alphaLeaf_(0.72), // The leaf absorbs 75% of the incident radiation

    // Leaf reflectance (albedo) - ρ_leaf [-]
    rhoLeaf_(0.23), // The leaf reflects approximately 25% of the incident radiation

    // Transmittance through canopy - τ [-]
	taoLeaf_(0.05),
	
	// Leaf emissivity - ε_leaf [-]
    epsilonLeaf_(0.95) // Typical emissivity value for Pachira macrocarpa
	
	// Solar elevation angle - β [degrees]
	// betaDeg_(90.0) // Default solar elevation angle in degrees
	
{
    read(dict);
}

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::GLeaf::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        return true;
    }

    return false;
}

bool Foam::functionObjects::GLeaf::execute()
{
    const auto& T = mesh_.lookupObject<volScalarField>("T");
	const volScalarField& TLeaf = mesh_.lookupObject<volScalarField>("TLeaf");
	const volScalarField& rho_ = mesh_.lookupObject<volScalarField>("rho");
	//const auto& LAD = mesh_.lookupObject<volScalarField>("LAD");
	
	const dimensionedScalar StefanBoltzmann
    (
        dimPower/(sqr(dimLength)*pow4(dimTemperature)),
        5.670374e-8
    );

    const dimensionedScalar G_
    (
        "G_constant",
        dimPower/sqr(dimLength),
        100.0
    );
    
    Info << "Using constant G: " << G_.value() << " W/m²" << nl;
    Info << "Calculating the incident radiation on the leaf (GLeaf)" << nl;

    tmp<volScalarField> tGLeaf
    (
        volScalarField::New
        (
            "GLeaf",
            T.mesh(),
            dimPower/sqr(dimLength)  // [W/m²]
        )
    );
    volScalarField& GLeaf = tGLeaf.ref();


	// Extinction coefficient - kappa [1/m]
	// volScalarField betaRad = betaDeg_ * (constant::mathematical::pi / 180.0); // Conversion from degrees to radians
	// volScalarField kappa = 1.0 / (2.0 * sin(betaRad)); 
	
	// Transmittance through canopy - τ [-]
	// volScalarField tao = exp (-kappa * LAD);

	// Absorbed shortwave radiation - Rs,absorbed [W/m²]
	dimensionedScalar Rs_absorbed = G_ * alphaLeaf_;

	// Absorbed longwave radiation - RL,absorbed [W/m²]
	volScalarField RL_absorbed = epsilonLeaf_ * StefanBoltzmann * pow(T, 4);

	// Reflected shortwave radiation - Rs,reflected [W/m²]
	dimensionedScalar Rs_reflected = G_ * rhoLeaf_;
	
	// Transmitted shortwave radiation through canopy - Rs,transmitted [W/m²]
	dimensionedScalar Rs_transmitted  = G_ * taoLeaf_;

	// Outgoing longwave radiation - RL_outgoing [W/m²]
	volScalarField RL_outgoing = epsilonLeaf_ * StefanBoltzmann * pow(TLeaf, 4);

	// Net radiation - R_n [W/m²]
	volScalarField R_n = Rs_absorbed + RL_absorbed - Rs_reflected - Rs_transmitted  - RL_outgoing;

	// Isothermal net radiation - R_ni [W/m²]
	volScalarField R_ni = R_n + epsilonLeaf_ * StefanBoltzmann * (pow(TLeaf, 4) - pow(T, 4));

	// Conductance to radiative heat transfer - g_R [m/s]
	volScalarField g_R = 4 * epsilonLeaf_ * StefanBoltzmann * pow(T, 3) / (rho_ * Cp0_);
	
	// Radiative resistance - r_R [s/m]
	volScalarField r_R = 1 / g_R;
	
	// Sensible heat loss - H [W/m²]
	volScalarField H = (rho_ * Cp0_ * (TLeaf - T) / r_R);

	// Leaf energy balance equation - GLeaf [W/m²]
	GLeaf = R_ni - H;

    word fieldName("GLeaf");
    return store(fieldName, tGLeaf);    
}

bool Foam::functionObjects::GLeaf::write()
{
    return writeObject("GLeaf");
}

// ************************************************************************* //
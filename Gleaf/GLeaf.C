/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
  Copyright (C) 2019-2023 OpenFOAM Foundation
  Copyright (C) 2021-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
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
#include "surfaceMesh.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "fvMesh.H"
#include "fvc.H"

// * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(GLeaf, 0);
    addToRunTimeSelectionTable(functionObject, GLeaf, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField Foam::functionObjects::GLeaf::plantFilter
(
    const volScalarField& LAD
) const
{
    volScalarField plantFilter = LAD * dimensionedScalar(dimLength, 1.0);

    forAll(plantFilter, i)
    {
        if(plantFilter[i] != 0)
        {
            plantFilter[i] = 1;
        } 
		else 
		{
            plantFilter[i] = 0;
        }
    }

    return plantFilter;
}

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * //

/**
 * Constructor: Initializes the GLeaf function object.
 *
 * @param name       Function object name.
 * @param runTime    Simulation time.
 * @param dict       Dictionary with user-defined parameters.
 */
Foam::functionObjects::GLeaf::GLeaf
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
	Cp0_("Cp0", dimSpecificHeatCapacity, 0.0),   /**< Specific heat capacity - Cp₀ [J/(kg·K)] */
    alphaLeaf_("alphaLeaf", dimless, 0.0),       /**< Absorptance coefficient - α_leaf [-]: Fraction of incident radiation absorbed by the leaf */
    rhoLeaf_("rhoLeaf", dimless, 0.0),           /**< Reflectance coefficient (albedo) - ρ_leaf [-]: Fraction of incident radiation reflected by the leaf */
    taoLeaf_("taoLeaf", dimless, 0.0),           /**< References transmittance coefficient - τ_leaf [-]: Fraction of incident radiation transmitted through the leaf */    
    kappaLeaf_SW_("kappaLeaf_SW", dimless, 0.0), /**< Shortwave radiation extinction coefficient - k [-]: Governs the attenuation of shortwave radiation within the canopy, considering solar directionality and leaf optical properties */
    epsilonLeaf_("epsilonLeaf", dimless, 0.0),   /**< Leaf emissivity coefficient - ε_leaf [-]: Fraction of thermal radiation emitted by the leaf relative to a black body */
    epsilonSky_("epsilonSky", dimless, 0.0),     /**< Sky emissivity coefficient - ε_sky [-]: Fraction of thermal radiation emitted by the sky relative to a black body */
	C_lw_("C_lw", dimless, 0.0),                 /**< Longwave radiation empirical constant - C_lw [-]: Empirical constant for quantifying the net absorption of longwave radiation */
	H_("H", dimLength, 0.0),                     /**< Height of the tree - h [m]: Defines the total height of the vegetation canopy, determining the vertical extent over which radiation attenuation occurs. */
	h_("h", dimLength, 0.0)                      /**< High of the base of the canopy */
{
    read(dict);
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

/**
 * Reads parameters from the dictionary.
 *
 * @param dict   Dictionary containing the function object settings.
 * @return       True if read successfully, false otherwise.
 */
bool Foam::functionObjects::GLeaf::read(const dictionary& dict)
{
    IOdictionary leafDict
    (
        IOobject
        (
            "leafProperties",
            this->mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (fvMeshFunctionObject::read(dict))
    {
        if (leafDict.found("leafProperties"))
        {
            // Extract the subdictionary leafProperties
            const dictionary& gProps = leafDict.subDict("leafProperties");

            Cp0_.readIfPresent(gProps);
            alphaLeaf_.readIfPresent(gProps);
            rhoLeaf_.readIfPresent(gProps);
            taoLeaf_.readIfPresent(gProps);
            kappaLeaf_SW_.readIfPresent(gProps);
            epsilonLeaf_.readIfPresent(gProps);
            epsilonSky_.readIfPresent(gProps);
            C_lw_.readIfPresent(gProps);
            H_.readIfPresent(gProps);
            h_.readIfPresent(gProps);

            Info << "\n ✅ Values after reading leafProperties (GLeaf Library):\n"
                 << "Cp0 = " << Cp0_.value() << "\n"
                 << "alphaLeaf = " << alphaLeaf_.value() << "\n"
                 << "rhoLeaf = " << rhoLeaf_.value() << "\n"
                 << "taoLeaf = " << taoLeaf_.value() << "\n"
                 << "kappaLeaf_SW = " << kappaLeaf_SW_.value() << "\n"
                 << "epsilonLeaf = " << epsilonLeaf_.value() << "\n"
                 << "epsilonSky = " << epsilonSky_.value() << "\n"
                 << "C_lw = " << C_lw_.value() << "\n"
                 << "H = " << H_.value() << "\n"
                 << "h = " << h_.value() << "\n";
        }
        else
        {
            Info << "\n❌ Subdictionary 'leafProperties' not found in leafProperties (GLeaf Library).\n";
        }

        return true;
    }

    return false;
}

/**
 * Executes the GLeaf calculation.
 *
 * Computes incident radiation balance for leaf surfaces.
 * 
 * @return True if execution is successful.
 */

bool Foam::functionObjects::GLeaf::execute()
{
    const auto& T = mesh_.lookupObject<volScalarField>("T");
	//const auto& G = lookupObject<volScalarField>("G");
	const auto& LAD = mesh_.lookupObject<volScalarField>("LAD");
    const volScalarField& TLeaf = mesh_.lookupObject<volScalarField>("TLeaf");
    const volScalarField& rho_ = mesh_.lookupObject<volScalarField>("rho");

    /** Stefan-Boltzmann constant - σ [W/(m²·K⁴)] */
    const dimensionedScalar StefanBoltzmann
    (
        dimPower/(sqr(dimLength)*pow4(dimTemperature)),
        5.670374e-8
    );

    /** Incident radiation constant - G [W/m²] */
    const dimensionedScalar G_
    (
        "G_constant",
        dimPower/sqr(dimLength),
        99 // 235
    );

    Info << "Using constant G: " << G_.value() << " W/m²" << nl;
    Info << "Calculating the incident radiation on the leaf (GLeaf)" << nl;
	
    /** GLeaf output field [W/m²] */
    tmp<volScalarField> tGLeaf
    (
        volScalarField::New
        (
            "GLeaf",
            T.mesh(),
            dimPower/sqr(dimLength)
        )
    );
    volScalarField& GLeaf = tGLeaf.ref();
	
	/** qr_sw output field [W/m²] */
    tmp<volScalarField> tqr_sw
    (
        volScalarField::New
        (
            "qr_sw",
            T.mesh(),
            dimPower/sqr(dimLength)
        )
    );
    volScalarField& qr_sw = tqr_sw.ref();
	
	/** qr_lw output field [W/m²] */
    tmp<volScalarField> tqr_lw
    (
        volScalarField::New
        (
            "qr_lw",
            T.mesh(),
            dimPower/sqr(dimLength)
        )
    );
    volScalarField& qr_lw = tqr_lw.ref();
	
	/** Ga output field [W/m³] */
    tmp<volScalarField> tGa
    (
        volScalarField::New
        (
            "Ga",
            T.mesh(),
            dimensionSet(1, -1, -3, 0, 0, 0, 0)
        )
    );
    volScalarField& Ga = tGa.ref();
	
	
	Info<< "\nRADIATION MODEL 1\n" << endl;
	/** RADIATION MODEL 1 */
	
	/** Transmittance coefficient through the canopy - τ [-] */
	volScalarField y = mesh_.C().component(1);
	volScalarField tao = exp(- kappaLeaf_SW_ * LAD * (H_ - y));
	Info << "Current tao: min = " << min(tao).value() << ", max = " << max(tao).value() << endl;

	/** Transmitted shortwave radiation (Beer-Lambert law) - Rsₜ [W/m²] */
	volScalarField Rs_transmitted_M1 = G_ * tao;
	Info << "Current Rs_transmitted_M1 [W/m²]: min = " << min(Rs_transmitted_M1).value() << ", max = " << max(Rs_transmitted_M1).value() << endl;
	
	/** Reflected shortwave radiation - Rsᵣ [W/m²] */
    volScalarField Rs_reflected_M1 = Rs_transmitted_M1 * rhoLeaf_;
	Info << "Current Rs_reflected_M1 [W/m²]: min = " << min(Rs_reflected_M1).value() << ", max = " << max(Rs_reflected_M1).value() << endl;
	
    /** Absorbed sky longwave radiation - RLₐ_sky [W/m²] */
    volScalarField T_sky = T - dimensionedScalar("offset", dimensionSet(0,0,0,1,0), 15.0);
	volScalarField RL_absorbed_Sky_M1 = C_lw_ * epsilonSky_ * StefanBoltzmann * pow(T_sky, 4);
	const dimensionedScalar RL_validation(dimPower/sqr(dimLength), 525);
	dimensionedScalar RL_absorbed_Sky_validation = C_lw_ * RL_validation;
	Info << "Current RL_absorbed_Sky_M1 [W/m²]: min = " << min(RL_absorbed_Sky_M1).value() << ", max = " << max(RL_absorbed_Sky_M1).value() << endl;
	Info << "Current RL_absorbed_Sky_validation [W/m²] = " << RL_absorbed_Sky_validation.value() << endl;
	
    /** Density of net radiation - Ga [W/m³] */
	volScalarField Ga_M1 = kappaLeaf_SW_ * LAD * (Rs_transmitted_M1 - Rs_reflected_M1) + RL_absorbed_Sky_M1 / (H_ - h_);
	volScalarField Ga_validation = kappaLeaf_SW_ * LAD * (Rs_transmitted_M1 - Rs_reflected_M1) + RL_absorbed_Sky_validation / (H_ - h_);
	Info << "Current Ga_M1 [W/m³]: min = " << min(Ga_M1).value() << ", max = " << max(Ga_M1).value() << endl;
	Info << "Current Ga_validation [W/m³]: min = " << min(Ga_validation).value() << ", max = " << max(Ga_validation).value() << endl;
	
	/** Net radiative flux at the leaf surface [W/m²] */
	const dimensionedScalar LAD_0("LAD_0", dimensionSet(0,-1,0,0,0,0,0), 8.77);
	volScalarField q_rad = Rs_transmitted_M1 + RL_absorbed_Sky_validation;
	volScalarField q_rad_M1_verify = Ga_M1 / LAD_0;
	volScalarField q_rad_validation_verify = Ga_validation / LAD_0;
	Info << "Current q_rad [W/m²]: min = " << min(q_rad).value() << ", max = " << max(q_rad).value() << endl;
	Info << "Current q_rad_M1_verify [W/m²]: min = " << min(q_rad_M1_verify).value() << ", max = " << max(q_rad_M1_verify).value() << endl;
	Info << "Current q_rad_validation_verify [W/m²]: min = " << min(q_rad_validation_verify).value() << ", max = " << max(q_rad_validation_verify).value() << endl;
		
	volScalarField q_rad_M1
	(
		IOobject("q_rad_M1", this->mesh_.time().timeName(), this->mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
		this->mesh_,
		dimensionedScalar("zero", dimPower/sqr(dimLength), 0.0)
	);

	volScalarField q_rad_validation
	(
		IOobject("q_rad_validation", this->mesh_.time().timeName(), this->mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
		this->mesh_,
		dimensionedScalar("zero", dimPower/sqr(dimLength), 0.0)
	);

	forAll(mesh_.C(), i)
	{
		
		if (LAD[i] > SMALL)
		{
			q_rad_M1[i] = Ga_M1[i] / LAD[i];
			q_rad_validation[i] = Ga_validation[i] / LAD[i];
		}
		else
		{
			q_rad_M1[i] = 0.0;
			q_rad_validation[i] = 0.0;
		}
	}
	
	Info << "Current q_rad_M1 [W/m²]: min = " << min(q_rad_M1).value() << ", max = " << max(q_rad_M1).value() << endl;
	Info << "Current q_rad_validation [W/m²]: min = " << min(q_rad_validation).value() << ", max = " << max(q_rad_validation).value() << endl;
	
	
	Info<< "\nRADIATION MODEL 2\n" << endl;
	/** RADIATION MODEL 2 */
		
	/** Reflected shortwave radiation - Rsᵣ [W/m²] */
    dimensionedScalar Rs_transmitted_M2 = G_ * taoLeaf_;
	Info << "Rs_transmitted_M2 [W/m²] = " << Rs_transmitted_M2.value() << endl;
	
	/** Reflected shortwave radiation - Rsᵣ [W/m²] */
    dimensionedScalar Rs_reflected_M2 = G_ * rhoLeaf_;
	Info << "Rs_reflected_M2 [W/m²] = " << Rs_reflected_M2.value() << endl;
	
    /** Absorbed sky longwave radiation - RLₐ_sky [W/m²] */
	volScalarField RL_absorbed_Sky_M2 = epsilonSky_ * StefanBoltzmann * pow(T_sky, 4);
	Info << "Current RL_absorbed_Sky_M2 [W/m²]: min = " << min(RL_absorbed_Sky_M2).value() << ", max = " << max(RL_absorbed_Sky_M2).value() << endl;
	
	/** Emitted longwave radiation - RLₒ [W/m²] */
    volScalarField RL_emitted_M2 = epsilonLeaf_ * StefanBoltzmann * pow(TLeaf, 4);
	Info << "Current RL_emitted_M2 [W/m²]: min = " << min(RL_emitted_M2).value() << ", max = " << max(RL_emitted_M2).value() << endl;
	
	/** Net radiation - Rₙ [W/m²] */
    volScalarField R_n = G_ - Rs_transmitted_M2 - Rs_reflected_M2 + RL_absorbed_Sky_M2 - RL_emitted_M2;
	Info << "Current R_n: min [W/m²] = " << min(R_n).value() << ", max = " << max(R_n).value() << endl;

    /** Isothermal net radiation - Rₙᵢ [W/m²] */
    volScalarField R_ni = R_n + (epsilonLeaf_ * StefanBoltzmann * (pow(TLeaf, 4) - pow(T, 4)));
	Info << "Current R_ni [W/m²]: min = " << min(R_ni).value() << ", max = " << max(R_ni).value() << endl;

    /** Radiative heat transfer conductance - gR [m/s] */
    volScalarField g_R = 4 * epsilonLeaf_ * StefanBoltzmann * pow(T, 3) / (rho_ * Cp0_);
	Info << "Current g_R [m/s]: min = " << min(g_R).value() << ", max = " << max(g_R).value() << endl;

    /** Radiative heat transfer resistance - rR [s/m] */
    volScalarField r_R = 1 / g_R;
	Info << "Current r_R [s/m]: min = " << min(r_R).value() << ", max = " << max(r_R).value() << endl;

    /** Sensible heat loss - H_l [W/m²] */
    volScalarField H_l = rho_ * Cp0_ * (TLeaf - T) / r_R;
	Info << "Current H [W/m²]: min = " << min(H_l).value() << ", max = " << max(H_l).value() << endl;
	
	volScalarField q_rad_M2 = R_ni - H_l;
	Info << "Current q_rad_M2 [W/m²]: min = " << min(q_rad_M2).value() << ", max = " << max(q_rad_M2).value() << endl;


	Info<< "\nRESULTS\n" << endl;
	/** ✅ Leaf energy balance equation - GLeaf [W/m²] */
    //GLeaf = q_rad_M1;
	GLeaf = q_rad_validation;
	//GLeaf = q_rad_M2;
	
	Info << "Current GLeaf [W/m²]: min = " << min(GLeaf).value() << ", max = " << max(GLeaf).value() << endl;
	
	/** ✅ Density of net radiation - Ga [W/m³] */
	Ga = Ga_validation;
	Info << "Current Ga [W/m³]: min = " << min(Ga).value() << ", max = " << max(Ga).value() << endl;
	
	/** ✅ Shortwave net radiation - qr_sw [W/m²] */
	qr_sw = Rs_transmitted_M1 - Rs_reflected_M1;
	//qr_sw = G_ - Rs_transmitted_M2 - Rs_reflected_M2;
	Info << "Current qr_sw [W/m²]: min = " << min(qr_sw).value() << ", max = " << max(qr_sw).value() << endl;
	
	/** ✅ Longwave net radiation - qr_lw [W/m²] */	
	//qr_lw = RL_absorbed_Sky_M1;
	qr_lw = RL_absorbed_Sky_validation;
	//qr_lw = RL_absorbed_Sky_M2 - RL_emitted_M2;
	Info << "Current qr_lw [W/m²]: min = " << min(qr_lw).value() << ", max = " << max(qr_lw).value() << endl;
	
    Foam::word fieldNameGLeaf("GLeaf");
	Foam::word fieldNameqr_sw("qr_sw");
	Foam::word fieldNameqr_lw("qr_lw");
	Foam::word fieldNameGa("Ga");

	return store(fieldNameGLeaf, tGLeaf) && store(fieldNameqr_sw, tqr_sw) && store(fieldNameqr_lw, tqr_lw) && store(fieldNameGa, tGa);
}

/**
 * Writes the GLeaf field to a file.
 *
 * @return True if the write operation is successful.
 */
bool Foam::functionObjects::GLeaf::write()
{
    return writeObject("GLeaf") && writeObject("qr_sw") && writeObject("qr_lw") && writeObject("Ga");
}
// ************************************************************************* //
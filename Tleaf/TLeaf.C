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

#include "TLeaf.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(TLeaf, 0);
    addToRunTimeSelectionTable(functionObject, TLeaf, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::TLeaf::converged
(
    const volScalarField& phi
) const
{
    return
        max(mag(phi.primitiveField() - phi.prevIter().primitiveField()))
      < tolerance_;
}

Foam::volScalarField Foam::functionObjects::TLeaf::plantFilter
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

/**
 * Constructor: Initializes the TLeaf function object.
 *
 * @param name       Function object name.
 * @param runTime    Simulation time.
 * @param dict       Dictionary with user-defined parameters.
 */
 
Foam::functionObjects::TLeaf::TLeaf
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
	
	fvMeshFunctionObject(name, runTime, dict),
    Cp0_("Cp0", dimSpecificHeatCapacity, 0.0),      /**< Specific heat capacity - Cp₀ [J/(kg·K)] */
	L_("L", dimLength, 0.0),                        /**< Characteristic leaf size - L [m] */
    C_("C", sqrt(dimTime)/dimLength, 0.0),          /**< Proportionality factor - C [sqrt(s)/m] */
    r_sMin_("r_sMin", dimTime/dimLength, 0.0),      /**< Minimal stomatal resistance - rₛ,Min [s/m] */
    R_a_("R_a", dimGasConstant, 0.0),               /**< Gas constant of dry air - Rₐ [J/(kg·K)] */
    R_v_("R_v", dimGasConstant, 0.0),               /**< Gas constant of water vapour - Rᵥ [J/(kg·K)] */
    L_v_("L_v", dimensionSet(0,2,-2,0,0,0,0), 0.0), /**< Latent heat of vaporization - Lᵥ [J/kg] */
    tolerance_(0.0),                                /**< Convergence tolerance for TLeaf */
    maxIter_(0.0)                                   /**< Maximum number of iterations */
	
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/**
 * Reads function object parameters from the dictionary.
 *
 * @param dict   Dictionary containing user-defined parameters.
 * @return       True if reading is successful.
 */

bool Foam::functionObjects::TLeaf::read(const dictionary& dict)
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
            L_.readIfPresent(gProps);
            C_.readIfPresent(gProps);
            r_sMin_.readIfPresent(gProps);
            R_a_.readIfPresent(gProps);
            R_v_.readIfPresent(gProps);
            L_v_.readIfPresent(gProps);
			tolerance_ = gProps.getOrDefault("tolerance", 1e-4);
            maxIter_ = gProps.getOrDefault("maxIter", 100);

            Info << "\n ✅ Values after reading leafProperties (TLeaf Library):\n"
                 << "Cp0 = " << Cp0_.value() << "\n"
				 << "L = " << L_.value() << "\n"
                 << "C = " << C_.value() << "\n"
                 << "r_sMin = " << r_sMin_.value() << "\n"
                 << "R_a = " << R_a_.value() << "\n"
                 << "R_v = " << R_v_.value() << "\n"
                 << "L_v = " << L_v_.value() << "\n"
                 << "tolerance = " << tolerance_ << "\n"
                 << "maxIter = " << maxIter_ << "\n";
        }
        else
        {
            Info << "\n❌ Subdictionary 'leafProperties' not found in leafProperties (TLeaf Library).\n";
        }

        return true;
    }

    return false;
}

/**
 * Executes the TLeaf calculation.
 *
 * Computes the temperature of the leaf based on radiation, humidity, and aerodynamics.
 * 
 * @return True if execution is successful.
 */
bool Foam::functionObjects::TLeaf::execute()
{
    // Assign and build fields
    const auto& T = mesh_.lookupObject<volScalarField>("T");
	const auto& GLeaf = mesh_.lookupObject<volScalarField>("GLeaf");	
	const auto& LAD = mesh_.lookupObject<volScalarField>("LAD");
	const auto& p = mesh_.lookupObject<volScalarField>("p");
	const auto& rho = mesh_.lookupObject<volScalarField>("rho");	
	const auto& w = mesh_.lookupObject<volScalarField>("specHum");
	const auto& U = mesh_.lookupObject<volVectorField>("U");
	const auto& qr_sw = mesh_.lookupObject<volScalarField>("qr_sw");

    Info<< "\nCalculating the Temperature of the leaf (TLeaf)" << endl;

    /** TLeaf output field [K] */
    tmp<volScalarField> tTLeaf
    (
        volScalarField::New
        (
            "TLeaf",
            T.mesh(),
            dimTemperature
        )
    );
    volScalarField& TLeaf = tTLeaf.ref();    
    TLeaf = T; // Initial guess
    label i = 0;
    TLeaf.storePrevIter();
	
	/** g_vLeaf output field [mg/(s·m²)] */
    tmp<volScalarField> tg_vLeaf_mg
    (
        volScalarField::New
        (
            "g_vLeaf",
            T.mesh(),
            dimensionSet(1, -2, -1, 0, 0, 0, 0)
        )
    );
    volScalarField& g_vLeaf_mg = tg_vLeaf_mg.ref();
	
	/** r_s output field [s/m] */
    tmp<volScalarField> tr_sW
    (
        volScalarField::New
        (
            "r_s",
            T.mesh(),
            dimensionSet(0, -1, 1, 0, 0, 0, 0)
        )
    );
    volScalarField& r_sW = tr_sW.ref();
	
	/** r_a output field [s/m] */
    tmp<volScalarField> tr_aW
    (
        volScalarField::New
        (
            "r_a",
            T.mesh(),
            dimensionSet(0, -1, 1, 0, 0, 0, 0)
        )
    );
    volScalarField& r_aW = tr_aW.ref();
	
	/** qPlantLat output field [W/m²] */
    tmp<volScalarField> tqPlantLat
    (
        volScalarField::New
        (
            "qPlantLat",
            T.mesh(),
            dimPower/sqr(dimLength)
        )
    );
    volScalarField& qPlantLat = tqPlantLat.ref();

	/** Stefan-Boltzmann constant - σ [W/(m²·K⁴)] */
    const dimensionedScalar StefanBoltzmann
    (
        dimPower/(sqr(dimLength)*pow4(dimTemperature)),
        5.670374e-8
    );

    /** Pressure unit - p₀ [Pa] */
	const dimensionedScalar pressureUnit
    (
        dimPressure,
        1
    );

    do
    {
        TLeaf.storePrevIter();
        
		Info<< "\n¡¡¡¡NEW TLeaf ITERATION!!!!\n" << endl;
        /** Net radiative flux - qPlantRad [W/m²] */
        volScalarField qPlantRad = GLeaf;
		Info << "Current qPlantRad [W/m²]: min = " << min(qPlantRad).value() << ", max = " << max(qPlantRad).value() << ", mean: " << average(qPlantRad).value() << endl;
			
		
		/** Saturation vapour pressure at air (Tetens equation) - pVSatAir [Pa] */
		const dimensionedScalar A(dimless, 17.27); // Adjustment coefficient that describes how the saturation vapor pressure changes with temperature
		const dimensionedScalar B(dimTemperature, 273.15); // Constant to convert Kelvin to Kelvin Celsius degrees
		const dimensionedScalar C(dimTemperature, 237.3); // Constant that adjusts the equation for common atmospheric temperatures
		const dimensionedScalar E(dimless, 610.78); // Saturation vapor pressure at 0°C [Pa]
		volScalarField pVSatAir = pressureUnit * E * exp((A * (T - B)) / ((T - B) + C));
		volScalarField T_C = T - B;
		Info << "Current T [K]: min = " << min(T).value() << ", max = " << max(T).value() << ", mean: " << average(T).value() << endl;
		Info << "Current T [°C]: min = " << min(T_C).value() << ", max = " << max(T_C).value() << ", mean: " << average(T_C).value() << endl;
		Info << "Current pVSatAir [Pa]: min = " << min(pVSatAir).value() << ", max = " << max(pVSatAir).value() << ", mean: " << average(pVSatAir).value() << endl;

		
		/** Partial vapour pressure at air - pVAir [Pa] */
		const dimensionedScalar pF(dimDensity, 1); //Pressure unit adjustment factor
		volScalarField pVAir = w * (p * pF) / ((R_a_ / R_v_) + w);		
		Info << "Current pVAir [Pa]: min = " << min(pVAir).value() << ", max = " << max(pVAir).value() << ", mean: " << average(pVAir).value() << endl;
		Info << "Current p [Pa]: min = " << min(p).value() << ", max = " << max(p).value() << ", mean: " << average(p).value() << endl;
		Info << "Current w [kg/kg]: min = " << min(w).value() << ", max = " << max(w).value() << ", mean: " << average(w).value() << endl;


        /** Saturation vapour pressure at leaf temperature - pVLeaf ≈ pVSatLeaf [Pa] */
        volScalarField pVSatLeaf = pressureUnit * E * exp((A * (TLeaf - B)) / ((TLeaf - B) + C)); // The vapour pressure at the leaf is the vapour pressure within the leaf stomata which is close to the saturation vapour pressure at the leaf temperature,
		volScalarField TLeaf_C = TLeaf - B;
		Info << "Current TLeaf [K]: min = " << min(TLeaf).value() << ", max = " << max(TLeaf).value() << ", mean: " << average(TLeaf).value() << endl;
		Info << "Current TLeaf [°C]: min = " << min(TLeaf_C).value() << ", max = " << max(TLeaf_C).value() << ", mean: " << average(TLeaf_C).value() << endl;
		Info << "Current pVSatLeaf [Pa]: min = " << min(pVSatLeaf).value() << ", max = " << max(pVSatLeaf).value() << ", mean: " << average(pVSatLeaf).value() << endl;


		/** Boundary layer conductance - gₐW [m/s] */
        const dimensionedScalar Umin(dimVelocity, 0.03);
        const dimensionedScalar Umax(dimVelocity, 1000);
        volScalarField Umag(mag(U));
        Umag.clip(Umin, Umax);
        volScalarField g_aW = (1 / C_) * sqrt(Umag / L_);		
		
		
		/** Boundary layer resistance - rₐw [s/m] */
        r_aW = 1 / g_aW;
		Info << "Current mag(U) [m/s]: min = " << min(mag(U)).value() << ", max = " << max(mag(U)).value() << ", mean: " << average(mag(U)).value() << endl;
		Info << "Current r_aW [s/m]: min = " << min(r_aW).value() << ", max = " << max(r_aW).value() << ", mean: " << average(r_aW).value() << endl;		
		
		
        /** Stomatal resistance [rₛW = rₛ,Min * f₁(G₀) * f₂(D)] - rₛ [s/m] */
        const dimensionedScalar a_1(dimPower/sqr(dimLength), 169); // [W/m²]
        const dimensionedScalar a_2(dimPower/sqr(dimLength), 18); // [W/m²]
        const dimensionedScalar a_3(dimless/sqr(dimPressure), 0.005); // [1/kPa²]
		const dimensionedScalar D0(dimPressure, 1.2); // [kPa]        
		volScalarField D = (pVSatAir - pVAir) / 1000; // Vapor pressure deficit [kPa]				
        r_sW = r_sMin_ * (a_1 + qr_sw) / (a_2 + qr_sw) * (1 + a_3 * sqr(D - D0));
        //r_sW = r_aW * 0 + r_sMin_;
		Info << "Current r_sW [s/m]: min = " << min(r_sW).value() << ", max = " << max(r_sW).value() << ", mean: " << average(r_sW).value() << endl;
		Info << "Current D [kPa]: min = " << min(D).value() << ", max = " << max(D).value() << ", mean: " << average(D).value() << endl;
		

		/** Convective mass transfer coefficient - hcₘ [s/m] */
		volScalarField h_cm = (1 / (r_aW + r_sW)) * (R_a_ / R_v_) * (rho / (p * pF));
		Info << "Current h_cm [s/m]: min = " << min(h_cm).value() << ", max = " << max(h_cm).value() << ", mean: " << average(h_cm).value() << endl;
		Info << "Current rho [kg/m³]: min = " << min(rho).value() << ", max = " << max(rho).value() << ", mean: " << average(rho).value() << endl;
		
		
        /** Vapour mass flux from the leaf (TRANSPIRATION MODEL) - gvLeaf [kg/(s·m²)] */
		volScalarField g_vLeaf = h_cm * (pVSatLeaf - pVAir);
		g_vLeaf_mg = g_vLeaf * 1000;
		Info << "Current g_vLeaf [kg/(s·m²)] (transpiration): min = " << min(g_vLeaf).value() << ", max = " << max(g_vLeaf).value() << ", mean: " << average(g_vLeaf).value() << endl;
		Info << "Current g_vLeaf [mg/(s·m²)] (transpiration): min = " << min(g_vLeaf_mg).value() << ", max = " << max(g_vLeaf_mg).value() << ", mean: " << average(g_vLeaf_mg).value() << endl;
		
		
		/** Latent heat flux from the leaf - qPlantLat [W/m²] */
        qPlantLat = L_v_ * g_vLeaf;
		Info << "Current qPlantLat [W/m²]: min = " << min(qPlantLat).value() << ", max = " << max(qPlantLat).value() << ", mean: " << average(qPlantLat).value() << endl;
        
		
		/** Convective heat transfer coefficient - hcₕ [W/(m²·K)] */
		volScalarField h_ch = (2 * rho * Cp0_ ) / r_aW;
		Info << "Current h_ch [W/(m²·K)]: min = " << min(h_ch).value() << ", max = " << max(h_ch).value() << ", mean: " << average(h_ch).value() << endl;
		
		
		/** ✅ Leaf energy balance equation - TLeaf [K] */
		TLeaf = T + (qPlantRad - qPlantLat) / h_ch;
		

        // Reset TLeaf to T where there is not any plant
        TLeaf = (1 - plantFilter(LAD)) * T + plantFilter(LAD) * TLeaf;
		

    } while (!converged(TLeaf) && i++ < maxIter_);

    if (i < maxIter_) {
        Info
            << "The leaf temperature converged within\n" << i
            << " iterations" << endl;
    } else if (i == maxIter_) {
        WarningInFunction
            << "The leaf temperature did not converge within" << i
            << " iterations" << nl;
    }

    Foam::word fieldNameTLeaf("TLeaf");
	Foam::word fieldNameg_vLeaf_mg("g_vLeaf");
	Foam::word fieldNamer_s("r_s");
	Foam::word fieldNamer_a("r_a");
	Foam::word fieldNameqPlantLat("qPlantLat");

    return store(fieldNameTLeaf, tTLeaf) && store(fieldNameg_vLeaf_mg, tg_vLeaf_mg) && store(fieldNamer_s, tr_sW) && store(fieldNamer_a, tr_aW) && store(fieldNameqPlantLat, tqPlantLat);
}

/**
 * Writes the TLeaf field to a file.
 *
 * @return True if the write operation is successful.
 */ 
bool Foam::functionObjects::TLeaf::write()
{
    return writeObject("TLeaf") && writeObject("g_vLeaf") && writeObject("r_s") && writeObject("r_a") && writeObject("qPlantLat");
}
// ************************************************************************* //
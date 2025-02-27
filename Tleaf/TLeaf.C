/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

Foam::functionObjects::TLeaf::TLeaf
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    Cp0_("Cp0", dimensionSet(0,2,-2,-1,0,0,0), 1003.5),  // Specific heat capacity of air [J/(kg K)]   
    //rho0_("rho0", dimDensity, 1.225),  // Density of air [kg/m^3]
	r_("r", dimDensity, 1),  // Constant to adjust pressure units [kg/m^3]
	//p0("p0", dimPressure, 101325), // Pressure of the air [Pa]
	L_("L", dimLength, 0.1),  // Characteristic leaf size [m]
    C_("C", sqrt(dimTime)/dimLength, 130), // Proportionality factor [s^0.5/m]
    r_sMin_("r_sMin", dimTime/dimLength, 150), // Minimal stomatal resistance [s/m]
    R_a_("R_a", dimGasConstant, 287.042),  // Gas constants of dry air [J/(kg K)]
	R_v_("R_v", dimGasConstant, 461.524),  // Gas constants water vapour [J/(kg K)]
	L_v_("L_v", dimensionSet(0,2,-2,0,0,0,0), 2.5e+6),  // Latent heat of vaporization [J/kg]
    //G0_("G0", dimPower/sqr(dimLength), 500),
	g("g", dimAcceleration, 9.81),
    tolerance_(1e-4),
    maxIter_(100)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::TLeaf::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        Cp0_.readIfPresent(dict);        
        //rho0_.readIfPresent(dict);
		//p0.readIfPresent(dict);
		L_.readIfPresent(dict);
        C_.readIfPresent(dict);
        r_sMin_.readIfPresent(dict);
        R_a_.readIfPresent(dict);        
		R_v_.readIfPresent(dict);        
		L_v_.readIfPresent(dict);
        //G0_.readIfPresent(dict);
        tolerance_ = dict.getOrDefault("tolerance", 1e-4);
        maxIter_ = dict.getOrDefault("maxIter", 100);

        return true;
    }

    return false;
}


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

    // Calculate the temperature of the leaf - (TLeaf) [K] quantity
    Info<< "Calculating the Temperature of the leaf (TLeaf)" << endl;

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

    // Initial guess
    TLeaf = T;

    label i = 0;

    TLeaf.storePrevIter();

    const dimensionedScalar StefanBoltzmann
    (
        dimPower/(sqr(dimLength)*pow4(dimTemperature)),
        5.670374e-8
    );

    const dimensionedScalar pressureUnit
    (
        dimPressure,
        1
    );

    do
    {
        TLeaf.storePrevIter();
        
		Info<< "\n¡¡¡¡NUEVA ITERACIÓN PARA TLeaf!!!!\n" << endl;
        // Net radiative flux to the leaf - qPlantRad [W/m^2]
        volScalarField qPlantRad = GLeaf;		
        // Info << "Dimensions of qPlantRad: " << qPlantRad.dimensions() << endl;
		Info << "Current qPlantRad: min = " << min(qPlantRad).value() << ", max = " << max(qPlantRad).value() << endl;
			
		
		// Saturation vapour pressure at air (Tetens equation) - pVSatAir [Pa]
		const dimensionedScalar A(dimless, 17.27);
		const dimensionedScalar B(dimTemperature, 273.15);
		const dimensionedScalar C(dimTemperature, 237.3);
		const dimensionedScalar E(dimless, 0.61078);
		volScalarField pVSatAir = 1000 * pressureUnit * E * exp((A * (T - B)) / ((T - B) + C));
		// Info << "Dimensions of pVSatAir: " << pVSatAir.dimensions() << endl;
		Info << "Current pVSatAir: min = " << min(pVSatAir).value() << ", max = " << max(pVSatAir).value() << endl;
		Info << "Current T: min = " << min(T).value() << ", max = " << max(T).value() << endl;

		
		// Partial vapour pressure at air (Tetens equation) - pVAir [Pa]
        volScalarField pVAir = w * pVSatAir;
		// Info << "Dimensions of pVAir: " << pVAir.dimensions() << endl;
		Info << "Current pVAir: min = " << min(pVAir).value() << ", max = " << max(pVAir).value() << endl;
		Info << "Current w: min = " << min(w).value() << ", max = " << max(w).value() << endl;		


        // Saturation vapour pressure at TLeaf (Tetens equation) - pVSatLeaf [Pa]
        volScalarField pVSatLeaf = 1000 * pressureUnit * E * exp((A * (TLeaf - B)) / ((TLeaf - B) + C));
		// Info << "Dimensions of pVSatLeaf: " << pVSatLeaf.dimensions() << endl;
		Info << "Current TLeaf: min = " << min(TLeaf).value() << ", max = " << max(TLeaf).value() << endl;
		Info << "Current pVSatLeaf: min = " << min(pVSatLeaf).value() << ", max = " << max(pVSatLeaf).value() << endl;


		// Aerodynamic conductance to heat transfer - g_H [m/s]
        const dimensionedScalar Umin(dimVelocity, 0.5);
        const dimensionedScalar Umax(dimVelocity, 1000);
        volScalarField Umag(mag(U));
        Umag.clip(Umin, Umax);
        volScalarField g_H = (1 / C_) * sqrt(Umag / L_);		
		
		// Aerodinamic resistance - r_a [s/m]
        volScalarField r_a = 1 / g_H;
		// Info << "Dimensions of r_a: " << r_a.dimensions() << endl;
		Info << "Current mag(U): min = " << min(mag(U)).value() << ", max = " << max(mag(U)).value() << endl;
		Info << "Current r_a: min = " << min(r_a).value() << ", max = " << max(r_a).value() << endl;		
		
		
        // Stomatal resistance - r_s [r_s = r_sMin * f1(G0) * f2(D)] [s/m]
        const dimensionedScalar a_1(dimPower/sqr(dimLength), 169); // [W/m^2]
        const dimensionedScalar a_2(dimPower/sqr(dimLength), 18); // [W/m^2]
        const dimensionedScalar a_3(dimless/sqr(dimPressure), 0.005); // [1/kPa^2]
        const dimensionedScalar D0(dimPressure, 1.2); // [kPa]        
        volScalarField D = (pVSatAir - pVAir) / 1000; // [kPa]
        volScalarField r_s = r_sMin_ * (a_1 + qPlantRad) / (a_2 + qPlantRad) * (1 + a_3 * sqr(D - D0));
        // volScalarField r_s = r_a * 0 + r_sMin_;
		// Info << "Dimensions of r_s: " << r_s.dimensions() << endl;
		Info << "Current r_s: min = " << min(r_s).value() << ", max = " << max(r_s).value() << endl;
		Info << "Current D: min = " << min(D).value() << ", max = " << max(D).value() << endl;
		

		// Convective mass transfer coefficient - h_cm [s/m]
		volScalarField h_cm = (1/(r_a + r_s)) * (R_a_/R_v_) * (rho / (p * r_));
		// Info << "Dimensions of h_cm: " << h_cm.dimensions() << endl;
		// Info << "Dimensions of p: " << p.dimensions() << endl;
		Info << "Current p: min = " << min(p).value() << ", max = " << max(p).value() << endl;
		Info << "Current h_cm: min = " << min(h_cm).value() << ", max = " << max(h_cm).value() << endl;
		Info << "Current rho: min = " << min(rho).value() << ", max = " << max(rho).value() << endl;
		
		
        // Vapour mass flux from the leaf (TRANSPIRATION MODEL) - g_vLeaf [kg/(s m^2)]       
		volScalarField g_vLeaf = h_cm * (pVSatLeaf - pVAir);
		volScalarField g_vLeaf_mg = h_cm * (pVSatLeaf - pVAir)*1000;		
		// Info << "Dimensions of g_vLeaf: " << g_vLeaf.dimensions() << endl;
		Info << "Current g_vLeaf_kg (transpiration): min = " << min(g_vLeaf).value() << ", max = " << max(g_vLeaf).value() << endl;
		Info << "Current g_vLeaf_mg (transpiration): min = " << min(g_vLeaf_mg).value() << ", max = " << max(g_vLeaf_mg).value() << endl;
		
		
		// Latent heat flux from the leaf - qPlantLat [W/m^2]
        volScalarField qPlantLat = L_v_ * g_vLeaf;
		// Info << "Dimensions of qPlantLat: " << qPlantLat.dimensions() << endl;
		Info << "Current qPlantLat: min = " << min(qPlantLat).value() << ", max = " << max(qPlantLat).value() << endl;
        
		
		// Convective heat transfer coefficient - h_ch [W/(m^2 K)]
		volScalarField h_ch = (2.0 * rho * Cp0_ ) / r_a;
		// Info << "Dimensions of h_ch: " << h_ch.dimensions() << endl;
		Info << "Current h_ch: min = " << min(h_ch).value() << ", max = " << max(h_ch).value() << endl;
		
		
		// Leaf energy balance - TLeaf [K]
		TLeaf = T + (qPlantRad - qPlantLat) / h_ch;
		

        // reset TLeaf to T where there is not any plant
        TLeaf = (1 - plantFilter(LAD)) * T + plantFilter(LAD) * TLeaf;
		

    } while (!converged(TLeaf) && i++ < maxIter_);

    if (i < maxIter_) {
        Info
            << "The leaf temperature converged within " << i
            << " iterations" << endl;
    } else if (i == maxIter_) {
        WarningInFunction
            << "The leaf temperature did not converge within " << i
            << " iterations" << nl;
    }

    word fieldNameTLeaf = "TLeaf";

    return store(fieldNameTLeaf, tTLeaf);
}

bool Foam::functionObjects::TLeaf::write()
{
    return writeObject("TLeaf");
}


// ************************************************************************* //

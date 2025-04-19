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

#include "ESource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
	defineTypeNameAndDebug(ESource, 0);
	addRemovableToRunTimeSelectionTable
(
    option,
    ESource,
    dictionary
);

}
}

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * //

/**
 * Constructor: Initializes the ESource function object.
 *
 * @param name       Function object name.
 * @param modelType  Type of model.
 * @param dict       Dictionary with user-defined parameters.
 * @param mesh       Mesh reference.
 */
Foam::fv::ESource::ESource
(
    const Foam::word& name,
    const Foam::word& modelType,
    const Foam::dictionary& dict,
    const Foam::fvMesh& mesh

)
:
    fv::cellSetOption(name, modelType, dict, mesh),
	Cp0_("Cp0", dimSpecificHeatCapacity, 0.0),    /**< Specific heat capacity - Cp₀ [J/(kg·K)] */
	L_("L", dimLength, 0.0),                      /**< Characteristic leaf size - L [m] */
	C_("C", sqrt(dimTime)/dimLength, 0.0)         /**< Proportionality factor - C [sqrt(s)/m] */

{
    read(dict);
}


// * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * //

/**
 * Destructor: Destroys the ESource object.
 */
Foam::fv::ESource::~ESource()
{
    Info << "Destructor: ESource" << endl;
}


// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

/**
 * Reads parameters from the dictionary.
 *
 * @param dict   Dictionary containing source term settings.
 * @return True if data is read successfully, false otherwise.
 */

bool Foam::fv::ESource::read(const dictionary& dict)
{
	Info << " 🟡 access to ESource::read()" << endl;
    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }
	
    coeffs_.readEntry("fields", fieldNames_);
    applied_.resize(fieldNames_.size(), false);
	
	Info << " 🟡 return = true from read" << endl;

    // Extract the subdictionary leafProperties
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

    if (leafDict.found("leafProperties"))
    {
        const dictionary& gProps = leafDict.subDict("leafProperties");
		
		Cp0_.readIfPresent(gProps);
        L_.readIfPresent(gProps);
        C_.readIfPresent(gProps);

        Info << "🟡 Values after reading leafProperties:\n"
             << "Cp0 = " << Cp0_.value() << "\n"
			 << "L = " << L_.value() << "\n"
             << "C = " << C_.value() << "\n";
    }
    else
    {
        Info << "\n❌ Subdictionary 'leafProperties' not found in leafProperties.\n";
    }

    return true;
}

/**
 * Adds explicit source term to the temperature equation.
 *
 * @param rho     Density field.
 * @param eqn     Temperature equation matrix.
 * @param fieldi  Index of the field.
 */
void Foam::fv::ESource::addSup
(
    const volScalarField& rho,
	fvMatrix<scalar>& eqn,
    const label fieldi
)
    {
        
        Info << " 🟡 TSource::addSup() addSup() applying source in the scalar field T" << endl;
        
		//const volScalarField& h = eqn.psi();
		const dimensionedScalar timeUnit
    (
        dimTime,
        1
    );
        // Temperature source term - ST [W/(m³·s)]
		const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
		const volScalarField& TLeaf = mesh_.lookupObject<volScalarField>("TLeaf");
		const volScalarField& rho_ = mesh_.lookupObject<volScalarField>("rho");
		const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
		const volScalarField& LAD = mesh_.lookupObject<volScalarField>("LAD");

		/** Heat transfer conductance - gₐH [m/s] */
        const dimensionedScalar Umin(dimVelocity, 0.03);
        const dimensionedScalar Umax(dimVelocity, 1000);
        volScalarField Umag(mag(U));
        Umag.clip(Umin, Umax);
        volScalarField g_aH = (1 / C_) * sqrt(Umag / L_);		
		
		/** Heat transfer resistance - rₐH [s/m] */
        volScalarField r_aH = 1 / g_aH;
		
		/** Convective heat transfer coefficient - hcₕ [W/(m²·K)] */
		volScalarField h_ch = 2 * (rho_ * Cp0_ ) / r_aH;
		
		/** Sensible heat flux from the leaf - qPlantSen [W/m²] */
		volScalarField qPlantSen = h_ch * (TLeaf - T);
		Info << "Current qPlantSen: min = " << min(qPlantSen).value() << ", max = " << max(qPlantSen).value() << ", mean: " << average(qPlantSen).value() << endl;
		
		/** Temperature source term - ST [W/(m³·s)]*/
		volScalarField ST = LAD * qPlantSen / (rho_ * Cp0_);
		
		// Add source term to equation		
		eqn += ST;
		Info << "Current TSource: min = " << min(ST).value() << ", max = " << max(ST).value() << ", mean: " << average(ST).value() << endl;
    }

/**
 * Adds explicit contribution for incompressible flow (needed for fvOptions).
 *
 * @param eqn     Temperature equation matrix.
 * @param fieldi  Index of the field.
 */
void Foam::fv::ESource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)

/**
 * Writes the ESource field to a file.
 *
 * @return True if the write operation is successful.
 */
{
    return this->addSup(volScalarField::null(), eqn, fieldi);
}
// ************************************************************************* //
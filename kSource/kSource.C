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

#include "kSource.H"
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
	defineTypeNameAndDebug(kSource, 0);
	addRemovableToRunTimeSelectionTable
(
    option,
    kSource,
    dictionary
);

}
}

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * //

/**
 * Constructor: Initializes the kSource function object.
 *
 * @param name       Function object name.
 * @param modelType  Type of model.
 * @param dict       Dictionary with user-defined parameters.
 * @param mesh       Mesh reference.
 */

Foam::fv::kSource::kSource
(
    const Foam::word& name,
    const Foam::word& modelType,
    const Foam::dictionary& dict,
    const Foam::fvMesh& mesh

)
:
    fv::cellSetOption(name, modelType, dict, mesh),
	C_d_(0.0),  /**< Drag coefficient Cd - [-] */
    betaP_(0.0), /**< Production coefficient βp - [-] */
	betaD_(0.0)  /**< Dissipation coefficient βd - [-] */

{
    read(dict);
}

// * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * //

/**
 * Destructor: Destroys the kSource object.
 */

Foam::fv::kSource::~kSource()
{
    Info << "Destructor: kSource" << endl;
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

/**
 * Reads parameters from the dictionary.
 *
 * @param dict   Dictionary containing source term settings.
 * @return True if data is read successfully, false otherwise.
 */
 
bool Foam::fv::kSource::read(const dictionary& dict)
{
    Info << " 🔴 access to kSource::read()" << endl;
    if (!fv::cellSetOption::read(dict))
    {
        Info << "return = false from read" << endl;
		return false;
    }

    coeffs_.readEntry("fields", fieldNames_);
    applied_.resize(fieldNames_.size(), false);
	
	Info << " 🔴 return = true from read" << endl;

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

        if (gProps.found("CdCoeffs"))
		{
			const dictionary& CdDict = gProps.subDict("CdCoeffs");
			if (CdDict.found(this->name()))
			{
				CdDict.lookup(this->name()) >> C_d_;
				Info << "C_d from CdCoeffs[" << this->name() << "] = " << C_d_ << endl;
			}
			else
			{
				Info << "CdCoeffs: entry for " << this->name() << " not found. Using default.\n";
			}
		}
		else
		{
			gProps.readIfPresent("C_d", C_d_); // fallback if no CdCoeffs
		}
        gProps.readIfPresent("betaP", betaP_);
        gProps.readIfPresent("betaD", betaD_);

        Info << "🔴 Values after reading leafProperties:\n"
             << "C_d = " << C_d_ << "\n"
             << "betaP = " << betaP_ << "\n"
             << "betaD = " << betaD_ << "\n";
    }
    else
    {
        Info << "\n❌ Subdictionary 'leafProperties' not found in leafProperties.\n";
    }

    return true;
}

/**
 * Adds explicit source term to the turbulent kinetic energy (k) equation.
 *
 * @param rho     Density field.
 * @param eqn     Turbulent kinetic energy equation matrix.
 * @param fieldi  Index of the field.
 */
void Foam::fv::kSource::addSup
(
    const volScalarField& rho,
	fvMatrix<scalar>& eqn,
    const label fieldi
)
    {
        
        Info << "🔴 kSource::addSup() applying source in the scalar field k" << endl;

		const volScalarField& k = eqn.psi();
        const volScalarField& rho_ = mesh_.lookupObject<volScalarField>("rho");
		const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
		const volScalarField& LAD = mesh_.lookupObject<volScalarField>("LAD");
		
		// Turbulent Kinetic Energy, TKE source term - Sk [W/m³] * (1 / rho_)		
		volScalarField Sk_production = rho_ * C_d_ * LAD * betaP_ * pow(mag(U), 3) * (1 / rho_);
		
		fvMatrix<scalar> Sk_dissipation
		(
			- fvm::Sp(rho_ * C_d_ * LAD * betaD_ * mag(U) * (1 / rho_), k)
		);
		
		fvMatrix<scalar> Sk
		(
			Sk_production + Sk_dissipation
		);
		
		// Add source term to equation
		eqn += Sk;
		
		// This is only for printing the values to the console (debugging purposes)
		volScalarField Sk_P = rho_ * C_d_ * LAD * betaP_ * pow(mag(U), 3) * (1 / rho_);
		volScalarField Sk_D = rho_ * C_d_ * LAD * betaD_ * mag(U) * (1 / rho_);

		volScalarField kSource
		(
			IOobject
			(
				"kSource",
				mesh_.time().timeName(),
				mesh_,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			Sk_P - Sk_D * k
		);
		
		Info << "Current kSource: min = " << min(mag(kSource)).value() << ", max = " << max(mag(kSource)).value() << ", mean = " << average(mag(kSource)).value() << endl;
	}

/**
 * Adds explicit contribution for incompressible flow (needed for fvOptions).
 *
 * @param eqn     Turbulent kinetic energy equation matrix.
 * @param fieldi  Index of the field.
 */
void Foam::fv::kSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)

/**
 * Writes the kSource field to a file.
 *
 * @return True if the write operation is successful.
 */
{
    return this->addSup(volScalarField::null(), eqn, fieldi);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
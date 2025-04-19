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

#include "epsilonSource.H"
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
	defineTypeNameAndDebug(epsilonSource, 0);
	addRemovableToRunTimeSelectionTable
(
    option,
    epsilonSource,
    dictionary
);

}
}

// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * //

/**
 * Constructor: Initializes the epsilonSource function object.
 *
 * @param name       Function object name.
 * @param modelType  Type of model.
 * @param dict       Dictionary with user-defined parameters.
 * @param mesh       Mesh reference.
 */
Foam::fv::epsilonSource::epsilonSource
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
	betaD_(0.0), /**< Dissipation coefficient βd - [-] */
	C4_(0.0),    /**< Model coefficient 1 C₄ - [-] */
	C5_(0.0)     /**< Model coefficient 2 C₅ - [-] */

{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

/**
 * Destructor: Destroys the epsilonSource object.
 */
Foam::fv::epsilonSource::~epsilonSource()
{
    Info << "Destructor: epsilonSource" << endl;
}

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

/**
 * Reads parameters from the dictionary.
 *
 * @param dict   Dictionary containing source term settings.
 * @return True if data is read successfully, false otherwise.
 */

bool Foam::fv::epsilonSource::read(const dictionary& dict)
{
	Info << " 🟢 access to εSource::read()" << endl;
    if (!fv::cellSetOption::read(dict))
    {
        return false;
    }
	
    coeffs_.readEntry("fields", fieldNames_);
    applied_.resize(fieldNames_.size(), false);
	
	Info << " 🟢 return = true from read" << endl;

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
		gProps.readIfPresent("C4", C4_);
		gProps.readIfPresent("C5", C5_);

        Info << "🟢 Values after reading leafProperties:\n"
             << "C_d = " << C_d_ << "\n"
             << "betaP = " << betaP_ << "\n"
             << "betaD = " << betaD_ << "\n"
			 << "C4 = " << C4_ << "\n"
			 << "C5 = " << C5_ << "\n";
    }
    else
    {
        Info << "\n❌ Subdictionary 'leafProperties' not found in leafProperties.\n";
    }

    return true;
}

/**
 * Adds explicit source term to the dissipation rate (epsilon) equation.
 *
 * @param rho     Density field.
 * @param eqn     Dissipation rate equation matrix.
 * @param fieldi  Index of the field.
 */
void Foam::fv::epsilonSource::addSup
(
    const volScalarField& rho,
	fvMatrix<scalar>& eqn,
    const label fieldi
)
    {
        
        Info << "🟢 εSource::addSup() addSup() applying source in the scalar field ε" << endl;
        
		const volScalarField& epsilon = eqn.psi();
        const volScalarField& k = mesh_.lookupObject<volScalarField>("k");
		const volScalarField& rho_ = mesh_.lookupObject<volScalarField>("rho");
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
		const volScalarField& LAD = mesh_.lookupObject<volScalarField>("LAD");
		
		// Dissipation Rate, TDR source term - Sε [W/(m³·s)] * (1 / rho_)		
		fvMatrix<scalar> Sε_production
        (
			fvm::Sp(rho_ * C_d_ * LAD * betaP_ * C4_ * pow(mag(U),3)/k * (1 / rho_), epsilon)
		);
		
		fvMatrix<scalar> Sε_dissipation
        (
			fvm::Sp(- rho_ *C_d_ * LAD * betaD_ * C5_ * mag(U) * (1 / rho_), epsilon)
		);
		
        fvMatrix<scalar> Sε_total
        (
			Sε_production + Sε_dissipation
		);
		
		// Add source term to equation	
        eqn += Sε_total;

		// This is only for printing the values to the console (debugging purposes)
		volScalarField Sε_P = rho_ * C_d_ * LAD * betaP_ * C4_ * pow(mag(U),3)/k * (1 / rho_);
		volScalarField Sε_D = rho_ *C_d_ * LAD * betaD_ * C5_ * mag(U) * (1 / rho_);

		volScalarField εSource
		(
			IOobject
			(
				"εSource",
				mesh_.time().timeName(),
				mesh_,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			Sε_P * epsilon - Sε_D * epsilon
		);
		
		Info << "Current εSource: min = " << min(mag(εSource)).value() << ", max = " << max(mag(εSource)).value() << ", mean = " << average(mag(εSource)).value() << endl;
	}

/**
 * Adds explicit contribution for incompressible flow (needed for fvOptions).
 *
 * @param eqn     Dissipation rate equation matrix.
 * @param fieldi  Index of the field.
 */
void Foam::fv::epsilonSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)

/**
 * Writes the epsilonSource field to a file.
 *
 * @return True if the write operation is successful.
 */
{
    return this->addSup(volScalarField::null(), eqn, fieldi);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
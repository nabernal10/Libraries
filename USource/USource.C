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

#include "USource.H"
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
	defineTypeNameAndDebug(USource, 0);
	addRemovableToRunTimeSelectionTable
(
    option,
    USource,
    dictionary
);

}
}

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

/**
 * Constructor: Initializes the USource function object.
 *
 * @param name       Function object name.
 * @param modelType  Type of model.
 * @param dict       Dictionary with user-defined parameters.
 * @param mesh       Mesh reference.
 */

Foam::fv::USource::USource
(
    const Foam::word& name,
    const Foam::word& modelType,
    const Foam::dictionary& dict,
    const Foam::fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    C_d_(0.0)  /**< Drag coefficient Cd - [-] */
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

/**
 * Destructor: Destroys the USource object.
 */

Foam::fv::USource::~USource()
{
    Info << "Destructor: USource" << endl;
}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

/**
 * Reads parameters from the dictionary.
 *
 * @param dict   Dictionary containing source term settings.
 * @return True if data is read successfully, false otherwise.
 */

bool Foam::fv::USource::read(const dictionary& dict)
{
	Info << " 🔵 access to USource::read()" << endl;
    if (!fv::cellSetOption::read(dict))
    {
    	Info << "return = false from read" << endl;
        return false;
    }

    coeffs_.readEntry("fields", fieldNames_);
    applied_.resize(fieldNames_.size(), false);
    
    Info << " 🔵 return = true from read" << endl;

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

        gProps.readIfPresent("C_d", C_d_);

        Info << "🔵 Values after reading leafProperties:\n"
             << "C_d = " << C_d_ << "\n";
    }
    else
    {
        Info << "\n❌ Subdictionary 'leafProperties' not found in leafProperties.\n";
    }

    return true;
}

/**
 * Adds explicit source term to the momentum equation.
 *
 * @param rho     Density field.
 * @param eqn     Momentum equation matrix.
 * @param fieldi  Index of the field.
 */
void Foam::fv::USource::addSup
(
    const volScalarField& rho,
	fvMatrix<vector>& eqn,
    const label fieldi
)
{
		Info << "🔵 USource::addSup() applying source in the scalar field U" << endl;
		
		const volScalarField& rho_ = mesh_.lookupObject<volScalarField>("rho");		
		const volVectorField& U = eqn.psi();
		const volScalarField& LAD = mesh_.lookupObject<volScalarField>("LAD");

		// Source term for momentum equation [N/m³] * (1 / rho_)
		fvMatrix<vector> SU
		(
			fvm::Sp(rho_ * C_d_ * LAD * mag(U) * (1 / rho_), U)
		);
		
		// Add source term to equation		
		eqn -= SU;
}

/**
 * Adds explicit contribution for incompressible flow (needed for fvOptions).
 *
 * @param eqn     Momentum equation matrix.
 * @param fieldi  Index of the field.
 */
 
void Foam::fv::USource::addSup //(Para que fvOptions detecte USource)
(
    fvMatrix<vector>& eqn,
    const label fieldi
)

/**
 * Writes the USource field to a file.
 *
 * @return True if the write operation is successful.
 */
{
    return this->addSup(volScalarField::null(), eqn, fieldi);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
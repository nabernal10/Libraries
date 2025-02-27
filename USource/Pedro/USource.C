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

#include "USource.H"
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

defineTypeNameAndDebug(USource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    USource,
    dictionary
);

}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::USource::USource
(
    const Foam::word& name,
    const Foam::word& modelType,
    const Foam::dictionary& dict,
    const Foam::fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    C_d(1.2)  // C_d(0.2)
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::USource::~USource()
{
    Info << "Destructor: USource" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::fv::USource::read(const dictionary& dict)
{
		Info << "access to USource::read()" << endl;
    if (!fv::cellSetOption::read(dict))
    {
    		Info << "return = false from read" << endl;
        return false;
    }

    coeffs_.readEntry("fields", fieldNames_);
    applied_.resize(fieldNames_.size(), false);
    
    Info << "return = true from read" << endl;

    return true;
}


void Foam::fv::USource::addSup
(
    const volScalarField& rho,
	fvMatrix<vector>& eqn,
    const label fieldi
)
{
    Info<<"mySource::addSup()\n";
    scalar rho_ = 1; // constant rho value
		const volVectorField& U = eqn.psi();
		const volScalarField& LAD = mesh_.lookupObject<volScalarField>("LAD");
		volScalarField plantCells = pos(LAD); // Devuelve 1 donde LAD > 0, 0 en otro caso
		scalar numPlantCells = gSum(plantCells);

		Info << " 🔹 Número de celdas con LAD > 0: " << numPlantCells << endl;

		// Momentum source term - Su [N/m^3]

		fvMatrix<vector> SU
    (
        fvm::Sp(rho_*C_d*LAD*mag(U), U)
		);
				
    eqn -= SU;	
    }

void Foam::fv::USource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    return this->addSup(volScalarField::null(), eqn, fieldi);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //

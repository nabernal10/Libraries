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

#include "epsilonSource.H"
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

defineTypeNameAndDebug(epsilonSource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    epsilonSource,
    dictionary
);

}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::epsilonSource::epsilonSource
(
    const Foam::word& name,
    const Foam::word& modelType,
    const Foam::dictionary& dict,
    const Foam::fvMesh& mesh

)
:
    fv::cellSetOption(name, modelType, dict, mesh),
	C_d(0.7),
    betaP_(1.0),
	betaD_(5.1),
	C4_(0.9),
	C5_(0.9)

{
    read(dict);           // <<-- Add
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::epsilonSource::~epsilonSource()
{
    Info << "Destructor: epsilonSource" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
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

    return true;
}


void Foam::fv::epsilonSource::addSup
(
    const volScalarField& rho,
	fvMatrix<scalar>& eqn,
    const label fieldi
)
    {
        
        Info << "🟢 εSource::addSup() aplicando fuente en el campo escalar ε" << endl;
        
		const volScalarField& epsilon = eqn.psi();
        const volScalarField& k = mesh_.lookupObject<volScalarField>("k");
		const volScalarField& rho_ = mesh_.lookupObject<volScalarField>("rho");
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
		const volScalarField& LAD = mesh_.lookupObject<volScalarField>("LAD");
		
		// Dissipation Rate, TDR source term - Sε [W/(m^3 s)] * (1 / rho_)
		
		fvMatrix<scalar> Sε_production
        (
			fvm::Sp(rho_ * C_d * LAD * betaP_ * C4_ * pow(mag(U),3)/k * (1 / rho_), epsilon)
		);
		
		fvMatrix<scalar> Sε_dissipation
        (
			fvm::Sp(- rho_ *C_d * LAD * betaD_ * C5_ * mag(U) * (1 / rho_), epsilon)
		);
		
        fvMatrix<scalar> Sε
        (
			Sε_production + Sε_dissipation
		);
		
        eqn += Sε;
    }

void Foam::fv::epsilonSource::addSup //(Para que fvOptions detecte εSource)
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    return this->addSup(volScalarField::null(), eqn, fieldi);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
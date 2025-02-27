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

Class
    Foam::functionObjects::GLeaf

Description
    Calculates GLeaf from G and T.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    GLeaf1
    {
        // Mandatory entries
        type            GLeaf;
        libs            (fieldFunctionObjects);

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                                | Type   | Reqd | Deflt
      type      | Type name: GLeaf                           | word   | yes  | -
      libs      | Library name: fieldFunctionObjects         | word   | yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

SourceFiles
    GLeaf.C

\*---------------------------------------------------------------------------*/

#ifndef GLeaf_H
#define GLeaf_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class GLeaf Declaration
\*---------------------------------------------------------------------------*/

class GLeaf
:
    public fvMeshFunctionObject
{

public:

    //- Runtime type information
    TypeName("GLeaf");


    // Constructors

        //- Construct from Time and dictionary
        GLeaf
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~GLeaf() = default;


    // Member Functions

        //- Read the data needed for the GLeaf calculation
        virtual bool read(const dictionary&);

        //- Calculate the Incident radiation on the (GLeaf)
        virtual bool execute();

        //- Write the GLeaf field
        virtual bool write();
		
private:
	// - Specific heat capacity [J/(kg K)]
    dimensionedScalar Cp0_;
	
	// - Absorptance (α_leaf): fraction of radiation absorbed.
    dimensionedScalar alphaLeaf_;
	
	// - Reflectance or albedo (ρ_leaf): fraction of radiation reflected.
    dimensionedScalar rhoLeaf_;
	
	// - Extinction coefficient.
	dimensionedScalar kappa_;
	
	// - Extinction coefficient.
	dimensionedScalar taoLeaf_;
	
	// - Emissivity (ε_leaf): ability to emit thermal radiation.
    dimensionedScalar epsilonLeaf_;
	
	// Solar elevation angle - β [degrees]
	dimensionedScalar betaDeg_;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

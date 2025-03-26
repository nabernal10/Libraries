/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | www.openfoam.com
    \\  /    A nd           |
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

    OpenFOAM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
   Class: Foam::functionObjects::GLeaf

   Description:
       Computes the net radiation balance for a leaf surface.
       It calculates incident shortwave and longwave radiation.

   Usage:
       To use this function object in `system/controlDict.functions`:
       \verbatim
       GLeaf1
       {
           type            GLeaf;
           libs            (fieldFunctionObjects);

           Cp0             ...;    // Specific heat capacity [J/(kg·K)]
           alphaLeaf       ...;    // Absorptance coefficient [-]
           rhoLeaf         ...;    // Reflectance coefficient (albedo) [-]
           taoLeaf         ...;    // Transmittance coefficient [-]
           kappaLeafSW     ...;    // Shortwave extinction coefficient [-]
		   epsilonLeaf     ...;    // Leaf emissivity coefficient [-]
		   epsilonSky      ...;    // Sky emissivity coefficient [-]
 		   betaDeg         ...;    // Solar elevation angle [degrees]
		   C_lw            ...;    // Longwave radiation empirical constant [-]
		   H               ...;    // Tree height - H [m]
		   h               ...;    // Canopy base height - h [m]
       }
       \endverbatim

   Property Table:
   | Property      | Description                                 | Type   | Required |
   |---------------|---------------------------------------------|--------|----------|
   | type          | Type name: GLeaf                            | word   | Yes      |
   | libs          | Library name: fieldFunctionObjects          | word   | Yes      |
   | Cp0           | Specific heat capacity of air               | scalar | No       |
   | alphaLeaf     | Leaf absorptance (fraction absorbed)        | scalar | No       |
   | rhoLeaf       | Leaf reflectance (albedo)                   | scalar | No       |
   | taoLeaf       | Leaf transmittance (fraction transmitted)   | scalar | No       |
   | kappaLeafSW   | Shortwave radiation extinction              | scalar | No       |
   | epsilonLeaf   | Leaf emissivity                             | scalar | No       |
   | epsilonSky    | Sky emissivity                              | scalar | No       |
   | betaDeg       | Solar elevation angle                       | scalar | No       |
   | C_lw          | Longwave radiation empirical constant       | scalar | No       |
   | H             | Tree height                                 | scalar | No       |
   | h             | Canopy base height                          | scalar | No       |

   See also:
       - \link functionObject.H \endlink

   SourceFiles:
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

        /**
         * Constructor: Initializes the GLeaf function object
         *
         * @param name       Function object name
         * @param runTime    Simulation time
         * @param dict       Dictionary with user-defined parameters
         */
        GLeaf
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~GLeaf() = default;

    // Member Functions

        /**
         * Reads parameters from the dictionary
         *
         * @param dict   Dictionary containing the function object settings
         * @return       True if read successfully, false otherwise
         */
        virtual bool read(const dictionary&);

        /**
         * Computes the incident radiation on the leaf (GLeaf)
         *
         * @return True if execution is successful
         */
        virtual bool execute();

        /**
         * Writes the GLeaf field to a file
         *
         * @return True if write operation is successful
         */
        virtual bool write();
		
private:

    // Leaf properties

    /** Specific heat capacity - Cp₀ [J/(kg·K)] */
    dimensionedScalar Cp0_;
	
    /** Absorptance coefficient - α_leaf [-]: Fraction of radiation absorbed */
    dimensionedScalar alphaLeaf_;
	
    /** Reflectance coefficient or albedo - ρ_leaf [-]: Fraction of radiation reflected */
    dimensionedScalar rhoLeaf_;
	
    /** Transmittance coefficient - τ_leaf [-]: Fraction of radiation transmitted */
    dimensionedScalar taoLeaf_;
	
	/** Extinction coefficient - k_leafSW [-] */
    dimensionedScalar kappaLeaf_SW_;
	
    /** Emissivity coefficient - ε_leaf [-]: Leaf ability to emit thermal radiation */
    dimensionedScalar epsilonLeaf_;
	
	/** Emissivity coefficient - ε_sky [-]: Sky ability to emit thermal radiation */
    dimensionedScalar epsilonSky_;
	
    /** Solar elevation angle - β [degrees]: Ability to determine the impact of solar angle on radiation balance */
    dimensionedScalar betaDeg_;
	
    /** Long-wave radiation empirical constant - C_lw [-]: Empirical constant for quantifying the net absorption of long-wave radiation */
    dimensionedScalar C_lw_;
	
	/** Tree height - H [m] */
    dimensionedScalar H_;
	
	/** Canopy base height - h [m] */
    dimensionedScalar h_;

    /**
     * Computes a binary field based on leaf area density.
     *
     * @param leafAreaDensity Leaf area density field.
     * @return Field with values of 1 where a tree is present, 0 otherwise.
     */
    volScalarField plantFilter(const volScalarField&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
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
   Class: Foam::functionObjects::TLeaf

   Description:
       Computes leaf temperature (TLeaf) based on energy balance equations.

   Usage:
       To use this function object in `system/controlDict.functions`:
       \verbatim
       TLeaf1
       {
           type            TLeaf;
           libs            (fieldFunctionObjects);

           Cp0             ...;    // Specific heat capacity [J/(kg·K)]
           L               ...;    // Characteristic leaf size [m]
           C               ...;    // Proportionality factor [sqrt(s)/m]
           r_sMin          ...;    // Minimum stomatal resistance [s/m]
           R_a             ...;    // Gas constant of dry air [J/(kg·K)]
           R_v             ...;    // Gas constant of vapor [J/(kg·K)]
           L_v             ...;    // Latent heat of vaporization [J/kg]
		   tolerance       ...;    // Convergence tolerance for TLeaf
           maxIter         ...;    // Maximum iterations for temperature solution
       }
       \endverbatim

   Property Table:
   | Property   | Description                                | Type   | Required |
   |------------|--------------------------------------------|--------|----------|
   | type       | Type name: TLeaf                           | word   | Yes      |
   | libs       | Library name: fieldFunctionObjects         | word   | Yes      |
   | Cp0        | Specific heat capacity of air              | scalar | No       |
   | L          | Characteristic leaf size                   | scalar | No       |
   | C          | Proportionality factor                     | scalar | No       |
   | r_sMin     | Minimum stomatal resistance                | scalar | No       |
   | R_a        | Gas constant of dry air                    | scalar | No       |
   | R_v        | Gas constant of vapor                      | scalar | No       |
   | L_v        | Latent heat of vaporization                | scalar | No       |  
   | tolerance  | Convergence tolerance for TLeaf            | scalar | No       |
   | maxIter    | Maximum iterations for temperature solution| scalar | No       |

   See also:
       - \link functionObject.H \endlink

   SourceFiles:
       TLeaf.C

\*---------------------------------------------------------------------------*/

#ifndef TLeaf_H
#define TLeaf_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class TLeaf Declaration
\*---------------------------------------------------------------------------*/

class TLeaf
:
    public fvMeshFunctionObject
{
    // Private Data

    /** Specific heat capacity - Cp₀ [J/(kg·K)] */
    dimensionedScalar Cp0_;

    /** Characteristic leaf size - L [m] */
    dimensionedScalar L_;

    /** Proportionality factor - C [sqrt(s)/m] */
    dimensionedScalar C_;

    /** Minimum stomatal resistance - rₛ,Min [s/m] */
    dimensionedScalar r_sMin_;

    /** Gas constant of dry air - Rₐ [J/(kg·K)] */
    dimensionedScalar R_a_;

    /** Gas constant of vapor - Rᵥ [J/(kg·K)] */
    dimensionedScalar R_v_;

    /** Latent heat of vaporization - Lᵥ [J/kg] */
    dimensionedScalar L_v_;
		
    /** Tolerance for iterative solution - tolerance [-] */
    scalar tolerance_;

    /** Maximum iterations allowed for solving TLeaf - maxIter [-] */
    int maxIter_;


    // Private Member Functions

    /**
     * Checks if the leaf temperature iteration has converged.
     *
     * @param TLeaf Field of computed leaf temperature.
     * @return True if the iteration has converged.
     */
    bool converged(const volScalarField&) const;

    /**
     * Computes a binary field based on leaf area density.
     *
     * @param leafAreaDensity Leaf area density field.
     * @return Field with values of 1 where a tree is present, 0 otherwise.
     */
    volScalarField plantFilter(const volScalarField&) const;

public:

    //- Runtime type information
    TypeName("TLeaf");

    // Constructors

    /**
     * Constructor: Initializes the TLeaf function object.
     *
     * @param name       Function object name.
     * @param runTime    Simulation time.
     * @param dict       Dictionary with user-defined parameters.
     */
    TLeaf
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Destructor
    virtual ~TLeaf() = default;

    // Member Functions

    /**
     * Reads parameters from the dictionary.
     *
     * @param dict   Dictionary containing function object settings.
     * @return       True if read successfully, false otherwise.
     */
    virtual bool read(const dictionary&);

    /**
     * Computes the leaf temperature based on energy balance.
     *
     * @return True if execution is successful.
     */
    virtual bool execute();

    /**
     * Writes the TLeaf field to a file.
     *
     * @return True if the write operation is successful.
     */
    virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
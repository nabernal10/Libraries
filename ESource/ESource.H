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
    Class: Foam::fv::ESource

    Description:
        Computes a source term contribution to the energy equation.
        The source term is defined by user-provided coefficients and
        added explicitly to the energy equation.

    Usage:
        To use this function object in `system/controlDict.functions`:
        \verbatim
        ESource1
        {
            type            ESource;
            active          yes;

            ESourceCoeffs
            {
                fields      (T);
                Cp0         ...;    // Specific heat capacity [J/(kg·K)]
                L           ...;    // Characteristic leaf size [m]
                C           ...;    // Proportionality factor [sqrt(s)/m]
            }
        }
        \endverbatim

    Property Table:
    | Property   | Description                                | Type   | Required |
    |------------|--------------------------------------------|--------|----------|
    | type       | Type name: ESource                         | word   | Yes      |
    | fields     | Fields to apply the source term            | list   | Yes      |
    | Cp0        | Specific heat capacity                     | scalar | No       |
    | L          | Characteristic length                      | scalar | No       |
    | C          | Proportionality factor                     | scalar | No       |
    | active     | Whether the source term is active          | bool   | Yes      |

    See also:
        - \link cellSetOption.H \endlink

    SourceFiles:
        ESource.C

\*---------------------------------------------------------------------------*/

#ifndef ESource_H
#define ESource_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class ESource Declaration
\*---------------------------------------------------------------------------*/

class ESource
:
    public cellSetOption
{
    // Private Data

    /** Specific heat capacity [J/(kg·K)] */
    dimensionedScalar Cp0_;

    /** Characteristic length [m] */
    dimensionedScalar L_;

    /** Proportionality factor [sqrt(s)/m] */
    dimensionedScalar C_;

    // Private Member Functions   

public:

 
    //- Runtime type information
    TypeName("ESource");


    // Constructors

    /**
     * Constructor: Initializes the ESource object.
     *
     * @param name       Name of the source term.
     * @param modelType  Type of model.
     * @param dict       Dictionary containing model parameters.
     * @param mesh       Mesh reference.
     */
        ESource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~ESource();


    // Member Functions

    /**
     * Read data from dictionary.
     *
     * @param dict  Dictionary containing source term settings.
     * @return True if data is read successfully, false otherwise.
     */
        virtual bool read(const dictionary&);
		

    /**
     * Add explicit source term to the energy equation.
     *
     * @param rho     Density field.
     * @param eqn     Energy equation matrix.
     * @param fieldi  Index of the field.
     */
        virtual void addSup
        (
            const volScalarField& rho,
	        fvMatrix<scalar>& eqn,
            const label fieldi
        );
		
    /**
     * Add explicit contribution to the incompressible equation.
     *
     * @param eqn     Energy equation matrix.
     * @param fieldi  Index of the field.
     */
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
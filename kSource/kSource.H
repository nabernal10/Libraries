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
    Class: Foam::fv::kSource

    Description:
        Adds a source term to the turbulence kinetic energy equation (k).
        The source term is defined based on user-defined coefficients.

    Usage:
        To use this function object in `system/controlDict.functions`:
        \verbatim
        kSource1
        {
            type            kSource;
            active          yes;

            kSourceCoeffs
            {
                fields      (k);
                C_d         ...;    // Drag coefficient [-]
                betaP       ...;    // Production coefficient [-]
                betaD       ...;    // Dissipation coefficient [-]
            }
        }
        \endverbatim

    Property Table:
    | Property   | Description                                | Type   | Required |
    |------------|--------------------------------------------|--------|----------|
    | type       | Type name: kSource                         | word   | Yes      |
    | fields     | Fields to apply the source term            | list   | Yes      |
    | C_d        | Drag coefficient                           | scalar | No       |
    | betaP      | Production coefficient                     | scalar | No       |
    | betaD      | Dissipation coefficient                    | scalar | No       |
    | active     | Whether the source term is active          | bool   | Yes      |

    See also:
        - \link cellSetOption.H \endlink

    SourceFiles:
        kSource.C

\*---------------------------------------------------------------------------*/

#ifndef kSource_H
#define kSource_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class kSource Declaration
\*---------------------------------------------------------------------------*/

class kSource
:
    public cellSetOption
{
	// Private Data

		/** Drag coefficient [-] */
		scalar C_d_;

		/** Production coefficient [-] */
		scalar betaP_;

		/** Dissipation coefficient [-] */
		scalar betaD_;

    // Private Member Functions  

public:

 
    //- Runtime type information
    TypeName("kSource");


    // Constructors

    /**
     * Constructor: Initializes the kSource object.
     *
     * @param name       Name of the source term.
     * @param modelType  Type of model.
     * @param dict       Dictionary containing model parameters.
     * @param mesh       Mesh reference.
     */
        kSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~kSource();


    // Member Functions

    /**
     * Read data from dictionary.
     *
     * @param dict  Dictionary containing source term settings.
     * @return True if data is read successfully, false otherwise.
     */
        virtual bool read(const dictionary&);
		
    /**
     * Add explicit source term to the turbulence kinetic energy equation (k).
     *
     * @param rho     Density field.
     * @param eqn     Momentum equation matrix.
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
     * @param eqn     Turbulence kinetic energy equation.
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
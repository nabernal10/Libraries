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
    Class: Foam::fv::epsilonSource

    Description:
        Adds a source term to the turbulence dissipation equation (epsilon).
        The source term is defined based on user-defined coefficients.

    Usage:
        To use this function object in `system/controlDict.functions`:
        \verbatim
        epsilonSource1
        {
            type            epsilonSource;
            active          yes;

            epsilonSourceCoeffs
            {
                fields      (epsilon);
                C_d         ...;    // Drag coefficient [-]
                betaP       ...;    // Production coefficient [-]
                betaD       ...;    // Dissipation coefficient [-]
                C4          ...;    // Model coefficient 1 [-]
                C5          ...;    // Model coefficient 2 [-]
            }
        }
        \endverbatim

    Property Table:
    | Property   | Description                                | Type   | Required |
    |------------|--------------------------------------------|--------|----------|
    | type       | Type name: epsilonSource                   | word   | Yes      |
    | fields     | Fields to apply the source term            | list   | Yes      |
    | C_d        | Drag coefficient                           | scalar | No       |
    | betaP      | Production coefficient                     | scalar | No       |
    | betaD      | Dissipation coefficient                    | scalar | No       |
    | C4         | Model coefficient 1                        | scalar | No       |
    | C5         | Model coefficient 2                        | scalar | No       |
    | active     | Whether the source term is active          | bool   | Yes      |

    See also:
        - \link cellSetOption.H \endlink

    SourceFiles:
        epsilonSource.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonSource_H
#define epsilonSource_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class epsilonSource Declaration
\*---------------------------------------------------------------------------*/

class epsilonSource
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

    /** Model coefficient 1 [-] */
    scalar C4_;

    /** Model coefficient 2 [-] */
    scalar C5_;
	
	// Private Member Functions        

public:

 
    //- Runtime type information
    TypeName("epsilonSource");


    // Constructors

    /**
     * Constructor: Initializes the epsilonSource object.
     *
     * @param name       Name of the source term.
     * @param modelType  Type of model.
     * @param dict       Dictionary containing model parameters.
     * @param mesh       Mesh reference.
     */
        epsilonSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~epsilonSource();


    // Member Functions

    /**
     * Read data from dictionary.
     *
     * @param dict  Dictionary containing source term settings.
     * @return True if data is read successfully, false otherwise.
     */
        virtual bool read(const dictionary&);
		
    /**
     * Add explicit source term to the dissipation equation.
     *
     * @param rho     Density field.
     * @param eqn     Dissipation equation matrix.
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
     * @param eqn     Dissipation equation matrix.
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
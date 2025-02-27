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

Description
    Template for use with dynamic code generation of a source.
    The hook functions take the following arguments:

    \verbatim
    codeCorrect
    (
        GeometricField<Type, fvPatchField, volMesh>& fld
    )

    codeAddSup
    (
        fvMatrix<Type}>& eqn,
        const label fieldi
    )

    codeConstrain
    (
        fvMatrix<Type>& eqn,
        const label fieldi
    )
    \endverbatim

    where :
        fieldi is the index in the fields entry
        eqn is the fvMatrix

    \verbatim
    energySource
    {
        type            scalarCodedSource;

        active          true;
        selectionMode   all;

        scalarCodedSourceCoeffs
        {
            fields  (h);
            name    sourceTime;

            codeInclude
            #{
            #};

            codeCorrect
            #{
                Pout<< "**codeCorrect**" << endl;
            #};

            codeAddSup
            #{
                const Time& time = mesh().time();
                const scalarField& V = mesh_.V();
                scalarField& heSource = eqn.source();
                heSource -= 0.1*sqr(time.value())*V;
            #};

            codeConstrain
            #{
                Pout<< "**codeConstrain**" << endl;
            #};
        }

        sourceTimeCoeffs
        {
            // Dummy entry
        }
    }
    \endverbatim

SourceFiles
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
                         A templated CodedFvOption
\*---------------------------------------------------------------------------*/

class epsilonSource
:
    public cellSetOption
{
    // Private Data
		// Cofficients
		scalar C_d;
        scalar betaP_;
	    scalar betaD_;
		scalar C4_;
	    scalar C5_;
	// Private Member Functions

        

public:

 
    //- Runtime type information
    TypeName("epsilonSource");


    // Constructors

        //- Construct from patch and internal field
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

        //- Read the data needed for the source calculation
        virtual bool read(const dictionary&);
		
        //- Explicit/implicit matrix contributions
        virtual void addSup
        (
            const volScalarField& rho,
	        fvMatrix<scalar>& eqn,
            const label fieldi
        );		
		
		//- Agregar contribución explícita para flujo incompresible (necesario para fvOptions)
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


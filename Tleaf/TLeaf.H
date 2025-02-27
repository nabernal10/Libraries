/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Class
    Foam::functionObjects::TLeaf

Description
    Calculates TLeaf using a leaf energy balance given incident radiation and tree parameters.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    TLeaf1
    {
        // Mandatory entries
        type            TLeaf;
        libs            (fieldFunctionObjects);

        // Optional entries
        Cp0             <scalar>;
        rho0            <scalar>;
		r               <scalar>;
        p0              <scalar>;
        L               <scalar>;
        C               <scalar>;
        r_sMin          <scalar>;
        R_a             <scalar>;
        R_v             <scalar>;
        L_v             <scalar>;
        G0              <scalar>;
        tolerance       <scalar>;
        maxIter         <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                                | Type   | Reqd | Deflt
      type      | Type name: TLeaf                           | word   | yes  | -
      libs      | Library name: fieldFunctionObjects         | word   | yes  | -
      Cp0       | Specific heat capacity [J/(kg K)]          | scalar | no   | 1003.5
      rho0      | Air density [kg/m^3]                       | scalar | no   | 1.225
      r         | Constant to adjust pressure units [kg/m^3] | scalar | no   | 1
	  p0        | Absolute pressure [Pa]                     | scalar | no   | 101325
      L         | Characteristic leaf size [m]               | scalar | no   | 0.1
      C         | Proportionality factor [sqrt(s)/m]         | scalar | no   | 130
      r_sMin    | Minimum stomatal resistance [s/m]          | scalar | no   | 400
      R_a       | Gas constant of dry air [J/(kg K)]         | scalar | no   | 287.042
      R_v       | Gas constant of vapour [J/(kg K)]          | scalar | no   | 461.524
      L_v       | Latent heat of vapourization [J/kg]        | scalar | no   | 2.5e+6
      L_v       | Reference radiation [W/m^2]                | scalar | no   | 500
      tolerance | Residual control for the leaf temperature  | scalar | no   | 1e-4
      maxIter   | Maximum number of iterations allowed       | scalar | no   | 100
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

SourceFiles
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

        //- Specific heat capacity [J/(kg K)]
        dimensionedScalar Cp0_;

        //- Density of air [kg/m^3]
        dimensionedScalar rho0_;
		
		//- Constant to adjust pressure units [kg/m^3]
        dimensionedScalar r_;

        //- Absolute pressure [Pa]
        dimensionedScalar p0;

        //- Characteristic leaf size [m]
        dimensionedScalar L_;

        //- Proportionality factor [sqrt(s)/m]
        dimensionedScalar C_;

        //- Minimum stomatal resistance [s/m]
        dimensionedScalar r_sMin_;

        //- Gas constant of dry air [J/(kg K)]
        dimensionedScalar R_a_;

        //- Gas constant of vapour [J/(kg K)]
        dimensionedScalar R_v_;

        //- Latent heat of vapourization [J/kg]
        dimensionedScalar L_v_;

        //- Reference radiation [W/m^2]
        dimensionedScalar G0_;

        //- Gravity [m/s^2]
        dimensionedScalar g;
		
		//- Tolerance criteria for iterative process to find TLeaf
        scalar tolerance_;

        //- Maximum number of correctors for leaf temperature
        int maxIter_;


    // Private Member Functions

        //- Return true if the leaf temperature iteration has converged
        bool converged(const volScalarField&) const;

        //- Return a field based on leafAreaDensity that is 1 where there is a tree, 0 where there isn't
        volScalarField plantFilter(const volScalarField&) const;

public:

    //- Runtime type information
    TypeName("TLeaf");


    // Constructors

        //- Construct from Time and dictionary
        TLeaf
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~TLeaf() = default;


    // Member Functions

        //- Read the data needed for the TLeaf calculation
        virtual bool read(const dictionary&);

        //- Calculate the temperature of the leaf (TLeaf)
        virtual bool execute();

        //- Write the TLeaf field
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

// Solve the temperature equation

Obtener referencias a los campos en `createFields.H`
const volScalarField& LAD = mesh.lookupObject<volScalarField>("LAD");
const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");
const volScalarField& TLeaf = mesh.lookupObject<volScalarField>("TLeaf");

// Definir constantes relevantes
dimensionedScalar Cp0("Cp0", dimSpecificHeatCapacity, 1003.5);  // Calor especifico del aire [J/(kg K)]
dimensionedScalar L("L", dimLength, 0.1);  // Tamano caracteristico de la hoja [m]
dimensionedScalar C("C", sqrt(dimTime)/dimLength, 130); // Factor de proporcionalidad [s^0.5/m]

// Calcular la resistencia aerodinamica - r_a [s/m]
const dimensionedScalar Umin(dimVelocity, 0.5);
const dimensionedScalar Umax(dimVelocity, 1000);
volScalarField Umag(mag(U));
Umag.clip(Umin, Umax);
volScalarField r_a = C * sqrt(L / Umag);

// Coeficiente de transferencia de calor por conveccion - h_ch [W/(m^2 K)]
volScalarField h_ch = (2.0 * rho * Cp0) / r_a;

// Flujo de calor sensible desde la hoja - qPlantSen [W/m^3]
volScalarField qPlantSen = h_ch * (TLeaf - T);

// Termino fuente de la ecuacion de temperatura - ST [K/s]
volScalarField ST = LAD * qPlantSen / (rho * Cp0);

// Definir la ecuacion de temperatura con el termino fuente
alphat = turbulence->nut() / Prt;
alphat.correctBoundaryConditions();
volScalarField alphaEff("alphaEff", turbulence->nu() / Pr + alphat);

fvScalarMatrix TEqn
(
    fvm::div(phi, T)
  - fvm::laplacian(alphaEff, T)
 ==
    radiation->ST(rhoCpRef, const_cast<volScalarField&>(T))  // Correccion del tipo de dato
  + ST
);

TEqn.relax();
fvOptions.constrain(TEqn);
TEqn.solve();
radiation->correct();

// Mensajes de depuracion
Info << "Aplicando fvOptions.correct(T)" << endl;
fvOptions.correct(T);
Info << "Se ejecuto fvOptions.correct(T)" << endl;

// Actualizacion de rhok
rhok = 1.0 - beta * (T - TRef);

// Solve the momentum equation
MRF.correctBoundaryVelocity(U);

// Obtener referencias a los campos en `createFields.H`
// const volScalarField& LAD = mesh.lookupObject<volScalarField>("LAD");
// const volScalarField& rho = mesh.lookupObject<volScalarField>("rho");

// Definir el coeficiente de arrastre como una variable
dimensionedScalar Cd("C_d", dimless, 1.2);

// Definir el término fuente USource con signo negativo
fvMatrix<vector> USourceTerm(
    - fvm::Sp(rho * Cd * LAD * mag(U) / rho, U)  // Signo negativo para representar resistencia
);

tmp<fvVectorMatrix> tUEqn
(
    fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    USourceTerm
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

Info << "Aplicando restricciones de fvOptions antes de resolver UEqn" << endl;
fvOptions.constrain(UEqn);

if (simple.momentumPredictor())
{
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            (
              - ghf * fvc::snGrad(rhok)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    // Aplicar corrección de fvOptions después de resolver
    Info << "🔍 Aplicando fvOptions.correct(U)" << endl;
    fvOptions.correct(U);
    Info << "✅ Se ejecutó fvOptions.correct(U)" << endl;
}

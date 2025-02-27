// Solve the specific humidity equation

// Define el coeficiente de difusión efectivo para la ecuación de humedad
volScalarField DspecEff("DspecEff", turbulence->nu()/sigma_v + turbulence->nut()/sigma_vt);

// Solve the specific humidity equation
fvScalarMatrix specHumEqn
(                                                                                                                          
    fvm::div(phi, specHum)
 ==                                                                         
    fvm::laplacian(DspecEff, specHum)
	+ fvOptions(specHum)
);

specHumEqn.relax(); 
fvOptions.constrain(specHumEqn); 
specHumEqn.solve();
fvOptions.correct(specHum);

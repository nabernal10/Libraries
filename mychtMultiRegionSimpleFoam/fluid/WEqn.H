{   
    
	fvScalarMatrix specHumEqn                                                   
    (                                                                                                                          
        fvm::div(phi,specHum)                                                  
     ==                                                                         
        fvm::laplacian(turb.muEff(),specHum)                                          
      + fvOptions(rho, specHum)                                                 
    );                                                                         
    
	specHumEqn.relax(); 
	
    fvOptions.constrain(specHumEqn); 
	
    specHumEqn.solve();   
	
    fvOptions.correct(specHum);
}

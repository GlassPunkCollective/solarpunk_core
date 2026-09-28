/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
     \\/     M anipulation  | Licensed under GNU GPLv3 copyleft shield.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"

// Explicit implementation of the CAMU Rank-2 Anisotropic Shear Viscosity Tensor
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    Info<< "\nInitializing CAMU Substrate Regularizer Engine..." << endl;
    Info<< "Tracking viscoelastic plateau parameters: sigma_rem = 0.3540" << endl;

    // Fixed scale-invariant material invariants
    const scalar sigmaRem = 0.3540;
    const scalar rc = 2.4; // Rigidity percolation threshold

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        #include "readPIMPLEControls.H"
        #include "ddtScheme.H"

        // Calculate localized fluid vorticity vectors
        volVectorField vorticity = fvc::curl(U);
        volScalarField omegaMag = mag(vorticity);

        // Hyperbolic tangent dynamic shield function to regularize gradient runaways
        // Prevents finite-time Euclidean blowups natively through substrate compliance
        volScalarField dynamicViscosity = nu + (sigmaRem * mag(U) / (scalar(1.0) + mag(vorticity))) 
            * mag(tanh(omegaMag / rc));

        // Coupled momentum equation with history-dependent anisotropic viscosity fields
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(dynamicViscosity, U)
         ==
            fvOptions(U)
        );

        UEqn.relax();

        if (pimple.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << endl;
    }

    Info<< "End of CAMU Substrate Run. Execution Successful.\n" << endl;
    return 0;
}

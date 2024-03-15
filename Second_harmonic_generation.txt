#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fftw3.h" /* Fourier Transform Library */

/* Macros for real and imaginary parts */
#define REAL 0
#define IMAG 1

#define SIZE 15000 /* Number of mesh nodes */

int main()
{
    /* Vacuum impedance in ohms */
    double imp0 = 377.0;
    double e0 = 8.854 * pow(10.0, -12.0);
    double pi = M_PI;
    double mu0 = 4 * pi * pow(10.0, -7);
    double c = 1 / sqrt(e0 * mu0);

    /* Open files to save obtained data */
    FILE* videoFile;
    FILE* nodePlot1File;
    FILE* nodePlot2File;
    FILE* fftFile;
    FILE* fftInvFile;
    FILE* fftInvCheckFile;

    /* Open files */
    videoFile = fopen("Etvideo.txt", "w"); /* Write values of all nodes for various time instances */
    nodePlot1File = fopen("Enodos1.txt", "w"); /* " " for the time passing through the discontinuity */
    nodePlot2File = fopen("Enodos2.txt", "w"); /* " " for the time passed after the discontinuity */
    fftFile = fopen("datosTF.txt", "w"); /* Write FT for a point scanning all times */
    fftInvFile = fopen("datosTFI.txt", "w");
    fftInvCheckFile = fopen("datosTFIcheck.txt", "w");
    double ezt3[SIZE] = {0.}, hyt3[SIZE] = {0.}, dzt3[SIZE] = {0.}; /* Arrays for M field (A/m), E field (V/m), and D field (C/m^2) */

    int wall = 200; /* Discontinuity position */

    /* Parameters for the point source */
    double w = 2; /* Pulse frequency 10^15 rad/s */
    double tau = 30; /* Period in which the cosine^2 pulse is defined in fs (10^-15 s) */

    int numTime = 11000; /* Number of time divisions */
    double qTime = 0; /* Initial time */
    double deltatime = 0.004; /* Time step in fs */

    fprintf(archivoPlotNodos1, "%g 0\n%d 0\n%d 0\n", deltatime, SIZE, wall);
    fprintf(archivoPlotNodos2, "%g 0\n%d 0\n%d 0\n", deltatime, SIZE, wall);

    int j; /* j is the variable for the time loop */
    int mm; /* mm is the variable traversing the nodes of M and E fields */

    fprintf(archivoTF, "%d 0\n%f 0\n", numTime, deltatime); /* First row of the file is to define FT parameters in Mathematica */

    /* Variables for MATLAB video snapshots */
    int photoStep = 25;
    double photos = numTime / photoStep;
    fprintf(videoFile, "%d\n %g\n %g\n %d\n", SIZE, photos, deltatime, wall);

/*---------------------------------------------------------------------------------------------------*/
    /* Maxima and minima of transmitted and reflected waves */
    double maxTransf = 0, maxRefle = 0;
    double minTransf = 0, minRefle = 0;

    /* Positions of maxima and minima of transmitted and reflected waves */
    int pos1Trans, pos2Trans, pos1Refle, pos2Refle;

    /* Wavelength in vacuum and in dielectric */
    double lambda0, lambdaprima;

    /* Dispersion parameters */
    double f1 = 0.55, delta1 = 0.5 * f1;
    double w1 = 2 * pi * f1;

    double es = 3, einf = 2, G1 = 1;
    double b1 = G1 * (es - einf);

    double a1, c1, g1;
    double pzt1[SIZE] = {0.}, pzt2[SIZE] = {0.}, pzt3[SIZE] = {0.}, dzt1[SIZE] = {0.}, dzt2[SIZE] = {0.};

    /* Formulas: Artl Joseph, R. M., & Taflove, A. (1997) */
    a1 = 2 + 2 * delta1 * deltatime + pow(w1, 2) * pow(deltatime, 2) * (1 + b1);  /*17.a*/
    c1 = pow(w1, 2) * pow(deltatime, 2) * b1;                           /*17.b*/
    g1 = -2 + 2 * delta1 * deltatime - pow(w1, 2) * pow(deltatime, 2) * (1 + b1); /*17.c*/

/*----------------------------------------SECOND HARMONIC------------------------------------------*/
/* Weight of the second harmonic */
double factor = 8000000000; // Maximum value: 20000000000

double pzNLt3[SIZE] = {0.}; /* Array for second-order nonlinearity */

/*--------------------------------------FOURIER TRANSFORM-----------------------------------------*/
double absValue;
int fftMeasure = wall + 50;

/* Input array for FFT */
fftw_complex x[numTime]; /* equivalent to double x[n][2] */

/* Output array for FFT */
fftw_complex y[numTime];
/*--------------------------------------------------------------------------------------------------*/
    /* Initialize fields to zero for all nodes */
    for (mm = 0; mm < SIZE; mm++){
       ezt3[mm] = 0.0;
       dzt3[mm] = 0.0;
       pzt1[mm] = 0.0;
       pzt2[mm] = 0.0;
       pzt3[mm] = 0.0;
       dzt2[mm] = 0.0;
       dzt1[mm] = 0.0;
       pzNLt3[mm] = 0.0;
    }
    for (mm = 0; mm < SIZE - 1; mm++){
       hyt3[mm] = 0.0;
    }

    /* Time loop */
    for (j = 0; j < numTime ; j++){

        /* Leap-Frog method */
        /* M field loop in space */
        for (mm = 0; mm < SIZE - 1; mm++){
            hyt3[mm] = hyt3[mm] + (ezt3[mm + 1] - ezt3[mm]) / imp0;
        }

        /* Source at node 0 */
        if(qTime <= tau){
            ezt3[0] = pow((sin((pi * qTime) / tau)), 2) * cos(w * qTime);
        }
        else{
            ezt3[0] = ezt3[1]; /* Boundary condition at the vacuum wall


		/* Loop for E field in space */
		for (mm = 1; mm < wall; mm++){ /* Vacuum propagation */
			dzt3[mm] = dzt3[mm] + (hyt3[mm] - hyt3[mm - 1]) * (imp0) * e0;
			ezt3[mm] = dzt3[mm] / e0;
		}

		for (mm = wall; mm < SIZE; mm++){ /* Dielectric propagation without losses */
			pzt1[mm] = pzt2[mm];
			pzt2[mm] = pzt3[mm];

			dzt1[mm] = dzt2[mm];
			dzt2[mm] = dzt3[mm];

			dzt3[mm] = dzt3[mm] + (hyt3[mm] - hyt3[mm - 1]) * (imp0) * e0;

			/* Formulas: Artl Joseph, R. M., & Taflove, A. (1997) */
			pzt3[mm] = (4 * pzt2[mm] + g1 * pzt1[mm] + c1 * (dzt3[mm] + dzt1[mm])) / a1; /* 16.a */

			pzNLt3[mm] = factor * ((pzt3[mm] + pzt2[mm]) / 2) * ((pzt3[mm] + pzt2[mm]) / 2);

			ezt3[mm] = (dzt3[mm] - (pzt3[mm] + pzNLt3[mm])) / (e0 * einf);
		}

		x[j][REAL] = ezt3[fftMeasure]; /* Point where the FT is performed. After traversing much material */
		x[j][IMAG] = 0;
		fprintf(fftInvCheckFile, "%d %f\n", j, ezt3[fftMeasure]);

		for (mm = 0; mm < SIZE; mm++){

			/* Save all ez values for various times: when j is divisible by 20 */
			if (j % photoStep == 0){
				fprintf(videoFile, "%g\n", ezt3[mm]);
			}
			double photo1 = 21;
			/* Save all ez values for a specific time */
			if ((qTime > photo1 - deltatime / 2) & (qTime < photo1 + deltatime / 2)){
				fprintf(nodePlot1File, "%d %g \n", mm, ezt3[mm]);
			}
			double photo2 = 40;
			if ((qTime > photo2 - deltatime / 2) & (qTime < photo2 + deltatime / 2)){
				fprintf(nodePlot2File, "%d %g \n", mm, ezt3[mm]);

				/* Conditions to obtain the value and position of the maximum and minimum of the transmitted and reflected waves */
				if ((mm > wall) & (maxTransf < ezt3[mm])){
					maxTransf = ezt3[mm];
					pos1Trans = mm;
				}
				if ((mm > wall) & (minTransf > ezt3[mm])){
					minTransf = ezt3[mm];
					pos2Trans = mm;
				}
				if ((mm < wall) & (maxRefle < ezt3[mm])){
					maxRefle = ezt3[mm];
					pos1Refle = mm;
				}
				if ((mm < wall) & (minRefle > ezt3[mm])){
					minRefle = ezt3[mm];
					pos2Refle = mm;
				}
			}
		}
		/* Time step */
		qTime = qTime + deltatime;
		//printf ("%g %g %g %g\n",ezt3[fftMeasure],dzt3[fftMeasure], pzt3[fftMeasure], pzNLt3[fftMeasure]);
	}

/* Execute Fourier Transform (FT) */
fftw_plan p;
p = fftw_plan_dft_1d(numTime, x, y, FFTW_FORWARD, FFTW_ESTIMATE); // Creating a plan for forward FT
fftw_execute(p); // Executing the FT

/* Output the FT results */
for (j = 0; j < numTime; j++) {
    valorabs = sqrt(pow(y[j][REAL], 2) + pow(y[j][IMAG], 2)); // Calculate the absolute value
    fprintf(archivoTF, "%d %f\n", j, valorabs); // Write the result to a file
}

/* Clean up */
fftw_destroy_plan(p); // Destroy the plan

/* Execute inverse Fourier Transform (IFT) */
fftw_plan p2;
p2 = fftw_plan_dft_1d(numTime, y, x, FFTW_BACKWARD, FFTW_ESTIMATE); // Creating a plan for backward IFT
fftw_execute(p2); // Executing the IFT

/* Output the IFT results */
for (j = 0; j < numTime; j++) {
    fprintf(archivoTFI, "%d %f \n", j, x[j][REAL] / numTime); // Write the result to a file
}

/* Clean up */
fftw_destroy_plan(p2); // Destroy the plan

/* Print simulation completion message */
printf("Simulation completed for time: %g\n", qTime);

/* Print maximum and minimum of transmitted and reflected waves */
printf("Max and min of the transmitted: %g and %g\n", maxtransf, mintransf);
printf("Max and min of the reflected: %g and %g\n", maxrefle, minrefle);

/* The wavelength in vacuum will be twice the absolute distance between the maximum and minimum of the reflected in nm */
lambda0 = abs((pos1refle - pos2refle) * 2);
lambdaprima = abs((pos1trans - pos2trans) * 2);
printf("Wavelength in vacuum: %g\nand in the other medium: %g\n", lambda0, lambdaprima);

/* Print dispersion coefficients */
printf("Dispersion coefficients: %g %g %g", a1, c1, g1);

/* Close files when the temporal loop ends */
fclose(archivoVideo);
fclose(archivoPlotNodos1);
fclose(archivoPlotNodos2);
fclose(archivoTF);
fclose(archivoTFI);
fclose(archivoTFIcheck);

return 0;
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Number of nodes in the mesh */
#define SIZE 14000

int main() {
    /* Impedance of free space in ohms */
    double imp0 = 377.0;
    double e0 = 8.854 * pow(10.0, -12.0);
    double pi = M_PI;
    double mu0 = 4 * pi * pow(10.0, -7);
    double c = 1 / sqrt(e0 * mu0);

    /* Open files to store obtained data */
    FILE *archivoPlotNodos1;
    FILE *archivoPlotNodos2;
    FILE *archivoTF;
    FILE *archivoIntensidadSHG;

    /* Open files */
    archivoPlotNodos1 = fopen("Enodos1.txt", "w"); /* For the time it takes to cross the discontinuity */
    archivoPlotNodos2 = fopen("Enodos2.txt", "w"); /* For the time that has passed the discontinuity */
    archivoTF = fopen("datosTF.txt", "w"); /* Write the TF for a point sweeping all times */
    archivoIntensidadSHG = fopen("datosIntensidadSHG.txt", "w");

    /* Declare arrays for the M (A/m), E (V/m), and D (C/m^2) fields */
    double ezt3[SIZE] = {0.}, hyt3[SIZE] = {0.}, dzt3[SIZE] = {0.};
    double eyt3[SIZE] = {0.}, hzt3[SIZE] = {0.}, dyt3[SIZE] = {0.};

    int wall = 100; /* Position of the discontinuity */

    /* Parameters of the point source */
    double wz = 1; /* Pulse frequency 10^15 rad/s */
    double tau = 40; /* Period in which the cosine^2 pulse is defined in fs (10^-15 s) */

    int numTime = 25000; /* Number of time divisions */
    double qTime = 0; /* Initial time */
    double deltatime = 0.01; /* Time step in fs */

    double foto1, foto2;
    fprintf(archivoPlotNodos1, "%g 0 0\n%d 0 0\n%d 0 0\n", deltatime, SIZE, wall);
    fprintf(archivoPlotNodos2, "%g 0 0\n%d 0 0\n%d 0 0\n", deltatime, SIZE, wall);

    int j; /* j is the variable of the time loop */
    int mm; /* mm is the variable that goes through the nodes of the M and E fields */

    /* The first row of the file is to define the parameters of the TF in Mathematica */
    fprintf(archivoTF, "%d %f\n%f 0\n", numTime, deltatime, wz);

        /*----------TRANSMITTED AND REFLECTED--------*/
    /* Maximums and minimums of the transmitted and reflected */
    double maxtransfEz = 0, maxrefleEz = 0;
    double mintransfEz = 0, minrefleEz = 0;

    /* Positions of maximums and minimums of the transmitted and reflected */
    int pos1transEz, pos2transEz, pos1refleEz, pos2refleEz;

    /* Wavelength in vacuum and in the dielectric */
    double lambda0Ez, lambdaprimaEz;

    /* Parameters of dispersion: Electric field in z */
    double f1z = 0.7, delta1z = 0.5 * f1z;
    double w1z = 2 * pi * f1z;

    double esz = 2, einfz = 1, G1z = 1;
    double b1z = G1z * (esz - einfz);

    double a1z, c1z, g1z;
    double pzt1[SIZE] = {0.}, pzt2[SIZE] = {0.}, pzt3[SIZE] = {0.}, dzt1[SIZE] = {0.}, dzt2[SIZE] = {0.};

    /* Formulas: Artl Joseph, R. M., & Taflove, A. (1997) */
    a1z = 2 + 2 * delta1z * deltatime + pow(w1z, 2) * pow(deltatime, 2) * (1 + b1z);  /* 17.a */
    c1z = pow(w1z, 2) * pow(deltatime, 2) * b1z;                                      /* 17.b */
    g1z = -2 + 2 * delta1z * deltatime - pow(w1z, 2) * pow(deltatime, 2) * (1 + b1z);  /* 17.c */

    /* Parameters of dispersion: Electric field in y */
    double f1y = 0.5, delta1y = 0.2 * f1y;
    double w1y = 2 * pi * f1y;

    double esy = 1.6289, einfy = 1, G1y = 1; // esy = 1.62702
    double b1y = G1y * (esy - einfy);

    double a1y, c1y, g1y;
    double pyt1[SIZE] = {0.}, pyt2[SIZE] = {0.}, pyt3[SIZE] = {0.}, dyt1[SIZE] = {0.}, dyt2[SIZE] = {0.};

    /* Formulas: Artl Joseph, R. M., & Taflove, A. (1997) */
    a1y = 2 + 2 * delta1y * deltatime + pow(w1y, 2) * pow(deltatime, 2) * (1 + b1y);  /* 17.a */
    c1y = pow(w1y, 2) * pow(deltatime, 2) * b1y;                                      /* 17.b */
    g1y = -2 + 2 * delta1y * deltatime - pow(w1y, 2) * pow(deltatime, 2) * (1 + b1y);  /* 17.c */

    /*-------------SECOND HARMONIC------------*/
    /* Weight of the second harmonic */
    double factorez = 8000000;
    double factorey = 10000000000;

    /* Declare arrays for the second order nonlinearity */
    double pzNLt3[SIZE] = {0.};
    double pyNLt3[SIZE] = {0.};

    /*------FOURIER TRANSFORM--------*/

    int medidaTF = wall + 10;
    int mTF = 0;
    int sumamTF = 20;
    int mTFlimit = 130;
    fprintf(archivoIntensidadSHG, "%d\n%f\n%f\n%d\n%d\n%d\n", numTime, deltatime, wz, mTFlimit, sumamTF, medidaTF);

    /* FDTD */
    // Activate this loop if you want to perform the FT
    // for mTFlimit points of the mesh with a step of sumamTF
    // for (mTF = 0; mTF < mTFlimit; mTF++) {
    /* Initialize fields to zero for all nodes */
    for (mm = 0; mm < SIZE; mm++) {
        ezt3[mm] = 0.0;
        dzt3[mm] = 0.0;
        pzt1[mm] = 0.0;
        pzt2[mm] = 0.0;
        pzt3[mm] = 0.0;
        dzt2[mm] = 0.0;
        dzt1[mm] = 0.0;
        pzNLt3[mm] = 0.0;
        // Perpendicular component
        eyt3[mm] = 0.0;
        dyt3[mm] = 0.0;
        pyt1[mm] = 0.0;
        pyt2[mm] = 0.0;
        pyt3[mm] = 0.0;
        dyt2[mm] = 0.0;
        dyt1[mm] = 0.0;
        pyNLt3[mm] = 0.0;
    }
    for (mm = 0; mm < SIZE - 1; mm++) {
        hyt3[mm] = 0.0;
        // Perpendicular component
        hzt3[mm] = 0.0;
    }

    /* Temporal loop */
    for (j = 0; j < numTime; j++) {
        /* Leap-Frog method */
        /* Loop for field M in space */
        for (mm = 0; mm < SIZE - 1; mm++) {
            hyt3[mm] = hyt3[mm] + (ezt3[mm + 1] - ezt3[mm]) / imp0;
            hzt3[mm] = hzt3[mm] + (eyt3[mm + 1] - eyt3[mm]) / imp0;
        }

 
	        /* Source at node 0 */
        if (qTime <= tau) {
            ezt3[0] = pow((sin((pi * qTime) / tau)), 2) * cos(wz * qTime);
            // Not programming the source in y
            // eyt3[0]=pow((sin((pi*qTime)/tau)),2)*cos(wy*qTime);
        }
        else {
            ezt3[0] = ezt3[1]; /* Boundary condition at the vacuum wall */
            eyt3[0] = eyt3[1];
        }

        /* Loop for field E in space */
        for (mm = 1; mm < wall; mm++) { /* Vacuum propagation */
            dzt3[mm] = dzt3[mm] + (hyt3[mm] - hyt3[mm - 1]) * (imp0) * e0;
            ezt3[mm] = dzt3[mm] / e0;

            dyt3[mm] = dyt3[mm] + (hzt3[mm] - hzt3[mm - 1]) * (imp0) * e0;
            eyt3[mm] = dyt3[mm] / e0;
        }

        /* Dielectric propagation without losses */
        for (mm = wall; mm < SIZE; mm++) {
            pzt1[mm] = pzt2[mm];
            pzt2[mm] = pzt3[mm];

            dzt1[mm] = dzt2[mm];
            dzt2[mm] = dzt3[mm];
            dzt3[mm] = dzt3[mm] + (hyt3[mm] - hyt3[mm - 1]) * (imp0) * e0;

            /* Formulas: Artl Joseph, R. M., & Taflove, A. (1997) */
            pzt3[mm] = (4 * pzt2[mm] + g1z * pzt1[mm] + c1z * (dzt3[mm] + dzt1[mm])) / a1z; /* Equation 16.a */
            pzNLt3[mm] = factorez * ((pzt3[mm] + pyt2[mm]) / 2) * ((pzt3[mm] + pyt2[mm]) / 2);
            ezt3[mm] = (dzt3[mm] - (pzt3[mm] + pzNLt3[mm])) / (e0 * einfz);

            // Perpendicular component
            pyt1[mm] = pyt2[mm];
            pyt2[mm] = pyt3[mm];

            dyt1[mm] = dyt2[mm];
            dyt2[mm] = dyt3[mm];
            dyt3[mm] = dyt3[mm] + (hzt3[mm] - hzt3[mm - 1]) * (imp0) * e0;

            pyt3[mm] = (4 * pyt2[mm] + g1y * pyt1[mm] + c1y * (dyt3[mm] + dyt1[mm])) / a1y;
            pyNLt3[mm] = factorey * ((pyt3[mm] + pzt2[mm]) / 2) * ((pyt3[mm] + pzt2[mm]) / 2);
            eyt3[mm] = (dyt3[mm] - (pyt3[mm] + pyNLt3[mm])) / (e0 * einfy);

            //printf ("%d %g %g\n",mm,ezt3[mm],eyt3[mm]);
        }

        //fprintf (archivoTF,"%f %f\n",eyt3[medidaTF],ezt3[medidaTF]);
        fprintf(archivoIntensidadSHG, "%f\n", eyt3[medidaTF]);

        for (mm = 0; mm < SIZE; mm++) {
            /* Save all values of ez for a certain time */
            foto1 = 28;
            if ((qTime > foto1 - deltatime / 2) & (qTime < foto1 + deltatime / 2)) {
                fprintf(archivoPlotNodos1, "%d %g %g\n", mm, eyt3[mm], ezt3[mm]);
            }
            foto2 = 170;
            if ((qTime > foto2 - deltatime / 2) & (qTime < foto2 + deltatime / 2)) {
                fprintf(archivoPlotNodos2, "%d %g %g\n", mm, eyt3[mm], ezt3[mm]);

                /* Conditions to obtain the value and position of the maximum and minimum transmitted and reflected */
                if ((mm > wall) & (maxtransfEz < ezt3[mm])) {
                    maxtransfEz = ezt3[mm];
                    pos1transEz = mm;
                }
                if ((mm > wall) & (mintransfEz > ezt3[mm])) {
                    mintransfEz = ezt3[mm];
                    pos2transEz = mm;
                }
                if ((mm < wall) & (maxrefleEz < ezt3[mm])) {
                    maxrefleEz = ezt3[mm];
                    pos1refleEz = mm;
                }
                if ((mm < wall) & (minrefleEz > ezt3[mm])) {
                    minrefleEz = ezt3[mm];
                    pos2refleEz = mm;
                }
            }
        }
        /* Time step */
        qTime = qTime + deltatime;
    }
//medidaTF=medidaTF+sumamTF;
//qTime=0;
//}

    printf("Simulation completed for time: %g\n", qTime);
    printf("Max and min of the transmitted: %g and %g\n", maxtransfEz, mintransfEz);
    printf("Max and min of the reflected: %g and %g\n", maxrefleEz, minrefleEz);

    /* The wavelength in vacuum will be twice the distance between the absolute minimum and maximum of the reflected in nm */
    lambda0Ez = abs((pos1refleEz - pos2refleEz) * 2);
    lambdaprimaEz = abs((pos1transEz - pos2transEz) * 2);
    printf("Wavelength in vacuum: %g\nand in the other medium: %g\n", lambda0Ez, lambdaprimaEz);
    printf("Dispersion coefficients in z: %g %g %g \n\n", a1z, c1z, g1z);
    printf("Dispersion coefficients in y: %g %g %g \n\n", a1y, c1y, g1y);

    /* Close files when the time loop ends */
    fclose(archivoPlotNodos1);
    fclose(archivoPlotNodos2);
    fclose(archivoTF);

    return 0;
}


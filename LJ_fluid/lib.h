#ifndef __LIB_H__
#define __LIB_H__

#define RHO 0.860   // densita' - unita ridotte
#define T 0.849     // temperatura - unita ridotte
#define N 256       // particelle per riempire il quadratino


void generate_FCC(double L, double X[][3]);
void tuttoZero(double v[][3]);
void oneStepVerlet(double X[][3], double V[][3], double a[][3], double L, double dt);
void calculateForce(double X[][3], double a[][3], double L);
void rescaleVelocities(double V[][3]);
double U_pot(double X[][3], double L);
double E_cinetica(double V[][3]);
double mean(double X[], int size);
void F_distribuzione(const double X[][3], double g[], double dr, double L);


#endif

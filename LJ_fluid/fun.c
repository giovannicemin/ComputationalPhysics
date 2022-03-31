/*
 * file dove inserire tutte le cose in piu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"

void generate_FCC(double L, double X[][3])
{
  /*
    This function fills the matrix X[N][3] with the coordinates of
    particles in a FCC lattice within a box of side L whose center is the
    center of the coordinate system.
    The box in filled uniformly if N = 4 n^3 with n integer hence for
    N = 4     32    108    256    500    864   1372   2048   2916
    4000   5324 6912   8788  10976  13500  16384 ...
  */

  int i,j,k,m,n,p,c;

  /* position within the primary cell */
  const double rFCC[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
                             {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
  double b, rCell[3];

  for (c = 1; ; c++)
    if (4*c*c*c >= N){
      break;}

  b = L / (double)c;            /* side length of the primary cell */
  p = 0;                        /* particles placed so far */
  for (i = 0; i < c; i++)
    {
      rCell[0] = i;
      for (j = 0; j < c; j++)
        {
          rCell[1] = j;
          for (k = 0; k < c; k++)
            {
              rCell[2] = k;
              for (m = 0; m < 4; m++) /* 4 particles in cell */
                if (p < N)
                  {

                    /* add the com to each bead, and project to the real cell */
                    for(n=0;n<3;n++)
                      {
                        X[p][n] = b * (rCell[n] + rFCC[m][n]);
                        X[p][n] -= L * rint(X[p][n]/L);
                      }
                ++p;
                  }
            }
        }
    }
}


void tuttoZero(double v[][3])
{
    for(int i=0; i<N; i++){
        v[i][0] = v[i][1] = v[i][2] = 0;
    }
}

void oneStepVerlet(double X[][3], double V[][3], double a[][3], double L, double dt)
{
    double v_temp[N][3], dx[3], r;

    for(int i=0; i<N; i++){ // cicla attraverso le particelle che si muovono
        for(int k=0; k<3; k++){
            v_temp[i][k] = V[i][k] + dt*0.5*a[i][k];
            X[i][k] = X[i][k] + dt*v_temp[i][k];

            // effetto pacman
            X[i][k] -= L*rint(X[i][k]/L);
            //printf("%lf ", X[i][k]);
        }
        //printf("\n");
    }

    calculateForce(X, a, L);

    for(int i=0; i<N; i++){
        for(int k=0; k<3; k++){
            V[i][k] = v_temp[i][k] + dt*0.5*a[i][k];
        }
    }
}

void calculateForce(double X[][3], double a[][3], double L)
{
    double dx[3], r_sq;
    tuttoZero(a);

    // calcolo le forze sapendo che le particelle interagiscono mediante LJ
    for(int i=0; i<N; i++){ // cicla attraverso le particelle su cui calcolo la forza

        for(int j=0; j<N; j++) // cicla attraverso le particelle che interagiscono con la i-esima
          if(j != i){

            for(int k=0; k<3; k++){ // calcolo la distanza fra i-esima e j-esima
                dx[k] = X[i][k] - X[j][k];
                dx[k] -= L*rint(dx[k]/L);
            }
            r_sq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
            //printf("%d %d %lf \n", i, j, r_sq);

            for(int k=0; k<3; k++) // calcolo l'accelerazione DI UNA SOLA PARTICELLA
              if(r_sq < 9){        // CUT-OFF
                a[i][k] += 24.0*( 2.0*pow(r_sq, -7) - pow(r_sq, -4) )*dx[k];
                //printf(" %lf ", a_temp[k]);
            }
            //printf("\n");
        }
    }
}

void rescaleVelocities(double V[][3])
{
    double k=0, alpha, temperatura;

    // come prima cosa calcolo la temperatura attuale e quindi alpha
    k = E_cinetica(V);

    temperatura = (2.0*k)/(3.0*N);
    alpha = sqrt(T/temperatura);
    //printf("%lf %lf \n", temperatura, alpha);

    // a questo punto riscalo tutte le velocita'
    for(int i=0; i<N; i++){
        for(int m=0; m<3; m++){
            V[i][m] *= alpha;
        }
    }
}

double U_pot(double X[][3], double L)
{
    double dx[3], tot = 0, r_sq;

    for(int i=0; i<N-1; i++){
        for(int j=i+1; j<N; j++){
            for(int k=0; k<3; k++){ // calcolo la distanza fra i-esima e j-esima
                dx[k] = X[i][k] - X[j][k];
                dx[k] -= L*rint(dx[k]/L);
            }
            r_sq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

            if(r_sq<9){ tot += 4*( pow(r_sq, -6) - pow(r_sq, -3) );}
        }
    }
    return tot;
}

double E_cinetica(double V[][3])
{
    double tot = 0;

    for(int i=0; i<N; i++){
        tot += pow(V[i][0],2) + pow(V[i][1],2) + pow(V[i][2],2);
    }
    tot *= 0.5;
    return tot;
}

void F_distribuzione(const double X[][3], double g[], double dr, double L)
{
    double r, dx[3];
    int pos;

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++)
          if(j!=i){
            for(int k=0; k<3; k++){ // calcolo la distanza fra i-esima e j-esima
                dx[k] = X[i][k] - X[j][k];
                dx[k] -= L*rint(dx[k]/L);
            }
            r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

            if(r<(L/2)){g[(int) floor(r/dr)]++; }

      }
  }
  //printf("\n");
}

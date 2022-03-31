/* Cemin Giovanni 11.03.2021
**
** Programmino per i liquido di Lennard-Jones
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"

int main(int argc, char *argv[])
{
    double X[N][3];                  // posizioni delle particelle
    double V[N][3];                  // velocita' delle particelle
    double a[N][3];                  // accelerazioni delle particelle
    double L = pow(N/RHO, 1.0/3);    // lunghezza della scatola

    double dt = 0.005;               // incremento del tempo per Verlet

    //FILE * file = fopen("potenziale.txt", "w");
    FILE * file2 = fopen("struttura.txt", "w");

    // come prima cosa creo il mio cubetto si liquido
    generate_FCC(L, X);
    tuttoZero(V);
    calculateForce(X, a, L);

    // adesso devo far partire il tutto

    // ==================== FASE DI EQUILIBRAZIONE =======================
    // Quello che faccio e' equilibrare il mio sistema ogni 50 passaggi di Verlet
    // vediamo se questa soluzione va bene
    double tempo = 0, max;
    int cont = 0;

    while(tempo < 10){
        printf("Tempo: %lf \n", tempo);

        oneStepVerlet(X, V, a, L, dt);

        max = 0;
        for(int i=0; i<N; i++){
            //printf("%d \n", i); //printf("%lf %lf %lf \n", X[i][0], X[i][1], X[i][2]); //printf("%lf %lf %lf \n", a[i][0], a[i][1], a[i][2]);
        }

        if(cont == 50){
            // qui riscalo le mie velocita in modo da metchare la temperatura
            rescaleVelocities(V);
            cont = 0;
        }

        // vedere come si comporta il potenziale
        //fprintf(file, "%lf %lf \n", tempo, U_pot(X, L));

        tempo += dt;
        cont ++;
    }

    // ====================== FASE DI PRODUZIONE ==============================
    // in questa sezione, dopo avere un sistema per bene si prendono i dati che
    // ci servono

    double vmPotenziale = 0, vmCinetica = 0;

    // variabili per proprieta' strutturali
    int divisioni = 1000;
    double dr = (L/2)/divisioni, g[divisioni + 1];

    for(int i=0; i<divisioni; i++){
        g[i] = 0;
    }

    cont = 0;
    while(tempo < 90){
        printf("Tempo: %lf \n", tempo);

        if(fabs(rint(tempo)-tempo) < 0.0001){
            // valor medio potenziale
            vmPotenziale += U_pot(X, L);

            // valor medio energia cinetica
            vmCinetica += E_cinetica(V);

            // proprieta' strutturali
            F_distribuzione(X, g, dr, L);

            cont ++;
        }

        oneStepVerlet(X, V, a, L, dt);
        //fprintf(file, "%lf %lf \n", tempo, U_pot(X, L));

        tempo += dt;
    }

    // rinormalizzazione
    vmPotenziale /= cont;
    vmCinetica /= cont;

    // correzione all'energia potenziale
    vmPotenziale -= N*8.0*M_PI*RHO/(81.0); // 81 = 3*c^3 ove c e' il cutoff =3

    printf("Valor medio potenziale = %lf \nValor medio per particella = %lf \n",
           vmPotenziale, vmPotenziale/N);

    printf("Valor medio energia cinetica = %lf \nTemperatura media = %lf \n",
           vmCinetica, vmCinetica*(2.0/3.0)/N);

    for(int i=0; i<divisioni; i++){
        g[i] /= (cont * 4*M_PI*N*RHO * pow(dr*i + dr/2.0, 2)*dr);
        fprintf(file2, "%lf %lf \n", dr*i + dr/2.0, g[i]);
    }

    return 0;
}


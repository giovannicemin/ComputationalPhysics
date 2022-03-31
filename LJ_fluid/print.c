/*
** Programmino che serve per fare i grafici con GNU Plot
** molto bene, cominciamo
*/

#include <stdio.h>
#include <stdlib.h>
#define N_POINTS 1000

int main(int argc, char *argv[])
{
    // faccio tutti e 5 i grafici insieme
    char * commandsForGnuplot[] = {"set title \"grafico E = 5\" ",
                                   "set xrange [0:5]",
                                   "set yrange [-7:4]",
                                   "set grid",
                                   "plot 'energia.txt' lc 7 w lines"};
                                   //"replot 'energia_0-1_Integral.txt' lc 2 w lines"};


    int numCommands = 5;
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

    for (int i=0; i < numCommands; i++)
    {
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
    }

}

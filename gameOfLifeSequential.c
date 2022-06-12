#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*
#define N_ITER 10
#define gridRows 14
#define gridColumns 5
*/

#define ALIVE '1'
#define DEAD '0'

double what_time_is_it()
{
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    return now.tv_sec + now.tv_nsec*1e-9;
}

int main(int argc, char *argv[])
{
    double t1 = what_time_is_it();

    int gridRows = atoi(argv[1]);
    int gridColumns = atoi(argv[2]);
    int n_iter = atoi(argv[3]);

    /*
    //stack allocation
    char grid0[gridRows][gridColumns] = {{ALIVE, DEAD, ALIVE, DEAD, DEAD},
                                        {ALIVE, DEAD, ALIVE, ALIVE, DEAD},
                                        {ALIVE, ALIVE, ALIVE, ALIVE, DEAD},
                                        {ALIVE, DEAD, ALIVE, ALIVE, ALIVE},
                                        {ALIVE, DEAD, ALIVE, DEAD, DEAD},
                                        {ALIVE, DEAD, ALIVE, DEAD, DEAD},
                                        {ALIVE, DEAD, ALIVE, DEAD, DEAD},
                                        {ALIVE, DEAD, ALIVE, DEAD, DEAD},
                                        {ALIVE, ALIVE, ALIVE, DEAD, ALIVE},
                                        {ALIVE, DEAD, ALIVE, DEAD, DEAD},
                                        {ALIVE, DEAD, ALIVE, DEAD, DEAD},
                                        {ALIVE, DEAD, ALIVE, DEAD, DEAD},
                                        {ALIVE, DEAD, ALIVE, DEAD, DEAD},
                                        {ALIVE, DEAD, ALIVE, DEAD, DEAD}}; //uncomment to initialize
                                        
    char grid1[gridRows][gridColumns];
    char (*grid)[gridRows][gridColumns] = &grid0;
    char (*newGrid)[gridRows][gridColumns] = &grid1;
    */

    //heap + random allocation
    char (*grid)[gridRows][gridColumns] = malloc((size_t)gridRows*(size_t)gridColumns);
    char (*newGrid)[gridRows][gridColumns] = malloc((size_t)gridRows*(size_t)gridColumns);
    
    //inizializzazione casuale
    time_t t;
    srand((unsigned) time(&t));
    for (int row=0; row<gridRows; row++) {
        for (int col=0; col<gridColumns; col++) {
            if (rand() % 2 == 0) {
                (*grid)[row][col] = DEAD;
            } else {
                (*grid)[row][col] = ALIVE;
            }
        }
    }

    /* //stampa prima griglia
    for (int row=0; row<gridRows; row++) {
        for (int col=0; col<gridColumns; col++) {
            printf("%c", (*grid)[row][col]);
        }
        printf("\n");
    }
    */

    for (int iter = 0; iter < n_iter; iter++) {
        for (int row=0; row<gridRows; row++) { //per ogni cella
            for (int col=0; col<gridColumns; col++) {
                int aliveNeighbours = 0; //conta i vicini ALIVE
                for (int row1 = row - 1; row1 <= row + 1; row1++) { //righe sopra-corrente-sotto
                    for (int col1 = col - 1; col1 <= col + 1; col1++) { //colonne sinistra-corrente-destra
                        if (!(row1 == row && col1 == col)) { //non contare la cella stessa come aliveNeighbour
                            if (row1 >= 0 &&  col1 >= 0) { //check bound inferiori griglia
                                if (row1 < gridRows && col1 < gridColumns) { //check bound superiori griglia
                                    if ((*grid)[row1][col1] == ALIVE) {
                                        aliveNeighbours++;
                                    }
                                }
                            }
                        }
                    }
                }

                //calcola il nuovo stato delle celle
                (*newGrid)[row][col] = (*grid)[row][col];
                if ((*grid)[row][col] == ALIVE) {
                    if (aliveNeighbours < 2 || aliveNeighbours > 3) {
                        (*newGrid)[row][col] = DEAD;
                    }
                } else {
                    if (aliveNeighbours == 3) {
                        (*newGrid)[row][col] = ALIVE;
                    }
                }
            }
        }

        /*
        //stampa griglia aggiornata
        printf("\n");
        for (int row=0; row<gridRows; row++) {
            for (int col=0; col<gridColumns; col++) {
                printf("%c", (*newGrid)[row][col]);
            }
            printf("\n");
        }
        */

        //swappa le griglie
        char (*temp)[gridRows][gridColumns] = grid;
        grid = newGrid;
        newGrid = temp;
    }

    
    printf("time taken %.6lfs\n", what_time_is_it() - t1);
    return 0;
}

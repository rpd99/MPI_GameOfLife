#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#define MASTER 0
#define FIRST_WORKER_RANK 1

#define ALIVE '1'
#define DEAD '0'

struct workPortion {
    int firstRow, rows;
    int topGhost, bottomGhost; //presenza o meno di ghost row, possibile cambio nome in isXGhost
};

int totalRows, totalCols, nIters;

void initializeGrid (int gridRows, int gridCols, char (*grid) [][gridCols])
{
    time_t t;
    srand((unsigned) time(&t));
    for (int row = 0; row < gridRows; row++) {
        for (int col = 0; col < gridCols; col++) {
            if (rand() % 2 == 0) {
                (*grid)[row][col] = DEAD;
            } else {
                (*grid)[row][col] = ALIVE;
            }
        }
    }
}

/* Data una griglia M*N determina la miglior partizione della griglia per nWorkers worker */
void defineWorkPortions (int gridRows, int gridCols, struct workPortion portions[], int nPortions)
{
    int rowsPerPortion = gridRows / nPortions;
    int remainderRows = gridRows % nPortions;

    for (int i = 0; i < nPortions; i++) {
        portions[i].topGhost = portions[i].bottomGhost = true;
        portions[i].firstRow = rowsPerPortion * i; //determina prima riga della porzione i-esima
        portions[i].rows = rowsPerPortion; //determina quante righe per la porzione i-esima
    }

    portions[0].topGhost = portions[nPortions - 1].bottomGhost = false; //prima e ultima porzione hanno solo una ghost row
    portions[nPortions - 1].rows += remainderRows; //assegna le righe restanti all'ultima porzione
}

void subGameOfLife (int gridRows, int gridCols, char (*currGrid) [][gridCols], char (*nextGrid) [][gridCols], int firstGameRow, int gameRows)
{
    //per ogni cella
    for (int row = firstGameRow; row < firstGameRow + gameRows; row++) {
        for (int col = 0; col < gridCols; col++) {
            int aliveNeighbours = 0; //conta i vicini ALIVE
            for (int row1 = row - 1; row1 <= row + 1; row1++) { //righe sopra-corrente-sotto
                for (int col1 = col - 1; col1 <= col + 1; col1++) { //colonne sinistra-corrente-destra
                    if (!(row1 == row && col1 == col)) { //non contare la cella stessa come aliveNeighbour
                        if (row1 >= 0 &&  col1 >= 0) { //check bound inferiore/sinistro griglia
                            if (row1 < gridRows && col1 < gridCols) { //check bound superiore/destro griglia
                                if ((*currGrid)[row1][col1] == ALIVE) {
                                    aliveNeighbours++;
                                }
                            }
                        }
                    }
                }
            }

            //calcola il nuovo stato delle celle
            (*nextGrid)[row][col] = (*currGrid)[row][col];
            if ((*currGrid)[row][col] == ALIVE) {
                if (aliveNeighbours < 2 || aliveNeighbours > 3) {
                    (*nextGrid)[row][col] = DEAD;
                }
            } else {
                if (aliveNeighbours == 3) {
                    (*nextGrid)[row][col] = ALIVE;
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);
    double t1 = MPI_Wtime();

    int myRank, nProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    //usage?
    totalRows = atoi(argv[1]);
    totalCols = atoi(argv[2]);
    nIters = atoi(argv[3]);

    int nWorkers = nProcs - FIRST_WORKER_RANK;
    struct workPortion p[nWorkers];
    defineWorkPortions(totalRows, totalCols, p, nWorkers);
    
    //printf("args: nRows %d nCols %d nIters %d\n", totalRows, totalCols, nIters);

    if (myRank == MASTER) {
        /*//Stack allocation
        char grid0[14][5] = {{ALIVE, DEAD, ALIVE, DEAD, DEAD},
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
        char (*grid)[14][5] = &grid0;*/

        //heap allocation
        char (*grid)[totalRows][totalCols] = malloc((size_t)totalRows*(size_t)totalCols);
        initializeGrid(totalRows, totalCols, grid);

        /*//debug - stampa griglia
        printf("Sono il master e stampo l'intera griglia:\n");
        for (int row = 0; row < totalRows; row++) {
            for (int col = 0; col < totalCols; col++) {
                printf("%c", (*grid)[row][col]);
            }
            printf("\n");
        }*/

        int counts[nProcs];
        int displs[nProcs];

        //per ora solo splitTypeRows1D
        counts[0] = displs[0] = 0;
        for (int i = 0; i < nWorkers; i++) {
            counts[i+1] = (p[i].rows + p[i].topGhost + p[i].bottomGhost) * totalCols;
            displs[i+1] = (p[i].firstRow - p[i].topGhost) * totalCols;
        }

        /* //debug - counts/displs
        printf("Master counts:\n");
        for (int i = 0; i < nProcs; i++) {
            printf("%d, ", counts[i]);
        }
        printf("\n");

        printf("Master displs:\n");
        for (int i = 0; i < nProcs; i++) {
            printf("%d, ", displs[i]);
        }
        printf("\n");
        */

        MPI_Scatterv(grid, counts, displs, MPI_CHAR, NULL, 0, MPI_CHAR, MASTER, MPI_COMM_WORLD);

        counts[0] = displs[0] = 0;
        for (int i = 0; i < nWorkers; i++) {
            counts[i+1] = p[i].rows * totalCols;
            displs[i+1] = p[i].firstRow * totalCols;
        }

        MPI_Request gather_res_req;
        for(int i = 0; i < nIters; i++) {
            MPI_Igatherv(NULL, 0, MPI_CHAR, grid, counts, displs, MPI_CHAR, MASTER, MPI_COMM_WORLD, &gather_res_req);
            MPI_Wait(&gather_res_req, MPI_STATUS_IGNORE);

            /*
            //debug - stampa risultato iterazione
            printf("I'm the master process and I've received iteration %d:\n", i);
            for (int row = 0; row < totalRows; row++) {
                for (int col = 0; col < totalCols; col++) {
                    printf("%c", (*grid)[row][col]);
                }
                printf("\n");
            }
            fflush(stdout);
            */
        }
    } else { //worker code
        struct workPortion portion = p[myRank - FIRST_WORKER_RANK];

        MPI_Request gatherToMaster_req = MPI_REQUEST_NULL;
        MPI_Request neigh_rcv_reqs[2], neigh_snd_reqs[2];
        int n_neigh_rcv = 0;
        int n_neigh_snd = 0;

        /* /debug - worker's portion
        printf("I'm worker %d and this is my portion: firstRow: %d, firstCol: %d, nRows: %d, nCols: %d, left: %d, right: %d, bottom: %d, top: %d\n",
               myRank, portion.firstRow, portion.firstCol, portion.rows, portion.cols,
               portion.leftGhost, portion.rightGhost, portion.bottomGhost, portion.topGhost);
        */

        //aggiungi le ghost rows/cols se presenti
        int gridRows, gridCols;
        gridRows = portion.rows + portion.topGhost + portion.bottomGhost;
        gridCols = totalCols;

        //alloca lo spazio per le due griglie
        char (*currGrid)[gridRows][gridCols] = malloc((size_t)gridRows*(size_t)gridCols);
        char (*nextGrid)[gridRows][gridCols] = malloc((size_t)gridRows*(size_t)gridCols);

        //printf("I'm worker %d and this is my grid size: %d %d\n", myRank, gridRows, gridCols);

        //scatterv per ricevere porzione dal master
        MPI_Scatterv(NULL, NULL, NULL, NULL, currGrid, gridRows*gridCols, MPI_CHAR, MASTER, MPI_COMM_WORLD);
        
        //printf("I'm worker %d and I've received my portion from the master\n", myRank);
        //fflush(stdout);

        /*//debug - stampa porzione ricevuta
        printf("rcvside Process %d portion %d %d\n", myRank, gridRows, gridCols);
        for (int row = 0; row < gridRows; row++) {
            for (int col = 0; col < gridCols; col++) {
                printf("%c", (*currGrid)[row][col]);
            }
            printf("\n");
        }
        printf("\n");
        */
        
        for (int iter = 0; iter < nIters; iter++) {
            //Calcolo parziale in attesa delle ghost rows
            subGameOfLife(gridRows, gridCols, currGrid, nextGrid, portion.topGhost ? 1 : 0, portion.rows);

            MPI_Waitall(n_neigh_rcv, neigh_rcv_reqs, MPI_STATUSES_IGNORE);

            /*//debug - stampa porzioni dopo scambio di righe
            if (iter == 0) {
                printf("Process %d portion rows %d after row exchange\n", myRank, gridRows);
                for (int row = 0; row < gridRows; row++) {
                    for (int col = 0; col < gridCols; col++) {
                        printf("%c", (*nextGrid)[row][col]);
                    }
                    printf("\n");
                }
                printf("\n");
            }
            fflush(stdout);
            */

            //Calcola le righe e colonne restanti che dipendono dalle ghost rows
            if (portion.topGhost) {
                subGameOfLife(gridRows, gridCols, currGrid, nextGrid, 1, 1);
            }
            if (portion.bottomGhost) {
                subGameOfLife(gridRows, gridCols, currGrid, nextGrid, gridRows - 2, 1);
            }
            
            MPI_Waitall(n_neigh_snd, neigh_snd_reqs, MPI_STATUSES_IGNORE);

            /*//debug - stampa porzione processo 1 dopo calcolo
            if (rank == 1) {
                printf("calcolo1 Process %d portion %d %d\n", rank, portionRows, portionCols);
                for (int row = 0; row < portionRows; row++) {
                    for (int col = 0; col < TOTAL_COLS; col++) {
                        printf("%c", (*nextGrid)[row][col]);
                    }
                    printf("\n");
                }
                printf("\n");
            }
            */

            //controlla che l'invio del risultato dell'iterazione precedente sia completo
            //quindi il master ha ricevuto tutti i risultati fino a iter - 2 e sta completando solo iter - 1
            //a questo punto si può mandare iter
            //inoltre la griglia mandata in iter - 1 è riutilizzabile per computare iter + 1

            /* Collective operations can (but are not required to) complete as soon as the caller's
            participation in the collective communication is finished. A blocking operation is complete
            as soon as the call returns. A nonblocking (immediate) call requires a separate completion
            call (cf. Section 3.7). The completion of a collective operation indicates that the caller is free
            to modify locations in the communication buffer. It does not indicate that other processes
            in the group have completed or even started the operation (unless otherwise implied by the
            description of the operation). Thus, a collective communication operation may, or may not,
            have the effect of synchronizing all participating MPI processes. */
            MPI_Wait(&gatherToMaster_req, MPI_STATUS_IGNORE); //qui vuol dire che posso riutilizzare il buffer d'invio
            //igatherv del risultato al master, per ora solo splitTypeRows1D
            MPI_Igatherv(&(*nextGrid)[portion.topGhost ? 1 : 0][0], portion.rows * gridCols, MPI_CHAR,
                        NULL, NULL, NULL, NULL, MASTER, MPI_COMM_WORLD, &gatherToMaster_req);
            
            //printf("I'm worker %d and I've sent iteration %d to the master\n", myRank, iter);
            //fflush(stdout);

            n_neigh_rcv = n_neigh_snd = 0;
            if (portion.topGhost) {
                MPI_Irecv(&(*nextGrid)[0][0], gridCols, MPI_CHAR, myRank - 1, iter, MPI_COMM_WORLD, &neigh_rcv_reqs[n_neigh_rcv++]);
                MPI_Isend(&(*nextGrid)[1][0], gridCols, MPI_CHAR, myRank - 1, iter, MPI_COMM_WORLD, &neigh_snd_reqs[n_neigh_snd++]);
            }
            if (portion.bottomGhost) {
                MPI_Irecv(&(*nextGrid)[gridRows - 1][0], gridCols, MPI_CHAR, myRank + 1, iter, MPI_COMM_WORLD, &neigh_rcv_reqs[n_neigh_rcv++]);
                MPI_Isend(&(*nextGrid)[gridRows - 2][0], gridCols, MPI_CHAR, myRank + 1, iter, MPI_COMM_WORLD, &neigh_snd_reqs[n_neigh_snd++]);
            }

            //inverti griglie
            char (*temp)[][gridRows] = currGrid;
            currGrid = nextGrid;
            nextGrid = temp;
        }
    }

    printf("MPI time taken %.6lfs\n", MPI_Wtime() - t1);
    MPI_Finalize();
    return 0;
}

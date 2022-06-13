#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#define MASTER 0
#define ALIVE '1'
#define DEAD '0'

//commentare
struct work_portion {
    int first_row, rows;
    int top_ghost, bottom_ghost;
};

//commentare
int n_iters, total_rows, total_cols;

//commentare
void init_grid (int grid_rows, int grid_cols, char (*grid) [][grid_cols])
{
    time_t t;
    srand((unsigned) time(&t));
    for (int row = 0; row < grid_rows; row++) {
        for (int col = 0; col < grid_cols; col++) {
            if (rand() % 2 == 0) {
                (*grid)[row][col] = DEAD;
            } else {
                (*grid)[row][col] = ALIVE;
            }
        }
    }
}

/* Data una griglia M*N determina la miglior partizione della griglia per n_workers worker */
void define_work_portions (int grid_rows, int grid_cols, struct work_portion portions[], int n_portions)
{
    int rowsPerPortion = grid_rows / n_portions;
    int remainderRows = grid_rows % n_portions;

    for (int i = 0; i < n_portions; i++) {
        portions[i].top_ghost = portions[i].bottom_ghost = true;
        portions[i].first_row = rowsPerPortion * i; //determina prima riga della porzione i-esima
        portions[i].rows = rowsPerPortion; //determina quante righe per la porzione i-esima
    }

    portions[0].top_ghost = portions[n_portions - 1].bottom_ghost = false; //prima e ultima porzione hanno solo una ghost row
    portions[n_portions - 1].rows += remainderRows; //assegna le righe restanti all'ultima porzione
}

//commentare
void sub_game_of_life (int grid_rows, int grid_cols, char (*curr_grid) [][grid_cols], char (*next_grid) [][grid_cols], int first_game_row, int game_rows)
{
    //per ogni cella
    for (int row = first_game_row; row < first_game_row + game_rows; row++) {
        for (int col = 0; col < grid_cols; col++) {
            int alive_neighbours = 0; //conta i vicini ALIVE
            for (int row1 = row - 1; row1 <= row + 1; row1++) { //righe sopra-corrente-sotto
                for (int col1 = col - 1; col1 <= col + 1; col1++) { //colonne sinistra-corrente-destra
                    if (!(row1 == row && col1 == col)) { //non contare la cella stessa come aliveNeighbour
                        if (row1 >= 0 &&  col1 >= 0) { //check bound inferiore/sinistro griglia
                            if (row1 < grid_rows && col1 < grid_cols) { //check bound superiore/destro griglia
                                if ((*curr_grid)[row1][col1] == ALIVE) {
                                    alive_neighbours++;
                                }
                            }
                        }
                    }
                }
            }

            //calcola il nuovo stato delle celle
            (*next_grid)[row][col] = (*curr_grid)[row][col];
            if ((*curr_grid)[row][col] == ALIVE) {
                if (alive_neighbours < 2 || alive_neighbours > 3) {
                    (*next_grid)[row][col] = DEAD;
                }
            } else {
                if (alive_neighbours == 3) {
                    (*next_grid)[row][col] = ALIVE;
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    MPI_Init(NULL, NULL);
    double t1 = MPI_Wtime();

    int rank, n_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 4) {
        printf("Usage example: mpirun -np <n_procs> %s <n_rows> <n_cols> <n_iters>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    total_rows = atoi(argv[1]);
    total_cols = atoi(argv[2]);
    n_iters = atoi(argv[3]);

    int n_workers = n_procs - 1; //il master non lavora
    struct work_portion portions[n_workers];
    define_work_portions(total_rows, total_cols, portions, n_workers);
    
    //printf("args: nRows %d nCols %d n_iters %d\n", total_rows, total_cols, n_iters);

    if (rank == MASTER) {
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
        char (*grid)[total_rows][total_cols] = malloc((size_t)total_rows*(size_t)total_cols);
        init_grid(total_rows, total_cols, grid);

        /*//debug - stampa griglia
        printf("Sono il master e stampo l'intera griglia:\n");
        for (int row = 0; row < total_rows; row++) {
            for (int col = 0; col < total_cols; col++) {
                printf("%c", (*grid)[row][col]);
            }
            printf("\n");
        }*/

        int counts[n_procs], displs[n_procs];

        counts[0] = displs[0] = 0;
        for (int i = 0; i < n_workers; i++) {
            counts[i+1] = (portions[i].rows + portions[i].top_ghost + portions[i].bottom_ghost) * total_cols;
            displs[i+1] = (portions[i].first_row - portions[i].top_ghost) * total_cols;
        }

        /* //debug - counts/displs
        printf("Master counts:\n");
        for (int i = 0; i < n_procs; i++) {
            printf("%d, ", counts[i]);
        }
        printf("\n");

        printf("Master displs:\n");
        for (int i = 0; i < n_procs; i++) {
            printf("%d, ", displs[i]);
        }
        printf("\n");
        */

        MPI_Scatterv(grid, counts, displs, MPI_CHAR, NULL, 0, MPI_CHAR, MASTER, MPI_COMM_WORLD);

        counts[0] = displs[0] = 0;
        for (int i = 0; i < n_workers; i++) {
            counts[i+1] = portions[i].rows * total_cols;
            displs[i+1] = portions[i].first_row * total_cols;
        }

        MPI_Request gather_res_req;
        for(int i = 0; i < n_iters; i++) {
            MPI_Igatherv(NULL, 0, MPI_CHAR, grid, counts, displs, MPI_CHAR, MASTER, MPI_COMM_WORLD, &gather_res_req);
            MPI_Wait(&gather_res_req, MPI_STATUS_IGNORE);

            /* //debug - stampa risultato iterazione
            printf("I'm the master process and I've received iteration %d:\n", i);
            for (int row = 0; row < total_rows; row++) {
                for (int col = 0; col < total_cols; col++) {
                    printf("%c", (*grid)[row][col]);
                }
                printf("\n");
            }
            fflush(stdout);
            */
        }
    } else { //worker code
        struct work_portion portion = portions[rank - 1];

        MPI_Request gather_to_master_req = MPI_REQUEST_NULL;
        MPI_Request neigh_rcv_reqs[2], neigh_snd_reqs[2];
        int n_neigh_rcv = 0, n_neigh_snd = 0;

        //aggiungi le ghost rows/cols se presenti
        int grid_rows, grid_cols;
        grid_rows = portion.rows + portion.top_ghost + portion.bottom_ghost;
        grid_cols = total_cols;

        //alloca lo spazio per le due griglie
        char (*curr_grid)[grid_rows][grid_cols] = malloc((size_t)grid_rows*(size_t)grid_cols);
        char (*next_grid)[grid_rows][grid_cols] = malloc((size_t)grid_rows*(size_t)grid_cols);

        //scatterv per ricevere porzione dal master
        MPI_Scatterv(NULL, NULL, NULL, NULL, curr_grid, grid_rows*grid_cols, MPI_CHAR, MASTER, MPI_COMM_WORLD);
        
        for (int iter = 0; iter < n_iters; iter++) {
            //Calcolo parziale in attesa delle ghost rows
            sub_game_of_life(grid_rows, grid_cols, curr_grid, next_grid, portion.top_ghost ? 1 : 0, portion.rows);

            MPI_Waitall(n_neigh_rcv, neigh_rcv_reqs, MPI_STATUSES_IGNORE);

            //Calcola le righe e colonne restanti che dipendono dalle ghost rows
            if (portion.top_ghost) {
                sub_game_of_life(grid_rows, grid_cols, curr_grid, next_grid, 1, 1);
            }
            if (portion.bottom_ghost) {
                sub_game_of_life(grid_rows, grid_cols, curr_grid, next_grid, grid_rows - 2, 1);
            }
            
            MPI_Waitall(n_neigh_snd, neigh_snd_reqs, MPI_STATUSES_IGNORE);

            //wait della gather - qui vuol dire che posso riutilizzare il buffer d'invio
            MPI_Wait(&gather_to_master_req, MPI_STATUS_IGNORE); 
            //igatherv del risultato al master
            MPI_Igatherv(&(*next_grid)[portion.top_ghost ? 1 : 0][0], portion.rows * grid_cols, MPI_CHAR,
                        NULL, NULL, NULL, NULL, MASTER, MPI_COMM_WORLD, &gather_to_master_req);

            n_neigh_rcv = n_neigh_snd = 0;
            if (portion.top_ghost) {
                MPI_Irecv(&(*next_grid)[0][0], grid_cols, MPI_CHAR, rank - 1, iter, MPI_COMM_WORLD, &neigh_rcv_reqs[n_neigh_rcv++]);
                MPI_Isend(&(*next_grid)[1][0], grid_cols, MPI_CHAR, rank - 1, iter, MPI_COMM_WORLD, &neigh_snd_reqs[n_neigh_snd++]);
            }
            if (portion.bottom_ghost) {
                MPI_Irecv(&(*next_grid)[grid_rows - 1][0], grid_cols, MPI_CHAR, rank + 1, iter, MPI_COMM_WORLD, &neigh_rcv_reqs[n_neigh_rcv++]);
                MPI_Isend(&(*next_grid)[grid_rows - 2][0], grid_cols, MPI_CHAR, rank + 1, iter, MPI_COMM_WORLD, &neigh_snd_reqs[n_neigh_snd++]);
            }

            //inverti griglie
            char (*temp)[][grid_rows] = curr_grid;
            curr_grid = next_grid;
            next_grid = temp;
        }
    }

    printf("MPI time taken %.6lfs\n", MPI_Wtime() - t1);
    MPI_Finalize();
    return 0;
}

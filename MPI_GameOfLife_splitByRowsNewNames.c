#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#define MASTER 0
#define ALIVE '1'
#define DEAD '0'

//commentare
struct work_portion {
    int first_row, rows;
    bool has_top_ghost, has_bottom_ghost;
};

//commentare
int n_iters, tot_rows, tot_cols;

//commentare
void init_grid(int grid_rows, int grid_cols, char (*grid) [][grid_cols])
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
void define_work_portions(int grid_rows, int grid_cols, struct work_portion portions[], int n_portions)
{
    int rows_per_portion = grid_rows / n_portions;
    int remainder_rows = grid_rows % n_portions;

    for (int i = 0; i < n_portions; i++) {
        portions[i].has_top_ghost = portions[i].has_bottom_ghost = true;

        if (i == 0)
            portions[i].has_top_ghost = false;

        if (i == n_portions - 1)
            portions[i].has_bottom_ghost = false;

        portions[i].first_row = rows_per_portion * i - portions[i].has_top_ghost; //determina prima riga della porzione i-esima
        portions[i].rows = rows_per_portion + portions[i].has_top_ghost + portions[i].has_bottom_ghost; //determina quante righe per la porzione i-esima
    }

    portions[n_portions - 1].rows += remainder_rows; //assegna le righe restanti all'ultima porzione
}


int count_alive_neighbours(int grid_rows, int grid_cols, char (*curr_grid)[][grid_cols], int row, int col)
{
    int alive_neighbours = 0; //conta i vicini ALIVE
    for (int row1 = row - 1; row1 <= row + 1; row1++) { //righe sopra-corrente-sotto
        for (int col1 = col - 1; col1 <= col + 1; col1++) { //colonne sinistra-corrente-destra
            if (!(row1 == row && col1 == col)) { //non contare la cella stessa come aliveNeighbour
                if (row1 >= 0 &&  col1 >= 0) { //check bound superiore/sinistro griglia
                    if (row1 < grid_rows && col1 < grid_cols) { //check bound inferiore/destro griglia
                        if ((*curr_grid)[row1][col1] == ALIVE) {
                            alive_neighbours++;
                        }
                    }
                }
            }
        }
    }

    return alive_neighbours;
}

//commentare
void sub_game_of_life(int grid_rows, int grid_cols, char (*curr_grid)[][grid_cols], char (*next_grid) [][grid_cols], int first_game_row, int game_rows)
{
    //per ogni cella
    for (int row = first_game_row; row < first_game_row + game_rows; row++) {
        for (int col = 0; col < grid_cols; col++) {
            int alive_neighbours = count_alive_neighbours(grid_rows, grid_cols, curr_grid, row, col);

            //aggiorna la cella
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

    tot_rows = atoi(argv[1]);
    tot_cols = atoi(argv[2]);
    n_iters = atoi(argv[3]);

    int n_workers = n_procs - 1; //il master non lavora
    struct work_portion portions[n_workers];
    define_work_portions(tot_rows, tot_cols, portions, n_workers);
    
    //printf("args: nRows %d nCols %d n_iters %d\n", tot_rows, tot_cols, n_iters);

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
        char (*grid)[tot_rows][tot_cols] = malloc((size_t)tot_rows*tot_cols);
        init_grid(tot_rows, tot_cols, grid);

        /*//debug - stampa griglia
        printf("Sono il master e stampo l'intera griglia:\n");
        for (int row = 0; row < tot_rows; row++) {
            for (int col = 0; col < tot_cols; col++) {
                printf("%c", (*grid)[row][col]);
            }
            printf("\n");
        }*/

        int counts[n_procs], displs[n_procs];

        counts[0] = displs[0] = 0;
        for (int i = 0; i < n_workers; i++) {
            counts[i+1] = portions[i].rows * tot_cols;
            displs[i+1] = portions[i].first_row * tot_cols;
        }

        //Invia porzioni ai worker
        MPI_Scatterv(grid, counts, displs, MPI_CHAR, NULL, 0, MPI_CHAR, MASTER, MPI_COMM_WORLD);

        counts[0] = displs[0] = 0;
        for (int i = 0; i < n_workers; i++) {
            counts[i+1] = (portions[i].rows - portions[i].has_top_ghost - portions[i].has_bottom_ghost) * tot_cols;
            displs[i+1] = (portions[i].first_row + portions[i].has_top_ghost) * tot_cols;
        }

        //Ricevi risultati dai worker
        MPI_Request gather_res_req;
        for(int i = 0; i < n_iters; i++) {
            MPI_Igatherv(NULL, 0, MPI_CHAR, grid, counts, displs, MPI_CHAR, MASTER, MPI_COMM_WORLD, &gather_res_req);
            MPI_Wait(&gather_res_req, MPI_STATUS_IGNORE);

            /* //debug - stampa risultato iterazione
            printf("I'm the master process and I've received iteration %d:\n", i);
            for (int row = 0; row < tot_rows; row++) {
                for (int col = 0; col < tot_cols; col++) {
                    printf("%c", (*grid)[row][col]);
                }
                printf("\n");
            }
            fflush(stdout);
            */
        }
        free(grid);
    } else { //worker code
        struct work_portion portion = portions[rank - 1];

        MPI_Request gather_to_master_req = MPI_REQUEST_NULL;
        MPI_Request neigh_rcv_reqs[2], neigh_snd_reqs[2];
        int n_neigh_rcv = 0, n_neigh_snd = 0;

        //alloca lo spazio per le due griglie
        char (*curr_grid)[portion.rows][tot_cols];
        curr_grid = malloc((size_t)portion.rows*tot_cols);
        char (*next_grid)[portion.rows][tot_cols];
        next_grid = malloc((size_t)portion.rows*tot_cols);

        //scatterv per ricevere porzione dal master
        MPI_Scatterv(NULL, NULL, NULL, NULL, curr_grid, portion.rows*tot_cols, MPI_CHAR, MASTER, MPI_COMM_WORLD);
        
        for (int iter = 0; iter < n_iters; iter++) {
            //Calcolo parziale in attesa delle ghost rows
            sub_game_of_life(portion.rows, tot_cols,
                             curr_grid, next_grid,
                             portion.has_top_ghost ? 1 : 0,
                             portion.rows - portion.has_top_ghost - portion.has_bottom_ghost);

            //Aspetta le ghost rows
            MPI_Waitall(n_neigh_rcv, neigh_rcv_reqs, MPI_STATUSES_IGNORE);

            //Calcola le righe restanti che dipendono dalle ghost rows
            if (portion.has_top_ghost) {
                sub_game_of_life(portion.rows, tot_cols, curr_grid, next_grid, 1, 1);
            }
            if (portion.has_bottom_ghost) {
                sub_game_of_life(portion.rows, tot_cols, curr_grid, next_grid, portion.rows - 2, 1);
            }
            
            //Wait della gather - qui vuol dire che si può riutilizzare il buffer d'invio
            MPI_Wait(&gather_to_master_req, MPI_STATUS_IGNORE); 

            //Invio del risultato al master
            int first_row_to_send = portion.has_top_ghost ? 1 : 0;
            int rows_to_send = (portion.rows - portion.has_top_ghost - portion.has_bottom_ghost) * tot_cols;
            MPI_Igatherv(&(*next_grid)[first_row_to_send][0], rows_to_send, MPI_CHAR,
                         NULL, NULL, NULL, NULL, MASTER, MPI_COMM_WORLD, &gather_to_master_req);

            //Aspetta le send delle ghost row ai vicini
            MPI_Waitall(n_neigh_snd, neigh_snd_reqs, MPI_STATUSES_IGNORE);

            if (iter != n_iters - 1) {
                //Scambio ghost rows tra i vicini
                n_neigh_rcv = n_neigh_snd = 0;
                if (portion.has_top_ghost) {
                    MPI_Irecv(&(*next_grid)[0][0], tot_cols, MPI_CHAR, rank - 1, iter, MPI_COMM_WORLD, &neigh_rcv_reqs[n_neigh_rcv++]);
                    MPI_Isend(&(*next_grid)[1][0], tot_cols, MPI_CHAR, rank - 1, iter, MPI_COMM_WORLD, &neigh_snd_reqs[n_neigh_snd++]);
                }
                if (portion.has_bottom_ghost) {
                    MPI_Irecv(&(*next_grid)[portion.rows - 1][0], tot_cols, MPI_CHAR, rank + 1, iter, MPI_COMM_WORLD, &neigh_rcv_reqs[n_neigh_rcv++]);
                    MPI_Isend(&(*next_grid)[portion.rows - 2][0], tot_cols, MPI_CHAR, rank + 1, iter, MPI_COMM_WORLD, &neigh_snd_reqs[n_neigh_snd++]);
                }

                //Scambia puntatori griglie
                char (*temp)[][tot_cols] = curr_grid;
                curr_grid = next_grid;
                next_grid = temp;
            }
        }
        free(curr_grid);
        free(next_grid);
    }

    printf("MPI time taken %.6lfs\n", MPI_Wtime() - t1);
    MPI_Finalize();
    return 0;
}

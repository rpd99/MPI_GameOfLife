#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MASTER 0
#define ALIVE '1'
#define DEAD '0'

/**
 * Description
 * @member start_row
 * @member n_rows
 * @member has_top_ghost
 * @member has_bottom_ghost
 */
struct work_portion {
    int start_row, n_rows;
    bool has_top_ghost, has_bottom_ghost;
};

// commentare
int n_iters, tot_rows, tot_cols;

/**
 * Description
 * @param grid_rows
 * @param grid_cols
 * @param grid
 */
void init_grid(int grid_rows, int grid_cols, char (*grid)[][grid_cols]) {
    time_t t;
    srand((unsigned)time(&t));
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

/**
 * Data una griglia M*N partiziona la griglia per n_workers worker
 * @param grid_rows
 * @param grid_cols
 * @param portions
 * @param n_portions
 */
void define_work_portions(int grid_rows, int grid_cols, struct work_portion portions[], int n_portions) {
    int rows_per_portion = grid_rows / n_portions;
    int remainder_rows = grid_rows % n_portions;

    for (int i = 0; i < n_portions; i++) {
        portions[i].has_top_ghost = i > 0;                                                                 // la prima porzione non ha la top ghost row
        portions[i].has_bottom_ghost = i < n_portions - 1;                                                 // l'ultima porzione non ha la bottom ghost row
        portions[i].start_row = rows_per_portion * i - portions[i].has_top_ghost;                          // determina prima riga della porzione i-esima
        portions[i].n_rows = rows_per_portion + portions[i].has_top_ghost + portions[i].has_bottom_ghost;  // determina quante righe per la porzione i-esima

        // Distribuisci equamente le righe restanti tra i primi (grid_rows % n_portions) worker
        if (i < remainder_rows) {
            portions[i].n_rows += 1;
            portions[i].start_row += i;
        } else {
            portions[i].start_row += remainder_rows;
        }
    }
}

/**
 * Data una griglia M*N partiziona la griglia per n_workers worker
 * @param grid_rows
 * @param grid_cols
 * @param grid
 * @param row
 * @param col
 */
int count_alive_neighbours(int grid_rows, int grid_cols, char (*grid)[][grid_cols], int row, int col) {
    int alive_neighbours = 0;                        // conta i vicini ALIVE
    for (int r = row - 1; r <= row + 1; r++) {       // righe sopra-corrente-sotto
        for (int c = col - 1; c <= col + 1; c++) {   // colonne sinistra-corrente-destra
            if (!(r == row && c == col) &&           // non contare la cella stessa come aliveNeighbour
                (r >= 0 && c >= 0) &&                // check bound superiore/sinistro griglia
                (r < grid_rows && c < grid_cols) &&  // check bound inferiore/destro griglia
                ((*grid)[r][c] == ALIVE)) {
                alive_neighbours++;
            }
        }
    }

    return alive_neighbours;
}

/**
 * Data una griglia M*N partiziona la griglia per n_workers worker
 * @param grid_rows
 * @param grid_cols
 * @param curr_grid
 * @param next_grid
 * @param game_start_row
 * @param game_n_rows
 */
void sub_game_of_life(int grid_rows, int grid_cols, char (*curr_grid)[][grid_cols], char (*next_grid)[][grid_cols], int game_start_row, int game_n_rows) {
    // per ogni cella
    for (int row = game_start_row; row < game_start_row + game_n_rows; row++) {
        for (int col = 0; col < grid_cols; col++) {
            int alive_neighbours = count_alive_neighbours(grid_rows, grid_cols, curr_grid, row, col);

            // aggiorna la cella
            (*next_grid)[row][col] = (*curr_grid)[row][col];
            if ((*curr_grid)[row][col] == ALIVE) {
                if (alive_neighbours < 2 || alive_neighbours > 3)
                    (*next_grid)[row][col] = DEAD;
            } else {
                if (alive_neighbours == 3)
                    (*next_grid)[row][col] = ALIVE;
            }
        }
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(NULL, NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();

    if (argc != 4) {
        printf("Usage example: mpirun -np <n_procs> %s <n_rows> <n_cols> <n_iters>\n", argv[0]);
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    int rank, n_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    tot_rows = atoi(argv[1]);
    tot_cols = atoi(argv[2]);
    n_iters = atoi(argv[3]);

    int n_workers = n_procs - 1;  // il master non lavora
    struct work_portion portions[n_workers];
    define_work_portions(tot_rows, tot_cols, portions, n_workers);

    int scatter_counts[n_procs], scatter_displs[n_procs];
    int gather_counts[n_procs], gather_displs[n_procs];

    scatter_counts[0] = scatter_displs[0] = 0;
    gather_counts[0] = gather_displs[0] = 0;

    for (int i = 0; i < n_workers; i++) {
        // Scatter
        scatter_counts[i + 1] = portions[i].n_rows * tot_cols;
        scatter_displs[i + 1] = portions[i].start_row * tot_cols;
        // Gather
        gather_counts[i + 1] = (portions[i].n_rows - portions[i].has_top_ghost - portions[i].has_bottom_ghost) * tot_cols;
        gather_displs[i + 1] = (portions[i].start_row + portions[i].has_top_ghost) * tot_cols;
    }

    if (rank == MASTER) {  // Master code
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

        // heap allocation
        char(*grid)[tot_rows][tot_cols] = malloc((size_t)tot_rows * tot_cols);
        init_grid(tot_rows, tot_cols, grid);

        /*//debug - stampa griglia
        printf("Sono il master e stampo l'intera griglia:\n");
        for (int row = 0; row < tot_rows; row++) {
            for (int col = 0; col < tot_cols; col++) {
                printf("%c", (*grid)[row][col]);
            }
            printf("\n");
        }*/

        // Invia porzioni ai worker
        MPI_Scatterv(grid, scatter_counts, scatter_displs, MPI_CHAR, NULL, 0, MPI_CHAR, MASTER, MPI_COMM_WORLD);

        // Ricevi risultati dai worker
        MPI_Request gather_results_req;
        for (int i = 0; i < n_iters; i++) {
            MPI_Igatherv(NULL, 0, MPI_CHAR, grid, gather_counts, gather_displs, MPI_CHAR, MASTER, MPI_COMM_WORLD, &gather_results_req);
            MPI_Wait(&gather_results_req, MPI_STATUS_IGNORE);
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
    } else {  // worker code
        struct work_portion portion = portions[rank - 1];

        MPI_Request gather_to_master_req = MPI_REQUEST_NULL;
        MPI_Request neigh_rcv_reqs[2], neigh_snd_reqs[2];
        int n_neigh_rcv = 0, n_neigh_snd = 0;

        // alloca lo spazio per le due griglie
        char(*curr_grid)[portion.n_rows][tot_cols];
        curr_grid = malloc((size_t)portion.n_rows * tot_cols);
        char(*next_grid)[portion.n_rows][tot_cols];
        next_grid = malloc((size_t)portion.n_rows * tot_cols);

        // Indici delle top/ghost rows (se presenti nella porzione)
        int top_ghost_row_index = 0;
        int bot_ghost_row_index = portion.n_rows - 1;

        // Scatterv per ricevere porzione dal master
        MPI_Scatterv(NULL, NULL, NULL, NULL, curr_grid, scatter_counts[rank], MPI_CHAR, MASTER, MPI_COMM_WORLD);

        for (int iter = 0; iter < n_iters; iter++) {
            // Scambio ghost row tra i vicini
            n_neigh_rcv = n_neigh_snd = 0;
            if (portion.has_top_ghost) {
                MPI_Irecv(&(*curr_grid)[top_ghost_row_index][0], tot_cols, MPI_CHAR, rank - 1, iter, MPI_COMM_WORLD, &neigh_rcv_reqs[n_neigh_rcv++]);
                MPI_Isend(&(*curr_grid)[top_ghost_row_index + 1][0], tot_cols, MPI_CHAR, rank - 1, iter, MPI_COMM_WORLD, &neigh_snd_reqs[n_neigh_snd++]);
            }
            if (portion.has_bottom_ghost) {
                MPI_Irecv(&(*curr_grid)[bot_ghost_row_index][0], tot_cols, MPI_CHAR, rank + 1, iter, MPI_COMM_WORLD, &neigh_rcv_reqs[n_neigh_rcv++]);
                MPI_Isend(&(*curr_grid)[bot_ghost_row_index - 1][0], tot_cols, MPI_CHAR, rank + 1, iter, MPI_COMM_WORLD, &neigh_snd_reqs[n_neigh_snd++]);
            }

            // Calcolo parziale in attesa delle ghost rows
            int game_start_row = portion.has_top_ghost ? 2 : 0;
            int game_n_rows = portion.n_rows - (portion.has_top_ghost ? 2 : 0) - (portion.has_bottom_ghost ? 2 : 0);
            sub_game_of_life(portion.n_rows, tot_cols, curr_grid, next_grid, game_start_row, game_n_rows);

            // Aspetta le ghost rows
            MPI_Waitall(n_neigh_rcv, neigh_rcv_reqs, MPI_STATUSES_IGNORE);

            // Calcola le righe restanti che dipendono dalle ghost rows
            if (portion.has_top_ghost)
                sub_game_of_life(portion.n_rows, tot_cols, curr_grid, next_grid, top_ghost_row_index + 1, 1);

            if (portion.has_bottom_ghost)
                sub_game_of_life(portion.n_rows, tot_cols, curr_grid, next_grid, bot_ghost_row_index - 1, 1);

            // Wait della gather - qui vuol dire che si può riutilizzare il buffer d'invio
            MPI_Wait(&gather_to_master_req, MPI_STATUS_IGNORE);

            // Invio del risultato al master
            int first_row_to_send = portion.has_top_ghost ? 1 : 0;
            MPI_Igatherv(&(*next_grid)[first_row_to_send][0], gather_counts[rank], MPI_CHAR,
                         NULL, NULL, NULL, NULL, MASTER, MPI_COMM_WORLD, &gather_to_master_req);

            // Aspetta le send delle ghost row ai vicini
            MPI_Waitall(n_neigh_snd, neigh_snd_reqs, MPI_STATUSES_IGNORE);

            // Scambia puntatori griglie
            char(*temp)[][tot_cols] = curr_grid;
            curr_grid = next_grid;
            next_grid = temp;
        }
        MPI_Wait(&gather_to_master_req, MPI_STATUS_IGNORE);  // aspetta l'invio dell'ultimo risultato
        free(curr_grid);
        free(next_grid);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("MPI time taken %.6lfs\n", MPI_Wtime() - t1);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SOFTENING 1e-9f //Utilizzato per prevenire divergenze numeriche

typedef struct //Struttura dati che rappresenta la particella
{
    float x;
    float y;
    float z;
    float vx;
    float vy;
    float vz;
} Body;

int rank;              //rank del processo
int process;           //numero di processi
int iteration = 20;    //numero di iterazioni
int bodies;            //numero di particelle
const float dt = 0.1f; //time step
MPI_Datatype MPIbody;  //tipo di dato MPI
Body *collection;      //array di particelle

void bodyForce(Body *p, int start, int lenght);     //Si occupa di calcolare le variazioni della velocitÃ  e della posizione delle particelle
void randomizeBodies(Body *collection, int bodies); //Inizializza le particelle
void printBodies(Body *body);                       //Stampa la posizione delle particelle

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);                      //Inizializzo MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);        //identifica il processo corrente
    MPI_Comm_size(MPI_COMM_WORLD, &process);     //determina il numero di processi coinvolti nella computazione
    MPI_Type_contiguous(6, MPI_FLOAT, &MPIbody); //viene creato un tipo contiguo
    MPI_Type_commit(&MPIbody);                   // la commit viene fatta per poter utilizzare l'oggetto nelle comunicazioni

    bodies = atoi(argv[1]); //ottengo il numero di particelle

    collection = malloc(sizeof(Body) * bodies); //alloco l'insiemi delle particelle
    int part = bodies / process;                //calcolo il numero di particelle da assegnare ad ogni processo
    int rest = bodies % process;                //calcolo il resto
    int sum = 0;                                //sum viene utilizzato per la suddivisione in caso di resto
    int *dspl, *sc;                             //array di interi per gestire il resto
    dspl = malloc(sizeof(int) * process);       //alloco l'array che descrive lo spostamento
    sc = malloc(sizeof(int) * process);         //alloco l'array che descrive quanti elementi inviare
    double times;                               //utilizzato per calcolare il tempo di esecuzione

    times = MPI_Wtime(); //calcolo il tempo di esecuzione

    for (int i = 0; i < process; i++) //calcolo lo spostamento e quanti elementi inviare ad ogni processo
    {
        sc[i] = part; //al processo assegno un numero di particellle uguale alla parte calcolata
        if (rest > 0) //se c'Ã¨ resto gli assegno una particella in piÃ¹ e riduco il resto
        {
            sc[i]++;
            rest--;
        }
        dspl[i] = sum; //asseggno al processo lo scostamento
        sum += sc[i];  //aggiorno lo spostamento
    }

    randomizeBodies(collection, bodies); //vengono inizializzate le particelle

    for (int iter = 0; iter < iteration; iter++) //ad ogni iterazione i processi si occupano di calcolare le nuove posizioni e velocitÃ 
    {
        MPI_Bcast(collection, bodies, MPIbody, 0, MPI_COMM_WORLD); //invio l'array di particelle a tutti i processi

        bodyForce(collection, dspl[rank], sc[rank]); //ogni processo calcola la sua porzione

        MPI_Gatherv(&collection[dspl[rank]], sc[rank], MPIbody, collection, sc, dspl, MPIbody, 0, MPI_COMM_WORLD); //il master raccoglie i risultati aggiornando l'array
    }

    if (rank == 0) //il master stampa le posizioni finali ed il tempo impiegato
    {
        double timee = MPI_Wtime();
        double end = timee - times;
        printBodies(collection);
        printf("\nTime: %f\n", end);
        fflush(stdout);
    }

    MPI_Type_free(&MPIbody); //libera la memoria allocata dal tipo contiguo
    free(collection);        //dealloco l'array delle particelle
    free(sc);                //dealloco gli array dei contatori
    free(dspl);
    MPI_Finalize(); //termino l'esecuzione di MPI
    return 0;
}

void bodyForce(Body *p, int start, int lenght)
{
    for (int i = start; i < start + lenght; i++)
    {
        float Fx = 0.0f;
        float Fy = 0.0f;
        float Fz = 0.0f;

        for (int j = 0; j < bodies; j++)
        {
            float dx = p[j].x - p[i].x;
            float dy = p[j].y - p[i].y;
            float dz = p[j].z - p[i].z;
            float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
            float invDist = 1.0f / sqrtf(distSqr);
            float invDist3 = invDist * invDist * invDist;

            Fx += dx * invDist3;
            Fy += dy * invDist3;
            Fz += dz * invDist3;
        }
        p[i].vx += dt * Fx;
        p[i].vy += dt * Fy;
        p[i].vz += dt * Fz;
    }

    for (int l = start; l < start + lenght; l++)
    {
        p[l].x += p[l].vx * dt;
        p[l].y += p[l].vy * dt;
        p[l].z += p[l].vz * dt;
    }
}

void randomizeBodies(Body *collection, int bodies)
{
    for (int i = 0; i < bodies; i++)
    {
        collection[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        collection[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        collection[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        collection[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        collection[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        collection[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    }
}

void printBodies(Body *body)
{
    for (int i = 0; i < bodies; i++)
    {
        printf("\n %d\n", i);
        printf("x= %f\ty= %f\tz= %f\tvx= %f\tvy= %f\tvz= %f\t\n\n", body[i].x, body[i].y, body[i].z, body[i].vx, body[i].vy, body[i].vz);
        fflush(stdout);
    }
}

/*
Compilare mpicc nbody.c -lm -std=c99 -o nbody
Lanciare mpirun -np numero_processi nbody numero_particelle
*/

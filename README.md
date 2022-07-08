# N-Body
 A solution to the n-body problem using MPI
 
 
Nel problema N-body, abbiamo bisogno di trovare le posizioni e le velocità di un gruppo di particelle che interagiscono in un periodo di tempo. Ad esempio, come un astrofisico potrebbe voler conoscere le posizioni e le velocità di un gruppo di stelle. Allo stesso modo, un chimico potrebbe voler conoscere le posizioni e le velocità di un gruppo di molecole o atomi. 
Un N-body solver è un programma che trova la soluzione ad un problema N-body attraverso la simulazione del comportamento delle particelle. Gli input al problema sono la massa, la posizione e la velocità di ogni particella all’inizio della simulazione. Gli output tipicamente sono la posizione e la velocità di ogni particella in una sequenza di tempo specificata dall’utente, oppure semplicemente la posizione e la velocità di ogni particella alla fine del tempo specificato dall’utente.   


Configurazione: 
8 istanze EC2 t2.large Ubuntu Server 18.04 LTS (HVM), SSD Volume Type 
16 Processori  

Per la compilazione del sorgente utilizzare: 
mpicc nbody.c -lm -std=c99 -o nbody 
Per l'esecuzione: 
mpirun -np numero_processi nbody numero_particelle 
  



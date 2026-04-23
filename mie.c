/* ========================================================================== */                                                                        
/*   mie.c                                                                    */
/*   11 de Junio del 2016                                                     */
/*   Luis Adri�n Padilla Salas                                                */
/*                                                                            */
/*   Este programa es una rutina de Din�mica Molecular (con OMP)              */
/*   para un potencial efectivo tipo Mie(n, m)                                */
/*   con termostato y barostato de Berendsen                                  */
/*   correcciones de largo alcance y n dimensiones (1, 2 y 3).                */
/* ========================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Semilla para la generacion de numeros pseudo-aleatorios */
double seed = 0.0;
FILE *movie, *gr, *perfil;
                           
/* Declaracion de las funciones */
void crear_posiciones_uniformemente(int *nat, double *rho, int dofx, double rx[], double ry[], 
double rz[], double boxx, double boxy, double boxz, int ndivx, int ndivy,
int ndivz);

void crear_posiciones_centradas_rectangularmente_tercios(int *nat, double *rho, int dofx, double rx[], double ry[], // Coloca las partículas centradas en el centro de la caja de simulación
double rz[], double boxx, double boxy, double boxz, int ndivx, int ndivy, 
int ndivz);

void crear_posiciones_centradas_rectangularmente_mitad(int *nat, double *rho, int dofx, double rx[], double ry[],
double rz[], double boxx, double boxy, double boxz, int ndivx, int ndivy, 
int ndivz);

void crear_posiciones_mas_espacio_líquido(int *nat, double *rho, int dofx, double rx[], double ry[],
double rz[], double boxx, double boxy, double boxz, int ndivx, int ndivy, 
int ndivz);

void crear_velocidades(int nat, int dofx, double amasa, double temp,
double vx[], double vy[], double vz[]);

double gauss(void);

double aleatorio_uniforme(void);

void leer_configuracion(int *nat, double *rho, int dofx, double rx[], 
double ry[], double rz[], double vx[], double vy[], double vz[], double *boxx, 
double *boxy, double *boxz);

void guardar_configuracion(int nat, double boxx, double boxy, double boxz,
double rx[], double ry[], double rz[], double vx[], double vy[], double vz[]);

void fuerza(int nat, int dofx, double coef, double urc, double rcut, 
double *epi, double rx[], double ry[], double rz[], double boxx, double boxy, 
double boxz, double fx[], double fy[], double fz[], double expn, double expm, 
double *wxxr, double *wyyr, double *wzzr);

void lrc(double coef, double expn, double expm, double pi, double rho, 
double rcut, int dofx, double *ulr, double *plr);

void berendsen(int ncolectivo, int nat, int dofx, double time, double taut, 
double taup, double temp, double tempi, double presion, double presilr,
double vx[], double vy[], double vz[], double rx[], double ry[], double rz[],
double *boxx, double *boxy, double *boxz, double *rho);

void pelicula_gro(int iconf, int nat, int dofx, double boxx, double boxy, 
double boxz, double rx[], double ry[], double rz[]); 

void distribuciones(int dofx, int nat, double rho, double deltar, double rx[], 
double ry[], double rz[], double boxx, double boxy, double boxz, double pi,
int ndist, int nparticulas[], int npares[], int nconf, int iconf);

/* Inicio del programa principal */
int main(void)
{
double amasa, pi, coef, ntot, ladomenor, time, timeh, rho, temp, rcut;
double presion, taup, taut, sumrho, promrho;
double boxx, boxy, boxz, deltar, ulr, plr, urc, snrc, smrc, epi, eki, tempi;
double eti, presi, virk, virr, wxxr, wyyr, wzzr, eto, error;
double wxxv, wyyv, wzzv, pxx, pyy, pzz;
double epilr, etilr, presilr, sumek, sumep, sumt, sump, promek, promep;
double promp, promt, expn, expm;
int nat, dofx, nmovie, nprint, i, iconf, nperfil, ndist;
int ncolectivo, nopcion, nconfequi, nconf, ndivx, ndivy, ndivz, num_atomos;
int MAXNAT, MAXNG;
char trash[100];
FILE *in, *resumen, *todo, *presiones;       

amasa = 1.0;
pi = acos(-1.0);

/* Creamos los archivos de salida de la DM */
in = fopen( "in.dat", "r" );
resumen = fopen("resumen.dat", "w");
todo = fopen("todo.dat", "w");
presiones = fopen("presiones.dat", "w");

/* Leemos las variables del archivo de entrada in.dat    */
fscanf(in, "%d %s", &ncolectivo, trash);
fscanf(in, "%d %s", &nopcion, trash);
fscanf(in, "%d %s", &dofx, trash); // Número de dimensiones a trabajar
fscanf(in, "%lf %s", &expn, trash);
fscanf(in, "%lf %s", &expm, trash);
fscanf(in, "%lf %s", &time, trash);
fscanf(in, "%lf %lf %lf %s", &boxx, &boxy, &boxz, trash);
fscanf(in, "%d %d %d %s", &ndivx, &ndivy, &ndivz, trash); // Aquí se almacena la entrada de la distribución de partículas
fscanf(in, "%lf %s", &temp, trash);
fscanf(in, "%lf %s", &taut, trash);
fscanf(in, "%lf %s", &presion, trash);
fscanf(in, "%lf %s", &taup, trash);
fscanf(in, "%lf %s", &rcut, trash);
fscanf(in, "%d %s", &nconfequi, trash);
fscanf(in, "%d %s", &nconf, trash);
fscanf(in, "%d %s", &nperfil, trash);
fscanf(in, "%lf %s", &deltar, trash); // Perfil de densidad
fscanf(in, "%d %s", &nmovie, trash);
fscanf(in, "%d %s", &nprint, trash);   

ndist = round((nconf - nconfequi) / nperfil);

printf("\nPotencial Mie (%lf, %lf)\n", expn, expm);
fprintf(resumen,"Potencial Mie (%lf, %lf)\n", expn, expm);

if (ncolectivo == 1){
   printf("Colectivo NVE \n");
   fprintf(resumen, "Colectivo NVE \n");
   } 
else if (ncolectivo == 2){
   printf("Colectivo NVT \n");
   fprintf(resumen, "Colectivo NVT \n");
   }   
else {
   printf("Colectivo NPT \n");
   fprintf(resumen, "Colectivo NPT \n");
}
   
if (nopcion == 1) {
   printf("Corrida nueva \n");
   fprintf(resumen, "Corrida nueva \n");
   }
else if (nopcion == 2) {
   printf("Continuacion de corrida \n");
   fprintf(resumen, "Continuacion de corrida \n");
   }
      
printf("Movimientos a realizar: %d \n", nconf);   
fprintf(resumen, "Movimientos a realizar: %d \n", nconf);
printf("Movimientos para equilibrar: %d \n", nconfequi);
fprintf(resumen, "Movimientos para equilibrar: %d \n", nconfequi);
printf("Movimientos para promediar: %d \n", nconf-nconfequi);   
fprintf(resumen, "Movimientos para promediar: %d \n", nconf-nconfequi);

/* Se comprueba que la distancia de corte sea menor que la mitad del lado mas
peque�o de la caja */
if (dofx == 3){
   if (boxx < boxy){
   ladomenor = boxx;
   }
   else {
   ladomenor = boxy;
   }
   if (boxz < ladomenor){
   ladomenor = boxz;
   }
   MAXNAT = 1 + ndivx * ndivy * ndivz;
   MAXNG = 1 + (int)(boxx / deltar);
}
else if (dofx == 2){
   if (boxx < boxy){
   ladomenor = boxx;
   }
   else {
   ladomenor = boxy;
   }
   boxz = 0.0;
   ndivz = 0;
   MAXNAT = 1 + ndivx * ndivy;
   MAXNG = 1 + (int)(boxx / deltar);
}
else if (dofx == 1){
   ladomenor = boxx;
   boxy = 0.0;
   ndivy = 0;
   boxz = 0.0;
   ndivz = 0;
   MAXNAT = 1 + ndivx;
   MAXNG = 1 + (int)(boxx / deltar);
}
   
double rx[MAXNAT], ry[MAXNAT], rz[MAXNAT];
double vx[MAXNAT], vy[MAXNAT], vz[MAXNAT];
double fx[MAXNAT], fy[MAXNAT], fz[MAXNAT];
int nparticulas[MAXNG], npares[MAXNG];
    
for(i=1 ; i <= MAXNG ; i++){
   nparticulas[i] = 0;
   npares[i] = 0;
}

if (rcut >= (0.5*ladomenor)){
   printf("ERROR: La distancia de corte es mayor que L/2... \n");
   exit(0);
}

/* Creamos las posiciones y velocidades iniciales si es primer corrida */   
if (nopcion == 1){
   crear_posiciones_mas_espacio_líquido(&nat, &rho, dofx, rx, ry, rz, boxx, boxy, boxz, ndivx, 
   ndivy, ndivz);
   crear_velocidades(nat, dofx, amasa, temp, vx, vy, vz);
}
/* Leemos el archivo de configuraciones si es continuacion de corrida */
else if (nopcion == 2){
   leer_configuracion(&nat, &rho, dofx, rx, ry, rz, vx, vy, vz, &boxx, &boxy, 
   &boxz);
}

printf("Numero de particulas: %d \n", nat);
printf("Dimensiones de la caja: %lf %lf %lf \n", boxx, boxy, boxz);

// Nuevo: Se especifican las condiciones de la minibox en el archivo de resumen
//fprintf(resumen, "Dimensiones de la minibox: %lf %lf %lf",  )

fprintf(resumen,"Numero de particulas: %d \n", nat);
fprintf(resumen,"Dimensiones de la caja: %lf %lf %lf \n", boxx, boxy, boxz);
fprintf(resumen,"Densidad inicial reducida: %lf \n", rho);

/* Calculamos las correcciones de largo alcance */
coef = (1.0 * expn/(expn-expm)) * pow(1.0 * expn/expm, 1.0 * expm/(expn-expm));
printf("Coef: %lf\n",coef);
fprintf(resumen, "Coef: %lf\n",coef);
lrc(coef, expn, expm, pi, rho, rcut, dofx, &ulr, &plr); 

printf("Ulr: %lf \n", ulr);
printf("Plr: %lf \n", plr);
fprintf(resumen, "Ulr: %.6lf \n", ulr);
fprintf(resumen, "Plr: %.6lf \n", plr);

printf("    iconf   rho      eki      epi      etot     tempi    presi    error\n");

/* Calculamos el potencial en la rcut */
snrc = pow((1.0 / rcut), expn);
smrc = pow((1.0 / rcut), expm);
urc = coef * (snrc - smrc);

/* Calculamos la fuerza inicial sobre las part�culas a t=0 */
fuerza(nat, dofx, coef, urc, rcut, &epi, rx, ry, rz, boxx, boxy, boxz, fx, fy,
fz, expn, expm, &wxxr, &wyyr, &wzzr);

/* Iniciamos variables para promediar */
sumek = 0.0;
sumep = 0.0;
sumt = 0.0;
sump = 0.0;
sumrho = 0.0;

/* Movemos las particulas con el algoritmo de Verlet de velocidades */
timeh = 0.5 * time;

for(iconf=1 ; iconf<=nconf ; iconf++){

   /* Grabamos las posiciones iniciales en un archivo de visualizacion gro */
   if (iconf == 1){
   pelicula_gro(iconf, nat, dofx, boxx, boxy, boxz, rx, ry, rz); 
   }

   /* Calculamos las velocidades al tiempo dt/2 */
   for(i=1 ; i <= nat ; i++){
      vx[i] = vx[i] + timeh * fx[i];
      vy[i] = vy[i] + timeh * fy[i];
      vz[i] = vz[i] + timeh * fz[i];
   }
   
   /* Calculamos las posiciones al tiempo dt */
   for(i=1 ; i <= nat ; i++){
      rx[i] = rx[i] + time * vx[i];
      ry[i] = ry[i] + time * vy[i];
      rz[i] = rz[i] + time * vz[i];
   }
   
   /* Aplicamos condiciones periodicas */
   if (dofx == 3){
      for(i=1 ; i <= nat ; i++){
         rx[i] = rx[i] - round(rx[i] / boxx) * boxx;
         ry[i] = ry[i] - round(ry[i] / boxy) * boxy;
         rz[i] = rz[i] - round(rz[i] / boxz) * boxz;
      }
   }
   else if (dofx == 2){
      for(i=1 ; i <= nat ; i++){
         rx[i] = rx[i] - round(rx[i] / boxx) * boxx;
         ry[i] = ry[i] - round(ry[i] / boxy) * boxy;
         rz[i] = 0.0;
      }
   }
   else if (dofx == 1){
      for(i=1 ; i <= nat ; i++){
         rx[i] = rx[i] - round(rx[i] / boxx) * boxx;
         ry[i] = 0.0;
         rz[i] = 0.0;
      }
   }
   
   /* Calculamos la fuerza al tiempo dt */
   fuerza(nat, dofx, coef, urc, rcut, &epi, rx, ry, rz, boxx, boxy, boxz, fx,
   fy, fz, expn, expm, &wxxr, &wyyr, &wzzr);

   /* Calculamos las velocidades al tiempo dt */
   for(i=1 ; i <= nat ; i++){
      vx[i] = vx[i] + timeh * fx[i];
      vy[i] = vy[i] + timeh * fy[i];
      vz[i] = vz[i] + timeh * fz[i];
   }

   /* Calculamos la temperatura y energia cinetica instantanea */
   eki = 0.0;
   wxxv = 0.0;
   wyyv = 0.0;
   wzzv = 0.0;
   for(i=1 ; i <= nat ; i++){
      eki = eki + (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
      wxxv = wxxv + vx[i] * vx[i];
      wyyv = wyyv + vy[i] * vy[i];
      wzzv = wzzv + vz[i] * vz[i];
   }
   virk = eki;
   virr = wxxr + wyyr + wzzr;
   eki = 0.5 * eki / nat;
   epi = epi / nat;
   eti = eki + epi;
   tempi = 2.0 * eki / dofx;
   
   /* Calculamos el error relativo en la energ�a total */
   if (iconf == 1){
      eto = eti;
   }
   else if (iconf == nconfequi + 1){
      eto = eti;
   }
   error = (eti - eto) / eto;
   error = fabs(error);
   
   /* Calculamos la presion como la ideal mas la de exceso debido a las 
   interacciones y las presiones en cada direccion */   
   presi = (virk + virr) * rho / (nat * dofx);
   pxx = (wxxr + wxxv) * rho /(nat);
   pyy = (wyyr + wyyv) * rho /(nat);
   pzz = (wzzr + wzzv) * rho /(nat);
   fprintf(presiones,"%9d %10.6lf %10.6lf %10.6lf\n",iconf, pxx, pyy, pzz);
    
   /* Calculamos las correcciones de largo alcance para NPT en cada paso */   
   if(ncolectivo == 3){
      lrc(coef, expn, expm, pi, rho, rcut, dofx, &ulr, &plr); 
   }
      
   /* Sumamos las correcciones de largo alcance */
   epilr = epi + ulr;
   etilr = eti + ulr;
   presilr = presi + plr;
   
   /* Imprimimos los resultados en pantalla */
   if (fmod(iconf, nprint) == 0){
   printf("%9d %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf\n", iconf, rho,
   eki, epilr, etilr, tempi, presilr, error);
   }
   /* Imprimimos todos los resultados en el archivo todo.dat */
   fprintf(todo,"%9d %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf\n",
   iconf, rho, eki, epilr, etilr, tempi, presilr, error);
   
   /* Acumulamos las propiedades para promediar */
   if (iconf > nconfequi){
   sumek = sumek + eki;
   sumep = sumep + epilr;
   sumt = sumt + tempi;
   sump = sump + presilr;
   sumrho = sumrho + rho;
   }
   
   /* Aplicamos el termostato y barostato de Berendsen para el siguiente 
   movimiento en caso de usar el colectivo NVT o NPT */
   if (ncolectivo == 2 || ncolectivo == 3){
      berendsen(ncolectivo, nat, dofx, time, taut, taup, temp, tempi, presion,
      presilr, vx, vy, vz, rx, ry, rz, &boxx, &boxy, &boxz, &rho);
   }
      
   /* Grabamos las nuevas posiciones en un archivo de visualizacion gro */
   if (fmod(iconf, nmovie) == 0){
   pelicula_gro(iconf, nat, dofx, boxx, boxy, boxz, rx, ry, rz); 
   }
   
   /* Calculamos rho(r) y g(r) varias veces para promediarlas */
   if (iconf > nconfequi){
      if (fmod(iconf-nconfequi, nperfil) == 0){
         distribuciones(dofx, nat, rho, deltar, rx, ry, rz, boxx, boxy, boxz, 
         pi, ndist, nparticulas, npares, nconf, iconf);
      }
   }

}

/* Calculamos los promedios */
promek = sumek / (nconf - nconfequi);
promep = sumep / (nconf - nconfequi);
promt = sumt / (nconf - nconfequi);
promp = sump / (nconf - nconfequi);
promrho = sumrho / (nconf - nconfequi);

/* Imprimimos los promedios en pantalla y archivo de salida */
printf("\nPROMEDIOS\n\n");

printf("Densidad: %.6lf \n", promrho);
printf("Energia Cinetica: %.6lf \n", promek);
printf("Energia Potencial: %.6lf \n", promep);
printf("Temperatura: %.6lf \n", promt);
printf("Presion: %.6lf \n", promp);

fprintf(resumen, "\nPROMEDIOS\n\n");
fprintf(resumen, "Densidad: %.6lf \n", promrho);
fprintf(resumen, "Energia Cinetica: %.6lf \n", promek);
fprintf(resumen, "Energia Potencial: %.6lf \n", promep);
fprintf(resumen, "Temperatura: %.6lf \n", promt);
fprintf(resumen, "Presion: %.6lf \n", promp);

/* Guardamos la configuraci�n final de las particulas */
guardar_configuracion(nat, boxx, boxy, boxz, rx, ry, rz, vx, vy, vz);  

return 0;
}

/* Fin del programa principal */

/* Funciones del programa principal */
/* Creamos las posiciones iniciales en la caja */
void crear_posiciones_uniformemente(int *nat, // El número total de partículas ingresado
   double *rho, int dofx, double rx[], double ry[], 
double rz[], double boxx, double boxy, double boxz, int ndivx, int ndivy, 
int ndivz) {
   double volumen, area, longitud, deltax, deltay, deltaz;
   int i, j, k, n;

   printf("Creando configuracion inicial...\n");
   if (dofx == 3){
      /* Para 3D */
      /* Calculamos la densidad y el numero de particulas */
      volumen = boxx * boxy * boxz;
      *nat = ndivx * ndivy * ndivz;   
      *rho = (*nat) / volumen;
      printf("Densidad volumetrica reducida: %lf \n", *rho);
      /* Calculamos la separacion entre particulas */
      deltax = boxx / ndivx;
      deltay = boxy / ndivy;         
      deltaz = boxz / ndivz;
      /* Creamos las posiciones iniciales en la caja */  
      n = 0;
      for(i=1; i<=ndivx ; i++){
         for(j=1; j<=ndivy ; j++){
            for(k=1; k<=ndivz ; k++){
               n += 1;
               rx[n] = (i-1) * deltax + 0.5;
               ry[n] = (j-1) * deltay + 0.5;
               rz[n] = (k-1) * deltaz + 0.5;
            }
         }
      }
      
      /* Centramos las posiciones simetricamente respecto al origen */
      for(i=1; i<=(*nat) ; i++){
      rx[i] = rx[i] - 0.5 * boxx;
      ry[i] = ry[i] - 0.5 * boxy;
      rz[i] = rz[i] - 0.5 * boxz;
      }   
   }
   else if (dofx == 2){
      /* Para 2D */
      /* Calculamos la densidad y el numero de particulas */
      area = boxx * boxy;
      *nat = ndivx * ndivy;   
      *rho = (*nat) / area;
      printf("Densidad superficial reducida: %lf \n", *rho);
      /* Calculamos la separacion entre particulas */
      deltax = boxx / ndivx;
      deltay = boxy / ndivy;         
      /* Creamos las posiciones iniciales en un cuadrado */
      n = 0;
      for(i=1; i<=ndivx ; i++){
         for(j=1; j<=ndivy ; j++){
            n += 1;
            rx[n] = (i-1) * deltax + 0.5;
            ry[n] = (j-1) * deltay + 0.5;
         }
      }
      
      /* Centramos las posiciones simetricamente respecto al origen */
      for(i=1; i<=(*nat) ; i++){
      rx[i] = rx[i] - 0.5 * boxx;
      ry[i] = ry[i] - 0.5 * boxy;
      rz[i] = 0.0;
      }   
   }
   else if (dofx == 1){
      /* Para 1D */
      /* Calculamos la densidad y el numero de particulas */
      longitud = boxx;
      *nat = ndivx;   
      *rho = (*nat) / longitud;
      printf("Densidad lineal reducida: %lf \n", *rho);
      /* Calculamos la separacion entre particulas */
      deltax = boxx / ndivx;         
      /* Creamos las posiciones iniciales en una linea */
      n = 0;
      for(i=1; i<=ndivx ; i++){
         n += 1;
         rx[n] = (i-1) * deltax + 0.5;
      }
      
      /* Centramos las posiciones simetricamente respecto al origen */
      for(i=1; i<=(*nat) ; i++){
      rx[i] = rx[i] - 0.5 * boxx;
      ry[i] = 0.0;
      rz[i] = 0.0;
      }   
   }

   return;
}

void crear_posiciones_centradas_rectangularmente_tercios(int *nat, double *rho, int dofx, double rx[], double ry[], // Coloca las partículas centradas en el centro de la caja de simulación
double rz[], double boxx, double boxy, double boxz, int ndivx, int ndivy, 
int ndivz) { 
   double volumen_caja, longitud, area, deltax, deltay, deltaz;
   double volumen_minibox, minibox_rho, area_minibox, longitud_minibox, minibox_x;
   int i, j, k, n; 

   printf("Creando configuración inicial rectangularmente centrada...\n");
   if (dofx == 3) {
      // Para 3D
      // Calculamos el volumen total de la caja
      volumen_caja = boxx * boxy * boxz; 
      *nat = ndivx * ndivy * ndivz; // Se obtiene el número de átomos
      *rho = (*nat) / volumen_caja;    
      printf("Densidad volumétrica reducida global: %lf \n", *rho);

      // En esta version se genera un volumen interno para distribuir las partículas.
      minibox_x = boxx / 3.0;
      volumen_minibox = minibox_x * boxy * boxz;           
      minibox_rho = (*nat) / volumen_minibox;
      printf("Densidad volumetrica reducida de la mini caja centrada: %lf \n", minibox_rho);

      // Calculamos la distancia entre partículas
      deltax = minibox_x / ndivx;
      deltay = boxy / ndivy;
      deltaz = boxz / ndivz;

      // Creamos las posiciones iniciales de las particulas dentro de la mini-box
      n = 0;
      for(i=1; i<=ndivx ; i++){
         for(j=1; j<=ndivy ; j++){
            for(k=1; k<=ndivz ; k++){
               n += 1;
               rx[n] = (i-1) * deltax + 0.5;
               ry[n] = (j-1) * deltay + 0.5;
               rz[n] = (k-1) * deltaz + 0.5;
            }
         }
      }

      // Centramos las posiciones simetricamente respecto al centro de la caja
      for(i=1; i<=(*nat) ; i++){
         rx[i] = rx[i] - 0.5 * minibox_x;
         ry[i] = ry[i] - 0.5 * boxy;
         rz[i] = rz[i] - 0.5 * boxz;
      }  
   }

   else if (dofx == 2) { 
      // Para 2D
      // calculamos la densidad y número de partículas
      area = boxx * boxy;
      *nat = ndivx * ndivy;   
      *rho = (*nat) / area;
      printf("Densidad superficial reducida: %lf \n", *rho);

      minibox_x = boxx / 3.0;
      area_minibox = minibox_x * boxy;
      minibox_rho = (*nat) / area_minibox;
      printf("Densidad superficial reducida: %lf \n", minibox_rho);
      /* Calculamos la separacion entre particulas */
      deltax = minibox_x / ndivx;
      deltay = boxy / ndivy;

      /* Creamos las posiciones iniciales en un cuadrado */
      n = 0;
      for(i=1; i<=ndivx ; i++){
         for(j=1; j<=ndivy ; j++){
            n += 1;
            rx[n] = (i-1) * deltax + 0.5;
            ry[n] = (j-1) * deltay + 0.5;
         }
      }
      
      /* Centramos las posiciones simetricamente respecto al origen */
      for(i=1; i<=(*nat) ; i++){
         rx[i] = rx[i] - 0.5 * boxx;
         ry[i] = ry[i] - 0.5 * boxy;
         rz[i] = 0.0;
      }   
   }
   
   else if (dofx == 1){
      /* Para 1D */
      /* Calculamos la densidad y el numero de particulas */
      longitud = boxx;
      *nat = ndivx;   
      *rho = (*nat) / longitud;
      printf("Densidad lineal reducida: %lf \n", *rho);

      // Calculamos la sección de generación de partículas
      longitud_minibox = boxx / 3.0;

      /* Calculamos la separacion entre particulas */
      deltax = longitud_minibox / ndivx;

      /* Creamos las posiciones iniciales en una linea */
      n = 0;
      for(i=1; i<=ndivx ; i++){
         n += 1;
         rx[n] = (i-1) * deltax + 0.5;
      }
      
      /* Centramos las posiciones simetricamente respecto al origen */
      for(i=1; i<=(*nat) ; i++){
         rx[i] = rx[i] - 0.5 * boxx;
         ry[i] = 0.0;
         rz[i] = 0.0;
      }   
   }

   return;
}

void crear_posiciones_centradas_rectangularmente_mitad(int *nat, double *rho, int dofx, double rx[], double ry[],
double rz[], double boxx, double boxy, double boxz, int ndivx, int ndivy, 
int ndivz) {  // Esta función permite tener relaciones de aspecto 2:1 o menores a 3:1 sin que las energías potenciales se disparen
   double volumen_caja, longitud, area, deltax, deltay, deltaz;
   double volumen_minibox, minibox_rho, area_minibox, longitud_minibox, minibox_x;
   int i, j, k, n; 

   printf("Creando configuración centrada (dividida en 2)... \n");
   if (dofx == 3) {
      // Para 3D
      // Calculamos el volumen total de la caja
      volumen_caja = boxx * boxy * boxz; 
      *nat = ndivx * ndivy * ndivz; // Se obtiene el número de átomos
      *rho = (*nat) / volumen_caja;    
      printf("Densidad volumétrica reducida global: %lf \n", *rho);

      // En esta version se genera un volumen interno para distribuir las partículas.
      minibox_x = boxx / 2.0; // La caja se divide en 4 partes sobe el eje x
      volumen_minibox = minibox_x * boxy * boxz;  // La minibox debe tener dimensiones iguales en los 3 ejes
      minibox_rho = (*nat) / volumen_minibox;
      printf("Densidad volumetrica reducida de la mini caja centrada: %lf \n", minibox_rho);

      // Calculamos la distancia entre partículas
      deltax = minibox_x / ndivx;
      deltay = boxy / ndivy;
      deltaz = boxz / ndivz;

      // Creamos las posiciones iniciales de las particulas dentro de la mini-box
      n = 0;
      for(i=1; i<=ndivx ; i++){
         for(j=1; j<=ndivy ; j++){
            for(k=1; k<=ndivz ; k++){
               n += 1;
               rx[n] = (i-1) * deltax + 0.5;
               ry[n] = (j-1) * deltay + 0.5;
               rz[n] = (k-1) * deltaz + 0.5;
            }
         }
      }

      // Centramos las posiciones simetricamente respecto al centro de la caja
      for(i=1; i<=(*nat) ; i++){
         rx[i] = rx[i] - 0.5 * minibox_x;
         ry[i] = ry[i] - 0.5 * boxy;
         rz[i] = rz[i] - 0.5 * boxz;
      }  
   }

   else if (dofx == 2) { 
      // Para 2D
      // calculamos la densidad y número de partículas
      area = boxx * boxy;
      *nat = ndivx * ndivy;   
      *rho = (*nat) / area;
      printf("Densidad superficial reducida: %lf \n", *rho);

      minibox_x = boxx / 4.0;
      area_minibox = (2.0 * minibox_x) * boxy;
      minibox_rho = (*nat) / area_minibox;
      printf("Densidad superficial reducida: %lf \n", minibox_rho);
      /* Calculamos la separacion entre particulas */
      deltax = minibox_x / ndivx;
      deltay = boxy / ndivy;

      /* Creamos las posiciones iniciales en un cuadrado */
      n = 0;
      for(i=1; i<=ndivx ; i++){
         for(j=1; j<=ndivy ; j++){
            n += 1;
            rx[n] = (i-1) * deltax + 0.5;
            ry[n] = (j-1) * deltay + 0.5;
         }
      }
      
      /* Centramos las posiciones simetricamente respecto al origen */
      for(i=1; i<=(*nat) ; i++){
         rx[i] = rx[i] - 0.5 * boxx;
         ry[i] = ry[i] - 0.5 * boxy;
         rz[i] = 0.0;
      }   
   }
   
   else if (dofx == 1){
      /* Para 1D */
      /* Calculamos la densidad y el numero de particulas */
      longitud = boxx;
      *nat = ndivx;   
      *rho = (*nat) / longitud;
      printf("Densidad lineal reducida: %lf \n", *rho);

      // Calculamos la sección de generación de partículas
      longitud_minibox = boxx / 4.0;

      /* Calculamos la separacion entre particulas */
      deltax = longitud_minibox / ndivx;

      /* Creamos las posiciones iniciales en una linea */
      n = 0;
      for(i=1; i<=ndivx ; i++){
         n += 1;
         rx[n] = (i-1) * deltax + 0.5;
      }
      
      /* Centramos las posiciones simetricamente respecto al origen */
      for(i=1; i<=(*nat) ; i++){
         rx[i] = rx[i] - 0.5 * boxx;
         ry[i] = 0.0;
         rz[i] = 0.0;
      }   
   }

   return;
}

void crear_posiciones_mas_espacio_líquido(int *nat, double *rho, int dofx, double rx[], double ry[],
double rz[], double boxx, double boxy, double boxz, int ndivx, int ndivy, 
int ndivz) {
   double volumen_caja, longitud, area, deltax, deltay, deltaz, divisiones_box, secciones_para_minibox;
   double volumen_minibox, minibox_rho, area_minibox, longitud_minibox, minibox_x;
   int i, j, k, n; 

   divisiones_box = 6.0; // <-- Señala el número de divisiones que quieres para la caja 
   secciones_para_minibox = 4.0; // Cuántas de esas divisiones corresponden a la minibox

   printf("Creando configuración centrada (dividida en %lf)... \n", divisiones_box);
   if (dofx == 3) {
      // Para 3D
      // Calculamos el volumen total de la caja
      volumen_caja = boxx * boxy * boxz; 
      *nat = ndivx * ndivy * ndivz; // Se obtiene el número de átomos
      *rho = (*nat) / volumen_caja;    
      printf("Densidad volumétrica reducida global: %lf \n", *rho);

      // En esta version se genera un volumen interno para distribuir las partículas.
      minibox_x = boxx / divisiones_box; // La caja se divide en partes sobe el eje x
      volumen_minibox = (minibox_x * secciones_para_minibox) * boxy * boxz;  // Se toma la sección de la caja en la que se colocarán las partículas (minibox)
      minibox_rho = (*nat) / volumen_minibox;
      printf("Densidad volumetrica reducida de la mini caja centrada: %lf \n", minibox_rho);

      // Calculamos la distancia entre partículas
      deltax = minibox_x / ndivx;
      deltay = boxy / ndivy;
      deltaz = boxz / ndivz;

      // Creamos las posiciones iniciales de las particulas dentro de la mini-box
      n = 0;
      for(i=1; i<=ndivx ; i++){
         for(j=1; j<=ndivy ; j++){
            for(k=1; k<=ndivz ; k++){
               n += 1;
               rx[n] = (i-1) * deltax + 0.5;
               ry[n] = (j-1) * deltay + 0.5;
               rz[n] = (k-1) * deltaz + 0.5;
            }
         }
      }

      // Centramos las posiciones simetricamente respecto al centro de la caja
      for(i=1; i<=(*nat) ; i++){
         rx[i] = rx[i] - 0.5 * minibox_x;
         ry[i] = ry[i] - 0.5 * boxy;
         rz[i] = rz[i] - 0.5 * boxz;
      }  
   }

   return;
}

/* Creamos las velocidades iniciales de acuerdo a una distribucion gaussiana */
void crear_velocidades(int nat, int dofx, double amasa, double temp,
double vx[], double vy[], double vz[])
{
int i;
double rtemp, sumx, sumy, sumz, uk, tempi, dof, prueba;

/* Asignamos velocidades a las part�culas de acuerdo a la equipartici�n de
la energ�a y con distribuci�n gaussiana */
if (dofx == 3){
   for(i=1; i<=nat ; i++){
      rtemp = sqrt(temp / amasa);
      vx[i] = rtemp * gauss();
      vy[i] = rtemp * gauss(); 
      vz[i] = rtemp * gauss();
   }
}
else if (dofx == 2){
   for(i=1; i<=nat ; i++){
      rtemp = sqrt(temp / amasa);
      vx[i] = rtemp * gauss();
      vy[i] = rtemp * gauss(); 
      vz[i] = 0.0;
   }
}
else if (dofx == 1){
   for(i=1; i<=nat ; i++){
      rtemp = sqrt(temp / amasa);
      vx[i] = rtemp * gauss();
      vy[i] = 0.0; 
      vz[i] = 0.0;
   }
} 
   
/* Removemos el momento del centro de masa */
sumx = 0.0;
sumy = 0.0;
sumz = 0.0;
for(i=1; i<=nat ; i++){
   sumx = sumx + vx[i];
   sumy = sumy + vy[i];
   sumz = sumz + vz[i];
}
   
sumx = sumx / nat;
sumy = sumy / nat;
sumz = sumz / nat;
uk = 0.0;
for(i=1; i<=nat ; i++){
   vx[i] = vx[i] - sumx;
   vy[i] = vy[i] - sumy;
   vz[i] = vz[i] - sumz;
   uk = uk + amasa * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);      
}

uk = 0.5 * uk;
dof = (double)(dofx * nat);
tempi = 2.0 * (uk / dof);
printf("Temperatura Inicial: %lf \n", tempi);

return; 
}

/* Esta funci�n convierte aleatorios uniformemente distribuidos en aleatorios 
gaussianamente distribuidos con distribucion normal */
double gauss(void)
{
double a1 = 3.949846138, a3 = 0.252408784, a5 = 0.076542912, a7 = 0.008355968;
double a9 = 0.029899776, sum, r, r2, rgauss;
int i;

sum = 0.0;
for(i=1; i<=12 ; i++){
   sum = sum + aleatorio_uniforme();
}
r = (sum - 6.0) / 4.0;
r2 = r * r;
rgauss = ((((a9 * r2 + a7) * r2 + a5) * r2 + a3) * r2 + a1) * r;

return rgauss;
}

/* Esta funci�n crea numeros aleatorios uniformemente distribuidos */
double aleatorio_uniforme(void)
{
int l = 5, c = 1;
int m = 2147483647;
double rrand;

seed = fmod( (seed) * l + c, m);
rrand = 1.0 * (seed) / m;

return rrand;
}

/* Leemos el archivo de configuraciones iniciales */
void leer_configuracion(int *nat, double *rho, int dofx, double rx[], 
double ry[], double rz[], double vx[], double vy[], double vz[], double *boxx, 
double *boxy, double *boxz)
{
int i, fnat;
double fboxx, fboxy, fboxz, frho;
FILE *xyz;

/* Leemos el numero de particulas, tama�o de la caja y calculamos rho */
xyz = fopen("xyz.dat", "r");
printf("Leyendo la configuracion inicial...\n");
fscanf(xyz, "%d", &fnat);
fscanf(xyz, "%lf %lf %lf", &fboxx, &fboxy, &fboxz);

if (dofx == 3){
   frho = fnat / (fboxx * fboxy * fboxz);
   printf("Densidad volumetrica reducida: %lf \n", frho);
}
else if (dofx == 2){
   frho = fnat / (fboxx * fboxy);
   printf("Densidad superficial reducida: %lf \n", frho);
}
else if (dofx == 1){
   frho = fnat / (fboxx);
   printf("Densidad lineal reducida: %lf \n", frho);
}

*nat = fnat;
*boxx = fboxx;
*boxy = fboxy;
*boxz = fboxz;
*rho = frho;

/* Leemos las posiciones y velocidades */
for(i=1 ; i<=fnat ; i++){
fscanf(xyz, "%lf %lf %lf %lf %lf %lf", &rx[i], &ry[i], &rz[i], &vx[i], &vy[i], 
&vz[i]);
}
fclose(xyz); 

return;
}

/* Guardamos un archivo con la configuracion final */
void guardar_configuracion(int nat, double boxx, double boxy, double boxz,
double rx[], double ry[], double rz[], double vx[], double vy[], double vz[])
{
int i;
FILE *xyz;

xyz = fopen("xyz.dat", "w");

fprintf(xyz, "%12d\n", nat);
fprintf(xyz, "%21.16lf %25.16lf %25.16lf\n", boxx, boxy, boxz);

for(i=1 ; i<=nat ; i++){
fprintf(xyz, "%21.16lf %25.16lf %25.16lf %25.16lf %25.16lf %25.16lf\n", rx[i], 
ry[i], rz[i], vx[i], vy[i], vz[i]);
}
fclose(xyz);

return;
}

/* Calculamos la fuerza entre particulas */
void fuerza(int nat, int dofx, double coef, double urc, double rcut, 
double *epi, double rx[], double ry[], double rz[], double boxx, double boxy, 
double boxz, double fx[], double fy[], double fz[], double expn, double expm, 
double *wxxr, double *wyyr, double *wzzr)
{
double dx, dy, dz, rij, uij, duij, sn, sm, fxij, fyij, fzij;
double epj, wxj, wyj, wzj;
int i, j;

for(i=1 ; i<= nat ; i++){
   fx[i] = 0.0;
   fy[i] = 0.0;
   fz[i] = 0.0;
}

*wxxr = 0.0;
*wyyr = 0.0;
*wzzr = 0.0;
*epi = 0.0;

if (dofx == 3){
   /* Esta etiqueta paraleliza el programa */
   #pragma omp parallel for schedule(static) shared(boxx,boxy,boxz,rx,ry,rz,fx,fy,fz) private(i,j,dx,dy,dz,rij,sn,sm,uij,duij,fxij,fyij,fzij,epj,wxj,wyj,wzj)
   /* Esta etiqueta con la gráfica */
   //#pragma acc parallel loop copyin(rx, ry, rz) copy(fx, fy, fz) reduction(+:epi, wxxr, wyyr, wzzr)
   for(i=1 ; i<=nat-1 ; i++){
      epj = 0.0;
      wxj = 0.0;
      wyj = 0.0;
      wzj = 0.0;
      for(j=i+1 ; j<=nat; j++){
         /* Calculamos la distancia entre particulas */
         dx = rx[i] - rx[j];
         dy = ry[i] - ry[j];
         dz = rz[i] - rz[j];
         /* Aplicamos la condicion de imagen minima */
         dx = dx - round(dx/boxx) * boxx;
         dy = dy - round(dy/boxy) * boxy;
         dz = dz - round(dz/boxz) * boxz;
         rij = sqrt(dx * dx + dy * dy + dz * dz);     
         /* Calculamos el potencial y su derivada para las particulas interiores
         al radio de corte */
         if (rij <= rcut){
            sn = pow((1.0 / rij), expn);
            sm = pow((1.0 / rij), expm);
            uij = coef * (sn - sm) - urc;
            duij = coef * (expm * sm - expn * sn) / rij;
            /* Acumulamos la fuerza, potencial y el virial */
            fxij = - (duij * dx / rij);
            fyij = - (duij * dy / rij);
            fzij = - (duij * dz / rij);             
            epj = epj + uij;
            wxj = wxj + fxij * dx;
            wyj = wyj + fyij * dy;
            wzj = wzj + fzij * dz;
            /* Acumulamos la fuerza para la particula i */
            fx[i] = fx[i] + fxij;
            fy[i] = fy[i] + fyij;
            fz[i] = fz[i] + fzij;
            /* Acumulamos la fuerza para la particula j Fji = -Fij */
            fx[j] = fx[j] - fxij;
            fy[j] = fy[j] - fyij;
            fz[j] = fz[j] - fzij;
         }     
      }
      /* Esta etiqueta sincroniza la escritura entre los cpu */
      #pragma omp critical(orden)
      {
      *epi = *epi + epj;
      *wxxr = *wxxr + wxj;
      *wyyr = *wyyr + wyj;
      *wzzr = *wzzr + wzj;
      }
   }
}
else if (dofx == 2){
   /* Esta etiqueta paraleliza el programa */
   #pragma omp parallel for schedule(static) shared(boxx,boxy,rx,ry,fx,fy) private(i,j,dx,dy,rij,sn,sm,uij,duij,fxij,fyij,epj,wxj,wyj)
   for(i=1 ; i<=nat-1 ; i++){
      epj = 0.0;
      wxj = 0.0;
      wyj = 0.0;
      for(j=i+1 ; j<=nat; j++){
         /* Calculamos la distancia entre particulas */
         dx = rx[i] - rx[j];
         dy = ry[i] - ry[j];
         /* Aplicamos la condicion de imagen minima */
         dx = dx - round(dx/boxx) * boxx;
         dy = dy - round(dy/boxy) * boxy;
         rij = sqrt(dx * dx + dy * dy);     
         /* Calculamos el potencial y su derivada para las particulas interiores
         al radio de corte */
         if (rij <= rcut){
            sn = pow((1.0 / rij), expn);
            sm = pow((1.0 / rij), expm);
            uij = coef * (sn - sm) - urc;
            duij = coef * (expm * sm - expn * sn) / rij;
            /* Acumulamos la fuerza, potencial y el virial */
            fxij = - (duij * dx / rij);
            fyij = - (duij * dy / rij);           
            epj = epj + uij;
            wxj = wxj + fxij * dx;
            wyj = wyj + fyij * dy;
            /* Acumulamos la fuerza para la particula i */
            fx[i] = fx[i] + fxij;
            fy[i] = fy[i] + fyij;
            /* Acumulamos la fuerza para la particula j Fji = -Fij */
            fx[j] = fx[j] - fxij;
            fy[j] = fy[j] - fyij;
         }     
      }
      /* Esta etiqueta sincroniza la escritura entre los cpu */
      #pragma omp critical(orden)
      {
      *epi = *epi + epj;
      *wxxr = *wxxr + wxj;
      *wyyr = *wyyr + wyj;
      }
   }
}
else if (dofx == 1){
   /* Esta etiqueta paraleliza el programa */
   #pragma omp parallel for schedule(static) shared(boxx,rx,fx) private(i,j,dx,rij,sn,sm,uij,duij,fxij,epj,wxj)
   for(i=1 ; i<=nat-1 ; i++){
      epj = 0.0;
      wxj = 0.0;
      for(j=i+1 ; j<=nat; j++){
         /* Calculamos la distancia entre particulas */
         dx = rx[i] - rx[j];
         /* Aplicamos la condicion de imagen minima */
         dx = dx - round(dx/boxx) * boxx;
         rij = sqrt(dx * dx);    
         /* Calculamos el potencial y su derivada para las particulas interiores
         al radio de corte */
         if (rij <= rcut){
            sn = pow((1.0 / rij), expn);
            sm = pow((1.0 / rij), expm);
            uij = coef * (sn - sm) - urc;
            duij = coef * (expm * sm - expn * sn) / rij;
            /* Acumulamos la fuerza, potencial y el virial */
            fxij = - (duij * dx / rij);            
            epj = epj + uij;
            wxj = wxj + fxij * dx;
            /* Acumulamos la fuerza para la particula i */
            fx[i] = fx[i] + fxij;
            /* Acumulamos la fuerza para la particula j Fji = -Fij */
            fx[j] = fx[j] - fxij; 
         }     
      }
      /* Esta etiqueta sincroniza la escritura entre los cpu */
      #pragma omp critical(orden)
      {
      *epi = *epi + epj;
      *wxxr = *wxxr + wxj;     
      }
   }
}

return;
}

/* Calculamos las correcciones de largo alcance */
void lrc(double coef, double expn, double expm, double pi, double rho, 
double rcut, int dofx, double *ulr, double *plr)
{
double u, p;

if (dofx == 3){
   u = 2.0 * coef * pi * rho * (1.0 / ((expn-3.0) * pow(rcut,(expn-3.0))) - 
   1.0 / ((expm-3.0) * pow(rcut,(expm-3.0))));
   p = (2.0 / 3.0) * coef * pi * pow(rho,2.0) * (expn / ((expn-3.0) * 
   pow(rcut,(expn-3.0))) - expm / ((expm-3.0) * pow(rcut,(expm-3.0))));
}
else if (dofx == 2){
   u = coef * pi * rho * (1.0 / ((expn-2.0) * pow(rcut,(expn-2.0))) - 
   1.0 / ((expm-2.0) * pow(rcut,(expm-2.0))));
   p = (1.0 / 3.0) * coef * pi * pow(rho,2.0) * (expn / ((expn-2.0) * 
   pow(rcut,(expn-2.0))) - expm / ((expm-2.0) * pow(rcut,(expm-2.0))));
}
else if (dofx == 1){
   u = 0.5 * coef * rho * (1.0 / ((expn-1.0) * pow(rcut,(expn-1.0))) - 
   1.0 / ((expm-1.0) * pow(rcut,(expm-1.0))));
   p = (1.0 / 6.0) * coef * pow(rho,2.0) * (expn / ((expn-1.0) * 
   pow(rcut,(expn-1.0))) - expm / ((expm-1.0) * pow(rcut,(expm-1.0))));
}
   
*ulr = 2.0 * u;
*plr = p;

return;
}

/* Aplicamos el termostato y barostato de Berendsen */
void berendsen(int ncolectivo, int nat, int dofx, double time, double taut, 
double taup, double temp, double tempi, double presion, double presilr,
double vx[], double vy[], double vz[], double rx[], double ry[], double rz[],
double *boxx, double *boxy, double *boxz, double *rho)
{
double factor, pcmx, pcmy, pcmz, factor2, bx, by, bz;
int i;

/* Aplicamos el termostato tanto para NVT como NPT */
factor = sqrt(1.0 + time/taut * (temp/tempi - 1.0));
//factor = sqrt(temp / tempi);  escalamiento de velocidades
for(i=1 ; i <= nat ; i++){
   vx[i] = vx[i] * factor;
   vy[i] = vy[i] * factor;
   vz[i] = vz[i] * factor;
}  
   
/* Removemos el momento del centro de masa */
pcmx = 0.0;
pcmy = 0.0;
pcmz = 0.0;
for(i=1; i<=nat ; i++){
   pcmx = pcmx + vx[i];
   pcmy = pcmy + vy[i];
   pcmz = pcmz + vz[i];
}
   
pcmx = pcmx / nat;
pcmy = pcmy / nat;
pcmz = pcmz / nat;
for(i=1; i<=nat ; i++){
   vx[i] = vx[i] - pcmx;
   vy[i] = vy[i] - pcmy;
   vz[i] = vz[i] - pcmz;     
}
   
/* Aplicamos el barostato de Berendsen en caso de usar el colectivo NPT */
if (ncolectivo == 3){
   factor2 = pow(1.0 - time/taup * (presion - presilr), 1.0/dofx);  
   bx = (*boxx) * factor2;
   by = (*boxy) * factor2;
   bz = (*boxz) * factor2;
    
   for(i=1 ; i <= nat ; i++){
      rx[i] = rx[i] * factor2;
      ry[i] = ry[i] * factor2;
      rz[i] = rz[i] * factor2;
   } 
      
   if(dofx == 3){
      *rho = nat / (bx * by * bz);
   }
   else if(dofx == 2){
      *rho = nat / (bx * by);
   }
   else if(dofx == 1){
      *rho = nat / bx;
   }
   *boxx = bx;
   *boxy = by;
   *boxz = bz;
}

return;
}

/* Creamos una pelicula con las posiciones de las particulas */
void pelicula_gro(int iconf, int nat, int dofx, double boxx, double boxy, 
double boxz, double rx[], double ry[], double rz[])
{
int i;
char marca;
double sr = 0.3405, x, y, z, lx, ly, lz;

if (iconf == 1){
movie = fopen("movie.gro", "w");
}

fprintf(movie, " CONFIGURACION %12d\n", iconf);
fprintf(movie, "%12d\n", nat);

for(i=1 ; i <= nat ; i++){
   marca = 'O';
   if (i == 1){ // Se señalan partículas para distinguir movimiento (la 1 y la 10)
   marca = 'N';   
   }
   else if (i == 10){
   marca = 'N';   
   }
   fprintf(movie, "%5dSOL %3c %7d %7.3lf %7.3lf %7.3lf\n", i, marca, i, rx[i] * 
   sr, ry[i] * sr, rz[i] * sr);
}

fprintf(movie, "%12.5lf %11.5lf %11.5lf \n", boxx * sr, boxy * sr, boxz * sr);

return;
}

/* Creamos el perfil de densidad por particula y la funcion de 
distribucion entre pares */
void distribuciones(int dofx, int nat, double rho, double deltar, double rx[], 
double ry[], double rz[], double boxx, double boxy, double boxz, double pi,
int ndist, int nparticulas[], int npares[], int nconf, int iconf)
{
int i, j, bin, binp;
double dx, dy, dz, rij, factor, denominador, rdf, r, rg;
double nperfil, areabin, nbin;

/* Contamos el numero de pares y el numero de particulas y 
hacemos el histograma */
if (dofx == 3){
   for(i=1 ; i<=nat-1 ; i++){
      for(j=i+1 ; j<=nat ; j++){
         /* Calculamos la distancia entre la particula i y la particula j */
         dx = rx[i] - rx[j];
         dy = ry[i] - ry[j];
         dz = rz[i] - rz[j];
         /* Aplicamos la condicion de imagen minima */
         dx = dx - round(dx / boxx) * boxx;
         dy = dy - round(dy / boxy) * boxy;
         dz = dz - round(dz / boxz) * boxz;
         rij = sqrt(dx * dx + dy * dy + dz * dz);
         if (rij < 0.5*boxx){
            /* bin es el numero de subconjunto para la distancia rij */
            bin = (int)(round(rij / deltar)); 
            /* Calculamos el numero de pares dentro del subconjunto ibin */
            npares[bin] = npares[bin] + 2;        
         }
      }
   }
}    
else if (dofx == 2){
   for(i=1 ; i<=nat-1 ; i++){
      for(j=i+1 ; j<=nat ; j++){
         /* Calculamos la distancia entre la particula i y la particula j */
         dx = rx[i] - rx[j];
         dy = ry[i] - ry[j];
         /* Aplicamos la condicion de imagen minima */
         dx = dx - round(dx / boxx) * boxx;
         dy = dy - round(dy / boxy) * boxy;
         rij = sqrt(dx * dx + dy * dy);
         if (rij < 0.5*boxx){
            /* ibin es el numero de subconjunto para la distancia rij */
            bin = (int)(round(rij / deltar)); 
            /* Calculamos el numero de pares dentro del subconjunto ibin */
            npares[bin] = npares[bin] + 2;        
         }
      }
   }
}    
else if (dofx == 1){
   for(i=1 ; i<=nat-1 ; i++){
      for(j=i+1 ; j<=nat ; j++){
         /* Calculamos la distancia entre la particula i y la particula j */
         dx = rx[i] - rx[j];
         /* Aplicamos la condicion de imagen minima */
         dx = dx - round(dx / boxx) * boxx;
         rij = sqrt(dx * dx);
         if (rij < 0.5*boxx){
            /* ibin es el numero de subconjunto para la distancia rij */
            bin = (int)(round(rij / deltar)); 
            /* Calculamos el numero de pares dentro del subconjunto ibin */
            npares[bin] = npares[bin] + 2;        
         }
      }
   }
}    
       
for(i=1 ; i<=nat ; i++){
   binp = (int)(round((rx[i] + 0.5 * boxx) / deltar));
   nparticulas[binp] = nparticulas[binp] + 1; 
}

   /* Escribimos las distribuciones promedio en los archivos de salida */
   if (iconf == nconf){
      gr = fopen("gr_dm.dat", "w");
      perfil = fopen("perfil_dm.dat", "w"); 
      nbin = (int)(boxx / deltar);
      
      if (dofx == 3){
         for(i=1 ; i<=nbin-1 ; i++){
            r = i * deltar;
            areabin = deltar * boxy * boxz; 
            nperfil = nparticulas[i] / (areabin * ndist);
            fprintf(perfil,"%9.5lf %12.5lf\n", r, nperfil);

            rg = r + deltar;
            factor = (4.0/3.0) * pi * rho;
            denominador = factor * nat * (pow(rg, 3) - pow(r, 3));
            rdf = npares[i] / (denominador * ndist);
            fprintf(gr,"%9.5lf %12.5lf\n", r, rdf);  
         }
      }
      else if (dofx == 2){
         for(i=1 ; i<=nbin-1 ; i++){
            r = i * deltar;
            areabin = deltar * boxy; 
            nperfil = nparticulas[i] / (areabin * ndist);
            fprintf(perfil,"%9.5lf %12.5lf\n", r, nperfil);

            rg = r + deltar;   
            factor = pi * rho;
            denominador = factor * nat * (pow(rg, 2) - pow(r, 2));
            rdf = npares[i] / (denominador * ndist);
            fprintf(gr,"%9.5lf %12.5lf\n", r, rdf);   
         }
      }
      else if (dofx == 1){
         for(i=1 ; i<=nbin-1 ; i++){
            r = i * deltar;
            areabin = deltar; 
            nperfil = nparticulas[i] / (areabin * ndist);
            fprintf(perfil,"%9.5lf %12.5lf\n", r, nperfil);
            
            rg = r + deltar;         
            factor = 2.0 * rho;
            denominador = factor * nat * deltar;
            rdf = npares[i] / (denominador * ndist);
            fprintf(gr,"%9.5lf %12.5lf\n", r, rdf);   
         }
      }
   }
      
return;         
}         

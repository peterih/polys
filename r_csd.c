/**********************************************************************
  "r_csd.c"

  Programme de conversion de donnees de la Cambridge Database en donnees 
  au format Polys 

  Written by C. Durand, INRA Nantes 1994.
  Implementation in POLYS by Soren Balling Engelsen, INRA-94.
**********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "polys.h"
#include "extern.h"


char    liste[200];
void    symetrie(int sform, int centre , int nbre);

double puiss(double X , int Y)
/**********************************************************************
**********************************************************************/
{ 
   int l;
   float val;

   val = X;
   for (l=1; l<=Y; l++)
      val /= 10.0;
   return val; 
} /* End of puiss */

/**********************************************************************
        Reads CSD format database files.
        ====================================
1 Enregistrement:  
2 Enregistrement:  
3 Enregistrement:  
4 Enregistrement:  
5 Enregistrement:  
6 Enregistrement:  
7 Enregistrement:  
**********************************************************************/
void read_CSD(char *csd_id, BOOLEAN ref)

{  FILE    *fp;
   int     sform = 1;
   int     i, k, s, l, longueur, sup, lignes, centre, nbre, liaisons, NSAT;
   char    spacegrp, record[100], chaine[11], connect[800], fname[15];
   char  *unused ;

   /* open database file */
   if ((fp = fileopen("CSD.DAT", "rt")) != NULL)
   {  printf("CSD database file CSD.DAT opened for read-only access\n\n");

      strcpy(chaine, csd_id);
      /* find the specified CSD id */
      do 
      { unused = fgets(record, 100, fp);
         if (feof(fp)) 
        { printf("\n\nCSD entry %s not found!\n\n", csd_id);
          goto fin;
        }
      } while (strncmp(chaine, record + 1, strlen(chaine)));  
      strcpy(chaine, "          "); 

      /************** Enregistrement 1 ******************************/ 
      M.id = strdup("        ");
      strncpy(M.id, record + 1, 8); 
      if (!(record[56]-48))
      { printf("\n\nCSD entry %s contains NO unit cell information!\n\n",M.id);
        goto fin;
      }

      if (!(record[58]-48))
      { printf("\n\nCSD entry %s contains NO coordinates!\n\n", M.id);
        goto fin;
      }

      strncpy(chaine, record + 44, 3);
      M.nat = atoi(chaine);
      printf("No. of atoms in CSD entry %s: %d\n\n", M.id, M.nat); 

      strncpy(chaine, record + 26, 3);
      s = atoi(chaine);

      strncpy(chaine, record + 29, 3);
      s += atoi(chaine);

      strncpy(chaine, record + 32, 3);
      s += atoi(chaine);

      strncpy(chaine, record + 35, 3);
      s += atoi(chaine);

      lignes = s / 80;
      sup = s % 80;
      if (sup > 0) lignes++;

      strncpy(chaine, record + 38, 3);
      nbre = atoi(chaine);
      centre = record[59] - 49;

      strncpy(chaine, record + 47, 3);
      NSAT = atoi(chaine);

      strncpy(chaine, record + 53, 3);
      liaisons = atoi(chaine);

      /************** Enregistrement 2 ******************************/ 
      unused = fgets(record, 100, fp);

      strncpy(chaine, record, 6);
      CELL.a = (double) atof(chaine); 
      strncpy(chaine, record + 6, 6);
      CELL.b = (double) atof(chaine);
      strncpy(chaine, record + 12, 6);
      CELL.c = (double) atof(chaine);
      strncpy(chaine, record + 18, 6);
      CELL.alpha = (double) atof(chaine);
      strncpy(chaine, record + 24, 6);
      CELL.beta = (double) atof(chaine);
      strncpy(chaine, record + 30, 6);
      CELL.gamma = (double) atof(chaine);

      strncpy(chaine, record + 36, 6);
      CELL.a     = puiss(CELL.a,     chaine[0] - 48);
      CELL.b     = puiss(CELL.b,     chaine[1] - 48);
      CELL.c     = puiss(CELL.c,     chaine[2] - 48);
      CELL.alpha = puiss(CELL.alpha, chaine[3] - 48);
      CELL.beta  = puiss(CELL.beta,  chaine[4] - 48);
      CELL.gamma = puiss(CELL.gamma, chaine[5] - 48);

      spacegrp = record[63];  
      printf("Space group symmetry: %c\n\n", spacegrp);
      switch ( record[63])
      { case 'P': spacegrp = 1;
                  break;
        case 'I': spacegrp = 2;
                  break;
        case 'R': spacegrp = 3;
                  break;
        case 'F': spacegrp = 4;
                  break;
        case 'A': spacegrp = 5;
                  break;
        case 'B': spacegrp = 6;
                  break;
        case 'C': spacegrp = 7;
                  break;
        case 'H': spacegrp = 8;
                  break;
        default:
         printf("\n\nCSD entry %s contains NO space group information!\n\n",
                               M.id);
         goto fin; 
      }

      /************** Enregistrement 3 ******************************/ 
      for (l=0; l<lignes; l++)
      unused = fgets(record, 100, fp);

      /************** Enregistrement 4 ******************************/ 
      for(i=0; i<(nbre/6 + 1); i++)
      { unused = fgets(record, 100, fp);
        strncat(liste, record, 75);
      }

      /************** Enregistrement 5 ******************************/ 
      unused = fgets(record, 100, fp);

      /************** Enregistrement 6 ******************************/ 
      i = 0; 
      while (i < M.nat)
      { unused = fgets(record, 100, fp);
        for (k=0; k<=2; k++)
        { if ((i+1) > M.nat) break; 
          strncpy(A[i].lab, record + 27*k, 5);
          switch (A[i].lab[0])
          {
            case 'C': A[i].ztyp = 600;
                      if (A[i].lab[1] == 'L') A[i].ztyp = 1700;
                      break;       
            case 'B': A[i].ztyp = 500;
                      if (A[i].lab[1] == 'R') A[i].ztyp = 3500;
                      break;
            case 'H': A[i].ztyp = 100;
                      break;
            case 'O': A[i].ztyp = 800;
                      break;
            case 'N': A[i].ztyp = 700;
                      break;
            case 'P': A[i].ztyp = 1500;
                      break;
            case 'S': A[i].ztyp = 1600;
                      if (A[i].lab[1] == 'N') A[i].ztyp = 5000;
                      break;
            case 'I': A[i].ztyp = 5300;
                      break;
            case 'F': A[i].ztyp = 900;
                      break;
            default : A[i].ztyp = 9900;
          }

          strncpy(chaine, record + 5 + 27*k, 7);
          A[i].pos.x = atof(chaine)/100000.0;

          strncpy(chaine, record + 12 + 27*k, 7);
          A[i].pos.y = atof(chaine)/100000.0;

          strncpy(chaine, record + 19 + 27*k, 7);
          A[i].pos.z = atof(chaine)/100000.0;

          A[i].res = 1;
          A[i].chg = 0.0;
          strcpy(chaine, "         ");
          i ++;
        } 
      }


      /************** Enregistrement 7 ******************************/ 
      longueur = M.nat < 100 ? 2 : 3;
      for(i=0; i<(liaisons/(82/longueur) + 1); i++)
      { unused = fgets(record, 100, fp);
        strncat(connect, record, 84 - 2*longueur);
      }
      for (k=1; k<(M.nat+1); k++)
      { s = 1;
        for (i=0; i<liaisons; i++)
        { strncpy(chaine, connect + i*longueur, longueur);
          l = atoi(chaine);
          strcpy(chaine, "         ");
          if (k == l) s++;
        }
        strncpy(chaine, connect + (k-1)*longueur, longueur);
        l = atoi(chaine);  
        if (l == 0) s--;
        A[k-1].ztyp += s;
      }

      fclose(fp);

      /* Affichage */
      M.titel = strdup("        ");
      strcpy(M.titel, M.id); 
      printf("%s\n", M.id);
      printf("%s\n", M.titel);
      printf("%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
               CELL.a, CELL.b, CELL.c, CELL.alpha, CELL.beta, CELL.gamma);
      printf("%5d\n", M.nat);
      for ( i=0; i<M.nat; i++)
         printf("%5d %3d %15.5f %15.5f %15.5f %7.4f %6d %s\n",
                  i, A[i].res, A[i].pos.x, A[i].pos.y, A[i].pos.z,
                     A[i].chg, A[i].ztyp, A[i].lab);

     frac2carte(CELL, FALSE);

     /* print symmetry neighbors */
     symetrie(sform, centre , nbre);
     printf("%d symmetry operations in the unit cell\n", 
             nbre + nbre*!centre);
     printf("Space group type no. %d\n", spacegrp);
 
     /* LOOP OUT */
     fin:;
  }

} /* End of read_csd */

void symetrie(int sform, int centre , int nbre)
/**********************************************************************
**********************************************************************/
{
   int     i, j, k, point, matrice[3][3], transl[3];
   char    symet[10];

   /* prospection par symetrie */
   for (k=0; k<nbre; k++)
   {
      /* lecture de la matrice */
      for (i=0; i<3; i++)
      {
         for (j=0; j<3; j++)
         {  matrice[i][j] = liste[i*5 + k*15 + j] - 49;
            MAT[nmat][i][j] = (double) matrice[i][j];
         }
         transl[i] = liste[i*5 + k*15 + 4] - 48;
      }

      /* traduction de la matrice */
      strcpy(symet, "");
      if (sform) 
        strcpy(symet, "(");
      for (j=0; j<3; j++)
      {  point = 0;
         if (matrice[0][j])
         {  if (matrice[0][j] == -1) 
               strcat(symet,"-x");
            else 
               strcat(symet,"x");
            point = 1;
         }
         if (matrice[1][j])
         {  if (matrice[1][j] == -1) 
               strcat(symet,"-y");
            else 
            {  if (point) 
                  strcat(symet,"+");
               strcat(symet,"y");           
            }
            point = 1;
         }
         if (matrice[2][j])
         {  if (matrice[2][j] == -1) 
               strcat(symet,"-z");
            else 
            {  if (point) 
                  strcat(symet,"+");
               strcat(symet,"z");           
            }
            point = 1;
         }     
         printf("%s", symet);
         strcpy(symet, "");

         if (transl[j])
            if (sform)   
            {  MAT[nmat][j][3] = (double) transl[j] / 12.00;
               printf("+%6.4f", transl[j] / 12.00);
            }
            else 
            {  MAT[nmat][j][3] = (double) transl[j] / 12.00;
               switch (transl[j]) 
               {  case 9:  printf("+3/4");
                           break;
                  case 8:  printf("+2/3");
                           break;
                  default: printf("+1/%d", 12 / transl[j]); 
             }
          }  
          if (j != 2) 
          {  putchar(','); 
          }
      }
      if (sform) 
         putchar(')');
      printf("\n");
      MAT[nmat][3][0] = 0.0;
      MAT[nmat][3][1] = 0.0;
      MAT[nmat][3][2] = 0.0;
      MAT[nmat][3][3] = 1.0;
      printf(" MATRIX #%d\n", nmat);
      matprint(MAT[nmat]);
      nmat++;


      /* traduction de la matrice avec centre de symetrie */

      if (!centre)
      {  if (sform) 
            strcpy(symet, "(");
         for (i=0; i<3; i++)
            for (j=0; j<3; j++)
               MAT[nmat][i][j] = (double) -matrice[i][j];
         for (j=0; j<3; j++)
         {  point = 0;
            if (matrice[0][j])
            {  if (matrice[0][j] == 1) 
                  strcat(symet,"-x");
               else 
                  strcat(symet,"x");
               point = 1;
            }
            if (matrice[1][j])
            {  if (matrice[1][j] == 1) 
                  strcat(symet,"-y");
               else 
               {  if (point) 
                     strcat(symet,"+");
                  strcat(symet,"y");           
               }
               point = 1;
            }
            if (matrice[2][j])
            {  if (matrice[2][j] == 1) 
                  strcat(symet,"-z");
               else 
               {  if (point) 
                     strcat(symet,"+");
                  strcat(symet,"z");           
               }
               point = 1;
            }
            printf("%s", symet);
            strcpy(symet, "");

            if (transl[j])
               if (sform)
               {  MAT[nmat][j][3] = (double) transl[j] / 12.00;
                  printf("+%6.4f", transl[j] / 12.00);
               }
               else
               {  MAT[nmat][j][3] = (double) transl[j] / 12.00;
                  switch (transl[j])
                  {  case 9:  printf("+3/4");
                              break;
                     case 8:  printf("+2/3");
                              break;
                     default: printf("+1/%d", 12 / transl[j]);
                  }
               }

            if (j != 2) 
            {   putchar(','); 
            }
         }
         
         if (sform)    
            putchar(')'); 
         printf("\n");

         MAT[nmat][3][0] = 0.0;
         MAT[nmat][3][1] = 0.0;
         MAT[nmat][3][2] = 0.0;
         MAT[nmat][3][3] = 1.0;
         printf(" MATRIX #%d\n", nmat);
         matprint(MAT[nmat]);
         nmat++;
      }
   }
   return;

} /* End of symmetrie */

/* End of file */

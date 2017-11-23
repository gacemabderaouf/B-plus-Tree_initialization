// Taille d'un noeud = 1024 octets [127 valeurs(entier), (ordre+1) fils (entier), nb (entier:nombre de valeurs dans le bloc)]
//  On est donc avec des arbres d'ordre 127.

/**definitions du programme d'insertion standard**/
#define ordre 127
#define initial 500
#define itiration 100
#define pas 1000000

/**definitions du programme bulk loading**/
#define Ordre 127
#define NBvalsleaf 254
#define Nbleaves 393701
#define NIL -1
#define SIZE 1024




/**************************************************************************************************/
/**************************************************************************************************/
/**************************************************************************************************/

/*bulk loading*/

/**< L'entete du fichier */
struct Tentete
{
    int teteListe;
    int teteTree;
    int offsetsController;
    int TreeLevel;
}enteteFile;
typedef struct Tentete Tentete;
/**< La structure de la feuille */
struct leaf
{
    int valuesTable[NBvalsleaf];
    int suivant;
    int nb;
}bufferleaf;
typedef struct leaf leaf;
/**< La structure d'un noeud interne */
struct internalNode
{
    int Tvals[Ordre-1];
    int Tfils[Ordre];
    int degre;
}bufferInternalNode;
typedef struct internalNode internalNode;
/**< Structure d'une liste pour la manipulation de l'arbre */
struct Tlist
{
    internalNode node;
    int ptr;
    struct Tlist *suivant;
};
typedef struct Tlist Tlist;
/**********************************************************************************************/
/**< Variables globales: */
FILE *myFile;
/**< Les  compteurs */
int nbrLectures=0;
int nbrEcritures=0;
int insertedIninternalNodes=0;
int cptInStats1=0,cptInStats2=0;
/**********************************************************************************************/

/**< Lecture d'une feuille */void lireDireLeaf(int i)
{
    nbrLectures++;
    fseek(myFile,sizeof(Tentete)+i*SIZE,SEEK_SET);
    fread(&bufferleaf,SIZE,1,myFile);
}
/**< Ecriture d'une feuille */void EcrireDireLeaf(int i)
{
    fseek(myFile,sizeof(Tentete)+i*SIZE,SEEK_SET);
    fwrite(&bufferleaf,SIZE,1,myFile);
}
/**< ecriture d'un noeud interne */void ecrireDireInternalNode(internalNode *buff,int i)
{
    nbrEcritures++;
    fseek(myFile,sizeof(Tentete)+i*SIZE,SEEK_SET);
    fwrite(buff,SIZE,1,myFile);
}
/**< lecture d'un noeud interne */void lireDireInternalNode(internalNode *buff,int i)
{
    fseek(myFile,sizeof(Tentete)+i*SIZE,SEEK_SET);
    fread(buff,SIZE,1,myFile);
}
/**< Affectation a l'entete */void affEntete(int i,int val)
{
    if (i==1) enteteFile.teteListe=val;
    else if (i==2) enteteFile.teteTree=val;
    else if (i==3) enteteFile.offsetsController=val;
    else enteteFile.TreeLevel=val;
}
/**< Extraire une information de l'entete */int Entete (int i)
{
    if (i==1) return enteteFile.teteListe;
    else if (i==2) return enteteFile.teteTree;
    else if (i==3) return enteteFile.offsetsController;
    else return enteteFile.TreeLevel;
}
/**< alloeur un noeud interne */int allocBlocInternalNode()
{
    affEntete(3,Entete(3)+1);
    return Entete(3);
}
/**< Ecrire l'entete du fichier */void ecrireEntete()
{
    fseek(myFile,0,SEEK_SET);
    fwrite(&enteteFile,sizeof(Tentete),1,myFile);
}
/**< Lecture directe de l'entete */void lireEntete(FILE *pf)
{
    fseek(pf,0,SEEK_SET);
    fread(&enteteFile,sizeof(Tentete),1,pf);
}
/**< Creeer un fichier LOF constituant la base de données */void createTestFile()
{
    myFile=fopen("treeFile.bin","wb+");
    int i,j;
    for (i=0;i<Nbleaves;i++)
    {
        for (j=0;j<NBvalsleaf;j++) bufferleaf.valuesTable[j]=i*NBvalsleaf+j;
        if (i!=Nbleaves-1) bufferleaf.suivant=(i+1);
        else bufferleaf.suivant=NIL;
        bufferleaf.nb=254;
        EcrireDireLeaf(i);
    }
    affEntete(1,0);
    affEntete(2,NIL);
    affEntete(3,Nbleaves);
    affEntete(4,0);
}
/**< Eclater un noeud interne */void eclatement(int ptr,int *fd,int *val,Tlist *cursor)
{
    int i,j=-1,tempo=*val;
    (*val)=(cursor)->node.Tvals[(Ordre-1)/2];
    (cursor)->node.degre=(Ordre/2)+1;
    ecrireDireInternalNode(&(cursor->node),ptr);
    insertedIninternalNodes+=(cursor->node.degre)-1;
    for (i=Ordre-2;i>(Ordre-1)/2;i--)
    {
        j++;
        (cursor)->node.Tvals[j]=(cursor)->node.Tvals[i];
        (cursor)->node.Tfils[j]=(cursor)->node.Tfils[i];
    }
    (cursor)->node.Tfils[j+1]=(cursor)->node.Tfils[Ordre-1];
    (cursor)->node.Tfils[j+2]=*fd;
    (cursor)->node.Tvals[j+1]=tempo;
    *fd=allocBlocInternalNode();
    (cursor)->ptr=*fd;
}

/**********************************************************************************/
/**********************************************************************************/

/*insertion standard*/

typedef int TVal;	// type de la valeur
#include <stdio.h>
#include <stdlib.h>
#include<time.h>

/*******variables globales*******/
int nbLect;
int nbEcr;
int indexVal,indexBloq;
float moyCharge,moyChargeIdex,moyChargeFich;


/*******Definition des structures*******/


///Definition de la structure entête

typedef struct blocE {
	int racine;
	int nbnoeud;
	int nbvaleurs;
	int prof;
} Entet;


///Definition de la structure TNoeud (un bloc du fichier B+Arbre)

typedef struct noeud {
   int nb;
   TVal cle[ordre];
   int fils[(ordre+1)];
   int inutil[(1024-8*(ordre+1))/4];
} TNoeud;

typedef struct feuille {
    int nb;
    int suiv;
    int prec;
    TVal cle[253];
} TFeuille;

///Definition de la structure TBloq (du fichier TOF)

typedef struct bloq {
    int nb;
    TVal cle[255];
}TBloq;


///Definition de la structure TPile utilise dans le parcours du B+Arbre
typedef struct pile TPile;
struct pile {
   int nmBloc;
   TNoeud bloc;
   TPile* suiv;
} ;

/*///liste pour les statistiques, on la notera file

typedef struct list Tlist;
struct list {
    int bloq;
    Tlist* suiv;
};*///??????????????

///Definition de structure BTree un Fichier en B+Arbre

typedef struct btree {	 // Un Fichier organise en B+Arbre
   FILE *arb;		         // Le fichier contenant les noeuds
   Entet entete;        // les caracteristiques du B+Arbre :(racine,nb de blocs,nb de valeurs)
   TPile *pile;	         // Pile de parcours (garde la trace des blocs visites)
} BTree;
#define TailleEntete sizeof(Entet)
#define TailleBloc sizeof(TNoeud)
/*
*
*
*/

/*******Entetes des fonctions*******/
void ins_noeud( TVal v, int j, int fd, TNoeud *n );
int rech_noeud( TVal v, int *i, const TNoeud *n );
void Affich_noeud(const TNoeud *n);
void Init_pile(BTree *bt);
void Empiler(BTree *bt, int i, TNoeud n);
void Depiler(BTree *bt, int *i, TNoeud *n);
BTree *Ouvrir_btree(char *name);
void Fermer_btree(BTree *bt);
void LireBloc_btree(BTree *bt, TNoeud *buffer, int i);
void EcrireBloc_btree(BTree *bt, const TNoeud *buffer, int i);
int AllocBloc_btree(BTree *bt);
void affich_entete_btree(BTree *bt);
int Rech_btree(BTree *bt, TVal c, int *i, int *j, TFeuille *buff);
void Ins_btree(BTree *bt, TVal v);
void Stat(BTree *bt);
void EcrireBloc_bloq(FILE *f, const TBloq *buffer, int i);
void LireBloc_bloq(FILE *f, TBloq *buffer, int i);
int rech_feuille( TVal v, int *i, const TFeuille *n );
/*
*
*
*/

/*******Corps des fonctions*******/


///Insertion ordonnee de (v:cle, j:position, fd:fils droit) dans le noeud

void ins_noeud( TVal v, int j, int fd, TNoeud *n )
{
   if (n->nb < ordre) {
      int jj ;

  /* decalage a droite des cles et des fils */
      for (jj=n->nb; jj>j; jj--) {
    n->cle[jj]    = n->cle[jj-1];
    n->fils[jj+1] = n->fils[jj];
      }

  // insertion de v dans cle et fd dans fils
      n->cle[j]    = v;
      n->fils[j+1] = fd;
      n->nb++;
   }
}

///Insertion ordonnee de (v:cle, j:position ) dans la feuille

void ins_feuille( TVal v, int j, TFeuille *n )
{
   if (n->nb < 253) {
      int jj ;

  /* decalage a droite des cles  */
      for (jj=n->nb; jj>j; jj--) {
    n->cle[jj]    = n->cle[jj-1];
      }

  // insertion de v dans cle
      n->cle[j]    = v;
      n->nb++;
   }
}

// recherche de v et retourne sa pos dans i

int rech_noeud( TVal v, int *i, const TNoeud *n )
{
   int bi = 0, bs = n->nb-1, trouv = 0;
   while ( !trouv && bi <= bs )  {
  *i = (bi +bs) / 2;
  if ( v < n->cle[*i] )
    bs = *i - 1;
  else
     if ( v > n->cle[*i] )
    bi = *i + 1;
     else
    trouv = 1;
   }
   if ( !trouv )
  *i = bi;
  else
    (*i)++; // chemin à prendre fils[i]

   return trouv;

}

// recherche de v dans une feuille et retourne sa pos dans i

int rech_feuille( TVal v, int *i, const TFeuille *n )
{
   int bi = 0, bs = n->nb-1, trouv = 0;
   while ( !trouv && bi <= bs )  {
  *i = (bi +bs) / 2;
  if ( v < n->cle[*i] )
    bs = *i - 1;
  else
     if ( v > n->cle[*i] )
    bi = *i + 1;
     else
    trouv = 1;
   }
   if ( !trouv )
  *i = bi;


   return trouv;

}



// Initialiser la pile p

void Init_pile(BTree *bt)
{
    TPile *q;

    while((bt->pile)!=NULL){

      q=(bt->pile)->suiv;
      free(bt->pile);
      bt->pile=q;
  }

}


// Empiler le num de bloc (i) et son contenu (n)

void Empiler(BTree *bt, int i, TNoeud n)
{
    TPile *q;

    q=(TPile*) malloc(sizeof(TPile));
    q->nmBloc=i;
    q->bloc=n;
    q->suiv=bt->pile;
    bt->pile=q;
}


// Depiler le num de bloc dans i et son contenu dans n

void Depiler(BTree *bt, int *i, TNoeud *n)
{
    TPile* Q;

   *i = (bt->pile)->nmBloc;
   *n = (bt->pile)->bloc;
    Q=bt->pile;
    (bt->pile)=(bt->pile)->suiv;
    free(Q);
}

// Crée ou ouvre un fichier B-Arbre

BTree *Ouvrir_btree(char *name)
{
   BTree *bt = malloc( sizeof(BTree) );


      bt->arb = fopen(name,"wb+");
      bt->entete.racine = -1;     // Racine = -1
      bt->entete.nbnoeud = 0;     // Nb de noeud = 0
      bt->entete.prof = 0;      // Hauteur = 0
      bt->entete.nbvaleurs = 0;// Nb de valeurs initiales
      bt->pile=NULL;
      fwrite( &(bt->entete),sizeof(Entet),1, bt->arb );
   return bt;
}


// Ferme un fichier B-Arbre

void Fermer_btree(BTree *bt)
{
   fseek(bt->arb, 0, 0);        // Sauver l'entete
   fwrite( &(bt->entete), sizeof(Entet), 1, bt->arb ); //  et
   fclose(bt->arb);         // Fermer le fichier
   free(bt);
}


// lire le bloc i dans la var buffer

void LireBloc_btree(BTree *bt, TNoeud *buffer, int i)
{
   fseek(bt->arb, TailleEntete + i*TailleBloc, 0);
   fread(buffer, TailleBloc, 1, bt->arb);
   nbLect++;
}

// lire le bloc i dans la var buffer

void LireBloc_feuille(BTree *bt, TFeuille *buffer, int i)
{
   fseek(bt->arb, TailleEntete + i*TailleBloc, 0);
   fread(buffer, TailleBloc, 1, bt->arb);
   nbLect++;
}

// lire le bloc i du fichier TOF dans la var buffer

void LireBloc_bloq(FILE *f, TBloq *buffer, int i)
{
   fseek(f,  i*sizeof(TBloq), 0);
   fread(buffer, sizeof(TBloq), 1, f);
   nbLect++;
}



// ecrire la var buffer dans le bloc i

void EcrireBloc_btree(BTree *bt, const TNoeud *buffer, int i)
{
   fseek(bt->arb, TailleEntete + i*TailleBloc, 0);
   fwrite(buffer, TailleBloc, 1, bt->arb);
   nbEcr++;
}

// ecrire la var buffer dans le bloc i

void EcrireBloc_feuille(BTree *bt, const TFeuille *buffer, int i)
{
   fseek(bt->arb, TailleEntete + i*TailleBloc, 0);
   fwrite(buffer, TailleBloc, 1, bt->arb);
   nbEcr++;
}

// ecrire la var buffer dans le fichier TOF

void EcrireBloc_bloq(FILE *f, const TBloq *buffer, int i)
{
   fseek(f, i*sizeof(TBloq), 0);
   fwrite(buffer, sizeof(TBloq), 1, f);
}

// Allocation d'un nouveau bloc au fichier :
// on rajoute un bloc à la fin du fichier (s'il n'y a pas de blocs déjà libéré)


int AllocBloc_btree(BTree *bt)
{
   return bt->entete.nbnoeud+1;
}




// Affichage des caracteristiques du fichier (l'entete)

void affich_entete_btree(BTree *bt)
{

   printf("Racine:\t\t@%d\n", bt->entete.racine);
   printf("Nb bloc:\t%d (%d dans l'index)\n", bt->entete.nbnoeud,indexBloq);
   printf("Hauteur:\t%d\n", bt->entete.prof);
   printf("Nb valeurs:\t%d\n", bt->entete.nbvaleurs);
}

// RECHERCHE de c dans le B-Arbre, retourne en plus du booleen,
// le num de bloc i, le depl j et le contenu du dernier bloc dans buf
// de plus les blocs de la branche visitee sont sauvegardés dans la pile

int Rech_btree(BTree *bt, TVal c, int *i, int *j, TFeuille *buff)
{
    TNoeud buffer,*buf=&buffer;

    *i = bt->entete.racine;
   *j = -1;
   Init_pile (bt );
   int trouv = 0;
   if(*i!=-1){
    int k=bt->entete.prof,h;

    for(h=2;h<=k;h++){
     LireBloc_btree(bt, buf, *i);
     trouv = rech_noeud(c, j, buf);
     Empiler( bt, *i, *buf );
     *i = buf->fils[*j];
    }
    LireBloc_feuille(bt,buff,*i);
    trouv = rech_feuille(c, j, buff);
 }

   return trouv;
}
#define m printf("\n");
// INSERTTION d'une cle v dans le B-Arbre

void Ins_btree(BTree *bt, TVal v)
{
   int i,j;
   TNoeud buf;
   buf.nb = 0;
   TFeuille buff;
   buff.nb=0;
   int k,inutil;

   if ( !Rech_btree(bt,v,&i,&j,&buff) )
    {

         if (i == -1) // si arbre vide ...
         {

                bt->entete.racine = AllocBloc_btree(bt);  // nouvelle racine
                ins_feuille(v, 0, &buff);
                buff.prec=0; buff.suiv=0;
                EcrireBloc_feuille(bt, &buff, bt->entete.racine);
                bt->entete.nbvaleurs++;
                bt->entete.nbnoeud++;
                bt->entete.prof++;

         }
         else
         {   // si arbre non vide ...

                if(buff.nb<253)// il y a de la place dans le bloque
                {
                    ins_feuille(v,j, &buff);
                    EcrireBloc_feuille(bt,&buff,i);
                    bt->entete.nbvaleurs++;

                }
                else
                { // un eclatement s'impose
                    TVal c[(ordre+1)]; // tableau temporaire de cles
                    int f[(ordre+2)];  // tableau temporaire de fils
                    int fd,suiv;
                    TVal cle[254]; // tableau temporaire de cles

                    k=0;
                    while (k<253 && buff.cle[k] < v)
                    {
                        cle[k] = buff.cle[k];
                        k++;
                    }
                    cle[k] = v;
                    k++;
                    while (k < 254)
                    {
                       cle[k] = buff.cle[k-1];
                       k++;
                    }
                      // 1ere moitie de la sequence dans i ...
                    for (k=0; k < (254 / 2); k++)
                    {
                        buff.cle[k] = cle[k];
                    }
                    buff.nb = 254/2;

                      // 2e moitie dans un nouveau bloc ...

                    fd = AllocBloc_btree(bt); // un nouveau noeud pour l'eclatement
                    suiv=buff.suiv;
                    buff.suiv=fd*(-1);
                    EcrireBloc_feuille(bt, &buff, i); //ecreture de la 1ere moitie

                    while (k < 254)
                    {
                        buf.cle[k-254/2] =cle[k];
                        k++;
                    }
                    buff.prec=i;
                    buff.suiv=suiv;
                    buff.nb = 254/2;
                    EcrireBloc_feuille(bt, &buff, fd);
                    if(suiv!=0) {
                            LireBloc_feuille(bt,&buff,suiv*(-1));
                            buff.prec=fd;
                            EcrireBloc_feuille(bt,&buff,suiv*(-1));
                    }
                    bt->entete.nbnoeud++;
                    bt->entete.nbvaleurs++;
                          // insertion du sep dans l'index (prepa de v,buf,i,j)


                    if((bt->entete.prof)>1)
                    {

                        Depiler(bt,&i,&buf);
                        inutil=rech_noeud(v,&j,&buf);
                    }
                    else
                    {
                       i=-1;
                       buf.nb=0;
                    }
                    v=cle[254/2];
                    if (i == -1)
                    {  // si index vide ...
                        k=bt->entete.racine;
                        bt->entete.racine = AllocBloc_btree(bt);  // nouvelle racine
                        ins_noeud(v,0, fd, &buf);
                        buf.fils[0] =k;
                        EcrireBloc_btree(bt, &buf, bt->entete.racine);
                        bt->entete.nbvaleurs++;
                        bt->entete.nbnoeud++;
                        bt->entete.prof++;
                    }
                    else
                    {   // si index non vide ...
                        int continu = 1;
                        while (continu)   // Boucle Principale...
                           // Si le bloc n'est pas plein ...
                        if (buf.nb < ordre)
                         {
                              ins_noeud(v,j, fd, &buf); // inserer dans le bloc
                              EcrireBloc_btree(bt,&buf,i);// et
                              continu = 0;  // sortir de la boucle
                         }
                        else
                         {
                           // Sinon : le bloc est plein ...
                            // Eclatement du noeud courant


                            // Construction de la seq ordonnee dans c et f
                              k = 0;
                              f[0] = buf.fils[0];
                              while (k<ordre && buf.cle[k] < v)
                              {
                                  c[k] = buf.cle[k];
                                  f[k+1] = buf.fils[k+1];
                                  k++;
                              }
                              c[k] = v;
                              f[k+1] = fd;
                              k++;
                              while (k < (ordre+1))
                               {
                                  c[k] = buf.cle[k-1];
                                  f[k+1] = buf.fils[k];
                                  k++;
                                }

                                // 1ere moitie de la sequence dans i ...
                                buf.fils[0] = f[0];
                                for (k=0; k < (ordre / 2); k++)
                                {
                                    buf.cle[k] = c[k];
                                    buf.fils[k+1] =f[k+1];
                                }
                                buf.nb = k++;
                                EcrireBloc_btree(bt, &buf, i);

                                // 2e moitie dans un nouveau bloc ...
                                fd = AllocBloc_btree(bt); // un nouveau noeud pour l'eclatement
                                buf.fils[0] = f[k];
                                while (k < (ordre+1))
                                {
                                    buf.cle[k-(ordre/2)-1] =c[k];
                                    buf.fils[k-(ordre/2)] = f[k+1];
                                    k++;
                                }
                                buf.nb = k - (ordre/2) - 1;
                                EcrireBloc_btree(bt, &buf, fd);

                                    // le separateur dans le pere ...
                                v = c[ordre / 2];
                                if (bt->pile!=NULL)
                                {
                                    // si le pere existe,le depiler
                                    k=i;
                                    Depiler( bt, &i,&buf );
                                    j=0;
                                    while(buf.fils[j]!=k )j++;
                                }
                                else
                                {
                                     // si le pere n'existe pas : donc nouvelle racine
                                    bt->entete.prof++;
                                    bt->entete.nbnoeud++;
                                    bt->entete.racine = AllocBloc_btree(bt);
                                    buf.cle[0] = v;
                                    buf.fils[0] = i; // le 1er fils est l'ancienne racine
                                    buf.fils[1] = fd;// le 2e = noeud de l'eclatement
                                    buf.nb = 1;
                                    EcrireBloc_btree(bt, &buf, bt->entete.racine);
                                    continu = 0;  // sortir de la boucle
                                }
                                bt->entete.nbnoeud++;

                          } //  Fin sinon (le bloc est plein...) et Fin boucle (while (continu)...)
                          bt->entete.nbvaleurs++;
                    }
               }


         } // Fin : si arbre non vide.
    }


} // Fin Insertion

/*void Stat(BTree *bt)
{
    TNoeud buf;
    TFeuille buff;
    int k,i=bt->entete.racine,racine=1;
    Tlist *tete=NULL,*prec=NULL,*visit,*tempo;


    for(k=2;k<=bt->entete.prof;k++)
    {
        if(racine){
                printf("kkkkkkkkkkkkkkk");
            LireBloc_btree(bt,&buf,i);
            int h;
            for(h=0;h<=buf.nb;h++)
            {
                prec=tete;
                tete=(Tlist*)malloc(sizeof(Tlist));
                printf("tete liste %p",tete);m//??????
                tete->bloq=buf.fils[h];
                tete->suiv=prec;
            }

            moyCharge=moyCharge+((float)buf.nb/ordre)*100;
            indexVal+=buf.nb;
            indexBloq++;
            moyChargeIdex=moyChargeIdex+((float)buf.nb/ordre)*100;
            racine=0;
            visit=tete;
            tete=NULL;
            printf("fin stat racine");m//??????????????
        }

        while(visit!=NULL)
        {

            i=visit->bloq;
            LireBloc_btree(bt,&buf,i);
             int h;
            for(h=0;h<=buf.nb;h++)
            {
                prec=tete;
                tete=(Tlist*)malloc(sizeof(Tlist));
                tete->bloq=buf.fils[h];
                tete->suiv=prec;
            }
            tempo=visit;
            visit=visit->suiv;
            free(tempo);
            moyCharge=moyCharge+((float)buf.nb/ordre)*100;
            indexVal+=buf.nb;
            indexBloq++;
            moyChargeIdex=moyChargeIdex+((float)buf.nb/ordre)*100;
        }
        visit=tete;
        tete=NULL;

    }
    if (bt->entete.prof>1){
        tete=visit;
        i=tete->bloq;
        while(tete!=NULL)
        {
            tempo=tete;
            tete=tete->suiv;
            free(tempo);
        }

        while(i!=-1)
        {
            LireBloc_feuille(bt,&buff,i);
            k=i;
            i=buff.suiv;
        }
    }
    else k=i;
    while(k!=-1)
    {
        LireBloc_feuille(bt,&buff,k);
        moyCharge=moyCharge+((float)buff.nb/253)*100;
        moyChargeFich=moyChargeFich+((float)buff.nb/253)*100;
        k=buff.prec;
    }

}*/
void Stat(BTree *bt)
{
    TNoeud buf;
    int k;


    for(k=1;k<=bt->entete.nbnoeud;k++)
    {
        LireBloc_btree(bt,&buf,k);
        if(buf.cle[0]>0)
        {
            if(k!=bt->entete.racine){
                moyCharge=moyCharge+((float)buf.nb/ordre)*100;
                indexVal+=buf.nb;
                indexBloq++;
                moyChargeIdex=moyChargeIdex+((float)buf.nb/ordre)*100;
            }
            else {
                    indexVal+=buf.nb;
                    indexBloq++;
            }
        }
        else
        {
            moyCharge=moyCharge+((float)buf.nb/253)*100;
            moyChargeFich=moyChargeFich+((float)buf.nb/253)*100;
        }

    }

}
/***************************************************************************************/

int main()
{
    int Root=allocBlocInternalNode(),fd=-1,fg=-1;/**< Allouer une recine au depart */
    int testFileHead;
    int stop,val,createNewBloc=0;
    clock_t end;
    double time_spent=0;
    Tlist *tempoListHead=malloc(sizeof(Tlist)),*cursor=NULL,*q;/**< On initilaize une liste pour manipuler une liste pour garder la branche toute a droite */
    Tlist *usedInStats=NULL;
    int nbinsertedVals=0,repere=10;
    FILE *statisticsFile=fopen("stats.txt","w+");/**< creation d'un fichier contenant des statistics */
    createTestFile();/**< Lancement de la creation d'un fichier LOF */
    printf("La creation du fichier original des enregistrements a ete eu lieu\n");
    clock_t begin = clock();/**< ORIGINE DES TEMPS */
    /**
     * Pour manipuler la derniere branche on utilise une liste, (servire presque comme une pile) a chaque fois qu'on ajoute
     * une racine, on ajout un maillon a la fin de la liste, on savegarde aussi la tete de la liste, pour y acceder a chaque ittération
     * a la fin, on ecrit toute la branche a droite
     */
    tempoListHead->node.degre=1;
    tempoListHead->node.Tfils[0]=0;
    tempoListHead->ptr=Root;
    tempoListHead->suivant=NULL;
    affEntete(2,Root);
    affEntete(4,1);
    lireDireLeaf(0);
    testFileHead=bufferleaf.suivant;
    while (testFileHead!=-1)/**< On sort de la boucle que si on join toutes les feuilles à l'arbre */
    {
        cursor=tempoListHead;/**< cursor est un pointeur pour parcourire la liste */
        stop=0;
        lireDireLeaf(testFileHead);
        fd=testFileHead;
        val=bufferleaf.valuesTable[0];
        while (!stop)
        {
            if (!createNewBloc)
            {

                if ((cursor)->node.degre<Ordre)
                {
                    cursor->node.degre++;
                    cursor->node.Tfils[cursor->node.degre-1]=fd;
                    cursor->node.Tvals[cursor->node.degre-2]=val;
                    stop=1;
                }
                else
                {
                    fg=cursor->ptr;
                    eclatement(cursor->ptr,&fd,&val,cursor);
                    if (cursor->suivant!=NULL) cursor=cursor->suivant;
                    else createNewBloc=1;
                }
            }
            else
            {
                Root=allocBlocInternalNode();
                affEntete(2,Root);
                affEntete(4,Entete(4)+1);
                bufferInternalNode.degre=2;
                bufferInternalNode.Tfils[0]=fg;
                bufferInternalNode.Tfils[1]=fd;
                bufferInternalNode.Tvals[0]=val;
                q=malloc(sizeof(Tlist));
                q->node=bufferInternalNode;
                q->ptr=Root;
                q->suivant=NULL;
                cursor->suivant=q;
                stop=1;
                createNewBloc=0;
            }
        }
        testFileHead=bufferleaf.suivant;
        nbinsertedVals+=254;
        if (nbinsertedVals>repere){
            end=clock();
            time_spent=((double)(end-begin)/CLOCKS_PER_SEC)*1000;
            cptInStats1=0;cptInStats2=0;usedInStats=tempoListHead;
            while (usedInStats!=NULL)
            {
                /**< Faire des states sur la pile qui n'est pas encore ecrite */
                cptInStats1++;
                cptInStats2+=tempoListHead->node.degre-1;
                usedInStats=usedInStats->suivant;
            }
            fprintf(statisticsFile,"%i;%f;%i;%i;%f;%f\n",nbinsertedVals,time_spent,nbrEcritures+cptInStats1,Entete(4)+1,((float)(insertedIninternalNodes+cptInStats2)/(nbrEcritures+cptInStats1)/126)*100,((float)(insertedIninternalNodes+cptInStats2+nbinsertedVals)/(nbrEcritures+cptInStats1+nbinsertedVals/127)/127)*100);
            repere+=10000;
        }
    }
    ecrireEntete();
    while (tempoListHead!=NULL)/**< Dans cette boucle on ecrit toute la liste qui contient la derniere branche */
    {
        ecrireDireInternalNode(&(tempoListHead->node),tempoListHead->ptr);
        tempoListHead=tempoListHead->suivant;
    }
    fclose(myFile);
    fclose(statisticsFile);

    /**********************************************
    ***********************************************/
    int k=0,vall=rand()%50+1,i=0,nbVal=initial,nbiter;
    TBloq buf;
    FILE *tof,*stat;
    char name[20]="fich_btree";
    clock_t ins_begin,ins_end;
    double ins_spent;

    srand(time(NULL));
    stat=fopen("fich_stat.txt","a+");
    fprintf(stat,"nbrVal--hauteur--temps--nbrLect--nbrEcr--nbrNoeudFich--nbrNoeudIndex--nbValEn+--moyCharg--moyChargFich--moyChargindex \n");

    for(nbiter=0;nbiter<itiration;nbiter++){
            i=0;
            //printf("construction du fichier TOF de %d valeurs ...",k);m
            tof=fopen("fich_tof","ab+");

            while(k<=nbVal)
            {
                int j,nb=0;

                buf.cle[0]=vall;k++;nb++;
                for(j=1;j<255;j++)
                {
                    k++;
                    if(k>nbVal) break;
                    buf.cle[j]=buf.cle[j-1]+rand()%50+1;
                    nb++;


                }
                buf.nb=nb;
                EcrireBloc_bloq(tof,&buf,i);
                i++;
                vall=buf.cle[254]+rand()%50+1;
            }
            //printf("fin de construction du fichier TOF de %d valeurs ...",k-1);m m


            BTree* bt=Ouvrir_btree(name);
            k=1;i=0;

            nbLect=0; nbEcr=0; moyCharge=0; indexVal=0;
            indexBloq=0; moyChargeFich=0; moyChargeIdex=0;
            ins_begin=clock();
            while(k<=nbVal)
            {
                int j;
                LireBloc_bloq(tof,&buf,i);
                for(j=0;j<buf.nb;j++)
                {
                    Ins_btree(bt,buf.cle[j]);
                    k++;
                }
                i++;
            }
            ins_end=clock();
            ins_spent=(double)(ins_end-ins_begin)/CLOCKS_PER_SEC;

            //k=ordre;



            //printf("il sagit d'un arbre d'ordre %d",k);m
            //printf("temps necessaire d'insertion: %f s",ins_spent);m
            //printf("le nombre d'acces disque est %d : %d en lecture %d en ecriture ",nbLect+nbEcr,nbLect,nbEcr);m
            k=nbVal;
            fprintf(stat,"%d;%d;%f;%d;%d;",k,bt->entete.prof,ins_spent,nbLect,nbEcr);
            Stat(bt);

            //affich_entete_btree(bt);
            //printf("taille du fichier d'origine %.2f Ko --- taille du B+arbre %.2f Ko \n il y a %d plus de valeurs dans le B+arbre que dans le fichier d'origine",(float)(nbVal/255+1)/1000,(float)bt->entete.nbnoeud/1000,indexVal);m
            //if (indexBloq>1) printf("moy de charge %.2f %% --> dans l'index %.2f %% dans les bloques du fichier %.2f %%",moyCharge/(bt->entete.nbnoeud-1),moyChargeIdex/(indexBloq-1),moyChargeFich/(bt->entete.nbnoeud-indexBloq));
            //else printf("moy de charge %.2f %% ",moyCharge/(bt->entete.nbnoeud-1));

            if (indexBloq>1) fprintf(stat,"%d;%d;%d;%f;%f;%f;\n",bt->entete.nbnoeud-indexBloq, indexBloq,indexVal, moyCharge/(bt->entete.nbnoeud-1),moyChargeFich/(bt->entete.nbnoeud-indexBloq),moyChargeIdex/(indexBloq-1));
            else             fprintf(stat,"%d;%d;%d;%f;%f;\n",bt->entete.nbnoeud-indexBloq, indexBloq,indexVal, moyCharge/(bt->entete.nbnoeud-1),moyChargeFich/(bt->entete.nbnoeud-indexBloq));
            Fermer_btree(bt);;
            fclose(tof);
            nbVal=nbVal+pas;
            printf("%d-",nbiter);
    }

        fclose(stat);



    return 0;
}

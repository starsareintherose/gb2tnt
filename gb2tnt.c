/*
Copyright (C) 2008 Pablo A. Goloboff
Copyright (C) 2023 Guoyi Zhang;

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc., 59
Temple Place - Suite 330, Boston, MA 02111-1307, USA.

Email: pablogolo@csnat.edu.ar
Mail: Pablo A. Goloboff, INSUE, Instituto Miguel Lillo, Miguel Lillo 205, 4000 S.M. de Tucuman, Argentina. 
Email: GuoyiZhang@malacology.net
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
int verb = 0 ;
int use_translation = 0 ; 
int dargc , laquiero ;
int showedskipped = 0 ; 
char ** dargv ;
int fileargs ;   // args 1 to fileargs-1 are file names!
FILE * inpf , * opsf , * curinput , * opsf , * notsfile ;
int laschar ;
#define NUCLEAR 1
#define MITOCH  2
#define PLASTID 3
#define CHLORO  4
#define TRUE    1
#define FALSE   0
#define MAXUSERACC 50000 
int genometype = 0 ;
int use_string_matching = 0 ;
double string_similarity ; 

int stringis ( char * , char * ) ;  //  are both strings the same, case aside ??? 
int rd ( void ) ;                  //  return next char. from input file; save it in laschar
void outer ( int , char * txt ) ;  // if first arg. is TRUE, then output txt and exit 
void gonln ( void ) ;              //  read all the way to the ENTER
void gotostring ( char * , int ) ;  // find string of first arg, with margin of sec arg characters
void rdto ( int , char * ) ;        //  saves to string every byte it finds before first arg
void stornexline ( void ) ;         // puts in stringsp everything to the ENTER
void * mymalloc ( unsigned long int ) ;  // allocs mem or exits
void process (void);
void save_to_not ( void );
void setopts ( void );
void makelower ( char * txt );
void parsit ( void );
void parse_translation ( void );
void spew_name ( void );
void output_translation ( void );
void effect_complementation ( void );

#define MAXNUMCHUNKS 200
typedef struct { int from ; int to ; } Chunktyp ;
Chunktyp chunk[MAXNUMCHUNKS] ;
int numchunks = 0 ;

#define MAXSEQLENGTH 25000 
char bytestring[MAXSEQLENGTH] ;
char complement_string[MAXSEQLENGTH] ;

unsigned long int accepted = 0 , rejected = 0 ; 
char stringsp[160] , headerline[160];
int prodnumber = 0 , genenumber = 0 ;
#define MAX_USER_DEFINITIONS 60 
char * prodname[MAX_USER_DEFINITIONS] , * genename[MAX_USER_DEFINITIONS ] , 
     accnumber [ 100 ] ;
char taxname [ 100 ] ;
char taxonomy [ 1000 ] ;
char ** useraccname ;
int useraccnumber = 0 ;

unsigned long int dafsize , bytesread = 0 ;

#ifdef LINUX

int getch ( void )
{
    return ( getc ( stdin ) ) ;
}
#endif
    
void * mymalloc ( unsigned long int size )
{
    void * pt = malloc ( size ) ;
    outer ( pt == NULL , "Not enough memory!" ) ;
    return pt ; 
}

void dildit ( void )
{
    static int prv = 0 ; 
    double fract ;
    unsigned long int ifract ; 
    ++ bytesread ; 
    fract = ( ( double ) bytesread / dafsize ) * 100 ;
    ifract = fract ; 
    if ( ifract == prv ) return ;
    prv = ifract ; 
    fprintf ( stderr , "\rParsing ... %lu%%" , ifract ) ; 
}

void undild ( void )
{
    fprintf ( stderr , "\r                             \r" ) ; 
}

int rd ( void )
{
    int i ;
    if ( feof ( curinput ) ) {
        fclose ( curinput ) ;
        if ( inpf == curinput ) {
             undild () ; 
             if ( !useraccnumber )
                  fprintf ( stderr , "\nAccepted %lu accessions, rejected %lu" , accepted , rejected ) ; 
             fprintf ( stderr , "\nFinished parsing %s (input)\n" , dargv[1] ) ; }
        else fprintf ( stderr , "\nFinished parsing %s (options)\n" , dargv[2] ) ;
        if ( inpf == curinput ) {
            getc ( stdin ) ;
            exit ( 1 ) ; }}
    if ( curinput == inpf ) dildit () ; 
    i = getc ( curinput ) ;
    if ( i == 13 ) 
        if ( ( i = getc ( curinput ) ) != 10  ) ungetc ( i , curinput ) ; 
    return ( laschar = i ) ; 
}

void openit ( void )
{
   struct stat buf ;
   curinput = inpf = fopen ( dargv [ 1 ] , "rb" ) ;
   if ( inpf == NULL ) {
       fprintf ( stderr , "Error trying to open input file %s: " , dargv[1] ) ;
       outer ( 1 , "cannot open file" ) ; }
   fstat ( fileno ( inpf ) , &buf ) ;
   dafsize = buf.st_size ; 
}
    
void outer ( int doit , char * txt )
{
    if ( !doit ) return ;
    fprintf ( stderr , "%s\n" , txt ) ;
    getc ( stdin ) ; 
    exit ( 0 ) ; 
}

void dohelp ( void )
{
    fprintf ( stderr , "\n\n"
"Usage is: \n\n"

"  Give input file name (1st arg) and options file (2nd arg)\n"
"  Output can be redirected with \">\"\n\n"

"  Inside options file: \n\n"

"    gene \"name(s)\"\n"
"    product \"name(s)\"\n"
"    protein\n"
"    genome \"type\"\n"
"    stringmatch S    (S=similarity=1-(E/L)\n"
"                      where E=edit cost, and L=length)\n\n"

"  List of skipped sequences goes to file gb2tnt.not\n\n"

"  To see details of skipped sequences, use gb2tnt.not as\n"
"  options file:\n\n"

"         \"gb2tnt input gb2tnt.not\"\n\n"

"  Alternatively, string similarity can be given as 3d and\n"
"  4th argument (which overrides values in options file)\n\n"

  ) ;
    getc ( stdin ) ; 
    exit ( 0 ) ; 
}

int main ( int argc , char ** argv )
{
    dargc = argc ;
    dargv = argv ;
    setopts () ;
    process () ; 
    return 0;
}    

void gonln ( void )
{
    while ( laschar != 10 ) rd () ;
}    

void gotostring ( char * string , int shift )
{
    char * cp ; 
    int i ;
    int somenonwhite ; 
    strcpy ( stringsp , " " ) ; 
    while ( strcmp ( stringsp , string ) && !feof ( inpf ) ) {
         gonln () ;
         somenonwhite = 0 ;
         if ( shift > 0 ) {
            for ( i = 0 ; i < shift && !somenonwhite ; ++ i )
                if ( rd () != 32 ) somenonwhite = 1 ;
            if ( somenonwhite ) continue ; }
         else {
            while ( isspace ( rd () ) ) ;
            ungetc ( laschar , inpf ) ; }
         * ( cp = stringsp ) = rd () ;
         if ( * cp ++ != * string ) continue ; 
         while ( !isspace( * cp = rd () ) ) {
             cp ++ ;
             if ( cp - stringsp > 98 ) break ; }
         * cp = '\0' ; }
    if ( feof ( inpf ) ) {
        fprintf ( stderr , "\nDone!" ) ; 
        exit ( 1 ) ; }
}

void rdaccnumber ( void )
{
    char * cp = accnumber ; 
    while ( isspace ( rd () ) ) ;
    while ( !isspace ( laschar ) ) {
         * cp ++ = laschar ;
         rd () ; }
    * cp = '\0' ;
}

void rdto ( int donde , char * stor )
{
    char * cp ;
    int a = 0 ; 
    cp = stor ;
    while ( isspace( laschar ) ) rd () ;
    while ( laschar != donde ) {
        if ( a != '_' || laschar != '_' ) 
           * cp ++ = laschar ;
        a = laschar ; 
        rd () ;
        if ( laschar == 32 ) laschar = '_' ;
        if ( laschar == 10 && donde != 10 ) laschar = '_' ; }
    * cp = '\0' ;
}    

void stornexline ( void )
{
    char * cp = stringsp ; 
    rd () ;
    while ( isspace ( laschar ) ) rd () ;
    while ( laschar != 10 && laschar != 13 ) {
        * cp ++ = laschar ;
        rd () ; }
    * cp = '\0' ; 
}

typedef struct 
      { int up , diag , lef ; 
        int min ; } Stringcomptyp ; 
Stringcomptyp ** cellcost ; 
int gapcost = 1 , gapextcost = 1 ; 
int suscost = 1 ;
int mademem = 0 ;

void ** loray ( int wid , int hei , int size ) 
{
  int a ; 
  void ** pp ; 
  pp = ( void ** ) malloc ( wid * sizeof ( void * ) ) ; 
  outer ( ( pp == NULL ) , "Not enough RAM") ;  
  for ( a = 0 ; a < wid ; ++ a ) 
    if ( ( pp [ a ] = ( void * ) malloc ( hei * size ) ) == NULL ) 
        outer ( 1 , "Not enough RAM") ;  
  return pp ; 
}

double doneedwunsch ( char * ap , char * bp ) 
{
    int wid , hei , i , j , dacos ; 
    char * app , * bpp ; 
    char * abecs , * bbecs , * anp , * bnp ;
    double val ; 
    int HIGH = 10000000 ; 
    wid = strlen ( ap ) ; 
    hei = strlen ( bp ) ; 
    if ( ! mademem ) {
         cellcost = ( Stringcomptyp ** ) loray ( 100 , 100 , sizeof ( Stringcomptyp ) ) ;
         mademem = 1 ; }
    outer ( ( hei >= 99 || wid >= 99 ) , "String is too long to use string-matching" ) ; 
    cellcost[0][0].min = cellcost[0][0].diag = 0 ; 
    cellcost[0][0].up = cellcost[0][0].lef = HIGH ; 
    bpp = bp ; 
    for ( j = 0 ; j < hei ; ++ j ) {
       app = ap ; 
       for ( i = 0 ; i < wid ; ++ i ) {
          if ( !i && !j ) {
             continue ; }
          dacos = 0 ; 
          if ( * app != * bpp ) dacos = suscost ; 
          if ( j ) {
            if ( cellcost[i][j-1].min == cellcost[i][j-1].up )  
               cellcost[i][j].up = cellcost[i][j-1].min + gapextcost ; 
            else 
               cellcost[i][j].up = cellcost[i][j-1].min + gapcost ; }
          else cellcost[i][j].up = cellcost[i][j].diag = HIGH ; 
          if ( i ) {
            if ( cellcost[i-1][j].min == cellcost[i-1][j].lef ) 
               cellcost[i][j].lef = cellcost[i-1][j].min + gapextcost ; 
            else 
               cellcost[i][j].lef = cellcost[i-1][j].min + gapcost ; }
          else cellcost[i][j].lef = cellcost[i][j].diag = HIGH ; 
          if ( i && j ) cellcost[i][j].diag = cellcost[i-1][j-1].min + dacos ; 
          dacos = cellcost[i][j].diag ; 
          if ( dacos > cellcost[i][j].up ) dacos = cellcost[i][j].up ; 
          if ( dacos > cellcost[i][j].lef ) dacos = cellcost[i][j].lef ; 
          cellcost[i][j].min = dacos ; 
          ++ app ; }
       ++ bpp ; }
    dacos = cellcost[wid-1][hei-1].min ;
    val = ( double ) dacos / ( double ) hei ;
    val = 1 - val ; 
    return val ; 
} 

int stringis ( char * a , char * b )
{

    if ( !strcmp ( b , "\"?\"" ) ) return 1 ; 
    if ( use_string_matching ) {
        if ( doneedwunsch ( a , b ) >= string_similarity ) return 1 ;
        return 0 ; }
    while ( * a && * b )
        if ( tolower ( * a ++ ) != tolower ( * b ++ ) ) return 0 ;
    return 1 ;
}

int rdliteral ( void )
{
    char * cp = stringsp ; 
    rd () ; 
    while ( isspace ( laschar ) && !feof ( curinput ) ) rd () ;
    if ( laschar == '\"' ) {
       * cp ++ = laschar ; 
       while ( rd () != '\"' && !feof ( curinput ) ) * cp ++ = laschar ;
       * cp ++ = laschar ;
       * cp = '\0' ;
       return 1 ; }
    else {
        * cp = laschar ;
        while ( !isspace ( laschar ) && !feof ( curinput ) ) * ++ cp = rd () ;
        * cp = '\0' ;
        return 0 ; }
}

int istrunc ( char * a )   // is a a truncation of stringsp ?? 
{
    char * b = stringsp ; 
    while ( 1 ) {
        if ( ! * a ) return 1 ;
        if ( * a ++ != * b ++ ) return 0 ; }
}

int isamatch ( char * a , char * b )
{
    while ( 1 ) {
        if ( ! * a ) return 1 ;
        if ( * a ++ != * b ++ ) return 0 ; }
}

void process (void)
{
    char * cp ;
    int i , mygenometype , showed_headerline ;
    int showacc_only ;
    int found_translation ; 
    while ( !feof ( inpf ) ) {
        gotostring ( "ACCESSION" , 0 ) ;
        rdaccnumber () ;
        showacc_only = 0 ; 
        if ( useraccnumber ) {
             for ( i = 0 ; i < useraccnumber && !showacc_only ; ++ i ) 
                if ( !strcmp ( accnumber , useraccname[i] ) ) showacc_only = 1 ; }
        laquiero = found_translation = 0 ;
        gotostring ( "ORGANISM" , -1 ) ;
        rdto ( 10 , taxname ) ;
        rdto ( '.' , taxonomy ) ; 
        gotostring ( "FEATURES" , 0 ) ;
        mygenometype = NUCLEAR ;
        /**** Empieza accession...   ****/
        while ( 1 ) {
          stornexline ( ) ;
          if ( istrunc ( "ORIGIN" ) ) { laquiero = -1 ; break ; } 
          if ( istrunc ( "/organelle=" ) ) {
              if ( istrunc ( "/organelle=\"mitoc" ) ) 
                   mygenometype = MITOCH ; 
              else if ( istrunc ( "/organelle=\"chlorop" ) ) 
                   mygenometype = CHLORO ; 
              else if ( istrunc ( "/organelle=\"plastid" ) ) 
                   mygenometype = PLASTID ; }
          if ( istrunc ( "tRNA" ) ) break ; 
          if ( istrunc ( "rRNA" ) ) break ; 
          if ( istrunc ( "CDS" ) ) {

/***  linea agregada por el caso de Marcos... ****/
if ( use_translation ) laquiero = FALSE ; 

              break ; }}  /***** End initial parsing... ******/
        if ( genometype && mygenometype != genometype )
             laquiero = -1 ; 
        while ( laquiero == FALSE || ( use_translation && !found_translation ) ) {
           if ( !istrunc ( "gene" ) ) {
                if ( istrunc ( "ORIGIN" ) ) break ;
                if ( istrunc ( "tRNA" ) || 
                     istrunc ( "rRNA" ) || 
                     istrunc ( "CDS" ) ) {
                     showed_headerline = 0 ; 
                     strcpy ( headerline , stringsp ) ; 
                     cp = stringsp + 10 ;
                     while ( isspace ( * cp ) ) ++ cp ;
                     if ( isamatch ( "join(" , cp ) ) {
                       while ( * cp != ')' && * cp ) ++ cp ;
                       if ( !*cp ) {
                         stornexline () ;
                         cp = stringsp ;
                         while ( isspace ( * cp ) ) ++ cp ; }
                       strcat ( headerline , cp ) ; }}
                stornexline () ;
                if ( showacc_only ) {
                     if ( showacc_only ++ == 1 ) {
                           fprintf ( stdout , "\n%s , %s " , accnumber , taxname ) ;
                           if ( mygenometype == NUCLEAR ) fprintf ( stdout , ", nuclear " ) ; 
                           if ( mygenometype == CHLORO ) fprintf ( stdout , ", chloro " ) ; 
                           if ( mygenometype == MITOCH ) fprintf ( stdout , ", mitoch " ) ; 
                           if ( mygenometype == PLASTID ) fprintf ( stdout , ", plastid " ) ; }
                     if ( !showed_headerline ) 
                        if ( isamatch ( "tRNA" , headerline ) || 
                             isamatch ( "rRNA" , headerline ) || 
                             isamatch ( "CDS" , headerline ) ) {
                                fprintf ( stdout , "\n   %s , " , headerline ) ;
                                ++ showed_headerline ; }
                     if ( ( istrunc ( "/product=" ) || istrunc ( "/gene=" ) ) )  {
                        if ( istrunc ( "/product=" ) ) cp = stringsp + 9 ;
                        else cp = stringsp + 6 ; 
                     fprintf ( stdout , ", %s" , cp ) ; }}
                if (

/***  linea agregada por el problema de Marcos.... ****/
laquiero &&

                istrunc ( "/translation=" ) && use_translation ) {
                     parse_translation () ;
                     found_translation = 1 ; }
                if ( istrunc ( "/product=" ) && prodnumber ) {
                    cp = stringsp + 9 ;
                    makelower ( cp ) ; 
                    for ( i = 0 ; i < prodnumber && !laquiero ; ++ i ) 
                        if ( stringis ( cp , prodname[i] ) ) laquiero = 1 ; }
                if ( istrunc ( "/gene=" ) && genenumber ) {
                    cp = stringsp + 6 ;
                    makelower ( cp ) ; 
                    for ( i = 0 ; i < genenumber && !laquiero ; ++ i ) 
                        if ( stringis ( cp , genename[i] ) ) laquiero = 1 ; }}
           else stornexline () ; }

        if ( useraccnumber ) {
            // if ( showacc_only == 2 ) fprintf ( stdout , "\n" ) ; 
            continue ; }
        if ( laquiero == TRUE && ( !use_translation || found_translation ) ) {
            if ( use_translation ) output_translation () ; 
            else parsit () ; }
        else 
            save_to_not () ; 
           }     // ---  while ( !feof ( inpf ) ) 
}

void save_to_not ( void )
{
       ++ rejected ; 
       if ( !showedskipped )
            fprintf ( notsfile , "accession " ) ;
       showedskipped = 1 ; 
       fprintf ( notsfile , "\"%s\" " , accnumber ) ;
}

void setopts ( void )
{
    int i ; 
    if ( dargc > 1 ) 
       if ( !strcmp ( dargv[1] , "--help" ) ) dohelp () ;
    outer ( dargc <  3 , "Specify input file and option file (\"--help\" to get help)" ) ; 
    outer ( ( curinput = opsf = fopen ( dargv[ 2 ] , "rb" ) ) == NULL , "Cannot open file with options" ) ;
    if ( dargc > 4 ) 
        if ( !strcmp ( dargv[3], "stringmatch" ) ) {
            use_string_matching = 1 ;
            string_similarity = atof ( dargv[4] ) ; }
    while ( !feof ( opsf ) ) {
        outer ( rdliteral () , "Unexpected literal string in options file" ) ;
        if ( !strcmp ( stringsp , "protein" ) ) use_translation = 1 ;
        else if ( !strcmp ( stringsp , "stringmatch" ) ) {
            rdliteral () ;
            if ( !use_string_matching ) {
                use_string_matching = 1 ;
                string_similarity = atof ( stringsp ) ; }}
        else if ( !strcmp ( stringsp , "gene" ) ) {
             while ( !feof ( curinput ) ) {
                 while ( isspace ( rd () ) && !feof ( curinput ) ) ;
                 ungetc ( laschar , opsf ) ; 
                 if ( laschar != '\"' ) break ;
                 outer ( genenumber == MAX_USER_DEFINITIONS  , "Cannot define so many gene names" ) ; 
                 rdliteral () ;
                 makelower ( stringsp ) ; 
                 genename [ genenumber ] = mymalloc ( ( strlen ( stringsp ) + 1 ) * sizeof ( char ) ) ;
                 strcpy ( genename [ genenumber ++ ] , stringsp ) ; }}
        else if ( !strcmp ( stringsp , "accession" ) ) {
             useraccname = mymalloc ( MAXUSERACC * sizeof ( char * ) ) ;
             while ( !feof ( curinput ) ) {
                 while ( isspace ( rd () ) && !feof ( curinput ) ) ;
                 ungetc ( laschar , opsf ) ; 
                 if ( laschar != '\"' ) break ;
                 rdliteral () ;
                 if ( useraccnumber == MAXUSERACC ) 
                     fprintf ( stderr , "Cannot define so many accessions, will skip from %s on" , stringsp ) ; 
                 else {
                    useraccname [ useraccnumber ] = mymalloc ( ( strlen ( stringsp ) + 1 ) *  sizeof ( char ) ) ;
                    strcpy ( useraccname [ useraccnumber ] , stringsp + 1 ) ;
                    i = strlen ( useraccname [ useraccnumber ] ) ;
                    useraccname [ useraccnumber ] [ i - 1 ] = '\0' ;
                    ++ useraccnumber ; }}}
        else if ( !strcmp ( stringsp , "product" ) ) {        
             while ( !feof ( curinput ) ) {
                 while ( isspace ( rd () ) && !feof ( curinput ) ) ;
                 ungetc ( laschar , opsf ) ; 
                 if ( laschar != '\"' ) break ;
                 outer ( prodnumber == MAX_USER_DEFINITIONS , "Cannot define so many product names" ) ; 
                 rdliteral () ;
                 makelower ( stringsp ) ; 
                 prodname [ prodnumber ] = malloc ( ( strlen ( stringsp ) + 1 ) *  sizeof ( char ) ) ;
                 strcpy ( prodname [ prodnumber ++ ] , stringsp ) ; }}
       else if ( !strcmp ( stringsp , "genome" ) ) {
           outer ( !rdliteral () , "Syntax error after \"genome\"" ) ;
           if ( !strcmp ( stringsp , "\"mitochondrial\"" ) ) genometype = MITOCH ;
           else if ( !strcmp ( stringsp , "\"nuclear\"" ) ) genometype = NUCLEAR ;
           else if ( !strcmp ( stringsp , "\"plastid\"" ) ) genometype = PLASTID ;
           else if ( !strcmp ( stringsp , "\"chloroplast\"" ) ) genometype = CHLORO ;
           else outer ( 1 , "Unrecognized genome option" ) ; }}
    if ( use_string_matching ) 
       fprintf ( stderr , "\nUsing string similarity of %.3f\n" , string_similarity ) ; 
    if ( !useraccnumber ) 
        outer ( ( notsfile = fopen ( "gb2tnt.not" , "wb" ) ) == NULL , "Cannot open file for skipped accessions" ) ; 
    openit () ;
}

void makelower ( char * txt )
{
    char * cp = txt ;
    while ( * cp ) {
         * cp = tolower ( *cp ) ;
         ++ cp ; }
    return ; 
}

int wrong_location ; 

char * storchunk ( char * cp )
{
     if ( * cp == '<' ) ++ cp ;
     if ( !isdigit ( * cp ) ) { wrong_location = 1 ; * cp = '0' ; return cp ; }
     chunk[numchunks].from = atoi ( cp ) ;
     while ( * cp ++ != '.' && * cp ) ;
     if ( * cp != '.' ) { wrong_location = 1 ; * cp = '0' ; return cp ; }
     ++ cp ;
     if ( * cp == '>' ) ++ cp ;
     if ( !isdigit ( * cp ) ) { wrong_location = 1 ; * cp = '0' ; return cp ; }
     chunk[numchunks].to = atoi ( cp ) ;
     while ( isdigit ( * cp ) ) ++ cp ; 
     ++ numchunks ;
     return cp ; 
}

int species_read = 0 ; 

void parsit ( void )
{
    char * cp = headerline + 15 ;
    int i , atchunk , atpos ;
    char * bytept = bytestring , now ;
    int complement_it = 0 , didntmatch ; 
    numchunks = 0 ; 
    wrong_location = 0 ; 
    while ( isspace ( * cp ) ) ++ cp ;
    if ( * cp == '<' || isdigit ( * cp ) ) storchunk ( cp ) ;
    else {
       i = cp [ 4 ] ;
       cp [ 4 ] = '\0' ;
       didntmatch = 0 ; 
       if ( strcmp ( "join" , cp ) ) didntmatch = 1 ; 
       cp [ 4 ] = i ;
       i = cp [ 10 ] ; 
       cp [ 10 ] = '\0' ;
       if ( didntmatch ) 
          if ( strcmp ( "complement" , cp ) ) {
              fprintf ( stderr , "OOPS!!\nFound unrecognized location specifier: %s\nFor accession %s\n" , cp , accnumber ) ;
              save_to_not () ; return ; }
       if ( !strcmp ( "complement" , cp ) ) {
           complement_it = 1 ;
           cp += 11 ; }
       else { cp [ 10 ] = i ; cp += 5 ; }
       while ( * cp != ')' && !wrong_location ) {
             cp = storchunk ( cp ) ;
             while ( * cp == ',' || isspace ( * cp ) ) ++ cp ; }}
    if ( wrong_location ) { save_to_not () ; return ; }
    gotostring ( "ORIGIN" , 0 ) ;
    rdliteral () ; 
    atchunk = atpos = 0 ;
    ++ species_read ; 
    while ( 1 ) {
        now = laschar ; 
        while ( laschar == 10 || laschar == 13 || laschar == 32 || ( laschar >= '0' && laschar <= '9' ) ) {
            now = laschar ; 
            rd () ; }
        if ( now == '/' || laschar == '/' ) return ; 
        ++ atpos ;
        if ( chunk[atchunk].from <= atpos && chunk[atchunk].to >= atpos ) {
            * bytept ++ = laschar ;
            outer ( ( bytept - bytestring >= MAXSEQLENGTH ) , "OOPS -- sequence is too long!\nChange MAXSEQLENGTH and re-compile" ) ; }
        laschar = 32 ; 
        if ( atpos == chunk[atchunk].to ) 
            if ( ++ atchunk == numchunks ) break ; }
    * bytept = '\0' ;
    ++ accepted ;
    spew_name () ; 
    if ( complement_it ) effect_complementation () ; 
    fprintf ( stdout , "%s" , bytestring ) ; 
    fprintf ( stdout , "\n" ) ; 
    fflush ( stdout ) ; 
}

void parse_translation ( void )
{
    char * cp , * bytept ; 
    cp = stringsp + 14 ;
    bytept = bytestring ;
    while ( * cp != '\"' && * cp ) {
           if ( !isspace ( * cp ) ) * bytept ++ = * cp ++ ;
           if ( * cp == 10 || * cp == 13 || * cp == '\0' ) {
                 stornexline () ;
                 cp = stringsp ;
                 while ( isspace ( * cp ) ) ++ cp ; }}
    * bytept = '\0' ;
}    

void spew_name ( void )
{
    char * cp = taxonomy , * here , * begg ;
    int numsemicols = 0 ; 
    while ( * cp != ';' && * cp ) ++ cp ;
    begg = here = ( cp += 2 ) ;
    while ( numsemicols < 2 && * cp ) {
        if ( * cp == ';' ) {
            ++ cp ;
            ++ numsemicols ; }
        * here ++ = * cp ++ ; }
   * here = '\0' ; 
    fprintf ( stdout , ">%s____%s_@%s\n" , taxname , accnumber , begg ) ;
}    

void output_translation ( void )
{
    ++ accepted ;             
    spew_name () ; 
    fprintf ( stdout , "%s" , bytestring ) ; 
    fprintf ( stdout , "\n" ) ; 
    fflush ( stdout ) ; 
}

char tmpmask[256] ;
char makeit[9] ;
char antimask[16] ; 

void effect_complementation ( void )
{
    int i , j , k , l ;
    int bit ; 
    int x , y ; 
    for ( i = 0 ; i < 256 ; ++ i ) tmpmask[i] = 0 ;
    tmpmask [ 'a' ] = tmpmask [ 'A' ] = 1 ; // tmpmask [ '0' ] ; 
    tmpmask [ 'g' ] = tmpmask [ 'G' ] = 2 ; // tmpmask [ '1' ] ; 
    tmpmask [ 'c' ] = tmpmask [ 'C' ] = 4 ; // tmpmask [ '2' ] ; 
    tmpmask [ 't' ] = tmpmask [ 'T' ] = 8 ; // tmpmask [ '3' ] ; 
    tmpmask [ 'R' ] = tmpmask [ 'r' ] = tmpmask [ 'a' ] | tmpmask [ 'g' ] ; 
    tmpmask [ 'Y' ] = tmpmask [ 'y' ] = tmpmask [ 't' ] | tmpmask [ 'c' ] ; 
    tmpmask [ 'W' ] = tmpmask [ 'w' ] = tmpmask [ 'a' ] | tmpmask [ 't' ] ; 
    tmpmask [ 'S' ] = tmpmask [ 's' ] = tmpmask [ 'c' ] | tmpmask [ 'g' ] ; 
    tmpmask [ 'M' ] = tmpmask [ 'm' ] = tmpmask [ 'a' ] | tmpmask [ 'c' ] ; 
    tmpmask [ 'K' ] = tmpmask [ 'k' ] = tmpmask [ 'g' ] | tmpmask [ 't' ] ; 
    tmpmask [ 'B' ] = tmpmask [ 'b' ] = tmpmask [ 'c' ] | tmpmask [ 'g' ] | tmpmask [ 't' ] ; 
    tmpmask [ 'D' ] = tmpmask [ 'd' ] = tmpmask [ 'a' ] | tmpmask [ 'g' ] | tmpmask [ 't' ] ; 
    tmpmask [ 'H' ] = tmpmask [ 'h' ] = tmpmask [ 'a' ] | tmpmask [ 'c' ] | tmpmask [ 't' ] ; 
    tmpmask [ 'V' ] = tmpmask [ 'v' ] = tmpmask [ 'a' ] | tmpmask [ 'c' ] | tmpmask [ 'g' ] ; 
    tmpmask [ 'N' ] = tmpmask [ 'n' ] = 1 | 2 | 4 | 8 ;
    makeit[tmpmask['a']] = tmpmask['t'] ;
    makeit[tmpmask['c']] = tmpmask['g'] ;
    makeit[tmpmask['g']] = tmpmask['c'] ;
    makeit[tmpmask['t']] = tmpmask['a'] ;
    for ( i = 0 ; i < 256 ; ++ i )
        if ( tmpmask[i] )
             antimask[ tmpmask[i] ] = i ; 
    j = strlen ( bytestring ) ;
    for ( i = 0 ; i < j ; ++ i )
       complement_string[i] = bytestring[i] ; 
    for ( i = 0 , k = j - 1 ; i < j ; ++ i , -- k ) {
        x = tmpmask [ complement_string [i]] ;
        y = 0 ;
        for ( l = 0 , bit = 1 ; l < 4 ; ++ l , bit <<= 1 ) 
            if ( ( bit & x ) ) 
                y |= makeit[ bit ] ;
        bytestring [ k ] = antimask [ y ] ; }
	return;
}
    

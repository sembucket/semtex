/*****************************************************************************
 * SLIT: reproduce specified columns of text file on stdout.                 *
 *                                                                           *
 * A limited form of awk?                                                    *
 *                                                                           *
 * Usage: slit [-c (ch1[,ch2,ch3, .....] [filename]                          *
 * Can be used as a filter.                                                  *
 *                                                                           *
 * Maximum number of columns which input file can contain:	MAXCOL       *
 * Maximum number of columns which can be output:		MAXCOL       *
 * Maximum number of digits in each column specifier:	        NUMDIG       *
 * Maximum length of string in each column:		        MAXSTR       *
 *                                                                           *
******************************************************************************/

/*------------------*
 * RCS Information: *
 *------------------*/
static char
  rcsid[] = "$Id$";

#include <stdio.h>
#include <ctype.h>
#include <nrutil.h>

#define MAXCOL	32
#define	MAXSTR	32
#define NUMDIG	8

int  parse   (char *strin, char *strout, int *pos, char sep);
void slit    (FILE *fp, int n, int *col);
int  getline (FILE *fp, char coltext[][MAXSTR], int *nwords);





main(int argc, char *argv[])
/*===========================================================================*
 * This does adminstration for the routines which do all the work.           *
 *===========================================================================*/
{
  static char usage[] =
    "usage: slit [options] [input]\n"
    "[options]:\n"
    "-h                  ... print this message\n"
    "-c col1[,col2,...]  ... reproduce these columns of file [Default: 1]\n";

  int	i=0, ncol=1, colnum[MAXCOL];
  char	c, num[NUMDIG];
  FILE *fp;
  
  colnum[0] = 0;

  while (--argc && (*++argv)[0] == '-') {
    switch (c = *++argv[0]) {
    case 'h':
      fputs(usage, stderr);
      exit(0);
      break;
    case 'c':
      ncol = 0;
      if (*++argv[0])
	while (parse(*argv, num, &i, ',') != 0 && ncol < MAXCOL) {
	  colnum[ncol] = atoi(num) - 1;
	  if (colnum[ncol] < 0) {
	    (void)fprintf(stderr, "slit: bad column\n");
	    exit(1);
	  }
	  ncol++;
	}
      else {
	++argv;	argc--;
	while (parse(*argv, num, &i, ',') != 0 && ncol < MAXCOL) {
	  colnum[ncol] = atoi(num) - 1;
	  if (colnum[ncol] < 0) {
	    (void)fprintf(stderr, "slit: bad column\n");
	    exit(1);
	  }
	  ncol++;
	}
      }
      break;
    default:
      message("slit", "unknown option", WARNING);
      fprintf(stderr, usage);
      exit(1);
      break;
    }
  }
 
  if (argc == 1) {
    if ((fp = fopen(*argv, "r")) == (FILE*) NULL)
      message("pdf", "couldn't open input file", ERROR);
    else { 
      slit(fp, ncol, colnum);
      fclose(fp);
    }
  }
  else 
    slit(stdin, ncol, colnum);

  exit(0);
}





int parse(char *strin, char *strout, int *pos, char sep)
/*===========================================================================*
 * Parse strin into strout until a non-solid character or separator          *
 * character occurs in strin or end of strin is reached.                     *
 *                                                                           *
 * Pos is the index of the first character in strin to be examined.          *
 * Update pos to point at the next non-separator character in strin.         *
 *                                                                           *
 * Parse returns the number of characters parsed from strin.                 *
 *===========================================================================*/
{
  int	i = 0;
  
  while (    (strin[*pos] != sep ) 
	 && !isspace(strin[*pos] ) 
	 &&  (strin[*pos] != '\0') ) {
    strout[i] = strin[*pos];
    (*pos)++;
    i++;
  }

  if (strin[*pos] == sep) (*pos)++;
  strout[i] = '\0';

  return i;
}





void slit(FILE *fp, int n, int *col)
/*===========================================================================*
 * Slit into columns, print up.                                              *
 *===========================================================================*/
{
  int	i, nwords;
  char	coltext[MAXCOL][MAXSTR];
  
  while (getline(fp, coltext, &nwords) != EOF) {
    if (nwords > 0) {
      (void)printf("%s", coltext[col[0]]);
      for (i=1; i<n; i++) 
	(void)printf(" %s", coltext[col[i]]);
    }
    (void)printf("\n");
  }
}





int getline(FILE *fp, char coltext[][MAXSTR], int *nwords)
/*===========================================================================*
 * The parsing of each input line is done here.                              *
 *===========================================================================*/
{
  int	c, i, rec, inword;
  
  rec    = -1;
  inword =  0;
  *nwords = 0;

  while ((c=getc(fp)) != EOF && c != '\n') {
    if (!isspace(c)) {
      if (!inword) {
	inword = 1;
	(*nwords)++;
	rec++;
	i = 0;
	coltext[rec][i] = c;
	i++;
      } else {
	coltext[rec][i] = c;
	i++;
      }
    } else
      inword = 0;

    if (rec >= 0) coltext[rec][i] = '\0';
  }

  return c;
}

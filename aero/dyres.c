/*****************************************************************************
 * Dyres: update a session.cou file to reflect new state-variable data.
 *
 * See the usage prompt below.  Dyres does rely on the fact that Aero reads
 * through .cou files to get the last state-variable data, so we can append
 * to the end of the file.
 *****************************************************************************/

static char
RCSid[] = "$Id$";

static char usage[] =
  "usage: dyres session[.fld]\n"
  "   or: dyres session.fld session.sta session.cou\n\n"
  "   or: dyres -h (...generates this message)\n";

#include <stdio.h>
#include <string.h>

#define NFIELD  14


int main (int argc, char *argv[])
/* ------------------------------------------------------------------------- *
 * This does the driving: no subroutines used.
 * ------------------------------------------------------------------------- */
{
  char   f_name[FILENAME_MAX], line[BUFSIZ];
  FILE  *fld_fp, *cou_fp, *sta_fp;
  double dtime;
  double step, time,
         xpos, xvel, xacc, xpf, xvf, xtf,
         ypos, yvel, yacc, ypf, yvf, ytf;
  int    found=0, nline=0;


  switch (argc) {
  case 2:
    if (strstr(argv[1], "-h")) {
      fprintf(stderr, usage);
      exit(0);
    } else {
      strcpy(f_name, argv[1]);
      if (!(fld_fp = fopen(f_name, "r")))
	if (!(fld_fp = fopen( strcat(f_name, ".fld"), "r" ))) {
	  fprintf(stderr, "dyres: Unable to open field file %s\n", f_name);
	  exit(1);
	}
      strcpy(f_name + strlen(f_name)-3, "sta");
      if (!(sta_fp = fopen(f_name, "r"))) {
	fprintf(stderr, "dyres: Unable to open .sta file %s\n", f_name);
	exit(1);
      }
      strcpy(f_name + strlen(f_name)-3, "cou");
      if (!(cou_fp = fopen(f_name, "a+"))) {
	fprintf(stderr, "dyres: Unable to open .cou file %s\n", f_name);
	exit(1);
      }
    }
    break;
  case 4:
    if (!(fld_fp = fopen(argv[1], "r"))) {
	fprintf(stderr, "dyres: Unable to open field file %s\n", f_name);
	exit(1);
      }
    if (!(sta_fp = fopen(argv[2], "r"))) {
	fprintf(stderr, "dyres: Unable to open .sta file %s\n", f_name);
	exit(1);
      }
    if (!(cou_fp = fopen(argv[3], "a+"))) {
	fprintf(stderr, "dyres: Unable to open .cou file %s\n", f_name);
	exit(1);
      }
    break;
  default:
    fprintf(stderr, usage);
    exit(1);
    break;
  }

  /* -- Extract dump time from field file (5th line of file---non-general?). */

  fgets(line, BUFSIZ, fld_fp);
  fgets(line, BUFSIZ, fld_fp);
  fgets(line, BUFSIZ, fld_fp);
  fgets(line, BUFSIZ, fld_fp);
  fgets(line, BUFSIZ, fld_fp);

  if (strstr(line, "Time") && !(strstr(line, "Step"))) {
    sscanf(line, "%lf", &dtime);
    fprintf(stderr, "Dump time: %.1f\n", dtime);
  } else {
    fprintf(stderr, "dyres: unable to find dump time in field file\n");
    exit(1);
  }
  
  fclose(fld_fp);

  /* -- Get state-variable info from state-variable file. */

  while (!feof(sta_fp) && !found) {
    fgets(line, BUFSIZ, sta_fp);
    nline++;
    if (line[0] == '#') continue;
    if (sscanf(line,
	       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	       &step, &time,
	       &xpos, &xvel, &xacc, &xpf, &xvf, &xtf,
	       &ypos, &yvel, &yacc, &ypf, &yvf, &ytf) != NFIELD) {
      fprintf(stderr, "dyres: terminated at incomplete line (%d)"
	      " in state-variable file\n", nline);
      exit(1);
    }
    found = (time == dtime);
  }
  
  fclose(sta_fp);

  if (!found) {
    fprintf(stderr, "dyres: unable to find dump time (%.1f)"
	    " in state-variable file\n", dtime);
    exit(1);
  }

  /* -- Now append state-variable information to coupling file. *

  fprintf(cou_fp, "x-state %14g %14g\n", xpos, xvel);
  fprintf(cou_fp, "y-state %14g %14g\n", ypos, yvel);

  fclose(cou_fp);

  return 0;
}

#ifndef STATS_H
#define STATS_H

/* ------------------------------------------------------------------------- *
 * stat_t                                                                    * 
 *                                                                           * 
 * This data structure is used to extract statistics for the primitive vari- * 
 * ables and Reynolds stresses.  Statistics are accumulated in physical      * 
 * space if IO_STAT > 0, where IO_STAT gives the frequency in time steps for * 
 * sending data through the statistics module, i.e.                          * 
 *                                                                           * 
 *     IO_STAT = 0      no statistics                                        * 
 *             = 1      accumulate every time step                           * 
 *             = 2      accumulate every other time step                     * 
 *                                                                           * 
 * and so forth.                                                             * 
 *                                                                           * 
 * Output                                                                    * 
 * ------                                                                    * 
 * At the very least, average  values of the primitive variables [uvwp] are  * 
 * accumulated each IO_STAT time step.  In addition, the independent comp-   * 
 * onents of the "Reynolds stress" tensor are accumulated.                   * 
 *                                                                           * 
 * The following naming convention is used for the Reynolds stresses:        * 
 *                                                                           * 
 *     | u u    u v    u w |       |  A      B      D  |                     * 
 *     | v u    v v    v w |   =   |  *      C      E  |                     * 
 *     | w u    w v    w w |       |  *      *      F  |                     * 
 *                                                                           * 
 * The '*' quantities are not computed because of symmetry.  Only the fields * 
 * [ABC] are generated for a 2D flow.  Note that these are not the true      * 
 * Reynolds stresses, but the raw data (in addition to the mean values of    * 
 * [uvw]) that are needed to compute the Reynolds stresses.   For example,   * 
 * the mean value of (u'v') can be computed as:                              * 
 *                                                                           * 
 *     <u'v'> = <(u-U)(v-V)> = <uv> - <U><V> = B - U V,                      * 
 *                                                                           * 
 * where B is defined above and U,V are the mean values of U(t) and V(t).    * 
 *                                                                           * 
 * Statistics are written to the file "session.avg", containing the values   * 
 * in physical space of the fields [uvpABC] or [uvwpABCDEF].  The value of   * 
 * "Step" in this file is the number of fields that have been accumulated in * 
 * the averages.                                                             * 
 *                                                                           * 
 * Initialization                                                            * 
 * --------------                                                            * 
 * When a stat_t is allocated, the statistics are initialized either to zero * 
 * or the values stored in the file "session.avg" if present.  This way,     *
 * multiple runs in the same directory will maintain running averages.       * 
 * ------------------------------------------------------------------------- */

/* $Id$ */

#include "speclib/field.h"

/* Macros */

#define STAT_NAME(s)  (s->name)
#define STAT_NFLD(s)  (s->nfld)
#define STAT_NAVG(s)  (s->navg)
#define STAT_SRC(s,i) (s->src[i])
#define STAT_AVG(s,i) (s->avg[i])

/* Declarations */

typedef struct stat {
  char  *name;        /* Name of the session file */
  int    nfld;        /* Number of fields being tracked */
  int    navg;        /* Number of averages accumulated */

  Field  **src;       /* Fields supplying the source data */
  Field  **avg;       /* Fields containing the averages */

  struct domain *d;   /* Computational domain */
} stat_t;

/* Prototypes */

stat_t *stat_alloc  (struct domain *d);
void    stat_free   (stat_t *s);
void    stat_update (stat_t *s);
void    stat_write  (stat_t *s);

#endif

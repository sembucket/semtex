#ifndef FAMILY_H
#define FAMILY_H

/* NB: This is a different type of "family" than what's defined as part of *
 *     the Element declaration.  This family is a group of elements that   *
 *     have the same geometric shape, i.e. the same map transforms each    *
 *     one to the bi-unit square.  It's stored as a linked list.           */

typedef struct family {             /* ........ FAMILY Definition ........ */
  int            id             ;   /* ID number of the parent element     */
  int            members        ;   /* Number of members                   */
  int            set            ;   /* Set flag (on = 1, off = 0)          */
  int            nr, ns         ;   /*                                     */
  double         **xr, **xs     ;   /* ----------------------------------- */
  double         **yr, **ys     ;   /*                                     */
  double         **jac          ;   /*       F A M I L Y   D A T A         */
  double         **rx, **ry     ;   /*             (Geometry)              */
  double         **sx, **sy     ;   /*                                     */
  double         **mass         ;   /* ----------------------------------- */
  struct element *parent        ;   /* Defining element for the family     */
  struct family  *next          ;   /* Pointer to the next family          */
} Family;

#endif

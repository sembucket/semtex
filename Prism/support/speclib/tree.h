#ifndef TREE_H
#define TREE_H

/* Tree code (named symbol tables)
 *
 * $Id$
 *
 * Copyright (c) 1998 R. D. Henderson and Caltech
 * 
 * ------------------------------------------------------------------------- */

typedef void (*PFV)(void*);  /* Pointer to a function returning void */

typedef enum { 
  Black, 
  Red 
} Color;

typedef struct  node  {      /* ...... Node ...... */
  char *        name    ;    /* Name of this node  */
  Color         color   ;    /* Node color         */
  void *        other   ;    /* Stored information */
  struct node * left    ;    /* Left  sub-tree     */
  struct node * right   ;    /* Right sub-tree     */
  struct node * parent  ;    /* Parent             */
} Node;

typedef struct tree {        /* ...... Tree ...... */
  struct node * root    ;    /* Root of the tree   */
  PFV           show    ;    /* Show a node        */
  PFV           kill    ;    /* Kill a node        */
} Tree;

/* Prototypes */

Node*  tree_search            (Node *x, char *key);
Node*  tree_minimum           (Node *x);
Node*  tree_maximum           (Node *x);
Node*  tree_successor         (Node *x);

void   tree_insert            (Tree *T, Node *x);
void   tree_walk              (Tree *T);

Node*  create_node            (char *key);
Tree*  create_tree            (PFV show, PFV kill);

#endif 

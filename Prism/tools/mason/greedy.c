/*
 * Bandwidth optimizer based on a Greedy-type algorithm
 * 
 * RCS Information
 * ------------------------
 * $Author$
 * $Date$
 * $RCSfile$
 * $Revision$
 * ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mason.h"


/* Internal variables */

typedef struct Queue {         /* .......... Element Queue .......... */
  int            id       ;    /* Queue entry ID number               */
  int            bw       ;    /* The bandwidth of this queue entry   */
  Element       *element  ;    /* Element stored in this queue entry  */
  struct Queue  *prev     ;    /* Previous element in the queue       */
  struct Queue  *next     ;    /* Next element in the queue           */
} Queue, *QueueP;

static int (*indexer)(Element *U, Domain *Omega);

/* Internal functions */

static int  
  optimize_xx   (Domain *omega, int start),
  bandwidth_est (Element *U, Domain *Omega),
  cross_bw      (Domain *omega, Edge *edge, Node pos, Node *min, Node *max);

static void 
  reset         (Domain *omega),
  save_ops      (Queue  *queue);

static QueueP
  copy          (Queue *queue, int id, Queue *dest),
  delete        (Queue *queue, int id),
  insert        (Queue *queue, int id, Element *U),
  lookup        (Queue *queue, int id),
  randomize     (Queue *queue),
  qindex        (Queue *queue, Domain *omega),
  showQueue     (Queue *queue, FILE *fp);


/* ------------------------------------------------------------------------ */

void greedy (Domain *omega)
{
  Element *U;
  int bw, best, id, i;

  indexer = bandwidth_est;

  if (rflag) {               /* ...keep the random numbering */
    optimize_xx (omega, -1);
    return;                  /* ...done, with numbering random */
  }
   
  /* Save the bandwidth of the initial numbering */
    
  best = optimize_xx (omega, id = 1);

  /* Search through all boundary elements */
  
#ifdef TEST
  printf ("Bandwidth [init] = %d\n", bw);
#endif

  for (U = omega->U; U; U = U->next) {
    for (i = 0; i < U->edges; i++)
      if (strchr("EPMS", U->elist[i].type) == NULL) {
	bw = optimize_xx (omega, U->id);
	
	if (omega->verbose)
	  fprintf (stderr, "mason: id = %3d, bw = %3d\n", U->id, bw);
	
	if (bw < best) {
	  best = bw;
	  id   = U->id;
	}
	
	break;
      }
  }
  
#ifdef TEST
  printf ("Bandwidth [best] = %d from node %d\n", best, id);
  exit (1);
#endif
  
  optimize_xx (omega, opstart = id);

  return;
}


/* Do an optimization pass starting from element "start" */

static int optimize_xx (Domain *omega, int start)
{
  int     bw;
  Queue   *todo, *done, *qmax;
  Element *U;
  
  /* Set up the queues of elements to number */

  reset  (omega);
  done = (QueueP) NULL;
  todo = (QueueP) NULL;
  for (U = omega->U; U ; U = U->next)
    todo = insert (todo, U->id, U);
  
  /* Randomize the queue to get rid of bias */
  
  todo = randomize (todo);

  if (start < 0) {
  
    /* If start < 0, return the "random" numbering bandwidth */

    while (todo) {
      linkNodesE (omega, todo->element);
      done  = copy   (todo, todo->id, done);
      todo  = delete (todo, todo->id);
    }

    /* Uniform numbering, no optimization */

  } else if (start == 0) {
    for (U = omega->U; U ; U = U->next) {
      linkNodesE (omega, U);
      done  = copy   (todo, U->id, done);
      todo  = delete (todo, U->id);
    }
  
  /* Otherwise, perform a normal optimization pass */

  } else {
  
    linkNodesE    (omega, lookup(todo, start)->element);
    done = copy   (todo,  start, done);
    todo = delete (todo,  start);
    
    while (todo) {
      qmax = qindex (todo,  omega);
      linkNodesE    (omega, qmax->element);
      done = copy   (todo,  qmax->id, done);
      todo = delete (todo,  qmax->id);
    }
  }
  

  /* Save the final bandwidth */

  bw = qindex (done, omega)->bw;
  

  /* Save the optimizer results in a log file */

  if (omega->verbose) save_ops (done);  


  /* Clean up */

  while (done)
    done = delete (done, done->id);

  return bw;
}


/* ------------------------------------------------------------------------ *
 *                    P R I V A T E     F U N C T I O N S                   *
 * ------------------------------------------------------------------------ */

/* Write the optimizer results to a file */

static void save_ops (Queue *queue)
{
  FILE *fp;

  if (!(fp = fopen (".mason.opt", "w")))
    fprintf (stderr, "optimizer: warning: unable to open a log file\n");
  else  {
    showQueue (queue, fp);
    fclose    (fp);
  }

  return;
}

/* Reset all nodes to NULL */

static void reset (Domain *omega)
{
  Element *U = omega->U;
  Edge    *edge;
  int i;

  while (U) {
    for (edge = U->elist; edge; edge = edge->next) {
      const int dir = edge->dir;
      const int np  = edge->np;
      int       *p  = edge->nodes;

      for (i = 0; i < np; i++, p += dir)
 	  *p = (Node) NULL;

      *(edge-> right ->node) = (Node) NULL;
      *(edge-> left  ->node) = (Node) NULL;
    }
    U = U->next;
  }

  omega->nodes = 0;
  return;
}

/* ------------------------------------------------------------------------ *
 * Queue Functions                                                          *
 *                                                                          *
 * The following function manage the simple queues for number the elements. *
 * Each function modifies a queue and returns a pointer to the new queue.   *
 * The beginning and end of the queue are signaled by a NULL pointer.       *
 * ------------------------------------------------------------------------ */

static QueueP insert (Queue *queue, int id, Element *U)
{
  QueueP new   = (QueueP) calloc(1,sizeof(Queue));
  QueueP q     =  queue;

  new->id      = id;
  new->element = U;

  /* Search for the end of this queue and append to it */

  if (queue) {
    while (q->next) q = q->next;

    q  ->next = new;
    new->prev = q;
  }

  /* This is the first entry */
  
  else
    queue = new;
  
  return queue;
}

/* Delete and entry from the queue */
    
static QueueP delete (Queue *queue, int id)
{
  QueueP p, n, q;

  if ((q=lookup(queue, id))) {

    p = q->prev;
    n = q->next;

    free (q);

    if (n) 
      n->prev = p;    /* Is this the last element in the queue? */
    if (p) 
      p->next = n;    /* Is this an entry from the middle? */
    else
      queue   = n;    /* Is this the first element in the queue? */
  }    
  
  return queue;       /* Return the modified queue */
}

/* Search for an entry in the queue */

static QueueP lookup (Queue *queue, int id)
{
  QueueP q;

  if (!(q = queue))     /* If the queue is empty, return a NULL */
    return queue;

  while (q)
    if (q->id == id)
      return q;
    else
      q = q->next;

  sprintf   (error_buf, "lookup: Queue ID = %d not found", id);
  error_msg (error_buf);
  return (QueueP) NULL;
}

/* Copy from one queue to another */

static QueueP copy (Queue *orig, int id, Queue *dest)
{
  QueueP q;

  q     = lookup (orig, id);
  dest  = insert (dest, id, q->element);

  return dest;
}

/* Show the entries in a queue */

static QueueP showQueue (Queue *queue, FILE *fp)
{
  QueueP q = queue;

  fputs ("# Queue ID   Element ID   Bandwidth\n" 
	 "# ---------------------------------\n", fp);

  while (q) {
    fprintf (fp, "%8d     %8d    %8d\n", q->id, q->element->id, q->bw);
    q = q->next;
  }

  return (QueueP) NULL;
}


/* Compute the queueing index (bandwidth), return the maximum */

#define MAX(a,b) ((a > b) ? (a) : (b))
#define MIN(a,b) ((a > b) ? (b) : (a))

static QueueP qindex (Queue *queue, Domain *omega)
{
  Queue *q, *qmax;

  if (!queue) return queue;
  
  qmax     = queue;    /* Start up the search for the maximum */
  qmax->bw = 0;

  for (q = queue; q ; q = q->next)
    if ((q->bw = (*indexer)(q->element, omega)) > qmax->bw)
      qmax = q;

  return qmax;
}

/* Estimate what an element's bandwidth will be if numbered */

static int bandwidth_est (Element *U, Domain *omega)
{
  Edge *e;
  int   p;
  
  Node min = 0xfffff;
  Node max = -min;
  int  pos = omega->nodes + 1;

  for (e = U->elist; e ; e = e->next) {
    const int n = (e->np-1)*(e->dir);

    if (!(p = *(e->right->node)))
      p = pos++;
    max = MAX (p, max);
    min = MIN (p, min);

    if (e->type == 'S') 
      pos = cross_bw (omega, e, pos, &min, &max);
    else {
      if (!(p = *(e->nodes + n)))
	p = (pos += e->np);
      max = MAX (p, max);
      min = MIN (p, min);
    }

    if (!(p = *(e->left->node)))
      p = pos++;
    max = MAX (p, max);
    min = MIN (p, min);
  }

  return max-min;
}


#if 0

/* Estimate an element's bandwidth relative to its neighbors */

static int bandwidth_zero (Element *U, Domain *omega)
{
  Edge *e;

  int min = 0xfffff;
  int max = -min;
  int p;

  for (e = U->elist; e ; e = e->next) {
    const int n = (e->np-1)*(e->dir);

    p   = *(e->right->node);
    max = MAX (p, max);
    min = MIN (p, min);

    p   = *(e->nodes + n);
    max = MAX (p, max);
    min = MIN (p, min);
    
    p   = *(e->left->node);
    max = MAX (p, max);
    min = MIN (p, min);
  }
  
  return max-min;
}

#endif

static Node cross_bw 
  (Domain *omega, Edge *edge, Node pos, Node *pmin, Node *pmax)
{
  SegmentP master;
  int      b, id, n, p;

  Node min = 0xfffff;
  Node max = -min;

  PatchP   patch = findPatch   (edge->bc.p.patch,   omega);
  SegmentP slave = findSegment (edge->bc.p.segment, patch->slaves);

  for (b = 0; b < (*slave).branches; b++) {

    if (!(master = findSegment (id = slave->branch_ID[b], patch->masters))) {
      sprintf   (error_buf, "cross_bw: failed searching for segment %d.%dm", 
		 patch->id, id);
      error_msg (error_buf);
    }
      
    n = (master->edge->np-1)*(master->edge->dir);

    if (!(p = *master->edge->right->node))
      p = pos++;
    max = MAX (p, max);
    min = MIN (p, min);

    if (!(p = *(master->edge->nodes + n)))
      p = (pos += master->edge->np);
    max = MAX (p, max);
    min = MIN (p, min);
    
    if (!(p = *master->edge->left->node))
      p = pos++;
    max = MAX (p, max);
    min = MIN (p, min);
  }

  *pmax = MAX(*pmax,max);
  *pmin = MIN(*pmin,min);

  return pos;
}

#undef MAX
#undef MIN

/* Randomize a queue */

static int qcomp (const void *, const void *);

static QueueP randomize (Queue *queue)
{
  Queue *q, *q_sort;
  int    n, nq;
  
  for (nq = 0, q = queue; q; q = q->next)
       nq++;

  q_sort = (Queue*) calloc (nq, sizeof(Queue));
  for (n = 0, q = queue; q; n++) {
    memcpy (q_sort + n, q, sizeof(Queue));
    q_sort[n].bw = rand();
    queue = q;
    q = q->next;
    free (queue);
  }
  
  qsort (q_sort, nq, sizeof(Queue), qcomp);

  for (n = 0; n < nq; n++) {
    q_sort[n].next = q_sort + n + 1;
    q_sort[n].prev = q_sort + n - 1;
  }

  q_sort[ 0  ].prev = NULL;
  q_sort[nq-1].next = NULL;

  return q_sort;
}

static int qcomp (const void *p1, const void *p2)
{
  const int bw1 = ((Queue*) p1)->bw;
  const int bw2 = ((Queue*) p2)->bw;

  if (bw1 > bw2)
    return  1;
  if (bw1 < bw2)
    return -1;

  return 0;
}

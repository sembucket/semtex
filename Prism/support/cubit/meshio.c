/*
 * Mesh I/O 
 *
 * $Id$
 * ------------------------------------------------------------------------ */

#include <stdlib.h>
#include <assert.h>
#include "cubit.h"
#include "mesh.h"

int Mesh_import (Mesh *mesh, FILE *fp)
{
  Element *list = ReadMesh(fp);
  Element *elmt = list;

  while (elmt) {
    Element *next = elmt->next;
    Mesh_add(mesh,elmt);
    elmt = next;
  }

  mesh->bc = ReadBCs(fp, 1, mesh->head);
  mesh->nr = mesh->head->nr;
  mesh->ns = mesh->head->ns;
  mesh->nz = mesh->head->nz;

  ReadKeys(fp, mesh->head);

  Mesh_connect(mesh);

  return 0;
}

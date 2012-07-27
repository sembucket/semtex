///////////////////////////////////////////////////////////////////////////////
// image.C: routines for image processing/dumping.
//
// Copyright (c) 1999 <--> $Date$, Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "sview.h"


int writetiff (char* filename   ,
	       char* description,
	       int   compression)
// ----------------------------------------------------------------------------
// This code is adapted from a routine that comes with the GLUT
// distribution.  It uses libtiff to write a TIFF file corresponding
// to the current window.  Note that libtiff uses stdio.
//
// Legal values for compression are: COMPRESSION_LZW, COMPRESSION_NONE
// and COMPRESSION_PACKBITS.
//
// OpenGL's default 4 byte pack alignment would leave extra bytes at
// the end of each image row so that each full row contained a number
// of bytes divisible by 4.  Ie, an RGB row with 3 pixels and 8-bit
// componets would be laid out like "RGBRGBRGBxxx" where the last
// three "xxx" bytes exist just to pad the row out to 12 bytes (12 is
// divisible by 4). To make sure the rows are packed as tight as
// possible (no row padding), set the pack alignment to 1.
// ----------------------------------------------------------------------------
{
  TIFF         *file;
  GLubyte      *image, *p;
  GLint        *vp;
  register int i;

  file = TIFFOpen (filename, "w");
  if (file == NULL) { return 1; }

  glPixelStorei (GL_PACK_ALIGNMENT, 1);

  vp = (GLint*) malloc (4 * sizeof (GLint));
  glGetIntegerv (GL_VIEWPORT, vp);

  image = (GLubyte *) malloc (vp[2] * vp[3] * sizeof(GLubyte) * 3);

  glutSwapBuffers();
  glReadPixels   (vp[0], vp[1], vp[2], vp[3], GL_RGB, GL_UNSIGNED_BYTE, image);

  TIFFSetField (file, TIFFTAG_IMAGEWIDTH,       (uint32) vp[2]     );
  TIFFSetField (file, TIFFTAG_IMAGELENGTH,      (uint32) vp[3]     );
  TIFFSetField (file, TIFFTAG_BITSPERSAMPLE,    8                  );
  TIFFSetField (file, TIFFTAG_COMPRESSION,      compression        );
  TIFFSetField (file, TIFFTAG_PHOTOMETRIC,      PHOTOMETRIC_RGB    );
  TIFFSetField (file, TIFFTAG_ORIENTATION,      ORIENTATION_TOPLEFT);
  TIFFSetField (file, TIFFTAG_SAMPLESPERPIXEL,  3                  );
  TIFFSetField (file, TIFFTAG_PLANARCONFIG,     PLANARCONFIG_CONTIG);
  TIFFSetField (file, TIFFTAG_ROWSPERSTRIP,     1                  );
  TIFFSetField (file, TIFFTAG_IMAGEDESCRIPTION, description        );

  p = image;
  for (i = vp[3] - 1; i >= 0; i--) {
    //    if (TIFFWriteScanline (file, p, i, 0) < 0) {
    if (TIFFWriteScanline (file, p, vp[3] - i - 1, 0) < 0) {
      free (vp); free (image); TIFFClose (file); return 1;
    }
    p += vp[2] * sizeof (GLubyte) * 3;
  }
  free (vp); free (image); TIFFClose (file); return 0;
}

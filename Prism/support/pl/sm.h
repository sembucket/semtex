#ifndef SM_H
#define SM_H

#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for SM-callable functions */

void sm_defvar   (const char *token, const char *expr);
void sm_device   (const char *);
void sm_limits   (double, double, double, double);
void sm_box      (int, int, int, int);
void sm_gflush   (void);
void sm_erase    (void);
void sm_graphics (void);
void sm_conn     (const float [], const float [], int);
void sm_expand   (double);
void sm_angle    (double);
void sm_ptype    (const float *, int);
void sm_points   (const float [], const float [], int);
void sm_shade    (int, const float [], const float [], int);
void sm_relocate (double, double);
void sm_label    (const char *label);
void sm_xlabel   (const char *label);
void sm_ylabel   (const char *label);
void sm_levels   (const float [], int);
void sm_contour  (void);
void sm_defimage (const float **, double, double, double, double, int, int);
void sm_delimage (void);
void sm_dot      (void);
void sm_alpha    (void);
void sm_curs     (float *x, float *y, int *n);
void sm_line     (double,double,double,double);
void sm_draw     (double,double);
void sm_ltype    (int);
void sm_connect  (float [], float [], int);

#ifdef __cplusplus
}
#endif

#endif

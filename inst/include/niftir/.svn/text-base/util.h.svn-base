#ifdef __cplusplus
extern "C" {
#endif

#include "niftir/nifti1_io.h"

#include <Rdefines.h>
#include <Rinternals.h>

SEXP mat44_to_SEXP(mat44 *mat);

SEXP mat33_to_SEXP(mat33 *mat);

SEXP float_to_SEXP(float val);

SEXP short_to_SEXP(short val);

SEXP int_to_SEXP(int val);

SEXP pchar_to_SEXP(const char* val);

/****************************************
 * Functions to convert from R => C type
 ****************************************/

void SEXP_to_mat44(SEXP val, mat44 *mat);

void SEXP_to_mat33(SEXP val, mat33 *mat);

void SEXP_to_float(SEXP val_sexp, float *val);

void SEXP_to_short(SEXP val_sexp, short *val);

void SEXP_to_int(SEXP val_sexp, int *val);

void SEXP_to_pchar(SEXP val_sexp, char* val, int max_num);


#ifdef __cplusplus
}
#endif

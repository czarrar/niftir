#include "niftir/util.h"

/****************************************
 * Functions to convert from C to R type
 ****************************************/

SEXP mat44_to_SEXP(mat44 *mat) {
    SEXP ret_val;
    int c,r;
    PROTECT(ret_val = NEW_NUMERIC(4*4));
    for(r=0;r<4;++r)
        for(c=0;c<4;++c)
            NUMERIC_POINTER(ret_val)[r+c*4]=(double)(mat->m[r][c]);

    SEXP dim;
    PROTECT(dim=NEW_INTEGER(2));
    INTEGER_POINTER(dim)[0]=4;
    INTEGER_POINTER(dim)[1]=4;
    SET_DIM(ret_val,dim);

    UNPROTECT(2);
    return ret_val;
}

SEXP mat33_to_SEXP(mat33 *mat) {
    SEXP ret_val;
    int c,r;
    PROTECT(ret_val = NEW_NUMERIC(4*4));
    for(r=0;r<3;++r)
        for(c=0;c<3;++c)
            NUMERIC_POINTER(ret_val)[r+c*4]=(double)mat->m[r][c];
    
    SEXP dim;
    PROTECT(dim=NEW_INTEGER(2));
    INTEGER_POINTER(dim)[0]=3;
    INTEGER_POINTER(dim)[1]=3;
    SET_DIM(ret_val,dim);

    UNPROTECT(2);
    return ret_val;
}

SEXP float_to_SEXP(float val) {
    SEXP ret_val;
    PROTECT(ret_val=NEW_NUMERIC(1));
    NUMERIC_POINTER(ret_val)[0]=val;
    UNPROTECT(1);
    return ret_val;
}

SEXP short_to_SEXP(short val) {
    SEXP ret_val;
    PROTECT(ret_val=NEW_NUMERIC(1));
    NUMERIC_POINTER(ret_val)[0]=val;
    UNPROTECT(1);
    return ret_val;
}

SEXP int_to_SEXP(int val) {
    SEXP ret_val;
    PROTECT(ret_val=NEW_INTEGER(1));
    INTEGER_POINTER(ret_val)[0]=val;
    UNPROTECT(1);
    return ret_val;
}

SEXP pchar_to_SEXP(const char* val) {
    SEXP ret_val;
    PROTECT(ret_val=NEW_CHARACTER(1));
    if(val!=NULL)
        SET_STRING_ELT(ret_val, 0, mkChar(val));
    else
        SET_STRING_ELT(ret_val, 0, mkChar(""));
    UNPROTECT(1);
    return ret_val;
}

/****************************************
 * Functions to convert from R => C type
 ****************************************/
 
void SEXP_to_mat44(SEXP val, mat44 *mat) {
    int c,r;
    PROTECT(val = AS_NUMERIC(val));
    if (LENGTH(val)==16) {
        for(r=0;r<4;++r)
    	    for(c=0;c<4;++c)
    	        (mat->m[r][c])=(float)NUMERIC_POINTER(val)[r+c*4];
    } else {
        error("matrix must be 4x4\n");
    }
    UNPROTECT(1);    
}

void SEXP_to_mat33(SEXP val, mat33 *mat) {
    int c,r;
    PROTECT(val = AS_NUMERIC(val));
    if(LENGTH(val)==9) {
        for(r=0;r<3;++r)
    	    for(c=0;c<3;++c)
    	        (mat->m[r][c])=(float)NUMERIC_POINTER(val)[r+c*4];
    } else {
        error("matrix must be 3x3\n");
    }
    UNPROTECT(1);    
}

void SEXP_to_float(SEXP val_sexp, float *val) {
    PROTECT(val_sexp=AS_NUMERIC(val_sexp));
    *val=(float)NUMERIC_POINTER(val_sexp)[0];
    UNPROTECT(1);
}

void SEXP_to_short(SEXP val_sexp, short *val) {
    PROTECT(val_sexp=AS_NUMERIC(val_sexp));
    *val=(short)NUMERIC_POINTER(val_sexp)[0];
    UNPROTECT(1);
}

void SEXP_to_int(SEXP val_sexp, int *val) {
    PROTECT(val_sexp=AS_INTEGER(val_sexp));
    *val=(int)INTEGER_POINTER(val_sexp)[0];
    UNPROTECT(1);
}

void SEXP_to_pchar(SEXP val_sexp, char* val, int max_num) {
    PROTECT(val_sexp=AS_CHARACTER(val_sexp));
    const char *pcstring = CHAR(CHARACTER_POINTER(val_sexp)[0]);
    if(strlen(pcstring)<max_num)
        strcpy(val,pcstring);
    else
        error("character string too long\n");
    UNPROTECT(1);
}

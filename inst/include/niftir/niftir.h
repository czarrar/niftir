#include "niftir/util.h"

#ifdef __cplusplus
extern "C" {
#endif

#define GET_ELEMENT(x, i)	VECTOR_ELT(x, i)
#define GET_LIST_ELEMENT(x, s)  getListElement(x, s)
#define NIRTS_FTYPE NIFTI_FTYPE_NIFTI1_1

extern char *gni_version;

#define NUM_HEADER_ATTRS 31

const char *header_attributes[] = {
    "toffset", 			/* 1 */
    "descrip", 			/* 2 */
    "fname",   			/* 3 */
    "iname",   			/* 4 */
    "slice.duration", 	/* 5 */
    "qform.code", 		/* 6 */
    "sform.code", 		/* 7 */
    "quatern.b", 		/* 8 */
    "quatern.c", 		/* 9 */
    "quatern.d", 		/* 10 */
    "qoffset.x", 		/* 11 */
    "qoffset.y", 		/* 12 */
    "qoffset.z", 		/* 13 */
    "qfac",			    /* 14 */
    "pixdim",			/* 15 */
    "nifti.type",		/* 16 */
    "sizeof.hdr",	    /* 17 */
    "datatype",         /* 18 */
    "scl.slope",        /* 19 nifti1: Data scaling: slope.  analyze 7.5: float funused1; */
    "scl.inter",        /* 20 nifti1: Data scaling: offset. analyze 7.5: float funused2; */
    "qto.xyz", 			/* 21 */
    "qto.ijk", 			/* 22 */
    "sto.xyz", 			/* 23 */
    "sto.ijk", 			/* 24 */    
    "dim",              /* 25 */
    "xyz.units",        /* 26 */
    "time.units",       /* 27 */
    "intent.code",      /* 28 */
    "intent.name",      /* 29 */
    "cal.min",          /* 30 */
    "cal.max",          /* 31 */
    NULL
};

SEXP getListElement(SEXP list, const char *str);


/****************************************
 * Helper functions
 ****************************************/

SEXP getHeaderAttributes(void);

//SEXP getXformString(SEXP xx) {
//    SEXP Rstr;
//    
//    const char *Cstr = nifti_xform_string(SEXP_to_int(xx));
//    SEXP_to_pchar(Rstr, Cstr, 100);
//    
//    return Rstr;
//}
 

/****************************************
 * Niftilib wrapper functions
 ****************************************/

nifti_image *niftir_image_read(SEXP file, SEXP read_data);

SEXP niftir_image_free(SEXP nim);

SEXP get_nifti_header(nifti_image *pnim);

void set_nifti_header(nifti_image *pnim, SEXP header);

nifti_image *create_nifti_image(SEXP header, SEXP Rdim, SEXP Rdatatype, SEXP outfile);

SEXP read_nifti(SEXP file, SEXP read_data);

SEXP write_nifti(SEXP header, SEXP data, SEXP outfile);


/****************************************
 * BigNifti Read/Write Functions
 ****************************************/

// Reading Functions

SEXP read_bignifti_header(SEXP file);

SEXP read_bignifti_data(SEXP nim_addr, SEXP big_addr);

SEXP read_partial_bignifti_data(SEXP nim_addr, SEXP big_addr, \
                                SEXP rowIndices, SEXP colIndices, \
                                SEXP totalVoxs);

// Writing Functions

SEXP write_bignifti(SEXP header, SEXP big_addr, SEXP indices, SEXP outfile);


#ifdef __cplusplus
}
#endif

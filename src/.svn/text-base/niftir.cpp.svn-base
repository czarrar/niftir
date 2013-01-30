#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/isna.hpp"
#include "bigmemory/util.h"

#include "niftir/niftir.h"
#include "niftir/nifti1.h"

#include <utility>
#include <vector>
#include <algorithm>
#include <map>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

// Masking function: pass two big matrices, pass the indices of first to mask

template<typename bT, typename dT, typename MatrixAccessorType>
void NiftiDataToBigMatrix(void *pnim_data, MatrixAccessorType m, index_type ncols, index_type nrows) {
    index_type i=0;
    index_type j=0;
    index_type k=0;
    bT *pColumn;
        
    for (i=0; i < nrows; ++i)
  	    for (j=0; j < ncols; ++j, ++k)       
            m[j][i] = (bT)((dT)pnim_data)[k];
}

template<typename bT, typename MatrixAccessorType>
SEXP NiftiToBigMatrix(nifti_image *pnim, BigMatrix *pMat, MatrixAccessorType m) {
    // Loop through columns and rows; assign the elements in pnim that correspond
    switch (pnim->datatype) {
        case NIFTI_TYPE_UINT8:
            NiftiDataToBigMatrix<bT, unsigned char*>(pnim->data, m, pMat->ncol(), pMat->nrow());
            break;
        case NIFTI_TYPE_INT8:
            NiftiDataToBigMatrix<bT, signed char*>(pnim->data, m, pMat->ncol(), pMat->nrow());
            break;            
        case NIFTI_TYPE_UINT16:
            NiftiDataToBigMatrix<bT, unsigned short*>(pnim->data, m, pMat->ncol(), pMat->nrow());
            break;            
        case NIFTI_TYPE_INT16:
            NiftiDataToBigMatrix<bT, short*>(pnim->data, m, pMat->ncol(), pMat->nrow());
            break;            
        case NIFTI_TYPE_UINT32:
            NiftiDataToBigMatrix<bT, unsigned int*>(pnim->data, m, pMat->ncol(), pMat->nrow());
            break;            
        case NIFTI_TYPE_INT32:
            NiftiDataToBigMatrix<bT, int*>(pnim->data, m, pMat->ncol(), pMat->nrow());
            break;
        case NIFTI_TYPE_FLOAT32:
            NiftiDataToBigMatrix<bT, float*>(pnim->data, m, pMat->ncol(), pMat->nrow());
            break;
        case NIFTI_TYPE_FLOAT64:
            NiftiDataToBigMatrix<bT, double*>(pnim->data, m, pMat->ncol(), pMat->nrow());
            break;
    	default:
    	    warning("unsupported data format (identifier %d)", pnim->datatype);
            break;
    }
    
    return R_NilValue;
}

template<typename CType, typename NType, typename BMAccessorType>
void BigMatrixToNiftiData(nifti_image *pnim, BigMatrix *pMat, SEXP indices) {
    BMAccessorType m( *pMat );
    
    int *pCols = INTEGER_DATA(indices);
    index_type numCols = GET_LENGTH(indices);
    if (numCols != pMat->ncol())
        error("indices must have the same length as number of columns in big matrix");
    index_type numRows = pMat->nrow();
    
    index_type i=0;
    index_type j=0;
    CType *pColumn;
    
    index_type numVoxs = (index_type)(pnim->dim[1]*pnim->dim[2]*pnim->dim[3]);
    
    for (i = 0; i < numCols; ++i) {
        pColumn = m[i];
        for (j = 0; j < numRows; ++j) {
            ((NType*)pnim->data)[(j*numVoxs)+(static_cast<index_type>(pCols[i])-1)] = (NType)(pColumn[j]);
        }
    }
    
    return;
}

template<typename CType, typename BMAccessorType>
SEXP BigMatrixToNifti(nifti_image *pnim, BigMatrix *pMat, SEXP indices) {
    switch (pnim->datatype) {
        case NIFTI_TYPE_UINT8:
            BigMatrixToNiftiData<CType, unsigned char, BMAccessorType>(pnim, pMat, indices);
            break;
        case NIFTI_TYPE_INT8:
            BigMatrixToNiftiData<CType, signed char, BMAccessorType>(pnim, pMat, indices);
            break;            
        case NIFTI_TYPE_UINT16:
            BigMatrixToNiftiData<CType, unsigned short, BMAccessorType>(pnim, pMat, indices);
            break;            
        case NIFTI_TYPE_INT16:
            BigMatrixToNiftiData<CType, short, BMAccessorType>(pnim, pMat, indices);
            break;            
        case NIFTI_TYPE_UINT32:
            BigMatrixToNiftiData<CType, unsigned int, BMAccessorType>(pnim, pMat, indices);
            break;            
        case NIFTI_TYPE_INT32:
            BigMatrixToNiftiData<CType, int, BMAccessorType>(pnim, pMat, indices);
            break;
        case NIFTI_TYPE_FLOAT32:
            BigMatrixToNiftiData<CType, float, BMAccessorType>(pnim, pMat, indices);
            break;
        case NIFTI_TYPE_FLOAT64:
            BigMatrixToNiftiData<CType, double, BMAccessorType>(pnim, pMat, indices);
            break;
    	default:
    	    warning("unsupported data format (identifier %d)", pnim->datatype);
            break;
    }
    
    return R_NilValue;
}

extern "C" {

//#include "International.h"

SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;
    
    int len_list = length(list);
    for (i = 0; i < len_list; i++)
    	if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
    	    elmt = VECTOR_ELT(list, i);
    	    break;
    	}
    return elmt;
}


/****************************************
 * Helper functions
 ****************************************/

SEXP getHeaderAttributes(void) {
    SEXP header;
    
    PROTECT(header=NEW_CHARACTER(NUM_HEADER_ATTRS));
    
    for(size_t i = 0; i < NUM_HEADER_ATTRS; ++i)
        SET_STRING_ELT(header, i, mkChar(header_attributes[i]));
    
    UNPROTECT(1);
    return header;
}

SEXP niftir_xform_string(SEXP Rxx) {
    int Cxx;
    SEXP_to_int(Rxx, &Cxx);
    char *Cstr = nifti_xform_string(Cxx);
    return pchar_to_SEXP(Cstr);
}

SEXP niftir_datatype_string(SEXP Rxx) {
    int Cxx;
    SEXP_to_int(Rxx, &Cxx);
    char *Cstr = nifti_datatype_string(Cxx);
    return pchar_to_SEXP(Cstr);
}
 

/****************************************
 * Niftilib wrapper functions
 ****************************************/

nifti_image *niftir_image_read(SEXP file, SEXP read_data) {
    nifti_image *pnim;
    
    PROTECT(read_data = AS_INTEGER(read_data));
    PROTECT(file = AS_CHARACTER(file));

    if(!isString(file) || length(file) != 1)
        error("niftir_image_read: file is not a single string\n");
    if(length(read_data) != 1)
        error("niftir_image_read: read_data is not a single integer\n");
    
    int *piread_data = INTEGER_POINTER(read_data);
    const char *pcfilename  = CHAR(STRING_ELT(file , 0));
    
    pnim = nifti_image_read(pcfilename, piread_data[0]) ;
    
    if (pnim==NULL) {
        error("Rnifti_image_read: Cannot open file \"%s\"",pcfilename);
        UNPROTECT(2);
    } else {
        UNPROTECT(2);
    }
    
    return pnim;
}

SEXP niftir_image_free(SEXP nim) {
    nifti_image *pnim = (nifti_image *)R_ExternalPtrAddr(nim);

    if(pnim!=NULL) {
        nifti_image_free(pnim);
        R_ClearExternalPtr(nim);
    }
    else {
        Rprintf("niftir_image_free: not a nifti pointer.\n");
    }
    
    return R_NilValue;
}

SEXP get_nifti_header(nifti_image *pnim) {
    SEXP header;
    
    PROTECT(header = NEW_LIST(NUM_HEADER_ATTRS));
    
    // toffset
    SET_ELEMENT(header, 0, float_to_SEXP(pnim->toffset));
    
    // descrip
    SET_ELEMENT(header, 1, pchar_to_SEXP(pnim->descrip));
    
    // fname
    SET_ELEMENT(header, 2, pchar_to_SEXP(pnim->fname));
    
    // iname
    SET_ELEMENT(header, 3, pchar_to_SEXP(pnim->iname));
    
    // slice duration
    SET_ELEMENT(header, 4, float_to_SEXP(pnim->slice_duration));
    
    // qform code
    SET_ELEMENT(header, 5, int_to_SEXP(pnim->qform_code));
    
    // sform code
    SET_ELEMENT(header, 6, int_to_SEXP(pnim->sform_code));
    
    // quatern_b
    SET_ELEMENT(header, 7, float_to_SEXP(pnim->quatern_b));
    
    // quatern_c
    SET_ELEMENT(header, 8, float_to_SEXP(pnim->quatern_c));
    
    // quatern_d
    SET_ELEMENT(header, 9, float_to_SEXP(pnim->quatern_d));
    
    // qoffset_x
    SET_ELEMENT(header, 10, float_to_SEXP(pnim->qoffset_x));
    
    // qoffset_y
    SET_ELEMENT(header, 11, float_to_SEXP(pnim->qoffset_y));
    
    // qoffset_z
    SET_ELEMENT(header, 12, float_to_SEXP(pnim->qoffset_z));
    
    // qfac
    SET_ELEMENT(header, 13, float_to_SEXP(pnim->qfac));
    
    // pixdim
    SEXP pixdim = R_NilValue;
    if (pnim->dim[0] > 0 && pnim->dim[0] < 8) {
        PROTECT(pixdim = NEW_NUMERIC(pnim->dim[0]));
        for (int i = 0; i < pnim->dim[0]; ++i)
            NUMERIC_POINTER(pixdim)[i] = pnim->pixdim[i+1];
        UNPROTECT(1);
    } else {
        error("number of dimensions (dim[0]) > 7");
    }
    SET_ELEMENT(header, 14, pixdim);
    
    // nifti type
    SET_ELEMENT(header, 15, int_to_SEXP(pnim->nifti_type));
    
    // size of header
    struct nifti_1_header hdr = nifti_convert_nim2nhdr(pnim);
    SET_ELEMENT(header, 16, int_to_SEXP(hdr.sizeof_hdr));
    
    // datatype
    SET_ELEMENT(header, 17, int_to_SEXP(pnim->datatype));
    
    // scl_slope: nifti1: Data scaling: slope.  analyze 7.5: float funused1
    SET_ELEMENT(header, 18, float_to_SEXP(pnim->scl_slope));
    
    // scl_inter nifti1: Data scaling: offset  analyze 7.5: float funused1
    SET_ELEMENT(header, 19, float_to_SEXP(pnim->scl_inter));
    
    // qto_xyz
    SET_ELEMENT(header, 20, mat44_to_SEXP(&(pnim->qto_xyz)));
    
    // qto_ij
    SET_ELEMENT(header, 21, mat44_to_SEXP(&(pnim->qto_ijk)));
    
    // sto_xyz
    SET_ELEMENT(header, 22, mat44_to_SEXP(&(pnim->sto_xyz)));
    
    // sto_ijk
    SET_ELEMENT(header, 23, mat44_to_SEXP(&(pnim->sto_ijk)));
    
    // dim
    SEXP dim = R_NilValue;
    if (pnim->dim[0] > 0 && pnim->dim[0] < 8) {
        PROTECT(dim = NEW_INTEGER(pnim->dim[0]));
        for (int i = 0; i < pnim->dim[0]; ++i)
            INTEGER_POINTER(dim)[i] = pnim->dim[i+1];
        UNPROTECT(1);
    } else {
        error("number of dimensions (dim[0]) > 7");
    }
    SET_ELEMENT(header, 24, dim);
    
    // xyz units
    SET_ELEMENT(header, 25, int_to_SEXP(pnim->xyz_units));
    
    // time units
    SET_ELEMENT(header, 26, int_to_SEXP(pnim->time_units));    
    
    // intent code
    SET_ELEMENT(header, 27, int_to_SEXP(pnim->intent_code));
    
    // intent name
    SET_ELEMENT(header, 28, pchar_to_SEXP(pnim->intent_name));
    
    // calibration param min
    SET_ELEMENT(header, 29, float_to_SEXP(pnim->cal_min));
    
    // calibration param max
    SET_ELEMENT(header, 30, float_to_SEXP(pnim->cal_max));
    
    // Set list element names
    SEXP names;
    PROTECT(names = allocVector(STRSXP, NUM_HEADER_ATTRS));
    for (int i = 0; header_attributes[i]!=NULL; ++i)
        SET_STRING_ELT(names, i, mkChar(header_attributes[i]));
    SET_NAMES(header, names);
    UNPROTECT(1);
            
    UNPROTECT(1);
    
    return header;
}

SEXP read_nifti(SEXP file, SEXP read_data) {
    // niftilib struct (temporary)
    nifti_image *pnim;
    
    // list to return to user
    SEXP nim;
    
    
    /***************
     * Get nifti image
     ***************/
    
     pnim = niftir_image_read(file, read_data);
    
    
    /******************
     * Get header info in list
     ******************/
     
     SEXP header = R_NilValue;
     PROTECT(header = get_nifti_header(pnim));
     
     if (pnim->data == NULL) {
         UNPROTECT(1);
         return header;
     }
    
    
    /******************
     * Create nifti R list (holds header and data)
     *****************/
     
     PROTECT(nim = NEW_LIST(2));
     SET_ELEMENT(nim, 0, header);
     
    
    /***************
     * Convert nifti image to R data structure
     ***************/
    
    // Allocate
    SEXP arr;
    PROTECT(arr = NEW_NUMERIC((double)pnim->nvox));
    
    // Convert/save nifti image to double
    switch (pnim->datatype) {
        case NIFTI_TYPE_UINT8:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                DOUBLE_DATA(arr)[i] = (double)((unsigned char*)pnim->data)[i];
            break;
        case NIFTI_TYPE_INT8:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                DOUBLE_DATA(arr)[i] = (double)((signed char*)pnim->data)[i];
            break;
        case NIFTI_TYPE_UINT16:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                DOUBLE_DATA(arr)[i] = (double)((unsigned short*)pnim->data)[i];
            break;
        case NIFTI_TYPE_INT16:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                DOUBLE_DATA(arr)[i] = (double)((short*)pnim->data)[i];
            break;
        case NIFTI_TYPE_UINT32:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                DOUBLE_DATA(arr)[i] = (double)((unsigned int*)pnim->data)[i];
            break;
        case NIFTI_TYPE_INT32:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                DOUBLE_DATA(arr)[i] = (double)((int*)pnim->data)[i];
            break;
        case NIFTI_TYPE_FLOAT32:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                DOUBLE_DATA(arr)[i] = (double)((float*)pnim->data)[i];
            break;
        case NIFTI_TYPE_FLOAT64:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                DOUBLE_DATA(arr)[i] = (double)((double*)pnim->data)[i];
            break;
    	default:
    	    warning("unsupported data format (identifier %d)", pnim->datatype);
            break;
    }
    
    // set image dimensions (based on dim from header)
    SET_DIM(arr, GET_ELEMENT(header, 24));
    
    // Add data to list
    SET_ELEMENT(nim, 1, arr);
    
    // Give names to list
    SEXP names;
    PROTECT(names = allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("header"));
    SET_STRING_ELT(names, 1, mkChar("image"));
    SET_NAMES(nim, names);
    UNPROTECT(1);
    
    // Remove niftilib object
    nifti_image_free(pnim);
    
    // Unprotect nim and arr
    UNPROTECT(3);
    
    return nim;
}

nifti_image *create_nifti_image(SEXP header, SEXP Rdim, SEXP Rdatatype, SEXP outfile) {
    nifti_image *pnim;
    int dim[8], datatype;
    
    if (!IS_LIST(header))
        error("header must be a list");
    
    // Get dimensions from data
    int lendim = length(Rdim);
    if (lendim>1 && lendim<8) {
        dim[0] = lendim;
        for(int i = 0; i < lendim; ++i)
            dim[i+1] = INTEGER_POINTER(Rdim)[i];
        for(int i = lendim; i < 8; ++i)
            dim[i+1] = 0;
    } else {
        error("number of dimensions %d not compatible", lendim);
    }
    UNPROTECT(1);
    
    // Get datatype from header or based on data
    SEXP_to_int(Rdatatype, &datatype);
            
    // Create niftilib object
    pnim = nifti_make_new_nim(dim, datatype, 1);
    
    // Set output file
    PROTECT(outfile = AS_CHARACTER(outfile));
    const char *pcoutfile = CHAR(STRING_ELT(outfile, 0));
    int fname_result = nifti_set_filenames(pnim, pcoutfile, 1, 1);
    UNPROTECT(1);
    if (fname_result != 0)
        error("file already exists, can't write output\n");
    
    set_nifti_header(pnim, header);
    
    return pnim;
}

void set_nifti_header(nifti_image *pnim, SEXP header) {
    /* Go through each element in header and copy over to pnim */
    SEXP tmpRElement;
    
    // toffset
    PROTECT(tmpRElement = GET_ELEMENT(header, 0));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->toffset));
    UNPROTECT(1);
    
    // descrip
    PROTECT(tmpRElement = GET_ELEMENT(header, 1));
    if (tmpRElement != R_NilValue)
        SEXP_to_pchar(tmpRElement, pnim->descrip, 80);
    UNPROTECT(1);
        
    //header and image names (fname, iname)
    // already done
    
    // slice duration
    PROTECT(tmpRElement = GET_ELEMENT(header, 4));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->slice_duration));
    UNPROTECT(1);
    
    // qform code
    PROTECT(tmpRElement = GET_ELEMENT(header, 5));
    if (tmpRElement != R_NilValue)
        SEXP_to_int(tmpRElement, &(pnim->qform_code));
    UNPROTECT(1);
    
    // sform code
    PROTECT(tmpRElement = GET_ELEMENT(header, 6));
    if (tmpRElement != R_NilValue)
        SEXP_to_int(tmpRElement, &(pnim->sform_code));
    UNPROTECT(1);
    
    // quatern_b
    PROTECT(tmpRElement = GET_ELEMENT(header, 7));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->quatern_b));
    UNPROTECT(1);
    
    // quatern_c
    PROTECT(tmpRElement = GET_ELEMENT(header, 8));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->quatern_c));
    UNPROTECT(1);
    
    // quatern_d
    PROTECT(tmpRElement = GET_ELEMENT(header, 9));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->quatern_d));
    UNPROTECT(1);
    
    // qoffset_x
    PROTECT(tmpRElement = GET_ELEMENT(header, 10));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->qoffset_x));
    UNPROTECT(1);
    
    // qoffset_y
    PROTECT(tmpRElement = GET_ELEMENT(header, 11));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->qoffset_y));
    UNPROTECT(1);
    
    // qoffset_z
    PROTECT(tmpRElement = GET_ELEMENT(header, 12));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->qoffset_z));
    UNPROTECT(1);
    
    // qfac
    PROTECT(tmpRElement = GET_ELEMENT(header, 13));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->qfac));
    UNPROTECT(1);
    
    // pixdim
    PROTECT(tmpRElement = GET_ELEMENT(header, 14));
    pnim->pixdim[0] = pnim->ndim;
    if (tmpRElement != R_NilValue) {
        tmpRElement = AS_NUMERIC(tmpRElement);
        for(int i = 0; i < pnim->dim[0]; ++i) {
            pnim->pixdim[i+1] = (float)(NUMERIC_POINTER(tmpRElement)[i]);
        }
    }
    UNPROTECT(1);
    
    // nifti type
    PROTECT(tmpRElement = GET_ELEMENT(header, 15));
    if (tmpRElement != R_NilValue)
        SEXP_to_int(tmpRElement, &(pnim->nifti_type));
    else
        pnim->nifti_type = NIRTS_FTYPE;
    UNPROTECT(1);
    
    // size of header
    // nothing for 16
    
    // datatype 17
    //pnim->datatype = datatype;
    
    // scl_slope: nifti1: Data scaling: slope.  analyze 7.5: float funused1
    PROTECT(tmpRElement = GET_ELEMENT(header, 18));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->scl_slope));
    UNPROTECT(1);
    
    // scl_inter nifti1: Data scaling: offset  analyze 7.5: float funused1
    PROTECT(tmpRElement = GET_ELEMENT(header, 19));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->scl_inter));
    UNPROTECT(1);
    
    // qto_xyz
    PROTECT(tmpRElement = GET_ELEMENT(header, 20));
    if (tmpRElement != R_NilValue)
        SEXP_to_mat44(tmpRElement, &(pnim->qto_xyz));
    UNPROTECT(1);
    
    // qto_ij
    PROTECT(tmpRElement = GET_ELEMENT(header, 21));
    if (tmpRElement != R_NilValue)
        SEXP_to_mat44(tmpRElement, &(pnim->qto_ijk));
    UNPROTECT(1);
    
    // sto_xyz
    PROTECT(tmpRElement = GET_ELEMENT(header, 22));
    if (tmpRElement != R_NilValue)
        SEXP_to_mat44(tmpRElement, &(pnim->sto_xyz));
    UNPROTECT(1);
    
    // sto_ijk
    PROTECT(tmpRElement = GET_ELEMENT(header, 23));
    if (tmpRElement != R_NilValue)
        SEXP_to_mat44(tmpRElement, &(pnim->sto_ijk));
    UNPROTECT(1);
    
    // dim
    // already set
    
    // xyz unit
    PROTECT(tmpRElement = GET_ELEMENT(header, 25));
    if (tmpRElement != R_NilValue)
        SEXP_to_int(tmpRElement, &(pnim->xyz_units));
    UNPROTECT(1);
    
    // time unit
    PROTECT(tmpRElement = GET_ELEMENT(header, 26));
    if (tmpRElement != R_NilValue)
        SEXP_to_int(tmpRElement, &(pnim->time_units));
    UNPROTECT(1);
    
    // intent code
    PROTECT(tmpRElement = GET_ELEMENT(header, 27));
    if (tmpRElement != R_NilValue)
        SEXP_to_int(tmpRElement, &(pnim->intent_code));
    UNPROTECT(1);
    
    // intent name
    PROTECT(tmpRElement = GET_ELEMENT(header, 28));
    if (tmpRElement != R_NilValue)
        SEXP_to_pchar(tmpRElement, pnim->intent_name, 16);
    UNPROTECT(1);
    
    // calibration param min
    PROTECT(tmpRElement = GET_ELEMENT(header, 29));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->cal_min));
    UNPROTECT(1);
    
    // calibration param max
    PROTECT(tmpRElement = GET_ELEMENT(header, 30));
    if (tmpRElement != R_NilValue)
        SEXP_to_float(tmpRElement, &(pnim->cal_max));
    UNPROTECT(1);
    
    nifti_update_dims_from_array(pnim);
    
    return;
}

SEXP write_nifti(SEXP header, SEXP data, SEXP outfile) {    
    SEXP Rdim, Rdatatype;
    
    if (!IS_VECTOR(data))
        error("data must be vector type");
    
    // Get dim
    PROTECT(Rdim = GET_DIM(data));
    
    // Get datatype
    PROTECT(Rdatatype = GET_LIST_ELEMENT(header, "datatype"));
    if (Rdatatype == R_NilValue) {
        int datatype;
        switch (TYPEOF(data)) {
            case REALSXP:
                datatype = DT_FLOAT32;
                break;
            case INTSXP:
                datatype = DT_INT32;
                break;
            case LGLSXP:
                datatype = DT_BINARY;
                break;
            default:
                datatype = DT_FLOAT32;
                break;
        }
        Rdatatype = int_to_SEXP(datatype);
    }
    UNPROTECT(1);
    
    nifti_image *pnim = create_nifti_image(header, Rdim, Rdatatype, outfile);
    
    /* Save R data to pnim */
    // Convert/save nifti image to double
    switch (pnim->datatype) {
        case NIFTI_TYPE_UINT8:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                ((unsigned char*)pnim->data)[i] = (unsigned char)(DOUBLE_DATA(data)[i]);
            break;
        case NIFTI_TYPE_INT8:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                ((signed char*)pnim->data)[i] = (signed char)(DOUBLE_DATA(data)[i]);
            break;
        case NIFTI_TYPE_UINT16:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                ((unsigned short*)pnim->data)[i] = (unsigned short)(DOUBLE_DATA(data)[i]);
            break;
        case NIFTI_TYPE_INT16:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                ((short*)pnim->data)[i] = (short)(DOUBLE_DATA(data)[i]);
            break;
        case NIFTI_TYPE_UINT32:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                ((unsigned int*)pnim->data)[i] = (unsigned int)(DOUBLE_DATA(data)[i]);
            break;
        case NIFTI_TYPE_INT32:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                ((int*)pnim->data)[i] = (int)(DOUBLE_DATA(data)[i]);
            break;
        case NIFTI_TYPE_FLOAT32:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                ((float*)pnim->data)[i] = (float)(DOUBLE_DATA(data)[i]);
            break;
        case NIFTI_TYPE_FLOAT64:
            for(size_t i = 0 ; i < pnim->nvox; ++i)
                ((double*)pnim->data)[i] = (double)(DOUBLE_DATA(data)[i]);
            break;
    	default:
    	    warning("unsupported data format (identifier %d)", pnim->datatype);
            break;
    }
    
    if (!nifti_nim_is_valid(pnim, 1))
        error("data seems invalid");
    
    if (pnim!=NULL)
        nifti_image_write(pnim);
    else
        error("pnim was NULL");
    
    nifti_image_free(pnim);
    
    return R_NilValue;
}


/****************************************
 * big nifti functions
 ****************************************/

SEXP read_bignifti_header(SEXP file) {
    SEXP nim, header, ptr;
    nifti_image *pnim;

    // Read in file
    pnim = niftir_image_read(file, int_to_SEXP(0));
    
    // Convert header to list
    PROTECT(header = get_nifti_header(pnim));

    // Get nifti pointer    
    ptr = R_MakeExternalPtr(pnim, install("NIFTI_TYPE_TAG"), R_NilValue);
    PROTECT(ptr);
    R_RegisterCFinalizer(ptr, (R_CFinalizer_t) niftir_image_free);

    // Create nifti R list (holds header and nifti address)
    PROTECT(nim = NEW_LIST(2));
    SET_ELEMENT(nim, 0, header);
    SET_ELEMENT(nim, 1, ptr);

    // Give names to list
    SEXP names;
    PROTECT(names = allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("header"));
    SET_STRING_ELT(names, 1, mkChar("address"));
    SET_NAMES(nim, names);
    
    UNPROTECT(4);
    return nim;
}

SEXP read_bignifti_data(SEXP nim_addr, SEXP big_addr) {
    // Get nifti object pointer
    nifti_image *pnim = (nifti_image *)R_ExternalPtrAddr(nim_addr);
    
    // Get nifti data
    if (nifti_image_load(pnim) < 0) {
       nifti_image_free(pnim);
       error("Could not load nifti data");
    }
    
    // Save data to big matrix object
    BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(big_addr));    
    
    if (pMat->separated_columns()) {
        switch (pMat->matrix_type()) {
            case 1:
                return NiftiToBigMatrix<char>(pnim, pMat, SepMatrixAccessor<char>(*pMat));
                break;
            case 2:
                return NiftiToBigMatrix<short>(pnim, pMat, SepMatrixAccessor<short>(*pMat));
                break;
            case 4:
                return NiftiToBigMatrix<int>(pnim, pMat, SepMatrixAccessor<int>(*pMat));
                break;
            case 8:
                return NiftiToBigMatrix<double>(pnim, pMat, SepMatrixAccessor<double>(*pMat));
                break;
        }
    }
    else {
        switch (pMat->matrix_type()) {
            case 1:
                return NiftiToBigMatrix<char>(pnim, pMat, MatrixAccessor<char>(*pMat));
                break;
            case 2:
                return NiftiToBigMatrix<short>(pnim, pMat, MatrixAccessor<short>(*pMat));
                break;
            case 4:
                return NiftiToBigMatrix<int>(pnim, pMat, MatrixAccessor<int>(*pMat));
                break;
            case 8:
                return NiftiToBigMatrix<double>(pnim, pMat, MatrixAccessor<double>(*pMat));
                break;
        }
    }
    
    error("failed to identify big matrix type");
}

SEXP write_bignifti(SEXP header, SEXP big_addr, SEXP indices, SEXP outfile) {    
    SEXP Rdim, Rdatatype;
    
    // Load big matrix
    BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(big_addr));    
    
    // Get dim
    PROTECT(Rdim = GET_LIST_ELEMENT(header, "dim"));
    if (Rdim == R_NilValue)
        error("header must have a proper dim (dimension) attribute");
    if (GET_LENGTH(Rdim) != 4)
        error("header must have a 4D dim (dimension) attribute");
    
    // Get datatype
    PROTECT(Rdatatype = GET_LIST_ELEMENT(header, "datatype"));
    if (Rdatatype == R_NilValue) {
        int datatype;
        switch(pMat->matrix_type()) {
            case 1:     // char
                datatype = DT_INT8;
                break;
            case 2:     // short
                datatype = DT_INT16;
                break;
            case 4:     // int
                datatype = DT_INT32;
                break;
            case 8:     // double
                datatype = DT_FLOAT64;
                break;
            default:
                error("unrecognized big matrix data type");
                break;
        }
        Rdatatype = int_to_SEXP(datatype);
    }
    UNPROTECT(1);
    
    // Load nifti object
    nifti_image *pnim = create_nifti_image(header, Rdim, Rdatatype, outfile);
    
    if (pMat->separated_columns()) {
        switch (pMat->matrix_type()) {
            case 1:
                BigMatrixToNifti<char, SepMatrixAccessor<char> >(pnim, pMat, indices);
                break;
            case 2:
                BigMatrixToNifti<short, SepMatrixAccessor<short> >(pnim, pMat, indices);
                break;
            case 4:
                BigMatrixToNifti<int, SepMatrixAccessor<int> >(pnim, pMat, indices);
                break;
            case 8:
                BigMatrixToNifti<double, SepMatrixAccessor<double> >(pnim, pMat, indices);
                break;
        }
    }
    else {
        switch (pMat->matrix_type()) {
            case 1:
                BigMatrixToNifti<char, MatrixAccessor<char> >(pnim, pMat, indices);
                break;
            case 2:
                BigMatrixToNifti<short, MatrixAccessor<short> >(pnim, pMat, indices);
                break;
            case 4:
                BigMatrixToNifti<int, MatrixAccessor<int> >(pnim, pMat, indices);
                break;
            case 8:
                BigMatrixToNifti<double, MatrixAccessor<double> >(pnim, pMat, indices);
                break;
        }
    }
    
    if (!nifti_nim_is_valid(pnim, 1))
        error("data seems invalid");
    
    if (pnim!=NULL)
        nifti_image_write(pnim);
    else
        error("pnim was NULL");
        
    nifti_image_free(pnim);
    
    return R_NilValue;
}

} // End Extern C

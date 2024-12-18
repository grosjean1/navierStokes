//------------------------------------------------------------------------------
// GB_mex_Matrix_extractElement: MATLAB interface for x = A(i,j)
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// x = A (i,j), where i and j are zero-based.  If i and j arrays, then
// x (k) = A (i (k), j (k)) is done for all k.

// I and J and zero-based

#include "GB_mex.h"

#define FREE_ALL                        \
{                                       \
    GB_MATRIX_FREE (&A) ;               \
    GB_FREE_MEMORY (Xtemp) ;            \
    GB_mx_put_global (malloc_debug) ;   \
}

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{

    bool malloc_debug = GB_mx_get_global ( ) ;
    GrB_Matrix A = NULL ;
    void *Y = NULL ;
    void *Xtemp = NULL ;

    // check inputs
    if (nargout > 1 || nargin < 3 || nargin > 4)
    {
        mexErrMsgTxt ("Usage: x = GB_mex_Matrix_extractElement "
            "(A, I, J, xclass)");
    }

    #define GET_DEEP_COPY ;
    #define FREE_DEEP_COPY ;

    // get A (shallow copy)
    A = GB_mx_mxArray_to_Matrix (pargin [0], "A input", false) ;
    if (A == NULL)
    {
        FREE_ALL ;
        mexErrMsgTxt ("A failed") ;
    }
    mxClassID aclass = GB_mx_Type_to_classID (A->type) ;

    // get I
    GrB_Index *I, ni ; 
    if (!GB_mx_mxArray_to_indices (&I, pargin [1], &ni))
    {
        FREE_ALL ;
        mexErrMsgTxt ("I failed") ;
    }

    // get J
    GrB_Index *J, nj ; 
    if (!GB_mx_mxArray_to_indices (&J, pargin [2], &nj))
    {
        FREE_ALL ;
        mexErrMsgTxt ("J failed") ;
    }

    if (ni != nj)
    {
        FREE_ALL ;
        mexErrMsgTxt ("I and J must be the same size") ;
    }

    // get xclass, default is class (A), and the corresponding xtype
    mxClassID xclass ;
    GrB_Type xtype ;

    if (A->type == Complex)
    {
        // input argument xclass is ignored
        xtype = Complex ;
        xclass = mxDOUBLE_CLASS ;
        // create Xtemp
        GB_CALLOC_MEMORY (Xtemp, ni, sizeof (double complex)) ;
    }
    else
    {
        xclass = GB_mx_string_to_classID (aclass, PARGIN (3)) ;
        xtype = GB_mx_classID_to_Type (xclass) ;
        if (xtype == NULL)
        {
            FREE_ALL ;
            mexErrMsgTxt ("X must be numeric") ;
        }
        // create Y
        pargout [0] = mxCreateNumericMatrix (ni, 1, xclass, mxREAL) ;
        Y = mxGetData (pargout [0]) ;
    }

    size_t s = sizeof (double complex) ;

    // x = A (i,j)
    switch (xtype->code)
    {
        case GB_BOOL_code   :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                bool *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break ;

        case GB_INT8_code   :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                int8_t *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break ;

        case GB_UINT8_code  :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                uint8_t *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break ;

        case GB_INT16_code  :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                int16_t *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break ;

        case GB_UINT16_code :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                uint16_t *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break ;

        case GB_INT32_code  :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                int32_t *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break ;

        case GB_UINT32_code :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                uint32_t *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break ;

        case GB_INT64_code  :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                int64_t *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break ;

        case GB_UINT64_code :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                uint64_t *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break ;

        case GB_FP32_code   :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                float *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break ;

        case GB_FP64_code   :

            for (int64_t k = 0 ; k < ni ; k++)
            {
                double *X = Y ;
                METHOD (GrB_Matrix_extractElement (&X [k], A, I [k], J [k])) ;
            }
            break;

        case GB_UDT_code   :
            {
                // user-defined complex type
                for (int64_t k = 0 ; k < ni ; k++)
                {
                    METHOD (GrB_Matrix_extractElement (Xtemp +(k*s),
                        A, I [k], J [k])) ;
                }
            }
            break;

        default              :
        
            FREE_ALL ;
            mexErrMsgTxt ("unsupported class") ;
    }

    if (A->type == Complex)
    {
        // create the MATLAB complex X
        pargout [0] = mxCreateNumericMatrix (ni, 1, mxDOUBLE_CLASS, mxCOMPLEX) ;
        GB_mx_complex_split (ni, Xtemp, pargout [0]) ;
    }

    FREE_ALL ;

}


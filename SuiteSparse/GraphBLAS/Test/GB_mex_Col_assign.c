//------------------------------------------------------------------------------
// GB_mex_assign: C<Mask>(I,J) = accum (C (I,J), A)
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

// This function is a wrapper for GrB_Matrix_assign, GrB_Matrix_assign_T
// GrB_Vector_assign, and GrB_Vector_assign_T.  For these uses, the Mask must
// always be the same size as C.

// This mexFunction does not call GrB_Row_assign or GrB_Col_assign, since
// the Mask is a single row or column in these cases, and C is not modified
// outside that single row (for GrB_Row_assign) or column (for GrB_Col_assign).

// This function does the same thing as the MATLAB mimic GB_spec_assign.m.

//------------------------------------------------------------------------------

#include "GB_mex.h"

#define FREE_ALL                        \
{                                       \
    GB_MATRIX_FREE (&A) ;               \
    GB_MATRIX_FREE (&Mask) ;            \
    GB_MATRIX_FREE (&C) ;               \
    GrB_free (&desc) ;                  \
    GB_mx_put_global (malloc_debug) ;   \
}

#define GET_DEEP_COPY \
    C = GB_mx_mxArray_to_Matrix (pargin [0], "C input", true) ;

#define FREE_DEEP_COPY GB_MATRIX_FREE (&C) ;

GrB_Matrix C = NULL ;
GrB_Matrix Mask = NULL ;
GrB_Matrix A = NULL ;
GrB_Descriptor desc = NULL ;
GrB_BinaryOp accum = NULL ;
GrB_Index *I = NULL, ni ;
GrB_Index *J = NULL, nj ;
bool malloc_debug = false ;
GrB_Info info = GrB_SUCCESS ;

//------------------------------------------------------------------------------
// assign: perform a single assignment
//------------------------------------------------------------------------------

#define OK(method)                      \
{                                       \
    info = method ;                     \
    if (info != GrB_SUCCESS)            \
    {                                   \
        GB_MATRIX_FREE (&mask) ;        \
        GB_MATRIX_FREE (&u) ;           \
        return (info) ;                 \
    }                                   \
}

GrB_Info assign ( )
{
    bool at = (desc != NULL && desc->in0 == GrB_TRAN) ;
    GrB_Info info ;
    GrB_Matrix mask = NULL, u = NULL ;

    // printf ("\n--- assign:\n") ;
    ASSERT_OK (GB_check (C, "C", 0)) ;
    ASSERT_OK_OR_NULL (GB_check (Mask, "Mask", 0)) ;
    ASSERT_OK (GB_check (A, "A", 0)) ;
    ASSERT_OK_OR_NULL (GB_check (accum, "accum", 0)) ;
    ASSERT_OK_OR_NULL (GB_check (desc, "desc", 0)) ;

    /*
    if (I == NULL)
    {
        printf ("I is NULL\n") ;
    }
    else if (I == GrB_ALL)
    {
        printf ("I is ALL\n") ;
    }
    else
    {
        for (int64_t k = 0 ; k < ni ; k++) printf ("I [%lld] = %lld\n", k, I [k]) ;
    }
    if (J == NULL)
    {
        printf ("J is NULL\n") ;
    }
    else if (J == GrB_ALL)
    {
        printf ("J is ALL\n") ;
    }
    else
    {
        for (int64_t k = 0 ; k < nj ; k++) printf ("J [%lld] = %lld\n", k, J [k]) ;
    }
    */

    if (A->nrows == 1 && A->ncols == 1 && NNZ (A) == 1)
    {
        // scalar expansion to matrix or vector
        void *Ax = A->x ;

        if (ni == 1 && nj == 1 && Mask == NULL && I != GrB_ALL && J != GrB_ALL
            && GB_op_is_second (accum, C->type) && A->type->code != GB_UDT_code
            && desc == NULL)
        {
            // printf ("setElement\n") ;
            // test GrB_Matrix_setElement
            #define ASSIGN(type)                                        \
            {                                                           \
                type x = ((type *) Ax) [0] ;                            \
                OK (GrB_Matrix_setElement (C, x, I [0], J [0])) ;       \
            } break ;

            switch (A->type->code)
            {
                case GB_BOOL_code   : ASSIGN (bool) ;
                case GB_INT8_code   : ASSIGN (int8_t) ;
                case GB_UINT8_code  : ASSIGN (uint8_t) ;
                case GB_INT16_code  : ASSIGN (int16_t) ;
                case GB_UINT16_code : ASSIGN (uint16_t) ;
                case GB_INT32_code  : ASSIGN (int32_t) ;
                case GB_UINT32_code : ASSIGN (uint32_t) ;
                case GB_INT64_code  : ASSIGN (int64_t) ;
                case GB_UINT64_code : ASSIGN (uint64_t) ;
                case GB_FP32_code   : ASSIGN (float) ;
                case GB_FP64_code   : ASSIGN (double) ;
                case GB_UDT_code    :
                default:
                    FREE_ALL ;
                    mexErrMsgTxt ("unsupported class") ;
            }
            #undef ASSIGN

            ASSERT_OK (GB_check (C, "C after setElement", 0)) ;

        }
        if (C->ncols == 1)
        {

            // test GrB_Vector_assign_scalar functions
            // printf ("scalar assign to vector\n") ;
            #define ASSIGN(type)                                        \
            {                                                           \
                type x = ((type *) Ax) [0] ;                            \
                OK (GrB_assign ((GrB_Vector) C, (GrB_Vector) Mask,      \
                    accum, x, I, ni, desc)) ;      \
            } break ;

            switch (A->type->code)
            {
                case GB_BOOL_code   : ASSIGN (bool) ;
                case GB_INT8_code   : ASSIGN (int8_t) ;
                case GB_UINT8_code  : ASSIGN (uint8_t) ;
                case GB_INT16_code  : ASSIGN (int16_t) ;
                case GB_UINT16_code : ASSIGN (uint16_t) ;
                case GB_INT32_code  : ASSIGN (int32_t) ;
                case GB_UINT32_code : ASSIGN (uint32_t) ;
                case GB_INT64_code  : ASSIGN (int64_t) ;
                case GB_UINT64_code : ASSIGN (uint64_t) ;
                case GB_FP32_code   : ASSIGN (float) ;
                case GB_FP64_code   : ASSIGN (double) ;
                case GB_UDT_code    :
                {
                    OK (GrB_assign ((GrB_Vector) C, (GrB_Vector) Mask,
                        accum, Ax, I, ni, desc)) ;
                }
                break ;
                default:
                    FREE_ALL ;
                    mexErrMsgTxt ("unsupported class") ;
            }
            #undef ASSIGN

        }
        else
        {

            // test Matrix_assign_scalar functions
            // printf ("scalar assign to matrix\n") ;
            #define ASSIGN(type)                                            \
            {                                                               \
                type x = ((type *) Ax) [0] ;                                \
                OK (GrB_assign (C, Mask, accum, x, I, ni, J, nj,desc)) ;    \
            } break ;

            switch (A->type->code)
            {
                case GB_BOOL_code   : ASSIGN (bool) ;
                case GB_INT8_code   : ASSIGN (int8_t) ;
                case GB_UINT8_code  : ASSIGN (uint8_t) ;
                case GB_INT16_code  : ASSIGN (int16_t) ;
                case GB_UINT16_code : ASSIGN (uint16_t) ;
                case GB_INT32_code  : ASSIGN (int32_t) ;
                case GB_UINT32_code : ASSIGN (uint32_t) ;
                case GB_INT64_code  : ASSIGN (int64_t) ;
                case GB_UINT64_code : ASSIGN (uint64_t) ;
                case GB_FP32_code   : ASSIGN (float) ;
                case GB_FP64_code   : ASSIGN (double) ;
                case GB_UDT_code    :
                {
                    OK (GrB_assign (C, Mask, accum, Ax, I, ni, J, nj, desc)) ;
                }
                break ;

                default:
                    FREE_ALL ;
                    mexErrMsgTxt ("unsupported class") ;
            }
            #undef ASSIGN

        }

    }
    else if (C->ncols == 1 && A->ncols == 1 &&
        (Mask == NULL || Mask->ncols == 1) && !at)
    {
        // test GrB_Vector_assign
        // printf ("vector assign\n") ;
        OK (GrB_assign ((GrB_Vector) C, (GrB_Vector) Mask, accum,
            (GrB_Vector) A, I, ni, desc)) ;
    }
    else
    {
        // standard submatrix assignment
        // printf ("submatrix assign\n") ;
        OK (GrB_assign (C, Mask, accum, A, I, ni, J, nj, desc)) ;
    }
    return (info) ;
}

//------------------------------------------------------------------------------
// many_assign: do a sequence of assignments
//------------------------------------------------------------------------------

// The list of assignments is in a struct array

GrB_Info many_assign
(
    int nwork,
    int fA,
    int fI,
    int fJ,
    int faccum,
    int fMask,
    int fdesc,
    mxClassID cclass,
    const mxArray *pargin [ ]
)
{

#define CATCH 9999999999

    GrB_Info info = GrB_SUCCESS ;

    for (int64_t k = 0 ; k < nwork ; k++)
    {

        //----------------------------------------------------------------------
        // get the kth work to do
        //----------------------------------------------------------------------

        // each struct has fields A, I, J, and optionally Mask, accum, and desc

        mxArray *p ;

        // if (k == CATCH) GB_check (C, "C start", 3) ;

        // [ turn off malloc debugging
        bool save = GB_thread_local.malloc_debug ;
        GB_thread_local.malloc_debug = false ;

        // get Mask (shallow copy)
        Mask = NULL ;
        if (fMask >= 0)
        {
            p = mxGetFieldByNumber (pargin [1], k, fMask) ;
            Mask = GB_mx_mxArray_to_Matrix (p, "Mask", false) ;
            if (Mask == NULL && !mxIsEmpty (p))
            {
                FREE_ALL ;
                mexErrMsgTxt ("Mask failed") ;
            }
        }
        // if (k == CATCH) GB_check (Mask, "Mask", 3) ;

        // get A (shallow copy)
        p = mxGetFieldByNumber (pargin [1], k, fA) ;
        A = GB_mx_mxArray_to_Matrix (p, "A", false) ;
        if (A == NULL)
        {
            FREE_ALL ;
            mexErrMsgTxt ("A failed") ;
        }
        // if (k == CATCH) GB_check (A, "A", 3) ;

        // get accum; default: NOP, default class is class(C)
        accum = NULL ;
        if (faccum >= 0)
        {
            p = mxGetFieldByNumber (pargin [1], k, faccum) ;
            if (!GB_mx_mxArray_to_BinaryOp (&accum, p, "accum",
                GB_NOP_opcode, cclass,
                C->type == Complex, A->type == Complex))
            {
                FREE_ALL ;
                mexErrMsgTxt ("accum failed") ;
            }
        }
        // if (k == CATCH) GB_check (accum, "accum", 3) ;

        // get I
        p = mxGetFieldByNumber (pargin [1], k, fI) ;
        if (!GB_mx_mxArray_to_indices (&I, p, &ni))
        {
            FREE_ALL ;
            mexErrMsgTxt ("I failed") ;
        }

        /*
        if (k == CATCH)
        {
            printf ("I: ") ;
            if (I == GrB_ALL) printf ("all") ;
            else for (int kk = 0 ; kk < ni ; kk++) printf ("%lld ", I [kk]) ;
            printf ("\n") ;
        }
        */

        // get J
        p = mxGetFieldByNumber (pargin [1], k, fJ) ;
        if (!GB_mx_mxArray_to_indices (&J, p, &nj))
        {
            FREE_ALL ;
            mexErrMsgTxt ("J failed") ;
        }

        /*
        if (k == CATCH)
        {
            printf ("J: ") ;
            if (J == GrB_ALL) printf ("all") ;
            else for (int kk = 0 ; kk < ni ; kk++) printf ("%lld ", J [kk]) ;
            printf ("\n") ;
        }
        */

        /*
        printf ("many assign: fI %d fJ %d ni %lld nj %lld\n", fI, fJ, ni, nj) ;
        for (int kk = 0 ; kk < nj ; kk++) printf ("J[%d]=%lld\n", kk, J[kk]) ;
        */

        // get desc
        desc = NULL ;
        if (fdesc > 0)
        {
            p = mxGetFieldByNumber (pargin [1], k, fdesc) ;
            if (!GB_mx_mxArray_to_Descriptor (&desc, p, "desc"))
            {
                FREE_ALL ;
                mexErrMsgTxt ("desc failed") ;
            }
        }
        // restore malloc debugging to test the method
        GB_thread_local.malloc_debug = save ;   // ]

        // GB_check (desc, "desc", 3) ;

        //----------------------------------------------------------------------
        // C<Mask>(I,J) = A
        //----------------------------------------------------------------------

        // GB_thread_local.line = (k == CATCH) ? -911 : 0 ;

        info = assign ( ) ;

        // GB_thread_local.line = 0 ;

        // if (k == CATCH) GB_check (C, "C done", 3) ;

        GB_MATRIX_FREE (&A) ;
        GB_MATRIX_FREE (&Mask) ;
        GrB_free (&desc) ;

        if (info != GrB_SUCCESS)
        {
            return (info) ;
        }
    }
    return (info) ;
}

//------------------------------------------------------------------------------
// GB_mex_assign mexFunction
//------------------------------------------------------------------------------

void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    malloc_debug = GB_mx_get_global ( ) ;
    A = NULL ;
    C = NULL ;
    Mask = NULL ;
    desc = NULL ;

    // printf ("\n========================= GB_mex_assign:\n") ;

    if (nargout > 1 || ! (nargin == 2 || nargin == 6 || nargin == 7))
    {
        mexErrMsgTxt ("Usage: C = GB_mex_assign "
       "(C, Mask, accum, A, I, J, desc) or (C, Work)");
    }

    //--------------------------------------------------------------------------
    // get C (make a deep copy)
    //--------------------------------------------------------------------------

    GET_DEEP_COPY ;
    if (C == NULL)
    {
        FREE_ALL ;
        mexErrMsgTxt ("C failed") ;
    }
    mxClassID cclass = GB_mx_Type_to_classID (C->type) ;

    if (nargin == 2)
    {

        //----------------------------------------------------------------------
        // get a list of work to do: a struct array of length nwork
        //----------------------------------------------------------------------

        // each entry is a struct with fields:
        // Mask, accum, A, I, J, desc

        if (!mxIsStruct (pargin [1]))
        {
            FREE_ALL ;
            mexErrMsgTxt ("2nd argument must be a struct") ;
        }

        int nwork = mxGetNumberOfElements (pargin [1]) ;
        int nf = mxGetNumberOfFields (pargin [1]) ;
        for (int f = 0 ; f < nf ; f++)
        {
            mxArray *p ;
            for (int k = 0 ; k < nwork ; k++)
            {
                p = mxGetFieldByNumber (pargin [1], k, f) ;
            }
        }

        int fA = mxGetFieldNumber (pargin [1], "A") ;
        int fI = mxGetFieldNumber (pargin [1], "I") ;
        int fJ = mxGetFieldNumber (pargin [1], "J") ;
        int faccum = mxGetFieldNumber (pargin [1], "accum") ;
        int fMask = mxGetFieldNumber (pargin [1], "Mask") ;
        int fdesc = mxGetFieldNumber (pargin [1], "desc") ;

        if (fA < 0 || fI < 0 || fJ < 0) mexErrMsgTxt ("A,I,J required") ;

        // GB_thread_local.line = 0 ;

        METHOD (many_assign (nwork, fA, fI, fJ, faccum, fMask, fdesc, cclass,
            pargin)) ;

    }
    else
    {

        //----------------------------------------------------------------------
        // C<Mask>(I,J) = A, with a single assignment
        //----------------------------------------------------------------------

        // get Mask (shallow copy)
        Mask = GB_mx_mxArray_to_Matrix (pargin [1], "Mask", false) ;
        if (Mask == NULL && !mxIsEmpty (pargin [1]))
        {
            FREE_ALL ;
            mexErrMsgTxt ("Mask failed") ;
        }

        // get A (shallow copy)
        A = GB_mx_mxArray_to_Matrix (pargin [3], "A", false) ;
        if (A == NULL)
        {
            FREE_ALL ;
            mexErrMsgTxt ("A failed") ;
        }

        // get accum; default: NOP, default class is class(C)
        accum = NULL ;
        if (!GB_mx_mxArray_to_BinaryOp (&accum, pargin [2], "accum",
            GB_NOP_opcode, cclass, C->type == Complex, A->type == Complex))
        {
            FREE_ALL ;
            mexErrMsgTxt ("accum failed") ;
        }

        // get I
        if (!GB_mx_mxArray_to_indices (&I, pargin [4], &ni))
        {
            FREE_ALL ;
            mexErrMsgTxt ("I failed") ;
        }

        // get J
        if (!GB_mx_mxArray_to_indices (&J, pargin [5], &nj))
        {
            FREE_ALL ;
            mexErrMsgTxt ("J failed") ;
        }

        // get desc
        if (!GB_mx_mxArray_to_Descriptor (&desc, PARGIN (6), "desc"))
        {
            FREE_ALL ;
            mexErrMsgTxt ("desc failed") ;
        }

        // C<Mask>(I,J) = A

        // GB_thread_local.line = 0 ;

        METHOD (assign ( )) ;
    }

    //--------------------------------------------------------------------------
    // return C to MATLAB as a struct
    //--------------------------------------------------------------------------

    ASSERT_OK (GB_check (C, "Final C before wait", 0)) ;
    GrB_wait ( ) ;
    pargout [0] = GB_mx_Matrix_to_mxArray (&C, "C assign result", true) ;
    FREE_ALL ;
}


//------------------------------------------------------------------------------
// GrB_Vector_select: select entries from a vector
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

GrB_Info GxB_Vector_select          // w<mask> = accum (w, select(u,k))
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GxB_SelectOp op,          // operator to apply to the entries
    const GrB_Vector u,             // first input:  vector u
    const void *k,                  // optional input for select operator
    const GrB_Descriptor desc       // descriptor for w and mask
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    WHERE ("GxB_Vector_select (w, mask, accum, op, u, k, desc)") ;
    RETURN_IF_NULL_OR_UNINITIALIZED (w) ;
    RETURN_IF_UNINITIALIZED (mask) ;
    RETURN_IF_NULL_OR_UNINITIALIZED (u) ;

    // get the descriptor
    GET_DESCRIPTOR (info, desc, C_replace, Mask_comp, ignore0, ignore1) ;

    //--------------------------------------------------------------------------
    // select the entries; do not transpose; assemble pending entries
    //--------------------------------------------------------------------------

    return (GB_select (
        (GrB_Matrix) w,     C_replace,      // w and its descriptor
        (GrB_Matrix) mask,  Mask_comp,      // mask and its descriptor
        accum,                              // optional accum for Z=accum(C,T)
        op,                                 // operator to select the entries
        (GrB_Matrix) u,                     // first input: u
        k,                                  // optional input for select op
        false)) ;                           // u, not transposed
}


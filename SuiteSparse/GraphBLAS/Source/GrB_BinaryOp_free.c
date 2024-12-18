//------------------------------------------------------------------------------
// GrB_BinaryOp_free: free a binary operator
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

GrB_Info GrB_BinaryOp_free          // free a user-created binary operator
(
    GrB_BinaryOp *binaryop          // handle of binary operator to free
)
{

    if (binaryop != NULL)
    {
        // only free a user-defined operator, not a built-in one
        GrB_BinaryOp op = *binaryop ;
        if (op != NULL && op->opcode == GB_USER_opcode)
        {
            if (op->magic == MAGIC)
            {
                op->magic = FREED ;         // to help detect dangling pointers
                GB_FREE_MEMORY (*binaryop) ;
            }
            (*binaryop) = NULL ;
        }
    }

    return (GrB_SUCCESS) ;
}


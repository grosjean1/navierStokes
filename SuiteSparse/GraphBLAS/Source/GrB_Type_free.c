//------------------------------------------------------------------------------
// GrB_Type_free:  free a user-defined type
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

#include "GB.h"

GrB_Info GrB_Type_free          // free a user-defined type
(
    GrB_Type *type              // handle of user-defined type to free
)
{

    if (type != NULL)
    {
        // only free a user-defined type, not a built-in one
        GrB_Type t = *type ;
        if (t != NULL && t->code == GB_UDT_code)
        {
            if (t->magic == MAGIC)
            {
                t->magic = FREED ;       // to help detect dangling pointers
                GB_FREE_MEMORY (*type) ;
            }
            (*type) = NULL ;
        }
    }

    return (GrB_SUCCESS) ;
}


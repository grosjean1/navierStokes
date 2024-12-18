//------------------------------------------------------------------------------
// GraphBLAS.h: definitions for the GraphBLAS package
//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017, All Rights Reserved.
// http://suitesparse.com   See GraphBLAS/Doc/License.txt for license.

//------------------------------------------------------------------------------

// SuiteSparse:GraphBLAS is an full implementation of the GraphBLAS standard,
// which defines a set of sparse matrix operations on an extended algebra of
// semirings, using an almost unlimited variety of operators and types.  When
// applied to sparse adjacency matrices, these algebraic operations are
// equivalent to computations on graphs.  GraphBLAS provides a powerful and
// expressive framework creating graph algorithms based on the elegant
// mathematics of sparse matrix operations on a semiring.

// This GraphBLAS.h file contains GraphBLAS definitions for user applications
// to #include.  Functions and variables with GB_ or _opaque in their name need
// to be defined in this file and are thus technically visible to the user, but
// they must not be not accessed in user code.  They are not guaranteed to be
// present in all implementations of GraphBLAS.

// The GraphBLAS API Specification 1.1.0 is provisional, but this
// implementation fully conforms to that specificatin.  It does include
// functions and features that are extensions to the spec.  These are cataloged
// here and tagged with "SPEC".

// All functions and definitions that are extensions to the spec are given
// names of the form GxB_* for functions and built-in objects, or GXB_ for
// macros, so it is clear which are in the spec and which are extensions.

#ifndef GRAPHBLAS_H
#define GRAPHBLAS_H

//------------------------------------------------------------------------------
// GraphBLAS version
//------------------------------------------------------------------------------

// SPEC: the following macros are extensions to the spec

// There are two version numbers that user codes can check against with
// compile-time #if tests:  the version of this GraphBLAS implementation,
// and the version of the GraphBLAS specification it conforms to.  User code
// can use tests like this:
//
//      #if GXB >= GXB_VERSION (2,0,3)
//      ... use features in GraphBLAS specification 2.0.3 ...
//      #else
//      ... only use features in early specifications
//      #endif
//
//      #if GXB_IMPLEMENTATION > GXB_VERSION (1,4,0)
//      ... use features from version 1.4.0 of a GraphBLAS package
//      #endif

// X_GRAPHBLAS: names this particular implementation:
#define GXB_SUITESPARSE_GRAPHBLAS

// GXB_VERSION: a single integer for comparing spec and version levels
#define GXB_VERSION(major,minor,sub) \
    (((major)*1000ULL + (minor))*1000ULL + (sub))

// This implementation conforms to the GraphBLAS provisional release 1.1.0
#define GXB_MAJOR 1
#define GXB_MINOR 1
#define GXB_SUB   0
#define GXB GXB_VERSION(GXB_MAJOR, GXB_MINOR, GXB_SUB)

// The version of this implementation:
#define GXB_IMPLEMENTATION_MAJOR 1
#define GXB_IMPLEMENTATION_MINOR 1
#define GXB_IMPLEMENTATION_SUB   0
#define GXB_IMPLEMENTATION \
        GXB_VERSION (GXB_IMPLEMENTATION_MAJOR, \
                     GXB_IMPLEMENTATION_MINOR, \
                     GXB_IMPLEMENTATION_SUB)

// The 'about' string the describes this particular implementation of GraphBLAS:
#define GXB_ABOUT \
"SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017, All Rights Reserved.\n" \
"http://suitesparse.com  Dept of Computer Sci. & Eng, Texas A&M University\n"

// and its date:
#define GXB_DATE "Dec 1, 2017"

// The 'spec' string describes the GraphBLAS spec:
#define GXB_SPEC \
"GraphBLAS C API, provisional release, by Aydin Buluc, Timothy\n"   \
"Mattson, Scott McMillan, Jose' Moreira, Carl Yang.  Based on\n"    \
"\"GraphBLAS Mathematics\" by Jeremy Kepner.\n"

// and its date:
#define GXB_SPEC_DATE "Oct 10, 2017"

// The GraphBLAS license for this particular implementation of GraphBLAS:
#define GXB_LICENSE \
"SuiteSparse:GraphBLAS, Copyright 2017, Timothy A. Davis\n"                  \
"\n"                                                                         \
"Licensed under the Apache License, Version 2.0 (the \"License\");\n"        \
"you may not use SuiteSparse:GraphBLAS except in compliance with the\n"      \
"License.  You may obtain a copy of the License at\n"                        \
"\n"                                                                         \
"    http://www.apache.org/licenses/LICENSE-2.0  \n"                         \
"\n"                                                                         \
"Unless required by applicable law or agreed to in writing, software\n"      \
"distributed under the License is distributed on an \"AS IS\" BASIS,\n"      \
"WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.\n" \
"See the License for the specific language governing permissions and\n"      \
"limitations under the License.\n"

//------------------------------------------------------------------------------
// include files required by GraphBLAS
//------------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <inttypes.h>
#include <stddef.h>
#include <limits.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

//------------------------------------------------------------------------------
// GraphBLAS error and informational codes
//------------------------------------------------------------------------------

// All GraphBLAS functions return a code that indicates if it was successful
// or not.  If more information is required, the GrB_error function can be
// called, which returns a string that provides more information on the last
// return value from GraphBLAS.

typedef enum
{

    GrB_SUCCESS,                // all is well

    //--------------------------------------------------------------------------
    // informational codes, not an error:
    //--------------------------------------------------------------------------

    // The GraphBLAS spec lists this as an 'error' code, but it just means that
    // A(i,j) is not present in the matrix, having been requested by
    // GrB_*_extractElement.  The function cannot return the proper value
    // because the value of 'implicit zeros' depends on the semiring.  For the
    // conventational plus-times semiring, the implied 'zero' actually has the
    // value of zero.  For the max-plus semiring, it has the value -infinity.
    // A matrix does not keep track of its semiring, and the user can change
    // the semiring used to operate on the matrix.  How mathematically
    // well-defined that change of semiring is depends the user; GraphBLAS will
    // not change the explicit values in the matrix if the semiring changes.
    // As a result, GraphBLAS needs to return not a value, but an indication
    // that the value of A(i,j) is implicit.  The user application can use this
    // indicator (GrB_NO_VALUE) to use the semiring's addititive identity, or
    // it can take other action, as it chooses.  In either case, it is safe to
    // ask for values that are not there, which is why this return condition is
    // not really an 'error' code but an informational code.

    GrB_NO_VALUE,               // A(i,j) requested but not there

    //--------------------------------------------------------------------------
    // API errors:
    //--------------------------------------------------------------------------

    // In non-blocking mode, these errors are caught right away.

    GrB_UNINITIALIZED_OBJECT,   // object has not been initialized
    GrB_INVALID_OBJECT,         // object is corrupted
    GrB_NULL_POINTER,           // input pointer is NULL
    GrB_INVALID_VALUE,          // generic error code; some value is bad
    GrB_INVALID_INDEX,          // a row or column index is out of bounds;
                                // used for indices passed as scalars, not
                                // in a list.
    GrB_DOMAIN_MISMATCH,        // object domains are not compatible
    GrB_DIMENSION_MISMATCH,     // matrix dimensions do not match
    GrB_OUTPUT_NOT_EMPTY,       // output matrix already has values in it

    //--------------------------------------------------------------------------
    // execution errors:
    //--------------------------------------------------------------------------

    // In non-blocking mode, these errors can be deferred.

    GrB_OUT_OF_MEMORY,          // out of memory
    GrB_INSUFFICIENT_SPACE,     // output array not large enough
    GrB_INDEX_OUT_OF_BOUNDS,    // a row or column index is out of bounds;
                                // used for indices in a list of indices.
    GrB_PANIC                   // SuiteSparse:GraphBLAS never panics

}
GrB_Info ;

//==============================================================================
//=== GraphBLAS context methods ================================================
//==============================================================================

// GrB_init must called before any other GraphBLAS operation.  GrB_finalize
// must be called as the last GraphBLAS operation.

// GrB_init defines the mode that GraphBLAS will use:  blocking or
// non-blocking.  With blocking mode, all operations finish before returning to
// the user application.  With non-blocking mode, operations can be left
// pending, and are computed only when needed.

// The GrB_wait ( ) function forces all pending operations to complete.
// Blocking mode is as if GrB_wait is called whenever a GraphBLAS method or
// operation returns to the user.

// The non-blocking mode is unpredictable if user-defined functions have side
// effects or if they rely on global variables, which are not under the control
// of GraphBLAS.  Suppose the user application creates a user-defined operator
// that accesses a global variable.  That operator is then used in a GraphBLAS
// operation, which is left pending.  If the user application then changes the
// global variable before pending operations complete, the pending operations
// will be eventually computed with this different value.

// The non-blocking mode can have side effects if user-defined functions have
// side effects or if they rely on global variables, which are not under the
// control of GraphBLAS.  Suppose the user creates a user-defined operator that
// accesses a global variable.  That operator is then used in a GraphBLAS
// operation, which is left pending.  If the user then changes the global
// variable before pending operations complete, the pending operations will be
// eventually computed with this different value.

// Worse yet, a user-defined operator might be freed before it is needed to
// finish a pending operation.  This causes undefined behavior.  To avoid this,
// call GrB_wait before modifying any global variables relied upon by
// user-defined operators, or before freeing any user-defined types, operators,
// monoids, or semirings.

typedef enum
{
    GrB_NONBLOCKING,    // methods may return with pending computations
    GrB_BLOCKING        // no computations are ever left pending
}
GrB_Mode ;

GrB_Info GrB_init           // start up GraphBLAS
(
    const GrB_Mode mode     // blocking or non-blocking mode
) ;

// In non-blocking mode, GraphBLAS operations need not complete until their
// results are required.  GrB_wait ensures all pending operations are finished.

GrB_Info GrB_wait ( ) ;     // finish all pending computations

// GrB_finalize does not call GrB_wait; any pending computations are abandoned.

GrB_Info GrB_finalize ( ) ;     // finish GraphBLAS

//==============================================================================
//=== GraphBLAS sequence termination ===========================================
//==============================================================================

// Each GraphBLAS method and operation returns a GrB_Info error code.
// GrB_error returns additional information on the error in a thread-safe
// null-terminated string.  The string returned by GrB_error is statically
// allocated in thread local storage and must not be free'd.

const char *GrB_error ( ) ;     // return a string describing the last error

//==============================================================================
//=== GraphBLAS types, operators, monoids, and semirings =======================
//==============================================================================

//------------------------------------------------------------------------------
// GraphBLAS types
//------------------------------------------------------------------------------

// A GraphBLAS GrB_Type defines the type of scalar values that a matrix
// contains, and the type of scalar operands for a unary or binary operator.
// There are eleven built-in types, and a user application can define any types
// of its own as well.  The built-in types correspond to built-in types in C
// (#include <stdbool.h> #include <stdint.h>), and the classes in MATLAB, as
// listed in the "extern GrB_Type ..." table below The user application can
// also define new types based on any typedef in the C language whose values
// are held in a contiguous region of memory.

#define GB_LEN 128

typedef struct
{
    int64_t magic ;         // for detecting uninitialized objects
    size_t size ;           // size of the type
    int code ;              // the type code
    char name [GB_LEN] ;    // name of the type
}
GB_Type_opaque ;            // CONTENT NOT USER-ACCESSIBLE

// The GrB_Type handle user-accessible, but GB_Type_opaque is not:
typedef GB_Type_opaque *GrB_Type ;

// GraphBLAS predefined types and the counterparts in pure C and in MATLAB
extern GrB_Type
    GrB_BOOL   ,        // in C: bool       in MATLAB: logical
    GrB_INT8   ,        // in C: int8_t     in MATLAB: int8
    GrB_UINT8  ,        // in C: uint8_t    in MATLAB: uint8
    GrB_INT16  ,        // in C: int16_t    in MATLAB: int16
    GrB_UINT16 ,        // in C: uint16_t   in MATLAB: uint16
    GrB_INT32  ,        // in C: int32_t    in MATLAB: int32
    GrB_UINT32 ,        // in C: uint32_t   in MATLAB: uint32
    GrB_INT64  ,        // in C: int64_t    in MATLAB: int64
    GrB_UINT64 ,        // in C: uint64_t   in MATLAB: uint64
    GrB_FP32   ,        // in C: float      in MATLAB: single
    GrB_FP64   ;        // in C: double     in MATLAB: double

// The user-callable function has the following signature.
// It is actually implemented as a macro.

/*

GrB_Info GrB_Type_new           // create a new GraphBLAS type
(
    GrB_Type *type,             // handle of user type to create
    <ctype>                     // a C type
) ;

*/

// GB_STR: convert the content of x into a string "x"
#define GB_XSTR(x) GB_STR(x)
#define GB_STR(x) #x

// GrB_Type_new is user-callable; GB_Type_new should not be called directly.
#define GrB_Type_new(utype, ctype) \
    GB_Type_new (utype, sizeof (ctype), GB_STR(ctype))

// This function is not user-callable; use GrB_Type_new instead

GrB_Info GB_Type_new            // create a new GraphBLAS type
(
    GrB_Type *type,             // handle of user type to create
    const size_t size,          // size of the user type
    const char *name            // name of the type
) ;

// SPEC: GxB_Type_size is an extension to the spec

GrB_Info GxB_Type_size          // determine the size of the type
(
    size_t *size,               // the sizeof the type
    GrB_Type type               // type to determine the sizeof
) ;

GrB_Info GrB_Type_free          // free a user-defined type
(
    GrB_Type *type              // handle of user-defined type to free
) ;


//------------------------------------------------------------------------------
// GraphBLAS unary and binary operators
//------------------------------------------------------------------------------

// GraphBLAS defines built-in unary and binary operators, and the user may also
// define new ones via function pointers.  When a user function z=f(x,y) or
// z=f(x) is called by GraphBLAS, the pointers x, y, and z are guaranteed to be
// non-NULL and to point to unique valid space of the expected type.

// GraphBLAS provides 256 built-in binary operators z=f(x,y) and 45 built-in
// unary operators z=f(x) that operate on the 11 built-in types.  Built-in
// types are statically allocated and need not be freed when the application
// finishes.

//------------------------------------------------------------------------------
// unary operators
//------------------------------------------------------------------------------

// GrB_UnaryOp: a function z=f(x).  The function f must have the signature:

//      void f (void *z, const void *x) ;

// The pointers are void * but they are always of pointers to objects of type
// ztype and xtype, respectively.  The function must typecast its arguments as
// needed from void* to ztype* and xtype*.

typedef struct
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Type xtype ;        // type of x
    GrB_Type ztype ;        // type of z
    void *function ;        // a pointer to the unary function
    char name [GB_LEN] ;    // name of the unary operator
    int opcode ;            // operator opcode
}
GB_UnaryOp_opaque ;         // CONTENT NOT USER-ACCESSIBLE

// The GrB_UnaryOp handle (user-accesible)
typedef GB_UnaryOp_opaque *GrB_UnaryOp ;

//------------------------------------------------------------------------------
// built-in unary operators, z = f(x)
//------------------------------------------------------------------------------

// There are 67 unary operators: 6 kinds * 11 types and GrB_LNOT.

// SPEC: ONE and ABS unary operators are extensions to the spec
// as are the LNOT_TYPE operators.

extern GrB_UnaryOp
    // For these three functions z=f(x), z and x have the same type.
    // The suffix in the name is the type of x and z.
    // z = x             z = -x             z = 1/x             z = ! (x != 0)
    // identity          additive           multiplicative      logical
    //                   inverse            inverse             negation
    GrB_IDENTITY_BOOL,   GrB_AINV_BOOL,     GrB_MINV_BOOL,      GxB_LNOT_BOOL,
    GrB_IDENTITY_INT8,   GrB_AINV_INT8,     GrB_MINV_INT8,      GxB_LNOT_INT8,
    GrB_IDENTITY_UINT8,  GrB_AINV_UINT8,    GrB_MINV_UINT8,     GxB_LNOT_UINT8,
    GrB_IDENTITY_INT16,  GrB_AINV_INT16,    GrB_MINV_INT16,     GxB_LNOT_INT16,
    GrB_IDENTITY_UINT16, GrB_AINV_UINT16,   GrB_MINV_UINT16,    GxB_LNOT_UINT16,
    GrB_IDENTITY_INT32,  GrB_AINV_INT32,    GrB_MINV_INT32,     GxB_LNOT_INT32,
    GrB_IDENTITY_UINT32, GrB_AINV_UINT32,   GrB_MINV_UINT32,    GxB_LNOT_UINT32,
    GrB_IDENTITY_INT64,  GrB_AINV_INT64,    GrB_MINV_INT64,     GxB_LNOT_INT64,
    GrB_IDENTITY_UINT64, GrB_AINV_UINT64,   GrB_MINV_UINT64,    GxB_LNOT_UINT64,
    GrB_IDENTITY_FP32,   GrB_AINV_FP32,     GrB_MINV_FP32,      GxB_LNOT_FP32,
    GrB_IDENTITY_FP64,   GrB_AINV_FP64,     GrB_MINV_FP64,      GxB_LNOT_FP64,

    // z = 1             z = abs(x)
    // one               absolute value
    //
    GxB_ONE_BOOL,        GxB_ABS_BOOL,
    GxB_ONE_INT8,        GxB_ABS_INT8,
    GxB_ONE_UINT8,       GxB_ABS_UINT8,
    GxB_ONE_INT16,       GxB_ABS_INT16,
    GxB_ONE_UINT16,      GxB_ABS_UINT16,
    GxB_ONE_INT32,       GxB_ABS_INT32,
    GxB_ONE_UINT32,      GxB_ABS_UINT32,
    GxB_ONE_INT64,       GxB_ABS_INT64,
    GxB_ONE_UINT64,      GxB_ABS_UINT64,
    GxB_ONE_FP32,        GxB_ABS_FP32,
    GxB_ONE_FP64,        GxB_ABS_FP64,

    // Boolean negation, z = !x, where both x and x are boolean.  There is no
    // suffix since z and x are only boolean.  This operator is identical to
    // GxB_LNOT_BOOL; it just has a different name.
    GrB_LNOT ;

//------------------------------------------------------------------------------
// methods for unary operators
//------------------------------------------------------------------------------

// The user-callable function GrB_UnaryOp_new has the following signature.  It
// is actually implemented as a macro so that the name of the unary function
// can be kept by GraphBLAS.

/*

GrB_Info GrB_UnaryOp_new            // create a new user-defined unary operator
(
    GrB_UnaryOp *unaryop,           // handle for the new unary operator
    void *function,                 // pointer to the unary function
    const GrB_Type ztype,           // type of output z
    const GrB_Type xtype            // type of input x
) ;

*/

#define GrB_UnaryOp_new(op,f,z,x) GB_UnaryOp_new (op,f,z,x, GB_STR(f))

// This function is NOT user-callable:

GrB_Info GB_UnaryOp_new             // create a new user-defined unary operator
(
    GrB_UnaryOp *unaryop,           // handle for the new unary operator
    void *function,                 // pointer to the unary function
    const GrB_Type ztype,           // type of output z
    const GrB_Type xtype,           // type of input x
    const char *name                // name of the underlying function
) ;

// SPEC: GxB_UnaryOp_ztype is an extension to the spec

GrB_Info GxB_UnaryOp_ztype          // return the type of z
(
    GrB_Type *ztype,                // return type of output z
    const GrB_UnaryOp unaryop       // unary operator
) ;

// SPEC: GxB_UnaryOp_xtype is an extension to the spec

GrB_Info GxB_UnaryOp_xtype          // return the type of x
(
    GrB_Type *xtype,                // return type of input x
    const GrB_UnaryOp unaryop       // unary operator
) ;

GrB_Info GrB_UnaryOp_free           // free a user-created unary operator
(
    GrB_UnaryOp *unaryop            // handle of unary operator to free
) ;

//------------------------------------------------------------------------------
// binary operators
//------------------------------------------------------------------------------

// GrB_BinaryOp: a function z=f(x,y).  The function f must have the signature:

//      void f (void *z, const void *x, const void *y) ;

// The pointers are void * but they are always of pointers to objects of type
// ztype, xtype, and ytype, respectively.  See Demo/usercomplex.c for examples.

typedef struct
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Type xtype ;        // type of x
    GrB_Type ytype ;        // type of y
    GrB_Type ztype ;        // type of z
    void *function ;        // a pointer to the binary function
    char name [GB_LEN] ;    // name of the binary operator
    int opcode ;            // operator opcode
}
GB_BinaryOp_opaque ;        // CONTENT NOT USER-ACCESSIBLE

// The GrB_BinaryOp handle (user-accesible)
typedef GB_BinaryOp_opaque *GrB_BinaryOp ;

//------------------------------------------------------------------------------
// built-in binary operators, z = f(x,y)
//------------------------------------------------------------------------------

// There are three sets of built-in binary operators.  For the first set of
// 17 kinds of operators, x,y,z all have the same type, and they are available
// for all 11 types, for a total of 17*11 = 187 operators.  All of them have
// a "_TYPE" suffix that denotes the type of x,y,z:

//      8 general: FIRST, SECOND, MIN, MAX, PLUS, MINUS, TIMES, DIV
//      6 comparison: ISEQ, ISNE, ISGT, ISLT, ISGE, ISLE
//      3 logical: LOR, LAND, LXOR

// For the second set, there are 6 comparison operators where, x,y have the
// same type but z is always boolean, for a total of 6*11 = 66 operators.
// All of them have a "_TYPE" suffix that denotes the type of x,y (not z):

//      6 comparison: EQ, NE, GT, LT, GE, LE

// The final set of operators is for boolean x,y,z only, and they have no
// suffix:

//      3 logical: LOR, LAND, LXOR

// Thus there are 187+66+3 = 256 built-in binary operators.  Some are redundant
// but are included to keep the name space of operators uniform.

// For eight binary operators z=f(x,y), x, y, and z are all the same type:
// FIRST, SECOND, MIN, MAX, PLUS, MINUS, TIMES, and DIV, for all 11 types.

extern GrB_BinaryOp
    // z = x            z = y               z = min(x,y)        z = max (x,y)
    GrB_FIRST_BOOL,     GrB_SECOND_BOOL,    GrB_MIN_BOOL,       GrB_MAX_BOOL,
    GrB_FIRST_INT8,     GrB_SECOND_INT8,    GrB_MIN_INT8,       GrB_MAX_INT8,
    GrB_FIRST_UINT8,    GrB_SECOND_UINT8,   GrB_MIN_UINT8,      GrB_MAX_UINT8,
    GrB_FIRST_INT16,    GrB_SECOND_INT16,   GrB_MIN_INT16,      GrB_MAX_INT16,
    GrB_FIRST_UINT16,   GrB_SECOND_UINT16,  GrB_MIN_UINT16,     GrB_MAX_UINT16,
    GrB_FIRST_INT32,    GrB_SECOND_INT32,   GrB_MIN_INT32,      GrB_MAX_INT32,
    GrB_FIRST_UINT32,   GrB_SECOND_UINT32,  GrB_MIN_UINT32,     GrB_MAX_UINT32,
    GrB_FIRST_INT64,    GrB_SECOND_INT64,   GrB_MIN_INT64,      GrB_MAX_INT64,
    GrB_FIRST_UINT64,   GrB_SECOND_UINT64,  GrB_MIN_UINT64,     GrB_MAX_UINT64,
    GrB_FIRST_FP32,     GrB_SECOND_FP32,    GrB_MIN_FP32,       GrB_MAX_FP32,
    GrB_FIRST_FP64,     GrB_SECOND_FP64,    GrB_MIN_FP64,       GrB_MAX_FP64,

    // z = x+y          z = x-y             z = x*y             z = x/y
    GrB_PLUS_BOOL,      GrB_MINUS_BOOL,     GrB_TIMES_BOOL,     GrB_DIV_BOOL,
    GrB_PLUS_INT8,      GrB_MINUS_INT8,     GrB_TIMES_INT8,     GrB_DIV_INT8,
    GrB_PLUS_UINT8,     GrB_MINUS_UINT8,    GrB_TIMES_UINT8,    GrB_DIV_UINT8,
    GrB_PLUS_INT16,     GrB_MINUS_INT16,    GrB_TIMES_INT16,    GrB_DIV_INT16,
    GrB_PLUS_UINT16,    GrB_MINUS_UINT16,   GrB_TIMES_UINT16,   GrB_DIV_UINT16,
    GrB_PLUS_INT32,     GrB_MINUS_INT32,    GrB_TIMES_INT32,    GrB_DIV_INT32,
    GrB_PLUS_UINT32,    GrB_MINUS_UINT32,   GrB_TIMES_UINT32,   GrB_DIV_UINT32,
    GrB_PLUS_INT64,     GrB_MINUS_INT64,    GrB_TIMES_INT64,    GrB_DIV_INT64,
    GrB_PLUS_UINT64,    GrB_MINUS_UINT64,   GrB_TIMES_UINT64,   GrB_DIV_UINT64,
    GrB_PLUS_FP32,      GrB_MINUS_FP32,     GrB_TIMES_FP32,     GrB_DIV_FP32,
    GrB_PLUS_FP64,      GrB_MINUS_FP64,     GrB_TIMES_FP64,     GrB_DIV_FP64 ;

// Six comparison operators z=f(x,y) return the same type as their inputs.
// Each of them compute z = (x OP y), where x, y, and z all have the same type.
// The value z is either 1 for true or 0 for false, but it is a value with the
// same type as x and y.  Z is not bool (unless x and y are also bool).  These
// operators compute the same thing as the 6 sets of EQ, NE, GT, LT, GE, and LE
// operators.  They just return their result z as the same type as x and y,
// instead of returning a value z that is boolean.  Since their ztype is
// non-boolean, they can be used as multiply operators in a semring with
// non-boolean monoids (PLUS, for example).

extern GrB_BinaryOp
    // z = (x == y)     z = (x != y)        z = (x > y)         z = (x < y)
    GxB_ISEQ_BOOL,      GxB_ISNE_BOOL,      GxB_ISGT_BOOL,      GxB_ISLT_BOOL,
    GxB_ISEQ_INT8,      GxB_ISNE_INT8,      GxB_ISGT_INT8,      GxB_ISLT_INT8,
    GxB_ISEQ_UINT8,     GxB_ISNE_UINT8,     GxB_ISGT_UINT8,     GxB_ISLT_UINT8,
    GxB_ISEQ_INT16,     GxB_ISNE_INT16,     GxB_ISGT_INT16,     GxB_ISLT_INT16,
    GxB_ISEQ_UINT16,    GxB_ISNE_UINT16,    GxB_ISGT_UINT16,    GxB_ISLT_UINT16,
    GxB_ISEQ_INT32,     GxB_ISNE_INT32,     GxB_ISGT_INT32,     GxB_ISLT_INT32,
    GxB_ISEQ_UINT32,    GxB_ISNE_UINT32,    GxB_ISGT_UINT32,    GxB_ISLT_UINT32,
    GxB_ISEQ_INT64,     GxB_ISNE_INT64,     GxB_ISGT_INT64,     GxB_ISLT_INT64,
    GxB_ISEQ_UINT64,    GxB_ISNE_UINT64,    GxB_ISGT_UINT64,    GxB_ISLT_UINT64,
    GxB_ISEQ_FP32,      GxB_ISNE_FP32,      GxB_ISGT_FP32,      GxB_ISLT_FP32,
    GxB_ISEQ_FP64,      GxB_ISNE_FP64,      GxB_ISGT_FP64,      GxB_ISLT_FP64,

    // z = (x >= y)     z = (x <= y)
    GxB_ISGE_BOOL,      GxB_ISLE_BOOL,
    GxB_ISGE_INT8,      GxB_ISLE_INT8,
    GxB_ISGE_UINT8,     GxB_ISLE_UINT8,
    GxB_ISGE_INT16,     GxB_ISLE_INT16,
    GxB_ISGE_UINT16,    GxB_ISLE_UINT16,
    GxB_ISGE_INT32,     GxB_ISLE_INT32,
    GxB_ISGE_UINT32,    GxB_ISLE_UINT32,
    GxB_ISGE_INT64,     GxB_ISLE_INT64,
    GxB_ISGE_UINT64,    GxB_ISLE_UINT64,
    GxB_ISGE_FP32,      GxB_ISLE_FP32,
    GxB_ISGE_FP64,      GxB_ISLE_FP64 ;

// Six comparison operators z=f(x,y) return their result as boolean, but where
// x and y have the same type (any one of the 11 built-in types).  The suffix
// in their names refers to the type of x and y since z is always boolean.  If
// used as multiply operators in a semiring, they can only be combined with
// boolean monoids.  The _BOOL versions of these operators give the same
// results as their IS*_BOOL counterparts.

extern GrB_BinaryOp
    // z = (x == y)     z = (x != y)        z = (x > y)         z = (x < y)
    GrB_EQ_BOOL,        GrB_NE_BOOL,        GrB_GT_BOOL,        GrB_LT_BOOL,
    GrB_EQ_INT8,        GrB_NE_INT8,        GrB_GT_INT8,        GrB_LT_INT8,
    GrB_EQ_UINT8,       GrB_NE_UINT8,       GrB_GT_UINT8,       GrB_LT_UINT8,
    GrB_EQ_INT16,       GrB_NE_INT16,       GrB_GT_INT16,       GrB_LT_INT16,
    GrB_EQ_UINT16,      GrB_NE_UINT16,      GrB_GT_UINT16,      GrB_LT_UINT16,
    GrB_EQ_INT32,       GrB_NE_INT32,       GrB_GT_INT32,       GrB_LT_INT32,
    GrB_EQ_UINT32,      GrB_NE_UINT32,      GrB_GT_UINT32,      GrB_LT_UINT32,
    GrB_EQ_INT64,       GrB_NE_INT64,       GrB_GT_INT64,       GrB_LT_INT64,
    GrB_EQ_UINT64,      GrB_NE_UINT64,      GrB_GT_UINT64,      GrB_LT_UINT64,
    GrB_EQ_FP32,        GrB_NE_FP32,        GrB_GT_FP32,        GrB_LT_FP32,
    GrB_EQ_FP64,        GrB_NE_FP64,        GrB_GT_FP64,        GrB_LT_FP64,

    // z = (x >= y)     z = (x <= y)
    GrB_GE_BOOL,        GrB_LE_BOOL,
    GrB_GE_INT8,        GrB_LE_INT8,
    GrB_GE_UINT8,       GrB_LE_UINT8,
    GrB_GE_INT16,       GrB_LE_INT16,
    GrB_GE_UINT16,      GrB_LE_UINT16,
    GrB_GE_INT32,       GrB_LE_INT32,
    GrB_GE_UINT32,      GrB_LE_UINT32,
    GrB_GE_INT64,       GrB_LE_INT64,
    GrB_GE_UINT64,      GrB_LE_UINT64,
    GrB_GE_FP32,        GrB_LE_FP32,
    GrB_GE_FP64,        GrB_LE_FP64 ;

// Three binary operators operate on each of the types, converting them
// internally to boolean and returning a value 1 or 0 in the same type, for
// true or false.  Each of them compute z = ((x != 0) OP (y != 0)), where x, y,
// and z all the same type.  These operators are useful as multiply operators
// when combined with non-boolean monoids of the same type.

extern GrB_BinaryOp
    // z = (x || y)     z = (x && y)        z = (x != y)
    GxB_LOR_BOOL,       GxB_LAND_BOOL,      GxB_LXOR_BOOL,
    GxB_LOR_INT8,       GxB_LAND_INT8,      GxB_LXOR_INT8,
    GxB_LOR_UINT8,      GxB_LAND_UINT8,     GxB_LXOR_UINT8,
    GxB_LOR_INT16,      GxB_LAND_INT16,     GxB_LXOR_INT16,
    GxB_LOR_UINT16,     GxB_LAND_UINT16,    GxB_LXOR_UINT16,
    GxB_LOR_INT32,      GxB_LAND_INT32,     GxB_LXOR_INT32,
    GxB_LOR_UINT32,     GxB_LAND_UINT32,    GxB_LXOR_UINT32,
    GxB_LOR_INT64,      GxB_LAND_INT64,     GxB_LXOR_INT64,
    GxB_LOR_UINT64,     GxB_LAND_UINT64,    GxB_LXOR_UINT64,
    GxB_LOR_FP32,       GxB_LAND_FP32,      GxB_LXOR_FP32,
    GxB_LOR_FP64,       GxB_LAND_FP64,      GxB_LXOR_FP64 ;

// Finally, three binary operate only on boolean types: LOR, LAND, LXOR.  The
// naming convention differs (_BOOL is not appended to the name).  They are
// the same as GxB_LOR_BOOL, GxB_LAND_BOOL, and GxB_LXOR_BOOL; they just
// have a simpler name.

extern GrB_BinaryOp
    // z = (x || y)     z = (x && y)        z = (x != y)
    GrB_LOR,            GrB_LAND,           GrB_LXOR ;

// Some of the boolean operators compute the same thing but have unique names.
// For example, x*y and x&&y give the same results for boolean x and y.
// Operations such as x < y when x and y are boolean are treated as if true=1
// and false=0.  Below is the truth table for all 17 binary operators with
// boolean inputs.  This table is defined by how C typecasts boolean values for
// non-boolean operations.  For example, if x, y, and z are boolean, x = true,
// and y = true, then z = x + y = true + true = true.  DIV (x/y) is defined
// below.

//                                                     is  is  is  is  is  is
//  x y    1st 2nd min max +   -   *   /   or  and xor eq  ne  gt  lt  ge  le
//  0 0    0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   1   1
//  0 1    0   1   0   1   1   1   0   0   1   0   1   0   1   0   1   0   1
//  1 0    1   0   0   1   1   1   0   1   1   0   1   0   1   1   0   1   0
//  1 1    1   1   1   1   1   0   1   1   1   1   0   1   0   0   0   1   1

// SPEC: the definition of divide-by-zero is an extension to the spec

// GraphBLAS includes a GrB_DIV_BOOL operator in its specification, but does
// not define what boolean "division" means.  SuiteSparse:GraphBLAS makes the
// following interpretation.

// GraphBLAS does not generate exceptions for divide-by-zero, so the results
// are defined just as they are in MATLAB.  Floating-point divide-by-zero
// follows the IEEE 754 standard: 1/0 is +Inf, -1/0 is -Inf, and 0/0 is NaN.
// For integer division by zero, if x is positive, x/0 is the largest integer,
// -x/0 is the integer minimum (zero for unsigned integers), and 0/0 is zero.
// For example, for int8, 1/0 is 127, and -1/0 is -128.  For uint8, 1/0 is 255
// and 0/0 is zero.

// Boolean division is not in MATLAB.  For SuiteSparse:GraphBLAS, boolean
// division is treated as if it were an unsigned integer type with true=1 and
// false=0, and with the max and min value being 1 and 0.  As a result,
// GrB_IDENTITY_BOOL, GrB_AINV_BOOL, and GrB_MINV_BOOL all give the same result
// (z = x).

// With this convention for boolean "division", there are 10 unique binary
// operators that are purely boolean; 13 *_BOOL operators are redundant but are
// included in GraphBLAS so that the name space of operators is complete:

//      z = x           FIRST, DIV
//      z = y           SECOND
//      z = (x && y)    AND, MIN, TIMES
//      z = (x || y)    OR, MAX, PLUS
//      z = (x != y)    XOR, MINUS, NE, ISNE
//      z = (x == y)    EQ, ISEQ
//      z = (x >  y)    GT, ISGT
//      z = (x <  y)    LT, ISLT
//      z = (x >= y)    GE, ISGE
//      z = (x >= y)    LE, ISLE

// Three more that have no_BOOL suffix are also redundant with the operators
// of the form GxB_*_BOOL.
// (GrB_LOR, GrB_LAND, and GrB_LXOR).

// There are thus 256 built-in binary operators with unique names, 16 of which
// are redundant, giving 240 built-in binary operators that compute unique
// results.

//------------------------------------------------------------------------------
// methods for binary operators
//------------------------------------------------------------------------------

// The user-callable function GxB_BinaryOp_new has the following signature.  It
// is implemented as a macro so that the name of the select function can be
// kept by GraphBLAS.

/*

GrB_Info GrB_BinaryOp_new
(
    GrB_BinaryOp *binaryop,         // handle for the new binary operator
    void *function,                 // pointer to the binary function
    const GrB_Type ztype,           // type of output z
    const GrB_Type xtype,           // type of input x
    const GrB_Type ytype            // type of input y
) ;

*/

#define GrB_BinaryOp_new(op,f,z,x,y) GB_BinaryOp_new (op,f,z,x,y, GB_STR(f))

// This function is NOT user-callable:

GrB_Info GB_BinaryOp_new
(
    GrB_BinaryOp *binaryop,         // handle for the new binary operator
    void *function,                 // pointer to the binary function
    const GrB_Type ztype,           // type of output z
    const GrB_Type xtype,           // type of input x
    const GrB_Type ytype,           // type of input y
    const char *name                // name of the underlying function
) ;


// SPEC: GxB_BinaryOp_ztype is an extension to the spec

GrB_Info GxB_BinaryOp_ztype         // return the type of z
(
    GrB_Type *ztype,                // return type of output z
    const GrB_BinaryOp binaryop     // binary operator to query
) ;

// SPEC: GxB_BinaryOp_xtype is an extension to the spec

GrB_Info GxB_BinaryOp_xtype         // return the type of x
(
    GrB_Type *xtype,                // return type of input x
    const GrB_BinaryOp binaryop     // binary operator to query
) ;

// SPEC: GxB_BinaryOp_ytype is an extension to the spec

GrB_Info GxB_BinaryOp_ytype         // return the type of y
(
    GrB_Type *ytype,                // return type of input y
    const GrB_BinaryOp binaryop     // binary operator to query
) ;

GrB_Info GrB_BinaryOp_free          // free a user-created binary operator
(
    GrB_BinaryOp *binaryop          // handle of binary operator to free
) ;

//------------------------------------------------------------------------------
// Select operators
//------------------------------------------------------------------------------

// SPEC: GxB_SelectOp and all related functions are an extenstion to the spec.

// GxB_SelectOp is an operator used by GxB_select to select entries from an
// input matrix A that are kept in the output C.  If an entry A(i,j) in the
// matrix A, of size nrows-by-ncols, has the value aij, then it calls the
// select function as result = f (i, j, nrows, ncols, aij, k).  If the function
// returns true, the entry is kept in the output C.  If f returns false, the
// entry is not kept in C.  The type of x for the GxB_SelectOp operator may be
// any of the 11 built-in types, or any user-defined type.  It may also be
// GrB_NULL, to indicate that the function is type-generic and does not depend
// at all on the value aij.  In this case, x is passed to f as a NULL pointer.
// The k parameter is a const void * pointer to a user-defined object that is
// passed to GxB_select.  It is not used by GxB_select except to pass it to the
// function f.

// GxB_SelectOp:  a function z=f(i,j,m,n,x,k) for the GxB_Select operation.
// The function f must have the signature:

//      bool f (const GrB_Index i, const GrB_Index j,
//              const GrB_Index nrows, const GrB_Index ncols,
//              const void *x, const void *k) ;

typedef struct
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Type xtype ;        // type of x, or NULL if generic
    void *function ;        // a pointer to the select function
    char name [GB_LEN] ;    // name of the select operator
    int opcode ;            // operator opcode
}
GB_SelectOp_opaque ;        // CONTENT NOT USER-ACCESSIBLE

// The GxB_SelectOp handle (user-accesible)
typedef GB_SelectOp_opaque *GxB_SelectOp ;

//------------------------------------------------------------------------------
// built-in select operators
//------------------------------------------------------------------------------

// GxB_select (C, Mask, accum, op, A, k, desc) always returns a matrix C of the
// same size as A (or A' if GrB_TRAN is in the descriptor).

extern GxB_SelectOp
    GxB_TRIL,       // C=tril(A,k):   returns true if ((j-i) <= k)
    GxB_TRIU,       // C=triu(A,k):   returns true if ((j-i) >= k)
    GxB_DIAG,       // C=diag(A,k):   returns true if ((j-i) == k)
    GxB_OFFDIAG,    // C=A-diag(A,k): returns true if ((j-i) != k)
    GxB_NONZERO ;   // C=A(A~=0):     returns true if aij is nonzero,
                    //                for any built-in or user-defined type

// For GxB_TRIL, GxB_TRIU, GxB_DIAG, and GxB_OFFDIAG, the parameter k is a
// const void * pointer to a single scalar value of type GrB_Index.  These
// select operators do not depend on the values of A, but just their position.

// For GxB_NONZERO, the result depends only on the value of A(i,j), and the k
// parameter may be GrB_NULL.  It works on all built-in types and all
// user-defined types.  When applied to user-defined types the operator it
// returns true if all bits in the user-defined value are zero, which can be
// tested regardless of how the user-defined type is defined.

//------------------------------------------------------------------------------
// select operators
//------------------------------------------------------------------------------

// The user-callable function GxB_SelectOp_new has the following signature.  It
// is implemented as a macro so that the name of the select function can be
// kept by GraphBLAS.

/*

GrB_Info GxB_SelectOp_new       // create a new user-defined select operator
(
    GxB_SelectOp *selectop,     // handle for the new select operator
    void *function,             // pointer to the select function
    const GrB_Type xtype        // type of input x, or NULL if type-generic
) ;

*/

#define GxB_SelectOp_new(op,f,x) GB_SelectOp_new (op,f,x, GB_STR(f))

// This function is NOT user-callable:

GrB_Info GB_SelectOp_new        // create a new user-defined select operator
(
    GxB_SelectOp *selectop,     // handle for the new select operator
    void *function,             // pointer to the select function
    const GrB_Type xtype,       // type of input x
    const char *name            // name of the underlying function
) ;

GrB_Info GxB_SelectOp_xtype     // return the type of x
(
    GrB_Type *xtype,            // return type of input x
    const GxB_SelectOp selectop // select operator
) ;

GrB_Info GxB_SelectOp_free      // free a user-created select operator
(
    GxB_SelectOp *selectop      // handle of select operator to free
) ;

//------------------------------------------------------------------------------
// GraphBLAS Monoid
//------------------------------------------------------------------------------

// A monoid is an associative operator z=op(x,y) where all three types of z, x,
// and y are identical.  The monoid also has an identity element, such that
// op(x,identity) = op(identity,x) = x.

typedef struct
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_BinaryOp op ;       // binary operator of the monoid
    void *identity ;        // identity of the monoid; size is op->ztype->size
    bool identity_is_zero ; // true if all bits of identity are zero
    bool user_defined ;     // true if monoid is user-defined
}
GB_Monoid_opaque ;          // CONTENT NOT USER-ACCESSIBLE

// The GrB_Monoid handle (user-accesible)
typedef GB_Monoid_opaque *GrB_Monoid ;

// Create a new Monoid with a specific type of identity, which must match
// the binary_op type.  The binary_op's three types must all be the same.

GrB_Info GrB_Monoid_BOOL_new        // create a new boolean monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const bool identity             // identity value of the monoid
) ;

GrB_Info GrB_Monoid_INT8_new        // create a new int8 monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const int8_t identity           // identity value of the monoid
) ;

GrB_Info GrB_Monoid_UINT8_new       // create a new uint8 monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const uint8_t identity          // identity value of the monoid
) ;

GrB_Info GrB_Monoid_INT16_new       // create a new int16 monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const int16_t identity          // identity value of the monoid
) ;

GrB_Info GrB_Monoid_UINT16_new      // create a new uint16 monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const uint16_t identity         // identity value of the monoid
) ;

GrB_Info GrB_Monoid_INT32_new       // create a new int32 monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const int32_t identity          // identity value of the monoid
) ;

GrB_Info GrB_Monoid_UINT32_new      // create a new uint32 monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const uint32_t identity         // identity value of the monoid
) ;

GrB_Info GrB_Monoid_INT64_new       // create a new int64 monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const int64_t identity          // identity value of the monoid
) ;

GrB_Info GrB_Monoid_UINT64_new      // create a new uint64 monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const uint64_t identity         // identity value of the monoid
) ;

GrB_Info GrB_Monoid_FP32_new        // create a new float monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const float identity            // identity value of the monoid
) ;

GrB_Info GrB_Monoid_FP64_new        // create a new double monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const double identity           // identity value of the monoid
) ;

GrB_Info GrB_Monoid_UDT_new         // create a monoid with a user-defined type
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const void *identity            // identity value of the monoid
) ;

// Type-generic method for creating a new monoid:

/*

GrB_Info GrB_Monoid_new             // create a monoid
(
    GrB_Monoid *monoid,             // handle of monoid to create
    const GrB_BinaryOp op,          // binary operator of the monoid
    const <type> identity           // identity value of the monoid
) ;

*/

#define GrB_Monoid_new(monoid,op,identity)          \
    _Generic                                        \
    (                                               \
        (identity),                                 \
        const bool     : GrB_Monoid_BOOL_new   ,    \
              bool     : GrB_Monoid_BOOL_new   ,    \
        const int8_t   : GrB_Monoid_INT8_new   ,    \
              int8_t   : GrB_Monoid_INT8_new   ,    \
        const uint8_t  : GrB_Monoid_UINT8_new  ,    \
              uint8_t  : GrB_Monoid_UINT8_new  ,    \
        const int16_t  : GrB_Monoid_INT16_new  ,    \
              int16_t  : GrB_Monoid_INT16_new  ,    \
        const uint16_t : GrB_Monoid_UINT16_new ,    \
              uint16_t : GrB_Monoid_UINT16_new ,    \
        const int32_t  : GrB_Monoid_INT32_new  ,    \
              int32_t  : GrB_Monoid_INT32_new  ,    \
        const uint32_t : GrB_Monoid_UINT32_new ,    \
              uint32_t : GrB_Monoid_UINT32_new ,    \
        const int64_t  : GrB_Monoid_INT64_new  ,    \
              int64_t  : GrB_Monoid_INT64_new  ,    \
        const uint64_t : GrB_Monoid_UINT64_new ,    \
              uint64_t : GrB_Monoid_UINT64_new ,    \
        const float    : GrB_Monoid_FP32_new   ,    \
              float    : GrB_Monoid_FP32_new   ,    \
        const double   : GrB_Monoid_FP64_new   ,    \
              double   : GrB_Monoid_FP64_new   ,    \
        const void *   : GrB_Monoid_UDT_new    ,    \
              void *   : GrB_Monoid_UDT_new         \
    )                                               \
    (monoid, op, identity) ;

// SPEC: GxB_Monoid_operator is an extension to the spec

GrB_Info GxB_Monoid_operator        // return the monoid operator
(
    GrB_BinaryOp *op,               // returns the binary op of the monoid
    const GrB_Monoid monoid         // monoid to query
) ;

// SPEC: GxB_Monoid_identity is an extension to the spec

GrB_Info GxB_Monoid_identity        // return the monoid identity
(
    void *identity,                 // returns the identity of the monoid
    const GrB_Monoid monoid         // monoid to query
) ;

GrB_Info GrB_Monoid_free            // free a user-created monoid
(
    GrB_Monoid *monoid              // handle of monoid to free
) ;

//------------------------------------------------------------------------------
// GraphBLAS Semiring
//------------------------------------------------------------------------------

// A semiring defines all the operators required to define the multiplication
// of two sparse matrices in GraphBLAS, C=A*B.  The "add" operator is a
// commutative and associative monoid, and the binary "multiply" operator
// defines a function z=fmult(x,y) where the type of z matches the exactly with
// the monoid type.

typedef struct
{
    int64_t magic ;                 // for detecting uninitialized objects
    GrB_Monoid add ;                // add operator of the semiring
    GrB_BinaryOp multiply ;         // multiply operator of the semiring
    bool user_defined ;             // true if semiring is user-defined
}
GB_Semiring_opaque ;                // CONTENT NOT USER-ACCESSIBLE

// The GrB_Semiring handle (user-accesible)
typedef GB_Semiring_opaque *GrB_Semiring ;

GrB_Info GrB_Semiring_new           // create a semiring
(
    GrB_Semiring *semiring,         // handle of semiring to create
    const GrB_Monoid add,           // add monoid of the semiring
    const GrB_BinaryOp multiply     // multiply operator of the semiring
) ;

// SPEC: GxB_Semiring_add is an extension to the spec

GrB_Info GxB_Semiring_add           // return the add monoid of a semiring
(
    GrB_Monoid *add,                // returns add monoid of the semiring
    const GrB_Semiring semiring     // semiring to query
) ;

// SPEC: GxB_Semiring_multiply is an extension to the spec

GrB_Info GxB_Semiring_multiply      // return multiply operator of a semiring
(
    GrB_BinaryOp *multiply,         // returns multiply operator of the semiring
    const GrB_Semiring semiring     // semiring to query
) ;

GrB_Info GrB_Semiring_free          // free a user-created semiring
(
    GrB_Semiring *semiring          // handle of semiring to free
) ;

//==============================================================================
//=== GraphBLAS Matrix and Vector objects ======================================
//==============================================================================

// Sparse matrices and vectors are the primary objects in GraphBLAS.  All other
// objects exist to support them, and all the operations do their work on them.

// A sparse matrix is nrows-by-ncols and stored in a compressed sparse column
// form.  The row indices are kept sorted.  Also present is a list of pending
// tuples, held in (i,j,x) form in an unsorted format.  These are pending
// updates to the matrix, having been put there by the setElement method and/or
// assign operations.  The row and column indices of a matrix are of type
// GrB_Index, and they range from 0 to the dimesion minus 1.  That is, they are
// zero-based.

// Like all GraphBLAS objects, the GrB_Vector and GrB_Matrix are opaque to
// the user; their structure may change in future releases.

// GrB_Index: row or column index, or matrix dimension.  This typedef is used
// for row and column indices, or matrix and vector dimensions.

typedef uint64_t GrB_Index ;

// The GraphBLAS GrB_Matrix object; content not user-accessible

typedef struct
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Type type ;         // the type of each numerical entry
    int64_t nrows ;         // number of rows
    int64_t ncols ;         // number of columns
    int64_t nzmax ;         // size of i and x arrays
    int64_t *p ;            // column pointers, array of size ncols+1
    int64_t *i ;            // row indices, array of size nzmax
    void *x ;               // values, size nzmax; each size A->type->size
    bool p_shallow ;        // true if p is a shallow copy
    bool i_shallow ;        // true if i is a shallow copy
    bool x_shallow ;        // true if x is a shallow copy
    int64_t npending ;      // number of pending tuples to add to the matrix
    int64_t max_npending ;  // size of ipending, jpending, and xpending arrays
    bool sorted_pending ;   // true if pending tuples are in sorted order
    int64_t *ipending ;     // row indices of pending tuples
    int64_t *jpending ;     // col indices of pending tuples; NULL if ncols <= 1
    void *xpending ;        // values of pending tuples
    GrB_BinaryOp operator_pending ; // operator to assemble duplications
    int64_t nzombies ;      // number of zombines marked for deletion
    void *queue_next ;      // next matrix in the matrix queue
    void *queue_prev ;      // prev matrix in the matrix queue
    bool enqueued ;         // true if the matrix is in the queue
}
GB_Matrix_opaque ;          // CONTENT NOT USER-ACCESSIBLE

// The GrB_Matrix handle (user-accesible)
typedef GB_Matrix_opaque *GrB_Matrix ;

// The GraphBLAS GrB_Vector object; content not user-accessible.  The content
// is exactly the same as a GrB_Matrix (SuiteSparse:GraphBLAS requires these
// to objects to be identical in size and content).

typedef struct
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Type type ;         // the type of each numerical entry
    int64_t nrows ;         // number of rows
    int64_t ncols ;         // always 1
    int64_t nzmax ;         // size of i and x arrays
    int64_t *p ;            // column pointers, array of size ncols+1 == 2
    int64_t *i ;            // row indices, array of size nzmax
    void *x ;               // values, size nzmax; each size A->type->size
    bool p_shallow ;        // true if p is a shallow copy
    bool i_shallow ;        // true if i is a shallow copy
    bool x_shallow ;        // true if x is a shallow copy
    int64_t npending ;      // number of pending tuples to add to the matrix
    int64_t max_npending ;  // size of ipending, jpending, and xpending arrays
    bool sorted_pending ;   // true if pending tuples are in sorted order
    int64_t *ipending ;     // row indices of pending tuples
    int64_t *jpending ;     // always NULL
    void *xpending ;        // values of pending tuples
    GrB_BinaryOp operator_pending ; // operator to assemble duplications
    int64_t nzombies ;      // number of zombines marked for deletion
    void *queue_next ;      // next matrix in the matrix queue
    void *queue_prev ;      // prev matrix in the matrix queue
    bool enqueued ;         // true if the matrix is in the queue
}
GB_Vector_opaque ;          // CONTENT NOT USER-ACCESSIBLE

// The GrB_Vector handle (user-accesible)
typedef GB_Vector_opaque *GrB_Vector ;

//==============================================================================
//=== GraphBLAS Vector methods =================================================
//==============================================================================

// These methods create, free, copy, and clear a vector.  The size, nvals,
// and type methods return basic information about a vector.

GrB_Info GrB_Vector_new     // create a new vector with no entries
(
    GrB_Vector *v,          // handle of vector to create
    const GrB_Type type,    // type of vector to create
    const GrB_Index n       // vector dimension is n-by-1
) ;

GrB_Info GrB_Vector_dup     // make an exact copy of a vector
(
    GrB_Vector *w,          // handle of output vector to create
    const GrB_Vector u      // input vector to copy
) ;

GrB_Info GrB_Vector_clear   // clear a vector of all entries;
(                           // type and dimension remain unchanged.
    GrB_Vector v            // vector to clear
) ;

GrB_Info GrB_Vector_size    // get the dimension of a vector
(
    GrB_Index *n,           // vector dimension is n-by-1
    const GrB_Vector v      // vector to query
) ;

GrB_Info GrB_Vector_nvals   // get the number of entries in a vector
(
    GrB_Index *nvals,       // vector has nvals entries
    const GrB_Vector v      // vector to query
) ;

// SPEC: GxB_Vector_type is an extension to the spec

GrB_Info GxB_Vector_type    // get the type of a vector
(
    GrB_Type *type,         // returns the type of the vector
    const GrB_Vector v      // vector to query
) ;


//------------------------------------------------------------------------------
// GrB_Vector_build
//------------------------------------------------------------------------------

// GrB_Vector_build:  w = sparse (I,1,X) in MATLAB notation, but using any
// associative operator to assemble duplicate entries.

// Build a vector w from a set of (i,x) tuples.  The type and dimension of the
// vector is already defined in w (via GrB_Vector_new), which must initially
// have no entries.  I [0..nvals-1] is the list of row indices, and X
// [0..nvals-1] is the list of numerical values.  The kth tuple is (I[k],X[k]),
// and tuples can appear in any order.  Values are typecasted from X into the
// type of the dup operator, as needed (user-defined types cannot be cast).
// Duplicates are assembled together with the dup operator.  If two tuples
// (i,x1) and (i,x2) have the same row index, then w(i) = dup (x1,x2).  All
// three types of dup must be the same.  The types of C, X, and dup must be
// compatible.

// SPEC: extension: well-defined behavior of a non-associative dup operator.

// The GraphBLAS spec requires dup to be associative and does not define the
// order in which duplicates are assembled.  Currently this implementation
// assembles duplicates in the order they appear in I and X.  For example, if
// (i,x1), (i,x2), and (i,x3) appear in that order in I and X, then w(i) =
// dup(dup(x1,x2),x3).  This means that using the non-associative FIRST
// operator as dup means that w(i) is set equal to the first entry in the list,
// x1, and SECOND gives the last one, x3.  SuiteSparse:GraphBLAS guarantees
// this ordering.  However, per the spec, this order of assembly is not
// guaranteed in all implementations.  Thus dup must be associative and results
// are not guaranteed in all implementations if it is not.

GrB_Info GrB_Vector_build_BOOL      // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const bool *X,                  // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_INT8      // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const int8_t *X,                // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_UINT8     // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const uint8_t *X,               // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_INT16     // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const int16_t *X,               // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_UINT16    // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const uint16_t *X,              // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_INT32     // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const int32_t *X,               // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_UINT32    // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const uint32_t *X,              // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_INT64     // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const int64_t *X,               // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_UINT64    // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const uint64_t *X,              // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_FP32      // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const float *X,                 // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_FP64      // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const double *X,                // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Vector_build_UDT       // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const void *X,                  // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

// Type-generic version:  X can be a pointer to any supported C type or void *
// for a user-defined type.

/*

GrB_Info GrB_Vector_build           // build a vector from (I,X) tuples
(
    GrB_Vector w,                   // vector to build
    const GrB_Index *I,             // array of row indices of tuples
    const <type> *X,                // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

*/

#define GrB_Vector_build(w,I,X,nvals,dup)               \
    _Generic                                            \
    (                                                   \
        (X),                                            \
        const bool      *: GrB_Vector_build_BOOL   ,    \
              bool      *: GrB_Vector_build_BOOL   ,    \
        const int8_t    *: GrB_Vector_build_INT8   ,    \
              int8_t    *: GrB_Vector_build_INT8   ,    \
        const uint8_t   *: GrB_Vector_build_UINT8  ,    \
              uint8_t   *: GrB_Vector_build_UINT8  ,    \
        const int16_t   *: GrB_Vector_build_INT16  ,    \
              int16_t   *: GrB_Vector_build_INT16  ,    \
        const uint16_t  *: GrB_Vector_build_UINT16 ,    \
              uint16_t  *: GrB_Vector_build_UINT16 ,    \
        const int32_t   *: GrB_Vector_build_INT32  ,    \
              int32_t   *: GrB_Vector_build_INT32  ,    \
        const uint32_t  *: GrB_Vector_build_UINT32 ,    \
              uint32_t  *: GrB_Vector_build_UINT32 ,    \
        const int64_t   *: GrB_Vector_build_INT64  ,    \
              int64_t   *: GrB_Vector_build_INT64  ,    \
        const uint64_t  *: GrB_Vector_build_UINT64 ,    \
              uint64_t  *: GrB_Vector_build_UINT64 ,    \
        const float     *: GrB_Vector_build_FP32   ,    \
              float     *: GrB_Vector_build_FP32   ,    \
        const double    *: GrB_Vector_build_FP64   ,    \
              double    *: GrB_Vector_build_FP64   ,    \
        const void      *: GrB_Vector_build_UDT    ,    \
              void      *: GrB_Vector_build_UDT         \
    )                                                   \
    (w, I, ((const void *) (X)), nvals, dup)

//------------------------------------------------------------------------------
// GrB_Vector_setElement
//------------------------------------------------------------------------------

// Set a single scalar in a vector, w(i) = x, typecasting from the type of x to
// the type of w as needed.

GrB_Info GrB_Vector_setElement_BOOL     // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const bool x,                       // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_INT8     // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const int8_t x,                     // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_UINT8    // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const uint8_t x,                    // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_INT16    // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const int16_t x,                    // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_UINT16   // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const uint16_t x,                   // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_INT32    // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const int32_t x,                    // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_UINT32   // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const uint32_t x,                   // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_INT64    // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const int64_t x,                    // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_UINT64   // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const uint64_t x,                   // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_FP32     // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const float x,                      // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_FP64     // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const double x,                     // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

GrB_Info GrB_Vector_setElement_UDT      // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const void *x,                      // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

// Type-generic version:  x can be any supported C type or void * for a
// user-defined type.

/*

GrB_Info GrB_Vector_setElement          // w(i) = x
(
    GrB_Vector w,                       // vector to modify
    const <type> x,                     // scalar to assign to w(i)
    const GrB_Index i                   // row index
) ;

*/

#define GrB_Vector_setElement(w,x,i)                \
    _Generic                                        \
    (                                               \
        (x),                                        \
        const bool      : GrB_Vector_setElement_BOOL   ,  \
              bool      : GrB_Vector_setElement_BOOL   ,  \
        const int8_t    : GrB_Vector_setElement_INT8   ,  \
              int8_t    : GrB_Vector_setElement_INT8   ,  \
        const uint8_t   : GrB_Vector_setElement_UINT8  ,  \
              uint8_t   : GrB_Vector_setElement_UINT8  ,  \
        const int16_t   : GrB_Vector_setElement_INT16  ,  \
              int16_t   : GrB_Vector_setElement_INT16  ,  \
        const uint16_t  : GrB_Vector_setElement_UINT16 ,  \
              uint16_t  : GrB_Vector_setElement_UINT16 ,  \
        const int32_t   : GrB_Vector_setElement_INT32  ,  \
              int32_t   : GrB_Vector_setElement_INT32  ,  \
        const uint32_t  : GrB_Vector_setElement_UINT32 ,  \
              uint32_t  : GrB_Vector_setElement_UINT32 ,  \
        const int64_t   : GrB_Vector_setElement_INT64  ,  \
              int64_t   : GrB_Vector_setElement_INT64  ,  \
        const uint64_t  : GrB_Vector_setElement_UINT64 ,  \
              uint64_t  : GrB_Vector_setElement_UINT64 ,  \
        const float     : GrB_Vector_setElement_FP32   ,  \
              float     : GrB_Vector_setElement_FP32   ,  \
        const double    : GrB_Vector_setElement_FP64   ,  \
              double    : GrB_Vector_setElement_FP64   ,  \
        const void *    : GrB_Vector_setElement_UDT    ,  \
              void *    : GrB_Vector_setElement_UDT       \
    )                                               \
    (w, x, i)

//------------------------------------------------------------------------------
// GrB_Vector_extractElement
//------------------------------------------------------------------------------

// Extract a single entry from a vector, x = v(i), typecasting from the type of
// v to the type of x as needed.

// Returns GrB_SUCCESS if v(i) is present, and sets x to its value.
// Returns GrB_NO_VALUE if v(i) is not present, and x is unmodified.

GrB_Info GrB_Vector_extractElement_BOOL     // x = v(i)
(
    bool *x,                        // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_INT8     // x = v(i)
(
    int8_t *x,                      // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_UINT8    // x = v(i)
(
    uint8_t *x,                     // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_INT16    // x = v(i)
(
    int16_t *x,                     // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_UINT16   // x = v(i)
(
    uint16_t *x,                    // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_INT32    // x = v(i)
(
    int32_t *x,                     // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_UINT32   // x = v(i)
(
    uint32_t *x,                    // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_INT64    // x = v(i)
(
    int64_t *x,                     // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_UINT64   // x = v(i)
(
    uint64_t *x,                    // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_FP32     // x = v(i)
(
    float *x,                       // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_FP64     // x = v(i)
(
    double *x,                      // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

GrB_Info GrB_Vector_extractElement_UDT      // x = v(i)
(
    void *x,                        // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

// Type-generic version:  x can be a pointer to any supported C type or void *
// for a user-defined type.

/*

GrB_Info GrB_Vector_extractElement  // x = v(i)
(
    <type> *x,                      // scalar extracted
    const GrB_Vector v,             // vector to extract an entry from
    const GrB_Index i               // row index
) ;

*/

#define GrB_Vector_extractElement(x,v,i)                \
    _Generic                                            \
    (                                                   \
        (x),                                            \
        bool     *: GrB_Vector_extractElement_BOOL   ,  \
        int8_t   *: GrB_Vector_extractElement_INT8   ,  \
        uint8_t  *: GrB_Vector_extractElement_UINT8  ,  \
        int16_t  *: GrB_Vector_extractElement_INT16  ,  \
        uint16_t *: GrB_Vector_extractElement_UINT16 ,  \
        int32_t  *: GrB_Vector_extractElement_INT32  ,  \
        uint32_t *: GrB_Vector_extractElement_UINT32 ,  \
        int64_t  *: GrB_Vector_extractElement_INT64  ,  \
        uint64_t *: GrB_Vector_extractElement_UINT64 ,  \
        float    *: GrB_Vector_extractElement_FP32   ,  \
        double   *: GrB_Vector_extractElement_FP64   ,  \
        void     *: GrB_Vector_extractElement_UDT       \
    )                                                   \
    (x, v, i)

//------------------------------------------------------------------------------
// GrB_Vector_extractTuples
//------------------------------------------------------------------------------

// Extracts all tuples from a vector, like [I,~,X] = find (v) in MATLAB.  If
// any parameter I and/or X is NULL, then that component is not extracted.  The
// size of the I and X arrays (those that are not NULL) is given by nvals,
// which must be at least as large as GrB_nvals (&nvals, v).  The values in the
// vector are typecasted to the type of X, as needed.

// If any parameter I and/or X is NULL, that component is not extracted.  So to
// extract just the row indices, pass I as non-NULL, and X as NULL.  This is
// like [I,~,~] = find (v) in MATLAB.

GrB_Info GrB_Vector_extractTuples_BOOL      // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    bool *X,                    // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_INT8      // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    int8_t *X,                  // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_UINT8     // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    uint8_t *X,                 // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_INT16     // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    int16_t *X,                 // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_UINT16    // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    uint16_t *X,                // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_INT32     // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    int32_t *X,                 // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_UINT32    // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    uint32_t *X,                // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_INT64     // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    int64_t *X,                 // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_UINT64    // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    uint64_t *X,                // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_FP32      // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    float *X,                   // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_FP64      // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    double *X,                  // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

GrB_Info GrB_Vector_extractTuples_UDT       // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    void *X,                    // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

// Type-generic version:  X can be a pointer to any supported C type or void *
// for a user-defined type.

/*

GrB_Info GrB_Vector_extractTuples           // [I,~,X] = find (v)
(
    GrB_Index *I,               // array for returning row indices of tuples
    <type> *X,                  // array for returning values of tuples
    GrB_Index *nvals,           // I, X size on input; # tuples on output
    const GrB_Vector v          // vector to extract tuples from
) ;

*/

#define GrB_Vector_extractTuples(I,X,nvals,v)           \
    _Generic                                            \
    (                                                   \
        (X),                                            \
        bool     *: GrB_Vector_extractTuples_BOOL   ,   \
        int8_t   *: GrB_Vector_extractTuples_INT8   ,   \
        uint8_t  *: GrB_Vector_extractTuples_UINT8  ,   \
        int16_t  *: GrB_Vector_extractTuples_INT16  ,   \
        uint16_t *: GrB_Vector_extractTuples_UINT16 ,   \
        int32_t  *: GrB_Vector_extractTuples_INT32  ,   \
        uint32_t *: GrB_Vector_extractTuples_UINT32 ,   \
        int64_t  *: GrB_Vector_extractTuples_INT64  ,   \
        uint64_t *: GrB_Vector_extractTuples_UINT64 ,   \
        float    *: GrB_Vector_extractTuples_FP32   ,   \
        double   *: GrB_Vector_extractTuples_FP64   ,   \
        void     *: GrB_Vector_extractTuples_UDT        \
    )                                                   \
    (I, X, nvals, v)

//------------------------------------------------------------------------------
// GrB_Vector_free
//------------------------------------------------------------------------------

GrB_Info GrB_Vector_free    // free a vector
(
    GrB_Vector *v           // handle of vector to free
) ;

//==============================================================================
//=== GraphBLAS Matrix methods =================================================
//==============================================================================

// These methods create, free, copy, and clear a matrix.  The nrows, ncols,
// nvals, and type methods return basic information about a matrix.

GrB_Info GrB_Matrix_new     // create a new matrix with no entries
(
    GrB_Matrix *A,          // handle of matrix to create
    const GrB_Type type,    // type of matrix to create
    const GrB_Index nrows,  // matrix dimension is nrows-by-ncols
    const GrB_Index ncols
) ;

GrB_Info GrB_Matrix_dup     // make an exact copy of a matrix
(
    GrB_Matrix *C,          // handle of output matrix to create
    const GrB_Matrix A      // input matrix to copy
) ;

GrB_Info GrB_Matrix_clear   // clear a matrix of all entries;
(                           // type and dimensions remain unchanged
    GrB_Matrix A            // matrix to clear
) ;

GrB_Info GrB_Matrix_nrows   // get the number of rows of a matrix
(
    GrB_Index *nrows,       // matrix has nrows rows
    const GrB_Matrix A      // matrix to query
) ;

GrB_Info GrB_Matrix_ncols   // get the number of columns of a matrix
(
    GrB_Index *ncols,       // matrix has ncols columns
    const GrB_Matrix A      // matrix to query
) ;

GrB_Info GrB_Matrix_nvals   // get the number of entries in a matrix
(
    GrB_Index *nvals,       // matrix has nvals entries
    const GrB_Matrix A      // matrix to query
) ;

// SPEC: GxB_Matrix_type is an extension to the spec

GrB_Info GxB_Matrix_type    // get the type of a matrix
(
    GrB_Type *type,         // returns the type of the matrix
    const GrB_Matrix A      // matrix to query
) ;

//------------------------------------------------------------------------------
// GrB_Matrix_build
//------------------------------------------------------------------------------

// GrB_Matrix_build:  C = sparse (I,J,X) in MATLAB notation, but using any
// associative operator to assemble duplicate entries.

// Builds a matrix C from a set of (i,j,x) tuples.  The type and dimension of
// the matrix is already defined in C (via GrB_Matrix_new), which must
// initially have no entries.  I [0..nvals-1] is the list of row indices, J
// [0..nvals-1] is the list of column indices, and X [0..nvals-1] is the list
// of numerical values.  The kth triplet is (I[k],J[k],X[k]), and tuples can
// appear in any order.  Values are typecasted from X into the type of C, as
// needed (user-defined types cannot be cast).  Duplicates are assembled
// together with the dup operator.  If two tuples (i,j,x1) and (i,j,x2) have
// the same row index, then C(i,j) = dup(x1,x2).  All three types of dup must
// be the same, and dup, C, and X must be compatible.

// SPEC: extension: well-defined behavior of a non-associative dup operator.

// The dup operator must be associative in general, and the GraphBLAS spec
// states the order of assembly is not defined.  However, SuiteSparse:GraphBLAS
// does guarantee an ordering; see the description of GrB_Vector_build for more
// details.

GrB_Info GrB_Matrix_build_BOOL      // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const bool *X,                  // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_INT8      // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const int8_t *X,                // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_UINT8     // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const uint8_t *X,               // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_INT16     // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const int16_t *X,               // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_UINT16    // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const uint16_t *X,              // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_INT32     // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const int32_t *X,               // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_UINT32    // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const uint32_t *X,              // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_INT64     // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const int64_t *X,               // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_UINT64    // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const uint64_t *X,              // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_FP32      // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const float *X,                 // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_FP64      // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const double *X,                // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

GrB_Info GrB_Matrix_build_UDT       // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const void *X,                  // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

// Type-generic version:  X can be a pointer to any supported C type or void *
// for a user-defined type.

/*

GrB_Info GrB_Matrix_build           // build a matrix from (I,J,X) tuples
(
    GrB_Matrix C,                   // matrix to build
    const GrB_Index *I,             // array of row indices of tuples
    const GrB_Index *J,             // array of column indices of tuples
    const <type> *X,                // array of values of tuples
    const GrB_Index nvals,          // number of tuples
    const GrB_BinaryOp dup          // binary function to assemble duplicates
) ;

*/

#define GrB_Matrix_build(C,I,J,X,nvals,dup)             \
    _Generic                                            \
    (                                                   \
        (X),                                            \
        const bool      *: GrB_Matrix_build_BOOL   ,    \
              bool      *: GrB_Matrix_build_BOOL   ,    \
        const int8_t    *: GrB_Matrix_build_INT8   ,    \
              int8_t    *: GrB_Matrix_build_INT8   ,    \
        const uint8_t   *: GrB_Matrix_build_UINT8  ,    \
              uint8_t   *: GrB_Matrix_build_UINT8  ,    \
        const int16_t   *: GrB_Matrix_build_INT16  ,    \
              int16_t   *: GrB_Matrix_build_INT16  ,    \
        const uint16_t  *: GrB_Matrix_build_UINT16 ,    \
              uint16_t  *: GrB_Matrix_build_UINT16 ,    \
        const int32_t   *: GrB_Matrix_build_INT32  ,    \
              int32_t   *: GrB_Matrix_build_INT32  ,    \
        const uint32_t  *: GrB_Matrix_build_UINT32 ,    \
              uint32_t  *: GrB_Matrix_build_UINT32 ,    \
        const int64_t   *: GrB_Matrix_build_INT64  ,    \
              int64_t   *: GrB_Matrix_build_INT64  ,    \
        const uint64_t  *: GrB_Matrix_build_UINT64 ,    \
              uint64_t  *: GrB_Matrix_build_UINT64 ,    \
        const float     *: GrB_Matrix_build_FP32   ,    \
              float     *: GrB_Matrix_build_FP32   ,    \
        const double    *: GrB_Matrix_build_FP64   ,    \
              double    *: GrB_Matrix_build_FP64   ,    \
        const void      *: GrB_Matrix_build_UDT    ,    \
              void      *: GrB_Matrix_build_UDT         \
    )                                                   \
    (C, I, J, ((const void *) (X)), nvals, dup)

//------------------------------------------------------------------------------
// GrB_Matrix_setElement
//------------------------------------------------------------------------------

// Set a single entry in a matrix, C(i,j) = x in MATLAB notation, typecasting
// from the type of x to the type of C, as needed.

GrB_Info GrB_Matrix_setElement_BOOL     // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const bool x,                       // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_INT8     // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const int8_t x,                     // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_UINT8    // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const uint8_t x,                    // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_INT16    // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const int16_t x,                    // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_UINT16   // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const uint16_t x,                   // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_INT32    // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const int32_t x,                    // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_UINT32   // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const uint32_t x,                   // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_INT64    // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const int64_t x,                    // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_UINT64   // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const uint64_t x,                   // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_FP32     // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const float x,                      // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_FP64     // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const double x,                     // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_setElement_UDT      // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const void *x,                      // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

// Type-generic version:  x can be any supported C type or void * for a
// user-defined type.

/*

GrB_Info GrB_Matrix_setElement          // C (i,j) = x
(
    GrB_Matrix C,                       // matrix to modify
    const <type> x,                     // scalar to assign to C(i,j)
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

*/

#define GrB_Matrix_setElement(C,x,i,j)                    \
    _Generic                                              \
    (                                                     \
        (x),                                              \
        const bool      : GrB_Matrix_setElement_BOOL   ,  \
              bool      : GrB_Matrix_setElement_BOOL   ,  \
        const int8_t    : GrB_Matrix_setElement_INT8   ,  \
              int8_t    : GrB_Matrix_setElement_INT8   ,  \
        const uint8_t   : GrB_Matrix_setElement_UINT8  ,  \
              uint8_t   : GrB_Matrix_setElement_UINT8  ,  \
        const int16_t   : GrB_Matrix_setElement_INT16  ,  \
              int16_t   : GrB_Matrix_setElement_INT16  ,  \
        const uint16_t  : GrB_Matrix_setElement_UINT16 ,  \
              uint16_t  : GrB_Matrix_setElement_UINT16 ,  \
        const int32_t   : GrB_Matrix_setElement_INT32  ,  \
              int32_t   : GrB_Matrix_setElement_INT32  ,  \
        const uint32_t  : GrB_Matrix_setElement_UINT32 ,  \
              uint32_t  : GrB_Matrix_setElement_UINT32 ,  \
        const int64_t   : GrB_Matrix_setElement_INT64  ,  \
              int64_t   : GrB_Matrix_setElement_INT64  ,  \
        const uint64_t  : GrB_Matrix_setElement_UINT64 ,  \
              uint64_t  : GrB_Matrix_setElement_UINT64 ,  \
        const float     : GrB_Matrix_setElement_FP32   ,  \
              float     : GrB_Matrix_setElement_FP32   ,  \
        const double    : GrB_Matrix_setElement_FP64   ,  \
              double    : GrB_Matrix_setElement_FP64   ,  \
        const void *    : GrB_Matrix_setElement_UDT    ,  \
              void *    : GrB_Matrix_setElement_UDT       \
    )                                                     \
    (C, x, i, j)

//------------------------------------------------------------------------------
// GrB_Matrix_extractElement
//------------------------------------------------------------------------------

// Extract a single entry from a matrix, x = A(i,j), typecasting from the type
// of A to the type of x, as needed.

// Returns GrB_SUCCESS if A(i,j) is present, and sets x to its value.
// Returns GrB_NO_VALUE if A(i,j) is not present, and x is unmodified.

GrB_Info GrB_Matrix_extractElement_BOOL     // x = A(i,j)
(
    bool *x,                            // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_INT8     // x = A(i,j)
(
    int8_t *x,                          // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_UINT8    // x = A(i,j)
(
    uint8_t *x,                         // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_INT16    // x = A(i,j)
(
    int16_t *x,                         // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_UINT16   // x = A(i,j)
(
    uint16_t *x,                        // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_INT32    // x = A(i,j)
(
    int32_t *x,                         // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_UINT32   // x = A(i,j)
(
    uint32_t *x,                        // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_INT64    // x = A(i,j)
(
    int64_t *x,                         // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_UINT64   // x = A(i,j)
(
    uint64_t *x,                        // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_FP32     // x = A(i,j)
(
    float *x,                           // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_FP64     // x = A(i,j)
(
    double *x,                          // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

GrB_Info GrB_Matrix_extractElement_UDT      // x = A(i,j)
(
    void *x,                            // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

// Type-generic version:  x can be a pointer to any supported C type or void *
// for a user-defined type.

/*

GrB_Info GrB_Matrix_extractElement      // x = A(i,j)
(
    <type> *x,                          // extracted scalar
    const GrB_Matrix A,                 // matrix to extract a scalar from
    const GrB_Index i,                  // row index
    const GrB_Index j                   // column index
) ;

*/

#define GrB_Matrix_extractElement(x,A,i,j)              \
    _Generic                                            \
    (                                                   \
        (x),                                            \
        bool     *: GrB_Matrix_extractElement_BOOL   ,  \
        int8_t   *: GrB_Matrix_extractElement_INT8   ,  \
        uint8_t  *: GrB_Matrix_extractElement_UINT8  ,  \
        int16_t  *: GrB_Matrix_extractElement_INT16  ,  \
        uint16_t *: GrB_Matrix_extractElement_UINT16 ,  \
        int32_t  *: GrB_Matrix_extractElement_INT32  ,  \
        uint32_t *: GrB_Matrix_extractElement_UINT32 ,  \
        int64_t  *: GrB_Matrix_extractElement_INT64  ,  \
        uint64_t *: GrB_Matrix_extractElement_UINT64 ,  \
        float    *: GrB_Matrix_extractElement_FP32   ,  \
        double   *: GrB_Matrix_extractElement_FP64   ,  \
        void     *: GrB_Matrix_extractElement_UDT       \
    )                                                   \
    (x, A, i, j)

//------------------------------------------------------------------------------
// GrB_Matrix_extractTuples
//------------------------------------------------------------------------------

// Extracts all tuples from a matrix, like [I,J,X] = find (A) in MATLAB.  If
// any parameter I, J and/or X is NULL, then that component is not extracted.
// The size of the I, J, and X arrays (those that are not NULL) is given by
// nvals, which must be at least as large as GrB_nvals (&nvals, A).  The values
// in the matrix are typecasted to the type of X, as needed.

// If any parameter I, J, and/or X is NULL, that component is not extracted.
// So to extract just the row and col indices, pass I and J as non-NULL,
// and X as NULL.  This is like [I,J,~] = find (A).

GrB_Info GrB_Matrix_extractTuples_BOOL      // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    bool *X,                    // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_INT8      // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    int8_t *X,                  // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_UINT8     // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    uint8_t *X,                 // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_INT16     // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    int16_t *X,                 // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_UINT16    // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    uint16_t *X,                // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_INT32     // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    int32_t *X,                 // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_UINT32    // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    uint32_t *X,                // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_INT64     // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    int64_t *X,                 // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_UINT64    // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    uint64_t *X,                // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_FP32      // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    float *X,                   // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_FP64      // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    double *X,                  // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

GrB_Info GrB_Matrix_extractTuples_UDT       // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    void *X,                    // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

// Type-generic version:  X can be a pointer to any supported C type or void *
// for a user-defined type.

/*

GrB_Info GrB_Matrix_extractTuples           // [I,J,X] = find (A)
(
    GrB_Index *I,               // array for returning row indices of tuples
    GrB_Index *J,               // array for returning col indices of tuples
    <type> *X,                  // array for returning values of tuples
    GrB_Index *nvals,           // I,J,X size on input; # tuples on output
    const GrB_Matrix A          // matrix to extract tuples from
) ;

*/

#define GrB_Matrix_extractTuples(I,J,X,nvals,A)         \
    _Generic                                            \
    (                                                   \
        (X),                                            \
        bool     *: GrB_Matrix_extractTuples_BOOL   ,   \
        int8_t   *: GrB_Matrix_extractTuples_INT8   ,   \
        uint8_t  *: GrB_Matrix_extractTuples_UINT8  ,   \
        int16_t  *: GrB_Matrix_extractTuples_INT16  ,   \
        uint16_t *: GrB_Matrix_extractTuples_UINT16 ,   \
        int32_t  *: GrB_Matrix_extractTuples_INT32  ,   \
        uint32_t *: GrB_Matrix_extractTuples_UINT32 ,   \
        int64_t  *: GrB_Matrix_extractTuples_INT64  ,   \
        uint64_t *: GrB_Matrix_extractTuples_UINT64 ,   \
        float    *: GrB_Matrix_extractTuples_FP32   ,   \
        double   *: GrB_Matrix_extractTuples_FP64   ,   \
        void     *: GrB_Matrix_extractTuples_UDT        \
    )                                                   \
    (I, J, X, nvals, A)

//------------------------------------------------------------------------------
// GrB_Matrix_free
//------------------------------------------------------------------------------

GrB_Info GrB_Matrix_free    // free a matrix
(
    GrB_Matrix *A           // handle of matrix to free
) ;

//==============================================================================
//=== GraphBLAS Descriptor =====================================================
//==============================================================================

// The Descriptor is used to modify the behavior of GraphBLAS operations.
//
// OUTP: can be DEFAULT or REPLACE.  If REPLACE, then C is cleared after taking
//      part in the accum operation but before the mask.  In other words,
//      C<Mask> = accum (C,T) is split into Z = accum(C,T) ; C=0 ; C<Mask> = Z.
//
// MASK: can be DEFAULT or SCMP.  If DEFAULT, the mask is used normally,
//      where Mask(i,j)=1 means C(i,j) can be modified by C<Mask>=Z, and
//      Mask(i,j)=0 means it cannot be modified even if Z(i,j) is has been
//      computed and differs from C(i,j).  If SCMP, this is the same as
//      taking the logical complement of the Mask.
//
// INP0: can be DEFAULT or TRAN.  If DEFAULT, the first input is used as-is.
//      If TRAN, it is transposed.  Only matrices are transposed this way.
//
// INP1: the same as INP0 but for the second input

typedef enum
{
    GrB_OUTP,       // descriptor for output of a method
    GrB_MASK,       // descriptor for the mask input of a method
    GrB_INP0,       // descriptor for the first input of a method
    GrB_INP1        // descriptor for the second input of a method
}
GrB_Desc_Field ;

// SPEC: GxB_DEFAULT is an extension to the spec

typedef enum
{
    GxB_DEFAULT,    // default behavior of the method
    GrB_REPLACE,    // clear the output before assigning new values to it
    GrB_SCMP,       // use the structural complement of the input
    GrB_TRAN        // use the transpose of the input
}
GrB_Desc_Value ;

typedef struct
{
    int64_t magic ;         // for detecting uninitialized objects
    GrB_Desc_Value out ;    // output descriptor
    GrB_Desc_Value mask ;   // mask descriptor
    GrB_Desc_Value in0 ;    // first input descriptor (A for C=A*B, for example)
    GrB_Desc_Value in1 ;    // second input descriptor (B for C=A*B)
}
GB_Descriptor_opaque ;      // CONTENT NOT USER-ACCESSIBLE

// The GrB_Descriptor handle (user-accesible)
typedef GB_Descriptor_opaque *GrB_Descriptor ;

GrB_Info GrB_Descriptor_new     // create a new descriptor
(
    GrB_Descriptor *descriptor  // handle of descriptor to create
) ;

GrB_Info GrB_Descriptor_set     // set a parameter in a descriptor
(
    GrB_Descriptor desc,        // descriptor to modify
    const GrB_Desc_Field field, // parameter to change
    const GrB_Desc_Value val    // value to change it to
) ;

// SPEC: GxB_Descriptor_get is an extension to the spec

GrB_Info GxB_Descriptor_get     // get a parameter from a descriptor
(
    GrB_Desc_Value *val,        // value of the parameter
    const GrB_Descriptor desc,  // descriptor to query; NULL means defaults
    const GrB_Desc_Field field  // parameter to query
) ;

GrB_Info GrB_Descriptor_free    // free a descriptor
(
    GrB_Descriptor *descriptor  // handle of descriptor to free
) ;

//==============================================================================
//=== GrB_free =================================================================
//==============================================================================

// GrB_free: free a GraphBLAS object.  Each GraphBLAS object has a specific
// GrB_*_new and GrB_*_free method.  There is no generic GrB_new, but the
// generic GrB_free method can free any GraphBLAS object.  It is safe to free
// an object twice, and it is also safe to (attempt to) free a built-in object.
// In that case, GrB_free silently does nothing and returns GrB_SUCCESS.  By
// the GraphBLAS spec, GrB_*_free functions can return GrB_SUCCESS or
// GrB_PANIC; in this implementation they never panic.

#define GrB_free(object)                         \
    _Generic                                     \
    (                                            \
        (object),                                \
        GrB_Type       *: GrB_Type_free       ,  \
        GrB_UnaryOp    *: GrB_UnaryOp_free    ,  \
        GrB_BinaryOp   *: GrB_BinaryOp_free   ,  \
        GxB_SelectOp   *: GxB_SelectOp_free   ,  \
        GrB_Monoid     *: GrB_Monoid_free     ,  \
        GrB_Semiring   *: GrB_Semiring_free   ,  \
        GrB_Vector     *: GrB_Vector_free     ,  \
        GrB_Matrix     *: GrB_Matrix_free     ,  \
        GrB_Descriptor *: GrB_Descriptor_free    \
    )                                            \
    (object)

//==============================================================================
//=== GraphBLAS operations =====================================================
//==============================================================================

// Each GraphBLAS operation can be modified by an optional Mask, an optional
// accum operator, and a descriptor.

// The primary computation of an operation computes a matrix or vector T.  If
// accum is NULL, Z=T.  Otherwise, Z=accum(C,T) is computed, where accum is a
// binary operator applied in an element-wise add manner.  Next, C is
// optionally cleared if the OUTP:REPLACE descriptor is enabled.  Finally,
// C<Mask>=Z is computed.  If there is no Mask, C=Z, or if an empty Mask
// (Mask==NULL) is complemented via the descriptor, C is not modified at all.
// Otherwise C(Mask)=Z(Mask) is computed using MATLAB-style logical index, if
// the Mask is not complemented.  Otherwise C(~Mask)=Z(~Mask) is computed.
// This description is terse; see the User Guide for more details.

// GrB_NULL is used for the accum argument when no accum operation is desired,
// for the Mask argument when no Mask is desired, and for the descriptor
// argument when the default descriptor is desired.

#define GrB_NULL NULL

// An object that has been freed is a GrB_INVALID_HANDLE, a NULL pointer.

#define GrB_INVALID_HANDLE NULL

//------------------------------------------------------------------------------
// matrix and vector multiplication over a semiring
//------------------------------------------------------------------------------

// Each of these methods compute a matrix multiplication over a semiring.  The
// inputs are typecasted into the inputs of the semiring's multiply operator.
// The result T=A*B has the type of the multiplier output, which is also the 3
// types of the 'add' operator.  The 'add' operator is a commutatitive and
// associative monoid.

GrB_Info GrB_mxm                    // C<Mask> = accum (C, A*B)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Semiring semiring,    // defines '+' and '*' for A*B
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Matrix B,             // second input: matrix B
    const GrB_Descriptor desc       // descriptor for C, Mask, A, and B
) ;

GrB_Info GrB_vxm                    // w'<Mask> = accum (w, u'*A)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_Semiring semiring,    // defines '+' and '*' for u'*A
    const GrB_Vector u,             // first input:  vector u
    const GrB_Matrix A,             // second input: matrix A
    const GrB_Descriptor desc       // descriptor for w, mask, and A
) ;

GrB_Info GrB_mxv                    // w<Mask> = accum (w, A*u)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_Semiring semiring,    // defines '+' and '*' for A*B
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Vector u,             // second input: vector u
    const GrB_Descriptor desc       // descriptor for w, mask, and A
) ;

//------------------------------------------------------------------------------
// element-wise matrix and vector operations: using set intersection
//------------------------------------------------------------------------------

// GrB_eWiseMult computes C<Mask> = accum (C, A.*B), where ".*" is MATLAB
// notation, and where pairs of elements in two matrices (or vectors) are
// pairwise "multiplied" with C(i,j) = mult (A(i,j),B(i,j)).  The
// "multiplication" operator can be any binary operator.  This is not matrix
// multiplication in the conventional linear algebra sense; see GrB_mxm and
// related methods for that operation.  The pattern of the result T=A.*B is the
// set intersection (not union) of A and B.  Entries outside of the
// intersection are not computed.  This is primary difference with
// GrB_eWiseAdd.

// The input matrices A and/or B may be transposed first, via the descriptor.

// For a semiring, the mult operator is the semiring's multiply operator; note
// that this differs from the eWiseAdd methods which use the semiring's add
// operator instead. For a monoid, the mult operator is the monoid operator.

GrB_Info GrB_eWiseMult_Vector_Semiring       // w<Mask> = accum (w, u.*v)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_Semiring semiring,    // defines '.*' for t=u.*v
    const GrB_Vector u,             // first input:  vector u
    const GrB_Vector v,             // second input: vector v
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_eWiseMult_Vector_Monoid         // w<Mask> = accum (w, u.*v)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_Monoid monoid,        // defines '.*' for t=u.*v
    const GrB_Vector u,             // first input:  vector u
    const GrB_Vector v,             // second input: vector v
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_eWiseMult_Vector_BinaryOp       // w<Mask> = accum (w, u.*v)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_BinaryOp mult,        // defines '.*' for t=u.*v
    const GrB_Vector u,             // first input:  vector u
    const GrB_Vector v,             // second input: vector v
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_eWiseMult_Matrix_Semiring       // C<Mask> = accum (C, A.*B)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Semiring semiring,    // defines '.*' for T=A.*B
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Matrix B,             // second input: matrix B
    const GrB_Descriptor desc       // descriptor for C, Mask, A, and B
) ;

GrB_Info GrB_eWiseMult_Matrix_Monoid         // C<Mask> = accum (C, A.*B)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Monoid monoid,        // defines '.*' for T=A.*B
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Matrix B,             // second input: matrix B
    const GrB_Descriptor desc       // descriptor for C, Mask, A, and B
) ;

GrB_Info GrB_eWiseMult_Matrix_BinaryOp       // C<Mask> = accum (C, A.*B)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_BinaryOp mult,        // defines '.*' for T=A.*B
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Matrix B,             // second input: matrix B
    const GrB_Descriptor desc       // descriptor for C, Mask, A, and B
) ;

// All 6 of the above type-specific functions are captured in a single
// type-generic function, GrB_eWiseMult:

#define GrB_eWiseMult(C,Mask,accum,op,A,B,desc)                         \
    _Generic                                                            \
    (                                                                   \
        (C),                                                            \
        GrB_Matrix :                                                    \
            _Generic                                                    \
            (                                                           \
                (op),                                                   \
                const GrB_Semiring : GrB_eWiseMult_Matrix_Semiring ,    \
                      GrB_Semiring : GrB_eWiseMult_Matrix_Semiring ,    \
                const GrB_Monoid   : GrB_eWiseMult_Matrix_Monoid   ,    \
                      GrB_Monoid   : GrB_eWiseMult_Matrix_Monoid   ,    \
                const GrB_BinaryOp : GrB_eWiseMult_Matrix_BinaryOp ,    \
                      GrB_BinaryOp : GrB_eWiseMult_Matrix_BinaryOp      \
            ),                                                          \
        GrB_Vector :                                                    \
            _Generic                                                    \
            (                                                           \
                (op),                                                   \
                const GrB_Semiring : GrB_eWiseMult_Vector_Semiring ,    \
                      GrB_Semiring : GrB_eWiseMult_Vector_Semiring ,    \
                const GrB_Monoid   : GrB_eWiseMult_Vector_Monoid   ,    \
                      GrB_Monoid   : GrB_eWiseMult_Vector_Monoid   ,    \
                const GrB_BinaryOp : GrB_eWiseMult_Vector_BinaryOp ,    \
                      GrB_BinaryOp : GrB_eWiseMult_Vector_BinaryOp      \
            )                                                           \
    )                                                                   \
    (C, Mask, accum, op, A, B, desc)

//------------------------------------------------------------------------------
// element-wise matrix and vector operations: using set union
//------------------------------------------------------------------------------

// GrB_eWiseAdd computes C<Mask> = accum (C, A+B), where pairs of elements in
// two matrices (or two vectors) are pairwise "added".  The "add" operator can
// be any binary operator.  With the plus operator, this is the same matrix
// addition in conventional linear algebra.  The pattern of the result T=A+B is
// the set union (not intersection) of A and B.  Entries outside of the union
// are not computed.  That is, if both A(i,j) and B(i,j) are present in the
// pattern of A and B, then T(i,j) = A(i,j) "+" B(i,j).  If only A(i,j) is
// present then T(i,j) = A (i,j) and the "+" operator is not used.  Likewise,
// if only B(i,j) is in the pattern of B but A(i,j) is not in the pattern of A,
// then T(i,j) = B(i,j).  This is primary difference between GrB_eWiseAdd and
// GrB_eWiseMult; the same set of binary operators can be used in both methods,
// and the action they take on entries in the intersection of the pattern of A
// and B is identical.

// The input matrices A and/or B may be transposed first, via the descriptor.

// For a semiring, the mult operator is the semiring's add operator; note that
// this differs from the eWiseMult methods which use the semiring's multiply
// operator instead. For a monoid, the mult operator is the monoid operator.

GrB_Info GrB_eWiseAdd_Vector_Semiring       // w<Mask> = accum (w, u+v)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_Semiring semiring,    // defines '+' for t=u+v
    const GrB_Vector u,             // first input:  vector u
    const GrB_Vector v,             // second input: vector v
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_eWiseAdd_Vector_Monoid         // w<Mask> = accum (w, u+v)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_Monoid monoid,        // defines '+' for t=u+v
    const GrB_Vector u,             // first input:  vector u
    const GrB_Vector v,             // second input: vector v
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_eWiseAdd_Vector_BinaryOp       // w<Mask> = accum (w, u+v)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_BinaryOp add,         // defines '+' for t=u+v
    const GrB_Vector u,             // first input:  vector u
    const GrB_Vector v,             // second input: vector v
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_eWiseAdd_Matrix_Semiring       // C<Mask> = accum (C, A+B)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Semiring semiring,    // defines '+' for T=A+B
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Matrix B,             // second input: matrix B
    const GrB_Descriptor desc       // descriptor for C, Mask, A, and B
) ;

GrB_Info GrB_eWiseAdd_Matrix_Monoid         // C<Mask> = accum (C, A+B)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Monoid monoid,        // defines '+' for T=A+B
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Matrix B,             // second input: matrix B
    const GrB_Descriptor desc       // descriptor for C, Mask, A, and B
) ;

GrB_Info GrB_eWiseAdd_Matrix_BinaryOp       // C<Mask> = accum (C, A+B)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_BinaryOp add,         // defines '+' for T=A+B
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Matrix B,             // second input: matrix B
    const GrB_Descriptor desc       // descriptor for C, Mask, A, and B
) ;

#define GrB_eWiseAdd(C,Mask,accum,op,A,B,desc)                          \
    _Generic                                                            \
    (                                                                   \
        (C),                                                            \
        GrB_Matrix :                                                    \
            _Generic                                                    \
            (                                                           \
                (op),                                                   \
                const GrB_Semiring : GrB_eWiseAdd_Matrix_Semiring ,     \
                      GrB_Semiring : GrB_eWiseAdd_Matrix_Semiring ,     \
                const GrB_Monoid   : GrB_eWiseAdd_Matrix_Monoid   ,     \
                      GrB_Monoid   : GrB_eWiseAdd_Matrix_Monoid   ,     \
                const GrB_BinaryOp : GrB_eWiseAdd_Matrix_BinaryOp ,     \
                      GrB_BinaryOp : GrB_eWiseAdd_Matrix_BinaryOp       \
            ),                                                          \
        GrB_Vector :                                                    \
            _Generic                                                    \
            (                                                           \
                (op),                                                   \
                const GrB_Semiring : GrB_eWiseAdd_Vector_Semiring ,     \
                      GrB_Semiring : GrB_eWiseAdd_Vector_Semiring ,     \
                const GrB_Monoid   : GrB_eWiseAdd_Vector_Monoid   ,     \
                      GrB_Monoid   : GrB_eWiseAdd_Vector_Monoid   ,     \
                const GrB_BinaryOp : GrB_eWiseAdd_Vector_BinaryOp ,     \
                      GrB_BinaryOp : GrB_eWiseAdd_Vector_BinaryOp       \
            )                                                           \
    )                                                                   \
    (C, Mask, accum, op, A, B, desc)

//------------------------------------------------------------------------------
// matrix and vector extract
//------------------------------------------------------------------------------

// Extract entries from a matrix or vector; T = A(I,J) in MATLAB notation.
// This (like most GraphBLAS methods) is then followed by C<Mask>=accum(C,T).

// The input matrix A may be transposed first, via the descriptor.

// To extract all rows of a matrix or vector, as in A (:,J) in MATLAB, use
// I=GrB_ALL as the input argument.  For all columns of a matrix, use
// J=GrB_ALL.  GrB_ALL is a predefined pointer that is not NULL so that
// out-of-memory conditions can be (I=NULL) distinguished from a request for
// all rows (I=GrB_ALL).  The pointer GrB_ALL should never dereferenced, and it
// must not be freed or modified.

// Each of these can be used with their generic name, GrB_extract.

extern const uint64_t *GrB_ALL ;

GrB_Info GrB_Vector_extract         // w<mask> = accum (w, u(I))
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_Vector u,             // first input:  vector u
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Matrix_extract         // C<Mask> = accum (C, A(I,J))
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C, Mask, and A
) ;

GrB_Info GrB_Col_extract            // w<mask> = accum (w, A(I,j))
(
    GrB_Vector w,                   // input/output matrix for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index j,              // column index
    const GrB_Descriptor desc       // descriptor for w, mask, and A
) ;

//------------------------------------------------------------------------------
// GrB_extract: generic matrix/vector extraction
//------------------------------------------------------------------------------

// GrB_extract is a generic interface to the following functions:

// GrB_Vector_extract (w,mask,acc,u,I,ni,d)      // w<m>    = acc (w, u(I))
// GrB_Col_extract    (w,mask,acc,A,I,ni,j,d)    // w<m>    = acc (w, A(I,j))
// GrB_Matrix_extract (C,Mask,acc,A,I,ni,J,nj,d) // C<Mask> = acc (C, A(I,J))

#define GrB_extract(arg1,Mask,accum,arg4,...) \
    _Generic                                                \
    (                                                       \
        (arg1),                                             \
        GrB_Vector :                                        \
            _Generic                                        \
            (                                               \
                (arg4),                                     \
                const GrB_Vector : GrB_Vector_extract ,     \
                      GrB_Vector : GrB_Vector_extract ,     \
                const GrB_Matrix : GrB_Col_extract    ,     \
                      GrB_Matrix : GrB_Col_extract          \
            ),                                              \
        GrB_Matrix : GrB_Matrix_extract                     \
    )                                                       \
    (arg1, Mask, accum, arg4, __VA_ARGS__)

//------------------------------------------------------------------------------
// matrix and vector subassign: C(I,J)<Mask> = accum (C(I,J), A)
//------------------------------------------------------------------------------

// Assign entries in a matrix or vector; C(I,J) = A in MATLAB notation.
// Each of these can be used with their generic name, GxB_subassign.

// SPEC: The GxB_*_subassign functions are extensions to the spec.

// Each GxB_subassign function is very similar to its corresponding GrB_assign
// function in the spec, but they differ in three ways:

// (1) the mask in the GxB_subassign functions has the same dimensions as
//      w(I) for vectors and C(I,J) for matrices.  In GrB_assign, the mask is
//      the same size as w or C, respectively (except for GrB_Row_asssign and
//      GrB_Col_assign, in which case the mask is the same size as a row or
//      column of C, respectively).  The two masks are related.  If M is the
//      mask for GrB_assign, then M(I,J) is the mask for GxB_subassign.  If
//      there is no mask, or if I and J are both GrB_ALL, then the two masks
//      are the same.

//      For GrB_Row_assign and GrB_Col_assign, the mask vector is the same
//      size as a row or column of C, respectively.  For the corresponding
//      GxB_Row_subassign and GxB_Col_subassign operations, the mask is the
//      same size as the subrow C(i,J) or subcolumn C(I,j), respectively.

// (2) They differ in how C is affected in areas outside the C(I,J) submatrix.
//      In GxB_subassign, C(I,J) is the only part of C that can be modified,
//      and no part of C outside the submatrix is ever modified.  In
//      GrB_assign, it is possible to modify C outside the submatrix, but only
//      in one specific manner.  Suppose the mask M is present (or, suppose it
//      is not present but GrB_SCMP is true).  After (optionally) complementing
//      the mask, the value of M(i,j) can be 0 for some entry outside the
//      C(I,J) submatrix.  If the GrB_REPLACE descriptor is true, the
//      GrB_assign deletes this entry.  This case does not occur if GrB_REPLACE
//      is false.  With GrB_assign, it is not possible to change entries
//      outside the submatrix C(I,J), except to delete them in this
//      circumstance.

// (3) They differ in how duplicate indices are treated in I and J.  For both
//      assign and subassign operations, results are not defined for
//      GrB_Matrix_*assign, GrB_Vector_*assign, GrB_Row_*assign, and
//      GrB_Col_*assign when duplicate indices appear in I and J.  The scalar
//      expansion operations, GrB_*_assign_TYPE, are well-defined if duplicate
//      indices appear (the results are the same as if duplicates are removed
//      first from I and J).  However, the subassign scalar expansion
//      operations, GxB_*_subassign_TYPE are not well-defined if duplicate
//      indices appear in I and J.

// GxB_subassign and GrB_assign are identical if GrB_REPLACE is set to its
// default value of false, or if the masks happen to be the same.  The two
// masks can be the same in two cases:  either there is no mask (and GrB_SCMP
// is false), or I and J are both GrB_ALL.  In this case, the two algorithms
// are identical and have the same performance.

// GxB_subassign is much faster than GrB_assign, when the latter must examine
// the entire matrix C to delete entries (when GrB__REPLACE is true), and it
// must deal with a much larger Mask matrix.  However, both methods have
// specific uses.  Consider using C(I,J)+=F for many submatrices F (for
// example, when assembling a finite-element matrix).  If the Mask is meant as
// a specification for which entries of C should appear in the final result,
// then use GrB_assign.  If the Mask is meant to control which entries of the
// submatrix C(I,J) are modified by the finite-element F, then use
// GxB_subassign.  This is particularly useful is the Mask is a "template" that
// follows along with the finite-element F, independent of where it is applied
// C.  Using GrB_assign would be very difficult in this case since a new Mask,
// the same size as C, would need to be constructed for each finite-element F.

// In GraphBLAS notation, the two methods can be described as follows:

// matrix and vector subassign: C(I,J)<Mask> = accum (C(I,J), A)
// matrix and vector    assign: C<Mask>(I,J) = accum (C(I,J), A)

// This notation does not include the details of the GrB_SCMP and GrB_REPLACE
// descriptors, but it does illustrate the difference in the Mask.  In the
// subassign, Mask is the same size as C(I,J) and A.  If I[0]=i and J[0]=j,
// Then Mask(0,0) controls how C(i,j) is modified by the subassign, from the
// value A(0,0).  In the assign, Mask is the same size as C, and Mask(i,j)
// controls how C(i,j) is modified.

// The GxB_subassign and GrB_assign functions have the same signatures; they
// differ only in how they consider the Mask and the GrB_REPLACE descriptor,
// and in how duplicate indices are treated for scalar expansion.

GrB_Info GxB_Vector_subassign       // w(I)<mask> = accum (w(I),u)
(
    GrB_Vector w,                   // input/output matrix for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w(I),t)
    const GrB_Vector u,             // first input:  vector u
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Matrix_subassign       // C(I,J)<Mask> = accum (C(I,J),A)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),T)
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J), Mask, and A
) ;

GrB_Info GxB_Col_subassign          // C(I,j)<mask> = accum (C(I,j),u)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Vector mask,          // optional mask for C(I,j), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(C(I,j),t)
    const GrB_Vector u,             // input vector
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index j,              // column index
    const GrB_Descriptor desc       // descriptor for C(I,j) and mask
) ;

GrB_Info GxB_Row_subassign          // C(i,J)<mask'> = accum (C(i,J),u')
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Vector mask,          // optional mask for C(i,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(C(i,J),t)
    const GrB_Vector u,             // input vector
    const GrB_Index i,              // row index
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(i,J) and mask
) ;

//------------------------------------------------------------------------------
// GxB_Vector_subassign_[SCALAR]:  scalar expansion assignment to subvector
//------------------------------------------------------------------------------

// Assigns a single scalar to a subvector, w(I)<mask> = accum(w(I),x).  The
// scalar x is implicitly expanded into a vector u of size ni-by-1, with each
// entry in u equal to x, and then w(I)<mask> = accum(w(I),u) is done.

// Each of these can be used with their generic name, GxB_subassign.

GrB_Info GxB_Vector_subassign_BOOL  // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w(I),x)
    const bool x,                   // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_INT8  // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const int8_t x,                 // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_UINT8 // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const uint8_t x,                // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_INT16 // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const int16_t x,                // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_UINT16   // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const uint16_t x,               // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_INT32    // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const int32_t x,                // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_UINT32   // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const uint32_t x,               // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_INT64    // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const int64_t x,                // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_UINT64   // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const uint64_t x,               // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_FP32     // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const float x,                  // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_FP64     // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const double x,                 // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

GrB_Info GxB_Vector_subassign_UDT      // w(I)<mask> = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w(I), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const void *x,                  // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w(I) and mask
) ;

//------------------------------------------------------------------------------
// GxB_Matrix_subassign_[SCALAR]:  scalar expansion assignment to submatrix
//------------------------------------------------------------------------------

// Assigns a single scalar to a submatrix, C(I,J)<Mask> = accum(C(I,J),x).  The
// scalar x is implicitly expanded into a matrix A of size ni-by-nj, with each
// entry in A equal to x, and then C(I,J)<Mask> = accum(C(I,J),A) is done.

// Each of these can be used with their generic name, GxB_subassign.

GrB_Info GxB_Matrix_subassign_BOOL  // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const bool x,                   // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_INT8  // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const int8_t x,                 // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_UINT8 // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const uint8_t x,                // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_INT16 // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const int16_t x,                // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_UINT16   // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const uint16_t x,               // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_INT32    // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const int32_t x,                // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_UINT32   // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const uint32_t x,               // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_INT64    // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const int64_t x,                // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_UINT64   // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const uint64_t x,               // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_FP32     // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const float x,                  // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_FP64     // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const double x,                 // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

GrB_Info GxB_Matrix_subassign_UDT      // C(I,J)<Mask> = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C(I,J), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const void *x,                  // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(I,J) and Mask
) ;

//------------------------------------------------------------------------------
// GxB_subassign: generic submatrix/subvector assignment
//------------------------------------------------------------------------------

// GxB_subassign is a generic function that provides access to all specific
// GxB_*_subassign* functions:

// GxB_Vector_subassign   (w,mask,acc,u,I,ni,d)     // w(I)<mask>   =acc(w(I),u)
// GxB_Matrix_subassign   (C,Mask,acc,A,I,ni,J,nj,d)// C(I,J)<Mask> =acc(C(I,J),A)
// GxB_Col_subassign      (C,mask,acc,u,I,ni,j,d)   // C(I,j)<mask> =acc(C(I,j),u)
// GxB_Row_subassign      (C,mask,acc,u,i,J,nj,d)   // C(i,J)<mask'>=acc(C(i,J),u')
// GxB_Vector_subassign_T (w,mask,acc,x,I,ni,d)     // w(I)<mask>   =acc(w(I),x)
// GxB_Matrix_subassign_T (C,Mask,acc,x,I,ni,J,nj,d)// C(I,J)<Mask> =acc(C(I,J),x)

#define GxB_subassign(arg1,Mask,accum,arg4,arg5,...)               \
    _Generic                                                    \
    (                                                           \
        (arg1),                                                 \
        GrB_Vector :                                            \
            _Generic                                            \
            (                                                   \
                (arg4),                                         \
                const bool       : GxB_Vector_subassign_BOOL   ,   \
                      bool       : GxB_Vector_subassign_BOOL   ,   \
                const int8_t     : GxB_Vector_subassign_INT8   ,   \
                      int8_t     : GxB_Vector_subassign_INT8   ,   \
                const uint8_t    : GxB_Vector_subassign_UINT8  ,   \
                      uint8_t    : GxB_Vector_subassign_UINT8  ,   \
                const int16_t    : GxB_Vector_subassign_INT16  ,   \
                      int16_t    : GxB_Vector_subassign_INT16  ,   \
                const uint16_t   : GxB_Vector_subassign_UINT16 ,   \
                      uint16_t   : GxB_Vector_subassign_UINT16 ,   \
                const int32_t    : GxB_Vector_subassign_INT32  ,   \
                      int32_t    : GxB_Vector_subassign_INT32  ,   \
                const uint32_t   : GxB_Vector_subassign_UINT32 ,   \
                      uint32_t   : GxB_Vector_subassign_UINT32 ,   \
                const int64_t    : GxB_Vector_subassign_INT64  ,   \
                      int64_t    : GxB_Vector_subassign_INT64  ,   \
                const uint64_t   : GxB_Vector_subassign_UINT64 ,   \
                      uint64_t   : GxB_Vector_subassign_UINT64 ,   \
                const float      : GxB_Vector_subassign_FP32   ,   \
                      float      : GxB_Vector_subassign_FP32   ,   \
                const double     : GxB_Vector_subassign_FP64   ,   \
                      double     : GxB_Vector_subassign_FP64   ,   \
                const void *     : GxB_Vector_subassign_UDT    ,   \
                      void *     : GxB_Vector_subassign_UDT    ,   \
                default          : GxB_Vector_subassign            \
            ),                                                  \
        default :                                               \
            _Generic                                            \
            (                                                   \
                (arg4),                                         \
                const bool       : GxB_Matrix_subassign_BOOL   ,   \
                      bool       : GxB_Matrix_subassign_BOOL   ,   \
                const int8_t     : GxB_Matrix_subassign_INT8   ,   \
                      int8_t     : GxB_Matrix_subassign_INT8   ,   \
                const uint8_t    : GxB_Matrix_subassign_UINT8  ,   \
                      uint8_t    : GxB_Matrix_subassign_UINT8  ,   \
                const int16_t    : GxB_Matrix_subassign_INT16  ,   \
                      int16_t    : GxB_Matrix_subassign_INT16  ,   \
                const uint16_t   : GxB_Matrix_subassign_UINT16 ,   \
                      uint16_t   : GxB_Matrix_subassign_UINT16 ,   \
                const int32_t    : GxB_Matrix_subassign_INT32  ,   \
                      int32_t    : GxB_Matrix_subassign_INT32  ,   \
                const uint32_t   : GxB_Matrix_subassign_UINT32 ,   \
                      uint32_t   : GxB_Matrix_subassign_UINT32 ,   \
                const int64_t    : GxB_Matrix_subassign_INT64  ,   \
                      int64_t    : GxB_Matrix_subassign_INT64  ,   \
                const uint64_t   : GxB_Matrix_subassign_UINT64 ,   \
                      uint64_t   : GxB_Matrix_subassign_UINT64 ,   \
                const float      : GxB_Matrix_subassign_FP32   ,   \
                      float      : GxB_Matrix_subassign_FP32   ,   \
                const double     : GxB_Matrix_subassign_FP64   ,   \
                      double     : GxB_Matrix_subassign_FP64   ,   \
                const void *     : GxB_Matrix_subassign_UDT    ,   \
                      void *     : GxB_Matrix_subassign_UDT    ,   \
                const GrB_Vector :                              \
                    _Generic                                    \
                    (                                           \
                        (arg5),                                 \
                        const GrB_Index *: GxB_Col_subassign ,     \
                              GrB_Index *: GxB_Col_subassign ,     \
                        default          : GxB_Row_subassign       \
                    ),                                          \
                GrB_Vector :                                    \
                    _Generic                                    \
                    (                                           \
                        (arg5),                                 \
                        const GrB_Index *: GxB_Col_subassign ,     \
                              GrB_Index *: GxB_Col_subassign ,     \
                        default          : GxB_Row_subassign       \
                    ),                                          \
                default    : GxB_Matrix_subassign                  \
            )                                                   \
    )                                                           \
    (arg1, Mask, accum, arg4, arg5, __VA_ARGS__)

//------------------------------------------------------------------------------
// matrix and vector assign: C<Mask>(I,J) = accum (C(I,J), A)
//------------------------------------------------------------------------------

// Assign entries in a matrix or vector; C(I,J) = A in MATLAB notation.
// Each of these can be used with their generic name, GrB_assign.

GrB_Info GrB_Vector_assign          // w<mask>(I) = accum (w(I),u)
(
    GrB_Vector w,                   // input/output matrix for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w(I),t)
    const GrB_Vector u,             // first input:  vector u
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Matrix_assign          // C<Mask>(I,J) = accum (C(I,J),A)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),T)
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C, Mask, and A
) ;

GrB_Info GrB_Col_assign             // C<mask>(I,j) = accum (C(I,j),u)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Vector mask,          // optional mask for C(:,j), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(C(I,j),t)
    const GrB_Vector u,             // input vector
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index j,              // column index
    const GrB_Descriptor desc       // descriptor for C(:,j) and mask
) ;

GrB_Info GrB_Row_assign             // C<mask'>(i,J) = accum (C(i,J),u')
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Vector mask,          // optional mask for C(i,:), unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(C(i,J),t)
    const GrB_Vector u,             // input vector
    const GrB_Index i,              // row index
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C(i,:) and mask
) ;

//------------------------------------------------------------------------------
// GrB_Vector_assign_[SCALAR]:  scalar expansion assignment to subvector
//------------------------------------------------------------------------------

// Assigns a single scalar to a subvector, w<mask>(I) = accum(w(I),x).  The
// scalar x is implicitly expanded into a vector u of size ni-by-1, with each
// entry in u equal to x, and then w<mask>(I) = accum(w(I),u) is done.

// Each of these can be used with their generic name, GrB_assign.

GrB_Info GrB_Vector_assign_BOOL     // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w(I),x)
    const bool x,                   // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_INT8     // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const int8_t x,                 // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_UINT8    // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const uint8_t x,                // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_INT16    // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const int16_t x,                // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_UINT16   // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const uint16_t x,               // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_INT32    // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const int32_t x,                // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_UINT32   // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const uint32_t x,               // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_INT64    // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const int64_t x,                // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_UINT64   // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const uint64_t x,               // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_FP32     // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const float x,                  // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_FP64     // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const double x,                 // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Vector_assign_UDT      // w<mask>(I) = accum (w(I),x)
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(w(I),x)
    const void *x,                  // scalar to assign to w(I)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

//------------------------------------------------------------------------------
// GrB_Matrix_assign_[SCALAR]:  scalar expansion assignment to submatrix
//------------------------------------------------------------------------------

// Assigns a single scalar to a submatrix, C<Mask>(I,J) = accum(C(I,J),x).  The
// scalar x is implicitly expanded into a matrix A of size ni-by-nj, with each
// entry in A equal to x, and then C<Mask>(I,J) = accum(C(I,J),A) is done.

// Each of these can be used with their generic name, GrB_assign.

GrB_Info GrB_Matrix_assign_BOOL     // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const bool x,                   // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_INT8     // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const int8_t x,                 // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_UINT8    // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const uint8_t x,                // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_INT16    // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const int16_t x,                // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_UINT16   // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const uint16_t x,               // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_INT32    // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const int32_t x,                // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_UINT32   // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const uint32_t x,               // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_INT64    // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const int64_t x,                // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_UINT64   // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const uint64_t x,               // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_FP32     // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const float x,                  // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_FP64     // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const double x,                 // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

GrB_Info GrB_Matrix_assign_UDT      // C<Mask>(I,J) = accum (C(I,J),x)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C(I,J),x)
    const void *x,                  // scalar to assign to C(I,J)
    const GrB_Index *I,             // row indices
    const GrB_Index ni,             // number of row indices
    const GrB_Index *J,             // column indices
    const GrB_Index nj,             // number of column indices
    const GrB_Descriptor desc       // descriptor for C and Mask
) ;

//------------------------------------------------------------------------------
// GrB_assign: generic submatrix/subvector assignment
//------------------------------------------------------------------------------

// GrB_assign is a generic function that provides access to all specific
// GrB_*_assign* functions:

// GrB_Vector_assign   (w,mask,acc,u,I,ni,d)     // w<mask>(I)   =acc(w(I),u)
// GrB_Matrix_assign   (C,Mask,acc,A,I,ni,J,nj,d)// C<Mask>(I,J) =acc(C(I,J),A)
// GrB_Col_assign      (C,mask,acc,u,I,ni,j,d)   // C<mask>(I,j) =acc(C(I,j),u)
// GrB_Row_assign      (C,mask,acc,u,i,J,nj,d)   // C<mask'>(i,J)=acc(C(i,J),u')
// GrB_Vector_assign_T (w,mask,acc,x,I,ni,d)     // w<mask>(I)   =acc(w(I),x)
// GrB_Matrix_assign_T (C,Mask,acc,x,I,ni,J,nj,d)// C<Mask>(I,J) =acc(C(I,J),x)

#define GrB_assign(arg1,Mask,accum,arg4,arg5,...)               \
    _Generic                                                    \
    (                                                           \
        (arg1),                                                 \
        GrB_Vector :                                            \
            _Generic                                            \
            (                                                   \
                (arg4),                                         \
                const bool       : GrB_Vector_assign_BOOL   ,   \
                      bool       : GrB_Vector_assign_BOOL   ,   \
                const int8_t     : GrB_Vector_assign_INT8   ,   \
                      int8_t     : GrB_Vector_assign_INT8   ,   \
                const uint8_t    : GrB_Vector_assign_UINT8  ,   \
                      uint8_t    : GrB_Vector_assign_UINT8  ,   \
                const int16_t    : GrB_Vector_assign_INT16  ,   \
                      int16_t    : GrB_Vector_assign_INT16  ,   \
                const uint16_t   : GrB_Vector_assign_UINT16 ,   \
                      uint16_t   : GrB_Vector_assign_UINT16 ,   \
                const int32_t    : GrB_Vector_assign_INT32  ,   \
                      int32_t    : GrB_Vector_assign_INT32  ,   \
                const uint32_t   : GrB_Vector_assign_UINT32 ,   \
                      uint32_t   : GrB_Vector_assign_UINT32 ,   \
                const int64_t    : GrB_Vector_assign_INT64  ,   \
                      int64_t    : GrB_Vector_assign_INT64  ,   \
                const uint64_t   : GrB_Vector_assign_UINT64 ,   \
                      uint64_t   : GrB_Vector_assign_UINT64 ,   \
                const float      : GrB_Vector_assign_FP32   ,   \
                      float      : GrB_Vector_assign_FP32   ,   \
                const double     : GrB_Vector_assign_FP64   ,   \
                      double     : GrB_Vector_assign_FP64   ,   \
                const void *     : GrB_Vector_assign_UDT    ,   \
                      void *     : GrB_Vector_assign_UDT    ,   \
                default          : GrB_Vector_assign            \
            ),                                                  \
        default :                                               \
            _Generic                                            \
            (                                                   \
                (arg4),                                         \
                const bool       : GrB_Matrix_assign_BOOL   ,   \
                      bool       : GrB_Matrix_assign_BOOL   ,   \
                const int8_t     : GrB_Matrix_assign_INT8   ,   \
                      int8_t     : GrB_Matrix_assign_INT8   ,   \
                const uint8_t    : GrB_Matrix_assign_UINT8  ,   \
                      uint8_t    : GrB_Matrix_assign_UINT8  ,   \
                const int16_t    : GrB_Matrix_assign_INT16  ,   \
                      int16_t    : GrB_Matrix_assign_INT16  ,   \
                const uint16_t   : GrB_Matrix_assign_UINT16 ,   \
                      uint16_t   : GrB_Matrix_assign_UINT16 ,   \
                const int32_t    : GrB_Matrix_assign_INT32  ,   \
                      int32_t    : GrB_Matrix_assign_INT32  ,   \
                const uint32_t   : GrB_Matrix_assign_UINT32 ,   \
                      uint32_t   : GrB_Matrix_assign_UINT32 ,   \
                const int64_t    : GrB_Matrix_assign_INT64  ,   \
                      int64_t    : GrB_Matrix_assign_INT64  ,   \
                const uint64_t   : GrB_Matrix_assign_UINT64 ,   \
                      uint64_t   : GrB_Matrix_assign_UINT64 ,   \
                const float      : GrB_Matrix_assign_FP32   ,   \
                      float      : GrB_Matrix_assign_FP32   ,   \
                const double     : GrB_Matrix_assign_FP64   ,   \
                      double     : GrB_Matrix_assign_FP64   ,   \
                const void *     : GrB_Matrix_assign_UDT    ,   \
                      void *     : GrB_Matrix_assign_UDT    ,   \
                const GrB_Vector :                              \
                    _Generic                                    \
                    (                                           \
                        (arg5),                                 \
                        const GrB_Index *: GrB_Col_assign ,     \
                              GrB_Index *: GrB_Col_assign ,     \
                        default          : GrB_Row_assign       \
                    ),                                          \
                GrB_Vector :                                    \
                    _Generic                                    \
                    (                                           \
                        (arg5),                                 \
                        const GrB_Index *: GrB_Col_assign ,     \
                              GrB_Index *: GrB_Col_assign ,     \
                        default          : GrB_Row_assign       \
                    ),                                          \
                default    : GrB_Matrix_assign                  \
            )                                                   \
    )                                                           \
    (arg1, Mask, accum, arg4, arg5, __VA_ARGS__)


//------------------------------------------------------------------------------
// matrix and vector apply
//------------------------------------------------------------------------------

// Apply a unary operator to the entries in a matrix or vector,
// C<Mask> = accum (C, op (A)).

// The input matrix A may be optionally transposed first, via the descriptor.

GrB_Info GrB_Vector_apply           // w<mask> = accum (w, op(u))
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_UnaryOp op,           // operator to apply to the entries
    const GrB_Vector u,             // first input:  vector u
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GrB_Matrix_apply           // C<Mask> = accum (C, op(A)) or op(A')
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_UnaryOp op,           // operator to apply to the entries
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Descriptor desc       // descriptor for C, mask, and A
) ;

//------------------------------------------------------------------------------
// GrB_apply: generic matrix/vector apply
//------------------------------------------------------------------------------

// GrB_apply is a generic function for applying a unary operator to a matrix
// or vector and provides access to these functions:

// GrB_Vector_apply (w,mask,acc,op,u,d)  // w<mask> = accum (w, op(u))
// GrB_Matrix_apply (C,Mask,acc,op,A,d)  // C<Mask> = accum (C, op(A))

#define GrB_apply(C,Mask,accum,op,A,desc)       \
    _Generic                                    \
    (                                           \
        (C),                                    \
        GrB_Vector   : GrB_Vector_apply ,       \
        GrB_Matrix   : GrB_Matrix_apply         \
    )                                           \
    (C, Mask, accum, op, A, desc)

//------------------------------------------------------------------------------
// matrix and vector selection
//------------------------------------------------------------------------------

// Select a subset of entries from a matrix or vector.
// C<Mask> = accum (C, op (A,k)), where the entries of op(A,k) are a subset of
// the entries of A.

// The input matrix A may be optionally transposed first, via the descriptor.

GrB_Info GxB_Vector_select          // w<mask> = accum (w, op(u,k))
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GxB_SelectOp op,          // operator to apply to the entries
    const GrB_Vector u,             // first input:  vector u
    const void *k,                  // optional input for the select operator
    const GrB_Descriptor desc       // descriptor for w and mask
) ;

GrB_Info GxB_Matrix_select          // C<Mask> = accum (C, op(A,k)) or op(A',k)
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GxB_SelectOp op,          // operator to apply to the entries
    const GrB_Matrix A,             // first input:  matrix A
    const void *k,                  // optional input for the select operator
    const GrB_Descriptor desc       // descriptor for C, mask, and A
) ;

//------------------------------------------------------------------------------
// GxB_select: generic matrix/vector select
//------------------------------------------------------------------------------

// GrB_select is a generic function for applying a select operator to a matrix
// or vector and provides access to these functions:

// GrB_Vector_select (w,mask,acc,op,u,k,d)  // w<mask> = accum (w, op(u,k))
// GrB_Matrix_select (C,Mask,acc,op,A,k,d)  // C<Mask> = accum (C, op(A,k))

#define GxB_select(C,Mask,accum,op,A,k,desc)    \
    _Generic                                    \
    (                                           \
        (C),                                    \
        GrB_Vector   : GxB_Vector_select ,      \
        GrB_Matrix   : GxB_Matrix_select        \
    )                                           \
    (C, Mask, accum, op, A, k, desc)


//------------------------------------------------------------------------------
// matrix and vector reduction
//------------------------------------------------------------------------------

// Reduce the entries in a matrix to a vector.  By default these methods
// compute a column vector t such that t(i) = sum (A (i,:)), and where "sum" is
// a commutative and associative monoid with an identity value.  A can be
// transposed, which reduces down the columns instead of the rows.  This
// behavior is the transpose of the MATLAB convention, where r=sum(A) produces
// a row vector and sums each column.

GrB_Info GrB_Matrix_reduce_Monoid   // w<mask> = accum (w,reduce(A))
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_Monoid monoid,        // reduce operator for t=reduce(A)
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Descriptor desc       // descriptor for w, mask, and A
) ;

GrB_Info GrB_Matrix_reduce_BinaryOp // w<mask> = accum (w,reduce(A))
(
    GrB_Vector w,                   // input/output vector for results
    const GrB_Vector mask,          // optional mask for w, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for z=accum(w,t)
    const GrB_BinaryOp op,          // reduce operator for t=reduce(A)
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Descriptor desc       // descriptor for w, mask, and A
) ;

//------------------------------------------------------------------------------
// reduce a vector to a scalar
//------------------------------------------------------------------------------

// Reduce entries in a vector to a scalar, c = accum (c, reduce_to_scalar(u))

// All entries in the vector are "summed" to a single scalar t using the reduce
// monoid, which must be associative (otherwise the results are undefined).
// The result is either assigned to the output scalar c (if accum is NULL), or
// it accumulated in the result c via c = accum(c,t).  If the vector has no
// entries, the result t is the identity value of the monoid.  Unlike most
// other GraphBLAS operations, this operation uses an accum operator but no
// mask.

// Like all GraphBLAS operations, these take a last argument of a GraphBLAS
// descriptor.  However, it is unused in the current GraphBLAS spec.  It may be
// used in the future.

GrB_Info GrB_Vector_reduce_BOOL     // c = accum (c, reduce_to_scalar (u))
(
    bool *c,                        // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_INT8     // c = accum (c, reduce_to_scalar (u))
(
    int8_t *c,                      // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_UINT8    // c = accum (c, reduce_to_scalar (u))
(
    uint8_t *c,                     // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_INT16    // c = accum (c, reduce_to_scalar (u))
(
    int16_t *c,                     // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_UINT16   // c = accum (c, reduce_to_scalar (u))
(
    uint16_t *c,                    // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_INT32    // c = accum (c, reduce_to_scalar (u))
(
    int32_t *c,                     // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_UINT32   // c = accum (c, reduce_to_scalar (u))
(
    uint32_t *c,                    // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_INT64    // c = accum (c, reduce_to_scalar (u))
(
    int64_t *c,                     // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_UINT64   // c = accum (c, reduce_to_scalar (u))
(
    uint64_t *c,                    // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_FP32     // c = accum (c, reduce_to_scalar (u))
(
    float *c,                       // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_FP64     // c = accum (c, reduce_to_scalar (u))
(
    double *c,                      // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Vector_reduce_UDT      // c = accum (c, reduce_to_scalar (u))
(
    void *c,                        // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Vector u,             // vector to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

//------------------------------------------------------------------------------
// reduce a matrix to a scalar
//------------------------------------------------------------------------------

// Reduce entries in a matrix to a scalar, c = accum (c, reduce_to_scalar(A))

// All entries in the matrix are "summed" to a single scalar t using the reduce
// monoid, which must be associative (otherwise the results are undefined).
// The result is either assigned to the output scalar c (if accum is NULL), or
// it accumulated in the result c via c = accum(c,t).  If the matrix has no
// entries, the result t is the identity value of the monoid.  Unlike most
// other GraphBLAS operations, this operation uses an accum operator but no
// mask.

GrB_Info GrB_Matrix_reduce_BOOL     // c = accum (c, reduce_to_scalar (A))
(
    bool *c,                        // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_INT8     // c = accum (c, reduce_to_scalar (A))
(
    int8_t *c,                      // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_UINT8    // c = accum (c, reduce_to_scalar (A))
(
    uint8_t *c,                     // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_INT16    // c = accum (c, reduce_to_scalar (A))
(
    int16_t *c,                     // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_UINT16   // c = accum (c, reduce_to_scalar (A))
(
    uint16_t *c,                    // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_INT32    // c = accum (c, reduce_to_scalar (A))
(
    int32_t *c,                     // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_UINT32   // c = accum (c, reduce_to_scalar (A))
(
    uint32_t *c,                    // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_INT64    // c = accum (c, reduce_to_scalar (A))
(
    int64_t *c,                     // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_UINT64   // c = accum (c, reduce_to_scalar (A))
(
    uint64_t *c,                    // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_FP32     // c = accum (c, reduce_to_scalar (A))
(
    float *c,                       // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_FP64     // c = accum (c, reduce_to_scalar (A))
(
    double *c,                      // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

GrB_Info GrB_Matrix_reduce_UDT      // c = accum (c, reduce_to_scalar (A))
(
    void *c,                        // result scalar
    const GrB_BinaryOp accum,       // optional accum for c=accum(c,t)
    const GrB_Monoid monoid,        // monoid to do the reduction
    const GrB_Matrix A,             // matrix to reduce
    const GrB_Descriptor desc       // descriptor (currently unused)
) ;

//------------------------------------------------------------------------------
// GrB_reduce: generic matrix/vector reduction to a vector or scalar
//------------------------------------------------------------------------------

// GrB_reduce is a generic function that provides access to all GrB_*reduce*
// functions:

// reduce matrix to vector:
// GrB_Matrix_reduce_Monoid   (w,mask,acc,mo,A,d) // w<mask> = acc (w,reduce(A))
// GrB_Matrix_reduce_BinaryOp (w,mask,acc,op,A,d) // w<mask> = acc (w,reduce(A))
// GrB_Vector_reduce_[SCALAR] (c, acc,monoid,u, d)
// GrB_Matrix_reduce_[SCALAR] (c, acc,monoid,A, d)

#define GrB_reduce(arg1,arg2,arg3,arg4,...)                 \
    _Generic                                                \
    (                                                       \
        (arg4),                                             \
        const GrB_Vector :                                  \
            _Generic                                        \
            (                                               \
                (arg1),                                     \
                bool     * : GrB_Vector_reduce_BOOL   ,     \
                int8_t   * : GrB_Vector_reduce_INT8   ,     \
                uint8_t  * : GrB_Vector_reduce_UINT8  ,     \
                int16_t  * : GrB_Vector_reduce_INT16  ,     \
                uint16_t * : GrB_Vector_reduce_UINT16 ,     \
                int32_t  * : GrB_Vector_reduce_INT32  ,     \
                uint32_t * : GrB_Vector_reduce_UINT32 ,     \
                int64_t  * : GrB_Vector_reduce_INT64  ,     \
                uint64_t * : GrB_Vector_reduce_UINT64 ,     \
                float    * : GrB_Vector_reduce_FP32   ,     \
                double   * : GrB_Vector_reduce_FP64   ,     \
                default    : GrB_Vector_reduce_UDT          \
            ),                                              \
        GrB_Vector :                                        \
            _Generic                                        \
            (                                               \
                (arg1),                                     \
                bool     * : GrB_Vector_reduce_BOOL   ,     \
                int8_t   * : GrB_Vector_reduce_INT8   ,     \
                uint8_t  * : GrB_Vector_reduce_UINT8  ,     \
                int16_t  * : GrB_Vector_reduce_INT16  ,     \
                uint16_t * : GrB_Vector_reduce_UINT16 ,     \
                int32_t  * : GrB_Vector_reduce_INT32  ,     \
                uint32_t * : GrB_Vector_reduce_UINT32 ,     \
                int64_t  * : GrB_Vector_reduce_INT64  ,     \
                uint64_t * : GrB_Vector_reduce_UINT64 ,     \
                float    * : GrB_Vector_reduce_FP32   ,     \
                double   * : GrB_Vector_reduce_FP64   ,     \
                default    : GrB_Vector_reduce_UDT          \
            ),                                              \
        const GrB_Matrix :                                  \
            _Generic                                        \
            (                                               \
                (arg1),                                     \
                bool     * : GrB_Matrix_reduce_BOOL   ,     \
                int8_t   * : GrB_Matrix_reduce_INT8   ,     \
                uint8_t  * : GrB_Matrix_reduce_UINT8  ,     \
                int16_t  * : GrB_Matrix_reduce_INT16  ,     \
                uint16_t * : GrB_Matrix_reduce_UINT16 ,     \
                int32_t  * : GrB_Matrix_reduce_INT32  ,     \
                uint32_t * : GrB_Matrix_reduce_UINT32 ,     \
                int64_t  * : GrB_Matrix_reduce_INT64  ,     \
                uint64_t * : GrB_Matrix_reduce_UINT64 ,     \
                float    * : GrB_Matrix_reduce_FP32   ,     \
                double   * : GrB_Matrix_reduce_FP64   ,     \
                default    : GrB_Matrix_reduce_UDT          \
            ),                                              \
        GrB_Matrix :                                        \
            _Generic                                        \
            (                                               \
                (arg1),                                     \
                bool     * : GrB_Matrix_reduce_BOOL   ,     \
                int8_t   * : GrB_Matrix_reduce_INT8   ,     \
                uint8_t  * : GrB_Matrix_reduce_UINT8  ,     \
                int16_t  * : GrB_Matrix_reduce_INT16  ,     \
                uint16_t * : GrB_Matrix_reduce_UINT16 ,     \
                int32_t  * : GrB_Matrix_reduce_INT32  ,     \
                uint32_t * : GrB_Matrix_reduce_UINT32 ,     \
                int64_t  * : GrB_Matrix_reduce_INT64  ,     \
                uint64_t * : GrB_Matrix_reduce_UINT64 ,     \
                float    * : GrB_Matrix_reduce_FP32   ,     \
                double   * : GrB_Matrix_reduce_FP64   ,     \
                default    : GrB_Matrix_reduce_UDT          \
            ),                                              \
        const GrB_Monoid   : GrB_Matrix_reduce_Monoid   ,   \
              GrB_Monoid   : GrB_Matrix_reduce_Monoid   ,   \
        const GrB_BinaryOp : GrB_Matrix_reduce_BinaryOp ,   \
              GrB_BinaryOp : GrB_Matrix_reduce_BinaryOp     \
    )                                                       \
    (arg1, arg2, arg3, arg4, __VA_ARGS__)

//------------------------------------------------------------------------------
// matrix transpose
//------------------------------------------------------------------------------

// T = A' is computed by default, but A can also be transposed via the
// descriptor.  In this case A is not transposed at all, and T = A.  The result
// is then passed through the Mask and accum, like almost all other GraphBLAS
// operations.  This makes GrB_transpose a direct interface to the accum/mask
// operation, C<Mask> = accum (C,A), or C<Mask> = accum (C,A') by default.

GrB_Info GrB_transpose              // C<Mask> = accum (C, A')
(
    GrB_Matrix C,                   // input/output matrix for results
    const GrB_Matrix Mask,          // optional mask for C, unused if NULL
    const GrB_BinaryOp accum,       // optional accum for Z=accum(C,T)
    const GrB_Matrix A,             // first input:  matrix A
    const GrB_Descriptor desc       // descriptor for C, Mask, and A
) ;



//==============================================================================
// additional predefined objects
//==============================================================================

// SPEC: predefined monoids and semirings are extensions to the spec

//------------------------------------------------------------------------------
// built-in monoids
//------------------------------------------------------------------------------

// 44 unique monoids can be constructed using built-in types and operators, all
// of which are defined below.  Four operators (min, max, plus, times) are
// available for each of the 10 non-Boolean types, and four purely Boolean
// monoids are available.

extern GrB_Monoid

    // MIN monoids:
    GxB_MIN_INT8_MONOID,          // identity: INT8_MAX
    GxB_MIN_UINT8_MONOID,         // identity: UINT8_MAX
    GxB_MIN_INT16_MONOID,         // identity: INT16_MAX
    GxB_MIN_UINT16_MONOID,        // identity: UINT16_MAX
    GxB_MIN_INT32_MONOID,         // identity: INT32_MAX
    GxB_MIN_UINT32_MONOID,        // identity: UINT32_MAX
    GxB_MIN_INT64_MONOID,         // identity: INT64_MAX
    GxB_MIN_UINT64_MONOID,        // identity: UINT64_MAX
    GxB_MIN_FP32_MONOID,          // identity: INFINITY
    GxB_MIN_FP64_MONOID,          // identity: INFINITY

    // MAX monoids:
    GxB_MAX_INT8_MONOID,          // identity: INT8_MIN
    GxB_MAX_UINT8_MONOID,         // identity: 0
    GxB_MAX_INT16_MONOID,         // identity: INT16_MIN
    GxB_MAX_UINT16_MONOID,        // identity: 0
    GxB_MAX_INT32_MONOID,         // identity: INT32_MIN
    GxB_MAX_UINT32_MONOID,        // identity: 0
    GxB_MAX_INT64_MONOID,         // identity: INT64_MIN
    GxB_MAX_UINT64_MONOID,        // identity: 0
    GxB_MAX_FP32_MONOID,          // identity: -INFINITY
    GxB_MAX_FP64_MONOID,          // identity: -INFINITY

    // PLUS monoids:
    GxB_PLUS_INT8_MONOID,         // identity: 0
    GxB_PLUS_UINT8_MONOID,        // identity: 0
    GxB_PLUS_INT16_MONOID,        // identity: 0
    GxB_PLUS_UINT16_MONOID,       // identity: 0
    GxB_PLUS_INT32_MONOID,        // identity: 0
    GxB_PLUS_UINT32_MONOID,       // identity: 0
    GxB_PLUS_INT64_MONOID,        // identity: 0
    GxB_PLUS_UINT64_MONOID,       // identity: 0
    GxB_PLUS_FP32_MONOID,         // identity: 0
    GxB_PLUS_FP64_MONOID,         // identity: 0

    // TIMES monoids:
    GxB_TIMES_INT8_MONOID,        // identity: 1
    GxB_TIMES_UINT8_MONOID,       // identity: 1
    GxB_TIMES_INT16_MONOID,       // identity: 1
    GxB_TIMES_UINT16_MONOID,      // identity: 1
    GxB_TIMES_INT32_MONOID,       // identity: 1
    GxB_TIMES_UINT32_MONOID,      // identity: 1
    GxB_TIMES_INT64_MONOID,       // identity: 1
    GxB_TIMES_UINT64_MONOID,      // identity: 1
    GxB_TIMES_FP32_MONOID,        // identity: 1
    GxB_TIMES_FP64_MONOID,        // identity: 1

    // Boolean monoids:
    GxB_LOR_BOOL_MONOID,          // identity: false
    GxB_LAND_BOOL_MONOID,         // identity: true
    GxB_LXOR_BOOL_MONOID,         // identity: false
    GxB_EQ_BOOL_MONOID ;          // identity: true

//------------------------------------------------------------------------------
// built-in semirings
//------------------------------------------------------------------------------

// Using built-in types and operators, 960 unique semirings can be built.  This
// count excludes redundant Boolean operators (for example GxB_TIMES_BOOL and
// GxB_LAND_BOOL are different operators but they are redundant since they
// always return the same result):

// 680 semirings with a multiply operator TxT -> T where T is non-Boolean, from
// the complete cross product of:

//      4 add monoids (MIN, MAX, PLUS, TIMES)
//      17 multiply operators:
//          (FIRST, SECOND, MIN, MAX, PLUS, MINUS, TIMES, DIV,
//           ISEQ, ISNE, ISGT, ISLT, ISGE, ISLE,
//           LOR, LAND, LXOR)
//      10 non-Boolean types, T

// 240 semirings with a comparison operator TxT -> bool, where T is
// non-Boolean, from the complete cross product of:

//      4 Boolean add monoids: (LAND, LOR, LXOR, EQ)
//      6 multiply operators: (EQ, NE, GT, LT, GE, LE)
//      10 non-Boolean types, T

// 40 semirings with purely Boolean types, bool x bool -> bool, from the
// complete cross product of:

//      4 Boolean add monoids (LAND, LOR, LXOR, EQ)
//      10 multiply operators:
//          (FIRST, SECOND, LOR, LAND, LXOR, EQ, GT, LT, GE, LE)

// In the names below, each semiring has a name of the form GxB_add_mult_T
// where add is the additive monoid, mult is the multiply operator, and T is
// the type.  The type T is always the type of x and y for the z=mult(x,y)
// operator.  The monoid's three types and the ztype of the mult operator are
// always the same.  This is the type T for the first set, and Boolean for
// the second and third sets of semirngs.

extern GrB_Semiring

//------------------------------------------------------------------------------
// 680 non-Boolean semirings where all types are the same, given by suffix _T
//------------------------------------------------------------------------------

// semirings with multiply op: z = FIRST (x,y), all types x,y,z the same:
GxB_MIN_FIRST_INT8     , GxB_MAX_FIRST_INT8     , GxB_PLUS_FIRST_INT8    , GxB_TIMES_FIRST_INT8   ,
GxB_MIN_FIRST_UINT8    , GxB_MAX_FIRST_UINT8    , GxB_PLUS_FIRST_UINT8   , GxB_TIMES_FIRST_UINT8  ,
GxB_MIN_FIRST_INT16    , GxB_MAX_FIRST_INT16    , GxB_PLUS_FIRST_INT16   , GxB_TIMES_FIRST_INT16  ,
GxB_MIN_FIRST_UINT16   , GxB_MAX_FIRST_UINT16   , GxB_PLUS_FIRST_UINT16  , GxB_TIMES_FIRST_UINT16 ,
GxB_MIN_FIRST_INT32    , GxB_MAX_FIRST_INT32    , GxB_PLUS_FIRST_INT32   , GxB_TIMES_FIRST_INT32  ,
GxB_MIN_FIRST_UINT32   , GxB_MAX_FIRST_UINT32   , GxB_PLUS_FIRST_UINT32  , GxB_TIMES_FIRST_UINT32 ,
GxB_MIN_FIRST_INT64    , GxB_MAX_FIRST_INT64    , GxB_PLUS_FIRST_INT64   , GxB_TIMES_FIRST_INT64  ,
GxB_MIN_FIRST_UINT64   , GxB_MAX_FIRST_UINT64   , GxB_PLUS_FIRST_UINT64  , GxB_TIMES_FIRST_UINT64 ,
GxB_MIN_FIRST_FP32     , GxB_MAX_FIRST_FP32     , GxB_PLUS_FIRST_FP32    , GxB_TIMES_FIRST_FP32   ,
GxB_MIN_FIRST_FP64     , GxB_MAX_FIRST_FP64     , GxB_PLUS_FIRST_FP64    , GxB_TIMES_FIRST_FP64   ,

// semirings with multiply op: z = SECOND (x,y), all types x,y,z the same:
GxB_MIN_SECOND_INT8    , GxB_MAX_SECOND_INT8    , GxB_PLUS_SECOND_INT8   , GxB_TIMES_SECOND_INT8  ,
GxB_MIN_SECOND_UINT8   , GxB_MAX_SECOND_UINT8   , GxB_PLUS_SECOND_UINT8  , GxB_TIMES_SECOND_UINT8 ,
GxB_MIN_SECOND_INT16   , GxB_MAX_SECOND_INT16   , GxB_PLUS_SECOND_INT16  , GxB_TIMES_SECOND_INT16 ,
GxB_MIN_SECOND_UINT16  , GxB_MAX_SECOND_UINT16  , GxB_PLUS_SECOND_UINT16 , GxB_TIMES_SECOND_UINT16,
GxB_MIN_SECOND_INT32   , GxB_MAX_SECOND_INT32   , GxB_PLUS_SECOND_INT32  , GxB_TIMES_SECOND_INT32 ,
GxB_MIN_SECOND_UINT32  , GxB_MAX_SECOND_UINT32  , GxB_PLUS_SECOND_UINT32 , GxB_TIMES_SECOND_UINT32,
GxB_MIN_SECOND_INT64   , GxB_MAX_SECOND_INT64   , GxB_PLUS_SECOND_INT64  , GxB_TIMES_SECOND_INT64 ,
GxB_MIN_SECOND_UINT64  , GxB_MAX_SECOND_UINT64  , GxB_PLUS_SECOND_UINT64 , GxB_TIMES_SECOND_UINT64,
GxB_MIN_SECOND_FP32    , GxB_MAX_SECOND_FP32    , GxB_PLUS_SECOND_FP32   , GxB_TIMES_SECOND_FP32  ,
GxB_MIN_SECOND_FP64    , GxB_MAX_SECOND_FP64    , GxB_PLUS_SECOND_FP64   , GxB_TIMES_SECOND_FP64  ,

// semirings with multiply op: z = MIN (x,y), all types x,y,z the same:
GxB_MIN_MIN_INT8       , GxB_MAX_MIN_INT8       , GxB_PLUS_MIN_INT8      , GxB_TIMES_MIN_INT8     ,
GxB_MIN_MIN_UINT8      , GxB_MAX_MIN_UINT8      , GxB_PLUS_MIN_UINT8     , GxB_TIMES_MIN_UINT8    ,
GxB_MIN_MIN_INT16      , GxB_MAX_MIN_INT16      , GxB_PLUS_MIN_INT16     , GxB_TIMES_MIN_INT16    ,
GxB_MIN_MIN_UINT16     , GxB_MAX_MIN_UINT16     , GxB_PLUS_MIN_UINT16    , GxB_TIMES_MIN_UINT16   ,
GxB_MIN_MIN_INT32      , GxB_MAX_MIN_INT32      , GxB_PLUS_MIN_INT32     , GxB_TIMES_MIN_INT32    ,
GxB_MIN_MIN_UINT32     , GxB_MAX_MIN_UINT32     , GxB_PLUS_MIN_UINT32    , GxB_TIMES_MIN_UINT32   ,
GxB_MIN_MIN_INT64      , GxB_MAX_MIN_INT64      , GxB_PLUS_MIN_INT64     , GxB_TIMES_MIN_INT64    ,
GxB_MIN_MIN_UINT64     , GxB_MAX_MIN_UINT64     , GxB_PLUS_MIN_UINT64    , GxB_TIMES_MIN_UINT64   ,
GxB_MIN_MIN_FP32       , GxB_MAX_MIN_FP32       , GxB_PLUS_MIN_FP32      , GxB_TIMES_MIN_FP32     ,
GxB_MIN_MIN_FP64       , GxB_MAX_MIN_FP64       , GxB_PLUS_MIN_FP64      , GxB_TIMES_MIN_FP64     ,

// semirings with multiply op: z = MAX (x,y), all types x,y,z the same:
GxB_MIN_MAX_INT8       , GxB_MAX_MAX_INT8       , GxB_PLUS_MAX_INT8      , GxB_TIMES_MAX_INT8     ,
GxB_MIN_MAX_UINT8      , GxB_MAX_MAX_UINT8      , GxB_PLUS_MAX_UINT8     , GxB_TIMES_MAX_UINT8    ,
GxB_MIN_MAX_INT16      , GxB_MAX_MAX_INT16      , GxB_PLUS_MAX_INT16     , GxB_TIMES_MAX_INT16    ,
GxB_MIN_MAX_UINT16     , GxB_MAX_MAX_UINT16     , GxB_PLUS_MAX_UINT16    , GxB_TIMES_MAX_UINT16   ,
GxB_MIN_MAX_INT32      , GxB_MAX_MAX_INT32      , GxB_PLUS_MAX_INT32     , GxB_TIMES_MAX_INT32    ,
GxB_MIN_MAX_UINT32     , GxB_MAX_MAX_UINT32     , GxB_PLUS_MAX_UINT32    , GxB_TIMES_MAX_UINT32   ,
GxB_MIN_MAX_INT64      , GxB_MAX_MAX_INT64      , GxB_PLUS_MAX_INT64     , GxB_TIMES_MAX_INT64    ,
GxB_MIN_MAX_UINT64     , GxB_MAX_MAX_UINT64     , GxB_PLUS_MAX_UINT64    , GxB_TIMES_MAX_UINT64   ,
GxB_MIN_MAX_FP32       , GxB_MAX_MAX_FP32       , GxB_PLUS_MAX_FP32      , GxB_TIMES_MAX_FP32     ,
GxB_MIN_MAX_FP64       , GxB_MAX_MAX_FP64       , GxB_PLUS_MAX_FP64      , GxB_TIMES_MAX_FP64     ,

// semirings with multiply op: z = PLUS (x,y), all types x,y,z the same:
GxB_MIN_PLUS_INT8      , GxB_MAX_PLUS_INT8      , GxB_PLUS_PLUS_INT8     , GxB_TIMES_PLUS_INT8    ,
GxB_MIN_PLUS_UINT8     , GxB_MAX_PLUS_UINT8     , GxB_PLUS_PLUS_UINT8    , GxB_TIMES_PLUS_UINT8   ,
GxB_MIN_PLUS_INT16     , GxB_MAX_PLUS_INT16     , GxB_PLUS_PLUS_INT16    , GxB_TIMES_PLUS_INT16   ,
GxB_MIN_PLUS_UINT16    , GxB_MAX_PLUS_UINT16    , GxB_PLUS_PLUS_UINT16   , GxB_TIMES_PLUS_UINT16  ,
GxB_MIN_PLUS_INT32     , GxB_MAX_PLUS_INT32     , GxB_PLUS_PLUS_INT32    , GxB_TIMES_PLUS_INT32   ,
GxB_MIN_PLUS_UINT32    , GxB_MAX_PLUS_UINT32    , GxB_PLUS_PLUS_UINT32   , GxB_TIMES_PLUS_UINT32  ,
GxB_MIN_PLUS_INT64     , GxB_MAX_PLUS_INT64     , GxB_PLUS_PLUS_INT64    , GxB_TIMES_PLUS_INT64   ,
GxB_MIN_PLUS_UINT64    , GxB_MAX_PLUS_UINT64    , GxB_PLUS_PLUS_UINT64   , GxB_TIMES_PLUS_UINT64  ,
GxB_MIN_PLUS_FP32      , GxB_MAX_PLUS_FP32      , GxB_PLUS_PLUS_FP32     , GxB_TIMES_PLUS_FP32    ,
GxB_MIN_PLUS_FP64      , GxB_MAX_PLUS_FP64      , GxB_PLUS_PLUS_FP64     , GxB_TIMES_PLUS_FP64    ,

// semirings with multiply op: z = MINUS (x,y), all types x,y,z the same:
GxB_MIN_MINUS_INT8     , GxB_MAX_MINUS_INT8     , GxB_PLUS_MINUS_INT8    , GxB_TIMES_MINUS_INT8   ,
GxB_MIN_MINUS_UINT8    , GxB_MAX_MINUS_UINT8    , GxB_PLUS_MINUS_UINT8   , GxB_TIMES_MINUS_UINT8  ,
GxB_MIN_MINUS_INT16    , GxB_MAX_MINUS_INT16    , GxB_PLUS_MINUS_INT16   , GxB_TIMES_MINUS_INT16  ,
GxB_MIN_MINUS_UINT16   , GxB_MAX_MINUS_UINT16   , GxB_PLUS_MINUS_UINT16  , GxB_TIMES_MINUS_UINT16 ,
GxB_MIN_MINUS_INT32    , GxB_MAX_MINUS_INT32    , GxB_PLUS_MINUS_INT32   , GxB_TIMES_MINUS_INT32  ,
GxB_MIN_MINUS_UINT32   , GxB_MAX_MINUS_UINT32   , GxB_PLUS_MINUS_UINT32  , GxB_TIMES_MINUS_UINT32 ,
GxB_MIN_MINUS_INT64    , GxB_MAX_MINUS_INT64    , GxB_PLUS_MINUS_INT64   , GxB_TIMES_MINUS_INT64  ,
GxB_MIN_MINUS_UINT64   , GxB_MAX_MINUS_UINT64   , GxB_PLUS_MINUS_UINT64  , GxB_TIMES_MINUS_UINT64 ,
GxB_MIN_MINUS_FP32     , GxB_MAX_MINUS_FP32     , GxB_PLUS_MINUS_FP32    , GxB_TIMES_MINUS_FP32   ,
GxB_MIN_MINUS_FP64     , GxB_MAX_MINUS_FP64     , GxB_PLUS_MINUS_FP64    , GxB_TIMES_MINUS_FP64   ,

// semirings with multiply op: z = TIMES (x,y), all types x,y,z the same:
GxB_MIN_TIMES_INT8     , GxB_MAX_TIMES_INT8     , GxB_PLUS_TIMES_INT8    , GxB_TIMES_TIMES_INT8   ,
GxB_MIN_TIMES_UINT8    , GxB_MAX_TIMES_UINT8    , GxB_PLUS_TIMES_UINT8   , GxB_TIMES_TIMES_UINT8  ,
GxB_MIN_TIMES_INT16    , GxB_MAX_TIMES_INT16    , GxB_PLUS_TIMES_INT16   , GxB_TIMES_TIMES_INT16  ,
GxB_MIN_TIMES_UINT16   , GxB_MAX_TIMES_UINT16   , GxB_PLUS_TIMES_UINT16  , GxB_TIMES_TIMES_UINT16 ,
GxB_MIN_TIMES_INT32    , GxB_MAX_TIMES_INT32    , GxB_PLUS_TIMES_INT32   , GxB_TIMES_TIMES_INT32  ,
GxB_MIN_TIMES_UINT32   , GxB_MAX_TIMES_UINT32   , GxB_PLUS_TIMES_UINT32  , GxB_TIMES_TIMES_UINT32 ,
GxB_MIN_TIMES_INT64    , GxB_MAX_TIMES_INT64    , GxB_PLUS_TIMES_INT64   , GxB_TIMES_TIMES_INT64  ,
GxB_MIN_TIMES_UINT64   , GxB_MAX_TIMES_UINT64   , GxB_PLUS_TIMES_UINT64  , GxB_TIMES_TIMES_UINT64 ,
GxB_MIN_TIMES_FP32     , GxB_MAX_TIMES_FP32     , GxB_PLUS_TIMES_FP32    , GxB_TIMES_TIMES_FP32   ,
GxB_MIN_TIMES_FP64     , GxB_MAX_TIMES_FP64     , GxB_PLUS_TIMES_FP64    , GxB_TIMES_TIMES_FP64   ,

// semirings with multiply op: z = DIV (x,y), all types x,y,z the same:
GxB_MIN_DIV_INT8       , GxB_MAX_DIV_INT8       , GxB_PLUS_DIV_INT8      , GxB_TIMES_DIV_INT8     ,
GxB_MIN_DIV_UINT8      , GxB_MAX_DIV_UINT8      , GxB_PLUS_DIV_UINT8     , GxB_TIMES_DIV_UINT8    ,
GxB_MIN_DIV_INT16      , GxB_MAX_DIV_INT16      , GxB_PLUS_DIV_INT16     , GxB_TIMES_DIV_INT16    ,
GxB_MIN_DIV_UINT16     , GxB_MAX_DIV_UINT16     , GxB_PLUS_DIV_UINT16    , GxB_TIMES_DIV_UINT16   ,
GxB_MIN_DIV_INT32      , GxB_MAX_DIV_INT32      , GxB_PLUS_DIV_INT32     , GxB_TIMES_DIV_INT32    ,
GxB_MIN_DIV_UINT32     , GxB_MAX_DIV_UINT32     , GxB_PLUS_DIV_UINT32    , GxB_TIMES_DIV_UINT32   ,
GxB_MIN_DIV_INT64      , GxB_MAX_DIV_INT64      , GxB_PLUS_DIV_INT64     , GxB_TIMES_DIV_INT64    ,
GxB_MIN_DIV_UINT64     , GxB_MAX_DIV_UINT64     , GxB_PLUS_DIV_UINT64    , GxB_TIMES_DIV_UINT64   ,
GxB_MIN_DIV_FP32       , GxB_MAX_DIV_FP32       , GxB_PLUS_DIV_FP32      , GxB_TIMES_DIV_FP32     ,
GxB_MIN_DIV_FP64       , GxB_MAX_DIV_FP64       , GxB_PLUS_DIV_FP64      , GxB_TIMES_DIV_FP64     ,

// semirings with multiply op: z = ISEQ (x,y), all types x,y,z the same:
GxB_MIN_ISEQ_INT8      , GxB_MAX_ISEQ_INT8      , GxB_PLUS_ISEQ_INT8     , GxB_TIMES_ISEQ_INT8    ,
GxB_MIN_ISEQ_UINT8     , GxB_MAX_ISEQ_UINT8     , GxB_PLUS_ISEQ_UINT8    , GxB_TIMES_ISEQ_UINT8   ,
GxB_MIN_ISEQ_INT16     , GxB_MAX_ISEQ_INT16     , GxB_PLUS_ISEQ_INT16    , GxB_TIMES_ISEQ_INT16   ,
GxB_MIN_ISEQ_UINT16    , GxB_MAX_ISEQ_UINT16    , GxB_PLUS_ISEQ_UINT16   , GxB_TIMES_ISEQ_UINT16  ,
GxB_MIN_ISEQ_INT32     , GxB_MAX_ISEQ_INT32     , GxB_PLUS_ISEQ_INT32    , GxB_TIMES_ISEQ_INT32   ,
GxB_MIN_ISEQ_UINT32    , GxB_MAX_ISEQ_UINT32    , GxB_PLUS_ISEQ_UINT32   , GxB_TIMES_ISEQ_UINT32  ,
GxB_MIN_ISEQ_INT64     , GxB_MAX_ISEQ_INT64     , GxB_PLUS_ISEQ_INT64    , GxB_TIMES_ISEQ_INT64   ,
GxB_MIN_ISEQ_UINT64    , GxB_MAX_ISEQ_UINT64    , GxB_PLUS_ISEQ_UINT64   , GxB_TIMES_ISEQ_UINT64  ,
GxB_MIN_ISEQ_FP32      , GxB_MAX_ISEQ_FP32      , GxB_PLUS_ISEQ_FP32     , GxB_TIMES_ISEQ_FP32    ,
GxB_MIN_ISEQ_FP64      , GxB_MAX_ISEQ_FP64      , GxB_PLUS_ISEQ_FP64     , GxB_TIMES_ISEQ_FP64    ,

// semirings with multiply op: z = ISNE (x,y), all types x,y,z the same:
GxB_MIN_ISNE_INT8      , GxB_MAX_ISNE_INT8      , GxB_PLUS_ISNE_INT8     , GxB_TIMES_ISNE_INT8    ,
GxB_MIN_ISNE_UINT8     , GxB_MAX_ISNE_UINT8     , GxB_PLUS_ISNE_UINT8    , GxB_TIMES_ISNE_UINT8   ,
GxB_MIN_ISNE_INT16     , GxB_MAX_ISNE_INT16     , GxB_PLUS_ISNE_INT16    , GxB_TIMES_ISNE_INT16   ,
GxB_MIN_ISNE_UINT16    , GxB_MAX_ISNE_UINT16    , GxB_PLUS_ISNE_UINT16   , GxB_TIMES_ISNE_UINT16  ,
GxB_MIN_ISNE_INT32     , GxB_MAX_ISNE_INT32     , GxB_PLUS_ISNE_INT32    , GxB_TIMES_ISNE_INT32   ,
GxB_MIN_ISNE_UINT32    , GxB_MAX_ISNE_UINT32    , GxB_PLUS_ISNE_UINT32   , GxB_TIMES_ISNE_UINT32  ,
GxB_MIN_ISNE_INT64     , GxB_MAX_ISNE_INT64     , GxB_PLUS_ISNE_INT64    , GxB_TIMES_ISNE_INT64   ,
GxB_MIN_ISNE_UINT64    , GxB_MAX_ISNE_UINT64    , GxB_PLUS_ISNE_UINT64   , GxB_TIMES_ISNE_UINT64  ,
GxB_MIN_ISNE_FP32      , GxB_MAX_ISNE_FP32      , GxB_PLUS_ISNE_FP32     , GxB_TIMES_ISNE_FP32    ,
GxB_MIN_ISNE_FP64      , GxB_MAX_ISNE_FP64      , GxB_PLUS_ISNE_FP64     , GxB_TIMES_ISNE_FP64    ,

// semirings with multiply op: z = ISGT (x,y), all types x,y,z the same:
GxB_MIN_ISGT_INT8      , GxB_MAX_ISGT_INT8      , GxB_PLUS_ISGT_INT8     , GxB_TIMES_ISGT_INT8    ,
GxB_MIN_ISGT_UINT8     , GxB_MAX_ISGT_UINT8     , GxB_PLUS_ISGT_UINT8    , GxB_TIMES_ISGT_UINT8   ,
GxB_MIN_ISGT_INT16     , GxB_MAX_ISGT_INT16     , GxB_PLUS_ISGT_INT16    , GxB_TIMES_ISGT_INT16   ,
GxB_MIN_ISGT_UINT16    , GxB_MAX_ISGT_UINT16    , GxB_PLUS_ISGT_UINT16   , GxB_TIMES_ISGT_UINT16  ,
GxB_MIN_ISGT_INT32     , GxB_MAX_ISGT_INT32     , GxB_PLUS_ISGT_INT32    , GxB_TIMES_ISGT_INT32   ,
GxB_MIN_ISGT_UINT32    , GxB_MAX_ISGT_UINT32    , GxB_PLUS_ISGT_UINT32   , GxB_TIMES_ISGT_UINT32  ,
GxB_MIN_ISGT_INT64     , GxB_MAX_ISGT_INT64     , GxB_PLUS_ISGT_INT64    , GxB_TIMES_ISGT_INT64   ,
GxB_MIN_ISGT_UINT64    , GxB_MAX_ISGT_UINT64    , GxB_PLUS_ISGT_UINT64   , GxB_TIMES_ISGT_UINT64  ,
GxB_MIN_ISGT_FP32      , GxB_MAX_ISGT_FP32      , GxB_PLUS_ISGT_FP32     , GxB_TIMES_ISGT_FP32    ,
GxB_MIN_ISGT_FP64      , GxB_MAX_ISGT_FP64      , GxB_PLUS_ISGT_FP64     , GxB_TIMES_ISGT_FP64    ,

// semirings with multiply op: z = ISLT (x,y), all types x,y,z the same:
GxB_MIN_ISLT_INT8      , GxB_MAX_ISLT_INT8      , GxB_PLUS_ISLT_INT8     , GxB_TIMES_ISLT_INT8    ,
GxB_MIN_ISLT_UINT8     , GxB_MAX_ISLT_UINT8     , GxB_PLUS_ISLT_UINT8    , GxB_TIMES_ISLT_UINT8   ,
GxB_MIN_ISLT_INT16     , GxB_MAX_ISLT_INT16     , GxB_PLUS_ISLT_INT16    , GxB_TIMES_ISLT_INT16   ,
GxB_MIN_ISLT_UINT16    , GxB_MAX_ISLT_UINT16    , GxB_PLUS_ISLT_UINT16   , GxB_TIMES_ISLT_UINT16  ,
GxB_MIN_ISLT_INT32     , GxB_MAX_ISLT_INT32     , GxB_PLUS_ISLT_INT32    , GxB_TIMES_ISLT_INT32   ,
GxB_MIN_ISLT_UINT32    , GxB_MAX_ISLT_UINT32    , GxB_PLUS_ISLT_UINT32   , GxB_TIMES_ISLT_UINT32  ,
GxB_MIN_ISLT_INT64     , GxB_MAX_ISLT_INT64     , GxB_PLUS_ISLT_INT64    , GxB_TIMES_ISLT_INT64   ,
GxB_MIN_ISLT_UINT64    , GxB_MAX_ISLT_UINT64    , GxB_PLUS_ISLT_UINT64   , GxB_TIMES_ISLT_UINT64  ,
GxB_MIN_ISLT_FP32      , GxB_MAX_ISLT_FP32      , GxB_PLUS_ISLT_FP32     , GxB_TIMES_ISLT_FP32    ,
GxB_MIN_ISLT_FP64      , GxB_MAX_ISLT_FP64      , GxB_PLUS_ISLT_FP64     , GxB_TIMES_ISLT_FP64    ,

// semirings with multiply op: z = ISGE (x,y), all types x,y,z the same:
GxB_MIN_ISGE_INT8      , GxB_MAX_ISGE_INT8      , GxB_PLUS_ISGE_INT8     , GxB_TIMES_ISGE_INT8    ,
GxB_MIN_ISGE_UINT8     , GxB_MAX_ISGE_UINT8     , GxB_PLUS_ISGE_UINT8    , GxB_TIMES_ISGE_UINT8   ,
GxB_MIN_ISGE_INT16     , GxB_MAX_ISGE_INT16     , GxB_PLUS_ISGE_INT16    , GxB_TIMES_ISGE_INT16   ,
GxB_MIN_ISGE_UINT16    , GxB_MAX_ISGE_UINT16    , GxB_PLUS_ISGE_UINT16   , GxB_TIMES_ISGE_UINT16  ,
GxB_MIN_ISGE_INT32     , GxB_MAX_ISGE_INT32     , GxB_PLUS_ISGE_INT32    , GxB_TIMES_ISGE_INT32   ,
GxB_MIN_ISGE_UINT32    , GxB_MAX_ISGE_UINT32    , GxB_PLUS_ISGE_UINT32   , GxB_TIMES_ISGE_UINT32  ,
GxB_MIN_ISGE_INT64     , GxB_MAX_ISGE_INT64     , GxB_PLUS_ISGE_INT64    , GxB_TIMES_ISGE_INT64   ,
GxB_MIN_ISGE_UINT64    , GxB_MAX_ISGE_UINT64    , GxB_PLUS_ISGE_UINT64   , GxB_TIMES_ISGE_UINT64  ,
GxB_MIN_ISGE_FP32      , GxB_MAX_ISGE_FP32      , GxB_PLUS_ISGE_FP32     , GxB_TIMES_ISGE_FP32    ,
GxB_MIN_ISGE_FP64      , GxB_MAX_ISGE_FP64      , GxB_PLUS_ISGE_FP64     , GxB_TIMES_ISGE_FP64    ,

// semirings with multiply op: z = ISLE (x,y), all types x,y,z the same:
GxB_MIN_ISLE_INT8      , GxB_MAX_ISLE_INT8      , GxB_PLUS_ISLE_INT8     , GxB_TIMES_ISLE_INT8    ,
GxB_MIN_ISLE_UINT8     , GxB_MAX_ISLE_UINT8     , GxB_PLUS_ISLE_UINT8    , GxB_TIMES_ISLE_UINT8   ,
GxB_MIN_ISLE_INT16     , GxB_MAX_ISLE_INT16     , GxB_PLUS_ISLE_INT16    , GxB_TIMES_ISLE_INT16   ,
GxB_MIN_ISLE_UINT16    , GxB_MAX_ISLE_UINT16    , GxB_PLUS_ISLE_UINT16   , GxB_TIMES_ISLE_UINT16  ,
GxB_MIN_ISLE_INT32     , GxB_MAX_ISLE_INT32     , GxB_PLUS_ISLE_INT32    , GxB_TIMES_ISLE_INT32   ,
GxB_MIN_ISLE_UINT32    , GxB_MAX_ISLE_UINT32    , GxB_PLUS_ISLE_UINT32   , GxB_TIMES_ISLE_UINT32  ,
GxB_MIN_ISLE_INT64     , GxB_MAX_ISLE_INT64     , GxB_PLUS_ISLE_INT64    , GxB_TIMES_ISLE_INT64   ,
GxB_MIN_ISLE_UINT64    , GxB_MAX_ISLE_UINT64    , GxB_PLUS_ISLE_UINT64   , GxB_TIMES_ISLE_UINT64  ,
GxB_MIN_ISLE_FP32      , GxB_MAX_ISLE_FP32      , GxB_PLUS_ISLE_FP32     , GxB_TIMES_ISLE_FP32    ,
GxB_MIN_ISLE_FP64      , GxB_MAX_ISLE_FP64      , GxB_PLUS_ISLE_FP64     , GxB_TIMES_ISLE_FP64    ,

// semirings with multiply op: z = LOR (x,y), all types x,y,z the same:
GxB_MIN_LOR_INT8       , GxB_MAX_LOR_INT8       , GxB_PLUS_LOR_INT8      , GxB_TIMES_LOR_INT8     ,
GxB_MIN_LOR_UINT8      , GxB_MAX_LOR_UINT8      , GxB_PLUS_LOR_UINT8     , GxB_TIMES_LOR_UINT8    ,
GxB_MIN_LOR_INT16      , GxB_MAX_LOR_INT16      , GxB_PLUS_LOR_INT16     , GxB_TIMES_LOR_INT16    ,
GxB_MIN_LOR_UINT16     , GxB_MAX_LOR_UINT16     , GxB_PLUS_LOR_UINT16    , GxB_TIMES_LOR_UINT16   ,
GxB_MIN_LOR_INT32      , GxB_MAX_LOR_INT32      , GxB_PLUS_LOR_INT32     , GxB_TIMES_LOR_INT32    ,
GxB_MIN_LOR_UINT32     , GxB_MAX_LOR_UINT32     , GxB_PLUS_LOR_UINT32    , GxB_TIMES_LOR_UINT32   ,
GxB_MIN_LOR_INT64      , GxB_MAX_LOR_INT64      , GxB_PLUS_LOR_INT64     , GxB_TIMES_LOR_INT64    ,
GxB_MIN_LOR_UINT64     , GxB_MAX_LOR_UINT64     , GxB_PLUS_LOR_UINT64    , GxB_TIMES_LOR_UINT64   ,
GxB_MIN_LOR_FP32       , GxB_MAX_LOR_FP32       , GxB_PLUS_LOR_FP32      , GxB_TIMES_LOR_FP32     ,
GxB_MIN_LOR_FP64       , GxB_MAX_LOR_FP64       , GxB_PLUS_LOR_FP64      , GxB_TIMES_LOR_FP64     ,

// semirings with multiply op: z = LAND (x,y), all types x,y,z the same:
GxB_MIN_LAND_INT8      , GxB_MAX_LAND_INT8      , GxB_PLUS_LAND_INT8     , GxB_TIMES_LAND_INT8    ,
GxB_MIN_LAND_UINT8     , GxB_MAX_LAND_UINT8     , GxB_PLUS_LAND_UINT8    , GxB_TIMES_LAND_UINT8   ,
GxB_MIN_LAND_INT16     , GxB_MAX_LAND_INT16     , GxB_PLUS_LAND_INT16    , GxB_TIMES_LAND_INT16   ,
GxB_MIN_LAND_UINT16    , GxB_MAX_LAND_UINT16    , GxB_PLUS_LAND_UINT16   , GxB_TIMES_LAND_UINT16  ,
GxB_MIN_LAND_INT32     , GxB_MAX_LAND_INT32     , GxB_PLUS_LAND_INT32    , GxB_TIMES_LAND_INT32   ,
GxB_MIN_LAND_UINT32    , GxB_MAX_LAND_UINT32    , GxB_PLUS_LAND_UINT32   , GxB_TIMES_LAND_UINT32  ,
GxB_MIN_LAND_INT64     , GxB_MAX_LAND_INT64     , GxB_PLUS_LAND_INT64    , GxB_TIMES_LAND_INT64   ,
GxB_MIN_LAND_UINT64    , GxB_MAX_LAND_UINT64    , GxB_PLUS_LAND_UINT64   , GxB_TIMES_LAND_UINT64  ,
GxB_MIN_LAND_FP32      , GxB_MAX_LAND_FP32      , GxB_PLUS_LAND_FP32     , GxB_TIMES_LAND_FP32    ,
GxB_MIN_LAND_FP64      , GxB_MAX_LAND_FP64      , GxB_PLUS_LAND_FP64     , GxB_TIMES_LAND_FP64    ,

// semirings with multiply op: z = LXOR (x,y), all types x,y,z the same:
GxB_MIN_LXOR_INT8      , GxB_MAX_LXOR_INT8      , GxB_PLUS_LXOR_INT8     , GxB_TIMES_LXOR_INT8    ,
GxB_MIN_LXOR_UINT8     , GxB_MAX_LXOR_UINT8     , GxB_PLUS_LXOR_UINT8    , GxB_TIMES_LXOR_UINT8   ,
GxB_MIN_LXOR_INT16     , GxB_MAX_LXOR_INT16     , GxB_PLUS_LXOR_INT16    , GxB_TIMES_LXOR_INT16   ,
GxB_MIN_LXOR_UINT16    , GxB_MAX_LXOR_UINT16    , GxB_PLUS_LXOR_UINT16   , GxB_TIMES_LXOR_UINT16  ,
GxB_MIN_LXOR_INT32     , GxB_MAX_LXOR_INT32     , GxB_PLUS_LXOR_INT32    , GxB_TIMES_LXOR_INT32   ,
GxB_MIN_LXOR_UINT32    , GxB_MAX_LXOR_UINT32    , GxB_PLUS_LXOR_UINT32   , GxB_TIMES_LXOR_UINT32  ,
GxB_MIN_LXOR_INT64     , GxB_MAX_LXOR_INT64     , GxB_PLUS_LXOR_INT64    , GxB_TIMES_LXOR_INT64   ,
GxB_MIN_LXOR_UINT64    , GxB_MAX_LXOR_UINT64    , GxB_PLUS_LXOR_UINT64   , GxB_TIMES_LXOR_UINT64  ,
GxB_MIN_LXOR_FP32      , GxB_MAX_LXOR_FP32      , GxB_PLUS_LXOR_FP32     , GxB_TIMES_LXOR_FP32    ,
GxB_MIN_LXOR_FP64      , GxB_MAX_LXOR_FP64      , GxB_PLUS_LXOR_FP64     , GxB_TIMES_LXOR_FP64    ,

//------------------------------------------------------------------------------
// 240 semirings with comparison ops of the form TxT->bool, and Boolean monoids
//------------------------------------------------------------------------------

// semirings with multiply op: z = EQ (x,y), where z is Boolean and x,y are given by the suffix:
GxB_LOR_EQ_INT8        , GxB_LAND_EQ_INT8       , GxB_LXOR_EQ_INT8       , GxB_EQ_EQ_INT8         ,
GxB_LOR_EQ_UINT8       , GxB_LAND_EQ_UINT8      , GxB_LXOR_EQ_UINT8      , GxB_EQ_EQ_UINT8        ,
GxB_LOR_EQ_INT16       , GxB_LAND_EQ_INT16      , GxB_LXOR_EQ_INT16      , GxB_EQ_EQ_INT16        ,
GxB_LOR_EQ_UINT16      , GxB_LAND_EQ_UINT16     , GxB_LXOR_EQ_UINT16     , GxB_EQ_EQ_UINT16       ,
GxB_LOR_EQ_INT32       , GxB_LAND_EQ_INT32      , GxB_LXOR_EQ_INT32      , GxB_EQ_EQ_INT32        ,
GxB_LOR_EQ_UINT32      , GxB_LAND_EQ_UINT32     , GxB_LXOR_EQ_UINT32     , GxB_EQ_EQ_UINT32       ,
GxB_LOR_EQ_INT64       , GxB_LAND_EQ_INT64      , GxB_LXOR_EQ_INT64      , GxB_EQ_EQ_INT64        ,
GxB_LOR_EQ_UINT64      , GxB_LAND_EQ_UINT64     , GxB_LXOR_EQ_UINT64     , GxB_EQ_EQ_UINT64       ,
GxB_LOR_EQ_FP32        , GxB_LAND_EQ_FP32       , GxB_LXOR_EQ_FP32       , GxB_EQ_EQ_FP32         ,
GxB_LOR_EQ_FP64        , GxB_LAND_EQ_FP64       , GxB_LXOR_EQ_FP64       , GxB_EQ_EQ_FP64         ,

// semirings with multiply op: z = NE (x,y), where z is Boolean and x,y are given by the suffix:
GxB_LOR_NE_INT8        , GxB_LAND_NE_INT8       , GxB_LXOR_NE_INT8       , GxB_EQ_NE_INT8         ,
GxB_LOR_NE_UINT8       , GxB_LAND_NE_UINT8      , GxB_LXOR_NE_UINT8      , GxB_EQ_NE_UINT8        ,
GxB_LOR_NE_INT16       , GxB_LAND_NE_INT16      , GxB_LXOR_NE_INT16      , GxB_EQ_NE_INT16        ,
GxB_LOR_NE_UINT16      , GxB_LAND_NE_UINT16     , GxB_LXOR_NE_UINT16     , GxB_EQ_NE_UINT16       ,
GxB_LOR_NE_INT32       , GxB_LAND_NE_INT32      , GxB_LXOR_NE_INT32      , GxB_EQ_NE_INT32        ,
GxB_LOR_NE_UINT32      , GxB_LAND_NE_UINT32     , GxB_LXOR_NE_UINT32     , GxB_EQ_NE_UINT32       ,
GxB_LOR_NE_INT64       , GxB_LAND_NE_INT64      , GxB_LXOR_NE_INT64      , GxB_EQ_NE_INT64        ,
GxB_LOR_NE_UINT64      , GxB_LAND_NE_UINT64     , GxB_LXOR_NE_UINT64     , GxB_EQ_NE_UINT64       ,
GxB_LOR_NE_FP32        , GxB_LAND_NE_FP32       , GxB_LXOR_NE_FP32       , GxB_EQ_NE_FP32         ,
GxB_LOR_NE_FP64        , GxB_LAND_NE_FP64       , GxB_LXOR_NE_FP64       , GxB_EQ_NE_FP64         ,

// semirings with multiply op: z = GT (x,y), where z is Boolean and x,y are given by the suffix:
GxB_LOR_GT_INT8        , GxB_LAND_GT_INT8       , GxB_LXOR_GT_INT8       , GxB_EQ_GT_INT8         ,
GxB_LOR_GT_UINT8       , GxB_LAND_GT_UINT8      , GxB_LXOR_GT_UINT8      , GxB_EQ_GT_UINT8        ,
GxB_LOR_GT_INT16       , GxB_LAND_GT_INT16      , GxB_LXOR_GT_INT16      , GxB_EQ_GT_INT16        ,
GxB_LOR_GT_UINT16      , GxB_LAND_GT_UINT16     , GxB_LXOR_GT_UINT16     , GxB_EQ_GT_UINT16       ,
GxB_LOR_GT_INT32       , GxB_LAND_GT_INT32      , GxB_LXOR_GT_INT32      , GxB_EQ_GT_INT32        ,
GxB_LOR_GT_UINT32      , GxB_LAND_GT_UINT32     , GxB_LXOR_GT_UINT32     , GxB_EQ_GT_UINT32       ,
GxB_LOR_GT_INT64       , GxB_LAND_GT_INT64      , GxB_LXOR_GT_INT64      , GxB_EQ_GT_INT64        ,
GxB_LOR_GT_UINT64      , GxB_LAND_GT_UINT64     , GxB_LXOR_GT_UINT64     , GxB_EQ_GT_UINT64       ,
GxB_LOR_GT_FP32        , GxB_LAND_GT_FP32       , GxB_LXOR_GT_FP32       , GxB_EQ_GT_FP32         ,
GxB_LOR_GT_FP64        , GxB_LAND_GT_FP64       , GxB_LXOR_GT_FP64       , GxB_EQ_GT_FP64         ,

// semirings with multiply op: z = LT (x,y), where z is Boolean and x,y are given by the suffix:
GxB_LOR_LT_INT8        , GxB_LAND_LT_INT8       , GxB_LXOR_LT_INT8       , GxB_EQ_LT_INT8         ,
GxB_LOR_LT_UINT8       , GxB_LAND_LT_UINT8      , GxB_LXOR_LT_UINT8      , GxB_EQ_LT_UINT8        ,
GxB_LOR_LT_INT16       , GxB_LAND_LT_INT16      , GxB_LXOR_LT_INT16      , GxB_EQ_LT_INT16        ,
GxB_LOR_LT_UINT16      , GxB_LAND_LT_UINT16     , GxB_LXOR_LT_UINT16     , GxB_EQ_LT_UINT16       ,
GxB_LOR_LT_INT32       , GxB_LAND_LT_INT32      , GxB_LXOR_LT_INT32      , GxB_EQ_LT_INT32        ,
GxB_LOR_LT_UINT32      , GxB_LAND_LT_UINT32     , GxB_LXOR_LT_UINT32     , GxB_EQ_LT_UINT32       ,
GxB_LOR_LT_INT64       , GxB_LAND_LT_INT64      , GxB_LXOR_LT_INT64      , GxB_EQ_LT_INT64        ,
GxB_LOR_LT_UINT64      , GxB_LAND_LT_UINT64     , GxB_LXOR_LT_UINT64     , GxB_EQ_LT_UINT64       ,
GxB_LOR_LT_FP32        , GxB_LAND_LT_FP32       , GxB_LXOR_LT_FP32       , GxB_EQ_LT_FP32         ,
GxB_LOR_LT_FP64        , GxB_LAND_LT_FP64       , GxB_LXOR_LT_FP64       , GxB_EQ_LT_FP64         ,

// semirings with multiply op: z = GE (x,y), where z is Boolean and x,y are given by the suffix:
GxB_LOR_GE_INT8        , GxB_LAND_GE_INT8       , GxB_LXOR_GE_INT8       , GxB_EQ_GE_INT8         ,
GxB_LOR_GE_UINT8       , GxB_LAND_GE_UINT8      , GxB_LXOR_GE_UINT8      , GxB_EQ_GE_UINT8        ,
GxB_LOR_GE_INT16       , GxB_LAND_GE_INT16      , GxB_LXOR_GE_INT16      , GxB_EQ_GE_INT16        ,
GxB_LOR_GE_UINT16      , GxB_LAND_GE_UINT16     , GxB_LXOR_GE_UINT16     , GxB_EQ_GE_UINT16       ,
GxB_LOR_GE_INT32       , GxB_LAND_GE_INT32      , GxB_LXOR_GE_INT32      , GxB_EQ_GE_INT32        ,
GxB_LOR_GE_UINT32      , GxB_LAND_GE_UINT32     , GxB_LXOR_GE_UINT32     , GxB_EQ_GE_UINT32       ,
GxB_LOR_GE_INT64       , GxB_LAND_GE_INT64      , GxB_LXOR_GE_INT64      , GxB_EQ_GE_INT64        ,
GxB_LOR_GE_UINT64      , GxB_LAND_GE_UINT64     , GxB_LXOR_GE_UINT64     , GxB_EQ_GE_UINT64       ,
GxB_LOR_GE_FP32        , GxB_LAND_GE_FP32       , GxB_LXOR_GE_FP32       , GxB_EQ_GE_FP32         ,
GxB_LOR_GE_FP64        , GxB_LAND_GE_FP64       , GxB_LXOR_GE_FP64       , GxB_EQ_GE_FP64         ,

// semirings with multiply op: z = LE (x,y), where z is Boolean and x,y are given by the suffix:
GxB_LOR_LE_INT8        , GxB_LAND_LE_INT8       , GxB_LXOR_LE_INT8       , GxB_EQ_LE_INT8         ,
GxB_LOR_LE_UINT8       , GxB_LAND_LE_UINT8      , GxB_LXOR_LE_UINT8      , GxB_EQ_LE_UINT8        ,
GxB_LOR_LE_INT16       , GxB_LAND_LE_INT16      , GxB_LXOR_LE_INT16      , GxB_EQ_LE_INT16        ,
GxB_LOR_LE_UINT16      , GxB_LAND_LE_UINT16     , GxB_LXOR_LE_UINT16     , GxB_EQ_LE_UINT16       ,
GxB_LOR_LE_INT32       , GxB_LAND_LE_INT32      , GxB_LXOR_LE_INT32      , GxB_EQ_LE_INT32        ,
GxB_LOR_LE_UINT32      , GxB_LAND_LE_UINT32     , GxB_LXOR_LE_UINT32     , GxB_EQ_LE_UINT32       ,
GxB_LOR_LE_INT64       , GxB_LAND_LE_INT64      , GxB_LXOR_LE_INT64      , GxB_EQ_LE_INT64        ,
GxB_LOR_LE_UINT64      , GxB_LAND_LE_UINT64     , GxB_LXOR_LE_UINT64     , GxB_EQ_LE_UINT64       ,
GxB_LOR_LE_FP32        , GxB_LAND_LE_FP32       , GxB_LXOR_LE_FP32       , GxB_EQ_LE_FP32         ,
GxB_LOR_LE_FP64        , GxB_LAND_LE_FP64       , GxB_LXOR_LE_FP64       , GxB_EQ_LE_FP64         ,

//------------------------------------------------------------------------------
// 40 purely Boolean semirings
//------------------------------------------------------------------------------

// purely boolean semirings (in the form GxB_(add monoid)_(multipy operator)_BOOL:
GxB_LOR_FIRST_BOOL     , GxB_LAND_FIRST_BOOL    , GxB_LXOR_FIRST_BOOL    , GxB_EQ_FIRST_BOOL      , 
GxB_LOR_SECOND_BOOL    , GxB_LAND_SECOND_BOOL   , GxB_LXOR_SECOND_BOOL   , GxB_EQ_SECOND_BOOL     , 
GxB_LOR_LOR_BOOL       , GxB_LAND_LOR_BOOL      , GxB_LXOR_LOR_BOOL      , GxB_EQ_LOR_BOOL        , 
GxB_LOR_LAND_BOOL      , GxB_LAND_LAND_BOOL     , GxB_LXOR_LAND_BOOL     , GxB_EQ_LAND_BOOL       , 
GxB_LOR_LXOR_BOOL      , GxB_LAND_LXOR_BOOL     , GxB_LXOR_LXOR_BOOL     , GxB_EQ_LXOR_BOOL       , 
GxB_LOR_EQ_BOOL        , GxB_LAND_EQ_BOOL       , GxB_LXOR_EQ_BOOL       , GxB_EQ_EQ_BOOL         , 
GxB_LOR_GT_BOOL        , GxB_LAND_GT_BOOL       , GxB_LXOR_GT_BOOL       , GxB_EQ_GT_BOOL         , 
GxB_LOR_LT_BOOL        , GxB_LAND_LT_BOOL       , GxB_LXOR_LT_BOOL       , GxB_EQ_LT_BOOL         , 
GxB_LOR_GE_BOOL        , GxB_LAND_GE_BOOL       , GxB_LXOR_GE_BOOL       , GxB_EQ_GE_BOOL         , 
GxB_LOR_LE_BOOL        , GxB_LAND_LE_BOOL       , GxB_LXOR_LE_BOOL       , GxB_EQ_LE_BOOL         ; 

#endif


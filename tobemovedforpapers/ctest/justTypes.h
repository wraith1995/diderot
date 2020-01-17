/*---------- begin lib-h-head.in ----------*/
/*! \file justTypes.h
 *
 * C interface to library generated from justTypes.diderot.
 *
 * Command: /home/teocollin/gitcode/diderot/bin/diderotc --debug --log --json --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree --double --namespace=justTypes justTypes.diderot
 * Version: master:2016-07-29
 */

#ifndef _JUSTTYPES_H_
#define _JUSTTYPES_H_

#define DIDEROT_DOUBLE_PRECISION
#define DIDEROT_INT
#define DIDEROT_TARGET_SEQUENTIAL

#include "diderot/config.h"

#if defined(HAVE_STDBOOL_H)
#  include <stdbool.h>
#elif !defined(__bool_true_false_are_defined)
#  ifndef HAVE__BOOL
#    ifdef __cplusplus
       typedef bool _Bool;
#    else
#      define _Bool signed char
#    endif
#  endif
#  define bool _Bool
#  define false 0
#  define true 1
#  define __bool_true_false_are_defined 1
#endif

#include <stdint.h>
#include <string.h>
#include "teem/nrrd.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct justTypes_struct_world justTypes_world_t;
/*---------- end lib-h-head.in ----------*/

/**** Functions etc. for input variables ****/
bool justTypes_input_set_a (justTypes_world_t *wrld, void *v);
bool justTypes_input_set_b (justTypes_world_t *wrld, void *v);
bool justTypes_input_set_c (justTypes_world_t *wrld, void *v);
/*---------- begin lib-h-body.in ----------*/

/***** World query operations *****/

//! Return the total number of strands (active+stable) in the world
uint32_t justTypes_num_strands (justTypes_world_t *wrld);

//! Return the total number of active strands
uint32_t justTypes_num_active_strands (justTypes_world_t *wrld);

//! Return the total number of stable strands
uint32_t justTypes_num_stable_strands (justTypes_world_t *wrld);

//! Return true if there are any recorded error conditions
bool justTypes_any_errors (justTypes_world_t *wrld);

//! Return the pending error message (if any).  This call clears the pending error
//! state.
char *justTypes_get_errors (justTypes_world_t *wrld);

/***** Program running operations *****/

//! Allocate the program's world
//! \return the new world or NULL if there are any errors
justTypes_world_t *justTypes_new_world ();

//! Initialize the execution state for the world.  This includes allocating processor
//! and GPU resources for parallel execution.
//! \param wrld the world-state of the Diderot program
//! \return true if there are any errors
bool justTypes_init_world (justTypes_world_t *wrld);

//! Initiaize the globals and create the initial set of strands
//! \param wrld the world-state of the Diderot program
//! \return true if there are any errors
bool justTypes_create_strands (justTypes_world_t *wrld);

//! Run the Diderot program
//! \param wrld the world-state of the Diderot program
//! \param maxNSteps the limit on the number of super steps; 0 means unlimited
//! \return the number of steps taken.
uint32_t justTypes_run (justTypes_world_t *wrld, uint32_t maxNSteps);

//! shutdown and deallocate the world
//! \param wrld the world-state of the Diderot program
void justTypes_shutdown (justTypes_world_t *wrld);

/***** Runtime options *****/

//! Set verbose mode
//! \param wrld the world-state of the Diderot program
//! \param mode the mode value to set; true means verbose
void justTypes_set_verbose (justTypes_world_t *wrld, bool mode);

//! Get verbose mode
//! \return true if there are any errors
bool justTypes_get_verbose (justTypes_world_t *wrld);

//! Register a callback function for Diderot print() calls
//! \param wrld the world-state of the Diderot program
//! \param pr the printing callback function
//! \param data an opaque data value that will be passed to the
//!             printer
//! \return true if there are any errors
bool justTypes_set_printer_cb (justTypes_world_t *wrld, bool (*pr)(void *, char *), void *data);
/*---------- end lib-h-body.in ----------*/

/**** Getters for output values ****/
bool justTypes_output_get_pos (justTypes_world_t *wrld, Nrrd *data);
bool justTypes_output_get__pos (justTypes_world_t *wrld, Nrrd *data);
/*---------- begin lib-h-foot.in ----------*/

#ifdef __cplusplus
}
#endif

#endif /* !_JUSTTYPES_H_ */
/*---------- end lib-h-foot.in ----------*/

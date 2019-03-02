/*---------- begin lib-h-head.in ----------*/
/*! \file 15.h
 *
 * C interface to library generated from 15.diderot.
 *
 * Command: /home/teocollin/installs/diderot/diderot/bin/diderotc --log --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree 15.diderot
 * Version: master:2016-07-29
 */

#ifndef _15_H_
#define _15_H_

#define DIDEROT_SINGLE_PRECISION
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

typedef struct Diderot_struct_world Diderot_world_t;
/*---------- end lib-h-head.in ----------*/

/**** Functions etc. for input variables ****/
void Diderot_input_get_atest (Diderot_world_t *wrld, float v[15]);
bool Diderot_input_set_atest (Diderot_world_t *wrld, float v[15]);
/*---------- begin lib-h-body.in ----------*/

/***** World query operations *****/

//! Return the total number of strands (active+stable) in the world
uint32_t Diderot_num_strands (Diderot_world_t *wrld);

//! Return the total number of active strands
uint32_t Diderot_num_active_strands (Diderot_world_t *wrld);

//! Return the total number of stable strands
uint32_t Diderot_num_stable_strands (Diderot_world_t *wrld);

//! Return true if there are any recorded error conditions
bool Diderot_any_errors (Diderot_world_t *wrld);

//! Return the pending error message (if any).  This call clears the pending error
//! state.
char *Diderot_get_errors (Diderot_world_t *wrld);

/***** Program running operations *****/

//! Allocate the program's world
//! \return the new world or NULL if there are any errors
Diderot_world_t *Diderot_new_world ();

//! Initialize the execution state for the world.  This includes allocating processor
//! and GPU resources for parallel execution.
//! \param wrld the world-state of the Diderot program
//! \return true if there are any errors
bool Diderot_init_world (Diderot_world_t *wrld);

//! Initiaize the globals and create the initial set of strands
//! \param wrld the world-state of the Diderot program
//! \return true if there are any errors
bool Diderot_create_strands (Diderot_world_t *wrld);

//! Run the Diderot program
//! \param wrld the world-state of the Diderot program
//! \param maxNSteps the limit on the number of super steps; 0 means unlimited
//! \return the number of steps taken.
uint32_t Diderot_run (Diderot_world_t *wrld, uint32_t maxNSteps);

//! shutdown and deallocate the world
//! \param wrld the world-state of the Diderot program
void Diderot_shutdown (Diderot_world_t *wrld);

/***** Runtime options *****/

//! Set verbose mode
//! \param wrld the world-state of the Diderot program
//! \param mode the mode value to set; true means verbose
void Diderot_set_verbose (Diderot_world_t *wrld, bool mode);

//! Get verbose mode
//! \return true if there are any errors
bool Diderot_get_verbose (Diderot_world_t *wrld);

//! Register a callback function for Diderot print() calls
//! \param wrld the world-state of the Diderot program
//! \param pr the printing callback function
//! \param data an opaque data value that will be passed to the
//!             printer
//! \return true if there are any errors
bool Diderot_set_printer_cb (Diderot_world_t *wrld, bool (*pr)(void *, char *), void *data);
/*---------- end lib-h-body.in ----------*/

/**** Getters for output values ****/
bool Diderot_output_get_umm (Diderot_world_t *wrld, Nrrd *data);
bool Diderot_output_get_temp (Diderot_world_t *wrld, Nrrd *data);
/*---------- begin lib-h-foot.in ----------*/

#ifdef __cplusplus
}
#endif

#endif /* !_15_H_ */
/*---------- end lib-h-foot.in ----------*/

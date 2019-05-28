/*---------- begin lib-h-head.in ----------*/
/*! \file evalProg.h
 *
 * C interface to library generated from evalProg.diderot.
 *
 * Command: /home/teocollin/gitcode/diderot/bin/diderotc --debug --log --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree --double --namespace=evalProg --target=parallel evalProg.diderot
 * Version: master:2016-07-29
 */

#ifndef _EVALPROG_H_
#define _EVALPROG_H_

#define DIDEROT_DOUBLE_PRECISION
#define DIDEROT_INT
#define DIDEROT_TARGET_PARALLEL

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

typedef struct evalProg_struct_world evalProg_world_t;
/*---------- end lib-h-head.in ----------*/

/**** Functions etc. for input variables ****/
bool evalProg_input_set_meshData (evalProg_world_t *wrld, void *v);
bool evalProg_input_set_space (evalProg_world_t *wrld, void *v);
bool evalProg_input_set_data (evalProg_world_t *wrld, void *v);
bool evalProg_input_set_by_name_ipos (evalProg_world_t *wrld, const char *s);
bool evalProg_input_set_ipos (evalProg_world_t *wrld, Nrrd *nin);
/*---------- begin lib-h-body.in ----------*/

/***** World query operations *****/

//! Return the total number of strands (active+stable) in the world
uint32_t evalProg_num_strands (evalProg_world_t *wrld);

//! Return the total number of active strands
uint32_t evalProg_num_active_strands (evalProg_world_t *wrld);

//! Return the total number of stable strands
uint32_t evalProg_num_stable_strands (evalProg_world_t *wrld);

//! Return true if there are any recorded error conditions
bool evalProg_any_errors (evalProg_world_t *wrld);

//! Return the pending error message (if any).  This call clears the pending error
//! state.
char *evalProg_get_errors (evalProg_world_t *wrld);

/***** Program running operations *****/

//! Allocate the program's world
//! \return the new world or NULL if there are any errors
evalProg_world_t *evalProg_new_world ();

//! Initialize the execution state for the world.  This includes allocating processor
//! and GPU resources for parallel execution.
//! \param wrld the world-state of the Diderot program
//! \return true if there are any errors
bool evalProg_init_world (evalProg_world_t *wrld);

//! Initiaize the globals and create the initial set of strands
//! \param wrld the world-state of the Diderot program
//! \return true if there are any errors
bool evalProg_create_strands (evalProg_world_t *wrld);

//! Run the Diderot program
//! \param wrld the world-state of the Diderot program
//! \param maxNSteps the limit on the number of super steps; 0 means unlimited
//! \return the number of steps taken.
uint32_t evalProg_run (evalProg_world_t *wrld, uint32_t maxNSteps);

//! shutdown and deallocate the world
//! \param wrld the world-state of the Diderot program
void evalProg_shutdown (evalProg_world_t *wrld);

/***** Runtime options *****/

//! Set verbose mode
//! \param wrld the world-state of the Diderot program
//! \param mode the mode value to set; true means verbose
void evalProg_set_verbose (evalProg_world_t *wrld, bool mode);

//! Get verbose mode
//! \return true if there are any errors
bool evalProg_get_verbose (evalProg_world_t *wrld);

//! Register a callback function for Diderot print() calls
//! \param wrld the world-state of the Diderot program
//! \param pr the printing callback function
//! \param data an opaque data value that will be passed to the
//!             printer
//! \return true if there are any errors
bool evalProg_set_printer_cb (evalProg_world_t *wrld, bool (*pr)(void *, char *), void *data);
/*---------- end lib-h-body.in ----------*/
/*---------- begin lib-h-par-extras.in ----------*/
//! get the number of hardware cores
//! \param wrld the world-state of the Diderot program
//! \return the number of cores on the system
uint32_t evalProg_get_num_cores (evalProg_world_t *wrld);

//! set the number of workers.  The value should be between 0 and the number of
//! hardware cores.  Note that this function should be called after evalProg_init_world
//! and before evalProg_create_strands.
//! \param wrld the world-state of the Diderot program
//! \param nWorkers the requested number of workers; 0 means set the number of
//!        workers to the number of cores.
//! \return true if there are any errors
bool evalProg_set_num_workers (evalProg_world_t *wrld, uint32_t nWorkers);

//! get the number of workers.
//! \param wrld the world-state of the Diderot program
//! \return the number of workers
uint32_t evalProg_get_num_workers (evalProg_world_t *wrld);
/*---------- end lib-h-par-extras.in ----------*/

/**** Getters for output values ****/
bool evalProg_output_get_normal (evalProg_world_t *wrld, Nrrd *data);
/*---------- begin lib-h-foot.in ----------*/

#ifdef __cplusplus
}
#endif

#endif /* !_EVALPROG_H_ */
/*---------- end lib-h-foot.in ----------*/

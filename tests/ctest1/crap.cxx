/*---------- begin cxx-head.in ----------*/
/*! \file crap.cxx
 *
 * Generated from crap.diderot.
 *
 * Command: ../bin/diderotc --exec --log --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree crap.diderot
 * Version: master:2016-07-29
 */
/*---------- end cxx-head.in ----------*/

#define DIDEROT_STRAND_HAS_CONSTR
#define DIDEROT_HAS_CONSTS
#define DIDEROT_NO_INPUTS
#define DIDEROT_HAS_KILL_ALL
/*---------- begin exec-incl.in ----------*/
#define DIDEROT_STANDALONE_EXEC
#define DIDEROT_SINGLE_PRECISION
#define DIDEROT_INT
#define DIDEROT_TARGET_SEQUENTIAL
#include "diderot/diderot.hxx"

#ifdef DIDEROT_ENABLE_LOGGING
#define IF_LOGGING(...)         __VA_ARGS__
#else
#define IF_LOGGING(...)
#endif
/*---------- end exec-incl.in ----------*/

// ***** Begin synthesized types *****

namespace diderot {
    template <>
    struct dynseq_traits<int32_t> {
        using value_type = int32_t;
        using base_type = int32_t;
        static const __details::load_fn_ptr<base_type> *load_fn_tbl;
        static const uint32_t values_per_elem = 1;
    };
    const __details::load_fn_ptr< dynseq_traits< int32_t >::base_type > *dynseq_traits< int32_t >::load_fn_tbl = nrrdILoad;
    template <>
    struct dynseq_traits<float> {
        using value_type = float;
        using base_type = float;
        static const __details::load_fn_ptr<base_type> *load_fn_tbl;
        static const uint32_t values_per_elem = 1;
    };
    const __details::load_fn_ptr< dynseq_traits< float >::base_type > *dynseq_traits< float >::load_fn_tbl = nrrdFLoad;
} // namespace diderot
// ***** End synthesized types *****

/*---------- begin namespace-open.in ----------*/
namespace Diderot {

static std::string ProgramName = "crap";

struct world;
struct gg_strand;
/*---------- end namespace-open.in ----------*/

/*---------- begin nrrd-save-helper.in ----------*/
/* helper function for saving output to nrrd file */
inline bool nrrd_save_helper (std::string const &file, Nrrd *nin)
{
    if (nrrdSave (file.c_str(), nin, nullptr)) {
        std::cerr << "Error saving \"" << file << "\":\n" << biffGetDone(NRRD) << std::endl;
        return true;
    }
    else {
        return false;
    }
}

/* Helper function for saving dynamic sequence output to a nrrd file.
 * Dynamic sequence output is represented by two nrrds: one for lengths
 * and one for the actual sequence data.  If the sum of the lengths is
 * zero, then the second nrrd will be empty and should not be output,
 * since teem barfs on empty nrrds.
 */
inline bool dynseq_save_helper (std::string const &prefix, Nrrd *lens, Nrrd *data)
{
    if (nrrd_save_helper(prefix + "-len.nrrd", lens)) {
        return true;
    }
  // check for an empty data file; in this case nrrdEmpty will have been called
  // on the data nrrd, so we can just check for a 0-dimension object.
    if (data->dim == 0) {
        std::cerr << "Warning: all sequences in output are empty, so no '"
            << prefix << "-data.nrrd' file produced\n";
        return false;
    }
  // write the data
    return nrrd_save_helper(prefix + "-data.nrrd", data);
}

#ifdef DIDEROT_EXEC_SNAPSHOT
// version of dynseq_save_helper for snapshots
inline bool dynseq_save_helper (
    std::string const &prefix,
    std::string const &suffix,
    Nrrd *lens,
    Nrrd *data)
{
    if (nrrd_save_helper(prefix + "-len" + suffix + ".nrrd", lens)) {
        return true;
    }
  // check for an empty data file; in this case nrrdEmpty will have been called
  // on the data nrrd, so we can just check for a 0-dimension object.
    std::string dataFile = prefix + "-data" + suffix + ".nrrd";
    if (data->dim == 0) {
        std::cerr << "Warning: all sequences in snapshot are empty, so no '"
            << dataFile << "' file produced\n";
        return false;
    }
  // write the data
    return nrrd_save_helper(dataFile, data);
}
#endif // !DIDEROT_EXEC_SNAPSHOT
/*---------- end nrrd-save-helper.in ----------*/

struct globals {
    diderot::array< int32_t, 2 > gv_ll;
    diderot::dynseq< int32_t > gv_a1;
    diderot::dynseq< float > gv_a2;
    ~globals () { }
};
struct gg_strand {
    float sv_result;
    int32_t sv_i;
    int32_t sv_j;
};
/*---------- begin seq-sarr.in ----------*/
// forward declarations of strand methods
#ifdef DIDEROT_HAS_START_METHOD
static diderot::strand_status gg_start (gg_strand *self);
#endif // DIDEROT_HAS_START_METHOD
static diderot::strand_status gg_update (world *wrld, globals *glob, gg_strand *self);
#ifdef DIDEROT_HAS_STABILIZE_METHOD
static void gg_stabilize (gg_strand *self);
#endif // DIDEROT_HAS_STABILIZE_METHOD

// if we have both communication and "die", then we need to track when strands die
// so that we can rebuild the list of strands use to construct the kd-tree
#if defined(DIDEROT_HAS_STRAND_COMMUNICATION) && !defined(DIDEROT_HAS_STRAND_DIE)
#  define TRACK_STRAND_DEATH
#endif

// strand_array for SEQUENTIAL/NO BSP/SINGLE STATE/DIRECT ACCESS
//
struct strand_array {
    typedef gg_strand strand_t;
    typedef uint32_t index_t;
    typedef index_t sid_t;              // strand ID (index into strand-state storage)

    uint8_t             *_status;       // the array of status information for the strands
    char                *_storage;      // points to array of gg_strand structs
    uint32_t            _nItems;        // number of items in the _storage and _status arrays
    uint32_t            _nStable;       // number of stable strands
    uint32_t            _nActive;       // number of active strands
    uint32_t            _nFresh;        // number of fresh strands (new strands from create_strands)
#ifdef TRACK_STRAND_DEATH
    bool                _died;          // a strand died in the current superstep.
#endif

    strand_array () : _status(nullptr), _storage(nullptr), _nItems(0) { }
    ~strand_array ();

    uint32_t in_state_index () const { return 0; /* dummy */ }

    uint32_t num_active () const { return this->_nActive; }
    uint32_t num_stable () const { return this->_nStable; }
    uint32_t num_alive () const { return this->_nActive+this->_nStable; }

  // return the ID of a strand, which is the same as the ix index
    sid_t id (index_t ix) const
    {
        assert (ix < this->_nItems);
        return ix;
    }
  // return a pointer to the strand with the given ID
    gg_strand *id_to_strand (sid_t id) const
    {
        assert (id < this->_nItems);
        return reinterpret_cast<gg_strand *>(this->_storage + id * sizeof(gg_strand));
    }

  // return a strand's status
    diderot::strand_status status (index_t ix) const
    {
        assert (ix < this->_nItems);
        return static_cast<diderot::strand_status>(this->_status[ix]);
    }
  // return a pointer to the given strand
    gg_strand *strand (index_t ix) const
    {
        return this->id_to_strand(this->id(ix));
    }
  // return a pointer to the local state of strand ix
    gg_strand *local_state (index_t ix) const
    {
        return this->strand(ix);
    }
  // return a pointer to the local state of strand with the given ID
    gg_strand *id_to_local_state (sid_t id) const
    {
        return this->id_to_strand(id);
    }

  // is an index valid for the strand array?
    bool validIndex (index_t ix) const { return (ix < this->_nItems); }

  // is a given strand alive?
    bool isAlive (index_t ix) const
    {
#ifdef DIDEROT_HAS_STRAND_DIE
        return aliveSts(this->status(ix));
#else
        return true;
#endif
    }

  // allocate space for nItems
    bool alloc (uint32_t nItems)
    {
        this->_storage = static_cast<char *>(std::malloc (nItems * sizeof(gg_strand)));
        if (this->_storage == nullptr) {
            return true;
        }
        this->_status = static_cast<uint8_t *>(std::malloc (nItems * sizeof(uint8_t)));
        if (this->_status == nullptr) {
            std::free (this->_storage);
            return true;
        }
        this->_nItems = nItems;
        this->_nActive = 0;
        this->_nStable = 0;
        this->_nFresh = 0;
        return false;
    }

  // initialize the first nStrands locations as new active strands
    void create_strands (uint32_t nStrands)
    {
        assert (this->_nActive == 0);
        assert (this->_nItems == nStrands);
        for (index_t ix = 0;  ix < nStrands;  ix++) {
            this->_status[ix] = diderot::kActive;
            new (this->strand(ix)) gg_strand;
        }
        this->_nActive = nStrands;
        this->_nFresh = nStrands;
#ifdef TRACK_STRAND_DEATH
        this->_died = false;
#endif
    }

  // swap in and out states (NOP for this version)
    void swap () { }

#ifdef DIDEROT_HAS_START_METHOD
  // invoke strand's start method
    diderot::strand_status strand_start (index_t ix)
    {
        return gg_start(this->strand(ix));
    }
#endif // DIDEROT_HAS_START_METHOD

  // invoke strand's update method
    diderot::strand_status strand_update (world *wrld, globals *glob, index_t ix)
    {
        return gg_update(wrld, glob, this->strand(ix));
    }

  // invoke strand's stabilize method
    index_t strand_stabilize (index_t ix)
    {
#ifdef DIDEROT_HAS_STABILIZE_METHOD
        gg_stabilize (this->strand(ix));
#endif // DIDEROT_HAS_STABILIZE_METHOD
        this->_status[ix] = diderot::kStable;
        this->_nActive--;
        this->_nStable++;
      // skip to next active strand
        do {
            ix++;
        } while ((ix < this->_nItems) && notActiveSts(this->status(ix)));
        return ix;
    }

  // mark the given strand as dead
    index_t kill (index_t ix)
    {
#ifdef TRACK_STRAND_DEATH
        this->_died = true;
#endif
        this->_status[ix] = diderot::kDead;
        this->_nActive--;
      // skip to next active strand
        do {
            ix++;
        } while ((ix < this->_nItems) && notActiveSts(this->status(ix)));
        return ix;
    }

  // finish the local-phase of a superstep (NOP)
#ifdef TRACK_STRAND_DEATH
    bool finish_step ()
    {
        bool res = this->_died;
        this->_died = false;
        return res;
    }
#else
    bool finish_step () { return false; }
#endif

  // finish a kill_all operation (NOP)
    void finish_kill_all () { }

  // finish a stabilize_all operation (NOP)
    void finish_stabilize_all () { }

    index_t begin_alive () const
    {
        index_t ix = 0;
#ifdef DIDEROT_HAS_STRAND_DIE
        while ((ix < this->_nItems) && notAliveSts(this->status(ix))) {
            ix++;
        }
#endif
        return ix;
    }
    index_t end_alive () const { return this->_nItems; }
    index_t next_alive (index_t &ix) const
    {
        ix++;
#ifdef DIDEROT_HAS_STRAND_DIE
        while ((ix < this->_nItems) && notAliveSts(this->status(ix))) {
            ix++;
        }
#endif
        return ix;
    }

  // iterator over active strands
    index_t begin_active () const
    {
        index_t ix = 0;
        while ((ix < this->_nItems) && notActiveSts(this->status(ix))) {
            ix++;
        }
        return ix;
    }
    index_t end_active () const { return this->_nItems; }
    index_t next_active (index_t &ix) const
    {
        do {
            ix++;
        } while ((ix < this->_nItems) && notActiveSts(this->status(ix)));
        return ix;
    }

  // iterator over stable strands
    index_t begin_stable () const
    {
        index_t ix = 0;
        while ((ix < this->_nItems) && (this->status(ix) != diderot::kStable)) {
            ix++;
        }
        return ix;
    }
    index_t end_stable () const { return this->_nItems; }
    index_t next_stable (index_t &ix) const
    {
        do {
            ix++;
        } while ((ix < this->_nItems) && (this->status(ix) != diderot::kStable));
        return ix;
    }

  // iterator over fresh strands; since the only new strands were created by create_strand
  // we iterate over all of them
    index_t begin_fresh () const { return 0; }
    index_t end_fresh () const { return this->_nFresh; }
    index_t next_fresh (index_t &ix) const { return ++ix; }

}; // struct strand_array

strand_array::~strand_array ()
{
  // run destructors to reclaim any dynamic memory attached to the strand state
    for (auto ix = this->begin_alive();  ix != this->end_alive();  ix = this->next_alive(ix)) {
        this->strand(ix)->~gg_strand();
    }
    if (this->_status != nullptr) std::free (this->_status);
    if (this->_storage != nullptr) std::free (this->_storage);
}
/*---------- end seq-sarr.in ----------*/

struct world : public diderot::world_base {
    strand_array _strands;
    globals *_globals;
    world ();
    ~world ();
    bool init ();
    bool alloc (int32_t base[1], uint32_t size[1]);
    bool create_strands ();
    uint32_t run (uint32_t max_nsteps);
    void swap_state ();
    void kill_all ();
};
// ***** Begin synthesized operations *****

// ***** End synthesized operations *****

static std::string OutputFile = "result.nrrd";
static void register_outputs (diderot::options *opts)
{
    opts->add("o,output", "specify output file", &OutputFile, true);
}
static void init_consts (world *wrld)
{
    globals *glob = wrld->_globals;
    glob->gv_ll[0] = 2;
    glob->gv_ll[1] = 2;
}
static bool init_globals (world *wrld)
{
    diderot::dynseq< int32_t > l_accum_0;
    diderot::dynseq< int32_t > l_accum_2;
    diderot::dynseq< float > l_accum_4;
    diderot::dynseq< float > l_accum_6;
    globals *glob = wrld->_globals;
    l_accum_0 = {};
    for (int32_t i_j_1 = 0; i_j_1 <= 10; ++i_j_1) {
        l_accum_0 = diderot::dynseq< int32_t >::append(l_accum_0, 1 + i_j_1);
    }
    l_accum_2 = {};
    for (auto it_0 = l_accum_0.cbegin(); it_0 != l_accum_0.cend(); ++it_0) {
        auto i_k_3 = *it_0;
        l_accum_2 = diderot::dynseq< int32_t >::append(l_accum_2, i_k_3 + 2);
    }
    glob->gv_a1 = l_accum_2;
    l_accum_4 = {};
    for (int32_t i_j_5 = 0; i_j_5 <= 10; ++i_j_5) {
        l_accum_4 = diderot::dynseq< float >::append(l_accum_4, 0.1e1f + static_cast<float>(i_j_5));
    }
    l_accum_6 = {};
    for (auto it_1 = l_accum_4.cbegin(); it_1 != l_accum_4.cend(); ++it_1) {
        auto i_k_7 = *it_1;
        l_accum_6 = diderot::dynseq< float >::append(l_accum_6, i_k_7 + 0.2e1f);
    }
    glob->gv_a2 = l_accum_6;
    return false;
}
static void gg_init (gg_strand *self, int32_t p_i_8, int32_t p_j_9)
{
    self->sv_result = 0.1e1f * static_cast<float>(p_i_8);
    self->sv_i = p_i_8;
    self->sv_j = p_j_9;
}
static diderot::strand_status gg_update (world *wrld, globals *glob, gg_strand *self)
{
    float l__t_10 = self->sv_result * static_cast<float>(self->sv_i * 10);
    wrld->print() << glob->gv_ll[0] << " ! " << glob->gv_ll[1] << "bah\n" << std::flush;
    wrld->print() << glob->gv_a1 << std::flush;
    wrld->print() << l__t_10 << "\n" << "numCell:" << 2 << "\n" << glob->gv_a1 << glob->gv_a2 << "\n" << 0.2e1f << "\n" << self->sv_j << "\n" << std::flush;
    self->sv_result = l__t_10;
    return diderot::kStabilize;
}
bool output_get_result (world *wrld, Nrrd *nData)
{
    // Compute sizes of nrrd file
    size_t sizes[1];
    sizes[0] = wrld->_strands.num_stable();
    // Allocate nData nrrd
    if (nrrdMaybeAlloc_nva(nData, nrrdTypeFloat, 1, sizes) != 0) {
        char *msg = biffGetDone(NRRD);
        biffMsgAdd(wrld->_errors, msg);
        std::free(msg);
        return true;
    }
    // copy data to output nrrd
    char *cp = reinterpret_cast<char *>(nData->data);
    for (auto ix = wrld->_strands.begin_stable(); ix != wrld->_strands.end_stable(); ix = wrld->_strands.next_stable(
        ix)) {
        memcpy(cp, &wrld->_strands.strand(ix)->sv_result, 1 * sizeof(float));
        cp += 1 * sizeof(float);
    }
    nData->axis[0].kind = nrrdKindList;
    return false;
}
static bool write_output (world *wrld)
{
    Nrrd *nData;
    nData = nrrdNew();
    if (output_get_result(wrld, nData)) {
        wrld->error("Error getting nrrd data for \'result\'");
        return true;
    }
    else if (nrrd_save_helper(OutputFile, nData)) {
        return true;
    }
    nrrdNuke(nData);
    return false;
}
/*---------- begin world-methods.in ----------*/
// Allocate the program's world
//
world::world ()
    : diderot::world_base (ProgramName, false, 1)
{
#ifndef DIDEROT_NO_GLOBALS
    this->_globals = new globals;
#endif

#ifdef DIDEROT_HAS_STRAND_COMMUNICATION
    this->_tree = nullptr;
#endif
} // world constructor

// shutdown and deallocate the world
//
world::~world ()
{
#ifndef DIDEROT_NO_GLOBALS
    delete this->_globals;
#endif

#ifdef DIDEROT_HAS_STRAND_COMMUNICATION
    delete this->_tree;
#endif

} // world destructor

// Initialize the program's world
//
bool world::init ()
{
    assert (this->_stage == diderot::POST_NEW);

#if !defined(DIDEROT_STANDALONE_EXEC) && !defined(DIDEROT_NO_INPUTS)
  // initialize the defined flags for the input globals
    init_defined_inputs (this);
#endif

#ifdef DIDEROT_TARGET_PARALLEL
  // get CPU info
    if (this->_sched->get_cpu_info (this)) {
        return true;
    }
#endif

#ifdef DIDEROT_HAS_CONSTS
    init_consts (this);
#endif

    this->_stage = diderot::POST_INIT;

    return false;

}

// allocate the initial strands and initialize the rest of the world structure.
//
bool world::alloc (int32_t base[1], uint32_t size[1])
{
    size_t numStrands = 1;
    for (uint32_t i = 0;  i < 1;  i++) {
        numStrands *= size[i];
        this->_base[i] = base[i];
        this->_size[i] = size[i];
    }

    if (this->_verbose) {
        std::cerr << "world::alloc: " << size[0];
        for (uint32_t i = 1;  i < 1;  i++) {
            std::cerr << " x " << size[i];
        }
        std::cerr << std::endl;
    }

#ifdef DIDEROT_TARGET_PARALLEL
  // determine the block size based on the initial number of strands and the
  // number of workers
    this->_strands.set_block_size (this->_sched->_numWorkers, numStrands);
#endif

  // allocate the strand array
    if (this->_strands.alloc (numStrands)) {
        biffMsgAdd (this->_errors, "unable to allocate strand-state array\n");
        return true;
    }

  // initialize strand state pointers etc.
    this->_strands.create_strands (numStrands);

#ifdef DIDEROT_HAS_STRAND_COMMUNICATION
    this->_tree = new diderot::kdtree<0, float, strand_array> (&this->_strands);
#endif

    return false;

} // world::alloc

// swap input and output states
//
inline void world::swap_state ()
{
    this->_strands.swap ();
}

#ifdef DIDEROT_HAS_KILL_ALL
void world::kill_all ()
{
    if (this->_strands.num_active() > 0) {
        for (auto ix = this->_strands.begin_active();
            ix != this->_strands.end_active();
            )
        {
            assert (this->_strands.status(ix) == diderot::kActive);
            ix = this->_strands.kill (ix);
        }
        this->_strands.finish_kill_all();
    }
    assert (this->_strands.num_active() == 0);
}
#endif

#ifdef DIDEROT_HAS_STABILIZE_ALL
void world::stabilize_all ()
{
#ifndef DIDEROT_NO_GLOBALS
    globals *glob = this->_globals;
#endif

    if (this->_strands.num_active() > 0) {
        for (auto ix = this->_strands.begin_active();
            ix != this->_strands.end_active();
            )
        {
            assert (this->_strands.status(ix) == diderot::kActive);
            this->_strands._status[ix] = diderot::kStable;
            ix = this->_strands.strand_stabilize (ix);
        }
        this->_strands.finish_stabilize_all();
    }
    assert (this->_strands.num_active() == 0);
}
#endif
/*---------- end world-methods.in ----------*/

bool world::create_strands ()
{
    if (init_globals(this)) {
        return true;
    }
    globals *glob = this->_globals;
    diderot::dynseq< int32_t > seq_2 = glob->gv_a1;
    int lo_3 = 5;
    int hi_4 = 11;
    int32_t base[1] = {0,};
    uint32_t size[1] = {static_cast<uint32_t>(seq_2.length() * (hi_4 - lo_3 + 1)),};
    if (this->alloc(base, size)) {
        return true;
    }
    uint32_t ix = 0;
    for (auto it_5 = seq_2.cbegin(); it_5 != seq_2.cend(); ++it_5) {
        auto i_j_11 = *it_5;
        for (int i_k_12 = lo_3; i_k_12 <= hi_4; i_k_12++) {
            gg_init(this->_strands.strand(ix), i_j_11, i_k_12);
            ++ix;
        }
    }
    this->swap_state();
    this->_stage = diderot::POST_CREATE;
    return false;
}
/*---------- begin seq-run-nobsp.in ----------*/
//! Run the Diderot program (sequential version without BSP semantics)
//! \param max_nsteps the limit on the number of super steps; 0 means unlimited
//! \return the number of steps taken, or 0 if done or there is an error.
uint32_t world::run (uint32_t max_nsteps)
{
    if (this->_stage == diderot::POST_CREATE) {
#ifdef DIDEROT_HAS_GLOBAL_START
        this->global_start();
#endif
        this->_stage = diderot::RUNNING;
    }
    else if (this->_stage == diderot::DONE) {
        return 0;
    }
    assert (this->_stage == diderot::RUNNING);

#ifndef DIDEROT_NO_GLOBALS
    globals *glob = this->_globals;
#endif

    if (max_nsteps == 0) {
        max_nsteps = 0xffffffff;  // essentially unlimited
    }

    double t0 = airTime();

    if (this->_verbose) {
        std::cerr << "run with " << this->_strands.num_alive() << " strands ..." << std::endl;
    }

#ifdef DIDEROT_HAS_START_METHOD
    this->run_start_methods();
#endif

  // iterate until all strands are stable
    uint32_t maxSteps = 0;
    for (auto ix = this->_strands.begin_active();
         ix != this->_strands.end_active();
         )
    {
        diderot::strand_status sts = this->_strands.status(ix);
        uint32_t nSteps = 0;
        while ((! sts) && (nSteps < max_nsteps)) {
            nSteps++;
            sts = this->_strands.strand_update(this, glob, ix);
        }
        switch (sts) {
          case diderot::kStabilize:
          // stabilize the strand's state.
            ix = this->_strands.strand_stabilize (ix);
            break;
#ifdef DIDEROT_HAS_STRAND_DIE
          case diderot::kDie:
            ix = this->_strands.kill (ix);
            break;
#endif
          default:
            assert (sts == this->_strands.status(ix));
	    ix = this->_strands.next_active(ix);
            break;
        }
        if (maxSteps < nSteps) maxSteps = nSteps;
    }

    this->_run_time += airTime() - t0;

    if (this->_strands.num_active() == 0) {
        this->_stage = diderot::DONE;
    }

    return maxSteps;

} // world::run
/*---------- end seq-run-nobsp.in ----------*/

/*---------- begin namespace-close.in ----------*/

} // namespace Diderot
/*---------- end namespace-close.in ----------*/

/*---------- begin exit-with-error.in ----------*/
/* helper function that reports an error, deletes the world, and then exits */
#ifdef HAVE_FUNC_ATTRIBUTE_NORETURN
void exit_with_error (Diderot::world *wrld, std::string const &msg) __attribute__ ((noreturn));
#endif
void exit_with_error (Diderot::world *wrld, std::string const &msg)
{
    std::cerr << msg << ":\n" << wrld->get_errors() << std::endl;
    delete wrld;
    exit (1);
}
/*---------- end exit-with-error.in ----------*/

/*---------- begin seq-main.in ----------*/
using namespace Diderot;

//! Main function for standalone sequential C target
//
int main (int argc, const char **argv)
{
    bool        timingFlg = false;      //! true if timing computation
    uint32_t    stepLimit = 0;          //! limit on number of execution steps (0 means unlimited)
    std::string printFile = "-";        //! file to direct printed output into
#ifdef DIDEROT_EXEC_SNAPSHOT
    uint32_t    snapshotPeriod = 0;     //! supersteps per snapshot; 0 means no snapshots
#endif
    uint32_t    nSteps = 0;             //! number of supersteps taken

  // create the world
    world *wrld = new (std::nothrow) world();
    if (wrld == nullptr) {
        exit_with_error (wrld, "Error: unable to create world");
    }

#ifndef DIDEROT_NO_INPUTS
  // initialize the default values for the inputs
    cmd_line_inputs inputs;
    init_defaults (&inputs);
#endif

  // handle command-line options
    {
        diderot::options *opts = new diderot::options ();
        opts->add ("l,limit", "specify limit on number of super-steps (0 means unlimited)",
            &stepLimit, true);
#ifdef DIDEROT_EXEC_SNAPSHOT
        opts->add ("s,snapshot",
            "specify number of super-steps per snapshot (0 means no snapshots)",
            &snapshotPeriod, true);
#endif
        opts->add ("print", "specify where to direct printed output", &printFile, true);
        opts->addFlag ("v,verbose", "enable runtime-system messages", &(wrld->_verbose));
        opts->addFlag ("t,timing", "enable execution timing", &timingFlg);
#ifndef DIDEROT_NO_INPUTS
      // register options for setting global inputs
        register_inputs (&inputs, opts);
#endif
        register_outputs (opts);
        opts->process (argc, argv);
        delete opts;
    }

  // redirect printing (if necessary)
    if (printFile.compare("-") != 0) {
        wrld->_printTo = new std::ofstream (printFile);
        if (wrld->_printTo->fail()) {
            exit_with_error (wrld, "Error opening print file");
        }
        diderot::__details::config_ostream (*wrld->_printTo);
    }
    else {
        diderot::__details::config_ostream (std::cout);
    }

  // initialize scheduler stuff
    if (wrld->_verbose) {
        std::cerr << "initializing world ..." << std::endl;
    }
    if (wrld->init()) {
        exit_with_error (wrld, "Error initializing world");
    }

#ifndef DIDEROT_NO_INPUTS
  // initialize the input globals
    if (init_inputs (wrld, &inputs)) {
        exit_with_error (wrld, "Error initializing inputs");
    }
#endif

  // run the generated global initialization code
    if (wrld->_verbose) {
        std::cerr << "initializing globals and creating strands ...\n";
    }
    if (wrld->create_strands()) {
        exit_with_error (wrld, "Error in global initialization");
    }

#ifdef DIDEROT_EXEC_SNAPSHOT

    if (snapshotPeriod > 0) {
     // write initial state as snapshot 0
        if (write_snapshot (wrld, "-0000")) {
            exit_with_error (wrld, "Error generating snapshot");
        }
     // run the program for `snapshotPeriod` steps at a time with a snapshot after each run
        while (true) {
            uint32_t n, limit;
          // determine a step limit for the next run
            if (stepLimit > 0) {
                if (stepLimit <= nSteps) {
                    break;
                }
                limit = std::min(stepLimit - nSteps, snapshotPeriod);
            }
            else {
                limit = snapshotPeriod;
            }
          // run the program for upto limit steps
            if ((n = wrld->run (limit)) == 0) {
                break;
            }
            nSteps += n;
            if (wrld->_errors->errNum > 0) {
                break;
            }
            else if (wrld->_strands.num_alive() == 0) {
                wrld->error("no alive strands, so no snapshot at step %d", nSteps);
                break;
            }
          // write a snapshot with the step count as a suffix
            std::string suffix = std::to_string(nSteps);
            if (suffix.length() < 4) {
                suffix = std::string("0000").substr(0, 4 - suffix.length()) + suffix;
            }
            suffix = "-" + suffix;
            if (write_snapshot (wrld, suffix)) {
                exit_with_error (wrld, "Error generating snapshot");
	    }
        }
    }
    else {
        nSteps = wrld->run (stepLimit);
    }

#else // !DIDEROT_EXEC_SNAPSHOT

    nSteps = wrld->run (stepLimit);

#endif // DIDEROT_EXEC_SNAPSHOT

    if (wrld->_errors->errNum > 0) {
        exit_with_error (wrld, "Error during execution");
    }

    if ((stepLimit != 0) && (wrld->_strands.num_active() > 0)) {
#ifdef DIDEROT_STRAND_ARRAY
        if (wrld->_verbose) {
            std::cerr << "Step limit expired; "
                << wrld->_strands.num_active() << " active strands remaining" << std::endl;
        }
#else
      // step limit expired, so kill remaining strands
        if (wrld->_verbose) {
            std::cerr << "Step limit expired. Killing remaining "
                << wrld->_strands.num_active() << " active strands" << std::endl;
        }
        wrld->kill_all();
#endif
    }

    if (wrld->_verbose) {
        std::cerr << "done: " << nSteps << " steps, in " << wrld->_run_time << " seconds";
#ifndef DIDEROT_STRAND_ARRAY
        std::cerr << "; " << wrld->_strands.num_stable() << " stable strands" << std::endl;
#else
        std::cerr << std::endl;
#endif
    }
    else if (timingFlg) {
        std::cout << "usr=" << wrld->_run_time << std::endl;
    }

  // output the final strand states
    if (wrld->_strands.num_stable() > 0) {
        if (write_output (wrld)) {
            exit_with_error (wrld, "Error generating output");
        }
    }
    else {
        std::cerr << "Error: no stable strands at termination, so no output\n";
        delete wrld;
        return 1;
    }

    delete wrld;

    return 0;

} // main
/*---------- end seq-main.in ----------*/


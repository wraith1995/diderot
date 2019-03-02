/*---------- begin cxx-head.in ----------*/
/*! \file 15.cxx
 *
 * Generated from 15.diderot.
 *
 * Command: /home/teocollin/installs/diderot/diderot/bin/diderotc --log --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree 15.diderot
 * Version: master:2016-07-29
 */
/*---------- end cxx-head.in ----------*/

#define DIDEROT_STRAND_HAS_CONSTR
/*---------- begin lib-cxx-incl.in ----------*/
#include "15.h"
#include "diderot/diderot.hxx"

#ifdef DIDEROT_ENABLE_LOGGING
#define IF_LOGGING(...)         __VA_ARGS__
#else
#define IF_LOGGING(...)
#endif

static std::string ProgramName = "15";
/*---------- end lib-cxx-incl.in ----------*/

// ***** Begin synthesized types *****

namespace Diderot {
    typedef float vec15 __attribute__ ((vector_size (32)));
    typedef float vec7 __attribute__ ((vector_size (32)));
    struct tensor_ref_15 : public diderot::tensor_ref<float,15> {
        tensor_ref_15 (const float *src);
        tensor_ref_15 (struct tensor_15 const & ten);
        tensor_ref_15 (tensor_ref_15 const & ten);
    };
    struct tensor_15 : public diderot::tensor<float,15> {
        tensor_15 ()
            : diderot::tensor<float,15>()
        { }
        tensor_15 (std::initializer_list< float > const & il)
            : diderot::tensor<float,15>(il)
        { }
        tensor_15 (const float *src)
            : diderot::tensor<float,15>(src)
        { }
        tensor_15 (tensor_15 const & ten)
            : diderot::tensor<float,15>(ten._data)
        { }
        ~tensor_15 () { }
        tensor_15 & operator= (tensor_15 const & src);
        tensor_15 & operator= (tensor_ref_15 const & src);
        tensor_15 & operator= (std::initializer_list< float > const & il);
        tensor_15 & operator= (const float *src);
    };
    inline tensor_ref_15::tensor_ref_15 (const float *src)
        : diderot::tensor_ref<float,15>(src)
    { }
    inline tensor_ref_15::tensor_ref_15 (struct tensor_15 const & ten)
        : diderot::tensor_ref<float,15>(ten._data)
    { }
    inline tensor_ref_15::tensor_ref_15 (tensor_ref_15 const & ten)
        : diderot::tensor_ref<float,15>(ten._data)
    { }
    inline tensor_15 & tensor_15::operator= (tensor_15 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_15 & tensor_15::operator= (tensor_ref_15 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_15 & tensor_15::operator= (std::initializer_list< float > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_15 & tensor_15::operator= (const float *src)
    {
        this->copy(src);
        return *this;
    }
} // namespace Diderot
namespace diderot {
    template <>
    struct dynseq_traits<int32_t> {
        using value_type = int32_t;
        using base_type = int32_t;
        static const __details::load_fn_ptr<base_type> *load_fn_tbl;
        static const uint32_t values_per_elem = 1;
    };
    const __details::load_fn_ptr< dynseq_traits< int32_t >::base_type > *dynseq_traits< int32_t >::load_fn_tbl = nrrdILoad;
} // namespace diderot
// ***** End synthesized types *****

/*---------- begin namespace-open.in ----------*/
namespace Diderot {

static std::string ProgramName = "15";

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

typedef struct {
    bool gv_atest;
} defined_inputs;
struct globals {
    tensor_15 gv_atest;
    diderot::dynseq< int32_t > gv_a;
    ~globals () { }
};
struct gg_strand {
    tensor_15 sv_umm;
    float sv_temp;
};
/*---------- begin seq-sarr.in ----------*/
// forward declarations of strand methods
#ifdef DIDEROT_HAS_START_METHOD
static diderot::strand_status gg_start (gg_strand *self);
#endif // DIDEROT_HAS_START_METHOD
static diderot::strand_status gg_update (world *wrld, gg_strand *self);
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
    diderot::strand_status strand_update (world *wrld, index_t ix)
    {
        return gg_update(wrld, this->strand(ix));
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
    defined_inputs _definedInp;
    world ();
    ~world ();
    bool init ();
    bool alloc (int32_t base[1], uint32_t size[1]);
    bool create_strands ();
    uint32_t run (uint32_t max_nsteps);
    void swap_state ();
};
// ***** Begin synthesized operations *****

inline float vdot7 (vec7 u, vec7 v)
{
    vec7 w = u * v;
    return w[0] + w[1] + w[2] + w[3] + w[4] + w[5] + w[6];
}
inline float vdot8 (vec8 u, vec8 v)
{
    vec8 w = u * v;
    return w[0] + w[1] + w[2] + w[3] + w[4] + w[5] + w[6] + w[7];
}
inline vec15 vload15 (const float *vp)
{
    return __extension__ (vec15){vp[0], vp[1], vp[2], vp[3], vp[4], vp[5], vp[6], vp[7]};
}
inline vec7 vload7 (const float *vp)
{
    return __extension__ (vec7){vp[0], vp[1], vp[2], vp[3], vp[4], vp[5], vp[6], 0.e0f};
}
// ***** End synthesized operations *****

extern "C" void Diderot_input_get_atest (Diderot_world_t *cWrld, float v[15])
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    std::memcpy(v, wrld->_globals->gv_atest.addr(0), sizeof(tensor_15));
}
extern "C" bool Diderot_input_set_atest (Diderot_world_t *cWrld, float v[15])
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_atest = true;
    wrld->_globals->gv_atest = v;
    return false;
}
static bool check_defined (world *wrld)
{
    return false;
}
static void init_defined_inputs (world *wrld)
{
    wrld->_definedInp.gv_atest = false;
}
static void init_defaults (globals *glob)
{
    glob->gv_atest[0] = 0.e0f;
    glob->gv_atest[1] = 0.1e1f;
    glob->gv_atest[2] = 0.2e1f;
    glob->gv_atest[3] = 0.3e1f;
    glob->gv_atest[4] = 0.4e1f;
    glob->gv_atest[5] = 0.5e1f;
    glob->gv_atest[6] = 0.6e1f;
    glob->gv_atest[7] = 0.8e1f;
    glob->gv_atest[8] = 0.9e1f;
    glob->gv_atest[9] = 0.1e2f;
    glob->gv_atest[10] = 0.11e2f;
    glob->gv_atest[11] = 0.12e2f;
    glob->gv_atest[12] = 0.13e2f;
    glob->gv_atest[13] = 0.14e2f;
    glob->gv_atest[14] = 0.15e2f;
}
static bool init_globals (world *wrld)
{
    globals *glob = wrld->_globals;
    diderot::array< int32_t, 1 > t_2;
    t_2[0] = 1;
    diderot::dynseq< int32_t > l__t_1 = diderot::dynseq< int32_t >(1, t_2.data());
    glob->gv_a = l__t_1;
    return false;
}
static void gg_init (globals *glob, gg_strand *self, int32_t p_i_3)
{
    self->sv_umm = tensor_ref_15(glob->gv_atest);
    self->sv_temp = 0.e0f;
}
static diderot::strand_status gg_update (world *wrld, gg_strand *self)
{
    float l_temp_4 = vdot7(vload7(tensor_ref_15(self->sv_umm).addr(8)), vload7(tensor_ref_15(self->sv_umm).addr(8))) + vdot8(
        vload15(tensor_ref_15(self->sv_umm).addr(0)), vload15(tensor_ref_15(self->sv_umm).addr(0)));
    wrld->print() << l_temp_4 << std::flush;
    self->sv_temp = l_temp_4;
    return diderot::kStabilize;
}
extern "C" bool Diderot_output_get_umm (Diderot_world_t *cWrld, Nrrd *nData)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    // Compute sizes of nrrd file
    size_t sizes[2];
    sizes[0] = 15;
    sizes[1] = wrld->_strands.num_stable();
    // Allocate nData nrrd
    if (nrrdMaybeAlloc_nva(nData, nrrdTypeFloat, 2, sizes) != 0) {
        char *msg = biffGetDone(NRRD);
        biffMsgAdd(wrld->_errors, msg);
        std::free(msg);
        return true;
    }
    // copy data to output nrrd
    char *cp = reinterpret_cast<char *>(nData->data);
    for (auto ix = wrld->_strands.begin_stable(); ix != wrld->_strands.end_stable(); ix = wrld->_strands.next_stable(
        ix)) {
        memcpy(cp, &wrld->_strands.strand(ix)->sv_umm, 15 * sizeof(float));
        cp += 15 * sizeof(float);
    }
    nData->axis[0].kind = nrrdKindList;
    nData->axis[1].kind = nrrdKindList;
    return false;
}
extern "C" bool Diderot_output_get_temp (Diderot_world_t *cWrld, Nrrd *nData)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
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
        memcpy(cp, &wrld->_strands.strand(ix)->sv_temp, 1 * sizeof(float));
        cp += 1 * sizeof(float);
    }
    nData->axis[0].kind = nrrdKindList;
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
    diderot::dynseq< int32_t > seq_0 = glob->gv_a;
    int32_t base[1] = {0,};
    uint32_t size[1] = {static_cast<uint32_t>(seq_0.length()),};
    if (this->alloc(base, size)) {
        return true;
    }
    uint32_t ix = 0;
    for (auto it_1 = seq_0.cbegin(); it_1 != seq_0.cend(); ++it_1) {
        auto i_j_5 = *it_1;
        gg_init(this->_globals, this->_strands.strand(ix), i_j_5);
        ++ix;
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
            sts = this->_strands.strand_update(this, ix);
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

/*---------- begin c-wrappers.in ----------*/
extern "C" uint32_t Diderot_num_strands (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    return w->_strands.num_alive();
}

extern "C" uint32_t Diderot_num_active_strands (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    return w->_strands.num_active();
}

extern "C" uint32_t Diderot_num_stable_strands (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    return w->_strands.num_stable();
}

extern "C" bool Diderot_any_errors (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    return (w->_errors->errNum > 0);
}

extern "C" char *Diderot_get_errors (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    char *msg = biffMsgStrGet (w->_errors);
    biffMsgClear (w->_errors);
    return msg;
}

extern "C" Diderot_world_t *Diderot_new_world ()
{
    Diderot::world *w = new (std::nothrow) Diderot::world();
    return reinterpret_cast<Diderot_world_t *>(w);
}

extern "C" bool Diderot_init_world (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);

    if (w->_stage != diderot::POST_NEW) {
        w->error ("multiple calls to Diderot_init_world");
        return true;
    }

    if (w->init()) {
        return true;
    }

#ifndef DIDEROT_NO_INPUTS
    if (w != nullptr) {
        init_defined_inputs (w);
        init_defaults (w->_globals);
    }
#endif

    return false;
}

extern "C" bool Diderot_create_strands (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);

    if (w->_stage < diderot::POST_INIT) {
        w->error ("must call Diderot_init_world before Diderot_create_strands");
        return true;
    }
    else if (w->_stage > diderot::POST_INIT) {
        w->error ("multiple calls to Diderot_create_strands");
        return true;
    }

#ifdef DIDEROT_TARGET_PARALLEL
    if (w->_sched->create_workers (w)) {
        return true;
    }
#endif

#ifndef DIDEROT_NO_INPUTS
    if (check_defined(w)) {
        return true;
    }
#endif

    return static_cast<bool>(w->create_strands());
}

extern "C" uint32_t Diderot_run (Diderot_world_t *wrld, uint32_t maxNSteps)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);

    if (w->_stage < diderot::POST_CREATE) {
        w->error ("attempt to run uninitialized program");
        return 0;
    }
    else if (w->_stage == diderot::DONE) {
        return 0;
    }

    return w->run(maxNSteps);
}

extern "C" void Diderot_shutdown (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    delete w;
}

extern "C" void Diderot_set_verbose (Diderot_world_t *wrld, bool mode)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    w->_verbose = (mode ? true : false);
}

extern "C" bool Diderot_get_verbose (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    return static_cast<bool>(w->_verbose);
}

extern "C" bool Diderot_set_printer_cb (Diderot_world_t *wrld, bool (*pr)(void *, char *), void *data)
{
  /* FIXME: implement printer callback */
    return true;
}

#ifdef DIDEROT_TARGET_PARALLEL

extern "C" uint32_t Diderot_get_num_cores (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    return w->_sched->_numHWCores;
}

extern "C" bool Diderot_set_num_workers (Diderot_world_t *wrld, uint32_t nw)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    if (w->_sched->_numHWCores < nw) {
        w->_sched->_numWorkers = w->_sched->_numHWCores;
        return true;
    }
    else if (nw > 0) {
        w->_sched->_numWorkers = nw;
    }
    else {
        w->_sched->_numWorkers = w->_sched->_numHWCores;
    }
    return false;
}

extern "C" uint32_t Diderot_get_num_workers (Diderot_world_t *wrld)
{
    Diderot::world *w = reinterpret_cast<Diderot::world *>(wrld);
    return w->_sched->_numWorkers;
}

#endif /* DIDEROT_TARGET_PARALLEL */
/*---------- end c-wrappers.in ----------*/


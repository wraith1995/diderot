/*---------- begin cxx-head.in ----------*/
/*! \file qd.cxx
 *
 * Generated from qd.diderot.
 *
 * Command: ../bin/diderotc --exec --log --dump-all qd.diderot
 * Version: master:2016-07-29
 */
/*---------- end cxx-head.in ----------*/

#define DIDEROT_STRAND_HAS_CONSTR
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

namespace Diderot {
    typedef float vec2 __attribute__ ((vector_size (8)));
    typedef float vec3 __attribute__ ((vector_size (16)));
    struct tensor_ref_2 : public diderot::tensor_ref<float,2> {
        tensor_ref_2 (const float *src);
        tensor_ref_2 (struct tensor_2 const & ten);
        tensor_ref_2 (tensor_ref_2 const & ten);
    };
    struct tensor_ref_3 : public diderot::tensor_ref<float,3> {
        tensor_ref_3 (const float *src);
        tensor_ref_3 (struct tensor_3 const & ten);
        tensor_ref_3 (tensor_ref_3 const & ten);
    };
    struct tensor_ref_2_2 : public diderot::tensor_ref<float,4> {
        tensor_ref_2_2 (const float *src);
        tensor_ref_2_2 (struct tensor_2_2 const & ten);
        tensor_ref_2_2 (tensor_ref_2_2 const & ten);
        tensor_ref_2 last (uint32_t i)
        {
            return &this->_data[i];
        }
    };
    struct tensor_ref_3_3 : public diderot::tensor_ref<float,9> {
        tensor_ref_3_3 (const float *src);
        tensor_ref_3_3 (struct tensor_3_3 const & ten);
        tensor_ref_3_3 (tensor_ref_3_3 const & ten);
        tensor_ref_3 last (uint32_t i)
        {
            return &this->_data[i];
        }
    };
    struct tensor_3 : public diderot::tensor<float,3> {
        tensor_3 ()
            : diderot::tensor<float,3>()
        { }
        tensor_3 (std::initializer_list< float > const & il)
            : diderot::tensor<float,3>(il)
        { }
        tensor_3 (const float *src)
            : diderot::tensor<float,3>(src)
        { }
        tensor_3 (tensor_3 const & ten)
            : diderot::tensor<float,3>(ten._data)
        { }
        ~tensor_3 () { }
        tensor_3 & operator= (tensor_3 const & src);
        tensor_3 & operator= (tensor_ref_3 const & src);
        tensor_3 & operator= (std::initializer_list< float > const & il);
        tensor_3 & operator= (const float *src);
    };
    struct tensor_2 : public diderot::tensor<float,2> {
        tensor_2 ()
            : diderot::tensor<float,2>()
        { }
        tensor_2 (std::initializer_list< float > const & il)
            : diderot::tensor<float,2>(il)
        { }
        tensor_2 (const float *src)
            : diderot::tensor<float,2>(src)
        { }
        tensor_2 (tensor_2 const & ten)
            : diderot::tensor<float,2>(ten._data)
        { }
        ~tensor_2 () { }
        tensor_2 & operator= (tensor_2 const & src);
        tensor_2 & operator= (tensor_ref_2 const & src);
        tensor_2 & operator= (std::initializer_list< float > const & il);
        tensor_2 & operator= (const float *src);
    };
    struct tensor_3_3 : public diderot::tensor<float,9> {
        tensor_3_3 ()
            : diderot::tensor<float,9>()
        { }
        tensor_3_3 (std::initializer_list< float > const & il)
            : diderot::tensor<float,9>(il)
        { }
        tensor_3_3 (const float *src)
            : diderot::tensor<float,9>(src)
        { }
        tensor_3_3 (tensor_3_3 const & ten)
            : diderot::tensor<float,9>(ten._data)
        { }
        ~tensor_3_3 () { }
        tensor_3_3 & operator= (tensor_3_3 const & src);
        tensor_3_3 & operator= (tensor_ref_3_3 const & src);
        tensor_3_3 & operator= (std::initializer_list< float > const & il);
        tensor_3_3 & operator= (const float *src);
        tensor_ref_3 last (uint32_t i)
        {
            return &this->_data[i];
        }
    };
    struct tensor_2_2 : public diderot::tensor<float,4> {
        tensor_2_2 ()
            : diderot::tensor<float,4>()
        { }
        tensor_2_2 (std::initializer_list< float > const & il)
            : diderot::tensor<float,4>(il)
        { }
        tensor_2_2 (const float *src)
            : diderot::tensor<float,4>(src)
        { }
        tensor_2_2 (tensor_2_2 const & ten)
            : diderot::tensor<float,4>(ten._data)
        { }
        ~tensor_2_2 () { }
        tensor_2_2 & operator= (tensor_2_2 const & src);
        tensor_2_2 & operator= (tensor_ref_2_2 const & src);
        tensor_2_2 & operator= (std::initializer_list< float > const & il);
        tensor_2_2 & operator= (const float *src);
        tensor_ref_2 last (uint32_t i)
        {
            return &this->_data[i];
        }
    };
    inline tensor_ref_2::tensor_ref_2 (const float *src)
        : diderot::tensor_ref<float,2>(src)
    { }
    inline tensor_ref_2::tensor_ref_2 (struct tensor_2 const & ten)
        : diderot::tensor_ref<float,2>(ten._data)
    { }
    inline tensor_ref_2::tensor_ref_2 (tensor_ref_2 const & ten)
        : diderot::tensor_ref<float,2>(ten._data)
    { }
    inline tensor_ref_3::tensor_ref_3 (const float *src)
        : diderot::tensor_ref<float,3>(src)
    { }
    inline tensor_ref_3::tensor_ref_3 (struct tensor_3 const & ten)
        : diderot::tensor_ref<float,3>(ten._data)
    { }
    inline tensor_ref_3::tensor_ref_3 (tensor_ref_3 const & ten)
        : diderot::tensor_ref<float,3>(ten._data)
    { }
    inline tensor_ref_2_2::tensor_ref_2_2 (const float *src)
        : diderot::tensor_ref<float,4>(src)
    { }
    inline tensor_ref_2_2::tensor_ref_2_2 (struct tensor_2_2 const & ten)
        : diderot::tensor_ref<float,4>(ten._data)
    { }
    inline tensor_ref_2_2::tensor_ref_2_2 (tensor_ref_2_2 const & ten)
        : diderot::tensor_ref<float,4>(ten._data)
    { }
    inline tensor_ref_3_3::tensor_ref_3_3 (const float *src)
        : diderot::tensor_ref<float,9>(src)
    { }
    inline tensor_ref_3_3::tensor_ref_3_3 (struct tensor_3_3 const & ten)
        : diderot::tensor_ref<float,9>(ten._data)
    { }
    inline tensor_ref_3_3::tensor_ref_3_3 (tensor_ref_3_3 const & ten)
        : diderot::tensor_ref<float,9>(ten._data)
    { }
    inline tensor_3 & tensor_3::operator= (tensor_3 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_3 & tensor_3::operator= (tensor_ref_3 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_3 & tensor_3::operator= (std::initializer_list< float > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_3 & tensor_3::operator= (const float *src)
    {
        this->copy(src);
        return *this;
    }
    inline tensor_2 & tensor_2::operator= (tensor_2 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_2 & tensor_2::operator= (tensor_ref_2 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_2 & tensor_2::operator= (std::initializer_list< float > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_2 & tensor_2::operator= (const float *src)
    {
        this->copy(src);
        return *this;
    }
    inline tensor_3_3 & tensor_3_3::operator= (tensor_3_3 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_3_3 & tensor_3_3::operator= (tensor_ref_3_3 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_3_3 & tensor_3_3::operator= (std::initializer_list< float > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_3_3 & tensor_3_3::operator= (const float *src)
    {
        this->copy(src);
        return *this;
    }
    inline tensor_2_2 & tensor_2_2::operator= (tensor_2_2 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_2_2 & tensor_2_2::operator= (tensor_ref_2_2 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_2_2 & tensor_2_2::operator= (std::initializer_list< float > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_2_2 & tensor_2_2::operator= (const float *src)
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

static std::string ProgramName = "qd";

struct world;
struct dump_strand;
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
    tensor_2 gv_v1;
    tensor_2 gv_v2;
    tensor_3 gv_w1;
    tensor_3 gv_w2;
    tensor_3 gv_w3;
    diderot::dynseq< int32_t > gv_itter;
    ~globals () { }
};
struct dump_strand {
    tensor_2_2 sv_t1;
    tensor_3_3 sv_t2;
};
/*---------- begin seq-sarr.in ----------*/
// forward declarations of strand methods
#ifdef DIDEROT_HAS_START_METHOD
static diderot::strand_status dump_start (dump_strand *self);
#endif // DIDEROT_HAS_START_METHOD
static diderot::strand_status dump_update (world *wrld, dump_strand *self);
#ifdef DIDEROT_HAS_STABILIZE_METHOD
static void dump_stabilize (dump_strand *self);
#endif // DIDEROT_HAS_STABILIZE_METHOD

// if we have both communication and "die", then we need to track when strands die
// so that we can rebuild the list of strands use to construct the kd-tree
#if defined(DIDEROT_HAS_STRAND_COMMUNICATION) && !defined(DIDEROT_HAS_STRAND_DIE)
#  define TRACK_STRAND_DEATH
#endif

// strand_array for SEQUENTIAL/NO BSP/SINGLE STATE/DIRECT ACCESS
//
struct strand_array {
    typedef dump_strand strand_t;
    typedef uint32_t index_t;
    typedef index_t sid_t;              // strand ID (index into strand-state storage)

    uint8_t             *_status;       // the array of status information for the strands
    char                *_storage;      // points to array of dump_strand structs
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
    dump_strand *id_to_strand (sid_t id) const
    {
        assert (id < this->_nItems);
        return reinterpret_cast<dump_strand *>(this->_storage + id * sizeof(dump_strand));
    }

  // return a strand's status
    diderot::strand_status status (index_t ix) const
    {
        assert (ix < this->_nItems);
        return static_cast<diderot::strand_status>(this->_status[ix]);
    }
  // return a pointer to the given strand
    dump_strand *strand (index_t ix) const
    {
        return this->id_to_strand(this->id(ix));
    }
  // return a pointer to the local state of strand ix
    dump_strand *local_state (index_t ix) const
    {
        return this->strand(ix);
    }
  // return a pointer to the local state of strand with the given ID
    dump_strand *id_to_local_state (sid_t id) const
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
        this->_storage = static_cast<char *>(std::malloc (nItems * sizeof(dump_strand)));
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
            new (this->strand(ix)) dump_strand;
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
        return dump_start(this->strand(ix));
    }
#endif // DIDEROT_HAS_START_METHOD

  // invoke strand's update method
    diderot::strand_status strand_update (world *wrld, index_t ix)
    {
        return dump_update(wrld, this->strand(ix));
    }

  // invoke strand's stabilize method
    index_t strand_stabilize (index_t ix)
    {
#ifdef DIDEROT_HAS_STABILIZE_METHOD
        dump_stabilize (this->strand(ix));
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
        this->strand(ix)->~dump_strand();
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

inline vec2 vload2 (const float *vp)
{
    return __extension__ (vec2){vp[0], vp[1]};
}
inline vec2 vcons2 (float r0, float r1)
{
    return __extension__ (vec2){r0, r1};
}
inline vec3 vload3 (const float *vp)
{
    return __extension__ (vec3){vp[0], vp[1], vp[2], 0.e0f};
}
static std::ostream& operator<< (std::ostream& outs, tensor_ref_2_2 const & ten)
{
    return outs << "[[" << ten._data[0] << "," << ten._data[1] << "],[" << ten._data[2] << "," << ten._data[3] << "]]";
}
static std::ostream& operator<< (std::ostream& outs, tensor_ref_3_3 const & ten)
{
    return outs << "[[" << ten._data[0] << "," << ten._data[1] << "," << ten._data[2] << "],[" << ten._data[3] << "," << ten._data[4] << "," << ten._data[5] << "],[" << ten._data[6] << "," << ten._data[7] << "," << ten._data[8] << "]]";
}
inline vec3 vcons3 (float r0, float r1, float r2)
{
    return __extension__ (vec3){r0, r1, r2, 0.e0f};
}
inline float vdot2 (vec2 u, vec2 v)
{
    vec2 w = u * v;
    return w[0] + w[1];
}
inline float vdot3 (vec3 u, vec3 v)
{
    vec3 w = u * v;
    return w[0] + w[1] + w[2];
}
// ***** End synthesized operations *****

typedef struct {
    tensor_2 gv_v1;
    tensor_2 gv_v2;
    tensor_3 gv_w1;
    tensor_3 gv_w2;
    tensor_3 gv_w3;
} cmd_line_inputs;
static void init_defaults (cmd_line_inputs *inp)
{
    inp->gv_v1[0] = 0.1e1f;
    inp->gv_v1[1] = 0.e0f;
    inp->gv_v2[0] = 0.2e1f;
    inp->gv_v2[1] = 0.1e1f;
    inp->gv_w1[0] = 0.1e1f;
    inp->gv_w1[1] = 0.e0f;
    inp->gv_w1[2] = 0.e0f;
    inp->gv_w2[0] = 0.1e1f;
    inp->gv_w2[1] = 0.1e1f;
    inp->gv_w2[2] = 0.e0f;
    inp->gv_w3[0] = 0.1e1f;
    inp->gv_w3[1] = 0.1e1f;
    inp->gv_w3[2] = 0.1e1f;
}
static void register_inputs (cmd_line_inputs *inp, diderot::options *opts)
{
    opts->add("v1", "", 2, inp->gv_v1._data, true);
    opts->add("v2", "", 2, inp->gv_v2._data, true);
    opts->add("w1", "", 3, inp->gv_w1._data, true);
    opts->add("w2", "", 3, inp->gv_w2._data, true);
    opts->add("w3", "", 3, inp->gv_w3._data, true);
}
static bool init_inputs (world *wrld, cmd_line_inputs *inp)
{
    globals *glob = wrld->_globals;
    glob->gv_v1 = inp->gv_v1;
    glob->gv_v2 = inp->gv_v2;
    glob->gv_w1 = inp->gv_w1;
    glob->gv_w2 = inp->gv_w2;
    glob->gv_w3 = inp->gv_w3;
    return false;
}
static std::string OutPrefix_t1 = "t1.nrrd";
static std::string OutPrefix_t2 = "t2.nrrd";
static void register_outputs (diderot::options *opts)
{
    opts->add("o-t1,output-t1", "specify output file for t1", &OutPrefix_t1, true);
    opts->add("o-t2,output-t2", "specify output file for t2", &OutPrefix_t2, true);
}
static bool init_globals (world *wrld)
{
    globals *glob = wrld->_globals;
    diderot::array< int32_t, 1 > t_6;
    t_6[0] = 1;
    diderot::dynseq< int32_t > l__t_5 = diderot::dynseq< int32_t >(1, t_6.data());
    glob->gv_itter = l__t_5;
    return false;
}
static void dump_init (globals *glob, dump_strand *self, int32_t p_i_7)
{
    self->sv_t1[0] = tensor_ref_2(glob->gv_v1)[0];
    self->sv_t1[1] = tensor_ref_2(glob->gv_v1)[1];
    self->sv_t1[2] = tensor_ref_2(glob->gv_v2)[0];
    self->sv_t1[3] = tensor_ref_2(glob->gv_v2)[1];
    self->sv_t2[0] = tensor_ref_3(glob->gv_w1)[0];
    self->sv_t2[1] = tensor_ref_3(glob->gv_w1)[1];
    self->sv_t2[2] = tensor_ref_3(glob->gv_w1)[2];
    self->sv_t2[3] = tensor_ref_3(glob->gv_w2)[0];
    self->sv_t2[4] = tensor_ref_3(glob->gv_w2)[1];
    self->sv_t2[5] = tensor_ref_3(glob->gv_w2)[2];
    self->sv_t2[6] = tensor_ref_3(glob->gv_w3)[0];
    self->sv_t2[7] = tensor_ref_3(glob->gv_w3)[1];
    self->sv_t2[8] = tensor_ref_3(glob->gv_w3)[2];
}
static diderot::strand_status dump_update (world *wrld, dump_strand *self)
{
    wrld->print() << tensor_ref_2_2(self->sv_t1) << "\n" << std::flush;
    float l_r_8 = tensor_ref_2_2(self->sv_t1)[0];
    float l_r_9 = tensor_ref_2_2(self->sv_t1)[1];
    float l_r_10 = 0.e0f * l_r_8 + -0.1e1f * l_r_9;
    float l_r_11 = tensor_ref_2_2(self->sv_t1)[2];
    float l_r_12 = 0.e0f * l_r_11;
    float l_r_13 = tensor_ref_2_2(self->sv_t1)[3];
    float l_r_14 = l_r_12 + -0.1e1f * l_r_13;
    float l_r_15 = 0.1e1f * l_r_8 + 0.e0f * l_r_9;
    float l_r_16 = 0.e0f * l_r_13;
    float l_r_17 = 0.1e1f * l_r_11 + l_r_16;
    float l_op1_e3_l_5_18 = vdot2(vload2(tensor_ref_2_2(self->sv_t1).last(0).addr(0)),
        vcons2(l_r_12 + 0.1e1f * l_r_13, -0.1e1f * l_r_11 + l_r_16));
    tensor_2_2 t_19 = {(0.e0f * l_r_10 + -0.1e1f * l_r_14) / 0.1e1f / l_op1_e3_l_5_18,(0.e0f * l_r_15 + -0.1e1f * l_r_17) / 0.1e1f / l_op1_e3_l_5_18,(0.1e1f * l_r_10 + 0.e0f * l_r_14) / 0.1e1f / l_op1_e3_l_5_18,(0.1e1f * l_r_15 + 0.e0f * l_r_17) / 0.1e1f / l_op1_e3_l_5_18,};
    wrld->print() << tensor_ref_2_2(t_19) << "\n" << std::flush;
    wrld->print() << tensor_ref_3_3(self->sv_t2) << "\n" << std::flush;
    float l_r_20 = tensor_ref_3_3(self->sv_t2)[0];
    float l_r_21 = 0.e0f * l_r_20;
    float l_r_22 = tensor_ref_3_3(self->sv_t2)[1];
    float l_r_23 = 0.e0f * l_r_22;
    float l_r_24 = tensor_ref_3_3(self->sv_t2)[2];
    float l_r_25 = 0.e0f * l_r_24;
    float l_r_26 = l_r_21 + l_r_23;
    float l_r_27 = l_r_26 + l_r_25;
    float l_r_28 = tensor_ref_3_3(self->sv_t2)[3];
    float l_r_29 = 0.e0f * l_r_28;
    float l_r_30 = tensor_ref_3_3(self->sv_t2)[4];
    float l_r_31 = 0.e0f * l_r_30;
    float l_r_32 = tensor_ref_3_3(self->sv_t2)[5];
    float l_r_33 = 0.e0f * l_r_32;
    float l_r_34 = l_r_29 + l_r_31;
    float l_r_35 = l_r_34 + l_r_33;
    float l_r_36 = tensor_ref_3_3(self->sv_t2)[6];
    float l_r_37 = 0.e0f * l_r_36;
    float l_r_38 = tensor_ref_3_3(self->sv_t2)[7];
    float l_r_39 = 0.e0f * l_r_38;
    float l_r_40 = tensor_ref_3_3(self->sv_t2)[8];
    float l_r_41 = 0.e0f * l_r_40;
    float l_r_42 = l_r_37 + l_r_39;
    float l_r_43 = l_r_42 + l_r_41;
    float l_r_44 = l_r_26 + -0.1e1f * l_r_24;
    float l_r_45 = l_r_34 + -0.1e1f * l_r_32;
    float l_r_46 = l_r_42 + -0.1e1f * l_r_40;
    float l_r_47 = l_r_21 + 0.1e1f * l_r_22 + l_r_25;
    float l_r_48 = l_r_29 + 0.1e1f * l_r_30 + l_r_33;
    float l_r_49 = l_r_37 + 0.1e1f * l_r_38 + l_r_41;
    float l_r_50 = l_r_26 + 0.1e1f * l_r_24;
    float l_r_51 = l_r_34 + 0.1e1f * l_r_32;
    float l_r_52 = l_r_42 + 0.1e1f * l_r_40;
    float l_r_53 = -0.1e1f * l_r_20 + l_r_23 + l_r_25;
    float l_r_54 = -0.1e1f * l_r_28 + l_r_31 + l_r_33;
    float l_r_55 = -0.1e1f * l_r_36 + l_r_39 + l_r_41;
    float l_r_56 = l_r_21 + -0.1e1f * l_r_22 + l_r_25;
    float l_r_57 = l_r_29 + -0.1e1f * l_r_30 + l_r_33;
    float l_r_58 = l_r_37 + -0.1e1f * l_r_38 + l_r_41;
    float l_r_59 = 0.1e1f * l_r_20 + l_r_23 + l_r_25;
    float l_r_60 = 0.1e1f * l_r_28 + l_r_31 + l_r_33;
    float l_r_61 = 0.1e1f * l_r_36 + l_r_39 + l_r_41;
    float l_r_62 = l_r_20 * l_r_35 + l_r_22 * l_r_51 + l_r_24 * l_r_57;
    float l_r_63 = l_r_20 * l_r_43 + l_r_22 * l_r_52 + l_r_24 * l_r_58;
    float l_r_64 = l_r_20 * l_r_45 + l_r_22 * l_r_35 + l_r_24 * l_r_60;
    float l_r_65 = l_r_20 * l_r_46 + l_r_22 * l_r_43 + l_r_24 * l_r_61;
    float l_r_66 = l_r_20 * l_r_48 + l_r_22 * l_r_54 + l_r_24 * l_r_35;
    float l_r_67 = l_r_20 * l_r_49 + l_r_22 * l_r_55 + l_r_24 * l_r_43;
    float l_r_68 = l_r_28 * l_r_27 + l_r_30 * l_r_50 + l_r_32 * l_r_56;
    float l_r_69 = l_r_28 * l_r_43 + l_r_30 * l_r_52 + l_r_32 * l_r_58;
    float l_r_70 = l_r_28 * l_r_44 + l_r_30 * l_r_27 + l_r_32 * l_r_59;
    float l_r_71 = l_r_28 * l_r_46 + l_r_30 * l_r_43 + l_r_32 * l_r_61;
    float l_r_72 = l_r_28 * l_r_47 + l_r_30 * l_r_53 + l_r_32 * l_r_27;
    float l_r_73 = l_r_28 * l_r_49 + l_r_30 * l_r_55 + l_r_32 * l_r_43;
    float l_r_74 = l_r_36 * l_r_27 + l_r_38 * l_r_50 + l_r_40 * l_r_56;
    float l_r_75 = l_r_36 * l_r_35 + l_r_38 * l_r_51 + l_r_40 * l_r_57;
    float l_r_76 = l_r_36 * l_r_44 + l_r_38 * l_r_27 + l_r_40 * l_r_59;
    float l_r_77 = l_r_36 * l_r_45 + l_r_38 * l_r_35 + l_r_40 * l_r_60;
    float l_r_78 = l_r_36 * l_r_47 + l_r_38 * l_r_53 + l_r_40 * l_r_27;
    float l_r_79 = l_r_36 * l_r_48 + l_r_38 * l_r_54 + l_r_40 * l_r_35;
    float l_r_80 = 0.e0f * (l_r_20 * l_r_27 + l_r_22 * l_r_50 + l_r_24 * l_r_56);
    float l_r_81 = 0.e0f * l_r_63;
    float l_r_82 = 0.e0f * l_r_68;
    float l_r_83 = 0.e0f * (l_r_28 * l_r_35 + l_r_30 * l_r_51 + l_r_32 * l_r_57);
    float l_r_84 = 0.e0f * l_r_74;
    float l_r_85 = 0.e0f * (l_r_36 * l_r_43 + l_r_38 * l_r_52 + l_r_40 * l_r_58);
    float l_r_86 = l_r_80 + 0.e0f * l_r_62;
    float l_r_87 = 0.e0f * (l_r_20 * l_r_44 + l_r_22 * l_r_27 + l_r_24 * l_r_59);
    float l_r_88 = 0.e0f * l_r_65;
    float l_r_89 = 0.e0f * l_r_70;
    float l_r_90 = 0.e0f * (l_r_28 * l_r_45 + l_r_30 * l_r_35 + l_r_32 * l_r_60);
    float l_r_91 = 0.e0f * l_r_76;
    float l_r_92 = 0.e0f * (l_r_36 * l_r_46 + l_r_38 * l_r_43 + l_r_40 * l_r_61);
    float l_r_93 = l_r_87 + 0.e0f * l_r_64;
    float l_r_94 = 0.e0f * (l_r_20 * l_r_47 + l_r_22 * l_r_53 + l_r_24 * l_r_27);
    float l_r_95 = 0.e0f * l_r_67;
    float l_r_96 = 0.e0f * l_r_72;
    float l_r_97 = 0.e0f * (l_r_28 * l_r_48 + l_r_30 * l_r_54 + l_r_32 * l_r_35);
    float l_r_98 = 0.e0f * l_r_78;
    float l_r_99 = 0.e0f * (l_r_36 * l_r_49 + l_r_38 * l_r_55 + l_r_40 * l_r_43);
    float l_r_100 = l_r_94 + 0.e0f * l_r_66;
    float l_r_101 = 0.e0f * l_r_69;
    float l_r_102 = 0.e0f * l_r_75;
    float l_r_103 = 0.e0f * l_r_71;
    float l_r_104 = 0.e0f * l_r_77;
    float l_r_105 = 0.e0f * l_r_73;
    float l_r_106 = 0.e0f * l_r_79;
    tensor_ref_3 l_projParam_107 = tensor_ref_3_3(self->sv_t2).last(3);
    float l_op1_e3_l_7_108 = vdot3(vload3(tensor_ref_3_3(self->sv_t2).last(0).addr(0)),
        vcons3(vdot3(vload3(l_projParam_107.addr(0)), vcons3(l_r_43, l_r_52, l_r_58)),
            vdot3(vload3(l_projParam_107.addr(0)), vcons3(l_r_46, l_r_43, l_r_61)),
            vdot3(vload3(l_projParam_107.addr(0)), vcons3(l_r_49, l_r_55, l_r_43))));
    tensor_3_3 t_109 = {(l_r_86 + l_r_81 + l_r_82 + l_r_83 + 0.1e1f * l_r_69 + l_r_84 + -0.1e1f * l_r_75 + l_r_85) / 0.2e1f / l_op1_e3_l_7_108,(l_r_93 + l_r_88 + l_r_89 + l_r_90 + 0.1e1f * l_r_71 + l_r_91 + -0.1e1f * l_r_77 + l_r_92) / 0.2e1f / l_op1_e3_l_7_108,(l_r_100 + l_r_95 + l_r_96 + l_r_97 + 0.1e1f * l_r_73 + l_r_98 + -0.1e1f * l_r_79 + l_r_99) / 0.2e1f / l_op1_e3_l_7_108,(l_r_86 + -0.1e1f * l_r_63 + l_r_82 + l_r_83 + l_r_101 + 0.1e1f * l_r_74 + l_r_102 + l_r_85) / 0.2e1f / l_op1_e3_l_7_108,(l_r_93 + -0.1e1f * l_r_65 + l_r_89 + l_r_90 + l_r_103 + 0.1e1f * l_r_76 + l_r_104 + l_r_92) / 0.2e1f / l_op1_e3_l_7_108,(l_r_100 + -0.1e1f * l_r_67 + l_r_96 + l_r_97 + l_r_105 + 0.1e1f * l_r_78 + l_r_106 + l_r_99) / 0.2e1f / l_op1_e3_l_7_108,(l_r_80 + 0.1e1f * l_r_62 + l_r_81 + -0.1e1f * l_r_68 + l_r_83 + l_r_101 + l_r_84 + l_r_102 + l_r_85) / 0.2e1f / l_op1_e3_l_7_108,(l_r_87 + 0.1e1f * l_r_64 + l_r_88 + -0.1e1f * l_r_70 + l_r_90 + l_r_103 + l_r_91 + l_r_104 + l_r_92) / 0.2e1f / l_op1_e3_l_7_108,(l_r_94 + 0.1e1f * l_r_66 + l_r_95 + -0.1e1f * l_r_72 + l_r_97 + l_r_105 + l_r_98 + l_r_106 + l_r_99) / 0.2e1f / l_op1_e3_l_7_108,};
    wrld->print() << tensor_ref_3_3(t_109) << "\n" << std::flush;
    return diderot::kStabilize;
}
bool output_get_t1 (world *wrld, Nrrd *nData)
{
    // Compute sizes of nrrd file
    size_t sizes[2];
    sizes[0] = 4;
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
        memcpy(cp, &wrld->_strands.strand(ix)->sv_t1, 4 * sizeof(float));
        cp += 4 * sizeof(float);
    }
    nData->axis[0].kind = nrrdKind2DMatrix;
    nData->axis[1].kind = nrrdKindList;
    return false;
}
bool output_get_t2 (world *wrld, Nrrd *nData)
{
    // Compute sizes of nrrd file
    size_t sizes[2];
    sizes[0] = 9;
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
        memcpy(cp, &wrld->_strands.strand(ix)->sv_t2, 9 * sizeof(float));
        cp += 9 * sizeof(float);
    }
    nData->axis[0].kind = nrrdKind3DMatrix;
    nData->axis[1].kind = nrrdKindList;
    return false;
}
static bool write_output (world *wrld)
{
    Nrrd *nData;
    nData = nrrdNew();
    if (output_get_t1(wrld, nData)) {
        wrld->error("Error getting nrrd data for \'t1\'");
        return true;
    }
    else if (nrrd_save_helper(OutPrefix_t1, nData)) {
        return true;
    }
    nrrdNuke(nData);
    nData = nrrdNew();
    if (output_get_t2(wrld, nData)) {
        wrld->error("Error getting nrrd data for \'t2\'");
        return true;
    }
    else if (nrrd_save_helper(OutPrefix_t2, nData)) {
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
    diderot::dynseq< int32_t > seq_0 = glob->gv_itter;
    int32_t base[1] = {0,};
    uint32_t size[1] = {static_cast<uint32_t>(seq_0.length()),};
    if (this->alloc(base, size)) {
        return true;
    }
    uint32_t ix = 0;
    for (auto it_1 = seq_0.cbegin(); it_1 != seq_0.cend(); ++it_1) {
        auto i_i_110 = *it_1;
        dump_init(this->_globals, this->_strands.strand(ix), i_i_110);
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


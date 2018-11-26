/*---------- begin cxx-head.in ----------*/
/*! \file umm.cxx
 *
 * Generated from test/umm.diderot.
 *
 * Command: bin/diderotc --log --exec --dump-pt --dump-ast --dump-simple test/umm.diderot
 * Version: master:2016-07-29
 */
/*---------- end cxx-head.in ----------*/

#define DIDEROT_STRAND_HAS_CONSTR
#define DIDEROT_NO_INPUTS
#define DIDEROT_STRAND_ARRAY
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
    typedef float vec5 __attribute__ ((vector_size (32)));
    typedef float vec4 __attribute__ ((vector_size (16)));
    struct tensor_ref_4 : public diderot::tensor_ref<float,4> {
        tensor_ref_4 (const float *src);
        tensor_ref_4 (struct tensor_4 const & ten);
        tensor_ref_4 (tensor_ref_4 const & ten);
    };
    struct tensor_4 : public diderot::tensor<float,4> {
        tensor_4 ()
            : diderot::tensor<float,4>()
        { }
        tensor_4 (std::initializer_list< float > const & il)
            : diderot::tensor<float,4>(il)
        { }
        tensor_4 (const float *src)
            : diderot::tensor<float,4>(src)
        { }
        tensor_4 (tensor_4 const & ten)
            : diderot::tensor<float,4>(ten._data)
        { }
        ~tensor_4 () { }
        tensor_4 & operator= (tensor_4 const & src);
        tensor_4 & operator= (tensor_ref_4 const & src);
        tensor_4 & operator= (std::initializer_list< float > const & il);
        tensor_4 & operator= (const float *src);
    };
    inline tensor_ref_4::tensor_ref_4 (const float *src)
        : diderot::tensor_ref<float,4>(src)
    { }
    inline tensor_ref_4::tensor_ref_4 (struct tensor_4 const & ten)
        : diderot::tensor_ref<float,4>(ten._data)
    { }
    inline tensor_ref_4::tensor_ref_4 (tensor_ref_4 const & ten)
        : diderot::tensor_ref<float,4>(ten._data)
    { }
    inline tensor_4 & tensor_4::operator= (tensor_4 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_4 & tensor_4::operator= (tensor_ref_4 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_4 & tensor_4::operator= (std::initializer_list< float > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_4 & tensor_4::operator= (const float *src)
    {
        this->copy(src);
        return *this;
    }
} // namespace Diderot
// ***** End synthesized types *****

/*---------- begin namespace-open.in ----------*/
namespace Diderot {

static std::string ProgramName = "umm";

struct world;
struct f_strand;
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
    tensor_4 gv_q0;
    tensor_4 gv_q1;
};
struct f_strand {
    float sv_b;
    tensor_4 sv_z1;
    float sv_a;
    tensor_4 sv_frm;
    float sv_q;
};
/*---------- begin seq-sarr.in ----------*/
// forward declarations of strand methods
#ifdef DIDEROT_HAS_START_METHOD
static diderot::strand_status f_start (f_strand *self);
#endif // DIDEROT_HAS_START_METHOD
static diderot::strand_status f_update (world *wrld, globals *glob, f_strand *self);
#ifdef DIDEROT_HAS_STABILIZE_METHOD
static void f_stabilize (f_strand *self);
#endif // DIDEROT_HAS_STABILIZE_METHOD

// if we have both communication and "die", then we need to track when strands die
// so that we can rebuild the list of strands use to construct the kd-tree
#if defined(DIDEROT_HAS_STRAND_COMMUNICATION) && !defined(DIDEROT_HAS_STRAND_DIE)
#  define TRACK_STRAND_DEATH
#endif

// strand_array for SEQUENTIAL/NO BSP/SINGLE STATE/DIRECT ACCESS
//
struct strand_array {
    typedef f_strand strand_t;
    typedef uint32_t index_t;
    typedef index_t sid_t;              // strand ID (index into strand-state storage)

    uint8_t             *_status;       // the array of status information for the strands
    char                *_storage;      // points to array of f_strand structs
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
    f_strand *id_to_strand (sid_t id) const
    {
        assert (id < this->_nItems);
        return reinterpret_cast<f_strand *>(this->_storage + id * sizeof(f_strand));
    }

  // return a strand's status
    diderot::strand_status status (index_t ix) const
    {
        assert (ix < this->_nItems);
        return static_cast<diderot::strand_status>(this->_status[ix]);
    }
  // return a pointer to the given strand
    f_strand *strand (index_t ix) const
    {
        return this->id_to_strand(this->id(ix));
    }
  // return a pointer to the local state of strand ix
    f_strand *local_state (index_t ix) const
    {
        return this->strand(ix);
    }
  // return a pointer to the local state of strand with the given ID
    f_strand *id_to_local_state (sid_t id) const
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
        this->_storage = static_cast<char *>(std::malloc (nItems * sizeof(f_strand)));
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
            new (this->strand(ix)) f_strand;
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
        return f_start(this->strand(ix));
    }
#endif // DIDEROT_HAS_START_METHOD

  // invoke strand's update method
    diderot::strand_status strand_update (world *wrld, globals *glob, index_t ix)
    {
        return f_update(wrld, glob, this->strand(ix));
    }

  // invoke strand's stabilize method
    index_t strand_stabilize (index_t ix)
    {
#ifdef DIDEROT_HAS_STABILIZE_METHOD
        f_stabilize (this->strand(ix));
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
        this->strand(ix)->~f_strand();
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
};
// ***** Begin synthesized operations *****

inline float vdot5 (vec5 u, vec5 v)
{
    vec5 w = u * v;
    return w[0] + w[1] + w[2] + w[3] + w[4];
}
static std::ostream& operator<< (std::ostream& outs, tensor_ref_4 const & ten)
{
    return outs << "[" << ten._data[0] << "," << ten._data[1] << "," << ten._data[2] << "," << ten._data[3] << "]";
}
inline vec4 vload4 (const float *vp)
{
    return __extension__ (vec4){vp[0], vp[1], vp[2], vp[3]};
}
inline vec4 vcons4 (float r0, float r1, float r2, float r3)
{
    return __extension__ (vec4){r0, r1, r2, r3};
}
inline vec5 vcons5 (float r0, float r1, float r2, float r3, float r4)
{
    return __extension__ (vec5){r0, r1, r2, r3, r4, 0.e0f, 0.e0f, 0.e0f};
}
inline void vpack4 (tensor_4 &dst, vec4 v0)
{
    dst._data[0] = v0[0];
    dst._data[1] = v0[1];
    dst._data[2] = v0[2];
    dst._data[3] = v0[3];
}
inline vec4 vscale4 (float s, vec4 v)
{
    return __extension__ (vec4){s, s, s, s} * v;
}
inline float vdot4 (vec4 u, vec4 v)
{
    vec4 w = u * v;
    return w[0] + w[1] + w[2] + w[3];
}
// ***** End synthesized operations *****

static std::string OutPrefix_a = "a.nrrd";
static std::string OutPrefix_frm = "frm.nrrd";
static std::string OutPrefix_q = "q.nrrd";
static void register_outputs (diderot::options *opts)
{
    opts->add("o-a,output-a", "specify output file for a", &OutPrefix_a, true);
    opts->add("o-frm,output-frm", "specify output file for frm", &OutPrefix_frm, true);
    opts->add("o-q,output-q", "specify output file for q", &OutPrefix_q, true);
}
static bool init_globals (world *wrld)
{
    globals *glob = wrld->_globals;
    glob->gv_q0[0] = 0.4e1f;
    glob->gv_q0[1] = 0.3e1f;
    glob->gv_q0[2] = 0.2e1f;
    glob->gv_q0[3] = 0.1e1f;
    glob->gv_q1[0] = 0.1e1f;
    glob->gv_q1[1] = 0.2e1f;
    glob->gv_q1[2] = 0.3e1f;
    glob->gv_q1[3] = 0.1e1f;
    return false;
}
static void f_init (globals *glob, f_strand *self, int32_t p_i_2)
{
    float l__t_10;
    float l_b_14;
    float l_b_17;
    float l_b_19;
    float l_b_20;
    float l_z_24;
    vec4 v_3 = vcons4(0.e0f, 0.e0f, 0.e0f, 0.e0f);
    float l_op1_e3_l_12_4 = 0.1e1f / std::sqrt(
        vdot4(vload4(tensor_ref_4(glob->gv_q0).addr(0)), vload4(tensor_ref_4(glob->gv_q0).addr(0))));
    float l_op1_e3_l_13_5 = 0.1e1f / std::sqrt(
        vdot4(vload4(tensor_ref_4(glob->gv_q1).addr(0)), vload4(tensor_ref_4(glob->gv_q1).addr(0))));
    float l_result_6 = l_op1_e3_l_12_4 * l_op1_e3_l_13_5 * vdot4(vload4(tensor_ref_4(glob->gv_q0).addr(0)),
        vload4(tensor_ref_4(glob->gv_q1).addr(0)));
    vec4 v_7 = v_3;
    float l__t_8 = 0.1e1f - l_result_6 * l_result_6;
    if (l_result_6 < 0.e0f) {
        vec4 v_9 = v_7 - vscale4(l_op1_e3_l_12_4, vload4(tensor_ref_4(glob->gv_q0).addr(0))) - vscale4(l_op1_e3_l_13_5,
            vload4(tensor_ref_4(glob->gv_q1).addr(0)));
        l__t_10 = 0.314159265358979323846264338327950288e1f - 0.2e1f * std::asin(std::sqrt(vdot4(v_9, v_9)) / 0.2e1f);
    }
    else {
        vec4 v_11 = vscale4(l_op1_e3_l_12_4, vload4(tensor_ref_4(glob->gv_q0).addr(0))) - vscale4(l_op1_e3_l_13_5,
            vload4(tensor_ref_4(glob->gv_q1).addr(0)));
        l__t_10 = 0.2e1f * std::asin(std::sqrt(vdot4(v_11, v_11)) / 0.2e1f);
    }
    float l_op1_e3_l_7_12 = 0.1e1f - 0.5e0f;
    float l_r_13 = l__t_10 * l__t_10;
    if (l_r_13 * l_op1_e3_l_7_12 * l_op1_e3_l_7_12 + 0.1e1f == 0.1e1f) {
        l_b_14 = 0.1e1f;
    }
    else {
        float l_op1_e3_l_7_15 = l__t_10 * l_op1_e3_l_7_12;
        l_b_14 = std::sin(l_op1_e3_l_7_15) / l_op1_e3_l_7_15;
    }
    bool l__t_16 = l_r_13 + 0.1e1f == 0.1e1f;
    if (l__t_16) {
        l_b_17 = 0.1e1f;
    }
    else {
        l_b_17 = std::sin(l__t_10) / l__t_10;
    }
    float l_r_18 = 0.5e0f * l__t_10;
    if (l_r_18 * 0.5e0f * l__t_10 + 0.1e1f == 0.1e1f) {
        l_b_19 = 0.1e1f;
    }
    else {
        l_b_19 = std::sin(l_r_18) / l_r_18;
    }
    if (l__t_16) {
        l_b_20 = 0.1e1f;
    }
    else {
        l_b_20 = std::sin(l__t_10) / l__t_10;
    }
    float l_r_21 = l_op1_e3_l_7_12 * (l_b_14 / l_b_17) * l_op1_e3_l_12_4;
    float l_r_22 = 0.5e0f * (l_b_19 / l_b_20) * l_op1_e3_l_13_5;
    vec5 v_23 = vcons5(0.1e1f, 0.e0f, 0.e0f, 0.e0f, 0.e0f);
    l_z_24 = 0.2e1f;
    vec4 v_26 = vcons4(l_r_21 * tensor_ref_4(glob->gv_q0)[0], l_r_21 * tensor_ref_4(glob->gv_q0)[1],
        l_r_21 * tensor_ref_4(glob->gv_q0)[2], l_r_21 * tensor_ref_4(glob->gv_q0)[3]) + vcons4(
        l_r_22 * tensor_ref_4(glob->gv_q1)[0], l_r_22 * tensor_ref_4(glob->gv_q1)[1],
        l_r_22 * tensor_ref_4(glob->gv_q1)[2], l_r_22 * tensor_ref_4(glob->gv_q1)[3]);
    vec5 v_27 = v_23;
    for (int32_t i_j_25 = 0; i_j_25 <= 5; ++i_j_25) {
        l_z_24 = l_z_24 + std::sqrt(vdot5(v_27, v_27)) + static_cast<float>(i_j_25);
    }
    self->sv_b = 0.1e1f;
    vpack4(self->sv_z1, v_7);
    self->sv_a = l__t_8;
    vpack4(self->sv_frm, v_26);
    self->sv_q = l_z_24;
}
static diderot::strand_status f_update (world *wrld, globals *glob, f_strand *self)
{
    float l_z_30;
    vec4 v_31;
    l_z_30 = 0.2e1f;
    v_31 = vload4(tensor_ref_4(self->sv_z1).addr(0));
    for (int32_t i_j_32 = 0; i_j_32 <= 5; ++i_j_32) {
        l_z_30 = l_z_30 + std::abs(self->sv_b) + static_cast<float>(i_j_32);
        v_31 = v_31 + vload4(tensor_ref_4(glob->gv_q1).addr(0));
    }
    tensor_4 _arg_33;
    vpack4(_arg_33, v_31);
    wrld->print() << self->sv_q << " " << l_z_30 << " " << tensor_ref_4(_arg_33) << " " << tensor_ref_4(self->sv_frm) << "\n" << std::flush;
    vpack4(self->sv_z1, v_31);
    return diderot::kStabilize;
}
bool output_get_a (world *wrld, Nrrd *nData)
{
    // Compute sizes of nrrd file
    size_t sizes[1];
    sizes[0] = wrld->_size[0];
    // Allocate nData nrrd
    if (nrrdMaybeAlloc_nva(nData, nrrdTypeFloat, 1, sizes) != 0) {
        char *msg = biffGetDone(NRRD);
        biffMsgAdd(wrld->_errors, msg);
        std::free(msg);
        return true;
    }
    // copy data to output nrrd
    char *cp = reinterpret_cast<char *>(nData->data);
    for (auto ix = wrld->_strands.begin_alive(); ix != wrld->_strands.end_alive(); ix = wrld->_strands.next_alive(ix)) {
        memcpy(cp, &wrld->_strands.strand(ix)->sv_a, 1 * sizeof(float));
        cp += 1 * sizeof(float);
    }
    nData->axis[0].kind = nrrdKindSpace;
    return false;
}
bool output_get_frm (world *wrld, Nrrd *nData)
{
    // Compute sizes of nrrd file
    size_t sizes[2];
    sizes[0] = 4;
    sizes[1] = wrld->_size[0];
    // Allocate nData nrrd
    if (nrrdMaybeAlloc_nva(nData, nrrdTypeFloat, 2, sizes) != 0) {
        char *msg = biffGetDone(NRRD);
        biffMsgAdd(wrld->_errors, msg);
        std::free(msg);
        return true;
    }
    // copy data to output nrrd
    char *cp = reinterpret_cast<char *>(nData->data);
    for (auto ix = wrld->_strands.begin_alive(); ix != wrld->_strands.end_alive(); ix = wrld->_strands.next_alive(ix)) {
        memcpy(cp, &wrld->_strands.strand(ix)->sv_frm, 4 * sizeof(float));
        cp += 4 * sizeof(float);
    }
    nData->axis[0].kind = nrrdKind4Vector;
    nData->axis[1].kind = nrrdKindSpace;
    return false;
}
bool output_get_q (world *wrld, Nrrd *nData)
{
    // Compute sizes of nrrd file
    size_t sizes[1];
    sizes[0] = wrld->_size[0];
    // Allocate nData nrrd
    if (nrrdMaybeAlloc_nva(nData, nrrdTypeFloat, 1, sizes) != 0) {
        char *msg = biffGetDone(NRRD);
        biffMsgAdd(wrld->_errors, msg);
        std::free(msg);
        return true;
    }
    // copy data to output nrrd
    char *cp = reinterpret_cast<char *>(nData->data);
    for (auto ix = wrld->_strands.begin_alive(); ix != wrld->_strands.end_alive(); ix = wrld->_strands.next_alive(ix)) {
        memcpy(cp, &wrld->_strands.strand(ix)->sv_q, 1 * sizeof(float));
        cp += 1 * sizeof(float);
    }
    nData->axis[0].kind = nrrdKindSpace;
    return false;
}
static bool write_output (world *wrld)
{
    Nrrd *nData;
    nData = nrrdNew();
    if (output_get_a(wrld, nData)) {
        wrld->error("Error getting nrrd data for \'a\'");
        return true;
    }
    else if (nrrd_save_helper(OutPrefix_a, nData)) {
        return true;
    }
    nrrdNuke(nData);
    nData = nrrdNew();
    if (output_get_frm(wrld, nData)) {
        wrld->error("Error getting nrrd data for \'frm\'");
        return true;
    }
    else if (nrrd_save_helper(OutPrefix_frm, nData)) {
        return true;
    }
    nrrdNuke(nData);
    nData = nrrdNew();
    if (output_get_q(wrld, nData)) {
        wrld->error("Error getting nrrd data for \'q\'");
        return true;
    }
    else if (nrrd_save_helper(OutPrefix_q, nData)) {
        return true;
    }
    nrrdNuke(nData);
    return false;
}
/*---------- begin world-methods.in ----------*/
// Allocate the program's world
//
world::world ()
    : diderot::world_base (ProgramName, true, 1)
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
    int lo_0 = 0;
    int hi_1 = 1;
    int32_t base[1] = {lo_0,};
    uint32_t size[1] = {static_cast<uint32_t>(hi_1 - lo_0 + 1),};
    if (this->alloc(base, size)) {
        return true;
    }
    uint32_t ix = 0;
    for (int i_i_35 = lo_0; i_i_35 <= hi_1; i_i_35++) {
        f_init(this->_globals, this->_strands.strand(ix), i_i_35);
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


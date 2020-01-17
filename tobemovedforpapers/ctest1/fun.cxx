/*---------- begin cxx-head.in ----------*/
/*! \file fun.cxx
 *
 * Generated from fun.diderot.
 *
 * Command: ../bin/diderotc --log --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree fun.diderot
 * Version: master:2016-07-29
 */
/*---------- end cxx-head.in ----------*/

#define DIDEROT_STRAND_HAS_CONSTR
#define DIDEROT_STRAND_ARRAY
/*---------- begin lib-cxx-incl.in ----------*/
#include "fun.h"
#include "diderot/diderot.hxx"

#ifdef DIDEROT_ENABLE_LOGGING
#define IF_LOGGING(...)         __VA_ARGS__
#else
#define IF_LOGGING(...)
#endif

static std::string ProgramName = "fun";
/*---------- end lib-cxx-incl.in ----------*/

// ***** Begin synthesized types *****

namespace Diderot {
    typedef float vec3 __attribute__ ((vector_size (16)));
    typedef float vec6 __attribute__ ((vector_size (32)));
    struct tensor_ref_3 : public diderot::tensor_ref<float,3> {
        tensor_ref_3 (const float *src);
        tensor_ref_3 (struct tensor_3 const & ten);
        tensor_ref_3 (tensor_ref_3 const & ten);
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
    inline tensor_ref_3::tensor_ref_3 (const float *src)
        : diderot::tensor_ref<float,3>(src)
    { }
    inline tensor_ref_3::tensor_ref_3 (struct tensor_3 const & ten)
        : diderot::tensor_ref<float,3>(ten._data)
    { }
    inline tensor_ref_3::tensor_ref_3 (tensor_ref_3 const & ten)
        : diderot::tensor_ref<float,3>(ten._data)
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
} // namespace Diderot
// ***** End synthesized types *****

/*---------- begin namespace-open.in ----------*/
namespace Diderot {

static std::string ProgramName = "fun";

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

typedef struct {
    bool gv_A;
} defined_inputs;
struct globals {
    diderot::image3d< float, float, 3 > gv_A;
    int32_t gv_shape;
    ~globals ()
    {
        this->gv_A.unregister_global();
    }
};
struct f_strand {
    tensor_3 sv_val;
    tensor_3_3 sv_result;
};
/*---------- begin seq-sarr.in ----------*/
// forward declarations of strand methods
#ifdef DIDEROT_HAS_START_METHOD
static diderot::strand_status f_start (f_strand *self);
#endif // DIDEROT_HAS_START_METHOD
static diderot::strand_status f_update (globals *glob, f_strand *self);
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
    diderot::strand_status strand_update (globals *glob, index_t ix)
    {
        return f_update(glob, this->strand(ix));
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
    defined_inputs _definedInp;
    world ();
    ~world ();
    bool init ();
    bool alloc (int32_t base[3], uint32_t size[3]);
    bool create_strands ();
    uint32_t run (uint32_t max_nsteps);
    void swap_state ();
};
// ***** Begin synthesized operations *****

inline vec3 vfloor3 (vec3 v)
{
    return __extension__ (vec3){diderot::floor(v[0]), diderot::floor(v[1]), diderot::floor(v[2]), 0.e0f};
}
inline float vdot6 (vec6 u, vec6 v)
{
    vec6 w = u * v;
    return w[0] + w[1] + w[2] + w[3] + w[4] + w[5];
}
inline diderot::array< int, 3 > vtoi3 (vec3 v0)
{
    diderot::array< int, 3 > res = {int32_t(v0[0]),int32_t(v0[1]),int32_t(v0[2]),};
    return res;
}
inline vec3 vload3 (const float *vp)
{
    return __extension__ (vec3){vp[0], vp[1], vp[2], 0.e0f};
}
inline vec3 vcons3 (float r0, float r1, float r2)
{
    return __extension__ (vec3){r0, r1, r2, 0.e0f};
}
template <typename TY, const int VOXSZ>
inline bool inside3Ds3 (vec3 x0, diderot::image3d< float, TY, VOXSZ > img)
{
    return 2 < x0[0] && x0[2] < img.size(2) - 3 && 2 < x0[2] && x0[1] < img.size(1) - 3 && 2 < x0[1] && x0[0] < img.size(
        0) - 3;
}
inline void vpack3 (tensor_3 &dst, vec3 v0)
{
    dst._data[0] = v0[0];
    dst._data[1] = v0[1];
    dst._data[2] = v0[2];
}
template <typename TY, const int VOXSZ>
inline tensor_ref_3_3 world2image (diderot::image3d< float, TY, VOXSZ > const & img)
{
    return tensor_ref_3_3(img.world2image());
}
inline vec6 vcons6 (float r0, float r1, float r2, float r3, float r4, float r5)
{
    return __extension__ (vec6){r0, r1, r2, r3, r4, r5, 0.e0f, 0.e0f};
}
inline float vdot3 (vec3 u, vec3 v)
{
    vec3 w = u * v;
    return w[0] + w[1] + w[2];
}
template <typename TY, const int VOXSZ>
inline tensor_ref_3 translate (diderot::image3d< float, TY, VOXSZ > const & img)
{
    return tensor_ref_3(img.translate());
}
// ***** End synthesized operations *****

extern "C" bool Diderot_input_set_by_name_A (Diderot_world_t *cWrld, const char *s)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    if (wrld->_globals->gv_A.load(wrld, s)) {
        return true;
    }
    return false;
}
extern "C" bool Diderot_input_set_A (Diderot_world_t *cWrld, Nrrd *nin)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    if (wrld->_globals->gv_A.load(wrld, nin)) {
        return true;
    }
    return false;
}
static bool check_defined (world *wrld)
{
    if (!wrld->_definedInp.gv_A) {
        biffMsgAdd(wrld->_errors, "undefined input \"A\"\n");
        return true;
    }
    return false;
}
static void init_defined_inputs (world *wrld)
{
    wrld->_definedInp.gv_A = false;
}
static void init_defaults (globals *glob)
{
}
static bool init_globals (world *wrld)
{
    globals *glob = wrld->_globals;
    glob->gv_shape = 10;
    glob->gv_A.register_global();
    return false;
}
static void f_init (globals *glob, f_strand *self, int32_t p_i_0, int32_t p_j_1, int32_t p_k_2)
{
    float l_r_3 = 0.4e1f * static_cast<float>(1 / glob->gv_shape);
    vpack3(self->sv_val,
        vcons3(l_r_3 * static_cast<float>(p_i_0), l_r_3 * static_cast<float>(p_j_1), l_r_3 * static_cast<float>(p_k_2)) + vcons3(
            -0.2e1f, -0.2e1f, -0.2e1f));
    self->sv_result[0] = 0.e0f;
    self->sv_result[1] = 0.e0f;
    self->sv_result[2] = 0.e0f;
    self->sv_result[3] = 0.e0f;
    self->sv_result[4] = 0.e0f;
    self->sv_result[5] = 0.e0f;
    self->sv_result[6] = 0.e0f;
    self->sv_result[7] = 0.e0f;
    self->sv_result[8] = 0.e0f;
}
static diderot::strand_status f_update (globals *glob, f_strand *self)
{
    tensor_3_3 l_result_440;
    tensor_ref_3_3 l_Mtransform_5 = world2image(glob->gv_A);
    vec3 v_6 = vcons3(vdot3(vload3(l_Mtransform_5.last(0).addr(0)), vload3(tensor_ref_3(self->sv_val).addr(0))),
        vdot3(vload3(l_Mtransform_5.last(3).addr(0)), vload3(tensor_ref_3(self->sv_val).addr(0))),
        vdot3(vload3(l_Mtransform_5.last(6).addr(0)), vload3(tensor_ref_3(self->sv_val).addr(0)))) + vload3(
        translate(glob->gv_A).addr(0));
    vec3 v_7 = v_6;
    if (inside3Ds3(v_6, glob->gv_A)) {
        vec3 v_8 = vfloor3(v_7);
        vec3 v_9 = v_7 - v_8;
        diderot::array< int32_t, 3 > l_n_10 = vtoi3(v_8);
        vec3 v_11 = vcons3(l_Mtransform_5[0], l_Mtransform_5[3], l_Mtransform_5[6]);
        vec3 v_12 = vcons3(l_Mtransform_5[1], l_Mtransform_5[4], l_Mtransform_5[7]);
        vec3 v_13 = vcons3(l_Mtransform_5[2], l_Mtransform_5[5], l_Mtransform_5[8]);
        int32_t l_idx_14 = l_n_10[0] + -2;
        int32_t l_idx_15 = l_n_10[1] + -2;
        int32_t l_idx_16 = l_n_10[2] + -2;
        int32_t l_sz_17 = glob->gv_A.size(0);
        int32_t l_sz_18 = glob->gv_A.size(1);
        int32_t l_mulRes_19 = l_sz_18 * l_idx_16;
        int32_t l_mulRes_20 = l_sz_17 * (l_idx_15 + l_mulRes_19);
        int32_t l_offp_21 = 3 * (l_idx_14 + l_mulRes_20);
        int32_t l_addRes_22 = l_idx_14 + 1;
        int32_t l_offp_23 = 3 * (l_addRes_22 + l_mulRes_20);
        int32_t l_addRes_24 = l_idx_14 + 2;
        int32_t l_offp_25 = 3 * (l_addRes_24 + l_mulRes_20);
        int32_t l_addRes_26 = l_idx_14 + 3;
        int32_t l_offp_27 = 3 * (l_addRes_26 + l_mulRes_20);
        int32_t l_addRes_28 = l_idx_14 + 4;
        int32_t l_offp_29 = 3 * (l_addRes_28 + l_mulRes_20);
        int32_t l_addRes_30 = l_idx_14 + 5;
        int32_t l_offp_31 = 3 * (l_addRes_30 + l_mulRes_20);
        vec6 v_32 = vcons6(glob->gv_A[l_offp_21], glob->gv_A[l_offp_23], glob->gv_A[l_offp_25], glob->gv_A[l_offp_27],
            glob->gv_A[l_offp_29], glob->gv_A[l_offp_31]);
        int32_t l_addRes_33 = l_idx_15 + 1;
        int32_t l_mulRes_34 = l_sz_17 * (l_addRes_33 + l_mulRes_19);
        int32_t l_offp_35 = 3 * (l_idx_14 + l_mulRes_34);
        int32_t l_offp_36 = 3 * (l_addRes_22 + l_mulRes_34);
        int32_t l_offp_37 = 3 * (l_addRes_24 + l_mulRes_34);
        int32_t l_offp_38 = 3 * (l_addRes_26 + l_mulRes_34);
        int32_t l_offp_39 = 3 * (l_addRes_28 + l_mulRes_34);
        int32_t l_offp_40 = 3 * (l_addRes_30 + l_mulRes_34);
        vec6 v_41 = vcons6(glob->gv_A[l_offp_35], glob->gv_A[l_offp_36], glob->gv_A[l_offp_37], glob->gv_A[l_offp_38],
            glob->gv_A[l_offp_39], glob->gv_A[l_offp_40]);
        int32_t l_addRes_42 = l_idx_15 + 2;
        int32_t l_mulRes_43 = l_sz_17 * (l_addRes_42 + l_mulRes_19);
        int32_t l_offp_44 = 3 * (l_idx_14 + l_mulRes_43);
        int32_t l_offp_45 = 3 * (l_addRes_22 + l_mulRes_43);
        int32_t l_offp_46 = 3 * (l_addRes_24 + l_mulRes_43);
        int32_t l_offp_47 = 3 * (l_addRes_26 + l_mulRes_43);
        int32_t l_offp_48 = 3 * (l_addRes_28 + l_mulRes_43);
        int32_t l_offp_49 = 3 * (l_addRes_30 + l_mulRes_43);
        vec6 v_50 = vcons6(glob->gv_A[l_offp_44], glob->gv_A[l_offp_45], glob->gv_A[l_offp_46], glob->gv_A[l_offp_47],
            glob->gv_A[l_offp_48], glob->gv_A[l_offp_49]);
        int32_t l_addRes_51 = l_idx_15 + 3;
        int32_t l_mulRes_52 = l_sz_17 * (l_addRes_51 + l_mulRes_19);
        int32_t l_offp_53 = 3 * (l_idx_14 + l_mulRes_52);
        int32_t l_offp_54 = 3 * (l_addRes_22 + l_mulRes_52);
        int32_t l_offp_55 = 3 * (l_addRes_24 + l_mulRes_52);
        int32_t l_offp_56 = 3 * (l_addRes_26 + l_mulRes_52);
        int32_t l_offp_57 = 3 * (l_addRes_28 + l_mulRes_52);
        int32_t l_offp_58 = 3 * (l_addRes_30 + l_mulRes_52);
        vec6 v_59 = vcons6(glob->gv_A[l_offp_53], glob->gv_A[l_offp_54], glob->gv_A[l_offp_55], glob->gv_A[l_offp_56],
            glob->gv_A[l_offp_57], glob->gv_A[l_offp_58]);
        int32_t l_addRes_60 = l_idx_15 + 4;
        int32_t l_mulRes_61 = l_sz_17 * (l_addRes_60 + l_mulRes_19);
        int32_t l_offp_62 = 3 * (l_idx_14 + l_mulRes_61);
        int32_t l_offp_63 = 3 * (l_addRes_22 + l_mulRes_61);
        int32_t l_offp_64 = 3 * (l_addRes_24 + l_mulRes_61);
        int32_t l_offp_65 = 3 * (l_addRes_26 + l_mulRes_61);
        int32_t l_offp_66 = 3 * (l_addRes_28 + l_mulRes_61);
        int32_t l_offp_67 = 3 * (l_addRes_30 + l_mulRes_61);
        vec6 v_68 = vcons6(glob->gv_A[l_offp_62], glob->gv_A[l_offp_63], glob->gv_A[l_offp_64], glob->gv_A[l_offp_65],
            glob->gv_A[l_offp_66], glob->gv_A[l_offp_67]);
        int32_t l_addRes_69 = l_idx_15 + 5;
        int32_t l_mulRes_70 = l_sz_17 * (l_addRes_69 + l_mulRes_19);
        int32_t l_offp_71 = 3 * (l_idx_14 + l_mulRes_70);
        int32_t l_offp_72 = 3 * (l_addRes_22 + l_mulRes_70);
        int32_t l_offp_73 = 3 * (l_addRes_24 + l_mulRes_70);
        int32_t l_offp_74 = 3 * (l_addRes_26 + l_mulRes_70);
        int32_t l_offp_75 = 3 * (l_addRes_28 + l_mulRes_70);
        int32_t l_offp_76 = 3 * (l_addRes_30 + l_mulRes_70);
        vec6 v_77 = vcons6(glob->gv_A[l_offp_71], glob->gv_A[l_offp_72], glob->gv_A[l_offp_73], glob->gv_A[l_offp_74],
            glob->gv_A[l_offp_75], glob->gv_A[l_offp_76]);
        int32_t l_mulRes_78 = l_sz_18 * (l_idx_16 + 1);
        int32_t l_mulRes_79 = l_sz_17 * (l_idx_15 + l_mulRes_78);
        int32_t l_offp_80 = 3 * (l_idx_14 + l_mulRes_79);
        int32_t l_offp_81 = 3 * (l_addRes_22 + l_mulRes_79);
        int32_t l_offp_82 = 3 * (l_addRes_24 + l_mulRes_79);
        int32_t l_offp_83 = 3 * (l_addRes_26 + l_mulRes_79);
        int32_t l_offp_84 = 3 * (l_addRes_28 + l_mulRes_79);
        int32_t l_offp_85 = 3 * (l_addRes_30 + l_mulRes_79);
        vec6 v_86 = vcons6(glob->gv_A[l_offp_80], glob->gv_A[l_offp_81], glob->gv_A[l_offp_82], glob->gv_A[l_offp_83],
            glob->gv_A[l_offp_84], glob->gv_A[l_offp_85]);
        int32_t l_mulRes_87 = l_sz_17 * (l_addRes_33 + l_mulRes_78);
        int32_t l_offp_88 = 3 * (l_idx_14 + l_mulRes_87);
        int32_t l_offp_89 = 3 * (l_addRes_22 + l_mulRes_87);
        int32_t l_offp_90 = 3 * (l_addRes_24 + l_mulRes_87);
        int32_t l_offp_91 = 3 * (l_addRes_26 + l_mulRes_87);
        int32_t l_offp_92 = 3 * (l_addRes_28 + l_mulRes_87);
        int32_t l_offp_93 = 3 * (l_addRes_30 + l_mulRes_87);
        vec6 v_94 = vcons6(glob->gv_A[l_offp_88], glob->gv_A[l_offp_89], glob->gv_A[l_offp_90], glob->gv_A[l_offp_91],
            glob->gv_A[l_offp_92], glob->gv_A[l_offp_93]);
        int32_t l_mulRes_95 = l_sz_17 * (l_addRes_42 + l_mulRes_78);
        int32_t l_offp_96 = 3 * (l_idx_14 + l_mulRes_95);
        int32_t l_offp_97 = 3 * (l_addRes_22 + l_mulRes_95);
        int32_t l_offp_98 = 3 * (l_addRes_24 + l_mulRes_95);
        int32_t l_offp_99 = 3 * (l_addRes_26 + l_mulRes_95);
        int32_t l_offp_100 = 3 * (l_addRes_28 + l_mulRes_95);
        int32_t l_offp_101 = 3 * (l_addRes_30 + l_mulRes_95);
        vec6 v_102 = vcons6(glob->gv_A[l_offp_96], glob->gv_A[l_offp_97], glob->gv_A[l_offp_98], glob->gv_A[l_offp_99],
            glob->gv_A[l_offp_100], glob->gv_A[l_offp_101]);
        int32_t l_mulRes_103 = l_sz_17 * (l_addRes_51 + l_mulRes_78);
        int32_t l_offp_104 = 3 * (l_idx_14 + l_mulRes_103);
        int32_t l_offp_105 = 3 * (l_addRes_22 + l_mulRes_103);
        int32_t l_offp_106 = 3 * (l_addRes_24 + l_mulRes_103);
        int32_t l_offp_107 = 3 * (l_addRes_26 + l_mulRes_103);
        int32_t l_offp_108 = 3 * (l_addRes_28 + l_mulRes_103);
        int32_t l_offp_109 = 3 * (l_addRes_30 + l_mulRes_103);
        vec6 v_110 = vcons6(glob->gv_A[l_offp_104], glob->gv_A[l_offp_105], glob->gv_A[l_offp_106],
            glob->gv_A[l_offp_107], glob->gv_A[l_offp_108], glob->gv_A[l_offp_109]);
        int32_t l_mulRes_111 = l_sz_17 * (l_addRes_60 + l_mulRes_78);
        int32_t l_offp_112 = 3 * (l_idx_14 + l_mulRes_111);
        int32_t l_offp_113 = 3 * (l_addRes_22 + l_mulRes_111);
        int32_t l_offp_114 = 3 * (l_addRes_24 + l_mulRes_111);
        int32_t l_offp_115 = 3 * (l_addRes_26 + l_mulRes_111);
        int32_t l_offp_116 = 3 * (l_addRes_28 + l_mulRes_111);
        int32_t l_offp_117 = 3 * (l_addRes_30 + l_mulRes_111);
        vec6 v_118 = vcons6(glob->gv_A[l_offp_112], glob->gv_A[l_offp_113], glob->gv_A[l_offp_114],
            glob->gv_A[l_offp_115], glob->gv_A[l_offp_116], glob->gv_A[l_offp_117]);
        int32_t l_mulRes_119 = l_sz_17 * (l_addRes_69 + l_mulRes_78);
        int32_t l_offp_120 = 3 * (l_idx_14 + l_mulRes_119);
        int32_t l_offp_121 = 3 * (l_addRes_22 + l_mulRes_119);
        int32_t l_offp_122 = 3 * (l_addRes_24 + l_mulRes_119);
        int32_t l_offp_123 = 3 * (l_addRes_26 + l_mulRes_119);
        int32_t l_offp_124 = 3 * (l_addRes_28 + l_mulRes_119);
        int32_t l_offp_125 = 3 * (l_addRes_30 + l_mulRes_119);
        vec6 v_126 = vcons6(glob->gv_A[l_offp_120], glob->gv_A[l_offp_121], glob->gv_A[l_offp_122],
            glob->gv_A[l_offp_123], glob->gv_A[l_offp_124], glob->gv_A[l_offp_125]);
        int32_t l_mulRes_127 = l_sz_18 * (l_idx_16 + 2);
        int32_t l_mulRes_128 = l_sz_17 * (l_idx_15 + l_mulRes_127);
        int32_t l_offp_129 = 3 * (l_idx_14 + l_mulRes_128);
        int32_t l_offp_130 = 3 * (l_addRes_22 + l_mulRes_128);
        int32_t l_offp_131 = 3 * (l_addRes_24 + l_mulRes_128);
        int32_t l_offp_132 = 3 * (l_addRes_26 + l_mulRes_128);
        int32_t l_offp_133 = 3 * (l_addRes_28 + l_mulRes_128);
        int32_t l_offp_134 = 3 * (l_addRes_30 + l_mulRes_128);
        vec6 v_135 = vcons6(glob->gv_A[l_offp_129], glob->gv_A[l_offp_130], glob->gv_A[l_offp_131],
            glob->gv_A[l_offp_132], glob->gv_A[l_offp_133], glob->gv_A[l_offp_134]);
        int32_t l_mulRes_136 = l_sz_17 * (l_addRes_33 + l_mulRes_127);
        int32_t l_offp_137 = 3 * (l_idx_14 + l_mulRes_136);
        int32_t l_offp_138 = 3 * (l_addRes_22 + l_mulRes_136);
        int32_t l_offp_139 = 3 * (l_addRes_24 + l_mulRes_136);
        int32_t l_offp_140 = 3 * (l_addRes_26 + l_mulRes_136);
        int32_t l_offp_141 = 3 * (l_addRes_28 + l_mulRes_136);
        int32_t l_offp_142 = 3 * (l_addRes_30 + l_mulRes_136);
        vec6 v_143 = vcons6(glob->gv_A[l_offp_137], glob->gv_A[l_offp_138], glob->gv_A[l_offp_139],
            glob->gv_A[l_offp_140], glob->gv_A[l_offp_141], glob->gv_A[l_offp_142]);
        int32_t l_mulRes_144 = l_sz_17 * (l_addRes_42 + l_mulRes_127);
        int32_t l_offp_145 = 3 * (l_idx_14 + l_mulRes_144);
        int32_t l_offp_146 = 3 * (l_addRes_22 + l_mulRes_144);
        int32_t l_offp_147 = 3 * (l_addRes_24 + l_mulRes_144);
        int32_t l_offp_148 = 3 * (l_addRes_26 + l_mulRes_144);
        int32_t l_offp_149 = 3 * (l_addRes_28 + l_mulRes_144);
        int32_t l_offp_150 = 3 * (l_addRes_30 + l_mulRes_144);
        vec6 v_151 = vcons6(glob->gv_A[l_offp_145], glob->gv_A[l_offp_146], glob->gv_A[l_offp_147],
            glob->gv_A[l_offp_148], glob->gv_A[l_offp_149], glob->gv_A[l_offp_150]);
        int32_t l_mulRes_152 = l_sz_17 * (l_addRes_51 + l_mulRes_127);
        int32_t l_offp_153 = 3 * (l_idx_14 + l_mulRes_152);
        int32_t l_offp_154 = 3 * (l_addRes_22 + l_mulRes_152);
        int32_t l_offp_155 = 3 * (l_addRes_24 + l_mulRes_152);
        int32_t l_offp_156 = 3 * (l_addRes_26 + l_mulRes_152);
        int32_t l_offp_157 = 3 * (l_addRes_28 + l_mulRes_152);
        int32_t l_offp_158 = 3 * (l_addRes_30 + l_mulRes_152);
        vec6 v_159 = vcons6(glob->gv_A[l_offp_153], glob->gv_A[l_offp_154], glob->gv_A[l_offp_155],
            glob->gv_A[l_offp_156], glob->gv_A[l_offp_157], glob->gv_A[l_offp_158]);
        int32_t l_mulRes_160 = l_sz_17 * (l_addRes_60 + l_mulRes_127);
        int32_t l_offp_161 = 3 * (l_idx_14 + l_mulRes_160);
        int32_t l_offp_162 = 3 * (l_addRes_22 + l_mulRes_160);
        int32_t l_offp_163 = 3 * (l_addRes_24 + l_mulRes_160);
        int32_t l_offp_164 = 3 * (l_addRes_26 + l_mulRes_160);
        int32_t l_offp_165 = 3 * (l_addRes_28 + l_mulRes_160);
        int32_t l_offp_166 = 3 * (l_addRes_30 + l_mulRes_160);
        vec6 v_167 = vcons6(glob->gv_A[l_offp_161], glob->gv_A[l_offp_162], glob->gv_A[l_offp_163],
            glob->gv_A[l_offp_164], glob->gv_A[l_offp_165], glob->gv_A[l_offp_166]);
        int32_t l_mulRes_168 = l_sz_17 * (l_addRes_69 + l_mulRes_127);
        int32_t l_offp_169 = 3 * (l_idx_14 + l_mulRes_168);
        int32_t l_offp_170 = 3 * (l_addRes_22 + l_mulRes_168);
        int32_t l_offp_171 = 3 * (l_addRes_24 + l_mulRes_168);
        int32_t l_offp_172 = 3 * (l_addRes_26 + l_mulRes_168);
        int32_t l_offp_173 = 3 * (l_addRes_28 + l_mulRes_168);
        int32_t l_offp_174 = 3 * (l_addRes_30 + l_mulRes_168);
        vec6 v_175 = vcons6(glob->gv_A[l_offp_169], glob->gv_A[l_offp_170], glob->gv_A[l_offp_171],
            glob->gv_A[l_offp_172], glob->gv_A[l_offp_173], glob->gv_A[l_offp_174]);
        int32_t l_mulRes_176 = l_sz_18 * (l_idx_16 + 3);
        int32_t l_mulRes_177 = l_sz_17 * (l_idx_15 + l_mulRes_176);
        int32_t l_offp_178 = 3 * (l_idx_14 + l_mulRes_177);
        int32_t l_offp_179 = 3 * (l_addRes_22 + l_mulRes_177);
        int32_t l_offp_180 = 3 * (l_addRes_24 + l_mulRes_177);
        int32_t l_offp_181 = 3 * (l_addRes_26 + l_mulRes_177);
        int32_t l_offp_182 = 3 * (l_addRes_28 + l_mulRes_177);
        int32_t l_offp_183 = 3 * (l_addRes_30 + l_mulRes_177);
        vec6 v_184 = vcons6(glob->gv_A[l_offp_178], glob->gv_A[l_offp_179], glob->gv_A[l_offp_180],
            glob->gv_A[l_offp_181], glob->gv_A[l_offp_182], glob->gv_A[l_offp_183]);
        int32_t l_mulRes_185 = l_sz_17 * (l_addRes_33 + l_mulRes_176);
        int32_t l_offp_186 = 3 * (l_idx_14 + l_mulRes_185);
        int32_t l_offp_187 = 3 * (l_addRes_22 + l_mulRes_185);
        int32_t l_offp_188 = 3 * (l_addRes_24 + l_mulRes_185);
        int32_t l_offp_189 = 3 * (l_addRes_26 + l_mulRes_185);
        int32_t l_offp_190 = 3 * (l_addRes_28 + l_mulRes_185);
        int32_t l_offp_191 = 3 * (l_addRes_30 + l_mulRes_185);
        vec6 v_192 = vcons6(glob->gv_A[l_offp_186], glob->gv_A[l_offp_187], glob->gv_A[l_offp_188],
            glob->gv_A[l_offp_189], glob->gv_A[l_offp_190], glob->gv_A[l_offp_191]);
        int32_t l_mulRes_193 = l_sz_17 * (l_addRes_42 + l_mulRes_176);
        int32_t l_offp_194 = 3 * (l_idx_14 + l_mulRes_193);
        int32_t l_offp_195 = 3 * (l_addRes_22 + l_mulRes_193);
        int32_t l_offp_196 = 3 * (l_addRes_24 + l_mulRes_193);
        int32_t l_offp_197 = 3 * (l_addRes_26 + l_mulRes_193);
        int32_t l_offp_198 = 3 * (l_addRes_28 + l_mulRes_193);
        int32_t l_offp_199 = 3 * (l_addRes_30 + l_mulRes_193);
        vec6 v_200 = vcons6(glob->gv_A[l_offp_194], glob->gv_A[l_offp_195], glob->gv_A[l_offp_196],
            glob->gv_A[l_offp_197], glob->gv_A[l_offp_198], glob->gv_A[l_offp_199]);
        int32_t l_mulRes_201 = l_sz_17 * (l_addRes_51 + l_mulRes_176);
        int32_t l_offp_202 = 3 * (l_idx_14 + l_mulRes_201);
        int32_t l_offp_203 = 3 * (l_addRes_22 + l_mulRes_201);
        int32_t l_offp_204 = 3 * (l_addRes_24 + l_mulRes_201);
        int32_t l_offp_205 = 3 * (l_addRes_26 + l_mulRes_201);
        int32_t l_offp_206 = 3 * (l_addRes_28 + l_mulRes_201);
        int32_t l_offp_207 = 3 * (l_addRes_30 + l_mulRes_201);
        vec6 v_208 = vcons6(glob->gv_A[l_offp_202], glob->gv_A[l_offp_203], glob->gv_A[l_offp_204],
            glob->gv_A[l_offp_205], glob->gv_A[l_offp_206], glob->gv_A[l_offp_207]);
        int32_t l_mulRes_209 = l_sz_17 * (l_addRes_60 + l_mulRes_176);
        int32_t l_offp_210 = 3 * (l_idx_14 + l_mulRes_209);
        int32_t l_offp_211 = 3 * (l_addRes_22 + l_mulRes_209);
        int32_t l_offp_212 = 3 * (l_addRes_24 + l_mulRes_209);
        int32_t l_offp_213 = 3 * (l_addRes_26 + l_mulRes_209);
        int32_t l_offp_214 = 3 * (l_addRes_28 + l_mulRes_209);
        int32_t l_offp_215 = 3 * (l_addRes_30 + l_mulRes_209);
        vec6 v_216 = vcons6(glob->gv_A[l_offp_210], glob->gv_A[l_offp_211], glob->gv_A[l_offp_212],
            glob->gv_A[l_offp_213], glob->gv_A[l_offp_214], glob->gv_A[l_offp_215]);
        int32_t l_mulRes_217 = l_sz_17 * (l_addRes_69 + l_mulRes_176);
        int32_t l_offp_218 = 3 * (l_idx_14 + l_mulRes_217);
        int32_t l_offp_219 = 3 * (l_addRes_22 + l_mulRes_217);
        int32_t l_offp_220 = 3 * (l_addRes_24 + l_mulRes_217);
        int32_t l_offp_221 = 3 * (l_addRes_26 + l_mulRes_217);
        int32_t l_offp_222 = 3 * (l_addRes_28 + l_mulRes_217);
        int32_t l_offp_223 = 3 * (l_addRes_30 + l_mulRes_217);
        vec6 v_224 = vcons6(glob->gv_A[l_offp_218], glob->gv_A[l_offp_219], glob->gv_A[l_offp_220],
            glob->gv_A[l_offp_221], glob->gv_A[l_offp_222], glob->gv_A[l_offp_223]);
        int32_t l_mulRes_225 = l_sz_18 * (l_idx_16 + 4);
        int32_t l_mulRes_226 = l_sz_17 * (l_idx_15 + l_mulRes_225);
        int32_t l_offp_227 = 3 * (l_idx_14 + l_mulRes_226);
        int32_t l_offp_228 = 3 * (l_addRes_22 + l_mulRes_226);
        int32_t l_offp_229 = 3 * (l_addRes_24 + l_mulRes_226);
        int32_t l_offp_230 = 3 * (l_addRes_26 + l_mulRes_226);
        int32_t l_offp_231 = 3 * (l_addRes_28 + l_mulRes_226);
        int32_t l_offp_232 = 3 * (l_addRes_30 + l_mulRes_226);
        vec6 v_233 = vcons6(glob->gv_A[l_offp_227], glob->gv_A[l_offp_228], glob->gv_A[l_offp_229],
            glob->gv_A[l_offp_230], glob->gv_A[l_offp_231], glob->gv_A[l_offp_232]);
        int32_t l_mulRes_234 = l_sz_17 * (l_addRes_33 + l_mulRes_225);
        int32_t l_offp_235 = 3 * (l_idx_14 + l_mulRes_234);
        int32_t l_offp_236 = 3 * (l_addRes_22 + l_mulRes_234);
        int32_t l_offp_237 = 3 * (l_addRes_24 + l_mulRes_234);
        int32_t l_offp_238 = 3 * (l_addRes_26 + l_mulRes_234);
        int32_t l_offp_239 = 3 * (l_addRes_28 + l_mulRes_234);
        int32_t l_offp_240 = 3 * (l_addRes_30 + l_mulRes_234);
        vec6 v_241 = vcons6(glob->gv_A[l_offp_235], glob->gv_A[l_offp_236], glob->gv_A[l_offp_237],
            glob->gv_A[l_offp_238], glob->gv_A[l_offp_239], glob->gv_A[l_offp_240]);
        int32_t l_mulRes_242 = l_sz_17 * (l_addRes_42 + l_mulRes_225);
        int32_t l_offp_243 = 3 * (l_idx_14 + l_mulRes_242);
        int32_t l_offp_244 = 3 * (l_addRes_22 + l_mulRes_242);
        int32_t l_offp_245 = 3 * (l_addRes_24 + l_mulRes_242);
        int32_t l_offp_246 = 3 * (l_addRes_26 + l_mulRes_242);
        int32_t l_offp_247 = 3 * (l_addRes_28 + l_mulRes_242);
        int32_t l_offp_248 = 3 * (l_addRes_30 + l_mulRes_242);
        vec6 v_249 = vcons6(glob->gv_A[l_offp_243], glob->gv_A[l_offp_244], glob->gv_A[l_offp_245],
            glob->gv_A[l_offp_246], glob->gv_A[l_offp_247], glob->gv_A[l_offp_248]);
        int32_t l_mulRes_250 = l_sz_17 * (l_addRes_51 + l_mulRes_225);
        int32_t l_offp_251 = 3 * (l_idx_14 + l_mulRes_250);
        int32_t l_offp_252 = 3 * (l_addRes_22 + l_mulRes_250);
        int32_t l_offp_253 = 3 * (l_addRes_24 + l_mulRes_250);
        int32_t l_offp_254 = 3 * (l_addRes_26 + l_mulRes_250);
        int32_t l_offp_255 = 3 * (l_addRes_28 + l_mulRes_250);
        int32_t l_offp_256 = 3 * (l_addRes_30 + l_mulRes_250);
        vec6 v_257 = vcons6(glob->gv_A[l_offp_251], glob->gv_A[l_offp_252], glob->gv_A[l_offp_253],
            glob->gv_A[l_offp_254], glob->gv_A[l_offp_255], glob->gv_A[l_offp_256]);
        int32_t l_mulRes_258 = l_sz_17 * (l_addRes_60 + l_mulRes_225);
        int32_t l_offp_259 = 3 * (l_idx_14 + l_mulRes_258);
        int32_t l_offp_260 = 3 * (l_addRes_22 + l_mulRes_258);
        int32_t l_offp_261 = 3 * (l_addRes_24 + l_mulRes_258);
        int32_t l_offp_262 = 3 * (l_addRes_26 + l_mulRes_258);
        int32_t l_offp_263 = 3 * (l_addRes_28 + l_mulRes_258);
        int32_t l_offp_264 = 3 * (l_addRes_30 + l_mulRes_258);
        vec6 v_265 = vcons6(glob->gv_A[l_offp_259], glob->gv_A[l_offp_260], glob->gv_A[l_offp_261],
            glob->gv_A[l_offp_262], glob->gv_A[l_offp_263], glob->gv_A[l_offp_264]);
        int32_t l_mulRes_266 = l_sz_17 * (l_addRes_69 + l_mulRes_225);
        int32_t l_offp_267 = 3 * (l_idx_14 + l_mulRes_266);
        int32_t l_offp_268 = 3 * (l_addRes_22 + l_mulRes_266);
        int32_t l_offp_269 = 3 * (l_addRes_24 + l_mulRes_266);
        int32_t l_offp_270 = 3 * (l_addRes_26 + l_mulRes_266);
        int32_t l_offp_271 = 3 * (l_addRes_28 + l_mulRes_266);
        int32_t l_offp_272 = 3 * (l_addRes_30 + l_mulRes_266);
        vec6 v_273 = vcons6(glob->gv_A[l_offp_267], glob->gv_A[l_offp_268], glob->gv_A[l_offp_269],
            glob->gv_A[l_offp_270], glob->gv_A[l_offp_271], glob->gv_A[l_offp_272]);
        int32_t l_mulRes_274 = l_sz_18 * (l_idx_16 + 5);
        int32_t l_mulRes_275 = l_sz_17 * (l_idx_15 + l_mulRes_274);
        int32_t l_offp_276 = 3 * (l_idx_14 + l_mulRes_275);
        int32_t l_offp_277 = 3 * (l_addRes_22 + l_mulRes_275);
        int32_t l_offp_278 = 3 * (l_addRes_24 + l_mulRes_275);
        int32_t l_offp_279 = 3 * (l_addRes_26 + l_mulRes_275);
        int32_t l_offp_280 = 3 * (l_addRes_28 + l_mulRes_275);
        int32_t l_offp_281 = 3 * (l_addRes_30 + l_mulRes_275);
        vec6 v_282 = vcons6(glob->gv_A[l_offp_276], glob->gv_A[l_offp_277], glob->gv_A[l_offp_278],
            glob->gv_A[l_offp_279], glob->gv_A[l_offp_280], glob->gv_A[l_offp_281]);
        int32_t l_mulRes_283 = l_sz_17 * (l_addRes_33 + l_mulRes_274);
        int32_t l_offp_284 = 3 * (l_idx_14 + l_mulRes_283);
        int32_t l_offp_285 = 3 * (l_addRes_22 + l_mulRes_283);
        int32_t l_offp_286 = 3 * (l_addRes_24 + l_mulRes_283);
        int32_t l_offp_287 = 3 * (l_addRes_26 + l_mulRes_283);
        int32_t l_offp_288 = 3 * (l_addRes_28 + l_mulRes_283);
        int32_t l_offp_289 = 3 * (l_addRes_30 + l_mulRes_283);
        vec6 v_290 = vcons6(glob->gv_A[l_offp_284], glob->gv_A[l_offp_285], glob->gv_A[l_offp_286],
            glob->gv_A[l_offp_287], glob->gv_A[l_offp_288], glob->gv_A[l_offp_289]);
        int32_t l_mulRes_291 = l_sz_17 * (l_addRes_42 + l_mulRes_274);
        int32_t l_offp_292 = 3 * (l_idx_14 + l_mulRes_291);
        int32_t l_offp_293 = 3 * (l_addRes_22 + l_mulRes_291);
        int32_t l_offp_294 = 3 * (l_addRes_24 + l_mulRes_291);
        int32_t l_offp_295 = 3 * (l_addRes_26 + l_mulRes_291);
        int32_t l_offp_296 = 3 * (l_addRes_28 + l_mulRes_291);
        int32_t l_offp_297 = 3 * (l_addRes_30 + l_mulRes_291);
        vec6 v_298 = vcons6(glob->gv_A[l_offp_292], glob->gv_A[l_offp_293], glob->gv_A[l_offp_294],
            glob->gv_A[l_offp_295], glob->gv_A[l_offp_296], glob->gv_A[l_offp_297]);
        int32_t l_mulRes_299 = l_sz_17 * (l_addRes_51 + l_mulRes_274);
        int32_t l_offp_300 = 3 * (l_idx_14 + l_mulRes_299);
        int32_t l_offp_301 = 3 * (l_addRes_22 + l_mulRes_299);
        int32_t l_offp_302 = 3 * (l_addRes_24 + l_mulRes_299);
        int32_t l_offp_303 = 3 * (l_addRes_26 + l_mulRes_299);
        int32_t l_offp_304 = 3 * (l_addRes_28 + l_mulRes_299);
        int32_t l_offp_305 = 3 * (l_addRes_30 + l_mulRes_299);
        vec6 v_306 = vcons6(glob->gv_A[l_offp_300], glob->gv_A[l_offp_301], glob->gv_A[l_offp_302],
            glob->gv_A[l_offp_303], glob->gv_A[l_offp_304], glob->gv_A[l_offp_305]);
        int32_t l_mulRes_307 = l_sz_17 * (l_addRes_60 + l_mulRes_274);
        int32_t l_offp_308 = 3 * (l_idx_14 + l_mulRes_307);
        int32_t l_offp_309 = 3 * (l_addRes_22 + l_mulRes_307);
        int32_t l_offp_310 = 3 * (l_addRes_24 + l_mulRes_307);
        int32_t l_offp_311 = 3 * (l_addRes_26 + l_mulRes_307);
        int32_t l_offp_312 = 3 * (l_addRes_28 + l_mulRes_307);
        int32_t l_offp_313 = 3 * (l_addRes_30 + l_mulRes_307);
        vec6 v_314 = vcons6(glob->gv_A[l_offp_308], glob->gv_A[l_offp_309], glob->gv_A[l_offp_310],
            glob->gv_A[l_offp_311], glob->gv_A[l_offp_312], glob->gv_A[l_offp_313]);
        int32_t l_mulRes_315 = l_sz_17 * (l_addRes_69 + l_mulRes_274);
        int32_t l_offp_316 = 3 * (l_idx_14 + l_mulRes_315);
        int32_t l_offp_317 = 3 * (l_addRes_22 + l_mulRes_315);
        int32_t l_offp_318 = 3 * (l_addRes_24 + l_mulRes_315);
        int32_t l_offp_319 = 3 * (l_addRes_26 + l_mulRes_315);
        int32_t l_offp_320 = 3 * (l_addRes_28 + l_mulRes_315);
        int32_t l_offp_321 = 3 * (l_addRes_30 + l_mulRes_315);
        vec6 v_322 = vcons6(glob->gv_A[l_offp_316], glob->gv_A[l_offp_317], glob->gv_A[l_offp_318],
            glob->gv_A[l_offp_319], glob->gv_A[l_offp_320], glob->gv_A[l_offp_321]);
        vec6 v_323 = vcons6(glob->gv_A[l_offp_21 + 1], glob->gv_A[l_offp_23 + 1], glob->gv_A[l_offp_25 + 1],
            glob->gv_A[l_offp_27 + 1], glob->gv_A[l_offp_29 + 1], glob->gv_A[l_offp_31 + 1]);
        vec6 v_324 = vcons6(glob->gv_A[l_offp_35 + 1], glob->gv_A[l_offp_36 + 1], glob->gv_A[l_offp_37 + 1],
            glob->gv_A[l_offp_38 + 1], glob->gv_A[l_offp_39 + 1], glob->gv_A[l_offp_40 + 1]);
        vec6 v_325 = vcons6(glob->gv_A[l_offp_44 + 1], glob->gv_A[l_offp_45 + 1], glob->gv_A[l_offp_46 + 1],
            glob->gv_A[l_offp_47 + 1], glob->gv_A[l_offp_48 + 1], glob->gv_A[l_offp_49 + 1]);
        vec6 v_326 = vcons6(glob->gv_A[l_offp_53 + 1], glob->gv_A[l_offp_54 + 1], glob->gv_A[l_offp_55 + 1],
            glob->gv_A[l_offp_56 + 1], glob->gv_A[l_offp_57 + 1], glob->gv_A[l_offp_58 + 1]);
        vec6 v_327 = vcons6(glob->gv_A[l_offp_62 + 1], glob->gv_A[l_offp_63 + 1], glob->gv_A[l_offp_64 + 1],
            glob->gv_A[l_offp_65 + 1], glob->gv_A[l_offp_66 + 1], glob->gv_A[l_offp_67 + 1]);
        vec6 v_328 = vcons6(glob->gv_A[l_offp_71 + 1], glob->gv_A[l_offp_72 + 1], glob->gv_A[l_offp_73 + 1],
            glob->gv_A[l_offp_74 + 1], glob->gv_A[l_offp_75 + 1], glob->gv_A[l_offp_76 + 1]);
        vec6 v_329 = vcons6(glob->gv_A[l_offp_80 + 1], glob->gv_A[l_offp_81 + 1], glob->gv_A[l_offp_82 + 1],
            glob->gv_A[l_offp_83 + 1], glob->gv_A[l_offp_84 + 1], glob->gv_A[l_offp_85 + 1]);
        vec6 v_330 = vcons6(glob->gv_A[l_offp_88 + 1], glob->gv_A[l_offp_89 + 1], glob->gv_A[l_offp_90 + 1],
            glob->gv_A[l_offp_91 + 1], glob->gv_A[l_offp_92 + 1], glob->gv_A[l_offp_93 + 1]);
        vec6 v_331 = vcons6(glob->gv_A[l_offp_96 + 1], glob->gv_A[l_offp_97 + 1], glob->gv_A[l_offp_98 + 1],
            glob->gv_A[l_offp_99 + 1], glob->gv_A[l_offp_100 + 1], glob->gv_A[l_offp_101 + 1]);
        vec6 v_332 = vcons6(glob->gv_A[l_offp_104 + 1], glob->gv_A[l_offp_105 + 1], glob->gv_A[l_offp_106 + 1],
            glob->gv_A[l_offp_107 + 1], glob->gv_A[l_offp_108 + 1], glob->gv_A[l_offp_109 + 1]);
        vec6 v_333 = vcons6(glob->gv_A[l_offp_112 + 1], glob->gv_A[l_offp_113 + 1], glob->gv_A[l_offp_114 + 1],
            glob->gv_A[l_offp_115 + 1], glob->gv_A[l_offp_116 + 1], glob->gv_A[l_offp_117 + 1]);
        vec6 v_334 = vcons6(glob->gv_A[l_offp_120 + 1], glob->gv_A[l_offp_121 + 1], glob->gv_A[l_offp_122 + 1],
            glob->gv_A[l_offp_123 + 1], glob->gv_A[l_offp_124 + 1], glob->gv_A[l_offp_125 + 1]);
        vec6 v_335 = vcons6(glob->gv_A[l_offp_129 + 1], glob->gv_A[l_offp_130 + 1], glob->gv_A[l_offp_131 + 1],
            glob->gv_A[l_offp_132 + 1], glob->gv_A[l_offp_133 + 1], glob->gv_A[l_offp_134 + 1]);
        vec6 v_336 = vcons6(glob->gv_A[l_offp_137 + 1], glob->gv_A[l_offp_138 + 1], glob->gv_A[l_offp_139 + 1],
            glob->gv_A[l_offp_140 + 1], glob->gv_A[l_offp_141 + 1], glob->gv_A[l_offp_142 + 1]);
        vec6 v_337 = vcons6(glob->gv_A[l_offp_145 + 1], glob->gv_A[l_offp_146 + 1], glob->gv_A[l_offp_147 + 1],
            glob->gv_A[l_offp_148 + 1], glob->gv_A[l_offp_149 + 1], glob->gv_A[l_offp_150 + 1]);
        vec6 v_338 = vcons6(glob->gv_A[l_offp_153 + 1], glob->gv_A[l_offp_154 + 1], glob->gv_A[l_offp_155 + 1],
            glob->gv_A[l_offp_156 + 1], glob->gv_A[l_offp_157 + 1], glob->gv_A[l_offp_158 + 1]);
        vec6 v_339 = vcons6(glob->gv_A[l_offp_161 + 1], glob->gv_A[l_offp_162 + 1], glob->gv_A[l_offp_163 + 1],
            glob->gv_A[l_offp_164 + 1], glob->gv_A[l_offp_165 + 1], glob->gv_A[l_offp_166 + 1]);
        vec6 v_340 = vcons6(glob->gv_A[l_offp_169 + 1], glob->gv_A[l_offp_170 + 1], glob->gv_A[l_offp_171 + 1],
            glob->gv_A[l_offp_172 + 1], glob->gv_A[l_offp_173 + 1], glob->gv_A[l_offp_174 + 1]);
        vec6 v_341 = vcons6(glob->gv_A[l_offp_178 + 1], glob->gv_A[l_offp_179 + 1], glob->gv_A[l_offp_180 + 1],
            glob->gv_A[l_offp_181 + 1], glob->gv_A[l_offp_182 + 1], glob->gv_A[l_offp_183 + 1]);
        vec6 v_342 = vcons6(glob->gv_A[l_offp_186 + 1], glob->gv_A[l_offp_187 + 1], glob->gv_A[l_offp_188 + 1],
            glob->gv_A[l_offp_189 + 1], glob->gv_A[l_offp_190 + 1], glob->gv_A[l_offp_191 + 1]);
        vec6 v_343 = vcons6(glob->gv_A[l_offp_194 + 1], glob->gv_A[l_offp_195 + 1], glob->gv_A[l_offp_196 + 1],
            glob->gv_A[l_offp_197 + 1], glob->gv_A[l_offp_198 + 1], glob->gv_A[l_offp_199 + 1]);
        vec6 v_344 = vcons6(glob->gv_A[l_offp_202 + 1], glob->gv_A[l_offp_203 + 1], glob->gv_A[l_offp_204 + 1],
            glob->gv_A[l_offp_205 + 1], glob->gv_A[l_offp_206 + 1], glob->gv_A[l_offp_207 + 1]);
        vec6 v_345 = vcons6(glob->gv_A[l_offp_210 + 1], glob->gv_A[l_offp_211 + 1], glob->gv_A[l_offp_212 + 1],
            glob->gv_A[l_offp_213 + 1], glob->gv_A[l_offp_214 + 1], glob->gv_A[l_offp_215 + 1]);
        vec6 v_346 = vcons6(glob->gv_A[l_offp_218 + 1], glob->gv_A[l_offp_219 + 1], glob->gv_A[l_offp_220 + 1],
            glob->gv_A[l_offp_221 + 1], glob->gv_A[l_offp_222 + 1], glob->gv_A[l_offp_223 + 1]);
        vec6 v_347 = vcons6(glob->gv_A[l_offp_227 + 1], glob->gv_A[l_offp_228 + 1], glob->gv_A[l_offp_229 + 1],
            glob->gv_A[l_offp_230 + 1], glob->gv_A[l_offp_231 + 1], glob->gv_A[l_offp_232 + 1]);
        vec6 v_348 = vcons6(glob->gv_A[l_offp_235 + 1], glob->gv_A[l_offp_236 + 1], glob->gv_A[l_offp_237 + 1],
            glob->gv_A[l_offp_238 + 1], glob->gv_A[l_offp_239 + 1], glob->gv_A[l_offp_240 + 1]);
        vec6 v_349 = vcons6(glob->gv_A[l_offp_243 + 1], glob->gv_A[l_offp_244 + 1], glob->gv_A[l_offp_245 + 1],
            glob->gv_A[l_offp_246 + 1], glob->gv_A[l_offp_247 + 1], glob->gv_A[l_offp_248 + 1]);
        vec6 v_350 = vcons6(glob->gv_A[l_offp_251 + 1], glob->gv_A[l_offp_252 + 1], glob->gv_A[l_offp_253 + 1],
            glob->gv_A[l_offp_254 + 1], glob->gv_A[l_offp_255 + 1], glob->gv_A[l_offp_256 + 1]);
        vec6 v_351 = vcons6(glob->gv_A[l_offp_259 + 1], glob->gv_A[l_offp_260 + 1], glob->gv_A[l_offp_261 + 1],
            glob->gv_A[l_offp_262 + 1], glob->gv_A[l_offp_263 + 1], glob->gv_A[l_offp_264 + 1]);
        vec6 v_352 = vcons6(glob->gv_A[l_offp_267 + 1], glob->gv_A[l_offp_268 + 1], glob->gv_A[l_offp_269 + 1],
            glob->gv_A[l_offp_270 + 1], glob->gv_A[l_offp_271 + 1], glob->gv_A[l_offp_272 + 1]);
        vec6 v_353 = vcons6(glob->gv_A[l_offp_276 + 1], glob->gv_A[l_offp_277 + 1], glob->gv_A[l_offp_278 + 1],
            glob->gv_A[l_offp_279 + 1], glob->gv_A[l_offp_280 + 1], glob->gv_A[l_offp_281 + 1]);
        vec6 v_354 = vcons6(glob->gv_A[l_offp_284 + 1], glob->gv_A[l_offp_285 + 1], glob->gv_A[l_offp_286 + 1],
            glob->gv_A[l_offp_287 + 1], glob->gv_A[l_offp_288 + 1], glob->gv_A[l_offp_289 + 1]);
        vec6 v_355 = vcons6(glob->gv_A[l_offp_292 + 1], glob->gv_A[l_offp_293 + 1], glob->gv_A[l_offp_294 + 1],
            glob->gv_A[l_offp_295 + 1], glob->gv_A[l_offp_296 + 1], glob->gv_A[l_offp_297 + 1]);
        vec6 v_356 = vcons6(glob->gv_A[l_offp_300 + 1], glob->gv_A[l_offp_301 + 1], glob->gv_A[l_offp_302 + 1],
            glob->gv_A[l_offp_303 + 1], glob->gv_A[l_offp_304 + 1], glob->gv_A[l_offp_305 + 1]);
        vec6 v_357 = vcons6(glob->gv_A[l_offp_308 + 1], glob->gv_A[l_offp_309 + 1], glob->gv_A[l_offp_310 + 1],
            glob->gv_A[l_offp_311 + 1], glob->gv_A[l_offp_312 + 1], glob->gv_A[l_offp_313 + 1]);
        vec6 v_358 = vcons6(glob->gv_A[l_offp_316 + 1], glob->gv_A[l_offp_317 + 1], glob->gv_A[l_offp_318 + 1],
            glob->gv_A[l_offp_319 + 1], glob->gv_A[l_offp_320 + 1], glob->gv_A[l_offp_321 + 1]);
        vec6 v_359 = vcons6(glob->gv_A[l_offp_21 + 2], glob->gv_A[l_offp_23 + 2], glob->gv_A[l_offp_25 + 2],
            glob->gv_A[l_offp_27 + 2], glob->gv_A[l_offp_29 + 2], glob->gv_A[l_offp_31 + 2]);
        vec6 v_360 = vcons6(glob->gv_A[l_offp_35 + 2], glob->gv_A[l_offp_36 + 2], glob->gv_A[l_offp_37 + 2],
            glob->gv_A[l_offp_38 + 2], glob->gv_A[l_offp_39 + 2], glob->gv_A[l_offp_40 + 2]);
        vec6 v_361 = vcons6(glob->gv_A[l_offp_44 + 2], glob->gv_A[l_offp_45 + 2], glob->gv_A[l_offp_46 + 2],
            glob->gv_A[l_offp_47 + 2], glob->gv_A[l_offp_48 + 2], glob->gv_A[l_offp_49 + 2]);
        vec6 v_362 = vcons6(glob->gv_A[l_offp_53 + 2], glob->gv_A[l_offp_54 + 2], glob->gv_A[l_offp_55 + 2],
            glob->gv_A[l_offp_56 + 2], glob->gv_A[l_offp_57 + 2], glob->gv_A[l_offp_58 + 2]);
        vec6 v_363 = vcons6(glob->gv_A[l_offp_62 + 2], glob->gv_A[l_offp_63 + 2], glob->gv_A[l_offp_64 + 2],
            glob->gv_A[l_offp_65 + 2], glob->gv_A[l_offp_66 + 2], glob->gv_A[l_offp_67 + 2]);
        vec6 v_364 = vcons6(glob->gv_A[l_offp_71 + 2], glob->gv_A[l_offp_72 + 2], glob->gv_A[l_offp_73 + 2],
            glob->gv_A[l_offp_74 + 2], glob->gv_A[l_offp_75 + 2], glob->gv_A[l_offp_76 + 2]);
        vec6 v_365 = vcons6(glob->gv_A[l_offp_80 + 2], glob->gv_A[l_offp_81 + 2], glob->gv_A[l_offp_82 + 2],
            glob->gv_A[l_offp_83 + 2], glob->gv_A[l_offp_84 + 2], glob->gv_A[l_offp_85 + 2]);
        vec6 v_366 = vcons6(glob->gv_A[l_offp_88 + 2], glob->gv_A[l_offp_89 + 2], glob->gv_A[l_offp_90 + 2],
            glob->gv_A[l_offp_91 + 2], glob->gv_A[l_offp_92 + 2], glob->gv_A[l_offp_93 + 2]);
        vec6 v_367 = vcons6(glob->gv_A[l_offp_96 + 2], glob->gv_A[l_offp_97 + 2], glob->gv_A[l_offp_98 + 2],
            glob->gv_A[l_offp_99 + 2], glob->gv_A[l_offp_100 + 2], glob->gv_A[l_offp_101 + 2]);
        vec6 v_368 = vcons6(glob->gv_A[l_offp_104 + 2], glob->gv_A[l_offp_105 + 2], glob->gv_A[l_offp_106 + 2],
            glob->gv_A[l_offp_107 + 2], glob->gv_A[l_offp_108 + 2], glob->gv_A[l_offp_109 + 2]);
        vec6 v_369 = vcons6(glob->gv_A[l_offp_112 + 2], glob->gv_A[l_offp_113 + 2], glob->gv_A[l_offp_114 + 2],
            glob->gv_A[l_offp_115 + 2], glob->gv_A[l_offp_116 + 2], glob->gv_A[l_offp_117 + 2]);
        vec6 v_370 = vcons6(glob->gv_A[l_offp_120 + 2], glob->gv_A[l_offp_121 + 2], glob->gv_A[l_offp_122 + 2],
            glob->gv_A[l_offp_123 + 2], glob->gv_A[l_offp_124 + 2], glob->gv_A[l_offp_125 + 2]);
        vec6 v_371 = vcons6(glob->gv_A[l_offp_129 + 2], glob->gv_A[l_offp_130 + 2], glob->gv_A[l_offp_131 + 2],
            glob->gv_A[l_offp_132 + 2], glob->gv_A[l_offp_133 + 2], glob->gv_A[l_offp_134 + 2]);
        vec6 v_372 = vcons6(glob->gv_A[l_offp_137 + 2], glob->gv_A[l_offp_138 + 2], glob->gv_A[l_offp_139 + 2],
            glob->gv_A[l_offp_140 + 2], glob->gv_A[l_offp_141 + 2], glob->gv_A[l_offp_142 + 2]);
        vec6 v_373 = vcons6(glob->gv_A[l_offp_145 + 2], glob->gv_A[l_offp_146 + 2], glob->gv_A[l_offp_147 + 2],
            glob->gv_A[l_offp_148 + 2], glob->gv_A[l_offp_149 + 2], glob->gv_A[l_offp_150 + 2]);
        vec6 v_374 = vcons6(glob->gv_A[l_offp_153 + 2], glob->gv_A[l_offp_154 + 2], glob->gv_A[l_offp_155 + 2],
            glob->gv_A[l_offp_156 + 2], glob->gv_A[l_offp_157 + 2], glob->gv_A[l_offp_158 + 2]);
        vec6 v_375 = vcons6(glob->gv_A[l_offp_161 + 2], glob->gv_A[l_offp_162 + 2], glob->gv_A[l_offp_163 + 2],
            glob->gv_A[l_offp_164 + 2], glob->gv_A[l_offp_165 + 2], glob->gv_A[l_offp_166 + 2]);
        vec6 v_376 = vcons6(glob->gv_A[l_offp_169 + 2], glob->gv_A[l_offp_170 + 2], glob->gv_A[l_offp_171 + 2],
            glob->gv_A[l_offp_172 + 2], glob->gv_A[l_offp_173 + 2], glob->gv_A[l_offp_174 + 2]);
        vec6 v_377 = vcons6(glob->gv_A[l_offp_178 + 2], glob->gv_A[l_offp_179 + 2], glob->gv_A[l_offp_180 + 2],
            glob->gv_A[l_offp_181 + 2], glob->gv_A[l_offp_182 + 2], glob->gv_A[l_offp_183 + 2]);
        vec6 v_378 = vcons6(glob->gv_A[l_offp_186 + 2], glob->gv_A[l_offp_187 + 2], glob->gv_A[l_offp_188 + 2],
            glob->gv_A[l_offp_189 + 2], glob->gv_A[l_offp_190 + 2], glob->gv_A[l_offp_191 + 2]);
        vec6 v_379 = vcons6(glob->gv_A[l_offp_194 + 2], glob->gv_A[l_offp_195 + 2], glob->gv_A[l_offp_196 + 2],
            glob->gv_A[l_offp_197 + 2], glob->gv_A[l_offp_198 + 2], glob->gv_A[l_offp_199 + 2]);
        vec6 v_380 = vcons6(glob->gv_A[l_offp_202 + 2], glob->gv_A[l_offp_203 + 2], glob->gv_A[l_offp_204 + 2],
            glob->gv_A[l_offp_205 + 2], glob->gv_A[l_offp_206 + 2], glob->gv_A[l_offp_207 + 2]);
        vec6 v_381 = vcons6(glob->gv_A[l_offp_210 + 2], glob->gv_A[l_offp_211 + 2], glob->gv_A[l_offp_212 + 2],
            glob->gv_A[l_offp_213 + 2], glob->gv_A[l_offp_214 + 2], glob->gv_A[l_offp_215 + 2]);
        vec6 v_382 = vcons6(glob->gv_A[l_offp_218 + 2], glob->gv_A[l_offp_219 + 2], glob->gv_A[l_offp_220 + 2],
            glob->gv_A[l_offp_221 + 2], glob->gv_A[l_offp_222 + 2], glob->gv_A[l_offp_223 + 2]);
        vec6 v_383 = vcons6(glob->gv_A[l_offp_227 + 2], glob->gv_A[l_offp_228 + 2], glob->gv_A[l_offp_229 + 2],
            glob->gv_A[l_offp_230 + 2], glob->gv_A[l_offp_231 + 2], glob->gv_A[l_offp_232 + 2]);
        vec6 v_384 = vcons6(glob->gv_A[l_offp_235 + 2], glob->gv_A[l_offp_236 + 2], glob->gv_A[l_offp_237 + 2],
            glob->gv_A[l_offp_238 + 2], glob->gv_A[l_offp_239 + 2], glob->gv_A[l_offp_240 + 2]);
        vec6 v_385 = vcons6(glob->gv_A[l_offp_243 + 2], glob->gv_A[l_offp_244 + 2], glob->gv_A[l_offp_245 + 2],
            glob->gv_A[l_offp_246 + 2], glob->gv_A[l_offp_247 + 2], glob->gv_A[l_offp_248 + 2]);
        vec6 v_386 = vcons6(glob->gv_A[l_offp_251 + 2], glob->gv_A[l_offp_252 + 2], glob->gv_A[l_offp_253 + 2],
            glob->gv_A[l_offp_254 + 2], glob->gv_A[l_offp_255 + 2], glob->gv_A[l_offp_256 + 2]);
        vec6 v_387 = vcons6(glob->gv_A[l_offp_259 + 2], glob->gv_A[l_offp_260 + 2], glob->gv_A[l_offp_261 + 2],
            glob->gv_A[l_offp_262 + 2], glob->gv_A[l_offp_263 + 2], glob->gv_A[l_offp_264 + 2]);
        vec6 v_388 = vcons6(glob->gv_A[l_offp_267 + 2], glob->gv_A[l_offp_268 + 2], glob->gv_A[l_offp_269 + 2],
            glob->gv_A[l_offp_270 + 2], glob->gv_A[l_offp_271 + 2], glob->gv_A[l_offp_272 + 2]);
        vec6 v_389 = vcons6(glob->gv_A[l_offp_276 + 2], glob->gv_A[l_offp_277 + 2], glob->gv_A[l_offp_278 + 2],
            glob->gv_A[l_offp_279 + 2], glob->gv_A[l_offp_280 + 2], glob->gv_A[l_offp_281 + 2]);
        vec6 v_390 = vcons6(glob->gv_A[l_offp_284 + 2], glob->gv_A[l_offp_285 + 2], glob->gv_A[l_offp_286 + 2],
            glob->gv_A[l_offp_287 + 2], glob->gv_A[l_offp_288 + 2], glob->gv_A[l_offp_289 + 2]);
        vec6 v_391 = vcons6(glob->gv_A[l_offp_292 + 2], glob->gv_A[l_offp_293 + 2], glob->gv_A[l_offp_294 + 2],
            glob->gv_A[l_offp_295 + 2], glob->gv_A[l_offp_296 + 2], glob->gv_A[l_offp_297 + 2]);
        vec6 v_392 = vcons6(glob->gv_A[l_offp_300 + 2], glob->gv_A[l_offp_301 + 2], glob->gv_A[l_offp_302 + 2],
            glob->gv_A[l_offp_303 + 2], glob->gv_A[l_offp_304 + 2], glob->gv_A[l_offp_305 + 2]);
        vec6 v_393 = vcons6(glob->gv_A[l_offp_308 + 2], glob->gv_A[l_offp_309 + 2], glob->gv_A[l_offp_310 + 2],
            glob->gv_A[l_offp_311 + 2], glob->gv_A[l_offp_312 + 2], glob->gv_A[l_offp_313 + 2]);
        vec6 v_394 = vcons6(glob->gv_A[l_offp_316 + 2], glob->gv_A[l_offp_317 + 2], glob->gv_A[l_offp_318 + 2],
            glob->gv_A[l_offp_319 + 2], glob->gv_A[l_offp_320 + 2], glob->gv_A[l_offp_321 + 2]);
        float l_vZ__395 = v_9[2];
        vec6 v_396 = vcons6(l_vZ__395 + 0.2e1f, l_vZ__395 + 0.1e1f, l_vZ__395, l_vZ__395 - 0.1e1f, l_vZ__395 - 0.2e1f,
            l_vZ__395 - 0.3e1f);
        vec6 v_397 = vcons6(0.961875e1f, 0.1875e-1f, 0.8625e0f, 0.8625e0f, 0.1875e-1f, 0.961875e1f);
        vec6 v_398 = vcons6(-0.23625e2f, 0.4375e1f, 0.e0f, 0.e0f, -0.4375e1f, 0.23625e2f);
        vec6 v_399 = vcons6(0.2334375e2f, -0.1065625e2f, -0.14375e1f, -0.14375e1f, -0.1065625e2f, 0.2334375e2f);
        vec6 v_400 = vcons6(-0.12e2f, 0.1e2f, 0.e0f, 0.e0f, -0.1e2f, 0.12e2f);
        vec6 v_401 = vcons6(0.340625e1f, -0.459375e1f, 0.11875e1f, 0.11875e1f, -0.459375e1f, 0.340625e1f);
        vec6 v_402 = vcons6(-0.508333333333e0f, 0.104166666667e1f, -0.583333333333e0f, 0.583333333333e0f,
            -0.104166666667e1f, 0.508333333333e0f);
        vec6 v_403 = vcons6(0.3125e-1f, -0.9375e-1f, 0.625e-1f, 0.625e-1f, -0.9375e-1f, 0.3125e-1f);
        vec6 v_404 = v_397 + v_396 * (v_398 + v_396 * (v_399 + v_396 * (v_400 + v_396 * (v_401 + v_396 * (v_402 + v_396 * v_403)))));
        vec6 v_405 = vcons6(0.466875e2f, -0.213125e2f, -0.2875e1f, -0.2875e1f, -0.213125e2f, 0.466875e2f);
        vec6 v_406 = vcons6(-0.36e2f, 0.3e2f, 0.e0f, 0.e0f, -0.3e2f, 0.36e2f);
        vec6 v_407 = vcons6(0.13625e2f, -0.18375e2f, 0.475e1f, 0.475e1f, -0.18375e2f, 0.13625e2f);
        vec6 v_408 = vcons6(-0.254166666667e1f, 0.520833333333e1f, -0.291666666667e1f, 0.291666666667e1f,
            -0.520833333333e1f, 0.254166666667e1f);
        vec6 v_409 = vcons6(0.1875e0f, -0.5625e0f, 0.375e0f, 0.375e0f, -0.5625e0f, 0.1875e0f);
        vec6 v_410 = v_398 + v_396 * (v_405 + v_396 * (v_406 + v_396 * (v_407 + v_396 * (v_408 + v_396 * v_409))));
        float l_vY__411 = v_9[1];
        vec6 v_412 = vcons6(l_vY__411 + 0.2e1f, l_vY__411 + 0.1e1f, l_vY__411, l_vY__411 - 0.1e1f, l_vY__411 - 0.2e1f,
            l_vY__411 - 0.3e1f);
        vec6 v_413 = v_397 + v_412 * (v_398 + v_412 * (v_399 + v_412 * (v_400 + v_412 * (v_401 + v_412 * (v_402 + v_412 * v_403)))));
        vec6 v_414 = v_398 + v_412 * (v_405 + v_412 * (v_406 + v_412 * (v_407 + v_412 * (v_408 + v_412 * v_409))));
        float l_vX__415 = v_9[0];
        vec6 v_416 = vcons6(l_vX__415 + 0.2e1f, l_vX__415 + 0.1e1f, l_vX__415, l_vX__415 - 0.1e1f, l_vX__415 - 0.2e1f,
            l_vX__415 - 0.3e1f);
        vec6 v_417 = v_397 + v_416 * (v_398 + v_416 * (v_399 + v_416 * (v_400 + v_416 * (v_401 + v_416 * (v_402 + v_416 * v_403)))));
        vec6 v_418 = v_398 + v_416 * (v_405 + v_416 * (v_406 + v_416 * (v_407 + v_416 * (v_408 + v_416 * v_409))));
        vec6 v_419 = vcons6(vdot6(v_417, v_32), vdot6(v_417, v_41), vdot6(v_417, v_50), vdot6(v_417, v_59),
            vdot6(v_417, v_68), vdot6(v_417, v_77));
        vec6 v_420 = vcons6(vdot6(v_417, v_86), vdot6(v_417, v_94), vdot6(v_417, v_102), vdot6(v_417, v_110),
            vdot6(v_417, v_118), vdot6(v_417, v_126));
        vec6 v_421 = vcons6(vdot6(v_417, v_135), vdot6(v_417, v_143), vdot6(v_417, v_151), vdot6(v_417, v_159),
            vdot6(v_417, v_167), vdot6(v_417, v_175));
        vec6 v_422 = vcons6(vdot6(v_417, v_184), vdot6(v_417, v_192), vdot6(v_417, v_200), vdot6(v_417, v_208),
            vdot6(v_417, v_216), vdot6(v_417, v_224));
        vec6 v_423 = vcons6(vdot6(v_417, v_233), vdot6(v_417, v_241), vdot6(v_417, v_249), vdot6(v_417, v_257),
            vdot6(v_417, v_265), vdot6(v_417, v_273));
        vec6 v_424 = vcons6(vdot6(v_417, v_282), vdot6(v_417, v_290), vdot6(v_417, v_298), vdot6(v_417, v_306),
            vdot6(v_417, v_314), vdot6(v_417, v_322));
        vec3 v_425 = vcons3(
            vdot6(v_404,
                vcons6(
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_32), vdot6(v_418, v_41), vdot6(v_418, v_50), vdot6(v_418, v_59),
                            vdot6(v_418, v_68), vdot6(v_418, v_77))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_86), vdot6(v_418, v_94), vdot6(v_418, v_102), vdot6(v_418, v_110),
                            vdot6(v_418, v_118), vdot6(v_418, v_126))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_135), vdot6(v_418, v_143), vdot6(v_418, v_151), vdot6(v_418, v_159),
                            vdot6(v_418, v_167), vdot6(v_418, v_175))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_184), vdot6(v_418, v_192), vdot6(v_418, v_200), vdot6(v_418, v_208),
                            vdot6(v_418, v_216), vdot6(v_418, v_224))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_233), vdot6(v_418, v_241), vdot6(v_418, v_249), vdot6(v_418, v_257),
                            vdot6(v_418, v_265), vdot6(v_418, v_273))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_282), vdot6(v_418, v_290), vdot6(v_418, v_298), vdot6(v_418, v_306),
                            vdot6(v_418, v_314), vdot6(v_418, v_322))))),
            vdot6(v_404,
                vcons6(vdot6(v_414, v_419), vdot6(v_414, v_420), vdot6(v_414, v_421), vdot6(v_414, v_422),
                    vdot6(v_414, v_423), vdot6(v_414, v_424))),
            vdot6(v_410,
                vcons6(vdot6(v_413, v_419), vdot6(v_413, v_420), vdot6(v_413, v_421), vdot6(v_413, v_422),
                    vdot6(v_413, v_423), vdot6(v_413, v_424))));
        vec6 v_426 = vcons6(vdot6(v_417, v_323), vdot6(v_417, v_324), vdot6(v_417, v_325), vdot6(v_417, v_326),
            vdot6(v_417, v_327), vdot6(v_417, v_328));
        vec6 v_427 = vcons6(vdot6(v_417, v_329), vdot6(v_417, v_330), vdot6(v_417, v_331), vdot6(v_417, v_332),
            vdot6(v_417, v_333), vdot6(v_417, v_334));
        vec6 v_428 = vcons6(vdot6(v_417, v_335), vdot6(v_417, v_336), vdot6(v_417, v_337), vdot6(v_417, v_338),
            vdot6(v_417, v_339), vdot6(v_417, v_340));
        vec6 v_429 = vcons6(vdot6(v_417, v_341), vdot6(v_417, v_342), vdot6(v_417, v_343), vdot6(v_417, v_344),
            vdot6(v_417, v_345), vdot6(v_417, v_346));
        vec6 v_430 = vcons6(vdot6(v_417, v_347), vdot6(v_417, v_348), vdot6(v_417, v_349), vdot6(v_417, v_350),
            vdot6(v_417, v_351), vdot6(v_417, v_352));
        vec6 v_431 = vcons6(vdot6(v_417, v_353), vdot6(v_417, v_354), vdot6(v_417, v_355), vdot6(v_417, v_356),
            vdot6(v_417, v_357), vdot6(v_417, v_358));
        vec3 v_432 = vcons3(
            vdot6(v_404,
                vcons6(
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_323), vdot6(v_418, v_324), vdot6(v_418, v_325), vdot6(v_418, v_326),
                            vdot6(v_418, v_327), vdot6(v_418, v_328))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_329), vdot6(v_418, v_330), vdot6(v_418, v_331), vdot6(v_418, v_332),
                            vdot6(v_418, v_333), vdot6(v_418, v_334))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_335), vdot6(v_418, v_336), vdot6(v_418, v_337), vdot6(v_418, v_338),
                            vdot6(v_418, v_339), vdot6(v_418, v_340))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_341), vdot6(v_418, v_342), vdot6(v_418, v_343), vdot6(v_418, v_344),
                            vdot6(v_418, v_345), vdot6(v_418, v_346))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_347), vdot6(v_418, v_348), vdot6(v_418, v_349), vdot6(v_418, v_350),
                            vdot6(v_418, v_351), vdot6(v_418, v_352))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_353), vdot6(v_418, v_354), vdot6(v_418, v_355), vdot6(v_418, v_356),
                            vdot6(v_418, v_357), vdot6(v_418, v_358))))),
            vdot6(v_404,
                vcons6(vdot6(v_414, v_426), vdot6(v_414, v_427), vdot6(v_414, v_428), vdot6(v_414, v_429),
                    vdot6(v_414, v_430), vdot6(v_414, v_431))),
            vdot6(v_410,
                vcons6(vdot6(v_413, v_426), vdot6(v_413, v_427), vdot6(v_413, v_428), vdot6(v_413, v_429),
                    vdot6(v_413, v_430), vdot6(v_413, v_431))));
        vec6 v_433 = vcons6(vdot6(v_417, v_359), vdot6(v_417, v_360), vdot6(v_417, v_361), vdot6(v_417, v_362),
            vdot6(v_417, v_363), vdot6(v_417, v_364));
        vec6 v_434 = vcons6(vdot6(v_417, v_365), vdot6(v_417, v_366), vdot6(v_417, v_367), vdot6(v_417, v_368),
            vdot6(v_417, v_369), vdot6(v_417, v_370));
        vec6 v_435 = vcons6(vdot6(v_417, v_371), vdot6(v_417, v_372), vdot6(v_417, v_373), vdot6(v_417, v_374),
            vdot6(v_417, v_375), vdot6(v_417, v_376));
        vec6 v_436 = vcons6(vdot6(v_417, v_377), vdot6(v_417, v_378), vdot6(v_417, v_379), vdot6(v_417, v_380),
            vdot6(v_417, v_381), vdot6(v_417, v_382));
        vec6 v_437 = vcons6(vdot6(v_417, v_383), vdot6(v_417, v_384), vdot6(v_417, v_385), vdot6(v_417, v_386),
            vdot6(v_417, v_387), vdot6(v_417, v_388));
        vec6 v_438 = vcons6(vdot6(v_417, v_389), vdot6(v_417, v_390), vdot6(v_417, v_391), vdot6(v_417, v_392),
            vdot6(v_417, v_393), vdot6(v_417, v_394));
        vec3 v_439 = vcons3(
            vdot6(v_404,
                vcons6(
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_359), vdot6(v_418, v_360), vdot6(v_418, v_361), vdot6(v_418, v_362),
                            vdot6(v_418, v_363), vdot6(v_418, v_364))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_365), vdot6(v_418, v_366), vdot6(v_418, v_367), vdot6(v_418, v_368),
                            vdot6(v_418, v_369), vdot6(v_418, v_370))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_371), vdot6(v_418, v_372), vdot6(v_418, v_373), vdot6(v_418, v_374),
                            vdot6(v_418, v_375), vdot6(v_418, v_376))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_377), vdot6(v_418, v_378), vdot6(v_418, v_379), vdot6(v_418, v_380),
                            vdot6(v_418, v_381), vdot6(v_418, v_382))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_383), vdot6(v_418, v_384), vdot6(v_418, v_385), vdot6(v_418, v_386),
                            vdot6(v_418, v_387), vdot6(v_418, v_388))),
                    vdot6(v_413,
                        vcons6(vdot6(v_418, v_389), vdot6(v_418, v_390), vdot6(v_418, v_391), vdot6(v_418, v_392),
                            vdot6(v_418, v_393), vdot6(v_418, v_394))))),
            vdot6(v_404,
                vcons6(vdot6(v_414, v_433), vdot6(v_414, v_434), vdot6(v_414, v_435), vdot6(v_414, v_436),
                    vdot6(v_414, v_437), vdot6(v_414, v_438))),
            vdot6(v_410,
                vcons6(vdot6(v_413, v_433), vdot6(v_413, v_434), vdot6(v_413, v_435), vdot6(v_413, v_436),
                    vdot6(v_413, v_437), vdot6(v_413, v_438))));
        l_result_440[0] = vdot3(v_425, v_11);
        l_result_440[1] = vdot3(v_425, v_12);
        l_result_440[2] = vdot3(v_425, v_13);
        l_result_440[3] = vdot3(v_432, v_11);
        l_result_440[4] = vdot3(v_432, v_12);
        l_result_440[5] = vdot3(v_432, v_13);
        l_result_440[6] = vdot3(v_439, v_11);
        l_result_440[7] = vdot3(v_439, v_12);
        l_result_440[8] = vdot3(v_439, v_13);
    }
    else {
        l_result_440 = tensor_ref_3_3(self->sv_result);
    }
    self->sv_result = l_result_440;
    return diderot::kStabilize;
}
extern "C" bool Diderot_output_get_result (Diderot_world_t *cWrld, Nrrd *nData)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    // Compute sizes of nrrd file
    size_t sizes[4];
    sizes[0] = 9;
    sizes[1] = wrld->_size[2];
    sizes[2] = wrld->_size[1];
    sizes[3] = wrld->_size[0];
    // Allocate nData nrrd
    if (nrrdMaybeAlloc_nva(nData, nrrdTypeFloat, 4, sizes) != 0) {
        char *msg = biffGetDone(NRRD);
        biffMsgAdd(wrld->_errors, msg);
        std::free(msg);
        return true;
    }
    // copy data to output nrrd
    char *cp = reinterpret_cast<char *>(nData->data);
    for (auto ix = wrld->_strands.begin_alive(); ix != wrld->_strands.end_alive(); ix = wrld->_strands.next_alive(ix)) {
        memcpy(cp, &wrld->_strands.strand(ix)->sv_result, 9 * sizeof(float));
        cp += 9 * sizeof(float);
    }
    nData->axis[0].kind = nrrdKind3DMatrix;
    nData->axis[1].kind = nrrdKindSpace;
    nData->axis[2].kind = nrrdKindSpace;
    nData->axis[3].kind = nrrdKindSpace;
    return false;
}
/*---------- begin world-methods.in ----------*/
// Allocate the program's world
//
world::world ()
    : diderot::world_base (ProgramName, true, 3)
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
bool world::alloc (int32_t base[3], uint32_t size[3])
{
    size_t numStrands = 1;
    for (uint32_t i = 0;  i < 3;  i++) {
        numStrands *= size[i];
        this->_base[i] = base[i];
        this->_size[i] = size[i];
    }

    if (this->_verbose) {
        std::cerr << "world::alloc: " << size[0];
        for (uint32_t i = 1;  i < 3;  i++) {
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
    int32_t l__t_441 = glob->gv_shape - 1;
    int lo_0 = 0;
    int lo_1 = 0;
    int lo_2 = 0;
    int32_t base[3] = {lo_0,lo_1,lo_2,};
    uint32_t size[3] = {static_cast<uint32_t>(l__t_441 - lo_0 + 1),static_cast<uint32_t>(l__t_441 - lo_1 + 1),static_cast<uint32_t>(l__t_441 - lo_2 + 1),};
    if (this->alloc(base, size)) {
        return true;
    }
    uint32_t ix = 0;
    for (int i_i_442 = lo_0; i_i_442 <= l__t_441; i_i_442++) {
        for (int i_j_443 = lo_1; i_j_443 <= l__t_441; i_j_443++) {
            for (int i_k_444 = lo_2; i_k_444 <= l__t_441; i_k_444++) {
                f_init(this->_globals, this->_strands.strand(ix), i_i_442, i_j_443, i_k_444);
                ++ix;
            }
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
            sts = this->_strands.strand_update(glob, ix);
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


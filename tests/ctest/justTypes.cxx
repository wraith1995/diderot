/*---------- begin cxx-head.in ----------*/
/*! \file justTypes.cxx
 *
 * Generated from justTypes.diderot.
 *
 * Command: /home/teocollin/gitcode/diderot/bin/diderotc --debug --log --json --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree --double --namespace=justTypes justTypes.diderot
 * Version: master:2016-07-29
 */
/*---------- end cxx-head.in ----------*/


/**** User mandated includes ****/

#include <spatialindex/capi/sidx_api.h> 

#include <spatialindex/capi/sidx_impl.h> 


/**** Diderot defs and includes ****/

#define DIDEROT_STRAND_HAS_CONSTR
/*---------- begin lib-cxx-incl.in ----------*/
#include "justTypes.h"
#include "diderot/diderot.hxx"

#ifdef DIDEROT_ENABLE_LOGGING
#define IF_LOGGING(...)         __VA_ARGS__
#else
#define IF_LOGGING(...)
#endif

static std::string ProgramName = "justTypes";
/*---------- end lib-cxx-incl.in ----------*/

// ***** Begin synthesized types *****

namespace justTypes {
    typedef double vec3 __attribute__ ((vector_size (32)));
    typedef double vec2 __attribute__ ((vector_size (16)));
    struct tensor_ref_3 : public diderot::tensor_ref<double,3> {
        tensor_ref_3 (const double *src);
        tensor_ref_3 (struct tensor_3 const & ten);
        tensor_ref_3 (tensor_ref_3 const & ten);
    };
    struct tensor_3 : public diderot::tensor<double,3> {
        tensor_3 ()
            : diderot::tensor<double,3>()
        { }
        tensor_3 (std::initializer_list< double > const & il)
            : diderot::tensor<double,3>(il)
        { }
        tensor_3 (const double *src)
            : diderot::tensor<double,3>(src)
        { }
        tensor_3 (tensor_3 const & ten)
            : diderot::tensor<double,3>(ten._data)
        { }
        ~tensor_3 () { }
        tensor_3 & operator= (tensor_3 const & src);
        tensor_3 & operator= (tensor_ref_3 const & src);
        tensor_3 & operator= (std::initializer_list< double > const & il);
        tensor_3 & operator= (const double *src);
    };
    struct msh {
        int32_t *indexMap;
        double *coordMap;
        int32_t dim;
        int32_t mapDim;
        int32_t numCells;
        void *index;
        int32_t *con;
        msh operator= (std::string file)
        {
            // No something with the nrrd
        }
    };
    struct fns {
        int32_t *indexMap;
        msh mesh;
        fns operator= (std::string file)
        {
            // No something with the nrrd
        }
        fns loadFem (msh mesh)
        {
            fns space = *this;
            space.mesh = mesh;
            return space;
        }
    };
    struct mesh_cell_msh {
        int32_t cell;
        msh mesh;
        std::ostream & operator<< (std::ostream & os)
        {
            return os << this->cell;
        }
    };
    struct FUNC {
        double *coordMap;
        fns space;
        FUNC operator= (std::string file)
        {
            // No something with the nrrd
        }
        FUNC loadFem (fns space)
        {
            FUNC func = *this;
            func.space = space;
            return func;
        }
    };
    inline tensor_ref_3::tensor_ref_3 (const double *src)
        : diderot::tensor_ref<double,3>(src)
    { }
    inline tensor_ref_3::tensor_ref_3 (struct tensor_3 const & ten)
        : diderot::tensor_ref<double,3>(ten._data)
    { }
    inline tensor_ref_3::tensor_ref_3 (tensor_ref_3 const & ten)
        : diderot::tensor_ref<double,3>(ten._data)
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
    inline tensor_3 & tensor_3::operator= (std::initializer_list< double > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_3 & tensor_3::operator= (const double *src)
    {
        this->copy(src);
        return *this;
    }
    std::ostream & operator<< (std::ostream & os, const mesh_cell_msh &  cell)
    {
        return os << cell.cell;
    }
    mesh_cell_msh makeFem (msh mesh, int32_t cellInt)
    {
        mesh_cell_msh cell;
        cell.cell = cellInt;
        cell.mesh = mesh;
        return cell;
    }
    char * copy_to (const mesh_cell_msh dumb, char *cp)
    {
        size_t nbytes = sizeof(int32_t);
        (std::memcpy)(cp, &(dumb.cell), nbytes);
        return cp + nbytes;
    }
} // namespace justTypes
namespace diderot {
    template <>
    struct dynseq_traits<justTypes::mesh_cell_msh> {
        using value_type = justTypes::mesh_cell_msh;
        using base_type = int32_t;
        static const __details::load_fn_ptr<base_type> *load_fn_tbl;
        static const uint32_t values_per_elem = 1;
    };
    const __details::load_fn_ptr< dynseq_traits< justTypes::mesh_cell_msh >::base_type > *dynseq_traits< justTypes::mesh_cell_msh >::load_fn_tbl = nrrdILoad;
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
namespace justTypes {

static std::string ProgramName = "justTypes";

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
    bool gv_a;
    bool gv_0b036D_intermedateGlobal;
    bool gv_0c036F_intermedateGlobal;
} defined_inputs;
struct globals {
    msh gv_a;
    fns gv_0b036D_intermedateGlobal;
    FUNC gv_0c036F_intermedateGlobal;
    diderot::dynseq< mesh_cell_msh > gv_0cell_a;
    ~globals () { }
};
struct gg_strand {
    tensor_3 sv_pos;
    tensor_3 sv__pos;
    mesh_cell_msh sv_j;
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

inline vec2 vcons2 (double r0, double r1)
{
    return __extension__ (vec2){r0, r1};
}
static std::ostream& operator<< (std::ostream& outs, tensor_ref_3 const & ten)
{
    return outs << "[" << ten._data[0] << "," << ten._data[1] << "," << ten._data[2] << "]";
}
inline vec3 vcons3 (double r0, double r1, double r2)
{
    return __extension__ (vec3){r0, r1, r2, 0.e0};
}
inline void vpack3 (tensor_3 &dst, vec3 v0)
{
    dst._data[0] = v0[0];
    dst._data[1] = v0[1];
    dst._data[2] = v0[2];
}
inline vec3 vscale3 (double s, vec3 v)
{
    return __extension__ (vec3){s, s, s, 0.e0} * v;
}
inline double vdot3 (vec3 u, vec3 v)
{
    vec3 w = u * v;
    return w[0] + w[1] + w[2];
}
// ***** End synthesized operations *****

extern "C" bool justTypes_input_set_a (justTypes_world_t *cWrld, void *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_a = true;
    std::memcpy(&wrld->_globals->gv_a, v, sizeof(msh));
    return false;
}
extern "C" bool justTypes_input_set_b (justTypes_world_t *cWrld, void *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_0b036D_intermedateGlobal = true;
    std::memcpy(&wrld->_globals->gv_0b036D_intermedateGlobal, v, sizeof(fns));
    return false;
}
extern "C" bool justTypes_input_set_c (justTypes_world_t *cWrld, void *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_0c036F_intermedateGlobal = true;
    std::memcpy(&wrld->_globals->gv_0c036F_intermedateGlobal, v, sizeof(FUNC));
    return false;
}
static bool check_defined (world *wrld)
{
    if (!wrld->_definedInp.gv_a) {
        biffMsgAdd(wrld->_errors, "undefined input \"a\"\n");
        return true;
    }
    if (!wrld->_definedInp.gv_0b036D_intermedateGlobal) {
        biffMsgAdd(wrld->_errors, "undefined input \"b\"\n");
        return true;
    }
    if (!wrld->_definedInp.gv_0c036F_intermedateGlobal) {
        biffMsgAdd(wrld->_errors, "undefined input \"c\"\n");
        return true;
    }
    return false;
}
static void init_defined_inputs (world *wrld)
{
    wrld->_definedInp.gv_a = false;
    wrld->_definedInp.gv_0b036D_intermedateGlobal = false;
    wrld->_definedInp.gv_0c036F_intermedateGlobal = false;
}
static void init_defaults (globals *glob)
{
}
static bool init_globals (world *wrld)
{
    diderot::dynseq< mesh_cell_msh > l__t_0;
    globals *glob = wrld->_globals;
    l__t_0 = {};
    int32_t hi_0 = glob->gv_a.numCells - 1;
    for (int32_t i__t_1 = 0; i__t_1 <= hi_0; ++i__t_1) {
        l__t_0 = diderot::dynseq< mesh_cell_msh >::append(l__t_0, makeFem(glob->gv_a, i__t_1));
    }
    glob->gv_0cell_a = l__t_0;
    return false;
}
static void gg_init (gg_strand *self, mesh_cell_msh p_j_2)
{
    vec3 v_3 = vcons3(0.25e0, 0.25e0, 0.25e0);
    vpack3(self->sv_pos, v_3);
    vpack3(self->sv__pos, v_3);
    self->sv_j = p_j_2;
}
static diderot::strand_status gg_update (world *wrld, gg_strand *self)
{
    int32_t l_face_166;
    double l_time_167;
    int32_t l_face_176;
    double l_time_177;
    int32_t l_face_184;
    double l_time_185;
    int32_t l_face_194;
    double l_time_195;
    vec2 v_196;
    int32_t l_face_213;
    double l_time_214;
    int32_t l_face_220;
    double l_time_221;
    int32_t l_face_227;
    double l_time_228;
    int32_t l_face_234;
    double l_time_235;
    vec2 v_236;
    vec3 v_6 = vcons3(0.1e0, -0.1e1, 0.1e0);
    msh l__t_7 = (self->sv_j.mesh);
    vec3 v_8 = vcons3(0.e0, 0.e0, 0.e0);
    vec3 v_9 = vcons3(0.e0, 0.1e1, 0.e0);
    vec3 v_10 = vcons3(-0.e0, -0.e0, -0.e0);
    vec3 v_11 = vcons3(0.57735026919e0, -0.57735026919e0, 0.57735026919e0);
    int32_t l_mulRes_12 = self->sv_j.cell * 4;
    int32_t t_13 = l__t_7.indexMap[l_mulRes_12];
    int32_t l_mulRes_14 = 3 * t_13;
    double l_dof_load_15 = l__t_7.coordMap[l_mulRes_14];
    double l_dof_load_16 = l__t_7.coordMap[1 + l_mulRes_14];
    double l_dof_load_17 = l__t_7.coordMap[2 + l_mulRes_14];
    int32_t t_18 = l__t_7.indexMap[l_mulRes_12 + 1];
    int32_t l_mulRes_19 = 3 * t_18;
    double l_dof_load_20 = l__t_7.coordMap[l_mulRes_19];
    double l_dof_load_21 = l__t_7.coordMap[1 + l_mulRes_19];
    double l_dof_load_22 = l__t_7.coordMap[2 + l_mulRes_19];
    int32_t t_23 = l__t_7.indexMap[l_mulRes_12 + 2];
    int32_t l_mulRes_24 = 3 * t_23;
    double l_dof_load_25 = l__t_7.coordMap[l_mulRes_24];
    double l_dof_load_26 = l__t_7.coordMap[1 + l_mulRes_24];
    double l_dof_load_27 = l__t_7.coordMap[2 + l_mulRes_24];
    int32_t t_28 = l__t_7.indexMap[l_mulRes_12 + 3];
    int32_t l_mulRes_29 = 3 * t_28;
    double l_dof_load_30 = l__t_7.coordMap[l_mulRes_29];
    double l_dof_load_31 = l__t_7.coordMap[1 + l_mulRes_29];
    double l_dof_load_32 = l__t_7.coordMap[2 + l_mulRes_29];
    double l_prod_33 = 0.1e1 * 0.1e1;
    double l_prod_34 = 0.1e1 * l_prod_33;
    double l_basisEval_35 = -0.1e1 * l_prod_34;
    double l_basisEval_36 = 0.1e1 * l_prod_34;
    double l_r_37 = l_dof_load_15 * l_basisEval_35;
    double l_r_38 = l_dof_load_25 * 0.e0;
    double l_r_39 = l_dof_load_30 * 0.e0;
    double l_r_40 = l_r_37 + l_dof_load_20 * l_basisEval_36 + l_r_38 + l_r_39;
    double l_r_41 = l_r_37 + l_dof_load_20 * 0.e0;
    double l_r_42 = l_r_41 + l_dof_load_25 * l_basisEval_36 + l_r_39;
    double l_r_43 = l_r_41 + l_r_38 + l_dof_load_30 * l_basisEval_36;
    double l_r_44 = l_dof_load_16 * l_basisEval_35;
    double l_r_45 = l_dof_load_26 * 0.e0;
    double l_r_46 = l_dof_load_31 * 0.e0;
    double l_r_47 = l_r_44 + l_dof_load_21 * l_basisEval_36 + l_r_45 + l_r_46;
    double l_r_48 = l_r_44 + l_dof_load_21 * 0.e0;
    double l_r_49 = l_r_48 + l_dof_load_26 * l_basisEval_36 + l_r_46;
    double l_r_50 = l_r_48 + l_r_45 + l_dof_load_31 * l_basisEval_36;
    double l_r_51 = l_dof_load_17 * l_basisEval_35;
    double l_r_52 = l_dof_load_27 * 0.e0;
    double l_r_53 = l_dof_load_32 * 0.e0;
    double l_r_54 = l_r_51 + l_dof_load_22 * l_basisEval_36 + l_r_52 + l_r_53;
    double l_r_55 = l_r_51 + l_dof_load_22 * 0.e0;
    double l_r_56 = l_r_55 + l_dof_load_27 * l_basisEval_36 + l_r_53;
    double l_r_57 = l_r_55 + l_r_52 + l_dof_load_32 * l_basisEval_36;
    double l_r_58 = 0.e0 * l_r_40;
    double l_r_59 = 0.e0 * l_r_47;
    double l_r_60 = 0.e0 * l_r_54;
    double l_r_61 = l_r_58 + l_r_59;
    double l_r_62 = l_r_61 + l_r_60;
    double l_r_63 = 0.e0 * l_r_42;
    double l_r_64 = 0.e0 * l_r_49;
    double l_r_65 = 0.e0 * l_r_56;
    double l_r_66 = l_r_63 + l_r_64;
    double l_r_67 = l_r_66 + l_r_65;
    double l_r_68 = 0.e0 * l_r_43;
    double l_r_69 = 0.e0 * l_r_50;
    double l_r_70 = 0.e0 * l_r_57;
    double l_r_71 = l_r_68 + l_r_69;
    double l_r_72 = l_r_71 + l_r_70;
    double l_r_73 = l_r_61 + -0.1e1 * l_r_54;
    double l_r_74 = l_r_66 + -0.1e1 * l_r_56;
    double l_r_75 = l_r_71 + -0.1e1 * l_r_57;
    double l_r_76 = l_r_58 + 0.1e1 * l_r_47 + l_r_60;
    double l_r_77 = l_r_63 + 0.1e1 * l_r_49 + l_r_65;
    double l_r_78 = l_r_68 + 0.1e1 * l_r_50 + l_r_70;
    double l_r_79 = l_r_61 + 0.1e1 * l_r_54;
    double l_r_80 = l_r_66 + 0.1e1 * l_r_56;
    double l_r_81 = l_r_71 + 0.1e1 * l_r_57;
    double l_r_82 = -0.1e1 * l_r_40 + l_r_59 + l_r_60;
    double l_r_83 = -0.1e1 * l_r_42 + l_r_64 + l_r_65;
    double l_r_84 = -0.1e1 * l_r_43 + l_r_69 + l_r_70;
    double l_r_85 = l_r_58 + -0.1e1 * l_r_47 + l_r_60;
    double l_r_86 = l_r_63 + -0.1e1 * l_r_49 + l_r_65;
    double l_r_87 = l_r_68 + -0.1e1 * l_r_50 + l_r_70;
    double l_r_88 = 0.1e1 * l_r_40 + l_r_59 + l_r_60;
    double l_r_89 = 0.1e1 * l_r_42 + l_r_64 + l_r_65;
    double l_r_90 = 0.1e1 * l_r_43 + l_r_69 + l_r_70;
    double l_r_91 = l_r_40 * l_r_67 + l_r_47 * l_r_80 + l_r_54 * l_r_86;
    double l_r_92 = l_r_40 * l_r_72 + l_r_47 * l_r_81 + l_r_54 * l_r_87;
    double l_r_93 = l_r_40 * l_r_74 + l_r_47 * l_r_67 + l_r_54 * l_r_89;
    double l_r_94 = l_r_40 * l_r_75 + l_r_47 * l_r_72 + l_r_54 * l_r_90;
    double l_r_95 = l_r_40 * l_r_77 + l_r_47 * l_r_83 + l_r_54 * l_r_67;
    double l_r_96 = l_r_40 * l_r_78 + l_r_47 * l_r_84 + l_r_54 * l_r_72;
    double l_r_97 = l_r_42 * l_r_62 + l_r_49 * l_r_79 + l_r_56 * l_r_85;
    double l_r_98 = l_r_42 * l_r_72 + l_r_49 * l_r_81 + l_r_56 * l_r_87;
    double l_r_99 = l_r_42 * l_r_73 + l_r_49 * l_r_62 + l_r_56 * l_r_88;
    double l_r_100 = l_r_42 * l_r_75 + l_r_49 * l_r_72 + l_r_56 * l_r_90;
    double l_r_101 = l_r_42 * l_r_76 + l_r_49 * l_r_82 + l_r_56 * l_r_62;
    double l_r_102 = l_r_42 * l_r_78 + l_r_49 * l_r_84 + l_r_56 * l_r_72;
    double l_r_103 = l_r_43 * l_r_62 + l_r_50 * l_r_79 + l_r_57 * l_r_85;
    double l_r_104 = l_r_43 * l_r_67 + l_r_50 * l_r_80 + l_r_57 * l_r_86;
    double l_r_105 = l_r_43 * l_r_73 + l_r_50 * l_r_62 + l_r_57 * l_r_88;
    double l_r_106 = l_r_43 * l_r_74 + l_r_50 * l_r_67 + l_r_57 * l_r_89;
    double l_r_107 = l_r_43 * l_r_76 + l_r_50 * l_r_82 + l_r_57 * l_r_62;
    double l_r_108 = l_r_43 * l_r_77 + l_r_50 * l_r_83 + l_r_57 * l_r_67;
    vec3 v_109 = vcons3(l_r_42, l_r_49, l_r_56);
    double l_r_110 = 0.e0 * (l_r_40 * l_r_62 + l_r_47 * l_r_79 + l_r_54 * l_r_85);
    double l_r_111 = 0.e0 * l_r_92;
    double l_r_112 = 0.e0 * l_r_97;
    double l_r_113 = 0.e0 * (l_r_42 * l_r_67 + l_r_49 * l_r_80 + l_r_56 * l_r_86);
    double l_r_114 = 0.e0 * l_r_103;
    double l_r_115 = 0.e0 * (l_r_43 * l_r_72 + l_r_50 * l_r_81 + l_r_57 * l_r_87);
    double l_r_116 = l_r_110 + 0.e0 * l_r_91;
    double l_r_117 = 0.e0 * (l_r_40 * l_r_73 + l_r_47 * l_r_62 + l_r_54 * l_r_88);
    double l_r_118 = 0.e0 * l_r_94;
    double l_r_119 = 0.e0 * l_r_99;
    double l_r_120 = 0.e0 * (l_r_42 * l_r_74 + l_r_49 * l_r_67 + l_r_56 * l_r_89);
    double l_r_121 = 0.e0 * l_r_105;
    double l_r_122 = 0.e0 * (l_r_43 * l_r_75 + l_r_50 * l_r_72 + l_r_57 * l_r_90);
    double l_r_123 = l_r_117 + 0.e0 * l_r_93;
    double l_r_124 = 0.e0 * (l_r_40 * l_r_76 + l_r_47 * l_r_82 + l_r_54 * l_r_62);
    double l_r_125 = 0.e0 * l_r_96;
    double l_r_126 = 0.e0 * l_r_101;
    double l_r_127 = 0.e0 * (l_r_42 * l_r_77 + l_r_49 * l_r_83 + l_r_56 * l_r_67);
    double l_r_128 = 0.e0 * l_r_107;
    double l_r_129 = 0.e0 * (l_r_43 * l_r_78 + l_r_50 * l_r_84 + l_r_57 * l_r_72);
    double l_r_130 = l_r_124 + 0.e0 * l_r_95;
    double l_r_131 = 0.e0 * l_r_98;
    double l_r_132 = 0.e0 * l_r_104;
    double l_r_133 = 0.e0 * l_r_100;
    double l_r_134 = 0.e0 * l_r_106;
    double l_r_135 = 0.e0 * l_r_102;
    double l_r_136 = 0.e0 * l_r_108;
    double l_op1_e3_l_28_137 = 0.2e1 * vdot3(vcons3(l_r_40, l_r_47, l_r_54),
        vcons3(vdot3(v_109, vcons3(l_r_72, l_r_81, l_r_87)), vdot3(v_109, vcons3(l_r_75, l_r_72, l_r_90)),
            vdot3(v_109, vcons3(l_r_78, l_r_84, l_r_72))));
    double l_prod_138 = v_10[0] * l_prod_33;
    double l_prod_139 = 0.1e1 * (v_10[1] * 0.1e1);
    double l_prod_140 = 0.1e1 * (0.1e1 * v_10[2]);
    double l_sum_141 = l_basisEval_36 + (-0.1e1 * l_prod_140 + (-0.1e1 * l_prod_139 + -0.1e1 * l_prod_138));
    double l_basisEval_142 = 0.1e1 * l_prod_138;
    double l_basisEval_143 = 0.1e1 * l_prod_139;
    double l_basisEval_144 = 0.1e1 * l_prod_140;
    vec3 v_145 = vcons3(
        (l_r_116 + l_r_111 + l_r_112 + l_r_113 + 0.1e1 * l_r_98 + l_r_114 + -0.1e1 * l_r_104 + l_r_115) / l_op1_e3_l_28_137,
        (l_r_123 + l_r_118 + l_r_119 + l_r_120 + 0.1e1 * l_r_100 + l_r_121 + -0.1e1 * l_r_106 + l_r_122) / l_op1_e3_l_28_137,
        (l_r_130 + l_r_125 + l_r_126 + l_r_127 + 0.1e1 * l_r_102 + l_r_128 + -0.1e1 * l_r_108 + l_r_129) / l_op1_e3_l_28_137);
    vec3 v_146 = vcons3(
        (l_r_116 + -0.1e1 * l_r_92 + l_r_112 + l_r_113 + l_r_131 + 0.1e1 * l_r_103 + l_r_132 + l_r_115) / l_op1_e3_l_28_137,
        (l_r_123 + -0.1e1 * l_r_94 + l_r_119 + l_r_120 + l_r_133 + 0.1e1 * l_r_105 + l_r_134 + l_r_122) / l_op1_e3_l_28_137,
        (l_r_130 + -0.1e1 * l_r_96 + l_r_126 + l_r_127 + l_r_135 + 0.1e1 * l_r_107 + l_r_136 + l_r_129) / l_op1_e3_l_28_137);
    vec3 v_147 = vcons3(
        (l_r_110 + 0.1e1 * l_r_91 + l_r_111 + -0.1e1 * l_r_97 + l_r_113 + l_r_131 + l_r_114 + l_r_132 + l_r_115) / l_op1_e3_l_28_137,
        (l_r_117 + 0.1e1 * l_r_93 + l_r_118 + -0.1e1 * l_r_99 + l_r_120 + l_r_133 + l_r_121 + l_r_134 + l_r_122) / l_op1_e3_l_28_137,
        (l_r_124 + 0.1e1 * l_r_95 + l_r_125 + -0.1e1 * l_r_101 + l_r_127 + l_r_135 + l_r_128 + l_r_136 + l_r_129) / l_op1_e3_l_28_137);
    vec3 v_148 = v_6 - vcons3(
        l_dof_load_15 * l_sum_141 + l_dof_load_20 * l_basisEval_142 + l_dof_load_25 * l_basisEval_143 + l_dof_load_30 * l_basisEval_144,
        l_dof_load_16 * l_sum_141 + l_dof_load_21 * l_basisEval_142 + l_dof_load_26 * l_basisEval_143 + l_dof_load_31 * l_basisEval_144,
        l_dof_load_17 * l_sum_141 + l_dof_load_22 * l_basisEval_142 + l_dof_load_27 * l_basisEval_143 + l_dof_load_32 * l_basisEval_144);
    vec3 v_149 = vcons3(vdot3(v_145, v_148), vdot3(v_146, v_148), vdot3(v_147, v_148));
    vec3 v_150 = vcons3(vdot3(v_145, v_9), vdot3(v_146, v_9), vdot3(v_147, v_9));
    double l_op1_e3_l_49_151 = vdot3(v_11, v_150);
    double l__t_152 = (0.57735026919e0 - vdot3(v_11, v_149)) / l_op1_e3_l_49_151;
    vec3 v_153 = v_11;
    vec3 v_154 = v_145;
    vec3 v_155 = v_146;
    vec3 v_156 = v_147;
    vec3 v_157 = v_149;
    vec3 v_158 = v_150;
    vec3 v_159 = v_6;
    vec3 v_160 = v_8;
    vec3 v_161 = v_9;
    if (l__t_152 > -0.e0 && HUGE_VAL > l__t_152) {
        int32_t l_face_164;
        double l__t_165;
        vec3 v_162 = vscale3(l__t_152, v_158);
        vec3 v_163 = vcons3(0.1e-8, 0.1e-8, 0.1e-8) + v_162 + v_157;
        if (0.1e1 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_162 + v_157) && (v_163[0] > -0.e0 && (v_163[1] > -0.e0 && v_163[2] > -0.e0))) {
            l_face_164 = 0;
            l__t_165 = l__t_152;
        }
        else {
            l_face_164 = -1;
            l__t_165 = HUGE_VAL;
        }
        l_face_166 = l_face_164;
        l_time_167 = l__t_165;
    }
    else {
        l_face_166 = -1;
        l_time_167 = HUGE_VAL;
    }
    vec3 v_168 = vcons3(0.1e1, 0.e0, 0.e0);
    double l_op1_e3_l_49_169 = vdot3(v_168, v_158);
    double l__t_170 = (0.e0 - vdot3(v_168, v_157)) / l_op1_e3_l_49_169;
    vec3 v_171 = v_168;
    if (l__t_170 > -0.e0 && l_time_167 > l__t_170) {
        int32_t l_ilit_174;
        double l__t_175;
        vec3 v_172 = vscale3(l__t_170, v_158);
        vec3 v_173 = vcons3(0.1e-8, 0.1e-8, 0.1e-8) + v_172 + v_157;
        if (0.1e1 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_172 + v_157) && (v_173[0] > -0.e0 && (v_173[1] > -0.e0 && v_173[2] > -0.e0))) {
            l_ilit_174 = 1;
            l__t_175 = l__t_170;
        }
        else {
            l_ilit_174 = l_face_166;
            l__t_175 = l_time_167;
        }
        l_face_176 = l_ilit_174;
        l_time_177 = l__t_175;
    }
    else {
        l_face_176 = l_face_166;
        l_time_177 = l_time_167;
    }
    double l_op1_e3_l_49_178 = vdot3(v_161, v_158);
    double l__t_179 = (0.e0 - vdot3(v_161, v_157)) / l_op1_e3_l_49_178;
    if (l__t_179 > -0.e0 && l_time_177 > l__t_179) {
        int32_t l_ilit_182;
        double l__t_183;
        vec3 v_180 = vscale3(l__t_179, v_158);
        vec3 v_181 = vcons3(0.1e-8, 0.1e-8, 0.1e-8) + v_180 + v_157;
        if (0.1e1 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_180 + v_157) && (v_181[0] > -0.e0 && (v_181[1] > -0.e0 && v_181[2] > -0.e0))) {
            l_ilit_182 = 2;
            l__t_183 = l__t_179;
        }
        else {
            l_ilit_182 = l_face_176;
            l__t_183 = l_time_177;
        }
        l_face_184 = l_ilit_182;
        l_time_185 = l__t_183;
    }
    else {
        l_face_184 = l_face_176;
        l_time_185 = l_time_177;
    }
    vec3 v_186 = vcons3(0.e0, 0.e0, 0.1e1);
    double l_op1_e3_l_49_187 = vdot3(v_186, v_158);
    double l__t_188 = (0.e0 - vdot3(v_186, v_157)) / l_op1_e3_l_49_187;
    vec3 v_189 = v_186;
    if (l__t_188 > -0.e0 && l_time_185 > l__t_188) {
        int32_t l_ilit_192;
        double l__t_193;
        vec3 v_190 = vscale3(l__t_188, v_158);
        vec3 v_191 = vcons3(0.1e-8, 0.1e-8, 0.1e-8) + v_190 + v_157;
        if (0.1e1 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_190 + v_157) && (v_191[0] > -0.e0 && (v_191[1] > -0.e0 && v_191[2] > -0.e0))) {
            l_ilit_192 = 3;
            l__t_193 = l__t_188;
        }
        else {
            l_ilit_192 = l_face_184;
            l__t_193 = l_time_185;
        }
        l_face_194 = l_ilit_192;
        l_time_195 = l__t_193;
    }
    else {
        l_face_194 = l_face_184;
        l_time_195 = l_time_185;
    }
    if (l_face_194 != -1) {
        v_196 = vcons2(l_time_195, static_cast<double>(l_face_194));
    }
    else {
        v_196 = vcons2(-0.1e1, -0.1e1);
    }
    double l_t5_197 = v_196[0];
    double l_prod_198 = v_160[0] * l_prod_33;
    double l_prod_199 = 0.1e1 * (v_160[1] * 0.1e1);
    double l_prod_200 = 0.1e1 * (0.1e1 * v_160[2]);
    double l_sum_201 = l_basisEval_36 + (-0.1e1 * l_prod_200 + (-0.1e1 * l_prod_199 + -0.1e1 * l_prod_198));
    double l_basisEval_202 = 0.1e1 * l_prod_198;
    double l_basisEval_203 = 0.1e1 * l_prod_199;
    double l_basisEval_204 = 0.1e1 * l_prod_200;
    vec3 v_205 = v_159 - vcons3(
        l_dof_load_15 * l_sum_201 + l_dof_load_20 * l_basisEval_202 + l_dof_load_25 * l_basisEval_203 + l_dof_load_30 * l_basisEval_204,
        l_dof_load_16 * l_sum_201 + l_dof_load_21 * l_basisEval_202 + l_dof_load_26 * l_basisEval_203 + l_dof_load_31 * l_basisEval_204,
        l_dof_load_17 * l_sum_201 + l_dof_load_22 * l_basisEval_202 + l_dof_load_27 * l_basisEval_203 + l_dof_load_32 * l_basisEval_204);
    vec3 v_206 = vcons3(vdot3(v_154, v_205), vdot3(v_155, v_205), vdot3(v_156, v_205));
    double l__t_207 = (0.57735026919e0 - vdot3(v_153, v_206)) / l_op1_e3_l_49_151;
    vec3 v_208 = v_206;
    if (l__t_207 > -0.e0 && HUGE_VAL > l__t_207) {
        int32_t l_face_211;
        double l__t_212;
        vec3 v_209 = vscale3(l__t_207, v_158);
        vec3 v_210 = vcons3(0.1e-8, 0.1e-8, 0.1e-8) + v_209 + v_208;
        if (0.1e1 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_209 + v_208) && (v_210[0] > -0.e0 && (v_210[1] > -0.e0 && v_210[2] > -0.e0))) {
            l_face_211 = 0;
            l__t_212 = l__t_207;
        }
        else {
            l_face_211 = -1;
            l__t_212 = HUGE_VAL;
        }
        l_face_213 = l_face_211;
        l_time_214 = l__t_212;
    }
    else {
        l_face_213 = -1;
        l_time_214 = HUGE_VAL;
    }
    double l__t_215 = (0.e0 - vdot3(v_171, v_208)) / l_op1_e3_l_49_169;
    if (l__t_215 > -0.e0 && l_time_214 > l__t_215) {
        int32_t l_ilit_218;
        double l__t_219;
        vec3 v_216 = vscale3(l__t_215, v_158);
        vec3 v_217 = vcons3(0.1e-8, 0.1e-8, 0.1e-8) + v_216 + v_208;
        if (0.1e1 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_216 + v_208) && (v_217[0] > -0.e0 && (v_217[1] > -0.e0 && v_217[2] > -0.e0))) {
            l_ilit_218 = 1;
            l__t_219 = l__t_215;
        }
        else {
            l_ilit_218 = l_face_213;
            l__t_219 = l_time_214;
        }
        l_face_220 = l_ilit_218;
        l_time_221 = l__t_219;
    }
    else {
        l_face_220 = l_face_213;
        l_time_221 = l_time_214;
    }
    double l__t_222 = (0.e0 - vdot3(v_161, v_208)) / l_op1_e3_l_49_178;
    if (l__t_222 > -0.e0 && l_time_221 > l__t_222) {
        int32_t l_ilit_225;
        double l__t_226;
        vec3 v_223 = vscale3(l__t_222, v_158);
        vec3 v_224 = vcons3(0.1e-8, 0.1e-8, 0.1e-8) + v_223 + v_208;
        if (0.1e1 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_223 + v_208) && (v_224[0] > -0.e0 && (v_224[1] > -0.e0 && v_224[2] > -0.e0))) {
            l_ilit_225 = 2;
            l__t_226 = l__t_222;
        }
        else {
            l_ilit_225 = l_face_220;
            l__t_226 = l_time_221;
        }
        l_face_227 = l_ilit_225;
        l_time_228 = l__t_226;
    }
    else {
        l_face_227 = l_face_220;
        l_time_228 = l_time_221;
    }
    double l__t_229 = (0.e0 - vdot3(v_189, v_208)) / l_op1_e3_l_49_187;
    if (l__t_229 > -0.e0 && l_time_228 > l__t_229) {
        int32_t l_ilit_232;
        double l__t_233;
        vec3 v_230 = vscale3(l__t_229, v_158);
        vec3 v_231 = vcons3(0.1e-8, 0.1e-8, 0.1e-8) + v_230 + v_208;
        if (0.1e1 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_230 + v_208) && (v_231[0] > -0.e0 && (v_231[1] > -0.e0 && v_231[2] > -0.e0))) {
            l_ilit_232 = 3;
            l__t_233 = l__t_229;
        }
        else {
            l_ilit_232 = l_face_227;
            l__t_233 = l_time_228;
        }
        l_face_234 = l_ilit_232;
        l_time_235 = l__t_233;
    }
    else {
        l_face_234 = l_face_227;
        l_time_235 = l_time_228;
    }
    if (l_face_234 != -1) {
        v_236 = vcons2(l_time_235, static_cast<double>(l_face_234));
    }
    else {
        v_236 = vcons2(-0.1e1, -0.1e1);
    }
    double l_t6_237 = v_236[0];
    wrld->print() << "\n" << std::flush;
    tensor_3 _arg_238;
    vpack3(_arg_238, v_159 + vscale3(l_t5_197, v_161));
    wrld->print() << "We have time " << l_t5_197 << " and result" << tensor_ref_3(_arg_238) << "\n" << std::flush;
    tensor_3 _arg_239;
    vpack3(_arg_239, v_159 + vscale3(l_t6_237, v_161));
    wrld->print() << "We have time " << l_t6_237 << " and result" << tensor_ref_3(_arg_239) << "\n" << std::flush;
    wrld->print() << "\n" << std::flush;
    return diderot::kStabilize;
}
extern "C" bool justTypes_output_get_pos (justTypes_world_t *cWrld, Nrrd *nData)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    // Compute sizes of nrrd file
    size_t sizes[2];
    sizes[0] = 3;
    sizes[1] = wrld->_strands.num_stable();
    // Allocate nData nrrd
    if (nrrdMaybeAlloc_nva(nData, nrrdTypeDouble, 2, sizes) != 0) {
        char *msg = biffGetDone(NRRD);
        biffMsgAdd(wrld->_errors, msg);
        std::free(msg);
        return true;
    }
    // copy data to output nrrd
    char *cp = reinterpret_cast<char *>(nData->data);
    for (auto ix = wrld->_strands.begin_stable(); ix != wrld->_strands.end_stable(); ix = wrld->_strands.next_stable(
        ix)) {
        memcpy(cp, &wrld->_strands.strand(ix)->sv_pos, 3 * sizeof(double));
        cp += 3 * sizeof(double);
    }
    nData->axis[0].kind = nrrdKind3Vector;
    nData->axis[1].kind = nrrdKindList;
    return false;
}
extern "C" bool justTypes_output_get__pos (justTypes_world_t *cWrld, Nrrd *nData)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    // Compute sizes of nrrd file
    size_t sizes[2];
    sizes[0] = 3;
    sizes[1] = wrld->_strands.num_stable();
    // Allocate nData nrrd
    if (nrrdMaybeAlloc_nva(nData, nrrdTypeDouble, 2, sizes) != 0) {
        char *msg = biffGetDone(NRRD);
        biffMsgAdd(wrld->_errors, msg);
        std::free(msg);
        return true;
    }
    // copy data to output nrrd
    char *cp = reinterpret_cast<char *>(nData->data);
    for (auto ix = wrld->_strands.begin_stable(); ix != wrld->_strands.end_stable(); ix = wrld->_strands.next_stable(
        ix)) {
        memcpy(cp, &wrld->_strands.strand(ix)->sv__pos, 3 * sizeof(double));
        cp += 3 * sizeof(double);
    }
    nData->axis[0].kind = nrrdKind3Vector;
    nData->axis[1].kind = nrrdKindList;
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
    this->_tree = new diderot::kdtree<0, double, strand_array> (&this->_strands);
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
    diderot::dynseq< mesh_cell_msh > seq_1 = glob->gv_0cell_a;
    int32_t base[1] = {0,};
    uint32_t size[1] = {static_cast<uint32_t>(seq_1.length()),};
    if (this->alloc(base, size)) {
        return true;
    }
    uint32_t ix = 0;
    for (auto it_2 = seq_1.cbegin(); it_2 != seq_1.cend(); ++it_2) {
        auto i_k_240 = *it_2;
        gg_init(this->_strands.strand(ix), i_k_240);
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

} // namespace justTypes
/*---------- end namespace-close.in ----------*/

/*---------- begin c-wrappers.in ----------*/
extern "C" uint32_t justTypes_num_strands (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
    return w->_strands.num_alive();
}

extern "C" uint32_t justTypes_num_active_strands (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
    return w->_strands.num_active();
}

extern "C" uint32_t justTypes_num_stable_strands (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
    return w->_strands.num_stable();
}

extern "C" bool justTypes_any_errors (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
    return (w->_errors->errNum > 0);
}

extern "C" char *justTypes_get_errors (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
    char *msg = biffMsgStrGet (w->_errors);
    biffMsgClear (w->_errors);
    return msg;
}

extern "C" justTypes_world_t *justTypes_new_world ()
{
    justTypes::world *w = new (std::nothrow) justTypes::world();
    return reinterpret_cast<justTypes_world_t *>(w);
}

extern "C" bool justTypes_init_world (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);

    if (w->_stage != diderot::POST_NEW) {
        w->error ("multiple calls to justTypes_init_world");
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

extern "C" bool justTypes_create_strands (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);

    if (w->_stage < diderot::POST_INIT) {
        w->error ("must call justTypes_init_world before justTypes_create_strands");
        return true;
    }
    else if (w->_stage > diderot::POST_INIT) {
        w->error ("multiple calls to justTypes_create_strands");
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

extern "C" uint32_t justTypes_run (justTypes_world_t *wrld, uint32_t maxNSteps)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);

    if (w->_stage < diderot::POST_CREATE) {
        w->error ("attempt to run uninitialized program");
        return 0;
    }
    else if (w->_stage == diderot::DONE) {
        return 0;
    }

    return w->run(maxNSteps);
}

extern "C" void justTypes_shutdown (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
    delete w;
}

extern "C" void justTypes_set_verbose (justTypes_world_t *wrld, bool mode)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
    w->_verbose = (mode ? true : false);
}

extern "C" bool justTypes_get_verbose (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
    return static_cast<bool>(w->_verbose);
}

extern "C" bool justTypes_set_printer_cb (justTypes_world_t *wrld, bool (*pr)(void *, char *), void *data)
{
  /* FIXME: implement printer callback */
    return true;
}

#ifdef DIDEROT_TARGET_PARALLEL

extern "C" uint32_t justTypes_get_num_cores (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
    return w->_sched->_numHWCores;
}

extern "C" bool justTypes_set_num_workers (justTypes_world_t *wrld, uint32_t nw)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
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

extern "C" uint32_t justTypes_get_num_workers (justTypes_world_t *wrld)
{
    justTypes::world *w = reinterpret_cast<justTypes::world *>(wrld);
    return w->_sched->_numWorkers;
}

#endif /* DIDEROT_TARGET_PARALLEL */
/*---------- end c-wrappers.in ----------*/


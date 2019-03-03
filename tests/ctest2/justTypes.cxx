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

#define DIDEROT_HAS_STABILIZE_METHOD
#define DIDEROT_STRAND_HAS_CONSTR
#define DIDEROT_STRAND_ARRAY
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
    typedef double vec2 __attribute__ ((vector_size (16)));
    typedef double vec3 __attribute__ ((vector_size (32)));
    typedef double vec4 __attribute__ ((vector_size (32)));
    struct tensor_ref_2 : public diderot::tensor_ref<double,2> {
        tensor_ref_2 (const double *src);
        tensor_ref_2 (struct tensor_2 const & ten);
        tensor_ref_2 (tensor_ref_2 const & ten);
    };
    struct tensor_ref_4 : public diderot::tensor_ref<double,4> {
        tensor_ref_4 (const double *src);
        tensor_ref_4 (struct tensor_4 const & ten);
        tensor_ref_4 (tensor_ref_4 const & ten);
    };
    struct tensor_ref_3 : public diderot::tensor_ref<double,3> {
        tensor_ref_3 (const double *src);
        tensor_ref_3 (struct tensor_3 const & ten);
        tensor_ref_3 (tensor_ref_3 const & ten);
    };
    struct tensor_ref_7 : public diderot::tensor_ref<double,7> {
        tensor_ref_7 (const double *src);
        tensor_ref_7 (struct tensor_7 const & ten);
        tensor_ref_7 (tensor_ref_7 const & ten);
    };
    struct tensor_ref_3_2 : public diderot::tensor_ref<double,6> {
        tensor_ref_3_2 (const double *src);
        tensor_ref_3_2 (struct tensor_3_2 const & ten);
        tensor_ref_3_2 (tensor_ref_3_2 const & ten);
        tensor_ref_2 last (uint32_t i)
        {
            return &this->_data[i];
        }
    };
    struct tensor_2 : public diderot::tensor<double,2> {
        tensor_2 ()
            : diderot::tensor<double,2>()
        { }
        tensor_2 (std::initializer_list< double > const & il)
            : diderot::tensor<double,2>(il)
        { }
        tensor_2 (const double *src)
            : diderot::tensor<double,2>(src)
        { }
        tensor_2 (tensor_2 const & ten)
            : diderot::tensor<double,2>(ten._data)
        { }
        ~tensor_2 () { }
        tensor_2 & operator= (tensor_2 const & src);
        tensor_2 & operator= (tensor_ref_2 const & src);
        tensor_2 & operator= (std::initializer_list< double > const & il);
        tensor_2 & operator= (const double *src);
    };
    struct tensor_7 : public diderot::tensor<double,7> {
        tensor_7 ()
            : diderot::tensor<double,7>()
        { }
        tensor_7 (std::initializer_list< double > const & il)
            : diderot::tensor<double,7>(il)
        { }
        tensor_7 (const double *src)
            : diderot::tensor<double,7>(src)
        { }
        tensor_7 (tensor_7 const & ten)
            : diderot::tensor<double,7>(ten._data)
        { }
        ~tensor_7 () { }
        tensor_7 & operator= (tensor_7 const & src);
        tensor_7 & operator= (tensor_ref_7 const & src);
        tensor_7 & operator= (std::initializer_list< double > const & il);
        tensor_7 & operator= (const double *src);
    };
    struct tensor_4 : public diderot::tensor<double,4> {
        tensor_4 ()
            : diderot::tensor<double,4>()
        { }
        tensor_4 (std::initializer_list< double > const & il)
            : diderot::tensor<double,4>(il)
        { }
        tensor_4 (const double *src)
            : diderot::tensor<double,4>(src)
        { }
        tensor_4 (tensor_4 const & ten)
            : diderot::tensor<double,4>(ten._data)
        { }
        ~tensor_4 () { }
        tensor_4 & operator= (tensor_4 const & src);
        tensor_4 & operator= (tensor_ref_4 const & src);
        tensor_4 & operator= (std::initializer_list< double > const & il);
        tensor_4 & operator= (const double *src);
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
    struct tensor_3_2 : public diderot::tensor<double,6> {
        tensor_3_2 ()
            : diderot::tensor<double,6>()
        { }
        tensor_3_2 (std::initializer_list< double > const & il)
            : diderot::tensor<double,6>(il)
        { }
        tensor_3_2 (const double *src)
            : diderot::tensor<double,6>(src)
        { }
        tensor_3_2 (tensor_3_2 const & ten)
            : diderot::tensor<double,6>(ten._data)
        { }
        ~tensor_3_2 () { }
        tensor_3_2 & operator= (tensor_3_2 const & src);
        tensor_3_2 & operator= (tensor_ref_3_2 const & src);
        tensor_3_2 & operator= (std::initializer_list< double > const & il);
        tensor_3_2 & operator= (const double *src);
        tensor_ref_2 last (uint32_t i)
        {
            return &this->_data[i];
        }
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
    struct mesh_cell_msh {
        int32_t cell;
        msh mesh;
        std::ostream & operator<< (std::ostream & os)
        {
            return os << this->cell;
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
    struct mesh_pos_msh {
        msh mesh;
        int32_t cell;
        tensor_3 refPos;
        tensor_3 worldPos;
        bool wpc;
        bool valid;
        mesh_pos_msh (msh mesh, int32_t cell, tensor_ref_3 refPos, tensor_ref_3 worldPos, bool wpc, bool valid)
            : mesh(mesh), cell(cell), wpc(wpc), valid(valid)
        {
            mesh = mesh;
            cell = cell;
            refPos = refPos;
            worldPos = worldPos;
            wpc = wpc;
            valid = valid;
        }
        mesh_pos_msh ()
            : cell(-1), wpc(false), valid(false)
        { }
        mesh_pos_msh (msh mesh)
            : cell(-1), wpc(false), valid(false), mesh(mesh)
        { }
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
    inline tensor_ref_2::tensor_ref_2 (const double *src)
        : diderot::tensor_ref<double,2>(src)
    { }
    inline tensor_ref_2::tensor_ref_2 (struct tensor_2 const & ten)
        : diderot::tensor_ref<double,2>(ten._data)
    { }
    inline tensor_ref_2::tensor_ref_2 (tensor_ref_2 const & ten)
        : diderot::tensor_ref<double,2>(ten._data)
    { }
    inline tensor_ref_4::tensor_ref_4 (const double *src)
        : diderot::tensor_ref<double,4>(src)
    { }
    inline tensor_ref_4::tensor_ref_4 (struct tensor_4 const & ten)
        : diderot::tensor_ref<double,4>(ten._data)
    { }
    inline tensor_ref_4::tensor_ref_4 (tensor_ref_4 const & ten)
        : diderot::tensor_ref<double,4>(ten._data)
    { }
    inline tensor_ref_3::tensor_ref_3 (const double *src)
        : diderot::tensor_ref<double,3>(src)
    { }
    inline tensor_ref_3::tensor_ref_3 (struct tensor_3 const & ten)
        : diderot::tensor_ref<double,3>(ten._data)
    { }
    inline tensor_ref_3::tensor_ref_3 (tensor_ref_3 const & ten)
        : diderot::tensor_ref<double,3>(ten._data)
    { }
    inline tensor_ref_7::tensor_ref_7 (const double *src)
        : diderot::tensor_ref<double,7>(src)
    { }
    inline tensor_ref_7::tensor_ref_7 (struct tensor_7 const & ten)
        : diderot::tensor_ref<double,7>(ten._data)
    { }
    inline tensor_ref_7::tensor_ref_7 (tensor_ref_7 const & ten)
        : diderot::tensor_ref<double,7>(ten._data)
    { }
    inline tensor_ref_3_2::tensor_ref_3_2 (const double *src)
        : diderot::tensor_ref<double,6>(src)
    { }
    inline tensor_ref_3_2::tensor_ref_3_2 (struct tensor_3_2 const & ten)
        : diderot::tensor_ref<double,6>(ten._data)
    { }
    inline tensor_ref_3_2::tensor_ref_3_2 (tensor_ref_3_2 const & ten)
        : diderot::tensor_ref<double,6>(ten._data)
    { }
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
    inline tensor_2 & tensor_2::operator= (std::initializer_list< double > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_2 & tensor_2::operator= (const double *src)
    {
        this->copy(src);
        return *this;
    }
    inline tensor_7 & tensor_7::operator= (tensor_7 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_7 & tensor_7::operator= (tensor_ref_7 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_7 & tensor_7::operator= (std::initializer_list< double > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_7 & tensor_7::operator= (const double *src)
    {
        this->copy(src);
        return *this;
    }
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
    inline tensor_4 & tensor_4::operator= (std::initializer_list< double > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_4 & tensor_4::operator= (const double *src)
    {
        this->copy(src);
        return *this;
    }
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
    inline tensor_3_2 & tensor_3_2::operator= (tensor_3_2 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_3_2 & tensor_3_2::operator= (tensor_ref_3_2 const & src)
    {
        this->copy(src._data);
        return *this;
    }
    inline tensor_3_2 & tensor_3_2::operator= (std::initializer_list< double > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_3_2 & tensor_3_2::operator= (const double *src)
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
    std::ostream & operator<< (std::ostream & os, const mesh_pos_msh &  pos)
    {
        return os << pos.cell;
    }
    mesh_pos_msh allBuild (msh mesh, int32_t cell, tensor_ref_3 refPos, tensor_ref_3 worldPos, bool wpc, bool valid)
    {
        mesh_pos_msh pos;
        pos.mesh = mesh;
        pos.cell = cell;
        pos.refPos = refPos._data;
        pos.worldPos = worldPos._data;
        pos.wpc = wpc;
        pos.valid = valid;
        return pos;
    }
    mesh_pos_msh invalidBuild (msh mesh)
    {
        mesh_pos_msh pos;
        pos.mesh = mesh;
        return pos;
    }
    mesh_pos_msh invalidBuild (msh mesh, tensor_ref_3 refPos)
    {
        mesh_pos_msh pos;
        pos.mesh = mesh;
        pos.refPos = refPos;
        return pos;
    }
    mesh_pos_msh refBuild (msh mesh, int32_t cell, tensor_ref_3 refPos)
    {
        mesh_pos_msh pos;
        pos.mesh = mesh;
        pos.cell = cell;
        pos.refPos = refPos;
        pos.wpc = false;
        pos.valid = true;
        return pos;
    }
    char * copy_to  (const mesh_pos_msh dumb, char *cp)
    {
        size_t nbytes1 = sizeof(tensor_3);
        size_t nbytes2 = sizeof(int32_t);
        (std::memcpy)(cp, &(dumb.cell), nbytes2);
        cp += nbytes2;
        if (dumb.valid) {
            (std::memcpy)(cp, dumb.refPos._data, nbytes1);
            cp += nbytes1;
            (std::memcpy)(cp, dumb.worldPos._data, nbytes1);
            cp += nbytes1;
        }
        else {
            const double dumb[6] = {};
            size_t nbytes3 = sizeof(double[6]);
            (std::memcpy)(cp, &dumb, nbytes3);
            cp += nbytes3;
        }
        return cp;
    }
    mesh_pos_msh makeFem (msh mesh, tensor_ref_7 data)
    {
        double cell = data[0];
        tensor_ref_3 ref = &data[1];
        tensor_ref_3 world = &data[4];
        int32_t cellInt = static_cast<int32_t>(cell);
        bool test = 0 < cellInt;
        return allBuild(mesh, cellInt, ref, world, test, test);
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
struct raycast_strand;
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
    bool gv_0b043B_intermedateGlobal;
    bool gv_0c043D_intermedateGlobal;
    bool gv_isoval;
    bool gv_thick;
    bool gv_camEye;
    bool gv_camAt;
    bool gv_camUp;
    bool gv_camFOV;
    bool gv_iresU;
    bool gv_iresV;
    bool gv_camNear;
    bool gv_camFar;
    bool gv_refStep;
    bool gv_rayStep;
    bool gv_lightVsp;
    bool gv_phongKa;
    bool gv_phongKd;
    bool gv_debug;
    bool gv_su;
    bool gv_sv;
} defined_inputs;
struct globals {
    msh gv_a;
    fns gv_0b043B_intermedateGlobal;
    FUNC gv_0c043D_intermedateGlobal;
    double gv_isoval;
    double gv_thick;
    tensor_3 gv_camEye;
    tensor_3 gv_camAt;
    tensor_3 gv_camUp;
    double gv_camFOV;
    int32_t gv_iresU;
    int32_t gv_iresV;
    double gv_camNear;
    double gv_camFar;
    double gv_refStep;
    double gv_rayStep;
    tensor_3 gv_lightVsp;
    double gv_phongKa;
    double gv_phongKd;
    bool gv_debug;
    int32_t gv_su;
    int32_t gv_sv;
    double gv_camDist;
    tensor_3 gv_camN;
    tensor_3 gv_camU;
    tensor_3 gv_camV;
    double gv_camVmax;
    double gv_camUmax;
    tensor_3 gv_light;
    msh gv__t;
    fns gv__tX;
    FUNC gv_c;
    diderot::image1d< double, double, 3 > gv_I;
    ~globals ()
    {
        this->gv_I.unregister_global();
    }
};
struct raycast_strand {
    double sv_rayU;
    double sv_rayV;
    double sv_rayN;
    tensor_3 sv_rayVec;
    double sv_transp;
    tensor_3 sv_rgb;
    tensor_4 sv_rgba;
    double sv_gray;
    int32_t sv_ui;
    int32_t sv_vi;
};
/*---------- begin seq-sarr.in ----------*/
// forward declarations of strand methods
#ifdef DIDEROT_HAS_START_METHOD
static diderot::strand_status raycast_start (raycast_strand *self);
#endif // DIDEROT_HAS_START_METHOD
static diderot::strand_status raycast_update (world *wrld, globals *glob, raycast_strand *self);
#ifdef DIDEROT_HAS_STABILIZE_METHOD
static void raycast_stabilize (world *wrld, globals *glob, raycast_strand *self);
#endif // DIDEROT_HAS_STABILIZE_METHOD

// if we have both communication and "die", then we need to track when strands die
// so that we can rebuild the list of strands use to construct the kd-tree
#if defined(DIDEROT_HAS_STRAND_COMMUNICATION) && !defined(DIDEROT_HAS_STRAND_DIE)
#  define TRACK_STRAND_DEATH
#endif

// strand_array for SEQUENTIAL/NO BSP/SINGLE STATE/DIRECT ACCESS
//
struct strand_array {
    typedef raycast_strand strand_t;
    typedef uint32_t index_t;
    typedef index_t sid_t;              // strand ID (index into strand-state storage)

    uint8_t             *_status;       // the array of status information for the strands
    char                *_storage;      // points to array of raycast_strand structs
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
    raycast_strand *id_to_strand (sid_t id) const
    {
        assert (id < this->_nItems);
        return reinterpret_cast<raycast_strand *>(this->_storage + id * sizeof(raycast_strand));
    }

  // return a strand's status
    diderot::strand_status status (index_t ix) const
    {
        assert (ix < this->_nItems);
        return static_cast<diderot::strand_status>(this->_status[ix]);
    }
  // return a pointer to the given strand
    raycast_strand *strand (index_t ix) const
    {
        return this->id_to_strand(this->id(ix));
    }
  // return a pointer to the local state of strand ix
    raycast_strand *local_state (index_t ix) const
    {
        return this->strand(ix);
    }
  // return a pointer to the local state of strand with the given ID
    raycast_strand *id_to_local_state (sid_t id) const
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
        this->_storage = static_cast<char *>(std::malloc (nItems * sizeof(raycast_strand)));
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
            new (this->strand(ix)) raycast_strand;
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
        return raycast_start(this->strand(ix));
    }
#endif // DIDEROT_HAS_START_METHOD

  // invoke strand's update method
    diderot::strand_status strand_update (world *wrld, globals *glob, index_t ix)
    {
        return raycast_update(wrld, glob, this->strand(ix));
    }

  // invoke strand's stabilize method
    index_t strand_stabilize (world *wrld, globals *glob, index_t ix)
    {
#ifdef DIDEROT_HAS_STABILIZE_METHOD
        raycast_stabilize (wrld, glob, this->strand(ix));
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
        this->strand(ix)->~raycast_strand();
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
    bool alloc (int32_t base[2], uint32_t size[2]);
    bool create_strands ();
    uint32_t run (uint32_t max_nsteps);
    void swap_state ();
};
// ***** Begin synthesized operations *****

inline double clamp (double lo, double hi, double x)
{
    return x <= lo ? lo : hi <= x ? hi : x;
}
inline vec2 vload2 (const double *vp)
{
    return __extension__ (vec2){vp[0], vp[1]};
}
inline vec2 vcons2 (double r0, double r1)
{
    return __extension__ (vec2){r0, r1};
}
static std::ostream& operator<< (std::ostream& outs, tensor_ref_3 const & ten)
{
    return outs << "[" << ten._data[0] << "," << ten._data[1] << "," << ten._data[2] << "]";
}
inline vec3 vload3 (const double *vp)
{
    return __extension__ (vec3){vp[0], vp[1], vp[2], 0.e0};
}
static std::ostream& operator<< (std::ostream& outs, tensor_ref_4 const & ten)
{
    return outs << "[" << ten._data[0] << "," << ten._data[1] << "," << ten._data[2] << "," << ten._data[3] << "]";
}
inline vec3 vcons3 (double r0, double r1, double r2)
{
    return __extension__ (vec3){r0, r1, r2, 0.e0};
}
inline vec4 vload4 (const double *vp)
{
    return __extension__ (vec4){vp[0], vp[1], vp[2], vp[3]};
}
template <typename TY, const int VOXSZ>
inline double world2image (diderot::image1d< double, TY, VOXSZ > const & img)
{
    return img.world2image();
}
inline vec4 vcons4 (double r0, double r1, double r2, double r3)
{
    return __extension__ (vec4){r0, r1, r2, r3};
}
inline void vpack3 (tensor_3 &dst, vec3 v0)
{
    dst._data[0] = v0[0];
    dst._data[1] = v0[1];
    dst._data[2] = v0[2];
}
template <typename TY, const int VOXSZ>
inline double translate (diderot::image1d< double, TY, VOXSZ > const & img)
{
    return img.translate();
}
inline vec3 vscale3 (double s, vec3 v)
{
    return __extension__ (vec3){s, s, s, 0.e0} * v;
}
diderot::dynseq< int32_t > mesh_geom_msh (void *index, msh *mesh, const double *data)
{
    //diderot::dynseq<int32_t> myFunction(void * index, msh * mesh, double * data)

    try {

    Index * idx = reinterpret_cast<Index*>(index);

    

    SpatialIndex::Region* r = 0;

    //float to double conversion here..

    double dataP[2] = {static_cast<double>(data[0]), static_cast<double>(data[1])};

    r = new SpatialIndex::Region(dataP, dataP, 2);

      

    int64_t nResultCount;

    

    

    IdVisitor* visitor = new IdVisitor;

    //do the actually query

    idx->index().intersectsWithQuery(*r, *visitor);

    nResultCount = visitor->GetResultCount();

    

    std::vector<uint64_t>& results = visitor->GetResults();

    uint32_t vecSize = results.size();

    diderot::dynseq<int32_t> result = diderot::dynseq<int32_t>(vecSize);

    //copy data to a dynamic sequence.

    for(auto it = results.cbegin(); it != results.cend(); ++it){

      result.append((static_cast<int32_t>(*it)));

    

     }

    //clean up and return.

    delete r;

    delete visitor;

    return(result);

    } catch (Tools::IllegalArgumentException) {

    diderot::dynseq<int32_t> result = diderot::dynseq<int32_t>();

    return(result);

    }

}
inline void vpack4 (tensor_4 &dst, vec4 v0)
{
    dst._data[0] = v0[0];
    dst._data[1] = v0[1];
    dst._data[2] = v0[2];
    dst._data[3] = v0[3];
}
inline double vdot2 (vec2 u, vec2 v)
{
    vec2 w = u * v;
    return w[0] + w[1];
}
inline double vdot3 (vec3 u, vec3 v)
{
    vec3 w = u * v;
    return w[0] + w[1] + w[2];
}
inline double vdot4 (vec4 u, vec4 v)
{
    vec4 w = u * v;
    return w[0] + w[1] + w[2] + w[3];
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
    wrld->_definedInp.gv_0b043B_intermedateGlobal = true;
    std::memcpy(&wrld->_globals->gv_0b043B_intermedateGlobal, v, sizeof(fns));
    return false;
}
extern "C" bool justTypes_input_set_c (justTypes_world_t *cWrld, void *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_0c043D_intermedateGlobal = true;
    std::memcpy(&wrld->_globals->gv_0c043D_intermedateGlobal, v, sizeof(FUNC));
    return false;
}
extern "C" void justTypes_input_get_isoval (justTypes_world_t *cWrld, double *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_isoval;
}
extern "C" bool justTypes_input_set_isoval (justTypes_world_t *cWrld, double v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_isoval = true;
    wrld->_globals->gv_isoval = v;
    return false;
}
extern "C" void justTypes_input_get_thick (justTypes_world_t *cWrld, double *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_thick;
}
extern "C" bool justTypes_input_set_thick (justTypes_world_t *cWrld, double v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_thick = true;
    wrld->_globals->gv_thick = v;
    return false;
}
extern "C" void justTypes_input_get_camEye (justTypes_world_t *cWrld, double v[3])
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    std::memcpy(v, wrld->_globals->gv_camEye.addr(0), sizeof(tensor_3));
}
extern "C" bool justTypes_input_set_camEye (justTypes_world_t *cWrld, double v[3])
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_camEye = true;
    wrld->_globals->gv_camEye = v;
    return false;
}
extern "C" void justTypes_input_get_camAt (justTypes_world_t *cWrld, double v[3])
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    std::memcpy(v, wrld->_globals->gv_camAt.addr(0), sizeof(tensor_3));
}
extern "C" bool justTypes_input_set_camAt (justTypes_world_t *cWrld, double v[3])
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_camAt = true;
    wrld->_globals->gv_camAt = v;
    return false;
}
extern "C" void justTypes_input_get_camUp (justTypes_world_t *cWrld, double v[3])
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    std::memcpy(v, wrld->_globals->gv_camUp.addr(0), sizeof(tensor_3));
}
extern "C" bool justTypes_input_set_camUp (justTypes_world_t *cWrld, double v[3])
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_camUp = true;
    wrld->_globals->gv_camUp = v;
    return false;
}
extern "C" void justTypes_input_get_camFOV (justTypes_world_t *cWrld, double *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_camFOV;
}
extern "C" bool justTypes_input_set_camFOV (justTypes_world_t *cWrld, double v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_camFOV = true;
    wrld->_globals->gv_camFOV = v;
    return false;
}
extern "C" void justTypes_input_get_iresU (justTypes_world_t *cWrld, int32_t *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_iresU;
}
extern "C" bool justTypes_input_set_iresU (justTypes_world_t *cWrld, int32_t v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_iresU = true;
    wrld->_globals->gv_iresU = v;
    return false;
}
extern "C" void justTypes_input_get_iresV (justTypes_world_t *cWrld, int32_t *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_iresV;
}
extern "C" bool justTypes_input_set_iresV (justTypes_world_t *cWrld, int32_t v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_iresV = true;
    wrld->_globals->gv_iresV = v;
    return false;
}
extern "C" void justTypes_input_get_camNear (justTypes_world_t *cWrld, double *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_camNear;
}
extern "C" bool justTypes_input_set_camNear (justTypes_world_t *cWrld, double v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_camNear = true;
    wrld->_globals->gv_camNear = v;
    return false;
}
extern "C" void justTypes_input_get_camFar (justTypes_world_t *cWrld, double *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_camFar;
}
extern "C" bool justTypes_input_set_camFar (justTypes_world_t *cWrld, double v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_camFar = true;
    wrld->_globals->gv_camFar = v;
    return false;
}
extern "C" void justTypes_input_get_refStep (justTypes_world_t *cWrld, double *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_refStep;
}
extern "C" bool justTypes_input_set_refStep (justTypes_world_t *cWrld, double v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_refStep = true;
    wrld->_globals->gv_refStep = v;
    return false;
}
extern "C" void justTypes_input_get_rayStep (justTypes_world_t *cWrld, double *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_rayStep;
}
extern "C" bool justTypes_input_set_rayStep (justTypes_world_t *cWrld, double v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_rayStep = true;
    wrld->_globals->gv_rayStep = v;
    return false;
}
extern "C" void justTypes_input_get_lightVsp (justTypes_world_t *cWrld, double v[3])
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    std::memcpy(v, wrld->_globals->gv_lightVsp.addr(0), sizeof(tensor_3));
}
extern "C" bool justTypes_input_set_lightVsp (justTypes_world_t *cWrld, double v[3])
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_lightVsp = true;
    wrld->_globals->gv_lightVsp = v;
    return false;
}
extern "C" void justTypes_input_get_phongKa (justTypes_world_t *cWrld, double *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_phongKa;
}
extern "C" bool justTypes_input_set_phongKa (justTypes_world_t *cWrld, double v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_phongKa = true;
    wrld->_globals->gv_phongKa = v;
    return false;
}
extern "C" void justTypes_input_get_phongKd (justTypes_world_t *cWrld, double *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_phongKd;
}
extern "C" bool justTypes_input_set_phongKd (justTypes_world_t *cWrld, double v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_phongKd = true;
    wrld->_globals->gv_phongKd = v;
    return false;
}
extern "C" void justTypes_input_get_debug (justTypes_world_t *cWrld, bool *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_debug;
}
extern "C" bool justTypes_input_set_debug (justTypes_world_t *cWrld, bool v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_debug = true;
    wrld->_globals->gv_debug = v;
    return false;
}
extern "C" void justTypes_input_get_su (justTypes_world_t *cWrld, int32_t *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_su;
}
extern "C" bool justTypes_input_set_su (justTypes_world_t *cWrld, int32_t v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_su = true;
    wrld->_globals->gv_su = v;
    return false;
}
extern "C" void justTypes_input_get_sv (justTypes_world_t *cWrld, int32_t *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    *v = wrld->_globals->gv_sv;
}
extern "C" bool justTypes_input_set_sv (justTypes_world_t *cWrld, int32_t v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_sv = true;
    wrld->_globals->gv_sv = v;
    return false;
}
static bool check_defined (world *wrld)
{
    if (!wrld->_definedInp.gv_a) {
        biffMsgAdd(wrld->_errors, "undefined input \"a\"\n");
        return true;
    }
    if (!wrld->_definedInp.gv_0b043B_intermedateGlobal) {
        biffMsgAdd(wrld->_errors, "undefined input \"b\"\n");
        return true;
    }
    if (!wrld->_definedInp.gv_0c043D_intermedateGlobal) {
        biffMsgAdd(wrld->_errors, "undefined input \"c\"\n");
        return true;
    }
    return false;
}
static void init_defined_inputs (world *wrld)
{
    wrld->_definedInp.gv_a = false;
    wrld->_definedInp.gv_0b043B_intermedateGlobal = false;
    wrld->_definedInp.gv_0c043D_intermedateGlobal = false;
    wrld->_definedInp.gv_isoval = false;
    wrld->_definedInp.gv_thick = false;
    wrld->_definedInp.gv_camEye = false;
    wrld->_definedInp.gv_camAt = false;
    wrld->_definedInp.gv_camUp = false;
    wrld->_definedInp.gv_camFOV = false;
    wrld->_definedInp.gv_iresU = false;
    wrld->_definedInp.gv_iresV = false;
    wrld->_definedInp.gv_camNear = false;
    wrld->_definedInp.gv_camFar = false;
    wrld->_definedInp.gv_refStep = false;
    wrld->_definedInp.gv_rayStep = false;
    wrld->_definedInp.gv_lightVsp = false;
    wrld->_definedInp.gv_phongKa = false;
    wrld->_definedInp.gv_phongKd = false;
    wrld->_definedInp.gv_debug = false;
    wrld->_definedInp.gv_su = false;
    wrld->_definedInp.gv_sv = false;
}
static void init_defaults (globals *glob)
{
    glob->gv_isoval = 0.2e0;
    glob->gv_thick = 0.7e-1;
    glob->gv_camEye[0] = 0.3e1;
    glob->gv_camEye[1] = -0.3e1;
    glob->gv_camEye[2] = 0.3e1;
    glob->gv_camAt[0] = 0.5e0;
    glob->gv_camAt[1] = 0.5e0;
    glob->gv_camAt[2] = 0.5e0;
    vec3 v_2 = vcons3(0.e0, 0.e0, 0.1e1);
    vpack3(glob->gv_camUp, v_2);
    glob->gv_camFOV = 0.137e2;
    glob->gv_iresU = 500;
    glob->gv_iresV = 500;
    glob->gv_camNear = -0.2e1;
    glob->gv_camFar = 0.5e1;
    glob->gv_refStep = 0.1e1;
    glob->gv_rayStep = 0.1e-1;
    vpack3(glob->gv_lightVsp, v_2);
    glob->gv_phongKa = 0.1e0;
    glob->gv_phongKd = 0.9e0;
    glob->gv_debug = false;
    glob->gv_su = 250;
    glob->gv_sv = 250;
}
mesh_pos_msh fn_findPos (world *wrld, msh p_mesh_5, tensor_ref_3 p_pos_6)
{
    vec3 v_8;
    int32_t l_cellInt_153;
    int32_t l_newtonInt_154;
    vec3 v_155;
    int32_t l_numCell_7 = p_mesh_5.numCells - 1;
    v_8 = vcons3(
        0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e0,
        0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e0,
        0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e0);
    diderot::dynseq< int32_t > t_10 = mesh_geom_msh(p_mesh_5.index, &p_mesh_5, p_pos_6._data);
    int32_t l_cellInt_11 = 0;
    for (auto it_0 = t_10.cbegin(); it_0 != t_10.cend(); ++it_0) {
        auto i_cellItter_9 = *it_0;
        vec3 v_12;
        v_12 = v_8;
        for (int32_t i_newtonItter_13 = 0; i_newtonItter_13 <= 2; ++i_newtonItter_13) {
            int32_t l_mulRes_14 = i_cellItter_9 * 4;
            int32_t t_15 = p_mesh_5.indexMap[l_mulRes_14];
            int32_t l_mulRes_16 = 3 * t_15;
            double l_dof_load_17 = p_mesh_5.coordMap[l_mulRes_16];
            double l_dof_load_18 = p_mesh_5.coordMap[1 + l_mulRes_16];
            double l_dof_load_19 = p_mesh_5.coordMap[2 + l_mulRes_16];
            int32_t t_20 = p_mesh_5.indexMap[l_mulRes_14 + 1];
            int32_t l_mulRes_21 = 3 * t_20;
            double l_dof_load_22 = p_mesh_5.coordMap[l_mulRes_21];
            double l_dof_load_23 = p_mesh_5.coordMap[1 + l_mulRes_21];
            double l_dof_load_24 = p_mesh_5.coordMap[2 + l_mulRes_21];
            int32_t t_25 = p_mesh_5.indexMap[l_mulRes_14 + 2];
            int32_t l_mulRes_26 = 3 * t_25;
            double l_dof_load_27 = p_mesh_5.coordMap[l_mulRes_26];
            double l_dof_load_28 = p_mesh_5.coordMap[1 + l_mulRes_26];
            double l_dof_load_29 = p_mesh_5.coordMap[2 + l_mulRes_26];
            int32_t t_30 = p_mesh_5.indexMap[l_mulRes_14 + 3];
            int32_t l_mulRes_31 = 3 * t_30;
            double l_dof_load_32 = p_mesh_5.coordMap[l_mulRes_31];
            double l_dof_load_33 = p_mesh_5.coordMap[1 + l_mulRes_31];
            double l_dof_load_34 = p_mesh_5.coordMap[2 + l_mulRes_31];
            double l_prod_35 = 0.1e1 * 0.1e1;
            double l_prod_36 = 0.1e1 * l_prod_35;
            double l_basisEval_37 = -0.1e1 * l_prod_36;
            double l_basisEval_38 = 0.1e1 * l_prod_36;
            double l_r_39 = l_dof_load_17 * l_basisEval_37;
            double l_r_40 = l_dof_load_27 * 0.e0;
            double l_r_41 = l_dof_load_32 * 0.e0;
            double l_r_42 = l_r_39 + l_dof_load_22 * l_basisEval_38 + l_r_40 + l_r_41;
            double l_r_43 = l_r_39 + l_dof_load_22 * 0.e0;
            double l_r_44 = l_r_43 + l_dof_load_27 * l_basisEval_38 + l_r_41;
            double l_r_45 = l_r_43 + l_r_40 + l_dof_load_32 * l_basisEval_38;
            double l_r_46 = l_dof_load_18 * l_basisEval_37;
            double l_r_47 = l_dof_load_28 * 0.e0;
            double l_r_48 = l_dof_load_33 * 0.e0;
            double l_r_49 = l_r_46 + l_dof_load_23 * l_basisEval_38 + l_r_47 + l_r_48;
            double l_r_50 = l_r_46 + l_dof_load_23 * 0.e0;
            double l_r_51 = l_r_50 + l_dof_load_28 * l_basisEval_38 + l_r_48;
            double l_r_52 = l_r_50 + l_r_47 + l_dof_load_33 * l_basisEval_38;
            double l_r_53 = l_dof_load_19 * l_basisEval_37;
            double l_r_54 = l_dof_load_29 * 0.e0;
            double l_r_55 = l_dof_load_34 * 0.e0;
            double l_r_56 = l_r_53 + l_dof_load_24 * l_basisEval_38 + l_r_54 + l_r_55;
            double l_r_57 = l_r_53 + l_dof_load_24 * 0.e0;
            double l_r_58 = l_r_57 + l_dof_load_29 * l_basisEval_38 + l_r_55;
            double l_r_59 = l_r_57 + l_r_54 + l_dof_load_34 * l_basisEval_38;
            double l_r_60 = 0.e0 * l_r_42;
            double l_r_61 = 0.e0 * l_r_49;
            double l_r_62 = 0.e0 * l_r_56;
            double l_r_63 = l_r_60 + l_r_61;
            double l_r_64 = l_r_63 + l_r_62;
            double l_r_65 = 0.e0 * l_r_44;
            double l_r_66 = 0.e0 * l_r_51;
            double l_r_67 = 0.e0 * l_r_58;
            double l_r_68 = l_r_65 + l_r_66;
            double l_r_69 = l_r_68 + l_r_67;
            double l_r_70 = 0.e0 * l_r_45;
            double l_r_71 = 0.e0 * l_r_52;
            double l_r_72 = 0.e0 * l_r_59;
            double l_r_73 = l_r_70 + l_r_71;
            double l_r_74 = l_r_73 + l_r_72;
            double l_r_75 = l_r_63 + -0.1e1 * l_r_56;
            double l_r_76 = l_r_68 + -0.1e1 * l_r_58;
            double l_r_77 = l_r_73 + -0.1e1 * l_r_59;
            double l_r_78 = l_r_60 + 0.1e1 * l_r_49 + l_r_62;
            double l_r_79 = l_r_65 + 0.1e1 * l_r_51 + l_r_67;
            double l_r_80 = l_r_70 + 0.1e1 * l_r_52 + l_r_72;
            double l_r_81 = l_r_63 + 0.1e1 * l_r_56;
            double l_r_82 = l_r_68 + 0.1e1 * l_r_58;
            double l_r_83 = l_r_73 + 0.1e1 * l_r_59;
            double l_r_84 = -0.1e1 * l_r_42 + l_r_61 + l_r_62;
            double l_r_85 = -0.1e1 * l_r_44 + l_r_66 + l_r_67;
            double l_r_86 = -0.1e1 * l_r_45 + l_r_71 + l_r_72;
            double l_r_87 = l_r_60 + -0.1e1 * l_r_49 + l_r_62;
            double l_r_88 = l_r_65 + -0.1e1 * l_r_51 + l_r_67;
            double l_r_89 = l_r_70 + -0.1e1 * l_r_52 + l_r_72;
            double l_r_90 = 0.1e1 * l_r_42 + l_r_61 + l_r_62;
            double l_r_91 = 0.1e1 * l_r_44 + l_r_66 + l_r_67;
            double l_r_92 = 0.1e1 * l_r_45 + l_r_71 + l_r_72;
            double l_r_93 = l_r_42 * l_r_69 + l_r_49 * l_r_82 + l_r_56 * l_r_88;
            double l_r_94 = l_r_42 * l_r_74 + l_r_49 * l_r_83 + l_r_56 * l_r_89;
            double l_r_95 = l_r_42 * l_r_76 + l_r_49 * l_r_69 + l_r_56 * l_r_91;
            double l_r_96 = l_r_42 * l_r_77 + l_r_49 * l_r_74 + l_r_56 * l_r_92;
            double l_r_97 = l_r_42 * l_r_79 + l_r_49 * l_r_85 + l_r_56 * l_r_69;
            double l_r_98 = l_r_42 * l_r_80 + l_r_49 * l_r_86 + l_r_56 * l_r_74;
            double l_r_99 = l_r_44 * l_r_64 + l_r_51 * l_r_81 + l_r_58 * l_r_87;
            double l_r_100 = l_r_44 * l_r_74 + l_r_51 * l_r_83 + l_r_58 * l_r_89;
            double l_r_101 = l_r_44 * l_r_75 + l_r_51 * l_r_64 + l_r_58 * l_r_90;
            double l_r_102 = l_r_44 * l_r_77 + l_r_51 * l_r_74 + l_r_58 * l_r_92;
            double l_r_103 = l_r_44 * l_r_78 + l_r_51 * l_r_84 + l_r_58 * l_r_64;
            double l_r_104 = l_r_44 * l_r_80 + l_r_51 * l_r_86 + l_r_58 * l_r_74;
            double l_r_105 = l_r_45 * l_r_64 + l_r_52 * l_r_81 + l_r_59 * l_r_87;
            double l_r_106 = l_r_45 * l_r_69 + l_r_52 * l_r_82 + l_r_59 * l_r_88;
            double l_r_107 = l_r_45 * l_r_75 + l_r_52 * l_r_64 + l_r_59 * l_r_90;
            double l_r_108 = l_r_45 * l_r_76 + l_r_52 * l_r_69 + l_r_59 * l_r_91;
            double l_r_109 = l_r_45 * l_r_78 + l_r_52 * l_r_84 + l_r_59 * l_r_64;
            double l_r_110 = l_r_45 * l_r_79 + l_r_52 * l_r_85 + l_r_59 * l_r_69;
            vec3 v_111 = vcons3(l_r_44, l_r_51, l_r_58);
            double l_r_112 = 0.e0 * (l_r_42 * l_r_64 + l_r_49 * l_r_81 + l_r_56 * l_r_87);
            double l_r_113 = 0.e0 * l_r_94;
            double l_r_114 = 0.e0 * l_r_99;
            double l_r_115 = 0.e0 * (l_r_44 * l_r_69 + l_r_51 * l_r_82 + l_r_58 * l_r_88);
            double l_r_116 = 0.e0 * l_r_105;
            double l_r_117 = 0.e0 * (l_r_45 * l_r_74 + l_r_52 * l_r_83 + l_r_59 * l_r_89);
            double l_r_118 = l_r_112 + 0.e0 * l_r_93;
            double l_r_119 = 0.e0 * (l_r_42 * l_r_75 + l_r_49 * l_r_64 + l_r_56 * l_r_90);
            double l_r_120 = 0.e0 * l_r_96;
            double l_r_121 = 0.e0 * l_r_101;
            double l_r_122 = 0.e0 * (l_r_44 * l_r_76 + l_r_51 * l_r_69 + l_r_58 * l_r_91);
            double l_r_123 = 0.e0 * l_r_107;
            double l_r_124 = 0.e0 * (l_r_45 * l_r_77 + l_r_52 * l_r_74 + l_r_59 * l_r_92);
            double l_r_125 = l_r_119 + 0.e0 * l_r_95;
            double l_r_126 = 0.e0 * (l_r_42 * l_r_78 + l_r_49 * l_r_84 + l_r_56 * l_r_64);
            double l_r_127 = 0.e0 * l_r_98;
            double l_r_128 = 0.e0 * l_r_103;
            double l_r_129 = 0.e0 * (l_r_44 * l_r_79 + l_r_51 * l_r_85 + l_r_58 * l_r_69);
            double l_r_130 = 0.e0 * l_r_109;
            double l_r_131 = 0.e0 * (l_r_45 * l_r_80 + l_r_52 * l_r_86 + l_r_59 * l_r_74);
            double l_r_132 = l_r_126 + 0.e0 * l_r_97;
            double l_r_133 = 0.e0 * l_r_100;
            double l_r_134 = 0.e0 * l_r_106;
            double l_r_135 = 0.e0 * l_r_102;
            double l_r_136 = 0.e0 * l_r_108;
            double l_r_137 = 0.e0 * l_r_104;
            double l_r_138 = 0.e0 * l_r_110;
            double l_op1_e3_l_21_139 = 0.2e1 * vdot3(vcons3(l_r_42, l_r_49, l_r_56),
                vcons3(vdot3(v_111, vcons3(l_r_74, l_r_83, l_r_89)), vdot3(v_111, vcons3(l_r_77, l_r_74, l_r_92)),
                    vdot3(v_111, vcons3(l_r_80, l_r_86, l_r_74))));
            double l_prod_140 = v_12[0] * l_prod_35;
            double l_prod_141 = 0.1e1 * (v_12[1] * 0.1e1);
            double l_prod_142 = 0.1e1 * (0.1e1 * v_12[2]);
            double l_sum_143 = l_basisEval_38 + (-0.1e1 * l_prod_142 + (-0.1e1 * l_prod_141 + -0.1e1 * l_prod_140));
            double l_basisEval_144 = 0.1e1 * l_prod_140;
            double l_basisEval_145 = 0.1e1 * l_prod_141;
            double l_basisEval_146 = 0.1e1 * l_prod_142;
            vec3 v_147 = vcons3(
                l_dof_load_17 * l_sum_143 + l_dof_load_22 * l_basisEval_144 + l_dof_load_27 * l_basisEval_145 + l_dof_load_32 * l_basisEval_146,
                l_dof_load_18 * l_sum_143 + l_dof_load_23 * l_basisEval_144 + l_dof_load_28 * l_basisEval_145 + l_dof_load_33 * l_basisEval_146,
                l_dof_load_19 * l_sum_143 + l_dof_load_24 * l_basisEval_144 + l_dof_load_29 * l_basisEval_145 + l_dof_load_34 * l_basisEval_146) - vload3(
                p_pos_6.addr(0));
            vec3 v_148 = vcons3(
                vdot3(
                    vcons3(
                        (l_r_118 + l_r_113 + l_r_114 + l_r_115 + 0.1e1 * l_r_100 + l_r_116 + -0.1e1 * l_r_106 + l_r_117) / l_op1_e3_l_21_139,
                        (l_r_125 + l_r_120 + l_r_121 + l_r_122 + 0.1e1 * l_r_102 + l_r_123 + -0.1e1 * l_r_108 + l_r_124) / l_op1_e3_l_21_139,
                        (l_r_132 + l_r_127 + l_r_128 + l_r_129 + 0.1e1 * l_r_104 + l_r_130 + -0.1e1 * l_r_110 + l_r_131) / l_op1_e3_l_21_139),
                    v_147),
                vdot3(
                    vcons3(
                        (l_r_118 + -0.1e1 * l_r_94 + l_r_114 + l_r_115 + l_r_133 + 0.1e1 * l_r_105 + l_r_134 + l_r_117) / l_op1_e3_l_21_139,
                        (l_r_125 + -0.1e1 * l_r_96 + l_r_121 + l_r_122 + l_r_135 + 0.1e1 * l_r_107 + l_r_136 + l_r_124) / l_op1_e3_l_21_139,
                        (l_r_132 + -0.1e1 * l_r_98 + l_r_128 + l_r_129 + l_r_137 + 0.1e1 * l_r_109 + l_r_138 + l_r_131) / l_op1_e3_l_21_139),
                    v_147),
                vdot3(
                    vcons3(
                        (l_r_112 + 0.1e1 * l_r_93 + l_r_113 + -0.1e1 * l_r_99 + l_r_115 + l_r_133 + l_r_116 + l_r_134 + l_r_117) / l_op1_e3_l_21_139,
                        (l_r_119 + 0.1e1 * l_r_95 + l_r_120 + -0.1e1 * l_r_101 + l_r_122 + l_r_135 + l_r_123 + l_r_136 + l_r_124) / l_op1_e3_l_21_139,
                        (l_r_126 + 0.1e1 * l_r_97 + l_r_127 + -0.1e1 * l_r_103 + l_r_129 + l_r_137 + l_r_130 + l_r_138 + l_r_131) / l_op1_e3_l_21_139),
                    v_147));
            vec3 v_149 = v_12 - v_148;
            vec3 v_150 = v_149;
            if (0.1e-7 * 0.1e-7 >= vdot3(v_148, v_148)) {
                vec3 v_151 = vcons3(0.1e-4, 0.1e-4, 0.1e-4) + v_150;
                if (0.1e1 + 0.1e-4 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_150) && (v_151[0] > -0.e0 && (v_151[1] > -0.e0 && v_151[2] > -0.e0))) {
                    tensor_3 _arg_152;
                    vpack3(_arg_152, v_150);
                    return allBuild(p_mesh_5, i_cellItter_9, _arg_152, p_pos_6, true, true);
                }
            }
            v_12 = v_150;
        }
        v_8 = v_12;
    }
    l_cellInt_153 = l_cellInt_11;
    l_newtonInt_154 = 0;
    v_155 = v_8;
    int32_t hi_1 = 3 * l_numCell_7;
    for (int32_t i_itter_156 = 0; i_itter_156 <= hi_1; ++i_itter_156) {
        int32_t l_cellInt_298;
        int32_t l_newtonInt_299;
        int32_t l_mulRes_157 = l_cellInt_153 * 4;
        int32_t t_158 = p_mesh_5.indexMap[l_mulRes_157];
        int32_t l_mulRes_159 = 3 * t_158;
        double l_dof_load_160 = p_mesh_5.coordMap[l_mulRes_159];
        double l_dof_load_161 = p_mesh_5.coordMap[1 + l_mulRes_159];
        double l_dof_load_162 = p_mesh_5.coordMap[2 + l_mulRes_159];
        int32_t t_163 = p_mesh_5.indexMap[l_mulRes_157 + 1];
        int32_t l_mulRes_164 = 3 * t_163;
        double l_dof_load_165 = p_mesh_5.coordMap[l_mulRes_164];
        double l_dof_load_166 = p_mesh_5.coordMap[1 + l_mulRes_164];
        double l_dof_load_167 = p_mesh_5.coordMap[2 + l_mulRes_164];
        int32_t t_168 = p_mesh_5.indexMap[l_mulRes_157 + 2];
        int32_t l_mulRes_169 = 3 * t_168;
        double l_dof_load_170 = p_mesh_5.coordMap[l_mulRes_169];
        double l_dof_load_171 = p_mesh_5.coordMap[1 + l_mulRes_169];
        double l_dof_load_172 = p_mesh_5.coordMap[2 + l_mulRes_169];
        int32_t t_173 = p_mesh_5.indexMap[l_mulRes_157 + 3];
        int32_t l_mulRes_174 = 3 * t_173;
        double l_dof_load_175 = p_mesh_5.coordMap[l_mulRes_174];
        double l_dof_load_176 = p_mesh_5.coordMap[1 + l_mulRes_174];
        double l_dof_load_177 = p_mesh_5.coordMap[2 + l_mulRes_174];
        double l_prod_178 = 0.1e1 * 0.1e1;
        double l_prod_179 = 0.1e1 * l_prod_178;
        double l_basisEval_180 = -0.1e1 * l_prod_179;
        double l_basisEval_181 = 0.1e1 * l_prod_179;
        double l_r_182 = l_dof_load_160 * l_basisEval_180;
        double l_r_183 = l_dof_load_170 * 0.e0;
        double l_r_184 = l_dof_load_175 * 0.e0;
        double l_r_185 = l_r_182 + l_dof_load_165 * l_basisEval_181 + l_r_183 + l_r_184;
        double l_r_186 = l_r_182 + l_dof_load_165 * 0.e0;
        double l_r_187 = l_r_186 + l_dof_load_170 * l_basisEval_181 + l_r_184;
        double l_r_188 = l_r_186 + l_r_183 + l_dof_load_175 * l_basisEval_181;
        double l_r_189 = l_dof_load_161 * l_basisEval_180;
        double l_r_190 = l_dof_load_171 * 0.e0;
        double l_r_191 = l_dof_load_176 * 0.e0;
        double l_r_192 = l_r_189 + l_dof_load_166 * l_basisEval_181 + l_r_190 + l_r_191;
        double l_r_193 = l_r_189 + l_dof_load_166 * 0.e0;
        double l_r_194 = l_r_193 + l_dof_load_171 * l_basisEval_181 + l_r_191;
        double l_r_195 = l_r_193 + l_r_190 + l_dof_load_176 * l_basisEval_181;
        double l_r_196 = l_dof_load_162 * l_basisEval_180;
        double l_r_197 = l_dof_load_172 * 0.e0;
        double l_r_198 = l_dof_load_177 * 0.e0;
        double l_r_199 = l_r_196 + l_dof_load_167 * l_basisEval_181 + l_r_197 + l_r_198;
        double l_r_200 = l_r_196 + l_dof_load_167 * 0.e0;
        double l_r_201 = l_r_200 + l_dof_load_172 * l_basisEval_181 + l_r_198;
        double l_r_202 = l_r_200 + l_r_197 + l_dof_load_177 * l_basisEval_181;
        double l_r_203 = 0.e0 * l_r_185;
        double l_r_204 = 0.e0 * l_r_192;
        double l_r_205 = 0.e0 * l_r_199;
        double l_r_206 = l_r_203 + l_r_204;
        double l_r_207 = l_r_206 + l_r_205;
        double l_r_208 = 0.e0 * l_r_187;
        double l_r_209 = 0.e0 * l_r_194;
        double l_r_210 = 0.e0 * l_r_201;
        double l_r_211 = l_r_208 + l_r_209;
        double l_r_212 = l_r_211 + l_r_210;
        double l_r_213 = 0.e0 * l_r_188;
        double l_r_214 = 0.e0 * l_r_195;
        double l_r_215 = 0.e0 * l_r_202;
        double l_r_216 = l_r_213 + l_r_214;
        double l_r_217 = l_r_216 + l_r_215;
        double l_r_218 = l_r_206 + -0.1e1 * l_r_199;
        double l_r_219 = l_r_211 + -0.1e1 * l_r_201;
        double l_r_220 = l_r_216 + -0.1e1 * l_r_202;
        double l_r_221 = l_r_203 + 0.1e1 * l_r_192 + l_r_205;
        double l_r_222 = l_r_208 + 0.1e1 * l_r_194 + l_r_210;
        double l_r_223 = l_r_213 + 0.1e1 * l_r_195 + l_r_215;
        double l_r_224 = l_r_206 + 0.1e1 * l_r_199;
        double l_r_225 = l_r_211 + 0.1e1 * l_r_201;
        double l_r_226 = l_r_216 + 0.1e1 * l_r_202;
        double l_r_227 = -0.1e1 * l_r_185 + l_r_204 + l_r_205;
        double l_r_228 = -0.1e1 * l_r_187 + l_r_209 + l_r_210;
        double l_r_229 = -0.1e1 * l_r_188 + l_r_214 + l_r_215;
        double l_r_230 = l_r_203 + -0.1e1 * l_r_192 + l_r_205;
        double l_r_231 = l_r_208 + -0.1e1 * l_r_194 + l_r_210;
        double l_r_232 = l_r_213 + -0.1e1 * l_r_195 + l_r_215;
        double l_r_233 = 0.1e1 * l_r_185 + l_r_204 + l_r_205;
        double l_r_234 = 0.1e1 * l_r_187 + l_r_209 + l_r_210;
        double l_r_235 = 0.1e1 * l_r_188 + l_r_214 + l_r_215;
        double l_r_236 = l_r_185 * l_r_212 + l_r_192 * l_r_225 + l_r_199 * l_r_231;
        double l_r_237 = l_r_185 * l_r_217 + l_r_192 * l_r_226 + l_r_199 * l_r_232;
        double l_r_238 = l_r_185 * l_r_219 + l_r_192 * l_r_212 + l_r_199 * l_r_234;
        double l_r_239 = l_r_185 * l_r_220 + l_r_192 * l_r_217 + l_r_199 * l_r_235;
        double l_r_240 = l_r_185 * l_r_222 + l_r_192 * l_r_228 + l_r_199 * l_r_212;
        double l_r_241 = l_r_185 * l_r_223 + l_r_192 * l_r_229 + l_r_199 * l_r_217;
        double l_r_242 = l_r_187 * l_r_207 + l_r_194 * l_r_224 + l_r_201 * l_r_230;
        double l_r_243 = l_r_187 * l_r_217 + l_r_194 * l_r_226 + l_r_201 * l_r_232;
        double l_r_244 = l_r_187 * l_r_218 + l_r_194 * l_r_207 + l_r_201 * l_r_233;
        double l_r_245 = l_r_187 * l_r_220 + l_r_194 * l_r_217 + l_r_201 * l_r_235;
        double l_r_246 = l_r_187 * l_r_221 + l_r_194 * l_r_227 + l_r_201 * l_r_207;
        double l_r_247 = l_r_187 * l_r_223 + l_r_194 * l_r_229 + l_r_201 * l_r_217;
        double l_r_248 = l_r_188 * l_r_207 + l_r_195 * l_r_224 + l_r_202 * l_r_230;
        double l_r_249 = l_r_188 * l_r_212 + l_r_195 * l_r_225 + l_r_202 * l_r_231;
        double l_r_250 = l_r_188 * l_r_218 + l_r_195 * l_r_207 + l_r_202 * l_r_233;
        double l_r_251 = l_r_188 * l_r_219 + l_r_195 * l_r_212 + l_r_202 * l_r_234;
        double l_r_252 = l_r_188 * l_r_221 + l_r_195 * l_r_227 + l_r_202 * l_r_207;
        double l_r_253 = l_r_188 * l_r_222 + l_r_195 * l_r_228 + l_r_202 * l_r_212;
        vec3 v_254 = vcons3(l_r_187, l_r_194, l_r_201);
        double l_r_255 = 0.e0 * (l_r_185 * l_r_207 + l_r_192 * l_r_224 + l_r_199 * l_r_230);
        double l_r_256 = 0.e0 * l_r_237;
        double l_r_257 = 0.e0 * l_r_242;
        double l_r_258 = 0.e0 * (l_r_187 * l_r_212 + l_r_194 * l_r_225 + l_r_201 * l_r_231);
        double l_r_259 = 0.e0 * l_r_248;
        double l_r_260 = 0.e0 * (l_r_188 * l_r_217 + l_r_195 * l_r_226 + l_r_202 * l_r_232);
        double l_r_261 = l_r_255 + 0.e0 * l_r_236;
        double l_r_262 = 0.e0 * (l_r_185 * l_r_218 + l_r_192 * l_r_207 + l_r_199 * l_r_233);
        double l_r_263 = 0.e0 * l_r_239;
        double l_r_264 = 0.e0 * l_r_244;
        double l_r_265 = 0.e0 * (l_r_187 * l_r_219 + l_r_194 * l_r_212 + l_r_201 * l_r_234);
        double l_r_266 = 0.e0 * l_r_250;
        double l_r_267 = 0.e0 * (l_r_188 * l_r_220 + l_r_195 * l_r_217 + l_r_202 * l_r_235);
        double l_r_268 = l_r_262 + 0.e0 * l_r_238;
        double l_r_269 = 0.e0 * (l_r_185 * l_r_221 + l_r_192 * l_r_227 + l_r_199 * l_r_207);
        double l_r_270 = 0.e0 * l_r_241;
        double l_r_271 = 0.e0 * l_r_246;
        double l_r_272 = 0.e0 * (l_r_187 * l_r_222 + l_r_194 * l_r_228 + l_r_201 * l_r_212);
        double l_r_273 = 0.e0 * l_r_252;
        double l_r_274 = 0.e0 * (l_r_188 * l_r_223 + l_r_195 * l_r_229 + l_r_202 * l_r_217);
        double l_r_275 = l_r_269 + 0.e0 * l_r_240;
        double l_r_276 = 0.e0 * l_r_243;
        double l_r_277 = 0.e0 * l_r_249;
        double l_r_278 = 0.e0 * l_r_245;
        double l_r_279 = 0.e0 * l_r_251;
        double l_r_280 = 0.e0 * l_r_247;
        double l_r_281 = 0.e0 * l_r_253;
        double l_op1_e3_l_21_282 = 0.2e1 * vdot3(vcons3(l_r_185, l_r_192, l_r_199),
            vcons3(vdot3(v_254, vcons3(l_r_217, l_r_226, l_r_232)), vdot3(v_254, vcons3(l_r_220, l_r_217, l_r_235)),
                vdot3(v_254, vcons3(l_r_223, l_r_229, l_r_217))));
        double l_prod_283 = v_155[0] * l_prod_178;
        double l_prod_284 = 0.1e1 * (v_155[1] * 0.1e1);
        double l_prod_285 = 0.1e1 * (0.1e1 * v_155[2]);
        double l_sum_286 = l_basisEval_181 + (-0.1e1 * l_prod_285 + (-0.1e1 * l_prod_284 + -0.1e1 * l_prod_283));
        double l_basisEval_287 = 0.1e1 * l_prod_283;
        double l_basisEval_288 = 0.1e1 * l_prod_284;
        double l_basisEval_289 = 0.1e1 * l_prod_285;
        vec3 v_290 = vcons3(
            l_dof_load_160 * l_sum_286 + l_dof_load_165 * l_basisEval_287 + l_dof_load_170 * l_basisEval_288 + l_dof_load_175 * l_basisEval_289,
            l_dof_load_161 * l_sum_286 + l_dof_load_166 * l_basisEval_287 + l_dof_load_171 * l_basisEval_288 + l_dof_load_176 * l_basisEval_289,
            l_dof_load_162 * l_sum_286 + l_dof_load_167 * l_basisEval_287 + l_dof_load_172 * l_basisEval_288 + l_dof_load_177 * l_basisEval_289) - vload3(
            p_pos_6.addr(0));
        vec3 v_291 = vcons3(
            vdot3(
                vcons3(
                    (l_r_261 + l_r_256 + l_r_257 + l_r_258 + 0.1e1 * l_r_243 + l_r_259 + -0.1e1 * l_r_249 + l_r_260) / l_op1_e3_l_21_282,
                    (l_r_268 + l_r_263 + l_r_264 + l_r_265 + 0.1e1 * l_r_245 + l_r_266 + -0.1e1 * l_r_251 + l_r_267) / l_op1_e3_l_21_282,
                    (l_r_275 + l_r_270 + l_r_271 + l_r_272 + 0.1e1 * l_r_247 + l_r_273 + -0.1e1 * l_r_253 + l_r_274) / l_op1_e3_l_21_282),
                v_290),
            vdot3(
                vcons3(
                    (l_r_261 + -0.1e1 * l_r_237 + l_r_257 + l_r_258 + l_r_276 + 0.1e1 * l_r_248 + l_r_277 + l_r_260) / l_op1_e3_l_21_282,
                    (l_r_268 + -0.1e1 * l_r_239 + l_r_264 + l_r_265 + l_r_278 + 0.1e1 * l_r_250 + l_r_279 + l_r_267) / l_op1_e3_l_21_282,
                    (l_r_275 + -0.1e1 * l_r_241 + l_r_271 + l_r_272 + l_r_280 + 0.1e1 * l_r_252 + l_r_281 + l_r_274) / l_op1_e3_l_21_282),
                v_290),
            vdot3(
                vcons3(
                    (l_r_255 + 0.1e1 * l_r_236 + l_r_256 + -0.1e1 * l_r_242 + l_r_258 + l_r_276 + l_r_259 + l_r_277 + l_r_260) / l_op1_e3_l_21_282,
                    (l_r_262 + 0.1e1 * l_r_238 + l_r_263 + -0.1e1 * l_r_244 + l_r_265 + l_r_278 + l_r_266 + l_r_279 + l_r_267) / l_op1_e3_l_21_282,
                    (l_r_269 + 0.1e1 * l_r_240 + l_r_270 + -0.1e1 * l_r_246 + l_r_272 + l_r_280 + l_r_273 + l_r_281 + l_r_274) / l_op1_e3_l_21_282),
                v_290));
        vec3 v_292 = v_155 - v_291;
        vec3 v_293 = v_292;
        if (0.1e-7 * 0.1e-7 >= vdot3(v_291, v_291)) {
            vec3 v_294 = vcons3(0.1e-4, 0.1e-4, 0.1e-4) + v_293;
            if (0.1e1 + 0.1e-4 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_293) && (v_294[0] > -0.e0 && (v_294[1] > -0.e0 && v_294[2] > -0.e0))) {
                tensor_3 _arg_295;
                vpack3(_arg_295, v_293);
                return allBuild(p_mesh_5, l_cellInt_153, _arg_295, p_pos_6, true, true);
            }
        }
        int32_t l_newtonInt_296 = l_newtonInt_154 + 1;
        if (l_newtonInt_296 >= 2) {
            int32_t l_cellInt_297;
            if (l_cellInt_153 >= l_numCell_7) {
                return invalidBuild(p_mesh_5);
            }
            else {
                l_cellInt_297 = l_cellInt_153 + 1;
            }
            l_cellInt_298 = l_cellInt_297;
            l_newtonInt_299 = 0;
        }
        else {
            l_cellInt_298 = l_cellInt_153;
            l_newtonInt_299 = l_newtonInt_296;
        }
        l_cellInt_153 = l_cellInt_298;
        l_newtonInt_154 = l_newtonInt_299;
        v_155 = v_293;
    }
    wrld->print() << "Bad end 2" << "\n" << std::flush;
    return invalidBuild(p_mesh_5);
}
static bool init_globals (world *wrld)
{
    diderot::dynseq< mesh_cell_msh > l__t_300;
    diderot::image1d< double, double, 3 > l_I_304;
    globals *glob = wrld->_globals;
    l__t_300 = {};
    int32_t hi_2 = glob->gv_a.numCells - 1;
    for (int32_t i__t_301 = 0; i__t_301 <= hi_2; ++i__t_301) {
        l__t_300 = diderot::dynseq< mesh_cell_msh >::append(l__t_300, makeFem(glob->gv_a, i__t_301));
    }
    FUNC l_c_302 = glob->gv_0c043D_intermedateGlobal.loadFem(glob->gv_0b043B_intermedateGlobal.loadFem(glob->gv_a));
    fns l__t_303 = l_c_302.space;
    {
        diderot::nrrd_proxy proxy = diderot::nrrd_proxy(311, 3);
        if (l_I_304.load(wrld, "cmap.nrrd", &proxy)) {
            return true;
        }
    }
    vec3 v_305 = vload3(tensor_ref_3(glob->gv_camAt).addr(0)) - vload3(tensor_ref_3(glob->gv_camEye).addr(0));
    double l_camDist_306 = std::sqrt(vdot3(v_305, v_305));
    glob->gv_camDist = l_camDist_306;
    double l_op1_e3_l_9_307 = 0.1e1 / l_camDist_306;
    vec3 v_308 = vscale3(l_op1_e3_l_9_307, v_305);
    vpack3(glob->gv_camN, v_308);
    double l_r_310 = tensor_ref_3(glob->gv_camUp)[0];
    double l_r_311 = 0.e0 * l_r_310;
    double l_r_312 = tensor_ref_3(glob->gv_camUp)[1];
    double l_r_313 = 0.e0 * l_r_312;
    double l_r_314 = tensor_ref_3(glob->gv_camUp)[2];
    double l_r_315 = 0.e0 * l_r_314;
    double l_r_316 = l_r_311 + l_r_313;
    double l_r_317 = l_r_316 + l_r_315;
    vec3 v_318 = vcons3(vdot3(v_305, vcons3(l_r_317, l_r_316 + 0.1e1 * l_r_314, l_r_311 + -0.1e1 * l_r_312 + l_r_315)),
        vdot3(v_305, vcons3(l_r_316 + -0.1e1 * l_r_314, l_r_317, 0.1e1 * l_r_310 + l_r_313 + l_r_315)),
        vdot3(v_305, vcons3(l_r_311 + 0.1e1 * l_r_312 + l_r_315, -0.1e1 * l_r_310 + l_r_313 + l_r_315, l_r_317)));
    double l_op1_e3_l_38_319 = 0.1e1 / std::sqrt(l_op1_e3_l_9_307 * l_op1_e3_l_9_307 * vdot3(v_318, v_318));
    vec3 v_320 = vscale3(l_op1_e3_l_38_319, vscale3(l_op1_e3_l_9_307, v_318));
    vpack3(glob->gv_camU, v_320);
    double l_r_322 = v_305[0];
    double l_r_323 = 0.e0 * l_r_322;
    double l_r_324 = v_305[1];
    double l_r_325 = 0.e0 * l_r_324;
    double l_r_326 = v_305[2];
    double l_r_327 = 0.e0 * l_r_326;
    double l_r_328 = l_r_323 + l_r_325;
    double l_r_329 = l_r_328 + l_r_327;
    vec3 v_330 = vscale3(l_op1_e3_l_9_307,
        vcons3(vdot3(v_318, vcons3(l_r_329, l_r_328 + 0.1e1 * l_r_326, l_r_323 + -0.1e1 * l_r_324 + l_r_327)),
            vdot3(v_318, vcons3(l_r_328 + -0.1e1 * l_r_326, l_r_329, 0.1e1 * l_r_322 + l_r_325 + l_r_327)),
            vdot3(v_318, vcons3(l_r_323 + 0.1e1 * l_r_324 + l_r_327, -0.1e1 * l_r_322 + l_r_325 + l_r_327, l_r_329))));
    double l_r_331 = l_op1_e3_l_38_319 * l_op1_e3_l_9_307;
    double l_r_332 = l_r_331 * v_330[0];
    double l_r_333 = l_r_331 * v_330[1];
    double l_r_334 = l_r_331 * v_330[2];
    glob->gv_camV[0] = l_r_332;
    glob->gv_camV[1] = l_r_333;
    glob->gv_camV[2] = l_r_334;
    double l_op1_e3_l_10_336 = std::tan(glob->gv_camFOV * 0.314159265358979323846264338327950288e1 / 0.36e3);
    glob->gv_camVmax = l_op1_e3_l_10_336 * l_camDist_306;
    glob->gv_camUmax = static_cast<double>(glob->gv_iresU) * l_op1_e3_l_10_336 * l_camDist_306 / static_cast<double>(glob->gv_iresV);
    double l_r_337 = tensor_ref_3(glob->gv_lightVsp)[0];
    double l_r_338 = tensor_ref_3(glob->gv_lightVsp)[1];
    double l_r_339 = tensor_ref_3(glob->gv_lightVsp)[2];
    vpack3(glob->gv_light,
        vscale3(
            0.1e1 / std::sqrt(
                vdot3(vload3(tensor_ref_3(glob->gv_lightVsp).addr(0)), vload3(tensor_ref_3(glob->gv_lightVsp).addr(0)))),
            vcons3(v_320[0] * l_r_337 + l_r_332 * l_r_338 + v_308[0] * l_r_339,
                v_320[1] * l_r_337 + l_r_333 * l_r_338 + v_308[1] * l_r_339,
                v_320[2] * l_r_337 + l_r_334 * l_r_338 + v_308[2] * l_r_339)));
    glob->gv__t = l__t_303.mesh;
    glob->gv__tX = l__t_303;
    glob->gv_c = l_c_302;
    glob->gv_I = l_I_304;
    glob->gv_I.register_global();
    return false;
}
static void raycast_init (globals *glob, raycast_strand *self, int32_t p_ui_341, int32_t p_vi_342)
{
    double l_op1_e3_l_9_343 = -glob->gv_camUmax;
    double l_rayU_344 = l_op1_e3_l_9_343 + (static_cast<double>(p_ui_341) - (-0.5e0)) / (static_cast<double>(glob->gv_iresU) - 0.5e0 - (-0.5e0)) * (glob->gv_camUmax - l_op1_e3_l_9_343);
    double l_rayV_345 = glob->gv_camVmax + (static_cast<double>(p_vi_342) - (-0.5e0)) / (static_cast<double>(glob->gv_iresV) - 0.5e0 - (-0.5e0)) * (-glob->gv_camVmax - glob->gv_camVmax);
    vec3 v_346 = vscale3(l_rayU_344, vload3(tensor_ref_3(glob->gv_camU).addr(0))) + vscale3(l_rayV_345,
        vload3(tensor_ref_3(glob->gv_camV).addr(0)));
    double l_r_347 = 0.1e1 / glob->gv_camDist;
    self->sv_rayU = l_rayU_344;
    self->sv_rayV = l_rayV_345;
    self->sv_rayN = glob->gv_camNear;
    vpack3(self->sv_rayVec,
        vload3(tensor_ref_3(glob->gv_camN).addr(0)) + vcons3(l_r_347 * v_346[0], l_r_347 * v_346[1],
            l_r_347 * v_346[2]));
    self->sv_transp = 0.1e1;
    self->sv_rgb[0] = 0.e0;
    self->sv_rgb[1] = 0.e0;
    self->sv_rgb[2] = 0.e0;
    self->sv_rgba[0] = 0.e0;
    self->sv_rgba[1] = 0.e0;
    self->sv_rgba[2] = 0.e0;
    self->sv_rgba[3] = 0.e0;
    self->sv_gray = 0.e0;
    self->sv_ui = p_ui_341;
    self->sv_vi = p_vi_342;
}
static diderot::strand_status raycast_update (world *wrld, globals *glob, raycast_strand *self)
{
    vec3 v_582;
    double l_transp_583;
    double l_transp_585;
    if (glob->gv_debug) {
        bool l__t_351;
        if (self->sv_ui != glob->gv_su) {
            l__t_351 = true;
        }
        else {
            l__t_351 = self->sv_vi != glob->gv_sv;
        }
        if (l__t_351) {
            return diderot::kStabilize;
        }
    }
    vec3 v_352 = vload3(tensor_ref_3(glob->gv_camEye).addr(0)) + vscale3(self->sv_rayN,
        vload3(tensor_ref_3(self->sv_rayVec).addr(0)));
    double l__t_353 = v_352[0];
    vec3 v_354 = v_352;
    if (0.e0 < l__t_353) {
        vec3 v_580;
        double l_transp_581;
        if (l__t_353 < 0.1e1) {
            vec3 v_578;
            double l_transp_579;
            double l__t_355 = v_354[1];
            if (0.e0 < l__t_355) {
                vec3 v_576;
                double l_transp_577;
                if (l__t_355 < 0.1e1) {
                    vec3 v_574;
                    double l_transp_575;
                    double l__t_356 = v_354[2];
                    if (0.e0 < l__t_356) {
                        vec3 v_572;
                        double l_transp_573;
                        if (l__t_356 < 0.1e1) {
                            vec3 v_569;
                            double l_transp_570;
                            tensor_3 _arg_357;
                            vpack3(_arg_357, v_354);
                            mesh_pos_msh l_p_358 = fn_findPos(wrld, glob->gv_a, _arg_357);
                            if (l_p_358.valid) {
                                vec3 v_567;
                                double l_transp_568;
                                tensor_3 _arg_359;
                                vpack3(_arg_359, v_354);
                                mesh_pos_msh l_callFindPos_360 = fn_findPos(wrld, glob->gv__t, _arg_359);
                                int32_t l_intPos_361 = l_callFindPos_360.cell;
                                tensor_ref_3 l_refPos_362 = l_callFindPos_360.refPos;
                                int32_t l_mulRes_363 = l_intPos_361 * 10;
                                int32_t t_364 = glob->gv__tX.indexMap[l_mulRes_363];
                                int32_t t_365 = glob->gv__tX.indexMap[l_mulRes_363 + 1];
                                int32_t t_366 = glob->gv__tX.indexMap[l_mulRes_363 + 2];
                                int32_t t_367 = glob->gv__tX.indexMap[l_mulRes_363 + 3];
                                int32_t t_368 = glob->gv__tX.indexMap[l_mulRes_363 + 4];
                                int32_t t_369 = glob->gv__tX.indexMap[l_mulRes_363 + 5];
                                int32_t t_370 = glob->gv__tX.indexMap[l_mulRes_363 + 6];
                                int32_t t_371 = glob->gv__tX.indexMap[l_mulRes_363 + 7];
                                int32_t t_372 = glob->gv__tX.indexMap[l_mulRes_363 + 8];
                                int32_t t_373 = glob->gv__tX.indexMap[l_mulRes_363 + 9];
                                double t_374 = glob->gv_c.coordMap[1 * t_373];
                                double t_375 = glob->gv_c.coordMap[1 * t_372];
                                double t_376 = glob->gv_c.coordMap[1 * t_371];
                                double t_377 = glob->gv_c.coordMap[1 * t_370];
                                double t_378 = glob->gv_c.coordMap[1 * t_369];
                                double t_379 = glob->gv_c.coordMap[1 * t_368];
                                double t_380 = glob->gv_c.coordMap[1 * t_367];
                                double t_381 = glob->gv_c.coordMap[1 * t_366];
                                double t_382 = glob->gv_c.coordMap[1 * t_365];
                                double t_383 = glob->gv_c.coordMap[1 * t_364];
                                vec4 v_384 = vcons4(t_383, t_382, t_381, t_380);
                                vec4 v_385 = vcons4(t_379, t_378, t_377, t_376);
                                vec2 v_386 = vcons2(t_375, t_374);
                                double l_varAcc_387 = l_refPos_362[0];
                                double l_varAcc_388 = l_refPos_362[1];
                                double l_varAcc_389 = l_refPos_362[2];
                                double l_prod_390 = 0.1e1 * 0.1e1;
                                double l_prod_391 = l_varAcc_387 * l_varAcc_387 * l_prod_390;
                                double l_prod_392 = l_varAcc_388 * 0.1e1;
                                double l_prod_393 = l_varAcc_387 * l_prod_392;
                                double l_prod_394 = 0.1e1 * l_varAcc_389;
                                double l_prod_395 = l_varAcc_387 * l_prod_394;
                                double l_prod_396 = l_varAcc_387 * l_prod_390;
                                double l_prod_397 = 0.1e1 * (l_varAcc_388 * l_varAcc_388 * 0.1e1);
                                double l_prod_398 = 0.1e1 * (l_varAcc_388 * l_varAcc_389);
                                double l_prod_399 = 0.1e1 * l_prod_392;
                                double l_prod_400 = 0.1e1 * (0.1e1 * (l_varAcc_389 * l_varAcc_389));
                                double l_prod_401 = 0.1e1 * l_prod_394;
                                double l_prod_402 = 0.1e1 * l_prod_390;
                                double l_mult_403 = 0.1e1 * l_prod_402;
                                double l_mult_404 = 0.2e1 * l_prod_400;
                                double l_mult_405 = 0.4e1 * l_prod_398;
                                double l_mult_406 = 0.2e1 * l_prod_397;
                                double l_mult_407 = 0.4e1 * l_prod_395;
                                double l_mult_408 = 0.4e1 * l_prod_393;
                                double l_mult_409 = 0.2e1 * l_prod_391;
                                double l_mult_410 = 0.4e1 * l_prod_401;
                                double l_mult_411 = -0.4e1 * l_prod_398;
                                double l_mult_412 = -0.4e1 * l_prod_395;
                                double l_mult_413 = 0.4e1 * l_prod_399;
                                double l_mult_414 = -0.4e1 * l_prod_393;
                                double l_mult_415 = 0.4e1 * l_prod_396;
                                double l_compositionl_416 = vdot4(v_385,
                                    vcons4(l_mult_405, l_mult_407, l_mult_408,
                                        l_mult_410 + (-0.4e1 * l_prod_400 + (l_mult_411 + l_mult_412)))) + (vdot2(
                                    v_386,
                                    vcons2(l_mult_413 + (l_mult_411 + (-0.4e1 * l_prod_397 + l_mult_414)),
                                        l_mult_415 + (l_mult_412 + (l_mult_414 + -0.4e1 * l_prod_391)))) + vdot4(v_384,
                                    vcons4(
                                        l_mult_403 + (-0.3e1 * l_prod_401 + (l_mult_404 + (-0.3e1 * l_prod_399 + (l_mult_405 + (l_mult_406 + (-0.3e1 * l_prod_396 + (l_mult_407 + (l_mult_408 + l_mult_409)))))))),
                                        -0.1e1 * l_prod_396 + l_mult_409, -0.1e1 * l_prod_399 + l_mult_406,
                                        -0.1e1 * l_prod_401 + l_mult_404)));
                                double l_sum_417 = -0.3e1 * l_prod_402 + (l_mult_410 + (l_mult_413 + l_mult_415));
                                double l_mult_418 = -0.1e1 * l_prod_402;
                                double l_basisEval_419 = -0.4e1 * l_prod_401;
                                double l_basisEval_420 = -0.4e1 * l_prod_399;
                                double l_mult_421 = 0.4e1 * l_prod_402;
                                double l_mult_422 = -0.4e1 * l_prod_396;
                                double l_vdot_423 = vdot4(v_385, vcons4(0.e0, l_mult_410, l_mult_413, l_basisEval_419)) + (vdot2(
                                    v_386,
                                    vcons2(l_basisEval_420,
                                        l_mult_421 + (l_basisEval_419 + (l_basisEval_420 + -0.8e1 * l_prod_396)))) + vdot4(
                                    v_384, vcons4(l_sum_417, l_mult_418 + l_mult_415, 0.e0, 0.e0)));
                                double l_vdot_424 = vdot4(v_385, vcons4(l_mult_410, 0.e0, l_mult_415, l_basisEval_419)) + (vdot2(
                                    v_386,
                                    vcons2(l_mult_421 + (l_basisEval_419 + (-0.8e1 * l_prod_399 + l_mult_422)),
                                        l_mult_422)) + vdot4(v_384,
                                    vcons4(l_sum_417, 0.e0, l_mult_418 + l_mult_413, 0.e0)));
                                double l_vdot_425 = vdot4(v_385,
                                    vcons4(l_mult_413, l_mult_415, 0.e0,
                                        l_mult_421 + (-0.8e1 * l_prod_401 + (l_basisEval_420 + l_mult_422)))) + (vdot2(
                                    v_386, vcons2(l_basisEval_420, l_mult_422)) + vdot4(v_384,
                                    vcons4(l_sum_417, 0.e0, 0.e0, l_mult_418 + l_mult_410)));
                                int32_t l_mulRes_426 = l_intPos_361 * 4;
                                int32_t t_427 = glob->gv__t.indexMap[l_mulRes_426];
                                int32_t l_mulRes_428 = 3 * t_427;
                                int32_t t_429 = glob->gv__t.indexMap[l_mulRes_426 + 1];
                                int32_t l_mulRes_430 = 3 * t_429;
                                double l_dof_load_431 = glob->gv__t.coordMap[l_mulRes_430];
                                double l_dof_load_432 = glob->gv__t.coordMap[1 + l_mulRes_430];
                                double l_dof_load_433 = glob->gv__t.coordMap[2 + l_mulRes_430];
                                int32_t t_434 = glob->gv__t.indexMap[l_mulRes_426 + 2];
                                int32_t l_mulRes_435 = 3 * t_434;
                                double l_dof_load_436 = glob->gv__t.coordMap[l_mulRes_435];
                                double l_dof_load_437 = glob->gv__t.coordMap[1 + l_mulRes_435];
                                double l_dof_load_438 = glob->gv__t.coordMap[2 + l_mulRes_435];
                                int32_t t_439 = glob->gv__t.indexMap[l_mulRes_426 + 3];
                                int32_t l_mulRes_440 = 3 * t_439;
                                double l_dof_load_441 = glob->gv__t.coordMap[l_mulRes_440];
                                double l_dof_load_442 = glob->gv__t.coordMap[1 + l_mulRes_440];
                                double l_dof_load_443 = glob->gv__t.coordMap[2 + l_mulRes_440];
                                double t_444 = glob->gv__t.coordMap[l_mulRes_428];
                                double l_r_445 = t_444 * l_mult_418;
                                double l_r_446 = l_dof_load_436 * 0.e0;
                                double l_r_447 = l_dof_load_441 * 0.e0;
                                double l_r_448 = l_r_445 + l_dof_load_431 * l_mult_403 + l_r_446 + l_r_447;
                                double l_r_449 = l_r_445 + l_dof_load_431 * 0.e0;
                                double l_r_450 = l_r_449 + l_dof_load_436 * l_mult_403 + l_r_447;
                                double l_r_451 = l_r_449 + l_r_446 + l_dof_load_441 * l_mult_403;
                                double t_452 = glob->gv__t.coordMap[1 + l_mulRes_428];
                                double l_r_453 = t_452 * l_mult_418;
                                double l_r_454 = l_dof_load_437 * 0.e0;
                                double l_r_455 = l_dof_load_442 * 0.e0;
                                double l_r_456 = l_r_453 + l_dof_load_432 * l_mult_403 + l_r_454 + l_r_455;
                                double l_r_457 = l_r_453 + l_dof_load_432 * 0.e0;
                                double l_r_458 = l_r_457 + l_dof_load_437 * l_mult_403 + l_r_455;
                                double l_r_459 = l_r_457 + l_r_454 + l_dof_load_442 * l_mult_403;
                                double t_460 = glob->gv__t.coordMap[2 + l_mulRes_428];
                                double l_r_461 = t_460 * l_mult_418;
                                double l_r_462 = l_dof_load_438 * 0.e0;
                                double l_r_463 = l_dof_load_443 * 0.e0;
                                double l_r_464 = l_r_461 + l_dof_load_433 * l_mult_403 + l_r_462 + l_r_463;
                                double l_r_465 = l_r_461 + l_dof_load_433 * 0.e0;
                                double l_r_466 = l_r_465 + l_dof_load_438 * l_mult_403 + l_r_463;
                                double l_r_467 = l_r_465 + l_r_462 + l_dof_load_443 * l_mult_403;
                                double l_r_468 = 0.e0 * l_r_448;
                                double l_r_469 = 0.e0 * l_r_456;
                                double l_r_470 = 0.e0 * l_r_464;
                                double l_r_471 = l_r_468 + l_r_469;
                                double l_r_472 = l_r_471 + l_r_470;
                                double l_r_473 = 0.e0 * l_r_450;
                                double l_r_474 = 0.e0 * l_r_458;
                                double l_r_475 = 0.e0 * l_r_466;
                                double l_r_476 = l_r_473 + l_r_474;
                                double l_r_477 = l_r_476 + l_r_475;
                                double l_r_478 = 0.e0 * l_r_451;
                                double l_r_479 = 0.e0 * l_r_459;
                                double l_r_480 = 0.e0 * l_r_467;
                                double l_r_481 = l_r_478 + l_r_479;
                                double l_r_482 = l_r_481 + l_r_480;
                                double l_r_483 = l_r_471 + -0.1e1 * l_r_464;
                                double l_r_484 = l_r_476 + -0.1e1 * l_r_466;
                                double l_r_485 = l_r_481 + -0.1e1 * l_r_467;
                                double l_r_486 = l_r_468 + 0.1e1 * l_r_456 + l_r_470;
                                double l_r_487 = l_r_473 + 0.1e1 * l_r_458 + l_r_475;
                                double l_r_488 = l_r_478 + 0.1e1 * l_r_459 + l_r_480;
                                double l_r_489 = l_r_471 + 0.1e1 * l_r_464;
                                double l_r_490 = l_r_476 + 0.1e1 * l_r_466;
                                double l_r_491 = l_r_481 + 0.1e1 * l_r_467;
                                double l_r_492 = -0.1e1 * l_r_448 + l_r_469 + l_r_470;
                                double l_r_493 = -0.1e1 * l_r_450 + l_r_474 + l_r_475;
                                double l_r_494 = -0.1e1 * l_r_451 + l_r_479 + l_r_480;
                                double l_r_495 = l_r_468 + -0.1e1 * l_r_456 + l_r_470;
                                double l_r_496 = l_r_473 + -0.1e1 * l_r_458 + l_r_475;
                                double l_r_497 = l_r_478 + -0.1e1 * l_r_459 + l_r_480;
                                double l_r_498 = 0.1e1 * l_r_448 + l_r_469 + l_r_470;
                                double l_r_499 = 0.1e1 * l_r_450 + l_r_474 + l_r_475;
                                double l_r_500 = 0.1e1 * l_r_451 + l_r_479 + l_r_480;
                                double l_r_501 = l_r_448 * l_r_477 + l_r_456 * l_r_490 + l_r_464 * l_r_496;
                                double l_r_502 = l_r_448 * l_r_482 + l_r_456 * l_r_491 + l_r_464 * l_r_497;
                                double l_r_503 = l_r_448 * l_r_484 + l_r_456 * l_r_477 + l_r_464 * l_r_499;
                                double l_r_504 = l_r_448 * l_r_485 + l_r_456 * l_r_482 + l_r_464 * l_r_500;
                                double l_r_505 = l_r_448 * l_r_487 + l_r_456 * l_r_493 + l_r_464 * l_r_477;
                                double l_r_506 = l_r_448 * l_r_488 + l_r_456 * l_r_494 + l_r_464 * l_r_482;
                                double l_r_507 = l_r_450 * l_r_472 + l_r_458 * l_r_489 + l_r_466 * l_r_495;
                                double l_r_508 = l_r_450 * l_r_482 + l_r_458 * l_r_491 + l_r_466 * l_r_497;
                                double l_r_509 = l_r_450 * l_r_483 + l_r_458 * l_r_472 + l_r_466 * l_r_498;
                                double l_r_510 = l_r_450 * l_r_485 + l_r_458 * l_r_482 + l_r_466 * l_r_500;
                                double l_r_511 = l_r_450 * l_r_486 + l_r_458 * l_r_492 + l_r_466 * l_r_472;
                                double l_r_512 = l_r_450 * l_r_488 + l_r_458 * l_r_494 + l_r_466 * l_r_482;
                                double l_r_513 = l_r_451 * l_r_472 + l_r_459 * l_r_489 + l_r_467 * l_r_495;
                                double l_r_514 = l_r_451 * l_r_477 + l_r_459 * l_r_490 + l_r_467 * l_r_496;
                                double l_r_515 = l_r_451 * l_r_483 + l_r_459 * l_r_472 + l_r_467 * l_r_498;
                                double l_r_516 = l_r_451 * l_r_484 + l_r_459 * l_r_477 + l_r_467 * l_r_499;
                                double l_r_517 = l_r_451 * l_r_486 + l_r_459 * l_r_492 + l_r_467 * l_r_472;
                                double l_r_518 = l_r_451 * l_r_487 + l_r_459 * l_r_493 + l_r_467 * l_r_477;
                                vec3 v_519 = vcons3(l_r_450, l_r_458, l_r_466);
                                double l_r_520 = 0.e0 * (l_r_448 * l_r_472 + l_r_456 * l_r_489 + l_r_464 * l_r_495);
                                double l_r_521 = 0.e0 * l_r_502;
                                double l_r_522 = 0.e0 * l_r_507;
                                double l_r_523 = 0.e0 * (l_r_450 * l_r_477 + l_r_458 * l_r_490 + l_r_466 * l_r_496);
                                double l_r_524 = 0.e0 * l_r_513;
                                double l_r_525 = 0.e0 * (l_r_451 * l_r_482 + l_r_459 * l_r_491 + l_r_467 * l_r_497);
                                double l_r_526 = l_r_520 + 0.e0 * l_r_501;
                                double l_r_527 = 0.e0 * (l_r_448 * l_r_483 + l_r_456 * l_r_472 + l_r_464 * l_r_498);
                                double l_r_528 = 0.e0 * l_r_504;
                                double l_r_529 = 0.e0 * l_r_509;
                                double l_r_530 = 0.e0 * (l_r_450 * l_r_484 + l_r_458 * l_r_477 + l_r_466 * l_r_499);
                                double l_r_531 = 0.e0 * l_r_515;
                                double l_r_532 = 0.e0 * (l_r_451 * l_r_485 + l_r_459 * l_r_482 + l_r_467 * l_r_500);
                                double l_r_533 = l_r_527 + 0.e0 * l_r_503;
                                double l_r_534 = 0.e0 * (l_r_448 * l_r_486 + l_r_456 * l_r_492 + l_r_464 * l_r_472);
                                double l_r_535 = 0.e0 * l_r_506;
                                double l_r_536 = 0.e0 * l_r_511;
                                double l_r_537 = 0.e0 * (l_r_450 * l_r_487 + l_r_458 * l_r_493 + l_r_466 * l_r_477);
                                double l_r_538 = 0.e0 * l_r_517;
                                double l_r_539 = 0.e0 * (l_r_451 * l_r_488 + l_r_459 * l_r_494 + l_r_467 * l_r_482);
                                double l_r_540 = l_r_534 + 0.e0 * l_r_505;
                                double l_r_541 = 0.e0 * l_r_508;
                                double l_r_542 = 0.e0 * l_r_514;
                                double l_r_543 = 0.e0 * l_r_510;
                                double l_r_544 = 0.e0 * l_r_516;
                                double l_r_545 = 0.e0 * l_r_512;
                                double l_r_546 = 0.e0 * l_r_518;
                                double l_op1_e3_l_36_547 = 0.2e1 * vdot3(vcons3(l_r_448, l_r_456, l_r_464),
                                    vcons3(vdot3(v_519, vcons3(l_r_482, l_r_491, l_r_497)),
                                        vdot3(v_519, vcons3(l_r_485, l_r_482, l_r_500)),
                                        vdot3(v_519, vcons3(l_r_488, l_r_494, l_r_482))));
                                vec3 v_548 = -vcons3(
                                    l_vdot_423 * ((l_r_526 + l_r_521 + l_r_522 + l_r_523 + 0.1e1 * l_r_508 + l_r_524 + -0.1e1 * l_r_514 + l_r_525) / l_op1_e3_l_36_547) + l_vdot_424 * ((l_r_526 + -0.1e1 * l_r_502 + l_r_522 + l_r_523 + l_r_541 + 0.1e1 * l_r_513 + l_r_542 + l_r_525) / l_op1_e3_l_36_547) + l_vdot_425 * ((l_r_520 + 0.1e1 * l_r_501 + l_r_521 + -0.1e1 * l_r_507 + l_r_523 + l_r_541 + l_r_524 + l_r_542 + l_r_525) / l_op1_e3_l_36_547),
                                    l_vdot_423 * ((l_r_533 + l_r_528 + l_r_529 + l_r_530 + 0.1e1 * l_r_510 + l_r_531 + -0.1e1 * l_r_516 + l_r_532) / l_op1_e3_l_36_547) + l_vdot_424 * ((l_r_533 + -0.1e1 * l_r_504 + l_r_529 + l_r_530 + l_r_543 + 0.1e1 * l_r_515 + l_r_544 + l_r_532) / l_op1_e3_l_36_547) + l_vdot_425 * ((l_r_527 + 0.1e1 * l_r_503 + l_r_528 + -0.1e1 * l_r_509 + l_r_530 + l_r_543 + l_r_531 + l_r_544 + l_r_532) / l_op1_e3_l_36_547),
                                    l_vdot_423 * ((l_r_540 + l_r_535 + l_r_536 + l_r_537 + 0.1e1 * l_r_512 + l_r_538 + -0.1e1 * l_r_518 + l_r_539) / l_op1_e3_l_36_547) + l_vdot_424 * ((l_r_540 + -0.1e1 * l_r_506 + l_r_536 + l_r_537 + l_r_545 + 0.1e1 * l_r_517 + l_r_546 + l_r_539) / l_op1_e3_l_36_547) + l_vdot_425 * ((l_r_534 + 0.1e1 * l_r_505 + l_r_535 + -0.1e1 * l_r_511 + l_r_537 + l_r_545 + l_r_538 + l_r_546 + l_r_539) / l_op1_e3_l_36_547));
                                double l_op1_e3_l_57_549 = std::sqrt(vdot3(v_548, v_548));
                                double l_a_550 = 0.1e1 * clamp(0.e0, 0.1e1,
                                    0.13e1 * (0.1e1 - std::abs(l_compositionl_416 - glob->gv_isoval) / (glob->gv_thick * l_op1_e3_l_57_549)));
                                vec3 v_551 = v_548;
                                if (l_a_550 > 0.e0) {
                                    tensor_3_2 l_voxels_561;
                                    double l_imgPos_552 = world2image(glob->gv_I) * l_compositionl_416 + translate(
                                        glob->gv_I);
                                    double l_nd_553 = std::floor(l_imgPos_552);
                                    double l_f_554 = l_imgPos_552 - l_nd_553;
                                    int32_t l_n_555 = std::lround(l_nd_553);
                                    double l__t_556 = std::pow(0.1e1 - l_a_550,
                                        glob->gv_rayStep * std::sqrt(
                                            vdot3(vload3(tensor_ref_3(self->sv_rayVec).addr(0)),
                                                vload3(tensor_ref_3(self->sv_rayVec).addr(0)))) / glob->gv_refStep);
                                    double l_op1_e3_l_19_557 = (self->sv_rayN - glob->gv_camNear) / (glob->gv_camFar - glob->gv_camNear) * (0.7e0 - 0.11e1);
                                    double l_op1_e3_l_20_558 = glob->gv_phongKd * std::max(0.e0,
                                        0.1e1 / l_op1_e3_l_57_549 * vdot3(v_551,
                                            vload3(tensor_ref_3(glob->gv_light).addr(0))));
                                    if (glob->gv_I.inside(l_n_555, 2)) {
                                        int32_t l_offp_559 = 3 * l_n_555;
                                        int32_t l_offp_560 = 3 * (l_n_555 + 1);
                                        l_voxels_561[0] = glob->gv_I[l_offp_559];
                                        l_voxels_561[1] = glob->gv_I[l_offp_560];
                                        l_voxels_561[2] = glob->gv_I[l_offp_559 + 1];
                                        l_voxels_561[3] = glob->gv_I[l_offp_560 + 1];
                                        l_voxels_561[4] = glob->gv_I[l_offp_559 + 2];
                                        l_voxels_561[5] = glob->gv_I[l_offp_560 + 2];
                                    }
                                    else {
                                        int32_t l_offp_562 = 3 * glob->gv_I.clamp(0, l_n_555);
                                        int32_t l_offp_563 = 3 * glob->gv_I.clamp(0, l_n_555 + 1);
                                        l_voxels_561[0] = glob->gv_I[l_offp_562];
                                        l_voxels_561[1] = glob->gv_I[l_offp_563];
                                        l_voxels_561[2] = glob->gv_I[l_offp_562 + 1];
                                        l_voxels_561[3] = glob->gv_I[l_offp_563 + 1];
                                        l_voxels_561[4] = glob->gv_I[l_offp_562 + 2];
                                        l_voxels_561[5] = glob->gv_I[l_offp_563 + 2];
                                    }
                                    vec2 v_564 = vcons2(0.1e1, 0.1e1) + vcons2(l_f_554, l_f_554 - 0.1e1) * vcons2(
                                        -0.1e1, 0.1e1);
                                    double l_op1_e3_l_22_565 = 0.1e1 - l__t_556;
                                    double l_r_566 = self->sv_transp * l_op1_e3_l_22_565 * (0.11e1 + l_op1_e3_l_19_557) * (glob->gv_phongKa + l_op1_e3_l_20_558);
                                    v_567 = vload3(tensor_ref_3(self->sv_rgb).addr(0)) + vcons3(
                                        l_r_566 * vdot2(vload2(l_voxels_561.last(0).addr(0)), v_564),
                                        l_r_566 * vdot2(vload2(l_voxels_561.last(2).addr(0)), v_564),
                                        l_r_566 * vdot2(vload2(l_voxels_561.last(4).addr(0)), v_564));
                                    l_transp_568 = self->sv_transp * (0.1e1 - l_op1_e3_l_22_565);
                                }
                                else {
                                    v_567 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                                    l_transp_568 = self->sv_transp;
                                }
                                wrld->print() << "yay!\n" << std::flush;
                                v_569 = v_567;
                                l_transp_570 = l_transp_568;
                            }
                            else {
                                tensor_3 _arg_571;
                                vpack3(_arg_571, v_354);
                                wrld->print() << tensor_ref_3(_arg_571) << "\n" << std::flush;
                                wrld->print() << "wait what\?\n" << std::flush;
                                v_569 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                                l_transp_570 = self->sv_transp;
                            }
                            v_572 = v_569;
                            l_transp_573 = l_transp_570;
                        }
                        else {
                            v_572 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                            l_transp_573 = self->sv_transp;
                        }
                        v_574 = v_572;
                        l_transp_575 = l_transp_573;
                    }
                    else {
                        v_574 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                        l_transp_575 = self->sv_transp;
                    }
                    v_576 = v_574;
                    l_transp_577 = l_transp_575;
                }
                else {
                    v_576 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                    l_transp_577 = self->sv_transp;
                }
                v_578 = v_576;
                l_transp_579 = l_transp_577;
            }
            else {
                v_578 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                l_transp_579 = self->sv_transp;
            }
            v_580 = v_578;
            l_transp_581 = l_transp_579;
        }
        else {
            v_580 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
            l_transp_581 = self->sv_transp;
        }
        v_582 = v_580;
        l_transp_583 = l_transp_581;
    }
    else {
        v_582 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
        l_transp_583 = self->sv_transp;
    }
    if (l_transp_583 < 0.1e-1) {
        self->sv_transp = 0.e0;
        vpack3(self->sv_rgb, v_582);
        return diderot::kStabilize;
    }
    else {
        l_transp_585 = l_transp_583;
    }
    if (self->sv_rayN > glob->gv_camFar) {
        self->sv_transp = l_transp_585;
        vpack3(self->sv_rgb, v_582);
        return diderot::kStabilize;
    }
    self->sv_rayN = self->sv_rayN + glob->gv_rayStep;
    self->sv_transp = l_transp_585;
    vpack3(self->sv_rgb, v_582);
    return diderot::kActive;
}
static void raycast_stabilize (world *wrld, globals *glob, raycast_strand *self)
{
    vec4 v_589;
    double l_a_588 = 0.1e1 - self->sv_transp;
    if (l_a_588 > 0.e0) {
        v_589 = vcons4(tensor_ref_3(self->sv_rgb)[0] / l_a_588, tensor_ref_3(self->sv_rgb)[1] / l_a_588,
            tensor_ref_3(self->sv_rgb)[2] / l_a_588, l_a_588);
    }
    else {
        v_589 = vload4(tensor_ref_4(self->sv_rgba).addr(0));
    }
    if (self->sv_ui == glob->gv_su) {
        if (self->sv_vi == glob->gv_sv) {
            if (glob->gv_debug) {
                tensor_4 _arg_590;
                vpack4(_arg_590, v_589);
                wrld->print() << l_a_588 << tensor_ref_4(_arg_590) << std::flush;
            }
        }
    }
    vpack4(self->sv_rgba, v_589);
}
extern "C" bool justTypes_output_get_rgba (justTypes_world_t *cWrld, Nrrd *nData)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    // Compute sizes of nrrd file
    size_t sizes[3];
    sizes[0] = 4;
    sizes[1] = wrld->_size[1];
    sizes[2] = wrld->_size[0];
    // Allocate nData nrrd
    if (nrrdMaybeAlloc_nva(nData, nrrdTypeDouble, 3, sizes) != 0) {
        char *msg = biffGetDone(NRRD);
        biffMsgAdd(wrld->_errors, msg);
        std::free(msg);
        return true;
    }
    // copy data to output nrrd
    char *cp = reinterpret_cast<char *>(nData->data);
    for (auto ix = wrld->_strands.begin_alive(); ix != wrld->_strands.end_alive(); ix = wrld->_strands.next_alive(ix)) {
        memcpy(cp, &wrld->_strands.strand(ix)->sv_rgba, 4 * sizeof(double));
        cp += 4 * sizeof(double);
    }
    nData->axis[0].kind = nrrdKind4Vector;
    nData->axis[1].kind = nrrdKindSpace;
    nData->axis[2].kind = nrrdKindSpace;
    return false;
}
extern "C" bool justTypes_output_get_gray (justTypes_world_t *cWrld, Nrrd *nData)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    // Compute sizes of nrrd file
    size_t sizes[2];
    sizes[0] = wrld->_size[1];
    sizes[1] = wrld->_size[0];
    // Allocate nData nrrd
    if (nrrdMaybeAlloc_nva(nData, nrrdTypeDouble, 2, sizes) != 0) {
        char *msg = biffGetDone(NRRD);
        biffMsgAdd(wrld->_errors, msg);
        std::free(msg);
        return true;
    }
    // copy data to output nrrd
    char *cp = reinterpret_cast<char *>(nData->data);
    for (auto ix = wrld->_strands.begin_alive(); ix != wrld->_strands.end_alive(); ix = wrld->_strands.next_alive(ix)) {
        memcpy(cp, &wrld->_strands.strand(ix)->sv_gray, 1 * sizeof(double));
        cp += 1 * sizeof(double);
    }
    nData->axis[0].kind = nrrdKindSpace;
    nData->axis[1].kind = nrrdKindSpace;
    return false;
}
/*---------- begin world-methods.in ----------*/
// Allocate the program's world
//
world::world ()
    : diderot::world_base (ProgramName, true, 2)
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
bool world::alloc (int32_t base[2], uint32_t size[2])
{
    size_t numStrands = 1;
    for (uint32_t i = 0;  i < 2;  i++) {
        numStrands *= size[i];
        this->_base[i] = base[i];
        this->_size[i] = size[i];
    }

    if (this->_verbose) {
        std::cerr << "world::alloc: " << size[0];
        for (uint32_t i = 1;  i < 2;  i++) {
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
            ix = this->_strands.strand_stabilize (this, glob, ix);
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
    int lo_3 = 0;
    int hi_4 = glob->gv_iresV - 1;
    int lo_5 = 0;
    int hi_6 = glob->gv_iresU - 1;
    int32_t base[2] = {lo_3,lo_5,};
    uint32_t size[2] = {static_cast<uint32_t>(hi_4 - lo_3 + 1),static_cast<uint32_t>(hi_6 - lo_5 + 1),};
    if (this->alloc(base, size)) {
        return true;
    }
    uint32_t ix = 0;
    for (int i_vi_592 = lo_3; i_vi_592 <= hi_4; i_vi_592++) {
        for (int i_ui_593 = lo_5; i_ui_593 <= hi_6; i_ui_593++) {
            raycast_init(this->_globals, this->_strands.strand(ix), i_ui_593, i_vi_592);
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
            ix = this->_strands.strand_stabilize (this, glob, ix);
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


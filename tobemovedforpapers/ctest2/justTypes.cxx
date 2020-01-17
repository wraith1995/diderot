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
    struct tensor_ref_4 : public diderot::tensor_ref<double,4> {
        tensor_ref_4 (const double *src);
        tensor_ref_4 (struct tensor_4 const & ten);
        tensor_ref_4 (tensor_ref_4 const & ten);
    };
    struct tensor_ref_7 : public diderot::tensor_ref<double,7> {
        tensor_ref_7 (const double *src);
        tensor_ref_7 (struct tensor_7 const & ten);
        tensor_ref_7 (tensor_ref_7 const & ten);
    };
    struct tensor_ref_3 : public diderot::tensor_ref<double,3> {
        tensor_ref_3 (const double *src);
        tensor_ref_3 (struct tensor_3 const & ten);
        tensor_ref_3 (tensor_ref_3 const & ten);
    };
    struct tensor_ref_2 : public diderot::tensor_ref<double,2> {
        tensor_ref_2 (const double *src);
        tensor_ref_2 (struct tensor_2 const & ten);
        tensor_ref_2 (tensor_ref_2 const & ten);
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
    inline tensor_ref_4::tensor_ref_4 (const double *src)
        : diderot::tensor_ref<double,4>(src)
    { }
    inline tensor_ref_4::tensor_ref_4 (struct tensor_4 const & ten)
        : diderot::tensor_ref<double,4>(ten._data)
    { }
    inline tensor_ref_4::tensor_ref_4 (tensor_ref_4 const & ten)
        : diderot::tensor_ref<double,4>(ten._data)
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
    inline tensor_ref_3::tensor_ref_3 (const double *src)
        : diderot::tensor_ref<double,3>(src)
    { }
    inline tensor_ref_3::tensor_ref_3 (struct tensor_3 const & ten)
        : diderot::tensor_ref<double,3>(ten._data)
    { }
    inline tensor_ref_3::tensor_ref_3 (tensor_ref_3 const & ten)
        : diderot::tensor_ref<double,3>(ten._data)
    { }
    inline tensor_ref_2::tensor_ref_2 (const double *src)
        : diderot::tensor_ref<double,2>(src)
    { }
    inline tensor_ref_2::tensor_ref_2 (struct tensor_2 const & ten)
        : diderot::tensor_ref<double,2>(ten._data)
    { }
    inline tensor_ref_2::tensor_ref_2 (tensor_ref_2 const & ten)
        : diderot::tensor_ref<double,2>(ten._data)
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

    double dataP[3] = {static_cast<double>(data[0]), static_cast<double>(data[1]), static_cast<double>(data[2])};

    r = new SpatialIndex::Region(dataP, dataP, 3);

      

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
    vec3 v_7;
    v_7 = vcons3(
        0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e0,
        0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e0,
        0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e0);
    diderot::dynseq< int32_t > t_9 = mesh_geom_msh(p_mesh_5.index, &p_mesh_5, p_pos_6._data);
    for (auto it_0 = t_9.cbegin(); it_0 != t_9.cend(); ++it_0) {
        auto i_cellItter_8 = *it_0;
        vec3 v_10;
        v_10 = v_7;
        for (int32_t i_newtonItter_11 = 0; i_newtonItter_11 <= 2; ++i_newtonItter_11) {
            int32_t l_mulRes_12 = i_cellItter_8 * 4;
            int32_t t_13 = p_mesh_5.indexMap[l_mulRes_12];
            int32_t l_mulRes_14 = 3 * t_13;
            double l_dof_load_15 = p_mesh_5.coordMap[l_mulRes_14];
            double l_dof_load_16 = p_mesh_5.coordMap[1 + l_mulRes_14];
            double l_dof_load_17 = p_mesh_5.coordMap[2 + l_mulRes_14];
            int32_t t_18 = p_mesh_5.indexMap[l_mulRes_12 + 1];
            int32_t l_mulRes_19 = 3 * t_18;
            double l_dof_load_20 = p_mesh_5.coordMap[l_mulRes_19];
            double l_dof_load_21 = p_mesh_5.coordMap[1 + l_mulRes_19];
            double l_dof_load_22 = p_mesh_5.coordMap[2 + l_mulRes_19];
            int32_t t_23 = p_mesh_5.indexMap[l_mulRes_12 + 2];
            int32_t l_mulRes_24 = 3 * t_23;
            double l_dof_load_25 = p_mesh_5.coordMap[l_mulRes_24];
            double l_dof_load_26 = p_mesh_5.coordMap[1 + l_mulRes_24];
            double l_dof_load_27 = p_mesh_5.coordMap[2 + l_mulRes_24];
            int32_t t_28 = p_mesh_5.indexMap[l_mulRes_12 + 3];
            int32_t l_mulRes_29 = 3 * t_28;
            double l_dof_load_30 = p_mesh_5.coordMap[l_mulRes_29];
            double l_dof_load_31 = p_mesh_5.coordMap[1 + l_mulRes_29];
            double l_dof_load_32 = p_mesh_5.coordMap[2 + l_mulRes_29];
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
            double l_op1_e3_l_21_137 = 0.2e1 * vdot3(vcons3(l_r_40, l_r_47, l_r_54),
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
                l_dof_load_15 * l_sum_141 + l_dof_load_20 * l_basisEval_142 + l_dof_load_25 * l_basisEval_143 + l_dof_load_30 * l_basisEval_144,
                l_dof_load_16 * l_sum_141 + l_dof_load_21 * l_basisEval_142 + l_dof_load_26 * l_basisEval_143 + l_dof_load_31 * l_basisEval_144,
                l_dof_load_17 * l_sum_141 + l_dof_load_22 * l_basisEval_142 + l_dof_load_27 * l_basisEval_143 + l_dof_load_32 * l_basisEval_144) - vload3(
                p_pos_6.addr(0));
            vec3 v_146 = vcons3(
                vdot3(
                    vcons3(
                        (l_r_116 + l_r_111 + l_r_112 + l_r_113 + 0.1e1 * l_r_98 + l_r_114 + -0.1e1 * l_r_104 + l_r_115) / l_op1_e3_l_21_137,
                        (l_r_123 + l_r_118 + l_r_119 + l_r_120 + 0.1e1 * l_r_100 + l_r_121 + -0.1e1 * l_r_106 + l_r_122) / l_op1_e3_l_21_137,
                        (l_r_130 + l_r_125 + l_r_126 + l_r_127 + 0.1e1 * l_r_102 + l_r_128 + -0.1e1 * l_r_108 + l_r_129) / l_op1_e3_l_21_137),
                    v_145),
                vdot3(
                    vcons3(
                        (l_r_116 + -0.1e1 * l_r_92 + l_r_112 + l_r_113 + l_r_131 + 0.1e1 * l_r_103 + l_r_132 + l_r_115) / l_op1_e3_l_21_137,
                        (l_r_123 + -0.1e1 * l_r_94 + l_r_119 + l_r_120 + l_r_133 + 0.1e1 * l_r_105 + l_r_134 + l_r_122) / l_op1_e3_l_21_137,
                        (l_r_130 + -0.1e1 * l_r_96 + l_r_126 + l_r_127 + l_r_135 + 0.1e1 * l_r_107 + l_r_136 + l_r_129) / l_op1_e3_l_21_137),
                    v_145),
                vdot3(
                    vcons3(
                        (l_r_110 + 0.1e1 * l_r_91 + l_r_111 + -0.1e1 * l_r_97 + l_r_113 + l_r_131 + l_r_114 + l_r_132 + l_r_115) / l_op1_e3_l_21_137,
                        (l_r_117 + 0.1e1 * l_r_93 + l_r_118 + -0.1e1 * l_r_99 + l_r_120 + l_r_133 + l_r_121 + l_r_134 + l_r_122) / l_op1_e3_l_21_137,
                        (l_r_124 + 0.1e1 * l_r_95 + l_r_125 + -0.1e1 * l_r_101 + l_r_127 + l_r_135 + l_r_128 + l_r_136 + l_r_129) / l_op1_e3_l_21_137),
                    v_145));
            vec3 v_147 = v_10 - v_146;
            vec3 v_148 = v_147;
            if (0.1e-7 * 0.1e-7 >= vdot3(v_146, v_146)) {
                vec3 v_149 = vcons3(0.1e-4, 0.1e-4, 0.1e-4) + v_148;
                if (0.1e1 + 0.1e-4 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_148) && (v_149[0] > -0.e0 && (v_149[1] > -0.e0 && v_149[2] > -0.e0))) {
                    tensor_3 _arg_150;
                    vpack3(_arg_150, v_148);
                    return allBuild(p_mesh_5, i_cellItter_8, _arg_150, p_pos_6, true, true);
                }
            }
            v_10 = v_148;
        }
        v_7 = v_10;
    }
    wrld->print() << "Bad end 2" << "\n" << std::flush;
    return invalidBuild(p_mesh_5);
}
static bool init_globals (world *wrld)
{
    diderot::dynseq< mesh_cell_msh > l__t_151;
    diderot::image1d< double, double, 3 > l_I_155;
    globals *glob = wrld->_globals;
    l__t_151 = {};
    int32_t hi_1 = glob->gv_a.numCells - 1;
    for (int32_t i__t_152 = 0; i__t_152 <= hi_1; ++i__t_152) {
        l__t_151 = diderot::dynseq< mesh_cell_msh >::append(l__t_151, makeFem(glob->gv_a, i__t_152));
    }
    FUNC l_c_153 = glob->gv_0c043D_intermedateGlobal.loadFem(glob->gv_0b043B_intermedateGlobal.loadFem(glob->gv_a));
    fns l__t_154 = l_c_153.space;
    {
        diderot::nrrd_proxy proxy = diderot::nrrd_proxy(311, 3);
        if (l_I_155.load(wrld, "cmap.nrrd", &proxy)) {
            return true;
        }
    }
    vec3 v_156 = vload3(tensor_ref_3(glob->gv_camAt).addr(0)) - vload3(tensor_ref_3(glob->gv_camEye).addr(0));
    double l_camDist_157 = std::sqrt(vdot3(v_156, v_156));
    glob->gv_camDist = l_camDist_157;
    double l_op1_e3_l_9_158 = 0.1e1 / l_camDist_157;
    vec3 v_159 = vscale3(l_op1_e3_l_9_158, v_156);
    vpack3(glob->gv_camN, v_159);
    double l_r_161 = tensor_ref_3(glob->gv_camUp)[0];
    double l_r_162 = 0.e0 * l_r_161;
    double l_r_163 = tensor_ref_3(glob->gv_camUp)[1];
    double l_r_164 = 0.e0 * l_r_163;
    double l_r_165 = tensor_ref_3(glob->gv_camUp)[2];
    double l_r_166 = 0.e0 * l_r_165;
    double l_r_167 = l_r_162 + l_r_164;
    double l_r_168 = l_r_167 + l_r_166;
    vec3 v_169 = vcons3(vdot3(v_156, vcons3(l_r_168, l_r_167 + 0.1e1 * l_r_165, l_r_162 + -0.1e1 * l_r_163 + l_r_166)),
        vdot3(v_156, vcons3(l_r_167 + -0.1e1 * l_r_165, l_r_168, 0.1e1 * l_r_161 + l_r_164 + l_r_166)),
        vdot3(v_156, vcons3(l_r_162 + 0.1e1 * l_r_163 + l_r_166, -0.1e1 * l_r_161 + l_r_164 + l_r_166, l_r_168)));
    double l_op1_e3_l_38_170 = 0.1e1 / std::sqrt(l_op1_e3_l_9_158 * l_op1_e3_l_9_158 * vdot3(v_169, v_169));
    vec3 v_171 = vscale3(l_op1_e3_l_38_170, vscale3(l_op1_e3_l_9_158, v_169));
    vpack3(glob->gv_camU, v_171);
    double l_r_173 = v_156[0];
    double l_r_174 = 0.e0 * l_r_173;
    double l_r_175 = v_156[1];
    double l_r_176 = 0.e0 * l_r_175;
    double l_r_177 = v_156[2];
    double l_r_178 = 0.e0 * l_r_177;
    double l_r_179 = l_r_174 + l_r_176;
    double l_r_180 = l_r_179 + l_r_178;
    vec3 v_181 = vscale3(l_op1_e3_l_9_158,
        vcons3(vdot3(v_169, vcons3(l_r_180, l_r_179 + 0.1e1 * l_r_177, l_r_174 + -0.1e1 * l_r_175 + l_r_178)),
            vdot3(v_169, vcons3(l_r_179 + -0.1e1 * l_r_177, l_r_180, 0.1e1 * l_r_173 + l_r_176 + l_r_178)),
            vdot3(v_169, vcons3(l_r_174 + 0.1e1 * l_r_175 + l_r_178, -0.1e1 * l_r_173 + l_r_176 + l_r_178, l_r_180))));
    double l_r_182 = l_op1_e3_l_38_170 * l_op1_e3_l_9_158;
    double l_r_183 = l_r_182 * v_181[0];
    double l_r_184 = l_r_182 * v_181[1];
    double l_r_185 = l_r_182 * v_181[2];
    glob->gv_camV[0] = l_r_183;
    glob->gv_camV[1] = l_r_184;
    glob->gv_camV[2] = l_r_185;
    double l_op1_e3_l_10_187 = std::tan(glob->gv_camFOV * 0.314159265358979323846264338327950288e1 / 0.36e3);
    glob->gv_camVmax = l_op1_e3_l_10_187 * l_camDist_157;
    glob->gv_camUmax = static_cast<double>(glob->gv_iresU) * l_op1_e3_l_10_187 * l_camDist_157 / static_cast<double>(glob->gv_iresV);
    double l_r_188 = tensor_ref_3(glob->gv_lightVsp)[0];
    double l_r_189 = tensor_ref_3(glob->gv_lightVsp)[1];
    double l_r_190 = tensor_ref_3(glob->gv_lightVsp)[2];
    vpack3(glob->gv_light,
        vscale3(
            0.1e1 / std::sqrt(
                vdot3(vload3(tensor_ref_3(glob->gv_lightVsp).addr(0)), vload3(tensor_ref_3(glob->gv_lightVsp).addr(0)))),
            vcons3(v_171[0] * l_r_188 + l_r_183 * l_r_189 + v_159[0] * l_r_190,
                v_171[1] * l_r_188 + l_r_184 * l_r_189 + v_159[1] * l_r_190,
                v_171[2] * l_r_188 + l_r_185 * l_r_189 + v_159[2] * l_r_190)));
    glob->gv__t = l__t_154.mesh;
    glob->gv__tX = l__t_154;
    glob->gv_c = l_c_153;
    glob->gv_I = l_I_155;
    glob->gv_I.register_global();
    return false;
}
static void raycast_init (globals *glob, raycast_strand *self, int32_t p_ui_192, int32_t p_vi_193)
{
    double l_op1_e3_l_9_194 = -glob->gv_camUmax;
    double l_rayU_195 = l_op1_e3_l_9_194 + (static_cast<double>(p_ui_192) - (-0.5e0)) / (static_cast<double>(glob->gv_iresU) - 0.5e0 - (-0.5e0)) * (glob->gv_camUmax - l_op1_e3_l_9_194);
    double l_rayV_196 = glob->gv_camVmax + (static_cast<double>(p_vi_193) - (-0.5e0)) / (static_cast<double>(glob->gv_iresV) - 0.5e0 - (-0.5e0)) * (-glob->gv_camVmax - glob->gv_camVmax);
    vec3 v_197 = vscale3(l_rayU_195, vload3(tensor_ref_3(glob->gv_camU).addr(0))) + vscale3(l_rayV_196,
        vload3(tensor_ref_3(glob->gv_camV).addr(0)));
    double l_r_198 = 0.1e1 / glob->gv_camDist;
    self->sv_rayU = l_rayU_195;
    self->sv_rayV = l_rayV_196;
    self->sv_rayN = glob->gv_camNear;
    vpack3(self->sv_rayVec,
        vload3(tensor_ref_3(glob->gv_camN).addr(0)) + vcons3(l_r_198 * v_197[0], l_r_198 * v_197[1],
            l_r_198 * v_197[2]));
    self->sv_transp = 0.1e1;
    self->sv_rgb[0] = 0.e0;
    self->sv_rgb[1] = 0.e0;
    self->sv_rgb[2] = 0.e0;
    self->sv_rgba[0] = 0.e0;
    self->sv_rgba[1] = 0.e0;
    self->sv_rgba[2] = 0.e0;
    self->sv_rgba[3] = 0.e0;
    self->sv_gray = 0.e0;
    self->sv_ui = p_ui_192;
    self->sv_vi = p_vi_193;
}
static diderot::strand_status raycast_update (world *wrld, globals *glob, raycast_strand *self)
{
    vec3 v_433;
    double l_transp_434;
    double l_transp_436;
    if (glob->gv_debug) {
        bool l__t_202;
        if (self->sv_ui != glob->gv_su) {
            l__t_202 = true;
        }
        else {
            l__t_202 = self->sv_vi != glob->gv_sv;
        }
        if (l__t_202) {
            return diderot::kStabilize;
        }
    }
    vec3 v_203 = vload3(tensor_ref_3(glob->gv_camEye).addr(0)) + vscale3(self->sv_rayN,
        vload3(tensor_ref_3(self->sv_rayVec).addr(0)));
    double l__t_204 = v_203[0];
    vec3 v_205 = v_203;
    if (0.e0 < l__t_204) {
        vec3 v_431;
        double l_transp_432;
        if (l__t_204 < 0.1e1) {
            vec3 v_429;
            double l_transp_430;
            double l__t_206 = v_205[1];
            if (0.e0 < l__t_206) {
                vec3 v_427;
                double l_transp_428;
                if (l__t_206 < 0.1e1) {
                    vec3 v_425;
                    double l_transp_426;
                    double l__t_207 = v_205[2];
                    if (0.e0 < l__t_207) {
                        vec3 v_423;
                        double l_transp_424;
                        if (l__t_207 < 0.1e1) {
                            vec3 v_420;
                            double l_transp_421;
                            tensor_3 _arg_208;
                            vpack3(_arg_208, v_205);
                            mesh_pos_msh l_p_209 = fn_findPos(wrld, glob->gv_a, _arg_208);
                            if (l_p_209.valid) {
                                vec3 v_418;
                                double l_transp_419;
                                tensor_3 _arg_210;
                                vpack3(_arg_210, v_205);
                                mesh_pos_msh l_callFindPos_211 = fn_findPos(wrld, glob->gv__t, _arg_210);
                                int32_t l_intPos_212 = l_callFindPos_211.cell;
                                tensor_ref_3 l_refPos_213 = l_callFindPos_211.refPos;
                                int32_t l_mulRes_214 = l_intPos_212 * 10;
                                int32_t t_215 = glob->gv__tX.indexMap[l_mulRes_214];
                                int32_t t_216 = glob->gv__tX.indexMap[l_mulRes_214 + 1];
                                int32_t t_217 = glob->gv__tX.indexMap[l_mulRes_214 + 2];
                                int32_t t_218 = glob->gv__tX.indexMap[l_mulRes_214 + 3];
                                int32_t t_219 = glob->gv__tX.indexMap[l_mulRes_214 + 4];
                                int32_t t_220 = glob->gv__tX.indexMap[l_mulRes_214 + 5];
                                int32_t t_221 = glob->gv__tX.indexMap[l_mulRes_214 + 6];
                                int32_t t_222 = glob->gv__tX.indexMap[l_mulRes_214 + 7];
                                int32_t t_223 = glob->gv__tX.indexMap[l_mulRes_214 + 8];
                                int32_t t_224 = glob->gv__tX.indexMap[l_mulRes_214 + 9];
                                double t_225 = glob->gv_c.coordMap[1 * t_224];
                                double t_226 = glob->gv_c.coordMap[1 * t_223];
                                double t_227 = glob->gv_c.coordMap[1 * t_222];
                                double t_228 = glob->gv_c.coordMap[1 * t_221];
                                double t_229 = glob->gv_c.coordMap[1 * t_220];
                                double t_230 = glob->gv_c.coordMap[1 * t_219];
                                double t_231 = glob->gv_c.coordMap[1 * t_218];
                                double t_232 = glob->gv_c.coordMap[1 * t_217];
                                double t_233 = glob->gv_c.coordMap[1 * t_216];
                                double t_234 = glob->gv_c.coordMap[1 * t_215];
                                vec4 v_235 = vcons4(t_234, t_233, t_232, t_231);
                                vec4 v_236 = vcons4(t_230, t_229, t_228, t_227);
                                vec2 v_237 = vcons2(t_226, t_225);
                                double l_varAcc_238 = l_refPos_213[0];
                                double l_varAcc_239 = l_refPos_213[1];
                                double l_varAcc_240 = l_refPos_213[2];
                                double l_prod_241 = 0.1e1 * 0.1e1;
                                double l_prod_242 = l_varAcc_238 * l_varAcc_238 * l_prod_241;
                                double l_prod_243 = l_varAcc_239 * 0.1e1;
                                double l_prod_244 = l_varAcc_238 * l_prod_243;
                                double l_prod_245 = 0.1e1 * l_varAcc_240;
                                double l_prod_246 = l_varAcc_238 * l_prod_245;
                                double l_prod_247 = l_varAcc_238 * l_prod_241;
                                double l_prod_248 = 0.1e1 * (l_varAcc_239 * l_varAcc_239 * 0.1e1);
                                double l_prod_249 = 0.1e1 * (l_varAcc_239 * l_varAcc_240);
                                double l_prod_250 = 0.1e1 * l_prod_243;
                                double l_prod_251 = 0.1e1 * (0.1e1 * (l_varAcc_240 * l_varAcc_240));
                                double l_prod_252 = 0.1e1 * l_prod_245;
                                double l_prod_253 = 0.1e1 * l_prod_241;
                                double l_mult_254 = 0.1e1 * l_prod_253;
                                double l_mult_255 = 0.2e1 * l_prod_251;
                                double l_mult_256 = 0.4e1 * l_prod_249;
                                double l_mult_257 = 0.2e1 * l_prod_248;
                                double l_mult_258 = 0.4e1 * l_prod_246;
                                double l_mult_259 = 0.4e1 * l_prod_244;
                                double l_mult_260 = 0.2e1 * l_prod_242;
                                double l_mult_261 = 0.4e1 * l_prod_252;
                                double l_mult_262 = -0.4e1 * l_prod_249;
                                double l_mult_263 = -0.4e1 * l_prod_246;
                                double l_mult_264 = 0.4e1 * l_prod_250;
                                double l_mult_265 = -0.4e1 * l_prod_244;
                                double l_mult_266 = 0.4e1 * l_prod_247;
                                double l_compositionl_267 = vdot4(v_236,
                                    vcons4(l_mult_256, l_mult_258, l_mult_259,
                                        l_mult_261 + (-0.4e1 * l_prod_251 + (l_mult_262 + l_mult_263)))) + (vdot2(
                                    v_237,
                                    vcons2(l_mult_264 + (l_mult_262 + (-0.4e1 * l_prod_248 + l_mult_265)),
                                        l_mult_266 + (l_mult_263 + (l_mult_265 + -0.4e1 * l_prod_242)))) + vdot4(v_235,
                                    vcons4(
                                        l_mult_254 + (-0.3e1 * l_prod_252 + (l_mult_255 + (-0.3e1 * l_prod_250 + (l_mult_256 + (l_mult_257 + (-0.3e1 * l_prod_247 + (l_mult_258 + (l_mult_259 + l_mult_260)))))))),
                                        -0.1e1 * l_prod_247 + l_mult_260, -0.1e1 * l_prod_250 + l_mult_257,
                                        -0.1e1 * l_prod_252 + l_mult_255)));
                                double l_sum_268 = -0.3e1 * l_prod_253 + (l_mult_261 + (l_mult_264 + l_mult_266));
                                double l_mult_269 = -0.1e1 * l_prod_253;
                                double l_basisEval_270 = -0.4e1 * l_prod_252;
                                double l_basisEval_271 = -0.4e1 * l_prod_250;
                                double l_mult_272 = 0.4e1 * l_prod_253;
                                double l_mult_273 = -0.4e1 * l_prod_247;
                                double l_vdot_274 = vdot4(v_236, vcons4(0.e0, l_mult_261, l_mult_264, l_basisEval_270)) + (vdot2(
                                    v_237,
                                    vcons2(l_basisEval_271,
                                        l_mult_272 + (l_basisEval_270 + (l_basisEval_271 + -0.8e1 * l_prod_247)))) + vdot4(
                                    v_235, vcons4(l_sum_268, l_mult_269 + l_mult_266, 0.e0, 0.e0)));
                                double l_vdot_275 = vdot4(v_236, vcons4(l_mult_261, 0.e0, l_mult_266, l_basisEval_270)) + (vdot2(
                                    v_237,
                                    vcons2(l_mult_272 + (l_basisEval_270 + (-0.8e1 * l_prod_250 + l_mult_273)),
                                        l_mult_273)) + vdot4(v_235,
                                    vcons4(l_sum_268, 0.e0, l_mult_269 + l_mult_264, 0.e0)));
                                double l_vdot_276 = vdot4(v_236,
                                    vcons4(l_mult_264, l_mult_266, 0.e0,
                                        l_mult_272 + (-0.8e1 * l_prod_252 + (l_basisEval_271 + l_mult_273)))) + (vdot2(
                                    v_237, vcons2(l_basisEval_271, l_mult_273)) + vdot4(v_235,
                                    vcons4(l_sum_268, 0.e0, 0.e0, l_mult_269 + l_mult_261)));
                                int32_t l_mulRes_277 = l_intPos_212 * 4;
                                int32_t t_278 = glob->gv__t.indexMap[l_mulRes_277];
                                int32_t l_mulRes_279 = 3 * t_278;
                                int32_t t_280 = glob->gv__t.indexMap[l_mulRes_277 + 1];
                                int32_t l_mulRes_281 = 3 * t_280;
                                double l_dof_load_282 = glob->gv__t.coordMap[l_mulRes_281];
                                double l_dof_load_283 = glob->gv__t.coordMap[1 + l_mulRes_281];
                                double l_dof_load_284 = glob->gv__t.coordMap[2 + l_mulRes_281];
                                int32_t t_285 = glob->gv__t.indexMap[l_mulRes_277 + 2];
                                int32_t l_mulRes_286 = 3 * t_285;
                                double l_dof_load_287 = glob->gv__t.coordMap[l_mulRes_286];
                                double l_dof_load_288 = glob->gv__t.coordMap[1 + l_mulRes_286];
                                double l_dof_load_289 = glob->gv__t.coordMap[2 + l_mulRes_286];
                                int32_t t_290 = glob->gv__t.indexMap[l_mulRes_277 + 3];
                                int32_t l_mulRes_291 = 3 * t_290;
                                double l_dof_load_292 = glob->gv__t.coordMap[l_mulRes_291];
                                double l_dof_load_293 = glob->gv__t.coordMap[1 + l_mulRes_291];
                                double l_dof_load_294 = glob->gv__t.coordMap[2 + l_mulRes_291];
                                double t_295 = glob->gv__t.coordMap[l_mulRes_279];
                                double l_r_296 = t_295 * l_mult_269;
                                double l_r_297 = l_dof_load_287 * 0.e0;
                                double l_r_298 = l_dof_load_292 * 0.e0;
                                double l_r_299 = l_r_296 + l_dof_load_282 * l_mult_254 + l_r_297 + l_r_298;
                                double l_r_300 = l_r_296 + l_dof_load_282 * 0.e0;
                                double l_r_301 = l_r_300 + l_dof_load_287 * l_mult_254 + l_r_298;
                                double l_r_302 = l_r_300 + l_r_297 + l_dof_load_292 * l_mult_254;
                                double t_303 = glob->gv__t.coordMap[1 + l_mulRes_279];
                                double l_r_304 = t_303 * l_mult_269;
                                double l_r_305 = l_dof_load_288 * 0.e0;
                                double l_r_306 = l_dof_load_293 * 0.e0;
                                double l_r_307 = l_r_304 + l_dof_load_283 * l_mult_254 + l_r_305 + l_r_306;
                                double l_r_308 = l_r_304 + l_dof_load_283 * 0.e0;
                                double l_r_309 = l_r_308 + l_dof_load_288 * l_mult_254 + l_r_306;
                                double l_r_310 = l_r_308 + l_r_305 + l_dof_load_293 * l_mult_254;
                                double t_311 = glob->gv__t.coordMap[2 + l_mulRes_279];
                                double l_r_312 = t_311 * l_mult_269;
                                double l_r_313 = l_dof_load_289 * 0.e0;
                                double l_r_314 = l_dof_load_294 * 0.e0;
                                double l_r_315 = l_r_312 + l_dof_load_284 * l_mult_254 + l_r_313 + l_r_314;
                                double l_r_316 = l_r_312 + l_dof_load_284 * 0.e0;
                                double l_r_317 = l_r_316 + l_dof_load_289 * l_mult_254 + l_r_314;
                                double l_r_318 = l_r_316 + l_r_313 + l_dof_load_294 * l_mult_254;
                                double l_r_319 = 0.e0 * l_r_299;
                                double l_r_320 = 0.e0 * l_r_307;
                                double l_r_321 = 0.e0 * l_r_315;
                                double l_r_322 = l_r_319 + l_r_320;
                                double l_r_323 = l_r_322 + l_r_321;
                                double l_r_324 = 0.e0 * l_r_301;
                                double l_r_325 = 0.e0 * l_r_309;
                                double l_r_326 = 0.e0 * l_r_317;
                                double l_r_327 = l_r_324 + l_r_325;
                                double l_r_328 = l_r_327 + l_r_326;
                                double l_r_329 = 0.e0 * l_r_302;
                                double l_r_330 = 0.e0 * l_r_310;
                                double l_r_331 = 0.e0 * l_r_318;
                                double l_r_332 = l_r_329 + l_r_330;
                                double l_r_333 = l_r_332 + l_r_331;
                                double l_r_334 = l_r_322 + -0.1e1 * l_r_315;
                                double l_r_335 = l_r_327 + -0.1e1 * l_r_317;
                                double l_r_336 = l_r_332 + -0.1e1 * l_r_318;
                                double l_r_337 = l_r_319 + 0.1e1 * l_r_307 + l_r_321;
                                double l_r_338 = l_r_324 + 0.1e1 * l_r_309 + l_r_326;
                                double l_r_339 = l_r_329 + 0.1e1 * l_r_310 + l_r_331;
                                double l_r_340 = l_r_322 + 0.1e1 * l_r_315;
                                double l_r_341 = l_r_327 + 0.1e1 * l_r_317;
                                double l_r_342 = l_r_332 + 0.1e1 * l_r_318;
                                double l_r_343 = -0.1e1 * l_r_299 + l_r_320 + l_r_321;
                                double l_r_344 = -0.1e1 * l_r_301 + l_r_325 + l_r_326;
                                double l_r_345 = -0.1e1 * l_r_302 + l_r_330 + l_r_331;
                                double l_r_346 = l_r_319 + -0.1e1 * l_r_307 + l_r_321;
                                double l_r_347 = l_r_324 + -0.1e1 * l_r_309 + l_r_326;
                                double l_r_348 = l_r_329 + -0.1e1 * l_r_310 + l_r_331;
                                double l_r_349 = 0.1e1 * l_r_299 + l_r_320 + l_r_321;
                                double l_r_350 = 0.1e1 * l_r_301 + l_r_325 + l_r_326;
                                double l_r_351 = 0.1e1 * l_r_302 + l_r_330 + l_r_331;
                                double l_r_352 = l_r_299 * l_r_328 + l_r_307 * l_r_341 + l_r_315 * l_r_347;
                                double l_r_353 = l_r_299 * l_r_333 + l_r_307 * l_r_342 + l_r_315 * l_r_348;
                                double l_r_354 = l_r_299 * l_r_335 + l_r_307 * l_r_328 + l_r_315 * l_r_350;
                                double l_r_355 = l_r_299 * l_r_336 + l_r_307 * l_r_333 + l_r_315 * l_r_351;
                                double l_r_356 = l_r_299 * l_r_338 + l_r_307 * l_r_344 + l_r_315 * l_r_328;
                                double l_r_357 = l_r_299 * l_r_339 + l_r_307 * l_r_345 + l_r_315 * l_r_333;
                                double l_r_358 = l_r_301 * l_r_323 + l_r_309 * l_r_340 + l_r_317 * l_r_346;
                                double l_r_359 = l_r_301 * l_r_333 + l_r_309 * l_r_342 + l_r_317 * l_r_348;
                                double l_r_360 = l_r_301 * l_r_334 + l_r_309 * l_r_323 + l_r_317 * l_r_349;
                                double l_r_361 = l_r_301 * l_r_336 + l_r_309 * l_r_333 + l_r_317 * l_r_351;
                                double l_r_362 = l_r_301 * l_r_337 + l_r_309 * l_r_343 + l_r_317 * l_r_323;
                                double l_r_363 = l_r_301 * l_r_339 + l_r_309 * l_r_345 + l_r_317 * l_r_333;
                                double l_r_364 = l_r_302 * l_r_323 + l_r_310 * l_r_340 + l_r_318 * l_r_346;
                                double l_r_365 = l_r_302 * l_r_328 + l_r_310 * l_r_341 + l_r_318 * l_r_347;
                                double l_r_366 = l_r_302 * l_r_334 + l_r_310 * l_r_323 + l_r_318 * l_r_349;
                                double l_r_367 = l_r_302 * l_r_335 + l_r_310 * l_r_328 + l_r_318 * l_r_350;
                                double l_r_368 = l_r_302 * l_r_337 + l_r_310 * l_r_343 + l_r_318 * l_r_323;
                                double l_r_369 = l_r_302 * l_r_338 + l_r_310 * l_r_344 + l_r_318 * l_r_328;
                                vec3 v_370 = vcons3(l_r_301, l_r_309, l_r_317);
                                double l_r_371 = 0.e0 * (l_r_299 * l_r_323 + l_r_307 * l_r_340 + l_r_315 * l_r_346);
                                double l_r_372 = 0.e0 * l_r_353;
                                double l_r_373 = 0.e0 * l_r_358;
                                double l_r_374 = 0.e0 * (l_r_301 * l_r_328 + l_r_309 * l_r_341 + l_r_317 * l_r_347);
                                double l_r_375 = 0.e0 * l_r_364;
                                double l_r_376 = 0.e0 * (l_r_302 * l_r_333 + l_r_310 * l_r_342 + l_r_318 * l_r_348);
                                double l_r_377 = l_r_371 + 0.e0 * l_r_352;
                                double l_r_378 = 0.e0 * (l_r_299 * l_r_334 + l_r_307 * l_r_323 + l_r_315 * l_r_349);
                                double l_r_379 = 0.e0 * l_r_355;
                                double l_r_380 = 0.e0 * l_r_360;
                                double l_r_381 = 0.e0 * (l_r_301 * l_r_335 + l_r_309 * l_r_328 + l_r_317 * l_r_350);
                                double l_r_382 = 0.e0 * l_r_366;
                                double l_r_383 = 0.e0 * (l_r_302 * l_r_336 + l_r_310 * l_r_333 + l_r_318 * l_r_351);
                                double l_r_384 = l_r_378 + 0.e0 * l_r_354;
                                double l_r_385 = 0.e0 * (l_r_299 * l_r_337 + l_r_307 * l_r_343 + l_r_315 * l_r_323);
                                double l_r_386 = 0.e0 * l_r_357;
                                double l_r_387 = 0.e0 * l_r_362;
                                double l_r_388 = 0.e0 * (l_r_301 * l_r_338 + l_r_309 * l_r_344 + l_r_317 * l_r_328);
                                double l_r_389 = 0.e0 * l_r_368;
                                double l_r_390 = 0.e0 * (l_r_302 * l_r_339 + l_r_310 * l_r_345 + l_r_318 * l_r_333);
                                double l_r_391 = l_r_385 + 0.e0 * l_r_356;
                                double l_r_392 = 0.e0 * l_r_359;
                                double l_r_393 = 0.e0 * l_r_365;
                                double l_r_394 = 0.e0 * l_r_361;
                                double l_r_395 = 0.e0 * l_r_367;
                                double l_r_396 = 0.e0 * l_r_363;
                                double l_r_397 = 0.e0 * l_r_369;
                                double l_op1_e3_l_36_398 = 0.2e1 * vdot3(vcons3(l_r_299, l_r_307, l_r_315),
                                    vcons3(vdot3(v_370, vcons3(l_r_333, l_r_342, l_r_348)),
                                        vdot3(v_370, vcons3(l_r_336, l_r_333, l_r_351)),
                                        vdot3(v_370, vcons3(l_r_339, l_r_345, l_r_333))));
                                vec3 v_399 = -vcons3(
                                    l_vdot_274 * ((l_r_377 + l_r_372 + l_r_373 + l_r_374 + 0.1e1 * l_r_359 + l_r_375 + -0.1e1 * l_r_365 + l_r_376) / l_op1_e3_l_36_398) + l_vdot_275 * ((l_r_377 + -0.1e1 * l_r_353 + l_r_373 + l_r_374 + l_r_392 + 0.1e1 * l_r_364 + l_r_393 + l_r_376) / l_op1_e3_l_36_398) + l_vdot_276 * ((l_r_371 + 0.1e1 * l_r_352 + l_r_372 + -0.1e1 * l_r_358 + l_r_374 + l_r_392 + l_r_375 + l_r_393 + l_r_376) / l_op1_e3_l_36_398),
                                    l_vdot_274 * ((l_r_384 + l_r_379 + l_r_380 + l_r_381 + 0.1e1 * l_r_361 + l_r_382 + -0.1e1 * l_r_367 + l_r_383) / l_op1_e3_l_36_398) + l_vdot_275 * ((l_r_384 + -0.1e1 * l_r_355 + l_r_380 + l_r_381 + l_r_394 + 0.1e1 * l_r_366 + l_r_395 + l_r_383) / l_op1_e3_l_36_398) + l_vdot_276 * ((l_r_378 + 0.1e1 * l_r_354 + l_r_379 + -0.1e1 * l_r_360 + l_r_381 + l_r_394 + l_r_382 + l_r_395 + l_r_383) / l_op1_e3_l_36_398),
                                    l_vdot_274 * ((l_r_391 + l_r_386 + l_r_387 + l_r_388 + 0.1e1 * l_r_363 + l_r_389 + -0.1e1 * l_r_369 + l_r_390) / l_op1_e3_l_36_398) + l_vdot_275 * ((l_r_391 + -0.1e1 * l_r_357 + l_r_387 + l_r_388 + l_r_396 + 0.1e1 * l_r_368 + l_r_397 + l_r_390) / l_op1_e3_l_36_398) + l_vdot_276 * ((l_r_385 + 0.1e1 * l_r_356 + l_r_386 + -0.1e1 * l_r_362 + l_r_388 + l_r_396 + l_r_389 + l_r_397 + l_r_390) / l_op1_e3_l_36_398));
                                double l_op1_e3_l_57_400 = std::sqrt(vdot3(v_399, v_399));
                                double l_a_401 = 0.1e1 * clamp(0.e0, 0.1e1,
                                    0.13e1 * (0.1e1 - std::abs(l_compositionl_267 - glob->gv_isoval) / (glob->gv_thick * l_op1_e3_l_57_400)));
                                vec3 v_402 = v_399;
                                if (l_a_401 > 0.e0) {
                                    tensor_3_2 l_voxels_412;
                                    double l_imgPos_403 = world2image(glob->gv_I) * l_compositionl_267 + translate(
                                        glob->gv_I);
                                    double l_nd_404 = std::floor(l_imgPos_403);
                                    double l_f_405 = l_imgPos_403 - l_nd_404;
                                    int32_t l_n_406 = std::lround(l_nd_404);
                                    double l__t_407 = std::pow(0.1e1 - l_a_401,
                                        glob->gv_rayStep * std::sqrt(
                                            vdot3(vload3(tensor_ref_3(self->sv_rayVec).addr(0)),
                                                vload3(tensor_ref_3(self->sv_rayVec).addr(0)))) / glob->gv_refStep);
                                    double l_op1_e3_l_19_408 = (self->sv_rayN - glob->gv_camNear) / (glob->gv_camFar - glob->gv_camNear) * (0.7e0 - 0.11e1);
                                    double l_op1_e3_l_20_409 = glob->gv_phongKd * std::max(0.e0,
                                        0.1e1 / l_op1_e3_l_57_400 * vdot3(v_402,
                                            vload3(tensor_ref_3(glob->gv_light).addr(0))));
                                    if (glob->gv_I.inside(l_n_406, 2)) {
                                        int32_t l_offp_410 = 3 * l_n_406;
                                        int32_t l_offp_411 = 3 * (l_n_406 + 1);
                                        l_voxels_412[0] = glob->gv_I[l_offp_410];
                                        l_voxels_412[1] = glob->gv_I[l_offp_411];
                                        l_voxels_412[2] = glob->gv_I[l_offp_410 + 1];
                                        l_voxels_412[3] = glob->gv_I[l_offp_411 + 1];
                                        l_voxels_412[4] = glob->gv_I[l_offp_410 + 2];
                                        l_voxels_412[5] = glob->gv_I[l_offp_411 + 2];
                                    }
                                    else {
                                        int32_t l_offp_413 = 3 * glob->gv_I.clamp(0, l_n_406);
                                        int32_t l_offp_414 = 3 * glob->gv_I.clamp(0, l_n_406 + 1);
                                        l_voxels_412[0] = glob->gv_I[l_offp_413];
                                        l_voxels_412[1] = glob->gv_I[l_offp_414];
                                        l_voxels_412[2] = glob->gv_I[l_offp_413 + 1];
                                        l_voxels_412[3] = glob->gv_I[l_offp_414 + 1];
                                        l_voxels_412[4] = glob->gv_I[l_offp_413 + 2];
                                        l_voxels_412[5] = glob->gv_I[l_offp_414 + 2];
                                    }
                                    vec2 v_415 = vcons2(0.1e1, 0.1e1) + vcons2(l_f_405, l_f_405 - 0.1e1) * vcons2(
                                        -0.1e1, 0.1e1);
                                    double l_op1_e3_l_22_416 = 0.1e1 - l__t_407;
                                    double l_r_417 = self->sv_transp * l_op1_e3_l_22_416 * (0.11e1 + l_op1_e3_l_19_408) * (glob->gv_phongKa + l_op1_e3_l_20_409);
                                    v_418 = vload3(tensor_ref_3(self->sv_rgb).addr(0)) + vcons3(
                                        l_r_417 * vdot2(vload2(l_voxels_412.last(0).addr(0)), v_415),
                                        l_r_417 * vdot2(vload2(l_voxels_412.last(2).addr(0)), v_415),
                                        l_r_417 * vdot2(vload2(l_voxels_412.last(4).addr(0)), v_415));
                                    l_transp_419 = self->sv_transp * (0.1e1 - l_op1_e3_l_22_416);
                                }
                                else {
                                    v_418 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                                    l_transp_419 = self->sv_transp;
                                }
                                v_420 = v_418;
                                l_transp_421 = l_transp_419;
                            }
                            else {
                                tensor_3 _arg_422;
                                vpack3(_arg_422, v_205);
                                wrld->print() << tensor_ref_3(_arg_422) << "\n" << std::flush;
                                wrld->print() << "wait what\?\n" << std::flush;
                                v_420 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                                l_transp_421 = self->sv_transp;
                            }
                            v_423 = v_420;
                            l_transp_424 = l_transp_421;
                        }
                        else {
                            v_423 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                            l_transp_424 = self->sv_transp;
                        }
                        v_425 = v_423;
                        l_transp_426 = l_transp_424;
                    }
                    else {
                        v_425 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                        l_transp_426 = self->sv_transp;
                    }
                    v_427 = v_425;
                    l_transp_428 = l_transp_426;
                }
                else {
                    v_427 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                    l_transp_428 = self->sv_transp;
                }
                v_429 = v_427;
                l_transp_430 = l_transp_428;
            }
            else {
                v_429 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
                l_transp_430 = self->sv_transp;
            }
            v_431 = v_429;
            l_transp_432 = l_transp_430;
        }
        else {
            v_431 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
            l_transp_432 = self->sv_transp;
        }
        v_433 = v_431;
        l_transp_434 = l_transp_432;
    }
    else {
        v_433 = vload3(tensor_ref_3(self->sv_rgb).addr(0));
        l_transp_434 = self->sv_transp;
    }
    if (l_transp_434 < 0.1e-1) {
        self->sv_transp = 0.e0;
        vpack3(self->sv_rgb, v_433);
        return diderot::kStabilize;
    }
    else {
        l_transp_436 = l_transp_434;
    }
    if (self->sv_rayN > glob->gv_camFar) {
        self->sv_transp = l_transp_436;
        vpack3(self->sv_rgb, v_433);
        return diderot::kStabilize;
    }
    self->sv_rayN = self->sv_rayN + glob->gv_rayStep;
    self->sv_transp = l_transp_436;
    vpack3(self->sv_rgb, v_433);
    return diderot::kActive;
}
static void raycast_stabilize (world *wrld, globals *glob, raycast_strand *self)
{
    vec4 v_440;
    double l_a_439 = 0.1e1 - self->sv_transp;
    if (l_a_439 > 0.e0) {
        v_440 = vcons4(tensor_ref_3(self->sv_rgb)[0] / l_a_439, tensor_ref_3(self->sv_rgb)[1] / l_a_439,
            tensor_ref_3(self->sv_rgb)[2] / l_a_439, l_a_439);
    }
    else {
        v_440 = vload4(tensor_ref_4(self->sv_rgba).addr(0));
    }
    if (self->sv_ui == glob->gv_su) {
        if (self->sv_vi == glob->gv_sv) {
            if (glob->gv_debug) {
                tensor_4 _arg_441;
                vpack4(_arg_441, v_440);
                wrld->print() << l_a_439 << tensor_ref_4(_arg_441) << std::flush;
            }
        }
    }
    vpack4(self->sv_rgba, v_440);
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
    int lo_2 = 0;
    int hi_3 = glob->gv_iresV - 1;
    int lo_4 = 0;
    int hi_5 = glob->gv_iresU - 1;
    int32_t base[2] = {lo_2,lo_4,};
    uint32_t size[2] = {static_cast<uint32_t>(hi_3 - lo_2 + 1),static_cast<uint32_t>(hi_5 - lo_4 + 1),};
    if (this->alloc(base, size)) {
        return true;
    }
    uint32_t ix = 0;
    for (int i_vi_443 = lo_2; i_vi_443 <= hi_3; i_vi_443++) {
        for (int i_ui_444 = lo_4; i_ui_444 <= hi_5; i_ui_444++) {
            raycast_init(this->_globals, this->_strands.strand(ix), i_ui_444, i_vi_443);
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


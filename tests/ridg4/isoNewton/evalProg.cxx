/*---------- begin cxx-head.in ----------*/
/*! \file evalProg.cxx
 *
 * Generated from evalProg.diderot.
 *
 * Command: /home/teocollin/gitcode/diderot/bin/diderotc --debug --log --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree --double --namespace=evalProg --target=parallel evalProg.diderot
 * Version: master:2016-07-29
 */
/*---------- end cxx-head.in ----------*/


/**** User mandated includes ****/

#include <spatialindex/capi/sidx_api.h> 

#include <spatialindex/capi/sidx_impl.h> 


/**** Diderot defs and includes ****/

#define DIDEROT_STRAND_HAS_CONSTR
/*---------- begin lib-cxx-incl.in ----------*/
#include "evalProg.h"
#include "diderot/diderot.hxx"

#ifdef DIDEROT_ENABLE_LOGGING
#define IF_LOGGING(...)         __VA_ARGS__
#else
#define IF_LOGGING(...)
#endif

static std::string ProgramName = "evalProg";
/*---------- end lib-cxx-incl.in ----------*/

// ***** Begin synthesized types *****

namespace evalProg {
    typedef double vec4 __attribute__ ((vector_size (32)));
    typedef double vec3 __attribute__ ((vector_size (32)));
    typedef double vec1 __attribute__ ((vector_size (8)));
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
    struct tensor_ref_3_3 : public diderot::tensor_ref<double,9> {
        tensor_ref_3_3 (const double *src);
        tensor_ref_3_3 (struct tensor_3_3 const & ten);
        tensor_ref_3_3 (tensor_ref_3_3 const & ten);
        tensor_ref_3 last (uint32_t i)
        {
            return &this->_data[i];
        }
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
    struct tensor_3_3 : public diderot::tensor<double,9> {
        tensor_3_3 ()
            : diderot::tensor<double,9>()
        { }
        tensor_3_3 (std::initializer_list< double > const & il)
            : diderot::tensor<double,9>(il)
        { }
        tensor_3_3 (const double *src)
            : diderot::tensor<double,9>(src)
        { }
        tensor_3_3 (tensor_3_3 const & ten)
            : diderot::tensor<double,9>(ten._data)
        { }
        ~tensor_3_3 () { }
        tensor_3_3 & operator= (tensor_3_3 const & src);
        tensor_3_3 & operator= (tensor_ref_3_3 const & src);
        tensor_3_3 & operator= (std::initializer_list< double > const & il);
        tensor_3_3 & operator= (const double *src);
        tensor_ref_3 last (uint32_t i)
        {
            return &this->_data[i];
        }
    };
    struct mesh_t {
        int32_t *indexMap;
        double *coordMap;
        int32_t dim;
        int32_t mapDim;
        int32_t numCells;
        void *index;
        int32_t *con;
        mesh_t operator= (std::string file)
        {
            // No something with the nrrd
        }
    };
    struct mesh_pos_mesh_t {
        mesh_t mesh;
        int32_t cell;
        tensor_3 refPos;
        tensor_3 worldPos;
        bool wpc;
        bool valid;
        int32_t face;
        mesh_pos_mesh_t (mesh_t mesh, int32_t cell, tensor_ref_3 refPos, tensor_ref_3 worldPos, bool wpc, bool valid, int32_t face)
            : mesh(mesh), cell(cell), wpc(wpc), valid(valid), face(face)
        {
            mesh = mesh;
            cell = cell;
            refPos = refPos;
            worldPos = worldPos;
            wpc = wpc;
            valid = valid;
            face = face;
        }
        mesh_pos_mesh_t ()
            : cell(-1), face(-1), wpc(false), valid(false)
        { }
        mesh_pos_mesh_t (mesh_t mesh)
            : cell(-1), face(-1), wpc(false), valid(false), mesh(mesh)
        { }
    };
    struct fns_t {
        int32_t *indexMap;
        mesh_t mesh;
        fns_t operator= (std::string file)
        {
            // No something with the nrrd
        }
        fns_t loadFem (mesh_t mesh)
        {
            fns_t space = *this;
            space.mesh = mesh;
            return space;
        }
    };
    struct mesh_cell_mesh_t {
        int32_t cell;
        mesh_t mesh;
        std::ostream & operator<< (std::ostream & os)
        {
            return os << this->cell;
        }
    };
    struct func_t {
        double *coordMap;
        fns_t space;
        func_t operator= (std::string file)
        {
            // No something with the nrrd
        }
        func_t loadFem (fns_t space)
        {
            func_t func = *this;
            func.space = space;
            return func;
        }
    };
    struct func_cell_func_t {
        int32_t cell;
        func_t func;
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
    inline tensor_ref_7::tensor_ref_7 (const double *src)
        : diderot::tensor_ref<double,7>(src)
    { }
    inline tensor_ref_7::tensor_ref_7 (struct tensor_7 const & ten)
        : diderot::tensor_ref<double,7>(ten._data)
    { }
    inline tensor_ref_7::tensor_ref_7 (tensor_ref_7 const & ten)
        : diderot::tensor_ref<double,7>(ten._data)
    { }
    inline tensor_ref_3_3::tensor_ref_3_3 (const double *src)
        : diderot::tensor_ref<double,9>(src)
    { }
    inline tensor_ref_3_3::tensor_ref_3_3 (struct tensor_3_3 const & ten)
        : diderot::tensor_ref<double,9>(ten._data)
    { }
    inline tensor_ref_3_3::tensor_ref_3_3 (tensor_ref_3_3 const & ten)
        : diderot::tensor_ref<double,9>(ten._data)
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
    inline tensor_3_3 & tensor_3_3::operator= (std::initializer_list< double > const & il)
    {
        this->copy(il);
        return *this;
    }
    inline tensor_3_3 & tensor_3_3::operator= (const double *src)
    {
        this->copy(src);
        return *this;
    }
    std::ostream & operator<< (std::ostream & os, const mesh_pos_mesh_t &  pos)
    {
        return os << pos.cell;
    }
    mesh_pos_mesh_t allBuild (mesh_t mesh, int32_t cell, tensor_ref_3 refPos, tensor_ref_3 worldPos, bool wpc, bool valid)
    {
        mesh_pos_mesh_t pos;
        pos.mesh = mesh;
        pos.cell = cell;
        pos.refPos = refPos._data;
        pos.worldPos = worldPos._data;
        pos.wpc = wpc;
        pos.valid = valid;
        return pos;
    }
    mesh_pos_mesh_t allBuild (mesh_t mesh, int32_t cell, tensor_ref_3 refPos, tensor_ref_3 worldPos, bool wpc, bool valid, int32_t face)
    {
        mesh_pos_mesh_t pos;
        pos.mesh = mesh;
        pos.cell = cell;
        pos.face = face;
        pos.refPos = refPos._data;
        pos.worldPos = worldPos._data;
        pos.wpc = wpc;
        pos.valid = valid;
        return pos;
    }
    mesh_pos_mesh_t invalidBuild (mesh_t mesh)
    {
        mesh_pos_mesh_t pos;
        pos.mesh = mesh;
        return pos;
    }
    mesh_pos_mesh_t invalidBuild (mesh_t mesh, tensor_ref_3 refPos)
    {
        mesh_pos_mesh_t pos;
        pos.mesh = mesh;
        pos.refPos = refPos;
        return pos;
    }
    mesh_pos_mesh_t refBuild (mesh_t mesh, int32_t cell, tensor_ref_3 refPos)
    {
        mesh_pos_mesh_t pos;
        pos.mesh = mesh;
        pos.cell = cell;
        pos.refPos = refPos;
        pos.wpc = false;
        pos.valid = true;
        return pos;
    }
    mesh_pos_mesh_t refBuild (mesh_t mesh, int32_t cell, tensor_ref_3 refPos, int32_t face)
    {
        mesh_pos_mesh_t pos;
        pos.mesh = mesh;
        pos.cell = cell;
        pos.refPos = refPos;
        pos.wpc = false;
        pos.valid = true;
        pos.face = face;
        return pos;
    }
    char * copy_to  (const mesh_pos_mesh_t dumb, char *cp)
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
    mesh_pos_mesh_t makeFem (mesh_t mesh, tensor_ref_7 data)
    {
        double cell = data[0];
        tensor_ref_3 ref = &data[1];
        tensor_ref_3 world = &data[4];
        int32_t cellInt = static_cast<int32_t>(cell);
        bool test = 0 < cellInt;
        return allBuild(mesh, cellInt, ref, world, test, test);
    }
    std::ostream & operator<< (std::ostream & os, const mesh_cell_mesh_t &  cell)
    {
        return os << cell.cell;
    }
    mesh_cell_mesh_t makeFem (mesh_t mesh, int32_t cellInt)
    {
        mesh_cell_mesh_t cell;
        cell.cell = cellInt;
        cell.mesh = mesh;
        return cell;
    }
    char * copy_to (const mesh_cell_mesh_t dumb, char *cp)
    {
        size_t nbytes = sizeof(int32_t);
        (std::memcpy)(cp, &(dumb.cell), nbytes);
        return cp + nbytes;
    }
    std::ostream & operator<< (std::ostream & os, const func_cell_func_t &  cell)
    {
        return os << cell.cell;
    }
    func_cell_func_t makeFem (func_t func, int32_t cellInt)
    {
        func_cell_func_t cell;
        cell.cell = cellInt;
        cell.func = func;
        return cell;
    }
    char * copy_to (const func_cell_func_t dumb, char *cp)
    {
        size_t nbytes = sizeof(int32_t);
        (std::memcpy)(cp, &(dumb.cell), nbytes);
        return cp + nbytes;
    }
} // namespace evalProg
namespace diderot {
    template <>
    struct dynseq_traits<evalProg::mesh_cell_mesh_t> {
        using value_type = evalProg::mesh_cell_mesh_t;
        using base_type = int32_t;
        static const __details::load_fn_ptr<base_type> *load_fn_tbl;
        static const uint32_t values_per_elem = 1;
    };
    const __details::load_fn_ptr< dynseq_traits< evalProg::mesh_cell_mesh_t >::base_type > *dynseq_traits< evalProg::mesh_cell_mesh_t >::load_fn_tbl = nrrdILoad;
    template <>
    struct dynseq_traits<evalProg::tensor_3> {
        using value_type = evalProg::tensor_3;
        using base_type = double;
        static const __details::load_fn_ptr<base_type> *load_fn_tbl;
        static const uint32_t values_per_elem = 3;
    };
    const __details::load_fn_ptr< dynseq_traits< evalProg::tensor_3 >::base_type > *dynseq_traits< evalProg::tensor_3 >::load_fn_tbl = nrrdDLoad;
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
namespace evalProg {

static std::string ProgramName = "evalProg";

struct world;
struct normal_strand;
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
    bool gv_meshData;
    bool gv_0space0396_intermedateGlobal;
    bool gv_0data0398_intermedateGlobal;
    bool gv_ipos;
} defined_inputs;
struct globals {
    mesh_t gv_meshData;
    fns_t gv_0space0396_intermedateGlobal;
    func_t gv_0data0398_intermedateGlobal;
    diderot::dynseq< tensor_3 > gv_ipos;
    func_t gv_data;
    ~globals () { }
};
struct normal_strand {
    tensor_3 sv_normal;
    mesh_pos_mesh_t sv_pos0;
    tensor_3 sv_xp;
};
/*---------- begin par-sarr.in ----------*/
// forward declaration of worker_cache type
struct worker_cache;
// forward declarations of strand methods
#ifdef DIDEROT_HAS_START_METHOD
static diderot::strand_status normal_start (normal_strand *self);
#endif // DIDEROT_HAS_START_METHOD
static diderot::strand_status normal_update (world *wrld, globals *glob, normal_strand *self);
#ifdef DIDEROT_HAS_STABILIZE_METHOD
static void normal_stabilize (normal_strand *self);
#endif // DIDEROT_HAS_STABILIZE_METHOD

// strand_array for PARALLEL_TARGET/NO BSP/SINGLE STATE/DIRECT ACCESS
//
struct strand_array {
    typedef normal_strand strand_t;
    typedef uint32_t index_t;
    typedef index_t sid_t;              // strand ID (index into strand-state storage)

    // scheduling block of strands
    //
    struct CACHE_ALIGN sched_block {
        index_t         _start;         // first index in block
        index_t         _stop;          // last index in block + 1
        uint32_t        _nStable;       // number of stable strands in the block
        uint32_t        _nDead;         // number of dead strands in the block

      // return the number of strands in the block
        uint32_t num_strands () const { return this->_stop - this->_start; }
      // return the number of active strands in the block
        uint32_t num_active () const
        {
#ifdef DIDEROT_HAS_STRAND_DIE
            return this->num_strands() - (this->_nStable + this->_nDead);
#else
            return this->num_strands() - this->_nStable;
#endif
        }
    };

    uint8_t             *_status;       // the array of status information for the strands
    char                *_storage;      // points to array of normal_strand structs
    sched_block         *_schedBlks;    // blocks of strands for parallel scheduling
    uint32_t            _nItems;        // number of items in the _storage and _status arrays
    uint32_t            _nStable;       // global number of stable strands
    uint32_t            _nActive;       // global number of active strands
    uint32_t            _nFresh;        // number of fresh strands (new strands from create_strands)
    uint32_t            _nSchedBlks;    // number of scheduling blocks
    uint32_t            _schedBlkSz;    // size of scheduling blocks
    atomic_uint32_t     _nextSchedBlk CACHE_ALIGN;
                                        // next block to schedule
    std::vector<worker_cache *> _workers;

    strand_array ()
        : _status(nullptr), _storage(nullptr), _schedBlks(nullptr), _nItems(0),
          _nStable(0), _nActive(0), _nFresh(0), _nSchedBlks(0), _schedBlkSz(0), _nextSchedBlk(0)
    { }
    ~strand_array ();

    uint32_t in_state_index () const { return 0; /* dummy */ }

    uint32_t num_active () const { return this->_nActive; }
    uint32_t num_stable () const { return this->_nStable; }
    uint32_t num_alive () const { return this->_nActive+this->_nStable; }

  // return the ID of a strand, which is just the value of the argument
    sid_t id (index_t ix) const
    {
        assert (ix < this->_nItems);
        return ix;
    }
  // return a pointer to the strand with the given ID
    normal_strand *id_to_strand (sid_t id) const
    {
        assert (id < this->_nItems);
        return reinterpret_cast<normal_strand *>(this->_storage + id * sizeof(normal_strand));
    }

  // return a strand's status
    diderot::strand_status status (index_t ix) const
    {
        return static_cast<diderot::strand_status>(this->_status[ix]);
    }
  // return a pointer to the given strand
    normal_strand *strand (index_t ix) const
    {
        return this->id_to_strand(this->id(ix));
    }
  // return a pointer to the local state of strand ix
    normal_strand *local_state (index_t ix) const
    {
        return this->strand(ix);
    }
  // return a pointer to the local state of strand with the given ID
    normal_strand *id_to_local_state (sid_t id) const
    {
        return this->id_to_strand(id);
    }

  // set the scheduling block size based on the number of workers and the number of
  // strands.  This should be called before alloc.
    void set_block_size (uint32_t nWorkers, uint32_t nStrands)
    {
        this->_schedBlkSz = diderot::sched_block_size (nWorkers, nStrands);
    }

  // allocate space for nItems organized into blkSz sized blocks of strands
    bool alloc (uint32_t nItems);

  // initialize the first nStrands locations as new active strands
    void create_strands (uint32_t nStrands);

  // swap in and out states (NOP for this version)
    void swap () { }

  // invoke strand's stabilize method (single-thread version)
    index_t strand_stabilize (index_t ix)
    {
#ifdef DIDEROT_HAS_STABILIZE_METHOD
        normal_stabilize (this->strand(ix));
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

  // mark the given strand as dead (single-thread version)
    index_t kill (index_t ix)
    {
        this->_status[ix] = diderot::kDead;
        this->_nActive--;
      // skip to next active strand
        do {
            ix++;
        } while ((ix < this->_nItems) && notActiveSts(this->status(ix)));
        return ix;
    }

  // prepare to run the workers
    void prepare_run ()
    {
        this->_nextSchedBlk = 0;
    }

#ifdef DIDEROT_BSP
  // finish the local-phase of a superstep; note that this function is only used
  // when BSP is enabled.
    bool finish_step ();
#endif

  // finish a kill_all operation (NOP)
    void finish_kill_all () { }

  // finish a stabilize_all operation (NOP)
    void finish_stabilize_all () { }

  // iterator over all alive strands (single-threaded version)
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

  // iterator over all active strands (single-threaded version)
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
        this->strand(ix)->~normal_strand();
    }
    if (this->_status != nullptr) std::free (this->_status);
    if (this->_storage != nullptr) std::free (this->_storage);
    if (this->_schedBlks != nullptr) std::free (this->_schedBlks);
}

bool strand_array::alloc (uint32_t nItems)
{
    if (this->_schedBlkSz == 0) {
        std::cerr << "Internal error: strand_array block size is 0\n";
        return true;
    }
    this->_storage = static_cast<char *>(std::malloc (nItems * sizeof(normal_strand)));
    if (this->_storage == nullptr) {
        return true;
    }
    this->_status = static_cast<uint8_t *>(std::malloc (nItems * sizeof(uint8_t)));
    if (this->_status == nullptr) {
        std::free (this->_storage);
        return true;
    }
    this->_nSchedBlks = (nItems + this->_schedBlkSz - 1) / this->_schedBlkSz;
    this->_schedBlks =
        static_cast<sched_block *>(std::malloc (this->_nSchedBlks * sizeof(sched_block)));
    if (this->_schedBlks == nullptr) {
        std::free (this->_storage);
        std::free (this->_status);
        return true;
    }
    this->_nItems = nItems;
    this->_nActive = 0;
    this->_nStable = 0;
    this->_nFresh = 0;
    return false;
}

void strand_array::create_strands (uint32_t nStrands)
{
    assert (this->_nActive == 0);
    assert (this->_nItems == nStrands);
    for (uint32_t ix = 0;  ix < nStrands;  ix++) {
#ifdef DIDEROT_HAS_START_METHOD
        this->_status[ix] = diderot::kNew;
#else
        this->_status[ix] = diderot::kActive;
#endif
        new(this->strand(ix)) normal_strand;
    }
    this->_nActive = nStrands;
    this->_nFresh = nStrands;
  // initialize the scheduling blocks
    for (uint32_t ix = 0, i = 0;  i < this->_nSchedBlks;  i++) {
        this->_schedBlks[i]._start = ix;
        ix += this->_schedBlkSz;
        if (ix < nStrands) {
            this->_schedBlks[i]._stop = ix;
        }
        else {
            this->_schedBlks[i]._stop = nStrands;
        }
        this->_schedBlks[i]._nDead = 0;
        this->_schedBlks[i]._nStable = 0;
    }
}

// a local copy of strand state for workers
struct worker_cache {
    typedef strand_array::strand_t strand_t;
    typedef strand_array::index_t index_t;
    typedef strand_array::sid_t sid_t;
    typedef strand_array::sched_block sched_block;

    uint8_t             *_status;       // the array of status information for the strands
    char                *_storage;      // points to array of normal_strand structs
    sched_block         *_schedBlks;    // blocks of strands for parallel scheduling
    atomic_uint32_t     *_nextBlkPtr;   // pointer to _nextSchedBlk
    uint32_t            _nStabilizing;  // count of strands run by this worker that stabilized in
                                        // the current superstep
#ifdef DIDEROT_HAS_STRAND_DIE
    uint32_t            _nDying;        // count of strands run by this worker that died in
                                        // the current superstep
#endif
    uint32_t            _nSchedBlks;    // number of scheduling blocks
    uint32_t            _schedBlkSz;    // size of scheduling blocks
#ifndef NDEBUG
    uint32_t        _nItems;            // number of items in the _storage and _status arrays
#endif

    void init (strand_array &sarr)
    {
        this->_status = sarr._status;
        this->_storage = sarr._storage;
        this->_schedBlks = sarr._schedBlks;
        this->_nextBlkPtr = &sarr._nextSchedBlk;
        this->_nSchedBlks = sarr._nSchedBlks;
        this->_schedBlkSz = sarr._schedBlkSz;
#ifndef NDEBUG
        this->_nItems = sarr._nItems;
#endif
        sarr._workers.push_back (this);
    }

  // refresh those parts of the cache that might change between steps
    void refresh ()
    {
        // this target does not support dynamic strands, so nothing can change
    }

  // return the ID of a strand, which is the value of the _idx array
    sid_t id (index_t ix) const
    {
        assert (ix < this->_nItems);
        return ix;
    }
  // return a pointer to the strand with the given ID
    normal_strand *id_to_strand (sid_t id) const
    {
        return reinterpret_cast<normal_strand *>(this->_storage + id * sizeof(normal_strand));
    }
  // return a strand's status
    diderot::strand_status status (index_t ix) const
    {
        return static_cast<diderot::strand_status>(this->_status[ix]);
    }
  // return a pointer to the given strand
    normal_strand *strand (index_t ix) const
    {
        return this->id_to_strand(this->id(ix));
    }
  // return a pointer to the local state of strand ix
    normal_strand *local_state (index_t ix) const
    {
        return this->strand(ix);
    }
  // return a pointer to the local state of strand with the given ID
    normal_strand *id_to_local_state (sid_t id) const
    {
        return this->id_to_strand(id);
    }

#ifdef DIDEROT_HAS_START_METHOD
  // invoke strand's start method
    diderot::strand_status strand_start (index_t ix)
    {
        return normal_start(this->strand(ix));
    }

    void run_start_methods (sched_block *bp);
#endif // DIDEROT_HAS_START_METHOD

  // invoke strand's update method
    diderot::strand_status strand_update (world *wrld, globals *glob, index_t ix)
    {
        return normal_update(wrld, glob, this->strand(ix));
    }

  // invoke strand's stabilize method (multithread version)
    index_t strand_stabilize (sched_block *bp, index_t ix)
    {
#ifdef DIDEROT_HAS_STABILIZE_METHOD
        normal_stabilize (this->strand(ix));
#endif // DIDEROT_HAS_STABILIZE_METHOD
        this->_status[ix] = diderot::kStable;
        bp->_nStable++;
      // skip to next active strand
        do {
            ix++;
        } while ((ix < bp->_stop) && notActiveSts(this->status(ix)));
        return ix;
    }

  // mark the given strand as dead (multithread version)
    index_t kill (sched_block *bp, index_t ix)
    {
        this->_status[ix] = diderot::kDead;
        bp->_nDead++;
      // skip to next active strand
        do {
            ix++;
        } while ((ix < bp->_stop) && notActiveSts(this->status(ix)));
        return ix;
    }

  // iterator over alive strands in a scheduling block
    index_t begin_alive (const sched_block *bp) const
    {
        index_t ix = bp->_start;
#ifdef DIDEROT_HAS_STRAND_DIE
        while ((ix < bp->_stop) && notAliveSts(this->status(ix))) {
            ix++;
        }
#endif
        return ix;
    }
    index_t end_alive (const sched_block *bp) const { return bp->_stop; }
    index_t next_alive (const sched_block *bp, index_t &ix) const
    {
#ifdef DIDEROT_HAS_STRAND_DIE
        do {
            ix++;
        } while ((ix < bp->_stop) && notAliveSts(this->status(ix)));
#endif
        return ix;
    }

  // iterator over active strands in a scheduling block
    index_t begin_active (const sched_block *bp) const
    {
        index_t ix = bp->_start;
        while ((ix < bp->_stop) && notActiveSts(this->status(ix))) {
            ix++;
        }
        return ix;
    }
    index_t end_active (const sched_block *bp) const { return bp->_stop; }
    index_t next_active (const sched_block *bp, index_t &ix) const
    {
        do {
            ix++;
        } while ((ix < bp->_stop) && notActiveSts(this->status(ix)));
        return ix;
    }

  // iterator over fresh strands in a scheduling block
    index_t begin_fresh (const sched_block *bp) const
    {
        index_t ix = bp->_start;
        while ((ix < bp->_stop) && (this->status(ix) != diderot::kNew)) {
            ix++;
        }
        return ix;
    }
    index_t end_fresh (const sched_block *bp) const { return bp->_stop; }
    index_t next_fresh (const sched_block *bp, index_t &ix) const
    {
        do {
            ix++;
        } while ((ix < bp->_stop) && (this->status(ix) != diderot::kNew));
        return ix;
    }

  // swap in and out states (NOP for this version)
    void swap () { }

  // get a block of strands
    sched_block *get_block ();

}; // struct worker_cache

strand_array::sched_block *worker_cache::get_block ()
{
    do {
        uint32_t blkId = this->_nextBlkPtr->fetch_add(1);
        if (blkId < this->_nSchedBlks) {
            strand_array::sched_block *bp = &this->_schedBlks[blkId];
            if (bp->num_active() > 0) {
                return bp;
            } // else skip stable block
        }
        else {  // no more blocks
            return nullptr;
        }
    } while (true);

}

#ifdef DIDEROT_BSP
// finish the update phase of a superstep.  Return true if there are any dead strands.
bool strand_array::finish_step ()
{
    int32_t nStabilizing = 0;
#ifdef DIDEROT_HAS_STRAND_DIE
    int32_t nDying = 0;
#endif

    for (auto it = this->_workers.begin();  it != this->_workers.end();  ++it) {
        worker_cache *wp = *it;
        nStabilizing += wp->_nStabilizing;
#ifdef DIDEROT_HAS_STRAND_DIE
        nDying += wp->_nDying;
#endif
    }

#ifdef DIDEROT_HAS_STRAND_DIE
    if (nDying > 0) {
      /* FIXME: compact dead strands */
/*
      // check to see if we need to compact dead strands?
        if ((this->_nStrands - this->_nActive) / this->_schedBlkSz > ??) {
        }
*/
    }
#endif

  // reset scheduler for next superstep
    this->_nextSchedBlk = 0;

  // update global count of stable strands
    this->_nStable += nStabilizing;
  // update global count of active strands
#ifdef DIDEROT_HAS_STRAND_DIE
    this->_nActive -= (nStabilizing + nDying);

    return (nDying > 0);
#else
    this->_nActive -= nStabilizing;

    return false;
#endif

}
#endif // DIDEROT_BSP
/*---------- end par-sarr.in ----------*/

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

diderot::dynseq< int32_t > mesh_geom_mesh_t (void *index, mesh_t *mesh, const double *data)
{
    //diderot::dynseq<int32_t> myFunction(void * index, mesh_t * mesh, double * data)

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
static std::ostream& operator<< (std::ostream& outs, tensor_ref_3 const & ten)
{
    return outs << "[" << ten._data[0] << "," << ten._data[1] << "," << ten._data[2] << "]";
}
inline vec3 vload3 (const double *vp)
{
    return __extension__ (vec3){vp[0], vp[1], vp[2], 0.e0};
}
inline vec3 vcons3 (double r0, double r1, double r2)
{
    return __extension__ (vec3){r0, r1, r2, 0.e0};
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
inline vec3 vscale3 (double s, vec3 v)
{
    return __extension__ (vec3){s, s, s, 0.e0} * v;
}
/*---------- begin eigenvecs3x3.in ----------*/
inline double dot3 (tensor_ref_3 const & a, tensor_ref_3 const & b)
{
    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

inline void cross3 (tensor_ref_3 const & v1, tensor_ref_3 const & v2, tensor_3 & res)
{
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

inline void nullspace1 (
    tensor_ref_3 const & r0,
    tensor_ref_3 const & r1,
    tensor_ref_3 const & r2,
    tensor_3 & res)
{
    tensor_3 crs;

    cross3(r0, r1, res);
    cross3(r1, r2, crs);

  /* ret += crs or ret -= crs; whichever makes res longer */
    if (dot3(res, crs) > 0.0) {
	res[0] += crs[0];
	res[1] += crs[1];
	res[2] += crs[2];
    } else {
	res[0] -= crs[0];
	res[1] -= crs[1];
	res[2] -= crs[2];
    }

    cross3(r0, r2, crs);
  /* ret += crs or ret -= crs; whichever makes res longer */
    if (dot3(res, crs) > 0.0) {
	res[0] += crs[0];
	res[1] += crs[1];
	res[2] += crs[2];
    } else {
	res[0] -= crs[0];
	res[1] -= crs[1];
	res[2] -= crs[2];
    }
}

/*
** All vectors are in the same 1D space, we have to find two
** mutually vectors perpendicular to that span
*/
static void nullspace2 (
    tensor_ref_3 const & r0,
    tensor_ref_3 const & r1,
    tensor_ref_3 const & r2,
    tensor_3 *rets)  // will point to either evec[0] or evec[1]
{
    tensor_3 sqr, sum;
    int idx;

    sum = r0;
    if (dot3(sum, r1) > 0) {
	sum[0] += r1[0];
	sum[1] += r1[1];
	sum[2] += r1[2];
    } else {
	sum[0] -= r1[0];
	sum[1] -= r1[1];
	sum[2] -= r1[2];
    }
    if (dot3(sum, r2) > 0) {
	sum[0] += r2[0];
	sum[1] += r2[1];
	sum[2] += r2[2];
    } else {
	sum[0] -= r2[0];
	sum[1] -= r2[1];
	sum[2] -= r2[2];
    }
  // find largest component, to get most stable expression for a perpendicular vector
    sqr[0] = sum[0]*sum[0];
    sqr[1] = sum[1]*sum[1];
    sqr[2] = sum[2]*sum[2];
    idx = 0;
    if (sqr[0] < sqr[1]) {
	idx = 1;
    }
    if (sqr[idx] < sqr[2]) {
	idx = 2;
    }

    if (0 == idx) {
      rets[0] = {sum[1] - sum[2], -sum[0], sum[0]};
    } else if (1 == idx) {
      rets[0] = {-sum[1], sum[0] - sum[2], sum[1]};
    } else {
      rets[0] = {-sum[2], sum[2], sum[0] - sum[1]};
    }

    cross3(rets[0], sum, rets[1]);
    return;
}

static int eigenvecs (tensor_ref_3_3 const &mat, diderot::array<double, 3> &eval, diderot::array<tensor_3,3> &evec)
{
    double len, dot;
    int roots;

  /* copy the given matrix elements */
    double M00 = mat[0];
    double M01 = mat[1];
    double M02 = mat[2];
    double M11 = mat[4];
    double M12 = mat[5];
    double M22 = mat[8];

  /*
  ** subtract out the eigenvalue mean (will add back to evals later);
  ** helps with numerical stability
  */
    double mean = (M00 + M11 + M22) / double(3);
    M00 -= mean;
    M11 -= mean;
    M22 -= mean;

  /*
  ** divide out L2 norm of eigenvalues (will multiply back later);
  ** this too seems to help with stability
  */
    double norm = std::sqrt(M00*M00 + 2*M01*M01 + 2*M02*M02 + M11*M11 + 2*M12*M12 + M22*M22);
    double rnorm = (norm > diderot::__details::EPSILON) ? double(1) / norm : double(1);
    M00 *= rnorm;
    M01 *= rnorm;
    M02 *= rnorm;
    M11 *= rnorm;
    M12 *= rnorm;
    M22 *= rnorm;

  /* this code is a mix of prior Teem code and ideas from Eberly's
   * "Eigensystems for 3 x 3 Symmetric Matrices (Revisited)"
   */
    double Q = (M01*M01 + M02*M02 + M12*M12 - M00*M11 - M00*M22 - M11*M22)/double(3);
    double QQQ = Q*Q*Q;
    double R = double(0.5)*(M00*M11*M22 + M02*(2*M01*M12 - M02*M11) - M00*M12*M12 - M01*M01*M22);
    double D = QQQ - R*R;
    if (D > diderot::__details::EPSILON) {
      /* three distinct roots- this is the most common case */
        double theta = std::atan2(std::sqrt(D), R) / double(3);
        double mm = std::sqrt(Q);
        double ss = std::sin(theta);
        double cc = std::cos(theta);
        eval[0] = 2*mm*cc;
        eval[1] = mm*(-cc + std::sqrt(3.0)*ss);
        eval[2] = mm*(-cc - std::sqrt(3.0)*ss);
        roots = diderot::__details::ROOT_THREE;
    }
  /* else D is near enough to zero */
    else if (std::abs(R) > diderot::__details::EPSILON) {
      /* one double root and one single root */
        double U = std::cbrt(R); /* cube root function */
        if (U > 0) {
            eval[0] = 2*U;
            eval[1] = -U;
            eval[2] = -U;
        } else {
            eval[0] = -U;
            eval[1] = -U;
            eval[2] = 2*U;
        }
        roots = diderot::__details::ROOT_SINGLE_DOUBLE;
    }
    else {
      /* a triple root! */
        eval[0] = eval[1] = eval[2] = 0.0;
        roots = diderot::__details::ROOT_TRIPLE;
    }
/* END #include "teigen-evals-A.c" */

    tensor_3 ev = tensor_3 { eval[0], eval[1], eval[2] };
    if (diderot::__details::ROOT_THREE == roots) {
         nullspace1 (
	    tensor_3 { M00 - eval[0], M01, M02 },
	    tensor_3 { M01, M11 - eval[0], M12 },
	    tensor_3 { M02, M12, M22 - eval[0] },
	    evec[0]);
        nullspace1 (
	    tensor_3 { M00 - eval[1], M01, M02 },
	    tensor_3 { M01, M11 - eval[1], M12 },
	    tensor_3 { M02, M12, M22 - eval[1] },
	    evec[1]);
        nullspace1 (
	    tensor_3 { M00 - eval[2], M01, M02 },
	    tensor_3 { M01, M11 - eval[2], M12 },
	    tensor_3 { M02, M12, M22 - eval[2] },
	    evec[2]);
    }
    else if (diderot::__details::ROOT_SINGLE_DOUBLE == roots) {
        if (eval[1] == eval[2]) {
          /* one big (eval[0]) , two small (eval[1,2]) */
            nullspace1 (
		tensor_3 { M00 - eval[0], M01, M02 },
		tensor_3 { M01, M11 - eval[0], M12 },
		tensor_3 { M02, M12, M22 - eval[0] },
		evec[0]);
            nullspace2 (
		tensor_3 { M00 - eval[1], M01, M02 },
		tensor_3 { M01, M11 - eval[1], M12 },
		tensor_3 { M02, M12, M22 - eval[1] },
		&evec[1]);
        }
        else {
          /* two big (eval[0,1]), one small (eval[2]) */
            nullspace2 (
		tensor_3 { M00 - eval[0], M01, M02 },
		tensor_3 { M01, M11 - eval[0], M12 },
		tensor_3 { M02, M12, M22 - eval[0] },
		&evec[0]);
            nullspace1 (
		tensor_3 { M00 - eval[2], M01, M02 },
		tensor_3 { M01, M11 - eval[2], M12 },
		tensor_3 { M02, M12, M22 - eval[2] },
		evec[2]);
        }
    }
    else {
      /* ROOT_TRIPLE == roots; use any basis for eigenvectors */
        evec[0] = tensor_3 { 1, 0, 0 };
        evec[1] = tensor_3 { 0, 1, 0 };
        evec[2] = tensor_3 { 0, 0, 1 };
    }
  /* we always make sure it's really orthonormal; keeping fixed the
   * eigenvector associated with the largest-magnitude eigenvalue
   */
    if (std::abs(eval[0]) > std::abs(eval[2])) {
      /* normalize evec[0] but don't move it */
        diderot::__details::normalize3(evec[0]._data);
      // compute evec[1] -= scale3(dot3(evec[1], evec[0]), evec[0]);
	double s = dot3(evec[1], evec[0]);
	evec[1][0] -= s*evec[0][0];
	evec[1][1] -= s*evec[0][1];
	evec[1][2] -= s*evec[0][2];
        diderot::__details::normalize3(evec[1]._data);
        cross3(evec[0], evec[1], evec[2]);
    }
    else {
      /* normalize evec[2] but don't move it */
        diderot::__details::normalize3(evec[2]._data);
      // compute evec[1] -= scale3(dot3(evec[1], evec[2]), evec[2]);
	double s = dot3(evec[1], evec[2]);
	evec[1][0] -= s*evec[2][0];
	evec[1][1] -= s*evec[2][1];
	evec[1][2] -= s*evec[2][2];
        diderot::__details::normalize3(evec[1]._data);
        cross3(evec[1], evec[2], evec[0]);
    }
    /* note that the right-handedness check has been folded into
       the code above to enforce orthogonality.  Indeed, some work
       could be removed by never really bothering to find all three
       eigenvectors; just find two and then use the cross-product.
       The next iteration of the code will do that */

  /* multiply back by eigenvalue L2 norm */
    eval[0] /= rnorm;
    eval[1] /= rnorm;
    eval[2] /= rnorm;

  /* add back in the eigenvalue mean */
    eval[0] += mean;
    eval[1] += mean;
    eval[2] += mean;

    return roots;
}
/*---------- end eigenvecs3x3.in ----------*/

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

extern "C" bool evalProg_input_set_meshData (evalProg_world_t *cWrld, void *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_meshData = true;
    std::memcpy(&wrld->_globals->gv_meshData, v, sizeof(mesh_t));
    return false;
}
extern "C" bool evalProg_input_set_space (evalProg_world_t *cWrld, void *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_0space0396_intermedateGlobal = true;
    std::memcpy(&wrld->_globals->gv_0space0396_intermedateGlobal, v, sizeof(fns_t));
    return false;
}
extern "C" bool evalProg_input_set_data (evalProg_world_t *cWrld, void *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_0data0398_intermedateGlobal = true;
    std::memcpy(&wrld->_globals->gv_0data0398_intermedateGlobal, v, sizeof(func_t));
    return false;
}
extern "C" bool evalProg_input_set_by_name_ipos (evalProg_world_t *cWrld, const char *s)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_ipos = true;
    if (wrld->_globals->gv_ipos.load(wrld, s)) {
        return true;
    }
    return false;
}
extern "C" bool evalProg_input_set_ipos (evalProg_world_t *cWrld, Nrrd *nin)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_ipos = true;
    if (wrld->_globals->gv_ipos.load(wrld, nin)) {
        return true;
    }
    return false;
}
static bool check_defined (world *wrld)
{
    if (!wrld->_definedInp.gv_meshData) {
        biffMsgAdd(wrld->_errors, "undefined input \"meshData\"\n");
        return true;
    }
    if (!wrld->_definedInp.gv_0space0396_intermedateGlobal) {
        biffMsgAdd(wrld->_errors, "undefined input \"space\"\n");
        return true;
    }
    if (!wrld->_definedInp.gv_0data0398_intermedateGlobal) {
        biffMsgAdd(wrld->_errors, "undefined input \"data\"\n");
        return true;
    }
    if (!wrld->_definedInp.gv_ipos) {
        biffMsgAdd(wrld->_errors, "undefined input \"ipos\"\n");
        return true;
    }
    return false;
}
static void init_defined_inputs (world *wrld)
{
    wrld->_definedInp.gv_meshData = false;
    wrld->_definedInp.gv_0space0396_intermedateGlobal = false;
    wrld->_definedInp.gv_0data0398_intermedateGlobal = false;
    wrld->_definedInp.gv_ipos = false;
}
static void init_defaults (globals *glob)
{
}
mesh_pos_mesh_t fn_findPos (mesh_t p_mesh_0, tensor_ref_3 p_pos_1)
{
    vec3 v_2 = vcons3(
        0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e0,
        0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e0,
        0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e0);
    diderot::dynseq< int32_t > t_4 = mesh_geom_mesh_t(p_mesh_0.index, &p_mesh_0, p_pos_1._data);
    vec3 v_5 = v_2;
    for (auto it_0 = t_4.cbegin(); it_0 != t_4.cend(); ++it_0) {
        auto i_cellItter_3 = *it_0;
        vec3 v_6;
        v_6 = v_5;
        for (int32_t i_newtonItter_7 = 0; i_newtonItter_7 <= 16; ++i_newtonItter_7) {
            int32_t l_mulRes_8 = i_cellItter_3 * 20;
            int32_t t_9 = p_mesh_0.indexMap[l_mulRes_8];
            int32_t l_mulRes_10 = 3 * t_9;
            double l_dof_load_11 = p_mesh_0.coordMap[l_mulRes_10];
            double l_dof_load_12 = p_mesh_0.coordMap[1 + l_mulRes_10];
            double l_dof_load_13 = p_mesh_0.coordMap[2 + l_mulRes_10];
            int32_t t_14 = p_mesh_0.indexMap[l_mulRes_8 + 1];
            int32_t l_mulRes_15 = 3 * t_14;
            double l_dof_load_16 = p_mesh_0.coordMap[l_mulRes_15];
            double l_dof_load_17 = p_mesh_0.coordMap[1 + l_mulRes_15];
            double l_dof_load_18 = p_mesh_0.coordMap[2 + l_mulRes_15];
            int32_t t_19 = p_mesh_0.indexMap[l_mulRes_8 + 2];
            int32_t l_mulRes_20 = 3 * t_19;
            double l_dof_load_21 = p_mesh_0.coordMap[l_mulRes_20];
            double l_dof_load_22 = p_mesh_0.coordMap[1 + l_mulRes_20];
            double l_dof_load_23 = p_mesh_0.coordMap[2 + l_mulRes_20];
            int32_t t_24 = p_mesh_0.indexMap[l_mulRes_8 + 3];
            int32_t l_mulRes_25 = 3 * t_24;
            double l_dof_load_26 = p_mesh_0.coordMap[l_mulRes_25];
            double l_dof_load_27 = p_mesh_0.coordMap[1 + l_mulRes_25];
            double l_dof_load_28 = p_mesh_0.coordMap[2 + l_mulRes_25];
            int32_t t_29 = p_mesh_0.indexMap[l_mulRes_8 + 4];
            int32_t l_mulRes_30 = 3 * t_29;
            double l_dof_load_31 = p_mesh_0.coordMap[l_mulRes_30];
            double l_dof_load_32 = p_mesh_0.coordMap[1 + l_mulRes_30];
            double l_dof_load_33 = p_mesh_0.coordMap[2 + l_mulRes_30];
            int32_t t_34 = p_mesh_0.indexMap[l_mulRes_8 + 5];
            int32_t l_mulRes_35 = 3 * t_34;
            double l_dof_load_36 = p_mesh_0.coordMap[l_mulRes_35];
            double l_dof_load_37 = p_mesh_0.coordMap[1 + l_mulRes_35];
            double l_dof_load_38 = p_mesh_0.coordMap[2 + l_mulRes_35];
            int32_t t_39 = p_mesh_0.indexMap[l_mulRes_8 + 6];
            int32_t l_mulRes_40 = 3 * t_39;
            double l_dof_load_41 = p_mesh_0.coordMap[l_mulRes_40];
            double l_dof_load_42 = p_mesh_0.coordMap[1 + l_mulRes_40];
            double l_dof_load_43 = p_mesh_0.coordMap[2 + l_mulRes_40];
            int32_t t_44 = p_mesh_0.indexMap[l_mulRes_8 + 7];
            int32_t l_mulRes_45 = 3 * t_44;
            double l_dof_load_46 = p_mesh_0.coordMap[l_mulRes_45];
            double l_dof_load_47 = p_mesh_0.coordMap[1 + l_mulRes_45];
            double l_dof_load_48 = p_mesh_0.coordMap[2 + l_mulRes_45];
            int32_t t_49 = p_mesh_0.indexMap[l_mulRes_8 + 8];
            int32_t l_mulRes_50 = 3 * t_49;
            double l_dof_load_51 = p_mesh_0.coordMap[l_mulRes_50];
            double l_dof_load_52 = p_mesh_0.coordMap[1 + l_mulRes_50];
            double l_dof_load_53 = p_mesh_0.coordMap[2 + l_mulRes_50];
            int32_t t_54 = p_mesh_0.indexMap[l_mulRes_8 + 9];
            int32_t l_mulRes_55 = 3 * t_54;
            double l_dof_load_56 = p_mesh_0.coordMap[l_mulRes_55];
            double l_dof_load_57 = p_mesh_0.coordMap[1 + l_mulRes_55];
            double l_dof_load_58 = p_mesh_0.coordMap[2 + l_mulRes_55];
            int32_t t_59 = p_mesh_0.indexMap[l_mulRes_8 + 10];
            int32_t l_mulRes_60 = 3 * t_59;
            double l_dof_load_61 = p_mesh_0.coordMap[l_mulRes_60];
            double l_dof_load_62 = p_mesh_0.coordMap[1 + l_mulRes_60];
            double l_dof_load_63 = p_mesh_0.coordMap[2 + l_mulRes_60];
            int32_t t_64 = p_mesh_0.indexMap[l_mulRes_8 + 11];
            int32_t l_mulRes_65 = 3 * t_64;
            double l_dof_load_66 = p_mesh_0.coordMap[l_mulRes_65];
            double l_dof_load_67 = p_mesh_0.coordMap[1 + l_mulRes_65];
            double l_dof_load_68 = p_mesh_0.coordMap[2 + l_mulRes_65];
            int32_t t_69 = p_mesh_0.indexMap[l_mulRes_8 + 12];
            int32_t l_mulRes_70 = 3 * t_69;
            double l_dof_load_71 = p_mesh_0.coordMap[l_mulRes_70];
            double l_dof_load_72 = p_mesh_0.coordMap[1 + l_mulRes_70];
            double l_dof_load_73 = p_mesh_0.coordMap[2 + l_mulRes_70];
            int32_t t_74 = p_mesh_0.indexMap[l_mulRes_8 + 13];
            int32_t l_mulRes_75 = 3 * t_74;
            double l_dof_load_76 = p_mesh_0.coordMap[l_mulRes_75];
            double l_dof_load_77 = p_mesh_0.coordMap[1 + l_mulRes_75];
            double l_dof_load_78 = p_mesh_0.coordMap[2 + l_mulRes_75];
            int32_t t_79 = p_mesh_0.indexMap[l_mulRes_8 + 14];
            int32_t l_mulRes_80 = 3 * t_79;
            double l_dof_load_81 = p_mesh_0.coordMap[l_mulRes_80];
            double l_dof_load_82 = p_mesh_0.coordMap[1 + l_mulRes_80];
            double l_dof_load_83 = p_mesh_0.coordMap[2 + l_mulRes_80];
            int32_t t_84 = p_mesh_0.indexMap[l_mulRes_8 + 15];
            int32_t l_mulRes_85 = 3 * t_84;
            double l_dof_load_86 = p_mesh_0.coordMap[l_mulRes_85];
            double l_dof_load_87 = p_mesh_0.coordMap[1 + l_mulRes_85];
            double l_dof_load_88 = p_mesh_0.coordMap[2 + l_mulRes_85];
            int32_t t_89 = p_mesh_0.indexMap[l_mulRes_8 + 16];
            int32_t l_mulRes_90 = 3 * t_89;
            double l_dof_load_91 = p_mesh_0.coordMap[l_mulRes_90];
            double l_dof_load_92 = p_mesh_0.coordMap[1 + l_mulRes_90];
            double l_dof_load_93 = p_mesh_0.coordMap[2 + l_mulRes_90];
            int32_t t_94 = p_mesh_0.indexMap[l_mulRes_8 + 17];
            int32_t l_mulRes_95 = 3 * t_94;
            double l_dof_load_96 = p_mesh_0.coordMap[l_mulRes_95];
            double l_dof_load_97 = p_mesh_0.coordMap[1 + l_mulRes_95];
            double l_dof_load_98 = p_mesh_0.coordMap[2 + l_mulRes_95];
            int32_t t_99 = p_mesh_0.indexMap[l_mulRes_8 + 18];
            int32_t l_mulRes_100 = 3 * t_99;
            double l_dof_load_101 = p_mesh_0.coordMap[l_mulRes_100];
            double l_dof_load_102 = p_mesh_0.coordMap[1 + l_mulRes_100];
            double l_dof_load_103 = p_mesh_0.coordMap[2 + l_mulRes_100];
            int32_t t_104 = p_mesh_0.indexMap[l_mulRes_8 + 19];
            int32_t l_mulRes_105 = 3 * t_104;
            double l_dof_load_106 = p_mesh_0.coordMap[l_mulRes_105];
            double l_dof_load_107 = p_mesh_0.coordMap[1 + l_mulRes_105];
            double l_dof_load_108 = p_mesh_0.coordMap[2 + l_mulRes_105];
            double l_varAcc_109 = v_5[0];
            double l_varAcc_110 = v_5[1];
            double l_varAcc_111 = v_5[2];
            double l_prod_112 = 0.1e1 * 0.1e1;
            double l_prod_113 = l_varAcc_109 * l_varAcc_109 * l_prod_112;
            double l_prod_114 = l_varAcc_110 * 0.1e1;
            double l_prod_115 = l_varAcc_109 * l_prod_114;
            double l_prod_116 = 0.1e1 * l_varAcc_111;
            double l_prod_117 = l_varAcc_109 * l_prod_116;
            double l_prod_118 = l_varAcc_109 * l_prod_112;
            double l_prod_119 = 0.1e1 * (l_varAcc_110 * l_varAcc_110 * 0.1e1);
            double l_prod_120 = 0.1e1 * (l_varAcc_110 * l_varAcc_111);
            double l_prod_121 = 0.1e1 * l_prod_114;
            double l_prod_122 = 0.1e1 * (0.1e1 * (l_varAcc_111 * l_varAcc_111));
            double l_prod_123 = 0.1e1 * l_prod_116;
            double l_prod_124 = 0.1e1 * l_prod_112;
            double l_mult_125 = -0.135e2 * l_prod_122;
            double l_mult_126 = -0.27e2 * l_prod_120;
            double l_mult_127 = -0.135e2 * l_prod_119;
            double l_mult_128 = -0.27e2 * l_prod_117;
            double l_mult_129 = -0.27e2 * l_prod_115;
            double l_mult_130 = -0.135e2 * l_prod_113;
            double l_sum_131 = -0.55e1 * l_prod_124 + (0.18e2 * l_prod_123 + (l_mult_125 + (0.18e2 * l_prod_121 + (l_mult_126 + (l_mult_127 + (0.18e2 * l_prod_118 + (l_mult_128 + (l_mult_129 + l_mult_130))))))));
            double l_mult_132 = 0.1e1 * l_prod_124;
            double l_mult_133 = 0.135e2 * l_prod_113;
            double l_sum_134 = l_mult_132 + (-0.9e1 * l_prod_118 + l_mult_133);
            double l_mult_135 = -0.45e1 * l_prod_123;
            double l_mult_136 = 0.27e2 * l_prod_117;
            double l_sum_137 = l_mult_135 + l_mult_136;
            double l_mult_138 = 0.135e2 * l_prod_122;
            double l_sum_139 = l_mult_135 + l_mult_138;
            double l_mult_140 = -0.45e1 * l_prod_121;
            double l_mult_141 = 0.27e2 * l_prod_115;
            double l_sum_142 = l_mult_140 + l_mult_141;
            double l_mult_143 = 0.135e2 * l_prod_119;
            double l_sum_144 = l_mult_140 + l_mult_143;
            double l_mult_145 = -0.225e2 * l_prod_123;
            double l_mult_146 = 0.27e2 * l_prod_120;
            double l_sum_147 = l_mult_145 + (0.27e2 * l_prod_122 + (l_mult_146 + l_mult_136));
            double l_mult_148 = 0.45e1 * l_prod_123;
            double l_sum_149 = l_mult_148 + l_mult_125;
            double l_mult_150 = -0.225e2 * l_prod_121;
            double l_sum_151 = l_mult_150 + (l_mult_146 + (0.27e2 * l_prod_119 + l_mult_141));
            double l_mult_152 = 0.45e1 * l_prod_121;
            double l_sum_153 = l_mult_152 + l_mult_127;
            double l_mult_154 = 0.9e1 * l_prod_124;
            double l_mult_155 = 0.54e2 * l_prod_117;
            double l_mult_156 = 0.54e2 * l_prod_115;
            double l_sum_157 = l_mult_154 + (l_mult_145 + (l_mult_138 + (l_mult_150 + (l_mult_146 + (l_mult_143 + (-0.45e2 * l_prod_118 + (l_mult_155 + (l_mult_156 + 0.405e2 * l_prod_113))))))));
            double l_mult_158 = -0.45e1 * l_prod_124;
            double l_sum_159 = l_mult_158 + (l_mult_148 + (l_mult_152 + (0.36e2 * l_prod_118 + (l_mult_128 + (l_mult_129 + -0.405e2 * l_prod_113)))));
            double l_mult_160 = 0.27e2 * l_prod_123;
            double l_mult_161 = -0.27e2 * l_prod_122;
            double l_mult_162 = -0.54e2 * l_prod_117;
            double l_sum_163 = l_mult_160 + (l_mult_161 + (l_mult_126 + l_mult_162));
            double l_mult_164 = 0.27e2 * l_prod_121;
            double l_mult_165 = -0.27e2 * l_prod_119;
            double l_mult_166 = -0.54e2 * l_prod_115;
            double l_sum_167 = l_mult_164 + (l_mult_126 + (l_mult_165 + l_mult_166));
            double l_sum_168 = l_mult_132 + (-0.9e1 * l_prod_121 + l_mult_143);
            double l_sum_169 = l_mult_135 + l_mult_146;
            double l_mult_170 = -0.45e1 * l_prod_118;
            double l_sum_171 = l_mult_170 + l_mult_133;
            double l_sum_172 = l_mult_170 + l_mult_141;
            double l_mult_173 = 0.54e2 * l_prod_120;
            double l_mult_174 = -0.225e2 * l_prod_118;
            double l_sum_175 = l_mult_154 + (l_mult_145 + (l_mult_138 + (-0.45e2 * l_prod_121 + (l_mult_173 + (0.405e2 * l_prod_119 + (l_mult_174 + (l_mult_136 + (l_mult_156 + l_mult_133))))))));
            double l_mult_176 = 0.45e1 * l_prod_118;
            double l_sum_177 = l_mult_158 + (l_mult_148 + (0.36e2 * l_prod_121 + (l_mult_126 + (-0.405e2 * l_prod_119 + (l_mult_176 + l_mult_129)))));
            double l_sum_178 = l_mult_174 + (l_mult_136 + (l_mult_141 + 0.27e2 * l_prod_113));
            double l_sum_179 = l_mult_176 + l_mult_130;
            double l_mult_180 = -0.54e2 * l_prod_120;
            double l_sum_181 = l_mult_160 + (l_mult_161 + (l_mult_180 + l_mult_128));
            double l_mult_182 = 0.27e2 * l_prod_118;
            double l_mult_183 = -0.27e2 * l_prod_113;
            double l_sum_184 = l_mult_182 + (l_mult_128 + (l_mult_166 + l_mult_183));
            double l_sum_185 = l_mult_132 + (-0.9e1 * l_prod_123 + l_mult_138);
            double l_sum_186 = l_mult_140 + l_mult_146;
            double l_sum_187 = l_mult_170 + l_mult_136;
            double l_sum_188 = l_mult_154 + (-0.45e2 * l_prod_123 + (0.405e2 * l_prod_122 + (l_mult_150 + (l_mult_173 + (l_mult_143 + (l_mult_174 + (l_mult_155 + (l_mult_141 + l_mult_133))))))));
            double l_sum_189 = l_mult_158 + (0.36e2 * l_prod_123 + (-0.405e2 * l_prod_122 + (l_mult_152 + (l_mult_126 + (l_mult_176 + l_mult_128)))));
            double l_sum_190 = l_mult_164 + (l_mult_180 + (l_mult_165 + l_mult_129));
            double l_sum_191 = l_mult_182 + (l_mult_162 + (l_mult_129 + l_mult_183));
            double l_r_192 = l_dof_load_11 * l_sum_131;
            double l_r_193 = l_dof_load_21 * 0.e0;
            double l_r_194 = l_dof_load_26 * 0.e0;
            double l_r_195 = l_dof_load_61 * l_sum_147;
            double l_r_196 = l_dof_load_66 * l_sum_149;
            double l_r_197 = l_dof_load_71 * l_sum_151;
            double l_r_198 = l_dof_load_76 * l_sum_153;
            double l_r_199 = l_r_192 + l_dof_load_16 * l_sum_134 + l_r_193 + l_r_194 + l_dof_load_31 * 0.e0 + l_dof_load_36 * 0.e0 + l_dof_load_41 * l_sum_137 + l_dof_load_46 * l_sum_139 + l_dof_load_51 * l_sum_142 + l_dof_load_56 * l_sum_144 + l_r_195 + l_r_196 + l_r_197 + l_r_198 + l_dof_load_81 * l_sum_157 + l_dof_load_86 * l_sum_159 + l_dof_load_91 * l_mult_146 + l_dof_load_96 * l_mult_126 + l_dof_load_101 * l_sum_163 + l_dof_load_106 * l_sum_167;
            double l_r_200 = l_dof_load_81 * l_sum_178;
            double l_r_201 = l_dof_load_86 * l_sum_179;
            double l_r_202 = l_r_192 + l_dof_load_16 * 0.e0;
            double l_r_203 = l_r_202 + l_dof_load_21 * l_sum_168 + l_r_194 + l_dof_load_31 * l_sum_169 + l_dof_load_36 * l_sum_139 + l_dof_load_41 * 0.e0 + l_dof_load_46 * 0.e0 + l_dof_load_51 * l_sum_171 + l_dof_load_56 * l_sum_172 + l_r_195 + l_r_196 + l_dof_load_71 * l_sum_175 + l_dof_load_76 * l_sum_177 + l_r_200 + l_r_201 + l_dof_load_91 * l_mult_136 + l_dof_load_96 * l_sum_181 + l_dof_load_101 * l_mult_128 + l_dof_load_106 * l_sum_184;
            double l_r_204 = l_r_202 + l_r_193 + l_dof_load_26 * l_sum_185 + l_dof_load_31 * l_sum_144 + l_dof_load_36 * l_sum_186 + l_dof_load_41 * l_sum_171 + l_dof_load_46 * l_sum_187 + l_dof_load_51 * 0.e0 + l_dof_load_56 * 0.e0 + l_dof_load_61 * l_sum_188 + l_dof_load_66 * l_sum_189 + l_r_197 + l_r_198 + l_r_200 + l_r_201 + l_dof_load_91 * l_mult_141 + l_dof_load_96 * l_sum_190 + l_dof_load_101 * l_sum_191 + l_dof_load_106 * l_mult_129;
            double l_r_205 = l_dof_load_12 * l_sum_131;
            double l_r_206 = l_dof_load_22 * 0.e0;
            double l_r_207 = l_dof_load_27 * 0.e0;
            double l_r_208 = l_dof_load_62 * l_sum_147;
            double l_r_209 = l_dof_load_67 * l_sum_149;
            double l_r_210 = l_dof_load_72 * l_sum_151;
            double l_r_211 = l_dof_load_77 * l_sum_153;
            double l_r_212 = l_r_205 + l_dof_load_17 * l_sum_134 + l_r_206 + l_r_207 + l_dof_load_32 * 0.e0 + l_dof_load_37 * 0.e0 + l_dof_load_42 * l_sum_137 + l_dof_load_47 * l_sum_139 + l_dof_load_52 * l_sum_142 + l_dof_load_57 * l_sum_144 + l_r_208 + l_r_209 + l_r_210 + l_r_211 + l_dof_load_82 * l_sum_157 + l_dof_load_87 * l_sum_159 + l_dof_load_92 * l_mult_146 + l_dof_load_97 * l_mult_126 + l_dof_load_102 * l_sum_163 + l_dof_load_107 * l_sum_167;
            double l_r_213 = l_dof_load_82 * l_sum_178;
            double l_r_214 = l_dof_load_87 * l_sum_179;
            double l_r_215 = l_r_205 + l_dof_load_17 * 0.e0;
            double l_r_216 = l_r_215 + l_dof_load_22 * l_sum_168 + l_r_207 + l_dof_load_32 * l_sum_169 + l_dof_load_37 * l_sum_139 + l_dof_load_42 * 0.e0 + l_dof_load_47 * 0.e0 + l_dof_load_52 * l_sum_171 + l_dof_load_57 * l_sum_172 + l_r_208 + l_r_209 + l_dof_load_72 * l_sum_175 + l_dof_load_77 * l_sum_177 + l_r_213 + l_r_214 + l_dof_load_92 * l_mult_136 + l_dof_load_97 * l_sum_181 + l_dof_load_102 * l_mult_128 + l_dof_load_107 * l_sum_184;
            double l_r_217 = l_r_215 + l_r_206 + l_dof_load_27 * l_sum_185 + l_dof_load_32 * l_sum_144 + l_dof_load_37 * l_sum_186 + l_dof_load_42 * l_sum_171 + l_dof_load_47 * l_sum_187 + l_dof_load_52 * 0.e0 + l_dof_load_57 * 0.e0 + l_dof_load_62 * l_sum_188 + l_dof_load_67 * l_sum_189 + l_r_210 + l_r_211 + l_r_213 + l_r_214 + l_dof_load_92 * l_mult_141 + l_dof_load_97 * l_sum_190 + l_dof_load_102 * l_sum_191 + l_dof_load_107 * l_mult_129;
            double l_r_218 = l_dof_load_13 * l_sum_131;
            double l_r_219 = l_dof_load_23 * 0.e0;
            double l_r_220 = l_dof_load_28 * 0.e0;
            double l_r_221 = l_dof_load_63 * l_sum_147;
            double l_r_222 = l_dof_load_68 * l_sum_149;
            double l_r_223 = l_dof_load_73 * l_sum_151;
            double l_r_224 = l_dof_load_78 * l_sum_153;
            double l_r_225 = l_r_218 + l_dof_load_18 * l_sum_134 + l_r_219 + l_r_220 + l_dof_load_33 * 0.e0 + l_dof_load_38 * 0.e0 + l_dof_load_43 * l_sum_137 + l_dof_load_48 * l_sum_139 + l_dof_load_53 * l_sum_142 + l_dof_load_58 * l_sum_144 + l_r_221 + l_r_222 + l_r_223 + l_r_224 + l_dof_load_83 * l_sum_157 + l_dof_load_88 * l_sum_159 + l_dof_load_93 * l_mult_146 + l_dof_load_98 * l_mult_126 + l_dof_load_103 * l_sum_163 + l_dof_load_108 * l_sum_167;
            double l_r_226 = l_dof_load_83 * l_sum_178;
            double l_r_227 = l_dof_load_88 * l_sum_179;
            double l_r_228 = l_r_218 + l_dof_load_18 * 0.e0;
            double l_r_229 = l_r_228 + l_dof_load_23 * l_sum_168 + l_r_220 + l_dof_load_33 * l_sum_169 + l_dof_load_38 * l_sum_139 + l_dof_load_43 * 0.e0 + l_dof_load_48 * 0.e0 + l_dof_load_53 * l_sum_171 + l_dof_load_58 * l_sum_172 + l_r_221 + l_r_222 + l_dof_load_73 * l_sum_175 + l_dof_load_78 * l_sum_177 + l_r_226 + l_r_227 + l_dof_load_93 * l_mult_136 + l_dof_load_98 * l_sum_181 + l_dof_load_103 * l_mult_128 + l_dof_load_108 * l_sum_184;
            double l_r_230 = l_r_228 + l_r_219 + l_dof_load_28 * l_sum_185 + l_dof_load_33 * l_sum_144 + l_dof_load_38 * l_sum_186 + l_dof_load_43 * l_sum_171 + l_dof_load_48 * l_sum_187 + l_dof_load_53 * 0.e0 + l_dof_load_58 * 0.e0 + l_dof_load_63 * l_sum_188 + l_dof_load_68 * l_sum_189 + l_r_223 + l_r_224 + l_r_226 + l_r_227 + l_dof_load_93 * l_mult_141 + l_dof_load_98 * l_sum_190 + l_dof_load_103 * l_sum_191 + l_dof_load_108 * l_mult_129;
            double l_r_231 = 0.e0 * l_r_199;
            double l_r_232 = 0.e0 * l_r_212;
            double l_r_233 = 0.e0 * l_r_225;
            double l_r_234 = l_r_231 + l_r_232;
            double l_r_235 = l_r_234 + l_r_233;
            double l_r_236 = 0.e0 * l_r_203;
            double l_r_237 = 0.e0 * l_r_216;
            double l_r_238 = 0.e0 * l_r_229;
            double l_r_239 = l_r_236 + l_r_237;
            double l_r_240 = l_r_239 + l_r_238;
            double l_r_241 = 0.e0 * l_r_204;
            double l_r_242 = 0.e0 * l_r_217;
            double l_r_243 = 0.e0 * l_r_230;
            double l_r_244 = l_r_241 + l_r_242;
            double l_r_245 = l_r_244 + l_r_243;
            double l_r_246 = l_r_234 + -0.1e1 * l_r_225;
            double l_r_247 = l_r_239 + -0.1e1 * l_r_229;
            double l_r_248 = l_r_244 + -0.1e1 * l_r_230;
            double l_r_249 = l_r_231 + 0.1e1 * l_r_212 + l_r_233;
            double l_r_250 = l_r_236 + 0.1e1 * l_r_216 + l_r_238;
            double l_r_251 = l_r_241 + 0.1e1 * l_r_217 + l_r_243;
            double l_r_252 = l_r_234 + 0.1e1 * l_r_225;
            double l_r_253 = l_r_239 + 0.1e1 * l_r_229;
            double l_r_254 = l_r_244 + 0.1e1 * l_r_230;
            double l_r_255 = -0.1e1 * l_r_199 + l_r_232 + l_r_233;
            double l_r_256 = -0.1e1 * l_r_203 + l_r_237 + l_r_238;
            double l_r_257 = -0.1e1 * l_r_204 + l_r_242 + l_r_243;
            double l_r_258 = l_r_231 + -0.1e1 * l_r_212 + l_r_233;
            double l_r_259 = l_r_236 + -0.1e1 * l_r_216 + l_r_238;
            double l_r_260 = l_r_241 + -0.1e1 * l_r_217 + l_r_243;
            double l_r_261 = 0.1e1 * l_r_199 + l_r_232 + l_r_233;
            double l_r_262 = 0.1e1 * l_r_203 + l_r_237 + l_r_238;
            double l_r_263 = 0.1e1 * l_r_204 + l_r_242 + l_r_243;
            double l_r_264 = l_r_199 * l_r_240 + l_r_212 * l_r_253 + l_r_225 * l_r_259;
            double l_r_265 = l_r_199 * l_r_245 + l_r_212 * l_r_254 + l_r_225 * l_r_260;
            double l_r_266 = l_r_199 * l_r_247 + l_r_212 * l_r_240 + l_r_225 * l_r_262;
            double l_r_267 = l_r_199 * l_r_248 + l_r_212 * l_r_245 + l_r_225 * l_r_263;
            double l_r_268 = l_r_199 * l_r_250 + l_r_212 * l_r_256 + l_r_225 * l_r_240;
            double l_r_269 = l_r_199 * l_r_251 + l_r_212 * l_r_257 + l_r_225 * l_r_245;
            double l_r_270 = l_r_203 * l_r_235 + l_r_216 * l_r_252 + l_r_229 * l_r_258;
            double l_r_271 = l_r_203 * l_r_245 + l_r_216 * l_r_254 + l_r_229 * l_r_260;
            double l_r_272 = l_r_203 * l_r_246 + l_r_216 * l_r_235 + l_r_229 * l_r_261;
            double l_r_273 = l_r_203 * l_r_248 + l_r_216 * l_r_245 + l_r_229 * l_r_263;
            double l_r_274 = l_r_203 * l_r_249 + l_r_216 * l_r_255 + l_r_229 * l_r_235;
            double l_r_275 = l_r_203 * l_r_251 + l_r_216 * l_r_257 + l_r_229 * l_r_245;
            double l_r_276 = l_r_204 * l_r_235 + l_r_217 * l_r_252 + l_r_230 * l_r_258;
            double l_r_277 = l_r_204 * l_r_240 + l_r_217 * l_r_253 + l_r_230 * l_r_259;
            double l_r_278 = l_r_204 * l_r_246 + l_r_217 * l_r_235 + l_r_230 * l_r_261;
            double l_r_279 = l_r_204 * l_r_247 + l_r_217 * l_r_240 + l_r_230 * l_r_262;
            double l_r_280 = l_r_204 * l_r_249 + l_r_217 * l_r_255 + l_r_230 * l_r_235;
            double l_r_281 = l_r_204 * l_r_250 + l_r_217 * l_r_256 + l_r_230 * l_r_240;
            vec3 v_282 = vcons3(l_r_203, l_r_216, l_r_229);
            double l_r_283 = 0.e0 * (l_r_199 * l_r_235 + l_r_212 * l_r_252 + l_r_225 * l_r_258);
            double l_r_284 = 0.e0 * l_r_265;
            double l_r_285 = 0.e0 * l_r_270;
            double l_r_286 = 0.e0 * (l_r_203 * l_r_240 + l_r_216 * l_r_253 + l_r_229 * l_r_259);
            double l_r_287 = 0.e0 * l_r_276;
            double l_r_288 = 0.e0 * (l_r_204 * l_r_245 + l_r_217 * l_r_254 + l_r_230 * l_r_260);
            double l_r_289 = l_r_283 + 0.e0 * l_r_264;
            double l_r_290 = 0.e0 * (l_r_199 * l_r_246 + l_r_212 * l_r_235 + l_r_225 * l_r_261);
            double l_r_291 = 0.e0 * l_r_267;
            double l_r_292 = 0.e0 * l_r_272;
            double l_r_293 = 0.e0 * (l_r_203 * l_r_247 + l_r_216 * l_r_240 + l_r_229 * l_r_262);
            double l_r_294 = 0.e0 * l_r_278;
            double l_r_295 = 0.e0 * (l_r_204 * l_r_248 + l_r_217 * l_r_245 + l_r_230 * l_r_263);
            double l_r_296 = l_r_290 + 0.e0 * l_r_266;
            double l_r_297 = 0.e0 * (l_r_199 * l_r_249 + l_r_212 * l_r_255 + l_r_225 * l_r_235);
            double l_r_298 = 0.e0 * l_r_269;
            double l_r_299 = 0.e0 * l_r_274;
            double l_r_300 = 0.e0 * (l_r_203 * l_r_250 + l_r_216 * l_r_256 + l_r_229 * l_r_240);
            double l_r_301 = 0.e0 * l_r_280;
            double l_r_302 = 0.e0 * (l_r_204 * l_r_251 + l_r_217 * l_r_257 + l_r_230 * l_r_245);
            double l_r_303 = l_r_297 + 0.e0 * l_r_268;
            double l_r_304 = 0.e0 * l_r_271;
            double l_r_305 = 0.e0 * l_r_277;
            double l_r_306 = 0.e0 * l_r_273;
            double l_r_307 = 0.e0 * l_r_279;
            double l_r_308 = 0.e0 * l_r_275;
            double l_r_309 = 0.e0 * l_r_281;
            double l_op1_e3_l_21_310 = 0.2e1 * vdot3(vcons3(l_r_199, l_r_212, l_r_225),
                vcons3(vdot3(v_282, vcons3(l_r_245, l_r_254, l_r_260)),
                    vdot3(v_282, vcons3(l_r_248, l_r_245, l_r_263)), vdot3(v_282, vcons3(l_r_251, l_r_257, l_r_245))));
            double l_varAcc_311 = v_6[0];
            double l_varAcc_312 = v_6[1];
            double l_varAcc_313 = v_6[2];
            double l_prod2_314 = l_varAcc_311 * l_varAcc_311;
            double l_prod_315 = l_prod2_314 * l_varAcc_311 * l_prod_112;
            double l_prod_316 = l_varAcc_312 * 0.1e1;
            double l_prod_317 = l_prod2_314 * l_prod_316;
            double l_prod_318 = 0.1e1 * l_varAcc_313;
            double l_prod_319 = l_prod2_314 * l_prod_318;
            double l_prod_320 = l_prod2_314 * l_prod_112;
            double l_prod2_321 = l_varAcc_312 * l_varAcc_312;
            double l_prod_322 = l_prod2_321 * 0.1e1;
            double l_prod_323 = l_varAcc_311 * l_prod_322;
            double l_prod_324 = l_varAcc_312 * l_varAcc_313;
            double l_prod_325 = l_varAcc_311 * l_prod_324;
            double l_prod_326 = l_varAcc_311 * l_prod_316;
            double l_prod2_327 = l_varAcc_313 * l_varAcc_313;
            double l_prod_328 = 0.1e1 * l_prod2_327;
            double l_prod_329 = l_varAcc_311 * l_prod_328;
            double l_prod_330 = l_varAcc_311 * l_prod_318;
            double l_prod_331 = l_varAcc_311 * l_prod_112;
            double l_prod_332 = 0.1e1 * (l_prod2_321 * l_varAcc_312 * 0.1e1);
            double l_prod_333 = 0.1e1 * (l_prod2_321 * l_varAcc_313);
            double l_prod_334 = 0.1e1 * l_prod_322;
            double l_prod_335 = 0.1e1 * (l_varAcc_312 * l_prod2_327);
            double l_prod_336 = 0.1e1 * l_prod_324;
            double l_prod_337 = 0.1e1 * l_prod_316;
            double l_prod_338 = 0.1e1 * (0.1e1 * (l_prod2_327 * l_varAcc_313));
            double l_prod_339 = 0.1e1 * l_prod_328;
            double l_prod_340 = 0.1e1 * l_prod_318;
            double l_mult_341 = -0.135e2 * l_prod_335;
            double l_mult_342 = -0.135e2 * l_prod_333;
            double l_mult_343 = -0.135e2 * l_prod_329;
            double l_mult_344 = -0.27e2 * l_prod_325;
            double l_mult_345 = -0.135e2 * l_prod_323;
            double l_mult_346 = -0.135e2 * l_prod_319;
            double l_mult_347 = -0.135e2 * l_prod_317;
            double l_sum_348 = l_mult_132 + (-0.55e1 * l_prod_340 + (0.9e1 * l_prod_339 + (-0.45e1 * l_prod_338 + (-0.55e1 * l_prod_337 + (0.18e2 * l_prod_336 + (l_mult_341 + (0.9e1 * l_prod_334 + (l_mult_342 + (-0.45e1 * l_prod_332 + (-0.55e1 * l_prod_331 + (0.18e2 * l_prod_330 + (l_mult_343 + (0.18e2 * l_prod_326 + (l_mult_344 + (l_mult_345 + (0.9e1 * l_prod_320 + (l_mult_346 + (l_mult_347 + -0.45e1 * l_prod_315))))))))))))))))));
            double l_sum_349 = 0.1e1 * l_prod_331 + (-0.45e1 * l_prod_320 + 0.45e1 * l_prod_315);
            double l_sum_350 = 0.1e1 * l_prod_337 + (-0.45e1 * l_prod_334 + 0.45e1 * l_prod_332);
            double l_sum_351 = 0.1e1 * l_prod_340 + (-0.45e1 * l_prod_339 + 0.45e1 * l_prod_338);
            double l_mult_352 = -0.45e1 * l_prod_336;
            double l_mult_353 = 0.135e2 * l_prod_333;
            double l_sum_354 = l_mult_352 + l_mult_353;
            double l_mult_355 = 0.135e2 * l_prod_335;
            double l_sum_356 = l_mult_352 + l_mult_355;
            double l_mult_357 = -0.45e1 * l_prod_330;
            double l_mult_358 = 0.135e2 * l_prod_319;
            double l_sum_359 = l_mult_357 + l_mult_358;
            double l_mult_360 = 0.135e2 * l_prod_329;
            double l_sum_361 = l_mult_357 + l_mult_360;
            double l_mult_362 = -0.45e1 * l_prod_326;
            double l_mult_363 = 0.135e2 * l_prod_317;
            double l_sum_364 = l_mult_362 + l_mult_363;
            double l_mult_365 = 0.135e2 * l_prod_323;
            double l_sum_366 = l_mult_362 + l_mult_365;
            double l_mult_367 = -0.225e2 * l_prod_336;
            double l_mult_368 = -0.225e2 * l_prod_330;
            double l_mult_369 = 0.27e2 * l_prod_325;
            double l_sum_370 = 0.9e1 * l_prod_340 + (-0.225e2 * l_prod_339 + (0.135e2 * l_prod_338 + (l_mult_367 + (0.27e2 * l_prod_335 + (l_mult_353 + (l_mult_368 + (0.27e2 * l_prod_329 + (l_mult_369 + l_mult_358))))))));
            double l_mult_371 = 0.45e1 * l_prod_336;
            double l_mult_372 = 0.45e1 * l_prod_330;
            double l_sum_373 = -0.45e1 * l_prod_340 + (0.18e2 * l_prod_339 + (-0.135e2 * l_prod_338 + (l_mult_371 + (l_mult_341 + (l_mult_372 + l_mult_343)))));
            double l_mult_374 = -0.225e2 * l_prod_326;
            double l_sum_375 = 0.9e1 * l_prod_337 + (l_mult_367 + (l_mult_355 + (-0.225e2 * l_prod_334 + (0.27e2 * l_prod_333 + (0.135e2 * l_prod_332 + (l_mult_374 + (l_mult_369 + (0.27e2 * l_prod_323 + l_mult_363))))))));
            double l_mult_376 = 0.45e1 * l_prod_326;
            double l_sum_377 = -0.45e1 * l_prod_337 + (l_mult_371 + (0.18e2 * l_prod_334 + (l_mult_342 + (-0.135e2 * l_prod_332 + (l_mult_376 + l_mult_345)))));
            double l_sum_378 = 0.9e1 * l_prod_331 + (l_mult_368 + (l_mult_360 + (l_mult_374 + (l_mult_369 + (l_mult_365 + (-0.225e2 * l_prod_320 + (0.27e2 * l_prod_319 + (0.27e2 * l_prod_317 + 0.135e2 * l_prod_315))))))));
            double l_sum_379 = -0.45e1 * l_prod_331 + (l_mult_372 + (l_mult_376 + (0.18e2 * l_prod_320 + (l_mult_346 + (l_mult_347 + -0.135e2 * l_prod_315)))));
            double l_sum_380 = 0.27e2 * l_prod_336 + (-0.27e2 * l_prod_335 + (-0.27e2 * l_prod_333 + l_mult_344));
            double l_sum_381 = 0.27e2 * l_prod_330 + (-0.27e2 * l_prod_329 + (l_mult_344 + -0.27e2 * l_prod_319));
            double l_sum_382 = 0.27e2 * l_prod_326 + (l_mult_344 + (-0.27e2 * l_prod_323 + -0.27e2 * l_prod_317));
            vec3 v_383 = vcons3(
                l_dof_load_11 * l_sum_348 + l_dof_load_16 * l_sum_349 + l_dof_load_21 * l_sum_350 + l_dof_load_26 * l_sum_351 + l_dof_load_31 * l_sum_354 + l_dof_load_36 * l_sum_356 + l_dof_load_41 * l_sum_359 + l_dof_load_46 * l_sum_361 + l_dof_load_51 * l_sum_364 + l_dof_load_56 * l_sum_366 + l_dof_load_61 * l_sum_370 + l_dof_load_66 * l_sum_373 + l_dof_load_71 * l_sum_375 + l_dof_load_76 * l_sum_377 + l_dof_load_81 * l_sum_378 + l_dof_load_86 * l_sum_379 + l_dof_load_91 * l_mult_369 + l_dof_load_96 * l_sum_380 + l_dof_load_101 * l_sum_381 + l_dof_load_106 * l_sum_382,
                l_dof_load_12 * l_sum_348 + l_dof_load_17 * l_sum_349 + l_dof_load_22 * l_sum_350 + l_dof_load_27 * l_sum_351 + l_dof_load_32 * l_sum_354 + l_dof_load_37 * l_sum_356 + l_dof_load_42 * l_sum_359 + l_dof_load_47 * l_sum_361 + l_dof_load_52 * l_sum_364 + l_dof_load_57 * l_sum_366 + l_dof_load_62 * l_sum_370 + l_dof_load_67 * l_sum_373 + l_dof_load_72 * l_sum_375 + l_dof_load_77 * l_sum_377 + l_dof_load_82 * l_sum_378 + l_dof_load_87 * l_sum_379 + l_dof_load_92 * l_mult_369 + l_dof_load_97 * l_sum_380 + l_dof_load_102 * l_sum_381 + l_dof_load_107 * l_sum_382,
                l_dof_load_13 * l_sum_348 + l_dof_load_18 * l_sum_349 + l_dof_load_23 * l_sum_350 + l_dof_load_28 * l_sum_351 + l_dof_load_33 * l_sum_354 + l_dof_load_38 * l_sum_356 + l_dof_load_43 * l_sum_359 + l_dof_load_48 * l_sum_361 + l_dof_load_53 * l_sum_364 + l_dof_load_58 * l_sum_366 + l_dof_load_63 * l_sum_370 + l_dof_load_68 * l_sum_373 + l_dof_load_73 * l_sum_375 + l_dof_load_78 * l_sum_377 + l_dof_load_83 * l_sum_378 + l_dof_load_88 * l_sum_379 + l_dof_load_93 * l_mult_369 + l_dof_load_98 * l_sum_380 + l_dof_load_103 * l_sum_381 + l_dof_load_108 * l_sum_382) - vload3(
                p_pos_1.addr(0));
            vec3 v_384 = vcons3(
                vdot3(
                    vcons3(
                        (l_r_289 + l_r_284 + l_r_285 + l_r_286 + 0.1e1 * l_r_271 + l_r_287 + -0.1e1 * l_r_277 + l_r_288) / l_op1_e3_l_21_310,
                        (l_r_296 + l_r_291 + l_r_292 + l_r_293 + 0.1e1 * l_r_273 + l_r_294 + -0.1e1 * l_r_279 + l_r_295) / l_op1_e3_l_21_310,
                        (l_r_303 + l_r_298 + l_r_299 + l_r_300 + 0.1e1 * l_r_275 + l_r_301 + -0.1e1 * l_r_281 + l_r_302) / l_op1_e3_l_21_310),
                    v_383),
                vdot3(
                    vcons3(
                        (l_r_289 + -0.1e1 * l_r_265 + l_r_285 + l_r_286 + l_r_304 + 0.1e1 * l_r_276 + l_r_305 + l_r_288) / l_op1_e3_l_21_310,
                        (l_r_296 + -0.1e1 * l_r_267 + l_r_292 + l_r_293 + l_r_306 + 0.1e1 * l_r_278 + l_r_307 + l_r_295) / l_op1_e3_l_21_310,
                        (l_r_303 + -0.1e1 * l_r_269 + l_r_299 + l_r_300 + l_r_308 + 0.1e1 * l_r_280 + l_r_309 + l_r_302) / l_op1_e3_l_21_310),
                    v_383),
                vdot3(
                    vcons3(
                        (l_r_283 + 0.1e1 * l_r_264 + l_r_284 + -0.1e1 * l_r_270 + l_r_286 + l_r_304 + l_r_287 + l_r_305 + l_r_288) / l_op1_e3_l_21_310,
                        (l_r_290 + 0.1e1 * l_r_266 + l_r_291 + -0.1e1 * l_r_272 + l_r_293 + l_r_306 + l_r_294 + l_r_307 + l_r_295) / l_op1_e3_l_21_310,
                        (l_r_297 + 0.1e1 * l_r_268 + l_r_298 + -0.1e1 * l_r_274 + l_r_300 + l_r_308 + l_r_301 + l_r_309 + l_r_302) / l_op1_e3_l_21_310),
                    v_383));
            vec3 v_385 = v_6 - v_384;
            vec3 v_386 = v_385;
            if (0.1e-7 * 0.1e-7 >= vdot3(v_384, v_384)) {
                vec3 v_387 = vcons3(0.1e-8, 0.1e-8, 0.1e-8) + v_386;
                if (0.1e1 + 0.1e-8 > vdot3(vcons3(0.1e1, 0.1e1, 0.1e1), v_386) && (v_387[0] > -0.e0 && (v_387[1] > -0.e0 && v_387[2] > -0.e0))) {
                    tensor_3 _arg_388;
                    vpack3(_arg_388, v_386);
                    return allBuild(p_mesh_0, i_cellItter_3, _arg_388, p_pos_1, true, true);
                }
            }
            v_6 = v_386;
        }
    }
    return invalidBuild(p_mesh_0);
}
static bool init_globals (world *wrld)
{
    diderot::dynseq< mesh_cell_mesh_t > l__t_389;
    globals *glob = wrld->_globals;
    l__t_389 = {};
    int32_t hi_1 = glob->gv_meshData.numCells - 1;
    for (int32_t i__t_390 = 0; i__t_390 <= hi_1; ++i__t_390) {
        l__t_389 = diderot::dynseq< mesh_cell_mesh_t >::append(l__t_389, makeFem(glob->gv_meshData, i__t_390));
    }
    glob->gv_data = glob->gv_0data0398_intermedateGlobal.loadFem(
        glob->gv_0space0396_intermedateGlobal.loadFem(glob->gv_meshData));
    return false;
}
static void normal_init (normal_strand *self, mesh_pos_mesh_t p_pos0_391, tensor_ref_3 p_xp_392)
{
    self->sv_normal[0] = 0.e0;
    self->sv_normal[1] = 0.e0;
    self->sv_normal[2] = 0.e0;
    self->sv_pos0 = p_pos0_391;
    self->sv_xp = p_xp_392;
}
static diderot::strand_status normal_update (world *wrld, globals *glob, normal_strand *self)
{
    vec3 v_3265;
    if (self->sv_pos0.valid) {
        diderot::array< double, 3 > l_evals_3263;
        diderot::array< tensor_3, 3 > l__t_3264;
        mesh_cell_mesh_t l__t_394 = makeFem(self->sv_pos0.mesh, self->sv_pos0.cell);
        func_cell_func_t l__t_395 = makeFem(glob->gv_data, l__t_394.cell);
        tensor_ref_3 l_evalPoint_396 = self->sv_pos0.refPos;
        func_t l__t_397 = (l__t_395.func);
        fns_t l__t_398 = l__t_397.space;
        int32_t l__t_399 = l__t_395.cell;
        mesh_t l__t_400 = l__t_398.mesh;
        int32_t l_mulRes_401 = l__t_399 * 20;
        int32_t t_402 = l__t_400.indexMap[l_mulRes_401];
        int32_t l_mulRes_403 = 3 * t_402;
        double l_dof_load_404 = l__t_400.coordMap[l_mulRes_403];
        double l_dof_load_405 = l__t_400.coordMap[1 + l_mulRes_403];
        double l_dof_load_406 = l__t_400.coordMap[2 + l_mulRes_403];
        int32_t t_407 = l__t_400.indexMap[l_mulRes_401 + 1];
        int32_t l_mulRes_408 = 3 * t_407;
        double l_dof_load_409 = l__t_400.coordMap[l_mulRes_408];
        double l_dof_load_410 = l__t_400.coordMap[1 + l_mulRes_408];
        double l_dof_load_411 = l__t_400.coordMap[2 + l_mulRes_408];
        int32_t t_412 = l__t_400.indexMap[l_mulRes_401 + 2];
        int32_t l_mulRes_413 = 3 * t_412;
        double l_dof_load_414 = l__t_400.coordMap[l_mulRes_413];
        double l_dof_load_415 = l__t_400.coordMap[1 + l_mulRes_413];
        double l_dof_load_416 = l__t_400.coordMap[2 + l_mulRes_413];
        int32_t t_417 = l__t_400.indexMap[l_mulRes_401 + 3];
        int32_t l_mulRes_418 = 3 * t_417;
        double l_dof_load_419 = l__t_400.coordMap[l_mulRes_418];
        double l_dof_load_420 = l__t_400.coordMap[1 + l_mulRes_418];
        double l_dof_load_421 = l__t_400.coordMap[2 + l_mulRes_418];
        int32_t t_422 = l__t_400.indexMap[l_mulRes_401 + 4];
        int32_t l_mulRes_423 = 3 * t_422;
        double l_dof_load_424 = l__t_400.coordMap[l_mulRes_423];
        double l_dof_load_425 = l__t_400.coordMap[1 + l_mulRes_423];
        double l_dof_load_426 = l__t_400.coordMap[2 + l_mulRes_423];
        int32_t t_427 = l__t_400.indexMap[l_mulRes_401 + 5];
        int32_t l_mulRes_428 = 3 * t_427;
        double l_dof_load_429 = l__t_400.coordMap[l_mulRes_428];
        double l_dof_load_430 = l__t_400.coordMap[1 + l_mulRes_428];
        double l_dof_load_431 = l__t_400.coordMap[2 + l_mulRes_428];
        int32_t t_432 = l__t_400.indexMap[l_mulRes_401 + 6];
        int32_t l_mulRes_433 = 3 * t_432;
        double l_dof_load_434 = l__t_400.coordMap[l_mulRes_433];
        double l_dof_load_435 = l__t_400.coordMap[1 + l_mulRes_433];
        double l_dof_load_436 = l__t_400.coordMap[2 + l_mulRes_433];
        int32_t t_437 = l__t_400.indexMap[l_mulRes_401 + 7];
        int32_t l_mulRes_438 = 3 * t_437;
        double l_dof_load_439 = l__t_400.coordMap[l_mulRes_438];
        double l_dof_load_440 = l__t_400.coordMap[1 + l_mulRes_438];
        double l_dof_load_441 = l__t_400.coordMap[2 + l_mulRes_438];
        int32_t t_442 = l__t_400.indexMap[l_mulRes_401 + 8];
        int32_t l_mulRes_443 = 3 * t_442;
        double l_dof_load_444 = l__t_400.coordMap[l_mulRes_443];
        double l_dof_load_445 = l__t_400.coordMap[1 + l_mulRes_443];
        double l_dof_load_446 = l__t_400.coordMap[2 + l_mulRes_443];
        int32_t t_447 = l__t_400.indexMap[l_mulRes_401 + 9];
        int32_t l_mulRes_448 = 3 * t_447;
        double l_dof_load_449 = l__t_400.coordMap[l_mulRes_448];
        double l_dof_load_450 = l__t_400.coordMap[1 + l_mulRes_448];
        double l_dof_load_451 = l__t_400.coordMap[2 + l_mulRes_448];
        int32_t t_452 = l__t_400.indexMap[l_mulRes_401 + 10];
        int32_t l_mulRes_453 = 3 * t_452;
        double l_dof_load_454 = l__t_400.coordMap[l_mulRes_453];
        double l_dof_load_455 = l__t_400.coordMap[1 + l_mulRes_453];
        double l_dof_load_456 = l__t_400.coordMap[2 + l_mulRes_453];
        int32_t t_457 = l__t_400.indexMap[l_mulRes_401 + 11];
        int32_t l_mulRes_458 = 3 * t_457;
        double l_dof_load_459 = l__t_400.coordMap[l_mulRes_458];
        double l_dof_load_460 = l__t_400.coordMap[1 + l_mulRes_458];
        double l_dof_load_461 = l__t_400.coordMap[2 + l_mulRes_458];
        int32_t t_462 = l__t_400.indexMap[l_mulRes_401 + 12];
        int32_t l_mulRes_463 = 3 * t_462;
        double l_dof_load_464 = l__t_400.coordMap[l_mulRes_463];
        double l_dof_load_465 = l__t_400.coordMap[1 + l_mulRes_463];
        double l_dof_load_466 = l__t_400.coordMap[2 + l_mulRes_463];
        int32_t t_467 = l__t_400.indexMap[l_mulRes_401 + 13];
        int32_t l_mulRes_468 = 3 * t_467;
        double l_dof_load_469 = l__t_400.coordMap[l_mulRes_468];
        double l_dof_load_470 = l__t_400.coordMap[1 + l_mulRes_468];
        double l_dof_load_471 = l__t_400.coordMap[2 + l_mulRes_468];
        int32_t t_472 = l__t_400.indexMap[l_mulRes_401 + 14];
        int32_t l_mulRes_473 = 3 * t_472;
        double l_dof_load_474 = l__t_400.coordMap[l_mulRes_473];
        double l_dof_load_475 = l__t_400.coordMap[1 + l_mulRes_473];
        double l_dof_load_476 = l__t_400.coordMap[2 + l_mulRes_473];
        int32_t t_477 = l__t_400.indexMap[l_mulRes_401 + 15];
        int32_t l_mulRes_478 = 3 * t_477;
        double l_dof_load_479 = l__t_400.coordMap[l_mulRes_478];
        double l_dof_load_480 = l__t_400.coordMap[1 + l_mulRes_478];
        double l_dof_load_481 = l__t_400.coordMap[2 + l_mulRes_478];
        int32_t t_482 = l__t_400.indexMap[l_mulRes_401 + 16];
        int32_t l_mulRes_483 = 3 * t_482;
        double l_dof_load_484 = l__t_400.coordMap[l_mulRes_483];
        double l_dof_load_485 = l__t_400.coordMap[1 + l_mulRes_483];
        double l_dof_load_486 = l__t_400.coordMap[2 + l_mulRes_483];
        int32_t t_487 = l__t_400.indexMap[l_mulRes_401 + 17];
        int32_t l_mulRes_488 = 3 * t_487;
        double l_dof_load_489 = l__t_400.coordMap[l_mulRes_488];
        double l_dof_load_490 = l__t_400.coordMap[1 + l_mulRes_488];
        double l_dof_load_491 = l__t_400.coordMap[2 + l_mulRes_488];
        int32_t t_492 = l__t_400.indexMap[l_mulRes_401 + 18];
        int32_t l_mulRes_493 = 3 * t_492;
        double l_dof_load_494 = l__t_400.coordMap[l_mulRes_493];
        double l_dof_load_495 = l__t_400.coordMap[1 + l_mulRes_493];
        double l_dof_load_496 = l__t_400.coordMap[2 + l_mulRes_493];
        int32_t t_497 = l__t_400.indexMap[l_mulRes_401 + 19];
        int32_t l_mulRes_498 = 3 * t_497;
        double l_dof_load_499 = l__t_400.coordMap[l_mulRes_498];
        double l_dof_load_500 = l__t_400.coordMap[1 + l_mulRes_498];
        double l_dof_load_501 = l__t_400.coordMap[2 + l_mulRes_498];
        double l_varAcc_502 = l_evalPoint_396[0];
        double l_varAcc_503 = l_evalPoint_396[1];
        double l_varAcc_504 = l_evalPoint_396[2];
        double l_prod2_505 = l_varAcc_502 * l_varAcc_502;
        double l_prod_506 = 0.1e1 * 0.1e1;
        double l_prod_507 = l_prod2_505 * l_prod_506;
        double l_prod_508 = l_varAcc_503 * 0.1e1;
        double l_prod_509 = l_varAcc_502 * l_prod_508;
        double l_prod_510 = 0.1e1 * l_varAcc_504;
        double l_prod_511 = l_varAcc_502 * l_prod_510;
        double l_prod_512 = l_varAcc_502 * l_prod_506;
        double l_prod2_513 = l_varAcc_503 * l_varAcc_503;
        double l_prod_514 = l_prod2_513 * 0.1e1;
        double l_prod_515 = 0.1e1 * l_prod_514;
        double l_prod_516 = l_varAcc_503 * l_varAcc_504;
        double l_prod_517 = 0.1e1 * l_prod_516;
        double l_prod_518 = 0.1e1 * l_prod_508;
        double l_prod2_519 = l_varAcc_504 * l_varAcc_504;
        double l_prod_520 = 0.1e1 * l_prod2_519;
        double l_prod_521 = 0.1e1 * l_prod_520;
        double l_prod_522 = 0.1e1 * l_prod_510;
        double l_prod_523 = 0.1e1 * l_prod_506;
        double l_mult_524 = -0.135e2 * l_prod_521;
        double l_mult_525 = -0.27e2 * l_prod_517;
        double l_mult_526 = -0.135e2 * l_prod_515;
        double l_mult_527 = -0.27e2 * l_prod_511;
        double l_mult_528 = -0.27e2 * l_prod_509;
        double l_mult_529 = -0.135e2 * l_prod_507;
        double l_sum_530 = -0.55e1 * l_prod_523 + (0.18e2 * l_prod_522 + (l_mult_524 + (0.18e2 * l_prod_518 + (l_mult_525 + (l_mult_526 + (0.18e2 * l_prod_512 + (l_mult_527 + (l_mult_528 + l_mult_529))))))));
        double l_mult_531 = 0.1e1 * l_prod_523;
        double l_mult_532 = 0.135e2 * l_prod_507;
        double l_sum_533 = l_mult_531 + (-0.9e1 * l_prod_512 + l_mult_532);
        double l_mult_534 = -0.45e1 * l_prod_522;
        double l_mult_535 = 0.27e2 * l_prod_511;
        double l_sum_536 = l_mult_534 + l_mult_535;
        double l_mult_537 = 0.135e2 * l_prod_521;
        double l_sum_538 = l_mult_534 + l_mult_537;
        double l_mult_539 = -0.45e1 * l_prod_518;
        double l_mult_540 = 0.27e2 * l_prod_509;
        double l_sum_541 = l_mult_539 + l_mult_540;
        double l_mult_542 = 0.135e2 * l_prod_515;
        double l_sum_543 = l_mult_539 + l_mult_542;
        double l_mult_544 = -0.225e2 * l_prod_522;
        double l_mult_545 = 0.27e2 * l_prod_517;
        double l_sum_546 = l_mult_544 + (0.27e2 * l_prod_521 + (l_mult_545 + l_mult_535));
        double l_mult_547 = 0.45e1 * l_prod_522;
        double l_sum_548 = l_mult_547 + l_mult_524;
        double l_mult_549 = -0.225e2 * l_prod_518;
        double l_sum_550 = l_mult_549 + (l_mult_545 + (0.27e2 * l_prod_515 + l_mult_540));
        double l_mult_551 = 0.45e1 * l_prod_518;
        double l_sum_552 = l_mult_551 + l_mult_526;
        double l_mult_553 = 0.9e1 * l_prod_523;
        double l_mult_554 = 0.54e2 * l_prod_511;
        double l_mult_555 = 0.54e2 * l_prod_509;
        double l_sum_556 = l_mult_553 + (l_mult_544 + (l_mult_537 + (l_mult_549 + (l_mult_545 + (l_mult_542 + (-0.45e2 * l_prod_512 + (l_mult_554 + (l_mult_555 + 0.405e2 * l_prod_507))))))));
        double l_mult_557 = -0.45e1 * l_prod_523;
        double l_mult_558 = 0.36e2 * l_prod_512;
        double l_sum_559 = l_mult_557 + (l_mult_547 + (l_mult_551 + (l_mult_558 + (l_mult_527 + (l_mult_528 + -0.405e2 * l_prod_507)))));
        double l_mult_560 = 0.27e2 * l_prod_522;
        double l_mult_561 = -0.27e2 * l_prod_521;
        double l_mult_562 = -0.54e2 * l_prod_511;
        double l_sum_563 = l_mult_560 + (l_mult_561 + (l_mult_525 + l_mult_562));
        double l_mult_564 = 0.27e2 * l_prod_518;
        double l_mult_565 = -0.27e2 * l_prod_515;
        double l_mult_566 = -0.54e2 * l_prod_509;
        double l_sum_567 = l_mult_564 + (l_mult_525 + (l_mult_565 + l_mult_566));
        double l_sum_568 = l_mult_531 + (-0.9e1 * l_prod_518 + l_mult_542);
        double l_sum_569 = l_mult_534 + l_mult_545;
        double l_mult_570 = -0.45e1 * l_prod_512;
        double l_sum_571 = l_mult_570 + l_mult_532;
        double l_sum_572 = l_mult_570 + l_mult_540;
        double l_mult_573 = 0.54e2 * l_prod_517;
        double l_mult_574 = -0.225e2 * l_prod_512;
        double l_sum_575 = l_mult_553 + (l_mult_544 + (l_mult_537 + (-0.45e2 * l_prod_518 + (l_mult_573 + (0.405e2 * l_prod_515 + (l_mult_574 + (l_mult_535 + (l_mult_555 + l_mult_532))))))));
        double l_mult_576 = 0.36e2 * l_prod_518;
        double l_mult_577 = 0.45e1 * l_prod_512;
        double l_sum_578 = l_mult_557 + (l_mult_547 + (l_mult_576 + (l_mult_525 + (-0.405e2 * l_prod_515 + (l_mult_577 + l_mult_528)))));
        double l_sum_579 = l_mult_574 + (l_mult_535 + (l_mult_540 + 0.27e2 * l_prod_507));
        double l_sum_580 = l_mult_577 + l_mult_529;
        double l_mult_581 = -0.54e2 * l_prod_517;
        double l_sum_582 = l_mult_560 + (l_mult_561 + (l_mult_581 + l_mult_527));
        double l_mult_583 = 0.27e2 * l_prod_512;
        double l_mult_584 = -0.27e2 * l_prod_507;
        double l_sum_585 = l_mult_583 + (l_mult_527 + (l_mult_566 + l_mult_584));
        double l_sum_586 = l_mult_531 + (-0.9e1 * l_prod_522 + l_mult_537);
        double l_sum_587 = l_mult_539 + l_mult_545;
        double l_sum_588 = l_mult_570 + l_mult_535;
        double l_sum_589 = l_mult_553 + (-0.45e2 * l_prod_522 + (0.405e2 * l_prod_521 + (l_mult_549 + (l_mult_573 + (l_mult_542 + (l_mult_574 + (l_mult_554 + (l_mult_540 + l_mult_532))))))));
        double l_mult_590 = 0.36e2 * l_prod_522;
        double l_sum_591 = l_mult_557 + (l_mult_590 + (-0.405e2 * l_prod_521 + (l_mult_551 + (l_mult_525 + (l_mult_577 + l_mult_527)))));
        double l_sum_592 = l_mult_564 + (l_mult_581 + (l_mult_565 + l_mult_528));
        double l_sum_593 = l_mult_583 + (l_mult_562 + (l_mult_528 + l_mult_584));
        double l_r_594 = l_dof_load_404 * l_sum_530;
        double l_r_595 = l_dof_load_414 * 0.e0;
        double l_r_596 = l_dof_load_419 * 0.e0;
        double l_r_597 = l_dof_load_424 * 0.e0;
        double l_r_598 = l_dof_load_429 * 0.e0;
        double l_r_599 = l_dof_load_454 * l_sum_546;
        double l_r_600 = l_dof_load_459 * l_sum_548;
        double l_r_601 = l_dof_load_464 * l_sum_550;
        double l_r_602 = l_dof_load_469 * l_sum_552;
        double l_r_603 = l_r_594 + l_dof_load_409 * l_sum_533 + l_r_595 + l_r_596 + l_r_597 + l_r_598 + l_dof_load_434 * l_sum_536 + l_dof_load_439 * l_sum_538 + l_dof_load_444 * l_sum_541 + l_dof_load_449 * l_sum_543 + l_r_599 + l_r_600 + l_r_601 + l_r_602 + l_dof_load_474 * l_sum_556 + l_dof_load_479 * l_sum_559 + l_dof_load_484 * l_mult_545 + l_dof_load_489 * l_mult_525 + l_dof_load_494 * l_sum_563 + l_dof_load_499 * l_sum_567;
        double l_r_604 = l_dof_load_409 * 0.e0;
        double l_r_605 = l_dof_load_434 * 0.e0;
        double l_r_606 = l_dof_load_439 * 0.e0;
        double l_r_607 = l_dof_load_474 * l_sum_579;
        double l_r_608 = l_dof_load_479 * l_sum_580;
        double l_r_609 = l_r_594 + l_r_604;
        double l_r_610 = l_r_609 + l_dof_load_414 * l_sum_568 + l_r_596 + l_dof_load_424 * l_sum_569 + l_dof_load_429 * l_sum_538 + l_r_605 + l_r_606 + l_dof_load_444 * l_sum_571 + l_dof_load_449 * l_sum_572 + l_r_599 + l_r_600 + l_dof_load_464 * l_sum_575 + l_dof_load_469 * l_sum_578 + l_r_607 + l_r_608 + l_dof_load_484 * l_mult_535 + l_dof_load_489 * l_sum_582 + l_dof_load_494 * l_mult_527 + l_dof_load_499 * l_sum_585;
        double l_r_611 = l_dof_load_444 * 0.e0;
        double l_r_612 = l_dof_load_449 * 0.e0;
        double l_r_613 = l_r_609 + l_r_595 + l_dof_load_419 * l_sum_586 + l_dof_load_424 * l_sum_543 + l_dof_load_429 * l_sum_587 + l_dof_load_434 * l_sum_571 + l_dof_load_439 * l_sum_588 + l_r_611 + l_r_612 + l_dof_load_454 * l_sum_589 + l_dof_load_459 * l_sum_591 + l_r_601 + l_r_602 + l_r_607 + l_r_608 + l_dof_load_484 * l_mult_540 + l_dof_load_489 * l_sum_592 + l_dof_load_494 * l_sum_593 + l_dof_load_499 * l_mult_528;
        double l_r_614 = l_dof_load_405 * l_sum_530;
        double l_r_615 = l_dof_load_415 * 0.e0;
        double l_r_616 = l_dof_load_420 * 0.e0;
        double l_r_617 = l_dof_load_425 * 0.e0;
        double l_r_618 = l_dof_load_430 * 0.e0;
        double l_r_619 = l_dof_load_455 * l_sum_546;
        double l_r_620 = l_dof_load_460 * l_sum_548;
        double l_r_621 = l_dof_load_465 * l_sum_550;
        double l_r_622 = l_dof_load_470 * l_sum_552;
        double l_r_623 = l_r_614 + l_dof_load_410 * l_sum_533 + l_r_615 + l_r_616 + l_r_617 + l_r_618 + l_dof_load_435 * l_sum_536 + l_dof_load_440 * l_sum_538 + l_dof_load_445 * l_sum_541 + l_dof_load_450 * l_sum_543 + l_r_619 + l_r_620 + l_r_621 + l_r_622 + l_dof_load_475 * l_sum_556 + l_dof_load_480 * l_sum_559 + l_dof_load_485 * l_mult_545 + l_dof_load_490 * l_mult_525 + l_dof_load_495 * l_sum_563 + l_dof_load_500 * l_sum_567;
        double l_r_624 = l_dof_load_410 * 0.e0;
        double l_r_625 = l_dof_load_435 * 0.e0;
        double l_r_626 = l_dof_load_440 * 0.e0;
        double l_r_627 = l_dof_load_475 * l_sum_579;
        double l_r_628 = l_dof_load_480 * l_sum_580;
        double l_r_629 = l_r_614 + l_r_624;
        double l_r_630 = l_r_629 + l_dof_load_415 * l_sum_568 + l_r_616 + l_dof_load_425 * l_sum_569 + l_dof_load_430 * l_sum_538 + l_r_625 + l_r_626 + l_dof_load_445 * l_sum_571 + l_dof_load_450 * l_sum_572 + l_r_619 + l_r_620 + l_dof_load_465 * l_sum_575 + l_dof_load_470 * l_sum_578 + l_r_627 + l_r_628 + l_dof_load_485 * l_mult_535 + l_dof_load_490 * l_sum_582 + l_dof_load_495 * l_mult_527 + l_dof_load_500 * l_sum_585;
        double l_r_631 = l_dof_load_445 * 0.e0;
        double l_r_632 = l_dof_load_450 * 0.e0;
        double l_r_633 = l_r_629 + l_r_615 + l_dof_load_420 * l_sum_586 + l_dof_load_425 * l_sum_543 + l_dof_load_430 * l_sum_587 + l_dof_load_435 * l_sum_571 + l_dof_load_440 * l_sum_588 + l_r_631 + l_r_632 + l_dof_load_455 * l_sum_589 + l_dof_load_460 * l_sum_591 + l_r_621 + l_r_622 + l_r_627 + l_r_628 + l_dof_load_485 * l_mult_540 + l_dof_load_490 * l_sum_592 + l_dof_load_495 * l_sum_593 + l_dof_load_500 * l_mult_528;
        double l_r_634 = l_dof_load_406 * l_sum_530;
        double l_r_635 = l_dof_load_416 * 0.e0;
        double l_r_636 = l_dof_load_421 * 0.e0;
        double l_r_637 = l_dof_load_426 * 0.e0;
        double l_r_638 = l_dof_load_431 * 0.e0;
        double l_r_639 = l_dof_load_456 * l_sum_546;
        double l_r_640 = l_dof_load_461 * l_sum_548;
        double l_r_641 = l_dof_load_466 * l_sum_550;
        double l_r_642 = l_dof_load_471 * l_sum_552;
        double l_r_643 = l_r_634 + l_dof_load_411 * l_sum_533 + l_r_635 + l_r_636 + l_r_637 + l_r_638 + l_dof_load_436 * l_sum_536 + l_dof_load_441 * l_sum_538 + l_dof_load_446 * l_sum_541 + l_dof_load_451 * l_sum_543 + l_r_639 + l_r_640 + l_r_641 + l_r_642 + l_dof_load_476 * l_sum_556 + l_dof_load_481 * l_sum_559 + l_dof_load_486 * l_mult_545 + l_dof_load_491 * l_mult_525 + l_dof_load_496 * l_sum_563 + l_dof_load_501 * l_sum_567;
        double l_r_644 = l_dof_load_411 * 0.e0;
        double l_r_645 = l_dof_load_436 * 0.e0;
        double l_r_646 = l_dof_load_441 * 0.e0;
        double l_r_647 = l_dof_load_476 * l_sum_579;
        double l_r_648 = l_dof_load_481 * l_sum_580;
        double l_r_649 = l_r_634 + l_r_644;
        double l_r_650 = l_r_649 + l_dof_load_416 * l_sum_568 + l_r_636 + l_dof_load_426 * l_sum_569 + l_dof_load_431 * l_sum_538 + l_r_645 + l_r_646 + l_dof_load_446 * l_sum_571 + l_dof_load_451 * l_sum_572 + l_r_639 + l_r_640 + l_dof_load_466 * l_sum_575 + l_dof_load_471 * l_sum_578 + l_r_647 + l_r_648 + l_dof_load_486 * l_mult_535 + l_dof_load_491 * l_sum_582 + l_dof_load_496 * l_mult_527 + l_dof_load_501 * l_sum_585;
        double l_r_651 = l_dof_load_446 * 0.e0;
        double l_r_652 = l_dof_load_451 * 0.e0;
        double l_r_653 = l_r_649 + l_r_635 + l_dof_load_421 * l_sum_586 + l_dof_load_426 * l_sum_543 + l_dof_load_431 * l_sum_587 + l_dof_load_436 * l_sum_571 + l_dof_load_441 * l_sum_588 + l_r_651 + l_r_652 + l_dof_load_456 * l_sum_589 + l_dof_load_461 * l_sum_591 + l_r_641 + l_r_642 + l_r_647 + l_r_648 + l_dof_load_486 * l_mult_540 + l_dof_load_491 * l_sum_592 + l_dof_load_496 * l_sum_593 + l_dof_load_501 * l_mult_528;
        double l_r_654 = 0.e0 * l_r_603;
        double l_r_655 = 0.e0 * l_r_623;
        double l_r_656 = 0.e0 * l_r_643;
        double l_r_657 = l_r_654 + l_r_655;
        double l_r_658 = l_r_657 + l_r_656;
        double l_r_659 = 0.e0 * l_r_610;
        double l_r_660 = 0.e0 * l_r_630;
        double l_r_661 = 0.e0 * l_r_650;
        double l_r_662 = l_r_659 + l_r_660;
        double l_r_663 = l_r_662 + l_r_661;
        double l_r_664 = 0.e0 * l_r_613;
        double l_r_665 = 0.e0 * l_r_633;
        double l_r_666 = 0.e0 * l_r_653;
        double l_r_667 = l_r_664 + l_r_665;
        double l_r_668 = l_r_667 + l_r_666;
        double l_r_669 = l_r_657 + -0.1e1 * l_r_643;
        double l_r_670 = l_r_662 + -0.1e1 * l_r_650;
        double l_r_671 = l_r_667 + -0.1e1 * l_r_653;
        double l_r_672 = l_r_654 + 0.1e1 * l_r_623 + l_r_656;
        double l_r_673 = l_r_659 + 0.1e1 * l_r_630 + l_r_661;
        double l_r_674 = l_r_664 + 0.1e1 * l_r_633 + l_r_666;
        double l_r_675 = l_r_657 + 0.1e1 * l_r_643;
        double l_r_676 = l_r_662 + 0.1e1 * l_r_650;
        double l_r_677 = l_r_667 + 0.1e1 * l_r_653;
        double l_r_678 = -0.1e1 * l_r_603 + l_r_655 + l_r_656;
        double l_r_679 = -0.1e1 * l_r_610 + l_r_660 + l_r_661;
        double l_r_680 = -0.1e1 * l_r_613 + l_r_665 + l_r_666;
        double l_r_681 = l_r_654 + -0.1e1 * l_r_623 + l_r_656;
        double l_r_682 = l_r_659 + -0.1e1 * l_r_630 + l_r_661;
        double l_r_683 = l_r_664 + -0.1e1 * l_r_633 + l_r_666;
        double l_r_684 = 0.1e1 * l_r_603 + l_r_655 + l_r_656;
        double l_r_685 = 0.1e1 * l_r_610 + l_r_660 + l_r_661;
        double l_r_686 = 0.1e1 * l_r_613 + l_r_665 + l_r_666;
        double l_r_687 = l_r_603 * l_r_663 + l_r_623 * l_r_676 + l_r_643 * l_r_682;
        double l_r_688 = l_r_603 * l_r_668 + l_r_623 * l_r_677 + l_r_643 * l_r_683;
        double l_r_689 = l_r_603 * l_r_670 + l_r_623 * l_r_663 + l_r_643 * l_r_685;
        double l_r_690 = l_r_603 * l_r_671 + l_r_623 * l_r_668 + l_r_643 * l_r_686;
        double l_r_691 = l_r_603 * l_r_673 + l_r_623 * l_r_679 + l_r_643 * l_r_663;
        double l_r_692 = l_r_603 * l_r_674 + l_r_623 * l_r_680 + l_r_643 * l_r_668;
        double l_r_693 = l_r_610 * l_r_658 + l_r_630 * l_r_675 + l_r_650 * l_r_681;
        double l_r_694 = l_r_610 * l_r_668 + l_r_630 * l_r_677 + l_r_650 * l_r_683;
        double l_r_695 = l_r_610 * l_r_669 + l_r_630 * l_r_658 + l_r_650 * l_r_684;
        double l_r_696 = l_r_610 * l_r_671 + l_r_630 * l_r_668 + l_r_650 * l_r_686;
        double l_r_697 = l_r_610 * l_r_672 + l_r_630 * l_r_678 + l_r_650 * l_r_658;
        double l_r_698 = l_r_610 * l_r_674 + l_r_630 * l_r_680 + l_r_650 * l_r_668;
        double l_r_699 = l_r_613 * l_r_658 + l_r_633 * l_r_675 + l_r_653 * l_r_681;
        double l_r_700 = l_r_613 * l_r_663 + l_r_633 * l_r_676 + l_r_653 * l_r_682;
        double l_r_701 = l_r_613 * l_r_669 + l_r_633 * l_r_658 + l_r_653 * l_r_684;
        double l_r_702 = l_r_613 * l_r_670 + l_r_633 * l_r_663 + l_r_653 * l_r_685;
        double l_r_703 = l_r_613 * l_r_672 + l_r_633 * l_r_678 + l_r_653 * l_r_658;
        double l_r_704 = l_r_613 * l_r_673 + l_r_633 * l_r_679 + l_r_653 * l_r_663;
        vec3 v_705 = vcons3(l_r_610, l_r_630, l_r_650);
        double l_vdot_706 = vdot3(v_705, vcons3(l_r_668, l_r_677, l_r_683));
        double l_vdot_707 = vdot3(v_705, vcons3(l_r_671, l_r_668, l_r_686));
        double l_vdot_708 = vdot3(v_705, vcons3(l_r_674, l_r_680, l_r_668));
        double l_op1_e3_l_16_709 = vdot3(vcons3(l_r_603, l_r_623, l_r_643), vcons3(l_vdot_706, l_vdot_707, l_vdot_708));
        double l_r_710 = 0.e0 * (l_r_603 * l_r_658 + l_r_623 * l_r_675 + l_r_643 * l_r_681);
        double l_r_711 = 0.e0 * l_r_688;
        double l_r_712 = 0.e0 * l_r_693;
        double l_r_713 = 0.e0 * (l_r_610 * l_r_663 + l_r_630 * l_r_676 + l_r_650 * l_r_682);
        double l_r_714 = 0.e0 * l_r_699;
        double l_r_715 = 0.e0 * (l_r_613 * l_r_668 + l_r_633 * l_r_677 + l_r_653 * l_r_683);
        double l_r_716 = l_r_710 + 0.e0 * l_r_687;
        double l_r_717 = l_r_716 + l_r_711 + l_r_712 + l_r_713 + 0.1e1 * l_r_694 + l_r_714 + -0.1e1 * l_r_700 + l_r_715;
        double l_r_718 = 0.e0 * (l_r_603 * l_r_669 + l_r_623 * l_r_658 + l_r_643 * l_r_684);
        double l_r_719 = 0.e0 * l_r_690;
        double l_r_720 = 0.e0 * l_r_695;
        double l_r_721 = 0.e0 * (l_r_610 * l_r_670 + l_r_630 * l_r_663 + l_r_650 * l_r_685);
        double l_r_722 = 0.e0 * l_r_701;
        double l_r_723 = 0.e0 * (l_r_613 * l_r_671 + l_r_633 * l_r_668 + l_r_653 * l_r_686);
        double l_r_724 = l_r_718 + 0.e0 * l_r_689;
        double l_r_725 = l_r_724 + l_r_719 + l_r_720 + l_r_721 + 0.1e1 * l_r_696 + l_r_722 + -0.1e1 * l_r_702 + l_r_723;
        double l_r_726 = 0.e0 * (l_r_603 * l_r_672 + l_r_623 * l_r_678 + l_r_643 * l_r_658);
        double l_r_727 = 0.e0 * l_r_692;
        double l_r_728 = 0.e0 * l_r_697;
        double l_r_729 = 0.e0 * (l_r_610 * l_r_673 + l_r_630 * l_r_679 + l_r_650 * l_r_663);
        double l_r_730 = 0.e0 * l_r_703;
        double l_r_731 = 0.e0 * (l_r_613 * l_r_674 + l_r_633 * l_r_680 + l_r_653 * l_r_668);
        double l_r_732 = l_r_726 + 0.e0 * l_r_691;
        double l_r_733 = l_r_732 + l_r_727 + l_r_728 + l_r_729 + 0.1e1 * l_r_698 + l_r_730 + -0.1e1 * l_r_704 + l_r_731;
        double l_r_734 = 0.e0 * l_r_694;
        double l_r_735 = 0.e0 * l_r_700;
        double l_r_736 = l_r_716 + -0.1e1 * l_r_688 + l_r_712 + l_r_713 + l_r_734 + 0.1e1 * l_r_699 + l_r_735 + l_r_715;
        double l_r_737 = 0.e0 * l_r_696;
        double l_r_738 = 0.e0 * l_r_702;
        double l_r_739 = l_r_724 + -0.1e1 * l_r_690 + l_r_720 + l_r_721 + l_r_737 + 0.1e1 * l_r_701 + l_r_738 + l_r_723;
        double l_r_740 = 0.e0 * l_r_698;
        double l_r_741 = 0.e0 * l_r_704;
        double l_r_742 = l_r_732 + -0.1e1 * l_r_692 + l_r_728 + l_r_729 + l_r_740 + 0.1e1 * l_r_703 + l_r_741 + l_r_731;
        double l_r_743 = l_r_710 + 0.1e1 * l_r_687 + l_r_711 + -0.1e1 * l_r_693 + l_r_713 + l_r_734 + l_r_714 + l_r_735 + l_r_715;
        double l_r_744 = l_r_718 + 0.1e1 * l_r_689 + l_r_719 + -0.1e1 * l_r_695 + l_r_721 + l_r_737 + l_r_722 + l_r_738 + l_r_723;
        double l_r_745 = l_r_726 + 0.1e1 * l_r_691 + l_r_727 + -0.1e1 * l_r_697 + l_r_729 + l_r_740 + l_r_730 + l_r_741 + l_r_731;
        double l_op1_e3_l_18_746 = 0.2e1 * l_op1_e3_l_16_709;
        int32_t l_mulRes_747 = l__t_399 * 84;
        int32_t t_748 = l__t_398.indexMap[l_mulRes_747];
        int32_t t_749 = l__t_398.indexMap[l_mulRes_747 + 1];
        int32_t t_750 = l__t_398.indexMap[l_mulRes_747 + 2];
        int32_t t_751 = l__t_398.indexMap[l_mulRes_747 + 3];
        int32_t t_752 = l__t_398.indexMap[l_mulRes_747 + 4];
        int32_t t_753 = l__t_398.indexMap[l_mulRes_747 + 5];
        int32_t t_754 = l__t_398.indexMap[l_mulRes_747 + 6];
        int32_t t_755 = l__t_398.indexMap[l_mulRes_747 + 7];
        int32_t t_756 = l__t_398.indexMap[l_mulRes_747 + 8];
        int32_t t_757 = l__t_398.indexMap[l_mulRes_747 + 9];
        int32_t t_758 = l__t_398.indexMap[l_mulRes_747 + 10];
        int32_t t_759 = l__t_398.indexMap[l_mulRes_747 + 11];
        int32_t t_760 = l__t_398.indexMap[l_mulRes_747 + 12];
        int32_t t_761 = l__t_398.indexMap[l_mulRes_747 + 13];
        int32_t t_762 = l__t_398.indexMap[l_mulRes_747 + 14];
        int32_t t_763 = l__t_398.indexMap[l_mulRes_747 + 15];
        int32_t t_764 = l__t_398.indexMap[l_mulRes_747 + 16];
        int32_t t_765 = l__t_398.indexMap[l_mulRes_747 + 17];
        int32_t t_766 = l__t_398.indexMap[l_mulRes_747 + 18];
        int32_t t_767 = l__t_398.indexMap[l_mulRes_747 + 19];
        int32_t t_768 = l__t_398.indexMap[l_mulRes_747 + 20];
        int32_t t_769 = l__t_398.indexMap[l_mulRes_747 + 21];
        int32_t t_770 = l__t_398.indexMap[l_mulRes_747 + 22];
        int32_t t_771 = l__t_398.indexMap[l_mulRes_747 + 23];
        int32_t t_772 = l__t_398.indexMap[l_mulRes_747 + 24];
        int32_t t_773 = l__t_398.indexMap[l_mulRes_747 + 25];
        int32_t t_774 = l__t_398.indexMap[l_mulRes_747 + 26];
        int32_t t_775 = l__t_398.indexMap[l_mulRes_747 + 27];
        int32_t t_776 = l__t_398.indexMap[l_mulRes_747 + 28];
        int32_t t_777 = l__t_398.indexMap[l_mulRes_747 + 29];
        int32_t t_778 = l__t_398.indexMap[l_mulRes_747 + 30];
        int32_t t_779 = l__t_398.indexMap[l_mulRes_747 + 31];
        int32_t t_780 = l__t_398.indexMap[l_mulRes_747 + 32];
        int32_t t_781 = l__t_398.indexMap[l_mulRes_747 + 33];
        int32_t t_782 = l__t_398.indexMap[l_mulRes_747 + 34];
        int32_t t_783 = l__t_398.indexMap[l_mulRes_747 + 35];
        int32_t t_784 = l__t_398.indexMap[l_mulRes_747 + 36];
        int32_t t_785 = l__t_398.indexMap[l_mulRes_747 + 37];
        int32_t t_786 = l__t_398.indexMap[l_mulRes_747 + 38];
        int32_t t_787 = l__t_398.indexMap[l_mulRes_747 + 39];
        int32_t t_788 = l__t_398.indexMap[l_mulRes_747 + 40];
        int32_t t_789 = l__t_398.indexMap[l_mulRes_747 + 41];
        int32_t t_790 = l__t_398.indexMap[l_mulRes_747 + 42];
        int32_t t_791 = l__t_398.indexMap[l_mulRes_747 + 43];
        int32_t t_792 = l__t_398.indexMap[l_mulRes_747 + 44];
        int32_t t_793 = l__t_398.indexMap[l_mulRes_747 + 45];
        int32_t t_794 = l__t_398.indexMap[l_mulRes_747 + 46];
        int32_t t_795 = l__t_398.indexMap[l_mulRes_747 + 47];
        int32_t t_796 = l__t_398.indexMap[l_mulRes_747 + 48];
        int32_t t_797 = l__t_398.indexMap[l_mulRes_747 + 49];
        int32_t t_798 = l__t_398.indexMap[l_mulRes_747 + 50];
        int32_t t_799 = l__t_398.indexMap[l_mulRes_747 + 51];
        int32_t t_800 = l__t_398.indexMap[l_mulRes_747 + 52];
        int32_t t_801 = l__t_398.indexMap[l_mulRes_747 + 53];
        int32_t t_802 = l__t_398.indexMap[l_mulRes_747 + 54];
        int32_t t_803 = l__t_398.indexMap[l_mulRes_747 + 55];
        int32_t t_804 = l__t_398.indexMap[l_mulRes_747 + 56];
        int32_t t_805 = l__t_398.indexMap[l_mulRes_747 + 57];
        int32_t t_806 = l__t_398.indexMap[l_mulRes_747 + 58];
        int32_t t_807 = l__t_398.indexMap[l_mulRes_747 + 59];
        int32_t t_808 = l__t_398.indexMap[l_mulRes_747 + 60];
        int32_t t_809 = l__t_398.indexMap[l_mulRes_747 + 61];
        int32_t t_810 = l__t_398.indexMap[l_mulRes_747 + 62];
        int32_t t_811 = l__t_398.indexMap[l_mulRes_747 + 63];
        int32_t t_812 = l__t_398.indexMap[l_mulRes_747 + 64];
        int32_t t_813 = l__t_398.indexMap[l_mulRes_747 + 65];
        int32_t t_814 = l__t_398.indexMap[l_mulRes_747 + 66];
        int32_t t_815 = l__t_398.indexMap[l_mulRes_747 + 67];
        int32_t t_816 = l__t_398.indexMap[l_mulRes_747 + 68];
        int32_t t_817 = l__t_398.indexMap[l_mulRes_747 + 69];
        int32_t t_818 = l__t_398.indexMap[l_mulRes_747 + 70];
        int32_t t_819 = l__t_398.indexMap[l_mulRes_747 + 71];
        int32_t t_820 = l__t_398.indexMap[l_mulRes_747 + 72];
        int32_t t_821 = l__t_398.indexMap[l_mulRes_747 + 73];
        int32_t t_822 = l__t_398.indexMap[l_mulRes_747 + 74];
        int32_t t_823 = l__t_398.indexMap[l_mulRes_747 + 75];
        int32_t t_824 = l__t_398.indexMap[l_mulRes_747 + 76];
        int32_t t_825 = l__t_398.indexMap[l_mulRes_747 + 77];
        int32_t t_826 = l__t_398.indexMap[l_mulRes_747 + 78];
        int32_t t_827 = l__t_398.indexMap[l_mulRes_747 + 79];
        int32_t t_828 = l__t_398.indexMap[l_mulRes_747 + 80];
        int32_t t_829 = l__t_398.indexMap[l_mulRes_747 + 81];
        int32_t t_830 = l__t_398.indexMap[l_mulRes_747 + 82];
        int32_t t_831 = l__t_398.indexMap[l_mulRes_747 + 83];
        double t_832 = l__t_397.coordMap[1 * t_831];
        double t_833 = l__t_397.coordMap[1 * t_830];
        double t_834 = l__t_397.coordMap[1 * t_829];
        double t_835 = l__t_397.coordMap[1 * t_828];
        double t_836 = l__t_397.coordMap[1 * t_827];
        double t_837 = l__t_397.coordMap[1 * t_826];
        double t_838 = l__t_397.coordMap[1 * t_825];
        double t_839 = l__t_397.coordMap[1 * t_824];
        double t_840 = l__t_397.coordMap[1 * t_823];
        double t_841 = l__t_397.coordMap[1 * t_822];
        double t_842 = l__t_397.coordMap[1 * t_821];
        double t_843 = l__t_397.coordMap[1 * t_820];
        double t_844 = l__t_397.coordMap[1 * t_819];
        double t_845 = l__t_397.coordMap[1 * t_818];
        double t_846 = l__t_397.coordMap[1 * t_817];
        double t_847 = l__t_397.coordMap[1 * t_816];
        double t_848 = l__t_397.coordMap[1 * t_815];
        double t_849 = l__t_397.coordMap[1 * t_814];
        double t_850 = l__t_397.coordMap[1 * t_813];
        double t_851 = l__t_397.coordMap[1 * t_812];
        double t_852 = l__t_397.coordMap[1 * t_811];
        double t_853 = l__t_397.coordMap[1 * t_810];
        double t_854 = l__t_397.coordMap[1 * t_809];
        double t_855 = l__t_397.coordMap[1 * t_808];
        double t_856 = l__t_397.coordMap[1 * t_807];
        double t_857 = l__t_397.coordMap[1 * t_806];
        double t_858 = l__t_397.coordMap[1 * t_805];
        double t_859 = l__t_397.coordMap[1 * t_804];
        double t_860 = l__t_397.coordMap[1 * t_803];
        double t_861 = l__t_397.coordMap[1 * t_802];
        double t_862 = l__t_397.coordMap[1 * t_801];
        double t_863 = l__t_397.coordMap[1 * t_800];
        double t_864 = l__t_397.coordMap[1 * t_799];
        double t_865 = l__t_397.coordMap[1 * t_798];
        double t_866 = l__t_397.coordMap[1 * t_797];
        double t_867 = l__t_397.coordMap[1 * t_796];
        double t_868 = l__t_397.coordMap[1 * t_795];
        double t_869 = l__t_397.coordMap[1 * t_794];
        double t_870 = l__t_397.coordMap[1 * t_793];
        double t_871 = l__t_397.coordMap[1 * t_792];
        double t_872 = l__t_397.coordMap[1 * t_791];
        double t_873 = l__t_397.coordMap[1 * t_790];
        double t_874 = l__t_397.coordMap[1 * t_789];
        double t_875 = l__t_397.coordMap[1 * t_788];
        double t_876 = l__t_397.coordMap[1 * t_787];
        double t_877 = l__t_397.coordMap[1 * t_786];
        double t_878 = l__t_397.coordMap[1 * t_785];
        double t_879 = l__t_397.coordMap[1 * t_784];
        double t_880 = l__t_397.coordMap[1 * t_783];
        double t_881 = l__t_397.coordMap[1 * t_782];
        double t_882 = l__t_397.coordMap[1 * t_781];
        double t_883 = l__t_397.coordMap[1 * t_780];
        double t_884 = l__t_397.coordMap[1 * t_779];
        double t_885 = l__t_397.coordMap[1 * t_778];
        double t_886 = l__t_397.coordMap[1 * t_777];
        double t_887 = l__t_397.coordMap[1 * t_776];
        double t_888 = l__t_397.coordMap[1 * t_775];
        double t_889 = l__t_397.coordMap[1 * t_774];
        double t_890 = l__t_397.coordMap[1 * t_773];
        double t_891 = l__t_397.coordMap[1 * t_772];
        double t_892 = l__t_397.coordMap[1 * t_771];
        double t_893 = l__t_397.coordMap[1 * t_770];
        double t_894 = l__t_397.coordMap[1 * t_769];
        double t_895 = l__t_397.coordMap[1 * t_768];
        double t_896 = l__t_397.coordMap[1 * t_767];
        double t_897 = l__t_397.coordMap[1 * t_766];
        double t_898 = l__t_397.coordMap[1 * t_765];
        double t_899 = l__t_397.coordMap[1 * t_764];
        double t_900 = l__t_397.coordMap[1 * t_763];
        double t_901 = l__t_397.coordMap[1 * t_762];
        double t_902 = l__t_397.coordMap[1 * t_761];
        double t_903 = l__t_397.coordMap[1 * t_760];
        double t_904 = l__t_397.coordMap[1 * t_759];
        double t_905 = l__t_397.coordMap[1 * t_758];
        double t_906 = l__t_397.coordMap[1 * t_757];
        double t_907 = l__t_397.coordMap[1 * t_756];
        double t_908 = l__t_397.coordMap[1 * t_755];
        double t_909 = l__t_397.coordMap[1 * t_754];
        double t_910 = l__t_397.coordMap[1 * t_753];
        double t_911 = l__t_397.coordMap[1 * t_752];
        double t_912 = l__t_397.coordMap[1 * t_751];
        double t_913 = l__t_397.coordMap[1 * t_750];
        double t_914 = l__t_397.coordMap[1 * t_749];
        double t_915 = l__t_397.coordMap[1 * t_748];
        vec4 v_916 = vcons4(t_915, t_914, t_913, t_912);
        vec4 v_917 = vcons4(t_911, t_910, t_909, t_908);
        vec4 v_918 = vcons4(t_907, t_906, t_905, t_904);
        vec4 v_919 = vcons4(t_903, t_902, t_901, t_900);
        vec4 v_920 = vcons4(t_899, t_898, t_897, t_896);
        vec4 v_921 = vcons4(t_895, t_894, t_893, t_892);
        vec4 v_922 = vcons4(t_891, t_890, t_889, t_888);
        vec4 v_923 = vcons4(t_887, t_886, t_885, t_884);
        vec4 v_924 = vcons4(t_883, t_882, t_881, t_880);
        vec4 v_925 = vcons4(t_879, t_878, t_877, t_876);
        vec4 v_926 = vcons4(t_875, t_874, t_873, t_872);
        vec4 v_927 = vcons4(t_871, t_870, t_869, t_868);
        vec4 v_928 = vcons4(t_867, t_866, t_865, t_864);
        vec4 v_929 = vcons4(t_863, t_862, t_861, t_860);
        vec4 v_930 = vcons4(t_859, t_858, t_857, t_856);
        vec4 v_931 = vcons4(t_855, t_854, t_853, t_852);
        vec4 v_932 = vcons4(t_851, t_850, t_849, t_848);
        vec4 v_933 = vcons4(t_847, t_846, t_845, t_844);
        vec4 v_934 = vcons4(t_843, t_842, t_841, t_840);
        vec4 v_935 = vcons4(t_839, t_838, t_837, t_836);
        vec4 v_936 = vcons4(t_835, t_834, t_833, t_832);
        double l_prod3_937 = l_prod2_505 * l_varAcc_502;
        double l_prod4_938 = l_prod3_937 * l_varAcc_502;
        double l_prod_939 = l_prod4_938 * l_prod_506;
        double l_prod_940 = l_prod3_937 * l_prod_508;
        double l_prod_941 = l_prod3_937 * l_prod_510;
        double l_prod_942 = l_prod3_937 * l_prod_506;
        double l_prod_943 = l_prod2_505 * l_prod_514;
        double l_prod_944 = l_prod2_505 * l_prod_516;
        double l_prod_945 = l_prod2_505 * l_prod_508;
        double l_prod_946 = l_prod2_505 * l_prod_520;
        double l_prod_947 = l_prod2_505 * l_prod_510;
        double l_prod3_948 = l_prod2_513 * l_varAcc_503;
        double l_prod_949 = l_prod3_948 * 0.1e1;
        double l_prod_950 = l_varAcc_502 * l_prod_949;
        double l_prod_951 = l_prod2_513 * l_varAcc_504;
        double l_prod_952 = l_varAcc_502 * l_prod_951;
        double l_prod_953 = l_varAcc_502 * l_prod_514;
        double l_prod_954 = l_varAcc_503 * l_prod2_519;
        double l_prod_955 = l_varAcc_502 * l_prod_954;
        double l_prod_956 = l_varAcc_502 * l_prod_516;
        double l_prod3_957 = l_prod2_519 * l_varAcc_504;
        double l_prod_958 = 0.1e1 * l_prod3_957;
        double l_prod_959 = l_varAcc_502 * l_prod_958;
        double l_prod_960 = l_varAcc_502 * l_prod_520;
        double l_prod4_961 = l_prod3_948 * l_varAcc_503;
        double l_prod_962 = l_prod4_961 * 0.1e1;
        double l_prod_963 = 0.1e1 * l_prod_962;
        double l_prod_964 = l_prod3_948 * l_varAcc_504;
        double l_prod_965 = 0.1e1 * l_prod_964;
        double l_prod_966 = 0.1e1 * l_prod_949;
        double l_prod_967 = l_prod2_513 * l_prod2_519;
        double l_prod_968 = 0.1e1 * l_prod_967;
        double l_prod_969 = 0.1e1 * l_prod_951;
        double l_prod_970 = l_varAcc_503 * l_prod3_957;
        double l_prod_971 = 0.1e1 * l_prod_970;
        double l_prod_972 = 0.1e1 * l_prod_954;
        double l_prod4_973 = l_prod3_957 * l_varAcc_504;
        double l_prod_974 = 0.1e1 * l_prod4_973;
        double l_prod_975 = 0.1e1 * l_prod_974;
        double l_prod_976 = 0.1e1 * l_prod_958;
        double l_mult_977 = 0.1944e4 * l_prod_975;
        double l_mult_978 = 0.7776e4 * l_prod_971;
        double l_mult_979 = 0.11664e5 * l_prod_968;
        double l_mult_980 = 0.7776e4 * l_prod_965;
        double l_mult_981 = 0.1944e4 * l_prod_963;
        double l_mult_982 = 0.7776e4 * l_prod_959;
        double l_mult_983 = 0.23328e5 * l_prod_955;
        double l_mult_984 = 0.23328e5 * l_prod_952;
        double l_mult_985 = 0.7776e4 * l_prod_950;
        double l_mult_986 = 0.11664e5 * l_prod_946;
        double l_mult_987 = 0.23328e5 * l_prod_944;
        double l_mult_988 = 0.11664e5 * l_prod_943;
        double l_mult_989 = 0.7776e4 * l_prod_941;
        double l_mult_990 = 0.7776e4 * l_prod_940;
        double l_mult_991 = 0.1944e4 * l_prod_939;
        double l_sum_992 = l_mult_990 + l_mult_991;
        double l_sum_993 = 0.1624e3 * l_prod_523 + (-0.1323e4 * l_prod_522 + (0.3780e4 * l_prod_521 + (-0.4536e4 * l_prod_976 + (l_mult_977 + (-0.1323e4 * l_prod_518 + (0.7560e4 * l_prod_517 + (-0.13608e5 * l_prod_972 + (l_mult_978 + (0.3780e4 * l_prod_515 + (-0.13608e5 * l_prod_969 + (l_mult_979 + (-0.4536e4 * l_prod_966 + (l_mult_980 + (l_mult_981 + (-0.1323e4 * l_prod_512 + (0.7560e4 * l_prod_511 + (-0.13608e5 * l_prod_960 + (l_mult_982 + (0.7560e4 * l_prod_509 + (-0.27216e5 * l_prod_956 + (l_mult_983 + (-0.13608e5 * l_prod_953 + (l_mult_984 + (l_mult_985 + (0.3780e4 * l_prod_507 + (-0.13608e5 * l_prod_947 + (l_mult_986 + (-0.13608e5 * l_prod_945 + (l_mult_987 + (l_mult_988 + (-0.4536e4 * l_prod_942 + (l_mult_989 + l_sum_992))))))))))))))))))))))))))))))));
        double l_mult_994 = 0.274e2 * l_prod_523;
        double l_mult_995 = -0.180e3 * l_prod_522;
        double l_mult_996 = 0.2268e4 * l_prod_511;
        double l_mult_997 = -0.7776e4 * l_prod_947;
        double l_mult_998 = -0.99e2 * l_prod_522;
        double l_mult_999 = 0.594e3 * l_prod_521;
        double l_mult_1000 = 0.972e3 * l_prod_511;
        double l_mult_1001 = -0.5832e4 * l_prod_960;
        double l_mult_1002 = -0.1944e4 * l_prod_947;
        double l_sum_1003 = l_mult_1002 + l_mult_986;
        double l_mult_1004 = -0.72e2 * l_prod_522;
        double l_mult_1005 = 0.648e3 * l_prod_521;
        double l_mult_1006 = -0.1296e4 * l_prod_976;
        double l_mult_1007 = 0.432e3 * l_prod_511;
        double l_mult_1008 = -0.3888e4 * l_prod_960;
        double l_sum_1009 = l_mult_1007 + (l_mult_1008 + l_mult_982);
        double l_mult_1010 = -0.54e2 * l_prod_522;
        double l_mult_1011 = -0.1944e4 * l_prod_976;
        double l_sum_1012 = l_mult_1010 + (l_mult_999 + (l_mult_1011 + l_mult_977));
        double l_mult_1013 = -0.180e3 * l_prod_518;
        double l_mult_1014 = 0.2268e4 * l_prod_509;
        double l_mult_1015 = -0.7776e4 * l_prod_945;
        double l_mult_1016 = -0.99e2 * l_prod_518;
        double l_mult_1017 = 0.594e3 * l_prod_515;
        double l_mult_1018 = 0.972e3 * l_prod_509;
        double l_mult_1019 = -0.5832e4 * l_prod_953;
        double l_mult_1020 = -0.1944e4 * l_prod_945;
        double l_sum_1021 = l_mult_1020 + l_mult_988;
        double l_mult_1022 = -0.72e2 * l_prod_518;
        double l_mult_1023 = 0.648e3 * l_prod_515;
        double l_mult_1024 = -0.1296e4 * l_prod_966;
        double l_mult_1025 = 0.432e3 * l_prod_509;
        double l_mult_1026 = -0.3888e4 * l_prod_953;
        double l_sum_1027 = l_mult_1025 + (l_mult_1026 + l_mult_985);
        double l_mult_1028 = -0.54e2 * l_prod_518;
        double l_mult_1029 = -0.1944e4 * l_prod_966;
        double l_sum_1030 = l_mult_1028 + (l_mult_1017 + (l_mult_1029 + l_mult_981));
        double l_mult_1031 = 0.2088e4 * l_prod_522;
        double l_mult_1032 = -0.10044e5 * l_prod_521;
        double l_mult_1033 = 0.15552e5 * l_prod_976;
        double l_mult_1034 = -0.7776e4 * l_prod_975;
        double l_mult_1035 = -0.10044e5 * l_prod_517;
        double l_mult_1036 = 0.31104e5 * l_prod_972;
        double l_mult_1037 = -0.23328e5 * l_prod_971;
        double l_mult_1038 = 0.15552e5 * l_prod_969;
        double l_mult_1039 = -0.23328e5 * l_prod_968;
        double l_mult_1040 = -0.7776e4 * l_prod_965;
        double l_mult_1041 = -0.10044e5 * l_prod_511;
        double l_mult_1042 = 0.31104e5 * l_prod_960;
        double l_mult_1043 = -0.23328e5 * l_prod_959;
        double l_mult_1044 = 0.31104e5 * l_prod_956;
        double l_mult_1045 = -0.46656e5 * l_prod_955;
        double l_mult_1046 = -0.23328e5 * l_prod_952;
        double l_mult_1047 = 0.15552e5 * l_prod_947;
        double l_mult_1048 = -0.23328e5 * l_prod_946;
        double l_mult_1049 = -0.23328e5 * l_prod_944;
        double l_mult_1050 = -0.7776e4 * l_prod_941;
        double l_sum_1051 = l_mult_1031 + (l_mult_1032 + (l_mult_1033 + (l_mult_1034 + (l_mult_1035 + (l_mult_1036 + (l_mult_1037 + (l_mult_1038 + (l_mult_1039 + (l_mult_1040 + (l_mult_1041 + (l_mult_1042 + (l_mult_1043 + (l_mult_1044 + (l_mult_1045 + (l_mult_1046 + (l_mult_1047 + (l_mult_1048 + (l_mult_1049 + l_mult_1050))))))))))))))))));
        double l_mult_1052 = -0.1071e4 * l_prod_522;
        double l_mult_1053 = 0.9342e4 * l_prod_521;
        double l_mult_1054 = 0.11664e5 * l_prod_975;
        double l_mult_1055 = 0.2916e4 * l_prod_517;
        double l_mult_1056 = -0.21384e5 * l_prod_972;
        double l_mult_1057 = 0.23328e5 * l_prod_971;
        double l_mult_1058 = -0.1944e4 * l_prod_969;
        double l_mult_1059 = 0.2916e4 * l_prod_511;
        double l_mult_1060 = -0.21384e5 * l_prod_960;
        double l_mult_1061 = 0.23328e5 * l_prod_959;
        double l_mult_1062 = -0.3888e4 * l_prod_956;
        double l_sum_1063 = l_mult_1052 + (l_mult_1053 + (-0.19440e5 * l_prod_976 + (l_mult_1054 + (l_mult_1055 + (l_mult_1056 + (l_mult_1057 + (l_mult_1058 + (l_mult_979 + (l_mult_1059 + (l_mult_1060 + (l_mult_1061 + (l_mult_1062 + (l_mult_983 + l_sum_1003)))))))))))));
        double l_mult_1064 = 0.360e3 * l_prod_522;
        double l_mult_1065 = -0.3672e4 * l_prod_521;
        double l_mult_1066 = 0.10368e5 * l_prod_976;
        double l_mult_1067 = -0.432e3 * l_prod_517;
        double l_mult_1068 = 0.3888e4 * l_prod_972;
        double l_mult_1069 = -0.7776e4 * l_prod_971;
        double l_mult_1070 = -0.432e3 * l_prod_511;
        double l_mult_1071 = 0.3888e4 * l_prod_960;
        double l_mult_1072 = -0.7776e4 * l_prod_959;
        double l_sum_1073 = l_mult_1070 + (l_mult_1071 + l_mult_1072);
        double l_sum_1074 = l_mult_1064 + (l_mult_1065 + (l_mult_1066 + (l_mult_1034 + (l_mult_1067 + (l_mult_1068 + (l_mult_1069 + l_sum_1073))))));
        double l_mult_1075 = 0.2088e4 * l_prod_518;
        double l_mult_1076 = 0.15552e5 * l_prod_972;
        double l_mult_1077 = -0.10044e5 * l_prod_515;
        double l_mult_1078 = 0.31104e5 * l_prod_969;
        double l_mult_1079 = 0.15552e5 * l_prod_966;
        double l_mult_1080 = -0.23328e5 * l_prod_965;
        double l_mult_1081 = -0.7776e4 * l_prod_963;
        double l_mult_1082 = -0.10044e5 * l_prod_509;
        double l_mult_1083 = -0.23328e5 * l_prod_955;
        double l_mult_1084 = 0.31104e5 * l_prod_953;
        double l_mult_1085 = -0.46656e5 * l_prod_952;
        double l_mult_1086 = -0.23328e5 * l_prod_950;
        double l_mult_1087 = 0.15552e5 * l_prod_945;
        double l_mult_1088 = -0.23328e5 * l_prod_943;
        double l_mult_1089 = -0.7776e4 * l_prod_940;
        double l_sum_1090 = l_mult_1088 + l_mult_1089;
        double l_sum_1091 = l_mult_1075 + (l_mult_1035 + (l_mult_1076 + (l_mult_1069 + (l_mult_1077 + (l_mult_1078 + (l_mult_1039 + (l_mult_1079 + (l_mult_1080 + (l_mult_1081 + (l_mult_1082 + (l_mult_1044 + (l_mult_1083 + (l_mult_1084 + (l_mult_1085 + (l_mult_1086 + (l_mult_1087 + (l_mult_1049 + l_sum_1090)))))))))))))))));
        double l_mult_1092 = -0.1071e4 * l_prod_518;
        double l_mult_1093 = -0.1944e4 * l_prod_972;
        double l_mult_1094 = 0.9342e4 * l_prod_515;
        double l_mult_1095 = -0.21384e5 * l_prod_969;
        double l_mult_1096 = 0.23328e5 * l_prod_965;
        double l_mult_1097 = 0.11664e5 * l_prod_963;
        double l_mult_1098 = 0.2916e4 * l_prod_509;
        double l_mult_1099 = -0.21384e5 * l_prod_953;
        double l_mult_1100 = 0.23328e5 * l_prod_950;
        double l_sum_1101 = l_mult_1100 + l_sum_1021;
        double l_sum_1102 = l_mult_1092 + (l_mult_1055 + (l_mult_1093 + (l_mult_1094 + (l_mult_1095 + (l_mult_979 + (-0.19440e5 * l_prod_966 + (l_mult_1096 + (l_mult_1097 + (l_mult_1098 + (l_mult_1062 + (l_mult_1099 + (l_mult_984 + l_sum_1101))))))))))));
        double l_mult_1103 = 0.360e3 * l_prod_518;
        double l_mult_1104 = -0.3672e4 * l_prod_515;
        double l_mult_1105 = 0.3888e4 * l_prod_969;
        double l_mult_1106 = 0.10368e5 * l_prod_966;
        double l_mult_1107 = -0.432e3 * l_prod_509;
        double l_mult_1108 = 0.3888e4 * l_prod_953;
        double l_mult_1109 = -0.7776e4 * l_prod_950;
        double l_sum_1110 = l_mult_1107 + (l_mult_1108 + l_mult_1109);
        double l_sum_1111 = l_mult_1081 + l_sum_1110;
        double l_sum_1112 = l_mult_1103 + (l_mult_1067 + (l_mult_1104 + (l_mult_1105 + (l_mult_1106 + (l_mult_1040 + l_sum_1111)))));
        double l_mult_1113 = -0.6264e3 * l_prod_523;
        double l_mult_1114 = 0.4176e4 * l_prod_522;
        double l_mult_1115 = -0.3888e4 * l_prod_975;
        double l_mult_1116 = 0.4176e4 * l_prod_518;
        double l_mult_1117 = -0.20088e5 * l_prod_517;
        double l_mult_1118 = -0.15552e5 * l_prod_971;
        double l_mult_1119 = -0.15552e5 * l_prod_965;
        double l_mult_1120 = -0.3888e4 * l_prod_963;
        double l_mult_1121 = -0.30132e5 * l_prod_511;
        double l_mult_1122 = 0.46656e5 * l_prod_960;
        double l_mult_1123 = -0.30132e5 * l_prod_509;
        double l_mult_1124 = 0.93312e5 * l_prod_956;
        double l_mult_1125 = -0.69984e5 * l_prod_955;
        double l_mult_1126 = 0.46656e5 * l_prod_953;
        double l_mult_1127 = -0.69984e5 * l_prod_952;
        double l_mult_1128 = -0.46656e5 * l_prod_946;
        double l_mult_1129 = -0.93312e5 * l_prod_944;
        double l_mult_1130 = -0.46656e5 * l_prod_943;
        double l_mult_1131 = -0.38880e5 * l_prod_941;
        double l_mult_1132 = -0.38880e5 * l_prod_940;
        double l_mult_1133 = -0.11664e5 * l_prod_939;
        double l_mult_1134 = 0.1053e4 * l_prod_523;
        double l_mult_1135 = -0.5220e4 * l_prod_522;
        double l_mult_1136 = -0.7128e4 * l_prod_976;
        double l_mult_1137 = -0.5220e4 * l_prod_518;
        double l_mult_1138 = 0.18684e5 * l_prod_517;
        double l_mult_1139 = -0.7128e4 * l_prod_966;
        double l_mult_1140 = 0.47304e5 * l_prod_511;
        double l_mult_1141 = -0.58320e5 * l_prod_960;
        double l_mult_1142 = 0.47304e5 * l_prod_509;
        double l_mult_1143 = -0.116640e6 * l_prod_956;
        double l_mult_1144 = 0.69984e5 * l_prod_955;
        double l_mult_1145 = -0.58320e5 * l_prod_953;
        double l_mult_1146 = 0.69984e5 * l_prod_952;
        double l_mult_1147 = 0.69984e5 * l_prod_946;
        double l_mult_1148 = 0.139968e6 * l_prod_944;
        double l_mult_1149 = 0.69984e5 * l_prod_943;
        double l_mult_1150 = 0.77760e5 * l_prod_941;
        double l_mult_1151 = 0.77760e5 * l_prod_940;
        double l_mult_1152 = 0.29160e5 * l_prod_939;
        double l_mult_1153 = -0.1016e4 * l_prod_523;
        double l_mult_1154 = 0.3384e4 * l_prod_522;
        double l_mult_1155 = 0.1296e4 * l_prod_976;
        double l_mult_1156 = 0.3384e4 * l_prod_518;
        double l_mult_1157 = -0.7344e4 * l_prod_517;
        double l_mult_1158 = 0.1296e4 * l_prod_966;
        double l_mult_1159 = -0.36720e5 * l_prod_511;
        double l_mult_1160 = -0.36720e5 * l_prod_509;
        double l_mult_1161 = 0.62208e5 * l_prod_956;
        double l_mult_1162 = -0.77760e5 * l_prod_941;
        double l_mult_1163 = -0.77760e5 * l_prod_940;
        double l_mult_1164 = 0.594e3 * l_prod_523;
        double l_mult_1165 = -0.1197e4 * l_prod_522;
        double l_mult_1166 = -0.1197e4 * l_prod_518;
        double l_mult_1167 = 0.1188e4 * l_prod_517;
        double l_mult_1168 = 0.14256e5 * l_prod_511;
        double l_mult_1169 = 0.14256e5 * l_prod_509;
        double l_mult_1170 = -0.11664e5 * l_prod_956;
        double l_mult_1171 = 0.38880e5 * l_prod_941;
        double l_mult_1172 = 0.38880e5 * l_prod_940;
        double l_mult_1173 = -0.1944e3 * l_prod_523;
        double l_mult_1174 = 0.180e3 * l_prod_522;
        double l_mult_1175 = 0.180e3 * l_prod_518;
        double l_mult_1176 = -0.2268e4 * l_prod_511;
        double l_mult_1177 = -0.2268e4 * l_prod_509;
        double l_mult_1178 = 0.7776e4 * l_prod_947;
        double l_mult_1179 = 0.7776e4 * l_prod_945;
        double l_mult_1180 = 0.20736e5 * l_prod_942;
        double l_mult_1181 = 0.648e3 * l_prod_517;
        double l_mult_1182 = -0.3888e4 * l_prod_969;
        double l_mult_1183 = 0.432e3 * l_prod_517;
        double l_sum_1184 = l_mult_1183 + (l_mult_1182 + l_mult_980);
        double l_mult_1185 = -0.3888e4 * l_prod_972;
        double l_sum_1186 = l_mult_1062 + l_mult_983;
        double l_mult_1187 = 0.324e3 * l_prod_517;
        double l_sum_1188 = l_mult_1058 + l_mult_979;
        double l_sum_1189 = l_mult_1187 + (l_mult_1093 + l_sum_1188);
        double l_sum_1190 = l_mult_1183 + (l_mult_1185 + l_mult_978);
        double l_mult_1191 = 0.12852e5 * l_prod_517;
        double l_mult_1192 = -0.34992e5 * l_prod_972;
        double l_mult_1193 = -0.34992e5 * l_prod_969;
        double l_mult_1194 = 0.46656e5 * l_prod_968;
        double l_mult_1195 = -0.34992e5 * l_prod_956;
        double l_mult_1196 = 0.46656e5 * l_prod_955;
        double l_mult_1197 = 0.46656e5 * l_prod_952;
        double l_mult_1198 = -0.3240e4 * l_prod_517;
        double l_mult_1199 = 0.23328e5 * l_prod_969;
        double l_mult_1200 = 0.3888e4 * l_prod_956;
        double l_mult_1201 = 0.23328e5 * l_prod_972;
        double l_sum_1202 = l_mult_1200 + l_mult_1083;
        double l_mult_1203 = -0.6156e4 * l_prod_522;
        double l_mult_1204 = 0.25704e5 * l_prod_521;
        double l_mult_1205 = -0.34992e5 * l_prod_976;
        double l_mult_1206 = 0.15552e5 * l_prod_975;
        double l_mult_1207 = 0.25704e5 * l_prod_517;
        double l_mult_1208 = -0.69984e5 * l_prod_972;
        double l_mult_1209 = 0.46656e5 * l_prod_971;
        double l_mult_1210 = 0.15552e5 * l_prod_965;
        double l_mult_1211 = 0.38556e5 * l_prod_511;
        double l_mult_1212 = -0.104976e6 * l_prod_960;
        double l_mult_1213 = 0.69984e5 * l_prod_959;
        double l_mult_1214 = -0.104976e6 * l_prod_956;
        double l_mult_1215 = 0.139968e6 * l_prod_955;
        double l_mult_1216 = -0.69984e5 * l_prod_947;
        double l_mult_1217 = 0.93312e5 * l_prod_946;
        double l_mult_1218 = 0.93312e5 * l_prod_944;
        double l_mult_1219 = 0.6984e4 * l_prod_522;
        double l_mult_1220 = -0.22464e5 * l_prod_521;
        double l_mult_1221 = 0.23328e5 * l_prod_976;
        double l_mult_1222 = -0.22464e5 * l_prod_517;
        double l_mult_1223 = 0.46656e5 * l_prod_972;
        double l_mult_1224 = -0.57672e5 * l_prod_511;
        double l_mult_1225 = 0.128304e6 * l_prod_960;
        double l_mult_1226 = -0.69984e5 * l_prod_959;
        double l_mult_1227 = 0.128304e6 * l_prod_956;
        double l_mult_1228 = -0.139968e6 * l_prod_955;
        double l_mult_1229 = -0.139968e6 * l_prod_946;
        double l_mult_1230 = -0.139968e6 * l_prod_944;
        double l_mult_1231 = -0.4032e4 * l_prod_522;
        double l_mult_1232 = 0.7992e4 * l_prod_521;
        double l_mult_1233 = -0.3888e4 * l_prod_976;
        double l_mult_1234 = 0.7992e4 * l_prod_517;
        double l_mult_1235 = -0.7776e4 * l_prod_972;
        double l_mult_1236 = 0.42120e5 * l_prod_511;
        double l_mult_1237 = -0.66096e5 * l_prod_960;
        double l_mult_1238 = -0.66096e5 * l_prod_956;
        double l_mult_1239 = 0.1296e4 * l_prod_522;
        double l_mult_1240 = -0.1188e4 * l_prod_521;
        double l_mult_1241 = -0.1188e4 * l_prod_517;
        double l_mult_1242 = -0.15228e5 * l_prod_511;
        double l_mult_1243 = 0.11664e5 * l_prod_960;
        double l_mult_1244 = 0.11664e5 * l_prod_956;
        double l_mult_1245 = 0.46656e5 * l_prod_947;
        double l_mult_1246 = 0.2664e4 * l_prod_522;
        double l_mult_1247 = 0.42768e5 * l_prod_976;
        double l_mult_1248 = -0.23328e5 * l_prod_975;
        double l_mult_1249 = -0.6480e4 * l_prod_517;
        double l_mult_1250 = -0.46656e5 * l_prod_971;
        double l_mult_1251 = -0.9720e4 * l_prod_511;
        double l_mult_1252 = 0.69984e5 * l_prod_960;
        double l_sum_1253 = l_mult_1178 + l_mult_1128;
        double l_mult_1254 = -0.2214e4 * l_prod_522;
        double l_mult_1255 = 0.17496e5 * l_prod_521;
        double l_mult_1256 = -0.27216e5 * l_prod_976;
        double l_mult_1257 = 0.4212e4 * l_prod_517;
        double l_mult_1258 = -0.29160e5 * l_prod_972;
        double l_mult_1259 = 0.11664e5 * l_prod_511;
        double l_mult_1260 = -0.81648e5 * l_prod_960;
        double l_mult_1261 = -0.11664e5 * l_prod_947;
        double l_mult_1262 = 0.720e3 * l_prod_522;
        double l_mult_1263 = -0.4968e4 * l_prod_521;
        double l_mult_1264 = 0.3888e4 * l_prod_976;
        double l_mult_1265 = -0.648e3 * l_prod_517;
        double l_mult_1266 = -0.5832e4 * l_prod_511;
        double l_mult_1267 = 0.38880e5 * l_prod_960;
        double l_mult_1268 = -0.792e3 * l_prod_522;
        double l_mult_1269 = -0.22032e5 * l_prod_976;
        double l_mult_1270 = 0.864e3 * l_prod_517;
        double l_mult_1271 = 0.15552e5 * l_prod_971;
        double l_mult_1272 = 0.1296e4 * l_prod_511;
        double l_mult_1273 = -0.11664e5 * l_prod_960;
        double l_mult_1274 = 0.504e3 * l_prod_522;
        double l_mult_1275 = 0.12960e5 * l_prod_976;
        double l_mult_1276 = -0.1296e4 * l_prod_511;
        double l_mult_1277 = 0.108e3 * l_prod_522;
        double l_sum_1278 = l_mult_1277 + (l_mult_1240 + (l_mult_1264 + l_mult_1115));
        double l_mult_1279 = -0.6156e4 * l_prod_518;
        double l_mult_1280 = 0.25704e5 * l_prod_515;
        double l_mult_1281 = -0.69984e5 * l_prod_969;
        double l_mult_1282 = -0.34992e5 * l_prod_966;
        double l_mult_1283 = 0.46656e5 * l_prod_965;
        double l_mult_1284 = 0.15552e5 * l_prod_963;
        double l_mult_1285 = 0.38556e5 * l_prod_509;
        double l_mult_1286 = -0.104976e6 * l_prod_953;
        double l_mult_1287 = 0.139968e6 * l_prod_952;
        double l_mult_1288 = 0.69984e5 * l_prod_950;
        double l_mult_1289 = -0.69984e5 * l_prod_945;
        double l_mult_1290 = 0.93312e5 * l_prod_943;
        double l_mult_1291 = 0.6984e4 * l_prod_518;
        double l_mult_1292 = -0.22464e5 * l_prod_515;
        double l_mult_1293 = 0.46656e5 * l_prod_969;
        double l_mult_1294 = 0.23328e5 * l_prod_966;
        double l_mult_1295 = -0.57672e5 * l_prod_509;
        double l_mult_1296 = 0.128304e6 * l_prod_953;
        double l_mult_1297 = -0.139968e6 * l_prod_952;
        double l_mult_1298 = -0.69984e5 * l_prod_950;
        double l_mult_1299 = -0.139968e6 * l_prod_943;
        double l_mult_1300 = -0.4032e4 * l_prod_518;
        double l_mult_1301 = 0.7992e4 * l_prod_515;
        double l_mult_1302 = -0.7776e4 * l_prod_969;
        double l_mult_1303 = -0.3888e4 * l_prod_966;
        double l_mult_1304 = 0.42120e5 * l_prod_509;
        double l_mult_1305 = -0.66096e5 * l_prod_953;
        double l_mult_1306 = 0.1296e4 * l_prod_518;
        double l_mult_1307 = -0.1188e4 * l_prod_515;
        double l_mult_1308 = -0.15228e5 * l_prod_509;
        double l_mult_1309 = 0.11664e5 * l_prod_953;
        double l_mult_1310 = 0.46656e5 * l_prod_945;
        double l_mult_1311 = 0.2664e4 * l_prod_518;
        double l_mult_1312 = 0.42768e5 * l_prod_966;
        double l_mult_1313 = -0.46656e5 * l_prod_965;
        double l_mult_1314 = -0.23328e5 * l_prod_963;
        double l_mult_1315 = -0.9720e4 * l_prod_509;
        double l_mult_1316 = 0.69984e5 * l_prod_953;
        double l_sum_1317 = l_mult_1179 + l_mult_1130;
        double l_mult_1318 = -0.2214e4 * l_prod_518;
        double l_mult_1319 = 0.17496e5 * l_prod_515;
        double l_mult_1320 = -0.29160e5 * l_prod_969;
        double l_mult_1321 = -0.27216e5 * l_prod_966;
        double l_mult_1322 = 0.11664e5 * l_prod_509;
        double l_mult_1323 = -0.81648e5 * l_prod_953;
        double l_mult_1324 = -0.11664e5 * l_prod_945;
        double l_mult_1325 = 0.720e3 * l_prod_518;
        double l_mult_1326 = -0.4968e4 * l_prod_515;
        double l_mult_1327 = 0.3888e4 * l_prod_966;
        double l_mult_1328 = -0.5832e4 * l_prod_509;
        double l_mult_1329 = 0.38880e5 * l_prod_953;
        double l_mult_1330 = -0.792e3 * l_prod_518;
        double l_mult_1331 = -0.22032e5 * l_prod_966;
        double l_mult_1332 = 0.1296e4 * l_prod_509;
        double l_mult_1333 = -0.11664e5 * l_prod_953;
        double l_mult_1334 = 0.504e3 * l_prod_518;
        double l_mult_1335 = 0.12960e5 * l_prod_966;
        double l_mult_1336 = -0.1296e4 * l_prod_509;
        double l_mult_1337 = 0.108e3 * l_prod_518;
        double l_sum_1338 = l_mult_1337 + (l_mult_1307 + (l_mult_1327 + l_mult_1120));
        double l_mult_1339 = -0.31968e5 * l_prod_517;
        double l_mult_1340 = 0.77760e5 * l_prod_972;
        double l_mult_1341 = 0.77760e5 * l_prod_969;
        double l_mult_1342 = 0.116640e6 * l_prod_956;
        double l_mult_1343 = 0.26568e5 * l_prod_517;
        double l_mult_1344 = -0.50544e5 * l_prod_972;
        double l_mult_1345 = -0.50544e5 * l_prod_969;
        double l_mult_1346 = -0.139968e6 * l_prod_956;
        double l_mult_1347 = -0.8640e4 * l_prod_517;
        double l_mult_1348 = 0.7776e4 * l_prod_972;
        double l_mult_1349 = 0.7776e4 * l_prod_969;
        double l_mult_1350 = 0.69984e5 * l_prod_956;
        double l_mult_1351 = 0.7128e4 * l_prod_517;
        double l_mult_1352 = -0.4536e4 * l_prod_517;
        double l_mult_1353 = -0.864e3 * l_prod_517;
        double l_mult_1354 = 0.72e1 * l_prod_523;
        double l_mult_1355 = -0.180e3 * l_prod_512;
        double l_sum_1356 = l_mult_1354 + (l_mult_1355 + (0.1134e4 * l_prod_507 + (-0.2592e4 * l_prod_942 + l_mult_991)));
        double l_mult_1357 = 0.45e1 * l_prod_523;
        double l_mult_1358 = -0.99e2 * l_prod_512;
        double l_mult_1359 = 0.1188e4 * l_prod_509;
        double l_mult_1360 = 0.486e3 * l_prod_507;
        double l_mult_1361 = -0.5832e4 * l_prod_945;
        double l_mult_1362 = -0.648e3 * l_prod_942;
        double l_sum_1363 = l_mult_1362 + l_mult_990;
        double l_mult_1364 = 0.4e1 * l_prod_523;
        double l_mult_1365 = 0.216e3 * l_prod_515;
        double l_mult_1366 = -0.72e2 * l_prod_512;
        double l_mult_1367 = 0.216e3 * l_prod_507;
        double l_mult_1368 = -0.3888e4 * l_prod_945;
        double l_sum_1369 = l_mult_1367 + (l_mult_1368 + l_mult_988);
        double l_mult_1370 = 0.486e3 * l_prod_515;
        double l_mult_1371 = -0.648e3 * l_prod_966;
        double l_mult_1372 = -0.54e2 * l_prod_512;
        double l_sum_1373 = l_mult_1372 + (l_mult_1359 + (l_mult_1019 + l_mult_985));
        double l_sum_1374 = l_mult_1354 + (l_mult_1013 + (0.1134e4 * l_prod_515 + (-0.2592e4 * l_prod_966 + l_mult_981)));
        double l_mult_1375 = -0.3132e3 * l_prod_523;
        double l_mult_1376 = -0.5022e4 * l_prod_521;
        double l_mult_1377 = 0.5184e4 * l_prod_976;
        double l_mult_1378 = -0.1944e4 * l_prod_975;
        double l_mult_1379 = -0.34992e5 * l_prod_968;
        double l_mult_1380 = 0.20736e5 * l_prod_966;
        double l_mult_1381 = -0.31104e5 * l_prod_965;
        double l_mult_1382 = -0.9720e4 * l_prod_963;
        double l_mult_1383 = 0.2088e4 * l_prod_512;
        double l_mult_1384 = 0.15552e5 * l_prod_960;
        double l_mult_1385 = -0.20088e5 * l_prod_509;
        double l_mult_1386 = -0.31104e5 * l_prod_950;
        double l_mult_1387 = -0.5022e4 * l_prod_507;
        double l_mult_1388 = -0.11664e5 * l_prod_946;
        double l_mult_1389 = 0.31104e5 * l_prod_945;
        double l_mult_1390 = -0.46656e5 * l_prod_944;
        double l_mult_1391 = -0.34992e5 * l_prod_943;
        double l_mult_1392 = 0.5184e4 * l_prod_942;
        double l_mult_1393 = -0.15552e5 * l_prod_940;
        double l_mult_1394 = -0.1944e4 * l_prod_939;
        double l_sum_1395 = l_mult_1375 + (l_mult_1031 + (l_mult_1376 + (l_mult_1377 + (l_mult_1378 + (l_mult_1116 + (l_mult_1117 + (l_mult_1036 + (l_mult_1118 + (-0.15066e5 * l_prod_515 + (l_mult_1293 + (l_mult_1379 + (l_mult_1380 + (l_mult_1381 + (l_mult_1382 + (l_mult_1383 + (l_mult_1041 + (l_mult_1384 + (l_mult_1072 + (l_mult_1385 + (l_mult_1161 + (l_mult_1045 + (l_mult_1126 + (l_mult_1127 + (l_mult_1386 + (l_mult_1387 + (l_mult_1047 + (l_mult_1388 + (l_mult_1389 + (l_mult_1390 + (l_mult_1391 + (l_mult_1392 + (l_mult_1050 + (l_mult_1393 + l_mult_1394)))))))))))))))))))))))))))))))));
        double l_mult_1396 = 0.2565e3 * l_prod_523;
        double l_mult_1397 = 0.1458e4 * l_prod_521;
        double l_mult_1398 = -0.648e3 * l_prod_976;
        double l_mult_1399 = -0.58320e5 * l_prod_969;
        double l_mult_1400 = 0.34992e5 * l_prod_968;
        double l_mult_1401 = 0.19440e5 * l_prod_963;
        double l_mult_1402 = -0.1071e4 * l_prod_512;
        double l_mult_1403 = -0.1944e4 * l_prod_960;
        double l_mult_1404 = 0.18684e5 * l_prod_509;
        double l_mult_1405 = -0.42768e5 * l_prod_956;
        double l_mult_1406 = 0.46656e5 * l_prod_950;
        double l_mult_1407 = 0.1458e4 * l_prod_507;
        double l_mult_1408 = -0.21384e5 * l_prod_945;
        double l_mult_1409 = 0.34992e5 * l_prod_943;
        double l_sum_1410 = l_mult_1396 + (l_mult_1052 + (l_mult_1397 + (l_mult_1398 + (l_mult_1137 + (l_mult_1138 + (l_mult_1056 + (l_mult_978 + (0.23652e5 * l_prod_515 + (l_mult_1399 + (l_mult_1400 + (-0.37584e5 * l_prod_966 + (l_mult_1283 + (l_mult_1401 + (l_mult_1402 + (l_mult_1059 + (l_mult_1403 + (l_mult_1404 + (l_mult_1405 + (l_mult_983 + (l_mult_1145 + (l_mult_1146 + (l_mult_1406 + (l_mult_1407 + (l_mult_1002 + (l_mult_1408 + (l_mult_987 + (l_mult_1409 + l_sum_1363)))))))))))))))))))))))))));
        double l_mult_1411 = -0.148e3 * l_prod_523;
        double l_mult_1412 = -0.216e3 * l_prod_521;
        double l_mult_1413 = -0.11664e5 * l_prod_968;
        double l_mult_1414 = -0.19440e5 * l_prod_963;
        double l_mult_1415 = 0.360e3 * l_prod_512;
        double l_mult_1416 = -0.7344e4 * l_prod_509;
        double l_mult_1417 = 0.7776e4 * l_prod_956;
        double l_mult_1418 = -0.216e3 * l_prod_507;
        double l_mult_1419 = 0.3888e4 * l_prod_945;
        double l_mult_1420 = -0.11664e5 * l_prod_943;
        double l_sum_1421 = l_mult_1418 + (l_mult_1419 + l_mult_1420);
        double l_sum_1422 = l_mult_1411 + (l_mult_1064 + (l_mult_1412 + (l_mult_1156 + (l_mult_1157 + (l_mult_1068 + (-0.18360e5 * l_prod_515 + (l_mult_1078 + (l_mult_1413 + (0.33696e5 * l_prod_966 + (l_mult_1381 + (l_mult_1414 + (l_mult_1415 + (l_mult_1070 + (l_mult_1416 + (l_mult_1417 + (l_mult_1084 + (l_mult_1046 + (l_mult_1386 + l_sum_1421))))))))))))))))));
        double l_mult_1423 = 0.495e2 * l_prod_523;
        double l_mult_1424 = -0.5832e4 * l_prod_969;
        double l_mult_1425 = 0.9720e4 * l_prod_963;
        double l_sum_1426 = l_mult_1423 + (l_mult_1010 + (l_mult_1166 + (l_mult_1167 + (0.7128e4 * l_prod_515 + (l_mult_1424 + (-0.14904e5 * l_prod_966 + (l_mult_980 + (l_mult_1425 + l_sum_1373))))))));
        double l_mult_1427 = -0.72e1 * l_prod_523;
        double l_mult_1428 = 0.2592e4 * l_prod_966;
        double l_mult_1429 = -0.1944e4 * l_prod_963;
        double l_sum_1430 = l_mult_1427 + (l_mult_1175 + (-0.1134e4 * l_prod_515 + (l_mult_1428 + l_mult_1429)));
        double l_mult_1431 = -0.5022e4 * l_prod_515;
        double l_mult_1432 = 0.5184e4 * l_prod_966;
        double l_mult_1433 = 0.4176e4 * l_prod_512;
        double l_mult_1434 = -0.20088e5 * l_prod_511;
        double l_mult_1435 = -0.15552e5 * l_prod_959;
        double l_mult_1436 = -0.15552e5 * l_prod_950;
        double l_mult_1437 = -0.34992e5 * l_prod_946;
        double l_mult_1438 = -0.69984e5 * l_prod_944;
        double l_mult_1439 = -0.31104e5 * l_prod_941;
        double l_mult_1440 = -0.31104e5 * l_prod_940;
        double l_mult_1441 = -0.9720e4 * l_prod_939;
        double l_sum_1442 = l_mult_1375 + (l_mult_1031 + (l_mult_1376 + (l_mult_1377 + (l_mult_1378 + (l_mult_1075 + (l_mult_1035 + (l_mult_1076 + (l_mult_1069 + (l_mult_1431 + (l_mult_1038 + (l_mult_1413 + (l_mult_1432 + (l_mult_1040 + (l_mult_1429 + (l_mult_1433 + (l_mult_1434 + (l_mult_1042 + (l_mult_1435 + (l_mult_1385 + (l_mult_1161 + (l_mult_1045 + (l_mult_1084 + (l_mult_1085 + (l_mult_1436 + (-0.15066e5 * l_prod_507 + (l_mult_1245 + (l_mult_1437 + (l_mult_1310 + (l_mult_1438 + (l_mult_1391 + (l_mult_1180 + (l_mult_1439 + (l_mult_1440 + l_mult_1441)))))))))))))))))))))))))))))))));
        double l_mult_1443 = 0.1458e4 * l_prod_515;
        double l_mult_1444 = -0.5220e4 * l_prod_512;
        double l_mult_1445 = 0.18684e5 * l_prod_511;
        double l_mult_1446 = -0.58320e5 * l_prod_947;
        double l_mult_1447 = 0.34992e5 * l_prod_946;
        double l_mult_1448 = -0.58320e5 * l_prod_945;
        double l_mult_1449 = 0.69984e5 * l_prod_944;
        double l_mult_1450 = 0.46656e5 * l_prod_941;
        double l_mult_1451 = 0.46656e5 * l_prod_940;
        double l_mult_1452 = 0.19440e5 * l_prod_939;
        double l_sum_1453 = l_mult_1396 + (l_mult_1052 + (l_mult_1397 + (l_mult_1398 + (l_mult_1092 + (l_mult_1055 + (l_mult_1093 + (l_mult_1443 + (l_mult_1058 + (l_mult_1371 + (l_mult_1444 + (l_mult_1445 + (l_mult_1060 + (l_mult_982 + (l_mult_1404 + (l_mult_1405 + (l_mult_983 + (l_mult_1099 + (l_mult_984 + (l_mult_985 + (0.23652e5 * l_prod_507 + (l_mult_1446 + (l_mult_1447 + (l_mult_1448 + (l_mult_1449 + (l_mult_1409 + (-0.37584e5 * l_prod_942 + (l_mult_1450 + (l_mult_1451 + l_mult_1452))))))))))))))))))))))))))));
        double l_mult_1454 = -0.216e3 * l_prod_515;
        double l_mult_1455 = 0.3384e4 * l_prod_512;
        double l_mult_1456 = -0.7344e4 * l_prod_511;
        double l_mult_1457 = 0.31104e5 * l_prod_947;
        double l_mult_1458 = -0.19440e5 * l_prod_939;
        double l_sum_1459 = l_mult_1411 + (l_mult_1064 + (l_mult_1412 + (l_mult_1103 + (l_mult_1067 + (l_mult_1454 + (l_mult_1455 + (l_mult_1456 + (l_mult_1071 + (l_mult_1416 + (l_mult_1417 + (l_mult_1108 + (-0.18360e5 * l_prod_507 + (l_mult_1457 + (l_mult_1388 + (l_mult_1389 + (l_mult_1049 + (l_mult_1420 + (0.33696e5 * l_prod_942 + (l_mult_1439 + (l_mult_1440 + l_mult_1458))))))))))))))))))));
        double l_mult_1460 = -0.1197e4 * l_prod_512;
        double l_mult_1461 = 0.1188e4 * l_prod_511;
        double l_mult_1462 = -0.5832e4 * l_prod_947;
        double l_mult_1463 = 0.9720e4 * l_prod_939;
        double l_sum_1464 = l_mult_1423 + (l_mult_1010 + (l_mult_1028 + (l_mult_1460 + (l_mult_1461 + (l_mult_1359 + (0.7128e4 * l_prod_507 + (l_mult_1462 + (l_mult_1361 + (-0.14904e5 * l_prod_942 + (l_mult_989 + (l_mult_990 + l_mult_1463)))))))))));
        double l_mult_1465 = 0.180e3 * l_prod_512;
        double l_mult_1466 = 0.2592e4 * l_prod_942;
        double l_sum_1467 = l_mult_1427 + (l_mult_1465 + (-0.1134e4 * l_prod_507 + (l_mult_1466 + l_mult_1394)));
        double l_mult_1468 = -0.36e2 * l_prod_522;
        double l_mult_1469 = 0.648e3 * l_prod_511;
        double l_mult_1470 = -0.7776e4 * l_prod_956;
        double l_sum_1471 = l_mult_1007 + (l_mult_1470 + l_mult_984);
        double l_mult_1472 = 0.216e3 * l_prod_521;
        double l_mult_1473 = -0.27e2 * l_prod_522;
        double l_mult_1474 = 0.162e3 * l_prod_521;
        double l_mult_1475 = 0.324e3 * l_prod_511;
        double l_sum_1476 = l_mult_1475 + (l_mult_1403 + l_sum_1186);
        double l_mult_1477 = 0.324e3 * l_prod_521;
        double l_mult_1478 = -0.3078e4 * l_prod_522;
        double l_mult_1479 = 0.12852e5 * l_prod_521;
        double l_mult_1480 = -0.17496e5 * l_prod_976;
        double l_mult_1481 = 0.7776e4 * l_prod_975;
        double l_mult_1482 = -0.52488e5 * l_prod_969;
        double l_mult_1483 = 0.69984e5 * l_prod_968;
        double l_mult_1484 = 0.31104e5 * l_prod_965;
        double l_mult_1485 = 0.12852e5 * l_prod_511;
        double l_mult_1486 = -0.34992e5 * l_prod_960;
        double l_mult_1487 = -0.69984e5 * l_prod_956;
        double l_mult_1488 = 0.93312e5 * l_prod_955;
        double l_mult_1489 = 0.23328e5 * l_prod_946;
        double l_mult_1490 = 0.46656e5 * l_prod_944;
        double l_mult_1491 = 0.1332e4 * l_prod_522;
        double l_mult_1492 = -0.3240e4 * l_prod_521;
        double l_mult_1493 = 0.1944e4 * l_prod_976;
        double l_mult_1494 = 0.64152e5 * l_prod_969;
        double l_mult_1495 = -0.69984e5 * l_prod_968;
        double l_mult_1496 = -0.3240e4 * l_prod_511;
        double l_mult_1497 = 0.46656e5 * l_prod_956;
        double l_mult_1498 = 0.1944e4 * l_prod_947;
        double l_mult_1499 = -0.396e3 * l_prod_522;
        double l_mult_1500 = 0.432e3 * l_prod_521;
        double l_mult_1501 = -0.33048e5 * l_prod_969;
        double l_mult_1502 = 0.23328e5 * l_prod_968;
        double l_mult_1503 = 0.54e2 * l_prod_522;
        double l_mult_1504 = 0.5832e4 * l_prod_969;
        double l_mult_1505 = -0.11232e5 * l_prod_521;
        double l_mult_1506 = 0.21384e5 * l_prod_976;
        double l_mult_1507 = -0.11664e5 * l_prod_975;
        double l_mult_1508 = 0.23328e5 * l_prod_960;
        double l_sum_1509 = l_mult_1498 + l_mult_1388;
        double l_mult_1510 = -0.297e3 * l_prod_522;
        double l_mult_1511 = 0.2106e4 * l_prod_521;
        double l_mult_1512 = 0.1944e4 * l_prod_969;
        double l_mult_1513 = 0.3996e4 * l_prod_521;
        double l_mult_1514 = -0.11016e5 * l_prod_976;
        double l_mult_1515 = -0.324e3 * l_prod_521;
        double l_mult_1516 = 0.648e3 * l_prod_976;
        double l_sum_1517 = l_mult_1503 + (-0.594e3 * l_prod_521 + (l_mult_1493 + l_mult_1378));
        double l_mult_1518 = 0.25704e5 * l_prod_511;
        double l_mult_1519 = -0.69984e5 * l_prod_960;
        double l_mult_1520 = 0.46656e5 * l_prod_959;
        double l_mult_1521 = -0.52488e5 * l_prod_947;
        double l_mult_1522 = 0.31104e5 * l_prod_941;
        double l_mult_1523 = -0.22464e5 * l_prod_511;
        double l_mult_1524 = 0.64152e5 * l_prod_947;
        double l_mult_1525 = -0.69984e5 * l_prod_946;
        double l_mult_1526 = -0.46656e5 * l_prod_941;
        double l_mult_1527 = 0.7992e4 * l_prod_511;
        double l_mult_1528 = -0.7776e4 * l_prod_960;
        double l_mult_1529 = -0.33048e5 * l_prod_947;
        double l_mult_1530 = -0.1188e4 * l_prod_511;
        double l_mult_1531 = 0.5832e4 * l_prod_947;
        double l_mult_1532 = -0.6480e4 * l_prod_511;
        double l_mult_1533 = -0.46656e5 * l_prod_959;
        double l_sum_1534 = l_mult_1417 + (l_mult_1045 + (l_mult_1531 + l_mult_1437));
        double l_mult_1535 = 0.4212e4 * l_prod_511;
        double l_mult_1536 = -0.29160e5 * l_prod_960;
        double l_sum_1537 = l_mult_1462 + l_mult_1447;
        double l_mult_1538 = -0.648e3 * l_prod_511;
        double l_mult_1539 = 0.864e3 * l_prod_511;
        double l_mult_1540 = 0.15552e5 * l_prod_959;
        double l_sum_1541 = l_mult_1539 + (l_mult_1528 + l_mult_1540);
        double l_mult_1542 = 0.540e3 * l_prod_523;
        double l_mult_1543 = 0.19278e5 * l_prod_515;
        double l_mult_1544 = -0.23328e5 * l_prod_966;
        double l_mult_1545 = -0.6156e4 * l_prod_512;
        double l_mult_1546 = 0.62208e5 * l_prod_950;
        double l_mult_1547 = 0.19278e5 * l_prod_507;
        double l_mult_1548 = -0.104976e6 * l_prod_945;
        double l_mult_1549 = 0.104976e6 * l_prod_943;
        double l_mult_1550 = -0.23328e5 * l_prod_942;
        double l_mult_1551 = 0.62208e5 * l_prod_940;
        double l_mult_1552 = -0.360e3 * l_prod_523;
        double l_mult_1553 = -0.1620e4 * l_prod_521;
        double l_mult_1554 = -0.4860e4 * l_prod_515;
        double l_mult_1555 = 0.6984e4 * l_prod_512;
        double l_mult_1556 = -0.44928e5 * l_prod_509;
        double l_mult_1557 = -0.28836e5 * l_prod_507;
        double l_mult_1558 = 0.128304e6 * l_prod_945;
        double l_mult_1559 = -0.104976e6 * l_prod_943;
        double l_mult_1560 = 0.41472e5 * l_prod_942;
        double l_mult_1561 = 0.180e3 * l_prod_523;
        double l_mult_1562 = -0.4032e4 * l_prod_512;
        double l_mult_1563 = 0.15984e5 * l_prod_509;
        double l_mult_1564 = -0.15552e5 * l_prod_956;
        double l_mult_1565 = 0.21060e5 * l_prod_507;
        double l_mult_1566 = -0.66096e5 * l_prod_945;
        double l_mult_1567 = -0.36288e5 * l_prod_942;
        double l_mult_1568 = -0.54e2 * l_prod_523;
        double l_mult_1569 = 0.1296e4 * l_prod_512;
        double l_mult_1570 = -0.2376e4 * l_prod_509;
        double l_mult_1571 = -0.7614e4 * l_prod_507;
        double l_mult_1572 = 0.11664e5 * l_prod_945;
        double l_mult_1573 = 0.15552e5 * l_prod_942;
        double l_mult_1574 = -0.28836e5 * l_prod_515;
        double l_mult_1575 = 0.41472e5 * l_prod_966;
        double l_mult_1576 = 0.2664e4 * l_prod_512;
        double l_mult_1577 = -0.4860e4 * l_prod_507;
        double l_mult_1578 = 0.69984e5 * l_prod_945;
        double l_sum_1579 = l_mult_1466 + l_mult_1440;
        double l_mult_1580 = 0.135e3 * l_prod_523;
        double l_mult_1581 = 0.5832e4 * l_prod_515;
        double l_mult_1582 = -0.2214e4 * l_prod_512;
        double l_mult_1583 = -0.58320e5 * l_prod_956;
        double l_mult_1584 = 0.5832e4 * l_prod_507;
        double l_mult_1585 = -0.81648e5 * l_prod_945;
        double l_mult_1586 = -0.3888e4 * l_prod_942;
        double l_mult_1587 = -0.36e2 * l_prod_523;
        double l_mult_1588 = -0.648e3 * l_prod_515;
        double l_mult_1589 = 0.720e3 * l_prod_512;
        double l_mult_1590 = -0.9936e4 * l_prod_509;
        double l_mult_1591 = -0.2916e4 * l_prod_507;
        double l_mult_1592 = 0.38880e5 * l_prod_945;
        double l_mult_1593 = 0.21060e5 * l_prod_515;
        double l_mult_1594 = -0.36288e5 * l_prod_966;
        double l_mult_1595 = -0.792e3 * l_prod_512;
        double l_mult_1596 = 0.648e3 * l_prod_507;
        double l_mult_1597 = -0.2916e4 * l_prod_515;
        double l_mult_1598 = 0.504e3 * l_prod_512;
        double l_mult_1599 = -0.648e3 * l_prod_507;
        double l_mult_1600 = -0.7614e4 * l_prod_515;
        double l_mult_1601 = 0.108e3 * l_prod_512;
        double l_mult_1602 = -0.31968e5 * l_prod_511;
        double l_mult_1603 = 0.77760e5 * l_prod_960;
        double l_mult_1604 = 0.155520e6 * l_prod_956;
        double l_mult_1605 = -0.1620e4 * l_prod_522;
        double l_mult_1606 = 0.3564e4 * l_prod_521;
        double l_mult_1607 = 0.26568e5 * l_prod_511;
        double l_mult_1608 = -0.50544e5 * l_prod_960;
        double l_mult_1609 = -0.101088e6 * l_prod_956;
        double l_mult_1610 = 0.432e3 * l_prod_522;
        double l_mult_1611 = -0.432e3 * l_prod_521;
        double l_mult_1612 = -0.8640e4 * l_prod_511;
        double l_mult_1613 = 0.7776e4 * l_prod_960;
        double l_mult_1614 = 0.15552e5 * l_prod_956;
        double l_mult_1615 = 0.7128e4 * l_prod_511;
        double l_mult_1616 = 0.324e3 * l_prod_522;
        double l_mult_1617 = -0.4536e4 * l_prod_511;
        double l_mult_1618 = -0.864e3 * l_prod_511;
        double l_mult_1619 = -0.23328e5 * l_prod_976;
        double l_mult_1620 = -0.2268e4 * l_prod_521;
        double l_sum_1621 = l_mult_1618 + (l_mult_1613 + l_mult_1435);
        double l_sum_1622 = l_mult_1362 + l_mult_989;
        double l_mult_1623 = -0.3888e4 * l_prod_947;
        double l_sum_1624 = l_mult_1367 + (l_mult_1623 + l_mult_986);
        double l_mult_1625 = 0.486e3 * l_prod_521;
        double l_sum_1626 = l_mult_1372 + (l_mult_1461 + (l_mult_1001 + l_mult_982));
        double l_sum_1627 = l_mult_1354 + (l_mult_995 + (0.1134e4 * l_prod_521 + (-0.2592e4 * l_prod_976 + l_mult_977)));
        double l_mult_1628 = 0.20736e5 * l_prod_976;
        double l_mult_1629 = -0.9720e4 * l_prod_975;
        double l_mult_1630 = -0.31104e5 * l_prod_971;
        double l_mult_1631 = -0.31104e5 * l_prod_959;
        double l_mult_1632 = 0.15552e5 * l_prod_953;
        double l_mult_1633 = -0.15552e5 * l_prod_941;
        double l_sum_1634 = l_mult_1375 + (l_mult_1114 + (-0.15066e5 * l_prod_521 + (l_mult_1628 + (l_mult_1629 + (l_mult_1075 + (l_mult_1117 + (l_mult_1223 + (l_mult_1630 + (l_mult_1431 + (l_mult_1078 + (l_mult_1379 + (l_mult_1432 + (l_mult_1119 + (l_mult_1429 + (l_mult_1383 + (l_mult_1434 + (l_mult_1122 + (l_mult_1631 + (l_mult_1082 + (l_mult_1161 + (l_mult_1125 + (l_mult_1632 + (l_mult_1085 + (l_mult_1109 + (l_mult_1387 + (l_mult_1457 + (l_mult_1437 + (l_mult_1087 + (l_mult_1390 + (l_mult_1420 + (l_mult_1392 + (l_mult_1633 + (l_mult_1089 + l_mult_1394)))))))))))))))))))))))))))))))));
        double l_mult_1635 = 0.19440e5 * l_prod_975;
        double l_mult_1636 = -0.58320e5 * l_prod_972;
        double l_mult_1637 = -0.1944e4 * l_prod_953;
        double l_mult_1638 = -0.21384e5 * l_prod_947;
        double l_sum_1639 = l_mult_1396 + (l_mult_1135 + (0.23652e5 * l_prod_521 + (-0.37584e5 * l_prod_976 + (l_mult_1635 + (l_mult_1092 + (l_mult_1138 + (l_mult_1636 + (l_mult_1209 + (l_mult_1443 + (l_mult_1095 + (l_mult_1400 + (l_mult_1371 + (l_mult_980 + (l_mult_1402 + (l_mult_1445 + (l_mult_1141 + (l_mult_1520 + (l_mult_1098 + (l_mult_1405 + (l_mult_1144 + (l_mult_1637 + (l_mult_984 + (l_mult_1407 + (l_mult_1638 + (l_mult_1447 + (l_mult_1020 + (l_mult_987 + l_sum_1622)))))))))))))))))))))))))));
        double l_mult_1640 = -0.19440e5 * l_prod_975;
        double l_mult_1641 = 0.3888e4 * l_prod_947;
        double l_sum_1642 = l_mult_1418 + (l_mult_1641 + l_mult_1388);
        double l_sum_1643 = l_mult_1411 + (l_mult_1154 + (-0.18360e5 * l_prod_521 + (0.33696e5 * l_prod_976 + (l_mult_1640 + (l_mult_1103 + (l_mult_1157 + (l_mult_1036 + (l_mult_1630 + (l_mult_1454 + (l_mult_1105 + (l_mult_1413 + (l_mult_1415 + (l_mult_1456 + (l_mult_1042 + (l_mult_1631 + (l_mult_1107 + (l_mult_1417 + (l_mult_1083 + l_sum_1642))))))))))))))))));
        double l_mult_1644 = 0.9720e4 * l_prod_975;
        double l_mult_1645 = -0.5832e4 * l_prod_972;
        double l_sum_1646 = l_mult_1423 + (l_mult_1165 + (0.7128e4 * l_prod_521 + (-0.14904e5 * l_prod_976 + (l_mult_1644 + (l_mult_1028 + (l_mult_1167 + (l_mult_1645 + (l_mult_978 + l_sum_1626))))))));
        double l_mult_1647 = 0.2592e4 * l_prod_976;
        double l_sum_1648 = l_mult_1427 + (l_mult_1174 + (-0.1134e4 * l_prod_521 + (l_mult_1647 + l_mult_1378)));
        double l_mult_1649 = -0.36e2 * l_prod_518;
        double l_mult_1650 = 0.648e3 * l_prod_509;
        double l_mult_1651 = 0.324e3 * l_prod_515;
        double l_sum_1652 = l_mult_1020 + l_mult_987;
        double l_mult_1653 = -0.27e2 * l_prod_518;
        double l_mult_1654 = 0.162e3 * l_prod_515;
        double l_mult_1655 = 0.324e3 * l_prod_509;
        double l_sum_1656 = l_mult_1637 + l_mult_984;
        double l_sum_1657 = l_mult_1655 + (l_mult_1062 + l_sum_1656);
        double l_sum_1658 = l_mult_1371 + l_mult_980;
        double l_sum_1659 = l_mult_1025 + (l_mult_1470 + l_mult_983);
        double l_sum_1660 = l_mult_1365 + (l_mult_1182 + l_mult_979);
        double l_sum_1661 = l_mult_1028 + (l_mult_1167 + (l_mult_1645 + l_mult_978));
        double l_mult_1662 = -0.3078e4 * l_prod_518;
        double l_mult_1663 = -0.52488e5 * l_prod_972;
        double l_mult_1664 = 0.31104e5 * l_prod_971;
        double l_mult_1665 = 0.12852e5 * l_prod_515;
        double l_mult_1666 = -0.17496e5 * l_prod_966;
        double l_mult_1667 = 0.7776e4 * l_prod_963;
        double l_mult_1668 = 0.12852e5 * l_prod_509;
        double l_mult_1669 = -0.34992e5 * l_prod_953;
        double l_mult_1670 = 0.93312e5 * l_prod_952;
        double l_mult_1671 = 0.23328e5 * l_prod_943;
        double l_mult_1672 = 0.1332e4 * l_prod_518;
        double l_mult_1673 = 0.5832e4 * l_prod_972;
        double l_mult_1674 = -0.11232e5 * l_prod_515;
        double l_mult_1675 = 0.21384e5 * l_prod_966;
        double l_mult_1676 = -0.11664e5 * l_prod_963;
        double l_mult_1677 = -0.3240e4 * l_prod_509;
        double l_mult_1678 = 0.23328e5 * l_prod_953;
        double l_mult_1679 = 0.1944e4 * l_prod_945;
        double l_sum_1680 = l_mult_1679 + l_mult_1420;
        double l_mult_1681 = -0.396e3 * l_prod_518;
        double l_mult_1682 = 0.3996e4 * l_prod_515;
        double l_mult_1683 = -0.11016e5 * l_prod_966;
        double l_mult_1684 = 0.54e2 * l_prod_518;
        double l_mult_1685 = 0.1944e4 * l_prod_966;
        double l_sum_1686 = l_mult_1684 + (-0.594e3 * l_prod_515 + (l_mult_1685 + l_mult_1429));
        double l_mult_1687 = 0.64152e5 * l_prod_972;
        double l_mult_1688 = -0.3240e4 * l_prod_515;
        double l_mult_1689 = -0.297e3 * l_prod_518;
        double l_mult_1690 = 0.2106e4 * l_prod_515;
        double l_mult_1691 = -0.324e3 * l_prod_515;
        double l_mult_1692 = 0.648e3 * l_prod_966;
        double l_mult_1693 = -0.33048e5 * l_prod_972;
        double l_mult_1694 = 0.432e3 * l_prod_515;
        double l_mult_1695 = 0.1944e4 * l_prod_972;
        double l_mult_1696 = 0.19278e5 * l_prod_521;
        double l_mult_1697 = 0.62208e5 * l_prod_959;
        double l_mult_1698 = 0.25704e5 * l_prod_509;
        double l_mult_1699 = 0.15552e5 * l_prod_950;
        double l_mult_1700 = -0.104976e6 * l_prod_947;
        double l_mult_1701 = 0.104976e6 * l_prod_946;
        double l_mult_1702 = -0.52488e5 * l_prod_945;
        double l_mult_1703 = 0.62208e5 * l_prod_941;
        double l_mult_1704 = 0.31104e5 * l_prod_940;
        double l_mult_1705 = -0.4860e4 * l_prod_521;
        double l_mult_1706 = -0.1620e4 * l_prod_515;
        double l_mult_1707 = -0.44928e5 * l_prod_511;
        double l_mult_1708 = -0.22464e5 * l_prod_509;
        double l_mult_1709 = 0.128304e6 * l_prod_947;
        double l_mult_1710 = -0.104976e6 * l_prod_946;
        double l_mult_1711 = 0.64152e5 * l_prod_945;
        double l_mult_1712 = -0.46656e5 * l_prod_940;
        double l_mult_1713 = 0.15984e5 * l_prod_511;
        double l_mult_1714 = 0.7992e4 * l_prod_509;
        double l_mult_1715 = -0.66096e5 * l_prod_947;
        double l_mult_1716 = -0.33048e5 * l_prod_945;
        double l_mult_1717 = -0.2376e4 * l_prod_511;
        double l_mult_1718 = -0.1188e4 * l_prod_509;
        double l_mult_1719 = 0.11664e5 * l_prod_947;
        double l_mult_1720 = 0.5832e4 * l_prod_945;
        double l_mult_1721 = -0.28836e5 * l_prod_521;
        double l_mult_1722 = 0.41472e5 * l_prod_976;
        double l_mult_1723 = -0.6480e4 * l_prod_509;
        double l_mult_1724 = 0.69984e5 * l_prod_947;
        double l_sum_1725 = l_mult_1466 + l_mult_1439;
        double l_mult_1726 = 0.5832e4 * l_prod_521;
        double l_mult_1727 = 0.4212e4 * l_prod_509;
        double l_mult_1728 = -0.81648e5 * l_prod_947;
        double l_mult_1729 = -0.648e3 * l_prod_521;
        double l_mult_1730 = -0.9936e4 * l_prod_511;
        double l_mult_1731 = -0.648e3 * l_prod_509;
        double l_mult_1732 = 0.38880e5 * l_prod_947;
        double l_mult_1733 = 0.21060e5 * l_prod_521;
        double l_mult_1734 = -0.36288e5 * l_prod_976;
        double l_mult_1735 = 0.864e3 * l_prod_509;
        double l_mult_1736 = -0.2916e4 * l_prod_521;
        double l_mult_1737 = -0.7614e4 * l_prod_521;
        double l_mult_1738 = -0.69984e5 * l_prod_953;
        double l_mult_1739 = -0.69984e5 * l_prod_943;
        double l_mult_1740 = -0.7776e4 * l_prod_953;
        double l_mult_1741 = -0.46656e5 * l_prod_950;
        double l_sum_1742 = l_mult_1720 + l_mult_1391;
        double l_mult_1743 = -0.29160e5 * l_prod_953;
        double l_sum_1744 = l_mult_1361 + l_mult_1409;
        double l_sum_1745 = l_mult_1735 + (l_mult_1740 + l_mult_1699);
        double l_mult_1746 = -0.31968e5 * l_prod_509;
        double l_mult_1747 = 0.77760e5 * l_prod_953;
        double l_mult_1748 = -0.1620e4 * l_prod_518;
        double l_mult_1749 = 0.3564e4 * l_prod_515;
        double l_mult_1750 = 0.26568e5 * l_prod_509;
        double l_mult_1751 = -0.50544e5 * l_prod_953;
        double l_mult_1752 = 0.432e3 * l_prod_518;
        double l_mult_1753 = -0.432e3 * l_prod_515;
        double l_mult_1754 = -0.8640e4 * l_prod_509;
        double l_mult_1755 = 0.7776e4 * l_prod_953;
        double l_mult_1756 = 0.7128e4 * l_prod_509;
        double l_mult_1757 = 0.324e3 * l_prod_518;
        double l_mult_1758 = -0.2268e4 * l_prod_515;
        double l_mult_1759 = -0.4536e4 * l_prod_509;
        double l_mult_1760 = -0.864e3 * l_prod_509;
        double l_sum_1761 = l_mult_1760 + (l_mult_1755 + l_mult_1436);
        double l_mult_1762 = 0.2268e4 * l_prod_517;
        double l_mult_1763 = 0.972e3 * l_prod_517;
        double l_mult_1764 = 0.594e3 * l_prod_507;
        double l_mult_1765 = -0.1944e4 * l_prod_942;
        double l_sum_1766 = l_mult_1372 + (l_mult_1764 + (l_mult_1765 + l_mult_991));
        double l_mult_1767 = -0.1296e4 * l_prod_942;
        double l_sum_1768 = l_mult_1764 + (l_mult_1361 + l_mult_988);
        double l_mult_1769 = -0.30132e5 * l_prod_517;
        double l_mult_1770 = -0.46656e5 * l_prod_968;
        double l_mult_1771 = -0.38880e5 * l_prod_965;
        double l_mult_1772 = -0.93312e5 * l_prod_952;
        double l_mult_1773 = -0.38880e5 * l_prod_950;
        double l_mult_1774 = -0.10044e5 * l_prod_507;
        double l_mult_1775 = 0.10368e5 * l_prod_942;
        double l_mult_1776 = -0.23328e5 * l_prod_940;
        double l_mult_1777 = -0.3888e4 * l_prod_939;
        double l_mult_1778 = 0.47304e5 * l_prod_517;
        double l_mult_1779 = 0.77760e5 * l_prod_965;
        double l_mult_1780 = 0.29160e5 * l_prod_963;
        double l_mult_1781 = 0.77760e5 * l_prod_950;
        double l_mult_1782 = 0.9342e4 * l_prod_507;
        double l_mult_1783 = -0.7128e4 * l_prod_942;
        double l_mult_1784 = 0.23328e5 * l_prod_940;
        double l_mult_1785 = -0.36720e5 * l_prod_517;
        double l_mult_1786 = -0.77760e5 * l_prod_965;
        double l_mult_1787 = -0.77760e5 * l_prod_950;
        double l_mult_1788 = -0.3672e4 * l_prod_507;
        double l_mult_1789 = 0.1296e4 * l_prod_942;
        double l_mult_1790 = 0.14256e5 * l_prod_517;
        double l_mult_1791 = 0.38880e5 * l_prod_965;
        double l_mult_1792 = 0.38880e5 * l_prod_950;
        double l_mult_1793 = -0.2268e4 * l_prod_517;
        double l_mult_1794 = -0.23328e5 * l_prod_941;
        double l_mult_1795 = -0.7776e4 * l_prod_939;
        double l_sum_1796 = l_mult_1776 + l_mult_1795;
        double l_sum_1797 = l_mult_1383 + (l_mult_1041 + (l_mult_1384 + (l_mult_1072 + (l_mult_1082 + (l_mult_1044 + (l_mult_1083 + (l_mult_1632 + (l_mult_1046 + (l_mult_1109 + (l_mult_1774 + (l_mult_1457 + (l_mult_1048 + (l_mult_1389 + (l_mult_1390 + (l_mult_1088 + (l_mult_1573 + (l_mult_1794 + l_sum_1796)))))))))))))))));
        double l_mult_1798 = 0.23328e5 * l_prod_941;
        double l_mult_1799 = 0.11664e5 * l_prod_939;
        double l_sum_1800 = l_mult_1784 + l_mult_1799;
        double l_sum_1801 = l_mult_1402 + (l_mult_1059 + (l_mult_1403 + (l_mult_1098 + (l_mult_1062 + (l_mult_1637 + (l_mult_1782 + (l_mult_1638 + (l_mult_986 + (l_mult_1408 + (l_mult_987 + (l_mult_988 + (-0.19440e5 * l_prod_942 + (l_mult_1798 + l_sum_1800)))))))))))));
        double l_sum_1802 = l_mult_1089 + l_mult_1795;
        double l_sum_1803 = l_mult_1415 + (l_mult_1070 + (l_mult_1107 + (l_mult_1788 + (l_mult_1641 + (l_mult_1419 + (l_mult_1775 + (l_mult_1050 + l_sum_1802)))))));
        double l_sum_1804 = l_mult_1007 + (l_mult_1623 + l_mult_989);
        double l_sum_1805 = l_mult_1623 + l_mult_987;
        double l_sum_1806 = l_mult_1475 + (l_mult_1403 + l_sum_1003);
        double l_mult_1807 = 0.38556e5 * l_prod_517;
        double l_mult_1808 = -0.104976e6 * l_prod_972;
        double l_mult_1809 = 0.69984e5 * l_prod_971;
        double l_mult_1810 = 0.93312e5 * l_prod_968;
        double l_mult_1811 = -0.34992e5 * l_prod_947;
        double l_mult_1812 = 0.46656e5 * l_prod_946;
        double l_mult_1813 = 0.15552e5 * l_prod_941;
        double l_mult_1814 = -0.57672e5 * l_prod_517;
        double l_mult_1815 = 0.128304e6 * l_prod_972;
        double l_mult_1816 = -0.69984e5 * l_prod_971;
        double l_mult_1817 = -0.139968e6 * l_prod_968;
        double l_mult_1818 = 0.23328e5 * l_prod_947;
        double l_mult_1819 = 0.42120e5 * l_prod_517;
        double l_mult_1820 = -0.66096e5 * l_prod_972;
        double l_mult_1821 = -0.15228e5 * l_prod_517;
        double l_mult_1822 = 0.11664e5 * l_prod_972;
        double l_mult_1823 = -0.9720e4 * l_prod_517;
        double l_mult_1824 = 0.69984e5 * l_prod_972;
        double l_sum_1825 = l_mult_1641 + l_mult_1048;
        double l_sum_1826 = l_mult_1244 + (l_mult_1125 + l_sum_1825);
        double l_mult_1827 = 0.11664e5 * l_prod_517;
        double l_mult_1828 = -0.81648e5 * l_prod_972;
        double l_mult_1829 = -0.11664e5 * l_prod_969;
        double l_mult_1830 = -0.5832e4 * l_prod_517;
        double l_mult_1831 = 0.38880e5 * l_prod_972;
        double l_mult_1832 = 0.1296e4 * l_prod_517;
        double l_mult_1833 = -0.11664e5 * l_prod_972;
        double l_mult_1834 = -0.1296e4 * l_prod_517;
        double l_mult_1835 = 0.25704e5 * l_prod_507;
        double l_mult_1836 = -0.34992e5 * l_prod_942;
        double l_mult_1837 = 0.69984e5 * l_prod_940;
        double l_mult_1838 = 0.15552e5 * l_prod_939;
        double l_mult_1839 = -0.22464e5 * l_prod_507;
        double l_mult_1840 = 0.42768e5 * l_prod_942;
        double l_mult_1841 = -0.69984e5 * l_prod_940;
        double l_mult_1842 = -0.23328e5 * l_prod_939;
        double l_mult_1843 = 0.7992e4 * l_prod_507;
        double l_mult_1844 = -0.22032e5 * l_prod_942;
        double l_mult_1845 = -0.1188e4 * l_prod_507;
        double l_mult_1846 = 0.3888e4 * l_prod_942;
        double l_sum_1847 = l_mult_1601 + (l_mult_1845 + (l_mult_1846 + l_mult_1777));
        double l_mult_1848 = 0.23328e5 * l_prod_942;
        double l_mult_1849 = 0.17496e5 * l_prod_507;
        double l_mult_1850 = -0.29160e5 * l_prod_947;
        double l_mult_1851 = -0.27216e5 * l_prod_942;
        double l_mult_1852 = -0.4968e4 * l_prod_507;
        double l_mult_1853 = 0.12960e5 * l_prod_942;
        double l_mult_1854 = 0.77760e5 * l_prod_947;
        double l_mult_1855 = -0.50544e5 * l_prod_947;
        double l_mult_1856 = -0.36e2 * l_prod_512;
        double l_mult_1857 = 0.324e3 * l_prod_507;
        double l_mult_1858 = -0.27e2 * l_prod_512;
        double l_sum_1859 = 0.162e3 * l_prod_507 + (l_mult_1002 + l_sum_1652);
        double l_mult_1860 = 0.62208e5 * l_prod_971;
        double l_mult_1861 = -0.104976e6 * l_prod_969;
        double l_mult_1862 = 0.104976e6 * l_prod_968;
        double l_mult_1863 = 0.62208e5 * l_prod_965;
        double l_mult_1864 = -0.3078e4 * l_prod_512;
        double l_mult_1865 = -0.52488e5 * l_prod_960;
        double l_mult_1866 = 0.31104e5 * l_prod_959;
        double l_mult_1867 = -0.52488e5 * l_prod_953;
        double l_mult_1868 = 0.31104e5 * l_prod_950;
        double l_mult_1869 = -0.34992e5 * l_prod_945;
        double l_mult_1870 = 0.15552e5 * l_prod_940;
        double l_mult_1871 = -0.44928e5 * l_prod_517;
        double l_mult_1872 = 0.128304e6 * l_prod_969;
        double l_mult_1873 = -0.104976e6 * l_prod_968;
        double l_mult_1874 = 0.1332e4 * l_prod_512;
        double l_mult_1875 = 0.5832e4 * l_prod_960;
        double l_mult_1876 = 0.64152e5 * l_prod_953;
        double l_mult_1877 = -0.1620e4 * l_prod_507;
        double l_mult_1878 = 0.23328e5 * l_prod_945;
        double l_mult_1879 = 0.648e3 * l_prod_942;
        double l_sum_1880 = l_mult_1879 + l_mult_1089;
        double l_mult_1881 = 0.15984e5 * l_prod_517;
        double l_mult_1882 = -0.66096e5 * l_prod_969;
        double l_mult_1883 = -0.396e3 * l_prod_512;
        double l_mult_1884 = -0.33048e5 * l_prod_953;
        double l_mult_1885 = -0.2376e4 * l_prod_517;
        double l_mult_1886 = 0.11664e5 * l_prod_969;
        double l_mult_1887 = 0.54e2 * l_prod_512;
        double l_mult_1888 = 0.5832e4 * l_prod_953;
        double l_sum_1889 = l_mult_1887 + (l_mult_1718 + (l_mult_1888 + l_mult_1109));
        double l_mult_1890 = 0.69984e5 * l_prod_969;
        double l_mult_1891 = 0.64152e5 * l_prod_960;
        double l_sum_1892 = l_mult_1879 + l_mult_1050;
        double l_mult_1893 = -0.81648e5 * l_prod_969;
        double l_mult_1894 = -0.297e3 * l_prod_512;
        double l_mult_1895 = -0.9936e4 * l_prod_517;
        double l_mult_1896 = 0.38880e5 * l_prod_969;
        double l_mult_1897 = 0.1944e4 * l_prod_953;
        double l_mult_1898 = -0.33048e5 * l_prod_960;
        double l_mult_1899 = 0.1944e4 * l_prod_960;
        double l_sum_1900 = l_mult_1887 + (l_mult_1530 + (l_mult_1875 + l_mult_1072));
        double l_mult_1901 = 0.12852e5 * l_prod_507;
        double l_mult_1902 = -0.17496e5 * l_prod_942;
        double l_mult_1903 = 0.7776e4 * l_prod_939;
        double l_mult_1904 = -0.11232e5 * l_prod_507;
        double l_mult_1905 = 0.21384e5 * l_prod_942;
        double l_mult_1906 = 0.3996e4 * l_prod_507;
        double l_mult_1907 = -0.11016e5 * l_prod_942;
        double l_mult_1908 = 0.1944e4 * l_prod_942;
        double l_sum_1909 = l_mult_1887 + (-0.594e3 * l_prod_507 + (l_mult_1908 + l_mult_1394));
        double l_mult_1910 = -0.3240e4 * l_prod_507;
        double l_sum_1911 = l_mult_1419 + (l_mult_1390 + (l_mult_1908 + l_mult_1794));
        double l_mult_1912 = 0.2106e4 * l_prod_507;
        double l_sum_1913 = l_mult_1765 + l_mult_1798;
        double l_mult_1914 = -0.324e3 * l_prod_507;
        double l_mult_1915 = 0.432e3 * l_prod_507;
        double l_sum_1916 = l_mult_1908 + l_mult_1776;
        double l_mult_1917 = -0.29160e5 * l_prod_945;
        double l_sum_1918 = l_mult_1765 + l_mult_1784;
        double l_mult_1919 = 0.77760e5 * l_prod_945;
        double l_mult_1920 = -0.1620e4 * l_prod_512;
        double l_mult_1921 = -0.50544e5 * l_prod_945;
        double l_mult_1922 = 0.432e3 * l_prod_512;
        double l_mult_1923 = 0.3564e4 * l_prod_507;
        double l_mult_1924 = 0.324e3 * l_prod_512;
        double l_mult_1925 = -0.2268e4 * l_prod_507;
        double l_mult_1926 = -0.432e3 * l_prod_507;
        double l_sum_1927 = l_mult_1764 + (l_mult_1462 + l_mult_986);
        double l_mult_1928 = -0.38880e5 * l_prod_971;
        double l_mult_1929 = -0.38880e5 * l_prod_959;
        double l_mult_1930 = -0.93312e5 * l_prod_955;
        double l_mult_1931 = 0.29160e5 * l_prod_975;
        double l_mult_1932 = 0.77760e5 * l_prod_971;
        double l_mult_1933 = 0.77760e5 * l_prod_959;
        double l_mult_1934 = -0.77760e5 * l_prod_971;
        double l_mult_1935 = -0.77760e5 * l_prod_959;
        double l_mult_1936 = 0.38880e5 * l_prod_971;
        double l_mult_1937 = 0.38880e5 * l_prod_959;
        double l_sum_1938 = l_mult_1025 + (l_mult_1368 + l_mult_990);
        double l_sum_1939 = l_mult_1655 + (l_mult_1637 + l_sum_1021);
        double l_sum_1940 = l_mult_1368 + l_mult_987;
        double l_mult_1941 = 0.69984e5 * l_prod_965;
        double l_mult_1942 = 0.46656e5 * l_prod_943;
        double l_mult_1943 = -0.69984e5 * l_prod_965;
        double l_sum_1944 = l_mult_1419 + l_mult_1088;
        double l_mult_1945 = 0.69984e5 * l_prod_941;
        double l_mult_1946 = -0.69984e5 * l_prod_941;
        double l_sum_1947 = l_mult_1942 + l_mult_1784;
        double l_sum_1948 = l_mult_1088 + l_mult_1776;
        double l_sum_1949 = l_mult_1086 + l_sum_1944;
        double l_vdot_1950 = vdot4(v_917, vcons4(0.e0, 0.e0, 0.e0, 0.e0)) + (vdot4(v_918,
            vcons4(0.e0, l_mult_995 + (l_mult_996 + (l_mult_997 + l_mult_989)),
                l_mult_998 + (l_mult_999 + (l_mult_1000 + (l_mult_1001 + l_sum_1003))),
                l_mult_1004 + (l_mult_1005 + (l_mult_1006 + l_sum_1009)))) + (vdot4(v_919,
            vcons4(l_sum_1012, 0.e0, l_mult_1013 + (l_mult_1014 + (l_mult_1015 + l_mult_990)),
                l_mult_1016 + (l_mult_1017 + (l_mult_1018 + (l_mult_1019 + l_sum_1021))))) + (vdot4(v_920,
            vcons4(l_mult_1022 + (l_mult_1023 + (l_mult_1024 + l_sum_1027)), l_sum_1030, 0.e0, l_sum_1051)) + (vdot4(
            v_921, vcons4(l_sum_1063, l_sum_1074, l_sum_1012, 0.e0)) + (vdot4(v_922,
            vcons4(l_sum_1091, l_sum_1102, l_sum_1112, l_sum_1030)) + (vdot4(v_923,
            vcons4(0.e0,
                l_mult_1113 + (l_mult_1114 + (l_mult_1032 + (l_mult_1066 + (l_mult_1115 + (l_mult_1116 + (l_mult_1117 + (l_mult_1036 + (l_mult_1118 + (l_mult_1077 + (l_mult_1078 + (l_mult_1039 + (l_mult_1106 + (l_mult_1119 + (l_mult_1120 + (0.6264e4 * l_prod_512 + (l_mult_1121 + (l_mult_1122 + (l_mult_1043 + (l_mult_1123 + (l_mult_1124 + (l_mult_1125 + (l_mult_1126 + (l_mult_1127 + (l_mult_1086 + (-0.20088e5 * l_prod_507 + (0.62208e5 * l_prod_947 + (l_mult_1128 + (0.62208e5 * l_prod_945 + (l_mult_1129 + (l_mult_1130 + (0.25920e5 * l_prod_942 + (l_mult_1131 + (l_mult_1132 + l_mult_1133))))))))))))))))))))))))))))))))),
                l_mult_1134 + (l_mult_1135 + (l_mult_1053 + (l_mult_1136 + (l_mult_977 + (l_mult_1137 + (l_mult_1138 + (l_mult_1056 + (l_mult_978 + (l_mult_1094 + (l_mult_1095 + (l_mult_979 + (l_mult_1139 + (l_mult_980 + (l_mult_981 + (-0.12447e5 * l_prod_512 + (l_mult_1140 + (l_mult_1141 + (l_mult_1061 + (l_mult_1142 + (l_mult_1143 + (l_mult_1144 + (l_mult_1145 + (l_mult_1146 + (l_mult_1100 + (0.44388e5 * l_prod_507 + (-0.112752e6 * l_prod_947 + (l_mult_1147 + (-0.112752e6 * l_prod_945 + (l_mult_1148 + (l_mult_1149 + (-0.61560e5 * l_prod_942 + (l_mult_1150 + (l_mult_1151 + l_mult_1152))))))))))))))))))))))))))))))))),
                l_mult_1153 + (l_mult_1154 + (l_mult_1065 + (l_mult_1155 + (l_mult_1156 + (l_mult_1157 + (l_mult_1068 + (l_mult_1104 + (l_mult_1105 + (l_mult_1158 + (0.13392e5 * l_prod_512 + (l_mult_1159 + (l_mult_1042 + (l_mult_1072 + (l_mult_1160 + (l_mult_1161 + (l_mult_1083 + (l_mult_1084 + (l_mult_1046 + (l_mult_1109 + (-0.52272e5 * l_prod_507 + (0.101088e6 * l_prod_947 + (l_mult_1128 + (0.101088e6 * l_prod_945 + (l_mult_1129 + (l_mult_1130 + (0.77760e5 * l_prod_942 + (l_mult_1162 + (l_mult_1163 + -0.38880e5 * l_prod_939)))))))))))))))))))))))))))))) + (vdot4(
            v_924,
            vcons4(
                l_mult_1164 + (l_mult_1165 + (l_mult_999 + (l_mult_1166 + (l_mult_1167 + (l_mult_1017 + (-0.8289e4 * l_prod_512 + (l_mult_1168 + (l_mult_1001 + (l_mult_1169 + (l_mult_1170 + (l_mult_1019 + (0.34668e5 * l_prod_507 + (-0.44712e5 * l_prod_947 + (l_mult_986 + (-0.44712e5 * l_prod_945 + (l_mult_987 + (l_mult_988 + (-0.55080e5 * l_prod_942 + (l_mult_1171 + (l_mult_1172 + l_mult_1152)))))))))))))))))))),
                l_mult_1173 + (l_mult_1174 + (l_mult_1175 + (0.2808e4 * l_prod_512 + (l_mult_1176 + (l_mult_1177 + (-0.12312e5 * l_prod_507 + (l_mult_1178 + (l_mult_1179 + (l_mult_1180 + (l_mult_1050 + (l_mult_1089 + l_mult_1133))))))))))),
                l_mult_1167 + (l_mult_1170 + l_mult_987), l_mult_1181 + (l_mult_1182 + (l_mult_1062 + l_mult_984)))) + (vdot4(
            v_925, vcons4(l_sum_1184, 0.e0, l_mult_1181 + (l_mult_1185 + l_sum_1186), l_sum_1189)) + (vdot4(v_926,
            vcons4(0.e0, l_sum_1190, 0.e0, 0.e0)) + (vdot4(v_927,
            vcons4(
                l_mult_1191 + (l_mult_1192 + (l_mult_1057 + (l_mult_1193 + (l_mult_1194 + (l_mult_1096 + (l_mult_1195 + (l_mult_1196 + (l_mult_1197 + l_mult_987)))))))),
                l_mult_1198 + (l_mult_1068 + (l_mult_1199 + (l_mult_1039 + (l_mult_1080 + (l_mult_1200 + l_mult_1046))))),
                l_sum_1184, 0.e0)) + (vdot4(v_928,
            vcons4(l_mult_1198 + (l_mult_1201 + (l_mult_1037 + (l_mult_1105 + (l_mult_1039 + l_sum_1202)))),
                l_sum_1189, 0.e0, l_sum_1190)) + (vdot4(v_929,
            vcons4(0.e0, 0.e0,
                l_mult_1203 + (l_mult_1204 + (l_mult_1205 + (l_mult_1206 + (l_mult_1207 + (l_mult_1208 + (l_mult_1209 + (l_mult_1193 + (l_mult_1194 + (l_mult_1210 + (l_mult_1211 + (l_mult_1212 + (l_mult_1213 + (l_mult_1214 + (l_mult_1215 + (l_mult_1146 + (l_mult_1216 + (l_mult_1217 + (l_mult_1218 + l_mult_1171)))))))))))))))))),
                l_mult_1219 + (l_mult_1220 + (l_mult_1221 + (l_mult_1034 + (l_mult_1222 + (l_mult_1223 + (l_mult_1037 + (l_mult_1199 + (l_mult_1039 + (l_mult_1040 + (l_mult_1224 + (l_mult_1225 + (l_mult_1226 + (l_mult_1227 + (l_mult_1228 + (l_mult_1127 + (0.124416e6 * l_prod_947 + (l_mult_1229 + (l_mult_1230 + l_mult_1162)))))))))))))))))))) + (vdot4(
            v_930,
            vcons4(
                l_mult_1231 + (l_mult_1232 + (l_mult_1233 + (l_mult_1234 + (l_mult_1235 + (l_mult_1182 + (l_mult_1236 + (l_mult_1237 + (l_mult_1061 + (l_mult_1238 + (l_mult_1196 + (l_mult_984 + (-0.108864e6 * l_prod_947 + (l_mult_1217 + (l_mult_1218 + l_mult_1150)))))))))))))),
                l_mult_1239 + (l_mult_1240 + (l_mult_1241 + (l_mult_1242 + (l_mult_1243 + (l_mult_1244 + (l_mult_1245 + (l_mult_1048 + (l_mult_1049 + l_mult_1131)))))))),
                l_mult_1246 + (l_mult_1220 + (l_mult_1247 + (l_mult_1248 + (l_mult_1249 + (l_mult_1223 + (l_mult_1250 + (l_mult_1105 + (l_mult_1039 + (l_mult_1251 + (l_mult_1252 + (l_mult_1226 + (l_mult_1244 + (l_mult_1125 + l_sum_1253))))))))))))),
                l_mult_1254 + (l_mult_1255 + (l_mult_1256 + (l_mult_1054 + (l_mult_1257 + (l_mult_1258 + (l_mult_1057 + (l_mult_1058 + (l_mult_979 + (l_mult_1259 + (l_mult_1260 + (l_mult_1213 + (l_mult_1170 + (l_mult_1144 + (l_mult_1261 + l_mult_1147)))))))))))))))) + (vdot4(
            v_931,
            vcons4(
                l_mult_1262 + (l_mult_1263 + (l_mult_1264 + (l_mult_1265 + (l_mult_1068 + (l_mult_1266 + (l_mult_1267 + (l_mult_1043 + (l_mult_1200 + (l_mult_1083 + l_sum_1253))))))))),
                l_mult_1268 + (l_mult_1232 + (l_mult_1269 + (l_mult_1206 + (l_mult_1270 + (l_mult_1235 + (l_mult_1271 + (l_mult_1272 + (l_mult_1273 + l_mult_1061)))))))),
                l_mult_1274 + (l_mult_1263 + (l_mult_1275 + (l_mult_1034 + (l_mult_1067 + (l_mult_1068 + (l_mult_1069 + (l_mult_1276 + (l_mult_1243 + l_mult_1043)))))))),
                l_sum_1278)) + (vdot4(v_932,
            vcons4(
                l_mult_1279 + (l_mult_1207 + (l_mult_1192 + (l_mult_1271 + (l_mult_1280 + (l_mult_1281 + (l_mult_1194 + (l_mult_1282 + (l_mult_1283 + (l_mult_1284 + (l_mult_1285 + (l_mult_1214 + (l_mult_1144 + (l_mult_1286 + (l_mult_1287 + (l_mult_1288 + (l_mult_1289 + (l_mult_1218 + (l_mult_1290 + l_mult_1172)))))))))))))))))),
                l_mult_1291 + (l_mult_1222 + (l_mult_1201 + (l_mult_1069 + (l_mult_1292 + (l_mult_1293 + (l_mult_1039 + (l_mult_1294 + (l_mult_1080 + (l_mult_1081 + (l_mult_1295 + (l_mult_1227 + (l_mult_1125 + (l_mult_1296 + (l_mult_1297 + (l_mult_1298 + (0.124416e6 * l_prod_945 + (l_mult_1230 + (l_mult_1299 + l_mult_1163)))))))))))))))))),
                l_mult_1300 + (l_mult_1234 + (l_mult_1185 + (l_mult_1301 + (l_mult_1302 + (l_mult_1303 + (l_mult_1304 + (l_mult_1238 + (l_mult_983 + (l_mult_1305 + (l_mult_1197 + (l_mult_1100 + (-0.108864e6 * l_prod_945 + (l_mult_1218 + (l_mult_1290 + l_mult_1151)))))))))))))),
                l_mult_1306 + (l_mult_1241 + (l_mult_1307 + (l_mult_1308 + (l_mult_1244 + (l_mult_1309 + (l_mult_1310 + (l_mult_1049 + (l_mult_1088 + l_mult_1132)))))))))) + (vdot4(
            v_933,
            vcons4(
                l_mult_1311 + (l_mult_1249 + (l_mult_1068 + (l_mult_1292 + (l_mult_1293 + (l_mult_1039 + (l_mult_1312 + (l_mult_1313 + (l_mult_1314 + (l_mult_1315 + (l_mult_1244 + (l_mult_1316 + (l_mult_1127 + (l_mult_1298 + l_sum_1317))))))))))))),
                l_mult_1318 + (l_mult_1257 + (l_mult_1093 + (l_mult_1319 + (l_mult_1320 + (l_mult_979 + (l_mult_1321 + (l_mult_1096 + (l_mult_1097 + (l_mult_1322 + (l_mult_1170 + (l_mult_1323 + (l_mult_1146 + (l_mult_1288 + (l_mult_1324 + l_mult_1149)))))))))))))),
                l_mult_1325 + (l_mult_1265 + (l_mult_1326 + (l_mult_1105 + (l_mult_1327 + (l_mult_1328 + (l_mult_1200 + (l_mult_1329 + (l_mult_1046 + (l_mult_1086 + l_sum_1317))))))))),
                l_mult_1330 + (l_mult_1270 + (l_mult_1301 + (l_mult_1302 + (l_mult_1331 + (l_mult_1210 + (l_mult_1284 + (l_mult_1332 + (l_mult_1333 + l_mult_1100)))))))))) + (vdot4(
            v_934,
            vcons4(
                l_mult_1334 + (l_mult_1067 + (l_mult_1326 + (l_mult_1105 + (l_mult_1335 + (l_mult_1040 + (l_mult_1081 + (l_mult_1336 + (l_mult_1309 + l_mult_1086)))))))),
                l_sum_1338,
                l_mult_1339 + (l_mult_1340 + (l_mult_1250 + (l_mult_1341 + (-0.93312e5 * l_prod_968 + (l_mult_1313 + (l_mult_1342 + (l_mult_1228 + (l_mult_1297 + l_mult_1129)))))))),
                l_mult_1343 + (l_mult_1344 + (l_mult_1057 + (l_mult_1345 + (l_mult_1194 + (l_mult_1096 + (l_mult_1346 + (l_mult_1215 + (l_mult_1287 + l_mult_1148)))))))))) + (vdot4(
            v_935,
            vcons4(
                l_mult_1347 + (l_mult_1348 + (l_mult_1349 + (l_mult_1350 + (l_mult_1045 + (l_mult_1085 + l_mult_1129))))),
                l_mult_1351 + (l_mult_1235 + (l_mult_1345 + (l_mult_1194 + (l_mult_1283 + (l_mult_1170 + l_mult_1146))))),
                l_mult_1352 + (l_mult_1068 + (l_mult_1078 + (l_mult_1039 + (l_mult_1080 + (l_mult_1244 + l_mult_1127))))),
                l_mult_1353 + (l_mult_1349 + l_mult_1119))) + (vdot4(v_936,
            vcons4(
                l_mult_1351 + (l_mult_1344 + (l_mult_1209 + (l_mult_1302 + (l_mult_1194 + (l_mult_1170 + l_mult_1144))))),
                l_mult_1352 + (l_mult_1036 + (l_mult_1037 + (l_mult_1105 + (l_mult_1039 + (l_mult_1244 + l_mult_1125))))),
                l_mult_1265 + (l_mult_1068 + (l_mult_1105 + l_mult_1039)), l_mult_1353 + (l_mult_1348 + l_mult_1118))) + vdot4(
            v_916,
            vcons4(l_sum_993,
                l_mult_994 + (-0.405e3 * l_prod_512 + (0.1836e4 * l_prod_507 + (-0.3240e4 * l_prod_942 + l_mult_991))),
                0.e0, 0.e0)))))))))))))))))))));
        double l_vdot_1951 = vdot4(v_917, vcons4(0.e0, 0.e0, 0.e0, 0.e0)) + (vdot4(v_918,
            vcons4(0.e0, 0.e0, 0.e0, 0.e0)) + (vdot4(v_919,
            vcons4(0.e0, 0.e0, l_sum_1356,
                l_mult_1357 + (l_mult_1028 + (l_mult_1358 + (l_mult_1359 + (l_mult_1360 + (l_mult_1361 + l_sum_1363))))))) + (vdot4(
            v_920,
            vcons4(
                l_mult_1364 + (l_mult_1022 + (l_mult_1365 + (l_mult_1366 + (l_mult_1332 + (l_mult_1026 + l_sum_1369))))),
                l_mult_1357 + (l_mult_1016 + (l_mult_1370 + (l_mult_1371 + l_sum_1373))), l_sum_1374, l_sum_1051)) + (vdot4(
            v_921, vcons4(l_sum_1063, l_sum_1074, l_sum_1012, 0.e0)) + (vdot4(v_922,
            vcons4(l_sum_1395, l_sum_1410, l_sum_1422, l_sum_1426)) + (vdot4(v_923,
            vcons4(l_sum_1430, l_sum_1442, l_sum_1453, l_sum_1459)) + (vdot4(v_924,
            vcons4(l_sum_1464, l_sum_1467, l_mult_1010 + (l_mult_1461 + (l_mult_1462 + l_mult_989)),
                l_mult_1468 + (l_mult_1183 + (l_mult_1469 + (l_mult_1470 + (l_mult_1002 + l_mult_987)))))) + (vdot4(
            v_925,
            vcons4(l_mult_1468 + (l_mult_1181 + (l_mult_1058 + l_sum_1471)),
                l_mult_1010 + (l_mult_1167 + (l_mult_1424 + l_mult_980)),
                l_mult_1468 + (l_mult_1472 + (l_mult_1469 + (l_mult_1008 + l_sum_1003))),
                l_mult_1473 + (l_mult_1474 + (l_mult_1187 + (l_mult_1093 + l_sum_1476))))) + (vdot4(v_926,
            vcons4(l_mult_1468 + (l_mult_1472 + (l_mult_1181 + (l_mult_1185 + l_sum_1188))),
                l_mult_1468 + (l_mult_1477 + (l_mult_1398 + l_sum_1009)),
                l_mult_1468 + (l_mult_1477 + (l_mult_1398 + l_sum_1190)), l_sum_1012)) + (vdot4(v_927,
            vcons4(
                l_mult_1478 + (l_mult_1479 + (l_mult_1480 + (l_mult_1481 + (l_mult_1207 + (l_mult_1208 + (l_mult_1209 + (l_mult_1482 + (l_mult_1483 + (l_mult_1484 + (l_mult_1485 + (l_mult_1486 + (l_mult_1061 + (l_mult_1487 + (l_mult_1488 + (l_mult_1146 + (-0.17496e5 * l_prod_947 + (l_mult_1489 + (l_mult_1490 + l_mult_989)))))))))))))))))),
                l_mult_1491 + (l_mult_1492 + (l_mult_1493 + (l_mult_1222 + (l_mult_1223 + (l_mult_1037 + (l_mult_1494 + (l_mult_1495 + (l_mult_1313 + (l_mult_1496 + (l_mult_1071 + (l_mult_1497 + (l_mult_1045 + (l_mult_1127 + (l_mult_1498 + l_mult_1049)))))))))))))),
                l_mult_1499 + (l_mult_1500 + (l_mult_1234 + (l_mult_1235 + (l_mult_1501 + (l_mult_1502 + (l_mult_1484 + l_sum_1471)))))),
                l_mult_1503 + (l_mult_1241 + (l_mult_1504 + l_mult_1040)))) + (vdot4(v_928,
            vcons4(
                l_mult_1491 + (l_mult_1505 + (l_mult_1506 + (l_mult_1507 + (l_mult_1249 + (l_mult_1223 + (l_mult_1250 + (l_mult_1504 + (l_mult_1379 + (l_mult_1496 + (l_mult_1508 + (l_mult_1043 + (l_mult_1417 + (l_mult_1045 + l_sum_1509))))))))))))),
                l_mult_1510 + (l_mult_1511 + (l_mult_1011 + (l_mult_1257 + (l_mult_1258 + (l_mult_1057 + (l_mult_1424 + (l_mult_1400 + l_sum_1476))))))),
                l_mult_590 + (l_mult_1412 + (l_mult_1265 + (l_mult_1068 + (l_mult_1512 + l_mult_1413)))),
                l_mult_1499 + (l_mult_1513 + (l_mult_1514 + (l_mult_1481 + (l_mult_1270 + (l_mult_1235 + (l_mult_1271 + l_sum_1009)))))))) + (vdot4(
            v_929,
            vcons4(l_mult_590 + (l_mult_1515 + (l_mult_1516 + (l_mult_1067 + (l_mult_1068 + l_mult_1069)))),
                l_sum_1517,
                l_mult_1478 + (l_mult_1479 + (l_mult_1480 + (l_mult_1481 + (l_mult_1191 + (l_mult_1192 + (l_mult_1057 + (-0.17496e5 * l_prod_969 + (l_mult_1502 + (l_mult_980 + (l_mult_1518 + (l_mult_1519 + (l_mult_1520 + (l_mult_1487 + (l_mult_1488 + (l_mult_1197 + (l_mult_1521 + (l_mult_1147 + (l_mult_1449 + l_mult_1522)))))))))))))))))),
                l_mult_1491 + (l_mult_1492 + (l_mult_1493 + (l_mult_1198 + (l_mult_1068 + (l_mult_1512 + (l_mult_1523 + (l_mult_1122 + (l_mult_1043 + (l_mult_1497 + (l_mult_1045 + (l_mult_1046 + (l_mult_1524 + (l_mult_1525 + (l_mult_1438 + l_mult_1526)))))))))))))))) + (vdot4(
            v_930,
            vcons4(
                l_mult_1499 + (l_mult_1500 + (l_mult_1183 + (l_mult_1527 + (l_mult_1528 + (l_mult_1470 + (l_mult_1529 + (l_mult_1489 + (l_mult_987 + l_mult_1522)))))))),
                l_mult_1503 + (l_mult_1530 + (l_mult_1531 + l_mult_1050)),
                l_mult_1491 + (l_mult_1505 + (l_mult_1506 + (l_mult_1507 + (l_mult_1198 + (l_mult_1201 + (l_mult_1037 + (l_mult_1512 + (l_mult_1413 + (l_mult_1532 + (l_mult_1122 + (l_mult_1533 + l_sum_1534))))))))))),
                l_mult_1510 + (l_mult_1511 + (l_mult_1011 + (l_mult_1187 + (l_mult_1093 + (l_mult_1535 + (l_mult_1536 + (l_mult_1061 + (l_mult_1062 + (l_mult_983 + l_sum_1537))))))))))) + (vdot4(
            v_931,
            vcons4(l_mult_590 + (l_mult_1412 + (l_mult_1538 + (l_mult_1071 + l_sum_1509))),
                l_mult_1499 + (l_mult_1513 + (l_mult_1514 + (l_mult_1481 + (l_mult_1183 + (l_mult_1185 + (l_mult_978 + l_sum_1541)))))),
                l_mult_590 + (l_mult_1515 + (l_mult_1516 + l_sum_1073)), l_sum_1517)) + (vdot4(v_932,
            vcons4(
                l_mult_1542 + (l_mult_1478 + (0.6426e4 * l_prod_521 + (-0.5832e4 * l_prod_976 + (l_mult_977 + (l_mult_1279 + (l_mult_1207 + (l_mult_1192 + (l_mult_1271 + (l_mult_1543 + (l_mult_1482 + (l_mult_1400 + (l_mult_1544 + (l_mult_1484 + (l_mult_1425 + (l_mult_1545 + (l_mult_1518 + (l_mult_1486 + (l_mult_1540 + (0.51408e5 * l_prod_509 + (l_mult_1346 + (l_mult_1488 + (l_mult_1286 + (l_mult_1287 + (l_mult_1546 + (l_mult_1547 + (l_mult_1521 + (l_mult_1447 + (l_mult_1548 + (l_mult_1148 + (l_mult_1549 + (l_mult_1550 + (l_mult_1522 + (l_mult_1551 + l_mult_1463))))))))))))))))))))))))))))))))),
                l_mult_1552 + (l_mult_1491 + (l_mult_1553 + (l_mult_1516 + (l_mult_1311 + (l_mult_1249 + (l_mult_1068 + (l_mult_1554 + (l_mult_1504 + (l_mult_1428 + (l_mult_1555 + (l_mult_1523 + (l_mult_1508 + (l_mult_1072 + (l_mult_1556 + (l_mult_1124 + (l_mult_1045 + (l_mult_1316 + (l_mult_1127 + (l_mult_1386 + (l_mult_1557 + (l_mult_1524 + (l_mult_1437 + (l_mult_1558 + (l_mult_1230 + (l_mult_1559 + (l_mult_1560 + (l_mult_1526 + (-0.93312e5 * l_prod_940 + l_mult_1458)))))))))))))))))))))))))))),
                l_mult_1561 + (l_mult_1499 + (l_mult_1472 + (l_mult_1330 + (l_mult_1270 + (l_mult_1023 + (l_mult_1562 + (l_mult_1527 + (l_mult_1008 + (l_mult_1563 + (l_mult_1564 + (l_mult_1333 + (l_mult_1565 + (l_mult_1529 + (l_mult_986 + (l_mult_1566 + (l_mult_1490 + (l_mult_1409 + (l_mult_1567 + (l_mult_1522 + (l_mult_1551 + l_mult_1452)))))))))))))))))))),
                l_mult_1568 + (l_mult_1503 + (l_mult_1337 + (l_mult_1569 + (l_mult_1530 + (l_mult_1570 + (l_mult_1571 + (l_mult_1531 + (l_mult_1572 + (l_mult_1573 + (l_mult_1050 + (l_mult_1393 + l_mult_1441))))))))))))) + (vdot4(
            v_933,
            vcons4(
                l_mult_1552 + (l_mult_1491 + (l_mult_1553 + (l_mult_1516 + (l_mult_1291 + (l_mult_1222 + (l_mult_1201 + (l_mult_1069 + (l_mult_1574 + (l_mult_1494 + (l_mult_1379 + (l_mult_1575 + (l_mult_1313 + (l_mult_1414 + (l_mult_1576 + (l_mult_1532 + (l_mult_1071 + (l_mult_1556 + (l_mult_1124 + (l_mult_1045 + (l_mult_1296 + (l_mult_1297 + (-0.93312e5 * l_prod_950 + (l_mult_1577 + (l_mult_1531 + (l_mult_1578 + (l_mult_1438 + (l_mult_1559 + l_sum_1579))))))))))))))))))))))))))),
                l_mult_1580 + (l_mult_1510 + (l_mult_1474 + (l_mult_1318 + (l_mult_1257 + (l_mult_1093 + (l_mult_1581 + (l_mult_1424 + (l_mult_1303 + (l_mult_1582 + (l_mult_1535 + (l_mult_1403 + (0.34992e5 * l_prod_509 + (l_mult_1583 + (l_mult_983 + (l_mult_1323 + (l_mult_1146 + (l_mult_1406 + (l_mult_1584 + (l_mult_1462 + (l_mult_1585 + (l_mult_1449 + (l_mult_1549 + (l_mult_1586 + l_mult_1451))))))))))))))))))))))),
                l_mult_1587 + (l_mult_590 + (l_mult_1334 + (l_mult_1067 + (l_mult_1588 + (l_mult_1589 + (l_mult_1538 + (l_mult_1590 + (l_mult_1417 + (l_mult_1309 + (l_mult_1591 + (l_mult_1498 + (l_mult_1592 + (l_mult_1049 + (l_mult_1391 + l_sum_1579)))))))))))))),
                l_mult_1561 + (l_mult_1499 + (l_mult_1472 + (l_mult_1300 + (l_mult_1234 + (l_mult_1185 + (l_mult_1593 + (l_mult_1501 + (l_mult_979 + (l_mult_1594 + (l_mult_1484 + (l_mult_1401 + (l_mult_1595 + (l_mult_1539 + (l_mult_1563 + (l_mult_1564 + (l_mult_1305 + (l_mult_1197 + (l_mult_1546 + (l_mult_1596 + (l_mult_1324 + l_mult_1409)))))))))))))))))))))) + (vdot4(
            v_934,
            vcons4(
                l_mult_1587 + (l_mult_590 + (l_mult_1325 + (l_mult_1265 + (l_mult_1597 + (l_mult_1512 + (l_mult_1428 + (l_mult_1598 + (l_mult_1070 + (l_mult_1590 + (l_mult_1417 + (l_mult_1329 + (l_mult_1046 + (l_mult_1386 + (l_mult_1599 + (l_mult_1572 + l_mult_1391))))))))))))))),
                l_mult_1568 + (l_mult_1503 + (l_mult_1306 + (l_mult_1241 + (l_mult_1600 + (l_mult_1504 + (l_mult_1079 + (l_mult_1040 + (l_mult_1382 + (l_mult_1601 + (l_mult_1570 + (l_mult_1309 + l_mult_1436))))))))))),
                0.4320e4 * l_prod_522 + (-0.15984e5 * l_prod_521 + (0.19440e5 * l_prod_976 + (l_mult_1034 + (l_mult_1339 + (l_mult_1340 + (l_mult_1250 + (0.58320e5 * l_prod_969 + (l_mult_1495 + (l_mult_1381 + (l_mult_1602 + (l_mult_1603 + (l_mult_1533 + (l_mult_1604 + (-0.186624e6 * l_prod_955 + (l_mult_1297 + (0.58320e5 * l_prod_947 + (l_mult_1525 + (l_mult_1230 + l_mult_1439)))))))))))))))))),
                l_mult_1605 + (l_mult_1606 + (l_mult_1011 + (l_mult_1351 + (l_mult_1235 + (l_mult_1424 + (l_mult_1607 + (l_mult_1608 + (l_mult_1061 + (l_mult_1609 + (l_mult_1488 + (l_mult_1146 + (l_mult_1216 + (l_mult_1147 + (l_mult_1148 + l_mult_1450)))))))))))))))) + (vdot4(
            v_935,
            vcons4(
                l_mult_1610 + (l_mult_1611 + (l_mult_1353 + (l_mult_1612 + (l_mult_1613 + (l_mult_1614 + (0.34992e5 * l_prod_947 + (l_mult_1048 + (l_mult_1390 + l_mult_1439)))))))),
                l_mult_1605 + (l_mult_1606 + (l_mult_1011 + (l_mult_1343 + (l_mult_1344 + (l_mult_1057 + (l_mult_1281 + (l_mult_1483 + (l_mult_1283 + (l_mult_1615 + (l_mult_1528 + (l_mult_1609 + (l_mult_1488 + (l_mult_1287 + (l_mult_1462 + l_mult_1449)))))))))))))),
                l_mult_1616 + (l_mult_1515 + (l_mult_1352 + (l_mult_1068 + (l_mult_1504 + (l_mult_1617 + (l_mult_1071 + (l_mult_1161 + (l_mult_1045 + (l_mult_1127 + (l_mult_1531 + l_mult_1438)))))))))),
                l_mult_1610 + (l_mult_1611 + (l_mult_1347 + (l_mult_1348 + (0.34992e5 * l_prod_969 + (l_mult_1039 + (l_mult_1381 + (l_mult_1618 + (l_mult_1614 + l_mult_1085)))))))))) + (vdot4(
            v_936,
            vcons4(
                l_mult_1605 + (0.13284e5 * l_prod_521 + (l_mult_1619 + (l_mult_1054 + (l_mult_1351 + (l_mult_1344 + (l_mult_1209 + (l_mult_1424 + (l_mult_1400 + (l_mult_1615 + (l_mult_1608 + (l_mult_1520 + (l_mult_1564 + (l_mult_1488 + l_sum_1537))))))))))))),
                l_mult_1616 + (l_mult_1620 + (l_mult_1493 + (l_mult_1265 + (l_mult_1068 + (l_mult_1617 + (l_mult_1042 + (l_mult_1043 + l_sum_1534))))))),
                l_mult_1616 + (l_mult_1620 + (l_mult_1493 + (l_mult_1352 + (l_mult_1036 + (l_mult_1037 + (l_mult_1504 + (l_mult_1379 + (l_mult_1538 + (l_mult_1071 + (l_mult_1417 + l_mult_1045)))))))))),
                l_mult_1610 + (-0.4320e4 * l_prod_521 + (0.11664e5 * l_prod_976 + (l_mult_1034 + (l_mult_1353 + (l_mult_1348 + (l_mult_1118 + l_sum_1621)))))))) + vdot4(
            v_916, vcons4(l_sum_993, 0.e0, 0.e0, 0.e0)))))))))))))))))))));
        double l_vdot_1952 = vdot4(v_917, vcons4(0.e0, 0.e0, 0.e0, 0.e0)) + (vdot4(v_918,
            vcons4(0.e0, l_sum_1356,
                l_mult_1357 + (l_mult_1010 + (l_mult_1358 + (l_mult_1461 + (l_mult_1360 + (l_mult_1462 + l_sum_1622))))),
                l_mult_1364 + (l_mult_1004 + (l_mult_1472 + (l_mult_1366 + (l_mult_1272 + (l_mult_1008 + l_sum_1624))))))) + (vdot4(
            v_919,
            vcons4(l_mult_1357 + (l_mult_998 + (l_mult_1625 + (l_mult_1398 + l_sum_1626))), l_sum_1627, 0.e0, 0.e0)) + (vdot4(
            v_920, vcons4(0.e0, 0.e0, 0.e0, l_sum_1634)) + (vdot4(v_921,
            vcons4(l_sum_1639, l_sum_1643, l_sum_1646, l_sum_1648)) + (vdot4(v_922,
            vcons4(l_sum_1091, l_sum_1102, l_sum_1112, l_sum_1030)) + (vdot4(v_923,
            vcons4(0.e0, l_sum_1442, l_sum_1453, l_sum_1459)) + (vdot4(v_924,
            vcons4(l_sum_1464, l_sum_1467, l_mult_1028 + (l_mult_1359 + (l_mult_1361 + l_mult_990)),
                l_mult_1649 + (l_mult_1365 + (l_mult_1650 + (l_mult_1026 + l_sum_1021))))) + (vdot4(v_925,
            vcons4(l_mult_1649 + (l_mult_1651 + (l_mult_1371 + l_sum_1027)), l_sum_1030,
                l_mult_1649 + (l_mult_1183 + (l_mult_1650 + (l_mult_1470 + l_sum_1652))),
                l_mult_1653 + (l_mult_1187 + (l_mult_1654 + (l_mult_1058 + l_sum_1657))))) + (vdot4(v_926,
            vcons4(l_mult_1649 + (l_mult_1183 + (l_mult_1651 + (l_mult_1182 + l_sum_1658))),
                l_mult_1649 + (l_mult_1181 + (l_mult_1093 + l_sum_1659)),
                l_mult_1649 + (l_mult_1181 + (l_mult_1093 + l_sum_1660)), l_sum_1661)) + (vdot4(v_927,
            vcons4(
                l_mult_1662 + (l_mult_1207 + (l_mult_1663 + (l_mult_1664 + (l_mult_1665 + (l_mult_1281 + (l_mult_1483 + (l_mult_1666 + (l_mult_1283 + (l_mult_1667 + (l_mult_1668 + (l_mult_1487 + (l_mult_1144 + (l_mult_1669 + (l_mult_1670 + (l_mult_1100 + (-0.17496e5 * l_prod_945 + (l_mult_1490 + (l_mult_1671 + l_mult_990)))))))))))))))))),
                l_mult_1672 + (l_mult_1249 + (l_mult_1673 + (l_mult_1674 + (l_mult_1293 + (l_mult_1379 + (l_mult_1675 + (l_mult_1313 + (l_mult_1676 + (l_mult_1677 + (l_mult_1417 + (l_mult_1678 + (l_mult_1085 + (l_mult_1086 + l_sum_1680))))))))))))),
                l_mult_1681 + (l_mult_1270 + (l_mult_1682 + (l_mult_1302 + (l_mult_1683 + (l_mult_1210 + (l_mult_1667 + l_sum_1027)))))),
                l_sum_1686)) + (vdot4(v_928,
            vcons4(
                l_mult_1672 + (l_mult_1222 + (l_mult_1687 + (l_mult_1250 + (l_mult_1688 + (l_mult_1293 + (l_mult_1495 + (l_mult_1685 + (l_mult_1080 + (l_mult_1677 + (l_mult_1497 + (l_mult_1125 + (l_mult_1108 + (l_mult_1085 + (l_mult_1679 + l_mult_1049)))))))))))))),
                l_mult_1689 + (l_mult_1257 + (l_mult_1645 + (l_mult_1690 + (l_mult_1320 + (l_mult_1400 + (l_mult_1029 + (l_mult_1096 + l_sum_1657))))))),
                l_mult_576 + (l_mult_1067 + (l_mult_1691 + (l_mult_1105 + (l_mult_1692 + l_mult_1040)))),
                l_mult_1681 + (l_mult_1234 + (l_mult_1693 + (l_mult_1664 + (l_mult_1694 + (l_mult_1302 + (l_mult_1502 + l_sum_1659)))))))) + (vdot4(
            v_929,
            vcons4(l_mult_576 + (l_mult_1265 + (l_mult_1695 + (l_mult_1454 + (l_mult_1105 + l_mult_1413)))),
                l_mult_1684 + (l_mult_1241 + (l_mult_1673 + l_mult_1069)),
                l_mult_1542 + (l_mult_1203 + (l_mult_1696 + (l_mult_1619 + (l_mult_1644 + (l_mult_1662 + (l_mult_1207 + (l_mult_1663 + (l_mult_1664 + (0.6426e4 * l_prod_515 + (l_mult_1193 + (l_mult_1400 + (-0.5832e4 * l_prod_966 + (l_mult_1210 + (l_mult_981 + (l_mult_1545 + (0.51408e5 * l_prod_511 + (l_mult_1212 + (l_mult_1697 + (l_mult_1698 + (l_mult_1346 + (l_mult_1215 + (l_mult_1669 + (l_mult_1670 + (l_mult_1699 + (l_mult_1547 + (l_mult_1700 + (l_mult_1701 + (l_mult_1702 + (l_mult_1148 + (l_mult_1409 + (l_mult_1550 + (l_mult_1703 + (l_mult_1704 + l_mult_1463))))))))))))))))))))))))))))))))),
                l_mult_1552 + (l_mult_1246 + (l_mult_1705 + (l_mult_1647 + (l_mult_1672 + (l_mult_1249 + (l_mult_1673 + (l_mult_1706 + (l_mult_1105 + (l_mult_1692 + (l_mult_1555 + (l_mult_1707 + (l_mult_1252 + (l_mult_1631 + (l_mult_1708 + (l_mult_1124 + (l_mult_1125 + (l_mult_1678 + (l_mult_1085 + (l_mult_1109 + (l_mult_1557 + (l_mult_1709 + (l_mult_1710 + (l_mult_1711 + (l_mult_1230 + (l_mult_1391 + (l_mult_1560 + (-0.93312e5 * l_prod_941 + (l_mult_1712 + l_mult_1458)))))))))))))))))))))))))))))) + (vdot4(
            v_930,
            vcons4(
                l_mult_1561 + (l_mult_1268 + (l_mult_1005 + (l_mult_1681 + (l_mult_1270 + (l_mult_1365 + (l_mult_1562 + (l_mult_1713 + (l_mult_1273 + (l_mult_1714 + (l_mult_1564 + (l_mult_1026 + (l_mult_1565 + (l_mult_1715 + (l_mult_1447 + (l_mult_1716 + (l_mult_1490 + (l_mult_988 + (l_mult_1567 + (l_mult_1703 + (l_mult_1704 + l_mult_1452)))))))))))))))))))),
                l_mult_1568 + (l_mult_1277 + (l_mult_1684 + (l_mult_1569 + (l_mult_1717 + (l_mult_1718 + (l_mult_1571 + (l_mult_1719 + (l_mult_1720 + (l_mult_1573 + (l_mult_1633 + (l_mult_1089 + l_mult_1441))))))))))),
                l_mult_1552 + (l_mult_1219 + (l_mult_1721 + (l_mult_1722 + (l_mult_1640 + (l_mult_1672 + (l_mult_1222 + (l_mult_1687 + (l_mult_1250 + (l_mult_1706 + (l_mult_1199 + (l_mult_1379 + (l_mult_1692 + (l_mult_1040 + (l_mult_1576 + (l_mult_1707 + (l_mult_1225 + (-0.93312e5 * l_prod_959 + (l_mult_1723 + (l_mult_1124 + (l_mult_1228 + (l_mult_1108 + (l_mult_1085 + (l_mult_1577 + (l_mult_1724 + (l_mult_1710 + (l_mult_1720 + (l_mult_1438 + l_sum_1725))))))))))))))))))))))))))),
                l_mult_1580 + (l_mult_1254 + (l_mult_1726 + (l_mult_1233 + (l_mult_1689 + (l_mult_1257 + (l_mult_1645 + (l_mult_1654 + (l_mult_1058 + (l_mult_1582 + (0.34992e5 * l_prod_511 + (l_mult_1260 + (l_mult_1520 + (l_mult_1727 + (l_mult_1583 + (l_mult_1144 + (l_mult_1637 + (l_mult_984 + (l_mult_1584 + (l_mult_1728 + (l_mult_1701 + (l_mult_1361 + (l_mult_1449 + (l_mult_1586 + l_mult_1450))))))))))))))))))))))))) + (vdot4(
            v_931,
            vcons4(
                l_mult_1587 + (l_mult_1274 + (l_mult_1729 + (l_mult_576 + (l_mult_1067 + (l_mult_1589 + (l_mult_1730 + (l_mult_1243 + (l_mult_1731 + (l_mult_1417 + (l_mult_1591 + (l_mult_1732 + (l_mult_1437 + (l_mult_1679 + (l_mult_1049 + l_sum_1725)))))))))))))),
                l_mult_1561 + (l_mult_1231 + (l_mult_1733 + (l_mult_1734 + (l_mult_1635 + (l_mult_1681 + (l_mult_1234 + (l_mult_1693 + (l_mult_1664 + (l_mult_1365 + (l_mult_1182 + (l_mult_979 + (l_mult_1595 + (l_mult_1713 + (l_mult_1237 + (l_mult_1697 + (l_mult_1735 + (l_mult_1564 + (l_mult_1196 + (l_mult_1596 + (l_mult_1261 + l_mult_1447)))))))))))))))))))),
                l_mult_1587 + (l_mult_1262 + (l_mult_1736 + (l_mult_1647 + (l_mult_576 + (l_mult_1265 + (l_mult_1695 + (l_mult_1598 + (l_mult_1730 + (l_mult_1267 + (l_mult_1631 + (l_mult_1107 + (l_mult_1417 + (l_mult_1083 + (l_mult_1599 + (l_mult_1719 + l_mult_1437))))))))))))))),
                l_mult_1568 + (l_mult_1239 + (l_mult_1737 + (l_mult_1033 + (l_mult_1629 + (l_mult_1684 + (l_mult_1241 + (l_mult_1673 + (l_mult_1069 + (l_mult_1601 + (l_mult_1717 + (l_mult_1243 + l_mult_1435))))))))))))) + (vdot4(
            v_932,
            vcons4(
                l_mult_1662 + (l_mult_1191 + (-0.17496e5 * l_prod_972 + (l_mult_978 + (l_mult_1665 + (l_mult_1193 + (l_mult_1502 + (l_mult_1666 + (l_mult_1096 + (l_mult_1667 + (l_mult_1698 + (l_mult_1487 + (l_mult_1196 + (l_mult_1738 + (l_mult_1670 + (l_mult_1406 + (l_mult_1702 + (l_mult_1449 + (l_mult_1149 + l_mult_1704)))))))))))))))))),
                l_mult_1672 + (l_mult_1198 + (l_mult_1695 + (l_mult_1688 + (l_mult_1105 + (l_mult_1685 + (l_mult_1708 + (l_mult_1497 + (l_mult_1083 + (l_mult_1126 + (l_mult_1085 + (l_mult_1086 + (l_mult_1711 + (l_mult_1438 + (l_mult_1739 + l_mult_1712)))))))))))))),
                l_mult_1681 + (l_mult_1183 + (l_mult_1694 + (l_mult_1714 + (l_mult_1470 + (l_mult_1740 + (l_mult_1716 + (l_mult_987 + (l_mult_1671 + l_mult_1704)))))))),
                l_mult_1684 + (l_mult_1718 + (l_mult_1720 + l_mult_1089)))) + (vdot4(v_933,
            vcons4(
                l_mult_1672 + (l_mult_1198 + (l_mult_1695 + (l_mult_1674 + (l_mult_1199 + (l_mult_1413 + (l_mult_1675 + (l_mult_1080 + (l_mult_1676 + (l_mult_1723 + (l_mult_1417 + (l_mult_1126 + (l_mult_1085 + (l_mult_1741 + l_sum_1742))))))))))))),
                l_mult_1689 + (l_mult_1187 + (l_mult_1690 + (l_mult_1058 + (l_mult_1029 + (l_mult_1727 + (l_mult_1062 + (l_mult_1743 + (l_mult_984 + (l_mult_1100 + l_sum_1744))))))))),
                l_mult_576 + (l_mult_1454 + (l_mult_1731 + (l_mult_1108 + l_sum_1680))),
                l_mult_1681 + (l_mult_1183 + (l_mult_1682 + (l_mult_1182 + (l_mult_1683 + (l_mult_980 + (l_mult_1667 + l_sum_1745)))))))) + (vdot4(
            v_934,
            vcons4(l_mult_576 + (l_mult_1691 + (l_mult_1692 + l_sum_1110)), l_sum_1686,
                0.4320e4 * l_prod_518 + (l_mult_1339 + (0.58320e5 * l_prod_972 + (l_mult_1630 + (-0.15984e5 * l_prod_515 + (l_mult_1341 + (l_mult_1495 + (0.19440e5 * l_prod_966 + (l_mult_1313 + (l_mult_1081 + (l_mult_1746 + (l_mult_1604 + (l_mult_1228 + (l_mult_1747 + (-0.186624e6 * l_prod_952 + (l_mult_1741 + (0.58320e5 * l_prod_945 + (l_mult_1230 + (l_mult_1739 + l_mult_1440)))))))))))))))))),
                l_mult_1748 + (l_mult_1351 + (l_mult_1645 + (l_mult_1749 + (l_mult_1302 + (l_mult_1029 + (l_mult_1750 + (l_mult_1609 + (l_mult_1144 + (l_mult_1751 + (l_mult_1670 + (l_mult_1100 + (l_mult_1289 + (l_mult_1148 + (l_mult_1149 + l_mult_1451)))))))))))))))) + (vdot4(
            v_935,
            vcons4(
                l_mult_1752 + (l_mult_1353 + (l_mult_1753 + (l_mult_1754 + (l_mult_1614 + (l_mult_1755 + (0.34992e5 * l_prod_945 + (l_mult_1390 + (l_mult_1088 + l_mult_1440)))))))),
                l_mult_1748 + (l_mult_1351 + (l_mult_1645 + (0.13284e5 * l_prod_515 + (l_mult_1345 + (l_mult_1400 + (l_mult_1544 + (l_mult_1283 + (l_mult_1097 + (l_mult_1756 + (l_mult_1564 + (l_mult_1751 + (l_mult_1670 + (l_mult_1406 + l_sum_1744))))))))))))),
                l_mult_1757 + (l_mult_1265 + (l_mult_1758 + (l_mult_1105 + (l_mult_1685 + (l_mult_1759 + (l_mult_1417 + (l_mult_1084 + (l_mult_1085 + (l_mult_1086 + l_sum_1742))))))))),
                l_mult_1752 + (l_mult_1353 + (-0.4320e4 * l_prod_515 + (l_mult_1349 + (0.11664e5 * l_prod_966 + (l_mult_1119 + (l_mult_1081 + l_sum_1761)))))))) + (vdot4(
            v_936,
            vcons4(
                l_mult_1748 + (l_mult_1343 + (l_mult_1208 + (l_mult_1209 + (l_mult_1749 + (l_mult_1345 + (l_mult_1483 + (l_mult_1029 + (l_mult_1096 + (l_mult_1756 + (l_mult_1609 + (l_mult_1215 + (l_mult_1740 + (l_mult_1670 + (l_mult_1361 + l_mult_1449)))))))))))))),
                l_mult_1757 + (l_mult_1352 + (l_mult_1673 + (l_mult_1691 + (l_mult_1105 + (l_mult_1759 + (l_mult_1161 + (l_mult_1125 + (l_mult_1108 + (l_mult_1085 + (l_mult_1720 + l_mult_1438)))))))))),
                l_mult_1757 + (l_mult_1352 + (l_mult_1673 + (l_mult_1758 + (l_mult_1078 + (l_mult_1379 + (l_mult_1685 + (l_mult_1080 + (l_mult_1731 + (l_mult_1417 + (l_mult_1108 + l_mult_1085)))))))))),
                l_mult_1752 + (l_mult_1347 + (0.34992e5 * l_prod_972 + (l_mult_1630 + (l_mult_1753 + (l_mult_1349 + (l_mult_1039 + (l_mult_1760 + (l_mult_1614 + l_mult_1045)))))))))) + vdot4(
            v_916, vcons4(l_sum_993, 0.e0, 0.e0, 0.e0)))))))))))))))))))));
        double l_vdot_1953 = vdot4(v_917,
            vcons4(l_mult_995 + (l_mult_1762 + (l_mult_1302 + l_mult_980)),
                l_mult_998 + (l_mult_999 + (l_mult_1763 + (l_mult_1645 + l_sum_1188))),
                l_mult_1004 + (l_mult_1005 + (l_mult_1006 + l_sum_1190)), l_sum_1012)) + (vdot4(v_918,
            vcons4(0.e0, 0.e0, 0.e0, 0.e0)) + (vdot4(v_919, vcons4(0.e0, 0.e0, 0.e0, l_sum_1766)) + (vdot4(v_920,
            vcons4(l_mult_1366 + (l_mult_1025 + (l_mult_1596 + (l_mult_1368 + (l_mult_1767 + l_mult_990)))),
                l_mult_1358 + (l_mult_1018 + (l_mult_1637 + l_sum_1768)),
                l_mult_1355 + (l_mult_1014 + (l_mult_1740 + l_mult_985)), l_sum_1051)) + (vdot4(v_921,
            vcons4(l_sum_1063, l_sum_1074, l_sum_1012, 0.e0)) + (vdot4(v_922,
            vcons4(
                l_mult_1113 + (l_mult_1114 + (l_mult_1032 + (l_mult_1066 + (l_mult_1115 + (0.6264e4 * l_prod_518 + (l_mult_1769 + (l_mult_1223 + (l_mult_1037 + (-0.20088e5 * l_prod_515 + (0.62208e5 * l_prod_969 + (l_mult_1770 + (0.25920e5 * l_prod_966 + (l_mult_1771 + (l_mult_1676 + (l_mult_1433 + (l_mult_1434 + (l_mult_1042 + (l_mult_1435 + (l_mult_1123 + (l_mult_1124 + (l_mult_1125 + (0.62208e5 * l_prod_953 + (l_mult_1772 + (l_mult_1773 + (l_mult_1774 + (l_mult_1457 + (l_mult_1048 + (l_mult_1310 + (l_mult_1438 + (l_mult_1130 + (l_mult_1775 + (l_mult_1633 + (l_mult_1776 + l_mult_1777))))))))))))))))))))))))))))))))),
                l_mult_1134 + (l_mult_1135 + (l_mult_1053 + (l_mult_1136 + (l_mult_977 + (-0.12447e5 * l_prod_518 + (l_mult_1778 + (l_mult_1636 + (l_mult_1057 + (0.44388e5 * l_prod_515 + (-0.112752e6 * l_prod_969 + (l_mult_1483 + (-0.61560e5 * l_prod_966 + (l_mult_1779 + (l_mult_1780 + (l_mult_1444 + (l_mult_1445 + (l_mult_1060 + (l_mult_982 + (l_mult_1142 + (l_mult_1143 + (l_mult_1144 + (-0.112752e6 * l_prod_953 + (l_mult_1287 + (l_mult_1781 + (l_mult_1782 + (l_mult_1638 + (l_mult_986 + (l_mult_1448 + (l_mult_1449 + (l_mult_1149 + (l_mult_1783 + (l_mult_989 + (l_mult_1784 + l_mult_991))))))))))))))))))))))))))))))))),
                l_mult_1153 + (l_mult_1154 + (l_mult_1065 + (l_mult_1155 + (0.13392e5 * l_prod_518 + (l_mult_1785 + (l_mult_1036 + (l_mult_1069 + (-0.52272e5 * l_prod_515 + (0.101088e6 * l_prod_969 + (l_mult_1770 + (0.77760e5 * l_prod_966 + (l_mult_1786 + (-0.38880e5 * l_prod_963 + (l_mult_1455 + (l_mult_1456 + (l_mult_1071 + (l_mult_1160 + (l_mult_1161 + (l_mult_1083 + (0.101088e6 * l_prod_953 + (l_mult_1772 + (l_mult_1787 + (l_mult_1788 + (l_mult_1641 + (l_mult_1389 + (l_mult_1049 + (l_mult_1130 + (l_mult_1789 + l_mult_1089)))))))))))))))))))))))))))),
                l_mult_1164 + (l_mult_1165 + (l_mult_999 + (-0.8289e4 * l_prod_518 + (l_mult_1790 + (l_mult_1645 + (0.34668e5 * l_prod_515 + (-0.44712e5 * l_prod_969 + (l_mult_979 + (-0.55080e5 * l_prod_966 + (l_mult_1791 + (l_mult_1780 + (l_mult_1460 + (l_mult_1461 + (l_mult_1169 + (l_mult_1170 + (-0.44712e5 * l_prod_953 + (l_mult_984 + (l_mult_1792 + l_sum_1768)))))))))))))))))))) + (vdot4(
            v_923,
            vcons4(
                l_mult_1173 + (l_mult_1174 + (0.2808e4 * l_prod_518 + (l_mult_1793 + (-0.12312e5 * l_prod_515 + (l_mult_1349 + (l_mult_1380 + (l_mult_1040 + (l_mult_1676 + (l_mult_1465 + (l_mult_1177 + (l_mult_1755 + l_mult_1109))))))))))),
                l_sum_1797, l_sum_1801, l_sum_1803)) + (vdot4(v_924, vcons4(l_sum_1766, 0.e0, 0.e0, l_sum_1804)) + (vdot4(
            v_925,
            vcons4(l_mult_1469 + (l_mult_1062 + l_sum_1805), l_mult_1461 + (l_mult_1170 + l_mult_984), 0.e0,
                l_sum_1806)) + (vdot4(v_926, vcons4(l_mult_1469 + (l_mult_1008 + l_sum_1186), 0.e0, l_sum_1009, 0.e0)) + (vdot4(
            v_927,
            vcons4(
                l_mult_1203 + (l_mult_1204 + (l_mult_1205 + (l_mult_1206 + (l_mult_1807 + (l_mult_1808 + (l_mult_1809 + (l_mult_1281 + (l_mult_1810 + (l_mult_1791 + (l_mult_1518 + (l_mult_1519 + (l_mult_1520 + (l_mult_1214 + (l_mult_1215 + (l_mult_1670 + (l_mult_1811 + (l_mult_1812 + (l_mult_1449 + l_mult_1813)))))))))))))))))),
                l_mult_1219 + (l_mult_1220 + (l_mult_1221 + (l_mult_1034 + (l_mult_1814 + (l_mult_1815 + (l_mult_1816 + (0.124416e6 * l_prod_969 + (l_mult_1817 + (l_mult_1786 + (l_mult_1523 + (l_mult_1122 + (l_mult_1043 + (l_mult_1227 + (l_mult_1228 + (l_mult_1297 + (l_mult_1818 + (l_mult_1048 + (l_mult_1438 + l_mult_1050)))))))))))))))))),
                l_mult_1231 + (l_mult_1232 + (l_mult_1233 + (l_mult_1819 + (l_mult_1820 + (l_mult_1057 + (-0.108864e6 * l_prod_969 + (l_mult_1810 + (l_mult_1779 + (l_mult_1527 + (l_mult_1528 + (l_mult_1238 + (l_mult_1196 + (l_mult_1670 + l_sum_1805))))))))))))),
                l_mult_1239 + (l_mult_1240 + (l_mult_1821 + (l_mult_1822 + (l_mult_1293 + (l_mult_1039 + (l_mult_1771 + (l_mult_1530 + (l_mult_1244 + l_mult_1046)))))))))) + (vdot4(
            v_928,
            vcons4(
                l_mult_1246 + (l_mult_1220 + (l_mult_1247 + (l_mult_1248 + (l_mult_1823 + (l_mult_1824 + (l_mult_1816 + (l_mult_1349 + (l_mult_1770 + (l_mult_1532 + (l_mult_1122 + (l_mult_1533 + l_sum_1826))))))))))),
                l_mult_1254 + (l_mult_1255 + (l_mult_1256 + (l_mult_1054 + (l_mult_1827 + (l_mult_1828 + (l_mult_1809 + (l_mult_1829 + (l_mult_1483 + (l_mult_1535 + (l_mult_1536 + (l_mult_1061 + (l_mult_1170 + (l_mult_1144 + l_sum_1003))))))))))))),
                l_mult_1262 + (l_mult_1263 + (l_mult_1264 + (l_mult_1830 + (l_mult_1831 + (l_mult_1037 + (l_mult_1349 + (l_mult_1770 + (l_mult_1538 + (l_mult_1071 + l_sum_1202))))))))),
                l_mult_1268 + (l_mult_1232 + (l_mult_1269 + (l_mult_1206 + (l_mult_1832 + (l_mult_1833 + (l_mult_1057 + l_sum_1541)))))))) + (vdot4(
            v_929,
            vcons4(
                l_mult_1274 + (l_mult_1263 + (l_mult_1275 + (l_mult_1034 + (l_mult_1834 + (l_mult_1822 + (l_mult_1037 + l_sum_1073)))))),
                l_sum_1278,
                l_mult_1485 + (l_mult_1486 + (l_mult_1061 + (l_mult_1195 + (l_mult_1196 + (l_mult_984 + (l_mult_1811 + (l_mult_1812 + (l_mult_1490 + l_mult_1798)))))))),
                l_mult_1496 + (l_mult_1071 + (l_mult_1200 + (l_mult_1818 + (l_mult_1048 + (l_mult_1049 + l_mult_1794))))))) + (vdot4(
            v_930,
            vcons4(l_sum_1804, 0.e0,
                l_mult_1496 + (l_mult_1508 + (l_mult_1043 + (l_mult_1200 + (l_mult_1083 + l_sum_1825)))), l_sum_1806)) + (vdot4(
            v_931, vcons4(0.e0, l_sum_1009, 0.e0, 0.e0)) + (vdot4(v_932,
            vcons4(
                l_mult_1545 + (l_mult_1518 + (l_mult_1486 + (l_mult_1540 + (l_mult_1285 + (l_mult_1214 + (l_mult_1144 + (l_mult_1738 + (l_mult_1670 + (l_mult_1792 + (l_mult_1835 + (l_mult_1216 + (l_mult_1812 + (l_mult_1548 + (l_mult_1148 + (l_mult_1290 + (l_mult_1836 + (l_mult_1450 + (l_mult_1837 + l_mult_1838)))))))))))))))))),
                l_mult_1576 + (l_mult_1532 + (l_mult_1071 + (l_mult_1315 + (l_mult_1244 + (l_mult_1755 + (l_mult_1839 + (l_mult_1245 + (l_mult_1048 + (l_mult_1578 + (l_mult_1438 + (l_mult_1130 + (l_mult_1840 + (l_mult_1526 + (l_mult_1841 + l_mult_1842)))))))))))))),
                l_mult_1595 + (l_mult_1539 + (l_mult_1332 + (l_mult_1843 + (l_mult_997 + (l_mult_1324 + (l_mult_1844 + (l_mult_1813 + (l_mult_1784 + l_mult_1838)))))))),
                l_sum_1847)) + (vdot4(v_933,
            vcons4(
                l_mult_1555 + (l_mult_1523 + (l_mult_1508 + (l_mult_1072 + (l_mult_1295 + (l_mult_1227 + (l_mult_1125 + (0.124416e6 * l_prod_953 + (l_mult_1297 + (l_mult_1787 + (l_mult_1839 + (l_mult_1245 + (l_mult_1048 + (l_mult_1558 + (l_mult_1230 + (l_mult_1299 + (l_mult_1848 + (l_mult_1794 + (l_mult_1841 + l_mult_1795)))))))))))))))))),
                l_mult_1582 + (l_mult_1535 + (l_mult_1403 + (l_mult_1322 + (l_mult_1170 + (l_mult_1333 + (l_mult_1849 + (l_mult_1850 + (l_mult_986 + (l_mult_1585 + (l_mult_1449 + (l_mult_1149 + (l_mult_1851 + (l_mult_1798 + (l_mult_1837 + l_mult_1799)))))))))))))),
                l_mult_1598 + (l_mult_1070 + (l_mult_1336 + (l_mult_1852 + (l_mult_1641 + (l_mult_1572 + (l_mult_1853 + (l_mult_1050 + l_sum_1796))))))),
                l_mult_1562 + (l_mult_1527 + (l_mult_1008 + (l_mult_1304 + (l_mult_1238 + (l_mult_983 + (-0.108864e6 * l_prod_953 + (l_mult_1670 + (l_mult_1781 + (l_mult_1843 + (l_mult_997 + (l_mult_1566 + (l_mult_1490 + (l_mult_1290 + (l_mult_1586 + l_mult_1784)))))))))))))))) + (vdot4(
            v_934,
            vcons4(
                l_mult_1589 + (l_mult_1538 + (l_mult_1328 + (l_mult_1200 + (l_mult_1755 + (l_mult_1852 + (l_mult_1641 + (l_mult_1592 + (l_mult_1049 + (l_mult_1130 + (l_mult_1846 + l_mult_1776)))))))))),
                l_mult_1569 + (l_mult_1530 + (l_mult_1308 + (l_mult_1244 + (l_mult_1126 + (l_mult_1046 + (l_mult_1773 + (l_mult_1845 + (l_mult_1572 + l_mult_1088)))))))),
                l_mult_1602 + (l_mult_1603 + (l_mult_1533 + (l_mult_1342 + (l_mult_1228 + (l_mult_1772 + (l_mult_1854 + (-0.93312e5 * l_prod_946 + (l_mult_1230 + l_mult_1526)))))))),
                l_mult_1615 + (l_mult_1528 + (l_mult_1170 + (l_mult_1855 + (l_mult_1812 + (l_mult_1449 + l_mult_1450))))))) + (vdot4(
            v_935,
            vcons4(l_mult_1618 + (l_mult_1178 + l_mult_1633),
                l_mult_1607 + (l_mult_1608 + (l_mult_1061 + (l_mult_1346 + (l_mult_1215 + (l_mult_1287 + (l_mult_1855 + (l_mult_1812 + (l_mult_1148 + l_mult_1798)))))))),
                l_mult_1617 + (l_mult_1071 + (l_mult_1244 + (l_mult_1457 + (l_mult_1048 + (l_mult_1438 + l_mult_1794))))),
                l_mult_1612 + (l_mult_1613 + (l_mult_1350 + (l_mult_1045 + (l_mult_1772 + (l_mult_1178 + l_mult_1390))))))) + (vdot4(
            v_936,
            vcons4(
                l_mult_1615 + (l_mult_1608 + (l_mult_1520 + (l_mult_1170 + (l_mult_1144 + (l_mult_997 + l_mult_1812))))),
                l_mult_1538 + (l_mult_1071 + l_sum_1825), l_mult_1617 + (l_mult_1042 + (l_mult_1043 + l_sum_1826)),
                l_sum_1621)) + vdot4(v_916,
            vcons4(l_sum_993, 0.e0,
                l_mult_994 + (-0.405e3 * l_prod_518 + (0.1836e4 * l_prod_515 + (-0.3240e4 * l_prod_966 + l_mult_981))),
                0.e0)))))))))))))))))))));
        double l_vdot_1954 = vdot4(v_917,
            vcons4(l_sum_1374,
                l_mult_1357 + (l_mult_1010 + (l_mult_1016 + (l_mult_1167 + (l_mult_1370 + (l_mult_1424 + l_sum_1658))))),
                l_mult_1364 + (l_mult_1004 + (l_mult_1472 + (l_mult_1022 + (l_mult_1832 + (l_mult_1185 + l_sum_1660))))),
                l_mult_1357 + (l_mult_998 + (l_mult_1625 + (l_mult_1398 + l_sum_1661))))) + (vdot4(v_918,
            vcons4(l_sum_1627, 0.e0, 0.e0, 0.e0)) + (vdot4(v_919, vcons4(0.e0, 0.e0, 0.e0, 0.e0)) + (vdot4(v_920,
            vcons4(0.e0, 0.e0, 0.e0, l_sum_1634)) + (vdot4(v_921,
            vcons4(l_sum_1639, l_sum_1643, l_sum_1646, l_sum_1648)) + (vdot4(v_922,
            vcons4(l_sum_1395, l_sum_1410, l_sum_1422, l_sum_1426)) + (vdot4(v_923,
            vcons4(l_sum_1430, l_sum_1797, l_sum_1801, l_sum_1803)) + (vdot4(v_924,
            vcons4(l_sum_1766, 0.e0, l_sum_1766,
                l_mult_1856 + (l_mult_1025 + (l_mult_1857 + (l_mult_1368 + l_sum_1363))))) + (vdot4(v_925,
            vcons4(l_mult_1856 + (l_mult_1650 + (l_mult_1637 + l_sum_1369)), l_sum_1373,
                l_mult_1856 + (l_mult_1007 + (l_mult_1857 + (l_mult_1623 + l_sum_1622))),
                l_mult_1858 + (l_mult_1475 + (l_mult_1655 + (l_mult_1062 + l_sum_1859))))) + (vdot4(v_926,
            vcons4(l_mult_1856 + (l_mult_1007 + (l_mult_1650 + (l_mult_1470 + l_sum_1656))),
                l_mult_1856 + (l_mult_1469 + (l_mult_1403 + l_sum_1624)),
                l_mult_1856 + (l_mult_1469 + (l_mult_1403 + l_sum_1659)), l_sum_1626)) + (vdot4(v_927,
            vcons4(
                l_mult_1542 + (l_mult_1203 + (l_mult_1696 + (l_mult_1619 + (l_mult_1644 + (l_mult_1279 + (0.51408e5 * l_prod_517 + (l_mult_1808 + (l_mult_1860 + (l_mult_1543 + (l_mult_1861 + (l_mult_1862 + (l_mult_1544 + (l_mult_1863 + (l_mult_1425 + (l_mult_1864 + (l_mult_1518 + (l_mult_1865 + (l_mult_1866 + (l_mult_1698 + (l_mult_1346 + (l_mult_1215 + (l_mult_1867 + (l_mult_1287 + (l_mult_1868 + (0.6426e4 * l_prod_507 + (l_mult_1811 + (l_mult_1447 + (l_mult_1869 + (l_mult_1218 + (l_mult_1409 + (-0.5832e4 * l_prod_942 + (l_mult_1813 + (l_mult_1870 + l_mult_991))))))))))))))))))))))))))))))))),
                l_mult_1552 + (l_mult_1246 + (l_mult_1705 + (l_mult_1647 + (l_mult_1291 + (l_mult_1871 + (l_mult_1824 + (l_mult_1630 + (l_mult_1574 + (l_mult_1872 + (l_mult_1873 + (l_mult_1575 + (-0.93312e5 * l_prod_965 + (l_mult_1414 + (l_mult_1874 + (l_mult_1532 + (l_mult_1875 + (l_mult_1708 + (l_mult_1124 + (l_mult_1125 + (l_mult_1876 + (l_mult_1297 + (l_mult_1741 + (l_mult_1877 + (l_mult_1641 + (l_mult_1878 + (l_mult_1390 + (l_mult_1391 + l_sum_1880))))))))))))))))))))))))))),
                l_mult_1561 + (l_mult_1268 + (l_mult_1005 + (l_mult_1300 + (l_mult_1881 + (l_mult_1833 + (l_mult_1593 + (l_mult_1882 + (l_mult_1400 + (l_mult_1594 + (l_mult_1863 + (l_mult_1401 + (l_mult_1883 + (l_mult_1539 + (l_mult_1714 + (l_mult_1564 + (l_mult_1884 + (l_mult_1197 + (l_mult_1868 + l_sum_1369)))))))))))))))))),
                l_mult_1568 + (l_mult_1277 + (l_mult_1306 + (l_mult_1885 + (l_mult_1600 + (l_mult_1886 + (l_mult_1079 + (l_mult_1119 + (l_mult_1382 + l_sum_1889)))))))))) + (vdot4(
            v_928,
            vcons4(
                l_mult_1552 + (l_mult_1219 + (l_mult_1721 + (l_mult_1722 + (l_mult_1640 + (l_mult_1311 + (l_mult_1871 + (l_mult_1815 + (-0.93312e5 * l_prod_971 + (l_mult_1554 + (l_mult_1890 + (l_mult_1873 + (l_mult_1428 + (l_mult_1381 + (l_mult_1874 + (l_mult_1523 + (l_mult_1891 + (l_mult_1533 + (l_mult_1723 + (l_mult_1124 + (l_mult_1228 + (l_mult_1888 + (l_mult_1127 + (l_mult_1877 + (l_mult_1818 + (l_mult_1437 + (l_mult_1419 + (l_mult_1390 + l_sum_1892))))))))))))))))))))))))))),
                l_mult_1580 + (l_mult_1254 + (l_mult_1726 + (l_mult_1233 + (l_mult_1318 + (0.34992e5 * l_prod_517 + (l_mult_1828 + (l_mult_1209 + (l_mult_1581 + (l_mult_1893 + (l_mult_1862 + (l_mult_1303 + (l_mult_1283 + (l_mult_1894 + (l_mult_1535 + (l_mult_1001 + (l_mult_1727 + (l_mult_1583 + (l_mult_1144 + (l_mult_1019 + (l_mult_1146 + l_sum_1859)))))))))))))))))))),
                l_mult_1587 + (l_mult_1274 + (l_mult_1729 + (l_mult_1325 + (l_mult_1895 + (l_mult_1822 + (l_mult_1597 + (l_mult_1896 + (l_mult_1379 + (l_mult_1428 + (l_mult_1381 + (l_mult_558 + (l_mult_1070 + (l_mult_1731 + (l_mult_1417 + (l_mult_1897 + l_mult_1046))))))))))))))),
                l_mult_1561 + (l_mult_1231 + (l_mult_1733 + (l_mult_1734 + (l_mult_1635 + (l_mult_1330 + (l_mult_1881 + (l_mult_1820 + (l_mult_1860 + (l_mult_1023 + (l_mult_1829 + (l_mult_1400 + (l_mult_1883 + (l_mult_1527 + (l_mult_1898 + (l_mult_1866 + (l_mult_1735 + (l_mult_1564 + (l_mult_1196 + l_sum_1624)))))))))))))))))))) + (vdot4(
            v_929,
            vcons4(
                l_mult_1587 + (l_mult_1262 + (l_mult_1736 + (l_mult_1647 + (l_mult_1334 + (l_mult_1895 + (l_mult_1831 + (l_mult_1630 + (l_mult_1588 + (l_mult_1886 + (l_mult_1379 + (l_mult_558 + (l_mult_1538 + (l_mult_1899 + (l_mult_1107 + (l_mult_1417 + l_mult_1083))))))))))))))),
                l_mult_1568 + (l_mult_1239 + (l_mult_1737 + (l_mult_1033 + (l_mult_1629 + (l_mult_1337 + (l_mult_1885 + (l_mult_1822 + (l_mult_1118 + l_sum_1900)))))))),
                l_mult_1864 + (l_mult_1518 + (l_mult_1865 + (l_mult_1866 + (l_mult_1668 + (l_mult_1487 + (l_mult_1144 + (-0.17496e5 * l_prod_953 + (l_mult_1197 + (l_mult_985 + (l_mult_1901 + (l_mult_1216 + (l_mult_1147 + (l_mult_1869 + (l_mult_1218 + (l_mult_1671 + (l_mult_1902 + (l_mult_1450 + (l_mult_1784 + l_mult_1903)))))))))))))))))),
                l_mult_1874 + (l_mult_1532 + (l_mult_1875 + (l_mult_1677 + (l_mult_1417 + (l_mult_1897 + (l_mult_1904 + (l_mult_1245 + (l_mult_1437 + (l_mult_1878 + (l_mult_1390 + (l_mult_1420 + (l_mult_1905 + (l_mult_1526 + (l_mult_1776 + l_mult_1133)))))))))))))))) + (vdot4(
            v_930,
            vcons4(
                l_mult_1883 + (l_mult_1539 + (l_mult_1025 + (l_mult_1906 + (l_mult_997 + (l_mult_1368 + (l_mult_1907 + (l_mult_1813 + (l_mult_990 + l_mult_1903)))))))),
                l_sum_1909,
                l_mult_1874 + (l_mult_1523 + (l_mult_1891 + (l_mult_1533 + (l_mult_1677 + (l_mult_1497 + (l_mult_1125 + (l_mult_1897 + (l_mult_1046 + (l_mult_1910 + (l_mult_1245 + (l_mult_1525 + l_sum_1911))))))))))),
                l_mult_1894 + (l_mult_1535 + (l_mult_1001 + (l_mult_1655 + (l_mult_1062 + (l_mult_1912 + (l_mult_1850 + (l_mult_1447 + (l_mult_1020 + (l_mult_987 + l_sum_1913))))))))))) + (vdot4(
            v_931,
            vcons4(l_mult_558 + (l_mult_1070 + (l_mult_1914 + (l_mult_1641 + l_sum_1892))),
                l_mult_1883 + (l_mult_1527 + (l_mult_1898 + (l_mult_1866 + (l_mult_1025 + (l_mult_1470 + (l_mult_983 + (l_mult_1915 + (l_mult_997 + l_mult_1489)))))))),
                l_mult_558 + (l_mult_1538 + (l_mult_1899 + l_sum_1642)), l_sum_1900)) + (vdot4(v_932,
            vcons4(
                l_mult_1864 + (l_mult_1485 + (-0.17496e5 * l_prod_960 + (l_mult_982 + (l_mult_1698 + (l_mult_1487 + (l_mult_1196 + (l_mult_1867 + (l_mult_1146 + (l_mult_1868 + (l_mult_1901 + (l_mult_1811 + (l_mult_1489 + (l_mult_1289 + (l_mult_1218 + (l_mult_1149 + (l_mult_1902 + (l_mult_1798 + (l_mult_1451 + l_mult_1903)))))))))))))))))),
                l_mult_1874 + (l_mult_1496 + (l_mult_1899 + (l_mult_1723 + (l_mult_1417 + (l_mult_1888 + (l_mult_1904 + (l_mult_1818 + (l_mult_1388 + (l_mult_1310 + (l_mult_1390 + (l_mult_1391 + (l_mult_1905 + (l_mult_1794 + (l_mult_1712 + l_mult_1133)))))))))))))),
                l_mult_1883 + (l_mult_1007 + (l_mult_1735 + (l_mult_1906 + (l_mult_1623 + (l_mult_1015 + (l_mult_1907 + (l_mult_989 + (l_mult_1870 + l_mult_1903)))))))),
                l_sum_1909)) + (vdot4(v_933,
            vcons4(
                l_mult_1874 + (l_mult_1496 + (l_mult_1899 + (l_mult_1708 + (l_mult_1497 + (l_mult_1083 + (l_mult_1876 + (l_mult_1127 + (l_mult_1741 + (l_mult_1910 + (l_mult_1641 + (l_mult_1310 + (l_mult_1390 + (l_mult_1739 + l_sum_1916))))))))))))),
                l_mult_1894 + (l_mult_1475 + (l_mult_1727 + (l_mult_1062 + (l_mult_1019 + (l_mult_1912 + (l_mult_1002 + (l_mult_1917 + (l_mult_987 + (l_mult_1409 + l_sum_1918))))))))),
                l_mult_558 + (l_mult_1107 + (l_mult_1914 + (l_mult_1419 + l_sum_1880))),
                l_mult_1883 + (l_mult_1007 + (l_mult_1714 + (l_mult_1470 + (l_mult_1884 + (l_mult_984 + (l_mult_1868 + (l_mult_1915 + (l_mult_1015 + l_mult_1671)))))))))) + (vdot4(
            v_934,
            vcons4(l_mult_558 + (l_mult_1731 + (l_mult_1897 + l_sum_1421)), l_sum_1889,
                0.4320e4 * l_prod_512 + (l_mult_1602 + (0.58320e5 * l_prod_960 + (l_mult_1631 + (l_mult_1746 + (l_mult_1604 + (l_mult_1228 + (0.58320e5 * l_prod_953 + (l_mult_1297 + (l_mult_1386 + (-0.15984e5 * l_prod_507 + (l_mult_1854 + (l_mult_1525 + (l_mult_1919 + (-0.186624e6 * l_prod_944 + (l_mult_1739 + (0.19440e5 * l_prod_942 + (l_mult_1526 + (l_mult_1712 + l_mult_1795)))))))))))))))))),
                l_mult_1920 + (l_mult_1615 + (l_mult_1001 + (l_mult_1756 + (l_mult_1564 + (l_mult_1019 + (0.13284e5 * l_prod_507 + (l_mult_1855 + (l_mult_1447 + (l_mult_1921 + (l_mult_1218 + (l_mult_1409 + (l_mult_1550 + (l_mult_1450 + (l_mult_1451 + l_mult_1799)))))))))))))))) + (vdot4(
            v_935,
            vcons4(
                l_mult_1922 + (l_mult_1618 + (l_mult_1760 + (-0.4320e4 * l_prod_507 + (l_mult_1178 + (l_mult_1179 + (0.11664e5 * l_prod_942 + (l_mult_1633 + (l_mult_1393 + l_mult_1795)))))))),
                l_mult_1920 + (l_mult_1615 + (l_mult_1001 + (l_mult_1750 + (l_mult_1609 + (l_mult_1144 + (l_mult_1738 + (l_mult_1287 + (l_mult_1406 + (l_mult_1923 + (l_mult_997 + (l_mult_1921 + (l_mult_1218 + (l_mult_1149 + l_sum_1918))))))))))))),
                l_mult_1924 + (l_mult_1538 + (l_mult_1759 + (l_mult_1417 + (l_mult_1888 + (l_mult_1925 + (l_mult_1641 + (l_mult_1389 + (l_mult_1390 + (l_mult_1391 + l_sum_1916))))))))),
                l_mult_1922 + (l_mult_1618 + (l_mult_1754 + (l_mult_1614 + (0.34992e5 * l_prod_953 + (l_mult_1085 + (l_mult_1386 + (l_mult_1926 + (l_mult_1179 + l_mult_1088)))))))))) + (vdot4(
            v_936,
            vcons4(
                l_mult_1920 + (l_mult_1607 + (l_mult_1519 + (l_mult_1520 + (l_mult_1756 + (l_mult_1609 + (l_mult_1215 + (l_mult_1019 + (l_mult_1146 + (l_mult_1923 + (l_mult_1855 + (l_mult_1147 + (l_mult_1015 + (l_mult_1218 + l_sum_1913))))))))))))),
                l_mult_1924 + (l_mult_1617 + (l_mult_1875 + (l_mult_1731 + (l_mult_1417 + (l_mult_1925 + (l_mult_1457 + (l_mult_1437 + l_sum_1911))))))),
                l_mult_1924 + (l_mult_1617 + (l_mult_1875 + (l_mult_1759 + (l_mult_1161 + (l_mult_1125 + (l_mult_1888 + (l_mult_1127 + (l_mult_1914 + (l_mult_1641 + (l_mult_1419 + l_mult_1390)))))))))),
                l_mult_1922 + (l_mult_1612 + (0.34992e5 * l_prod_960 + (l_mult_1631 + (l_mult_1760 + (l_mult_1614 + (l_mult_1045 + (l_mult_1926 + (l_mult_1178 + l_mult_1048)))))))))) + vdot4(
            v_916, vcons4(l_sum_993, 0.e0, 0.e0, 0.e0)))))))))))))))))))));
        double l_vdot_1955 = vdot4(v_917,
            vcons4(0.e0, l_sum_1030,
                l_mult_1022 + (l_mult_1183 + (l_mult_1023 + (l_mult_1182 + (l_mult_1024 + l_mult_980)))),
                l_mult_1016 + (l_mult_1763 + (l_mult_1093 + (l_mult_1017 + (l_mult_1424 + l_mult_979)))))) + (vdot4(
            v_918,
            vcons4(l_mult_1013 + (l_mult_1762 + (l_mult_1235 + l_mult_978)), 0.e0, l_sum_1766,
                l_mult_1366 + (l_mult_1007 + (l_mult_1596 + (l_mult_1623 + (l_mult_1767 + l_mult_989)))))) + (vdot4(
            v_919,
            vcons4(l_mult_1358 + (l_mult_1000 + (l_mult_1403 + l_sum_1927)),
                l_mult_1355 + (l_mult_996 + (l_mult_1528 + l_mult_982)), 0.e0, 0.e0)) + (vdot4(v_920,
            vcons4(0.e0, 0.e0, 0.e0,
                l_mult_1113 + (0.6264e4 * l_prod_522 + (-0.20088e5 * l_prod_521 + (0.25920e5 * l_prod_976 + (l_mult_1507 + (l_mult_1116 + (l_mult_1769 + (0.62208e5 * l_prod_972 + (l_mult_1928 + (l_mult_1077 + (l_mult_1293 + (l_mult_1770 + (l_mult_1106 + (l_mult_1080 + (l_mult_1120 + (l_mult_1433 + (l_mult_1121 + (0.62208e5 * l_prod_960 + (l_mult_1929 + (l_mult_1385 + (l_mult_1124 + (l_mult_1930 + (l_mult_1084 + (l_mult_1127 + (l_mult_1436 + (l_mult_1774 + (l_mult_1245 + (l_mult_1128 + (l_mult_1389 + (l_mult_1438 + (l_mult_1088 + (l_mult_1775 + (l_mult_1794 + (l_mult_1393 + l_mult_1777))))))))))))))))))))))))))))))))))) + (vdot4(
            v_921,
            vcons4(
                l_mult_1134 + (-0.12447e5 * l_prod_522 + (0.44388e5 * l_prod_521 + (-0.61560e5 * l_prod_976 + (l_mult_1931 + (l_mult_1137 + (l_mult_1778 + (-0.112752e6 * l_prod_972 + (l_mult_1932 + (l_mult_1094 + (l_mult_1399 + (l_mult_1483 + (l_mult_1139 + (l_mult_1096 + (l_mult_981 + (l_mult_1444 + (l_mult_1140 + (-0.112752e6 * l_prod_960 + (l_mult_1933 + (l_mult_1404 + (l_mult_1143 + (l_mult_1215 + (l_mult_1099 + (l_mult_1146 + (l_mult_985 + (l_mult_1782 + (l_mult_1446 + (l_mult_1147 + (l_mult_1408 + (l_mult_1449 + (l_mult_988 + (l_mult_1783 + (l_mult_1798 + l_sum_992)))))))))))))))))))))))))))))))),
                l_mult_1153 + (0.13392e5 * l_prod_522 + (-0.52272e5 * l_prod_521 + (0.77760e5 * l_prod_976 + (-0.38880e5 * l_prod_975 + (l_mult_1156 + (l_mult_1785 + (0.101088e6 * l_prod_972 + (l_mult_1934 + (l_mult_1104 + (l_mult_1078 + (l_mult_1770 + (l_mult_1158 + (l_mult_1040 + (l_mult_1455 + (l_mult_1159 + (0.101088e6 * l_prod_960 + (l_mult_1935 + (l_mult_1416 + (l_mult_1161 + (l_mult_1930 + (l_mult_1108 + (l_mult_1046 + (l_mult_1788 + (l_mult_1457 + (l_mult_1128 + (l_mult_1419 + (l_mult_1049 + (l_mult_1789 + l_mult_1050)))))))))))))))))))))))))))),
                l_mult_1164 + (-0.8289e4 * l_prod_522 + (0.34668e5 * l_prod_521 + (-0.55080e5 * l_prod_976 + (l_mult_1931 + (l_mult_1166 + (l_mult_1790 + (-0.44712e5 * l_prod_972 + (l_mult_1936 + (l_mult_1017 + (l_mult_1424 + (l_mult_979 + (l_mult_1460 + (l_mult_1168 + (-0.44712e5 * l_prod_960 + (l_mult_1937 + (l_mult_1359 + (l_mult_1170 + (l_mult_983 + l_sum_1927)))))))))))))))))),
                l_mult_1173 + (0.2808e4 * l_prod_522 + (-0.12312e5 * l_prod_521 + (l_mult_1628 + (l_mult_1507 + (l_mult_1175 + (l_mult_1793 + (l_mult_1348 + (l_mult_1069 + (l_mult_1465 + (l_mult_1176 + (l_mult_1613 + l_mult_1072))))))))))))) + (vdot4(
            v_922, vcons4(l_sum_1091, l_sum_1102, l_sum_1112, l_sum_1030)) + (vdot4(v_923,
            vcons4(0.e0, l_sum_1797, l_sum_1801, l_sum_1803)) + (vdot4(v_924, vcons4(l_sum_1766, 0.e0, 0.e0, 0.e0)) + (vdot4(
            v_925, vcons4(0.e0, 0.e0, l_sum_1938, l_sum_1939)) + (vdot4(v_926,
            vcons4(l_sum_1027, l_mult_1650 + (l_mult_1062 + l_sum_1940),
                l_mult_1650 + (l_mult_1062 + (l_mult_1026 + l_mult_984)), l_mult_1359 + (l_mult_1170 + l_mult_983))) + (vdot4(
            v_927,
            vcons4(
                l_mult_1279 + (l_mult_1807 + (l_mult_1208 + (l_mult_1936 + (l_mult_1280 + (l_mult_1861 + (l_mult_1810 + (l_mult_1282 + (l_mult_1941 + (l_mult_1284 + (l_mult_1698 + (l_mult_1214 + (l_mult_1488 + (l_mult_1738 + (l_mult_1287 + (l_mult_1406 + (l_mult_1869 + (l_mult_1449 + (l_mult_1942 + l_mult_1870)))))))))))))))))),
                l_mult_1311 + (l_mult_1823 + (l_mult_1348 + (l_mult_1292 + (l_mult_1890 + (l_mult_1770 + (l_mult_1312 + (l_mult_1943 + (l_mult_1314 + (l_mult_1723 + (l_mult_1244 + (l_mult_1126 + (l_mult_1127 + (l_mult_1741 + l_sum_1944))))))))))))),
                l_mult_1330 + (l_mult_1832 + (l_mult_1301 + (l_mult_1829 + (l_mult_1331 + (l_mult_1096 + (l_mult_1284 + l_sum_1745)))))),
                l_sum_1338)) + (vdot4(v_928,
            vcons4(
                l_mult_1291 + (l_mult_1814 + (0.124416e6 * l_prod_972 + (l_mult_1934 + (l_mult_1292 + (l_mult_1872 + (l_mult_1817 + (l_mult_1294 + (l_mult_1943 + (l_mult_1081 + (l_mult_1708 + (l_mult_1227 + (l_mult_1228 + (l_mult_1126 + (l_mult_1297 + (l_mult_1086 + (l_mult_1878 + (l_mult_1438 + l_sum_1090))))))))))))))))),
                l_mult_1318 + (l_mult_1827 + (l_mult_1833 + (l_mult_1319 + (l_mult_1893 + (l_mult_1483 + (l_mult_1321 + (l_mult_1941 + (l_mult_1097 + (l_mult_1727 + (l_mult_1170 + (l_mult_1743 + (l_mult_1146 + l_sum_1101)))))))))))),
                l_mult_1334 + (l_mult_1834 + (l_mult_1326 + (l_mult_1886 + (l_mult_1335 + (l_mult_1080 + l_sum_1111))))),
                l_mult_1300 + (l_mult_1819 + (-0.108864e6 * l_prod_972 + (l_mult_1932 + (l_mult_1301 + (l_mult_1882 + (l_mult_1810 + (l_mult_1303 + (l_mult_1096 + (l_mult_1714 + (l_mult_1238 + (l_mult_1488 + (l_mult_1740 + (l_mult_1197 + l_sum_1940))))))))))))))) + (vdot4(
            v_929,
            vcons4(
                l_mult_1325 + (l_mult_1830 + (l_mult_1348 + (l_mult_1326 + (l_mult_1896 + (l_mult_1770 + (l_mult_1327 + (l_mult_1080 + (l_mult_1731 + (l_mult_1200 + (l_mult_1108 + l_mult_1046)))))))))),
                l_mult_1306 + (l_mult_1821 + (l_mult_1223 + (l_mult_1928 + (l_mult_1307 + (l_mult_1886 + (l_mult_1039 + (l_mult_1718 + (l_mult_1244 + l_mult_1083)))))))),
                l_mult_1545 + (l_mult_1211 + (l_mult_1519 + (l_mult_1937 + (l_mult_1698 + (l_mult_1214 + (l_mult_1488 + (l_mult_1669 + (l_mult_1146 + (l_mult_1699 + (l_mult_1835 + (l_mult_1700 + (l_mult_1217 + (l_mult_1289 + (l_mult_1148 + (l_mult_1942 + (l_mult_1836 + (l_mult_1945 + (l_mult_1451 + l_mult_1838)))))))))))))))))),
                l_mult_1576 + (l_mult_1251 + (l_mult_1613 + (l_mult_1723 + (l_mult_1244 + (l_mult_1108 + (l_mult_1839 + (l_mult_1724 + (l_mult_1128 + (l_mult_1310 + (l_mult_1438 + (l_mult_1088 + (l_mult_1840 + (l_mult_1946 + (l_mult_1712 + l_mult_1842)))))))))))))))) + (vdot4(
            v_930,
            vcons4(
                l_mult_1595 + (l_mult_1272 + (l_mult_1735 + (l_mult_1843 + (l_mult_1261 + (l_mult_1015 + (l_mult_1844 + (l_mult_1798 + (l_mult_1870 + l_mult_1838)))))))),
                l_sum_1847,
                l_mult_1555 + (l_mult_1224 + (0.124416e6 * l_prod_960 + (l_mult_1935 + (l_mult_1708 + (l_mult_1227 + (l_mult_1228 + (l_mult_1678 + (l_mult_1127 + (l_mult_1109 + (l_mult_1839 + (l_mult_1709 + (l_mult_1229 + (l_mult_1310 + (l_mult_1230 + (l_mult_1088 + (l_mult_1848 + (l_mult_1946 + l_sum_1796))))))))))))))))),
                l_mult_1582 + (l_mult_1259 + (l_mult_1273 + (l_mult_1727 + (l_mult_1170 + (l_mult_1637 + (l_mult_1849 + (l_mult_1728 + (l_mult_1147 + (l_mult_1917 + (l_mult_1449 + (l_mult_988 + (l_mult_1851 + (l_mult_1945 + l_sum_1800))))))))))))))) + (vdot4(
            v_931,
            vcons4(
                l_mult_1598 + (l_mult_1276 + (l_mult_1107 + (l_mult_1852 + (l_mult_1719 + (l_mult_1419 + (l_mult_1853 + (l_mult_1794 + l_sum_1802))))))),
                l_mult_1562 + (l_mult_1236 + (-0.108864e6 * l_prod_960 + (l_mult_1933 + (l_mult_1714 + (l_mult_1238 + (l_mult_1488 + (l_mult_1026 + (l_mult_984 + (l_mult_1843 + (l_mult_1715 + (l_mult_1217 + (l_mult_1015 + (l_mult_1490 + (l_mult_1586 + l_mult_1798)))))))))))))),
                l_mult_1589 + (l_mult_1266 + (l_mult_1613 + (l_mult_1731 + (l_mult_1200 + (l_mult_1852 + (l_mult_1732 + (l_mult_1128 + (l_mult_1419 + (l_mult_1049 + (l_mult_1846 + l_mult_1794)))))))))),
                l_mult_1569 + (l_mult_1242 + (l_mult_1122 + (l_mult_1929 + (l_mult_1718 + (l_mult_1244 + (l_mult_1083 + (l_mult_1845 + (l_mult_1719 + l_mult_1048)))))))))) + (vdot4(
            v_932,
            vcons4(
                l_mult_1668 + (l_mult_1195 + (l_mult_983 + (l_mult_1669 + (l_mult_1197 + (l_mult_1100 + (l_mult_1869 + (l_mult_1490 + l_sum_1947))))))),
                l_mult_1677 + (l_mult_1200 + (l_mult_1108 + (l_mult_1878 + (l_mult_1049 + l_sum_1948)))), l_sum_1938,
                0.e0)) + (vdot4(v_933,
            vcons4(l_mult_1677 + (l_mult_1200 + (l_mult_1678 + (l_mult_1046 + l_sum_1949))), l_sum_1939, 0.e0,
                l_sum_1027)) + (vdot4(v_934,
            vcons4(0.e0, 0.e0,
                l_mult_1746 + (l_mult_1342 + (l_mult_1930 + (l_mult_1747 + (l_mult_1297 + (l_mult_1741 + (l_mult_1919 + (l_mult_1230 + (-0.93312e5 * l_prod_943 + l_mult_1712)))))))),
                l_mult_1756 + (l_mult_1170 + (l_mult_1740 + (l_mult_1921 + (l_mult_1449 + (l_mult_1942 + l_mult_1451))))))) + (vdot4(
            v_935,
            vcons4(l_mult_1760 + (l_mult_1179 + l_mult_1393),
                l_mult_1756 + (l_mult_1170 + (l_mult_1751 + (l_mult_1146 + (l_mult_1406 + (l_mult_1015 + l_mult_1942))))),
                l_mult_1731 + (l_mult_1108 + l_sum_1944), l_sum_1761)) + (vdot4(v_936,
            vcons4(
                l_mult_1750 + (l_mult_1346 + (l_mult_1215 + (l_mult_1751 + (l_mult_1287 + (l_mult_1100 + (l_mult_1921 + (l_mult_1148 + l_sum_1947))))))),
                l_mult_1759 + (l_mult_1244 + (l_mult_1108 + (l_mult_1389 + (l_mult_1438 + l_sum_1948)))),
                l_mult_1759 + (l_mult_1244 + (l_mult_1084 + (l_mult_1127 + l_sum_1949))),
                l_mult_1754 + (l_mult_1350 + (l_mult_1930 + (l_mult_1755 + (l_mult_1085 + (l_mult_1179 + l_mult_1390))))))) + vdot4(
            v_916,
            vcons4(l_sum_993, 0.e0, 0.e0,
                l_mult_994 + (-0.405e3 * l_prod_522 + (0.1836e4 * l_prod_521 + (-0.3240e4 * l_prod_976 + l_mult_977))))))))))))))))))))))));
        double l_r_1956 = l_r_717 / l_op1_e3_l_18_746;
        double l_r_1957 = l_r_725 / l_op1_e3_l_18_746;
        double l_r_1958 = l_r_733 / l_op1_e3_l_18_746;
        double l_r_1959 = l_r_736 / l_op1_e3_l_18_746;
        double l_r_1960 = l_r_739 / l_op1_e3_l_18_746;
        double l_r_1961 = l_r_742 / l_op1_e3_l_18_746;
        double l_r_1962 = l_r_743 / l_op1_e3_l_18_746;
        double l_r_1963 = l_r_744 / l_op1_e3_l_18_746;
        double l_r_1964 = l_r_745 / l_op1_e3_l_18_746;
        double l_r_1965 = l_vdot_1950 * l_r_1956 + l_vdot_1951 * l_r_1959 + l_vdot_1952 * l_r_1962;
        double l_r_1966 = l_vdot_1950 * l_r_1957 + l_vdot_1951 * l_r_1960 + l_vdot_1952 * l_r_1963;
        double l_r_1967 = l_vdot_1950 * l_r_1958 + l_vdot_1951 * l_r_1961 + l_vdot_1952 * l_r_1964;
        double l_r_1968 = l_vdot_1951 * l_r_1956 + l_vdot_1953 * l_r_1959 + l_vdot_1954 * l_r_1962;
        double l_r_1969 = l_vdot_1951 * l_r_1957 + l_vdot_1953 * l_r_1960 + l_vdot_1954 * l_r_1963;
        double l_r_1970 = l_vdot_1951 * l_r_1958 + l_vdot_1953 * l_r_1961 + l_vdot_1954 * l_r_1964;
        double l_r_1971 = l_vdot_1952 * l_r_1956 + l_vdot_1954 * l_r_1959 + l_vdot_1955 * l_r_1962;
        double l_r_1972 = l_vdot_1952 * l_r_1957 + l_vdot_1954 * l_r_1960 + l_vdot_1955 * l_r_1963;
        double l_r_1973 = l_vdot_1952 * l_r_1958 + l_vdot_1954 * l_r_1961 + l_vdot_1955 * l_r_1964;
        double l_prod_1974 = l_prod4_938 * l_varAcc_502 * l_prod_506;
        double l_prod_1975 = l_prod4_938 * l_prod_508;
        double l_prod_1976 = l_prod4_938 * l_prod_510;
        double l_prod_1977 = l_prod3_937 * l_prod_514;
        double l_prod_1978 = l_prod3_937 * l_prod_516;
        double l_prod_1979 = l_prod3_937 * l_prod_520;
        double l_prod_1980 = l_prod2_505 * l_prod_949;
        double l_prod_1981 = l_prod2_505 * l_prod_951;
        double l_prod_1982 = l_prod2_505 * l_prod_954;
        double l_prod_1983 = l_prod2_505 * l_prod_958;
        double l_prod_1984 = l_varAcc_502 * l_prod_962;
        double l_prod_1985 = l_varAcc_502 * l_prod_964;
        double l_prod_1986 = l_varAcc_502 * l_prod_967;
        double l_prod_1987 = l_varAcc_502 * l_prod_970;
        double l_prod_1988 = l_varAcc_502 * l_prod_974;
        double l_prod_1989 = 0.1e1 * (l_prod4_961 * l_varAcc_503 * 0.1e1);
        double l_prod_1990 = 0.1e1 * (l_prod4_961 * l_varAcc_504);
        double l_prod_1991 = 0.1e1 * (l_prod3_948 * l_prod2_519);
        double l_prod_1992 = 0.1e1 * (l_prod2_513 * l_prod3_957);
        double l_prod_1993 = 0.1e1 * (l_varAcc_503 * l_prod4_973);
        double l_prod_1994 = 0.1e1 * (0.1e1 * (l_prod4_973 * l_varAcc_504));
        double l_mult_1995 = 0.3888e3 * l_prod_1994;
        double l_mult_1996 = 0.1944e4 * l_prod_1993;
        double l_mult_1997 = 0.3888e4 * l_prod_1992;
        double l_mult_1998 = 0.3888e4 * l_prod_1991;
        double l_mult_1999 = 0.1944e4 * l_prod_1990;
        double l_mult_2000 = 0.3888e3 * l_prod_1989;
        double l_mult_2001 = 0.1944e4 * l_prod_1988;
        double l_mult_2002 = 0.7776e4 * l_prod_1987;
        double l_mult_2003 = 0.11664e5 * l_prod_1986;
        double l_mult_2004 = 0.7776e4 * l_prod_1985;
        double l_mult_2005 = 0.1944e4 * l_prod_1984;
        double l_mult_2006 = 0.3888e4 * l_prod_1983;
        double l_mult_2007 = 0.11664e5 * l_prod_1982;
        double l_mult_2008 = 0.11664e5 * l_prod_1981;
        double l_mult_2009 = 0.3888e4 * l_prod_1980;
        double l_mult_2010 = 0.3888e4 * l_prod_1979;
        double l_mult_2011 = 0.7776e4 * l_prod_1978;
        double l_mult_2012 = 0.3888e4 * l_prod_1977;
        double l_mult_2013 = 0.1944e4 * l_prod_1976;
        double l_mult_2014 = 0.1944e4 * l_prod_1975;
        double l_mult_2015 = 0.3888e3 * l_prod_1974;
        double l_sum_2016 = -0.147e2 * l_prod_523 + (0.1624e3 * l_prod_522 + (-0.6615e3 * l_prod_521 + (0.1260e4 * l_prod_976 + (-0.1134e4 * l_prod_975 + (l_mult_1995 + (0.1624e3 * l_prod_518 + (-0.1323e4 * l_prod_517 + (0.3780e4 * l_prod_972 + (-0.4536e4 * l_prod_971 + (l_mult_1996 + (-0.6615e3 * l_prod_515 + (0.3780e4 * l_prod_969 + (-0.6804e4 * l_prod_968 + (l_mult_1997 + (0.1260e4 * l_prod_966 + (-0.4536e4 * l_prod_965 + (l_mult_1998 + (-0.1134e4 * l_prod_963 + (l_mult_1999 + (l_mult_2000 + (0.1624e3 * l_prod_512 + (-0.1323e4 * l_prod_511 + (0.3780e4 * l_prod_960 + (-0.4536e4 * l_prod_959 + (l_mult_2001 + (-0.1323e4 * l_prod_509 + (0.7560e4 * l_prod_956 + (-0.13608e5 * l_prod_955 + (l_mult_2002 + (0.3780e4 * l_prod_953 + (-0.13608e5 * l_prod_952 + (l_mult_2003 + (-0.4536e4 * l_prod_950 + (l_mult_2004 + (l_mult_2005 + (-0.6615e3 * l_prod_507 + (0.3780e4 * l_prod_947 + (-0.6804e4 * l_prod_946 + (l_mult_2006 + (0.3780e4 * l_prod_945 + (-0.13608e5 * l_prod_944 + (l_mult_2007 + (-0.6804e4 * l_prod_943 + (l_mult_2008 + (l_mult_2009 + (0.1260e4 * l_prod_942 + (-0.4536e4 * l_prod_941 + (l_mult_2010 + (-0.4536e4 * l_prod_940 + (l_mult_2011 + (l_mult_2012 + (-0.1134e4 * l_prod_939 + (l_mult_2013 + (l_mult_2014 + l_mult_2015))))))))))))))))))))))))))))))))))))))))))))))))))))));
        double l_mult_2017 = -0.1e1 * l_prod_523;
        double l_mult_2018 = 0.72e1 * l_prod_522;
        double l_mult_2019 = -0.180e3 * l_prod_511;
        double l_mult_2020 = -0.99e2 * l_prod_511;
        double l_mult_2021 = 0.594e3 * l_prod_960;
        double l_mult_2022 = -0.2916e4 * l_prod_946;
        double l_mult_2023 = -0.648e3 * l_prod_941;
        double l_sum_2024 = l_mult_2023 + l_mult_2010;
        double l_mult_2025 = 0.4e1 * l_prod_522;
        double l_mult_2026 = -0.36e2 * l_prod_521;
        double l_mult_2027 = 0.72e2 * l_prod_976;
        double l_mult_2028 = -0.72e2 * l_prod_511;
        double l_mult_2029 = 0.648e3 * l_prod_960;
        double l_mult_2030 = 0.216e3 * l_prod_947;
        double l_mult_2031 = -0.1944e4 * l_prod_946;
        double l_sum_2032 = l_mult_2030 + (l_mult_2031 + l_mult_2006);
        double l_mult_2033 = -0.495e2 * l_prod_521;
        double l_mult_2034 = 0.162e3 * l_prod_976;
        double l_mult_2035 = -0.162e3 * l_prod_975;
        double l_mult_2036 = -0.1944e4 * l_prod_959;
        double l_sum_2037 = l_mult_562 + (l_mult_2021 + (l_mult_2036 + l_mult_2001));
        double l_sum_2038 = l_mult_2018 + (-0.90e2 * l_prod_521 + (0.378e3 * l_prod_976 + (-0.648e3 * l_prod_975 + l_mult_1995)));
        double l_mult_2039 = 0.72e1 * l_prod_518;
        double l_mult_2040 = -0.180e3 * l_prod_509;
        double l_mult_2041 = -0.99e2 * l_prod_509;
        double l_mult_2042 = 0.594e3 * l_prod_953;
        double l_mult_2043 = -0.2916e4 * l_prod_943;
        double l_mult_2044 = -0.648e3 * l_prod_940;
        double l_sum_2045 = l_mult_2044 + l_mult_2012;
        double l_mult_2046 = 0.4e1 * l_prod_518;
        double l_mult_2047 = -0.36e2 * l_prod_515;
        double l_mult_2048 = 0.72e2 * l_prod_966;
        double l_mult_2049 = -0.72e2 * l_prod_509;
        double l_mult_2050 = 0.648e3 * l_prod_953;
        double l_mult_2051 = 0.216e3 * l_prod_945;
        double l_mult_2052 = -0.1944e4 * l_prod_943;
        double l_sum_2053 = l_mult_2051 + (l_mult_2052 + l_mult_2009);
        double l_mult_2054 = -0.495e2 * l_prod_515;
        double l_mult_2055 = 0.162e3 * l_prod_966;
        double l_mult_2056 = -0.162e3 * l_prod_963;
        double l_mult_2057 = -0.1944e4 * l_prod_950;
        double l_sum_2058 = l_mult_566 + (l_mult_2042 + (l_mult_2057 + l_mult_2005));
        double l_sum_2059 = l_mult_2039 + (-0.90e2 * l_prod_515 + (0.378e3 * l_prod_966 + (-0.648e3 * l_prod_963 + l_mult_2000)));
        double l_mult_2060 = -0.3132e3 * l_prod_522;
        double l_mult_2061 = 0.5184e4 * l_prod_975;
        double l_mult_2062 = -0.1944e4 * l_prod_1994;
        double l_mult_2063 = 0.2088e4 * l_prod_517;
        double l_mult_2064 = -0.10044e5 * l_prod_972;
        double l_mult_2065 = -0.7776e4 * l_prod_1993;
        double l_mult_2066 = -0.5022e4 * l_prod_969;
        double l_mult_2067 = 0.15552e5 * l_prod_968;
        double l_mult_2068 = -0.11664e5 * l_prod_1992;
        double l_mult_2069 = 0.5184e4 * l_prod_965;
        double l_mult_2070 = -0.7776e4 * l_prod_1991;
        double l_mult_2071 = -0.1944e4 * l_prod_1990;
        double l_mult_2072 = 0.2088e4 * l_prod_511;
        double l_mult_2073 = -0.10044e5 * l_prod_960;
        double l_mult_2074 = -0.7776e4 * l_prod_1988;
        double l_mult_2075 = -0.10044e5 * l_prod_956;
        double l_mult_2076 = 0.31104e5 * l_prod_955;
        double l_mult_2077 = -0.23328e5 * l_prod_1987;
        double l_mult_2078 = 0.15552e5 * l_prod_952;
        double l_mult_2079 = -0.23328e5 * l_prod_1986;
        double l_mult_2080 = -0.7776e4 * l_prod_1985;
        double l_mult_2081 = -0.5022e4 * l_prod_947;
        double l_mult_2082 = 0.15552e5 * l_prod_946;
        double l_mult_2083 = -0.11664e5 * l_prod_1983;
        double l_mult_2084 = 0.15552e5 * l_prod_944;
        double l_mult_2085 = -0.23328e5 * l_prod_1982;
        double l_mult_2086 = -0.11664e5 * l_prod_1981;
        double l_mult_2087 = 0.5184e4 * l_prod_941;
        double l_mult_2088 = -0.7776e4 * l_prod_1979;
        double l_mult_2089 = -0.7776e4 * l_prod_1978;
        double l_mult_2090 = -0.1944e4 * l_prod_1976;
        double l_sum_2091 = l_mult_2060 + (0.2088e4 * l_prod_521 + (-0.5022e4 * l_prod_976 + (l_mult_2061 + (l_mult_2062 + (l_mult_2063 + (l_mult_2064 + (l_mult_1271 + (l_mult_2065 + (l_mult_2066 + (l_mult_2067 + (l_mult_2068 + (l_mult_2069 + (l_mult_2070 + (l_mult_2071 + (l_mult_2072 + (l_mult_2073 + (l_mult_1540 + (l_mult_2074 + (l_mult_2075 + (l_mult_2076 + (l_mult_2077 + (l_mult_2078 + (l_mult_2079 + (l_mult_2080 + (l_mult_2081 + (l_mult_2082 + (l_mult_2083 + (l_mult_2084 + (l_mult_2085 + (l_mult_2086 + (l_mult_2087 + (l_mult_2088 + (l_mult_2089 + l_mult_2090)))))))))))))))))))))))))))))))));
        double l_mult_2092 = 0.2565e3 * l_prod_522;
        double l_mult_2093 = 0.3888e4 * l_prod_1994;
        double l_mult_2094 = -0.1071e4 * l_prod_517;
        double l_mult_2095 = 0.9342e4 * l_prod_972;
        double l_mult_2096 = 0.11664e5 * l_prod_1993;
        double l_mult_2097 = 0.1458e4 * l_prod_969;
        double l_mult_2098 = -0.10692e5 * l_prod_968;
        double l_mult_2099 = 0.11664e5 * l_prod_1992;
        double l_mult_2100 = -0.648e3 * l_prod_965;
        double l_mult_2101 = -0.1071e4 * l_prod_511;
        double l_mult_2102 = 0.9342e4 * l_prod_960;
        double l_mult_2103 = 0.11664e5 * l_prod_1988;
        double l_mult_2104 = 0.2916e4 * l_prod_956;
        double l_mult_2105 = -0.21384e5 * l_prod_955;
        double l_mult_2106 = 0.23328e5 * l_prod_1987;
        double l_mult_2107 = -0.1944e4 * l_prod_952;
        double l_mult_2108 = 0.1458e4 * l_prod_947;
        double l_mult_2109 = -0.10692e5 * l_prod_946;
        double l_mult_2110 = 0.11664e5 * l_prod_1983;
        double l_mult_2111 = -0.1944e4 * l_prod_944;
        double l_sum_2112 = l_mult_2092 + (-0.2610e4 * l_prod_521 + (0.7884e4 * l_prod_976 + (-0.9396e4 * l_prod_975 + (l_mult_2093 + (l_mult_2094 + (l_mult_2095 + (-0.19440e5 * l_prod_971 + (l_mult_2096 + (l_mult_2097 + (l_mult_2098 + (l_mult_2099 + (l_mult_2100 + (l_mult_1998 + (l_mult_2101 + (l_mult_2102 + (-0.19440e5 * l_prod_959 + (l_mult_2103 + (l_mult_2104 + (l_mult_2105 + (l_mult_2106 + (l_mult_2107 + (l_mult_2003 + (l_mult_2108 + (l_mult_2109 + (l_mult_2110 + (l_mult_2111 + (l_mult_2007 + l_sum_2024)))))))))))))))))))))))))));
        double l_mult_2113 = -0.148e3 * l_prod_522;
        double l_mult_2114 = -0.3888e4 * l_prod_1994;
        double l_mult_2115 = 0.360e3 * l_prod_517;
        double l_mult_2116 = -0.3672e4 * l_prod_972;
        double l_mult_2117 = 0.10368e5 * l_prod_971;
        double l_mult_2118 = -0.216e3 * l_prod_969;
        double l_mult_2119 = 0.1944e4 * l_prod_968;
        double l_mult_2120 = -0.3888e4 * l_prod_1992;
        double l_mult_2121 = 0.360e3 * l_prod_511;
        double l_mult_2122 = -0.3672e4 * l_prod_960;
        double l_mult_2123 = 0.10368e5 * l_prod_959;
        double l_mult_2124 = -0.432e3 * l_prod_956;
        double l_mult_2125 = 0.3888e4 * l_prod_955;
        double l_mult_2126 = -0.7776e4 * l_prod_1987;
        double l_mult_2127 = -0.216e3 * l_prod_947;
        double l_mult_2128 = 0.1944e4 * l_prod_946;
        double l_mult_2129 = -0.3888e4 * l_prod_1983;
        double l_sum_2130 = l_mult_2127 + (l_mult_2128 + l_mult_2129);
        double l_sum_2131 = l_mult_2113 + (0.1692e4 * l_prod_521 + (-0.6120e4 * l_prod_976 + (0.8424e4 * l_prod_975 + (l_mult_2114 + (l_mult_2115 + (l_mult_2116 + (l_mult_2117 + (l_mult_2065 + (l_mult_2118 + (l_mult_2119 + (l_mult_2120 + (l_mult_2121 + (l_mult_2122 + (l_mult_2123 + (l_mult_2074 + (l_mult_2124 + (l_mult_2125 + (l_mult_2126 + l_sum_2130))))))))))))))))));
        double l_mult_2132 = 0.495e2 * l_prod_522;
        double l_mult_2133 = 0.1944e4 * l_prod_1994;
        double l_mult_2134 = 0.594e3 * l_prod_972;
        double l_mult_2135 = -0.1944e4 * l_prod_971;
        double l_sum_2136 = l_mult_2132 + (-0.5985e3 * l_prod_521 + (0.2376e4 * l_prod_976 + (-0.3726e4 * l_prod_975 + (l_mult_2133 + (l_mult_581 + (l_mult_2134 + (l_mult_2135 + (l_mult_1996 + l_sum_2037))))))));
        double l_mult_2137 = -0.72e1 * l_prod_522;
        double l_mult_2138 = 0.648e3 * l_prod_975;
        double l_mult_2139 = -0.3888e3 * l_prod_1994;
        double l_sum_2140 = l_mult_2137 + (0.90e2 * l_prod_521 + (-0.378e3 * l_prod_976 + (l_mult_2138 + l_mult_2139)));
        double l_mult_2141 = -0.3132e3 * l_prod_518;
        double l_mult_2142 = -0.5022e4 * l_prod_972;
        double l_mult_2143 = 0.5184e4 * l_prod_971;
        double l_mult_2144 = -0.1944e4 * l_prod_1993;
        double l_mult_2145 = -0.10044e5 * l_prod_969;
        double l_mult_2146 = -0.7776e4 * l_prod_1992;
        double l_mult_2147 = -0.11664e5 * l_prod_1991;
        double l_mult_2148 = 0.5184e4 * l_prod_963;
        double l_mult_2149 = -0.7776e4 * l_prod_1990;
        double l_mult_2150 = -0.1944e4 * l_prod_1989;
        double l_mult_2151 = 0.2088e4 * l_prod_509;
        double l_mult_2152 = 0.15552e5 * l_prod_955;
        double l_mult_2153 = -0.10044e5 * l_prod_953;
        double l_mult_2154 = 0.31104e5 * l_prod_952;
        double l_mult_2155 = -0.23328e5 * l_prod_1985;
        double l_mult_2156 = -0.7776e4 * l_prod_1984;
        double l_mult_2157 = -0.5022e4 * l_prod_945;
        double l_mult_2158 = -0.11664e5 * l_prod_1982;
        double l_mult_2159 = 0.15552e5 * l_prod_943;
        double l_mult_2160 = -0.23328e5 * l_prod_1981;
        double l_mult_2161 = -0.11664e5 * l_prod_1980;
        double l_mult_2162 = 0.5184e4 * l_prod_940;
        double l_mult_2163 = -0.7776e4 * l_prod_1977;
        double l_mult_2164 = -0.1944e4 * l_prod_1975;
        double l_sum_2165 = l_mult_2141 + (l_mult_2063 + (l_mult_2142 + (l_mult_2143 + (l_mult_2144 + (0.2088e4 * l_prod_515 + (l_mult_2145 + (l_mult_2067 + (l_mult_2146 + (-0.5022e4 * l_prod_966 + (l_mult_1210 + (l_mult_2147 + (l_mult_2148 + (l_mult_2149 + (l_mult_2150 + (l_mult_2151 + (l_mult_2075 + (l_mult_2152 + (l_mult_2126 + (l_mult_2153 + (l_mult_2154 + (l_mult_2079 + (l_mult_1699 + (l_mult_2155 + (l_mult_2156 + (l_mult_2157 + (l_mult_2084 + (l_mult_2158 + (l_mult_2159 + (l_mult_2160 + (l_mult_2161 + (l_mult_2162 + (l_mult_2089 + (l_mult_2163 + l_mult_2164)))))))))))))))))))))))))))))))));
        double l_mult_2166 = 0.2565e3 * l_prod_518;
        double l_mult_2167 = 0.1458e4 * l_prod_972;
        double l_mult_2168 = -0.648e3 * l_prod_971;
        double l_mult_2169 = 0.9342e4 * l_prod_969;
        double l_mult_2170 = 0.11664e5 * l_prod_1991;
        double l_mult_2171 = 0.11664e5 * l_prod_1990;
        double l_mult_2172 = 0.3888e4 * l_prod_1989;
        double l_mult_2173 = -0.1071e4 * l_prod_509;
        double l_mult_2174 = -0.1944e4 * l_prod_955;
        double l_mult_2175 = 0.9342e4 * l_prod_953;
        double l_mult_2176 = -0.21384e5 * l_prod_952;
        double l_mult_2177 = 0.23328e5 * l_prod_1985;
        double l_mult_2178 = 0.11664e5 * l_prod_1984;
        double l_mult_2179 = 0.1458e4 * l_prod_945;
        double l_mult_2180 = -0.10692e5 * l_prod_943;
        double l_mult_2181 = 0.11664e5 * l_prod_1980;
        double l_sum_2182 = l_mult_2166 + (l_mult_2094 + (l_mult_2167 + (l_mult_2168 + (-0.2610e4 * l_prod_515 + (l_mult_2169 + (l_mult_2098 + (l_mult_1997 + (0.7884e4 * l_prod_966 + (-0.19440e5 * l_prod_965 + (l_mult_2170 + (-0.9396e4 * l_prod_963 + (l_mult_2171 + (l_mult_2172 + (l_mult_2173 + (l_mult_2104 + (l_mult_2174 + (l_mult_2175 + (l_mult_2176 + (l_mult_2003 + (-0.19440e5 * l_prod_950 + (l_mult_2177 + (l_mult_2178 + (l_mult_2179 + (l_mult_2111 + (l_mult_2180 + (l_mult_2008 + (l_mult_2181 + l_sum_2045)))))))))))))))))))))))))));
        double l_mult_2183 = -0.148e3 * l_prod_518;
        double l_mult_2184 = -0.216e3 * l_prod_972;
        double l_mult_2185 = -0.3672e4 * l_prod_969;
        double l_mult_2186 = 0.10368e5 * l_prod_965;
        double l_mult_2187 = -0.3888e4 * l_prod_1991;
        double l_mult_2188 = -0.3888e4 * l_prod_1989;
        double l_mult_2189 = 0.360e3 * l_prod_509;
        double l_mult_2190 = -0.3672e4 * l_prod_953;
        double l_mult_2191 = 0.3888e4 * l_prod_952;
        double l_mult_2192 = 0.10368e5 * l_prod_950;
        double l_mult_2193 = -0.216e3 * l_prod_945;
        double l_mult_2194 = 0.1944e4 * l_prod_943;
        double l_mult_2195 = -0.3888e4 * l_prod_1980;
        double l_sum_2196 = l_mult_2193 + (l_mult_2194 + l_mult_2195);
        double l_sum_2197 = l_mult_2183 + (l_mult_2115 + (l_mult_2184 + (0.1692e4 * l_prod_515 + (l_mult_2185 + (l_mult_2119 + (-0.6120e4 * l_prod_966 + (l_mult_2186 + (l_mult_2187 + (0.8424e4 * l_prod_963 + (l_mult_2149 + (l_mult_2188 + (l_mult_2189 + (l_mult_2124 + (l_mult_2190 + (l_mult_2191 + (l_mult_2192 + (l_mult_2080 + (l_mult_2156 + l_sum_2196))))))))))))))))));
        double l_mult_2198 = 0.495e2 * l_prod_518;
        double l_mult_2199 = 0.594e3 * l_prod_969;
        double l_mult_2200 = -0.1944e4 * l_prod_965;
        double l_mult_2201 = 0.1944e4 * l_prod_1989;
        double l_sum_2202 = l_mult_2198 + (l_mult_581 + (-0.5985e3 * l_prod_515 + (l_mult_2199 + (0.2376e4 * l_prod_966 + (l_mult_2200 + (-0.3726e4 * l_prod_963 + (l_mult_1999 + (l_mult_2201 + l_sum_2058))))))));
        double l_mult_2203 = -0.72e1 * l_prod_518;
        double l_mult_2204 = 0.648e3 * l_prod_963;
        double l_mult_2205 = -0.3888e3 * l_prod_1989;
        double l_sum_2206 = l_mult_2203 + (0.90e2 * l_prod_515 + (-0.378e3 * l_prod_966 + (l_mult_2204 + l_mult_2205)));
        double l_mult_2207 = 0.36e2 * l_prod_523;
        double l_mult_2208 = 0.1044e4 * l_prod_521;
        double l_mult_2209 = -0.1674e4 * l_prod_976;
        double l_mult_2210 = 0.1296e4 * l_prod_975;
        double l_mult_2211 = 0.1044e4 * l_prod_515;
        double l_mult_2212 = -0.1674e4 * l_prod_966;
        double l_mult_2213 = 0.1296e4 * l_prod_963;
        double l_mult_2214 = 0.4176e4 * l_prod_511;
        double l_mult_2215 = -0.3888e4 * l_prod_1988;
        double l_mult_2216 = 0.4176e4 * l_prod_509;
        double l_mult_2217 = -0.20088e5 * l_prod_956;
        double l_mult_2218 = -0.15552e5 * l_prod_1987;
        double l_mult_2219 = -0.15552e5 * l_prod_1985;
        double l_mult_2220 = -0.3888e4 * l_prod_1984;
        double l_mult_2221 = -0.34992e5 * l_prod_1982;
        double l_mult_2222 = -0.34992e5 * l_prod_1981;
        double l_mult_2223 = -0.15552e5 * l_prod_1979;
        double l_mult_2224 = -0.31104e5 * l_prod_1978;
        double l_mult_2225 = -0.15552e5 * l_prod_1977;
        double l_mult_2226 = -0.9720e4 * l_prod_1976;
        double l_mult_2227 = -0.9720e4 * l_prod_1975;
        double l_mult_2228 = -0.23328e4 * l_prod_1974;
        double l_mult_2229 = -0.45e2 * l_prod_523;
        double l_mult_2230 = -0.5355e3 * l_prod_521;
        double l_mult_2231 = 0.486e3 * l_prod_976;
        double l_mult_2232 = -0.5355e3 * l_prod_515;
        double l_mult_2233 = -0.972e3 * l_prod_968;
        double l_mult_2234 = 0.486e3 * l_prod_966;
        double l_mult_2235 = -0.5220e4 * l_prod_511;
        double l_mult_2236 = -0.5220e4 * l_prod_509;
        double l_mult_2237 = 0.18684e5 * l_prod_956;
        double l_mult_2238 = -0.29160e5 * l_prod_946;
        double l_mult_2239 = 0.34992e5 * l_prod_1982;
        double l_mult_2240 = -0.29160e5 * l_prod_943;
        double l_mult_2241 = 0.34992e5 * l_prod_1981;
        double l_mult_2242 = 0.23328e5 * l_prod_1979;
        double l_mult_2243 = 0.46656e5 * l_prod_1978;
        double l_mult_2244 = 0.23328e5 * l_prod_1977;
        double l_mult_2245 = 0.19440e5 * l_prod_1976;
        double l_mult_2246 = 0.19440e5 * l_prod_1975;
        double l_mult_2247 = 0.5832e4 * l_prod_1974;
        double l_mult_2248 = 0.40e2 * l_prod_523;
        double l_mult_2249 = 0.180e3 * l_prod_521;
        double l_mult_2250 = -0.72e2 * l_prod_976;
        double l_mult_2251 = 0.180e3 * l_prod_515;
        double l_mult_2252 = -0.72e2 * l_prod_966;
        double l_mult_2253 = 0.3384e4 * l_prod_511;
        double l_mult_2254 = 0.3384e4 * l_prod_509;
        double l_mult_2255 = -0.7344e4 * l_prod_956;
        double l_mult_2256 = 0.31104e5 * l_prod_944;
        double l_mult_2257 = -0.19440e5 * l_prod_1976;
        double l_mult_2258 = -0.19440e5 * l_prod_1975;
        double l_mult_2259 = -0.225e2 * l_prod_523;
        double l_mult_2260 = -0.1197e4 * l_prod_511;
        double l_mult_2261 = -0.1197e4 * l_prod_509;
        double l_mult_2262 = 0.1188e4 * l_prod_956;
        double l_mult_2263 = -0.5832e4 * l_prod_944;
        double l_mult_2264 = 0.9720e4 * l_prod_1976;
        double l_mult_2265 = 0.9720e4 * l_prod_1975;
        double l_mult_2266 = 0.180e3 * l_prod_511;
        double l_mult_2267 = 0.180e3 * l_prod_509;
        double l_mult_2268 = 0.2592e4 * l_prod_941;
        double l_mult_2269 = 0.2592e4 * l_prod_940;
        double l_mult_2270 = 0.5184e4 * l_prod_939;
        double l_mult_2271 = -0.36e2 * l_prod_517;
        double l_mult_2272 = 0.216e3 * l_prod_969;
        double l_mult_2273 = 0.648e3 * l_prod_956;
        double l_mult_2274 = -0.3888e4 * l_prod_952;
        double l_mult_2275 = 0.324e3 * l_prod_969;
        double l_mult_2276 = 0.432e3 * l_prod_956;
        double l_sum_2277 = l_mult_2276 + (l_mult_2274 + l_mult_2004);
        double l_mult_2278 = 0.216e3 * l_prod_972;
        double l_mult_2279 = -0.3888e4 * l_prod_955;
        double l_sum_2280 = l_mult_2111 + l_mult_2007;
        double l_mult_2281 = 0.162e3 * l_prod_972;
        double l_mult_2282 = 0.162e3 * l_prod_969;
        double l_mult_2283 = 0.324e3 * l_prod_956;
        double l_sum_2284 = l_mult_2107 + l_mult_2003;
        double l_sum_2285 = l_mult_2283 + (l_mult_2174 + l_sum_2284);
        double l_mult_2286 = -0.1944e4 * l_prod_968;
        double l_sum_2287 = l_mult_2100 + l_mult_1998;
        double l_mult_2288 = 0.324e3 * l_prod_972;
        double l_sum_2289 = l_mult_2276 + (l_mult_2279 + l_mult_2002);
        double l_sum_2290 = l_mult_2272 + (l_mult_2286 + l_mult_1997);
        double l_sum_2291 = l_mult_581 + (l_mult_2134 + (l_mult_2135 + l_mult_1996));
        double l_mult_2292 = -0.3078e4 * l_prod_517;
        double l_mult_2293 = 0.12852e5 * l_prod_972;
        double l_mult_2294 = -0.17496e5 * l_prod_971;
        double l_mult_2295 = 0.7776e4 * l_prod_1993;
        double l_mult_2296 = 0.12852e5 * l_prod_969;
        double l_mult_2297 = 0.23328e5 * l_prod_1992;
        double l_mult_2298 = -0.17496e5 * l_prod_965;
        double l_mult_2299 = 0.23328e5 * l_prod_1991;
        double l_mult_2300 = 0.7776e4 * l_prod_1990;
        double l_mult_2301 = 0.12852e5 * l_prod_956;
        double l_mult_2302 = -0.34992e5 * l_prod_955;
        double l_mult_2303 = -0.34992e5 * l_prod_952;
        double l_mult_2304 = 0.46656e5 * l_prod_1986;
        double l_mult_2305 = 0.23328e5 * l_prod_1982;
        double l_mult_2306 = 0.23328e5 * l_prod_1981;
        double l_mult_2307 = 0.1332e4 * l_prod_517;
        double l_mult_2308 = -0.3240e4 * l_prod_972;
        double l_mult_2309 = 0.1944e4 * l_prod_971;
        double l_mult_2310 = -0.11232e5 * l_prod_969;
        double l_mult_2311 = 0.21384e5 * l_prod_965;
        double l_mult_2312 = -0.23328e5 * l_prod_1991;
        double l_mult_2313 = -0.11664e5 * l_prod_1990;
        double l_mult_2314 = -0.3240e4 * l_prod_956;
        double l_mult_2315 = 0.1944e4 * l_prod_944;
        double l_mult_2316 = -0.396e3 * l_prod_517;
        double l_mult_2317 = 0.432e3 * l_prod_972;
        double l_mult_2318 = 0.3996e4 * l_prod_969;
        double l_mult_2319 = -0.3888e4 * l_prod_968;
        double l_mult_2320 = -0.11016e5 * l_prod_965;
        double l_mult_2321 = 0.7776e4 * l_prod_1991;
        double l_mult_2322 = -0.594e3 * l_prod_969;
        double l_mult_2323 = 0.1944e4 * l_prod_965;
        double l_mult_2324 = -0.11232e5 * l_prod_972;
        double l_mult_2325 = 0.21384e5 * l_prod_971;
        double l_mult_2326 = -0.11664e5 * l_prod_1993;
        double l_mult_2327 = -0.3240e4 * l_prod_969;
        double l_mult_2328 = -0.23328e5 * l_prod_1992;
        double l_mult_2329 = -0.297e3 * l_prod_517;
        double l_mult_2330 = 0.2106e4 * l_prod_972;
        double l_mult_2331 = 0.2106e4 * l_prod_969;
        double l_mult_2332 = 0.36e2 * l_prod_517;
        double l_mult_2333 = -0.324e3 * l_prod_969;
        double l_mult_2334 = 0.648e3 * l_prod_965;
        double l_mult_2335 = 0.3996e4 * l_prod_972;
        double l_mult_2336 = -0.11016e5 * l_prod_971;
        double l_mult_2337 = 0.432e3 * l_prod_969;
        double l_mult_2338 = 0.7776e4 * l_prod_1992;
        double l_mult_2339 = -0.324e3 * l_prod_972;
        double l_mult_2340 = 0.648e3 * l_prod_971;
        double l_mult_2341 = -0.594e3 * l_prod_972;
        double l_mult_2342 = 0.540e3 * l_prod_522;
        double l_mult_2343 = -0.3078e4 * l_prod_521;
        double l_mult_2344 = 0.6426e4 * l_prod_976;
        double l_mult_2345 = -0.5832e4 * l_prod_975;
        double l_mult_2346 = -0.17496e5 * l_prod_968;
        double l_mult_2347 = -0.6156e4 * l_prod_511;
        double l_mult_2348 = 0.15552e5 * l_prod_1988;
        double l_mult_2349 = 0.25704e5 * l_prod_956;
        double l_mult_2350 = 0.46656e5 * l_prod_1987;
        double l_mult_2351 = 0.15552e5 * l_prod_1985;
        double l_mult_2352 = -0.52488e5 * l_prod_946;
        double l_mult_2353 = 0.34992e5 * l_prod_1983;
        double l_mult_2354 = -0.52488e5 * l_prod_944;
        double l_mult_2355 = 0.69984e5 * l_prod_1982;
        double l_mult_2356 = 0.31104e5 * l_prod_1979;
        double l_mult_2357 = 0.31104e5 * l_prod_1978;
        double l_mult_2358 = -0.360e3 * l_prod_522;
        double l_mult_2359 = 0.1332e4 * l_prod_521;
        double l_mult_2360 = -0.1620e4 * l_prod_976;
        double l_mult_2361 = -0.1620e4 * l_prod_969;
        double l_mult_2362 = 0.6984e4 * l_prod_511;
        double l_mult_2363 = -0.22464e5 * l_prod_960;
        double l_mult_2364 = -0.22464e5 * l_prod_956;
        double l_mult_2365 = 0.64152e5 * l_prod_946;
        double l_mult_2366 = -0.34992e5 * l_prod_1983;
        double l_mult_2367 = 0.64152e5 * l_prod_944;
        double l_mult_2368 = -0.69984e5 * l_prod_1982;
        double l_mult_2369 = -0.46656e5 * l_prod_1978;
        double l_mult_2370 = -0.396e3 * l_prod_521;
        double l_mult_2371 = 0.216e3 * l_prod_976;
        double l_mult_2372 = -0.4032e4 * l_prod_511;
        double l_mult_2373 = 0.7992e4 * l_prod_960;
        double l_mult_2374 = -0.3888e4 * l_prod_959;
        double l_mult_2375 = 0.7992e4 * l_prod_956;
        double l_mult_2376 = -0.7776e4 * l_prod_955;
        double l_mult_2377 = -0.33048e5 * l_prod_946;
        double l_mult_2378 = -0.33048e5 * l_prod_944;
        double l_mult_2379 = 0.54e2 * l_prod_521;
        double l_mult_2380 = -0.1188e4 * l_prod_960;
        double l_mult_2381 = -0.1188e4 * l_prod_956;
        double l_mult_2382 = 0.5832e4 * l_prod_946;
        double l_mult_2383 = 0.5832e4 * l_prod_944;
        double l_mult_2384 = 0.3492e4 * l_prod_521;
        double l_mult_2385 = -0.9612e4 * l_prod_976;
        double l_mult_2386 = 0.10368e5 * l_prod_975;
        double l_mult_2387 = 0.2664e4 * l_prod_511;
        double l_mult_2388 = -0.6480e4 * l_prod_956;
        double l_mult_2389 = -0.46656e5 * l_prod_1987;
        double l_sum_2390 = l_mult_2268 + l_mult_2223;
        double l_mult_2391 = 0.135e3 * l_prod_522;
        double l_mult_2392 = -0.1107e4 * l_prod_521;
        double l_mult_2393 = -0.972e3 * l_prod_975;
        double l_mult_2394 = -0.2214e4 * l_prod_511;
        double l_mult_2395 = 0.4212e4 * l_prod_956;
        double l_mult_2396 = -0.29160e5 * l_prod_955;
        double l_mult_2397 = -0.40824e5 * l_prod_946;
        double l_mult_2398 = -0.3888e4 * l_prod_941;
        double l_mult_2399 = 0.252e3 * l_prod_521;
        double l_mult_2400 = -0.216e3 * l_prod_976;
        double l_mult_2401 = 0.720e3 * l_prod_511;
        double l_mult_2402 = -0.4968e4 * l_prod_960;
        double l_mult_2403 = 0.3888e4 * l_prod_959;
        double l_mult_2404 = -0.648e3 * l_prod_956;
        double l_mult_2405 = 0.19440e5 * l_prod_946;
        double l_mult_2406 = -0.2016e4 * l_prod_521;
        double l_mult_2407 = 0.7020e4 * l_prod_976;
        double l_mult_2408 = -0.9072e4 * l_prod_975;
        double l_mult_2409 = -0.792e3 * l_prod_511;
        double l_mult_2410 = 0.864e3 * l_prod_956;
        double l_mult_2411 = 0.15552e5 * l_prod_1987;
        double l_mult_2412 = 0.648e3 * l_prod_947;
        double l_mult_2413 = -0.5832e4 * l_prod_946;
        double l_mult_2414 = 0.360e3 * l_prod_521;
        double l_mult_2415 = -0.972e3 * l_prod_976;
        double l_mult_2416 = 0.504e3 * l_prod_511;
        double l_mult_2417 = -0.2538e4 * l_prod_976;
        double l_mult_2418 = 0.3888e4 * l_prod_975;
        double l_mult_2419 = 0.108e3 * l_prod_511;
        double l_mult_2420 = 0.540e3 * l_prod_518;
        double l_mult_2421 = -0.3078e4 * l_prod_515;
        double l_mult_2422 = 0.6426e4 * l_prod_966;
        double l_mult_2423 = -0.5832e4 * l_prod_963;
        double l_mult_2424 = -0.6156e4 * l_prod_509;
        double l_mult_2425 = 0.46656e5 * l_prod_1985;
        double l_mult_2426 = 0.15552e5 * l_prod_1984;
        double l_mult_2427 = -0.52488e5 * l_prod_943;
        double l_mult_2428 = 0.69984e5 * l_prod_1981;
        double l_mult_2429 = 0.34992e5 * l_prod_1980;
        double l_mult_2430 = 0.31104e5 * l_prod_1977;
        double l_mult_2431 = -0.360e3 * l_prod_518;
        double l_mult_2432 = -0.1620e4 * l_prod_972;
        double l_mult_2433 = 0.1332e4 * l_prod_515;
        double l_mult_2434 = -0.1620e4 * l_prod_966;
        double l_mult_2435 = 0.6984e4 * l_prod_509;
        double l_mult_2436 = -0.22464e5 * l_prod_953;
        double l_mult_2437 = 0.64152e5 * l_prod_943;
        double l_mult_2438 = -0.69984e5 * l_prod_1981;
        double l_mult_2439 = -0.34992e5 * l_prod_1980;
        double l_mult_2440 = -0.396e3 * l_prod_515;
        double l_mult_2441 = 0.216e3 * l_prod_966;
        double l_mult_2442 = -0.4032e4 * l_prod_509;
        double l_mult_2443 = 0.7992e4 * l_prod_953;
        double l_mult_2444 = -0.7776e4 * l_prod_952;
        double l_mult_2445 = -0.3888e4 * l_prod_950;
        double l_mult_2446 = -0.33048e5 * l_prod_943;
        double l_mult_2447 = 0.54e2 * l_prod_515;
        double l_mult_2448 = -0.1188e4 * l_prod_953;
        double l_mult_2449 = 0.5832e4 * l_prod_943;
        double l_mult_2450 = 0.3492e4 * l_prod_515;
        double l_mult_2451 = -0.9612e4 * l_prod_966;
        double l_mult_2452 = 0.10368e5 * l_prod_963;
        double l_mult_2453 = 0.2664e4 * l_prod_509;
        double l_mult_2454 = -0.46656e5 * l_prod_1985;
        double l_sum_2455 = l_mult_2269 + l_mult_2225;
        double l_mult_2456 = 0.135e3 * l_prod_518;
        double l_mult_2457 = -0.1107e4 * l_prod_515;
        double l_mult_2458 = -0.972e3 * l_prod_963;
        double l_mult_2459 = -0.2214e4 * l_prod_509;
        double l_mult_2460 = -0.29160e5 * l_prod_952;
        double l_mult_2461 = -0.40824e5 * l_prod_943;
        double l_mult_2462 = -0.3888e4 * l_prod_940;
        double l_mult_2463 = 0.252e3 * l_prod_515;
        double l_mult_2464 = -0.216e3 * l_prod_966;
        double l_mult_2465 = 0.720e3 * l_prod_509;
        double l_mult_2466 = -0.4968e4 * l_prod_953;
        double l_mult_2467 = 0.3888e4 * l_prod_950;
        double l_mult_2468 = 0.19440e5 * l_prod_943;
        double l_mult_2469 = -0.2016e4 * l_prod_515;
        double l_mult_2470 = 0.7020e4 * l_prod_966;
        double l_mult_2471 = -0.9072e4 * l_prod_963;
        double l_mult_2472 = -0.792e3 * l_prod_509;
        double l_mult_2473 = 0.648e3 * l_prod_945;
        double l_mult_2474 = -0.5832e4 * l_prod_943;
        double l_mult_2475 = 0.360e3 * l_prod_515;
        double l_mult_2476 = -0.972e3 * l_prod_966;
        double l_mult_2477 = 0.504e3 * l_prod_509;
        double l_mult_2478 = -0.2538e4 * l_prod_966;
        double l_mult_2479 = 0.3888e4 * l_prod_963;
        double l_mult_2480 = 0.108e3 * l_prod_509;
        double l_mult_2481 = -0.31968e5 * l_prod_956;
        double l_mult_2482 = 0.77760e5 * l_prod_955;
        double l_mult_2483 = 0.77760e5 * l_prod_952;
        double l_mult_2484 = -0.1620e4 * l_prod_517;
        double l_mult_2485 = 0.3564e4 * l_prod_972;
        double l_mult_2486 = 0.3564e4 * l_prod_969;
        double l_mult_2487 = 0.26568e5 * l_prod_956;
        double l_mult_2488 = -0.50544e5 * l_prod_955;
        double l_mult_2489 = -0.50544e5 * l_prod_952;
        double l_mult_2490 = -0.432e3 * l_prod_972;
        double l_mult_2491 = -0.432e3 * l_prod_969;
        double l_mult_2492 = -0.8640e4 * l_prod_956;
        double l_mult_2493 = 0.7776e4 * l_prod_955;
        double l_mult_2494 = 0.7776e4 * l_prod_952;
        double l_mult_2495 = -0.25272e5 * l_prod_968;
        double l_mult_2496 = 0.7128e4 * l_prod_956;
        double l_mult_2497 = -0.2268e4 * l_prod_969;
        double l_mult_2498 = -0.4536e4 * l_prod_956;
        double l_mult_2499 = 0.3888e4 * l_prod_968;
        double l_mult_2500 = -0.864e3 * l_prod_956;
        double l_mult_2501 = -0.2268e4 * l_prod_972;
        double l_mult_2502 = -0.180e3 * l_prod_517;
        double l_mult_2503 = -0.99e2 * l_prod_517;
        double l_mult_2504 = -0.2916e4 * l_prod_968;
        double l_mult_2505 = -0.72e2 * l_prod_517;
        double l_mult_2506 = 0.648e3 * l_prod_972;
        double l_mult_2507 = 0.72e1 * l_prod_512;
        double l_sum_2508 = l_mult_2507 + (-0.90e2 * l_prod_507 + (0.378e3 * l_prod_942 + (-0.648e3 * l_prod_939 + l_mult_2015)));
        double l_mult_2509 = -0.495e2 * l_prod_507;
        double l_mult_2510 = 0.594e3 * l_prod_945;
        double l_mult_2511 = 0.162e3 * l_prod_942;
        double l_mult_2512 = -0.1944e4 * l_prod_940;
        double l_mult_2513 = -0.162e3 * l_prod_939;
        double l_sum_2514 = l_mult_2513 + l_mult_2014;
        double l_mult_2515 = 0.4e1 * l_prod_512;
        double l_mult_2516 = 0.216e3 * l_prod_953;
        double l_mult_2517 = -0.36e2 * l_prod_507;
        double l_mult_2518 = 0.72e2 * l_prod_942;
        double l_mult_2519 = -0.648e3 * l_prod_950;
        double l_sum_2520 = l_mult_584 + (l_mult_2510 + (l_mult_2043 + l_mult_2009));
        double l_mult_2521 = 0.4176e4 * l_prod_517;
        double l_mult_2522 = -0.3888e4 * l_prod_1993;
        double l_mult_2523 = -0.15552e5 * l_prod_1991;
        double l_mult_2524 = -0.9720e4 * l_prod_1990;
        double l_mult_2525 = -0.23328e4 * l_prod_1989;
        double l_mult_2526 = -0.3132e3 * l_prod_512;
        double l_mult_2527 = -0.5022e4 * l_prod_960;
        double l_mult_2528 = 0.5184e4 * l_prod_959;
        double l_mult_2529 = -0.1944e4 * l_prod_1988;
        double l_mult_2530 = -0.34992e5 * l_prod_1986;
        double l_mult_2531 = -0.31104e5 * l_prod_1985;
        double l_mult_2532 = -0.9720e4 * l_prod_1984;
        double l_mult_2533 = 0.1044e4 * l_prod_507;
        double l_mult_2534 = -0.10044e5 * l_prod_945;
        double l_mult_2535 = -0.15552e5 * l_prod_1980;
        double l_mult_2536 = -0.1674e4 * l_prod_942;
        double l_mult_2537 = -0.3888e4 * l_prod_1979;
        double l_mult_2538 = 0.10368e5 * l_prod_940;
        double l_mult_2539 = -0.15552e5 * l_prod_1978;
        double l_mult_2540 = -0.11664e5 * l_prod_1977;
        double l_mult_2541 = 0.1296e4 * l_prod_939;
        double l_mult_2542 = -0.3888e4 * l_prod_1975;
        double l_mult_2543 = -0.3888e3 * l_prod_1974;
        double l_mult_2544 = -0.5220e4 * l_prod_517;
        double l_mult_2545 = -0.29160e5 * l_prod_968;
        double l_mult_2546 = 0.19440e5 * l_prod_1990;
        double l_mult_2547 = 0.5832e4 * l_prod_1989;
        double l_mult_2548 = 0.2565e3 * l_prod_512;
        double l_mult_2549 = 0.1458e4 * l_prod_960;
        double l_mult_2550 = -0.648e3 * l_prod_959;
        double l_mult_2551 = 0.34992e5 * l_prod_1986;
        double l_mult_2552 = 0.19440e5 * l_prod_1984;
        double l_mult_2553 = -0.5355e3 * l_prod_507;
        double l_mult_2554 = -0.972e3 * l_prod_946;
        double l_mult_2555 = 0.9342e4 * l_prod_945;
        double l_mult_2556 = -0.21384e5 * l_prod_944;
        double l_mult_2557 = 0.23328e5 * l_prod_1980;
        double l_mult_2558 = 0.486e3 * l_prod_942;
        double l_mult_2559 = 0.11664e5 * l_prod_1977;
        double l_mult_2560 = 0.3384e4 * l_prod_517;
        double l_mult_2561 = -0.19440e5 * l_prod_1990;
        double l_mult_2562 = -0.148e3 * l_prod_512;
        double l_mult_2563 = -0.216e3 * l_prod_960;
        double l_mult_2564 = -0.11664e5 * l_prod_1986;
        double l_mult_2565 = -0.19440e5 * l_prod_1984;
        double l_mult_2566 = 0.180e3 * l_prod_507;
        double l_mult_2567 = -0.3672e4 * l_prod_945;
        double l_mult_2568 = 0.3888e4 * l_prod_944;
        double l_mult_2569 = -0.72e2 * l_prod_942;
        double l_mult_2570 = -0.3888e4 * l_prod_1977;
        double l_mult_2571 = -0.1197e4 * l_prod_517;
        double l_mult_2572 = 0.9720e4 * l_prod_1990;
        double l_mult_2573 = 0.495e2 * l_prod_512;
        double l_mult_2574 = -0.5832e4 * l_prod_952;
        double l_mult_2575 = 0.9720e4 * l_prod_1984;
        double l_mult_2576 = 0.180e3 * l_prod_517;
        double l_mult_2577 = 0.2592e4 * l_prod_965;
        double l_mult_2578 = -0.72e1 * l_prod_512;
        double l_mult_2579 = 0.2592e4 * l_prod_950;
        double l_mult_2580 = -0.1944e4 * l_prod_1984;
        double l_mult_2581 = -0.5022e4 * l_prod_953;
        double l_mult_2582 = 0.5184e4 * l_prod_950;
        double l_mult_2583 = -0.10044e5 * l_prod_947;
        double l_mult_2584 = -0.7776e4 * l_prod_1983;
        double l_mult_2585 = -0.7776e4 * l_prod_1980;
        double l_mult_2586 = -0.11664e5 * l_prod_1979;
        double l_mult_2587 = -0.23328e5 * l_prod_1978;
        double l_mult_2588 = -0.7776e4 * l_prod_1976;
        double l_mult_2589 = -0.7776e4 * l_prod_1975;
        double l_mult_2590 = -0.1944e4 * l_prod_1974;
        double l_sum_2591 = l_mult_2526 + (l_mult_2072 + (l_mult_2527 + (l_mult_2528 + (l_mult_2529 + (l_mult_2151 + (l_mult_2075 + (l_mult_2152 + (l_mult_2126 + (l_mult_2581 + (l_mult_2078 + (l_mult_2564 + (l_mult_2582 + (l_mult_2080 + (l_mult_2580 + (0.2088e4 * l_prod_507 + (l_mult_2583 + (l_mult_2082 + (l_mult_2584 + (l_mult_2534 + (l_mult_2256 + (l_mult_2085 + (l_mult_2159 + (l_mult_2160 + (l_mult_2585 + (-0.5022e4 * l_prod_942 + (l_mult_1813 + (l_mult_2586 + (l_mult_1870 + (l_mult_2587 + (l_mult_2540 + (l_mult_2270 + (l_mult_2588 + (l_mult_2589 + l_mult_2590)))))))))))))))))))))))))))))))));
        double l_mult_2592 = 0.1458e4 * l_prod_953;
        double l_mult_2593 = 0.9342e4 * l_prod_947;
        double l_mult_2594 = 0.11664e5 * l_prod_1979;
        double l_mult_2595 = 0.23328e5 * l_prod_1978;
        double l_mult_2596 = 0.11664e5 * l_prod_1976;
        double l_mult_2597 = 0.11664e5 * l_prod_1975;
        double l_mult_2598 = 0.3888e4 * l_prod_1974;
        double l_sum_2599 = l_mult_2548 + (l_mult_2101 + (l_mult_2549 + (l_mult_2550 + (l_mult_2173 + (l_mult_2104 + (l_mult_2174 + (l_mult_2592 + (l_mult_2107 + (l_mult_2519 + (-0.2610e4 * l_prod_507 + (l_mult_2593 + (l_mult_2109 + (l_mult_2006 + (l_mult_2555 + (l_mult_2556 + (l_mult_2007 + (l_mult_2180 + (l_mult_2008 + (l_mult_2009 + (0.7884e4 * l_prod_942 + (-0.19440e5 * l_prod_941 + (l_mult_2594 + (-0.19440e5 * l_prod_940 + (l_mult_2595 + (l_mult_2559 + (-0.9396e4 * l_prod_939 + (l_mult_2596 + (l_mult_2597 + l_mult_2598))))))))))))))))))))))))))));
        double l_mult_2600 = -0.216e3 * l_prod_953;
        double l_mult_2601 = -0.3672e4 * l_prod_947;
        double l_mult_2602 = 0.10368e5 * l_prod_941;
        double l_mult_2603 = -0.3888e4 * l_prod_1974;
        double l_sum_2604 = l_mult_2562 + (l_mult_2121 + (l_mult_2563 + (l_mult_2189 + (l_mult_2124 + (l_mult_2600 + (0.1692e4 * l_prod_507 + (l_mult_2601 + (l_mult_2128 + (l_mult_2567 + (l_mult_2568 + (l_mult_2194 + (-0.6120e4 * l_prod_942 + (l_mult_2602 + (l_mult_2537 + (l_mult_2538 + (l_mult_2089 + (l_mult_2570 + (0.8424e4 * l_prod_939 + (l_mult_2588 + (l_mult_2589 + l_mult_2603))))))))))))))))))));
        double l_mult_2605 = 0.594e3 * l_prod_947;
        double l_mult_2606 = -0.1944e4 * l_prod_941;
        double l_mult_2607 = 0.1944e4 * l_prod_1974;
        double l_sum_2608 = l_mult_2573 + (l_mult_562 + (l_mult_566 + (-0.5985e3 * l_prod_507 + (l_mult_2605 + (l_mult_2510 + (0.2376e4 * l_prod_942 + (l_mult_2606 + (l_mult_2512 + (-0.3726e4 * l_prod_939 + (l_mult_2013 + (l_mult_2014 + l_mult_2607)))))))))));
        double l_mult_2609 = 0.648e3 * l_prod_939;
        double l_sum_2610 = l_mult_2578 + (0.90e2 * l_prod_507 + (-0.378e3 * l_prod_942 + (l_mult_2609 + l_mult_2543)));
        double l_mult_2611 = -0.36e2 * l_prod_511;
        double l_mult_2612 = 0.324e3 * l_prod_947;
        double l_mult_2613 = -0.3888e4 * l_prod_944;
        double l_sum_2614 = l_mult_2030 + (l_mult_2613 + l_mult_2008);
        double l_mult_2615 = 0.216e3 * l_prod_960;
        double l_mult_2616 = 0.162e3 * l_prod_960;
        double l_sum_2617 = 0.162e3 * l_prod_947 + (l_mult_2554 + l_sum_2280);
        double l_mult_2618 = 0.324e3 * l_prod_960;
        double l_mult_2619 = -0.6156e4 * l_prod_517;
        double l_mult_2620 = 0.15552e5 * l_prod_1993;
        double l_mult_2621 = -0.52488e5 * l_prod_968;
        double l_mult_2622 = 0.34992e5 * l_prod_1992;
        double l_mult_2623 = 0.31104e5 * l_prod_1991;
        double l_mult_2624 = -0.3078e4 * l_prod_511;
        double l_mult_2625 = 0.12852e5 * l_prod_960;
        double l_mult_2626 = -0.17496e5 * l_prod_959;
        double l_mult_2627 = 0.7776e4 * l_prod_1988;
        double l_mult_2628 = -0.52488e5 * l_prod_952;
        double l_mult_2629 = 0.69984e5 * l_prod_1986;
        double l_mult_2630 = 0.31104e5 * l_prod_1985;
        double l_mult_2631 = -0.17496e5 * l_prod_946;
        double l_mult_2632 = -0.34992e5 * l_prod_944;
        double l_mult_2633 = 0.46656e5 * l_prod_1982;
        double l_mult_2634 = 0.7776e4 * l_prod_1979;
        double l_mult_2635 = 0.15552e5 * l_prod_1978;
        double l_mult_2636 = 0.6984e4 * l_prod_517;
        double l_mult_2637 = -0.22464e5 * l_prod_972;
        double l_mult_2638 = 0.64152e5 * l_prod_968;
        double l_mult_2639 = -0.34992e5 * l_prod_1992;
        double l_mult_2640 = 0.1332e4 * l_prod_511;
        double l_mult_2641 = -0.3240e4 * l_prod_960;
        double l_mult_2642 = 0.1944e4 * l_prod_959;
        double l_mult_2643 = 0.64152e5 * l_prod_952;
        double l_mult_2644 = -0.69984e5 * l_prod_1986;
        double l_mult_2645 = -0.1620e4 * l_prod_947;
        double l_mult_2646 = 0.648e3 * l_prod_941;
        double l_mult_2647 = -0.4032e4 * l_prod_517;
        double l_mult_2648 = 0.7992e4 * l_prod_972;
        double l_mult_2649 = -0.3888e4 * l_prod_971;
        double l_mult_2650 = -0.33048e5 * l_prod_968;
        double l_mult_2651 = -0.396e3 * l_prod_511;
        double l_mult_2652 = 0.432e3 * l_prod_960;
        double l_mult_2653 = -0.33048e5 * l_prod_952;
        double l_mult_2654 = 0.23328e5 * l_prod_1986;
        double l_mult_2655 = -0.1188e4 * l_prod_972;
        double l_mult_2656 = 0.5832e4 * l_prod_968;
        double l_mult_2657 = 0.5832e4 * l_prod_952;
        double l_mult_2658 = 0.2664e4 * l_prod_517;
        double l_mult_2659 = -0.11232e5 * l_prod_960;
        double l_mult_2660 = 0.21384e5 * l_prod_959;
        double l_mult_2661 = -0.11664e5 * l_prod_1988;
        double l_sum_2662 = l_mult_2646 + l_mult_2537;
        double l_mult_2663 = -0.2214e4 * l_prod_517;
        double l_mult_2664 = -0.40824e5 * l_prod_968;
        double l_mult_2665 = -0.3888e4 * l_prod_965;
        double l_mult_2666 = -0.297e3 * l_prod_511;
        double l_mult_2667 = 0.2106e4 * l_prod_960;
        double l_mult_2668 = 0.720e3 * l_prod_517;
        double l_mult_2669 = -0.4968e4 * l_prod_972;
        double l_mult_2670 = 0.3888e4 * l_prod_971;
        double l_mult_2671 = 0.19440e5 * l_prod_968;
        double l_mult_2672 = 0.36e2 * l_prod_511;
        double l_mult_2673 = 0.1944e4 * l_prod_952;
        double l_mult_2674 = -0.792e3 * l_prod_517;
        double l_mult_2675 = 0.648e3 * l_prod_969;
        double l_mult_2676 = -0.5832e4 * l_prod_968;
        double l_mult_2677 = 0.3996e4 * l_prod_960;
        double l_mult_2678 = -0.11016e5 * l_prod_959;
        double l_mult_2679 = 0.504e3 * l_prod_517;
        double l_mult_2680 = -0.324e3 * l_prod_960;
        double l_mult_2681 = 0.648e3 * l_prod_959;
        double l_mult_2682 = 0.108e3 * l_prod_517;
        double l_sum_2683 = l_mult_554 + (-0.594e3 * l_prod_960 + (l_mult_2642 + l_mult_2529));
        double l_mult_2684 = 0.12852e5 * l_prod_947;
        double l_mult_2685 = 0.23328e5 * l_prod_1983;
        double l_mult_2686 = -0.17496e5 * l_prod_941;
        double l_mult_2687 = 0.7776e4 * l_prod_1976;
        double l_mult_2688 = -0.11232e5 * l_prod_947;
        double l_mult_2689 = 0.21384e5 * l_prod_941;
        double l_mult_2690 = -0.23328e5 * l_prod_1979;
        double l_mult_2691 = -0.11664e5 * l_prod_1976;
        double l_mult_2692 = 0.3996e4 * l_prod_947;
        double l_mult_2693 = -0.3888e4 * l_prod_946;
        double l_mult_2694 = -0.11016e5 * l_prod_941;
        double l_mult_2695 = -0.594e3 * l_prod_947;
        double l_mult_2696 = 0.1944e4 * l_prod_941;
        double l_mult_2697 = -0.3240e4 * l_prod_947;
        double l_mult_2698 = -0.23328e5 * l_prod_1983;
        double l_sum_2699 = l_mult_2568 + (l_mult_2085 + (l_mult_2696 + l_mult_2586));
        double l_mult_2700 = 0.2106e4 * l_prod_947;
        double l_sum_2701 = l_mult_2606 + l_mult_2594;
        double l_mult_2702 = -0.324e3 * l_prod_947;
        double l_mult_2703 = 0.432e3 * l_prod_947;
        double l_mult_2704 = 0.7776e4 * l_prod_1983;
        double l_mult_2705 = 0.540e3 * l_prod_512;
        double l_mult_2706 = -0.3078e4 * l_prod_507;
        double l_mult_2707 = 0.31104e5 * l_prod_1980;
        double l_mult_2708 = 0.6426e4 * l_prod_942;
        double l_mult_2709 = 0.34992e5 * l_prod_1977;
        double l_mult_2710 = -0.5832e4 * l_prod_939;
        double l_mult_2711 = 0.15552e5 * l_prod_1975;
        double l_mult_2712 = -0.360e3 * l_prod_512;
        double l_mult_2713 = -0.1620e4 * l_prod_960;
        double l_mult_2714 = 0.3492e4 * l_prod_507;
        double l_mult_2715 = -0.22464e5 * l_prod_945;
        double l_mult_2716 = -0.9612e4 * l_prod_942;
        double l_mult_2717 = -0.34992e5 * l_prod_1977;
        double l_mult_2718 = 0.10368e5 * l_prod_939;
        double l_mult_2719 = -0.2016e4 * l_prod_507;
        double l_mult_2720 = 0.7992e4 * l_prod_945;
        double l_mult_2721 = -0.7776e4 * l_prod_944;
        double l_mult_2722 = 0.7020e4 * l_prod_942;
        double l_mult_2723 = -0.9072e4 * l_prod_939;
        double l_mult_2724 = -0.1188e4 * l_prod_945;
        double l_mult_2725 = -0.2538e4 * l_prod_942;
        double l_mult_2726 = 0.3888e4 * l_prod_940;
        double l_mult_2727 = 0.3888e4 * l_prod_939;
        double l_mult_2728 = 0.1332e4 * l_prod_507;
        double l_mult_2729 = -0.1620e4 * l_prod_942;
        double l_sum_2730 = l_mult_2609 + l_mult_2589;
        double l_mult_2731 = 0.135e3 * l_prod_512;
        double l_mult_2732 = -0.1107e4 * l_prod_507;
        double l_mult_2733 = -0.29160e5 * l_prod_944;
        double l_mult_2734 = -0.972e3 * l_prod_939;
        double l_mult_2735 = 0.360e3 * l_prod_507;
        double l_mult_2736 = -0.4968e4 * l_prod_945;
        double l_mult_2737 = -0.972e3 * l_prod_942;
        double l_mult_2738 = -0.396e3 * l_prod_507;
        double l_mult_2739 = 0.216e3 * l_prod_942;
        double l_mult_2740 = 0.252e3 * l_prod_507;
        double l_mult_2741 = -0.216e3 * l_prod_942;
        double l_mult_2742 = 0.54e2 * l_prod_507;
        double l_mult_2743 = 0.77760e5 * l_prod_944;
        double l_mult_2744 = -0.1620e4 * l_prod_511;
        double l_mult_2745 = 0.3564e4 * l_prod_960;
        double l_mult_2746 = -0.25272e5 * l_prod_946;
        double l_mult_2747 = -0.50544e5 * l_prod_944;
        double l_mult_2748 = -0.432e3 * l_prod_960;
        double l_mult_2749 = 0.3888e4 * l_prod_946;
        double l_mult_2750 = 0.7776e4 * l_prod_944;
        double l_mult_2751 = 0.3564e4 * l_prod_947;
        double l_mult_2752 = -0.2268e4 * l_prod_947;
        double l_mult_2753 = -0.432e3 * l_prod_947;
        double l_mult_2754 = -0.2268e4 * l_prod_960;
        double l_sum_2755 = l_mult_2513 + l_mult_2013;
        double l_sum_2756 = l_mult_584 + (l_mult_2605 + (l_mult_2022 + l_mult_2006));
        double l_mult_2757 = -0.23328e4 * l_prod_1994;
        double l_mult_2758 = -0.9720e4 * l_prod_1993;
        double l_mult_2759 = -0.15552e5 * l_prod_1992;
        double l_mult_2760 = -0.3888e4 * l_prod_1990;
        double l_mult_2761 = -0.9720e4 * l_prod_1988;
        double l_mult_2762 = -0.31104e5 * l_prod_1987;
        double l_mult_2763 = -0.15552e5 * l_prod_1983;
        double l_mult_2764 = -0.3888e4 * l_prod_1976;
        double l_mult_2765 = 0.5832e4 * l_prod_1994;
        double l_mult_2766 = 0.19440e5 * l_prod_1993;
        double l_mult_2767 = 0.19440e5 * l_prod_1988;
        double l_mult_2768 = -0.972e3 * l_prod_943;
        double l_mult_2769 = -0.19440e5 * l_prod_1993;
        double l_mult_2770 = -0.19440e5 * l_prod_1988;
        double l_mult_2771 = 0.9720e4 * l_prod_1993;
        double l_mult_2772 = 0.9720e4 * l_prod_1988;
        double l_mult_2773 = -0.5832e4 * l_prod_955;
        double l_mult_2774 = 0.2592e4 * l_prod_971;
        double l_mult_2775 = 0.2592e4 * l_prod_959;
        double l_mult_2776 = -0.36e2 * l_prod_509;
        double l_mult_2777 = 0.324e3 * l_prod_945;
        double l_mult_2778 = 0.324e3 * l_prod_953;
        double l_mult_2779 = 0.162e3 * l_prod_953;
        double l_sum_2780 = 0.162e3 * l_prod_945 + (l_mult_2111 + (l_mult_2768 + l_mult_2008));
        double l_sum_2781 = l_mult_2051 + (l_mult_2613 + l_mult_2007);
        double l_mult_2782 = 0.31104e5 * l_prod_1992;
        double l_mult_2783 = 0.34992e5 * l_prod_1991;
        double l_mult_2784 = 0.15552e5 * l_prod_1990;
        double l_mult_2785 = -0.3078e4 * l_prod_509;
        double l_mult_2786 = -0.52488e5 * l_prod_955;
        double l_mult_2787 = 0.31104e5 * l_prod_1987;
        double l_mult_2788 = 0.12852e5 * l_prod_953;
        double l_mult_2789 = -0.17496e5 * l_prod_950;
        double l_mult_2790 = 0.7776e4 * l_prod_1984;
        double l_mult_2791 = -0.17496e5 * l_prod_943;
        double l_mult_2792 = 0.46656e5 * l_prod_1981;
        double l_mult_2793 = 0.7776e4 * l_prod_1977;
        double l_mult_2794 = -0.22464e5 * l_prod_969;
        double l_mult_2795 = -0.34992e5 * l_prod_1991;
        double l_mult_2796 = 0.1332e4 * l_prod_509;
        double l_mult_2797 = 0.5832e4 * l_prod_955;
        double l_mult_2798 = -0.11232e5 * l_prod_953;
        double l_mult_2799 = 0.21384e5 * l_prod_950;
        double l_mult_2800 = -0.11664e5 * l_prod_1984;
        double l_mult_2801 = -0.1620e4 * l_prod_945;
        double l_mult_2802 = 0.648e3 * l_prod_940;
        double l_sum_2803 = l_mult_2802 + l_mult_2570;
        double l_mult_2804 = 0.7992e4 * l_prod_969;
        double l_mult_2805 = -0.396e3 * l_prod_509;
        double l_mult_2806 = 0.3996e4 * l_prod_953;
        double l_mult_2807 = -0.11016e5 * l_prod_950;
        double l_mult_2808 = -0.1188e4 * l_prod_969;
        double l_mult_2809 = 0.3888e4 * l_prod_965;
        double l_mult_2810 = 0.1944e4 * l_prod_950;
        double l_sum_2811 = l_mult_555 + (-0.594e3 * l_prod_953 + (l_mult_2810 + l_mult_2580));
        double l_mult_2812 = 0.64152e5 * l_prod_955;
        double l_mult_2813 = -0.3240e4 * l_prod_953;
        double l_mult_2814 = -0.297e3 * l_prod_509;
        double l_mult_2815 = 0.2106e4 * l_prod_953;
        double l_mult_2816 = -0.4968e4 * l_prod_969;
        double l_mult_2817 = 0.36e2 * l_prod_509;
        double l_mult_2818 = -0.324e3 * l_prod_953;
        double l_mult_2819 = 0.648e3 * l_prod_950;
        double l_mult_2820 = -0.33048e5 * l_prod_955;
        double l_mult_2821 = 0.432e3 * l_prod_953;
        double l_mult_2822 = 0.1944e4 * l_prod_955;
        double l_mult_2823 = 0.31104e5 * l_prod_1983;
        double l_mult_2824 = 0.12852e5 * l_prod_945;
        double l_mult_2825 = 0.7776e4 * l_prod_1980;
        double l_mult_2826 = 0.34992e5 * l_prod_1979;
        double l_mult_2827 = -0.17496e5 * l_prod_940;
        double l_mult_2828 = 0.15552e5 * l_prod_1976;
        double l_mult_2829 = 0.7776e4 * l_prod_1975;
        double l_mult_2830 = -0.1620e4 * l_prod_953;
        double l_mult_2831 = -0.22464e5 * l_prod_947;
        double l_mult_2832 = -0.11232e5 * l_prod_945;
        double l_mult_2833 = -0.34992e5 * l_prod_1979;
        double l_mult_2834 = 0.21384e5 * l_prod_940;
        double l_mult_2835 = -0.11664e5 * l_prod_1975;
        double l_mult_2836 = 0.7992e4 * l_prod_947;
        double l_mult_2837 = 0.3996e4 * l_prod_945;
        double l_mult_2838 = -0.11016e5 * l_prod_940;
        double l_mult_2839 = -0.1188e4 * l_prod_947;
        double l_mult_2840 = -0.594e3 * l_prod_945;
        double l_mult_2841 = 0.3888e4 * l_prod_941;
        double l_mult_2842 = 0.1944e4 * l_prod_940;
        double l_mult_2843 = -0.3240e4 * l_prod_945;
        double l_sum_2844 = l_mult_2609 + l_mult_2588;
        double l_mult_2845 = 0.2106e4 * l_prod_945;
        double l_mult_2846 = -0.4968e4 * l_prod_947;
        double l_mult_2847 = -0.324e3 * l_prod_945;
        double l_mult_2848 = 0.432e3 * l_prod_945;
        double l_mult_2849 = -0.23328e5 * l_prod_1977;
        double l_mult_2850 = -0.3888e4 * l_prod_943;
        double l_mult_2851 = -0.23328e5 * l_prod_1980;
        double l_sum_2852 = l_mult_2842 + l_mult_2540;
        double l_sum_2853 = l_mult_2512 + l_mult_2559;
        double l_mult_2854 = -0.1620e4 * l_prod_509;
        double l_mult_2855 = 0.3564e4 * l_prod_953;
        double l_mult_2856 = -0.25272e5 * l_prod_943;
        double l_mult_2857 = -0.432e3 * l_prod_953;
        double l_mult_2858 = 0.3888e4 * l_prod_943;
        double l_mult_2859 = 0.3564e4 * l_prod_945;
        double l_mult_2860 = -0.2268e4 * l_prod_953;
        double l_mult_2861 = -0.2268e4 * l_prod_945;
        double l_mult_2862 = -0.432e3 * l_prod_945;
        double l_vdot_2863 = vdot4(v_917, vcons4(0.e0, 0.e0, 0.e0, 0.e0)) + (vdot4(v_918,
            vcons4(0.e0,
                l_mult_2018 + (l_mult_2019 + (0.1134e4 * l_prod_947 + (-0.2592e4 * l_prod_941 + l_mult_2013))),
                l_mult_547 + (l_mult_561 + (l_mult_2020 + (l_mult_2021 + (0.486e3 * l_prod_947 + (l_mult_2022 + l_sum_2024))))),
                l_mult_2025 + (l_mult_2026 + (l_mult_2027 + (l_mult_2028 + (l_mult_2029 + (-0.1296e4 * l_prod_959 + l_sum_2032))))))) + (vdot4(
            v_919,
            vcons4(l_mult_547 + (l_mult_2033 + (l_mult_2034 + (l_mult_2035 + l_sum_2037))), l_sum_2038,
                l_mult_2039 + (l_mult_2040 + (0.1134e4 * l_prod_945 + (-0.2592e4 * l_prod_940 + l_mult_2014))),
                l_mult_551 + (l_mult_565 + (l_mult_2041 + (l_mult_2042 + (0.486e3 * l_prod_945 + (l_mult_2043 + l_sum_2045))))))) + (vdot4(
            v_920,
            vcons4(
                l_mult_2046 + (l_mult_2047 + (l_mult_2048 + (l_mult_2049 + (l_mult_2050 + (-0.1296e4 * l_prod_950 + l_sum_2053))))),
                l_mult_551 + (l_mult_2054 + (l_mult_2055 + (l_mult_2056 + l_sum_2058))), l_sum_2059, l_sum_2091)) + (vdot4(
            v_921, vcons4(l_sum_2112, l_sum_2131, l_sum_2136, l_sum_2140)) + (vdot4(v_922,
            vcons4(l_sum_2165, l_sum_2182, l_sum_2197, l_sum_2202)) + (vdot4(v_923,
            vcons4(l_sum_2206,
                l_mult_2207 + (l_mult_2060 + (l_mult_2208 + (l_mult_2209 + (l_mult_2210 + (l_mult_2139 + (l_mult_2141 + (l_mult_2063 + (l_mult_2142 + (l_mult_2143 + (l_mult_2144 + (l_mult_2211 + (l_mult_2066 + (0.7776e4 * l_prod_968 + (l_mult_2120 + (l_mult_2212 + (l_mult_2069 + (l_mult_2187 + (l_mult_2213 + (l_mult_2071 + (l_mult_2205 + (-0.6264e3 * l_prod_512 + (l_mult_2214 + (l_mult_2073 + (l_mult_2123 + (l_mult_2215 + (l_mult_2216 + (l_mult_2217 + (l_mult_2076 + (l_mult_2218 + (l_mult_2153 + (l_mult_2154 + (l_mult_2079 + (l_mult_2192 + (l_mult_2219 + (l_mult_2220 + (0.3132e4 * l_prod_507 + (-0.15066e5 * l_prod_947 + (l_mult_1489 + (l_mult_2083 + (-0.15066e5 * l_prod_945 + (l_mult_1490 + (l_mult_2221 + (l_mult_1671 + (l_mult_2222 + (l_mult_2161 + (-0.6696e4 * l_prod_942 + (0.20736e5 * l_prod_941 + (l_mult_2223 + (0.20736e5 * l_prod_940 + (l_mult_2224 + (l_mult_2225 + (0.6480e4 * l_prod_939 + (l_mult_2226 + (l_mult_2227 + l_mult_2228)))))))))))))))))))))))))))))))))))))))))))))))))))))),
                l_mult_2229 + (l_mult_2092 + (l_mult_2230 + (l_mult_2231 + (l_mult_2035 + (l_mult_2166 + (l_mult_2094 + (l_mult_2167 + (l_mult_2168 + (l_mult_2232 + (l_mult_2097 + (l_mult_2233 + (l_mult_2234 + (l_mult_2100 + (l_mult_2056 + (0.1053e4 * l_prod_512 + (l_mult_2235 + (l_mult_2102 + (-0.7128e4 * l_prod_959 + (l_mult_2001 + (l_mult_2236 + (l_mult_2237 + (l_mult_2105 + (l_mult_2002 + (l_mult_2175 + (l_mult_2176 + (l_mult_2003 + (-0.7128e4 * l_prod_950 + (l_mult_2004 + (l_mult_2005 + (-0.62235e4 * l_prod_507 + (0.23652e5 * l_prod_947 + (l_mult_2238 + (l_mult_2110 + (0.23652e5 * l_prod_945 + (-0.58320e5 * l_prod_944 + (l_mult_2239 + (l_mult_2240 + (l_mult_2241 + (l_mult_2181 + (0.14796e5 * l_prod_942 + (-0.37584e5 * l_prod_941 + (l_mult_2242 + (-0.37584e5 * l_prod_940 + (l_mult_2243 + (l_mult_2244 + (-0.15390e5 * l_prod_939 + (l_mult_2245 + (l_mult_2246 + l_mult_2247)))))))))))))))))))))))))))))))))))))))))))))))),
                l_mult_2248 + (l_mult_2113 + (l_mult_2249 + (l_mult_2250 + (l_mult_2183 + (l_mult_2115 + (l_mult_2184 + (l_mult_2251 + (l_mult_2118 + (l_mult_2252 + (-0.1016e4 * l_prod_512 + (l_mult_2253 + (l_mult_2122 + (0.1296e4 * l_prod_959 + (l_mult_2254 + (l_mult_2255 + (l_mult_2125 + (l_mult_2190 + (l_mult_2191 + (0.1296e4 * l_prod_950 + (0.6696e4 * l_prod_507 + (-0.18360e5 * l_prod_947 + (l_mult_2082 + (l_mult_2129 + (-0.18360e5 * l_prod_945 + (l_mult_2256 + (l_mult_2158 + (l_mult_2159 + (l_mult_2086 + (l_mult_2195 + (-0.17424e5 * l_prod_942 + (0.33696e5 * l_prod_941 + (l_mult_2223 + (0.33696e5 * l_prod_940 + (l_mult_2224 + (l_mult_2225 + (l_mult_1452 + (l_mult_2257 + (l_mult_2258 + -0.7776e4 * l_prod_1974)))))))))))))))))))))))))))))))))))))))) + (vdot4(
            v_924,
            vcons4(
                l_mult_2259 + (l_mult_2132 + (l_mult_561 + (l_mult_2198 + (l_mult_581 + (l_mult_565 + (0.594e3 * l_prod_512 + (l_mult_2260 + (l_mult_2021 + (l_mult_2261 + (l_mult_2262 + (l_mult_2042 + (-0.41445e4 * l_prod_507 + (0.7128e4 * l_prod_947 + (l_mult_2022 + (0.7128e4 * l_prod_945 + (l_mult_2263 + (l_mult_2043 + (0.11556e5 * l_prod_942 + (-0.14904e5 * l_prod_941 + (l_mult_2010 + (-0.14904e5 * l_prod_940 + (l_mult_2011 + (l_mult_2012 + (-0.13770e5 * l_prod_939 + (l_mult_2264 + (l_mult_2265 + l_mult_2247)))))))))))))))))))))))))),
                l_mult_1354 + (l_mult_2137 + (l_mult_2203 + (-0.1944e3 * l_prod_512 + (l_mult_2266 + (l_mult_2267 + (0.1404e4 * l_prod_507 + (-0.1134e4 * l_prod_947 + (-0.1134e4 * l_prod_945 + (-0.4104e4 * l_prod_942 + (l_mult_2268 + (l_mult_2269 + (l_mult_2270 + (l_mult_2090 + (l_mult_2164 + l_mult_2228)))))))))))))),
                l_mult_581 + (l_mult_2262 + (l_mult_2263 + l_mult_2011)),
                l_mult_2271 + (l_mult_2272 + (l_mult_2273 + (l_mult_2274 + (l_mult_2111 + l_mult_2008)))))) + (vdot4(
            v_925,
            vcons4(l_mult_2271 + (l_mult_2275 + (l_mult_2100 + l_sum_2277)),
                l_mult_581 + (l_mult_2199 + (l_mult_2200 + l_mult_1999)),
                l_mult_2271 + (l_mult_2278 + (l_mult_2273 + (l_mult_2279 + l_sum_2280))),
                l_mult_525 + (l_mult_2281 + (l_mult_2282 + (l_mult_2233 + l_sum_2285))))) + (vdot4(v_926,
            vcons4(l_mult_2271 + (l_mult_2278 + (l_mult_2275 + (l_mult_2286 + l_sum_2287))),
                l_mult_2271 + (l_mult_2288 + (l_mult_2168 + l_sum_2289)),
                l_mult_2271 + (l_mult_2288 + (l_mult_2168 + l_sum_2290)), l_sum_2291)) + (vdot4(v_927,
            vcons4(
                l_mult_2292 + (l_mult_2293 + (l_mult_2294 + (l_mult_2295 + (l_mult_2296 + (l_mult_1379 + (l_mult_2297 + (l_mult_2298 + (l_mult_2299 + (l_mult_2300 + (l_mult_2301 + (l_mult_2302 + (l_mult_2106 + (l_mult_2303 + (l_mult_2304 + (l_mult_2177 + (-0.17496e5 * l_prod_944 + (l_mult_2305 + (l_mult_2306 + l_mult_2011)))))))))))))))))),
                l_mult_2307 + (l_mult_2308 + (l_mult_2309 + (l_mult_2310 + (l_mult_1502 + (l_mult_2068 + (l_mult_2311 + (l_mult_2312 + (l_mult_2313 + (l_mult_2314 + (l_mult_2125 + (l_mult_984 + (l_mult_2079 + (l_mult_2155 + (l_mult_2315 + l_mult_2086)))))))))))))),
                l_mult_2316 + (l_mult_2317 + (l_mult_2318 + (l_mult_2319 + (l_mult_2320 + (l_mult_2321 + (l_mult_2300 + l_sum_2277)))))),
                l_mult_573 + (l_mult_2322 + (l_mult_2323 + l_mult_2071)))) + (vdot4(v_928,
            vcons4(
                l_mult_2307 + (l_mult_2324 + (l_mult_2325 + (l_mult_2326 + (l_mult_2327 + (l_mult_1502 + (l_mult_2328 + (l_mult_2323 + (l_mult_2147 + (l_mult_2314 + (l_mult_983 + (l_mult_2077 + (l_mult_2191 + (l_mult_2079 + (l_mult_2315 + l_mult_2158)))))))))))))),
                l_mult_2329 + (l_mult_2330 + (l_mult_2135 + (l_mult_2331 + (-0.14580e5 * l_prod_968 + (l_mult_2099 + (l_mult_2200 + (l_mult_2170 + l_sum_2285))))))),
                l_mult_2332 + (l_mult_2184 + (l_mult_2333 + (l_mult_2119 + (l_mult_2334 + l_mult_2187)))),
                l_mult_2316 + (l_mult_2335 + (l_mult_2336 + (l_mult_2295 + (l_mult_2337 + (l_mult_2319 + (l_mult_2338 + l_sum_2289)))))))) + (vdot4(
            v_929,
            vcons4(l_mult_2332 + (l_mult_2339 + (l_mult_2340 + (l_mult_2118 + (l_mult_2119 + l_mult_2120)))),
                l_mult_573 + (l_mult_2341 + (l_mult_2309 + l_mult_2144)),
                l_mult_2342 + (l_mult_2343 + (l_mult_2344 + (l_mult_2345 + (l_mult_2133 + (l_mult_2292 + (l_mult_2293 + (l_mult_2294 + (l_mult_2295 + (0.6426e4 * l_prod_969 + (l_mult_2346 + (l_mult_2099 + (-0.5832e4 * l_prod_965 + (l_mult_2321 + (l_mult_1999 + (l_mult_2347 + (0.25704e5 * l_prod_960 + (-0.34992e5 * l_prod_959 + (l_mult_2348 + (l_mult_2349 + (l_mult_1125 + (l_mult_2350 + (l_mult_2303 + (l_mult_2304 + (l_mult_2351 + (0.19278e5 * l_prod_947 + (l_mult_2352 + (l_mult_2353 + (l_mult_2354 + (l_mult_2355 + (l_mult_2241 + (l_mult_1794 + (l_mult_2356 + (l_mult_2357 + l_mult_2264))))))))))))))))))))))))))))))))),
                l_mult_2358 + (l_mult_2359 + (l_mult_2360 + (l_mult_2138 + (l_mult_2307 + (l_mult_2308 + (l_mult_2309 + (l_mult_2361 + (l_mult_2119 + (l_mult_2334 + (l_mult_2362 + (l_mult_2363 + (l_mult_1061 + (l_mult_2074 + (l_mult_2364 + (l_mult_1196 + (l_mult_2077 + (l_mult_984 + (l_mult_2079 + (l_mult_2080 + (-0.28836e5 * l_prod_947 + (l_mult_2365 + (l_mult_2366 + (l_mult_2367 + (l_mult_2368 + (l_mult_2222 + (0.41472e5 * l_prod_941 + (-0.46656e5 * l_prod_1979 + (l_mult_2369 + l_mult_2257)))))))))))))))))))))))))))))) + (vdot4(
            v_930,
            vcons4(
                l_mult_1174 + (l_mult_2370 + (l_mult_2371 + (l_mult_2316 + (l_mult_2317 + (l_mult_2272 + (l_mult_2372 + (l_mult_2373 + (l_mult_2374 + (l_mult_2375 + (l_mult_2376 + (l_mult_2274 + (0.21060e5 * l_prod_947 + (l_mult_2377 + (l_mult_2110 + (l_mult_2378 + (l_mult_2305 + (l_mult_2008 + (-0.36288e5 * l_prod_941 + (l_mult_2356 + (l_mult_2357 + l_mult_2245)))))))))))))))))))),
                l_mult_1010 + (l_mult_2379 + (l_mult_573 + (l_mult_1272 + (l_mult_2380 + (l_mult_2381 + (-0.7614e4 * l_prod_947 + (l_mult_2382 + (l_mult_2383 + (l_mult_1813 + (l_mult_2088 + (l_mult_2089 + l_mult_2226))))))))))),
                l_mult_2358 + (l_mult_2384 + (l_mult_2385 + (l_mult_2386 + (l_mult_2114 + (l_mult_2307 + (l_mult_2324 + (l_mult_2325 + (l_mult_2326 + (l_mult_2361 + (l_mult_979 + (l_mult_2068 + (l_mult_2334 + (l_mult_2187 + (l_mult_2387 + (l_mult_2363 + (0.42768e5 * l_prod_959 + (-0.23328e5 * l_prod_1988 + (l_mult_2388 + (l_mult_1196 + (l_mult_2389 + (l_mult_2191 + (l_mult_2079 + (-0.4860e4 * l_prod_947 + (l_mult_1447 + (l_mult_2366 + (l_mult_2383 + (l_mult_2221 + l_sum_2390))))))))))))))))))))))))))),
                l_mult_2391 + (l_mult_2392 + (l_mult_1493 + (l_mult_2393 + (l_mult_2329 + (l_mult_2330 + (l_mult_2135 + (l_mult_2282 + (l_mult_2233 + (l_mult_2394 + (0.17496e5 * l_prod_960 + (-0.27216e5 * l_prod_959 + (l_mult_2103 + (l_mult_2395 + (l_mult_2396 + (l_mult_2106 + (l_mult_2107 + (l_mult_2003 + (l_mult_1531 + (l_mult_2397 + (l_mult_2353 + (l_mult_2263 + (l_mult_2239 + (l_mult_2398 + l_mult_2242))))))))))))))))))))))))) + (vdot4(
            v_931,
            vcons4(
                l_mult_1468 + (l_mult_2399 + (l_mult_2400 + (l_mult_2332 + (l_mult_2184 + (l_mult_2401 + (l_mult_2402 + (l_mult_2403 + (l_mult_2404 + (l_mult_2125 + (-0.2916e4 * l_prod_947 + (l_mult_2405 + (l_mult_2083 + (l_mult_2315 + (l_mult_2158 + l_sum_2390)))))))))))))),
                l_mult_1174 + (l_mult_2406 + (l_mult_2407 + (l_mult_2408 + (l_mult_2093 + (l_mult_2316 + (l_mult_2335 + (l_mult_2336 + (l_mult_2295 + (l_mult_2272 + (l_mult_2286 + (l_mult_1997 + (l_mult_2409 + (l_mult_2373 + (-0.22032e5 * l_prod_959 + (l_mult_2348 + (l_mult_2410 + (l_mult_2376 + (l_mult_2411 + (l_mult_2412 + (l_mult_2413 + l_mult_2110)))))))))))))))))))),
                l_mult_1468 + (l_mult_2414 + (l_mult_2415 + (l_mult_2138 + (l_mult_2332 + (l_mult_2339 + (l_mult_2340 + (l_mult_2416 + (l_mult_2402 + (0.12960e5 * l_prod_959 + (l_mult_2074 + (l_mult_2124 + (l_mult_2125 + (l_mult_2126 + (-0.648e3 * l_prod_947 + (l_mult_2382 + l_mult_2083))))))))))))))),
                l_mult_1010 + (l_mult_1005 + (l_mult_2417 + (l_mult_2418 + (l_mult_2062 + (l_mult_573 + (l_mult_2341 + (l_mult_2309 + (l_mult_2144 + (l_mult_2419 + (l_mult_2380 + (l_mult_2403 + l_mult_2215))))))))))))) + (vdot4(
            v_932,
            vcons4(
                l_mult_2420 + (l_mult_2292 + (0.6426e4 * l_prod_972 + (-0.5832e4 * l_prod_971 + (l_mult_1996 + (l_mult_2421 + (l_mult_2296 + (l_mult_2346 + (l_mult_2338 + (l_mult_2422 + (l_mult_2298 + (l_mult_2170 + (l_mult_2423 + (l_mult_2300 + (l_mult_2201 + (l_mult_2424 + (l_mult_2349 + (l_mult_2302 + (l_mult_2411 + (0.25704e5 * l_prod_953 + (l_mult_1127 + (l_mult_2304 + (-0.34992e5 * l_prod_950 + (l_mult_2425 + (l_mult_2426 + (0.19278e5 * l_prod_945 + (l_mult_2354 + (l_mult_2239 + (l_mult_2427 + (l_mult_2428 + (l_mult_2429 + (l_mult_1776 + (l_mult_2357 + (l_mult_2430 + l_mult_2265))))))))))))))))))))))))))))))))),
                l_mult_2431 + (l_mult_2307 + (l_mult_2432 + (l_mult_2340 + (l_mult_2433 + (l_mult_2327 + (l_mult_2119 + (l_mult_2434 + (l_mult_2323 + (l_mult_2204 + (l_mult_2435 + (l_mult_2364 + (l_mult_983 + (l_mult_2126 + (l_mult_2436 + (l_mult_1197 + (l_mult_2079 + (l_mult_1100 + (l_mult_2155 + (l_mult_2156 + (-0.28836e5 * l_prod_945 + (l_mult_2367 + (l_mult_2221 + (l_mult_2437 + (l_mult_2438 + (l_mult_2439 + (0.41472e5 * l_prod_940 + (l_mult_2369 + (-0.46656e5 * l_prod_1977 + l_mult_2258)))))))))))))))))))))))))))),
                l_mult_1175 + (l_mult_2316 + (l_mult_2278 + (l_mult_2440 + (l_mult_2337 + (l_mult_2441 + (l_mult_2442 + (l_mult_2375 + (l_mult_2279 + (l_mult_2443 + (l_mult_2444 + (l_mult_2445 + (0.21060e5 * l_prod_945 + (l_mult_2378 + (l_mult_2007 + (l_mult_2446 + (l_mult_2306 + (l_mult_2181 + (-0.36288e5 * l_prod_940 + (l_mult_2357 + (l_mult_2430 + l_mult_2246)))))))))))))))))))),
                l_mult_1028 + (l_mult_573 + (l_mult_2447 + (l_mult_1332 + (l_mult_2381 + (l_mult_2448 + (-0.7614e4 * l_prod_945 + (l_mult_2383 + (l_mult_2449 + (l_mult_1870 + (l_mult_2089 + (l_mult_2163 + l_mult_2227))))))))))))) + (vdot4(
            v_933,
            vcons4(
                l_mult_2431 + (l_mult_2307 + (l_mult_2432 + (l_mult_2340 + (l_mult_2450 + (l_mult_2310 + (l_mult_979 + (l_mult_2120 + (l_mult_2451 + (l_mult_2311 + (l_mult_2147 + (l_mult_2452 + (l_mult_2313 + (l_mult_2188 + (l_mult_2453 + (l_mult_2388 + (l_mult_2125 + (l_mult_2436 + (l_mult_1197 + (l_mult_2079 + (0.42768e5 * l_prod_950 + (l_mult_2454 + (-0.23328e5 * l_prod_1984 + (-0.4860e4 * l_prod_945 + (l_mult_2383 + (l_mult_1409 + (l_mult_2222 + (l_mult_2439 + l_sum_2455))))))))))))))))))))))))))),
                l_mult_2456 + (l_mult_2329 + (l_mult_2281 + (l_mult_2457 + (l_mult_2331 + (l_mult_2233 + (l_mult_1685 + (l_mult_2200 + (l_mult_2458 + (l_mult_2459 + (l_mult_2395 + (l_mult_2174 + (0.17496e5 * l_prod_953 + (l_mult_2460 + (l_mult_2003 + (-0.27216e5 * l_prod_950 + (l_mult_2177 + (l_mult_2178 + (l_mult_1720 + (l_mult_2263 + (l_mult_2461 + (l_mult_2241 + (l_mult_2429 + (l_mult_2462 + l_mult_2244))))))))))))))))))))))),
                l_mult_1649 + (l_mult_2332 + (l_mult_2463 + (l_mult_2118 + (l_mult_2464 + (l_mult_2465 + (l_mult_2404 + (l_mult_2466 + (l_mult_2191 + (l_mult_2467 + (-0.2916e4 * l_prod_945 + (l_mult_2315 + (l_mult_2468 + (l_mult_2086 + (l_mult_2161 + l_sum_2455)))))))))))))),
                l_mult_1175 + (l_mult_2316 + (l_mult_2278 + (l_mult_2469 + (l_mult_2318 + (l_mult_2286 + (l_mult_2470 + (l_mult_2320 + (l_mult_1998 + (l_mult_2471 + (l_mult_2300 + (l_mult_2172 + (l_mult_2472 + (l_mult_2410 + (l_mult_2443 + (l_mult_2444 + (-0.22032e5 * l_prod_950 + (l_mult_2351 + (l_mult_2426 + (l_mult_2473 + (l_mult_2474 + l_mult_2181)))))))))))))))))))))) + (vdot4(
            v_934,
            vcons4(
                l_mult_1649 + (l_mult_2332 + (l_mult_2475 + (l_mult_2333 + (l_mult_2476 + (l_mult_2334 + (l_mult_2204 + (l_mult_2477 + (l_mult_2124 + (l_mult_2466 + (l_mult_2191 + (0.12960e5 * l_prod_950 + (l_mult_2080 + (l_mult_2156 + (-0.648e3 * l_prod_945 + (l_mult_2449 + l_mult_2161))))))))))))))),
                l_mult_1028 + (l_mult_573 + (l_mult_1023 + (l_mult_2322 + (l_mult_2478 + (l_mult_2323 + (l_mult_2479 + (l_mult_2071 + (l_mult_2150 + (l_mult_2480 + (l_mult_2448 + (l_mult_2467 + l_mult_2220))))))))))),
                0.4320e4 * l_prod_517 + (-0.15984e5 * l_prod_972 + (0.19440e5 * l_prod_971 + (l_mult_2065 + (-0.15984e5 * l_prod_969 + (0.38880e5 * l_prod_968 + (l_mult_2328 + (0.19440e5 * l_prod_965 + (l_mult_2312 + (l_mult_2149 + (l_mult_2481 + (l_mult_2482 + (l_mult_2389 + (l_mult_2483 + (-0.93312e5 * l_prod_1986 + (l_mult_2454 + (0.58320e5 * l_prod_944 + (l_mult_2368 + (l_mult_2438 + l_mult_2224)))))))))))))))))),
                l_mult_2484 + (l_mult_2485 + (l_mult_2135 + (l_mult_2486 + (l_mult_2319 + (l_mult_2200 + (l_mult_2487 + (l_mult_2488 + (l_mult_2106 + (l_mult_2489 + (l_mult_2304 + (l_mult_2177 + (l_mult_1438 + (l_mult_2355 + (l_mult_2428 + l_mult_2243)))))))))))))))) + (vdot4(
            v_935,
            vcons4(
                l_mult_1183 + (l_mult_2490 + (l_mult_2491 + (l_mult_2492 + (l_mult_2493 + (l_mult_2494 + (0.34992e5 * l_prod_944 + (l_mult_2085 + (l_mult_2160 + l_mult_2224)))))))),
                l_mult_2484 + (l_mult_2485 + (l_mult_2135 + (0.13284e5 * l_prod_969 + (l_mult_2495 + (l_mult_2099 + (l_mult_1080 + (l_mult_2299 + (l_mult_2171 + (l_mult_2496 + (l_mult_2376 + (l_mult_2489 + (l_mult_2304 + (l_mult_2425 + (l_mult_2263 + l_mult_2241)))))))))))))),
                l_mult_1187 + (l_mult_2339 + (l_mult_2497 + (l_mult_2119 + (l_mult_2323 + (l_mult_2498 + (l_mult_2125 + (l_mult_2154 + (l_mult_2079 + (l_mult_2155 + (l_mult_2383 + l_mult_2222)))))))))),
                l_mult_1183 + (l_mult_2490 + (-0.4320e4 * l_prod_969 + (l_mult_2499 + (0.11664e5 * l_prod_965 + (l_mult_2070 + (l_mult_2149 + (l_mult_2500 + (l_mult_2494 + l_mult_2219)))))))))) + (vdot4(
            v_936,
            vcons4(
                l_mult_2484 + (0.13284e5 * l_prod_972 + (l_mult_1037 + (l_mult_2096 + (l_mult_2486 + (l_mult_2495 + (l_mult_2297 + (l_mult_2200 + (l_mult_2170 + (l_mult_2496 + (l_mult_2488 + (l_mult_2350 + (l_mult_2444 + (l_mult_2304 + (l_mult_2263 + l_mult_2239)))))))))))))),
                l_mult_1187 + (l_mult_2501 + (l_mult_2309 + (l_mult_2333 + (l_mult_2119 + (l_mult_2498 + (l_mult_2076 + (l_mult_2077 + (l_mult_2191 + (l_mult_2079 + (l_mult_2383 + l_mult_2221)))))))))),
                l_mult_1187 + (l_mult_2501 + (l_mult_2309 + (l_mult_2497 + (l_mult_2067 + (l_mult_2068 + (l_mult_2323 + (l_mult_2147 + (l_mult_2404 + (l_mult_2125 + (l_mult_2191 + l_mult_2079)))))))))),
                l_mult_1183 + (-0.4320e4 * l_prod_972 + (0.11664e5 * l_prod_971 + (l_mult_2065 + (l_mult_2491 + (l_mult_2499 + (l_mult_2146 + (l_mult_2500 + (l_mult_2493 + l_mult_2218)))))))))) + vdot4(
            v_916,
            vcons4(l_sum_2016,
                l_mult_2017 + (0.274e2 * l_prod_512 + (-0.2025e3 * l_prod_507 + (0.612e3 * l_prod_942 + (-0.810e3 * l_prod_939 + l_mult_2015)))),
                0.e0, 0.e0)))))))))))))))))))));
        double l_vdot_2864 = vdot4(v_917,
            vcons4(l_mult_2018 + (l_mult_2502 + (0.1134e4 * l_prod_969 + (-0.2592e4 * l_prod_965 + l_mult_1999))),
                l_mult_547 + (l_mult_561 + (l_mult_2503 + (l_mult_2134 + (0.486e3 * l_prod_969 + (l_mult_2504 + l_sum_2287))))),
                l_mult_2025 + (l_mult_2026 + (l_mult_2027 + (l_mult_2505 + (l_mult_2506 + (-0.1296e4 * l_prod_971 + l_sum_2290))))),
                l_mult_547 + (l_mult_2033 + (l_mult_2034 + (l_mult_2035 + l_sum_2291))))) + (vdot4(v_918,
            vcons4(l_sum_2038, 0.e0, 0.e0, 0.e0)) + (vdot4(v_919,
            vcons4(0.e0, 0.e0, l_sum_2508,
                l_mult_577 + (l_mult_566 + (l_mult_2509 + (l_mult_2510 + (l_mult_2511 + (l_mult_2512 + l_sum_2514))))))) + (vdot4(
            v_920,
            vcons4(
                l_mult_2515 + (l_mult_2049 + (l_mult_2516 + (l_mult_2517 + (l_mult_2473 + (l_mult_2052 + (l_mult_2518 + (-0.1296e4 * l_prod_940 + l_mult_2012))))))),
                l_mult_577 + (l_mult_2041 + (0.486e3 * l_prod_953 + (l_mult_2519 + l_sum_2520))),
                l_mult_2507 + (l_mult_2040 + (0.1134e4 * l_prod_953 + (-0.2592e4 * l_prod_950 + l_mult_2005))),
                l_sum_2091)) + (vdot4(v_921, vcons4(l_sum_2112, l_sum_2131, l_sum_2136, l_sum_2140)) + (vdot4(v_922,
            vcons4(
                l_mult_2207 + (l_mult_2060 + (l_mult_2208 + (l_mult_2209 + (l_mult_2210 + (l_mult_2139 + (-0.6264e3 * l_prod_518 + (l_mult_2521 + (l_mult_2064 + (l_mult_2117 + (l_mult_2522 + (0.3132e4 * l_prod_515 + (-0.15066e5 * l_prod_969 + (l_mult_1502 + (l_mult_2068 + (-0.6696e4 * l_prod_966 + (0.20736e5 * l_prod_965 + (l_mult_2523 + (0.6480e4 * l_prod_963 + (l_mult_2524 + (l_mult_2525 + (l_mult_2526 + (l_mult_2072 + (l_mult_2527 + (l_mult_2528 + (l_mult_2529 + (l_mult_2216 + (l_mult_2217 + (l_mult_2076 + (l_mult_2218 + (-0.15066e5 * l_prod_953 + (l_mult_1197 + (l_mult_2530 + (0.20736e5 * l_prod_950 + (l_mult_2531 + (l_mult_2532 + (l_mult_2533 + (l_mult_2081 + (0.7776e4 * l_prod_946 + (l_mult_2129 + (l_mult_2534 + (l_mult_2256 + (l_mult_2085 + (l_mult_1671 + (l_mult_2222 + (l_mult_2535 + (l_mult_2536 + (l_mult_2087 + (l_mult_2537 + (l_mult_2538 + (l_mult_2539 + (l_mult_2540 + (l_mult_2541 + (l_mult_2090 + (l_mult_2542 + l_mult_2543)))))))))))))))))))))))))))))))))))))))))))))))))))))),
                l_mult_2229 + (l_mult_2092 + (l_mult_2230 + (l_mult_2231 + (l_mult_2035 + (0.1053e4 * l_prod_518 + (l_mult_2544 + (l_mult_2095 + (-0.7128e4 * l_prod_971 + (l_mult_1996 + (-0.62235e4 * l_prod_515 + (0.23652e5 * l_prod_969 + (l_mult_2545 + (l_mult_2099 + (0.14796e5 * l_prod_966 + (-0.37584e5 * l_prod_965 + (l_mult_2299 + (-0.15390e5 * l_prod_963 + (l_mult_2546 + (l_mult_2547 + (l_mult_2548 + (l_mult_2101 + (l_mult_2549 + (l_mult_2550 + (l_mult_2236 + (l_mult_2237 + (l_mult_2105 + (l_mult_2002 + (0.23652e5 * l_prod_953 + (-0.58320e5 * l_prod_952 + (l_mult_2551 + (-0.37584e5 * l_prod_950 + (l_mult_2425 + (l_mult_2552 + (l_mult_2553 + (l_mult_2108 + (l_mult_2554 + (l_mult_2555 + (l_mult_2556 + (l_mult_2007 + (l_mult_2240 + (l_mult_2241 + (l_mult_2557 + (l_mult_2558 + (l_mult_2023 + (-0.7128e4 * l_prod_940 + (l_mult_2011 + (l_mult_2559 + l_sum_2514))))))))))))))))))))))))))))))))))))))))))))))),
                l_mult_2248 + (l_mult_2113 + (l_mult_2249 + (l_mult_2250 + (-0.1016e4 * l_prod_518 + (l_mult_2560 + (l_mult_2116 + (0.1296e4 * l_prod_971 + (0.6696e4 * l_prod_515 + (-0.18360e5 * l_prod_969 + (l_mult_2067 + (l_mult_2120 + (-0.17424e5 * l_prod_966 + (0.33696e5 * l_prod_965 + (l_mult_2523 + (l_mult_1401 + (l_mult_2561 + (-0.7776e4 * l_prod_1989 + (l_mult_2562 + (l_mult_2121 + (l_mult_2563 + (l_mult_2254 + (l_mult_2255 + (l_mult_2125 + (-0.18360e5 * l_prod_953 + (l_mult_2154 + (l_mult_2564 + (0.33696e5 * l_prod_950 + (l_mult_2531 + (l_mult_2565 + (l_mult_2566 + (l_mult_2127 + (l_mult_2567 + (l_mult_2568 + (l_mult_2159 + (l_mult_2086 + (l_mult_2535 + (l_mult_2569 + (0.1296e4 * l_prod_940 + l_mult_2570)))))))))))))))))))))))))))))))))))))),
                l_mult_2259 + (l_mult_2132 + (l_mult_561 + (0.594e3 * l_prod_518 + (l_mult_2571 + (l_mult_2134 + (-0.41445e4 * l_prod_515 + (0.7128e4 * l_prod_969 + (l_mult_2504 + (0.11556e5 * l_prod_966 + (-0.14904e5 * l_prod_965 + (l_mult_1998 + (-0.13770e5 * l_prod_963 + (l_mult_2572 + (l_mult_2547 + (l_mult_2573 + (l_mult_562 + (l_mult_2261 + (l_mult_2262 + (0.7128e4 * l_prod_953 + (l_mult_2574 + (-0.14904e5 * l_prod_950 + (l_mult_2004 + (l_mult_2575 + l_sum_2520))))))))))))))))))))))))) + (vdot4(
            v_923,
            vcons4(
                l_mult_1354 + (l_mult_2137 + (-0.1944e3 * l_prod_518 + (l_mult_2576 + (0.1404e4 * l_prod_515 + (-0.1134e4 * l_prod_969 + (-0.4104e4 * l_prod_966 + (l_mult_2577 + (l_mult_2148 + (l_mult_2071 + (l_mult_2525 + (l_mult_2578 + (l_mult_2267 + (-0.1134e4 * l_prod_953 + (l_mult_2579 + l_mult_2580)))))))))))))),
                l_sum_2591, l_sum_2599, l_sum_2604)) + (vdot4(v_924,
            vcons4(l_sum_2608, l_sum_2610, l_mult_562 + (l_mult_2605 + (l_mult_2606 + l_mult_2013)),
                l_mult_2611 + (l_mult_2276 + (l_mult_2612 + (l_mult_2613 + (l_mult_2023 + l_mult_2011)))))) + (vdot4(
            v_925,
            vcons4(l_mult_2611 + (l_mult_2273 + (l_mult_2107 + l_sum_2614)),
                l_mult_562 + (l_mult_2262 + (l_mult_2574 + l_mult_2004)),
                l_mult_2611 + (l_mult_2615 + (l_mult_2612 + (l_mult_2031 + l_sum_2024))),
                l_mult_527 + (l_mult_2616 + (l_mult_2283 + (l_mult_2174 + l_sum_2617))))) + (vdot4(v_926,
            vcons4(l_mult_2611 + (l_mult_2615 + (l_mult_2273 + (l_mult_2279 + l_sum_2284))),
                l_mult_2611 + (l_mult_2618 + (l_mult_2550 + l_sum_2032)),
                l_mult_2611 + (l_mult_2618 + (l_mult_2550 + l_sum_2289)), l_sum_2037)) + (vdot4(v_927,
            vcons4(
                l_mult_2342 + (l_mult_2343 + (l_mult_2344 + (l_mult_2345 + (l_mult_2133 + (l_mult_2619 + (0.25704e5 * l_prod_972 + (-0.34992e5 * l_prod_971 + (l_mult_2620 + (0.19278e5 * l_prod_969 + (l_mult_2621 + (l_mult_2622 + (l_mult_1080 + (l_mult_2623 + (l_mult_2572 + (l_mult_2624 + (l_mult_2625 + (l_mult_2626 + (l_mult_2627 + (l_mult_2349 + (l_mult_1125 + (l_mult_2350 + (l_mult_2628 + (l_mult_2629 + (l_mult_2630 + (0.6426e4 * l_prod_947 + (l_mult_2631 + (l_mult_2110 + (l_mult_2632 + (l_mult_2633 + (l_mult_2241 + (-0.5832e4 * l_prod_941 + (l_mult_2634 + (l_mult_2635 + l_mult_2013))))))))))))))))))))))))))))))))),
                l_mult_2358 + (l_mult_2359 + (l_mult_2360 + (l_mult_2138 + (l_mult_2636 + (l_mult_2637 + (l_mult_1057 + (l_mult_2065 + (-0.28836e5 * l_prod_969 + (l_mult_2638 + (l_mult_2639 + (0.41472e5 * l_prod_965 + (-0.46656e5 * l_prod_1991 + (l_mult_2561 + (l_mult_2640 + (l_mult_2641 + (l_mult_2642 + (l_mult_2364 + (l_mult_1196 + (l_mult_2077 + (l_mult_2643 + (l_mult_2644 + (l_mult_2454 + (l_mult_2645 + (l_mult_2128 + (l_mult_987 + (l_mult_2085 + (l_mult_2222 + (l_mult_2646 + l_mult_2089)))))))))))))))))))))))))))),
                l_mult_1174 + (l_mult_2370 + (l_mult_2371 + (l_mult_2647 + (l_mult_2648 + (l_mult_2649 + (0.21060e5 * l_prod_969 + (l_mult_2650 + (l_mult_2099 + (-0.36288e5 * l_prod_965 + (l_mult_2623 + (l_mult_2546 + (l_mult_2651 + (l_mult_2652 + (l_mult_2375 + (l_mult_2376 + (l_mult_2653 + (l_mult_2654 + (l_mult_2630 + l_sum_2614)))))))))))))))))),
                l_mult_1010 + (l_mult_2379 + (l_mult_1832 + (l_mult_2655 + (-0.7614e4 * l_prod_969 + (l_mult_2656 + (l_mult_1210 + (l_mult_2070 + (l_mult_2524 + (l_mult_554 + (l_mult_2381 + (l_mult_2657 + l_mult_2080))))))))))))) + (vdot4(
            v_928,
            vcons4(
                l_mult_2358 + (l_mult_2384 + (l_mult_2385 + (l_mult_2386 + (l_mult_2114 + (l_mult_2658 + (l_mult_2637 + (0.42768e5 * l_prod_971 + (-0.23328e5 * l_prod_1993 + (-0.4860e4 * l_prod_969 + (l_mult_1400 + (l_mult_2639 + (l_mult_2577 + (l_mult_2523 + (l_mult_2640 + (l_mult_2659 + (l_mult_2660 + (l_mult_2661 + (l_mult_2388 + (l_mult_1196 + (l_mult_2389 + (l_mult_2657 + (l_mult_2530 + (l_mult_2645 + (l_mult_986 + (l_mult_2083 + (l_mult_2568 + (l_mult_2085 + l_sum_2662))))))))))))))))))))))))))),
                l_mult_2391 + (l_mult_2392 + (l_mult_1493 + (l_mult_2393 + (l_mult_2663 + (0.17496e5 * l_prod_972 + (-0.27216e5 * l_prod_971 + (l_mult_2096 + (l_mult_1504 + (l_mult_2664 + (l_mult_2622 + (l_mult_2665 + (l_mult_2299 + (l_mult_2666 + (l_mult_2667 + (l_mult_2036 + (l_mult_2395 + (l_mult_2396 + (l_mult_2106 + (l_mult_2574 + (l_mult_2551 + l_sum_2617)))))))))))))))))))),
                l_mult_1468 + (l_mult_2399 + (l_mult_2400 + (l_mult_2668 + (l_mult_2669 + (l_mult_2670 + (-0.2916e4 * l_prod_969 + (l_mult_2671 + (l_mult_2068 + (l_mult_2577 + (l_mult_2523 + (l_mult_2672 + (l_mult_2563 + (l_mult_2404 + (l_mult_2125 + (l_mult_2673 + l_mult_2564))))))))))))))),
                l_mult_1174 + (l_mult_2406 + (l_mult_2407 + (l_mult_2408 + (l_mult_2093 + (l_mult_2674 + (l_mult_2648 + (-0.22032e5 * l_prod_971 + (l_mult_2620 + (l_mult_2675 + (l_mult_2676 + (l_mult_2099 + (l_mult_2651 + (l_mult_2677 + (l_mult_2678 + (l_mult_2627 + (l_mult_2410 + (l_mult_2376 + (l_mult_2411 + l_sum_2032)))))))))))))))))))) + (vdot4(
            v_929,
            vcons4(
                l_mult_1468 + (l_mult_2414 + (l_mult_2415 + (l_mult_2138 + (l_mult_2679 + (l_mult_2669 + (0.12960e5 * l_prod_971 + (l_mult_2065 + (-0.648e3 * l_prod_969 + (l_mult_2656 + (l_mult_2068 + (l_mult_2672 + (l_mult_2680 + (l_mult_2681 + (l_mult_2124 + (l_mult_2125 + l_mult_2126))))))))))))))),
                l_mult_1010 + (l_mult_1005 + (l_mult_2417 + (l_mult_2418 + (l_mult_2062 + (l_mult_2682 + (l_mult_2655 + (l_mult_2670 + (l_mult_2522 + l_sum_2683)))))))),
                l_mult_2624 + (l_mult_2625 + (l_mult_2626 + (l_mult_2627 + (l_mult_2301 + (l_mult_2302 + (l_mult_2106 + (-0.17496e5 * l_prod_952 + (l_mult_2654 + (l_mult_2004 + (l_mult_2684 + (l_mult_1437 + (l_mult_2685 + (l_mult_2632 + (l_mult_2633 + (l_mult_2306 + (l_mult_2686 + (l_mult_2242 + (l_mult_2595 + l_mult_2687)))))))))))))))))),
                l_mult_2640 + (l_mult_2641 + (l_mult_2642 + (l_mult_2314 + (l_mult_2125 + (l_mult_2673 + (l_mult_2688 + (l_mult_1489 + (l_mult_2083 + (l_mult_987 + (l_mult_2085 + (l_mult_2086 + (l_mult_2689 + (l_mult_2690 + (l_mult_2587 + l_mult_2691)))))))))))))))) + (vdot4(
            v_930,
            vcons4(
                l_mult_2651 + (l_mult_2652 + (l_mult_2276 + (l_mult_2692 + (l_mult_2693 + (l_mult_2613 + (l_mult_2694 + (l_mult_2634 + (l_mult_2011 + l_mult_2687)))))))),
                l_mult_554 + (l_mult_2695 + (l_mult_2696 + l_mult_2090)),
                l_mult_2640 + (l_mult_2659 + (l_mult_2660 + (l_mult_2661 + (l_mult_2314 + (l_mult_983 + (l_mult_2077 + (l_mult_2673 + (l_mult_2564 + (l_mult_2697 + (l_mult_1489 + (l_mult_2698 + l_sum_2699))))))))))),
                l_mult_2666 + (l_mult_2667 + (l_mult_2036 + (l_mult_2283 + (l_mult_2174 + (l_mult_2700 + (-0.14580e5 * l_prod_946 + (l_mult_2110 + (l_mult_2111 + (l_mult_2007 + l_sum_2701))))))))))) + (vdot4(
            v_931,
            vcons4(l_mult_2672 + (l_mult_2563 + (l_mult_2702 + (l_mult_2128 + l_sum_2662))),
                l_mult_2651 + (l_mult_2677 + (l_mult_2678 + (l_mult_2627 + (l_mult_2276 + (l_mult_2279 + (l_mult_2002 + (l_mult_2703 + (l_mult_2693 + l_mult_2704)))))))),
                l_mult_2672 + (l_mult_2680 + (l_mult_2681 + l_sum_2130)), l_sum_2683)) + (vdot4(v_932,
            vcons4(
                l_mult_2705 + (l_mult_2624 + (0.6426e4 * l_prod_960 + (-0.5832e4 * l_prod_959 + (l_mult_2001 + (l_mult_2424 + (l_mult_2349 + (l_mult_2302 + (l_mult_2411 + (0.19278e5 * l_prod_953 + (l_mult_2628 + (l_mult_2551 + (l_mult_1086 + (l_mult_2630 + (l_mult_2575 + (l_mult_2706 + (l_mult_2684 + (l_mult_2631 + (l_mult_2704 + (0.25704e5 * l_prod_945 + (l_mult_1438 + (l_mult_2633 + (l_mult_2427 + (l_mult_2428 + (l_mult_2707 + (l_mult_2708 + (l_mult_2686 + (l_mult_2594 + (-0.34992e5 * l_prod_940 + (l_mult_2243 + (l_mult_2709 + (l_mult_2710 + (l_mult_2687 + (l_mult_2711 + l_mult_2607))))))))))))))))))))))))))))))))),
                l_mult_2712 + (l_mult_2640 + (l_mult_2713 + (l_mult_2681 + (l_mult_2453 + (l_mult_2388 + (l_mult_2125 + (-0.4860e4 * l_prod_953 + (l_mult_2657 + (l_mult_2579 + (l_mult_2714 + (l_mult_2688 + (l_mult_986 + (l_mult_2129 + (l_mult_2715 + (l_mult_1490 + (l_mult_2085 + (l_mult_1409 + (l_mult_2222 + (l_mult_2535 + (l_mult_2716 + (l_mult_2689 + (l_mult_2586 + (0.42768e5 * l_prod_940 + (l_mult_2369 + (l_mult_2717 + (l_mult_2718 + (l_mult_2691 + (-0.23328e5 * l_prod_1975 + l_mult_2603)))))))))))))))))))))))))))),
                l_mult_1465 + (l_mult_2651 + (l_mult_2615 + (l_mult_2472 + (l_mult_2410 + (l_mult_2050 + (l_mult_2719 + (l_mult_2692 + (l_mult_2031 + (l_mult_2720 + (l_mult_2721 + (l_mult_2474 + (l_mult_2722 + (l_mult_2694 + (l_mult_2010 + (-0.22032e5 * l_prod_940 + (l_mult_2635 + (l_mult_2559 + (l_mult_2723 + (l_mult_2687 + (l_mult_2711 + l_mult_2598)))))))))))))))))))),
                l_mult_1372 + (l_mult_554 + (l_mult_2480 + (l_mult_1596 + (l_mult_2695 + (l_mult_2724 + (l_mult_2725 + (l_mult_2696 + (l_mult_2726 + (l_mult_2727 + (l_mult_2090 + (l_mult_2542 + l_mult_2590))))))))))))) + (vdot4(
            v_933,
            vcons4(
                l_mult_2712 + (l_mult_2640 + (l_mult_2713 + (l_mult_2681 + (l_mult_2435 + (l_mult_2364 + (l_mult_983 + (l_mult_2126 + (-0.28836e5 * l_prod_953 + (l_mult_2643 + (l_mult_2530 + (0.41472e5 * l_prod_950 + (l_mult_2454 + (l_mult_2565 + (l_mult_2728 + (l_mult_2697 + (l_mult_2128 + (l_mult_2715 + (l_mult_1490 + (l_mult_2085 + (l_mult_2437 + (l_mult_2438 + (-0.46656e5 * l_prod_1980 + (l_mult_2729 + (l_mult_2696 + (l_mult_1784 + (l_mult_2587 + (l_mult_2717 + l_sum_2730))))))))))))))))))))))))))),
                l_mult_2731 + (l_mult_2666 + (l_mult_2616 + (l_mult_2459 + (l_mult_2395 + (l_mult_2174 + (l_mult_1888 + (l_mult_2574 + (l_mult_2445 + (l_mult_2732 + (l_mult_2700 + (l_mult_2554 + (0.17496e5 * l_prod_945 + (l_mult_2733 + (l_mult_2007 + (l_mult_2461 + (l_mult_2241 + (l_mult_2557 + (l_mult_1908 + (l_mult_2606 + (-0.27216e5 * l_prod_940 + (l_mult_2595 + (l_mult_2709 + (l_mult_2734 + l_mult_2597))))))))))))))))))))))),
                l_mult_1856 + (l_mult_2672 + (l_mult_2477 + (l_mult_2124 + (-0.648e3 * l_prod_953 + (l_mult_2735 + (l_mult_2702 + (l_mult_2736 + (l_mult_2568 + (l_mult_2449 + (l_mult_2737 + (l_mult_2646 + (0.12960e5 * l_prod_940 + (l_mult_2089 + (l_mult_2540 + l_sum_2730)))))))))))))),
                l_mult_1465 + (l_mult_2651 + (l_mult_2615 + (l_mult_2442 + (l_mult_2375 + (l_mult_2279 + (0.21060e5 * l_prod_953 + (l_mult_2653 + (l_mult_2003 + (-0.36288e5 * l_prod_950 + (l_mult_2630 + (l_mult_2552 + (l_mult_2738 + (l_mult_2703 + (l_mult_2720 + (l_mult_2721 + (l_mult_2446 + (l_mult_2306 + (l_mult_2707 + (l_mult_2739 + (l_mult_2462 + l_mult_2559)))))))))))))))))))))) + (vdot4(
            v_934,
            vcons4(
                l_mult_1856 + (l_mult_2672 + (l_mult_2465 + (l_mult_2404 + (-0.2916e4 * l_prod_953 + (l_mult_2673 + (l_mult_2579 + (l_mult_2740 + (l_mult_2127 + (l_mult_2736 + (l_mult_2568 + (l_mult_2468 + (l_mult_2086 + (l_mult_2535 + (l_mult_2741 + (l_mult_2726 + l_mult_2540))))))))))))))),
                l_mult_1372 + (l_mult_554 + (l_mult_1332 + (l_mult_2381 + (-0.7614e4 * l_prod_953 + (l_mult_2657 + (l_mult_1699 + (l_mult_2080 + (l_mult_2532 + (l_mult_2742 + (l_mult_2724 + (l_mult_2449 + l_mult_2585))))))))))),
                0.4320e4 * l_prod_511 + (-0.15984e5 * l_prod_960 + (0.19440e5 * l_prod_959 + (l_mult_2074 + (l_mult_2481 + (l_mult_2482 + (l_mult_2389 + (0.58320e5 * l_prod_952 + (l_mult_2644 + (l_mult_2531 + (-0.15984e5 * l_prod_947 + (0.38880e5 * l_prod_946 + (l_mult_2698 + (l_mult_2743 + (-0.93312e5 * l_prod_1982 + (l_mult_2438 + (0.19440e5 * l_prod_941 + (l_mult_2690 + (l_mult_2369 + l_mult_2588)))))))))))))))))),
                l_mult_2744 + (l_mult_2745 + (l_mult_2036 + (l_mult_2496 + (l_mult_2376 + (l_mult_2574 + (0.13284e5 * l_prod_947 + (l_mult_2746 + (l_mult_2110 + (l_mult_2747 + (l_mult_2633 + (l_mult_2241 + (l_mult_1794 + (l_mult_2242 + (l_mult_2243 + l_mult_2596)))))))))))))))) + (vdot4(
            v_935,
            vcons4(
                l_mult_1007 + (l_mult_2748 + (l_mult_2500 + (-0.4320e4 * l_prod_947 + (l_mult_2749 + (l_mult_2750 + (0.11664e5 * l_prod_941 + (l_mult_2088 + (l_mult_2539 + l_mult_2588)))))))),
                l_mult_2744 + (l_mult_2745 + (l_mult_2036 + (l_mult_2487 + (l_mult_2488 + (l_mult_2106 + (l_mult_1127 + (l_mult_2629 + (l_mult_2425 + (l_mult_2751 + (l_mult_2693 + (l_mult_2747 + (l_mult_2633 + (l_mult_2428 + (l_mult_2606 + l_mult_2595)))))))))))))),
                l_mult_1475 + (l_mult_2680 + (l_mult_2498 + (l_mult_2125 + (l_mult_2657 + (l_mult_2752 + (l_mult_2128 + (l_mult_2256 + (l_mult_2085 + (l_mult_2222 + (l_mult_2696 + l_mult_2587)))))))))),
                l_mult_1007 + (l_mult_2748 + (l_mult_2492 + (l_mult_2493 + (0.34992e5 * l_prod_952 + (l_mult_2079 + (l_mult_2531 + (l_mult_2753 + (l_mult_2750 + l_mult_2160)))))))))) + (vdot4(
            v_936,
            vcons4(
                l_mult_2744 + (0.13284e5 * l_prod_960 + (l_mult_1043 + (l_mult_2103 + (l_mult_2496 + (l_mult_2488 + (l_mult_2350 + (l_mult_2574 + (l_mult_2551 + (l_mult_2751 + (l_mult_2746 + (l_mult_2685 + (l_mult_2721 + (l_mult_2633 + l_sum_2701))))))))))))),
                l_mult_1475 + (l_mult_2754 + (l_mult_2642 + (l_mult_2404 + (l_mult_2125 + (l_mult_2752 + (l_mult_2082 + (l_mult_2083 + l_sum_2699))))))),
                l_mult_1475 + (l_mult_2754 + (l_mult_2642 + (l_mult_2498 + (l_mult_2076 + (l_mult_2077 + (l_mult_2657 + (l_mult_2530 + (l_mult_2702 + (l_mult_2128 + (l_mult_2568 + l_mult_2085)))))))))),
                l_mult_1007 + (-0.4320e4 * l_prod_960 + (0.11664e5 * l_prod_959 + (l_mult_2074 + (l_mult_2500 + (l_mult_2493 + (l_mult_2218 + (l_mult_2753 + (l_mult_2749 + l_mult_2584)))))))))) + vdot4(
            v_916,
            vcons4(l_sum_2016, 0.e0,
                l_mult_2017 + (0.274e2 * l_prod_518 + (-0.2025e3 * l_prod_515 + (0.612e3 * l_prod_966 + (-0.810e3 * l_prod_963 + l_mult_2000)))),
                0.e0)))))))))))))))))))));
        double l_vdot_2865 = vdot4(v_917,
            vcons4(l_sum_2059,
                l_mult_551 + (l_mult_581 + (l_mult_2054 + (l_mult_2199 + (l_mult_2055 + (l_mult_2200 + (l_mult_2056 + l_mult_1999)))))),
                l_mult_2046 + (l_mult_2505 + (l_mult_2278 + (l_mult_2047 + (l_mult_2675 + (l_mult_2286 + (l_mult_2048 + (-0.1296e4 * l_prod_965 + l_mult_1998))))))),
                l_mult_551 + (l_mult_2503 + (0.486e3 * l_prod_972 + (l_mult_2168 + (l_mult_565 + (l_mult_2199 + (l_mult_2504 + l_mult_1997)))))))) + (vdot4(
            v_918,
            vcons4(l_mult_2039 + (l_mult_2502 + (0.1134e4 * l_prod_972 + (-0.2592e4 * l_prod_971 + l_mult_1996))),
                l_sum_2508,
                l_mult_577 + (l_mult_562 + (l_mult_2509 + (l_mult_2605 + (l_mult_2511 + (l_mult_2606 + l_sum_2755))))),
                l_mult_2515 + (l_mult_2028 + (l_mult_2615 + (l_mult_2517 + (l_mult_2412 + (l_mult_2031 + (l_mult_2518 + (-0.1296e4 * l_prod_941 + l_mult_2010))))))))) + (vdot4(
            v_919,
            vcons4(l_mult_577 + (l_mult_2020 + (0.486e3 * l_prod_960 + (l_mult_2550 + l_sum_2756))),
                l_mult_2507 + (l_mult_2019 + (0.1134e4 * l_prod_960 + (-0.2592e4 * l_prod_959 + l_mult_2001))), 0.e0,
                0.e0)) + (vdot4(v_920,
            vcons4(0.e0, 0.e0, 0.e0,
                l_mult_2207 + (-0.6264e3 * l_prod_522 + (0.3132e4 * l_prod_521 + (-0.6696e4 * l_prod_976 + (0.6480e4 * l_prod_975 + (l_mult_2757 + (l_mult_2141 + (l_mult_2521 + (-0.15066e5 * l_prod_972 + (0.20736e5 * l_prod_971 + (l_mult_2758 + (l_mult_2211 + (l_mult_2145 + (l_mult_1502 + (l_mult_2759 + (l_mult_2212 + (l_mult_2186 + (l_mult_2147 + (l_mult_2213 + (l_mult_2760 + (l_mult_2205 + (l_mult_2526 + (l_mult_2214 + (-0.15066e5 * l_prod_960 + (0.20736e5 * l_prod_959 + (l_mult_2761 + (l_mult_2151 + (l_mult_2217 + (l_mult_1196 + (l_mult_2762 + (l_mult_2581 + (l_mult_2154 + (l_mult_2530 + (l_mult_2582 + (l_mult_2219 + (l_mult_2580 + (l_mult_2533 + (l_mult_2583 + (l_mult_1489 + (l_mult_2763 + (l_mult_2157 + (l_mult_2256 + (l_mult_2221 + (0.7776e4 * l_prod_943 + (l_mult_2160 + (l_mult_2195 + (l_mult_2536 + (l_mult_2602 + (l_mult_2586 + (l_mult_2162 + (l_mult_2539 + (l_mult_2570 + (l_mult_2541 + (l_mult_2764 + (l_mult_2164 + l_mult_2543)))))))))))))))))))))))))))))))))))))))))))))))))))))))) + (vdot4(
            v_921,
            vcons4(
                l_mult_2229 + (0.1053e4 * l_prod_522 + (-0.62235e4 * l_prod_521 + (0.14796e5 * l_prod_976 + (-0.15390e5 * l_prod_975 + (l_mult_2765 + (l_mult_2166 + (l_mult_2544 + (0.23652e5 * l_prod_972 + (-0.37584e5 * l_prod_971 + (l_mult_2766 + (l_mult_2232 + (l_mult_2169 + (l_mult_2545 + (l_mult_2297 + (l_mult_2234 + (-0.7128e4 * l_prod_965 + (l_mult_2170 + (l_mult_2056 + (l_mult_1999 + (l_mult_2548 + (l_mult_2235 + (0.23652e5 * l_prod_960 + (-0.37584e5 * l_prod_959 + (l_mult_2767 + (l_mult_2173 + (l_mult_2237 + (-0.58320e5 * l_prod_955 + (l_mult_2350 + (l_mult_2592 + (l_mult_2176 + (l_mult_2551 + (l_mult_2519 + (l_mult_2004 + (l_mult_2553 + (l_mult_2593 + (l_mult_2238 + (l_mult_2685 + (l_mult_2179 + (l_mult_2556 + (l_mult_2239 + (l_mult_2768 + (l_mult_2008 + (l_mult_2558 + (-0.7128e4 * l_prod_941 + (l_mult_2594 + (l_mult_2044 + (l_mult_2011 + l_sum_2755))))))))))))))))))))))))))))))))))))))))))))))),
                l_mult_2248 + (-0.1016e4 * l_prod_522 + (0.6696e4 * l_prod_521 + (-0.17424e5 * l_prod_976 + (l_mult_1635 + (-0.7776e4 * l_prod_1994 + (l_mult_2183 + (l_mult_2560 + (-0.18360e5 * l_prod_972 + (0.33696e5 * l_prod_971 + (l_mult_2769 + (l_mult_2251 + (l_mult_2185 + (l_mult_2067 + (l_mult_2759 + (l_mult_2252 + (0.1296e4 * l_prod_965 + (l_mult_2187 + (l_mult_2562 + (l_mult_2253 + (-0.18360e5 * l_prod_960 + (0.33696e5 * l_prod_959 + (l_mult_2770 + (l_mult_2189 + (l_mult_2255 + (l_mult_2076 + (l_mult_2762 + (l_mult_2600 + (l_mult_2191 + (l_mult_2564 + (l_mult_2566 + (l_mult_2601 + (l_mult_2082 + (l_mult_2763 + (l_mult_2193 + (l_mult_2568 + (l_mult_2158 + (l_mult_2569 + (0.1296e4 * l_prod_941 + l_mult_2537)))))))))))))))))))))))))))))))))))))),
                l_mult_2259 + (0.594e3 * l_prod_522 + (-0.41445e4 * l_prod_521 + (0.11556e5 * l_prod_976 + (-0.13770e5 * l_prod_975 + (l_mult_2765 + (l_mult_2198 + (l_mult_2571 + (0.7128e4 * l_prod_972 + (-0.14904e5 * l_prod_971 + (l_mult_2771 + (l_mult_565 + (l_mult_2199 + (l_mult_2504 + (l_mult_1997 + (l_mult_2573 + (l_mult_2260 + (0.7128e4 * l_prod_960 + (-0.14904e5 * l_prod_959 + (l_mult_2772 + (l_mult_566 + (l_mult_2262 + (l_mult_2773 + (l_mult_2002 + l_sum_2756))))))))))))))))))))))),
                l_mult_1354 + (-0.1944e3 * l_prod_522 + (0.1404e4 * l_prod_521 + (-0.4104e4 * l_prod_976 + (l_mult_2061 + (l_mult_2757 + (l_mult_2203 + (l_mult_2576 + (-0.1134e4 * l_prod_972 + (l_mult_2774 + (l_mult_2144 + (l_mult_2578 + (l_mult_2266 + (-0.1134e4 * l_prod_960 + (l_mult_2775 + l_mult_2529)))))))))))))))) + (vdot4(
            v_922, vcons4(l_sum_2165, l_sum_2182, l_sum_2197, l_sum_2202)) + (vdot4(v_923,
            vcons4(l_sum_2206, l_sum_2591, l_sum_2599, l_sum_2604)) + (vdot4(v_924,
            vcons4(l_sum_2608, l_sum_2610, l_mult_566 + (l_mult_2510 + (l_mult_2512 + l_mult_2014)),
                l_mult_2776 + (l_mult_2516 + (l_mult_2777 + (l_mult_2052 + l_sum_2045))))) + (vdot4(v_925,
            vcons4(l_mult_2776 + (l_mult_2778 + (l_mult_2519 + l_sum_2053)), l_sum_2058,
                l_mult_2776 + (l_mult_2276 + (l_mult_2777 + (l_mult_2613 + (l_mult_2044 + l_mult_2011)))),
                l_mult_528 + (l_mult_2283 + (l_mult_2779 + (l_mult_2107 + l_sum_2780))))) + (vdot4(v_926,
            vcons4(l_mult_2776 + (l_mult_2276 + (l_mult_2778 + (l_mult_2274 + (l_mult_2519 + l_mult_2004)))),
                l_mult_2776 + (l_mult_2273 + (l_mult_2174 + l_sum_2781)),
                l_mult_2776 + (l_mult_2273 + (l_mult_2174 + (l_mult_2516 + (l_mult_2274 + l_mult_2003)))),
                l_mult_566 + (l_mult_2262 + (l_mult_2773 + l_mult_2002)))) + (vdot4(v_927,
            vcons4(
                l_mult_2420 + (l_mult_2619 + (0.19278e5 * l_prod_972 + (l_mult_1037 + (l_mult_2771 + (l_mult_2421 + (0.25704e5 * l_prod_969 + (l_mult_2621 + (l_mult_2782 + (l_mult_2422 + (-0.34992e5 * l_prod_965 + (l_mult_2783 + (l_mult_2423 + (l_mult_2784 + (l_mult_2201 + (l_mult_2785 + (l_mult_2349 + (l_mult_2786 + (l_mult_2787 + (l_mult_2788 + (l_mult_1127 + (l_mult_2629 + (l_mult_2789 + (l_mult_2425 + (l_mult_2790 + (0.6426e4 * l_prod_945 + (l_mult_2632 + (l_mult_2239 + (l_mult_2791 + (l_mult_2792 + (l_mult_2181 + (-0.5832e4 * l_prod_940 + (l_mult_2635 + (l_mult_2793 + l_mult_2014))))))))))))))))))))))))))))))))),
                l_mult_2431 + (l_mult_2658 + (-0.4860e4 * l_prod_972 + (l_mult_2774 + (l_mult_2450 + (l_mult_2794 + (l_mult_1400 + (l_mult_2759 + (l_mult_2451 + (0.42768e5 * l_prod_965 + (l_mult_2795 + (l_mult_2452 + (-0.23328e5 * l_prod_1990 + (l_mult_2188 + (l_mult_2796 + (l_mult_2388 + (l_mult_2797 + (l_mult_2798 + (l_mult_1197 + (l_mult_2530 + (l_mult_2799 + (l_mult_2454 + (l_mult_2800 + (l_mult_2801 + (l_mult_2568 + (l_mult_988 + (l_mult_2160 + (l_mult_2161 + l_sum_2803))))))))))))))))))))))))))),
                l_mult_1175 + (l_mult_2674 + (l_mult_2506 + (l_mult_2469 + (l_mult_2804 + (l_mult_2676 + (l_mult_2470 + (-0.22032e5 * l_prod_965 + (l_mult_2170 + (l_mult_2471 + (l_mult_2784 + (l_mult_2172 + (l_mult_2805 + (l_mult_2410 + (l_mult_2806 + (l_mult_2444 + (l_mult_2807 + (l_mult_2351 + (l_mult_2790 + l_sum_2053)))))))))))))))))),
                l_mult_1028 + (l_mult_2682 + (l_mult_1023 + (l_mult_2808 + (l_mult_2478 + (l_mult_2809 + (l_mult_2479 + (l_mult_2760 + (l_mult_2150 + l_sum_2811)))))))))) + (vdot4(
            v_928,
            vcons4(
                l_mult_2431 + (l_mult_2636 + (-0.28836e5 * l_prod_972 + (0.41472e5 * l_prod_971 + (l_mult_2769 + (l_mult_2433 + (l_mult_2794 + (l_mult_2638 + (-0.46656e5 * l_prod_1992 + (l_mult_2434 + (l_mult_1096 + (l_mult_2795 + (l_mult_2204 + (l_mult_2149 + (l_mult_2796 + (l_mult_2364 + (l_mult_2812 + (l_mult_2389 + (l_mult_2813 + (l_mult_1197 + (l_mult_2644 + (l_mult_2810 + (l_mult_2155 + (l_mult_2801 + (l_mult_987 + (l_mult_2221 + (l_mult_2194 + (l_mult_2160 + (l_mult_2802 + l_mult_2089)))))))))))))))))))))))))))),
                l_mult_2456 + (l_mult_2663 + (l_mult_1673 + (l_mult_2649 + (l_mult_2457 + (0.17496e5 * l_prod_969 + (l_mult_2664 + (l_mult_2297 + (l_mult_1685 + (-0.27216e5 * l_prod_965 + (l_mult_2783 + (l_mult_2458 + (l_mult_2171 + (l_mult_2814 + (l_mult_2395 + (l_mult_2773 + (l_mult_2815 + (l_mult_2460 + (l_mult_2551 + (l_mult_2057 + (l_mult_2177 + l_sum_2780)))))))))))))))))))),
                l_mult_1649 + (l_mult_2679 + (-0.648e3 * l_prod_972 + (l_mult_2475 + (l_mult_2816 + (l_mult_2656 + (l_mult_2476 + (0.12960e5 * l_prod_965 + (l_mult_2147 + (l_mult_2204 + (l_mult_2149 + (l_mult_2817 + (l_mult_2124 + (l_mult_2818 + (l_mult_2191 + (l_mult_2819 + l_mult_2080))))))))))))))),
                l_mult_1175 + (l_mult_2647 + (0.21060e5 * l_prod_972 + (-0.36288e5 * l_prod_971 + (l_mult_2766 + (l_mult_2440 + (l_mult_2804 + (l_mult_2650 + (l_mult_2782 + (l_mult_2441 + (l_mult_2665 + (l_mult_2170 + (l_mult_2805 + (l_mult_2375 + (l_mult_2820 + (l_mult_2787 + (l_mult_2821 + (l_mult_2444 + (l_mult_2654 + l_sum_2781)))))))))))))))))))) + (vdot4(
            v_929,
            vcons4(
                l_mult_1649 + (l_mult_2668 + (-0.2916e4 * l_prod_972 + (l_mult_2774 + (l_mult_2463 + (l_mult_2816 + (l_mult_2671 + (l_mult_2759 + (l_mult_2464 + (l_mult_2809 + (l_mult_2147 + (l_mult_2817 + (l_mult_2404 + (l_mult_2822 + (l_mult_2600 + (l_mult_2191 + l_mult_2564))))))))))))))),
                l_mult_1028 + (l_mult_1832 + (-0.7614e4 * l_prod_972 + (l_mult_1271 + (l_mult_2758 + (l_mult_2447 + (l_mult_2808 + (l_mult_2656 + (l_mult_2146 + (l_mult_555 + (l_mult_2381 + (l_mult_2797 + l_mult_2126))))))))))),
                l_mult_2705 + (l_mult_2347 + (0.19278e5 * l_prod_960 + (l_mult_1043 + (l_mult_2772 + (l_mult_2785 + (l_mult_2349 + (l_mult_2786 + (l_mult_2787 + (0.6426e4 * l_prod_953 + (l_mult_2303 + (l_mult_2551 + (-0.5832e4 * l_prod_950 + (l_mult_2351 + (l_mult_2005 + (l_mult_2706 + (0.25704e5 * l_prod_947 + (l_mult_2352 + (l_mult_2823 + (l_mult_2824 + (l_mult_1438 + (l_mult_2355 + (l_mult_2791 + (l_mult_2792 + (l_mult_2825 + (l_mult_2708 + (-0.34992e5 * l_prod_941 + (l_mult_2826 + (l_mult_2827 + (l_mult_2243 + (l_mult_2559 + (l_mult_2710 + (l_mult_2828 + (l_mult_2829 + l_mult_2607))))))))))))))))))))))))))))))))),
                l_mult_2712 + (l_mult_2387 + (-0.4860e4 * l_prod_960 + (l_mult_2775 + (l_mult_2796 + (l_mult_2388 + (l_mult_2797 + (l_mult_2830 + (l_mult_2191 + (l_mult_2819 + (l_mult_2714 + (l_mult_2831 + (l_mult_1447 + (l_mult_2763 + (l_mult_2832 + (l_mult_1490 + (l_mult_2221 + (l_mult_988 + (l_mult_2160 + (l_mult_2195 + (l_mult_2716 + (0.42768e5 * l_prod_941 + (l_mult_2833 + (l_mult_2834 + (l_mult_2369 + (l_mult_2540 + (l_mult_2718 + (-0.23328e5 * l_prod_1976 + (l_mult_2835 + l_mult_2603)))))))))))))))))))))))))))))) + (vdot4(
            v_930,
            vcons4(
                l_mult_1465 + (l_mult_2409 + (l_mult_2029 + (l_mult_2805 + (l_mult_2410 + (l_mult_2516 + (l_mult_2719 + (l_mult_2836 + (l_mult_2413 + (l_mult_2837 + (l_mult_2721 + (l_mult_2052 + (l_mult_2722 + (-0.22032e5 * l_prod_941 + (l_mult_2594 + (l_mult_2838 + (l_mult_2635 + (l_mult_2012 + (l_mult_2723 + (l_mult_2828 + (l_mult_2829 + l_mult_2598)))))))))))))))))))),
                l_mult_1372 + (l_mult_2419 + (l_mult_555 + (l_mult_1596 + (l_mult_2839 + (l_mult_2840 + (l_mult_2725 + (l_mult_2841 + (l_mult_2842 + (l_mult_2727 + (l_mult_2764 + (l_mult_2164 + l_mult_2590))))))))))),
                l_mult_2712 + (l_mult_2362 + (-0.28836e5 * l_prod_960 + (0.41472e5 * l_prod_959 + (l_mult_2770 + (l_mult_2796 + (l_mult_2364 + (l_mult_2812 + (l_mult_2389 + (l_mult_2830 + (l_mult_984 + (l_mult_2530 + (l_mult_2819 + (l_mult_2080 + (l_mult_2728 + (l_mult_2831 + (l_mult_2365 + (-0.46656e5 * l_prod_1983 + (l_mult_2843 + (l_mult_1490 + (l_mult_2368 + (l_mult_2194 + (l_mult_2160 + (l_mult_2729 + (l_mult_1798 + (l_mult_2833 + (l_mult_2842 + (l_mult_2587 + l_sum_2844))))))))))))))))))))))))))),
                l_mult_2731 + (l_mult_2394 + (l_mult_1875 + (l_mult_2374 + (l_mult_2814 + (l_mult_2395 + (l_mult_2773 + (l_mult_2779 + (l_mult_2107 + (l_mult_2732 + (0.17496e5 * l_prod_947 + (l_mult_2397 + (l_mult_2685 + (l_mult_2845 + (l_mult_2733 + (l_mult_2239 + (l_mult_2768 + (l_mult_2008 + (l_mult_1908 + (-0.27216e5 * l_prod_941 + (l_mult_2826 + (l_mult_2512 + (l_mult_2595 + (l_mult_2734 + l_mult_2596))))))))))))))))))))))))) + (vdot4(
            v_931,
            vcons4(
                l_mult_1856 + (l_mult_2416 + (-0.648e3 * l_prod_960 + (l_mult_2817 + (l_mult_2124 + (l_mult_2735 + (l_mult_2846 + (l_mult_2382 + (l_mult_2847 + (l_mult_2568 + (l_mult_2737 + (0.12960e5 * l_prod_941 + (l_mult_2586 + (l_mult_2802 + (l_mult_2089 + l_sum_2844)))))))))))))),
                l_mult_1465 + (l_mult_2372 + (0.21060e5 * l_prod_960 + (-0.36288e5 * l_prod_959 + (l_mult_2767 + (l_mult_2805 + (l_mult_2375 + (l_mult_2820 + (l_mult_2787 + (l_mult_2516 + (l_mult_2274 + (l_mult_2003 + (l_mult_2738 + (l_mult_2836 + (l_mult_2377 + (l_mult_2823 + (l_mult_2848 + (l_mult_2721 + (l_mult_2305 + (l_mult_2739 + (l_mult_2398 + l_mult_2594)))))))))))))))))))),
                l_mult_1856 + (l_mult_2401 + (-0.2916e4 * l_prod_960 + (l_mult_2775 + (l_mult_2817 + (l_mult_2404 + (l_mult_2822 + (l_mult_2740 + (l_mult_2846 + (l_mult_2405 + (l_mult_2763 + (l_mult_2193 + (l_mult_2568 + (l_mult_2158 + (l_mult_2741 + (l_mult_2841 + l_mult_2586))))))))))))))),
                l_mult_1372 + (l_mult_1272 + (-0.7614e4 * l_prod_960 + (l_mult_1540 + (l_mult_2761 + (l_mult_555 + (l_mult_2381 + (l_mult_2797 + (l_mult_2126 + (l_mult_2742 + (l_mult_2839 + (l_mult_2382 + l_mult_2584))))))))))))) + (vdot4(
            v_932,
            vcons4(
                l_mult_2785 + (l_mult_2301 + (-0.17496e5 * l_prod_955 + (l_mult_2002 + (l_mult_2788 + (l_mult_2303 + (l_mult_2654 + (l_mult_2789 + (l_mult_2177 + (l_mult_2790 + (l_mult_2824 + (l_mult_2632 + (l_mult_2305 + (l_mult_1391 + (l_mult_2792 + (l_mult_2557 + (l_mult_2827 + (l_mult_2595 + (l_mult_2244 + l_mult_2829)))))))))))))))))),
                l_mult_2796 + (l_mult_2314 + (l_mult_2822 + (l_mult_2813 + (l_mult_2191 + (l_mult_2810 + (l_mult_2832 + (l_mult_987 + (l_mult_2158 + (l_mult_1671 + (l_mult_2160 + (l_mult_2161 + (l_mult_2834 + (l_mult_2587 + (l_mult_2849 + l_mult_2835)))))))))))))),
                l_mult_2805 + (l_mult_2276 + (l_mult_2821 + (l_mult_2837 + (l_mult_2613 + (l_mult_2850 + (l_mult_2838 + (l_mult_2011 + (l_mult_2793 + l_mult_2829)))))))),
                l_mult_555 + (l_mult_2840 + (l_mult_2842 + l_mult_2164)))) + (vdot4(v_933,
            vcons4(
                l_mult_2796 + (l_mult_2314 + (l_mult_2822 + (l_mult_2798 + (l_mult_984 + (l_mult_2564 + (l_mult_2799 + (l_mult_2155 + (l_mult_2800 + (l_mult_2843 + (l_mult_2568 + (l_mult_1671 + (l_mult_2160 + (l_mult_2851 + l_sum_2852))))))))))))),
                l_mult_2814 + (l_mult_2283 + (l_mult_2815 + (l_mult_2107 + (l_mult_2057 + (l_mult_2845 + (l_mult_2111 + (-0.14580e5 * l_prod_943 + (l_mult_2008 + (l_mult_2181 + l_sum_2853))))))))),
                l_mult_2817 + (l_mult_2600 + (l_mult_2847 + (l_mult_2194 + l_sum_2803))),
                l_mult_2805 + (l_mult_2276 + (l_mult_2806 + (l_mult_2274 + (l_mult_2807 + (l_mult_2004 + (l_mult_2790 + (l_mult_2848 + (l_mult_2850 + l_mult_2825)))))))))) + (vdot4(
            v_934,
            vcons4(l_mult_2817 + (l_mult_2818 + (l_mult_2819 + l_sum_2196)), l_sum_2811,
                0.4320e4 * l_prod_509 + (l_mult_2481 + (0.58320e5 * l_prod_955 + (l_mult_2762 + (-0.15984e5 * l_prod_953 + (l_mult_2483 + (l_mult_2644 + (0.19440e5 * l_prod_950 + (l_mult_2454 + (l_mult_2156 + (-0.15984e5 * l_prod_945 + (l_mult_2743 + (l_mult_2368 + (0.38880e5 * l_prod_943 + (-0.93312e5 * l_prod_1981 + (l_mult_2851 + (0.19440e5 * l_prod_940 + (l_mult_2369 + (l_mult_2849 + l_mult_2589)))))))))))))))))),
                l_mult_2854 + (l_mult_2496 + (l_mult_2773 + (l_mult_2855 + (l_mult_2444 + (l_mult_2057 + (0.13284e5 * l_prod_945 + (l_mult_2747 + (l_mult_2239 + (l_mult_2856 + (l_mult_2792 + (l_mult_2181 + (l_mult_1776 + (l_mult_2243 + (l_mult_2244 + l_mult_2597)))))))))))))))) + (vdot4(
            v_935,
            vcons4(
                l_mult_1025 + (l_mult_2500 + (l_mult_2857 + (-0.4320e4 * l_prod_945 + (l_mult_2750 + (l_mult_2858 + (0.11664e5 * l_prod_940 + (l_mult_2539 + (l_mult_2163 + l_mult_2589)))))))),
                l_mult_2854 + (l_mult_2496 + (l_mult_2773 + (0.13284e5 * l_prod_953 + (l_mult_2489 + (l_mult_2551 + (l_mult_1086 + (l_mult_2425 + (l_mult_2178 + (l_mult_2859 + (l_mult_2721 + (l_mult_2856 + (l_mult_2792 + (l_mult_2557 + l_sum_2853))))))))))))),
                l_mult_1655 + (l_mult_2404 + (l_mult_2860 + (l_mult_2191 + (l_mult_2810 + (l_mult_2861 + (l_mult_2568 + (l_mult_2159 + (l_mult_2160 + (l_mult_2161 + l_sum_2852))))))))),
                l_mult_1025 + (l_mult_2500 + (-0.4320e4 * l_prod_953 + (l_mult_2494 + (0.11664e5 * l_prod_950 + (l_mult_2219 + (l_mult_2156 + (l_mult_2862 + (l_mult_2858 + l_mult_2585)))))))))) + (vdot4(
            v_936,
            vcons4(
                l_mult_2854 + (l_mult_2487 + (l_mult_1125 + (l_mult_2350 + (l_mult_2855 + (l_mult_2489 + (l_mult_2629 + (l_mult_2057 + (l_mult_2177 + (l_mult_2859 + (l_mult_2747 + (l_mult_2355 + (l_mult_2850 + (l_mult_2792 + (l_mult_2512 + l_mult_2595)))))))))))))),
                l_mult_1655 + (l_mult_2498 + (l_mult_2797 + (l_mult_2818 + (l_mult_2191 + (l_mult_2861 + (l_mult_2256 + (l_mult_2221 + (l_mult_2194 + (l_mult_2160 + (l_mult_2842 + l_mult_2587)))))))))),
                l_mult_1655 + (l_mult_2498 + (l_mult_2797 + (l_mult_2860 + (l_mult_2154 + (l_mult_2530 + (l_mult_2810 + (l_mult_2155 + (l_mult_2847 + (l_mult_2568 + (l_mult_2194 + l_mult_2160)))))))))),
                l_mult_1025 + (l_mult_2492 + (0.34992e5 * l_prod_955 + (l_mult_2762 + (l_mult_2857 + (l_mult_2494 + (l_mult_2079 + (l_mult_2862 + (l_mult_2750 + l_mult_2085)))))))))) + vdot4(
            v_916,
            vcons4(l_sum_2016, 0.e0, 0.e0,
                l_mult_2017 + (0.274e2 * l_prod_522 + (-0.2025e3 * l_prod_521 + (0.612e3 * l_prod_976 + (-0.810e3 * l_prod_975 + l_mult_1995)))))))))))))))))))))))));
        double l_sum_2866 = l_mult_1653 + l_mult_1858;
        double l_sum_2867 = 0.18e2 * l_prod_523 + (l_mult_1473 + l_sum_2866);
        double l_mult_2868 = -0.9e1 * l_prod_523;
        double l_sum_2869 = l_mult_2868 + l_mult_583;
        double l_sum_2870 = l_mult_2229 + (l_mult_1503 + (l_mult_1684 + 0.81e2 * l_prod_512));
        double l_sum_2871 = l_mult_2207 + (l_mult_1473 + (l_mult_1653 + -0.81e2 * l_prod_512));
        double l_sum_2872 = l_mult_557 + l_mult_583;
        double l_sum_2873 = l_mult_557 + l_mult_564;
        double l_sum_2874 = l_mult_2259 + (l_mult_560 + (l_mult_1684 + l_mult_583));
        double l_sum_2875 = l_mult_1357 + l_mult_1653;
        double l_sum_2876 = l_mult_2259 + (l_mult_560 + (l_mult_564 + l_mult_1887));
        double l_sum_2877 = l_mult_1357 + l_mult_1858;
        double l_mult_2878 = 0.27e2 * l_prod_523;
        double l_sum_2879 = l_mult_2878 + (l_mult_1473 + (l_mult_1028 + l_mult_1372));
        double l_sum_2880 = l_mult_557 + l_mult_560;
        double l_sum_2881 = l_mult_2259 + (l_mult_1503 + (l_mult_564 + l_mult_583));
        double l_sum_2882 = l_mult_1357 + l_mult_1473;
        double l_sum_2883 = l_mult_2878 + (l_mult_1010 + (l_mult_1653 + l_mult_1372));
        double l_sum_2884 = l_mult_2868 + l_mult_564;
        double l_sum_2885 = l_mult_2229 + (l_mult_1503 + (0.81e2 * l_prod_518 + l_mult_1887));
        double l_sum_2886 = l_mult_2207 + (l_mult_1473 + (-0.81e2 * l_prod_518 + l_mult_1858));
        double l_sum_2887 = l_mult_2878 + (l_mult_1010 + (l_mult_1028 + l_mult_1858));
        double l_sum_2888 = l_mult_2868 + l_mult_560;
        double l_sum_2889 = l_mult_2229 + (0.81e2 * l_prod_522 + (l_mult_1684 + l_mult_1887));
        double l_sum_2890 = l_mult_2207 + (-0.81e2 * l_prod_522 + l_sum_2866);
        double l_r_2891 = l_dof_load_404 * l_sum_2867;
        double l_r_2892 = l_dof_load_454 * l_mult_560;
        double l_r_2893 = l_dof_load_459 * 0.e0;
        double l_r_2894 = l_dof_load_464 * l_mult_564;
        double l_r_2895 = l_dof_load_469 * 0.e0;
        double l_r_2896 = l_dof_load_484 * 0.e0;
        double l_r_2897 = l_r_2891 + l_dof_load_409 * l_sum_2869 + l_r_595 + l_r_596 + l_r_597 + l_r_598 + l_dof_load_434 * l_mult_560 + l_r_606 + l_dof_load_444 * l_mult_564 + l_r_612 + l_r_2892 + l_r_2893 + l_r_2894 + l_r_2895 + l_dof_load_474 * l_sum_2870 + l_dof_load_479 * l_sum_2871 + l_r_2896 + l_dof_load_489 * 0.e0 + l_dof_load_494 * l_mult_1010 + l_dof_load_499 * l_mult_1028;
        double l_r_2898 = l_dof_load_464 * l_sum_2874;
        double l_r_2899 = l_dof_load_469 * l_sum_2875;
        double l_r_2900 = l_dof_load_474 * l_sum_2876;
        double l_r_2901 = l_dof_load_479 * l_sum_2877;
        double l_r_2902 = l_r_2891 + l_r_604;
        double l_r_2903 = l_r_2902 + l_r_595;
        double l_r_2904 = l_r_2903 + l_r_596;
        double l_r_2905 = l_r_2904 + l_r_597 + l_r_598;
        double l_r_2906 = l_r_2905 + l_r_605 + l_r_606 + l_dof_load_444 * l_sum_2872 + l_dof_load_449 * l_sum_2873 + l_r_2892 + l_r_2893 + l_r_2898 + l_r_2899 + l_r_2900 + l_r_2901 + l_dof_load_484 * l_mult_560 + l_dof_load_489 * l_mult_1473 + l_dof_load_494 * l_mult_1473 + l_dof_load_499 * l_sum_2879;
        double l_r_2907 = l_dof_load_454 * l_sum_2881;
        double l_r_2908 = l_dof_load_459 * l_sum_2882;
        double l_r_2909 = l_r_2905 + l_dof_load_434 * l_sum_2872 + l_dof_load_439 * l_sum_2880 + l_r_611 + l_r_612 + l_r_2907 + l_r_2908 + l_r_2894 + l_r_2895 + l_r_2900 + l_r_2901 + l_dof_load_484 * l_mult_564 + l_dof_load_489 * l_mult_1653 + l_dof_load_494 * l_sum_2883 + l_dof_load_499 * l_mult_1653;
        double l_r_2910 = l_dof_load_474 * l_mult_583;
        double l_r_2911 = l_dof_load_479 * 0.e0;
        double l_r_2912 = l_r_2902 + l_dof_load_414 * l_sum_2884 + l_r_596 + l_dof_load_424 * l_mult_560 + l_r_598 + l_r_605 + l_r_606 + l_r_611 + l_dof_load_449 * l_mult_583 + l_r_2892 + l_r_2893 + l_dof_load_464 * l_sum_2885 + l_dof_load_469 * l_sum_2886 + l_r_2910 + l_r_2911 + l_r_2896 + l_dof_load_489 * l_mult_1010 + l_dof_load_494 * 0.e0 + l_dof_load_499 * l_mult_1372;
        double l_r_2913 = l_r_2904 + l_dof_load_424 * l_sum_2873 + l_dof_load_429 * l_sum_2880 + l_r_605 + l_r_606 + l_r_611 + l_r_612 + l_r_2907 + l_r_2908 + l_r_2898 + l_r_2899 + l_r_2910 + l_r_2911 + l_dof_load_484 * l_mult_583 + l_dof_load_489 * l_sum_2887 + l_dof_load_494 * l_mult_1858 + l_dof_load_499 * l_mult_1858;
        double l_r_2914 = l_r_2903 + l_dof_load_419 * l_sum_2888 + l_r_597 + l_dof_load_429 * l_mult_564 + l_r_605 + l_dof_load_439 * l_mult_583 + l_r_611 + l_r_612 + l_dof_load_454 * l_sum_2889 + l_dof_load_459 * l_sum_2890 + l_r_2894 + l_r_2895 + l_r_2910 + l_r_2911 + l_r_2896 + l_dof_load_489 * l_mult_1028 + l_dof_load_494 * l_mult_1372 + l_dof_load_499 * 0.e0;
        double l_r_2915 = l_dof_load_405 * l_sum_2867;
        double l_r_2916 = l_dof_load_455 * l_mult_560;
        double l_r_2917 = l_dof_load_460 * 0.e0;
        double l_r_2918 = l_dof_load_465 * l_mult_564;
        double l_r_2919 = l_dof_load_470 * 0.e0;
        double l_r_2920 = l_dof_load_485 * 0.e0;
        double l_r_2921 = l_r_2915 + l_dof_load_410 * l_sum_2869 + l_r_615 + l_r_616 + l_r_617 + l_r_618 + l_dof_load_435 * l_mult_560 + l_r_626 + l_dof_load_445 * l_mult_564 + l_r_632 + l_r_2916 + l_r_2917 + l_r_2918 + l_r_2919 + l_dof_load_475 * l_sum_2870 + l_dof_load_480 * l_sum_2871 + l_r_2920 + l_dof_load_490 * 0.e0 + l_dof_load_495 * l_mult_1010 + l_dof_load_500 * l_mult_1028;
        double l_r_2922 = l_dof_load_465 * l_sum_2874;
        double l_r_2923 = l_dof_load_470 * l_sum_2875;
        double l_r_2924 = l_dof_load_475 * l_sum_2876;
        double l_r_2925 = l_dof_load_480 * l_sum_2877;
        double l_r_2926 = l_r_2915 + l_r_624;
        double l_r_2927 = l_r_2926 + l_r_615;
        double l_r_2928 = l_r_2927 + l_r_616;
        double l_r_2929 = l_r_2928 + l_r_617 + l_r_618;
        double l_r_2930 = l_r_2929 + l_r_625 + l_r_626 + l_dof_load_445 * l_sum_2872 + l_dof_load_450 * l_sum_2873 + l_r_2916 + l_r_2917 + l_r_2922 + l_r_2923 + l_r_2924 + l_r_2925 + l_dof_load_485 * l_mult_560 + l_dof_load_490 * l_mult_1473 + l_dof_load_495 * l_mult_1473 + l_dof_load_500 * l_sum_2879;
        double l_r_2931 = l_dof_load_455 * l_sum_2881;
        double l_r_2932 = l_dof_load_460 * l_sum_2882;
        double l_r_2933 = l_r_2929 + l_dof_load_435 * l_sum_2872 + l_dof_load_440 * l_sum_2880 + l_r_631 + l_r_632 + l_r_2931 + l_r_2932 + l_r_2918 + l_r_2919 + l_r_2924 + l_r_2925 + l_dof_load_485 * l_mult_564 + l_dof_load_490 * l_mult_1653 + l_dof_load_495 * l_sum_2883 + l_dof_load_500 * l_mult_1653;
        double l_r_2934 = l_dof_load_475 * l_mult_583;
        double l_r_2935 = l_dof_load_480 * 0.e0;
        double l_r_2936 = l_r_2926 + l_dof_load_415 * l_sum_2884 + l_r_616 + l_dof_load_425 * l_mult_560 + l_r_618 + l_r_625 + l_r_626 + l_r_631 + l_dof_load_450 * l_mult_583 + l_r_2916 + l_r_2917 + l_dof_load_465 * l_sum_2885 + l_dof_load_470 * l_sum_2886 + l_r_2934 + l_r_2935 + l_r_2920 + l_dof_load_490 * l_mult_1010 + l_dof_load_495 * 0.e0 + l_dof_load_500 * l_mult_1372;
        double l_r_2937 = l_r_2928 + l_dof_load_425 * l_sum_2873 + l_dof_load_430 * l_sum_2880 + l_r_625 + l_r_626 + l_r_631 + l_r_632 + l_r_2931 + l_r_2932 + l_r_2922 + l_r_2923 + l_r_2934 + l_r_2935 + l_dof_load_485 * l_mult_583 + l_dof_load_490 * l_sum_2887 + l_dof_load_495 * l_mult_1858 + l_dof_load_500 * l_mult_1858;
        double l_r_2938 = l_r_2927 + l_dof_load_420 * l_sum_2888 + l_r_617 + l_dof_load_430 * l_mult_564 + l_r_625 + l_dof_load_440 * l_mult_583 + l_r_631 + l_r_632 + l_dof_load_455 * l_sum_2889 + l_dof_load_460 * l_sum_2890 + l_r_2918 + l_r_2919 + l_r_2934 + l_r_2935 + l_r_2920 + l_dof_load_490 * l_mult_1028 + l_dof_load_495 * l_mult_1372 + l_dof_load_500 * 0.e0;
        double l_r_2939 = l_dof_load_406 * l_sum_2867;
        double l_r_2940 = l_dof_load_456 * l_mult_560;
        double l_r_2941 = l_dof_load_461 * 0.e0;
        double l_r_2942 = l_dof_load_466 * l_mult_564;
        double l_r_2943 = l_dof_load_471 * 0.e0;
        double l_r_2944 = l_dof_load_486 * 0.e0;
        double l_r_2945 = l_r_2939 + l_dof_load_411 * l_sum_2869 + l_r_635 + l_r_636 + l_r_637 + l_r_638 + l_dof_load_436 * l_mult_560 + l_r_646 + l_dof_load_446 * l_mult_564 + l_r_652 + l_r_2940 + l_r_2941 + l_r_2942 + l_r_2943 + l_dof_load_476 * l_sum_2870 + l_dof_load_481 * l_sum_2871 + l_r_2944 + l_dof_load_491 * 0.e0 + l_dof_load_496 * l_mult_1010 + l_dof_load_501 * l_mult_1028;
        double l_r_2946 = l_dof_load_466 * l_sum_2874;
        double l_r_2947 = l_dof_load_471 * l_sum_2875;
        double l_r_2948 = l_dof_load_476 * l_sum_2876;
        double l_r_2949 = l_dof_load_481 * l_sum_2877;
        double l_r_2950 = l_r_2939 + l_r_644;
        double l_r_2951 = l_r_2950 + l_r_635;
        double l_r_2952 = l_r_2951 + l_r_636;
        double l_r_2953 = l_r_2952 + l_r_637 + l_r_638;
        double l_r_2954 = l_r_2953 + l_r_645 + l_r_646 + l_dof_load_446 * l_sum_2872 + l_dof_load_451 * l_sum_2873 + l_r_2940 + l_r_2941 + l_r_2946 + l_r_2947 + l_r_2948 + l_r_2949 + l_dof_load_486 * l_mult_560 + l_dof_load_491 * l_mult_1473 + l_dof_load_496 * l_mult_1473 + l_dof_load_501 * l_sum_2879;
        double l_r_2955 = l_dof_load_456 * l_sum_2881;
        double l_r_2956 = l_dof_load_461 * l_sum_2882;
        double l_r_2957 = l_r_2953 + l_dof_load_436 * l_sum_2872 + l_dof_load_441 * l_sum_2880 + l_r_651 + l_r_652 + l_r_2955 + l_r_2956 + l_r_2942 + l_r_2943 + l_r_2948 + l_r_2949 + l_dof_load_486 * l_mult_564 + l_dof_load_491 * l_mult_1653 + l_dof_load_496 * l_sum_2883 + l_dof_load_501 * l_mult_1653;
        double l_r_2958 = l_dof_load_476 * l_mult_583;
        double l_r_2959 = l_dof_load_481 * 0.e0;
        double l_r_2960 = l_r_2950 + l_dof_load_416 * l_sum_2884 + l_r_636 + l_dof_load_426 * l_mult_560 + l_r_638 + l_r_645 + l_r_646 + l_r_651 + l_dof_load_451 * l_mult_583 + l_r_2940 + l_r_2941 + l_dof_load_466 * l_sum_2885 + l_dof_load_471 * l_sum_2886 + l_r_2958 + l_r_2959 + l_r_2944 + l_dof_load_491 * l_mult_1010 + l_dof_load_496 * 0.e0 + l_dof_load_501 * l_mult_1372;
        double l_r_2961 = l_r_2952 + l_dof_load_426 * l_sum_2873 + l_dof_load_431 * l_sum_2880 + l_r_645 + l_r_646 + l_r_651 + l_r_652 + l_r_2955 + l_r_2956 + l_r_2946 + l_r_2947 + l_r_2958 + l_r_2959 + l_dof_load_486 * l_mult_583 + l_dof_load_491 * l_sum_2887 + l_dof_load_496 * l_mult_1858 + l_dof_load_501 * l_mult_1858;
        double l_r_2962 = l_r_2951 + l_dof_load_421 * l_sum_2888 + l_r_637 + l_dof_load_431 * l_mult_564 + l_r_645 + l_dof_load_441 * l_mult_583 + l_r_651 + l_r_652 + l_dof_load_456 * l_sum_2889 + l_dof_load_461 * l_sum_2890 + l_r_2942 + l_r_2943 + l_r_2958 + l_r_2959 + l_r_2944 + l_dof_load_491 * l_mult_1028 + l_dof_load_496 * l_mult_1372 + l_dof_load_501 * 0.e0;
        double l_r_2963 = 0.e0 * l_r_2897;
        double l_r_2964 = 0.e0 * l_r_2921;
        double l_r_2965 = 0.e0 * l_r_2945;
        double l_r_2966 = l_r_2963 + l_r_2964;
        double l_r_2967 = l_r_2966 + l_r_2965;
        double l_r_2968 = 0.e0 * l_r_2906;
        double l_r_2969 = 0.e0 * l_r_2930;
        double l_r_2970 = 0.e0 * l_r_2954;
        double l_r_2971 = l_r_2968 + l_r_2969;
        double l_r_2972 = l_r_2971 + l_r_2970;
        double l_r_2973 = 0.e0 * l_r_2909;
        double l_r_2974 = 0.e0 * l_r_2933;
        double l_r_2975 = 0.e0 * l_r_2957;
        double l_r_2976 = l_r_2973 + l_r_2974;
        double l_r_2977 = l_r_2976 + l_r_2975;
        double l_r_2978 = 0.e0 * l_r_2912;
        double l_r_2979 = 0.e0 * l_r_2936;
        double l_r_2980 = 0.e0 * l_r_2960;
        double l_r_2981 = l_r_2978 + l_r_2979;
        double l_r_2982 = l_r_2981 + l_r_2980;
        double l_r_2983 = 0.e0 * l_r_2913;
        double l_r_2984 = 0.e0 * l_r_2937;
        double l_r_2985 = 0.e0 * l_r_2961;
        double l_r_2986 = l_r_2983 + l_r_2984;
        double l_r_2987 = l_r_2986 + l_r_2985;
        double l_r_2988 = 0.e0 * l_r_2914;
        double l_r_2989 = 0.e0 * l_r_2938;
        double l_r_2990 = 0.e0 * l_r_2962;
        double l_r_2991 = l_r_2988 + l_r_2989;
        double l_r_2992 = l_r_2991 + l_r_2990;
        double l_r_2993 = l_r_2966 + -0.1e1 * l_r_2945;
        double l_r_2994 = l_r_2971 + -0.1e1 * l_r_2954;
        double l_r_2995 = l_r_2976 + -0.1e1 * l_r_2957;
        double l_r_2996 = l_r_2981 + -0.1e1 * l_r_2960;
        double l_r_2997 = l_r_2986 + -0.1e1 * l_r_2961;
        double l_r_2998 = l_r_2991 + -0.1e1 * l_r_2962;
        double l_r_2999 = l_r_2963 + 0.1e1 * l_r_2921 + l_r_2965;
        double l_r_3000 = l_r_2968 + 0.1e1 * l_r_2930 + l_r_2970;
        double l_r_3001 = l_r_2973 + 0.1e1 * l_r_2933 + l_r_2975;
        double l_r_3002 = l_r_2978 + 0.1e1 * l_r_2936 + l_r_2980;
        double l_r_3003 = l_r_2983 + 0.1e1 * l_r_2937 + l_r_2985;
        double l_r_3004 = l_r_2988 + 0.1e1 * l_r_2938 + l_r_2990;
        double l_r_3005 = l_r_2966 + 0.1e1 * l_r_2945;
        double l_r_3006 = l_r_2971 + 0.1e1 * l_r_2954;
        double l_r_3007 = l_r_2976 + 0.1e1 * l_r_2957;
        double l_r_3008 = l_r_2981 + 0.1e1 * l_r_2960;
        double l_r_3009 = l_r_2986 + 0.1e1 * l_r_2961;
        double l_r_3010 = l_r_2991 + 0.1e1 * l_r_2962;
        double l_r_3011 = -0.1e1 * l_r_2897 + l_r_2964 + l_r_2965;
        double l_r_3012 = -0.1e1 * l_r_2906 + l_r_2969 + l_r_2970;
        double l_r_3013 = -0.1e1 * l_r_2909 + l_r_2974 + l_r_2975;
        double l_r_3014 = -0.1e1 * l_r_2912 + l_r_2979 + l_r_2980;
        double l_r_3015 = -0.1e1 * l_r_2913 + l_r_2984 + l_r_2985;
        double l_r_3016 = -0.1e1 * l_r_2914 + l_r_2989 + l_r_2990;
        double l_r_3017 = l_r_2963 + -0.1e1 * l_r_2921 + l_r_2965;
        double l_r_3018 = l_r_2968 + -0.1e1 * l_r_2930 + l_r_2970;
        double l_r_3019 = l_r_2973 + -0.1e1 * l_r_2933 + l_r_2975;
        double l_r_3020 = l_r_2978 + -0.1e1 * l_r_2936 + l_r_2980;
        double l_r_3021 = l_r_2983 + -0.1e1 * l_r_2937 + l_r_2985;
        double l_r_3022 = l_r_2988 + -0.1e1 * l_r_2938 + l_r_2990;
        double l_r_3023 = 0.1e1 * l_r_2897 + l_r_2964 + l_r_2965;
        double l_r_3024 = 0.1e1 * l_r_2906 + l_r_2969 + l_r_2970;
        double l_r_3025 = 0.1e1 * l_r_2909 + l_r_2974 + l_r_2975;
        double l_r_3026 = 0.1e1 * l_r_2912 + l_r_2979 + l_r_2980;
        double l_r_3027 = 0.1e1 * l_r_2913 + l_r_2984 + l_r_2985;
        double l_r_3028 = 0.1e1 * l_r_2914 + l_r_2989 + l_r_2990;
        double l_r_3029 = l_r_658 * l_r_2906 + l_r_675 * l_r_2930 + l_r_681 * l_r_2954;
        double l_r_3030 = l_r_658 * l_r_2909 + l_r_675 * l_r_2933 + l_r_681 * l_r_2957;
        double l_r_3031 = l_r_658 * l_r_2913 + l_r_675 * l_r_2937 + l_r_681 * l_r_2961;
        double l_r_3032 = l_r_663 * l_r_2906 + l_r_676 * l_r_2930 + l_r_682 * l_r_2954;
        double l_r_3033 = l_r_663 * l_r_2909 + l_r_676 * l_r_2933 + l_r_682 * l_r_2957;
        double l_r_3034 = l_r_663 * l_r_2913 + l_r_676 * l_r_2937 + l_r_682 * l_r_2961;
        double l_r_3035 = l_r_668 * l_r_2906 + l_r_677 * l_r_2930 + l_r_683 * l_r_2954;
        double l_r_3036 = l_r_668 * l_r_2909 + l_r_677 * l_r_2933 + l_r_683 * l_r_2957;
        double l_r_3037 = l_r_668 * l_r_2913 + l_r_677 * l_r_2937 + l_r_683 * l_r_2961;
        double l_r_3038 = l_r_669 * l_r_2906 + l_r_658 * l_r_2930 + l_r_684 * l_r_2954;
        double l_r_3039 = l_r_669 * l_r_2909 + l_r_658 * l_r_2933 + l_r_684 * l_r_2957;
        double l_r_3040 = l_r_669 * l_r_2913 + l_r_658 * l_r_2937 + l_r_684 * l_r_2961;
        double l_r_3041 = l_r_670 * l_r_2906 + l_r_663 * l_r_2930 + l_r_685 * l_r_2954;
        double l_r_3042 = l_r_670 * l_r_2909 + l_r_663 * l_r_2933 + l_r_685 * l_r_2957;
        double l_r_3043 = l_r_670 * l_r_2913 + l_r_663 * l_r_2937 + l_r_685 * l_r_2961;
        double l_r_3044 = l_r_671 * l_r_2906 + l_r_668 * l_r_2930 + l_r_686 * l_r_2954;
        double l_r_3045 = l_r_671 * l_r_2909 + l_r_668 * l_r_2933 + l_r_686 * l_r_2957;
        double l_r_3046 = l_r_671 * l_r_2913 + l_r_668 * l_r_2937 + l_r_686 * l_r_2961;
        double l_r_3047 = l_r_672 * l_r_2906 + l_r_678 * l_r_2930 + l_r_658 * l_r_2954;
        double l_r_3048 = l_r_672 * l_r_2909 + l_r_678 * l_r_2933 + l_r_658 * l_r_2957;
        double l_r_3049 = l_r_672 * l_r_2913 + l_r_678 * l_r_2937 + l_r_658 * l_r_2961;
        double l_r_3050 = l_r_673 * l_r_2906 + l_r_679 * l_r_2930 + l_r_663 * l_r_2954;
        double l_r_3051 = l_r_673 * l_r_2909 + l_r_679 * l_r_2933 + l_r_663 * l_r_2957;
        double l_r_3052 = l_r_673 * l_r_2913 + l_r_679 * l_r_2937 + l_r_663 * l_r_2961;
        double l_r_3053 = l_r_674 * l_r_2906 + l_r_680 * l_r_2930 + l_r_668 * l_r_2954;
        double l_r_3054 = l_r_674 * l_r_2909 + l_r_680 * l_r_2933 + l_r_668 * l_r_2957;
        double l_r_3055 = l_r_674 * l_r_2913 + l_r_680 * l_r_2937 + l_r_668 * l_r_2961;
        double l_r_3056 = l_r_603 * l_r_2972 + l_r_623 * l_r_3006 + l_r_643 * l_r_3018;
        double l_r_3057 = l_r_603 * l_r_2977 + l_r_623 * l_r_3007 + l_r_643 * l_r_3019;
        double l_r_3058 = l_r_603 * l_r_2987 + l_r_623 * l_r_3009 + l_r_643 * l_r_3021;
        double l_r_3059 = l_r_603 * l_r_2994 + l_r_623 * l_r_2972 + l_r_643 * l_r_3024;
        double l_r_3060 = l_r_603 * l_r_2995 + l_r_623 * l_r_2977 + l_r_643 * l_r_3025;
        double l_r_3061 = l_r_603 * l_r_2997 + l_r_623 * l_r_2987 + l_r_643 * l_r_3027;
        double l_r_3062 = l_r_603 * l_r_3000 + l_r_623 * l_r_3012 + l_r_643 * l_r_2972;
        double l_r_3063 = l_r_603 * l_r_3001 + l_r_623 * l_r_3013 + l_r_643 * l_r_2977;
        double l_r_3064 = l_r_603 * l_r_3003 + l_r_623 * l_r_3015 + l_r_643 * l_r_2987;
        double l_r_3065 = l_r_610 * l_r_2972 + l_r_630 * l_r_3006 + l_r_650 * l_r_3018;
        double l_r_3066 = l_r_610 * l_r_2977 + l_r_630 * l_r_3007 + l_r_650 * l_r_3019;
        double l_r_3067 = l_r_610 * l_r_2987 + l_r_630 * l_r_3009 + l_r_650 * l_r_3021;
        double l_r_3068 = l_r_610 * l_r_2994 + l_r_630 * l_r_2972 + l_r_650 * l_r_3024;
        double l_r_3069 = l_r_610 * l_r_2995 + l_r_630 * l_r_2977 + l_r_650 * l_r_3025;
        double l_r_3070 = l_r_610 * l_r_2997 + l_r_630 * l_r_2987 + l_r_650 * l_r_3027;
        double l_r_3071 = l_r_610 * l_r_3000 + l_r_630 * l_r_3012 + l_r_650 * l_r_2972;
        double l_r_3072 = l_r_610 * l_r_3001 + l_r_630 * l_r_3013 + l_r_650 * l_r_2977;
        double l_r_3073 = l_r_610 * l_r_3003 + l_r_630 * l_r_3015 + l_r_650 * l_r_2987;
        double l_r_3074 = l_r_613 * l_r_2972 + l_r_633 * l_r_3006 + l_r_653 * l_r_3018;
        double l_r_3075 = l_r_613 * l_r_2977 + l_r_633 * l_r_3007 + l_r_653 * l_r_3019;
        double l_r_3076 = l_r_613 * l_r_2987 + l_r_633 * l_r_3009 + l_r_653 * l_r_3021;
        double l_r_3077 = l_r_613 * l_r_2994 + l_r_633 * l_r_2972 + l_r_653 * l_r_3024;
        double l_r_3078 = l_r_613 * l_r_2995 + l_r_633 * l_r_2977 + l_r_653 * l_r_3025;
        double l_r_3079 = l_r_613 * l_r_2997 + l_r_633 * l_r_2987 + l_r_653 * l_r_3027;
        double l_r_3080 = l_r_613 * l_r_3000 + l_r_633 * l_r_3012 + l_r_653 * l_r_2972;
        double l_r_3081 = l_r_613 * l_r_3001 + l_r_633 * l_r_3013 + l_r_653 * l_r_2977;
        double l_r_3082 = l_r_613 * l_r_3003 + l_r_633 * l_r_3015 + l_r_653 * l_r_2987;
        vec3 v_3083 = vcons3(l_r_658 * l_r_2897 + l_r_675 * l_r_2921 + l_r_681 * l_r_2945, l_r_3029, l_r_3030) + vcons3(
            l_r_603 * l_r_2967 + l_r_623 * l_r_3005 + l_r_643 * l_r_3017, l_r_3056, l_r_3057);
        vec3 v_3084 = vcons3(l_r_3029, l_r_658 * l_r_2912 + l_r_675 * l_r_2936 + l_r_681 * l_r_2960, l_r_3031) + vcons3(
            l_r_610 * l_r_2967 + l_r_630 * l_r_3005 + l_r_650 * l_r_3017, l_r_3065, l_r_3066);
        vec3 v_3085 = vcons3(l_r_3030, l_r_3031, l_r_658 * l_r_2914 + l_r_675 * l_r_2938 + l_r_681 * l_r_2962) + vcons3(
            l_r_613 * l_r_2967 + l_r_633 * l_r_3005 + l_r_653 * l_r_3017, l_r_3074, l_r_3075);
        vec3 v_3086 = vcons3(l_r_663 * l_r_2897 + l_r_676 * l_r_2921 + l_r_682 * l_r_2945, l_r_3032, l_r_3033) + vcons3(
            l_r_3056, l_r_603 * l_r_2982 + l_r_623 * l_r_3008 + l_r_643 * l_r_3020, l_r_3058);
        vec3 v_3087 = vcons3(l_r_3032, l_r_663 * l_r_2912 + l_r_676 * l_r_2936 + l_r_682 * l_r_2960, l_r_3034) + vcons3(
            l_r_3065, l_r_610 * l_r_2982 + l_r_630 * l_r_3008 + l_r_650 * l_r_3020, l_r_3067);
        vec3 v_3088 = vcons3(l_r_3033, l_r_3034, l_r_663 * l_r_2914 + l_r_676 * l_r_2938 + l_r_682 * l_r_2962) + vcons3(
            l_r_3074, l_r_613 * l_r_2982 + l_r_633 * l_r_3008 + l_r_653 * l_r_3020, l_r_3076);
        vec3 v_3089 = vcons3(l_r_668 * l_r_2897 + l_r_677 * l_r_2921 + l_r_683 * l_r_2945, l_r_3035, l_r_3036) + vcons3(
            l_r_3057, l_r_3058, l_r_603 * l_r_2992 + l_r_623 * l_r_3010 + l_r_643 * l_r_3022);
        vec3 v_3090 = vcons3(l_r_3035, l_r_668 * l_r_2912 + l_r_677 * l_r_2936 + l_r_683 * l_r_2960, l_r_3037) + vcons3(
            l_r_3066, l_r_3067, l_r_610 * l_r_2992 + l_r_630 * l_r_3010 + l_r_650 * l_r_3022);
        vec3 v_3091 = vcons3(l_r_3036, l_r_3037, l_r_668 * l_r_2914 + l_r_677 * l_r_2938 + l_r_683 * l_r_2962) + vcons3(
            l_r_3075, l_r_3076, l_r_613 * l_r_2992 + l_r_633 * l_r_3010 + l_r_653 * l_r_3022);
        vec3 v_3092 = vcons3(l_r_669 * l_r_2897 + l_r_658 * l_r_2921 + l_r_684 * l_r_2945, l_r_3038, l_r_3039) + vcons3(
            l_r_603 * l_r_2993 + l_r_623 * l_r_2967 + l_r_643 * l_r_3023, l_r_3059, l_r_3060);
        vec3 v_3093 = vcons3(l_r_3038, l_r_669 * l_r_2912 + l_r_658 * l_r_2936 + l_r_684 * l_r_2960, l_r_3040) + vcons3(
            l_r_610 * l_r_2993 + l_r_630 * l_r_2967 + l_r_650 * l_r_3023, l_r_3068, l_r_3069);
        vec3 v_3094 = vcons3(l_r_3039, l_r_3040, l_r_669 * l_r_2914 + l_r_658 * l_r_2938 + l_r_684 * l_r_2962) + vcons3(
            l_r_613 * l_r_2993 + l_r_633 * l_r_2967 + l_r_653 * l_r_3023, l_r_3077, l_r_3078);
        vec3 v_3095 = vcons3(l_r_670 * l_r_2897 + l_r_663 * l_r_2921 + l_r_685 * l_r_2945, l_r_3041, l_r_3042) + vcons3(
            l_r_3059, l_r_603 * l_r_2996 + l_r_623 * l_r_2982 + l_r_643 * l_r_3026, l_r_3061);
        vec3 v_3096 = vcons3(l_r_3041, l_r_670 * l_r_2912 + l_r_663 * l_r_2936 + l_r_685 * l_r_2960, l_r_3043) + vcons3(
            l_r_3068, l_r_610 * l_r_2996 + l_r_630 * l_r_2982 + l_r_650 * l_r_3026, l_r_3070);
        vec3 v_3097 = vcons3(l_r_3042, l_r_3043, l_r_670 * l_r_2914 + l_r_663 * l_r_2938 + l_r_685 * l_r_2962) + vcons3(
            l_r_3077, l_r_613 * l_r_2996 + l_r_633 * l_r_2982 + l_r_653 * l_r_3026, l_r_3079);
        vec3 v_3098 = vcons3(l_r_671 * l_r_2897 + l_r_668 * l_r_2921 + l_r_686 * l_r_2945, l_r_3044, l_r_3045) + vcons3(
            l_r_3060, l_r_3061, l_r_603 * l_r_2998 + l_r_623 * l_r_2992 + l_r_643 * l_r_3028);
        vec3 v_3099 = vcons3(l_r_3044, l_r_671 * l_r_2912 + l_r_668 * l_r_2936 + l_r_686 * l_r_2960, l_r_3046) + vcons3(
            l_r_3069, l_r_3070, l_r_610 * l_r_2998 + l_r_630 * l_r_2992 + l_r_650 * l_r_3028);
        vec3 v_3100 = vcons3(l_r_3045, l_r_3046, l_r_671 * l_r_2914 + l_r_668 * l_r_2938 + l_r_686 * l_r_2962) + vcons3(
            l_r_3078, l_r_3079, l_r_613 * l_r_2998 + l_r_633 * l_r_2992 + l_r_653 * l_r_3028);
        vec3 v_3101 = vcons3(l_r_672 * l_r_2897 + l_r_678 * l_r_2921 + l_r_658 * l_r_2945, l_r_3047, l_r_3048) + vcons3(
            l_r_603 * l_r_2999 + l_r_623 * l_r_3011 + l_r_643 * l_r_2967, l_r_3062, l_r_3063);
        vec3 v_3102 = vcons3(l_r_3047, l_r_672 * l_r_2912 + l_r_678 * l_r_2936 + l_r_658 * l_r_2960, l_r_3049) + vcons3(
            l_r_610 * l_r_2999 + l_r_630 * l_r_3011 + l_r_650 * l_r_2967, l_r_3071, l_r_3072);
        vec3 v_3103 = vcons3(l_r_3048, l_r_3049, l_r_672 * l_r_2914 + l_r_678 * l_r_2938 + l_r_658 * l_r_2962) + vcons3(
            l_r_613 * l_r_2999 + l_r_633 * l_r_3011 + l_r_653 * l_r_2967, l_r_3080, l_r_3081);
        vec3 v_3104 = vcons3(l_r_673 * l_r_2897 + l_r_679 * l_r_2921 + l_r_663 * l_r_2945, l_r_3050, l_r_3051) + vcons3(
            l_r_3062, l_r_603 * l_r_3002 + l_r_623 * l_r_3014 + l_r_643 * l_r_2982, l_r_3064);
        vec3 v_3105 = vcons3(l_r_3050, l_r_673 * l_r_2912 + l_r_679 * l_r_2936 + l_r_663 * l_r_2960, l_r_3052) + vcons3(
            l_r_3071, l_r_610 * l_r_3002 + l_r_630 * l_r_3014 + l_r_650 * l_r_2982, l_r_3073);
        vec3 v_3106 = vcons3(l_r_3051, l_r_3052, l_r_673 * l_r_2914 + l_r_679 * l_r_2938 + l_r_663 * l_r_2962) + vcons3(
            l_r_3080, l_r_613 * l_r_3002 + l_r_633 * l_r_3014 + l_r_653 * l_r_2982, l_r_3082);
        vec3 v_3107 = vcons3(l_r_674 * l_r_2897 + l_r_680 * l_r_2921 + l_r_668 * l_r_2945, l_r_3053, l_r_3054) + vcons3(
            l_r_3063, l_r_3064, l_r_603 * l_r_3004 + l_r_623 * l_r_3016 + l_r_643 * l_r_2992);
        vec3 v_3108 = vcons3(l_r_3053, l_r_674 * l_r_2912 + l_r_680 * l_r_2936 + l_r_668 * l_r_2960, l_r_3055) + vcons3(
            l_r_3072, l_r_3073, l_r_610 * l_r_3004 + l_r_630 * l_r_3016 + l_r_650 * l_r_2992);
        vec3 v_3109 = vcons3(l_r_3054, l_r_3055, l_r_674 * l_r_2914 + l_r_680 * l_r_2938 + l_r_668 * l_r_2962) + vcons3(
            l_r_3081, l_r_3082, l_r_613 * l_r_3004 + l_r_633 * l_r_3016 + l_r_653 * l_r_2992);
        double l_r_3110 = 0.e0 * v_3083[0];
        double l_r_3111 = v_3086[0];
        double l_r_3112 = v_3089[0];
        double l_r_3113 = 0.e0 * l_r_3112;
        double l_r_3114 = v_3084[0];
        double l_r_3115 = 0.e0 * l_r_3114;
        double l_r_3116 = 0.e0 * v_3087[0];
        double l_r_3117 = v_3090[0];
        double l_r_3118 = v_3085[0];
        double l_r_3119 = 0.e0 * l_r_3118;
        double l_r_3120 = v_3088[0];
        double l_r_3121 = 0.e0 * v_3091[0];
        double l_r_3122 = l_r_3110 + 0.e0 * l_r_3111;
        double l_r_3123 = 0.e0 * v_3083[1];
        double l_r_3124 = v_3086[1];
        double l_r_3125 = v_3089[1];
        double l_r_3126 = 0.e0 * l_r_3125;
        double l_r_3127 = v_3084[1];
        double l_r_3128 = 0.e0 * l_r_3127;
        double l_r_3129 = 0.e0 * v_3087[1];
        double l_r_3130 = v_3090[1];
        double l_r_3131 = v_3085[1];
        double l_r_3132 = 0.e0 * l_r_3131;
        double l_r_3133 = v_3088[1];
        double l_r_3134 = 0.e0 * v_3091[1];
        double l_r_3135 = l_r_3123 + 0.e0 * l_r_3124;
        double l_r_3136 = 0.e0 * v_3083[2];
        double l_r_3137 = v_3086[2];
        double l_r_3138 = v_3089[2];
        double l_r_3139 = 0.e0 * l_r_3138;
        double l_r_3140 = v_3084[2];
        double l_r_3141 = 0.e0 * l_r_3140;
        double l_r_3142 = 0.e0 * v_3087[2];
        double l_r_3143 = v_3090[2];
        double l_r_3144 = v_3085[2];
        double l_r_3145 = 0.e0 * l_r_3144;
        double l_r_3146 = v_3088[2];
        double l_r_3147 = 0.e0 * v_3091[2];
        double l_r_3148 = l_r_3136 + 0.e0 * l_r_3137;
        double l_r_3149 = 0.e0 * v_3092[0];
        double l_r_3150 = v_3095[0];
        double l_r_3151 = v_3098[0];
        double l_r_3152 = 0.e0 * l_r_3151;
        double l_r_3153 = v_3093[0];
        double l_r_3154 = 0.e0 * l_r_3153;
        double l_r_3155 = 0.e0 * v_3096[0];
        double l_r_3156 = v_3099[0];
        double l_r_3157 = v_3094[0];
        double l_r_3158 = 0.e0 * l_r_3157;
        double l_r_3159 = v_3097[0];
        double l_r_3160 = 0.e0 * v_3100[0];
        double l_r_3161 = l_r_3149 + 0.e0 * l_r_3150;
        double l_r_3162 = 0.e0 * v_3092[1];
        double l_r_3163 = v_3095[1];
        double l_r_3164 = v_3098[1];
        double l_r_3165 = 0.e0 * l_r_3164;
        double l_r_3166 = v_3093[1];
        double l_r_3167 = 0.e0 * l_r_3166;
        double l_r_3168 = 0.e0 * v_3096[1];
        double l_r_3169 = v_3099[1];
        double l_r_3170 = v_3094[1];
        double l_r_3171 = 0.e0 * l_r_3170;
        double l_r_3172 = v_3097[1];
        double l_r_3173 = 0.e0 * v_3100[1];
        double l_r_3174 = l_r_3162 + 0.e0 * l_r_3163;
        double l_r_3175 = 0.e0 * v_3092[2];
        double l_r_3176 = v_3095[2];
        double l_r_3177 = v_3098[2];
        double l_r_3178 = 0.e0 * l_r_3177;
        double l_r_3179 = v_3093[2];
        double l_r_3180 = 0.e0 * l_r_3179;
        double l_r_3181 = 0.e0 * v_3096[2];
        double l_r_3182 = v_3099[2];
        double l_r_3183 = v_3094[2];
        double l_r_3184 = 0.e0 * l_r_3183;
        double l_r_3185 = v_3097[2];
        double l_r_3186 = 0.e0 * v_3100[2];
        double l_r_3187 = l_r_3175 + 0.e0 * l_r_3176;
        double l_r_3188 = 0.e0 * v_3101[0];
        double l_r_3189 = v_3104[0];
        double l_r_3190 = v_3107[0];
        double l_r_3191 = 0.e0 * l_r_3190;
        double l_r_3192 = v_3102[0];
        double l_r_3193 = 0.e0 * l_r_3192;
        double l_r_3194 = 0.e0 * v_3105[0];
        double l_r_3195 = v_3108[0];
        double l_r_3196 = v_3103[0];
        double l_r_3197 = 0.e0 * l_r_3196;
        double l_r_3198 = v_3106[0];
        double l_r_3199 = 0.e0 * v_3109[0];
        double l_r_3200 = l_r_3188 + 0.e0 * l_r_3189;
        double l_r_3201 = 0.e0 * v_3101[1];
        double l_r_3202 = v_3104[1];
        double l_r_3203 = v_3107[1];
        double l_r_3204 = 0.e0 * l_r_3203;
        double l_r_3205 = v_3102[1];
        double l_r_3206 = 0.e0 * l_r_3205;
        double l_r_3207 = 0.e0 * v_3105[1];
        double l_r_3208 = v_3108[1];
        double l_r_3209 = v_3103[1];
        double l_r_3210 = 0.e0 * l_r_3209;
        double l_r_3211 = v_3106[1];
        double l_r_3212 = 0.e0 * v_3109[1];
        double l_r_3213 = l_r_3201 + 0.e0 * l_r_3202;
        double l_r_3214 = 0.e0 * v_3101[2];
        double l_r_3215 = v_3104[2];
        double l_r_3216 = v_3107[2];
        double l_r_3217 = 0.e0 * l_r_3216;
        double l_r_3218 = v_3102[2];
        double l_r_3219 = 0.e0 * l_r_3218;
        double l_r_3220 = 0.e0 * v_3105[2];
        double l_r_3221 = v_3108[2];
        double l_r_3222 = v_3103[2];
        double l_r_3223 = 0.e0 * l_r_3222;
        double l_r_3224 = v_3106[2];
        double l_r_3225 = 0.e0 * v_3109[2];
        double l_r_3226 = l_r_3214 + 0.e0 * l_r_3215;
        double l_r_3227 = 0.e0 * l_r_3117;
        double l_r_3228 = 0.e0 * l_r_3120;
        double l_r_3229 = 0.e0 * l_r_3130;
        double l_r_3230 = 0.e0 * l_r_3133;
        double l_r_3231 = 0.e0 * l_r_3143;
        double l_r_3232 = 0.e0 * l_r_3146;
        double l_r_3233 = 0.e0 * l_r_3156;
        double l_r_3234 = 0.e0 * l_r_3159;
        double l_r_3235 = 0.e0 * l_r_3169;
        double l_r_3236 = 0.e0 * l_r_3172;
        double l_r_3237 = 0.e0 * l_r_3182;
        double l_r_3238 = 0.e0 * l_r_3185;
        double l_r_3239 = 0.e0 * l_r_3195;
        double l_r_3240 = 0.e0 * l_r_3198;
        double l_r_3241 = 0.e0 * l_r_3208;
        double l_r_3242 = 0.e0 * l_r_3211;
        double l_r_3243 = 0.e0 * l_r_3221;
        double l_r_3244 = 0.e0 * l_r_3224;
        vec3 v_3245 = vcons3(l_vdot_706 * l_r_2897 + l_vdot_707 * l_r_2921 + l_vdot_708 * l_r_2945,
            l_vdot_706 * l_r_2906 + l_vdot_707 * l_r_2930 + l_vdot_708 * l_r_2954,
            l_vdot_706 * l_r_2909 + l_vdot_707 * l_r_2933 + l_vdot_708 * l_r_2957) + vcons3(
            l_r_603 * l_r_3117 + l_r_623 * l_r_3156 + l_r_643 * l_r_3195,
            l_r_603 * l_r_3130 + l_r_623 * l_r_3169 + l_r_643 * l_r_3208,
            l_r_603 * l_r_3143 + l_r_623 * l_r_3182 + l_r_643 * l_r_3221);
        double l_r_3246 = v_3245[0];
        double l_r_3247 = v_3245[1];
        double l_r_3248 = v_3245[2];
        vec3 v_3249 = vscale3(l_op1_e3_l_16_709,
            vcons3(
                l_r_3122 + l_r_3113 + l_r_3115 + l_r_3116 + 0.1e1 * l_r_3117 + l_r_3119 + -0.1e1 * l_r_3120 + l_r_3121,
                l_r_3135 + l_r_3126 + l_r_3128 + l_r_3129 + 0.1e1 * l_r_3130 + l_r_3132 + -0.1e1 * l_r_3133 + l_r_3134,
                l_r_3148 + l_r_3139 + l_r_3141 + l_r_3142 + 0.1e1 * l_r_3143 + l_r_3145 + -0.1e1 * l_r_3146 + l_r_3147)) - vcons3(
            l_r_717 * l_r_3246, l_r_717 * l_r_3247, l_r_717 * l_r_3248);
        vec3 v_3250 = vscale3(l_op1_e3_l_16_709,
            vcons3(
                l_r_3161 + l_r_3152 + l_r_3154 + l_r_3155 + 0.1e1 * l_r_3156 + l_r_3158 + -0.1e1 * l_r_3159 + l_r_3160,
                l_r_3174 + l_r_3165 + l_r_3167 + l_r_3168 + 0.1e1 * l_r_3169 + l_r_3171 + -0.1e1 * l_r_3172 + l_r_3173,
                l_r_3187 + l_r_3178 + l_r_3180 + l_r_3181 + 0.1e1 * l_r_3182 + l_r_3184 + -0.1e1 * l_r_3185 + l_r_3186)) - vcons3(
            l_r_725 * l_r_3246, l_r_725 * l_r_3247, l_r_725 * l_r_3248);
        vec3 v_3251 = vscale3(l_op1_e3_l_16_709,
            vcons3(
                l_r_3200 + l_r_3191 + l_r_3193 + l_r_3194 + 0.1e1 * l_r_3195 + l_r_3197 + -0.1e1 * l_r_3198 + l_r_3199,
                l_r_3213 + l_r_3204 + l_r_3206 + l_r_3207 + 0.1e1 * l_r_3208 + l_r_3210 + -0.1e1 * l_r_3211 + l_r_3212,
                l_r_3226 + l_r_3217 + l_r_3219 + l_r_3220 + 0.1e1 * l_r_3221 + l_r_3223 + -0.1e1 * l_r_3224 + l_r_3225)) - vcons3(
            l_r_733 * l_r_3246, l_r_733 * l_r_3247, l_r_733 * l_r_3248);
        vec3 v_3252 = vscale3(l_op1_e3_l_16_709,
            vcons3(
                l_r_3122 + -0.1e1 * l_r_3112 + l_r_3115 + l_r_3116 + l_r_3227 + 0.1e1 * l_r_3118 + l_r_3228 + l_r_3121,
                l_r_3135 + -0.1e1 * l_r_3125 + l_r_3128 + l_r_3129 + l_r_3229 + 0.1e1 * l_r_3131 + l_r_3230 + l_r_3134,
                l_r_3148 + -0.1e1 * l_r_3138 + l_r_3141 + l_r_3142 + l_r_3231 + 0.1e1 * l_r_3144 + l_r_3232 + l_r_3147)) - vcons3(
            l_r_736 * l_r_3246, l_r_736 * l_r_3247, l_r_736 * l_r_3248);
        vec3 v_3253 = vscale3(l_op1_e3_l_16_709,
            vcons3(
                l_r_3161 + -0.1e1 * l_r_3151 + l_r_3154 + l_r_3155 + l_r_3233 + 0.1e1 * l_r_3157 + l_r_3234 + l_r_3160,
                l_r_3174 + -0.1e1 * l_r_3164 + l_r_3167 + l_r_3168 + l_r_3235 + 0.1e1 * l_r_3170 + l_r_3236 + l_r_3173,
                l_r_3187 + -0.1e1 * l_r_3177 + l_r_3180 + l_r_3181 + l_r_3237 + 0.1e1 * l_r_3183 + l_r_3238 + l_r_3186)) - vcons3(
            l_r_739 * l_r_3246, l_r_739 * l_r_3247, l_r_739 * l_r_3248);
        vec3 v_3254 = vscale3(l_op1_e3_l_16_709,
            vcons3(
                l_r_3200 + -0.1e1 * l_r_3190 + l_r_3193 + l_r_3194 + l_r_3239 + 0.1e1 * l_r_3196 + l_r_3240 + l_r_3199,
                l_r_3213 + -0.1e1 * l_r_3203 + l_r_3206 + l_r_3207 + l_r_3241 + 0.1e1 * l_r_3209 + l_r_3242 + l_r_3212,
                l_r_3226 + -0.1e1 * l_r_3216 + l_r_3219 + l_r_3220 + l_r_3243 + 0.1e1 * l_r_3222 + l_r_3244 + l_r_3225)) - vcons3(
            l_r_742 * l_r_3246, l_r_742 * l_r_3247, l_r_742 * l_r_3248);
        vec3 v_3255 = vscale3(l_op1_e3_l_16_709,
            vcons3(
                l_r_3110 + 0.1e1 * l_r_3111 + l_r_3113 + -0.1e1 * l_r_3114 + l_r_3116 + l_r_3227 + l_r_3119 + l_r_3228 + l_r_3121,
                l_r_3123 + 0.1e1 * l_r_3124 + l_r_3126 + -0.1e1 * l_r_3127 + l_r_3129 + l_r_3229 + l_r_3132 + l_r_3230 + l_r_3134,
                l_r_3136 + 0.1e1 * l_r_3137 + l_r_3139 + -0.1e1 * l_r_3140 + l_r_3142 + l_r_3231 + l_r_3145 + l_r_3232 + l_r_3147)) - vcons3(
            l_r_743 * l_r_3246, l_r_743 * l_r_3247, l_r_743 * l_r_3248);
        vec3 v_3256 = vscale3(l_op1_e3_l_16_709,
            vcons3(
                l_r_3149 + 0.1e1 * l_r_3150 + l_r_3152 + -0.1e1 * l_r_3153 + l_r_3155 + l_r_3233 + l_r_3158 + l_r_3234 + l_r_3160,
                l_r_3162 + 0.1e1 * l_r_3163 + l_r_3165 + -0.1e1 * l_r_3166 + l_r_3168 + l_r_3235 + l_r_3171 + l_r_3236 + l_r_3173,
                l_r_3175 + 0.1e1 * l_r_3176 + l_r_3178 + -0.1e1 * l_r_3179 + l_r_3181 + l_r_3237 + l_r_3184 + l_r_3238 + l_r_3186)) - vcons3(
            l_r_744 * l_r_3246, l_r_744 * l_r_3247, l_r_744 * l_r_3248);
        vec3 v_3257 = vscale3(l_op1_e3_l_16_709,
            vcons3(
                l_r_3188 + 0.1e1 * l_r_3189 + l_r_3191 + -0.1e1 * l_r_3192 + l_r_3194 + l_r_3239 + l_r_3197 + l_r_3240 + l_r_3199,
                l_r_3201 + 0.1e1 * l_r_3202 + l_r_3204 + -0.1e1 * l_r_3205 + l_r_3207 + l_r_3241 + l_r_3210 + l_r_3242 + l_r_3212,
                l_r_3214 + 0.1e1 * l_r_3215 + l_r_3217 + -0.1e1 * l_r_3218 + l_r_3220 + l_r_3243 + l_r_3223 + l_r_3244 + l_r_3225)) - vcons3(
            l_r_745 * l_r_3246, l_r_745 * l_r_3247, l_r_745 * l_r_3248);
        double l_op1_e3_l_90_3258 = l_op1_e3_l_18_746 * l_op1_e3_l_16_709;
        vec3 v_3259 = vcons3(l_r_1956 * l_r_1965 + l_r_1959 * l_r_1968 + l_r_1962 * l_r_1971,
            l_r_1956 * l_r_1966 + l_r_1959 * l_r_1969 + l_r_1962 * l_r_1972,
            l_r_1956 * l_r_1967 + l_r_1959 * l_r_1970 + l_r_1962 * l_r_1973) + vcons3(
            l_vdot_2863 * (v_3249[0] / l_op1_e3_l_90_3258) + l_vdot_2864 * (v_3252[0] / l_op1_e3_l_90_3258) + l_vdot_2865 * (v_3255[0] / l_op1_e3_l_90_3258),
            l_vdot_2863 * (v_3249[1] / l_op1_e3_l_90_3258) + l_vdot_2864 * (v_3252[1] / l_op1_e3_l_90_3258) + l_vdot_2865 * (v_3255[1] / l_op1_e3_l_90_3258),
            l_vdot_2863 * (v_3249[2] / l_op1_e3_l_90_3258) + l_vdot_2864 * (v_3252[2] / l_op1_e3_l_90_3258) + l_vdot_2865 * (v_3255[2] / l_op1_e3_l_90_3258));
        vec3 v_3260 = vcons3(l_r_1957 * l_r_1965 + l_r_1960 * l_r_1968 + l_r_1963 * l_r_1971,
            l_r_1957 * l_r_1966 + l_r_1960 * l_r_1969 + l_r_1963 * l_r_1972,
            l_r_1957 * l_r_1967 + l_r_1960 * l_r_1970 + l_r_1963 * l_r_1973) + vcons3(
            l_vdot_2863 * (v_3250[0] / l_op1_e3_l_90_3258) + l_vdot_2864 * (v_3253[0] / l_op1_e3_l_90_3258) + l_vdot_2865 * (v_3256[0] / l_op1_e3_l_90_3258),
            l_vdot_2863 * (v_3250[1] / l_op1_e3_l_90_3258) + l_vdot_2864 * (v_3253[1] / l_op1_e3_l_90_3258) + l_vdot_2865 * (v_3256[1] / l_op1_e3_l_90_3258),
            l_vdot_2863 * (v_3250[2] / l_op1_e3_l_90_3258) + l_vdot_2864 * (v_3253[2] / l_op1_e3_l_90_3258) + l_vdot_2865 * (v_3256[2] / l_op1_e3_l_90_3258));
        vec3 v_3261 = vcons3(l_r_1958 * l_r_1965 + l_r_1961 * l_r_1968 + l_r_1964 * l_r_1971,
            l_r_1958 * l_r_1966 + l_r_1961 * l_r_1969 + l_r_1964 * l_r_1972,
            l_r_1958 * l_r_1967 + l_r_1961 * l_r_1970 + l_r_1964 * l_r_1973) + vcons3(
            l_vdot_2863 * (v_3251[0] / l_op1_e3_l_90_3258) + l_vdot_2864 * (v_3254[0] / l_op1_e3_l_90_3258) + l_vdot_2865 * (v_3257[0] / l_op1_e3_l_90_3258),
            l_vdot_2863 * (v_3251[1] / l_op1_e3_l_90_3258) + l_vdot_2864 * (v_3254[1] / l_op1_e3_l_90_3258) + l_vdot_2865 * (v_3257[1] / l_op1_e3_l_90_3258),
            l_vdot_2863 * (v_3251[2] / l_op1_e3_l_90_3258) + l_vdot_2864 * (v_3254[2] / l_op1_e3_l_90_3258) + l_vdot_2865 * (v_3257[2] / l_op1_e3_l_90_3258));
        tensor_3_3 t_3262 = {v_3259[0],v_3259[1],v_3259[2],v_3260[0],v_3260[1],v_3260[2],v_3261[0],v_3261[1],v_3261[2],};
        eigenvecs(tensor_ref_3_3(t_3262), l_evals_3263, l__t_3264);
        pthread_mutex_lock(&wrld->_sched->_prLock);
        wrld->print() << l__t_394 << "," << std::flush;
        pthread_mutex_unlock(&wrld->_sched->_prLock);
        v_3265 = vload3(l__t_3264[2].addr(0));
    }
    else {
        pthread_mutex_lock(&wrld->_sched->_prLock);
        wrld->print() << "Error at input pos:" << tensor_ref_3(self->sv_xp) << "\n" << std::flush;
        pthread_mutex_unlock(&wrld->_sched->_prLock);
        v_3265 = vload3(tensor_ref_3(self->sv_normal).addr(0));
    }
    vpack3(self->sv_normal, v_3265);
    return diderot::kStabilize;
}
extern "C" bool evalProg_output_get_normal (evalProg_world_t *cWrld, Nrrd *nData)
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
        memcpy(cp, &wrld->_strands.strand(ix)->sv_normal, 3 * sizeof(double));
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
    diderot::dynseq< tensor_3 > seq_2 = glob->gv_ipos;
    int32_t base[1] = {0,};
    uint32_t size[1] = {static_cast<uint32_t>(seq_2.length()),};
    if (this->alloc(base, size)) {
        return true;
    }
    uint32_t ix = 0;
    for (auto it_3 = seq_2.cbegin(); it_3 != seq_2.cend(); ++it_3) {
        auto i_x_3267 = *it_3;
        mesh_pos_mesh_t l__t_3268 = fn_findPos(glob->gv_meshData, i_x_3267);
        normal_init(this->_strands.strand(ix), l__t_3268, i_x_3267);
        ++ix;
    }
    this->swap_state();
    this->_stage = diderot::POST_CREATE;
    return false;
}
/*---------- begin par-worker-nobsp.in ----------*/
struct CACHE_ALIGN worker_arg {
    world       *_wrld;         //!< world pointer
    uint32_t    _id;            //!< worker ID
    uint32_t    _maxNSteps;     //!< maximum number of steps to take; 0 == infinity
    uint32_t    _nSteps;        //!< max number of steps taken by a strand in call to run
    uint32_t    _nStable;       //!< number of strands that stabilized in call to run
    uint32_t    _nDead;         //!< number of strands that died in call to run
    worker_cache _strands;
};

/* Worker task for when we do not need super-step synchronization */
static void worker (void *arg)
{
    worker_arg *myArg = reinterpret_cast<worker_arg *>(arg);
    world *wrld = myArg->_wrld;
#ifndef DIDEROT_NO_GLOBALS
    globals *glob = wrld->_globals;
#endif

  // iterate until there is no more work to do
    uint32_t numDead = 0;
    uint32_t numStabilized = 0;
    uint32_t maxSteps = 0;
    uint32_t maxNSteps = myArg->_maxNSteps;
    strand_array::sched_block *blk;
    IF_LOGGING ( LogGetStrandBlock(wrld, myArg->_id+1); )
    while ((blk = myArg->_strands.get_block()) != nullptr) {
        IF_LOGGING ( LogGotStrandBlock(wrld, myArg->_id+1); )
        uint32_t nStable = blk->_nStable;
#ifdef DIDEROT_HAS_STRAND_DIE
        uint32_t nDead = blk->_nDead;
#endif
      // update the strands
        for (auto ix = myArg->_strands.begin_active(blk);
            ix != myArg->_strands.end_active(blk);
        ) {
          // run the strand to completion, or until the step limit is exceeded
            normal_strand *self = myArg->_strands.strand(ix);
            diderot::strand_status sts = myArg->_strands.status(ix);
#ifdef DIDEROT_HAS_START_METHOD
            if (sts == diderot::kNew) {
                IF_LOGGING ( LogStrandStart(wrld, myArg->_id+1, ix); )
                sts = normal_start(self);
            }
#endif
            uint32_t nSteps = 0;
            while ((! sts) && (nSteps < maxNSteps)) {
                nSteps++;
                sts = normal_update(wrld, glob, self);
            }
            switch (sts) {
              case diderot::kStabilize:
              // stabilize the strand's state.
                IF_LOGGING ( LogStrandStabilize(wrld, myArg->_id+1, ix); )
                ix = myArg->_strands.strand_stabilize (blk, ix);
                break;
#ifdef DIDEROT_HAS_STRAND_DIE
              case diderot::kDie:
                IF_LOGGING ( LogStrandDie(wrld, myArg->_id+1, ix); )
                ix = myArg->_strands.kill (blk, ix);
                break;
#endif
              default:
                assert (sts == myArg->_strands.status(ix));
                ix = myArg->_strands.next_active(blk, ix);
                break;
            }
            if (maxSteps < nSteps) maxSteps = nSteps;
        }
        numStabilized += (blk->_nStable - nStable);
#ifdef DIDEROT_HAS_STRAND_DIE
        numDead += (blk->_nDead - nDead);
#endif
        IF_LOGGING ( LogGetStrandBlock(wrld, myArg->_id+1); )
    }
    IF_LOGGING ( LogNoStrandBlock(wrld, myArg->_id+1); )

  // update global counts of active and stable strands
    myArg->_nSteps = maxSteps;
    myArg->_nStable = numStabilized;
    myArg->_nDead = numDead;

}
/*---------- end par-worker-nobsp.in ----------*/
/*---------- begin par-run.in ----------*/
//! Run the Diderot program (parallel version)
//! \param maxNSteps the limit on the number of super steps; 0 means unlimited
//! \return the number of steps taken, or 0 if done or there is an error.
uint32_t world::run (uint32_t maxNSteps)
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

    diderot::scheduler *sched = this->_sched;

    if (maxNSteps == 0) {
        maxNSteps = 0xffffffff;  // essentially unlimited
    }

  // set task pointer
    sched->_task = worker;

  // initialize per-worker info
    this->_strands._workers.clear();
    worker_arg *args = new worker_arg[sched->_numWorkers];
    for (int i = 0;  i < sched->_numWorkers;  i++) {
        worker_arg *p = &args[i];
        p->_wrld = this;
        p->_id = i;
        p->_maxNSteps = maxNSteps;
        p->_nSteps = 0;
#ifndef DIDEROT_BSP
        p->_nStable = 0;
        p->_nDead = 0;
#endif
        p->_strands.init (this->_strands);
        sched->_info[i]._data = p;
    }

    double t0 = airTime();

  // Start worker threads
    if (this->_verbose) {
        std::cerr << "run with " << this->_strands.num_active() << " active strands / "
            << sched->_numWorkers << " workers ..." << std::endl;
    }
    this->_strands.prepare_run ();
    sched->_gate.release_workers (IF_LOGGING( this ));

  // wait for the computation to finish
    sched->_gate.controller_wait (IF_LOGGING( this ));

  // get max # steps and update global counts of active and stable strands when no-bsp
    uint32_t nSteps = 0;
    for (uint32_t i = 0;  i < sched->_numWorkers;  i++) {
        nSteps = std::max (nSteps, args[i]._nSteps);
#ifndef DIDEROT_BSP
      // if there is no BSP, then the controller updates #active and #stable
        this->_strands._nActive -= args[i]._nStable + args[i]._nDead;
        this->_strands._nStable += args[i]._nStable;
#endif
    }
    delete[] args;

    t0 = airTime() - t0;
    if (this->_verbose) {
        std::cerr << nSteps << " steps done in " << t0 << " seconds" << std::endl;
    }
    this->_run_time += t0;

    if (this->_strands.num_active() == 0) {
        this->_stage = diderot::DONE;
    }

    return nSteps;

} // world::run
/*---------- end par-run.in ----------*/

/*---------- begin namespace-close.in ----------*/

} // namespace evalProg
/*---------- end namespace-close.in ----------*/

/*---------- begin c-wrappers.in ----------*/
extern "C" uint32_t evalProg_num_strands (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
    return w->_strands.num_alive();
}

extern "C" uint32_t evalProg_num_active_strands (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
    return w->_strands.num_active();
}

extern "C" uint32_t evalProg_num_stable_strands (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
    return w->_strands.num_stable();
}

extern "C" bool evalProg_any_errors (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
    return (w->_errors->errNum > 0);
}

extern "C" char *evalProg_get_errors (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
    char *msg = biffMsgStrGet (w->_errors);
    biffMsgClear (w->_errors);
    return msg;
}

extern "C" evalProg_world_t *evalProg_new_world ()
{
    evalProg::world *w = new (std::nothrow) evalProg::world();
    return reinterpret_cast<evalProg_world_t *>(w);
}

extern "C" bool evalProg_init_world (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);

    if (w->_stage != diderot::POST_NEW) {
        w->error ("multiple calls to evalProg_init_world");
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

extern "C" bool evalProg_create_strands (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);

    if (w->_stage < diderot::POST_INIT) {
        w->error ("must call evalProg_init_world before evalProg_create_strands");
        return true;
    }
    else if (w->_stage > diderot::POST_INIT) {
        w->error ("multiple calls to evalProg_create_strands");
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

extern "C" uint32_t evalProg_run (evalProg_world_t *wrld, uint32_t maxNSteps)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);

    if (w->_stage < diderot::POST_CREATE) {
        w->error ("attempt to run uninitialized program");
        return 0;
    }
    else if (w->_stage == diderot::DONE) {
        return 0;
    }

    return w->run(maxNSteps);
}

extern "C" void evalProg_shutdown (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
    delete w;
}

extern "C" void evalProg_set_verbose (evalProg_world_t *wrld, bool mode)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
    w->_verbose = (mode ? true : false);
}

extern "C" bool evalProg_get_verbose (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
    return static_cast<bool>(w->_verbose);
}

extern "C" bool evalProg_set_printer_cb (evalProg_world_t *wrld, bool (*pr)(void *, char *), void *data)
{
  /* FIXME: implement printer callback */
    return true;
}

#ifdef DIDEROT_TARGET_PARALLEL

extern "C" uint32_t evalProg_get_num_cores (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
    return w->_sched->_numHWCores;
}

extern "C" bool evalProg_set_num_workers (evalProg_world_t *wrld, uint32_t nw)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
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

extern "C" uint32_t evalProg_get_num_workers (evalProg_world_t *wrld)
{
    evalProg::world *w = reinterpret_cast<evalProg::world *>(wrld);
    return w->_sched->_numWorkers;
}

#endif /* DIDEROT_TARGET_PARALLEL */
/*---------- end c-wrappers.in ----------*/


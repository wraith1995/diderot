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
    typedef double vec1 __attribute__ ((vector_size (8)));
    typedef double vec4 __attribute__ ((vector_size (32)));
    typedef double vec3 __attribute__ ((vector_size (32)));
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
    bool gv_0space0392_intermedateGlobal;
    bool gv_0data0394_intermedateGlobal;
    bool gv_ipos;
} defined_inputs;
struct globals {
    mesh_t gv_meshData;
    fns_t gv_0space0392_intermedateGlobal;
    func_t gv_0data0394_intermedateGlobal;
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
    wrld->_definedInp.gv_0space0392_intermedateGlobal = true;
    std::memcpy(&wrld->_globals->gv_0space0392_intermedateGlobal, v, sizeof(fns_t));
    return false;
}
extern "C" bool evalProg_input_set_data (evalProg_world_t *cWrld, void *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_0data0394_intermedateGlobal = true;
    std::memcpy(&wrld->_globals->gv_0data0394_intermedateGlobal, v, sizeof(func_t));
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
    if (!wrld->_definedInp.gv_0space0392_intermedateGlobal) {
        biffMsgAdd(wrld->_errors, "undefined input \"space\"\n");
        return true;
    }
    if (!wrld->_definedInp.gv_0data0394_intermedateGlobal) {
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
    wrld->_definedInp.gv_0space0392_intermedateGlobal = false;
    wrld->_definedInp.gv_0data0394_intermedateGlobal = false;
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
    glob->gv_data = glob->gv_0data0394_intermedateGlobal.loadFem(
        glob->gv_0space0392_intermedateGlobal.loadFem(glob->gv_meshData));
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
    vec3 v_1888;
    if (self->sv_pos0.valid) {
        mesh_cell_mesh_t l__t_394 = makeFem(self->sv_pos0.mesh, self->sv_pos0.cell);
        func_cell_func_t l__t_395 = makeFem(glob->gv_data, l__t_394.cell);
        tensor_ref_3 l_evalPoint_396 = self->sv_pos0.refPos;
        func_t l__t_397 = (l__t_395.func);
        fns_t l__t_398 = l__t_397.space;
        int32_t l__t_399 = l__t_395.cell;
        mesh_t l__t_400 = l__t_398.mesh;
        int32_t l_mulRes_401 = l__t_399 * 84;
        int32_t t_402 = l__t_398.indexMap[l_mulRes_401];
        int32_t t_403 = l__t_398.indexMap[l_mulRes_401 + 1];
        int32_t t_404 = l__t_398.indexMap[l_mulRes_401 + 2];
        int32_t t_405 = l__t_398.indexMap[l_mulRes_401 + 3];
        int32_t t_406 = l__t_398.indexMap[l_mulRes_401 + 4];
        int32_t t_407 = l__t_398.indexMap[l_mulRes_401 + 5];
        int32_t t_408 = l__t_398.indexMap[l_mulRes_401 + 6];
        int32_t t_409 = l__t_398.indexMap[l_mulRes_401 + 7];
        int32_t t_410 = l__t_398.indexMap[l_mulRes_401 + 8];
        int32_t t_411 = l__t_398.indexMap[l_mulRes_401 + 9];
        int32_t t_412 = l__t_398.indexMap[l_mulRes_401 + 10];
        int32_t t_413 = l__t_398.indexMap[l_mulRes_401 + 11];
        int32_t t_414 = l__t_398.indexMap[l_mulRes_401 + 12];
        int32_t t_415 = l__t_398.indexMap[l_mulRes_401 + 13];
        int32_t t_416 = l__t_398.indexMap[l_mulRes_401 + 14];
        int32_t t_417 = l__t_398.indexMap[l_mulRes_401 + 15];
        int32_t t_418 = l__t_398.indexMap[l_mulRes_401 + 16];
        int32_t t_419 = l__t_398.indexMap[l_mulRes_401 + 17];
        int32_t t_420 = l__t_398.indexMap[l_mulRes_401 + 18];
        int32_t t_421 = l__t_398.indexMap[l_mulRes_401 + 19];
        int32_t t_422 = l__t_398.indexMap[l_mulRes_401 + 20];
        int32_t t_423 = l__t_398.indexMap[l_mulRes_401 + 21];
        int32_t t_424 = l__t_398.indexMap[l_mulRes_401 + 22];
        int32_t t_425 = l__t_398.indexMap[l_mulRes_401 + 23];
        int32_t t_426 = l__t_398.indexMap[l_mulRes_401 + 24];
        int32_t t_427 = l__t_398.indexMap[l_mulRes_401 + 25];
        int32_t t_428 = l__t_398.indexMap[l_mulRes_401 + 26];
        int32_t t_429 = l__t_398.indexMap[l_mulRes_401 + 27];
        int32_t t_430 = l__t_398.indexMap[l_mulRes_401 + 28];
        int32_t t_431 = l__t_398.indexMap[l_mulRes_401 + 29];
        int32_t t_432 = l__t_398.indexMap[l_mulRes_401 + 30];
        int32_t t_433 = l__t_398.indexMap[l_mulRes_401 + 31];
        int32_t t_434 = l__t_398.indexMap[l_mulRes_401 + 32];
        int32_t t_435 = l__t_398.indexMap[l_mulRes_401 + 33];
        int32_t t_436 = l__t_398.indexMap[l_mulRes_401 + 34];
        int32_t t_437 = l__t_398.indexMap[l_mulRes_401 + 35];
        int32_t t_438 = l__t_398.indexMap[l_mulRes_401 + 36];
        int32_t t_439 = l__t_398.indexMap[l_mulRes_401 + 37];
        int32_t t_440 = l__t_398.indexMap[l_mulRes_401 + 38];
        int32_t t_441 = l__t_398.indexMap[l_mulRes_401 + 39];
        int32_t t_442 = l__t_398.indexMap[l_mulRes_401 + 40];
        int32_t t_443 = l__t_398.indexMap[l_mulRes_401 + 41];
        int32_t t_444 = l__t_398.indexMap[l_mulRes_401 + 42];
        int32_t t_445 = l__t_398.indexMap[l_mulRes_401 + 43];
        int32_t t_446 = l__t_398.indexMap[l_mulRes_401 + 44];
        int32_t t_447 = l__t_398.indexMap[l_mulRes_401 + 45];
        int32_t t_448 = l__t_398.indexMap[l_mulRes_401 + 46];
        int32_t t_449 = l__t_398.indexMap[l_mulRes_401 + 47];
        int32_t t_450 = l__t_398.indexMap[l_mulRes_401 + 48];
        int32_t t_451 = l__t_398.indexMap[l_mulRes_401 + 49];
        int32_t t_452 = l__t_398.indexMap[l_mulRes_401 + 50];
        int32_t t_453 = l__t_398.indexMap[l_mulRes_401 + 51];
        int32_t t_454 = l__t_398.indexMap[l_mulRes_401 + 52];
        int32_t t_455 = l__t_398.indexMap[l_mulRes_401 + 53];
        int32_t t_456 = l__t_398.indexMap[l_mulRes_401 + 54];
        int32_t t_457 = l__t_398.indexMap[l_mulRes_401 + 55];
        int32_t t_458 = l__t_398.indexMap[l_mulRes_401 + 56];
        int32_t t_459 = l__t_398.indexMap[l_mulRes_401 + 57];
        int32_t t_460 = l__t_398.indexMap[l_mulRes_401 + 58];
        int32_t t_461 = l__t_398.indexMap[l_mulRes_401 + 59];
        int32_t t_462 = l__t_398.indexMap[l_mulRes_401 + 60];
        int32_t t_463 = l__t_398.indexMap[l_mulRes_401 + 61];
        int32_t t_464 = l__t_398.indexMap[l_mulRes_401 + 62];
        int32_t t_465 = l__t_398.indexMap[l_mulRes_401 + 63];
        int32_t t_466 = l__t_398.indexMap[l_mulRes_401 + 64];
        int32_t t_467 = l__t_398.indexMap[l_mulRes_401 + 65];
        int32_t t_468 = l__t_398.indexMap[l_mulRes_401 + 66];
        int32_t t_469 = l__t_398.indexMap[l_mulRes_401 + 67];
        int32_t t_470 = l__t_398.indexMap[l_mulRes_401 + 68];
        int32_t t_471 = l__t_398.indexMap[l_mulRes_401 + 69];
        int32_t t_472 = l__t_398.indexMap[l_mulRes_401 + 70];
        int32_t t_473 = l__t_398.indexMap[l_mulRes_401 + 71];
        int32_t t_474 = l__t_398.indexMap[l_mulRes_401 + 72];
        int32_t t_475 = l__t_398.indexMap[l_mulRes_401 + 73];
        int32_t t_476 = l__t_398.indexMap[l_mulRes_401 + 74];
        int32_t t_477 = l__t_398.indexMap[l_mulRes_401 + 75];
        int32_t t_478 = l__t_398.indexMap[l_mulRes_401 + 76];
        int32_t t_479 = l__t_398.indexMap[l_mulRes_401 + 77];
        int32_t t_480 = l__t_398.indexMap[l_mulRes_401 + 78];
        int32_t t_481 = l__t_398.indexMap[l_mulRes_401 + 79];
        int32_t t_482 = l__t_398.indexMap[l_mulRes_401 + 80];
        int32_t t_483 = l__t_398.indexMap[l_mulRes_401 + 81];
        int32_t t_484 = l__t_398.indexMap[l_mulRes_401 + 82];
        int32_t t_485 = l__t_398.indexMap[l_mulRes_401 + 83];
        double t_486 = l__t_397.coordMap[1 * t_485];
        double t_487 = l__t_397.coordMap[1 * t_484];
        double t_488 = l__t_397.coordMap[1 * t_483];
        double t_489 = l__t_397.coordMap[1 * t_482];
        double t_490 = l__t_397.coordMap[1 * t_481];
        double t_491 = l__t_397.coordMap[1 * t_480];
        double t_492 = l__t_397.coordMap[1 * t_479];
        double t_493 = l__t_397.coordMap[1 * t_478];
        double t_494 = l__t_397.coordMap[1 * t_477];
        double t_495 = l__t_397.coordMap[1 * t_476];
        double t_496 = l__t_397.coordMap[1 * t_475];
        double t_497 = l__t_397.coordMap[1 * t_474];
        double t_498 = l__t_397.coordMap[1 * t_473];
        double t_499 = l__t_397.coordMap[1 * t_472];
        double t_500 = l__t_397.coordMap[1 * t_471];
        double t_501 = l__t_397.coordMap[1 * t_470];
        double t_502 = l__t_397.coordMap[1 * t_469];
        double t_503 = l__t_397.coordMap[1 * t_468];
        double t_504 = l__t_397.coordMap[1 * t_467];
        double t_505 = l__t_397.coordMap[1 * t_466];
        double t_506 = l__t_397.coordMap[1 * t_465];
        double t_507 = l__t_397.coordMap[1 * t_464];
        double t_508 = l__t_397.coordMap[1 * t_463];
        double t_509 = l__t_397.coordMap[1 * t_462];
        double t_510 = l__t_397.coordMap[1 * t_461];
        double t_511 = l__t_397.coordMap[1 * t_460];
        double t_512 = l__t_397.coordMap[1 * t_459];
        double t_513 = l__t_397.coordMap[1 * t_458];
        double t_514 = l__t_397.coordMap[1 * t_457];
        double t_515 = l__t_397.coordMap[1 * t_456];
        double t_516 = l__t_397.coordMap[1 * t_455];
        double t_517 = l__t_397.coordMap[1 * t_454];
        double t_518 = l__t_397.coordMap[1 * t_453];
        double t_519 = l__t_397.coordMap[1 * t_452];
        double t_520 = l__t_397.coordMap[1 * t_451];
        double t_521 = l__t_397.coordMap[1 * t_450];
        double t_522 = l__t_397.coordMap[1 * t_449];
        double t_523 = l__t_397.coordMap[1 * t_448];
        double t_524 = l__t_397.coordMap[1 * t_447];
        double t_525 = l__t_397.coordMap[1 * t_446];
        double t_526 = l__t_397.coordMap[1 * t_445];
        double t_527 = l__t_397.coordMap[1 * t_444];
        double t_528 = l__t_397.coordMap[1 * t_443];
        double t_529 = l__t_397.coordMap[1 * t_442];
        double t_530 = l__t_397.coordMap[1 * t_441];
        double t_531 = l__t_397.coordMap[1 * t_440];
        double t_532 = l__t_397.coordMap[1 * t_439];
        double t_533 = l__t_397.coordMap[1 * t_438];
        double t_534 = l__t_397.coordMap[1 * t_437];
        double t_535 = l__t_397.coordMap[1 * t_436];
        double t_536 = l__t_397.coordMap[1 * t_435];
        double t_537 = l__t_397.coordMap[1 * t_434];
        double t_538 = l__t_397.coordMap[1 * t_433];
        double t_539 = l__t_397.coordMap[1 * t_432];
        double t_540 = l__t_397.coordMap[1 * t_431];
        double t_541 = l__t_397.coordMap[1 * t_430];
        double t_542 = l__t_397.coordMap[1 * t_429];
        double t_543 = l__t_397.coordMap[1 * t_428];
        double t_544 = l__t_397.coordMap[1 * t_427];
        double t_545 = l__t_397.coordMap[1 * t_426];
        double t_546 = l__t_397.coordMap[1 * t_425];
        double t_547 = l__t_397.coordMap[1 * t_424];
        double t_548 = l__t_397.coordMap[1 * t_423];
        double t_549 = l__t_397.coordMap[1 * t_422];
        double t_550 = l__t_397.coordMap[1 * t_421];
        double t_551 = l__t_397.coordMap[1 * t_420];
        double t_552 = l__t_397.coordMap[1 * t_419];
        double t_553 = l__t_397.coordMap[1 * t_418];
        double t_554 = l__t_397.coordMap[1 * t_417];
        double t_555 = l__t_397.coordMap[1 * t_416];
        double t_556 = l__t_397.coordMap[1 * t_415];
        double t_557 = l__t_397.coordMap[1 * t_414];
        double t_558 = l__t_397.coordMap[1 * t_413];
        double t_559 = l__t_397.coordMap[1 * t_412];
        double t_560 = l__t_397.coordMap[1 * t_411];
        double t_561 = l__t_397.coordMap[1 * t_410];
        double t_562 = l__t_397.coordMap[1 * t_409];
        double t_563 = l__t_397.coordMap[1 * t_408];
        double t_564 = l__t_397.coordMap[1 * t_407];
        double t_565 = l__t_397.coordMap[1 * t_406];
        double t_566 = l__t_397.coordMap[1 * t_405];
        double t_567 = l__t_397.coordMap[1 * t_404];
        double t_568 = l__t_397.coordMap[1 * t_403];
        double t_569 = l__t_397.coordMap[1 * t_402];
        vec4 v_570 = vcons4(t_569, t_568, t_567, t_566);
        vec4 v_571 = vcons4(t_565, t_564, t_563, t_562);
        vec4 v_572 = vcons4(t_561, t_560, t_559, t_558);
        vec4 v_573 = vcons4(t_557, t_556, t_555, t_554);
        vec4 v_574 = vcons4(t_553, t_552, t_551, t_550);
        vec4 v_575 = vcons4(t_549, t_548, t_547, t_546);
        vec4 v_576 = vcons4(t_545, t_544, t_543, t_542);
        vec4 v_577 = vcons4(t_541, t_540, t_539, t_538);
        vec4 v_578 = vcons4(t_537, t_536, t_535, t_534);
        vec4 v_579 = vcons4(t_533, t_532, t_531, t_530);
        vec4 v_580 = vcons4(t_529, t_528, t_527, t_526);
        vec4 v_581 = vcons4(t_525, t_524, t_523, t_522);
        vec4 v_582 = vcons4(t_521, t_520, t_519, t_518);
        vec4 v_583 = vcons4(t_517, t_516, t_515, t_514);
        vec4 v_584 = vcons4(t_513, t_512, t_511, t_510);
        vec4 v_585 = vcons4(t_509, t_508, t_507, t_506);
        vec4 v_586 = vcons4(t_505, t_504, t_503, t_502);
        vec4 v_587 = vcons4(t_501, t_500, t_499, t_498);
        vec4 v_588 = vcons4(t_497, t_496, t_495, t_494);
        vec4 v_589 = vcons4(t_493, t_492, t_491, t_490);
        vec4 v_590 = vcons4(t_489, t_488, t_487, t_486);
        double l_varAcc_591 = l_evalPoint_396[0];
        double l_varAcc_592 = l_evalPoint_396[1];
        double l_varAcc_593 = l_evalPoint_396[2];
        double l_prod2_594 = l_varAcc_591 * l_varAcc_591;
        double l_prod3_595 = l_prod2_594 * l_varAcc_591;
        double l_prod4_596 = l_prod3_595 * l_varAcc_591;
        double l_prod_597 = 0.1e1 * 0.1e1;
        double l_prod_598 = l_prod4_596 * l_varAcc_591 * l_prod_597;
        double l_prod_599 = l_varAcc_592 * 0.1e1;
        double l_prod_600 = l_prod4_596 * l_prod_599;
        double l_prod_601 = 0.1e1 * l_varAcc_593;
        double l_prod_602 = l_prod4_596 * l_prod_601;
        double l_prod_603 = l_prod4_596 * l_prod_597;
        double l_prod2_604 = l_varAcc_592 * l_varAcc_592;
        double l_prod_605 = l_prod2_604 * 0.1e1;
        double l_prod_606 = l_prod3_595 * l_prod_605;
        double l_prod_607 = l_varAcc_592 * l_varAcc_593;
        double l_prod_608 = l_prod3_595 * l_prod_607;
        double l_prod_609 = l_prod3_595 * l_prod_599;
        double l_prod2_610 = l_varAcc_593 * l_varAcc_593;
        double l_prod_611 = 0.1e1 * l_prod2_610;
        double l_prod_612 = l_prod3_595 * l_prod_611;
        double l_prod_613 = l_prod3_595 * l_prod_601;
        double l_prod_614 = l_prod3_595 * l_prod_597;
        double l_prod3_615 = l_prod2_604 * l_varAcc_592;
        double l_prod_616 = l_prod3_615 * 0.1e1;
        double l_prod_617 = l_prod2_594 * l_prod_616;
        double l_prod_618 = l_prod2_604 * l_varAcc_593;
        double l_prod_619 = l_prod2_594 * l_prod_618;
        double l_prod_620 = l_prod2_594 * l_prod_605;
        double l_prod_621 = l_varAcc_592 * l_prod2_610;
        double l_prod_622 = l_prod2_594 * l_prod_621;
        double l_prod_623 = l_prod2_594 * l_prod_607;
        double l_prod_624 = l_prod2_594 * l_prod_599;
        double l_prod3_625 = l_prod2_610 * l_varAcc_593;
        double l_prod_626 = 0.1e1 * l_prod3_625;
        double l_prod_627 = l_prod2_594 * l_prod_626;
        double l_prod_628 = l_prod2_594 * l_prod_611;
        double l_prod_629 = l_prod2_594 * l_prod_601;
        double l_prod_630 = l_prod2_594 * l_prod_597;
        double l_prod4_631 = l_prod3_615 * l_varAcc_592;
        double l_prod_632 = l_prod4_631 * 0.1e1;
        double l_prod_633 = l_varAcc_591 * l_prod_632;
        double l_prod_634 = l_prod3_615 * l_varAcc_593;
        double l_prod_635 = l_varAcc_591 * l_prod_634;
        double l_prod_636 = l_varAcc_591 * l_prod_616;
        double l_prod_637 = l_prod2_604 * l_prod2_610;
        double l_prod_638 = l_varAcc_591 * l_prod_637;
        double l_prod_639 = l_varAcc_591 * l_prod_618;
        double l_prod_640 = l_varAcc_591 * l_prod_605;
        double l_prod_641 = l_varAcc_592 * l_prod3_625;
        double l_prod_642 = l_varAcc_591 * l_prod_641;
        double l_prod_643 = l_varAcc_591 * l_prod_621;
        double l_prod_644 = l_varAcc_591 * l_prod_607;
        double l_prod_645 = l_varAcc_591 * l_prod_599;
        double l_prod4_646 = l_prod3_625 * l_varAcc_593;
        double l_prod_647 = 0.1e1 * l_prod4_646;
        double l_prod_648 = l_varAcc_591 * l_prod_647;
        double l_prod_649 = l_varAcc_591 * l_prod_626;
        double l_prod_650 = l_varAcc_591 * l_prod_611;
        double l_prod_651 = l_varAcc_591 * l_prod_601;
        double l_prod_652 = l_varAcc_591 * l_prod_597;
        double l_prod_653 = 0.1e1 * (l_prod4_631 * l_varAcc_592 * 0.1e1);
        double l_prod_654 = 0.1e1 * (l_prod4_631 * l_varAcc_593);
        double l_prod_655 = 0.1e1 * l_prod_632;
        double l_prod_656 = 0.1e1 * (l_prod3_615 * l_prod2_610);
        double l_prod_657 = 0.1e1 * l_prod_634;
        double l_prod_658 = 0.1e1 * l_prod_616;
        double l_prod_659 = 0.1e1 * (l_prod2_604 * l_prod3_625);
        double l_prod_660 = 0.1e1 * l_prod_637;
        double l_prod_661 = 0.1e1 * l_prod_618;
        double l_prod_662 = 0.1e1 * l_prod_605;
        double l_prod_663 = 0.1e1 * (l_varAcc_592 * l_prod4_646);
        double l_prod_664 = 0.1e1 * l_prod_641;
        double l_prod_665 = 0.1e1 * l_prod_621;
        double l_prod_666 = 0.1e1 * l_prod_607;
        double l_prod_667 = 0.1e1 * l_prod_599;
        double l_prod_668 = 0.1e1 * (0.1e1 * (l_prod4_646 * l_varAcc_593));
        double l_prod_669 = 0.1e1 * l_prod_647;
        double l_prod_670 = 0.1e1 * l_prod_626;
        double l_prod_671 = 0.1e1 * l_prod_611;
        double l_prod_672 = 0.1e1 * l_prod_601;
        double l_prod_673 = 0.1e1 * l_prod_597;
        double l_mult_674 = 0.3888e3 * l_prod_668;
        double l_mult_675 = 0.1944e4 * l_prod_663;
        double l_mult_676 = 0.3888e4 * l_prod_659;
        double l_mult_677 = 0.3888e4 * l_prod_656;
        double l_mult_678 = 0.1944e4 * l_prod_654;
        double l_mult_679 = 0.3888e3 * l_prod_653;
        double l_mult_680 = 0.1944e4 * l_prod_648;
        double l_mult_681 = 0.7776e4 * l_prod_642;
        double l_mult_682 = 0.11664e5 * l_prod_638;
        double l_mult_683 = 0.7776e4 * l_prod_635;
        double l_mult_684 = 0.1944e4 * l_prod_633;
        double l_mult_685 = 0.3888e4 * l_prod_627;
        double l_mult_686 = 0.11664e5 * l_prod_622;
        double l_mult_687 = 0.11664e5 * l_prod_619;
        double l_mult_688 = 0.3888e4 * l_prod_617;
        double l_mult_689 = 0.3888e4 * l_prod_612;
        double l_mult_690 = 0.7776e4 * l_prod_608;
        double l_mult_691 = 0.3888e4 * l_prod_606;
        double l_mult_692 = 0.1944e4 * l_prod_602;
        double l_mult_693 = 0.1944e4 * l_prod_600;
        double l_mult_694 = 0.3888e3 * l_prod_598;
        double l_sum_695 = -0.147e2 * l_prod_673 + (0.1624e3 * l_prod_672 + (-0.6615e3 * l_prod_671 + (0.1260e4 * l_prod_670 + (-0.1134e4 * l_prod_669 + (l_mult_674 + (0.1624e3 * l_prod_667 + (-0.1323e4 * l_prod_666 + (0.3780e4 * l_prod_665 + (-0.4536e4 * l_prod_664 + (l_mult_675 + (-0.6615e3 * l_prod_662 + (0.3780e4 * l_prod_661 + (-0.6804e4 * l_prod_660 + (l_mult_676 + (0.1260e4 * l_prod_658 + (-0.4536e4 * l_prod_657 + (l_mult_677 + (-0.1134e4 * l_prod_655 + (l_mult_678 + (l_mult_679 + (0.1624e3 * l_prod_652 + (-0.1323e4 * l_prod_651 + (0.3780e4 * l_prod_650 + (-0.4536e4 * l_prod_649 + (l_mult_680 + (-0.1323e4 * l_prod_645 + (0.7560e4 * l_prod_644 + (-0.13608e5 * l_prod_643 + (l_mult_681 + (0.3780e4 * l_prod_640 + (-0.13608e5 * l_prod_639 + (l_mult_682 + (-0.4536e4 * l_prod_636 + (l_mult_683 + (l_mult_684 + (-0.6615e3 * l_prod_630 + (0.3780e4 * l_prod_629 + (-0.6804e4 * l_prod_628 + (l_mult_685 + (0.3780e4 * l_prod_624 + (-0.13608e5 * l_prod_623 + (l_mult_686 + (-0.6804e4 * l_prod_620 + (l_mult_687 + (l_mult_688 + (0.1260e4 * l_prod_614 + (-0.4536e4 * l_prod_613 + (l_mult_689 + (-0.4536e4 * l_prod_609 + (l_mult_690 + (l_mult_691 + (-0.1134e4 * l_prod_603 + (l_mult_692 + (l_mult_693 + l_mult_694))))))))))))))))))))))))))))))))))))))))))))))))))))));
        double l_mult_696 = -0.1e1 * l_prod_673;
        double l_mult_697 = 0.72e1 * l_prod_672;
        double l_mult_698 = -0.180e3 * l_prod_651;
        double l_mult_699 = 0.45e1 * l_prod_672;
        double l_mult_700 = -0.27e2 * l_prod_671;
        double l_mult_701 = -0.99e2 * l_prod_651;
        double l_mult_702 = 0.594e3 * l_prod_650;
        double l_mult_703 = -0.2916e4 * l_prod_628;
        double l_mult_704 = -0.648e3 * l_prod_613;
        double l_sum_705 = l_mult_704 + l_mult_689;
        double l_mult_706 = 0.4e1 * l_prod_672;
        double l_mult_707 = -0.36e2 * l_prod_671;
        double l_mult_708 = 0.72e2 * l_prod_670;
        double l_mult_709 = -0.72e2 * l_prod_651;
        double l_mult_710 = 0.648e3 * l_prod_650;
        double l_mult_711 = 0.216e3 * l_prod_629;
        double l_mult_712 = -0.1944e4 * l_prod_628;
        double l_sum_713 = l_mult_711 + (l_mult_712 + l_mult_685);
        double l_mult_714 = -0.495e2 * l_prod_671;
        double l_mult_715 = 0.162e3 * l_prod_670;
        double l_mult_716 = -0.162e3 * l_prod_669;
        double l_mult_717 = -0.54e2 * l_prod_651;
        double l_mult_718 = -0.1944e4 * l_prod_649;
        double l_sum_719 = l_mult_717 + (l_mult_702 + (l_mult_718 + l_mult_680));
        double l_sum_720 = l_mult_697 + (-0.90e2 * l_prod_671 + (0.378e3 * l_prod_670 + (-0.648e3 * l_prod_669 + l_mult_674)));
        double l_mult_721 = 0.72e1 * l_prod_667;
        double l_mult_722 = -0.180e3 * l_prod_645;
        double l_mult_723 = 0.45e1 * l_prod_667;
        double l_mult_724 = -0.27e2 * l_prod_662;
        double l_mult_725 = -0.99e2 * l_prod_645;
        double l_mult_726 = 0.594e3 * l_prod_640;
        double l_mult_727 = -0.2916e4 * l_prod_620;
        double l_mult_728 = -0.648e3 * l_prod_609;
        double l_sum_729 = l_mult_728 + l_mult_691;
        double l_mult_730 = 0.4e1 * l_prod_667;
        double l_mult_731 = -0.36e2 * l_prod_662;
        double l_mult_732 = 0.72e2 * l_prod_658;
        double l_mult_733 = -0.72e2 * l_prod_645;
        double l_mult_734 = 0.648e3 * l_prod_640;
        double l_mult_735 = 0.216e3 * l_prod_624;
        double l_mult_736 = -0.1944e4 * l_prod_620;
        double l_sum_737 = l_mult_735 + (l_mult_736 + l_mult_688);
        double l_mult_738 = -0.495e2 * l_prod_662;
        double l_mult_739 = 0.162e3 * l_prod_658;
        double l_mult_740 = -0.162e3 * l_prod_655;
        double l_mult_741 = -0.54e2 * l_prod_645;
        double l_mult_742 = -0.1944e4 * l_prod_636;
        double l_sum_743 = l_mult_741 + (l_mult_726 + (l_mult_742 + l_mult_684));
        double l_sum_744 = l_mult_721 + (-0.90e2 * l_prod_662 + (0.378e3 * l_prod_658 + (-0.648e3 * l_prod_655 + l_mult_679)));
        double l_mult_745 = -0.3132e3 * l_prod_672;
        double l_mult_746 = 0.5184e4 * l_prod_669;
        double l_mult_747 = -0.1944e4 * l_prod_668;
        double l_mult_748 = 0.2088e4 * l_prod_666;
        double l_mult_749 = -0.10044e5 * l_prod_665;
        double l_mult_750 = 0.15552e5 * l_prod_664;
        double l_mult_751 = -0.7776e4 * l_prod_663;
        double l_mult_752 = -0.5022e4 * l_prod_661;
        double l_mult_753 = 0.15552e5 * l_prod_660;
        double l_mult_754 = -0.11664e5 * l_prod_659;
        double l_mult_755 = 0.5184e4 * l_prod_657;
        double l_mult_756 = -0.7776e4 * l_prod_656;
        double l_mult_757 = -0.1944e4 * l_prod_654;
        double l_mult_758 = 0.2088e4 * l_prod_651;
        double l_mult_759 = -0.10044e5 * l_prod_650;
        double l_mult_760 = 0.15552e5 * l_prod_649;
        double l_mult_761 = -0.7776e4 * l_prod_648;
        double l_mult_762 = -0.10044e5 * l_prod_644;
        double l_mult_763 = 0.31104e5 * l_prod_643;
        double l_mult_764 = -0.23328e5 * l_prod_642;
        double l_mult_765 = 0.15552e5 * l_prod_639;
        double l_mult_766 = -0.23328e5 * l_prod_638;
        double l_mult_767 = -0.7776e4 * l_prod_635;
        double l_mult_768 = -0.5022e4 * l_prod_629;
        double l_mult_769 = 0.15552e5 * l_prod_628;
        double l_mult_770 = -0.11664e5 * l_prod_627;
        double l_mult_771 = 0.15552e5 * l_prod_623;
        double l_mult_772 = -0.23328e5 * l_prod_622;
        double l_mult_773 = -0.11664e5 * l_prod_619;
        double l_mult_774 = 0.5184e4 * l_prod_613;
        double l_mult_775 = -0.7776e4 * l_prod_612;
        double l_mult_776 = -0.7776e4 * l_prod_608;
        double l_mult_777 = -0.1944e4 * l_prod_602;
        double l_sum_778 = l_mult_745 + (0.2088e4 * l_prod_671 + (-0.5022e4 * l_prod_670 + (l_mult_746 + (l_mult_747 + (l_mult_748 + (l_mult_749 + (l_mult_750 + (l_mult_751 + (l_mult_752 + (l_mult_753 + (l_mult_754 + (l_mult_755 + (l_mult_756 + (l_mult_757 + (l_mult_758 + (l_mult_759 + (l_mult_760 + (l_mult_761 + (l_mult_762 + (l_mult_763 + (l_mult_764 + (l_mult_765 + (l_mult_766 + (l_mult_767 + (l_mult_768 + (l_mult_769 + (l_mult_770 + (l_mult_771 + (l_mult_772 + (l_mult_773 + (l_mult_774 + (l_mult_775 + (l_mult_776 + l_mult_777)))))))))))))))))))))))))))))))));
        double l_mult_779 = 0.2565e3 * l_prod_672;
        double l_mult_780 = 0.3888e4 * l_prod_668;
        double l_mult_781 = -0.1071e4 * l_prod_666;
        double l_mult_782 = 0.9342e4 * l_prod_665;
        double l_mult_783 = 0.11664e5 * l_prod_663;
        double l_mult_784 = 0.1458e4 * l_prod_661;
        double l_mult_785 = -0.10692e5 * l_prod_660;
        double l_mult_786 = 0.11664e5 * l_prod_659;
        double l_mult_787 = -0.648e3 * l_prod_657;
        double l_mult_788 = -0.1071e4 * l_prod_651;
        double l_mult_789 = 0.9342e4 * l_prod_650;
        double l_mult_790 = 0.11664e5 * l_prod_648;
        double l_mult_791 = 0.2916e4 * l_prod_644;
        double l_mult_792 = -0.21384e5 * l_prod_643;
        double l_mult_793 = 0.23328e5 * l_prod_642;
        double l_mult_794 = -0.1944e4 * l_prod_639;
        double l_mult_795 = 0.1458e4 * l_prod_629;
        double l_mult_796 = -0.10692e5 * l_prod_628;
        double l_mult_797 = 0.11664e5 * l_prod_627;
        double l_mult_798 = -0.1944e4 * l_prod_623;
        double l_sum_799 = l_mult_779 + (-0.2610e4 * l_prod_671 + (0.7884e4 * l_prod_670 + (-0.9396e4 * l_prod_669 + (l_mult_780 + (l_mult_781 + (l_mult_782 + (-0.19440e5 * l_prod_664 + (l_mult_783 + (l_mult_784 + (l_mult_785 + (l_mult_786 + (l_mult_787 + (l_mult_677 + (l_mult_788 + (l_mult_789 + (-0.19440e5 * l_prod_649 + (l_mult_790 + (l_mult_791 + (l_mult_792 + (l_mult_793 + (l_mult_794 + (l_mult_682 + (l_mult_795 + (l_mult_796 + (l_mult_797 + (l_mult_798 + (l_mult_686 + l_sum_705)))))))))))))))))))))))))));
        double l_mult_800 = -0.148e3 * l_prod_672;
        double l_mult_801 = -0.3888e4 * l_prod_668;
        double l_mult_802 = 0.360e3 * l_prod_666;
        double l_mult_803 = -0.3672e4 * l_prod_665;
        double l_mult_804 = 0.10368e5 * l_prod_664;
        double l_mult_805 = -0.216e3 * l_prod_661;
        double l_mult_806 = 0.1944e4 * l_prod_660;
        double l_mult_807 = -0.3888e4 * l_prod_659;
        double l_mult_808 = 0.360e3 * l_prod_651;
        double l_mult_809 = -0.3672e4 * l_prod_650;
        double l_mult_810 = 0.10368e5 * l_prod_649;
        double l_mult_811 = -0.432e3 * l_prod_644;
        double l_mult_812 = 0.3888e4 * l_prod_643;
        double l_mult_813 = -0.7776e4 * l_prod_642;
        double l_mult_814 = -0.216e3 * l_prod_629;
        double l_mult_815 = 0.1944e4 * l_prod_628;
        double l_mult_816 = -0.3888e4 * l_prod_627;
        double l_sum_817 = l_mult_814 + (l_mult_815 + l_mult_816);
        double l_sum_818 = l_mult_800 + (0.1692e4 * l_prod_671 + (-0.6120e4 * l_prod_670 + (0.8424e4 * l_prod_669 + (l_mult_801 + (l_mult_802 + (l_mult_803 + (l_mult_804 + (l_mult_751 + (l_mult_805 + (l_mult_806 + (l_mult_807 + (l_mult_808 + (l_mult_809 + (l_mult_810 + (l_mult_761 + (l_mult_811 + (l_mult_812 + (l_mult_813 + l_sum_817))))))))))))))))));
        double l_mult_819 = 0.495e2 * l_prod_672;
        double l_mult_820 = 0.1944e4 * l_prod_668;
        double l_mult_821 = -0.54e2 * l_prod_666;
        double l_mult_822 = 0.594e3 * l_prod_665;
        double l_mult_823 = -0.1944e4 * l_prod_664;
        double l_sum_824 = l_mult_819 + (-0.5985e3 * l_prod_671 + (0.2376e4 * l_prod_670 + (-0.3726e4 * l_prod_669 + (l_mult_820 + (l_mult_821 + (l_mult_822 + (l_mult_823 + (l_mult_675 + l_sum_719))))))));
        double l_mult_825 = -0.72e1 * l_prod_672;
        double l_mult_826 = 0.648e3 * l_prod_669;
        double l_mult_827 = -0.3888e3 * l_prod_668;
        double l_sum_828 = l_mult_825 + (0.90e2 * l_prod_671 + (-0.378e3 * l_prod_670 + (l_mult_826 + l_mult_827)));
        double l_mult_829 = -0.3132e3 * l_prod_667;
        double l_mult_830 = -0.5022e4 * l_prod_665;
        double l_mult_831 = 0.5184e4 * l_prod_664;
        double l_mult_832 = -0.1944e4 * l_prod_663;
        double l_mult_833 = -0.10044e5 * l_prod_661;
        double l_mult_834 = -0.7776e4 * l_prod_659;
        double l_mult_835 = 0.15552e5 * l_prod_657;
        double l_mult_836 = -0.11664e5 * l_prod_656;
        double l_mult_837 = 0.5184e4 * l_prod_655;
        double l_mult_838 = -0.7776e4 * l_prod_654;
        double l_mult_839 = -0.1944e4 * l_prod_653;
        double l_mult_840 = 0.2088e4 * l_prod_645;
        double l_mult_841 = 0.15552e5 * l_prod_643;
        double l_mult_842 = -0.10044e5 * l_prod_640;
        double l_mult_843 = 0.31104e5 * l_prod_639;
        double l_mult_844 = 0.15552e5 * l_prod_636;
        double l_mult_845 = -0.23328e5 * l_prod_635;
        double l_mult_846 = -0.7776e4 * l_prod_633;
        double l_mult_847 = -0.5022e4 * l_prod_624;
        double l_mult_848 = -0.11664e5 * l_prod_622;
        double l_mult_849 = 0.15552e5 * l_prod_620;
        double l_mult_850 = -0.23328e5 * l_prod_619;
        double l_mult_851 = -0.11664e5 * l_prod_617;
        double l_mult_852 = 0.5184e4 * l_prod_609;
        double l_mult_853 = -0.7776e4 * l_prod_606;
        double l_mult_854 = -0.1944e4 * l_prod_600;
        double l_sum_855 = l_mult_829 + (l_mult_748 + (l_mult_830 + (l_mult_831 + (l_mult_832 + (0.2088e4 * l_prod_662 + (l_mult_833 + (l_mult_753 + (l_mult_834 + (-0.5022e4 * l_prod_658 + (l_mult_835 + (l_mult_836 + (l_mult_837 + (l_mult_838 + (l_mult_839 + (l_mult_840 + (l_mult_762 + (l_mult_841 + (l_mult_813 + (l_mult_842 + (l_mult_843 + (l_mult_766 + (l_mult_844 + (l_mult_845 + (l_mult_846 + (l_mult_847 + (l_mult_771 + (l_mult_848 + (l_mult_849 + (l_mult_850 + (l_mult_851 + (l_mult_852 + (l_mult_776 + (l_mult_853 + l_mult_854)))))))))))))))))))))))))))))))));
        double l_mult_856 = 0.2565e3 * l_prod_667;
        double l_mult_857 = 0.1458e4 * l_prod_665;
        double l_mult_858 = -0.648e3 * l_prod_664;
        double l_mult_859 = 0.9342e4 * l_prod_661;
        double l_mult_860 = 0.11664e5 * l_prod_656;
        double l_mult_861 = 0.11664e5 * l_prod_654;
        double l_mult_862 = 0.3888e4 * l_prod_653;
        double l_mult_863 = -0.1071e4 * l_prod_645;
        double l_mult_864 = -0.1944e4 * l_prod_643;
        double l_mult_865 = 0.9342e4 * l_prod_640;
        double l_mult_866 = -0.21384e5 * l_prod_639;
        double l_mult_867 = 0.23328e5 * l_prod_635;
        double l_mult_868 = 0.11664e5 * l_prod_633;
        double l_mult_869 = 0.1458e4 * l_prod_624;
        double l_mult_870 = -0.10692e5 * l_prod_620;
        double l_mult_871 = 0.11664e5 * l_prod_617;
        double l_sum_872 = l_mult_856 + (l_mult_781 + (l_mult_857 + (l_mult_858 + (-0.2610e4 * l_prod_662 + (l_mult_859 + (l_mult_785 + (l_mult_676 + (0.7884e4 * l_prod_658 + (-0.19440e5 * l_prod_657 + (l_mult_860 + (-0.9396e4 * l_prod_655 + (l_mult_861 + (l_mult_862 + (l_mult_863 + (l_mult_791 + (l_mult_864 + (l_mult_865 + (l_mult_866 + (l_mult_682 + (-0.19440e5 * l_prod_636 + (l_mult_867 + (l_mult_868 + (l_mult_869 + (l_mult_798 + (l_mult_870 + (l_mult_687 + (l_mult_871 + l_sum_729)))))))))))))))))))))))))));
        double l_mult_873 = -0.148e3 * l_prod_667;
        double l_mult_874 = -0.216e3 * l_prod_665;
        double l_mult_875 = -0.3672e4 * l_prod_661;
        double l_mult_876 = 0.10368e5 * l_prod_657;
        double l_mult_877 = -0.3888e4 * l_prod_656;
        double l_mult_878 = -0.3888e4 * l_prod_653;
        double l_mult_879 = 0.360e3 * l_prod_645;
        double l_mult_880 = -0.3672e4 * l_prod_640;
        double l_mult_881 = 0.3888e4 * l_prod_639;
        double l_mult_882 = 0.10368e5 * l_prod_636;
        double l_mult_883 = -0.216e3 * l_prod_624;
        double l_mult_884 = 0.1944e4 * l_prod_620;
        double l_mult_885 = -0.3888e4 * l_prod_617;
        double l_sum_886 = l_mult_883 + (l_mult_884 + l_mult_885);
        double l_sum_887 = l_mult_873 + (l_mult_802 + (l_mult_874 + (0.1692e4 * l_prod_662 + (l_mult_875 + (l_mult_806 + (-0.6120e4 * l_prod_658 + (l_mult_876 + (l_mult_877 + (0.8424e4 * l_prod_655 + (l_mult_838 + (l_mult_878 + (l_mult_879 + (l_mult_811 + (l_mult_880 + (l_mult_881 + (l_mult_882 + (l_mult_767 + (l_mult_846 + l_sum_886))))))))))))))))));
        double l_mult_888 = 0.495e2 * l_prod_667;
        double l_mult_889 = 0.594e3 * l_prod_661;
        double l_mult_890 = -0.1944e4 * l_prod_657;
        double l_mult_891 = 0.1944e4 * l_prod_653;
        double l_sum_892 = l_mult_888 + (l_mult_821 + (-0.5985e3 * l_prod_662 + (l_mult_889 + (0.2376e4 * l_prod_658 + (l_mult_890 + (-0.3726e4 * l_prod_655 + (l_mult_678 + (l_mult_891 + l_sum_743))))))));
        double l_mult_893 = -0.72e1 * l_prod_667;
        double l_mult_894 = 0.648e3 * l_prod_655;
        double l_mult_895 = -0.3888e3 * l_prod_653;
        double l_sum_896 = l_mult_893 + (0.90e2 * l_prod_662 + (-0.378e3 * l_prod_658 + (l_mult_894 + l_mult_895)));
        double l_mult_897 = 0.36e2 * l_prod_673;
        double l_mult_898 = 0.1044e4 * l_prod_671;
        double l_mult_899 = -0.1674e4 * l_prod_670;
        double l_mult_900 = 0.1296e4 * l_prod_669;
        double l_mult_901 = 0.1044e4 * l_prod_662;
        double l_mult_902 = -0.1674e4 * l_prod_658;
        double l_mult_903 = 0.1296e4 * l_prod_655;
        double l_mult_904 = 0.4176e4 * l_prod_651;
        double l_mult_905 = -0.3888e4 * l_prod_648;
        double l_mult_906 = 0.4176e4 * l_prod_645;
        double l_mult_907 = -0.20088e5 * l_prod_644;
        double l_mult_908 = -0.15552e5 * l_prod_642;
        double l_mult_909 = -0.15552e5 * l_prod_635;
        double l_mult_910 = -0.3888e4 * l_prod_633;
        double l_mult_911 = 0.23328e5 * l_prod_628;
        double l_mult_912 = 0.46656e5 * l_prod_623;
        double l_mult_913 = -0.34992e5 * l_prod_622;
        double l_mult_914 = 0.23328e5 * l_prod_620;
        double l_mult_915 = -0.34992e5 * l_prod_619;
        double l_mult_916 = -0.15552e5 * l_prod_612;
        double l_mult_917 = -0.31104e5 * l_prod_608;
        double l_mult_918 = -0.15552e5 * l_prod_606;
        double l_mult_919 = -0.9720e4 * l_prod_602;
        double l_mult_920 = -0.9720e4 * l_prod_600;
        double l_mult_921 = -0.23328e4 * l_prod_598;
        double l_mult_922 = -0.45e2 * l_prod_673;
        double l_mult_923 = -0.5355e3 * l_prod_671;
        double l_mult_924 = 0.486e3 * l_prod_670;
        double l_mult_925 = -0.5355e3 * l_prod_662;
        double l_mult_926 = -0.972e3 * l_prod_660;
        double l_mult_927 = 0.486e3 * l_prod_658;
        double l_mult_928 = -0.5220e4 * l_prod_651;
        double l_mult_929 = -0.5220e4 * l_prod_645;
        double l_mult_930 = 0.18684e5 * l_prod_644;
        double l_mult_931 = -0.29160e5 * l_prod_628;
        double l_mult_932 = 0.34992e5 * l_prod_622;
        double l_mult_933 = -0.29160e5 * l_prod_620;
        double l_mult_934 = 0.34992e5 * l_prod_619;
        double l_mult_935 = 0.23328e5 * l_prod_612;
        double l_mult_936 = 0.46656e5 * l_prod_608;
        double l_mult_937 = 0.23328e5 * l_prod_606;
        double l_mult_938 = 0.19440e5 * l_prod_602;
        double l_mult_939 = 0.19440e5 * l_prod_600;
        double l_mult_940 = 0.5832e4 * l_prod_598;
        double l_mult_941 = 0.40e2 * l_prod_673;
        double l_mult_942 = 0.180e3 * l_prod_671;
        double l_mult_943 = -0.72e2 * l_prod_670;
        double l_mult_944 = 0.180e3 * l_prod_662;
        double l_mult_945 = -0.72e2 * l_prod_658;
        double l_mult_946 = 0.3384e4 * l_prod_651;
        double l_mult_947 = 0.3384e4 * l_prod_645;
        double l_mult_948 = -0.7344e4 * l_prod_644;
        double l_mult_949 = 0.31104e5 * l_prod_623;
        double l_mult_950 = -0.19440e5 * l_prod_602;
        double l_mult_951 = -0.19440e5 * l_prod_600;
        double l_mult_952 = -0.225e2 * l_prod_673;
        double l_mult_953 = -0.1197e4 * l_prod_651;
        double l_mult_954 = -0.1197e4 * l_prod_645;
        double l_mult_955 = 0.1188e4 * l_prod_644;
        double l_mult_956 = -0.5832e4 * l_prod_623;
        double l_mult_957 = 0.9720e4 * l_prod_602;
        double l_mult_958 = 0.9720e4 * l_prod_600;
        double l_mult_959 = 0.72e1 * l_prod_673;
        double l_mult_960 = 0.180e3 * l_prod_651;
        double l_mult_961 = 0.180e3 * l_prod_645;
        double l_mult_962 = 0.2592e4 * l_prod_613;
        double l_mult_963 = 0.2592e4 * l_prod_609;
        double l_mult_964 = 0.5184e4 * l_prod_603;
        double l_mult_965 = -0.36e2 * l_prod_666;
        double l_mult_966 = 0.216e3 * l_prod_661;
        double l_mult_967 = 0.648e3 * l_prod_644;
        double l_mult_968 = -0.3888e4 * l_prod_639;
        double l_mult_969 = 0.324e3 * l_prod_661;
        double l_mult_970 = 0.432e3 * l_prod_644;
        double l_sum_971 = l_mult_970 + (l_mult_968 + l_mult_683);
        double l_mult_972 = 0.216e3 * l_prod_665;
        double l_mult_973 = -0.3888e4 * l_prod_643;
        double l_sum_974 = l_mult_798 + l_mult_686;
        double l_mult_975 = -0.27e2 * l_prod_666;
        double l_mult_976 = 0.162e3 * l_prod_665;
        double l_mult_977 = 0.162e3 * l_prod_661;
        double l_mult_978 = 0.324e3 * l_prod_644;
        double l_sum_979 = l_mult_794 + l_mult_682;
        double l_sum_980 = l_mult_978 + (l_mult_864 + l_sum_979);
        double l_mult_981 = -0.1944e4 * l_prod_660;
        double l_sum_982 = l_mult_787 + l_mult_677;
        double l_mult_983 = 0.324e3 * l_prod_665;
        double l_sum_984 = l_mult_970 + (l_mult_973 + l_mult_681);
        double l_sum_985 = l_mult_966 + (l_mult_981 + l_mult_676);
        double l_sum_986 = l_mult_821 + (l_mult_822 + (l_mult_823 + l_mult_675));
        double l_mult_987 = -0.3078e4 * l_prod_666;
        double l_mult_988 = 0.12852e5 * l_prod_665;
        double l_mult_989 = -0.17496e5 * l_prod_664;
        double l_mult_990 = 0.7776e4 * l_prod_663;
        double l_mult_991 = 0.12852e5 * l_prod_661;
        double l_mult_992 = 0.23328e5 * l_prod_659;
        double l_mult_993 = -0.17496e5 * l_prod_657;
        double l_mult_994 = 0.23328e5 * l_prod_656;
        double l_mult_995 = 0.7776e4 * l_prod_654;
        double l_mult_996 = 0.12852e5 * l_prod_644;
        double l_mult_997 = -0.34992e5 * l_prod_643;
        double l_mult_998 = -0.34992e5 * l_prod_639;
        double l_mult_999 = 0.46656e5 * l_prod_638;
        double l_mult_1000 = 0.23328e5 * l_prod_622;
        double l_mult_1001 = 0.23328e5 * l_prod_619;
        double l_mult_1002 = 0.1332e4 * l_prod_666;
        double l_mult_1003 = -0.3240e4 * l_prod_665;
        double l_mult_1004 = 0.1944e4 * l_prod_664;
        double l_mult_1005 = -0.11232e5 * l_prod_661;
        double l_mult_1006 = 0.23328e5 * l_prod_660;
        double l_mult_1007 = 0.21384e5 * l_prod_657;
        double l_mult_1008 = -0.23328e5 * l_prod_656;
        double l_mult_1009 = -0.11664e5 * l_prod_654;
        double l_mult_1010 = -0.3240e4 * l_prod_644;
        double l_mult_1011 = 0.23328e5 * l_prod_639;
        double l_mult_1012 = 0.1944e4 * l_prod_623;
        double l_mult_1013 = -0.396e3 * l_prod_666;
        double l_mult_1014 = 0.432e3 * l_prod_665;
        double l_mult_1015 = 0.3996e4 * l_prod_661;
        double l_mult_1016 = -0.3888e4 * l_prod_660;
        double l_mult_1017 = -0.11016e5 * l_prod_657;
        double l_mult_1018 = 0.7776e4 * l_prod_656;
        double l_mult_1019 = 0.54e2 * l_prod_666;
        double l_mult_1020 = -0.594e3 * l_prod_661;
        double l_mult_1021 = 0.1944e4 * l_prod_657;
        double l_mult_1022 = -0.11232e5 * l_prod_665;
        double l_mult_1023 = 0.21384e5 * l_prod_664;
        double l_mult_1024 = -0.11664e5 * l_prod_663;
        double l_mult_1025 = -0.3240e4 * l_prod_661;
        double l_mult_1026 = -0.23328e5 * l_prod_659;
        double l_mult_1027 = 0.23328e5 * l_prod_643;
        double l_mult_1028 = -0.297e3 * l_prod_666;
        double l_mult_1029 = 0.2106e4 * l_prod_665;
        double l_mult_1030 = 0.2106e4 * l_prod_661;
        double l_mult_1031 = 0.36e2 * l_prod_666;
        double l_mult_1032 = -0.324e3 * l_prod_661;
        double l_mult_1033 = 0.648e3 * l_prod_657;
        double l_mult_1034 = 0.3996e4 * l_prod_665;
        double l_mult_1035 = -0.11016e5 * l_prod_664;
        double l_mult_1036 = 0.432e3 * l_prod_661;
        double l_mult_1037 = 0.7776e4 * l_prod_659;
        double l_mult_1038 = -0.324e3 * l_prod_665;
        double l_mult_1039 = 0.648e3 * l_prod_664;
        double l_mult_1040 = -0.594e3 * l_prod_665;
        double l_mult_1041 = 0.540e3 * l_prod_672;
        double l_mult_1042 = -0.3078e4 * l_prod_671;
        double l_mult_1043 = 0.6426e4 * l_prod_670;
        double l_mult_1044 = -0.5832e4 * l_prod_669;
        double l_mult_1045 = -0.17496e5 * l_prod_660;
        double l_mult_1046 = -0.6156e4 * l_prod_651;
        double l_mult_1047 = 0.15552e5 * l_prod_648;
        double l_mult_1048 = 0.25704e5 * l_prod_644;
        double l_mult_1049 = -0.69984e5 * l_prod_643;
        double l_mult_1050 = 0.46656e5 * l_prod_642;
        double l_mult_1051 = 0.15552e5 * l_prod_635;
        double l_mult_1052 = -0.52488e5 * l_prod_628;
        double l_mult_1053 = 0.34992e5 * l_prod_627;
        double l_mult_1054 = -0.52488e5 * l_prod_623;
        double l_mult_1055 = 0.69984e5 * l_prod_622;
        double l_mult_1056 = -0.23328e5 * l_prod_613;
        double l_mult_1057 = 0.31104e5 * l_prod_612;
        double l_mult_1058 = 0.31104e5 * l_prod_608;
        double l_mult_1059 = -0.360e3 * l_prod_672;
        double l_mult_1060 = 0.1332e4 * l_prod_671;
        double l_mult_1061 = -0.1620e4 * l_prod_670;
        double l_mult_1062 = -0.1620e4 * l_prod_661;
        double l_mult_1063 = 0.6984e4 * l_prod_651;
        double l_mult_1064 = -0.22464e5 * l_prod_650;
        double l_mult_1065 = -0.22464e5 * l_prod_644;
        double l_mult_1066 = 0.46656e5 * l_prod_643;
        double l_mult_1067 = 0.64152e5 * l_prod_628;
        double l_mult_1068 = -0.34992e5 * l_prod_627;
        double l_mult_1069 = 0.64152e5 * l_prod_623;
        double l_mult_1070 = -0.69984e5 * l_prod_622;
        double l_mult_1071 = -0.46656e5 * l_prod_608;
        double l_mult_1072 = 0.180e3 * l_prod_672;
        double l_mult_1073 = -0.396e3 * l_prod_671;
        double l_mult_1074 = 0.216e3 * l_prod_670;
        double l_mult_1075 = -0.4032e4 * l_prod_651;
        double l_mult_1076 = 0.7992e4 * l_prod_650;
        double l_mult_1077 = -0.3888e4 * l_prod_649;
        double l_mult_1078 = 0.7992e4 * l_prod_644;
        double l_mult_1079 = -0.7776e4 * l_prod_643;
        double l_mult_1080 = -0.33048e5 * l_prod_628;
        double l_mult_1081 = -0.33048e5 * l_prod_623;
        double l_mult_1082 = -0.54e2 * l_prod_672;
        double l_mult_1083 = 0.54e2 * l_prod_671;
        double l_mult_1084 = 0.1296e4 * l_prod_651;
        double l_mult_1085 = -0.1188e4 * l_prod_650;
        double l_mult_1086 = -0.1188e4 * l_prod_644;
        double l_mult_1087 = 0.5832e4 * l_prod_628;
        double l_mult_1088 = 0.5832e4 * l_prod_623;
        double l_mult_1089 = 0.15552e5 * l_prod_613;
        double l_mult_1090 = 0.3492e4 * l_prod_671;
        double l_mult_1091 = -0.9612e4 * l_prod_670;
        double l_mult_1092 = 0.10368e5 * l_prod_669;
        double l_mult_1093 = 0.11664e5 * l_prod_660;
        double l_mult_1094 = 0.2664e4 * l_prod_651;
        double l_mult_1095 = -0.6480e4 * l_prod_644;
        double l_mult_1096 = -0.46656e5 * l_prod_642;
        double l_mult_1097 = 0.34992e5 * l_prod_628;
        double l_sum_1098 = l_mult_962 + l_mult_916;
        double l_mult_1099 = 0.135e3 * l_prod_672;
        double l_mult_1100 = -0.1107e4 * l_prod_671;
        double l_mult_1101 = 0.1944e4 * l_prod_670;
        double l_mult_1102 = -0.972e3 * l_prod_669;
        double l_mult_1103 = -0.2214e4 * l_prod_651;
        double l_mult_1104 = 0.4212e4 * l_prod_644;
        double l_mult_1105 = -0.29160e5 * l_prod_643;
        double l_mult_1106 = -0.40824e5 * l_prod_628;
        double l_mult_1107 = -0.3888e4 * l_prod_613;
        double l_mult_1108 = -0.36e2 * l_prod_672;
        double l_mult_1109 = 0.252e3 * l_prod_671;
        double l_mult_1110 = -0.216e3 * l_prod_670;
        double l_mult_1111 = 0.720e3 * l_prod_651;
        double l_mult_1112 = -0.4968e4 * l_prod_650;
        double l_mult_1113 = 0.3888e4 * l_prod_649;
        double l_mult_1114 = -0.648e3 * l_prod_644;
        double l_mult_1115 = 0.19440e5 * l_prod_628;
        double l_mult_1116 = -0.2016e4 * l_prod_671;
        double l_mult_1117 = 0.7020e4 * l_prod_670;
        double l_mult_1118 = -0.9072e4 * l_prod_669;
        double l_mult_1119 = -0.792e3 * l_prod_651;
        double l_mult_1120 = 0.864e3 * l_prod_644;
        double l_mult_1121 = 0.15552e5 * l_prod_642;
        double l_mult_1122 = 0.648e3 * l_prod_629;
        double l_mult_1123 = -0.5832e4 * l_prod_628;
        double l_mult_1124 = 0.360e3 * l_prod_671;
        double l_mult_1125 = -0.972e3 * l_prod_670;
        double l_mult_1126 = 0.504e3 * l_prod_651;
        double l_mult_1127 = 0.648e3 * l_prod_671;
        double l_mult_1128 = -0.2538e4 * l_prod_670;
        double l_mult_1129 = 0.3888e4 * l_prod_669;
        double l_mult_1130 = 0.108e3 * l_prod_651;
        double l_mult_1131 = 0.540e3 * l_prod_667;
        double l_mult_1132 = -0.3078e4 * l_prod_662;
        double l_mult_1133 = 0.6426e4 * l_prod_658;
        double l_mult_1134 = -0.5832e4 * l_prod_655;
        double l_mult_1135 = -0.6156e4 * l_prod_645;
        double l_mult_1136 = -0.69984e5 * l_prod_639;
        double l_mult_1137 = 0.46656e5 * l_prod_635;
        double l_mult_1138 = 0.15552e5 * l_prod_633;
        double l_mult_1139 = -0.52488e5 * l_prod_620;
        double l_mult_1140 = 0.69984e5 * l_prod_619;
        double l_mult_1141 = 0.34992e5 * l_prod_617;
        double l_mult_1142 = -0.23328e5 * l_prod_609;
        double l_mult_1143 = 0.31104e5 * l_prod_606;
        double l_mult_1144 = -0.360e3 * l_prod_667;
        double l_mult_1145 = -0.1620e4 * l_prod_665;
        double l_mult_1146 = 0.1332e4 * l_prod_662;
        double l_mult_1147 = -0.1620e4 * l_prod_658;
        double l_mult_1148 = 0.6984e4 * l_prod_645;
        double l_mult_1149 = -0.22464e5 * l_prod_640;
        double l_mult_1150 = 0.46656e5 * l_prod_639;
        double l_mult_1151 = 0.64152e5 * l_prod_620;
        double l_mult_1152 = -0.69984e5 * l_prod_619;
        double l_mult_1153 = -0.34992e5 * l_prod_617;
        double l_mult_1154 = 0.180e3 * l_prod_667;
        double l_mult_1155 = -0.396e3 * l_prod_662;
        double l_mult_1156 = 0.216e3 * l_prod_658;
        double l_mult_1157 = -0.4032e4 * l_prod_645;
        double l_mult_1158 = 0.7992e4 * l_prod_640;
        double l_mult_1159 = -0.7776e4 * l_prod_639;
        double l_mult_1160 = -0.3888e4 * l_prod_636;
        double l_mult_1161 = -0.33048e5 * l_prod_620;
        double l_mult_1162 = -0.54e2 * l_prod_667;
        double l_mult_1163 = 0.54e2 * l_prod_662;
        double l_mult_1164 = 0.1296e4 * l_prod_645;
        double l_mult_1165 = -0.1188e4 * l_prod_640;
        double l_mult_1166 = 0.5832e4 * l_prod_620;
        double l_mult_1167 = 0.15552e5 * l_prod_609;
        double l_mult_1168 = 0.3492e4 * l_prod_662;
        double l_mult_1169 = -0.9612e4 * l_prod_658;
        double l_mult_1170 = 0.10368e5 * l_prod_655;
        double l_mult_1171 = 0.2664e4 * l_prod_645;
        double l_mult_1172 = -0.46656e5 * l_prod_635;
        double l_mult_1173 = 0.34992e5 * l_prod_620;
        double l_sum_1174 = l_mult_963 + l_mult_918;
        double l_mult_1175 = 0.135e3 * l_prod_667;
        double l_mult_1176 = -0.1107e4 * l_prod_662;
        double l_mult_1177 = 0.1944e4 * l_prod_658;
        double l_mult_1178 = -0.972e3 * l_prod_655;
        double l_mult_1179 = -0.2214e4 * l_prod_645;
        double l_mult_1180 = -0.29160e5 * l_prod_639;
        double l_mult_1181 = -0.40824e5 * l_prod_620;
        double l_mult_1182 = -0.3888e4 * l_prod_609;
        double l_mult_1183 = -0.36e2 * l_prod_667;
        double l_mult_1184 = 0.252e3 * l_prod_662;
        double l_mult_1185 = -0.216e3 * l_prod_658;
        double l_mult_1186 = 0.720e3 * l_prod_645;
        double l_mult_1187 = -0.4968e4 * l_prod_640;
        double l_mult_1188 = 0.3888e4 * l_prod_636;
        double l_mult_1189 = 0.19440e5 * l_prod_620;
        double l_mult_1190 = -0.2016e4 * l_prod_662;
        double l_mult_1191 = 0.7020e4 * l_prod_658;
        double l_mult_1192 = -0.9072e4 * l_prod_655;
        double l_mult_1193 = -0.792e3 * l_prod_645;
        double l_mult_1194 = 0.648e3 * l_prod_624;
        double l_mult_1195 = -0.5832e4 * l_prod_620;
        double l_mult_1196 = 0.360e3 * l_prod_662;
        double l_mult_1197 = -0.972e3 * l_prod_658;
        double l_mult_1198 = 0.504e3 * l_prod_645;
        double l_mult_1199 = 0.648e3 * l_prod_662;
        double l_mult_1200 = -0.2538e4 * l_prod_658;
        double l_mult_1201 = 0.3888e4 * l_prod_655;
        double l_mult_1202 = 0.108e3 * l_prod_645;
        double l_mult_1203 = -0.31968e5 * l_prod_644;
        double l_mult_1204 = 0.77760e5 * l_prod_643;
        double l_mult_1205 = 0.77760e5 * l_prod_639;
        double l_mult_1206 = -0.1620e4 * l_prod_666;
        double l_mult_1207 = 0.3564e4 * l_prod_665;
        double l_mult_1208 = 0.3564e4 * l_prod_661;
        double l_mult_1209 = 0.26568e5 * l_prod_644;
        double l_mult_1210 = -0.50544e5 * l_prod_643;
        double l_mult_1211 = -0.50544e5 * l_prod_639;
        double l_mult_1212 = -0.69984e5 * l_prod_623;
        double l_mult_1213 = 0.432e3 * l_prod_666;
        double l_mult_1214 = -0.432e3 * l_prod_665;
        double l_mult_1215 = -0.432e3 * l_prod_661;
        double l_mult_1216 = -0.8640e4 * l_prod_644;
        double l_mult_1217 = 0.7776e4 * l_prod_643;
        double l_mult_1218 = 0.7776e4 * l_prod_639;
        double l_mult_1219 = -0.25272e5 * l_prod_660;
        double l_mult_1220 = -0.23328e5 * l_prod_657;
        double l_mult_1221 = 0.7128e4 * l_prod_644;
        double l_mult_1222 = 0.324e3 * l_prod_666;
        double l_mult_1223 = -0.2268e4 * l_prod_661;
        double l_mult_1224 = -0.4536e4 * l_prod_644;
        double l_mult_1225 = 0.3888e4 * l_prod_660;
        double l_mult_1226 = -0.864e3 * l_prod_644;
        double l_mult_1227 = -0.23328e5 * l_prod_664;
        double l_mult_1228 = -0.2268e4 * l_prod_665;
        double l_mult_1229 = -0.180e3 * l_prod_666;
        double l_mult_1230 = -0.99e2 * l_prod_666;
        double l_mult_1231 = -0.2916e4 * l_prod_660;
        double l_mult_1232 = -0.72e2 * l_prod_666;
        double l_mult_1233 = 0.648e3 * l_prod_665;
        double l_mult_1234 = 0.72e1 * l_prod_652;
        double l_sum_1235 = l_mult_1234 + (-0.90e2 * l_prod_630 + (0.378e3 * l_prod_614 + (-0.648e3 * l_prod_603 + l_mult_694)));
        double l_mult_1236 = 0.45e1 * l_prod_652;
        double l_mult_1237 = -0.495e2 * l_prod_630;
        double l_mult_1238 = 0.594e3 * l_prod_624;
        double l_mult_1239 = 0.162e3 * l_prod_614;
        double l_mult_1240 = -0.1944e4 * l_prod_609;
        double l_mult_1241 = -0.162e3 * l_prod_603;
        double l_sum_1242 = l_mult_1241 + l_mult_693;
        double l_mult_1243 = 0.4e1 * l_prod_652;
        double l_mult_1244 = 0.216e3 * l_prod_640;
        double l_mult_1245 = -0.36e2 * l_prod_630;
        double l_mult_1246 = 0.72e2 * l_prod_614;
        double l_mult_1247 = -0.648e3 * l_prod_636;
        double l_mult_1248 = -0.27e2 * l_prod_630;
        double l_sum_1249 = l_mult_1248 + (l_mult_1238 + (l_mult_727 + l_mult_688));
        double l_mult_1250 = 0.4176e4 * l_prod_666;
        double l_mult_1251 = -0.3888e4 * l_prod_663;
        double l_mult_1252 = -0.15552e5 * l_prod_656;
        double l_mult_1253 = -0.9720e4 * l_prod_654;
        double l_mult_1254 = -0.23328e4 * l_prod_653;
        double l_mult_1255 = -0.3132e3 * l_prod_652;
        double l_mult_1256 = -0.5022e4 * l_prod_650;
        double l_mult_1257 = 0.5184e4 * l_prod_649;
        double l_mult_1258 = -0.1944e4 * l_prod_648;
        double l_mult_1259 = -0.34992e5 * l_prod_638;
        double l_mult_1260 = -0.31104e5 * l_prod_635;
        double l_mult_1261 = -0.9720e4 * l_prod_633;
        double l_mult_1262 = 0.1044e4 * l_prod_630;
        double l_mult_1263 = -0.10044e5 * l_prod_624;
        double l_mult_1264 = -0.15552e5 * l_prod_617;
        double l_mult_1265 = -0.1674e4 * l_prod_614;
        double l_mult_1266 = -0.3888e4 * l_prod_612;
        double l_mult_1267 = 0.10368e5 * l_prod_609;
        double l_mult_1268 = -0.15552e5 * l_prod_608;
        double l_mult_1269 = -0.11664e5 * l_prod_606;
        double l_mult_1270 = 0.1296e4 * l_prod_603;
        double l_mult_1271 = -0.3888e4 * l_prod_600;
        double l_mult_1272 = -0.3888e3 * l_prod_598;
        double l_mult_1273 = -0.5220e4 * l_prod_666;
        double l_mult_1274 = -0.29160e5 * l_prod_660;
        double l_mult_1275 = 0.19440e5 * l_prod_654;
        double l_mult_1276 = 0.5832e4 * l_prod_653;
        double l_mult_1277 = 0.2565e3 * l_prod_652;
        double l_mult_1278 = 0.1458e4 * l_prod_650;
        double l_mult_1279 = -0.648e3 * l_prod_649;
        double l_mult_1280 = 0.34992e5 * l_prod_638;
        double l_mult_1281 = 0.19440e5 * l_prod_633;
        double l_mult_1282 = -0.5355e3 * l_prod_630;
        double l_mult_1283 = -0.972e3 * l_prod_628;
        double l_mult_1284 = 0.9342e4 * l_prod_624;
        double l_mult_1285 = -0.21384e5 * l_prod_623;
        double l_mult_1286 = 0.23328e5 * l_prod_617;
        double l_mult_1287 = 0.486e3 * l_prod_614;
        double l_mult_1288 = 0.11664e5 * l_prod_606;
        double l_mult_1289 = 0.3384e4 * l_prod_666;
        double l_mult_1290 = -0.19440e5 * l_prod_654;
        double l_mult_1291 = -0.148e3 * l_prod_652;
        double l_mult_1292 = -0.216e3 * l_prod_650;
        double l_mult_1293 = -0.11664e5 * l_prod_638;
        double l_mult_1294 = -0.19440e5 * l_prod_633;
        double l_mult_1295 = 0.180e3 * l_prod_630;
        double l_mult_1296 = -0.3672e4 * l_prod_624;
        double l_mult_1297 = 0.3888e4 * l_prod_623;
        double l_mult_1298 = -0.72e2 * l_prod_614;
        double l_mult_1299 = -0.3888e4 * l_prod_606;
        double l_mult_1300 = -0.1197e4 * l_prod_666;
        double l_mult_1301 = 0.9720e4 * l_prod_654;
        double l_mult_1302 = 0.495e2 * l_prod_652;
        double l_mult_1303 = -0.5832e4 * l_prod_639;
        double l_mult_1304 = 0.9720e4 * l_prod_633;
        double l_mult_1305 = 0.180e3 * l_prod_666;
        double l_mult_1306 = 0.2592e4 * l_prod_657;
        double l_mult_1307 = -0.72e1 * l_prod_652;
        double l_mult_1308 = 0.2592e4 * l_prod_636;
        double l_mult_1309 = -0.1944e4 * l_prod_633;
        double l_mult_1310 = -0.5022e4 * l_prod_640;
        double l_mult_1311 = 0.5184e4 * l_prod_636;
        double l_mult_1312 = -0.10044e5 * l_prod_629;
        double l_mult_1313 = -0.7776e4 * l_prod_627;
        double l_mult_1314 = -0.7776e4 * l_prod_617;
        double l_mult_1315 = -0.11664e5 * l_prod_612;
        double l_mult_1316 = -0.23328e5 * l_prod_608;
        double l_mult_1317 = -0.7776e4 * l_prod_602;
        double l_mult_1318 = -0.7776e4 * l_prod_600;
        double l_mult_1319 = -0.1944e4 * l_prod_598;
        double l_sum_1320 = l_mult_1255 + (l_mult_758 + (l_mult_1256 + (l_mult_1257 + (l_mult_1258 + (l_mult_840 + (l_mult_762 + (l_mult_841 + (l_mult_813 + (l_mult_1310 + (l_mult_765 + (l_mult_1293 + (l_mult_1311 + (l_mult_767 + (l_mult_1309 + (0.2088e4 * l_prod_630 + (l_mult_1312 + (l_mult_769 + (l_mult_1313 + (l_mult_1263 + (l_mult_949 + (l_mult_772 + (l_mult_849 + (l_mult_850 + (l_mult_1314 + (-0.5022e4 * l_prod_614 + (l_mult_1089 + (l_mult_1315 + (l_mult_1167 + (l_mult_1316 + (l_mult_1269 + (l_mult_964 + (l_mult_1317 + (l_mult_1318 + l_mult_1319)))))))))))))))))))))))))))))))));
        double l_mult_1321 = 0.1458e4 * l_prod_640;
        double l_mult_1322 = 0.9342e4 * l_prod_629;
        double l_mult_1323 = 0.11664e5 * l_prod_612;
        double l_mult_1324 = 0.23328e5 * l_prod_608;
        double l_mult_1325 = 0.11664e5 * l_prod_602;
        double l_mult_1326 = 0.11664e5 * l_prod_600;
        double l_mult_1327 = 0.3888e4 * l_prod_598;
        double l_sum_1328 = l_mult_1277 + (l_mult_788 + (l_mult_1278 + (l_mult_1279 + (l_mult_863 + (l_mult_791 + (l_mult_864 + (l_mult_1321 + (l_mult_794 + (l_mult_1247 + (-0.2610e4 * l_prod_630 + (l_mult_1322 + (l_mult_796 + (l_mult_685 + (l_mult_1284 + (l_mult_1285 + (l_mult_686 + (l_mult_870 + (l_mult_687 + (l_mult_688 + (0.7884e4 * l_prod_614 + (-0.19440e5 * l_prod_613 + (l_mult_1323 + (-0.19440e5 * l_prod_609 + (l_mult_1324 + (l_mult_1288 + (-0.9396e4 * l_prod_603 + (l_mult_1325 + (l_mult_1326 + l_mult_1327))))))))))))))))))))))))))));
        double l_mult_1329 = -0.216e3 * l_prod_640;
        double l_mult_1330 = -0.3672e4 * l_prod_629;
        double l_mult_1331 = 0.10368e5 * l_prod_613;
        double l_mult_1332 = -0.3888e4 * l_prod_598;
        double l_sum_1333 = l_mult_1291 + (l_mult_808 + (l_mult_1292 + (l_mult_879 + (l_mult_811 + (l_mult_1329 + (0.1692e4 * l_prod_630 + (l_mult_1330 + (l_mult_815 + (l_mult_1296 + (l_mult_1297 + (l_mult_884 + (-0.6120e4 * l_prod_614 + (l_mult_1331 + (l_mult_1266 + (l_mult_1267 + (l_mult_776 + (l_mult_1299 + (0.8424e4 * l_prod_603 + (l_mult_1317 + (l_mult_1318 + l_mult_1332))))))))))))))))))));
        double l_mult_1334 = 0.594e3 * l_prod_629;
        double l_mult_1335 = -0.1944e4 * l_prod_613;
        double l_mult_1336 = 0.1944e4 * l_prod_598;
        double l_sum_1337 = l_mult_1302 + (l_mult_717 + (l_mult_741 + (-0.5985e3 * l_prod_630 + (l_mult_1334 + (l_mult_1238 + (0.2376e4 * l_prod_614 + (l_mult_1335 + (l_mult_1240 + (-0.3726e4 * l_prod_603 + (l_mult_692 + (l_mult_693 + l_mult_1336)))))))))));
        double l_mult_1338 = 0.648e3 * l_prod_603;
        double l_sum_1339 = l_mult_1307 + (0.90e2 * l_prod_630 + (-0.378e3 * l_prod_614 + (l_mult_1338 + l_mult_1272)));
        double l_mult_1340 = -0.36e2 * l_prod_651;
        double l_mult_1341 = 0.324e3 * l_prod_629;
        double l_mult_1342 = -0.3888e4 * l_prod_623;
        double l_sum_1343 = l_mult_711 + (l_mult_1342 + l_mult_687);
        double l_mult_1344 = 0.216e3 * l_prod_650;
        double l_mult_1345 = -0.27e2 * l_prod_651;
        double l_mult_1346 = 0.162e3 * l_prod_650;
        double l_sum_1347 = 0.162e3 * l_prod_629 + (l_mult_1283 + l_sum_974);
        double l_mult_1348 = 0.324e3 * l_prod_650;
        double l_mult_1349 = -0.6156e4 * l_prod_666;
        double l_mult_1350 = 0.15552e5 * l_prod_663;
        double l_mult_1351 = -0.52488e5 * l_prod_660;
        double l_mult_1352 = 0.34992e5 * l_prod_659;
        double l_mult_1353 = 0.31104e5 * l_prod_656;
        double l_mult_1354 = -0.3078e4 * l_prod_651;
        double l_mult_1355 = 0.12852e5 * l_prod_650;
        double l_mult_1356 = -0.17496e5 * l_prod_649;
        double l_mult_1357 = 0.7776e4 * l_prod_648;
        double l_mult_1358 = -0.52488e5 * l_prod_639;
        double l_mult_1359 = 0.69984e5 * l_prod_638;
        double l_mult_1360 = 0.31104e5 * l_prod_635;
        double l_mult_1361 = -0.17496e5 * l_prod_628;
        double l_mult_1362 = -0.34992e5 * l_prod_623;
        double l_mult_1363 = 0.46656e5 * l_prod_622;
        double l_mult_1364 = 0.7776e4 * l_prod_612;
        double l_mult_1365 = 0.15552e5 * l_prod_608;
        double l_mult_1366 = 0.6984e4 * l_prod_666;
        double l_mult_1367 = -0.22464e5 * l_prod_665;
        double l_mult_1368 = 0.64152e5 * l_prod_660;
        double l_mult_1369 = -0.34992e5 * l_prod_659;
        double l_mult_1370 = 0.1332e4 * l_prod_651;
        double l_mult_1371 = -0.3240e4 * l_prod_650;
        double l_mult_1372 = 0.1944e4 * l_prod_649;
        double l_mult_1373 = 0.64152e5 * l_prod_639;
        double l_mult_1374 = -0.69984e5 * l_prod_638;
        double l_mult_1375 = -0.1620e4 * l_prod_629;
        double l_mult_1376 = 0.23328e5 * l_prod_623;
        double l_mult_1377 = 0.648e3 * l_prod_613;
        double l_mult_1378 = -0.4032e4 * l_prod_666;
        double l_mult_1379 = 0.7992e4 * l_prod_665;
        double l_mult_1380 = -0.3888e4 * l_prod_664;
        double l_mult_1381 = -0.33048e5 * l_prod_660;
        double l_mult_1382 = -0.396e3 * l_prod_651;
        double l_mult_1383 = 0.432e3 * l_prod_650;
        double l_mult_1384 = -0.33048e5 * l_prod_639;
        double l_mult_1385 = 0.23328e5 * l_prod_638;
        double l_mult_1386 = 0.1296e4 * l_prod_666;
        double l_mult_1387 = -0.1188e4 * l_prod_665;
        double l_mult_1388 = 0.5832e4 * l_prod_660;
        double l_mult_1389 = 0.54e2 * l_prod_651;
        double l_mult_1390 = 0.5832e4 * l_prod_639;
        double l_mult_1391 = 0.2664e4 * l_prod_666;
        double l_mult_1392 = 0.34992e5 * l_prod_660;
        double l_mult_1393 = -0.11232e5 * l_prod_650;
        double l_mult_1394 = 0.21384e5 * l_prod_649;
        double l_mult_1395 = -0.11664e5 * l_prod_648;
        double l_mult_1396 = 0.11664e5 * l_prod_628;
        double l_sum_1397 = l_mult_1377 + l_mult_1266;
        double l_mult_1398 = -0.2214e4 * l_prod_666;
        double l_mult_1399 = -0.40824e5 * l_prod_660;
        double l_mult_1400 = -0.3888e4 * l_prod_657;
        double l_mult_1401 = -0.297e3 * l_prod_651;
        double l_mult_1402 = 0.2106e4 * l_prod_650;
        double l_mult_1403 = 0.720e3 * l_prod_666;
        double l_mult_1404 = -0.4968e4 * l_prod_665;
        double l_mult_1405 = 0.3888e4 * l_prod_664;
        double l_mult_1406 = 0.19440e5 * l_prod_660;
        double l_mult_1407 = 0.36e2 * l_prod_651;
        double l_mult_1408 = 0.1944e4 * l_prod_639;
        double l_mult_1409 = -0.792e3 * l_prod_666;
        double l_mult_1410 = 0.648e3 * l_prod_661;
        double l_mult_1411 = -0.5832e4 * l_prod_660;
        double l_mult_1412 = 0.3996e4 * l_prod_650;
        double l_mult_1413 = -0.11016e5 * l_prod_649;
        double l_mult_1414 = 0.504e3 * l_prod_666;
        double l_mult_1415 = -0.324e3 * l_prod_650;
        double l_mult_1416 = 0.648e3 * l_prod_649;
        double l_mult_1417 = 0.108e3 * l_prod_666;
        double l_sum_1418 = l_mult_1389 + (-0.594e3 * l_prod_650 + (l_mult_1372 + l_mult_1258));
        double l_mult_1419 = 0.12852e5 * l_prod_629;
        double l_mult_1420 = 0.23328e5 * l_prod_627;
        double l_mult_1421 = -0.17496e5 * l_prod_613;
        double l_mult_1422 = 0.7776e4 * l_prod_602;
        double l_mult_1423 = -0.11232e5 * l_prod_629;
        double l_mult_1424 = 0.21384e5 * l_prod_613;
        double l_mult_1425 = -0.23328e5 * l_prod_612;
        double l_mult_1426 = -0.11664e5 * l_prod_602;
        double l_mult_1427 = 0.3996e4 * l_prod_629;
        double l_mult_1428 = -0.3888e4 * l_prod_628;
        double l_mult_1429 = -0.11016e5 * l_prod_613;
        double l_mult_1430 = -0.594e3 * l_prod_629;
        double l_mult_1431 = 0.1944e4 * l_prod_613;
        double l_mult_1432 = -0.3240e4 * l_prod_629;
        double l_mult_1433 = -0.23328e5 * l_prod_627;
        double l_sum_1434 = l_mult_1297 + (l_mult_772 + (l_mult_1431 + l_mult_1315));
        double l_mult_1435 = 0.2106e4 * l_prod_629;
        double l_sum_1436 = l_mult_1335 + l_mult_1323;
        double l_mult_1437 = -0.324e3 * l_prod_629;
        double l_mult_1438 = 0.432e3 * l_prod_629;
        double l_mult_1439 = 0.7776e4 * l_prod_627;
        double l_mult_1440 = 0.540e3 * l_prod_652;
        double l_mult_1441 = -0.23328e5 * l_prod_636;
        double l_mult_1442 = -0.3078e4 * l_prod_630;
        double l_mult_1443 = 0.31104e5 * l_prod_617;
        double l_mult_1444 = 0.6426e4 * l_prod_614;
        double l_mult_1445 = 0.34992e5 * l_prod_606;
        double l_mult_1446 = -0.5832e4 * l_prod_603;
        double l_mult_1447 = 0.15552e5 * l_prod_600;
        double l_mult_1448 = -0.360e3 * l_prod_652;
        double l_mult_1449 = -0.1620e4 * l_prod_650;
        double l_mult_1450 = 0.3492e4 * l_prod_630;
        double l_mult_1451 = -0.22464e5 * l_prod_624;
        double l_mult_1452 = -0.9612e4 * l_prod_614;
        double l_mult_1453 = -0.34992e5 * l_prod_606;
        double l_mult_1454 = 0.10368e5 * l_prod_603;
        double l_mult_1455 = 0.180e3 * l_prod_652;
        double l_mult_1456 = -0.2016e4 * l_prod_630;
        double l_mult_1457 = 0.7992e4 * l_prod_624;
        double l_mult_1458 = -0.7776e4 * l_prod_623;
        double l_mult_1459 = 0.7020e4 * l_prod_614;
        double l_mult_1460 = -0.9072e4 * l_prod_603;
        double l_mult_1461 = -0.54e2 * l_prod_652;
        double l_mult_1462 = 0.648e3 * l_prod_630;
        double l_mult_1463 = -0.1188e4 * l_prod_624;
        double l_mult_1464 = -0.2538e4 * l_prod_614;
        double l_mult_1465 = 0.3888e4 * l_prod_609;
        double l_mult_1466 = 0.3888e4 * l_prod_603;
        double l_mult_1467 = 0.1332e4 * l_prod_630;
        double l_mult_1468 = -0.1620e4 * l_prod_614;
        double l_sum_1469 = l_mult_1338 + l_mult_1318;
        double l_mult_1470 = 0.135e3 * l_prod_652;
        double l_mult_1471 = -0.1107e4 * l_prod_630;
        double l_mult_1472 = -0.29160e5 * l_prod_623;
        double l_mult_1473 = 0.1944e4 * l_prod_614;
        double l_mult_1474 = -0.972e3 * l_prod_603;
        double l_mult_1475 = -0.36e2 * l_prod_652;
        double l_mult_1476 = 0.360e3 * l_prod_630;
        double l_mult_1477 = -0.4968e4 * l_prod_624;
        double l_mult_1478 = -0.972e3 * l_prod_614;
        double l_mult_1479 = -0.396e3 * l_prod_630;
        double l_mult_1480 = 0.216e3 * l_prod_614;
        double l_mult_1481 = 0.252e3 * l_prod_630;
        double l_mult_1482 = -0.216e3 * l_prod_614;
        double l_mult_1483 = 0.54e2 * l_prod_630;
        double l_mult_1484 = 0.77760e5 * l_prod_623;
        double l_mult_1485 = -0.1620e4 * l_prod_651;
        double l_mult_1486 = 0.3564e4 * l_prod_650;
        double l_mult_1487 = -0.25272e5 * l_prod_628;
        double l_mult_1488 = -0.50544e5 * l_prod_623;
        double l_mult_1489 = 0.432e3 * l_prod_651;
        double l_mult_1490 = -0.432e3 * l_prod_650;
        double l_mult_1491 = 0.3888e4 * l_prod_628;
        double l_mult_1492 = 0.7776e4 * l_prod_623;
        double l_mult_1493 = 0.3564e4 * l_prod_629;
        double l_mult_1494 = 0.324e3 * l_prod_651;
        double l_mult_1495 = -0.2268e4 * l_prod_629;
        double l_mult_1496 = -0.432e3 * l_prod_629;
        double l_mult_1497 = -0.23328e5 * l_prod_649;
        double l_mult_1498 = -0.2268e4 * l_prod_650;
        double l_sum_1499 = l_mult_1241 + l_mult_692;
        double l_sum_1500 = l_mult_1248 + (l_mult_1334 + (l_mult_703 + l_mult_685));
        double l_mult_1501 = -0.23328e4 * l_prod_668;
        double l_mult_1502 = -0.9720e4 * l_prod_663;
        double l_mult_1503 = -0.15552e5 * l_prod_659;
        double l_mult_1504 = -0.3888e4 * l_prod_654;
        double l_mult_1505 = -0.9720e4 * l_prod_648;
        double l_mult_1506 = -0.31104e5 * l_prod_642;
        double l_mult_1507 = -0.15552e5 * l_prod_627;
        double l_mult_1508 = -0.3888e4 * l_prod_602;
        double l_mult_1509 = 0.5832e4 * l_prod_668;
        double l_mult_1510 = 0.19440e5 * l_prod_663;
        double l_mult_1511 = 0.19440e5 * l_prod_648;
        double l_mult_1512 = -0.972e3 * l_prod_620;
        double l_mult_1513 = -0.19440e5 * l_prod_663;
        double l_mult_1514 = -0.19440e5 * l_prod_648;
        double l_mult_1515 = 0.9720e4 * l_prod_663;
        double l_mult_1516 = 0.9720e4 * l_prod_648;
        double l_mult_1517 = -0.5832e4 * l_prod_643;
        double l_mult_1518 = 0.2592e4 * l_prod_664;
        double l_mult_1519 = 0.2592e4 * l_prod_649;
        double l_mult_1520 = -0.36e2 * l_prod_645;
        double l_mult_1521 = 0.324e3 * l_prod_624;
        double l_mult_1522 = 0.324e3 * l_prod_640;
        double l_mult_1523 = -0.27e2 * l_prod_645;
        double l_mult_1524 = 0.162e3 * l_prod_640;
        double l_sum_1525 = 0.162e3 * l_prod_624 + (l_mult_798 + (l_mult_1512 + l_mult_687));
        double l_sum_1526 = l_mult_735 + (l_mult_1342 + l_mult_686);
        double l_mult_1527 = 0.31104e5 * l_prod_659;
        double l_mult_1528 = 0.34992e5 * l_prod_656;
        double l_mult_1529 = 0.15552e5 * l_prod_654;
        double l_mult_1530 = -0.3078e4 * l_prod_645;
        double l_mult_1531 = -0.52488e5 * l_prod_643;
        double l_mult_1532 = 0.31104e5 * l_prod_642;
        double l_mult_1533 = 0.12852e5 * l_prod_640;
        double l_mult_1534 = -0.17496e5 * l_prod_636;
        double l_mult_1535 = 0.7776e4 * l_prod_633;
        double l_mult_1536 = -0.17496e5 * l_prod_620;
        double l_mult_1537 = 0.46656e5 * l_prod_619;
        double l_mult_1538 = 0.7776e4 * l_prod_606;
        double l_mult_1539 = -0.22464e5 * l_prod_661;
        double l_mult_1540 = -0.34992e5 * l_prod_656;
        double l_mult_1541 = 0.1332e4 * l_prod_645;
        double l_mult_1542 = 0.5832e4 * l_prod_643;
        double l_mult_1543 = -0.11232e5 * l_prod_640;
        double l_mult_1544 = 0.21384e5 * l_prod_636;
        double l_mult_1545 = -0.11664e5 * l_prod_633;
        double l_mult_1546 = -0.1620e4 * l_prod_624;
        double l_mult_1547 = 0.11664e5 * l_prod_620;
        double l_mult_1548 = 0.648e3 * l_prod_609;
        double l_sum_1549 = l_mult_1548 + l_mult_1299;
        double l_mult_1550 = 0.7992e4 * l_prod_661;
        double l_mult_1551 = -0.396e3 * l_prod_645;
        double l_mult_1552 = 0.3996e4 * l_prod_640;
        double l_mult_1553 = -0.11016e5 * l_prod_636;
        double l_mult_1554 = -0.1188e4 * l_prod_661;
        double l_mult_1555 = 0.3888e4 * l_prod_657;
        double l_mult_1556 = 0.54e2 * l_prod_645;
        double l_mult_1557 = 0.1944e4 * l_prod_636;
        double l_sum_1558 = l_mult_1556 + (-0.594e3 * l_prod_640 + (l_mult_1557 + l_mult_1309));
        double l_mult_1559 = 0.64152e5 * l_prod_643;
        double l_mult_1560 = -0.3240e4 * l_prod_640;
        double l_mult_1561 = -0.297e3 * l_prod_645;
        double l_mult_1562 = 0.2106e4 * l_prod_640;
        double l_mult_1563 = -0.4968e4 * l_prod_661;
        double l_mult_1564 = 0.36e2 * l_prod_645;
        double l_mult_1565 = -0.324e3 * l_prod_640;
        double l_mult_1566 = 0.648e3 * l_prod_636;
        double l_mult_1567 = -0.33048e5 * l_prod_643;
        double l_mult_1568 = 0.432e3 * l_prod_640;
        double l_mult_1569 = 0.1944e4 * l_prod_643;
        double l_mult_1570 = 0.31104e5 * l_prod_627;
        double l_mult_1571 = 0.12852e5 * l_prod_624;
        double l_mult_1572 = 0.7776e4 * l_prod_617;
        double l_mult_1573 = 0.34992e5 * l_prod_612;
        double l_mult_1574 = -0.17496e5 * l_prod_609;
        double l_mult_1575 = 0.15552e5 * l_prod_602;
        double l_mult_1576 = 0.7776e4 * l_prod_600;
        double l_mult_1577 = -0.1620e4 * l_prod_640;
        double l_mult_1578 = -0.22464e5 * l_prod_629;
        double l_mult_1579 = -0.11232e5 * l_prod_624;
        double l_mult_1580 = -0.34992e5 * l_prod_612;
        double l_mult_1581 = 0.21384e5 * l_prod_609;
        double l_mult_1582 = -0.11664e5 * l_prod_600;
        double l_mult_1583 = 0.7992e4 * l_prod_629;
        double l_mult_1584 = 0.3996e4 * l_prod_624;
        double l_mult_1585 = -0.11016e5 * l_prod_609;
        double l_mult_1586 = -0.1188e4 * l_prod_629;
        double l_mult_1587 = -0.594e3 * l_prod_624;
        double l_mult_1588 = 0.3888e4 * l_prod_613;
        double l_mult_1589 = 0.1944e4 * l_prod_609;
        double l_mult_1590 = -0.3240e4 * l_prod_624;
        double l_sum_1591 = l_mult_1338 + l_mult_1317;
        double l_mult_1592 = 0.2106e4 * l_prod_624;
        double l_mult_1593 = -0.4968e4 * l_prod_629;
        double l_mult_1594 = -0.324e3 * l_prod_624;
        double l_mult_1595 = 0.432e3 * l_prod_624;
        double l_mult_1596 = -0.23328e5 * l_prod_606;
        double l_mult_1597 = -0.3888e4 * l_prod_620;
        double l_mult_1598 = -0.23328e5 * l_prod_617;
        double l_sum_1599 = l_mult_1589 + l_mult_1269;
        double l_sum_1600 = l_mult_1240 + l_mult_1288;
        double l_mult_1601 = -0.1620e4 * l_prod_645;
        double l_mult_1602 = 0.3564e4 * l_prod_640;
        double l_mult_1603 = -0.25272e5 * l_prod_620;
        double l_mult_1604 = 0.432e3 * l_prod_645;
        double l_mult_1605 = -0.432e3 * l_prod_640;
        double l_mult_1606 = 0.3888e4 * l_prod_620;
        double l_mult_1607 = 0.3564e4 * l_prod_624;
        double l_mult_1608 = 0.324e3 * l_prod_645;
        double l_mult_1609 = -0.2268e4 * l_prod_640;
        double l_mult_1610 = -0.2268e4 * l_prod_624;
        double l_mult_1611 = -0.432e3 * l_prod_624;
        double l_vdot_1612 = vdot4(v_571, vcons4(0.e0, 0.e0, 0.e0, 0.e0)) + (vdot4(v_572,
            vcons4(0.e0, l_mult_697 + (l_mult_698 + (0.1134e4 * l_prod_629 + (-0.2592e4 * l_prod_613 + l_mult_692))),
                l_mult_699 + (l_mult_700 + (l_mult_701 + (l_mult_702 + (0.486e3 * l_prod_629 + (l_mult_703 + l_sum_705))))),
                l_mult_706 + (l_mult_707 + (l_mult_708 + (l_mult_709 + (l_mult_710 + (-0.1296e4 * l_prod_649 + l_sum_713))))))) + (vdot4(
            v_573,
            vcons4(l_mult_699 + (l_mult_714 + (l_mult_715 + (l_mult_716 + l_sum_719))), l_sum_720,
                l_mult_721 + (l_mult_722 + (0.1134e4 * l_prod_624 + (-0.2592e4 * l_prod_609 + l_mult_693))),
                l_mult_723 + (l_mult_724 + (l_mult_725 + (l_mult_726 + (0.486e3 * l_prod_624 + (l_mult_727 + l_sum_729))))))) + (vdot4(
            v_574,
            vcons4(
                l_mult_730 + (l_mult_731 + (l_mult_732 + (l_mult_733 + (l_mult_734 + (-0.1296e4 * l_prod_636 + l_sum_737))))),
                l_mult_723 + (l_mult_738 + (l_mult_739 + (l_mult_740 + l_sum_743))), l_sum_744, l_sum_778)) + (vdot4(
            v_575, vcons4(l_sum_799, l_sum_818, l_sum_824, l_sum_828)) + (vdot4(v_576,
            vcons4(l_sum_855, l_sum_872, l_sum_887, l_sum_892)) + (vdot4(v_577,
            vcons4(l_sum_896,
                l_mult_897 + (l_mult_745 + (l_mult_898 + (l_mult_899 + (l_mult_900 + (l_mult_827 + (l_mult_829 + (l_mult_748 + (l_mult_830 + (l_mult_831 + (l_mult_832 + (l_mult_901 + (l_mult_752 + (0.7776e4 * l_prod_660 + (l_mult_807 + (l_mult_902 + (l_mult_755 + (l_mult_877 + (l_mult_903 + (l_mult_757 + (l_mult_895 + (-0.6264e3 * l_prod_652 + (l_mult_904 + (l_mult_759 + (l_mult_810 + (l_mult_905 + (l_mult_906 + (l_mult_907 + (l_mult_763 + (l_mult_908 + (l_mult_842 + (l_mult_843 + (l_mult_766 + (l_mult_882 + (l_mult_909 + (l_mult_910 + (0.3132e4 * l_prod_630 + (-0.15066e5 * l_prod_629 + (l_mult_911 + (l_mult_770 + (-0.15066e5 * l_prod_624 + (l_mult_912 + (l_mult_913 + (l_mult_914 + (l_mult_915 + (l_mult_851 + (-0.6696e4 * l_prod_614 + (0.20736e5 * l_prod_613 + (l_mult_916 + (0.20736e5 * l_prod_609 + (l_mult_917 + (l_mult_918 + (0.6480e4 * l_prod_603 + (l_mult_919 + (l_mult_920 + l_mult_921)))))))))))))))))))))))))))))))))))))))))))))))))))))),
                l_mult_922 + (l_mult_779 + (l_mult_923 + (l_mult_924 + (l_mult_716 + (l_mult_856 + (l_mult_781 + (l_mult_857 + (l_mult_858 + (l_mult_925 + (l_mult_784 + (l_mult_926 + (l_mult_927 + (l_mult_787 + (l_mult_740 + (0.1053e4 * l_prod_652 + (l_mult_928 + (l_mult_789 + (-0.7128e4 * l_prod_649 + (l_mult_680 + (l_mult_929 + (l_mult_930 + (l_mult_792 + (l_mult_681 + (l_mult_865 + (l_mult_866 + (l_mult_682 + (-0.7128e4 * l_prod_636 + (l_mult_683 + (l_mult_684 + (-0.62235e4 * l_prod_630 + (0.23652e5 * l_prod_629 + (l_mult_931 + (l_mult_797 + (0.23652e5 * l_prod_624 + (-0.58320e5 * l_prod_623 + (l_mult_932 + (l_mult_933 + (l_mult_934 + (l_mult_871 + (0.14796e5 * l_prod_614 + (-0.37584e5 * l_prod_613 + (l_mult_935 + (-0.37584e5 * l_prod_609 + (l_mult_936 + (l_mult_937 + (-0.15390e5 * l_prod_603 + (l_mult_938 + (l_mult_939 + l_mult_940)))))))))))))))))))))))))))))))))))))))))))))))),
                l_mult_941 + (l_mult_800 + (l_mult_942 + (l_mult_943 + (l_mult_873 + (l_mult_802 + (l_mult_874 + (l_mult_944 + (l_mult_805 + (l_mult_945 + (-0.1016e4 * l_prod_652 + (l_mult_946 + (l_mult_809 + (0.1296e4 * l_prod_649 + (l_mult_947 + (l_mult_948 + (l_mult_812 + (l_mult_880 + (l_mult_881 + (0.1296e4 * l_prod_636 + (0.6696e4 * l_prod_630 + (-0.18360e5 * l_prod_629 + (l_mult_769 + (l_mult_816 + (-0.18360e5 * l_prod_624 + (l_mult_949 + (l_mult_848 + (l_mult_849 + (l_mult_773 + (l_mult_885 + (-0.17424e5 * l_prod_614 + (0.33696e5 * l_prod_613 + (l_mult_916 + (0.33696e5 * l_prod_609 + (l_mult_917 + (l_mult_918 + (0.19440e5 * l_prod_603 + (l_mult_950 + (l_mult_951 + -0.7776e4 * l_prod_598)))))))))))))))))))))))))))))))))))))))) + (vdot4(
            v_578,
            vcons4(
                l_mult_952 + (l_mult_819 + (l_mult_700 + (l_mult_888 + (l_mult_821 + (l_mult_724 + (0.594e3 * l_prod_652 + (l_mult_953 + (l_mult_702 + (l_mult_954 + (l_mult_955 + (l_mult_726 + (-0.41445e4 * l_prod_630 + (0.7128e4 * l_prod_629 + (l_mult_703 + (0.7128e4 * l_prod_624 + (l_mult_956 + (l_mult_727 + (0.11556e5 * l_prod_614 + (-0.14904e5 * l_prod_613 + (l_mult_689 + (-0.14904e5 * l_prod_609 + (l_mult_690 + (l_mult_691 + (-0.13770e5 * l_prod_603 + (l_mult_957 + (l_mult_958 + l_mult_940)))))))))))))))))))))))))),
                l_mult_959 + (l_mult_825 + (l_mult_893 + (-0.1944e3 * l_prod_652 + (l_mult_960 + (l_mult_961 + (0.1404e4 * l_prod_630 + (-0.1134e4 * l_prod_629 + (-0.1134e4 * l_prod_624 + (-0.4104e4 * l_prod_614 + (l_mult_962 + (l_mult_963 + (l_mult_964 + (l_mult_777 + (l_mult_854 + l_mult_921)))))))))))))),
                l_mult_821 + (l_mult_955 + (l_mult_956 + l_mult_690)),
                l_mult_965 + (l_mult_966 + (l_mult_967 + (l_mult_968 + (l_mult_798 + l_mult_687)))))) + (vdot4(v_579,
            vcons4(l_mult_965 + (l_mult_969 + (l_mult_787 + l_sum_971)),
                l_mult_821 + (l_mult_889 + (l_mult_890 + l_mult_678)),
                l_mult_965 + (l_mult_972 + (l_mult_967 + (l_mult_973 + l_sum_974))),
                l_mult_975 + (l_mult_976 + (l_mult_977 + (l_mult_926 + l_sum_980))))) + (vdot4(v_580,
            vcons4(l_mult_965 + (l_mult_972 + (l_mult_969 + (l_mult_981 + l_sum_982))),
                l_mult_965 + (l_mult_983 + (l_mult_858 + l_sum_984)),
                l_mult_965 + (l_mult_983 + (l_mult_858 + l_sum_985)), l_sum_986)) + (vdot4(v_581,
            vcons4(
                l_mult_987 + (l_mult_988 + (l_mult_989 + (l_mult_990 + (l_mult_991 + (-0.34992e5 * l_prod_660 + (l_mult_992 + (l_mult_993 + (l_mult_994 + (l_mult_995 + (l_mult_996 + (l_mult_997 + (l_mult_793 + (l_mult_998 + (l_mult_999 + (l_mult_867 + (-0.17496e5 * l_prod_623 + (l_mult_1000 + (l_mult_1001 + l_mult_690)))))))))))))))))),
                l_mult_1002 + (l_mult_1003 + (l_mult_1004 + (l_mult_1005 + (l_mult_1006 + (l_mult_754 + (l_mult_1007 + (l_mult_1008 + (l_mult_1009 + (l_mult_1010 + (l_mult_812 + (l_mult_1011 + (l_mult_766 + (l_mult_845 + (l_mult_1012 + l_mult_773)))))))))))))),
                l_mult_1013 + (l_mult_1014 + (l_mult_1015 + (l_mult_1016 + (l_mult_1017 + (l_mult_1018 + (l_mult_995 + l_sum_971)))))),
                l_mult_1019 + (l_mult_1020 + (l_mult_1021 + l_mult_757)))) + (vdot4(v_582,
            vcons4(
                l_mult_1002 + (l_mult_1022 + (l_mult_1023 + (l_mult_1024 + (l_mult_1025 + (l_mult_1006 + (l_mult_1026 + (l_mult_1021 + (l_mult_836 + (l_mult_1010 + (l_mult_1027 + (l_mult_764 + (l_mult_881 + (l_mult_766 + (l_mult_1012 + l_mult_848)))))))))))))),
                l_mult_1028 + (l_mult_1029 + (l_mult_823 + (l_mult_1030 + (-0.14580e5 * l_prod_660 + (l_mult_786 + (l_mult_890 + (l_mult_860 + l_sum_980))))))),
                l_mult_1031 + (l_mult_874 + (l_mult_1032 + (l_mult_806 + (l_mult_1033 + l_mult_877)))),
                l_mult_1013 + (l_mult_1034 + (l_mult_1035 + (l_mult_990 + (l_mult_1036 + (l_mult_1016 + (l_mult_1037 + l_sum_984)))))))) + (vdot4(
            v_583,
            vcons4(l_mult_1031 + (l_mult_1038 + (l_mult_1039 + (l_mult_805 + (l_mult_806 + l_mult_807)))),
                l_mult_1019 + (l_mult_1040 + (l_mult_1004 + l_mult_832)),
                l_mult_1041 + (l_mult_1042 + (l_mult_1043 + (l_mult_1044 + (l_mult_820 + (l_mult_987 + (l_mult_988 + (l_mult_989 + (l_mult_990 + (0.6426e4 * l_prod_661 + (l_mult_1045 + (l_mult_786 + (-0.5832e4 * l_prod_657 + (l_mult_1018 + (l_mult_678 + (l_mult_1046 + (0.25704e5 * l_prod_650 + (-0.34992e5 * l_prod_649 + (l_mult_1047 + (l_mult_1048 + (l_mult_1049 + (l_mult_1050 + (l_mult_998 + (l_mult_999 + (l_mult_1051 + (0.19278e5 * l_prod_629 + (l_mult_1052 + (l_mult_1053 + (l_mult_1054 + (l_mult_1055 + (l_mult_934 + (l_mult_1056 + (l_mult_1057 + (l_mult_1058 + l_mult_957))))))))))))))))))))))))))))))))),
                l_mult_1059 + (l_mult_1060 + (l_mult_1061 + (l_mult_826 + (l_mult_1002 + (l_mult_1003 + (l_mult_1004 + (l_mult_1062 + (l_mult_806 + (l_mult_1033 + (l_mult_1063 + (l_mult_1064 + (0.23328e5 * l_prod_649 + (l_mult_761 + (l_mult_1065 + (l_mult_1066 + (l_mult_764 + (l_mult_1011 + (l_mult_766 + (l_mult_767 + (-0.28836e5 * l_prod_629 + (l_mult_1067 + (l_mult_1068 + (l_mult_1069 + (l_mult_1070 + (l_mult_915 + (0.41472e5 * l_prod_613 + (-0.46656e5 * l_prod_612 + (l_mult_1071 + l_mult_950)))))))))))))))))))))))))))))) + (vdot4(
            v_584,
            vcons4(
                l_mult_1072 + (l_mult_1073 + (l_mult_1074 + (l_mult_1013 + (l_mult_1014 + (l_mult_966 + (l_mult_1075 + (l_mult_1076 + (l_mult_1077 + (l_mult_1078 + (l_mult_1079 + (l_mult_968 + (0.21060e5 * l_prod_629 + (l_mult_1080 + (l_mult_797 + (l_mult_1081 + (l_mult_1000 + (l_mult_687 + (-0.36288e5 * l_prod_613 + (l_mult_1057 + (l_mult_1058 + l_mult_938)))))))))))))))))))),
                l_mult_1082 + (l_mult_1083 + (l_mult_1019 + (l_mult_1084 + (l_mult_1085 + (l_mult_1086 + (-0.7614e4 * l_prod_629 + (l_mult_1087 + (l_mult_1088 + (l_mult_1089 + (l_mult_775 + (l_mult_776 + l_mult_919))))))))))),
                l_mult_1059 + (l_mult_1090 + (l_mult_1091 + (l_mult_1092 + (l_mult_801 + (l_mult_1002 + (l_mult_1022 + (l_mult_1023 + (l_mult_1024 + (l_mult_1062 + (l_mult_1093 + (l_mult_754 + (l_mult_1033 + (l_mult_877 + (l_mult_1094 + (l_mult_1064 + (0.42768e5 * l_prod_649 + (-0.23328e5 * l_prod_648 + (l_mult_1095 + (l_mult_1066 + (l_mult_1096 + (l_mult_881 + (l_mult_766 + (-0.4860e4 * l_prod_629 + (l_mult_1097 + (l_mult_1068 + (l_mult_1088 + (l_mult_913 + l_sum_1098))))))))))))))))))))))))))),
                l_mult_1099 + (l_mult_1100 + (l_mult_1101 + (l_mult_1102 + (l_mult_1028 + (l_mult_1029 + (l_mult_823 + (l_mult_977 + (l_mult_926 + (l_mult_1103 + (0.17496e5 * l_prod_650 + (-0.27216e5 * l_prod_649 + (l_mult_790 + (l_mult_1104 + (l_mult_1105 + (l_mult_793 + (l_mult_794 + (l_mult_682 + (0.5832e4 * l_prod_629 + (l_mult_1106 + (l_mult_1053 + (l_mult_956 + (l_mult_932 + (l_mult_1107 + l_mult_935))))))))))))))))))))))))) + (vdot4(
            v_585,
            vcons4(
                l_mult_1108 + (l_mult_1109 + (l_mult_1110 + (l_mult_1031 + (l_mult_874 + (l_mult_1111 + (l_mult_1112 + (l_mult_1113 + (l_mult_1114 + (l_mult_812 + (-0.2916e4 * l_prod_629 + (l_mult_1115 + (l_mult_770 + (l_mult_1012 + (l_mult_848 + l_sum_1098)))))))))))))),
                l_mult_1072 + (l_mult_1116 + (l_mult_1117 + (l_mult_1118 + (l_mult_780 + (l_mult_1013 + (l_mult_1034 + (l_mult_1035 + (l_mult_990 + (l_mult_966 + (l_mult_981 + (l_mult_676 + (l_mult_1119 + (l_mult_1076 + (-0.22032e5 * l_prod_649 + (l_mult_1047 + (l_mult_1120 + (l_mult_1079 + (l_mult_1121 + (l_mult_1122 + (l_mult_1123 + l_mult_797)))))))))))))))))))),
                l_mult_1108 + (l_mult_1124 + (l_mult_1125 + (l_mult_826 + (l_mult_1031 + (l_mult_1038 + (l_mult_1039 + (l_mult_1126 + (l_mult_1112 + (0.12960e5 * l_prod_649 + (l_mult_761 + (l_mult_811 + (l_mult_812 + (l_mult_813 + (-0.648e3 * l_prod_629 + (l_mult_1087 + l_mult_770))))))))))))))),
                l_mult_1082 + (l_mult_1127 + (l_mult_1128 + (l_mult_1129 + (l_mult_747 + (l_mult_1019 + (l_mult_1040 + (l_mult_1004 + (l_mult_832 + (l_mult_1130 + (l_mult_1085 + (l_mult_1113 + l_mult_905))))))))))))) + (vdot4(
            v_586,
            vcons4(
                l_mult_1131 + (l_mult_987 + (0.6426e4 * l_prod_665 + (-0.5832e4 * l_prod_664 + (l_mult_675 + (l_mult_1132 + (l_mult_991 + (l_mult_1045 + (l_mult_1037 + (l_mult_1133 + (l_mult_993 + (l_mult_860 + (l_mult_1134 + (l_mult_995 + (l_mult_891 + (l_mult_1135 + (l_mult_1048 + (l_mult_997 + (l_mult_1121 + (0.25704e5 * l_prod_640 + (l_mult_1136 + (l_mult_999 + (-0.34992e5 * l_prod_636 + (l_mult_1137 + (l_mult_1138 + (0.19278e5 * l_prod_624 + (l_mult_1054 + (l_mult_932 + (l_mult_1139 + (l_mult_1140 + (l_mult_1141 + (l_mult_1142 + (l_mult_1058 + (l_mult_1143 + l_mult_958))))))))))))))))))))))))))))))))),
                l_mult_1144 + (l_mult_1002 + (l_mult_1145 + (l_mult_1039 + (l_mult_1146 + (l_mult_1025 + (l_mult_806 + (l_mult_1147 + (l_mult_1021 + (l_mult_894 + (l_mult_1148 + (l_mult_1065 + (l_mult_1027 + (l_mult_813 + (l_mult_1149 + (l_mult_1150 + (l_mult_766 + (0.23328e5 * l_prod_636 + (l_mult_845 + (l_mult_846 + (-0.28836e5 * l_prod_624 + (l_mult_1069 + (l_mult_913 + (l_mult_1151 + (l_mult_1152 + (l_mult_1153 + (0.41472e5 * l_prod_609 + (l_mult_1071 + (-0.46656e5 * l_prod_606 + l_mult_951)))))))))))))))))))))))))))),
                l_mult_1154 + (l_mult_1013 + (l_mult_972 + (l_mult_1155 + (l_mult_1036 + (l_mult_1156 + (l_mult_1157 + (l_mult_1078 + (l_mult_973 + (l_mult_1158 + (l_mult_1159 + (l_mult_1160 + (0.21060e5 * l_prod_624 + (l_mult_1081 + (l_mult_686 + (l_mult_1161 + (l_mult_1001 + (l_mult_871 + (-0.36288e5 * l_prod_609 + (l_mult_1058 + (l_mult_1143 + l_mult_939)))))))))))))))))))),
                l_mult_1162 + (l_mult_1019 + (l_mult_1163 + (l_mult_1164 + (l_mult_1086 + (l_mult_1165 + (-0.7614e4 * l_prod_624 + (l_mult_1088 + (l_mult_1166 + (l_mult_1167 + (l_mult_776 + (l_mult_853 + l_mult_920))))))))))))) + (vdot4(
            v_587,
            vcons4(
                l_mult_1144 + (l_mult_1002 + (l_mult_1145 + (l_mult_1039 + (l_mult_1168 + (l_mult_1005 + (l_mult_1093 + (l_mult_807 + (l_mult_1169 + (l_mult_1007 + (l_mult_836 + (l_mult_1170 + (l_mult_1009 + (l_mult_878 + (l_mult_1171 + (l_mult_1095 + (l_mult_812 + (l_mult_1149 + (l_mult_1150 + (l_mult_766 + (0.42768e5 * l_prod_636 + (l_mult_1172 + (-0.23328e5 * l_prod_633 + (-0.4860e4 * l_prod_624 + (l_mult_1088 + (l_mult_1173 + (l_mult_915 + (l_mult_1153 + l_sum_1174))))))))))))))))))))))))))),
                l_mult_1175 + (l_mult_1028 + (l_mult_976 + (l_mult_1176 + (l_mult_1030 + (l_mult_926 + (l_mult_1177 + (l_mult_890 + (l_mult_1178 + (l_mult_1179 + (l_mult_1104 + (l_mult_864 + (0.17496e5 * l_prod_640 + (l_mult_1180 + (l_mult_682 + (-0.27216e5 * l_prod_636 + (l_mult_867 + (l_mult_868 + (0.5832e4 * l_prod_624 + (l_mult_956 + (l_mult_1181 + (l_mult_934 + (l_mult_1141 + (l_mult_1182 + l_mult_937))))))))))))))))))))))),
                l_mult_1183 + (l_mult_1031 + (l_mult_1184 + (l_mult_805 + (l_mult_1185 + (l_mult_1186 + (l_mult_1114 + (l_mult_1187 + (l_mult_881 + (l_mult_1188 + (-0.2916e4 * l_prod_624 + (l_mult_1012 + (l_mult_1189 + (l_mult_773 + (l_mult_851 + l_sum_1174)))))))))))))),
                l_mult_1154 + (l_mult_1013 + (l_mult_972 + (l_mult_1190 + (l_mult_1015 + (l_mult_981 + (l_mult_1191 + (l_mult_1017 + (l_mult_677 + (l_mult_1192 + (l_mult_995 + (l_mult_862 + (l_mult_1193 + (l_mult_1120 + (l_mult_1158 + (l_mult_1159 + (-0.22032e5 * l_prod_636 + (l_mult_1051 + (l_mult_1138 + (l_mult_1194 + (l_mult_1195 + l_mult_871)))))))))))))))))))))) + (vdot4(
            v_588,
            vcons4(
                l_mult_1183 + (l_mult_1031 + (l_mult_1196 + (l_mult_1032 + (l_mult_1197 + (l_mult_1033 + (l_mult_894 + (l_mult_1198 + (l_mult_811 + (l_mult_1187 + (l_mult_881 + (0.12960e5 * l_prod_636 + (l_mult_767 + (l_mult_846 + (-0.648e3 * l_prod_624 + (l_mult_1166 + l_mult_851))))))))))))))),
                l_mult_1162 + (l_mult_1019 + (l_mult_1199 + (l_mult_1020 + (l_mult_1200 + (l_mult_1021 + (l_mult_1201 + (l_mult_757 + (l_mult_839 + (l_mult_1202 + (l_mult_1165 + (l_mult_1188 + l_mult_910))))))))))),
                0.4320e4 * l_prod_666 + (-0.15984e5 * l_prod_665 + (0.19440e5 * l_prod_664 + (l_mult_751 + (-0.15984e5 * l_prod_661 + (0.38880e5 * l_prod_660 + (l_mult_1026 + (0.19440e5 * l_prod_657 + (l_mult_1008 + (l_mult_838 + (l_mult_1203 + (l_mult_1204 + (l_mult_1096 + (l_mult_1205 + (-0.93312e5 * l_prod_638 + (l_mult_1172 + (0.58320e5 * l_prod_623 + (l_mult_1070 + (l_mult_1152 + l_mult_917)))))))))))))))))),
                l_mult_1206 + (l_mult_1207 + (l_mult_823 + (l_mult_1208 + (l_mult_1016 + (l_mult_890 + (l_mult_1209 + (l_mult_1210 + (l_mult_793 + (l_mult_1211 + (l_mult_999 + (l_mult_867 + (l_mult_1212 + (l_mult_1055 + (l_mult_1140 + l_mult_936)))))))))))))))) + (vdot4(
            v_589,
            vcons4(
                l_mult_1213 + (l_mult_1214 + (l_mult_1215 + (l_mult_1216 + (l_mult_1217 + (l_mult_1218 + (0.34992e5 * l_prod_623 + (l_mult_772 + (l_mult_850 + l_mult_917)))))))),
                l_mult_1206 + (l_mult_1207 + (l_mult_823 + (0.13284e5 * l_prod_661 + (l_mult_1219 + (l_mult_786 + (l_mult_1220 + (l_mult_994 + (l_mult_861 + (l_mult_1221 + (l_mult_1079 + (l_mult_1211 + (l_mult_999 + (l_mult_1137 + (l_mult_956 + l_mult_934)))))))))))))),
                l_mult_1222 + (l_mult_1038 + (l_mult_1223 + (l_mult_806 + (l_mult_1021 + (l_mult_1224 + (l_mult_812 + (l_mult_843 + (l_mult_766 + (l_mult_845 + (l_mult_1088 + l_mult_915)))))))))),
                l_mult_1213 + (l_mult_1214 + (-0.4320e4 * l_prod_661 + (l_mult_1225 + (0.11664e5 * l_prod_657 + (l_mult_756 + (l_mult_838 + (l_mult_1226 + (l_mult_1218 + l_mult_909)))))))))) + (vdot4(
            v_590,
            vcons4(
                l_mult_1206 + (0.13284e5 * l_prod_665 + (l_mult_1227 + (l_mult_783 + (l_mult_1208 + (l_mult_1219 + (l_mult_992 + (l_mult_890 + (l_mult_860 + (l_mult_1221 + (l_mult_1210 + (l_mult_1050 + (l_mult_1159 + (l_mult_999 + (l_mult_956 + l_mult_932)))))))))))))),
                l_mult_1222 + (l_mult_1228 + (l_mult_1004 + (l_mult_1032 + (l_mult_806 + (l_mult_1224 + (l_mult_763 + (l_mult_764 + (l_mult_881 + (l_mult_766 + (l_mult_1088 + l_mult_913)))))))))),
                l_mult_1222 + (l_mult_1228 + (l_mult_1004 + (l_mult_1223 + (l_mult_753 + (l_mult_754 + (l_mult_1021 + (l_mult_836 + (l_mult_1114 + (l_mult_812 + (l_mult_881 + l_mult_766)))))))))),
                l_mult_1213 + (-0.4320e4 * l_prod_665 + (0.11664e5 * l_prod_664 + (l_mult_751 + (l_mult_1215 + (l_mult_1225 + (l_mult_834 + (l_mult_1226 + (l_mult_1217 + l_mult_908)))))))))) + vdot4(
            v_570,
            vcons4(l_sum_695,
                l_mult_696 + (0.274e2 * l_prod_652 + (-0.2025e3 * l_prod_630 + (0.612e3 * l_prod_614 + (-0.810e3 * l_prod_603 + l_mult_694)))),
                0.e0, 0.e0)))))))))))))))))))));
        double l_vdot_1613 = vdot4(v_571,
            vcons4(l_mult_697 + (l_mult_1229 + (0.1134e4 * l_prod_661 + (-0.2592e4 * l_prod_657 + l_mult_678))),
                l_mult_699 + (l_mult_700 + (l_mult_1230 + (l_mult_822 + (0.486e3 * l_prod_661 + (l_mult_1231 + l_sum_982))))),
                l_mult_706 + (l_mult_707 + (l_mult_708 + (l_mult_1232 + (l_mult_1233 + (-0.1296e4 * l_prod_664 + l_sum_985))))),
                l_mult_699 + (l_mult_714 + (l_mult_715 + (l_mult_716 + l_sum_986))))) + (vdot4(v_572,
            vcons4(l_sum_720, 0.e0, 0.e0, 0.e0)) + (vdot4(v_573,
            vcons4(0.e0, 0.e0, l_sum_1235,
                l_mult_1236 + (l_mult_741 + (l_mult_1237 + (l_mult_1238 + (l_mult_1239 + (l_mult_1240 + l_sum_1242))))))) + (vdot4(
            v_574,
            vcons4(
                l_mult_1243 + (l_mult_733 + (l_mult_1244 + (l_mult_1245 + (l_mult_1194 + (l_mult_736 + (l_mult_1246 + (-0.1296e4 * l_prod_609 + l_mult_691))))))),
                l_mult_1236 + (l_mult_725 + (0.486e3 * l_prod_640 + (l_mult_1247 + l_sum_1249))),
                l_mult_1234 + (l_mult_722 + (0.1134e4 * l_prod_640 + (-0.2592e4 * l_prod_636 + l_mult_684))),
                l_sum_778)) + (vdot4(v_575, vcons4(l_sum_799, l_sum_818, l_sum_824, l_sum_828)) + (vdot4(v_576,
            vcons4(
                l_mult_897 + (l_mult_745 + (l_mult_898 + (l_mult_899 + (l_mult_900 + (l_mult_827 + (-0.6264e3 * l_prod_667 + (l_mult_1250 + (l_mult_749 + (l_mult_804 + (l_mult_1251 + (0.3132e4 * l_prod_662 + (-0.15066e5 * l_prod_661 + (l_mult_1006 + (l_mult_754 + (-0.6696e4 * l_prod_658 + (0.20736e5 * l_prod_657 + (l_mult_1252 + (0.6480e4 * l_prod_655 + (l_mult_1253 + (l_mult_1254 + (l_mult_1255 + (l_mult_758 + (l_mult_1256 + (l_mult_1257 + (l_mult_1258 + (l_mult_906 + (l_mult_907 + (l_mult_763 + (l_mult_908 + (-0.15066e5 * l_prod_640 + (l_mult_1150 + (l_mult_1259 + (0.20736e5 * l_prod_636 + (l_mult_1260 + (l_mult_1261 + (l_mult_1262 + (l_mult_768 + (0.7776e4 * l_prod_628 + (l_mult_816 + (l_mult_1263 + (l_mult_949 + (l_mult_772 + (l_mult_914 + (l_mult_915 + (l_mult_1264 + (l_mult_1265 + (l_mult_774 + (l_mult_1266 + (l_mult_1267 + (l_mult_1268 + (l_mult_1269 + (l_mult_1270 + (l_mult_777 + (l_mult_1271 + l_mult_1272)))))))))))))))))))))))))))))))))))))))))))))))))))))),
                l_mult_922 + (l_mult_779 + (l_mult_923 + (l_mult_924 + (l_mult_716 + (0.1053e4 * l_prod_667 + (l_mult_1273 + (l_mult_782 + (-0.7128e4 * l_prod_664 + (l_mult_675 + (-0.62235e4 * l_prod_662 + (0.23652e5 * l_prod_661 + (l_mult_1274 + (l_mult_786 + (0.14796e5 * l_prod_658 + (-0.37584e5 * l_prod_657 + (l_mult_994 + (-0.15390e5 * l_prod_655 + (l_mult_1275 + (l_mult_1276 + (l_mult_1277 + (l_mult_788 + (l_mult_1278 + (l_mult_1279 + (l_mult_929 + (l_mult_930 + (l_mult_792 + (l_mult_681 + (0.23652e5 * l_prod_640 + (-0.58320e5 * l_prod_639 + (l_mult_1280 + (-0.37584e5 * l_prod_636 + (l_mult_1137 + (l_mult_1281 + (l_mult_1282 + (l_mult_795 + (l_mult_1283 + (l_mult_1284 + (l_mult_1285 + (l_mult_686 + (l_mult_933 + (l_mult_934 + (l_mult_1286 + (l_mult_1287 + (l_mult_704 + (-0.7128e4 * l_prod_609 + (l_mult_690 + (l_mult_1288 + l_sum_1242))))))))))))))))))))))))))))))))))))))))))))))),
                l_mult_941 + (l_mult_800 + (l_mult_942 + (l_mult_943 + (-0.1016e4 * l_prod_667 + (l_mult_1289 + (l_mult_803 + (0.1296e4 * l_prod_664 + (0.6696e4 * l_prod_662 + (-0.18360e5 * l_prod_661 + (l_mult_753 + (l_mult_807 + (-0.17424e5 * l_prod_658 + (0.33696e5 * l_prod_657 + (l_mult_1252 + (0.19440e5 * l_prod_655 + (l_mult_1290 + (-0.7776e4 * l_prod_653 + (l_mult_1291 + (l_mult_808 + (l_mult_1292 + (l_mult_947 + (l_mult_948 + (l_mult_812 + (-0.18360e5 * l_prod_640 + (l_mult_843 + (l_mult_1293 + (0.33696e5 * l_prod_636 + (l_mult_1260 + (l_mult_1294 + (l_mult_1295 + (l_mult_814 + (l_mult_1296 + (l_mult_1297 + (l_mult_849 + (l_mult_773 + (l_mult_1264 + (l_mult_1298 + (0.1296e4 * l_prod_609 + l_mult_1299)))))))))))))))))))))))))))))))))))))),
                l_mult_952 + (l_mult_819 + (l_mult_700 + (0.594e3 * l_prod_667 + (l_mult_1300 + (l_mult_822 + (-0.41445e4 * l_prod_662 + (0.7128e4 * l_prod_661 + (l_mult_1231 + (0.11556e5 * l_prod_658 + (-0.14904e5 * l_prod_657 + (l_mult_677 + (-0.13770e5 * l_prod_655 + (l_mult_1301 + (l_mult_1276 + (l_mult_1302 + (l_mult_717 + (l_mult_954 + (l_mult_955 + (0.7128e4 * l_prod_640 + (l_mult_1303 + (-0.14904e5 * l_prod_636 + (l_mult_683 + (l_mult_1304 + l_sum_1249))))))))))))))))))))))))) + (vdot4(
            v_577,
            vcons4(
                l_mult_959 + (l_mult_825 + (-0.1944e3 * l_prod_667 + (l_mult_1305 + (0.1404e4 * l_prod_662 + (-0.1134e4 * l_prod_661 + (-0.4104e4 * l_prod_658 + (l_mult_1306 + (l_mult_837 + (l_mult_757 + (l_mult_1254 + (l_mult_1307 + (l_mult_961 + (-0.1134e4 * l_prod_640 + (l_mult_1308 + l_mult_1309)))))))))))))),
                l_sum_1320, l_sum_1328, l_sum_1333)) + (vdot4(v_578,
            vcons4(l_sum_1337, l_sum_1339, l_mult_717 + (l_mult_1334 + (l_mult_1335 + l_mult_692)),
                l_mult_1340 + (l_mult_970 + (l_mult_1341 + (l_mult_1342 + (l_mult_704 + l_mult_690)))))) + (vdot4(
            v_579,
            vcons4(l_mult_1340 + (l_mult_967 + (l_mult_794 + l_sum_1343)),
                l_mult_717 + (l_mult_955 + (l_mult_1303 + l_mult_683)),
                l_mult_1340 + (l_mult_1344 + (l_mult_1341 + (l_mult_712 + l_sum_705))),
                l_mult_1345 + (l_mult_1346 + (l_mult_978 + (l_mult_864 + l_sum_1347))))) + (vdot4(v_580,
            vcons4(l_mult_1340 + (l_mult_1344 + (l_mult_967 + (l_mult_973 + l_sum_979))),
                l_mult_1340 + (l_mult_1348 + (l_mult_1279 + l_sum_713)),
                l_mult_1340 + (l_mult_1348 + (l_mult_1279 + l_sum_984)), l_sum_719)) + (vdot4(v_581,
            vcons4(
                l_mult_1041 + (l_mult_1042 + (l_mult_1043 + (l_mult_1044 + (l_mult_820 + (l_mult_1349 + (0.25704e5 * l_prod_665 + (-0.34992e5 * l_prod_664 + (l_mult_1350 + (0.19278e5 * l_prod_661 + (l_mult_1351 + (l_mult_1352 + (l_mult_1220 + (l_mult_1353 + (l_mult_1301 + (l_mult_1354 + (l_mult_1355 + (l_mult_1356 + (l_mult_1357 + (l_mult_1048 + (l_mult_1049 + (l_mult_1050 + (l_mult_1358 + (l_mult_1359 + (l_mult_1360 + (0.6426e4 * l_prod_629 + (l_mult_1361 + (l_mult_797 + (l_mult_1362 + (l_mult_1363 + (l_mult_934 + (-0.5832e4 * l_prod_613 + (l_mult_1364 + (l_mult_1365 + l_mult_692))))))))))))))))))))))))))))))))),
                l_mult_1059 + (l_mult_1060 + (l_mult_1061 + (l_mult_826 + (l_mult_1366 + (l_mult_1367 + (0.23328e5 * l_prod_664 + (l_mult_751 + (-0.28836e5 * l_prod_661 + (l_mult_1368 + (l_mult_1369 + (0.41472e5 * l_prod_657 + (-0.46656e5 * l_prod_656 + (l_mult_1290 + (l_mult_1370 + (l_mult_1371 + (l_mult_1372 + (l_mult_1065 + (l_mult_1066 + (l_mult_764 + (l_mult_1373 + (l_mult_1374 + (l_mult_1172 + (l_mult_1375 + (l_mult_815 + (l_mult_1376 + (l_mult_772 + (l_mult_915 + (l_mult_1377 + l_mult_776)))))))))))))))))))))))))))),
                l_mult_1072 + (l_mult_1073 + (l_mult_1074 + (l_mult_1378 + (l_mult_1379 + (l_mult_1380 + (0.21060e5 * l_prod_661 + (l_mult_1381 + (l_mult_786 + (-0.36288e5 * l_prod_657 + (l_mult_1353 + (l_mult_1275 + (l_mult_1382 + (l_mult_1383 + (l_mult_1078 + (l_mult_1079 + (l_mult_1384 + (l_mult_1385 + (l_mult_1360 + l_sum_1343)))))))))))))))))),
                l_mult_1082 + (l_mult_1083 + (l_mult_1386 + (l_mult_1387 + (-0.7614e4 * l_prod_661 + (l_mult_1388 + (l_mult_835 + (l_mult_756 + (l_mult_1253 + (l_mult_1389 + (l_mult_1086 + (l_mult_1390 + l_mult_767))))))))))))) + (vdot4(
            v_582,
            vcons4(
                l_mult_1059 + (l_mult_1090 + (l_mult_1091 + (l_mult_1092 + (l_mult_801 + (l_mult_1391 + (l_mult_1367 + (0.42768e5 * l_prod_664 + (-0.23328e5 * l_prod_663 + (-0.4860e4 * l_prod_661 + (l_mult_1392 + (l_mult_1369 + (l_mult_1306 + (l_mult_1252 + (l_mult_1370 + (l_mult_1393 + (l_mult_1394 + (l_mult_1395 + (l_mult_1095 + (l_mult_1066 + (l_mult_1096 + (l_mult_1390 + (l_mult_1259 + (l_mult_1375 + (l_mult_1396 + (l_mult_770 + (l_mult_1297 + (l_mult_772 + l_sum_1397))))))))))))))))))))))))))),
                l_mult_1099 + (l_mult_1100 + (l_mult_1101 + (l_mult_1102 + (l_mult_1398 + (0.17496e5 * l_prod_665 + (-0.27216e5 * l_prod_664 + (l_mult_783 + (0.5832e4 * l_prod_661 + (l_mult_1399 + (l_mult_1352 + (l_mult_1400 + (l_mult_994 + (l_mult_1401 + (l_mult_1402 + (l_mult_718 + (l_mult_1104 + (l_mult_1105 + (l_mult_793 + (l_mult_1303 + (l_mult_1280 + l_sum_1347)))))))))))))))))))),
                l_mult_1108 + (l_mult_1109 + (l_mult_1110 + (l_mult_1403 + (l_mult_1404 + (l_mult_1405 + (-0.2916e4 * l_prod_661 + (l_mult_1406 + (l_mult_754 + (l_mult_1306 + (l_mult_1252 + (l_mult_1407 + (l_mult_1292 + (l_mult_1114 + (l_mult_812 + (l_mult_1408 + l_mult_1293))))))))))))))),
                l_mult_1072 + (l_mult_1116 + (l_mult_1117 + (l_mult_1118 + (l_mult_780 + (l_mult_1409 + (l_mult_1379 + (-0.22032e5 * l_prod_664 + (l_mult_1350 + (l_mult_1410 + (l_mult_1411 + (l_mult_786 + (l_mult_1382 + (l_mult_1412 + (l_mult_1413 + (l_mult_1357 + (l_mult_1120 + (l_mult_1079 + (l_mult_1121 + l_sum_713)))))))))))))))))))) + (vdot4(
            v_583,
            vcons4(
                l_mult_1108 + (l_mult_1124 + (l_mult_1125 + (l_mult_826 + (l_mult_1414 + (l_mult_1404 + (0.12960e5 * l_prod_664 + (l_mult_751 + (-0.648e3 * l_prod_661 + (l_mult_1388 + (l_mult_754 + (l_mult_1407 + (l_mult_1415 + (l_mult_1416 + (l_mult_811 + (l_mult_812 + l_mult_813))))))))))))))),
                l_mult_1082 + (l_mult_1127 + (l_mult_1128 + (l_mult_1129 + (l_mult_747 + (l_mult_1417 + (l_mult_1387 + (l_mult_1405 + (l_mult_1251 + l_sum_1418)))))))),
                l_mult_1354 + (l_mult_1355 + (l_mult_1356 + (l_mult_1357 + (l_mult_996 + (l_mult_997 + (l_mult_793 + (-0.17496e5 * l_prod_639 + (l_mult_1385 + (l_mult_683 + (l_mult_1419 + (-0.34992e5 * l_prod_628 + (l_mult_1420 + (l_mult_1362 + (l_mult_1363 + (l_mult_1001 + (l_mult_1421 + (l_mult_935 + (l_mult_1324 + l_mult_1422)))))))))))))))))),
                l_mult_1370 + (l_mult_1371 + (l_mult_1372 + (l_mult_1010 + (l_mult_812 + (l_mult_1408 + (l_mult_1423 + (l_mult_911 + (l_mult_770 + (l_mult_1376 + (l_mult_772 + (l_mult_773 + (l_mult_1424 + (l_mult_1425 + (l_mult_1316 + l_mult_1426)))))))))))))))) + (vdot4(
            v_584,
            vcons4(
                l_mult_1382 + (l_mult_1383 + (l_mult_970 + (l_mult_1427 + (l_mult_1428 + (l_mult_1342 + (l_mult_1429 + (l_mult_1364 + (l_mult_690 + l_mult_1422)))))))),
                l_mult_1389 + (l_mult_1430 + (l_mult_1431 + l_mult_777)),
                l_mult_1370 + (l_mult_1393 + (l_mult_1394 + (l_mult_1395 + (l_mult_1010 + (l_mult_1027 + (l_mult_764 + (l_mult_1408 + (l_mult_1293 + (l_mult_1432 + (l_mult_911 + (l_mult_1433 + l_sum_1434))))))))))),
                l_mult_1401 + (l_mult_1402 + (l_mult_718 + (l_mult_978 + (l_mult_864 + (l_mult_1435 + (-0.14580e5 * l_prod_628 + (l_mult_797 + (l_mult_798 + (l_mult_686 + l_sum_1436))))))))))) + (vdot4(
            v_585,
            vcons4(l_mult_1407 + (l_mult_1292 + (l_mult_1437 + (l_mult_815 + l_sum_1397))),
                l_mult_1382 + (l_mult_1412 + (l_mult_1413 + (l_mult_1357 + (l_mult_970 + (l_mult_973 + (l_mult_681 + (l_mult_1438 + (l_mult_1428 + l_mult_1439)))))))),
                l_mult_1407 + (l_mult_1415 + (l_mult_1416 + l_sum_817)), l_sum_1418)) + (vdot4(v_586,
            vcons4(
                l_mult_1440 + (l_mult_1354 + (0.6426e4 * l_prod_650 + (-0.5832e4 * l_prod_649 + (l_mult_680 + (l_mult_1135 + (l_mult_1048 + (l_mult_997 + (l_mult_1121 + (0.19278e5 * l_prod_640 + (l_mult_1358 + (l_mult_1280 + (l_mult_1441 + (l_mult_1360 + (l_mult_1304 + (l_mult_1442 + (l_mult_1419 + (l_mult_1361 + (l_mult_1439 + (0.25704e5 * l_prod_624 + (l_mult_1212 + (l_mult_1363 + (l_mult_1139 + (l_mult_1140 + (l_mult_1443 + (l_mult_1444 + (l_mult_1421 + (l_mult_1323 + (-0.34992e5 * l_prod_609 + (l_mult_936 + (l_mult_1445 + (l_mult_1446 + (l_mult_1422 + (l_mult_1447 + l_mult_1336))))))))))))))))))))))))))))))))),
                l_mult_1448 + (l_mult_1370 + (l_mult_1449 + (l_mult_1416 + (l_mult_1171 + (l_mult_1095 + (l_mult_812 + (-0.4860e4 * l_prod_640 + (l_mult_1390 + (l_mult_1308 + (l_mult_1450 + (l_mult_1423 + (l_mult_1396 + (l_mult_816 + (l_mult_1451 + (l_mult_912 + (l_mult_772 + (l_mult_1173 + (l_mult_915 + (l_mult_1264 + (l_mult_1452 + (l_mult_1424 + (l_mult_1315 + (0.42768e5 * l_prod_609 + (l_mult_1071 + (l_mult_1453 + (l_mult_1454 + (l_mult_1426 + (-0.23328e5 * l_prod_600 + l_mult_1332)))))))))))))))))))))))))))),
                l_mult_1455 + (l_mult_1382 + (l_mult_1344 + (l_mult_1193 + (l_mult_1120 + (l_mult_734 + (l_mult_1456 + (l_mult_1427 + (l_mult_712 + (l_mult_1457 + (l_mult_1458 + (l_mult_1195 + (l_mult_1459 + (l_mult_1429 + (l_mult_689 + (-0.22032e5 * l_prod_609 + (l_mult_1365 + (l_mult_1288 + (l_mult_1460 + (l_mult_1422 + (l_mult_1447 + l_mult_1327)))))))))))))))))))),
                l_mult_1461 + (l_mult_1389 + (l_mult_1202 + (l_mult_1462 + (l_mult_1430 + (l_mult_1463 + (l_mult_1464 + (l_mult_1431 + (l_mult_1465 + (l_mult_1466 + (l_mult_777 + (l_mult_1271 + l_mult_1319))))))))))))) + (vdot4(
            v_587,
            vcons4(
                l_mult_1448 + (l_mult_1370 + (l_mult_1449 + (l_mult_1416 + (l_mult_1148 + (l_mult_1065 + (l_mult_1027 + (l_mult_813 + (-0.28836e5 * l_prod_640 + (l_mult_1373 + (l_mult_1259 + (0.41472e5 * l_prod_636 + (l_mult_1172 + (l_mult_1294 + (l_mult_1467 + (l_mult_1432 + (l_mult_815 + (l_mult_1451 + (l_mult_912 + (l_mult_772 + (l_mult_1151 + (l_mult_1152 + (-0.46656e5 * l_prod_617 + (l_mult_1468 + (l_mult_1431 + (0.23328e5 * l_prod_609 + (l_mult_1316 + (l_mult_1453 + l_sum_1469))))))))))))))))))))))))))),
                l_mult_1470 + (l_mult_1401 + (l_mult_1346 + (l_mult_1179 + (l_mult_1104 + (l_mult_864 + (0.5832e4 * l_prod_640 + (l_mult_1303 + (l_mult_1160 + (l_mult_1471 + (l_mult_1435 + (l_mult_1283 + (0.17496e5 * l_prod_624 + (l_mult_1472 + (l_mult_686 + (l_mult_1181 + (l_mult_934 + (l_mult_1286 + (l_mult_1473 + (l_mult_1335 + (-0.27216e5 * l_prod_609 + (l_mult_1324 + (l_mult_1445 + (l_mult_1474 + l_mult_1326))))))))))))))))))))))),
                l_mult_1475 + (l_mult_1407 + (l_mult_1198 + (l_mult_811 + (-0.648e3 * l_prod_640 + (l_mult_1476 + (l_mult_1437 + (l_mult_1477 + (l_mult_1297 + (l_mult_1166 + (l_mult_1478 + (l_mult_1377 + (0.12960e5 * l_prod_609 + (l_mult_776 + (l_mult_1269 + l_sum_1469)))))))))))))),
                l_mult_1455 + (l_mult_1382 + (l_mult_1344 + (l_mult_1157 + (l_mult_1078 + (l_mult_973 + (0.21060e5 * l_prod_640 + (l_mult_1384 + (l_mult_682 + (-0.36288e5 * l_prod_636 + (l_mult_1360 + (l_mult_1281 + (l_mult_1479 + (l_mult_1438 + (l_mult_1457 + (l_mult_1458 + (l_mult_1161 + (l_mult_1001 + (l_mult_1443 + (l_mult_1480 + (l_mult_1182 + l_mult_1288)))))))))))))))))))))) + (vdot4(
            v_588,
            vcons4(
                l_mult_1475 + (l_mult_1407 + (l_mult_1186 + (l_mult_1114 + (-0.2916e4 * l_prod_640 + (l_mult_1408 + (l_mult_1308 + (l_mult_1481 + (l_mult_814 + (l_mult_1477 + (l_mult_1297 + (l_mult_1189 + (l_mult_773 + (l_mult_1264 + (l_mult_1482 + (l_mult_1465 + l_mult_1269))))))))))))))),
                l_mult_1461 + (l_mult_1389 + (l_mult_1164 + (l_mult_1086 + (-0.7614e4 * l_prod_640 + (l_mult_1390 + (l_mult_844 + (l_mult_767 + (l_mult_1261 + (l_mult_1483 + (l_mult_1463 + (l_mult_1166 + l_mult_1314))))))))))),
                0.4320e4 * l_prod_651 + (-0.15984e5 * l_prod_650 + (0.19440e5 * l_prod_649 + (l_mult_761 + (l_mult_1203 + (l_mult_1204 + (l_mult_1096 + (0.58320e5 * l_prod_639 + (l_mult_1374 + (l_mult_1260 + (-0.15984e5 * l_prod_629 + (0.38880e5 * l_prod_628 + (l_mult_1433 + (l_mult_1484 + (-0.93312e5 * l_prod_622 + (l_mult_1152 + (0.19440e5 * l_prod_613 + (l_mult_1425 + (l_mult_1071 + l_mult_1317)))))))))))))))))),
                l_mult_1485 + (l_mult_1486 + (l_mult_718 + (l_mult_1221 + (l_mult_1079 + (l_mult_1303 + (0.13284e5 * l_prod_629 + (l_mult_1487 + (l_mult_797 + (l_mult_1488 + (l_mult_1363 + (l_mult_934 + (l_mult_1056 + (l_mult_935 + (l_mult_936 + l_mult_1325)))))))))))))))) + (vdot4(
            v_589,
            vcons4(
                l_mult_1489 + (l_mult_1490 + (l_mult_1226 + (-0.4320e4 * l_prod_629 + (l_mult_1491 + (l_mult_1492 + (0.11664e5 * l_prod_613 + (l_mult_775 + (l_mult_1268 + l_mult_1317)))))))),
                l_mult_1485 + (l_mult_1486 + (l_mult_718 + (l_mult_1209 + (l_mult_1210 + (l_mult_793 + (l_mult_1136 + (l_mult_1359 + (l_mult_1137 + (l_mult_1493 + (l_mult_1428 + (l_mult_1488 + (l_mult_1363 + (l_mult_1140 + (l_mult_1335 + l_mult_1324)))))))))))))),
                l_mult_1494 + (l_mult_1415 + (l_mult_1224 + (l_mult_812 + (l_mult_1390 + (l_mult_1495 + (l_mult_815 + (l_mult_949 + (l_mult_772 + (l_mult_915 + (l_mult_1431 + l_mult_1316)))))))))),
                l_mult_1489 + (l_mult_1490 + (l_mult_1216 + (l_mult_1217 + (0.34992e5 * l_prod_639 + (l_mult_766 + (l_mult_1260 + (l_mult_1496 + (l_mult_1492 + l_mult_850)))))))))) + (vdot4(
            v_590,
            vcons4(
                l_mult_1485 + (0.13284e5 * l_prod_650 + (l_mult_1497 + (l_mult_790 + (l_mult_1221 + (l_mult_1210 + (l_mult_1050 + (l_mult_1303 + (l_mult_1280 + (l_mult_1493 + (l_mult_1487 + (l_mult_1420 + (l_mult_1458 + (l_mult_1363 + l_sum_1436))))))))))))),
                l_mult_1494 + (l_mult_1498 + (l_mult_1372 + (l_mult_1114 + (l_mult_812 + (l_mult_1495 + (l_mult_769 + (l_mult_770 + l_sum_1434))))))),
                l_mult_1494 + (l_mult_1498 + (l_mult_1372 + (l_mult_1224 + (l_mult_763 + (l_mult_764 + (l_mult_1390 + (l_mult_1259 + (l_mult_1437 + (l_mult_815 + (l_mult_1297 + l_mult_772)))))))))),
                l_mult_1489 + (-0.4320e4 * l_prod_650 + (0.11664e5 * l_prod_649 + (l_mult_761 + (l_mult_1226 + (l_mult_1217 + (l_mult_908 + (l_mult_1496 + (l_mult_1491 + l_mult_1313)))))))))) + vdot4(
            v_570,
            vcons4(l_sum_695, 0.e0,
                l_mult_696 + (0.274e2 * l_prod_667 + (-0.2025e3 * l_prod_662 + (0.612e3 * l_prod_658 + (-0.810e3 * l_prod_655 + l_mult_679)))),
                0.e0)))))))))))))))))))));
        double l_vdot_1614 = vdot4(v_571,
            vcons4(l_sum_744,
                l_mult_723 + (l_mult_821 + (l_mult_738 + (l_mult_889 + (l_mult_739 + (l_mult_890 + (l_mult_740 + l_mult_678)))))),
                l_mult_730 + (l_mult_1232 + (l_mult_972 + (l_mult_731 + (l_mult_1410 + (l_mult_981 + (l_mult_732 + (-0.1296e4 * l_prod_657 + l_mult_677))))))),
                l_mult_723 + (l_mult_1230 + (0.486e3 * l_prod_665 + (l_mult_858 + (l_mult_724 + (l_mult_889 + (l_mult_1231 + l_mult_676)))))))) + (vdot4(
            v_572,
            vcons4(l_mult_721 + (l_mult_1229 + (0.1134e4 * l_prod_665 + (-0.2592e4 * l_prod_664 + l_mult_675))),
                l_sum_1235,
                l_mult_1236 + (l_mult_717 + (l_mult_1237 + (l_mult_1334 + (l_mult_1239 + (l_mult_1335 + l_sum_1499))))),
                l_mult_1243 + (l_mult_709 + (l_mult_1344 + (l_mult_1245 + (l_mult_1122 + (l_mult_712 + (l_mult_1246 + (-0.1296e4 * l_prod_613 + l_mult_689))))))))) + (vdot4(
            v_573,
            vcons4(l_mult_1236 + (l_mult_701 + (0.486e3 * l_prod_650 + (l_mult_1279 + l_sum_1500))),
                l_mult_1234 + (l_mult_698 + (0.1134e4 * l_prod_650 + (-0.2592e4 * l_prod_649 + l_mult_680))), 0.e0,
                0.e0)) + (vdot4(v_574,
            vcons4(0.e0, 0.e0, 0.e0,
                l_mult_897 + (-0.6264e3 * l_prod_672 + (0.3132e4 * l_prod_671 + (-0.6696e4 * l_prod_670 + (0.6480e4 * l_prod_669 + (l_mult_1501 + (l_mult_829 + (l_mult_1250 + (-0.15066e5 * l_prod_665 + (0.20736e5 * l_prod_664 + (l_mult_1502 + (l_mult_901 + (l_mult_833 + (l_mult_1006 + (l_mult_1503 + (l_mult_902 + (l_mult_876 + (l_mult_836 + (l_mult_903 + (l_mult_1504 + (l_mult_895 + (l_mult_1255 + (l_mult_904 + (-0.15066e5 * l_prod_650 + (0.20736e5 * l_prod_649 + (l_mult_1505 + (l_mult_840 + (l_mult_907 + (l_mult_1066 + (l_mult_1506 + (l_mult_1310 + (l_mult_843 + (l_mult_1259 + (l_mult_1311 + (l_mult_909 + (l_mult_1309 + (l_mult_1262 + (l_mult_1312 + (l_mult_911 + (l_mult_1507 + (l_mult_847 + (l_mult_949 + (l_mult_913 + (0.7776e4 * l_prod_620 + (l_mult_850 + (l_mult_885 + (l_mult_1265 + (l_mult_1331 + (l_mult_1315 + (l_mult_852 + (l_mult_1268 + (l_mult_1299 + (l_mult_1270 + (l_mult_1508 + (l_mult_854 + l_mult_1272)))))))))))))))))))))))))))))))))))))))))))))))))))))))) + (vdot4(
            v_575,
            vcons4(
                l_mult_922 + (0.1053e4 * l_prod_672 + (-0.62235e4 * l_prod_671 + (0.14796e5 * l_prod_670 + (-0.15390e5 * l_prod_669 + (l_mult_1509 + (l_mult_856 + (l_mult_1273 + (0.23652e5 * l_prod_665 + (-0.37584e5 * l_prod_664 + (l_mult_1510 + (l_mult_925 + (l_mult_859 + (l_mult_1274 + (l_mult_992 + (l_mult_927 + (-0.7128e4 * l_prod_657 + (l_mult_860 + (l_mult_740 + (l_mult_678 + (l_mult_1277 + (l_mult_928 + (0.23652e5 * l_prod_650 + (-0.37584e5 * l_prod_649 + (l_mult_1511 + (l_mult_863 + (l_mult_930 + (-0.58320e5 * l_prod_643 + (l_mult_1050 + (l_mult_1321 + (l_mult_866 + (l_mult_1280 + (l_mult_1247 + (l_mult_683 + (l_mult_1282 + (l_mult_1322 + (l_mult_931 + (l_mult_1420 + (l_mult_869 + (l_mult_1285 + (l_mult_932 + (l_mult_1512 + (l_mult_687 + (l_mult_1287 + (-0.7128e4 * l_prod_613 + (l_mult_1323 + (l_mult_728 + (l_mult_690 + l_sum_1499))))))))))))))))))))))))))))))))))))))))))))))),
                l_mult_941 + (-0.1016e4 * l_prod_672 + (0.6696e4 * l_prod_671 + (-0.17424e5 * l_prod_670 + (0.19440e5 * l_prod_669 + (-0.7776e4 * l_prod_668 + (l_mult_873 + (l_mult_1289 + (-0.18360e5 * l_prod_665 + (0.33696e5 * l_prod_664 + (l_mult_1513 + (l_mult_944 + (l_mult_875 + (l_mult_753 + (l_mult_1503 + (l_mult_945 + (0.1296e4 * l_prod_657 + (l_mult_877 + (l_mult_1291 + (l_mult_946 + (-0.18360e5 * l_prod_650 + (0.33696e5 * l_prod_649 + (l_mult_1514 + (l_mult_879 + (l_mult_948 + (l_mult_763 + (l_mult_1506 + (l_mult_1329 + (l_mult_881 + (l_mult_1293 + (l_mult_1295 + (l_mult_1330 + (l_mult_769 + (l_mult_1507 + (l_mult_883 + (l_mult_1297 + (l_mult_848 + (l_mult_1298 + (0.1296e4 * l_prod_613 + l_mult_1266)))))))))))))))))))))))))))))))))))))),
                l_mult_952 + (0.594e3 * l_prod_672 + (-0.41445e4 * l_prod_671 + (0.11556e5 * l_prod_670 + (-0.13770e5 * l_prod_669 + (l_mult_1509 + (l_mult_888 + (l_mult_1300 + (0.7128e4 * l_prod_665 + (-0.14904e5 * l_prod_664 + (l_mult_1515 + (l_mult_724 + (l_mult_889 + (l_mult_1231 + (l_mult_676 + (l_mult_1302 + (l_mult_953 + (0.7128e4 * l_prod_650 + (-0.14904e5 * l_prod_649 + (l_mult_1516 + (l_mult_741 + (l_mult_955 + (l_mult_1517 + (l_mult_681 + l_sum_1500))))))))))))))))))))))),
                l_mult_959 + (-0.1944e3 * l_prod_672 + (0.1404e4 * l_prod_671 + (-0.4104e4 * l_prod_670 + (l_mult_746 + (l_mult_1501 + (l_mult_893 + (l_mult_1305 + (-0.1134e4 * l_prod_665 + (l_mult_1518 + (l_mult_832 + (l_mult_1307 + (l_mult_960 + (-0.1134e4 * l_prod_650 + (l_mult_1519 + l_mult_1258)))))))))))))))) + (vdot4(
            v_576, vcons4(l_sum_855, l_sum_872, l_sum_887, l_sum_892)) + (vdot4(v_577,
            vcons4(l_sum_896, l_sum_1320, l_sum_1328, l_sum_1333)) + (vdot4(v_578,
            vcons4(l_sum_1337, l_sum_1339, l_mult_741 + (l_mult_1238 + (l_mult_1240 + l_mult_693)),
                l_mult_1520 + (l_mult_1244 + (l_mult_1521 + (l_mult_736 + l_sum_729))))) + (vdot4(v_579,
            vcons4(l_mult_1520 + (l_mult_1522 + (l_mult_1247 + l_sum_737)), l_sum_743,
                l_mult_1520 + (l_mult_970 + (l_mult_1521 + (l_mult_1342 + (l_mult_728 + l_mult_690)))),
                l_mult_1523 + (l_mult_978 + (l_mult_1524 + (l_mult_794 + l_sum_1525))))) + (vdot4(v_580,
            vcons4(l_mult_1520 + (l_mult_970 + (l_mult_1522 + (l_mult_968 + (l_mult_1247 + l_mult_683)))),
                l_mult_1520 + (l_mult_967 + (l_mult_864 + l_sum_1526)),
                l_mult_1520 + (l_mult_967 + (l_mult_864 + (l_mult_1244 + (l_mult_968 + l_mult_682)))),
                l_mult_741 + (l_mult_955 + (l_mult_1517 + l_mult_681)))) + (vdot4(v_581,
            vcons4(
                l_mult_1131 + (l_mult_1349 + (0.19278e5 * l_prod_665 + (l_mult_1227 + (l_mult_1515 + (l_mult_1132 + (0.25704e5 * l_prod_661 + (l_mult_1351 + (l_mult_1527 + (l_mult_1133 + (-0.34992e5 * l_prod_657 + (l_mult_1528 + (l_mult_1134 + (l_mult_1529 + (l_mult_891 + (l_mult_1530 + (l_mult_1048 + (l_mult_1531 + (l_mult_1532 + (l_mult_1533 + (l_mult_1136 + (l_mult_1359 + (l_mult_1534 + (l_mult_1137 + (l_mult_1535 + (0.6426e4 * l_prod_624 + (l_mult_1362 + (l_mult_932 + (l_mult_1536 + (l_mult_1537 + (l_mult_871 + (-0.5832e4 * l_prod_609 + (l_mult_1365 + (l_mult_1538 + l_mult_693))))))))))))))))))))))))))))))))),
                l_mult_1144 + (l_mult_1391 + (-0.4860e4 * l_prod_665 + (l_mult_1518 + (l_mult_1168 + (l_mult_1539 + (l_mult_1392 + (l_mult_1503 + (l_mult_1169 + (0.42768e5 * l_prod_657 + (l_mult_1540 + (l_mult_1170 + (-0.23328e5 * l_prod_654 + (l_mult_878 + (l_mult_1541 + (l_mult_1095 + (l_mult_1542 + (l_mult_1543 + (l_mult_1150 + (l_mult_1259 + (l_mult_1544 + (l_mult_1172 + (l_mult_1545 + (l_mult_1546 + (l_mult_1297 + (l_mult_1547 + (l_mult_850 + (l_mult_851 + l_sum_1549))))))))))))))))))))))))))),
                l_mult_1154 + (l_mult_1409 + (l_mult_1233 + (l_mult_1190 + (l_mult_1550 + (l_mult_1411 + (l_mult_1191 + (-0.22032e5 * l_prod_657 + (l_mult_860 + (l_mult_1192 + (l_mult_1529 + (l_mult_862 + (l_mult_1551 + (l_mult_1120 + (l_mult_1552 + (l_mult_1159 + (l_mult_1553 + (l_mult_1051 + (l_mult_1535 + l_sum_737)))))))))))))))))),
                l_mult_1162 + (l_mult_1417 + (l_mult_1199 + (l_mult_1554 + (l_mult_1200 + (l_mult_1555 + (l_mult_1201 + (l_mult_1504 + (l_mult_839 + l_sum_1558)))))))))) + (vdot4(
            v_582,
            vcons4(
                l_mult_1144 + (l_mult_1366 + (-0.28836e5 * l_prod_665 + (0.41472e5 * l_prod_664 + (l_mult_1513 + (l_mult_1146 + (l_mult_1539 + (l_mult_1368 + (-0.46656e5 * l_prod_659 + (l_mult_1147 + (0.23328e5 * l_prod_657 + (l_mult_1540 + (l_mult_894 + (l_mult_838 + (l_mult_1541 + (l_mult_1065 + (l_mult_1559 + (l_mult_1096 + (l_mult_1560 + (l_mult_1150 + (l_mult_1374 + (l_mult_1557 + (l_mult_845 + (l_mult_1546 + (l_mult_1376 + (l_mult_913 + (l_mult_884 + (l_mult_850 + (l_mult_1548 + l_mult_776)))))))))))))))))))))))))))),
                l_mult_1175 + (l_mult_1398 + (0.5832e4 * l_prod_665 + (l_mult_1380 + (l_mult_1176 + (0.17496e5 * l_prod_661 + (l_mult_1399 + (l_mult_992 + (l_mult_1177 + (-0.27216e5 * l_prod_657 + (l_mult_1528 + (l_mult_1178 + (l_mult_861 + (l_mult_1561 + (l_mult_1104 + (l_mult_1517 + (l_mult_1562 + (l_mult_1180 + (l_mult_1280 + (l_mult_742 + (l_mult_867 + l_sum_1525)))))))))))))))))))),
                l_mult_1183 + (l_mult_1414 + (-0.648e3 * l_prod_665 + (l_mult_1196 + (l_mult_1563 + (l_mult_1388 + (l_mult_1197 + (0.12960e5 * l_prod_657 + (l_mult_836 + (l_mult_894 + (l_mult_838 + (l_mult_1564 + (l_mult_811 + (l_mult_1565 + (l_mult_881 + (l_mult_1566 + l_mult_767))))))))))))))),
                l_mult_1154 + (l_mult_1378 + (0.21060e5 * l_prod_665 + (-0.36288e5 * l_prod_664 + (l_mult_1510 + (l_mult_1155 + (l_mult_1550 + (l_mult_1381 + (l_mult_1527 + (l_mult_1156 + (l_mult_1400 + (l_mult_860 + (l_mult_1551 + (l_mult_1078 + (l_mult_1567 + (l_mult_1532 + (l_mult_1568 + (l_mult_1159 + (l_mult_1385 + l_sum_1526)))))))))))))))))))) + (vdot4(
            v_583,
            vcons4(
                l_mult_1183 + (l_mult_1403 + (-0.2916e4 * l_prod_665 + (l_mult_1518 + (l_mult_1184 + (l_mult_1563 + (l_mult_1406 + (l_mult_1503 + (l_mult_1185 + (l_mult_1555 + (l_mult_836 + (l_mult_1564 + (l_mult_1114 + (l_mult_1569 + (l_mult_1329 + (l_mult_881 + l_mult_1293))))))))))))))),
                l_mult_1162 + (l_mult_1386 + (-0.7614e4 * l_prod_665 + (l_mult_750 + (l_mult_1502 + (l_mult_1163 + (l_mult_1554 + (l_mult_1388 + (l_mult_834 + (l_mult_1556 + (l_mult_1086 + (l_mult_1542 + l_mult_813))))))))))),
                l_mult_1440 + (l_mult_1046 + (0.19278e5 * l_prod_650 + (l_mult_1497 + (l_mult_1516 + (l_mult_1530 + (l_mult_1048 + (l_mult_1531 + (l_mult_1532 + (0.6426e4 * l_prod_640 + (l_mult_998 + (l_mult_1280 + (-0.5832e4 * l_prod_636 + (l_mult_1051 + (l_mult_684 + (l_mult_1442 + (0.25704e5 * l_prod_629 + (l_mult_1052 + (l_mult_1570 + (l_mult_1571 + (l_mult_1212 + (l_mult_1055 + (l_mult_1536 + (l_mult_1537 + (l_mult_1572 + (l_mult_1444 + (-0.34992e5 * l_prod_613 + (l_mult_1573 + (l_mult_1574 + (l_mult_936 + (l_mult_1288 + (l_mult_1446 + (l_mult_1575 + (l_mult_1576 + l_mult_1336))))))))))))))))))))))))))))))))),
                l_mult_1448 + (l_mult_1094 + (-0.4860e4 * l_prod_650 + (l_mult_1519 + (l_mult_1541 + (l_mult_1095 + (l_mult_1542 + (l_mult_1577 + (l_mult_881 + (l_mult_1566 + (l_mult_1450 + (l_mult_1578 + (l_mult_1097 + (l_mult_1507 + (l_mult_1579 + (l_mult_912 + (l_mult_913 + (l_mult_1547 + (l_mult_850 + (l_mult_885 + (l_mult_1452 + (0.42768e5 * l_prod_613 + (l_mult_1580 + (l_mult_1581 + (l_mult_1071 + (l_mult_1269 + (l_mult_1454 + (-0.23328e5 * l_prod_602 + (l_mult_1582 + l_mult_1332)))))))))))))))))))))))))))))) + (vdot4(
            v_584,
            vcons4(
                l_mult_1455 + (l_mult_1119 + (l_mult_710 + (l_mult_1551 + (l_mult_1120 + (l_mult_1244 + (l_mult_1456 + (l_mult_1583 + (l_mult_1123 + (l_mult_1584 + (l_mult_1458 + (l_mult_736 + (l_mult_1459 + (-0.22032e5 * l_prod_613 + (l_mult_1323 + (l_mult_1585 + (l_mult_1365 + (l_mult_691 + (l_mult_1460 + (l_mult_1575 + (l_mult_1576 + l_mult_1327)))))))))))))))))))),
                l_mult_1461 + (l_mult_1130 + (l_mult_1556 + (l_mult_1462 + (l_mult_1586 + (l_mult_1587 + (l_mult_1464 + (l_mult_1588 + (l_mult_1589 + (l_mult_1466 + (l_mult_1508 + (l_mult_854 + l_mult_1319))))))))))),
                l_mult_1448 + (l_mult_1063 + (-0.28836e5 * l_prod_650 + (0.41472e5 * l_prod_649 + (l_mult_1514 + (l_mult_1541 + (l_mult_1065 + (l_mult_1559 + (l_mult_1096 + (l_mult_1577 + (l_mult_1011 + (l_mult_1259 + (l_mult_1566 + (l_mult_767 + (l_mult_1467 + (l_mult_1578 + (l_mult_1067 + (-0.46656e5 * l_prod_627 + (l_mult_1590 + (l_mult_912 + (l_mult_1070 + (l_mult_884 + (l_mult_850 + (l_mult_1468 + (0.23328e5 * l_prod_613 + (l_mult_1580 + (l_mult_1589 + (l_mult_1316 + l_sum_1591))))))))))))))))))))))))))),
                l_mult_1470 + (l_mult_1103 + (0.5832e4 * l_prod_650 + (l_mult_1077 + (l_mult_1561 + (l_mult_1104 + (l_mult_1517 + (l_mult_1524 + (l_mult_794 + (l_mult_1471 + (0.17496e5 * l_prod_629 + (l_mult_1106 + (l_mult_1420 + (l_mult_1592 + (l_mult_1472 + (l_mult_932 + (l_mult_1512 + (l_mult_687 + (l_mult_1473 + (-0.27216e5 * l_prod_613 + (l_mult_1573 + (l_mult_1240 + (l_mult_1324 + (l_mult_1474 + l_mult_1325))))))))))))))))))))))))) + (vdot4(
            v_585,
            vcons4(
                l_mult_1475 + (l_mult_1126 + (-0.648e3 * l_prod_650 + (l_mult_1564 + (l_mult_811 + (l_mult_1476 + (l_mult_1593 + (l_mult_1087 + (l_mult_1594 + (l_mult_1297 + (l_mult_1478 + (0.12960e5 * l_prod_613 + (l_mult_1315 + (l_mult_1548 + (l_mult_776 + l_sum_1591)))))))))))))),
                l_mult_1455 + (l_mult_1075 + (0.21060e5 * l_prod_650 + (-0.36288e5 * l_prod_649 + (l_mult_1511 + (l_mult_1551 + (l_mult_1078 + (l_mult_1567 + (l_mult_1532 + (l_mult_1244 + (l_mult_968 + (l_mult_682 + (l_mult_1479 + (l_mult_1583 + (l_mult_1080 + (l_mult_1570 + (l_mult_1595 + (l_mult_1458 + (l_mult_1000 + (l_mult_1480 + (l_mult_1107 + l_mult_1323)))))))))))))))))))),
                l_mult_1475 + (l_mult_1111 + (-0.2916e4 * l_prod_650 + (l_mult_1519 + (l_mult_1564 + (l_mult_1114 + (l_mult_1569 + (l_mult_1481 + (l_mult_1593 + (l_mult_1115 + (l_mult_1507 + (l_mult_883 + (l_mult_1297 + (l_mult_848 + (l_mult_1482 + (l_mult_1588 + l_mult_1315))))))))))))))),
                l_mult_1461 + (l_mult_1084 + (-0.7614e4 * l_prod_650 + (l_mult_760 + (l_mult_1505 + (l_mult_1556 + (l_mult_1086 + (l_mult_1542 + (l_mult_813 + (l_mult_1483 + (l_mult_1586 + (l_mult_1087 + l_mult_1313))))))))))))) + (vdot4(
            v_586,
            vcons4(
                l_mult_1530 + (l_mult_996 + (-0.17496e5 * l_prod_643 + (l_mult_681 + (l_mult_1533 + (l_mult_998 + (l_mult_1385 + (l_mult_1534 + (l_mult_867 + (l_mult_1535 + (l_mult_1571 + (l_mult_1362 + (l_mult_1000 + (-0.34992e5 * l_prod_620 + (l_mult_1537 + (l_mult_1286 + (l_mult_1574 + (l_mult_1324 + (l_mult_937 + l_mult_1576)))))))))))))))))),
                l_mult_1541 + (l_mult_1010 + (l_mult_1569 + (l_mult_1560 + (l_mult_881 + (l_mult_1557 + (l_mult_1579 + (l_mult_1376 + (l_mult_848 + (l_mult_914 + (l_mult_850 + (l_mult_851 + (l_mult_1581 + (l_mult_1316 + (l_mult_1596 + l_mult_1582)))))))))))))),
                l_mult_1551 + (l_mult_970 + (l_mult_1568 + (l_mult_1584 + (l_mult_1342 + (l_mult_1597 + (l_mult_1585 + (l_mult_690 + (l_mult_1538 + l_mult_1576)))))))),
                l_mult_1556 + (l_mult_1587 + (l_mult_1589 + l_mult_854)))) + (vdot4(v_587,
            vcons4(
                l_mult_1541 + (l_mult_1010 + (l_mult_1569 + (l_mult_1543 + (l_mult_1011 + (l_mult_1293 + (l_mult_1544 + (l_mult_845 + (l_mult_1545 + (l_mult_1590 + (l_mult_1297 + (l_mult_914 + (l_mult_850 + (l_mult_1598 + l_sum_1599))))))))))))),
                l_mult_1561 + (l_mult_978 + (l_mult_1562 + (l_mult_794 + (l_mult_742 + (l_mult_1592 + (l_mult_798 + (-0.14580e5 * l_prod_620 + (l_mult_687 + (l_mult_871 + l_sum_1600))))))))),
                l_mult_1564 + (l_mult_1329 + (l_mult_1594 + (l_mult_884 + l_sum_1549))),
                l_mult_1551 + (l_mult_970 + (l_mult_1552 + (l_mult_968 + (l_mult_1553 + (l_mult_683 + (l_mult_1535 + (l_mult_1595 + (l_mult_1597 + l_mult_1572)))))))))) + (vdot4(
            v_588,
            vcons4(l_mult_1564 + (l_mult_1565 + (l_mult_1566 + l_sum_886)), l_sum_1558,
                0.4320e4 * l_prod_645 + (l_mult_1203 + (0.58320e5 * l_prod_643 + (l_mult_1506 + (-0.15984e5 * l_prod_640 + (l_mult_1205 + (l_mult_1374 + (0.19440e5 * l_prod_636 + (l_mult_1172 + (l_mult_846 + (-0.15984e5 * l_prod_624 + (l_mult_1484 + (l_mult_1070 + (0.38880e5 * l_prod_620 + (-0.93312e5 * l_prod_619 + (l_mult_1598 + (0.19440e5 * l_prod_609 + (l_mult_1071 + (l_mult_1596 + l_mult_1318)))))))))))))))))),
                l_mult_1601 + (l_mult_1221 + (l_mult_1517 + (l_mult_1602 + (l_mult_1159 + (l_mult_742 + (0.13284e5 * l_prod_624 + (l_mult_1488 + (l_mult_932 + (l_mult_1603 + (l_mult_1537 + (l_mult_871 + (l_mult_1142 + (l_mult_936 + (l_mult_937 + l_mult_1326)))))))))))))))) + (vdot4(
            v_589,
            vcons4(
                l_mult_1604 + (l_mult_1226 + (l_mult_1605 + (-0.4320e4 * l_prod_624 + (l_mult_1492 + (l_mult_1606 + (0.11664e5 * l_prod_609 + (l_mult_1268 + (l_mult_853 + l_mult_1318)))))))),
                l_mult_1601 + (l_mult_1221 + (l_mult_1517 + (0.13284e5 * l_prod_640 + (l_mult_1211 + (l_mult_1280 + (l_mult_1441 + (l_mult_1137 + (l_mult_868 + (l_mult_1607 + (l_mult_1458 + (l_mult_1603 + (l_mult_1537 + (l_mult_1286 + l_sum_1600))))))))))))),
                l_mult_1608 + (l_mult_1114 + (l_mult_1609 + (l_mult_881 + (l_mult_1557 + (l_mult_1610 + (l_mult_1297 + (l_mult_849 + (l_mult_850 + (l_mult_851 + l_sum_1599))))))))),
                l_mult_1604 + (l_mult_1226 + (-0.4320e4 * l_prod_640 + (l_mult_1218 + (0.11664e5 * l_prod_636 + (l_mult_909 + (l_mult_846 + (l_mult_1611 + (l_mult_1606 + l_mult_1314)))))))))) + (vdot4(
            v_590,
            vcons4(
                l_mult_1601 + (l_mult_1209 + (l_mult_1049 + (l_mult_1050 + (l_mult_1602 + (l_mult_1211 + (l_mult_1359 + (l_mult_742 + (l_mult_867 + (l_mult_1607 + (l_mult_1488 + (l_mult_1055 + (l_mult_1597 + (l_mult_1537 + (l_mult_1240 + l_mult_1324)))))))))))))),
                l_mult_1608 + (l_mult_1224 + (l_mult_1542 + (l_mult_1565 + (l_mult_881 + (l_mult_1610 + (l_mult_949 + (l_mult_913 + (l_mult_884 + (l_mult_850 + (l_mult_1589 + l_mult_1316)))))))))),
                l_mult_1608 + (l_mult_1224 + (l_mult_1542 + (l_mult_1609 + (l_mult_843 + (l_mult_1259 + (l_mult_1557 + (l_mult_845 + (l_mult_1594 + (l_mult_1297 + (l_mult_884 + l_mult_850)))))))))),
                l_mult_1604 + (l_mult_1216 + (0.34992e5 * l_prod_643 + (l_mult_1506 + (l_mult_1605 + (l_mult_1218 + (l_mult_766 + (l_mult_1611 + (l_mult_1492 + l_mult_772)))))))))) + vdot4(
            v_570,
            vcons4(l_sum_695, 0.e0, 0.e0,
                l_mult_696 + (0.274e2 * l_prod_672 + (-0.2025e3 * l_prod_671 + (0.612e3 * l_prod_670 + (-0.810e3 * l_prod_669 + l_mult_674)))))))))))))))))))))))));
        int32_t l_mulRes_1615 = l__t_399 * 20;
        int32_t t_1616 = l__t_400.indexMap[l_mulRes_1615];
        int32_t l_mulRes_1617 = 3 * t_1616;
        int32_t t_1618 = l__t_400.indexMap[l_mulRes_1615 + 1];
        int32_t l_mulRes_1619 = 3 * t_1618;
        double l_dof_load_1620 = l__t_400.coordMap[l_mulRes_1619];
        double l_dof_load_1621 = l__t_400.coordMap[1 + l_mulRes_1619];
        double l_dof_load_1622 = l__t_400.coordMap[2 + l_mulRes_1619];
        int32_t t_1623 = l__t_400.indexMap[l_mulRes_1615 + 2];
        int32_t l_mulRes_1624 = 3 * t_1623;
        double l_dof_load_1625 = l__t_400.coordMap[l_mulRes_1624];
        double l_dof_load_1626 = l__t_400.coordMap[1 + l_mulRes_1624];
        double l_dof_load_1627 = l__t_400.coordMap[2 + l_mulRes_1624];
        int32_t t_1628 = l__t_400.indexMap[l_mulRes_1615 + 3];
        int32_t l_mulRes_1629 = 3 * t_1628;
        double l_dof_load_1630 = l__t_400.coordMap[l_mulRes_1629];
        double l_dof_load_1631 = l__t_400.coordMap[1 + l_mulRes_1629];
        double l_dof_load_1632 = l__t_400.coordMap[2 + l_mulRes_1629];
        int32_t t_1633 = l__t_400.indexMap[l_mulRes_1615 + 4];
        int32_t l_mulRes_1634 = 3 * t_1633;
        double l_dof_load_1635 = l__t_400.coordMap[l_mulRes_1634];
        double l_dof_load_1636 = l__t_400.coordMap[1 + l_mulRes_1634];
        double l_dof_load_1637 = l__t_400.coordMap[2 + l_mulRes_1634];
        int32_t t_1638 = l__t_400.indexMap[l_mulRes_1615 + 5];
        int32_t l_mulRes_1639 = 3 * t_1638;
        double l_dof_load_1640 = l__t_400.coordMap[l_mulRes_1639];
        double l_dof_load_1641 = l__t_400.coordMap[1 + l_mulRes_1639];
        double l_dof_load_1642 = l__t_400.coordMap[2 + l_mulRes_1639];
        int32_t t_1643 = l__t_400.indexMap[l_mulRes_1615 + 6];
        int32_t l_mulRes_1644 = 3 * t_1643;
        double l_dof_load_1645 = l__t_400.coordMap[l_mulRes_1644];
        double l_dof_load_1646 = l__t_400.coordMap[1 + l_mulRes_1644];
        double l_dof_load_1647 = l__t_400.coordMap[2 + l_mulRes_1644];
        int32_t t_1648 = l__t_400.indexMap[l_mulRes_1615 + 7];
        int32_t l_mulRes_1649 = 3 * t_1648;
        double l_dof_load_1650 = l__t_400.coordMap[l_mulRes_1649];
        double l_dof_load_1651 = l__t_400.coordMap[1 + l_mulRes_1649];
        double l_dof_load_1652 = l__t_400.coordMap[2 + l_mulRes_1649];
        int32_t t_1653 = l__t_400.indexMap[l_mulRes_1615 + 8];
        int32_t l_mulRes_1654 = 3 * t_1653;
        double l_dof_load_1655 = l__t_400.coordMap[l_mulRes_1654];
        double l_dof_load_1656 = l__t_400.coordMap[1 + l_mulRes_1654];
        double l_dof_load_1657 = l__t_400.coordMap[2 + l_mulRes_1654];
        int32_t t_1658 = l__t_400.indexMap[l_mulRes_1615 + 9];
        int32_t l_mulRes_1659 = 3 * t_1658;
        double l_dof_load_1660 = l__t_400.coordMap[l_mulRes_1659];
        double l_dof_load_1661 = l__t_400.coordMap[1 + l_mulRes_1659];
        double l_dof_load_1662 = l__t_400.coordMap[2 + l_mulRes_1659];
        int32_t t_1663 = l__t_400.indexMap[l_mulRes_1615 + 10];
        int32_t l_mulRes_1664 = 3 * t_1663;
        double l_dof_load_1665 = l__t_400.coordMap[l_mulRes_1664];
        double l_dof_load_1666 = l__t_400.coordMap[1 + l_mulRes_1664];
        double l_dof_load_1667 = l__t_400.coordMap[2 + l_mulRes_1664];
        int32_t t_1668 = l__t_400.indexMap[l_mulRes_1615 + 11];
        int32_t l_mulRes_1669 = 3 * t_1668;
        double l_dof_load_1670 = l__t_400.coordMap[l_mulRes_1669];
        double l_dof_load_1671 = l__t_400.coordMap[1 + l_mulRes_1669];
        double l_dof_load_1672 = l__t_400.coordMap[2 + l_mulRes_1669];
        int32_t t_1673 = l__t_400.indexMap[l_mulRes_1615 + 12];
        int32_t l_mulRes_1674 = 3 * t_1673;
        double l_dof_load_1675 = l__t_400.coordMap[l_mulRes_1674];
        double l_dof_load_1676 = l__t_400.coordMap[1 + l_mulRes_1674];
        double l_dof_load_1677 = l__t_400.coordMap[2 + l_mulRes_1674];
        int32_t t_1678 = l__t_400.indexMap[l_mulRes_1615 + 13];
        int32_t l_mulRes_1679 = 3 * t_1678;
        double l_dof_load_1680 = l__t_400.coordMap[l_mulRes_1679];
        double l_dof_load_1681 = l__t_400.coordMap[1 + l_mulRes_1679];
        double l_dof_load_1682 = l__t_400.coordMap[2 + l_mulRes_1679];
        int32_t t_1683 = l__t_400.indexMap[l_mulRes_1615 + 14];
        int32_t l_mulRes_1684 = 3 * t_1683;
        double l_dof_load_1685 = l__t_400.coordMap[l_mulRes_1684];
        double l_dof_load_1686 = l__t_400.coordMap[1 + l_mulRes_1684];
        double l_dof_load_1687 = l__t_400.coordMap[2 + l_mulRes_1684];
        int32_t t_1688 = l__t_400.indexMap[l_mulRes_1615 + 15];
        int32_t l_mulRes_1689 = 3 * t_1688;
        double l_dof_load_1690 = l__t_400.coordMap[l_mulRes_1689];
        double l_dof_load_1691 = l__t_400.coordMap[1 + l_mulRes_1689];
        double l_dof_load_1692 = l__t_400.coordMap[2 + l_mulRes_1689];
        int32_t t_1693 = l__t_400.indexMap[l_mulRes_1615 + 16];
        int32_t l_mulRes_1694 = 3 * t_1693;
        double l_dof_load_1695 = l__t_400.coordMap[l_mulRes_1694];
        double l_dof_load_1696 = l__t_400.coordMap[1 + l_mulRes_1694];
        double l_dof_load_1697 = l__t_400.coordMap[2 + l_mulRes_1694];
        int32_t t_1698 = l__t_400.indexMap[l_mulRes_1615 + 17];
        int32_t l_mulRes_1699 = 3 * t_1698;
        double l_dof_load_1700 = l__t_400.coordMap[l_mulRes_1699];
        double l_dof_load_1701 = l__t_400.coordMap[1 + l_mulRes_1699];
        double l_dof_load_1702 = l__t_400.coordMap[2 + l_mulRes_1699];
        int32_t t_1703 = l__t_400.indexMap[l_mulRes_1615 + 18];
        int32_t l_mulRes_1704 = 3 * t_1703;
        double l_dof_load_1705 = l__t_400.coordMap[l_mulRes_1704];
        double l_dof_load_1706 = l__t_400.coordMap[1 + l_mulRes_1704];
        double l_dof_load_1707 = l__t_400.coordMap[2 + l_mulRes_1704];
        int32_t t_1708 = l__t_400.indexMap[l_mulRes_1615 + 19];
        int32_t l_mulRes_1709 = 3 * t_1708;
        double l_dof_load_1710 = l__t_400.coordMap[l_mulRes_1709];
        double l_dof_load_1711 = l__t_400.coordMap[1 + l_mulRes_1709];
        double l_dof_load_1712 = l__t_400.coordMap[2 + l_mulRes_1709];
        double l_mult_1713 = -0.135e2 * l_prod_671;
        double l_mult_1714 = -0.135e2 * l_prod_662;
        double l_mult_1715 = -0.135e2 * l_prod_630;
        double l_sum_1716 = -0.55e1 * l_prod_673 + (0.18e2 * l_prod_672 + (l_mult_1713 + (0.18e2 * l_prod_667 + (l_mult_975 + (l_mult_1714 + (0.18e2 * l_prod_652 + (l_mult_1345 + (l_mult_1523 + l_mult_1715))))))));
        double l_mult_1717 = 0.1e1 * l_prod_673;
        double l_mult_1718 = 0.135e2 * l_prod_630;
        double l_sum_1719 = l_mult_1717 + (-0.9e1 * l_prod_652 + l_mult_1718);
        double l_mult_1720 = -0.45e1 * l_prod_672;
        double l_mult_1721 = 0.27e2 * l_prod_651;
        double l_sum_1722 = l_mult_1720 + l_mult_1721;
        double l_mult_1723 = 0.135e2 * l_prod_671;
        double l_sum_1724 = l_mult_1720 + l_mult_1723;
        double l_mult_1725 = -0.45e1 * l_prod_667;
        double l_mult_1726 = 0.27e2 * l_prod_645;
        double l_sum_1727 = l_mult_1725 + l_mult_1726;
        double l_mult_1728 = 0.135e2 * l_prod_662;
        double l_sum_1729 = l_mult_1725 + l_mult_1728;
        double l_mult_1730 = -0.225e2 * l_prod_672;
        double l_mult_1731 = 0.27e2 * l_prod_666;
        double l_sum_1732 = l_mult_1730 + (0.27e2 * l_prod_671 + (l_mult_1731 + l_mult_1721));
        double l_sum_1733 = l_mult_699 + l_mult_1713;
        double l_mult_1734 = -0.225e2 * l_prod_667;
        double l_sum_1735 = l_mult_1734 + (l_mult_1731 + (0.27e2 * l_prod_662 + l_mult_1726));
        double l_sum_1736 = l_mult_723 + l_mult_1714;
        double l_mult_1737 = 0.9e1 * l_prod_673;
        double l_sum_1738 = l_mult_1737 + (l_mult_1730 + (l_mult_1723 + (l_mult_1734 + (l_mult_1731 + (l_mult_1728 + (-0.45e2 * l_prod_652 + (l_mult_1389 + (l_mult_1556 + 0.405e2 * l_prod_630))))))));
        double l_mult_1739 = -0.45e1 * l_prod_673;
        double l_sum_1740 = l_mult_1739 + (l_mult_699 + (l_mult_723 + (0.36e2 * l_prod_652 + (l_mult_1345 + (l_mult_1523 + -0.405e2 * l_prod_630)))));
        double l_mult_1741 = 0.27e2 * l_prod_672;
        double l_sum_1742 = l_mult_1741 + (l_mult_700 + (l_mult_975 + l_mult_717));
        double l_mult_1743 = 0.27e2 * l_prod_667;
        double l_sum_1744 = l_mult_1743 + (l_mult_975 + (l_mult_724 + l_mult_741));
        double l_sum_1745 = l_mult_1717 + (-0.9e1 * l_prod_667 + l_mult_1728);
        double l_sum_1746 = l_mult_1720 + l_mult_1731;
        double l_mult_1747 = -0.45e1 * l_prod_652;
        double l_sum_1748 = l_mult_1747 + l_mult_1718;
        double l_sum_1749 = l_mult_1747 + l_mult_1726;
        double l_mult_1750 = -0.225e2 * l_prod_652;
        double l_sum_1751 = l_mult_1737 + (l_mult_1730 + (l_mult_1723 + (-0.45e2 * l_prod_667 + (l_mult_1019 + (0.405e2 * l_prod_662 + (l_mult_1750 + (l_mult_1721 + (l_mult_1556 + l_mult_1718))))))));
        double l_sum_1752 = l_mult_1739 + (l_mult_699 + (0.36e2 * l_prod_667 + (l_mult_975 + (-0.405e2 * l_prod_662 + (l_mult_1236 + l_mult_1523)))));
        double l_sum_1753 = l_mult_1750 + (l_mult_1721 + (l_mult_1726 + 0.27e2 * l_prod_630));
        double l_sum_1754 = l_mult_1236 + l_mult_1715;
        double l_sum_1755 = l_mult_1741 + (l_mult_700 + (l_mult_821 + l_mult_1345));
        double l_mult_1756 = 0.27e2 * l_prod_652;
        double l_sum_1757 = l_mult_1756 + (l_mult_1345 + (l_mult_741 + l_mult_1248));
        double l_sum_1758 = l_mult_1717 + (-0.9e1 * l_prod_672 + l_mult_1723);
        double l_sum_1759 = l_mult_1725 + l_mult_1731;
        double l_sum_1760 = l_mult_1747 + l_mult_1721;
        double l_sum_1761 = l_mult_1737 + (-0.45e2 * l_prod_672 + (0.405e2 * l_prod_671 + (l_mult_1734 + (l_mult_1019 + (l_mult_1728 + (l_mult_1750 + (l_mult_1389 + (l_mult_1726 + l_mult_1718))))))));
        double l_sum_1762 = l_mult_1739 + (0.36e2 * l_prod_672 + (-0.405e2 * l_prod_671 + (l_mult_723 + (l_mult_975 + (l_mult_1236 + l_mult_1345)))));
        double l_sum_1763 = l_mult_1743 + (l_mult_821 + (l_mult_724 + l_mult_1523));
        double l_sum_1764 = l_mult_1756 + (l_mult_717 + (l_mult_1523 + l_mult_1248));
        double t_1765 = l__t_400.coordMap[l_mulRes_1617];
        double l_r_1766 = t_1765 * l_sum_1716;
        double l_r_1767 = l_dof_load_1625 * 0.e0;
        double l_r_1768 = l_dof_load_1630 * 0.e0;
        double l_r_1769 = l_dof_load_1665 * l_sum_1732;
        double l_r_1770 = l_dof_load_1670 * l_sum_1733;
        double l_r_1771 = l_dof_load_1675 * l_sum_1735;
        double l_r_1772 = l_dof_load_1680 * l_sum_1736;
        double l_r_1773 = l_r_1766 + l_dof_load_1620 * l_sum_1719 + l_r_1767 + l_r_1768 + l_dof_load_1635 * 0.e0 + l_dof_load_1640 * 0.e0 + l_dof_load_1645 * l_sum_1722 + l_dof_load_1650 * l_sum_1724 + l_dof_load_1655 * l_sum_1727 + l_dof_load_1660 * l_sum_1729 + l_r_1769 + l_r_1770 + l_r_1771 + l_r_1772 + l_dof_load_1685 * l_sum_1738 + l_dof_load_1690 * l_sum_1740 + l_dof_load_1695 * l_mult_1731 + l_dof_load_1700 * l_mult_975 + l_dof_load_1705 * l_sum_1742 + l_dof_load_1710 * l_sum_1744;
        double l_r_1774 = l_dof_load_1685 * l_sum_1753;
        double l_r_1775 = l_dof_load_1690 * l_sum_1754;
        double l_r_1776 = l_r_1766 + l_dof_load_1620 * 0.e0;
        double l_r_1777 = l_r_1776 + l_dof_load_1625 * l_sum_1745 + l_r_1768 + l_dof_load_1635 * l_sum_1746 + l_dof_load_1640 * l_sum_1724 + l_dof_load_1645 * 0.e0 + l_dof_load_1650 * 0.e0 + l_dof_load_1655 * l_sum_1748 + l_dof_load_1660 * l_sum_1749 + l_r_1769 + l_r_1770 + l_dof_load_1675 * l_sum_1751 + l_dof_load_1680 * l_sum_1752 + l_r_1774 + l_r_1775 + l_dof_load_1695 * l_mult_1721 + l_dof_load_1700 * l_sum_1755 + l_dof_load_1705 * l_mult_1345 + l_dof_load_1710 * l_sum_1757;
        double l_r_1778 = l_r_1776 + l_r_1767 + l_dof_load_1630 * l_sum_1758 + l_dof_load_1635 * l_sum_1729 + l_dof_load_1640 * l_sum_1759 + l_dof_load_1645 * l_sum_1748 + l_dof_load_1650 * l_sum_1760 + l_dof_load_1655 * 0.e0 + l_dof_load_1660 * 0.e0 + l_dof_load_1665 * l_sum_1761 + l_dof_load_1670 * l_sum_1762 + l_r_1771 + l_r_1772 + l_r_1774 + l_r_1775 + l_dof_load_1695 * l_mult_1726 + l_dof_load_1700 * l_sum_1763 + l_dof_load_1705 * l_sum_1764 + l_dof_load_1710 * l_mult_1523;
        double t_1779 = l__t_400.coordMap[1 + l_mulRes_1617];
        double l_r_1780 = t_1779 * l_sum_1716;
        double l_r_1781 = l_dof_load_1626 * 0.e0;
        double l_r_1782 = l_dof_load_1631 * 0.e0;
        double l_r_1783 = l_dof_load_1666 * l_sum_1732;
        double l_r_1784 = l_dof_load_1671 * l_sum_1733;
        double l_r_1785 = l_dof_load_1676 * l_sum_1735;
        double l_r_1786 = l_dof_load_1681 * l_sum_1736;
        double l_r_1787 = l_r_1780 + l_dof_load_1621 * l_sum_1719 + l_r_1781 + l_r_1782 + l_dof_load_1636 * 0.e0 + l_dof_load_1641 * 0.e0 + l_dof_load_1646 * l_sum_1722 + l_dof_load_1651 * l_sum_1724 + l_dof_load_1656 * l_sum_1727 + l_dof_load_1661 * l_sum_1729 + l_r_1783 + l_r_1784 + l_r_1785 + l_r_1786 + l_dof_load_1686 * l_sum_1738 + l_dof_load_1691 * l_sum_1740 + l_dof_load_1696 * l_mult_1731 + l_dof_load_1701 * l_mult_975 + l_dof_load_1706 * l_sum_1742 + l_dof_load_1711 * l_sum_1744;
        double l_r_1788 = l_dof_load_1686 * l_sum_1753;
        double l_r_1789 = l_dof_load_1691 * l_sum_1754;
        double l_r_1790 = l_r_1780 + l_dof_load_1621 * 0.e0;
        double l_r_1791 = l_r_1790 + l_dof_load_1626 * l_sum_1745 + l_r_1782 + l_dof_load_1636 * l_sum_1746 + l_dof_load_1641 * l_sum_1724 + l_dof_load_1646 * 0.e0 + l_dof_load_1651 * 0.e0 + l_dof_load_1656 * l_sum_1748 + l_dof_load_1661 * l_sum_1749 + l_r_1783 + l_r_1784 + l_dof_load_1676 * l_sum_1751 + l_dof_load_1681 * l_sum_1752 + l_r_1788 + l_r_1789 + l_dof_load_1696 * l_mult_1721 + l_dof_load_1701 * l_sum_1755 + l_dof_load_1706 * l_mult_1345 + l_dof_load_1711 * l_sum_1757;
        double l_r_1792 = l_r_1790 + l_r_1781 + l_dof_load_1631 * l_sum_1758 + l_dof_load_1636 * l_sum_1729 + l_dof_load_1641 * l_sum_1759 + l_dof_load_1646 * l_sum_1748 + l_dof_load_1651 * l_sum_1760 + l_dof_load_1656 * 0.e0 + l_dof_load_1661 * 0.e0 + l_dof_load_1666 * l_sum_1761 + l_dof_load_1671 * l_sum_1762 + l_r_1785 + l_r_1786 + l_r_1788 + l_r_1789 + l_dof_load_1696 * l_mult_1726 + l_dof_load_1701 * l_sum_1763 + l_dof_load_1706 * l_sum_1764 + l_dof_load_1711 * l_mult_1523;
        double t_1793 = l__t_400.coordMap[2 + l_mulRes_1617];
        double l_r_1794 = t_1793 * l_sum_1716;
        double l_r_1795 = l_dof_load_1627 * 0.e0;
        double l_r_1796 = l_dof_load_1632 * 0.e0;
        double l_r_1797 = l_dof_load_1667 * l_sum_1732;
        double l_r_1798 = l_dof_load_1672 * l_sum_1733;
        double l_r_1799 = l_dof_load_1677 * l_sum_1735;
        double l_r_1800 = l_dof_load_1682 * l_sum_1736;
        double l_r_1801 = l_r_1794 + l_dof_load_1622 * l_sum_1719 + l_r_1795 + l_r_1796 + l_dof_load_1637 * 0.e0 + l_dof_load_1642 * 0.e0 + l_dof_load_1647 * l_sum_1722 + l_dof_load_1652 * l_sum_1724 + l_dof_load_1657 * l_sum_1727 + l_dof_load_1662 * l_sum_1729 + l_r_1797 + l_r_1798 + l_r_1799 + l_r_1800 + l_dof_load_1687 * l_sum_1738 + l_dof_load_1692 * l_sum_1740 + l_dof_load_1697 * l_mult_1731 + l_dof_load_1702 * l_mult_975 + l_dof_load_1707 * l_sum_1742 + l_dof_load_1712 * l_sum_1744;
        double l_r_1802 = l_dof_load_1687 * l_sum_1753;
        double l_r_1803 = l_dof_load_1692 * l_sum_1754;
        double l_r_1804 = l_r_1794 + l_dof_load_1622 * 0.e0;
        double l_r_1805 = l_r_1804 + l_dof_load_1627 * l_sum_1745 + l_r_1796 + l_dof_load_1637 * l_sum_1746 + l_dof_load_1642 * l_sum_1724 + l_dof_load_1647 * 0.e0 + l_dof_load_1652 * 0.e0 + l_dof_load_1657 * l_sum_1748 + l_dof_load_1662 * l_sum_1749 + l_r_1797 + l_r_1798 + l_dof_load_1677 * l_sum_1751 + l_dof_load_1682 * l_sum_1752 + l_r_1802 + l_r_1803 + l_dof_load_1697 * l_mult_1721 + l_dof_load_1702 * l_sum_1755 + l_dof_load_1707 * l_mult_1345 + l_dof_load_1712 * l_sum_1757;
        double l_r_1806 = l_r_1804 + l_r_1795 + l_dof_load_1632 * l_sum_1758 + l_dof_load_1637 * l_sum_1729 + l_dof_load_1642 * l_sum_1759 + l_dof_load_1647 * l_sum_1748 + l_dof_load_1652 * l_sum_1760 + l_dof_load_1657 * 0.e0 + l_dof_load_1662 * 0.e0 + l_dof_load_1667 * l_sum_1761 + l_dof_load_1672 * l_sum_1762 + l_r_1799 + l_r_1800 + l_r_1802 + l_r_1803 + l_dof_load_1697 * l_mult_1726 + l_dof_load_1702 * l_sum_1763 + l_dof_load_1707 * l_sum_1764 + l_dof_load_1712 * l_mult_1523;
        double l_r_1807 = 0.e0 * l_r_1773;
        double l_r_1808 = 0.e0 * l_r_1787;
        double l_r_1809 = 0.e0 * l_r_1801;
        double l_r_1810 = l_r_1807 + l_r_1808;
        double l_r_1811 = l_r_1810 + l_r_1809;
        double l_r_1812 = 0.e0 * l_r_1777;
        double l_r_1813 = 0.e0 * l_r_1791;
        double l_r_1814 = 0.e0 * l_r_1805;
        double l_r_1815 = l_r_1812 + l_r_1813;
        double l_r_1816 = l_r_1815 + l_r_1814;
        double l_r_1817 = 0.e0 * l_r_1778;
        double l_r_1818 = 0.e0 * l_r_1792;
        double l_r_1819 = 0.e0 * l_r_1806;
        double l_r_1820 = l_r_1817 + l_r_1818;
        double l_r_1821 = l_r_1820 + l_r_1819;
        double l_r_1822 = l_r_1810 + -0.1e1 * l_r_1801;
        double l_r_1823 = l_r_1815 + -0.1e1 * l_r_1805;
        double l_r_1824 = l_r_1820 + -0.1e1 * l_r_1806;
        double l_r_1825 = l_r_1807 + 0.1e1 * l_r_1787 + l_r_1809;
        double l_r_1826 = l_r_1812 + 0.1e1 * l_r_1791 + l_r_1814;
        double l_r_1827 = l_r_1817 + 0.1e1 * l_r_1792 + l_r_1819;
        double l_r_1828 = l_r_1810 + 0.1e1 * l_r_1801;
        double l_r_1829 = l_r_1815 + 0.1e1 * l_r_1805;
        double l_r_1830 = l_r_1820 + 0.1e1 * l_r_1806;
        double l_r_1831 = -0.1e1 * l_r_1773 + l_r_1808 + l_r_1809;
        double l_r_1832 = -0.1e1 * l_r_1777 + l_r_1813 + l_r_1814;
        double l_r_1833 = -0.1e1 * l_r_1778 + l_r_1818 + l_r_1819;
        double l_r_1834 = l_r_1807 + -0.1e1 * l_r_1787 + l_r_1809;
        double l_r_1835 = l_r_1812 + -0.1e1 * l_r_1791 + l_r_1814;
        double l_r_1836 = l_r_1817 + -0.1e1 * l_r_1792 + l_r_1819;
        double l_r_1837 = 0.1e1 * l_r_1773 + l_r_1808 + l_r_1809;
        double l_r_1838 = 0.1e1 * l_r_1777 + l_r_1813 + l_r_1814;
        double l_r_1839 = 0.1e1 * l_r_1778 + l_r_1818 + l_r_1819;
        double l_r_1840 = l_r_1773 * l_r_1816 + l_r_1787 * l_r_1829 + l_r_1801 * l_r_1835;
        double l_r_1841 = l_r_1773 * l_r_1821 + l_r_1787 * l_r_1830 + l_r_1801 * l_r_1836;
        double l_r_1842 = l_r_1773 * l_r_1823 + l_r_1787 * l_r_1816 + l_r_1801 * l_r_1838;
        double l_r_1843 = l_r_1773 * l_r_1824 + l_r_1787 * l_r_1821 + l_r_1801 * l_r_1839;
        double l_r_1844 = l_r_1773 * l_r_1826 + l_r_1787 * l_r_1832 + l_r_1801 * l_r_1816;
        double l_r_1845 = l_r_1773 * l_r_1827 + l_r_1787 * l_r_1833 + l_r_1801 * l_r_1821;
        double l_r_1846 = l_r_1777 * l_r_1811 + l_r_1791 * l_r_1828 + l_r_1805 * l_r_1834;
        double l_r_1847 = l_r_1777 * l_r_1821 + l_r_1791 * l_r_1830 + l_r_1805 * l_r_1836;
        double l_r_1848 = l_r_1777 * l_r_1822 + l_r_1791 * l_r_1811 + l_r_1805 * l_r_1837;
        double l_r_1849 = l_r_1777 * l_r_1824 + l_r_1791 * l_r_1821 + l_r_1805 * l_r_1839;
        double l_r_1850 = l_r_1777 * l_r_1825 + l_r_1791 * l_r_1831 + l_r_1805 * l_r_1811;
        double l_r_1851 = l_r_1777 * l_r_1827 + l_r_1791 * l_r_1833 + l_r_1805 * l_r_1821;
        double l_r_1852 = l_r_1778 * l_r_1811 + l_r_1792 * l_r_1828 + l_r_1806 * l_r_1834;
        double l_r_1853 = l_r_1778 * l_r_1816 + l_r_1792 * l_r_1829 + l_r_1806 * l_r_1835;
        double l_r_1854 = l_r_1778 * l_r_1822 + l_r_1792 * l_r_1811 + l_r_1806 * l_r_1837;
        double l_r_1855 = l_r_1778 * l_r_1823 + l_r_1792 * l_r_1816 + l_r_1806 * l_r_1838;
        double l_r_1856 = l_r_1778 * l_r_1825 + l_r_1792 * l_r_1831 + l_r_1806 * l_r_1811;
        double l_r_1857 = l_r_1778 * l_r_1826 + l_r_1792 * l_r_1832 + l_r_1806 * l_r_1816;
        vec3 v_1858 = vcons3(l_r_1777, l_r_1791, l_r_1805);
        double l_r_1859 = 0.e0 * (l_r_1773 * l_r_1811 + l_r_1787 * l_r_1828 + l_r_1801 * l_r_1834);
        double l_r_1860 = 0.e0 * l_r_1841;
        double l_r_1861 = 0.e0 * l_r_1846;
        double l_r_1862 = 0.e0 * (l_r_1777 * l_r_1816 + l_r_1791 * l_r_1829 + l_r_1805 * l_r_1835);
        double l_r_1863 = 0.e0 * l_r_1852;
        double l_r_1864 = 0.e0 * (l_r_1778 * l_r_1821 + l_r_1792 * l_r_1830 + l_r_1806 * l_r_1836);
        double l_r_1865 = l_r_1859 + 0.e0 * l_r_1840;
        double l_r_1866 = 0.e0 * (l_r_1773 * l_r_1822 + l_r_1787 * l_r_1811 + l_r_1801 * l_r_1837);
        double l_r_1867 = 0.e0 * l_r_1843;
        double l_r_1868 = 0.e0 * l_r_1848;
        double l_r_1869 = 0.e0 * (l_r_1777 * l_r_1823 + l_r_1791 * l_r_1816 + l_r_1805 * l_r_1838);
        double l_r_1870 = 0.e0 * l_r_1854;
        double l_r_1871 = 0.e0 * (l_r_1778 * l_r_1824 + l_r_1792 * l_r_1821 + l_r_1806 * l_r_1839);
        double l_r_1872 = l_r_1866 + 0.e0 * l_r_1842;
        double l_r_1873 = 0.e0 * (l_r_1773 * l_r_1825 + l_r_1787 * l_r_1831 + l_r_1801 * l_r_1811);
        double l_r_1874 = 0.e0 * l_r_1845;
        double l_r_1875 = 0.e0 * l_r_1850;
        double l_r_1876 = 0.e0 * (l_r_1777 * l_r_1826 + l_r_1791 * l_r_1832 + l_r_1805 * l_r_1816);
        double l_r_1877 = 0.e0 * l_r_1856;
        double l_r_1878 = 0.e0 * (l_r_1778 * l_r_1827 + l_r_1792 * l_r_1833 + l_r_1806 * l_r_1821);
        double l_r_1879 = l_r_1873 + 0.e0 * l_r_1844;
        double l_r_1880 = 0.e0 * l_r_1847;
        double l_r_1881 = 0.e0 * l_r_1853;
        double l_r_1882 = 0.e0 * l_r_1849;
        double l_r_1883 = 0.e0 * l_r_1855;
        double l_r_1884 = 0.e0 * l_r_1851;
        double l_r_1885 = 0.e0 * l_r_1857;
        double l_op1_e3_l_26_1886 = 0.2e1 * vdot3(vcons3(l_r_1773, l_r_1787, l_r_1801),
            vcons3(vdot3(v_1858, vcons3(l_r_1821, l_r_1830, l_r_1836)),
                vdot3(v_1858, vcons3(l_r_1824, l_r_1821, l_r_1839)),
                vdot3(v_1858, vcons3(l_r_1827, l_r_1833, l_r_1821))));
        vec3 v_1887 = -vcons3(
            l_vdot_1612 * ((l_r_1865 + l_r_1860 + l_r_1861 + l_r_1862 + 0.1e1 * l_r_1847 + l_r_1863 + -0.1e1 * l_r_1853 + l_r_1864) / l_op1_e3_l_26_1886) + l_vdot_1613 * ((l_r_1865 + -0.1e1 * l_r_1841 + l_r_1861 + l_r_1862 + l_r_1880 + 0.1e1 * l_r_1852 + l_r_1881 + l_r_1864) / l_op1_e3_l_26_1886) + l_vdot_1614 * ((l_r_1859 + 0.1e1 * l_r_1840 + l_r_1860 + -0.1e1 * l_r_1846 + l_r_1862 + l_r_1880 + l_r_1863 + l_r_1881 + l_r_1864) / l_op1_e3_l_26_1886),
            l_vdot_1612 * ((l_r_1872 + l_r_1867 + l_r_1868 + l_r_1869 + 0.1e1 * l_r_1849 + l_r_1870 + -0.1e1 * l_r_1855 + l_r_1871) / l_op1_e3_l_26_1886) + l_vdot_1613 * ((l_r_1872 + -0.1e1 * l_r_1843 + l_r_1868 + l_r_1869 + l_r_1882 + 0.1e1 * l_r_1854 + l_r_1883 + l_r_1871) / l_op1_e3_l_26_1886) + l_vdot_1614 * ((l_r_1866 + 0.1e1 * l_r_1842 + l_r_1867 + -0.1e1 * l_r_1848 + l_r_1869 + l_r_1882 + l_r_1870 + l_r_1883 + l_r_1871) / l_op1_e3_l_26_1886),
            l_vdot_1612 * ((l_r_1879 + l_r_1874 + l_r_1875 + l_r_1876 + 0.1e1 * l_r_1851 + l_r_1877 + -0.1e1 * l_r_1857 + l_r_1878) / l_op1_e3_l_26_1886) + l_vdot_1613 * ((l_r_1879 + -0.1e1 * l_r_1845 + l_r_1875 + l_r_1876 + l_r_1884 + 0.1e1 * l_r_1856 + l_r_1885 + l_r_1878) / l_op1_e3_l_26_1886) + l_vdot_1614 * ((l_r_1873 + 0.1e1 * l_r_1844 + l_r_1874 + -0.1e1 * l_r_1850 + l_r_1876 + l_r_1884 + l_r_1877 + l_r_1885 + l_r_1878) / l_op1_e3_l_26_1886));
        pthread_mutex_lock(&wrld->_sched->_prLock);
        wrld->print() << l__t_394 << "," << std::flush;
        pthread_mutex_unlock(&wrld->_sched->_prLock);
        v_1888 = vscale3(0.1e1 / std::sqrt(vdot3(v_1887, v_1887)), v_1887);
    }
    else {
        pthread_mutex_lock(&wrld->_sched->_prLock);
        wrld->print() << "Error at input pos:" << tensor_ref_3(self->sv_xp) << "\n" << std::flush;
        pthread_mutex_unlock(&wrld->_sched->_prLock);
        v_1888 = vload3(tensor_ref_3(self->sv_normal).addr(0));
    }
    vpack3(self->sv_normal, v_1888);
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
        auto i_x_1890 = *it_3;
        mesh_pos_mesh_t l__t_1891 = fn_findPos(glob->gv_meshData, i_x_1890);
        normal_init(this->_strands.strand(ix), l__t_1891, i_x_1890);
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


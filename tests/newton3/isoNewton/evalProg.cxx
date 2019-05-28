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
    struct mesh_cell_mesh_t {
        int32_t cell;
        mesh_t mesh;
        std::ostream & operator<< (std::ostream & os)
        {
            return os << this->cell;
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
    bool gv_0space0391_intermedateGlobal;
    bool gv_0data0393_intermedateGlobal;
    bool gv_ipos;
} defined_inputs;
struct globals {
    mesh_t gv_meshData;
    fns_t gv_0space0391_intermedateGlobal;
    func_t gv_0data0393_intermedateGlobal;
    diderot::dynseq< tensor_3 > gv_ipos;
    func_t gv_data;
    ~globals () { }
};
struct normal_strand {
    tensor_3 sv_normal;
    mesh_pos_mesh_t sv_pos0;
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
    wrld->_definedInp.gv_0space0391_intermedateGlobal = true;
    std::memcpy(&wrld->_globals->gv_0space0391_intermedateGlobal, v, sizeof(fns_t));
    return false;
}
extern "C" bool evalProg_input_set_data (evalProg_world_t *cWrld, void *v)
{
    world *wrld = reinterpret_cast<world *>(cWrld);
    wrld->_definedInp.gv_0data0393_intermedateGlobal = true;
    std::memcpy(&wrld->_globals->gv_0data0393_intermedateGlobal, v, sizeof(func_t));
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
    if (!wrld->_definedInp.gv_0space0391_intermedateGlobal) {
        biffMsgAdd(wrld->_errors, "undefined input \"space\"\n");
        return true;
    }
    if (!wrld->_definedInp.gv_0data0393_intermedateGlobal) {
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
    wrld->_definedInp.gv_0space0391_intermedateGlobal = false;
    wrld->_definedInp.gv_0data0393_intermedateGlobal = false;
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
    glob->gv_data = glob->gv_0data0393_intermedateGlobal.loadFem(
        glob->gv_0space0391_intermedateGlobal.loadFem(glob->gv_meshData));
    return false;
}
static void normal_init (normal_strand *self, mesh_pos_mesh_t p_pos0_391)
{
    self->sv_normal[0] = 0.e0;
    self->sv_normal[1] = 0.e0;
    self->sv_normal[2] = 0.e0;
    self->sv_pos0 = p_pos0_391;
}
static diderot::strand_status normal_update (world *wrld, globals *glob, normal_strand *self)
{
    vec3 v_770;
    if (self->sv_pos0.valid) {
        mesh_cell_mesh_t l__t_393 = makeFem(self->sv_pos0.mesh, self->sv_pos0.cell);
        func_cell_func_t l__t_394 = makeFem(glob->gv_data, l__t_393.cell);
        tensor_ref_3 l_evalPoint_395 = self->sv_pos0.refPos;
        func_t l__t_396 = (l__t_394.func);
        fns_t l__t_397 = l__t_396.space;
        mesh_t l__t_398 = l__t_397.mesh;
        int32_t l_mulRes_399 = l__t_394.cell * 20;
        int32_t l_addRes_400 = l_mulRes_399 + 1;
        int32_t l_addRes_401 = l_mulRes_399 + 2;
        int32_t l_addRes_402 = l_mulRes_399 + 3;
        int32_t l_addRes_403 = l_mulRes_399 + 4;
        int32_t l_addRes_404 = l_mulRes_399 + 5;
        int32_t l_addRes_405 = l_mulRes_399 + 6;
        int32_t l_addRes_406 = l_mulRes_399 + 7;
        int32_t l_addRes_407 = l_mulRes_399 + 8;
        int32_t l_addRes_408 = l_mulRes_399 + 9;
        int32_t l_addRes_409 = l_mulRes_399 + 10;
        int32_t l_addRes_410 = l_mulRes_399 + 11;
        int32_t l_addRes_411 = l_mulRes_399 + 12;
        int32_t l_addRes_412 = l_mulRes_399 + 13;
        int32_t l_addRes_413 = l_mulRes_399 + 14;
        int32_t l_addRes_414 = l_mulRes_399 + 15;
        int32_t l_addRes_415 = l_mulRes_399 + 16;
        int32_t l_addRes_416 = l_mulRes_399 + 17;
        int32_t l_addRes_417 = l_mulRes_399 + 18;
        int32_t l_addRes_418 = l_mulRes_399 + 19;
        int32_t t_419 = l__t_397.indexMap[l_mulRes_399];
        int32_t t_420 = l__t_397.indexMap[l_addRes_400];
        int32_t t_421 = l__t_397.indexMap[l_addRes_401];
        int32_t t_422 = l__t_397.indexMap[l_addRes_402];
        int32_t t_423 = l__t_397.indexMap[l_addRes_403];
        int32_t t_424 = l__t_397.indexMap[l_addRes_404];
        int32_t t_425 = l__t_397.indexMap[l_addRes_405];
        int32_t t_426 = l__t_397.indexMap[l_addRes_406];
        int32_t t_427 = l__t_397.indexMap[l_addRes_407];
        int32_t t_428 = l__t_397.indexMap[l_addRes_408];
        int32_t t_429 = l__t_397.indexMap[l_addRes_409];
        int32_t t_430 = l__t_397.indexMap[l_addRes_410];
        int32_t t_431 = l__t_397.indexMap[l_addRes_411];
        int32_t t_432 = l__t_397.indexMap[l_addRes_412];
        int32_t t_433 = l__t_397.indexMap[l_addRes_413];
        int32_t t_434 = l__t_397.indexMap[l_addRes_414];
        int32_t t_435 = l__t_397.indexMap[l_addRes_415];
        int32_t t_436 = l__t_397.indexMap[l_addRes_416];
        int32_t t_437 = l__t_397.indexMap[l_addRes_417];
        int32_t t_438 = l__t_397.indexMap[l_addRes_418];
        double t_439 = l__t_396.coordMap[1 * t_438];
        double t_440 = l__t_396.coordMap[1 * t_437];
        double t_441 = l__t_396.coordMap[1 * t_436];
        double t_442 = l__t_396.coordMap[1 * t_435];
        double t_443 = l__t_396.coordMap[1 * t_434];
        double t_444 = l__t_396.coordMap[1 * t_433];
        double t_445 = l__t_396.coordMap[1 * t_432];
        double t_446 = l__t_396.coordMap[1 * t_431];
        double t_447 = l__t_396.coordMap[1 * t_430];
        double t_448 = l__t_396.coordMap[1 * t_429];
        double t_449 = l__t_396.coordMap[1 * t_428];
        double t_450 = l__t_396.coordMap[1 * t_427];
        double t_451 = l__t_396.coordMap[1 * t_426];
        double t_452 = l__t_396.coordMap[1 * t_425];
        double t_453 = l__t_396.coordMap[1 * t_424];
        double t_454 = l__t_396.coordMap[1 * t_423];
        double t_455 = l__t_396.coordMap[1 * t_422];
        double t_456 = l__t_396.coordMap[1 * t_421];
        double t_457 = l__t_396.coordMap[1 * t_420];
        double t_458 = l__t_396.coordMap[1 * t_419];
        vec4 v_459 = vcons4(t_458, t_457, t_456, t_455);
        vec4 v_460 = vcons4(t_454, t_453, t_452, t_451);
        vec4 v_461 = vcons4(t_450, t_449, t_448, t_447);
        vec4 v_462 = vcons4(t_446, t_445, t_444, t_443);
        vec4 v_463 = vcons4(t_442, t_441, t_440, t_439);
        double l_varAcc_464 = l_evalPoint_395[0];
        double l_varAcc_465 = l_evalPoint_395[1];
        double l_varAcc_466 = l_evalPoint_395[2];
        double l_prod_467 = 0.1e1 * 0.1e1;
        double l_prod_468 = l_varAcc_464 * l_varAcc_464 * l_prod_467;
        double l_prod_469 = l_varAcc_465 * 0.1e1;
        double l_prod_470 = l_varAcc_464 * l_prod_469;
        double l_prod_471 = 0.1e1 * l_varAcc_466;
        double l_prod_472 = l_varAcc_464 * l_prod_471;
        double l_prod_473 = l_varAcc_464 * l_prod_467;
        double l_prod_474 = 0.1e1 * (l_varAcc_465 * l_varAcc_465 * 0.1e1);
        double l_prod_475 = 0.1e1 * (l_varAcc_465 * l_varAcc_466);
        double l_prod_476 = 0.1e1 * l_prod_469;
        double l_prod_477 = 0.1e1 * (0.1e1 * (l_varAcc_466 * l_varAcc_466));
        double l_prod_478 = 0.1e1 * l_prod_471;
        double l_prod_479 = 0.1e1 * l_prod_467;
        double l_mult_480 = -0.135e2 * l_prod_477;
        double l_mult_481 = -0.27e2 * l_prod_475;
        double l_mult_482 = -0.135e2 * l_prod_474;
        double l_mult_483 = -0.27e2 * l_prod_472;
        double l_mult_484 = -0.27e2 * l_prod_470;
        double l_mult_485 = -0.135e2 * l_prod_468;
        double l_sum_486 = -0.55e1 * l_prod_479 + (0.18e2 * l_prod_478 + (l_mult_480 + (0.18e2 * l_prod_476 + (l_mult_481 + (l_mult_482 + (0.18e2 * l_prod_473 + (l_mult_483 + (l_mult_484 + l_mult_485))))))));
        double l_mult_487 = 0.1e1 * l_prod_479;
        double l_mult_488 = 0.135e2 * l_prod_468;
        double l_sum_489 = l_mult_487 + (-0.9e1 * l_prod_473 + l_mult_488);
        double l_mult_490 = -0.45e1 * l_prod_478;
        double l_mult_491 = 0.27e2 * l_prod_472;
        double l_sum_492 = l_mult_490 + l_mult_491;
        double l_mult_493 = 0.135e2 * l_prod_477;
        double l_sum_494 = l_mult_490 + l_mult_493;
        double l_mult_495 = -0.45e1 * l_prod_476;
        double l_mult_496 = 0.27e2 * l_prod_470;
        double l_sum_497 = l_mult_495 + l_mult_496;
        double l_mult_498 = 0.135e2 * l_prod_474;
        double l_sum_499 = l_mult_495 + l_mult_498;
        double l_mult_500 = -0.225e2 * l_prod_478;
        double l_mult_501 = 0.27e2 * l_prod_475;
        double l_sum_502 = l_mult_500 + (0.27e2 * l_prod_477 + (l_mult_501 + l_mult_491));
        double l_mult_503 = 0.45e1 * l_prod_478;
        double l_sum_504 = l_mult_503 + l_mult_480;
        double l_mult_505 = -0.225e2 * l_prod_476;
        double l_sum_506 = l_mult_505 + (l_mult_501 + (0.27e2 * l_prod_474 + l_mult_496));
        double l_mult_507 = 0.45e1 * l_prod_476;
        double l_sum_508 = l_mult_507 + l_mult_482;
        double l_mult_509 = 0.9e1 * l_prod_479;
        double l_mult_510 = 0.54e2 * l_prod_472;
        double l_mult_511 = 0.54e2 * l_prod_470;
        double l_sum_512 = l_mult_509 + (l_mult_500 + (l_mult_493 + (l_mult_505 + (l_mult_501 + (l_mult_498 + (-0.45e2 * l_prod_473 + (l_mult_510 + (l_mult_511 + 0.405e2 * l_prod_468))))))));
        double l_mult_513 = -0.45e1 * l_prod_479;
        double l_sum_514 = l_mult_513 + (l_mult_503 + (l_mult_507 + (0.36e2 * l_prod_473 + (l_mult_483 + (l_mult_484 + -0.405e2 * l_prod_468)))));
        double l_mult_515 = 0.27e2 * l_prod_478;
        double l_mult_516 = -0.27e2 * l_prod_477;
        double l_mult_517 = -0.54e2 * l_prod_472;
        double l_sum_518 = l_mult_515 + (l_mult_516 + (l_mult_481 + l_mult_517));
        double l_mult_519 = 0.27e2 * l_prod_476;
        double l_mult_520 = -0.27e2 * l_prod_474;
        double l_mult_521 = -0.54e2 * l_prod_470;
        double l_sum_522 = l_mult_519 + (l_mult_481 + (l_mult_520 + l_mult_521));
        double l_sum_523 = l_mult_487 + (-0.9e1 * l_prod_476 + l_mult_498);
        double l_sum_524 = l_mult_490 + l_mult_501;
        double l_mult_525 = -0.45e1 * l_prod_473;
        double l_sum_526 = l_mult_525 + l_mult_488;
        double l_sum_527 = l_mult_525 + l_mult_496;
        double l_mult_528 = 0.54e2 * l_prod_475;
        double l_mult_529 = -0.225e2 * l_prod_473;
        double l_sum_530 = l_mult_509 + (l_mult_500 + (l_mult_493 + (-0.45e2 * l_prod_476 + (l_mult_528 + (0.405e2 * l_prod_474 + (l_mult_529 + (l_mult_491 + (l_mult_511 + l_mult_488))))))));
        double l_mult_531 = 0.45e1 * l_prod_473;
        double l_sum_532 = l_mult_513 + (l_mult_503 + (0.36e2 * l_prod_476 + (l_mult_481 + (-0.405e2 * l_prod_474 + (l_mult_531 + l_mult_484)))));
        double l_sum_533 = l_mult_529 + (l_mult_491 + (l_mult_496 + 0.27e2 * l_prod_468));
        double l_sum_534 = l_mult_531 + l_mult_485;
        double l_mult_535 = -0.54e2 * l_prod_475;
        double l_sum_536 = l_mult_515 + (l_mult_516 + (l_mult_535 + l_mult_483));
        double l_mult_537 = 0.27e2 * l_prod_473;
        double l_mult_538 = -0.27e2 * l_prod_468;
        double l_sum_539 = l_mult_537 + (l_mult_483 + (l_mult_521 + l_mult_538));
        double l_sum_540 = l_mult_487 + (-0.9e1 * l_prod_478 + l_mult_493);
        double l_sum_541 = l_mult_495 + l_mult_501;
        double l_sum_542 = l_mult_525 + l_mult_491;
        double l_sum_543 = l_mult_509 + (-0.45e2 * l_prod_478 + (0.405e2 * l_prod_477 + (l_mult_505 + (l_mult_528 + (l_mult_498 + (l_mult_529 + (l_mult_510 + (l_mult_496 + l_mult_488))))))));
        double l_sum_544 = l_mult_513 + (0.36e2 * l_prod_478 + (-0.405e2 * l_prod_477 + (l_mult_507 + (l_mult_481 + (l_mult_531 + l_mult_483)))));
        double l_sum_545 = l_mult_519 + (l_mult_535 + (l_mult_520 + l_mult_484));
        double l_sum_546 = l_mult_537 + (l_mult_517 + (l_mult_484 + l_mult_538));
        double l_vdot_547 = vdot4(v_460, vcons4(0.e0, 0.e0, l_sum_492, l_sum_494)) + (vdot4(v_461,
            vcons4(l_sum_497, l_sum_499, l_sum_502, l_sum_504)) + (vdot4(v_462,
            vcons4(l_sum_506, l_sum_508, l_sum_512, l_sum_514)) + (vdot4(v_463,
            vcons4(l_mult_501, l_mult_481, l_sum_518, l_sum_522)) + vdot4(v_459,
            vcons4(l_sum_486, l_sum_489, 0.e0, 0.e0)))));
        double l_vdot_548 = vdot4(v_460, vcons4(l_sum_524, l_sum_494, 0.e0, 0.e0)) + (vdot4(v_461,
            vcons4(l_sum_526, l_sum_527, l_sum_502, l_sum_504)) + (vdot4(v_462,
            vcons4(l_sum_530, l_sum_532, l_sum_533, l_sum_534)) + (vdot4(v_463,
            vcons4(l_mult_491, l_sum_536, l_mult_483, l_sum_539)) + vdot4(v_459,
            vcons4(l_sum_486, 0.e0, l_sum_523, 0.e0)))));
        double l_vdot_549 = vdot4(v_460, vcons4(l_sum_499, l_sum_541, l_sum_526, l_sum_542)) + (vdot4(v_461,
            vcons4(0.e0, 0.e0, l_sum_543, l_sum_544)) + (vdot4(v_462,
            vcons4(l_sum_506, l_sum_508, l_sum_533, l_sum_534)) + (vdot4(v_463,
            vcons4(l_mult_496, l_sum_545, l_sum_546, l_mult_484)) + vdot4(v_459,
            vcons4(l_sum_486, 0.e0, 0.e0, l_sum_540)))));
        int32_t t_550 = l__t_398.indexMap[l_mulRes_399];
        int32_t l_mulRes_551 = 3 * t_550;
        int32_t t_552 = l__t_398.indexMap[l_addRes_400];
        int32_t l_mulRes_553 = 3 * t_552;
        double l_dof_load_554 = l__t_398.coordMap[l_mulRes_553];
        double l_dof_load_555 = l__t_398.coordMap[1 + l_mulRes_553];
        double l_dof_load_556 = l__t_398.coordMap[2 + l_mulRes_553];
        int32_t t_557 = l__t_398.indexMap[l_addRes_401];
        int32_t l_mulRes_558 = 3 * t_557;
        double l_dof_load_559 = l__t_398.coordMap[l_mulRes_558];
        double l_dof_load_560 = l__t_398.coordMap[1 + l_mulRes_558];
        double l_dof_load_561 = l__t_398.coordMap[2 + l_mulRes_558];
        int32_t t_562 = l__t_398.indexMap[l_addRes_402];
        int32_t l_mulRes_563 = 3 * t_562;
        double l_dof_load_564 = l__t_398.coordMap[l_mulRes_563];
        double l_dof_load_565 = l__t_398.coordMap[1 + l_mulRes_563];
        double l_dof_load_566 = l__t_398.coordMap[2 + l_mulRes_563];
        int32_t t_567 = l__t_398.indexMap[l_addRes_403];
        int32_t l_mulRes_568 = 3 * t_567;
        double l_dof_load_569 = l__t_398.coordMap[l_mulRes_568];
        double l_dof_load_570 = l__t_398.coordMap[1 + l_mulRes_568];
        double l_dof_load_571 = l__t_398.coordMap[2 + l_mulRes_568];
        int32_t t_572 = l__t_398.indexMap[l_addRes_404];
        int32_t l_mulRes_573 = 3 * t_572;
        double l_dof_load_574 = l__t_398.coordMap[l_mulRes_573];
        double l_dof_load_575 = l__t_398.coordMap[1 + l_mulRes_573];
        double l_dof_load_576 = l__t_398.coordMap[2 + l_mulRes_573];
        int32_t t_577 = l__t_398.indexMap[l_addRes_405];
        int32_t l_mulRes_578 = 3 * t_577;
        double l_dof_load_579 = l__t_398.coordMap[l_mulRes_578];
        double l_dof_load_580 = l__t_398.coordMap[1 + l_mulRes_578];
        double l_dof_load_581 = l__t_398.coordMap[2 + l_mulRes_578];
        int32_t t_582 = l__t_398.indexMap[l_addRes_406];
        int32_t l_mulRes_583 = 3 * t_582;
        double l_dof_load_584 = l__t_398.coordMap[l_mulRes_583];
        double l_dof_load_585 = l__t_398.coordMap[1 + l_mulRes_583];
        double l_dof_load_586 = l__t_398.coordMap[2 + l_mulRes_583];
        int32_t t_587 = l__t_398.indexMap[l_addRes_407];
        int32_t l_mulRes_588 = 3 * t_587;
        double l_dof_load_589 = l__t_398.coordMap[l_mulRes_588];
        double l_dof_load_590 = l__t_398.coordMap[1 + l_mulRes_588];
        double l_dof_load_591 = l__t_398.coordMap[2 + l_mulRes_588];
        int32_t t_592 = l__t_398.indexMap[l_addRes_408];
        int32_t l_mulRes_593 = 3 * t_592;
        double l_dof_load_594 = l__t_398.coordMap[l_mulRes_593];
        double l_dof_load_595 = l__t_398.coordMap[1 + l_mulRes_593];
        double l_dof_load_596 = l__t_398.coordMap[2 + l_mulRes_593];
        int32_t t_597 = l__t_398.indexMap[l_addRes_409];
        int32_t l_mulRes_598 = 3 * t_597;
        double l_dof_load_599 = l__t_398.coordMap[l_mulRes_598];
        double l_dof_load_600 = l__t_398.coordMap[1 + l_mulRes_598];
        double l_dof_load_601 = l__t_398.coordMap[2 + l_mulRes_598];
        int32_t t_602 = l__t_398.indexMap[l_addRes_410];
        int32_t l_mulRes_603 = 3 * t_602;
        double l_dof_load_604 = l__t_398.coordMap[l_mulRes_603];
        double l_dof_load_605 = l__t_398.coordMap[1 + l_mulRes_603];
        double l_dof_load_606 = l__t_398.coordMap[2 + l_mulRes_603];
        int32_t t_607 = l__t_398.indexMap[l_addRes_411];
        int32_t l_mulRes_608 = 3 * t_607;
        double l_dof_load_609 = l__t_398.coordMap[l_mulRes_608];
        double l_dof_load_610 = l__t_398.coordMap[1 + l_mulRes_608];
        double l_dof_load_611 = l__t_398.coordMap[2 + l_mulRes_608];
        int32_t t_612 = l__t_398.indexMap[l_addRes_412];
        int32_t l_mulRes_613 = 3 * t_612;
        double l_dof_load_614 = l__t_398.coordMap[l_mulRes_613];
        double l_dof_load_615 = l__t_398.coordMap[1 + l_mulRes_613];
        double l_dof_load_616 = l__t_398.coordMap[2 + l_mulRes_613];
        int32_t t_617 = l__t_398.indexMap[l_addRes_413];
        int32_t l_mulRes_618 = 3 * t_617;
        double l_dof_load_619 = l__t_398.coordMap[l_mulRes_618];
        double l_dof_load_620 = l__t_398.coordMap[1 + l_mulRes_618];
        double l_dof_load_621 = l__t_398.coordMap[2 + l_mulRes_618];
        int32_t t_622 = l__t_398.indexMap[l_addRes_414];
        int32_t l_mulRes_623 = 3 * t_622;
        double l_dof_load_624 = l__t_398.coordMap[l_mulRes_623];
        double l_dof_load_625 = l__t_398.coordMap[1 + l_mulRes_623];
        double l_dof_load_626 = l__t_398.coordMap[2 + l_mulRes_623];
        int32_t t_627 = l__t_398.indexMap[l_addRes_415];
        int32_t l_mulRes_628 = 3 * t_627;
        double l_dof_load_629 = l__t_398.coordMap[l_mulRes_628];
        double l_dof_load_630 = l__t_398.coordMap[1 + l_mulRes_628];
        double l_dof_load_631 = l__t_398.coordMap[2 + l_mulRes_628];
        int32_t t_632 = l__t_398.indexMap[l_addRes_416];
        int32_t l_mulRes_633 = 3 * t_632;
        double l_dof_load_634 = l__t_398.coordMap[l_mulRes_633];
        double l_dof_load_635 = l__t_398.coordMap[1 + l_mulRes_633];
        double l_dof_load_636 = l__t_398.coordMap[2 + l_mulRes_633];
        int32_t t_637 = l__t_398.indexMap[l_addRes_417];
        int32_t l_mulRes_638 = 3 * t_637;
        double l_dof_load_639 = l__t_398.coordMap[l_mulRes_638];
        double l_dof_load_640 = l__t_398.coordMap[1 + l_mulRes_638];
        double l_dof_load_641 = l__t_398.coordMap[2 + l_mulRes_638];
        int32_t t_642 = l__t_398.indexMap[l_addRes_418];
        int32_t l_mulRes_643 = 3 * t_642;
        double l_dof_load_644 = l__t_398.coordMap[l_mulRes_643];
        double l_dof_load_645 = l__t_398.coordMap[1 + l_mulRes_643];
        double l_dof_load_646 = l__t_398.coordMap[2 + l_mulRes_643];
        double t_647 = l__t_398.coordMap[l_mulRes_551];
        double l_r_648 = t_647 * l_sum_486;
        double l_r_649 = l_dof_load_559 * 0.e0;
        double l_r_650 = l_dof_load_564 * 0.e0;
        double l_r_651 = l_dof_load_599 * l_sum_502;
        double l_r_652 = l_dof_load_604 * l_sum_504;
        double l_r_653 = l_dof_load_609 * l_sum_506;
        double l_r_654 = l_dof_load_614 * l_sum_508;
        double l_r_655 = l_r_648 + l_dof_load_554 * l_sum_489 + l_r_649 + l_r_650 + l_dof_load_569 * 0.e0 + l_dof_load_574 * 0.e0 + l_dof_load_579 * l_sum_492 + l_dof_load_584 * l_sum_494 + l_dof_load_589 * l_sum_497 + l_dof_load_594 * l_sum_499 + l_r_651 + l_r_652 + l_r_653 + l_r_654 + l_dof_load_619 * l_sum_512 + l_dof_load_624 * l_sum_514 + l_dof_load_629 * l_mult_501 + l_dof_load_634 * l_mult_481 + l_dof_load_639 * l_sum_518 + l_dof_load_644 * l_sum_522;
        double l_r_656 = l_dof_load_619 * l_sum_533;
        double l_r_657 = l_dof_load_624 * l_sum_534;
        double l_r_658 = l_r_648 + l_dof_load_554 * 0.e0;
        double l_r_659 = l_r_658 + l_dof_load_559 * l_sum_523 + l_r_650 + l_dof_load_569 * l_sum_524 + l_dof_load_574 * l_sum_494 + l_dof_load_579 * 0.e0 + l_dof_load_584 * 0.e0 + l_dof_load_589 * l_sum_526 + l_dof_load_594 * l_sum_527 + l_r_651 + l_r_652 + l_dof_load_609 * l_sum_530 + l_dof_load_614 * l_sum_532 + l_r_656 + l_r_657 + l_dof_load_629 * l_mult_491 + l_dof_load_634 * l_sum_536 + l_dof_load_639 * l_mult_483 + l_dof_load_644 * l_sum_539;
        double l_r_660 = l_r_658 + l_r_649 + l_dof_load_564 * l_sum_540 + l_dof_load_569 * l_sum_499 + l_dof_load_574 * l_sum_541 + l_dof_load_579 * l_sum_526 + l_dof_load_584 * l_sum_542 + l_dof_load_589 * 0.e0 + l_dof_load_594 * 0.e0 + l_dof_load_599 * l_sum_543 + l_dof_load_604 * l_sum_544 + l_r_653 + l_r_654 + l_r_656 + l_r_657 + l_dof_load_629 * l_mult_496 + l_dof_load_634 * l_sum_545 + l_dof_load_639 * l_sum_546 + l_dof_load_644 * l_mult_484;
        double t_661 = l__t_398.coordMap[1 + l_mulRes_551];
        double l_r_662 = t_661 * l_sum_486;
        double l_r_663 = l_dof_load_560 * 0.e0;
        double l_r_664 = l_dof_load_565 * 0.e0;
        double l_r_665 = l_dof_load_600 * l_sum_502;
        double l_r_666 = l_dof_load_605 * l_sum_504;
        double l_r_667 = l_dof_load_610 * l_sum_506;
        double l_r_668 = l_dof_load_615 * l_sum_508;
        double l_r_669 = l_r_662 + l_dof_load_555 * l_sum_489 + l_r_663 + l_r_664 + l_dof_load_570 * 0.e0 + l_dof_load_575 * 0.e0 + l_dof_load_580 * l_sum_492 + l_dof_load_585 * l_sum_494 + l_dof_load_590 * l_sum_497 + l_dof_load_595 * l_sum_499 + l_r_665 + l_r_666 + l_r_667 + l_r_668 + l_dof_load_620 * l_sum_512 + l_dof_load_625 * l_sum_514 + l_dof_load_630 * l_mult_501 + l_dof_load_635 * l_mult_481 + l_dof_load_640 * l_sum_518 + l_dof_load_645 * l_sum_522;
        double l_r_670 = l_dof_load_620 * l_sum_533;
        double l_r_671 = l_dof_load_625 * l_sum_534;
        double l_r_672 = l_r_662 + l_dof_load_555 * 0.e0;
        double l_r_673 = l_r_672 + l_dof_load_560 * l_sum_523 + l_r_664 + l_dof_load_570 * l_sum_524 + l_dof_load_575 * l_sum_494 + l_dof_load_580 * 0.e0 + l_dof_load_585 * 0.e0 + l_dof_load_590 * l_sum_526 + l_dof_load_595 * l_sum_527 + l_r_665 + l_r_666 + l_dof_load_610 * l_sum_530 + l_dof_load_615 * l_sum_532 + l_r_670 + l_r_671 + l_dof_load_630 * l_mult_491 + l_dof_load_635 * l_sum_536 + l_dof_load_640 * l_mult_483 + l_dof_load_645 * l_sum_539;
        double l_r_674 = l_r_672 + l_r_663 + l_dof_load_565 * l_sum_540 + l_dof_load_570 * l_sum_499 + l_dof_load_575 * l_sum_541 + l_dof_load_580 * l_sum_526 + l_dof_load_585 * l_sum_542 + l_dof_load_590 * 0.e0 + l_dof_load_595 * 0.e0 + l_dof_load_600 * l_sum_543 + l_dof_load_605 * l_sum_544 + l_r_667 + l_r_668 + l_r_670 + l_r_671 + l_dof_load_630 * l_mult_496 + l_dof_load_635 * l_sum_545 + l_dof_load_640 * l_sum_546 + l_dof_load_645 * l_mult_484;
        double t_675 = l__t_398.coordMap[2 + l_mulRes_551];
        double l_r_676 = t_675 * l_sum_486;
        double l_r_677 = l_dof_load_561 * 0.e0;
        double l_r_678 = l_dof_load_566 * 0.e0;
        double l_r_679 = l_dof_load_601 * l_sum_502;
        double l_r_680 = l_dof_load_606 * l_sum_504;
        double l_r_681 = l_dof_load_611 * l_sum_506;
        double l_r_682 = l_dof_load_616 * l_sum_508;
        double l_r_683 = l_r_676 + l_dof_load_556 * l_sum_489 + l_r_677 + l_r_678 + l_dof_load_571 * 0.e0 + l_dof_load_576 * 0.e0 + l_dof_load_581 * l_sum_492 + l_dof_load_586 * l_sum_494 + l_dof_load_591 * l_sum_497 + l_dof_load_596 * l_sum_499 + l_r_679 + l_r_680 + l_r_681 + l_r_682 + l_dof_load_621 * l_sum_512 + l_dof_load_626 * l_sum_514 + l_dof_load_631 * l_mult_501 + l_dof_load_636 * l_mult_481 + l_dof_load_641 * l_sum_518 + l_dof_load_646 * l_sum_522;
        double l_r_684 = l_dof_load_621 * l_sum_533;
        double l_r_685 = l_dof_load_626 * l_sum_534;
        double l_r_686 = l_r_676 + l_dof_load_556 * 0.e0;
        double l_r_687 = l_r_686 + l_dof_load_561 * l_sum_523 + l_r_678 + l_dof_load_571 * l_sum_524 + l_dof_load_576 * l_sum_494 + l_dof_load_581 * 0.e0 + l_dof_load_586 * 0.e0 + l_dof_load_591 * l_sum_526 + l_dof_load_596 * l_sum_527 + l_r_679 + l_r_680 + l_dof_load_611 * l_sum_530 + l_dof_load_616 * l_sum_532 + l_r_684 + l_r_685 + l_dof_load_631 * l_mult_491 + l_dof_load_636 * l_sum_536 + l_dof_load_641 * l_mult_483 + l_dof_load_646 * l_sum_539;
        double l_r_688 = l_r_686 + l_r_677 + l_dof_load_566 * l_sum_540 + l_dof_load_571 * l_sum_499 + l_dof_load_576 * l_sum_541 + l_dof_load_581 * l_sum_526 + l_dof_load_586 * l_sum_542 + l_dof_load_591 * 0.e0 + l_dof_load_596 * 0.e0 + l_dof_load_601 * l_sum_543 + l_dof_load_606 * l_sum_544 + l_r_681 + l_r_682 + l_r_684 + l_r_685 + l_dof_load_631 * l_mult_496 + l_dof_load_636 * l_sum_545 + l_dof_load_641 * l_sum_546 + l_dof_load_646 * l_mult_484;
        double l_r_689 = 0.e0 * l_r_655;
        double l_r_690 = 0.e0 * l_r_669;
        double l_r_691 = 0.e0 * l_r_683;
        double l_r_692 = l_r_689 + l_r_690;
        double l_r_693 = l_r_692 + l_r_691;
        double l_r_694 = 0.e0 * l_r_659;
        double l_r_695 = 0.e0 * l_r_673;
        double l_r_696 = 0.e0 * l_r_687;
        double l_r_697 = l_r_694 + l_r_695;
        double l_r_698 = l_r_697 + l_r_696;
        double l_r_699 = 0.e0 * l_r_660;
        double l_r_700 = 0.e0 * l_r_674;
        double l_r_701 = 0.e0 * l_r_688;
        double l_r_702 = l_r_699 + l_r_700;
        double l_r_703 = l_r_702 + l_r_701;
        double l_r_704 = l_r_692 + -0.1e1 * l_r_683;
        double l_r_705 = l_r_697 + -0.1e1 * l_r_687;
        double l_r_706 = l_r_702 + -0.1e1 * l_r_688;
        double l_r_707 = l_r_689 + 0.1e1 * l_r_669 + l_r_691;
        double l_r_708 = l_r_694 + 0.1e1 * l_r_673 + l_r_696;
        double l_r_709 = l_r_699 + 0.1e1 * l_r_674 + l_r_701;
        double l_r_710 = l_r_692 + 0.1e1 * l_r_683;
        double l_r_711 = l_r_697 + 0.1e1 * l_r_687;
        double l_r_712 = l_r_702 + 0.1e1 * l_r_688;
        double l_r_713 = -0.1e1 * l_r_655 + l_r_690 + l_r_691;
        double l_r_714 = -0.1e1 * l_r_659 + l_r_695 + l_r_696;
        double l_r_715 = -0.1e1 * l_r_660 + l_r_700 + l_r_701;
        double l_r_716 = l_r_689 + -0.1e1 * l_r_669 + l_r_691;
        double l_r_717 = l_r_694 + -0.1e1 * l_r_673 + l_r_696;
        double l_r_718 = l_r_699 + -0.1e1 * l_r_674 + l_r_701;
        double l_r_719 = 0.1e1 * l_r_655 + l_r_690 + l_r_691;
        double l_r_720 = 0.1e1 * l_r_659 + l_r_695 + l_r_696;
        double l_r_721 = 0.1e1 * l_r_660 + l_r_700 + l_r_701;
        double l_r_722 = l_r_655 * l_r_698 + l_r_669 * l_r_711 + l_r_683 * l_r_717;
        double l_r_723 = l_r_655 * l_r_703 + l_r_669 * l_r_712 + l_r_683 * l_r_718;
        double l_r_724 = l_r_655 * l_r_705 + l_r_669 * l_r_698 + l_r_683 * l_r_720;
        double l_r_725 = l_r_655 * l_r_706 + l_r_669 * l_r_703 + l_r_683 * l_r_721;
        double l_r_726 = l_r_655 * l_r_708 + l_r_669 * l_r_714 + l_r_683 * l_r_698;
        double l_r_727 = l_r_655 * l_r_709 + l_r_669 * l_r_715 + l_r_683 * l_r_703;
        double l_r_728 = l_r_659 * l_r_693 + l_r_673 * l_r_710 + l_r_687 * l_r_716;
        double l_r_729 = l_r_659 * l_r_703 + l_r_673 * l_r_712 + l_r_687 * l_r_718;
        double l_r_730 = l_r_659 * l_r_704 + l_r_673 * l_r_693 + l_r_687 * l_r_719;
        double l_r_731 = l_r_659 * l_r_706 + l_r_673 * l_r_703 + l_r_687 * l_r_721;
        double l_r_732 = l_r_659 * l_r_707 + l_r_673 * l_r_713 + l_r_687 * l_r_693;
        double l_r_733 = l_r_659 * l_r_709 + l_r_673 * l_r_715 + l_r_687 * l_r_703;
        double l_r_734 = l_r_660 * l_r_693 + l_r_674 * l_r_710 + l_r_688 * l_r_716;
        double l_r_735 = l_r_660 * l_r_698 + l_r_674 * l_r_711 + l_r_688 * l_r_717;
        double l_r_736 = l_r_660 * l_r_704 + l_r_674 * l_r_693 + l_r_688 * l_r_719;
        double l_r_737 = l_r_660 * l_r_705 + l_r_674 * l_r_698 + l_r_688 * l_r_720;
        double l_r_738 = l_r_660 * l_r_707 + l_r_674 * l_r_713 + l_r_688 * l_r_693;
        double l_r_739 = l_r_660 * l_r_708 + l_r_674 * l_r_714 + l_r_688 * l_r_698;
        vec3 v_740 = vcons3(l_r_659, l_r_673, l_r_687);
        double l_r_741 = 0.e0 * (l_r_655 * l_r_693 + l_r_669 * l_r_710 + l_r_683 * l_r_716);
        double l_r_742 = 0.e0 * l_r_723;
        double l_r_743 = 0.e0 * l_r_728;
        double l_r_744 = 0.e0 * (l_r_659 * l_r_698 + l_r_673 * l_r_711 + l_r_687 * l_r_717);
        double l_r_745 = 0.e0 * l_r_734;
        double l_r_746 = 0.e0 * (l_r_660 * l_r_703 + l_r_674 * l_r_712 + l_r_688 * l_r_718);
        double l_r_747 = l_r_741 + 0.e0 * l_r_722;
        double l_r_748 = 0.e0 * (l_r_655 * l_r_704 + l_r_669 * l_r_693 + l_r_683 * l_r_719);
        double l_r_749 = 0.e0 * l_r_725;
        double l_r_750 = 0.e0 * l_r_730;
        double l_r_751 = 0.e0 * (l_r_659 * l_r_705 + l_r_673 * l_r_698 + l_r_687 * l_r_720);
        double l_r_752 = 0.e0 * l_r_736;
        double l_r_753 = 0.e0 * (l_r_660 * l_r_706 + l_r_674 * l_r_703 + l_r_688 * l_r_721);
        double l_r_754 = l_r_748 + 0.e0 * l_r_724;
        double l_r_755 = 0.e0 * (l_r_655 * l_r_707 + l_r_669 * l_r_713 + l_r_683 * l_r_693);
        double l_r_756 = 0.e0 * l_r_727;
        double l_r_757 = 0.e0 * l_r_732;
        double l_r_758 = 0.e0 * (l_r_659 * l_r_708 + l_r_673 * l_r_714 + l_r_687 * l_r_698);
        double l_r_759 = 0.e0 * l_r_738;
        double l_r_760 = 0.e0 * (l_r_660 * l_r_709 + l_r_674 * l_r_715 + l_r_688 * l_r_703);
        double l_r_761 = l_r_755 + 0.e0 * l_r_726;
        double l_r_762 = 0.e0 * l_r_729;
        double l_r_763 = 0.e0 * l_r_735;
        double l_r_764 = 0.e0 * l_r_731;
        double l_r_765 = 0.e0 * l_r_737;
        double l_r_766 = 0.e0 * l_r_733;
        double l_r_767 = 0.e0 * l_r_739;
        double l_op1_e3_l_26_768 = 0.2e1 * vdot3(vcons3(l_r_655, l_r_669, l_r_683),
            vcons3(vdot3(v_740, vcons3(l_r_703, l_r_712, l_r_718)), vdot3(v_740, vcons3(l_r_706, l_r_703, l_r_721)),
                vdot3(v_740, vcons3(l_r_709, l_r_715, l_r_703))));
        vec3 v_769 = -vcons3(
            l_vdot_547 * ((l_r_747 + l_r_742 + l_r_743 + l_r_744 + 0.1e1 * l_r_729 + l_r_745 + -0.1e1 * l_r_735 + l_r_746) / l_op1_e3_l_26_768) + l_vdot_548 * ((l_r_747 + -0.1e1 * l_r_723 + l_r_743 + l_r_744 + l_r_762 + 0.1e1 * l_r_734 + l_r_763 + l_r_746) / l_op1_e3_l_26_768) + l_vdot_549 * ((l_r_741 + 0.1e1 * l_r_722 + l_r_742 + -0.1e1 * l_r_728 + l_r_744 + l_r_762 + l_r_745 + l_r_763 + l_r_746) / l_op1_e3_l_26_768),
            l_vdot_547 * ((l_r_754 + l_r_749 + l_r_750 + l_r_751 + 0.1e1 * l_r_731 + l_r_752 + -0.1e1 * l_r_737 + l_r_753) / l_op1_e3_l_26_768) + l_vdot_548 * ((l_r_754 + -0.1e1 * l_r_725 + l_r_750 + l_r_751 + l_r_764 + 0.1e1 * l_r_736 + l_r_765 + l_r_753) / l_op1_e3_l_26_768) + l_vdot_549 * ((l_r_748 + 0.1e1 * l_r_724 + l_r_749 + -0.1e1 * l_r_730 + l_r_751 + l_r_764 + l_r_752 + l_r_765 + l_r_753) / l_op1_e3_l_26_768),
            l_vdot_547 * ((l_r_761 + l_r_756 + l_r_757 + l_r_758 + 0.1e1 * l_r_733 + l_r_759 + -0.1e1 * l_r_739 + l_r_760) / l_op1_e3_l_26_768) + l_vdot_548 * ((l_r_761 + -0.1e1 * l_r_727 + l_r_757 + l_r_758 + l_r_766 + 0.1e1 * l_r_738 + l_r_767 + l_r_760) / l_op1_e3_l_26_768) + l_vdot_549 * ((l_r_755 + 0.1e1 * l_r_726 + l_r_756 + -0.1e1 * l_r_732 + l_r_758 + l_r_766 + l_r_759 + l_r_767 + l_r_760) / l_op1_e3_l_26_768));
        pthread_mutex_lock(&wrld->_sched->_prLock);
        wrld->print() << l__t_393 << "," << std::flush;
        pthread_mutex_unlock(&wrld->_sched->_prLock);
        v_770 = vscale3(0.1e1 / std::sqrt(vdot3(v_769, v_769)), v_769);
    }
    else {
        pthread_mutex_lock(&wrld->_sched->_prLock);
        wrld->print() << "Error at input pos\n" << std::flush;
        pthread_mutex_unlock(&wrld->_sched->_prLock);
        v_770 = vload3(tensor_ref_3(self->sv_normal).addr(0));
    }
    vpack3(self->sv_normal, v_770);
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
        auto i_x_772 = *it_3;
        mesh_pos_mesh_t l__t_773 = fn_findPos(glob->gv_meshData, i_x_772);
        normal_init(this->_strands.strand(ix), l__t_773);
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


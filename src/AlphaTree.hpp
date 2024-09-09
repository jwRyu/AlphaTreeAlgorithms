#pragma once
#include <AlphaTreeConfig.hpp>
#include <HeapQueue.hpp>
#include <HierarQueue.hpp>
#include <HybridQueue.hpp>
#include <LadderQueue.hpp>
#include <PixelDissimilarity.hpp>
#include <Trie.hpp>
#include <defines.hpp>
#include <radixsort_teeninga/sort/radix_sort_parallel.h>
#include <radixsort_teeninga/sort/sort_item.h>
#include <walltime.hpp>

using namespace pmt;

// Minimum number of pixels for TSE to be used. TSE does not work well with very small images (also saving memory is not
// so important on small images).
#define TSE_MINSIZE 10000

#define TSE_A 1.3901
#define TSE_SIGMA -2.1989
#define TSE_B -0.1906
#define TSE_M 0.05

#define IMGIDX_32BITS 0
#define IMGIDX_64BITS 1

#define PIXEL_8BIT 1
#define PIXEL_16BIT 2
#define PIXEL_32BIT 4
#define PIXEL_64BIT 8
#define PIXEL_FLOAT 16
#define PIXEL_DOUBLE 32

// algcode
#define UNIONFIND 0                           // p1
#define FLOOD_HIERARQUEUE 1                   // p1,2
#define FLOOD_HIERARQUEUE_CACHE 2             // p2
#define FLOOD_TRIE 3                          // p1,2
#define FLOOD_TRIE_CACHE 4                    // p2
#define FLOOD_HEAPQUEUE 5                     // p1,2
#define FLOOD_HEAPQUEUE_CACHE 6               // p1,2
#define FLOOD_HIERARQUEUE_HYPERGRAPH 7        // p1
#define FLOOD_TRIE_HYPERGRAPH 8               // p1
#define FLOOD_HIERARHEAPQUEUE_CACHE 9         // p2
#define FLOOD_HIERARHEAPQUEUE_CACHE_HISTEQ 10 // p2
#define FLOOD_HIERARQUEUE_PAR 11              // p1
#define PILOT_RANK 12                         // p1
#define FLOOD_HIERARQUEUE_CACHE_PAR 13        // p2
#define FLOOD_HIERARHEAPQUEUE 14              // p2
#define FLOOD_LADDERQUEUE 15                  // p2
#define FLOOD_HEAPQUEUE_NAIVE 16              // p2

#define ROOTIDX -1

#define CHKRNG(var, a, b) ((var >= a) && (var < b))
#define QUANTIZE_RANK(rank, binsize) (rank) / (binsize)

template <class Pixel> class AlphaNode {
  public:
    ImgIdx area = 0;
    double alpha = std::numeric_limits<Pixel>::infinity();
    double sumPix = 0.0;
    Pixel minPix = std::numeric_limits<Pixel>::max();
    Pixel maxPix = std::numeric_limits<Pixel>::min();
    ImgIdx parentIdx = ROOTIDX;
    ImgIdx _rootIdx = ROOTIDX;

    AlphaNode() = default;
    AlphaNode(Pixel pixelVal, double alpha_, ImgIdx parentidx_ = ROOTIDX);
    AlphaNode(double alpha_, ImgIdx parentidx_ = ROOTIDX);

    inline void set(ImgIdx area_in, double level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in);
    inline void add(AlphaNode *q);
    inline void add(const AlphaNode &q);
    inline void add(const Pixel &pix_val);
    inline void copy(AlphaNode *q);
    inline void connect_to_parent(AlphaNode *pPar, ImgIdx iPar);
    void print(AlphaNode *_node);
    void print(AlphaNode *_node, int heading);

    bool operator<(const AlphaNode &other) const { return alpha < other.alpha; }
};

template <class Pixel> class RankItem {
  public:
    Pixel alpha;
    ImgIdx dimgidx;

    void operator=(const RankItem &item);
    ImgIdx get_pidx0(ImgIdx _connectivity = 4);
    ImgIdx get_pidx1(ImgIdx _width, ImgIdx _connectivity = 4);
};

template <class Pixel> class AlphaTree {
  public:
    ImgIdx _maxSize = 0;
    ImgIdx _curSize = 0;
    ImgIdx _height = 0;
    ImgIdx _width = 0;
    ImgIdx _channel = 0;
    ImgIdx _connectivity = 0;
    AlphaNode<Pixel> *_node = nullptr;
    AlphaNode<Pixel> *_nodeIn = nullptr;
    ImgIdx *_parentAry = nullptr;
    ImgIdx num_node = 0;
    ImgIdx num_node_in = 0;
    ImgIdx _rootIdx = ROOTIDX;
    double nrmsd = 0.0;

    AlphaTree() : _maxSize(0), _curSize(0), _node(0), _parentAry(0) {}
    ~AlphaTree();
    void clear();

    void BuildAlphaTree(const Pixel *img, int height_in, int width_in, int channel_in, std::string dMetric,
                        int connectivity_in, int algorithm, int numthreads, int tse, double fparam1 = 0.0,
                        double fparam2 = 0.0, int iparam1 = 0);

    void AlphaFilter(double *outimg, double alpha);
    void AreaFilter(double *outimg, double area);

    void printTree() const;
    void printGraph(const uint8_t *isVisited, const bool *isRedundant, const Pixel *img) const;
    void printGraph(const uint8_t *isVisited, const uint8_t *edge, const Pixel *img) const;
    void printParentAry() const;
    // void printAll(const uint8_t *isVisited, const bool *isRedundant, const Pixel *img) const;
    void printAll(const uint8_t *isVisited, const uint8_t *edge, const Pixel *img) const;
    void printVisit(ImgIdx p, double q) const;
    void printIsAvailable(const uint8_t *isAvailable) const;
    static void printBinary(uint8_t num);

  private:
    PixelDissimilarity<Pixel> _pixelDissim;

    void Unionfind(const Pixel *img);
    void FloodHierarQueueNoCache(const Pixel *img, int tse);
    void FloodHeapQueueNoCache(const Pixel *img);
    void FloodHeapQueueNaiveNoCache(const Pixel *img);
    void FloodHeapQueue(const Pixel *img);
    void FloodHierarQueue(const Pixel *img);
    void FloodLadderQueue(const Pixel *img, int thres = 64);
    void FloodHierarHeapQueueNoCache(const Pixel *img, double a = 12.0, double r = 0.5, int listsize = 12);
    void FloodHierarHeapQueue(const Pixel *img, double a = 12.0, double r = 0.5, int listsize = 12);
    void FloodHierHeapQueueHisteq(const Pixel *img, int listsize = 12, int a = 0);
    void FloodTrieHypergraph(const Pixel *img);
    void FloodHierarQueueHypergraph(const Pixel *img);
    void FloodTrieNoCache(const Pixel *img);
    void FloodTrie(const Pixel *img);
    void FloodHierQueueParallel(const Pixel *img, int numthreads);
    void HybridParallel(const Pixel *img, int numthreads);

    ImgIdx mergePartition(Pixel *dimg, int64_t blksz_x, int64_t blksz_y, ImgIdx npartition_x, ImgIdx npartition_y,
                          ImgIdx *subtree_cur);
    void floodPartition(const Pixel *img, const Pixel *dimg, ImgIdx startPixelIndex, int blockIndex, ImgIdx blockArea,
                        Pixel maxdiff, HierarQueue *queue, const ImgIdx *subtree_start, ImgIdx *subtree_cur,
                        uint8_t *isVisited, const uint8_t *isAvailable);
    Pixel computePartitionDifferences(const Pixel *img, Pixel *dimg, ImgIdx startPixelIndex, ImgIdx blockWidth,
                                      ImgIdx blockHeight, ImgIdx *blockDiffHist);
    std::pair<ImgIdx, ImgIdx> computePartitionSize(int numthreads);
    void setBlockDimensions(ImgIdx npartition_x, ImgIdx npartition_y, ImgIdx blksz_x, ImgIdx blksz_y, ImgIdx *blocksize,
                            ImgIdx *blockWidths, ImgIdx *blockHeights, ImgIdx *startpidx, ImgIdx *subtree_start);

    void runFloodHHPQ(ImgIdx startingPixel, const Pixel *img, double a, double r, int listsize, ImgIdx imgSize,
                      ImgIdx nredges, ImgIdx dimgSize, uint64_t numLevels, const ImgIdx *dhist, const double *dimg,
                      const uint8_t *isAvailable);

    void sortAlphaNodes();
    void markRedundant(ImgIdx imgIdx, ImgIdx eIdx, uint8_t *edgeStatus, ImgIdx *queuedEdges,
                       uint8_t *numQueuedEdges) const;
    void registerEdge(ImgIdx imgIdx, ImgIdx edgeIdx, ImgIdx *queuedEdges, uint8_t *numQueuedEdges) const;
    void compute_dimg_hhpq(double *dimg, ImgIdx *dhist, const Pixel *img, double a);
    void compute_dimg_hhpq_par(double *dimg, ImgIdx *dhist, const Pixel *img, double a);
    Pixel abs_diff(Pixel p, Pixel q);
    uint8_t compute_incidedge_queue(Pixel d0, Pixel d1);
    void compute_dimg_par4(RankItem<double> *&rankitem, const Pixel *img, SortValue<double> *&vals);
    void compute_dimg_par4(RankItem<double> *&rankitem, const Pixel *img, SortValue<Pixel> *&vals);
    Pixel compute_dimg(double *dimg, const Pixel *img);
    Pixel compute_dimg1(Pixel *dimg, ImgIdx *dhist, const Pixel *img);
    void compute_dimg(ImgIdx &minidx, double &mindiff, Pixel *dimg, ImgIdx *dhist, const Pixel *img, double a);
    void compute_dimg(Pixel *dimg, ImgIdx *dhist, const Pixel *img, double a);
    Pixel compute_dimg(Pixel *dimg, ImgIdx *dhist, const Pixel *img);
    double compute_dimg(double *dimg, ImgIdx *dhist, const Pixel *img);
    void set_isAvailable(uint8_t *isAvailable);
    void set_isAvailable(uint8_t *isAvailable, int npartitions_hor, int npartitions_ver);
    uint8_t is_available(uint8_t isAvailable, uint8_t iNeighbour) const;
    void set_field(uint8_t *arr, ImgIdx idx, uint8_t in);
    uint8_t get_field(uint8_t *arr, ImgIdx idx);
    void connectPix2Node(ImgIdx pidx, Pixel pix_val, ImgIdx iNode, Pixel level);
    void connectPix2Node(ImgIdx pidx, Pixel pix_val, ImgIdx iNode);
    void connectPix2Node0(ImgIdx pidx, Pixel pix_val, ImgIdx iNode, Pixel level);
    ImgIdx NewAlphaNode();
    ImgIdx NewAlphaNode(Pixel level, AlphaNode<Pixel> *pCopy);
    ImgIdx NewAlphaNode1(double level, AlphaNode<Pixel> *pCopy);
    ImgIdx NewAlphaNode(Pixel level);
    uint8_t is_visited(uint8_t *isVisited, ImgIdx p);
    void visit(uint8_t *isVisited, ImgIdx p);
    ImgIdx TreeSizeEstimation(ImgIdx *dhist, int64_t numlevels, ImgIdx imgSize, ImgIdx nredges);
    ImgIdx TreeSizeEstimation(ImgIdx *dhist, int64_t numlevels, ImgIdx imgSize, ImgIdx nredges, double m);
    ImgIdx TreeSizeEstimation(ImgIdx *dhist, int64_t numlevels, ImgIdx imgSize, ImgIdx nredges, double m,
                              ImgIdx reserve);
    void remove_redundant_node(ImgIdx &prev_top, ImgIdx &stack_top);
    int get_bitdepth(uint64_t num);
    ImgIdx initialize_node(const Pixel *img, Pixel *dimg, Pixel maxpixval);
    void initialize_node1(const Pixel *img, RankItem<double> *rankitem, Pixel maxpixval);
    void initialize_node1(const Pixel *img, RankItem<double> *rankitem, Pixel maxpixval, int32_t *rank2rankitem);
    void initialize_node(const Pixel *img, RankItem<Pixel> *rankitem, Pixel maxpixval);
    void initialize_node_par(const Pixel *img, RankItem<Pixel> *rankitem, Pixel maxpixval);
    void initialize_node_par1(const Pixel *img, RankItem<double> *rankitem, Pixel maxpixval, int32_t *rank2rankitem);
    void init_hypergraph_nodes(Pixel *dimg);
    void init_hypergraph_nodes(ImgIdx *rank);
    void set_isAvailable_hypergraph(uint8_t *isAvailable);
    uint8_t push_neighbor(Trie<TrieIdx> *queue, uint8_t *isVisited, ImgIdx *rank, ImgIdx p);
    void set_isAvailable_par_hypergraph(uint8_t *isAvailable, int8_t npartition_x, int8_t npartition_y);
    void cumsum(ImgIdx *hist, ImgIdx size, ImgIdx &maxidx);
    void cumsum(ImgIdx *hist, ImgIdx size, uint32_t *histeqmap, int eqhistsize);
    uint8_t push_neighbor(HierarQueue *queue, uint8_t *isVisited, uint8_t *dimg, ImgIdx p);
    uint8_t push_neighbor(HierarQueue *queue, uint8_t *isVisited, uint16_t *dimg, ImgIdx p);
    uint8_t push_neighbor(HierarQueue *queue, uint8_t *isVisited, uint32_t *dimg, ImgIdx p);
    uint8_t push_neighbor(HierarQueue *queue, uint8_t *isVisited, uint64_t *dimg, ImgIdx p);
    void canonicalize(ImgIdx nidx);
    ImgIdx merge_subtrees(uint8_t *dimg, int64_t blksz_x, int64_t blksz_y, int16_t npartition_x, int16_t npartition_y,
                          ImgIdx *subtree_cur, int tse);
    ImgIdx merge_subtrees1(uint8_t *dimg, int64_t blksz_x, int64_t blksz_y, int16_t npartition_x, int16_t npartition_y,
                           ImgIdx *subtree_cur, int tse, ImgIdx *hypernode_level);
    int migrate_subtree(int blk, int numpartitions, ImgIdx &nidx, ImgIdx &nidx_lim, int &nidxblk, ImgIdx &blkts,
                        char *blkflooddone, ImgIdx *subtree_cur, ImgIdx *subtree_start, ImgIdx *subtree_nborderedges,
                        omp_lock_t *locks, int &numbusythr, int &numblkproc, int &outofmemory);
    ImgIdx parflood_node_alloc(ImgIdx *subtree_size, ImgIdx *subtree_start, ImgIdx *blkws, ImgIdx *blkhs,
                               int numpartitions, double sizemult);
    void set_isAvailable_par(uint8_t *isAvailable, int16_t npartition_x, int16_t npartition_y);
    ImgIdx find_root(ImgIdx p);
    ImgIdx find_root_in(ImgIdx p);
    void blockwise_tse(ImgIdx *subtree_size, ImgIdx *subtree_nborderedges, double *nrmsds, ImgIdx *dhist,
                       ImgIdx *subtree_max, ImgIdx *blkws, ImgIdx *blkhs, int8_t npartition_x, int8_t npartition_y,
                       ImgIdx numbins);
    void quantize_ranks_compute_histogram(uint8_t *qrank, ImgIdx *rank, const Pixel *img, ImgIdx *dhist, ImgIdx *blkws,
                                          ImgIdx *blkhs, ImgIdx *startpidx, int64_t binsize, ImgIdx numbins,
                                          int8_t npartition_x, int8_t npartition_y, ImgIdx *subtree_max);
    uint8_t pow_quantization(ImgIdx rank, uint64_t qint);
    void pow_quantize_ranks(uint8_t *qrank, ImgIdx *rank, int64_t dimgsize, int64_t qint);
    ImgIdx find_root(AlphaNode<Pixel> *pilottree, ImgIdx p, Pixel below_this_qlevel);
    ImgIdx descendroots(ImgIdx q, int64_t qlevel, AlphaNode<Pixel> *pilottree);
    void unionfind_refine_qlevel(int64_t qlevel, int64_t binsize, ImgIdx nredges, AlphaNode<Pixel> *pilottree,
                                 RankItem<double> *rankitem, int8_t *redundant_edge, int32_t *rank2rankitem);
    void compute_dhist_par(uint8_t *qrank, ImgIdx *dhist, ImgIdx *startpidx, int32_t numbins, int8_t npartition_x,
                           int8_t npartition_y, int64_t blksz_x, int64_t blksz_y, int64_t blksz_xn, int64_t blksz_yn);
    void compute_dhist_par_hypergraph(uint8_t *qrank, ImgIdx *dhist, ImgIdx *startpidx, int32_t numbins,
                                      int8_t npartition_x, int8_t npartition_y, int64_t blksz_x, int64_t blksz_y,
                                      int64_t blksz_xn, int64_t blksz_yn, ImgIdx *blkmaxpidx);
    void fix_subtreeidx(ImgIdx *subtreestart, ImgIdx *startpidx, ImgIdx *cursizes, int8_t npartition_x,
                        int8_t npartition_y, int numpartitions, int64_t blksz_x, int64_t blksz_y, int64_t blksz_xn,
                        int64_t blksz_yn);
    void merge_subtrees(uint8_t *qrank, ImgIdx *qindex, int64_t blksz_x, int64_t blksz_y, ImgIdx neighbor_offset,
                        ImgIdx shamt, ImgIdx npartition_x, ImgIdx npartition_y, int32_t numbins);
    void merge_subtrees(uint8_t *qrank, int64_t blksz_x, int64_t blksz_y, ImgIdx neighbor_offset, ImgIdx shamt,
                        ImgIdx npartition_x, ImgIdx npartition_y, int32_t numbins);
    void connect_pilotnode(AlphaNode<Pixel> *pilottree, ImgIdx nredges, ImgIdx imgSize);
    void set_qindex(ImgIdx *qindex, ImgIdx *dhist, int64_t numpartitions, int32_t numbins, ImgIdx npartition_x,
                    ImgIdx npartition_y, int64_t blksz_x, int64_t blksz_y, int64_t blksz_xn, int64_t blksz_yn);
    void set_qindex(ImgIdx *qindex, ImgIdx *dhist, int64_t numpartitions, int32_t numbins);
    void set_subtree_root(ImgIdx **subtreerootary, ImgIdx *strary, ImgIdx nonzero_nodeidx_start,
                          ImgIdx rootlevel_nodeidx_start);
    void find_redundant_nodes(uint8_t *is_redundant, ImgIdx *rank);
    void set_subblock_properties(ImgIdx *startpidx, ImgIdx *blkws, ImgIdx *blkhs, ImgIdx *blocksize,
                                 int8_t npartition_x, int8_t npartition_y, int64_t blksz_x, int64_t blksz_y,
                                 int64_t blksz_xn, int64_t blksz_yn);
    void memalloc_queues(HierarQueue ***queues, int64_t numpartitions, ImgIdx *blocksize, ImgIdx *subtree_max);
    void compute_dimg_and_rank2index(RankItem<double> *&rankitem, const Pixel *img, ImgIdx nredges,
                                     int32_t *rank2rankitem);
    void compute_difference_and_sort(RankItem<double> *&rankitem, const Pixel *img, ImgIdx nredges);
    void compute_difference_and_sort(ImgIdx *rank, RankItem<double> *&rankitem, const Pixel *img, ImgIdx nredges,
                                     int32_t *&rank2rankitem);
    ImgIdx NewAlphaNode(ImgIdx &size, ImgIdx &maxsize);
    ImgIdx NewAlphaNode(AlphaNode<Pixel> *tree, ImgIdx &size, ImgIdx &maxsize, Pixel level, AlphaNode<Pixel> *pCopy);
    void remove_redundant_node(AlphaNode<Pixel> *tree, ImgIdx &size, ImgIdx &prev_top, ImgIdx &stack_top);
    void connectPix2Node(AlphaNode<Pixel> *tree, ImgIdx pidx, Pixel pix_val, ImgIdx iNode, ImgIdx *pAry);
    ImgIdx find_root1(ImgIdx p, ImgIdx qlevel);
    ImgIdx get_level_root(ImgIdx p, ImgIdx nodeidx);
    ImgIdx get_level_root(ImgIdx p);
    ImgIdx get_level_root(ImgIdx p, Pixel alpha);
    ImgIdx get_level_root(ImgIdx p, AlphaNode<Pixel> *tree);
    void swap(ImgIdx &x, ImgIdx &y);
    void swap(AlphaNode<Pixel> **x, AlphaNode<Pixel> **y);
    ImgIdx get_nearest_common_ancestor(ImgIdx x, ImgIdx y);
    Pixel get_nearest_common_ancestor_level(ImgIdx x, ImgIdx y);
    Pixel connect(ImgIdx x, ImgIdx y, ImgIdx newidx, Pixel alpha);
    Pixel connect(ImgIdx x, ImgIdx y, Pixel alpha, ImgIdx newidx);
    void canonicalize();
    void merge_subtrees(ImgIdx *rank, RankItem<Pixel> *rankitem, int64_t blksz_x, int64_t blksz_y,
                        ImgIdx neighbor_offset, ImgIdx shamt, ImgIdx npartition_x, ImgIdx npartition_y);
    void set_subimgsizes(ImgIdx **subimgsizes, int8_t npartition_x, int8_t npartition_y, int64_t blksz,
                         int64_t blksz_lastcol, int64_t blksz_lastrow, int64_t blksz_last);
};
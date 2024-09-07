#pragma once
#include "AlphaTreeConfig.h"
#include "HeapQueue.h"
#include "HierarQueue.h"
#include "HybridQueue.h"
#include "LadderQueue.hpp"
#include "Trie.h"
#include "defines.h"
#include "radixsort_teeninga/sort/radix_sort_parallel.h"
#include "radixsort_teeninga/sort/sort_item.h"
#include "walltime.h"
#include <PixelDissimilarity.hpp>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <omp.h>

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

    void BuildAlphaTree(Pixel *img, int height_in, int width_in, int channel_in, std::string dMetric,
                        int connectivity_in, int algorithm, int numthreads, int tse, double fparam1 = 0.0,
                        double fparam2 = 0.0, int iparam1 = 0);

    void AlphaFilter(double *outimg, double alpha);
    void AreaFilter(double *outimg, double area);

    void printTree() const;
    void printGraph(_uint8 *isVisited, bool *isRedundant, Pixel *img) const;
    void printGraph(_uint8 *isVisited, _uint8 *edge, Pixel *img) const;
    void printParentAry() const;
    void printAll(_uint8 *isVisited, bool *isRedundant, Pixel *img) const;
    void printAll(_uint8 *isVisited, _uint8 *edge, Pixel *img) const;
    void printVisit(ImgIdx p, double q) const;

  private:
    PixelDissimilarity<Pixel> _pixelDissim;

    void Unionfind(Pixel *img);
    void FloodHierarQueueNoCache(Pixel *img, int tse);
    void FloodHeapQueueNoCache(Pixel *img);
    void FloodHeapQueueNaiveNoCache(Pixel *img);
    void FloodHeapQueue(Pixel *img);
    void FloodHierarQueue(Pixel *img);
    void FloodLadderQueue(Pixel *img, int thres = 64);
    void FloodHierarHeapQueueNoCache(Pixel *img, double a = 12.0, double r = 0.5, int listsize = 12);
    void FloodHierarHeapQueue(Pixel *img, double a = 12.0, double r = 0.5, int listsize = 12);
    void FloodHierarHeapQueuePar(Pixel *img, double a = 12.0, double r = 0.5, int listsize = 12);
    void FloodHierHeapQueueHisteq(Pixel *img, int listsize = 12, int a = 0);
    void Flood_Hierarqueue_par(Pixel *img, int numthreads);
    void FloodTrieHypergraph(Pixel *img);

    void sortAlphaNodes();
    void markRedundant(ImgIdx imgIdx, ImgIdx eIdx, _uint8 *edgeStatus, ImgIdx *queuedEdges, _uint8 *numQueuedEdges);
    void queueEdge(ImgIdx imgIdx, ImgIdx edgeIdx, ImgIdx *queuedEdges, _uint8 *numQueuedEdges);
    void compute_dimg_hhpq(double *dimg, ImgIdx *dhist, Pixel *img, double a);
    void compute_dimg_hhpq_par(double *dimg, ImgIdx *dhist, Pixel *img, double a);
    Pixel abs_diff(Pixel p, Pixel q);
    _uint8 compute_incidedge_queue(Pixel d0, Pixel d1);
    void compute_dimg_par4(RankItem<double> *&rankitem, Pixel *img, SortValue<double> *&vals);
    void compute_dimg_par4(RankItem<double> *&rankitem, Pixel *img, SortValue<Pixel> *&vals);
    Pixel compute_dimg(double *dimg, Pixel *img);
    Pixel compute_dimg1(Pixel *dimg, ImgIdx *dhist, Pixel *img);
    void compute_dimg(ImgIdx &minidx, double &mindiff, Pixel *dimg, ImgIdx *dhist, Pixel *img, double a);
    void compute_dimg(Pixel *dimg, ImgIdx *dhist, Pixel *img, double a);
    Pixel compute_dimg(Pixel *dimg, ImgIdx *dhist, Pixel *img);
    double compute_dimg(double *dimg, ImgIdx *dhist, Pixel *img);
    void set_isAvailable(_uint8 *isAvailable);
    void set_isAvailable(_uint8 *isAvailable, int npartitions_hor, int npartitions_ver);
    _uint8 is_available(_uint8 isAvailable, _uint8 iNeighbour);
    void set_field(_uint8 *arr, ImgIdx idx, _uint8 in);
    _uint8 get_field(_uint8 *arr, ImgIdx idx);
    void connectPix2Node(ImgIdx pidx, Pixel pix_val, ImgIdx iNode, Pixel level);
    void connectPix2Node(ImgIdx pidx, Pixel pix_val, ImgIdx iNode);
    void connectPix2Node0(ImgIdx pidx, Pixel pix_val, ImgIdx iNode, Pixel level);
    ImgIdx NewAlphaNode();
    ImgIdx NewAlphaNode(Pixel level, AlphaNode<Pixel> *pCopy);
    ImgIdx NewAlphaNode1(double level, AlphaNode<Pixel> *pCopy);
    ImgIdx NewAlphaNode(Pixel level);
    _uint8 is_visited(_uint8 *isVisited, ImgIdx p);
    void visit(_uint8 *isVisited, ImgIdx p);
    ImgIdx TreeSizeEstimation(ImgIdx *dhist, _int64 numlevels, ImgIdx imgSize, ImgIdx nredges);
    ImgIdx TreeSizeEstimation(ImgIdx *dhist, _int64 numlevels, ImgIdx imgSize, ImgIdx nredges, double m);
    ImgIdx TreeSizeEstimation(ImgIdx *dhist, _int64 numlevels, ImgIdx imgSize, ImgIdx nredges, double m,
                              ImgIdx reserve);
    void remove_redundant_node(ImgIdx &prev_top, ImgIdx &stack_top);
    int get_bitdepth(_uint64 num);
    ImgIdx initialize_node(Pixel *img, Pixel *dimg, Pixel maxpixval);
    void initialize_node1(Pixel *img, RankItem<double> *rankitem, Pixel maxpixval);
    void initialize_node1(Pixel *img, RankItem<double> *rankitem, Pixel maxpixval, _int32 *rank2rankitem);
    void initialize_node(Pixel *img, RankItem<Pixel> *rankitem, Pixel maxpixval);
    void initialize_node_par(Pixel *img, RankItem<Pixel> *rankitem, Pixel maxpixval);
    void initialize_node_par1(Pixel *img, RankItem<double> *rankitem, Pixel maxpixval, _int32 *rank2rankitem);
    void init_hypergraph_nodes(Pixel *dimg);
    void init_hypergraph_nodes(ImgIdx *rank);
    void set_isAvailable_hypergraph(_uint8 *isAvailable);
    _uint8 push_neighbor(Trie<TrieIdx> *queue, _uint8 *isVisited, ImgIdx *rank, ImgIdx p);
    void set_isAvailable_par_hypergraph(_uint8 *isAvailable, _int8 npartition_x, _int8 npartition_y);
    void cumsum(ImgIdx *hist, ImgIdx size, ImgIdx &maxidx);
    void cumsum(ImgIdx *hist, ImgIdx size, _uint32 *histeqmap, int eqhistsize);
    _uint8 push_neighbor(HierarQueue *queue, _uint8 *isVisited, _uint8 *dimg, ImgIdx p);
    _uint8 push_neighbor(HierarQueue *queue, _uint8 *isVisited, _uint16 *dimg, ImgIdx p);
    _uint8 push_neighbor(HierarQueue *queue, _uint8 *isVisited, _uint32 *dimg, ImgIdx p);
    _uint8 push_neighbor(HierarQueue *queue, _uint8 *isVisited, _uint64 *dimg, ImgIdx p);
    void FloodHierarQueueHypergraph(Pixel *img);
    void canonicalize(ImgIdx nidx);
    ImgIdx merge_subtrees(Pixel *dimg, _int64 blksz_x, _int64 blksz_y, ImgIdx npartition_x, ImgIdx npartition_y,
                          ImgIdx *subtree_cur, int tse, ImgIdx *nrbnode = NULL);
    ImgIdx merge_subtrees(Pixel *dimg, _int64 blksz_x, _int64 blksz_y, ImgIdx npartition_x, ImgIdx npartition_y,
                          ImgIdx *subtree_cur, ImgIdx *subtree_start = NULL, ImgIdx *blkhs = NULL,
                          ImgIdx *blkws = NULL);
    ImgIdx merge_subtrees(_uint8 *dimg, _int64 blksz_x, _int64 blksz_y, _int16 npartition_x, _int16 npartition_y,
                          ImgIdx *subtree_cur, int tse);
    ImgIdx merge_subtrees1(_uint8 *dimg, _int64 blksz_x, _int64 blksz_y, _int16 npartition_x, _int16 npartition_y,
                           ImgIdx *subtree_cur, int tse, ImgIdx *hypernode_level);
    int migrate_subtree(int blk, int numpartitions, ImgIdx &nidx, ImgIdx &nidx_lim, int &nidxblk, ImgIdx &blkts,
                        char *blkflooddone, ImgIdx *subtree_cur, ImgIdx *subtree_start, ImgIdx *subtree_nborderedges,
                        omp_lock_t *locks, int &numbusythr, int &numblkproc, int &outofmemory);
    ImgIdx parflood_node_alloc(ImgIdx *subtree_size, ImgIdx *subtree_start, ImgIdx *blkws, ImgIdx *blkhs,
                               int numpartitions, double sizemult);
    void set_isAvailable_par(_uint8 *isAvailable, _int16 npartition_x, _int16 npartition_y);
    ImgIdx find_root(ImgIdx p);
    ImgIdx find_root_in(ImgIdx p);
    void blockwise_tse(ImgIdx *subtree_size, ImgIdx *subtree_nborderedges, double *nrmsds, ImgIdx *dhist,
                       ImgIdx *subtree_max, ImgIdx *blkws, ImgIdx *blkhs, _int8 npartition_x, _int8 npartition_y,
                       ImgIdx numbins);
    void quantize_ranks_compute_histogram(_uint8 *qrank, ImgIdx *rank, Pixel *img, ImgIdx *dhist, ImgIdx *blkws,
                                          ImgIdx *blkhs, ImgIdx *startpidx, _int64 binsize, ImgIdx numbins,
                                          _int8 npartition_x, _int8 npartition_y, ImgIdx *subtree_max);
    _uint8 pow_quantization(ImgIdx rank, _uint64 qint);
    void pow_quantize_ranks(_uint8 *qrank, ImgIdx *rank, _int64 dimgsize, _int64 qint);
    ImgIdx find_root(AlphaNode<Pixel> *pilottree, ImgIdx p, Pixel below_this_qlevel);
    ImgIdx descendroots(ImgIdx q, _int64 qlevel, AlphaNode<Pixel> *pilottree);
    void unionfind_refine_qlevel(_int64 qlevel, _int64 binsize, ImgIdx nredges, AlphaNode<Pixel> *pilottree,
                                 RankItem<double> *rankitem, _int8 *redundant_edge, _int32 *rank2rankitem);
    void compute_dhist_par(_uint8 *qrank, ImgIdx *dhist, ImgIdx *startpidx, _int32 numbins, _int8 npartition_x,
                           _int8 npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn);
    void compute_dhist_par_hypergraph(_uint8 *qrank, ImgIdx *dhist, ImgIdx *startpidx, _int32 numbins,
                                      _int8 npartition_x, _int8 npartition_y, _int64 blksz_x, _int64 blksz_y,
                                      _int64 blksz_xn, _int64 blksz_yn, ImgIdx *blkmaxpidx);
    void fix_subtreeidx(ImgIdx *subtreestart, ImgIdx *startpidx, ImgIdx *cursizes, _int8 npartition_x,
                        _int8 npartition_y, int numpartitions, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn,
                        _int64 blksz_yn);
    void merge_subtrees(_uint8 *qrank, ImgIdx *qindex, _int64 blksz_x, _int64 blksz_y, ImgIdx neighbor_offset,
                        ImgIdx shamt, ImgIdx npartition_x, ImgIdx npartition_y, _int32 numbins);
    void merge_subtrees(_uint8 *qrank, _int64 blksz_x, _int64 blksz_y, ImgIdx neighbor_offset, ImgIdx shamt,
                        ImgIdx npartition_x, ImgIdx npartition_y, _int32 numbins);
    void connect_pilotnode(AlphaNode<Pixel> *pilottree, ImgIdx nredges, ImgIdx imgSize);
    void set_qindex(ImgIdx *qindex, ImgIdx *dhist, _int64 numpartitions, _int32 numbins, ImgIdx npartition_x,
                    ImgIdx npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn);
    void set_qindex(ImgIdx *qindex, ImgIdx *dhist, _int64 numpartitions, _int32 numbins);
    void set_subtree_root(ImgIdx **subtreerootary, ImgIdx *strary, ImgIdx nonzero_nodeidx_start,
                          ImgIdx rootlevel_nodeidx_start);
    void find_redundant_nodes(_uint8 *is_redundant, ImgIdx *rank);
    void set_subblock_properties(ImgIdx *startpidx, ImgIdx *blkws, ImgIdx *blkhs, ImgIdx *blocksize, _int8 npartition_x,
                                 _int8 npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn);
    void memalloc_queues(HierarQueue ***queues, _int64 numpartitions, ImgIdx *blocksize, ImgIdx *subtree_max);
    void compute_dimg_and_rank2index(RankItem<double> *&rankitem, Pixel *img, ImgIdx nredges, _int32 *rank2rankitem);
    void compute_difference_and_sort(RankItem<double> *&rankitem, Pixel *img, ImgIdx nredges);
    void compute_difference_and_sort(ImgIdx *rank, RankItem<double> *&rankitem, Pixel *img, ImgIdx nredges,
                                     _int32 *&rank2rankitem);
    void HybridParallel(Pixel *img, int numthreads);
    ImgIdx NewAlphaNode(ImgIdx &size, ImgIdx &maxsize);
    ImgIdx NewAlphaNode(AlphaNode<Pixel> *tree, ImgIdx &size, ImgIdx &maxsize, Pixel level, AlphaNode<Pixel> *pCopy);
    void remove_redundant_node(AlphaNode<Pixel> *tree, ImgIdx &size, ImgIdx &prev_top, ImgIdx &stack_top);
    void connectPix2Node(AlphaNode<Pixel> *tree, ImgIdx pidx, Pixel pix_val, ImgIdx iNode, ImgIdx *pAry);
    ImgIdx find_root1(ImgIdx p, ImgIdx qlevel);
    void FloodTrieNoCache(Pixel *img);
    void FloodTrie(Pixel *img);
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
    void merge_subtrees(ImgIdx *rank, RankItem<Pixel> *rankitem, _int64 blksz_x, _int64 blksz_y, ImgIdx neighbor_offset,
                        ImgIdx shamt, ImgIdx npartition_x, ImgIdx npartition_y);
    void set_subimgsizes(ImgIdx **subimgsizes, _int8 npartition_x, _int8 npartition_y, _int64 blksz,
                         _int64 blksz_lastcol, _int64 blksz_lastrow, _int64 blksz_last);
};
#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <cfloat>
#include <cassert>
#include "defines.h"
#include "HierarQueue.h"
#include "Trie.h"
#include "HybridQueue.h"
#include "HeapQueue.h"
#include "radixsort_teeninga/sort/radix_sort_parallel.h"
#include "radixsort_teeninga/sort/sort_item.h"
#include "walltime.h"
#include "LadderQueue.hpp"

using namespace pmt; // TODO: remove

#define CHKRNG(var,a,b) ( (var >= a) && (var < b) )
#define QUANTIZE_RANK(rank, binsize) (rank) / (binsize)

template<class Pixel>
class AlphaNode
{
public:
	Imgidx area;
	double alpha;
	double sumPix;
	Pixel minPix;
	Pixel maxPix;
	Imgidx parentidx;
	Imgidx rootidx;

	void set(Imgidx area_in, double level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in);
	void add(AlphaNode* q);
	void add(Pixel pix_val);
	void copy(AlphaNode* q);
	void connect_to_parent(AlphaNode* pPar, Imgidx iPar);
	void print(AlphaNode* node);
	void print(AlphaNode* node, int heading);
};

template<class Pixel>
class RankItem
{
public:
	Pixel alpha;
	Imgidx dimgidx;

	void operator=(const RankItem& item);
	Imgidx get_pidx0(Imgidx connectivity = 4);
	Imgidx get_pidx1(Imgidx width, Imgidx connectivity = 4);
};

template<class Pixel>
class AlphaTree
{
public:
	Imgidx maxSize;
	Imgidx curSize;
	Imgidx height, width, channel, connectivity;
	AlphaNode<Pixel> *node, *node_in;
	Imgidx num_node, num_node_in;
	Imgidx rootidx;
	Imgidx* parentAry;
	double nrmsd;

	AlphaTree() : maxSize(0), curSize(0), node(0), parentAry(0) {}
	~AlphaTree();

	inline void clear() { Free(node); Free(parentAry); node = NULL; parentAry = NULL; curSize = 0; }

	void BuildAlphaTree(Pixel *img, int height_in, int width_in, int channel_in, int connectivity_in, const char* algorithm, int numthreads, int tse, double fparam1 = 0.0, double fparam2 = 0.0, int iparam1 = 0);

    void AlphaFilter(double *outimg, double alpha);

    void AreaFilter(double *outimg, double area);

    void print_tree();

  	static void getAlgorithmName(char* out, int algorithmCode);

	static void getAlgorithmDescriptionFromCode(char* out, int algorithmCode);

  	static void getAlgorithmDescription(char* out, const char* algorithmName);

	static int constexpr NUM_ALGORITHMS = 16;

private:

	///< Root or subtree roots are denoted by having NULLINDEX as parent
	static Imgidx constexpr NULLINDEX = -1;

	static int constexpr FLOOD_HIERARHEAPQUEUE_CACHE = 0;
	static int constexpr FLOOD_HIERARHEAPQUEUE_CACHE_HISTEQ = 1;
	static int constexpr FLOOD_HIERARHEAPQUEUE = 2;
	static int constexpr PILOT_RANK = 3;
	static int constexpr UNIONFIND = 4;
	static int constexpr FLOOD_HIERARQUEUE_CACHE = 5;
	static int constexpr FLOOD_HIERARQUEUE = 8;
	static int constexpr FLOOD_HIERARQUEUE_PAR = 7;
	static int constexpr FLOOD_TRIE_CACHE = 9;
	static int constexpr FLOOD_TRIE = 10;
	static int constexpr FLOOD_HEAPQUEUE_CACHE = 11;
	static int constexpr FLOOD_HEAPQUEUE = 12;
	static int constexpr FLOOD_HEAPQUEUE_NAIVE = 13;
	static int constexpr FLOOD_HIERARQUEUE_HYPERGRAPH = 14;
	static int constexpr FLOOD_TRIE_HYPERGRAPH = 15;
	static int constexpr FLOOD_LADDERQUEUE = 16;

	static char constexpr UNIONFIND_COMMAND[] = "UnionFind";
	static char constexpr UNIONFIND_COMMAND_ALT[] = "UF";
	static char constexpr FLOOD_HIERARQUEUE_CACHE_COMMAND[] = "FloodHierQueue";
	static char constexpr FLOOD_HIERARQUEUE_CACHE_COMMAND_ALT[] = "FHIQ";
	static char constexpr FLOOD_HIERARQUEUE_COMMAND[] = "FloodHierQueueNoCache";
	static char constexpr FLOOD_TRIE_CACHE_COMMAND[] = "FloodTrie";
	static char constexpr FLOOD_TRIE_CACHE_COMMAND_ALT[] = "FT";
	static char constexpr FLOOD_TRIE_COMMAND[] = "FloodTrieNoCache";
	static char constexpr FLOOD_HEAPQUEUE_CACHE_COMMAND[] = "FloodHeapQueue";
	static char constexpr FLOOD_HEAPQUEUE_CACHE_COMMAND_ALT[] = "FHEQ";
	static char constexpr FLOOD_HEAPQUEUE_COMMAND[] = "FloodHeapQueueNoCache";
	static char constexpr FLOOD_HEAPQUEUE_NAIVE_COMMAND[] = "FloodHeapQueueNaive";
	static char constexpr FLOOD_HIERARQUEUE_HYPERGRAPH_COMMAND[] = "FloodHierQueueHyperGraph";
	static char constexpr FLOOD_TRIE_HYPERGRAPH_COMMAND[] = "FloodTrieHyperGraph";
	static char constexpr FLOOD_HIERARHEAPQUEUE_CACHE_COMMAND[] = "FloodHiearHeapQueue";
	static char constexpr FLOOD_HIERARHEAPQUEUE_CACHE_COMMAND_ALT[] = "FHHQ";
	static char constexpr FLOOD_HIERARHEAPQUEUE_COMMAND[] = "FloodHiearHeapQueueNoCache";
	static char constexpr FLOOD_HIERARHEAPQUEUE_CACHE_HISTEQ_COMMAND[] = "FloodHiearHeapQueueHistogramEqualization";
	static char constexpr FLOOD_HIERARQUEUE_PAR_COMMAND[] = "FloodHierQueueParallel";
	static char constexpr FLOOD_HIERARQUEUE_PAR_COMMAND_ALT[] = "FHIQP";
	static char constexpr PILOT_RANK_COMMAND[] = "HybridAlgorithm";
	static char constexpr PILOT_RANK_COMMAND_ALT[] = "HA";
	static char constexpr FLOOD_LADDERQUEUE_COMMAND[] = "FloodLadderQueue";

	///< Tree size estimation parameters
	// Minimum number of pixels for TSE to be used. TSE does not work well with very small images.
	static Imgidx constexpr TSE_MINSIZE = 10000;
	static double constexpr A = 1.3901;
	static double constexpr SIGMA = -2.1989;
	static double constexpr B = -0.1906;
	static double constexpr M = 0.05;

	static int parseAlgorithmCode(const char* algorithmName);
    Pixel abs_diff(Pixel p, Pixel q);
    _uint8 compute_incidedge_queue(Pixel d0, Pixel d1);
    void compute_dimg_par4(RankItem<double> *&rankitem, Pixel *img, SortValue<double> *&vals);
    void compute_dimg_par4(RankItem<double> *&rankitem, Pixel *img, SortValue<Pixel> *&vals);
    Pixel compute_dimg(double *dimg, Pixel *img);
    Pixel compute_dimg1(Pixel *dimg, Imgidx *dhist, Pixel *img);
    void compute_dimg(Imgidx &minidx, double &mindiff, Pixel *dimg, Imgidx *dhist, Pixel *img, double a);
    void compute_dimg(Pixel *dimg, Imgidx *dhist, Pixel *img, double a);
    Pixel compute_dimg(Pixel *dimg, Imgidx *dhist, Pixel *img);
    double compute_dimg(double *dimg, Imgidx *dhist, Pixel *img);
    void set_isAvailable(_uint8 *isAvailable);
    void set_isAvailable(_uint8 *isAvailable, int npartitions_hor, int npartitions_ver);
    _uint8 is_available(_uint8 isAvailable, _uint8 iNeighbour);
    void set_field(_uint8 *arr, Imgidx idx, _uint8 in);
    _uint8 get_field(_uint8 *arr, Imgidx idx);
    void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx iNode, Pixel level);
    void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx iNode);
    void connectPix2Node0(Imgidx pidx, Pixel pix_val, Imgidx iNode, Pixel level);
    Imgidx NewAlphaNode();
    Imgidx NewAlphaNode(Pixel level, AlphaNode<Pixel> *pCopy);
    Imgidx NewAlphaNode1(double level, AlphaNode<Pixel> *pCopy);
    Imgidx NewAlphaNode(Pixel level);
    _uint8 is_visited(_uint8 *isVisited, Imgidx p);
    void visit(_uint8 *isVisited, Imgidx p);
    Imgidx TreeSizeEstimation(Imgidx *dhist, _int64 numlevels, Imgidx imgsize, Imgidx nredges);
    Imgidx TreeSizeEstimation(Imgidx *dhist, _int64 numlevels, Imgidx imgsize, Imgidx nredges, double m);
    Imgidx TreeSizeEstimation(Imgidx *dhist, _int64 numlevels, Imgidx imgsize, Imgidx nredges, double m, Imgidx reserve);
    void remove_redundant_node(Imgidx &prev_top, Imgidx &stack_top);
    void Flood_HierarQueue(Pixel *img, int tse);
    void Flood_HeapQueue(Pixel *img);
    void Flood_HeapQueue_Naive(Pixel *img);
    void Flood_HeapQueue_Cache(Pixel *img);
    void Flood_HierarQueue_Cache(Pixel *img);
	void Flood_LadderQueue(Pixel *img, int thres = 64);
	int get_bitdepth(_uint64 num);
	void Flood_HierarHeapQueue(Pixel* img, double a = 12.0, double r = 0.5, int listsize = 12);
	void Flood_HierarHeapQueue_Cache(Pixel* img, double a = 12.0, double r = 0.5, int listsize = 12);
	void Flood_HierarHeapQueue_Cache_histeq(Pixel* img, int listsize = 12, int a = 0);
	Imgidx initialize_node(Pixel *img, Pixel *dimg, Pixel maxpixval);
	void initialize_node1(Pixel *img, RankItem<double> *rankitem, Pixel maxpixval);
	void initialize_node1(Pixel *img, RankItem<double> *rankitem, Pixel maxpixval, _int32* rank2rankitem);
	void initialize_node(Pixel *img, RankItem<Pixel> *rankitem, Pixel maxpixval);
	void initialize_node_par(Pixel *img, RankItem<Pixel> *rankitem, Pixel maxpixval);
	void initialize_node_par1(Pixel *img, RankItem<double> *rankitem, Pixel maxpixval, _int32* rank2rankitem);
	void init_hypergraph_nodes(Pixel *dimg);
	void init_hypergraph_nodes(Imgidx *rank);
	void set_isAvailable_hypergraph(_uint8 *isAvailable);
	_uint8 push_neighbor(Trie<trieidx> *queue, _uint8 *isVisited, Imgidx *rank, Imgidx p);
	void Flood_Trie_Hypergraph(Pixel *img);
	void set_isAvailable_par_hypergraph(_uint8* isAvailable, _int8 npartition_x, _int8 npartition_y);
	void cumsum(Imgidx *hist, Imgidx size, Imgidx &maxidx);
	void cumsum(Imgidx *hist, Imgidx size, _uint32 *histeqmap, int eqhistsize);
	_uint8 push_neighbor(HierarQueue *queue, _uint8 *isVisited, _uint8 *dimg, Imgidx p);
	_uint8 push_neighbor(HierarQueue *queue, _uint8 *isVisited, _uint16 *dimg, Imgidx p);
	_uint8 push_neighbor(HierarQueue *queue, _uint8 *isVisited, _uint32 *dimg, Imgidx p);
	_uint8 push_neighbor(HierarQueue *queue, _uint8 *isVisited, _uint64 *dimg, Imgidx p);
	void Flood_HierarQueue_Hypergraph(Pixel* img);
	void canonicalize(Imgidx nidx);
	Imgidx merge_subtrees(Pixel *dimg, _int64 blksz_x, _int64 blksz_y, Imgidx npartition_x, Imgidx npartition_y, Imgidx* subtree_cur, int tse, Imgidx *nrbnode = NULL);
	Imgidx merge_subtrees(Pixel *dimg, _int64 blksz_x, _int64 blksz_y, Imgidx npartition_x, Imgidx npartition_y, Imgidx* subtree_cur, Imgidx* subtree_start = NULL, Imgidx *blkhs = NULL, Imgidx *blkws = NULL);
	Imgidx merge_subtrees(_uint8 *dimg, _int64 blksz_x, _int64 blksz_y, _int16 npartition_x, _int16 npartition_y, Imgidx* subtree_cur, int tse);
	Imgidx merge_subtrees1(_uint8 *dimg, _int64 blksz_x, _int64 blksz_y, _int16 npartition_x, _int16 npartition_y, Imgidx* subtree_cur, int tse, Imgidx* hypernode_level);
	int migrate_subtree(int blk, int numpartitions, Imgidx & nidx, Imgidx & nidx_lim, int & nidxblk, Imgidx & blkts, char *blkflooddone, Imgidx *subtree_cur, Imgidx *subtree_start, Imgidx *subtree_nborderedges, omp_lock_t *locks, int &numbusythr, int &numblkproc, int &outofmemory);
	Imgidx parflood_node_alloc(Imgidx *subtree_size, Imgidx *subtree_start, Imgidx *blkws, Imgidx *blkhs, int numpartitions, double sizemult);
	void set_isAvailable_par(_uint8* isAvailable, _int16 npartition_x, _int16 npartition_y);
	void Flood_Hierarqueue_par(Pixel *img, int numthreads);
	Imgidx find_root(Imgidx p);
	Imgidx find_root_in(Imgidx p);
	void Unionfind(Pixel* img);
	void blockwise_tse(Imgidx *subtree_size, Imgidx *subtree_nborderedges, double *nrmsds, Imgidx *dhist, Imgidx *subtree_max, Imgidx *blkws, Imgidx *blkhs, _int8 npartition_x, _int8 npartition_y, Imgidx numbins);
	void quantize_ranks_compute_histogram(_uint8 *qrank, Imgidx* rank, Pixel* img, Imgidx *dhist, Imgidx *blkws, Imgidx *blkhs, Imgidx *startpidx, _int64 binsize, Imgidx numbins, _int8 npartition_x, _int8 npartition_y, Imgidx* subtree_max);
	_uint8 pow_quantization(Imgidx rank, _uint64 qint);
	void pow_quantize_ranks(_uint8 *qrank, Imgidx *rank, _int64 dimgsize, _int64 qint);
	Imgidx find_root(AlphaNode<Pixel> *pilottree, Imgidx p, Pixel below_this_qlevel);
	Imgidx descendroots(Imgidx q, _int64 qlevel, AlphaNode<Pixel> *pilottree);
	void unionfind_refine_qlevel(_int64 qlevel, _int64 binsize, Imgidx nredges, AlphaNode<Pixel>* pilottree, RankItem<double>* rankitem, _int8* redundant_edge, _int32* rank2rankitem);
	void compute_dhist_par(_uint8 *qrank, Imgidx *dhist, Imgidx *startpidx, _int32 numbins, _int8 npartition_x, _int8 npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn);
	void compute_dhist_par_hypergraph(_uint8 *qrank, Imgidx *dhist, Imgidx *startpidx, _int32 numbins, _int8 npartition_x, _int8 npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn, Imgidx *blkmaxpidx);
	void fix_subtreeidx(Imgidx *subtreestart, Imgidx *startpidx, Imgidx *cursizes, _int8 npartition_x, _int8 npartition_y, int numpartitions, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn);
	void merge_subtrees(_uint8 *qrank, Imgidx *qindex, _int64 blksz_x, _int64 blksz_y, Imgidx neighbor_offset, Imgidx shamt, Imgidx npartition_x, Imgidx npartition_y, _int32 numbins);
	void merge_subtrees(_uint8 *qrank, _int64 blksz_x, _int64 blksz_y, Imgidx neighbor_offset, Imgidx shamt, Imgidx npartition_x, Imgidx npartition_y, _int32 numbins);
	void connect_pilotnode(AlphaNode<Pixel> *pilottree, Imgidx nredges, Imgidx imgsize);
	void set_qindex(Imgidx *qindex, Imgidx *dhist, _int64 numpartitions, _int32 numbins, Imgidx npartition_x, Imgidx npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn);
	void set_qindex(Imgidx *qindex, Imgidx *dhist, _int64 numpartitions, _int32 numbins);
	void set_subtree_root(Imgidx** subtreerootary, Imgidx *strary, Imgidx nonzero_nodeidx_start, Imgidx rootlevel_nodeidx_start);
	void find_redundant_nodes(_uint8* is_redundant, Imgidx *rank);
	void set_subblock_properties(Imgidx* startpidx, Imgidx *blkws, Imgidx *blkhs, Imgidx *blocksize, _int8 npartition_x, _int8 npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn);
	void memalloc_queues(HierarQueue*** queues, _int64 numpartitions, Imgidx* blocksize, Imgidx* subtree_max);
	void compute_dimg_and_rank2index(RankItem<double>*& rankitem, Pixel* img, Imgidx nredges, _int32* rank2rankitem);
	void compute_difference_and_sort(RankItem<double>*& rankitem, Pixel* img, Imgidx nredges);
	void print_all_trees(AlphaNode<Pixel>* pilottree);
	void compute_difference_and_sort(Imgidx* rank, RankItem<double>*& rankitem, Pixel* img, Imgidx nredges, _int32*& rank2rankitem);
	void Pilot_Rank(Pixel* img, int numthreads);
	Imgidx NewAlphaNode(Imgidx &size, Imgidx &maxsize);
	Imgidx NewAlphaNode(AlphaNode<Pixel> *tree, Imgidx &size, Imgidx &maxsize, Pixel level, AlphaNode<Pixel> *pCopy);
	void remove_redundant_node(AlphaNode<Pixel> *tree, Imgidx &size, Imgidx& prev_top, Imgidx& stack_top);
	void connectPix2Node(AlphaNode<Pixel> *tree, Imgidx pidx, Pixel pix_val, Imgidx iNode, Imgidx *pAry);
	Imgidx find_root1(Imgidx p, Imgidx qlevel);
	void Flood_Trie(Pixel* img);
	void Flood_Trie_Cache(Pixel* img);
	Imgidx get_level_root(Imgidx p, Imgidx nodeidx);
	Imgidx get_level_root(Imgidx p);
	Imgidx get_level_root(Imgidx p, Pixel alpha);
	Imgidx get_level_root(Imgidx p, AlphaNode<Pixel> *tree);
	void swap(Imgidx& x, Imgidx& y);
	void swap(AlphaNode<Pixel>** x, AlphaNode<Pixel>** y);
	Imgidx get_nearest_common_ancestor(Imgidx x, Imgidx y);
	Pixel get_nearest_common_ancestor_level(Imgidx x, Imgidx y);
	Pixel connect(Imgidx x, Imgidx y, Imgidx newidx, Pixel alpha);
	Pixel connect(Imgidx x, Imgidx y, Pixel alpha, Imgidx newidx);
	void canonicalize();
	void merge_subtrees(Imgidx *rank, RankItem<Pixel>* rankitem, _int64 blksz_x, _int64 blksz_y, Imgidx neighbor_offset, Imgidx shamt, Imgidx npartition_x, Imgidx npartition_y);
	void set_subimgsizes(Imgidx** subimgsizes, _int8 npartition_x, _int8 npartition_y, _int64 blksz, _int64 blksz_lastcol, _int64 blksz_lastrow, _int64 blksz_last);

};
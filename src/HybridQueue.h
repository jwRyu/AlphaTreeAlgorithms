#pragma once
#include "allocator.h"
#include "defines.h"
#include "Trie.h"
#include "HierarQueue.h"
#include "HeapQueue.h"
#include "walltime.h"//tmp

#include <iostream>
#include <cfloat>
#include <cmath>

#define LISTSIZE_DEFAULT 12
#define HEAPSIZE_DEFAULT 128

#define PROFILE			0

template<class Pixel>
class HierarHeapQueue_HEQ
{
	HQentry<Pixel> *list;
	HeapQueue_naive_quad<Pixel> **hqueue;
	HQentry<Pixel> **storage;
	Imgidx *storage_cursize;
	Imgidx *qsizes;
	_uint32 *histeqmap;

	Imgidx thr_hqueue, curthr, numlevels;
	double a;
	Imgidx queue_minlev;

	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

public:
	void initHQ(Imgidx *dhist, _uint32* histeqmap_in, Imgidx numlevels_in, Imgidx size, double a_in, int listsize);
	HierarHeapQueue_HEQ(Imgidx *dhist, _uint32 *histeqmap_in, Imgidx numlevels_in, double a_in, Imgidx size);
	HierarHeapQueue_HEQ(Imgidx *dhist, _uint32 *histeqmap_in, Imgidx numlevels_in, Imgidx size, double a_in, int listsize);
	~HierarHeapQueue_HEQ();

    void push_1stitem(Imgidx idx, Pixel alpha);
    void end_pushes(_uint8 *isVisited);

    inline void start_pushes() { emptytop = 1; }
	inline Pixel get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline Pixel top_alpha() { return list[0].alpha; }
	
	
	void push(Imgidx idx, Pixel alpha);
	void push_queue(Imgidx idx, Pixel alpha);
	Imgidx pop(_uint8 *isVisited);
    int check_queue_level(_uint8 *isVisited);
    void pop_queue(_uint8 *isVisited);
};

template<class Pixel>
class HierarHeapQueue
{
	HeapQueue_naive_quad<Pixel> **hqueue;
	HQentry<Pixel> **storage;
	Imgidx *storage_cursize;
	Imgidx *qsizes;
	Imgidx thr_hqueue, curthr, numlevels;
	double a;
	Imgidx queue_minlev;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;

public:

	HierarHeapQueue(Imgidx *dhist, Imgidx numlevels_in, Imgidx size, double a_in, int listsize, int connectivity, double r);
	~HierarHeapQueue();
	inline Imgidx top() { return hqueue[queue_minlev]->top(); }
	inline Pixel top_alpha() { return hqueue[queue_minlev]->top_alpha(); }
	inline void push_1stitem(Imgidx idx, Pixel alpha) {push_queue(idx, alpha);}
	inline void end_pushes(_uint8* isVisited) {pop(isVisited);}
	inline void push(Imgidx idx, Pixel alpha) {push_queue(idx, alpha);}
	void push_queue(Imgidx idx, Pixel alpha);

	Imgidx pop(_uint8 *isVisited);
	int check_queue_level(_uint8* isVisited);
	void pop_queue(_uint8* isVisited);
};

template<class Pixel>
class HierarHeapQueue_cache
{
	HQentry<Pixel> *list;
	HeapQueue_naive_quad<Pixel> **hqueue;
	HQentry<Pixel> **storage;
	Imgidx *storage_cursize;
	Imgidx *qsizes;

	Imgidx thr_hqueue, curthr, numlevels;
	double a;
	Imgidx queue_minlev;

	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

	Imgidx totalsize;

public:
#if PROFILE
	double t0;
	double tconv;
	double tcache;
	double tqueue;

	Imgidx num_cache;
	Imgidx num_cache_ovfl;
	Imgidx num_hq;
	Imgidx num_store;
	Imgidx num_conv;
#endif

	void initHQ(Imgidx *dhist, Imgidx numlevels_in, Imgidx size, double a_in, int listsize, int connectivity, double r);
	HierarHeapQueue_cache(Imgidx *dhist, Imgidx numlevels_in, double a_in, Imgidx size, Imgidx connectivity = 4, double r = 0.2);
	HierarHeapQueue_cache(Imgidx *dhist, Imgidx numlevels_in, Imgidx size, double a_in, int listsize, Imgidx connectivity = 4, double r = 0.2);
	~HierarHeapQueue_cache();
	inline void start_pushes() {emptytop = 1;}
	inline Pixel get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline Pixel top_alpha() { return list[0].alpha; }
	void push_1stitem(Imgidx idx, Pixel alpha);
	void end_pushes(_uint8 *isVisited);
	void push(Imgidx idx, Pixel alpha);
	void push_queue(Imgidx idx, Pixel alpha);
	Imgidx pop(_uint8 *isVisited);
	int check_queue_level(_uint8 *isVisited);
	void pop_queue(_uint8 *isVisited);
};

template<class Pixel>
class Cache_Heapqueue
{
	HQentry<Pixel> *list;
	HeapQueue_naive<Pixel> *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

	void initHQ(Imgidx size, size_t listsize);
	
public:
	double qtime;//tmp

	Cache_Heapqueue(Imgidx size);
	Cache_Heapqueue(Imgidx size, size_t listsize);
	~Cache_Heapqueue();

	inline void start_pushes() { emptytop = 1; }
	inline Imgidx get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline Pixel top_alpha() { return list[0].alpha; }
	inline void end_pushes() { if(emptytop) pop(); }
	inline void push_queue(Imgidx idx, Pixel alpha) { hqueue->push(idx, alpha);}
	inline void pop_queue() { hqueue->pop(); }
	void push_1stitem(Imgidx idx, Pixel alpha);
	void push(Imgidx idx, Pixel alpha);
	Imgidx pop();
};

template<class Pixel>
class Cache_Quad_Heapqueue
{
	//MinList1 *list, *list_end, *head, *tail;
	HQentry<Pixel> *list;
	HeapQueue_naive_quad<Pixel> *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

	void initHQ(Imgidx size, size_t listsize);

public:

	double qtime;//tmp

	Cache_Quad_Heapqueue(Imgidx size);
	Cache_Quad_Heapqueue(Imgidx size, size_t listsize);
	~Cache_Quad_Heapqueue();

	inline void start_pushes() { emptytop = 1; }
	inline Pixel get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline Pixel top_alpha() { return list[0].alpha; }
	inline void end_pushes() { if(emptytop) pop(); }
	inline void push_queue(Imgidx idx, Pixel alpha) { hqueue->push(idx, alpha); }
	inline void pop_queue() { hqueue->pop(); }

	void push_1stitem(Imgidx idx, Pixel alpha);
	void push(Imgidx idx, Pixel alpha);
	Imgidx pop();
};

template<class Pixel>
class CirCache_Hierqueue
{
	HQentry<_int32> *list;
	HierarQueue *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list, liststart, mask;
	Imgidx maxSize_queue, mask_field;
	int emptytop;
	void initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize);

public:
	double qtime;//tmp

	CirCache_Hierqueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels);
	CirCache_Hierqueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize);
	~CirCache_Hierqueue();

    inline void start_pushes() { emptytop = 1; }
	inline _int32 get_minlev() { return list[liststart].alpha; }
	inline Imgidx top() { return list[liststart].pidx; }
	inline _int32 top_alpha() { return list[liststart].alpha; }
	inline void end_pushes() { if(emptytop) pop(); }
	inline void push_queue(Imgidx idx, _int32 alpha) { hqueue->push(idx, alpha); }
	inline void pop_queue() { hqueue->pop(); hqueue->find_minlev(); }

	void push_1stitem(Imgidx idx, _int32 alpha);
    void push(Imgidx idx, _int32 alpha);
	Imgidx pop();
};

template<class Pixel>
class HierarQueueCache
{
	HQentry<_int32> *list;
	HierarQueue *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;
	void initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize = LISTSIZE_DEFAULT);
	
public:
	double qtime;//tmp

	HierarQueueCache(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels);
	HierarQueueCache(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize);
	~HierarQueueCache();

	inline void start_pushes() { emptytop = 1; }
	inline _int32 get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline _int32 top_alpha() { return list[0].alpha; }
	inline void end_pushes() { if(emptytop) pop(); }
	inline void push_queue(Imgidx idx, _int32 alpha) { hqueue->push(idx, alpha); }
	inline void pop_queue() { hqueue->pop();hqueue->find_minlev(); }
	
	void push_1stitem(Imgidx idx, _int32 alpha);
	void push(Imgidx idx, _int32 alpha);
	Imgidx pop();
};

template<class Pixel>
class Cache_Hierqueue_l1
{
	HQentry<_int32> *list;
	HQueue_l1idx *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

	void initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize);

public:

	double qtime;//tmp

	Cache_Hierqueue_l1(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels);
	Cache_Hierqueue_l1(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize);
	~Cache_Hierqueue_l1();

	inline void start_pushes() { emptytop = 1; }
	inline _int32 get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline _int32 top_alpha() { return list[0].alpha; }
	inline void end_pushes() { if(emptytop) pop(); }
	inline void push_queue(Imgidx idx, _int32 alpha) { hqueue->push(idx, alpha); }
	inline void pop_queue() { hqueue->pop();hqueue->find_minlev(); }

	void push_1stitem(Imgidx idx, _int32 alpha);
	void push(Imgidx idx, _int32 alpha);
	Imgidx pop();
};

template<class Pixel>
class Cache_Hierqueue_l2
{
	HQentry<_int32> *list;
	HQueue_l2idx *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

	void initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize);
public:

	double qtime;//tmp

	Cache_Hierqueue_l2(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels);
	Cache_Hierqueue_l2(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize);
	~Cache_Hierqueue_l2();
	
	inline void start_pushes() { emptytop = 1; }
	inline _int32 get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline _int32 top_alpha() { return list[0].alpha; }
	inline void end_pushes() { if(emptytop) pop(); }
	inline void push_queue(Imgidx idx, _int32 alpha) { hqueue->push(idx, alpha); }
	inline void pop_queue() { hqueue->pop();hqueue->find_minlev(); }

	void push_1stitem(Imgidx idx, _int32 alpha);
	void push(Imgidx idx, _int32 alpha);
	Imgidx pop();
};

class Trie_Cache
{
	Imgidx *list;
	Trie<_uint64> *trie;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;

	void initHQ(Imgidx size, size_t listsize);
public:
	Trie_Cache(Imgidx size);
	Trie_Cache(Imgidx size, size_t listsize);
	~Trie_Cache();

	inline Imgidx get_minlev() { return list[0]; }
	inline Imgidx top() { return list[0]; }
	inline void push_queue(Imgidx idx) { trie->push(idx); }
	inline void pop_queue() { trie->pop(); }

	void push(Imgidx idx);
	void pop();
};

class HybridQueue_HQueue_Rank
{
	Imgidx *list;
	HQueue_l1idx_rank *queue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;

	void initHQ(Imgidx size, size_t listsize);

public:
	HybridQueue_HQueue_Rank(Imgidx size);
	HybridQueue_HQueue_Rank(Imgidx size, size_t listsize);
	~HybridQueue_HQueue_Rank();
	
	inline Imgidx get_minlev() { return list[0]; }
	inline Imgidx top() { return list[0]; }
	inline void push_queue(Imgidx idx) { queue->push(idx); }
	void push(Imgidx idx);
	void pop();
	inline void pop_queue() { queue->pop(); }
};

class HybridQueue_HQueue_Rank1
{
	Imgidx *list;
	HQueue_l1idx_rank *queue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	_int16 l0, mask;
	
	void initHQ(Imgidx size, size_t listsize);

public:
	HybridQueue_HQueue_Rank1(Imgidx size);
	HybridQueue_HQueue_Rank1(Imgidx size, size_t listsize);
	~HybridQueue_HQueue_Rank1();

	inline Imgidx get_minlev() { return list[l0]; }
	inline Imgidx top() { return list[l0]; }
	inline void push_queue(Imgidx idx) { queue->push(idx); }
	inline void pop_queue() { queue->pop(); }

	void push(Imgidx idx);
	void pop();
};

class HybridQueue_HQueue
{
	Imgidx *list;
	_int64 *levels;
	HQueue_l1idx *queue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 minlevnotfixed;

	void initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize);
	
public:
	HybridQueue_HQueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels);
	HybridQueue_HQueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize);
	~HybridQueue_HQueue();

	inline Imgidx top(){return list[0];}
	inline _int64 get_minlev() { return levels[0]; }
	inline void find_minlev() {	queue->find_minlev(); }
	inline void push_queue(Imgidx idx, _int64 level) { if (queue->push(idx, level)) minlevnotfixed = 0; }
	inline void pop_queue() { minlevnotfixed = 1;queue->pop(); }

	void push(Imgidx idx, _int64 level);
	Imgidx pop();	
};

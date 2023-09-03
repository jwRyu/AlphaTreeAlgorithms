#pragma once

#include <cstdio>
#include "defines.h"
#include "allocator.h"

//Do not use beyond 20-bit image
class HierarQueue
{
public:
	Imgidx *queue;
	Imgidx *bottom, *cur;
	_int32 numlevel;
	_int64 qsize;
	_int64 min_level, max_level;

	void print();
	HierarQueue(_uint64 qsize_in, _int32 numlevels);
	HierarQueue(_uint64 qsize_in);
	void reset_queue();
	Imgidx set_queue(Imgidx *dhist);
	Imgidx set_queue(Imgidx *dhist, _int32 maxpix);

	HierarQueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels);
	HierarQueue(Imgidx *dhist, _int32 numlevels);
	HierarQueue(_int32 numlevels, Imgidx binsize);
	~HierarQueue();

	_int8 push(Imgidx pidx, _int64 level);
	void find_minlev();

	inline Imgidx pop() { return queue[bottom[min_level]++]; }
	inline Imgidx top() { return queue[bottom[min_level]]; }
	inline _int64 get_minlev() { return min_level; }
};

class HQueue_l1idx
{
public:
	Imgidx *queue;
	Imgidx *bottom, *cur;
	_uint64 *seeker;

	_int64 qsize, seekersize;
	_int32 min_level;
	HQueue_l1idx(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels);
	~HQueue_l1idx();

	int push(Imgidx pidx, _int32 level);
	Imgidx pop();
	void find_minlev();

	inline Imgidx top() { return queue[cur[min_level] - 1]; }
	inline _int32 get_minlev() { return min_level; }
};

class HQueue_l2idx
{
	Imgidx *queue;
	Imgidx *bottom, *cur;
	_uint64 *seeker,*seeker2;
public:
	_int64 qsize;
	_int64 min_level;
	HQueue_l2idx(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels);
	~HQueue_l2idx();

	void push(Imgidx pidx, _int64 level);
	Imgidx pop();
	void find_minlev();

	inline Imgidx top() { return queue[cur[min_level] - 1]; }
	inline _int64 get_minlev() { return min_level; }
};

class HQueue_l1idx_rank
{
	struct hqueue_word
	{
		_int64 qword[64];
		_int64 seeker;
	};

	hqueue_word *queue;
public:
	_int64 qsize, seekersize;
	_int64 min_level;
	HQueue_l1idx_rank(_int64 qsize_in);
	~HQueue_l1idx_rank();

	void push(Imgidx pidx);
	void pop();
	void find_minlev();

	inline Imgidx top() { return min_level; }
	inline _int64 get_minlev() { return min_level; }
};
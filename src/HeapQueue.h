#pragma once
#include <cfloat>
#include "defines.h"
#define PROFILE 0
#if PROFILE
#include <vector>
#endif

template<class Pixel>
class HQentry
{
public:
	HQentry() {}
	~HQentry() {}
	Imgidx pidx;
	Pixel alpha;
	inline void operator=(const HQentry& item)
	{
		this->pidx = item.pidx;
		this->alpha = item.alpha;
	}
};

// Heap-based priority queue
// Use push buffer to minimize push() and pop() call and time
template<class Pixel> //
class HeapQueue
{
	Imgidx maxsize;
	HQentry<Pixel> *arr;
	Pixel pop_level;
	Pixel max_level;
	HQentry<Pixel> pushlist[7]; // Connectivity - 1
	_uint8 flags[7];
	_int8 pushlistidx;

#if PROFILE
	std::vector<_uint64> num_memmove_push;
	std::vector<_uint64> num_memmove_pop;
	std::vector<_uint64> num_items_push;
	std::vector<_uint64> num_items_pop;

	_uint64 num_memmove_push_i;
	_uint64 num_memmove_pop_i;
#endif


public:
	Imgidx cursize;
	double qtime, timing;

	HeapQueue(Imgidx maxsize_in);
	~HeapQueue();
	int minidx_pushlist();
	void find_minlev();
	Imgidx pop();
	void push(Imgidx pidx, Pixel alpha);
	void push_run(Imgidx pidx, Pixel alpha);
	inline Pixel get_minlev() { return arr[1].alpha; }
	inline Imgidx top() { return arr[1].pidx; }
};

// Heap-based priority queue
template<class Pixel>
class HeapQueue_naive
{
	Imgidx cursize;
	Imgidx maxsize;
	HQentry<Pixel> *arr;
	Pixel pop_level;

#if PROFILE
	std::vector<_uint64> num_memmove_push;
	std::vector<_uint64> num_memmove_pop;
	std::vector<_uint64> num_items_push;
	std::vector<_uint64> num_items_pop;

	_uint64 num_memmove_push_i;
	_uint64 num_memmove_pop_i;
#endif

public:
	double qtime; //tmp
	HeapQueue_naive(Imgidx maxsize_in);
	~HeapQueue_naive();
	Imgidx pop();
	void push(Imgidx pidx, Pixel alpha);
	inline Pixel get_minlev() { return arr[1].alpha; }
	inline Imgidx top() { return arr[1].pidx; }
	inline Pixel top_alpha() { return arr[1].alpha; }
};

// Heap-based priority queue
template<class Pixel>
class HeapQueue_naive_quad
{
	Imgidx cursize;
	Imgidx maxsize;
	HQentry<Pixel> *arr;
	Pixel pop_level;

public:
	double qtime; //tmp
	inline Imgidx get_cursize(){return cursize;}
	HeapQueue_naive_quad(Imgidx maxsize_in);
	~HeapQueue_naive_quad();

#if PROFILE
	_uint64 pop();
	_uint64 push(Imgidx pidx, Pixel alpha);
#else
	Imgidx pop();
	void push(Imgidx pidx, Pixel alpha);
#endif

	inline Pixel get_minlev() { return arr[1].alpha; }
	inline Imgidx top() { return arr[1].pidx; }
	inline Pixel top_alpha() {return arr[1].alpha;}
};

//Heap-based priority queue
class HeapQueue_rank
{
	Imgidx *arr;
public:
	Imgidx cursize;
	Imgidx maxsize;
	Imgidx min_level;
	HeapQueue_rank(Imgidx maxsize_in);
	~HeapQueue_rank();
	Imgidx pop();
	void push(Imgidx pidx);

	inline Imgidx get_minlev() { return arr[1]; }
	inline void find_minlev() {};
	inline Imgidx top() { return arr[1]; }
};

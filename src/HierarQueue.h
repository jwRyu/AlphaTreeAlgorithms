#pragma once

#include <cstdio>
#include "defines.h"
#include "allocator.h"

#define HQUEUE_DEBUG 0
#if HQUEUE_DEBUG
#include <iostream>
using namespace std;
#endif

//Do not use beyond 20-bit image
template <class Imgidx>
class HierarQueue//: public PriorityQueue<Imgidx>
{
public://tmp
	Imgidx *queue;
	Imgidx *bottom, *cur;
	_int32 numlevel;
#if HQUEUE_DEBUG
	Imgidx *cnt;
	Imgidx *curmax;
#endif

	_int64 qsize;
	_int64 min_level, max_level;

	void print()
	{
		for(_int32 i = min_level;i < numlevel;i++)
		{
			if(cur[i] == bottom[i])
				continue;
			printf("Level %d(%d-%d): ", i, (int)bottom[i], (int)cur[i] - 1);
			for(Imgidx j = cur[i] - 1;j != bottom[i];j--)
				printf("%d ", (int)queue[j]);
			printf("%d ", (int)queue[bottom[i]]);
			printf("\n");
		}
	}

	HierarQueue(_uint64 qsize_in, _int32 numlevels)
	{
		//tmp
		queue = (Imgidx*)Malloc((size_t)qsize_in * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));

	#if HQUEUE_DEBUG
		cnt = (Imgidx*)Calloc((size_t)(numlevels) * sizeof(Imgidx));
		curmax = (Imgidx*)Calloc((size_t)(numlevels + 1) * sizeof(Imgidx));
	#endif
		this->numlevel = numlevels;
		max_level = 0;

		qsize = qsize_in;
		min_level = numlevels - 1;

		bottom[numlevels] = 0;
		cur[numlevels] = 1;
#if HQUEUE_DEBUG
		//cout << "  queue init - minlev = " << min_level << " / qsize = " << qsize_in << endl;
#endif
	}


	HierarQueue(_uint64 qsize_in)
	{
		_int32 numlevels = (_int32)1 << 20;
		//tmp
		queue = (Imgidx*)Malloc((size_t)qsize_in * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc(((size_t)numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc(((size_t)numlevels + 1) * sizeof(Imgidx));

	#if HQUEUE_DEBUG
		cnt = (Imgidx*)Calloc((size_t)(numlevels) * sizeof(Imgidx));
		curmax = (Imgidx*)Calloc((size_t)(numlevels + 1) * sizeof(Imgidx));
	#endif
		this->numlevel = numlevels;
		max_level = 0;

		qsize = qsize_in;
		min_level = numlevels - 1;

		bottom[numlevels] = 0;
		cur[numlevels] = 1;
#if HQUEUE_DEBUG
		//cout << "  queue init - minlev = " << min_level << " / qsize = " << qsize_in << endl;
#endif
	}

	void reset_queue()
	{
		for (_int32 i = 0; i < numlevel; i++)
			cur[i] = bottom[i];
		min_level = numlevel;
	}

	Imgidx set_queue(Imgidx *dhist)
	{
		Imgidx sum_hist = 0;
		for (_int32 i = 0; i < numlevel; i++)
		{
			bottom[i] = cur[i] = sum_hist;

			#if HQUEUE_DEBUG
				curmax[i] = cur[i];
			#endif

			if(dhist[i])
			{
				max_level = i;
				sum_hist += dhist[i];
			}
		}
		return sum_hist;
	}

	Imgidx set_queue(Imgidx *dhist, _int32 maxpix)
	{
		Imgidx sum_hist = 0;
		_int32 numlevels = (_int32)1 << 20;

		for(int i = 0;i < numlevels;i++)
			bottom[i] = cur[i] = 0;
		for(int i = 0;i < qsize;i++)
			queue[i] = 0;

		max_level = 0;
		for (_int32 i = 0; i <= maxpix; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			if(dhist[i])
			{
				max_level = i;
				sum_hist += dhist[i];
			}
		}

		min_level = maxpix + 1;
		bottom[maxpix + 1] = 0;
		cur[maxpix + 1] = 1;

		return sum_hist;
	}

	HierarQueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
	{
		//tmp
		queue = (Imgidx*)Malloc((size_t)qsize_in * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));

#if HQUEUE_DEBUG
		cnt = (Imgidx*)Calloc((size_t)(numlevels) * sizeof(Imgidx));
		curmax = (Imgidx*)Calloc((size_t)(numlevels + 1) * sizeof(Imgidx));
#endif

		this->numlevel = numlevels;
		max_level = 0;

		qsize = qsize_in;
		min_level = numlevels - 1;

		Imgidx sum_hist = 0;
		for (_int32 i = 0; i < numlevels; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			if(dhist[i])
			{
				max_level = i;
				sum_hist += dhist[i];
			}
		}

		bottom[numlevels] = 0;
		cur[numlevels] = 1;
	}


	HierarQueue(Imgidx *dhist, _int32 numlevels)
	{
		_uint64 dsum = 0;
		max_level = 0;
		for(int i = 0;i < numlevels;i++)
		{
			if(dhist[i])
			{
				dsum += dhist[i];
				max_level = i;
			}
		}

		//tmp

		queue = (Imgidx*)Malloc((size_t)(dsum + 1) * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));

#if HQUEUE_DEBUG
		cnt = (Imgidx*)Calloc((size_t)(numlevels) * sizeof(Imgidx));
#endif

		this->numlevel = numlevels;

		qsize = (dsum + 1);
		min_level = numlevels - 1;

		Imgidx sum_hist = 0;
		for (_int32 i = 0; i < numlevels; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i];
		}
		bottom[numlevels] = 0;
		cur[numlevels] = 1;
	}

	HierarQueue(_int32 numlevels, Imgidx binsize)
	{
		qsize = numlevels * binsize;
		//tmp
		queue = (Imgidx*)Malloc((size_t)qsize * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
#if HQUEUE_DEBUG
		cnt = (Imgidx*)Calloc((size_t)(numlevels) * sizeof(Imgidx));
#endif
		min_level = numlevels - 1;
		this->numlevel = numlevels;

		Imgidx sum_hist = 0;
		for (_int32 i = 0; i < numlevels; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += binsize;
		}
		bottom[numlevels] = 0;
		cur[numlevels] = 1;
	}

	~HierarQueue()
	{
//tmptmp

#if HQUEUE_DEBUG
		//for(int i = 0;i < numlevel;i++)
		//	cout << cnt[i] << " ";
		//cout << endl;
		Free(cnt);
		Free(curmax);
#endif

		Free(queue);
		Free(bottom);
		Free(cur);
	}

#if HQUEUE_DEBUG
	Imgidx numfreeslots()
	{
		Imgidx numfreeslots = 0;
		for(int i = 0;i < 1;i++)
		{
			if(i == numlevel - 1)
				numfreeslots += qsize - curmax[i];
			else
				numfreeslots += bottom[i + 1] - curmax[i];
		}
		return numfreeslots;
	}
#endif

	inline _int8 push(Imgidx pidx, _int64 level)
	{
#if HQUEUE_DEBUG
	cnt[level]++;
	cout << "  push " << pidx << " at " << level << " cur[" << level << "] = " << cur[level]+1 << "(minlev= " << min_level << ")" << endl;
#endif
		queue[cur[level]++] = pidx;

		//#if HQUEUE_DEBUG
		//	if(cur[level] > curmax[level])
		//		curmax[level] = cur[level];
		//#endif


		if (level < min_level)
		{
#if HQUEUE_DEBUG
			//cout << "  minlev updated from " << min_level << " to " << level << endl;
#endif
			min_level = level;
			return 1;
		}
		else
			return 0;
		//min_level = min(level, min_level);
	}

//pop operated as stack for the classic alpha tree implementation (i.e. non-hypergraph)... but in max tree it has to be changed.
	inline Imgidx pop()
	{
#if HQUEUE_DEBUG
		//cout << "  pop " << queue[bottom[min_level]] << " at " << min_level << " bottom[" << min_level << "] = " << bottom[min_level]+1 << "(minlev= " << min_level << ")" << endl;
#endif
		return queue[bottom[min_level]++];
		//return queue[bottom[min_level]++]; //alternative implementation
	}

//inline Imgidx pop()
//{
//#if HQUEUE_DEBUG
//cout << "pop " << queue[cur[min_level] - 1] << " at " << min_level << endl;
//#endif
//	return queue[--cur[min_level]];
//}
	inline Imgidx top()
	{
		return queue[bottom[min_level]];
	}
	inline _int64 get_minlev()
	{
		return min_level;
	}
	inline void find_minlev()
	{
		while (bottom[min_level] == cur[min_level])
			min_level++;
	}
};

template <class Imgidx>
class HQueue_l1idx
{
public://tmp
	Imgidx *queue;
	Imgidx *bottom, *cur;
	_uint64 *seeker;
	//move public: here

	_int64 qsize, seekersize;
	_int32 min_level;
	HQueue_l1idx(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
	{
		//tmp
		queue = (Imgidx*)Malloc((size_t)qsize_in * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		seekersize = (numlevels + 1 + 63) >> 6;
		seeker = (_uint64 *)Malloc((size_t)(seekersize) * sizeof(_uint64));

		qsize = qsize_in;
		min_level = numlevels - 1;

		Imgidx sum_hist = 0;
		for (_int32 i = 0; i < numlevels; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i];
		}
		for (_int64 i = 0; i < seekersize; i++)
			seeker[i] = 0;
		seeker[numlevels >> 6] |= (_uint64)1 << (numlevels & 63);
		bottom[numlevels] = 0;
		cur[numlevels] = 1;
	}
	~HQueue_l1idx()
	{
		Free(queue);
		Free(bottom);
		Free(cur);
		Free(seeker);
	}

	inline int push(Imgidx pidx, _int32 level)
	{
		_int64 qidx = cur[level]++;
		queue[qidx] = pidx;
		seeker[level >> 6] |= (_uint64)1 << (level & 63);
		if (level <= min_level)
		{
			min_level = level;
			return 1;
		}
		return 0;
		//		min_level = min(level, min_level);
	}

//stack -> queue. should change later
	inline Imgidx pop()
	{
		Imgidx popidx = --cur[min_level];

		if(bottom[min_level] == cur[min_level])
			seeker[min_level >> 6] &= ~((_uint64)1 << (min_level & 63));
		return queue[popidx];
	}

	inline Imgidx top()
	{
		return queue[cur[min_level] - 1];
	}

	inline _int32 get_minlev()
	{
		return min_level;
	}

	inline void find_minlev()
	{
		Imgidx qidx, widx;
		_uint64 w;

		for (qidx = min_level >> 6; !seeker[qidx]; qidx++)
			;

		w = seeker[qidx];

		if (w & 0xffffffff)
			widx = 0;
		else
		{
			widx = 32;
			w >>= 32;
		}

		while (!(w&(_uint64)1))
		{
			w >>= 1;
			widx++;
		}

		min_level = ((qidx << 6) + widx);
//
// 		while (bottom[min_level] == cur[min_level])
// 			min_level++;
//
//
// 		if ( != min_level)
// 			qidx = qidx;
	}
};

template <class Imgidx>
class HQueue_l2idx
{
	Imgidx *queue;
	Imgidx *bottom, *cur;
	_uint64 *seeker,*seeker2;
public:
	_int64 qsize;
	_int64 min_level;
	HQueue_l2idx(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
	{
		_int64 seekersize, seeker2size;
		//tmp
		queue = (Imgidx*)Malloc((size_t)qsize_in * sizeof(Imgidx));
		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));
		cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx));

		seekersize = (numlevels + 1 + 63) >> 6;
		seeker2size = (seekersize + 63) >> 6;

		seeker = (_uint64 *)Malloc((size_t)(seekersize) * sizeof(_uint64));
		seeker2 = (_uint64 *)Malloc((size_t)(seeker2size) * sizeof(_uint64));

		qsize = qsize_in;
		min_level = numlevels - 1;

		Imgidx sum_hist = 0;
		for (_int32 i = 0; i < numlevels; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i];
		}
		for (_int64 i = 0; i < seekersize; i++)
			seeker[i] = 0;
		seeker[numlevels >> 6] |= (_uint64)1 << (numlevels & 63);
		for (_int64 i = 0; i < seeker2size; i++)
			seeker2[i] = 0;
		seeker2[numlevels >> 12] |= (_uint64)1 << ((numlevels >> 6) & 63);
		bottom[numlevels] = 0;
		cur[numlevels] = 1;
	}
	~HQueue_l2idx()
	{
		Free(queue);
		Free(bottom);
		Free(cur);
		Free(seeker);
		Free(seeker2);
	}

	inline void push(Imgidx pidx, _int64 level)
	{
		_int64 qidx = cur[level]++;
		queue[qidx] = pidx;
		seeker[level >> 6] |= (_uint64)1 << (level & 63);
		seeker2[level >> 12] |= (_uint64)1 << ((level >> 6) & 63);
		if (level < min_level)
		{
			min_level = level;
			//return 1;
		}
		//		else
					//return 0;
		//		min_level = min(level, min_level);
	}

//stack -> queue. should change later
	inline Imgidx pop()
	{
		Imgidx popidx = --cur[min_level];

		if (bottom[min_level] == cur[min_level])
		{
			seeker[min_level >> 6] &= ~((_uint64)1 << (min_level & 63));
			if (!seeker[min_level >> 6])
				seeker2[min_level >> 12] &= ~((_uint64)1 << ((min_level >> 6) & 63));
		}
		return queue[popidx];
	}

	inline Imgidx top()
	{
		return queue[cur[min_level] - 1];
	}

	inline _int64 get_minlev()
	{
		return min_level;
	}

	inline void find_minlev()
	{
		Imgidx qidx, widx;
		_uint64 w;

		for (qidx = min_level >> 12; !seeker2[qidx]; qidx++)
			;

		w = seeker2[qidx];
		if (w & 0xffffffff)
			widx = 0;
		else
		{
			widx = 32;
			w >>= 32;
		}

		while (!(w&(_uint64)1))
		{
			w >>= 1;
			widx++;
		}

		qidx = ((qidx << 6) + widx);

		w = seeker[qidx];
		if (w & 0xffffffff)
			widx = 0;
		else
		{
			widx = 32;
			w >>= 32;
		}

		while (!(w&(_uint64)1))
		{
			w >>= 1;
			widx++;
		}

		min_level = ((qidx << 6) + widx);
		//
		// 		while (bottom[min_level] == cur[min_level])
		// 			min_level++;
		//
		//
		// 		if ( != min_level)
		// 			qidx = qidx;
	}
};

template <class Imgidx>
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
//	double jumpdist, jumpnum;
	HQueue_l1idx_rank(_int64 qsize_in)
	{
// 		jumpnum = 0;
// 		jumpdist = 0;
		//tmp
		qsize = (qsize_in + (1<<12)) >> 12;
		queue = (hqueue_word *)Calloc((size_t)qsize * sizeof(hqueue_word));

		min_level = qsize_in;
	}
	~HQueue_l1idx_rank()
	{
		Free(queue);
	}

	inline void push(Imgidx pidx)
	{
		_int64 qidx = pidx >> 12;
		_int64 widx = (pidx >> 6) & 63;
		_int64 bitpos = pidx & 63;

		queue[qidx].qword[widx] |= (_int64)1 << (bitpos);
		queue[qidx].seeker |= (_int64)1 << (widx);
		min_level = min_level < pidx ? min_level : pidx;
	}

//stack -> queue. should change later
	inline void pop()
	{
		_int64 qidx = min_level >> 12;
		_int64 widx = (min_level >> 6) & 63;
		_int64 bitpos = min_level & 63;
		_int64 w, skr;

		queue[qidx].qword[widx] &= ~((_int64)1 << (bitpos));
		if(!queue[qidx].qword[widx])
			queue[qidx].seeker &= ~((_int64)1 << (widx));

		//find_min_level
		for (qidx = min_level >> 12; !queue[qidx].seeker; qidx++)
			;

		skr = queue[qidx].seeker;
		widx = (skr & 0xffffffff) ? 0 : 32;
		widx += ((skr >> widx) & 0xffff) ? 0 : 16;
		widx += ((skr >> widx) & 0xff) ? 0 : 8;
		widx += ((skr >> widx) & 0xf) ? 0 : 4;
		widx += ((skr >> widx) & 0x3) ? 0 : 2;
		widx += ((skr >> widx) & 0x1) ? 0 : 1;

		w = queue[qidx].qword[widx];
		bitpos = (w & 0xffffffff) ? 0 : 32;
		bitpos += ((w >> bitpos) & 0xffff) ? 0 : 16;
		bitpos += ((w >> bitpos) & 0xff) ? 0 : 8;
		bitpos += ((w >> bitpos) & 0xf) ? 0 : 4;
		bitpos += ((w >> bitpos) & 0x3) ? 0 : 2;
		bitpos += ((w >> bitpos) & 0x1) ? 0 : 1;

// 		jumpdist += (double)((qidx << 12) + (widx << 6) + bitpos) - (double)min_level;
// 		jumpnum++;

		min_level = (qidx << 12) + (widx << 6) + bitpos;
	}

	inline Imgidx top()
	{
		return min_level;
	}

	inline _int64 get_minlev()
	{
		return min_level;
	}

	inline void find_minlev()
	{
		_int64 qidx, widx, bitpos, w, skr;

		for (qidx = min_level >> 12; !queue[qidx].seeker; qidx++)
			;

		skr = queue[qidx].seeker;
		widx = (skr & 0xffffffff) ? 0 : 32;
		widx += ((skr >> widx) & 0xffff) ? 0 : 16;
		widx += ((skr >> widx) & 0xff) ? 0 : 8;
		widx += ((skr >> widx) & 0xf) ? 0 : 4;
		widx += ((skr >> widx) & 0x3) ? 0 : 2;
		widx += ((skr >> widx) & 0x1) ? 0 : 1;

		w = queue[qidx].qword[widx];
		bitpos = (w & 0xffffffff) ? 0 : 32;
		bitpos += ((w >> bitpos) & 0xffff) ? 0 : 16;
		bitpos += ((w >> bitpos) & 0xff) ? 0 : 8;
		bitpos += ((w >> bitpos) & 0xf) ? 0 : 4;
		bitpos += ((w >> bitpos) & 0x3) ? 0 : 2;
		bitpos += ((w >> bitpos) & 0x1) ? 0 : 1;

		min_level = (qidx << 12) + (widx << 6) + bitpos;
 	}
};




/*
template <class Imgidx>
class HQueue_ubr
{
	_int8 *queue;
//	Imgidx *bottom, *cur;
public:
	_uint64 qsize;
	Imgidx min_level;
	HQueue_ubr(_uint64 qsize)
	{
		//tmp
		queue = (_int8*)Malloc((size_t)(qsize + 1) * sizeof(_int8));
		min_level = (Imgidx)(qsize);

		for (_uint64 i = 0; i < qsize; i++)
			queue[i] = 0;
		//queue[qsize] = 1;
//		bottom = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx), 1);
	//	cur = (Imgidx*)Malloc((size_t)(numlevels + 1) * sizeof(Imgidx), 1);


		// 		qsize = qsize_in;
		// 		min_level = numlevels - 1;
		//
		// 		Imgidx sum_hist = 0;
		// 		for (_int32 i = 0; i < numlevels; i++)
		// 		{
		// 			bottom[i] = cur[i] = sum_hist;
		// 			sum_hist += dhist[i];
		// 		}
		// 		bottom[numlevels] = 0;
		// 		cur[numlevels] = 1;
	}
	~HQueue_ubr()
	{
		Free(queue);
		//Free(bottom, 1);
		//Free(cur, 1);
	}

	inline void push(Imgidx pidx)
	{
		min_level = min(pidx, min_level);
		queue[pidx] = 1;
	}

	inline Imgidx top()
	{
		return min_level;
	}

	inline void pop()
	{
		queue[min_level] = 0;
		//std::cout << "popping from" << min_level << std::endl;
		while (!queue[min_level])
			min_level++;
	}

	inline void find_minlev()
	{

	}
};
*/

/*

struct neighbouridx {
	//  -  3  -
	//  2  p  1
	//  -  0  -
	_uint8 neighbour;
	_uint32 pidx;
};
template <>
class HQueue <neighbouridx>
{
	neighbouridx *queue;
	_uint64 *bottom, *cur;
	_uint64 qsize;
	_uint64 min_level, max_level;
public:
	HQueue(_uint64 qsize, _uint64 *dhist, _uint32 dhistsize, _uint8 neighbours)
	{
		_uint64 nn = neighbours >> 1;
		int shamt;
		queue = (neighbouridx*)Malloc((size_t)qsize * nn * sizeof(neighbouridx));
		bottom = (_uint64*)Malloc((size_t)(dhistsize + 1) * sizeof(_uint64));
		cur = (_uint64*)Malloc((size_t)(dhistsize + 1) * sizeof(_uint64));

		this->qsize = qsize;
		min_level = max_level = dhistsize;

		for (shamt = -1; nn; nn >>= 1)
			shamt++;

		_uint64 sum_hist = 0;
		for (_uint64 i = 0; i < dhistsize; i++)
		{
			bottom[i] = cur[i] = sum_hist;
			sum_hist += dhist[i] << shamt;
		}
		bottom[dhistsize] = 0;
		cur[dhistsize] = 1;
}
	~HQueue()
	{
		Free(queue);
		Free(bottom);
		Free(cur);
	}

	inline void hqueue_push(_uint64 pidx, _uint64 level)
	{
		min_level = min(level, min_level);
#if HQUEUE_DEBUG
		assert(level < max_level);
		assert(cur[level] < qsize);
#endif
		queue[cur[level]++] = pidx;
	}

	inline T hqueue_pop()
	{
		return queue[--cur[min_level]];
	}

	inline void hqueue_find_min_level()
	{
		while (bottom[min_level] == cur[min_level])
			min_level++;
	}
};
inline void hqueue_push(HQueue<neighidx>* hqueue, _uint32 idx, _uint8 neighbor, _uint32 level)
{
	hqueue->min_level = min(level, hqueue->min_level);
#if HQUEUE_DEBUG
	assert(level < hqueue->max_level);
	assert(hqueue->cur[level] < hqueue->qsize);
#endif
	hqueue->queue[hqueue->cur[level]].pidx = idx;
	hqueue->queue[hqueue->cur[level]++].neighbour = neighbor;
}
*/

/*
void hqueue_new(HQueue<_uint32>** hqueue, _uint64 qsize, _uint32 *dhist, _uint32 dhistsize);
void hqueue_new(HQueue<neighidx>** hqueue, _uint64 qsize, _uint32 *dhist, _uint32 dhistsize, _uint8 neighbours);
*/

#pragma once

#define HEAPQ_STAT 0
#define HEAPQ_DEBUG 0
#if (HEAPQ_STAT || HEAPQ_DEBUG)
#include <iostream>
#include <fstream>
using namespace std;
#endif

#include <cfloat>
#include "walltime.h"//tmp

template<class Imgidx, class Pixel>
class HQentry
{
public:
	HQentry() {}
	~HQentry() {}
	Imgidx pidx;
	Pixel alpha;
	uint32 moves;
	uint32 cache_moves;
	uint8 pure_cache;

	inline void operator=(const HQentry& item)
	{
		this->pidx = item.pidx;
		this->alpha = item.alpha;
		this->moves = item.moves + 1;
		this->cache_moves = item.cache_moves;
		this->pure_cache = item.pure_cache;
	}
};

//Heap-based priority queue
//uses push buffer to minimize push() and pop() call and time
template<class Imgidx, class Pixel>
class HeapQueue
{
public:
	Imgidx maxsize;
	HQentry<Imgidx, Pixel> *arr;
	Pixel pop_level;
	Pixel max_level;
	HQentry<Imgidx, Pixel> pushlist[7];//Connectivity - 1
	uint8 flags[7];
	int8 pushlistidx;

	uint64 numcmp;

#if HEAPQ_STAT
	ofstream f;
	Imgidx popcnt;
#endif

		Imgidx cursize;
		//tmp
		double qtime, timing;

//	Pixel min_level;
	HeapQueue(Imgidx maxsize_in) : maxsize(maxsize_in), numcmp(0), cursize(0), qtime(0), timing(0)
	{
		arr = new HQentry<Imgidx, Pixel>[maxsize + 1];
		pushlistidx = 0;
#if HEAPQ_STAT
		this->f.open("D:/RUG/2019/TTMA_ISMM/qstat.dat", ofstream::out);
		popcnt = 0;
#endif
		//max_level = -1;
		flags[0] = flags[1] = flags[2] = flags[3] = flags[4] = flags[5] = flags[6] = 1;
	}

	~HeapQueue()
	{
		delete[] arr;

#if HEAPQ_STAT
		f.close();
#endif
	}

	void printQ()
	{
		printf("Q content: ");
		for(int i = 1;i <= cursize;i++)
		{
			printf("%d(%f) ", (int)arr[i].pidx, (double)log2(1 + arr[i].alpha));
		}
		printf("\n");
	}

	inline uint8 is_empty(){return cursize < 1;}

	inline int minidx_pushlist()
	{
		int minidx = -1;
		for (int i = 0; i < pushlistidx; i++)
		{
			numcmp++;
			if (flags[i] && (minidx < 0 || pushlist[i].alpha < pushlist[minidx].alpha))
				minidx = i;
		}
		flags[minidx] = 0;
		return minidx;
	}

	inline void find_minlev()
	{
		//timing = 1;
		//double tt = get_wall_time();//tmp

		if(!pushlistidx)
		{
			pop();
#if HEAPQ_DEBUG
		cout << "find_minlev: no item in the list - new minlev is now " << arr[1].pidx << " at " << log2(1 + arr[1].alpha) << " cursize: " << cursize << endl;
#endif
			return;
		}
		int minidx = minidx_pushlist();
		numcmp++;
		if (pushlist[minidx].alpha <= arr[1].alpha)
		{
#if HEAPQ_DEBUG
			cout << "pop: " << arr[1].pidx << " at " << log2(1 + arr[1].alpha) << " cursize: " << cursize << endl;
#endif
			arr[1].moves = pushlist[minidx].moves + 1;
			arr[1].alpha = pushlist[minidx].alpha;
			arr[1].pidx = pushlist[minidx].pidx;
		}
		else
		{
			pop();
			push_run(pushlist[minidx].pidx, pushlist[minidx].alpha, pushlist[minidx].moves);
		}

		for (int i = 1; i < pushlistidx; i++)
		{
			minidx = minidx_pushlist();
			pushlist[minidx].moves++;
			push_run(pushlist[minidx].pidx, pushlist[minidx].alpha, pushlist[minidx].moves);
		}
#if HEAPQ_DEBUG
		//cout << "find_minlev: new minlev is now " << arr[1].pidx << " at " << log2(1 + arr[1].alpha) << endl;
#endif
//		min_level = arr[1].alpha;
		pushlistidx = 0;
		flags[0] = flags[1] = flags[2] = flags[3] = flags[4] = flags[5] = flags[6] = 1;

		//qtime += get_wall_time() - tt;//tmp
		//timing = 0;
	}

	inline Pixel get_minlev() { return arr[1].alpha; }
	inline Pixel top_alpha() { return arr[1].alpha; }

	inline Imgidx top() { return arr[1].pidx; }
	inline Imgidx top_moves() { return arr[1].moves; }
	inline uint32 top_cache_moves() { return arr[1].cache_moves; }

	//inline Imgidx pop()
	//{
		//pop_level = arr[1].alpha;
//		pushlistidx = 0;
//#if HEAPQ_DEBUG
//		cout << "pop: " << arr[1].pidx << " at " << arr[1].alpha << endl;
//#endif
//		return arr[1].pidx;
//	}

	inline Imgidx pop()
	{
		//double tt = 0;
		//if(!timing)			tt = get_wall_time(); //tmp

		Imgidx outval = arr[1].pidx;
		Imgidx current = 1, next, next0, next1, curidx;
		Pixel curalpha;

		uint32 curmove;

#if HEAPQ_DEBUG
		cout << "pop: " << arr[1].pidx << " at " << log2(1 + arr[1].alpha) << " cursize: " << cursize << endl;
#endif
		// 		ulong val_out = queue->array[1];
		// 		ulong current = 1, moved;
		// 		value val_cur;


		curidx = arr[cursize].pidx;
		curalpha = arr[cursize].alpha;
		curmove = arr[cursize].moves;
//		if (curidx < 0)
	//		curidx = curidx;
		cursize--;

#if HEAPQ_STAT
		f << "1" << std::endl << (int)arr[1].alpha << std::endl;
#endif
		while (1)
		{
			next0 = current << 1;
			next1 = next0 + 1;
			if (next0 > cursize)
				break;
				
			numcmp++;
			if (next1 <= cursize && arr[next1].alpha < arr[next0].alpha)
				next = next1;
			else
				next = next0;

			numcmp++;
			if (curalpha < arr[next].alpha)
				break;

			arr[current] = arr[next];
			current = next;
		}
		arr[current].alpha = curalpha;
		arr[current].pidx = curidx;
		arr[current].moves = curmove + 1;
//
// 		if (cursize)
// 			min_level = arr[1].alpha;
// 		else
// 			min_level = (Pixel)-1;

		//validate();

		//if(!timing)		qtime += get_wall_time() - tt; //tmp

		return outval;
	}

	inline void push(Imgidx pidx, Pixel alpha, uint32 moves = 0)
	{
		//double tt = get_wall_time(); //tmp
		pushlist[pushlistidx].moves = moves;
		pushlist[pushlistidx].pidx = pidx;
		pushlist[pushlistidx++].alpha = alpha;
		//qtime += get_wall_time() - tt; //tmp
#if HEAPQ_DEBUG
		//cout << "push: " << pidx << " at " << log2(1 + alpha) << endl;
#endif
	}

// 	void push(Imgidx pidx)
// 	{
// 		push_run(pidx, (Pixel)pidx);
// 	}

	inline void push_run(Imgidx pidx, Pixel alpha, uint32 moves = 0)
	{
		Imgidx current, next;

#if HEAPQ_STAT
		f << "0" << std::endl << (int)alpha << std::endl;
#endif
		//		value val_cur = tree[pixpos].gval;
		cursize++;
		current = cursize;

	//	if (pidx < 0)
	//		pidx = pidx;

		next = current >> 1;
		numcmp++;
		while (next && (arr[next].alpha > alpha))
		{
			numcmp++;
			arr[current] = arr[next];
			current = next;
			next = next >> 1;
		}

		arr[current].pidx = pidx;
		arr[current].alpha = alpha;
		arr[current].moves = moves + 1;
		arr[current].cache_moves = 0;
#if HEAPQ_DEBUG
		cout << "pushrun: " << pidx << " at " << log2(1 + arr[1].alpha) << " cursize: " << cursize << endl;
#endif

// 		if (current == 1)
// 			min_level = alpha;
		//validate();
	}
};

//Heap-based priority queue
template<class Imgidx, class Pixel>
class HeapQueue_naive
{
	Imgidx cursize;
	Imgidx maxsize;
	HQentry<Imgidx, Pixel> *arr;
	Pixel pop_level;
	//HQentry<Imgidx, Pixel> pushlist[7];//Connectivity - 1
	//int8 pushlistidx;

#if HEAPQ_STAT
	ofstream f;
	Imgidx popcnt;
#endif

public:
	double qtime; //tmp
//	Pixel min_level;
	HeapQueue_naive(Imgidx maxsize_in) : cursize(0), maxsize(maxsize_in), qtime(0)
	{
		arr = new HQentry<Imgidx, Pixel>[maxsize + 1];
#if HEAPQ_STAT
		this->f.open("D:/RUG/2019/TTMA_ISMM/qstat.dat", ofstream::out);
		popcnt = 0;
#endif

		if((Pixel)-1 > (Pixel)1)
			arr[1].alpha = (Pixel)-1;
		else if(sizeof(Pixel) == 8)
			arr[1].alpha = DBL_MAX;
		else
			arr[1].alpha = FLT_MAX;
	}

	~HeapQueue_naive()
	{
		delete[] arr;

#if HEAPQ_STAT
		f.close();
#endif
	}

	inline Pixel get_minlev() { return arr[1].alpha; }

/*
	inline void validate()
	{
		Imgidx i, j;

		for (i = cursize; i > 1; i--)
		{
			j = i >> 1;
			if (arr[i].alpha < arr[j].alpha)
			{
				std::cout << "invalidate queue" << std::endl;
				std::cin >> i;
			}
		}
	}*/

	inline Imgidx top() { return arr[1].pidx; }
	inline Pixel top_alpha() { return arr[1].alpha; }

	inline Imgidx pop()
	{
		//double tt = get_wall_time(); //tmp

		Imgidx outval = arr[1].pidx;
		Imgidx current = 1, next, next0, next1, curidx;
		Pixel curalpha;

#if HEAPQ_DEBUG
		cout << "pop: " << arr[1].pidx << " at " << arr[1].alpha << " cursize: " << cursize << endl;
#endif
		// 		ulong val_out = queue->array[1];
		// 		ulong current = 1, moved;
		// 		value val_cur;


		curidx = arr[cursize].pidx;
		curalpha = arr[cursize].alpha;
//		if (curidx < 0)
	//		curidx = curidx;
		cursize--;

#if HEAPQ_STAT
		f << "1" << std::endl << (int)arr[1].alpha << std::endl;
#endif
		while (1)
		{
			next0 = current << 1;
			next1 = next0 + 1;
			if (next0 > cursize)
				break;
			if (next1 <= cursize && arr[next1].alpha < arr[next0].alpha)
				next = next1;
			else
				next = next0;

			if (curalpha < arr[next].alpha)
				break;

			arr[current] = arr[next];
			current = next;
		}
		arr[current].alpha = curalpha;
		arr[current].pidx = curidx;
//
// 		if (cursize)
// 			min_level = arr[1].alpha;
// 		else
// 			min_level = (Pixel)-1;

		//validate();

		//qtime += get_wall_time() - tt; //tmp
		return outval;
	}

// 	void push(Imgidx pidx)
// 	{
// 		push_run(pidx, (Pixel)pidx);
// 	}

	inline void push(Imgidx pidx, Pixel alpha)
	{
		//double tt = get_wall_time(); //tmp

		Imgidx current, next;

#if HEAPQ_STAT
		f << "0" << std::endl << (int)alpha << std::endl;
#endif
		//		value val_cur = tree[pixpos].gval;
		cursize++;
		current = cursize;

	//	if (pidx < 0)
	//		pidx = pidx;

		next = current >> 1;
		while (next && (arr[next].alpha > alpha))
		{
			arr[current] = arr[next];
			current = next;
			next = next >> 1;
		}

		arr[current].pidx = pidx;
		arr[current].alpha = alpha;
#if HEAPQ_DEBUG
		cout << "pushrun: " << pidx << " at " << alpha << " cursize: " << cursize << endl;
#endif

// 		if (current == 1)
// 			min_level = alpha;
		//validate();
		//qtime += get_wall_time() - tt; //tmp
	}
};

//Heap-based priority queue
template<class Imgidx, class Pixel>
class HeapQueue_naive_quad
{
public:
	uint64 numcmp;
	Imgidx cursize;
	Imgidx cursize_max;
	Imgidx maxsize;
	HQentry<Imgidx, Pixel> *arr;
	Pixel pop_level;
	//HQentry<Imgidx, Pixel> pushlist[7];//Connectivity - 1
	//int8 pushlistidx;

#if HEAPQ_STAT
	ofstream f;
	Imgidx popcnt;
#endif

	double qtime; //tmp
//	Pixel min_level;
	inline Imgidx get_cursize(){return cursize;}
	HeapQueue_naive_quad(Imgidx maxsize_in) : numcmp(0), cursize(0), cursize_max(0), maxsize(maxsize_in), qtime(0)
	{
		//arr = new HQentry<Imgidx, Pixel>[maxsize + 2];
		arr = (HQentry<Imgidx, Pixel>*)Calloc(sizeof(HQentry<Imgidx, Pixel>) * (maxsize + 2));
#if HEAPQ_STAT
		this->f.open("D:/RUG/2019/TTMA_ISMM/qstat.dat", ofstream::out);
		popcnt = 0;
#endif
		if((Pixel)-1 > (Pixel)1)
		 	arr[1].alpha = (Pixel)-1;
		else if(sizeof(Pixel) == 8)
			arr[1].alpha = DBL_MAX;
		else
			arr[1].alpha = FLT_MAX;
	}

	~HeapQueue_naive_quad()
	{
		//delete[] arr;
		Free(arr);

#if HEAPQ_STAT
		f.close();
#endif
	}

	inline Pixel get_minlev() { return arr[1].alpha; }
	inline Imgidx top() { return arr[1].pidx; }
	inline Pixel top_alpha() {return arr[1].alpha;}
	inline uint32 top_moves() { return arr[1].moves; }
	inline uint32 top_cache_moves() { return arr[1].cache_moves; }
//	{
//		 Pixel ret = arr[1].alpha;
//		 return ret;
//	}

	inline Imgidx pop()
	{
		//double tt = get_wall_time(); //tmp

		Imgidx outval = arr[1].pidx;
		Imgidx current = 1, next, next0, curidx;
		Pixel curalpha;
		uint32 curmove, curcachemove;

#if HEAPQ_DEBUG
		cout << "pop: " << arr[1].pidx << " at " << arr[1].alpha << endl;
#endif
		// 		ulong val_out = queue->array[1];
		// 		ulong current = 1, moved;
		// 		value val_cur;


		curidx = arr[cursize].pidx;
		curalpha = arr[cursize].alpha;
		curmove = arr[cursize].moves;
		curcachemove = arr[cursize].cache_moves;
//		if (curidx < 0)
	//		curidx = curidx;
		cursize--;
		if(cursize == 0)
			return 0;

#if HEAPQ_STAT
		f << "1" << std::endl << (int)arr[1].alpha << std::endl;
#endif
		while (1)
		{
			next0 = (current << 2) - 2;
			next = next0;
			//next1 = next0 + 1;
			if(next0 + 3 <= cursize)
			{
				if (arr[next0 + 1].alpha < arr[next].alpha)	{next = next0 + 1; numcmp++;}
				if (arr[next0 + 2].alpha < arr[next].alpha)	{next = next0 + 2; numcmp++;}
				if (arr[next0 + 3].alpha < arr[next].alpha)	{next = next0 + 3; numcmp++;}
			}
			else
			{
				if (next0 > cursize)	break;
				if (next0 == cursize)			    	goto MIN_NEXT_FOUND;
				if (arr[next0 + 1].alpha < arr[next].alpha)	{next = next0 + 1; numcmp++;}
				if (next0 + 1 == cursize)				goto MIN_NEXT_FOUND;
				if (arr[next0 + 2].alpha < arr[next].alpha)	{next = next0 + 2; numcmp++;}
				if (next0 + 2 == cursize)				goto MIN_NEXT_FOUND;
				if (arr[next0 + 3].alpha < arr[next].alpha)	{next = next0 + 3; numcmp++;}
			}

MIN_NEXT_FOUND:
			numcmp++;
			if (curalpha < arr[next].alpha)
				break;

			arr[current] = arr[next];
			current = next;
		}
		arr[current].alpha = curalpha;
		arr[current].pidx = curidx;
		arr[current].moves = curmove + 1;
		arr[current].cache_moves = curcachemove;
//
// 		if (cursize)
// 			min_level = arr[1].alpha;
// 		else
// 			min_level = (Pixel)-1;

		//validate();

		//qtime += get_wall_time() - tt; //tmp
		return outval;
	}

// 	void push(Imgidx pidx)
// 	{
// 		push_run(pidx, (Pixel)pidx);
// 	}

	inline void push(Imgidx pidx, Pixel alpha, uint32 moves = 0, uint32 cache_moves = 0)
	{
		//double tt = get_wall_time(); //tmp

		Imgidx current, next;

#if HEAPQ_STAT
		f << "0" << std::endl << (int)alpha << std::endl;
#endif
		//		value val_cur = tree[pixpos].gval;
		cursize++;
		if(cursize > cursize_max) 
			cursize_max = cursize;

		current = cursize;

	//	if (pidx < 0)
	//		pidx = pidx;

		next = (current + 2) >> 2;
		numcmp++;
		while (next && (arr[next].alpha > alpha))
		{
			numcmp++;
			arr[current] = arr[next];
			current = next;
			next = (next + 2) >> 2;
		}

		arr[current].pidx = pidx;
		arr[current].alpha = alpha;
		arr[current].moves = moves;
		arr[current].cache_moves = cache_moves;
#if HEAPQ_DEBUG
		cout << "push: " << pidx << " at " << alpha << " at slot[" << current << "]" << endl;
#endif

// 		if (current == 1)
// 			min_level = alpha;
		//validate();
		//qtime += get_wall_time() - tt; //tmp
	}
};

//Heap-based priority queue
template<class Imgidx>
class HeapQueue_rank
{

	Imgidx *arr;
#if HEAPQ_STAT
	ofstream f;
	Imgidx *depth;
	Imgidx cur_depth;
	Imgidx popcnt;
#endif
public:
	Imgidx cursize;
	Imgidx maxsize;
	Imgidx min_level;
	HeapQueue_rank(Imgidx maxsize_in) : cursize(0), maxsize(maxsize_in)
	{
		arr = new Imgidx[maxsize];
		arr--;
#if HEAPQ_STAT
		popcnt = 0;
		f.open("D:/RUG/2019/TTMA_ISMM/qstat.dat", std::ofstream::app);
		f << -1 << '\n' << maxsize << endl;
		depth = (Imgidx*)Calloc(maxsize * sizeof(Imgidx));
#endif

	}
	~HeapQueue_rank()
	{
		delete[] (arr + 1);

#if HEAPQ_STAT
		f.close();
#endif
	}

	inline Imgidx get_minlev() { return arr[1]; }
	inline void find_minlev() {};
	/*
	inline void validate()
	{
		Imgidx i, j;

		for (i = cursize; i > 1; i--)
		{
			j = i >> 1;
			if (arr[i].alpha < arr[j].alpha)
			{
				std::cout << "invalidate queue" << std::endl;
				std::cin >> i;
			}
		}
	}
	*/
	inline Imgidx top() { return arr[1]; }
	inline Imgidx pop()
	{
		Imgidx outval;
		Imgidx current = 1, next, next0, next1, curidx;
		//Imgidx curalpha;

		// 		ulong val_out = queue->array[1];
		// 		ulong current = 1, moved;
		// 		value val_cur;


		outval = arr[current];
		curidx = arr[cursize--];
#if HEAPQ_STAT
		Imgidx curdepth = depth[cursize + 1];
		cur_depth = cur_depth > depth[1] ? cur_depth : depth[1];
		f << '1' << std::endl << arr[1] << std::endl << cur_depth++ << std::endl;
		if (popcnt++ % 3000 == 0)
			cout << "pop " << popcnt << '/' << maxsize << endl;
#endif
		while (1)
		{
			next0 = current << 1;
			next1 = next0 + 1;
			if (next0 > cursize)
				break;
			if (next1 <= cursize && arr[next1] < arr[next0])
				next = next1;
			else
				next = next0;

			if (curidx < arr[next])
				break;

			arr[current] = arr[next];
#if HEAPQ_STAT
			depth[current] = depth[next];
#endif
			current = next;
		}
		arr[current] = curidx;
#if HEAPQ_STAT
		depth[current] = curdepth;
		f << arr[1] - min_level << endl;
#endif

		if (cursize)
			min_level = arr[1];
		else
			min_level = (Imgidx)-1;

		//validate();

		return outval;
	}
	inline void push(Imgidx pidx)
	{
		Imgidx current, next;

#if HEAPQ_STAT
		f << '0' << '\n' << pidx << endl;
#endif
		//		value val_cur = tree[pixpos].gval;
		cursize++;
		current = cursize;

		//	if (pidx < 0)
		//		pidx = pidx;

		next = current >> 1;
		while (next && (arr[next] > pidx))
		{
			arr[current] = arr[next];
#if HEAPQ_STAT
			depth[current] = depth[next];
#endif
			current = next;
			next = next >> 1;
		}

		arr[current] = pidx;
#if HEAPQ_STAT
		depth[current] = 0;
#endif

		if (current == 1)
		{
			min_level = pidx;
#if HEAPQ_STAT
			if (cursize > 2)
			{
				if (arr[2] < arr[3])
					depth[2] = depth[2] > cur_depth? depth[2] : cur_depth;
				else
					depth[3] = depth[3] > cur_depth ? depth[3] : cur_depth;
			}
			else if (cursize > 1)
				depth[2] = depth[2] > cur_depth ? depth[2] : cur_depth;
			cur_depth = 0;
#endif
		}
		//validate();
	}

	//for hybrid queueing
// 	inline Imgidx get_arr(Imgidx arridx) { return arr[arridx]; }
// 	inline void set_arr(Imgidx arridx, Imgidx val) { arr[arridx] = val; }
};

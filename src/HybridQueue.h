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
using namespace std;

#define LISTSIZE_DEFAULT 12
#define HEAPSIZE_DEFAULT 128

#define TRACK_QUEUEING	0
#define QUEUE_DEBUG		0

#define PROFILE			0

//tmp
#if (TRACK_QUEUEING || QUEUE_DEBUG)
#include <iostream>
#include <fstream>
#include <cfloat>
using namespace std;
#endif

template<class Imgidx, class Pixel>
class HierarHeapQueue_HEQ
{
	HQentry<Imgidx, Pixel> *list;
	HeapQueue_naive_quad<Imgidx, Pixel> **hqueue;
	HQentry<Imgidx, Pixel> **storage;
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
	void initHQ(Imgidx *dhist, _uint32* histeqmap_in, Imgidx numlevels_in, Imgidx size, double a_in, int listsize)
	{
		/*		cnt = 0;//tmp*/
		//Imgidx i;

		list = (HQentry<Imgidx, Pixel>*)Malloc(listsize * sizeof(HQentry<Imgidx, Pixel>));
		maxSize_list = listsize - 1;
		curSize_list = -1;

		this->numlevels = numlevels_in;
		this->a = a_in;
		this->queue_minlev = numlevels;
		this->histeqmap = histeqmap_in; //do not free dhist outside
		this->qsizes = dhist; //do not free dhist outside

		storage_cursize = (Imgidx*)Calloc(numlevels * sizeof(Imgidx));

		Imgidx cumsum = 0;
		Imgidx thr_nonredundantnodes = (Imgidx)(size * 0.6);
		for(int level = 0;level < numlevels;level++)
		{
			cumsum += qsizes[level];
			if(cumsum > thr_nonredundantnodes)
			{
				thr_hqueue = curthr = level;
				break;
			}
		}

		hqueue = (HeapQueue_naive_quad<Imgidx, Pixel>**)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Imgidx, Pixel>*));
		for(int level = 0;level < thr_hqueue;level++)
			hqueue[level] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[level]);
		//hqueue[numlevels] = new HeapQueue_naive_quad<Imgidx, Pixel>(1);
		//hqueue[numlevels]->push(0,(Pixel)-1);

		storage = (HQentry<Imgidx, Pixel>**)Calloc((numlevels - thr_hqueue) * sizeof(HQentry<Imgidx, Pixel>*));
		storage -= thr_hqueue;
		for(int level = thr_hqueue;level < numlevels;level++)
			storage[level] = (HQentry<Imgidx, Pixel>*)Malloc(qsizes[level] * sizeof(HQentry<Imgidx, Pixel>));
	}
	HierarHeapQueue_HEQ(Imgidx *dhist, _uint32 *histeqmap_in, Imgidx numlevels_in, double a_in, Imgidx size)
	{
		initHQ(dhist, histeqmap_in, numlevels_in, size, a_in, 12);
	}
	HierarHeapQueue_HEQ(Imgidx *dhist, _uint32 *histeqmap_in, Imgidx numlevels_in, Imgidx size, double a_in, int listsize)
	{
		initHQ(dhist, histeqmap_in, numlevels_in, size, a_in, listsize);
	}

	~HierarHeapQueue_HEQ()
	{
		Free(list);
		Free(qsizes);
		Free(histeqmap);
		Free(storage_cursize);

		for(int level = 0;level < numlevels;level++)
			if(hqueue[level]) delete hqueue[level];
		Free(hqueue);

		for(int level = thr_hqueue;level < numlevels;level++)
			if(storage[level]) Free(storage[level]);
		Free(storage + thr_hqueue);
	}

	inline void start_pushes() { emptytop = 1; }
	inline Pixel get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline Pixel top_alpha() { return list[0].alpha; }
	inline void push_1stitem(Imgidx idx, Pixel alpha)
	{
		list[0].pidx = idx;
		list[0].alpha = alpha;
		curSize_list++;
	}

	inline void end_pushes(_uint8 *isVisited)
	{
		if(emptytop)
			pop(isVisited);
	}

	void push(Imgidx idx, Pixel alpha)
	{
#if QUEUE_DEBUG
		printf("push: %d at %.2f\n", (int)idx, log2((double)alpha));
		//cout << "push: " << idx << " at " << (int)alpha << endl;
#endif
		//printf("Q: push %d at level %d\n", (int)idx, (int)alpha);

//#if TRACK_QUEUEING
		//tmp
//		f << '0' << '\n' << idx << endl;
//#endif

		if(emptytop && alpha < list[0].alpha)
		{
			emptytop = 0;
			list[0].pidx = idx;
			list[0].alpha = alpha;
			return;
		}

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;

		bool push2list = (queue_minlev < curthr) ?  alpha < hqueue[queue_minlev]->top_alpha()
																						 :  (Imgidx)histeqmap[(int)(a * log2(1 + (double)alpha))] < queue_minlev;

		if (push2list)
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				int i;
				for (i = curSize_list; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
				curSize_list++;
			}
			else if (alpha < list[curSize_list].alpha)// push to the full list
			{
				push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

				int i;
				for (i = curSize_list - 1; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
			}
			else
				push_queue(idx, alpha); // push to the queue
		}
		else
			push_queue(idx, alpha); // push to the queue


#if QUEUE_DEBUG
		printf("List(%d): ", (int)(curSize_list + 1));
		for(int i = 0;i <= curSize_list;i++)
		{
			printf("%d-%.3f ", (int)list[i].pidx, log2((double)list[i].alpha));
		}
		printf("\n");
#endif

			//qtime += get_cpu_time() - tt; //tmp
	}

	void push_queue(Imgidx idx, Pixel alpha)
	{
		int level = (int)histeqmap[(int)(a * log2(1 + (double)alpha))];

		//hidx = (int)(log2(1 + pow((double)dimg[dimgidx++],a)));
		if(level < queue_minlev)
			queue_minlev = level;

		if(level < curthr)
			hqueue[level]->push(idx, alpha);
		else
		{
			Imgidx cur = storage_cursize[level]++;
			storage[level][cur].pidx = idx;
			storage[level][cur].alpha = alpha;
		}
	}

	Imgidx pop(_uint8 *isVisited)
	{
		Imgidx ret = top();
		//double tt = get_cpu_time(); //tmp

		_int8 i;
		//Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;
#if QUEUE_DEBUG
		printf("pop: %d at %.2f\n", (int)list[0].pidx, log2((double)list[0].alpha));
		//cout << "pop: " << list[0].pidx << " at " << log2((double)list[0].alpha) << endl;
#endif
		//tmp

		//printf("Q: pop %d at level %d\n", (int)list[0].pidx, (int)list[0].alpha);

//#if TRACK_QUEUEING
//		f << '1' << '\n' << list[0] << endl;
//#endif

		if (curSize_list == 0)
		{
#if QUEUE_DEBUG
			if(queue_minlev == numlevels) //this should never happen
			{
				cout << "Error on HybridQueue.h: Trying to empty a queue." << endl;
			}
#endif
			while(!check_queue_level(isVisited))
				queue_minlev++;
			list[0].pidx = hqueue[queue_minlev]->top();
			list[0].alpha = hqueue[queue_minlev]->top_alpha();

			pop_queue(isVisited);
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if QUEUE_DEBUG
		printf("List(%d): ", (int)(curSize_list + 1));
		for(i = 0;i <= curSize_list;i++)
		{
			printf("%d-%.3f ", (int)list[i].pidx, log2((double)list[i].alpha));
		}
		printf("\n");
#endif



#if TRACK_QUEUEING
		//f << list[0] << endl;
#endif

		//qtime += get_cpu_time() - tt; //tmp
		return ret;
	}

	int check_queue_level(_uint8 *isVisited)
	{
		if(queue_minlev < curthr)
			return hqueue[queue_minlev]->get_cursize();
		else
		{
			//while(curthr < queue_minlev)
			//{
			//	printf("Piep ");
			//	hqueue[curthr] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[curthr]);

			//	Free(storage[curthr]);
			//	storage[curthr] = 0;
			//	curthr++;
			//}
			curthr++;

			hqueue[queue_minlev] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[queue_minlev]);

			HQentry<Imgidx, Pixel>* store = storage[queue_minlev];
			Imgidx cur = storage_cursize[queue_minlev];
			HeapQueue_naive_quad<Imgidx, Pixel> *pQ = hqueue[queue_minlev];
			for(Imgidx p = 0;p < cur;p++)
			{
				if(!isVisited[store[p].pidx])
					pQ->push(store[p].pidx, store[p].alpha);
			}

			Free(storage[queue_minlev]);
			storage[queue_minlev] = 0;

			return pQ->get_cursize();
		}
	}

	void pop_queue(_uint8 *isVisited)
	{
		hqueue[queue_minlev]->pop(); 

		if(!hqueue[queue_minlev]->get_cursize())
		{
			do
			{
				queue_minlev++;
			}while(queue_minlev < numlevels && !check_queue_level(isVisited));
		}
	}
};

template<class Imgidx, class Pixel>
class HierarHeapQueue
{
	HeapQueue_naive_quad<Imgidx, Pixel> **hqueue;
	HQentry<Imgidx, Pixel> **storage;
	Imgidx *storage_cursize;
	Imgidx *qsizes;

	Imgidx thr_hqueue, curthr, numlevels;
	double a;
	Imgidx queue_minlev;

	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;

public:
	#if TRACK_QUEUEING
		Imgidx *in_size;
		ofstream f;
		Imgidx numproc;
		Imgidx numqueue;
	#endif

	HierarHeapQueue(Imgidx *dhist, Imgidx numlevels_in, Imgidx size, double a_in, int listsize, int connectivity, double r)
	{
		#if TRACK_QUEUEING
				f.open("../../groupmeeting/next/queuelog.dat", std::ofstream::out);
				f << size << endl;
				numproc = 0;
				numqueue = 0;
		#endif

		this->numlevels = numlevels_in;
		this->a = a_in;
		this->queue_minlev = numlevels;

		Imgidx cumsum = 0;
		qsizes = dhist; //do not free dhist outside
		if(r >= 1)
		{
			thr_hqueue = curthr = numlevels;
			hqueue = (HeapQueue_naive_quad<Imgidx, Pixel>**)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Imgidx, Pixel>*));
			for(int level = 0;level < thr_hqueue;level++)
				hqueue[level] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[level]);
			storage = 0;
			storage_cursize = 0;
		}
		else
		{
			storage_cursize = (Imgidx*)Calloc(numlevels * sizeof(Imgidx));
			Imgidx thr_nonredundantnodes = (Imgidx)(size * r);
			for(int level = 0;level < numlevels;level++)
			{
				cumsum += qsizes[level];
				if(cumsum > thr_nonredundantnodes)
				{
					thr_hqueue = curthr = level;
					break;
				}
			}

			hqueue = (HeapQueue_naive_quad<Imgidx, Pixel>**)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Imgidx, Pixel>*));
			for(int level = 0;level < thr_hqueue;level++)
				hqueue[level] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[level]);
			//hqueue[numlevels] = new HeapQueue_naive_quad<Imgidx, Pixel>(1);
			//hqueue[numlevels]->push(0,(Pixel)-1);

			storage = (HQentry<Imgidx, Pixel>**)Calloc((numlevels - thr_hqueue) * sizeof(HQentry<Imgidx, Pixel>*));
			storage -= thr_hqueue;
			for(int level = thr_hqueue;level < numlevels;level++)
				storage[level] = (HQentry<Imgidx, Pixel>*)Malloc(qsizes[level] * sizeof(HQentry<Imgidx, Pixel>));
		}
	}

	~HierarHeapQueue()
	{
		Free(qsizes);

		for(int level = 0;level < numlevels;level++)
			if(hqueue[level]) delete hqueue[level];
		Free(hqueue);

		if(storage)
		{
			Free(storage_cursize);
			for(int level = thr_hqueue;level < numlevels;level++)
				if(storage[level]) Free(storage[level]);
			Free(storage + thr_hqueue);
		}
	}

	inline Imgidx top() { return hqueue[queue_minlev]->top(); }
	inline Pixel top_alpha() { return hqueue[queue_minlev]->top_alpha(); }
	inline void push_1stitem(Imgidx idx, Pixel alpha)
	{
		push_queue(idx, alpha);

		#if TRACK_QUEUEING
			//tmp
			f << '0' << '\n' << idx << endl << alpha << endl;
			numproc++;
		#endif
	}

	inline void end_pushes(_uint8* isVisited)
	{
		pop(isVisited);
	}

	void push(Imgidx idx, Pixel alpha)
	{
#if QUEUE_DEBUG
		printf("push: %d at %.2f\n", (int)idx, log2((double)alpha));
		//cout << "push: " << idx << " at " << (int)alpha << endl;
#endif
		//printf("Q: push %d at level %d\n", (int)idx, (int)alpha);

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		push_queue(idx, alpha); // push to the queue
	}

	void push_queue(Imgidx idx, Pixel alpha)
	{

	#if TRACK_QUEUEING
		numqueue++;
	#endif
		int level = (int)(a * log2(1 + (double)alpha));

		//hidx = (int)(log2(1 + pow((double)dimg[dimgidx++],a)));
		if(level < queue_minlev)
			queue_minlev = level;

		if(level < curthr)
			hqueue[level]->push(idx, alpha);
		else
		{
			Imgidx cur = storage_cursize[level]++;
			storage[level][cur].pidx = idx;
			storage[level][cur].alpha = alpha;
		}
	}


	#if TRACK_QUEUEING
	inline void mark(int m)
	{
		if(m == 0)
			f << m << endl;
		else
		{
			f << '1' << endl << m << endl;
		}
	}
	inline void redundant(Pixel alpha)
	{
		f << '2' << endl << alpha << endl;
	}
	#endif

	Imgidx pop(_uint8 *isVisited)
	{
		Imgidx ret = top();
		//double tt = get_cpu_time(); //tmp

		//Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;

#if QUEUE_DEBUG
		printf("pop: %d at %.2f\n", (int)top_alpha(), log2((double)top_alpha()));
#endif
		pop_queue(isVisited);

#if TRACK_QUEUEING
		//f << list[0] << endl;
#endif

		//qtime += get_cpu_time() - tt; //tmp
		return ret;
	}

	int check_queue_level(_uint8* isVisited)
	{
		if(queue_minlev < curthr)
			return hqueue[queue_minlev]->get_cursize();
		else
		{
			while(curthr < queue_minlev)
			{
				hqueue[curthr] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[curthr]);

				Free(storage[curthr]);
				storage[curthr] = 0;
				curthr++;
			}
			curthr++;

			hqueue[queue_minlev] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[queue_minlev]);

			HQentry<Imgidx, Pixel>* store = storage[queue_minlev];
			Imgidx cur = storage_cursize[queue_minlev];
			HeapQueue_naive_quad<Imgidx, Pixel> *pQ = hqueue[queue_minlev];
			for(Imgidx p = 0;p < cur;p++)
			{
				if(!isVisited[store[p].pidx])
					pQ->push(store[p].pidx, store[p].alpha);
			}

			Free(storage[queue_minlev]);
			storage[queue_minlev] = 0;

			return pQ->get_cursize();
		}
	}

	void pop_queue(_uint8* isVisited)
	{
		hqueue[queue_minlev]->pop();

		if(!hqueue[queue_minlev]->get_cursize())
		{
			do
			{
				queue_minlev++;
			}while(queue_minlev < numlevels && !check_queue_level(isVisited));
		}
	}
};


template<class Imgidx, class Pixel>
class HierarHeapQueue_cache
{
	HQentry<Imgidx, Pixel> *list;
	HeapQueue_naive_quad<Imgidx, Pixel> **hqueue;
	HQentry<Imgidx, Pixel> **storage;
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

	#if TRACK_QUEUEING
		Imgidx *in_size;
		ofstream f;
		Imgidx numproc;
		Imgidx numqueue;
	#endif

	void initHQ(Imgidx *dhist, Imgidx numlevels_in, Imgidx size, double a_in, int listsize, int connectivity, double r)
	{
		/*		cnt = 0;//tmp*/
		//Imgidx i;
		totalsize = size;

#if PROFILE
		t0 = get_cpu_time();
		tconv = tcache = tqueue = 0;
		num_cache =
		num_cache_ovfl =
		num_hq =
		num_store =
		num_conv = 0;
#endif

#if TRACK_QUEUEING
		f.open("../../groupmeeting/next/queuelog.dat", std::ofstream::out);
		f << size << endl;
		numproc = 0;
		numqueue = 0;
#endif

		list = (HQentry<Imgidx, Pixel>*)Malloc((listsize + 1) * sizeof(HQentry<Imgidx, Pixel>));
		list = list - 1;
		list[-1].alpha = 0;
		maxSize_list = listsize - 1;
		curSize_list = -1;

		this->numlevels = numlevels_in;
		this->a = a_in;
		this->queue_minlev = numlevels;

		Imgidx cumsum = 0;
		qsizes = dhist; //do not free dhist outside
		if(r >= 1)
		{
			thr_hqueue = curthr = numlevels;
			hqueue = (HeapQueue_naive_quad<Imgidx, Pixel>**)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Imgidx, Pixel>*));
			for(int level = 0;level < thr_hqueue;level++)
				hqueue[level] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[level]);
			storage = 0;
			storage_cursize = 0;
		}
		else
		{
			storage_cursize = (Imgidx*)Calloc(numlevels * sizeof(Imgidx));
			Imgidx thr_nonredundantnodes = (Imgidx)(size * r);
			for(int level = 0;level < numlevels;level++)
			{
				cumsum += qsizes[level];
				if(cumsum > thr_nonredundantnodes)
				{
					thr_hqueue = curthr = level;
					break;
				}
			}

			hqueue = (HeapQueue_naive_quad<Imgidx, Pixel>**)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Imgidx, Pixel>*));
			for(int level = 0;level < thr_hqueue;level++)
				hqueue[level] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[level]);
			//hqueue[numlevels] = new HeapQueue_naive_quad<Imgidx, Pixel>(1);
			//hqueue[numlevels]->push(0,(Pixel)-1);

			storage = (HQentry<Imgidx, Pixel>**)Calloc((numlevels - thr_hqueue) * sizeof(HQentry<Imgidx, Pixel>*));
			storage -= thr_hqueue;
			for(int level = thr_hqueue;level < numlevels;level++)
				storage[level] = (HQentry<Imgidx, Pixel>*)Malloc(qsizes[level] * sizeof(HQentry<Imgidx, Pixel>));
		}		
	}

	HierarHeapQueue_cache(Imgidx *dhist, Imgidx numlevels_in, double a_in, Imgidx size, Imgidx connectivity = 4, double r = 0.2)
	{
		initHQ(dhist, numlevels_in, size, a_in, 12, (int)connectivity, r);
	}
	HierarHeapQueue_cache(Imgidx *dhist, Imgidx numlevels_in, Imgidx size, double a_in, int listsize, Imgidx connectivity = 4, double r = 0.2)
	{
		initHQ(dhist, numlevels_in, size, a_in, listsize, (int)connectivity, r);
	}

	~HierarHeapQueue_cache()
	{
		Free(list + 1);
		Free(qsizes);

		for(int level = 0;level < numlevels;level++)
			if(hqueue[level]) delete hqueue[level];
		Free(hqueue);

		if(storage)
		{
			Free(storage_cursize);
			for(int level = thr_hqueue;level < numlevels;level++)
				if(storage[level]) Free(storage[level]);
			Free(storage + thr_hqueue);
		}

#if PROFILE
		double t1 = get_cpu_time() - t0;
		double sz = (double)totalsize;
		printf("HQ Profile - Total: %f, Cache: %f, Queue: %f, Conv: %f\n", t1, tcache, tqueue, tconv);
		printf("Cached: %d(%f), C. ovfl: %d(%f), Heap Queued: %d(%f), Stored: %d(%f), Converted: %d(%f), \n",
		 (int)num_cache, (double)num_cache / sz, 
		 (int)num_cache_ovfl, (double)num_cache_ovfl / sz, 
		 (int)num_hq, (double)num_hq / sz, 
		 (int)num_store, (double)num_store / sz, 
		 (int)num_conv, (double)num_conv / sz);
		printf("Initial: %d queues + %d storages (%d total) End: %d queues + %d storages (%d total)\n",
		 (int)thr_hqueue, (int)(numlevels - thr_hqueue), (int)numlevels, (int)curthr, (int)(numlevels - curthr), (int)numlevels);

#endif
	}

	inline void start_pushes()
	{
	#if TRACK_QUEUEING
			f << '1' << endl;
	#endif
		 emptytop = 1;
	}
	inline Pixel get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline Pixel top_alpha() { return list[0].alpha; }
	inline void push_1stitem(Imgidx idx, Pixel alpha)
	{
		list[0].pidx = idx;
		list[0].alpha = alpha;
		curSize_list++;

		#if TRACK_QUEUEING
				//tmp
				f << '0' << '\n' << idx << endl << alpha << endl;
				numproc++;
		#endif
	}

	inline void end_pushes(_uint8 *isVisited)
	{
		if(emptytop)
			pop(isVisited);
	}

	void push(Imgidx idx, Pixel alpha)
	{
#if QUEUE_DEBUG
		printf("push: %d at %.2f\n", (int)idx, log2((double)alpha));
		//cout << "push: " << idx << " at " << (int)alpha << endl;
#endif
		//printf("Q: push %d at level %d\n", (int)idx, (int)alpha);

#if TRACK_QUEUEING
		//tmp
		numproc++;
		f << '0' << '\n' << idx << endl << alpha << endl;
#endif

		if(emptytop && alpha < list[0].alpha)
		{
#if PROFILE
			num_cache++;
#endif
			emptytop = 0;
			list[0].pidx = idx;
			list[0].alpha = alpha;
			return;
		}

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;

		bool push2list = (queue_minlev < curthr) ?  alpha < hqueue[queue_minlev]->top_alpha()
																						 :  (int)(a * log2(1 + (double)alpha)) < queue_minlev;

		if (push2list)
		{
#if PROFILE
		num_cache++;
		double t1 = get_cpu_time(), t2, tq = 0;
#endif
			if (curSize_list < maxSize_list) //spare room in the list
			{
				int i;
				for (i = curSize_list; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
				curSize_list++;
			}
			else if (alpha < list[curSize_list].alpha)// push to the full list
			{
#if PROFILE
				num_cache_ovfl++;
				t2 = get_cpu_time();
#endif
				push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

#if PROFILE
				tq = get_cpu_time() - t2;
#endif
				int i;
				for (i = curSize_list - 1; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
			}
			else
			{
#if PROFILE
				num_cache_ovfl++;
				t2 = get_cpu_time();
#endif
				push_queue(idx, alpha); // push to the queue
#if PROFILE
				tq = get_cpu_time() - t2;
#endif
			}

#if PROFILE
			tcache += get_cpu_time() - t1 - tq;
			tqueue += tq;
#endif
		}
		else
		{
#if PROFILE
			double t1 = get_cpu_time();
#endif
			push_queue(idx, alpha); // push to the queue
#if PROFILE
			tqueue += get_cpu_time() - t1;
#endif
		}


#if QUEUE_DEBUG
		printf("List(%d): ", (int)(curSize_list + 1));
		for(int i = 0;i <= curSize_list;i++)
		{
			printf("%d-%.3f ", (int)list[i].pidx, log2((double)list[i].alpha));
		}
		printf("\n");
#endif

			//qtime += get_cpu_time() - tt; //tmp
	}

	void push_queue(Imgidx idx, Pixel alpha)
	{

	#if TRACK_QUEUEING
		numqueue++;
	#endif
		int level = (int)(a * log2(1 + (double)alpha));

		//hidx = (int)(log2(1 + pow((double)dimg[dimgidx++],a)));
		if(level < queue_minlev)
			queue_minlev = level;

		if(level < curthr)
		{
#if PROFILE
			num_hq++;
#endif
			hqueue[level]->push(idx, alpha);
		}
		else
		{
#if PROFILE
			num_store++;
#endif
			Imgidx cur = storage_cursize[level]++;
			storage[level][cur].pidx = idx;
			storage[level][cur].alpha = alpha;
		}
	}


	#if TRACK_QUEUEING
	inline void mark(int m)
	{
		if(m == 0)
			f << m << endl;
		else
		{
			f << '1' << endl << m << endl;
		}
	}
	inline void redundant(Pixel alpha)
	{
		f << '2' << endl << alpha << endl;
	}
	#endif

	Imgidx pop(_uint8 *isVisited)
	{
		Imgidx ret = top();
		//double tt = get_cpu_time(); //tmp

		_int8 i;
		//Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;

#if QUEUE_DEBUG
		printf("pop: %d at %.2f\n", (int)list[0].pidx, log2((double)list[0].alpha));
		//cout << "pop: " << list[0].pidx << " at " << log2((double)list[0].alpha) << endl;
#endif
		//tmp

		//printf("Q: pop %d at level %d\n", (int)list[0].pidx, (int)list[0].alpha);
		if (curSize_list == 0)
		{
#if QUEUE_DEBUG
			if(queue_minlev == numlevels) //this should never happen
			{
				cout << "Error on HybridQueue.h: Trying to empty a queue." << endl;
			}
#endif
			while(!check_queue_level(isVisited))
				queue_minlev++;
			list[0].pidx = hqueue[queue_minlev]->top();
			list[0].alpha = hqueue[queue_minlev]->top_alpha();

			pop_queue(isVisited);
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if QUEUE_DEBUG
		printf("List(%d): ", (int)(curSize_list + 1));
		for(i = 0;i <= curSize_list;i++)
		{
			printf("%d-%.3f ", (int)list[i].pidx, log2((double)list[i].alpha));
		}
		printf("\n");
#endif



#if TRACK_QUEUEING
		//f << list[0] << endl;
#endif

		//qtime += get_cpu_time() - tt; //tmp
		return ret;
	}

	int check_queue_level(_uint8 *isVisited)
	{
		if(queue_minlev < curthr)
			return hqueue[queue_minlev]->get_cursize();
		else
		{			
#if PROFILE
			double t1 = get_cpu_time();
#endif
			while(curthr < queue_minlev)
			{
				hqueue[curthr] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[curthr]);

				Free(storage[curthr]);
				storage[curthr] = 0;
				curthr++;
			}
			curthr++;

			hqueue[queue_minlev] = new HeapQueue_naive_quad<Imgidx, Pixel>(qsizes[queue_minlev]);

			HQentry<Imgidx, Pixel>* store = storage[queue_minlev];
			Imgidx cur = storage_cursize[queue_minlev];
			HeapQueue_naive_quad<Imgidx, Pixel> *pQ = hqueue[queue_minlev];
			for(Imgidx p = 0;p < cur;p++)
			{
				if(!isVisited[store[p].pidx])
					pQ->push(store[p].pidx, store[p].alpha);
			}
#if PROFILE
			num_conv += cur;
#endif

			Free(storage[queue_minlev]);
			storage[queue_minlev] = 0;
	
#if PROFILE
			tconv += get_cpu_time() - t1;
#endif
			return pQ->get_cursize();
		}
	}

	void pop_queue(_uint8 *isVisited)
	{
			
#if PROFILE
		double t1 = get_cpu_time();
#endif
		hqueue[queue_minlev]->pop(); 
			
#if PROFILE
		tqueue += get_cpu_time() - t1;
#endif

		if(!hqueue[queue_minlev]->get_cursize())
		{
			do
			{
				queue_minlev++;
			}while(queue_minlev < numlevels && !check_queue_level(isVisited));
		}
	}
};

template<class Imgidx, class Pixel>//, class Qidx>
class Cache_Heapqueue
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	HQentry<Imgidx, Pixel> *list;
	HeapQueue_naive<Imgidx, Pixel> *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
	//	_int32 cnt;
	void initHQ(Imgidx size, size_t listsize)
	{
		/*		cnt = 0;//tmp*/
		//Imgidx i;
		this->maxSize_queue = size;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		//
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
		hqueue = new HeapQueue_naive<Imgidx, Pixel>(size);
		list = (HQentry<Imgidx, Pixel>*)Malloc((listsize + 1) * sizeof(HQentry<Imgidx, Pixel>));
		list[0].pidx = 0;
		list[0].alpha = 0;
		list++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		//minidx_queue = size;


		//tmp
//#if TRACK_QUEUEING
//		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
//		f << -1 << '\n' << size << endl;
//#endif
		qtime = 0;//tmp
	}
public:

	double qtime;//tmp

	Cache_Heapqueue(Imgidx size)
	{
		initHQ(size, 12);
	}
	Cache_Heapqueue(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}

	inline void start_pushes() { emptytop = 1; }
	inline Imgidx get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline Pixel top_alpha() { return list[0].alpha; }
	inline void push_1stitem(Imgidx idx, Pixel alpha)
	{
		list[0].pidx = idx;
		list[0].alpha = alpha;
		curSize_list++;
	}

	inline void end_pushes()
	{
		if(emptytop)
			pop();
	}

	inline void push(Imgidx idx, Pixel alpha)
	{
		//double tt = get_cpu_time(); //tmp

		//MinList1<Imgidx> *p, *q;
		_int16 i;
#if QUEUE_DEBUG
		cout << "push: " << idx << " at " << (int)alpha << endl;
#endif
		//printf("Q: push %d at level %d\n", (int)idx, (int)alpha);

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		if(emptytop && alpha < list[0].alpha)
		{
			emptytop = 0;
			list[0].pidx = idx;
			list[0].alpha = alpha;
			return;
		}

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if (alpha < hqueue->top_alpha())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
				curSize_list++;
			}
			else if (alpha < list[curSize_list].alpha)// push to the full list
			{
				push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

				for (i = curSize_list - 1; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
			}
			else
				push_queue(idx, alpha); // push to the queue
		}
		else
			push_queue(idx, alpha); // push to the queue


			//qtime += get_cpu_time() - tt; //tmp
	}
	inline void push_queue(Imgidx idx, Pixel alpha)
	{
		hqueue->push(idx, alpha);
	}
	inline Imgidx pop()
	{
		Imgidx ret = top();
		//double tt = get_cpu_time(); //tmp

		_int8 i;
		//Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;
#if QUEUE_DEBUG
		cout << "pop: " << list[0].pidx << " at " << (int)list[0].alpha << endl;
#endif

		//tmp

		//printf("Q: pop %d at level %d\n", (int)list[0].pidx, (int)list[0].alpha);

		if (curSize_list == 0)
		{
			list[0].pidx = hqueue->top();
			list[0].alpha = hqueue->top_alpha();

			pop_queue();
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		//f << list[0] << endl;
#endif

		//qtime += get_cpu_time() - tt; //tmp
		return ret;
	}
	inline void pop_queue()
	{
		hqueue->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	~Cache_Heapqueue()
	{
		delete hqueue;
		Free(list - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};

template<class Imgidx, class Pixel>//, class Qidx>
class Cache_Quad_Heapqueue
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	HQentry<Imgidx, Pixel> *list;
	HeapQueue_naive_quad<Imgidx, Pixel> *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
	//	_int32 cnt;
	void initHQ(Imgidx size, size_t listsize)
	{
		/*		cnt = 0;//tmp*/
		//Imgidx i;
		this->maxSize_queue = size;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		//
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
		hqueue = new HeapQueue_naive_quad<Imgidx, Pixel>(size);
		list = (HQentry<Imgidx, Pixel>*)Malloc((listsize + 1) * sizeof(HQentry<Imgidx, Pixel>));
		list[0].pidx = 0;
		list[0].alpha = 0;
		list++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		//minidx_queue = size;


		//tmp
//#if TRACK_QUEUEING
//		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
//		f << -1 << '\n' << size << endl;
//#endif
		qtime = 0;//tmp
	}
public:

	double qtime;//tmp

	Cache_Quad_Heapqueue(Imgidx size)
	{
		initHQ(size, 12);
	}
	Cache_Quad_Heapqueue(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}

	inline void start_pushes() { emptytop = 1; }
	inline Pixel get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline Pixel top_alpha() { return list[0].alpha; }
	inline void push_1stitem(Imgidx idx, Pixel alpha)
	{
		list[0].pidx = idx;
		list[0].alpha = alpha;
		curSize_list++;
	}

	inline void end_pushes()
	{
		if(emptytop)
			pop();
	}

	inline void push(Imgidx idx, Pixel alpha)
	{
		//double tt = get_cpu_time(); //tmp

		//MinList1<Imgidx> *p, *q;
		_int16 i;
#if QUEUE_DEBUG
		cout << "push: " << idx << " at " << (int)alpha << endl;
#endif
		//printf("Q: push %d at level %d\n", (int)idx, (int)alpha);

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		if(emptytop && alpha < list[0].alpha)
		{
			emptytop = 0;
			list[0].pidx = idx;
			list[0].alpha = alpha;
			return;
		}

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if (alpha < hqueue->top_alpha())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
				curSize_list++;
			}
			else if (alpha < list[curSize_list].alpha)// push to the full list
			{
				push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

				for (i = curSize_list - 1; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
			}
			else
				push_queue(idx, alpha); // push to the queue
		}
		else
			push_queue(idx, alpha); // push to the queue


			//qtime += get_cpu_time() - tt; //tmp
	}
	inline void push_queue(Imgidx idx, Pixel alpha)
	{
		hqueue->push(idx, alpha);
	}
	inline Imgidx pop()
	{
		Imgidx ret = top();
		//double tt = get_cpu_time(); //tmp

		_int8 i;
		//Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;
#if QUEUE_DEBUG
		cout << "pop: " << list[0].pidx << " at " << (int)list[0].alpha << endl;
#endif

		//tmp

		//printf("Q: pop %d at level %d\n", (int)list[0].pidx, (int)list[0].alpha);

#if TRACK_QUEUEING
		//f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 0)
		{
			list[0].pidx = hqueue->top();
			list[0].alpha = hqueue->top_alpha();

			pop_queue();
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		//f << list[0] << endl;
#endif

		//qtime += get_cpu_time() - tt; //tmp
		return ret;
	}
	inline void pop_queue()
	{
		hqueue->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	~Cache_Quad_Heapqueue()
	{
		delete hqueue;
		Free(list - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};

template<class Imgidx, class Pixel>//, class Qidx>
class CirCache_Hierqueue
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	HQentry<Imgidx, _int32> *list;
	HierarQueue<Imgidx> *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list, liststart, mask;
	Imgidx maxSize_queue, mask_field;
	int emptytop;

#if TRACK_QUEUEING
	Imgidx *in_size;
	ofstream f;
#endif

	//	_int32 cnt;
	//void initHQ(Imgidx size, size_t listsize)
	void initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
	{
		int shamt = 0;
		for(int lsize = listsize;lsize;lsize>>=1)
			shamt++;
		listsize = 1 << (shamt - 1);
		mask = listsize - 1;
		/*		cnt = 0;//tmp*/
		//Imgidx i;
		this->maxSize_queue = qsize_in;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		//
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
		hqueue = new HierarQueue<Imgidx>(qsize_in, dhist, numlevels);
		list = (HQentry<Imgidx, _int32>*)Malloc(listsize * sizeof(HQentry<Imgidx, _int32>));
		maxSize_list = listsize;
		curSize_list = 0;
		liststart = 0;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		//minidx_queue = size;


		//tmp
//#if TRACK_QUEUEING
//		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
//		f << -1 << '\n' << size << endl;
//#endif
		qtime = 0;//tmp
	}
public:

	double qtime;//tmp

	CirCache_Hierqueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
	{
		initHQ(qsize_in, dhist, numlevels, 16);
	}
	CirCache_Hierqueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
	{
		initHQ(qsize_in, dhist, numlevels, listsize);
	}

	inline void start_pushes() { emptytop = 1; }
	inline _int32 get_minlev() { return list[liststart].alpha; }
	inline Imgidx top() { return list[liststart].pidx; }
	inline _int32 top_alpha() { return list[liststart].alpha; }
	inline void push_1stitem(Imgidx idx, _int32 alpha)
	{
		list[0].pidx = idx;
		list[0].alpha = alpha;
		curSize_list++;
	}

	inline void end_pushes()
	{
		if(emptytop) pop();

	}

	inline void push(Imgidx idx, _int32 alpha)
	{
		//double tt = get_cpu_time(); //tmp

		//MinList1<Imgidx> *p, *q;
		_int16 i, j, k;
#if QUEUE_DEBUG
		cout << "chierQ - push: " << idx << " at " << (int)alpha << endl;
#endif
		//printf("Q: push %d at level %d\n", (int)idx, (int)alpha);

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		if(emptytop && alpha < list[liststart].alpha)
		{
			emptytop = 0;
			list[liststart].pidx = idx;
			list[liststart].alpha = alpha;
			return;
		}

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if ((_int64)alpha < hqueue->get_minlev())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				j = (liststart - 1) & mask;
				for (i = (liststart + curSize_list - 1) & mask; i != j && alpha < list[i].alpha; i = (i - 1) & mask)
					list[(i + 1) & mask] = list[i];
				list[(i + 1) & mask].pidx = idx;
				list[(i + 1) & mask].alpha = alpha;
				curSize_list++;
			}
			else if (alpha < list[j = ((liststart + curSize_list - 1) & mask)].alpha)// push to the full list
			{
				push_queue(list[j].pidx, list[j].alpha);

				k = (liststart - 1) & mask;
				for (i = (j - 1) & mask; i != k && alpha < list[i].alpha; i = (i - 1) & mask)
					list[(i + 1) & mask] = list[i];
				list[(i + 1) & mask].pidx = idx;
				list[(i + 1) & mask].alpha = alpha;
			}
			else
				push_queue(idx, alpha); // push to the queue
		}
		else
			push_queue(idx, alpha); // push to the queue


#if QUEUE_DEBUG
		printf("List(%d): ", (int)curSize_list);
		i = liststart;
		for(j = 0;j < curSize_list;j++)
		{
			printf("(%d, %d): %d-%d ", (int)liststart, (int)curSize_list, (int)list[i].pidx, (int)list[i].alpha);
			i = (i + 1) & mask;
		}
		printf("\n");
#endif


			//qtime += get_cpu_time() - tt; //tmp
	}
	inline void push_queue(Imgidx idx, _int32 alpha)
	{
		hqueue->push(idx, alpha);
	}
	inline Imgidx pop()
	{
		Imgidx ret = top();
		//double tt = get_cpu_time(); //tmp

		//_int8 i;
		//Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;
#if QUEUE_DEBUG
		cout << "chierQ - pop: " << list[liststart].pidx << " at " << (int)list[liststart].alpha << endl;
#endif

		//tmp

		//printf("Q: pop %d at level %d\n", (int)list[0].pidx, (int)list[0].alpha);

#if TRACK_QUEUEING
		f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 1)
		{
			list[liststart].pidx = hqueue->top();
			list[liststart].alpha = hqueue->get_minlev();

			pop_queue();
		}
		else
		{
			liststart = (liststart + 1) & mask;
			//for (i = 0; i < curSize_list; i++)
			//	list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		//f << list[0] << endl;
#endif

		//qtime += get_cpu_time() - tt; //tmp
		return ret;
	}
	inline void pop_queue()
	{
		hqueue->pop();
		hqueue->find_minlev();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	~CirCache_Hierqueue()
	{
		delete hqueue;
		Free(list);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};

template<class Imgidx, class Pixel>//, class Qidx>
class HierarQueueCache
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	HQentry<Imgidx, _int32> *list;
	HierarQueue<Imgidx> *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif

	//	_int32 cnt;
	//void initHQ(Imgidx size, size_t listsize)
	void initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize = LISTSIZE_DEFAULT)
	{
		/*		cnt = 0;//tmp*/
		//Imgidx i;
		this->maxSize_queue = qsize_in;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		//
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
		hqueue = new HierarQueue<Imgidx>(qsize_in, dhist, numlevels);
		list = (HQentry<Imgidx, _int32>*)Malloc((listsize + 1) * sizeof(HQentry<Imgidx, _int32>));
		list[0].pidx = 0;
		list[0].alpha = 0;
		list++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		//minidx_queue = size;


		//tmp
//#if TRACK_QUEUEING
//		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
//		f << -1 << '\n' << size << endl;
//#endif
		qtime = 0;//tmp
	}
public:

	double qtime;//tmp

	HierarQueueCache(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
	{
		initHQ(qsize_in, dhist, numlevels, 12);
	}
	HierarQueueCache(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
	{
		initHQ(qsize_in, dhist, numlevels, listsize);
	}

	inline void start_pushes() { emptytop = 1; }
	inline _int32 get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline _int32 top_alpha() { return list[0].alpha; }
	inline void push_1stitem(Imgidx idx, _int32 alpha)
	{
		list[0].pidx = idx;
		list[0].alpha = alpha;
		curSize_list++;
	}

	inline void end_pushes()
	{
		if(emptytop) pop();

	}

	inline void push(Imgidx idx, _int32 alpha)
	{
		//double tt = get_cpu_time(); //tmp

		//MinList1<Imgidx> *p, *q;
		_int16 i;
#if QUEUE_DEBUG
//		cout << "chierQ - push: " << idx << " at " << (int)alpha << endl;
#endif
		//printf("Q: push %d at level %d\n", (int)idx, (int)alpha);

#if TRACK_QUEUEING
		//tmp
//		f << '0' << '\n' << idx << endl;
#endif

		if(emptytop && alpha < list[0].alpha)
		{
			emptytop = 0;
			list[0].pidx = idx;
			list[0].alpha = alpha;
			return;
		}

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if ((_int64)alpha < hqueue->get_minlev())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
				curSize_list++;
			}
			else if (alpha < list[curSize_list].alpha)// push to the full list
			{
				push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

				for (i = curSize_list - 1; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
			}
			else
				push_queue(idx, alpha); // push to the queue
		}
		else
			push_queue(idx, alpha); // push to the queue


#if QUEUE_DEBUG
		printf("List(%d): ", (int)(curSize_list + 1));
		for(i = 0;i <= curSize_list;i++)
		{
			printf("%d-%d ", (int)list[i].pidx, (int)list[i].alpha);
		}
		printf("\n");
#endif

			//qtime += get_cpu_time() - tt; //tmp
	}
	inline void push_queue(Imgidx idx, _int32 alpha)
	{
		hqueue->push(idx, alpha);
	}
	inline Imgidx pop()
	{
		Imgidx ret = top();
		//double tt = get_cpu_time(); //tmp

		_int8 i;
		//Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;
#if QUEUE_DEBUG
//		cout << "chierQ - pop: " << list[0].pidx << " at " << (int)list[0].alpha << endl;
#endif
		//tmp

		//printf("Q: pop %d at level %d\n", (int)list[0].pidx, (int)list[0].alpha);

#if TRACK_QUEUEING
//		f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 0)
		{
			list[0].pidx = hqueue->top();
			list[0].alpha = hqueue->get_minlev();

			pop_queue();
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		//f << list[0] << endl;
#endif

		//qtime += get_cpu_time() - tt; //tmp
		return ret;
	}
	inline void pop_queue()
	{
		hqueue->pop();
		hqueue->find_minlev();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	~HierarQueueCache()
	{
		delete hqueue;
		Free(list - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};

template<class Imgidx, class Pixel>//, class Qidx>
class Cache_Hierqueue_l1
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	HQentry<Imgidx, _int32> *list;
	HQueue_l1idx<Imgidx> *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif

	//	_int32 cnt;
	//void initHQ(Imgidx size, size_t listsize)
	void initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
	{
		/*		cnt = 0;//tmp*/
		//Imgidx i;
		this->maxSize_queue = qsize_in;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		//
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
		hqueue = new HQueue_l1idx<Imgidx>(qsize_in, dhist, numlevels);
		list = (HQentry<Imgidx, _int32>*)Malloc((listsize + 1) * sizeof(HQentry<Imgidx, _int32>));
		list[0].pidx = 0;
		list[0].alpha = 0;
		list++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		//minidx_queue = size;


		//tmp
//#if TRACK_QUEUEING
//		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
//		f << -1 << '\n' << size << endl;
//#endif
		qtime = 0;//tmp
	}
public:

	double qtime;//tmp

	Cache_Hierqueue_l1(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
	{
		initHQ(qsize_in, dhist, numlevels, 12);
	}
	Cache_Hierqueue_l1(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
	{
		initHQ(qsize_in, dhist, numlevels, listsize);
	}

	inline void start_pushes() { emptytop = 1; }
	inline _int32 get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline _int32 top_alpha() { return list[0].alpha; }
	inline void push_1stitem(Imgidx idx, _int32 alpha)
	{
		list[0].pidx = idx;
		list[0].alpha = alpha;
		curSize_list++;
	}

	inline void end_pushes()
	{
		if(emptytop) pop();

	}

	inline void push(Imgidx idx, _int32 alpha)
	{
		//double tt = get_cpu_time(); //tmp

		//MinList1<Imgidx> *p, *q;
		_int16 i;
#if QUEUE_DEBUG
		cout << "chierQ - push: " << idx << " at " << (int)alpha << endl;
#endif
		//printf("Q: push %d at level %d\n", (int)idx, (int)alpha);

#if TRACK_QUEUEING
		//tmp
		//f << '0' << '\n' << idx << endl;
#endif

		if(emptytop && alpha < list[0].alpha)
		{
			emptytop = 0;
			list[0].pidx = idx;
			list[0].alpha = alpha;
			return;
		}

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if ((_int64)alpha < hqueue->get_minlev())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
				curSize_list++;
			}
			else if (alpha < list[curSize_list].alpha)// push to the full list
			{
				push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

				for (i = curSize_list - 1; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
			}
			else
				push_queue(idx, alpha); // push to the queue
		}
		else
			push_queue(idx, alpha); // push to the queue


			//qtime += get_cpu_time() - tt; //tmp
	}
	inline void push_queue(Imgidx idx, _int32 alpha)
	{
		hqueue->push(idx, alpha);
	}
	inline Imgidx pop()
	{
		Imgidx ret = top();
		//double tt = get_cpu_time(); //tmp

		_int8 i;
		//Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;
#if QUEUE_DEBUG
		cout << "chierQ - pop: " << list[0].pidx << " at " << (int)list[0].alpha << endl;
#endif

		//tmp

		//printf("Q: pop %d at level %d\n", (int)list[0].pidx, (int)list[0].alpha);

#if TRACK_QUEUEING
		//f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 0)
		{
			list[0].pidx = hqueue->top();
			list[0].alpha = hqueue->get_minlev();

			pop_queue();
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		//f << list[0] << endl;
#endif

		//qtime += get_cpu_time() - tt; //tmp
		return ret;
	}
	inline void pop_queue()
	{
		hqueue->pop();
		hqueue->find_minlev();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	~Cache_Hierqueue_l1()
	{
		delete hqueue;
		Free(list - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};

template<class Imgidx, class Pixel>//, class Qidx>
class Cache_Hierqueue_l2
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	HQentry<Imgidx, _int32> *list;
	HQueue_l2idx<Imgidx> *hqueue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;
	int emptytop;

#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif

	//	_int32 cnt;
	//void initHQ(Imgidx size, size_t listsize)
	void initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
	{
		/*		cnt = 0;//tmp*/
		//Imgidx i;
		this->maxSize_queue = qsize_in;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		//
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
		hqueue = new HQueue_l2idx<Imgidx>(qsize_in, dhist, numlevels);
		list = (HQentry<Imgidx, _int32>*)Malloc((listsize + 1) * sizeof(HQentry<Imgidx, _int32>));
		list[0].pidx = 0;
		list[0].alpha = 0;
		list++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		//minidx_queue = size;


		//tmp
//#if TRACK_QUEUEING
//		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
//		f << -1 << '\n' << size << endl;
//#endif
		qtime = 0;//tmp
	}
public:

	double qtime;//tmp

	Cache_Hierqueue_l2(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
	{
		initHQ(qsize_in, dhist, numlevels, 12);
	}
	Cache_Hierqueue_l2(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
	{
		initHQ(qsize_in, dhist, numlevels, listsize);
	}

	inline void start_pushes() { emptytop = 1; }
	inline _int32 get_minlev() { return list[0].alpha; }
	inline Imgidx top() { return list[0].pidx; }
	inline _int32 top_alpha() { return list[0].alpha; }
	inline void push_1stitem(Imgidx idx, _int32 alpha)
	{
		list[0].pidx = idx;
		list[0].alpha = alpha;
		curSize_list++;
	}

	inline void end_pushes()
	{
		if(emptytop) pop();

	}

	inline void push(Imgidx idx, _int32 alpha)
	{
		//double tt = get_cpu_time(); //tmp

		//MinList1<Imgidx> *p, *q;
		_int16 i;
#if QUEUE_DEBUG
		cout << "chierQ - push: " << idx << " at " << (int)alpha << endl;
#endif
		//printf("Q: push %d at level %d\n", (int)idx, (int)alpha);

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		if(emptytop && alpha < list[0].alpha)
		{
			emptytop = 0;
			list[0].pidx = idx;
			list[0].alpha = alpha;
			return;
		}

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if ((_int64)alpha < hqueue->get_minlev())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
				curSize_list++;
			}
			else if (alpha < list[curSize_list].alpha)// push to the full list
			{
				push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

				for (i = curSize_list - 1; alpha < list[i].alpha; i--)
					list[i + 1] = list[i];
				list[i + 1].pidx = idx;
				list[i + 1].alpha = alpha;
			}
			else
				push_queue(idx, alpha); // push to the queue
		}
		else
			push_queue(idx, alpha); // push to the queue


			//qtime += get_cpu_time() - tt; //tmp
	}
	inline void push_queue(Imgidx idx, _int32 alpha)
	{
		hqueue->push(idx, alpha);
	}
	inline Imgidx pop()
	{
		Imgidx ret = top();
		//double tt = get_cpu_time(); //tmp

		_int8 i;
		//Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;
#if QUEUE_DEBUG
		cout << "chierQ - pop: " << list[0].pidx << " at " << (int)list[0].alpha << endl;
#endif

		//tmp

		//printf("Q: pop %d at level %d\n", (int)list[0].pidx, (int)list[0].alpha);

#if TRACK_QUEUEING
		f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 0)
		{
			list[0].pidx = hqueue->top();
			list[0].alpha = hqueue->get_minlev();

			pop_queue();
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		//f << list[0] << endl;
#endif

		//qtime += get_cpu_time() - tt; //tmp
		return ret;
	}
	inline void pop_queue()
	{
		hqueue->pop();
		hqueue->find_minlev();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	~Cache_Hierqueue_l2()
	{
		delete hqueue;
		Free(list - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};

template<class Imgidx, class Trieidx>//, class Qidx>
class Trie_Cache
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	Imgidx *list;
	Trie<Imgidx, Trieidx> *trie;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;


#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
	//	_int32 cnt;
	void initHQ(Imgidx size, size_t listsize)
	{
		/*		cnt = 0;//tmp*/
		//Imgidx i;
		this->maxSize_queue = size;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		//
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
		trie = new Trie<Imgidx, _int64>(size);
		list = (Imgidx*)Malloc((listsize + 1) * sizeof(Imgidx));
		list[0] = 0;
		list++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		minidx_queue = size;


		//tmp
#if TRACK_QUEUEING
		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
		f << -1 << '\n' << size << endl;
#endif
	}
public:
	Trie_Cache(Imgidx size)
	{
		initHQ(size, LISTSIZE_DEFAULT);
	}
	Trie_Cache(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}

	inline Imgidx get_minlev() { return list[0]; }
	inline Imgidx top() { return list[0]; }
	inline void push(Imgidx idx)
	{
		//MinList1<Imgidx> *p, *q;
		_int16 i;

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if (idx < trie->top())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; idx < list[i]; i--)
					list[i + 1] = list[i];
				list[i + 1] = idx;
				curSize_list++;
			}
			else if (idx < list[curSize_list])// push to the full list
			{
				push_queue(list[curSize_list]);

				for (i = curSize_list - 1; idx < list[i]; i--)
					list[i + 1] = list[i];
				list[i + 1] = idx;
			}
			else
				push_queue(idx); // push to the queue
		}
		else
			push_queue(idx); // push to the queue
	}
	inline void push_queue(Imgidx idx)
	{
		trie->push(idx);
	}
	inline void pop()
	{
		_int8 i;
		//Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;


		//tmp

#if TRACK_QUEUEING
		f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 0)
		{
			list[0] = trie->top();

			pop_queue();
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		f << list[0] << endl;
#endif
	}
	inline void pop_queue()
	{
		trie->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	// 	_int8 checklist()//tmp
	// 	{
	// 		MinList<Imgidx> *p;
	// 		if (head)
	// 		{
	// 			for (p = head; p; p = p->next)
	// 			{
	// 				if (p->idx < 0 || p->next && p->idx > p->next->idx)
	// 					return 1;
	// 			}
	// 			if (tail->next)
	// 				return 1;
	// 		}
	// 		int n = 0;
	// 		for (int i = 0; i < maxSize_list; i++)
	// 		{
	// 			if (list[i].idx != -1)
	// 				n++;
	// 		}
	// 		if (n != curSize_list)
	// 			return 1;
	// 		return 0;
	// 	}

	~Trie_Cache()
	{
		delete trie;
		Free(list - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};

template<class Imgidx>//, class Qidx>
class HybridQueue_HQueue_Rank
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	Imgidx *list;
	HQueue_l1idx_rank<Imgidx> *queue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;


#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
	//	_int32 cnt;
	void initHQ(Imgidx size, size_t listsize)
	{
		/*		cnt = 0;//tmp*/
		Imgidx i;
		this->maxSize_queue = size;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		//
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
		queue = new HQueue_l1idx_rank<Imgidx>(size);
		list = (Imgidx*)Malloc((listsize + 1) * sizeof(Imgidx));
		list[0] = 0;
		list++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		minidx_queue = size;


		//tmp
#if TRACK_QUEUEING
		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
		f << -1 << '\n' << size << endl;
#endif
	}
public:
	HybridQueue_HQueue_Rank(Imgidx size)
	{
		initHQ(size, LISTSIZE_DEFAULT);
	}
	HybridQueue_HQueue_Rank(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}

// 	double get_jumpnum() { return queue->jumpnum; }
// 	double get_jumpdist() { return queue->jumpdist; }

	inline Imgidx get_minlev() { return list[0]; }
	inline Imgidx top() { return list[0]; }
	inline void push(Imgidx idx)
	{
		//MinList1<Imgidx> *p, *q;
		_int16 i;

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		if (idx < queue->top())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; idx < list[i]; i--)
					list[i + 1] = list[i];
				list[i + 1] = idx;
				curSize_list++;
			}
			else if (idx < list[curSize_list])// push to the full list
			{
				push_queue(list[curSize_list]);

				for (i = curSize_list - 1; idx < list[i]; i--)
					list[i + 1] = list[i];
				list[i + 1] = idx;
			}
			else
				push_queue(idx); // push to the queue
		}
		else
			push_queue(idx); // push to the queue
	}
	inline void push_queue(Imgidx idx)
	{
		queue->push(idx);
	}
	inline void pop()
	{
		_int8 i;
		Imgidx idx;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;


		//tmp

#if TRACK_QUEUEING
		f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 0)
		{
			list[0] = queue->top();

			pop_queue();
		}
		else
		{
			for (i = 0; i < curSize_list; i++)
				list[i] = list[i + 1];
			curSize_list--;
		}

#if TRACK_QUEUEING
		f << list[0] << endl;
#endif
	}
	inline void pop_queue()
	{
		queue->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	// 	_int8 checklist()//tmp
	// 	{
	// 		MinList<Imgidx> *p;
	// 		if (head)
	// 		{
	// 			for (p = head; p; p = p->next)
	// 			{
	// 				if (p->idx < 0 || p->next && p->idx > p->next->idx)
	// 					return 1;
	// 			}
	// 			if (tail->next)
	// 				return 1;
	// 		}
	// 		int n = 0;
	// 		for (int i = 0; i < maxSize_list; i++)
	// 		{
	// 			if (list[i].idx != -1)
	// 				n++;
	// 		}
	// 		if (n != curSize_list)
	// 			return 1;
	// 		return 0;
	// 	}

	~HybridQueue_HQueue_Rank()
	{
		delete queue;
		Free(list - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};

template<class Imgidx>//, class Qidx>
class HybridQueue_HQueue_Rank1
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	Imgidx *list;
	HQueue_l1idx_rank<Imgidx> *queue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 shamt, nbit;

	_int16 l0, mask;


#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
	//	_int32 cnt;
	void initHQ(Imgidx size, size_t listsize)
	{
		/*		cnt = 0;//tmp*/
		//Imgidx i;
		this->maxSize_queue = size;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		//
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
		queue = new HQueue_l1idx_rank<Imgidx>(size);

		if (listsize > 256)
			listsize = 256;
		else if (listsize > 128)
			listsize = 256;
		else if (listsize > 64)
			listsize = 128;
		else if (listsize > 32)
			listsize = 64;
		else if (listsize > 16)
			listsize = 32;
		else if (listsize > 8)
			listsize = 16;
		else if (listsize > 4)
			listsize = 8;
		else
			listsize = 4;

		mask = listsize - 1;
		list = (Imgidx*)Malloc((listsize) * sizeof(Imgidx));
		maxSize_list = listsize;
		curSize_list = 0;
		l0 = listsize >> 1;
		list[l0] = size;

		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		minidx_queue = size;


		//tmp
#if TRACK_QUEUEING
		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
		f << -1 << '\n' << size << endl;
#endif
	}
public:
	HybridQueue_HQueue_Rank1(Imgidx size)
	{
		initHQ(size, LISTSIZE_DEFAULT);
	}
	HybridQueue_HQueue_Rank1(Imgidx size, size_t listsize)
	{
		initHQ(size, listsize);
	}

	// 	double get_jumpnum() { return queue->jumpnum; }
	// 	double get_jumpdist() { return queue->jumpdist; }

	inline Imgidx get_minlev() { return list[l0]; }
	inline Imgidx top() { return list[l0]; }
	inline void push(Imgidx idx)
	{
		//MinList1<Imgidx> *p, *q;
		_int16 i, j, lm;

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;
		lm = (l0 - 1) & mask;
		if (idx < queue->top())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				if (idx < list[l0])
				{
					list[lm] = idx;
					l0 = lm;
				}
				else
				{
					i = (l0 + curSize_list) & mask;
					j = (l0 + curSize_list - 1) & mask;
					while (idx < list[j])
					{
						list[i] = list[j];
						i = j;
						j = (j - 1) & mask;
					}
					list[i] = idx;
				}
				curSize_list++;
			}
			else if (idx < list[curSize_list])// push to the full list
			{
				push_queue(list[lm]);

				if (idx < list[l0])
				{
					list[lm] = idx;
					l0 = lm;
				}
				else
				{
					i = lm;
					j = (lm - 1) & mask;
					while (idx < list[j])
					{
						list[i] = list[j];
						i = j;
						j = (j - 1) & mask;
					}
					list[i] = idx;
				}
			}
			else
				push_queue(idx); // push to the queue
		}
		else
			push_queue(idx); // push to the queue
	}
	inline void push_queue(Imgidx idx)
	{
		queue->push(idx);
	}
	inline void pop()
	{
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;


		//tmp

#if TRACK_QUEUEING
		f << '1' << '\n' << list[0] << endl;
#endif

		if (curSize_list == 1)
		{
			list[l0] = queue->top();

			pop_queue();
		}
		else
		{
			l0 = (l0 + 1) & mask;
			curSize_list--;
		}

#if TRACK_QUEUEING
		f << list[0] << endl;
#endif
	}
	inline void pop_queue()
	{
		queue->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	// 	_int8 checklist()//tmp
	// 	{
	// 		MinList<Imgidx> *p;
	// 		if (head)
	// 		{
	// 			for (p = head; p; p = p->next)
	// 			{
	// 				if (p->idx < 0 || p->next && p->idx > p->next->idx)
	// 					return 1;
	// 			}
	// 			if (tail->next)
	// 				return 1;
	// 		}
	// 		int n = 0;
	// 		for (int i = 0; i < maxSize_list; i++)
	// 		{
	// 			if (list[i].idx != -1)
	// 				n++;
	// 		}
	// 		if (n != curSize_list)
	// 			return 1;
	// 		return 0;
	// 	}

	~HybridQueue_HQueue_Rank1()
	{
		delete queue;
		Free(list);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};


// template<class Imgidx>//, class Qidx>
// class HybridQueue_HH //Heap + Hierarchical
// {
// 	//MinList1<Imgidx> *list, *list_end, *head, *tail;
// 	HeapQueue_rank *heapqueue;				//prior queue
// 	HQueue_l1idx_rank<Imgidx> *hierarqueue; //secondary queue
// 	Imgidx minidx_hierarqueue, maxidx_heapqueue;
// //	_int16 curSize_list, maxSize_list; heapqueue->cursize
// 	Imgidx maxSize_queue, mask_field;
// 	_int8 shamt, nbit;
//
//
// #if TRACK_QUEUEING
// 	Imgidx *in_size;
//
// 	ofstream f;
// #endif
// 	//	_int32 cnt;
// 	void initHQ(Imgidx size, size_t heapsize)
// 	{
// 		/*		cnt = 0;//tmp*/
// 		Imgidx i;
// 		this->maxSize_queue = size;
// 		/*		shamt = 2;*/
// 		// 		nbit = sizeof(Qidx) * 8;
// 		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
// 		// 			shamt++;
// 		// 		mask_field = (1 << shamt) - 1;
// 		// 		qsize = (size + mask_field) >> shamt;
// 		//
// 		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
// 		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
// 		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
// 		hierarqueue = new HQueue_l1idx_rank<Imgidx>(size);
// 		heapqueue = new HeapQueue_rank<Imgidx>(heapsize);
// 		maxidx_heapqueue = 0;
// // 		list = (Imgidx*)Malloc((heapsize) * sizeof(Imgidx));
// //		list[0] = 0;
// //		list++;
// //		maxSize_list = heapsize - 1;
// //		curSize_list = -1;
// 		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
// 		//list_end = list + listsize;
// 		//maxSize_list = listsize;
// 		//head = tail = 0;
//
//
// 		//		for (i = 0; i < size; i++)
// 		//			queue[i] = -1;
// 		//		queue[size] = 0;
// 		//for (i = 0; i < listsize; i++)
// 			//list[i].idx = -1;
// 		//curSize_list = 0;
// 		//		minidx_queue = size >> shamt;
// 		minidx_hierarqueue = size;
//
//
// 		//tmp
// #if TRACK_QUEUEING
// 		f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
// 		f << -1 << '\n' << size << endl;
// #endif
// 	}
// public:
// 	HybridQueue_HQueue_Rank(Imgidx size)
// 	{
// 		initHQ(size, HEAPSIZE_DEFAULT);
// 	}
// 	HybridQueue_HQueue_Rank(Imgidx size, size_t heapsize)
// 	{
// 		initHQ(size, heapsize);
// 	}
//
// 	// 	double get_jumpnum() { return queue->jumpnum; }
// 	// 	double get_jumpdist() { return queue->jumpdist; }
//
// 	inline Imgidx get_minlev() { return heapqueue->get_minlev(); }
// 	inline Imgidx top() { return heapqueue->top(); }
// 	inline void push(Imgidx idx)
// 	{
// 		//MinList1<Imgidx> *p, *q;
//
// #if TRACK_QUEUEING
// 		//tmp
// 		f << '0' << '\n' << idx << endl;
// #endif
//
// 		// 		cnt++;//tmp
// 		//
// 		// 		if (cnt == 786)//tmp
// 		// 			idx = idx;
// 		if (idx < hierarqueue->top())
// 		{
// 			if (heapqueue->cursize < heapqueue->maxsize) //spare room in the list
// 			{
// 				maxidx_heapqueue = maxidx_heapqueue > idx ? maxidx_heapqueue: idx;
// 				heapqueue->push(idx);
// 			}
// 			else if (idx < )// push to the full list
// 			{
// 				for (Imgidx i = heapqueue->cursize; heapqueue->get_arr(i) != maxidx_heapqueue; i--)
// 					;
//
// 			}
// 			else
// 				push_queue(idx); // push to the queue
// 		}
// 		else
// 			push_queue(idx); // push to the queue
// 	}
// 	inline void push_queue(Imgidx idx)
// 	{
// 		hierarqueue->push(idx);
// 	}
// 	inline void pop()
// 	{
// 		_int8 i;
// 		Imgidx idx;
// 		// 		cnt++;//tmp
// 		// 		if (cnt == 776)//tmp
// 		// 			cnt = cnt;
//
//
// 		//tmp
//
// #if TRACK_QUEUEING
// 		f << '1' << '\n' << list[0] << endl;
// #endif
//
// 		if (curSize_list == 0)
// 		{
// 			list[0] = hierarqueue->top();
//
// 			pop_queue();
// 		}
// 		else
// 		{
// 			for (i = 0; i < curSize_list; i++)
// 				list[i] = list[i + 1];
// 			curSize_list--;
// 		}
//
// #if TRACK_QUEUEING
// 		f << list[0] << endl;
// #endif
// 	}
// 	inline void pop_queue()
// 	{
// 		hierarqueue->pop();
// 		// 		queue[minidx_queue] = -1;
// 		// 		while (queue[++minidx_queue] == -1)
// 		// 			;
// 	}
//
// 	// 	_int8 checklist()//tmp
// 	// 	{
// 	// 		MinList<Imgidx> *p;
// 	// 		if (head)
// 	// 		{
// 	// 			for (p = head; p; p = p->next)
// 	// 			{
// 	// 				if (p->idx < 0 || p->next && p->idx > p->next->idx)
// 	// 					return 1;
// 	// 			}
// 	// 			if (tail->next)
// 	// 				return 1;
// 	// 		}
// 	// 		int n = 0;
// 	// 		for (int i = 0; i < maxSize_list; i++)
// 	// 		{
// 	// 			if (list[i].idx != -1)
// 	// 				n++;
// 	// 		}
// 	// 		if (n != curSize_list)
// 	// 			return 1;
// 	// 		return 0;
// 	// 	}
//
// 	~HybridQueue_HQueue_Rank()
// 	{
// 		delete hierarqueue;
// 		Free(list - 1);
// 		//Free(queue);
//
// #if TRACK_QUEUEING
// 		//tmp
// 		f.close();
// #endif
// 	}
// };



template<class Imgidx>//, class Qidx>
class HybridQueue_HQueue
{
	//MinList1<Imgidx> *list, *list_end, *head, *tail;
	Imgidx *list;
	_int64 *levels;
	HQueue_l1idx<Imgidx> *queue;
	Imgidx minidx_queue;
	_int16 curSize_list, maxSize_list;
	Imgidx maxSize_queue, mask_field;
	_int8 minlevnotfixed;

#if TRACK_QUEUEING
	Imgidx *in_size;

	ofstream f;
#endif
	//	_int32 cnt;
	void initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
	{
		/*		cnt = 0;//tmp*/
		Imgidx i;
		this->maxSize_queue = qsize_in;
		/*		shamt = 2;*/
		// 		nbit = sizeof(Qidx) * 8;
		// 		for (_int8 nbyte = sizeof(Qidx); nbyte; nbyte >>= 1)
		// 			shamt++;
		// 		mask_field = (1 << shamt) - 1;
		// 		qsize = (size + mask_field) >> shamt;
		//
		// 		queue = (Qidx*)Malloc(qsize * sizeof(Qidx*));
		//queue = (_int8*)Malloc((size + 1) * sizeof(_int8));
		//trie = (Trie<Imgidx, _int64>*)Malloc(size * sizeof(Trie<Imgidx, _int64>*));
		queue = new HQueue_l1idx<Imgidx>(qsize_in, dhist, numlevels);
		list = (Imgidx*)Malloc((listsize) * sizeof(Imgidx));
		levels = (_int64*)Malloc((listsize + 1) * sizeof(_int64));
		levels[0] = 0;
		levels++;
		maxSize_list = listsize - 1;
		curSize_list = -1;
		//list = (MinList1<Imgidx>*)Malloc(listsize * sizeof(MinList1<Imgidx>));
		//list_end = list + listsize;
		//maxSize_list = listsize;
		//head = tail = 0;


		//		for (i = 0; i < size; i++)
		//			queue[i] = -1;
		//		queue[size] = 0;
		//for (i = 0; i < listsize; i++)
			//list[i].idx = -1;
		//curSize_list = 0;
		//		minidx_queue = size >> shamt;
		minidx_queue = qsize_in;

		minlevnotfixed = 0;
		//tmp
#if TRACK_QUEUEING
		//f.open("D:/RUG/2019/TTMA_ISMM/queuelog.dat", std::ofstream::app);
		//f << -1 << '\n' << size << endl;
#endif
	}
public:
	HybridQueue_HQueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
	{
		initHQ(qsize_in, dhist, numlevels, LISTSIZE_DEFAULT);
	}
	HybridQueue_HQueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
	{
		initHQ(qsize_in, dhist, numlevels, listsize);
	}
	inline Imgidx top(){return list[0];}
	inline _int64 get_minlev() { return levels[0]; }
	inline void find_minlev()
	{
//		if (minlevnotfixed)
			queue->find_minlev();

	}
	inline void push(Imgidx idx, _int64 level)
	{
		//MinList1<Imgidx> *p, *q;
		_int16 i;

//		cout << "- pushing " << (int)idx << " at level " << (int)level << endl;

#if TRACK_QUEUEING
		//tmp
		f << '0' << '\n' << idx << endl;
#endif

		// 		cnt++;//tmp
		//
		// 		if (cnt == 786)//tmp
		// 			idx = idx;

		if (curSize_list == -1) //should be run only the first time
		{
			list[0] = idx;
			levels[0] = level;
			curSize_list = 0;
			push_queue(idx, level);
			return;
		}

		if (level < queue->get_minlev())
		{
			if (curSize_list < maxSize_list) //spare room in the list
			{
				for (i = curSize_list; level < levels[i]; i--)
				{
					list[i + 1] = list[i];
					levels[i + 1] = levels[i];
				}
				list[i + 1] = idx;
				levels[i + 1] = level;
				curSize_list++;
			}
			else if (level < levels[curSize_list])// push to the full list
			{
				push_queue(list[curSize_list], levels[curSize_list]);

				for (i = curSize_list - 1; level < levels[i]; i--)
				{
					list[i + 1] = list[i];
					levels[i + 1] = levels[i];
				}
				list[i + 1] = idx;
				levels[i + 1] = level;
			}
			else
				push_queue(idx, level); // push to the queue
		}
		else
			push_queue(idx, level); // push to the queue
	}
	inline void push_queue(Imgidx idx, _int64 level)
	{
		if (queue->push(idx, level))
			minlevnotfixed = 0;
	}
	inline Imgidx pop()
	{
		Imgidx ret;
		// 		cnt++;//tmp
		// 		if (cnt == 776)//tmp
		// 			cnt = cnt;

//		cout << "- popping " << (int)list[0] << " at level " << (int)levels[0] << endl;

		//tmp

#if TRACK_QUEUEING
		f << '1' << '\n' << list[0] << endl;
#endif

		ret = list[0];

		if (curSize_list == 0)
		{
			list[0] = queue->top();
			levels[0] = queue->get_minlev();

			pop_queue();
		}
		else
		{
			for (_int8 i = 0; i < curSize_list; i++)
			{
				list[i] = list[i + 1];
				levels[i] = levels[i + 1];
			}
			curSize_list--;
		}
		return ret;

#if TRACK_QUEUEING
		f << list[0] << endl;
#endif
	}
	inline void pop_queue()
	{
		minlevnotfixed = 1;
		queue->pop();
		// 		queue[minidx_queue] = -1;
		// 		while (queue[++minidx_queue] == -1)
		// 			;
	}

	// 	_int8 checklist()//tmp
	// 	{
	// 		MinList<Imgidx> *p;
	// 		if (head)
	// 		{
	// 			for (p = head; p; p = p->next)
	// 			{
	// 				if (p->idx < 0 || p->next && p->idx > p->next->idx)
	// 					return 1;
	// 			}
	// 			if (tail->next)
	// 				return 1;
	// 		}
	// 		int n = 0;
	// 		for (int i = 0; i < maxSize_list; i++)
	// 		{
	// 			if (list[i].idx != -1)
	// 				n++;
	// 		}
	// 		if (n != curSize_list)
	// 			return 1;
	// 		return 0;
	// 	}

	~HybridQueue_HQueue()
	{
		delete queue;
		Free(list);
		Free(levels - 1);
		//Free(queue);

#if TRACK_QUEUEING
		//tmp
		f.close();
#endif
	}
};

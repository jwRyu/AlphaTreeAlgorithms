#include "HybridQueue.h"

template<class Pixel>
void HierarHeapQueue_HEQ<Pixel>::initHQ(Imgidx *dhist, _uint32* histeqmap_in, Imgidx numlevels_in, Imgidx size, double a_in, int listsize)
{
    list = (HQentry<Pixel>*)Malloc(listsize * sizeof(HQentry<Pixel>));
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

    hqueue = (HeapQueue_naive_quad<Pixel>**)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Pixel>*));
    for(int level = 0;level < thr_hqueue;level++)
        hqueue[level] = new HeapQueue_naive_quad<Pixel>(qsizes[level]);

    storage = (HQentry<Pixel>**)Calloc((numlevels - thr_hqueue) * sizeof(HQentry<Pixel>*));
    storage -= thr_hqueue;
    for(int level = thr_hqueue;level < numlevels;level++)
        storage[level] = (HQentry<Pixel>*)Malloc(qsizes[level] * sizeof(HQentry<Pixel>));
}

template<class Pixel>
HierarHeapQueue_HEQ<Pixel>::HierarHeapQueue_HEQ(Imgidx *dhist, _uint32 *histeqmap_in, Imgidx numlevels_in, double a_in, Imgidx size)
{
    initHQ(dhist, histeqmap_in, numlevels_in, size, a_in, 12);
}

template<class Pixel>
HierarHeapQueue_HEQ<Pixel>::HierarHeapQueue_HEQ(Imgidx *dhist, _uint32 *histeqmap_in, Imgidx numlevels_in, Imgidx size, double a_in, int listsize)
{
    initHQ(dhist, histeqmap_in, numlevels_in, size, a_in, listsize);
}

template<class Pixel>
HierarHeapQueue_HEQ<Pixel>::~HierarHeapQueue_HEQ()
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

template<class Pixel>
void HierarHeapQueue_HEQ<Pixel>::push_1stitem(Imgidx idx, Pixel alpha)
{
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template<class Pixel>
void HierarHeapQueue_HEQ<Pixel>::end_pushes(_uint8 *isVisited)
{
    if(emptytop)
        pop(isVisited);
}

template<class Pixel>
void HierarHeapQueue_HEQ<Pixel>::push(Imgidx idx, Pixel alpha)
{
    if(emptytop && alpha < list[0].alpha)
    {
        emptytop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }

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
}

template<class Pixel>
void HierarHeapQueue_HEQ<Pixel>::push_queue(Imgidx idx, Pixel alpha)
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

template<class Pixel>
Imgidx HierarHeapQueue_HEQ<Pixel>::pop(_uint8 *isVisited)
{
    Imgidx ret = top();
    if (curSize_list == 0)
    {
        while(!check_queue_level(isVisited))
            queue_minlev++;
        list[0].pidx = hqueue[queue_minlev]->top();
        list[0].alpha = hqueue[queue_minlev]->top_alpha();

        pop_queue(isVisited);
    }
    else
    {
        for (int i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }
    return ret;
}

template<class Pixel>
int HierarHeapQueue_HEQ<Pixel>::check_queue_level(_uint8 *isVisited)
{
    if(queue_minlev < curthr)
        return hqueue[queue_minlev]->get_cursize();
    else
    {
        //while(curthr < queue_minlev)
        //{
        //	printf("Piep ");
        //	hqueue[curthr] = new HeapQueue_naive_quad<Pixel>(qsizes[curthr]);

        //	Free(storage[curthr]);
        //	storage[curthr] = 0;
        //	curthr++;
        //}
        curthr++;

        hqueue[queue_minlev] = new HeapQueue_naive_quad<Pixel>(qsizes[queue_minlev]);

        HQentry<Pixel>* store = storage[queue_minlev];
        Imgidx cur = storage_cursize[queue_minlev];
        HeapQueue_naive_quad<Pixel> *pQ = hqueue[queue_minlev];
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

template<class Pixel>
void HierarHeapQueue_HEQ<Pixel>::pop_queue(_uint8 *isVisited)
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

template class HierarHeapQueue_HEQ<_uint8>;
template class HierarHeapQueue_HEQ<_uint16>;
template class HierarHeapQueue_HEQ<_uint32>;
template class HierarHeapQueue_HEQ<_uint64>;


template<class Pixel>
HierarHeapQueue<Pixel>::HierarHeapQueue(Imgidx *dhist, Imgidx numlevels_in, Imgidx size, double a_in, int listsize, int connectivity, double r)
{
    this->numlevels = numlevels_in;
    this->a = a_in;
    this->queue_minlev = numlevels;

    Imgidx cumsum = 0;
    qsizes = dhist; //do not free dhist outside
    if(r >= 1)
    {
        thr_hqueue = curthr = numlevels;
        hqueue = (HeapQueue_naive_quad<Pixel>**)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Pixel>*));
        for(int level = 0;level < thr_hqueue;level++)
            hqueue[level] = new HeapQueue_naive_quad<Pixel>(qsizes[level]);
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

        hqueue = (HeapQueue_naive_quad<Pixel>**)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Pixel>*));
        for(int level = 0;level < thr_hqueue;level++)
            hqueue[level] = new HeapQueue_naive_quad<Pixel>(qsizes[level]);
        //hqueue[numlevels] = new HeapQueue_naive_quad<Pixel>(1);
        //hqueue[numlevels]->push(0,(Pixel)-1);

        storage = (HQentry<Pixel>**)Calloc((numlevels - thr_hqueue) * sizeof(HQentry<Pixel>*));
        storage -= thr_hqueue;
        for(int level = thr_hqueue;level < numlevels;level++)
            storage[level] = (HQentry<Pixel>*)Malloc(qsizes[level] * sizeof(HQentry<Pixel>));
    }
}

template<class Pixel>
HierarHeapQueue<Pixel>::~HierarHeapQueue()
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

template<class Pixel>
void HierarHeapQueue<Pixel>::push_queue(Imgidx idx, Pixel alpha)
{
    int level = (int)(a * log2(1 + (double)alpha));
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

template<class Pixel>
Imgidx HierarHeapQueue<Pixel>::pop(_uint8 *isVisited)
{
    Imgidx ret = top();

    pop_queue(isVisited);
    return ret;
}

template<class Pixel>
int HierarHeapQueue<Pixel>::check_queue_level(_uint8* isVisited)
{
    if(queue_minlev < curthr)
        return hqueue[queue_minlev]->get_cursize();
    else
    {
        while(curthr < queue_minlev)
        {
            hqueue[curthr] = new HeapQueue_naive_quad<Pixel>(qsizes[curthr]);

            Free(storage[curthr]);
            storage[curthr] = 0;
            curthr++;
        }
        curthr++;

        hqueue[queue_minlev] = new HeapQueue_naive_quad<Pixel>(qsizes[queue_minlev]);

        HQentry<Pixel>* store = storage[queue_minlev];
        Imgidx cur = storage_cursize[queue_minlev];
        HeapQueue_naive_quad<Pixel> *pQ = hqueue[queue_minlev];
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

template<class Pixel>
void HierarHeapQueue<Pixel>::pop_queue(_uint8* isVisited)
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

template class HierarHeapQueue<_uint8>;
template class HierarHeapQueue<_uint16>;
template class HierarHeapQueue<_uint32>;
template class HierarHeapQueue<_uint64>;

//////////////////////////////////////////////////////////////////////////////////////////////////////
// HierarHeapQueue_cache

template<class Pixel>
void HierarHeapQueue_cache<Pixel>::initHQ(Imgidx *dhist, Imgidx numlevels_in, Imgidx size, double a_in, int listsize, int connectivity, double r)
{
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

    list = (HQentry<Pixel>*)Malloc((listsize) * sizeof(HQentry<Pixel>));
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
        hqueue = (HeapQueue_naive_quad<Pixel>**)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Pixel>*));
        for(int level = 0;level < thr_hqueue;level++)
            hqueue[level] = new HeapQueue_naive_quad<Pixel>(qsizes[level]);
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

        hqueue = (HeapQueue_naive_quad<Pixel>**)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Pixel>*));
        for(int level = 0;level < thr_hqueue;level++)
            hqueue[level] = new HeapQueue_naive_quad<Pixel>(qsizes[level]);

        storage = (HQentry<Pixel>**)Calloc((numlevels - thr_hqueue) * sizeof(HQentry<Pixel>*));
        storage -= thr_hqueue;
        for(int level = thr_hqueue;level < numlevels;level++)
            storage[level] = (HQentry<Pixel>*)Malloc(qsizes[level] * sizeof(HQentry<Pixel>));
    }		
}

template<class Pixel>
HierarHeapQueue_cache<Pixel>::HierarHeapQueue_cache(Imgidx *dhist, Imgidx numlevels_in, double a_in, Imgidx size, Imgidx connectivity, double r)
{
    initHQ(dhist, numlevels_in, size, a_in, 12, (int)connectivity, r);
}

template<class Pixel>
HierarHeapQueue_cache<Pixel>::HierarHeapQueue_cache(Imgidx *dhist, Imgidx numlevels_in, Imgidx size, double a_in, int listsize, Imgidx connectivity, double r)
{
    initHQ(dhist, numlevels_in, size, a_in, listsize, (int)connectivity, r);
}

template<class Pixel>
HierarHeapQueue_cache<Pixel>::~HierarHeapQueue_cache()
{
    Free(list);
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

template<class Pixel>
void HierarHeapQueue_cache<Pixel>::push_1stitem(Imgidx idx, Pixel alpha)
{
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template<class Pixel>
void HierarHeapQueue_cache<Pixel>::end_pushes(_uint8 *isVisited)
{
    if(emptytop)
        pop(isVisited);
}

template<class Pixel>
void HierarHeapQueue_cache<Pixel>::push(Imgidx idx, Pixel alpha)
{
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
}

template<class Pixel>
void HierarHeapQueue_cache<Pixel>::push_queue(Imgidx idx, Pixel alpha)
{
    int level = (int)(a * log2(1 + (double)alpha));

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

template<class Pixel>
Imgidx HierarHeapQueue_cache<Pixel>::pop(_uint8 *isVisited)
{
    Imgidx ret = top();
    if (curSize_list == 0)
    {
        while(!check_queue_level(isVisited))
            queue_minlev++;
        list[0].pidx = hqueue[queue_minlev]->top();
        list[0].alpha = hqueue[queue_minlev]->top_alpha();

        pop_queue(isVisited);
    }
    else
    {
        for (int i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }

    return ret;
}

template<class Pixel>
int HierarHeapQueue_cache<Pixel>::check_queue_level(_uint8 *isVisited)
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
            hqueue[curthr] = new HeapQueue_naive_quad<Pixel>(qsizes[curthr]);

            Free(storage[curthr]);
            storage[curthr] = 0;
            curthr++;
        }
        curthr++;

        hqueue[queue_minlev] = new HeapQueue_naive_quad<Pixel>(qsizes[queue_minlev]);

        HQentry<Pixel>* store = storage[queue_minlev];
        Imgidx cur = storage_cursize[queue_minlev];
        HeapQueue_naive_quad<Pixel> *pQ = hqueue[queue_minlev];
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

template<class Pixel>
void HierarHeapQueue_cache<Pixel>::pop_queue(_uint8 *isVisited)
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

template class HierarHeapQueue_cache<_uint8>;
template class HierarHeapQueue_cache<_uint16>;
template class HierarHeapQueue_cache<_uint32>;
template class HierarHeapQueue_cache<_uint64>;

// HierarHeapQueue_cache end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Cache_Heapqueue

template<class Pixel>
void Cache_Heapqueue<Pixel>::initHQ(Imgidx size, size_t listsize)
{
    this->maxSize_queue = size;
    hqueue = new HeapQueue_naive<Pixel>(size);
    list = (HQentry<Pixel>*)Malloc((listsize + 1) * sizeof(HQentry<Pixel>));
    list[0].pidx = 0;
    list[0].alpha = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    qtime = 0;//tmp
}

template<class Pixel>
Cache_Heapqueue<Pixel>::Cache_Heapqueue(Imgidx size)
{
    initHQ(size, 12);
}

template<class Pixel>
Cache_Heapqueue<Pixel>::Cache_Heapqueue(Imgidx size, size_t listsize)
{
    initHQ(size, listsize);
}

template<class Pixel>
Cache_Heapqueue<Pixel>::~Cache_Heapqueue()
{
    delete hqueue;
    Free(list - 1);
}

template<class Pixel>
void Cache_Heapqueue<Pixel>::push_1stitem(Imgidx idx, Pixel alpha)
{
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template<class Pixel>
void Cache_Heapqueue<Pixel>::push(Imgidx idx, Pixel alpha)
{
    _int16 i;
    if(emptytop && alpha < list[0].alpha)
    {
        emptytop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }

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
        else if (alpha < list[curSize_list].alpha) // push to the full list
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
}

template<class Pixel>
Imgidx Cache_Heapqueue<Pixel>::pop()
{
    Imgidx ret = top();
    _int8 i;
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

    return ret;
}

template class Cache_Heapqueue<_uint8>;
template class Cache_Heapqueue<_uint16>;
template class Cache_Heapqueue<_uint32>;
template class Cache_Heapqueue<_uint64>;
template class Cache_Heapqueue<float>;
template class Cache_Heapqueue<double>;

// Cache_Heapqueue end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Cache_Quad_Heapqueue

template<class Pixel>
void Cache_Quad_Heapqueue<Pixel>::initHQ(Imgidx size, size_t listsize)
{
    this->maxSize_queue = size;
    hqueue = new HeapQueue_naive_quad<Pixel>(size);
    list = (HQentry<Pixel>*)Malloc((listsize + 1) * sizeof(HQentry<Pixel>));
    list[0].pidx = 0;
    list[0].alpha = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    qtime = 0;//tmp
}

template<class Pixel>
Cache_Quad_Heapqueue<Pixel>::Cache_Quad_Heapqueue(Imgidx size)
{
    initHQ(size, 12);
}

template<class Pixel>
Cache_Quad_Heapqueue<Pixel>::Cache_Quad_Heapqueue(Imgidx size, size_t listsize)
{
    initHQ(size, listsize);
}

template<class Pixel>
Cache_Quad_Heapqueue<Pixel>::~Cache_Quad_Heapqueue()
{
    delete hqueue;
    Free(list - 1);
}

template<class Pixel>
void Cache_Quad_Heapqueue<Pixel>::push_1stitem(Imgidx idx, Pixel alpha)
{
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template<class Pixel>
void Cache_Quad_Heapqueue<Pixel>::push(Imgidx idx, Pixel alpha)
{
    _int16 i;

    if(emptytop && alpha < list[0].alpha)
    {
        emptytop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }
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
}

template<class Pixel>
Imgidx Cache_Quad_Heapqueue<Pixel>::pop()
{
    Imgidx ret = top();
    _int8 i;

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
    return ret;
}

template class Cache_Quad_Heapqueue<_uint8>;
template class Cache_Quad_Heapqueue<_uint16>;
template class Cache_Quad_Heapqueue<_uint32>;
template class Cache_Quad_Heapqueue<_uint64>;
template class Cache_Quad_Heapqueue<float>;
template class Cache_Quad_Heapqueue<double>;

// Cache_Quad_Heapqueue end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// CirCache_Hierqueue

template<class Pixel>
void CirCache_Hierqueue<Pixel>::initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
{
    int shamt = 0;
    for(int lsize = listsize;lsize;lsize>>=1)
        shamt++;
    listsize = 1 << (shamt - 1);
    mask = listsize - 1;
    this->maxSize_queue = qsize_in;
    hqueue = new HierarQueue(qsize_in, dhist, numlevels);
    list = (HQentry<_int32>*)Malloc(listsize * sizeof(HQentry<_int32>));
    maxSize_list = listsize;
    curSize_list = 0;
    liststart = 0;
    qtime = 0;//tmp
}

template<class Pixel>
CirCache_Hierqueue<Pixel>::CirCache_Hierqueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
{
    initHQ(qsize_in, dhist, numlevels, 16);
}

template<class Pixel>
CirCache_Hierqueue<Pixel>::CirCache_Hierqueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
{
    initHQ(qsize_in, dhist, numlevels, listsize);
}

template<class Pixel>
CirCache_Hierqueue<Pixel>::~CirCache_Hierqueue()
{
    delete hqueue;
    Free(list);
}

template<class Pixel>
void CirCache_Hierqueue<Pixel>::push_1stitem(Imgidx idx, _int32 alpha)
{
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template<class Pixel>
void CirCache_Hierqueue<Pixel>::push(Imgidx idx, _int32 alpha)
{
    _int16 i, j, k;

    if(emptytop && alpha < list[liststart].alpha)
    {
        emptytop = 0;
        list[liststart].pidx = idx;
        list[liststart].alpha = alpha;
        return;
    }
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
}

template<class Pixel>
Imgidx CirCache_Hierqueue<Pixel>::pop()
{
    Imgidx ret = top();

    if (curSize_list == 1)
    {
        list[liststart].pidx = hqueue->top();
        list[liststart].alpha = hqueue->get_minlev();

        pop_queue();
    }
    else
    {
        liststart = (liststart + 1) & mask;
        curSize_list--;
    }

    return ret;
}

template class CirCache_Hierqueue<_uint8>;
template class CirCache_Hierqueue<_uint16>;
template class CirCache_Hierqueue<_uint32>;
template class CirCache_Hierqueue<_uint64>;
template class CirCache_Hierqueue<float>;
template class CirCache_Hierqueue<double>;

// CirCache_Hierqueue end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HierarQueueCache

template<class Pixel>
void HierarQueueCache<Pixel>::initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
{
    this->maxSize_queue = qsize_in;
    hqueue = new HierarQueue(qsize_in, dhist, numlevels);
    list = (HQentry<_int32>*)Malloc((listsize + 1) * sizeof(HQentry<_int32>));
    list[0].pidx = 0;
    list[0].alpha = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    qtime = 0;//tmp
}

template<class Pixel>
HierarQueueCache<Pixel>::HierarQueueCache(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
{
    initHQ(qsize_in, dhist, numlevels, 12);
}

template<class Pixel>
HierarQueueCache<Pixel>::HierarQueueCache(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
{
    initHQ(qsize_in, dhist, numlevels, listsize);
}

template<class Pixel>
HierarQueueCache<Pixel>::~HierarQueueCache()
{
    delete hqueue;
    Free(list - 1);
}

template<class Pixel>
void HierarQueueCache<Pixel>::push_1stitem(Imgidx idx, _int32 alpha)
{
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template<class Pixel>
void HierarQueueCache<Pixel>::push(Imgidx idx, _int32 alpha)
{
    _int16 i;

    if(emptytop && alpha < list[0].alpha)
    {
        emptytop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }
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
        else if (alpha < list[curSize_list].alpha) // push to the full list
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
}

template<class Pixel>
Imgidx HierarQueueCache<Pixel>::pop()
{
    Imgidx ret = top();
    _int8 i;

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
    return ret;
}

template class HierarQueueCache<_uint8>;
template class HierarQueueCache<_uint16>;
template class HierarQueueCache<_uint32>;
template class HierarQueueCache<_uint64>;

// HierarQueueCache end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Cache_Hierqueue_l1

template<class Pixel>
void Cache_Hierqueue_l1<Pixel>::initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
{
    this->maxSize_queue = qsize_in;
    hqueue = new HQueue_l1idx(qsize_in, dhist, numlevels);
    list = (HQentry<_int32>*)Malloc((listsize + 1) * sizeof(HQentry<_int32>));
    list[0].pidx = 0;
    list[0].alpha = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;

    qtime = 0;//tmp
}

template<class Pixel>
Cache_Hierqueue_l1<Pixel>::Cache_Hierqueue_l1(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
{
    initHQ(qsize_in, dhist, numlevels, 12);
}

template<class Pixel>
Cache_Hierqueue_l1<Pixel>::Cache_Hierqueue_l1(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
{
    initHQ(qsize_in, dhist, numlevels, listsize);
}

template<class Pixel>
Cache_Hierqueue_l1<Pixel>::~Cache_Hierqueue_l1()
{
    delete hqueue;
    Free(list - 1);
}

template<class Pixel>
void Cache_Hierqueue_l1<Pixel>::push_1stitem(Imgidx idx, _int32 alpha)
{
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template<class Pixel>
void Cache_Hierqueue_l1<Pixel>::push(Imgidx idx, _int32 alpha)
{
    _int16 i;

    if(emptytop && alpha < list[0].alpha)
    {
        emptytop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }

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
}

template<class Pixel>
Imgidx Cache_Hierqueue_l1<Pixel>::pop()
{
    Imgidx ret = top();
    _int8 i;

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

    return ret;
}

template class Cache_Hierqueue_l1<_uint8>;
template class Cache_Hierqueue_l1<_uint16>;
template class Cache_Hierqueue_l1<_uint32>;
template class Cache_Hierqueue_l1<_uint64>;

// Cache_Hierqueue_l1 end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Cache_Hierqueue_l2

template<class Pixel>
void Cache_Hierqueue_l2<Pixel>::initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
{
    this->maxSize_queue = qsize_in;
    hqueue = new HQueue_l2idx(qsize_in, dhist, numlevels);
    list = (HQentry<_int32>*)Malloc((listsize + 1) * sizeof(HQentry<_int32>));
    list[0].pidx = 0;
    list[0].alpha = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    qtime = 0;//tmp
}

template<class Pixel>
Cache_Hierqueue_l2<Pixel>::Cache_Hierqueue_l2(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
{
    initHQ(qsize_in, dhist, numlevels, 12);
}

template<class Pixel>
Cache_Hierqueue_l2<Pixel>::Cache_Hierqueue_l2(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
{
    initHQ(qsize_in, dhist, numlevels, listsize);
}

template<class Pixel>
Cache_Hierqueue_l2<Pixel>::~Cache_Hierqueue_l2()
{
    delete hqueue;
    Free(list - 1);
}

template<class Pixel>
void Cache_Hierqueue_l2<Pixel>::push_1stitem(Imgidx idx, _int32 alpha)
{
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template<class Pixel>
void Cache_Hierqueue_l2<Pixel>::push(Imgidx idx, _int32 alpha)
{
    _int16 i;

    if(emptytop && alpha < list[0].alpha)
    {
        emptytop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }

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
}

template<class Pixel>
Imgidx Cache_Hierqueue_l2<Pixel>::pop()
{
    Imgidx ret = top();
    _int8 i;

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
    return ret;
}

template class Cache_Hierqueue_l2<_uint8>;
template class Cache_Hierqueue_l2<_uint16>;
template class Cache_Hierqueue_l2<_uint32>;
template class Cache_Hierqueue_l2<_uint64>;

// Cache_Hierqueue_l2 end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Trie_Cache

void Trie_Cache::initHQ(Imgidx size, size_t listsize)
{
    this->maxSize_queue = size;
    trie = new Trie<_uint64>(size);
    list = (Imgidx*)Malloc((listsize + 1) * sizeof(Imgidx));
    list[0] = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    minidx_queue = size;
}

Trie_Cache::Trie_Cache(Imgidx size)
{
    initHQ(size, LISTSIZE_DEFAULT);
}

Trie_Cache::Trie_Cache(Imgidx size, size_t listsize)
{
    initHQ(size, listsize);
}

Trie_Cache::~Trie_Cache()
{
    delete trie;
    Free(list - 1);
}

void Trie_Cache::push(Imgidx idx)
{
    _int16 i;
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

void Trie_Cache::pop()
{
    _int8 i;
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
}

// Trie_Cache end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HybridQueue_HQueue_Rank

void HybridQueue_HQueue_Rank::initHQ(Imgidx size, size_t listsize)
{
    this->maxSize_queue = size;
    queue = new HQueue_l1idx_rank(size);
    list = (Imgidx*)Malloc((listsize + 1) * sizeof(Imgidx));
    list[0] = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    minidx_queue = size;
}

HybridQueue_HQueue_Rank::HybridQueue_HQueue_Rank(Imgidx size)
{
    initHQ(size, LISTSIZE_DEFAULT);
}
HybridQueue_HQueue_Rank::HybridQueue_HQueue_Rank(Imgidx size, size_t listsize)
{
    initHQ(size, listsize);
}
HybridQueue_HQueue_Rank::~HybridQueue_HQueue_Rank()
{
    delete queue;
    Free(list - 1);
}
void HybridQueue_HQueue_Rank::push(Imgidx idx)
{
    _int16 i;
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
void HybridQueue_HQueue_Rank::pop()
{
    if (curSize_list == 0)
    {
        list[0] = queue->top();

        pop_queue();
    }
    else
    {
        for (int i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }
}

// HybridQueue_HQueue_Rank end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HybridQueue_HQueue_Rank1

void HybridQueue_HQueue_Rank1::initHQ(Imgidx size, size_t listsize)
{
    this->maxSize_queue = size;
    queue = new HQueue_l1idx_rank(size);

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
    minidx_queue = size;
}

HybridQueue_HQueue_Rank1::HybridQueue_HQueue_Rank1(Imgidx size)
{
    initHQ(size, LISTSIZE_DEFAULT);
}

HybridQueue_HQueue_Rank1::HybridQueue_HQueue_Rank1(Imgidx size, size_t listsize)
{
    initHQ(size, listsize);
}

HybridQueue_HQueue_Rank1::~HybridQueue_HQueue_Rank1()
{
    delete queue;
    Free(list);
}

void HybridQueue_HQueue_Rank1::push(Imgidx idx)
{
    _int16 i, j, lm;
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

void HybridQueue_HQueue_Rank1::pop()
{
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
}

// HybridQueue_HQueue_Rank1 end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HybridQueue_HQueue

void HybridQueue_HQueue::initHQ(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
{
    this->maxSize_queue = qsize_in;
    queue = new HQueue_l1idx(qsize_in, dhist, numlevels);
    list = (Imgidx*)Malloc((listsize) * sizeof(Imgidx));
    levels = (_int64*)Malloc((listsize + 1) * sizeof(_int64));
    levels[0] = 0;
    levels++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    minidx_queue = qsize_in;
    minlevnotfixed = 0;
}

HybridQueue_HQueue::HybridQueue_HQueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels)
{
    initHQ(qsize_in, dhist, numlevels, LISTSIZE_DEFAULT);
}

HybridQueue_HQueue::HybridQueue_HQueue(_uint64 qsize_in, Imgidx *dhist, _int32 numlevels, size_t listsize)
{
    initHQ(qsize_in, dhist, numlevels, listsize);
}

HybridQueue_HQueue::~HybridQueue_HQueue()
{
    delete queue;
    Free(list);
    Free(levels - 1);
}

void HybridQueue_HQueue::push(Imgidx idx, _int64 level)
{
    _int16 i;
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

Imgidx HybridQueue_HQueue::pop()
{
    Imgidx ret;

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
}

// HybridQueue_HQueue end
//////////////////////////////////////////////////////////////////////////////////////////////////////

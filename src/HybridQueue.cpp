#include "HybridQueue.hpp"

#if PROFILE
#include <fstream>
#endif

template <class Pixel>
void HierarHeapQueue_HEQ<Pixel>::initHQ(ImgIdx *dhist, uint32_t *histeqmap_in, ImgIdx numlevels_in, ImgIdx size,
                                        double a_in, int listsize) {
    list = (HQentry<Pixel> *)Malloc(listsize * sizeof(HQentry<Pixel>));
    maxSize_list = listsize - 1;
    curSize_list = -1;

    this->numlevels = numlevels_in;
    this->a = a_in;
    this->queue_minlev = numlevels;
    this->histeqmap = histeqmap_in; // do not free dhist outside
    this->qsizes = dhist;           // do not free dhist outside

    storage_cursize = (ImgIdx *)Calloc(numlevels * sizeof(ImgIdx));

    ImgIdx cumsum = 0;
    ImgIdx thr_nonredundantnodes = (ImgIdx)(size * 0.6);
    for (int level = 0; level < numlevels; level++) {
        cumsum += qsizes[level];
        if (cumsum > thr_nonredundantnodes) {
            thr_hqueue = curthr = level;
            break;
        }
    }

    hqueue = (HeapQueue_naive_quad<Pixel> **)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Pixel> *));
    for (int level = 0; level < thr_hqueue; level++)
        hqueue[level] = new HeapQueue_naive_quad<Pixel>(qsizes[level]);

    storage = (HQentry<Pixel> **)Calloc((numlevels - thr_hqueue) * sizeof(HQentry<Pixel> *));
    storage -= thr_hqueue;
    for (int level = thr_hqueue; level < numlevels; level++)
        storage[level] = (HQentry<Pixel> *)Malloc(qsizes[level] * sizeof(HQentry<Pixel>));
}

template <class Pixel>
HierarHeapQueue_HEQ<Pixel>::HierarHeapQueue_HEQ(ImgIdx *dhist, uint32_t *histeqmap_in, ImgIdx numlevels_in, double a_in,
                                                ImgIdx size) {
    initHQ(dhist, histeqmap_in, numlevels_in, size, a_in, 12);
}

template <class Pixel>
HierarHeapQueue_HEQ<Pixel>::HierarHeapQueue_HEQ(ImgIdx *dhist, uint32_t *histeqmap_in, ImgIdx numlevels_in, ImgIdx size,
                                                double a_in, int listsize) {
    initHQ(dhist, histeqmap_in, numlevels_in, size, a_in, listsize);
}

template <class Pixel> HierarHeapQueue_HEQ<Pixel>::~HierarHeapQueue_HEQ() {
    Free(list);
    Free(qsizes);
    Free(histeqmap);
    Free(storage_cursize);

    for (int level = 0; level < numlevels; level++)
        if (hqueue[level])
            delete hqueue[level];
    Free(hqueue);

    for (int level = thr_hqueue; level < numlevels; level++)
        if (storage[level])
            Free(storage[level]);
    Free(storage + thr_hqueue);
}

template <class Pixel> void HierarHeapQueue_HEQ<Pixel>::push_1stitem(ImgIdx idx, Pixel alpha) {
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template <class Pixel> void HierarHeapQueue_HEQ<Pixel>::endPushes(uint8_t *isVisited) {
    if (_emptyTop)
        pop(isVisited);
}

template <class Pixel> void HierarHeapQueue_HEQ<Pixel>::push(ImgIdx idx, Pixel alpha) {
    if (_emptyTop && alpha < list[0].alpha) {
        _emptyTop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }

    bool pushToList = (queue_minlev < curthr) ? alpha < hqueue[queue_minlev]->top_alpha()
                                              : (ImgIdx)histeqmap[(int)(a * log2(1 + (double)alpha))] < queue_minlev;

    if (pushToList) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            int i;
            for (i = curSize_list; alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
            curSize_list++;
        } else if (alpha < list[curSize_list].alpha) // push to the full list
        {
            push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

            int i;
            for (i = curSize_list - 1; alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
        } else
            push_queue(idx, alpha); // push to the queue
    } else
        push_queue(idx, alpha); // push to the queue
}

template <class Pixel> void HierarHeapQueue_HEQ<Pixel>::push_queue(ImgIdx idx, Pixel alpha) {
    int level = (int)histeqmap[(int)(a * log2(1 + (double)alpha))];

    // hidx = (int)(log2(1 + pow((double)dimg[dimgidx++],a)));
    if (level < queue_minlev)
        queue_minlev = level;

    if (level < curthr)
        hqueue[level]->push(idx, alpha);
    else {
        ImgIdx cur = storage_cursize[level]++;
        storage[level][cur].pidx = idx;
        storage[level][cur].alpha = alpha;
    }
}

template <class Pixel> ImgIdx HierarHeapQueue_HEQ<Pixel>::pop(uint8_t *isVisited) {
    ImgIdx ret = top();
    if (curSize_list == 0) {
        while (!check_queue_level(isVisited))
            queue_minlev++;
        list[0].pidx = hqueue[queue_minlev]->top();
        list[0].alpha = hqueue[queue_minlev]->top_alpha();

        pop_queue(isVisited);
    } else {
        for (int i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }
    return ret;
}

template <class Pixel> int HierarHeapQueue_HEQ<Pixel>::check_queue_level(uint8_t *isVisited) {
    if (queue_minlev < curthr)
        return hqueue[queue_minlev]->get_cursize();
    else {
        // while(curthr < queue_minlev)
        //{
        //	printf("Piep ");
        //	hqueue[curthr] = new HeapQueue_naive_quad<Pixel>(qsizes[curthr]);

        //	Free(storage[curthr]);
        //	storage[curthr] = 0;
        //	curthr++;
        //}
        curthr++;

        hqueue[queue_minlev] = new HeapQueue_naive_quad<Pixel>(qsizes[queue_minlev]);

        HQentry<Pixel> *store = storage[queue_minlev];
        ImgIdx cur = storage_cursize[queue_minlev];
        HeapQueue_naive_quad<Pixel> *pQ = hqueue[queue_minlev];
        for (ImgIdx p = 0; p < cur; p++) {
            if (!isVisited[store[p].pidx])
                pQ->push(store[p].pidx, store[p].alpha);
        }

        Free(storage[queue_minlev]);
        storage[queue_minlev] = 0;

        return pQ->get_cursize();
    }
}

template <class Pixel> void HierarHeapQueue_HEQ<Pixel>::pop_queue(uint8_t *isVisited) {
    hqueue[queue_minlev]->pop();

    if (!hqueue[queue_minlev]->get_cursize()) {
        do {
            queue_minlev++;
        } while (queue_minlev < numlevels && !check_queue_level(isVisited));
    }
}

template class HierarHeapQueue_HEQ<uint8_t>;
template class HierarHeapQueue_HEQ<uint16_t>;
template class HierarHeapQueue_HEQ<uint32_t>;
template class HierarHeapQueue_HEQ<uint64_t>;

template <class Pixel>
HierarHeapQueue<Pixel>::HierarHeapQueue(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize,
                                        int connectivity, double r) {
    this->numlevels = numlevels_in;
    this->a = a_in;
    this->queue_minlev = numlevels;

    ImgIdx cumsum = 0;
    qsizes = dhist; // do not free dhist outside
    if (r >= 1) {
        thr_hqueue = curthr = numlevels;
        hqueue = (HeapQueue_naive_quad<Pixel> **)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Pixel> *));
        for (int level = 0; level < thr_hqueue; level++)
            hqueue[level] = new HeapQueue_naive_quad<Pixel>(qsizes[level]);
        storage = 0;
        storage_cursize = 0;
    } else {
        storage_cursize = (ImgIdx *)Calloc(numlevels * sizeof(ImgIdx));
        ImgIdx thr_nonredundantnodes = (ImgIdx)(size * r);
        for (int level = 0; level < numlevels; level++) {
            cumsum += qsizes[level];
            if (cumsum > thr_nonredundantnodes) {
                thr_hqueue = curthr = level;
                break;
            }
        }

        hqueue = (HeapQueue_naive_quad<Pixel> **)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Pixel> *));
        for (int level = 0; level < thr_hqueue; level++)
            hqueue[level] = new HeapQueue_naive_quad<Pixel>(qsizes[level]);
        // hqueue[numlevels] = new HeapQueue_naive_quad<Pixel>(1);
        // hqueue[numlevels]->push(0,(Pixel)-1);

        storage = (HQentry<Pixel> **)Calloc((numlevels - thr_hqueue) * sizeof(HQentry<Pixel> *));
        storage -= thr_hqueue;
        for (int level = thr_hqueue; level < numlevels; level++)
            storage[level] = (HQentry<Pixel> *)Malloc(qsizes[level] * sizeof(HQentry<Pixel>));
    }
}

template <class Pixel> HierarHeapQueue<Pixel>::~HierarHeapQueue() {
    Free(qsizes);

    for (int level = 0; level < numlevels; level++)
        if (hqueue[level])
            delete hqueue[level];
    Free(hqueue);

    if (storage) {
        Free(storage_cursize);
        for (int level = thr_hqueue; level < numlevels; level++)
            if (storage[level])
                Free(storage[level]);
        Free(storage + thr_hqueue);
    }
}

template <class Pixel> void HierarHeapQueue<Pixel>::push_queue(ImgIdx idx, Pixel alpha) {
    int level = (int)(a * log2(1 + (double)alpha));
    if (level < queue_minlev)
        queue_minlev = level;

    if (level < curthr)
        hqueue[level]->push(idx, alpha);
    else {
        ImgIdx cur = storage_cursize[level]++;
        storage[level][cur].pidx = idx;
        storage[level][cur].alpha = alpha;
    }
}

template <class Pixel> ImgIdx HierarHeapQueue<Pixel>::pop(uint8_t *isVisited) {
    ImgIdx ret = top();

    pop_queue(isVisited);
    return ret;
}

template <class Pixel> int HierarHeapQueue<Pixel>::check_queue_level(uint8_t *isVisited) {
    if (queue_minlev < curthr)
        return hqueue[queue_minlev]->get_cursize();
    else {
        while (curthr < queue_minlev) {
            hqueue[curthr] = new HeapQueue_naive_quad<Pixel>(qsizes[curthr]);

            Free(storage[curthr]);
            storage[curthr] = 0;
            curthr++;
        }
        curthr++;

        hqueue[queue_minlev] = new HeapQueue_naive_quad<Pixel>(qsizes[queue_minlev]);

        HQentry<Pixel> *store = storage[queue_minlev];
        ImgIdx cur = storage_cursize[queue_minlev];
        HeapQueue_naive_quad<Pixel> *pQ = hqueue[queue_minlev];
        for (ImgIdx p = 0; p < cur; p++) {
            if (!isVisited[store[p].pidx])
                pQ->push(store[p].pidx, store[p].alpha);
        }

        Free(storage[queue_minlev]);
        storage[queue_minlev] = 0;

        return pQ->get_cursize();
    }
}

template <class Pixel> void HierarHeapQueue<Pixel>::pop_queue(uint8_t *isVisited) {
    hqueue[queue_minlev]->pop();

    if (!hqueue[queue_minlev]->get_cursize()) {
        do {
            queue_minlev++;
        } while (queue_minlev < numlevels && !check_queue_level(isVisited));
    }
}

template class HierarHeapQueue<uint8_t>;
template class HierarHeapQueue<uint16_t>;
template class HierarHeapQueue<uint32_t>;
template class HierarHeapQueue<uint64_t>;

//////////////////////////////////////////////////////////////////////////////////////////////////////
// HierarHeapQueue_cache

// hhpq

template <class Pixel>
void HierarHeapQueue_cache<Pixel>::initHQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize,
                                          int connectivity, double r) {
    maxSize = size;

    list = (HQentry<Pixel> *)Malloc((listsize) * sizeof(HQentry<Pixel>));
    maxSize_list = listsize - 1;
    curSize_list = -1;

    this->numlevels = numlevels_in;
    this->a = a_in;
    this->queue_minlev = numlevels;

    ImgIdx cumsum = 0;
    qsizes = dhist; // do not free dhist outside
    if (r >= 1) {
        thr_hqueue = curthr = numlevels;
        hqueue = (HeapQueue_naive_quad<Pixel> **)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Pixel> *));
        for (int level = 0; level < thr_hqueue; level++)
            hqueue[level] = new HeapQueue_naive_quad<Pixel>(qsizes[level]);
        storage = 0;
        storage_cursize = 0;
    } else {
        storage_cursize = (ImgIdx *)Calloc(numlevels * sizeof(ImgIdx));
        ImgIdx thr_nonredundantnodes = (ImgIdx)(size * r);
        for (int level = 0; level < numlevels; level++) {
            cumsum += qsizes[level];
            if (cumsum > thr_nonredundantnodes) {
                thr_hqueue = curthr = level;
                break;
            }
        }

        hqueue = (HeapQueue_naive_quad<Pixel> **)Calloc(numlevels * sizeof(HeapQueue_naive_quad<Pixel> *));
        for (int level = 0; level < thr_hqueue; level++)
            hqueue[level] = new HeapQueue_naive_quad<Pixel>(qsizes[level]);

        storage = (HQentry<Pixel> **)Calloc((numlevels - thr_hqueue) * sizeof(HQentry<Pixel> *));
        storage -= thr_hqueue;
        for (int level = thr_hqueue; level < numlevels; level++)
            storage[level] = (HQentry<Pixel> *)Malloc(qsizes[level] * sizeof(HQentry<Pixel>));
    }
}

template <class Pixel>
HierarHeapQueue_cache<Pixel>::HierarHeapQueue_cache(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in,
                                                    int listsize, ImgIdx connectivity, double r) {
    initHQ(dhist, numlevels_in, size, a_in, listsize, (int)connectivity, r);
}

template <class Pixel> HierarHeapQueue_cache<Pixel>::~HierarHeapQueue_cache() {
#if PROFILE
    std::ofstream outFile("MemmoveHHPQ.txt", std::ios::app);
    for (int i = 0; i < (int)num_memmove_push.size(); i++) {
        outFile << "0 " << num_memmove_push[i] << " " << num_items_push[i] << std::endl;
        // printf("%d %d\n", (int)num_memmove_push[i], (int)num_items_push[i]);
    }
    for (int i = 0; i < (int)num_memmove_pop.size(); i++) {
        outFile << "1 " << num_memmove_pop[i] << " " << num_items_pop[i] << std::endl;
        // printf("%d %d\n", (int)num_memmove_pop[i], (int)num_items_pop[i]);
    }
    outFile.close();
#endif

    Free(list);
    Free(qsizes);

    for (int level = 0; level < numlevels; level++)
        if (hqueue[level])
            delete hqueue[level];
    Free(hqueue);

    if (storage) {
        Free(storage_cursize);
        for (int level = thr_hqueue; level < numlevels; level++)
            if (storage[level])
                Free(storage[level]);
        Free(storage + thr_hqueue);
    }

#if PROFILE
    double t1 = get_cpu_time() - t0;
    double sz = (double)maxSize;
    printf("HQ Profile - Total: %f, Cache: %f, Queue: %f, Conv: %f\n", t1, tcache, tqueue, tconv);
    printf("Cached: %d(%f), C. ovfl: %d(%f), Heap Queued: %d(%f), Stored: %d(%f), Converted: %d(%f), \n",
           (int)num_cache, (double)num_cache / sz, (int)num_cache_ovfl, (double)num_cache_ovfl / sz, (int)num_hq,
           (double)num_hq / sz, (int)num_store, (double)num_store / sz, (int)num_conv, (double)num_conv / sz);
    printf("Initial: %d queues + %d storages (%d total) End: %d queues + %d storages (%d total)\n", (int)thr_hqueue,
           (int)(numlevels - thr_hqueue), (int)numlevels, (int)curthr, (int)(numlevels - curthr), (int)numlevels);

#endif
}

template <class Pixel> void HierarHeapQueue_cache<Pixel>::push_1stitem(ImgIdx idx, Pixel alpha) {
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;

#if PROFILE
    // printf("push %d, size %d\n", idx, curSize);
    curSize = 1;
    num_memmove_push.push_back(1);
    num_items_push.push_back(1);
#endif
}

template <class Pixel> void HierarHeapQueue_cache<Pixel>::endPushes(uint8_t *isVisited) {
    if (_emptyTop)
        pop(isVisited);
#if PROFILE
    else
        decrement_curSize();
#endif
}

template <class Pixel> void HierarHeapQueue_cache<Pixel>::push(ImgIdx idx, Pixel alpha) {
#if PROFILE
    // printf("push %d, size %d\n", idx, curSize);
    curSize++;
    num_memmove_push_i = 0;
#endif
    if (_emptyTop && alpha < list[0].alpha) {
#if PROFILE
        num_cache++;
        num_memmove_push.push_back(1);
        num_items_push.push_back(curSize);
#endif
        _emptyTop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }

    bool pushToList = (queue_minlev < curthr) ? alpha < hqueue[queue_minlev]->top_alpha()
                                              : (int)(a * log2(1 + (double)alpha)) < queue_minlev;

    if (pushToList) {
#if PROFILE
        num_cache++;
        double t1 = get_cpu_time(), t2, tq = 0;
#endif
        if (curSize_list < maxSize_list) // spare room in the list
        {
            int i;
            for (i = curSize_list; i >= 0 && alpha < list[i].alpha; i--) {
                list[i + 1] = list[i];
#if PROFILE
                num_memmove_push_i++;
#endif
            }
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
#if PROFILE
            num_memmove_push_i++;
#endif
            curSize_list++;
        } else if (alpha < list[curSize_list].alpha) // push to the full list
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
            for (i = curSize_list - 1; i >= 0 && alpha < list[i].alpha; i--) {
                list[i + 1] = list[i];
#if PROFILE
                num_memmove_push_i++;
#endif
            }
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
#if PROFILE
            num_memmove_push_i++;
#endif
        } else {
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
    } else {
#if PROFILE
        double t1 = get_cpu_time();
#endif
        push_queue(idx, alpha); // push to the queue
#if PROFILE
        tqueue += get_cpu_time() - t1;
#endif
    }
#if PROFILE
    num_memmove_push.push_back(num_memmove_push_i);
    num_items_push.push_back(curSize);
#endif
}

template <class Pixel> void HierarHeapQueue_cache<Pixel>::push_queue(ImgIdx idx, Pixel alpha) {
    int level = (int)(a * log2(1 + (double)alpha));

    if (level < queue_minlev)
        queue_minlev = level;

    if (level < curthr) {
#if PROFILE
        num_hq++;
        num_memmove_push_i += hqueue[level]->push(idx, alpha);
#else
        hqueue[level]->push(idx, alpha);
#endif
    } else {
#if PROFILE
        num_store++;
        num_memmove_push_i++;
#endif
        ImgIdx cur = storage_cursize[level]++;
        storage[level][cur].pidx = idx;
        storage[level][cur].alpha = alpha;
    }
}

template <class Pixel> ImgIdx HierarHeapQueue_cache<Pixel>::pop(uint8_t *isVisited) {
#if PROFILE
    // printf("pop %d, size %d\n", top(), curSize);
    curSize--;
    num_memmove_pop_i = 0;
#endif
    ImgIdx ret = top();
    if (curSize_list == 0) {
        while (!check_queue_level(isVisited))
            queue_minlev++;
        list[0].pidx = hqueue[queue_minlev]->top();
        list[0].alpha = hqueue[queue_minlev]->top_alpha();
#if PROFILE
        num_memmove_pop_i++;
#endif

        pop_queue(isVisited);
    } else {
        for (int i = 0; i < curSize_list; i++) {
            list[i] = list[i + 1];
#if PROFILE
            num_memmove_pop_i++;
#endif
        }
        curSize_list--;
    }
#if PROFILE
    num_memmove_pop.push_back(num_memmove_pop_i);
    num_items_pop.push_back(curSize);
#endif
    return ret;
}

template <class Pixel> int HierarHeapQueue_cache<Pixel>::check_queue_level(uint8_t *isVisited) {
#if PROFILE
#endif
    if (queue_minlev < curthr)
        return hqueue[queue_minlev]->get_cursize();
    else {
#if PROFILE
        double t1 = get_cpu_time();
#endif
        while (curthr < queue_minlev) {
            hqueue[curthr] = new HeapQueue_naive_quad<Pixel>(qsizes[curthr]);

            Free(storage[curthr]);
            storage[curthr] = 0;
            curthr++;
        }
        curthr++;

        hqueue[queue_minlev] = new HeapQueue_naive_quad<Pixel>(qsizes[queue_minlev]);

        HQentry<Pixel> *store = storage[queue_minlev];
        ImgIdx cur = storage_cursize[queue_minlev];
        HeapQueue_naive_quad<Pixel> *pQ = hqueue[queue_minlev];
        for (ImgIdx p = 0; p < cur; p++) {
            if (!isVisited[store[p].pidx])
#if PROFILE
                num_memmove_pop_i += pQ->push(store[p].pidx, store[p].alpha);
#else
                pQ->push(store[p].pidx, store[p].alpha);
#endif
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

template <class Pixel> void HierarHeapQueue_cache<Pixel>::pop_queue(uint8_t *isVisited) {

#if PROFILE
    double t1 = get_cpu_time();
    num_memmove_pop_i += hqueue[queue_minlev]->pop();
#else
    hqueue[queue_minlev]->pop();
#endif

#if PROFILE
    tqueue += get_cpu_time() - t1;
#endif

    if (!hqueue[queue_minlev]->get_cursize()) {
        do {
            queue_minlev++;
        } while (queue_minlev < numlevels && !check_queue_level(isVisited));
    }
}

template class HierarHeapQueue_cache<uint8_t>;
template class HierarHeapQueue_cache<uint16_t>;
template class HierarHeapQueue_cache<uint32_t>;
template class HierarHeapQueue_cache<uint64_t>;

// HierarHeapQueue_cache end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Cache_Heapqueue

template <class Pixel> void Cache_Heapqueue<Pixel>::initHQ(ImgIdx size, size_t listsize) {
    this->maxSize_queue = size;
    hqueue = new HeapQueue_naive<Pixel>(size);
    list = (HQentry<Pixel> *)Malloc((listsize + 1) * sizeof(HQentry<Pixel>));
    list[0].pidx = 0;
    list[0].alpha = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    qtime = 0; // tmp
}

template <class Pixel> Cache_Heapqueue<Pixel>::Cache_Heapqueue(ImgIdx size) { initHQ(size, 12); }

template <class Pixel> Cache_Heapqueue<Pixel>::Cache_Heapqueue(ImgIdx size, size_t listsize) { initHQ(size, listsize); }

template <class Pixel> Cache_Heapqueue<Pixel>::~Cache_Heapqueue() {
    delete hqueue;
    Free(list - 1);
}

template <class Pixel> void Cache_Heapqueue<Pixel>::push_1stitem(ImgIdx idx, Pixel alpha) {
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template <class Pixel> void Cache_Heapqueue<Pixel>::push(ImgIdx idx, Pixel alpha) {
    int16_t i;
    if (_emptyTop && alpha < list[0].alpha) {
        _emptyTop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }

    if (alpha < hqueue->top_alpha()) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            for (i = curSize_list; alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
            curSize_list++;
        } else if (alpha < list[curSize_list].alpha) // push to the full list
        {
            push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

            for (i = curSize_list - 1; alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
        } else
            push_queue(idx, alpha); // push to the queue
    } else
        push_queue(idx, alpha); // push to the queue
}

template <class Pixel> ImgIdx Cache_Heapqueue<Pixel>::pop() {
    ImgIdx ret = top();
    int8_t i;
    if (curSize_list == 0) {
        list[0].pidx = hqueue->top();
        list[0].alpha = hqueue->top_alpha();

        pop_queue();
    } else {
        for (i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }

    return ret;
}

template class Cache_Heapqueue<uint8_t>;
template class Cache_Heapqueue<uint16_t>;
template class Cache_Heapqueue<uint32_t>;
template class Cache_Heapqueue<uint64_t>;
template class Cache_Heapqueue<float>;
template class Cache_Heapqueue<double>;

// Cache_Heapqueue end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Cache_Quad_Heapqueue

template <class Pixel> void Cache_Quad_Heapqueue<Pixel>::initHQ(ImgIdx size, size_t listsize) {
    this->maxSize_queue = size;
    hqueue = new HeapQueue_naive_quad<Pixel>(size);
    list = (HQentry<Pixel> *)Malloc((listsize + 1) * sizeof(HQentry<Pixel>));
    list[0].pidx = 0;
    list[0].alpha = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    qtime = 0; // tmp
}

template <class Pixel> Cache_Quad_Heapqueue<Pixel>::Cache_Quad_Heapqueue(ImgIdx size) { initHQ(size, 12); }

template <class Pixel> Cache_Quad_Heapqueue<Pixel>::Cache_Quad_Heapqueue(ImgIdx size, size_t listsize) {
    initHQ(size, listsize);
}

template <class Pixel> Cache_Quad_Heapqueue<Pixel>::~Cache_Quad_Heapqueue() {
    delete hqueue;
    Free(list - 1);
}

template <class Pixel> void Cache_Quad_Heapqueue<Pixel>::push_1stitem(ImgIdx idx, Pixel alpha) {
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template <class Pixel> void Cache_Quad_Heapqueue<Pixel>::push(ImgIdx idx, Pixel alpha) {
    int16_t i;

    if (_emptyTop && alpha < list[0].alpha) {
        _emptyTop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }
    if (alpha < hqueue->top_alpha()) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            for (i = curSize_list; i >= 0 && alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
            curSize_list++;
        } else if (alpha < list[curSize_list].alpha) // push to the full list
        {
            push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

            for (i = curSize_list - 1; i >= 0 && alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
        } else
            push_queue(idx, alpha); // push to the queue
    } else
        push_queue(idx, alpha); // push to the queue
}

template <class Pixel> ImgIdx Cache_Quad_Heapqueue<Pixel>::pop() {
    ImgIdx ret = top();
    int8_t i;

    if (curSize_list == 0) {
        list[0].pidx = hqueue->top();
        list[0].alpha = hqueue->top_alpha();

        pop_queue();
    } else {
        for (i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }
    return ret;
}

template class Cache_Quad_Heapqueue<uint8_t>;
template class Cache_Quad_Heapqueue<uint16_t>;
template class Cache_Quad_Heapqueue<uint32_t>;
template class Cache_Quad_Heapqueue<uint64_t>;
template class Cache_Quad_Heapqueue<float>;
template class Cache_Quad_Heapqueue<double>;

// Cache_Quad_Heapqueue end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// CirCache_Hierqueue

template <class Pixel>
void CirCache_Hierqueue<Pixel>::initHQ(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize) {
    int shamt = 0;
    for (int lsize = listsize; lsize; lsize >>= 1)
        shamt++;
    listsize = 1 << (shamt - 1);
    mask = listsize - 1;
    this->maxSize_queue = qsize_in;
    hqueue = new HierarQueue(qsize_in, dhist, numlevels);
    list = (HQentry<int32_t> *)Malloc(listsize * sizeof(HQentry<int32_t>));
    maxSize_list = listsize;
    curSize_list = 0;
    liststart = 0;
    qtime = 0; // tmp
}

template <class Pixel>
CirCache_Hierqueue<Pixel>::CirCache_Hierqueue(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels) {
    initHQ(qsize_in, dhist, numlevels, 16);
}

template <class Pixel>
CirCache_Hierqueue<Pixel>::CirCache_Hierqueue(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize) {
    initHQ(qsize_in, dhist, numlevels, listsize);
}

template <class Pixel> CirCache_Hierqueue<Pixel>::~CirCache_Hierqueue() {
    delete hqueue;
    Free(list);
}

template <class Pixel> void CirCache_Hierqueue<Pixel>::push_1stitem(ImgIdx idx, int32_t alpha) {
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template <class Pixel> void CirCache_Hierqueue<Pixel>::push(ImgIdx idx, int32_t alpha) {
    int16_t i, j, k;

    if (_emptyTop && alpha < list[liststart].alpha) {
        _emptyTop = 0;
        list[liststart].pidx = idx;
        list[liststart].alpha = alpha;
        return;
    }
    if ((int64_t)alpha < hqueue->get_minlev()) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            j = (liststart - 1) & mask;
            for (i = (liststart + curSize_list - 1) & mask; i != j && alpha < list[i].alpha; i = (i - 1) & mask)
                list[(i + 1) & mask] = list[i];
            list[(i + 1) & mask].pidx = idx;
            list[(i + 1) & mask].alpha = alpha;
            curSize_list++;
        } else if (alpha < list[j = ((liststart + curSize_list - 1) & mask)].alpha) // push to the full list
        {
            push_queue(list[j].pidx, list[j].alpha);

            k = (liststart - 1) & mask;
            for (i = (j - 1) & mask; i != k && alpha < list[i].alpha; i = (i - 1) & mask)
                list[(i + 1) & mask] = list[i];
            list[(i + 1) & mask].pidx = idx;
            list[(i + 1) & mask].alpha = alpha;
        } else
            push_queue(idx, alpha); // push to the queue
    } else
        push_queue(idx, alpha); // push to the queue
}

template <class Pixel> ImgIdx CirCache_Hierqueue<Pixel>::pop() {
    ImgIdx ret = top();

    if (curSize_list == 1) {
        list[liststart].pidx = hqueue->top();
        list[liststart].alpha = hqueue->get_minlev();

        pop_queue();
    } else {
        liststart = (liststart + 1) & mask;
        curSize_list--;
    }

    return ret;
}

template class CirCache_Hierqueue<uint8_t>;
template class CirCache_Hierqueue<uint16_t>;
template class CirCache_Hierqueue<uint32_t>;
template class CirCache_Hierqueue<uint64_t>;
template class CirCache_Hierqueue<float>;
template class CirCache_Hierqueue<double>;

// CirCache_Hierqueue end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HierarQueueCache

template <class Pixel>
void HierarQueueCache<Pixel>::initHQ(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize) {
    this->maxSize_queue = qsize_in;
    hqueue = new HierarQueue(qsize_in, dhist, numlevels);
    list = (HQentry<int32_t> *)Malloc((listsize + 1) * sizeof(HQentry<int32_t>));
    list[0].pidx = 0;
    list[0].alpha = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    qtime = 0; // tmp
}

template <class Pixel> HierarQueueCache<Pixel>::HierarQueueCache(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels) {
    initHQ(qsize_in, dhist, numlevels, 12);
}

template <class Pixel>
HierarQueueCache<Pixel>::HierarQueueCache(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize) {
    initHQ(qsize_in, dhist, numlevels, listsize);
}

template <class Pixel> HierarQueueCache<Pixel>::~HierarQueueCache() {
    delete hqueue;
    Free(list - 1);
}

template <class Pixel> void HierarQueueCache<Pixel>::push_1stitem(ImgIdx idx, int32_t alpha) {
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template <class Pixel> void HierarQueueCache<Pixel>::push(ImgIdx idx, int32_t alpha) {
    int16_t i;

    if (_emptyTop && alpha < list[0].alpha) {
        _emptyTop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }
    if ((int64_t)alpha < hqueue->get_minlev()) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            for (i = curSize_list; alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
            curSize_list++;
        } else if (alpha < list[curSize_list].alpha) // push to the full list
        {
            push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

            for (i = curSize_list - 1; alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
        } else
            push_queue(idx, alpha); // push to the queue
    } else
        push_queue(idx, alpha); // push to the queue
}

template <class Pixel> ImgIdx HierarQueueCache<Pixel>::pop() {
    ImgIdx ret = top();
    int8_t i;

    if (curSize_list == 0) {
        list[0].pidx = hqueue->top();
        list[0].alpha = hqueue->get_minlev();

        pop_queue();
    } else {
        for (i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }
    return ret;
}

template class HierarQueueCache<uint8_t>;
template class HierarQueueCache<uint16_t>;
template class HierarQueueCache<uint32_t>;
template class HierarQueueCache<uint64_t>;

// HierarQueueCache end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Cache_Hierqueue_l1

template <class Pixel>
void Cache_Hierqueue_l1<Pixel>::initHQ(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize) {
    this->maxSize_queue = qsize_in;
    hqueue = new HQueue_l1idx(qsize_in, dhist, numlevels);
    list = (HQentry<int32_t> *)Malloc((listsize + 1) * sizeof(HQentry<int32_t>));
    list[0].pidx = 0;
    list[0].alpha = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;

    qtime = 0; // tmp
}

template <class Pixel>
Cache_Hierqueue_l1<Pixel>::Cache_Hierqueue_l1(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels) {
    initHQ(qsize_in, dhist, numlevels, 12);
}

template <class Pixel>
Cache_Hierqueue_l1<Pixel>::Cache_Hierqueue_l1(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize) {
    initHQ(qsize_in, dhist, numlevels, listsize);
}

template <class Pixel> Cache_Hierqueue_l1<Pixel>::~Cache_Hierqueue_l1() {
    delete hqueue;
    Free(list - 1);
}

template <class Pixel> void Cache_Hierqueue_l1<Pixel>::push_1stitem(ImgIdx idx, int32_t alpha) {
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template <class Pixel> void Cache_Hierqueue_l1<Pixel>::push(ImgIdx idx, int32_t alpha) {
    int16_t i;

    if (_emptyTop && alpha < list[0].alpha) {
        _emptyTop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }

    if ((int64_t)alpha < hqueue->get_minlev()) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            for (i = curSize_list; alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
            curSize_list++;
        } else if (alpha < list[curSize_list].alpha) // push to the full list
        {
            push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

            for (i = curSize_list - 1; alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
        } else
            push_queue(idx, alpha); // push to the queue
    } else
        push_queue(idx, alpha); // push to the queue
}

template <class Pixel> ImgIdx Cache_Hierqueue_l1<Pixel>::pop() {
    ImgIdx ret = top();
    int8_t i;

    if (curSize_list == 0) {
        list[0].pidx = hqueue->top();
        list[0].alpha = hqueue->get_minlev();

        pop_queue();
    } else {
        for (i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }

    return ret;
}

template class Cache_Hierqueue_l1<uint8_t>;
template class Cache_Hierqueue_l1<uint16_t>;
template class Cache_Hierqueue_l1<uint32_t>;
template class Cache_Hierqueue_l1<uint64_t>;

// Cache_Hierqueue_l1 end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Cache_Hierqueue_l2

template <class Pixel>
void Cache_Hierqueue_l2<Pixel>::initHQ(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize) {
    this->maxSize_queue = qsize_in;
    hqueue = new HQueue_l2idx(qsize_in, dhist, numlevels);
    list = (HQentry<int32_t> *)Malloc((listsize + 1) * sizeof(HQentry<int32_t>));
    list[0].pidx = 0;
    list[0].alpha = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    qtime = 0; // tmp
}

template <class Pixel>
Cache_Hierqueue_l2<Pixel>::Cache_Hierqueue_l2(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels) {
    initHQ(qsize_in, dhist, numlevels, 12);
}

template <class Pixel>
Cache_Hierqueue_l2<Pixel>::Cache_Hierqueue_l2(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize) {
    initHQ(qsize_in, dhist, numlevels, listsize);
}

template <class Pixel> Cache_Hierqueue_l2<Pixel>::~Cache_Hierqueue_l2() {
    delete hqueue;
    Free(list - 1);
}

template <class Pixel> void Cache_Hierqueue_l2<Pixel>::push_1stitem(ImgIdx idx, int32_t alpha) {
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template <class Pixel> void Cache_Hierqueue_l2<Pixel>::push(ImgIdx idx, int32_t alpha) {
    int16_t i;

    if (_emptyTop && alpha < list[0].alpha) {
        _emptyTop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }

    if ((int64_t)alpha < hqueue->get_minlev()) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            for (i = curSize_list; alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
            curSize_list++;
        } else if (alpha < list[curSize_list].alpha) // push to the full list
        {
            push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

            for (i = curSize_list - 1; alpha < list[i].alpha; i--)
                list[i + 1] = list[i];
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
        } else
            push_queue(idx, alpha); // push to the queue
    } else
        push_queue(idx, alpha); // push to the queue
}

template <class Pixel> ImgIdx Cache_Hierqueue_l2<Pixel>::pop() {
    ImgIdx ret = top();
    int8_t i;

    if (curSize_list == 0) {
        list[0].pidx = hqueue->top();
        list[0].alpha = hqueue->get_minlev();

        pop_queue();
    } else {
        for (i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }
    return ret;
}

template class Cache_Hierqueue_l2<uint8_t>;
template class Cache_Hierqueue_l2<uint16_t>;
template class Cache_Hierqueue_l2<uint32_t>;
template class Cache_Hierqueue_l2<uint64_t>;

// Cache_Hierqueue_l2 end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// Trie_Cache

void Trie_Cache::initHQ(ImgIdx size, size_t listsize) {
    this->maxSize_queue = size;
    trie = new Trie<uint64_t>(size);
    list = (ImgIdx *)Malloc((listsize + 1) * sizeof(ImgIdx));
    list[0] = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    minidx_queue = size;
}

Trie_Cache::Trie_Cache(ImgIdx size) { initHQ(size, LISTSIZE_DEFAULT); }

Trie_Cache::Trie_Cache(ImgIdx size, size_t listsize) { initHQ(size, listsize); }

Trie_Cache::~Trie_Cache() {
    delete trie;
    Free(list - 1);
}

void Trie_Cache::push(ImgIdx idx) {
    int16_t i;
    if (idx < trie->top()) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            for (i = curSize_list; idx < list[i]; i--)
                list[i + 1] = list[i];
            list[i + 1] = idx;
            curSize_list++;
        } else if (idx < list[curSize_list]) // push to the full list
        {
            push_queue(list[curSize_list]);

            for (i = curSize_list - 1; idx < list[i]; i--)
                list[i + 1] = list[i];
            list[i + 1] = idx;
        } else
            push_queue(idx); // push to the queue
    } else
        push_queue(idx); // push to the queue
}

void Trie_Cache::pop() {
    int8_t i;
    if (curSize_list == 0) {
        list[0] = trie->top();

        pop_queue();
    } else {
        for (i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }
}

// Trie_Cache end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HybridQueue_HQueue_Rank

void HybridQueue_HQueue_Rank::initHQ(ImgIdx size, size_t listsize) {
    this->maxSize_queue = size;
    queue = new HQueue_l1idx_rank(size);
    list = (ImgIdx *)Malloc((listsize + 1) * sizeof(ImgIdx));
    list[0] = 0;
    list++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    minidx_queue = size;
}

HybridQueue_HQueue_Rank::HybridQueue_HQueue_Rank(ImgIdx size) { initHQ(size, LISTSIZE_DEFAULT); }
HybridQueue_HQueue_Rank::HybridQueue_HQueue_Rank(ImgIdx size, size_t listsize) { initHQ(size, listsize); }
HybridQueue_HQueue_Rank::~HybridQueue_HQueue_Rank() {
    delete queue;
    Free(list - 1);
}
void HybridQueue_HQueue_Rank::push(ImgIdx idx) {
    int16_t i;
    if (idx < queue->top()) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            for (i = curSize_list; idx < list[i]; i--)
                list[i + 1] = list[i];
            list[i + 1] = idx;
            curSize_list++;
        } else if (idx < list[curSize_list]) // push to the full list
        {
            push_queue(list[curSize_list]);

            for (i = curSize_list - 1; idx < list[i]; i--)
                list[i + 1] = list[i];
            list[i + 1] = idx;
        } else
            push_queue(idx); // push to the queue
    } else
        push_queue(idx); // push to the queue
}
void HybridQueue_HQueue_Rank::pop() {
    if (curSize_list == 0) {
        list[0] = queue->top();

        pop_queue();
    } else {
        for (int i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    }
}

// HybridQueue_HQueue_Rank end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HybridQueue_HQueue_Rank1

void HybridQueue_HQueue_Rank1::initHQ(ImgIdx size, size_t listsize) {
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
    list = (ImgIdx *)Malloc((listsize) * sizeof(ImgIdx));
    maxSize_list = listsize;
    curSize_list = 0;
    l0 = listsize >> 1;
    list[l0] = size;
    minidx_queue = size;
}

HybridQueue_HQueue_Rank1::HybridQueue_HQueue_Rank1(ImgIdx size) { initHQ(size, LISTSIZE_DEFAULT); }

HybridQueue_HQueue_Rank1::HybridQueue_HQueue_Rank1(ImgIdx size, size_t listsize) { initHQ(size, listsize); }

HybridQueue_HQueue_Rank1::~HybridQueue_HQueue_Rank1() {
    delete queue;
    Free(list);
}

void HybridQueue_HQueue_Rank1::push(ImgIdx idx) {
    int16_t i, j, lm;
    lm = (l0 - 1) & mask;
    if (idx < queue->top()) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            if (idx < list[l0]) {
                list[lm] = idx;
                l0 = lm;
            } else {
                i = (l0 + curSize_list) & mask;
                j = (l0 + curSize_list - 1) & mask;
                while (idx < list[j]) {
                    list[i] = list[j];
                    i = j;
                    j = (j - 1) & mask;
                }
                list[i] = idx;
            }
            curSize_list++;
        } else if (idx < list[curSize_list]) // push to the full list
        {
            push_queue(list[lm]);

            if (idx < list[l0]) {
                list[lm] = idx;
                l0 = lm;
            } else {
                i = lm;
                j = (lm - 1) & mask;
                while (idx < list[j]) {
                    list[i] = list[j];
                    i = j;
                    j = (j - 1) & mask;
                }
                list[i] = idx;
            }
        } else
            push_queue(idx); // push to the queue
    } else
        push_queue(idx); // push to the queue
}

void HybridQueue_HQueue_Rank1::pop() {
    if (curSize_list == 1) {
        list[l0] = queue->top();

        pop_queue();
    } else {
        l0 = (l0 + 1) & mask;
        curSize_list--;
    }
}

// HybridQueue_HQueue_Rank1 end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HybridQueue_HQueue

void HybridQueue_HQueue::initHQ(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize) {
    this->maxSize_queue = qsize_in;
    queue = new HQueue_l1idx(qsize_in, dhist, numlevels);
    list = (ImgIdx *)Malloc((listsize) * sizeof(ImgIdx));
    levels = (int64_t *)Malloc((listsize + 1) * sizeof(int64_t));
    levels[0] = 0;
    levels++;
    maxSize_list = listsize - 1;
    curSize_list = -1;
    minidx_queue = qsize_in;
    minlevnotfixed = 0;
}

HybridQueue_HQueue::HybridQueue_HQueue(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels) {
    initHQ(qsize_in, dhist, numlevels, LISTSIZE_DEFAULT);
}

HybridQueue_HQueue::HybridQueue_HQueue(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels, size_t listsize) {
    initHQ(qsize_in, dhist, numlevels, listsize);
}

HybridQueue_HQueue::~HybridQueue_HQueue() {
    delete queue;
    Free(list);
    Free(levels - 1);
}

void HybridQueue_HQueue::push(ImgIdx idx, int64_t level) {
    int16_t i;
    if (curSize_list == -1) // should be run only the first time
    {
        list[0] = idx;
        levels[0] = level;
        curSize_list = 0;
        push_queue(idx, level);
        return;
    }

    if (level < queue->get_minlev()) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            for (i = curSize_list; level < levels[i]; i--) {
                list[i + 1] = list[i];
                levels[i + 1] = levels[i];
            }
            list[i + 1] = idx;
            levels[i + 1] = level;
            curSize_list++;
        } else if (level < levels[curSize_list]) // push to the full list
        {
            push_queue(list[curSize_list], levels[curSize_list]);

            for (i = curSize_list - 1; level < levels[i]; i--) {
                list[i + 1] = list[i];
                levels[i + 1] = levels[i];
            }
            list[i + 1] = idx;
            levels[i + 1] = level;
        } else
            push_queue(idx, level); // push to the queue
    } else
        push_queue(idx, level); // push to the queue
}

ImgIdx HybridQueue_HQueue::pop() {
    ImgIdx ret;

    ret = list[0];

    if (curSize_list == 0) {
        list[0] = queue->top();
        levels[0] = queue->get_minlev();

        pop_queue();
    } else {
        for (int8_t i = 0; i < curSize_list; i++) {
            list[i] = list[i + 1];
            levels[i] = levels[i + 1];
        }
        curSize_list--;
    }
    return ret;
}

// HybridQueue_HQueue end
//////////////////////////////////////////////////////////////////////////////////////////////////////

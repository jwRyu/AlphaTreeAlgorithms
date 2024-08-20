#include <HierarHeapQueue_cache.h>

#include <allocator.h>
#include <cmath>
#include <cstring>
#include <defines.h>

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
    qsizes = (ImgIdx *)Malloc((size_t)numlevels_in * sizeof(ImgIdx));
    std::memcpy(qsizes, dhist, (size_t)numlevels_in * sizeof(ImgIdx));
    // qsizes = dhist;
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
}

template <class Pixel> void HierarHeapQueue_cache<Pixel>::push_1stitem(ImgIdx idx, Pixel alpha) {
    list[0].pidx = idx;
    list[0].alpha = alpha;
    curSize_list++;
}

template <class Pixel> void HierarHeapQueue_cache<Pixel>::end_pushes(_uint8 *isVisited) {
    if (emptytop)
        pop(isVisited);
}

template <class Pixel> void HierarHeapQueue_cache<Pixel>::push(ImgIdx idx, Pixel alpha) {
    if (emptytop && alpha < list[0].alpha) {
        emptytop = 0;
        list[0].pidx = idx;
        list[0].alpha = alpha;
        return;
    }

    bool push2list = (queue_minlev < curthr) ? alpha < hqueue[queue_minlev]->top_alpha()
                                             : (int)(a * log2(1.0 + (double)alpha)) < queue_minlev;

    if (push2list) {
        // spare room in the list
        if (curSize_list < maxSize_list) {
            int i;
            for (i = curSize_list; alpha < list[i].alpha; i--) {
                list[i + 1] = list[i];
            }
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
            curSize_list++;
        } else if (alpha < list[curSize_list].alpha) { // push to the full list
            push_queue(list[curSize_list].pidx, list[curSize_list].alpha);

            int i;
            for (i = curSize_list - 1; alpha < list[i].alpha; i--) {
                list[i + 1] = list[i];
            }
            list[i + 1].pidx = idx;
            list[i + 1].alpha = alpha;
        } else {
            push_queue(idx, alpha); // push to the queue
        }

    } else {
        push_queue(idx, alpha); // push to the queue
    }
}

template <class Pixel> void HierarHeapQueue_cache<Pixel>::push_queue(ImgIdx idx, Pixel alpha) {
    int level = (int)(a * log2(1 + (double)alpha));

    if (level < queue_minlev)
        queue_minlev = level;

    if (level < curthr) {
        hqueue[level]->push(idx, alpha);
    } else {
        ImgIdx cur = storage_cursize[level]++;
        storage[level][cur].pidx = idx;
        storage[level][cur].alpha = alpha;
    }
}

template <class Pixel> ImgIdx HierarHeapQueue_cache<Pixel>::pop(_uint8 *isVisited) {
    ImgIdx ret = top();
    if (curSize_list == 0) {
        while (!check_queue_level(isVisited))
            queue_minlev++;
        list[0].pidx = hqueue[queue_minlev]->top();
        list[0].alpha = hqueue[queue_minlev]->top_alpha();

        pop_queue(isVisited);
    } else {
        for (int i = 0; i < curSize_list; i++) {
            list[i] = list[i + 1];
        }
        curSize_list--;
    }
    return ret;
}

template <class Pixel> int HierarHeapQueue_cache<Pixel>::check_queue_level(_uint8 *isVisited) {
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

template <class Pixel> void HierarHeapQueue_cache<Pixel>::pop_queue(_uint8 *isVisited) {
    hqueue[queue_minlev]->pop();

    if (!hqueue[queue_minlev]->get_cursize()) {
        do {
            queue_minlev++;
        } while (queue_minlev < numlevels && !check_queue_level(isVisited));
    }
}

template class HierarHeapQueue_cache<_uint8>;
template class HierarHeapQueue_cache<_uint16>;
template class HierarHeapQueue_cache<_uint32>;
template class HierarHeapQueue_cache<_uint64>;

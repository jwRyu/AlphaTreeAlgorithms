#include <HHPQ.hpp>
#include <allocator.h>
#include <cmath>

template <class Pixel> ImgIdx HHPQ<Pixel>::alphaToLevel(const double &alpha) const {
    return (ImgIdx)(a * log2(1.0 + alpha));
}

template <class Pixel>
void HHPQ<Pixel>::initHQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, int connectivity,
                         double r) {
    maxSize = size;

    list = (QItem<Pixel> *)Malloc((listsize) * sizeof(QItem<Pixel>));
    maxSize_list = listsize - 1;
    curSize_list = -1;

    this->numlevels = numlevels_in;
    this->a = a_in;
    this->queue_minlev = numlevels;

    ImgIdx cumsum = 0;
    qsizes = dhist; // do not free dhist outside
    if (r >= 1) {
        thr_hqueue = curthr = numlevels;
        hqueue = (QuadHeapQueue<Pixel> **)Calloc(numlevels * sizeof(QuadHeapQueue<Pixel> *));
        for (int level = 0; level < thr_hqueue; level++)
            hqueue[level] = new QuadHeapQueue<Pixel>(qsizes[level]);
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

        hqueue = (QuadHeapQueue<Pixel> **)Calloc(numlevels * sizeof(QuadHeapQueue<Pixel> *));
        for (int level = 0; level < thr_hqueue; level++)
            hqueue[level] = new QuadHeapQueue<Pixel>(qsizes[level]);

        storage = (QItem<Pixel> **)Calloc((numlevels - thr_hqueue) * sizeof(QItem<Pixel> *));
        storage -= thr_hqueue;
        for (int level = thr_hqueue; level < numlevels; level++)
            storage[level] = (QItem<Pixel> *)Malloc(qsizes[level] * sizeof(QItem<Pixel>));
    }
}

template <class Pixel> void QuadHeapQueue<Pixel>::print() {
    for (int i = 0; i < cursize; i++)
        arr[i + 1].print();
    printf("\n");
}

template <class Pixel> void HHPQ<Pixel>::print() {
    printf("---------- HHPQ<Pixel>::print START -------------\n");
    printf("Cache[%d / %d]: ", curSize_list + 1, maxSize_list + 1);
    size_t size = curSize_list;
    for (int i = -1; i < curSize_list; i++)
        list[i + 1].print();
    printf("\n");

    for (int level = 0; level < numlevels; level++) {
        if (level < curthr) {
            size += hqueue[level]->size();
            if (hqueue[level]->empty())
                continue;
            printf("Q: level[%d][%d / %d]: ", level, hqueue[level]->size(), hqueue[level]->sizeMax());
            hqueue[level]->print();
        } else {
            if (storage_cursize[level] == 0)
                continue;
            size += storage_cursize[level];
            printf("S: level[%d][%d / %d]: ", level, storage_cursize[level], qsizes[level]);
            for (int i = 0; i < storage_cursize[level]; i++)
                storage[level][i].print();
            printf("\n");
        }
    }

    printf("size = %d\n", (int)size);
    printf("---------- HHPQ<Pixel>::print END -------------\n");
    // std::getchar();
}

template <class Pixel>
HHPQ<Pixel>::HHPQ(ImgIdx *dhist, ImgIdx numlevels_in, ImgIdx size, double a_in, int listsize, ImgIdx connectivity,
                  double r) {
    initHQ(dhist, numlevels_in, size, a_in, listsize, (int)connectivity, r);
}

template <class Pixel> HHPQ<Pixel>::~HHPQ() {
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

template <class Pixel> void HHPQ<Pixel>::push_1stitem(ImgIdx idx) {
    list[0].index = idx;
    list[0].alpha = std::numeric_limits<Pixel>::max();
    curSize_list++;
}

template <class Pixel> void HHPQ<Pixel>::end_pushes(_uint8 *isVisited) {
    if (emptytop)
        pop(isVisited);
}

template <class Pixel> void HHPQ<Pixel>::push(const ImgIdx &idx, const Pixel &alpha) {
    const QItem<Pixel> newItem(idx, alpha);

    // printf("Pushing %d at %.2f\n", idx, (double)alpha);
    if (emptytop && newItem < list[0]) {
        emptytop = 0;
        list[0] = newItem;
        return;
    }

    bool pushToList =
        (queue_minlev < curthr) ? newItem < hqueue[queue_minlev]->top() : alphaToLevel(newItem.alpha) < queue_minlev;

    if (pushToList) {
        if (curSize_list < maxSize_list) // spare room in the list
        {
            int i = curSize_list++;
            for (; i >= 0 && alpha < list[i].alpha; i--)
                list[i + 1] = list[i];

            list[i + 1] = newItem;
        } else if (alpha < list[curSize_list].alpha) // push to the full list
        {
            push_queue(list[curSize_list]);
            int i;
            for (i = curSize_list - 1; i >= 0 && alpha < list[i].alpha; i--)
                list[i + 1] = list[i];

            list[i + 1] = newItem;
        } else {
            push_queue(newItem); // push to the queue
        }
    } else {
        push_queue(newItem); // push to the queue
    }
}

template <class Pixel> void HHPQ<Pixel>::push_queue(const QItem<Pixel> &item) {
    auto level = alphaToLevel(item.alpha);

    if (level < queue_minlev)
        queue_minlev = level;

    if (level < curthr) {
        hqueue[level]->push(item);
    } else {
        ImgIdx cur = storage_cursize[level]++;
        storage[level][cur] = item;
    }
}

template <class Pixel> ImgIdx HHPQ<Pixel>::pop(_uint8 *isVisited) {
    ImgIdx ret = top();
    if (curSize_list) {
        for (int i = 0; i < curSize_list; i++)
            list[i] = list[i + 1];
        curSize_list--;
    } else {
        while (!check_queue_level(isVisited))
            queue_minlev++;
        list[0] = hqueue[queue_minlev]->top();
        pop_queue(isVisited);
    }
    return ret;
}

template <class Pixel> int HHPQ<Pixel>::check_queue_level(_uint8 *isVisited) {
    if (queue_minlev < curthr)
        return hqueue[queue_minlev]->get_cursize();
    else {
        while (curthr < queue_minlev) {
            hqueue[curthr] = new QuadHeapQueue<Pixel>(qsizes[curthr]);
            Free(storage[curthr]);
            storage[curthr] = 0;
            curthr++;
        }
        curthr++;

        hqueue[queue_minlev] = new QuadHeapQueue<Pixel>(qsizes[queue_minlev]);

        QItem<Pixel> *store = storage[queue_minlev];
        const ImgIdx cur = storage_cursize[queue_minlev];
        QuadHeapQueue<Pixel> *pQ = hqueue[queue_minlev];
        for (ImgIdx p = 0; p < cur; p++) {
            if (!isVisited[store[p].index])
                pQ->push(store[p]);
        }
        Free(storage[queue_minlev]);
        storage[queue_minlev] = 0;
        return pQ->get_cursize();
    }
}

template <class Pixel> void HHPQ<Pixel>::pop_queue(_uint8 *isVisited) {
    hqueue[queue_minlev]->pop();
    if (!hqueue[queue_minlev]->get_cursize()) {
        do {
            queue_minlev++;
        } while (queue_minlev < numlevels && !check_queue_level(isVisited));
    }
}

template class HHPQ<_uint8>;
template class HHPQ<_uint16>;
template class HHPQ<_uint32>;
template class HHPQ<_uint64>;
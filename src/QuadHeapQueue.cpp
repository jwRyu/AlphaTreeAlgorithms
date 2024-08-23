#include <QuadHeapQueue.hpp>
#include <allocator.h>
#include <cfloat>

template <class Pixel> QuadHeapQueue<Pixel>::QuadHeapQueue(ImgIdx maxsize_in) : cursize(0), maxsize(maxsize_in) {
    arr = (QItem<Pixel> *)Malloc(sizeof(QItem<Pixel>) * (maxsize + 2));
    if ((Pixel)-1 > (Pixel)1)
        arr[1].alpha = (Pixel)-1;
    else if (sizeof(Pixel) == 8)
        arr[1].alpha = DBL_MAX;
    else
        arr[1].alpha = FLT_MAX;
}

template <class Pixel> QuadHeapQueue<Pixel>::~QuadHeapQueue() { Free(arr); }

template <class Pixel> ImgIdx QuadHeapQueue<Pixel>::pop() {
    const ImgIdx outval = arr[1].index;
    const QItem<Pixel> &lastItem = arr[cursize];
    cursize--;
    if (cursize == 0)
        return 0;
    ImgIdx current = 1;
    while (true) {
        const ImgIdx next0 = (current << 2) - 2;
        ImgIdx next = next0;
        if (next0 + 3 <= cursize) {
            if (arr[next0 + 1] < arr[next])
                next = next0 + 1;
            if (arr[next0 + 2] < arr[next])
                next = next0 + 2;
            if (arr[next0 + 3] < arr[next])
                next = next0 + 3;
        } else {
            if (next0 > cursize)
                break;
            if (next0 == cursize)
                goto MIN_NEXT_FOUND;
            if (arr[next0 + 1] < arr[next])
                next = next0 + 1;
            if (next0 + 1 == cursize)
                goto MIN_NEXT_FOUND;
            if (arr[next0 + 2] < arr[next])
                next = next0 + 2;
            if (next0 + 2 == cursize)
                goto MIN_NEXT_FOUND;
            if (arr[next0 + 3] < arr[next])
                next = next0 + 3;
        }

    MIN_NEXT_FOUND:
        if (lastItem < arr[next])
            break;

        arr[current] = arr[next];
        current = next;
    }
    arr[current] = lastItem;
    return outval;
}

template <class Pixel> void QuadHeapQueue<Pixel>::push(QItem<Pixel> item) {
    cursize++;
    ImgIdx current = cursize;
    ImgIdx next = (current + 2) >> 2;
    while (next && (arr[next].alpha > item.alpha)) {
        arr[current] = arr[next];
        current = next;
        next = (next + 2) >> 2;
    }
    arr[current] = item;
}

template class QuadHeapQueue<_uint8>;
template class QuadHeapQueue<_uint16>;
template class QuadHeapQueue<_uint32>;
template class QuadHeapQueue<_uint64>;
template class QuadHeapQueue<float>;
template class QuadHeapQueue<double>;
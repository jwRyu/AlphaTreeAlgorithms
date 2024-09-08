#include <QuadHeapQueue.hpp>
#include <allocator.hpp>
QuadHeapQueue::QuadHeapQueue(ImgIdx maxsize_in) : cursize(0), maxsize(maxsize_in) {
    arr = (QItem *)Malloc(sizeof(QItem) * (maxsize + 2));
    arr[1].alpha = DBL_MAX;
}

void QuadHeapQueue::print() {
    for (int i = 0; i < cursize; i++)
        arr[i + 1].print();
    printf("\n");
}

QuadHeapQueue::~QuadHeapQueue() { Free(arr); }

ImgIdx QuadHeapQueue::pop() {
    const ImgIdx outval = arr[1].index;
    const QItem &lastItem = arr[cursize];
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

void QuadHeapQueue::push(QItem item) {
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
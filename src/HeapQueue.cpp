#include <HeapQueue.hpp>
#include <allocator.hpp>

#if PROFILE
#include <fstream>
#endif
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HeapQueue

template <class Pixel>
HeapQueue<Pixel>::HeapQueue(ImgIdx maxsize_in) : maxsize(maxsize_in), cursize(0), qtime(0), timing(0) {
    arr = new HQentry<Pixel>[maxsize + 1];
    pushlistidx = 0;
    flags[0] = flags[1] = flags[2] = flags[3] = flags[4] = flags[5] = flags[6] = 1;
}

template <class Pixel> HeapQueue<Pixel>::~HeapQueue() {
#if PROFILE
    std::ofstream outFile("MemmoveHEAP.txt", std::ios::app);
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
    delete[] arr;
}

template <class Pixel> int HeapQueue<Pixel>::minidx_pushlist() {
    int minidx = -1;
    for (int i = 0; i < pushlistidx; i++) {
        if (flags[i] && (minidx < 0 || pushlist[i].alpha < pushlist[minidx].alpha))
            minidx = i;
    }
    flags[minidx] = 0;
    return minidx;
}

template <class Pixel> void HeapQueue<Pixel>::find_minlev() {
    if (!pushlistidx) {
        pop();
        return;
    }
    int minidx = minidx_pushlist();
    if (pushlist[minidx].alpha <= arr[1].alpha) {
        arr[1].alpha = pushlist[minidx].alpha;
        arr[1].pidx = pushlist[minidx].pidx;
    } else {
        pop();
        push_run(pushlist[minidx].pidx, pushlist[minidx].alpha);
    }

    for (int i = 1; i < pushlistidx; i++) {
        minidx = minidx_pushlist();
        push_run(pushlist[minidx].pidx, pushlist[minidx].alpha);
    }
    pushlistidx = 0;
    flags[0] = flags[1] = flags[2] = flags[3] = flags[4] = flags[5] = flags[6] = 1;
}

template <class Pixel> ImgIdx HeapQueue<Pixel>::pop() {
#if PROFILE
    printf("pop %d, size %d\n", arr[1].pidx, cursize);
    num_memmove_pop_i = 0;
#endif
    ImgIdx outval = arr[1].pidx;
    ImgIdx current = 1, next, next0, next1, curidx;
    Pixel curalpha;
    curidx = arr[cursize].pidx;
    curalpha = arr[cursize].alpha;
    cursize--;

    while (1) {
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
#if PROFILE
        num_memmove_pop_i++;
#endif
    }
    arr[current].alpha = curalpha;
    arr[current].pidx = curidx;
#if PROFILE
    num_memmove_pop_i++;
#endif

#if PROFILE
    num_memmove_pop.push_back(num_memmove_pop_i);
    num_items_pop.push_back(cursize);
#endif
    return outval;
}

template <class Pixel> void HeapQueue<Pixel>::push(ImgIdx pidx, Pixel alpha) {
    pushlist[pushlistidx].pidx = pidx;
    pushlist[pushlistidx++].alpha = alpha;
}

template <class Pixel> void HeapQueue<Pixel>::push_run(ImgIdx pidx, Pixel alpha) {
#if PROFILE
    // printf("push %d, size %d\n", pidx, cursize);
    num_memmove_push_i = 0;
#endif
    ImgIdx current, next;

    cursize++;
    current = cursize;
    next = current >> 1;
    while (next && (arr[next].alpha > alpha)) {
#if PROFILE
        num_memmove_push_i++;
#endif
        arr[current] = arr[next];
        current = next;
        next = next >> 1;
    }

#if PROFILE
    num_memmove_push_i++;
#endif
    arr[current].pidx = pidx;
    arr[current].alpha = alpha;
#if PROFILE
    num_memmove_push.push_back(num_memmove_push_i);
    num_items_push.push_back(cursize);
#endif
}

template class HeapQueue<uint8_t>;
template class HeapQueue<uint16_t>;
template class HeapQueue<uint32_t>;
template class HeapQueue<uint64_t>;
template class HeapQueue<float>;
template class HeapQueue<double>;

// HeapQueue end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HeapQueue_naive

template <class Pixel>
HeapQueue_naive<Pixel>::HeapQueue_naive(ImgIdx maxsize_in) : cursize(0), maxsize(maxsize_in), qtime(0) {
    arr = new HQentry<Pixel>[maxsize + 1];
    if ((Pixel)-1 > (Pixel)1)
        arr[1].alpha = (Pixel)-1;
    else if (sizeof(Pixel) == 8)
        arr[1].alpha = DBL_MAX;
    else
        arr[1].alpha = FLT_MAX;
}

template <class Pixel> HeapQueue_naive<Pixel>::~HeapQueue_naive() {
#if PROFILE
    std::ofstream outFile("MemmoveHEAP.txt", std::ios::app);
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
    delete[] arr;
}

template <class Pixel> ImgIdx HeapQueue_naive<Pixel>::pop() {
#if PROFILE
    // printf("pop %d, size %d\n", arr[1].pidx, cursize);
    num_memmove_pop_i = 0;
#endif
    ImgIdx outval = arr[1].pidx;
    ImgIdx current = 1, next, next0, next1, curidx;
    Pixel curalpha;
    curidx = arr[cursize].pidx;
    curalpha = arr[cursize].alpha;
    cursize--;

    while (1) {
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
#if PROFILE
        num_memmove_pop_i++;
#endif
    }
    arr[current].alpha = curalpha;
    arr[current].pidx = curidx;
#if PROFILE
    num_memmove_pop_i++;
#endif
#if PROFILE
    num_memmove_pop.push_back(num_memmove_pop_i);
    num_items_pop.push_back(cursize);
#endif
    return outval;
}

template <class Pixel> void HeapQueue_naive<Pixel>::push(ImgIdx pidx, Pixel alpha) {
#if PROFILE
    // printf("push %d, size %d\n", pidx, cursize);
    num_memmove_push_i = 0;
#endif
    ImgIdx current, next;
    cursize++;
    current = cursize;
    next = current >> 1;
    while (next && (arr[next].alpha > alpha)) {
        arr[current] = arr[next];
        current = next;
        next = next >> 1;
#if PROFILE
        num_memmove_push_i++;
#endif
    }

#if PROFILE
    num_memmove_push_i++;
#endif
    arr[current].pidx = pidx;
    arr[current].alpha = alpha;
#if PROFILE
    num_memmove_push.push_back(num_memmove_push_i);
    num_items_push.push_back(cursize);
#endif
}

template class HeapQueue_naive<uint8_t>;
template class HeapQueue_naive<uint16_t>;
template class HeapQueue_naive<uint32_t>;
template class HeapQueue_naive<uint64_t>;
template class HeapQueue_naive<float>;
template class HeapQueue_naive<double>;

// HeapQueue_naive end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HeapQueue_naive_quad

template <class Pixel>
HeapQueue_naive_quad<Pixel>::HeapQueue_naive_quad(ImgIdx maxsize_in) : cursize(0), maxsize(maxsize_in), qtime(0) {
    arr = (HQentry<Pixel> *)Malloc(sizeof(HQentry<Pixel>) * (maxsize + 2));
    if ((Pixel)-1 > (Pixel)1)
        arr[1].alpha = (Pixel)-1;
    else if (sizeof(Pixel) == 8)
        arr[1].alpha = DBL_MAX;
    else
        arr[1].alpha = FLT_MAX;
}

template <class Pixel> HeapQueue_naive_quad<Pixel>::~HeapQueue_naive_quad() { Free(arr); }

template <class Pixel>
#if PROFILE
uint64_t HeapQueue_naive_quad<Pixel>::pop()
#else
ImgIdx HeapQueue_naive_quad<Pixel>::pop()
#endif
{
#if PROFILE
    uint64_t nummove = 0;
#else
    ImgIdx outval = arr[1].pidx;
#endif
    ImgIdx current = 1, next, next0, curidx;
    Pixel curalpha;
    curidx = arr[cursize].pidx;
    curalpha = arr[cursize].alpha;
    cursize--;
    if (cursize == 0)
        return 0;
    while (1) {
        next0 = (current << 2) - 2;
        next = next0;
        if (next0 + 3 <= cursize) {
            if (arr[next0 + 1].alpha < arr[next].alpha)
                next = next0 + 1;
            if (arr[next0 + 2].alpha < arr[next].alpha)
                next = next0 + 2;
            if (arr[next0 + 3].alpha < arr[next].alpha)
                next = next0 + 3;
        } else {
            if (next0 > cursize)
                break;
            if (next0 == cursize)
                goto MIN_NEXT_FOUND;
            if (arr[next0 + 1].alpha < arr[next].alpha)
                next = next0 + 1;
            if (next0 + 1 == cursize)
                goto MIN_NEXT_FOUND;
            if (arr[next0 + 2].alpha < arr[next].alpha)
                next = next0 + 2;
            if (next0 + 2 == cursize)
                goto MIN_NEXT_FOUND;
            if (arr[next0 + 3].alpha < arr[next].alpha)
                next = next0 + 3;
        }

    MIN_NEXT_FOUND:
        if (curalpha < arr[next].alpha)
            break;

        arr[current] = arr[next];
#if PROFILE
        nummove++;
#endif
        current = next;
    }
    arr[current].alpha = curalpha;
    arr[current].pidx = curidx;
#if PROFILE
    nummove++;
    return nummove;
#else
    return outval;
#endif
}

template <class Pixel>
#if PROFILE
uint64_t HeapQueue_naive_quad<Pixel>::push(ImgIdx pidx, Pixel alpha)
#else
void HeapQueue_naive_quad<Pixel>::push(ImgIdx pidx, Pixel alpha)
#endif
{
#if PROFILE
    uint64_t nummove = 0;
#endif
    ImgIdx current, next;
    cursize++;
    current = cursize;

    next = (current + 2) >> 2;
    while (next && (arr[next].alpha > alpha)) {
        arr[current] = arr[next];
        current = next;
        next = (next + 2) >> 2;
#if PROFILE
        nummove++;
#endif
    }

    arr[current].pidx = pidx;
    arr[current].alpha = alpha;
#if PROFILE
    nummove++;
    return nummove;
#endif
}

template class HeapQueue_naive_quad<uint8_t>;
template class HeapQueue_naive_quad<uint16_t>;
template class HeapQueue_naive_quad<uint32_t>;
template class HeapQueue_naive_quad<uint64_t>;
template class HeapQueue_naive_quad<float>;
template class HeapQueue_naive_quad<double>;

// HeapQueue_naive_quad end
//////////////////////////////////////////////////////////////////////////////////////////////////////
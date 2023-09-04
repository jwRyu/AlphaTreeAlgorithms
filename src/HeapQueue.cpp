#include "HeapQueue.h"
#include "allocator.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////
// HeapQueue

template<class Pixel>
HeapQueue<Pixel>::HeapQueue(Imgidx maxsize_in) : maxsize(maxsize_in), cursize(0), qtime(0), timing(0)
{
    arr = new HQentry<Pixel>[maxsize + 1];
    pushlistidx = 0;
    flags[0] = flags[1] = flags[2] = flags[3] = flags[4] = flags[5] = flags[6] = 1;
}

template<class Pixel>
HeapQueue<Pixel>::~HeapQueue()
{
    delete[] arr;
}

template<class Pixel>
int HeapQueue<Pixel>::minidx_pushlist()
{
    int minidx = -1;
    for (int i = 0; i < pushlistidx; i++)
    {
        if (flags[i] && (minidx < 0 || pushlist[i].alpha < pushlist[minidx].alpha))
            minidx = i;
    }
    flags[minidx] = 0;
    return minidx;
}

template<class Pixel>
void HeapQueue<Pixel>::find_minlev()
{
    if(!pushlistidx)
    {
        pop();
        return;
    }
    int minidx = minidx_pushlist();
    if (pushlist[minidx].alpha <= arr[1].alpha)
    {
        arr[1].alpha = pushlist[minidx].alpha;
        arr[1].pidx = pushlist[minidx].pidx;
    }
    else
    {
        pop();
        push_run(pushlist[minidx].pidx, pushlist[minidx].alpha);
    }

    for (int i = 1; i < pushlistidx; i++)
    {
        minidx = minidx_pushlist();
        push_run(pushlist[minidx].pidx, pushlist[minidx].alpha);
    }
    pushlistidx = 0;
    flags[0] = flags[1] = flags[2] = flags[3] = flags[4] = flags[5] = flags[6] = 1;
}

template<class Pixel>
Imgidx HeapQueue<Pixel>::pop()
{
    Imgidx outval = arr[1].pidx;
    Imgidx current = 1, next, next0, next1, curidx;
    Pixel curalpha;
    curidx = arr[cursize].pidx;
    curalpha = arr[cursize].alpha;
    cursize--;

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

    return outval;
}

template<class Pixel>
void HeapQueue<Pixel>::push(Imgidx pidx, Pixel alpha)
{
    pushlist[pushlistidx].pidx = pidx;
    pushlist[pushlistidx++].alpha = alpha;
}

template<class Pixel>
void HeapQueue<Pixel>::push_run(Imgidx pidx, Pixel alpha)
{
    Imgidx current, next;

    cursize++;
    current = cursize;
    next = current >> 1;
    while (next && (arr[next].alpha > alpha))
    {
        arr[current] = arr[next];
        current = next;
        next = next >> 1;
    }

    arr[current].pidx = pidx;
    arr[current].alpha = alpha;
}

template class HeapQueue<_uint8>;
template class HeapQueue<_uint16>;
template class HeapQueue<_uint32>;
template class HeapQueue<_uint64>;
template class HeapQueue<float>;
template class HeapQueue<double>;

// HeapQueue end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HeapQueue_naive

template<class Pixel>
HeapQueue_naive<Pixel>::HeapQueue_naive(Imgidx maxsize_in) : cursize(0), maxsize(maxsize_in), qtime(0)
{
    arr = new HQentry<Pixel>[maxsize + 1];
    if((Pixel)-1 > (Pixel)1)
        arr[1].alpha = (Pixel)-1;
    else if(sizeof(Pixel) == 8)
        arr[1].alpha = DBL_MAX;
    else
        arr[1].alpha = FLT_MAX;
}

template<class Pixel>
HeapQueue_naive<Pixel>::~HeapQueue_naive()
{
    delete[] arr;
}

template<class Pixel>
Imgidx HeapQueue_naive<Pixel>::pop()
{
    Imgidx outval = arr[1].pidx;
    Imgidx current = 1, next, next0, next1, curidx;
    Pixel curalpha;
    curidx = arr[cursize].pidx;
    curalpha = arr[cursize].alpha;
    cursize--;

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

    return outval;
}

template<class Pixel>
void HeapQueue_naive<Pixel>::push(Imgidx pidx, Pixel alpha)
{
    Imgidx current, next;
    cursize++;
    current = cursize;
    next = current >> 1;
    while (next && (arr[next].alpha > alpha))
    {
        arr[current] = arr[next];
        current = next;
        next = next >> 1;
    }

    arr[current].pidx = pidx;
    arr[current].alpha = alpha;
}

template class HeapQueue_naive<_uint8>;
template class HeapQueue_naive<_uint16>;
template class HeapQueue_naive<_uint32>;
template class HeapQueue_naive<_uint64>;
template class HeapQueue_naive<float>;
template class HeapQueue_naive<double>;

// HeapQueue_naive end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HeapQueue_naive_quad

template<class Pixel>
HeapQueue_naive_quad<Pixel>::HeapQueue_naive_quad(Imgidx maxsize_in) : cursize(0), maxsize(maxsize_in), qtime(0)
{
    arr = (HQentry<Pixel>*)Malloc(sizeof(HQentry<Pixel>) * (maxsize + 2));
    if((Pixel)-1 > (Pixel)1)
        arr[1].alpha = (Pixel)-1;
    else if(sizeof(Pixel) == 8)
        arr[1].alpha = DBL_MAX;
    else
        arr[1].alpha = FLT_MAX;
}

template<class Pixel>
HeapQueue_naive_quad<Pixel>::~HeapQueue_naive_quad()
{
    Free(arr);
}

template<class Pixel>
Imgidx HeapQueue_naive_quad<Pixel>::pop()
{
    Imgidx outval = arr[1].pidx;
    Imgidx current = 1, next, next0, curidx;
    Pixel curalpha;
    curidx = arr[cursize].pidx;
    curalpha = arr[cursize].alpha;
    cursize--;
    if(cursize == 0)
        return 0;
    while (1)
    {
        next0 = (current << 2) - 2;
        next = next0;
        if(next0 + 3 <= cursize)
        {
            if (arr[next0 + 1].alpha < arr[next].alpha)	next = next0 + 1;
            if (arr[next0 + 2].alpha < arr[next].alpha)	next = next0 + 2;
            if (arr[next0 + 3].alpha < arr[next].alpha)	next = next0 + 3;
        }
        else
        {
            if (next0 > cursize)	break;
            if (next0 == cursize)			    	goto MIN_NEXT_FOUND;
            if (arr[next0 + 1].alpha < arr[next].alpha)	next = next0 + 1;
            if (next0 + 1 == cursize)				goto MIN_NEXT_FOUND;
            if (arr[next0 + 2].alpha < arr[next].alpha)	next = next0 + 2;
            if (next0 + 2 == cursize)				goto MIN_NEXT_FOUND;
            if (arr[next0 + 3].alpha < arr[next].alpha)	next = next0 + 3;
        }

MIN_NEXT_FOUND:
        if (curalpha < arr[next].alpha)
            break;

        arr[current] = arr[next];
        current = next;
    }
    arr[current].alpha = curalpha;
    arr[current].pidx = curidx;
    return outval;
}

template<class Pixel>
void HeapQueue_naive_quad<Pixel>::push(Imgidx pidx, Pixel alpha)
{
    Imgidx current, next;
    cursize++;
    current = cursize;

    next = (current + 2) >> 2;
    while (next && (arr[next].alpha > alpha))
    {
        arr[current] = arr[next];
        current = next;
        next = (next + 2) >> 2;
    }

    arr[current].pidx = pidx;
    arr[current].alpha = alpha;
}

template class HeapQueue_naive_quad<_uint8>;
template class HeapQueue_naive_quad<_uint16>;
template class HeapQueue_naive_quad<_uint32>;
template class HeapQueue_naive_quad<_uint64>;
template class HeapQueue_naive_quad<float>;
template class HeapQueue_naive_quad<double>;

// HeapQueue_naive_quad end
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
// HeapQueue_rank

HeapQueue_rank::HeapQueue_rank(Imgidx maxsize_in) : cursize(0), maxsize(maxsize_in)
{
    arr = new Imgidx[maxsize];
    arr--;
}

HeapQueue_rank::~HeapQueue_rank()
{
    delete[] (arr + 1);
}

Imgidx HeapQueue_rank::pop()
{
    Imgidx outval;
    Imgidx current = 1, next, next0, next1, curidx;
    outval = arr[current];
    curidx = arr[cursize--];

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
        current = next;
    }
    arr[current] = curidx;

    if (cursize)
        min_level = arr[1];
    else
        min_level = (Imgidx)-1;

    return outval;
}

void HeapQueue_rank::push(Imgidx pidx)
{
    Imgidx current, next;
    cursize++;
    current = cursize;
    next = current >> 1;
    while (next && (arr[next] > pidx))
    {
        arr[current] = arr[next];
        current = next;
        next = next >> 1;
    }

    arr[current] = pidx;

    if (current == 1)
    {
        min_level = pidx;
    }
}

// HeapQueue_rank end
//////////////////////////////////////////////////////////////////////////////////////////////////////

#include "HierarQueue.h"

void HierarQueue::print() {
    for (int32_t i = min_level; i < numlevel; i++) {
        if (cur[i] == bottom[i])
            continue;
        printf("Level %d(%d-%d): ", i, (int)bottom[i], (int)cur[i] - 1);
        for (ImgIdx j = cur[i] - 1; j != bottom[i]; j--)
            printf("%d ", (int)queue[j]);
        printf("%d ", (int)queue[bottom[i]]);
        printf("\n");
    }
}

HierarQueue::HierarQueue(uint64_t qsize_in, int32_t numlevels) {
    queue = (ImgIdx *)Malloc((size_t)qsize_in * sizeof(ImgIdx));
    bottom = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));
    cur = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));

    this->numlevel = numlevels;
    max_level = 0;

    qsize = qsize_in;
    min_level = numlevels - 1;

    bottom[numlevels] = 0;
    cur[numlevels] = 1;
}

HierarQueue::HierarQueue(uint64_t qsize_in) {
    int32_t numlevels = (int32_t)1 << 20;
    queue = (ImgIdx *)Malloc((size_t)qsize_in * sizeof(ImgIdx));
    bottom = (ImgIdx *)Malloc(((size_t)numlevels + 1) * sizeof(ImgIdx));
    cur = (ImgIdx *)Malloc(((size_t)numlevels + 1) * sizeof(ImgIdx));

    this->numlevel = numlevels;
    max_level = 0;

    qsize = qsize_in;
    min_level = numlevels - 1;

    bottom[numlevels] = 0;
    cur[numlevels] = 1;
}

void HierarQueue::reset_queue() {
    for (int32_t i = 0; i < numlevel; i++)
        cur[i] = bottom[i];
    min_level = numlevel;
}

ImgIdx HierarQueue::set_queue(ImgIdx *dhist) {
    ImgIdx sum_hist = 0;
    for (int32_t i = 0; i < numlevel; i++) {
        bottom[i] = cur[i] = sum_hist;

        if (dhist[i]) {
            max_level = i;
            sum_hist += dhist[i];
        }
    }
    return sum_hist;
}

ImgIdx HierarQueue::set_queue(ImgIdx *dhist, int32_t maxpix) {
    ImgIdx sum_hist = 0;
    int32_t numlevels = (int32_t)1 << 20;

    for (int i = 0; i < numlevels; i++)
        bottom[i] = cur[i] = 0;
    for (int i = 0; i < qsize; i++)
        queue[i] = 0;

    max_level = 0;
    for (int32_t i = 0; i <= maxpix; i++) {
        bottom[i] = cur[i] = sum_hist;
        if (dhist[i]) {
            max_level = i;
            sum_hist += dhist[i];
        }
    }

    min_level = maxpix + 1;
    bottom[maxpix + 1] = 0;
    cur[maxpix + 1] = 1;

    return sum_hist;
}

HierarQueue::HierarQueue(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels) {
    queue = (ImgIdx *)Malloc((size_t)qsize_in * sizeof(ImgIdx));
    bottom = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));
    cur = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));

    this->numlevel = numlevels;
    max_level = 0;

    qsize = qsize_in;
    min_level = numlevels - 1;

    ImgIdx sum_hist = 0;
    for (int32_t i = 0; i < numlevels; i++) {
        bottom[i] = cur[i] = sum_hist;
        if (dhist[i]) {
            max_level = i;
            sum_hist += dhist[i];
        }
    }

    bottom[numlevels] = 0;
    cur[numlevels] = 1;
}

HierarQueue::HierarQueue(ImgIdx *dhist, int32_t numlevels) {
    uint64_t dsum = 0;
    max_level = 0;
    for (int i = 0; i < numlevels; i++) {
        if (dhist[i]) {
            dsum += dhist[i];
            max_level = i;
        }
    }

    queue = (ImgIdx *)Malloc((size_t)(dsum + 1) * sizeof(ImgIdx));
    bottom = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));
    cur = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));

    this->numlevel = numlevels;

    qsize = (dsum + 1);
    min_level = numlevels - 1;

    ImgIdx sum_hist = 0;
    for (int32_t i = 0; i < numlevels; i++) {
        bottom[i] = cur[i] = sum_hist;
        sum_hist += dhist[i];
    }
    bottom[numlevels] = 0;
    cur[numlevels] = 1;
}

HierarQueue::HierarQueue(int32_t numlevels, ImgIdx binsize) {
    qsize = numlevels * binsize;
    // tmp
    queue = (ImgIdx *)Malloc((size_t)qsize * sizeof(ImgIdx));
    bottom = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));
    cur = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));
    min_level = numlevels - 1;
    this->numlevel = numlevels;

    ImgIdx sum_hist = 0;
    for (int32_t i = 0; i < numlevels; i++) {
        bottom[i] = cur[i] = sum_hist;
        sum_hist += binsize;
    }
    bottom[numlevels] = 0;
    cur[numlevels] = 1;
}

HierarQueue::~HierarQueue() {
    Free(queue);
    Free(bottom);
    Free(cur);
}

int8_t HierarQueue::push(ImgIdx pidx, int64_t level) {
    queue[cur[level]++] = pidx;

    if (level < min_level) {
        min_level = level;
        return 1;
    } else
        return 0;
}

void HierarQueue::find_minlev() {
    while (bottom[min_level] == cur[min_level])
        min_level++;
}

HQueue_l1idx::HQueue_l1idx(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels) {
    // tmp
    queue = (ImgIdx *)Malloc((size_t)qsize_in * sizeof(ImgIdx));
    bottom = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));
    cur = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));
    seekersize = (numlevels + 1 + 63) >> 6;
    seeker = (uint64_t *)Malloc((size_t)(seekersize) * sizeof(uint64_t));

    qsize = qsize_in;
    min_level = numlevels - 1;

    ImgIdx sum_hist = 0;
    for (int32_t i = 0; i < numlevels; i++) {
        bottom[i] = cur[i] = sum_hist;
        sum_hist += dhist[i];
    }
    for (int64_t i = 0; i < seekersize; i++)
        seeker[i] = 0;
    seeker[numlevels >> 6] |= (uint64_t)1 << (numlevels & 63);
    bottom[numlevels] = 0;
    cur[numlevels] = 1;
}

HQueue_l1idx::~HQueue_l1idx() {
    Free(queue);
    Free(bottom);
    Free(cur);
    Free(seeker);
}

int HQueue_l1idx::push(ImgIdx pidx, int32_t level) {
    int64_t qidx = cur[level]++;
    queue[qidx] = pidx;
    seeker[level >> 6] |= (uint64_t)1 << (level & 63);
    if (level <= min_level) {
        min_level = level;
        return 1;
    }
    return 0;
}

ImgIdx HQueue_l1idx::pop() {
    ImgIdx popidx = --cur[min_level];

    if (bottom[min_level] == cur[min_level])
        seeker[min_level >> 6] &= ~((uint64_t)1 << (min_level & 63));
    return queue[popidx];
}

void HQueue_l1idx::find_minlev() {
    ImgIdx qidx, widx;
    uint64_t w;

    for (qidx = min_level >> 6; !seeker[qidx]; qidx++)
        ;

    w = seeker[qidx];

    if (w & 0xffffffff)
        widx = 0;
    else {
        widx = 32;
        w >>= 32;
    }

    while (!(w & (uint64_t)1)) {
        w >>= 1;
        widx++;
    }

    min_level = ((qidx << 6) + widx);
}

HQueue_l2idx::HQueue_l2idx(uint64_t qsize_in, ImgIdx *dhist, int32_t numlevels) {
    int64_t seekersize, seeker2size;
    // tmp
    queue = (ImgIdx *)Malloc((size_t)qsize_in * sizeof(ImgIdx));
    bottom = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));
    cur = (ImgIdx *)Malloc((size_t)(numlevels + 1) * sizeof(ImgIdx));

    seekersize = (numlevels + 1 + 63) >> 6;
    seeker2size = (seekersize + 63) >> 6;

    seeker = (uint64_t *)Malloc((size_t)(seekersize) * sizeof(uint64_t));
    seeker2 = (uint64_t *)Malloc((size_t)(seeker2size) * sizeof(uint64_t));

    qsize = qsize_in;
    min_level = numlevels - 1;

    ImgIdx sum_hist = 0;
    for (int32_t i = 0; i < numlevels; i++) {
        bottom[i] = cur[i] = sum_hist;
        sum_hist += dhist[i];
    }
    for (int64_t i = 0; i < seekersize; i++)
        seeker[i] = 0;
    seeker[numlevels >> 6] |= (uint64_t)1 << (numlevels & 63);
    for (int64_t i = 0; i < seeker2size; i++)
        seeker2[i] = 0;
    seeker2[numlevels >> 12] |= (uint64_t)1 << ((numlevels >> 6) & 63);
    bottom[numlevels] = 0;
    cur[numlevels] = 1;
}

HQueue_l2idx::~HQueue_l2idx() {
    Free(queue);
    Free(bottom);
    Free(cur);
    Free(seeker);
    Free(seeker2);
}

void HQueue_l2idx::push(ImgIdx pidx, int64_t level) {
    int64_t qidx = cur[level]++;
    queue[qidx] = pidx;
    seeker[level >> 6] |= (uint64_t)1 << (level & 63);
    seeker2[level >> 12] |= (uint64_t)1 << ((level >> 6) & 63);
    if (level < min_level) {
        min_level = level;
    }
}

ImgIdx HQueue_l2idx::pop() {
    ImgIdx popidx = --cur[min_level];

    if (bottom[min_level] == cur[min_level]) {
        seeker[min_level >> 6] &= ~((uint64_t)1 << (min_level & 63));
        if (!seeker[min_level >> 6])
            seeker2[min_level >> 12] &= ~((uint64_t)1 << ((min_level >> 6) & 63));
    }
    return queue[popidx];
}

void HQueue_l2idx::find_minlev() {
    ImgIdx qidx, widx;
    uint64_t w;

    for (qidx = min_level >> 12; !seeker2[qidx]; qidx++)
        ;

    w = seeker2[qidx];
    if (w & 0xffffffff)
        widx = 0;
    else {
        widx = 32;
        w >>= 32;
    }

    while (!(w & (uint64_t)1)) {
        w >>= 1;
        widx++;
    }

    qidx = ((qidx << 6) + widx);

    w = seeker[qidx];
    if (w & 0xffffffff)
        widx = 0;
    else {
        widx = 32;
        w >>= 32;
    }

    while (!(w & (uint64_t)1)) {
        w >>= 1;
        widx++;
    }

    min_level = ((qidx << 6) + widx);
}

HQueue_l1idx_rank::HQueue_l1idx_rank(int64_t qsize_in) {
    qsize = (qsize_in + (1 << 12)) >> 12;
    queue = (hqueue_word *)Calloc((size_t)qsize * sizeof(hqueue_word));

    min_level = qsize_in;
}

HQueue_l1idx_rank::~HQueue_l1idx_rank() { Free(queue); }

void HQueue_l1idx_rank::push(ImgIdx pidx) {
    int64_t qidx = pidx >> 12;
    int64_t widx = (pidx >> 6) & 63;
    int64_t bitpos = pidx & 63;

    queue[qidx].qword[widx] |= (int64_t)1 << (bitpos);
    queue[qidx].seeker |= (int64_t)1 << (widx);
    min_level = min_level < pidx ? min_level : pidx;
}

void HQueue_l1idx_rank::pop() {
    int64_t qidx = min_level >> 12;
    int64_t widx = (min_level >> 6) & 63;
    int64_t bitpos = min_level & 63;
    int64_t w, skr;

    queue[qidx].qword[widx] &= ~((int64_t)1 << (bitpos));
    if (!queue[qidx].qword[widx])
        queue[qidx].seeker &= ~((int64_t)1 << (widx));

    for (qidx = min_level >> 12; !queue[qidx].seeker; qidx++)
        ;

    skr = queue[qidx].seeker;
    widx = (skr & 0xffffffff) ? 0 : 32;
    widx += ((skr >> widx) & 0xffff) ? 0 : 16;
    widx += ((skr >> widx) & 0xff) ? 0 : 8;
    widx += ((skr >> widx) & 0xf) ? 0 : 4;
    widx += ((skr >> widx) & 0x3) ? 0 : 2;
    widx += ((skr >> widx) & 0x1) ? 0 : 1;

    w = queue[qidx].qword[widx];
    bitpos = (w & 0xffffffff) ? 0 : 32;
    bitpos += ((w >> bitpos) & 0xffff) ? 0 : 16;
    bitpos += ((w >> bitpos) & 0xff) ? 0 : 8;
    bitpos += ((w >> bitpos) & 0xf) ? 0 : 4;
    bitpos += ((w >> bitpos) & 0x3) ? 0 : 2;
    bitpos += ((w >> bitpos) & 0x1) ? 0 : 1;

    min_level = (qidx << 12) + (widx << 6) + bitpos;
}

void HQueue_l1idx_rank::find_minlev() {
    int64_t qidx, widx, bitpos, w, skr;

    for (qidx = min_level >> 12; !queue[qidx].seeker; qidx++)
        ;

    skr = queue[qidx].seeker;
    widx = (skr & 0xffffffff) ? 0 : 32;
    widx += ((skr >> widx) & 0xffff) ? 0 : 16;
    widx += ((skr >> widx) & 0xff) ? 0 : 8;
    widx += ((skr >> widx) & 0xf) ? 0 : 4;
    widx += ((skr >> widx) & 0x3) ? 0 : 2;
    widx += ((skr >> widx) & 0x1) ? 0 : 1;

    w = queue[qidx].qword[widx];
    bitpos = (w & 0xffffffff) ? 0 : 32;
    bitpos += ((w >> bitpos) & 0xffff) ? 0 : 16;
    bitpos += ((w >> bitpos) & 0xff) ? 0 : 8;
    bitpos += ((w >> bitpos) & 0xf) ? 0 : 4;
    bitpos += ((w >> bitpos) & 0x3) ? 0 : 2;
    bitpos += ((w >> bitpos) & 0x1) ? 0 : 1;

    min_level = (qidx << 12) + (widx << 6) + bitpos;
}

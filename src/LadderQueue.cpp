
#include <cmath>
#include "assert.h"
#include "LadderQueue.hpp"

void LinkedEdgeList::add(LQNode* n) {
    nodeCount++;
    if (!head) {
        head = n;
        tail = n;
    } else {
        tail->next = n;
        n->prev = tail;
        tail = n;
    }
    assert(tail->next != tail);
    assert(n->info);
}

void LinkedEdgeList::add(Event* evt) {
    LQNode* n = new LQNode();
    n->info = evt;
    add(n);
}


void LinkedEdgeList::addInOrder(Event* evt) {
    LQNode* n = new LQNode();
    n->info = evt;
    addInOrder(n);
}

void LinkedEdgeList::addInOrder(LQNode* n) {
    LQNode* cur = tail;
    while (cur != nullptr) {
        if (cur->info->alpha <= n->info->alpha) {
            
            n->next = cur->next;
            if (n->next != nullptr) {
                n->next->prev = n;
            } else {
                tail = n;
            }
            n->prev = cur;
            cur->next = n;
            nodeCount++;
        
            return;
        } else {
            cur = cur->prev;
        }
    }

    n->next = head;
    head = n;
    if (n->next == nullptr)
        tail = n;
    else
        n->next->prev = n;
    nodeCount++;

    assert(tail->next != tail);
}


void LinkedEdgeList::addFirst(Event* evt) {
    nodeCount++;

    if (!head) {
        head = new LQNode();
        head->info = evt;
        tail = head;
    } else {
        LQNode* n = head;
        head = new LQNode();
        head->info = evt;
        head->next = n;
        n->prev = head;
    }
    assert(tail->next != tail);
}

Event* LinkedEdgeList::removeDeleteNode() {
    assert(nodeCount > 0);
    Event* evt = head->info;

    LQNode* n = head;
    head = head->next;
    if (head != nullptr)
        head->prev = nullptr;
    else
        tail = nullptr;
    nodeCount--;
    delete n;
    return evt;
}

LQNode* LinkedEdgeList::removeNode() {
    LQNode* n = head;
    head = head->next;
    if (head != nullptr)
        head->prev = nullptr;
    else
        tail = nullptr;
    --nodeCount;
    n->next = nullptr;
    n->prev = nullptr;
    return n;
}

void LinkedEdgeList::sort() {
    if (nodeCount <= 1)
        return;
    
    std::vector<LQNode*> temp;
    // temp.resize(nodeCount);
    while (!isEmpty())
        temp.push_back(removeNode());

    while (temp.size()) {
        LQNode *n = temp.back();
        temp.pop_back();
        addInOrder(n);
    }
}

void LinkedEdgeList::print() {
    LQNode *walker = head;
    if (head == NULL) {
        printf("EMPTY");
        return;
    }
        
    while (1) {
        printf("%d(%f)", walker->info->idx, walker->info->alpha);
        if (walker == tail)
            break;
        printf("->");
        walker = walker->next;
    }
}

Rung::Rung(double w, double bw) : bucketWidth(bw) {
    // Initialize the Rung with the specified parameters
    // Implement the constructor
    assert(w > 0 && bw > 0 && w > bw);
    bucketCount = (int)(std::ceil(w / bw));
    bucket.resize(bucketCount);
    for (int bidx = 0;bidx < bucketCount;bidx++)
        bucket[bidx] = nullptr;
}

Rung::~Rung() {
    for (int bidx = 0;bidx < bucketCount;bidx++) {
        if (bucket[bidx])
            delete bucket[bidx];
    }
    bucket.clear();
}

int Rung::getRungIdx(double alpha) {
    assert(alpha >= rStart && alpha <= rStart + bucketWidth * (double)bucketCount);
    int bIdx = (int)((alpha - rStart) / bucketWidth);
    if (bIdx >= bucketCount)
        bIdx = bucketCount - 1;
    return bIdx;
}

void Rung::checkRungExist(int bidx) {
    if (bucket[bidx] == nullptr) {
        bucket[bidx] = new LinkedEdgeList;
    }
}


void Rung::addInOrder(Event* evt) {
    int bidx = getRungIdx(evt->alpha);
    checkRungExist(bidx);
    bucket[bidx]->addInOrder(evt);
}

void Rung::add(Event* evt) {
    int bidx = getRungIdx(evt->alpha);
    checkRungExist(bidx);
    bucket[bidx]->add(evt);
}

void Rung::add(LQNode* n) {
    int bidx = getRungIdx(n->info->alpha);
    checkRungExist(bidx);
    bucket[bidx]->add(n);
}

bool Rung::print(LinkedEdgeList* bottom) {
    bool ret = false;
    printf("  Range = %f-%f, rStart = %f, rCur = %f, bucketWidth = %f, bucketCount = %d, minBucket = %d, maxBucket = %d\n",
     rStart, rStart + (double)bucketCount * bucketWidth, rStart, rCur(), bucketWidth, bucketCount, minBucket(), maxBucket());
    int bidx = 0;
    for (auto b : bucket) {
        if (b == bottom) {
            printf("    (*BOTTOM) Bucket[%d]: ", bidx++);
            ret = true;
        }
        else
            printf("    Bucket[%d]: ", bidx++);
        if (b && b->getNodeCount()) b->print();
        printf("\n");
    }
    return ret;
}

int Rung::minBucket() {
    int mb;
    for (mb = 0;mb < (int)bucket.size() && bucket[mb] == nullptr;mb++)
        ;
    return mb;
}

int Rung::maxBucket() {
    int mb;
    for (mb = (int)(bucketCount - 1);mb >= 0 && bucket[mb] == nullptr;mb--)
        ;
    return mb;
}

void LadderQueue::enqueue(Imgidx idx, double alpha) {
    Event* evt = new Event(idx, alpha);
    // print();
    ++size;
    const double ts = evt->alpha;

    if (ts >= topStart()) { // insert in top
        topInsert++;
        if (top->getNodeCount() == 0) {
            MinTS = MaxTS = ts;
        } else {
            if (ts > MaxTS)
                MaxTS = ts;
            if (ts < MinTS)
                MinTS = ts;
        }
        top->add(evt);
        return;
    }

    // if (bottom != nullptr && !bottom->isEmpty() && ts < bottom->getFirst()->alpha) {
    //     bottom->addFirst(evt);
    //     bottomInsert++;
    //     return;
    // }

    for (auto& rung : rungs) {
        if (ts >= rung->rCur()) {
            // insert in a rung
            if (rung->bucket[rung->minBucket()] == bottom)
                rung->addInOrder(evt);
            else
                rung->add(evt);
            rungInsert++;
            return;
        }
    }

    // Shouldn't reach here?
    // assert(false);

    if (bottom == NULL || bottom->getNodeCount() == 0 
    || rungs.back()->bucket[rungs.back()->minBucket()] == bottom) {
        // printf("*** replacing bottom ***\n");
        bottom = new LinkedEdgeList;    
    }
    bottomInsert++;
    bottom->addInOrder(evt);
        
    if (0 && bottom->getNodeCount() >= THRES && rungs.size() > 0
             && bottom->getFirst()->alpha != bottom->getLast()->alpha) {
        // Spawn a Rung
        double w = bottom->getLast()->alpha - bottom->getFirst()->alpha;
        double bw = w / (double)THRES;
        Rung* r = new Rung(w, bw);
        r->rStart = bottom->getFirst()->alpha;

        while (!bottom->isEmpty()) 
            r->add(bottom->removeNode());
        r->add(evt);
        rungs.push_back(r);
        if (rungs.size() > rungused)
            rungused = rungs.size();
        rungInsert++;
        delete bottom;
        bottom = nullptr;
    }

    // Rung* lastRung = rungs.back();
    // for (bottomIndex = 0;lastRung->bucket[bottomIndex] != bottom &&
    //          bottomIndex != (int)lastRung->bucket.size();bottomIndex++)
    //     ;
    // assert(bottomIndex != (int)lastRung->bucket.size());
}

/**
 * It returns event with the highest priority: element is first searched in
 * Bottom (sorted), if latter is empty, then the element with the highest
 * priority is searched in the last rung. Once the element to be returned is
 * found in the bucket, number of events of the bucket is compared to THRES: -
 * if < THRES, after sorting the bucket is moved into the bottom returning the
 * most priority event; - otherwise, a new rung is generated transferring
 * elements in the considered bucket and so on. If also the Ladder is empty,
 * elements in Top are moved in a rung (the first one) and process is
 * iterated.If the structure is empty an error comes out.
 * 
 * @return Event with the highest priority
 * 
 */
Imgidx LadderQueue::dequeue() {
    // print();
    --size;
    assert(inspectRung());
    Event* evt = bottom->removeDeleteNode();
    Imgidx ret = evt->idx;
    std::printf("dequeue %d at %f\n", evt->idx, evt->alpha);
    if (evt->idx == 22) {
        int aa = 10;
        aa = aa * 10;
    }
    // print();
    delete evt;
    return ret;
}

void LadderQueue::print() {
    printf(" Queue size = %d\n", size);
    printf(" Top: ");
    top->print();
    printf("\n");
    // printf("BottomIndex = %d\n", bottomIndex);
    int ridx = 0;
    bool bottomPrint = false;
    for (auto r : rungs) {
        printf("  Rung[%d]:\n", ridx++);
        bottomPrint |= r->print(bottom);

    }
    printf("\n");
    if (bottomPrint)
        return;
    printf(" Bottom: ");
    if (bottom) bottom->print();
    printf("\n");
}

/// @brief Inspect the bottom Rung, reorganize/repopulate as necessary and make the 
/// bottom Run ready and sorted for dequeue 
/// @return True if non-empty bottom rung is ready. False otherwise 
bool LadderQueue::inspectRung() {
    // printf("InspectRung() start\n");
    // print();
    if (bottom) {
        if (bottom->isEmpty()) {
            if (rungs.size()) {
                Rung* lastRung = rungs.back();

                if (bottom == lastRung->bucket[lastRung->minBucket()])
                    lastRung->bucket[lastRung->minBucket()] = nullptr;
            }
            delete bottom;
            bottom = nullptr;
        } else {
            return true;
        }
    }

    if (rungs.size() == 0) {
        if (top->getNodeCount() == 0) // Empty Queue
            return false;
        Rung* r = createNewRungFromTop();
        assert(r != nullptr);
        // print();
        while (!top->isEmpty())
            r->add(top->removeNode());
        // r->updateRCur();
        rungs.push_back(r);
    }
    
    Rung* lastRung = rungs.back();
    int k = lastRung->minBucket();
    while (k <= lastRung->maxBucket() &&
           (lastRung->bucket[k] == nullptr || lastRung->bucket[k]->getNodeCount() == 0)) {
        if (lastRung->bucket[k]) {
            if (bottom == lastRung->bucket[k])
                bottom = nullptr;
            delete lastRung->bucket[k];
            lastRung->bucket[k] = nullptr;
        }
        ++k;
    }

    if (k > lastRung->maxBucket()) { // Empty Rung
        rungs.pop_back();
        delete lastRung;
        return inspectRung();
    }
    // lastRung->rCur = lastRung->rStart + (lastRung->bucketWidth * k);
    // assert(k <= lastRung->maxBucket);
    LinkedEdgeList* bucket_k = lastRung->bucket[k];
    if (bucket_k->getNodeCount() > THRES) {
        // lastRung->minBucket = k;
        // lastRung->updateRCur();
        // print();
        // printf("Bucket overflow (");
        // bucket_k->print();
        // printf(")\n");
        // printf("BottomIdx = %d\n", bottomIndex);

        Rung* r = createNewRungFromBottom(bucket_k);
        assert(r != nullptr);
        while (!bucket_k->isEmpty()) {
            // auto h = bucket_k->getHead()->info;
            // printf("Moving (%d, %f) from bucket_k to the new rung\n", h->idx, h->alpha);
            // printf("Bucket before (");
            // bucket_k->print();
            // printf(")\n");
            // printf("Rung before:\n");
            // r->print();
            r->add(bucket_k->removeNode());
            // printf("Bucket after (");
            // bucket_k->print();
            // printf(")\n");
            // printf("Rung after:\n");
            // r->print();
        }

        delete bucket_k;
        lastRung->bucket[k] = nullptr;

        // lastRung->bucketCount--;
        // lastRung->minBucket = k + 1;
        // lastRung->updateRCur();
        rungs.push_back(r);
        if (rungs.size() > rungused)
            rungused = rungs.size();
        return inspectRung();
    } else {
        bottom = bucket_k;
        bottom->sort();        
        return true;
    }
}

LadderQueue::LadderQueue(int thr) : bottom(nullptr), THRES(thr) {
    top = new LinkedEdgeList;
}

LadderQueue::~LadderQueue() {
    delete top;
    for (auto r : rungs) {
        delete r;
    }
}

Rung* LadderQueue::createNewRungFromTop() {
    return createNewRung(top, true);
}

Rung* LadderQueue::createNewRungFromBottom(LinkedEdgeList* elist) {
    return createNewRung(elist, false);
}

Rung* LadderQueue::createNewRung(LinkedEdgeList* elist, const bool fromTop) {
    // assert(rungs.size() <= MAX_RUNGS);
    // assert(fromTop == (elist == top));
    // assert(!fromTop == (elist == bottom));

    // int rsize = rungs.size();

    double bw;
    double w;
    double rStart;

    if (elist == top) {
        double min = 0;
        // MinTS;
        // if (rsize > 0) {
        //     min = topStart;
        // } else {
        //     min = MinTS;
        // }
        w = getMaxTop() - min;
        // double nc = (double)THRES;
        bw = w / (double)THRES;
        // if (w % THRES != 0)
        //     bw++;
        // if (bw == 0)
        //     return nullptr;
        assert(bw != 0);
        rStart = min;
    } else if (elist == bottom) {

        // when creating a new rung here, check the max range

        Rung* lastRung = rungs.back();
        // int bottomIndex;
        // for (bottomIndex = 0;bottom != lastRung->bucket[bottomIndex];bottomIndex++)
        //     ;
        // rStart = bottom->getFirst()->alpha;        
        // w = lastRung->rCur - rStart; // TODO +1?
        w = lastRung->bucketWidth;
        assert(w > 0);
        bw = w / (double)THRES;
        rStart = lastRung->rCur();
        // if (w % THRES != 0)
        //     bw++;
        // if (bw == 0)
        //     return nullptr; // it is not possible to create a new rung with bw=0
    } else {
        Rung* lastRung = rungs.back();
        w = lastRung->bucketWidth;
        // bw = lastRung->bucketWidth / elist->getNodeCount();
        bw = w / (double)THRES;
        // if (w % THRES != 0)
        //     bw++;
        // if (bw == 0)
        //     return nullptr; // it is not possible to create a new rung with bw=0
        rStart = lastRung->rCur();
    }

    Rung* r = new Rung(w, bw);
    r->rStart = rStart;
    // r->rCur = r->rStart = rStart;
    // printf("Making a new Rung (%f - %f)\n", rStart, rStart + w);

    return r;
}
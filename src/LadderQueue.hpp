#pragma once

#include <iostream>
#include <vector>
#include <deque>
#include <cassert>
#include "defines.h"

// Define the Event class (replace with the actual Event class definition)
class Event {
	Imgidx idx;
	double alpha;

    Event(Imgidx idx, double alpha) : idx(idx), alpha(alpha) {}

    friend class LQNode;
    friend class LinkedEdgeList;
    friend class Rung;
    friend class LadderQueue;
};

struct LQNode {
    LQNode* next;
    LQNode* prev;
    Event* info;

    LQNode() : next(nullptr), prev(nullptr), info(nullptr) {}
};

class LinkedEdgeList {
private:
    LQNode* head;
    LQNode* tail;

    int nodeCount;

public:
    LinkedEdgeList() : head(nullptr), tail(nullptr), nodeCount(0) {}

    void add(Event* e);
    void add(LQNode* n);
    void addInOrder(Event* e); 
    void addInOrder(LQNode* e);
    void addFirst(Event* evt);
    Event* removeDeleteNode();
    LQNode* removeNode();
    void sort();
    bool isAllSameAlpha();

    inline LQNode* getHead() {return head;}
    inline Event* getFirst() {return head->info;}
    inline Event* getLast() {return tail->info;}
    inline int getNodeCount() {return nodeCount;}
    inline bool isEmpty() {return nodeCount == 0;}

    void print();
    // IMPLEMENTATION DONE
    /////////////////////////////////////////////

    ////////////////////////
    // TODO
    // std::string toString();
};

// Define the Rung class
class Rung {
public:
    Rung(double w, double bw);
    ~Rung();
    int getRungIdx(double alpha);
    void checkRungExist(int bidx);
    void addInOrder(Event *evt);
    void add(Event *evt);
    void add(LQNode *n);
    bool print(LinkedEdgeList* bottom);

private:
    double bucketWidth;
    double rStart;
    // double rCur;
    int bucketCount;
    std::vector<LinkedEdgeList*> bucket;
    // int minBucket;
    // int maxBucket;
    int minBucket();
    int maxBucket();
    inline double rCur() { return rStart + minBucket() * bucketWidth; }

    friend class LadderQueue;
};

class LadderQueue {
public:
    LadderQueue(int thr = 64);
    ~LadderQueue();

    // PUSH
    void enqueue(Imgidx idx, double alpha);

    // POP
    Imgidx dequeue();

    Imgidx getTopIdx();
    double getTopAlpha();

    void print();
    
    inline bool isEmpty() { return size == 0; }
    inline int getSize() { return size; }
    inline int getTopInsert() { return topInsert; }
    inline int getRungInsert() { return rungInsert; }
    inline int getBottomInsert() { return bottomInsert; }
    inline int getRungused() { return rungused; }

private:
    // Ladder Queue Structure
    LinkedEdgeList* top;
    std::deque<Rung*> rungs;
    LinkedEdgeList* bottom;

    // Operation
    int THRES;
    int size = 0;
    double MinTS = -1L;
    double MaxTS = -1L;

    // Statistics
    size_t rungused = 0;
    int topInsert = 0;
    int rungInsert = 0;
    int bottomInsert = 0;

    void smartSpawnStatsTopEnqueue(long ts_long, int occurrences, int newsize);
    bool inspectRung();
    Rung* createNewRungFromTop();
    Rung* createNewRungFromBottom(LinkedEdgeList* elist);
    Rung* createNewRung(LinkedEdgeList* elist, bool fromTop);
    inline double getMaxTop() { return MaxTS; }
    double topStart() { return rungs[0] ? (rungs[0]->rStart 
    + rungs[0]->bucketCount * rungs[0]->bucketWidth) : 0; }
};

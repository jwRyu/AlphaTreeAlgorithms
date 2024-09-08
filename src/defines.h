#pragma once
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <vector>

#define DELAYED_NODE_ALLOC 1
#define HQUEUE_COST_AMORTIZE 1

#define NULL_LEVELROOT 0xffffffff
#define NODE_CANDIDATE 0xfffffffe
//
// #define dimg_idx_v(pidx) ((pidx)<<1)
// #define dimg_idx_h(pidx) ((pidx)<<1)+1
//
// #define LEFT_AVAIL(pidx,width)			(((pidx) % (width)) != 0)
// #define RIGHT_AVAIL(pidx,width)			(((pidx) % (width)) != ((width) - 1))
// #define UP_AVAIL(pidx,width)				((pidx) > ((width) - 1))
// #define DOWN_AVAIL(pidx,width,imgsz)		((pidx) < (imgsz) - (width))

#define HIERARCHICAL_QUEUE 0
#define HIERARCHICAL_L1IDX_QUEUE 1
#define HIERARCHICAL_L2IDX_QUEUE 2
#define HEAP_QUEUE 3
#define HEAP_RANK_QUEUE 4
#define TRIE_QUEUE 5
#define TRIE_HYBRID_QUEUE 6

#define _max(a, b) (((a) > (b)) ? (a) : (b))
#define _min(a, b) (((a) > (b)) ? (b) : (a))
#define _clip(x, a, b) _min(_max(x, a), b)

typedef int32_t ImgIdx;
typedef int64_t TrieIdx;

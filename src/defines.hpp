#pragma once
#include <algorithm>
#include <allocator.hpp>
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
#include <map>
#include <omp.h>
#include <optional>
#include <random>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <unordered_map>
#include <vector>
#include <walltime.hpp>

#define DELAYED_NODE_ALLOC 1
#define HQUEUE_COST_AMORTIZE 1

#define NULL_LEVELROOT 0xffffffff
#define NODE_CANDIDATE 0xfffffffe

#define HIERARCHICAL_QUEUE 0
#define HIERARCHICAL_L1IDX_QUEUE 1
#define HIERARCHICAL_L2IDX_QUEUE 2
#define HEAP_QUEUE 3
#define HEAP_RANK_QUEUE 4
#define TRIE_QUEUE 5
#define TRIE_HYBRID_QUEUE 6

#define _max(a, b) (((a) > (b)) ? (a) : (b))
#define _min(a, b) (((a) > (b)) ? (b) : (a))
#define CLIP(x, a, b) _min(_max(x, a), b)

typedef int32_t ImgIdx;
typedef int64_t TrieIdx;

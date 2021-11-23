#pragma once
#include<sys/types.h>

#define DELAYED_NODE_ALLOC		1
#define HQUEUE_COST_AMORTIZE	1

#define NULL_LEVELROOT		0xffffffff
#define NODE_CANDIDATE		0xfffffffe
//
// #define dimg_idx_v(pidx) ((pidx)<<1)
// #define dimg_idx_h(pidx) ((pidx)<<1)+1
//
// #define LEFT_AVAIL(pidx,width)			(((pidx) % (width)) != 0)
// #define RIGHT_AVAIL(pidx,width)			(((pidx) % (width)) != ((width) - 1))
// #define UP_AVAIL(pidx,width)				((pidx) > ((width) - 1))
// #define DOWN_AVAIL(pidx,width,imgsz)		((pidx) < (imgsz) - (width))

#define HIERARCHICAL_QUEUE			0
#define HIERARCHICAL_L1IDX_QUEUE	1
#define HIERARCHICAL_L2IDX_QUEUE	2
#define HEAP_QUEUE					3
#define HEAP_RANK_QUEUE				4
#define TRIE_QUEUE					5
#define TRIE_HYBRID_QUEUE			6

#define _max(a,b) (((a)>(b))?(a):(b))
#define _min(a,b) (((a)>(b))?(b):(a))
#define _clip(x,a,b) _min(_max(x,a),b)

typedef u_int8_t _uint8;
typedef u_int16_t _uint16;
typedef u_int32_t _uint32;
typedef u_int64_t _uint64;
typedef int8_t _int8;
typedef int16_t _int16;
typedef int32_t _int32;
typedef int64_t _int64;


typedef _int64 trieidx;

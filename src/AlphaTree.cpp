#include "AlphaTree.h"

template<class Pixel>
void AlphaTree<Pixel>::BuildAlphaTree(Pixel *img, int height_in, int width_in, int channel_in, int connectivity_in, int algorithm, int numthreads, int tse, double fparam1, double fparam2, int iparam1)
{
    this->height = (Imgidx)height_in;
    this->width = (Imgidx)width_in;
    this->channel = (Imgidx)channel_in;
    this->connectivity = (Imgidx)connectivity_in;
    curSize = 0;

    if (connectivity != 4 && connectivity != 8)
    {
        std::cout << "connectivity should be 4 or 8\n" << std::endl;
        return;
    }

    switch (algorithm)
    {
        case(UNIONFIND):							Unionfind(img);													break;
        case(FLOOD_HIERARQUEUE):					Flood_HierarQueue(img, (HierarQueue*)0, tse);			break;
        case(FLOOD_HIERARQUEUE_CACHE):				Flood_HierarQueue_Cache(img);									break;
        case(FLOOD_TRIE):							Flood_Trie(img,(Trie<trieidx>*)0);						break;
        case(FLOOD_TRIE_CACHE):						Flood_Trie(img,(Trie_Cache<trieidx>*)0);				break;
        case(FLOOD_HEAPQUEUE): 						Flood_HeapQueue(img);											break;
        case(FLOOD_HEAPQUEUE_CACHE):				Flood_HeapQueue_Cache(img);										break;
        case(FLOOD_HIERARQUEUE_HYPERGRAPH):			Flood_HierarQueue_Hypergraph(img);								break;
        case(FLOOD_TRIE_HYPERGRAPH):				Flood_Trie_Hypergraph(img, (Trie<trieidx>*)0);			break;
        case(FLOOD_HIERARHEAPQUEUE):				Flood_HierarHeapQueue(img, fparam1, fparam2, iparam1);			break;
        case(FLOOD_HIERARHEAPQUEUE_CACHE):			Flood_HierarHeapQueue_Cache(img, fparam1, fparam2, iparam1);	break;
        case(FLOOD_HIERARHEAPQUEUE_CACHE_HISTEQ):	Flood_HierarHeapQueue_Cache_histeq(img);						break;
        case(PILOT_RANK):							Pilot_Rank(img, numthreads);									break;
        default: break;
    }
}

template class AlphaTree<_uint8>;
template class AlphaTree<_uint16>;
template class AlphaTree<_uint32>;
template class AlphaTree<_uint64>;
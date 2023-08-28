#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <cfloat>
#include <cassert>
#include "defines.h"
#include "HierarQueue.h"
#include "Trie.h"
#include "HybridQueue.h"
#include "HeapQueue.h"
#include "radixsort_teeninga/sort/radix_sort_parallel.h"
#include "radixsort_teeninga/sort/sort_item.h"

//for debugging
//#include <bitset>
#include "walltime.h"

using namespace pmt;

//Minimum number of pixels for TSE to be used. TSE does not work well with very small images (also saving memory is not so important on small images).
#define TSE_MINSIZE 10000

#define A		1.3901
#define SIGMA	-2.1989
#define B		-0.1906
#define M		0.05

#define	IMGIDX_32BITS		0
#define	IMGIDX_64BITS		1

#define	PIXEL_8BIT			1
#define	PIXEL_16BIT			2
#define	PIXEL_32BIT			4
#define	PIXEL_64BIT			8
#define	PIXEL_FLOAT			16
#define	PIXEL_DOUBLE		32

//algcode
#define UNIONFIND							0//p1
#define FLOOD_HIERARQUEUE					1//p1,2
#define FLOOD_HIERARQUEUE_CACHE				2//p2
#define FLOOD_TRIE							3//p1,2
#define FLOOD_TRIE_CACHE					4//p2
#define FLOOD_HEAPQUEUE						5//p1,2
#define FLOOD_HEAPQUEUE_CACHE				6//p1,2
#define FLOOD_HIERARQUEUE_HYPERGRAPH		7//p1
#define FLOOD_TRIE_HYPERGRAPH				8//p1
#define FLOOD_HIERARHEAPQUEUE_CACHE			9//p2
#define FLOOD_HIERARHEAPQUEUE_CACHE_HISTEQ	10//p2
#define FLOOD_HIERARQUEUE_PAR				11//p1
#define PILOT_RANK							12//p1
#define FLOOD_HIERARQUEUE_CACHE_PAR			13//p2
#define FLOOD_HIERARHEAPQUEUE				14//p2

#define ROOTIDX -1

#define CHKRNG(var,a,b) ( (var >= a) && (var < b) )
#define QUANTIZE_RANK(rank, binsize) (rank) / (binsize)

int nummv = 0, numedge = 0;
int memthod = 0;

extern _uint64 *qrecord;

template<class Imgidx, class Pixel>
class AlphaNode
{
public:
	Imgidx area;
	double alpha;  /* alpha of flat zone */
	double sumPix;
	Pixel minPix;
	Pixel maxPix;

	Imgidx parentidx;
	Imgidx rootidx;

	inline void set(Imgidx area_in, double level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in)
	{
		this->area = area_in;
		this->alpha = level;
		this->sumPix = sumPix_in;
		this->minPix = minPix_in;
		this->maxPix = maxPix_in;
	}

	inline void set(Imgidx area_in, double level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in, _int8 pix_type)
	{
		this->area = area_in;
		this->alpha = level;
		this->sumPix = sumPix_in;
		this->minPix = minPix_in;
		this->maxPix = maxPix_in;
	}

	inline void add(AlphaNode* q, _int8 pix_type)//, _int8 pix_type)
	{
		this->area += q->area;
		this->sumPix += (double)q->sumPix;
		this->minPix = _min(this->minPix, q->minPix);
		this->maxPix = _max(this->maxPix, q->maxPix);
	}
	inline void add(Pixel pix_val, _int8 pix_type)
	{
		this->area++;
		this->sumPix += (double)pix_val;
		this->minPix = _min(this->minPix, pix_val);
		this->maxPix = _max(this->maxPix, pix_val);
	}

	inline void copy(AlphaNode* q)
	{
		this->area = q->area;
		this->sumPix = q->sumPix;
		this->minPix = q->minPix;
		this->maxPix = q->maxPix;
	}

	inline void connect_to_parent(AlphaNode* pPar, Imgidx iPar, _int8 pix_type)
	{
		this->parentidx = iPar;
		pPar->add(this, pix_type);
	}

	void print(AlphaNode* node, _int8 pix_type)
	{
		double val;

		val = (double)this->sumPix;
		
		if (sizeof(Pixel) > 2)
		{
			printf("Node idx: %d  alpha: %f area: %d, sumpix: %.0f, min-max: %d-%d  rootidx: %d  parent: %d\n"
			,(int)(this - node)
			,(double)log2(this->alpha + 1)
			,(int)this->area
			,(double)log2(val + 1)
			,(int)log2(this->minPix + 1)
			,(int)log2(this->maxPix + 1)
			,(int)this->rootidx
			,(int)this->parentidx);			
		}
		else
		{
			printf("Node idx: %d  alpha: %f area: %d, sumpix: %.0f, min-max: %d-%d  rootidx: %d  parent: %d\n"
			,(int)(this - node)
			,(double)this->alpha
			,(int)this->area
			,(double)val
			,(int)this->minPix
			,(int)this->maxPix
			,(int)this->rootidx
			,(int)this->parentidx);
		}		
	}

	void print(AlphaNode* node, int heading, _int8 pix_type)
	{
		double val;

		val = (double)this->sumPix;

		printf("%d: Node idx: %d\talpha: %f\tparent: %d\n               area: %d, sumpix: %f\n               min-max: %d-%d    rootidx: %d\n"
		,heading
		,(int)(this - node)
		,(double)this->alpha
		,(int)this->parentidx
		,(int)this->area
		,(double)this->sumPix
		,(int)this->minPix
		,(int)this->maxPix
	  	,(int)this->rootidx);
	}
};

template<class Imgidx, class Pixel>
class RankItem
{
public:
	Pixel alpha;
	Imgidx dimgidx;
	//Imgidx p, q;

	inline void operator=(const RankItem& item)
	{
		this->alpha = item.alpha;
		this->dimgidx = item.dimgidx;
	}

	inline Imgidx get_pidx0(Imgidx connectivity = 4)
	{
		if (connectivity == 4)
			return (this->dimgidx >> 1);
		else if (connectivity == 8)
			return (this->dimgidx >> 2);
		else
		{
			return -1;
		}
	}
	inline Imgidx get_pidx1(Imgidx width, Imgidx connectivity = 4)
	{
		if (connectivity == 4)
			return (this->dimgidx >> 1) + width + (1 - width) * (this->dimgidx & 1);
		else if (connectivity == 8)
		{
				Imgidx neighboridx = (this->dimgidx & 2);
				return (this->dimgidx >> 2) + width * ((Imgidx)(neighboridx < 2) - (Imgidx)(neighboridx == 3)) + (Imgidx)(neighboridx > 0);
		}
		else
		{
			return -1;
		}
	}
};

template<class Imgidx, class Index, class Pixel>
class ATree
{
public:
	Imgidx maxSize;
	Imgidx curSize;
	Imgidx height, width, channel, connectivity;
	AlphaNode<Imgidx, Pixel> *node, *node_in;
	Imgidx num_node, num_node_in;
	Imgidx rootidx;
	Imgidx* parentAry;

	//double ratio_excessnodes;
	double nrmsd;
	_int8 pix_type;
	_int8 bit_depth;

	//for debugging
	double tijd0, tijd1, tijd2;
	int int0,int1,int2,int3;

	ATree(_int8 ptype, _int8 bit_depth_in) : maxSize(0), curSize(0), node(0), parentAry(0), pix_type(ptype), bit_depth(bit_depth_in) {}
	~ATree()
	{
		 if (node)
		 {
				Free(node);
				if (parentAry)
				{
					Free(parentAry);
				}
			}
	}
	inline void clear() { Free(node); Free(parentAry); node = NULL; parentAry = NULL; curSize = 0; }

	void BuildAlphaTree(Pixel *img, int height_in, int width_in, int channel_in, int connectivity_in, int algorithm, int numthreads, int tse, double fparam1 = 0.0, double fparam2 = 0.0, int iparam1 = 0)
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
			case(FLOOD_HIERARQUEUE):					Flood_HierarQueue(img, (HierarQueue<Imgidx>*)0, tse);			break;
			case(FLOOD_HIERARQUEUE_CACHE):				Flood_HierarQueue_Cache(img);									break;
			case(FLOOD_TRIE):							Flood_Trie(img,(Trie<Imgidx, trieidx>*)0);						break;
			case(FLOOD_TRIE_CACHE):						Flood_Trie(img,(Trie_Cache<Imgidx, trieidx>*)0);				break;
			case(FLOOD_HEAPQUEUE): 						Flood_HeapQueue(img);											break;
			case(FLOOD_HEAPQUEUE_CACHE):				Flood_HeapQueue_Cache(img);										break;
			case(FLOOD_HIERARQUEUE_HYPERGRAPH):			Flood_HierarQueue_Hypergraph(img);								break;
			case(FLOOD_TRIE_HYPERGRAPH):				Flood_Trie_Hypergraph(img, (Trie<Imgidx, trieidx>*)0);			break;
			case(FLOOD_HIERARHEAPQUEUE):				Flood_HierarHeapQueue(img, fparam1, fparam2, iparam1);			break;
			case(FLOOD_HIERARHEAPQUEUE_CACHE):			Flood_HierarHeapQueue_Cache(img, fparam1, fparam2, iparam1);	break;
			case(FLOOD_HIERARHEAPQUEUE_CACHE_HISTEQ):	Flood_HierarHeapQueue_Cache_histeq(img);						break;
			case(FLOOD_HIERARQUEUE_PAR):				Flood_Hierarqueue_par_tse(img, numthreads, tse);				break;
			case(PILOT_RANK):							Pilot_Rank(img, numthreads);									break;
			default: break;
		}
	}

	void AlphaFilter(double *outimg, double alpha)
	{
		Imgidx idx_lim, i, imgsize;
		AlphaNode<Imgidx, Pixel> *pNode;

		alpha = node[rootidx].alpha * alpha;

		imgsize = height * width;
		//val = 1;
		for (i = 0; i < imgsize; i++)
		{
			pNode = parentAry ? &node[parentAry[i]] : &node[i];
			while (pNode->parentidx != -1 && pNode->alpha < alpha)
				pNode = &node[pNode->parentidx];
			outimg[i] = (double)pNode->area;
		}
	}

	void AreaFilter(double *outimg, double area)
	{
		Imgidx i, imgsize;
		Imgidx iarea;
		AlphaNode<Imgidx, Pixel> *pNode;

		imgsize = height * width;
		iarea = (Imgidx)(area * (double)imgsize);
		iarea = _min(imgsize, _max(0, iarea));
		//val = 1;
		for (i = 0; i < imgsize; i++)
		{
			pNode = parentAry ? &node[parentAry[i]] : &node[i];
			while (pNode->parentidx != -1 && pNode->area < iarea)
				pNode = &node[pNode->parentidx];
			outimg[i] = (double)pNode->alpha;
		}
	}

	void print_tree()
	{
		for (int i = 0; i < maxSize; i++)
			if (node[i].area) node[i].print(node, pix_type);
	}

private:
	inline Pixel abs_diff(Pixel p, Pixel q)
	{
		if (p > q)
			return p - q;
		else
			return q - p;
	}

	inline _uint8 compute_incidedge_queue(Pixel d0, Pixel d1)
	{
		if (d0 <= d1)
			return 0x4;
		else
			return 0x1;
	}

#define _COMPUTE_RANK_PIXEL \
	rankitem[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))\
			: (Pixel)(abs((_int64)img[imgidx + width] - (_int64)img[imgidx]));\
	rankitem[contidx++].dimgidx = dimgidx++;\
	rankitem[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))\
			: (Pixel)(abs((_int64)img[imgidx + 1] - (_int64)img[imgidx]));\
	rankitem[contidx++].dimgidx = dimgidx++;

#define _COMPUTE_RANK_PIXEL_LASTCOL \
	rankitem[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + width] - (double)img[imgidx]))\
		: (Pixel)(abs((_int64)img[imgidx + width] - (_int64)img[imgidx]));\
	rankitem[contidx++].dimgidx = dimgidx;\
	dimgidx += 2;

#define _COMPUTE_RANK_PIXEL_LASTROW \
	dimgidx++;\
	rankitem[contidx].alpha = p64 ? (Pixel)(abs((double)img[imgidx + 1] - (double)img[imgidx]))\
		: (Pixel)(abs((_int64)img[imgidx + 1] - (_int64)img[imgidx]));\
	rankitem[contidx++].dimgidx = dimgidx++;

	template <typename Value>
	void compute_dimg_par4(RankItem<Imgidx, double>*& rankitem, Pixel* img, SortValue<Value>*& vals)
	{
		Imgidx contidx, dimgidx, imgidx, i, j;

		contidx = imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			if (channel == 1)
			{
				#pragma omp parallel for schedule(guided,1) private(imgidx, dimgidx, contidx, i, j)
				for (i = 0; i < height; i++)
				{
					double d;
					imgidx = i * width;
					dimgidx = imgidx << (connectivity >> 2);
					contidx = dimgidx - i;
					if (i < height - 1)
					{
						for (j = 0; j < width - 1; j++)
						{
							//caclulate histogram here
							vals[contidx].val_ = abs_diff(img[imgidx + width], img[imgidx]);
							d = (double)(vals[contidx].val_);
							rankitem[contidx].alpha = d;
							rankitem[contidx++].dimgidx = dimgidx++;

							vals[contidx].val_ = abs_diff(img[imgidx + 1], img[imgidx]);
							d = (double)(vals[contidx].val_);
							rankitem[contidx].alpha = d;
							rankitem[contidx++].dimgidx = dimgidx++;
							imgidx++;
						}
						vals[contidx].val_ = abs_diff(img[imgidx + width], img[imgidx]);
						d = (double)(vals[contidx].val_);
						rankitem[contidx].alpha = d;
						rankitem[contidx++].dimgidx = dimgidx;
						dimgidx += 2;
						imgidx++;
					}
					else
					{
						for (j = 0; j < width - 1; j++)
						{
							dimgidx++;
							vals[contidx].val_ = abs_diff(img[imgidx + 1], img[imgidx]);
							d = (double)(vals[contidx].val_);
							rankitem[contidx].alpha = d;
							rankitem[contidx++].dimgidx = dimgidx++;
							imgidx++;
						}
					}
				}
			}
			else
			{
				Imgidx chstride = channel, chstride2 = channel * 2;
				Imgidx imgsize = height * width;
				Imgidx dimgsize = height * width * (connectivity >> 1);
				double* dimg_3ch = (double*)Calloc(channel * dimgsize * sizeof(double));
				Imgidx linestride = width * (connectivity >> 1) - (connectivity >> 2);

				#pragma omp parallel for schedule(guided,1) private(imgidx, dimgidx, i, j)
				for (int hch = 0;hch < channel * height;hch++)
				{
					double d;
					int ch = hch / height;
					i = hch % height;
					Pixel *pimg = img + ch * imgsize + i * width;
					double *pdimg = dimg_3ch + ch + i * width * (connectivity >> 1) * channel;
					imgidx = dimgidx = 0;
					if (i < height - 1)
					{
						for (j = 0; j < width - 1; j++)
						{
							d = (double)pimg[imgidx + width] - (double)pimg[imgidx];
							pdimg[dimgidx] = d * d;
							dimgidx += chstride;
							d = (double)pimg[imgidx + 1] - (double)pimg[imgidx];
							pdimg[dimgidx] = d * d;
							dimgidx += chstride;
							imgidx++;
						}
						d = (double)pimg[imgidx + width] - (double)pimg[imgidx];
						pdimg[dimgidx] = d * d;
						dimgidx += chstride2;
						imgidx++;
					}
					else
					{
						for (j = 0; j < width - 1; j++)
						{
							dimgidx += chstride;
							d = (double)pimg[imgidx + 1] - (double)pimg[imgidx];
							pdimg[dimgidx] = d * d;
							dimgidx += chstride;
							imgidx++;
						}
					}
				}

				#pragma omp parallel for schedule(guided,1) private(imgidx, dimgidx, contidx, i, j)
				for (i = 0; i < height; i++)
				{
					double d;
					double *pdimg = dimg_3ch + i * width * (connectivity >> 1) * channel;
					contidx = i * linestride;
					imgidx = i * width;
					Imgidx dimgidx_3 = 0;
					dimgidx = i * width * (connectivity >> 1);
					if (i < height - 1)
					{
						for (j = 0; j < width - 1; j++)
						{
							d = pdimg[dimgidx_3] + pdimg[dimgidx_3 + 1] + pdimg[dimgidx_3 + 2]; //only for 3-ch
							rankitem[contidx].alpha = vals[contidx].val_ = sqrt(d/(double)channel);
							rankitem[contidx].dimgidx = dimgidx;
							contidx++;

							dimgidx_3 += chstride;
							dimgidx++;

							d = pdimg[dimgidx_3] + pdimg[dimgidx_3 + 1] + pdimg[dimgidx_3 + 2]; //only for 3-ch
							rankitem[contidx].alpha = vals[contidx].val_ = sqrt(d/(double)channel);
							rankitem[contidx].dimgidx = dimgidx;
							contidx++;

							dimgidx_3 += chstride;
							dimgidx++;

							imgidx++;
						}
						d = pdimg[dimgidx_3] + pdimg[dimgidx_3 + 1] + pdimg[dimgidx_3 + 2]; //only for 3-ch
						rankitem[contidx].alpha = vals[contidx].val_ = sqrt(d/(double)channel);
						rankitem[contidx].dimgidx = dimgidx;
						contidx++;
						dimgidx++;

						dimgidx_3 += chstride2;
						imgidx++;
					}
					else
					{
						for (j = 0; j < width - 1; j++)
						{
							dimgidx++;
							dimgidx_3 += chstride;

							d = pdimg[dimgidx_3] + pdimg[dimgidx_3 + 1] + pdimg[dimgidx_3 + 2]; //only for 3-ch
							rankitem[contidx].alpha = vals[contidx].val_ = sqrt(d/(double)channel);
							rankitem[contidx].dimgidx = dimgidx;
							contidx++;

							dimgidx_3 += chstride;
							dimgidx++;
							imgidx++;
						}
					}
				}
				Free(dimg_3ch);
			}
		}
		else if (connectivity == 8) //not implemented yet
		{
			if (channel == 1)
			{
				#pragma omp parallel for schedule(guided,1) private(imgidx, dimgidx, contidx, i, j)
				for (i = 0; i < height; i++)
				{
					double d;
					imgidx = i * width;
					dimgidx = imgidx << (connectivity >> 2);
					contidx = dimgidx - ((connectivity == 4) ? (i) : (3 * i));
					if (connectivity == 8 && i > 0)
						contidx -= width - 1;

					if (i < height - 1)
					{
						for (j = 0; j < width - 1; j++)
						{
							//caclulate histogram here
							vals[contidx].val_ = abs_diff(img[imgidx + width], img[imgidx]);
							d = (double)(vals[contidx].val_);
							rankitem[contidx].alpha = d;
							rankitem[contidx++].dimgidx = dimgidx++;

							vals[contidx].val_ = abs_diff(img[imgidx + width + 1], img[imgidx]);
							d = (double)(vals[contidx].val_);
							rankitem[contidx].alpha = d;
							rankitem[contidx++].dimgidx = dimgidx++;

							vals[contidx].val_ = abs_diff(img[imgidx + 1], img[imgidx]);
							d = (double)(vals[contidx].val_);
							rankitem[contidx].alpha = d;
							rankitem[contidx++].dimgidx = dimgidx++;

							if (i > 0)
							{
								vals[contidx].val_ = abs_diff(img[imgidx - width + 1], img[imgidx]);
								d = (double)(vals[contidx].val_);
								rankitem[contidx].alpha = d;
								rankitem[contidx++].dimgidx = dimgidx;
							}
							dimgidx++;

							imgidx++;
						}

						vals[contidx].val_ = abs_diff(img[imgidx + width], img[imgidx]);
						d = (double)(vals[contidx].val_);
						rankitem[contidx].alpha = d;
						rankitem[contidx++].dimgidx = dimgidx;
						dimgidx += 4;
						imgidx++;
					}
					else
					{
						for (j = 0; j < width - 1; j++)
						{
							dimgidx += 2;
							vals[contidx].val_ = abs_diff(img[imgidx + 1], img[imgidx]);
							d = (double)(vals[contidx].val_);
							rankitem[contidx].alpha = d;
							rankitem[contidx++].dimgidx = dimgidx++;

							vals[contidx].val_ = abs_diff(img[imgidx - width + 1], img[imgidx]);
							d = (double)(vals[contidx].val_);
							rankitem[contidx].alpha = d;
							rankitem[contidx++].dimgidx = dimgidx++;

							imgidx++;
						}
					}
				}
			}
		}
	}

	void compute_dimg1(Imgidx* rank, RankItem<Imgidx, double>*& rankitem, Pixel* img)
	{
		Imgidx contidx, dimgidx, imgidx, nbyte, i, j, nredges;
		Imgidx *hist, *h;
		Imgidx hsum;
		RankItem<Imgidx, double> *tmp, *r;
		Pixel hidx, h_offset, mask = (Pixel)0xffff, shamt;
		size_t hist_size = 65536;
		Pixel *incidents = (Pixel*)Malloc(connectivity * width * height * sizeof(Pixel));

		if (connectivity == 4)
			nredges = (width - 1) * height + width * (height - 1);
		else
			nredges = (width - 1) * height + width * (height - 1) + 2 * (width - 1) * (height - 1);

		tmp = (RankItem<Imgidx, double>*)Malloc(nredges * sizeof(RankItem<Imgidx, double>));
		i = (sizeof(Pixel) >> 1) * 65536;
		hist = (Imgidx*)Malloc(_max(hist_size * sizeof(Imgidx) * sizeof(Pixel), (sizeof(Pixel) * 65536 * sizeof(Imgidx)) >> 1));

		contidx = imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			if (channel == 1)
			{
				for (i = 0; i < height - 1; i++)
				{
					for (j = 0; j < width - 1; j++)
					{
							rankitem[contidx].alpha = abs((double)img[imgidx + width] - (double)img[imgidx]);
							rankitem[contidx++].dimgidx = dimgidx++;
							rankitem[contidx].alpha = abs((double)img[imgidx + 1] - (double)img[imgidx]);
							rankitem[contidx++].dimgidx = dimgidx++;
							imgidx++;
					}
					rankitem[contidx].alpha = abs((double)img[imgidx + width] - (double)img[imgidx]);
					rankitem[contidx++].dimgidx = dimgidx;
					dimgidx += 2;
					imgidx++;
				}
				for (j = 0; j < width - 1; j++)
				{
						dimgidx++;
						rankitem[contidx].alpha = abs((double)img[imgidx + 1] - (double)img[imgidx]);
						rankitem[contidx++].dimgidx = dimgidx++;
					imgidx++;
				}
			}
			else
			{
				double d;
				Pixel *pimg = img;
				Imgidx imgsize = height * width;
				double* dimg = (double*)Calloc(2 * imgsize * sizeof(double));
				contidx = 0;
				for (int ch = 0;ch < channel;ch++)
				{
					imgidx = dimgidx = 0;
					for (i = 0; i < height - 1; i++)
					{
						for (j = 0; j < width - 1; j++)
						{
							d = (double)pimg[imgidx + width] - (double)pimg[imgidx];
							if (ch < channel - 1)	dimg[dimgidx] += d * d;
							else
							{
								rankitem[contidx].alpha = sqrt((dimg[dimgidx] + d * d)/(double)channel);
								rankitem[contidx++].dimgidx = dimgidx;
							}
							dimgidx++;

							d = (double)pimg[imgidx + 1] - (double)pimg[imgidx];
							if (ch < channel - 1) dimg[dimgidx] += d * d;
							else
							{
									rankitem[contidx].alpha = sqrt((dimg[dimgidx] + d * d)/(double)channel);
									rankitem[contidx++].dimgidx = dimgidx;
							}
							dimgidx++;
							imgidx++;
						}
						d = (double)pimg[imgidx + width] - (double)pimg[imgidx];
						if (ch < channel - 1) dimg[dimgidx] += d * d;
						else
						{
								rankitem[contidx].alpha = sqrt((dimg[dimgidx] + d * d)/(double)channel);
								rankitem[contidx++].dimgidx = dimgidx;
						}
						dimgidx += 2;
						imgidx++;
					}
					for (j = 0; j < width - 1; j++)
					{
						dimgidx++;
						d = (double)pimg[imgidx + 1] - (double)pimg[imgidx];
						if (ch < channel - 1) dimg[dimgidx] += d * d;
						else
						{
								rankitem[contidx].alpha = sqrt((dimg[dimgidx] + d * d)/(double)channel);
								rankitem[contidx++].dimgidx = dimgidx;
						}
						dimgidx++;
						imgidx++;
					}
					pimg += imgsize;
				}
				Free(dimg);
			}
		}
		else if (connectivity == 8)
		{
		}


		for (i = 0; i < (Imgidx)(sizeof(Pixel) * 65536 >> 1); i++)
			hist[i] = 0;

		if (bit_depth <= 16) // count sort
		{
			for (i = 0; i < nredges; i++)
			{
				hidx = (Pixel)rankitem[i].alpha;
				hist[hidx]++;
			}

			h = hist;

			hsum = 0;
			for (i = 0; i < (Imgidx)(hist_size); i++)
			{
				hsum += h[i];
				h[i] = hsum;
			}
			for (i = nredges - 1; i >= 0; i--)
			{
				hidx = (Pixel)rankitem[i].alpha;
				j = --h[hidx];
				tmp[j] = rankitem[i]; //slow!
			}
			r = rankitem; rankitem = tmp; tmp = r; //swap p,q
		}
		else //Radix sort
		{
			shamt = (Pixel)16;
			for (i = 0; i < nredges; i++)
			{
				hidx = (Pixel)rankitem[i].alpha;
				hist[hidx & mask]++;

				h_offset = (Pixel)65536;
				for (nbyte = 2; nbyte < (Imgidx)(sizeof(Pixel)); nbyte += 2)
				{
					hidx = hidx >> shamt;
					hist[h_offset + (hidx & mask)]++;
					h_offset += (Pixel)65536;
				}
			}

			h = hist;
			shamt = 0;
			for (nbyte = 0; nbyte < (Imgidx)(sizeof(Pixel)); nbyte += 2)
			{
				hsum = 0;
				for (i = 0; i < (Imgidx)(hist_size); i++)
				{
					hsum += h[i];
					h[i] = hsum;
				}
				for (i = nredges - 1; i >= 0; i--)
				{
					hidx = ((Pixel)(rankitem[i].alpha) >> shamt) & mask;
					j = --h[hidx];
					tmp[j] = rankitem[i]; //slow!
				}
				r = rankitem; rankitem = tmp; tmp = r; //swap p,q
				shamt += 16;
				h += hist_size;
			}
		}
		for (i = 0; i < nredges; i++)
		{
			rank[rankitem[i].dimgidx] = i;
		}

		Free(tmp);
		Free(hist);
		Free(incidents);
	}

	void compute_dimg(Imgidx* rank, RankItem<Imgidx, Pixel>*& rankitem, Pixel* img)
	{
		Imgidx contidx, dimgidx, imgidx, nbyte, i, j, nredges;
		Imgidx *hist, *h;
		Imgidx hsum;
		RankItem<Imgidx, Pixel> *tmp, *r;
		Pixel hidx, h_offset, mask = (Pixel)0xffff, shamt;
		size_t hist_size = 65536;
		_int8 p64 = (pix_type == PIXEL_64BIT) || (pix_type == PIXEL_DOUBLE);
		Pixel *incidents = (Pixel*)Malloc(connectivity * width * height * sizeof(Pixel));

		if (connectivity == 4)
			nredges = (width - 1) * height + width * (height - 1);
		else
			nredges = (width - 1) * height + width * (height - 1) + 2 * (width - 1) * (height - 1);

		tmp = (RankItem<Imgidx, Pixel>*)Malloc(nredges * sizeof(RankItem<Imgidx, Pixel>));
		i = (sizeof(Pixel) >> 1) * 65536;
		hist = (Imgidx*)Malloc(_max(hist_size * sizeof(Imgidx) * sizeof(Pixel), (sizeof(Pixel) * 65536 * sizeof(Imgidx)) >> 1));


		contidx = imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					_COMPUTE_RANK_PIXEL
						imgidx++;
				}
				_COMPUTE_RANK_PIXEL_LASTCOL
					imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				_COMPUTE_RANK_PIXEL_LASTROW
				imgidx++;
			}
		}
		else if (connectivity == 8)
		{
			//   -  -  3
			//   -  p  2
			//   -  0  1
			//top,middle
			//#pragma omp parallel for schedule(guided,1) private(imgidx, dimgidx, contidx, i, j, runtime)
			for (i = 0; i < height; i++)
			{
				double d;
				imgidx = i * width;
				dimgidx = imgidx << (connectivity >> 2);
				contidx = dimgidx - ((connectivity == 4) ? (i) : (3 * i));
				if (connectivity == 8 && i > 0)
					contidx -= width - 1;

				if (i < height - 1)
				{
					for (j = 0; j < width - 1; j++)
					{
						//caclulate histogram here
						d = abs_diff(img[imgidx + width], img[imgidx]);
						rankitem[contidx].alpha = d;
						rankitem[contidx++].dimgidx = dimgidx++;

						d = abs_diff(img[imgidx + width + 1], img[imgidx]);
						rankitem[contidx].alpha = d;
						rankitem[contidx++].dimgidx = dimgidx++;

						d = abs_diff(img[imgidx + 1], img[imgidx]);
						rankitem[contidx].alpha = d;
						rankitem[contidx++].dimgidx = dimgidx++;

						if (i > 0)
						{
							d = abs_diff(img[imgidx - width + 1], img[imgidx]);
							rankitem[contidx].alpha = d;
							rankitem[contidx++].dimgidx = dimgidx;
						}
						dimgidx++;

						imgidx++;
					}

					d = abs_diff(img[imgidx + width], img[imgidx]);
					rankitem[contidx].alpha = d;
					rankitem[contidx++].dimgidx = dimgidx;
					dimgidx += 4;
					imgidx++;
				}
				else
				{
					for (j = 0; j < width - 1; j++)
					{
						dimgidx += 2;
						d = abs_diff(img[imgidx + 1], img[imgidx]);
						rankitem[contidx].alpha = d;
						rankitem[contidx++].dimgidx = dimgidx++;

						d = abs_diff(img[imgidx - width + 1], img[imgidx]);
						rankitem[contidx].alpha = d;
						rankitem[contidx++].dimgidx = dimgidx++;

						imgidx++;
					}
				}
			}
		}

		for (i = 0; i < (Imgidx)(sizeof(Pixel) * 65536 >> 1); i++)
			hist[i] = 0;

		if (bit_depth <= 16) // count sort
		{
			for (i = 0; i < nredges; i++)
			{
				hidx = rankitem[i].alpha;
				hist[hidx]++;
			}

			h = hist;

			hsum = 0;
			for (i = 0; i < (Imgidx)(hist_size); i++)
			{
				hsum += h[i];
				h[i] = hsum;
			}
			for (i = nredges - 1; i >= 0; i--)
			{
				hidx = rankitem[i].alpha;
				j = --h[hidx];
				tmp[j] = rankitem[i]; //slow!
			}
			r = rankitem; rankitem = tmp; tmp = r; //swap p,q
		}
		else //Radix sort
		{
			for (i = 0; i < nredges; i++)
			{
				hidx = rankitem[i].alpha;
				hist[hidx & mask]++;

				shamt = (Pixel)16;
				h_offset = (Pixel)65536;
				for (nbyte = 2; nbyte < (Imgidx)(sizeof(Pixel)); nbyte += 2)
				{
					hidx = hidx >> shamt;
					hist[h_offset + (hidx & mask)]++;
					h_offset += (Pixel)65536;
				}
			}

			h = hist;
			shamt = 0;
			for (nbyte = 0; nbyte < (Imgidx)(sizeof(Pixel)); nbyte += 2)
			{
				hsum = 0;
				for (i = 0; i < (Imgidx)(hist_size); i++)
				{
					hsum += h[i];
					h[i] = hsum;
				}
				for (i = nredges - 1; i >= 0; i--)
				{
					hidx = (rankitem[i].alpha >> shamt) & mask;
					j = --h[hidx];
					tmp[j] = rankitem[i]; //slow!
				}
				r = rankitem; rankitem = tmp; tmp = r; //swap p,q
				shamt += 16;
				h += hist_size;
			}
		}
		for (i = 0; i < nredges; i++)
		{
			rank[rankitem[i].dimgidx] = i;
		}

		Free(tmp);
		Free(hist);
		Free(incidents);
	}

	template <class Value>
	void compute_dimg(Value* dimg, Pixel* img)
	{
		Imgidx dimgidx, imgidx, stride_w = width, i, j;
		Imgidx imgsize = width * height;

		Pixel *pimg = img;

		imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			if (channel == 1)
			{
				for (i = 0; i < height - 1; i++)
				{
					for (j = 0; j < width - 1; j++)
					{
						dimg[dimgidx] = (Value)(abs_diff(img[imgidx + stride_w], img[imgidx]));
						dimgidx++;
						dimg[dimgidx] = (Value)(abs_diff(img[imgidx + 1], img[imgidx]));
						dimgidx++;
						imgidx++;
					}
					dimg[dimgidx] = (Value)(abs_diff(img[imgidx + stride_w], img[imgidx]));
					dimgidx+=2;
					imgidx++;
				}
				for (j = 0; j < width - 1; j++)
				{
					dimgidx++;
					dimg[dimgidx] = (Value)(abs_diff(img[imgidx + 1], img[imgidx]));
					dimgidx++;
					imgidx++;
				}
			}
			else
			{
				double d;
				for (int ch = 0;ch < channel;ch++)
				{
					imgidx = dimgidx = 0;
					for (i = 0; i < height - 1; i++)
					{
						for (j = 0; j < width - 1; j++)
						{
							d = (double)pimg[imgidx + stride_w] - (double)pimg[imgidx];
							if (ch == 0) 						   dimg[dimgidx] = d * d;
							else if (ch != channel - 1) dimg[dimgidx] += d * d;
							else 											 dimg[dimgidx] = sqrt(dimg[dimgidx] + d * d);
							dimgidx++;

							d = (double)pimg[imgidx + 1] - (double)pimg[imgidx];
							if (ch == 0) 						   dimg[dimgidx] = d * d;
							else if (ch != channel - 1) dimg[dimgidx] += d * d;
							else 											 dimg[dimgidx] = sqrt(dimg[dimgidx] + d * d);
							dimgidx++;
							imgidx++;
						}
						d = (double)pimg[imgidx + stride_w] - (double)pimg[imgidx];
						if (ch == 0) 						   dimg[dimgidx] = d * d;
						else if (ch != channel - 1) dimg[dimgidx] += d * d;
						else 											 dimg[dimgidx] = sqrt(dimg[dimgidx] + d * d);
						dimgidx+=2;
						imgidx++;
					}
					for (j = 0; j < width - 1; j++)
					{
						dimgidx++;
						d = (double)pimg[imgidx + 1] - (double)pimg[imgidx];
						if (ch == 0) 						   dimg[dimgidx] = d * d;
						else if (ch != channel - 1) dimg[dimgidx] += d * d;
						else 											 dimg[dimgidx] = sqrt(dimg[dimgidx] + d * d);
						dimgidx++;
						imgidx++;
					}
					pimg += imgsize;
				}
			}
		}
		else if (connectivity == 8)
		{
			if (channel == 1)
			{
				//   -  -  3
				//   -  p  2
				//   -  0  1
				//top,middle
				for (i = 0; i < height - 1; i++)
				{
					for (j = 0; j < width - 1; j++)
					{
						dimg[dimgidx++] = (Value)(abs_diff((_int64)img[imgidx + width], (_int64)img[imgidx]));//0
						dimg[dimgidx++] = (Value)(abs_diff((_int64)img[imgidx + width + 1], (_int64)img[imgidx]));//1
						dimg[dimgidx++] = (Value)(abs_diff((_int64)img[imgidx + 1], (_int64)img[imgidx]));//2
						if (i > 0)
						{
							dimg[dimgidx] = (Value)(abs_diff((_int64)img[imgidx - width + 1], (_int64)img[imgidx]));//3
						}
						dimgidx++;
						imgidx++;
					}
					dimg[dimgidx] = (Value)(abs_diff((_int64)img[imgidx + width], (_int64)img[imgidx]));//0
					dimgidx += 4;//skip 1,2,3
					imgidx++;
				}

				//bottom
				dimgidx += 2; //skip 0,1
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx++] = (Value)(abs_diff((_int64)img[imgidx + 1], (_int64)img[imgidx]));//2
					dimg[dimgidx] = (Value)(abs_diff((_int64)img[imgidx - width + 1], (_int64)img[imgidx]));//3
					dimgidx += 3;
					imgidx++;
				}
			}

		}
	}

	Pixel compute_dimg1(Pixel* dimg, Imgidx* dhist, Pixel* img)
	{
		Imgidx dimgidx, imgidx, stride_w = width, i, j;
		Pixel maxdiff = 0;

		imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + stride_w] - (_int64)img[imgidx]));
					maxdiff = _max(maxdiff,dimg[dimgidx]);
					dhist[dimg[dimgidx++]]++;

					dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + 1] - (_int64)img[imgidx]));
					maxdiff = _max(maxdiff,dimg[dimgidx]);
					dhist[dimg[dimgidx++]]++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + stride_w] - (_int64)img[imgidx]));
				maxdiff = _max(maxdiff,dimg[dimgidx]);
				dhist[dimg[dimgidx++]]++;
				dimgidx++;
				imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				dimgidx++;
				dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + 1] - (_int64)img[imgidx]));
				maxdiff = _max(maxdiff,dimg[dimgidx]);
				dhist[dimg[dimgidx++]]++;
				imgidx++;
			}
		}
		else if (connectivity == 8)
		{
			//   -  -  3
			//   -  p  2
			//   -  0  1
			//top,middle
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + width] - (_int64)img[imgidx]));//0
					maxdiff = _max(maxdiff,dimg[dimgidx]);
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + width + 1] - (_int64)img[imgidx]));//1
					maxdiff = _max(maxdiff,dimg[dimgidx]);
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + 1] - (_int64)img[imgidx]));//2
					maxdiff = _max(maxdiff,dimg[dimgidx]);
					dhist[dimg[dimgidx++]]++;
					if (i > 0)
					{
						dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx - width + 1] - (_int64)img[imgidx]));//3
						maxdiff = _max(maxdiff,dimg[dimgidx]);
						dhist[dimg[dimgidx]]++;
					}
					dimgidx++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + width] - (_int64)img[imgidx]));//0
				maxdiff = _max(maxdiff,dimg[dimgidx]);
				dhist[dimg[dimgidx]]++;
				dimgidx += 4;//skip 1,2,3
				imgidx++;
			}

			//bottom
			dimgidx += 2; //skip 0,1
			for (j = 0; j < width - 1; j++)
			{
				dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + 1] - (_int64)img[imgidx]));//2
				maxdiff = _max(maxdiff,dimg[dimgidx]);
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx - width + 1] - (_int64)img[imgidx]));//3
				maxdiff = _max(maxdiff,dimg[dimgidx]);
				dhist[dimg[dimgidx]]++;
				dimgidx += 3;
				imgidx++;
			}
		}

		return maxdiff;
	}

	void compute_dimg(Imgidx &minidx, double &mindiff, Pixel* dimg, Imgidx* dhist, Pixel* img, double a)
	{
		Imgidx dimgidx, imgidx, stride_w = width, i, j;
		int hidx;
		mindiff = (double)((Pixel)(-1));
		imgidx = dimgidx = 0;

		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + stride_w], img[imgidx]));
					if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
					hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
					dhist[hidx]++;
					dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));
					if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
					hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
					dhist[hidx]++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + stride_w], img[imgidx]));
				if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
				hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
				dhist[hidx]++;
				dimgidx++;
				imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				dimgidx++;
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));
				if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
				hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
				dhist[hidx]++;
				imgidx++;
			}
		}
		else if (connectivity == 8)
		{
			//   -  -  3
			//   -  p  2
			//   -  0  1
			//top,middle
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + width], img[imgidx]));//0
					if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
					hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
					dhist[hidx]++;
					dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + width + 1], img[imgidx]));//1
					if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
					hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
					dhist[hidx]++;
					dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));//2
					if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
					hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
					dhist[hidx]++;
					if (i > 0)
					{
						dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx - width + 1], img[imgidx]));//3
						if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
						hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
						dhist[hidx]++;
					}
					dimgidx++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + width], img[imgidx]));//0
				if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
				hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
				dhist[hidx]++;
				dimgidx += 4;//skip 1,2,3
				imgidx++;
			}

			//bottom
			dimgidx += 2; //skip 0,1
			for (j = 0; j < width - 1; j++)
			{
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));//2
				if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
				hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
				dhist[hidx]++;
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx - width + 1], img[imgidx]));//3
				if (dimg[dimgidx] < mindiff) {mindiff = dimg[dimgidx]; minidx = i * width + j;}
				hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
				dhist[hidx]++;
				dimgidx += 3;
				imgidx++;
			}
		}
	}


	void compute_dimg(Pixel* dimg, Imgidx* dhist, Pixel* img, double a)
	{
		Imgidx dimgidx, imgidx, stride_w = width, i, j;
		int hidx;
		imgidx = dimgidx = 0;

		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + stride_w], img[imgidx]));
					hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
					dhist[hidx]++;
					dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));
					hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
					dhist[hidx]++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + stride_w], img[imgidx]));
				hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
				dhist[hidx]++;
				dimgidx++;
				imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				dimgidx++;
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));
				hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
				dhist[hidx]++;
				imgidx++;
			}
		}
		else if (connectivity == 8)
		{
			//   -  -  3
			//   -  p  2
			//   -  0  1
			//top,middle
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + width], img[imgidx]));//0
					hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
					dhist[hidx]++;
					dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + width + 1], img[imgidx]));//1
					hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
					dhist[hidx]++;
					dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));//2
					hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
					dhist[hidx]++;
					if (i > 0)
					{
						dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx - width + 1], img[imgidx]));//3
						hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
						dhist[hidx]++;
					}
					dimgidx++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + width], img[imgidx]));//0
				hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
				dhist[hidx]++;
				dimgidx += 4;//skip 1,2,3
				imgidx++;
			}

			//bottom
			dimgidx += 2; //skip 0,1
			for (j = 0; j < width - 1; j++)
			{
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));//2
				hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
				dhist[hidx]++;
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx - width + 1], img[imgidx]));//3
				hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
				dhist[hidx]++;
				dimgidx += 3;
				imgidx++;
			}
		}
	}

	void compute_dimg(Pixel* dimg, Imgidx* dhist, Pixel* img)
	{
		Imgidx dimgidx, imgidx, stride_w = width, i, j;

		imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + stride_w] - (_int64)img[imgidx]));
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + 1] - (_int64)img[imgidx]));
					dhist[dimg[dimgidx++]]++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + stride_w] - (_int64)img[imgidx]));
				dhist[dimg[dimgidx++]]++;
				dimgidx++;
				imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				dimgidx++;
				dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + 1] - (_int64)img[imgidx]));
				dhist[dimg[dimgidx++]]++;
				imgidx++;
			}
		}
		else if (connectivity == 8)
		{
			//   -  -  3
			//   -  p  2
			//   -  0  1
			//top,middle
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + width] - (_int64)img[imgidx]));//0
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + width + 1] - (_int64)img[imgidx]));//1
					dhist[dimg[dimgidx++]]++;
					dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + 1] - (_int64)img[imgidx]));//2
					dhist[dimg[dimgidx++]]++;
					if (i > 0)
					{
						dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx - width + 1] - (_int64)img[imgidx]));//3
						dhist[dimg[dimgidx]]++;
					}
					dimgidx++;
					imgidx++;
				}
				dimg[dimgidx] = (Pixel)(abs((_int64)img[imgidx + width] - (_int64)img[imgidx]));//0
				dhist[dimg[dimgidx]]++;
				dimgidx += 4;//skip 1,2,3
				imgidx++;
			}

			//bottom
			dimgidx += 2; //skip 0,1
			for (j = 0; j < width - 1; j++)
			{
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));//2
				dhist[dimg[dimgidx++]]++;
				dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx - width + 1], img[imgidx]));//3
				dhist[dimg[dimgidx]]++;
				dimgidx += 3;
				imgidx++;
			}
		}
	}

	void compute_dimg(double* dimg, Imgidx* dhist, Pixel* img)
	{
		Imgidx dimgidx, imgidx, stride_w = width, i, j;
		Pixel d;

		imgidx = dimgidx = 0;
		if (connectivity == 4)
		{
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					d = abs_diff(img[imgidx + stride_w], img[imgidx]);
					dimg[dimgidx] = (double)d;
					dhist[d]++;
					dimgidx++;
					d = abs_diff(img[imgidx + 1], img[imgidx]);
					dimg[dimgidx] = (double)d;
					dhist[d]++;
					dimgidx++;
					imgidx++;
				}
				d = abs_diff(img[imgidx + stride_w], img[imgidx]);
				dimg[dimgidx] = (double)d;
				dhist[d]++;
				dimgidx+=2;
				imgidx++;
			}
			for (j = 0; j < width - 1; j++)
			{
				dimgidx++;
				d = abs_diff(img[imgidx + 1], img[imgidx]);
				dimg[dimgidx] = (double)d;
				dhist[d]++;
				dimgidx++;
				imgidx++;
			}
		}
		else if (connectivity == 8)
		{
			//   -  -  3
			//   -  p  2
			//   -  0  1
			//top,middle
			for (i = 0; i < height - 1; i++)
			{
				for (j = 0; j < width - 1; j++)
				{
					d = (Pixel)(abs_diff((_int64)img[imgidx + width], (_int64)img[imgidx]));//0
					dhist[d]++;
					dimg[dimgidx++] = (double)d;
					d = (Pixel)(abs_diff((_int64)img[imgidx + width + 1], (_int64)img[imgidx]));//1
					dhist[d]++;
					dimg[dimgidx++] = (double)d;
					d = (Pixel)(abs_diff((_int64)img[imgidx + 1], (_int64)img[imgidx]));//2
					dhist[d]++;
					dimg[dimgidx++] = (double)d;
					if (i > 0)
					{
						d = (Pixel)(abs_diff((_int64)img[imgidx - width + 1], (_int64)img[imgidx]));//3
						dhist[d]++;
						dimg[dimgidx] = (double)d;
					}
					dimgidx++;
					imgidx++;
				}
				d = (Pixel)(abs_diff((_int64)img[imgidx + width], (_int64)img[imgidx]));//0
				dhist[d]++;
				dimg[dimgidx] = (double)d;
				dimgidx += 4;//skip 1,2,3
				imgidx++;
			}

			//bottom
			dimgidx += 2; //skip 0,1
			for (j = 0; j < width - 1; j++)
			{
				d = (Pixel)(abs_diff((_int64)img[imgidx + 1], (_int64)img[imgidx]));//3
				dhist[d]++;
				dimg[dimgidx] = (double)d;
				dimgidx += 4;
				imgidx++;
			}
		}
	}

	void set_isAvailable(_uint8* isAvailable)
	{
		_int32 i, j, k;
		_int32 imgsize = width * height;

		if (connectivity == 4)
		{
			//		    Neighbour Index
			// 			       3
			// 			2    pixel    1
			// 			       0
			//
			//			Neighbour indices to bit field
			//			x x x x 3 2 1 0
			//         MSB			 LSB
			//			0: Neighbour pixel not available (corner of Image, or partition in later implementation)
			//			1: available
			for (i = 0; i < imgsize; i++)
				isAvailable[i] = 0xff;

			j = width * (height - 1);
			for (i = 0; i < width; i++)
			{
				isAvailable[i] &= 0x07;
				isAvailable[j] &= 0x0e;
				j++;
			}

			j = 0;
			k = width - 1;
			for (i = 0; i < height; i++)
			{
				isAvailable[j] &= 0xb;
				isAvailable[k] &= 0xd;
				j += width;
				k += width;
			}
		}
		else
		{
			//		    Neighbour Index
			// 			5      4      3
			// 			6    pixel    2
			// 			7      0      1
			//
			//			Neighbour indices to bit field
			//			7 6 5 4 3 2 1 0
			//      MSB   			 LSB
			//			0: Neighbour pixel not available (corner of Image, or partition in later implementation)
			//			1: available

						//initialize to all available
			for (i = 0; i < imgsize; i++)
				isAvailable[i] = 0b11111111;

			//top and bottom row
			for (i = 0; i < width; i++)
				isAvailable[i] &= 0b11000111;

			for (i = width * (height - 1); i < imgsize; i++)
				isAvailable[i] &= 0b01111100;

			//leftest and rightest column
			j = 0;
			k = width - 1;
			for (i = 0; i < height; i++)
			{
				isAvailable[j] &= 0b00011111;
				isAvailable[k] &= 0b11110001;
				j += width;
				k += width;
			}
		}
	}

	void set_isAvailable(_uint8* isAvailable, int npartitions_hor, int npartitions_ver)
	{
		_int32 i, j, k;
		Imgidx imgsize = width * height;
		Imgidx wstride = width / npartitions_ver;
		Imgidx hstride = height / npartitions_hor;

		set_isAvailable(isAvailable);

		if (connectivity == 4)
		{
			//hor partitions
			j = (hstride - 1) * width;
			for (i = 0; i < npartitions_hor - 1; i++)
			{
				k = j + width;
				for (; j < k; j++)
				{
					isAvailable[j] &= 0xe;
					isAvailable[j + width] &= 0x7;
				}

				j += (hstride - 1) * width;
			}

			//ver partitions

			for (i = 0; i < npartitions_ver - 1; i++)
			{
				j = (i + 1) * wstride - 1;
				for (; j < imgsize; j += width)
				{
					isAvailable[j] &= 0xd;
					isAvailable[j + 1] &= 0xb;
				}
			}
		}
		else
		{
		}
	}

	inline _uint8 is_available(_uint8 isAvailable, _uint8 iNeighbour)
	{
		return	(isAvailable >> iNeighbour) & 1;
	}

	inline void set_field(_uint8* arr, Imgidx idx, _uint8 in)
	{
		arr[idx] = in;
	}

	inline _uint8 get_field(_uint8* arr, Imgidx idx)
	{
		return arr[idx];
	}

	inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx iNode, Pixel level)
	{
		AlphaNode<Imgidx, Pixel>* pNode;
		pNode = node + iNode;
		parentAry[pidx] = iNode;
		if (pNode->area) //possibly unnecessary branch..
			pNode->add(pix_val, pix_type);
		else
			pNode->set(1, level, (double)pix_val, pix_val, pix_val);
	}

	//#if DELAYED_NODE_ALLOC
	inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx *levelroot, _int64 level)
	{
		AlphaNode<Imgidx, Pixel>* pNode;
		Imgidx iNode = levelroot[level];
		if (iNode >= NODE_CANDIDATE)
		{
			iNode = NewAlphaNode();
			levelroot[level] = iNode;
			parentAry[pidx] = iNode;
			pNode = node + iNode;

			pNode->set(1, level, (double)pix_val, pix_val, pix_val);
		}
		else
		{
			pNode = node + iNode;
			parentAry[pidx] = iNode;
			pNode->add(pix_val, pix_type);
		}
	}
	//#else
	inline void connectPix2Node(Imgidx pidx, Pixel pix_val, Imgidx iNode)
	{
		AlphaNode<Imgidx, Pixel> *pNode = &node[iNode];
		parentAry[pidx] = iNode;
		pNode->add(pix_val, pix_type);
	}

	inline void connectPix2Node0(Imgidx pidx, Pixel pix_val, Imgidx iNode, Pixel level)
	{
		AlphaNode<Imgidx, Pixel>* pNode;
		pNode = node + iNode;
		parentAry[pidx] = iNode;
		pNode->set(1, level, (double)pix_val, pix_val, pix_val, pix_type);
	}

	inline void connectNode2Node(Imgidx prev_top, Imgidx iPar, Pixel level)
	{
		AlphaNode<Imgidx, Pixel> *pPar, *pChild;
		pChild = node + prev_top;
		pPar = node + iPar;
		pChild->parentidx = iPar;
		if (pPar->area)
			pPar->add(pChild);
		else
		{
			pPar->alpha = level;
			pPar->copy(pChild);
		}
	}

	inline void connectNode2Node(Imgidx prev_top, Imgidx iPar)
	{
		AlphaNode<Imgidx, Pixel> *pPar, *pChild;
		pChild = node + prev_top;
		pPar = node + iPar;
		pChild->parentidx = iPar;
		if (pPar->area)
			pPar->add(pChild);
		else
			pPar->copy(pChild);
	}

	inline void connectNode2Node(Imgidx* levelroot, Imgidx prev_top, _int64 level)
	{
		AlphaNode<Imgidx, Pixel> *pPar, *pChild;
		Imgidx iPar = levelroot[level];
		if (iPar >= NODE_CANDIDATE)
		{
			iPar = NewAlphaNode();
			levelroot[level] = iPar;
			pPar = node + iPar;
			pChild = node + prev_top;
			pChild->parentidx = iPar;
			pPar->alpha = level;
			pPar->copy(pChild);
		}
		else
		{
			pPar = node + iPar;
			pChild = node + prev_top;
			pChild->parentidx = iPar;
			pPar->add(pChild);
		}
	}

	inline Imgidx NewAlphaNode()
	{
		if (curSize == maxSize)
		{
			std::cout << "Reallocating...\n";
			maxSize = _min((1 + (connectivity >> 1)) * height * width, maxSize + (Imgidx)(2 * height * width * 0.1));

			node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		}
		return curSize++;
	}

	inline Imgidx NewAlphaNode(Pixel level, AlphaNode<Imgidx, Pixel> *pCopy)
	{
		AlphaNode<Imgidx, Pixel> *pNew = node + curSize;

		if (curSize == maxSize)
		{
			std::cout << "Reallocating...\n";
			maxSize = _min((1 + (connectivity >> 1)) * height * width, maxSize + (Imgidx)(2 * height * width * 0.1));

			node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
			pNew = node + curSize;
		}
		pNew->alpha = level;
		pNew->copy(pCopy);
		return curSize++;
	}

	inline Imgidx NewAlphaNode1(double level, AlphaNode<Imgidx, Pixel> *pCopy)
	{
		AlphaNode<Imgidx, Pixel> *pNew = node + curSize;

		if (curSize == maxSize)
		{
			std::cout << "Reallocating...\n";
			maxSize = _min((1 + (connectivity >> 1)) * height * width, maxSize + (Imgidx)(2 * height * width * 0.1));

			node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
			pNew = node + curSize;
		}
		pNew->alpha = level;
		pNew->copy(pCopy);
		return curSize++;
	}

	//#else
	inline Imgidx NewAlphaNode(Pixel level) //Fix it later - no need to initialize
	{
		AlphaNode<Imgidx, Pixel> *pNew = node + curSize;

		if (curSize == maxSize)
		{
			std::cout << "Reallocating...\n";
			maxSize = _min((1 + (connectivity >> 1)) * height * width, maxSize + (Imgidx)(height * width * 0.1));

			node = (AlphaNode<Imgidx, Pixel>*)Realloc(node, maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
			pNew = node + curSize;
		}
		pNew->alpha = level;
		pNew->minPix = (_uint8)-1;
		pNew->maxPix = 0;
		pNew->sumPix = 0.0;
		//		pNew->parentidx = 0;
		pNew->area = 0;

		return curSize++;
	}
	//#endif
	inline _uint8 is_visited(_uint8* isVisited, Imgidx p)
	{
		return isVisited[p];
	}
	inline void visit(_uint8* isVisited, Imgidx p)
	{
		isVisited[p] = 1;
	}

#define _PUSH_NEIGHBORS_4N \
	{q = p << 1;\
	(is_available(isAv, 0) && !isVisited[p + width]) ? (void)queue->push(p + width, dimg[q]) : (void)0;\
	(is_available(isAv, 1) && !isVisited[p + 1]) 		 ? (void)queue->push(p + 1, dimg[q + 1]) : (void)0;\
	(is_available(isAv, 2) && !isVisited[p - 1]) 		 ? (void)queue->push(p - 1, dimg[q - 1]) : (void)0;\
	(is_available(isAv, 3) && !isVisited[p - width]) ? (void)queue->push(p - width, dimg[q - (width << 1)]) : (void)0;}

#define _PUSH_NEIGHBORS_8N \
	{q = p << 2;\
	(is_available(isAv, 0) && !isVisited[p + wstride1])?	(void)queue->push(p + wstride1, dimg[q]) 							   : (void)0;\
	(is_available(isAv, 1) && !isVisited[p + width])   ?	(void)queue->push(p + width, dimg[q + 1]) 						   : (void)0;\
	(is_available(isAv, 2) && !isVisited[p + wstride0])?	(void)queue->push(p + wstride0, dimg[q + 2]) 					   : (void)0;\
	(is_available(isAv, 3) && !isVisited[p + 1])	     ?	(void)queue->push(p + 1, dimg[q + 3]) 									 : (void)0;\
	(is_available(isAv, 4) && !isVisited[p - wstride1])?	(void)queue->push(p - wstride1, dimg[q - wstride_d + 4]) : (void)0;\
	(is_available(isAv, 5) && !isVisited[p - width])   ?	(void)queue->push(p - width, dimg[q - wstride_d + 1]) 	 : (void)0;\
	(is_available(isAv, 6) && !isVisited[p - wstride0])?	(void)queue->push(p - wstride0, dimg[q - wstride_d - 2]) : (void)0;\
	(is_available(isAv, 7) && !isVisited[p - 1])	     ?  (void)queue->push(p - 1, dimg[q - 1]) 									 : (void)0;}

#define _CREATE_NEW_NODE(minlev) \
	{Pixel pix_val = img[p];\
	current_level = minlev;\
	iNode = NewAlphaNode();\
	node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);\
	node[iNode].parentidx = stack_top;\
	node[iNode].rootidx = ROOTIDX;\
	stack_top = iNode;}

#define _CREATE_NEW_STACKTOP(minlev) \
	{iNode = NewAlphaNode(minlev, node + stack_top);\
	node[iNode].parentidx = node[stack_top].parentidx;\
	node[iNode].rootidx = ROOTIDX;\
	node[stack_top].parentidx = iNode;}

#define _CREATE_SINGLETON_NODE_0 \
	{iNode = NewAlphaNode(0, node + stack_top);\
	node[iNode].parentidx = stack_top;\
	node[iNode].rootidx = ROOTIDX;\
	parentAry[p] = iNode;\
	prev_top = iNode;}\

#define _CREATE_SINGLETON_NODE_1 \
	{Pixel pix_val = img[p];\
	iNode = NewAlphaNode();\
	node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);\
	node[stack_top].add(node + iNode, pix_type);\
	node[iNode].parentidx = stack_top;\
	node[iNode].rootidx = ROOTIDX;\
	parentAry[p] = iNode;}


	inline Imgidx TreeSizeEstimation(Imgidx *dhist, _int64 numlevels, Imgidx imgsize, Imgidx nredges)
	{
		return TreeSizeEstimation(dhist, numlevels, imgsize, nredges, M);
	}

	Imgidx TreeSizeEstimation(Imgidx *dhist, _int64 numlevels, Imgidx imgsize, Imgidx nredges, double m)
	{
		if (imgsize < TSE_MINSIZE)
			return 10 + nredges + imgsize;
		double tse_nrmsd = 0;
		for (_int64 p = 0; p < numlevels; p++)
			tse_nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
		tse_nrmsd = sqrt((tse_nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
		nrmsd = tse_nrmsd;
		Imgidx ret = _min(2 * imgsize, (Imgidx)(2 * imgsize * ((A * exp(SIGMA * tse_nrmsd) + B) + m)));
		double dret = _min((double)(2 * imgsize), (double)(2 * imgsize * ((A * exp(SIGMA * tse_nrmsd) + B) + m)));

		if (ret < 0)
		{
			printf("Warning: TSE returned 0< value\n");
			printf("nrmsd = %lf\n",tse_nrmsd);
			printf(" exp(SIGMA * nrmsd) = %lf\n", exp(SIGMA * tse_nrmsd));
			printf("A * exp(SIGMA * nrmsd) + B = %lf\n",A * exp(SIGMA * tse_nrmsd) + B);
			printf("2 * imgsize * ((A * exp(SIGMA * nrmsd) + B) + m) = %lf\n",2 * imgsize * ((A * exp(SIGMA * tse_nrmsd) + B) + m));
			printf("(Imgidx)(2 * imgsize * ((A * exp(SIGMA * nrmsd) + B) + m)) = %d\n", (int)(Imgidx)(2 * imgsize * ((A * exp(SIGMA * tse_nrmsd) + B) + m)));
			printf("nredges = %lf\n",(double)nredges);
			printf("imgsize = %lf\n",(double)imgsize);
			printf("1 + nredges + imgsize = %lf\n",(double)(1 + nredges + imgsize));
			printf("ret = %lf\n", (double)ret);
			printf("ret = %lf\n", (double)dret);
		}

		return ret;
	}

	Imgidx TreeSizeEstimation(Imgidx *dhist, _int64 numlevels, Imgidx imgsize, Imgidx nredges, double m, Imgidx reserve)
	{
		nrmsd = 0;
		for (_int64 p = 0; p < numlevels; p++)
			nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
		nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
		return _min(2.0 * imgsize, (Imgidx)(2.0 * (double)imgsize * ((A * exp(SIGMA * nrmsd) + B) + m))) + reserve;
	}

#define _SET_COMMON_MEMORY	\
	if (dhist) Free(dhist);\
	isVisited = (_uint8*)Calloc((size_t)((imgsize)));\
	isAvailable = (_uint8*)Malloc((size_t)(imgsize));\
	set_isAvailable(isAvailable);\
	parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(_int32));\
	node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

	inline void remove_redundant_node(Imgidx& prev_top, Imgidx& stack_top)
	{
		if (node[prev_top].parentidx == stack_top && node[prev_top].area == node[stack_top].area)
		{
			node[prev_top].parentidx = node[stack_top].parentidx;
			stack_top = prev_top;
			curSize--;
		}
	}

	template <class Queue>
	void Flood_HierarQueue(Pixel* img, Queue* queue, int tse)
	{
		if (sizeof(Pixel) > 2 || channel > 1)
		{
			printf("Error: Hierarchical queues do not work on >16 bits images or multispectral images\n");
			printf("Try Unionfind (algorithm code %d), flooding using Heapqueue (%d), trie queue (%d) or cached trie queue (%d) \n", UNIONFIND, FLOOD_HEAPQUEUE_CACHE, FLOOD_TRIE, FLOOD_TRIE_CACHE);
			return;
		}

		Imgidx imgsize, dimgsize, nredges, x0;
		_uint64 numlevels, max_level, current_level;
		Imgidx *dhist;
		Imgidx prev_top, stack_top, iNode;
		_uint8 *isVisited, *isAvailable, isAv;
		Imgidx p, q;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		max_level = (sizeof(Pixel)==8) ? 0xffffffffffffffff : (_int64)((Pixel)(-1));
		numlevels = max_level + 1;

		Pixel *dimg;
		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		memset(dhist, 0, (size_t)numlevels * sizeof(_int32));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		compute_dimg(dimg, dhist, img);//calculate pixel differences and make histogram

		//create hierarchical queue from dhist
		queue = new Queue(nredges + 1, dhist, numlevels); // +1 for the dummy node
		curSize = 0;

		if (!tse || imgsize < 10000 || sizeof(Pixel) > 1) //for small imags do not use TSE
			maxSize = 1 + imgsize + nredges;
		else
			maxSize = TreeSizeEstimation(dhist, numlevels, imgsize, nredges);

		if (dhist) Free(dhist);

		isVisited = (_uint8*)Calloc((size_t)((imgsize)));
		isAvailable = (_uint8*)Malloc((size_t)(imgsize));
		set_isAvailable(isAvailable);

		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(_int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

		stack_top = NewAlphaNode();/*dummy root*/
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
		current_level = max_level;
		x0 = 0; /*arbitrary starting point*/
		prev_top = stack_top;

		queue->push(x0, current_level);
		while (1) //flooding
		{
			while ((_uint64)queue->min_level <= current_level) //flood all levels below current_level
			{
				p = queue->pop();
				if (isVisited[p])
				{
					queue->find_minlev();
					continue;
				}
				isVisited[p] = 1;
				isAv = isAvailable[p];
				if (connectivity == 4)
				{
					q = p << 1;
					if (is_available(isAv, 0) && !isVisited[p + width]) queue->push(p + width, dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1]) 		 queue->push(p + 1, dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1]) 		 queue->push(p - 1, dimg[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width]) queue->push(p - width, dimg[q - (width << 1)]);
				}
				else if (connectivity == 8)
				{
					Imgidx width4 = width << 2;
					q = p << 2;
					if (is_available(isAv, 0) && !isVisited[p + width])				queue->push(p + width, 		dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + width + 1])		queue->push(p + width + 1, dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p + 1])						queue->push(p + 1, 				dimg[q + 2]);
					if (is_available(isAv, 3) && !isVisited[p - width + 1])		queue->push(p - width + 1, dimg[q + 3]);
					if (is_available(isAv, 4) && !isVisited[p - width])				queue->push(p - width, 		dimg[q - width4]);
					if (is_available(isAv, 5) && !isVisited[p - width - 1])		queue->push(p - width - 1, dimg[q - width4 - 3]);
					if (is_available(isAv, 6) && !isVisited[p - 1])						queue->push(p - 1, 				dimg[q - 2]);
					if (is_available(isAv, 7) && !isVisited[p + width - 1])		queue->push(p + width - 1, dimg[q + width4 - 1]);
				}
				else
				{
					//?
				}
				//else //later

				if (current_level > (_uint64)queue->min_level) //go to lower level
				{

					{//creat new node
						Pixel pix_val = img[p];
						current_level = queue->min_level;
						iNode = NewAlphaNode();
						node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
						node[iNode].parentidx = stack_top;
						node[iNode].rootidx = ROOTIDX;
						stack_top = iNode;
					}
					if (current_level)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						node[iNode].rootidx = ROOTIDX;
						parentAry[p] = iNode;
						prev_top = iNode;
					}
					else
						parentAry[p] = stack_top;
				}
				else
				{
					queue->find_minlev();

					if (current_level)
					{
						Pixel pix_val = img[p];
						iNode = NewAlphaNode();
						node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
						node[stack_top].add(node + iNode, pix_type);
						node[iNode].parentidx = stack_top;
						node[iNode].rootidx = ROOTIDX;
						parentAry[p] = iNode;
					}
					else
						connectPix2Node(p, img[p], stack_top);
				}
			}

			remove_redundant_node(prev_top, stack_top);

			//go to higher level
			iNode = node[stack_top].parentidx;
			if ((Pixel)queue->min_level < (Pixel)node[iNode].alpha) //new level from queue
			{
				iNode = NewAlphaNode(queue->min_level, node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[iNode].rootidx = ROOTIDX;
				node[stack_top].parentidx = iNode;
			}
			else //go to existing node
			{
				if (node[iNode].area == imgsize)	// root node found...done
					break;
				node[iNode].add(node + stack_top, pix_type);
			}

			if (node[iNode].area == imgsize)	// root node found...done
				break;

			prev_top = stack_top;
			stack_top = iNode;
			current_level = (_uint64)node[stack_top].alpha;
		}
		rootidx = (node[stack_top].area == imgsize) ? stack_top : iNode; //remove redundant root
		node[rootidx].parentidx = ROOTIDX;

		delete queue;
		Free(dimg);
		Free(isVisited);
		Free(isAvailable);
	}

	void Flood_HeapQueue(Pixel* img)
	{
		HeapQueue<Imgidx, double>* queue;

		Imgidx imgsize, dimgsize, nredges, x0;
		_uint64 numlevels, max_level;
		Imgidx *dhist;
		Imgidx stack_top, prev_top, iNode;
		_uint8 *isVisited, *isAvailable, isAv;
		Imgidx p, q;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;\
		max_level = (Pixel)(-1);
		numlevels = max_level + 1;
		double current_level;
		double *dimg;

		dimg = (double*)Malloc((size_t)dimgsize * sizeof(double));

		if (sizeof(Pixel) == 1 && channel == 1)
		{
			dhist = (Imgidx*)Calloc(numlevels * sizeof(Imgidx));
			compute_dimg(dimg, dhist, img);//calculate pixel differences and make histogram
			maxSize = TreeSizeEstimation(dhist, numlevels, imgsize, nredges);
		}
		else
		{
			dhist = 0;
			compute_dimg(dimg, img);
			maxSize = 1 + imgsize + nredges;
		}

		//create heap-based priority queue
		queue = new HeapQueue<Imgidx, double>(nredges);
		curSize = 0;

		_SET_COMMON_MEMORY

		x0 = 0; /*arbitrary starting point*/
		current_level = DBL_MAX;
		stack_top = NewAlphaNode();/*dummy root*/
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, current_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
		//current_level = max_level;
		prev_top = stack_top; /*to find redundant node*/

		queue->push_run(x0, current_level);
		while (1)
		{
			while (queue->get_minlev() <= current_level) //flood all levels below current_level
			{
				p = queue->top();
				if (isVisited[p])
				{
					queue->pop();
					continue;
				}
				isVisited[p] = 1;

				isAv = isAvailable[p];
				if (connectivity == 4)
				{
					q = p << 1;
					if (is_available(isAv, 0) && !isVisited[p + width])	queue->push(p + width, dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1])			queue->push(p + 1, dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1])			queue->push(p - 1, dimg[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])	queue->push(p - width, dimg[q - (width << 1)]);
				}
				else if (connectivity == 8)
				{
					Imgidx width4 = width << 2;
					q = p << 2;
					if (is_available(isAv, 0) && !isVisited[p + width]) 		 queue->push(p + width, 		dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + width + 1]) queue->push(p + width + 1, dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p + 1]) 		 		 queue->push(p + 1, 				dimg[q + 2]);
					if (is_available(isAv, 3) && !isVisited[p - width + 1]) queue->push(p - width + 1, dimg[q + 3]);
					if (is_available(isAv, 4) && !isVisited[p - width]) 		 queue->push(p - width, 		dimg[q - width4]);
					if (is_available(isAv, 5) && !isVisited[p - width - 1]) queue->push(p - width - 1, dimg[q - width4 - 3]);
					if (is_available(isAv, 6) && !isVisited[p - 1]) 				 queue->push(p - 1, 				dimg[q - 2]);
					if (is_available(isAv, 7) && !isVisited[p + width - 1]) queue->push(p + width - 1, dimg[q + width4 - 1]);
				}
				else
				{
					//?
				}

				queue->find_minlev();
				if (current_level > queue->get_minlev()) 
				{
						_CREATE_NEW_NODE(queue->get_minlev())
						if (current_level)
						{
							iNode = NewAlphaNode(0, node + stack_top);
							node[iNode].parentidx = stack_top;
							node[iNode].rootidx = ROOTIDX;
							parentAry[p] = iNode;
							prev_top = iNode;
						}
						else
							parentAry[p] = stack_top;
				}
				else
				{
					if (current_level)
						_CREATE_SINGLETON_NODE_1
					else
						connectPix2Node(p, img[p], stack_top);
					if (node[stack_top].area == imgsize)
						goto FLOOD_END;
				}
			}
			
			remove_redundant_node(prev_top, stack_top);

			if (node[stack_top].area == imgsize)	// root node found...done
				break;

			iNode = node[stack_top].parentidx;
			if (queue->get_minlev() < node[iNode].alpha)
			{
				iNode = NewAlphaNode1(queue->get_minlev(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[iNode].rootidx = ROOTIDX;
				node[stack_top].parentidx = iNode;
			}
			else
				node[iNode].add(node + stack_top, pix_type);

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
			if (node[stack_top].area == imgsize)	// root node found...done
				break;
		}
	FLOOD_END:
		node[stack_top].parentidx = ROOTIDX;

		delete queue;
		Free(dimg);
		Free(isVisited);
		Free(isAvailable);
	}

	void Flood_HeapQueue_Cache(Pixel* img)
	{
		Cache_Quad_Heapqueue<Imgidx, double>* queue;

		Imgidx imgsize, dimgsize, nredges, x0;
		_uint64 numlevels, max_level;
		Imgidx *dhist;
		Imgidx stack_top, prev_top, iNode;
		_uint8 *isVisited, *isAvailable, isAv;
		Imgidx p, q;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;\
		max_level = (Pixel)(-1);
		numlevels = max_level + 1;
		double current_level;
		double *dimg;

		dimg = (double*)Malloc((size_t)dimgsize * sizeof(double));

		if (sizeof(Pixel) == 1 && channel == 1)
		{
			dhist = (Imgidx*)Calloc(numlevels * sizeof(Imgidx));
			compute_dimg(dimg, dhist, img);//calculate pixel differences and make histogram
			maxSize = TreeSizeEstimation(dhist, numlevels, imgsize, nredges);
		}
		else
		{
			dhist = 0;
			compute_dimg(dimg, img);
			maxSize = 1 + imgsize + nredges;
		}

		//create heap-based priority queue
		queue = new Cache_Quad_Heapqueue<Imgidx, double>(nredges);
		curSize = 0;

		_SET_COMMON_MEMORY

		x0 = 0; /*arbitrary starting point*/
		current_level = DBL_MAX;
		stack_top = NewAlphaNode();/*dummy root*/
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, current_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
		prev_top = stack_top; /*to find redundant node*/


		queue->push_1stitem(x0, current_level);
		while (1)
		{
			while (queue->get_minlev() <= current_level) //flood all levels below current_level
			{
				p = queue->top();
				if (isVisited[p])
				{
					queue->pop();
					continue;
				}
				queue->start_pushes();
				isVisited[p] = 1;

				isAv = isAvailable[p];
				if (connectivity == 4)
				{
					q = p << 1;
					if (is_available(isAv, 0) && !isVisited[p + width])	queue->push(p + width, dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1])		queue->push(p + 1, dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1])		queue->push(p - 1, dimg[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])	queue->push(p - width, dimg[q - (width << 1)]);
				}
				else if (connectivity == 8)
				{
					Imgidx width4 = width << 2;
					q = p << 2;
					if (is_available(isAv, 0) && !isVisited[p + width])		queue->push(p + width, 		dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + width + 1])	queue->push(p + width + 1,	dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p + 1]) 		 	queue->push(p + 1, 			dimg[q + 2]);
					if (is_available(isAv, 3) && !isVisited[p - width + 1])	queue->push(p - width + 1,	dimg[q + 3]);
					if (is_available(isAv, 4) && !isVisited[p - width]) 		queue->push(p - width, 		dimg[q - width4]);
					if (is_available(isAv, 5) && !isVisited[p - width - 1])	queue->push(p - width - 1,	dimg[q - width4 - 3]);
					if (is_available(isAv, 6) && !isVisited[p - 1]) 			queue->push(p - 1, 			dimg[q - 2]);
					if (is_available(isAv, 7) && !isVisited[p + width - 1])	queue->push(p + width - 1,	dimg[q + width4 - 1]);
				}
				else
				{
					//?
				}


				queue->end_pushes();
				if (current_level > queue->get_minlev()) //remove typecasting later
				{
						_CREATE_NEW_NODE(queue->get_minlev())
						if (current_level)
						{
							iNode = NewAlphaNode(0, node + stack_top);
							node[iNode].parentidx = stack_top;
							node[iNode].rootidx = ROOTIDX;
							parentAry[p] = iNode;
							prev_top = iNode;
						}
						else
							parentAry[p] = stack_top;
				}
				else
				{
					if (current_level)
						_CREATE_SINGLETON_NODE_1
					else
						connectPix2Node(p, img[p], stack_top);
					if (node[stack_top].area == imgsize)
						goto FLOOD_END;
				}
			}
			remove_redundant_node(prev_top, stack_top);

			if (node[stack_top].area == imgsize)	// root node found...done
				break;

			iNode = node[stack_top].parentidx;
			if (queue->get_minlev() < node[iNode].alpha)
			{
				iNode = NewAlphaNode1(queue->get_minlev(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[iNode].rootidx = ROOTIDX;
				node[stack_top].parentidx = iNode;
			}
			else
				node[iNode].add(node + stack_top, pix_type);

			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
			if (node[stack_top].area == imgsize)	// root node found...done
				break;
		}
	FLOOD_END:
		node[stack_top].parentidx = ROOTIDX;

		delete queue;
		Free(dimg);
		Free(isVisited);
		Free(isAvailable);
	}

	void Flood_HierarQueue_Cache(Pixel* img)
	{
		HierarQueueCache<Imgidx, Pixel>* queue;
		Imgidx imgsize, dimgsize, nredges, x0;
		_uint64 numlevels, max_level, current_level;
		Imgidx *dhist;
		Imgidx prev_top, stack_top, iNode;
		_uint8 *isVisited, *isAvailable, isAv;
		Imgidx p, q;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		max_level = (sizeof(Pixel)==8) ? 0xffffffffffffffff : (_int64)((Pixel)(-1));
		numlevels = max_level + 1;

		Pixel *dimg;
		dhist = (Imgidx*)Calloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));

		max_level = compute_dimg1(dimg, dhist, img);//calculate pixel differences and make histogram
		numlevels = max_level + 1;

		//create hierarchical queue from dhist
		queue = new HierarQueueCache<Imgidx, Pixel>(nredges + 1, dhist, numlevels); // +1 for the dummy node
		curSize = 0;

		if (imgsize < 10000 || sizeof(Pixel) > 2) //for small imags do not use TSE
			maxSize = 1 + imgsize + dimgsize;
		else
	    	maxSize = TreeSizeEstimation(dhist, numlevels, imgsize, nredges);
		if (dhist)
			Free(dhist);

		isVisited = (_uint8*)Calloc((size_t)((imgsize)));
		isAvailable = (_uint8*)Malloc((size_t)(imgsize));
		set_isAvailable(isAvailable);
		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(_int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

		stack_top = NewAlphaNode();/*dummy root*/
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
		current_level = max_level;
		x0 = 0; /*arbitrary starting point*/
		prev_top = stack_top;

		queue->push_1stitem(x0, current_level);
		while (1) //flooding
		{
			while ((_int32)queue->top_alpha() <= (_int32)current_level) //flood all levels below current_level
			{
				p = queue->top();
				if (isVisited[p] == 1)
				{
					queue->pop();
					continue;
				}

				queue->start_pushes();
				isVisited[p] = 1;
				isAv = isAvailable[p];
				if (connectivity == 4)
				{
					q = p << 1;
					if (is_available(isAv, 0) && !isVisited[p + width])	queue->push(p + width, dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1]) 		queue->push(p + 1, dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1]) 		queue->push(p - 1, dimg[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])	queue->push(p - width, dimg[q - (width << 1)]);
				}
				else if (connectivity == 8)
				{
					Imgidx width4 = width << 2;
					q = p << 2;
					if (is_available(isAv, 0) && !isVisited[p + width]) 		queue->push(p + width, 		dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + width + 1])	queue->push(p + width + 1, dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p + 1]) 		 	queue->push(p + 1, 				dimg[q + 2]);
					if (is_available(isAv, 3) && !isVisited[p - width + 1])	queue->push(p - width + 1, dimg[q + 3]);
					if (is_available(isAv, 4) && !isVisited[p - width]) 		queue->push(p - width, 		dimg[q - width4]);
					if (is_available(isAv, 5) && !isVisited[p - width - 1])	queue->push(p - width - 1, dimg[q - width4 - 3]);
					if (is_available(isAv, 6) && !isVisited[p - 1]) 			queue->push(p - 1, 				dimg[q - 2]);
					if (is_available(isAv, 7) && !isVisited[p + width - 1])	queue->push(p + width - 1, dimg[q + width4 - 1]);
				}
				else
				{
					//?
				}

				queue->end_pushes();
				if (current_level > (_uint64)queue->top_alpha()) //go to lower level
				{
					_CREATE_NEW_NODE(queue->top_alpha())
					if (current_level)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						node[iNode].rootidx = ROOTIDX;
						parentAry[p] = iNode;
						prev_top = iNode;
					}
					else
					parentAry[p] = stack_top;
				}
				else
				{
					if (current_level)
					_CREATE_SINGLETON_NODE_1
					else
					connectPix2Node(p, img[p], stack_top);
					if (node[stack_top].area == imgsize)
						goto FLOOD_END;
				}
			}

			remove_redundant_node(prev_top, stack_top);

			if (node[stack_top].area == imgsize)	// root node found...done
				break;

			//go to higher level
			iNode = node[stack_top].parentidx;
			if ((_int32)queue->top_alpha() < (_int32)node[iNode].alpha)
			{
				iNode = NewAlphaNode1(queue->top_alpha(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[iNode].rootidx = ROOTIDX;
				node[stack_top].parentidx = iNode;
			}
			else //go to existing node
			{
				node[iNode].add(node + stack_top, pix_type);
			}

			prev_top = stack_top;
			stack_top = iNode;
			current_level = (_int32)node[stack_top].alpha;
			if (node[stack_top].area == imgsize)	// root node found...done
				break;
		}

FLOOD_END:
		rootidx = (node[stack_top].area == imgsize) ? stack_top : iNode; //remove redundant root
		node[rootidx].parentidx = ROOTIDX;

		delete queue;
		Free(dimg);
		Free(isVisited);
		Free(isAvailable);
	}

	int get_bitdepth(_uint64 num)
	{
		int ret = 0;
		while(num)
		{
			ret++;
			num >>= 1;
		}
		return ret;
	}

	void Flood_HierarHeapQueue(Pixel* img, double a = 12.0, double r = 0.5, int listsize = 12)
	{
	  	HierarHeapQueue<Imgidx, Pixel>* queue;

		Imgidx imgsize, dimgsize, nredges, x0;
		_uint64 numlevels, max_level;
		double current_level;
		Imgidx *dhist;
		Imgidx stack_top, iNode;
		_uint8 *isVisited, *isAvailable, isAv;
		Imgidx p, q;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (1 + (connectivity >> 1)) * width * height;

		max_level = (sizeof(Pixel)==8) ? 0xffffffffffffffff : (_int64)((Pixel)(-1));
		numlevels = (_uint64)(a * log2(1 + (double)max_level)) + 1;

		Pixel *dimg;
		dhist = (Imgidx*)Calloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));

		compute_dimg(x0, current_level, dimg, dhist, img, a);//calculate pixel differences and make histogram

		//create hierarchical queue from dhist
		queue = new HierarHeapQueue<Imgidx, Pixel>(dhist, numlevels, nredges, a, listsize, connectivity, r); // +1 for the dummy node
		curSize = 0;
		
		maxSize = 1 + imgsize + dimgsize; //Do not use TSE here, becasue dhist is a logged histogram (also this algorithm is for hdr)


		dhist = 0;
		isVisited = (_uint8*)Calloc((size_t)((imgsize)));
		isAvailable = (_uint8*)Malloc((size_t)(imgsize));
		set_isAvailable(isAvailable);
		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(_int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

		stack_top = NewAlphaNode();/*dummy root*/
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (double)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
		current_level = (double)max_level;

		bool firstpixel = 1;
		p = x0;
		Imgidx prev_top = stack_top;
		while (1) //flooding
		{
			while (firstpixel || (double)queue->top_alpha() <= (double)current_level) //flood all levels below current_level
			{
				if (!firstpixel)
				{
					p = queue->pop(isVisited);
					if (isVisited[p])
						continue;
				}
				else
					firstpixel = 0;
				
				isVisited[p] = 1;
				isAv = isAvailable[p];
				if (connectivity == 4)
				{
					q = p << 1;
					if (is_available(isAv, 0) && !isVisited[p + width])	queue->push(p + width, dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1])		queue->push(p + 1, dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1])		queue->push(p - 1, dimg[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])	queue->push(p - width, dimg[q - (width << 1)]);
				}
				else if (connectivity == 8)
				{
					Imgidx width4 = width << 2;
					q = p << 2;
					if (is_available(isAv, 0) && !isVisited[p + width])		{queue->push(p + width,		dimg[q]); 				};
					if (is_available(isAv, 1) && !isVisited[p + width + 1])	{queue->push(p + width + 1,	dimg[q + 1]); 			};
					if (is_available(isAv, 2) && !isVisited[p + 1])			{queue->push(p + 1,			dimg[q + 2]); 			};
					if (is_available(isAv, 3) && !isVisited[p - width + 1])	{queue->push(p - width + 1, dimg[q + 3]);			};
					if (is_available(isAv, 4) && !isVisited[p - width])		{queue->push(p - width,		dimg[q - width4]); 		};
					if (is_available(isAv, 5) && !isVisited[p - width - 1])	{queue->push(p - width - 1, dimg[q - width4 - 3]);	};
					if (is_available(isAv, 6) && !isVisited[p - 1])			{queue->push(p - 1,			dimg[q - 2]); 			};
					if (is_available(isAv, 7) && !isVisited[p + width - 1])	{queue->push(p + width - 1, dimg[q + width4 - 1]);	};
				}
				else
				{
					//?
				}
				
				if ((double)current_level > (double)queue->top_alpha()) //go to lower level
				{
					Pixel pix_val = img[p];
					current_level = queue->top_alpha();
					iNode = NewAlphaNode();
					node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
					node[iNode].parentidx = stack_top;
					node[iNode].rootidx = ROOTIDX;
					prev_top = stack_top;
					stack_top = iNode;
					if (current_level > 0)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						node[iNode].rootidx = ROOTIDX;
						parentAry[p] = iNode;
					}
					else
					parentAry[p] = stack_top;
				}
				else
				{
					if (current_level > 0)
						{
							Pixel pix_val = img[p];
							iNode = NewAlphaNode();
							node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
							node[stack_top].add(node + iNode, pix_type);
							node[iNode].parentidx = stack_top;
							node[iNode].rootidx = ROOTIDX;
							parentAry[p] = iNode;
						}
					else
					connectPix2Node(p, img[p], stack_top);
					if (node[stack_top].area == imgsize)
						goto FLOOD_END;
				}
	    	}

			if (node[prev_top].parentidx == stack_top && node[prev_top].area == node[stack_top].area)
			{
				node[prev_top].parentidx = node[stack_top].parentidx;
				stack_top = prev_top;
				curSize--;
			}

			if (node[stack_top].area == imgsize)	// root node found...done
				break;

			//go to higher level
			iNode = node[stack_top].parentidx;
			if ((double)queue->top_alpha() < (double)node[iNode].alpha)
			{
				iNode = NewAlphaNode1(queue->top_alpha(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[iNode].rootidx = ROOTIDX;
				node[stack_top].parentidx = iNode;
			}
			else //go to existing node
			{
			node[iNode].add(node + stack_top, pix_type);
			}
			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
			if (node[stack_top].area == imgsize)	// root node found...done
				break;
	  }
FLOOD_END:
	  rootidx = (node[stack_top].area == imgsize) ? stack_top : iNode; //remove redundant root
	  node[rootidx].parentidx = ROOTIDX;

	  delete queue;
	  Free(dimg);
	  Free(isVisited);
	  Free(isAvailable);
	}

	void Flood_HierarHeapQueue_Cache(Pixel* img, double a = 12.0, double r = 0.5, int listsize = 12)
	{
	  	HierarHeapQueue_cache<Imgidx, Pixel>* queue;

		Imgidx imgsize, dimgsize, nredges, x0;
		_uint64 numlevels, max_level;
		double current_level;
		Imgidx *dhist;
		Imgidx stack_top, iNode;
		_uint8 *isVisited, *isAvailable, isAv;
		Imgidx p, q;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (1 + (connectivity >> 1)) * width * height;

		//double a = 4.0;
		max_level = (sizeof(Pixel)==8) ? 0xffffffffffffffff : (_int64)((Pixel)(-1));
		numlevels = (_uint64)(a * log2(1 + (double)max_level)) + 1;

		Pixel *dimg;
		dhist = (Imgidx*)Calloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));

		compute_dimg(dimg, dhist, img, a);//calculate pixel differences and make histogram

		//create hierarchical queue from dhist
		queue = new HierarHeapQueue_cache<Imgidx, Pixel>(dhist, numlevels, nredges, a, listsize, connectivity, r); // +1 for the dummy node
		curSize = 0;
		maxSize = 1 + imgsize + dimgsize; //Do not use TSE here, becasue dhist is a logged histogram (also this algorithm is for hdr)

		dhist = 0;
		isVisited = (_uint8*)Calloc((size_t)((imgsize)));
		isAvailable = (_uint8*)Malloc((size_t)(imgsize));
		set_isAvailable(isAvailable);
		parentAry = (Imgidx*)Malloc((size_t)imgsize * sizeof(_int32));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));

		stack_top = NewAlphaNode();/*dummy root*/
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (double)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
		current_level = (double)max_level;
		x0 = 0; /*arbitrary starting point*/

		Imgidx prev_top = stack_top;
		queue->push_1stitem(x0, (Pixel)current_level);
		while (1) //flooding
		{
			while ((double)queue->top_alpha() <= (double)current_level) //flood all levels below current_level
			{
				p = queue->top();

				if (isVisited[p])
				{
					queue->pop(isVisited);
							//queue->mark(1);
					continue;
				}
				queue->start_pushes();
				isVisited[p] = 1;

				isAv = isAvailable[p];
				if (connectivity == 4)
				{
					q = p << 1;
					if (is_available(isAv, 0) && !isVisited[p + width])	queue->push(p + width, dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1])			queue->push(p + 1, dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1])			queue->push(p - 1, dimg[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])	queue->push(p - width, dimg[q - (width << 1)]);
				}
				else if (connectivity == 8)
				{
					Imgidx width4 = width << 2;
					q = p << 2;
					if (is_available(isAv, 0) && !isVisited[p + width]) 		 {queue->push(p + width, 		dimg[q]); 							};//printf("0:pushing %d at %.3f \n",(int)(p + width),log2((double)dimg[q] + 1));}
					if (is_available(isAv, 1) && !isVisited[p + width + 1]) {queue->push(p + width + 1, dimg[q + 1]); 					};//printf("1:pushing %d at %.3f \n",(int)(p + width + 1),log2((double)dimg[q + 1] + 1));}
					if (is_available(isAv, 2) && !isVisited[p + 1]) 		 		 {queue->push(p + 1, 				dimg[q + 2]); 					};//printf("2:pushing %d at %.3f \n",(int)(p + 1),log2((double)dimg[q + 2] + 1));}
					if (is_available(isAv, 3) && !isVisited[p - width + 1]) {queue->push(p - width + 1, dimg[q + 3]); 					};//printf("3:pushing %d at %.3f \n",(int)(p - width + 1),log2((double)dimg[q + 3] + 1));}
					if (is_available(isAv, 4) && !isVisited[p - width]) 		 {queue->push(p - width, 		dimg[q - width4]); 			};//printf("4:pushing %d at %.3f \n",(int)(p - width),log2((double)dimg[q - width4] + 1));}
					if (is_available(isAv, 5) && !isVisited[p - width - 1]) {queue->push(p - width - 1, dimg[q - width4 - 3]); };//printf("5:pushing %d at %.3f \n",(int)(p - width - 1),log2((double)dimg[q - width4 - 3] + 1));}
					if (is_available(isAv, 6) && !isVisited[p - 1]) 				 {queue->push(p - 1, 				dimg[q - 2]); 					};//printf("6:pushing %d at %.3f \n",(int)(p - 1),log2((double)dimg[q - 2] + 1));}
					if (is_available(isAv, 7) && !isVisited[p + width - 1]) {queue->push(p + width - 1, dimg[q + width4 - 1]); };//printf("7:pushing %d at %.3f \n",(int)(p + width - 1),log2((double)dimg[q + width4 - 1] + 1));}
				}
				else
				{
					//?
				}
				queue->end_pushes(isVisited);
				if ((double)current_level > (double)queue->top_alpha()) //go to lower level
				{
					Pixel pix_val = img[p];
					current_level = queue->top_alpha();
					iNode = NewAlphaNode();
					node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
					node[iNode].parentidx = stack_top;
					node[iNode].rootidx = ROOTIDX;
					prev_top = stack_top;
					stack_top = iNode;
					if (current_level > 0)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						node[iNode].rootidx = ROOTIDX;
						parentAry[p] = iNode;
						//prev_top = iNode;
					}
					else
					parentAry[p] = stack_top;
				}
				else
				{
					if (current_level > 0)
						{
							Pixel pix_val = img[p];
							iNode = NewAlphaNode();
							node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
							node[stack_top].add(node + iNode, pix_type);
							node[iNode].parentidx = stack_top;
							node[iNode].rootidx = ROOTIDX;
							parentAry[p] = iNode;
						}
					else
					connectPix2Node(p, img[p], stack_top);
					if (node[stack_top].area == imgsize)
						goto FLOOD_END;
				}
	    	}

			if (node[prev_top].parentidx == stack_top && node[prev_top].area == node[stack_top].area)
			{
				//queue->redundant(node[stack_top].alpha);
				node[prev_top].parentidx = node[stack_top].parentidx;
				stack_top = prev_top;
				curSize--;
			}

			if (node[stack_top].area == imgsize)	// root node found...done
				break;

			//go to higher level
			iNode = node[stack_top].parentidx;
			if ((double)queue->top_alpha() < (double)node[iNode].alpha)
			{
				iNode = NewAlphaNode1(queue->top_alpha(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[iNode].rootidx = ROOTIDX;
				node[stack_top].parentidx = iNode;
			}
			else //go to existing node
			{
			node[iNode].add(node + stack_top, pix_type);
			}
			prev_top = stack_top;
			stack_top = iNode;
			current_level = node[stack_top].alpha;
			if (node[stack_top].area == imgsize)	// root node found...done
				break;
	  }
FLOOD_END:
	  rootidx = (node[stack_top].area == imgsize) ? stack_top : iNode; //remove redundant root
	  node[rootidx].parentidx = ROOTIDX;

	  delete queue;
	  Free(dimg);
	  Free(isVisited);
	  Free(isAvailable);
	}

	void Flood_HierarHeapQueue_Cache_histeq(Pixel* img, int listsize = 12, int a = 0)
	{
		HierarHeapQueue_HEQ<Imgidx, Pixel>* queue;

		Imgidx imgsize, dimgsize, nredges, x0;
		_uint64 numlevels, max_level, current_level;
		Imgidx *dhist;
		Imgidx stack_top, iNode;
		_uint8 *isVisited, *isAvailable, isAv;
		Imgidx p, q;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;


		//double a = 4.0;
		max_level = (sizeof(Pixel)==8) ? 0xffffffffffffffff : (_int64)((Pixel)(-1));
		int bitdepth = get_bitdepth(max_level);

		if (a == 0) a = 1024;
		int eqhistsize = a;
		double coeff = (double)eqhistsize * 10 / (double)bitdepth;

		numlevels = (_uint64)(log2(1 + (double)max_level) * coeff) + 1;

		Pixel *dimg;
		dhist = (Imgidx*)Calloc((size_t)numlevels * sizeof(Imgidx));
		dimg = (Pixel*)Malloc((size_t)dimgsize * sizeof(Pixel));
		compute_dimg(dimg, dhist, img, coeff);//calculate pixel differences and make histogram
		Imgidx *eqhist = (Imgidx*)Calloc(eqhistsize * sizeof(Imgidx));

		_uint32* histeqmap = (_uint32*)Malloc(numlevels * sizeof(_uint32));

		cumsum(dhist, numlevels, histeqmap, eqhistsize);

		for (int ii = 0;ii < (int)numlevels;ii++)
		{
			int jj = histeqmap[ii];
			eqhist[jj] += dhist[ii];
		}

		queue = new HierarHeapQueue_HEQ<Imgidx, Pixel>(eqhist, histeqmap, eqhistsize, nredges, coeff, listsize); // +1 for the dummy node
		curSize = 0;

		maxSize = 1 + imgsize + dimgsize; //Do not use TSE here, becasue dhist is a logged histogram (also this algorithm is for hdr)

		_SET_COMMON_MEMORY
		stack_top = NewAlphaNode();/*dummy root*/
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0, pix_type);
		pNode->parentidx = stack_top;
		current_level = max_level;
		x0 = 0; /*arbitrary starting point*/

		queue->push_1stitem(x0, current_level);
		while (1) //flooding
		{
			while ((_uint64)queue->top_alpha() <= (_uint64)current_level) //flood all levels below current_level
			{
				p = queue->top();

				if (isVisited[p])
				{
					queue->pop(isVisited);
					continue;
				}

				queue->start_pushes();
				isVisited[p] = 1;
				isAv = isAvailable[p];
				if (connectivity == 4)
				{
					q = p << 1;
					if (is_available(isAv, 0) && !isVisited[p + width])	queue->push(p + width, dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1])		queue->push(p + 1, dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1])		queue->push(p - 1, dimg[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])	queue->push(p - width, dimg[q - (width << 1)]);
				}
				else if (connectivity == 8)
				{
					Imgidx width4 = width << 2;
					q = p << 2;
					if (is_available(isAv, 0) && !isVisited[p + width])		queue->push(p + width,		dimg[q]);
					if (is_available(isAv, 1) && !isVisited[p + width + 1])	queue->push(p + width + 1,	dimg[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p + 1])			queue->push(p + 1,			dimg[q + 2]);
					if (is_available(isAv, 3) && !isVisited[p - width + 1])	queue->push(p - width + 1,	dimg[q + 3]);
					if (is_available(isAv, 4) && !isVisited[p - width])		queue->push(p - width,		dimg[q - width4]);
					if (is_available(isAv, 5) && !isVisited[p - width - 1])	queue->push(p - width - 1,	dimg[q - width4 - 3]);
					if (is_available(isAv, 6) && !isVisited[p - 1])			queue->push(p - 1,			dimg[q - 2]);
					if (is_available(isAv, 7) && !isVisited[p + width - 1])	queue->push(p + width - 1,	dimg[q + width4 - 1]);
				}
				else
				{
					//?
				}

				queue->end_pushes(isVisited);
				if ((_uint64)current_level > (_uint64)queue->top_alpha()) //go to lower level
				{
					_CREATE_NEW_NODE(queue->top_alpha())
					if (current_level)
					{
						iNode = NewAlphaNode(0, node + stack_top);
						node[iNode].parentidx = stack_top;
						node[iNode].rootidx = ROOTIDX;
						parentAry[p] = iNode;
						//prev_top = iNode;
					}
					else
						parentAry[p] = stack_top;
				}
				else
				{
					if (current_level)
						_CREATE_SINGLETON_NODE_1
					else
						connectPix2Node(p, img[p], stack_top);
					if (node[stack_top].area == imgsize)
						goto FLOOD_END;
				}
			}
			if (node[stack_top].area == imgsize)	// root node found...done
				break;

			//go to higher level
			iNode = node[stack_top].parentidx;
			if ((double)queue->top_alpha() < (double)node[iNode].alpha)
			{
				iNode = NewAlphaNode1(queue->top_alpha(), node + stack_top);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[iNode].rootidx = ROOTIDX;
				node[stack_top].parentidx = iNode;
			}
			else //go to existing node
			{
				node[iNode].add(node + stack_top, pix_type);
			}

			stack_top = iNode;
			current_level = (_uint64)node[stack_top].alpha;
			if (node[stack_top].area == imgsize)	// root node found...done
				break;
		}
	FLOOD_END:
		rootidx = (node[stack_top].area == imgsize) ? stack_top : iNode; //remove redundant root
		node[rootidx].parentidx = ROOTIDX;

		delete queue;
		Free(dimg);
		Free(isVisited);
		Free(isAvailable);
	}

	Imgidx initialize_node(Pixel *img, Pixel *dimg, Pixel maxpixval)
	{
		Imgidx p, imgsize = width * height;
		Imgidx maxdiffidx = 0;
		Pixel maxdiffval = 0;

		for (p = 0; p < maxSize;p++)
		{
			if (p < imgsize)
				node[p].set(1, 0, (double)img[p], img[p], img[p]);
			else
			{
				Imgidx q = p - imgsize;
				if (maxdiffval < dimg[q])
				{
					maxdiffval = dimg[q];
					maxdiffidx = q;
				}
				node[p].set(0, dimg[q], 0.0, maxpixval, 0, pix_type);
			}
		}

		return maxdiffidx;
	}

	void initialize_node1(Pixel *img, RankItem<Imgidx, double> *rankitem, Pixel maxpixval)
	{
		Imgidx p, imgsize = width * height;

		for (p = 0; p < maxSize;p++)
		{
			if (p < imgsize)
				node[p].set(1, 0, (double)img[p], img[p], img[p]);
			else
				node[p].set(0, rankitem[p - imgsize].alpha, 0.0, maxpixval, 0, pix_type);
			node[p].rootidx = node[p].parentidx = ROOTIDX;
		}
	}

	void initialize_node1(Pixel *img, RankItem<Imgidx, double> *rankitem, Pixel maxpixval, Index* rank2rankitem)
	{
		Imgidx p, imgsize = width * height;

		for (p = 0; p < maxSize;p++)
		{
			Imgidx q = p;
			if (p < imgsize)
				node[p].set(1, 0, (double)img[p], img[p], img[p]);
			else
			{
				q = rank2rankitem[p - imgsize];
				node[p].set(0, rankitem[q].alpha, 0.0, maxpixval, 0, pix_type);
				q = q + imgsize;
			}
			node[q].rootidx = node[q].parentidx = ROOTIDX;
		}
	}

	void initialize_node(Pixel *img, RankItem<Imgidx, Pixel> *rankitem, Pixel maxpixval)
	{
		Imgidx p, imgsize = width * height;

		for (p = 0; p < maxSize;p++)
		{
			if (p < imgsize)
				node[p].set(1, 0, (double)img[p], img[p], img[p]);
			else
				node[p].set(0, rankitem[p - imgsize].alpha, 0.0, maxpixval, 0, pix_type);
			node[p].rootidx = node[p].parentidx = ROOTIDX;
		}

	}

	void initialize_node_par(Pixel *img, RankItem<Imgidx, Pixel> *rankitem, Pixel maxpixval)
	{
		Imgidx p, imgsize = width * height;

		#pragma omp parallel for schedule(guided, 1)
		for (p = 0; p < maxSize;p++)
		{
			if (p < imgsize)
			{
				node[p].set(1, 0, (double)img[p], img[p], img[p]);
				node[p].parentidx = node[p].rootidx = ROOTIDX;
			}
			else
			{
				node[p].set(0, rankitem[p - imgsize].alpha, 0.0, maxpixval, 0, pix_type);
				node[p].parentidx = node[p].rootidx = ROOTIDX;
			}
			//node[p].print(node);

			//node[p].thread = -1;
		}
	}

	void initialize_node_par1(Pixel *img, RankItem<Imgidx, double> *rankitem, Pixel maxpixval, Index* rank2rankitem)
	{
		Imgidx p, imgsize = width * height;

		#pragma omp parallel for schedule(guided, 1)
		for (p = 0; p < maxSize;p++)
		{
			if (p < imgsize)
			{
				node[p].set(1, 0, (double)img[p], img[p], img[p]);
				node[p].parentidx = node[p].rootidx = ROOTIDX;
			}
			else
			{
				if (rank2rankitem)
					node[p].set(0, rankitem[rank2rankitem[p - imgsize]].alpha, 0.0, maxpixval, 0, pix_type);
				else
					node[p].set(0, rankitem[p - imgsize].alpha, 0.0, maxpixval, 0, pix_type);
				node[p].parentidx = node[p].rootidx = ROOTIDX;
			}
		}
	}

	void init_hypergraph_nodes(Pixel *dimg)
	{
		if (connectivity == 4)
		{
			Imgidx imgsize = height * width;
			Imgidx p = 0;
			Imgidx wstride = width << 1;
			for (Imgidx y = 0;y < height;y++)
			{
				for (Imgidx x = 0;x < width;x++)
				{
					Imgidx q = p << 1;
					Pixel minalpha = (Pixel)(-1), minneighidx=0;

					if ((y < height - 1) && (dimg[q] < minalpha))		{minalpha = dimg[q];			minneighidx = q;}
					if ((x < width - 1) && (dimg[q + 1] < minalpha))	{minalpha = dimg[q + 1]; 		minneighidx = q + 1;}
					if ((x > 0) && (dimg[q - 1] < minalpha))			{minalpha = dimg[q - 1]; 		minneighidx = q - 1;}
					if ((y > 0) && (dimg[q - wstride] < minalpha))		{minalpha = dimg[q - wstride];	minneighidx = q - wstride;}

					node[p++].connect_to_parent(&node_in[minneighidx], minneighidx + imgsize, pix_type);
				}
			}
		}
		else
		{

		}
	}

	//connect pixels to one of incident edges with the minimum edge weight, so that
	//pixels are no longer needed to be inspected to make the min-tree implementation
	//easier
	void init_hypergraph_nodes(Imgidx *rank)
	{
		if (connectivity == 4)
		{
			Imgidx imgsize = height * width;
			for (Imgidx p = 0; p < imgsize;p++)
			{
				Imgidx q = p << 1;
				Imgidx y = p / width;
				Imgidx x = p % width;
				Imgidx minRank = 2 * imgsize;

				((y < height - 1) && (rank[q] < minRank)) 		? (minRank = rank[q])					: (Imgidx)0;
				((x < width - 1) && (rank[q + 1] < minRank)) 	? (minRank = rank[q + 1])				: (Imgidx)0;
				((x > 0) && (rank[q - 1] < minRank)) 			? (minRank = rank[q - 1])				: (Imgidx)0;
				((y > 0) && (rank[q - (width << 1)] < minRank)) ? (minRank = rank[q - (width << 1)])	: (Imgidx)0;

				node[p].connect_to_parent(&node_in[minRank], minRank + imgsize, pix_type);
			}
		}
		else
		{

		}
	}

	void set_isAvailable_hypergraph(_uint8 *isAvailable)
	{
		if (connectivity == 4)
		{
			Imgidx dimgidx;
			Imgidx width2 = 2*width;
			Imgidx dimgsize = width * height * 2;

			//first row
			for (dimgidx = 0;dimgidx < width2;)
			{
				isAvailable[dimgidx++] |= 0x20;//even neighbor 5
				isAvailable[dimgidx++] |= 0x30;//odd neighbor 4, 5
			}

			for (dimgidx-=3;dimgidx < dimgsize;dimgidx += width2)
			{
				isAvailable[dimgidx] |= 0x02;//odd neighbor 1
				isAvailable[dimgidx + 1] |= 0x09;//even neighbor 0, 3
			}

			for (dimgidx = 0;dimgidx < dimgsize;dimgidx += width2)
			{
				isAvailable[dimgidx] |= 0x12; //even neighbor 1, 4
				isAvailable[dimgidx + 1] |= 0x08; //odd neighbor 3
			}

			//last 2 rows
			for (dimgidx = width * (height - 2) * 2;dimgidx < dimgsize - width2;dimgidx += 2)
				isAvailable[dimgidx] |= 0x04;//even neighbor 2
			for (dimgidx+=1;dimgidx < dimgsize;dimgidx += 2)
				isAvailable[dimgidx] |= 0x05;//odd neighbor 0, 2

		}
		else//later
		{

		}
	}

	inline _uint8 push_neighbor(Trie<Imgidx, trieidx> *queue, _uint8 *isVisited, Imgidx *rank, Imgidx p)
	{
		isVisited[p] = 1;
		return queue->push(rank[p]);
	}

	#define _RANK_PUSH_NEIGHBORS_4N_HYPERGRAPH_EVEN	\
	((is_available(isAv,0) && !isVisited[p + 1]) ? push_neighbor(queue, isVisited, rank, p + 1) : 0) || \
	((is_available(isAv,1) && !isVisited[p - 1]) ? push_neighbor(queue, isVisited, rank, p - 1) : 0) || \
	((is_available(isAv,2) && !isVisited[p + wstride_d]) ? push_neighbor(queue, isVisited, rank, p + wstride_d) : 0) || \
	((is_available(isAv,3) && !isVisited[p + wstride_d + 1]) ? push_neighbor(queue, isVisited, rank, p + wstride_d + 1) : 0) || \
	((is_available(isAv,4) && !isVisited[p + wstride_d - 1]) ? push_neighbor(queue, isVisited, rank, p + wstride_d - 1) : 0) || \
	((is_available(isAv,5) && !isVisited[p - wstride_d]) ? push_neighbor(queue, isVisited, rank, p - wstride_d) : 0)

	#define _RANK_PUSH_NEIGHBORS_4N_HYPERGRAPH_ODD	\
	((is_available(isAv,0) && !isVisited[p + 1]) ? push_neighbor(queue, isVisited, rank, p + 1) : 0) || \
	((is_available(isAv,1) && !isVisited[p + 2]) ? push_neighbor(queue, isVisited, rank, p + 2) : 0) || \
	((is_available(isAv,2) && !isVisited[p - 1]) ? push_neighbor(queue, isVisited, rank, p - 1) : 0) || \
	((is_available(isAv,3) && !isVisited[p - 2]) ? push_neighbor(queue, isVisited, rank, p - 2) : 0) || \
	((is_available(isAv,4) && !isVisited[p - wstride_d - 1]) ? push_neighbor(queue, isVisited, rank, p - wstride_d - 1) : 0) || \
	((is_available(isAv,5) && !isVisited[p - wstride_d + 1]) ? push_neighbor(queue, isVisited, rank, p - wstride_d + 1) : 0)

	template <class Queue>
	void Flood_Trie_Hypergraph(Pixel *img, Queue *queue)
	{
			Imgidx imgsize, dimgsize, nredges;
			Imgidx current_rank = 0, next_rank = 0;
			RankItem<Imgidx, double> *rankitem, *pRank;
			Pixel maxpixval;
			Imgidx *rank;
			_int8 nbits;
			_uint8 *isVisited, *isAvailable, isAv;
			Imgidx p;
			imgsize = width * height;
			nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
			dimgsize = (connectivity >> 1) * width * height;
			maxSize = imgsize + nredges;
			num_node = maxSize;
			num_node_in = nredges;
			nbits = ((sizeof(Pixel) << 3) - 1);
			maxpixval = ~(1 << nbits);
			rankitem = (RankItem<Imgidx, double>*)Malloc(nredges * sizeof(RankItem<Imgidx, double>));
			parentAry = 0;
			rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
			node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
			node_in = node + imgsize;

			Imgidx wstride_d = width << 1;

			queue = new Queue(nredges);

			omp_set_num_threads(1);
			Index* rank2rankitem = (Index*)Calloc(nredges * sizeof(Index));
			compute_difference_and_sort(rank, rankitem, img, nredges, rank2rankitem);
			initialize_node1(img, rankitem, maxpixval, rank2rankitem);

			init_hypergraph_nodes(rank);

			isVisited = (_uint8*)Calloc(dimgsize * sizeof(_uint8));
			isAvailable = (_uint8*)Calloc(dimgsize * sizeof(_uint8));

			set_isAvailable_hypergraph(isAvailable);

			current_rank =  nredges - 1;
			queue->push(nredges - 1);
			isVisited[rank[nredges]-1] = 1;

			while (1)
			{
				while (1)
				{
					current_rank = queue->top();

					pRank = rankitem + rank2rankitem[current_rank];
					p = pRank->dimgidx;

					_uint8 gotolowerlevel = 0;
					if (connectivity == 4)
					{
						isAv = ~isAvailable[p];

						if (p & 1)	gotolowerlevel = _RANK_PUSH_NEIGHBORS_4N_HYPERGRAPH_ODD;
						else			gotolowerlevel = _RANK_PUSH_NEIGHBORS_4N_HYPERGRAPH_EVEN;
					}
					else
					{
						
					}

					if (!gotolowerlevel)
						break;
				}

				queue->pop();
				next_rank = queue->top();

				node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank + imgsize, pix_type);
				if (node_in[next_rank].area == imgsize)
					break;

				current_rank = next_rank;
			}

			rootidx = (node_in[current_rank].area == imgsize) ? current_rank + imgsize : next_rank + imgsize;
			node[rootidx].parentidx = ROOTIDX;

			delete queue;
			Free(rank2rankitem);
			Free(rank);
			Free(rankitem);
			Free(isVisited);
			Free(isAvailable);
	}

	/*
	6 - neighbor arrangement in isAvailable array
	(order of neighbors is designed to minimize cache miss)
	x - image edge (hypergraph node)
	◯ - neighbouring edge
	+ - image pixel (incident hypergraph edge)

	even-numbered edges
		  5
		  ◯
	1  ◯  +  ◯ 0
		  x
	4  ◯  +  ◯ 3
		  ◯
		  2

	odd-numbered edges
 		 4     5
 		 ◯     ◯
 	3 ◯  +  x  +  ◯ 1
 	     ◯     ◯
 		 2     0

	in isAvailable array:

								MSB           LSB
	isAvailable[i]	 X X # # # # # #
	Neighbor index   - - 5 4 3 2 1 0

	X : don't care
	# = 0 : neighbor not available (side, corner or etc.)
	# = 1 : neighbor available
	*/
	void set_isAvailable_par_hypergraph(_uint8* isAvailable, _int8 npartition_x, _int8 npartition_y)
	{
		Imgidx wstride = width / npartition_x * 2;
		Imgidx wres2 = (width % npartition_x) * 2;
		Imgidx hstride = height / npartition_y;
		Imgidx hres = height % npartition_y;

		if (connectivity == 4)
		{
			Imgidx dimgidx;
			Imgidx width2 = 2*width;
			Imgidx dimgsize = width * height * 2;

			//subimage first rows
			for (int y = 0;y < npartition_y;y++)
			{
				Imgidx nextrowidx = (y * hstride + 1) * width2;
				for (dimgidx = (y * hstride) * width2;dimgidx < nextrowidx;)
				{
					isAvailable[dimgidx++] |= 0x20;//even neighbor 5
					isAvailable[dimgidx++] |= 0x30;//odd neighbor 4, 5
				}
			}

			//subimage right edges
			for (int x = 0;x < npartition_x;x++)
			{
				if (x == 0)
				{
					dimgidx = width2 - 3;
					for (;dimgidx < dimgsize;dimgidx += width2)
					{
						isAvailable[dimgidx] |= 0x02;//odd neighbor 1
						isAvailable[dimgidx + 1] |= 0x09;//even neighbor 0, 3
					}
				}
				else
				{
					//edges on subblock borders belong to subblocks on the left
					dimgidx = width2 - 1 - x * wstride - wres2;
					for (;dimgidx < dimgsize;dimgidx += width2)
						isAvailable[dimgidx] |= 0x23;//odd neighbor 0, 1, 5
				}
			}

			//subimage left edges
			for (int x = 0;x < npartition_x;x++)
			{
				for (dimgidx = x * wstride;dimgidx < dimgsize;dimgidx += width2)
				{
					isAvailable[dimgidx] |= 0x12; //even neighbor 1, 4
					isAvailable[dimgidx + 1] |= 0x08; //odd neighbor 3
				}
			}

			//last 2 rows
			for (int y = 0;y < npartition_y;y++)
			{
				if (y == 0)
				{

					Imgidx subimgend = (height - 1) * width2;
					for (dimgidx = (height - 2) * width2;dimgidx < subimgend;dimgidx += 2)
						isAvailable[dimgidx] |= 0x04;//even neighbor 2
					subimgend = height * width2;
					for (dimgidx = (height - 1) * width2 + 1;dimgidx < subimgend;dimgidx += 2)
						isAvailable[dimgidx] |= 0x05;//odd neighbor 0, 2
				}
				else
				{
					dimgidx = (height - 1 - y * hstride - hres) * width2;
					Imgidx subimgend = dimgidx + width2;
					for (;dimgidx < subimgend;dimgidx += 2)
						isAvailable[dimgidx] |= 0x1c;//even neighbor 2, 3, 4
				}
			}

		}
		else
		{

		}
	}

	void cumsum(Imgidx *hist, Imgidx size, Imgidx &maxidx)
	{
		Imgidx sum = hist[0], hi;
		maxidx = 0;
		for (Imgidx i = 1;i < size;i++)
		{
			if (hist[i])
			{
				maxidx = i;
				sum += hist[i];
			}
			hist[i] = sum;
		}
	}

	void cumsum(Imgidx *hist, Imgidx size, _uint32 *histeqmap, int eqhistsize)
	{
		Imgidx sum = hist[0];
		for (Imgidx i = 1;i < size;i++)
		{
			sum += hist[i];
		}

		double coeff = (double)(eqhistsize - 1) / (double)sum;
		sum = 0;
		for (Imgidx i = 0;i < size;i++)
		{
			sum += hist[i];
			histeqmap[i] = (_uint16)((double)sum * coeff);
		}
	}

	inline _uint8 push_neighbor(HierarQueue<Imgidx> *queue, _uint8 *isVisited, _uint8 *dimg, Imgidx p)
	{
		isVisited[p] = 1;
		return queue->push(p, dimg[p]);
	}

	inline _uint8 push_neighbor(HierarQueue<Imgidx> *queue, _uint8 *isVisited, _uint16 *dimg, Imgidx p)
	{
		isVisited[p] = 1;
		return queue->push(p, dimg[p]);
	}

	inline _uint8 push_neighbor(HierarQueue<Imgidx> *queue, _uint8 *isVisited, _uint32 *dimg, Imgidx p)
	{
		isVisited[p] = 1;
		return queue->push(p, dimg[p]);
	}

	inline _uint8 push_neighbor(HierarQueue<Imgidx> *queue, _uint8 *isVisited, _uint64 *dimg, Imgidx p)
	{
		isVisited[p] = 1;
		return queue->push(p, dimg[p]);
	}

	#define _RANK_PUSH_NEIGHBORS_4N_HYPERGRAPH_PILOT_EVEN	\
	((is_available(isAv,0) && !isVisited[p + 1]) ? push_neighbor(queue, isVisited, dimg, p + 1) : 0) || \
	((is_available(isAv,1) && !isVisited[p - 1]) ? push_neighbor(queue, isVisited, dimg, p - 1) : 0) || \
	((is_available(isAv,2) && !isVisited[p + wstride_d]) ? push_neighbor(queue, isVisited, dimg, p + wstride_d) : 0) || \
	((is_available(isAv,3) && !isVisited[p + wstride_d + 1]) ? push_neighbor(queue, isVisited, dimg, p + wstride_d + 1) : 0) || \
	((is_available(isAv,4) && !isVisited[p + wstride_d - 1]) ? push_neighbor(queue, isVisited, dimg, p + wstride_d - 1) : 0) || \
	((is_available(isAv,5) && !isVisited[p - wstride_d]) ? push_neighbor(queue, isVisited, dimg, p - wstride_d) : 0)

	#define _RANK_PUSH_NEIGHBORS_4N_HYPERGRAPH_PILOT_ODD	\
	((is_available(isAv,0) && !isVisited[p + 1]) ? push_neighbor(queue, isVisited, dimg, p + 1) : 0) || \
	((is_available(isAv,1) && !isVisited[p + 2]) ? push_neighbor(queue, isVisited, dimg, p + 2) : 0) || \
	((is_available(isAv,2) && !isVisited[p - 1]) ? push_neighbor(queue, isVisited, dimg, p - 1) : 0) || \
	((is_available(isAv,3) && !isVisited[p - 2]) ? push_neighbor(queue, isVisited, dimg, p - 2) : 0) || \
	((is_available(isAv,4) && !isVisited[p - wstride_d - 1]) ? push_neighbor(queue, isVisited, dimg, p - wstride_d - 1) : 0) || \
	((is_available(isAv,5) && !isVisited[p - wstride_d + 1]) ? push_neighbor(queue, isVisited, dimg, p - wstride_d + 1) : 0)

	void Flood_HierarQueue_Hypergraph(Pixel* img)
	{
		if (sizeof(Pixel) > 2 || channel > 1)
		{
			printf("Error: Hierarchical queues do not work on >16 bits images or multispectral images\n");
			printf("Try Unionfind (algorithm code %d), flooding using Heapqueue (%d), trie queue (%d) or cached trie queue (%d) \n", UNIONFIND, FLOOD_HEAPQUEUE_CACHE, FLOOD_TRIE, FLOOD_TRIE_CACHE);
			return;
		}

		HierarQueue<Imgidx>* queue;

		Imgidx imgsize, dimgsize, nredges;
		_int64 numlevels, max_level, current_level;
		Imgidx *dhist;
		Pixel *dimg;
		Imgidx stack_top, iNode;
		_uint8 *isVisited, *isAvailable;
		Imgidx p, wstride_d = width << 1;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		max_level = (_int64)((Pixel)(-1));
		numlevels = max_level +1;


		dhist = (Imgidx*)Malloc((size_t)numlevels * sizeof(Imgidx));
		memset(dhist, 0, (size_t)numlevels * sizeof(_int32));
		dimg = (Pixel*)Calloc((size_t)dimgsize * sizeof(Pixel));

		compute_dimg(dimg, dhist, img);//calculate pixel differences and make histogram

		maxSize = imgsize + dimgsize;
		node = (AlphaNode<Imgidx, Pixel>*)Calloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		node_in = node + imgsize;

		Imgidx maxdiffidx = initialize_node(img, dimg, max_level);
		max_level = dimg[maxdiffidx];

		init_hypergraph_nodes(dimg);

		//create hierarchical queue from dhist
		queue = new HierarQueue<Imgidx>(nredges, dhist, numlevels);
		curSize = 0;

		Free(dhist);
		isVisited = (_uint8*)Calloc(dimgsize * sizeof(_uint8));
		isAvailable = (_uint8*)Calloc(dimgsize * sizeof(_uint8));

		omp_set_num_threads(1);
		_int8 npartition_x = 1, npartition_y = 1;
		set_isAvailable_par_hypergraph(isAvailable, npartition_x, npartition_y);

		stack_top = maxdiffidx + imgsize;//root
		AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
		pNode->parentidx = stack_top;
		current_level = max_level;
		queue->push(maxdiffidx, current_level);
		isVisited[maxdiffidx] = 1;
		while (1) //flooding
		{
			while (1) //flood all levels below current_level
			{
				p = queue->top();

				_uint8 gotolowerlevel = 0;
				if (connectivity == 4)
				{
					_uint8 isAv = ~isAvailable[p];

					if (p & 1)	gotolowerlevel = _RANK_PUSH_NEIGHBORS_4N_HYPERGRAPH_PILOT_ODD;
					else      gotolowerlevel = _RANK_PUSH_NEIGHBORS_4N_HYPERGRAPH_PILOT_EVEN;
				}
				else
				{
					//LATER
				}

				if (gotolowerlevel)
				{
					current_level = queue->min_level;
					iNode = queue->top() + imgsize;
					node[iNode].parentidx = stack_top;
					stack_top = iNode;
				}
				else
				{
					queue->pop();
					queue->find_minlev();

					if (queue->min_level == current_level)
					{
						iNode = queue->top() + imgsize;
						node[stack_top].add(node + iNode, pix_type);
						node[iNode].parentidx = stack_top;
					}
					else
						break;
				}
			}
			//go to higher level
			iNode = node[stack_top].parentidx;
			if ((_int64)queue->min_level < (_int64)node[iNode].alpha) //new level from queue
			{
				iNode = queue->top() + imgsize;
				node[iNode].add(node + stack_top,pix_type);
				node[iNode].parentidx = node[stack_top].parentidx;
				node[stack_top].parentidx = iNode;
			}
			else //go to existing node
			{
				if (node[iNode].area == imgsize)	// root node found...done
						break;
				node[iNode].add(node + stack_top, pix_type);
			}

			if (node[iNode].area == imgsize)	// root node found...done
					break;

			stack_top = iNode;
			current_level = node[stack_top].alpha;
		}

		rootidx = node[iNode].area == imgsize ? iNode : stack_top;
		node[rootidx].parentidx = ROOTIDX;

		delete queue;
		Free(dimg);
		Free(isVisited);
		Free(isAvailable);
	}

	void canonicalize(Imgidx nidx)
	{
		Imgidx p, q;

		p = get_level_root(nidx);//for 0-ccs
		if (p != nidx)
		 	node[nidx].parentidx = p;

		while(1)
		{
			q = node[p].parentidx;
			if (q == ROOTIDX)
				break;
			q = get_level_root(q);
			node[p].parentidx = q;
			p = q;
		}
	}

	inline Imgidx merge_subtrees(Pixel *dimg, _int64 blksz_x, _int64 blksz_y, Imgidx npartition_x, Imgidx npartition_y, Imgidx* subtree_cur, Imgidx* subtree_start = NULL, Imgidx *blkhs = NULL, Imgidx *blkws = NULL, Imgidx *nrbnode = NULL, Imgidx *subtree_nborderedges = NULL)
	{
		return merge_subtrees(dimg, blksz_x, blksz_y, npartition_x, npartition_y, subtree_cur, 0, subtree_start, blkhs, blkws);
	}

	//returns root node index
	Imgidx merge_subtrees(Pixel *dimg, _int64 blksz_x, _int64 blksz_y, Imgidx npartition_x, Imgidx npartition_y, Imgidx* subtree_cur, int tse, Imgidx* subtree_start = NULL, Imgidx *blkhs = NULL, Imgidx *blkws = NULL, Imgidx *nrbnode = NULL, Imgidx *subtree_nborderedges = NULL)
	{
		Imgidx numblk;
		_int64 blksz_x0 = blksz_x;
		_int64 blksz_y0 = blksz_y;


		Imgidx npartition_x0 = npartition_x;
		Imgidx npartition_y0 = npartition_y;
		Imgidx blkrow = width * blksz_y0;
		while(npartition_x > 1 || npartition_y > 1)
		{
				//merge horizontal borders
			if ((npartition_x == 1 || blksz_x >= blksz_y) && npartition_y > 1)
			{
				numblk = npartition_x * (npartition_y / 2);

				#pragma omp parallel for schedule(dynamic,1)
				for (int blk = 0;blk < numblk;blk++)
				{
					Imgidx x, y, r, p, p0, pn, dimgidx, bx, by, bidx;
					//Pixel qminlev;
					y = (1 + 2 * (blk / (int)npartition_x)) * blksz_y;
					x = (blk % (int)npartition_x) * blksz_x;

					p0 = (y - 1) * width + x;
					pn = (((blk % (int)npartition_x) == npartition_x - 1) ? y * width : p0 + blksz_x);

					by = _min((p0 / blkrow), npartition_y0 - 1);
						for (p = p0;p < pn;p++)
					{
						bx = _min(((p % width) / blksz_x0), npartition_x0 - 1);
						bidx = by * npartition_x0 + bx;

						dimgidx = p << 1;
						r = dimg[dimgidx];
						nrbnode[bidx]++;

						if (tse) connect(parentAry[p], parentAry[p + width], (Imgidx)subtree_cur[bidx]++, (Pixel)r);
						else    connect(p, p + width, (Imgidx)subtree_cur[bidx]++, (Pixel)r);
					}
				}
				npartition_y = (npartition_y + 1) / 2;
				blksz_y <<= 1;
				if (npartition_y == 1)
					blksz_y = height;
				else
				blksz_y = _min(blksz_y, height);

			}

			if ((npartition_y == 1 || blksz_x <= blksz_y) && npartition_x > 1)
			{
				numblk = npartition_y * (npartition_x / 2);

				#pragma omp parallel for schedule(dynamic,1)
				for (int blk = 0;blk < numblk;blk++)
				{
					Imgidx x, y, r, p, p0, pn, dimgidx, bx, by, bidx;;
					x = (1 + 2 * (blk / npartition_y)) * blksz_x;
					y = (blk % (int)npartition_y) * blksz_y;

					p0 = y * width + x - 1;
					pn = ((blk % (int)npartition_y) == npartition_y - 1) ? height * width : p0 + width * blksz_y;

					bx = _min(((p0 % width) / blksz_x0), npartition_x0 - 1);
					for (p = p0;p < pn;p += width)
					{
						by = _min((p / blkrow),npartition_y0 - 1);
						bidx = by * npartition_x0 + bx;

						dimgidx = (p << 1) + 1;
						r = dimg[dimgidx];
						nrbnode[bidx]++;

						if (tse) connect(parentAry[p], parentAry[p + 1], (Imgidx)subtree_cur[bidx]++, (Pixel)r);
						else    connect(p, p + 1, (Imgidx)subtree_cur[bidx]++, (Pixel)r);
					}
				}
				npartition_x = (npartition_x + 1) / 2;
				blksz_x <<= 1;
				if (npartition_x == 1)
					blksz_x = width;
				else
					blksz_x = _min(blksz_x, width);
			}
		}

		Imgidx p;
		if (tse)
		 	p = parentAry[0];
		else
			p = 0;
		while(node[p].parentidx != ROOTIDX)
			p = node[p].parentidx;

		return p;
	}

	Imgidx merge_subtrees(_uint8 *dimg, _int64 blksz_x, _int64 blksz_y, _int16 npartition_x, _int16 npartition_y, Imgidx* subtree_cur, int tse)
	{
		Imgidx numblk;
		_int64 blksz_x0 = blksz_x;
		_int64 blksz_y0 = blksz_y;

		Imgidx npartition_x0 = npartition_x;
		Imgidx npartition_y0 = npartition_y;
		Imgidx blkrow = width * blksz_y0;
		while(npartition_x > 1 || npartition_y > 1)
		{
			if ((npartition_x == 1 || blksz_x >= blksz_y) && npartition_y > 1)
			{
				numblk = npartition_x * (npartition_y / 2);


				#pragma omp parallel for
				for (int blk = 0;blk < numblk;blk++)
				{
					Imgidx x, y, r, p, p0, pn, dimgidx, bx, by, bidx;
					y = (1 + 2 * (blk / (int)npartition_x)) * blksz_y;
					x = (blk % (int)npartition_x) * blksz_x;

					p0 = (y - 1) * width + x;
					pn = (((blk % (int)npartition_x) == npartition_x - 1) ? y * width : p0 + blksz_x);

					by = _min((p0 / blkrow), npartition_y0 - 1);
				for (p = p0;p < pn;p++)
					{
						bx = _min(((p % width) / blksz_x0), npartition_x0 - 1);
						bidx = by * npartition_x0 + bx;

						dimgidx = p << 1;
						r = dimg[dimgidx];

						if (tse) connect(parentAry[p], parentAry[p + width], (Imgidx)subtree_cur[bidx]++, (Pixel)r);
						else    connect(p, p + width, (Imgidx)subtree_cur[bidx]++, (Pixel)r);
					}
				}
				npartition_y = (npartition_y + 1) / 2;
				blksz_y <<= 1;
				if (npartition_y == 1)
					blksz_y = height;
				else
				blksz_y = _min(blksz_y, height);
			}

			//merge vertical borders
			if ((npartition_y == 1 || blksz_x <= blksz_y) && npartition_x > 1)
			{
				numblk = npartition_y * (npartition_x / 2);

				#pragma omp parallel for
				for (int blk = 0;blk < numblk;blk++)
				{
					Imgidx x, y, r, p, p0, pn, dimgidx, bx, by, bidx;;
					x = (1 + 2 * (blk / npartition_y)) * blksz_x;
					y = (blk % (int)npartition_y) * blksz_y;

					p0 = y * width + x - 1;
					pn = ((blk % (int)npartition_y) == npartition_y - 1) ? height * width : p0 + width * blksz_y;

					bx = _min(((p0 % width) / blksz_x0), npartition_x0 - 1);
					for (p = p0;p < pn;p += width)
					{
						by = _min((p / blkrow),npartition_y0 - 1);
						bidx = by * npartition_x0 + bx;

						dimgidx = (p << 1) + 1;
						r = dimg[dimgidx];

						if (tse)	connect(parentAry[p], parentAry[p + 1], (Imgidx)subtree_cur[bidx]++, (Pixel)r);
						else		connect(p, p + 1, (Imgidx)subtree_cur[bidx]++, (Pixel)r);

					}
				}
				npartition_x = (npartition_x + 1) / 2;
				blksz_x <<= 1;
				if (npartition_x == 1)
					blksz_x = width;
				else
					blksz_x = _min(blksz_x, width);
			}
		}

		Imgidx p;
		if (tse)
		 	p = parentAry[0];
		else
			p = 0;
		while(node[p].parentidx != ROOTIDX)
			p = node[p].parentidx;

		return p;
	}

	Imgidx merge_subtrees1(_uint8 *dimg, _int64 blksz_x, _int64 blksz_y, _int16 npartition_x, _int16 npartition_y, Imgidx* subtree_cur, int tse, Imgidx* hypernode_level)
	{
		Imgidx numblk;
		_int64 blksz_x0 = blksz_x;
		_int64 blksz_y0 = blksz_y;

		Imgidx npartition_x0 = npartition_x;
		Imgidx npartition_y0 = npartition_y;
		Imgidx blkrow = width * blksz_y0;
		while(npartition_x > 1 || npartition_y > 1)
		{
			if ((npartition_x == 1 || blksz_x >= blksz_y) && npartition_y > 1)
			{
				numblk = npartition_x * (npartition_y / 2);

				#pragma omp parallel for schedule(dynamic,1)
				for (int blk = 0;blk < numblk;blk++)
				{
					Imgidx x, y, r, p, p0, pn, dimgidx, bx, by, bidx;
					y = (1 + 2 * (blk / (int)npartition_x)) * blksz_y;
					x = (blk % (int)npartition_x) * blksz_x;

					p0 = (y - 1) * width + x;
					pn = (((blk % (int)npartition_x) == npartition_x - 1) ? y * width : p0 + blksz_x);

					by = _min((p0 / blkrow), npartition_y0 - 1);
					for (p = p0;p < pn;p++)
					{
						bx = _min(((p % width) / blksz_x0), npartition_x0 - 1);
						bidx = by * npartition_x0 + bx;

						dimgidx = p << 1;
						r = dimg[dimgidx];

						if (tse) hypernode_level[dimgidx] = connect(parentAry[p], parentAry[p + width], (Imgidx)subtree_cur[bidx]++, (Pixel)r);
						else    hypernode_level[dimgidx] = connect(p, p + width, (Imgidx)subtree_cur[bidx]++, (Pixel)r);
					}
				}
				npartition_y = (npartition_y + 1) / 2;
				blksz_y <<= 1;
				if (npartition_y == 1)
					blksz_y = height;
				else
				blksz_y = _min(blksz_y, height);

			}

			if ((npartition_y == 1 || blksz_x <= blksz_y) && npartition_x > 1)
			{
				numblk = npartition_y * (npartition_x / 2);

				#pragma omp parallel for schedule(dynamic,1)
				for (int blk = 0;blk < numblk;blk++)
				{
					Imgidx x, y, r, p, p0, pn, dimgidx, bx, by, bidx;;
					x = (1 + 2 * (blk / npartition_y)) * blksz_x;
					y = (blk % (int)npartition_y) * blksz_y;

					p0 = y * width + x - 1;
					pn = ((blk % (int)npartition_y) == npartition_y - 1) ? height * width : p0 + width * blksz_y;

					bx = _min(((p0 % width) / blksz_x0), npartition_x0 - 1);
					for (p = p0;p < pn;p += width)
					{
						by = _min((p / blkrow),npartition_y0 - 1);
						bidx = by * npartition_x0 + bx;

						dimgidx = (p << 1) + 1;
						r = dimg[dimgidx];

						if (tse) hypernode_level[dimgidx] =  connect(parentAry[p], parentAry[p + 1], (Imgidx)subtree_cur[bidx]++, (Pixel)r);
						else    hypernode_level[dimgidx] =  connect(p, p + 1, (Imgidx)subtree_cur[bidx]++, (Pixel)r);

					}
				}
				npartition_x = (npartition_x + 1) / 2;
				blksz_x <<= 1;
				if (npartition_x == 1)
					blksz_x = width;
				else
					blksz_x = _min(blksz_x, width);
			}
		}

		Imgidx p;
		if (tse)
		 	p = parentAry[0];
		else
			p = 0;
		while(node[p].parentidx != ROOTIDX)
			p = node[p].parentidx;

		return p;
	}

	int migrate_subtree(int blk, int numpartitions, Imgidx & nidx, Imgidx & nidx_lim, int & nidxblk,
		Imgidx & blkts, char *blkflooddone, Imgidx *subtree_cur, Imgidx *subtree_start, Imgidx *subtree_nborderedges,
		omp_lock_t *locks, int &numbusythr, int &numblkproc, int &outofmemory)
	{
		if (omp_get_num_threads() == 1)
			outofmemory = 1;

		subtree_cur[nidxblk] = nidx;
		blkts += nidx - subtree_start[nidxblk];
		blkflooddone[nidxblk] = 2;

		omp_set_lock(locks + numpartitions);
		numbusythr--;
		omp_unset_lock(locks + numpartitions);
		omp_unset_lock(locks + nidxblk);

		while(1)
		{
			if (outofmemory || (!numbusythr))// && numblkproc >= omp_get_num_threads()))
			{
				outofmemory = 1;
				return 0;
			}

			for (int newblkidx = 0;newblkidx < numpartitions;newblkidx++)
			{
				if (blkflooddone[newblkidx] && (subtree_cur[newblkidx] + subtree_nborderedges[newblkidx] < subtree_start[newblkidx + 1]) && !omp_test_lock(locks + newblkidx))
				{
					omp_set_lock(locks + numpartitions);

					if (blkflooddone[newblkidx])
					{
						blkflooddone[newblkidx] = 0;
						numbusythr++;
						omp_unset_lock(locks + numpartitions);
					}
					else
					{
						omp_unset_lock(locks + numpartitions);
						continue;
					}
					nidxblk = newblkidx;
					nidx = subtree_cur[nidxblk];
					blkts -= subtree_cur[nidxblk] - subtree_start[nidxblk];
					nidx_lim = subtree_start[nidxblk + 1] - subtree_nborderedges[newblkidx];
					return 1;
				}
			}
		}
	}

	Imgidx parflood_node_alloc(Imgidx *subtree_size, Imgidx *subtree_start, Imgidx *blkws, Imgidx *blkhs, int numpartitions, double sizemult)
	{
		subtree_start[0] = 0;
		for (int blk = 0;blk < numpartitions;blk++)
		{
			Imgidx blkmaxsize = 1 + 3 * blkws[blk] * blkhs[blk];
			subtree_start[blk + 1] = _min(blkmaxsize, (Imgidx)((double)subtree_size[blk] * sizemult)) + subtree_start[blk];
		}
		maxSize = subtree_start[numpartitions];
		if (node)	Free(node);
		node = (AlphaNode<Imgidx, Pixel>*)Calloc((size_t)(maxSize) * sizeof(AlphaNode<Imgidx, Pixel>));

		return maxSize;
	}

	void set_isAvailable_par(_uint8* isAvailable, _int16 npartition_x, _int16 npartition_y)
	{
		_int32 i, j, k;
		Imgidx imgsize = width * height;
		Imgidx wstride = width / npartition_x;
		Imgidx hstride = height / npartition_y;

		set_isAvailable(isAvailable);

		if (connectivity == 4)
		{
			//hor partitions
			j = (hstride - 1) * width;
			for (i = 0; i < npartition_y - 1; i++)
			{
				k = j + width;
				for (; j < k; j++)
				{
					isAvailable[j] &= 0xe;
					isAvailable[j + width] &= 0x7;
				}
				j += (hstride - 1) * width;
			}

			//ver partitions
			for (i = 0; i < npartition_x - 1; i++)
			{
				j = (i + 1) * wstride - 1;
				for (; j < imgsize; j += width)
				{
					isAvailable[j] &= 0xd;
					isAvailable[j + 1] &= 0xb;
				}
			}
		}
		else
		{
			//		    Neighbour Index
			// 			6      5      4
			// 			7    pixel    3
			// 			0      1      2
			//
			//			Neighbour indices to bit field
			//			7 6 5 4 3 2 1 0
			//         MSB			 LSB
			//			0: Neighbour pixel not available (corner of Image, or partition in later implementation)
			//			1: available

			//initialize to all available
			for (i = 0; i < imgsize; i++)
				isAvailable[i] = 0xff;

			//four corners
			isAvailable[0] = 0x0e;
			isAvailable[width - 1] = 0x83;
			isAvailable[width*(height - 1)] = 0x38;
			isAvailable[width * height - 1] = 0xe0;

			//top and bottom row
			j = width * (height - 1) + 1;
			for (i = 1; i < width - 1; i++)
			{
				isAvailable[i] = 0x8f;
				isAvailable[j] = 0xf8;
				j++;
			}

			//leftest and rightest column
			j = width;
			k = (width << 1) - 1;
			for (i = 1; i < height - 1; i++)
			{
				isAvailable[j] = 0x3e;
				isAvailable[k] = 0xe3;
				j += width;
				k += width;
			}
		}
	}

	void Flood_Hierarqueue_par(Pixel *img, int numthreads)
	{
		if (sizeof(Pixel) > 2 || channel > 1)
		{
			printf("Error: Hierarchical queues do not work on >16 bits images or multispectral images\n");
			printf("Try Unionfind (algorithm code %d), flooding using Heapqueue (%d), trie queue (%d) or cached trie queue (%d) \n", UNIONFIND, FLOOD_HEAPQUEUE_CACHE, FLOOD_TRIE, FLOOD_TRIE_CACHE);
			return;
		}

		Imgidx imgsize, dimgsize;
		_int64 numlevels;
		Pixel *dimg;
		_uint8 *isVisited, *isAvailable;
		Imgidx p, q;
		imgsize = width * height;
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = (sizeof(Pixel) == 1) ? 256 : 65536;

		dimg = (Pixel*)Calloc((size_t)dimgsize * sizeof(Pixel));

		_int16 npartition_x, npartition_y;
		{
			_int16 optpart = 1;
			double optborderlength = (double)numthreads * (double)imgsize;
			for (int px = 2;px < numthreads; px++)
			{
				if (numthreads % px == 0)
				{
					int py = numthreads / px;

					if (((double)px * (double)height + (double)py * (double)width) < optborderlength)
					{
						optpart = px;
						optborderlength = ((double)px * (double)height + (double)py * (double)width);
					}
				}
			}
			npartition_x = (_int16)optpart;
			npartition_y = (_int16)numthreads / npartition_x;
		}

		isVisited = (_uint8*)Calloc((size_t)((imgsize)));
		isAvailable = (_uint8*)Malloc((size_t)(imgsize));
		set_isAvailable_par(isAvailable, npartition_x, npartition_y);

		_int64 blksz_x = width / npartition_x;
		_int64 blksz_y = height / npartition_y;
		_int64 blksz_xn = blksz_x + (width % npartition_x);
		_int64 blksz_yn = blksz_y + (height % npartition_y);
		_int64 numpartitions = (_int64)npartition_x * (_int64)npartition_y;

		p = q = 0;
		Imgidx *startpidx = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blocksize = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *subtree_start = (Imgidx*)Malloc((numpartitions + 1) * sizeof(Imgidx));
		Imgidx *subtree_cur = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blkws = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blkhs = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *dhist = (Imgidx*)Calloc((size_t)numlevels * (size_t)numpartitions * sizeof(Imgidx));

		for (_int16 y = 0; y < npartition_y; y++)
		{
			q = y * width * (Imgidx)blksz_y;
			bool lastrow = (y == npartition_y - 1);
			Imgidx blkh = lastrow ? blksz_yn : blksz_y;
			for (_int16 x = 0; x < npartition_x; x++)
			{
				startpidx[p] = q + (Imgidx)x * (Imgidx)blksz_x;
				bool lastcol = (x == npartition_x - 1);
				Imgidx blkw = lastcol ? blksz_xn : blksz_x;
				blocksize[p] = blkh * blkw * 2 - (Imgidx)lastrow * blkw - (Imgidx)lastcol * blkh;
				blkws[p] = blkw;
				blkhs[p] = blkh;
				p++;
			}
		}

		subtree_start[0] = startpidx[0] + imgsize;
		for (int blk = 1; blk < numpartitions;blk++)
		{
			subtree_start[blk] = subtree_start[blk - 1] + (blkws[blk - 1] * blkhs[blk - 1] * 2);
		}
		subtree_start[numpartitions] = subtree_start[numpartitions - 1] + (blkws[numpartitions - 1] * blkhs[numpartitions - 1] * 2);


		HierarQueue<Imgidx>** queues;
		queues = (HierarQueue<Imgidx>**)Calloc(numpartitions * sizeof(HierarQueue<Imgidx>*));
		for (int blk = 0;blk < numpartitions;blk++)
			queues[blk] = new HierarQueue<Imgidx>((_uint64)blocksize[blk] + 1);

		//singletons + inners + dummies
		node = (AlphaNode<Imgidx, Pixel>*)Calloc((size_t)(imgsize + dimgsize + numpartitions) * sizeof(AlphaNode<Imgidx, Pixel>));
		node_in = node + imgsize;

		#pragma omp parallel for private(p, q) schedule(dynamic,1)
		for (int blk = 0; blk < numpartitions; blk++)
		{
			Imgidx bwidth = blkws[blk];
			Imgidx bheight = blkhs[blk];
			Imgidx bareasum = bwidth * bheight;
			Imgidx *bhist = dhist + numlevels * blk;
			HierarQueue<Imgidx> *queue = queues[blk];
			Imgidx spidx = startpidx[blk];
			Imgidx nidx = subtree_start[blk];
			Imgidx iNode;
			Pixel maxdiff = 0;


			for (Imgidx i = 0;i < bheight;i++)
			{
				Pixel diff;
				p = spidx + i * width;
				bool notlastrow = i < bheight - 1;
				for (Imgidx j = 0;j < bwidth - 1;j++)
				{
					q = p << 1;
					node[p].set(1, 0, (double)img[p], img[p], img[p]);
					node[p].parentidx = node[p].rootidx = ROOTIDX;

					if (notlastrow)
					{
						diff = abs_diff(img[p], img[p + width]);
						dimg[q] = diff;
						bhist[diff]++;
						maxdiff = _max(maxdiff, diff);
					}
					diff = abs_diff(img[p], img[p + 1]);
					dimg[q + 1] = diff;
					bhist[diff]++;
					maxdiff = _max(maxdiff, diff);
					p++;
				}
				q = p << 1;
				node[p].set(1, 0, (double)img[p], img[p], img[p]);
				node[p].parentidx = node[p].rootidx = ROOTIDX;
				if (notlastrow)
				{
					diff = abs_diff(img[p], img[p + width]);
					dimg[q] = diff;
					bhist[diff]++;
					maxdiff = _max(maxdiff, diff);
				}
			}
			bhist[maxdiff]++;

			queue->set_queue(bhist, maxdiff);

			Imgidx stack_top = imgsize + dimgsize + blk;
			Imgidx prev_top = stack_top;
			AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
			pNode->set(0, maxdiff, (double)0.0, maxdiff, (Pixel)0, pix_type);
			pNode->parentidx = ROOTIDX;
			Pixel current_level = maxdiff;
			queue->push(startpidx[blk], current_level);
			while (1) //flooding
			{
				while ((_int64)queue->min_level <= (_int64)current_level) //flood all levels below current_level
				{
					p = queue->pop();
					if (is_visited(isVisited, p))
					{
						queue->find_minlev();
						continue;
					}

					isVisited[p] = 1;
					_uint8 isAv = isAvailable[p];

					if (connectivity == 4)
					{
						q = p << 1;
						(is_available(isAv, 0) && !isVisited[p + width]) ? (void)queue->push(p + width, dimg[q]) : (void)0;
						(is_available(isAv, 1) && !isVisited[p + 1]) 		 ? (void)queue->push(p + 1, dimg[q + 1]) : (void)0;
						(is_available(isAv, 2) && !isVisited[p - 1]) 		 ? (void)queue->push(p - 1, dimg[q - 1]) : (void)0;
						(is_available(isAv, 3) && !isVisited[p - width]) ? (void)queue->push(p - width, dimg[q - (width << 1)]) : (void)0;
					}
					//else 				  				 _PUSH_NEIGHBORS_8N //later!

					if ((_int64)current_level > (_int64)queue->min_level) //go to lower level
					{
						//_CREATE_NEW_NODE(queue->min_level)
						Pixel pix_val = node[p].minPix;
						current_level = queue->min_level;

						iNode = nidx++;
						node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
						node[iNode].parentidx = stack_top;
						node[iNode].rootidx = ROOTIDX;
						node[p].parentidx = iNode;
						stack_top = iNode;
					}
					else
					{
						queue->find_minlev();
						node[stack_top].add(node + p, pix_type);
						node[p].parentidx = stack_top;
					}
				}

				remove_redundant_node(node, nidx, prev_top, stack_top);

				//go to higher level
				iNode = node[stack_top].parentidx;
				if (iNode == ROOTIDX || (_int64)queue->min_level < (_int64)node[iNode].alpha) //new level from queue
				{
					iNode = nidx++;
					node[iNode].alpha = queue->min_level;
					node[iNode].copy(node + stack_top);
					node[iNode].parentidx = node[stack_top].parentidx;
					node[iNode].rootidx = ROOTIDX;
					node[stack_top].parentidx = iNode;
				}
				else //go to existing node
				{
					if (node[iNode].area == bareasum)
						break;
					node[iNode].add(node + stack_top, pix_type);
				}

				if (node[iNode].area == bareasum)
					break;

				prev_top = stack_top;
				stack_top = iNode;
				current_level = node[stack_top].alpha;
			}
			stack_top = (node[stack_top].area == bareasum) ? stack_top : iNode; //remove redundant root
			node[stack_top].parentidx = ROOTIDX;

			subtree_cur[blk] = nidx;
		}

		merge_subtrees(dimg, blksz_x, blksz_y, npartition_x, npartition_y, subtree_cur);

		Free(isVisited);
		Free(isAvailable);

		for (int blk = 0;blk < numpartitions;blk++)
			delete queues[blk];
		Free(queues);
		Free(dimg);
		Free(startpidx);
		Free(blocksize);
		Free(subtree_start);
		Free(subtree_cur);
		Free(blkws);
		Free(blkhs);
		Free(dhist);
	}

	void Flood_Hierarqueue_par_tse(Pixel *img, int numthreads, int tse)
	{
		if (sizeof(Pixel) > 2 || channel > 1)
		{
			printf("Error: Hierarchical queues do not work on >16 bits images or multispectral images\n");
			printf("Try Unionfind (algorithm code %d), flooding using Heapqueue (%d), trie queue (%d) or cached trie queue (%d) \n", UNIONFIND, FLOOD_HEAPQUEUE_CACHE, FLOOD_TRIE, FLOOD_TRIE_CACHE);
			return;
		}

		Imgidx imgsize, dimgsize;
		_int64 numlevels;
		Pixel *dimg;
		_uint8 *isVisited, *isAvailable;
		Imgidx p, q;
		imgsize = width * height;
//		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = (sizeof(Pixel) == 1) ? 256 : 65536;
		dimg = (Pixel*)Calloc((size_t)dimgsize * sizeof(Pixel));

		_int16 npartition_x, npartition_y;
		{
			_int16 optpart = 1;
			double optborderlength = (double)numthreads * (double)imgsize;
			for (int px = 2;px < numthreads; px++)
			{
				if (numthreads % px == 0)
				{
					int py = numthreads / px;

					if (((double)px * (double)height + (double)py * (double)width) < optborderlength)
					{
						optpart = px;
						optborderlength = ((double)px * (double)height + (double)py * (double)width);
					}
				}
			}
			npartition_x = (_int16)optpart;
			npartition_y = (_int16)numthreads / npartition_x;
		}

		parentAry = (Imgidx*)Malloc(imgsize * sizeof(Imgidx));
		isVisited = (_uint8*)Calloc((size_t)((imgsize)));
		isAvailable = (_uint8*)Malloc((size_t)(imgsize));
		set_isAvailable_par(isAvailable, npartition_x, npartition_y);

		_int64 blksz_x = width / npartition_x;
		_int64 blksz_y = height / npartition_y;
		_int64 blksz_xn = blksz_x + (width % npartition_x);
		_int64 blksz_yn = blksz_y + (height % npartition_y);
		_int64 numpartitions = (_int64)npartition_x * (_int64)npartition_y;

		Imgidx *startpidx = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blocksize = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *subtree_size = (Imgidx*)Calloc((numpartitions) * sizeof(Imgidx));
		Imgidx *subtree_start = (Imgidx*)Calloc((numpartitions + 1) * sizeof(Imgidx));
		Imgidx *subtree_nborderedges = (Imgidx*)Calloc((numpartitions) * sizeof(Imgidx));
		Imgidx *nrbnode = (Imgidx*)Calloc((numpartitions) * sizeof(Imgidx));
		Imgidx *subtree_cur = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *subtree_max = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blkws = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blkhs = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *dhist = (Imgidx*)Calloc((size_t)numlevels * (size_t)numpartitions * sizeof(Imgidx));
		char *blkflooddone = (char*)Calloc(numpartitions * sizeof(char));
		omp_lock_t *locks = (omp_lock_t*)Malloc((numpartitions + 1) * sizeof(omp_lock_t));

		omp_set_num_threads(_min(numpartitions, omp_get_num_procs()));

		for (p = 0;p < numpartitions + 1;p++)
			omp_init_lock(locks + p);

		p = q = 0;
		for (_int8 y = 0; y < npartition_y; y++)
		{
			q = y * width * (Imgidx)blksz_y;
			bool lastrow = (y == npartition_y - 1);
			Imgidx blkh = lastrow ? blksz_yn : blksz_y;
			for (_int8 x = 0; x < npartition_x; x++)
			{
				startpidx[p] = q + (Imgidx)x * (Imgidx)blksz_x;
				bool lastcol = (x == npartition_x - 1);
				Imgidx blkw = lastcol ? blksz_xn : blksz_x;
				blocksize[p] = blkh * blkw * 2 - (Imgidx)lastrow * blkw - (Imgidx)lastcol * blkh;
				blkws[p] = blkw;
				blkhs[p] = blkh;
				p++;
			}
		}

		#pragma omp parallel for private(p, q)
		for (int blk = 0; blk < numpartitions; blk++)
		{
			Imgidx bwidth = blkws[blk];
			Imgidx bheight = blkhs[blk];
			Imgidx *bhist = dhist + numlevels * blk;
			bool lastcol = (blk % npartition_x) == (npartition_x - 1);
			bool lastrow = (blk / npartition_x) == (npartition_y - 1);
			Imgidx spidx = startpidx[blk];
			Pixel maxdiff = 0;

			for (Imgidx i = 0;i < bheight;i++)
			{
				Pixel diff;
				p = spidx + i * width;
				bool blklastrow = (i == bheight - 1);
				for (Imgidx j = 0;j < bwidth - 1;j++)
				{
					q = p << 1;

					if (i < bheight - 1 || !lastrow)
					{
						diff = abs_diff(img[p], img[p + width]);
						dimg[q] = diff;
						bhist[diff]++;
						if (!blklastrow)
							maxdiff = _max(maxdiff, diff);
					}
					diff = abs_diff(img[p], img[p + 1]);
					dimg[q + 1] = diff;
					bhist[diff]++;
					maxdiff = _max(maxdiff, diff);
					p++;
				}
				q = p << 1;
				if (i < bheight - 1 || !lastrow)
				{
					diff = abs_diff(img[p], img[p + width]);
					dimg[q] = diff;
					bhist[diff]++;
					if (!blklastrow)
						maxdiff = _max(maxdiff, diff);
				}
				if (!lastcol)
				{
					diff = abs_diff(img[p], img[p + 1]);
					dimg[q + 1] = diff;
					bhist[diff]++;
					//maxdiff = _max(maxdiff, diff);
				}
			}
			bhist[maxdiff]++; //dummy for subroot
			subtree_max[blk] = maxdiff;


			Imgidx bhistsum = 0;
			for (int ii = 0;ii <= (int)maxdiff;ii++)
			{
				bhistsum += bhist[ii];
			}

			//Imgidx nrboderedges = (lastcol ? 0 : bheight) + (lastrow ? 0 : bwidth);
			Imgidx nrboderedges = (bheight + bwidth) * 2;
			subtree_nborderedges[blk] = nrboderedges;

			if (tse && imgsize > 1e5)
			{
				Imgidx est = TreeSizeEstimation(bhist, maxdiff + 1, bwidth * bheight, bhistsum, 0.5);
				if (est < 0)
				{
					printf("ERR: TSE yielded <0\n");
					Imgidx bhistsum = 0;
					for (int ii = 0;ii <= (int)maxdiff;ii++)
					{
						if (bhist[ii] < 0)
						{
							printf("bhist[%d] = %d\n", (int)ii, (int)bhist[ii]);
						}
						bhistsum += bhist[ii];
					}
				}
				//printf("blk %d: EST = %d\n", (int)blk, (int)est);
				subtree_size[blk] = est + nrboderedges;
			}
			else
				subtree_size[blk] = 1 + (1 + 2) * bwidth * bheight;//TreeSizeEstimation(bhist, maxdiff + 1, bwidth * bheight, bhistsum);

		}

		HierarQueue<Imgidx>** queues;
		queues = (HierarQueue<Imgidx>**)Calloc(numpartitions * sizeof(HierarQueue<Imgidx>*));
		for (int blk = 0;blk < numpartitions;blk++)
		{
			queues[blk] = new HierarQueue<Imgidx>((_uint64)blocksize[blk] + 1, (_int32)(subtree_max[blk] + 1));
		}

		double treesizemult_intv = 0.05;
		double treesizemult = 1.0 - treesizemult_intv;

		int flooddone = 0;
		while(!flooddone)
		{
			printf("Parallel flooding start -----------------------------------------\n");
			int numbusythr = 0;
			int outofmemory = 0;
			int numblkproc = 0;

			for (int blk = 0;blk < numpartitions;blk++)
				queues[blk]->reset_queue();
			for (int i = 0;i < imgsize;i++)
				isVisited[i] = 0;

			flooddone = 1;
			treesizemult = treesizemult + treesizemult_intv;
			maxSize = parflood_node_alloc(subtree_size, subtree_start, blkws, blkhs, numpartitions, treesizemult);
			curSize = maxSize;//doesnt do anything. only for debug

			for (p = 0;p < numpartitions;p++)
				omp_unset_lock(locks + p);
			omp_unset_lock(locks + numpartitions);
			numbusythr = 0;


			#pragma omp parallel for private(p, q)
			for (int blk = 0; blk < numpartitions; blk++)
			{
				if (outofmemory)
					continue;

				omp_set_lock(locks + blk);
				omp_set_lock(locks + numpartitions);
				numbusythr++;
				numblkproc++;
				omp_unset_lock(locks + numpartitions);

				Imgidx bwidth = blkws[blk];
				Imgidx bheight = blkhs[blk];
				Imgidx bareasum = bwidth * bheight;
				Imgidx *bhist = dhist + numlevels * blk;
				HierarQueue<Imgidx> *queue = queues[blk];
				Imgidx nidx = subtree_start[blk];
				Imgidx blkts = 0;
				int nidxblk = blk;
				Imgidx nidx_lim = subtree_start[blk + 1] - subtree_nborderedges[blk];//save room for nodes to be added in merge
				Imgidx iNode = 0;
				Pixel maxdiff = subtree_max[blk];

				queue->set_queue(bhist);

				Imgidx stack_top = nidx++;//imgsize + dimgsize + blk;
				Imgidx prev_top = stack_top;
				AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
				pNode->set(0, maxdiff, (double)0.0, (Pixel)-1, (Pixel)0, pix_type);
				pNode->parentidx = ROOTIDX;
				Pixel current_level = maxdiff;
				queue->push(startpidx[blk], current_level);
				if (outofmemory)
					continue;
				while (1) //flooding
				{
					while ((_int64)queue->min_level <= (_int64)current_level) //flood all levels below current_level
					{
						p = queue->pop();

						if (nidx < 0)
						{
							printf("thr%d: nidx < 0 (%d)\n", omp_get_thread_num(), (int)nidx);
						}
						//printf("thr%d: probing %d\n", omp_get_thread_num(), (int)p);
						if (is_visited(isVisited, p))
						{
							queue->find_minlev();
							continue;
						}

						isVisited[p] = 1;
						_uint8 isAv = isAvailable[p];

						if (connectivity == 4) _PUSH_NEIGHBORS_4N

						if ((_int64)current_level > (_int64)queue->min_level) //go to lower level
						{
							Pixel pix_val = img[p];
							current_level = queue->min_level;

							{
								if (nidx == nidx_lim)
								{
									if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
										break;
								}
								iNode = nidx++;
							}
							node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
							node[iNode].parentidx = stack_top;
							node[iNode].rootidx = ROOTIDX;
							stack_top = iNode;

							if (current_level)
							{
								{
									if (nidx == nidx_lim)
									{
										if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
											break;
									}
									iNode = nidx++;
								}
								node[iNode].copy(node + stack_top);
								node[iNode].alpha = 0;
								node[iNode].parentidx = stack_top;
								node[iNode].rootidx = ROOTIDX;
								prev_top = iNode;
							}
							parentAry[p] = iNode;
						}
						else
						{
							queue->find_minlev();

							if (current_level)
							{
								Pixel pix_val = img[p];
								{
									if (nidx == nidx_lim)
									{
										if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
											break;
									}
									iNode = nidx++;
								}
								node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
								node[stack_top].add(node + iNode, pix_type);
								node[iNode].parentidx = stack_top;
								node[iNode].rootidx = ROOTIDX;
								parentAry[p] = iNode;
							}
							else
							{
								parentAry[p] = stack_top;
								node[stack_top].add(img[p], pix_type);
							}
						}
						//if (stack_top == ROOTIDX) printf("SQEEEEAK\n");
					}

					if (outofmemory)
						break;

					remove_redundant_node(node, nidx, prev_top, stack_top);

					if (node[stack_top].area == bareasum)	// root node found...done
						break;

					//go to higher level
					iNode = node[stack_top].parentidx;
					if (iNode == ROOTIDX || (_int64)queue->min_level < (_int64)node[iNode].alpha) //new level from queue
					{
						//_CREATE_NEW_STACKTOP(queue->min_level)
						{
							if (nidx == nidx_lim)
							{
								if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
									break;
							}
							iNode = nidx++;
						}

						node[iNode].alpha = queue->min_level;
						node[iNode].copy(node + stack_top);
						node[iNode].parentidx = node[stack_top].parentidx;
						node[iNode].rootidx = ROOTIDX;
						node[stack_top].parentidx = iNode;
					}
					else //go to existing node
					{
						if (node[iNode].area == bareasum)	// root node found...done
							break;
						node[iNode].add(node + stack_top, pix_type);
					}

					if (node[iNode].area == bareasum)	// root node found...done
						break;

					prev_top = stack_top;
					stack_top = iNode;
					current_level = node[stack_top].alpha;
				}

				if (!outofmemory)
				{
					stack_top = (node[stack_top].area == bareasum) ? stack_top : iNode; //remove redundant root
					node[stack_top].parentidx = ROOTIDX;

					subtree_cur[nidxblk] = nidx;

					blkts += nidx - subtree_start[nidxblk];

					if (nidx == nidx_lim)//this should be really rare
					{
						blkflooddone[nidxblk] = 2; // flood done (at least for the native block), no free memory
					}
					else
					{
						blkflooddone[nidxblk] = 1; // flood done AND free memory available
						//if (blk != nidxblk)//not at the native block
						//printf("thr%d: subtree for blk %d: blkflooddone[%d] = 1\n", omp_get_thread_num(), (int)blk, (int)nidxblk);
					}

					omp_set_lock(locks + numpartitions);
					numbusythr--;
					//printf("th%d-- (%d/%d), blkflooddone[nidxblk] = %d\n", omp_get_thread_num(), numbusythr, omp_get_num_threads(), (int)blkflooddone[nidxblk]);
					omp_unset_lock(locks + numpartitions);
				}

				if (outofmemory)
				{
					printf("thr%d: subtree for blk %d: memory overflow - releasing lock %d\n", omp_get_thread_num(), (int)blk, (int)nidxblk);
					flooddone = 0;
				}
				else
				{
					//printf("thr%d: subtree for blk %d: releasing lock %d\n", omp_get_thread_num(), (int)blk, (int)nidxblk);
				}
				omp_unset_lock(locks + nidxblk);
				//printf("blk%d - area %d/%d\n", (int)blk, (int)(node[stack_top].area), (int)bareasum);
				//aa[omp_get_thread_num()]++;
			}//flood_end
		}

		rootidx = merge_subtrees(dimg, blksz_x, blksz_y, npartition_x, npartition_y, subtree_cur, 1, subtree_start, blkhs, blkws, nrbnode, subtree_nborderedges);

		Free(isVisited);
		Free(isAvailable);

		for (p = 0;p < numpartitions;p++)
			omp_destroy_lock(locks + p);
		Free(locks);
		Free(blkflooddone);

		for (int blk = 0;blk < numpartitions;blk++)
			delete queues[blk];
		Free(queues);
		Free(dimg);
		Free(startpidx);
		Free(blocksize);
		Free(subtree_start);
		Free(subtree_size);
		Free(subtree_nborderedges);
		Free(nrbnode);
		Free(subtree_cur);
		Free(subtree_max);
		Free(blkws);
		Free(blkhs);
		Free(dhist);
	}

	void Flood_Hierarqueue_cache_par_tse(Pixel *img, int numthreads, int tse)
	{
		if (sizeof(Pixel) > 2 || channel > 1)
		{
			printf("Error: Hierarchical queues do not work on >16 bits images or multispectral images\n");
			printf("Try Unionfind (algorithm code %d), flooding using Heapqueue (%d), trie queue (%d) or cached trie queue (%d) \n", UNIONFIND, FLOOD_HEAPQUEUE_CACHE, FLOOD_TRIE, FLOOD_TRIE_CACHE);
			return;
		}

		Imgidx imgsize, dimgsize, nredges;
		_int64 numlevels;
		Pixel *dimg;
		_uint8 *isVisited, *isAvailable;
		Imgidx p, q;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		numlevels = (sizeof(Pixel) == 1) ? 256 : 65536;

		double memuse_before = memuse;
		dimg = (Pixel*)Calloc((size_t)dimgsize * sizeof(Pixel));

		_int16 npartition_x, npartition_y;
		{
			_int16 optpart = 1;
			double optborderlength = (double)numthreads * (double)imgsize;
			for (int px = 2;px < numthreads; px++)
			{
				if (numthreads % px == 0)
				{
					int py = numthreads / px;

					if (((double)px * (double)height + (double)py * (double)width) < optborderlength)
					{
						optpart = px;
						optborderlength = ((double)px * (double)height + (double)py * (double)width);
					}
				}
			}
			npartition_x = (_int16)optpart;
			npartition_y = (_int16)numthreads / npartition_x;
		}

		parentAry = (Imgidx*)Malloc(imgsize * sizeof(Imgidx));
		isVisited = (_uint8*)Calloc((size_t)((imgsize)));
		isAvailable = (_uint8*)Malloc((size_t)(imgsize));
		set_isAvailable_par(isAvailable, npartition_x, npartition_y);

		_int64 blksz_x = width / npartition_x;
		_int64 blksz_y = height / npartition_y;
		_int64 blksz_xn = blksz_x + (width % npartition_x);
		_int64 blksz_yn = blksz_y + (height % npartition_y);
		_int64 numpartitions = (_int64)npartition_x * (_int64)npartition_y;

		Imgidx *startpidx = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blocksize = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *subtree_size = (Imgidx*)Calloc((numpartitions) * sizeof(Imgidx));
		Imgidx *subtree_start = (Imgidx*)Calloc((numpartitions + 1) * sizeof(Imgidx));
		Imgidx *subtree_nborderedges = (Imgidx*)Calloc((numpartitions) * sizeof(Imgidx));
		Imgidx *nrbnode = (Imgidx*)Calloc((numpartitions) * sizeof(Imgidx));
		Imgidx *subtree_cur = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *subtree_max = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blkws = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blkhs = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *dhist = (Imgidx*)Calloc((size_t)numlevels * (size_t)numpartitions * sizeof(Imgidx));
		char *blkflooddone = (char*)Calloc(numpartitions * sizeof(char));
		omp_lock_t *locks = (omp_lock_t*)Malloc((numpartitions + 1) * sizeof(omp_lock_t));

		omp_set_num_threads(_min(numpartitions, omp_get_num_procs()));

		for (p = 0;p < numpartitions + 1;p++)
			omp_init_lock(locks + p);

		p = q = 0;
		for (_int8 y = 0; y < npartition_y; y++)
		{
			q = y * width * (Imgidx)blksz_y;
			bool lastrow = (y == npartition_y - 1);
			Imgidx blkh = lastrow ? blksz_yn : blksz_y;
			for (_int8 x = 0; x < npartition_x; x++)
			{
				startpidx[p] = q + (Imgidx)x * (Imgidx)blksz_x;
				bool lastcol = (x == npartition_x - 1);
				Imgidx blkw = lastcol ? blksz_xn : blksz_x;
				blocksize[p] = blkh * blkw * 2 - (Imgidx)lastrow * blkw - (Imgidx)lastcol * blkh;
				blkws[p] = blkw;
				blkhs[p] = blkh;
				p++;
			}
		}

		#pragma omp parallel for private(p, q)
		for (int blk = 0; blk < numpartitions; blk++)
		{
			Imgidx bwidth = blkws[blk];
			Imgidx bheight = blkhs[blk];
			Imgidx *bhist = dhist + numlevels * blk;
			bool lastcol = (blk % npartition_x) == (npartition_x - 1);
			bool lastrow = (blk / npartition_x) == (npartition_y - 1);
			Imgidx spidx = startpidx[blk];
			Pixel maxdiff = 0;

			for (Imgidx i = 0;i < bheight;i++)
			{
				Pixel diff;
				p = spidx + i * width;
				bool blklastrow = (i == bheight - 1);
				for (Imgidx j = 0;j < bwidth - 1;j++)
				{
					q = p << 1;

					if (i < bheight - 1 || !lastrow)
					{
						diff = abs_diff(img[p], img[p + width]);
						dimg[q] = diff;
						bhist[diff]++;
						if (!blklastrow)
							maxdiff = _max(maxdiff, diff);
					}
					diff = abs_diff(img[p], img[p + 1]);
					dimg[q + 1] = diff;
					bhist[diff]++;
					maxdiff = _max(maxdiff, diff);
					p++;
				}
				q = p << 1;
				if (i < bheight - 1 || !lastrow)
				{
					diff = abs_diff(img[p], img[p + width]);
					dimg[q] = diff;
					bhist[diff]++;
					if (!blklastrow)
						maxdiff = _max(maxdiff, diff);
				}
				if (!lastcol)
				{
					diff = abs_diff(img[p], img[p + 1]);
					dimg[q + 1] = diff;
					bhist[diff]++;
					//maxdiff = _max(maxdiff, diff);
				}
			}
			bhist[maxdiff]++; //dummy for subroot
			subtree_max[blk] = maxdiff;


			Imgidx bhistsum = 0;
			for (int ii = 0;ii <= (int)maxdiff;ii++)
			{
				bhistsum += bhist[ii];
			}

			//Imgidx nrboderedges = (lastcol ? 0 : bheight) + (lastrow ? 0 : bwidth);
			Imgidx nrboderedges = (bheight + bwidth) * 2;
			subtree_nborderedges[blk] = nrboderedges;

			if (tse && imgsize > 1e5)
			{
				Imgidx est = TreeSizeEstimation(bhist, maxdiff + 1, bwidth * bheight, bhistsum, 0.5);
				if (est < 0)
				{
					printf("ERR: TSE yielded <0\n");
					Imgidx bhistsum = 0;
					for (int ii = 0;ii <= (int)maxdiff;ii++)
					{
						if (bhist[ii] < 0)
						{
							printf("bhist[%d] = %d\n", (int)ii, (int)bhist[ii]);
						}
						bhistsum += bhist[ii];
					}
				}
				//printf("blk %d: EST = %d\n", (int)blk, (int)est);
				subtree_size[blk] = est + nrboderedges;
			}
			else
				subtree_size[blk] = 1 + (1 + 2) * bwidth * bheight;//TreeSizeEstimation(bhist, maxdiff + 1, bwidth * bheight, bhistsum);
		}

		memuse_before = (double)memuse;

		HierarQueueCache<Imgidx,Pixel>** queues;
		queues = (HierarQueueCache<Imgidx,Pixel>**)Calloc(numpartitions * sizeof(HierarQueueCache<Imgidx,Pixel>*));
		for (int blk = 0;blk < numpartitions;blk++)
		{
			queues[blk] = new HierarQueueCache<Imgidx,Pixel>((_uint64)blocksize[blk] + 1, (_int32)(subtree_max[blk] + 1));
		}

		double treesizemult_intv = 0.05;
		double treesizemult = 1.0 - treesizemult_intv;

		int flooddone = 0;
		while(!flooddone)
		{
			printf("Parallel flooding start -----------------------------------------\n");
			int numbusythr = 0;
			int outofmemory = 0;
			int numblkproc = 0;

			for (int blk = 0;blk < numpartitions;blk++)
				queues[blk]->reset_queue();
			//	delete queues[blk];
			//Free(queues);
			//queues = (HierarQueue<Imgidx>**)Calloc(numpartitions * sizeof(HierarQueue<Imgidx>*));
			//for (int blk = 0;blk < numpartitions;blk++)
			//	queues[blk] = new HierarQueue<Imgidx>((_uint64)blocksize[blk] + 1);
			for (int i = 0;i < imgsize;i++)
				isVisited[i] = 0;

			flooddone = 1;
			treesizemult = treesizemult + treesizemult_intv;
			maxSize = parflood_node_alloc(subtree_size, subtree_start, blkws, blkhs, numpartitions, treesizemult);
			curSize = maxSize;//doesnt do anything. only for debug

			//for (int blkblk = 0; blkblk < numpartitions; blkblk++)
			{
			//	printf("blk %d: ESTSUB = %d\n", (int)blkblk, (int)subtree_start[blkblk + 1] - (int)subtree_start[blkblk] - (int)subtree_nborderedges[blkblk]);
			}

			for (p = 0;p < numpartitions;p++)
				omp_unset_lock(locks + p);
			omp_unset_lock(locks + numpartitions);
			numbusythr = 0;


			#pragma omp parallel for private(p, q)
			for (int blk = 0; blk < numpartitions; blk++)
			{
				if (outofmemory)
					continue;

				omp_set_lock(locks + blk);
				omp_set_lock(locks + numpartitions);
				numbusythr++;
				numblkproc++;
				omp_unset_lock(locks + numpartitions);

				Imgidx bwidth = blkws[blk];
				Imgidx bheight = blkhs[blk];
				Imgidx bareasum = bwidth * bheight;
				Imgidx *bhist = dhist + numlevels * blk;
				HierarQueueCache<Imgidx,Pixel> *queue = queues[blk];
				Imgidx nidx = subtree_start[blk], startingidx = nidx;
				Imgidx blkts = 0;
				int nidxblk = blk;
				Imgidx nidx_lim = subtree_start[blk + 1] - subtree_nborderedges[blk];//save room for nodes to be added in merge
				Imgidx iNode = 0;
				Pixel maxdiff = subtree_max[blk];


				queue->set_queue(bhist);

				Imgidx stack_top = nidx++;//imgsize + dimgsize + blk;
				Imgidx prev_top = stack_top;
				AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
				pNode->set(0, maxdiff, (double)0.0, (Pixel)-1, (Pixel)0, pix_type);
				pNode->parentidx = ROOTIDX;
				Pixel current_level = maxdiff;
				queue->push(startpidx[blk], current_level);
				if (outofmemory)
					continue;
				while (1) //flooding
				{
					while ((_int64)queue->min_level <= (_int64)current_level) //flood all levels below current_level
					{
						p = queue->pop();

						if (nidx < 0)
						{
							printf("thr%d: nidx < 0 (%d)\n", omp_get_thread_num(), (int)nidx);
						}

						if (is_visited(isVisited, p))
						{
							queue->find_minlev();
							continue;
						}

						isVisited[p] = 1;
						_uint8 isAv = isAvailable[p];

						if (connectivity == 4) _PUSH_NEIGHBORS_4N

						if ((_int64)current_level > (_int64)queue->min_level) //go to lower level
						{
							Pixel pix_val = img[p];
							current_level = queue->min_level;

							{
								if (nidx == nidx_lim)
								{
									if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
										break;
								}
								iNode = nidx++;
							}
							node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
							node[iNode].parentidx = stack_top;
							node[iNode].rootidx = ROOTIDX;
							stack_top = iNode;

							if (current_level)
							{
								{
									if (nidx == nidx_lim)
									{
										if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
											break;
									}
									iNode = nidx++;
								}
								node[iNode].copy(node + stack_top);
								node[iNode].alpha = 0;
								node[iNode].parentidx = stack_top;
								node[iNode].rootidx = ROOTIDX;
								prev_top = iNode;
							}
							parentAry[p] = iNode;
						}
						else
						{
							queue->find_minlev();

							if (current_level)
							{
								Pixel pix_val = img[p];
								{
									if (nidx == nidx_lim)
									{
										if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
											break;
									}
									iNode = nidx++;
								}
								node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
								node[stack_top].add(node + iNode, pix_type);
								node[iNode].parentidx = stack_top;
								node[iNode].rootidx = ROOTIDX;
								parentAry[p] = iNode;
							}
							else
							{
								parentAry[p] = stack_top;
								node[stack_top].add(img[p], pix_type);
							}
						}
						//if (stack_top == ROOTIDX) printf("SQEEEEAK\n");
					}

					if (outofmemory)
						break;

					remove_redundant_node(node, nidx, prev_top, stack_top);

					if (node[stack_top].area == bareasum)	// root node found...done
						break;

					//go to higher level
					iNode = node[stack_top].parentidx;
					if (iNode == ROOTIDX || (_int64)queue->min_level < (_int64)node[iNode].alpha) //new level from queue
					{
						//_CREATE_NEW_STACKTOP(queue->min_level)
						{
							if (nidx == nidx_lim)
							{
								if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
									break;
							}
							iNode = nidx++;
						}

						node[iNode].alpha = queue->min_level;
						node[iNode].copy(node + stack_top);
						node[iNode].parentidx = node[stack_top].parentidx;
						node[iNode].rootidx = ROOTIDX;
						node[stack_top].parentidx = iNode;
					}
					else //go to existing node
					{
						if (node[iNode].area == bareasum)	// root node found...done
							break;
						node[iNode].add(node + stack_top, pix_type);
					}

					if (node[iNode].area == bareasum)	// root node found...done
						break;

					prev_top = stack_top;
					stack_top = iNode;
					current_level = node[stack_top].alpha;
				}

				if (!outofmemory)
				{
					stack_top = (node[stack_top].area == bareasum) ? stack_top : iNode; //remove redundant root
					node[stack_top].parentidx = ROOTIDX;

					subtree_cur[nidxblk] = nidx;

					blkts += nidx - subtree_start[nidxblk];

					if (nidx == nidx_lim)//this should be really rare
					{
						blkflooddone[nidxblk] = 2; // flood done (at least for the native block), no free memory
					}
					else
					{
						blkflooddone[nidxblk] = 1; // flood done AND free memory available
					}

					omp_set_lock(locks + numpartitions);
					numbusythr--;
					omp_unset_lock(locks + numpartitions);
				}

				if (outofmemory)
				{

					printf("thr%d: subtree for blk %d: memory overflow - releasing lock %d\n", omp_get_thread_num(), (int)blk, (int)nidxblk);
					flooddone = 0;
				}
				else
				{
					//printf("thr%d: subtree for blk %d: releasing lock %d\n", omp_get_thread_num(), (int)blk, (int)nidxblk);
				}
				omp_unset_lock(locks + nidxblk);
			}//flood_end
		}

		rootidx = merge_subtrees(dimg, blksz_x, blksz_y, npartition_x, npartition_y, subtree_cur, 1, subtree_start, blkhs, blkws, nrbnode, subtree_nborderedges);

		//tmp
		Free(isVisited);
		Free(isAvailable);

		for (p = 0;p < numpartitions;p++)
			omp_destroy_lock(locks + p);
		Free(locks);
		Free(blkflooddone);

		for (int blk = 0;blk < numpartitions;blk++)
			delete queues[blk];
		Free(queues);
		Free(dimg);
		Free(startpidx);
		Free(blocksize);
		Free(subtree_start);
		Free(subtree_size);
		Free(subtree_nborderedges);
		Free(nrbnode);
		Free(subtree_cur);
		Free(subtree_max);
		Free(blkws);
		Free(blkhs);
		Free(dhist);
	}

	void quantize_dimg(_uint8 *qimg, Pixel *dimg, _int64 dimgsize, _int64 binsize)
	{
		#pragma omp parallel for
		for (Imgidx i = 0;i < dimgsize;i++)
			qimg[i] = QUANTIZE_RANK(dimg[i],binsize);
	}

	//Find subtree root and do path compression
	inline Imgidx find_root(Imgidx p)
	{
		if (p == ROOTIDX)
			return ROOTIDX;

		Imgidx r, q;

		for (r = p;node[r].rootidx != ROOTIDX;r = node[r].rootidx)
			;

		while(p != r)
		{
			q = node[p].rootidx;
			node[p].rootidx = r;
			p = q;
		}

		return r;
	}

	inline Imgidx find_root_in(Imgidx p)
	{
		if (p == ROOTIDX)
			return ROOTIDX;

		Imgidx r, q;

		for (r = p;node_in[r].rootidx != ROOTIDX;r = node_in[r].rootidx)
				;

		while(p != r)
		{
			q = node_in[p].rootidx;
			node_in[p].rootidx = r;
			p = q;
			//int2++;
		}

		return r;
	}

	void Unionfind(Pixel* img)
	{
		Imgidx imgsize, nredges;
		RankItem<Imgidx, double> *rankitem, *pRank;

		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		curSize = maxSize = imgsize + nredges; //to be compatible with flooding algorithms

		omp_set_num_threads(1);
		rankitem = (RankItem<Imgidx, double>*)Malloc(nredges * sizeof(RankItem<Imgidx, double>));
		compute_difference_and_sort(rankitem, img, nredges);

		//initialize_node(img, rankitem, maxpixval);
		maxSize = imgsize + nredges;
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		for (Imgidx p = 0; p < imgsize;p++)
		{
			node[p].set(1, 0, (double)img[p], img[p], img[p]);
			node[p].rootidx = node[p].parentidx = ROOTIDX;
		}

		bool unionbyrank = 0;//(sizeof(Pixel) <= 2);
		Imgidx *treedepth = 0;

		if (unionbyrank) treedepth = (Imgidx*)Calloc(maxSize * sizeof(Imgidx));

		Imgidx curSize = 0;
		for (Imgidx r = 0; r < nredges; r++)
		{
			pRank = rankitem + r;

			Imgidx x, x0;
			Imgidx y, y0;
			Imgidx z;

			Imgidx nodeaddr = curSize + imgsize;

			x0 = pRank->get_pidx0(connectivity);
			y0 = pRank->get_pidx1(width,connectivity);
			x = find_root(x0);
			y = find_root(y0);

			if (x == y) //already connected, nothing to do
				continue;

			if (x < y) {z = x;x = y;y = z;}

			if (!unionbyrank || node[x].alpha != pRank->alpha)
			{
				curSize++;
				//node_in[r].set(0, pRank->alpha, 0.0, maxpixval, 0, pix_type);
				node[nodeaddr].copy(node + x);
				node[nodeaddr].alpha = pRank->alpha;
				node[nodeaddr].parentidx = node[nodeaddr].rootidx = ROOTIDX;
	//			node_in[r].parentidx = node[x].parentidx;
				node[x].parentidx = node[x].rootidx = nodeaddr;
				node[nodeaddr].add(node + y, pix_type);
				node[y].parentidx = node[y].rootidx = nodeaddr;
				if (unionbyrank) treedepth[nodeaddr] = _max(treedepth[x],treedepth[y]) + 1;

				if (node[nodeaddr].area == imgsize)
				{
					rootidx = nodeaddr;
					break;
				}
			}
			else
			{
				if (node[x].alpha == node[y].alpha && treedepth[x] < treedepth[y]) {z = x;x = y;y = z;}
				node[x].add(node + y, pix_type);
				node[y].parentidx = node[y].rootidx = x;
				treedepth[x] = _max(treedepth[x],treedepth[y] + 1);

				if (node[x].area == imgsize)
				{
					rootidx = x;
					break;
				}
			}
		}
		if (unionbyrank) Free(treedepth);
		Free(rankitem);
//		Free(isVisited);
//		Free(isAvailable);
	}

	void blockwise_tse(Imgidx *subtree_size, Imgidx *subtree_nborderedges, double *nrmsds, Imgidx *dhist, Imgidx *subtree_max, Imgidx *blkws, Imgidx *blkhs, _int8 npartition_x, _int8 npartition_y, Imgidx numbins)
	{
		_int64 numpartitions = (_int64)npartition_x * (_int64)npartition_y;
		Pixel maxdiff;

		for (int blk = 0; blk < numpartitions; blk++)
		{
			maxdiff = subtree_max[blk];
			Imgidx *bhist = dhist + numbins * blk;
			Imgidx bhistsum = 0;
			Imgidx bwidth = blkws[blk];
			Imgidx bheight = blkhs[blk];
			//bool lastcol = (blk % npartition_x) == (npartition_x - 1);
			//bool lastrow = (blk / npartition_x) == (npartition_y - 1);

			for (int ii = 0;ii <= (int)maxdiff;ii++)
				bhistsum += bhist[ii];

			double nrmsd_blk = 0;
			for (int ii = 0; ii < (int)(maxdiff + 1); ii++)
				nrmsd_blk += ((double)dhist[ii]) * ((double)dhist[ii]);
			nrmsd_blk = sqrt((nrmsd_blk - (double)bhistsum) / ((double)bhistsum * ((double)bhistsum - 1.0)));
			nrmsd_blk = ((A * exp(SIGMA * nrmsd_blk) + B) + M);
			nrmsds[blk] = nrmsd_blk;

			//Imgidx nrboderedges = (lastcol ? 0 : bheight) + (lastrow ? 0 : bwidth);
			Imgidx nrboderedges = 2 * (bheight + bwidth);
			subtree_nborderedges[blk] = nrboderedges;
			subtree_size[blk] = TreeSizeEstimation(bhist, maxdiff + 1, bwidth * bheight, bhistsum, 0.5) + nrboderedges;
	//subtree_size[blk] = TreeSizeEstimation(bhist, maxdiff + 1, bwidth * bheight, bhistsum, 0.5); //no need to add room for borders (already computed in quantization)
			//printf("Subblock %d: bhist sum = %d, size estimate: %d\n", (int)blk, (int)bhistsum, (int)subtree_size[blk]);
		}
	}

	//blockwise quantization and histogram computation (for pilot_rank)
	void quantize_ranks_compute_histogram(_uint8 *qrank, Imgidx* rank, Pixel* img, Imgidx *dhist, Imgidx *blkws, Imgidx *blkhs,
		Imgidx *startpidx, _int64 binsize, Imgidx numbins, _int8 npartition_x, _int8 npartition_y, Imgidx* subtree_max)
	{
		_int64 numpartitions = (_int64)npartition_x * (_int64)npartition_y;

		#pragma omp parallel for
		for (int blk = 0; blk < numpartitions; blk++)
		{
			Imgidx p, q;
			Imgidx bwidth = blkws[blk];
			Imgidx bheight = blkhs[blk];
			Imgidx *bhist = dhist + numbins * blk;
			bool lastcol = (blk % npartition_x) == (npartition_x - 1);
			bool lastrow = (blk / npartition_x) == (npartition_y - 1);
			Imgidx spidx = startpidx[blk];
			Pixel maxdiff = 0;

			for (Imgidx i = 0;i < bheight;i++)
			{
				Imgidx r;
				_uint8 qr;
				p = spidx + i * width;
				bool blklastrow = (i == bheight - 1);
				for (Imgidx j = 0;j < bwidth - 1;j++)
				{
					q = p << 1;
					if (i < bheight - 1 || !lastrow)
					{
						r = rank[q];
						qr = QUANTIZE_RANK(r, binsize);
						qrank[q] = qr;
						bhist[qr]++;
						if (!blklastrow)
							maxdiff = _max(maxdiff, qr);
					}
					r = rank[q + 1];
					qr = QUANTIZE_RANK(r, binsize);
					qrank[q + 1] = qr;
					bhist[qr]++;
					maxdiff = _max(maxdiff, qr);
					p++;
				}

				q = p << 1;
				if (i < bheight - 1 || !lastrow)
				{
					r = rank[q];
					qr = QUANTIZE_RANK(r, binsize);
					qrank[q] = qr;
					bhist[qr]++;
					if (!blklastrow)
						maxdiff = _max(maxdiff, qr);
				}
				if (!lastcol)
				{
					r = rank[q + 1];
					qr = QUANTIZE_RANK(r, binsize);
					qrank[q + 1] = qr;
					bhist[qr]++;
				}
			}
			bhist[maxdiff]++; //dummy for subroot
			subtree_max[blk] = maxdiff;
		}

		//print_value_edge("qranks",qrank,0);
	}

	void quantize_ranks(_uint8 *qrank, Imgidx *rank, _int64 dimgsize, _int64 binsize)
	{
		#pragma omp parallel for
		for (Imgidx i = 0;i < dimgsize;i++)
			qrank[i] = QUANTIZE_RANK(rank[i],binsize);
	}

	inline _uint8 pow_quantization(Imgidx rank, _uint64 qint)
	{
		return (_uint8)(((double)rank * (double)rank) / (double)qint);
	}

	void pow_quantize_ranks(_uint8 *qrank, Imgidx *rank, _int64 dimgsize, _int64 qint)
	{
		#pragma omp parallel for
		for (Imgidx i = 0;i < dimgsize;i++)
			qrank[i] = pow_quantization(rank[i],qint);
	}

	inline Imgidx find_root(AlphaNode<Imgidx, Pixel> *pilottree, Imgidx p, Pixel below_this_qlevel)
	{
		Imgidx q = parentAry[p], r;

		//int cnt = 0;

		while(pilottree[(r = pilottree[q].parentidx)].alpha < below_this_qlevel) //fix qlevel (also sum)
		{
			q = r;
		}

		return q;
	}

	Imgidx descendroots(Imgidx q, _int64 qlevel, AlphaNode<Imgidx, Pixel> *pilottree)
	{
		Imgidx c = pilottree[q].parentidx;
		while((_int64)pilottree[c].alpha < qlevel)
		{
			//int1++;
			q = c;
			c = pilottree[c].parentidx;
		}
		return q;
	}

	//Hybrid_Pilot_Rank
	void unionfind_refine_qlevel(_int64 qlevel, _int64 binsize, Imgidx nredges, AlphaNode<Imgidx, Pixel> *pilottree, RankItem<Imgidx, double> *rankitem, _int8 *redundant_edge, Index *rank2rankitem)
	{
		RankItem<Imgidx, double> *pRank;
		Imgidx rank_start = (Imgidx)qlevel * (Imgidx)binsize;
		Imgidx rank_end = min(nredges - 1,((Imgidx)qlevel + 1) * (Imgidx)binsize - 1);
		Imgidx imgsize = width * height;

		if (qlevel)
		{
			for (Imgidx r = rank_start; r <= rank_end; r++)
			{
				Imgidx ridx = r;
				if (rank2rankitem)
					ridx = (Imgidx)rank2rankitem[r];

				pRank = rankitem + ridx;

				//look for ancestor at qlevel
				Imgidx x, x0;
				Imgidx y, y0;

				Imgidx nodeaddr = r + imgsize;

				if (redundant_edge[rankitem[ridx].dimgidx])
				{
					//printf("[skip]qlevel %d / rank %d(%d) - refining %d-%d\n", (int)qlevel, (int)r, (int)pRank->alpha, (int)pRank->get_pidx0(connectivity), (int)pRank->get_pidx1(width,connectivity));
					continue;
				}

				//printf("qlevel %d / rank %d(%d) - refining %d-%d\n", (int)qlevel, (int)r, (int)pRank->alpha, (int)pRank->get_pidx0(connectivity), (int)pRank->get_pidx1(width,connectivity));

				//non-zero level
				//if (qlevel)
				{
					x0 = find_root(pilottree, pRank->get_pidx0(connectivity), qlevel);
					y0 = find_root(pilottree, pRank->get_pidx1(width,connectivity), qlevel);

					if (x0 == y0) //already connected, nothing to do
					{
						continue;
					}

					x = find_root(pilottree[x0].rootidx);
					y = find_root(pilottree[y0].rootidx);
				}



				if ((x!= ROOTIDX) && (x == y)) //already connected, nothing to do
				{
					continue;
				}

				//add the new node to the refined tree
				{
					if (x == ROOTIDX)
					{
						{
							node_in[r].copy(pilottree + x0);
							node_in[r].parentidx = pilottree[x0].parentidx;
							pilottree[x0].rootidx = nodeaddr;
						}
					}
					else
					{
						node_in[r].copy(node + x);
						node_in[r].parentidx = node[x].parentidx;
						node[x].parentidx = node[x].rootidx = nodeaddr;
					}

					//attach to y
					if (y == ROOTIDX)
					{
						//if (qlevel)
						{
							node_in[r].add(pilottree + y0, pix_type);
							pilottree[y0].rootidx = nodeaddr;
						}
					}
					else
					{
						node_in[r].add(node + y, pix_type);
						node[y].parentidx = node[y].rootidx = nodeaddr;
					}
				}
			}
		}
		else
		{
			for (Imgidx r = rank_start; r <= rank_end; r++)
			{
				if (rank2rankitem)
					pRank = rankitem + rank2rankitem[r];
				else
					pRank = rankitem + r;


				//look for ancestor at qlevel
				Imgidx x, x0;
				Imgidx y, y0;

				Imgidx nodeaddr = r + imgsize;

				//else
				{
					//find subtree roots from two edge incidents
					x0 = pRank->get_pidx0(connectivity);
					y0 = pRank->get_pidx1(width,connectivity);

					x = find_root(node[x0].rootidx);
					y = find_root(node[y0].rootidx);
				}


				if ((x!= ROOTIDX) && (x == y)) //already connected, nothing to do
					continue;

				{
					//attach to x
					if (x == ROOTIDX)
					{
						//else
						{
							node_in[r].copy(node + x0);
							node_in[r].parentidx = parentAry[x0];
							node[x0].parentidx = node[x0].rootidx = nodeaddr;
						}
					}
					else
					{
						node_in[r].copy(node + x);
						node_in[r].parentidx = node[x].parentidx;
						node[x].parentidx = node[x].rootidx = nodeaddr;
					}

					//attach to y
					if (y == ROOTIDX)
					{
						//else
						{
							node_in[r].add(node + y0, pix_type);
							node[y0].parentidx = node[y0].rootidx = nodeaddr;
						}
					}
					else
					{
						node_in[r].add(node + y, pix_type);
						node[y].parentidx = node[y].rootidx = nodeaddr;
					}
				}
			}
		}


	}

	//compute edge histogram of the quantized rank image.
	//In the hypergraph implementation, edges on the subblock borders are also counted.
	void compute_dhist_par(_uint8 *qrank, Imgidx *dhist, Imgidx *startpidx, _int32 numbins, _int8 npartition_x, _int8 npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn)
	{
		if (connectivity == 4)
		{
			//gdhist = edge histogram for edges on the subblock borders
			Imgidx *gdhist = dhist + (int)npartition_x * (int)npartition_y * numbins;
			//Imgidx blk = 0;
			//for (Imgidx y = 0;y < (Imgidx)npartition_y;y++)
			//for (Imgidx x = 0;x < (Imgidx)npartition_x;x++)
		  #pragma omp parallel for
			for (Imgidx blk = 0;blk < (int)npartition_x * (int)npartition_y;blk++)
			{
				Imgidx x = blk % npartition_x;
				Imgidx y = blk / npartition_y;
				Imgidx lastcol = (x == npartition_x - 1);
				Imgidx lastrow = (y == npartition_y - 1);
				Imgidx xn = lastcol ? blksz_xn : blksz_x;
				Imgidx yn = lastrow ? blksz_yn : blksz_y;
				Imgidx p0 = startpidx[blk] << 1, p;
				Imgidx *pdhist = dhist + blk * numbins;

				//Imgidx maxval = qrank[p0];
				//Imgidx maxpidx;

				for (Imgidx i = 0;i < yn - 1;i++)
				{
					p = p0 + width * (i << 1);
					for (Imgidx j = 0;j < xn - 1;j++)
					{
						pdhist[qrank[p++]]++;
						pdhist[qrank[p++]]++;
					}
					pdhist[qrank[p++]]++;
					if (!lastcol)
						gdhist[qrank[p++]]++;
				}

				//the last row of the subblock is...
				if (lastrow)
				{//the last row of the image
					p = p0 + width * ((yn - 1) << 1) + 1;
					for (Imgidx j = 0;j < xn - 1;j++)
					{
						pdhist[qrank[p]]++;
						p += 2;
					}
					if (!lastcol)
						gdhist[qrank[p++]]++;
				}
				else
				{
					p = p0 + width * ((yn - 1) << 1);
					for (Imgidx j = 0;j < xn - 1;j++)
					{
						gdhist[qrank[p++]]++;
						pdhist[qrank[p++]]++;
					}
					gdhist[qrank[p++]]++;
					if (!lastcol)
						gdhist[qrank[p]]++;
				}
				//blk++;
			}
		}
		else
		{
			//L A T e r
		}
	}

	void compute_dhist_par_hypergraph(_uint8 *qrank, Imgidx *dhist, Imgidx *startpidx, _int32 numbins, _int8 npartition_x, _int8 npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn, Imgidx *blkmaxpidx)
	{
		if (connectivity == 4)
		{
			#pragma omp parallel for
			for (Imgidx blk = 0;blk < (int)npartition_x * (int)npartition_y;blk++)
			{
				Imgidx x = blk % npartition_x;
				Imgidx y = blk / npartition_y;
				Imgidx lastcol = (x == npartition_x - 1);
				Imgidx lastrow = (y == npartition_y - 1);
				Imgidx xn = lastcol ? blksz_xn : blksz_x;
				Imgidx yn = lastrow ? blksz_yn : blksz_y;
				Imgidx p0 = startpidx[blk] << 1, p;
				Imgidx *pdhist = dhist + blk * numbins;

				Imgidx maxval = qrank[p0];
				Imgidx maxpidx = p0;

				for (Imgidx i = 0;i < yn - 1;i++) //for subimage rows (except for the last)
				{
					p = p0 + width * (i << 1);
					for (Imgidx j = 0;j < xn - 1;j++) //for subimage cols (except for the last)
					{
						if (maxval < qrank[p]) {maxval = qrank[p]; maxpidx = p;}
						pdhist[qrank[p++]]++;
						if (maxval < qrank[p]) {maxval = qrank[p]; maxpidx = p;}
						pdhist[qrank[p++]]++;
					}
					if (maxval < qrank[p]) {maxval = qrank[p]; maxpidx = p;}
					pdhist[qrank[p++]]++;
					if (!lastcol)
					{
						if (maxval < qrank[p]) {maxval = qrank[p]; maxpidx = p;}
						pdhist[qrank[p++]]++;
					}
				}

				//the last row of the subblock is...
				if (lastrow)
				{//the last row of the image
					p = p0 + width * ((yn - 1) << 1) + 1;
					for (Imgidx j = 0;j < xn - 1;j++)
					{
						if (maxval < qrank[p]) {maxval = qrank[p]; maxpidx = p;}
						pdhist[qrank[p]]++;
						p += 2;
					}
					if (!lastcol)
					{
						if (maxval < qrank[p]) {maxval = qrank[p]; maxpidx = p;}
						pdhist[qrank[p++]]++;
					}
				}
				else
				{
					p = p0 + width * ((yn - 1) << 1);
					for (Imgidx j = 0;j < xn - 1;j++)
					{
						if (maxval < qrank[p]) {maxval = qrank[p]; maxpidx = p;}
						pdhist[qrank[p++]]++;
						if (maxval < qrank[p]) {maxval = qrank[p]; maxpidx = p;}
						pdhist[qrank[p++]]++;
					}

					if (maxval < qrank[p]) {maxval = qrank[p]; maxpidx = p;}
					pdhist[qrank[p++]]++;
					if (!lastcol)
					{
						if (maxval < qrank[p]) {maxval = qrank[p]; maxpidx = p;}
						pdhist[qrank[p]]++;
					}
				}

				blkmaxpidx[blk] = maxpidx;
				//blk++;
			}
		}
		else
		{
			//L A T e r
		}
	}

	void create_queues(HierarQueue<Imgidx> ***queues, Imgidx *dhist, _int8 npartition_x, _int8 npartition_y, _int32 numbins)
	{
			*queues = new HierarQueue<Imgidx>*[(int)npartition_x * (int)npartition_y];
			HierarQueue<Imgidx> **Q = *queues;

			//queues = new Trie_Cache<Imgidx, trieidx>*[npartition * npartition];//(nredges, listsize);
			Imgidx p = 0;
			Imgidx *dh = dhist;
			for (Imgidx y = 0;y < npartition_y; y++)
			{
				for (Imgidx x = 0;x < npartition_x; x++)
				{
					Q[p++] = new HierarQueue<Imgidx>(dh, numbins);
					//Q[p++]->push(nredges);
					dh += numbins;
				}
			}
	}

	//obsolete code (subtree nodes are now indexed based on their level)
	void fix_subtreeidx(Imgidx *subtreestart, Imgidx *startpidx, Imgidx *cursizes, _int8 npartition_x, _int8 npartition_y, int numpartitions, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn)
	{
		for (int b = 0;b < numpartitions;b++)
		{
			Imgidx poffset = subtreestart[b];

			//fix parentAry
			Imgidx bx = ((b % npartition_x) == npartition_x - 1) ? blksz_xn : blksz_x;
			Imgidx by = ((b / npartition_x) == npartition_y - 1) ? blksz_yn : blksz_y;
			Imgidx pidx = startpidx[b];
			for (Imgidx p = 0;p < by;p++)
			{
				for (Imgidx q = 0;q < bx;q++)
				{
					parentAry[pidx++] += poffset;
				}
				pidx += width - bx;
			}

			//fix node parentidxs
			AlphaNode<Imgidx, Pixel> *ptree = node + poffset;
			for (Imgidx p = 0;p < cursizes[b];p++)
				if (ptree[p].parentidx != ROOTIDX)
					ptree[p].parentidx += poffset;
		}
	}

	//The one with the pilottree indicing (slow?)
	void merge_subtrees(_uint8 *qrank, Imgidx *qindex, _int64 blksz_x, _int64 blksz_y, Imgidx neighbor_offset, Imgidx shamt, Imgidx npartition_x, Imgidx npartition_y, _int32 numbins)
	{
		Imgidx x, y, r, p, q, dimgidx;
		Imgidx imgsize = height * width;

		//merging border(hor)
		for (y = blksz_y;y <= blksz_y * (npartition_y - 1); y += blksz_y)
		{
			q = y * width;
			for (p = q - width;p < q;p++)
			{
				dimgidx = (p << shamt) + neighbor_offset;
				r = qrank[dimgidx];

				connect(parentAry[p], parentAry[p + width], (Pixel)r, (Imgidx)qindex[r]++);
			}
		}
		neighbor_offset = (connectivity == 4) ? 1 : 3;

		//merging border(ver)
		for (x = blksz_x - 1;x <= blksz_x * (npartition_x - 1); x += blksz_x)
		{
			q = x + imgsize;
			for (p = x;p < q;p += width)
			{
				dimgidx = (p << shamt) + neighbor_offset;
				r = qrank[dimgidx];

				connect(parentAry[p], parentAry[p + 1], (Pixel)r, (Imgidx)qindex[r]++);
			}
		}

		//reset rootidxs
		for (p = 0;p < maxSize;p++)
			if (node[p].area) node[p].rootidx = ROOTIDX;

		//Look for the root
		for (p = parentAry[0];node[p].area != imgsize;p = node[p].parentidx)
			;
		node[p].rootidx = node[p].parentidx = ROOTIDX;
		rootidx = p;
	}

	void merge_subtrees(_uint8 *qrank, _int64 blksz_x, _int64 blksz_y, Imgidx neighbor_offset, Imgidx shamt, Imgidx npartition_x, Imgidx npartition_y, _int32 numbins)
	{
		Imgidx x, y, r, p, q, dimgidx;
		Imgidx imgsize = height * width;

		//merging border(hor)
		for (y = blksz_y;y <= blksz_y * (npartition_y - 1); y += blksz_y)
		{
			q = y * width;
			for (p = q - width;p < q;p++)
			{
				dimgidx = (p << shamt) + neighbor_offset;
				r = qrank[dimgidx];

				connect(parentAry[p], parentAry[p + width], (Imgidx)curSize++, (Pixel)r);

			}
		}
		neighbor_offset = (connectivity == 4) ? 1 : 3;

		//merging border(ver)
		for (x = blksz_x - 1;x <= blksz_x * (npartition_x - 1); x += blksz_x)
		{
			q = x + imgsize;
			for (p = x;p < q;p += width)
			{
				dimgidx = (p << shamt) + neighbor_offset;
				r = qrank[dimgidx];

				connect(parentAry[p], parentAry[p + 1], (Imgidx)curSize++, (Pixel)r);

			}
		}

		//reset rootidxs
		for (p = 0;p < curSize;p++)
			node[p].rootidx = ROOTIDX;

		//Make sure that the root has the highest
		for (p = parentAry[0];node[p].parentidx != ROOTIDX;p = node[p].parentidx)
			;
		if (node[p].alpha != (Pixel)(numbins - 1))
		{
			q = curSize++;
			node[q].copy(node + p);
			node[q].alpha = numbins - 1;
			node[p].parentidx = q;
			node[q].parentidx = node[q].rootidx = ROOTIDX;
		}
	}

	void connect_pilotnode(AlphaNode<Imgidx, Pixel> *pilottree, Imgidx nredges, Imgidx imgsize)
	{
		//curSize = maxSize;
		Imgidx *rootindexcand = (Imgidx*)Calloc(omp_get_max_threads() * sizeof(Imgidx));
		for (int i = 0;i < omp_get_max_threads();i++)
			rootindexcand[i] = ROOTIDX;

		#pragma omp parallel for schedule(guided,1)
		for (Imgidx p = 0; p < maxSize; p++)
		{
			Imgidx q,r,s;
			//printf("p: %d\n",(int)p);
			if (p < imgsize)
			{
				q = parentAry[p];
				if (node[p].parentidx == ROOTIDX && pilottree[q].area == 1)
					node[p].parentidx = pilottree[q].rootidx;
			}
			else
			{
				if (node[p].rootidx == ROOTIDX && node[p].area)
				{
					if (node[p].area == imgsize)
					{
						if (rootindexcand[omp_get_thread_num()] == ROOTIDX)
						{
							rootindexcand[omp_get_thread_num()] = p;
						}
						continue;
					}
					q = node[p].parentidx;
					for (r = q;pilottree[r].rootidx == ROOTIDX;r = pilottree[r].parentidx)
						;

					s = pilottree[r].rootidx;
					while(q != r)
					{
						pilottree[q].rootidx = s;
						q = pilottree[q].parentidx;
					}
					node[p].parentidx = s;
				}
			}
		}

		rootidx = ROOTIDX;
		for (int i = 0;i < omp_get_max_threads();i++)
		{
			if (rootidx == ROOTIDX) rootidx = rootindexcand[i];
			else if (rootindexcand[i] != ROOTIDX)
				rootidx = _min(rootidx, rootindexcand[i]);
		}

		Free(rootindexcand);
	}

	void set_qindex(Imgidx *qindex, Imgidx *dhist, _int64 numpartitions, _int32 numbins, Imgidx npartition_x, Imgidx npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn)
	{

		//add rooms for singleton nodes
		{
			int p = 0;
			for (int x = 0; x < npartition_x;x++)
			{
				Imgidx sx = (x == npartition_x - 1) ? blksz_xn : blksz_x;
				for (int y = 0;y < npartition_y;y++)
				{
					Imgidx sy = (y == npartition_y - 1) ? blksz_yn : blksz_y;
					Imgidx blksize = sx * sy;
					dhist[p] += blksize;
					p += numbins;
				}
			}
		}

		//compute cumulative distribution
		for (int p = 0;p < numpartitions;p++)
		{
			Imgidx *pdhist = dhist + p * numbins;
			for (int q = 0;q < numbins;q++)
			{
				pdhist[q + numbins] += pdhist[q];
			}
		}

		Imgidx *cdhist = dhist + numpartitions * numbins;//histogram of the edges on the subblock borders
		for (int p = 0;p < numbins - 1;p++)
			cdhist[p + 1] += cdhist[p];
		for (int p = 0;p < numbins;p++)
		{
			Imgidx hp = ((p > 0) ? cdhist[p - 1] : 0);
			for (int q = 0;q < numpartitions + 1;q++)
				qindex[p + numbins * q] = hp + ((q > 0) ? dhist[(q - 1) * numbins + p] : 0);
		}
	}

	void set_qindex(Imgidx *qindex, Imgidx *dhist, _int64 numpartitions, _int32 numbins)
	{
		Imgidx imgsize = width * height;
		for (int p = 0;p < numpartitions;p++)
		{
			Imgidx *pdhist = dhist + p * numbins;
			for (int q = 0;q < numbins;q++)
			{
				pdhist[q + numbins] += pdhist[q];
			}
		}
		Imgidx *cdhist = dhist + numpartitions * numbins;
		for (int p = 0;p < numbins - 1;p++)
			cdhist[p + 1] += cdhist[p];
		for (int p = 0;p < numbins;p++)
		{
			Imgidx hp = imgsize + ((p > 0) ? cdhist[p - 1] : 0);
			for (int q = 0;q < numpartitions + 1;q++)
				qindex[p + numbins * q] = hp + ((q > 0) ? dhist[(q - 1) * numbins + p] : 0);
		}
	}

	void set_subtree_root(Imgidx** subtreerootary, Imgidx *strary, Imgidx nonzero_nodeidx_start, Imgidx rootlevel_nodeidx_start)
	{
			Imgidx imgsize = width * height;
			if (node[rootidx].alpha < 2)
			{
				cout << "Pilottree root node level lower than 2" << endl;
				return;
			}

		{

			for (Imgidx p = 1;p < (int)node[rootidx].alpha;p++)
				subtreerootary[p] = strary + (p - 1) * imgsize;

			for (Imgidx p = 0;p <= rootidx;p++)
				node[p].rootidx = ROOTIDX;

			for (Imgidx p = 0;p < imgsize;p++)
			{
				Imgidx q = parentAry[p];
				node[node[q].parentidx].rootidx = 1;
			}

			_int64 areasum = 0;
			for (Imgidx p = nonzero_nodeidx_start;p < rootlevel_nodeidx_start;p++)
			{
				if (node[p].rootidx != ROOTIDX)
				{
					node[node[p].parentidx].rootidx = 0;
					areasum += node[p].area;
					node[p].rootidx = areasum;
				}
			}

			Imgidx *pixelindex = (Imgidx*)Malloc(areasum * sizeof(Imgidx));
			for (Imgidx p = 0;p < imgsize;p++) //1-CCs
			{
				Imgidx q = node[parentAry[p]].parentidx;
				if (q < rootlevel_nodeidx_start)
					pixelindex[--node[q].rootidx] = p;

			}

			for (Imgidx p = nonzero_nodeidx_start;p < rootlevel_nodeidx_start;p++)
			{
				if (node[p].rootidx == ROOTIDX)
					continue;
				{

					Imgidx q = node[p].parentidx;
					if (q >= rootlevel_nodeidx_start)
						continue;
					Imgidx r = node[p].rootidx;
					Imgidx s = r + node[p].area;
					Imgidx t = node[q].rootidx;

					//node[p].rootidx = ROOTIDX;
					while(r < s)
					{
						pixelindex[--t] = pixelindex[r++];

					}

					node[q].rootidx = t;
				}
			}

			for (Imgidx p = 0;p < (int)((node[rootidx].alpha - 1) * imgsize);p++)
				strary[p] = -1;

			for (Imgidx p = 0;p < imgsize;p++)
				strary[p] = parentAry[p];
			for (Imgidx p = nonzero_nodeidx_start;p < rootlevel_nodeidx_start;p++)
			{
				if (!node[p].area)
					continue;

				Imgidx r = node[p].rootidx;
				Imgidx s = node[p].rootidx + node[p].area;
				node[p].rootidx = ROOTIDX;
				Imgidx level = node[p].alpha;
				while(r < s)
				{
					subtreerootary[level][pixelindex[r++]] = p;
				}
			}

			for (Imgidx level = 2;level < (int)node[rootidx].alpha;level++)
			{
				Imgidx *pstrary_prev = &subtreerootary[level - 1][0];
				Imgidx *pstrary = &subtreerootary[level][0];
				for (Imgidx p = 0;p < imgsize;p++)
				{
					if (pstrary[p] == -1)
						pstrary[p] = pstrary_prev[p];
				}
			}

			Free(pixelindex);
		}
	}

	//void init_hypergraph_nodes(_uint8* is_redundant, Imgidx *rank, RankItem<Imgidx, Pixel>* rankitem)
	void find_redundant_nodes(_uint8* is_redundant, Imgidx *rank)
	{
		if (connectivity == 4)
		{
			Imgidx imgsize = height * width;
			Imgidx imgsize2 = imgsize * 2;
			Imgidx width2 = width << 1;
			for (Imgidx p = 0; p < imgsize;p++)
			{
				Imgidx q = p << 1;
				Imgidx y = p / width;
				Imgidx x = p % width;
				//_int8 isAv = isAvailable[p];
				Imgidx maxRank = -1;

				((y < height - 1) && (rank[q] > maxRank)) 			? (maxRank = rank[q])  							 : (Imgidx)0;
				((x < width - 1) && (rank[q + 1] > maxRank)) 		? (maxRank = rank[q + 1])  					 : (Imgidx)0;
				((x > 0) && (rank[q - 1] > maxRank)) 						? (maxRank = rank[q - 1])   				 : (Imgidx)0;
				((y > 0) && (rank[q - width2] > maxRank)) 		  ? (maxRank = rank[q - width2]) 			 : (Imgidx)0;

				//is_redundant[rankitem[minRank].dimgidx] = 1;
				is_redundant[maxRank] = 1;
				//node[p].connect_to_parent(&node_in[minRank], minRank + imgsize, pix_type);
			}
		}
		else//later!
		{

		}
	}

	void set_subblock_properties(Imgidx* startpidx, Imgidx *blkws, Imgidx *blkhs, Imgidx *blocksize, _int8 npartition_x, _int8 npartition_y, _int64 blksz_x, _int64 blksz_y, _int64 blksz_xn, _int64 blksz_yn)
	{
		Imgidx p = 0, q;
		for (_int8 y = 0; y < npartition_y; y++)
		{
			q = y * width * (Imgidx)blksz_y;
			bool lastrow = (y == npartition_y - 1);
			Imgidx blkh = lastrow ? blksz_yn : blksz_y;
			for (_int8 x = 0; x < npartition_x; x++)
			{
				startpidx[p] = q + (Imgidx)x * (Imgidx)blksz_x;
				bool lastcol = (x == npartition_x - 1);
				Imgidx blkw = lastcol ? blksz_xn : blksz_x;
				blocksize[p] = blkh * blkw * 2 - (Imgidx)lastrow * blkw - (Imgidx)lastcol * blkh;
				blkws[p] = blkw;
				blkhs[p] = blkh;
				p++;
			}
		}
	}

	void memalloc_queues(HierarQueue<Imgidx>*** queues, _int64 numpartitions, Imgidx* blocksize, Imgidx* subtree_max)
	{
		//preallocate hqueues
		*queues = (HierarQueue<Imgidx>**)Calloc(numpartitions * sizeof(HierarQueue<Imgidx>*));
		for (int blk = 0;blk < numpartitions;blk++)
		{
			//printf("blk%d max = %d\n", blk, (int)subtree_max[blk]);
			(*queues)[blk] = new HierarQueue<Imgidx>((_uint64)blocksize[blk] + 1, (_int32)(subtree_max[blk] + 1));
			//queues[blk] = new HierarQueue<Imgidx>((_uint64)blocksize[blk] + 1, (_int32)(1<<20));
		}
	}

	void compute_dimg_and_rank2index(RankItem<Imgidx, double>*& rankitem, Pixel* img, Imgidx nredges, Index* rank2rankitem)
	{
		if (channel == 1)
		{
			SortValue<Pixel>* vals;// = new pmt::SortValue<Value>[N];
			vals = (SortValue<Pixel>*)Malloc(nredges * sizeof(SortValue<Pixel>));
			compute_dimg_par4(rankitem, img, vals);

			SortPair<Pixel,Index>* sort_space = (SortPair<Pixel,Index>*)Calloc(2 * nredges * sizeof(SortPair<Pixel,Index>));//new SortPair[2 * N];

			rank_to_index((SortValue<Pixel>*)vals, (Index)nredges, rank2rankitem, 0U, (uint_fast8_t)(sizeof(Pixel) << 3), sort_space, omp_get_max_threads());

			Free(vals);
			Free(sort_space);
		}
		else
		{
			SortValue<double>* vals;// = new pmt::SortValue<Value>[N];
			vals = (SortValue<double>*)Malloc(nredges * sizeof(SortValue<double>));
			compute_dimg_par4(rankitem, img, vals);

			SortPair<_uint64,Index>* sort_space = (SortPair<_uint64,Index>*)Calloc(2 * nredges * sizeof(SortPair<_uint64,Index>));//new SortPair[2 * N];

			rank_to_index((SortValue<_uint64>*)vals, (Index)nredges, rank2rankitem, 0U, (uint_fast8_t)(sizeof(double) << 3), sort_space, omp_get_max_threads());

			Free(vals);
			Free(sort_space);
		}
	}

	void compute_difference_and_sort(RankItem<Imgidx, double>*& rankitem, Pixel* img, Imgidx nredges)
	{
		Index* rank2rankitem = (Index*)Calloc(nredges * sizeof(Index));

		compute_dimg_and_rank2index(rankitem, img, nredges, rank2rankitem);

		RankItem<Imgidx, double> *sorteditem = (RankItem<Imgidx, double>*)Malloc(nredges * sizeof(RankItem<Imgidx, double>)), *tmp;
		#pragma omp parallel for schedule(guided,1)
		for (Imgidx i = 0; i < nredges; i++)
		{
			sorteditem[i] = rankitem[rank2rankitem[i]];
		}

		tmp = rankitem;
		rankitem = sorteditem;
		Free(tmp);
		Free(rank2rankitem);
	}

	void print_all_trees(AlphaNode<Imgidx, Pixel>* pilottree)
	{
		printf("Pilotree Start =======================\n");
		AlphaNode<Imgidx, Pixel> *tmp = node;
		node = pilottree;
		print_tree();
		node = tmp;
		printf("Pilotree End =======================\n");
		printf("Refined tree Start =======================\n");
		print_tree();
		printf("Refined tree End =======================\n");
	}

	void compute_difference_and_sort(Imgidx* rank, RankItem<Imgidx, double>*& rankitem, Pixel* img, Imgidx nredges, Index*& rank2rankitem)
	{
		compute_dimg_and_rank2index(rankitem, img, nredges, rank2rankitem);

		#pragma omp parallel for schedule(guided,1)
		for (Imgidx i = 0; i < nredges; i++)
		{
			rank[rankitem[rank2rankitem[i]].dimgidx] = i;
		}
	}

	void Pilot_Rank(Pixel* img, int numthreads)
	{
		Imgidx imgsize, dimgsize, nredges;
		RankItem<Imgidx, double> *rankitem;
		Imgidx p, q;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		rankitem = (RankItem<Imgidx, double>*)Malloc(nredges * sizeof(RankItem<Imgidx, double>));

		Imgidx *rank = (Imgidx*)Calloc((size_t)dimgsize * sizeof(Imgidx));
		_uint8 *qrank = (_uint8*)Malloc((size_t)dimgsize);
		parentAry = (Imgidx*)Calloc((size_t)imgsize * sizeof(Imgidx));


		_int64 numpartitions;
		if (numthreads > 2)
		 	numpartitions = _min(imgsize/2, _min(256, numthreads * 4));
		else
			numpartitions = 2;
		//_int64 numpartitions = numthreads;
		int npartition_x, npartition_y;
		{
			int optpart = 1;
			double optborderlength = (double)numpartitions * (double)imgsize;
			for (int px = 2;px < numpartitions; px++)
			{
				if (numpartitions % px == 0)
				{
					int py = numpartitions / px;

					if (((double)px * (double)height + (double)py * (double)width) < optborderlength)
					{
						optpart = px;
						optborderlength = ((double)px * (double)height + (double)py * (double)width);
					}
				}
			}
			npartition_x = (int)optpart;
			npartition_y = (int)numpartitions / npartition_x;
		}

		_int32 numbins = numpartitions; //number of levels in Quantization (= number of threads used)

		_uint8 *isVisited, *isAvailable;
		isVisited = (_uint8*)Calloc((size_t)((imgsize)));
		isAvailable = (_uint8*)Malloc((size_t)(imgsize));
		set_isAvailable_par(isAvailable, npartition_x, npartition_y);

		Imgidx binsize = nredges / (_int64)numbins;
		numbins = (nredges + binsize - 1) / binsize;
		_int64 numlevels = numbins;//for compatibility
		_int64 blksz_x = width / npartition_x;
		_int64 blksz_y = height / npartition_y;
		_int64 blksz_xn = blksz_x + (width % npartition_x);
		_int64 blksz_yn = blksz_y + (height % npartition_y);
		numpartitions = (_int64)npartition_x * (_int64)npartition_y;

		Imgidx *startpidx = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blocksize = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *subtree_size = (Imgidx*)Calloc((numpartitions) * sizeof(Imgidx));
		Imgidx *subtree_start = (Imgidx*)Calloc((numpartitions + 1) * sizeof(Imgidx));
		Imgidx *subtree_nborderedges = (Imgidx*)Calloc((numpartitions ) * sizeof(Imgidx));
		Imgidx *subtree_cur = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *subtree_max = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blkws = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *blkhs = (Imgidx*)Malloc(numpartitions * sizeof(Imgidx));
		Imgidx *dhist = (Imgidx*)Calloc((size_t)numbins * (size_t)numpartitions * sizeof(Imgidx));
		char *blkflooddone = (char*)Calloc(numpartitions * sizeof(char));
		omp_lock_t *locks = (omp_lock_t*)Malloc((numpartitions + 1) * sizeof(omp_lock_t));

		for (p = 0;p < numpartitions;p++)
			omp_init_lock(locks + p);

		double *nrmsds = (double*)Malloc(numpartitions * sizeof(double));

		omp_set_num_threads(_min(numthreads, omp_get_num_procs()));

		Index* rank2rankitem = (Index*)Calloc(nredges * sizeof(Index));

		compute_difference_and_sort(rank, rankitem, img, nredges, rank2rankitem);

		set_subblock_properties(startpidx, blkws, blkhs, blocksize, npartition_x, npartition_y, blksz_x, blksz_y, blksz_xn, blksz_yn);
		quantize_ranks_compute_histogram(qrank, rank, img, dhist, blkws, blkhs, startpidx, binsize, numbins, npartition_x, npartition_y, subtree_max);

		blockwise_tse(subtree_size, subtree_nborderedges, nrmsds, dhist, subtree_max, blkws, blkhs, npartition_x, npartition_y, numbins);

		HierarQueue<Imgidx>** queues;
		memalloc_queues(&queues, numpartitions, blocksize, subtree_max);

		Imgidx *hypernode_level = (Imgidx*)Malloc(dimgsize * sizeof(Imgidx));
		_int8 *levelroots = (_int8*)Calloc(numbins * numpartitions * sizeof(_int8));
		_int8 *redundant_edge = (_int8*)Calloc(dimgsize * sizeof(_int8));

		double treesizemult_intv = 0.05;
		double treesizemult = 1.0 - treesizemult_intv;

		int flooddone = 0;
		while(!flooddone)
		{
			int numbusythr = 0;
			int outofmemory = 0;
			int numblkproc = 0;

			Imgidx shamt = connectivity >> 2;
			Imgidx wstride_d = width << shamt;

			//reset queue, isvisited array, hypernode levels
			for (int blk = 0;blk < numpartitions;blk++)
				queues[blk]->reset_queue();

			#pragma omp parallel for private(p, q)
			for (int i = 0;i < imgsize;i++)
				isVisited[i] = 0;

			#pragma omp parallel for private(p, q)
			for (Imgidx i = 0; i < dimgsize;i++)
			{
				hypernode_level[i] = ROOTIDX;
				redundant_edge[i] = 0;
			}

			flooddone = 1; //reset this flag when memory overflows on all threads
			treesizemult = treesizemult + treesizemult_intv;

			//(re)allocate node array (expand size when reallocate)
			maxSize = parflood_node_alloc(subtree_size, subtree_start, blkws, blkhs, numpartitions, treesizemult);

			for (p = 0;p < numpartitions;p++)
				omp_unset_lock(locks + p);
			omp_unset_lock(locks + numpartitions);
			numbusythr = 0;


			#pragma omp parallel for private(p, q) schedule(dynamic,1)
			for (int blk = 0; blk < numpartitions; blk++) //flooding is somehow slower than the one without tse...
			{
				if (outofmemory)
					continue;

				omp_set_lock(locks + blk);
				omp_set_lock(locks + numpartitions);
				numbusythr++;
				numblkproc++;
				omp_unset_lock(locks + numpartitions);

				_uint8 connected_neighbor; //marks neighbors that are already visited
				Imgidx lsbclearmask = ~1; //mask for clearing 1st bit
				Imgidx bwidth = blkws[blk];
				Imgidx bheight = blkhs[blk];
				//Imgidx blksize = blocksize[blk];
				Imgidx bareasum = bwidth * bheight;
				Imgidx *bhist = dhist + numlevels * blk;
				HierarQueue<Imgidx> *queue = queues[blk];
				//bool lastcol = (blk % npartition_x) == (npartition_x - 1);
				//bool lastrow = (blk / npartition_x) == (npartition_y - 1);
				//Imgidx spidx = startpidx[blk];
				Imgidx nidx = subtree_start[blk];
				Imgidx blkts = 0;
				int nidxblk = blk;
				Imgidx nidx_lim = subtree_start[blk + 1] - subtree_nborderedges[blk + 1];//save room for nodes to be added in merge
				Imgidx iNode = 0;
				Pixel maxdiff = subtree_max[blk];

				_int8 *plr = levelroots + blk * numbins;
				for (p = 0;p < numbins;p++)
					plr[p] = 0;

				//printf("th%d: subtree for blk %d: setting queue\n", omp_get_thread_num(), (int)blk);
				queue->set_queue(bhist);
				//queue->set_queue(bhist, maxdiff);

				Imgidx stack_top = nidx++;//imgsize + dimgsize + blk;
				//printf("blk%d - Dummy %d\n", (int)blk, (int)(stack_top));
				Imgidx prev_top = stack_top;
				AlphaNode<Imgidx, Pixel> *pNode = node + stack_top;
				pNode->set(0, maxdiff, (double)0.0, (Pixel)-1, (Pixel)0, pix_type);
				pNode->parentidx = stack_top;
				pNode->rootidx = ROOTIDX;
				Pixel current_level = maxdiff;

				Imgidx x0 = startpidx[blk]; /*starting point*/
				queue->push(((x0 << shamt) << 1) & lsbclearmask, current_level);
				prev_top = stack_top; /*to find redundant node*/
				int firstpix = 1;

				if (outofmemory)
					continue;
				while (1) //flooding
				{
					while ((_int64)queue->min_level <= (_int64)current_level) //flood all levels below current_level
					{
						Imgidx qitem = queue->pop();
						Imgidx didx = qitem >> 1;
						//the pixel which pushed this item into the queue, where it is easier to find the level of the edge (level = stack_top)

						p = didx >> shamt;

						_uint8 isAv;
						if (firstpix)
						{
							isAv = isAvailable[p];
							firstpix = 0;
						}
						else
						{
							if ((qitem & 1)) //pixel at another end of the edge?
							{
								if (connectivity == 4)
								{
									if (didx & 1)//horizontal edge
									{
											p++;
											isAv = isAvailable[p];
											isAv &= ~0x4; //do not check the neighbor which pushed this pixel (p) into the queue
									}
									else//vertical edge
									{
											p += width;
											isAv = isAvailable[p];
											isAv &= ~0x8;
									}
								}
								else
								{
									isAv = 0;
								}
							}
							else
							{
								if (didx & 1)//horizontal edge
								{
										isAv = isAvailable[p];
										isAv &= ~0x2;
								}
								else//vertical edge
								{
										isAv = isAvailable[p];
										isAv &= ~0x1;
								}
							}
						}

						//printf("thr%d: probing %d\n", omp_get_thread_num(), (int)p);
						if (isVisited[p])
						{
							queue->find_minlev();
							continue;
						}

						isVisited[p] = 1;
						connected_neighbor = 0;
						if (connectivity == 4)
						{
							q = p << shamt;
							Imgidx q1;
							if (is_available(isAv, 0))
							{
								if (isVisited[p + width]) //neighbor alread visited - which means this edge might be on lower level than the alpha value of the edge corresponds to
								{
									//if (!plr[qrank[q]]) redundant_edge[q] = 1;
									//else
									connected_neighbor |= 0x1;
								}
								else//push new neighbor
									queue->push((q << 1) | 1, qrank[q]);
							}
							if (is_available(isAv, 1))
							{
								q1 = q + 1;
								if (isVisited[p + 1])
								{
									//if (!plr[qrank[q1]]) redundant_edge[q1] = 1;
									//else
									connected_neighbor |= 0x2;
								}
								else
									queue->push(((q1) << 1) | 1, qrank[q1]);
							}
							if (is_available(isAv, 2))
							{
								q1 = q - 1;
								if (isVisited[p - 1])
								{
									//if (!plr[qrank[q1]]) redundant_edge[q1] = 1;
									//else
									connected_neighbor |= 0x4;
								}
								else
									queue->push(((q1) << 1) & lsbclearmask, qrank[q1]);
							}

							if (is_available(isAv, 3))
							{
								q1 = q - wstride_d;
								if (isVisited[p - width])
								{
									//if (!plr[qrank[q1]]) redundant_edge[q1] = 1;
									//else
									connected_neighbor |= 0x8;
								}
								else
									queue->push(((q1) << 1) & lsbclearmask, qrank[q1]);
							}
						}

						if ((_int64)current_level > (_int64)queue->min_level) //go to lower level
						{
							//plr[queue->min_level] = 1;
							Pixel pix_val = img[p];
							current_level = queue->min_level;

							{
								if (nidx == nidx_lim)
								{
									if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
										break;
								}
								iNode = nidx++;
							}
							node[iNode].set(1, current_level, (double)pix_val, pix_val, pix_val, pix_type);
							node[iNode].parentidx = stack_top;
							node[iNode].rootidx = ROOTIDX;
							stack_top = iNode;

							if (current_level)
							{
								{
									if (nidx == nidx_lim)
									{
										if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
											break;
									}
									iNode = nidx++;
								}
								node[iNode].copy(node + stack_top);
								node[iNode].alpha = 0;
								node[iNode].parentidx = stack_top;
								node[iNode].rootidx = ROOTIDX;
								prev_top = iNode;
							}
							parentAry[p] = iNode;
						}
						else
						{
							queue->find_minlev();

							if (current_level)
							{
								Pixel pix_val = img[p];
								{
									if (nidx == nidx_lim)
									{
										if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
											break;
									}
									iNode = nidx++;
								}
								node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val, pix_type);
								node[stack_top].add(node + iNode, pix_type);
								node[iNode].parentidx = stack_top;
								node[iNode].rootidx = ROOTIDX;
								parentAry[p] = iNode;
							}
							else
							{
								parentAry[p] = stack_top;
								node[stack_top].add(img[p], pix_type);
							}
						}
						//if (stack_top == ROOTIDX) printf("SQEEEEAK\n");

						if (connected_neighbor)
						{
							// mark from the leaf to the stack_top node to help finding hypernode levelroots
							Imgidx squirrel, leaf1 = parentAry[p], leaf2;
							for (squirrel = leaf1;node[squirrel].parentidx != squirrel;squirrel = node[squirrel].parentidx)
								node[squirrel].rootidx = p;
							node[squirrel].rootidx = p;

							q = p << shamt;
							if (connected_neighbor & 0x1)
							{
								leaf2 = parentAry[p + width];
								if (qrank[q] > get_nearest_common_ancestor_level(p, leaf2))
									redundant_edge[q] = 1;
							}
							if (connected_neighbor & 0x2)
							{
								leaf2 = parentAry[p + 1];
								//printf("q = %d\n", (int)q);
								if (qrank[q + 1] > get_nearest_common_ancestor_level(p, leaf2))
									redundant_edge[q + 1] = 1;
							}
							if (connected_neighbor & 0x4)
							{
								leaf2 = parentAry[p - 1];
								if (qrank[q - 1] > get_nearest_common_ancestor_level(p, leaf2))
									redundant_edge[q - 1] = 1;
							}
							if (connected_neighbor & 0x8)
							{
								leaf2 = parentAry[p - width];
								if (qrank[q - wstride_d] > get_nearest_common_ancestor_level(p, leaf2))
									redundant_edge[q - wstride_d] = 1;
							}

							//cleanup pawprints
							for (squirrel = leaf1;node[squirrel].parentidx != squirrel;squirrel = node[squirrel].parentidx)
								node[squirrel].rootidx = ROOTIDX;
							node[squirrel].rootidx = ROOTIDX;
						}
					}

					if (outofmemory)
						break;

					//remove_redundant_node(node, nidx, prev_top, stack_top);
					if (node[prev_top].parentidx == stack_top && node[prev_top].area == node[stack_top].area)
					{
						//plr[(int)(node[prev_top].alpha)] = 0;
						node[prev_top].parentidx = node[stack_top].parentidx;
						stack_top = prev_top;
						//curSize--;
					}

					if (node[stack_top].area == bareasum)	// root node found...done
						break;

					//go to higher level
					iNode = node[stack_top].parentidx;
					if ((_int64)queue->min_level < (_int64)node[iNode].alpha) //new level from queue
					{
						//plr[queue->min_level] = 1;
						//_CREATE_NEW_STACKTOP(queue->min_level)
						{
							if (nidx == nidx_lim)
							{
								if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,	blkflooddone, subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr, numblkproc, outofmemory))
									break;
							}
							iNode = nidx++;
						}
						node[iNode].alpha = queue->min_level;
						node[iNode].copy(node + stack_top);
						node[iNode].parentidx = node[stack_top].parentidx;
						node[iNode].rootidx = ROOTIDX;
						node[stack_top].parentidx = iNode;
					}
					else //go to existing node
					{
						if (node[iNode].area == bareasum)	// root node found...done
							break;
						node[iNode].add(node + stack_top, pix_type);
					}


					if (node[iNode].area == bareasum)	// root node found...done
						break;

					prev_top = stack_top;
					stack_top = iNode;
					current_level = node[stack_top].alpha;
				}

				if (!outofmemory)
				{
					stack_top = (node[stack_top].area == bareasum) ? stack_top : iNode; //remove redundant root
					node[stack_top].parentidx = ROOTIDX;

					subtree_cur[nidxblk] = nidx;

					blkts += nidx - subtree_start[nidxblk];

					if (nidx == nidx_lim)//this should be really rare
					{
						blkflooddone[nidxblk] = 2; // flood done (at least for the native block), no free memory
					}
					else
					{
						blkflooddone[nidxblk] = 1; // flood done AND free memory available
					}

					omp_set_lock(locks + numpartitions);
					numbusythr--;
					//printf("th%d-- (%d/%d)\n", omp_get_thread_num(), numbusythr, omp_get_num_threads());
					omp_unset_lock(locks + numpartitions);
					omp_unset_lock(locks + nidxblk);
				}

				if (outofmemory)
				{
					printf("thr%d: subtree for blk %d: memory overflow - releasing lock %d\n", omp_get_thread_num(), (int)blk, (int)nidxblk);
					flooddone = 0;
				}
				else
				{
					//printf("thr%d: subtree for blk %d: releasing lock %d\n", omp_get_thread_num(), (int)blk, (int)nidxblk);
				}
				//omp_unset_lock(locks + nidxblk);
			}//flood_end
		}

		if (numpartitions > 1)
			rootidx = merge_subtrees1(qrank, blksz_x, blksz_y, npartition_x, npartition_y, subtree_cur, 1, hypernode_level);
		else
		{
			rootidx = parentAry[0];
			while(node[rootidx].parentidx != ROOTIDX)
				rootidx = node[rootidx].parentidx;
		}

		////////////////////////////////////////////////////////////////////////
		// Initialize refined tree
		////////////////////////////////////////////////////////////////////////

		#pragma omp parallel for
		for (p = 0;p < maxSize;p++)
			node[p].rootidx = ROOTIDX;

		maxSize = imgsize + nredges;
		AlphaNode<Imgidx, Pixel> *pilottree = node;
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)(maxSize + 1) * sizeof(AlphaNode<Imgidx, Pixel>));
		node_in = node + imgsize;

		//print_value_pixel("image",img);
		//print_value_edge("qrank",qrank);
		//print_all_trees(pilottree);
		//getchar();

		//pilot_rank2
		Pixel maxpixval = (Pixel)(~0) >> (((sizeof(Pixel) << 3)) - bit_depth);
		initialize_node_par1(img, rankitem, maxpixval, rank2rankitem);

		double *levelperthread = (double*)Calloc((numbins + 1) * sizeof(double));
		double *timeperthread = (double*)Calloc((numbins + 1) * sizeof(double));

		#pragma omp parallel for schedule(dynamic,1)
		for (_int64 qlevel = 0; qlevel < numbins; qlevel++) //this part is to be parallelised
		{
			if (qlevel > (_int64)pilottree[rootidx].alpha)
			{
				continue;
			}
			unionfind_refine_qlevel(qlevel, binsize, nredges, pilottree, rankitem, redundant_edge, rank2rankitem);

			levelperthread[omp_get_thread_num()]++;
		}

		connect_pilotnode(pilottree, nredges, imgsize);
		Free(parentAry);
		parentAry = 0;
		node[rootidx].parentidx = ROOTIDX;

		Free(levelperthread);
		Free(timeperthread);
		if (rank2rankitem)		Free(rank2rankitem);

		for (p = 0;p < numpartitions;p++)
				omp_destroy_lock(locks + p);
		Free(locks);
		for (int blk = 0;blk < numpartitions;blk++)
			delete queues[blk];
		Free(queues);
		Free(nrmsds);
		Free(blkflooddone);
		Free(blkws);
		Free(blkhs);
		Free(dhist);
		Free(startpidx);
		Free(blocksize);
		Free(subtree_start);
		Free(subtree_nborderedges);
		Free(subtree_size);
		Free(subtree_cur);
		Free(subtree_max);
		Free(hypernode_level);
		Free(levelroots);
		Free(redundant_edge);
		Free(pilottree);
		Free(rank);
		Free(rankitem);
		Free(isVisited);
		Free(isAvailable);
		Free(qrank);

	}

	inline Imgidx NewAlphaNode(Imgidx &size, Imgidx &maxsize)
	{
		if (size == maxsize)
		{
//std::cout << "Reallocating...\n";
		}
		return size++;
	}

	inline Imgidx NewAlphaNode(AlphaNode<Imgidx, Pixel> *tree, Imgidx &size, Imgidx &maxsize, Pixel level, AlphaNode<Imgidx, Pixel> *pCopy)
	{
		AlphaNode<Imgidx, Pixel> *pNew = tree + size;

			if (size == maxsize)
		{
		//	std::cout << "Reallocating...\n";

		}
		pNew->alpha = level;
		pNew->copy(pCopy);
		return size++;
	}

	inline void remove_redundant_node(AlphaNode<Imgidx, Pixel> *tree, Imgidx &size, Imgidx& prev_top, Imgidx& stack_top)
	{
		if (tree[prev_top].parentidx == stack_top && tree[prev_top].area == tree[stack_top].area)
		{
			tree[prev_top].parentidx = tree[stack_top].parentidx;
			stack_top = prev_top;
			size--;
		}
	}

	//for parallel pilottree
	inline void connectPix2Node(AlphaNode<Imgidx, Pixel> *tree, Imgidx pidx, Pixel pix_val, Imgidx iNode, Imgidx *pAry)
	{
		AlphaNode<Imgidx, Pixel> *pNode = &tree[iNode];
		pAry[pidx] = iNode;
		pNode->add(pix_val, pix_type);
	}

	inline Imgidx find_root1(Imgidx p, Imgidx qlevel)//, int &cnt)
	{
		Imgidx q, r;

		if (p == ROOTIDX)
			return ROOTIDX;

	//		if (!path_compression)
	//		{
	//		while((q = node[p].parentidx) != ROOTIDX)
	//			p = q;
	//			return p;
	//		}
	//		else
		{
			r = p;
			while((q = node[p].rootidx) != ROOTIDX)
				p = q;

			while(r != p)
			{
				q = node[r].rootidx;
				node[r].rootidx = p;
				//checkcoherence(node + r, qlevel);
				//cnt++;
				r = q;
			}

			return p;
		}
	}

	void create_queues(HierarQueue<Imgidx> ***queues, Imgidx nredges, Imgidx binsize, _int32 numbins)
	{
		Imgidx *qh = (Imgidx*)Calloc((2 * (numbins - 2) + 3) * sizeof(Imgidx));
		*queues = new HierarQueue<Imgidx>*[numbins];
		HierarQueue<Imgidx> **Q = *queues;
		qh[numbins - 2] = binsize;
		qh[numbins - 1] = nredges - binsize;
		Q[0] = new HierarQueue<Imgidx>(nredges, qh + numbins - 2, numbins);
		//Q[0] = new HierarQueue<Imgidx>(numbins, nredges);
		for (_int32 q = 1;q < numbins;q++)
		{
			qh[numbins - 2] = binsize * q;
			qh[numbins - 1] = (q == numbins - 1) ? nredges - qh[0] : binsize;
			qh[numbins] = nredges - qh[numbins - 1] - qh[numbins - 2];
			Q[q] = new HierarQueue<Imgidx>(nredges, qh + numbins - 1 - q, numbins);
		}
		Free(qh);
	}

	template <class Value>
	void print_value_pixel(const char *name, Value* arr)
	{
		Imgidx x, y;
		Imgidx p = 0;
		printf("%s: \n",name);
		for (y = 0;y < height;y++)
		{
			for (x = 0;x < width;x++)
			{
				printf("%2d ",(int)arr[p++]);
			}
			printf("\n");
		}
	}

	template <class Value>
	void print_value_edge(const char *name, Value* arr, Imgidx indexing = 0)
	{
		Imgidx x, y;
		Imgidx p;
		printf("%s: \n",name);
		if (indexing)
		{
			for (y = 0;y < height;y++)
			{
				Imgidx p = y * 2 * width - y + (y < height - 1);
				for (x = 0;x < width - 2;x++)
				{
						printf(". %d ", (int)arr[p]);
						p += 1 + (y < height - 1);
		//						cout << rank[p++] << '/' << rank[p++] << ' ';
				}
				printf(". %d .\n", (int)arr[p]);

				if (y < height - 1)
				{
					p = y * 2 * width - y;
					for (x = 0;x < width;x++)
					{
						printf("%d   ", (int)arr[p]);
						p += 2;
						//						cout << rank[p++] << '/' << rank[p++] << ' ';
					}
					printf("\n");
				}
			}
		}
		else
		{
			for (y = 0;y < height;y++)
			{
				p = y * 2 * width;
				for (x = 0;x < width - 2;x++)
				{
						printf(". %-2d ", (int)arr[p+1]);
						p += 2;
		//						cout << rank[p++] << '/' << rank[p++] << ' ';
				}
				printf(". %-2d .\n", (int)arr[p+1]);


				if (y < height - 1)
				{
					p = y * 2 * width;
					for (x = 0;x < width;x++)
					{
						printf("%-2d   ", (int)arr[p]);
						p += 2;
						//						cout << rank[p++] << '/' << rank[p++] << ' ';
					}
					printf("\n");
				}
			}
		}
	}

	template<class Queue>
	void Flood_Trie(Pixel* img, Queue *queue)
	{
		Imgidx imgsize, dimgsize, nredges;
		Imgidx current_rank = 0, next_rank = 0;
		RankItem<Imgidx, double> *rankitem, *pRank;
		//AlphaNode<Imgidx, Pixel> *pNode;
		Pixel maxpixval;
		Imgidx *rank, top_rank;
		_int8 nbits;
		//Imgidx *dhist;
		Imgidx prev_top = 0;
		_uint8 *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;
		Imgidx p, q;
		imgsize = width * height;
		nredges = width * (height - 1) + (width - 1) * height + ((connectivity == 8) ? ((width - 1) * (height - 1) * 2) : 0);
		dimgsize = (connectivity >> 1) * width * height;
		maxSize = imgsize + nredges;
		num_node = maxSize;
		num_node_in = nredges;
		nbits = ((sizeof(Pixel) << 3) - 1);
		maxpixval = ~(1 << nbits);
		rankitem = (RankItem<Imgidx, double>*)Malloc(nredges * sizeof(RankItem<Imgidx, double>));
		parentAry = 0;
		rank = (Imgidx*)Malloc((size_t)dimgsize * sizeof(Imgidx));
		node = (AlphaNode<Imgidx, Pixel>*)Malloc((size_t)maxSize * sizeof(AlphaNode<Imgidx, Pixel>));
		node_in = node + imgsize;
		isVisited = (_uint8*)Calloc((size_t)((imgsize)));
		isAvailable = (_uint8*)Malloc((size_t)(imgsize));

		set_isAvailable(isAvailable);

		queue = new Queue(nredges);

		omp_set_num_threads(1);
		Index* rank2rankitem = (Index*)Calloc(nredges * sizeof(Index));
		compute_difference_and_sort(rank, rankitem, img, nredges, rank2rankitem);
		//compute_dimg(rank, rankitem, img);


		initialize_node1(img, rankitem, maxpixval, rank2rankitem);

		//manually visit the first pixel
		isVisited[0] = 1;
		if (connectivity == 4)
		{
			queue->push(rank[0]);
			queue->push(rank[1]);
		}
		else if (connectivity == 8)
		{
			queue->push(rank[0]);
			queue->push(rank[1]);
			queue->push(rank[2]);
		}
		//else later
		current_rank = queue->top();
		node[0].connect_to_parent(&node_in[current_rank], current_rank + imgsize, pix_type);
		prev_top = current_rank;

		while (1)
		{
			while (1)
			{
				top_rank = queue->top();
				pRank = rankitem + rank2rankitem[top_rank];
				if (isVisited[pRank->get_pidx0(connectivity)])
				{
					if (isVisited[pRank->get_pidx1(width,connectivity)])
						break;
					p = pRank->get_pidx1(width,connectivity);
				}
				else
					p = pRank->get_pidx0(connectivity);

				isVisited[p] = 1;
				isAv = isAvailable[p];
				if (connectivity == 4)
				{
					q = p << 1;
					if (is_available(isAv, 0) && !isVisited[p + width])	queue->push(rank[q]);
					if (is_available(isAv, 1) && !isVisited[p + 1])			queue->push(rank[q + 1]);
					if (is_available(isAv, 2) && !isVisited[p - 1])			queue->push(rank[q - 1]);
					if (is_available(isAv, 3) && !isVisited[p - width])	queue->push(rank[q - (width << 1)]);
				}
				else if (connectivity == 8)
				{
					Imgidx width4 = width << 2;
					q = p << 2;
					if (is_available(isAv, 0) && !isVisited[p + width]) 		 queue->push(rank[q]);// printf("0:pushing %d \n",(int)rank[q]);}
					if (is_available(isAv, 1) && !isVisited[p + width + 1]) queue->push(rank[q + 1]);// printf("1:pushing %d \n",(int)rank[q+1]);}
					if (is_available(isAv, 2) && !isVisited[p + 1]) 		 		 queue->push(rank[q + 2]);// printf("2:pushing %d \n",(int)rank[q+2]);}
					if (is_available(isAv, 3) && !isVisited[p - width + 1]) queue->push(rank[q + 3]);// printf("3:pushing %d \n",(int)rank[q+3]);}
					if (is_available(isAv, 4) && !isVisited[p - width]) 		 queue->push(rank[q - width4]);// printf("4:pushing %d \n",(int)rank[q-width4]);}
					if (is_available(isAv, 5) && !isVisited[p - width - 1]) queue->push(rank[q - width4 - 3]);// printf("5:pushing %d \n",(int)rank[q-width4-3]);}
					if (is_available(isAv, 6) && !isVisited[p - 1]) 				 queue->push(rank[q - 2]);//  printf("6:pushing %d \n",(int)rank[q-2]);}
					if (is_available(isAv, 7) && !isVisited[p + width - 1]) queue->push(rank[q + width4 - 1]);// printf("7:pushing %d \n",(int)rank[q+width4-1]);}
				}
				else
				{
					//?
				}
				//else later

				next_rank = queue->top();
				node[p].connect_to_parent(&node_in[next_rank], next_rank + imgsize, pix_type);
				if (current_rank == next_rank)
					break;
				current_rank = next_rank;
			}

			queue->pop();
			next_rank = queue->top();

			//remove redundant node
			if (node_in[prev_top].parentidx == current_rank + imgsize && node_in[prev_top].area == node_in[current_rank].area)
				current_rank = prev_top;

			node_in[current_rank].connect_to_parent(&node_in[next_rank], next_rank + imgsize, pix_type);
			if (node_in[next_rank].area == imgsize)
				break;

			prev_top = current_rank;
			current_rank = next_rank;
		}

		rootidx = (node_in[current_rank].area == imgsize) ? current_rank + imgsize : next_rank + imgsize;
		node[rootidx].parentidx = ROOTIDX;

		delete queue;
		Free(rank2rankitem);
		Free(rank);
		Free(rankitem);
		Free(isVisited);
		Free(isAvailable);
	}

	//the name is misleading. It doesn't find levelroot but the highest node below nodeidx.
	inline Imgidx get_level_root(Imgidx p, Imgidx nodeidx)
	{
		while(node[p].parentidx != ROOTIDX && nodeidx > node[p].parentidx)
			p = node[p].parentidx;

		return p;
	}

	inline Imgidx get_level_root(Imgidx p)
	{
		return get_level_root(p, node);
	}

	inline Imgidx get_level_root(Imgidx p, Pixel alpha)
	{
		while(node[p].parentidx != ROOTIDX && alpha >= node[node[p].parentidx].alpha)
			p = node[p].parentidx;

		return p;
	}

	inline Imgidx get_level_root(Imgidx p, AlphaNode<Imgidx, Pixel> *tree)
	{
		if (p == ROOTIDX)
			return ROOTIDX;
		Pixel a = tree[p].alpha;
		while(1)
		{
			Imgidx parent = tree[p].parentidx;
			if (parent == ROOTIDX || tree[parent].alpha > a)
				break;
			p = parent;
			//cntcnt++;
		}

		return p;
	}

	void swap(Imgidx& x, Imgidx& y)
	{
		Imgidx tmp = x;
		x = y;
		y = tmp;
	}

	void swap(AlphaNode<Imgidx, Pixel>** x, AlphaNode<Imgidx, Pixel>** y)
	{
		AlphaNode<Imgidx, Pixel>* tmp = *x;
		*x = *y;
		*y = tmp;
	}

	inline Imgidx get_nearest_common_ancestor(Imgidx x, Imgidx y)
	{
		Imgidx p;
		for (p = y;p != ROOTIDX;p = node[p].parentidx)
			node[p].rootidx = ROOTIDX;
		for (p = x;p != ROOTIDX;p = node[p].parentidx)
			node[p].rootidx = x;
		for (p = y;p != ROOTIDX && node[p].rootidx != x;p = node[p].parentidx)
			;
		return p;
	}

	inline Pixel get_nearest_common_ancestor_level(Imgidx x, Imgidx y)
	{
		Imgidx p;
		for (p = y;node[p].rootidx != x;p = node[p].parentidx)
			;
		return node[p].alpha;
	}

	Pixel connect(Imgidx x, Imgidx y, Imgidx newidx, Pixel alpha)
	{
		Imgidx x0, y0, z;
	//		bool compxy;
		Imgidx imgsize = height * width;
	//		Imgidx x1 = x, y1 = x;
		AlphaNode<Imgidx, Pixel> n0, n1, *p, *q;
		p = &n0;
		q = &n1;

		x = get_level_root(x, alpha);
		y = get_level_root(y, alpha);

		if (x == y)
		{
			return node[x].alpha;
		}

		z = get_level_root(get_nearest_common_ancestor(x, y));

		node[newidx].copy(node + x);
		node[newidx].add(node + y, pix_type);
		node[newidx].alpha = alpha;
		node[newidx].parentidx = ROOTIDX;

		x0 = x;
		y0 = y;
		x = get_level_root(node[x].parentidx);
		y = get_level_root(node[y].parentidx);

		node[x0].parentidx = newidx;
		node[y0].parentidx = newidx;

		if (x == y)
		{
			if (x == ROOTIDX && y == ROOTIDX)
			{
				x = newidx;
			}
			else
			{
				node[newidx].parentidx = x;
			}
			return node[x].alpha;
		}

		//compxy = parentAry ? node[x].alpha > node[y].alpha : x > y;

		//y always has bigger alpha, or the same alpha but bigger area (for shorter path to level roots)
		if (x == ROOTIDX ||
			(y != ROOTIDX && ((node[x].alpha > node[y].alpha) ||
			(node[x].alpha == node[y].alpha && node[x].area > node[y].area))))
		{
			q->copy(node + x0);
			swap(x,y);
		}
		else
			q->copy(node + y0);

		node[newidx].parentidx = x;
		p->copy(node + x);
		node[x].add(q, pix_type);

		while(1)
		{
			if (y == ROOTIDX)
			{
				while(node[x].parentidx != ROOTIDX)
				{
					x = (node[x].parentidx);
					node[x].add(q, pix_type);
					if (node[x].area == imgsize)
					{
						node[x].parentidx = ROOTIDX;
						break;
					}
					//node[x].print(node,pix_type);
				}
				break;
			}

			if (y == z)//y is a common ancestor
			{

				if (x != y)
				{
					while(1)
					{
						x = (node[x].parentidx);
						if (x == y)
							break;
						node[x].add(q, pix_type);
						if (node[x].area == imgsize)
						{
							node[x].parentidx = ROOTIDX;
							break;
						}
					}
				}
				break;
			}

			while(1)
			{
				x0 = get_level_root(node[x].parentidx);
				if (x0 != ROOTIDX && (node[x0].alpha < node[y].alpha))
				{
					x = x0;
					x0 = get_level_root(node[x0].parentidx);
					p->copy(node + x);
					node[x].add(q, pix_type);
					if (node[x].area == imgsize)
					{
						node[x].parentidx = ROOTIDX;
						break;
					}
				}
				else
					break;
			}

			x0 = get_level_root(node[x].parentidx);
			node[x].parentidx = y;
			q->copy(node + y);
			node[y].add(p, pix_type);

			if (node[y].area == imgsize)
			{
				node[y].parentidx = ROOTIDX;
				break;
			}

			x = x0;
			swap(x,y);
			swap(&p,&q);
		}

		return alpha;
	}

	Pixel connect(Imgidx x, Imgidx y, Pixel alpha, Imgidx newidx)
	{
		Imgidx x0, y0, z;
//		bool compxy;
		Imgidx imgsize = height * width;
//		Imgidx x1 = x, y1 = x;
		AlphaNode<Imgidx, Pixel> n0, n1, *p, *q;
		p = &n0;
		q = &n1;

		x = get_level_root(x, newidx);
		y = get_level_root(y, newidx);
		z = get_nearest_common_ancestor(x, y);

		if (x == y)
		{
			return node[x].alpha;
		}

		node[newidx].copy(node + x);
		node[newidx].add(node + y, pix_type);
		node[newidx].alpha = alpha;
		node[newidx].parentidx = ROOTIDX;

		x0 = x;
		y0 = y;
		x = node[x].parentidx;
		y = node[y].parentidx;

		node[x0].parentidx = newidx;
		node[y0].parentidx = newidx;

		if (x == y)
		{
			if (x == ROOTIDX && y == ROOTIDX)
			{

			}
			else
			{
				node[newidx].parentidx = x;
			}
			return node[x].alpha;
		}

		//compxy = parentAry ? node[x].alpha > node[y].alpha : x > y;
		if (x == ROOTIDX || (y != ROOTIDX && (x > y)))
		{
			q->copy(node + x0);
			swap(x,y);
		}
		else
			q->copy(node + y0);

		node[newidx].parentidx = x;
		p->copy(node + x);
		node[x].add(q, pix_type);

		while(1)
		{
			if (y == ROOTIDX)
			{
				while(node[x].parentidx != ROOTIDX)
				{
					x = node[x].parentidx;
					node[x].add(q, pix_type);
					if (node[x].area == imgsize)
					{
						node[x].parentidx = ROOTIDX;
						break;
					}
				}
				break;
			}

			if (y == z)//y is a common ancestor
			{
				if (x != y)
				{
					while(1)
					{
						x = node[x].parentidx;
						if (x == y)
							break;
						node[x].add(q, pix_type);
						if (node[x].area == imgsize)
						{
							node[x].parentidx = ROOTIDX;
							break;
						}
					}
				}
				break;
			}

			while(1)
			{
				x0 = node[x].parentidx;
				if (x0 != ROOTIDX && (x0 < y))
				{
					x = x0;
					x0 = node[x0].parentidx;
					p->copy(node + x);
					node[x].add(q, pix_type);
					if (node[x].area == imgsize)
					{
						node[x].parentidx = ROOTIDX;
						break;
					}
				}
				else
					break;
			}

			x0 = node[x].parentidx;
			node[x].parentidx = y;
			q->copy(node + y);
			node[y].add(p, pix_type);

			if (node[y].area == imgsize)
			{
				node[y].parentidx = ROOTIDX;
				break;
			}

			x = x0;
			swap(x,y);
			swap(&p,&q);
		}

		return alpha;
	}

	void canonicalize()
	{
		Imgidx p;
		Imgidx numcan = 0;
		Imgidx imgsize = height * width;

		for (p = 0;p < imgsize;p++)
		{
			Imgidx q = parentAry[p];
			Imgidx r = q;

			//Canonicalize leaf nodes
			if (r != ROOTIDX && node[q].alpha == node[node[q].parentidx].alpha)
			{
				while(node[r].alpha == node[node[r].parentidx].alpha)
				{
					node[r].area = 0;
					r = node[r].parentidx;
				}
				parentAry[p] = r;
				numcan++;
			}
		}

		for (p = maxSize - 1;p >= 0; p--)
		{
			if (node[p].area && node[p].parentidx != ROOTIDX && node[node[p].parentidx].parentidx != ROOTIDX &&
				node[node[p].parentidx].alpha == node[node[node[p].parentidx].parentidx].alpha)
			{
				Imgidx q = node[p].parentidx;
				Imgidx r = node[q].parentidx;
				node[p].parentidx = r;
				numcan++;
			}
		}

		for (p = maxSize;p < maxSize;p++)
		{
			Imgidx q = p;
			if (p < imgsize)
				q = parentAry[p];

			if (node[q].area == 0 && (node[p].parentidx != ROOTIDX || node[p].rootidx != ROOTIDX))
			{
				node[q].parentidx = node[q].rootidx = ROOTIDX;
			}

		}
	}

	void merge_subtrees(Imgidx *rank, RankItem<Imgidx, Pixel>* rankitem, _int64 blksz_x, _int64 blksz_y, Imgidx neighbor_offset, Imgidx shamt, Imgidx npartition_x, Imgidx npartition_y)
	{
		Imgidx imgsize = height * width, numblk;
		while(npartition_x > 1 || npartition_y > 1)
		{
			if ((npartition_x == 1 || blksz_x >= blksz_y) && npartition_y > 1)
			{
				numblk = npartition_x * (npartition_y / 2);
				#pragma omp parallel for
				for (int blk = 0;blk < numblk;blk++)
				{
					Imgidx x, y, r, p, p0, pn, dimgidx;
					y = (1 + 2 * (blk / (int)npartition_x)) * blksz_y;
					x = (blk % (int)npartition_x) * blksz_x;


					p0 = (y - 1) * width + x;
					pn = (((blk % (int)npartition_x) == npartition_x - 1) ? y * width : p0 + blksz_x);
					for (p = p0;p < pn;p++)
					{
						dimgidx = p << 1;
						r = rank[dimgidx];
						canonicalize(p); canonicalize(p + width);
						connect(p, p + width, (Pixel)rankitem[r].alpha, (Imgidx)(r + imgsize));
					}
				}
				npartition_y = (npartition_y + 1) / 2;
				blksz_y <<= 1;
				if (npartition_y == 1)
					blksz_y = height;
				else
				blksz_y = _min(blksz_y, height);
			}

			if ((npartition_y == 1 || blksz_x <= blksz_y) && npartition_x > 1)
			{
				double t0 = get_cpu_time(), t;
				numblk = npartition_y * (npartition_x / 2);

				#pragma omp parallel for
				for (int blk = 0;blk < numblk;blk++)
				{
					Imgidx x, y, r, p, p0, pn, dimgidx;
					x = (1 + 2 * (blk / npartition_y)) * blksz_x;
					y = (blk % (int)npartition_y) * blksz_y;

					p0 = y * width + x - 1;
					pn = ((blk % (int)npartition_y) == npartition_y - 1) ? height * width : p0 + width * blksz_y;

					for (p = p0;p < pn;p += width)
					{
						dimgidx = (p << 1) + 1;
						r = rank[dimgidx];
						canonicalize(p); canonicalize(p + 1);
						connect(p, p + 1, (Pixel)rankitem[r].alpha, (Imgidx)(r + imgsize));
					}
				}
				npartition_x = (npartition_x + 1) / 2;
				blksz_x <<= 1;
				if (npartition_x == 1)
					blksz_x = width;
				else
					blksz_x = _min(blksz_x, width);
			}
		}
	}

	void create_queues(Trie_Cache<Imgidx, trieidx> ***queues, Imgidx nredges, _int8 npartition_x, _int8 npartition_y, _int8 listsize)
	{
			*queues = new Trie_Cache<Imgidx, trieidx>*[(int)npartition_x * (int)npartition_y];
			Trie_Cache<Imgidx, trieidx> **Q = *queues;

			Imgidx p = 0;
			for (Imgidx y = 0;y < npartition_y - 1; y++)
			{
				for (Imgidx x = 0;x < npartition_x - 1; x++)
				{
					Q[p] = new Trie_Cache<Imgidx, trieidx>(nredges + 1, listsize);
					Q[p++]->push(nredges);
				}
				Q[p] = new Trie_Cache<Imgidx, trieidx>(nredges + 1, listsize);
				Q[p++]->push(nredges);
			}
			for (Imgidx x = 0;x < npartition_x - 1; x++)
			{
				Q[p] = new Trie_Cache<Imgidx, trieidx>(nredges + 1, listsize);
				Q[p++]->push(nredges);
			}
			Q[p] = new Trie_Cache<Imgidx, trieidx>(nredges + 1, listsize);
			Q[p++]->push(nredges);
	}

	void set_subimgsizes(Imgidx** subimgsizes, _int8 npartition_x, _int8 npartition_y, _int64 blksz, _int64 blksz_lastcol, _int64 blksz_lastrow, _int64 blksz_last)
	{
		*subimgsizes = (Imgidx*)Malloc((int)npartition_x * (int)npartition_y * sizeof(Imgidx));
		Imgidx *sizes = *subimgsizes;
		Imgidx p = 0;
		for (Imgidx y = 0;y < npartition_y - 1; y++)
		{
			for (Imgidx x = 0;x < npartition_x - 1; x++)
				sizes[p++] = blksz;
			sizes[p++] = blksz_lastcol;
		}
		for (Imgidx x = 0;x < npartition_x - 1; x++)
			sizes[p++] = blksz_lastrow;
		sizes[p++] = blksz_last;
	}
};

class AlphaTree
{
	void *tree;
	_int8 imgidx, pix_type, bit_depth;
public:
	AlphaTree() : tree(0) {}
	~AlphaTree()
	{
		clear();
	}

	inline void clear()
	{
		if (tree)
		{
			if (imgidx == IMGIDX_32BITS)
			{
				if (pix_type == PIXEL_8BIT)			delete ((ATree<_int32, _uint32, _uint8>*)tree);
				else if (pix_type == PIXEL_16BIT)	delete ((ATree<_int32, _uint32, _uint16>*)tree);
				else if (pix_type == PIXEL_32BIT)	delete ((ATree<_int32, _uint32, _uint32>*)tree);
				else if (pix_type == PIXEL_64BIT)	delete ((ATree<_int32, _uint32, _uint64>*)tree);
				else if (pix_type == PIXEL_FLOAT)	delete ((ATree<_int32, _uint32, _uint32>*)tree);
				else								delete ((ATree<_int32, _uint32, _uint64>*)tree);
			}
			else
			{
				if (pix_type == PIXEL_8BIT)			delete ((ATree<_int64, _uint64, _uint8>*)tree);
				else if (pix_type == PIXEL_16BIT)	delete ((ATree<_int64, _uint64, _uint16>*)tree);
				else if (pix_type == PIXEL_32BIT)	delete ((ATree<_int64, _uint64, _uint32>*)tree);
				else if (pix_type == PIXEL_64BIT)	delete ((ATree<_int64, _uint64, _uint64>*)tree);
				else if (pix_type == PIXEL_FLOAT)	delete ((ATree<_int64, _uint64, _uint32>*)tree);
				else								delete ((ATree<_int64, _uint64, _uint64>*)tree);
			}
		}
		tree = 0;
	}

	#define BUILD_ALPHATREE(pixeltype, pixeltype_in_tree, float) void BuildAlphaTree(pixeltype *img, int height, int width, int channel, int connectivity, int algorithm, int numthreads, int tse, double fparam1, double fparam2, int iparam1) \
	{ \
		pix_type = sizeof(pixeltype); \
		if (float) \
			pix_type <<= 2; \
		bit_depth = sizeof(pixeltype) << 3; \
		/*determine index size based on the max alphatree size	*/\
		double number_of_node_max = (double)(1 + (connectivity >> 1)) * (double)width * (double)height;\
		double directional_edge = 2.0; /*some algorithms need to assign directions on edges*/\
		if (number_of_node_max * directional_edge < (double)0xefffffff) \
		{ \
			imgidx = IMGIDX_32BITS; \
			tree = new ATree<_int32, _uint32, pixeltype_in_tree>(pix_type, bit_depth); \
			((ATree<_int32, _uint32, pixeltype_in_tree>*)tree)->BuildAlphaTree((pixeltype_in_tree*)img, height, width, channel, connectivity, algorithm, numthreads, tse, fparam1, fparam2, iparam1); \
		} \
		else \
		{ \
			imgidx = IMGIDX_64BITS; \
			tree = new ATree<_int64, _uint64, pixeltype_in_tree>(pix_type, bit_depth); \
			((ATree<_int64, _uint64, pixeltype_in_tree>*)tree)->BuildAlphaTree((pixeltype_in_tree*)img, height, width, channel, connectivity, algorithm, numthreads, tse, fparam1, fparam2, iparam1); \
		} \
	}

	#define BUILD_ALPHATREE_DEFAULT(pixeltype, pixeltype_in_tree, float) void BuildAlphaTree(pixeltype *img, int height, int width, int channel) \
	{ \
		pix_type = sizeof(pixeltype); \
		if (float) \
			pix_type <<= 2; \
		bit_depth = sizeof(pixeltype) << 3; \
		int algorithm;\
		if (bit_depth <= 12)\
			algorithm = FLOOD_HIERARQUEUE;\
		else if (bit_depth < 32)\
			algorithm = FLOOD_HIERARQUEUE_CACHE;\
		else \
			algorithm = FLOOD_HIERARHEAPQUEUE_CACHE;	\
		if ((_int64)height * (_int64)width < (_int64)0xefffffff) \
		{ \
			imgidx = IMGIDX_32BITS; \
			tree = new ATree<_int32, _uint32, pixeltype_in_tree>(pix_type, bit_depth); \
			((ATree<_int32, _uint32, pixeltype_in_tree>*)tree)->BuildAlphaTree((pixeltype_in_tree*)img, height, width, channel, 4, algorithm, 0, 0, 0); \
		} \
		else \
		{ \
			imgidx = IMGIDX_64BITS; \
			tree = new ATree<_int64, _uint64, pixeltype_in_tree>(pix_type, bit_depth); \
			((ATree<_int64, _uint64, pixeltype_in_tree>*)tree)->BuildAlphaTree((pixeltype_in_tree*)img, height, width, channel, 4, algorithm, 0, 0, 0); \
		} \
	}

	BUILD_ALPHATREE(_uint8, _uint8, 0)
	BUILD_ALPHATREE(_uint16, _uint16, 0)
	BUILD_ALPHATREE(_uint32, _uint32, 0)
	BUILD_ALPHATREE(_uint64, _uint64, 0)
	BUILD_ALPHATREE(float, _uint32, 1)
	BUILD_ALPHATREE(double, _uint64, 1)
	BUILD_ALPHATREE_DEFAULT(_uint8, _uint8, 0)
	BUILD_ALPHATREE_DEFAULT(_uint16, _uint16, 0)
	BUILD_ALPHATREE_DEFAULT(_uint32, _uint32, 0)
	BUILD_ALPHATREE_DEFAULT(_uint64, _uint64, 0)
	BUILD_ALPHATREE_DEFAULT(float, _uint32, 1)
	BUILD_ALPHATREE_DEFAULT(double, _uint64, 1)

	void AreaFilter(double *outimg, double area)
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)			((ATree<_int32, _uint32, _uint8>*)tree)->AreaFilter((double*)outimg, area);
			else if (pix_type == PIXEL_16BIT)	((ATree<_int32, _uint32, _uint16>*)tree)->AreaFilter((double*)outimg, area);
			else if (pix_type == PIXEL_32BIT)	((ATree<_int32, _uint32, _uint32>*)tree)->AreaFilter((double*)outimg, area);
			else if (pix_type == PIXEL_64BIT)	((ATree<_int32, _uint32, _uint64>*)tree)->AreaFilter((double*)outimg, area);
			else if (pix_type == PIXEL_FLOAT)	((ATree<_int32, _uint32, _uint32>*)tree)->AreaFilter((double*)outimg, area);
			else								((ATree<_int32, _uint32, _uint64>*)tree)->AreaFilter((double*)outimg, area);
		}
		else
		{
			if (pix_type == PIXEL_8BIT)			((ATree<_int64, _uint64, _uint8>*)tree)->AreaFilter((double*)outimg, area);
			else if (pix_type == PIXEL_16BIT)	((ATree<_int64, _uint64, _uint16>*)tree)->AreaFilter((double*)outimg, area);
			else if (pix_type == PIXEL_32BIT)	((ATree<_int64, _uint64, _uint32>*)tree)->AreaFilter((double*)outimg, area);
			else if (pix_type == PIXEL_64BIT)	((ATree<_int64, _uint64, _uint64>*)tree)->AreaFilter((double*)outimg, area);
			else if (pix_type == PIXEL_FLOAT)	((ATree<_int64, _uint64, _uint32>*)tree)->AreaFilter((double*)outimg, area);
			else								((ATree<_int64, _uint64, _uint64>*)tree)->AreaFilter((double*)outimg, area);
		}
	}

	void AlphaFilter(double *outimg, double alpha)
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)			((ATree<_int32, _uint32, _uint8>*)tree)->AreaFilter((double*)outimg, alpha);
			else if (pix_type == PIXEL_16BIT)	((ATree<_int32, _uint32, _uint16>*)tree)->AreaFilter((double*)outimg, alpha);
			else if (pix_type == PIXEL_32BIT)	((ATree<_int32, _uint32, _uint32>*)tree)->AreaFilter((double*)outimg, alpha);
			else if (pix_type == PIXEL_64BIT)	((ATree<_int32, _uint32, _uint64>*)tree)->AreaFilter((double*)outimg, alpha);
			else if (pix_type == PIXEL_FLOAT)	((ATree<_int32, _uint32, _uint32>*)tree)->AreaFilter((double*)outimg, alpha);
			else								((ATree<_int32, _uint32, _uint64>*)tree)->AreaFilter((double*)outimg, alpha);
		}
		else
		{
			if (pix_type == PIXEL_8BIT)			((ATree<_int64, _uint64, _uint8>*)tree)->AreaFilter((double*)outimg, alpha);
			else if (pix_type == PIXEL_16BIT)	((ATree<_int64, _uint64, _uint16>*)tree)->AreaFilter((double*)outimg, alpha);
			else if (pix_type == PIXEL_32BIT)	((ATree<_int64, _uint64, _uint32>*)tree)->AreaFilter((double*)outimg, alpha);
			else if (pix_type == PIXEL_64BIT)	((ATree<_int64, _uint64, _uint64>*)tree)->AreaFilter((double*)outimg, alpha);
			else if (pix_type == PIXEL_FLOAT)	((ATree<_int64, _uint64, _uint32>*)tree)->AreaFilter((double*)outimg, alpha);
			else								((ATree<_int64, _uint64, _uint64>*)tree)->AreaFilter((double*)outimg, alpha);
		}
	}

	void print_tree()
	{
		if (imgidx == IMGIDX_32BITS)
		{
			if (pix_type == PIXEL_8BIT)			((ATree<_int32, _uint32, _uint8>*)tree)->print_tree();
			else if (pix_type == PIXEL_16BIT)	((ATree<_int32, _uint32, _uint16>*)tree)->print_tree();
			else if (pix_type == PIXEL_32BIT)	((ATree<_int32, _uint32, _uint32>*)tree)->print_tree();
			else if (pix_type == PIXEL_64BIT)	((ATree<_int32, _uint32, _uint64>*)tree)->print_tree();
			else if (pix_type == PIXEL_FLOAT)	((ATree<_int32, _uint32, _uint32>*)tree)->print_tree();
			else								((ATree<_int32, _uint32, _uint64>*)tree)->print_tree();
		}
		else
		{
			if (pix_type == PIXEL_8BIT)			((ATree<_int64, _uint64, _uint8>*)tree)->print_tree();
			else if (pix_type == PIXEL_16BIT)	((ATree<_int64, _uint64, _uint16>*)tree)->print_tree();
			else if (pix_type == PIXEL_32BIT)	((ATree<_int64, _uint64, _uint32>*)tree)->print_tree();
			else if (pix_type == PIXEL_64BIT)	((ATree<_int64, _uint64, _uint64>*)tree)->print_tree();
			else if (pix_type == PIXEL_FLOAT)	((ATree<_int64, _uint64, _uint32>*)tree)->print_tree();
			else								((ATree<_int64, _uint64, _uint64>*)tree)->print_tree();
		}
	}

	#define GET_MEMBER(returntype, methodname, member) returntype methodname()\
	{\
		if (imgidx == IMGIDX_32BITS)\
		{\
			if (pix_type == PIXEL_8BIT)			return (returntype)((ATree<_int32, _uint32, _uint8>*)tree)->member;\
			else if (pix_type == PIXEL_16BIT)	return (returntype)((ATree<_int32, _uint32, _uint16>*)tree)->member;\
			else if (pix_type == PIXEL_32BIT)	return (returntype)((ATree<_int32, _uint32, _uint32>*)tree)->member;\
			else if (pix_type == PIXEL_64BIT)	return (returntype)((ATree<_int32, _uint32, _uint64>*)tree)->member;\
			else if (pix_type == PIXEL_FLOAT)	return (returntype)((ATree<_int32, _uint32, _uint32>*)tree)->member;\
			else								return (returntype)((ATree<_int32, _uint32, _uint64>*)tree)->member;\
		}\
		else\
		{\
			if (pix_type == PIXEL_8BIT)			return (returntype)((ATree<_int64, _uint64, _uint8>*)tree)->member;\
			else if (pix_type == PIXEL_16BIT)	return (returntype)((ATree<_int64, _uint64, _uint16>*)tree)->member;\
			else if (pix_type == PIXEL_32BIT)	return (returntype)((ATree<_int64, _uint64, _uint32>*)tree)->member;\
			else if (pix_type == PIXEL_64BIT)	return (returntype)((ATree<_int64, _uint64, _uint64>*)tree)->member;\
			else if (pix_type == PIXEL_FLOAT)	return (returntype)((ATree<_int64, _uint64, _uint32>*)tree)->member;\
			else								return (returntype)((ATree<_int64, _uint64, _uint64>*)tree)->member;\
		}\
	}

	GET_MEMBER(_int64, get_maxsize, maxSize)
	GET_MEMBER(_int64, get_cursize, curSize)
	GET_MEMBER(double, get_nrmsd, nrmsd)
//	GET_MEMBER(auto, get_node, node)
};

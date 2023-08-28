#include <ctime>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <cstring>
#include "defines.h"
#include "AlphaTree.h"
#include "walltime.h"
#include "pgmio.h"

using namespace std;

#define RANDOMINPUT		 0
#define NUM_ALGORITHMS 12

#define DEBUG 0

#define OUTPUT_FNAME "./AlphaTree.dat"
#define OUTIMG_FNAME "./outimg.jpg"


_uint64 *qrecord;



#if DEBUG
void* buf;
_uint64 bufsize;
void save_buf(void* src, _uint64 size)
{
	memcpy(buf, src, size);
	bufsize = size;
}
_uint8 isChanged(void *src)
{
	_uint64 i;
	for (i = 0; i < bufsize; i++)
	{
		if (((_uint8*)buf)[i] != ((_uint8*)src)[i])
			return 1;
	}
	return 0;
}
#endif

void RandomizedHDRimage(_uint64* hdrimg, _uint8* ldrimg, _int64 imgsize)
{
	_uint64 pix;

	for (_int64 i = 0; i < imgsize; i++)
	{
		pix = ((_uint64)ldrimg[i]) << 56;
		pix |= ((_uint64)(rand() & 0xff) << 48);
		pix |= ((_uint64)(rand() & 0xff) << 40);
		pix |= ((_uint64)(rand() & 0xff) << 32);
		pix |= ((_uint64)(rand() & 0xff) << 24);
		pix |= ((_uint64)(rand() & 0xff) << 16);
		pix |= ((_uint64)(rand() & 0xff) << 8);
		pix |= ((_uint64)(rand() & 0xff));
		hdrimg[i] = pix;
	}
}

void RandomizedHDRimage(_uint64* hdrimg, _int64 imgsize)
{
	_uint64 pix;

	for (_int64 i = 0; i < imgsize; i++)
	{
		//pix = ((_uint64)ldrimg[i]) << 56;
		pix = ((_uint64)(rand() & 0xff) << 56); //tmp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		pix |= ((_uint64)(rand() & 0xff) << 48);
		pix |= ((_uint64)(rand() & 0xff) << 40);
		pix |= ((_uint64)(rand() & 0xff) << 32);
		pix |= ((_uint64)(rand() & 0xff) << 24);
		pix |= ((_uint64)(rand() & 0xff) << 16);
		pix |= ((_uint64)(rand() & 0xff) << 8);
		pix |= ((_uint64)(rand() & 0xff));

		hdrimg[i] = pix;
	}
}

void RandomizedHDRimage(_uint32* hdrimg, _uint8* ldrimg, _int64 imgsize)
{
	_uint32 pix;

	for (_int64 i = 0; i < imgsize; i++)
	{
		pix = ((_uint64)ldrimg[i]) << 24;
		pix |= ((_uint64)(rand() & 0xff) << 16);
		pix |= ((_uint64)(rand() & 0xff) << 8);
		pix |= ((_uint64)(rand() & 0xff));
		hdrimg[i] = pix;
	}
}

/*
void DeleteAlphaTree(AlphaTree* tree)
{
	Free(tree->parentAry);
	Free(tree->node);
	Free(tree);
}*/


void Randomizedimage(_uint8*& img, _int64 imgsize, int bit_depth, int ch)
{
	_uint64 pix;
	int shamt = 8 - bit_depth;

	if(img == 0) img = new _uint8[imgsize * ch];

	for (_int64 i = 0; i < imgsize * ch; i++)
	{
		pix = ((_uint8)(rand() & 0xff));

		img[i] = pix >> shamt;
	}
}

void Randomizedimage(_uint16*& img, _int64 imgsize, int bit_depth, int ch)
{
	_uint64 pix;
	int shamt = 16 - bit_depth;

	if(img == 0) img = new _uint16[imgsize * ch];

	for (_int64 i = 0; i < imgsize * ch; i++)
	{
		pix = ((_uint64)(rand() & 0xff) << 8);
		pix |= ((_uint64)(rand() & 0xff));

		img[i] = pix >> shamt;
	}
}

void Randomizedimage(_uint32*& img, _int64 imgsize, int bit_depth, int ch)
{
	_uint64 pix;
	int shamt = 32 - bit_depth;

	if(img == 0) img = new _uint32[imgsize * ch];

	for (_int64 i = 0; i < imgsize * ch; i++)
	{
		pix = ((_uint64)(rand() & 0xff) << 24);
		pix |= ((_uint64)(rand() & 0xff) << 16);
		pix |= ((_uint64)(rand() & 0xff) << 8);
		pix |= ((_uint64)(rand() & 0xff));

		img[i] = pix >> shamt;
	}
}

void Randomizedimage(_uint64*& img, _int64 imgsize, int bit_depth, int ch)
{
	_uint64 pix;
	int shamt = 32 - bit_depth;

	if(img == 0) img = new _uint64[imgsize * ch];

	for (_int64 i = 0; i < imgsize * ch; i++)
	{
		pix = ((_uint64)(rand() & 0xff) << 56);
		pix = ((_uint64)(rand() & 0xff) << 48);
		pix = ((_uint64)(rand() & 0xff) << 40);
		pix = ((_uint64)(rand() & 0xff) << 32);
		pix = ((_uint64)(rand() & 0xff) << 24);
		pix |= ((_uint64)(rand() & 0xff) << 16);
		pix |= ((_uint64)(rand() & 0xff) << 8);
		pix |= ((_uint64)(rand() & 0xff));

		img[i] = pix >> shamt;
	}
}

template<class T>
T getRand() {
	T ret;
	switch (sizeof(T)) {
		case (1):
			ret = (T)(rand() & 0xff);
			break;
		case (2):
			ret = (T)(rand() & 0xff);
			ret |= (T)((rand()  & 0xff) << 8);
			break;
		case (4):
			ret = (T)(rand() & 0xff);
			ret |= (T)((rand() & 0xff) << 8);
			ret |= (T)((rand() & 0xff) << 16);
			ret |= (T)((rand() & 0xff) << 24);
			break;
		case (8):
			ret = (T)(rand() & 0xff);
			ret |= (T)((_uint64)(rand() & 0xff) << 8);
			ret |= (T)((_uint64)(rand() & 0xff) << 16);
			ret |= (T)((_uint64)(rand() & 0xff) << 24);
			ret |= (T)((_uint64)(rand() & 0xff) << 32);
			ret |= (T)((_uint64)(rand() & 0xff) << 40);
			ret |= (T)((_uint64)(rand() & 0xff) << 48);
			ret |= (T)((_uint64)(rand() & 0xff) << 56);
			break;
		default:
			ret = 0;
			break;	
	}
	return ret;
}

template<class T>
T* getRandomizedImage(_int64 imgsize, int bit_depth, int ch)
{
	if ((int)(sizeof(T) * 8) < bit_depth)
		return NULL;
	int shamt = (sizeof(T) * 8) - bit_depth;

	T* img = new T[imgsize * ch];

	for (_int64 i = 0; i < imgsize * ch; i++)
	{
		T pix = getRand<T>();
		img[i] = pix >> shamt;
	}
	return img;
}

void imageblur(_uint64** img, int height, int width, int nch)
{
	int pidx = 0;
	int up,down,left,right, num_neighbor;
	double pixsum;

	_uint64* blurimg = new _uint64[height * width * nch];

	for(int ch = 0;ch < nch;ch++)
	{
		int offset = ch * height * width;
		for(int h = 0;h < height;h++)
		{
			up = h > 0;
			down = h < height - 1;
			for(int w = 0;w < width;w++)
			{
				left = w > 0;
				right = w < width -1;

				num_neighbor = up + down + left + right;
				pixsum = 0;
				if(left) pixsum += (*img)[pidx-1+offset];
				if(right) pixsum += (*img)[pidx+1+offset];
				if(up) pixsum += (*img)[pidx-width+offset];
				if(down) pixsum += (*img)[pidx+width+offset];
				pixsum += (double)(*img)[pidx+offset] * 2.0;

				blurimg[pidx+offset] = (_uint64)(pixsum / (2.0 + num_neighbor));
				pidx++;
			}
		}
	}

	_uint64 *tmp = *img;
	*img = blurimg;
	delete[] tmp;
}

void imageblur(_uint32** img, int height, int width, int nch)
{
	int pidx = 0;
	int up,down,left,right, num_neighbor;
	double pixsum;

	_uint32* blurimg = new _uint32[height * width * nch];

	for(int ch = 0;ch < nch;ch++)
	{
		int offset = ch * height * width;
		for(int h = 0;h < height;h++)
		{
			up = h > 0;
			down = h < height - 1;
			for(int w = 0;w < width;w++)
			{
				left = w > 0;
				right = w < width -1;

				num_neighbor = up + down + left + right;
				pixsum = 0;
				if(left) pixsum += (*img)[pidx-1+offset];
				if(right) pixsum += (*img)[pidx+1+offset];
				if(up) pixsum += (*img)[pidx-width+offset];
				if(down) pixsum += (*img)[pidx+width+offset];
				pixsum += (double)(*img)[pidx+offset] * 2.0;

				blurimg[pidx+offset] = (_uint32)(pixsum / (2.0 + num_neighbor));
				pidx++;
			}
		}
	}

	_uint32 *tmp = *img;
	*img = blurimg;
	delete[] tmp;
}

void imageblur(_uint16** img, int height, int width, int nch)
{
	int pidx = 0;
	int up,down,left,right, num_neighbor;
	double pixsum;

	_uint16* blurimg = new _uint16[height * width * nch];

	for(int ch = 0;ch < nch;ch++)
	{
		int offset = ch * height * width;
		for(int h = 0;h < height;h++)
		{
			up = h > 0;
			down = h < height - 1;
			for(int w = 0;w < width;w++)
			{
				left = w > 0;
				right = w < width -1;

				num_neighbor = up + down + left + right;
				pixsum = 0;
				if(left) pixsum += (*img)[pidx-1+offset];
				if(right) pixsum += (*img)[pidx+1+offset];
				if(up) pixsum += (*img)[pidx-width+offset];
				if(down) pixsum += (*img)[pidx+width+offset];
				pixsum += (double)(*img)[pidx+offset] * 2.0;

				blurimg[pidx+offset] = (_uint16)(pixsum / (2.0 + num_neighbor));
				pidx++;
			}
		}
	}

	_uint16 *tmp = *img;
	*img = blurimg;
	delete[] tmp;
}

void imageblur(_uint8** img, int height, int width, int nch)
{
	int pidx = 0;
	int up,down,left,right, num_neighbor;
	double pixsum;

	_uint8* blurimg = new _uint8[height * width * nch];

	for(int ch = 0;ch < nch;ch++)
	{
		int offset = ch * height * width;
		for(int h = 0;h < height;h++)
		{
			up = h > 0;
			down = h < height - 1;
			for(int w = 0;w < width;w++)
			{
				left = w > 0;
				right = w < width -1;

				num_neighbor = up + down + left + right;
				pixsum = 0;
				if(left) pixsum += (*img)[pidx-1+offset];
				if(right) pixsum += (*img)[pidx+1+offset];
				if(up) pixsum += (*img)[pidx-width+offset];
				if(down) pixsum += (*img)[pidx+width+offset];
				pixsum += (double)(*img)[pidx+offset] * 2.0;

				blurimg[pidx+offset] = (_uint8)(pixsum / (2.0 + num_neighbor));
				pidx++;
			}
		}
	}

	_uint8 *tmp = *img;
	*img = blurimg;
	delete[] tmp;
}

void sort_h1(double* h1_arr, int obs)
{
	int *hist, *h, *h1, i, j, bitfield, shamt, hsum;
	double *tmp, *p, *end;
	int64_t mask = 0xffff, num, offset1, offset2, offset3;

	offset1 = 65536;
	offset2 = 65536 << 1;
	offset3 = offset1 + offset2;

	hist = new int[65536 << 2];
	tmp = new double[obs];

	for (i = 0; i < 65536 << 2; i++)
		hist[i] = 0;

	end = h1_arr + obs;
	for (p = h1_arr; p < end; p++)
	{
		num = *((int64_t*)(p));

		hist[num & mask]++;
		hist[offset1 + ((num >> 16) & mask)]++;
		hist[offset2 + ((num >> 32) & mask)]++;
		hist[offset3 + ((num >> 48) & mask)]++;
	}

	shamt = 0;
	for (bitfield = 0; bitfield < 4; bitfield++)
	{
		hsum = 0;
		h = hist;
		h1 = h + 65536;
		while (h != h1)
		{
			hsum += *h;
			*(h++) = hsum;
		}

		for (p = h1_arr + obs - 1; p >= h1_arr; p--)
		{
			num = *((int64_t*)p);
			j = --hist[(num >> shamt) & mask];
			tmp[j] = *p;
		}
		p = tmp; tmp = h1_arr; h1_arr = p;

		hist += 65536;
		shamt += 16;
	}
	hist -= 65536 << 2;

	delete[] hist;
	delete[] tmp;
}

void alg_name(char*dst, int alg)
{
	switch(alg)
	{
		case(UNIONFIND):
			strcpy(dst,"Unionfind (Berger) (UNIONFIND)");
			break;
		case(FLOOD_HIERARQUEUE):
			strcpy(dst,"Flood using Hierarchical queue (Salembier) (FLOOD_HIERARQUEUE)");
			break;
		case(FLOOD_HIERARQUEUE_CACHE):
			strcpy(dst,"Flood using Hierarchical queue with cache (Salembier) (FLOOD_HIERARQUEUE_CACHE)");
			break;
		case(FLOOD_HIERARHEAPQUEUE_CACHE):
			strcpy(dst,"Flood using hierarchical heap queue with cache (FLOOD_HIERARHEAPQUEUE_CACHE)");
			break;
		case(FLOOD_HIERARHEAPQUEUE_CACHE_HISTEQ):
			strcpy(dst,"Flood using hist-equalized hierarchical heap queue with cache (FLOOD_HIERARHEAPQUEUE_CACHE_HISTEQ)");
			break;
		case(FLOOD_HEAPQUEUE):
			strcpy(dst,"Flood using heap queue (Wilkinson) (FLOOD_HEAPQUEUE)");
			break;
		case(FLOOD_HEAPQUEUE_CACHE):
			strcpy(dst,"Flood using quad heap queue with cache (Wilkinson) (FLOOD_HEAPQUEUE_CACHE)");
			break;
		case(FLOOD_TRIE):
			strcpy(dst,"Flood using Trie (Teeninga) (FLOOD_TRIE)");
			break;
		case(FLOOD_TRIE_CACHE):
			strcpy(dst,"Flood using cached trie (FLOOD_TRIE_CACHE)");
			break;
		case(FLOOD_TRIE_HYPERGRAPH):
			strcpy(dst,"Flood using Trie and hypergraph (FLOOD_TRIE_HYPERGRAPH)");
			break;
		case(PILOT_RANK):
			strcpy(dst,"Pilot-tree parallel (Ugo) (PILOT_RANK)");
			break;
		case(FLOOD_HIERARQUEUE_HYPERGRAPH):
			strcpy(dst,"Flooding using HierarQueue, Hypergraph (FLOOD_HIERARQUEUE_HYPERGRAPH)");
			break;
		case(FLOOD_HIERARQUEUE_PAR):
			strcpy(dst,"Block-based parallel using Hierarqueue (FLOOD_HQUEUE_PAR)");
			break;
		case(FLOOD_HIERARHEAPQUEUE):
			strcpy(dst,"Flood using hierarchical heap queue (FLOOD_HIERARHEAPQUEUE)");
			break;
		default:
			strcpy(dst,"NO_NAME");
			break;
	}
}

bool is_par(int alg)
{
	switch(alg)
	{
		case(UNIONFIND):
			return 0;
			break;
		case(FLOOD_HIERARQUEUE):
			return 0;
			break;
		case(FLOOD_HIERARQUEUE_CACHE):
			return 0;
			break;
		case(FLOOD_HEAPQUEUE_CACHE):
			return 0;
			break;
		case(FLOOD_TRIE):
			return 0;
			break;
		case(FLOOD_TRIE_CACHE):
			return 0;
			break;
		case(FLOOD_TRIE_HYPERGRAPH):
			return 0;
			break;
		case(PILOT_RANK):
			return 1;
			break;
		case(FLOOD_HIERARQUEUE_HYPERGRAPH):
			return 0;
			break;
		case(FLOOD_HIERARQUEUE_PAR):
			return 1;
			break;
		default:
			return 0;
			break;
	}
}

_uint64 rand64()
{
	_uint64 ret, num;
	ret = 0;
	num = (_uint64)(rand()%256) << 56;
	ret |= num;
	num = (_uint64)(rand()%256) << 48;
	ret |= num;
	num = (_uint64)(rand()%256) << 40;
	ret |= num;
	num = (_uint64)(rand()%256) << 32;
	ret |= num;
	num = (_uint64)(rand()%256) << 24;
	ret |= num;
	num = (_uint64)(rand()%256) << 16;
	ret |= num;
	num = (_uint64)(rand()%256) << 8;
	ret |= num;
	num = (_uint64)(rand()%256);
	ret |= num;

	return ret;
}


template <class Pixel>
void adjust_bitdepth(Pixel *img, int imgsize, int bitdepth)
{
	int shamt = (sizeof(Pixel) * 8) - bitdepth;
	for(int i = 0;i < imgsize;i++)
		img[i] = img[i] >> shamt;
}

void adjust_bitdepth(_uint8 *img8, _uint16 *img, int imgsize, int bitdepth)
{
	if(bitdepth == 8)
	{
		int shamt = 16 - bitdepth;
		for(int i = 0;i < imgsize;i++)
			img8[i] = img[i] >> shamt;
	}
	else
	{
		int shamt = 16 - bitdepth;
		for(int i = 0;i < imgsize;i++)
			img[i] = img[i] >> shamt;
	}
}

template <class Pixel>
void imread(Pixel*& img, char *filename, int imgidx, int& height, int& width, int nch, int bitdepth)
{
	char fname[128];
	if(nch == 1)
	{
		if(imgidx + 1 < 10)
			sprintf(fname, "%s0%d_gray.pgm%c",filename, imgidx + 1,'\0');
		else
			sprintf(fname, "%s%d_gray.pgm%c",filename, imgidx + 1,'\0');
	}
	else
	{
		if(imgidx + 1 < 10)
			sprintf(fname, "%s0%d_$.pgm%c",filename, imgidx + 1,'\0');
		else
			sprintf(fname, "%s%d_$.pgm%c",filename, imgidx + 1,'\0');
	}

	printf("Reading pgm file...\n");
	char *p = strstr(fname,"$");
	int gmax;
	if(p) //color image
	{
		Pixel *g;
		*p = 'r';
		read_pgm(fname, height, width, gmax, &g);
		img = new Pixel[height * width * nch];
		memcpy(img, g, height * width * sizeof(Pixel));
		delete[] g;

		*p = 'g';
		read_pgm(fname, height, width, gmax, &g);
		memcpy(img + height * width, g, height * width * sizeof(Pixel));
		delete[] g;

		*p = 'b';
		read_pgm(fname, height, width, gmax, &g);
		memcpy(img + 2 * height * width, g, height * width * sizeof(Pixel));
		delete[] g;
	}
	else
	{
		read_pgm(fname, height, width, gmax, &img);
	}
	printf("imread: %s\n", fname);
}

struct AtreeInputParams
{
	char *name;
	int nchannels;
	int numthreads;
	int testimgsize;
	int algorithmcode;
	int bitdepth;
	int tse;
	int connectivity;
	int numitr;
	double fparam1;
	double fparam2;
	int iparam1;
	char *fnameheader;
};

void parse_input(AtreeInputParams& input, int argc, char **argv)
{
	input.name = argv[1];
	input.nchannels = atoi(argv[2]);
	input.numthreads = atoi(argv[3]);
	input.testimgsize = atoi(argv[4]);
	input.algorithmcode = atoi(argv[5]);
	input.bitdepth = atoi(argv[6]);
	input.tse = atoi(argv[7]);
	input.fnameheader = argv[8];
	if(argc >= 10)
		input.connectivity = atoi(argv[9]);
	else
		input.connectivity = 4;
	input.numitr = atoi(argv[10]);
	input.fparam1 = atof(argv[11]);
	input.fparam2 = atof(argv[12]);
	input.iparam1 = atoi(argv[13]);
}

//args: Filename, nchannels, numthreads, testimgsize, algorithmcode, bitdepth, tseflag
int main(int argc, char **argv)
{
	if(argc == 1)
	{
		printf("args: Filename, nchannels, numthreads, testimgsize, algorithmcode, bitdepth, tseflag, outfnameheader connectivity float_param1 float_param2 int_param1\n");
		printf("No Arg - running default tests\n");

		int width = 4;
		int height = 4;

		_uint32 *img = getRandomizedImage<_uint32>(width * height, 32, 1);
	
		AlphaTree tree;

		tree.BuildAlphaTree(img, height, width, 1, 4, FLOOD_HIERARHEAPQUEUE_CACHE, 1, 0, 0.0, 0.0, 0);

		// return -1;
	}

	AtreeInputParams input;
	parse_input(input, argc, argv);
	int randomimg = (strcmp(input.name,"rand") == 0);
	double meanrunspeed[16] = {0,}, maxmemuse[16] = {0,};
	int nthr[] = {1, 2, 4, 8, 16, 32, 48, 64, 96, 128, 192, 256, 480, 960, 1920};
	char algname[256];
	char outfname[128];

	int numimg = randomimg ? 1 : 10;
	//int numitr = 10;

	int imgidxstart = 0;
	int thrstart = 0;

	srand(time(NULL));

	//if(!randomimg)
	{
		sprintf(outfname, "%s_ch%d_nthr%d_alg%d_bit%d_tse%d_conn%d.txt", input.fnameheader, input.nchannels, input.numthreads, input.algorithmcode, input.bitdepth, input.tse, input.connectivity);
		printf("output file name = %s\n",outfname);
		fstream f(outfname);
		if(0 && f.good())
		{
			char buf[128];
			ifstream fin(outfname);

			//f.seekg(0, std::ios::beg);
			fin.getline(buf,128,' ');
			imgidxstart = atoi(buf);
			fin.getline(buf,128,'\n');
			thrstart = atoi(buf);
			int numthr = randomimg ? 1 : input.numthreads;
			for(int thridx = 0;thridx < numthr;thridx++)
			{
				fin.getline(buf,128,' ');
				meanrunspeed[thridx] = atof(buf);
				fin.getline(buf,128,'\n');
				maxmemuse[thridx] = atof(buf);
			}
			printf("saved progress loaded (%d/%d)\n",imgidxstart, numimg);
			fin.close();
		}
		else
		{

		}
	}


	for(int imgidx = imgidxstart;imgidx < numimg;imgidx++)
	{
		int height, width;
		_uint16 *img;
		_uint8 *img8 = 0;

		height = width = input.testimgsize;
		int image_number = imgidx;
		//if(imgidx >= 4) image_number += 2;
		//if(imgidx >= 8) image_number += 1;
		//if(imgidx >= 9) image_number += 1;


		printf("========================================================================================\n");
		printf("========== image [%d/%d]: imgsize = %d x %d (%d bits, %d channels) ================\n",
		imgidx + 1, numimg, (int)height, (int)width, (int)input.bitdepth, (int)input.nchannels);
		printf("========================================================================================\n");
		int thritr = randomimg ? 1 : input.numthreads;
		for(int thridx = thrstart;thridx < thritr;thridx++)
		{
			int numthreads = randomimg ? input.numthreads : nthr[thridx];
			alg_name(algname, (int)input.algorithmcode);
			printf("-----------------------------------------------------------------------------------\n");
			printf("%d Running %s (%d threads)\n", (int)input.algorithmcode, algname, (int)numthreads);
			printf("-----------------------------------------------------------------------------------\n");
			double minruntime = 0;

			for (int testrep = 0; testrep < input.numitr; testrep++)
			{
				double t = get_cpu_time();
				double runtime;// = get_cpu_time() - t;
				//t = omp_get_wtime();

				if(!randomimg)
					imread(img, input.name, image_number, height, width, input.nchannels, input.bitdepth);//skip image 5 (unknown bug)

				if(!randomimg && input.bitdepth < 16)
				{
					if(input.bitdepth == 8)
						img8 = new _uint8[height * width];
					adjust_bitdepth(img8, img, height * width, input.bitdepth);
					delete[] img;
					img = 0;
				}

				AlphaTree *tree = new AlphaTree;

				if(randomimg)
				{
					img8 = 0;
					img = 0;
					int bitdepth = input.bitdepth;
					_uint32 *img32 = 0;
					_uint64 *img64 = 0;
					if(bitdepth <= 8)		Randomizedimage(img8, height * width, bitdepth, 1);
					else if(bitdepth <= 16) Randomizedimage(img, height * width, bitdepth, 1);
					else if(bitdepth <= 32) Randomizedimage(img32, height * width, bitdepth, 1);
					else					Randomizedimage(img64, height * width, bitdepth, 1);

					int queueprofile = 0;
					if(queueprofile)
					{
						qrecord = 0;
						if(bitdepth <= 8)		tree->BuildAlphaTree(img8, height, width, input.nchannels, input.connectivity, input.algorithmcode, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
						else if(bitdepth <= 16) tree->BuildAlphaTree(img, height, width, input.nchannels, input.connectivity, input.algorithmcode, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
						else if(bitdepth <= 32) tree->BuildAlphaTree(img32, height, width, input.nchannels, input.connectivity, input.algorithmcode, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
						else					tree->BuildAlphaTree(img64, height, width, input.nchannels, input.connectivity, input.algorithmcode, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
						delete tree;
						tree = new AlphaTree;
					}

					t = get_cpu_time();
					if(bitdepth <= 8)		tree->BuildAlphaTree(img8, height, width, input.nchannels, input.connectivity, input.algorithmcode, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
					else if(bitdepth <= 16) tree->BuildAlphaTree(img, height, width, input.nchannels, input.connectivity, input.algorithmcode, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
					else if(bitdepth <= 32) tree->BuildAlphaTree(img32, height, width, input.nchannels, input.connectivity, input.algorithmcode, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
					else					tree->BuildAlphaTree(img64, height, width, input.nchannels, input.connectivity, input.algorithmcode, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
					runtime = get_cpu_time() - t;

					if(img32) delete[] img32;
					if(img64) delete[] img64;
				}
				else
				{
					t = get_cpu_time();
					if(img8) tree->BuildAlphaTree(img8, height, width, input.nchannels, input.connectivity, input.algorithmcode, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
					else 	 tree->BuildAlphaTree(img, height, width, input.nchannels, input.connectivity, input.algorithmcode, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
					runtime = get_cpu_time() - t;
				}

				printf("-------------------Run %d: %.3f------------------\n",(int)testrep, runtime);

				//cout << "Run " << testrep << ": " << runtime << endl;
				if(!testrep) minruntime = runtime;
				minruntime = (minruntime > runtime) ? runtime : minruntime;

				delete tree;
				if(img) delete[] img;
				if(img8) delete[] img8;
			}
			double imgsize = (double)height * (double)width;
			double runspeed = imgsize / (1e6 * minruntime);
			meanrunspeed[thridx] += runspeed;
			maxmemuse[thridx] += max_memuse / imgsize;

			if(0 && thridx < thritr - 1)
			{
				ofstream fout(outfname);
				fout << imgidx << ' ' << thridx + 1 << endl;
				for(int thridx = 0;thridx < input.numthreads;thridx++)
				{
					fout << meanrunspeed[thridx] << ' ' << maxmemuse[thridx] << endl;
				}
				fout.close();
			}
		}

		thrstart = 0;
	}
	printf("Summary ======================================\n");
	int numthrs = randomimg ? 1 : input.numthreads;
	for(int thridx = 0;thridx < numthrs;thridx++)
	{
		printf("#thr = %d: %.3fMpix/s / %.3fB/pix \n", nthr[thridx], meanrunspeed[thridx]/numimg, maxmemuse[thridx]/numimg);
	}
	printf("==============================================\n");


	return 0;
}

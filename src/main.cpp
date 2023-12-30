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


// tmp
#include "LadderQueue.hpp"
#include "HeapQueue.h"

#define OUTPUT_FNAME "./AlphaTree.dat"
#define OUTIMG_FNAME "./outimg.jpg"

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
	int shamt = 64 - bit_depth;

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
	char *name = nullptr;
	int nchannels = 1;
	int numthreads = 1;
	int randimgsize = 100;
	char *algorithmname = nullptr;
	int bitdepth = 8;
	int tse = 1;
	int connectivity = 4;
	int numitr = 1;
	double fparam1 = 0.0;
	double fparam2 = 0.0;
	int iparam1 = 0;
	char *fnameheader = nullptr;
};

void parse_input(AtreeInputParams& input, int argc, char **argv)
{
	if(argc >= 2) input.name = argv[1];
	if(argc >= 3) input.algorithmname = argv[2];
	if(argc >= 4) input.randimgsize = atoi(argv[3]);
	if(argc >= 5) input.nchannels = atoi(argv[4]);
	if(argc >= 6) input.bitdepth = atoi(argv[5]);
	if(argc >= 7) input.numthreads = atoi(argv[6]);
	if(argc >= 8) input.tse = atoi(argv[7]);
	if(argc >= 9) input.fnameheader = argv[8];
	if(argc >= 10) input.connectivity = atoi(argv[9]);
	if(argc >= 11) input.numitr = atoi(argv[10]);
	if(argc >= 12) input.fparam1 = atof(argv[11]);
	if(argc >= 13) input.fparam2 = atof(argv[12]);
	if(argc >= 14) input.iparam1 = atoi(argv[13]);
}

int main(int argc, char **argv)
{
	if(argc == 2 && strcmp(argv[1], "help") == 0)
	{
		printf("Arguments (all optional) : imageFileName algorithmName randImgSize numChannels bitDepth numThreads runTSE outFileHeader connectivity numIteration floatParam1 floatParam2 intParam1\n");
		printf("Only .pgm files can be used for imageFileName. Type \"rand\" for imageFileName to run on randomly-generated images.\n");
		printf("algorithmName list\n");
		for (int algorithmCode = 0;algorithmCode < AlphaTree<uint8_t>::NUM_ALGORITHMS;algorithmCode++)
		{
			char algorithmName[256];
			AlphaTree<uint8_t>::getAlgorithmName(algorithmName, algorithmCode);
			char algorithmDesc[256];
			AlphaTree<uint8_t>::getAlgorithmDescriptionFromCode(algorithmDesc, algorithmCode);
			printf("%d) %s: %s\n", algorithmCode, algorithmName, algorithmDesc);
		}

		return 0;
	}

	using namespace std;

	AtreeInputParams input;
	parse_input(input, argc, argv);
	bool randomimg = (input.name == nullptr) || (strcmp(input.name,"rand") == 0);
	double meanrunspeed[16] = {0,}, maxmemuse[16] = {0,};
	int nthr[] = {1, 2, 4, 8, 16, 32, 48, 64, 96, 128, 192, 256, 480, 960, 1920};
	
	char outfname[128];

	int numimg = randomimg ? 1 : 10;
	int imgidxstart = 0;
	int thrstart = 0;
	
	if(!randomimg && input.fnameheader)
	{
		sprintf(outfname, "%s_ch%d_nthr%d_%s_bit%d_tse%d_conn%d.txt", input.fnameheader, input.nchannels, input.numthreads, input.algorithmname, input.bitdepth, input.tse, input.connectivity);
		printf("output file name = %s\n",outfname);
		fstream f(outfname);
		if(0 && f.good())
		{
			char buf[128];
			ifstream fin(outfname);

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

		height = width = input.randimgsize;
		int image_number = imgidx;


		printf("========================================================================================\n");
		printf("========== image [%d/%d]: imgsize = %d x %d (%d bits, %d channels) ================\n",
		imgidx + 1, numimg, (int)height, (int)width, (int)input.bitdepth, (int)input.nchannels);
		printf("========================================================================================\n");
		int thritr = randomimg ? 1 : input.numthreads;
		for(int thridx = thrstart;thridx < thritr;thridx++)
		{
			int numthreads = randomimg ? input.numthreads : nthr[thridx];
			char algname[256];
			if (input.algorithmname)
				AlphaTree<uint8_t>::getAlgorithmDescription(algname, input.algorithmname);
			else 
				AlphaTree<uint8_t>::getAlgorithmDescription(algname, "0");
			printf("-----------------------------------------------------------------------------------\n");
			printf("Running %s (%d threads)\n", algname, (int)numthreads);
			printf("-----------------------------------------------------------------------------------\n");
			double minruntime = 0;

			for (int testrep = 0; testrep < input.numitr; testrep++)
			{
				double t = 0.0;
				double runtime = 0.0;// = get_cpu_time() - t;
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


					t = get_wall_time();
					if(bitdepth <= 8) {
						AlphaTree<_uint8> *tree = new AlphaTree<_uint8>;
						tree->BuildAlphaTree(img8, height, width, input.nchannels, input.connectivity, input.algorithmname, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
						delete tree;
					} else if(bitdepth <= 16) {
						AlphaTree<_uint16> *tree = new AlphaTree<_uint16>;
						tree->BuildAlphaTree(img, height, width, input.nchannels, input.connectivity, input.algorithmname, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
						delete tree;
					} else if(bitdepth <= 32) {
						AlphaTree<_uint32> *tree = new AlphaTree<_uint32>;
						tree->BuildAlphaTree(img32, height, width, input.nchannels, input.connectivity, input.algorithmname, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
						delete tree;
					} else {
						AlphaTree<_uint64> *tree = new AlphaTree<_uint64>;
						tree->BuildAlphaTree(img64, height, width, input.nchannels, input.connectivity, input.algorithmname, (int)numthreads, input.tse, input.fparam1, input.fparam2, input.iparam1);
						delete tree;
					}
					runtime = get_wall_time() - t;

					if(img32) delete[] img32;
					if(img64) delete[] img64;
				}
				printf("-------------------Run %d: %.3f------------------\n",(int)testrep, runtime);

				if(!testrep) minruntime = runtime;
				minruntime = (minruntime > runtime) ? runtime : minruntime;

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

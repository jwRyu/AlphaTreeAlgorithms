#pragma once

#include "allocator.h"

#define TRIE_DEBUG 0

#if TRIE_DEBUG
#include <iostream>
using namespace std;
#endif


template<class Imgidx, class Trieidx>
class Trie
{
	Imgidx minidx;
	Trieidx **trie;
	//Imgidx *levelsize;
	Imgidx triesize, mask_field, nbits;
	int8 shamt, numlevels;
	//delayed non-leaf node push

#if TRIE_DEBUG
	Imgidx cursize;
#endif


	// 	//tmptmptmp
	// 	ofstream f;


		//int32 curSize;//tmp
public:
	Trie(Imgidx triesize_in)
	{
#if TRIE_DEBUG
		cursize = 0;
#endif
		//curSize = 0;//tmp
		triesize = triesize_in;
		Imgidx size;
		shamt = 2;
		for (int8 nbyte = sizeof(Trieidx); nbyte; nbyte >>= 1)
			shamt++;
		nbits = 1 << shamt;
		mask_field = nbits - 1;
		numlevels = 1;
		for (size = (triesize + 1) >> shamt; size; size >>= shamt)
			numlevels++;

		trie = (Trieidx**)Malloc(sizeof(Trieidx*) * numlevels);
		//levelsize = (Imgidx*)Malloc(sizeof(Imgidx) * numlevels);
		size = triesize + 1;
		for (int16 i = 0; i < numlevels; i++)
		{
			size >>= shamt;
			trie[i] = (Trieidx*)Malloc(sizeof(Trieidx) * (size + 1));
			for (Imgidx j = 0; j < (size + 1); j++)
				trie[i][j] = 0;
			//levelsize[i] = lvlsz;
		}
		minidx = triesize;
		//push(triesize);
		//curSize = 1;
		//curSize = 0;
		//tmp!
// 		f.open("C:/Users/jwryu/Google Drive/RUG/2019/AlphaTree_Trie/trie0rrr.dat",std::ofstream::out);
// 		f << triesize << endl;
	}
	~Trie()
	{
		for (int16 i = 0; i < numlevels; i++)
			Free(trie[i]);
		Free(trie);

		// 		f.close();//tmptmp
	}

	inline Imgidx top() { return minidx; }
	inline Imgidx get_minlev() { return minidx >> 1; }
	inline Imgidx min_incidence() { return minidx & 1; }
	inline int8 push(Imgidx in, int8 incidence)
	{
		//curSize++; //tmp
		//tmp
/*		f << '0' << '\n' << in << endl;*/
		return push((in << 1) + incidence);
	}
	inline int8 push(Imgidx in)
	{
		Imgidx n = in, s_in, shamt1;
		Trieidx *p;
		int8 ret = 0;

#if TRIE_DEBUG
		cout << "Push " << in << endl;
#endif

		//		curSize++; //tmp
				//tmp
		/*		f << '0' << '\n' << in << endl;*/

				//n = (in << 1) + incidence;
		s_in = n >> shamt;

		if (n < minidx)
		{
			minidx = n;
			ret = 1;
		}
		p = &(trie[0][s_in]);
		*p = *p | ((Trieidx)1 << (n & mask_field));
		n = s_in;
		s_in >>= shamt;
		for (int16 i = 1; i < numlevels; i++)
		{
			p = &(trie[i][s_in]);
			shamt1 = n & mask_field;
			//if (((*p) >> shamt1) & 1)
				//break;
			*p = *p | ((Trieidx)1 << shamt1);
			n = s_in;
			s_in >>= shamt;
		}
#if TRIE_DEBUG
		cursize++;
		cout << "Push returns " << (int)ret << endl;
#endif

		return ret;
	}
	inline void pop()
	{
		Imgidx s_idx = minidx >> shamt, shamt1;
		Trieidx *p, tmp;
		int16 lvl;

		//curSize--;//tmp
// 		//tmp
// 		f << '1' << '\n' << (minidx>>1) << endl;


#if TRIE_DEBUG
		cout << "Pop " << minidx << " - Cursize: " << cursize << endl;
#endif

		shamt1 = minidx & mask_field;
		p = &(trie[0][s_idx]);
		*p = *p & (~((Trieidx)1 << shamt1++));
		for (lvl = 0; !*p;)
		{
			minidx = s_idx;
			s_idx >>= shamt;
			shamt1 = minidx & mask_field;
			p = &(trie[++lvl][s_idx]);
			*p = *p & (~((Trieidx)1 << shamt1));
		}
		for (tmp = *p >> shamt1; !(tmp & 1); tmp >>= 1)
			shamt1++;
		minidx = (minidx & ~mask_field) | shamt1;
		while (lvl)
		{
			tmp = trie[--lvl][minidx];
			p = &tmp;
			s_idx = minidx;
			minidx <<= shamt;
			for (shamt1 = 0; !(tmp & 1); shamt1++)
				tmp >>= 1;
			minidx |= shamt1;
		}
#if TRIE_DEBUG
		cursize--;
		cout << "Pop: Next min: " << minidx << endl;
#endif
	}
};

Alpha-Tree algorithms

by Jiwoo Ryu (You), November 2021
Parallel radix sorting written by Paul Teeninga.

*The test code in the current version only accept pgm files as an input.

This C++ code implements several serial and shared-memory parallel alpha-tree algorithms for both low dynamic range (LDR) images and high dynamic range (HDR) images. Use the makefile to compile (the code was developed on gcc 9.3.0, Ubuntu 20.04.2 LTS). An alpha-tree can be built by calling the following method:

void BuildAlphaTree(Pixel *img, int height, int width, int channel, int connectivity, int algorithm, int numthreads, int tse)

- img: The input image (in either 1-ch or rgb)
- height, width, channel: Sizes of the input image
- connectivity: Pixel connectivity. Can be either 4 or 8. Note that some in algorithms including all parallel algorithms, 8-connectivity hasn't been implemented.
- algorithm: the index number of alpha-tree algorithm. See below for a list of algorithms with their indices.
- numthreads: the number of threads used in parallel algorithms. Ignored if the provided algorithm index corresponds to a serial algorithm.
- tse: a flag to indicate whether the tree size estimation (TSE) should be used. For some algorithms that TSE is not applicable, (flooding algorithms using trie) this flag is ignored.

List of algorithms

0. Unionfind: Union-find based alpha-tree algorithm similar to the Ouzounis's algorithm [1].
1. Flood using Hierarchical queue: An alpha-tree version of the flooding algorithm by [2].
2. Flood using Hierarchical queue with cache: Alpha-tree version of the flooding algorithm by [2], where a cache was used in its priority queue [3].
3. Flood using trie: Alpha-tree version of the flooding algorithm using trie by [4].
4. Flood using trie: [4] using the prioirty queue cache [3].
5. Flood using heap queue: Alpha-tree version of the flooding algorithm using heap queue by [5].
6. Flood using quad heap queue with cache: [5] using the prioirty queue cache [3].
7. Flooding using Hierarchical queue, Hypergraph (experimental): Flooding algorithm using hierarchical queue, on the hypergraph [6].
8. Flood using Trie and hypergraph (experimental): Flooding algorithm using trie queue, on the hypergraph [6].
9. Flood using hierarchical heap queue with cache (proposed): Flooding algorithm using hierarchical heap queue with cache [3].
11. Block-based parallel algorithm for LDR images (proposed): Parallel flooding algorithm for LDR image in [6].
12. Pilot-tree hybrid algorithm for HDR images (proposed): Alpha-tree version of hybrid max-tree algorithm in [7], for HDR images [6].


You don't have to use all of those algorithms. The algorithm indices of the best choices for different situations are as follows:

                                            Dynamic range
                                 low (<= 16-bit)      high (> 16-bit)
--------------------------------------------------------------------
	     single-thread    |          2                9
shared-memory multi-thread    |         11               12


After calling the method BuildAlphaTree, the alpha tree is constructed as the array of AlphaNodes (AlphaNode *node) in the class ATree. You can validate the alpha-tree by calling a simple area-based filter:

void AreaFilter(double *outimg, double area)

where outimg is the preallocated memory space to put the output image in, and area is the parameter of the area filter used as the threshold in filtering.

List of sentinel-2 remote sensing images used in [6]

This is the list of images used to test alpha-tree algorithms in my research [6]. You can download them from Copernicus Open Access Hub (https://scihub.copernicus.eu/dhus/#/home). You will need to make a free account, and probably have to wait for a few hours after requesting for each image to download them.

Here is the list of product names of the images I used that you can use to search in the open access hub.

S2A_MSIL1C_20151206T075212_N0204_R092_T37PFP_20151206T075547.SAFE
S2A_MSIL1C_20151206T075212_N0204_R092_T37PGQ_20151206T075547.SAFE
S2A_MSIL1C_20151229T075332_N0201_R135_T37PCN_20151229T080343.SAFE
S2A_MSIL1C_20151229T075332_N0201_R135_T37PDQ_20151229T080343.SAFE
S2A_MSIL1C_20151208T100412_N0204_R122_T32SNE_20151208T100409.SAFE
S2A_MSIL1C_20151208T100412_N0204_R122_T32SPE_20151208T100409.SAFE
S2A_MSIL1C_20151211T083342_N0204_R021_T36RUU_20151211T083340.SAFE
S2A_MSIL1C_20151211T083342_N0204_R021_T36RVV_20151211T083340.SAFE
S2A_MSIL1C_20151211T101412_N0204_R022_T32SLE_20151211T101905.SAFE

[1] G. K. Ouzounis and P. Soille, The alpha-tree algorithm, theory, algorithms, and applications, ser. JRC Technical Reports, Joint Research Centre. European Commission, 2012
[2] P. Salembier and J. Serra, “Flat zones filtering, connected operators, and filters by reconstruction,” IEEE Transactions on Image Processing, vol. 7, no. 4, pp. 1153–1160, 1995.
[3] J.  You,  S.  C.  Trager,  and  M.  H.  F.  Wilkinson,  “Hierarchical Heap Priority Queue for Alpha-Tree Algorithm”, to be submitted, 2021
[4] P. Teeninga and M. H. F. Wilkinson, “Fast and memory efficient sequentialmax-tree construction using a trie-based priority queue,”Submitted toPattern Recognition Letters, 2020.
[5] M. H. F. Wilkinson, “A fast component-tree algorithm for high dynamic-range images and second generation connectivity,” inICIP 2011, 18thIEEE International Conference on Image Processing, 2011, pp. 1021–1024.
[6] J.  You,  S.  C.  Trager,  and  M.  H.  F.  Wilkinson,  “Efficient Alpha-Tree Algorithms”, to be submitted, 2021
[7] U. Moschini, A. Meijster, and M. H. F. Wilkinson, “A hybrid shared-memory parallel max-tree algorithm for extreme dynamic-range images,”IEEE transactions on pattern analysis and machine intelligence, vol. 40,no. 3, pp. 513–526, 2017

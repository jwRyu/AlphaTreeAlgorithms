#include <AlphaTree.hpp>
#include <HHPQ.hpp>
#include <cmath>

#define RUN_FANCY_VIS_FOR_DEBUG 1

template <class Pixel> AlphaTree<Pixel>::~AlphaTree() { clear(); }

template <class Pixel>
AlphaNode<Pixel>::AlphaNode(Pixel pixelVal, double alpha_, ImgIdx parentidx_)
    : area(1), alpha(alpha_), sumPix((double)pixelVal), minPix(pixelVal), maxPix(pixelVal), parentIdx(parentidx_),
      _rootIdx(ROOTIDX) {
    // printf("AlphaNode<Pixel>::AlphaNode pixelVal = %d\n", (int)pixelVal);
}

template <class Pixel>
AlphaNode<Pixel>::AlphaNode(double alpha_, ImgIdx parentidx_)
    : area(0), alpha(alpha_), sumPix(0.0), minPix(0.0), maxPix(0.0), parentIdx(parentidx_), _rootIdx(ROOTIDX) {}

template <class Pixel> void AlphaTree<Pixel>::clear() {
    if (_node)
        Free(_node);
    _node = nullptr;
    if (_parentAry)
        Free(_parentAry);
    _parentAry = nullptr;
}

template <class Pixel>
void AlphaNode<Pixel>::set(ImgIdx area_in, double level, double sumPix_in, Pixel minPix_in, Pixel maxPix_in) {
    this->area = area_in;
    this->alpha = level;
    this->sumPix = sumPix_in;
    this->minPix = minPix_in;
    this->maxPix = maxPix_in;
}

template <class Pixel> void AlphaNode<Pixel>::add(AlphaNode *q) {
    this->area += q->area;
    this->sumPix += (double)q->sumPix;
    this->minPix = _min(this->minPix, q->minPix);
    this->maxPix = _max(this->maxPix, q->maxPix);
}

template <class Pixel> void AlphaNode<Pixel>::add(const AlphaNode &q) {
    this->area += q.area;
    this->sumPix += (double)q.sumPix;
    this->minPix = _min(this->minPix, q.minPix);
    this->maxPix = _max(this->maxPix, q.maxPix);
}

template <class Pixel> void AlphaNode<Pixel>::add(const Pixel &pix_val) {
    this->area++;
    this->sumPix += (double)pix_val;
    this->minPix = _min(this->minPix, pix_val);
    this->maxPix = _max(this->maxPix, pix_val);
}

template <class Pixel> void AlphaNode<Pixel>::copy(AlphaNode *q) {
    this->area = q->area;
    this->sumPix = q->sumPix;
    this->minPix = q->minPix;
    this->maxPix = q->maxPix;
}

template <class Pixel> void AlphaNode<Pixel>::connect_to_parent(AlphaNode *pPar, ImgIdx iPar) {
    this->parentIdx = iPar;
    pPar->add(this);
}

template <class Pixel> void AlphaNode<Pixel>::print(AlphaNode *_node) {
    double val = (double)this->sumPix;
    printf("Node idx: %d  alpha: %f area: %d, sumpix: %.0f, min-max: %d-%d  _rootIdx: %d  parent: %d\n",
           (int)(this - _node), (double)this->alpha, (int)this->area, (double)val, (int)this->minPix, (int)this->maxPix,
           (int)this->_rootIdx, (int)this->parentIdx);
}

template <class Pixel> void AlphaNode<Pixel>::print(AlphaNode *_node, int heading) {
    printf("%d: Node idx: %d\talpha: %f\tparent: %d\n               area: %d, sumpix: %f\n               min-max: "
           "%d-%d    _rootIdx: %d\n",
           heading, (int)(this - _node), (double)this->alpha, (int)this->parentIdx, (int)this->area,
           (double)this->sumPix, (int)this->minPix, (int)this->maxPix, (int)this->_rootIdx);
}

template <class Pixel>
void AlphaTree<Pixel>::BuildAlphaTree(const Pixel *img, int height_in, int width_in, int channel_in,
                                      std::string dMetric, int connectivity_in, int algorithm, int numthreads, int tse,
                                      double fparam1, double fparam2, int iparam1) {
    this->_height = (ImgIdx)height_in;
    this->_width = (ImgIdx)width_in;
    this->_channel = (ImgIdx)channel_in;
    this->_connectivity = (ImgIdx)connectivity_in;
    _curSize = 0;

    if (_connectivity != 4 && _connectivity != 8) {
        std::cout << "_connectivity should be 4 or 8\n" << std::endl;
        return;
    }

    if (dMetric == "L1")
        _pixelDissim = PixelDissimilarity<Pixel>(img, _height * _width, _channel, &PixelDissimilarity<Pixel>::L1);
    else if (dMetric == "L2")
        _pixelDissim = PixelDissimilarity<Pixel>(img, _height * _width, _channel, &PixelDissimilarity<Pixel>::L2);
    else if (dMetric == "LInfinity")
        _pixelDissim =
            PixelDissimilarity<Pixel>(img, _height * _width, _channel, &PixelDissimilarity<Pixel>::LInfinity);
    else
        _pixelDissim = PixelDissimilarity<Pixel>(img, _height * _width, _channel, &PixelDissimilarity<Pixel>::L2);

    // TEMP
    // EdgeSortingTest(img, numthreads);
    // TEMP
    // printf("BuildAlphaTree(): Temporary EdgeSortingTest Finished\n");
    // return;

    // switch
    if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("UnionFind"))
        Unionfind(img);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodHierQueueNoCache"))
        FloodHierarQueueNoCache(img, tse);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodHierQueue"))
        FloodHierarQueue(img);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodTrieQueueNoCache"))
        FloodTrieNoCache(img);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodTrieQueue"))
        FloodTrie(img);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodHeapQueueNoCache"))
        FloodHeapQueueNoCache(img);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodHeapQueueNaiveNoCache"))
        FloodHeapQueueNaiveNoCache(img);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodHeapQueue"))
        FloodHeapQueue(img);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodHierQueueHypergraph"))
        FloodHierarQueueHypergraph(img);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodTrieQueueHypergraph"))
        FloodTrieHypergraph(img);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodHierHeapQueueNoCache"))
        FloodHierarHeapQueueNoCache(img, fparam1, fparam2, iparam1);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodHierHeapQueue"))
        FloodHierarHeapQueue(img, fparam1, fparam2, iparam1);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodHierHeapQueueHisteq"))
        FloodHierHeapQueueHisteq(img);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodLadderQueue"))
        FloodLadderQueue(img, iparam1);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("HybridParallel"))
        HybridParallel(img, numthreads);
    else if (algorithm == alphatreeConfig.getAlphaTreeAlgorithmCode("FloodHierQueueParallel"))
        FloodHierQueueParallel(img, numthreads);
    else {
        std::cerr << "[AlphaTree::BuildAlphaTree]"
                  << " Unknown algorithm code (" << algorithm << std::endl;
        return;
    }
}

template <class Pixel> void AlphaTree<Pixel>::AlphaFilter(double *outimg, double alpha) {
    ImgIdx i, imgSize;
    AlphaNode<Pixel> *pNode;

    alpha = _node[_rootIdx].alpha * alpha;

    imgSize = _height * _width;
    // val = 1;
    for (i = 0; i < imgSize; i++) {
        pNode = _parentAry ? &_node[_parentAry[i]] : &_node[i];
        while (pNode->parentIdx != -1 && pNode->alpha < alpha)
            pNode = &_node[pNode->parentIdx];
        outimg[i] = (double)pNode->area;
    }
}

template <class Pixel> void AlphaTree<Pixel>::AreaFilter(double *outimg, double area) {
    ImgIdx i, imgSize;
    ImgIdx iarea;
    AlphaNode<Pixel> *pNode;

    imgSize = _height * _width;
    iarea = (ImgIdx)(area * (double)imgSize);
    iarea = _min(imgSize, _max(0, iarea));
    // val = 1;
    for (i = 0; i < imgSize; i++) {
        pNode = _parentAry ? &_node[_parentAry[i]] : &_node[i];
        while (pNode->parentIdx != -1 && pNode->area < iarea)
            pNode = &_node[pNode->parentIdx];
        outimg[i] = (double)pNode->alpha;
    }
}

template <class Pixel> void AlphaTree<Pixel>::printTree() const {
    for (int i = 0; i < _maxSize; i++)
        if (_node[i].area > 0)
            _node[i].print(_node);
}

template <class Pixel>
void AlphaTree<Pixel>::printGraph(const uint8_t *isVisited, const uint8_t *edge, const Pixel *img) const {
    printf("Index |   Graph (Left)   |   Image (Right)\n");
    printf("------------------------------------------\n");
    for (int i = 0; i < _height - 1; i++) {
        printf("%3d | ", i * _width);
        for (int j = 0; j < _width - 1; j++) {
            int imgIdx = i * _width + j;
            int dimgIdx = imgIdx * 2 + 1;
            if (edge == 0 || edge[dimgIdx] == QItem::EDGE_STANDBY)
                printf("%d   ", (int)isVisited[imgIdx]);
            else if (edge[dimgIdx] == QItem::EDGE_ENQUEUED)
                printf("%d . ", (int)isVisited[imgIdx]);
            else if (edge[dimgIdx] == QItem::EDGE_DEQUEUED)
                printf("%d ^ ", (int)isVisited[imgIdx]);
            else if (edge[dimgIdx] == QItem::EDGE_CONNECTED)
                printf("%d - ", (int)isVisited[imgIdx]);
            else if (edge[dimgIdx] == QItem::EDGE_REDUNDANT)
                printf("%d x ", (int)isVisited[imgIdx]);
            else if (edge[dimgIdx] == QItem::EDGE_ESSENTIAL)
                printf("%d * ", (int)isVisited[imgIdx]);
            else
                printf("%d ? ", (int)isVisited[imgIdx]);
        }
        printf("%d   |   ", (int)isVisited[i * _width + _width - 1]);

        if (_channel == 1) {
            for (int j = 0; j < _width; j++) {
                int imgIdx = i * _width + j;
                printf("%2d ", (int)img[imgIdx]);
            }
        } else {
            const auto imgSize = _height * _width;
            for (int j = 0; j < _width; j++) {
                int imgIdx = i * _width + j;
                for (int ch = 0; ch < _channel - 1; ch++)
                    printf("%2d/", (int)img[imgIdx + ch * imgSize]);
                printf("%2d ", (int)img[imgIdx + (_channel - 1) * imgSize]);
            }
        }
        printf("\n    | ");

        for (int j = 0; j < _width; j++) {
            int imgIdx = i * _width + j;
            int dimgIdx = imgIdx * 2;
            if (edge == 0 || edge[dimgIdx] == QItem::EDGE_STANDBY)
                printf("    ");
            else if (edge[dimgIdx] == QItem::EDGE_ENQUEUED)
                printf(".   ");
            else if (edge[dimgIdx] == QItem::EDGE_DEQUEUED)
                printf("^   ");
            else if (edge[dimgIdx] == QItem::EDGE_CONNECTED)
                printf("|   ");
            else if (edge[dimgIdx] == QItem::EDGE_REDUNDANT)
                printf("x   ");
            else if (edge[dimgIdx] == QItem::EDGE_ESSENTIAL)
                printf("*   ");
            else
                printf("?   ");
        }
        printf("|\n");
    }
    printf("%3d | ", (_height - 1) * _width);
    for (int j = 0; j < _width - 1; j++) {
        int imgIdx = (_height - 1) * _width + j;
        int dimgIdx = imgIdx * 2 + 1;
        if (edge == 0 || edge[dimgIdx] == QItem::EDGE_STANDBY)
            printf("%d   ", (int)isVisited[imgIdx]);
        else if (edge[dimgIdx] == QItem::EDGE_ENQUEUED)
            printf("%d . ", (int)isVisited[imgIdx]);
        else if (edge[dimgIdx] == QItem::EDGE_DEQUEUED)
            printf("%d ^ ", (int)isVisited[imgIdx]);
        else if (edge[dimgIdx] == QItem::EDGE_CONNECTED)
            printf("%d - ", (int)isVisited[imgIdx]);
        else if (edge[dimgIdx] == QItem::EDGE_REDUNDANT)
            printf("%d x ", (int)isVisited[imgIdx]);
        else if (edge[dimgIdx] == QItem::EDGE_ESSENTIAL)
            printf("%d * ", (int)isVisited[imgIdx]);
        else
            printf("%d ? ", (int)isVisited[imgIdx]);
    }
    printf("%d   |   ", (int)isVisited[_width * _height - 1]);

    if (_channel == 1) {
        for (int j = 0; j < _width; j++) {
            int imgIdx = (_height - 1) * _width + j;
            printf("%2d ", (int)img[imgIdx]);
        }
    } else {
        const auto imgSize = _height * _width;
        for (int j = 0; j < _width; j++) {
            int imgIdx = (_height - 1) * _width + j;
            for (int ch = 0; ch < _channel - 1; ch++)
                printf("%2d/", (int)img[imgIdx + ch * imgSize]);
            printf("%2d ", (int)img[imgIdx + (_channel - 1) * imgSize]);
        }
    }
    printf("\n");
}

template <class Pixel>
void AlphaTree<Pixel>::printGraph(const uint8_t *isVisited, const bool *isRedundant, const Pixel *img) const {
    printf("Index |   Graph (Left)   |   Image (Right)\n");
    printf("------------------------------------------\n");
    for (int i = 0; i < _height - 1; i++) {
        printf("%3d | ", i * _width);
        for (int j = 0; j < _width - 1; j++) {
            int imgIdx = i * _width + j;
            int dimgIdx = imgIdx * 2 + 1;
            if (isRedundant == 0 || isRedundant[dimgIdx] == 0)
                printf("%d   ", (int)isVisited[imgIdx]);
            else
                printf("%d x ", (int)isVisited[imgIdx]);
        }
        printf("%d   |   ", (int)isVisited[i * _width + _width - 1]);

        for (int j = 0; j < _width; j++) {
            int imgIdx = i * _width + j;
            printf("%2d ", (int)img[imgIdx]);
        }
        printf("\n    | ");

        for (int j = 0; j < _width; j++) {
            int imgIdx = i * _width + j;
            int dimgIdx = imgIdx * 2;
            if (isRedundant == 0 || isRedundant[dimgIdx] == 0)
                printf("    ");
            else
                printf("x   ");
        }
        printf("|\n");
    }
    printf("%3d | ", (_height - 1) * _width);
    for (int j = 0; j < _width - 1; j++) {
        int imgIdx = (_height - 1) * _width + j;
        int dimgIdx = imgIdx * 2 + 1;
        if (isRedundant == 0 || isRedundant[dimgIdx] == 0)
            printf("%d   ", (int)isVisited[imgIdx]);
        else
            printf("%d x ", (int)isVisited[imgIdx]);
    }
    printf("%d   |   ", (int)isVisited[_width * _height - 1]);

    for (int j = 0; j < _width; j++) {
        int imgIdx = (_height - 1) * _width + j;
        printf("%2d ", (int)img[imgIdx]);
    }
    printf("\n");
}

template <class Pixel> void AlphaTree<Pixel>::printParentAry() const {
    printf("Parent Array\n");
    for (int i = 0; i < _height; i++) {
        for (int j = 0; j < _width; j++) {
            int imgIdx = i * _width + j;
            if (_parentAry[imgIdx] >= 0)
                printf("%3d ", _parentAry[imgIdx]);
            else
                printf("  . ");
        }
        printf("\n");
    }
}

// template <class Pixel>
// void AlphaTree<Pixel>::printAll(const uint8_t *isVisited, const bool *isRedundant, const Pixel *img) const {
//     // printTree();
//     printParentAry();
//     printGraph(isVisited, isRedundant, img);
// }

template <class Pixel>
void AlphaTree<Pixel>::printAll(const uint8_t *isVisited, const uint8_t *edge, const Pixel *img) const {
    // printTree();
    printParentAry();
    printGraph(isVisited, edge, img);
}

template <class Pixel> Pixel AlphaTree<Pixel>::abs_diff(Pixel p, Pixel q) {
    if (p > q)
        return p - q;
    else
        return q - p;
}

template <class Pixel> uint8_t AlphaTree<Pixel>::compute_incidedge_queue(Pixel d0, Pixel d1) {
    if (d0 <= d1)
        return 0x4;
    else
        return 0x1;
}

template <class Pixel>
void AlphaTree<Pixel>::computeDimgRadixf(RankItem<double> *&rankitem, const Pixel *img, SortValue<float> *&vals) {
    if (_connectivity == 4) {
#pragma omp parallel for schedule(guided, 1)
        for (ImgIdx i = 0; i < _height; i++) {
            double d;
            ImgIdx imgidx = i * _width;
            ImgIdx dimgidx = imgidx << (_connectivity >> 2);
            ImgIdx contidx = dimgidx - i;
            if (i < _height - 1) {
                for (ImgIdx j = 0; j < _width - 1; j++) {
                    // caclulate histogram here
                    vals[contidx].val_ = (float)_pixelDissim.computeDissimilarity(imgidx, imgidx + _width);
                    d = (double)(vals[contidx].val_);
                    rankitem[contidx].alpha = d;
                    rankitem[contidx++].dimgidx = dimgidx++;

                    vals[contidx].val_ = (float)_pixelDissim.computeDissimilarity(imgidx, imgidx + 1);
                    d = (double)(vals[contidx].val_);
                    rankitem[contidx].alpha = d;
                    rankitem[contidx++].dimgidx = dimgidx++;
                    imgidx++;
                }
                vals[contidx].val_ = (float)_pixelDissim.computeDissimilarity(imgidx, imgidx + _width);
                d = (double)(vals[contidx].val_);
                rankitem[contidx].alpha = d;
                rankitem[contidx++].dimgidx = dimgidx;
                dimgidx += 2;
                imgidx++;
            } else {
                for (ImgIdx j = 0; j < _width - 1; j++) {
                    dimgidx++;
                    vals[contidx].val_ = (float)_pixelDissim.computeDissimilarity(imgidx, imgidx + 1);
                    d = (double)(vals[contidx].val_);
                    rankitem[contidx].alpha = d;
                    rankitem[contidx++].dimgidx = dimgidx++;
                    imgidx++;
                }
            }
        }
    }
}

template <class Pixel>
void AlphaTree<Pixel>::compute_dimg_par4(RankItem<double> *&rankitem, const Pixel *img, SortValue<Pixel> *&vals) {
    ImgIdx contidx, dimgidx, imgidx, i, j;

    contidx = imgidx = dimgidx = 0;
    if (_connectivity == 4) {
        if (_channel == 1) {
#pragma omp parallel for schedule(guided, 1) private(imgidx, dimgidx, contidx, i, j)
            for (i = 0; i < _height; i++) {
                double d;
                imgidx = i * _width;
                dimgidx = imgidx << (_connectivity >> 2);
                contidx = dimgidx - i;
                if (i < _height - 1) {
                    for (j = 0; j < _width - 1; j++) {
                        // caclulate histogram here
                        vals[contidx].val_ = abs_diff(img[imgidx + _width], img[imgidx]);
                        d = (double)(vals[contidx].val_);
                        rankitem[contidx].alpha = d;
                        rankitem[contidx++].dimgidx = dimgidx++;

                        vals[contidx].val_ = abs_diff(img[imgidx + 1], img[imgidx]);
                        d = (double)(vals[contidx].val_);
                        rankitem[contidx].alpha = d;
                        rankitem[contidx++].dimgidx = dimgidx++;
                        imgidx++;
                    }
                    vals[contidx].val_ = abs_diff(img[imgidx + _width], img[imgidx]);
                    d = (double)(vals[contidx].val_);
                    rankitem[contidx].alpha = d;
                    rankitem[contidx++].dimgidx = dimgidx;
                    dimgidx += 2;
                    imgidx++;
                } else {
                    for (j = 0; j < _width - 1; j++) {
                        dimgidx++;
                        vals[contidx].val_ = abs_diff(img[imgidx + 1], img[imgidx]);
                        d = (double)(vals[contidx].val_);
                        rankitem[contidx].alpha = d;
                        rankitem[contidx++].dimgidx = dimgidx++;
                        imgidx++;
                    }
                }
            }
        } else {
            ImgIdx chstride = _channel, chstride2 = _channel * 2;
            ImgIdx imgSize = _height * _width;
            ImgIdx dimgSize = _height * _width * (_connectivity >> 1);
            double *dimg_3ch = (double *)Calloc(_channel * dimgSize * sizeof(double));
            ImgIdx linestride = _width * (_connectivity >> 1) - (_connectivity >> 2);

#pragma omp parallel for schedule(guided, 1) private(imgidx, dimgidx, i, j)
            for (int hch = 0; hch < _channel * _height; hch++) {
                double d;
                int ch = hch / _height;
                i = hch % _height;
                const Pixel *pimg = img + ch * imgSize + i * _width;
                double *pdimg = dimg_3ch + ch + i * _width * (_connectivity >> 1) * _channel;
                imgidx = dimgidx = 0;
                if (i < _height - 1) {
                    for (j = 0; j < _width - 1; j++) {
                        d = (double)pimg[imgidx + _width] - (double)pimg[imgidx];
                        pdimg[dimgidx] = d * d;
                        dimgidx += chstride;
                        d = (double)pimg[imgidx + 1] - (double)pimg[imgidx];
                        pdimg[dimgidx] = d * d;
                        dimgidx += chstride;
                        imgidx++;
                    }
                    d = (double)pimg[imgidx + _width] - (double)pimg[imgidx];
                    pdimg[dimgidx] = d * d;
                    dimgidx += chstride2;
                    imgidx++;
                } else {
                    for (j = 0; j < _width - 1; j++) {
                        dimgidx += chstride;
                        d = (double)pimg[imgidx + 1] - (double)pimg[imgidx];
                        pdimg[dimgidx] = d * d;
                        dimgidx += chstride;
                        imgidx++;
                    }
                }
            }

#pragma omp parallel for schedule(guided, 1) private(imgidx, dimgidx, contidx, i, j)
            for (i = 0; i < _height; i++) {
                double d;
                double *pdimg = dimg_3ch + i * _width * (_connectivity >> 1) * _channel;
                contidx = i * linestride;
                imgidx = i * _width;
                ImgIdx dimgidx_3 = 0;
                dimgidx = i * _width * (_connectivity >> 1);
                if (i < _height - 1) {
                    for (j = 0; j < _width - 1; j++) {
                        d = pdimg[dimgidx_3] + pdimg[dimgidx_3 + 1] + pdimg[dimgidx_3 + 2]; // only for 3-ch
                        rankitem[contidx].alpha = vals[contidx].val_ = sqrt(d / (double)_channel);
                        rankitem[contidx].dimgidx = dimgidx;
                        contidx++;

                        dimgidx_3 += chstride;
                        dimgidx++;

                        d = pdimg[dimgidx_3] + pdimg[dimgidx_3 + 1] + pdimg[dimgidx_3 + 2]; // only for 3-ch
                        rankitem[contidx].alpha = vals[contidx].val_ = sqrt(d / (double)_channel);
                        rankitem[contidx].dimgidx = dimgidx;
                        contidx++;

                        dimgidx_3 += chstride;
                        dimgidx++;

                        imgidx++;
                    }
                    d = pdimg[dimgidx_3] + pdimg[dimgidx_3 + 1] + pdimg[dimgidx_3 + 2]; // only for 3-ch
                    rankitem[contidx].alpha = vals[contidx].val_ = sqrt(d / (double)_channel);
                    rankitem[contidx].dimgidx = dimgidx;
                    contidx++;
                    dimgidx++;

                    dimgidx_3 += chstride2;
                    imgidx++;
                } else {
                    for (j = 0; j < _width - 1; j++) {
                        dimgidx++;
                        dimgidx_3 += chstride;

                        d = pdimg[dimgidx_3] + pdimg[dimgidx_3 + 1] + pdimg[dimgidx_3 + 2]; // only for 3-ch
                        rankitem[contidx].alpha = vals[contidx].val_ = sqrt(d / (double)_channel);
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
    } else if (_connectivity == 8) // not implemented yet
    {
        if (_channel == 1) {
#pragma omp parallel for schedule(guided, 1) private(imgidx, dimgidx, contidx, i, j)
            for (i = 0; i < _height; i++) {
                double d;
                imgidx = i * _width;
                dimgidx = imgidx << (_connectivity >> 2);
                contidx = dimgidx - ((_connectivity == 4) ? (i) : (3 * i));
                if (_connectivity == 8 && i > 0)
                    contidx -= _width - 1;

                if (i < _height - 1) {
                    for (j = 0; j < _width - 1; j++) {
                        // caclulate histogram here
                        vals[contidx].val_ = abs_diff(img[imgidx + _width], img[imgidx]);
                        d = (double)(vals[contidx].val_);
                        rankitem[contidx].alpha = d;
                        rankitem[contidx++].dimgidx = dimgidx++;

                        vals[contidx].val_ = abs_diff(img[imgidx + _width + 1], img[imgidx]);
                        d = (double)(vals[contidx].val_);
                        rankitem[contidx].alpha = d;
                        rankitem[contidx++].dimgidx = dimgidx++;

                        vals[contidx].val_ = abs_diff(img[imgidx + 1], img[imgidx]);
                        d = (double)(vals[contidx].val_);
                        rankitem[contidx].alpha = d;
                        rankitem[contidx++].dimgidx = dimgidx++;

                        if (i > 0) {
                            vals[contidx].val_ = abs_diff(img[imgidx - _width + 1], img[imgidx]);
                            d = (double)(vals[contidx].val_);
                            rankitem[contidx].alpha = d;
                            rankitem[contidx++].dimgidx = dimgidx;
                        }
                        dimgidx++;

                        imgidx++;
                    }

                    vals[contidx].val_ = abs_diff(img[imgidx + _width], img[imgidx]);
                    d = (double)(vals[contidx].val_);
                    rankitem[contidx].alpha = d;
                    rankitem[contidx++].dimgidx = dimgidx;
                    dimgidx += 4;
                    imgidx++;
                } else {
                    for (j = 0; j < _width - 1; j++) {
                        dimgidx += 2;
                        vals[contidx].val_ = abs_diff(img[imgidx + 1], img[imgidx]);
                        d = (double)(vals[contidx].val_);
                        rankitem[contidx].alpha = d;
                        rankitem[contidx++].dimgidx = dimgidx++;

                        vals[contidx].val_ = abs_diff(img[imgidx - _width + 1], img[imgidx]);
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

template <class Pixel> Pixel AlphaTree<Pixel>::compute_dimg(double *dimg, const Pixel *img) {
    using Value = double;
    ImgIdx dimgidx, imgidx, stride_w = _width, i, j;
    ImgIdx imgSize = _width * _height;

    const Pixel *pimg = img;
    Pixel dmax = 0;

    imgidx = dimgidx = 0;
    if (_connectivity == 4) {
        if (_channel == 1) {
            for (i = 0; i < _height - 1; i++) {
                for (j = 0; j < _width - 1; j++) {
                    dimg[dimgidx] = (Value)(abs_diff(img[imgidx + stride_w], img[imgidx]));
                    dmax = _max(dmax, dimg[dimgidx]);
                    dimgidx++;
                    dimg[dimgidx] = (Value)(abs_diff(img[imgidx + 1], img[imgidx]));
                    dmax = _max(dmax, dimg[dimgidx]);
                    dimgidx++;
                    imgidx++;
                }
                dimg[dimgidx] = (Value)(abs_diff(img[imgidx + stride_w], img[imgidx]));
                dmax = _max(dmax, dimg[dimgidx]);
                dimgidx += 2;
                imgidx++;
            }
            for (j = 0; j < _width - 1; j++) {
                dimgidx++;
                dimg[dimgidx] = (Value)(abs_diff(img[imgidx + 1], img[imgidx]));
                dmax = _max(dmax, dimg[dimgidx]);
                dimgidx++;
                imgidx++;
            }
        } else {
            double d;
            for (int ch = 0; ch < _channel; ch++) {
                imgidx = dimgidx = 0;
                for (i = 0; i < _height - 1; i++) {
                    for (j = 0; j < _width - 1; j++) {
                        d = (double)pimg[imgidx + stride_w] - (double)pimg[imgidx];
                        if (ch == 0)
                            dimg[dimgidx] = d * d;
                        else if (ch != _channel - 1)
                            dimg[dimgidx] += d * d;
                        else
                            dimg[dimgidx] = sqrt(dimg[dimgidx] + d * d);
                        dmax = _max(dmax, dimg[dimgidx]);
                        dimgidx++;

                        d = (double)pimg[imgidx + 1] - (double)pimg[imgidx];
                        if (ch == 0)
                            dimg[dimgidx] = d * d;
                        else if (ch != _channel - 1)
                            dimg[dimgidx] += d * d;
                        else
                            dimg[dimgidx] = sqrt(dimg[dimgidx] + d * d);
                        dmax = _max(dmax, dimg[dimgidx]);
                        dimgidx++;
                        imgidx++;
                    }
                    d = (double)pimg[imgidx + stride_w] - (double)pimg[imgidx];
                    if (ch == 0)
                        dimg[dimgidx] = d * d;
                    else if (ch != _channel - 1)
                        dimg[dimgidx] += d * d;
                    else
                        dimg[dimgidx] = sqrt(dimg[dimgidx] + d * d);
                    dmax = _max(dmax, dimg[dimgidx]);
                    dimgidx += 2;
                    imgidx++;
                }
                for (j = 0; j < _width - 1; j++) {
                    dimgidx++;
                    d = (double)pimg[imgidx + 1] - (double)pimg[imgidx];
                    if (ch == 0)
                        dimg[dimgidx] = d * d;
                    else if (ch != _channel - 1)
                        dimg[dimgidx] += d * d;
                    else
                        dimg[dimgidx] = sqrt(dimg[dimgidx] + d * d);
                    dmax = _max(dmax, dimg[dimgidx]);
                    dimgidx++;
                    imgidx++;
                }
                pimg += imgSize;
            }
        }
    } else if (_connectivity == 8) {
        if (_channel == 1) {
            //   -  -  3
            //   -  p  2
            //   -  0  1
            // top,middle
            for (i = 0; i < _height - 1; i++) {
                for (j = 0; j < _width - 1; j++) {
                    dimg[dimgidx++] = (Value)(abs_diff((int64_t)img[imgidx + _width], (int64_t)img[imgidx])); // 0
                    dmax = _max(dmax, dimg[dimgidx - 1]);
                    dimg[dimgidx++] = (Value)(abs_diff((int64_t)img[imgidx + _width + 1], (int64_t)img[imgidx])); // 1
                    dmax = _max(dmax, dimg[dimgidx - 1]);
                    dimg[dimgidx++] = (Value)(abs_diff((int64_t)img[imgidx + 1], (int64_t)img[imgidx])); // 2
                    dmax = _max(dmax, dimg[dimgidx - 1]);
                    if (i > 0) {
                        dimg[dimgidx] = (Value)(abs_diff((int64_t)img[imgidx - _width + 1], (int64_t)img[imgidx])); // 3
                        dmax = _max(dmax, dimg[dimgidx]);
                    }
                    dimgidx++;
                    imgidx++;
                }
                dimg[dimgidx] = (Value)(abs_diff((int64_t)img[imgidx + _width], (int64_t)img[imgidx])); // 0
                dmax = _max(dmax, dimg[dimgidx]);
                dimgidx += 4; // skip 1,2,3
                imgidx++;
            }

            // bottom
            dimgidx += 2; // skip 0,1
            for (j = 0; j < _width - 1; j++) {
                dimg[dimgidx++] = (Value)(abs_diff((int64_t)img[imgidx + 1], (int64_t)img[imgidx])); // 2
                dmax = _max(dmax, dimg[dimgidx - 1]);
                dimg[dimgidx] = (Value)(abs_diff((int64_t)img[imgidx - _width + 1], (int64_t)img[imgidx])); // 3
                dmax = _max(dmax, dimg[dimgidx]);
                dimgidx += 3;
                imgidx++;
            }
        }
    }

    return dmax;
}

template <class Pixel> Pixel AlphaTree<Pixel>::compute_dimg1(Pixel *dimg, ImgIdx *dhist, const Pixel *img) {
    ImgIdx dimgidx, imgidx, stride_w = _width, i, j;
    Pixel maxdiff = 0;

    imgidx = dimgidx = 0;
    if (_connectivity == 4) {
        for (i = 0; i < _height - 1; i++) {
            for (j = 0; j < _width - 1; j++) {
                dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + stride_w] - (int64_t)img[imgidx]));
                maxdiff = _max(maxdiff, dimg[dimgidx]);
                dhist[dimg[dimgidx++]]++;

                dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + 1] - (int64_t)img[imgidx]));
                maxdiff = _max(maxdiff, dimg[dimgidx]);
                dhist[dimg[dimgidx++]]++;
                imgidx++;
            }
            dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + stride_w] - (int64_t)img[imgidx]));
            maxdiff = _max(maxdiff, dimg[dimgidx]);
            dhist[dimg[dimgidx++]]++;
            dimgidx++;
            imgidx++;
        }
        for (j = 0; j < _width - 1; j++) {
            dimgidx++;
            dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + 1] - (int64_t)img[imgidx]));
            maxdiff = _max(maxdiff, dimg[dimgidx]);
            dhist[dimg[dimgidx++]]++;
            imgidx++;
        }
    } else if (_connectivity == 8) {
        //   -  -  3
        //   -  p  2
        //   -  0  1
        // top,middle
        for (i = 0; i < _height - 1; i++) {
            for (j = 0; j < _width - 1; j++) {
                dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + _width] - (int64_t)img[imgidx])); // 0
                maxdiff = _max(maxdiff, dimg[dimgidx]);
                dhist[dimg[dimgidx++]]++;
                dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + _width + 1] - (int64_t)img[imgidx])); // 1
                maxdiff = _max(maxdiff, dimg[dimgidx]);
                dhist[dimg[dimgidx++]]++;
                dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + 1] - (int64_t)img[imgidx])); // 2
                maxdiff = _max(maxdiff, dimg[dimgidx]);
                dhist[dimg[dimgidx++]]++;
                if (i > 0) {
                    dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx - _width + 1] - (int64_t)img[imgidx])); // 3
                    maxdiff = _max(maxdiff, dimg[dimgidx]);
                    dhist[dimg[dimgidx]]++;
                }
                dimgidx++;
                imgidx++;
            }
            dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + _width] - (int64_t)img[imgidx])); // 0
            maxdiff = _max(maxdiff, dimg[dimgidx]);
            dhist[dimg[dimgidx]]++;
            dimgidx += 4; // skip 1,2,3
            imgidx++;
        }

        // bottom
        dimgidx += 2; // skip 0,1
        for (j = 0; j < _width - 1; j++) {
            dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + 1] - (int64_t)img[imgidx])); // 2
            maxdiff = _max(maxdiff, dimg[dimgidx]);
            dhist[dimg[dimgidx++]]++;
            dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx - _width + 1] - (int64_t)img[imgidx])); // 3
            maxdiff = _max(maxdiff, dimg[dimgidx]);
            dhist[dimg[dimgidx]]++;
            dimgidx += 3;
            imgidx++;
        }
    }

    return maxdiff;
}

template <class Pixel>
void AlphaTree<Pixel>::compute_dimg(ImgIdx &minidx, double &mindiff, Pixel *dimg, ImgIdx *dhist, const Pixel *img,
                                    double a) {
    ImgIdx dimgidx, imgidx, stride_w = _width, i, j;
    int hidx;
    mindiff = (double)((Pixel)(-1));
    imgidx = dimgidx = 0;

    if (_connectivity == 4) {
        for (i = 0; i < _height - 1; i++) {
            for (j = 0; j < _width - 1; j++) {
                dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + stride_w], img[imgidx]));
                if (dimg[dimgidx] < mindiff) {
                    mindiff = dimg[dimgidx];
                    minidx = i * _width + j;
                }
                hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
                dhist[hidx]++;
                dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));
                if (dimg[dimgidx] < mindiff) {
                    mindiff = dimg[dimgidx];
                    minidx = i * _width + j;
                }
                hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
                dhist[hidx]++;
                imgidx++;
            }
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + stride_w], img[imgidx]));
            if (dimg[dimgidx] < mindiff) {
                mindiff = dimg[dimgidx];
                minidx = i * _width + j;
            }
            hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
            dhist[hidx]++;
            dimgidx++;
            imgidx++;
        }
        for (j = 0; j < _width - 1; j++) {
            dimgidx++;
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));
            if (dimg[dimgidx] < mindiff) {
                mindiff = dimg[dimgidx];
                minidx = i * _width + j;
            }
            hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
            dhist[hidx]++;
            imgidx++;
        }
    } else if (_connectivity == 8) {
        //   -  -  3
        //   -  p  2
        //   -  0  1
        // top,middle
        for (i = 0; i < _height - 1; i++) {
            for (j = 0; j < _width - 1; j++) {
                dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + _width], img[imgidx])); // 0
                if (dimg[dimgidx] < mindiff) {
                    mindiff = dimg[dimgidx];
                    minidx = i * _width + j;
                }
                hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
                dhist[hidx]++;
                dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + _width + 1], img[imgidx])); // 1
                if (dimg[dimgidx] < mindiff) {
                    mindiff = dimg[dimgidx];
                    minidx = i * _width + j;
                }
                hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
                dhist[hidx]++;
                dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx])); // 2
                if (dimg[dimgidx] < mindiff) {
                    mindiff = dimg[dimgidx];
                    minidx = i * _width + j;
                }
                hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
                dhist[hidx]++;
                if (i > 0) {
                    dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx - _width + 1], img[imgidx])); // 3
                    if (dimg[dimgidx] < mindiff) {
                        mindiff = dimg[dimgidx];
                        minidx = i * _width + j;
                    }
                    hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
                    dhist[hidx]++;
                }
                dimgidx++;
                imgidx++;
            }
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + _width], img[imgidx])); // 0
            if (dimg[dimgidx] < mindiff) {
                mindiff = dimg[dimgidx];
                minidx = i * _width + j;
            }
            hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
            dhist[hidx]++;
            dimgidx += 4; // skip 1,2,3
            imgidx++;
        }

        // bottom
        dimgidx += 2; // skip 0,1
        for (j = 0; j < _width - 1; j++) {
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx])); // 2
            if (dimg[dimgidx] < mindiff) {
                mindiff = dimg[dimgidx];
                minidx = i * _width + j;
            }
            hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
            dhist[hidx]++;
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx - _width + 1], img[imgidx])); // 3
            if (dimg[dimgidx] < mindiff) {
                mindiff = dimg[dimgidx];
                minidx = i * _width + j;
            }
            hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
            dhist[hidx]++;
            dimgidx += 3;
            imgidx++;
        }
    }
}

template <class Pixel>
void AlphaTree<Pixel>::compute_dimg_hhpq(double *dimg, ImgIdx *dhist, const Pixel *img, double a) {
    ImgIdx imgidx = 0;
    ImgIdx dimgidx = 0;
    if (_connectivity == 4) {
        for (ImgIdx i = 0; i < _height - 1; i++) {
            for (ImgIdx j = 0; j < _width - 1; j++) {
                dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + _width);
                dhist[HHPQ::alphaToLevel(dimg[dimgidx++], a)]++;
                dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + 1);
                dhist[HHPQ::alphaToLevel(dimg[dimgidx++], a)]++;
                imgidx++;
            }
            dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + _width);
            dhist[HHPQ::alphaToLevel(dimg[dimgidx++], a)]++;
            dimgidx++;
            imgidx++;
        }
        for (ImgIdx j = 0; j < _width - 1; j++) {
            dimgidx++;
            dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + 1);
            dhist[HHPQ::alphaToLevel(dimg[dimgidx++], a)]++;
            imgidx++;
        }
    } else if (_connectivity == 8) {
        //   -  -  3
        //   -  p  2
        //   -  0  1
        // top,middle
        for (ImgIdx i = 0; i < _height - 1; i++) {
            for (ImgIdx j = 0; j < _width - 1; j++) {
                dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + _width);
                dhist[HHPQ::alphaToLevel(dimg[dimgidx++], a)]++;
                dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + _width + 1);
                dhist[HHPQ::alphaToLevel(dimg[dimgidx++], a)]++;
                dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + 1);
                dhist[HHPQ::alphaToLevel(dimg[dimgidx++], a)]++;
                if (i > 0) {
                    dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx - _width + 1);
                    dhist[HHPQ::alphaToLevel(dimg[dimgidx], a)]++;
                }
                dimgidx++;
                imgidx++;
            }
            dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + _width);
            dhist[HHPQ::alphaToLevel(dimg[dimgidx], a)]++;
            dimgidx += 4;
            imgidx++;
        }

        // bottom
        dimgidx += 2; // skip 0,1
        for (ImgIdx j = 0; j < _width - 1; j++) {
            dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + 1);
            dhist[HHPQ::alphaToLevel(dimg[dimgidx++], a)]++;
            dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx - _width + 1);
            dhist[HHPQ::alphaToLevel(dimg[dimgidx], a)]++;
            dimgidx += 3;
            imgidx++;
        }
    }
}

// void setEdgeStatus(double *dimg, ImgIdx *dhist, const Pixel *img, double a, uint8_t *edgeStatus) {
// {
//     const auto width2 = _width * 2;
//         for (ImgIdx i = 0; i < _height; i++) {
//             const bool top = i > 0;
//             const bool bottom = i < _height - 1;
//             for (ImgIdx j = 0; j < _width; j++) {
//                 const bool left = j > 0;
//                 const bool right = j < _width - 1;
//                 const ImgIdx imgidx = i * _width + j;
//                 const ImgIdx dimgidx = imgidx * (_connectivity / 2);
//                 double minAlpha = INFINITY;
//                 ImgIdx minIdx = dimgidx;

//                 if (bottom && dimg[dimgidx] < minAlpha) {
//                     minAlpha = dimg[dimgidx];
//                     minIdx = dimgidx;
//                 }
//                 if (right && dimg[dimgidx + 1] < minAlpha) {
//                     minAlpha = dimg[dimgidx + 1];
//                     minIdx = dimgidx + 1;
//                 }
//                 if (left && dimg[dimgidx - 1] < minAlpha) {
//                     minAlpha = dimg[dimgidx - 1];
//                     minIdx = dimgidx - 1;
//                 }
//                 if (top && dimg[dimgidx] < minAlpha) {
//                     minAlpha = dimg[dimgidx - width2];
//                     minIdx = dimgidx - width2;
//                 }
//                 edgeStatus[minIdx] = QItem::EDGE_ESSENTIAL;

//                 //                 dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + _width);
//                 // #pragma omp atomic
//                 //                 dhist[HHPQ::alphaToLevel(dimg[dimgidx], a)]++;

//                 //                 dimg[dimgidx + 1] = _pixelDissim.computeDissimilarity(imgidx, imgidx + 1);
//                 // #pragma omp atomic
//                 //                 dhist[HHPQ::alphaToLevel(dimg[dimgidx + 1], a)]++;
//             }
//         }
// }

template <class Pixel>
void AlphaTree<Pixel>::compute_dimg_hhpq_par(double *dimg, ImgIdx *dhist, const Pixel *img, double a) {

    if (_connectivity == 4) {
#pragma omp parallel for
        for (ImgIdx i = 0; i < _height; i++) {
            for (ImgIdx j = 0; j < _width; j++) {
                const ImgIdx imgidx = i * _width + j;
                const ImgIdx dimgidx = imgidx * (_connectivity / 2);
                if (i < _height - 1) {
                    dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + _width);
#pragma omp atomic
                    dhist[HHPQ::alphaToLevel(dimg[dimgidx], a)]++;
                }
                if (j < _width - 1) {
                    dimg[dimgidx + 1] = _pixelDissim.computeDissimilarity(imgidx, imgidx + 1);
#pragma omp atomic
                    dhist[HHPQ::alphaToLevel(dimg[dimgidx + 1], a)]++;
                }
            }
        }
    } else if (_connectivity == 8) {
        //   -  -  3
        //   -  p  2
        //   -  0  1
        // top,middle
#pragma omp parallel for
        for (ImgIdx i = 0; i < _height; i++) {
            const bool top = i > 0;
            const bool bottom = i < _height - 1;
            for (ImgIdx j = 0; j < _width; j++) {
                const bool right = j < _width - 1;
                const ImgIdx imgidx = i * _width + j;
                const ImgIdx dimgidx = imgidx * (_connectivity / 2);
                if (bottom) {
                    dimg[dimgidx] = _pixelDissim.computeDissimilarity(imgidx, imgidx + _width);
#pragma omp atomic
                    dhist[HHPQ::alphaToLevel(dimg[dimgidx], a)]++;
                }
                if (bottom && right) {
                    dimg[dimgidx + 1] = _pixelDissim.computeDissimilarity(imgidx, imgidx + _width + 1);
#pragma omp atomic
                    dhist[HHPQ::alphaToLevel(dimg[dimgidx + 1], a)]++;
                }
                if (right) {
                    dimg[dimgidx + 2] = _pixelDissim.computeDissimilarity(imgidx, imgidx + 1);
#pragma omp atomic
                    dhist[HHPQ::alphaToLevel(dimg[dimgidx + 2], a)]++;
                }
                if (top && right) {
                    dimg[dimgidx + 3] = _pixelDissim.computeDissimilarity(imgidx, imgidx - _width + 1);
#pragma omp atomic
                    dhist[HHPQ::alphaToLevel(dimg[dimgidx + 3], a)]++;
                }
            }
        }
    }
}

template <class Pixel> void AlphaTree<Pixel>::compute_dimg(Pixel *dimg, ImgIdx *dhist, const Pixel *img, double a) {
    ImgIdx dimgidx, imgidx, stride_w = _width, i, j;
    int hidx;
    imgidx = dimgidx = 0;

    if (_connectivity == 4) {
        for (i = 0; i < _height - 1; i++) {
            for (j = 0; j < _width - 1; j++) {
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
        for (j = 0; j < _width - 1; j++) {
            dimgidx++;
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx]));
            hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
            dhist[hidx]++;
            imgidx++;
        }
    } else if (_connectivity == 8) {
        //   -  -  3
        //   -  p  2
        //   -  0  1
        // top,middle
        for (i = 0; i < _height - 1; i++) {
            for (j = 0; j < _width - 1; j++) {
                dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + _width], img[imgidx])); // 0
                hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
                dhist[hidx]++;
                dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + _width + 1], img[imgidx])); // 1
                hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
                dhist[hidx]++;
                dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx])); // 2
                hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
                dhist[hidx]++;
                if (i > 0) {
                    dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx - _width + 1], img[imgidx])); // 3
                    hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
                    dhist[hidx]++;
                }
                dimgidx++;
                imgidx++;
            }
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + _width], img[imgidx])); // 0
            hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
            dhist[hidx]++;
            dimgidx += 4; // skip 1,2,3
            imgidx++;
        }

        // bottom
        dimgidx += 2; // skip 0,1
        for (j = 0; j < _width - 1; j++) {
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx])); // 2
            hidx = (int)(a * log2(1 + (double)dimg[dimgidx++]));
            dhist[hidx]++;
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx - _width + 1], img[imgidx])); // 3
            hidx = (int)(a * log2(1 + (double)dimg[dimgidx]));
            dhist[hidx]++;
            dimgidx += 3;
            imgidx++;
        }
    }
}

template <class Pixel> Pixel AlphaTree<Pixel>::compute_dimg(Pixel *dimg, ImgIdx *dhist, const Pixel *img) {
    ImgIdx dimgidx, imgidx, stride_w = _width, i, j;
    Pixel dmax = 0;

    imgidx = dimgidx = 0;
    if (_connectivity == 4) {
        for (i = 0; i < _height - 1; i++) {
            for (j = 0; j < _width - 1; j++) {
                dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + stride_w] - (int64_t)img[imgidx]));
                dmax = _max(dmax, dimg[dimgidx]);
                dhist[dimg[dimgidx++]]++;
                dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + 1] - (int64_t)img[imgidx]));
                dmax = _max(dmax, dimg[dimgidx]);
                dhist[dimg[dimgidx++]]++;
                imgidx++;
            }
            dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + stride_w] - (int64_t)img[imgidx]));
            dhist[dimg[dimgidx++]]++;
            dmax = _max(dmax, dimg[dimgidx]);
            dimgidx++;
            imgidx++;
        }
        for (j = 0; j < _width - 1; j++) {
            dimgidx++;
            dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + 1] - (int64_t)img[imgidx]));
            dmax = _max(dmax, dimg[dimgidx]);
            dhist[dimg[dimgidx++]]++;
            imgidx++;
        }
    } else if (_connectivity == 8) {
        //   -  -  3
        //   -  p  2
        //   -  0  1
        // top,middle
        for (i = 0; i < _height - 1; i++) {
            for (j = 0; j < _width - 1; j++) {
                dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + _width] - (int64_t)img[imgidx])); // 0
                dmax = _max(dmax, dimg[dimgidx]);
                dhist[dimg[dimgidx++]]++;
                dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + _width + 1] - (int64_t)img[imgidx])); // 1
                dmax = _max(dmax, dimg[dimgidx]);
                dhist[dimg[dimgidx++]]++;
                dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + 1] - (int64_t)img[imgidx])); // 2
                dmax = _max(dmax, dimg[dimgidx]);
                dhist[dimg[dimgidx++]]++;
                if (i > 0) {
                    dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx - _width + 1] - (int64_t)img[imgidx])); // 3
                    dmax = _max(dmax, dimg[dimgidx]);
                    dhist[dimg[dimgidx]]++;
                }
                dimgidx++;
                imgidx++;
            }
            dimg[dimgidx] = (Pixel)(abs((int64_t)img[imgidx + _width] - (int64_t)img[imgidx])); // 0
            dmax = _max(dmax, dimg[dimgidx]);
            dhist[dimg[dimgidx]]++;
            dimgidx += 4; // skip 1,2,3
            imgidx++;
        }

        // bottom
        dimgidx += 2; // skip 0,1
        for (j = 0; j < _width - 1; j++) {
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx + 1], img[imgidx])); // 2
            dmax = _max(dmax, dimg[dimgidx]);
            dhist[dimg[dimgidx++]]++;
            dimg[dimgidx] = (Pixel)(abs_diff(img[imgidx - _width + 1], img[imgidx])); // 3
            dmax = _max(dmax, dimg[dimgidx]);
            dhist[dimg[dimgidx]]++;
            dimgidx += 3;
            imgidx++;
        }
    }
    return dmax;
}

template <class Pixel> double AlphaTree<Pixel>::compute_dimg(double *dimg, ImgIdx *dhist, const Pixel *img) {
    ImgIdx dimgidx, imgidx, stride_w = _width, i, j;
    Pixel d;
    double dmax = 0;

    imgidx = dimgidx = 0;
    if (_connectivity == 4) {
        for (i = 0; i < _height - 1; i++) {
            for (j = 0; j < _width - 1; j++) {
                d = abs_diff(img[imgidx + stride_w], img[imgidx]);
                dimg[dimgidx] = (double)d;
                dmax = _max(dmax, dimg[dimgidx]);
                dhist[d]++;
                dimgidx++;
                d = abs_diff(img[imgidx + 1], img[imgidx]);
                dimg[dimgidx] = (double)d;
                dmax = _max(dmax, dimg[dimgidx]);
                dhist[d]++;
                dimgidx++;
                imgidx++;
            }
            d = abs_diff(img[imgidx + stride_w], img[imgidx]);
            dimg[dimgidx] = (double)d;
            dmax = _max(dmax, dimg[dimgidx]);
            dhist[d]++;
            dimgidx += 2;
            imgidx++;
        }
        for (j = 0; j < _width - 1; j++) {
            dimgidx++;
            d = abs_diff(img[imgidx + 1], img[imgidx]);
            dimg[dimgidx] = (double)d;
            dmax = _max(dmax, dimg[dimgidx]);
            dhist[d]++;
            dimgidx++;
            imgidx++;
        }
    } else if (_connectivity == 8) {
        //   -  -  3
        //   -  p  2
        //   -  0  1
        // top,middle
        for (i = 0; i < _height - 1; i++) {
            for (j = 0; j < _width - 1; j++) {
                d = (Pixel)(abs_diff((int64_t)img[imgidx + _width], (int64_t)img[imgidx])); // 0
                dhist[d]++;
                dimg[dimgidx++] = (double)d;
                dmax = _max(dmax, dimg[dimgidx]);
                d = (Pixel)(abs_diff((int64_t)img[imgidx + _width + 1], (int64_t)img[imgidx])); // 1
                dhist[d]++;
                dimg[dimgidx++] = (double)d;
                dmax = _max(dmax, dimg[dimgidx]);
                d = (Pixel)(abs_diff((int64_t)img[imgidx + 1], (int64_t)img[imgidx])); // 2
                dhist[d]++;
                dimg[dimgidx++] = (double)d;
                dmax = _max(dmax, dimg[dimgidx]);
                if (i > 0) {
                    d = (Pixel)(abs_diff((int64_t)img[imgidx - _width + 1], (int64_t)img[imgidx])); // 3
                    dhist[d]++;
                    dimg[dimgidx] = (double)d;
                    dmax = _max(dmax, dimg[dimgidx]);
                }
                dimgidx++;
                imgidx++;
            }
            d = (Pixel)(abs_diff((int64_t)img[imgidx + _width], (int64_t)img[imgidx])); // 0
            dhist[d]++;
            dimg[dimgidx] = (double)d;
            dmax = _max(dmax, dimg[dimgidx]);
            dimgidx += 4; // skip 1,2,3
            imgidx++;
        }

        // bottom
        dimgidx += 2; // skip 0,1
        for (j = 0; j < _width - 1; j++) {
            d = (Pixel)(abs_diff((int64_t)img[imgidx + 1], (int64_t)img[imgidx])); // 3
            dhist[d]++;
            dimg[dimgidx] = (double)d;
            dmax = _max(dmax, dimg[dimgidx]);
            dimgidx += 4;
            imgidx++;
        }
    }
    return dmax;
}

template <class Pixel> void AlphaTree<Pixel>::set_isAvailable(uint8_t *isAvailable) {
    int32_t i, j, k;
    int32_t imgSize = _width * _height;

    if (_connectivity == 4) {
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
        for (i = 0; i < imgSize; i++)
            isAvailable[i] = 0xff;

        j = _width * (_height - 1);
        for (i = 0; i < _width; i++) {
            isAvailable[i] &= 0x07;
            isAvailable[j] &= 0x0e;
            j++;
        }

        j = 0;
        k = _width - 1;
        for (i = 0; i < _height; i++) {
            isAvailable[j] &= 0xb;
            isAvailable[k] &= 0xd;
            j += _width;
            k += _width;
        }
    } else {
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

        // initialize to all available
        for (i = 0; i < imgSize; i++)
            isAvailable[i] = 0b11111111;

        // top and bottom row
        for (i = 0; i < _width; i++)
            isAvailable[i] &= 0b11000111;

        for (i = _width * (_height - 1); i < imgSize; i++)
            isAvailable[i] &= 0b01111100;

        // leftest and rightest column
        j = 0;
        k = _width - 1;
        for (i = 0; i < _height; i++) {
            isAvailable[j] &= 0b00011111;
            isAvailable[k] &= 0b11110001;
            j += _width;
            k += _width;
        }
    }
}

template <class Pixel>
void AlphaTree<Pixel>::set_isAvailable(uint8_t *isAvailable, int npartitions_hor, int npartitions_ver) {
    int32_t i, j, k;
    ImgIdx imgSize = _width * _height;
    ImgIdx wstride = _width / npartitions_ver;
    ImgIdx hstride = _height / npartitions_hor;

    set_isAvailable(isAvailable);

    if (_connectivity == 4) {
        // hor partitions
        j = (hstride - 1) * _width;
        for (i = 0; i < npartitions_hor - 1; i++) {
            k = j + _width;
            for (; j < k; j++) {
                isAvailable[j] &= 0xe;
                isAvailable[j + _width] &= 0x7;
            }

            j += (hstride - 1) * _width;
        }

        // ver partitions

        for (i = 0; i < npartitions_ver - 1; i++) {
            j = (i + 1) * wstride - 1;
            for (; j < imgSize; j += _width) {
                isAvailable[j] &= 0xd;
                isAvailable[j + 1] &= 0xb;
            }
        }
    } else {
    }
}

template <class Pixel> uint8_t AlphaTree<Pixel>::is_available(uint8_t isAvailable, uint8_t iNeighbour) const {
    return (isAvailable >> iNeighbour) & 1;
}

template <class Pixel> void AlphaTree<Pixel>::set_field(uint8_t *arr, ImgIdx idx, uint8_t in) { arr[idx] = in; }

template <class Pixel> uint8_t AlphaTree<Pixel>::get_field(uint8_t *arr, ImgIdx idx) { return arr[idx]; }

template <class Pixel> void AlphaTree<Pixel>::connectPix2Node(ImgIdx pidx, Pixel pix_val, ImgIdx iNode, Pixel level) {
    AlphaNode<Pixel> *pNode;
    pNode = _node + iNode;
    _parentAry[pidx] = iNode;
    if (pNode->area) // possibly unnecessary branch..
        pNode->add(pix_val);
    else
        pNode->set(1, level, (double)pix_val, pix_val, pix_val);
}

template <class Pixel> void AlphaTree<Pixel>::connectPix2Node(ImgIdx pidx, Pixel pix_val, ImgIdx iNode) {
    AlphaNode<Pixel> *pNode = &_node[iNode];
    _parentAry[pidx] = iNode;
    pNode->add(pix_val);
}

template <class Pixel> void AlphaTree<Pixel>::connectPix2Node0(ImgIdx pidx, Pixel pix_val, ImgIdx iNode, Pixel level) {
    AlphaNode<Pixel> *pNode;
    pNode = _node + iNode;
    _parentAry[pidx] = iNode;
    pNode->set(1, level, (double)pix_val, pix_val, pix_val);
}

template <class Pixel> ImgIdx AlphaTree<Pixel>::NewAlphaNode() {
    if (_curSize == _maxSize) {
        std::cout << "Reallocating...\n";
        _maxSize = _min((1 + (_connectivity >> 1)) * _height * _width, _maxSize + (ImgIdx)(2 * _height * _width * 0.1));

        _node = (AlphaNode<Pixel> *)Realloc(_node, _maxSize * sizeof(AlphaNode<Pixel>));
    }
    return _curSize++;
}

template <class Pixel> ImgIdx AlphaTree<Pixel>::NewAlphaNode(Pixel level, AlphaNode<Pixel> *pCopy) {
    AlphaNode<Pixel> *pNew = _node + _curSize;

    if (_curSize == _maxSize) {
        std::cout << "Reallocating...\n";
        _maxSize = _min((1 + (_connectivity >> 1)) * _height * _width, _maxSize + (ImgIdx)(2 * _height * _width * 0.1));

        _node = (AlphaNode<Pixel> *)Realloc(_node, _maxSize * sizeof(AlphaNode<Pixel>));
        pNew = _node + _curSize;
    }
    pNew->alpha = level;
    pNew->copy(pCopy);
    return _curSize++;
}

template <class Pixel> ImgIdx AlphaTree<Pixel>::NewAlphaNode1(double level, AlphaNode<Pixel> *pCopy) {
    AlphaNode<Pixel> *pNew = _node + _curSize;

    if (_curSize == _maxSize) {
        std::cout << "Reallocating...\n";
        _maxSize = _min((1 + (_connectivity >> 1)) * _height * _width, _maxSize + (ImgIdx)(2 * _height * _width * 0.1));

        _node = (AlphaNode<Pixel> *)Realloc(_node, _maxSize * sizeof(AlphaNode<Pixel>));
        pNew = _node + _curSize;
    }
    pNew->alpha = level;
    pNew->copy(pCopy);
    return _curSize++;
}

template <class Pixel>
ImgIdx AlphaTree<Pixel>::NewAlphaNode(Pixel level) // Fix it later - no need to initialize
{
    AlphaNode<Pixel> *pNew = _node + _curSize;

    if (_curSize == _maxSize) {
        std::cout << "Reallocating...\n";
        _maxSize = _min((1 + (_connectivity >> 1)) * _height * _width, _maxSize + (ImgIdx)(_height * _width * 0.1));

        _node = (AlphaNode<Pixel> *)Realloc(_node, _maxSize * sizeof(AlphaNode<Pixel>));
        pNew = _node + _curSize;
    }
    pNew->alpha = level;
    pNew->minPix = (uint8_t)-1;
    pNew->maxPix = 0;
    pNew->sumPix = 0.0;
    //		pNew->parentIdx = 0;
    pNew->area = 0;

    return _curSize++;
}

template <class Pixel> uint8_t AlphaTree<Pixel>::is_visited(uint8_t *isVisited, ImgIdx p) { return isVisited[p]; }

template <class Pixel> void AlphaTree<Pixel>::visit(uint8_t *isVisited, ImgIdx p) { isVisited[p] = 1; }

template <class Pixel>
ImgIdx AlphaTree<Pixel>::TreeSizeEstimation(ImgIdx *dhist, int64_t numlevels, ImgIdx imgSize, ImgIdx nredges) {
    return TreeSizeEstimation(dhist, numlevels, imgSize, nredges, TSE_M);
}

template <class Pixel>
ImgIdx AlphaTree<Pixel>::TreeSizeEstimation(ImgIdx *dhist, int64_t numlevels, ImgIdx imgSize, ImgIdx nredges,
                                            double m) {
    if (imgSize < TSE_MINSIZE)
        return 10 + nredges + imgSize;
    double tse_nrmsd = 0;
    for (int64_t p = 0; p < numlevels; p++)
        tse_nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
    tse_nrmsd = sqrt((tse_nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
    nrmsd = tse_nrmsd;
    ImgIdx ret = _min(2 * imgSize, (ImgIdx)(2 * imgSize * ((TSE_A * exp(TSE_SIGMA * tse_nrmsd) + TSE_B) + m)));
    double dret =
        _min((double)(2 * imgSize), (double)(2 * imgSize * ((TSE_A * exp(TSE_SIGMA * tse_nrmsd) + TSE_B) + m)));

    if (ret < 0) {
        printf("Warning: TSE returned 0< value\n");
        printf("nrmsd = %lf\n", tse_nrmsd);
        printf(" exp(SIGMA * nrmsd) = %lf\n", exp(TSE_SIGMA * tse_nrmsd));
        printf("A * exp(SIGMA * nrmsd) + B = %lf\n", TSE_A * exp(TSE_SIGMA * tse_nrmsd) + TSE_B);
        printf("2 * imgSize * ((A * exp(SIGMA * nrmsd) + B) + m) = %lf\n",
               2 * imgSize * ((TSE_A * exp(TSE_SIGMA * tse_nrmsd) + TSE_B) + m));
        printf("(ImgIdx)(2 * imgSize * ((A * exp(SIGMA * nrmsd) + B) + m)) = %d\n",
               (int)(ImgIdx)(2 * imgSize * ((TSE_A * exp(TSE_SIGMA * tse_nrmsd) + TSE_B) + m)));
        printf("nredges = %lf\n", (double)nredges);
        printf("imgSize = %lf\n", (double)imgSize);
        printf("1 + nredges + imgSize = %lf\n", (double)(1 + nredges + imgSize));
        printf("ret = %lf\n", (double)ret);
        printf("ret = %lf\n", (double)dret);
    }

    return ret;
}

template <class Pixel>
ImgIdx AlphaTree<Pixel>::TreeSizeEstimation(ImgIdx *dhist, int64_t numlevels, ImgIdx imgSize, ImgIdx nredges, double m,
                                            ImgIdx reserve) {
    nrmsd = 0;
    for (int64_t p = 0; p < numlevels; p++)
        nrmsd += ((double)dhist[p]) * ((double)dhist[p]);
    nrmsd = sqrt((nrmsd - (double)nredges) / ((double)nredges * ((double)nredges - 1.0)));
    return _min(2.0 * imgSize, (ImgIdx)(2.0 * (double)imgSize * ((TSE_A * exp(TSE_SIGMA * nrmsd) + TSE_B) + m))) +
           reserve;
}

template <class Pixel> void AlphaTree<Pixel>::remove_redundant_node(ImgIdx &prevTop, ImgIdx &stackTop) {
    if (_node[prevTop].parentIdx == stackTop && _node[prevTop].area == _node[stackTop].area) {
        _node[prevTop].parentIdx = _node[stackTop].parentIdx;
        stackTop = prevTop;
        _curSize--;
    }
}

template <class Pixel> void AlphaTree<Pixel>::FloodHierarQueueNoCache(const Pixel *img, int tse) {
    if (sizeof(Pixel) > 2 || _channel > 1) {
        printf("Error: Hierarchical queues do not work on >16 bits images or multispectral images\n");
        printf("Try Unionfind (algorithm code %d), flooding using Heapqueue (%d), trie queue (%d) or cached trie queue "
               "(%d) \n",
               UNIONFIND, FLOOD_HEAPQUEUE_CACHE, FLOOD_TRIE, FLOOD_TRIE_CACHE);
        return;
    }

    ImgIdx imgSize, dimgSize, nredges, x0;
    uint64_t numlevels, max_level, currentLevel;
    ImgIdx *dhist;
    ImgIdx prevTop, stackTop, iNode;
    uint8_t *isVisited, *isAvailable, isAv;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    max_level = (sizeof(Pixel) == 8) ? 0xffffffffffffffff : (int64_t)((Pixel)(-1));
    numlevels = max_level + 1;

    Pixel *dimg;
    dhist = (ImgIdx *)Malloc((size_t)numlevels * sizeof(ImgIdx));
    memset(dhist, 0, (size_t)numlevels * sizeof(int32_t));
    dimg = (Pixel *)Malloc((size_t)dimgSize * sizeof(Pixel));
    compute_dimg(dimg, dhist, img); // calculate pixel differences and make histogram

    // create hierarchical queue from dhist
    HierarQueue *queue = new HierarQueue(nredges + 1, dhist, numlevels); // +1 for the dummy _node
    _curSize = 0;

    if (!tse || imgSize < 10000 || sizeof(Pixel) > 1) // for small imags do not use TSE
        _maxSize = 1 + imgSize + nredges;
    else
        _maxSize = TreeSizeEstimation(dhist, numlevels, imgSize, nredges);

    if (dhist)
        Free(dhist);

    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable(isAvailable);

    _parentAry = (ImgIdx *)Malloc((size_t)imgSize * sizeof(int32_t));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));

    stackTop = NewAlphaNode(); /*dummy root*/
    AlphaNode<Pixel> *pNode = _node + stackTop;
    pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0);
    pNode->parentIdx = stackTop;
    currentLevel = max_level;
    x0 = 0; /*arbitrary starting point*/
    prevTop = stackTop;

    queue->push(x0, currentLevel);
    while (1) // flooding
    {
        while ((uint64_t)queue->min_level <= currentLevel) // flood all levels below currentLevel
        {
            p = queue->pop();
            if (isVisited[p]) {
                queue->find_minlev();
                continue;
            }
            isVisited[p] = 1;
            isAv = isAvailable[p];
            if (_connectivity == 4) {
                q = p << 1;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - (_width << 1)]);
            } else if (_connectivity == 8) {
                ImgIdx width4 = _width << 2;
                q = p << 2;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])
                    queue->push(p + _width + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 2]);
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])
                    queue->push(p - _width + 1, dimg[q + 3]);
                if (is_available(isAv, 4) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - width4]);
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])
                    queue->push(p - _width - 1, dimg[q - width4 - 3]);
                if (is_available(isAv, 6) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 2]);
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])
                    queue->push(p + _width - 1, dimg[q + width4 - 1]);
            } else {
                //?
            }
            // else //later

            if (currentLevel > (uint64_t)queue->min_level) // go to lower level
            {

                { // creat new _node
                    Pixel pix_val = img[p];
                    currentLevel = queue->min_level;
                    iNode = NewAlphaNode();
                    _node[iNode].set(1, currentLevel, (double)pix_val, pix_val, pix_val);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    stackTop = iNode;
                }
                if (currentLevel) {
                    iNode = NewAlphaNode(0, _node + stackTop);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                    prevTop = iNode;
                } else
                    _parentAry[p] = stackTop;
            } else {
                queue->find_minlev();

                if (currentLevel) {
                    Pixel pix_val = img[p];
                    iNode = NewAlphaNode();
                    _node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val);
                    _node[stackTop].add(_node + iNode);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                } else
                    connectPix2Node(p, img[p], stackTop);
            }
        }

        remove_redundant_node(prevTop, stackTop);

        // go to higher level
        iNode = _node[stackTop].parentIdx;
        if ((Pixel)queue->min_level < (Pixel)_node[iNode].alpha) // new level from queue
        {
            iNode = NewAlphaNode(queue->min_level, _node + stackTop);
            _node[iNode].parentIdx = _node[stackTop].parentIdx;
            _node[iNode]._rootIdx = ROOTIDX;
            _node[stackTop].parentIdx = iNode;
        } else // go to existing _node
        {
            if (_node[iNode].area == imgSize) // root _node found...done
                break;
            _node[iNode].add(_node + stackTop);
        }

        if (_node[iNode].area == imgSize) // root _node found...done
            break;

        prevTop = stackTop;
        stackTop = iNode;
        currentLevel = (uint64_t)_node[stackTop].alpha;
    }
    _rootIdx = (_node[stackTop].area == imgSize) ? stackTop : iNode; // remove redundant root
    _node[_rootIdx].parentIdx = ROOTIDX;

    delete queue;
    Free(dimg);
    Free(isVisited);
    Free(isAvailable);
}

template <class Pixel> void AlphaTree<Pixel>::FloodHeapQueueNoCache(const Pixel *img) {
    HeapQueue<double> *queue;

    ImgIdx imgSize, dimgSize, nredges, x0;
    uint64_t numlevels, max_level;
    ImgIdx *dhist;
    ImgIdx stackTop, prevTop, iNode;
    uint8_t *isVisited, *isAvailable, isAv;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    max_level = (Pixel)(-1);
    numlevels = max_level + 1;
    double currentLevel;
    double *dimg;

    dimg = (double *)Malloc((size_t)dimgSize * sizeof(double));

    if (sizeof(Pixel) == 1 && _channel == 1) {
        dhist = (ImgIdx *)Calloc(numlevels * sizeof(ImgIdx));
        compute_dimg(dimg, dhist, img); // calculate pixel differences and make histogram
        _maxSize = TreeSizeEstimation(dhist, numlevels, imgSize, nredges);
    } else {
        dhist = 0;
        compute_dimg(dimg, img);
        _maxSize = 1 + imgSize + nredges;
    }

    // create heap-based priority queue
    queue = new HeapQueue<double>(nredges);
    _curSize = 0;

    if (dhist)
        Free(dhist);
    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable(isAvailable);
    _parentAry = (ImgIdx *)Malloc((size_t)imgSize * sizeof(int32_t));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));

    x0 = 0; /*arbitrary starting point*/
    currentLevel = DBL_MAX;
    stackTop = NewAlphaNode(); /*dummy root*/
    AlphaNode<Pixel> *pNode = _node + stackTop;
    pNode->set(0, currentLevel, (double)0.0, (Pixel)max_level, (Pixel)0);
    pNode->parentIdx = stackTop;
    // currentLevel = max_level;
    prevTop = stackTop; /*to find redundant _node*/

    queue->push_run(x0, currentLevel);
    while (1) {
        while (queue->get_minlev() <= currentLevel) // flood all levels below currentLevel
        {
            p = queue->top();
            if (isVisited[p]) {
                queue->pop();
                continue;
            }
            isVisited[p] = 1;

            isAv = isAvailable[p];
            if (_connectivity == 4) {
                q = p << 1;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - (_width << 1)]);
            } else if (_connectivity == 8) {
                ImgIdx width4 = _width << 2;
                q = p << 2;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])
                    queue->push(p + _width + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 2]);
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])
                    queue->push(p - _width + 1, dimg[q + 3]);
                if (is_available(isAv, 4) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - width4]);
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])
                    queue->push(p - _width - 1, dimg[q - width4 - 3]);
                if (is_available(isAv, 6) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 2]);
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])
                    queue->push(p + _width - 1, dimg[q + width4 - 1]);
            } else {
                //?
            }

            queue->find_minlev();
            if (currentLevel > queue->get_minlev()) {
                Pixel pix_val = img[p];
                currentLevel = queue->get_minlev();
                iNode = NewAlphaNode();
                _node[iNode].set(1, currentLevel, (double)pix_val, pix_val, pix_val);
                _node[iNode].parentIdx = stackTop;
                _node[iNode]._rootIdx = ROOTIDX;
                stackTop = iNode;
                if (currentLevel) {
                    iNode = NewAlphaNode(0, _node + stackTop);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                    prevTop = iNode;
                } else
                    _parentAry[p] = stackTop;
            } else {
                if (currentLevel) {
                    Pixel pix_val = img[p];
                    iNode = NewAlphaNode();
                    _node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val);
                    _node[stackTop].add(_node + iNode);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                } else
                    connectPix2Node(p, img[p], stackTop);
                if (_node[stackTop].area == imgSize)
                    goto FLOOD_END;
            }
        }

        remove_redundant_node(prevTop, stackTop);

        if (_node[stackTop].area == imgSize) // root _node found...done
            break;

        iNode = _node[stackTop].parentIdx;
        if (queue->get_minlev() < _node[iNode].alpha) {
            iNode = NewAlphaNode1(queue->get_minlev(), _node + stackTop);
            _node[iNode].parentIdx = _node[stackTop].parentIdx;
            _node[iNode]._rootIdx = ROOTIDX;
            _node[stackTop].parentIdx = iNode;
        } else
            _node[iNode].add(_node + stackTop);

        prevTop = stackTop;
        stackTop = iNode;
        currentLevel = _node[stackTop].alpha;
        if (_node[stackTop].area == imgSize) // root _node found...done
            break;
    }
FLOOD_END:
    _node[stackTop].parentIdx = ROOTIDX;

    delete queue;
    Free(dimg);
    Free(isVisited);
    Free(isAvailable);
}

template <class Pixel> void AlphaTree<Pixel>::FloodHeapQueueNaiveNoCache(const Pixel *img) {
    HeapQueue_naive<double> *queue;

    ImgIdx imgSize, dimgSize, nredges, x0;
    uint64_t numlevels, max_level;
    ImgIdx *dhist;
    ImgIdx stackTop, prevTop, iNode;
    uint8_t *isVisited, *isAvailable, isAv;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    max_level = (Pixel)(-1);
    numlevels = max_level + 1;
    double currentLevel;
    double *dimg;

    dimg = (double *)Malloc((size_t)dimgSize * sizeof(double));

    if (sizeof(Pixel) == 1 && _channel == 1) {
        dhist = (ImgIdx *)Calloc(numlevels * sizeof(ImgIdx));
        compute_dimg(dimg, dhist, img); // calculate pixel differences and make histogram
        _maxSize = TreeSizeEstimation(dhist, numlevels, imgSize, nredges);
    } else {
        dhist = 0;
        compute_dimg(dimg, img);
        _maxSize = 1 + imgSize + nredges;
    }

    // create heap-based priority queue
    queue = new HeapQueue_naive<double>(nredges + 1);
    _curSize = 0;

    if (dhist)
        Free(dhist);
    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable(isAvailable);
    _parentAry = (ImgIdx *)Malloc((size_t)imgSize * sizeof(int32_t));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));

    x0 = 0; /*arbitrary starting point*/
    currentLevel = DBL_MAX;
    stackTop = NewAlphaNode(); /*dummy root*/
    AlphaNode<Pixel> *pNode = _node + stackTop;
    pNode->set(0, currentLevel, (double)0.0, (Pixel)max_level, (Pixel)0);
    pNode->parentIdx = stackTop;
    // currentLevel = max_level;
    prevTop = stackTop; /*to find redundant _node*/

    queue->push(x0, currentLevel);
    while (1) {
        while (queue->get_minlev() <= currentLevel) // flood all levels below currentLevel
        {
            p = queue->top();
            if (isVisited[p]) {
                queue->pop();
                continue;
            }
            isVisited[p] = 1;

            isAv = isAvailable[p];
            if (_connectivity == 4) {
                q = p << 1;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - (_width << 1)]);
            } else if (_connectivity == 8) {
                ImgIdx width4 = _width << 2;
                q = p << 2;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])
                    queue->push(p + _width + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 2]);
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])
                    queue->push(p - _width + 1, dimg[q + 3]);
                if (is_available(isAv, 4) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - width4]);
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])
                    queue->push(p - _width - 1, dimg[q - width4 - 3]);
                if (is_available(isAv, 6) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 2]);
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])
                    queue->push(p + _width - 1, dimg[q + width4 - 1]);
            } else {
                //?
            }

            if (currentLevel > queue->get_minlev()) {
                Pixel pix_val = img[p];
                currentLevel = queue->get_minlev();
                iNode = NewAlphaNode();
                _node[iNode].set(1, currentLevel, (double)pix_val, pix_val, pix_val);
                _node[iNode].parentIdx = stackTop;
                _node[iNode]._rootIdx = ROOTIDX;
                stackTop = iNode;
                if (currentLevel) {
                    iNode = NewAlphaNode(0, _node + stackTop);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                    prevTop = iNode;
                } else
                    _parentAry[p] = stackTop;
            } else {
                if (currentLevel) {
                    Pixel pix_val = img[p];
                    iNode = NewAlphaNode();
                    _node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val);
                    _node[stackTop].add(_node + iNode);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                } else
                    connectPix2Node(p, img[p], stackTop);
                if (_node[stackTop].area == imgSize)
                    goto FLOOD_END;
            }
        }

        remove_redundant_node(prevTop, stackTop);

        if (_node[stackTop].area == imgSize) // root _node found...done
            break;

        iNode = _node[stackTop].parentIdx;
        if (queue->get_minlev() < _node[iNode].alpha) {
            iNode = NewAlphaNode1(queue->get_minlev(), _node + stackTop);
            _node[iNode].parentIdx = _node[stackTop].parentIdx;
            _node[iNode]._rootIdx = ROOTIDX;
            _node[stackTop].parentIdx = iNode;
        } else
            _node[iNode].add(_node + stackTop);

        prevTop = stackTop;
        stackTop = iNode;
        currentLevel = _node[stackTop].alpha;
        if (_node[stackTop].area == imgSize) // root _node found...done
            break;
    }
FLOOD_END:
    _node[stackTop].parentIdx = ROOTIDX;

    delete queue;
    Free(dimg);
    Free(isVisited);
    Free(isAvailable);
}

template <class Pixel> void AlphaTree<Pixel>::FloodHeapQueue(const Pixel *img) {
    Cache_Quad_Heapqueue<double> *queue;

    ImgIdx imgSize, dimgSize, nredges, x0;
    uint64_t numlevels, max_level;
    ImgIdx *dhist;
    ImgIdx stackTop, prevTop, iNode;
    uint8_t *isVisited, *isAvailable, isAv;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    max_level = (Pixel)(-1);
    numlevels = max_level + 1;
    double currentLevel;
    double *dimg;

    dimg = (double *)Malloc((size_t)dimgSize * sizeof(double));

    if (sizeof(Pixel) == 1 && _channel == 1) {
        dhist = (ImgIdx *)Calloc(numlevels * sizeof(ImgIdx));
        compute_dimg(dimg, dhist, img); // calculate pixel differences and make histogram
        _maxSize = TreeSizeEstimation(dhist, numlevels, imgSize, nredges);
    } else {
        dhist = 0;
        compute_dimg(dimg, img);
        _maxSize = 1 + imgSize + nredges;
    }

    // create heap-based priority queue
    queue = new Cache_Quad_Heapqueue<double>(nredges);
    _curSize = 0;

    if (dhist)
        Free(dhist);
    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable(isAvailable);
    _parentAry = (ImgIdx *)Malloc((size_t)imgSize * sizeof(int32_t));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));

    x0 = 0; /*arbitrary starting point*/
    currentLevel = DBL_MAX;
    stackTop = NewAlphaNode(); /*dummy root*/
    AlphaNode<Pixel> *pNode = _node + stackTop;
    pNode->set(0, currentLevel, (double)0.0, (Pixel)max_level, (Pixel)0);
    pNode->parentIdx = stackTop;
    prevTop = stackTop; /*to find redundant _node*/

    queue->push_1stitem(x0, currentLevel);
    while (1) {
        while (queue->get_minlev() <= currentLevel) // flood all levels below currentLevel
        {
            p = queue->top();
            if (isVisited[p]) {
                queue->pop();
                continue;
            }
            queue->startPushes();
            isVisited[p] = 1;

            isAv = isAvailable[p];
            if (_connectivity == 4) {
                q = p << 1;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - (_width << 1)]);
            } else if (_connectivity == 8) {
                ImgIdx width4 = _width << 2;
                q = p << 2;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])
                    queue->push(p + _width + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 2]);
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])
                    queue->push(p - _width + 1, dimg[q + 3]);
                if (is_available(isAv, 4) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - width4]);
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])
                    queue->push(p - _width - 1, dimg[q - width4 - 3]);
                if (is_available(isAv, 6) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 2]);
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])
                    queue->push(p + _width - 1, dimg[q + width4 - 1]);
            } else {
                //?
            }

            queue->endPushes();
            if (currentLevel > queue->get_minlev()) // remove typecasting later
            {
                Pixel pix_val = img[p];
                currentLevel = queue->get_minlev();
                iNode = NewAlphaNode();
                _node[iNode].set(1, currentLevel, (double)pix_val, pix_val, pix_val);
                _node[iNode].parentIdx = stackTop;
                _node[iNode]._rootIdx = ROOTIDX;
                stackTop = iNode;
                if (currentLevel) {
                    iNode = NewAlphaNode(0, _node + stackTop);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                    prevTop = iNode;
                } else
                    _parentAry[p] = stackTop;
            } else {
                if (currentLevel) {
                    Pixel pix_val = img[p];
                    iNode = NewAlphaNode();
                    _node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val);
                    _node[stackTop].add(_node + iNode);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                } else
                    connectPix2Node(p, img[p], stackTop);
                if (_node[stackTop].area == imgSize)
                    goto FLOOD_END;
            }
        }
        remove_redundant_node(prevTop, stackTop);

        if (_node[stackTop].area == imgSize) // root _node found...done
            break;

        iNode = _node[stackTop].parentIdx;
        if (queue->get_minlev() < _node[iNode].alpha) {
            iNode = NewAlphaNode1(queue->get_minlev(), _node + stackTop);
            _node[iNode].parentIdx = _node[stackTop].parentIdx;
            _node[iNode]._rootIdx = ROOTIDX;
            _node[stackTop].parentIdx = iNode;
        } else
            _node[iNode].add(_node + stackTop);

        prevTop = stackTop;
        stackTop = iNode;
        currentLevel = _node[stackTop].alpha;
        if (_node[stackTop].area == imgSize) // root _node found...done
            break;
    }
FLOOD_END:
    _node[stackTop].parentIdx = ROOTIDX;

    delete queue;
    Free(dimg);
    Free(isVisited);
    Free(isAvailable);
}

template <class Pixel> void AlphaTree<Pixel>::FloodHierarQueue(const Pixel *img) {
    HierarQueueCache<Pixel> *queue;
    ImgIdx imgSize, dimgSize, nredges, x0;
    uint64_t numlevels, max_level, currentLevel;
    ImgIdx *dhist;
    ImgIdx prevTop, stackTop, iNode;
    uint8_t *isVisited, *isAvailable, isAv;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    max_level = (sizeof(Pixel) == 8) ? 0xffffffffffffffff : (int64_t)((Pixel)(-1));
    numlevels = max_level + 1;

    Pixel *dimg;
    dhist = (ImgIdx *)Calloc((size_t)numlevels * sizeof(ImgIdx));
    dimg = (Pixel *)Malloc((size_t)dimgSize * sizeof(Pixel));

    max_level = compute_dimg1(dimg, dhist, img); // calculate pixel differences and make histogram
    numlevels = max_level + 1;

    // create hierarchical queue from dhist
    queue = new HierarQueueCache<Pixel>(nredges + 1, dhist, numlevels); // +1 for the dummy _node
    _curSize = 0;

    if (imgSize < 10000 || sizeof(Pixel) > 2) // for small imags do not use TSE
        _maxSize = 1 + imgSize + dimgSize;
    else
        _maxSize = TreeSizeEstimation(dhist, numlevels, imgSize, nredges);
    if (dhist)
        Free(dhist);

    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable(isAvailable);
    _parentAry = (ImgIdx *)Malloc((size_t)imgSize * sizeof(int32_t));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));

    stackTop = NewAlphaNode(); /*dummy root*/
    AlphaNode<Pixel> *pNode = _node + stackTop;
    pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0);
    pNode->parentIdx = stackTop;
    currentLevel = max_level;
    x0 = 0; /*arbitrary starting point*/
    prevTop = stackTop;

    queue->push_1stitem(x0, currentLevel);
    while (1) // flooding
    {
        while ((int32_t)queue->top_alpha() <= (int32_t)currentLevel) // flood all levels below currentLevel
        {
            p = queue->top();
            if (isVisited[p] == 1) {
                queue->pop();
                continue;
            }

            queue->startPushes();
            isVisited[p] = 1;
            isAv = isAvailable[p];
            if (_connectivity == 4) {
                q = p << 1;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - (_width << 1)]);
            } else if (_connectivity == 8) {
                ImgIdx width4 = _width << 2;
                q = p << 2;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])
                    queue->push(p + _width + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 2]);
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])
                    queue->push(p - _width + 1, dimg[q + 3]);
                if (is_available(isAv, 4) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - width4]);
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])
                    queue->push(p - _width - 1, dimg[q - width4 - 3]);
                if (is_available(isAv, 6) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 2]);
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])
                    queue->push(p + _width - 1, dimg[q + width4 - 1]);
            } else {
                //?
            }

            queue->endPushes();
            if (currentLevel > (uint64_t)queue->top_alpha()) // go to lower level
            {
                Pixel pix_val = img[p];
                currentLevel = queue->top_alpha();
                iNode = NewAlphaNode();
                _node[iNode].set(1, currentLevel, (double)pix_val, pix_val, pix_val);
                _node[iNode].parentIdx = stackTop;
                _node[iNode]._rootIdx = ROOTIDX;
                stackTop = iNode;
                if (currentLevel) {
                    iNode = NewAlphaNode(0, _node + stackTop);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                    prevTop = iNode;
                } else
                    _parentAry[p] = stackTop;
            } else {
                if (currentLevel) {
                    Pixel pix_val = img[p];
                    iNode = NewAlphaNode();
                    _node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val);
                    _node[stackTop].add(_node + iNode);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                } else
                    connectPix2Node(p, img[p], stackTop);
                if (_node[stackTop].area == imgSize)
                    goto FLOOD_END;
            }
        }

        remove_redundant_node(prevTop, stackTop);

        if (_node[stackTop].area == imgSize) // root _node found...done
            break;

        // go to higher level
        iNode = _node[stackTop].parentIdx;
        if ((int32_t)queue->top_alpha() < (int32_t)_node[iNode].alpha) {
            iNode = NewAlphaNode1(queue->top_alpha(), _node + stackTop);
            _node[iNode].parentIdx = _node[stackTop].parentIdx;
            _node[iNode]._rootIdx = ROOTIDX;
            _node[stackTop].parentIdx = iNode;
        } else // go to existing _node
        {
            _node[iNode].add(_node + stackTop);
        }

        prevTop = stackTop;
        stackTop = iNode;
        currentLevel = (int32_t)_node[stackTop].alpha;
        if (_node[stackTop].area == imgSize) // root _node found...done
            break;
    }

FLOOD_END:
    _rootIdx = (_node[stackTop].area == imgSize) ? stackTop : iNode; // remove redundant root
    _node[_rootIdx].parentIdx = ROOTIDX;

    delete queue;
    Free(dimg);
    Free(isVisited);
    Free(isAvailable);
}

template <class Pixel> void AlphaTree<Pixel>::FloodLadderQueue(const Pixel *img, int thres) {
    LadderQueue *queue;

    ImgIdx imgSize, dimgSize, nredges, x0;
    uint64_t numlevels, max_level;
    ImgIdx *dhist;
    ImgIdx stackTop, prevTop, iNode;
    uint8_t *isVisited, *isAvailable, isAv;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    max_level = (Pixel)(-1);
    numlevels = max_level + 1;
    double currentLevel;
    double *dimg;
    Pixel dmax;

    dimg = (double *)Malloc((size_t)dimgSize * sizeof(double));

    if (sizeof(Pixel) == 1 && _channel == 1) {
        dhist = (ImgIdx *)Calloc(numlevels * sizeof(ImgIdx));
        dmax = compute_dimg(dimg, dhist, img); // calculate pixel differences and make histogram
        _maxSize = TreeSizeEstimation(dhist, numlevels, imgSize, nredges);
    } else {
        dhist = 0;
        dmax = compute_dimg(dimg, img);
        _maxSize = 1 + imgSize + nredges;
    }

    // create heap-based priority queue
    queue = new LadderQueue(thres);
    _curSize = 0;

    if (dhist)
        Free(dhist);
    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable(isAvailable);
    _parentAry = (ImgIdx *)Malloc((size_t)imgSize * sizeof(int32_t));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));

    x0 = 0; /*arbitrary starting point*/
    currentLevel = dmax;
    stackTop = NewAlphaNode(); /*dummy root*/
    AlphaNode<Pixel> *pNode = _node + stackTop;
    pNode->set(0, currentLevel, (double)0.0, (Pixel)max_level, (Pixel)0);
    pNode->parentIdx = stackTop;
    // currentLevel = max_level;
    prevTop = stackTop; /*to find redundant _node*/

    queue->enqueue(x0, currentLevel);
    while (1) {
        while (queue->getTopAlpha() <= currentLevel) // flood all levels below currentLevel
        {
            p = queue->dequeue();
            if (isVisited[p])
                continue;
            isVisited[p] = 1;

            isAv = isAvailable[p];
            if (_connectivity == 4) {
                q = p << 1;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->enqueue(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])
                    queue->enqueue(p + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])
                    queue->enqueue(p - 1, dimg[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])
                    queue->enqueue(p - _width, dimg[q - (_width << 1)]);
            } else if (_connectivity == 8) {
                ImgIdx width4 = _width << 2;
                q = p << 2;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->enqueue(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])
                    queue->enqueue(p + _width + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p + 1])
                    queue->enqueue(p + 1, dimg[q + 2]);
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])
                    queue->enqueue(p - _width + 1, dimg[q + 3]);
                if (is_available(isAv, 4) && !isVisited[p - _width])
                    queue->enqueue(p - _width, dimg[q - width4]);
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])
                    queue->enqueue(p - _width - 1, dimg[q - width4 - 3]);
                if (is_available(isAv, 6) && !isVisited[p - 1])
                    queue->enqueue(p - 1, dimg[q - 2]);
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])
                    queue->enqueue(p + _width - 1, dimg[q + width4 - 1]);
            } else {
                //?
            }

            double top_level = queue->getTopAlpha();
            if (currentLevel > top_level) {
                Pixel pix_val = img[p];
                currentLevel = top_level;
                iNode = NewAlphaNode();
                _node[iNode].set(1, currentLevel, (double)pix_val, pix_val, pix_val);
                _node[iNode].parentIdx = stackTop;
                _node[iNode]._rootIdx = ROOTIDX;
                stackTop = iNode;
                if (currentLevel) {
                    iNode = NewAlphaNode(0, _node + stackTop);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                    prevTop = iNode;
                } else
                    _parentAry[p] = stackTop;
            } else {
                if (currentLevel) {
                    Pixel pix_val = img[p];
                    iNode = NewAlphaNode();
                    _node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val);
                    _node[stackTop].add(_node + iNode);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                } else
                    connectPix2Node(p, img[p], stackTop);
                if (_node[stackTop].area == imgSize)
                    goto FLOOD_END;
            }
        }

        remove_redundant_node(prevTop, stackTop);

        if (_node[stackTop].area == imgSize) // root _node found...done
            break;

        double top_level = queue->getTopAlpha();
        iNode = _node[stackTop].parentIdx;
        if (top_level < _node[iNode].alpha) {
            iNode = NewAlphaNode1(top_level, _node + stackTop);
            _node[iNode].parentIdx = _node[stackTop].parentIdx;
            _node[iNode]._rootIdx = ROOTIDX;
            _node[stackTop].parentIdx = iNode;
        } else
            _node[iNode].add(_node + stackTop);

        prevTop = stackTop;
        stackTop = iNode;
        currentLevel = _node[stackTop].alpha;
        if (_node[stackTop].area == imgSize) // root _node found...done
            break;
    }
FLOOD_END:
    _node[stackTop].parentIdx = ROOTIDX;

    delete queue;
    Free(dimg);
    Free(isVisited);
    Free(isAvailable);
}

template <class Pixel> int AlphaTree<Pixel>::get_bitdepth(uint64_t num) {
    int ret = 0;
    while (num) {
        ret++;
        num >>= 1;
    }
    return ret;
}

template <class Pixel>
void AlphaTree<Pixel>::FloodHierarHeapQueueNoCache(const Pixel *img, double a, double r, int listsize) {
    HierarHeapQueue<Pixel> *queue;

    ImgIdx imgSize, dimgSize, nredges, x0;
    uint64_t numlevels, max_level;
    double currentLevel;
    ImgIdx *dhist;
    ImgIdx stackTop, iNode;
    uint8_t *isVisited, *isAvailable, isAv;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (1 + (_connectivity >> 1)) * _width * _height;

    max_level = (sizeof(Pixel) == 8) ? 0xffffffffffffffff : (int64_t)((Pixel)(-1));
    numlevels = (uint64_t)(a * log2(1 + (double)max_level)) + 1;

    Pixel *dimg;
    dhist = (ImgIdx *)Calloc((size_t)numlevels * sizeof(ImgIdx));
    dimg = (Pixel *)Malloc((size_t)dimgSize * sizeof(Pixel));

    compute_dimg(x0, currentLevel, dimg, dhist, img, a); // calculate pixel differences and make histogram

    // create hierarchical queue from dhist
    queue =
        new HierarHeapQueue<Pixel>(dhist, numlevels, nredges, a, listsize, _connectivity, r); // +1 for the dummy _node
    _curSize = 0;

    _maxSize = 1 + imgSize +
               dimgSize; // Do not use TSE here, becasue dhist is a logged histogram (also this algorithm is for hdr)

    dhist = 0;
    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable(isAvailable);
    _parentAry = (ImgIdx *)Malloc((size_t)imgSize * sizeof(int32_t));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));

    stackTop = NewAlphaNode(); /*dummy root*/
    AlphaNode<Pixel> *pNode = _node + stackTop;
    pNode->set(0, (double)max_level, (double)0.0, (Pixel)max_level, (Pixel)0);
    pNode->parentIdx = stackTop;
    currentLevel = (double)max_level;

    bool firstpixel = 1;
    p = x0;
    ImgIdx prevTop = stackTop;
    while (true) // flooding
    {
        while (firstpixel || (double)queue->top_alpha() <= (double)currentLevel) // flood all levels below currentLevel
        {
            if (!firstpixel) {
                p = queue->pop(isVisited);
                if (isVisited[p])
                    continue;
            } else
                firstpixel = 0;

            isVisited[p] = 1;
            isAv = isAvailable[p];
            if (_connectivity == 4) {
                q = p << 1;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - (_width << 1)]);
            } else if (_connectivity == 8) {
                ImgIdx width4 = _width << 2;
                q = p << 2;
                if (is_available(isAv, 0) && !isVisited[p + _width]) {
                    queue->push(p + _width, dimg[q]);
                };
                if (is_available(isAv, 1) && !isVisited[p + _width + 1]) {
                    queue->push(p + _width + 1, dimg[q + 1]);
                };
                if (is_available(isAv, 2) && !isVisited[p + 1]) {
                    queue->push(p + 1, dimg[q + 2]);
                };
                if (is_available(isAv, 3) && !isVisited[p - _width + 1]) {
                    queue->push(p - _width + 1, dimg[q + 3]);
                };
                if (is_available(isAv, 4) && !isVisited[p - _width]) {
                    queue->push(p - _width, dimg[q - width4]);
                };
                if (is_available(isAv, 5) && !isVisited[p - _width - 1]) {
                    queue->push(p - _width - 1, dimg[q - width4 - 3]);
                };
                if (is_available(isAv, 6) && !isVisited[p - 1]) {
                    queue->push(p - 1, dimg[q - 2]);
                };
                if (is_available(isAv, 7) && !isVisited[p + _width - 1]) {
                    queue->push(p + _width - 1, dimg[q + width4 - 1]);
                };
            } else {
                //?
            }

            if ((double)currentLevel > (double)queue->top_alpha()) // go to lower level
            {
                Pixel pix_val = img[p];
                currentLevel = queue->top_alpha();
                iNode = NewAlphaNode();
                _node[iNode].set(1, currentLevel, (double)pix_val, pix_val, pix_val);
                _node[iNode].parentIdx = stackTop;
                _node[iNode]._rootIdx = ROOTIDX;
                prevTop = stackTop;
                stackTop = iNode;
                if (currentLevel > 0) {
                    iNode = NewAlphaNode(0, _node + stackTop);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                } else
                    _parentAry[p] = stackTop;
            } else {
                if (currentLevel > 0) {
                    Pixel pix_val = img[p];
                    iNode = NewAlphaNode();
                    _node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val);
                    _node[stackTop].add(_node + iNode);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                } else
                    connectPix2Node(p, img[p], stackTop);
                if (_node[stackTop].area == imgSize)
                    goto FLOOD_END;
            }
        }

        if (_node[prevTop].parentIdx == stackTop && _node[prevTop].area == _node[stackTop].area) {
            _node[prevTop].parentIdx = _node[stackTop].parentIdx;
            stackTop = prevTop;
            _curSize--;
        }

        if (_node[stackTop].area == imgSize) // root _node found...done
            break;

        // go to higher level
        ImgIdx parentNodeIdx = _node[stackTop].parentIdx;
        if ((double)queue->top_alpha() < (double)_node[parentNodeIdx].alpha) {
            parentNodeIdx = NewAlphaNode1(queue->top_alpha(), _node + stackTop);
            _node[parentNodeIdx].parentIdx = _node[stackTop].parentIdx;
            _node[parentNodeIdx]._rootIdx = ROOTIDX;
            _node[stackTop].parentIdx = parentNodeIdx;
        } else // go to existing _node
            _node[parentNodeIdx].add(_node + stackTop);

        prevTop = stackTop;
        stackTop = parentNodeIdx;
        currentLevel = _node[stackTop].alpha;
        if (_node[stackTop].area == imgSize) // root _node found...done
            break;
    }
FLOOD_END:
    _rootIdx = (_node[stackTop].area == imgSize) ? stackTop : iNode; // remove redundant root
    _node[_rootIdx].parentIdx = ROOTIDX;

    delete queue;
    Free(dimg);
    Free(isVisited);
    Free(isAvailable);
}

// hhpq
template <class Pixel> void AlphaTree<Pixel>::FloodHierarHeapQueue(const Pixel *img, double a, double r, int listsize) {
    assert(_connectivity == 4 || _connectivity == 8 || _connectivity == 12);
    clear();
    const ImgIdx imgSize = _width * _height;
    const ImgIdx nredges = _width * (_height - 1) + (_width - 1) * _height +
                           ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    const ImgIdx dimgSize = (1 + (_connectivity >> 1)) * _width * _height;
    const double alphaMax = _pixelDissim.maximumDissmilarity();
    const uint64_t numLevels = HHPQ::alphaToLevel((double)alphaMax, a) + 1;
    assert(numLevels < 10e3); // More than 10k levels is infeasible and inefficient

    ImgIdx *dhist = (ImgIdx *)Calloc((size_t)numLevels * sizeof(ImgIdx));
    double *dimg = (double *)Malloc((size_t)dimgSize * sizeof(double));

    compute_dimg_hhpq(dimg, dhist, img, a); // Calculate pixel differences and make histogram

    uint8_t *isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable(isAvailable);

    _curSize = 0;
    _maxSize = 1 + imgSize + dimgSize;
    _parentAry = (ImgIdx *)Malloc((size_t)imgSize * sizeof(ImgIdx));
    memset(_parentAry, ROOTIDX, imgSize * sizeof(ImgIdx));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));

    ImgIdx startingPixel = 0;
    runFloodHHPQ(startingPixel, img, a, r, listsize, imgSize, nredges, dimgSize, numLevels, dhist, dimg, isAvailable);

    Free(dimg);
    Free(dhist);
    Free(isAvailable);
}

template <class Pixel> void AlphaTree<Pixel>::printVisit(ImgIdx p, double q) const {
    printf("Visiting %d at %.2f\n", p, q);
}

template <class Pixel>
void AlphaTree<Pixel>::registerEdge(ImgIdx imgIdx, ImgIdx edgeIdx, ImgIdx *queuedEdges, uint8_t *numQueuedEdges) const {
    auto nqe = numQueuedEdges[imgIdx]++;
    queuedEdges[imgIdx * _connectivity + nqe] = edgeIdx;
}

template <class Pixel>
void AlphaTree<Pixel>::markRedundant(ImgIdx imgIdx, ImgIdx eIdx, uint8_t *edgeStatus, ImgIdx *queuedEdges,
                                     uint8_t *numQueuedEdges) const {
    auto nqe = numQueuedEdges[imgIdx];
    for (int i = 0; i < nqe; i++) {
        auto edgeIdx = queuedEdges[imgIdx * _connectivity + i];
        if (edgeIdx != eIdx)
            edgeStatus[edgeIdx] = QItem::EDGE_REDUNDANT;
    }
}

template <class Pixel>
void AlphaTree<Pixel>::runFloodHHPQ(ImgIdx startingPixel, const Pixel *img, double a, double r, int listsize,
                                    ImgIdx imgSize, ImgIdx nredges, ImgIdx dimgSize, uint64_t numLevels,
                                    const ImgIdx *dhist, const double *dimg, const uint8_t *isAvailable) {
    uint8_t *isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    HHPQ *queue = new HHPQ(dhist, numLevels, nredges, isVisited, a, listsize, r);

    ImgIdx stackTop = _curSize++; // Dummy root with maximum possible alpha
    double currentLevel = std::numeric_limits<double>::infinity();
    _node[stackTop] = AlphaNode<Pixel>(currentLevel);

    ImgIdx prevTop = stackTop;
    queue->push(startingPixel, currentLevel);

    while (_node[stackTop].area < imgSize) { // Main flooding loop
        while (_node[stackTop].area < imgSize && !queue->empty() &&
               queue->front().alpha <= currentLevel) { // Flood all levels below currentLevel
            const ImgIdx p = queue->front().index;
            currentLevel = queue->front().alpha;
            queue->pop();
            if (isVisited[p]) {

#if RUN_FANCY_VIS_FOR_DEBUG
                // printVisit(p, currentLevel);
                // queue->print();
                // printAll(isVisited, dimg);
                // std::getchar();
#endif

                continue;
            }

            isVisited[p] = 1;

            // printVisit(p, currentLevel);
            // printTree();
            // queue->print();
            // printAll(isVisited, edgeStatus, img);
            // std::getchar();

            auto isAv = isAvailable[p];
            if (_connectivity == 4) {
                const ImgIdx q = p << 1;
                const ImgIdx width2 = _width * 2;
                // clang-format off
                if (is_available(isAv, 0) && !isVisited[p + _width])    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])         queue->push(p + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])         queue->push(p - 1, dimg[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])    queue->push(p - _width, dimg[q - width2]);
                // clang-format on
            } else if (_connectivity == 8) {
                const ImgIdx width4 = _width << 2;
                const ImgIdx q = p << 2;
                // clang-format off
                if (is_available(isAv, 0) && !isVisited[p + _width])        queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])    queue->push(p + _width + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p + 1])             queue->push(p + 1, dimg[q + 2]);
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])    queue->push(p - _width + 1, dimg[q + 3]);
                if (is_available(isAv, 4) && !isVisited[p - _width])        queue->push(p - _width, dimg[q - width4]);
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])    queue->push(p - _width - 1, dimg[q - width4 - 3]);
                if (is_available(isAv, 6) && !isVisited[p - 1])             queue->push(p - 1, dimg[q - 2]);
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])    queue->push(p + _width - 1, dimg[q + width4 - 1]);
                // clang-format on
            } else {
                //?
            }

            if (currentLevel > queue->front().alpha) // go to lower level
            {
                currentLevel = queue->front().alpha;
                const ImgIdx newNodeIdx = _curSize++;
                _node[newNodeIdx] = AlphaNode<Pixel>(img[p], queue->front().alpha, stackTop);
                prevTop = stackTop;
                stackTop = newNodeIdx;

                if (currentLevel > 0) {
                    const ImgIdx singletonNodeIdx = _curSize++;
                    _node[singletonNodeIdx] = AlphaNode<Pixel>(img[p], 0.0, newNodeIdx);
                    _parentAry[p] = singletonNodeIdx;
                } else
                    _parentAry[p] = stackTop;
            } else {
                if (currentLevel > 0) {
                    const ImgIdx singletonNodeIdx = _curSize++;
                    _node[singletonNodeIdx] = AlphaNode<Pixel>(img[p], 0.0, stackTop);
                    _node[stackTop].add(_node[singletonNodeIdx]);
                    _parentAry[p] = singletonNodeIdx;
                } else
                    connectPix2Node(p, img[p], stackTop);
                if (_node[stackTop].area == imgSize)
                    break;
            }
        }

        if (_node[stackTop].area < imgSize && _node[prevTop].parentIdx == stackTop &&
            _node[prevTop].area == _node[stackTop].area) {
            _node[prevTop].parentIdx = _node[stackTop].parentIdx;
            stackTop = prevTop;
            _curSize--;
        }

        // go to higher level
        if (_node[stackTop].area < imgSize) {
            ImgIdx newParentIdx = _node[stackTop].parentIdx;
            if (!queue->empty() && (double)queue->front().alpha < (double)_node[newParentIdx].alpha) {
                newParentIdx = _curSize++;
                _node[newParentIdx] = AlphaNode<Pixel>(_node[stackTop]);
                _node[newParentIdx].alpha = queue->front().alpha;

                _node[newParentIdx].parentIdx = _node[stackTop].parentIdx;
                _node[stackTop].parentIdx = newParentIdx;
            } else // go to existing _node
                _node[newParentIdx].add(_node[stackTop]);

            prevTop = stackTop;
            stackTop = newParentIdx;
            currentLevel = _node[stackTop].alpha;
        }

        // queue->print();
        // printTree();
        // printAll(isVisited, edgeStatus, img);
        // std::getchar();
    }
    _rootIdx = _node[prevTop].area == imgSize ? prevTop : stackTop;
    _node[_rootIdx].parentIdx = ROOTIDX;

    // sortAlphaNodes();
    // _curSize--; // Remove dummy root

    // printParentAry();
    // printTree();
    delete queue;
    Free(isVisited);
}

template <class Pixel> void AlphaTree<Pixel>::sortAlphaNodes() {
    // Step 1: Capture the original indices
    std::vector<ImgIdx> original_indices(_curSize);
    for (ImgIdx i = 0; i < _curSize; ++i) {
        original_indices[i] = i;
    }

    AlphaNode<Pixel> *node = _node;
    // Step 2: Sort the array and keep track of the new indices
    std::sort(original_indices.begin(), original_indices.end(), [&node](ImgIdx a, ImgIdx b) {
        return node[a] < node[b]; // Sort based on the AlphaNode's alpha values
    });

    // Create a new sorted array
    std::vector<AlphaNode<Pixel>> sorted_nodes(_curSize);
    for (ImgIdx i = 0; i < _curSize; ++i) {
        sorted_nodes[i] = _node[original_indices[i]];
    }

    // Step 3: Create a mapping from old index to new index
    std::vector<ImgIdx> index_map(_curSize);
    for (ImgIdx i = 0; i < _curSize; ++i) {
        index_map[original_indices[i]] = i;
    }

    // Step 4: Update the parentIdx to reflect the new positions
    for (size_t i = 0; i < sorted_nodes.size(); ++i) {
        if (sorted_nodes[i].parentIdx != ROOTIDX) {
            sorted_nodes[i].parentIdx = index_map[sorted_nodes[i].parentIdx];
        }
    }

    for (int i = 0; i < _height * _width; ++i) {
        if (_parentAry[i] != ROOTIDX)
            _parentAry[i] = index_map[_parentAry[i]];
    }

    // Replace the original array with the sorted one
    for (size_t i = 0; i < sorted_nodes.size(); i++)
        _node[i] = sorted_nodes[i];
}

template <class Pixel> void AlphaTree<Pixel>::FloodHierHeapQueueHisteq(const Pixel *img, int listsize, int a) {
    HierarHeapQueue_HEQ<Pixel> *queue;

    ImgIdx imgSize, dimgSize, nredges, x0;
    uint64_t numlevels, max_level, currentLevel;
    ImgIdx *dhist;
    ImgIdx stackTop, iNode;
    uint8_t *isVisited, *isAvailable, isAv;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;

    // double a = 4.0;
    max_level = (sizeof(Pixel) == 8) ? 0xffffffffffffffff : (int64_t)((Pixel)(-1));
    int bitdepth = get_bitdepth(max_level);

    if (a == 0)
        a = 1024;
    int eqhistsize = a;
    double coeff = (double)eqhistsize * 10 / (double)bitdepth;

    numlevels = (uint64_t)(log2(1 + (double)max_level) * coeff) + 1;

    Pixel *dimg;
    dhist = (ImgIdx *)Calloc((size_t)numlevels * sizeof(ImgIdx));
    dimg = (Pixel *)Malloc((size_t)dimgSize * sizeof(Pixel));
    compute_dimg(dimg, dhist, img, coeff); // calculate pixel differences and make histogram
    ImgIdx *eqhist = (ImgIdx *)Calloc(eqhistsize * sizeof(ImgIdx));

    uint32_t *histeqmap = (uint32_t *)Malloc(numlevels * sizeof(uint32_t));

    cumsum(dhist, numlevels, histeqmap, eqhistsize);

    for (int ii = 0; ii < (int)numlevels; ii++) {
        int jj = histeqmap[ii];
        eqhist[jj] += dhist[ii];
    }

    queue = new HierarHeapQueue_HEQ<Pixel>(eqhist, histeqmap, eqhistsize, nredges, coeff,
                                           listsize); // +1 for the dummy _node
    _curSize = 0;

    _maxSize = 1 + imgSize +
               dimgSize; // Do not use TSE here, becasue dhist is a logged histogram (also this algorithm is for hdr)

    if (dhist)
        Free(dhist);
    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable(isAvailable);
    _parentAry = (ImgIdx *)Malloc((size_t)imgSize * sizeof(int32_t));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));

    stackTop = NewAlphaNode(); /*dummy root*/
    AlphaNode<Pixel> *pNode = _node + stackTop;
    pNode->set(0, (Pixel)max_level, (double)0.0, (Pixel)max_level, (Pixel)0);
    pNode->parentIdx = stackTop;
    currentLevel = max_level;
    x0 = 0; /*arbitrary starting point*/

    queue->push_1stitem(x0, currentLevel);
    while (1) // flooding
    {
        while ((uint64_t)queue->top_alpha() <= (uint64_t)currentLevel) // flood all levels below currentLevel
        {
            p = queue->top();

            if (isVisited[p]) {
                queue->pop(isVisited);
                continue;
            }

            queue->startPushes();
            isVisited[p] = 1;
            isAv = isAvailable[p];
            if (_connectivity == 4) {
                q = p << 1;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - (_width << 1)]);
            } else if (_connectivity == 8) {
                ImgIdx width4 = _width << 2;
                q = p << 2;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])
                    queue->push(p + _width + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p + 1])
                    queue->push(p + 1, dimg[q + 2]);
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])
                    queue->push(p - _width + 1, dimg[q + 3]);
                if (is_available(isAv, 4) && !isVisited[p - _width])
                    queue->push(p - _width, dimg[q - width4]);
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])
                    queue->push(p - _width - 1, dimg[q - width4 - 3]);
                if (is_available(isAv, 6) && !isVisited[p - 1])
                    queue->push(p - 1, dimg[q - 2]);
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])
                    queue->push(p + _width - 1, dimg[q + width4 - 1]);
            } else {
                //?
            }

            queue->endPushes(isVisited);
            if ((uint64_t)currentLevel > (uint64_t)queue->top_alpha()) // go to lower level
            {
                Pixel pix_val = img[p];
                currentLevel = queue->top_alpha();
                iNode = NewAlphaNode();
                _node[iNode].set(1, currentLevel, (double)pix_val, pix_val, pix_val);
                _node[iNode].parentIdx = stackTop;
                _node[iNode]._rootIdx = ROOTIDX;
                stackTop = iNode;
                if (currentLevel) {
                    iNode = NewAlphaNode(0, _node + stackTop);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                    // prevTop = iNode;
                } else
                    _parentAry[p] = stackTop;
            } else {
                if (currentLevel) {
                    Pixel pix_val = img[p];
                    iNode = NewAlphaNode();
                    _node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val);
                    _node[stackTop].add(_node + iNode);
                    _node[iNode].parentIdx = stackTop;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _parentAry[p] = iNode;
                } else
                    connectPix2Node(p, img[p], stackTop);
                if (_node[stackTop].area == imgSize)
                    goto FLOOD_END;
            }
        }
        if (_node[stackTop].area == imgSize) // root _node found...done
            break;

        // go to higher level
        iNode = _node[stackTop].parentIdx;
        if ((double)queue->top_alpha() < (double)_node[iNode].alpha) {
            iNode = NewAlphaNode1(queue->top_alpha(), _node + stackTop);
            _node[iNode].parentIdx = _node[stackTop].parentIdx;
            _node[iNode]._rootIdx = ROOTIDX;
            _node[stackTop].parentIdx = iNode;
        } else // go to existing _node
        {
            _node[iNode].add(_node + stackTop);
        }

        stackTop = iNode;
        currentLevel = (uint64_t)_node[stackTop].alpha;
        if (_node[stackTop].area == imgSize) // root _node found...done
            break;
    }
FLOOD_END:
    _rootIdx = (_node[stackTop].area == imgSize) ? stackTop : iNode; // remove redundant root
    _node[_rootIdx].parentIdx = ROOTIDX;

    delete queue;
    Free(dimg);
    Free(isVisited);
    Free(isAvailable);
}

template <class Pixel> ImgIdx AlphaTree<Pixel>::initialize_node(const Pixel *img, Pixel *dimg, Pixel maxpixval) {
    ImgIdx p, imgSize = _width * _height;
    ImgIdx maxdiffidx = 0;
    Pixel maxdiffval = 0;

    for (p = 0; p < _maxSize; p++) {
        if (p < imgSize)
            _node[p].set(1, 0, (double)img[p], img[p], img[p]);
        else {
            ImgIdx q = p - imgSize;
            if (maxdiffval < dimg[q]) {
                maxdiffval = dimg[q];
                maxdiffidx = q;
            }
            _node[p].set(0, dimg[q], 0.0, maxpixval, 0);
        }
    }

    return maxdiffidx;
}

template <class Pixel>
void AlphaTree<Pixel>::initialize_node1(const Pixel *img, RankItem<double> *rankitem, Pixel maxpixval) {
    ImgIdx p, imgSize = _width * _height;

    for (p = 0; p < _maxSize; p++) {
        if (p < imgSize)
            _node[p].set(1, 0, (double)img[p], img[p], img[p]);
        else
            _node[p].set(0, rankitem[p - imgSize].alpha, 0.0, maxpixval, 0);
        _node[p]._rootIdx = _node[p].parentIdx = ROOTIDX;
    }
}

template <class Pixel>
void AlphaTree<Pixel>::initialize_node1(const Pixel *img, RankItem<double> *rankitem, Pixel maxpixval,
                                        int32_t *rankToIndex) {
    ImgIdx p, imgSize = _width * _height;

    for (p = 0; p < _maxSize; p++) {
        ImgIdx q = p;
        if (p < imgSize)
            _node[p].set(1, 0, (double)img[p], img[p], img[p]);
        else {
            q = rankToIndex[p - imgSize];
            _node[p].set(0, rankitem[q].alpha, 0.0, maxpixval, 0);
            q = q + imgSize;
        }
        _node[q]._rootIdx = _node[q].parentIdx = ROOTIDX;
    }
}

template <class Pixel>
void AlphaTree<Pixel>::initialize_node(const Pixel *img, RankItem<Pixel> *rankitem, Pixel maxpixval) {
    ImgIdx p, imgSize = _width * _height;

    for (p = 0; p < _maxSize; p++) {
        if (p < imgSize)
            _node[p].set(1, 0, (double)img[p], img[p], img[p]);
        else
            _node[p].set(0, rankitem[p - imgSize].alpha, 0.0, maxpixval, 0);
        _node[p]._rootIdx = _node[p].parentIdx = ROOTIDX;
    }
}

template <class Pixel>
void AlphaTree<Pixel>::initialize_node_par(const Pixel *img, RankItem<Pixel> *rankitem, Pixel maxpixval) {
    ImgIdx p, imgSize = _width * _height;

#pragma omp parallel for schedule(guided, 1)
    for (p = 0; p < _maxSize; p++) {
        if (p < imgSize) {
            _node[p].set(1, 0, (double)img[p], img[p], img[p]);
            _node[p].parentIdx = _node[p]._rootIdx = ROOTIDX;
        } else {
            _node[p].set(0, rankitem[p - imgSize].alpha, 0.0, maxpixval, 0);
            _node[p].parentIdx = _node[p]._rootIdx = ROOTIDX;
        }
        // _node[p].print(_node);

        // _node[p].thread = -1;
    }
}

template <class Pixel>
void AlphaTree<Pixel>::initialize_node_par1(const Pixel *img, RankItem<double> *rankitem, Pixel maxpixval,
                                            int32_t *rankToIndex) {
    ImgIdx p, imgSize = _width * _height;

#pragma omp parallel for schedule(guided, 1)
    for (p = 0; p < _maxSize; p++) {
        if (p < imgSize) {
            _node[p].set(1, 0, (double)img[p], img[p], img[p]);
            _node[p].parentIdx = _node[p]._rootIdx = ROOTIDX;
        } else {
            if (rankToIndex)
                _node[p].set(0, rankitem[rankToIndex[p - imgSize]].alpha, 0.0, maxpixval, 0);
            else
                _node[p].set(0, rankitem[p - imgSize].alpha, 0.0, maxpixval, 0);
            _node[p].parentIdx = _node[p]._rootIdx = ROOTIDX;
        }
    }
}

template <class Pixel> void AlphaTree<Pixel>::init_hypergraph_nodes(Pixel *dimg) {
    if (_connectivity == 4) {
        ImgIdx imgSize = _height * _width;
        ImgIdx p = 0;
        ImgIdx wstride = _width << 1;
        for (ImgIdx y = 0; y < _height; y++) {
            for (ImgIdx x = 0; x < _width; x++) {
                ImgIdx q = p << 1;
                Pixel minalpha = (Pixel)(-1), minneighidx = 0;

                if ((y < _height - 1) && (dimg[q] < minalpha)) {
                    minalpha = dimg[q];
                    minneighidx = q;
                }
                if ((x < _width - 1) && (dimg[q + 1] < minalpha)) {
                    minalpha = dimg[q + 1];
                    minneighidx = q + 1;
                }
                if ((x > 0) && (dimg[q - 1] < minalpha)) {
                    minalpha = dimg[q - 1];
                    minneighidx = q - 1;
                }
                if ((y > 0) && (dimg[q - wstride] < minalpha)) {
                    minalpha = dimg[q - wstride];
                    minneighidx = q - wstride;
                }

                _node[p++].connect_to_parent(&_nodeIn[minneighidx], minneighidx + imgSize);
            }
        }
    } else {
    }
}

// connect pixels to one of incident edges with the minimum edge weight, so that
// pixels are no longer needed to be inspected to make the min-tree implementation
// easier
template <class Pixel> void AlphaTree<Pixel>::init_hypergraph_nodes(ImgIdx *rank) {
    if (_connectivity == 4) {
        ImgIdx imgSize = _height * _width;
        for (ImgIdx p = 0; p < imgSize; p++) {
            ImgIdx q = p << 1;
            ImgIdx y = p / _width;
            ImgIdx x = p % _width;
            ImgIdx minRank = 2 * imgSize;

            ((y < _height - 1) && (rank[q] < minRank)) ? (minRank = rank[q]) : (ImgIdx)0;
            ((x < _width - 1) && (rank[q + 1] < minRank)) ? (minRank = rank[q + 1]) : (ImgIdx)0;
            ((x > 0) && (rank[q - 1] < minRank)) ? (minRank = rank[q - 1]) : (ImgIdx)0;
            ((y > 0) && (rank[q - (_width << 1)] < minRank)) ? (minRank = rank[q - (_width << 1)]) : (ImgIdx)0;

            _node[p].connect_to_parent(&_nodeIn[minRank], minRank + imgSize);
        }
    } else {
    }
}

template <class Pixel> void AlphaTree<Pixel>::set_isAvailable_hypergraph(uint8_t *isAvailable) {
    if (_connectivity == 4) {
        ImgIdx dimgidx;
        ImgIdx width2 = 2 * _width;
        ImgIdx dimgSize = _width * _height * 2;

        // first row
        for (dimgidx = 0; dimgidx < width2;) {
            isAvailable[dimgidx++] |= 0x20; // even neighbor 5
            isAvailable[dimgidx++] |= 0x30; // odd neighbor 4, 5
        }

        for (dimgidx -= 3; dimgidx < dimgSize; dimgidx += width2) {
            isAvailable[dimgidx] |= 0x02;     // odd neighbor 1
            isAvailable[dimgidx + 1] |= 0x09; // even neighbor 0, 3
        }

        for (dimgidx = 0; dimgidx < dimgSize; dimgidx += width2) {
            isAvailable[dimgidx] |= 0x12;     // even neighbor 1, 4
            isAvailable[dimgidx + 1] |= 0x08; // odd neighbor 3
        }

        // last 2 rows
        for (dimgidx = _width * (_height - 2) * 2; dimgidx < dimgSize - width2; dimgidx += 2)
            isAvailable[dimgidx] |= 0x04; // even neighbor 2
        for (dimgidx += 1; dimgidx < dimgSize; dimgidx += 2)
            isAvailable[dimgidx] |= 0x05; // odd neighbor 0, 2

    } else // later
    {
    }
}

template <class Pixel>
uint8_t AlphaTree<Pixel>::push_neighbor(Trie<TrieIdx> *queue, uint8_t *isVisited, ImgIdx *rank, ImgIdx p) {
    isVisited[p] = 1;
    return queue->push(rank[p]);
}

template <class Pixel> void AlphaTree<Pixel>::FloodTrieHypergraph(const Pixel *img) {
    ImgIdx imgSize, dimgSize, nredges;
    ImgIdx current_rank = 0, next_rank = 0;
    RankItem<double> *rankitem, *pRank;
    Pixel maxpixval;
    ImgIdx *rank;
    int8_t nbits;
    uint8_t *isVisited, *isAvailable, isAv;
    ImgIdx p;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    _maxSize = imgSize + nredges;
    num_node = _maxSize;
    num_node_in = nredges;
    nbits = ((sizeof(Pixel) << 3) - 1);
    maxpixval = ~(1 << nbits);
    rankitem = (RankItem<double> *)Malloc(nredges * sizeof(RankItem<double>));
    _parentAry = 0;
    rank = (ImgIdx *)Malloc((size_t)dimgSize * sizeof(ImgIdx));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));
    _nodeIn = _node + imgSize;

    ImgIdx wstride_d = _width << 1;

    Trie<TrieIdx> *queue = new Trie<TrieIdx>(nredges);

    omp_set_num_threads(1);
    int32_t *rankToIndex = (int32_t *)Calloc(nredges * sizeof(int32_t));
    computeDifferenceAndSort(rank, rankitem, img, nredges, rankToIndex);
    initialize_node1(img, rankitem, maxpixval, rankToIndex);

    init_hypergraph_nodes(rank);

    isVisited = (uint8_t *)Calloc(dimgSize * sizeof(uint8_t));
    isAvailable = (uint8_t *)Calloc(dimgSize * sizeof(uint8_t));

    set_isAvailable_hypergraph(isAvailable);

    current_rank = nredges - 1;
    queue->push(nredges - 1);
    isVisited[rank[nredges] - 1] = 1;

    while (1) {
        while (1) {
            current_rank = queue->top();

            pRank = rankitem + rankToIndex[current_rank];
            p = pRank->dimgidx;

            uint8_t gotolowerlevel = 0;
            if (_connectivity == 4) {
                isAv = ~isAvailable[p];

                if (p & 1)
                    gotolowerlevel =
                        ((is_available(isAv, 0) && !isVisited[p + 1]) ? push_neighbor(queue, isVisited, rank, p + 1)
                                                                      : 0) ||
                        ((is_available(isAv, 1) && !isVisited[p + 2]) ? push_neighbor(queue, isVisited, rank, p + 2)
                                                                      : 0) ||
                        ((is_available(isAv, 2) && !isVisited[p - 1]) ? push_neighbor(queue, isVisited, rank, p - 1)
                                                                      : 0) ||
                        ((is_available(isAv, 3) && !isVisited[p - 2]) ? push_neighbor(queue, isVisited, rank, p - 2)
                                                                      : 0) ||
                        ((is_available(isAv, 4) && !isVisited[p - wstride_d - 1])
                             ? push_neighbor(queue, isVisited, rank, p - wstride_d - 1)
                             : 0) ||
                        ((is_available(isAv, 5) && !isVisited[p - wstride_d + 1])
                             ? push_neighbor(queue, isVisited, rank, p - wstride_d + 1)
                             : 0);
                else
                    gotolowerlevel =
                        ((is_available(isAv, 0) && !isVisited[p + 1]) ? push_neighbor(queue, isVisited, rank, p + 1)
                                                                      : 0) ||
                        ((is_available(isAv, 1) && !isVisited[p - 1]) ? push_neighbor(queue, isVisited, rank, p - 1)
                                                                      : 0) ||
                        ((is_available(isAv, 2) && !isVisited[p + wstride_d])
                             ? push_neighbor(queue, isVisited, rank, p + wstride_d)
                             : 0) ||
                        ((is_available(isAv, 3) && !isVisited[p + wstride_d + 1])
                             ? push_neighbor(queue, isVisited, rank, p + wstride_d + 1)
                             : 0) ||
                        ((is_available(isAv, 4) && !isVisited[p + wstride_d - 1])
                             ? push_neighbor(queue, isVisited, rank, p + wstride_d - 1)
                             : 0) ||
                        ((is_available(isAv, 5) && !isVisited[p - wstride_d])
                             ? push_neighbor(queue, isVisited, rank, p - wstride_d)
                             : 0);
            } else {
            }

            if (!gotolowerlevel)
                break;
        }

        queue->pop();
        next_rank = queue->top();

        _nodeIn[current_rank].connect_to_parent(&_nodeIn[next_rank], next_rank + imgSize);
        if (_nodeIn[next_rank].area == imgSize)
            break;

        current_rank = next_rank;
    }

    _rootIdx = (_nodeIn[current_rank].area == imgSize) ? current_rank + imgSize : next_rank + imgSize;
    _node[_rootIdx].parentIdx = ROOTIDX;

    delete queue;
    Free(rankToIndex);
    Free(rank);
    Free(rankitem);
    Free(isVisited);
    Free(isAvailable);
}

/*
6 - neighbor arrangement in isAvailable array
(order of neighbors is designed to minimize cache miss)
x - image edge (hypergraph _node)
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
template <class Pixel>
void AlphaTree<Pixel>::set_isAvailable_par_hypergraph(uint8_t *isAvailable, int8_t npartition_x, int8_t npartition_y) {
    ImgIdx wstride = _width / npartition_x * 2;
    ImgIdx wres2 = (_width % npartition_x) * 2;
    ImgIdx hstride = _height / npartition_y;
    ImgIdx hres = _height % npartition_y;

    if (_connectivity == 4) {
        ImgIdx dimgidx;
        ImgIdx width2 = 2 * _width;
        ImgIdx dimgSize = _width * _height * 2;

        // subimage first rows
        for (int y = 0; y < npartition_y; y++) {
            ImgIdx nextrowidx = (y * hstride + 1) * width2;
            for (dimgidx = (y * hstride) * width2; dimgidx < nextrowidx;) {
                isAvailable[dimgidx++] |= 0x20; // even neighbor 5
                isAvailable[dimgidx++] |= 0x30; // odd neighbor 4, 5
            }
        }

        // subimage right edges
        for (int x = 0; x < npartition_x; x++) {
            if (x == 0) {
                dimgidx = width2 - 3;
                for (; dimgidx < dimgSize; dimgidx += width2) {
                    isAvailable[dimgidx] |= 0x02;     // odd neighbor 1
                    isAvailable[dimgidx + 1] |= 0x09; // even neighbor 0, 3
                }
            } else {
                // edges on subblock borders belong to subblocks on the left
                dimgidx = width2 - 1 - x * wstride - wres2;
                for (; dimgidx < dimgSize; dimgidx += width2)
                    isAvailable[dimgidx] |= 0x23; // odd neighbor 0, 1, 5
            }
        }

        // subimage left edges
        for (int x = 0; x < npartition_x; x++) {
            for (dimgidx = x * wstride; dimgidx < dimgSize; dimgidx += width2) {
                isAvailable[dimgidx] |= 0x12;     // even neighbor 1, 4
                isAvailable[dimgidx + 1] |= 0x08; // odd neighbor 3
            }
        }

        // last 2 rows
        for (int y = 0; y < npartition_y; y++) {
            if (y == 0) {

                ImgIdx subimgend = (_height - 1) * width2;
                for (dimgidx = (_height - 2) * width2; dimgidx < subimgend; dimgidx += 2)
                    isAvailable[dimgidx] |= 0x04; // even neighbor 2
                subimgend = _height * width2;
                for (dimgidx = (_height - 1) * width2 + 1; dimgidx < subimgend; dimgidx += 2)
                    isAvailable[dimgidx] |= 0x05; // odd neighbor 0, 2
            } else {
                dimgidx = (_height - 1 - y * hstride - hres) * width2;
                ImgIdx subimgend = dimgidx + width2;
                for (; dimgidx < subimgend; dimgidx += 2)
                    isAvailable[dimgidx] |= 0x1c; // even neighbor 2, 3, 4
            }
        }

    } else {
    }
}

template <class Pixel> void AlphaTree<Pixel>::cumsum(ImgIdx *hist, ImgIdx size, ImgIdx &maxidx) {
    ImgIdx sum = hist[0];
    maxidx = 0;
    for (ImgIdx i = 1; i < size; i++) {
        if (hist[i]) {
            maxidx = i;
            sum += hist[i];
        }
        hist[i] = sum;
    }
}

template <class Pixel> void AlphaTree<Pixel>::cumsum(ImgIdx *hist, ImgIdx size, uint32_t *histeqmap, int eqhistsize) {
    ImgIdx sum = hist[0];
    for (ImgIdx i = 1; i < size; i++) {
        sum += hist[i];
    }

    double coeff = (double)(eqhistsize - 1) / (double)sum;
    sum = 0;
    for (ImgIdx i = 0; i < size; i++) {
        sum += hist[i];
        histeqmap[i] = (uint16_t)((double)sum * coeff);
    }
}

template <class Pixel>
uint8_t AlphaTree<Pixel>::push_neighbor(HierarQueue *queue, uint8_t *isVisited, uint8_t *dimg, ImgIdx p) {
    isVisited[p] = 1;
    return queue->push(p, dimg[p]);
}

template <class Pixel>
uint8_t AlphaTree<Pixel>::push_neighbor(HierarQueue *queue, uint8_t *isVisited, uint16_t *dimg, ImgIdx p) {
    isVisited[p] = 1;
    return queue->push(p, dimg[p]);
}

template <class Pixel>
uint8_t AlphaTree<Pixel>::push_neighbor(HierarQueue *queue, uint8_t *isVisited, uint32_t *dimg, ImgIdx p) {
    isVisited[p] = 1;
    return queue->push(p, dimg[p]);
}

template <class Pixel>
uint8_t AlphaTree<Pixel>::push_neighbor(HierarQueue *queue, uint8_t *isVisited, uint64_t *dimg, ImgIdx p) {
    isVisited[p] = 1;
    return queue->push(p, dimg[p]);
}

template <class Pixel> void AlphaTree<Pixel>::FloodHierarQueueHypergraph(const Pixel *img) {
    if (sizeof(Pixel) > 2 || _channel > 1) {
        printf("Error: Hierarchical queues do not work on >16 bits images or multispectral images\n");
        printf("Try Unionfind (algorithm code %d), flooding using Heapqueue (%d), trie queue (%d) or cached trie queue "
               "(%d) \n",
               UNIONFIND, FLOOD_HEAPQUEUE_CACHE, FLOOD_TRIE, FLOOD_TRIE_CACHE);
        return;
    }

    HierarQueue *queue;

    ImgIdx imgSize, dimgSize, nredges;
    int64_t numlevels, max_level, currentLevel;
    ImgIdx *dhist;
    Pixel *dimg;
    ImgIdx stackTop, iNode;
    uint8_t *isVisited, *isAvailable;
    ImgIdx p, wstride_d = _width << 1;

    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    max_level = (int64_t)((Pixel)(-1));
    numlevels = max_level + 1;

    dhist = (ImgIdx *)Malloc((size_t)numlevels * sizeof(ImgIdx));
    memset(dhist, 0, (size_t)numlevels * sizeof(int32_t));
    dimg = (Pixel *)Calloc((size_t)dimgSize * sizeof(Pixel));

    compute_dimg(dimg, dhist, img); // calculate pixel differences and make histogram

    _maxSize = imgSize + dimgSize;
    _node = (AlphaNode<Pixel> *)Calloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));
    _nodeIn = _node + imgSize;

    ImgIdx maxdiffidx = initialize_node(img, dimg, max_level);
    max_level = dimg[maxdiffidx];

    init_hypergraph_nodes(dimg);

    // create hierarchical queue from dhist
    queue = new HierarQueue(nredges, dhist, numlevels);
    _curSize = 0;

    Free(dhist);
    isVisited = (uint8_t *)Calloc(dimgSize * sizeof(uint8_t));
    isAvailable = (uint8_t *)Calloc(dimgSize * sizeof(uint8_t));

    omp_set_num_threads(1);
    int8_t npartition_x = 1, npartition_y = 1;
    set_isAvailable_par_hypergraph(isAvailable, npartition_x, npartition_y);

    stackTop = maxdiffidx + imgSize; // root
    AlphaNode<Pixel> *pNode = _node + stackTop;
    pNode->parentIdx = stackTop;
    currentLevel = max_level;
    queue->push(maxdiffidx, currentLevel);
    isVisited[maxdiffidx] = 1;
    while (1) // flooding
    {
        while (1) // flood all levels below currentLevel
        {
            p = queue->top();

            uint8_t gotolowerlevel = 0;
            if (_connectivity == 4) {
                uint8_t isAv = ~isAvailable[p];

                if (p & 1)
                    gotolowerlevel =
                        ((is_available(isAv, 0) && !isVisited[p + 1]) ? push_neighbor(queue, isVisited, dimg, p + 1)
                                                                      : 0) ||
                        ((is_available(isAv, 1) && !isVisited[p + 2]) ? push_neighbor(queue, isVisited, dimg, p + 2)
                                                                      : 0) ||
                        ((is_available(isAv, 2) && !isVisited[p - 1]) ? push_neighbor(queue, isVisited, dimg, p - 1)
                                                                      : 0) ||
                        ((is_available(isAv, 3) && !isVisited[p - 2]) ? push_neighbor(queue, isVisited, dimg, p - 2)
                                                                      : 0) ||
                        ((is_available(isAv, 4) && !isVisited[p - wstride_d - 1])
                             ? push_neighbor(queue, isVisited, dimg, p - wstride_d - 1)
                             : 0) ||
                        ((is_available(isAv, 5) && !isVisited[p - wstride_d + 1])
                             ? push_neighbor(queue, isVisited, dimg, p - wstride_d + 1)
                             : 0);
                else
                    gotolowerlevel =
                        ((is_available(isAv, 0) && !isVisited[p + 1]) ? push_neighbor(queue, isVisited, dimg, p + 1)
                                                                      : 0) ||
                        ((is_available(isAv, 1) && !isVisited[p - 1]) ? push_neighbor(queue, isVisited, dimg, p - 1)
                                                                      : 0) ||
                        ((is_available(isAv, 2) && !isVisited[p + wstride_d])
                             ? push_neighbor(queue, isVisited, dimg, p + wstride_d)
                             : 0) ||
                        ((is_available(isAv, 3) && !isVisited[p + wstride_d + 1])
                             ? push_neighbor(queue, isVisited, dimg, p + wstride_d + 1)
                             : 0) ||
                        ((is_available(isAv, 4) && !isVisited[p + wstride_d - 1])
                             ? push_neighbor(queue, isVisited, dimg, p + wstride_d - 1)
                             : 0) ||
                        ((is_available(isAv, 5) && !isVisited[p - wstride_d])
                             ? push_neighbor(queue, isVisited, dimg, p - wstride_d)
                             : 0);
            } else {
                // LATER
            }

            if (gotolowerlevel) {
                currentLevel = queue->min_level;
                iNode = queue->top() + imgSize;
                _node[iNode].parentIdx = stackTop;
                stackTop = iNode;
            } else {
                queue->pop();
                queue->find_minlev();

                if (queue->min_level == currentLevel) {
                    iNode = queue->top() + imgSize;
                    _node[stackTop].add(_node + iNode);
                    _node[iNode].parentIdx = stackTop;
                } else
                    break;
            }
        }
        // go to higher level
        iNode = _node[stackTop].parentIdx;
        if ((int64_t)queue->min_level < (int64_t)_node[iNode].alpha) // new level from queue
        {
            iNode = queue->top() + imgSize;
            _node[iNode].add(_node + stackTop);
            _node[iNode].parentIdx = _node[stackTop].parentIdx;
            _node[stackTop].parentIdx = iNode;
        } else // go to existing _node
        {
            if (_node[iNode].area == imgSize) // root _node found...done
                break;
            _node[iNode].add(_node + stackTop);
        }

        if (_node[iNode].area == imgSize) // root _node found...done
            break;

        stackTop = iNode;
        currentLevel = _node[stackTop].alpha;
    }

    _rootIdx = _node[iNode].area == imgSize ? iNode : stackTop;
    _node[_rootIdx].parentIdx = ROOTIDX;

    delete queue;
    Free(dimg);
    Free(isVisited);
    Free(isAvailable);
}

template <class Pixel> void AlphaTree<Pixel>::canonicalize(ImgIdx nidx) {
    ImgIdx p, q;

    p = get_level_root(nidx); // for 0-ccs
    if (p != nidx)
        _node[nidx].parentIdx = p;

    while (1) {
        q = _node[p].parentIdx;
        if (q == ROOTIDX)
            break;
        q = get_level_root(q);
        _node[p].parentIdx = q;
        p = q;
    }
}

// returns root _node index
template <class Pixel>
ImgIdx AlphaTree<Pixel>::mergePartition(Pixel *dimg, int64_t blksz_x, int64_t blksz_y, ImgIdx npartition_x,
                                        ImgIdx npartition_y, ImgIdx *subtree_cur) {
    const int64_t blksz_x0 = blksz_x;
    const int64_t blksz_y0 = blksz_y;
    const ImgIdx npartition_x0 = npartition_x;
    const ImgIdx npartition_y0 = npartition_y;
    const ImgIdx blkrow = _width * blksz_y0;
    const int shamt = _connectivity / 4;
    const ImgIdx width4 = _width << 2;
    while (npartition_x > 1 || npartition_y > 1) {
        // merge horizontal borders
        if ((npartition_x == 1 || blksz_x >= blksz_y) && npartition_y > 1) {
            const ImgIdx numblk = npartition_x * (npartition_y / 2);

#pragma omp parallel for schedule(dynamic, 1)
            for (int blk = 0; blk < numblk; blk++) {
                const ImgIdx y = (1 + 2 * (blk / (int)npartition_x)) * blksz_y;
                const ImgIdx x = (blk % (int)npartition_x) * blksz_x;
                const ImgIdx p0 = (y - 1) * _width + x;
                const ImgIdx pn = (((blk % (int)npartition_x) == npartition_x - 1) ? y * _width : p0 + blksz_x);

                const ImgIdx by = _min((p0 / blkrow), npartition_y0 - 1);
                for (ImgIdx p = p0; p < pn; p++) {
                    const ImgIdx bx = _min(((p % _width) / blksz_x0), npartition_x0 - 1);
                    const ImgIdx bidx = by * npartition_x0 + bx;
                    const ImgIdx q = p << shamt;

                    if (_connectivity == 4) {
                        connect(p, p + _width, (ImgIdx)subtree_cur[bidx]++, (Pixel)dimg[q]);
                    } else if (_connectivity == 8) {
                        if (p > p0)
                            connect(p, p + _width - 1, (ImgIdx)subtree_cur[bidx]++, (Pixel)dimg[q + width4 - 1]);
                        connect(p, p + _width, (ImgIdx)subtree_cur[bidx]++, (Pixel)dimg[q]);
                        if (p < pn - 1)
                            connect(p, p + _width + 1, (ImgIdx)subtree_cur[bidx]++, (Pixel)dimg[q + 1]);
                    }
                }
            }
            npartition_y = (npartition_y + 1) / 2;
            blksz_y <<= 1;
            if (npartition_y == 1)
                blksz_y = _height;
            else
                blksz_y = _min(blksz_y, _height);
        }

        if ((npartition_y == 1 || blksz_x <= blksz_y) && npartition_x > 1) {
            const ImgIdx numblk = npartition_y * (npartition_x / 2);

#pragma omp parallel for schedule(dynamic, 1)
            for (int blk = 0; blk < numblk; blk++) {
                const ImgIdx x = (1 + 2 * (blk / npartition_y)) * blksz_x;
                const ImgIdx y = (blk % (int)npartition_y) * blksz_y;
                const ImgIdx p0 = y * _width + x - 1;
                const ImgIdx pn =
                    ((blk % (int)npartition_y) == npartition_y - 1) ? _height * _width : p0 + _width * blksz_y;
                const ImgIdx bx = _min(((p0 % _width) / blksz_x0), npartition_x0 - 1);
                for (ImgIdx p = p0; p < pn; p += _width) {
                    const ImgIdx by = _min((p / blkrow), npartition_y0 - 1);
                    const ImgIdx bidx = by * npartition_x0 + bx;
                    const ImgIdx q = p << shamt;

                    if (_connectivity == 4) {
                        connect(p, p + 1, (ImgIdx)subtree_cur[bidx]++, (Pixel)dimg[q + 1]);
                    } else if (_connectivity == 8) {
                        if (p > p0)
                            connect(p, p - _width + 1, (ImgIdx)subtree_cur[bidx]++, (Pixel)dimg[q + 3]);
                        connect(p, p + 1, (ImgIdx)subtree_cur[bidx]++, (Pixel)dimg[q + 2]);
                        if (p < pn - _width)
                            connect(p, p + _width + 1, (ImgIdx)subtree_cur[bidx]++, (Pixel)dimg[q + 1]);
                    }
                }
            }
            npartition_x = (npartition_x + 1) / 2;
            blksz_x <<= 1;
            if (npartition_x == 1)
                blksz_x = _width;
            else
                blksz_x = _min(blksz_x, _width);
        }
    }

    // Find root
    ImgIdx p = 0;
    while (_node[p].parentIdx != ROOTIDX)
        p = _node[p].parentIdx;
    _rootIdx = p;
    return _rootIdx;
}

template <class Pixel>
ImgIdx AlphaTree<Pixel>::merge_subtrees(uint8_t *dimg, int64_t blksz_x, int64_t blksz_y, int16_t npartition_x,
                                        int16_t npartition_y, ImgIdx *subtree_cur, int tse) {
    ImgIdx numblk;
    int64_t blksz_x0 = blksz_x;
    int64_t blksz_y0 = blksz_y;

    ImgIdx npartition_x0 = npartition_x;
    ImgIdx npartition_y0 = npartition_y;
    ImgIdx blkrow = _width * blksz_y0;
    while (npartition_x > 1 || npartition_y > 1) {
        if ((npartition_x == 1 || blksz_x >= blksz_y) && npartition_y > 1) {
            numblk = npartition_x * (npartition_y / 2);

#pragma omp parallel for
            for (int blk = 0; blk < numblk; blk++) {
                ImgIdx x, y, r, p, p0, pn, dimgidx, bx, by, bidx;
                y = (1 + 2 * (blk / (int)npartition_x)) * blksz_y;
                x = (blk % (int)npartition_x) * blksz_x;

                p0 = (y - 1) * _width + x;
                pn = (((blk % (int)npartition_x) == npartition_x - 1) ? y * _width : p0 + blksz_x);

                by = _min((p0 / blkrow), npartition_y0 - 1);
                for (p = p0; p < pn; p++) {
                    bx = _min(((p % _width) / blksz_x0), npartition_x0 - 1);
                    bidx = by * npartition_x0 + bx;

                    dimgidx = p << 1;
                    r = dimg[dimgidx];

                    if (tse)
                        connect(_parentAry[p], _parentAry[p + _width], (ImgIdx)subtree_cur[bidx]++, (Pixel)r);
                    else
                        connect(p, p + _width, (ImgIdx)subtree_cur[bidx]++, (Pixel)r);
                }
            }
            npartition_y = (npartition_y + 1) / 2;
            blksz_y <<= 1;
            if (npartition_y == 1)
                blksz_y = _height;
            else
                blksz_y = _min(blksz_y, _height);
        }

        // merge vertical borders
        if ((npartition_y == 1 || blksz_x <= blksz_y) && npartition_x > 1) {
            numblk = npartition_y * (npartition_x / 2);

#pragma omp parallel for
            for (int blk = 0; blk < numblk; blk++) {
                ImgIdx x, y, r, p, p0, pn, dimgidx, bx, by, bidx;
                ;
                x = (1 + 2 * (blk / npartition_y)) * blksz_x;
                y = (blk % (int)npartition_y) * blksz_y;

                p0 = y * _width + x - 1;
                pn = ((blk % (int)npartition_y) == npartition_y - 1) ? _height * _width : p0 + _width * blksz_y;

                bx = _min(((p0 % _width) / blksz_x0), npartition_x0 - 1);
                for (p = p0; p < pn; p += _width) {
                    by = _min((p / blkrow), npartition_y0 - 1);
                    bidx = by * npartition_x0 + bx;

                    dimgidx = (p << 1) + 1;
                    r = dimg[dimgidx];

                    if (tse)
                        connect(_parentAry[p], _parentAry[p + 1], (ImgIdx)subtree_cur[bidx]++, (Pixel)r);
                    else
                        connect(p, p + 1, (ImgIdx)subtree_cur[bidx]++, (Pixel)r);
                }
            }
            npartition_x = (npartition_x + 1) / 2;
            blksz_x <<= 1;
            if (npartition_x == 1)
                blksz_x = _width;
            else
                blksz_x = _min(blksz_x, _width);
        }
    }

    ImgIdx p;
    if (tse)
        p = _parentAry[0];
    else
        p = 0;
    while (_node[p].parentIdx != ROOTIDX)
        p = _node[p].parentIdx;

    return p;
}

template <class Pixel>
ImgIdx AlphaTree<Pixel>::mergeSubtreeHybrid(uint8_t *dimg, ImgIdx blksz_x, ImgIdx blksz_y, ImgIdx npartition_x,
                                            ImgIdx npartition_y, ImgIdx *subtree_cur) {
    const ImgIdx blksz_x0 = blksz_x;
    const ImgIdx blksz_y0 = blksz_y;
    const ImgIdx npartition_x0 = npartition_x;
    const ImgIdx npartition_y0 = npartition_y;
    const ImgIdx blkrow = _width * blksz_y0;
    while (npartition_x > 1 || npartition_y > 1) {
        if ((npartition_x == 1 || blksz_x >= blksz_y) && npartition_y > 1) {
            const ImgIdx numblk = npartition_x * (npartition_y / 2); // 8N

#pragma omp parallel for schedule(dynamic, 1)
            for (int blk = 0; blk < numblk; blk++) {
                const ImgIdx y = (1 + 2 * (blk / npartition_x)) * blksz_y; // 8N?
                const ImgIdx x = (blk % npartition_x) * blksz_x;
                const ImgIdx p0 = (y - 1) * _width + x;
                const ImgIdx pn = (((blk % npartition_x) == npartition_x - 1) ? y * _width : p0 + blksz_x);
                const ImgIdx by = _min((p0 / blkrow), npartition_y0 - 1);

                for (ImgIdx p = p0; p < pn; p++) {
                    const ImgIdx bx = _min(((p % _width) / blksz_x0), npartition_x0 - 1);
                    const ImgIdx bidx = by * npartition_x0 + bx;
                    const ImgIdx dimgidx = p << 1; // 8N

                    connect(_parentAry[p], _parentAry[p + _width], (ImgIdx)subtree_cur[bidx]++, (Pixel)dimg[dimgidx]);
                }
            }
            npartition_y = (npartition_y + 1) / 2; // 8N
            blksz_y <<= 1;
            blksz_y = (npartition_y == 1) ? _height : _min(blksz_y, _height);
        }

        if ((npartition_y == 1 || blksz_x <= blksz_y) && npartition_x > 1) {
            const ImgIdx numblk = npartition_y * (npartition_x / 2);

#pragma omp parallel for schedule(dynamic, 1)
            for (int blk = 0; blk < numblk; blk++) {
                const ImgIdx x = (1 + 2 * (blk / npartition_y)) * blksz_x;
                const ImgIdx y = (blk % (int)npartition_y) * blksz_y;
                const ImgIdx p0 = y * _width + x - 1;
                const ImgIdx pn =
                    ((blk % (int)npartition_y) == npartition_y - 1) ? _height * _width : p0 + _width * blksz_y;

                const ImgIdx bx = _min(((p0 % _width) / blksz_x0), npartition_x0 - 1);
                for (ImgIdx p = p0; p < pn; p += _width) {
                    const ImgIdx by = _min((p / blkrow), npartition_y0 - 1);
                    const ImgIdx bidx = by * npartition_x0 + bx;
                    const ImgIdx dimgidx = (p << 1) + 1; // 8N

                    connect(_parentAry[p], _parentAry[p + 1], (ImgIdx)subtree_cur[bidx]++, (Pixel)dimg[dimgidx]);
                    // connect 8N
                    // connect 8N
                }
            }
            npartition_x = (npartition_x + 1) / 2;
            blksz_x <<= 1;
            blksz_x = (npartition_x == 1) ? _width : _min(blksz_x, _width);
        }
    }

    ImgIdx p = _parentAry[0];
    while (_node[p].parentIdx != ROOTIDX)
        p = _node[p].parentIdx;

    return p;
}

template <class Pixel>
ImgIdx AlphaTree<Pixel>::merge_subtrees1(uint8_t *dimg, int64_t blksz_x, int64_t blksz_y, int16_t npartition_x,
                                         int16_t npartition_y, ImgIdx *subtree_cur, int tse, ImgIdx *hypernode_level) {
    ImgIdx numblk;
    int64_t blksz_x0 = blksz_x;
    int64_t blksz_y0 = blksz_y;

    ImgIdx npartition_x0 = npartition_x;
    ImgIdx npartition_y0 = npartition_y;
    ImgIdx blkrow = _width * blksz_y0;
    while (npartition_x > 1 || npartition_y > 1) {
        if ((npartition_x == 1 || blksz_x >= blksz_y) && npartition_y > 1) {
            numblk = npartition_x * (npartition_y / 2);

#pragma omp parallel for schedule(dynamic, 1)
            for (int blk = 0; blk < numblk; blk++) {
                ImgIdx x, y, r, p, p0, pn, dimgidx, bx, by, bidx;
                y = (1 + 2 * (blk / (int)npartition_x)) * blksz_y;
                x = (blk % (int)npartition_x) * blksz_x;

                p0 = (y - 1) * _width + x;
                pn = (((blk % (int)npartition_x) == npartition_x - 1) ? y * _width : p0 + blksz_x);

                by = _min((p0 / blkrow), npartition_y0 - 1);
                for (p = p0; p < pn; p++) {
                    bx = _min(((p % _width) / blksz_x0), npartition_x0 - 1);
                    bidx = by * npartition_x0 + bx;

                    dimgidx = p << 1;
                    r = dimg[dimgidx];

                    if (tse)
                        hypernode_level[dimgidx] =
                            connect(_parentAry[p], _parentAry[p + _width], (ImgIdx)subtree_cur[bidx]++, (Pixel)r);
                    else
                        hypernode_level[dimgidx] = connect(p, p + _width, (ImgIdx)subtree_cur[bidx]++, (Pixel)r);
                }
            }
            npartition_y = (npartition_y + 1) / 2;
            blksz_y <<= 1;
            if (npartition_y == 1)
                blksz_y = _height;
            else
                blksz_y = _min(blksz_y, _height);
        }

        if ((npartition_y == 1 || blksz_x <= blksz_y) && npartition_x > 1) {
            numblk = npartition_y * (npartition_x / 2);

#pragma omp parallel for schedule(dynamic, 1)
            for (int blk = 0; blk < numblk; blk++) {
                ImgIdx x, y, r, p, p0, pn, dimgidx, bx, by, bidx;
                ;
                x = (1 + 2 * (blk / npartition_y)) * blksz_x;
                y = (blk % (int)npartition_y) * blksz_y;

                p0 = y * _width + x - 1;
                pn = ((blk % (int)npartition_y) == npartition_y - 1) ? _height * _width : p0 + _width * blksz_y;

                bx = _min(((p0 % _width) / blksz_x0), npartition_x0 - 1);
                for (p = p0; p < pn; p += _width) {
                    by = _min((p / blkrow), npartition_y0 - 1);
                    bidx = by * npartition_x0 + bx;

                    dimgidx = (p << 1) + 1;
                    r = dimg[dimgidx];

                    if (tse)
                        hypernode_level[dimgidx] =
                            connect(_parentAry[p], _parentAry[p + 1], (ImgIdx)subtree_cur[bidx]++, (Pixel)r);
                    else
                        hypernode_level[dimgidx] = connect(p, p + 1, (ImgIdx)subtree_cur[bidx]++, (Pixel)r);
                }
            }
            npartition_x = (npartition_x + 1) / 2;
            blksz_x <<= 1;
            if (npartition_x == 1)
                blksz_x = _width;
            else
                blksz_x = _min(blksz_x, _width);
        }
    }

    ImgIdx p;
    if (tse)
        p = _parentAry[0];
    else
        p = 0;
    while (_node[p].parentIdx != ROOTIDX)
        p = _node[p].parentIdx;

    return p;
}

template <class Pixel>
int AlphaTree<Pixel>::migrate_subtree(int blk, int numpartitions, ImgIdx &nidx, ImgIdx &nidx_lim, int &nidxblk,
                                      ImgIdx &blkts, char *blkflooddone, ImgIdx *subtree_cur, ImgIdx *subtree_start,
                                      ImgIdx *subtree_nborderedges, omp_lock_t *locks, int &numbusythr, int &numblkproc,
                                      int &outofmemory) {
    if (omp_get_num_threads() == 1)
        outofmemory = 1;

    subtree_cur[nidxblk] = nidx;
    blkts += nidx - subtree_start[nidxblk];
    blkflooddone[nidxblk] = 2;

    omp_set_lock(locks + numpartitions);
    numbusythr--;
    omp_unset_lock(locks + numpartitions);
    omp_unset_lock(locks + nidxblk);

    while (1) {
        if (outofmemory || (!numbusythr)) // && numblkproc >= omp_get_num_threads()))
        {
            outofmemory = 1;
            return 0;
        }

        for (int newblkidx = 0; newblkidx < numpartitions; newblkidx++) {
            if (blkflooddone[newblkidx] &&
                (subtree_cur[newblkidx] + subtree_nborderedges[newblkidx] < subtree_start[newblkidx + 1]) &&
                !omp_test_lock(locks + newblkidx)) {
                omp_set_lock(locks + numpartitions);

                if (blkflooddone[newblkidx]) {
                    blkflooddone[newblkidx] = 0;
                    numbusythr++;
                    omp_unset_lock(locks + numpartitions);
                } else {
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

template <class Pixel>
ImgIdx AlphaTree<Pixel>::parflood_node_alloc(ImgIdx *subtree_size, ImgIdx *subtree_start, ImgIdx *blockWidths,
                                             ImgIdx *blockHeights, int numpartitions, double sizemult) {
    subtree_start[0] = 0;
    for (int blk = 0; blk < numpartitions; blk++) {
        ImgIdx blkmaxsize = 1 + 3 * blockWidths[blk] * blockHeights[blk];
        subtree_start[blk + 1] = _min(blkmaxsize, (ImgIdx)((double)subtree_size[blk] * sizemult)) + subtree_start[blk];
    }
    _maxSize = subtree_start[numpartitions];
    if (_node)
        Free(_node);
    _node = (AlphaNode<Pixel> *)Calloc((size_t)(_maxSize) * sizeof(AlphaNode<Pixel>));

    return _maxSize;
}

template <class Pixel> void AlphaTree<Pixel>::printBinary(uint8_t num) {
    // Loop through each bit (starting from the most significant bit)
    for (int i = sizeof(num) * 8 - 1; i >= 0; i--) {
        // Use bitwise AND and right shift to extract the bit at position i
        unsigned int bit = (num >> i) & 1;
        printf("%u", bit);
    }
}

template <class Pixel> void AlphaTree<Pixel>::printIsAvailable(const uint8_t *isAvailable) const {
    printf("AlphaTree<Pixel>::printIsAvailable\n");
    for (int i = 0; i < _height; i++) {
        for (int j = 0; j < _width; j++) {
            printBinary(isAvailable[i * _width + j]);
            printf(" ");
        }
        printf("\n");
    }
    printf("\n");
}

template <class Pixel>
void AlphaTree<Pixel>::set_isAvailable_par(uint8_t *isAvailable, int16_t npartition_x, int16_t npartition_y) {
    int32_t i, j, k;
    const ImgIdx imgSize = _width * _height;
    const ImgIdx wstride = _width / npartition_x;
    const ImgIdx hstride = _height / npartition_y;

    set_isAvailable(isAvailable);

    if (_connectivity == 4) {
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
        const uint8_t maskBottom = 0b00001110;
        const uint8_t maskTop = 0b00000111;
        const uint8_t maskLeft = 0b00001011;
        const uint8_t maskRight = 0b00001101;

        // hor partitions
        j = (hstride - 1) * _width;
        for (i = 0; i < npartition_y - 1; i++) {
            k = j + _width;
            for (; j < k; j++) {
                isAvailable[j] &= maskBottom;
                isAvailable[j + _width] &= maskTop;
            }
            j += (hstride - 1) * _width;
        }

        // ver partitions
        for (i = 0; i < npartition_x - 1; i++) {
            j = (i + 1) * wstride - 1;
            for (; j < imgSize; j += _width) {
                isAvailable[j] &= maskRight;
                isAvailable[j + 1] &= maskLeft;
            }
        }
    } else {
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

        const uint8_t maskBottom = 0b01111100;
        const uint8_t maskTop = 0b11000111;
        const uint8_t maskLeft = 0b00011111;
        const uint8_t maskRight = 0b11110001;

        // hor partitions
        j = (hstride - 1) * _width;
        for (i = 0; i < npartition_y - 1; i++) {
            k = j + _width;
            for (; j < k; j++) {
                isAvailable[j] &= maskBottom;
                isAvailable[j + _width] &= maskTop;
            }
            j += (hstride - 1) * _width;
        }

        // ver partitions
        for (i = 0; i < npartition_x - 1; i++) {
            j = (i + 1) * wstride - 1;
            for (; j < imgSize; j += _width) {
                isAvailable[j] &= maskRight;
                isAvailable[j + 1] &= maskLeft;
            }
        }
    }
}

template <class Pixel> std::pair<ImgIdx, ImgIdx> AlphaTree<Pixel>::computePartitionSize(int numthreads) {
    const ImgIdx imgSize = _width * _height;
    ImgIdx optpart = 1;
    double optborderlength = (double)numthreads * (double)imgSize;
    for (int px = 2; px < numthreads; px++) {
        if (numthreads % px == 0) {
            int py = numthreads / px;

            if (((double)px * (double)_height + (double)py * (double)_width) < optborderlength) {
                optpart = px;
                optborderlength = ((double)px * (double)_height + (double)py * (double)_width);
            }
        }
    }
    ImgIdx npartition_x = (ImgIdx)optpart;
    ImgIdx npartition_y = (ImgIdx)numthreads / npartition_x;

    return {npartition_x, npartition_y};
}

template <class Pixel>
void AlphaTree<Pixel>::setBlockDimensions(ImgIdx npartition_x, ImgIdx npartition_y, ImgIdx blksz_x, ImgIdx blksz_y,
                                          ImgIdx *blocksize, ImgIdx *blockWidths, ImgIdx *blockHeights,
                                          ImgIdx *startpidx, ImgIdx *subtree_start) {
    const ImgIdx blksz_xn = blksz_x + (_width % npartition_x);
    const ImgIdx blksz_yn = blksz_y + (_height % npartition_y);
    const ImgIdx numpartitions = npartition_x * npartition_y;
    const ImgIdx imgSize = _width * _height;
    const ImgIdx numNeighbors = _connectivity / 2;
    ImgIdx p = 0;
    for (ImgIdx y = 0; y < npartition_y; y++) {
        const ImgIdx q = y * _width * (ImgIdx)blksz_y;
        const bool lastrow = (y == npartition_y - 1);
        const ImgIdx blockHeight = lastrow ? blksz_yn : blksz_y;
        for (ImgIdx x = 0; x < npartition_x; x++) {
            startpidx[p] = q + (ImgIdx)x * (ImgIdx)blksz_x;
            const bool lastcol = (x == npartition_x - 1);
            const ImgIdx blockWidth = lastcol ? blksz_xn : blksz_x;
            blocksize[p] =
                blockHeight * blockWidth * numNeighbors - (ImgIdx)lastrow * blockWidth - (ImgIdx)lastcol * blockHeight;
            blockWidths[p] = blockWidth;
            blockHeights[p] = blockHeight;
            p++;
        }
    }

    subtree_start[0] = startpidx[0] + imgSize;
    for (int blk = 1; blk < numpartitions; blk++) {
        subtree_start[blk] = subtree_start[blk - 1] + (blockWidths[blk - 1] * blockHeights[blk - 1] * numNeighbors);
    }
    subtree_start[numpartitions] = subtree_start[numpartitions - 1] +
                                   (blockWidths[numpartitions - 1] * blockHeights[numpartitions - 1] * numNeighbors);
}

// template <class Pixel>
// void AlphaTree<Pixel>::computePartitionDiffLogBucket(const Pixel *img, double *dimg, LogBucket &lb,
//                                                      ImgIdx startPixelIndex, ImgIdx blockWidth, ImgIdx blockHeight) {
//     const int shamt = _connectivity / 4;
//     // Pixel maxdiff = std::numeric_limits<Pixel>::max();
//     if (_connectivity == 4) {
//         for (ImgIdx i = 0; i < blockHeight; i++) {
//             const bool bottom = i < blockHeight - 1;
//             for (ImgIdx j = 0; j < blockWidth; j++) {
//                 const bool right = j < blockWidth - 1;
//                 const ImgIdx p = startPixelIndex + i * _width + j;
//                 const ImgIdx q = p << shamt;

//                 if (bottom) {
//                     const double diff = _pixelDissim.computeDissimilarity(img[p], img[p + _width]);
//                     dimg[q] = diff;
//                     lb.push(q, diff);
//                 }
//                 if (right) {
//                     const double diff = _pixelDissim.computeDissimilarity(img[p], img[p + 1]);
//                     dimg[q + 1] = diff;
//                     lb.push(q + 1, diff);
//                 }
//             }
//         }
//     } else if (_connectivity == 8) {
//         for (ImgIdx i = 0; i < blockHeight; i++) {
//             const bool top = i > 0;
//             const bool bottom = i < blockHeight - 1;
//             for (ImgIdx j = 0; j < blockWidth; j++) {
//                 const bool right = j < blockWidth - 1;
//                 const ImgIdx p = startPixelIndex + i * _width + j;
//                 const ImgIdx q = p << shamt;

//                 if (bottom) {
//                     const double diff = _pixelDissim.computeDissimilarity(img[p], img[p + _width]);
//                     dimg[q] = diff;
//                     lb.push(q, diff);
//                 }
//                 if (bottom && right) {
//                     const double diff = _pixelDissim.computeDissimilarity(img[p], img[p + _width + 1]);
//                     dimg[q + 1] = diff;
//                     lb.push(q + 1, diff);
//                 }
//                 if (right) {
//                     const double diff = _pixelDissim.computeDissimilarity(img[p], img[p + 1]);
//                     dimg[q + 2] = diff;
//                     lb.push(q + 2, diff);
//                 }
//                 if (top && right) {
//                     const double diff = _pixelDissim.computeDissimilarity(img[p], img[p - _width + 1]);
//                     dimg[q + 3] = diff;
//                     lb.push(q + 3, diff);
//                 }
//             }
//         }
//     } else {
//     }
// }

template <class Pixel>
Pixel AlphaTree<Pixel>::computePartitionDifferences(const Pixel *img, Pixel *dimg, ImgIdx startPixelIndex,
                                                    ImgIdx blockWidth, ImgIdx blockHeight, ImgIdx *blockDiffHist) {
    const int shamt = _connectivity / 4;
    Pixel maxdiff = std::numeric_limits<Pixel>::max();
    if (_connectivity == 4) {
        for (ImgIdx i = 0; i < blockHeight; i++) {
            const bool bottom = i < blockHeight - 1;
            for (ImgIdx j = 0; j < blockWidth; j++) {
                const bool right = j < blockWidth - 1;
                const ImgIdx p = startPixelIndex + i * _width + j;
                const ImgIdx q = p << shamt;
                _node[p] = AlphaNode<Pixel>(img[p], 0.0);

                if (bottom) {
                    const Pixel diff = abs_diff(img[p], img[p + _width]);
                    dimg[q] = diff;
                    blockDiffHist[diff]++;
                }
                if (right) {
                    const Pixel diff = abs_diff(img[p], img[p + 1]);
                    dimg[q + 1] = diff;
                    blockDiffHist[diff]++;
                }
            }
        }
    } else if (_connectivity == 8) {
        for (ImgIdx i = 0; i < blockHeight; i++) {
            const bool top = i > 0;
            const bool bottom = i < blockHeight - 1;
            for (ImgIdx j = 0; j < blockWidth; j++) {
                const bool right = j < blockWidth - 1;
                const ImgIdx p = startPixelIndex + i * _width + j;
                const ImgIdx q = p << shamt;
                _node[p] = AlphaNode<Pixel>(img[p], 0.0);

                if (bottom) {
                    const Pixel diff = abs_diff(img[p], img[p + _width]);
                    dimg[q] = diff;
                    blockDiffHist[diff]++;
                }
                if (bottom && right) {
                    const Pixel diff = abs_diff(img[p], img[p + _width + 1]);
                    dimg[q + 1] = diff;
                    blockDiffHist[diff]++;
                }
                if (right) {
                    const Pixel diff = abs_diff(img[p], img[p + 1]);
                    dimg[q + 2] = diff;
                    blockDiffHist[diff]++;
                }
                if (top && right) {
                    const Pixel diff = abs_diff(img[p], img[p - _width + 1]);
                    dimg[q + 3] = diff;
                    blockDiffHist[diff]++;
                }
            }
        }
    } else {
    }

    blockDiffHist[maxdiff]++; // Room for subtree dummy root
    return maxdiff;
}

template <class Pixel>
void AlphaTree<Pixel>::floodPartition(const Pixel *img, const Pixel *dimg, ImgIdx startPixelIndex, int blockIndex,
                                      ImgIdx blockArea, Pixel maxdiff, HierarQueue *queue, const ImgIdx *subtree_start,
                                      ImgIdx *subtree_cur, uint8_t *isVisited, const uint8_t *isAvailable) {
    const ImgIdx imgSize = _width * _height;
    const ImgIdx dimgSize = (_connectivity >> 1) * _width * _height;
    ImgIdx stackTop = imgSize + dimgSize + blockIndex;
    ImgIdx prevTop = stackTop;
    _node[stackTop] = AlphaNode<Pixel>(maxdiff);
    ImgIdx nodeIndex = subtree_start[blockIndex];
    Pixel currentLevel = maxdiff;
    queue->push(startPixelIndex, currentLevel);
    while (_node[stackTop].area < blockArea) // Main flood loop
    {
        while ((int64_t)queue->min_level <= (int64_t)currentLevel) { // Flood all levels below currentLevel
            const ImgIdx p = queue->pop();
            if (is_visited(isVisited, p)) {
                queue->find_minlev();
                continue;
            }

            // printf("Visiting %d at %d, size = %d\n", p, (int)currentLevel,
            //        (int)(nodeIndex - subtree_start[blockIndex]));
            // printGraph(isVisited, (uint8_t *)nullptr, img);
            // std::getchar();

            isVisited[p] = 1;
            uint8_t isAv = isAvailable[p];

            if (_connectivity == 4) {
                const ImgIdx q = p << 1;
                const ImgIdx width2 = _width << 1;
                // clang-format off
                if (is_available(isAv, 0) && !isVisited[p + _width])    queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])         queue->push(p + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])         queue->push(p - 1, dimg[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])    queue->push(p - _width, dimg[q - width2]);
                // clang-format on
            } else if (_connectivity == 8) // To do later
            {
                const ImgIdx width4 = _width << 2;
                const ImgIdx q = p << 2;
                // clang-format off
                if (is_available(isAv, 0) && !isVisited[p + _width])        queue->push(p + _width, dimg[q]);
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])    queue->push(p + _width + 1, dimg[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p + 1])             queue->push(p + 1, dimg[q + 2]);
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])    queue->push(p - _width + 1, dimg[q + 3]);
                if (is_available(isAv, 4) && !isVisited[p - _width])        queue->push(p - _width, dimg[q - width4]);
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])    queue->push(p - _width - 1, dimg[q - width4 - 3]);
                if (is_available(isAv, 6) && !isVisited[p - 1])             queue->push(p - 1, dimg[q - 2]);
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])    queue->push(p + _width - 1, dimg[q + width4 - 1]);
                // clang-format on
            }

            if ((int64_t)currentLevel > (int64_t)queue->min_level) { // Go to lower level
                currentLevel = queue->min_level;
                const Pixel pixVal = _node[p].minPix;
                const ImgIdx iNode = nodeIndex++;
                _node[iNode] = AlphaNode<Pixel>(pixVal, currentLevel, stackTop);
                _node[p].parentIdx = iNode;
                stackTop = iNode;
            } else {
                queue->find_minlev();
                _node[stackTop].add(_node[p]);
                _node[p].parentIdx = stackTop;
            }
        }

        // remove_redundant_node(_node, nodeIndex, prevTop, stackTop);
        if (_node[prevTop].parentIdx == stackTop && _node[prevTop].area == _node[stackTop].area) {
            _node[prevTop].parentIdx = _node[stackTop].parentIdx;
            stackTop = prevTop;
            nodeIndex--;
        }

        // Go to higher level
        if (_node[stackTop].area < blockArea) {
            ImgIdx newParentIdx = _node[stackTop].parentIdx;
            if (newParentIdx == ROOTIDX ||
                (int64_t)queue->min_level < (int64_t)_node[newParentIdx].alpha) { // New level from queue
                newParentIdx = nodeIndex++;
                _node[newParentIdx] = AlphaNode<Pixel>(_node[stackTop]);
                _node[newParentIdx].alpha = queue->min_level;
                _node[newParentIdx].parentIdx = _node[stackTop].parentIdx;
                _node[newParentIdx]._rootIdx = ROOTIDX;
                _node[stackTop].parentIdx = newParentIdx;
            } else
                _node[newParentIdx].add(_node + stackTop);

            prevTop = stackTop;
            stackTop = newParentIdx;
            currentLevel = _node[stackTop].alpha;
        }
    }
    ImgIdx subTreeRootIdx = (_node[prevTop].area == blockArea) ? prevTop : stackTop; // remove redundant root
    _node[subTreeRootIdx].parentIdx = ROOTIDX;

    subtree_cur[blockIndex] = nodeIndex;
}

template <class Pixel> void AlphaTree<Pixel>::FloodHierQueueParallel(const Pixel *img, int numthreads) {
    if (sizeof(Pixel) > 2 || _channel > 1) {
        printf("Error: Hierarchical queues do not work on >16 bits images or multispectral images\n");
        printf("Try Unionfind (algorithm code %d), flooding using Heapqueue (%d), trie queue (%d) or cached trie queue "
               "(%d) \n",
               UNIONFIND, FLOOD_HEAPQUEUE_CACHE, FLOOD_TRIE, FLOOD_TRIE_CACHE);
        return;
    }

    const ImgIdx imgSize = _width * _height;
    const ImgIdx dimgSize = (_connectivity >> 1) * _width * _height;
    const ImgIdx numlevels = (sizeof(Pixel) == 1) ? 256 : 65536;

    Pixel *dimg = (Pixel *)Calloc((size_t)dimgSize * sizeof(Pixel));

    const auto [npartition_x, npartition_y] = computePartitionSize(numthreads);

    uint8_t *isVisited = (uint8_t *)Calloc((size_t)(imgSize));
    uint8_t *isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable_par(isAvailable, npartition_x, npartition_y);
    // printIsAvailable(isAvailable);

    const ImgIdx blksz_x = _width / npartition_x;
    const ImgIdx blksz_y = _height / npartition_y;
    const ImgIdx numpartitions = npartition_x * npartition_y;

    ImgIdx *startpidx = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *blocksize = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *subtree_start = (ImgIdx *)Malloc((numpartitions + 1) * sizeof(ImgIdx));
    ImgIdx *subtree_cur = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *blockWidths = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *blockHeights = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *dhist = (ImgIdx *)Calloc((size_t)numlevels * (size_t)numpartitions * sizeof(ImgIdx));

    setBlockDimensions(npartition_x, npartition_y, blksz_x, blksz_y, blocksize, blockWidths, blockHeights, startpidx,
                       subtree_start);

    HierarQueue **queues;
    queues = (HierarQueue **)Calloc(numpartitions * sizeof(HierarQueue *));
    for (int blockIndex = 0; blockIndex < numpartitions; blockIndex++)
        queues[blockIndex] = new HierarQueue((uint64_t)blocksize[blockIndex] + 1);

    // singletons + inners + dummies
    _node = (AlphaNode<Pixel> *)Calloc((size_t)(imgSize + dimgSize + numpartitions) * sizeof(AlphaNode<Pixel>));
    _maxSize = imgSize + dimgSize + numpartitions;
    _curSize = imgSize + dimgSize + numpartitions;
    _nodeIn = _node + imgSize;

#pragma omp parallel for
    for (int blockIndex = 0; blockIndex < numpartitions; blockIndex++) {
        // printf("blockIndex = %d/%d\n", blockIndex, numpartitions);
        const ImgIdx blockWidth = blockWidths[blockIndex];
        const ImgIdx blockHeight = blockHeights[blockIndex];
        const ImgIdx blockArea = blockWidth * blockHeight;
        const ImgIdx startPixelIndex = startpidx[blockIndex];
        ImgIdx *blockDiffHist = dhist + numlevels * blockIndex;
        HierarQueue *queue = queues[blockIndex];
        const Pixel maxdiff =
            computePartitionDifferences(img, dimg, startPixelIndex, blockWidth, blockHeight, blockDiffHist);

        queue->set_queue(blockDiffHist, maxdiff);

        floodPartition(img, dimg, startPixelIndex, blockIndex, blockArea, maxdiff, queue, subtree_start, subtree_cur,
                       isVisited, isAvailable);
    }

    mergePartition(dimg, blksz_x, blksz_y, npartition_x, npartition_y, subtree_cur);

    // printTree();
    // printAll(isVisited, nullptr, img);

    Free(isVisited);
    Free(isAvailable);

    for (int blockIndex = 0; blockIndex < numpartitions; blockIndex++)
        delete queues[blockIndex];
    Free(queues);
    Free(dimg);
    Free(startpidx);
    Free(blocksize);
    Free(subtree_start);
    Free(subtree_cur);
    Free(blockWidths);
    Free(blockHeights);
    Free(dhist);
}

// Find subtree root and do path compression
template <class Pixel> ImgIdx AlphaTree<Pixel>::find_root(ImgIdx p) {
    if (p == ROOTIDX)
        return ROOTIDX;

    ImgIdx r, q;

    for (r = p; _node[r]._rootIdx != ROOTIDX; r = _node[r]._rootIdx)
        ;

    while (p != r) {
        q = _node[p]._rootIdx;
        _node[p]._rootIdx = r;
        p = q;
    }

    return r;
}

template <class Pixel> ImgIdx AlphaTree<Pixel>::find_root_in(ImgIdx p) {
    if (p == ROOTIDX)
        return ROOTIDX;

    ImgIdx r, q;

    for (r = p; _nodeIn[r]._rootIdx != ROOTIDX; r = _nodeIn[r]._rootIdx)
        ;

    while (p != r) {
        q = _nodeIn[p]._rootIdx;
        _nodeIn[p]._rootIdx = r;
        p = q;
        // int2++;
    }

    return r;
}

template <class Pixel> void AlphaTree<Pixel>::Unionfind(const Pixel *img) {
    ImgIdx imgSize, nredges;
    RankItem<double> *rankitem, *pRank;

    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    _curSize = _maxSize = imgSize + nredges; // to be compatible with flooding algorithms

    omp_set_num_threads(1);
    rankitem = (RankItem<double> *)Malloc(nredges * sizeof(RankItem<double>));
    compute_difference_and_sort(rankitem, img, nredges);

    // initialize_node(img, rankitem, maxpixval);
    _maxSize = imgSize + nredges;
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));
    for (ImgIdx p = 0; p < imgSize; p++) {
        _node[p].set(1, 0, (double)img[p], img[p], img[p]);
        _node[p]._rootIdx = _node[p].parentIdx = ROOTIDX;
    }

    bool unionbyrank = 0; //(sizeof(Pixel) <= 2);
    ImgIdx *treedepth = 0;

    if (unionbyrank)
        treedepth = (ImgIdx *)Calloc(_maxSize * sizeof(ImgIdx));

    ImgIdx _curSize = 0;
    for (ImgIdx r = 0; r < nredges; r++) {
        pRank = rankitem + r;

        ImgIdx x, x0;
        ImgIdx y, y0;
        ImgIdx z;

        ImgIdx nodeaddr = _curSize + imgSize;

        x0 = pRank->get_pidx0(_connectivity);
        y0 = pRank->get_pidx1(_width, _connectivity);
        x = find_root(x0);
        y = find_root(y0);

        if (x == y) // already connected, nothing to do
            continue;

        if (x < y) {
            z = x;
            x = y;
            y = z;
        }

        if (!unionbyrank || _node[x].alpha != pRank->alpha) {
            _curSize++;
            // _nodeIn[r].set(0, pRank->alpha, 0.0, maxpixval, 0);
            _node[nodeaddr].copy(_node + x);
            _node[nodeaddr].alpha = pRank->alpha;
            _node[nodeaddr].parentIdx = _node[nodeaddr]._rootIdx = ROOTIDX;
            //			_nodeIn[r].parentIdx = _node[x].parentIdx;
            _node[x].parentIdx = _node[x]._rootIdx = nodeaddr;
            _node[nodeaddr].add(_node + y);
            _node[y].parentIdx = _node[y]._rootIdx = nodeaddr;
            if (unionbyrank)
                treedepth[nodeaddr] = _max(treedepth[x], treedepth[y]) + 1;

            if (_node[nodeaddr].area == imgSize) {
                _rootIdx = nodeaddr;
                break;
            }
        } else {
            if (_node[x].alpha == _node[y].alpha && treedepth[x] < treedepth[y]) {
                z = x;
                x = y;
                y = z;
            }
            _node[x].add(_node + y);
            _node[y].parentIdx = _node[y]._rootIdx = x;
            treedepth[x] = _max(treedepth[x], treedepth[y] + 1);

            if (_node[x].area == imgSize) {
                _rootIdx = x;
                break;
            }
        }
    }
    if (unionbyrank)
        Free(treedepth);
    Free(rankitem);
    //		Free(isVisited);
    //		Free(isAvailable);
}

template <class Pixel>
void AlphaTree<Pixel>::blockwise_tse(ImgIdx *subtree_size, ImgIdx *subtree_nborderedges, double *nrmsds, ImgIdx *dhist,
                                     ImgIdx *subtree_max, ImgIdx *blockWidths, ImgIdx *blockHeights,
                                     int8_t npartition_x, int8_t npartition_y, ImgIdx numbins) {
    int64_t numpartitions = (int64_t)npartition_x * (int64_t)npartition_y;
    Pixel maxdiff;

    for (int blk = 0; blk < numpartitions; blk++) {
        maxdiff = subtree_max[blk];
        ImgIdx *blockDiffHist = dhist + numbins * blk;
        ImgIdx bhistsum = 0;
        ImgIdx blockWidth = blockWidths[blk];
        ImgIdx blockHeight = blockHeights[blk];
        // bool lastcol = (blk % npartition_x) == (npartition_x - 1);
        // bool lastrow = (blk / npartition_x) == (npartition_y - 1);

        for (int ii = 0; ii <= (int)maxdiff; ii++)
            bhistsum += blockDiffHist[ii];

        double nrmsd_blk = 0;
        for (int ii = 0; ii < (int)(maxdiff + 1); ii++)
            nrmsd_blk += ((double)dhist[ii]) * ((double)dhist[ii]);
        nrmsd_blk = sqrt((nrmsd_blk - (double)bhistsum) / ((double)bhistsum * ((double)bhistsum - 1.0)));
        nrmsd_blk = ((TSE_A * exp(TSE_SIGMA * nrmsd_blk) + TSE_B) + TSE_M);
        nrmsds[blk] = nrmsd_blk;

        // ImgIdx nrboderedges = (lastcol ? 0 : blockHeight) + (lastrow ? 0 : blockWidth);
        ImgIdx nrboderedges = 2 * (blockHeight + blockWidth);
        subtree_nborderedges[blk] = nrboderedges;
        subtree_size[blk] =
            TreeSizeEstimation(blockDiffHist, maxdiff + 1, blockWidth * blockHeight, bhistsum, 0.5) + nrboderedges;
        // subtree_size[blk] = TreeSizeEstimation(blockDiffHist, maxdiff + 1, blockWidth * blockHeight, bhistsum, 0.5);
        // //no need to add room for borders (already computed in quantization) printf("Subblock %d: blockDiffHist sum =
        // %d, size estimate: %d\n", (int)blk, (int)bhistsum, (int)subtree_size[blk]);
    }
}

// blockwise quantization and histogram computation (for pilot_rank)
template <class Pixel>
void AlphaTree<Pixel>::quantize_ranks_compute_histogram(uint8_t *qrank, ImgIdx *rank, const Pixel *img, ImgIdx *dhist,
                                                        ImgIdx *blockWidths, ImgIdx *blockHeights, ImgIdx *startpidx,
                                                        int64_t binsize, ImgIdx numbins, int8_t npartition_x,
                                                        int8_t npartition_y, ImgIdx *subtree_max) {
    int64_t numpartitions = (int64_t)npartition_x * (int64_t)npartition_y;

#pragma omp parallel for
    for (int blk = 0; blk < numpartitions; blk++) {
        ImgIdx p, q;
        ImgIdx blockWidth = blockWidths[blk];
        ImgIdx blockHeight = blockHeights[blk];
        ImgIdx *blockDiffHist = dhist + numbins * blk;
        bool lastcol = (blk % npartition_x) == (npartition_x - 1);
        bool lastrow = (blk / npartition_x) == (npartition_y - 1);
        ImgIdx startPixelIndex = startpidx[blk];
        Pixel maxdiff = 0;

        for (ImgIdx i = 0; i < blockHeight; i++) {
            ImgIdx r;
            uint8_t qr;
            p = startPixelIndex + i * _width;
            bool blklastrow = (i == blockHeight - 1);
            for (ImgIdx j = 0; j < blockWidth - 1; j++) {
                q = p << 1;
                if (i < blockHeight - 1 || !lastrow) {
                    r = rank[q];
                    qr = QUANTIZE_RANK(r, binsize);
                    qrank[q] = qr;
                    blockDiffHist[qr]++;
                    if (!blklastrow)
                        maxdiff = _max(maxdiff, qr);
                }
                r = rank[q + 1];
                qr = QUANTIZE_RANK(r, binsize);
                qrank[q + 1] = qr;
                blockDiffHist[qr]++;
                maxdiff = _max(maxdiff, qr);
                p++;
            }

            q = p << 1;
            if (i < blockHeight - 1 || !lastrow) {
                r = rank[q];
                qr = QUANTIZE_RANK(r, binsize);
                qrank[q] = qr;
                blockDiffHist[qr]++;
                if (!blklastrow)
                    maxdiff = _max(maxdiff, qr);
            }
            if (!lastcol) {
                r = rank[q + 1];
                qr = QUANTIZE_RANK(r, binsize);
                qrank[q + 1] = qr;
                blockDiffHist[qr]++;
            }
        }
        blockDiffHist[maxdiff]++; // dummy for subroot
        subtree_max[blk] = maxdiff;
    }
}

template <class Pixel> uint8_t AlphaTree<Pixel>::pow_quantization(ImgIdx rank, uint64_t qint) {
    return (uint8_t)(((double)rank * (double)rank) / (double)qint);
}

template <class Pixel>
void AlphaTree<Pixel>::pow_quantize_ranks(uint8_t *qrank, ImgIdx *rank, int64_t dimgSize, int64_t qint) {
#pragma omp parallel for
    for (ImgIdx i = 0; i < (ImgIdx)dimgSize; i++)
        qrank[i] = pow_quantization(rank[i], qint);
}

template <class Pixel>
ImgIdx AlphaTree<Pixel>::find_root(AlphaNode<Pixel> *pilottree, ImgIdx p, Pixel below_this_qlevel) {
    ImgIdx q = _parentAry[p], r;

    // int cnt = 0;

    while (pilottree[(r = pilottree[q].parentIdx)].alpha < below_this_qlevel) // fix qlevel (also sum)
    {
        q = r;
    }

    return q;
}

template <class Pixel> ImgIdx AlphaTree<Pixel>::descendroots(ImgIdx q, int64_t qlevel, AlphaNode<Pixel> *pilottree) {
    ImgIdx c = pilottree[q].parentIdx;
    while ((int64_t)pilottree[c].alpha < qlevel) {
        // int1++;
        q = c;
        c = pilottree[c].parentIdx;
    }
    return q;
}

// Hybrid_Pilot_Rank
template <class Pixel>
void AlphaTree<Pixel>::unionfind_refine_qlevel(int64_t qlevel, int64_t binsize, ImgIdx nredges,
                                               AlphaNode<Pixel> *pilottree, RankItem<double> *rankitem,
                                               int8_t *redundant_edge, int32_t *rankToIndex) {
    RankItem<double> *pRank;
    ImgIdx rank_start = (ImgIdx)qlevel * (ImgIdx)binsize;
    ImgIdx rank_end = std::min(nredges - 1, ((ImgIdx)qlevel + 1) * (ImgIdx)binsize - 1);
    ImgIdx imgSize = _width * _height;

    if (qlevel) {
        for (ImgIdx r = rank_start; r <= rank_end; r++) {
            ImgIdx ridx = r;
            if (rankToIndex)
                ridx = (ImgIdx)rankToIndex[r];

            pRank = rankitem + ridx;

            // look for ancestor at qlevel
            ImgIdx x, x0;
            ImgIdx y, y0;

            ImgIdx nodeaddr = r + imgSize;

            if (redundant_edge[rankitem[ridx].dimgidx]) {
                // printf("[skip]qlevel %d / rank %d(%d) - refining %d-%d\n", (int)qlevel, (int)r, (int)pRank->alpha,
                // (int)pRank->get_pidx0(_connectivity), (int)pRank->get_pidx1(_width,_connectivity));
                continue;
            }

            // printf("qlevel %d / rank %d(%d) - refining %d-%d\n", (int)qlevel, (int)r, (int)pRank->alpha,
            // (int)pRank->get_pidx0(_connectivity), (int)pRank->get_pidx1(_width,_connectivity));

            // non-zero level
            // if (qlevel)
            {
                x0 = find_root(pilottree, pRank->get_pidx0(_connectivity), qlevel);
                y0 = find_root(pilottree, pRank->get_pidx1(_width, _connectivity), qlevel);

                if (x0 == y0) // already connected, nothing to do
                {
                    continue;
                }

                x = find_root(pilottree[x0]._rootIdx);
                y = find_root(pilottree[y0]._rootIdx);
            }

            if ((x != ROOTIDX) && (x == y)) // already connected, nothing to do
            {
                continue;
            }

            // add the new _node to the refined tree
            {
                if (x == ROOTIDX) {
                    {
                        _nodeIn[r].copy(pilottree + x0);
                        _nodeIn[r].parentIdx = pilottree[x0].parentIdx;
                        pilottree[x0]._rootIdx = nodeaddr;
                    }
                } else {
                    _nodeIn[r].copy(_node + x);
                    _nodeIn[r].parentIdx = _node[x].parentIdx;
                    _node[x].parentIdx = _node[x]._rootIdx = nodeaddr;
                }

                // attach to y
                if (y == ROOTIDX) {
                    // if (qlevel)
                    {
                        _nodeIn[r].add(pilottree + y0);
                        pilottree[y0]._rootIdx = nodeaddr;
                    }
                } else {
                    _nodeIn[r].add(_node + y);
                    _node[y].parentIdx = _node[y]._rootIdx = nodeaddr;
                }
            }
        }
    } else {
        for (ImgIdx r = rank_start; r <= rank_end; r++) {
            if (rankToIndex)
                pRank = rankitem + rankToIndex[r];
            else
                pRank = rankitem + r;

            // look for ancestor at qlevel
            ImgIdx x, x0;
            ImgIdx y, y0;

            ImgIdx nodeaddr = r + imgSize;

            // else
            {
                // find subtree roots from two edge incidents
                x0 = pRank->get_pidx0(_connectivity);
                y0 = pRank->get_pidx1(_width, _connectivity);

                x = find_root(_node[x0]._rootIdx);
                y = find_root(_node[y0]._rootIdx);
            }

            if ((x != ROOTIDX) && (x == y)) // already connected, nothing to do
                continue;

            {
                // attach to x
                if (x == ROOTIDX) {
                    // else
                    {
                        _nodeIn[r].copy(_node + x0);
                        _nodeIn[r].parentIdx = _parentAry[x0];
                        _node[x0].parentIdx = _node[x0]._rootIdx = nodeaddr;
                    }
                } else {
                    _nodeIn[r].copy(_node + x);
                    _nodeIn[r].parentIdx = _node[x].parentIdx;
                    _node[x].parentIdx = _node[x]._rootIdx = nodeaddr;
                }

                // attach to y
                if (y == ROOTIDX) {
                    // else
                    {
                        _nodeIn[r].add(_node + y0);
                        _node[y0].parentIdx = _node[y0]._rootIdx = nodeaddr;
                    }
                } else {
                    _nodeIn[r].add(_node + y);
                    _node[y].parentIdx = _node[y]._rootIdx = nodeaddr;
                }
            }
        }
    }
}

// compute edge histogram of the quantized rank image.
// In the hypergraph implementation, edges on the subblock borders are also counted.
template <class Pixel>
void AlphaTree<Pixel>::compute_dhist_par(uint8_t *qrank, ImgIdx *dhist, ImgIdx *startpidx, int32_t numbins,
                                         int8_t npartition_x, int8_t npartition_y, int64_t blksz_x, int64_t blksz_y,
                                         int64_t blksz_xn, int64_t blksz_yn) {
    if (_connectivity == 4) {
        // gdhist = edge histogram for edges on the subblock borders
        ImgIdx *gdhist = dhist + (int)npartition_x * (int)npartition_y * numbins;
// ImgIdx blk = 0;
// for (ImgIdx y = 0;y < (ImgIdx)npartition_y;y++)
// for (ImgIdx x = 0;x < (ImgIdx)npartition_x;x++)
#pragma omp parallel for
        for (ImgIdx blk = 0; blk < (int)npartition_x * (int)npartition_y; blk++) {
            ImgIdx x = blk % npartition_x;
            ImgIdx y = blk / npartition_y;
            ImgIdx lastcol = (x == npartition_x - 1);
            ImgIdx lastrow = (y == npartition_y - 1);
            ImgIdx xn = lastcol ? blksz_xn : blksz_x;
            ImgIdx yn = lastrow ? blksz_yn : blksz_y;
            ImgIdx p0 = startpidx[blk] << 1, p;
            ImgIdx *pdhist = dhist + blk * numbins;

            // ImgIdx maxval = qrank[p0];
            // ImgIdx maxpidx;

            for (ImgIdx i = 0; i < yn - 1; i++) {
                p = p0 + _width * (i << 1);
                for (ImgIdx j = 0; j < xn - 1; j++) {
                    pdhist[qrank[p++]]++;
                    pdhist[qrank[p++]]++;
                }
                pdhist[qrank[p++]]++;
                if (!lastcol)
                    gdhist[qrank[p++]]++;
            }

            // the last row of the subblock is...
            if (lastrow) { // the last row of the image
                p = p0 + _width * ((yn - 1) << 1) + 1;
                for (ImgIdx j = 0; j < xn - 1; j++) {
                    pdhist[qrank[p]]++;
                    p += 2;
                }
                if (!lastcol)
                    gdhist[qrank[p++]]++;
            } else {
                p = p0 + _width * ((yn - 1) << 1);
                for (ImgIdx j = 0; j < xn - 1; j++) {
                    gdhist[qrank[p++]]++;
                    pdhist[qrank[p++]]++;
                }
                gdhist[qrank[p++]]++;
                if (!lastcol)
                    gdhist[qrank[p]]++;
            }
            // blk++;
        }
    } else {
        // 8-N later
    }
}

template <class Pixel>
void AlphaTree<Pixel>::compute_dhist_par_hypergraph(uint8_t *qrank, ImgIdx *dhist, ImgIdx *startpidx, int32_t numbins,
                                                    int8_t npartition_x, int8_t npartition_y, int64_t blksz_x,
                                                    int64_t blksz_y, int64_t blksz_xn, int64_t blksz_yn,
                                                    ImgIdx *blkmaxpidx) {
    if (_connectivity == 4) {
#pragma omp parallel for
        for (ImgIdx blk = 0; blk < (int)npartition_x * (int)npartition_y; blk++) {
            ImgIdx x = blk % npartition_x;
            ImgIdx y = blk / npartition_y;
            ImgIdx lastcol = (x == npartition_x - 1);
            ImgIdx lastrow = (y == npartition_y - 1);
            ImgIdx xn = lastcol ? blksz_xn : blksz_x;
            ImgIdx yn = lastrow ? blksz_yn : blksz_y;
            ImgIdx p0 = startpidx[blk] << 1, p;
            ImgIdx *pdhist = dhist + blk * numbins;

            ImgIdx maxval = qrank[p0];
            ImgIdx maxpidx = p0;

            for (ImgIdx i = 0; i < yn - 1; i++) // for subimage rows (except for the last)
            {
                p = p0 + _width * (i << 1);
                for (ImgIdx j = 0; j < xn - 1; j++) // for subimage cols (except for the last)
                {
                    if (maxval < qrank[p]) {
                        maxval = qrank[p];
                        maxpidx = p;
                    }
                    pdhist[qrank[p++]]++;
                    if (maxval < qrank[p]) {
                        maxval = qrank[p];
                        maxpidx = p;
                    }
                    pdhist[qrank[p++]]++;
                }
                if (maxval < qrank[p]) {
                    maxval = qrank[p];
                    maxpidx = p;
                }
                pdhist[qrank[p++]]++;
                if (!lastcol) {
                    if (maxval < qrank[p]) {
                        maxval = qrank[p];
                        maxpidx = p;
                    }
                    pdhist[qrank[p++]]++;
                }
            }

            // the last row of the subblock is...
            if (lastrow) { // the last row of the image
                p = p0 + _width * ((yn - 1) << 1) + 1;
                for (ImgIdx j = 0; j < xn - 1; j++) {
                    if (maxval < qrank[p]) {
                        maxval = qrank[p];
                        maxpidx = p;
                    }
                    pdhist[qrank[p]]++;
                    p += 2;
                }
                if (!lastcol) {
                    if (maxval < qrank[p]) {
                        maxval = qrank[p];
                        maxpidx = p;
                    }
                    pdhist[qrank[p++]]++;
                }
            } else {
                p = p0 + _width * ((yn - 1) << 1);
                for (ImgIdx j = 0; j < xn - 1; j++) {
                    if (maxval < qrank[p]) {
                        maxval = qrank[p];
                        maxpidx = p;
                    }
                    pdhist[qrank[p++]]++;
                    if (maxval < qrank[p]) {
                        maxval = qrank[p];
                        maxpidx = p;
                    }
                    pdhist[qrank[p++]]++;
                }

                if (maxval < qrank[p]) {
                    maxval = qrank[p];
                    maxpidx = p;
                }
                pdhist[qrank[p++]]++;
                if (!lastcol) {
                    if (maxval < qrank[p]) {
                        maxval = qrank[p];
                        maxpidx = p;
                    }
                    pdhist[qrank[p]]++;
                }
            }

            blkmaxpidx[blk] = maxpidx;
            // blk++;
        }
    } else {
        // L A T e r
    }
}

// obsolete code (subtree nodes are now indexed based on their level)
template <class Pixel>
void AlphaTree<Pixel>::fix_subtreeidx(ImgIdx *subtreestart, ImgIdx *startpidx, ImgIdx *cursizes, int8_t npartition_x,
                                      int8_t npartition_y, int numpartitions, int64_t blksz_x, int64_t blksz_y,
                                      int64_t blksz_xn, int64_t blksz_yn) {
    for (int b = 0; b < numpartitions; b++) {
        ImgIdx poffset = subtreestart[b];

        // fix _parentAry
        ImgIdx bx = ((b % npartition_x) == npartition_x - 1) ? blksz_xn : blksz_x;
        ImgIdx by = ((b / npartition_x) == npartition_y - 1) ? blksz_yn : blksz_y;
        ImgIdx pidx = startpidx[b];
        for (ImgIdx p = 0; p < by; p++) {
            for (ImgIdx q = 0; q < bx; q++) {
                _parentAry[pidx++] += poffset;
            }
            pidx += _width - bx;
        }

        // fix _node parentidxs
        AlphaNode<Pixel> *ptree = _node + poffset;
        for (ImgIdx p = 0; p < cursizes[b]; p++)
            if (ptree[p].parentIdx != ROOTIDX)
                ptree[p].parentIdx += poffset;
    }
}

// The one with the pilottree indicing (slow?)
template <class Pixel>
void AlphaTree<Pixel>::merge_subtrees(uint8_t *qrank, ImgIdx *qindex, int64_t blksz_x, int64_t blksz_y,
                                      ImgIdx neighbor_offset, ImgIdx shamt, ImgIdx npartition_x, ImgIdx npartition_y,
                                      int32_t numbins) {
    ImgIdx x, y, r, p, q, dimgidx;
    ImgIdx imgSize = _height * _width;

    // merging border(hor)
    for (y = blksz_y; y <= blksz_y * (npartition_y - 1); y += blksz_y) {
        q = y * _width;
        for (p = q - _width; p < q; p++) {
            dimgidx = (p << shamt) + neighbor_offset;
            r = qrank[dimgidx];

            connect(_parentAry[p], _parentAry[p + _width], (Pixel)r, (ImgIdx)qindex[r]++);
        }
    }
    neighbor_offset = (_connectivity == 4) ? 1 : 3;

    // merging border(ver)
    for (x = blksz_x - 1; x <= blksz_x * (npartition_x - 1); x += blksz_x) {
        q = x + imgSize;
        for (p = x; p < q; p += _width) {
            dimgidx = (p << shamt) + neighbor_offset;
            r = qrank[dimgidx];

            connect(_parentAry[p], _parentAry[p + 1], (Pixel)r, (ImgIdx)qindex[r]++);
        }
    }

    // reset rootidxs
    for (p = 0; p < _maxSize; p++)
        if (_node[p].area)
            _node[p]._rootIdx = ROOTIDX;

    // Look for the root
    for (p = _parentAry[0]; _node[p].area != imgSize; p = _node[p].parentIdx)
        ;
    _node[p]._rootIdx = _node[p].parentIdx = ROOTIDX;
    _rootIdx = p;
}

template <class Pixel>
void AlphaTree<Pixel>::merge_subtrees(uint8_t *qrank, int64_t blksz_x, int64_t blksz_y, ImgIdx neighbor_offset,
                                      ImgIdx shamt, ImgIdx npartition_x, ImgIdx npartition_y, int32_t numbins) {
    ImgIdx x, y, r, p, q, dimgidx;
    ImgIdx imgSize = _height * _width;

    // merging border(hor)
    for (y = blksz_y; y <= blksz_y * (npartition_y - 1); y += blksz_y) {
        q = y * _width;
        for (p = q - _width; p < q; p++) {
            dimgidx = (p << shamt) + neighbor_offset;
            r = qrank[dimgidx];

            connect(_parentAry[p], _parentAry[p + _width], (ImgIdx)_curSize++, (Pixel)r);
        }
    }
    neighbor_offset = (_connectivity == 4) ? 1 : 3;

    // merging border(ver)
    for (x = blksz_x - 1; x <= blksz_x * (npartition_x - 1); x += blksz_x) {
        q = x + imgSize;
        for (p = x; p < q; p += _width) {
            dimgidx = (p << shamt) + neighbor_offset;
            r = qrank[dimgidx];

            connect(_parentAry[p], _parentAry[p + 1], (ImgIdx)_curSize++, (Pixel)r);
        }
    }

    // reset rootidxs
    for (p = 0; p < _curSize; p++)
        _node[p]._rootIdx = ROOTIDX;

    // Make sure that the root has the highest
    for (p = _parentAry[0]; _node[p].parentIdx != ROOTIDX; p = _node[p].parentIdx)
        ;
    if (_node[p].alpha != (Pixel)(numbins - 1)) {
        q = _curSize++;
        _node[q].copy(_node + p);
        _node[q].alpha = numbins - 1;
        _node[p].parentIdx = q;
        _node[q].parentIdx = _node[q]._rootIdx = ROOTIDX;
    }
}

template <class Pixel>
void AlphaTree<Pixel>::connect_pilotnode(AlphaNode<Pixel> *pilottree, ImgIdx nredges, ImgIdx imgSize) {
    // _curSize = _maxSize;
    ImgIdx *rootindexcand = (ImgIdx *)Calloc(omp_get_max_threads() * sizeof(ImgIdx));
    for (int i = 0; i < omp_get_max_threads(); i++)
        rootindexcand[i] = ROOTIDX;

#pragma omp parallel for schedule(guided, 1)
    for (ImgIdx p = 0; p < _maxSize; p++) {
        ImgIdx q, r, s;
        // printf("p: %d\n",(int)p);
        if (p < imgSize) {
            q = _parentAry[p];
            if (_node[p].parentIdx == ROOTIDX && pilottree[q].area == 1)
                _node[p].parentIdx = pilottree[q]._rootIdx;
        } else {
            if (_node[p]._rootIdx == ROOTIDX && _node[p].area) {
                if (_node[p].area == imgSize) {
                    if (rootindexcand[omp_get_thread_num()] == ROOTIDX) {
                        rootindexcand[omp_get_thread_num()] = p;
                    }
                    continue;
                }
                q = _node[p].parentIdx;
                for (r = q; pilottree[r]._rootIdx == ROOTIDX; r = pilottree[r].parentIdx)
                    ;

                s = pilottree[r]._rootIdx;
                while (q != r) {
                    pilottree[q]._rootIdx = s;
                    q = pilottree[q].parentIdx;
                }
                _node[p].parentIdx = s;
            }
        }
    }

    _rootIdx = ROOTIDX;
    for (int i = 0; i < omp_get_max_threads(); i++) {
        if (_rootIdx == ROOTIDX)
            _rootIdx = rootindexcand[i];
        else if (rootindexcand[i] != ROOTIDX)
            _rootIdx = _min(_rootIdx, rootindexcand[i]);
    }

    Free(rootindexcand);
}

template <class Pixel>
void AlphaTree<Pixel>::set_qindex(ImgIdx *qindex, ImgIdx *dhist, int64_t numpartitions, int32_t numbins,
                                  ImgIdx npartition_x, ImgIdx npartition_y, int64_t blksz_x, int64_t blksz_y,
                                  int64_t blksz_xn, int64_t blksz_yn) {

    // add rooms for singleton nodes
    {
        int p = 0;
        for (int x = 0; x < npartition_x; x++) {
            ImgIdx sx = (x == npartition_x - 1) ? blksz_xn : blksz_x;
            for (int y = 0; y < npartition_y; y++) {
                ImgIdx sy = (y == npartition_y - 1) ? blksz_yn : blksz_y;
                ImgIdx blksize = sx * sy;
                dhist[p] += blksize;
                p += numbins;
            }
        }
    }

    // compute cumulative distribution
    for (int p = 0; p < numpartitions; p++) {
        ImgIdx *pdhist = dhist + p * numbins;
        for (int q = 0; q < numbins; q++) {
            pdhist[q + numbins] += pdhist[q];
        }
    }

    ImgIdx *cdhist = dhist + numpartitions * numbins; // histogram of the edges on the subblock borders
    for (int p = 0; p < numbins - 1; p++)
        cdhist[p + 1] += cdhist[p];
    for (int p = 0; p < numbins; p++) {
        ImgIdx hp = ((p > 0) ? cdhist[p - 1] : 0);
        for (int q = 0; q < numpartitions + 1; q++)
            qindex[p + numbins * q] = hp + ((q > 0) ? dhist[(q - 1) * numbins + p] : 0);
    }
}

template <class Pixel>
void AlphaTree<Pixel>::set_qindex(ImgIdx *qindex, ImgIdx *dhist, int64_t numpartitions, int32_t numbins) {
    ImgIdx imgSize = _width * _height;
    for (int p = 0; p < numpartitions; p++) {
        ImgIdx *pdhist = dhist + p * numbins;
        for (int q = 0; q < numbins; q++) {
            pdhist[q + numbins] += pdhist[q];
        }
    }
    ImgIdx *cdhist = dhist + numpartitions * numbins;
    for (int p = 0; p < numbins - 1; p++)
        cdhist[p + 1] += cdhist[p];
    for (int p = 0; p < numbins; p++) {
        ImgIdx hp = imgSize + ((p > 0) ? cdhist[p - 1] : 0);
        for (int q = 0; q < numpartitions + 1; q++)
            qindex[p + numbins * q] = hp + ((q > 0) ? dhist[(q - 1) * numbins + p] : 0);
    }
}

template <class Pixel>
void AlphaTree<Pixel>::set_subtree_root(ImgIdx **subtreerootary, ImgIdx *strary, ImgIdx nonzero_nodeidx_start,
                                        ImgIdx rootlevel_nodeidx_start) {
    ImgIdx imgSize = _width * _height;
    if (_node[_rootIdx].alpha < 2) {
        std::cout << "Pilottree root _node level lower than 2" << std::endl;
        return;
    }

    {

        for (ImgIdx p = 1; p < (int)_node[_rootIdx].alpha; p++)
            subtreerootary[p] = strary + (p - 1) * imgSize;

        for (ImgIdx p = 0; p <= _rootIdx; p++)
            _node[p]._rootIdx = ROOTIDX;

        for (ImgIdx p = 0; p < imgSize; p++) {
            ImgIdx q = _parentAry[p];
            _node[_node[q].parentIdx]._rootIdx = 1;
        }

        int64_t areasum = 0;
        for (ImgIdx p = nonzero_nodeidx_start; p < rootlevel_nodeidx_start; p++) {
            if (_node[p]._rootIdx != ROOTIDX) {
                _node[_node[p].parentIdx]._rootIdx = 0;
                areasum += _node[p].area;
                _node[p]._rootIdx = areasum;
            }
        }

        ImgIdx *pixelindex = (ImgIdx *)Malloc(areasum * sizeof(ImgIdx));
        for (ImgIdx p = 0; p < imgSize; p++) // 1-CCs
        {
            ImgIdx q = _node[_parentAry[p]].parentIdx;
            if (q < rootlevel_nodeidx_start)
                pixelindex[--_node[q]._rootIdx] = p;
        }

        for (ImgIdx p = nonzero_nodeidx_start; p < rootlevel_nodeidx_start; p++) {
            if (_node[p]._rootIdx == ROOTIDX)
                continue;
            {

                ImgIdx q = _node[p].parentIdx;
                if (q >= rootlevel_nodeidx_start)
                    continue;
                ImgIdx r = _node[p]._rootIdx;
                ImgIdx s = r + _node[p].area;
                ImgIdx t = _node[q]._rootIdx;

                // _node[p]._rootIdx = ROOTIDX;
                while (r < s) {
                    pixelindex[--t] = pixelindex[r++];
                }

                _node[q]._rootIdx = t;
            }
        }

        for (ImgIdx p = 0; p < (int)((_node[_rootIdx].alpha - 1) * imgSize); p++)
            strary[p] = -1;

        for (ImgIdx p = 0; p < imgSize; p++)
            strary[p] = _parentAry[p];
        for (ImgIdx p = nonzero_nodeidx_start; p < rootlevel_nodeidx_start; p++) {
            if (!_node[p].area)
                continue;

            ImgIdx r = _node[p]._rootIdx;
            ImgIdx s = _node[p]._rootIdx + _node[p].area;
            _node[p]._rootIdx = ROOTIDX;
            ImgIdx level = _node[p].alpha;
            while (r < s) {
                subtreerootary[level][pixelindex[r++]] = p;
            }
        }

        for (ImgIdx level = 2; level < (int)_node[_rootIdx].alpha; level++) {
            ImgIdx *pstrary_prev = &subtreerootary[level - 1][0];
            ImgIdx *pstrary = &subtreerootary[level][0];
            for (ImgIdx p = 0; p < imgSize; p++) {
                if (pstrary[p] == -1)
                    pstrary[p] = pstrary_prev[p];
            }
        }

        Free(pixelindex);
    }
}

// void init_hypergraph_nodes(uint8_t* is_redundant, ImgIdx *rank, RankItem<Pixel>* rankitem)
template <class Pixel> void AlphaTree<Pixel>::find_redundant_nodes(uint8_t *is_redundant, ImgIdx *rank) {
    if (_connectivity == 4) {
        ImgIdx imgSize = _height * _width;
        ImgIdx width2 = _width << 1;
        for (ImgIdx p = 0; p < imgSize; p++) {
            ImgIdx q = p << 1;
            ImgIdx y = p / _width;
            ImgIdx x = p % _width;
            // int8_t isAv = isAvailable[p];
            ImgIdx maxRank = -1;

            ((y < _height - 1) && (rank[q] > maxRank)) ? (maxRank = rank[q]) : (ImgIdx)0;
            ((x < _width - 1) && (rank[q + 1] > maxRank)) ? (maxRank = rank[q + 1]) : (ImgIdx)0;
            ((x > 0) && (rank[q - 1] > maxRank)) ? (maxRank = rank[q - 1]) : (ImgIdx)0;
            ((y > 0) && (rank[q - width2] > maxRank)) ? (maxRank = rank[q - width2]) : (ImgIdx)0;

            // is_redundant[rankitem[minRank].dimgidx] = 1;
            is_redundant[maxRank] = 1;
            // _node[p].connect_to_parent(&_nodeIn[minRank], minRank + imgSize);
        }
    } else // later!
    {
    }
}

template <class Pixel>
void AlphaTree<Pixel>::set_subblock_properties(ImgIdx *startpidx, ImgIdx *blockWidths, ImgIdx *blockHeights,
                                               ImgIdx *blocksize, int8_t npartition_x, int8_t npartition_y,
                                               int64_t blksz_x, int64_t blksz_y, int64_t blksz_xn, int64_t blksz_yn) {
    ImgIdx p = 0, q;
    for (int8_t y = 0; y < npartition_y; y++) {
        q = y * _width * (ImgIdx)blksz_y;
        bool lastrow = (y == npartition_y - 1);
        ImgIdx blockHeight = lastrow ? blksz_yn : blksz_y;
        for (int8_t x = 0; x < npartition_x; x++) {
            startpidx[p] = q + (ImgIdx)x * (ImgIdx)blksz_x;
            bool lastcol = (x == npartition_x - 1);
            ImgIdx blockWidth = lastcol ? blksz_xn : blksz_x;
            blocksize[p] = blockHeight * blockWidth * 2 - (ImgIdx)lastrow * blockWidth - (ImgIdx)lastcol * blockHeight;
            blockWidths[p] = blockWidth;
            blockHeights[p] = blockHeight;
            p++;
        }
    }
}

template <class Pixel>
void AlphaTree<Pixel>::memalloc_queues(HierarQueue ***queues, int64_t numpartitions, ImgIdx *blocksize,
                                       ImgIdx *subtree_max) {
    // preallocate hqueues
    *queues = (HierarQueue **)Calloc(numpartitions * sizeof(HierarQueue *));
    for (int blk = 0; blk < numpartitions; blk++) {
        // printf("blk%d max = %d\n", blk, (int)subtree_max[blk]);
        (*queues)[blk] = new HierarQueue((uint64_t)blocksize[blk] + 1, (int32_t)(subtree_max[blk] + 1));
        // queues[blk] = new HierarQueue((uint64_t)blocksize[blk] + 1, (int32_t)(1<<20));
    }
}

template <class Pixel>
void AlphaTree<Pixel>::computeDimgAndRankToIndex(RankItem<double> *&rankitem, const Pixel *img, ImgIdx nredges,
                                                 int32_t *rankToIndex) {
    SortValue<float> *vals; // = new pmt::SortValue<Value>[N];
    vals = (SortValue<float> *)Malloc(nredges * sizeof(SortValue<float>));
    computeDimgRadixf(rankitem, img, vals);

    SortPair<uint32_t, ImgIdx> *sort_space =
        (SortPair<uint32_t, ImgIdx> *)Calloc(2 * nredges * sizeof(SortPair<uint32_t, ImgIdx>)); // new SortPair[2 * N];

    rank_to_index((SortValue<uint32_t> *)vals, (ImgIdx)nredges, rankToIndex, 0U, (uint_fast8_t)(sizeof(float) << 3),
                  sort_space, omp_get_max_threads());

    Free(vals);
    Free(sort_space);
}

template <class Pixel>
void AlphaTree<Pixel>::compute_difference_and_sort(RankItem<double> *&rankitem, const Pixel *img, ImgIdx nredges) {
    int32_t *rankToIndex = (int32_t *)Calloc(nredges * sizeof(int32_t));

    computeDimgAndRankToIndex(rankitem, img, nredges, rankToIndex);

    RankItem<double> *sorteditem = (RankItem<double> *)Malloc(nredges * sizeof(RankItem<double>)), *tmp;
#pragma omp parallel for schedule(guided, 1)
    for (ImgIdx i = 0; i < nredges; i++) {
        sorteditem[i] = rankitem[rankToIndex[i]];
    }

    tmp = rankitem;
    rankitem = sorteditem;
    Free(tmp);
    Free(rankToIndex);
}
template <class Pixel>
void AlphaTree<Pixel>::computeDifferenceAndSort(ImgIdx *rank, RankItem<double> *&rankitem, const Pixel *img,
                                                ImgIdx nredges, int32_t *&rankToIndex) {
    computeDimgAndRankToIndex(rankitem, img, nredges, rankToIndex);

#pragma omp parallel for schedule(guided, 1)
    for (ImgIdx i = 0; i < nredges; i++) {
        rank[rankitem[rankToIndex[i]].dimgidx] = i;
    }
}

template <class Pixel>
std::pair<ImgIdx, uint8_t> AlphaTree<Pixel>::computeFloodingContextHybridPartition(bool &isFirstPixel,
                                                                                   const ImgIdx &qitem,
                                                                                   const uint8_t *isAvailable) {
    const int shamt = (_connectivity >> 2) + 1;
    ImgIdx p = qitem >> shamt;
    uint8_t isAv = isAvailable[p];

    if (isFirstPixel) {
        isFirstPixel = 0;
        return {p, isAv};
    }

    const bool isRedundant = qitem & 1;
    const bool edgeDirection = (qitem >> 1) & 1;

    if (isRedundant) {
        if (_connectivity == 4) {
            if (edgeDirection) // horizontal edge
            {
                p++;
                isAv = isAvailable[p];
                isAv &= ~0x4; // do not check the neighbor which pushed this pixel (p) into the
                              // queue
            } else            // vertical edge
            {
                p += _width;
                isAv = isAvailable[p];
                isAv &= ~0x8;
            }
        } else {
            isAv = 0;
        }
    } else {
        if (edgeDirection) // horizontal edge
        {
            isAv = isAvailable[p];
            isAv &= ~0x2;
        } else // vertical edge
        {
            isAv = isAvailable[p];
            isAv &= ~0x1;
        }
    }

    return {p, isAv};
}

template <class Pixel>
void AlphaTree<Pixel>::markRedundantEdges(int8_t *redundant_edge, const ImgIdx &p, const uint8_t &connected_neighbor,
                                          const uint8_t *qrank) {
    // mark from the leaf to the stackTop _node to help finding hypernode levelroots
    const ImgIdx leaf1 = _parentAry[p];
    ImgIdx squirrel;
    for (squirrel = leaf1; _node[squirrel].parentIdx != squirrel; squirrel = _node[squirrel].parentIdx)
        _node[squirrel]._rootIdx = p;
    _node[squirrel]._rootIdx = p;

    const ImgIdx shamt = _connectivity >> 2;
    const ImgIdx wstride_d = _width << shamt;
    const ImgIdx q = p << shamt;
    if (connected_neighbor & 0x1) {
        const ImgIdx leaf2 = _parentAry[p + _width];
        if (qrank[q] > get_nearest_common_ancestor_level(p, leaf2))
            redundant_edge[q] = 1;
    }
    if (connected_neighbor & 0x2) {
        const ImgIdx leaf2 = _parentAry[p + 1];
        if (qrank[q + 1] > get_nearest_common_ancestor_level(p, leaf2))
            redundant_edge[q + 1] = 1;
    }
    if (connected_neighbor & 0x4) {
        const ImgIdx leaf2 = _parentAry[p - 1];
        if (qrank[q - 1] > get_nearest_common_ancestor_level(p, leaf2))
            redundant_edge[q - 1] = 1;
    }
    if (connected_neighbor & 0x8) {
        const ImgIdx leaf2 = _parentAry[p - _width];
        if (qrank[q - wstride_d] > get_nearest_common_ancestor_level(p, leaf2))
            redundant_edge[q - wstride_d] = 1;
    }

    // Cleanup
    for (squirrel = leaf1; _node[squirrel].parentIdx != squirrel; squirrel = _node[squirrel].parentIdx)
        _node[squirrel]._rootIdx = ROOTIDX;
    _node[squirrel]._rootIdx = ROOTIDX;
}

template <class Pixel> void AlphaTree<Pixel>::HybridParallel(const Pixel *img, int numthreads) {
    // call clear
    if (_node)
        Free(_node);
    ImgIdx imgSize, dimgSize, nredges;
    RankItem<double> *rankitem;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    rankitem = (RankItem<double> *)Malloc(nredges * sizeof(RankItem<double>));

    ImgIdx *indexToRank = (ImgIdx *)Calloc((size_t)dimgSize * sizeof(ImgIdx));
    uint8_t *qrank = (uint8_t *)Malloc((size_t)dimgSize);
    _parentAry = (ImgIdx *)Calloc((size_t)imgSize * sizeof(ImgIdx));

    int64_t numpartitions;
    if (numthreads > 2)
        numpartitions = _min(imgSize / 2, _min(256, numthreads * 4));
    else
        numpartitions = 2;
    // int64_t numpartitions = numthreads;
    int npartition_x, npartition_y;
    {
        int optpart = 1;
        double optborderlength = (double)numpartitions * (double)imgSize;
        for (int px = 2; px < numpartitions; px++) {
            if (numpartitions % px == 0) {
                int py = numpartitions / px;

                if (((double)px * (double)_height + (double)py * (double)_width) < optborderlength) {
                    optpart = px;
                    optborderlength = ((double)px * (double)_height + (double)py * (double)_width);
                }
            }
        }
        npartition_x = (int)optpart;
        npartition_y = (int)numpartitions / npartition_x;
    }

    int32_t numbins = numpartitions; // number of levels in Quantization (= number of threads used)

    uint8_t *isVisited, *isAvailable;
    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable_par(isAvailable, npartition_x, npartition_y);

    ImgIdx binsize = nredges / (int64_t)numbins;
    numbins = (nredges + binsize - 1) / binsize;
    const ImgIdx numlevels = numbins; // for compatibility
    const ImgIdx blksz_x = _width / npartition_x;
    const ImgIdx blksz_y = _height / npartition_y;
    const ImgIdx blksz_xn = blksz_x + (_width % npartition_x);
    const ImgIdx blksz_yn = blksz_y + (_height % npartition_y);
    numpartitions = (int64_t)npartition_x * (int64_t)npartition_y;

    ImgIdx *startpidx = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *blocksize = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *subtree_size = (ImgIdx *)Calloc((numpartitions) * sizeof(ImgIdx));
    ImgIdx *subtree_start = (ImgIdx *)Calloc((numpartitions + 1) * sizeof(ImgIdx));
    ImgIdx *subtree_nborderedges = (ImgIdx *)Calloc((numpartitions) * sizeof(ImgIdx));
    ImgIdx *subtree_cur = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *subtree_max = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *blockWidths = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *blockHeights = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *dhist = (ImgIdx *)Calloc((size_t)numbins * (size_t)numpartitions * sizeof(ImgIdx));

    double *nrmsds = (double *)Malloc(numpartitions * sizeof(double));

    omp_set_num_threads(_min(numthreads, omp_get_num_procs()));

    int32_t *rankToIndex = (int32_t *)Calloc(nredges * sizeof(int32_t));

    double t0 = get_wall_time();

    // HP
    computeDifferenceAndSort(indexToRank, rankitem, img, nredges, rankToIndex);

    double t1 = get_wall_time();

    printf("computeDifferenceAndSort = %f\n", t1 - t0);

    set_subblock_properties(startpidx, blockWidths, blockHeights, blocksize, npartition_x, npartition_y, blksz_x,
                            blksz_y, blksz_xn, blksz_yn);
    quantize_ranks_compute_histogram(qrank, indexToRank, img, dhist, blockWidths, blockHeights, startpidx, binsize,
                                     numbins, npartition_x, npartition_y, subtree_max);

    // subtree_size allocation (no TSE)
    for (int blk = 0; blk < numpartitions; blk++)
        subtree_size[blk] = blockWidths[blk] * blockHeights[blk] * _connectivity;

    HierarQueue **queues;
    memalloc_queues(&queues, numpartitions, blocksize, subtree_max);

    ImgIdx *hypernode_level = (ImgIdx *)Malloc(dimgSize * sizeof(ImgIdx));
#pragma omp parallel for
    for (ImgIdx i = 0; i < dimgSize; i++)
        hypernode_level[i] = ROOTIDX;

    int8_t *levelroots = (int8_t *)Calloc(numbins * numpartitions * sizeof(int8_t));
    int8_t *redundant_edge = (int8_t *)Calloc(dimgSize * sizeof(int8_t));

    double treesizemult_intv = 0.05;
    double treesizemult = 1.0 - treesizemult_intv;

    {
        ImgIdx shamt = _connectivity >> 2;
        ImgIdx wstride_d = _width << shamt;

        treesizemult = treesizemult + treesizemult_intv;

        subtree_start[0] = 0;
        for (int blk = 0; blk < numpartitions; blk++)
            subtree_start[blk + 1] = subtree_size[blk] + subtree_start[blk];

        _maxSize = subtree_start[numpartitions];
        _node = (AlphaNode<Pixel> *)Calloc((size_t)(_maxSize) * sizeof(AlphaNode<Pixel>));

#pragma omp parallel for private(p, q) schedule(dynamic, 1)
        for (int blk = 0; blk < numpartitions; blk++) {
            const ImgIdx blockArea = blockWidths[blk] * blockHeights[blk];
            ImgIdx *blockDiffHist = dhist + numlevels * blk;
            HierarQueue *queue = queues[blk];
            ImgIdx nodeIdx = subtree_start[blk];
            Pixel maxdiff = subtree_max[blk];

            // printf("th%d: subtree for blk %d: setting queue\n", omp_get_thread_num(), (int)blk);
            queue->set_queue(blockDiffHist);
            // queue->set_queue(blockDiffHist, maxdiff);

            ImgIdx stackTop = nodeIdx++;
            ImgIdx prevTop = stackTop;
            Pixel currentLevel = maxdiff;
            _node[stackTop] = AlphaNode<Pixel>((double)currentLevel, stackTop);

            ImgIdx x0 = startpidx[blk]; /*starting point*/
            queue->push((x0 << shamt) << 1, currentLevel);
            prevTop = stackTop; /*to find redundant _node*/
            bool isFirstPixel = true;

            while (_node[stackTop].area < blockArea) // flooding
            {
                while ((int64_t)queue->min_level <= (int64_t)currentLevel) // flood all levels below currentLevel
                {
                    auto [p, isAv] = computeFloodingContextHybridPartition(isFirstPixel, queue->pop(), isAvailable);

                    if (isVisited[p]) {
                        queue->find_minlev();
                        continue;
                    }
                    isVisited[p] = 1;
                    uint8_t connected_neighbor = 0;
                    if (_connectivity == 4) {
                        q = p << shamt;
                        if (is_available(isAv, 0)) {
                            if (isVisited[p + _width])
                                connected_neighbor |= 0x1;
                            else
                                queue->push((q << 1) | 1, qrank[q]);
                        }
                        if (is_available(isAv, 1)) {
                            if (isVisited[p + 1])
                                connected_neighbor |= 0x2;
                            else
                                queue->push(((q + 1) << 1) | 1, qrank[q + 1]);
                        }
                        if (is_available(isAv, 2)) {
                            if (isVisited[p - 1])
                                connected_neighbor |= 0x4;
                            else
                                queue->push(((q - 1) << 1), qrank[q - 1]);
                        }
                        if (is_available(isAv, 3)) {
                            if (isVisited[p - _width])
                                connected_neighbor |= 0x8;
                            else
                                queue->push(((q - wstride_d) << 1), qrank[q - wstride_d]);
                        }
                    } else { // 8N
                    }

                    const auto pixVal = img[p];
                    if ((int64_t)currentLevel > (int64_t)queue->min_level) { // Go to lower level
                        currentLevel = queue->min_level;

                        const ImgIdx newStackTop = nodeIdx++;
                        _node[newStackTop] = AlphaNode<Pixel>(pixVal, (double)currentLevel, stackTop);
                        stackTop = newStackTop;

                        if (currentLevel) {
                            const ImgIdx singletonNodeIdx = nodeIdx++;
                            _node[singletonNodeIdx] = AlphaNode<Pixel>(pixVal, 0.0, stackTop);
                            prevTop = singletonNodeIdx;
                            _parentAry[p] = singletonNodeIdx;
                        } else
                            _parentAry[p] = stackTop;
                    } else {
                        queue->find_minlev();

                        if (currentLevel) {
                            const ImgIdx singletonNodeIdx = nodeIdx++;
                            _node[singletonNodeIdx] = AlphaNode<Pixel>(pixVal, 0.0, stackTop);
                            _node[stackTop].add(_node[singletonNodeIdx]);
                            _parentAry[p] = singletonNodeIdx;
                        } else {
                            _parentAry[p] = stackTop;
                            _node[stackTop].add(pixVal);
                        }
                    }

                    if (connected_neighbor)
                        markRedundantEdges(redundant_edge, p, connected_neighbor, qrank);
                }

                // Remove redundant node
                if (_node[prevTop].parentIdx == stackTop && _node[prevTop].area == _node[stackTop].area) {
                    _node[prevTop].parentIdx = _node[stackTop].parentIdx;
                    stackTop = prevTop;
                }

                if (_node[stackTop].area < blockArea) {
                    // go to higher level
                    ImgIdx nextStackTop = _node[stackTop].parentIdx;
                    if ((int64_t)queue->min_level < (int64_t)_node[nextStackTop].alpha) { // new level from queue
                        const ImgIdx queueTop = nodeIdx++;

                        _node[queueTop] = AlphaNode<Pixel>(queue->min_level, _node[stackTop].parentIdx);
                        _node[queueTop].add(_node[stackTop]);
                        _node[stackTop].parentIdx = queueTop;
                        nextStackTop = queueTop;
                    } else
                        _node[nextStackTop].add(_node + stackTop);

                    prevTop = stackTop;
                    stackTop = nextStackTop;
                    currentLevel = _node[stackTop].alpha;
                }
            }
            _node[stackTop].parentIdx = ROOTIDX;
            subtree_cur[blk] = nodeIdx;
        }
    }

    if (numpartitions > 1)
        _rootIdx = mergeSubtreeHybrid(qrank, blksz_x, blksz_y, npartition_x, npartition_y, subtree_cur);
    else {
        _rootIdx = _parentAry[0];
        while (_node[_rootIdx].parentIdx != ROOTIDX)
            _rootIdx = _node[_rootIdx].parentIdx;
    }

    ////////////////////////////////////////////////////////////////////////
    // Initialize refined tree
    ////////////////////////////////////////////////////////////////////////

#pragma omp parallel for
    for (p = 0; p < _maxSize; p++)
        _node[p]._rootIdx = ROOTIDX;

    _maxSize = imgSize + nredges;
    AlphaNode<Pixel> *pilottree = _node;
    _node = (AlphaNode<Pixel> *)Malloc((size_t)(_maxSize + 1) * sizeof(AlphaNode<Pixel>));
    _nodeIn = _node + imgSize;

    Pixel maxpixval = std::numeric_limits<Pixel>::max();
    initialize_node_par1(img, rankitem, maxpixval, rankToIndex);

    double *levelperthread = (double *)Calloc((numbins + 1) * sizeof(double));
    double *timeperthread = (double *)Calloc((numbins + 1) * sizeof(double));

#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t qlevel = 0; qlevel < numbins; qlevel++) // this part is to be parallelised
    {
        if (qlevel > (int64_t)pilottree[_rootIdx].alpha) {
            continue;
        }
        unionfind_refine_qlevel(qlevel, binsize, nredges, pilottree, rankitem, redundant_edge, rankToIndex);

        levelperthread[omp_get_thread_num()]++;
    }

    connect_pilotnode(pilottree, nredges, imgSize);

    _node[_rootIdx].parentIdx = ROOTIDX;

    // printTree();
    // std::getchar();

    Free(_parentAry);
    _parentAry = nullptr;
    Free(levelperthread);
    Free(timeperthread);
    if (rankToIndex)
        Free(rankToIndex);

    for (int blk = 0; blk < numpartitions; blk++)
        delete queues[blk];
    Free(queues);
    Free(nrmsds);
    Free(blockWidths);
    Free(blockHeights);
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
    Free(indexToRank);
    Free(rankitem);
    Free(isVisited);
    Free(isAvailable);
    Free(qrank);
}

template <class Pixel> void AlphaTree<Pixel>::HybridParallelOld(const Pixel *img, int numthreads) {
    ImgIdx imgSize, dimgSize, nredges;
    RankItem<double> *rankitem;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    rankitem = (RankItem<double> *)Malloc(nredges * sizeof(RankItem<double>));

    ImgIdx *rank = (ImgIdx *)Calloc((size_t)dimgSize * sizeof(ImgIdx));
    uint8_t *qrank = (uint8_t *)Malloc((size_t)dimgSize);
    _parentAry = (ImgIdx *)Calloc((size_t)imgSize * sizeof(ImgIdx));

    int64_t numpartitions;
    if (numthreads > 2)
        numpartitions = _min(imgSize / 2, _min(256, numthreads * 4));
    else
        numpartitions = 2;
    // int64_t numpartitions = numthreads;
    int npartition_x, npartition_y;
    {
        int optpart = 1;
        double optborderlength = (double)numpartitions * (double)imgSize;
        for (int px = 2; px < numpartitions; px++) {
            if (numpartitions % px == 0) {
                int py = numpartitions / px;

                if (((double)px * (double)_height + (double)py * (double)_width) < optborderlength) {
                    optpart = px;
                    optborderlength = ((double)px * (double)_height + (double)py * (double)_width);
                }
            }
        }
        npartition_x = (int)optpart;
        npartition_y = (int)numpartitions / npartition_x;
    }

    int32_t numbins = numpartitions; // number of levels in Quantization (= number of threads used)

    uint8_t *isVisited, *isAvailable;
    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));
    set_isAvailable_par(isAvailable, npartition_x, npartition_y);

    ImgIdx binsize = nredges / (int64_t)numbins;
    numbins = (nredges + binsize - 1) / binsize;
    int64_t numlevels = numbins; // for compatibility
    int64_t blksz_x = _width / npartition_x;
    int64_t blksz_y = _height / npartition_y;
    int64_t blksz_xn = blksz_x + (_width % npartition_x);
    int64_t blksz_yn = blksz_y + (_height % npartition_y);
    numpartitions = (int64_t)npartition_x * (int64_t)npartition_y;

    ImgIdx *startpidx = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *blocksize = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *subtree_size = (ImgIdx *)Calloc((numpartitions) * sizeof(ImgIdx));
    ImgIdx *subtree_start = (ImgIdx *)Calloc((numpartitions + 1) * sizeof(ImgIdx));
    ImgIdx *subtree_nborderedges = (ImgIdx *)Calloc((numpartitions) * sizeof(ImgIdx));
    ImgIdx *subtree_cur = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *subtree_max = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *blockWidths = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *blockHeights = (ImgIdx *)Malloc(numpartitions * sizeof(ImgIdx));
    ImgIdx *dhist = (ImgIdx *)Calloc((size_t)numbins * (size_t)numpartitions * sizeof(ImgIdx));
    char *blkflooddone = (char *)Calloc(numpartitions * sizeof(char));
    omp_lock_t *locks = (omp_lock_t *)Malloc((numpartitions + 1) * sizeof(omp_lock_t));

    for (p = 0; p < numpartitions; p++)
        omp_init_lock(locks + p);

    double *nrmsds = (double *)Malloc(numpartitions * sizeof(double));

    omp_set_num_threads(_min(numthreads, omp_get_num_procs()));

    int32_t *rankToIndex = (int32_t *)Calloc(nredges * sizeof(int32_t));

    computeDifferenceAndSort(rank, rankitem, img, nredges, rankToIndex);

    // printGraph(isVisited, (uint8_t *)nullptr, img);
    // for (int r = 0; r < nredges; r++) {
    //     printf("Rank[%d]: index = %d, ", r, rankToIndex[r]);
    //     rankitem[rankToIndex[r]].print();
    //     printf("\n");
    // }
    // return;

    set_subblock_properties(startpidx, blockWidths, blockHeights, blocksize, npartition_x, npartition_y, blksz_x,
                            blksz_y, blksz_xn, blksz_yn);
    quantize_ranks_compute_histogram(qrank, rank, img, dhist, blockWidths, blockHeights, startpidx, binsize, numbins,
                                     npartition_x, npartition_y, subtree_max);

    blockwise_tse(subtree_size, subtree_nborderedges, nrmsds, dhist, subtree_max, blockWidths, blockHeights,
                  npartition_x, npartition_y, numbins);

    HierarQueue **queues;
    memalloc_queues(&queues, numpartitions, blocksize, subtree_max);

    ImgIdx *hypernode_level = (ImgIdx *)Malloc(dimgSize * sizeof(ImgIdx));
    int8_t *levelroots = (int8_t *)Calloc(numbins * numpartitions * sizeof(int8_t));
    int8_t *redundant_edge = (int8_t *)Calloc(dimgSize * sizeof(int8_t));

    double treesizemult_intv = 0.05;
    double treesizemult = 1.0 - treesizemult_intv;

    int flooddone = 0;
    while (!flooddone) {
        int numbusythr = 0;
        int outofmemory = 0;
        int numblkproc = 0;

        ImgIdx shamt = _connectivity >> 2;
        ImgIdx wstride_d = _width << shamt;

        // reset queue, isvisited array, hypernode levels
        for (int blk = 0; blk < numpartitions; blk++)
            queues[blk]->reset_queue();

#pragma omp parallel for private(p, q)
        for (int i = 0; i < imgSize; i++)
            isVisited[i] = 0;

#pragma omp parallel for private(p, q)
        for (ImgIdx i = 0; i < dimgSize; i++) {
            hypernode_level[i] = ROOTIDX;
            redundant_edge[i] = 0;
        }

        flooddone = 1; // reset this flag when memory overflows on all threads
        treesizemult = treesizemult + treesizemult_intv;

        //(re)allocate _node array (expand size when reallocate)
        _maxSize =
            parflood_node_alloc(subtree_size, subtree_start, blockWidths, blockHeights, numpartitions, treesizemult);

        for (p = 0; p < numpartitions; p++)
            omp_unset_lock(locks + p);
        omp_unset_lock(locks + numpartitions);
        numbusythr = 0;

#pragma omp parallel for private(p, q) schedule(dynamic, 1)
        for (int blk = 0; blk < numpartitions; blk++) // flooding is somehow slower than the one without tse...
        {
            if (outofmemory)
                continue;

            omp_set_lock(locks + blk);
            omp_set_lock(locks + numpartitions);
            numbusythr++;
            numblkproc++;
            omp_unset_lock(locks + numpartitions);

            uint8_t connected_neighbor; // marks neighbors that are already visited
            ImgIdx lsbclearmask = ~1;   // mask for clearing 1st bit
            ImgIdx blockWidth = blockWidths[blk];
            ImgIdx blockHeight = blockHeights[blk];
            // ImgIdx blksize = blocksize[blk];
            ImgIdx blockArea = blockWidth * blockHeight;
            ImgIdx *blockDiffHist = dhist + numlevels * blk;
            HierarQueue *queue = queues[blk];
            // bool lastcol = (blk % npartition_x) == (npartition_x - 1);
            // bool lastrow = (blk / npartition_x) == (npartition_y - 1);
            // ImgIdx startPixelIndex = startpidx[blk];
            ImgIdx nidx = subtree_start[blk];
            ImgIdx blkts = 0;
            int nidxblk = blk;
            ImgIdx nidx_lim =
                subtree_start[blk + 1] - subtree_nborderedges[blk + 1]; // save room for nodes to be added in merge
            ImgIdx iNode = 0;
            Pixel maxdiff = subtree_max[blk];

            int8_t *plr = levelroots + blk * numbins;
            for (p = 0; p < numbins; p++)
                plr[p] = 0;

            // printf("th%d: subtree for blk %d: setting queue\n", omp_get_thread_num(), (int)blk);
            queue->set_queue(blockDiffHist);
            // queue->set_queue(blockDiffHist, maxdiff);

            ImgIdx stackTop = nidx++; // imgSize + dimgSize + blk;
            // printf("blk%d - Dummy %d\n", (int)blk, (int)(stackTop));
            ImgIdx prevTop = stackTop;
            AlphaNode<Pixel> *pNode = _node + stackTop;
            pNode->set(0, maxdiff, (double)0.0, (Pixel)-1, (Pixel)0);
            pNode->parentIdx = stackTop;
            pNode->_rootIdx = ROOTIDX;
            Pixel currentLevel = maxdiff;

            ImgIdx x0 = startpidx[blk]; /*starting point*/
            queue->push(((x0 << shamt) << 1) & lsbclearmask, currentLevel);
            prevTop = stackTop; /*to find redundant _node*/
            int firstpix = 1;

            if (outofmemory)
                continue;
            while (1) // flooding
            {
                while ((int64_t)queue->min_level <= (int64_t)currentLevel) // flood all levels below currentLevel
                {
                    ImgIdx qitem = queue->pop();
                    ImgIdx didx = qitem >> 1;
                    // the pixel which pushed this item into the queue, where it is easier to find the level of the edge
                    // (level = stackTop)

                    p = didx >> shamt;

                    uint8_t isAv;
                    if (firstpix) {
                        isAv = isAvailable[p];
                        firstpix = 0;
                    } else {
                        if ((qitem & 1)) // pixel at another end of the edge?
                        {
                            if (_connectivity == 4) {
                                if (didx & 1) // horizontal edge
                                {
                                    p++;
                                    isAv = isAvailable[p];
                                    isAv &= ~0x4; // do not check the neighbor which pushed this pixel (p) into the
                                                  // queue
                                } else            // vertical edge
                                {
                                    p += _width;
                                    isAv = isAvailable[p];
                                    isAv &= ~0x8;
                                }
                            } else {
                                isAv = 0;
                            }
                        } else {
                            if (didx & 1) // horizontal edge
                            {
                                isAv = isAvailable[p];
                                isAv &= ~0x2;
                            } else // vertical edge
                            {
                                isAv = isAvailable[p];
                                isAv &= ~0x1;
                            }
                        }
                    }

                    // printf("thr%d: probing %d\n", omp_get_thread_num(), (int)p);
                    if (isVisited[p]) {
                        queue->find_minlev();
                        continue;
                    }

                    isVisited[p] = 1;
                    connected_neighbor = 0;
                    if (_connectivity == 4) {
                        q = p << shamt;
                        ImgIdx q1;
                        if (is_available(isAv, 0)) {
                            if (isVisited[p + _width]) // neighbor alread visited - which means this edge might be on
                                                       // lower level than the alpha value of the edge corresponds to
                            {
                                // if (!plr[qrank[q]]) redundant_edge[q] = 1;
                                // else
                                connected_neighbor |= 0x1;
                            } else // push new neighbor
                                queue->push((q << 1) | 1, qrank[q]);
                        }
                        if (is_available(isAv, 1)) {
                            q1 = q + 1;
                            if (isVisited[p + 1]) {
                                // if (!plr[qrank[q1]]) redundant_edge[q1] = 1;
                                // else
                                connected_neighbor |= 0x2;
                            } else
                                queue->push(((q1) << 1) | 1, qrank[q1]);
                        }
                        if (is_available(isAv, 2)) {
                            q1 = q - 1;
                            if (isVisited[p - 1]) {
                                // if (!plr[qrank[q1]]) redundant_edge[q1] = 1;
                                // else
                                connected_neighbor |= 0x4;
                            } else
                                queue->push(((q1) << 1) & lsbclearmask, qrank[q1]);
                        }

                        if (is_available(isAv, 3)) {
                            q1 = q - wstride_d;
                            if (isVisited[p - _width]) {
                                // if (!plr[qrank[q1]]) redundant_edge[q1] = 1;
                                // else
                                connected_neighbor |= 0x8;
                            } else
                                queue->push(((q1) << 1) & lsbclearmask, qrank[q1]);
                        }
                    }

                    if ((int64_t)currentLevel > (int64_t)queue->min_level) // go to lower level
                    {
                        // plr[queue->min_level] = 1;
                        Pixel pix_val = img[p];
                        currentLevel = queue->min_level;

                        {
                            if (nidx == nidx_lim) {
                                if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts, blkflooddone,
                                                     subtree_cur, subtree_start, subtree_nborderedges, locks,
                                                     numbusythr, numblkproc, outofmemory))
                                    break;
                            }
                            iNode = nidx++;
                        }
                        _node[iNode].set(1, currentLevel, (double)pix_val, pix_val, pix_val);
                        _node[iNode].parentIdx = stackTop;
                        _node[iNode]._rootIdx = ROOTIDX;
                        stackTop = iNode;

                        if (currentLevel) {
                            {
                                if (nidx == nidx_lim) {
                                    if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,
                                                         blkflooddone, subtree_cur, subtree_start, subtree_nborderedges,
                                                         locks, numbusythr, numblkproc, outofmemory))
                                        break;
                                }
                                iNode = nidx++;
                            }
                            _node[iNode].copy(_node + stackTop);
                            _node[iNode].alpha = 0;
                            _node[iNode].parentIdx = stackTop;
                            _node[iNode]._rootIdx = ROOTIDX;
                            prevTop = iNode;
                        }
                        _parentAry[p] = iNode;
                    } else {
                        queue->find_minlev();

                        if (currentLevel) {
                            Pixel pix_val = img[p];
                            {
                                if (nidx == nidx_lim) {
                                    if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts,
                                                         blkflooddone, subtree_cur, subtree_start, subtree_nborderedges,
                                                         locks, numbusythr, numblkproc, outofmemory))
                                        break;
                                }
                                iNode = nidx++;
                            }
                            _node[iNode].set(1, 0, (double)pix_val, pix_val, pix_val);
                            _node[stackTop].add(_node + iNode);
                            _node[iNode].parentIdx = stackTop;
                            _node[iNode]._rootIdx = ROOTIDX;
                            _parentAry[p] = iNode;
                        } else {
                            _parentAry[p] = stackTop;
                            _node[stackTop].add(img[p]);
                        }
                    }
                    // if (stackTop == ROOTIDX) printf("SQEEEEAK\n");

                    if (connected_neighbor) {
                        // mark from the leaf to the stackTop _node to help finding hypernode levelroots
                        ImgIdx squirrel, leaf1 = _parentAry[p], leaf2;
                        for (squirrel = leaf1; _node[squirrel].parentIdx != squirrel;
                             squirrel = _node[squirrel].parentIdx)
                            _node[squirrel]._rootIdx = p;
                        _node[squirrel]._rootIdx = p;

                        q = p << shamt;
                        if (connected_neighbor & 0x1) {
                            leaf2 = _parentAry[p + _width];
                            if (qrank[q] > get_nearest_common_ancestor_level(p, leaf2))
                                redundant_edge[q] = 1;
                        }
                        if (connected_neighbor & 0x2) {
                            leaf2 = _parentAry[p + 1];
                            // printf("q = %d\n", (int)q);
                            if (qrank[q + 1] > get_nearest_common_ancestor_level(p, leaf2))
                                redundant_edge[q + 1] = 1;
                        }
                        if (connected_neighbor & 0x4) {
                            leaf2 = _parentAry[p - 1];
                            if (qrank[q - 1] > get_nearest_common_ancestor_level(p, leaf2))
                                redundant_edge[q - 1] = 1;
                        }
                        if (connected_neighbor & 0x8) {
                            leaf2 = _parentAry[p - _width];
                            if (qrank[q - wstride_d] > get_nearest_common_ancestor_level(p, leaf2))
                                redundant_edge[q - wstride_d] = 1;
                        }

                        // cleanup pawprints
                        for (squirrel = leaf1; _node[squirrel].parentIdx != squirrel;
                             squirrel = _node[squirrel].parentIdx)
                            _node[squirrel]._rootIdx = ROOTIDX;
                        _node[squirrel]._rootIdx = ROOTIDX;
                    }
                }

                if (outofmemory)
                    break;

                // remove_redundant_node(_node, nidx, prevTop, stackTop);
                if (_node[prevTop].parentIdx == stackTop && _node[prevTop].area == _node[stackTop].area) {
                    // plr[(int)(_node[prevTop].alpha)] = 0;
                    _node[prevTop].parentIdx = _node[stackTop].parentIdx;
                    stackTop = prevTop;
                    // _curSize--;
                }

                if (_node[stackTop].area == blockArea) // root _node found...done
                    break;

                // go to higher level
                iNode = _node[stackTop].parentIdx;
                if ((int64_t)queue->min_level < (int64_t)_node[iNode].alpha) // new level from queue
                {
                    {
                        if (nidx == nidx_lim) {
                            if (!migrate_subtree(blk, numpartitions, nidx, nidx_lim, nidxblk, blkts, blkflooddone,
                                                 subtree_cur, subtree_start, subtree_nborderedges, locks, numbusythr,
                                                 numblkproc, outofmemory))
                                break;
                        }
                        iNode = nidx++;
                    }
                    _node[iNode].alpha = queue->min_level;
                    _node[iNode].copy(_node + stackTop);
                    _node[iNode].parentIdx = _node[stackTop].parentIdx;
                    _node[iNode]._rootIdx = ROOTIDX;
                    _node[stackTop].parentIdx = iNode;
                } else // go to existing _node
                {
                    if (_node[iNode].area == blockArea) // root _node found...done
                        break;
                    _node[iNode].add(_node + stackTop);
                }

                if (_node[iNode].area == blockArea) // root _node found...done
                    break;

                prevTop = stackTop;
                stackTop = iNode;
                currentLevel = _node[stackTop].alpha;
            }

            if (!outofmemory) {
                stackTop = (_node[stackTop].area == blockArea) ? stackTop : iNode; // remove redundant root
                _node[stackTop].parentIdx = ROOTIDX;

                subtree_cur[nidxblk] = nidx;

                blkts += nidx - subtree_start[nidxblk];

                if (nidx == nidx_lim) // this should be really rare
                {
                    blkflooddone[nidxblk] = 2; // flood done (at least for the native block), no free memory
                } else {
                    blkflooddone[nidxblk] = 1; // flood done AND free memory available
                }

                omp_set_lock(locks + numpartitions);
                numbusythr--;
                // printf("th%d-- (%d/%d)\n", omp_get_thread_num(), numbusythr, omp_get_num_threads());
                omp_unset_lock(locks + numpartitions);
                omp_unset_lock(locks + nidxblk);
            }

            if (outofmemory) {
                printf("thr%d: subtree for blk %d: memory overflow - releasing lock %d\n", omp_get_thread_num(),
                       (int)blk, (int)nidxblk);
                flooddone = 0;
            } else {
                // printf("thr%d: subtree for blk %d: releasing lock %d\n", omp_get_thread_num(), (int)blk,
                // (int)nidxblk);
            }
            // omp_unset_lock(locks + nidxblk);
        } // flood_end
    }

    if (numpartitions > 1)
        _rootIdx =
            merge_subtrees1(qrank, blksz_x, blksz_y, npartition_x, npartition_y, subtree_cur, 1, hypernode_level);
    else {
        _rootIdx = _parentAry[0];
        while (_node[_rootIdx].parentIdx != ROOTIDX)
            _rootIdx = _node[_rootIdx].parentIdx;
    }

    ////////////////////////////////////////////////////////////////////////
    // Initialize refined tree
    ////////////////////////////////////////////////////////////////////////

#pragma omp parallel for
    for (p = 0; p < _maxSize; p++)
        _node[p]._rootIdx = ROOTIDX;

    _maxSize = imgSize + nredges;
    AlphaNode<Pixel> *pilottree = _node;
    _node = (AlphaNode<Pixel> *)Malloc((size_t)(_maxSize + 1) * sizeof(AlphaNode<Pixel>));
    _nodeIn = _node + imgSize;

    Pixel maxpixval = std::numeric_limits<Pixel>::max();
    initialize_node_par1(img, rankitem, maxpixval, rankToIndex);

    double *levelperthread = (double *)Calloc((numbins + 1) * sizeof(double));
    double *timeperthread = (double *)Calloc((numbins + 1) * sizeof(double));

#pragma omp parallel for schedule(dynamic, 1)
    for (int64_t qlevel = 0; qlevel < numbins; qlevel++) // this part is to be parallelised
    {
        if (qlevel > (int64_t)pilottree[_rootIdx].alpha) {
            continue;
        }
        unionfind_refine_qlevel(qlevel, binsize, nredges, pilottree, rankitem, redundant_edge, rankToIndex);

        levelperthread[omp_get_thread_num()]++;
    }

    connect_pilotnode(pilottree, nredges, imgSize);
    Free(_parentAry);
    _parentAry = 0;
    _node[_rootIdx].parentIdx = ROOTIDX;

    Free(levelperthread);
    Free(timeperthread);
    if (rankToIndex)
        Free(rankToIndex);

    for (p = 0; p < numpartitions; p++)
        omp_destroy_lock(locks + p);
    Free(locks);
    for (int blk = 0; blk < numpartitions; blk++)
        delete queues[blk];
    Free(queues);
    Free(nrmsds);
    Free(blkflooddone);
    Free(blockWidths);
    Free(blockHeights);
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

template <class Pixel> ImgIdx AlphaTree<Pixel>::NewAlphaNode(ImgIdx &size, ImgIdx &maxsize) {
    if (size == maxsize) {
        // std::cout << "Reallocating...\n";
    }
    return size++;
}

template <class Pixel>
ImgIdx AlphaTree<Pixel>::NewAlphaNode(AlphaNode<Pixel> *tree, ImgIdx &size, ImgIdx &maxsize, Pixel level,
                                      AlphaNode<Pixel> *pCopy) {
    AlphaNode<Pixel> *pNew = tree + size;

    if (size == maxsize) {
        //	std::cout << "Reallocating...\n";
    }
    pNew->alpha = level;
    pNew->copy(pCopy);
    return size++;
}

template <class Pixel>
void AlphaTree<Pixel>::remove_redundant_node(AlphaNode<Pixel> *tree, ImgIdx &size, ImgIdx &prevTop, ImgIdx &stackTop) {
    if (tree[prevTop].parentIdx == stackTop && tree[prevTop].area == tree[stackTop].area) {
        tree[prevTop].parentIdx = tree[stackTop].parentIdx;
        stackTop = prevTop;
        size--;
    }
}

// for parallel pilottree
template <class Pixel>
void AlphaTree<Pixel>::connectPix2Node(AlphaNode<Pixel> *tree, ImgIdx pidx, Pixel pix_val, ImgIdx iNode, ImgIdx *pAry) {
    AlphaNode<Pixel> *pNode = &tree[iNode];
    pAry[pidx] = iNode;
    pNode->add(pix_val);
}

template <class Pixel>
ImgIdx AlphaTree<Pixel>::find_root1(ImgIdx p, ImgIdx qlevel) //, int &cnt)
{
    ImgIdx q, r;

    if (p == ROOTIDX)
        return ROOTIDX;

    //		if (!path_compression)
    //		{
    //		while((q = _node[p].parentIdx) != ROOTIDX)
    //			p = q;
    //			return p;
    //		}
    //		else
    {
        r = p;
        while ((q = _node[p]._rootIdx) != ROOTIDX)
            p = q;

        while (r != p) {
            q = _node[r]._rootIdx;
            _node[r]._rootIdx = p;
            // checkcoherence(_node + r, qlevel);
            // cnt++;
            r = q;
        }

        return p;
    }
}

template <class Pixel> void AlphaTree<Pixel>::FloodTrie(const Pixel *img) {
    ImgIdx imgSize, dimgSize, nredges;
    ImgIdx current_rank = 0, next_rank = 0;
    RankItem<double> *rankitem, *pRank;
    // AlphaNode<Pixel> *pNode;
    Pixel maxpixval;
    ImgIdx *rank, top_rank;
    int8_t nbits;
    // ImgIdx *dhist;
    ImgIdx prevTop = 0;
    uint8_t *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    _maxSize = imgSize + nredges;
    num_node = _maxSize;
    num_node_in = nredges;
    nbits = ((sizeof(Pixel) << 3) - 1);
    maxpixval = ~(1 << nbits);
    rankitem = (RankItem<double> *)Malloc(nredges * sizeof(RankItem<double>));
    _parentAry = 0;
    rank = (ImgIdx *)Malloc((size_t)dimgSize * sizeof(ImgIdx));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));
    _nodeIn = _node + imgSize;
    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));

    set_isAvailable(isAvailable);

    Trie_Cache *queue = new Trie_Cache(nredges);

    omp_set_num_threads(1);
    int32_t *rankToIndex = (int32_t *)Calloc(nredges * sizeof(int32_t));
    computeDifferenceAndSort(rank, rankitem, img, nredges, rankToIndex);

    initialize_node1(img, rankitem, maxpixval, rankToIndex);

    // manually visit the first pixel
    isVisited[0] = 1;
    if (_connectivity == 4) {
        queue->push(rank[0]);
        queue->push(rank[1]);
    } else if (_connectivity == 8) {
        queue->push(rank[0]);
        queue->push(rank[1]);
        queue->push(rank[2]);
    }
    // else later
    current_rank = queue->top();
    _node[0].connect_to_parent(&_nodeIn[current_rank], current_rank + imgSize);
    prevTop = current_rank;

    while (1) {
        while (1) {
            top_rank = queue->top();
            pRank = rankitem + rankToIndex[top_rank];
            if (isVisited[pRank->get_pidx0(_connectivity)]) {
                if (isVisited[pRank->get_pidx1(_width, _connectivity)])
                    break;
                p = pRank->get_pidx1(_width, _connectivity);
            } else
                p = pRank->get_pidx0(_connectivity);

            isVisited[p] = 1;
            isAv = isAvailable[p];
            if (_connectivity == 4) {
                q = p << 1;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(rank[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])
                    queue->push(rank[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])
                    queue->push(rank[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])
                    queue->push(rank[q - (_width << 1)]);
            } else if (_connectivity == 8) {
                ImgIdx width4 = _width << 2;
                q = p << 2;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(rank[q]); // printf("0:pushing %d \n",(int)rank[q]);}
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])
                    queue->push(rank[q + 1]); // printf("1:pushing %d \n",(int)rank[q+1]);}
                if (is_available(isAv, 2) && !isVisited[p + 1])
                    queue->push(rank[q + 2]); // printf("2:pushing %d \n",(int)rank[q+2]);}
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])
                    queue->push(rank[q + 3]); // printf("3:pushing %d \n",(int)rank[q+3]);}
                if (is_available(isAv, 4) && !isVisited[p - _width])
                    queue->push(rank[q - width4]); // printf("4:pushing %d \n",(int)rank[q-width4]);}
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])
                    queue->push(rank[q - width4 - 3]); // printf("5:pushing %d \n",(int)rank[q-width4-3]);}
                if (is_available(isAv, 6) && !isVisited[p - 1])
                    queue->push(rank[q - 2]); //  printf("6:pushing %d \n",(int)rank[q-2]);}
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])
                    queue->push(rank[q + width4 - 1]); // printf("7:pushing %d \n",(int)rank[q+width4-1]);}
            } else {
                //?
            }
            // else later

            next_rank = queue->top();
            _node[p].connect_to_parent(&_nodeIn[next_rank], next_rank + imgSize);
            if (current_rank == next_rank)
                break;
            current_rank = next_rank;
        }

        queue->pop();
        next_rank = queue->top();

        // remove redundant _node
        if (_nodeIn[prevTop].parentIdx == current_rank + imgSize && _nodeIn[prevTop].area == _nodeIn[current_rank].area)
            current_rank = prevTop;

        _nodeIn[current_rank].connect_to_parent(&_nodeIn[next_rank], next_rank + imgSize);
        if (_nodeIn[next_rank].area == imgSize)
            break;

        prevTop = current_rank;
        current_rank = next_rank;
    }

    _rootIdx = (_nodeIn[current_rank].area == imgSize) ? current_rank + imgSize : next_rank + imgSize;
    _node[_rootIdx].parentIdx = ROOTIDX;

    delete queue;
    Free(rankToIndex);
    Free(rank);
    Free(rankitem);
    Free(isVisited);
    Free(isAvailable);
}

template <class Pixel> void AlphaTree<Pixel>::FloodTrieNoCache(const Pixel *img) {
    ImgIdx imgSize, dimgSize, nredges;
    ImgIdx current_rank = 0, next_rank = 0;
    RankItem<double> *rankitem, *pRank;
    // AlphaNode<Pixel> *pNode;
    Pixel maxpixval;
    ImgIdx *rank, top_rank;
    int8_t nbits;
    // ImgIdx *dhist;
    ImgIdx prevTop = 0;
    uint8_t *isVisited, /**isVisited_edges,*/ *isAvailable, isAv;
    ImgIdx p, q;
    imgSize = _width * _height;
    nredges = _width * (_height - 1) + (_width - 1) * _height +
              ((_connectivity == 8) ? ((_width - 1) * (_height - 1) * 2) : 0);
    dimgSize = (_connectivity >> 1) * _width * _height;
    _maxSize = imgSize + nredges;
    num_node = _maxSize;
    num_node_in = nredges;
    nbits = ((sizeof(Pixel) << 3) - 1);
    maxpixval = ~(1 << nbits);
    rankitem = (RankItem<double> *)Malloc(nredges * sizeof(RankItem<double>));
    _parentAry = 0;
    rank = (ImgIdx *)Malloc((size_t)dimgSize * sizeof(ImgIdx));
    _node = (AlphaNode<Pixel> *)Malloc((size_t)_maxSize * sizeof(AlphaNode<Pixel>));
    _nodeIn = _node + imgSize;
    isVisited = (uint8_t *)Calloc((size_t)((imgSize)));
    isAvailable = (uint8_t *)Malloc((size_t)(imgSize));

    set_isAvailable(isAvailable);

    Trie<TrieIdx> *queue = new Trie<TrieIdx>(nredges);

    omp_set_num_threads(1);
    int32_t *rankToIndex = (int32_t *)Calloc(nredges * sizeof(int32_t));
    computeDifferenceAndSort(rank, rankitem, img, nredges, rankToIndex);

    initialize_node1(img, rankitem, maxpixval, rankToIndex);

    // manually visit the first pixel
    isVisited[0] = 1;
    if (_connectivity == 4) {
        queue->push(rank[0]);
        queue->push(rank[1]);
    } else if (_connectivity == 8) {
        queue->push(rank[0]);
        queue->push(rank[1]);
        queue->push(rank[2]);
    }
    // else later
    current_rank = queue->top();
    _node[0].connect_to_parent(&_nodeIn[current_rank], current_rank + imgSize);
    prevTop = current_rank;

    while (1) {
        while (1) {
            top_rank = queue->top();
            pRank = rankitem + rankToIndex[top_rank];
            if (isVisited[pRank->get_pidx0(_connectivity)]) {
                if (isVisited[pRank->get_pidx1(_width, _connectivity)])
                    break;
                p = pRank->get_pidx1(_width, _connectivity);
            } else
                p = pRank->get_pidx0(_connectivity);

            isVisited[p] = 1;
            isAv = isAvailable[p];
            if (_connectivity == 4) {
                q = p << 1;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(rank[q]);
                if (is_available(isAv, 1) && !isVisited[p + 1])
                    queue->push(rank[q + 1]);
                if (is_available(isAv, 2) && !isVisited[p - 1])
                    queue->push(rank[q - 1]);
                if (is_available(isAv, 3) && !isVisited[p - _width])
                    queue->push(rank[q - (_width << 1)]);
            } else if (_connectivity == 8) {
                ImgIdx width4 = _width << 2;
                q = p << 2;
                if (is_available(isAv, 0) && !isVisited[p + _width])
                    queue->push(rank[q]); // printf("0:pushing %d \n",(int)rank[q]);}
                if (is_available(isAv, 1) && !isVisited[p + _width + 1])
                    queue->push(rank[q + 1]); // printf("1:pushing %d \n",(int)rank[q+1]);}
                if (is_available(isAv, 2) && !isVisited[p + 1])
                    queue->push(rank[q + 2]); // printf("2:pushing %d \n",(int)rank[q+2]);}
                if (is_available(isAv, 3) && !isVisited[p - _width + 1])
                    queue->push(rank[q + 3]); // printf("3:pushing %d \n",(int)rank[q+3]);}
                if (is_available(isAv, 4) && !isVisited[p - _width])
                    queue->push(rank[q - width4]); // printf("4:pushing %d \n",(int)rank[q-width4]);}
                if (is_available(isAv, 5) && !isVisited[p - _width - 1])
                    queue->push(rank[q - width4 - 3]); // printf("5:pushing %d \n",(int)rank[q-width4-3]);}
                if (is_available(isAv, 6) && !isVisited[p - 1])
                    queue->push(rank[q - 2]); //  printf("6:pushing %d \n",(int)rank[q-2]);}
                if (is_available(isAv, 7) && !isVisited[p + _width - 1])
                    queue->push(rank[q + width4 - 1]); // printf("7:pushing %d \n",(int)rank[q+width4-1]);}
            } else {
                //?
            }
            // else later

            next_rank = queue->top();
            _node[p].connect_to_parent(&_nodeIn[next_rank], next_rank + imgSize);
            if (current_rank == next_rank)
                break;
            current_rank = next_rank;
        }

        queue->pop();
        next_rank = queue->top();

        // remove redundant _node
        if (_nodeIn[prevTop].parentIdx == current_rank + imgSize && _nodeIn[prevTop].area == _nodeIn[current_rank].area)
            current_rank = prevTop;

        _nodeIn[current_rank].connect_to_parent(&_nodeIn[next_rank], next_rank + imgSize);
        if (_nodeIn[next_rank].area == imgSize)
            break;

        prevTop = current_rank;
        current_rank = next_rank;
    }

    _rootIdx = (_nodeIn[current_rank].area == imgSize) ? current_rank + imgSize : next_rank + imgSize;
    _node[_rootIdx].parentIdx = ROOTIDX;

    delete queue;
    Free(rankToIndex);
    Free(rank);
    Free(rankitem);
    Free(isVisited);
    Free(isAvailable);
}

// the name is misleading. It doesn't find levelroot but the highest _node below nodeidx.
template <class Pixel> ImgIdx AlphaTree<Pixel>::get_level_root(ImgIdx p, ImgIdx nodeidx) {
    while (_node[p].parentIdx != ROOTIDX && nodeidx > _node[p].parentIdx)
        p = _node[p].parentIdx;

    return p;
}

template <class Pixel> ImgIdx AlphaTree<Pixel>::get_level_root(ImgIdx p) { return get_level_root(p, _node); }

template <class Pixel> ImgIdx AlphaTree<Pixel>::get_level_root(ImgIdx p, Pixel alpha) {
    while (_node[p].parentIdx != ROOTIDX && alpha >= _node[_node[p].parentIdx].alpha)
        p = _node[p].parentIdx;

    return p;
}

template <class Pixel> ImgIdx AlphaTree<Pixel>::get_level_root(ImgIdx p, AlphaNode<Pixel> *tree) {
    if (p == ROOTIDX)
        return ROOTIDX;
    Pixel a = tree[p].alpha;
    while (1) {
        ImgIdx parent = tree[p].parentIdx;
        if (parent == ROOTIDX || tree[parent].alpha > a)
            break;
        p = parent;
        // cntcnt++;
    }

    return p;
}

template <class Pixel> void AlphaTree<Pixel>::swap(ImgIdx &x, ImgIdx &y) {
    ImgIdx tmp = x;
    x = y;
    y = tmp;
}

template <class Pixel> void AlphaTree<Pixel>::swap(AlphaNode<Pixel> **x, AlphaNode<Pixel> **y) {
    AlphaNode<Pixel> *tmp = *x;
    *x = *y;
    *y = tmp;
}

template <class Pixel> ImgIdx AlphaTree<Pixel>::get_nearest_common_ancestor(ImgIdx x, ImgIdx y) {
    ImgIdx p;
    for (p = y; p != ROOTIDX; p = _node[p].parentIdx)
        _node[p]._rootIdx = ROOTIDX;
    for (p = x; p != ROOTIDX; p = _node[p].parentIdx)
        _node[p]._rootIdx = x;
    for (p = y; p != ROOTIDX && _node[p]._rootIdx != x; p = _node[p].parentIdx)
        ;
    return p;
}

template <class Pixel> Pixel AlphaTree<Pixel>::get_nearest_common_ancestor_level(ImgIdx x, ImgIdx y) {
    ImgIdx p;
    for (p = y; _node[p]._rootIdx != x; p = _node[p].parentIdx)
        ;
    return _node[p].alpha;
}

template <class Pixel> Pixel AlphaTree<Pixel>::connect(ImgIdx x, ImgIdx y, ImgIdx newidx, Pixel alpha) {
    ImgIdx x0, y0, z;
    //		bool compxy;
    ImgIdx imgSize = _height * _width;
    //		ImgIdx x1 = x, y1 = x;
    AlphaNode<Pixel> n0, n1, *p, *q;
    p = &n0;
    q = &n1;

    x = get_level_root(x, alpha);
    y = get_level_root(y, alpha);

    if (x == y) {
        return _node[x].alpha;
    }

    z = get_level_root(get_nearest_common_ancestor(x, y));

    _node[newidx].copy(_node + x);
    _node[newidx].add(_node + y);
    _node[newidx].alpha = alpha;
    _node[newidx].parentIdx = ROOTIDX;

    x0 = x;
    y0 = y;
    x = get_level_root(_node[x].parentIdx);
    y = get_level_root(_node[y].parentIdx);

    _node[x0].parentIdx = newidx;
    _node[y0].parentIdx = newidx;

    if (x == y) {
        if (x == ROOTIDX && y == ROOTIDX) {
            x = newidx;
        } else {
            _node[newidx].parentIdx = x;
        }
        return _node[x].alpha;
    }

    // y always has bigger alpha, or the same alpha but bigger area (for shorter path to level roots)
    if (x == ROOTIDX || (y != ROOTIDX && ((_node[x].alpha > _node[y].alpha) ||
                                          (_node[x].alpha == _node[y].alpha && _node[x].area > _node[y].area)))) {
        q->copy(_node + x0);
        swap(x, y);
    } else
        q->copy(_node + y0);

    _node[newidx].parentIdx = x;
    p->copy(_node + x);
    _node[x].add(q);

    while (1) {
        if (y == ROOTIDX) {
            while (_node[x].parentIdx != ROOTIDX) {
                x = (_node[x].parentIdx);
                _node[x].add(q);
                if (_node[x].area == imgSize) {
                    _node[x].parentIdx = ROOTIDX;
                    break;
                }
            }
            break;
        }

        if (y == z) // y is a common ancestor
        {

            if (x != y) {
                while (1) {
                    x = (_node[x].parentIdx);
                    if (x == y)
                        break;
                    _node[x].add(q);
                    if (_node[x].area == imgSize) {
                        _node[x].parentIdx = ROOTIDX;
                        break;
                    }
                }
            }
            break;
        }

        while (1) {
            x0 = get_level_root(_node[x].parentIdx);
            if (x0 != ROOTIDX && (_node[x0].alpha < _node[y].alpha)) {
                x = x0;
                x0 = get_level_root(_node[x0].parentIdx);
                p->copy(_node + x);
                _node[x].add(q);
                if (_node[x].area == imgSize) {
                    _node[x].parentIdx = ROOTIDX;
                    break;
                }
            } else
                break;
        }

        x0 = get_level_root(_node[x].parentIdx);
        _node[x].parentIdx = y;
        q->copy(_node + y);
        _node[y].add(p);

        if (_node[y].area == imgSize) {
            _node[y].parentIdx = ROOTIDX;
            break;
        }

        x = x0;
        swap(x, y);
        swap(&p, &q);
    }

    return alpha;
}

template <class Pixel> Pixel AlphaTree<Pixel>::connect(ImgIdx x, ImgIdx y, Pixel alpha, ImgIdx newidx) {
    ImgIdx x0, y0, z;
    //		bool compxy;
    ImgIdx imgSize = _height * _width;
    //		ImgIdx x1 = x, y1 = x;
    AlphaNode<Pixel> n0, n1, *p, *q;
    p = &n0;
    q = &n1;

    x = get_level_root(x, newidx);
    y = get_level_root(y, newidx);
    z = get_nearest_common_ancestor(x, y);

    if (x == y) {
        return _node[x].alpha;
    }

    _node[newidx].copy(_node + x);
    _node[newidx].add(_node + y);
    _node[newidx].alpha = alpha;
    _node[newidx].parentIdx = ROOTIDX;

    x0 = x;
    y0 = y;
    x = _node[x].parentIdx;
    y = _node[y].parentIdx;

    _node[x0].parentIdx = newidx;
    _node[y0].parentIdx = newidx;

    if (x == y) {
        if (x == ROOTIDX && y == ROOTIDX) {

        } else {
            _node[newidx].parentIdx = x;
        }
        return _node[x].alpha;
    }

    // compxy = _parentAry ? _node[x].alpha > _node[y].alpha : x > y;
    if (x == ROOTIDX || (y != ROOTIDX && (x > y))) {
        q->copy(_node + x0);
        swap(x, y);
    } else
        q->copy(_node + y0);

    _node[newidx].parentIdx = x;
    p->copy(_node + x);
    _node[x].add(q);

    while (1) {
        if (y == ROOTIDX) {
            while (_node[x].parentIdx != ROOTIDX) {
                x = _node[x].parentIdx;
                _node[x].add(q);
                if (_node[x].area == imgSize) {
                    _node[x].parentIdx = ROOTIDX;
                    break;
                }
            }
            break;
        }

        if (y == z) // y is a common ancestor
        {
            if (x != y) {
                while (1) {
                    x = _node[x].parentIdx;
                    if (x == y)
                        break;
                    _node[x].add(q);
                    if (_node[x].area == imgSize) {
                        _node[x].parentIdx = ROOTIDX;
                        break;
                    }
                }
            }
            break;
        }

        while (1) {
            x0 = _node[x].parentIdx;
            if (x0 != ROOTIDX && (x0 < y)) {
                x = x0;
                x0 = _node[x0].parentIdx;
                p->copy(_node + x);
                _node[x].add(q);
                if (_node[x].area == imgSize) {
                    _node[x].parentIdx = ROOTIDX;
                    break;
                }
            } else
                break;
        }

        x0 = _node[x].parentIdx;
        _node[x].parentIdx = y;
        q->copy(_node + y);
        _node[y].add(p);

        if (_node[y].area == imgSize) {
            _node[y].parentIdx = ROOTIDX;
            break;
        }

        x = x0;
        swap(x, y);
        swap(&p, &q);
    }

    return alpha;
}

template <class Pixel> void AlphaTree<Pixel>::canonicalize() {
    ImgIdx p;
    ImgIdx numcan = 0;
    ImgIdx imgSize = _height * _width;

    for (p = 0; p < imgSize; p++) {
        ImgIdx q = _parentAry[p];
        ImgIdx r = q;

        // Canonicalize leaf nodes
        if (r != ROOTIDX && _node[q].alpha == _node[_node[q].parentIdx].alpha) {
            while (_node[r].alpha == _node[_node[r].parentIdx].alpha) {
                _node[r].area = 0;
                r = _node[r].parentIdx;
            }
            _parentAry[p] = r;
            numcan++;
        }
    }

    for (p = _maxSize - 1; p >= 0; p--) {
        if (_node[p].area && _node[p].parentIdx != ROOTIDX && _node[_node[p].parentIdx].parentIdx != ROOTIDX &&
            _node[_node[p].parentIdx].alpha == _node[_node[_node[p].parentIdx].parentIdx].alpha) {
            ImgIdx q = _node[p].parentIdx;
            ImgIdx r = _node[q].parentIdx;
            _node[p].parentIdx = r;
            numcan++;
        }
    }

    for (p = _maxSize; p < _maxSize; p++) {
        ImgIdx q = p;
        if (p < imgSize)
            q = _parentAry[p];

        if (_node[q].area == 0 && (_node[p].parentIdx != ROOTIDX || _node[p]._rootIdx != ROOTIDX)) {
            _node[q].parentIdx = _node[q]._rootIdx = ROOTIDX;
        }
    }
}

template <class Pixel>
void AlphaTree<Pixel>::merge_subtrees(ImgIdx *rank, RankItem<Pixel> *rankitem, int64_t blksz_x, int64_t blksz_y,
                                      ImgIdx neighbor_offset, ImgIdx shamt, ImgIdx npartition_x, ImgIdx npartition_y) {
    ImgIdx imgSize = _height * _width, numblk;
    while (npartition_x > 1 || npartition_y > 1) {
        if ((npartition_x == 1 || blksz_x >= blksz_y) && npartition_y > 1) {
            numblk = npartition_x * (npartition_y / 2);
#pragma omp parallel for
            for (int blk = 0; blk < numblk; blk++) {
                ImgIdx x, y, r, p, p0, pn, dimgidx;
                y = (1 + 2 * (blk / (int)npartition_x)) * blksz_y;
                x = (blk % (int)npartition_x) * blksz_x;

                p0 = (y - 1) * _width + x;
                pn = (((blk % (int)npartition_x) == npartition_x - 1) ? y * _width : p0 + blksz_x);
                for (p = p0; p < pn; p++) {
                    dimgidx = p << 1;
                    r = rank[dimgidx];
                    canonicalize(p);
                    canonicalize(p + _width);
                    connect(p, p + _width, (Pixel)rankitem[r].alpha, (ImgIdx)(r + imgSize));
                }
            }
            npartition_y = (npartition_y + 1) / 2;
            blksz_y <<= 1;
            if (npartition_y == 1)
                blksz_y = _height;
            else
                blksz_y = _min(blksz_y, _height);
        }

        if ((npartition_y == 1 || blksz_x <= blksz_y) && npartition_x > 1) {
            numblk = npartition_y * (npartition_x / 2);

#pragma omp parallel for
            for (int blk = 0; blk < numblk; blk++) {
                ImgIdx x, y, r, p, p0, pn, dimgidx;
                x = (1 + 2 * (blk / npartition_y)) * blksz_x;
                y = (blk % (int)npartition_y) * blksz_y;

                p0 = y * _width + x - 1;
                pn = ((blk % (int)npartition_y) == npartition_y - 1) ? _height * _width : p0 + _width * blksz_y;

                for (p = p0; p < pn; p += _width) {
                    dimgidx = (p << 1) + 1;
                    r = rank[dimgidx];
                    canonicalize(p);
                    canonicalize(p + 1);
                    connect(p, p + 1, (Pixel)rankitem[r].alpha, (ImgIdx)(r + imgSize));
                }
            }
            npartition_x = (npartition_x + 1) / 2;
            blksz_x <<= 1;
            if (npartition_x == 1)
                blksz_x = _width;
            else
                blksz_x = _min(blksz_x, _width);
        }
    }
}

template <class Pixel>
void AlphaTree<Pixel>::set_subimgsizes(ImgIdx **subimgsizes, int8_t npartition_x, int8_t npartition_y, int64_t blksz,
                                       int64_t blksz_lastcol, int64_t blksz_lastrow, int64_t blksz_last) {
    *subimgsizes = (ImgIdx *)Malloc((int)npartition_x * (int)npartition_y * sizeof(ImgIdx));
    ImgIdx *sizes = *subimgsizes;
    ImgIdx p = 0;
    for (ImgIdx y = 0; y < npartition_y - 1; y++) {
        for (ImgIdx x = 0; x < npartition_x - 1; x++)
            sizes[p++] = blksz;
        sizes[p++] = blksz_lastcol;
    }
    for (ImgIdx x = 0; x < npartition_x - 1; x++)
        sizes[p++] = blksz_lastrow;
    sizes[p++] = blksz_last;
}

template class AlphaNode<uint8_t>;
template class AlphaNode<uint16_t>;
template class AlphaNode<uint32_t>;
template class AlphaNode<uint64_t>;

template class AlphaTree<uint8_t>;
template class AlphaTree<uint16_t>;
template class AlphaTree<uint32_t>;
template class AlphaTree<uint64_t>;
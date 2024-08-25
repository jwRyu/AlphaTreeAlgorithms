#pragma once
#include <QItem.hpp>
#include <defines.h>

template <class Pixel> class QuadHeapQueue {
    ImgIdx cursize;
    ImgIdx maxsize;
    QItem<Pixel> *arr;
    Pixel pop_level;

  public:
    ImgIdx get_cursize() { return cursize; }
    ImgIdx size() { return cursize; }
    ImgIdx sizeMax() { return maxsize; }
    bool empty() { return cursize == 0; }

    QuadHeapQueue(ImgIdx maxsize_in);
    ~QuadHeapQueue();

    void print();

    ImgIdx pop();
    void push(QItem<Pixel> item);
    inline const QItem<Pixel> &top() const { return arr[1]; }
};
#pragma once

#include <stddef.h>

extern size_t memuse, max_memuse;

void *Malloc(size_t size);
void *Calloc(size_t size);
// void* Malloc(size_t size, int fake);
void *Realloc(void *ptr, size_t size);
void Free(void *ptr);
// void Free(void* ptr, int fake);

#pragma once
#include <sig/gs_array.h>
#include <sig/gs_vec.h>

class CubeSurface
{
private:
    int cubeIndex;
    GsArray<GsVec> surfaceVerts;

public:
    CubeSurface(int cubeIndex, GsArray<GsVec> surfaceVerts)
    {
        this->cubeIndex = cubeIndex;
        this->surfaceVerts = surfaceVerts;
    }

    int getCubeIndex() { return cubeIndex; }
    GsArray<GsVec> getVerts() { return surfaceVerts; }

    ~CubeSurface() { }
};


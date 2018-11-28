#pragma once
#include <sig/gs_array.h>
#include <sig/gs_vec.h>

class CubeSurface
{
private:
    int cubeIndex;
    GsVec* surfaceVerts;
    GsArray<GsVec> surfaceNormals;

public:
    CubeSurface(int cubeIndex, GsVec* surfaceVerts, GsArray<GsVec> surfaceNormals)
    {
        this->cubeIndex = cubeIndex;
        this->surfaceVerts = surfaceVerts;
        this->surfaceNormals = surfaceNormals;
    }

    int getCubeIndex() { return cubeIndex; }
    GsVec* getVerts() { return surfaceVerts; }
    GsArray<GsVec> getNormals() { return surfaceNormals; }

    ~CubeSurface()
    {
        delete[] surfaceVerts;
    }
};


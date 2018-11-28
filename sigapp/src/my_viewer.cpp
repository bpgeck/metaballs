
#include "my_viewer.h"

MyViewer::MyViewer(int x, int y, int w, int h, const char* l) : WsViewer(x, y, w, h, l)
{
    build_ui();
    build_scene();
}

void MyViewer::build_ui()
{
    UiPanel *p;
    UiManager* uim = WsWindow::uim();
    p = uim->add_panel("", UiPanel::HorizLeft);
    
    p->add(new UiButton("Exit", EvExit));
    p->add(m_fidelitySlider = new UiSlider("fidelity: ", EvFidelity, 0, 0, 150));
    m_fidelitySlider->separate();
    m_fidelitySlider->range(1, 10);
    m_fidelitySlider->value(2.5f);

    p->add(m_thresholdSlider = new UiSlider("threshold: ", EvThreshold, 0, 0, 150));
    m_thresholdSlider->separate();
    m_thresholdSlider->range(1.5f, 0.05f);
    m_thresholdSlider->value(0.5f);

    p->top()->separate();
}

SnManipulator* MyViewer::add_model(SnShape* s, GsVec p)
{
    SnManipulator* manip = new SnManipulator;
    GsMat m;
    m.translation(p);
    manip->initial_mat(m);

    SnGroup* g = new SnGroup;
    SnLines* l = new SnLines;
    l->color(GsColor::orange);
    g->add(s);
    g->add(l);
    manip->child(g);

    rootg()->add(manip);
    return manip;
}

void MyViewer::build_scene()
{
    m_metaballCenters = new SnManipulator*[m_numCenters]; // holder for metaball centers to perform metaball calculations later on
    m_endpoints = new GsVec[m_numCenters];

    srand((int)gs_time());
    for (int i = 0; i < m_numCenters; i++)
    {
        SnPrimitive *center;
        center = new SnPrimitive(GsPrimitive::Box, 0.1f, 0.1f, 0.1f);
        center->prim().material.diffuse = GsColor::lightblue;

        // Place metaball centers randomly between the maximum and minimum world bounds for each dimension
        int maxX = m_worldBounds[0][1] * 2;
        int minX = m_worldBounds[0][0];
        int maxY = m_worldBounds[1][1] * 2;
        int minY = m_worldBounds[1][0];
        int maxZ = m_worldBounds[2][1] * 2;
        int minZ = m_worldBounds[2][0];

        GsVec point(rand() % maxX + minX, rand() % maxY + minY, rand() % maxZ + minZ);

        m_metaballCenters[i] = add_model(center, point);
        m_endpoints[i] = point;
    }

    camera().eye = GsVec(10, 10, 10);
    camera().center = GsVec(0, 0, 0);
}

void MyViewer::calculate_metaballs()
{
    // Use member variable m_subdivisionVal and m_worldBounds to subdivide the world into evenly distributed points
    float subdivVal = 1 / m_fidelitySlider->value(); // the higher the fidelity, the higher the quality
    int numCellsX = (int)((m_worldBounds[0][1] - m_worldBounds[0][0]) / subdivVal);
    int numCellsY = (int)((m_worldBounds[1][1] - m_worldBounds[1][0]) / subdivVal);
    int numCellsZ = (int)((m_worldBounds[2][1] - m_worldBounds[2][0]) / subdivVal);
    m_numCells[0] = numCellsX;
    m_numCells[1] = numCellsY;
    m_numCells[2] = numCellsZ;

    // Create the 3D array representing the metaball values in a subdivided 3D space
    m_metaballValues = new float**[numCellsX];
    for (int i = 0; i < m_numCells[0]; i++)
    {
        m_metaballValues[i] = new float*[numCellsY];
        for (int j = 0; j < m_numCells[1]; j++)
        {
            m_metaballValues[i][j] = new float[numCellsY];
            for (int k = 0; k < m_numCells[2]; k++)
            {
                m_metaballValues[i][j][k] = 0; // initialize all points to be 0
            }
        }
    }

    // Calcate the metaball value for each metaball at each world point (sum all metaball calculations at that point)
    float max = -1;
     m_threshold = m_thresholdSlider->value();
    for (int i = 0; i < m_numCells[0]; i++)
    {
        for (int j = 0; j < m_numCells[1]; j++)
        {
            for (int k = 0; k < m_numCells[2]; k++)
            {
                // convert the array index into world coordinates
                float worldX = ((float)i / m_numCells[0]) * (m_worldBounds[0][1] - m_worldBounds[0][0]) + m_worldBounds[0][0];
                float worldY = ((float)j / m_numCells[1]) * (m_worldBounds[1][1] - m_worldBounds[1][0]) + m_worldBounds[1][0];
                float worldZ = ((float)k / m_numCells[2]) * (m_worldBounds[2][1] - m_worldBounds[2][0]) + m_worldBounds[2][0];
                GsVec worldLoc(worldX, worldY, worldZ);

                // Calculate metaball value at the world point for each metaball
                float metaballVal = 0;
                for (int l = 0; l < m_numCenters; l++)
                {
                    float metaballX = m_metaballCenters[l]->mat().e14;
                    float metaballY = m_metaballCenters[l]->mat().e24;
                    float metaballZ = m_metaballCenters[l]->mat().e34;
                    GsVec metaballLoc(metaballX, metaballY, metaballZ);

                    metaballVal += 1 / (worldLoc - metaballLoc).norm2();
                }

                // Evaluate value at that point against the threshold. If greater, mark it as being "inside" of a metaball, by changing value to be greater than zero
                if (metaballVal > m_threshold)
                {
                    m_metaballValues[i][j][k] = metaballVal;

                    // Debugging: draw cubes on points considered "inside" metaball to get general idea of metaball shape before surface creation
                    /*
                    SnPrimitive* temp = new SnPrimitive(GsPrimitive::Box, 0.01f, 0.01f, 0.01f);
                    temp->prim().material.diffuse = GsColor::red;
                    add_model(temp, GsVec(worldX, worldY, worldZ));
                    */
                }
            }
        }
    }
}

GsVec MyViewer::lerp(GsVec a, GsVec b, float fa, float fb)
{
    GsVec p = a + ((b - a)*(m_threshold - fa)/(fb - fa));
    return p;
}

GsArray<GsVec> MyViewer::get_verts_from_cube(int cubeIndex, const GsArray<GsVec> cubeVerts, const GsArray<float> cubeVertMetaVals)
{
    int edgeCode = m_edgeTable[cubeIndex];

    if (edgeCode == 0)
        return GsArray<GsVec>();

    GsArray<GsVec> surfaceVerts(12, 12);

    // Check which edges are occupied, and get the linear interpolation of the vertices bordering that edge
    // All non-occupied edges will still have a value within surfaceVerts, it will just be an invalid value
    if (edgeCode & 1)
        surfaceVerts[0] = lerp(cubeVerts[0], cubeVerts[1], cubeVertMetaVals[0], cubeVertMetaVals[1]);
    if (edgeCode & 2)
        surfaceVerts[1] = lerp(cubeVerts[1], cubeVerts[2], cubeVertMetaVals[1], cubeVertMetaVals[2]);
    if (edgeCode & 4)
        surfaceVerts[2] = lerp(cubeVerts[2], cubeVerts[3], cubeVertMetaVals[2], cubeVertMetaVals[3]);
    if (edgeCode & 8)
        surfaceVerts[3] = lerp(cubeVerts[3], cubeVerts[0], cubeVertMetaVals[3], cubeVertMetaVals[0]);
    if (edgeCode & 16)
        surfaceVerts[4] = lerp(cubeVerts[4], cubeVerts[5], cubeVertMetaVals[4], cubeVertMetaVals[5]);
    if (edgeCode & 32)
        surfaceVerts[5] = lerp(cubeVerts[5], cubeVerts[6], cubeVertMetaVals[5], cubeVertMetaVals[6]);
    if (edgeCode & 64)
        surfaceVerts[6] = lerp(cubeVerts[6], cubeVerts[7], cubeVertMetaVals[6], cubeVertMetaVals[7]);
    if (edgeCode & 128)
        surfaceVerts[7] = lerp(cubeVerts[7], cubeVerts[4], cubeVertMetaVals[7], cubeVertMetaVals[4]);
    if (edgeCode & 256)
        surfaceVerts[8] = lerp(cubeVerts[0], cubeVerts[4], cubeVertMetaVals[0], cubeVertMetaVals[4]);
    if (edgeCode & 512)
        surfaceVerts[9] = lerp(cubeVerts[1], cubeVerts[5], cubeVertMetaVals[1], cubeVertMetaVals[5]);
    if (edgeCode & 1024)
        surfaceVerts[10] = lerp(cubeVerts[2], cubeVerts[6], cubeVertMetaVals[2], cubeVertMetaVals[6]);
    if (edgeCode & 2048)
        surfaceVerts[11] = lerp(cubeVerts[3], cubeVerts[7], cubeVertMetaVals[3], cubeVertMetaVals[7]);

    return surfaceVerts;
}

GsVec MyViewer::get_normal_from_vec(GsVec vert)
{
    GsVec normal(0, 0, 0);

    // Calculate the normal for each vertex at each mutaball surface vertex
    for (int j = 0; j < m_numCenters; j++)
    {
        float metaballX = m_metaballCenters[j]->mat().e14;
        float metaballY = m_metaballCenters[j]->mat().e24;
        float metaballZ = m_metaballCenters[j]->mat().e34;
        GsVec center(metaballX, metaballY, metaballZ);

        float denom = gs_pow(gs_pow(center.x - vert.x, 2) + gs_pow(center.y - vert.y, 2) + gs_pow(center.z - vert.z, 2), 2);
        GsVec temp(2 * (vert.x - center.x), 2 * (vert.y - center.y), 2 * (vert.z - center.z));
        temp /= denom;
        normal += temp;
    }

    return normal;
}

// Function used to calculate all surfaces for a cube, and a
void MyViewer::cube_thread(
    int x, int y, int z, 
    int minX, int maxX, int minY, int maxY, int minZ, int maxZ,
    int numCellsX, int numCellsY, int numCellsZ,
    GsArray<float> cubeMetaballValues)
{
    // Each cube has 8 vertices. These are set as occupied or unoccupied within m_metaballValues
    // We perform marching cubes on a case-by-case basis, where the occupied vertices corresponds to the specific case we use
    int cubeIndex = 0; // instead of having a giant switch statement, have each bit in cubeIndex correspond to a different corner (1 for occupied, 0 if not)
    if (cubeMetaballValues[0])
        cubeIndex |= 1;
    if (cubeMetaballValues[1])
        cubeIndex |= 2;
    if (cubeMetaballValues[2])
        cubeIndex |= 4;
    if (cubeMetaballValues[3])
        cubeIndex |= 8;
    if (cubeMetaballValues[4])
        cubeIndex |= 16;
    if (cubeMetaballValues[5])
        cubeIndex |= 32;
    if (cubeMetaballValues[6])
        cubeIndex |= 64;
    if (cubeMetaballValues[7])
        cubeIndex |= 128;

    // Calculate the world coordinates for each vertex of the cube
    GsArray<GsVec> cubeVerts;
    cubeVerts.push(GsVec(
        ((float)x / numCellsX) * (maxX - minX) + minX,
        ((float)y / numCellsY) * (maxY - minY) + minY,
        ((float)z / numCellsZ) * (maxZ - minZ) + minZ
    )); // vert 0

    cubeVerts.push(GsVec(
        ((float)(x + 1) / numCellsX) * (maxX - minX) + minX,
        ((float)y / numCellsY) * (maxY - minY) + minY,
        ((float)z / numCellsZ) * (maxZ - minZ) + minZ
    )); // vert 1

    cubeVerts.push(GsVec(
        ((float)(x + 1) / numCellsX) * (maxX - minX) + minX,
        ((float)y / numCellsY) * (maxY - minY) + minY,
        ((float)(z + 1) / numCellsZ) * (maxZ - minZ) + minZ
    )); // vert 2

    cubeVerts.push(GsVec(
        ((float)x / numCellsX) * (maxX - minX) + minX,
        ((float)y / numCellsY) * (maxY - minY) + minY,
        ((float)(z + 1) / numCellsZ) * (maxZ - minZ) + minZ
    )); // vert 3

    cubeVerts.push(GsVec(
        ((float)x / numCellsX) * (maxX - minX) + minX,
        ((float)(y + 1) / numCellsY) * (maxY - minY) + minY,
        ((float)z / numCellsZ) * (maxZ - minZ) + minZ
    )); // vert 4

    cubeVerts.push(GsVec(
        ((float)(x + 1) / numCellsX) * (maxX - minX) + minX,
        ((float)(y + 1) / numCellsY) * (maxY - minY) + minY,
        ((float)z / numCellsZ) * (maxZ - minZ) + minZ
    )); // vert 5

    cubeVerts.push(GsVec(
        ((float)(x + 1) / numCellsX) * (maxX - minX) + minX,
        ((float)(y + 1) / numCellsY) * (maxY - minY) + minY,
        ((float)(z + 1) / numCellsZ) * (maxZ - minZ) + minZ
    )); // vert 6

    cubeVerts.push(GsVec(
        ((float)x / numCellsX) * (maxX - minX) + minX,
        ((float)(y + 1) / numCellsY) * (maxY - minY) + minY,
        ((float)(z + 1) / numCellsZ) * (maxZ - minZ) + minZ
    )); // vert 7

    // Get the verts that compose the surface for this particular cube
    GsArray<GsVec> surfaceVerts = get_verts_from_cube(cubeIndex, cubeVerts, cubeMetaballValues); // surfaceVerts always has a size of 12

    if (surfaceVerts.size() == 0) // this will be the case when cube is completely within or completely outside of surface
        return;

    // GsArray<GsVec> surfaceNormals = get_normals_from_verts(surfaceVerts); // surfaceNormals always has the same size as cubeVerts
    CubeSurface surface(cubeIndex, surfaceVerts);
    surfaceQueue.push(surface);
}

void MyViewer::clear_queue()
{
    std::queue<CubeSurface>().swap(surfaceQueue);
}

void MyViewer::marching_cubes()
{
    if (m_metaballSurface)
        rootg()->remove(m_metaballSurface); // remove old surface from scene (if it exists)

    clear_queue();

    // Initialize surface
    rootg()->add(m_metaballSurface = new SnModel);
    m_metaballSurface->color(GsColor(1.0f, 0.0f, 1.0f, 0.8f));

    // Loop through each index in space, and calculate the marching cube at that location
    for (int i = 0; i < m_numCells[0] - 1; i++)
    {
        for (int j = 0; j < m_numCells[1] - 1; j++)
        {
            for (int k = 0; k < m_numCells[2] - 1; k++)
            {
                // Push the cube metaball values at each vertex to be used later when performing linear interpolation
                GsArray<float> cubeVertMetaVals;
                cubeVertMetaVals.push(m_metaballValues[i][j][k]);
                cubeVertMetaVals.push(m_metaballValues[i + 1][j][k]);
                cubeVertMetaVals.push(m_metaballValues[i + 1][j][k + 1]);
                cubeVertMetaVals.push(m_metaballValues[i][j][k + 1]);
                cubeVertMetaVals.push(m_metaballValues[i][j + 1][k]);
                cubeVertMetaVals.push(m_metaballValues[i + 1][j + 1][k]);
                cubeVertMetaVals.push(m_metaballValues[i + 1][j + 1][k + 1]);
                cubeVertMetaVals.push(m_metaballValues[i][j + 1][k + 1]);

                // Open a thread to process this specific cube
                cube_thread(
                    i, j, k,
                    m_worldBounds[0][0], m_worldBounds[0][1], m_worldBounds[1][0], m_worldBounds[1][1], m_worldBounds[2][0], m_worldBounds[2][1],
                    m_numCells[0], m_numCells[1], m_numCells[2],
                    cubeVertMetaVals);
            }
        }
    }

    int cubeCount = 0;
    while (!surfaceQueue.empty())
    {
        CubeSurface tempCube = surfaceQueue.front();
        surfaceQueue.pop();

        int cubeIndex = tempCube.getCubeIndex();
        GsArray<GsVec> surfaceVerts = tempCube.getVerts();

        // I load verts and faces in an inefficient way.
        // Since, m_triTable stores duplicates of indices, I just load the same vertex into the model's vertex array multiple times
        //   then, the faces are just 0, 1, 2, 3, ... , V.size()-1 since I load some vertices multiple times (and only need to reference each vertex once)
        for (int i = 0; m_triTable[cubeIndex][i] != -1; i += 3)
        {
            // m_triTable contains the INDEX of the edge where the surface vertex can be found
            // Use this index to get the actual world coordinates of the surface point along that edge;
            // After getting the world coordinates of that point, calculate the normal of the vertex
            m_metaballSurface->model()->V.push(surfaceVerts[m_triTable[cubeIndex][i]]);
            m_metaballSurface->model()->N.push(get_normal_from_vec(surfaceVerts[m_triTable[cubeIndex][i]]));

            m_metaballSurface->model()->V.push(surfaceVerts[m_triTable[cubeIndex][i + 1]]);
            m_metaballSurface->model()->N.push(get_normal_from_vec(surfaceVerts[m_triTable[cubeIndex][i + 1]]));

            m_metaballSurface->model()->V.push(surfaceVerts[m_triTable[cubeIndex][i + 2]]);
            m_metaballSurface->model()->N.push(get_normal_from_vec(surfaceVerts[m_triTable[cubeIndex][i + 2]]));

            m_metaballSurface->model()->F.push(GsModel::Face(
                m_metaballSurface->model()->V.size() - 3,
                m_metaballSurface->model()->V.size() - 2,
                m_metaballSurface->model()->V.size() - 1
            ));
        }
    } 

    m_metaballSurface->model()->set_mode(GsModel::Smooth, GsModel::NoMtl); // set the model to be smooth shading
}

float MyViewer::vec_distance(GsVec a, GsVec b)
{
    return sqrtf(gs_pow(a.x - b.x, 2) + gs_pow(a.y - b.y, 2) + gs_pow(a.z - b.z, 2));
}

void MyViewer::animate()
{
    camera().eye = GsVec(10, 10, 10);
    camera().center = GsVec(0, 0, 0);

    int maxX = m_worldBounds[0][1] * 2;
    int minX = m_worldBounds[0][0];
    int maxY = m_worldBounds[1][1] * 2;
    int minY = m_worldBounds[1][0];
    int maxZ = m_worldBounds[2][1] * 2;
    int minZ = m_worldBounds[2][0];

    float fr = 30; // framerate
    float spf = 1 / fr; // seconds per frame
    std::chrono::duration<double> dt; // the timespan of the last frame
    while (m_animating)
    {
        // Do nothing while we wait for the next frame to open
        auto frameStart = std::chrono::system_clock::now();
        do
        {
            ws_check();
            auto frameCurrent = std::chrono::system_clock::now();
            dt = frameCurrent - frameStart;
        } while (dt.count() < spf);

        for (int i = 0; i < m_numCenters; i++)
        {
            // Get the current position of the metaball center
            float metaballX = m_metaballCenters[i]->mat().e14;
            float metaballY = m_metaballCenters[i]->mat().e24;
            float metaballZ = m_metaballCenters[i]->mat().e34;
            GsVec center(metaballX, metaballY, metaballZ);

            // If the metaball has reached the endpoint, choose another random endpoint in the space
            if (vec_distance(center, m_endpoints[i]) < 0.75f)
            {
                GsVec newEndpoint(rand() % maxX + minX, rand() % maxY + minY, rand() % maxZ + minZ);
                m_endpoints[i] = newEndpoint;
            }

            // Linearly interpolate the position of the metaball center
            float t = dt.count() * m_animationSpeed;
            center = center + (m_endpoints[i] - center) * t;

            m_metaballCenters[i]->mat().e14 = center.x;
            m_metaballCenters[i]->mat().e24 = center.y;
            m_metaballCenters[i]->mat().e34 = center.z;
        }


        calculate_metaballs();
        marching_cubes();
        render();
    }
}

int MyViewer::handle_keyboard(const GsEvent &e)
{
    int ret = WsViewer::handle_keyboard(e); // 1st let system check events
    if (ret) return ret;

    switch (e.key)
    {
    case GsEvent::KeyEsc:
        gs_exit();
        break;
    case GsEvent::KeySpace:
        calculate_metaballs();
        marching_cubes();
        render();
        break;
    case 97:
        m_animating = !m_animating;
        animate();
        break;
    default:
        gsout << "Key pressed: " << e.key << gsnl;
    }

    return 0;
}

int MyViewer::uievent(int e)
{
    switch (e)
    {
    case EvFidelity:
        break;
    case EvExit:
        gs_exit();
        break;
    }
    return WsViewer::uievent(e);
}

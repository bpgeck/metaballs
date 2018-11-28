# Metaballs

![Professional metaballs visualized](https://i.imgur.com/ZVUssmW.jpg)

## What are metaballs?
Metaballs are smooth, blob-like surfaces constructed in 3D space that are used in computer graphics.

## Creation process
I begin by randomly choosing points in 3D space to serve as my metaball "centers".

![Metaball centers](https://i.imgur.com/HwMGuat.png)

I then interpret all points in my universe as having a value of 0. Then I increase the value at that point by its distance to a metaball center: f(x,y,z) = 1.0 / ((x - x_m)^2 + (y - y_m)^2 + (z - z_m)^2) where (x, y, z) is the point, and (x_m, y_m, z_m) is a metaball center. This process is repeated at each point for each metaball center (if there are 4 metaball centers in the universe, the value at the point is increased 4 times using the above equation. The smaller the distance to a metaball center, the bigger the increase).

After all points have been evaluated, the points that have a final value larger than a specified threshold are marked as being "inside" a metaball. Which metaball is not known, only that the point is close enough to a metaball that it is considered to be "inside".

![Points inside metaballs](https://i.imgur.com/cS0j1WO.png)

The surfaces of the metaballs are then constructed using the marching cubes algorithm. Smoothing is performed by finding the normal at each vertex and interpolating normals across the reamining surface.

![Surfaces](https://i.imgur.com/CJWhlu4.png)

And just for fun, I throw in some multithreading so constructing surfaces isn't awful.

Thanks for reading!

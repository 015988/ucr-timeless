
Accurate Tridimensional Reconstruction with Unsynchronized Cameras Regardless of Time Information

In this work, we approach the problem of the tridimensional reconstruction of the trajectory of an object using unsynchronized cameras. While most methods from the literature try to measure, as accurately as possible, the timing difference between the cameras, we chose a very different path: we ignore most, or all, of the time information from the videos. The idea is quite simple. Each pair of camera and projected trajectory generate a surface on the volume. The intersection of these curves is the desired result. Since this intersection can be very complex, we considered a Monte Carlo approach for the reconstruction, using random tridimensional points to estimate the region of intersection. Therefore, any camera calibration schema can be used. These points are later used to compute a continuous curve, which is the final result of the method. We compared this method to a very simple reconstruction approach, which assumes the frames are synchronized, and obtained outstanding results.

Original article: http://www.ic.unicamp.br/~reltech/2014/14-04.pdf

The VolumetricReconstruction subfolder contains the extension of this method to a full 3d reconstruction. This step needs a model to fit into the reconstructed points. We included a simple model (ball) as example.

Also available is a variant of space carving, inspired by the quadtree concept, zipfile.

FAPESP 2012/15811-1

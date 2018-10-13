## A ParaView Plugin for Disco ##

This is a C++ reader plugin for ParaView. It reads Disco checkpoints (HDF5 files with a .h5 extension) and loads them into VTK a VTK object suitable for visualization with Paraview.  Currently the mesh is formed with Hexahedron cells in a `vtkUnstructuredGrid'.  

#### Install ####

You first need a working installation of ParaView (preferably built from source, as binary distribution does not include the required headers).  Set the `Paraview_DIR' environment variable to the root directory of your ParaView install.

In a separate directory (e.g. `build/' within this project), call `cmake' with the root directory of `DiscoPV' as an argument.  If `build/' is in this project's root directory, you can just call `cmake ..'.  Then call `make' from the same directory you called `cmake'.  This will build the `libDiscoReader.so' shared library.

Make sure the directory containing `libDiscoReader.so' is in the `PV_PLUGIN_PATH' environment variable. It should be automatically loaded when you launch ParaView!


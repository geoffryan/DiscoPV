## A ParaView Plugin for Disco ##

This is a C++ reader plugin for ParaView. It reads Disco checkpoints (HDF5 files with a .h5 extension) and loads them into VTK a VTK object suitable for visualization with Paraview.  Currently the mesh is formed with Hexahedron cells in a `vtkUnstructuredGrid'.  

#### Install ####

You first need a working installation of ParaView (preferably built from source, as binary distribution does not include the required headers).  Set the `Paraview_DIR' environment variable to the root directory of your ParaView install.

In a separate directory (e.g. `build/' within this project), call `cmake' with the root directory of `DiscoPV' as an argument.  If `build/' is in this project's root directory, you can just call `cmake ..'.  Then call `make' from the same directory you called `cmake'.  This will build the `libDiscoReader.so' shared library.

Make sure the directory containing `libDiscoReader.so' is in the `PV_PLUGIN_PATH' environment variable. It should be automatically loaded when you launch ParaView!

#### Current Status ####

Currently we map the Disco grid onto a vtkUnstructuredMesh formed (manually) of tetrahedra (VTK_TETRA).  This was the simplest way (and still not that simple...) to guarantee a water tight mesh for the curvilinear grid.  This means *each* Disco cell is decomposed into ~24 tetrahedra, and each grid function value is copied into each of these cells.  This creates a not-so-small memory footprint and a bit of a performance hurdle (ie. can't animate 3D on a laptop).  More efficient meshes (like ~5 hexahedra per Disco cell) have noticably faster performance but are not watertight, limiting the visualization options (no streamlines).  The Polyhedra mesh is also watertight, faster to load,  and seems quicker at Volume renderings but is significantly slower at streamlines than the Tetra mesh. Open to suggestions to make a more efficient mesh!


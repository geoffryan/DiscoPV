#ifndef vtkDiscoReader_h
#define vtkDiscoReader_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkUnstructuredGrid.h"

#define H5_USE_16_API
#include <vtk_hdf5.h>

enum Geom{CARTESIAN, CYLINDRICAL, SPHERICAL};

class vtkDiscoReader : public vtkUnstructuredGridAlgorithm
{
public:
    static vtkDiscoReader *New();
    vtkTypeMacro(vtkDiscoReader, vtkUnstructuredGridAlgorithm);

    int CanReadFile(const char *fname);
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);
    void SetMeshType(int type);
    int GetMeshType();
protected:
    vtkDiscoReader();
    ~vtkDiscoReader();
    int RequestInformation(vtkInformation*, vtkInformationVector**,
                            vtkInformationVector*) override;
    int RequestData(vtkInformation*, vtkInformationVector**,
                            vtkInformationVector*) override;

    char *FileName;
private:
    Geom geometry;
    int Nz;
    int Nr;
    int NrTot;
    int NzTot;
    int jA;  //Indices into global array setting our subextent
    int jB;
    int kA;
    int kB;
    int cja; //Indices into local arrays ignoring extraneous outer cells.
    int cjb;
    int cka;
    int ckb;
    int ghostLevel;
    std::vector<int> Np;
    std::vector<int> Index;
    std::vector<int> Id_phi0;
    std::vector<double> rjmh;
    std::vector<double> zkmh;
    std::vector<std::vector<double>> pimh;
    double phimax;
    std::vector<int> chunkSize;
    std::vector<std::vector<int>> ncellspercell;
    int ncells;
    int cellType;
    void SetGeometry();
    int LoadGridDims(int piece, int numPieces, int ghostLevel);
    void MakeGrid(vtkUnstructuredGrid *output);
    void LoadGrid();
    void AddCells(vtkUnstructuredGrid *output,
                        const std::vector<std::vector<vtkIdType>> &cellPoints);
    void AddCellsTetrahedral(vtkUnstructuredGrid *output,
                        const std::vector<std::vector<vtkIdType>> &cellPoints);
    void AddCellsPolyhedral(vtkUnstructuredGrid *output,
                        const std::vector<std::vector<vtkIdType>> &cellPoints);
    void TagGhostCells(vtkUnstructuredGrid *output);
    void AddData(vtkUnstructuredGrid *output);
    void AddDataScalar(vtkUnstructuredGrid *output, const char *name, int id);
    void AddDataVector(vtkUnstructuredGrid *output, const char *name, int id0,
                        int basis);
    int AddDiscoCell(vtkUnstructuredGrid *output, 
                                    vtkIdType sDLB, vtkIdType sDRB, 
                                    vtkIdType sULB, vtkIdType sURB, 
                                    vtkIdType sDLF, vtkIdType sDRF, 
                                    vtkIdType sULF, vtkIdType sURF);
    int AddFacesToStream(std::vector<vtkIdType> &faceStream,
                                        vtkIdType sBL, vtkIdType sBR,
                                        vtkIdType sFL, vtkIdType sFR);
    int which4(double, double, double, double);
    void order4(double x0, double x1, double x2, double x3,
                            int *i0, int *i1, int *i2, int *i3);
    void readPatch(const char *group, const char *dset, void *data, 
                                hid_t type, int dim, hsize_t *start, 
                                hsize_t *loc_size, hsize_t *glo_size);
    void readSimple(const char *group, const char *dset, void *data, 
                        hid_t type);
    void readString(const char *group, const char *dset, char *buf, 
                        int len);
    void getH5dims(const char *group, const char *dset, hsize_t *dims);
    void calcXYZ(double x1, double x2, double x3, double *xyz);
    void calcVFromContravariant(double x1, double x2, double x3, double *V, 
                                double *Vxyz);
    void calcVFromOrthonormal(double x1, double x2, double x3, double *V, 
                                double *Vxyz);
    void dimsCreate(int numPieces, int dims[]);
};

#endif

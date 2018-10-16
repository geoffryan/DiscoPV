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
    std::vector<int> Np;
    std::vector<int> Index;
    std::vector<int> Id_phi0;
    std::vector<double> rjmh;
    std::vector<double> zkmh;
    std::vector<std::vector<double>> pimh;
    double phimax;
    std::vector<std::vector<int>> ncellspercell;
    int ncells;
    int cellType;
    void SetGeometry();
    void MakeGrid(vtkUnstructuredGrid *output);
    void LoadGrid();
    void AddCells(vtkUnstructuredGrid *output,
                        const std::vector<std::vector<vtkIdType>> &cellPoints);
    void AddCellsTetrahedral(vtkUnstructuredGrid *output,
                        const std::vector<std::vector<vtkIdType>> &cellPoints);
    void AddCellsPolyhedral(vtkUnstructuredGrid *output,
                        const std::vector<std::vector<vtkIdType>> &cellPoints);
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
    void getH5dims(const char *group, const char *dset, hsize_t *dims);
    void calcXYZ(double x1, double x2, double x3, double *xyz);
    void calcVFromContravariant(double x1, double x2, double x3, double *V, 
                                double *Vxyz);
    void calcVFromOrthonormal(double x1, double x2, double x3, double *V, 
                                double *Vxyz);
};

#endif

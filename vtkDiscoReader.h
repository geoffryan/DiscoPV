#ifndef vtkDiscoReader_h
#define vtkDiscoReader_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkUnstructuredGrid.h"

#define H5_USE_16_API
#include <vtk_hdf5.h>

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
    void MakeGridHexahedron(vtkUnstructuredGrid *output);
    int which4(double, double, double, double);
    void readPatch(const char *group, const char *dset, void *data, 
                                hid_t type, int dim, int *start, 
                                int *loc_size, int *glo_size);
    void readSimple(const char *group, const char *dset, void *data, 
                        hid_t type);
    void getH5dims(const char *group, const char *dset, hsize_t *dims);
};

#endif

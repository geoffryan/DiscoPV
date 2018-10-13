#ifndef vtkDiscoReader_h
#define vtkDiscoReader_h

#include "vtkUnstructuredGridAlgorithm.h"

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
    int which4(double, double, double, double);
};

#endif

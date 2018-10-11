#ifndef vtkDiscoReader_h
#define vtkDiscoReader_h

#include "vtkStructuredGridAlgorithm.h"

class vtkDiscoReader : public vtkStructuredGridAlgorithm
{
public:
    static vtkDiscoReader *New();
    vtkTypeMacro(vtkDiscoReader, vtkStructuredGridAlgorithm);

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
};

#endif

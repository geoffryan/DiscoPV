#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPoints.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkDiscoReader.h"

int vtkDiscoReader::RequestInformation(vtkInformation*, vtkInformationVector**,
                                        vtkInformationVector* outVec)
{
    vtkInformation *outInfo = outVec->GetInformationObject(0);
    int extent[6];

    //read file for extent
    extent[0] = 0;
    extent[1] = 6;
    extent[2] = 0;
    extent[3] = 2;
    extent[4] = 0;
    extent[5] = 1;
    
    //store extent
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

    //store other things

    return 1;
}

int vtkDiscoReader::RequestData(vtkInformation*, vtkInformationVector**,
                                        vtkInformationVector* outVec)
{
    vtkInformation *outInfo = outVec->GetInformationObject(0);

    vtkStructuredGrid *outData = vtkStructuredGrid::SafeDownCast
                                (outInfo->Get(vtkDataObject::DATA_OBJECT()));

    //int extent[6];
    //outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);
    
    int dims[3] = {7, 3, 2};

    outData->SetDimensions(dims);

    vtkPoints* points = vtkPoints::New();
    points->Allocate(dims[0]*dims[1]*dims[2]);

    double zmin = -0.5;
    double zmax = 0.5;
    double rmin = 0.1;
    double rmax = 1.0;
    double phimin = 0.0;
    double phimax = M_PI;

    int i,j,k;
    for(k=0; k<dims[2]; k++)
        for(j=0; j<dims[1]; j++)
            for(i=0; i<dims[0]; i++)
            {
                double z = zmin + k*(zmax-zmin)/(dims[2]-1);
                double r = rmin + j*(rmax-rmin)/(dims[1]-1);
                double phi = phimin + i*(phimax-phimin)/(dims[0]-1);
                double x[3] = {phi, r, z};
                int s = dims[1]*dims[0]*k + dims[0]*j + i;
                points->InsertPoint(s, x);
            }

    vtkDoubleArray *rho = vtkDoubleArray::New();
    rho->SetName("rho");
    rho->SetNumberOfValues((dims[0]-1)*(dims[1]-1)*(dims[2]-1));

    for(k=0; k<dims[2]-1; k++)
        for(j=0; j<dims[1]-1; j++)
            for(i=0; i<dims[0]-1; i++)
            {
                int s = dims[1]*dims[0]*k + dims[0]*j + i;
                double x = (double)s;
                rho->SetValue(s, x);
            }
    outData->SetPoints(points);
    outData->GetCellData()->SetScalars(rho);
    points->Delete();
    rho->Delete();
    
    return 1;
}



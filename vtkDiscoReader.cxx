#include <iostream>
#include <fstream>
#include <cmath>
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPoints.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredGrid.h"
#include "vtkDiscoReader.h"

vtkStandardNewMacro(vtkDiscoReader);

vtkDiscoReader::vtkDiscoReader()
{
    ofstream myfile;
    myfile.open("/Users/geoff/Projects/DiscoPV/construct.txt");
    myfile << "line";
    myfile.close();

    this->FileName = NULL;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
}

vtkDiscoReader::~vtkDiscoReader()
{
    ofstream myfile;
    myfile.open("/Users/geoff/Projects/DiscoPV/destruct.txt");
    myfile << "line";
    myfile.close();

    this->SetFileName(0);
}

int vtkDiscoReader::CanReadFile(const char *fname)
{
    return 1;
}

int vtkDiscoReader::RequestInformation(vtkInformation*, vtkInformationVector**,
                                        vtkInformationVector* outVec)
{
    ofstream myfile;
    myfile.open("/Users/geoff/Projects/DiscoPV/info.txt");
    myfile << "line";
    myfile.close();

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
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
            0,6,0,2,0,1);
            //extent[0], extent[1], extent[2], extent[3], extent[4], extent[5]);

    //store other things

    return 1;
}

int vtkDiscoReader::RequestData(vtkInformation*, vtkInformationVector**,
                                        vtkInformationVector* outVec)
{
    ofstream myfile;
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt");
    myfile << "line";
    myfile.close();

    vtkInformation *outInfo = outVec->GetInformationObject(0);

    vtkStructuredGrid *outData = vtkStructuredGrid::SafeDownCast
                                (outInfo->Get(vtkDataObject::DATA_OBJECT()));

    //int extent[6];
    //outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);
    
    int dims[3] = {7, 3, 2};

    outData->SetDimensions(dims);
    
    outData->GetDimensions(dims);
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", std::ios_base::app);
    myfile << "\n" << dims[0] << dims[1] << dims[2];
    myfile.close();

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
                double x[3] = {r*cos(phi), r*sin(phi), z};
                int s = dims[1]*dims[0]*k + dims[0]*j + i;
                points->InsertPoint(s, x);
            }
    outData->SetPoints(points);
    points->Delete();


    vtkDoubleArray *rho = vtkDoubleArray::New();
    rho->SetName("rho");
    rho->SetNumberOfValues((dims[0]-1)*(dims[1]-1)*(dims[2]-1));

    for(k=0; k<dims[2]-1; k++)
        for(j=0; j<dims[1]-1; j++)
            for(i=0; i<dims[0]-1; i++)
            {
                int s = (dims[1]-1)*(dims[0]-1)*k + (dims[0]-1)*j + i;
                double x = (double)s;
                rho->SetValue(s, x);
            }
    outData->GetCellData()->SetScalars(rho);
    rho->Delete();

    myfile.open("/Users/geoff/Projects/DiscoPV/data2.txt");
    myfile << "line";
    myfile.close();
    
    return 1;
}



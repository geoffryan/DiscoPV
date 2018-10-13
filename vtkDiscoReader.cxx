#include <iostream>
#include <fstream>
#include <cmath>
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPoints.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#define H5_USE_16_API
#include <vtk_hdf5.h>

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

    return 1;
}

int vtkDiscoReader::RequestData(vtkInformation*, vtkInformationVector**,
                                        vtkInformationVector* outVec)
{
    ofstream myfile;
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt");
    myfile << "line\n";
    myfile.close();

    vtkInformation *outInfo = outVec->GetInformationObject(0);

    vtkUnstructuredGrid *outData = vtkUnstructuredGrid::SafeDownCast
                                (outInfo->Get(vtkDataObject::DATA_OBJECT()));

    int Nr = 32;
    int Nz = 32;
    int Np = 64;

    double zmin = -0.5;
    double zmax = 0.5;
    double rmin = 0.1;
    double rmax = 1.0;
    double phimin = 0.0;
    double phimax = 1.5*M_PI;

    double rjmh[Nr+1];
    double zkmh[Nz+1];
    double pimh[Nr*Nz][Np+1];

    int i,j,k;

    for(j=0; j<=Nr; j++)
        rjmh[j] = rmin + j*(rmax-rmin)/Nr;

    for(k=0; k<=Nz; k++)
        zkmh[k] = zmin + k*(zmax-zmin)/Nz;

    srand(42);
    rand();
    rand();
    
    double dphi = (phimax-phimin) / Np;

    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = Nr*k + j;
            double p0 = -rand()*dphi/RAND_MAX;
            for(i=0; i<=Np; i++)
                pimh[jk][i] = p0 + i*dphi;
        }
    
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
    myfile << "Made grid\n";
    myfile.close();


    vtkPoints* points = vtkPoints::New();

    vtkIdType npointsTot = 0;
    vtkIdType npoints[Nz+1][Nr+1];
    vtkIdType indexPoints[Nz+1][Nr+1];
    vtkIdType cellPoints[Nz*Nr][4*(Np+1)];

    for(k=0; k<=Nz; k++)
        for(j=0; j<=Nr; j++)
        {
            int na = 4;
            if((j==0 && k==0) || (j==0 && k==Nz)
                    || (j==Nr && k==0) || (j==Nr && k==Nz))
                na = 1;
            else if(j==0 || j==Nr || k==0 || k==Nz)
                na = 2;

            indexPoints[k][j] = npointsTot;
            npoints[k][j] = na*(Np+1);
            npointsTot += npoints[k][j];
        }
    points->SetNumberOfPoints(npointsTot);
    
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
    myfile << "Counted Points\n";
    myfile.close();

    for(k=0; k<=Nz; k++)
    {
        double z = zkmh[k];

        int kD = k-1;
        int kU = k;

        for(j=0; j<=Nr; j++)
        {
            int iDL = 0;
            int iDR = 0;
            int iUL = 0;
            int iUR = 0;

            double r = rjmh[j];

            int jL = j-1;
            int jR = j;

            int jkDL = Nr*kD + jL;
            int jkDR = Nr*kD + jR;
            int jkUL = Nr*kU + jL;
            int jkUR = Nr*kU + jR;

            double phiDL, phiDR, phiUL, phiUR;

            if(kD >= 0 && jL >= 0)
                phiDL = pimh[jkDL][0];
            else
                phiDL = HUGE;
            if(kD >= 0 && jR < Nr)
                phiDR = pimh[jkDR][0];
            else
                phiDR = HUGE;
            if(kU < Nz && jL >= 0)
                phiUL = pimh[jkUL][0];
            else
                phiUL = HUGE;
            if(kU < Nz && jR < Nr)
                phiUR = pimh[jkUR][0];
            else
                phiUR = HUGE;
         
            /*
            myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                            std::ios_base::app);
            myfile << "k " << k << " j " << j << "\n";
            myfile.close();
            */

            for(i=0; i<npoints[k][j]; i++)
            {
                vtkIdType s = indexPoints[k][j] + i;

                int jkmin, *imin, opp;
                double *phimin;

                int which4 = this->which4(phiDL, phiDR, phiUL, phiUR);
                if(which4 == 0)
                {
                    imin = &iDL;
                    jkmin = jkDL;
                    phimin = &phiDL;
                    opp = 3;
                }
                else if(which4 == 1)
                {
                    imin = &iDR;
                    jkmin = jkDR;
                    phimin = &phiDR;
                    opp = 2;
                }
                else if(which4 == 2)
                {
                    imin = &iUL;
                    jkmin = jkUL;
                    phimin = &phiUL;
                    opp = 1;
                }
                else
                {
                    imin = &iUR;
                    jkmin = jkUR;
                    phimin = &phiUR;
                    opp = 0;
                }
                /*

                myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                                std::ios_base::app);
                myfile << s << " " << *phimin << " ";
                myfile.close();
                */
                double x[3] = {r*cos(*phimin), r*sin(*phimin), z};

                points->InsertPoint(s, x);
                cellPoints[jkmin][4*(*imin)+opp] = s;

                (*imin)++;
                if(*imin <= Np)
                    *phimin = pimh[jkmin][*imin];
                else
                    *phimin = HUGE;
            }
            /*
            myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                            std::ios_base::app);
            myfile <<  "\n";
            myfile.close();
            */
        }
    }

    outData->SetPoints(points);
    points->Delete();
    
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
    myfile << "Made Points\n";
    myfile.close();

    int ncellsTot = 0;
    int ncells[Nz*Nr][Np];

    vtkCellArray *cells = vtkCellArray::New();


    for(k=0; k<Nz; k++)
    {
        for(j=0; j<Nr; j++)
        {
            int jk = Nr*k + j;

            myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                std::ios_base::app);
            myfile << "jk " << jk <<"\n";
            myfile.close();

            for(i=0; i<Np; i++)
            {
                vtkIdType sDLB = cellPoints[jk][4*i+0];
                vtkIdType sDRB = cellPoints[jk][4*i+1];
                vtkIdType sULB = cellPoints[jk][4*i+2];
                vtkIdType sURB = cellPoints[jk][4*i+3];
                vtkIdType sDLF = cellPoints[jk][4*(i+1)+0];
                vtkIdType sDRF = cellPoints[jk][4*(i+1)+1];
                vtkIdType sULF = cellPoints[jk][4*(i+1)+2];
                vtkIdType sURF = cellPoints[jk][4*(i+1)+3];

                int nDL = (int)(sDLF-sDLB);
                int nDR = (int)(sDRF-sDRB);
                int nUL = (int)(sULF-sULB);
                int nUR = (int)(sURF-sURB);

                int nc = std::max(std::max(nDL,nDR), std::max(nUL,nUR));

                ncells[jk][i] = nc;

                vtkIdType p[8] = {sDRB, sDLB, sULB, sURB, 
                            sDRB+1, sDLB+1, sULB+1, sURB+1};
                cells->InsertNextCell(8, p);

                /*
                myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
                myfile << "p " << p[0] << " " << p[1] << " " << p[2]
                            << " " << p[3] << " " << p[4] << " " << p[5]
                            << " " << p[6] << " " << p[7] <<"\n";
                myfile << "pF " << sDRF << " " << sDLF << " " << sULF << " " 
                        << sURF << "\n";
                myfile.close();
                */

                while(p[4] < sDRF || p[5] < sDLF || p[6] < sULF || p[7]<sURF)
                {
                    p[0] = p[4];
                    p[1] = p[5];
                    p[2] = p[6];
                    p[3] = p[7];
                    if(p[4] < sDRF)
                        p[4]++;
                    if(p[5] < sDLF)
                        p[5]++;
                    if(p[6] < sULF)
                        p[6]++;
                    if(p[7] < sURF)
                        p[7]++;
                    cells->InsertNextCell(8, p);

                    /*
                    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                        std::ios_base::app);
                    myfile << p[4] << " " << p[5] << " " << p[6] << " " 
                            << p[7] <<"\n";
                    myfile.close();
                    */
                }
            }
        }
    }

    int nc = cells->GetNumberOfCells();

    outData->SetCells(VTK_HEXAHEDRON, cells);
    cells->Delete();
    
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
    myfile << "Made Cells\n";
    myfile.close();

    vtkDoubleArray *arho = vtkDoubleArray::New();
    arho->SetName("rho");
    arho->SetNumberOfValues(nc);
    vtkDoubleArray *aP = vtkDoubleArray::New();
    aP->SetName("P");
    aP->SetNumberOfValues(nc);
    vtkDoubleArray *aom = vtkDoubleArray::New();
    aom->SetName("om");
    aom->SetNumberOfValues(nc);
    
    vtkDoubleArray *av = vtkDoubleArray::New();
    av->SetName("v");
    av->SetNumberOfComponents(3);
    av->SetNumberOfTuples(nc);

    int s = 0;
    for(k=0; k<Nz; k++)
    {
        double z = 0.5*(zkmh[k] + zkmh[k+1]);
        for(j=0; j<Nr; j++)
        {
            double r = 0.5*(rjmh[j]+rjmh[j+1]);

            int jk = Nr*k + j;

            for(i=0; i<Np; i++)
            {
                double phi = 0.5*(pimh[jk][i]+pimh[jk][i+1]);
                
                double vr = 0.1*(rand()*2.0/RAND_MAX-1);
                double vp = sqrt(1./r);
                double V[3] = {vr*cos(phi) - vp*sin(phi), vr*sin(phi)+vp*cos(phi), 0.0};

                double R = 0.2 + phi/(2*M_PI);
                double Z = -0.5 + phi/(2*M_PI);
                double sig = 0.1;
                double phase = -((r-R)*(r-R) +(z-Z)*(z-Z)) / (2*sig*sig);
                double rho = 1.0 + 1.0*exp(phase);
                double P = phi;
                double om = pow(r, -1.5);

                int c;
                for(c=0; c<ncells[jk][i]; c++)
                {
                    arho->SetValue(s, rho);
                    aP->SetValue(s, P);
                    aom->SetValue(s, om);
                    av->SetTuple(s, V);
                    s++;
                }


            }
        }
    }

    outData->GetCellData()->AddArray(arho);
    outData->GetCellData()->AddArray(aP);
    outData->GetCellData()->AddArray(aom);
    outData->GetCellData()->AddArray(av);

    arho->Delete();
    aP->Delete();
    aom->Delete();
    av->Delete();


    /*
            myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                            std::ios_base::app);
            myfile << points->GetNumberOfPoints() << "\n";
            myfile.close();


            for(i=0; i<dims[0]; i++)
            {
                double phi = p0 + i*dphi;

                int s = dims[0]*(dims[1]*0+0) + i;
                double x[3] = {rm*cos(phi), rm*sin(phi), zm};
                points->InsertPoint(s, x);

                s = dims[0]*(dims[1]*0+1) + i;
                x[0] = rp*cos(phi);
                x[1] = rp*sin(phi);
                x[2] = zm;
                points->InsertPoint(s, x);

                s = dims[0]*(dims[1]*1+0) + i;
                x[0] = rm*cos(phi);
                x[1] = rm*sin(phi);
                x[2] = zp;
                points->InsertPoint(s, x);

                s = dims[0]*(dims[1]*1+1) + i;
                x[0] = rp*cos(phi);
                x[1] = rp*sin(phi);
                x[2] = zp;
                points->InsertPoint(s, x);
            }

            myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                            std::ios_base::app);
            myfile << points->GetNumberOfPoints() << "\n";
            myfile.close();

            annulus->SetPoints(points);
            points->Delete();

            vtkDoubleArray *rho = vtkDoubleArray::New();
            rho->SetName("rho");
            rho->SetNumberOfValues((dims[0]-1)*(dims[1]-1)*(dims[2]-1));
            vtkDoubleArray *P = vtkDoubleArray::New();
            P->SetName("P");
            P->SetNumberOfValues((dims[0]-1)*(dims[1]-1)*(dims[2]-1));
            vtkDoubleArray *om = vtkDoubleArray::New();
            om->SetName("om");
            om->SetNumberOfValues((dims[0]-1)*(dims[1]-1)*(dims[2]-1));
            
            vtkDoubleArray *v = vtkDoubleArray::New();
            v->SetName("v");
            v->SetNumberOfComponents(3);
            v->SetNumberOfTuples((dims[0]-1)*(dims[1]-1)*(dims[2]-1));

            for(i=0; i<dims[0]-1; i++)
            {
                int s = i;
                double phi = p0+(i+0.5)*dphi;

                double vr = 0.1*(rand()*2.0/RAND_MAX-1);
                double vp = sqrt(1./r);
                double V[3] = {vr*cos(phi) - vp*sin(phi), vr*sin(phi)+vp*cos(phi), 0.0};

                rho->SetValue(s, 1.0 + 0.1*(rand()*2.0/RAND_MAX-1));
                P->SetValue(s, phi);
                om->SetValue(s, pow(r, -1.5));
                v->SetTuple(s, V);
            }
            annulus->GetCellData()->AddArray(rho);
            annulus->GetCellData()->AddArray(P);
            annulus->GetCellData()->AddArray(om);
            annulus->GetCellData()->AddArray(v);
            rho->Delete();
            P->Delete();
            om->Delete();
            v->Delete();

            myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                            std::ios_base::app);
            myfile << annulus->GetNumberOfPoints() << "\n";
            myfile.close();


            outData->SetBlock(jk, annulus);
            annulus->Delete();
        }
    */


    /* 
    int dims[3] = {7, 3, 2};

    outData->SetDimensions(dims);
    
    outData->GetDimensions(dims);
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", std::ios_base::app);
    myfile << "\n" << dims[0] << dims[1] << dims[2];
    myfile.close();

    vtkPoints* points = vtkPoints::New();
    points->Allocate(dims[0]*dims[1]*dims[2]);


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
    */

    myfile.open("/Users/geoff/Projects/DiscoPV/data2.txt");
    myfile << "line";
    myfile.close();
    
    return 1;
}

int vtkDiscoReader::which4(double x0, double x1, double x2, double x3)
{
    int i01, i23;
    double x01, x23;

    if(x0 < x1)
    {
        i01 = 0;
        x01 = x0;
    }
    else
    {
        i01 = 1;
        x01 = x1;
    }

    if(x2 < x3)
    {
        i23 = 2;
        x23 = x3;
    }
    else
    {
        i23 = 3;
        x23 = x3;
    }

    if(x01 < x23)
        return i01;
    else
        return i23;
}



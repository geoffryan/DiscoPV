#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPoints.h"
#include "vtkStreamingDemandDrivenPipeline.h"

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

    this->MakeGridHexahedron(outData);

    int Nz = this->Nz;
    int Nr = this->Nr;
    int nc = this->ncells;
    double rmin = this->rjmh[0];
    double rmax = this->rjmh[Nr];
    double zmin = this->zkmh[0];
    double zmax = this->zkmh[Nz];


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

    int i,j,k;

    int s = 0;
    for(k=0; k<Nz; k++)
    {
        double z = 0.5*(this->zkmh[k] + this->zkmh[k+1]);
        for(j=0; j<Nr; j++)
        {
            double r = 0.5*(this->rjmh[j]+this->rjmh[j+1]);

            int jk = Nr*k + j;

            for(i=0; i<this->Np[jk]; i++)
            {
                double phi = 0.5*(this->pimh[jk][i]+this->pimh[jk][i+1]);
                
                double vr = 0.1*(rand()*2.0/RAND_MAX-1);
                double vp = sqrt(1./r);
                double V[3] = {vr*cos(phi) - vp*sin(phi), vr*sin(phi)+vp*cos(phi), 0.0};

                double R = rmin + (rmax-rmin)*phi/(2*M_PI);
                double Z = zmin + (zmax-zmin)*phi/(2*M_PI);
                double sig = 0.1;
                double phase = -((r-R)*(r-R) +(z-Z)*(z-Z)) / (2*sig*sig);
                double rho = 1.0 + 1.0*exp(phase);
                double P = phi;
                double om = pow(r, -1.5);

                int c;
                for(c=0; c<this->ncellspercell[jk][i]; c++)
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

    myfile.open("/Users/geoff/Projects/DiscoPV/data2.txt");
    myfile << "line";
    myfile.close();
    
    return 1;
}

void vtkDiscoReader::MakeGridHexahedron(vtkUnstructuredGrid *output)
{
    double phimin = 0.0;
    double phimax;
    readSimple("Pars", "Phi_Max", &phimax, H5T_NATIVE_DOUBLE);
    hsize_t dims[3];
    getH5dims("Grid", "r_jph", dims);
    int Nr = dims[0] - 1;
    getH5dims("Grid", "z_kph", dims);
    int Nz = dims[0] - 1;

    this->Nz = Nz;
    this->Nr = Nr;
    this->phimax = phimax;

    this->rjmh.reserve(Nr+1);
    this->zkmh.reserve(Nz+1);
    this->Np.reserve(Nz*Nr);
    this->Index.reserve(Nz*Nr);
    this->Id_phi0.reserve(Nz*Nr);
    this->pimh.reserve(Nr*Nz);

    readSimple("Grid", "r_jph", this->rjmh.data(), H5T_NATIVE_DOUBLE);
    readSimple("Grid", "z_kph", this->zkmh.data(), H5T_NATIVE_DOUBLE);
    readSimple("Grid", "Np", this->Np.data(), H5T_NATIVE_INT);
    readSimple("Grid", "Index", this->Index.data(), H5T_NATIVE_INT);
    readSimple("Grid", "Id_phi0", this->Id_phi0.data(), H5T_NATIVE_INT);

    int i,j,k;

    srand(42);
    rand();
    rand();

    int jk;
    for(jk=0; jk<Nr*Nz; jk++)
        this->pimh.push_back(std::vector<double>(this->Np[jk]+1));

    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = Nr*k + j;
            double dphi = (phimax-phimin) / this->Np[jk];
            double p0 = -rand()*dphi/RAND_MAX;
            for(i=0; i<=this->Np[jk]; i++)
                this->pimh[jk][i] = p0 + i*dphi;
        }
    
    ofstream myfile;
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
    myfile << "Made grid\n";
    myfile.close();

    vtkPoints* points = vtkPoints::New();

    vtkIdType npointsTot = 0;
    vtkIdType npoints[Nz+1][Nr+1];
    vtkIdType indexPoints[Nz+1][Nr+1];
    std::vector<std::vector<vtkIdType>> cellPoints;
    cellPoints.reserve(Nr*Nz);

    for(jk=0; jk<Nr*Nz; jk++)
            cellPoints.push_back(std::vector<vtkIdType>(4*(this->Np[jk]+1)));

    for(k=0; k<=Nz; k++)
        for(j=0; j<=Nr; j++)
        {
            int jkDL = Nr*(k-1) + (j-1);
            int jkDR = Nr*(k-1) + j;
            int jkUL = Nr*k + (j-1);
            int jkUR = Nr*k + j;

            int np = 0;
            if(j>0 && k>0)
                np += this->Np[jkDL]+1;
            if(j<Nr && k>0)
                np += this->Np[jkDR]+1;
            if(j>0 && k<Nz)
                np += this->Np[jkUL]+1;
            if(j<Nr && k<Nz)
                np += this->Np[jkUR]+1;

            npoints[k][j] = np;
            indexPoints[k][j] = npointsTot;
            npointsTot += npoints[k][j];
        }
    points->SetNumberOfPoints(npointsTot);
    
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
    myfile << "Counted Points: " << npointsTot << "\n";
    myfile.close();

    for(k=0; k<=Nz; k++)
    {
        double z = this->zkmh[k];

        int kD = k-1;
        int kU = k;

        for(j=0; j<=Nr; j++)
        {
            int iDL = 0;
            int iDR = 0;
            int iUL = 0;
            int iUR = 0;

            double r = this->rjmh[j];

            int jL = j-1;
            int jR = j;

            int jkDL = Nr*kD + jL;
            int jkDR = Nr*kD + jR;
            int jkUL = Nr*kU + jL;
            int jkUR = Nr*kU + jR;

            double phiDL, phiDR, phiUL, phiUR;

            if(kD >= 0 && jL >= 0)
                phiDL = this->pimh[jkDL][0];
            else
                phiDL = HUGE;
            if(kD >= 0 && jR < Nr)
                phiDR = this->pimh[jkDR][0];
            else
                phiDR = HUGE;
            if(kU < Nz && jL >= 0)
                phiUL = this->pimh[jkUL][0];
            else
                phiUL = HUGE;
            if(kU < Nz && jR < Nr)
                phiUR = this->pimh[jkUR][0];
            else
                phiUR = HUGE;
    
            myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
            myfile << k << " " << j << "\n";
            myfile.close();

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

                double x[3] = {r*cos(*phimin), r*sin(*phimin), z};

                points->InsertPoint(s, x);
                cellPoints[jkmin][4*(*imin)+opp] = s;

                (*imin)++;
                if(*imin <= this->Np[jkmin])
                    *phimin = this->pimh[jkmin][*imin];
                else
                    *phimin = HUGE;
            }
        }
    }

    output->SetPoints(points);
    points->Delete();
    
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
    myfile << "Made Points\n";
    myfile.close();

    int ncellsTot = 0;
    this->ncellspercell.reserve(Nz*Nr);
    for(jk=0; jk<Nr*Nz; jk++)
        this->ncellspercell.push_back(std::vector<int>(this->Np[jk]));

    vtkCellArray *cells = vtkCellArray::New();

    for(k=0; k<Nz; k++)
        for(j=0; j<Nr; j++)
        {
            int jk = Nr*k + j;

            for(i=0; i<this->Np[jk]; i++)
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

                this->ncellspercell[jk][i] = nc;

                vtkIdType p[8] = {sDRB, sDLB, sULB, sURB, 
                            sDRB+1, sDLB+1, sULB+1, sURB+1};
                cells->InsertNextCell(8, p);

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
                }
            }
        }

    int nc = cells->GetNumberOfCells();
    this->ncells = nc;

    output->SetCells(VTK_HEXAHEDRON, cells);
    cells->Delete();
    
    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
    myfile << "Made Cells\n";
    myfile.close();
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


void vtkDiscoReader::getH5dims(const char *group, const char *dset, 
                                hsize_t *dims)
{
    hid_t h5fil = H5Fopen(this->FileName, H5F_ACC_RDWR , H5P_DEFAULT );
    hid_t h5grp = H5Gopen1( h5fil , group );
    hid_t h5dst = H5Dopen1( h5grp , dset );
    hid_t h5spc = H5Dget_space( h5dst );

    H5Sget_simple_extent_dims( h5spc , dims , NULL);

    H5Sclose( h5spc );
    H5Dclose( h5dst );
    H5Gclose( h5grp );
    H5Fclose( h5fil );
}

void vtkDiscoReader::readSimple(const char *group, const char *dset, 
                                void *data, hid_t type)
{
    hid_t h5fil = H5Fopen( this->FileName, H5F_ACC_RDWR, H5P_DEFAULT);
    hid_t h5grp = H5Gopen1( h5fil , group );
    hid_t h5dst = H5Dopen1( h5grp , dset );

    H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

    H5Dclose( h5dst );
    H5Gclose( h5grp );
    H5Fclose( h5fil );
}

void vtkDiscoReader::readPatch(const char *group, const char *dset, void *data, 
                                hid_t type, int dim, int *start, 
                                int *loc_size, int *glo_size)
{
    hid_t h5fil = H5Fopen( this->FileName , H5F_ACC_RDWR , H5P_DEFAULT );
    hid_t h5grp = H5Gopen1( h5fil , group );
    hid_t h5dst = H5Dopen1( h5grp , dset );

    hsize_t mdims[dim];
    hsize_t fdims[dim];

    hsize_t fstart[dim];
    hsize_t fstride[dim];
    hsize_t fcount[dim];
    hsize_t fblock[dim];

    int d;
    for( d=0 ; d<dim ; ++d )
    {
        mdims[d] = loc_size[d];
        fdims[d] = glo_size[d];

        fstart[d]  = start[d];
        fstride[d] = 1;
        fcount[d]  = loc_size[d];
        fblock[d]  = 1;
    }
    hid_t mspace = H5Screate_simple(dim,mdims,NULL);
    hid_t fspace = H5Screate_simple(dim,fdims,NULL);

    H5Sselect_hyperslab( fspace , H5S_SELECT_SET , fstart , fstride , fcount , 
                        fblock );

    H5Dread( h5dst , type , mspace , fspace , H5P_DEFAULT , data );

    H5Sclose( mspace );
    H5Sclose( fspace );
    H5Dclose( h5dst );
    H5Gclose( h5grp );
    H5Fclose( h5fil );
}

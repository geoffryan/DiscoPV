#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <vector>
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUnsignedCharArray.h"
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
    this->cellType = VTK_TETRA;
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

void vtkDiscoReader::SetMeshType(int type)
{
    if(type == 0)
        this->cellType = VTK_TETRA;
    else if(type == 1)
        this->cellType = VTK_POLYHEDRON;
    else if(type == 2)
        this->cellType = VTK_HEXAHEDRON;
    else
        this->cellType = VTK_TETRA;
    this->Modified();
}

int vtkDiscoReader::GetMeshType()
{
    if(this->cellType == VTK_TETRA)
        return 0;
    else if(this->cellType == VTK_POLYHEDRON)
        return 1;
    else if(this->cellType == VTK_HEXAHEDRON)
        return 2;
    else
        return -1;
}

void vtkDiscoReader::SetGeometry()
{
    char buf[256];
    buf[255] = '\0';

    this->readString("Opts", "GEOMETRY", buf, 256);

    std::cout << "Geometry: " << buf << "\n";

    if(     strncmp(buf, "cylindrical", 11) == 0)
        this->geometry = CYLINDRICAL;
    else if(strncmp(buf, "spherical", 9) == 0)
        this->geometry = SPHERICAL;
    else if(strncmp(buf, "cartesian", 9) == 0)
        this->geometry = CARTESIAN;
    else
        this->geometry = CYLINDRICAL;
}

int vtkDiscoReader::RequestInformation(vtkInformation*, vtkInformationVector**,
                                        vtkInformationVector* outVec)
{
    ofstream myfile;
    myfile.open("/Users/geoff/Projects/DiscoPV/info.txt");
    myfile << "line";
    myfile.close();

    vtkInformation *outInfo = outVec->GetInformationObject(0);

    outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

    double tRange[2] = {2.3, 3.1};
    double times[1];
    readSimple("Grid", "T", times, H5T_NATIVE_DOUBLE);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), times, 1);
    //outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), tRange, 2);

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

    int piece = outInfo->Get(
                    vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
    int numPieces = outInfo->Get(
                vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
    int ghostLevel = outInfo->Get(
            vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());

    double time = outInfo->Get(
                        vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

    std::cout << "Requested time: " << time << "\n";

    this->SetGeometry();
    int numPiecesUsed = this->LoadGridDims(piece, numPieces, ghostLevel);

    //Return Empty data if too many pieces requested.
    if(piece >= numPiecesUsed)
    {
        outData->Initialize();
        return 1;
    }

    this->LoadGrid();               // Read in Grid data into internal vars
    this->MakeGrid(outData);        // Assemble grid data into vtk object
    this->TagGhostCells(outData);    // Tag cells as ghosts.
    this->AddData(outData);         // Add hydro data.

    myfile.open("/Users/geoff/Projects/DiscoPV/data2.txt");
    myfile << "line";
    myfile.close();
    
    return 1;
}
    
int vtkDiscoReader::LoadGridDims(int piece, int numPieces, int ghostLevel)
{
    hsize_t hdims[1];
    std::cout << "Piece " << piece << " Getting r_jph dims\n";
    getH5dims("Grid", "r_jph", hdims);
    this->NrTot = hdims[0] - 1;
    std::cout << "Piece " << piece << " Getting z_kph dims\n";
    getH5dims("Grid", "z_kph", hdims);
    this->NzTot = hdims[0] - 1;

    if(numPieces > this->NrTot * this->NzTot)
        numPieces = this->NrTot * this->NzTot;

    int dims[2], rank[2];
    this->dimsCreate(numPieces, dims);
    rank[0] = piece % dims[0];
    rank[1] = piece / dims[0];


    this->jA = (rank[0]*this->NrTot)/dims[0];
    this->jB = ((rank[0]+1)*this->NrTot)/dims[0];
    if(rank[0] == dims[0]-1)
        this->jB = this->NrTot;
    
    this->kA = (rank[1]*this->NzTot)/dims[1];
    this->kB = ((rank[1]+1)*this->NzTot)/dims[1];
    if(rank[1] == dims[1]-1)
        this->kB = this->NzTot;

    this->jA = std::max(this->jA - (ghostLevel+1), 0);
    this->jB = std::min(this->jB + (ghostLevel+1), this->NrTot);
    this->kA = std::max(this->kA - (ghostLevel+1), 0);
    this->kB = std::min(this->kB + (ghostLevel+1), this->NzTot);

    if(piece >= numPieces)
    {
        this->jA = -1;
        this->jB = -1;
        this->kA = -1;
        this->kB = -1;
    }

    this->Nr = this->jB - this->jA;
    this->Nz = this->kB - this->kA;

    this->cja = this->jA==0 ? 0 : 1;
    this->cjb = this->jB==this->NrTot ? Nr : Nr-1;
    this->cka = this->kA==0 ? 0 : 1;
    this->ckb = this->kB==this->NzTot ? Nz : Nz-1;

    this->ghostLevel = ghostLevel;

    char fname[1024];
    sprintf(fname, "/Users/geoff/Projects/DiscoPV/dims.%03d.txt", piece);

    ofstream myfile;
    myfile.open(fname);
    myfile << "Piece "<<piece<<" of " << numPieces << "\n";
    myfile << "Piece "<<piece<<" FileName: " << this->FileName << "\n";
    myfile << "Piece "<<piece<<" Dims: " << dims[0] << " x " << dims[1] << "\n";
    myfile << "Piece "<<piece<<" Rank: " << rank[0] << " x " << rank[1] << "\n";
    myfile << "Piece "<<piece<<" NrTot: " << this->NrTot << "\n";
    myfile << "Piece "<<piece<<" NzTot: " << this->NzTot << "\n";
    myfile << "Piece "<<piece<<" Nr: " << this->Nr << "\n";
    myfile << "Piece "<<piece<<" Nz: " << this->Nz << "\n";
    myfile << "Piece "<<piece<<" ghostLevels: " << ghostLevel << "\n";
    myfile << "Piece "<<piece<<" jAB: " << this->jA << " " << this->jB << "\n";
    myfile << "Piece "<<piece<<" kAB: " << this->kA << " " << this->kB << "\n";
    myfile << "Piece "<<piece<<" cjab: " << this->cja << " " << this->cjb << "\n";
    myfile << "Piece "<<piece<<" ckab: " << this->cka << " " << this->ckb << "\n";
    myfile.close();
    
    std::cout << "Piece " << piece << " GridDims Loaded.\n";

    return numPieces;
}

void vtkDiscoReader::MakeGrid(vtkUnstructuredGrid *output)
{

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

    int i, j, k, jk;

    int pkA = this->kA==0 ? 0 : 1;
    int pkB = this->kB==this->NzTot ? Nz : Nz-1;
    int pjA = this->jA==0 ? 0 : 1;
    int pjB = this->jB==this->NrTot ? Nr : Nr-1;

    for(jk=0; jk<Nr*Nz; jk++)
        cellPoints.push_back(std::vector<vtkIdType>(4*(this->Np[jk]+1)));

    for(k=this->cka; k<=this->ckb; k++)
        for(j=this->cja; j<=this->cjb; j++)
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

    for(k=this->cka; k<=this->ckb; k++)
    {
        double z = this->zkmh[k];

        int kD = k-1;
        int kU = k;

        for(j=this->cja; j<=this->cjb; j++)
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

            for(i=0; i<npoints[k][j]; i++)
            {
                vtkIdType s = indexPoints[k][j] + i;

                int jkmin, opp;
                int *imin;
                double *phimin;

                int w4 = this->which4(phiDL, phiDR, phiUL, phiUR);
                if(w4 == 0)
                {
                    imin = &iDL;
                    jkmin = jkDL;
                    phimin = &phiDL;
                    opp = 3;
                }
                else if(w4 == 1)
                {
                    imin = &iDR;
                    jkmin = jkDR;
                    phimin = &phiDR;
                    opp = 2;
                }
                else if(w4 == 2)
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

                double x[3];
                this->calcXYZ(r, *phimin, z, x);

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

    if(this->cellType == VTK_TETRA)
        this->AddCellsTetrahedral(output, cellPoints);
    else if(this->cellType == VTK_POLYHEDRON)
        this->AddCellsPolyhedral(output, cellPoints);
    else
        this->AddCells(output, cellPoints);

    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
    myfile << "Added cells: " << this->ncells << "\n";
    myfile.close();

    myfile.open("/Users/geoff/Projects/DiscoPV/data.txt", 
                    std::ios_base::app);
    myfile << "Made Cells\n";
    myfile.close();
}

void vtkDiscoReader::AddCells(vtkUnstructuredGrid *output,
                        const std::vector<std::vector<vtkIdType>> &cellPoints)
{
    /*
     * Rules for Cells to make a water-tight mesh:
     *    -cells MUST have edges at neighbour's phi-boundaries.
     *    -adjacent faces MUST be tesselated the same way.
     *
     *    This mesh does not follow these rules.
     */

    int Nr = this->Nr;
    int Nz = this->Nz;

    int i,j,k,jk;

    int ncellsTot = 0;
    this->ncellspercell.reserve(Nz*Nr);

    for(jk=0; jk<Nr*Nz; jk++)
        this->ncellspercell.push_back(std::vector<int>(this->Np[jk]));

    output->Allocate();

    for(k=this->cka; k<this->ckb; k++)
    {
        for(j=this->cja; j<this->cjb; j++)
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

                int nadd = 0;
                nadd += this->AddDiscoCell(output, 
                                                sDLB, sDRB, sULB, sURB,
                                                sDLF, sDRF, sULF, sURF);
                this->ncellspercell[jk][i] = nadd;
            }
        }
    }

    int nc = output->GetNumberOfCells();
    this->ncells = nc;
}

void vtkDiscoReader::AddCellsTetrahedral(vtkUnstructuredGrid *output,
                        const std::vector<std::vector<vtkIdType>> &cellPoints)
{
    /*
     * Rules for Cells to make a water-tight mesh:
     *    -cells MUST have edges at neighbour's phi-boundaries.
     *    -adjacent faces MUST be tesselated the same way.
     */

    int Nr = this->Nr;
    int Nz = this->Nz;

    int i,j,k,jk;

    int ncellsTot = 0;
    this->ncellspercell.reserve(Nz*Nr);

    for(jk=0; jk<Nr*Nz; jk++)
        this->ncellspercell.push_back(std::vector<int>(this->Np[jk]));

    output->Allocate();


    for(k=this->cka; k<this->ckb; k++)
    {
        int kD = k-1;
        int kU = k+1;

        for(j=this->cja; j<this->cjb; j++)
        {
            int jk = Nr*k + j;

            int jL = j-1;
            int jR = j+1;

            int iL = 0;
            int iR = 0;
            int iD = 0;
            int iU = 0;

            int jkL = Nr*k + jL;
            int jkR = Nr*k + jR;
            int jkD = Nr*kD + j;
            int jkU = Nr*kU + j;

            double phiL, phiR, phiD, phiU;
                
            if(j>0)
            {
                while(this->pimh[jkL][iL] < this->pimh[jk][0])
                    iL++;
                phiL = this->pimh[jkL][iL];
            }
            else
                phiL = HUGE;
            if(j<Nr-1)
            {
                while(this->pimh[jkR][iR] < this->pimh[jk][0])
                    iR++;
                phiR = this->pimh[jkR][iR];
            }
            else
                phiR = HUGE;
            if(k>0)
            {
                while(this->pimh[jkD][iD] < this->pimh[jk][0])
                    iD++;
                phiD = this->pimh[jkD][iD];
            }
            else
                phiD = HUGE;
            if(k<Nz-1)
            {
                while(this->pimh[jkU][iU] < this->pimh[jk][0])
                    iU++;
                phiU = this->pimh[jkU][iU];
            }
            else
                phiU = HUGE;

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

                int nadd = 0;

                double phiF = this->pimh[jk][i+1];
                vtkIdType s0b, s1b, s2b;
                vtkIdType p[4];

                //Bottom-left Loop
                
                p[0] = sDRB;
                p[1] = sDLB;
                p[2] = sULB;
                
                while(phiL < phiF || phiD < phiF)
                {
                    s0b = phiD<phiF ? cellPoints[jkD][4*iD+3] : sDRF;
                    s2b = phiL<phiF ? cellPoints[jkL][4*iL+3] : sULF;
                    if(phiL < phiD)
                        s1b = phiL<phiF ? cellPoints[jkL][4*iL+1] : sDLF;
                    else
                        s1b = phiD<phiF ? cellPoints[jkD][4*iD+2] : sDLF;

                    while(p[0] < s0b || p[1] < s1b || p[2] < s2b)
                    {
                        //Add L/D points before R/U.
                        //DL
                        if(p[1] < s1b)
                        {
                            p[3] = p[1]+1;
                            output->InsertNextCell(VTK_TETRA, 4, p);
                            nadd++;
                            p[1]++;
                        }
                        //DR
                        if(p[0] < s0b)
                        {
                            p[3] = p[0]+1;
                            output->InsertNextCell(VTK_TETRA, 4, p);
                            nadd++;
                            p[0]++;
                        }
                        //UL
                        if(p[2] < s2b)
                        {
                            p[3] = p[2]+1;
                            output->InsertNextCell(VTK_TETRA, 4, p);
                            nadd++;
                            p[2]++;
                        }
                    }

                    if(phiL < phiD)
                    {
                        iL++;
                        if(iL <= this->Np[jkL])
                            phiL = this->pimh[jkL][iL];
                        else
                            phiL = HUGE;
                    }
                    else
                    {
                        iD++;
                        if(iD <= this->Np[jkD])
                            phiD = this->pimh[jkD][iD];
                        else
                            phiD = HUGE;
                    }
                }
                s0b = sDRF;
                s1b = sDLF;
                s2b = sULF;

                while(p[0] < s0b || p[1] < s1b || p[2] < s2b)
                {
                    //Add L/D points before R/U.
                    //DL
                    if(p[1] < s1b)
                    {
                        p[3] = p[1]+1;
                        output->InsertNextCell(VTK_TETRA, 4, p);
                        nadd++;
                        p[1]++;
                    }
                    //DR
                    if(p[0] < s0b)
                    {
                        p[3] = p[0]+1;
                        output->InsertNextCell(VTK_TETRA, 4, p);
                        nadd++;
                        p[0]++;
                    }
                    //UL
                    if(p[2] < s2b)
                    {
                        p[3] = p[2]+1;
                        output->InsertNextCell(VTK_TETRA, 4, p);
                        nadd++;
                        p[2]++;
                    }
                }

                //Upper-right Loop

                p[0] = sULB;
                p[1] = sURB;
                p[2] = sDRB;
                
                while(phiR < phiF || phiU < phiF)
                {
                    s0b = phiU<phiF ? cellPoints[jkU][4*iU+0] : sULF;
                    s2b = phiR<phiF ? cellPoints[jkR][4*iR+0] : sDRF;
                    if(phiR < phiU)
                        s1b = phiR<phiF ? cellPoints[jkR][4*iR+2] : sURF;
                    else
                        s1b = phiU<phiF ? cellPoints[jkU][4*iU+1] : sURF;

                    while(p[0] < s0b || p[1] < s1b || p[2] < s2b)
                    {
                        //Add L/D points before R/U.
                        //UL
                        if(p[0] < s0b)
                        {
                            p[3] = p[0]+1;
                            output->InsertNextCell(VTK_TETRA, 4, p);
                            nadd++;
                            p[0] = p[3];
                        }
                        //DR
                        if(p[2] < s2b)
                        {
                            p[3] = p[2]+1;
                            output->InsertNextCell(VTK_TETRA, 4, p);
                            nadd++;
                            p[2] = p[3];
                        }
                        //UR
                        if(p[1] < s1b)
                        {
                            p[3] = p[1]+1;
                            output->InsertNextCell(VTK_TETRA, 4, p);
                            nadd++;
                            p[1] = p[3];
                        }
                    }

                    if(phiR < phiF)
                    {
                        iR++;
                        if(iR <= this->Np[jkR])
                            phiR = this->pimh[jkR][iR];
                        else
                            phiR = HUGE;
                    }
                    else
                    {
                        iU++;
                        if(iU <= this->Np[jkU])
                            phiU = this->pimh[jkU][iU];
                        else
                            phiU = HUGE;
                    }
                }

                s0b = sULF;
                s1b = sURF;
                s2b = sDRF;

                while(p[0] < s0b || p[1] < s1b || p[2] < s2b)
                {
                    //Add L/D points before R/U.
                    //UL
                    if(p[0] < s0b)
                    {
                        p[3] = p[0]+1;
                        output->InsertNextCell(VTK_TETRA, 4, p);
                        nadd++;
                        p[0] = p[3];
                    }
                    //DR
                    if(p[2] < s2b)
                    {
                        p[3] = p[2]+1;
                        output->InsertNextCell(VTK_TETRA, 4, p);
                        nadd++;
                        p[2] = p[3];
                    }
                    //UR
                    if(p[1] < s1b)
                    {
                        p[3] = p[1]+1;
                        output->InsertNextCell(VTK_TETRA, 4, p);
                        nadd++;
                        p[1] = p[3];
                    }
                }
                
                this->ncellspercell[jk][i] = nadd;
            }
        }
    }

    int nc = output->GetNumberOfCells();
    this->ncells = nc;
}

void vtkDiscoReader::AddCellsPolyhedral(vtkUnstructuredGrid *output,
                        const std::vector<std::vector<vtkIdType>> &cellPoints)
{
    /*
     * Rules for Cells to make a water-tight mesh:
     *    -cells MUST have edges at neighbour's phi-boundaries.
     *    -adjacent faces MUST be tesselated the same way.
     */

    int Nr = this->Nr;
    int Nz = this->Nz;

    int i,j,k,jk;

    int ncellsTot = 0;
    this->ncellspercell.reserve(Nz*Nr);

    for(jk=0; jk<Nr*Nz; jk++)
        this->ncellspercell.push_back(std::vector<int>(this->Np[jk]));

    output->Allocate();

    std::vector<vtkIdType> faceStream;

    for(k=this->cka; k<this->ckb; k++)
    {
        int kD = k-1;
        int kU = k+1;

        for(j=this->cja; j<this->cjb; j++)
        {
            int jk = Nr*k + j;

            int jL = j-1;
            int jR = j+1;

            int iL = 0;
            int iR = 0;
            int iD = 0;
            int iU = 0;

            int jkL = Nr*k + jL;
            int jkR = Nr*k + jR;
            int jkD = Nr*kD + j;
            int jkU = Nr*kU + j;

            double phiL, phiR, phiD, phiU;
                
            if(j>0)
            {
                while(this->pimh[jkL][iL] < this->pimh[jk][0])
                    iL++;
                phiL = this->pimh[jkL][iL];
            }
            else
                phiL = HUGE;
            if(j<Nr-1)
            {
                while(this->pimh[jkR][iR] < this->pimh[jk][0])
                    iR++;
                phiR = this->pimh[jkR][iR];
            }
            else
                phiR = HUGE;
            if(k>0)
            {
                while(this->pimh[jkD][iD] < this->pimh[jk][0])
                    iD++;
                phiD = this->pimh[jkD][iD];
            }
            else
                phiD = HUGE;
            if(k<Nz-1)
            {
                while(this->pimh[jkU][iU] < this->pimh[jk][0])
                    iU++;
                phiU = this->pimh[jkU][iU];
            }
            else
                phiU = HUGE;

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

                int nf = 0;
                faceStream.clear();
                
                faceStream.push_back(4);
                faceStream.push_back(sDRB);
                faceStream.push_back(sDLB);
                faceStream.push_back(sULB);
                faceStream.push_back(sURB);
                faceStream.push_back(4);
                faceStream.push_back(sDRF);
                faceStream.push_back(sDLF);
                faceStream.push_back(sULF);
                faceStream.push_back(sURF);
                nf += 2;

                double phiF = this->pimh[jk][i+1];
                vtkIdType sLa, sRa, sLb, sRb;

                //Bottom Loop
                sLa = sDLB;
                sRa = sDRB;
                
                while(phiD < phiF)
                {
                    sLb = cellPoints[jkD][4*iD+2];
                    sRb = cellPoints[jkD][4*iD+3];
                    nf += this->AddFacesToStream(faceStream, 
                                                    sLa, sRa, sLb, sRb);
                    sLa = sLb;
                    sRa = sRb;

                    iD++;
                    if(iD <= this->Np[jkD])
                        phiD = this->pimh[jkD][iD];
                    else
                        phiD = HUGE;
                }
                sLb = sDLF;
                sRb = sDRF;
                nf += this->AddFacesToStream(faceStream, sLa, sRa, sLb, sRb);

                //Top Loop
                sLa = sULB;
                sRa = sURB;
                
                while(phiU < phiF)
                {
                    sLb = cellPoints[jkU][4*iU+0];
                    sRb = cellPoints[jkU][4*iU+1];
                    nf += this->AddFacesToStream(faceStream, 
                                                    sLa, sRa, sLb, sRb);
                    sLa = sLb;
                    sRa = sRb;

                    iU++;
                    if(iU <= this->Np[jkU])
                        phiU = this->pimh[jkU][iU];
                    else
                        phiU = HUGE;
                }
                sLb = sULF;
                sRb = sURF;
                nf += this->AddFacesToStream(faceStream, sLa, sRa, sLb, sRb);
               
                // Radial Inner Loop
                sLa = sULB;
                sRa = sDLB;
                
                while(phiL < phiF)
                {
                    sLb = cellPoints[jkL][4*iL+3];
                    sRb = cellPoints[jkL][4*iL+1];
                    nf += this->AddFacesToStream(faceStream, 
                                                    sLa, sRa, sLb, sRb);
                    sLa = sLb;
                    sRa = sRb;

                    iL++;
                    if(iL <= this->Np[jkL])
                        phiL = this->pimh[jkL][iL];
                    else
                        phiL = HUGE;
                }
                sLb = sULF;
                sRb = sDLF;
                nf += this->AddFacesToStream(faceStream, sLa, sRa, sLb, sRb);
                
                //Outer Radial Loop
                sLa = sURB;
                sRa = sDRB;
                
                while(phiR < phiF)
                {
                    sLb = cellPoints[jkR][4*iR+2];
                    sRb = cellPoints[jkR][4*iR+0];
                    nf += this->AddFacesToStream(faceStream, 
                                                    sLa, sRa, sLb, sRb);
                    sLa = sLb;
                    sRa = sRb;

                    iR++;
                    if(iR <= this->Np[jkR])
                        phiR = this->pimh[jkR][iR];
                    else
                        phiR = HUGE;
                }
                sLb = sURF;
                sRb = sDRF;
                nf += this->AddFacesToStream(faceStream, sLa, sRa, sLb, sRb);

                output->InsertNextCell(VTK_POLYHEDRON, nf, faceStream.data());
                this->ncellspercell[jk][i] = 1;
            }
        }
    }

    int nc = output->GetNumberOfCells();
    this->ncells = nc;
}

void vtkDiscoReader::TagGhostCells(vtkUnstructuredGrid *output)
{
    int i,j,k,c;

    int ghostLevel = this->ghostLevel;
    int cka = this->cka;
    int ckb = this->ckb;
    int cja = this->cja;
    int cjb = this->cjb;

    vtkUnsignedCharArray *ghosts = output->GetCellGhostArray();
    if(!ghosts)
        ghosts = output->AllocateCellGhostArray();

    vtkIdType s = 0;
    for(k=cka; k<ckb; k++)
        for(j=cja; j<cjb; j++)
        {
            int jk = k * this->Nr + j;

            for(i=0; i<this->Np[jk]; i++)
                for(c=0; c<this->ncellspercell[jk][i]; c++)
                {
                    unsigned char val = ghosts->GetValue(s);
                    if(     (jA > 0             && j < cja + ghostLevel)
                        ||  (jB < this->NrTot   && j >= cjb - ghostLevel)
                        ||  (kA > 0             && k < cka + ghostLevel)
                        ||  (kB < this->NzTot   && k >= ckb - ghostLevel))
                        val = val | vtkDataSetAttributes::DUPLICATECELL;
                    ghosts->SetValue(s, val);
                    s++;
                }
        }
}

void vtkDiscoReader::LoadGrid()
{
    int Nr = this->Nr;
    int Nz = this->Nz;

    this->rjmh.reserve(Nr+1);
    this->zkmh.reserve(Nz+1);
    this->Np.reserve(Nz*Nr);
    this->Index.reserve(Nz*Nr);
    this->Id_phi0.reserve(Nz*Nr);

    hsize_t start[2] = {0,0};
    hsize_t loc_size[2] = {0,0};
    hsize_t glo_size[2] = {0,0};

    start[0] = this->jA;
    loc_size[0] = this->Nr+1;
    glo_size[0] = this->NrTot+1;
    readPatch("Grid", "r_jph", this->rjmh.data(), H5T_NATIVE_DOUBLE,
                1, start, loc_size, glo_size);
    start[0] = this->kA;
    loc_size[0] = this->Nz+1;
    glo_size[0] = this->NzTot+1;
    readPatch("Grid", "z_kph", this->zkmh.data(), H5T_NATIVE_DOUBLE,
                1, start, loc_size, glo_size);

    start[0] = this->kA;
    start[1] = this->jA;
    loc_size[0] = this->Nz;
    loc_size[1] = this->Nr;
    glo_size[0] = this->NzTot;
    glo_size[1] = this->NrTot;
    readPatch("Grid", "Np", this->Np.data(), H5T_NATIVE_INT,
                2, start, loc_size, glo_size);
    readPatch("Grid", "Index", this->Index.data(), H5T_NATIVE_INT,
                2, start, loc_size, glo_size);
    readPatch("Grid", "Id_phi0", this->Id_phi0.data(), H5T_NATIVE_INT,
                2, start, loc_size, glo_size);

    //Compute "chunk" size. This is the number of contiguous cells stored
    //forwards of Index[jk].  It is at least Np[jk], but may be larger if we're
    //lucky.  Only computed for the cells we know we'll need to read.

    this->chunkSize.reserve(Nr*Nz);
    int i,j,k,jk;
    for(k=this->cka; k < this->ckb; k++)
        for(j=this->cja; j < this->cjb; j++)
        {
            jk = k*Nr + j;
            this->chunkSize[jk] = Np[jk];
            int next = this->Index[jk] + this->Np[jk];
            int k1, j1;
            bool end = false;
            for(j1=j+1; j1<this->cjb; j1++)
            {
                int jk1 = k*Nr+j1;
                if(this->Index[jk1] == next)
                {
                    this->chunkSize[jk] += this->Np[jk1];
                    next += this->Np[jk1];
                }
                else
                {
                    end = true;
                    break;
                }
            }
            if(!end)
                for(k1=k+1; k1<this->ckb; k1++)
                {
                    for(j1=this->cja; j1<this->cjb; j1++)
                    {
                        int jk1 = k*Nr+j1;
                        if(this->Index[jk1] == next)
                        {
                            this->chunkSize[jk] += this->Np[jk1];
                            next += this->Np[jk1];
                        }
                        else
                        {
                            end = true;
                            break;
                        }
                    }
                    if(end)
                        break;
                }

        }


    double phimax;
    readSimple("Pars", "Phi_Max", &phimax, H5T_NATIVE_DOUBLE);
    this->phimax = phimax;
    
    this->pimh.reserve(Nr*Nz);

    for(jk=0; jk<Nr*Nz; jk++)
        this->pimh.push_back(std::vector<double>(this->Np[jk]+1));

    hsize_t dims[2];
    this->getH5dims("Data", "Cells", dims);

    start[0] = 0;
    start[1] = dims[1]-1;
    loc_size[0] = 0;
    loc_size[1] = 1;
    glo_size[0] = dims[0];
    glo_size[1] = dims[1];
    std::vector<double> buf;

    for(jk=0; jk<Nr*Nz; jk++)
    {
        int sa = this->Index[jk];
        int s0 = this->Id_phi0[jk];
        int i0 = s0-sa;
        int n = this->Np[jk];

        start[0] = sa;
        loc_size[0] = n;
        buf.reserve(n);

        this->readPatch("Data", "Cells", buf.data(), H5T_NATIVE_DOUBLE, 
                            2, start, loc_size, glo_size);

        for(i=0; i<n-i0; i++)
            this->pimh[jk][i+1] = buf[i0 + i];

        for(i=n-i0; i<n; i++)
            this->pimh[jk][i+1] = buf[i - (n-i0)];

        this->pimh[jk][0] = this->pimh[jk][n] - phimax;
    }
}

void vtkDiscoReader::AddData(vtkUnstructuredGrid *output)
{
    int num_c, num_n;
    this->readSimple("Opts", "NUM_C", &num_c, H5T_NATIVE_INT);
    this->readSimple("Opts", "NUM_N", &num_n, H5T_NATIVE_INT);

    int q;
    this->AddDataScalar(output, "rho", 0);
    this->AddDataScalar(output, "P", 1);
    if(this->geometry == CYLINDRICAL)
    {
        this->AddDataScalar(output, "vr", 2);
        this->AddDataScalar(output, "om", 3);
        this->AddDataScalar(output, "vz", 4);
    }
    else if(this->geometry == SPHERICAL)
    {
        this->AddDataScalar(output, "vr", 2);
        this->AddDataScalar(output, "om", 3);
        this->AddDataScalar(output, "vt", 4);
    }
    else if(this->geometry == CARTESIAN)
    {
        this->AddDataScalar(output, "vx", 2);
        this->AddDataScalar(output, "vy", 3);
        this->AddDataScalar(output, "vz", 4);
    }
    this->AddDataVector(output, "V", 2, 1);

    if(num_c == 8)
    {
        if(this->geometry == CYLINDRICAL)
        {
            this->AddDataScalar(output, "Br", 5);
            this->AddDataScalar(output, "Bp", 6);
            this->AddDataScalar(output, "Bz", 7);
        }
        else if(this->geometry == SPHERICAL)
        {
            this->AddDataScalar(output, "Br", 5);
            this->AddDataScalar(output, "Bp", 6);
            this->AddDataScalar(output, "Bt", 7);
        }
        else if(this->geometry == CARTESIAN)
        {
            this->AddDataScalar(output, "Bx", 5);
            this->AddDataScalar(output, "By", 6);
            this->AddDataScalar(output, "Bz", 7);
        }
        this->AddDataVector(output, "B", 5, 2);
    }
    char qname[128];
    for(q=0; q<num_n; q++)
    {
        sprintf(qname, "q%d", q);
        this->AddDataScalar(output, qname, num_c+q);
    }
}

void vtkDiscoReader::AddDataScalar(vtkUnstructuredGrid *output, 
                                    const char *name, int id)
{
    hsize_t dims[2];

    this->getH5dims("Data", "Cells", dims);

    std::vector<double> buf;

    hsize_t start[2] = {0, (hsize_t)id};
    hsize_t loc_size[2] = {0, 1};
    hsize_t glo_size[2] = {dims[0], dims[1]};

    vtkDoubleArray *q = vtkDoubleArray::New();
    q->SetName(name);
    q->SetNumberOfValues(this->ncells);

    int i, j, k, s, c, b0;

    s = 0;
    b0 = -1;
    for(k=this->cka; k<this->ckb; k++)
        for(j=this->cja; j<this->cjb; j++)
        {
            int jk = k * this->Nr + j;

            int sa = this->Index[jk];
            int s0 = this->Id_phi0[jk];
            int i0 = s0-sa;
            int n = this->Np[jk];

            if(b0 >= loc_size[0] || b0 < 0)
            {
                start[0] = sa;
                loc_size[0] = this->chunkSize[jk];
                buf.reserve(this->chunkSize[jk]);
                b0 = 0;
                this->readPatch("Data", "Cells", buf.data(), H5T_NATIVE_DOUBLE, 
                                    2, start, loc_size, glo_size);
            }

            for(i=0; i<n; i++)
            {
                double val;
                if(i < n-i0)
                    val = buf[b0+i0+i];
                else 
                    val = buf[b0+i-(n-i0)];
                for(c=0; c<this->ncellspercell[jk][i]; c++)
                {
                    q->SetValue(s, val);
                    s++;
                }
            }
            b0 += n;
        }
    output->GetCellData()->AddArray(q);
    q->Delete();
}

void vtkDiscoReader::AddDataVector(vtkUnstructuredGrid *output, 
                                    const char *name, int id0, int basis)
{
    hsize_t dims[2];

    this->getH5dims("Data", "Cells", dims);

    std::vector<double> buf;

    hsize_t start[2] = {0, (hsize_t)id0};
    hsize_t loc_size[2] = {0, 3};
    hsize_t glo_size[2] = {dims[0], dims[1]};


    vtkDoubleArray *q = vtkDoubleArray::New();
    q->SetName(name);
    q->SetNumberOfComponents(3);
    q->SetNumberOfTuples(this->ncells);

    int i, j, k, s, c, b0;
    double *buf0;
    int Nr = this->Nr;

    s = 0;
    b0 = -1;
    for(k=this->cka; k<this->ckb; k++)
    {
        double z = 0.5*(this->zkmh[k] + this->zkmh[k+1]);

        for(j=this->cja; j<this->cjb; j++)
        {
            int jk = Nr*k + j;

            int sa = this->Index[jk];
            int s0 = this->Id_phi0[jk];
            int i0 = s0-sa;
            int n = this->Np[jk];

            double r = 0.5*(this->rjmh[j] + this->rjmh[j+1]);

            if(b0 >= loc_size[0]*loc_size[1] || b0 < 0)
            {
                start[0] = sa;
                loc_size[0] = this->chunkSize[jk];
                buf.reserve(loc_size[0]*loc_size[1]);
                this->readPatch("Data", "Cells", buf.data(), H5T_NATIVE_DOUBLE, 
                                    2, start, loc_size, glo_size);
                buf0 = buf.data();
                b0 = 0;
            }

            for(i=0; i<n; i++)
            {
                double *vec;
                if(i < n-i0)
                    vec = buf0 + b0 + 3*(i0+i);
                else
                    vec = buf0 + b0 + 3*(i-(n-i0));
                double phi = 0.5*(this->pimh[jk][i] + this->pimh[jk][i+1]);
                double Vxyz[3];
                if(basis == 1)
                    this->calcVFromContravariant(r, phi, z, vec, Vxyz);
                else if(basis == 2)
                    this->calcVFromOrthonormal(r, phi, z, vec, Vxyz);
                else
                    for(c=0; c<3; c++)
                        Vxyz[c] = vec[c];

                for(c=0; c<this->ncellspercell[jk][i]; c++)
                {
                    q->SetTuple(s, Vxyz);
                    s++;
                }
            }

            b0 += 3*n;
        }
    }
    output->GetCellData()->AddArray(q);
    q->Delete();
}

int vtkDiscoReader::AddDiscoCell(vtkUnstructuredGrid *output,
                                    vtkIdType sDLB, vtkIdType sDRB, 
                                    vtkIdType sULB, vtkIdType sURB, 
                                    vtkIdType sDLF, vtkIdType sDRF, 
                                    vtkIdType sULF, vtkIdType sURF)
{
    int nc = 0;

    if(this->cellType == VTK_HEXAHEDRON)
    {
        vtkIdType p[8] = {sDRB, sDLB, sULB, sURB, 
                    sDRB+1, sDLB+1, sULB+1, sURB+1};
        output->InsertNextCell(this->cellType, 8, p);
        nc++;

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
            output->InsertNextCell(this->cellType, 8, p);
            nc++;
        }
    }
    else if(this->cellType == VTK_POLYHEDRON)
    {
        std::vector<vtkIdType> faceStream;
        int nf = 0;

        //Phi Faces
        faceStream.push_back(4);
        faceStream.push_back(sDRB);
        faceStream.push_back(sDLB);
        faceStream.push_back(sULB);
        faceStream.push_back(sURB);
        faceStream.push_back(4);
        faceStream.push_back(sDRF);
        faceStream.push_back(sDLF);
        faceStream.push_back(sULF);
        faceStream.push_back(sURF);
        nf += 2;
        //Inner Radial Faces
        nf += this->AddFacesToStream(faceStream, sULB, sDLB, sULF, sDLF);
        //Outer Radial Faces
        nf += this->AddFacesToStream(faceStream, sURB, sDRB, sURF, sDRF);
        //Bottom Z Faces
        nf += this->AddFacesToStream(faceStream, sDLB, sDRB, sDLF, sDRF);
        //Top Z Faces
        nf += this->AddFacesToStream(faceStream, sULB, sURB, sULF, sURF);

        output->InsertNextCell(this->cellType, nf, faceStream.data());
        nc++;
    }

    return nc;
}

int vtkDiscoReader::AddFacesToStream(std::vector<vtkIdType> &faceStream,
                                        vtkIdType sBL, vtkIdType sBR,
                                        vtkIdType sFL, vtkIdType sFR)
{
    int nf = 0;
    vtkIdType p[3];
    p[0] = sBL;
    p[1] = sBR;
    p[2] = sBR+1;
    faceStream.push_back(3);
    faceStream.push_back(p[0]);
    faceStream.push_back(p[1]);
    faceStream.push_back(p[2]);
    nf++;
    p[1] = p[2];

    while(p[0] < sFL || p[1] < sFR)
    {
        if(p[0] < sFL)
        {
            p[2] = p[0]+1;
            faceStream.push_back(3);
            faceStream.push_back(p[0]);
            faceStream.push_back(p[1]);
            faceStream.push_back(p[2]);
            nf++;
            p[0] = p[2];
        }

        if(p[1] < sFR)
        {
            p[2] = p[1]+1;
            faceStream.push_back(3);
            faceStream.push_back(p[0]);
            faceStream.push_back(p[1]);
            faceStream.push_back(p[2]);
            nf++;
            p[1] = p[2];
        }
    }

    return nf;
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
        x23 = x2;
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

void vtkDiscoReader::order4(double x0, double x1, double x2, double x3,
                            int *i0, int *i1, int *i2, int *i3)
{
    int lo0, lo1, lo2, lo3;
    int lo01, lo23, mida, midb;
    int hi01, hi23, hia, hib;
    double xlo01, xlo23, xmida, xmidb;
    double xhi01, xhi23, xhia, xhib;

    if(x0 < x1)
    {
        lo01 = 0;
        xlo01 = x0;
        hi01 = 1;
        xhi01 = x1;
    }
    else
    {
        lo01 = 1;
        xlo01 = x1;
        hi01 = 0;
        xhi01 = x0;
    }

    if(x2 < x3)
    {
        lo23 = 2;
        xlo23 = x2;
        hi23 = 3;
        xhi23 = x3;
    }
    else
    {
        lo23 = 3;
        xlo23 = x3;
        hi23 = 2;
        xhi23 = x2;
    }

    if(xlo01 < xlo23)
    {
        lo0 = lo01;
        mida = lo23;
        xmida = xlo23;
        midb = hi01;
        xmidb = xhi01;
        hia = hi23;
        xhia = xhi23;
    }
    else
    {
        lo0 = lo23;
        mida = lo01;
        xmida = xlo01;
        midb = hi23;
        xmidb = xhi23;
        hia = hi01;
        xhia = xhi01;
    }

    if(xmida < xmidb)
    {
        lo1 = mida;
        hib = midb;
        xhib = xmidb;
    }
    else
    {
        lo1 = midb;
        hib = mida;
        xhib = xmida;
    }

    if(xhia < xhib)
    {
        lo2 = hia;
        lo3 = hib;
    }
    else
    {
        lo2 = hib;
        lo3 = hia;
    }

    *i0 = lo0;
    *i1 = lo1;
    *i2 = lo2;
    *i3 = lo3;
}



void vtkDiscoReader::getH5dims(const char *group, const char *dset, 
                                hsize_t *dims)
{
    hid_t h5fil = H5Fopen(this->FileName, H5F_ACC_RDONLY , H5P_DEFAULT );
    if(h5fil < 0)
    {
        std::cout<< "ERR - getH5dims - Could Not Read File: " << this->FileName << "\n";
        dims[0] = -1;
        return;
    }
    hid_t h5grp = H5Gopen1( h5fil , group );
    if(h5grp < 0)
    {
        std::cout<< "ERR - getH5dims - Could Not Open Group: " << group << "\n";
        dims[0] = -1;
        H5Fclose(h5fil);
        return;
    }
    hid_t h5dst = H5Dopen1( h5grp , dset );
    if(h5dst < 0)
    {
        std::cout<< "ERR - getH5dims - Could Not Open Dataset: " << dset << "\n";
        dims[0] = -1;
        H5Gclose(h5grp);
        H5Fclose(h5fil);
        return;
    }
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
    hid_t h5fil = H5Fopen( this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    if(h5fil < 0)
    {
        std::cout<< "ERR - readSimple - Could Not Read File: " << this->FileName << "\n";
        data = NULL;
        return;
    }
    hid_t h5grp = H5Gopen1( h5fil , group );
    if(h5grp < 0)
    {
        std::cout<< "ERR - readSimple - Could Not Open Group: " << group << "\n";
        data = NULL;
        H5Fclose(h5fil);
        return;
    }
    hid_t h5dst = H5Dopen1( h5grp , dset );
    if(h5dst < 0)
    {
        std::cout<< "ERR - readSimple - Could Not Open Dataset: " << dset << "\n";
        data = NULL;
        H5Gclose(h5grp);
        H5Fclose(h5fil);
        return;
    }

    H5Dread( h5dst , type , H5S_ALL , H5S_ALL , H5P_DEFAULT , data );

    H5Dclose( h5dst );
    H5Gclose( h5grp );
    H5Fclose( h5fil );
}

void vtkDiscoReader::readString(const char *group, const char *dset, 
                                char *buf, int len)
{
    hid_t strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, H5T_VARIABLE);

    hid_t h5fil = H5Fopen( this->FileName , H5F_ACC_RDONLY , H5P_DEFAULT );
    hid_t h5grp = H5Gopen1( h5fil , group );
    hid_t h5dst = H5Dopen1( h5grp , dset );

    hid_t space = H5Dget_space(h5dst);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(space, dims, NULL);

    char **rdata = (char **)malloc(dims[0] * sizeof(char *));
    H5Dread(h5dst, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
    strncpy(buf, rdata[0], len-1);
    buf[len-1] = '\0';

    H5Dvlen_reclaim(strtype, space, H5P_DEFAULT, rdata);
    free(rdata);
    H5Dclose(h5dst);
    H5Sclose(space);
    H5Tclose(strtype);
    H5Gclose(h5grp);
    H5Fclose(h5fil);
}

void vtkDiscoReader::readPatch(const char *group, const char *dset, void *data, 
                                hid_t type, int dim, hsize_t *start, 
                                hsize_t *loc_size, hsize_t *glo_size)
{
    hid_t h5fil = H5Fopen( this->FileName , H5F_ACC_RDONLY , H5P_DEFAULT );
    if(h5fil < 0)
    {
        std::cout<< "ERR - readPatch - Could Not Read File: " << this->FileName << "\n";
        data = NULL;
        return;
    }
    hid_t h5grp = H5Gopen1( h5fil , group );
    if(h5grp < 0)
    {
        std::cout<< "ERR - readPatch - Could Not Open Group: " << group << "\n";
        data = NULL;
        H5Fclose(h5fil);
        return;
    }
    hid_t h5dst = H5Dopen1( h5grp , dset );
    if(h5dst < 0)
    {
        std::cout<< "ERR - readPatch - Could Not Open Dataset: " << dset << "\n";
        data = NULL;
        H5Gclose(h5grp);
        H5Fclose(h5fil);
        return;
    }

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

void vtkDiscoReader::calcXYZ(double x1, double x2, double x3, double *xyz)
{
    if(this->geometry == CYLINDRICAL)
    {
        xyz[0] = x1 * cos(x2);
        xyz[1] = x1 * sin(x2);
        xyz[2] = x3;
    }
    else if(this->geometry == SPHERICAL)
    {
        xyz[0] = x1*sin(x3)*cos(x2);
        xyz[1] = x1*sin(x3)*sin(x2);
        xyz[2] = x1*cos(x3);
    }
    else
    {
        xyz[0] = x1;
        xyz[1] = x2;
        xyz[2] = x3;
    }
}

void vtkDiscoReader::calcVFromContravariant(double x1, double x2, double x3, 
                                            double *V, double *Vxyz)
{
    if(this->geometry == CYLINDRICAL)
    {
        double cp = cos(x2);
        double sp = sin(x2);
        Vxyz[0] = cp*V[0] - x1*sp*V[1];
        Vxyz[1] = sp*V[0] + x1*cp*V[1];
        Vxyz[2] = V[2];
    }
    else if(this->geometry == SPHERICAL)
    {
        double cp = cos(x2);
        double sp = sin(x2);
        double ct = cos(x3);
        double st = sin(x3);
        Vxyz[0] = st*cp*V[0] + ct*cp*x1*V[2] - sp*x1*st*V[1];
        Vxyz[1] = st*sp*V[0] + ct*sp*x1*V[2] + cp*x1*st*V[1];
        Vxyz[2] = ct*V[0] - st*x1*V[2];
    }
    else
    {
        Vxyz[0] = V[0];
        Vxyz[1] = V[1];
        Vxyz[2] = V[2];
    }
}

void vtkDiscoReader::calcVFromOrthonormal(double x1, double x2, double x3, 
                                            double *V, double *Vxyz)
{
    if(this->geometry == CYLINDRICAL)
    {
        double cp = cos(x2);
        double sp = sin(x2);
        Vxyz[0] = cp*V[0] - sp*V[1];
        Vxyz[1] = sp*V[0] + cp*V[1];
        Vxyz[2] = V[2];
    }
    else if(this->geometry == SPHERICAL)
    {
        double cp = cos(x2);
        double sp = sin(x2);
        double ct = cos(x3);
        double st = sin(x3);
        Vxyz[0] = st*cp*V[0] + ct*cp*V[2] - sp*V[1];
        Vxyz[1] = st*sp*V[0] + ct*sp*V[2] + cp*V[1];
        Vxyz[2] = ct*V[0] - st*V[2];
    }
    else
    {
        Vxyz[0] = V[0];
        Vxyz[1] = V[1];
        Vxyz[2] = V[2];
    }
}

void vtkDiscoReader::dimsCreate(int numPieces, int dims[])
{
    std::vector<int> factors;

    int sqrtn = (int) sqrt(numPieces);

    while(numPieces % 2 == 0)
    {
        factors.push_back(2);
        numPieces /= 2;
    }

    int p;
    for(p=3; p<sqrtn && numPieces>1; p += 2)
        while(numPieces % p == 0)
        {
            factors.push_back(p);
            numPieces /= p;
        }
    if(numPieces > 1)
        factors.push_back(numPieces);

    int i;
    dims[0] = 1;
    dims[1] = 1;
    for(i=factors.size()-1; i>=0; i--)
    {
        int f = factors[i];
        if(dims[0]*this->NzTot <= dims[1]*this->NrTot)
            dims[0] *= f;
        else
            dims[1] *= f;
    }
}

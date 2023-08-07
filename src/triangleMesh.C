/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | faSavageHutterFOAM
    \\  /    A nd           | Copyright (C) 2020 Matthias Rauter
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    triangleMesh

Description
    Application to generate a triangle mesh.

Author
    Matthias Rauter matthias@rauter.it

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "stringOps.H"


#include <iostream>
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("triangleMeshDict");
    const word polyMeshDir = polyMesh::meshSubDir;

    IOdictionary demMeshDict
    (
         IOobject
         (
              dictName,
              runTime.system(),
              runTime,
              IOobject::MUST_READ,
              IOobject::NO_WRITE
         )
    );


    if (!demMeshDict.headerOk())
    {
        FatalErrorIn(args.executable())
            << "Cannot open mesh description file\n    "
            << demMeshDict.objectPath()
            << nl
            << exit(FatalError);
    }

    label N;
    vector p1, p2, p3;
    scalar meshH;
    word bottomName;
    word topName;
    word edgesName;

    demMeshDict.lookup("p1") >> p1;
    demMeshDict.lookup("p2") >> p2;
    demMeshDict.lookup("p3") >> p3;
    demMeshDict.lookup("N") >> N;
    demMeshDict.lookup("thickness") >> meshH;
    bottomName = demMeshDict.lookupOrDefault<word>("bottom", "bottom");
    topName = demMeshDict.lookupOrDefault<word>("top", "top");
    edgesName = demMeshDict.lookupOrDefault<word>("edges", "edges");


    Info << "Meshing between " << p1 << ", " << p2 << ", " << p3 << endl << endl;



    label pointsN = 0.5*(N*N+3*N)+1;
    label facesN = N*N;
    label edgesN = 1.5*(N*N+N);


    Info << "Surface mesh has " << endl;
    Info << "   " << pointsN << " points" << endl;
    Info << "   " << facesN << " faces" << endl;
    Info << "   " << edgesN << " edges" << endl << endl;



    pointField rawPoints(pointsN);
    faceList rawFaces(facesN);
    for (label i=0; i<facesN; i++)
        rawFaces[i] = face(3);

    label pointI = 0;
    label faceI = 0;

    for (label j=0; j <= N; j++)
    {
        for (label i=0; i <= N-j; i++)
        {
            rawPoints[pointI++] = p1 + (p2-p1)*scalar(j)/N + (p3-p1)*scalar(i)/N;
            //Info << "point" << pointI << ": " << j << ", " << i << nl;
        }
    }
    pointI = 0;

    for (label j=0; j < N; j++)
    {
        label nextRow = N-j+1;
        for (label i=0; i < N-j; i++)
        {
            rawFaces[faceI][0] = pointI;
            rawFaces[faceI][1] = pointI+1;
            rawFaces[faceI++][2] = pointI+nextRow;
            if(i < N-j-1)
            {
                rawFaces[faceI][0] = pointI+1;
                rawFaces[faceI][1] = pointI+nextRow;
                rawFaces[faceI++][2] = pointI+nextRow+1;
            }
            pointI++;
        }
        pointI++;
    }
    pointI++;


    Info << "Finished generating points. Found" << endl;
    Info << "   " << pointI << " points" << endl;
    Info << "   " << faceI << " faces" << endl << endl;

    /*
    Info << "Faces are" << endl;
    for (label i=0; i < facesN; i++)
    {
        Info << rawFaces[i] << endl;
    }
    */

    Info << "Reconstructing edges..." << endl;

    std::map<std::pair<label, label>, label> parentMap;
    std::map<std::pair<label, label>, label> neighborMap;
    std::map<std::pair<label, label>, label> singleParentMap;


    for(label faceI=0; faceI<facesN; faceI++)
    {
        for (label edgeI=0; edgeI<rawFaces[faceI].size(); edgeI++)
        {
            label p1 = rawFaces[faceI][edgeI];
            label p2 = rawFaces[faceI][(edgeI+1)%rawFaces[faceI].size()];
            if (p1 > p2)
            {
                label p = p1;
                p1 = p2;
                p2 = p;
            }
            std::pair<label, label> edge(p1,p2);
            if (parentMap.count(edge))
            {
                neighborMap[edge]= faceI;
            }
            else if (neighborMap.count(edge))
            {
                Info << "Warning: Edge used by more than two faces!" << endl;
            }
            else
            {
                parentMap[edge]= faceI;
            }
        }
    }

    Info << "Number of owners... " << parentMap.size() << endl;
    Info << "Number of neightbors... " << neighborMap.size() << endl;

    for(auto it = parentMap.begin(); it != parentMap.end(); ++it)
    {
        if (!neighborMap.count(it->first))
        {   
            singleParentMap[it->first] = it->second;
        }
    }

    for(auto it = singleParentMap.begin(); it != singleParentMap.end(); ++it) {
        parentMap.erase(it->first);
    }

    const label edgeCount = parentMap.size()+singleParentMap.size();
    const label borderedgeCount = singleParentMap.size();

    Info << "Number of inner edges... " << edgeCount - borderedgeCount << endl;
    Info << "Number of borders... " << borderedgeCount << endl << endl;

    const label cellCountZ = 1;

    const label cellCount = facesN*cellCountZ;
    const label pointCount = pointsN*(cellCountZ+1);
    const scalar dz = meshH/cellCountZ;

    const label faceCount = cellCount*2+edgeCount;
    //back+front
    //x-paralell-faces
    //y-paralell-faces


    Info<< "Creating Mesh with "
        << cellCount << " cells, "
        << faceCount << " faces and "
        << pointCount << " points" << nl;


    pointField points(pointCount);
    faceList faces(faceCount);
    labelList owner(faceCount);
    labelList neightbor(faceCount);

    label minZboundary, maxZboundary, xyboundary;


    pointField cellCenters(facesN);
    for (label i=0; i<facesN; i++)
    {
        cellCenters[i] = (rawPoints[rawFaces[i][0]]+
                          rawPoints[rawFaces[i][1]]+
                          rawPoints[rawFaces[i][2]])/3.;
    }

    for (label i=0; i<pointsN; i++)
    {
        points[i] = rawPoints[i];
        points[i+pointsN] = rawPoints[i]+vector(0, 0, dz);
    }

    faceI = 0;

    for(auto it = parentMap.begin(); it != parentMap.end(); ++it)
    {
        /* Checking face orientation */
        label p1 = it->first.first;
        label p2 = it->first.second;
        vector faceCenter = (rawPoints[p1]+
                             rawPoints[p2])/2.;

        vector edgeNormal = (rawPoints[p2]-rawPoints[p1]) ^ (vector(0, 0, dz));

        if (((faceCenter-cellCenters[it->second]) & edgeNormal) < 0)
        {
            label p = p1;
            p1 = p2;
            p2 = p;
        }

        /* Internal faces*/
        face &f = faces[faceI];
        f = face(4);
        f[0] = p1;
        f[1] = p2;
        f[2] = p2+pointsN;
        f[3] = p1+pointsN;

        owner[faceI] = it->second;
        neightbor[faceI] = neighborMap[it->first];
        faceI++;
    }
   
    label internalFaces = faceI;

    for(auto it = singleParentMap.begin(); it != singleParentMap.end(); ++it)
    {
        /* Checking face orientation */
        label p1 = it->first.first;
        label p2 = it->first.second;
        vector faceCenter = (rawPoints[p1]+
                             rawPoints[p2])/2.;

        vector edgeNormal = (rawPoints[p2]-rawPoints[p1]) ^ (vector(0, 0, dz));

        if (((faceCenter-cellCenters[it->second]) & edgeNormal) < 0)
        {
            label p = p1;
            p1 = p2;
            p2 = p;
        }

        /* External (x-y-borders) faces*/
        face &f = faces[faceI];
        f = face(4);
        f[0] = p1;
        f[1] = p2;
        f[2] = p2+pointsN;
        f[3] = p1+pointsN;

        owner[faceI] = it->second;
        neightbor[faceI] = -1;
        faceI++;
    }

    label xyFaces = faceI;
    label rawFaceI = 0;

    for (; faceI<facesN+xyFaces; faceI++, rawFaceI++)
    {
        /* Lower faces*/
        face &f = faces[faceI];
        f = face(3);
        f[0] = rawFaces[rawFaceI][0];
        f[1] = rawFaces[rawFaceI][1];
        f[2] = rawFaces[rawFaceI][2];


        /* Upper faces*/
        face &g = faces[faceI+facesN];
        g = face(3);
        g[0] = rawFaces[rawFaceI][2]+pointsN;
        g[1] = rawFaces[rawFaceI][1]+pointsN;
        g[2] = rawFaces[rawFaceI][0]+pointsN;

        vector faceNormal = (rawPoints[f[1]] - rawPoints[f[0]]) ^ (rawPoints[f[2]] - rawPoints[f[0]]);

        if ((faceNormal & vector (0,0,1)) > 0)
        {
            f.flip();
            g.flip();
        }

        owner[faceI] = rawFaceI;
        neightbor[faceI] = -1;
        owner[faceI+facesN] = rawFaceI;
        neightbor[faceI+facesN] = -1;
    }

    word regionName  = polyMesh::defaultRegion;

    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime
        ),
        std::move(points),
        std::move(faces),
        std::move(owner),
        std::move(neightbor)
    );
    xyboundary = internalFaces;
    minZboundary = xyboundary + borderedgeCount;
    maxZboundary = minZboundary + facesN;

    
    label start[3] = {xyboundary, minZboundary, maxZboundary};

    label size[3] = {borderedgeCount, facesN, facesN};

    word name[3] = {edgesName, bottomName, topName};

    List<polyPatch*> patches(3);

    for (int i=0; i<3;i++)
    {
        if (i==0)
        {
                polyPatch *patch = new polyPatch(
                        name[i],
                        size[i],
                        start[i],
                        0,
                        polyBoundaryMesh(
                            IOobject
                            (
                                regionName,
                                runTime.constant(),
                                runTime
                            ),
                            mesh
                            ),
                        name[i]
                        );
                patches[i] = patch;
        }
     else
        {
                emptyPolyPatch *patch = new emptyPolyPatch(
                        name[i],
                        size[i],
                        start[i],
                        0,
                        polyBoundaryMesh(
                            IOobject
                            (
                                regionName,
                                runTime.constant(),
                                runTime
                            ),
                            mesh
                            ),
                        name[i]
                        );
                patches[i] = patch;
        }
    }

    mesh.addPatches(patches);

    Info<< nl << "Writing polyMesh to "  << mesh.instance() << endl;
    mesh.setInstance(mesh.instance());
    mesh.removeFiles();
    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing polyMesh."
            << exit(FatalError);
    }


    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //

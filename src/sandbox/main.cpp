#include "MobiusStrip.h"
#include "MAdLibInterface.h"
#include "P1.h"

using namespace dolfin;
using namespace MAd;

void to_madlibmesh(const Mesh * mesh, 
		   MAd::pMesh MAdMesh)
{
    
    std::map<int,int> MAdToSolverIds;
    std::map<int,int> SolverToMAdIds;
    
    MAdToSolverIds.clear();
    SolverToMAdIds.clear();
    
    uint nVerts = mesh->num_vertices();
    
    const double * xyz = mesh->coordinates();
    
//     for (int iV = 0; iV < nVerts; iV++) {
// 	MAdMesh->add_point(iV + 1, xyz[iV][0], xyz[iV][1], xyz[iV][2]);
// 	std::cout << xyz[iV][0] << "\t" << xyz[iV][1] << "\t" << xyz[iV][2] << std:: endl;
// 	SolverToMAdIds[iV] = iV + 1;
// 	MAdToSolverIds[iV + 1] = iV;
//     }

//   // --- Build the elements ---
//   int dim = solverMesh->getDim();
//   int nElems = solverMesh->nElements();
//   if (dim==3)
//     {
//       const int ** elements = solverMesh->getElements();
//       const int * elemGeoTags = solverMesh->getElemGeoTags();
//       for (int iC=0; iC < nElems; iC++) {
//         pGRegion geom = GM_regionByTag(MAdMesh->model,
//                                        elemGeoTags[iC]);
//         MAdMesh->add_tet(elements[iC][0], elements[iC][1],
//                          elements[iC][2], elements[iC][3],
//                          (pGEntity)geom); 
//       }
//     }
//   else if (dim==2)
//     {
//       const int ** elements = solverMesh->getElements();
//       const int * elemGeoTags = solverMesh->getElemGeoTags();
//       for (int iC=0; iC < nElems; iC++) {
//         pGFace geom = GM_faceByTag(MAdMesh->model,
//                                    elemGeoTags[iC]);
//         MAdMesh->add_triangle(elements[iC][0], elements[iC][1],
//                               elements[iC][2], (pGEntity)geom); 
//       }
//     }

//   /*
//     Here, the entities of the MAdLib mesh sould be classified
//     on their corresponding geometrical entities, like for boundary 
//     faces in 3D for instance. The implementation of this step 
//     is highly dependent on the implementation of Solver_mesh and 
//     Solver_model so it is up to the reader to add the right 
//     instructions here. 

//     Note that the geometrical entities have been created in the 
//     execution of 'exportToMAdModel'. Any mesh entity can be 
//     associated to a geometrical entity using the EN_setWhatIn(...) 
//     function of the MAdLib mesh interface.

//     Note that all the steps involving geometrical entities can be
//     replaced by appropriate constraints on boundary mesh entities
//     (see AdaptInterface.h) but no mesh modification will therefore
//     be applied on the boundaries, which can be problematic for some
//     computations.
//   */

//   MAdMesh->classify_unclassified_entities();
//   MAdMesh->destroyStandAloneEntities();
}

class SizeField : public Expression
{
public:
    void eval(Array<double>& values, const Array<double>& x) const
    {
	values[0] = (1.0 + sin(DOLFIN_PI*5.0*x[0])*sin(DOLFIN_PI*5.0*x[1]))/2.0*0.01 + 0.01;
    }
};

int main(void)
{
//    MobiusStrip mesh(150, 7, 4, 0.5);
//    plot(mesh, "Mobius Strip");

    UnitSquare mesh(3, 3);

    SizeField size;
    P1::FunctionSpace P1(mesh);

    Function size_P1(P1);
    size_P1.interpolate(size);

    plot(size_P1);

    pGModel MAdModel = NULL;
    GM_create(&MAdModel,"theModel");
//    pMesh MAdMesh = M_new(MAdModel);

//    to_madlibmesh(&mesh, MAdMesh);

    for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
    {
	std::vector<double> coordinates(mesh.geometry().dim());
	for (uint i = 0; i < mesh.geometry().dim(); ++i)
	{
	    coordinates[i] = vertex->x()[i];
	    std::cout << coordinates[i] << "\t";
	}
	std::cout << std::endl;	    
    }

    return 0;
}

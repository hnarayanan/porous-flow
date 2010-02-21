#include "MobiusStrip.h"
#include "MAdLibInterface.h"
#include "P1.h"

using namespace dolfin;
using namespace MAd;

void to_madlibmesh(const Mesh* mesh, MAd::pMesh MAdMesh)
{
    
    std::map<int,int> MAdToSolverIds;
    std::map<int,int> SolverToMAdIds;
    
    MAdToSolverIds.clear();
    SolverToMAdIds.clear();
    
    uint current_vertex = 0;

    for (VertexIterator vertex(*mesh); !vertex.end(); ++vertex)
    {
	uint gdim = mesh->geometry().dim();

 	SolverToMAdIds[current_vertex] = current_vertex + 1;
 	MAdToSolverIds[current_vertex + 1] = current_vertex;

	std::vector<double> coordinates(gdim);
	for (uint i = 0; i < gdim; ++i)
	{
	    coordinates[i] = vertex->x()[i];
	}
	if (gdim == 1)
	    MAdMesh->add_point(current_vertex++, coordinates[0], 0.0, 0.0);
	else if (gdim == 2)
	    MAdMesh->add_point(current_vertex++, coordinates[0], coordinates[1], 0.0);
	else if(gdim == 3)
	    MAdMesh->add_point(current_vertex++, coordinates[0], coordinates[1], coordinates[2]);
    }
    

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
    pMesh MAdMesh = M_new(MAdModel);

    to_madlibmesh(&mesh, MAdMesh);
    return 0;
}

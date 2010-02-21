#include "MAdLibInterface.h"
using namespace MAd;

// -------------------------------------------------------------------
// This is an example of a callback function that takes care of a 
// nodal solution when local mesh modifications are applied.
// This function will be registered by 'MAdLibInterface' and will 
// then be called during every local mesh modification.
// --------------------------------------------------------------------
void Solver_CBFunction (pPList before, pPList after, void *data,
                        operationType type, pEntity ppp) {
  
  // Data can point to the object of type 'MAdLibInterface' for instance,
  // depending on what pointer was given when registering the callback function
  // It is not used in this example
  MAdLibInterface * mi = static_cast<MAdLibInterface *>(data);
  
  // The data id used to identify the data attached to mesh entities
  pMeshDataId dataId = MD_lookupMeshDataId("SolutionTag");

  // Do the right manipulation on data according to the mesh modification
  // that is currently applied
  switch (type) {
  case MAd_ESPLIT:
    // Edge split case:
    //   - 'before' contains the split edge (not deleted yet)
    //   - 'after'  contains the two new edges
    //   - 'ppp'    contains the new vertex
    {
      // find the edge to be deleted
      void * temp = NULL;
      pEdge pE = (pEdge) PList_next(before,&temp);
     
      // get coordinates and data at old nodes
      double data0 = 0.;
      pVertex pV0 = E_vertex((pEdge)pE, 0);
      int gotit0 = EN_getDataDbl((pEntity)pV0, dataId,  &data0);
      
      double data1 = 0.;
      pVertex pV1 = E_vertex((pEdge)pE, 1);
      int gotit1 = EN_getDataDbl((pEntity)pV1, dataId,  &data1);
      
      if ( !gotit0 || !gotit1) {
        printf("Error: one of the nodes has no data attached to\n");
        throw;
      }

      // interpolate the data at the new vertex (here linear interpolation)
      double t = E_linearParams(pE,(pVertex)ppp);
      double newData = (1.-t) * data0 + t * data1;
      
      // attach this data to the new vertex
      EN_attachDataDbl(ppp, dataId, newData);
    }
    break;
  case MAd_ECOLLAPSE:
    // Edge collapse case:
    //   - 'before' contains the regions (3D) or faces (2D) of the cavity 
    //                       before the edge collapse (not deleted yet)
    //   - 'after'  contains the regions (3D) or faces (2D) of the cavity 
    //                       after the edge collapse
    //   - 'ppp'    contains the vertex to be deleted (not deleted yet)
    {
      // remove the data on deleted vertex
      EN_deleteData(ppp, dataId);
    }
    break;
  case MAd_FSWAP:
    // Face swap case:
    //   - 'before' contains the regions of the cavity before the face swap (not deleted yet)
    //   - 'after'  contains the regions of the cavity after the face swap
    //   - 'ppp'    contains the swapped face (not deleted yet)
    {
      // nothing to be done for nodal solutions
    }
    break;
  case MAd_ESWAP:
    // Edge swap case:
    //   - 'before' contains the regions (3D) or faces (2D) of the cavity 
    //                       before the edge swap (not deleted yet)
    //   - 'after'  contains the regions (3D) or faces (2D) of the cavity 
    //                       after the edge swap
    //   - 'ppp'    contains the swapped edge (not deleted yet)
    {
      // nothing to be done for nodal solutions
    }
    break;
  default:
    printf("Error: no callback function should be called with this operation: %d",type);
    throw;
  }
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
MAdLibInterface::MAdLibInterface()
{}

//-----------------------------------------------------------------------------
MAdLibInterface::~MAdLibInterface()
{}

//-----------------------------------------------------------------------------
// Main routine for adaptation
void MAdLibInterface::adaptMesh()
{
  //-----------------------------------------------------
  // Step 1: Prepare for adaptation
  //-----------------------------------------------------

  // 1. Delete mesh/solution dependent data in the solver
  solver->deleteData();

  // 2.A. Build the MAdLib geometrical model.
  pGModel MAdModel = NULL;
  GM_create(&MAdModel,"theModel");
  exportToMAdModel(solver->getModel(), MAdModel);

  // 2.B. Build the MAdLib mesh.
  pMesh MAdMesh = M_new(MAdModel);
  exportToMAdMesh(solver->getMesh(), MAdMesh);
  
  // 3. Transfer solution to the MAdLib mesh as an attached data
  attachSolutionToMesh(MAdMesh);
  solver->deallocateSolution();

  // 4. Delete the solver mesh.
  solver->deleteMesh();

  // 5. Build the size field used in adaptation
  PWLSField * sizeField = new PWLSField(MAdMesh);
  buildSizeField(sizeField);

  //-----------------------------------------------------
  // Step 2: Run the adaptation
  //-----------------------------------------------------

  // 6.A. Build the adaptation tool
  MeshAdapter * adapter = new MeshAdapter(MAdMesh,sizeField);
  
  // 6.B. Register the callback function(s) of the solver
  adapter->addCallback(Solver_CBFunction,(void*)this);

  // 6.C. Edit the adaptation parameters if necessary
  adapter->setEdgeLenSqBounds( 1.0/3.0, 3.0 );
  adapter->setNoSwapQuality( 0.1 );
  adapter->setSliverQuality( 0.02 );
  adapter->setSliverPermissionInESplit( true, 10. );
  adapter->setSliverPermissionInECollapse( true, 0.1 );

  // 6.D. Run the adaptation procedure
  adapter->run();

  // 6.E. Optional output
  adapter->printStatistics(std::cout);
  M_writeMsh(MAdMesh,"adapted_mesh.msh",2);

  // 6.F. Clean the adaptation objects
  delete adapter;
  delete sizeField;

  //-----------------------------------------------------
  // Step 3: Rebuild solver data and mesh
  //-----------------------------------------------------

  // 7. Rebuild the solver mesh
  importFromMAdModel(MAdModel, solver->getModel());
  importFromMAdMesh(MAdMesh, solver->getMesh());

  // 8. Get the solution from the MAdLib mesh
  solver->allocateSolution();
  getSolutionFromMesh(MAdMesh);

  // 9. Delete MAdLib mesh
  delete MAdMesh;
  delete MAdModel;

  // 10. Build mesh/solution dependent data in the solver
  solver->allocateAndComputeData();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Converts a MAdLib mesh into a 'Solver_mesh'
void MAdLibInterface::importFromMAdMesh(const MAd::pMesh MAdMesh, 
                                        Solver_mesh * solverMesh)
{
  MAdToSolverIds.clear();
  SolverToMAdIds.clear();

  // --- Mesh dimension ---
  int dim = M_dim(MAdMesh);

  // --- Mesh size ---
  int numVertices = M_numVertices(MAdMesh);
  int numElements;
  if ( dim == 3 ) numElements = M_numRegions(MAdMesh);
  if ( dim == 2 ) numElements = M_numFaces(MAdMesh);

  solverMesh->allocate(numVertices,numElements);

  // --- Build vertices ---
  int solver_Id = 0; // will allow consecutive ids for the solver mesh
  VIter vit = M_vertexIter(MAdMesh);
  while (pVertex pv = VIter_next(vit))
    {
      // get MAdLib mesh id
      int MAd_Id = EN_id((pEntity)pv);

      // get coordinates
      double xyz[3];
      V_coord(pv,xyz);
      
      // add node in solver mesh
      solverMesh->addNode(solver_Id,xyz[0],xyz[1],xyz[2]);
      
      // fill in id's tables
      MAdToSolverIds[MAd_Id] = solver_Id;
      SolverToMAdIds[solver_Id] = MAd_Id;

      solver_Id++;
    }
  VIter_delete(vit);

  // --- Build elements ---
  int solver_elem_Id = 0;
  if (dim==3) {
    RIter rit = M_regionIter(MAdMesh);
    while (pRegion pr = RIter_next(rit))
      {
        // get list of node id's in the solver mesh
        int nodes[4];
        pPList rVerts = R_vertices(pr);
        void * temp = NULL;
        int iN = 0;
        while ( pVertex pv = (pVertex)PList_next(rVerts,&temp) )
          {
            int MAd_Id = EN_id((pEntity)pv);
            nodes[iN++] = MAdToSolverIds[MAd_Id];
          }
        PList_delete(rVerts);

        // add the element to the solver mesh
        solverMesh->addElement(solver_elem_Id, nodes);
        solver_elem_Id++;
      }
    RIter_delete(rit);
  }
  else  if (dim==2) {
    FIter fit = M_faceIter(MAdMesh);
    while (pFace pf = FIter_next(fit))
      {
        // get list of node id's in the solver mesh
        int nodes[3];
        pPList fVerts = F_vertices(pf,1);
        void * temp = NULL;
        int iN = 0;
        while ( pVertex pv = (pVertex)PList_next(fVerts,&temp) )
          {
            int MAd_Id = EN_id((pEntity)pv);
            nodes[iN++] = MAdToSolverIds[MAd_Id];
          }
        PList_delete(fVerts);

        // add the element to the solver mesh
        solverMesh->addElement(solver_elem_Id, nodes);
        solver_elem_Id++;
      }
    FIter_delete(fit);
  }
}

//-----------------------------------------------------------------------------
// Converts a 'Solver_mesh' into a MAdLib mesh
void MAdLibInterface::exportToMAdMesh(const Solver_mesh * solverMesh, 
                                      MAd::pMesh MAdMesh)
{
  // --- Build the vertices ---
  MAdToSolverIds.clear();
  SolverToMAdIds.clear();
  int nVerts = solverMesh->nVertices();
  const double ** xyz = solverMesh->getCoordinates();
  for (int iV=0; iV < nVerts; iV++) {
    MAdMesh->add_point(iV+1,xyz[iV][0],xyz[iV][1],xyz[iV][2]);
    SolverToMAdIds[iV] = iV+1;
    MAdToSolverIds[iV+1] = iV;
  }

  // --- Build the elements ---
  int dim = solverMesh->getDim();
  int nElems = solverMesh->nElements();
  if (dim==3)
    {
      const int ** elements = solverMesh->getElements();
      const int * elemGeoTags = solverMesh->getElemGeoTags();
      for (int iC=0; iC < nElems; iC++) {
        pGRegion geom = GM_regionByTag(MAdMesh->model,
                                       elemGeoTags[iC]);
        MAdMesh->add_tet(elements[iC][0], elements[iC][1],
                         elements[iC][2], elements[iC][3],
                         (pGEntity)geom); 
      }
    }
  else if (dim==2)
    {
      const int ** elements = solverMesh->getElements();
      const int * elemGeoTags = solverMesh->getElemGeoTags();
      for (int iC=0; iC < nElems; iC++) {
        pGFace geom = GM_faceByTag(MAdMesh->model,
                                   elemGeoTags[iC]);
        MAdMesh->add_triangle(elements[iC][0], elements[iC][1],
                              elements[iC][2], (pGEntity)geom); 
      }
    }

  /*
    Here, the entities of the MAdLib mesh sould be classified
    on their corresponding geometrical entities, like for boundary 
    faces in 3D for instance. The implementation of this step 
    is highly dependent on the implementation of Solver_mesh and 
    Solver_model so it is up to the reader to add the right 
    instructions here. 

    Note that the geometrical entities have been created in the 
    execution of 'exportToMAdModel'. Any mesh entity can be 
    associated to a geometrical entity using the EN_setWhatIn(...) 
    function of the MAdLib mesh interface.

    Note that all the steps involving geometrical entities can be
    replaced by appropriate constraints on boundary mesh entities
    (see AdaptInterface.h) but no mesh modification will therefore
    be applied on the boundaries, which can be problematic for some
    computations.
  */

  MAdMesh->classify_unclassified_entities();
  MAdMesh->destroyStandAloneEntities();
}

//-----------------------------------------------------------------------------
// Create in MAdModel all geometrical entities listed in solverModel.
void MAdLibInterface::exportToMAdModel(const Solver_model * solverModel, 
                                       MAd::pGModel MAdModel)
{
  std::set<std::pair<int,int> > geometry = solverModel->getAllGeoEntities();
  std::set<std::pair<int,int> >::const_iterator geoIter = geometry.begin();
  for (; geoIter != geometry.end(); geoIter++) {
    int dim = (*geoIter).first;
    int id  = (*geoIter).second;
    GM_entityByTag(MAdModel,dim,id);
  }
}

//-----------------------------------------------------------------------------
// Build a field of prescribed edges lengths on the domain.
void MAdLibInterface::buildSizeField(MAd::PWLSField * sizeField)
{
  // First option: keep actual edges lengths
  sizeField->setCurrentSize();

  // Second option: compute it from solver functions
  VIter vit = M_vertexIter(sizeField->getMesh());
  while (pVertex pv = VIter_next(vit))
    {
      // get solver point id
      int MAd_Id = EN_id((pEntity)pv);
      int solver_Id = MAdToSolverIds[MAd_Id];
      
      // get the edge length prescribed by the solver
      double length = solver->prescribedEdgeLength(solver_Id);

      // fill in the size field
      sizeField->setSize((pEntity)pv, length);
    }
  VIter_delete(vit);
}

//-----------------------------------------------------------------------------
void MAdLibInterface::attachSolutionToMesh(MAd::pMesh MAdMesh)
{
  // Get the solution database. Here we assume that it is a nodal solution.
  const Solver_solution * solution = solver->getSolution();

  // The data id used to identify the data attached to mesh entities
  pMeshDataId dataId = MD_lookupMeshDataId("SolutionTag");

  VIter vit = M_vertexIter(MAdMesh);
  while (pVertex pv = VIter_next(vit))
    {
      // get solver point id
      int MAd_Id = EN_id((pEntity)pv);
      int solver_Id = MAdToSolverIds[MAd_Id];
      
      double data = (*solution)[solver_Id];
      
      // attach data to the mesh vertex
      EN_attachDataDbl((pEntity)pv,dataId,data);
    }
  VIter_delete(vit);
}

//-----------------------------------------------------------------------------
void MAdLibInterface::getSolutionFromMesh(MAd::pMesh MAdMesh)
{
  // Get the solution database. Here we assume that it is a nodal solution.
  Solver_solution * solution = solver->getSolution();

  // The data id used to identify the data attached to mesh entities
  pMeshDataId dataId = MD_lookupMeshDataId("SolutionTag");

  VIter vit = M_vertexIter(MAdMesh);
  while (pVertex pv = VIter_next(vit))
    {
      // get solver point id
      pPoint pp = V_point(pv);
      int MAdId = P_id(pp);
      int solver_Id = MAdToSolverIds[MAdId];
      
      // get attached data and delete it
      double data;
      EN_getDataDbl((pEntity)pv,dataId,&data);
      EN_deleteData((pEntity)pv,dataId);
      
      *(*solution)[solver_Id] = data;
    }
  VIter_delete(vit);
}

//-----------------------------------------------------------------------------
void MAdLibInterface::importFromMAdModel(const MAd::pGModel MAdModel, Solver_model * solverModel)
{
    // Do nothing, yet.
}

std::set<std::pair<int,int> > Solver_model::getAllGeoEntities() const
{
    // Do nothing, yet.
}

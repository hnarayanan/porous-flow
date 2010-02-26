// -*- C++ -*-
// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
// -------------------------------------------------------------------
// Author: Gaetan Compere
//
// This file provides an example of an interface to MAdLib as it 
// could be implemented in a physical solver requiring mesh adaptivity.
// -------------------------------------------------------------------

#ifndef __MADLIBINTERFACE_H
#define __MADLIBINTERFACE_H

#include "MAdLib/ModelInterface.h"
#include "MAdLib/MeshDataBaseInterface.h"
#include "MAdLib/AdaptInterface.h"
#include "MAdLib/PWLinearSField.h"
#include <utility>
#include <set>

//-----------------------------------------------------------------------------
// For this example, what is needed in the solver side ?
//-----------------------------------------------------------------------------

/*
Solver class containing the solver geometrical model if any.
  ( Note that all the steps involving geometrical entities can be
    replaced by appropriate constraints on boundary mesh entities
    (see AdaptInterface.h) but no mesh modification will therefore
    be applied on the boundaries, which can be problematic for some
    computations. )
*/
class Solver_model
{
public:
  void addGeoEntity(int dim, int id);
  std::set<std::pair<int,int> > getAllGeoEntities() const;
  // return a set of pairs(dimension,id) 
  // each pair representing a geometric entity.
};

/*
  Solver class containing the solver mesh
*/
class Solver_mesh
{
public:
  void allocate (int nNodes, int nElements) {}
  void addNode (int id, double x, double y, double z) {}
  void addElement (int id, int * nodes) {}
  int getDim() const {return -1;}
  int nVertices() const {return -1;}
  int nElements() const {return -1;}
  const double ** getCoordinates() const {return NULL;}
  const int ** getElements() const {return NULL;}
  const int * getElemGeoTags() const {return NULL;}
};

/*
  Solver solution. We assume a nodal solution but the current example can be
  easily extended to other discretizations.
*/
class Solver_solution
{
public:
  double * operator[](int i) {return NULL;}
  const double operator[](int i) const {return 0.;}
};

/*
  Solver class containing pointers to solver data, solution and mesh
*/
class Solver
{
public:
  Solver_model    * getModel()    {return model;}
  Solver_mesh     * getMesh()     {return mesh;}
  Solver_solution * getSolution() {return solution;}
  void deleteMesh() {}
  void deallocateSolution() {}
  void allocateSolution() {}
  // optional functions:
  void deleteData() {}
  void allocateAndComputeData() {}
  double prescribedEdgeLength(int node) {return 0.;}
private:
  Solver_model    * model;
  Solver_mesh     * mesh;
  Solver_solution * solution;
};

//-----------------------------------------------------------------------------
// Class interfacing MAdLib with 'Solver'
//-----------------------------------------------------------------------------
class MAdLibInterface {

public:
   
  MAdLibInterface();
  ~MAdLibInterface();

  void adaptMesh();

private:

  // Mesh to mesh conversion
  void importFromMAdMesh(const MAd::pMesh, Solver_mesh *);
  void exportToMAdMesh(const Solver_mesh *, MAd::pMesh);
  void importFromMAdModel(const MAd::pGModel, Solver_model *);
  void exportToMAdModel(const Solver_model *, MAd::pGModel);

  // Size field construction
  void buildSizeField(MAd::PWLSField *);

  // Solution to solution conversion
  void attachSolutionToMesh(MAd::pMesh);
  void getSolutionFromMesh(MAd::pMesh);

private:

  // The solver that needs mesh adaptivity
  Solver * solver;

  // Correspondancy tables between nodal id's in the solver
  // and in the MAdLib mesh
  std::map<int,int> MAdToSolverIds;
  std::map<int,int> SolverToMAdIds;
};

//-----------------------------------------------------------------------------

#endif


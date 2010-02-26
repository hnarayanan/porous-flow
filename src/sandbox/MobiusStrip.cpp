#include "MobiusStrip.h"

using namespace dolfin;

MobiusStrip::MobiusStrip(uint nl, uint nw, uint num_twists, double width) : Mesh()
{
    // Receive mesh according to parallel policy
    if (MPI::is_receiver()) { MeshPartitioning::partition(*this); return; }
    
    if (nl < 1 || nw < 1)
	error("Size of Mobius strip must be at least 1 in each dimension.");
    
    if (num_twists % 2 != 0)
	error("Only strips with an even numbers of twists allowed!");
    
    rename("mesh", "Mesh of the Mobius strip");
    
    // Open mesh for editing
    MeshEditor editor;
    editor.open(*this, CellType::triangle, 2, 3);
    
    // Create empty lists to store vertices and cells
    editor.init_vertices(nl*nw);
    editor.init_cells(nl*(nw - 1)*2);

    // Populate the list of vertices:
    uint vertex = 0;
    for (uint i = 0; i < nl; i++)
    {
	const double u = static_cast<double>(i)/static_cast<double>(nl)*2*DOLFIN_PI;
	for (uint j = 0; j < nw; j++)
	{
	    const double v = static_cast<double>(j)/(static_cast<double>(nw) - 1.0)*width;
	    editor.add_vertex(vertex++,
			      cos(u) + v*cos(num_twists*u/2.0)*cos(u),
			      sin(u) + v*cos(num_twists*u/2.0)*sin(u),
                              v*sin(num_twists*u/2.0));
	}
    }

    // Populate list of cells
    uint cell = 0;
    for (uint i = 0; i < nl - 1; i++)
    {
	for (uint j = 0; j < nw - 1; j++)
	{
	    editor.add_cell(cell++, i*nw + j, (i + 1)*nw + j + 1, i*nw + j + 1);
	    editor.add_cell(cell++, i*nw + j, (i + 1)*nw + j    , (i + 1)*nw + j + 1);
	}
    }

    for (uint j = 0; j < nw - 1; j++)
    {
	editor.add_cell(cell++, (nl - 1)*nw + j, j + 1, (nl - 1)*nw + j + 1);
	editor.add_cell(cell++, (nl - 1)*nw + j, j    , j + 1);
    }

    // Close mesh editor
    editor.close();

    // Broadcast mesh according to parallel policy
    if (MPI::is_broadcaster()) { MeshPartitioning::partition(*this); return; }
}

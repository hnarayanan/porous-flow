#ifndef __MOBIUS_STRIP_H
#define __MOBIUS_STRIP_H

#include <dolfin.h>
#include <math.h>

namespace dolfin
{
    class MobiusStrip : public Mesh
    {
    public:
	MobiusStrip(uint nl, uint nw, uint num_twists, double width);
    };
}

#endif

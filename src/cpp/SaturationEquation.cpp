#include "SaturationEquation.h"

/// Constructor
saturationequation_0_finite_element_0::saturationequation_0_finite_element_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
saturationequation_0_finite_element_0::~saturationequation_0_finite_element_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* saturationequation_0_finite_element_0::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return the cell shape
ufc::shape saturationequation_0_finite_element_0::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int saturationequation_0_finite_element_0::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int saturationequation_0_finite_element_0::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int saturationequation_0_finite_element_0::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void saturationequation_0_finite_element_0::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void saturationequation_0_finite_element_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void saturationequation_0_finite_element_0::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
    
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 1)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[1][1] = \
    {{0}};
    
    static const double dmats1[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void saturationequation_0_finite_element_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double saturationequation_0_finite_element_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[1][1][2] = {{{0.333333333333333, 0.333333333333333}}};
    static const double W[1][1] = {{1}};
    static const double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void saturationequation_0_finite_element_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void saturationequation_0_finite_element_0::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int saturationequation_0_finite_element_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* saturationequation_0_finite_element_0::create_sub_element(unsigned int i) const
{
    return new saturationequation_0_finite_element_0();
}


/// Constructor
saturationequation_0_finite_element_1::saturationequation_0_finite_element_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
saturationequation_0_finite_element_1::~saturationequation_0_finite_element_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* saturationequation_0_finite_element_1::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return the cell shape
ufc::shape saturationequation_0_finite_element_1::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int saturationequation_0_finite_element_1::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int saturationequation_0_finite_element_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int saturationequation_0_finite_element_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void saturationequation_0_finite_element_1::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void saturationequation_0_finite_element_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void saturationequation_0_finite_element_1::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
    
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 1)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[1][1] = \
    {{0}};
    
    static const double dmats1[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void saturationequation_0_finite_element_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double saturationequation_0_finite_element_1::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[1][1][2] = {{{0.333333333333333, 0.333333333333333}}};
    static const double W[1][1] = {{1}};
    static const double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void saturationequation_0_finite_element_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void saturationequation_0_finite_element_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int saturationequation_0_finite_element_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* saturationequation_0_finite_element_1::create_sub_element(unsigned int i) const
{
    return new saturationequation_0_finite_element_1();
}

/// Constructor
saturationequation_0_dof_map_0::saturationequation_0_dof_map_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
saturationequation_0_dof_map_0::~saturationequation_0_dof_map_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* saturationequation_0_dof_map_0::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool saturationequation_0_dof_map_0::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool saturationequation_0_dof_map_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void saturationequation_0_dof_map_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void saturationequation_0_dof_map_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int saturationequation_0_dof_map_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int saturationequation_0_dof_map_0::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int saturationequation_0_dof_map_0::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int saturationequation_0_dof_map_0::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int saturationequation_0_dof_map_0::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int saturationequation_0_dof_map_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void saturationequation_0_dof_map_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[2][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void saturationequation_0_dof_map_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void saturationequation_0_dof_map_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void saturationequation_0_dof_map_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.333333333333333*x[0][0] + 0.333333333333333*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[0][1] = 0.333333333333333*x[0][1] + 0.333333333333333*x[1][1] + 0.333333333333333*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int saturationequation_0_dof_map_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* saturationequation_0_dof_map_0::create_sub_dof_map(unsigned int i) const
{
    return new saturationequation_0_dof_map_0();
}


/// Constructor
saturationequation_0_dof_map_1::saturationequation_0_dof_map_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
saturationequation_0_dof_map_1::~saturationequation_0_dof_map_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* saturationequation_0_dof_map_1::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool saturationequation_0_dof_map_1::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool saturationequation_0_dof_map_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void saturationequation_0_dof_map_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void saturationequation_0_dof_map_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int saturationequation_0_dof_map_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int saturationequation_0_dof_map_1::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int saturationequation_0_dof_map_1::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int saturationequation_0_dof_map_1::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int saturationequation_0_dof_map_1::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int saturationequation_0_dof_map_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void saturationequation_0_dof_map_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[2][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void saturationequation_0_dof_map_1::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void saturationequation_0_dof_map_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void saturationequation_0_dof_map_1::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.333333333333333*x[0][0] + 0.333333333333333*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[0][1] = 0.333333333333333*x[0][1] + 0.333333333333333*x[1][1] + 0.333333333333333*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int saturationequation_0_dof_map_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* saturationequation_0_dof_map_1::create_sub_dof_map(unsigned int i) const
{
    return new saturationequation_0_dof_map_1();
}


/// Constructor
saturationequation_0_cell_integral_0_quadrature::saturationequation_0_cell_integral_0_quadrature() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
saturationequation_0_cell_integral_0_quadrature::~saturationequation_0_cell_integral_0_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void saturationequation_0_cell_integral_0_quadrature::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Set scale factor
    const double det = std::abs(detJ);
    
    
    // Array of quadrature weights
    static const double W1 = 0.5;
    // Quadrature points on the UFC reference element: (0.333333333333333, 0.333333333333333)
    
    // Value of basis functions at quadrature points.
    static const double FE0[1][1] = \
    {{1}};
    
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', False), ('ignore zero tables', False), ('non zero columns', False), ('remove zero terms', False), ('ignore ones', False)
    // Total number of operations to compute element tensor: 4
    
    // Loop quadrature points for integral
    // Number of operations to compute element tensor for following IP loop = 4
    // Only 1 integration point, omitting IP loop.
    
    // Number of operations for primary indices: 4
    for (unsigned int j = 0; j < 1; j++)
    {
      for (unsigned int k = 0; k < 1; k++)
      {
        // Number of operations to compute entry: 4
        A[j*1 + k] += FE0[0][j]*FE0[0][k]*W1*det;
      }// end loop over 'k'
    }// end loop over 'j'
}

/// Constructor
saturationequation_0_cell_integral_0::saturationequation_0_cell_integral_0() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
saturationequation_0_cell_integral_0::~saturationequation_0_cell_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void saturationequation_0_cell_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c) const
{
    // Reset values of the element tensor block
    A[0] = 0;
    
    // Add all contributions to element tensor
    integral_0_quadrature.tabulate_tensor(A, w, c);
}

/// Constructor
saturationequation_form_0::saturationequation_form_0() : ufc::form()
{
    // Do nothing
}

/// Destructor
saturationequation_form_0::~saturationequation_form_0()
{
    // Do nothing
}

/// Return a string identifying the form
const char* saturationequation_form_0::signature() const
{
    return "Form([Integral(Product(BasisFunction(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 0), BasisFunction(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 1)), Measure('cell', 0, None))])";
}

/// Return the rank of the global tensor (r)
unsigned int saturationequation_form_0::rank() const
{
    return 2;
}

/// Return the number of coefficients (n)
unsigned int saturationequation_form_0::num_coefficients() const
{
    return 0;
}

/// Return the number of cell integrals
unsigned int saturationequation_form_0::num_cell_integrals() const
{
    return 1;
}

/// Return the number of exterior facet integrals
unsigned int saturationequation_form_0::num_exterior_facet_integrals() const
{
    return 0;
}

/// Return the number of interior facet integrals
unsigned int saturationequation_form_0::num_interior_facet_integrals() const
{
    return 0;
}

/// Create a new finite element for argument function i
ufc::finite_element* saturationequation_form_0::create_finite_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new saturationequation_0_finite_element_0();
      break;
    case 1:
      return new saturationequation_0_finite_element_1();
      break;
    }
    return 0;
}

/// Create a new dof map for argument function i
ufc::dof_map* saturationequation_form_0::create_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new saturationequation_0_dof_map_0();
      break;
    case 1:
      return new saturationequation_0_dof_map_1();
      break;
    }
    return 0;
}

/// Create a new cell integral on sub domain i
ufc::cell_integral* saturationequation_form_0::create_cell_integral(unsigned int i) const
{
    return new saturationequation_0_cell_integral_0();
}

/// Create a new exterior facet integral on sub domain i
ufc::exterior_facet_integral* saturationequation_form_0::create_exterior_facet_integral(unsigned int i) const
{
    return 0;
}

/// Create a new interior facet integral on sub domain i
ufc::interior_facet_integral* saturationequation_form_0::create_interior_facet_integral(unsigned int i) const
{
    return 0;
}


/// Constructor
saturationequation_1_finite_element_0::saturationequation_1_finite_element_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
saturationequation_1_finite_element_0::~saturationequation_1_finite_element_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* saturationequation_1_finite_element_0::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return the cell shape
ufc::shape saturationequation_1_finite_element_0::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int saturationequation_1_finite_element_0::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int saturationequation_1_finite_element_0::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int saturationequation_1_finite_element_0::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void saturationequation_1_finite_element_0::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void saturationequation_1_finite_element_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void saturationequation_1_finite_element_0::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
    
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 1)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[1][1] = \
    {{0}};
    
    static const double dmats1[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void saturationequation_1_finite_element_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double saturationequation_1_finite_element_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[1][1][2] = {{{0.333333333333333, 0.333333333333333}}};
    static const double W[1][1] = {{1}};
    static const double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void saturationequation_1_finite_element_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void saturationequation_1_finite_element_0::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int saturationequation_1_finite_element_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* saturationequation_1_finite_element_0::create_sub_element(unsigned int i) const
{
    return new saturationequation_1_finite_element_0();
}


/// Constructor
saturationequation_1_finite_element_1_0::saturationequation_1_finite_element_1_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
saturationequation_1_finite_element_1_0::~saturationequation_1_finite_element_1_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* saturationequation_1_finite_element_1_0::signature() const
{
    return "FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1)";
}

/// Return the cell shape
ufc::shape saturationequation_1_finite_element_1_0::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int saturationequation_1_finite_element_1_0::space_dimension() const
{
    return 6;
}

/// Return the rank of the value space
unsigned int saturationequation_1_finite_element_1_0::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int saturationequation_1_finite_element_1_0::value_dimension(unsigned int i) const
{
    return 2;
}

/// Evaluate basis function i at given point in cell
void saturationequation_1_finite_element_1_0::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Reset values
    values[0] = 0;
    values[1] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    
    // Table(s) of coefficients
    static const double coefficients0[6][3] = \
    {{0.942809041582063, 0.577350269189626, -0.333333333333333},
    {-0.471404520791032, -0.288675134594813, 0.166666666666667},
    {0.471404520791032, -0.577350269189626, -0.666666666666667},
    {0.471404520791032, 0.288675134594813, 0.833333333333333},
    {-0.471404520791032, -0.288675134594813, 0.166666666666667},
    {0.942809041582063, 0.577350269189626, -0.333333333333333}};
    
    static const double coefficients1[6][3] = \
    {{-0.471404520791032, 0, -0.333333333333333},
    {0.942809041582063, 0, 0.666666666666667},
    {0.471404520791032, 0, 0.333333333333333},
    {-0.942809041582063, 0, -0.666666666666667},
    {-0.471404520791032, 0.866025403784439, 0.166666666666667},
    {-0.471404520791032, -0.866025403784439, 0.166666666666667}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff1_0 = coefficients1[dof][0];
    const double coeff1_1 = coefficients1[dof][1];
    const double coeff1_2 = coefficients1[dof][2];
    
    // Compute value(s)
    const double tmp0_0 = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2;
    const double tmp0_1 = coeff1_0*basisvalue0 + coeff1_1*basisvalue1 + coeff1_2*basisvalue2;
    // Using contravariant Piola transform to map values back to the physical element
    values[0] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
    values[1] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
}

/// Evaluate all basis functions at given point in cell
void saturationequation_1_finite_element_1_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void saturationequation_1_finite_element_1_0::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
    
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 1)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 2*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    
    // Table(s) of coefficients
    static const double coefficients0[6][3] = \
    {{0.942809041582063, 0.577350269189626, -0.333333333333333},
    {-0.471404520791032, -0.288675134594813, 0.166666666666667},
    {0.471404520791032, -0.577350269189626, -0.666666666666667},
    {0.471404520791032, 0.288675134594813, 0.833333333333333},
    {-0.471404520791032, -0.288675134594813, 0.166666666666667},
    {0.942809041582063, 0.577350269189626, -0.333333333333333}};
    
    static const double coefficients1[6][3] = \
    {{-0.471404520791032, 0, -0.333333333333333},
    {0.942809041582063, 0, 0.666666666666667},
    {0.471404520791032, 0, 0.333333333333333},
    {-0.942809041582063, 0, -0.666666666666667},
    {-0.471404520791032, 0.866025403784439, 0.166666666666667},
    {-0.471404520791032, -0.866025403784439, 0.166666666666667}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[3][3] = \
    {{0, 0, 0},
    {4.89897948556636, 0, 0},
    {0, 0, 0}};
    
    static const double dmats1[3][3] = \
    {{0, 0, 0},
    {2.44948974278318, 0, 0},
    {4.24264068711928, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [2*num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    double coeff1_0 = 0;
    double coeff1_1 = 0;
    double coeff1_2 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff1_0 = 0;
    double new_coeff1_1 = 0;
    double new_coeff1_2 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
      new_coeff1_0 = coefficients1[dof][0];
      new_coeff1_1 = coefficients1[dof][1];
      new_coeff1_2 = coefficients1[dof][2];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
        coeff1_0 = new_coeff1_0;
        coeff1_1 = new_coeff1_1;
        coeff1_2 = new_coeff1_2;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2];
          new_coeff1_0 = coeff1_0*dmats0[0][0] + coeff1_1*dmats0[1][0] + coeff1_2*dmats0[2][0];
          new_coeff1_1 = coeff1_0*dmats0[0][1] + coeff1_1*dmats0[1][1] + coeff1_2*dmats0[2][1];
          new_coeff1_2 = coeff1_0*dmats0[0][2] + coeff1_1*dmats0[1][2] + coeff1_2*dmats0[2][2];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2];
          new_coeff1_0 = coeff1_0*dmats1[0][0] + coeff1_1*dmats1[1][0] + coeff1_2*dmats1[2][0];
          new_coeff1_1 = coeff1_0*dmats1[0][1] + coeff1_1*dmats1[1][1] + coeff1_2*dmats1[2][1];
          new_coeff1_2 = coeff1_0*dmats1[0][2] + coeff1_1*dmats1[1][2] + coeff1_2*dmats1[2][2];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      // Correct values by the contravariant Piola transform
      const double tmp0_0 = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2;
      const double tmp0_1 = new_coeff1_0*basisvalue0 + new_coeff1_1*basisvalue1 + new_coeff1_2*basisvalue2;
      derivatives[deriv_num] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
      derivatives[num_derivatives + deriv_num] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
        values[num_derivatives + row] += transform[row][col]*derivatives[num_derivatives + col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void saturationequation_1_finite_element_1_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double saturationequation_1_finite_element_1_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[6][1][2] = {{{0.666666666666667, 0.333333333333333}}, {{0.333333333333333, 0.666666666666667}}, {{0, 0.333333333333333}}, {{0, 0.666666666666667}}, {{0.333333333333333, 0}}, {{0.666666666666667, 0}}};
    static const double W[6][1] = {{1}, {1}, {1}, {1}, {1}, {1}};
    static const double D[6][1][2] = {{{1, 1}}, {{1, 1}}, {{1, 0}}, {{1, 0}}, {{0, -1}}, {{0, -1}}};
    
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    const double Jinv_00 =  J_11 / detJ;
    const double Jinv_01 = -J_01 / detJ;
    const double Jinv_10 = -J_10 / detJ;
    const double Jinv_11 =  J_00 / detJ;
    
    double copyofvalues[2];
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[2];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Copy old values:
    copyofvalues[0] = values[0];
    copyofvalues[1] = values[1];
    // Do the inverse of div piola 
    values[0] = detJ*(Jinv_00*copyofvalues[0]+Jinv_01*copyofvalues[1]);
    values[1] = detJ*(Jinv_10*copyofvalues[0]+Jinv_11*copyofvalues[1]);
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 2; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void saturationequation_1_finite_element_1_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void saturationequation_1_finite_element_1_0::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    // Evaluate at vertices and use Piola mapping
    vertex_values[0] = (1.0/detJ)*(dof_values[2]*2*J_00 + dof_values[3]*J_00 + dof_values[4]*(-2*J_01) + dof_values[5]*J_01);
    vertex_values[2] = (1.0/detJ)*(dof_values[0]*2*J_00 + dof_values[1]*J_00 + dof_values[4]*(J_00 + J_01) + dof_values[5]*(2*J_00 - 2*J_01));
    vertex_values[4] = (1.0/detJ)*(dof_values[0]*J_01 + dof_values[1]*2*J_01 + dof_values[2]*(J_00 + J_01) + dof_values[3]*(2*J_00 - 2*J_01));
    vertex_values[1] = (1.0/detJ)*(dof_values[2]*2*J_10 + dof_values[3]*J_10 + dof_values[4]*(-2*J_11) + dof_values[5]*J_11);
    vertex_values[3] = (1.0/detJ)*(dof_values[0]*2*J_10 + dof_values[1]*J_10 + dof_values[4]*(J_10 + J_11) + dof_values[5]*(2*J_10 - 2*J_11));
    vertex_values[5] = (1.0/detJ)*(dof_values[0]*J_11 + dof_values[1]*2*J_11 + dof_values[2]*(J_10 + J_11) + dof_values[3]*(2*J_10 - 2*J_11));
}

/// Return the number of sub elements (for a mixed element)
unsigned int saturationequation_1_finite_element_1_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* saturationequation_1_finite_element_1_0::create_sub_element(unsigned int i) const
{
    return new saturationequation_1_finite_element_1_0();
}


/// Constructor
saturationequation_1_finite_element_1_1::saturationequation_1_finite_element_1_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
saturationequation_1_finite_element_1_1::~saturationequation_1_finite_element_1_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* saturationequation_1_finite_element_1_1::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return the cell shape
ufc::shape saturationequation_1_finite_element_1_1::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int saturationequation_1_finite_element_1_1::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int saturationequation_1_finite_element_1_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int saturationequation_1_finite_element_1_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void saturationequation_1_finite_element_1_1::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void saturationequation_1_finite_element_1_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void saturationequation_1_finite_element_1_1::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
    
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 1)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[1][1] = \
    {{0}};
    
    static const double dmats1[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void saturationequation_1_finite_element_1_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double saturationequation_1_finite_element_1_1::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[1][1][2] = {{{0.333333333333333, 0.333333333333333}}};
    static const double W[1][1] = {{1}};
    static const double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void saturationequation_1_finite_element_1_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void saturationequation_1_finite_element_1_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int saturationequation_1_finite_element_1_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* saturationequation_1_finite_element_1_1::create_sub_element(unsigned int i) const
{
    return new saturationequation_1_finite_element_1_1();
}


/// Constructor
saturationequation_1_finite_element_1::saturationequation_1_finite_element_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
saturationequation_1_finite_element_1::~saturationequation_1_finite_element_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* saturationequation_1_finite_element_1::signature() const
{
    return "MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) })";
}

/// Return the cell shape
ufc::shape saturationequation_1_finite_element_1::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int saturationequation_1_finite_element_1::space_dimension() const
{
    return 7;
}

/// Return the rank of the value space
unsigned int saturationequation_1_finite_element_1::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int saturationequation_1_finite_element_1::value_dimension(unsigned int i) const
{
    return 3;
}

/// Evaluate basis function i at given point in cell
void saturationequation_1_finite_element_1::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Reset values
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    
    if (0 <= i && i <= 5)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_1_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
      const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
      const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    
      // Table(s) of coefficients
      static const double coefficients0[6][3] =   \
      {{0.942809041582063, 0.577350269189626, -0.333333333333333},
      {-0.471404520791032, -0.288675134594813, 0.166666666666667},
      {0.471404520791032, -0.577350269189626, -0.666666666666667},
      {0.471404520791032, 0.288675134594813, 0.833333333333333},
      {-0.471404520791032, -0.288675134594813, 0.166666666666667},
      {0.942809041582063, 0.577350269189626, -0.333333333333333}};
    
      static const double coefficients1[6][3] =   \
      {{-0.471404520791032, 0, -0.333333333333333},
      {0.942809041582063, 0, 0.666666666666667},
      {0.471404520791032, 0, 0.333333333333333},
      {-0.942809041582063, 0, -0.666666666666667},
      {-0.471404520791032, 0.866025403784439, 0.166666666666667},
      {-0.471404520791032, -0.866025403784439, 0.166666666666667}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff1_0 =   coefficients1[dof][0];
      const double coeff1_1 =   coefficients1[dof][1];
      const double coeff1_2 =   coefficients1[dof][2];
    
      // Compute value(s)
      const double tmp0_0 = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2;
      const double tmp0_1 = coeff1_0*basisvalue0 + coeff1_1*basisvalue1 + coeff1_2*basisvalue2;
      // Using contravariant Piola transform to map values back to the physical element
      values[0] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
      values[1] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
    }
    
    if (6 <= i && i <= 6)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 6;
    
      // Generate scalings
      const double scalings_y_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
      // Table(s) of coefficients
      static const double coefficients0[1][1] =   \
      {{1.41421356237309}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
    
      // Compute value(s)
      values[2] = coeff0_0*basisvalue0;
    }
    
}

/// Evaluate all basis functions at given point in cell
void saturationequation_1_finite_element_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void saturationequation_1_finite_element_1::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
    
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 1)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 3*num_derivatives; j++)
      values[j] = 0;
    
    if (0 <= i && i <= 5)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_1_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
      const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
      const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    
      // Table(s) of coefficients
      static const double coefficients0[6][3] =   \
      {{0.942809041582063, 0.577350269189626, -0.333333333333333},
      {-0.471404520791032, -0.288675134594813, 0.166666666666667},
      {0.471404520791032, -0.577350269189626, -0.666666666666667},
      {0.471404520791032, 0.288675134594813, 0.833333333333333},
      {-0.471404520791032, -0.288675134594813, 0.166666666666667},
      {0.942809041582063, 0.577350269189626, -0.333333333333333}};
    
      static const double coefficients1[6][3] =   \
      {{-0.471404520791032, 0, -0.333333333333333},
      {0.942809041582063, 0, 0.666666666666667},
      {0.471404520791032, 0, 0.333333333333333},
      {-0.942809041582063, 0, -0.666666666666667},
      {-0.471404520791032, 0.866025403784439, 0.166666666666667},
      {-0.471404520791032, -0.866025403784439, 0.166666666666667}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      static const double dmats0[3][3] =   \
      {{0, 0, 0},
      {4.89897948556636, 0, 0},
      {0, 0, 0}};
    
      static const double dmats1[3][3] =   \
      {{0, 0, 0},
      {2.44948974278318, 0, 0},
      {4.24264068711928, 0, 0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [2*num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
      double coeff0_1 = 0;
      double coeff0_2 = 0;
      double coeff1_0 = 0;
      double coeff1_1 = 0;
      double coeff1_2 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
      double new_coeff0_1 = 0;
      double new_coeff0_2 = 0;
      double new_coeff1_0 = 0;
      double new_coeff1_1 = 0;
      double new_coeff1_2 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
        new_coeff0_1 = coefficients0[dof][1];
        new_coeff0_2 = coefficients0[dof][2];
        new_coeff1_0 = coefficients1[dof][0];
        new_coeff1_1 = coefficients1[dof][1];
        new_coeff1_2 = coefficients1[dof][2];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
          coeff0_1 = new_coeff0_1;
          coeff0_2 = new_coeff0_2;
          coeff1_0 = new_coeff1_0;
          coeff1_1 = new_coeff1_1;
          coeff1_2 = new_coeff1_2;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0];
            new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1];
            new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2];
            new_coeff1_0 = coeff1_0*dmats0[0][0] + coeff1_1*dmats0[1][0] + coeff1_2*dmats0[2][0];
            new_coeff1_1 = coeff1_0*dmats0[0][1] + coeff1_1*dmats0[1][1] + coeff1_2*dmats0[2][1];
            new_coeff1_2 = coeff1_0*dmats0[0][2] + coeff1_1*dmats0[1][2] + coeff1_2*dmats0[2][2];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0];
            new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1];
            new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2];
            new_coeff1_0 = coeff1_0*dmats1[0][0] + coeff1_1*dmats1[1][0] + coeff1_2*dmats1[2][0];
            new_coeff1_1 = coeff1_0*dmats1[0][1] + coeff1_1*dmats1[1][1] + coeff1_2*dmats1[2][1];
            new_coeff1_2 = coeff1_0*dmats1[0][2] + coeff1_1*dmats1[1][2] + coeff1_2*dmats1[2][2];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        // Correct values by the contravariant Piola transform
        const double tmp0_0 = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2;
        const double tmp0_1 = new_coeff1_0*basisvalue0 + new_coeff1_1*basisvalue1 + new_coeff1_2*basisvalue2;
        derivatives[deriv_num] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
        derivatives[num_derivatives + deriv_num] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[row] += transform[row][col]*derivatives[col];
          values[num_derivatives + row] += transform[row][col]*derivatives[num_derivatives + col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
    if (6 <= i && i <= 6)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 6;
    
      // Generate scalings
      const double scalings_y_0 = 1;
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
      // Table(s) of coefficients
      static const double coefficients0[1][1] =   \
      {{1.41421356237309}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      static const double dmats0[1][1] =   \
      {{0}};
    
      static const double dmats1[1][1] =   \
      {{0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        derivatives[deriv_num] = new_coeff0_0*basisvalue0;
      }
    
      // Transform derivatives back to physical element
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        for (unsigned int col = 0; col < num_derivatives; col++)
        {
          values[2*num_derivatives + row] += transform[row][col]*derivatives[col];
        }
      }
      // Delete pointer to array of derivatives on FIAT element
      delete [] derivatives;
    
      // Delete pointer to array of combinations of derivatives and transform
      for (unsigned int row = 0; row < num_derivatives; row++)
      {
        delete [] combinations[row];
        delete [] transform[row];
      }
    
      delete [] combinations;
      delete [] transform;
    }
    
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void saturationequation_1_finite_element_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double saturationequation_1_finite_element_1::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[7][1][2] = {{{0.666666666666667, 0.333333333333333}}, {{0.333333333333333, 0.666666666666667}}, {{0, 0.333333333333333}}, {{0, 0.666666666666667}}, {{0.333333333333333, 0}}, {{0.666666666666667, 0}}, {{0.333333333333333, 0.333333333333333}}};
    static const double W[7][1] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}};
    static const double D[7][1][3] = {{{1, 1, 0}}, {{1, 1, 0}}, {{1, 0, 0}}, {{1, 0, 0}}, {{0, -1, 0}}, {{0, -1, 0}}, {{0, 0, 1}}};
    
    static const unsigned int mappings[7] = {1, 1, 1, 1, 1, 1, 0};
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    const double Jinv_00 =  J_11 / detJ;
    const double Jinv_01 = -J_01 / detJ;
    const double Jinv_10 = -J_10 / detJ;
    const double Jinv_11 =  J_00 / detJ;
    
    double copyofvalues[3];
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[3];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    if (mappings[i] == 0) { 
      // Affine map: Do nothing
    } else if (mappings[i] == 1) {
       // Copy old values:
      copyofvalues[0] = values[0];
      copyofvalues[1] = values[1];
      // Do the inverse of div piola 
      values[0] = detJ*(Jinv_00*copyofvalues[0]+Jinv_01*copyofvalues[1]);
      values[1] = detJ*(Jinv_10*copyofvalues[0]+Jinv_11*copyofvalues[1]); 
    } else { 
       // Other mappings not applicable. 
    }
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 3; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void saturationequation_1_finite_element_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void saturationequation_1_finite_element_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    // Evaluate at vertices and use Piola mapping
    vertex_values[0] = (1.0/detJ)*(dof_values[2]*2*J_00 + dof_values[3]*J_00 + dof_values[4]*(-2*J_01) + dof_values[5]*J_01);
    vertex_values[3] = (1.0/detJ)*(dof_values[0]*2*J_00 + dof_values[1]*J_00 + dof_values[4]*(J_00 + J_01) + dof_values[5]*(2*J_00 - 2*J_01));
    vertex_values[6] = (1.0/detJ)*(dof_values[0]*J_01 + dof_values[1]*2*J_01 + dof_values[2]*(J_00 + J_01) + dof_values[3]*(2*J_00 - 2*J_01));
    vertex_values[1] = (1.0/detJ)*(dof_values[2]*2*J_10 + dof_values[3]*J_10 + dof_values[4]*(-2*J_11) + dof_values[5]*J_11);
    vertex_values[4] = (1.0/detJ)*(dof_values[0]*2*J_10 + dof_values[1]*J_10 + dof_values[4]*(J_10 + J_11) + dof_values[5]*(2*J_10 - 2*J_11));
    vertex_values[7] = (1.0/detJ)*(dof_values[0]*J_11 + dof_values[1]*2*J_11 + dof_values[2]*(J_10 + J_11) + dof_values[3]*(2*J_10 - 2*J_11));
    // Evaluate at vertices and use affine mapping
    vertex_values[2] = dof_values[6];
    vertex_values[5] = dof_values[6];
    vertex_values[8] = dof_values[6];
}

/// Return the number of sub elements (for a mixed element)
unsigned int saturationequation_1_finite_element_1::num_sub_elements() const
{
    return 2;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* saturationequation_1_finite_element_1::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new saturationequation_1_finite_element_1_0();
      break;
    case 1:
      return new saturationequation_1_finite_element_1_1();
      break;
    }
    return 0;
}


/// Constructor
saturationequation_1_finite_element_2::saturationequation_1_finite_element_2() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
saturationequation_1_finite_element_2::~saturationequation_1_finite_element_2()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* saturationequation_1_finite_element_2::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return the cell shape
ufc::shape saturationequation_1_finite_element_2::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int saturationequation_1_finite_element_2::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int saturationequation_1_finite_element_2::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int saturationequation_1_finite_element_2::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void saturationequation_1_finite_element_2::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void saturationequation_1_finite_element_2::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void saturationequation_1_finite_element_2::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
    
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 1)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[1][1] = \
    {{0}};
    
    static const double dmats1[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void saturationequation_1_finite_element_2::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double saturationequation_1_finite_element_2::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[1][1][2] = {{{0.333333333333333, 0.333333333333333}}};
    static const double W[1][1] = {{1}};
    static const double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void saturationequation_1_finite_element_2::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void saturationequation_1_finite_element_2::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int saturationequation_1_finite_element_2::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* saturationequation_1_finite_element_2::create_sub_element(unsigned int i) const
{
    return new saturationequation_1_finite_element_2();
}


/// Constructor
saturationequation_1_finite_element_3::saturationequation_1_finite_element_3() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
saturationequation_1_finite_element_3::~saturationequation_1_finite_element_3()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* saturationequation_1_finite_element_3::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return the cell shape
ufc::shape saturationequation_1_finite_element_3::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int saturationequation_1_finite_element_3::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int saturationequation_1_finite_element_3::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int saturationequation_1_finite_element_3::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void saturationequation_1_finite_element_3::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void saturationequation_1_finite_element_3::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void saturationequation_1_finite_element_3::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
    
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 1)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[1][1] = \
    {{0}};
    
    static const double dmats1[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void saturationequation_1_finite_element_3::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double saturationequation_1_finite_element_3::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[1][1][2] = {{{0.333333333333333, 0.333333333333333}}};
    static const double W[1][1] = {{1}};
    static const double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void saturationequation_1_finite_element_3::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void saturationequation_1_finite_element_3::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int saturationequation_1_finite_element_3::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* saturationequation_1_finite_element_3::create_sub_element(unsigned int i) const
{
    return new saturationequation_1_finite_element_3();
}


/// Constructor
saturationequation_1_finite_element_4::saturationequation_1_finite_element_4() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
saturationequation_1_finite_element_4::~saturationequation_1_finite_element_4()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* saturationequation_1_finite_element_4::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return the cell shape
ufc::shape saturationequation_1_finite_element_4::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int saturationequation_1_finite_element_4::space_dimension() const
{
    return 3;
}

/// Return the rank of the value space
unsigned int saturationequation_1_finite_element_4::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int saturationequation_1_finite_element_4::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void saturationequation_1_finite_element_4::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    
    // Table(s) of coefficients
    static const double coefficients0[3][3] = \
    {{0.471404520791032, -0.288675134594813, -0.166666666666667},
    {0.471404520791032, 0.288675134594813, -0.166666666666667},
    {0.471404520791032, 0, 0.333333333333333}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2;
}

/// Evaluate all basis functions at given point in cell
void saturationequation_1_finite_element_4::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void saturationequation_1_finite_element_4::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
    
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 1)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    
    // Table(s) of coefficients
    static const double coefficients0[3][3] = \
    {{0.471404520791032, -0.288675134594813, -0.166666666666667},
    {0.471404520791032, 0.288675134594813, -0.166666666666667},
    {0.471404520791032, 0, 0.333333333333333}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[3][3] = \
    {{0, 0, 0},
    {4.89897948556636, 0, 0},
    {0, 0, 0}};
    
    static const double dmats1[3][3] = \
    {{0, 0, 0},
    {2.44948974278318, 0, 0},
    {4.24264068711928, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void saturationequation_1_finite_element_4::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double saturationequation_1_finite_element_4::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[3][1][2] = {{{0, 0}}, {{1, 0}}, {{0, 1}}};
    static const double W[3][1] = {{1}, {1}, {1}};
    static const double D[3][1][1] = {{{1}}, {{1}}, {{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void saturationequation_1_finite_element_4::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void saturationequation_1_finite_element_4::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
}

/// Return the number of sub elements (for a mixed element)
unsigned int saturationequation_1_finite_element_4::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* saturationequation_1_finite_element_4::create_sub_element(unsigned int i) const
{
    return new saturationequation_1_finite_element_4();
}


/// Constructor
saturationequation_1_finite_element_5::saturationequation_1_finite_element_5() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
saturationequation_1_finite_element_5::~saturationequation_1_finite_element_5()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* saturationequation_1_finite_element_5::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return the cell shape
ufc::shape saturationequation_1_finite_element_5::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int saturationequation_1_finite_element_5::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int saturationequation_1_finite_element_5::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int saturationequation_1_finite_element_5::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void saturationequation_1_finite_element_5::evaluate_basis(unsigned int i,
                                   double* values,
                                   const double* coordinates,
                                   const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
}

/// Evaluate all basis functions at given point in cell
void saturationequation_1_finite_element_5::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void saturationequation_1_finite_element_5::evaluate_basis_derivatives(unsigned int i,
                                               unsigned int n,
                                               double* values,
                                               const double* coordinates,
                                               const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
    
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 1)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
    
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    
    // Table(s) of coefficients
    static const double coefficients0[1][1] = \
    {{1.41421356237309}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[1][1] = \
    {{0}};
    
    static const double dmats1[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void saturationequation_1_finite_element_5::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double saturationequation_1_finite_element_5::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[1][1][2] = {{{0.333333333333333, 0.333333333333333}}};
    static const double W[1][1] = {{1}};
    static const double D[1][1][1] = {{{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void saturationequation_1_finite_element_5::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void saturationequation_1_finite_element_5::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
}

/// Return the number of sub elements (for a mixed element)
unsigned int saturationequation_1_finite_element_5::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* saturationequation_1_finite_element_5::create_sub_element(unsigned int i) const
{
    return new saturationequation_1_finite_element_5();
}

/// Constructor
saturationequation_1_dof_map_0::saturationequation_1_dof_map_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
saturationequation_1_dof_map_0::~saturationequation_1_dof_map_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* saturationequation_1_dof_map_0::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool saturationequation_1_dof_map_0::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool saturationequation_1_dof_map_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void saturationequation_1_dof_map_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void saturationequation_1_dof_map_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int saturationequation_1_dof_map_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int saturationequation_1_dof_map_0::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int saturationequation_1_dof_map_0::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int saturationequation_1_dof_map_0::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int saturationequation_1_dof_map_0::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int saturationequation_1_dof_map_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void saturationequation_1_dof_map_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[2][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void saturationequation_1_dof_map_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void saturationequation_1_dof_map_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void saturationequation_1_dof_map_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.333333333333333*x[0][0] + 0.333333333333333*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[0][1] = 0.333333333333333*x[0][1] + 0.333333333333333*x[1][1] + 0.333333333333333*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int saturationequation_1_dof_map_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* saturationequation_1_dof_map_0::create_sub_dof_map(unsigned int i) const
{
    return new saturationequation_1_dof_map_0();
}


/// Constructor
saturationequation_1_dof_map_1_0::saturationequation_1_dof_map_1_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
saturationequation_1_dof_map_1_0::~saturationequation_1_dof_map_1_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* saturationequation_1_dof_map_1_0::signature() const
{
    return "FFC dof map for FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool saturationequation_1_dof_map_1_0::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return true;
      break;
    case 2:
      return false;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool saturationequation_1_dof_map_1_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 2*m.num_entities[1];
    return false;
}

/// Initialize dof map for given cell
void saturationequation_1_dof_map_1_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void saturationequation_1_dof_map_1_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int saturationequation_1_dof_map_1_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int saturationequation_1_dof_map_1_0::local_dimension(const ufc::cell& c) const
{
    return 6;
}

/// Return the maximum dimension of the local finite element function space
unsigned int saturationequation_1_dof_map_1_0::max_local_dimension() const
{
    return 6;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int saturationequation_1_dof_map_1_0::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int saturationequation_1_dof_map_1_0::num_facet_dofs() const
{
    return 2;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int saturationequation_1_dof_map_1_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void saturationequation_1_dof_map_1_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 2*c.entity_indices[1][0];
    dofs[1] = 2*c.entity_indices[1][0] + 1;
    dofs[2] = 2*c.entity_indices[1][1];
    dofs[3] = 2*c.entity_indices[1][1] + 1;
    dofs[4] = 2*c.entity_indices[1][2];
    dofs[5] = 2*c.entity_indices[1][2] + 1;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void saturationequation_1_dof_map_1_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 0;
      dofs[1] = 1;
      break;
    case 1:
      dofs[0] = 2;
      dofs[1] = 3;
      break;
    case 2:
      dofs[0] = 4;
      dofs[1] = 5;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void saturationequation_1_dof_map_1_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void saturationequation_1_dof_map_1_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.666666666666667*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[0][1] = 0.666666666666667*x[1][1] + 0.333333333333333*x[2][1];
    coordinates[1][0] = 0.333333333333333*x[1][0] + 0.666666666666667*x[2][0];
    coordinates[1][1] = 0.333333333333333*x[1][1] + 0.666666666666667*x[2][1];
    coordinates[2][0] = 0.666666666666667*x[0][0] + 0.333333333333333*x[2][0];
    coordinates[2][1] = 0.666666666666667*x[0][1] + 0.333333333333333*x[2][1];
    coordinates[3][0] = 0.333333333333333*x[0][0] + 0.666666666666667*x[2][0];
    coordinates[3][1] = 0.333333333333333*x[0][1] + 0.666666666666667*x[2][1];
    coordinates[4][0] = 0.666666666666667*x[0][0] + 0.333333333333333*x[1][0];
    coordinates[4][1] = 0.666666666666667*x[0][1] + 0.333333333333333*x[1][1];
    coordinates[5][0] = 0.333333333333333*x[0][0] + 0.666666666666667*x[1][0];
    coordinates[5][1] = 0.333333333333333*x[0][1] + 0.666666666666667*x[1][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int saturationequation_1_dof_map_1_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* saturationequation_1_dof_map_1_0::create_sub_dof_map(unsigned int i) const
{
    return new saturationequation_1_dof_map_1_0();
}


/// Constructor
saturationequation_1_dof_map_1_1::saturationequation_1_dof_map_1_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
saturationequation_1_dof_map_1_1::~saturationequation_1_dof_map_1_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* saturationequation_1_dof_map_1_1::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool saturationequation_1_dof_map_1_1::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool saturationequation_1_dof_map_1_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void saturationequation_1_dof_map_1_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void saturationequation_1_dof_map_1_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int saturationequation_1_dof_map_1_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int saturationequation_1_dof_map_1_1::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int saturationequation_1_dof_map_1_1::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int saturationequation_1_dof_map_1_1::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int saturationequation_1_dof_map_1_1::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int saturationequation_1_dof_map_1_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void saturationequation_1_dof_map_1_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[2][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void saturationequation_1_dof_map_1_1::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void saturationequation_1_dof_map_1_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void saturationequation_1_dof_map_1_1::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.333333333333333*x[0][0] + 0.333333333333333*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[0][1] = 0.333333333333333*x[0][1] + 0.333333333333333*x[1][1] + 0.333333333333333*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int saturationequation_1_dof_map_1_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* saturationequation_1_dof_map_1_1::create_sub_dof_map(unsigned int i) const
{
    return new saturationequation_1_dof_map_1_1();
}


/// Constructor
saturationequation_1_dof_map_1::saturationequation_1_dof_map_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
saturationequation_1_dof_map_1::~saturationequation_1_dof_map_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* saturationequation_1_dof_map_1::signature() const
{
    return "FFC dof map for MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) })";
}

/// Return true iff mesh entities of topological dimension d are needed
bool saturationequation_1_dof_map_1::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return true;
      break;
    case 2:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool saturationequation_1_dof_map_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 2*m.num_entities[1] + m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void saturationequation_1_dof_map_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void saturationequation_1_dof_map_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int saturationequation_1_dof_map_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int saturationequation_1_dof_map_1::local_dimension(const ufc::cell& c) const
{
    return 7;
}

/// Return the maximum dimension of the local finite element function space
unsigned int saturationequation_1_dof_map_1::max_local_dimension() const
{
    return 7;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int saturationequation_1_dof_map_1::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int saturationequation_1_dof_map_1::num_facet_dofs() const
{
    return 2;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int saturationequation_1_dof_map_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void saturationequation_1_dof_map_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 2*c.entity_indices[1][0];
    dofs[1] = 2*c.entity_indices[1][0] + 1;
    dofs[2] = 2*c.entity_indices[1][1];
    dofs[3] = 2*c.entity_indices[1][1] + 1;
    dofs[4] = 2*c.entity_indices[1][2];
    dofs[5] = 2*c.entity_indices[1][2] + 1;
    unsigned int offset = 2*m.num_entities[1];
    dofs[6] = offset + c.entity_indices[2][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void saturationequation_1_dof_map_1::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 0;
      dofs[1] = 1;
      break;
    case 1:
      dofs[0] = 2;
      dofs[1] = 3;
      break;
    case 2:
      dofs[0] = 4;
      dofs[1] = 5;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void saturationequation_1_dof_map_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void saturationequation_1_dof_map_1::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.666666666666667*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[0][1] = 0.666666666666667*x[1][1] + 0.333333333333333*x[2][1];
    coordinates[1][0] = 0.333333333333333*x[1][0] + 0.666666666666667*x[2][0];
    coordinates[1][1] = 0.333333333333333*x[1][1] + 0.666666666666667*x[2][1];
    coordinates[2][0] = 0.666666666666667*x[0][0] + 0.333333333333333*x[2][0];
    coordinates[2][1] = 0.666666666666667*x[0][1] + 0.333333333333333*x[2][1];
    coordinates[3][0] = 0.333333333333333*x[0][0] + 0.666666666666667*x[2][0];
    coordinates[3][1] = 0.333333333333333*x[0][1] + 0.666666666666667*x[2][1];
    coordinates[4][0] = 0.666666666666667*x[0][0] + 0.333333333333333*x[1][0];
    coordinates[4][1] = 0.666666666666667*x[0][1] + 0.333333333333333*x[1][1];
    coordinates[5][0] = 0.333333333333333*x[0][0] + 0.666666666666667*x[1][0];
    coordinates[5][1] = 0.333333333333333*x[0][1] + 0.666666666666667*x[1][1];
    coordinates[6][0] = 0.333333333333333*x[0][0] + 0.333333333333333*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[6][1] = 0.333333333333333*x[0][1] + 0.333333333333333*x[1][1] + 0.333333333333333*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int saturationequation_1_dof_map_1::num_sub_dof_maps() const
{
    return 2;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* saturationequation_1_dof_map_1::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new saturationequation_1_dof_map_1_0();
      break;
    case 1:
      return new saturationequation_1_dof_map_1_1();
      break;
    }
    return 0;
}


/// Constructor
saturationequation_1_dof_map_2::saturationequation_1_dof_map_2() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
saturationequation_1_dof_map_2::~saturationequation_1_dof_map_2()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* saturationequation_1_dof_map_2::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool saturationequation_1_dof_map_2::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool saturationequation_1_dof_map_2::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void saturationequation_1_dof_map_2::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void saturationequation_1_dof_map_2::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int saturationequation_1_dof_map_2::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int saturationequation_1_dof_map_2::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int saturationequation_1_dof_map_2::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int saturationequation_1_dof_map_2::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int saturationequation_1_dof_map_2::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int saturationequation_1_dof_map_2::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void saturationequation_1_dof_map_2::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[2][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void saturationequation_1_dof_map_2::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void saturationequation_1_dof_map_2::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void saturationequation_1_dof_map_2::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.333333333333333*x[0][0] + 0.333333333333333*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[0][1] = 0.333333333333333*x[0][1] + 0.333333333333333*x[1][1] + 0.333333333333333*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int saturationequation_1_dof_map_2::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* saturationequation_1_dof_map_2::create_sub_dof_map(unsigned int i) const
{
    return new saturationequation_1_dof_map_2();
}


/// Constructor
saturationequation_1_dof_map_3::saturationequation_1_dof_map_3() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
saturationequation_1_dof_map_3::~saturationequation_1_dof_map_3()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* saturationequation_1_dof_map_3::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool saturationequation_1_dof_map_3::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool saturationequation_1_dof_map_3::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void saturationequation_1_dof_map_3::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void saturationequation_1_dof_map_3::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int saturationequation_1_dof_map_3::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int saturationequation_1_dof_map_3::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int saturationequation_1_dof_map_3::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int saturationequation_1_dof_map_3::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int saturationequation_1_dof_map_3::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int saturationequation_1_dof_map_3::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void saturationequation_1_dof_map_3::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[2][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void saturationequation_1_dof_map_3::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void saturationequation_1_dof_map_3::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void saturationequation_1_dof_map_3::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.333333333333333*x[0][0] + 0.333333333333333*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[0][1] = 0.333333333333333*x[0][1] + 0.333333333333333*x[1][1] + 0.333333333333333*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int saturationequation_1_dof_map_3::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* saturationequation_1_dof_map_3::create_sub_dof_map(unsigned int i) const
{
    return new saturationequation_1_dof_map_3();
}


/// Constructor
saturationequation_1_dof_map_4::saturationequation_1_dof_map_4() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
saturationequation_1_dof_map_4::~saturationequation_1_dof_map_4()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* saturationequation_1_dof_map_4::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool saturationequation_1_dof_map_4::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool saturationequation_1_dof_map_4::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void saturationequation_1_dof_map_4::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void saturationequation_1_dof_map_4::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int saturationequation_1_dof_map_4::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int saturationequation_1_dof_map_4::local_dimension(const ufc::cell& c) const
{
    return 3;
}

/// Return the maximum dimension of the local finite element function space
unsigned int saturationequation_1_dof_map_4::max_local_dimension() const
{
    return 3;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int saturationequation_1_dof_map_4::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int saturationequation_1_dof_map_4::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int saturationequation_1_dof_map_4::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void saturationequation_1_dof_map_4::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[2][0];
    dofs[1] = 3*c.entity_indices[2][0] + 1;
    dofs[2] = 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void saturationequation_1_dof_map_4::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void saturationequation_1_dof_map_4::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void saturationequation_1_dof_map_4::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = x[0][0];
    coordinates[0][1] = x[0][1];
    coordinates[1][0] = x[1][0];
    coordinates[1][1] = x[1][1];
    coordinates[2][0] = x[2][0];
    coordinates[2][1] = x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int saturationequation_1_dof_map_4::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* saturationequation_1_dof_map_4::create_sub_dof_map(unsigned int i) const
{
    return new saturationequation_1_dof_map_4();
}


/// Constructor
saturationequation_1_dof_map_5::saturationequation_1_dof_map_5() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
saturationequation_1_dof_map_5::~saturationequation_1_dof_map_5()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* saturationequation_1_dof_map_5::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool saturationequation_1_dof_map_5::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return true;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool saturationequation_1_dof_map_5::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void saturationequation_1_dof_map_5::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void saturationequation_1_dof_map_5::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int saturationequation_1_dof_map_5::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int saturationequation_1_dof_map_5::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int saturationequation_1_dof_map_5::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int saturationequation_1_dof_map_5::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int saturationequation_1_dof_map_5::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int saturationequation_1_dof_map_5::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void saturationequation_1_dof_map_5::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[2][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void saturationequation_1_dof_map_5::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void saturationequation_1_dof_map_5::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void saturationequation_1_dof_map_5::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.333333333333333*x[0][0] + 0.333333333333333*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[0][1] = 0.333333333333333*x[0][1] + 0.333333333333333*x[1][1] + 0.333333333333333*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int saturationequation_1_dof_map_5::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* saturationequation_1_dof_map_5::create_sub_dof_map(unsigned int i) const
{
    return new saturationequation_1_dof_map_5();
}


/// Constructor
saturationequation_1_cell_integral_0_quadrature::saturationequation_1_cell_integral_0_quadrature() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
saturationequation_1_cell_integral_0_quadrature::~saturationequation_1_cell_integral_0_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void saturationequation_1_cell_integral_0_quadrature::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c) const
{
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    const double Jinv_00 =  J_11 / detJ;
    const double Jinv_01 = -J_01 / detJ;
    const double Jinv_10 = -J_10 / detJ;
    const double Jinv_11 =  J_00 / detJ;
    
    // Set scale factor
    const double det = std::abs(detJ);
    
    
    // Array of quadrature weights
    static const double W1 = 0.5;
    // Quadrature points on the UFC reference element: (0.333333333333333, 0.333333333333333)
    
    // Value of basis functions at quadrature points.
    static const double FE0[1][1] = \
    {{1}};
    
    static const double FE0_D01[1][1] = \
    {{0}};
    
    static const double FE1_C0[1][7] = \
    {{0.666666666666667, -0.333333333333333, 0.333333333333333, 0.333333333333333, -0.333333333333333, 0.666666666666667, 0}};
    
    static const double FE1_C1[1][7] = \
    {{-0.333333333333333, 0.666666666666667, 0.333333333333333, -0.666666666666667, -0.333333333333333, -0.333333333333333, 0}};
    
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', False), ('ignore zero tables', False), ('non zero columns', False), ('remove zero terms', False), ('ignore ones', False)
    // Total number of operations to compute element tensor: 61
    
    // Loop quadrature points for integral
    // Number of operations to compute element tensor for following IP loop = 61
    // Only 1 integration point, omitting IP loop.
    
    // Function declarations
    double F0 = 0;
    double F1 = 0;
    double F2 = 0;
    double F3 = 0;
    double F4 = 0;
    
    // Total number of operations to compute function values = 6
    for (unsigned int r = 0; r < 1; r++)
    {
      F0 += FE0[0][r]*w[1][r];
      F1 += FE0[0][r]*w[2][r];
      F4 += FE0[0][r]*w[4][r];
    }// end loop over 'r'
    
    // Total number of operations to compute function values = 28
    for (unsigned int r = 0; r < 7; r++)
    {
      F2 += FE1_C0[0][r]*w[0][r];
      F3 += FE1_C1[0][r]*w[0][r];
    }// end loop over 'r'
    
    // Number of operations for primary indices: 27
    for (unsigned int j = 0; j < 1; j++)
    {
      // Number of operations to compute entry: 27
      A[j] += (FE0[0][j]*F0 + ((Jinv_00*FE0_D01[0][j] + Jinv_10*FE0_D01[0][j])*F1*((1.0/detJ)*J_00*F2 + (1.0/detJ)*J_01*F3) + (Jinv_01*FE0_D01[0][j] + Jinv_11*FE0_D01[0][j])*F1*((1.0/detJ)*J_10*F2 + (1.0/detJ)*J_11*F3))*F4)*W1*det;
    }// end loop over 'j'
}

/// Constructor
saturationequation_1_cell_integral_0::saturationequation_1_cell_integral_0() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
saturationequation_1_cell_integral_0::~saturationequation_1_cell_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void saturationequation_1_cell_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c) const
{
    // Reset values of the element tensor block
    A[0] = 0;
    
    // Add all contributions to element tensor
    integral_0_quadrature.tabulate_tensor(A, w, c);
}

/// Constructor
saturationequation_1_exterior_facet_integral_0_quadrature::saturationequation_1_exterior_facet_integral_0_quadrature() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
saturationequation_1_exterior_facet_integral_0_quadrature::~saturationequation_1_exterior_facet_integral_0_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void saturationequation_1_exterior_facet_integral_0_quadrature::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c,
                                    unsigned int facet) const
{
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Vertices on edges
    static unsigned int edge_vertices[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    
    // Get vertices
    const unsigned int v0 = edge_vertices[facet][0];
    const unsigned int v1 = edge_vertices[facet][1];
    
    // Compute scale factor (length of edge scaled by length of reference interval)
    const double dx0 = x[v1][0] - x[v0][0];
    const double dx1 = x[v1][1] - x[v0][1];
    const double det = std::sqrt(dx0*dx0 + dx1*dx1);
    
    const bool direction = dx1*(x[facet][0] - x[v0][0]) - dx0*(x[facet][1] - x[v0][1]) < 0;
    
    // Compute facet normals from the facet scale factor constants
    const double n0 = direction ? dx1 / det : -dx1 / det;
    const double n1 = direction ? -dx0 / det : dx0 / det;
    
    
    // Array of quadrature weights
    static const double W1 = 1;
    // Quadrature points on the UFC reference element: (0.5)
    
    // Value of basis functions at quadrature points.
    static const double FE0_f0[1][1] = \
    {{1}};
    
    static const double FE1_f0_C0[1][7] = \
    {{1, -0.5, -0.5, 1, -0.5, 1, 0}};
    
    static const double FE1_f0_C1[1][7] = \
    {{-0.5, 1, 0.5, -1, 0.5, -1, 0}};
    
    static const double FE1_f1_C0[1][7] = \
    {{0, 0, 0.5, 0.5, 0, 0, 0}};
    
    static const double FE1_f1_C1[1][7] = \
    {{-0.5, 1, 0.5, -1, -1, 0.5, 0}};
    
    static const double FE1_f2_C0[1][7] = \
    {{1, -0.5, 1, -0.5, -0.5, 1, 0}};
    
    static const double FE1_f2_C1[1][7] = \
    {{0, 0, 0, 0, -0.5, -0.5, 0}};
    
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', False), ('ignore zero tables', False), ('non zero columns', False), ('remove zero terms', False), ('ignore ones', False)
    switch ( facet )
    {
    case 0:
      {
      // Total number of operations to compute element tensor (from this point): 66
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 66
      // Only 1 integration point, omitting IP loop.
      
      // Function declarations
      double F0 = 0;
      double F1 = 0;
      double F2 = 0;
      double F3 = 0;
      
      // Total number of operations to compute function values = 4
      for (unsigned int r = 0; r < 1; r++)
      {
        F0 += FE0_f0[0][r]*w[4][r];
        F3 += FE0_f0[0][r]*w[2][r];
      }// end loop over 'r'
      
      // Total number of operations to compute function values = 28
      for (unsigned int r = 0; r < 7; r++)
      {
        F1 += FE1_f0_C0[0][r]*w[0][r];
        F2 += FE1_f0_C1[0][r]*w[0][r];
      }// end loop over 'r'
      
      // Number of operations for primary indices: 34
      for (unsigned int j = 0; j < 1; j++)
      {
        // Number of operations to compute entry: 34
        A[j] += FE0_f0[0][j]*((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2)) + std::abs((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2))))/(2)*F3*F0*-1*W1*det;
      }// end loop over 'j'
      }
      break;
    case 1:
      {
      // Total number of operations to compute element tensor (from this point): 66
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 66
      // Only 1 integration point, omitting IP loop.
      
      // Function declarations
      double F0 = 0;
      double F1 = 0;
      double F2 = 0;
      double F3 = 0;
      
      // Total number of operations to compute function values = 4
      for (unsigned int r = 0; r < 1; r++)
      {
        F0 += FE0_f0[0][r]*w[4][r];
        F3 += FE0_f0[0][r]*w[2][r];
      }// end loop over 'r'
      
      // Total number of operations to compute function values = 28
      for (unsigned int r = 0; r < 7; r++)
      {
        F1 += FE1_f1_C0[0][r]*w[0][r];
        F2 += FE1_f1_C1[0][r]*w[0][r];
      }// end loop over 'r'
      
      // Number of operations for primary indices: 34
      for (unsigned int j = 0; j < 1; j++)
      {
        // Number of operations to compute entry: 34
        A[j] += FE0_f0[0][j]*((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2)) + std::abs((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2))))/(2)*F3*F0*-1*W1*det;
      }// end loop over 'j'
      }
      break;
    case 2:
      {
      // Total number of operations to compute element tensor (from this point): 66
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 66
      // Only 1 integration point, omitting IP loop.
      
      // Function declarations
      double F0 = 0;
      double F1 = 0;
      double F2 = 0;
      double F3 = 0;
      
      // Total number of operations to compute function values = 4
      for (unsigned int r = 0; r < 1; r++)
      {
        F0 += FE0_f0[0][r]*w[4][r];
        F3 += FE0_f0[0][r]*w[2][r];
      }// end loop over 'r'
      
      // Total number of operations to compute function values = 28
      for (unsigned int r = 0; r < 7; r++)
      {
        F1 += FE1_f2_C0[0][r]*w[0][r];
        F2 += FE1_f2_C1[0][r]*w[0][r];
      }// end loop over 'r'
      
      // Number of operations for primary indices: 34
      for (unsigned int j = 0; j < 1; j++)
      {
        // Number of operations to compute entry: 34
        A[j] += FE0_f0[0][j]*((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2)) + std::abs((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2))))/(2)*F3*F0*-1*W1*det;
      }// end loop over 'j'
      }
      break;
    }
}

/// Constructor
saturationequation_1_exterior_facet_integral_0::saturationequation_1_exterior_facet_integral_0() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
saturationequation_1_exterior_facet_integral_0::~saturationequation_1_exterior_facet_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void saturationequation_1_exterior_facet_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c,
                                    unsigned int facet) const
{
    // Reset values of the element tensor block
    A[0] = 0;
    
    // Add all contributions to element tensor
    integral_0_quadrature.tabulate_tensor(A, w, c, facet);
}

/// Constructor
saturationequation_1_exterior_facet_integral_1_quadrature::saturationequation_1_exterior_facet_integral_1_quadrature() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
saturationequation_1_exterior_facet_integral_1_quadrature::~saturationequation_1_exterior_facet_integral_1_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void saturationequation_1_exterior_facet_integral_1_quadrature::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c,
                                    unsigned int facet) const
{
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Vertices on edges
    static unsigned int edge_vertices[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    
    // Get vertices
    const unsigned int v0 = edge_vertices[facet][0];
    const unsigned int v1 = edge_vertices[facet][1];
    
    // Compute scale factor (length of edge scaled by length of reference interval)
    const double dx0 = x[v1][0] - x[v0][0];
    const double dx1 = x[v1][1] - x[v0][1];
    const double det = std::sqrt(dx0*dx0 + dx1*dx1);
    
    const bool direction = dx1*(x[facet][0] - x[v0][0]) - dx0*(x[facet][1] - x[v0][1]) < 0;
    
    // Compute facet normals from the facet scale factor constants
    const double n0 = direction ? dx1 / det : -dx1 / det;
    const double n1 = direction ? -dx0 / det : dx0 / det;
    
    
    // Array of quadrature weights
    static const double W2[2] = {0.5, 0.5};
    // Quadrature points on the UFC reference element: (0.211324865405187), (0.788675134594813)
    
    // Value of basis functions at quadrature points.
    static const double FE2_f0[2][1] = \
    {{1},
    {1}};
    
    static const double FE3_f0[2][3] = \
    {{0, 0.788675134594813, 0.211324865405187},
    {0, 0.211324865405187, 0.788675134594813}};
    
    static const double FE3_f1[2][3] = \
    {{0.788675134594813, 0, 0.211324865405187},
    {0.211324865405187, 0, 0.788675134594813}};
    
    static const double FE3_f2[2][3] = \
    {{0.788675134594813, 0.211324865405187, 0},
    {0.211324865405187, 0.788675134594813, 0}};
    
    static const double FE4_f0_C0[2][7] = \
    {{1.57735026918963, -0.788675134594813, -0.211324865405187, 0.422649730810374, -0.788675134594813, 1.57735026918963, 0},
    {0.422649730810374, -0.211324865405187, -0.788675134594813, 1.57735026918963, -0.211324865405187, 0.422649730810374, 0}};
    
    static const double FE4_f0_C1[2][7] = \
    {{-0.211324865405187, 0.422649730810374, 0.211324865405187, -0.422649730810374, 0.788675134594813, -1.57735026918963, 0},
    {-0.788675134594813, 1.57735026918963, 0.788675134594813, -1.57735026918963, 0.211324865405187, -0.422649730810374, 0}};
    
    static const double FE4_f1_C0[2][7] = \
    {{0, 0, 1.36602540378444, -0.366025403784439, 0, 0, 0},
    {0, 0, -0.366025403784439, 1.36602540378444, 0, 0, 0}};
    
    static const double FE4_f1_C1[2][7] = \
    {{-0.211324865405187, 0.422649730810374, 0.211324865405187, -0.422649730810374, -1.57735026918963, 0.788675134594813, 0},
    {-0.788675134594813, 1.57735026918963, 0.788675134594813, -1.57735026918963, -0.422649730810374, 0.211324865405187, 0}};
    
    static const double FE4_f2_C0[2][7] = \
    {{0.422649730810375, -0.211324865405187, 1.57735026918963, -0.788675134594813, -0.211324865405187, 0.422649730810375, 0},
    {1.57735026918963, -0.788675134594813, 0.422649730810374, -0.211324865405187, -0.788675134594813, 1.57735026918963, 0}};
    
    static const double FE4_f2_C1[2][7] = \
    {{0, 0, 0, 0, -1.36602540378444, 0.366025403784439, 0},
    {0, 0, 0, 0, 0.366025403784439, -1.36602540378444, 0}};
    
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', False), ('ignore zero tables', False), ('non zero columns', False), ('remove zero terms', False), ('ignore ones', False)
    switch ( facet )
    {
    case 0:
      {
      // Total number of operations to compute element tensor (from this point): 142
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 142
      for (unsigned int ip = 0; ip < 2; ip++)
      {
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        
        // Total number of operations to compute function values = 2
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE2_f0[ip][r]*w[4][r];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 3; r++)
        {
          F3 += FE3_f0[ip][r]*w[3][r];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 28
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE4_f0_C0[ip][r]*w[0][r];
          F2 += FE4_f0_C1[ip][r]*w[0][r];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 35
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 35
          A[j] += FE2_f0[ip][j]*((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2)) + -1*std::abs((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2))))/(2)*F3*F0*-1*W2[ip]*det;
        }// end loop over 'j'
      }// end loop over 'ip'
      }
      break;
    case 1:
      {
      // Total number of operations to compute element tensor (from this point): 142
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 142
      for (unsigned int ip = 0; ip < 2; ip++)
      {
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        
        // Total number of operations to compute function values = 2
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE2_f0[ip][r]*w[4][r];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 3; r++)
        {
          F3 += FE3_f1[ip][r]*w[3][r];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 28
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE4_f1_C0[ip][r]*w[0][r];
          F2 += FE4_f1_C1[ip][r]*w[0][r];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 35
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 35
          A[j] += FE2_f0[ip][j]*((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2)) + -1*std::abs((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2))))/(2)*F3*F0*-1*W2[ip]*det;
        }// end loop over 'j'
      }// end loop over 'ip'
      }
      break;
    case 2:
      {
      // Total number of operations to compute element tensor (from this point): 142
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 142
      for (unsigned int ip = 0; ip < 2; ip++)
      {
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        
        // Total number of operations to compute function values = 2
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE2_f0[ip][r]*w[4][r];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 3; r++)
        {
          F3 += FE3_f2[ip][r]*w[3][r];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 28
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE4_f2_C0[ip][r]*w[0][r];
          F2 += FE4_f2_C1[ip][r]*w[0][r];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 35
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 35
          A[j] += FE2_f0[ip][j]*((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2)) + -1*std::abs((n1*((1.0/detJ)*J_10*F1 + (1.0/detJ)*J_11*F2) + n0*((1.0/detJ)*J_00*F1 + (1.0/detJ)*J_01*F2))))/(2)*F3*F0*-1*W2[ip]*det;
        }// end loop over 'j'
      }// end loop over 'ip'
      }
      break;
    }
}

/// Constructor
saturationequation_1_exterior_facet_integral_1::saturationequation_1_exterior_facet_integral_1() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
saturationequation_1_exterior_facet_integral_1::~saturationequation_1_exterior_facet_integral_1()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void saturationequation_1_exterior_facet_integral_1::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c,
                                    unsigned int facet) const
{
    // Reset values of the element tensor block
    A[0] = 0;
    
    // Add all contributions to element tensor
    integral_1_quadrature.tabulate_tensor(A, w, c, facet);
}

/// Constructor
saturationequation_1_interior_facet_integral_0_quadrature::saturationequation_1_interior_facet_integral_0_quadrature() : ufc::interior_facet_integral()
{
    // Do nothing
}

/// Destructor
saturationequation_1_interior_facet_integral_0_quadrature::~saturationequation_1_interior_facet_integral_0_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local interior facet
void saturationequation_1_interior_facet_integral_0_quadrature::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c0,
                                    const ufc::cell& c1,
                                    unsigned int facet0,
                                    unsigned int facet1) const
{
    // Extract vertex coordinates
    const double * const * x0 = c0.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J0_00 = x0[1][0] - x0[0][0];
    const double J0_01 = x0[2][0] - x0[0][0];
    const double J0_10 = x0[1][1] - x0[0][1];
    const double J0_11 = x0[2][1] - x0[0][1];
    
    // Compute determinant of Jacobian
    double detJ0 = J0_00*J0_11 - J0_01*J0_10;
    
    // Compute inverse of Jacobian
    
    // Extract vertex coordinates
    const double * const * x1 = c1.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J1_00 = x1[1][0] - x1[0][0];
    const double J1_01 = x1[2][0] - x1[0][0];
    const double J1_10 = x1[1][1] - x1[0][1];
    const double J1_11 = x1[2][1] - x1[0][1];
    
    // Compute determinant of Jacobian
    double detJ1 = J1_00*J1_11 - J1_01*J1_10;
    
    // Compute inverse of Jacobian
    
    // Vertices on edges
    static unsigned int edge_vertices[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    
    // Get vertices
    const unsigned int v0 = edge_vertices[facet0][0];
    const unsigned int v1 = edge_vertices[facet0][1];
    
    // Compute scale factor (length of edge scaled by length of reference interval)
    const double dx0 = x0[v1][0] - x0[v0][0];
    const double dx1 = x0[v1][1] - x0[v0][1];
    const double det = std::sqrt(dx0*dx0 + dx1*dx1);
    
    const bool direction = dx1*(x0[facet0][0] - x0[v0][0]) - dx0*(x0[facet0][1] - x0[v0][1]) < 0;
    
    // Compute facet normals from the facet scale factor constants
    const double n00 = direction ? dx1 / det : -dx1 / det;
    const double n01 = direction ? -dx0 / det : dx0 / det;
    
    // Compute facet normals from the facet scale factor constants
    const double n10 = !direction ? dx1 / det : -dx1 / det;
    const double n11 = !direction ? -dx0 / det : dx0 / det;
    
    
    // Array of quadrature weights
    static const double W1 = 1;
    // Quadrature points on the UFC reference element: (0.5)
    
    // Value of basis functions at quadrature points.
    static const double FE0_f0[1][1] = \
    {{1}};
    
    static const double FE1_f0_C0[1][7] = \
    {{1, -0.5, -0.5, 1, -0.5, 1, 0}};
    
    static const double FE1_f0_C1[1][7] = \
    {{-0.5, 1, 0.5, -1, 0.5, -1, 0}};
    
    static const double FE1_f1_C0[1][7] = \
    {{0, 0, 0.5, 0.5, 0, 0, 0}};
    
    static const double FE1_f1_C1[1][7] = \
    {{-0.5, 1, 0.5, -1, -1, 0.5, 0}};
    
    static const double FE1_f2_C0[1][7] = \
    {{1, -0.5, 1, -0.5, -0.5, 1, 0}};
    
    static const double FE1_f2_C1[1][7] = \
    {{0, 0, 0, 0, -0.5, -0.5, 0}};
    
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', False), ('ignore zero tables', False), ('non zero columns', False), ('remove zero terms', False), ('ignore ones', False)
    switch ( facet0 )
    {
    case 0:
      switch ( facet1 )
      {
      case 0:
        {
        // Total number of operations to compute element tensor (from this point): 191
        
        // Loop quadrature points for integral
        // Number of operations to compute element tensor for following IP loop = 191
        // Only 1 integration point, omitting IP loop.
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        double F4 = 0;
        double F5 = 0;
        double F6 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE0_f0[0][r]*w[4][r];
          F3 += FE0_f0[0][r]*w[2][r];
          F6 += FE0_f0[0][r]*w[2][r + 1];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 56
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE1_f0_C0[0][r]*w[0][r];
          F2 += FE1_f0_C1[0][r]*w[0][r];
          F4 += FE1_f0_C0[0][r]*w[0][r + 7];
          F5 += FE1_f0_C1[0][r]*w[0][r + 7];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 129
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 65
          A[(j + 1)] += FE0_f0[0][j]*-1*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
          // Number of operations to compute entry: 64
          A[j] += FE0_f0[0][j]*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
        }// end loop over 'j'
        }
        break;
      case 1:
        {
        // Total number of operations to compute element tensor (from this point): 191
        
        // Loop quadrature points for integral
        // Number of operations to compute element tensor for following IP loop = 191
        // Only 1 integration point, omitting IP loop.
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        double F4 = 0;
        double F5 = 0;
        double F6 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE0_f0[0][r]*w[4][r];
          F3 += FE0_f0[0][r]*w[2][r];
          F6 += FE0_f0[0][r]*w[2][r + 1];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 56
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE1_f0_C0[0][r]*w[0][r];
          F2 += FE1_f0_C1[0][r]*w[0][r];
          F4 += FE1_f1_C0[0][r]*w[0][r + 7];
          F5 += FE1_f1_C1[0][r]*w[0][r + 7];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 129
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 65
          A[(j + 1)] += FE0_f0[0][j]*-1*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
          // Number of operations to compute entry: 64
          A[j] += FE0_f0[0][j]*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
        }// end loop over 'j'
        }
        break;
      case 2:
        {
        // Total number of operations to compute element tensor (from this point): 191
        
        // Loop quadrature points for integral
        // Number of operations to compute element tensor for following IP loop = 191
        // Only 1 integration point, omitting IP loop.
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        double F4 = 0;
        double F5 = 0;
        double F6 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE0_f0[0][r]*w[4][r];
          F3 += FE0_f0[0][r]*w[2][r];
          F6 += FE0_f0[0][r]*w[2][r + 1];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 56
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE1_f0_C0[0][r]*w[0][r];
          F2 += FE1_f0_C1[0][r]*w[0][r];
          F4 += FE1_f2_C0[0][r]*w[0][r + 7];
          F5 += FE1_f2_C1[0][r]*w[0][r + 7];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 129
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 65
          A[(j + 1)] += FE0_f0[0][j]*-1*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
          // Number of operations to compute entry: 64
          A[j] += FE0_f0[0][j]*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
        }// end loop over 'j'
        }
        break;
      }
      break;
    case 1:
      switch ( facet1 )
      {
      case 0:
        {
        // Total number of operations to compute element tensor (from this point): 191
        
        // Loop quadrature points for integral
        // Number of operations to compute element tensor for following IP loop = 191
        // Only 1 integration point, omitting IP loop.
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        double F4 = 0;
        double F5 = 0;
        double F6 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE0_f0[0][r]*w[4][r];
          F3 += FE0_f0[0][r]*w[2][r];
          F6 += FE0_f0[0][r]*w[2][r + 1];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 56
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE1_f1_C0[0][r]*w[0][r];
          F2 += FE1_f1_C1[0][r]*w[0][r];
          F4 += FE1_f0_C0[0][r]*w[0][r + 7];
          F5 += FE1_f0_C1[0][r]*w[0][r + 7];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 129
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 65
          A[(j + 1)] += FE0_f0[0][j]*-1*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
          // Number of operations to compute entry: 64
          A[j] += FE0_f0[0][j]*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
        }// end loop over 'j'
        }
        break;
      case 1:
        {
        // Total number of operations to compute element tensor (from this point): 191
        
        // Loop quadrature points for integral
        // Number of operations to compute element tensor for following IP loop = 191
        // Only 1 integration point, omitting IP loop.
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        double F4 = 0;
        double F5 = 0;
        double F6 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE0_f0[0][r]*w[4][r];
          F3 += FE0_f0[0][r]*w[2][r];
          F6 += FE0_f0[0][r]*w[2][r + 1];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 56
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE1_f1_C0[0][r]*w[0][r];
          F2 += FE1_f1_C1[0][r]*w[0][r];
          F4 += FE1_f1_C0[0][r]*w[0][r + 7];
          F5 += FE1_f1_C1[0][r]*w[0][r + 7];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 129
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 65
          A[(j + 1)] += FE0_f0[0][j]*-1*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
          // Number of operations to compute entry: 64
          A[j] += FE0_f0[0][j]*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
        }// end loop over 'j'
        }
        break;
      case 2:
        {
        // Total number of operations to compute element tensor (from this point): 191
        
        // Loop quadrature points for integral
        // Number of operations to compute element tensor for following IP loop = 191
        // Only 1 integration point, omitting IP loop.
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        double F4 = 0;
        double F5 = 0;
        double F6 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE0_f0[0][r]*w[4][r];
          F3 += FE0_f0[0][r]*w[2][r];
          F6 += FE0_f0[0][r]*w[2][r + 1];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 56
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE1_f1_C0[0][r]*w[0][r];
          F2 += FE1_f1_C1[0][r]*w[0][r];
          F4 += FE1_f2_C0[0][r]*w[0][r + 7];
          F5 += FE1_f2_C1[0][r]*w[0][r + 7];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 129
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 65
          A[(j + 1)] += FE0_f0[0][j]*-1*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
          // Number of operations to compute entry: 64
          A[j] += FE0_f0[0][j]*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
        }// end loop over 'j'
        }
        break;
      }
      break;
    case 2:
      switch ( facet1 )
      {
      case 0:
        {
        // Total number of operations to compute element tensor (from this point): 191
        
        // Loop quadrature points for integral
        // Number of operations to compute element tensor for following IP loop = 191
        // Only 1 integration point, omitting IP loop.
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        double F4 = 0;
        double F5 = 0;
        double F6 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE0_f0[0][r]*w[4][r];
          F3 += FE0_f0[0][r]*w[2][r];
          F6 += FE0_f0[0][r]*w[2][r + 1];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 56
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE1_f2_C0[0][r]*w[0][r];
          F2 += FE1_f2_C1[0][r]*w[0][r];
          F4 += FE1_f0_C0[0][r]*w[0][r + 7];
          F5 += FE1_f0_C1[0][r]*w[0][r + 7];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 129
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 65
          A[(j + 1)] += FE0_f0[0][j]*-1*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
          // Number of operations to compute entry: 64
          A[j] += FE0_f0[0][j]*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
        }// end loop over 'j'
        }
        break;
      case 1:
        {
        // Total number of operations to compute element tensor (from this point): 191
        
        // Loop quadrature points for integral
        // Number of operations to compute element tensor for following IP loop = 191
        // Only 1 integration point, omitting IP loop.
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        double F4 = 0;
        double F5 = 0;
        double F6 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE0_f0[0][r]*w[4][r];
          F3 += FE0_f0[0][r]*w[2][r];
          F6 += FE0_f0[0][r]*w[2][r + 1];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 56
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE1_f2_C0[0][r]*w[0][r];
          F2 += FE1_f2_C1[0][r]*w[0][r];
          F4 += FE1_f1_C0[0][r]*w[0][r + 7];
          F5 += FE1_f1_C1[0][r]*w[0][r + 7];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 129
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 65
          A[(j + 1)] += FE0_f0[0][j]*-1*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
          // Number of operations to compute entry: 64
          A[j] += FE0_f0[0][j]*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
        }// end loop over 'j'
        }
        break;
      case 2:
        {
        // Total number of operations to compute element tensor (from this point): 191
        
        // Loop quadrature points for integral
        // Number of operations to compute element tensor for following IP loop = 191
        // Only 1 integration point, omitting IP loop.
        
        // Function declarations
        double F0 = 0;
        double F1 = 0;
        double F2 = 0;
        double F3 = 0;
        double F4 = 0;
        double F5 = 0;
        double F6 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 1; r++)
        {
          F0 += FE0_f0[0][r]*w[4][r];
          F3 += FE0_f0[0][r]*w[2][r];
          F6 += FE0_f0[0][r]*w[2][r + 1];
        }// end loop over 'r'
        
        // Total number of operations to compute function values = 56
        for (unsigned int r = 0; r < 7; r++)
        {
          F1 += FE1_f2_C0[0][r]*w[0][r];
          F2 += FE1_f2_C1[0][r]*w[0][r];
          F4 += FE1_f2_C0[0][r]*w[0][r + 7];
          F5 += FE1_f2_C1[0][r]*w[0][r + 7];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 129
        for (unsigned int j = 0; j < 1; j++)
        {
          // Number of operations to compute entry: 65
          A[(j + 1)] += FE0_f0[0][j]*-1*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
          // Number of operations to compute entry: 64
          A[j] += FE0_f0[0][j]*(-1*(std::abs((((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10)) + (((1.0/detJ1)*J1_10*F4 + (1.0/detJ1)*J1_11*F5)*n11 + ((1.0/detJ1)*J1_00*F4 + (1.0/detJ1)*J1_01*F5)*n10))/(2)*F6 + (std::abs((((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00)) + (((1.0/detJ0)*J0_10*F1 + (1.0/detJ0)*J0_11*F2)*n01 + ((1.0/detJ0)*J0_00*F1 + (1.0/detJ0)*J0_01*F2)*n00))/(2)*F3)*F0*-1*W1*det;
        }// end loop over 'j'
        }
        break;
      }
      break;
    }
}

/// Constructor
saturationequation_1_interior_facet_integral_0::saturationequation_1_interior_facet_integral_0() : ufc::interior_facet_integral()
{
    // Do nothing
}

/// Destructor
saturationequation_1_interior_facet_integral_0::~saturationequation_1_interior_facet_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local interior facet
void saturationequation_1_interior_facet_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c0,
                                    const ufc::cell& c1,
                                    unsigned int facet0,
                                    unsigned int facet1) const
{
    // Reset values of the element tensor block
    for (unsigned int j = 0; j < 2; j++)
      A[j] = 0;
    
    // Add all contributions to element tensor
    integral_0_quadrature.tabulate_tensor(A, w, c0, c1, facet0, facet1);
}

/// Constructor
saturationequation_form_1::saturationequation_form_1() : ufc::form()
{
    // Do nothing
}

/// Destructor
saturationequation_form_1::~saturationequation_form_1()
{
    // Do nothing
}

/// Return a string identifying the form
const char* saturationequation_form_1::signature() const
{
    return "Form([Integral(Sum(Product(BasisFunction(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 0), Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 1)), Product(IndexSum(Product(Indexed(ComponentTensor(Product(Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 2), Indexed(ListTensor(Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(0),), {Index(0): 2}))), MultiIndex((Index(0),), {Index(0): 2})), MultiIndex((Index(1),), {Index(1): 2})), Indexed(ComponentTensor(SpatialDerivative(BasisFunction(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 0), MultiIndex((Index(2),), {Index(2): 2})), MultiIndex((Index(2),), {Index(2): 2})), MultiIndex((Index(1),), {Index(1): 2}))), MultiIndex((Index(1),), {Index(1): 2})), PositiveRestricted(Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 4)))), Measure('cell', 0, None)), Integral(Product(IntValue(-1, (), (), {}), Product(PositiveRestricted(Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 4)), Product(Sum(PositiveRestricted(BasisFunction(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 0)), Product(IntValue(-1, (), (), {}), NegativeRestricted(BasisFunction(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 0)))), Sum(Product(IntValue(-1, (), (), {}), Product(NegativeRestricted(Division(Sum(Abs(IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(3),), {Index(3): 2})), Indexed(ListTensor(Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(3),), {Index(3): 2}))), MultiIndex((Index(3),), {Index(3): 2}))), IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(4),), {Index(4): 2})), Indexed(ListTensor(Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(4),), {Index(4): 2}))), MultiIndex((Index(4),), {Index(4): 2}))), FloatValue(2.0, (), (), {}))), NegativeRestricted(Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 2)))), Product(PositiveRestricted(Division(Sum(Abs(IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(5),), {Index(5): 2})), Indexed(ListTensor(Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(5),), {Index(5): 2}))), MultiIndex((Index(5),), {Index(5): 2}))), IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(6),), {Index(6): 2})), Indexed(ListTensor(Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(6),), {Index(6): 2}))), MultiIndex((Index(6),), {Index(6): 2}))), FloatValue(2.0, (), (), {}))), PositiveRestricted(Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 2))))))), Measure('interior_facet', 0, None)), Integral(Product(IntValue(-1, (), (), {}), Product(Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 4), Product(BasisFunction(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 0), Product(Division(Sum(Abs(IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(7),), {Index(7): 2})), Indexed(ListTensor(Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(7),), {Index(7): 2}))), MultiIndex((Index(7),), {Index(7): 2}))), IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(8),), {Index(8): 2})), Indexed(ListTensor(Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(8),), {Index(8): 2}))), MultiIndex((Index(8),), {Index(8): 2}))), FloatValue(2.0, (), (), {})), Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 2))))), Measure('exterior_facet', 0, None)), Integral(Product(IntValue(-1, (), (), {}), Product(Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 4), Product(BasisFunction(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0), 0), Product(Division(Sum(IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(9),), {Index(9): 2})), Indexed(ListTensor(Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(9),), {Index(9): 2}))), MultiIndex((Index(9),), {Index(9): 2})), Product(IntValue(-1, (), (), {}), Abs(IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(10),), {Index(10): 2})), Indexed(ListTensor(Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(Function(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 1), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 0)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(10),), {Index(10): 2}))), MultiIndex((Index(10),), {Index(10): 2}))))), FloatValue(2.0, (), (), {})), Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1), 3))))), Measure('exterior_facet', 1, None))])";
}

/// Return the rank of the global tensor (r)
unsigned int saturationequation_form_1::rank() const
{
    return 1;
}

/// Return the number of coefficients (n)
unsigned int saturationequation_form_1::num_coefficients() const
{
    return 5;
}

/// Return the number of cell integrals
unsigned int saturationequation_form_1::num_cell_integrals() const
{
    return 1;
}

/// Return the number of exterior facet integrals
unsigned int saturationequation_form_1::num_exterior_facet_integrals() const
{
    return 2;
}

/// Return the number of interior facet integrals
unsigned int saturationequation_form_1::num_interior_facet_integrals() const
{
    return 1;
}

/// Create a new finite element for argument function i
ufc::finite_element* saturationequation_form_1::create_finite_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new saturationequation_1_finite_element_0();
      break;
    case 1:
      return new saturationequation_1_finite_element_1();
      break;
    case 2:
      return new saturationequation_1_finite_element_2();
      break;
    case 3:
      return new saturationequation_1_finite_element_3();
      break;
    case 4:
      return new saturationequation_1_finite_element_4();
      break;
    case 5:
      return new saturationequation_1_finite_element_5();
      break;
    }
    return 0;
}

/// Create a new dof map for argument function i
ufc::dof_map* saturationequation_form_1::create_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new saturationequation_1_dof_map_0();
      break;
    case 1:
      return new saturationequation_1_dof_map_1();
      break;
    case 2:
      return new saturationequation_1_dof_map_2();
      break;
    case 3:
      return new saturationequation_1_dof_map_3();
      break;
    case 4:
      return new saturationequation_1_dof_map_4();
      break;
    case 5:
      return new saturationequation_1_dof_map_5();
      break;
    }
    return 0;
}

/// Create a new cell integral on sub domain i
ufc::cell_integral* saturationequation_form_1::create_cell_integral(unsigned int i) const
{
    return new saturationequation_1_cell_integral_0();
}

/// Create a new exterior facet integral on sub domain i
ufc::exterior_facet_integral* saturationequation_form_1::create_exterior_facet_integral(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new saturationequation_1_exterior_facet_integral_0();
      break;
    case 1:
      return new saturationequation_1_exterior_facet_integral_1();
      break;
    }
    return 0;
}

/// Create a new interior facet integral on sub domain i
ufc::interior_facet_integral* saturationequation_form_1::create_interior_facet_integral(unsigned int i) const
{
    return new saturationequation_1_interior_facet_integral_0();
}


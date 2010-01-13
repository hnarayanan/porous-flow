#include "DarcyFlow.h"

/// Constructor
darcyflow_0_finite_element_0_0::darcyflow_0_finite_element_0_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_0_finite_element_0_0::~darcyflow_0_finite_element_0_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_0_finite_element_0_0::signature() const
{
    return "FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2)";
}

/// Return the cell shape
ufc::shape darcyflow_0_finite_element_0_0::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_0_finite_element_0_0::space_dimension() const
{
    return 12;
}

/// Return the rank of the value space
unsigned int darcyflow_0_finite_element_0_0::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_0_finite_element_0_0::value_dimension(unsigned int i) const
{
    return 2;
}

/// Evaluate basis function i at given point in cell
void darcyflow_0_finite_element_0_0::evaluate_basis(unsigned int i,
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
    const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
    const double psitilde_bs_1_0 = 1;
    const double psitilde_bs_1_1 = 2.5*y + 1.5;
    const double psitilde_bs_2_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
    const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
    const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
    // Table(s) of coefficients
    static const double coefficients0[12][6] = \
    {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
    {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
    {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
    {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
    {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
    {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
    {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
    {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
    {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
    {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
    {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
    {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
    static const double coefficients1[12][6] = \
    {{0, 0, 0.2, 0, 0, 0.163299316185545},
    {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
    {0, 0, 0.6, 0, 0, 0.489897948556636},
    {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
    {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
    {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
    {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
    {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
    {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
    {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
    {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
    {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff0_3 = coefficients0[dof][3];
    const double coeff0_4 = coefficients0[dof][4];
    const double coeff0_5 = coefficients0[dof][5];
    const double coeff1_0 = coefficients1[dof][0];
    const double coeff1_1 = coefficients1[dof][1];
    const double coeff1_2 = coefficients1[dof][2];
    const double coeff1_3 = coefficients1[dof][3];
    const double coeff1_4 = coefficients1[dof][4];
    const double coeff1_5 = coefficients1[dof][5];
    
    // Compute value(s)
    const double tmp0_0 = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3 + coeff0_4*basisvalue4 + coeff0_5*basisvalue5;
    const double tmp0_1 = coeff1_0*basisvalue0 + coeff1_1*basisvalue1 + coeff1_2*basisvalue2 + coeff1_3*basisvalue3 + coeff1_4*basisvalue4 + coeff1_5*basisvalue5;
    // Using contravariant Piola transform to map values back to the physical element
    values[0] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
    values[1] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
}

/// Evaluate all basis functions at given point in cell
void darcyflow_0_finite_element_0_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_0_finite_element_0_0::evaluate_basis_derivatives(unsigned int i,
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
    const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
    const double psitilde_bs_1_0 = 1;
    const double psitilde_bs_1_1 = 2.5*y + 1.5;
    const double psitilde_bs_2_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
    const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
    const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
    // Table(s) of coefficients
    static const double coefficients0[12][6] = \
    {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
    {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
    {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
    {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
    {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
    {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
    {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
    {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
    {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
    {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
    {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
    {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
    static const double coefficients1[12][6] = \
    {{0, 0, 0.2, 0, 0, 0.163299316185545},
    {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
    {0, 0, 0.6, 0, 0, 0.489897948556636},
    {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
    {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
    {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
    {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
    {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
    {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
    {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
    {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
    {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[6][6] = \
    {{0, 0, 0, 0, 0, 0},
    {4.89897948556635, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 9.48683298050514, 0, 0, 0, 0},
    {4, 0, 7.07106781186548, 0, 0, 0},
    {0, 0, 0, 0, 0, 0}};
    
    static const double dmats1[6][6] = \
    {{0, 0, 0, 0, 0, 0},
    {2.44948974278318, 0, 0, 0, 0, 0},
    {4.24264068711928, 0, 0, 0, 0, 0},
    {2.58198889747161, 4.74341649025257, -0.912870929175277, 0, 0, 0},
    {2, 6.12372435695795, 3.53553390593274, 0, 0, 0},
    {-2.3094010767585, 0, 8.16496580927726, 0, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [2*num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    double coeff0_3 = 0;
    double coeff0_4 = 0;
    double coeff0_5 = 0;
    double coeff1_0 = 0;
    double coeff1_1 = 0;
    double coeff1_2 = 0;
    double coeff1_3 = 0;
    double coeff1_4 = 0;
    double coeff1_5 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff0_3 = 0;
    double new_coeff0_4 = 0;
    double new_coeff0_5 = 0;
    double new_coeff1_0 = 0;
    double new_coeff1_1 = 0;
    double new_coeff1_2 = 0;
    double new_coeff1_3 = 0;
    double new_coeff1_4 = 0;
    double new_coeff1_5 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
      new_coeff0_3 = coefficients0[dof][3];
      new_coeff0_4 = coefficients0[dof][4];
      new_coeff0_5 = coefficients0[dof][5];
      new_coeff1_0 = coefficients1[dof][0];
      new_coeff1_1 = coefficients1[dof][1];
      new_coeff1_2 = coefficients1[dof][2];
      new_coeff1_3 = coefficients1[dof][3];
      new_coeff1_4 = coefficients1[dof][4];
      new_coeff1_5 = coefficients1[dof][5];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
        coeff0_3 = new_coeff0_3;
        coeff0_4 = new_coeff0_4;
        coeff0_5 = new_coeff0_5;
        coeff1_0 = new_coeff1_0;
        coeff1_1 = new_coeff1_1;
        coeff1_2 = new_coeff1_2;
        coeff1_3 = new_coeff1_3;
        coeff1_4 = new_coeff1_4;
        coeff1_5 = new_coeff1_5;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0] + coeff0_4*dmats0[4][0] + coeff0_5*dmats0[5][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1] + coeff0_4*dmats0[4][1] + coeff0_5*dmats0[5][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2] + coeff0_4*dmats0[4][2] + coeff0_5*dmats0[5][2];
          new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3] + coeff0_4*dmats0[4][3] + coeff0_5*dmats0[5][3];
          new_coeff0_4 = coeff0_0*dmats0[0][4] + coeff0_1*dmats0[1][4] + coeff0_2*dmats0[2][4] + coeff0_3*dmats0[3][4] + coeff0_4*dmats0[4][4] + coeff0_5*dmats0[5][4];
          new_coeff0_5 = coeff0_0*dmats0[0][5] + coeff0_1*dmats0[1][5] + coeff0_2*dmats0[2][5] + coeff0_3*dmats0[3][5] + coeff0_4*dmats0[4][5] + coeff0_5*dmats0[5][5];
          new_coeff1_0 = coeff1_0*dmats0[0][0] + coeff1_1*dmats0[1][0] + coeff1_2*dmats0[2][0] + coeff1_3*dmats0[3][0] + coeff1_4*dmats0[4][0] + coeff1_5*dmats0[5][0];
          new_coeff1_1 = coeff1_0*dmats0[0][1] + coeff1_1*dmats0[1][1] + coeff1_2*dmats0[2][1] + coeff1_3*dmats0[3][1] + coeff1_4*dmats0[4][1] + coeff1_5*dmats0[5][1];
          new_coeff1_2 = coeff1_0*dmats0[0][2] + coeff1_1*dmats0[1][2] + coeff1_2*dmats0[2][2] + coeff1_3*dmats0[3][2] + coeff1_4*dmats0[4][2] + coeff1_5*dmats0[5][2];
          new_coeff1_3 = coeff1_0*dmats0[0][3] + coeff1_1*dmats0[1][3] + coeff1_2*dmats0[2][3] + coeff1_3*dmats0[3][3] + coeff1_4*dmats0[4][3] + coeff1_5*dmats0[5][3];
          new_coeff1_4 = coeff1_0*dmats0[0][4] + coeff1_1*dmats0[1][4] + coeff1_2*dmats0[2][4] + coeff1_3*dmats0[3][4] + coeff1_4*dmats0[4][4] + coeff1_5*dmats0[5][4];
          new_coeff1_5 = coeff1_0*dmats0[0][5] + coeff1_1*dmats0[1][5] + coeff1_2*dmats0[2][5] + coeff1_3*dmats0[3][5] + coeff1_4*dmats0[4][5] + coeff1_5*dmats0[5][5];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0] + coeff0_4*dmats1[4][0] + coeff0_5*dmats1[5][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1] + coeff0_4*dmats1[4][1] + coeff0_5*dmats1[5][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2] + coeff0_4*dmats1[4][2] + coeff0_5*dmats1[5][2];
          new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3] + coeff0_4*dmats1[4][3] + coeff0_5*dmats1[5][3];
          new_coeff0_4 = coeff0_0*dmats1[0][4] + coeff0_1*dmats1[1][4] + coeff0_2*dmats1[2][4] + coeff0_3*dmats1[3][4] + coeff0_4*dmats1[4][4] + coeff0_5*dmats1[5][4];
          new_coeff0_5 = coeff0_0*dmats1[0][5] + coeff0_1*dmats1[1][5] + coeff0_2*dmats1[2][5] + coeff0_3*dmats1[3][5] + coeff0_4*dmats1[4][5] + coeff0_5*dmats1[5][5];
          new_coeff1_0 = coeff1_0*dmats1[0][0] + coeff1_1*dmats1[1][0] + coeff1_2*dmats1[2][0] + coeff1_3*dmats1[3][0] + coeff1_4*dmats1[4][0] + coeff1_5*dmats1[5][0];
          new_coeff1_1 = coeff1_0*dmats1[0][1] + coeff1_1*dmats1[1][1] + coeff1_2*dmats1[2][1] + coeff1_3*dmats1[3][1] + coeff1_4*dmats1[4][1] + coeff1_5*dmats1[5][1];
          new_coeff1_2 = coeff1_0*dmats1[0][2] + coeff1_1*dmats1[1][2] + coeff1_2*dmats1[2][2] + coeff1_3*dmats1[3][2] + coeff1_4*dmats1[4][2] + coeff1_5*dmats1[5][2];
          new_coeff1_3 = coeff1_0*dmats1[0][3] + coeff1_1*dmats1[1][3] + coeff1_2*dmats1[2][3] + coeff1_3*dmats1[3][3] + coeff1_4*dmats1[4][3] + coeff1_5*dmats1[5][3];
          new_coeff1_4 = coeff1_0*dmats1[0][4] + coeff1_1*dmats1[1][4] + coeff1_2*dmats1[2][4] + coeff1_3*dmats1[3][4] + coeff1_4*dmats1[4][4] + coeff1_5*dmats1[5][4];
          new_coeff1_5 = coeff1_0*dmats1[0][5] + coeff1_1*dmats1[1][5] + coeff1_2*dmats1[2][5] + coeff1_3*dmats1[3][5] + coeff1_4*dmats1[4][5] + coeff1_5*dmats1[5][5];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      // Correct values by the contravariant Piola transform
      const double tmp0_0 = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3 + new_coeff0_4*basisvalue4 + new_coeff0_5*basisvalue5;
      const double tmp0_1 = new_coeff1_0*basisvalue0 + new_coeff1_1*basisvalue1 + new_coeff1_2*basisvalue2 + new_coeff1_3*basisvalue3 + new_coeff1_4*basisvalue4 + new_coeff1_5*basisvalue5;
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
void darcyflow_0_finite_element_0_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_0_finite_element_0_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[12][9][2] = {{{0.75, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.75, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}};
    static const double W[12][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}};
    static const double D[12][9][2] = {{{1, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, -1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, -1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, -1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.866025403784438, 0.433012701892219}, {0.866025403784439, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784439, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}}, {{0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}}, {{-0.17542966950853, 0.148329095604429}, {-0.0180753884489578, 0.444221856552778}, {0.0340120331760423, 0.307261606416198}, {-0.396781833524282, 3.42776476771047e-16}, {0.0799728864560174, -2.52114035295885e-16}, {0.171602319815636, -3.84645013221013e-15}, {-0.0271005739041022, -0.14832909560443}, {0.426146468103823, -0.44422185655278}, {0.341273639592246, -0.307261606416206}}};
    
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
    static const unsigned int ns[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 9, 9};
    for (unsigned int j = 0; j < ns[i]; j++) {
      // Evaluate basis functions for affine mapping
      const double w0 = 1.0 - X[i][j][0] - X[i][j][1];
      const double w1 = X[i][j][0];
      const double w2 = X[i][j][1];
      
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
        result += values[k]*D[i][j][k];
      // Multiply by weights
      result *= W[i][j];
    
    } // End for
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void darcyflow_0_finite_element_0_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_0_finite_element_0_0::interpolate_vertex_values(double* vertex_values,
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
    vertex_values[0] = (1.0/detJ)*(dof_values[3]*3*J_00 + dof_values[4]*(-3*J_00) + dof_values[5]*J_00 + dof_values[6]*(-3*J_01) + dof_values[7]*3*J_01 + dof_values[8]*J_01);
    vertex_values[2] = (1.0/detJ)*(dof_values[0]*3*J_00 + dof_values[1]*(-3*J_00) + dof_values[2]*J_00 + dof_values[6]*(J_00 - J_01) + dof_values[7]*(-3*J_00 + 3*J_01) + dof_values[8]*(3*J_00 - 3*J_01));
    vertex_values[4] = (1.0/detJ)*(dof_values[0]*J_01 + dof_values[1]*(-3*J_01) + dof_values[2]*3*J_01 + dof_values[3]*(J_00 - J_01) + dof_values[4]*(-3*J_00 + 3*J_01) + dof_values[5]*(3*J_00 - 3*J_01));
    vertex_values[1] = (1.0/detJ)*(dof_values[3]*3*J_10 + dof_values[4]*(-3*J_10) + dof_values[5]*J_10 + dof_values[6]*(-3*J_11) + dof_values[7]*3*J_11 + dof_values[8]*J_11);
    vertex_values[3] = (1.0/detJ)*(dof_values[0]*3*J_10 + dof_values[1]*(-3*J_10) + dof_values[2]*J_10 + dof_values[6]*(J_10 - J_11) + dof_values[7]*(-3*J_10 + 3*J_11) + dof_values[8]*(3*J_10 - 3*J_11));
    vertex_values[5] = (1.0/detJ)*(dof_values[0]*J_11 + dof_values[1]*(-3*J_11) + dof_values[2]*3*J_11 + dof_values[3]*(J_10 - J_11) + dof_values[4]*(-3*J_10 + 3*J_11) + dof_values[5]*(3*J_10 - 3*J_11));
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_0_finite_element_0_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_0_finite_element_0_0::create_sub_element(unsigned int i) const
{
    return new darcyflow_0_finite_element_0_0();
}


/// Constructor
darcyflow_0_finite_element_0_1::darcyflow_0_finite_element_0_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_0_finite_element_0_1::~darcyflow_0_finite_element_0_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_0_finite_element_0_1::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return the cell shape
ufc::shape darcyflow_0_finite_element_0_1::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_0_finite_element_0_1::space_dimension() const
{
    return 3;
}

/// Return the rank of the value space
unsigned int darcyflow_0_finite_element_0_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_0_finite_element_0_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void darcyflow_0_finite_element_0_1::evaluate_basis(unsigned int i,
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
void darcyflow_0_finite_element_0_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_0_finite_element_0_1::evaluate_basis_derivatives(unsigned int i,
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
void darcyflow_0_finite_element_0_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_0_finite_element_0_1::evaluate_dof(unsigned int i,
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
void darcyflow_0_finite_element_0_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_0_finite_element_0_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_0_finite_element_0_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_0_finite_element_0_1::create_sub_element(unsigned int i) const
{
    return new darcyflow_0_finite_element_0_1();
}


/// Constructor
darcyflow_0_finite_element_0::darcyflow_0_finite_element_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_0_finite_element_0::~darcyflow_0_finite_element_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_0_finite_element_0::signature() const
{
    return "MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) })";
}

/// Return the cell shape
ufc::shape darcyflow_0_finite_element_0::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_0_finite_element_0::space_dimension() const
{
    return 15;
}

/// Return the rank of the value space
unsigned int darcyflow_0_finite_element_0::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_0_finite_element_0::value_dimension(unsigned int i) const
{
    return 3;
}

/// Evaluate basis function i at given point in cell
void darcyflow_0_finite_element_0::evaluate_basis(unsigned int i,
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
    
    if (0 <= i && i <= 11)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
      const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
      const double psitilde_bs_1_0 = 1;
      const double psitilde_bs_1_1 = 2.5*y + 1.5;
      const double psitilde_bs_2_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
      const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
      const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
      const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
      const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
      const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
      // Table(s) of coefficients
      static const double coefficients0[12][6] =   \
      {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
      {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
      {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
      {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
      {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
      {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
      {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
      {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
      {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
      {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
      {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
      {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
      static const double coefficients1[12][6] =   \
      {{0, 0, 0.2, 0, 0, 0.163299316185545},
      {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
      {0, 0, 0.6, 0, 0, 0.489897948556636},
      {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
      {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
      {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
      {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
      {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
      {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
      {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
      {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
      {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
      const double coeff0_4 =   coefficients0[dof][4];
      const double coeff0_5 =   coefficients0[dof][5];
      const double coeff1_0 =   coefficients1[dof][0];
      const double coeff1_1 =   coefficients1[dof][1];
      const double coeff1_2 =   coefficients1[dof][2];
      const double coeff1_3 =   coefficients1[dof][3];
      const double coeff1_4 =   coefficients1[dof][4];
      const double coeff1_5 =   coefficients1[dof][5];
    
      // Compute value(s)
      const double tmp0_0 = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3 + coeff0_4*basisvalue4 + coeff0_5*basisvalue5;
      const double tmp0_1 = coeff1_0*basisvalue0 + coeff1_1*basisvalue1 + coeff1_2*basisvalue2 + coeff1_3*basisvalue3 + coeff1_4*basisvalue4 + coeff1_5*basisvalue5;
      // Using contravariant Piola transform to map values back to the physical element
      values[0] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
      values[1] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
    }
    
    if (12 <= i && i <= 14)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 12;
    
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
      static const double coefficients0[3][3] =   \
      {{0.471404520791032, -0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0, 0.333333333333333}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
    
      // Compute value(s)
      values[2] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2;
    }
    
}

/// Evaluate all basis functions at given point in cell
void darcyflow_0_finite_element_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_0_finite_element_0::evaluate_basis_derivatives(unsigned int i,
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
    
    if (0 <= i && i <= 11)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
      const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
      const double psitilde_bs_1_0 = 1;
      const double psitilde_bs_1_1 = 2.5*y + 1.5;
      const double psitilde_bs_2_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
      const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
      const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
      const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
      const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
      const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
      // Table(s) of coefficients
      static const double coefficients0[12][6] =   \
      {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
      {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
      {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
      {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
      {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
      {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
      {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
      {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
      {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
      {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
      {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
      {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
      static const double coefficients1[12][6] =   \
      {{0, 0, 0.2, 0, 0, 0.163299316185545},
      {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
      {0, 0, 0.6, 0, 0, 0.489897948556636},
      {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
      {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
      {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
      {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
      {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
      {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
      {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
      {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
      {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      static const double dmats0[6][6] =   \
      {{0, 0, 0, 0, 0, 0},
      {4.89897948556635, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
      {0, 9.48683298050514, 0, 0, 0, 0},
      {4, 0, 7.07106781186548, 0, 0, 0},
      {0, 0, 0, 0, 0, 0}};
    
      static const double dmats1[6][6] =   \
      {{0, 0, 0, 0, 0, 0},
      {2.44948974278318, 0, 0, 0, 0, 0},
      {4.24264068711928, 0, 0, 0, 0, 0},
      {2.58198889747161, 4.74341649025257, -0.912870929175277, 0, 0, 0},
      {2, 6.12372435695795, 3.53553390593274, 0, 0, 0},
      {-2.3094010767585, 0, 8.16496580927726, 0, 0, 0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [2*num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
      double coeff0_1 = 0;
      double coeff0_2 = 0;
      double coeff0_3 = 0;
      double coeff0_4 = 0;
      double coeff0_5 = 0;
      double coeff1_0 = 0;
      double coeff1_1 = 0;
      double coeff1_2 = 0;
      double coeff1_3 = 0;
      double coeff1_4 = 0;
      double coeff1_5 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
      double new_coeff0_1 = 0;
      double new_coeff0_2 = 0;
      double new_coeff0_3 = 0;
      double new_coeff0_4 = 0;
      double new_coeff0_5 = 0;
      double new_coeff1_0 = 0;
      double new_coeff1_1 = 0;
      double new_coeff1_2 = 0;
      double new_coeff1_3 = 0;
      double new_coeff1_4 = 0;
      double new_coeff1_5 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
        new_coeff0_1 = coefficients0[dof][1];
        new_coeff0_2 = coefficients0[dof][2];
        new_coeff0_3 = coefficients0[dof][3];
        new_coeff0_4 = coefficients0[dof][4];
        new_coeff0_5 = coefficients0[dof][5];
        new_coeff1_0 = coefficients1[dof][0];
        new_coeff1_1 = coefficients1[dof][1];
        new_coeff1_2 = coefficients1[dof][2];
        new_coeff1_3 = coefficients1[dof][3];
        new_coeff1_4 = coefficients1[dof][4];
        new_coeff1_5 = coefficients1[dof][5];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
          coeff0_1 = new_coeff0_1;
          coeff0_2 = new_coeff0_2;
          coeff0_3 = new_coeff0_3;
          coeff0_4 = new_coeff0_4;
          coeff0_5 = new_coeff0_5;
          coeff1_0 = new_coeff1_0;
          coeff1_1 = new_coeff1_1;
          coeff1_2 = new_coeff1_2;
          coeff1_3 = new_coeff1_3;
          coeff1_4 = new_coeff1_4;
          coeff1_5 = new_coeff1_5;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0] + coeff0_4*dmats0[4][0] + coeff0_5*dmats0[5][0];
            new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1] + coeff0_4*dmats0[4][1] + coeff0_5*dmats0[5][1];
            new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2] + coeff0_4*dmats0[4][2] + coeff0_5*dmats0[5][2];
            new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3] + coeff0_4*dmats0[4][3] + coeff0_5*dmats0[5][3];
            new_coeff0_4 = coeff0_0*dmats0[0][4] + coeff0_1*dmats0[1][4] + coeff0_2*dmats0[2][4] + coeff0_3*dmats0[3][4] + coeff0_4*dmats0[4][4] + coeff0_5*dmats0[5][4];
            new_coeff0_5 = coeff0_0*dmats0[0][5] + coeff0_1*dmats0[1][5] + coeff0_2*dmats0[2][5] + coeff0_3*dmats0[3][5] + coeff0_4*dmats0[4][5] + coeff0_5*dmats0[5][5];
            new_coeff1_0 = coeff1_0*dmats0[0][0] + coeff1_1*dmats0[1][0] + coeff1_2*dmats0[2][0] + coeff1_3*dmats0[3][0] + coeff1_4*dmats0[4][0] + coeff1_5*dmats0[5][0];
            new_coeff1_1 = coeff1_0*dmats0[0][1] + coeff1_1*dmats0[1][1] + coeff1_2*dmats0[2][1] + coeff1_3*dmats0[3][1] + coeff1_4*dmats0[4][1] + coeff1_5*dmats0[5][1];
            new_coeff1_2 = coeff1_0*dmats0[0][2] + coeff1_1*dmats0[1][2] + coeff1_2*dmats0[2][2] + coeff1_3*dmats0[3][2] + coeff1_4*dmats0[4][2] + coeff1_5*dmats0[5][2];
            new_coeff1_3 = coeff1_0*dmats0[0][3] + coeff1_1*dmats0[1][3] + coeff1_2*dmats0[2][3] + coeff1_3*dmats0[3][3] + coeff1_4*dmats0[4][3] + coeff1_5*dmats0[5][3];
            new_coeff1_4 = coeff1_0*dmats0[0][4] + coeff1_1*dmats0[1][4] + coeff1_2*dmats0[2][4] + coeff1_3*dmats0[3][4] + coeff1_4*dmats0[4][4] + coeff1_5*dmats0[5][4];
            new_coeff1_5 = coeff1_0*dmats0[0][5] + coeff1_1*dmats0[1][5] + coeff1_2*dmats0[2][5] + coeff1_3*dmats0[3][5] + coeff1_4*dmats0[4][5] + coeff1_5*dmats0[5][5];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0] + coeff0_4*dmats1[4][0] + coeff0_5*dmats1[5][0];
            new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1] + coeff0_4*dmats1[4][1] + coeff0_5*dmats1[5][1];
            new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2] + coeff0_4*dmats1[4][2] + coeff0_5*dmats1[5][2];
            new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3] + coeff0_4*dmats1[4][3] + coeff0_5*dmats1[5][3];
            new_coeff0_4 = coeff0_0*dmats1[0][4] + coeff0_1*dmats1[1][4] + coeff0_2*dmats1[2][4] + coeff0_3*dmats1[3][4] + coeff0_4*dmats1[4][4] + coeff0_5*dmats1[5][4];
            new_coeff0_5 = coeff0_0*dmats1[0][5] + coeff0_1*dmats1[1][5] + coeff0_2*dmats1[2][5] + coeff0_3*dmats1[3][5] + coeff0_4*dmats1[4][5] + coeff0_5*dmats1[5][5];
            new_coeff1_0 = coeff1_0*dmats1[0][0] + coeff1_1*dmats1[1][0] + coeff1_2*dmats1[2][0] + coeff1_3*dmats1[3][0] + coeff1_4*dmats1[4][0] + coeff1_5*dmats1[5][0];
            new_coeff1_1 = coeff1_0*dmats1[0][1] + coeff1_1*dmats1[1][1] + coeff1_2*dmats1[2][1] + coeff1_3*dmats1[3][1] + coeff1_4*dmats1[4][1] + coeff1_5*dmats1[5][1];
            new_coeff1_2 = coeff1_0*dmats1[0][2] + coeff1_1*dmats1[1][2] + coeff1_2*dmats1[2][2] + coeff1_3*dmats1[3][2] + coeff1_4*dmats1[4][2] + coeff1_5*dmats1[5][2];
            new_coeff1_3 = coeff1_0*dmats1[0][3] + coeff1_1*dmats1[1][3] + coeff1_2*dmats1[2][3] + coeff1_3*dmats1[3][3] + coeff1_4*dmats1[4][3] + coeff1_5*dmats1[5][3];
            new_coeff1_4 = coeff1_0*dmats1[0][4] + coeff1_1*dmats1[1][4] + coeff1_2*dmats1[2][4] + coeff1_3*dmats1[3][4] + coeff1_4*dmats1[4][4] + coeff1_5*dmats1[5][4];
            new_coeff1_5 = coeff1_0*dmats1[0][5] + coeff1_1*dmats1[1][5] + coeff1_2*dmats1[2][5] + coeff1_3*dmats1[3][5] + coeff1_4*dmats1[4][5] + coeff1_5*dmats1[5][5];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        // Correct values by the contravariant Piola transform
        const double tmp0_0 = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3 + new_coeff0_4*basisvalue4 + new_coeff0_5*basisvalue5;
        const double tmp0_1 = new_coeff1_0*basisvalue0 + new_coeff1_1*basisvalue1 + new_coeff1_2*basisvalue2 + new_coeff1_3*basisvalue3 + new_coeff1_4*basisvalue4 + new_coeff1_5*basisvalue5;
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
    
    if (12 <= i && i <= 14)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 12;
    
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
      static const double coefficients0[3][3] =   \
      {{0.471404520791032, -0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0, 0.333333333333333}};
    
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
void darcyflow_0_finite_element_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_0_finite_element_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[15][9][2] = {{{0.75, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.75, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}};
    static const double W[15][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}};
    static const double D[15][9][3] = {{{1, 1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0.866025403784438, 0.433012701892219, 0}, {0.866025403784439, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784439, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}}, {{0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}}, {{-0.17542966950853, 0.148329095604429, 0}, {-0.0180753884489578, 0.444221856552778, 0}, {0.0340120331760423, 0.307261606416198, 0}, {-0.396781833524282, 3.42776476771047e-16, 0}, {0.0799728864560174, -2.52114035295885e-16, 0}, {0.171602319815636, -3.84645013221013e-15, 0}, {-0.0271005739041022, -0.14832909560443, 0}, {0.426146468103823, -0.44422185655278, 0}, {0.341273639592246, -0.307261606416206, 0}}, {{0, 0, 1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    
    static const unsigned int mappings[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
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
    static const unsigned int ns[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 9, 9, 1, 1, 1};
    for (unsigned int j = 0; j < ns[i]; j++) {
      // Evaluate basis functions for affine mapping
      const double w0 = 1.0 - X[i][j][0] - X[i][j][1];
      const double w1 = X[i][j][0];
      const double w2 = X[i][j][1];
      
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
        result += values[k]*D[i][j][k];
      // Multiply by weights
      result *= W[i][j];
    
    } // End for
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void darcyflow_0_finite_element_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_0_finite_element_0::interpolate_vertex_values(double* vertex_values,
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
    vertex_values[0] = (1.0/detJ)*(dof_values[3]*3*J_00 + dof_values[4]*(-3*J_00) + dof_values[5]*J_00 + dof_values[6]*(-3*J_01) + dof_values[7]*3*J_01 + dof_values[8]*J_01);
    vertex_values[3] = (1.0/detJ)*(dof_values[0]*3*J_00 + dof_values[1]*(-3*J_00) + dof_values[2]*J_00 + dof_values[6]*(J_00 - J_01) + dof_values[7]*(-3*J_00 + 3*J_01) + dof_values[8]*(3*J_00 - 3*J_01));
    vertex_values[6] = (1.0/detJ)*(dof_values[0]*J_01 + dof_values[1]*(-3*J_01) + dof_values[2]*3*J_01 + dof_values[3]*(J_00 - J_01) + dof_values[4]*(-3*J_00 + 3*J_01) + dof_values[5]*(3*J_00 - 3*J_01));
    vertex_values[1] = (1.0/detJ)*(dof_values[3]*3*J_10 + dof_values[4]*(-3*J_10) + dof_values[5]*J_10 + dof_values[6]*(-3*J_11) + dof_values[7]*3*J_11 + dof_values[8]*J_11);
    vertex_values[4] = (1.0/detJ)*(dof_values[0]*3*J_10 + dof_values[1]*(-3*J_10) + dof_values[2]*J_10 + dof_values[6]*(J_10 - J_11) + dof_values[7]*(-3*J_10 + 3*J_11) + dof_values[8]*(3*J_10 - 3*J_11));
    vertex_values[7] = (1.0/detJ)*(dof_values[0]*J_11 + dof_values[1]*(-3*J_11) + dof_values[2]*3*J_11 + dof_values[3]*(J_10 - J_11) + dof_values[4]*(-3*J_10 + 3*J_11) + dof_values[5]*(3*J_10 - 3*J_11));
    // Evaluate at vertices and use affine mapping
    vertex_values[2] = dof_values[12];
    vertex_values[5] = dof_values[13];
    vertex_values[8] = dof_values[14];
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_0_finite_element_0::num_sub_elements() const
{
    return 2;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_0_finite_element_0::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_0_finite_element_0_0();
      break;
    case 1:
      return new darcyflow_0_finite_element_0_1();
      break;
    }
    return 0;
}


/// Constructor
darcyflow_0_finite_element_1_0::darcyflow_0_finite_element_1_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_0_finite_element_1_0::~darcyflow_0_finite_element_1_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_0_finite_element_1_0::signature() const
{
    return "FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2)";
}

/// Return the cell shape
ufc::shape darcyflow_0_finite_element_1_0::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_0_finite_element_1_0::space_dimension() const
{
    return 12;
}

/// Return the rank of the value space
unsigned int darcyflow_0_finite_element_1_0::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_0_finite_element_1_0::value_dimension(unsigned int i) const
{
    return 2;
}

/// Evaluate basis function i at given point in cell
void darcyflow_0_finite_element_1_0::evaluate_basis(unsigned int i,
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
    const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
    const double psitilde_bs_1_0 = 1;
    const double psitilde_bs_1_1 = 2.5*y + 1.5;
    const double psitilde_bs_2_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
    const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
    const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
    // Table(s) of coefficients
    static const double coefficients0[12][6] = \
    {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
    {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
    {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
    {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
    {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
    {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
    {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
    {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
    {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
    {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
    {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
    {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
    static const double coefficients1[12][6] = \
    {{0, 0, 0.2, 0, 0, 0.163299316185545},
    {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
    {0, 0, 0.6, 0, 0, 0.489897948556636},
    {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
    {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
    {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
    {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
    {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
    {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
    {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
    {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
    {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff0_3 = coefficients0[dof][3];
    const double coeff0_4 = coefficients0[dof][4];
    const double coeff0_5 = coefficients0[dof][5];
    const double coeff1_0 = coefficients1[dof][0];
    const double coeff1_1 = coefficients1[dof][1];
    const double coeff1_2 = coefficients1[dof][2];
    const double coeff1_3 = coefficients1[dof][3];
    const double coeff1_4 = coefficients1[dof][4];
    const double coeff1_5 = coefficients1[dof][5];
    
    // Compute value(s)
    const double tmp0_0 = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3 + coeff0_4*basisvalue4 + coeff0_5*basisvalue5;
    const double tmp0_1 = coeff1_0*basisvalue0 + coeff1_1*basisvalue1 + coeff1_2*basisvalue2 + coeff1_3*basisvalue3 + coeff1_4*basisvalue4 + coeff1_5*basisvalue5;
    // Using contravariant Piola transform to map values back to the physical element
    values[0] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
    values[1] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
}

/// Evaluate all basis functions at given point in cell
void darcyflow_0_finite_element_1_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_0_finite_element_1_0::evaluate_basis_derivatives(unsigned int i,
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
    const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
    const double psitilde_bs_1_0 = 1;
    const double psitilde_bs_1_1 = 2.5*y + 1.5;
    const double psitilde_bs_2_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
    const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
    const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
    // Table(s) of coefficients
    static const double coefficients0[12][6] = \
    {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
    {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
    {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
    {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
    {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
    {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
    {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
    {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
    {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
    {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
    {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
    {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
    static const double coefficients1[12][6] = \
    {{0, 0, 0.2, 0, 0, 0.163299316185545},
    {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
    {0, 0, 0.6, 0, 0, 0.489897948556636},
    {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
    {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
    {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
    {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
    {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
    {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
    {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
    {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
    {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[6][6] = \
    {{0, 0, 0, 0, 0, 0},
    {4.89897948556635, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 9.48683298050514, 0, 0, 0, 0},
    {4, 0, 7.07106781186548, 0, 0, 0},
    {0, 0, 0, 0, 0, 0}};
    
    static const double dmats1[6][6] = \
    {{0, 0, 0, 0, 0, 0},
    {2.44948974278318, 0, 0, 0, 0, 0},
    {4.24264068711928, 0, 0, 0, 0, 0},
    {2.58198889747161, 4.74341649025257, -0.912870929175277, 0, 0, 0},
    {2, 6.12372435695795, 3.53553390593274, 0, 0, 0},
    {-2.3094010767585, 0, 8.16496580927726, 0, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [2*num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    double coeff0_3 = 0;
    double coeff0_4 = 0;
    double coeff0_5 = 0;
    double coeff1_0 = 0;
    double coeff1_1 = 0;
    double coeff1_2 = 0;
    double coeff1_3 = 0;
    double coeff1_4 = 0;
    double coeff1_5 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff0_3 = 0;
    double new_coeff0_4 = 0;
    double new_coeff0_5 = 0;
    double new_coeff1_0 = 0;
    double new_coeff1_1 = 0;
    double new_coeff1_2 = 0;
    double new_coeff1_3 = 0;
    double new_coeff1_4 = 0;
    double new_coeff1_5 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
      new_coeff0_3 = coefficients0[dof][3];
      new_coeff0_4 = coefficients0[dof][4];
      new_coeff0_5 = coefficients0[dof][5];
      new_coeff1_0 = coefficients1[dof][0];
      new_coeff1_1 = coefficients1[dof][1];
      new_coeff1_2 = coefficients1[dof][2];
      new_coeff1_3 = coefficients1[dof][3];
      new_coeff1_4 = coefficients1[dof][4];
      new_coeff1_5 = coefficients1[dof][5];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
        coeff0_3 = new_coeff0_3;
        coeff0_4 = new_coeff0_4;
        coeff0_5 = new_coeff0_5;
        coeff1_0 = new_coeff1_0;
        coeff1_1 = new_coeff1_1;
        coeff1_2 = new_coeff1_2;
        coeff1_3 = new_coeff1_3;
        coeff1_4 = new_coeff1_4;
        coeff1_5 = new_coeff1_5;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0] + coeff0_4*dmats0[4][0] + coeff0_5*dmats0[5][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1] + coeff0_4*dmats0[4][1] + coeff0_5*dmats0[5][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2] + coeff0_4*dmats0[4][2] + coeff0_5*dmats0[5][2];
          new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3] + coeff0_4*dmats0[4][3] + coeff0_5*dmats0[5][3];
          new_coeff0_4 = coeff0_0*dmats0[0][4] + coeff0_1*dmats0[1][4] + coeff0_2*dmats0[2][4] + coeff0_3*dmats0[3][4] + coeff0_4*dmats0[4][4] + coeff0_5*dmats0[5][4];
          new_coeff0_5 = coeff0_0*dmats0[0][5] + coeff0_1*dmats0[1][5] + coeff0_2*dmats0[2][5] + coeff0_3*dmats0[3][5] + coeff0_4*dmats0[4][5] + coeff0_5*dmats0[5][5];
          new_coeff1_0 = coeff1_0*dmats0[0][0] + coeff1_1*dmats0[1][0] + coeff1_2*dmats0[2][0] + coeff1_3*dmats0[3][0] + coeff1_4*dmats0[4][0] + coeff1_5*dmats0[5][0];
          new_coeff1_1 = coeff1_0*dmats0[0][1] + coeff1_1*dmats0[1][1] + coeff1_2*dmats0[2][1] + coeff1_3*dmats0[3][1] + coeff1_4*dmats0[4][1] + coeff1_5*dmats0[5][1];
          new_coeff1_2 = coeff1_0*dmats0[0][2] + coeff1_1*dmats0[1][2] + coeff1_2*dmats0[2][2] + coeff1_3*dmats0[3][2] + coeff1_4*dmats0[4][2] + coeff1_5*dmats0[5][2];
          new_coeff1_3 = coeff1_0*dmats0[0][3] + coeff1_1*dmats0[1][3] + coeff1_2*dmats0[2][3] + coeff1_3*dmats0[3][3] + coeff1_4*dmats0[4][3] + coeff1_5*dmats0[5][3];
          new_coeff1_4 = coeff1_0*dmats0[0][4] + coeff1_1*dmats0[1][4] + coeff1_2*dmats0[2][4] + coeff1_3*dmats0[3][4] + coeff1_4*dmats0[4][4] + coeff1_5*dmats0[5][4];
          new_coeff1_5 = coeff1_0*dmats0[0][5] + coeff1_1*dmats0[1][5] + coeff1_2*dmats0[2][5] + coeff1_3*dmats0[3][5] + coeff1_4*dmats0[4][5] + coeff1_5*dmats0[5][5];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0] + coeff0_4*dmats1[4][0] + coeff0_5*dmats1[5][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1] + coeff0_4*dmats1[4][1] + coeff0_5*dmats1[5][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2] + coeff0_4*dmats1[4][2] + coeff0_5*dmats1[5][2];
          new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3] + coeff0_4*dmats1[4][3] + coeff0_5*dmats1[5][3];
          new_coeff0_4 = coeff0_0*dmats1[0][4] + coeff0_1*dmats1[1][4] + coeff0_2*dmats1[2][4] + coeff0_3*dmats1[3][4] + coeff0_4*dmats1[4][4] + coeff0_5*dmats1[5][4];
          new_coeff0_5 = coeff0_0*dmats1[0][5] + coeff0_1*dmats1[1][5] + coeff0_2*dmats1[2][5] + coeff0_3*dmats1[3][5] + coeff0_4*dmats1[4][5] + coeff0_5*dmats1[5][5];
          new_coeff1_0 = coeff1_0*dmats1[0][0] + coeff1_1*dmats1[1][0] + coeff1_2*dmats1[2][0] + coeff1_3*dmats1[3][0] + coeff1_4*dmats1[4][0] + coeff1_5*dmats1[5][0];
          new_coeff1_1 = coeff1_0*dmats1[0][1] + coeff1_1*dmats1[1][1] + coeff1_2*dmats1[2][1] + coeff1_3*dmats1[3][1] + coeff1_4*dmats1[4][1] + coeff1_5*dmats1[5][1];
          new_coeff1_2 = coeff1_0*dmats1[0][2] + coeff1_1*dmats1[1][2] + coeff1_2*dmats1[2][2] + coeff1_3*dmats1[3][2] + coeff1_4*dmats1[4][2] + coeff1_5*dmats1[5][2];
          new_coeff1_3 = coeff1_0*dmats1[0][3] + coeff1_1*dmats1[1][3] + coeff1_2*dmats1[2][3] + coeff1_3*dmats1[3][3] + coeff1_4*dmats1[4][3] + coeff1_5*dmats1[5][3];
          new_coeff1_4 = coeff1_0*dmats1[0][4] + coeff1_1*dmats1[1][4] + coeff1_2*dmats1[2][4] + coeff1_3*dmats1[3][4] + coeff1_4*dmats1[4][4] + coeff1_5*dmats1[5][4];
          new_coeff1_5 = coeff1_0*dmats1[0][5] + coeff1_1*dmats1[1][5] + coeff1_2*dmats1[2][5] + coeff1_3*dmats1[3][5] + coeff1_4*dmats1[4][5] + coeff1_5*dmats1[5][5];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      // Correct values by the contravariant Piola transform
      const double tmp0_0 = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3 + new_coeff0_4*basisvalue4 + new_coeff0_5*basisvalue5;
      const double tmp0_1 = new_coeff1_0*basisvalue0 + new_coeff1_1*basisvalue1 + new_coeff1_2*basisvalue2 + new_coeff1_3*basisvalue3 + new_coeff1_4*basisvalue4 + new_coeff1_5*basisvalue5;
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
void darcyflow_0_finite_element_1_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_0_finite_element_1_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[12][9][2] = {{{0.75, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.75, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}};
    static const double W[12][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}};
    static const double D[12][9][2] = {{{1, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, -1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, -1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, -1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.866025403784438, 0.433012701892219}, {0.866025403784439, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784439, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}}, {{0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}}, {{-0.17542966950853, 0.148329095604429}, {-0.0180753884489578, 0.444221856552778}, {0.0340120331760423, 0.307261606416198}, {-0.396781833524282, 3.42776476771047e-16}, {0.0799728864560174, -2.52114035295885e-16}, {0.171602319815636, -3.84645013221013e-15}, {-0.0271005739041022, -0.14832909560443}, {0.426146468103823, -0.44422185655278}, {0.341273639592246, -0.307261606416206}}};
    
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
    static const unsigned int ns[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 9, 9};
    for (unsigned int j = 0; j < ns[i]; j++) {
      // Evaluate basis functions for affine mapping
      const double w0 = 1.0 - X[i][j][0] - X[i][j][1];
      const double w1 = X[i][j][0];
      const double w2 = X[i][j][1];
      
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
        result += values[k]*D[i][j][k];
      // Multiply by weights
      result *= W[i][j];
    
    } // End for
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void darcyflow_0_finite_element_1_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_0_finite_element_1_0::interpolate_vertex_values(double* vertex_values,
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
    vertex_values[0] = (1.0/detJ)*(dof_values[3]*3*J_00 + dof_values[4]*(-3*J_00) + dof_values[5]*J_00 + dof_values[6]*(-3*J_01) + dof_values[7]*3*J_01 + dof_values[8]*J_01);
    vertex_values[2] = (1.0/detJ)*(dof_values[0]*3*J_00 + dof_values[1]*(-3*J_00) + dof_values[2]*J_00 + dof_values[6]*(J_00 - J_01) + dof_values[7]*(-3*J_00 + 3*J_01) + dof_values[8]*(3*J_00 - 3*J_01));
    vertex_values[4] = (1.0/detJ)*(dof_values[0]*J_01 + dof_values[1]*(-3*J_01) + dof_values[2]*3*J_01 + dof_values[3]*(J_00 - J_01) + dof_values[4]*(-3*J_00 + 3*J_01) + dof_values[5]*(3*J_00 - 3*J_01));
    vertex_values[1] = (1.0/detJ)*(dof_values[3]*3*J_10 + dof_values[4]*(-3*J_10) + dof_values[5]*J_10 + dof_values[6]*(-3*J_11) + dof_values[7]*3*J_11 + dof_values[8]*J_11);
    vertex_values[3] = (1.0/detJ)*(dof_values[0]*3*J_10 + dof_values[1]*(-3*J_10) + dof_values[2]*J_10 + dof_values[6]*(J_10 - J_11) + dof_values[7]*(-3*J_10 + 3*J_11) + dof_values[8]*(3*J_10 - 3*J_11));
    vertex_values[5] = (1.0/detJ)*(dof_values[0]*J_11 + dof_values[1]*(-3*J_11) + dof_values[2]*3*J_11 + dof_values[3]*(J_10 - J_11) + dof_values[4]*(-3*J_10 + 3*J_11) + dof_values[5]*(3*J_10 - 3*J_11));
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_0_finite_element_1_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_0_finite_element_1_0::create_sub_element(unsigned int i) const
{
    return new darcyflow_0_finite_element_1_0();
}


/// Constructor
darcyflow_0_finite_element_1_1::darcyflow_0_finite_element_1_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_0_finite_element_1_1::~darcyflow_0_finite_element_1_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_0_finite_element_1_1::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return the cell shape
ufc::shape darcyflow_0_finite_element_1_1::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_0_finite_element_1_1::space_dimension() const
{
    return 3;
}

/// Return the rank of the value space
unsigned int darcyflow_0_finite_element_1_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_0_finite_element_1_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void darcyflow_0_finite_element_1_1::evaluate_basis(unsigned int i,
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
void darcyflow_0_finite_element_1_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_0_finite_element_1_1::evaluate_basis_derivatives(unsigned int i,
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
void darcyflow_0_finite_element_1_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_0_finite_element_1_1::evaluate_dof(unsigned int i,
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
void darcyflow_0_finite_element_1_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_0_finite_element_1_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_0_finite_element_1_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_0_finite_element_1_1::create_sub_element(unsigned int i) const
{
    return new darcyflow_0_finite_element_1_1();
}


/// Constructor
darcyflow_0_finite_element_1::darcyflow_0_finite_element_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_0_finite_element_1::~darcyflow_0_finite_element_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_0_finite_element_1::signature() const
{
    return "MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) })";
}

/// Return the cell shape
ufc::shape darcyflow_0_finite_element_1::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_0_finite_element_1::space_dimension() const
{
    return 15;
}

/// Return the rank of the value space
unsigned int darcyflow_0_finite_element_1::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_0_finite_element_1::value_dimension(unsigned int i) const
{
    return 3;
}

/// Evaluate basis function i at given point in cell
void darcyflow_0_finite_element_1::evaluate_basis(unsigned int i,
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
    
    if (0 <= i && i <= 11)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
      const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
      const double psitilde_bs_1_0 = 1;
      const double psitilde_bs_1_1 = 2.5*y + 1.5;
      const double psitilde_bs_2_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
      const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
      const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
      const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
      const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
      const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
      // Table(s) of coefficients
      static const double coefficients0[12][6] =   \
      {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
      {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
      {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
      {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
      {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
      {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
      {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
      {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
      {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
      {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
      {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
      {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
      static const double coefficients1[12][6] =   \
      {{0, 0, 0.2, 0, 0, 0.163299316185545},
      {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
      {0, 0, 0.6, 0, 0, 0.489897948556636},
      {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
      {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
      {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
      {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
      {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
      {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
      {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
      {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
      {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
      const double coeff0_4 =   coefficients0[dof][4];
      const double coeff0_5 =   coefficients0[dof][5];
      const double coeff1_0 =   coefficients1[dof][0];
      const double coeff1_1 =   coefficients1[dof][1];
      const double coeff1_2 =   coefficients1[dof][2];
      const double coeff1_3 =   coefficients1[dof][3];
      const double coeff1_4 =   coefficients1[dof][4];
      const double coeff1_5 =   coefficients1[dof][5];
    
      // Compute value(s)
      const double tmp0_0 = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3 + coeff0_4*basisvalue4 + coeff0_5*basisvalue5;
      const double tmp0_1 = coeff1_0*basisvalue0 + coeff1_1*basisvalue1 + coeff1_2*basisvalue2 + coeff1_3*basisvalue3 + coeff1_4*basisvalue4 + coeff1_5*basisvalue5;
      // Using contravariant Piola transform to map values back to the physical element
      values[0] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
      values[1] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
    }
    
    if (12 <= i && i <= 14)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 12;
    
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
      static const double coefficients0[3][3] =   \
      {{0.471404520791032, -0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0, 0.333333333333333}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
    
      // Compute value(s)
      values[2] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2;
    }
    
}

/// Evaluate all basis functions at given point in cell
void darcyflow_0_finite_element_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_0_finite_element_1::evaluate_basis_derivatives(unsigned int i,
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
    
    if (0 <= i && i <= 11)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
      const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
      const double psitilde_bs_1_0 = 1;
      const double psitilde_bs_1_1 = 2.5*y + 1.5;
      const double psitilde_bs_2_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
      const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
      const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
      const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
      const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
      const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
      // Table(s) of coefficients
      static const double coefficients0[12][6] =   \
      {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
      {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
      {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
      {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
      {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
      {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
      {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
      {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
      {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
      {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
      {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
      {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
      static const double coefficients1[12][6] =   \
      {{0, 0, 0.2, 0, 0, 0.163299316185545},
      {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
      {0, 0, 0.6, 0, 0, 0.489897948556636},
      {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
      {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
      {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
      {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
      {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
      {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
      {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
      {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
      {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      static const double dmats0[6][6] =   \
      {{0, 0, 0, 0, 0, 0},
      {4.89897948556635, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
      {0, 9.48683298050514, 0, 0, 0, 0},
      {4, 0, 7.07106781186548, 0, 0, 0},
      {0, 0, 0, 0, 0, 0}};
    
      static const double dmats1[6][6] =   \
      {{0, 0, 0, 0, 0, 0},
      {2.44948974278318, 0, 0, 0, 0, 0},
      {4.24264068711928, 0, 0, 0, 0, 0},
      {2.58198889747161, 4.74341649025257, -0.912870929175277, 0, 0, 0},
      {2, 6.12372435695795, 3.53553390593274, 0, 0, 0},
      {-2.3094010767585, 0, 8.16496580927726, 0, 0, 0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [2*num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
      double coeff0_1 = 0;
      double coeff0_2 = 0;
      double coeff0_3 = 0;
      double coeff0_4 = 0;
      double coeff0_5 = 0;
      double coeff1_0 = 0;
      double coeff1_1 = 0;
      double coeff1_2 = 0;
      double coeff1_3 = 0;
      double coeff1_4 = 0;
      double coeff1_5 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
      double new_coeff0_1 = 0;
      double new_coeff0_2 = 0;
      double new_coeff0_3 = 0;
      double new_coeff0_4 = 0;
      double new_coeff0_5 = 0;
      double new_coeff1_0 = 0;
      double new_coeff1_1 = 0;
      double new_coeff1_2 = 0;
      double new_coeff1_3 = 0;
      double new_coeff1_4 = 0;
      double new_coeff1_5 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
        new_coeff0_1 = coefficients0[dof][1];
        new_coeff0_2 = coefficients0[dof][2];
        new_coeff0_3 = coefficients0[dof][3];
        new_coeff0_4 = coefficients0[dof][4];
        new_coeff0_5 = coefficients0[dof][5];
        new_coeff1_0 = coefficients1[dof][0];
        new_coeff1_1 = coefficients1[dof][1];
        new_coeff1_2 = coefficients1[dof][2];
        new_coeff1_3 = coefficients1[dof][3];
        new_coeff1_4 = coefficients1[dof][4];
        new_coeff1_5 = coefficients1[dof][5];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
          coeff0_1 = new_coeff0_1;
          coeff0_2 = new_coeff0_2;
          coeff0_3 = new_coeff0_3;
          coeff0_4 = new_coeff0_4;
          coeff0_5 = new_coeff0_5;
          coeff1_0 = new_coeff1_0;
          coeff1_1 = new_coeff1_1;
          coeff1_2 = new_coeff1_2;
          coeff1_3 = new_coeff1_3;
          coeff1_4 = new_coeff1_4;
          coeff1_5 = new_coeff1_5;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0] + coeff0_4*dmats0[4][0] + coeff0_5*dmats0[5][0];
            new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1] + coeff0_4*dmats0[4][1] + coeff0_5*dmats0[5][1];
            new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2] + coeff0_4*dmats0[4][2] + coeff0_5*dmats0[5][2];
            new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3] + coeff0_4*dmats0[4][3] + coeff0_5*dmats0[5][3];
            new_coeff0_4 = coeff0_0*dmats0[0][4] + coeff0_1*dmats0[1][4] + coeff0_2*dmats0[2][4] + coeff0_3*dmats0[3][4] + coeff0_4*dmats0[4][4] + coeff0_5*dmats0[5][4];
            new_coeff0_5 = coeff0_0*dmats0[0][5] + coeff0_1*dmats0[1][5] + coeff0_2*dmats0[2][5] + coeff0_3*dmats0[3][5] + coeff0_4*dmats0[4][5] + coeff0_5*dmats0[5][5];
            new_coeff1_0 = coeff1_0*dmats0[0][0] + coeff1_1*dmats0[1][0] + coeff1_2*dmats0[2][0] + coeff1_3*dmats0[3][0] + coeff1_4*dmats0[4][0] + coeff1_5*dmats0[5][0];
            new_coeff1_1 = coeff1_0*dmats0[0][1] + coeff1_1*dmats0[1][1] + coeff1_2*dmats0[2][1] + coeff1_3*dmats0[3][1] + coeff1_4*dmats0[4][1] + coeff1_5*dmats0[5][1];
            new_coeff1_2 = coeff1_0*dmats0[0][2] + coeff1_1*dmats0[1][2] + coeff1_2*dmats0[2][2] + coeff1_3*dmats0[3][2] + coeff1_4*dmats0[4][2] + coeff1_5*dmats0[5][2];
            new_coeff1_3 = coeff1_0*dmats0[0][3] + coeff1_1*dmats0[1][3] + coeff1_2*dmats0[2][3] + coeff1_3*dmats0[3][3] + coeff1_4*dmats0[4][3] + coeff1_5*dmats0[5][3];
            new_coeff1_4 = coeff1_0*dmats0[0][4] + coeff1_1*dmats0[1][4] + coeff1_2*dmats0[2][4] + coeff1_3*dmats0[3][4] + coeff1_4*dmats0[4][4] + coeff1_5*dmats0[5][4];
            new_coeff1_5 = coeff1_0*dmats0[0][5] + coeff1_1*dmats0[1][5] + coeff1_2*dmats0[2][5] + coeff1_3*dmats0[3][5] + coeff1_4*dmats0[4][5] + coeff1_5*dmats0[5][5];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0] + coeff0_4*dmats1[4][0] + coeff0_5*dmats1[5][0];
            new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1] + coeff0_4*dmats1[4][1] + coeff0_5*dmats1[5][1];
            new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2] + coeff0_4*dmats1[4][2] + coeff0_5*dmats1[5][2];
            new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3] + coeff0_4*dmats1[4][3] + coeff0_5*dmats1[5][3];
            new_coeff0_4 = coeff0_0*dmats1[0][4] + coeff0_1*dmats1[1][4] + coeff0_2*dmats1[2][4] + coeff0_3*dmats1[3][4] + coeff0_4*dmats1[4][4] + coeff0_5*dmats1[5][4];
            new_coeff0_5 = coeff0_0*dmats1[0][5] + coeff0_1*dmats1[1][5] + coeff0_2*dmats1[2][5] + coeff0_3*dmats1[3][5] + coeff0_4*dmats1[4][5] + coeff0_5*dmats1[5][5];
            new_coeff1_0 = coeff1_0*dmats1[0][0] + coeff1_1*dmats1[1][0] + coeff1_2*dmats1[2][0] + coeff1_3*dmats1[3][0] + coeff1_4*dmats1[4][0] + coeff1_5*dmats1[5][0];
            new_coeff1_1 = coeff1_0*dmats1[0][1] + coeff1_1*dmats1[1][1] + coeff1_2*dmats1[2][1] + coeff1_3*dmats1[3][1] + coeff1_4*dmats1[4][1] + coeff1_5*dmats1[5][1];
            new_coeff1_2 = coeff1_0*dmats1[0][2] + coeff1_1*dmats1[1][2] + coeff1_2*dmats1[2][2] + coeff1_3*dmats1[3][2] + coeff1_4*dmats1[4][2] + coeff1_5*dmats1[5][2];
            new_coeff1_3 = coeff1_0*dmats1[0][3] + coeff1_1*dmats1[1][3] + coeff1_2*dmats1[2][3] + coeff1_3*dmats1[3][3] + coeff1_4*dmats1[4][3] + coeff1_5*dmats1[5][3];
            new_coeff1_4 = coeff1_0*dmats1[0][4] + coeff1_1*dmats1[1][4] + coeff1_2*dmats1[2][4] + coeff1_3*dmats1[3][4] + coeff1_4*dmats1[4][4] + coeff1_5*dmats1[5][4];
            new_coeff1_5 = coeff1_0*dmats1[0][5] + coeff1_1*dmats1[1][5] + coeff1_2*dmats1[2][5] + coeff1_3*dmats1[3][5] + coeff1_4*dmats1[4][5] + coeff1_5*dmats1[5][5];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        // Correct values by the contravariant Piola transform
        const double tmp0_0 = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3 + new_coeff0_4*basisvalue4 + new_coeff0_5*basisvalue5;
        const double tmp0_1 = new_coeff1_0*basisvalue0 + new_coeff1_1*basisvalue1 + new_coeff1_2*basisvalue2 + new_coeff1_3*basisvalue3 + new_coeff1_4*basisvalue4 + new_coeff1_5*basisvalue5;
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
    
    if (12 <= i && i <= 14)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 12;
    
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
      static const double coefficients0[3][3] =   \
      {{0.471404520791032, -0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0, 0.333333333333333}};
    
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
void darcyflow_0_finite_element_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_0_finite_element_1::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[15][9][2] = {{{0.75, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.75, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}};
    static const double W[15][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}};
    static const double D[15][9][3] = {{{1, 1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0.866025403784438, 0.433012701892219, 0}, {0.866025403784439, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784439, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}}, {{0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}}, {{-0.17542966950853, 0.148329095604429, 0}, {-0.0180753884489578, 0.444221856552778, 0}, {0.0340120331760423, 0.307261606416198, 0}, {-0.396781833524282, 3.42776476771047e-16, 0}, {0.0799728864560174, -2.52114035295885e-16, 0}, {0.171602319815636, -3.84645013221013e-15, 0}, {-0.0271005739041022, -0.14832909560443, 0}, {0.426146468103823, -0.44422185655278, 0}, {0.341273639592246, -0.307261606416206, 0}}, {{0, 0, 1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    
    static const unsigned int mappings[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
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
    static const unsigned int ns[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 9, 9, 1, 1, 1};
    for (unsigned int j = 0; j < ns[i]; j++) {
      // Evaluate basis functions for affine mapping
      const double w0 = 1.0 - X[i][j][0] - X[i][j][1];
      const double w1 = X[i][j][0];
      const double w2 = X[i][j][1];
      
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
        result += values[k]*D[i][j][k];
      // Multiply by weights
      result *= W[i][j];
    
    } // End for
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void darcyflow_0_finite_element_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_0_finite_element_1::interpolate_vertex_values(double* vertex_values,
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
    vertex_values[0] = (1.0/detJ)*(dof_values[3]*3*J_00 + dof_values[4]*(-3*J_00) + dof_values[5]*J_00 + dof_values[6]*(-3*J_01) + dof_values[7]*3*J_01 + dof_values[8]*J_01);
    vertex_values[3] = (1.0/detJ)*(dof_values[0]*3*J_00 + dof_values[1]*(-3*J_00) + dof_values[2]*J_00 + dof_values[6]*(J_00 - J_01) + dof_values[7]*(-3*J_00 + 3*J_01) + dof_values[8]*(3*J_00 - 3*J_01));
    vertex_values[6] = (1.0/detJ)*(dof_values[0]*J_01 + dof_values[1]*(-3*J_01) + dof_values[2]*3*J_01 + dof_values[3]*(J_00 - J_01) + dof_values[4]*(-3*J_00 + 3*J_01) + dof_values[5]*(3*J_00 - 3*J_01));
    vertex_values[1] = (1.0/detJ)*(dof_values[3]*3*J_10 + dof_values[4]*(-3*J_10) + dof_values[5]*J_10 + dof_values[6]*(-3*J_11) + dof_values[7]*3*J_11 + dof_values[8]*J_11);
    vertex_values[4] = (1.0/detJ)*(dof_values[0]*3*J_10 + dof_values[1]*(-3*J_10) + dof_values[2]*J_10 + dof_values[6]*(J_10 - J_11) + dof_values[7]*(-3*J_10 + 3*J_11) + dof_values[8]*(3*J_10 - 3*J_11));
    vertex_values[7] = (1.0/detJ)*(dof_values[0]*J_11 + dof_values[1]*(-3*J_11) + dof_values[2]*3*J_11 + dof_values[3]*(J_10 - J_11) + dof_values[4]*(-3*J_10 + 3*J_11) + dof_values[5]*(3*J_10 - 3*J_11));
    // Evaluate at vertices and use affine mapping
    vertex_values[2] = dof_values[12];
    vertex_values[5] = dof_values[13];
    vertex_values[8] = dof_values[14];
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_0_finite_element_1::num_sub_elements() const
{
    return 2;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_0_finite_element_1::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_0_finite_element_1_0();
      break;
    case 1:
      return new darcyflow_0_finite_element_1_1();
      break;
    }
    return 0;
}


/// Constructor
darcyflow_0_finite_element_2::darcyflow_0_finite_element_2() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_0_finite_element_2::~darcyflow_0_finite_element_2()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_0_finite_element_2::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return the cell shape
ufc::shape darcyflow_0_finite_element_2::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_0_finite_element_2::space_dimension() const
{
    return 3;
}

/// Return the rank of the value space
unsigned int darcyflow_0_finite_element_2::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_0_finite_element_2::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void darcyflow_0_finite_element_2::evaluate_basis(unsigned int i,
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
void darcyflow_0_finite_element_2::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_0_finite_element_2::evaluate_basis_derivatives(unsigned int i,
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
void darcyflow_0_finite_element_2::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_0_finite_element_2::evaluate_dof(unsigned int i,
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
void darcyflow_0_finite_element_2::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_0_finite_element_2::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_0_finite_element_2::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_0_finite_element_2::create_sub_element(unsigned int i) const
{
    return new darcyflow_0_finite_element_2();
}


/// Constructor
darcyflow_0_finite_element_3::darcyflow_0_finite_element_3() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_0_finite_element_3::~darcyflow_0_finite_element_3()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_0_finite_element_3::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return the cell shape
ufc::shape darcyflow_0_finite_element_3::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_0_finite_element_3::space_dimension() const
{
    return 3;
}

/// Return the rank of the value space
unsigned int darcyflow_0_finite_element_3::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_0_finite_element_3::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void darcyflow_0_finite_element_3::evaluate_basis(unsigned int i,
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
void darcyflow_0_finite_element_3::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_0_finite_element_3::evaluate_basis_derivatives(unsigned int i,
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
void darcyflow_0_finite_element_3::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_0_finite_element_3::evaluate_dof(unsigned int i,
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
void darcyflow_0_finite_element_3::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_0_finite_element_3::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_0_finite_element_3::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_0_finite_element_3::create_sub_element(unsigned int i) const
{
    return new darcyflow_0_finite_element_3();
}

/// Constructor
darcyflow_0_dof_map_0_0::darcyflow_0_dof_map_0_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_0_dof_map_0_0::~darcyflow_0_dof_map_0_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_0_dof_map_0_0::signature() const
{
    return "FFC dof map for FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_0_dof_map_0_0::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_0_dof_map_0_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[1] + 3*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_0_dof_map_0_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_0_dof_map_0_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_0_dof_map_0_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_0_dof_map_0_0::local_dimension(const ufc::cell& c) const
{
    return 12;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_0_dof_map_0_0::max_local_dimension() const
{
    return 12;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_0_dof_map_0_0::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_0_dof_map_0_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_0_dof_map_0_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_0_dof_map_0_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[1][0];
    dofs[1] = 3*c.entity_indices[1][0] + 1;
    dofs[2] = 3*c.entity_indices[1][0] + 2;
    dofs[3] = 3*c.entity_indices[1][1];
    dofs[4] = 3*c.entity_indices[1][1] + 1;
    dofs[5] = 3*c.entity_indices[1][1] + 2;
    dofs[6] = 3*c.entity_indices[1][2];
    dofs[7] = 3*c.entity_indices[1][2] + 1;
    dofs[8] = 3*c.entity_indices[1][2] + 2;
    unsigned int offset = 3*m.num_entities[1];
    dofs[9] = offset + 3*c.entity_indices[2][0];
    dofs[10] = offset + 3*c.entity_indices[2][0] + 1;
    dofs[11] = offset + 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_0_dof_map_0_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      break;
    case 1:
      dofs[0] = 3;
      dofs[1] = 4;
      dofs[2] = 5;
      break;
    case 2:
      dofs[0] = 6;
      dofs[1] = 7;
      dofs[2] = 8;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void darcyflow_0_dof_map_0_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_0_dof_map_0_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.75*x[1][0] + 0.25*x[2][0];
    coordinates[0][1] = 0.75*x[1][1] + 0.25*x[2][1];
    coordinates[1][0] = 0.5*x[1][0] + 0.5*x[2][0];
    coordinates[1][1] = 0.5*x[1][1] + 0.5*x[2][1];
    coordinates[2][0] = 0.25*x[1][0] + 0.75*x[2][0];
    coordinates[2][1] = 0.25*x[1][1] + 0.75*x[2][1];
    coordinates[3][0] = 0.75*x[0][0] + 0.25*x[2][0];
    coordinates[3][1] = 0.75*x[0][1] + 0.25*x[2][1];
    coordinates[4][0] = 0.5*x[0][0] + 0.5*x[2][0];
    coordinates[4][1] = 0.5*x[0][1] + 0.5*x[2][1];
    coordinates[5][0] = 0.25*x[0][0] + 0.75*x[2][0];
    coordinates[5][1] = 0.25*x[0][1] + 0.75*x[2][1];
    coordinates[6][0] = 0.75*x[0][0] + 0.25*x[1][0];
    coordinates[6][1] = 0.75*x[0][1] + 0.25*x[1][1];
    coordinates[7][0] = 0.5*x[0][0] + 0.5*x[1][0];
    coordinates[7][1] = 0.5*x[0][1] + 0.5*x[1][1];
    coordinates[8][0] = 0.25*x[0][0] + 0.75*x[1][0];
    coordinates[8][1] = 0.25*x[0][1] + 0.75*x[1][1];
    coordinates[9][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[9][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[10][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[10][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[11][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[11][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int darcyflow_0_dof_map_0_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_0_dof_map_0_0::create_sub_dof_map(unsigned int i) const
{
    return new darcyflow_0_dof_map_0_0();
}


/// Constructor
darcyflow_0_dof_map_0_1::darcyflow_0_dof_map_0_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_0_dof_map_0_1::~darcyflow_0_dof_map_0_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_0_dof_map_0_1::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_0_dof_map_0_1::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_0_dof_map_0_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_0_dof_map_0_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_0_dof_map_0_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_0_dof_map_0_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_0_dof_map_0_1::local_dimension(const ufc::cell& c) const
{
    return 3;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_0_dof_map_0_1::max_local_dimension() const
{
    return 3;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_0_dof_map_0_1::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_0_dof_map_0_1::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_0_dof_map_0_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_0_dof_map_0_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[2][0];
    dofs[1] = 3*c.entity_indices[2][0] + 1;
    dofs[2] = 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_0_dof_map_0_1::tabulate_facet_dofs(unsigned int* dofs,
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
void darcyflow_0_dof_map_0_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_0_dof_map_0_1::tabulate_coordinates(double** coordinates,
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
unsigned int darcyflow_0_dof_map_0_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_0_dof_map_0_1::create_sub_dof_map(unsigned int i) const
{
    return new darcyflow_0_dof_map_0_1();
}


/// Constructor
darcyflow_0_dof_map_0::darcyflow_0_dof_map_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_0_dof_map_0::~darcyflow_0_dof_map_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_0_dof_map_0::signature() const
{
    return "FFC dof map for MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) })";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_0_dof_map_0::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_0_dof_map_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[1] + 6*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_0_dof_map_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_0_dof_map_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_0_dof_map_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_0_dof_map_0::local_dimension(const ufc::cell& c) const
{
    return 15;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_0_dof_map_0::max_local_dimension() const
{
    return 15;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_0_dof_map_0::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_0_dof_map_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_0_dof_map_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_0_dof_map_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[1][0];
    dofs[1] = 3*c.entity_indices[1][0] + 1;
    dofs[2] = 3*c.entity_indices[1][0] + 2;
    dofs[3] = 3*c.entity_indices[1][1];
    dofs[4] = 3*c.entity_indices[1][1] + 1;
    dofs[5] = 3*c.entity_indices[1][1] + 2;
    dofs[6] = 3*c.entity_indices[1][2];
    dofs[7] = 3*c.entity_indices[1][2] + 1;
    dofs[8] = 3*c.entity_indices[1][2] + 2;
    unsigned int offset = 3*m.num_entities[1];
    dofs[9] = offset + 3*c.entity_indices[2][0];
    dofs[10] = offset + 3*c.entity_indices[2][0] + 1;
    dofs[11] = offset + 3*c.entity_indices[2][0] + 2;
    offset = offset + 3*m.num_entities[2];
    dofs[12] = offset + 3*c.entity_indices[2][0];
    dofs[13] = offset + 3*c.entity_indices[2][0] + 1;
    dofs[14] = offset + 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_0_dof_map_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      break;
    case 1:
      dofs[0] = 3;
      dofs[1] = 4;
      dofs[2] = 5;
      break;
    case 2:
      dofs[0] = 6;
      dofs[1] = 7;
      dofs[2] = 8;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void darcyflow_0_dof_map_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_0_dof_map_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.75*x[1][0] + 0.25*x[2][0];
    coordinates[0][1] = 0.75*x[1][1] + 0.25*x[2][1];
    coordinates[1][0] = 0.5*x[1][0] + 0.5*x[2][0];
    coordinates[1][1] = 0.5*x[1][1] + 0.5*x[2][1];
    coordinates[2][0] = 0.25*x[1][0] + 0.75*x[2][0];
    coordinates[2][1] = 0.25*x[1][1] + 0.75*x[2][1];
    coordinates[3][0] = 0.75*x[0][0] + 0.25*x[2][0];
    coordinates[3][1] = 0.75*x[0][1] + 0.25*x[2][1];
    coordinates[4][0] = 0.5*x[0][0] + 0.5*x[2][0];
    coordinates[4][1] = 0.5*x[0][1] + 0.5*x[2][1];
    coordinates[5][0] = 0.25*x[0][0] + 0.75*x[2][0];
    coordinates[5][1] = 0.25*x[0][1] + 0.75*x[2][1];
    coordinates[6][0] = 0.75*x[0][0] + 0.25*x[1][0];
    coordinates[6][1] = 0.75*x[0][1] + 0.25*x[1][1];
    coordinates[7][0] = 0.5*x[0][0] + 0.5*x[1][0];
    coordinates[7][1] = 0.5*x[0][1] + 0.5*x[1][1];
    coordinates[8][0] = 0.25*x[0][0] + 0.75*x[1][0];
    coordinates[8][1] = 0.25*x[0][1] + 0.75*x[1][1];
    coordinates[9][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[9][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[10][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[10][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[11][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[11][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[12][0] = x[0][0];
    coordinates[12][1] = x[0][1];
    coordinates[13][0] = x[1][0];
    coordinates[13][1] = x[1][1];
    coordinates[14][0] = x[2][0];
    coordinates[14][1] = x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int darcyflow_0_dof_map_0::num_sub_dof_maps() const
{
    return 2;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_0_dof_map_0::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_0_dof_map_0_0();
      break;
    case 1:
      return new darcyflow_0_dof_map_0_1();
      break;
    }
    return 0;
}


/// Constructor
darcyflow_0_dof_map_1_0::darcyflow_0_dof_map_1_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_0_dof_map_1_0::~darcyflow_0_dof_map_1_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_0_dof_map_1_0::signature() const
{
    return "FFC dof map for FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_0_dof_map_1_0::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_0_dof_map_1_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[1] + 3*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_0_dof_map_1_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_0_dof_map_1_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_0_dof_map_1_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_0_dof_map_1_0::local_dimension(const ufc::cell& c) const
{
    return 12;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_0_dof_map_1_0::max_local_dimension() const
{
    return 12;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_0_dof_map_1_0::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_0_dof_map_1_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_0_dof_map_1_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_0_dof_map_1_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[1][0];
    dofs[1] = 3*c.entity_indices[1][0] + 1;
    dofs[2] = 3*c.entity_indices[1][0] + 2;
    dofs[3] = 3*c.entity_indices[1][1];
    dofs[4] = 3*c.entity_indices[1][1] + 1;
    dofs[5] = 3*c.entity_indices[1][1] + 2;
    dofs[6] = 3*c.entity_indices[1][2];
    dofs[7] = 3*c.entity_indices[1][2] + 1;
    dofs[8] = 3*c.entity_indices[1][2] + 2;
    unsigned int offset = 3*m.num_entities[1];
    dofs[9] = offset + 3*c.entity_indices[2][0];
    dofs[10] = offset + 3*c.entity_indices[2][0] + 1;
    dofs[11] = offset + 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_0_dof_map_1_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      break;
    case 1:
      dofs[0] = 3;
      dofs[1] = 4;
      dofs[2] = 5;
      break;
    case 2:
      dofs[0] = 6;
      dofs[1] = 7;
      dofs[2] = 8;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void darcyflow_0_dof_map_1_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_0_dof_map_1_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.75*x[1][0] + 0.25*x[2][0];
    coordinates[0][1] = 0.75*x[1][1] + 0.25*x[2][1];
    coordinates[1][0] = 0.5*x[1][0] + 0.5*x[2][0];
    coordinates[1][1] = 0.5*x[1][1] + 0.5*x[2][1];
    coordinates[2][0] = 0.25*x[1][0] + 0.75*x[2][0];
    coordinates[2][1] = 0.25*x[1][1] + 0.75*x[2][1];
    coordinates[3][0] = 0.75*x[0][0] + 0.25*x[2][0];
    coordinates[3][1] = 0.75*x[0][1] + 0.25*x[2][1];
    coordinates[4][0] = 0.5*x[0][0] + 0.5*x[2][0];
    coordinates[4][1] = 0.5*x[0][1] + 0.5*x[2][1];
    coordinates[5][0] = 0.25*x[0][0] + 0.75*x[2][0];
    coordinates[5][1] = 0.25*x[0][1] + 0.75*x[2][1];
    coordinates[6][0] = 0.75*x[0][0] + 0.25*x[1][0];
    coordinates[6][1] = 0.75*x[0][1] + 0.25*x[1][1];
    coordinates[7][0] = 0.5*x[0][0] + 0.5*x[1][0];
    coordinates[7][1] = 0.5*x[0][1] + 0.5*x[1][1];
    coordinates[8][0] = 0.25*x[0][0] + 0.75*x[1][0];
    coordinates[8][1] = 0.25*x[0][1] + 0.75*x[1][1];
    coordinates[9][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[9][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[10][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[10][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[11][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[11][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int darcyflow_0_dof_map_1_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_0_dof_map_1_0::create_sub_dof_map(unsigned int i) const
{
    return new darcyflow_0_dof_map_1_0();
}


/// Constructor
darcyflow_0_dof_map_1_1::darcyflow_0_dof_map_1_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_0_dof_map_1_1::~darcyflow_0_dof_map_1_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_0_dof_map_1_1::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_0_dof_map_1_1::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_0_dof_map_1_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_0_dof_map_1_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_0_dof_map_1_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_0_dof_map_1_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_0_dof_map_1_1::local_dimension(const ufc::cell& c) const
{
    return 3;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_0_dof_map_1_1::max_local_dimension() const
{
    return 3;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_0_dof_map_1_1::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_0_dof_map_1_1::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_0_dof_map_1_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_0_dof_map_1_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[2][0];
    dofs[1] = 3*c.entity_indices[2][0] + 1;
    dofs[2] = 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_0_dof_map_1_1::tabulate_facet_dofs(unsigned int* dofs,
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
void darcyflow_0_dof_map_1_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_0_dof_map_1_1::tabulate_coordinates(double** coordinates,
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
unsigned int darcyflow_0_dof_map_1_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_0_dof_map_1_1::create_sub_dof_map(unsigned int i) const
{
    return new darcyflow_0_dof_map_1_1();
}


/// Constructor
darcyflow_0_dof_map_1::darcyflow_0_dof_map_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_0_dof_map_1::~darcyflow_0_dof_map_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_0_dof_map_1::signature() const
{
    return "FFC dof map for MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) })";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_0_dof_map_1::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_0_dof_map_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[1] + 6*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_0_dof_map_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_0_dof_map_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_0_dof_map_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_0_dof_map_1::local_dimension(const ufc::cell& c) const
{
    return 15;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_0_dof_map_1::max_local_dimension() const
{
    return 15;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_0_dof_map_1::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_0_dof_map_1::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_0_dof_map_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_0_dof_map_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[1][0];
    dofs[1] = 3*c.entity_indices[1][0] + 1;
    dofs[2] = 3*c.entity_indices[1][0] + 2;
    dofs[3] = 3*c.entity_indices[1][1];
    dofs[4] = 3*c.entity_indices[1][1] + 1;
    dofs[5] = 3*c.entity_indices[1][1] + 2;
    dofs[6] = 3*c.entity_indices[1][2];
    dofs[7] = 3*c.entity_indices[1][2] + 1;
    dofs[8] = 3*c.entity_indices[1][2] + 2;
    unsigned int offset = 3*m.num_entities[1];
    dofs[9] = offset + 3*c.entity_indices[2][0];
    dofs[10] = offset + 3*c.entity_indices[2][0] + 1;
    dofs[11] = offset + 3*c.entity_indices[2][0] + 2;
    offset = offset + 3*m.num_entities[2];
    dofs[12] = offset + 3*c.entity_indices[2][0];
    dofs[13] = offset + 3*c.entity_indices[2][0] + 1;
    dofs[14] = offset + 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_0_dof_map_1::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      break;
    case 1:
      dofs[0] = 3;
      dofs[1] = 4;
      dofs[2] = 5;
      break;
    case 2:
      dofs[0] = 6;
      dofs[1] = 7;
      dofs[2] = 8;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void darcyflow_0_dof_map_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_0_dof_map_1::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.75*x[1][0] + 0.25*x[2][0];
    coordinates[0][1] = 0.75*x[1][1] + 0.25*x[2][1];
    coordinates[1][0] = 0.5*x[1][0] + 0.5*x[2][0];
    coordinates[1][1] = 0.5*x[1][1] + 0.5*x[2][1];
    coordinates[2][0] = 0.25*x[1][0] + 0.75*x[2][0];
    coordinates[2][1] = 0.25*x[1][1] + 0.75*x[2][1];
    coordinates[3][0] = 0.75*x[0][0] + 0.25*x[2][0];
    coordinates[3][1] = 0.75*x[0][1] + 0.25*x[2][1];
    coordinates[4][0] = 0.5*x[0][0] + 0.5*x[2][0];
    coordinates[4][1] = 0.5*x[0][1] + 0.5*x[2][1];
    coordinates[5][0] = 0.25*x[0][0] + 0.75*x[2][0];
    coordinates[5][1] = 0.25*x[0][1] + 0.75*x[2][1];
    coordinates[6][0] = 0.75*x[0][0] + 0.25*x[1][0];
    coordinates[6][1] = 0.75*x[0][1] + 0.25*x[1][1];
    coordinates[7][0] = 0.5*x[0][0] + 0.5*x[1][0];
    coordinates[7][1] = 0.5*x[0][1] + 0.5*x[1][1];
    coordinates[8][0] = 0.25*x[0][0] + 0.75*x[1][0];
    coordinates[8][1] = 0.25*x[0][1] + 0.75*x[1][1];
    coordinates[9][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[9][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[10][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[10][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[11][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[11][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[12][0] = x[0][0];
    coordinates[12][1] = x[0][1];
    coordinates[13][0] = x[1][0];
    coordinates[13][1] = x[1][1];
    coordinates[14][0] = x[2][0];
    coordinates[14][1] = x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int darcyflow_0_dof_map_1::num_sub_dof_maps() const
{
    return 2;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_0_dof_map_1::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_0_dof_map_1_0();
      break;
    case 1:
      return new darcyflow_0_dof_map_1_1();
      break;
    }
    return 0;
}


/// Constructor
darcyflow_0_dof_map_2::darcyflow_0_dof_map_2() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_0_dof_map_2::~darcyflow_0_dof_map_2()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_0_dof_map_2::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_0_dof_map_2::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_0_dof_map_2::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_0_dof_map_2::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_0_dof_map_2::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_0_dof_map_2::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_0_dof_map_2::local_dimension(const ufc::cell& c) const
{
    return 3;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_0_dof_map_2::max_local_dimension() const
{
    return 3;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_0_dof_map_2::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_0_dof_map_2::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_0_dof_map_2::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_0_dof_map_2::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[2][0];
    dofs[1] = 3*c.entity_indices[2][0] + 1;
    dofs[2] = 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_0_dof_map_2::tabulate_facet_dofs(unsigned int* dofs,
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
void darcyflow_0_dof_map_2::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_0_dof_map_2::tabulate_coordinates(double** coordinates,
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
unsigned int darcyflow_0_dof_map_2::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_0_dof_map_2::create_sub_dof_map(unsigned int i) const
{
    return new darcyflow_0_dof_map_2();
}


/// Constructor
darcyflow_0_dof_map_3::darcyflow_0_dof_map_3() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_0_dof_map_3::~darcyflow_0_dof_map_3()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_0_dof_map_3::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_0_dof_map_3::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_0_dof_map_3::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_0_dof_map_3::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_0_dof_map_3::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_0_dof_map_3::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_0_dof_map_3::local_dimension(const ufc::cell& c) const
{
    return 3;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_0_dof_map_3::max_local_dimension() const
{
    return 3;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_0_dof_map_3::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_0_dof_map_3::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_0_dof_map_3::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_0_dof_map_3::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[2][0];
    dofs[1] = 3*c.entity_indices[2][0] + 1;
    dofs[2] = 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_0_dof_map_3::tabulate_facet_dofs(unsigned int* dofs,
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
void darcyflow_0_dof_map_3::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_0_dof_map_3::tabulate_coordinates(double** coordinates,
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
unsigned int darcyflow_0_dof_map_3::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_0_dof_map_3::create_sub_dof_map(unsigned int i) const
{
    return new darcyflow_0_dof_map_3();
}


/// Constructor
darcyflow_0_cell_integral_0_quadrature::darcyflow_0_cell_integral_0_quadrature() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_0_cell_integral_0_quadrature::~darcyflow_0_cell_integral_0_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void darcyflow_0_cell_integral_0_quadrature::tabulate_tensor(double* A,
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
    static const double W16[16] = {0.0235683681933823, 0.0353880678980859, 0.0225840492823699, 0.00542322591052525, 0.0441850885223617, 0.0663442161070497, 0.0423397245217463, 0.0101672595644788, 0.0441850885223617, 0.0663442161070497, 0.0423397245217463, 0.0101672595644788, 0.0235683681933823, 0.0353880678980859, 0.0225840492823699, 0.00542322591052525};
    // Quadrature points on the UFC reference element: (0.0654669945550145, 0.0571041961145177), (0.0502101232113698, 0.276843013638124), (0.028912084224389, 0.583590432368917), (0.00970378512694614, 0.860240135656219), (0.311164552244357, 0.0571041961145177), (0.238648659731443, 0.276843013638124), (0.137419104134574, 0.583590432368917), (0.0461220799064521, 0.860240135656219), (0.631731251641125, 0.0571041961145177), (0.484508326630433, 0.276843013638124), (0.278990463496509, 0.583590432368917), (0.0936377844373285, 0.860240135656219), (0.877428809330468, 0.0571041961145177), (0.672946863150506, 0.276843013638124), (0.387497483406694, 0.583590432368917), (0.130056079216834, 0.860240135656219)
    
    // Value of basis functions at quadrature points.
    static const double FE0[16][3] = \
    {{0.877428809330468, 0.0654669945550145, 0.0571041961145176},
    {0.672946863150506, 0.0502101232113698, 0.276843013638124},
    {0.387497483406694, 0.0289120842243891, 0.583590432368917},
    {0.130056079216834, 0.00970378512694621, 0.860240135656219},
    {0.631731251641125, 0.311164552244357, 0.0571041961145176},
    {0.484508326630433, 0.238648659731443, 0.276843013638124},
    {0.278990463496509, 0.137419104134574, 0.583590432368917},
    {0.0936377844373285, 0.0461220799064521, 0.860240135656219},
    {0.311164552244357, 0.631731251641125, 0.0571041961145176},
    {0.238648659731443, 0.484508326630433, 0.276843013638124},
    {0.137419104134574, 0.278990463496509, 0.583590432368917},
    {0.046122079906452, 0.0936377844373286, 0.860240135656219},
    {0.0654669945550145, 0.877428809330468, 0.0571041961145176},
    {0.0502101232113699, 0.672946863150506, 0.276843013638124},
    {0.028912084224389, 0.387497483406694, 0.583590432368917},
    {0.00970378512694601, 0.130056079216835, 0.860240135656219}};
    
    static const double FE1_C0[16][15] = \
    {{-0.170685419408646, 0.0632770454351407, -0.0568951398028823, 1.90956695348697, -1.81090106573336, 0.537438661988833, -0.0837472332962584, 0.224389606395399, -0.143833325915269, 0.304955814070681, -0.0686579473711441, -0.509482652297143, 0, 0, 0},
    {-0.135504030796704, 0.095727184619675, -0.0451680102655682, 0.564801177511215, 0.302308415306002, -0.127964406070977, -0.0551122218098253, 0.155392453885218, -0.125559819252447, 0.199647532384875, -0.0754897103884233, -0.18867814808636, 0, 0, 0},
    {-0.0817208009879782, 0.0930597127003873, -0.0272402669959928, -0.161166771276691, 0.822947593101012, 0.202671847366731, -0.0244055390678904, 0.0760513451317736, -0.0845555289160805, 0.0874389396839381, -0.0618218070765862, 0.0537851807981598, 0, 0, 0},
    {-0.0285463747060986, 0.0427174730797178, -0.00951545823536593, 0.334648085721877, -1.13529326105324, 1.7595884526666, -0.00597268364196135, 0.0214608255192889, -0.0320891492995035, 0.0210161747954732, -0.026304792549122, 0.0672194217104174, 0, 0, 0},
    {-0.352553385292486, -0.00505375558837244, -0.117517795097497, 0.35932828797004, -0.810753134457506, 0.104092029118977, -0.20691958031771, 0.531356955732914, -0.263151600072271, 1.05219647448152, -0.249878776901424, -1.69627960908635, 0, 0, 0},
    {-0.374226882444645, 0.275108785207995, -0.124742294148216, -0.193373497992773, 0.539292762769455, -0.360908880885765, -0.149521818457378, 0.42378593106297, -0.349447358135483, 0.715250385312202, -0.313831905261333, -0.470158416914786, 0, 0, 0},
    {-0.298953251316829, 0.382668961000559, -0.0996510837722759, -0.251463237600451, 0.533314186559004, 0.148447425935037, -0.0787221563513435, 0.257095396474964, -0.319882178737762, 0.338117679132613, -0.278928042748716, 0.397098477807007, 0, 0, 0},
    {-0.125602762189973, 0.196317351988756, -0.041867587396657, 0.409159548515142, -1.36148144458439, 1.76557777645597, -0.0241889399469618, 0.0902454672905814, -0.143281409639669, 0.091161943126126, -0.123346962202501, 0.335428751352995, 0, 0, 0},
    {0.499312490876999, -0.820308224417584, 0.166437496958998, -0.483384736364376, -0.0504157766387051, -0.189016114772586, 0.0861885635738531, -0.338814624106707, 0.579561424262146, 1.0839028833065, -0.304795887911834, -1.52261645577714, 0, 0, 0},
    {-0.045035068445964, 0.0820460334420768, -0.0150116894819879, -0.488533046433263, 0.528161921186228, -0.504671751110184, -0.00575894823295976, 0.0265295859479077, -0.0542878096949923, 0.833141611399128, -0.518025498610478, 0.175558421682723, 0, 0, 0},
    {-0.369957318157544, 0.618912409483447, -0.123319106052513, -0.139149325349192, 0.049211665612856, 0.130805602804731, -0.0610803332210388, 0.245479772494593, -0.43219609098902, 0.481219181340367, -0.526787115211975, 1.18089768552763, 0, 0, 0},
    {-0.228305145265995, 0.380769567287661, -0.0761017150886639, 0.532299389646885, -1.66855840821762, 1.77937444609225, -0.0379856095832482, 0.152072934255162, -0.266421250771412, 0.161959374586259, -0.245971710536827, 0.723202253594404, 0, 0, 0},
    {1.98700146466709, -2.00167698522372, 0.662333821555697, -0.22493671723455, 0.11495368480019, -0.204973512723965, 0.658664941416539, -1.97966370438878, 1.99067034480625, 0.385264820895113, -0.207757227496455, -0.0696121062113665, 0, 0, 0},
    {0.69830429429312, -0.393280508894289, 0.232768098097708, -0.182806989684644, 0.274115161456172, -0.49210067232822, 0.309024044447415, -0.850816186992535, 0.622048347943412, 0.498253638570492, -0.592690657752215, 1.44685485358008, 0, 0, 0},
    {-0.261566652340953, 0.691439580353324, -0.0871888841136485, 0.123312584151275, -0.403233299001872, 0.157986960409395, 0.0202793478894422, 0.0466301883347679, -0.369034884344046, 0.449900030147638, -0.68962283152653, 2.03906633544705, 0, 0, 0},
    {-0.288680735202969, 0.509915581165432, -0.0962269117343215, 0.646548383115959, -1.91308699052222, 1.79453396926834, -0.0409182002437069, 0.178063312221737, -0.343989446693585, 0.200338640410223, -0.336900413932389, 1.04940901655684, 0, 0, 0}};
    
    static const double FE1_C0_D01[16][15] = \
    {{0, 0.261867978220057, 0, -8.23526193295253, 13.2561278005652, -4.88993187850268, 0.0654669945550143, -0.130933989110028, -0.0654669945550149, -0.226784321576239, -0.130933989110027, 1.24214888615813, 0, 0, 0},
    {0, 0.200840492845478, 0, -4.91778018004221, 6.43808183862086, -1.41988141215591, 0.0502101232113698, -0.100420246422739, -0.0502101232113697, -0.173932968912771, -0.100420246422738, 0.9526701056737, 0, 0, 0},
    {0, 0.115648336897555, 0, -0.286695987180273, -3.07966301494679, 3.42418317057584, 0.0289120842243893, -0.0578241684487783, -0.0289120842243884, -0.100154397658705, -0.0578241684487774, 0.548568228310153, 0, 0, 0},
    {0, 0.0388151405077836, 0, 3.88999137714981, -11.6635373327763, 7.79295352588035, 0.00970378512694645, -0.0194075702538928, -0.00970378512694535, -0.0336148977312042, -0.0194075702538923, 0.184116377556096, 0, 0, 0},
    {0, 1.24465820897742, 0, -5.04119368299108, 9.81636199291445, -4.15283920543466, 0.311164552244358, -0.622329104488715, -0.311164552244352, -1.07790562800329, -0.622329104488706, 5.90393227319176, 0, 0, 0},
    {0, 0.954594638925761, 0, -2.46807920528127, 3.79994232733985, -0.8545658025957, 0.238648659731444, -0.477297319462887, -0.238648659731439, -0.82670320762615, -0.47729731946288, 4.5280399517872, 0, 0, 0},
    {0, 0.549676416538291, 0, 1.12389527165213, -4.59876129368938, 3.74970423030639, 0.137419104134575, -0.27483820826915, -0.137419104134571, -0.476033740583362, -0.274838208269145, 2.6073441785107, 0, 0, 0},
    {0, 0.184488319625805, 0, 4.36342920928339, -12.1733934596894, 7.90220841021887, 0.0461220799064525, -0.0922441598129048, -0.0461220799064505, -0.159771571497453, -0.0922441598129032, 0.875104937572045, 0, 0, 0},
    {1.32430381059195e-14, 2.52692500656447, 1.72647338950337e-14, -0.873826590833098, 5.3284282013597, -3.19113910724437, 0.631731251641128, -1.26346250328225, -0.631731251641113, -2.18838124914301, -1.26346250328224, 11.9862577457696, 0, 0, 0},
    {1.05321938408392e-14, 1.93803330652171, 1.39329205724286e-14, 0.728096464405605, 0.357906990753991, -0.116986801898738, 0.484508326630436, -0.969016653260869, -0.484508326630424, -1.67838607682817, -0.969016653260855, 9.1928991448139, 0, 0, 0},
    {0, 1.11596185398602, 0, 2.96432294335728, -6.58076032475645, 4.17441830839218, 0.27899046349651, -0.557980926993019, -0.278990463496503, -0.966451315206285, -0.557980926993011, 5.29347186069019, 0, 0, 0},
    {0, 0.374551137749309, 0, 4.98113336818478, -12.8386133231216, 8.04475552381149, 0.0936377844373292, -0.187275568874658, -0.0936377844373259, -0.32437080030727, -0.187275568874655, 1.77665204324296, 0, 0, 0},
    {1.88935637283063e-14, 3.50971523732183, 2.50388415057924e-14, 2.32024165912835, 1.88866239370892, -2.45404643417635, 0.877428809330472, -1.75485761866094, -0.87742880933045, -3.03950255557007, -1.75485761866091, 16.6480411328033, 0, 0, 0},
    {1.48450423703333e-14, 2.69178745260199, 1.98883507573745e-14, 3.17779743916655, -2.28023252052702, 0.448328807661474, 0.67294686315051, -1.34589372630102, -0.672946863150492, -2.33115631554155, -1.345893726301, 12.7682689909274, 0, 0, 0},
    {0, 1.54998993362676, 1.2652520484005e-14, 4.37491420218969, -8.09985860349904, 4.49993936812274, 0.387497483406696, -0.774994966813391, -0.387497483406686, -1.34233065813094, -0.774994966813379, 7.35224781089074, 0, 0, 0},
    {0, 0.520224316867331, 0, 5.45457120031836, -13.3484694500347, 8.15401040815001, 0.130056079216835, -0.26011215843367, -0.130056079216831, -0.450527474073519, -0.260112158433666, 2.46764060325891, 0, 0, 0}};
    
    static const double FE1_C0_D10[16][15] = \
    {{-2.21439606533984, 0.704680828017973, -0.738132021779954, -7.90657452129608, 4.80773928706665, -2.13228642799139, -1.11556083111041, 2.96925368400077, -1.83696725600937, 4.31798471998647, -0.983274403119016, -7.16120815229814, 0, 0, 0},
    {-2.39747852146357, 1.70569106886155, -0.799159507154526, -5.24830922095658, 1.82293706979804, -1.51884058945151, -0.972106370305028, 2.74337224776458, -2.22453165831306, 3.71534118528489, -1.45326578085351, -3.28143601042227, 0, 0, 0},
    {-2.65305498930733, 3.10306505568055, -0.88435166310244, -1.53746728428702, -2.34373855851117, -0.662492450220086, -0.771849146509137, 2.42804995612072, -2.76555750590063, 2.87407267038241, -2.10935669628905, 2.1345851696144, 0, 0, 0},
    {-2.88355457847664, 4.36333026160929, -0.961184859492204, 1.80927097018115, -6.10158460994874, 0.109831762349485, -0.591240938709048, 2.14366673691031, -3.2534984992598, 2.11534848617999, -2.70107270105854, 7.01919237724623, 0, 0, 0},
    {0.733974626932281, -1.26089963349678, 0.244658208977424, -4.71250627133463, 3.33355394093061, -1.39519375492337, 0.112926957336302, -0.470512123650032, 0.865705878573407, 1.76462080070531, -0.491879287740324, -2.4994247652645, 0, 0, 0},
    {-0.136216083222685, 0.198182776700951, -0.0454053610742281, -2.79860824619563, 0.69230585067761, -0.9535249798913, -0.0299136877046616, 0.105232736483552, -0.151707756592252, 1.75703046914475, -1.07638870781336, 0.293933835691236, 0, 0, 0},
    {-1.3509707503851, 2.23500889639906, -0.450323583461695, -0.126876025454617, -2.99478067797227, -0.336971390489534, -0.22931404695821, 0.908951677378124, -1.57198028688859, 1.74643464160844, -1.89234265646868, 4.19336111981495, 0, 0, 0},
    {-2.44653504112257, 4.07198390337324, -0.815511680374179, 2.28270880231473, -6.32009437862577, 0.219086646688001, -0.409149464811518, 1.63381060999723, -2.85289725668524, 1.73687846488124, -2.62823611149953, 7.71018093726218, 0, 0, 0},
    {4.58077501969351, -3.82543322867094, 1.52692500656451, -0.545139179176643, 1.41015374455001, -0.433493656733078, 1.71576045432015, -4.95844591520479, 4.39193957193787, -1.56680606271386, 0.149254111053221, 3.58290070731338, 0, 0, 0},
    {2.81409991956521, -1.76869455849098, 0.938033306521741, 0.397567423491237, -0.782852150716323, -0.215945979194339, 1.19938464679029, -3.33680260010231, 2.55274857929665, -0.798018138461319, -0.584669374015374, 4.95879302871793, 0, 0, 0},
    {0.347885561958113, 1.10243802150358, 0.115961853986046, 1.71355164625053, -3.84420883414388, 0.0877426875962636, 0.478542749851463, -1.07304735368896, -0.0146953339073127, 0.275181917739669, -1.60919993774481, 6.87948880199444, 0, 0, 0},
    {-1.87634658675205, 3.69185826712623, -0.625448862250672, 2.90041296121612, -6.60518860581103, 0.361633760280628, -0.171570942157136, 0.96859074656496, -2.3302245068456, 1.24308077845179, -2.53320470243777, 8.61172804293309, 0, 0, 0},
    {7.52914571196563, -5.79101369018569, 2.50971523732188, 2.64892907078481, -0.0640316015860407, 0.303599016334941, 2.94424824276686, -8.39821172285559, 7.09461270652064, -4.12016998199502, 0.640649226431912, 8.24468409434702, 0, 0, 0},
    {5.07536235780609, -3.27620285065158, 1.69178745260204, 2.84726839825218, -1.91348336983675, 0.349369630365873, 2.14157732939066, -5.97494211138334, 4.62557248101746, -2.75632885460146, -0.207792300975223, 8.53416287483144, 0, 0, 0},
    {1.64996980088034, 0.234381862222092, 0.549989933626791, 3.12414290508293, -4.49525095360498, 0.413263747326815, 1.02107784940239, -2.59214563243155, 1.17888188510473, -0.852456111034303, -1.39218589792443, 8.93826475219498, 0, 0, 0},
    {-1.43932704939798, 3.40051190889018, -0.479775683132647, 3.3738507933497, -6.82369837448806, 0.470888644619145, 0.0105205317403942, 0.458734619651876, -1.92962326427103, 0.864610757153044, -2.46036811287876, 9.30271660294904, 0, 0, 0}};
    
    static const double FE1_C1[16][15] = \
    {{-0.0505824176867465, 0.0590143996433656, -0.15174725306024, 0.0737656310409652, -0.198113679768677, 0.128564039706022, -1.90692311150999, 1.74059857743675, -0.51483161550127, 0.080309006824432, 0.261739654310967, 0.439870546085776, 0, 0, 0},
    {-0.123558905237645, 0.0258760841370745, -0.370676715712938, 0.209759063131611, -0.543077031500869, 0.284476557818972, -0.566936126133585, 1.00044751850615, -0.183464225194967, 0.298606106185617, 0.973204153539475, 1.63553300166645, 0, 0, 0},
    {0.097565153136164, -0.711229475708454, 0.292695459408488, 0.00706835093882763, 0.0834284512585069, -0.39732896348348, 0.393440423411939, 0.307191394545831, 0.0642761810266409, 0.3624610904637, 1.18131756661549, 1.98528115464889, 0, 0, 0},
    {0.619786046331443, -2.06642188658317, 1.85935813899433, -0.568020109434232, 1.75582626519991, -1.91112407589154, 0.349962130335548, 0.0217107082214378, 0.0730073495432117, 0.179322465614749, 0.584440052443978, 0.982189594846425, 0, 0, 0},
    {-0.0505824176867469, 0.115135845719967, -0.151747253060241, 0.0597352695218153, -0.170052956730378, 0.142594401225173, -0.372641843944438, -0.513538065166794, 0.176963036498415, 0.0317064088249769, 0.233678931272666, 0.173663153309207, 0, 0, 0},
    {-0.123558905237646, 0.234547653480176, -0.37067671571294, 0.157591170795837, -0.438741246829321, 0.336644450154749, 0.20380962815237, -0.413304552014765, 0.355206306368442, 0.117891226086926, 0.868868368867924, 0.645716838597509, 0, 0, 0},
    {0.0975651531361626, -0.45793484105066, 0.292695459408486, -0.0562553077256191, 0.2100757685874, -0.33400530481903, 0.510918247340365, -0.253620077643653, 0.380962511958802, 0.143101502207754, 1.0546702492866, 0.783799207720622, 0, 0, 0},
    {0.619786046331443, -1.94110797121721, 1.85935813899433, -0.599348588275721, 1.81848322288289, -1.87979559705005, 0.290610190718374, -0.0898539795240764, 0.181267019222918, 0.0707974314601332, 0.521783094760997, 0.387773502241409, 0, 0, 0},
    {-0.0505824176867473, 0.188358660400515, -0.151747253060241, 0.0414295658516789, -0.133441549390105, 0.16090010489531, 0.176963036498415, -0.550149472507067, -0.372641843944438, -0.0317064088249769, 0.197067523932393, -0.173663153309207, 0, 0, 0},
    {-0.123558905237648, 0.506805777945706, -0.370676715712942, 0.0895266396794566, -0.30261218459656, 0.404708981271133, 0.355206306368442, -0.549433614247528, 0.203809628152369, -0.117891226086926, 0.732739306635159, -0.645716838597509, 0, 0, 0},
    {0.0975651531361609, -0.12745607776631, 0.292695459408484, -0.138874998546704, 0.37531515022957, -0.25138561399794, 0.380962511958801, -0.418859459285826, 0.510918247340364, -0.143101502207754, 0.889430867644423, -0.783799207720622, 0, 0, 0},
    {0.619786046331442, -1.77760830677144, 1.85935813899433, -0.640223504387163, 1.90023305510577, -1.83892068093861, 0.181267019222917, -0.17160381174696, 0.290610190718374, -0.0707974314601336, 0.440033262538112, -0.387773502241409, 0, 0, 0},
    {-0.0505824176867476, 0.244480106477117, -0.151747253060242, 0.0273992043325288, -0.105380826351805, 0.174930466414461, -0.51483161550127, 1.64786572401987, -1.90692311150999, -0.080309006824432, 0.169006800894092, -0.439870546085776, 0, 0, 0},
    {-0.123558905237649, 0.715477347288807, -0.370676715712943, 0.0373587473436828, -0.198276399925013, 0.45687687360691, -0.183464225194968, 0.655646886930288, -0.566936126133587, -0.298606106185616, 0.628403521963608, -1.63553300166645, 0, 0, 0},
    {0.0975651531361595, 0.125838556891484, 0.292695459408482, -0.202198657211151, 0.501962467558464, -0.18806195533349, 0.0642761810266382, -0.111342621754132, 0.393440423411937, -0.3624610904637, 0.762783550315526, -1.98528115464889, 0, 0, 0},
    {0.619786046331441, -1.65229439140548, 1.85935813899433, -0.671551983228653, 1.96289001278875, -1.80759220209711, 0.07300734954321, -0.185353039367407, 0.349962130335547, -0.17932246561475, 0.37737630485513, -0.982189594846426, 0, 0, 0}};
    
    static const double FE1_C1_D01[16][15] = \
    {{-0.771583215541922, 0.805034409303901, -2.31474964662578, 1.14901202487239, -3.06960726528671, 1.93732083729531, 7.90657452129608, -4.74083689954269, 2.13228642799141, 1.3074517480012, 4.29802404974481, 7.16120815229814, 0, 0, 0},
    {0.107372054552499, -1.01390361625952, 0.32211616365749, 0.0655748085980085, -0.0237775626435214, -0.495063026807999, 5.24830922095658, -3.63600019321208, 1.51884058945152, 0.599105507972679, 2.13114961719603, 3.28143601042227, 0, 0, 0},
    {1.33436172947566, -3.55307512205377, 4.003085188427, -1.44686424606897, 4.22809022161361, -3.89058267183369, 1.53746728428702, -2.09368822664505, 0.66249245022008, -0.389720149437935, -0.89372849213795, -2.1345851696144, 0, 0, 0},
    {2.44096054262487, -5.84310594474195, 7.32288162787462, -2.81090446340804, 8.06276946944095, -6.95293770709146, -1.80927097018116, -0.702706194285433, -0.109831762349505, -1.28152333349536, -3.6218089268161, -7.01919237724623, 0, 0, 0},
    {-0.771583215541927, 1.78782464006128, -2.31474964662579, 0.903314467183053, -2.57821214990804, 2.18301839498466, 4.71250627133463, -1.3010710918919, 1.39519375492338, 0.456330441574142, 3.80662893436611, 2.4994247652645, 0, 0, 0},
    {0.107372054552495, -0.260149470179218, 0.322116163657485, -0.122863727922062, 0.353099510396618, -0.306624490287918, 2.79860824619563, -0.997860681931056, 0.9535249798913, -0.0536647307407026, 1.75427254415588, -0.293933835691237, 0, 0, 0},
    {1.33436172947566, -3.11904704241303, 4.003085188427, -1.55537126597915, 4.44510426143398, -3.7820756519235, 0.126876025454613, -0.574589947902456, 0.336971390489522, -0.765599492362593, -1.11074253195832, -4.19336111981495, 0, 0, 0},
    {2.44096054262487, -5.69743276562393, 7.32288162787462, -2.84732275818754, 8.13560605899996, -6.91651941231195, -2.28270880231474, -0.19285006737235, -0.219086646688023, -1.40768000726161, -3.69464551637511, -7.71018093726218, 0, 0, 0},
    {-0.771583215541934, 3.07009143764837, -2.31474964662579, 0.582747767786289, -1.93707875111451, 2.50358509438144, 0.54513917917664, 3.18686269966286, 0.433493656733067, -0.654145179565585, 3.16549553557257, -3.58290070731338, 0, 0, 0},
    {0.107372054552489, 0.723289197416755, 0.322116163657478, -0.368723394821048, 0.84481884419459, -0.0607648233889185, -0.397567423491243, 2.44417465465481, 0.215945979194323, -0.905347599942728, 1.26255321035789, -4.95879302871793, 0, 0, 0},
    {1.33436172947566, -2.55276160496528, 4.00308518842699, -1.69694262534109, 4.72824698015784, -3.64050429256156, -1.71355164625053, 1.40740908316463, -0.0877426875962843, -1.25601706698552, -1.3938852506822, -6.87948880199444, 0, 0, 0},
    {2.44096054262487, -5.50736994750042, 7.32288162787462, -2.89483846271842, 8.23063746806171, -6.86900370778107, -2.90041296121613, 0.47236979605992, -0.361633760280654, -1.57227923607143, -3.78967692543686, -8.61172804293309, 0, 0, 0},
    {-0.771583215541939, 4.05288166840575, -2.3147496466258, 0.337050210096951, -1.44568363573583, 2.74928265207079, -2.64892907078482, 6.62662850731366, -0.303599016334966, -1.50526648599264, 2.67410042019388, -8.24468409434702, 0, 0, 0},
    {0.107372054552485, 1.47704334349706, 0.322116163657472, -0.557161931341118, 1.22169591723473, 0.127673713131162, -2.84726839825219, 5.08231416593583, -0.349369630365899, -1.55811783865611, 0.885676137317739, -8.53416287483143, 0, 0, 0},
    {1.33436172947566, -2.11873352532454, 4.00308518842699, -1.80544964525127, 4.94526101997821, -3.53199727265137, -3.12414290508294, 2.92650736190722, -0.413263747326842, -1.63189640991018, -1.61089929050257, -8.93826475219498, 0, 0, 0},
    {2.44096054262487, -5.36169676838239, 7.32288162787462, -2.93125675749792, 8.30347405762072, -6.83258541300156, -3.37385079334971, 0.982225922973004, -0.470888644619172, -1.69843590983768, -3.86251351499588, -9.30271660294904, 0, 0, 0}};
    
    static const double FE1_C1_D10[16][15] = \
    {{0, 0.228416784458073, 0, -0.0571041961145164, 0.114208392229032, 0.0571041961145201, 8.21017353763104, -13.1055974286363, 4.78121549877622, -0.197814737991444, -0.114208392229036, -1.08347594204888, 0, 0, 0},
    {0, 1.10737205455251, 0, -0.276843013638119, 0.553686027276237, 0.276843013638135, 5.59767885132247, -10.5174738663024, 4.36610898770371, -0.959012330683431, -0.553686027276254, -5.25272686440917, 0, 0, 0},
    {-1.26794480176767e-14, 2.33436172947569, -1.74830731879172e-14, -0.583590432368907, 1.16718086473781, 0.583590432368939, 1.95073103161385, -6.90454725165471, 3.78663535530301, -2.02161655934811, -1.16718086473785, -11.0728499218094, 0, 0, 0},
    {-1.85496645174957e-14, 3.44096054262492, -2.53739577937815e-14, -0.860240135656206, 1.72048027131241, 0.860240135656253, -1.338382325562, -3.64611697675064, 3.26401903100018, -2.97995924333304, -1.72048027131246, -16.3219089801953, 0, 0, 0},
    {0, 0.228416784458073, 0, -0.0571041961145166, 0.114208392229033, 0.0571041961145198, 4.27901261460156, -5.24327558257733, 0.850054575746733, -0.197814737991444, -0.114208392229036, -1.08347594204888, 0, 0, 0},
    {0, 1.10737205455251, 0, -0.276843013638119, 0.553686027276237, 0.276843013638134, 2.5826622670013, -4.48744069766009, 1.35109240338254, -0.959012330683431, -0.553686027276254, -5.25272686440917, 0, 0, 0},
    {-1.26217224721478e-14, 2.33436172947569, -1.73098965513305e-14, -0.583590432368907, 1.16718086473781, 0.583590432368939, 0.214618713050888, -3.43232261452878, 2.05052303674005, -2.02161655934811, -1.16718086473785, -11.0728499218094, 0, 0, 0},
    {-1.85302900481285e-14, 3.44096054262492, -2.53158343856799e-14, -0.860240135656206, 1.72048027131241, 0.860240135656252, -1.92107504203409, -2.48073154380645, 2.68132631452809, -2.97995924333304, -1.72048027131246, -16.3219089801953, 0, 0, 0},
    {0, 0.228416784458073, 0, -0.0571041961145168, 0.114208392229033, 0.0571041961145193, -0.850054575746736, 5.01485879811926, -4.27901261460156, -0.197814737991444, -0.114208392229037, -1.08347594204888, 0, 0, 0},
    {0, 1.10737205455251, 0, -0.276843013638119, 0.553686027276238, 0.276843013638134, -1.35109240338254, 3.3800686431076, -2.58266226700131, -0.959012330683431, -0.553686027276255, -5.25272686440917, 0, 0, 0},
    {-1.25464067556422e-14, 2.33436172947569, -1.70839494018134e-14, -0.583590432368907, 1.16718086473781, 0.583590432368939, -2.05052303674006, 1.09796088505312, -0.214618713050903, -2.02161655934811, -1.16718086473785, -11.0728499218094, 0, 0, 0},
    {-1.85050117761097e-14, 3.44096054262492, -2.52399995696236e-14, -0.860240135656206, 1.72048027131241, 0.860240135656252, -2.68132631452812, -0.960228998818408, 1.92107504203407, -2.97995924333304, -1.72048027131246, -16.3219089801953, 0, 0, 0},
    {0, 0.228416784458073, 0, -0.057104196114517, 0.114208392229034, 0.057104196114519, -4.78121549877622, 12.8771806441782, -8.21017353763104, -0.197814737991444, -0.114208392229037, -1.08347594204888, 0, 0, 0},
    {0, 1.10737205455251, 0, -0.276843013638119, 0.553686027276238, 0.276843013638134, -4.36610898770371, 9.41010181174993, -5.59767885132248, -0.959012330683431, -0.553686027276255, -5.25272686440917, 0, 0, 0},
    {-1.24886812101133e-14, 2.33436172947569, -1.69107727652268e-14, -0.583590432368907, 1.16718086473781, 0.583590432368939, -3.78663535530303, 4.57018552217905, -1.95073103161387, -2.02161655934811, -1.16718086473785, -11.0728499218094, 0, 0, 0},
    {-1.84856373067425e-14, 3.44096054262492, -2.51818761615219e-14, -0.860240135656206, 1.72048027131241, 0.860240135656252, -3.26401903100021, 0.205156434125782, 1.33838232556197, -2.97995924333304, -1.72048027131246, -16.3219089801953, 0, 0, 0}};
    
    static const double FE1_C2[16][15] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.877428809330468, 0.0654669945550145, 0.0571041961145176},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.672946863150506, 0.0502101232113698, 0.276843013638124},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.387497483406694, 0.0289120842243891, 0.583590432368917},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.130056079216834, 0.00970378512694621, 0.860240135656219},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.631731251641125, 0.311164552244357, 0.0571041961145176},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.484508326630433, 0.238648659731443, 0.276843013638124},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.278990463496509, 0.137419104134574, 0.583590432368917},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0936377844373285, 0.0461220799064521, 0.860240135656219},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.311164552244357, 0.631731251641125, 0.0571041961145176},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.238648659731443, 0.484508326630433, 0.276843013638124},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.137419104134574, 0.278990463496509, 0.583590432368917},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.046122079906452, 0.0936377844373286, 0.860240135656219},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0654669945550145, 0.877428809330468, 0.0571041961145176},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0502101232113699, 0.672946863150506, 0.276843013638124},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.028912084224389, 0.387497483406694, 0.583590432368917},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00970378512694601, 0.130056079216835, 0.860240135656219}};
    
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', False), ('ignore zero tables', False), ('non zero columns', False), ('remove zero terms', False), ('ignore ones', False)
    // Total number of operations to compute element tensor: 342192
    
    // Loop quadrature points for integral
    // Number of operations to compute element tensor for following IP loop = 342192
    for (unsigned int ip = 0; ip < 16; ip++)
    {
      
      // Function declarations
      double F0 = 0;
      double F1 = 0;
      
      // Total number of operations to compute function values = 12
      for (unsigned int r = 0; r < 3; r++)
      {
        F0 += FE0[ip][r]*w[0][r];
        F1 += FE0[ip][r]*w[1][r];
      }// end loop over 'r'
      
      // Number of operations for primary indices: 21375
      for (unsigned int j = 0; j < 15; j++)
      {
        for (unsigned int k = 0; k < 15; k++)
        {
          // Number of operations to compute entry: 95
          A[j*15 + k] += (((Jinv_00*(1.0/detJ)*J_00*FE1_C0_D10[ip][k] + Jinv_00*(1.0/detJ)*J_01*FE1_C1_D10[ip][k] + Jinv_10*(1.0/detJ)*J_00*FE1_C0_D01[ip][k] + Jinv_10*(1.0/detJ)*J_01*FE1_C1_D01[ip][k]) + (Jinv_01*(1.0/detJ)*J_10*FE1_C0_D10[ip][k] + Jinv_01*(1.0/detJ)*J_11*FE1_C1_D10[ip][k] + Jinv_11*(1.0/detJ)*J_10*FE1_C0_D01[ip][k] + Jinv_11*(1.0/detJ)*J_11*FE1_C1_D01[ip][k]))*FE1_C2[ip][j] + (((Jinv_01*(1.0/detJ)*J_10*FE1_C0_D10[ip][j] + Jinv_01*(1.0/detJ)*J_11*FE1_C1_D10[ip][j] + Jinv_11*(1.0/detJ)*J_10*FE1_C0_D01[ip][j] + Jinv_11*(1.0/detJ)*J_11*FE1_C1_D01[ip][j]) + (Jinv_00*(1.0/detJ)*J_00*FE1_C0_D10[ip][j] + Jinv_00*(1.0/detJ)*J_01*FE1_C1_D10[ip][j] + Jinv_10*(1.0/detJ)*J_00*FE1_C0_D01[ip][j] + Jinv_10*(1.0/detJ)*J_01*FE1_C1_D01[ip][j]))*FE1_C2[ip][k]*-1 + (((1.0/detJ)*J_00*FE1_C0[ip][j] + (1.0/detJ)*J_01*FE1_C1[ip][j])*((1.0/detJ)*J_00*FE1_C0[ip][k] + (1.0/detJ)*J_01*FE1_C1[ip][k]) + ((1.0/detJ)*J_10*FE1_C0[ip][j] + (1.0/detJ)*J_11*FE1_C1[ip][j])*((1.0/detJ)*J_10*FE1_C0[ip][k] + (1.0/detJ)*J_11*FE1_C1[ip][k]))*1/(F0)*1/(F1)))*W16[ip]*det;
        }// end loop over 'k'
      }// end loop over 'j'
    }// end loop over 'ip'
}

/// Constructor
darcyflow_0_cell_integral_0::darcyflow_0_cell_integral_0() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_0_cell_integral_0::~darcyflow_0_cell_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void darcyflow_0_cell_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c) const
{
    // Reset values of the element tensor block
    for (unsigned int j = 0; j < 225; j++)
      A[j] = 0;
    
    // Add all contributions to element tensor
    integral_0_quadrature.tabulate_tensor(A, w, c);
}

/// Constructor
darcyflow_0_exterior_facet_integral_0::darcyflow_0_exterior_facet_integral_0() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_0_exterior_facet_integral_0::~darcyflow_0_exterior_facet_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void darcyflow_0_exterior_facet_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c,
                                    unsigned int facet) const
{
    // Reset values of the element tensor block
    for (unsigned int j = 0; j < 225; j++)
      A[j] = 0;
    
    // Add all contributions to element tensor
}

/// Constructor
darcyflow_0_exterior_facet_integral_1::darcyflow_0_exterior_facet_integral_1() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_0_exterior_facet_integral_1::~darcyflow_0_exterior_facet_integral_1()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void darcyflow_0_exterior_facet_integral_1::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c,
                                    unsigned int facet) const
{
    // Reset values of the element tensor block
    for (unsigned int j = 0; j < 225; j++)
      A[j] = 0;
    
    // Add all contributions to element tensor
}

/// Constructor
darcyflow_0_exterior_facet_integral_2_quadrature::darcyflow_0_exterior_facet_integral_2_quadrature() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_0_exterior_facet_integral_2_quadrature::~darcyflow_0_exterior_facet_integral_2_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void darcyflow_0_exterior_facet_integral_2_quadrature::tabulate_tensor(double* A,
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
    static const double W3[3] = {0.277777777777778, 0.444444444444444, 0.277777777777778};
    // Quadrature points on the UFC reference element: (0.112701665379258), (0.5), (0.887298334620742)
    
    // Value of basis functions at quadrature points.
    static const double FE0_f0_C0[3][15] = \
    {{2.06189500386223, -1.86189500386223, 0.687298334620743, -0.0372983346207422, 0.161895003862227, -0.311895003862227, 0.737298334620742, -2.16189500386223, 2.01189500386223, 0.173205080756888, -0.299999999999999, 0.948683298050514, 0, 0, 0},
    {0, 0.499999999999994, 0, 0.124999999999998, -0.249999999999996, -0.125000000000005, 0.125000000000001, -0.250000000000001, -0.124999999999997, 0.433012701892221, -0.749999999999997, 2.37170824512628, 0, 0, 0},
    {-0.261895003862224, 0.461895003862222, -0.0872983346207398, 0.737298334620741, -2.16189500386222, 2.01189500386222, -0.0372983346207414, 0.161895003862225, -0.311895003862224, 0.173205080756888, -0.299999999999999, 0.948683298050513, 0, 0, 0}};
    
    static const double FE0_f0_C1[3][15] = \
    {{-0.0872983346207427, 0.461895003862227, -0.261895003862226, 0.0372983346207424, -0.161895003862227, 0.311895003862227, -0.737298334620742, 2.16189500386223, -2.01189500386223, -0.173205080756888, 0.299999999999999, -0.948683298050514, 0, 0, 0},
    {0, 0.500000000000006, 0, -0.124999999999998, 0.249999999999996, 0.125000000000005, -0.125000000000001, 0.250000000000001, 0.124999999999997, -0.433012701892221, 0.749999999999997, -2.37170824512628, 0, 0, 0},
    {0.687298334620741, -1.86189500386222, 2.06189500386222, -0.737298334620741, 2.16189500386222, -2.01189500386222, 0.0372983346207413, -0.161895003862225, 0.311895003862224, -0.173205080756888, 0.299999999999999, -0.948683298050513, 0, 0, 0}};
    
    static const double FE0_f0_C2[3][15] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.887298334620742, 0.112701665379258},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.112701665379258, 0.887298334620742}};
    
    static const double FE0_f1_C0[3][15] = \
    {{0, 0, 0, 1.97459666924148, -1.4, 0.425403330758517, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0.425403330758517, -1.4, 1.97459666924148, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    static const double FE0_f1_C1[3][15] = \
    {{-0.0872983346207405, 0.0618950038622225, -0.261895003862223, 0.137298334620741, -0.361895003862223, 0.211895003862223, -2.01189500386222, 2.36189500386222, -0.737298334620741, 0.173205080756888, 0.500000000000001, 0.948683298050514, 0, 0, 0},
    {0, -0.500000000000006, 0, 0.124999999999998, -0.249999999999995, -0.125000000000005, 0.125000000000001, 0.749999999999999, -0.124999999999998, 0.433012701892221, 1.25, 2.37170824512628, 0, 0, 0},
    {0.687298334620743, -2.26189500386223, 2.06189500386223, -0.637298334620742, 1.96189500386223, -2.11189500386223, 0.311895003862225, 0.0381049961377748, 0.0372983346207426, 0.173205080756888, 0.500000000000001, 0.948683298050514, 0, 0, 0}};
    
    static const double FE0_f1_C2[3][15] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.887298334620742, 0, 0.112701665379258},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.5},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.112701665379258, 0, 0.887298334620742}};
    
    static const double FE0_f2_C0[3][15] = \
    {{-0.261895003862226, 0.0618950038622272, -0.0872983346207427, 2.01189500386223, -2.36189500386223, 0.737298334620744, -0.137298334620742, 0.361895003862225, -0.211895003862226, 0.519615242270663, -0.100000000000001, -0.948683298050514, 0, 0, 0},
    {0, -0.499999999999994, 0, -0.124999999999998, -0.750000000000004, 0.125000000000005, -0.125000000000001, 0.250000000000001, 0.124999999999998, 1.29903810567666, -0.250000000000003, -2.37170824512628, 0, 0, 0},
    {2.06189500386222, -2.26189500386222, 0.68729833462074, -0.311895003862224, -0.0381049961377768, -0.0372983346207398, 0.637298334620741, -1.96189500386223, 2.11189500386222, 0.519615242270663, -0.100000000000001, -0.948683298050513, 0, 0, 0}};
    
    static const double FE0_f2_C1[3][15] = \
    {{0, 0, 0, 0, 0, 0, -1.97459666924148, 1.4, -0.425403330758517, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, -0.425403330758517, 1.4, -1.97459666924148, 0, 0, 0, 0, 0, 0}};
    
    static const double FE0_f2_C2[3][15] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.887298334620742, 0.112701665379258, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.112701665379258, 0.887298334620742, 0}};
    
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', False), ('ignore zero tables', False), ('non zero columns', False), ('remove zero terms', False), ('ignore ones', False)
    switch ( facet )
    {
    case 0:
      {
      // Total number of operations to compute element tensor (from this point): 22275
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 22275
      for (unsigned int ip = 0; ip < 3; ip++)
      {
        
        // Number of operations for primary indices: 7425
        for (unsigned int j = 0; j < 15; j++)
        {
          for (unsigned int k = 0; k < 15; k++)
          {
            // Number of operations to compute entry: 33
            A[j*15 + k] += ((((1.0/detJ)*J_00*FE0_f0_C0[ip][k] + (1.0/detJ)*J_01*FE0_f0_C1[ip][k])*n0 + ((1.0/detJ)*J_10*FE0_f0_C0[ip][k] + (1.0/detJ)*J_11*FE0_f0_C1[ip][k])*n1)*FE0_f0_C2[ip][j]*-1 + (((1.0/detJ)*J_10*FE0_f0_C0[ip][j] + (1.0/detJ)*J_11*FE0_f0_C1[ip][j])*n1 + ((1.0/detJ)*J_00*FE0_f0_C0[ip][j] + (1.0/detJ)*J_01*FE0_f0_C1[ip][j])*n0)*FE0_f0_C2[ip][k])*W3[ip]*det;
          }// end loop over 'k'
        }// end loop over 'j'
      }// end loop over 'ip'
      }
      break;
    case 1:
      {
      // Total number of operations to compute element tensor (from this point): 22275
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 22275
      for (unsigned int ip = 0; ip < 3; ip++)
      {
        
        // Number of operations for primary indices: 7425
        for (unsigned int j = 0; j < 15; j++)
        {
          for (unsigned int k = 0; k < 15; k++)
          {
            // Number of operations to compute entry: 33
            A[j*15 + k] += ((((1.0/detJ)*J_00*FE0_f1_C0[ip][k] + (1.0/detJ)*J_01*FE0_f1_C1[ip][k])*n0 + ((1.0/detJ)*J_10*FE0_f1_C0[ip][k] + (1.0/detJ)*J_11*FE0_f1_C1[ip][k])*n1)*FE0_f1_C2[ip][j]*-1 + (((1.0/detJ)*J_00*FE0_f1_C0[ip][j] + (1.0/detJ)*J_01*FE0_f1_C1[ip][j])*n0 + ((1.0/detJ)*J_10*FE0_f1_C0[ip][j] + (1.0/detJ)*J_11*FE0_f1_C1[ip][j])*n1)*FE0_f1_C2[ip][k])*W3[ip]*det;
          }// end loop over 'k'
        }// end loop over 'j'
      }// end loop over 'ip'
      }
      break;
    case 2:
      {
      // Total number of operations to compute element tensor (from this point): 22275
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 22275
      for (unsigned int ip = 0; ip < 3; ip++)
      {
        
        // Number of operations for primary indices: 7425
        for (unsigned int j = 0; j < 15; j++)
        {
          for (unsigned int k = 0; k < 15; k++)
          {
            // Number of operations to compute entry: 33
            A[j*15 + k] += ((((1.0/detJ)*J_10*FE0_f2_C0[ip][j] + (1.0/detJ)*J_11*FE0_f2_C1[ip][j])*n1 + ((1.0/detJ)*J_00*FE0_f2_C0[ip][j] + (1.0/detJ)*J_01*FE0_f2_C1[ip][j])*n0)*FE0_f2_C2[ip][k] + (((1.0/detJ)*J_00*FE0_f2_C0[ip][k] + (1.0/detJ)*J_01*FE0_f2_C1[ip][k])*n0 + ((1.0/detJ)*J_10*FE0_f2_C0[ip][k] + (1.0/detJ)*J_11*FE0_f2_C1[ip][k])*n1)*FE0_f2_C2[ip][j]*-1)*W3[ip]*det;
          }// end loop over 'k'
        }// end loop over 'j'
      }// end loop over 'ip'
      }
      break;
    }
}

/// Constructor
darcyflow_0_exterior_facet_integral_2::darcyflow_0_exterior_facet_integral_2() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_0_exterior_facet_integral_2::~darcyflow_0_exterior_facet_integral_2()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void darcyflow_0_exterior_facet_integral_2::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c,
                                    unsigned int facet) const
{
    // Reset values of the element tensor block
    for (unsigned int j = 0; j < 225; j++)
      A[j] = 0;
    
    // Add all contributions to element tensor
    integral_2_quadrature.tabulate_tensor(A, w, c, facet);
}

/// Constructor
darcyflow_form_0::darcyflow_form_0() : ufc::form()
{
    // Do nothing
}

/// Destructor
darcyflow_form_0::~darcyflow_form_0()
{
    // Do nothing
}

/// Return a string identifying the form
const char* darcyflow_form_0::signature() const
{
    return "Form([Integral(Sum(Product(IndexSum(Indexed(ListTensor(Indexed(SpatialDerivative(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 1), MultiIndex((Index(0),), {Index(0): 2})), MultiIndex((FixedIndex(0),), {})), Indexed(SpatialDerivative(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 1), MultiIndex((Index(0),), {Index(0): 2})), MultiIndex((FixedIndex(1),), {}))), MultiIndex((Index(0),), {Index(0): 2})), MultiIndex((Index(0),), {Index(0): 2})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(2),), {FixedIndex(2): 3}))), Sum(Product(IndexSum(Product(Indexed(ListTensor(Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(1),), {Index(1): 2})), Indexed(ListTensor(Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 1), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 1), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(1),), {Index(1): 2}))), MultiIndex((Index(1),), {Index(1): 2})), Product(Division(FloatValue(1.0, (), (), {}), Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1), 0)), Division(FloatValue(1.0, (), (), {}), Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1), 1)))), Product(IntValue(-1, (), (), {}), Product(IndexSum(Indexed(ListTensor(Indexed(SpatialDerivative(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((Index(2),), {Index(2): 2})), MultiIndex((FixedIndex(0),), {})), Indexed(SpatialDerivative(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((Index(2),), {Index(2): 2})), MultiIndex((FixedIndex(1),), {}))), MultiIndex((Index(2),), {Index(2): 2})), MultiIndex((Index(2),), {Index(2): 2})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 1), MultiIndex((FixedIndex(2),), {FixedIndex(2): 3})))))), Measure('cell', 0, None)), Integral(Sum(Product(IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(3),), {Index(3): 2})), Indexed(ListTensor(Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(3),), {Index(3): 2}))), MultiIndex((Index(3),), {Index(3): 2})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 1), MultiIndex((FixedIndex(2),), {FixedIndex(2): 3}))), Product(IntValue(-1, (), (), {}), Product(IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(4),), {Index(4): 2})), Indexed(ListTensor(Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 1), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 1), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(4),), {Index(4): 2}))), MultiIndex((Index(4),), {Index(4): 2})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(2),), {FixedIndex(2): 3}))))), Measure('exterior_facet', 2, None))])";
}

/// Return the rank of the global tensor (r)
unsigned int darcyflow_form_0::rank() const
{
    return 2;
}

/// Return the number of coefficients (n)
unsigned int darcyflow_form_0::num_coefficients() const
{
    return 2;
}

/// Return the number of cell integrals
unsigned int darcyflow_form_0::num_cell_integrals() const
{
    return 1;
}

/// Return the number of exterior facet integrals
unsigned int darcyflow_form_0::num_exterior_facet_integrals() const
{
    return 3;
}

/// Return the number of interior facet integrals
unsigned int darcyflow_form_0::num_interior_facet_integrals() const
{
    return 0;
}

/// Create a new finite element for argument function i
ufc::finite_element* darcyflow_form_0::create_finite_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_0_finite_element_0();
      break;
    case 1:
      return new darcyflow_0_finite_element_1();
      break;
    case 2:
      return new darcyflow_0_finite_element_2();
      break;
    case 3:
      return new darcyflow_0_finite_element_3();
      break;
    }
    return 0;
}

/// Create a new dof map for argument function i
ufc::dof_map* darcyflow_form_0::create_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_0_dof_map_0();
      break;
    case 1:
      return new darcyflow_0_dof_map_1();
      break;
    case 2:
      return new darcyflow_0_dof_map_2();
      break;
    case 3:
      return new darcyflow_0_dof_map_3();
      break;
    }
    return 0;
}

/// Create a new cell integral on sub domain i
ufc::cell_integral* darcyflow_form_0::create_cell_integral(unsigned int i) const
{
    return new darcyflow_0_cell_integral_0();
}

/// Create a new exterior facet integral on sub domain i
ufc::exterior_facet_integral* darcyflow_form_0::create_exterior_facet_integral(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_0_exterior_facet_integral_0();
      break;
    case 1:
      return new darcyflow_0_exterior_facet_integral_1();
      break;
    case 2:
      return new darcyflow_0_exterior_facet_integral_2();
      break;
    }
    return 0;
}

/// Create a new interior facet integral on sub domain i
ufc::interior_facet_integral* darcyflow_form_0::create_interior_facet_integral(unsigned int i) const
{
    return 0;
}


/// Constructor
darcyflow_1_finite_element_0_0::darcyflow_1_finite_element_0_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_1_finite_element_0_0::~darcyflow_1_finite_element_0_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_1_finite_element_0_0::signature() const
{
    return "FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2)";
}

/// Return the cell shape
ufc::shape darcyflow_1_finite_element_0_0::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_1_finite_element_0_0::space_dimension() const
{
    return 12;
}

/// Return the rank of the value space
unsigned int darcyflow_1_finite_element_0_0::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_1_finite_element_0_0::value_dimension(unsigned int i) const
{
    return 2;
}

/// Evaluate basis function i at given point in cell
void darcyflow_1_finite_element_0_0::evaluate_basis(unsigned int i,
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
    const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
    const double psitilde_bs_1_0 = 1;
    const double psitilde_bs_1_1 = 2.5*y + 1.5;
    const double psitilde_bs_2_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
    const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
    const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
    // Table(s) of coefficients
    static const double coefficients0[12][6] = \
    {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
    {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
    {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
    {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
    {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
    {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
    {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
    {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
    {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
    {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
    {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
    {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
    static const double coefficients1[12][6] = \
    {{0, 0, 0.2, 0, 0, 0.163299316185545},
    {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
    {0, 0, 0.6, 0, 0, 0.489897948556636},
    {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
    {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
    {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
    {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
    {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
    {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
    {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
    {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
    {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff0_3 = coefficients0[dof][3];
    const double coeff0_4 = coefficients0[dof][4];
    const double coeff0_5 = coefficients0[dof][5];
    const double coeff1_0 = coefficients1[dof][0];
    const double coeff1_1 = coefficients1[dof][1];
    const double coeff1_2 = coefficients1[dof][2];
    const double coeff1_3 = coefficients1[dof][3];
    const double coeff1_4 = coefficients1[dof][4];
    const double coeff1_5 = coefficients1[dof][5];
    
    // Compute value(s)
    const double tmp0_0 = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3 + coeff0_4*basisvalue4 + coeff0_5*basisvalue5;
    const double tmp0_1 = coeff1_0*basisvalue0 + coeff1_1*basisvalue1 + coeff1_2*basisvalue2 + coeff1_3*basisvalue3 + coeff1_4*basisvalue4 + coeff1_5*basisvalue5;
    // Using contravariant Piola transform to map values back to the physical element
    values[0] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
    values[1] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
}

/// Evaluate all basis functions at given point in cell
void darcyflow_1_finite_element_0_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_1_finite_element_0_0::evaluate_basis_derivatives(unsigned int i,
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
    const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
    const double psitilde_bs_1_0 = 1;
    const double psitilde_bs_1_1 = 2.5*y + 1.5;
    const double psitilde_bs_2_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
    const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
    const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
    // Table(s) of coefficients
    static const double coefficients0[12][6] = \
    {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
    {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
    {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
    {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
    {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
    {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
    {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
    {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
    {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
    {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
    {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
    {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
    static const double coefficients1[12][6] = \
    {{0, 0, 0.2, 0, 0, 0.163299316185545},
    {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
    {0, 0, 0.6, 0, 0, 0.489897948556636},
    {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
    {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
    {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
    {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
    {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
    {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
    {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
    {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
    {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[6][6] = \
    {{0, 0, 0, 0, 0, 0},
    {4.89897948556635, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 9.48683298050514, 0, 0, 0, 0},
    {4, 0, 7.07106781186548, 0, 0, 0},
    {0, 0, 0, 0, 0, 0}};
    
    static const double dmats1[6][6] = \
    {{0, 0, 0, 0, 0, 0},
    {2.44948974278318, 0, 0, 0, 0, 0},
    {4.24264068711928, 0, 0, 0, 0, 0},
    {2.58198889747161, 4.74341649025257, -0.912870929175277, 0, 0, 0},
    {2, 6.12372435695795, 3.53553390593274, 0, 0, 0},
    {-2.3094010767585, 0, 8.16496580927726, 0, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [2*num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    double coeff0_3 = 0;
    double coeff0_4 = 0;
    double coeff0_5 = 0;
    double coeff1_0 = 0;
    double coeff1_1 = 0;
    double coeff1_2 = 0;
    double coeff1_3 = 0;
    double coeff1_4 = 0;
    double coeff1_5 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff0_3 = 0;
    double new_coeff0_4 = 0;
    double new_coeff0_5 = 0;
    double new_coeff1_0 = 0;
    double new_coeff1_1 = 0;
    double new_coeff1_2 = 0;
    double new_coeff1_3 = 0;
    double new_coeff1_4 = 0;
    double new_coeff1_5 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
      new_coeff0_3 = coefficients0[dof][3];
      new_coeff0_4 = coefficients0[dof][4];
      new_coeff0_5 = coefficients0[dof][5];
      new_coeff1_0 = coefficients1[dof][0];
      new_coeff1_1 = coefficients1[dof][1];
      new_coeff1_2 = coefficients1[dof][2];
      new_coeff1_3 = coefficients1[dof][3];
      new_coeff1_4 = coefficients1[dof][4];
      new_coeff1_5 = coefficients1[dof][5];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
        coeff0_3 = new_coeff0_3;
        coeff0_4 = new_coeff0_4;
        coeff0_5 = new_coeff0_5;
        coeff1_0 = new_coeff1_0;
        coeff1_1 = new_coeff1_1;
        coeff1_2 = new_coeff1_2;
        coeff1_3 = new_coeff1_3;
        coeff1_4 = new_coeff1_4;
        coeff1_5 = new_coeff1_5;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0] + coeff0_4*dmats0[4][0] + coeff0_5*dmats0[5][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1] + coeff0_4*dmats0[4][1] + coeff0_5*dmats0[5][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2] + coeff0_4*dmats0[4][2] + coeff0_5*dmats0[5][2];
          new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3] + coeff0_4*dmats0[4][3] + coeff0_5*dmats0[5][3];
          new_coeff0_4 = coeff0_0*dmats0[0][4] + coeff0_1*dmats0[1][4] + coeff0_2*dmats0[2][4] + coeff0_3*dmats0[3][4] + coeff0_4*dmats0[4][4] + coeff0_5*dmats0[5][4];
          new_coeff0_5 = coeff0_0*dmats0[0][5] + coeff0_1*dmats0[1][5] + coeff0_2*dmats0[2][5] + coeff0_3*dmats0[3][5] + coeff0_4*dmats0[4][5] + coeff0_5*dmats0[5][5];
          new_coeff1_0 = coeff1_0*dmats0[0][0] + coeff1_1*dmats0[1][0] + coeff1_2*dmats0[2][0] + coeff1_3*dmats0[3][0] + coeff1_4*dmats0[4][0] + coeff1_5*dmats0[5][0];
          new_coeff1_1 = coeff1_0*dmats0[0][1] + coeff1_1*dmats0[1][1] + coeff1_2*dmats0[2][1] + coeff1_3*dmats0[3][1] + coeff1_4*dmats0[4][1] + coeff1_5*dmats0[5][1];
          new_coeff1_2 = coeff1_0*dmats0[0][2] + coeff1_1*dmats0[1][2] + coeff1_2*dmats0[2][2] + coeff1_3*dmats0[3][2] + coeff1_4*dmats0[4][2] + coeff1_5*dmats0[5][2];
          new_coeff1_3 = coeff1_0*dmats0[0][3] + coeff1_1*dmats0[1][3] + coeff1_2*dmats0[2][3] + coeff1_3*dmats0[3][3] + coeff1_4*dmats0[4][3] + coeff1_5*dmats0[5][3];
          new_coeff1_4 = coeff1_0*dmats0[0][4] + coeff1_1*dmats0[1][4] + coeff1_2*dmats0[2][4] + coeff1_3*dmats0[3][4] + coeff1_4*dmats0[4][4] + coeff1_5*dmats0[5][4];
          new_coeff1_5 = coeff1_0*dmats0[0][5] + coeff1_1*dmats0[1][5] + coeff1_2*dmats0[2][5] + coeff1_3*dmats0[3][5] + coeff1_4*dmats0[4][5] + coeff1_5*dmats0[5][5];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0] + coeff0_4*dmats1[4][0] + coeff0_5*dmats1[5][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1] + coeff0_4*dmats1[4][1] + coeff0_5*dmats1[5][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2] + coeff0_4*dmats1[4][2] + coeff0_5*dmats1[5][2];
          new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3] + coeff0_4*dmats1[4][3] + coeff0_5*dmats1[5][3];
          new_coeff0_4 = coeff0_0*dmats1[0][4] + coeff0_1*dmats1[1][4] + coeff0_2*dmats1[2][4] + coeff0_3*dmats1[3][4] + coeff0_4*dmats1[4][4] + coeff0_5*dmats1[5][4];
          new_coeff0_5 = coeff0_0*dmats1[0][5] + coeff0_1*dmats1[1][5] + coeff0_2*dmats1[2][5] + coeff0_3*dmats1[3][5] + coeff0_4*dmats1[4][5] + coeff0_5*dmats1[5][5];
          new_coeff1_0 = coeff1_0*dmats1[0][0] + coeff1_1*dmats1[1][0] + coeff1_2*dmats1[2][0] + coeff1_3*dmats1[3][0] + coeff1_4*dmats1[4][0] + coeff1_5*dmats1[5][0];
          new_coeff1_1 = coeff1_0*dmats1[0][1] + coeff1_1*dmats1[1][1] + coeff1_2*dmats1[2][1] + coeff1_3*dmats1[3][1] + coeff1_4*dmats1[4][1] + coeff1_5*dmats1[5][1];
          new_coeff1_2 = coeff1_0*dmats1[0][2] + coeff1_1*dmats1[1][2] + coeff1_2*dmats1[2][2] + coeff1_3*dmats1[3][2] + coeff1_4*dmats1[4][2] + coeff1_5*dmats1[5][2];
          new_coeff1_3 = coeff1_0*dmats1[0][3] + coeff1_1*dmats1[1][3] + coeff1_2*dmats1[2][3] + coeff1_3*dmats1[3][3] + coeff1_4*dmats1[4][3] + coeff1_5*dmats1[5][3];
          new_coeff1_4 = coeff1_0*dmats1[0][4] + coeff1_1*dmats1[1][4] + coeff1_2*dmats1[2][4] + coeff1_3*dmats1[3][4] + coeff1_4*dmats1[4][4] + coeff1_5*dmats1[5][4];
          new_coeff1_5 = coeff1_0*dmats1[0][5] + coeff1_1*dmats1[1][5] + coeff1_2*dmats1[2][5] + coeff1_3*dmats1[3][5] + coeff1_4*dmats1[4][5] + coeff1_5*dmats1[5][5];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      // Correct values by the contravariant Piola transform
      const double tmp0_0 = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3 + new_coeff0_4*basisvalue4 + new_coeff0_5*basisvalue5;
      const double tmp0_1 = new_coeff1_0*basisvalue0 + new_coeff1_1*basisvalue1 + new_coeff1_2*basisvalue2 + new_coeff1_3*basisvalue3 + new_coeff1_4*basisvalue4 + new_coeff1_5*basisvalue5;
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
void darcyflow_1_finite_element_0_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_1_finite_element_0_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[12][9][2] = {{{0.75, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.75, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}};
    static const double W[12][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}};
    static const double D[12][9][2] = {{{1, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, -1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, -1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, -1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.866025403784438, 0.433012701892219}, {0.866025403784439, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784439, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}, {0.866025403784438, 0.433012701892219}}, {{0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}, {0, 0.75}}, {{-0.17542966950853, 0.148329095604429}, {-0.0180753884489578, 0.444221856552778}, {0.0340120331760423, 0.307261606416198}, {-0.396781833524282, 3.42776476771047e-16}, {0.0799728864560174, -2.52114035295885e-16}, {0.171602319815636, -3.84645013221013e-15}, {-0.0271005739041022, -0.14832909560443}, {0.426146468103823, -0.44422185655278}, {0.341273639592246, -0.307261606416206}}};
    
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
    static const unsigned int ns[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 9, 9};
    for (unsigned int j = 0; j < ns[i]; j++) {
      // Evaluate basis functions for affine mapping
      const double w0 = 1.0 - X[i][j][0] - X[i][j][1];
      const double w1 = X[i][j][0];
      const double w2 = X[i][j][1];
      
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
        result += values[k]*D[i][j][k];
      // Multiply by weights
      result *= W[i][j];
    
    } // End for
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void darcyflow_1_finite_element_0_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_1_finite_element_0_0::interpolate_vertex_values(double* vertex_values,
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
    vertex_values[0] = (1.0/detJ)*(dof_values[3]*3*J_00 + dof_values[4]*(-3*J_00) + dof_values[5]*J_00 + dof_values[6]*(-3*J_01) + dof_values[7]*3*J_01 + dof_values[8]*J_01);
    vertex_values[2] = (1.0/detJ)*(dof_values[0]*3*J_00 + dof_values[1]*(-3*J_00) + dof_values[2]*J_00 + dof_values[6]*(J_00 - J_01) + dof_values[7]*(-3*J_00 + 3*J_01) + dof_values[8]*(3*J_00 - 3*J_01));
    vertex_values[4] = (1.0/detJ)*(dof_values[0]*J_01 + dof_values[1]*(-3*J_01) + dof_values[2]*3*J_01 + dof_values[3]*(J_00 - J_01) + dof_values[4]*(-3*J_00 + 3*J_01) + dof_values[5]*(3*J_00 - 3*J_01));
    vertex_values[1] = (1.0/detJ)*(dof_values[3]*3*J_10 + dof_values[4]*(-3*J_10) + dof_values[5]*J_10 + dof_values[6]*(-3*J_11) + dof_values[7]*3*J_11 + dof_values[8]*J_11);
    vertex_values[3] = (1.0/detJ)*(dof_values[0]*3*J_10 + dof_values[1]*(-3*J_10) + dof_values[2]*J_10 + dof_values[6]*(J_10 - J_11) + dof_values[7]*(-3*J_10 + 3*J_11) + dof_values[8]*(3*J_10 - 3*J_11));
    vertex_values[5] = (1.0/detJ)*(dof_values[0]*J_11 + dof_values[1]*(-3*J_11) + dof_values[2]*3*J_11 + dof_values[3]*(J_10 - J_11) + dof_values[4]*(-3*J_10 + 3*J_11) + dof_values[5]*(3*J_10 - 3*J_11));
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_1_finite_element_0_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_1_finite_element_0_0::create_sub_element(unsigned int i) const
{
    return new darcyflow_1_finite_element_0_0();
}


/// Constructor
darcyflow_1_finite_element_0_1::darcyflow_1_finite_element_0_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_1_finite_element_0_1::~darcyflow_1_finite_element_0_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_1_finite_element_0_1::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return the cell shape
ufc::shape darcyflow_1_finite_element_0_1::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_1_finite_element_0_1::space_dimension() const
{
    return 3;
}

/// Return the rank of the value space
unsigned int darcyflow_1_finite_element_0_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_1_finite_element_0_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void darcyflow_1_finite_element_0_1::evaluate_basis(unsigned int i,
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
void darcyflow_1_finite_element_0_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_1_finite_element_0_1::evaluate_basis_derivatives(unsigned int i,
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
void darcyflow_1_finite_element_0_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_1_finite_element_0_1::evaluate_dof(unsigned int i,
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
void darcyflow_1_finite_element_0_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_1_finite_element_0_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_1_finite_element_0_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_1_finite_element_0_1::create_sub_element(unsigned int i) const
{
    return new darcyflow_1_finite_element_0_1();
}


/// Constructor
darcyflow_1_finite_element_0::darcyflow_1_finite_element_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_1_finite_element_0::~darcyflow_1_finite_element_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_1_finite_element_0::signature() const
{
    return "MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) })";
}

/// Return the cell shape
ufc::shape darcyflow_1_finite_element_0::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_1_finite_element_0::space_dimension() const
{
    return 15;
}

/// Return the rank of the value space
unsigned int darcyflow_1_finite_element_0::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_1_finite_element_0::value_dimension(unsigned int i) const
{
    return 3;
}

/// Evaluate basis function i at given point in cell
void darcyflow_1_finite_element_0::evaluate_basis(unsigned int i,
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
    
    if (0 <= i && i <= 11)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
      const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
      const double psitilde_bs_1_0 = 1;
      const double psitilde_bs_1_1 = 2.5*y + 1.5;
      const double psitilde_bs_2_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
      const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
      const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
      const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
      const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
      const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
      // Table(s) of coefficients
      static const double coefficients0[12][6] =   \
      {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
      {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
      {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
      {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
      {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
      {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
      {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
      {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
      {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
      {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
      {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
      {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
      static const double coefficients1[12][6] =   \
      {{0, 0, 0.2, 0, 0, 0.163299316185545},
      {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
      {0, 0, 0.6, 0, 0, 0.489897948556636},
      {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
      {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
      {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
      {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
      {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
      {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
      {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
      {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
      {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
      const double coeff0_4 =   coefficients0[dof][4];
      const double coeff0_5 =   coefficients0[dof][5];
      const double coeff1_0 =   coefficients1[dof][0];
      const double coeff1_1 =   coefficients1[dof][1];
      const double coeff1_2 =   coefficients1[dof][2];
      const double coeff1_3 =   coefficients1[dof][3];
      const double coeff1_4 =   coefficients1[dof][4];
      const double coeff1_5 =   coefficients1[dof][5];
    
      // Compute value(s)
      const double tmp0_0 = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3 + coeff0_4*basisvalue4 + coeff0_5*basisvalue5;
      const double tmp0_1 = coeff1_0*basisvalue0 + coeff1_1*basisvalue1 + coeff1_2*basisvalue2 + coeff1_3*basisvalue3 + coeff1_4*basisvalue4 + coeff1_5*basisvalue5;
      // Using contravariant Piola transform to map values back to the physical element
      values[0] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
      values[1] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
    }
    
    if (12 <= i && i <= 14)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 12;
    
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
      static const double coefficients0[3][3] =   \
      {{0.471404520791032, -0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0, 0.333333333333333}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
    
      // Compute value(s)
      values[2] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2;
    }
    
}

/// Evaluate all basis functions at given point in cell
void darcyflow_1_finite_element_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_1_finite_element_0::evaluate_basis_derivatives(unsigned int i,
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
    
    if (0 <= i && i <= 11)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i;
    
      // Generate scalings
      const double scalings_y_0 = 1;
      const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
      const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    
      // Compute psitilde_a
      const double psitilde_a_0 = 1;
      const double psitilde_a_1 = x;
      const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
      // Compute psitilde_bs
      const double psitilde_bs_0_0 = 1;
      const double psitilde_bs_0_1 = 1.5*y + 0.5;
      const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
      const double psitilde_bs_1_0 = 1;
      const double psitilde_bs_1_1 = 2.5*y + 1.5;
      const double psitilde_bs_2_0 = 1;
    
      // Compute basisvalues
      const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
      const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
      const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
      const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
      const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
      const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
      // Table(s) of coefficients
      static const double coefficients0[12][6] =   \
      {{0, 0.519615242270664, -0.299999999999999, 0.365148371670111, -0.282842712474619, 0.163299316185545},
      {0, -0.404145188432739, 0.499999999999998, -0.243432247780075, 0.377123616632824, -0.272165526975908},
      {0, 0.173205080756889, -0.0999999999999985, 0.121716123890038, -0.0942809041582056, 0.0544331053951811},
      {0, -0.490747728811182, -0.0500000000000009, 0.395577402642619, 0.30641293851417, 0.2993820796735},
      {0, 0.230940107675851, 0, -0.182574185835054, -0.518544972870134, -0.816496580927727},
      {0, -0.20207259421637, 0.449999999999998, 0.0912870929175265, 0.0707106781186539, 0.571547606649409},
      {0, 0.202072594216369, -0.0499999999999998, 0.152145154862546, -0.0707106781186546, 0.0272165526975908},
      {0, -0.577350269189626, 0.2, -0.426006433615129, 0.235702260395516, -0.108866210790363},
      {0, 0.490747728811183, -0.349999999999999, 0.334719340697602, -0.30641293851417, 0.190515868883136},
      {0.816496580927726, 0.1, -0.288675134594813, -0.316227766016838, 0.0816496580927729, 0},
      {-0.471404520791032, -0.173205080756887, -0.0333333333333321, 0.0608580619450192, -0.141421356237309, 0.108866210790363},
      {0, 0.547722557505166, 0.948683298050514, 0.577350269189626, 0.447213595499958, -0.516397779494322}};
    
      static const double coefficients1[12][6] =   \
      {{0, 0, 0.2, 0, 0, 0.163299316185545},
      {0, 0.230940107675853, -0.6, 0, 0.188561808316415, -0.489897948556636},
      {0, 0, 0.6, 0, 0, 0.489897948556636},
      {0, -0.0577350269189617, -0.2, 0, -0.0471404520791024, -0.163299316185545},
      {0, 0.115470053837923, 0.6, 0, 0.0942809041582047, 0.489897948556636},
      {0, 0.0577350269189647, -0.6, 0, 0.047140452079105, -0.489897948556636},
      {0, 0.288675134594813, 0.4, -0.486864495560148, -0.235702260395516, -0.217732421580727},
      {0, -0.115470053837925, -0.2, 0.973728991120295, -0.0942809041582061, 0.108866210790364},
      {0, -0.288675134594814, 0.4, -0.486864495560148, 0.235702260395515, -0.217732421580727},
      {0, -0.200000000000001, 0, 0, -0.163299316185546, 0},
      {0.942809041582063, -0.115470053837927, 0.266666666666667, 0, -0.0942809041582075, -0.32659863237109},
      {0, -1.09544511501033, 0, 0, -0.894427190999916, 0}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      static const double dmats0[6][6] =   \
      {{0, 0, 0, 0, 0, 0},
      {4.89897948556635, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0},
      {0, 9.48683298050514, 0, 0, 0, 0},
      {4, 0, 7.07106781186548, 0, 0, 0},
      {0, 0, 0, 0, 0, 0}};
    
      static const double dmats1[6][6] =   \
      {{0, 0, 0, 0, 0, 0},
      {2.44948974278318, 0, 0, 0, 0, 0},
      {4.24264068711928, 0, 0, 0, 0, 0},
      {2.58198889747161, 4.74341649025257, -0.912870929175277, 0, 0, 0},
      {2, 6.12372435695795, 3.53553390593274, 0, 0, 0},
      {-2.3094010767585, 0, 8.16496580927726, 0, 0, 0}};
    
      // Compute reference derivatives
      // Declare pointer to array of derivatives on FIAT element
      double *derivatives = new double [2*num_derivatives];
    
      // Declare coefficients
      double coeff0_0 = 0;
      double coeff0_1 = 0;
      double coeff0_2 = 0;
      double coeff0_3 = 0;
      double coeff0_4 = 0;
      double coeff0_5 = 0;
      double coeff1_0 = 0;
      double coeff1_1 = 0;
      double coeff1_2 = 0;
      double coeff1_3 = 0;
      double coeff1_4 = 0;
      double coeff1_5 = 0;
    
      // Declare new coefficients
      double new_coeff0_0 = 0;
      double new_coeff0_1 = 0;
      double new_coeff0_2 = 0;
      double new_coeff0_3 = 0;
      double new_coeff0_4 = 0;
      double new_coeff0_5 = 0;
      double new_coeff1_0 = 0;
      double new_coeff1_1 = 0;
      double new_coeff1_2 = 0;
      double new_coeff1_3 = 0;
      double new_coeff1_4 = 0;
      double new_coeff1_5 = 0;
    
      // Loop possible derivatives
      for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
      {
        // Get values from coefficients array
        new_coeff0_0 = coefficients0[dof][0];
        new_coeff0_1 = coefficients0[dof][1];
        new_coeff0_2 = coefficients0[dof][2];
        new_coeff0_3 = coefficients0[dof][3];
        new_coeff0_4 = coefficients0[dof][4];
        new_coeff0_5 = coefficients0[dof][5];
        new_coeff1_0 = coefficients1[dof][0];
        new_coeff1_1 = coefficients1[dof][1];
        new_coeff1_2 = coefficients1[dof][2];
        new_coeff1_3 = coefficients1[dof][3];
        new_coeff1_4 = coefficients1[dof][4];
        new_coeff1_5 = coefficients1[dof][5];
    
        // Loop derivative order
        for (unsigned int j = 0; j < n; j++)
        {
          // Update old coefficients
          coeff0_0 = new_coeff0_0;
          coeff0_1 = new_coeff0_1;
          coeff0_2 = new_coeff0_2;
          coeff0_3 = new_coeff0_3;
          coeff0_4 = new_coeff0_4;
          coeff0_5 = new_coeff0_5;
          coeff1_0 = new_coeff1_0;
          coeff1_1 = new_coeff1_1;
          coeff1_2 = new_coeff1_2;
          coeff1_3 = new_coeff1_3;
          coeff1_4 = new_coeff1_4;
          coeff1_5 = new_coeff1_5;
    
          if(combinations[deriv_num][j] == 0)
          {
            new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0] + coeff0_4*dmats0[4][0] + coeff0_5*dmats0[5][0];
            new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1] + coeff0_4*dmats0[4][1] + coeff0_5*dmats0[5][1];
            new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2] + coeff0_4*dmats0[4][2] + coeff0_5*dmats0[5][2];
            new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3] + coeff0_4*dmats0[4][3] + coeff0_5*dmats0[5][3];
            new_coeff0_4 = coeff0_0*dmats0[0][4] + coeff0_1*dmats0[1][4] + coeff0_2*dmats0[2][4] + coeff0_3*dmats0[3][4] + coeff0_4*dmats0[4][4] + coeff0_5*dmats0[5][4];
            new_coeff0_5 = coeff0_0*dmats0[0][5] + coeff0_1*dmats0[1][5] + coeff0_2*dmats0[2][5] + coeff0_3*dmats0[3][5] + coeff0_4*dmats0[4][5] + coeff0_5*dmats0[5][5];
            new_coeff1_0 = coeff1_0*dmats0[0][0] + coeff1_1*dmats0[1][0] + coeff1_2*dmats0[2][0] + coeff1_3*dmats0[3][0] + coeff1_4*dmats0[4][0] + coeff1_5*dmats0[5][0];
            new_coeff1_1 = coeff1_0*dmats0[0][1] + coeff1_1*dmats0[1][1] + coeff1_2*dmats0[2][1] + coeff1_3*dmats0[3][1] + coeff1_4*dmats0[4][1] + coeff1_5*dmats0[5][1];
            new_coeff1_2 = coeff1_0*dmats0[0][2] + coeff1_1*dmats0[1][2] + coeff1_2*dmats0[2][2] + coeff1_3*dmats0[3][2] + coeff1_4*dmats0[4][2] + coeff1_5*dmats0[5][2];
            new_coeff1_3 = coeff1_0*dmats0[0][3] + coeff1_1*dmats0[1][3] + coeff1_2*dmats0[2][3] + coeff1_3*dmats0[3][3] + coeff1_4*dmats0[4][3] + coeff1_5*dmats0[5][3];
            new_coeff1_4 = coeff1_0*dmats0[0][4] + coeff1_1*dmats0[1][4] + coeff1_2*dmats0[2][4] + coeff1_3*dmats0[3][4] + coeff1_4*dmats0[4][4] + coeff1_5*dmats0[5][4];
            new_coeff1_5 = coeff1_0*dmats0[0][5] + coeff1_1*dmats0[1][5] + coeff1_2*dmats0[2][5] + coeff1_3*dmats0[3][5] + coeff1_4*dmats0[4][5] + coeff1_5*dmats0[5][5];
          }
          if(combinations[deriv_num][j] == 1)
          {
            new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0] + coeff0_4*dmats1[4][0] + coeff0_5*dmats1[5][0];
            new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1] + coeff0_4*dmats1[4][1] + coeff0_5*dmats1[5][1];
            new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2] + coeff0_4*dmats1[4][2] + coeff0_5*dmats1[5][2];
            new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3] + coeff0_4*dmats1[4][3] + coeff0_5*dmats1[5][3];
            new_coeff0_4 = coeff0_0*dmats1[0][4] + coeff0_1*dmats1[1][4] + coeff0_2*dmats1[2][4] + coeff0_3*dmats1[3][4] + coeff0_4*dmats1[4][4] + coeff0_5*dmats1[5][4];
            new_coeff0_5 = coeff0_0*dmats1[0][5] + coeff0_1*dmats1[1][5] + coeff0_2*dmats1[2][5] + coeff0_3*dmats1[3][5] + coeff0_4*dmats1[4][5] + coeff0_5*dmats1[5][5];
            new_coeff1_0 = coeff1_0*dmats1[0][0] + coeff1_1*dmats1[1][0] + coeff1_2*dmats1[2][0] + coeff1_3*dmats1[3][0] + coeff1_4*dmats1[4][0] + coeff1_5*dmats1[5][0];
            new_coeff1_1 = coeff1_0*dmats1[0][1] + coeff1_1*dmats1[1][1] + coeff1_2*dmats1[2][1] + coeff1_3*dmats1[3][1] + coeff1_4*dmats1[4][1] + coeff1_5*dmats1[5][1];
            new_coeff1_2 = coeff1_0*dmats1[0][2] + coeff1_1*dmats1[1][2] + coeff1_2*dmats1[2][2] + coeff1_3*dmats1[3][2] + coeff1_4*dmats1[4][2] + coeff1_5*dmats1[5][2];
            new_coeff1_3 = coeff1_0*dmats1[0][3] + coeff1_1*dmats1[1][3] + coeff1_2*dmats1[2][3] + coeff1_3*dmats1[3][3] + coeff1_4*dmats1[4][3] + coeff1_5*dmats1[5][3];
            new_coeff1_4 = coeff1_0*dmats1[0][4] + coeff1_1*dmats1[1][4] + coeff1_2*dmats1[2][4] + coeff1_3*dmats1[3][4] + coeff1_4*dmats1[4][4] + coeff1_5*dmats1[5][4];
            new_coeff1_5 = coeff1_0*dmats1[0][5] + coeff1_1*dmats1[1][5] + coeff1_2*dmats1[2][5] + coeff1_3*dmats1[3][5] + coeff1_4*dmats1[4][5] + coeff1_5*dmats1[5][5];
          }
    
        }
        // Compute derivatives on reference element as dot product of coefficients and basisvalues
        // Correct values by the contravariant Piola transform
        const double tmp0_0 = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3 + new_coeff0_4*basisvalue4 + new_coeff0_5*basisvalue5;
        const double tmp0_1 = new_coeff1_0*basisvalue0 + new_coeff1_1*basisvalue1 + new_coeff1_2*basisvalue2 + new_coeff1_3*basisvalue3 + new_coeff1_4*basisvalue4 + new_coeff1_5*basisvalue5;
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
    
    if (12 <= i && i <= 14)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 12;
    
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
      static const double coefficients0[3][3] =   \
      {{0.471404520791032, -0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0.288675134594813, -0.166666666666667},
      {0.471404520791032, 0, 0.333333333333333}};
    
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
void darcyflow_1_finite_element_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_1_finite_element_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    static const double X[15][9][2] = {{{0.75, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.25}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.5}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 0.75}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.25, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.5, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.75, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0.102717654809626, 0.088587959512704}, {0.0665540678391645, 0.409466864440735}, {0.0239311322870806, 0.787659461760847}, {0.455706020243648, 0.088587959512704}, {0.295266567779633, 0.409466864440735}, {0.106170269119576, 0.787659461760847}, {0.80869438567767, 0.088587959512704}, {0.523979067720101, 0.409466864440735}, {0.188409405952072, 0.787659461760847}}, {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}, {{0, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}}};
    static const double W[15][9] = {{1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {0.223257681932177, 0.25471234039954, 0.0775855332238378, 0.357212291091484, 0.407539744639264, 0.124136853158141, 0.223257681932177, 0.25471234039954, 0.0775855332238378}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0, 0, 0, 0}};
    static const double D[15][9][3] = {{{1, 1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, -1, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0.866025403784438, 0.433012701892219, 0}, {0.866025403784439, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784439, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}, {0.866025403784438, 0.433012701892219, 0}}, {{0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}, {0, 0.75, 0}}, {{-0.17542966950853, 0.148329095604429, 0}, {-0.0180753884489578, 0.444221856552778, 0}, {0.0340120331760423, 0.307261606416198, 0}, {-0.396781833524282, 3.42776476771047e-16, 0}, {0.0799728864560174, -2.52114035295885e-16, 0}, {0.171602319815636, -3.84645013221013e-15, 0}, {-0.0271005739041022, -0.14832909560443, 0}, {0.426146468103823, -0.44422185655278, 0}, {0.341273639592246, -0.307261606416206, 0}}, {{0, 0, 1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, {{0, 0, 1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    
    static const unsigned int mappings[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
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
    static const unsigned int ns[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 9, 9, 1, 1, 1};
    for (unsigned int j = 0; j < ns[i]; j++) {
      // Evaluate basis functions for affine mapping
      const double w0 = 1.0 - X[i][j][0] - X[i][j][1];
      const double w1 = X[i][j][0];
      const double w2 = X[i][j][1];
      
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
        result += values[k]*D[i][j][k];
      // Multiply by weights
      result *= W[i][j];
    
    } // End for
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void darcyflow_1_finite_element_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_1_finite_element_0::interpolate_vertex_values(double* vertex_values,
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
    vertex_values[0] = (1.0/detJ)*(dof_values[3]*3*J_00 + dof_values[4]*(-3*J_00) + dof_values[5]*J_00 + dof_values[6]*(-3*J_01) + dof_values[7]*3*J_01 + dof_values[8]*J_01);
    vertex_values[3] = (1.0/detJ)*(dof_values[0]*3*J_00 + dof_values[1]*(-3*J_00) + dof_values[2]*J_00 + dof_values[6]*(J_00 - J_01) + dof_values[7]*(-3*J_00 + 3*J_01) + dof_values[8]*(3*J_00 - 3*J_01));
    vertex_values[6] = (1.0/detJ)*(dof_values[0]*J_01 + dof_values[1]*(-3*J_01) + dof_values[2]*3*J_01 + dof_values[3]*(J_00 - J_01) + dof_values[4]*(-3*J_00 + 3*J_01) + dof_values[5]*(3*J_00 - 3*J_01));
    vertex_values[1] = (1.0/detJ)*(dof_values[3]*3*J_10 + dof_values[4]*(-3*J_10) + dof_values[5]*J_10 + dof_values[6]*(-3*J_11) + dof_values[7]*3*J_11 + dof_values[8]*J_11);
    vertex_values[4] = (1.0/detJ)*(dof_values[0]*3*J_10 + dof_values[1]*(-3*J_10) + dof_values[2]*J_10 + dof_values[6]*(J_10 - J_11) + dof_values[7]*(-3*J_10 + 3*J_11) + dof_values[8]*(3*J_10 - 3*J_11));
    vertex_values[7] = (1.0/detJ)*(dof_values[0]*J_11 + dof_values[1]*(-3*J_11) + dof_values[2]*3*J_11 + dof_values[3]*(J_10 - J_11) + dof_values[4]*(-3*J_10 + 3*J_11) + dof_values[5]*(3*J_10 - 3*J_11));
    // Evaluate at vertices and use affine mapping
    vertex_values[2] = dof_values[12];
    vertex_values[5] = dof_values[13];
    vertex_values[8] = dof_values[14];
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_1_finite_element_0::num_sub_elements() const
{
    return 2;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_1_finite_element_0::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_1_finite_element_0_0();
      break;
    case 1:
      return new darcyflow_1_finite_element_0_1();
      break;
    }
    return 0;
}


/// Constructor
darcyflow_1_finite_element_1::darcyflow_1_finite_element_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_1_finite_element_1::~darcyflow_1_finite_element_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_1_finite_element_1::signature() const
{
    return "FiniteElement('Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return the cell shape
ufc::shape darcyflow_1_finite_element_1::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_1_finite_element_1::space_dimension() const
{
    return 3;
}

/// Return the rank of the value space
unsigned int darcyflow_1_finite_element_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_1_finite_element_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void darcyflow_1_finite_element_1::evaluate_basis(unsigned int i,
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
void darcyflow_1_finite_element_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_1_finite_element_1::evaluate_basis_derivatives(unsigned int i,
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
void darcyflow_1_finite_element_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_1_finite_element_1::evaluate_dof(unsigned int i,
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
void darcyflow_1_finite_element_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_1_finite_element_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_1_finite_element_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_1_finite_element_1::create_sub_element(unsigned int i) const
{
    return new darcyflow_1_finite_element_1();
}


/// Constructor
darcyflow_1_finite_element_2::darcyflow_1_finite_element_2() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
darcyflow_1_finite_element_2::~darcyflow_1_finite_element_2()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* darcyflow_1_finite_element_2::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return the cell shape
ufc::shape darcyflow_1_finite_element_2::cell_shape() const
{
    return ufc::triangle;
}

/// Return the dimension of the finite element function space
unsigned int darcyflow_1_finite_element_2::space_dimension() const
{
    return 3;
}

/// Return the rank of the value space
unsigned int darcyflow_1_finite_element_2::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int darcyflow_1_finite_element_2::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void darcyflow_1_finite_element_2::evaluate_basis(unsigned int i,
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
void darcyflow_1_finite_element_2::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void darcyflow_1_finite_element_2::evaluate_basis_derivatives(unsigned int i,
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
void darcyflow_1_finite_element_2::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double darcyflow_1_finite_element_2::evaluate_dof(unsigned int i,
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
void darcyflow_1_finite_element_2::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void darcyflow_1_finite_element_2::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
}

/// Return the number of sub elements (for a mixed element)
unsigned int darcyflow_1_finite_element_2::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* darcyflow_1_finite_element_2::create_sub_element(unsigned int i) const
{
    return new darcyflow_1_finite_element_2();
}

/// Constructor
darcyflow_1_dof_map_0_0::darcyflow_1_dof_map_0_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_1_dof_map_0_0::~darcyflow_1_dof_map_0_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_1_dof_map_0_0::signature() const
{
    return "FFC dof map for FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_1_dof_map_0_0::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_1_dof_map_0_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[1] + 3*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_1_dof_map_0_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_1_dof_map_0_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_1_dof_map_0_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_1_dof_map_0_0::local_dimension(const ufc::cell& c) const
{
    return 12;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_1_dof_map_0_0::max_local_dimension() const
{
    return 12;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_1_dof_map_0_0::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_1_dof_map_0_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_1_dof_map_0_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_1_dof_map_0_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[1][0];
    dofs[1] = 3*c.entity_indices[1][0] + 1;
    dofs[2] = 3*c.entity_indices[1][0] + 2;
    dofs[3] = 3*c.entity_indices[1][1];
    dofs[4] = 3*c.entity_indices[1][1] + 1;
    dofs[5] = 3*c.entity_indices[1][1] + 2;
    dofs[6] = 3*c.entity_indices[1][2];
    dofs[7] = 3*c.entity_indices[1][2] + 1;
    dofs[8] = 3*c.entity_indices[1][2] + 2;
    unsigned int offset = 3*m.num_entities[1];
    dofs[9] = offset + 3*c.entity_indices[2][0];
    dofs[10] = offset + 3*c.entity_indices[2][0] + 1;
    dofs[11] = offset + 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_1_dof_map_0_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      break;
    case 1:
      dofs[0] = 3;
      dofs[1] = 4;
      dofs[2] = 5;
      break;
    case 2:
      dofs[0] = 6;
      dofs[1] = 7;
      dofs[2] = 8;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void darcyflow_1_dof_map_0_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_1_dof_map_0_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.75*x[1][0] + 0.25*x[2][0];
    coordinates[0][1] = 0.75*x[1][1] + 0.25*x[2][1];
    coordinates[1][0] = 0.5*x[1][0] + 0.5*x[2][0];
    coordinates[1][1] = 0.5*x[1][1] + 0.5*x[2][1];
    coordinates[2][0] = 0.25*x[1][0] + 0.75*x[2][0];
    coordinates[2][1] = 0.25*x[1][1] + 0.75*x[2][1];
    coordinates[3][0] = 0.75*x[0][0] + 0.25*x[2][0];
    coordinates[3][1] = 0.75*x[0][1] + 0.25*x[2][1];
    coordinates[4][0] = 0.5*x[0][0] + 0.5*x[2][0];
    coordinates[4][1] = 0.5*x[0][1] + 0.5*x[2][1];
    coordinates[5][0] = 0.25*x[0][0] + 0.75*x[2][0];
    coordinates[5][1] = 0.25*x[0][1] + 0.75*x[2][1];
    coordinates[6][0] = 0.75*x[0][0] + 0.25*x[1][0];
    coordinates[6][1] = 0.75*x[0][1] + 0.25*x[1][1];
    coordinates[7][0] = 0.5*x[0][0] + 0.5*x[1][0];
    coordinates[7][1] = 0.5*x[0][1] + 0.5*x[1][1];
    coordinates[8][0] = 0.25*x[0][0] + 0.75*x[1][0];
    coordinates[8][1] = 0.25*x[0][1] + 0.75*x[1][1];
    coordinates[9][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[9][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[10][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[10][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[11][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[11][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int darcyflow_1_dof_map_0_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_1_dof_map_0_0::create_sub_dof_map(unsigned int i) const
{
    return new darcyflow_1_dof_map_0_0();
}


/// Constructor
darcyflow_1_dof_map_0_1::darcyflow_1_dof_map_0_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_1_dof_map_0_1::~darcyflow_1_dof_map_0_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_1_dof_map_0_1::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_1_dof_map_0_1::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_1_dof_map_0_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_1_dof_map_0_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_1_dof_map_0_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_1_dof_map_0_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_1_dof_map_0_1::local_dimension(const ufc::cell& c) const
{
    return 3;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_1_dof_map_0_1::max_local_dimension() const
{
    return 3;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_1_dof_map_0_1::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_1_dof_map_0_1::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_1_dof_map_0_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_1_dof_map_0_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[2][0];
    dofs[1] = 3*c.entity_indices[2][0] + 1;
    dofs[2] = 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_1_dof_map_0_1::tabulate_facet_dofs(unsigned int* dofs,
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
void darcyflow_1_dof_map_0_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_1_dof_map_0_1::tabulate_coordinates(double** coordinates,
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
unsigned int darcyflow_1_dof_map_0_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_1_dof_map_0_1::create_sub_dof_map(unsigned int i) const
{
    return new darcyflow_1_dof_map_0_1();
}


/// Constructor
darcyflow_1_dof_map_0::darcyflow_1_dof_map_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_1_dof_map_0::~darcyflow_1_dof_map_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_1_dof_map_0::signature() const
{
    return "FFC dof map for MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) })";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_1_dof_map_0::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_1_dof_map_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[1] + 6*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_1_dof_map_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_1_dof_map_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_1_dof_map_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_1_dof_map_0::local_dimension(const ufc::cell& c) const
{
    return 15;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_1_dof_map_0::max_local_dimension() const
{
    return 15;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_1_dof_map_0::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_1_dof_map_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_1_dof_map_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_1_dof_map_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[1][0];
    dofs[1] = 3*c.entity_indices[1][0] + 1;
    dofs[2] = 3*c.entity_indices[1][0] + 2;
    dofs[3] = 3*c.entity_indices[1][1];
    dofs[4] = 3*c.entity_indices[1][1] + 1;
    dofs[5] = 3*c.entity_indices[1][1] + 2;
    dofs[6] = 3*c.entity_indices[1][2];
    dofs[7] = 3*c.entity_indices[1][2] + 1;
    dofs[8] = 3*c.entity_indices[1][2] + 2;
    unsigned int offset = 3*m.num_entities[1];
    dofs[9] = offset + 3*c.entity_indices[2][0];
    dofs[10] = offset + 3*c.entity_indices[2][0] + 1;
    dofs[11] = offset + 3*c.entity_indices[2][0] + 2;
    offset = offset + 3*m.num_entities[2];
    dofs[12] = offset + 3*c.entity_indices[2][0];
    dofs[13] = offset + 3*c.entity_indices[2][0] + 1;
    dofs[14] = offset + 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_1_dof_map_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      break;
    case 1:
      dofs[0] = 3;
      dofs[1] = 4;
      dofs[2] = 5;
      break;
    case 2:
      dofs[0] = 6;
      dofs[1] = 7;
      dofs[2] = 8;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void darcyflow_1_dof_map_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_1_dof_map_0::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.75*x[1][0] + 0.25*x[2][0];
    coordinates[0][1] = 0.75*x[1][1] + 0.25*x[2][1];
    coordinates[1][0] = 0.5*x[1][0] + 0.5*x[2][0];
    coordinates[1][1] = 0.5*x[1][1] + 0.5*x[2][1];
    coordinates[2][0] = 0.25*x[1][0] + 0.75*x[2][0];
    coordinates[2][1] = 0.25*x[1][1] + 0.75*x[2][1];
    coordinates[3][0] = 0.75*x[0][0] + 0.25*x[2][0];
    coordinates[3][1] = 0.75*x[0][1] + 0.25*x[2][1];
    coordinates[4][0] = 0.5*x[0][0] + 0.5*x[2][0];
    coordinates[4][1] = 0.5*x[0][1] + 0.5*x[2][1];
    coordinates[5][0] = 0.25*x[0][0] + 0.75*x[2][0];
    coordinates[5][1] = 0.25*x[0][1] + 0.75*x[2][1];
    coordinates[6][0] = 0.75*x[0][0] + 0.25*x[1][0];
    coordinates[6][1] = 0.75*x[0][1] + 0.25*x[1][1];
    coordinates[7][0] = 0.5*x[0][0] + 0.5*x[1][0];
    coordinates[7][1] = 0.5*x[0][1] + 0.5*x[1][1];
    coordinates[8][0] = 0.25*x[0][0] + 0.75*x[1][0];
    coordinates[8][1] = 0.25*x[0][1] + 0.75*x[1][1];
    coordinates[9][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[9][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[10][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[10][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[11][0] = 0.80869438567767*x[0][0] + 0.102717654809626*x[1][0] + 0.088587959512704*x[2][0];
    coordinates[11][1] = 0.80869438567767*x[0][1] + 0.102717654809626*x[1][1] + 0.088587959512704*x[2][1];
    coordinates[12][0] = x[0][0];
    coordinates[12][1] = x[0][1];
    coordinates[13][0] = x[1][0];
    coordinates[13][1] = x[1][1];
    coordinates[14][0] = x[2][0];
    coordinates[14][1] = x[2][1];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int darcyflow_1_dof_map_0::num_sub_dof_maps() const
{
    return 2;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_1_dof_map_0::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_1_dof_map_0_0();
      break;
    case 1:
      return new darcyflow_1_dof_map_0_1();
      break;
    }
    return 0;
}


/// Constructor
darcyflow_1_dof_map_1::darcyflow_1_dof_map_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_1_dof_map_1::~darcyflow_1_dof_map_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_1_dof_map_1::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_1_dof_map_1::needs_mesh_entities(unsigned int d) const
{
    switch ( d )
    {
    case 0:
      return true;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    }
    return false;
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool darcyflow_1_dof_map_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_1_dof_map_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_1_dof_map_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_1_dof_map_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_1_dof_map_1::local_dimension(const ufc::cell& c) const
{
    return 3;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_1_dof_map_1::max_local_dimension() const
{
    return 3;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_1_dof_map_1::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_1_dof_map_1::num_facet_dofs() const
{
    return 2;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_1_dof_map_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_1_dof_map_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_1_dof_map_1::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void darcyflow_1_dof_map_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_1_dof_map_1::tabulate_coordinates(double** coordinates,
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
unsigned int darcyflow_1_dof_map_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_1_dof_map_1::create_sub_dof_map(unsigned int i) const
{
    return new darcyflow_1_dof_map_1();
}


/// Constructor
darcyflow_1_dof_map_2::darcyflow_1_dof_map_2() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
darcyflow_1_dof_map_2::~darcyflow_1_dof_map_2()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* darcyflow_1_dof_map_2::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool darcyflow_1_dof_map_2::needs_mesh_entities(unsigned int d) const
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
bool darcyflow_1_dof_map_2::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 3*m.num_entities[2];
    return false;
}

/// Initialize dof map for given cell
void darcyflow_1_dof_map_2::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void darcyflow_1_dof_map_2::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int darcyflow_1_dof_map_2::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int darcyflow_1_dof_map_2::local_dimension(const ufc::cell& c) const
{
    return 3;
}

/// Return the maximum dimension of the local finite element function space
unsigned int darcyflow_1_dof_map_2::max_local_dimension() const
{
    return 3;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int darcyflow_1_dof_map_2::geometric_dimension() const
{
    return 2;
}

/// Return the number of dofs on each cell facet
unsigned int darcyflow_1_dof_map_2::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int darcyflow_1_dof_map_2::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void darcyflow_1_dof_map_2::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = 3*c.entity_indices[2][0];
    dofs[1] = 3*c.entity_indices[2][0] + 1;
    dofs[2] = 3*c.entity_indices[2][0] + 2;
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void darcyflow_1_dof_map_2::tabulate_facet_dofs(unsigned int* dofs,
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
void darcyflow_1_dof_map_2::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void darcyflow_1_dof_map_2::tabulate_coordinates(double** coordinates,
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
unsigned int darcyflow_1_dof_map_2::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* darcyflow_1_dof_map_2::create_sub_dof_map(unsigned int i) const
{
    return new darcyflow_1_dof_map_2();
}


/// Constructor
darcyflow_1_cell_integral_0_quadrature::darcyflow_1_cell_integral_0_quadrature() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_1_cell_integral_0_quadrature::~darcyflow_1_cell_integral_0_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void darcyflow_1_cell_integral_0_quadrature::tabulate_tensor(double* A,
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
    static const double W4[4] = {0.159020690871988, 0.0909793091280113, 0.159020690871988, 0.0909793091280113};
    // Quadrature points on the UFC reference element: (0.178558728263616, 0.155051025721682), (0.0750311102226081, 0.644948974278318), (0.666390246014701, 0.155051025721682), (0.280019915499074, 0.644948974278318)
    
    // Value of basis functions at quadrature points.
    static const double FE0[4][3] = \
    {{0.666390246014701, 0.178558728263616, 0.155051025721682},
    {0.280019915499074, 0.0750311102226082, 0.644948974278318},
    {0.178558728263616, 0.666390246014701, 0.155051025721682},
    {0.0750311102226081, 0.280019915499074, 0.644948974278318}};
    
    static const double FE1_C2[4][15] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.666390246014701, 0.178558728263616, 0.155051025721682},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.280019915499074, 0.0750311102226082, 0.644948974278318},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.178558728263616, 0.666390246014701, 0.155051025721682},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0750311102226081, 0.280019915499074, 0.644948974278318}};
    
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', False), ('ignore zero tables', False), ('non zero columns', False), ('remove zero terms', False), ('ignore ones', False)
    // Total number of operations to compute element tensor: 264
    
    // Loop quadrature points for integral
    // Number of operations to compute element tensor for following IP loop = 264
    for (unsigned int ip = 0; ip < 4; ip++)
    {
      
      // Function declarations
      double F0 = 0;
      
      // Total number of operations to compute function values = 6
      for (unsigned int r = 0; r < 3; r++)
      {
        F0 += FE0[ip][r]*w[1][r];
      }// end loop over 'r'
      
      // Number of operations for primary indices: 60
      for (unsigned int j = 0; j < 15; j++)
      {
        // Number of operations to compute entry: 4
        A[j] += FE1_C2[ip][j]*F0*W4[ip]*det;
      }// end loop over 'j'
    }// end loop over 'ip'
}

/// Constructor
darcyflow_1_cell_integral_0::darcyflow_1_cell_integral_0() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_1_cell_integral_0::~darcyflow_1_cell_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void darcyflow_1_cell_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c) const
{
    // Reset values of the element tensor block
    for (unsigned int j = 0; j < 15; j++)
      A[j] = 0;
    
    // Add all contributions to element tensor
    integral_0_quadrature.tabulate_tensor(A, w, c);
}

/// Constructor
darcyflow_1_exterior_facet_integral_0::darcyflow_1_exterior_facet_integral_0() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_1_exterior_facet_integral_0::~darcyflow_1_exterior_facet_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void darcyflow_1_exterior_facet_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c,
                                    unsigned int facet) const
{
    // Reset values of the element tensor block
    for (unsigned int j = 0; j < 15; j++)
      A[j] = 0;
    
    // Add all contributions to element tensor
}

/// Constructor
darcyflow_1_exterior_facet_integral_1_quadrature::darcyflow_1_exterior_facet_integral_1_quadrature() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_1_exterior_facet_integral_1_quadrature::~darcyflow_1_exterior_facet_integral_1_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void darcyflow_1_exterior_facet_integral_1_quadrature::tabulate_tensor(double* A,
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
    static const double FE0_f0[2][3] = \
    {{0, 0.788675134594813, 0.211324865405187},
    {0, 0.211324865405187, 0.788675134594813}};
    
    static const double FE0_f1[2][3] = \
    {{0.788675134594813, 0, 0.211324865405187},
    {0.211324865405187, 0, 0.788675134594813}};
    
    static const double FE0_f2[2][3] = \
    {{0.788675134594813, 0.211324865405187, 0},
    {0.211324865405187, 0.788675134594813, 0}};
    
    static const double FE1_f0_C0[2][15] = \
    {{1.36602540378444, -1.03269207045111, 0.455341801261482, -0.0386751345948141, 0.199358737117775, -0.449358737117775, 0.538675134594813, -1.53269207045111, 1.28269207045111, 0.288675134594813, -0.499999999999998, 1.58113883008419, 0, 0, 0},
    {-0.366025403784437, 0.699358737117768, -0.122008467928144, 0.538675134594811, -1.5326920704511, 1.2826920704511, -0.0386751345948124, 0.199358737117772, -0.44935873711777, 0.288675134594814, -0.499999999999998, 1.58113883008419, 0, 0, 0}};
    
    static const double FE1_f0_C1[2][15] = \
    {{-0.122008467928146, 0.699358737117776, -0.366025403784439, 0.0386751345948142, -0.199358737117775, 0.449358737117775, -0.538675134594813, 1.53269207045111, -1.28269207045111, -0.288675134594814, 0.499999999999998, -1.58113883008419, 0, 0, 0},
    {0.455341801261479, -1.0326920704511, 1.36602540378444, -0.538675134594811, 1.5326920704511, -1.2826920704511, 0.0386751345948124, -0.199358737117772, 0.44935873711777, -0.288675134594814, 0.499999999999998, -1.58113883008419, 0, 0, 0}};
    
    static const double FE1_f1_C0[2][15] = \
    {{0, 0, 0, 1.24401693585629, -0.333333333333332, 0.0893163974770408, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0.0893163974770407, -0.333333333333335, 1.24401693585629, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    static const double FE1_f1_C1[2][15] = \
    {{-0.122008467928146, 0.0326920704511013, -0.366025403784439, 0.205341801261478, -0.532692070451102, 0.282692070451102, -1.2826920704511, 1.86602540378444, -0.538675134594811, 0.288675134594814, 0.833333333333336, 1.58113883008419, 0, 0, 0},
    {0.455341801261479, -1.69935873711778, 1.36602540378444, -0.372008467928147, 1.19935873711777, -1.44935873711778, 0.449358737117772, 0.133974596215561, 0.0386751345948144, 0.288675134594814, 0.833333333333335, 1.58113883008419, 0, 0, 0}};
    
    static const double FE1_f2_C0[2][15] = \
    {{-0.36602540378444, 0.032692070451109, -0.122008467928148, 1.28269207045111, -1.86602540378444, 0.538675134594816, -0.20534180126148, 0.532692070451106, -0.282692070451107, 0.866025403784439, -0.166666666666669, -1.58113883008419, 0, 0, 0},
    {1.36602540378444, -1.69935873711777, 0.455341801261477, -0.449358737117771, -0.133974596215562, -0.0386751345948098, 0.372008467928146, -1.19935873711777, 1.44935873711777, 0.866025403784438, -0.166666666666669, -1.58113883008419, 0, 0, 0}};
    
    static const double FE1_f2_C1[2][15] = \
    {{0, 0, 0, 0, 0, 0, -1.24401693585629, 0.333333333333333, -0.0893163974770408, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, -0.0893163974770408, 0.333333333333333, -1.24401693585629, 0, 0, 0, 0, 0, 0}};
    
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', False), ('ignore zero tables', False), ('non zero columns', False), ('remove zero terms', False), ('ignore ones', False)
    switch ( facet )
    {
    case 0:
      {
      // Total number of operations to compute element tensor (from this point): 552
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 552
      for (unsigned int ip = 0; ip < 2; ip++)
      {
        
        // Function declarations
        double F0 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 3; r++)
        {
          F0 += FE0_f0[ip][r]*w[0][r];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 270
        for (unsigned int j = 0; j < 15; j++)
        {
          // Number of operations to compute entry: 18
          A[j] += (((1.0/detJ)*J_00*FE1_f0_C0[ip][j] + (1.0/detJ)*J_01*FE1_f0_C1[ip][j])*n0 + ((1.0/detJ)*J_10*FE1_f0_C0[ip][j] + (1.0/detJ)*J_11*FE1_f0_C1[ip][j])*n1)*F0*-1*W2[ip]*det;
        }// end loop over 'j'
      }// end loop over 'ip'
      }
      break;
    case 1:
      {
      // Total number of operations to compute element tensor (from this point): 552
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 552
      for (unsigned int ip = 0; ip < 2; ip++)
      {
        
        // Function declarations
        double F0 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 3; r++)
        {
          F0 += FE0_f1[ip][r]*w[0][r];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 270
        for (unsigned int j = 0; j < 15; j++)
        {
          // Number of operations to compute entry: 18
          A[j] += (((1.0/detJ)*J_00*FE1_f1_C0[ip][j] + (1.0/detJ)*J_01*FE1_f1_C1[ip][j])*n0 + ((1.0/detJ)*J_10*FE1_f1_C0[ip][j] + (1.0/detJ)*J_11*FE1_f1_C1[ip][j])*n1)*F0*-1*W2[ip]*det;
        }// end loop over 'j'
      }// end loop over 'ip'
      }
      break;
    case 2:
      {
      // Total number of operations to compute element tensor (from this point): 552
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 552
      for (unsigned int ip = 0; ip < 2; ip++)
      {
        
        // Function declarations
        double F0 = 0;
        
        // Total number of operations to compute function values = 6
        for (unsigned int r = 0; r < 3; r++)
        {
          F0 += FE0_f2[ip][r]*w[0][r];
        }// end loop over 'r'
        
        // Number of operations for primary indices: 270
        for (unsigned int j = 0; j < 15; j++)
        {
          // Number of operations to compute entry: 18
          A[j] += (((1.0/detJ)*J_00*FE1_f2_C0[ip][j] + (1.0/detJ)*J_01*FE1_f2_C1[ip][j])*n0 + ((1.0/detJ)*J_10*FE1_f2_C0[ip][j] + (1.0/detJ)*J_11*FE1_f2_C1[ip][j])*n1)*F0*-1*W2[ip]*det;
        }// end loop over 'j'
      }// end loop over 'ip'
      }
      break;
    }
}

/// Constructor
darcyflow_1_exterior_facet_integral_1::darcyflow_1_exterior_facet_integral_1() : ufc::exterior_facet_integral()
{
    // Do nothing
}

/// Destructor
darcyflow_1_exterior_facet_integral_1::~darcyflow_1_exterior_facet_integral_1()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local exterior facet
void darcyflow_1_exterior_facet_integral_1::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c,
                                    unsigned int facet) const
{
    // Reset values of the element tensor block
    for (unsigned int j = 0; j < 15; j++)
      A[j] = 0;
    
    // Add all contributions to element tensor
    integral_1_quadrature.tabulate_tensor(A, w, c, facet);
}

/// Constructor
darcyflow_form_1::darcyflow_form_1() : ufc::form()
{
    // Do nothing
}

/// Destructor
darcyflow_form_1::~darcyflow_form_1()
{
    // Do nothing
}

/// Return a string identifying the form
const char* darcyflow_form_1::signature() const
{
    return "Form([Integral(Product(Function(FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1), 1), Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(2),), {FixedIndex(2): 3}))), Measure('cell', 0, None)), Integral(Product(IntValue(-1, (), (), {}), Product(Function(FiniteElement('Lagrange', Cell('triangle', 1, Space(2)), 1), 0), IndexSum(Product(Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((Index(0),), {Index(0): 2})), Indexed(ListTensor(Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 3})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Brezzi-Douglas-Marini', Cell('triangle', 1, Space(2)), 2), FiniteElement('Discontinuous Lagrange', Cell('triangle', 1, Space(2)), 1)], **{'value_shape': (3,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 3}))), MultiIndex((Index(0),), {Index(0): 2}))), MultiIndex((Index(0),), {Index(0): 2})))), Measure('exterior_facet', 1, None))])";
}

/// Return the rank of the global tensor (r)
unsigned int darcyflow_form_1::rank() const
{
    return 1;
}

/// Return the number of coefficients (n)
unsigned int darcyflow_form_1::num_coefficients() const
{
    return 2;
}

/// Return the number of cell integrals
unsigned int darcyflow_form_1::num_cell_integrals() const
{
    return 1;
}

/// Return the number of exterior facet integrals
unsigned int darcyflow_form_1::num_exterior_facet_integrals() const
{
    return 2;
}

/// Return the number of interior facet integrals
unsigned int darcyflow_form_1::num_interior_facet_integrals() const
{
    return 0;
}

/// Create a new finite element for argument function i
ufc::finite_element* darcyflow_form_1::create_finite_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_1_finite_element_0();
      break;
    case 1:
      return new darcyflow_1_finite_element_1();
      break;
    case 2:
      return new darcyflow_1_finite_element_2();
      break;
    }
    return 0;
}

/// Create a new dof map for argument function i
ufc::dof_map* darcyflow_form_1::create_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_1_dof_map_0();
      break;
    case 1:
      return new darcyflow_1_dof_map_1();
      break;
    case 2:
      return new darcyflow_1_dof_map_2();
      break;
    }
    return 0;
}

/// Create a new cell integral on sub domain i
ufc::cell_integral* darcyflow_form_1::create_cell_integral(unsigned int i) const
{
    return new darcyflow_1_cell_integral_0();
}

/// Create a new exterior facet integral on sub domain i
ufc::exterior_facet_integral* darcyflow_form_1::create_exterior_facet_integral(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new darcyflow_1_exterior_facet_integral_0();
      break;
    case 1:
      return new darcyflow_1_exterior_facet_integral_1();
      break;
    }
    return 0;
}

/// Create a new interior facet integral on sub domain i
ufc::interior_facet_integral* darcyflow_form_1::create_interior_facet_integral(unsigned int i) const
{
    return 0;
}


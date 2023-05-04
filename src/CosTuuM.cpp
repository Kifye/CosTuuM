/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CosTuuM is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * CosTuuM is distributed in the hope that it will be useful, but WITOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CosTuuM. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file CosTuuM.cpp
 *
 * @brief Main program entry point.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "CommandLineParser.hpp"
#include "Configuration.hpp"
#include "ConfigurationInfo.hpp"
#include "Matrix.hpp"
#include "ParameterFile.hpp"
#include "TMatrixCalculator.hpp"
#include "Utilities.hpp"

#include <cinttypes>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

template<typename StartType, typename EndType>
std::vector<EndType> parse_to_typed_array(ParameterFile& params, const std::string& key, const std::vector<StartType>& default_value) {
    const std::vector<StartType> double_vector =
      params.get_value<std::vector<StartType>>(key, default_value);

  std::vector<EndType> float_vector(double_vector.size());
  for (size_t i = 0; i < double_vector.size(); ++i) {
    float_vector[i] = double_vector[i];
  }

  return float_vector;
}

template <Quantity _quantity_>
std::vector<float_type> parse_phys_to_typed_array(ParameterFile& params, const std::string& key, const std::string& default_str) {
    std::vector<std::string> parts;
    std::string value = params.get_value<std::string>(key);
    Utilities::split_string(value, parts);
    std::vector<double> double_vector;
    for (auto& part : parts) {
        std::pair<double, std::string> valunit = Utilities::split_value(part);
        double_vector.push_back(
            UnitConverter::to_SI<_quantity_>(valunit.first, valunit.second));
    }

  std::vector<float_type> float_vector(double_vector.size());
  for (size_t i = 0; i < double_vector.size(); ++i) {
    float_vector[i] = double_vector[i];
  }

  return float_vector;
}
/**
 * @brief Main program entry point.
 *
 * For now, does the entire T-matrix calculation for the Mishchenko default
 * parameters and outputs the scattering and extinction factors.
 *
 * @param argc Number of command line arguments (currently ignored).
 * @param argv Command line arguments (currently ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  std::cout << "Program configuration:" << std::endl;
  for (auto it = ConfigurationInfo::begin(); it != ConfigurationInfo::end();
       ++it) {
    std::cout << it.get_key() << ": " << it.get_value() << "\n";
  }
  std::cout << std::endl;

  CommandLineParser parser("C++ TMatrix code");
  parser.add_option("params", 'p', "Name of the parameter file.",
                    COMMANDLINEOPTION_STRINGARGUMENT, "", true);
  parser.add_option("output", 'o', "Name of the output dump file.",
                    COMMANDLINEOPTION_STRINGARGUMENT, "", true);
  parser.add_option("tmatrix_type", 't', "The way of obtaining T-matrix: base or spheroidal",
                      COMMANDLINEOPTION_STRINGARGUMENT, "base", false);
  parser.parse_arguments(argc, argv);

  std::cout << "Command line options:" << std::endl;
  parser.print_contents(std::cout);
  std::cout << std::endl;

  ParameterFile params(parser.get_value<std::string>("params"));

  bool use_spheroidal = false;
  if (parser.get_value<std::string>("tmatrix_type")== "spheroidal") {
    use_spheroidal = true;
  }

  /// input parameters

  /// dust particle
  // size of the particle (in same units as the wavelength)
  const std::vector<float_type> axi = parse_phys_to_typed_array<QUANTITY_LENGTH>(params, "DustParticle:size", "10. micron");

    std::cerr << "read axi = ";
    for (const auto& item : axi) {
        std::cerr << item << " ";
    }
    std::cerr << "\n";
  // ratio between the equal surface area sphere radius and equal volume sphere
  // radius (is recomputed if not equal to 1)
  float_type ratio_of_radii = 1.;
  if (!params.get_value<bool>("DustParticle:is equal volume sphere radius",
                              true)) {
    ratio_of_radii = 0.1;
  }
  // refractive index
  const std::vector<std::complex<float_type>> mr = 
  parse_to_typed_array<std::complex<double>, std::complex<float_type>>(
    params, "DustParticle:refractive index", std::vector<std::complex<double>>(1, std::complex<double>(1.5, 0.02)));
    std::cerr << "read mr = ";
    for (const auto& item : mr) {
        std::cerr << item << " ";
    }
    std::cerr << "\n";
  // ratio of horizontal and vertical axis of the spheroidal particle
  std::vector<float_type> axis_ratio = parse_to_typed_array<double, float_type>(params, "DustParticle:axis ratio", std::vector<double>{1, 0.5});
  std::cerr << "read axis_ratio = ";
    for (const auto& item : axis_ratio) {
        std::cerr << item << " ";
    }
    std::cerr << "\n";
  /// radiation
  // wavelength of incoming radiation (in same units as the particle size)
  const float_type wavelength = params.get_physical_value<QUANTITY_LENGTH>(
      "Radiation:incoming wavelength", "6.283185307 micron");

  /// calculation
  // tolerance for the calculation
  const float_type tolerance =
      params.get_value<double>("Calculation:tolerance", 1.e-4);
  // number of Gauss-Legendre points to use as a multiplier of the maximum
  // order of spherical harmonics used to decompose the electric field, during
  // the first loop of the algorithm
  const uint_fast32_t ndgs =
      params.get_value<uint_fast32_t>("Calculation:Gauss Legendre factor", 2);
  // maximum number of iterations during the first loop
  const uint_fast32_t maximum_order =
      params.get_value<uint_fast32_t>("Calculation:maximum order", 200);
  // maximum number of iterations during the second loop
  const uint_fast32_t maximum_ngauss = params.get_value<uint_fast32_t>(
      "Calculation:maximum number of Gauss Legendre points", 500);

  /// scattering event
  // first and second Euler angles describing the orientation of the particle
  // in the fixed laboratory reference frame
  const float_type alpha = params.get_physical_value<QUANTITY_ANGLE>(
      "DustParticle:alpha angle", "145. degrees");
  const float_type beta = params.get_physical_value<QUANTITY_ANGLE>(
      "DustParticle:beta angle", "52. degrees");

  // zenith and azimuth angle of the incoming photon
  const float_type theta_in = params.get_physical_value<QUANTITY_ANGLE>(
      "Radiation:incoming theta", "56. degrees");
  const float_type phi_in = params.get_physical_value<QUANTITY_ANGLE>(
      "Radiation:incoming phi", "114. degrees");

  // zenith and azimuth angle of the scattered photon
  const float_type theta_out = params.get_physical_value<QUANTITY_ANGLE>(
      "Radiation:outgoing theta", "65. degrees");
  const float_type phi_out = params.get_physical_value<QUANTITY_ANGLE>(
      "Radiation:outgoing phi", "128. degrees");

  // write used parameters to file
  {
    std::ofstream pfile("parameters-usedvalues.param");
    params.print_contents(pfile);
    pfile.close();

    std::cout << "Input parameters:" << std::endl;
    params.print_contents(std::cout);
    std::cout << std::endl;
  }

    const uint_fast32_t nol = axi.size();
    if (nol != mr.size()) {
        std::cerr << "wrong size of mr: " << mr.size() << " when expected " << nol;
        exit(1);
    }
        if (nol != axis_ratio.size()) {
        std::cerr << "wrong size of axis_ratio: " << axis_ratio.size() << " when expected " << nol;
        exit(1);
    }
    
  TMatrix *active_Tmatrix = TMatrixCalculator::calculate_TMatrix(nol,
      ratio_of_radii, axis_ratio, axi, wavelength, maximum_order, tolerance,
      ndgs, mr, maximum_ngauss, use_spheroidal);

   OrientationDistribution distribution(2 * active_Tmatrix->get_nmax());
   distribution.initialise();
   TMatrix *average_Tmatrix =
   TMatrixCalculator::apply_orientation_distribution(
         *active_Tmatrix, distribution);

  /// compute a scattering event using the T-matrix

  Matrix<float_type> Z = active_Tmatrix->get_scattering_matrix(
      alpha, beta, theta_in, phi_in, theta_out, phi_out);

  ctm_warning("Z[0,:] = %g %g %g %g", double(Z(0, 0)), double(Z(0, 1)),
              double(Z(0, 2)), double(Z(0, 3)));
  ctm_warning("Z[1,:] = %g %g %g %g", double(Z(1, 0)), double(Z(1, 1)),
              double(Z(1, 2)), double(Z(1, 3)));
  ctm_warning("Z[2,:] = %g %g %g %g", double(Z(2, 0)), double(Z(2, 1)),
              double(Z(2, 2)), double(Z(2, 3)));
  ctm_warning("Z[3,:] = %g %g %g %g", double(Z(3, 0)), double(Z(3, 1)),
              double(Z(3, 2)), double(Z(3, 3)));

  Matrix<float_type> K =
      active_Tmatrix->get_extinction_matrix(alpha, beta, theta_in, phi_in);

  ctm_warning("K[0,:] = %g %g %g %g", double(K(0, 0)), double(K(0, 1)),
              double(K(0, 2)), double(K(0, 3)));
  ctm_warning("K[1,:] = %g %g %g %g", double(K(1, 0)), double(K(1, 1)),
              double(K(1, 2)), double(K(1, 3)));
  ctm_warning("K[2,:] = %g %g %g %g", double(K(2, 0)), double(K(2, 1)),
              double(K(2, 2)), double(K(2, 3)));
  ctm_warning("K[3,:] = %g %g %g %g", double(K(3, 0)), double(K(3, 1)),
              double(K(3, 2)), double(K(3, 3)));

  ctm_warning("Extinction coefficient = %g", double(active_Tmatrix->get_extinction_coefficient()));
  ctm_warning("Scattering coefficient = %g", double(active_Tmatrix->get_scattering_coefficient()));

  ctm_warning("Extinction cross-section = %g", double(active_Tmatrix->get_extinction_coefficient()));
  ctm_warning("Scattering cross-section = %g", double(active_Tmatrix->get_scattering_coefficient()));

  K.binary_dump(parser.get_value<std::string>("output"));

  Matrix<float_type> Z_tot = average_Tmatrix->get_scattering_matrix(
      alpha, beta, theta_in, phi_in, theta_out, phi_out);

  ctm_warning("Z_tot[0,:] = %g %g %g %g", double(Z_tot(0, 0)), double(Z_tot(0, 1)),
              double(Z_tot(0, 2)), double(Z_tot(0, 3)));
  ctm_warning("Z_tot[1,:] = %g %g %g %g", double(Z_tot(1, 0)), double(Z_tot(1, 1)),
              double(Z(1, 2)), double(Z(1, 3)));
  ctm_warning("Z_tot[2,:] = %g %g %g %g", double(Z_tot(2, 0)), double(Z_tot(2, 1)),
              double(Z_tot(2, 2)), double(Z_tot(2, 3)));
  ctm_warning("Z_tot[3,:] = %g %g %g %g", double(Z_tot(3, 0)), double(Z_tot(3, 1)),
              double(Z_tot(3, 2)), double(Z_tot(3, 3)));

  Matrix<float_type> K_tot =
      average_Tmatrix->get_extinction_matrix(alpha, beta, theta_in, phi_in);

  ctm_warning("K_tot[0,:] = %g %g %g %g", double(K_tot(0, 0)), double(K_tot(0, 1)),
              double(K_tot(0, 2)), double(K_tot(0, 3)));
  ctm_warning("K_tot[1,:] = %g %g %g %g", double(K_tot(1, 0)), double(K_tot(1, 1)),
              double(K_tot(1, 2)), double(K_tot(1, 3)));
  ctm_warning("K_tot[2,:] = %g %g %g %g", double(K_tot(2, 0)), double(K_tot(2, 1)),
              double(K_tot(2, 2)), double(K_tot(2, 3)));
  ctm_warning("K_tot[3,:] = %g %g %g %g", double(K_tot(3, 0)), double(K_tot(3, 1)),
              double(K_tot(3, 2)), double(K_tot(3, 3)));

  K_tot.binary_dump(parser.get_value<std::string>("output"));

  // clean up
  delete active_Tmatrix;
  delete average_Tmatrix;

  // done!
  return 0;
}

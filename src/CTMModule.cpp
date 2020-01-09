/**
 * @file CTMModule.cpp
 *
 * @brief CosTuuM Python module.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Configuration.hpp"
#include "DavisGreensteinOrientationDistribution.hpp"
#include "MishchenkoOrientationDistribution.hpp"
#include "SizeBasedAlignmentDistribution.hpp"
#include "TMatrixCalculator.hpp"
#include "TaskManager.hpp"

#include <Python.h>
/*! @brief Use the NumPy 1.7 API. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "structmember.h"

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/// T-matrix object

/**
 * @brief C struct wrapper around a T-matrix object.
 */
typedef struct {
  /*! @brief Python object members. */
  PyObject_HEAD;
  /*! @brief T-matrix object. */
  TMatrix *_Tmatrix;
} TmatrixObject;

/**
 * @brief Destructor for the T-matrix object.
 *
 * @param self T-matrix object that is being deallocated.
 */
static void TmatrixObject_dealloc(TmatrixObject *self) {
  // first free the wrapped T-matrix object
  delete self->_Tmatrix;
  // then free the Python object
  Py_TYPE(self)->tp_free(reinterpret_cast<PyObject *>(self));
}

/**
 * @brief Constructor for the T-matrix object.
 *
 * Constructs an empty T-matrix object.
 *
 * @param type Type for the object.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Pointer to the newly created object.
 */
static PyObject *TmatrixObject_new(PyTypeObject *type, PyObject *args,
                                   PyObject *kwargs) {
  // allocate the Python object
  TmatrixObject *self =
      reinterpret_cast<TmatrixObject *>(type->tp_alloc(type, 0));
  // set the wrapped T-matrix to a sensible initial value
  self->_Tmatrix = nullptr;
  return reinterpret_cast<PyObject *>(self);
}

/**
 * @brief init() function for the T-matrix object.
 *
 * Required arguments are:
 *  - particle_radius: Radius of the scattering particle (in m).
 *  - axis_ratio: Axis ratio of the scattering particle.
 *  - wavelength: Wavelength of the scattered radiation (in m).
 *  - refractive_index: Complex refractive index of the scattering particle.
 *
 * Additional optional arguments are:
 *  - tolerance: Tolerance for the iterative T-matrix calculation.
 *  - maximum_order: Maximum order to use in the spherical basis function
 *    expansion.
 *  - gauss_legendre_factor: Multiplicative factor expressing the number of
 *    Gauss-Legendre quadrature points used during the first phase of the
 *    T-matrix calculation as a multiple of the order for each iteration.
 *  - maximum_ngauss: Maximum number of Gauss-Legendre quadrature points to use
 *    during the T-matrix calculation.
 *  - is_equal_volume_radius: Is the input particle radius an equal volume
 *    radius?
 *
 * @param self T-matrix object that is being initialised.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return 0 on success, 1 on failure.
 */
static int TmatrixObject_init(TmatrixObject *self, PyObject *args,
                              PyObject *kwargs) {

  // required arguments
  float_type particle_radius;
  float_type axis_ratio;
  float_type wavelength;
  std::complex<float_type> mr;

  // optional arguments
  float_type cos2beta = 1. / 3.;
  float_type tolerance = 1.e-4;
  uint_fast32_t maximum_order = 200;
  uint_fast32_t ndgs = 2;
  uint_fast32_t maximum_ngauss = 500;
  float_type ratio_of_radii = 1.;

  /// parse arguments
  // list of keywords (in the expected order)
  // note that we need to use strdup because the Python API expects char*,
  // while C++ strings are const char*
  // not doing this results in compilation warnings
  static char *kwlist[] = {strdup("particle_radius"),
                           strdup("axis_ratio"),
                           strdup("wavelength"),
                           strdup("refractive_index"),
                           strdup("cos2beta"),
                           strdup("tolerance"),
                           strdup("maximum_order"),
                           strdup("gauss_legendre_factor"),
                           strdup("maximum_ngauss"),
                           strdup("is_equal_volume_radius"),
                           nullptr};

  // allocate temporary variables to store double precision arguments
  double particle_radius_d, axis_ratio_d, wavelength_d;
  double cos2beta_d = static_cast<double>(cos2beta);
  double tolerance_d = static_cast<double>(tolerance);
  // allocate a temporary variable to store the complex refractive index
  Py_complex mr_temp;
  // allocate a temporary variable to store the equal volume radius boolean
  int_fast32_t is_equal_volume_radius = true;
  // parse the keywords/positional arguments
  // d is a double
  // D is a complex double (Py_complex)
  // I is an unsigned integer
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "dddD|ddIIIp", kwlist, &particle_radius_d,
          &axis_ratio_d, &wavelength_d, &mr_temp, &cos2beta_d, &tolerance_d,
          &maximum_order, &ndgs, &maximum_ngauss, &is_equal_volume_radius)) {
    // PyArg_ParseTupleAndKeywords will return 0 if a required argument was
    // missing, if an argument of the wrong type was provided or if the number
    // of arguments does not match the expectation
    // we use ctm_warning so that we can gracefully exit and let Python handle
    // the exception
    // ctm_error would call abort, which would kill the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // exit code 1 signals to Python that something was wrong
    return 1;
  }
  // unpack double precision arguments
  particle_radius = particle_radius_d;
  axis_ratio = axis_ratio_d;
  wavelength = wavelength_d;
  cos2beta = cos2beta_d;
  tolerance = tolerance_d;
  // get the complex components from the Py_complex and store them in the
  // std::complex variable
  mr.real(mr_temp.real);
  mr.imag(mr_temp.imag);
  // set the ratio_of_radii if it is not correct
  if (!is_equal_volume_radius) {
    ratio_of_radii = 0.1;
  }

  // construct the T-matrix
  self->_Tmatrix = TMatrixCalculator::calculate_TMatrix(
      ratio_of_radii, axis_ratio, particle_radius, wavelength, maximum_order,
      tolerance, ndgs, mr, maximum_ngauss);

  // average it out over the orientation distribution
  OrientationDistribution *orientation;
  if (cos2beta == 0. || cos2beta == 1.) {
    orientation = new DavisGreensteinOrientationDistribution(
        2 * self->_Tmatrix->get_nmax(), axis_ratio);
  } else if (cos2beta == 1. / 3.) {
    orientation = new OrientationDistribution(2 * self->_Tmatrix->get_nmax());
    orientation->initialise();
  } else {
    orientation = new MishchenkoOrientationDistribution(
        2 * self->_Tmatrix->get_nmax(), cos2beta);
  }
  self->_Tmatrix = TMatrixCalculator::apply_orientation_distribution(
      *self->_Tmatrix, *orientation);
  delete orientation;

  // everything went well: return 0
  return 0;
}

/**
 * @brief Get the maximum order for a T-matrix object.
 *
 * @param self T-matrix object.
 * @param Py_UNUSED Additional arguments are not used but need to be present.
 * @return Integer Python object containing the maximum order.
 */
static PyObject *TmatrixObject_get_nmax(TmatrixObject *self,
                                        PyObject *Py_UNUSED(ignored)) {
  // wrap the returned value in a Python object
  return PyLong_FromUnsignedLong(self->_Tmatrix->get_nmax());
}

/**
 * @brief Get the absorption cross section appropriately averaged over the
 * outgoing angle @f$\theta{}@f$.
 *
 * @param self T-matrix object.
 * @param Py_UNUSED Additional arguments are not used but need to be present.
 * @return Integer Python object containing the maximum order.
 */
static PyObject *TmatrixObject_get_average_absorption_cross_section(
    TmatrixObject *self, PyObject *Py_UNUSED(ignored)) {
  // wrap the returned value in a Python object
  return PyFloat_FromDouble(static_cast<double>(
      self->_Tmatrix->get_average_absorption_cross_section(20)));
}

/**
 * @brief Get the absorption cross section for the given input angle(s).
 *
 * Required arguments are:
 *  - theta: Input zenith angle (in radians).
 *
 * Additional optional arguments are:
 *  - ngauss: Number of Gauss-Legendre quadrature points to use to compute the
 *    integral to subtract the scattering contribution to the extinction.
 *
 * @param self T-matrix object being used.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Pointer to an absorption cross section value or array that has the
 * same shape as the input angles.
 */
static PyObject *TmatrixObject_get_absorption_cross_section(TmatrixObject *self,
                                                            PyObject *args,
                                                            PyObject *kwargs) {

  // required arguments
  PyArrayObject *thetas;

  // optional arguments
  uint_fast32_t ngauss = 100;

  // list of keywords (see comment above)
  static char *kwlist[] = {strdup("theta"), strdup("ngauss"), nullptr};

  // parse positional and keyword arguments
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O&|I", kwlist,
                                   PyArray_Converter, &thetas, &ngauss)) {
    // again, we do not call ctm_error to avoid killing the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // this time, a nullptr return will signal an error to Python
    return nullptr;
  }

  // determine the size of the input arrays
  npy_intp thetasize;

  // get the number of dimensions for the theta array
  const npy_intp thetadim = PyArray_NDIM(thetas);
  // we only accept 0D (scalar) or 1D input
  if (thetadim > 1) {
    ctm_warning("Wrong shape for input array!");
    return nullptr;
  }
  // check if we are in the 0D or 1D case
  if (thetadim > 0) {
    // if 1D, simply get the size from the 1 dimension
    const npy_intp *thetadims = PyArray_DIMS(thetas);
    thetasize = thetadims[0];
  } else {
    // if 0D, set the size to 1
    thetasize = 1;
    // reshape the array into a 1D array with 1 element, so that we can
    // manipulate it in the same way as a 1D array
    npy_intp newdims[1] = {1};
    PyArray_Dims newdimsobj;
    newdimsobj.ptr = newdims;
    newdimsobj.len = 1;
    thetas = reinterpret_cast<PyArrayObject *>(
        PyArray_Newshape(thetas, &newdimsobj, NPY_ANYORDER));
  }

  // create an uninitialised double NumPy array to store the results
  // shape: thetasize
  npy_intp dims[1] = {thetasize};
  PyArrayObject *Cabs =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNew(1, dims, NPY_DOUBLE));

  const double *thetadata = reinterpret_cast<double *>(PyArray_DATA(thetas));
  std::vector<float_type> thetadata_f(thetasize);
  for (npy_intp i = 0; i < thetasize; ++i) {
    thetadata_f[i] = thetadata[i];
  }

  double *Cabsdata = reinterpret_cast<double *>(PyArray_DATA(Cabs));
  std::vector<float_type> Cabsdata_f(thetasize);
  self->_Tmatrix->get_absorption_cross_section(&thetadata_f[0], thetasize,
                                               ngauss, &Cabsdata_f[0]);
  for (npy_intp i = 0; i < thetasize; ++i) {
    Cabsdata[i] = static_cast<double>(Cabsdata_f[i]);
  }

  // tell Python we are done with the input objects (so that memory is properly
  // deallocated)
  Py_DECREF(thetas);

  return PyArray_Return(Cabs);
}

/**
 * @brief Get the absorption cross sections for the given input angle(s).
 *
 * Required arguments are:
 *  - theta: Input zenith angle (in radians).
 *
 * Additional optional arguments are:
 *  - ngauss: Number of Gauss-Legendre quadrature points to use to compute the
 *    integral to subtract the scattering contribution to the extinction.
 *
 * @param self T-matrix object being used.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Pointer to an array containing the cross section values. The array
 * has two elements per input angle theta.
 */
static PyObject *
TmatrixObject_get_absorption_cross_sections(TmatrixObject *self, PyObject *args,
                                            PyObject *kwargs) {

  // required arguments
  PyArrayObject *thetas;

  // optional arguments
  uint_fast32_t ngauss = 100;

  // list of keywords (see comment above)
  static char *kwlist[] = {strdup("theta"), strdup("ngauss"), nullptr};

  // parse positional and keyword arguments
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O&|I", kwlist,
                                   PyArray_Converter, &thetas, &ngauss)) {
    // again, we do not call ctm_error to avoid killing the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // this time, a nullptr return will signal an error to Python
    return nullptr;
  }

  // determine the size of the input arrays
  npy_intp thetasize;

  // get the number of dimensions for the theta array
  const npy_intp thetadim = PyArray_NDIM(thetas);
  // we only accept 0D (scalar) or 1D input
  if (thetadim > 1) {
    ctm_warning("Wrong shape for input array!");
    return nullptr;
  }
  // check if we are in the 0D or 1D case
  if (thetadim > 0) {
    // if 1D, simply get the size from the 1 dimension
    const npy_intp *thetadims = PyArray_DIMS(thetas);
    thetasize = thetadims[0];
  } else {
    // if 0D, set the size to 1
    thetasize = 1;
    // reshape the array into a 1D array with 1 element, so that we can
    // manipulate it in the same way as a 1D array
    npy_intp newdims[1] = {1};
    PyArray_Dims newdimsobj;
    newdimsobj.ptr = newdims;
    newdimsobj.len = 1;
    thetas = reinterpret_cast<PyArrayObject *>(
        PyArray_Newshape(thetas, &newdimsobj, NPY_ANYORDER));
  }

  // create an uninitialised double NumPy array to store the results
  // shape: thetasize x 2
  npy_intp dims[2] = {thetasize, 2};
  PyArrayObject *Cabs =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNew(2, dims, NPY_DOUBLE));

  const double *thetadata = reinterpret_cast<double *>(PyArray_DATA(thetas));
  std::vector<float_type> thetadata_f(thetasize);
  for (npy_intp i = 0; i < thetasize; ++i) {
    thetadata_f[i] = thetadata[i];
  }
  double *Cabsdata = reinterpret_cast<double *>(PyArray_DATA(Cabs));
  std::vector<float_type> Cabsdata_f(2 * thetasize);
  self->_Tmatrix->get_absorption_cross_sections(&thetadata_f[0], thetasize,
                                                ngauss, &Cabsdata_f[0]);
  for (npy_intp i = 0; i < 2 * thetasize; ++i) {
    Cabsdata[i] = static_cast<double>(Cabsdata_f[i]);
  }

  // tell Python we are done with the input objects (so that memory is properly
  // deallocated)
  Py_DECREF(thetas);

  return PyArray_Squeeze(Cabs);
}

/**
 * @brief Compute the extinction matrix using the given T-matrix object.
 *
 * Required arguments are:
 *  - theta: Input zenith angle (in radians).
 *  - phi: Input azimuth angle (in radians).
 *
 * Additional optional arguments are:
 *  - alpha: Azimuth angle of the particle rotation axis (in radians).
 *  - beta: Zenith angle of the particle rotation axis (in radians).
 *
 * @param self T-matrix object being used.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Pointer to a 4x4 double precision NumPy array containing the
 * components of the extinction matrix.
 */
static PyObject *TmatrixObject_get_extinction_matrix(TmatrixObject *self,
                                                     PyObject *args,
                                                     PyObject *kwargs) {

  // required arguments
  PyArrayObject *thetas;
  PyArrayObject *phis;

  // optional arguments
  float_type alpha = 0.;
  float_type beta = 0.;

  // list of keywords (see comment above)
  static char *kwlist[] = {strdup("theta"), strdup("phi"), strdup("alpha"),
                           strdup("beta"), nullptr};

  // allocate temporary variables to store double precision arguments
  double alpha_d = static_cast<double>(alpha);
  double beta_d = static_cast<double>(beta);
  // parse positional and keyword arguments
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "O&O&|dd", kwlist, PyArray_Converter, &thetas,
          PyArray_Converter, &phis, &alpha_d, &beta_d)) {
    // again, we do not call ctm_error to avoid killing the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // this time, a nullptr return will signal an error to Python
    return nullptr;
  }
  // unpack double precision values
  alpha = alpha_d;
  beta = beta_d;

  // determine the size of the input arrays
  npy_intp thetasize;
  npy_intp phisize;

  // get the number of dimensions for the theta array
  const npy_intp thetadim = PyArray_NDIM(thetas);
  // we only accept 0D (scalar) or 1D input
  if (thetadim > 1) {
    ctm_warning("Wrong shape for input array!");
    return nullptr;
  }
  // check if we are in the 0D or 1D case
  if (thetadim > 0) {
    // if 1D, simply get the size from the 1 dimension
    const npy_intp *thetadims = PyArray_DIMS(thetas);
    thetasize = thetadims[0];
  } else {
    // if 0D, set the size to 1
    thetasize = 1;
    // reshape the array into a 1D array with 1 element, so that we can
    // manipulate it in the same way as a 1D array
    npy_intp newdims[1] = {1};
    PyArray_Dims newdimsobj;
    newdimsobj.ptr = newdims;
    newdimsobj.len = 1;
    thetas = reinterpret_cast<PyArrayObject *>(
        PyArray_Newshape(thetas, &newdimsobj, NPY_ANYORDER));
  }
  // now repeat for the phi array
  const npy_intp phidim = PyArray_NDIM(phis);
  if (phidim > 1) {
    ctm_warning("Wrong shape for input array!");
    return nullptr;
  }
  if (phidim > 0) {
    const npy_intp *phidims = PyArray_DIMS(phis);
    phisize = phidims[0];
  } else {
    phisize = 1;
    npy_intp newdims[1] = {1};
    PyArray_Dims newdimsobj;
    newdimsobj.ptr = newdims;
    newdimsobj.len = 1;
    phis = reinterpret_cast<PyArrayObject *>(
        PyArray_Newshape(phis, &newdimsobj, NPY_ANYORDER));
  }

  // create an uninitialised double NumPy array to store the results
  // shape: thetasize x phisize x 4 x 4
  npy_intp dims[4] = {thetasize, phisize, 4, 4};
  PyArrayObject *numpy_array =
      reinterpret_cast<PyArrayObject *>(PyArray_SimpleNew(4, dims, NPY_DOUBLE));

  // loop over all theta elements
  for (npy_intp itheta = 0; itheta < thetasize; ++itheta) {
    // get the corresponding theta angle
    const float_type theta =
        *(reinterpret_cast<double *>(PyArray_GETPTR1(thetas, itheta)));
    // loop over all phi elements
    for (npy_intp iphi = 0; iphi < phisize; ++iphi) {
      // get the corresponding phi angle
      const float_type phi =
          *(reinterpret_cast<double *>(PyArray_GETPTR1(phis, iphi)));

      // get the extinction matrix
      Matrix<float_type> K =
          self->_Tmatrix->get_extinction_matrix(alpha, beta, theta, phi);

      // copy the elements from the extinction matrix into the array
      for (uint_fast8_t i = 0; i < 4; ++i) {
        for (uint_fast8_t j = 0; j < 4; ++j) {
          *reinterpret_cast<double *>(PyArray_GETPTR4(
              numpy_array, itheta, iphi, i, j)) = static_cast<double>(K(i, j));
        }
      }
    }
  }

  // tell Python we are done with the input objects (so that memory is properly
  // deallocated)
  Py_DECREF(thetas);
  Py_DECREF(phis);

  // return the array
  // we squeeze it so that the output dimensions match the input arrays:
  // a 1D theta and 0D phi array e.g. will have a thetasize x 4 x 4 output
  return PyArray_Squeeze(numpy_array);
}

/*! @brief List of T-matrix object methods that are exposed to Python. */
static PyMethodDef TmatrixObject_methods[] = {
    {"get_nmax", reinterpret_cast<PyCFunction>(TmatrixObject_get_nmax),
     METH_NOARGS, "Return the maximum order of the T-matrix."},
    {"get_extinction_matrix",
     reinterpret_cast<PyCFunction>(TmatrixObject_get_extinction_matrix),
     METH_VARARGS | METH_KEYWORDS, "Return the extinction matrix."},
    {"get_average_absorption_cross_section",
     reinterpret_cast<PyCFunction>(
         TmatrixObject_get_average_absorption_cross_section),
     METH_NOARGS,
     "Return the angular average of the absorption cross section."},
    {"get_absorption_cross_section",
     reinterpret_cast<PyCFunction>(TmatrixObject_get_absorption_cross_section),
     METH_VARARGS | METH_KEYWORDS,
     "Return the absorption cross section for the given input angle(s)."},
    {"get_absorption_cross_sections",
     reinterpret_cast<PyCFunction>(TmatrixObject_get_absorption_cross_sections),
     METH_VARARGS | METH_KEYWORDS,
     "Return both absorption cross sections for the given input angle(s)."},
    {nullptr}};

/*! @brief Python Object type for the T-matrix (is edited in the module
 *  initialisation function since C++ does not allow fancy unordered struct
 *  initialisation). */
static PyTypeObject TmatrixObjectType = {PyVarObject_HEAD_INIT(nullptr, 0)};

/// MishchenkoOrientationDistributionObject

/**
 * @brief C struct wrapper around a MishchenkoOrientationDistribution object.
 */
typedef struct {
  /*! @brief Python object members. */
  PyObject_HEAD;
  /*! @brief T-matrix object. */
  MishchenkoOrientationDistribution *_distribution;
} MishchenkoOrientationDistributionObject;

/**
 * @brief Destructor for the MishchenkoOrientationDistributionObject.
 *
 * @param self MishchenkoOrientationDistributionObject that is being
 * deallocated.
 */
static void MishchenkoOrientationDistributionObject_dealloc(
    MishchenkoOrientationDistributionObject *self) {
  // first free the wrapped T-matrix object
  delete self->_distribution;
  // then free the Python object
  Py_TYPE(self)->tp_free(reinterpret_cast<PyObject *>(self));
}

/**
 * @brief Constructor for the MishchenkoOrientationDistributionObject.
 *
 * Constructs an empty MishchenkoOrientationDistributionObject.
 *
 * @param type Type for the object.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Pointer to the newly created object.
 */
static PyObject *MishchenkoOrientationDistributionObject_new(PyTypeObject *type,
                                                             PyObject *args,
                                                             PyObject *kwargs) {
  // allocate the Python object
  MishchenkoOrientationDistributionObject *self =
      (MishchenkoOrientationDistributionObject *)type->tp_alloc(type, 0);
  // set the wrapped T-matrix to a sensible initial value
  self->_distribution = nullptr;
  return reinterpret_cast<PyObject *>(self);
}

/**
 * @brief init() function for the MishchenkoOrientationDistributionObject.
 *
 * Required arguments are:
 *  - cos2beta: Parameter for the object.
 *
 * Additional optional arguments are:
 *  - maximum_order: Maximum order to use in the spherical basis function
 *    expansion.
 *
 * @param self MishchenkoOrientationDistributionObject that is being
 * initialised.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return 0 on success, 1 on failure.
 */
static int MishchenkoOrientationDistributionObject_init(
    MishchenkoOrientationDistributionObject *self, PyObject *args,
    PyObject *kwargs) {

  // required arguments
  float_type cos2beta;

  // optional arguments
  uint_fast32_t maximum_order = 200;

  /// parse arguments
  // list of keywords (in the expected order)
  // note that we need to use strdup because the Python API expects char*,
  // while C++ strings are const char*
  // not doing this results in compilation warnings
  static char *kwlist[] = {strdup("axis_ratio"), strdup("cos2beta"),
                           strdup("maximum_order"), nullptr};

  // allocate temporary variables to store double precision arguments
  double cos2beta_d;
  // parse the keywords/positional arguments
  // d is a double
  // I is an unsigned integer
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d|I", kwlist, &cos2beta_d,
                                   &maximum_order)) {
    // PyArg_ParseTupleAndKeywords will return 0 if a required argument was
    // missing, if an argument of the wrong type was provided or if the number
    // of arguments does not match the expectation
    // we use ctm_warning so that we can gracefully exit and let Python handle
    // the exception
    // ctm_error would call abort, which would kill the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // exit code 1 signals to Python that something was wrong
    return 1;
  }
  // unpack double precision arguments
  cos2beta = cos2beta_d;

  self->_distribution =
      new MishchenkoOrientationDistribution(maximum_order, cos2beta);

  // everything went well: return 0
  return 0;
}

/**
 * @brief Compute the orientation distribution for the given angle.
 *
 * Required arguments are:
 *  - beta: Zenith angle of the particle rotation axis (in radians).
 *
 * @param self MishchenkoOrientationDistributionObject being used.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Value of the distribution, wrapped in a Python float object.
 */
static PyObject *MishchenkoOrientationDistributionObject_get_distribution(
    MishchenkoOrientationDistributionObject *self, PyObject *args,
    PyObject *kwargs) {

  // required arguments
  float_type beta;

  // list of keywords (see comment above)
  static char *kwlist[] = {strdup("beta"), nullptr};

  // allocate temporary variables to store double precision arguments
  double beta_d;
  // parse positional and keyword arguments
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "d", kwlist, &beta_d)) {
    // again, we do not call ctm_error to avoid killing the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // this time, a nullptr return will signal an error to Python
    return nullptr;
  }
  // unpack double precision arguments
  beta = beta_d;

  const float_type value =
      (*static_cast<OrientationDistribution *>(self->_distribution))(beta);

  // return the array
  return PyFloat_FromDouble(static_cast<double>(value));
}

/*! @brief List of MishchenkoOrientationDistributionObject methods that are
 *  exposed to Python. */
static PyMethodDef MishchenkoOrientationDistributionObject_methods[] = {
    {"get_distribution",
     reinterpret_cast<PyCFunction>(
         MishchenkoOrientationDistributionObject_get_distribution),
     METH_VARARGS | METH_KEYWORDS, "Return the distribution function."},
    {nullptr}};

/*! @brief Python Object type for the MishchenkoOrientationDistributionObject
 *  (is edited in the module initialisation function since C++ does not allow
 *  fancy unordered struct initialisation). */
static PyTypeObject MishchenkoOrientationDistributionObjectType = {
    PyVarObject_HEAD_INIT(nullptr, 0)};

/// CosTuuM module

/**
 * @brief Get the equal volume radius for the particle with the given equal area
 * radius and axis ratio.
 *
 * Required arguments are:
 *  - radius: Equal area radius (in input length units).
 *  - axis_ratio: Axis ratio.
 *
 * @param self MishchenkoOrientationDistributionObject being used.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Equal volume radius, in same length units as input.
 */
static PyObject *get_equal_volume_radius(PyObject *self, PyObject *args,
                                         PyObject *kwargs) {

  // required arguments
  float_type radius;
  float_type axis_ratio;

  // list of keywords (see comment above)
  static char *kwlist[] = {strdup("radius"), strdup("axis_ratio"), nullptr};

  // allocate temporary variables to store double precision arguments
  double radius_d, axis_ratio_d;
  // parse positional and keyword arguments
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dd", kwlist, &radius_d,
                                   &axis_ratio_d)) {
    // again, we do not call ctm_error to avoid killing the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // this time, a nullptr return will signal an error to Python
    return nullptr;
  }
  // unpack double precision arguments
  radius = radius_d;
  axis_ratio = axis_ratio_d;

  // return the array
  return PyFloat_FromDouble(static_cast<double>(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
          axis_ratio) *
      radius));
}

/**
 * @brief Unpack the given 0D or 1D NumPy array into a std::vector with the
 * same number of elements.
 *
 * @param numpy_array Input NumPy array.
 * @return Output std::vector.
 * @tparam DATA_TYPE Data type of the elements in the result vector.
 * @tparam INPUT_TYPE Data type to use when reading values from the NumPy array.
 */
template <typename DATA_TYPE, typename INPUT_TYPE = DATA_TYPE>
inline std::vector<DATA_TYPE> unpack_numpy_array(PyArrayObject *numpy_array) {

  // determine the size of the array
  npy_intp array_size;
  const npy_intp array_ndim = PyArray_NDIM(numpy_array);
  if (array_ndim > 1) {
    ctm_warning("Wrong shape for input array!");
    return std::vector<DATA_TYPE>();
  }
  if (array_ndim > 0) {
    const npy_intp *array_dims = PyArray_DIMS(numpy_array);
    array_size = array_dims[0];
  } else {
    array_size = 1;
    npy_intp newdims[1] = {1};
    PyArray_Dims newdimsobj;
    newdimsobj.ptr = newdims;
    newdimsobj.len = 1;
    numpy_array = reinterpret_cast<PyArrayObject *>(
        PyArray_Newshape(numpy_array, &newdimsobj, NPY_ANYORDER));
  }

  std::vector<DATA_TYPE> array(array_size);
  for (npy_intp i = 0; i < array_size; ++i) {
    array[i] =
        *(reinterpret_cast<INPUT_TYPE *>(PyArray_GETPTR1(numpy_array, i)));
  }

  return array;
}

/**
 * @brief Get a table of absorption coefficients for the given input dust grain
 * sizes, dust grain types, wavelengths and angles, using the given shape and
 * alignment distributions.
 *
 * Required arguments are:
 *  - types: Array of dust grain types.
 *  - sizes: Array of dust grain sizes (in m).
 *  - wavelengths: Array of incoming/outgoing photon wavelengths (in m).
 *  - thetas: Array of zenith angles (in radians).
 *  - shape_distribution: Shape distribution object (needs proper
 *  implementation).
 *  - alignment_distribution: Alignment distribution object (needs proper
 *  implementation).
 *
 * Additional optional arguments are:
 *  - to be added
 *
 * @param self Module object.
 * @param args Positional arguments.
 * @param kwargs Keyword arguments.
 * @return Nothing.
 */
static PyObject *get_table(PyObject *self, PyObject *args, PyObject *kwargs) {

  // parse arguments //

  // required arguments
  PyArrayObject *input_types;
  PyArrayObject *input_sizes;
  PyArrayObject *input_wavelengths;
  PyArrayObject *input_thetas;

  // list of keywords
  static char *kwlist[] = {strdup("sizes"), strdup("types"),
                           strdup("wavelengths"), strdup("thetas"), nullptr};

  // parse positional and keyword arguments
  if (!PyArg_ParseTupleAndKeywords(
          args, kwargs, "O&O&O&O&", kwlist, PyArray_Converter, &input_types,
          PyArray_Converter, &input_sizes, PyArray_Converter,
          &input_wavelengths, PyArray_Converter, &input_thetas)) {
    // again, we do not call ctm_error to avoid killing the Python interpreter
    ctm_warning("Wrong arguments provided!");
    // this time, a nullptr return will signal an error to Python
    return nullptr;
  }

  uint_fast32_t shape_distribution_type = 2;
  ShapeDistribution *shape_distribution;
  if (shape_distribution_type == 0) {
    shape_distribution = new ShapeDistribution();
    shape_distribution->evaluate(20u);
  } else if (shape_distribution_type == 1) {
    shape_distribution = new DraineHensleyShapeDistribution(20u);
  } else if (shape_distribution_type == 2) {
    shape_distribution = new SingleShapeShapeDistribution(1.00001);
  }
  SizeBasedAlignmentDistribution alignment_distribution(1.e-5, 0, 100);
  TaskManager task_manager(10, 100, 2, 1.e-4, 1e10, *shape_distribution,
                           alignment_distribution);

  std::vector<int_fast32_t> types =
      unpack_numpy_array<int_fast32_t>(input_types);
  for (uint_fast32_t i = 0; i < types.size(); ++i) {
    task_manager.add_composition(types[i]);
  }
  Py_DECREF(input_types);

  std::vector<float_type> sizes =
      unpack_numpy_array<float_type, double>(input_sizes);
  for (uint_fast32_t i = 0; i < sizes.size(); ++i) {
    task_manager.add_size(sizes[i]);
  }
  Py_DECREF(input_sizes);

  std::vector<float_type> wavelengths =
      unpack_numpy_array<float_type, double>(input_wavelengths);
  for (uint_fast32_t i = 0; i < wavelengths.size(); ++i) {
    task_manager.add_wavelength(wavelengths[i]);
  }
  Py_DECREF(input_wavelengths);

  QuickSched quicksched(4, true, "CTMModule_taskgraph.log");

  std::vector<float_type> thetas =
      unpack_numpy_array<float_type, double>(input_thetas);
  Py_DECREF(input_thetas);

  std::vector<Task *> tasks;
  std::vector<Resource *> resources;
  ResultKey *result_key = nullptr;
  std::vector<Result *> results;
  TMatrixAuxiliarySpaceManager *space_manager = nullptr;
  task_manager.generate_tasks(thetas, 20, quicksched, tasks, resources,
                              result_key, results, space_manager);

  Py_BEGIN_ALLOW_THREADS;
  quicksched.execute_tasks();
  Py_END_ALLOW_THREADS;

  npy_intp result_dims[5] = {
      static_cast<npy_intp>(result_key->composition_size()),
      static_cast<npy_intp>(result_key->size_size()),
      static_cast<npy_intp>(result_key->wavelength_size()),
      static_cast<npy_intp>(thetas.size()), 2};
  PyArrayObject *result_array = reinterpret_cast<PyArrayObject *>(
      PyArray_SimpleNew(5, result_dims, NPY_DOUBLE));

  for (npy_intp itype = 0; itype < result_dims[0]; ++itype) {
    for (npy_intp isize = 0; isize < result_dims[1]; ++isize) {
      for (npy_intp ilambda = 0; ilambda < result_dims[2]; ++ilambda) {
        for (npy_intp itheta = 0; itheta < result_dims[3]; ++itheta) {
          const npy_intp result_index =
              result_key->get_result_index(0, isize, ilambda);
          const AbsorptionCoefficientResult &result =
              *static_cast<AbsorptionCoefficientResult *>(
                  results[result_index]);
          npy_intp Qabs_array_index[5] = {itype, isize, ilambda, itheta, 0};
          npy_intp Qabspol_array_index[5] = {itype, isize, ilambda, itheta, 1};
          *reinterpret_cast<double *>(
              PyArray_GetPtr(result_array, Qabs_array_index)) =
              static_cast<double>(result.get_Qabs(itheta));
          *reinterpret_cast<double *>(
              PyArray_GetPtr(result_array, Qabspol_array_index)) =
              static_cast<double>(result.get_Qabspol(itheta));
        } // theta
      }   // lambda
    }     // size
  }       // type

  return PyArray_Squeeze(result_array);
}

/*! @brief Methods exposed by the CTMmodule. */
static PyMethodDef CTMmethods[] = {
    {"get_equal_volume_radius",
     reinterpret_cast<PyCFunction>(get_equal_volume_radius),
     METH_VARARGS | METH_KEYWORDS,
     "Get the equal volume radius for the particle with the given equal area "
     "radius and axis ratio."},
    {"get_table", reinterpret_cast<PyCFunction>(get_table),
     METH_VARARGS | METH_KEYWORDS,
     "Test to see if the task-based algorithm can be coupled to the Python "
     "module."},
    {nullptr, nullptr, 0, nullptr}};

/*! @brief Module definition for the CTM module. */
static struct PyModuleDef CTMmodule = {PyModuleDef_HEAD_INIT, "CosTuuM",
                                       "T-matrix module.", -1, CTMmethods};

/**
 * @brief CosTuuM initialisation function.
 *
 * @return Pointer to the initialised module object.
 */
PyMODINIT_FUNC PyInit_CosTuuM() {

  // we need to call this to enable NumPy array functionality
  import_array();

  // create the module object
  PyObject *m = PyModule_Create(&CTMmodule);

  // set the required fields in the TmatrixObjectType struct
  // we need to do this because C++ does not support fancy C struct
  // initialisation
  // without this, we would need to manually set all elements in the struct
  // this is a pain (there are many) and makes it very hard to make the code
  // portable (addition of a single element to the struct would break it)
  TmatrixObjectType.tp_name = "CosTuuM.TMatrix";
  TmatrixObjectType.tp_basicsize = sizeof(TmatrixObject);
  TmatrixObjectType.tp_dealloc =
      reinterpret_cast<destructor>(TmatrixObject_dealloc);
  TmatrixObjectType.tp_methods = TmatrixObject_methods;
  TmatrixObjectType.tp_init = reinterpret_cast<initproc>(TmatrixObject_init);
  TmatrixObjectType.tp_new = TmatrixObject_new;

  // finalize creation of the TmatrixObjectType
  PyType_Ready(&TmatrixObjectType);
  // add a T-matrix object to the module
  // any call to CTMmodule.TMatrix() will use this same object
  Py_INCREF(&TmatrixObjectType);
  PyModule_AddObject(m, "TMatrix",
                     reinterpret_cast<PyObject *>(&TmatrixObjectType));

  // set the required fields in the MishchenkoOrientationDistributionObjectType
  // struct
  MishchenkoOrientationDistributionObjectType.tp_name =
      "CosTuuM.MishchenkoOrientationDistribution";
  MishchenkoOrientationDistributionObjectType.tp_basicsize =
      sizeof(MishchenkoOrientationDistributionObject);
  MishchenkoOrientationDistributionObjectType.tp_dealloc =
      (destructor)MishchenkoOrientationDistributionObject_dealloc;
  MishchenkoOrientationDistributionObjectType.tp_methods =
      MishchenkoOrientationDistributionObject_methods;
  MishchenkoOrientationDistributionObjectType.tp_init =
      (initproc)MishchenkoOrientationDistributionObject_init;
  MishchenkoOrientationDistributionObjectType.tp_new =
      MishchenkoOrientationDistributionObject_new;

  // finalize creation of the MishchenkoOrientationDistributionObjectType
  PyType_Ready(&MishchenkoOrientationDistributionObjectType);
  // add a MishchenkoOrientationDistributionObject to the module
  // any call to CTMmodule.MishchenkoOrientationDistribution() will use this
  // same object
  Py_INCREF(&MishchenkoOrientationDistributionObjectType);
  PyModule_AddObject(m, "MishchenkoOrientationDistribution",
                     reinterpret_cast<PyObject *>(
                         &MishchenkoOrientationDistributionObjectType));

  // return the module object
  return m;
}

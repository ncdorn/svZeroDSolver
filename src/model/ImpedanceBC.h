// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/**
 * @file ImpedanceBC.h
 * @brief model::ImpedanceBC source file
 */
#ifndef SVZERODSOLVER_MODEL_IMPEDANCEBC_HPP_
#define SVZERODSOLVER_MODEL_IMPEDANCEBC_HPP_

#include "Block.h"
#include "Model.h"
#include "Parameter.h"
#include "SparseSystem.h"

/**
 * @brief Pressure reference boundary condition.
 *
 * Applies a predefined pressure at a boundary.
 *
 * \f[
 * \begin{circuitikz}
 * \draw (1,0) node[anchor=south]{$P$} to [short, *-] (1.2,0) ;
 * \draw [-latex] (1.4,0) -- (2.2,0) node[right] {$Q$};
 * \draw (1,0) to [short, l=, *-] (1,-1)
 * node[ground]{$\hat{P}$};
 * \end{circuitikz}
 * \f]
 *
 * ### Governing equations
 *
 * \f[
 * P=\hat{P}
 * \f]
 *
 * ### Local contributions
 *
 * \f[
 * \mathbf{y}^{e}=\left[\begin{array}{ll}P & Q\end{array}\right]^{T}
 * \f]
 *
 * \f[
 * \mathbf{F}^{e}=\left[\begin{array}{ll}1 & 0\end{array}\right]
 * \f]
 *
 * \f[
 * \mathbf{C}^{e}=\left[\hat{P}\right]
 * \f]
 *
 * ### Parameters
 *
 * Parameter sequence for constructing this block
 *
 * * `0` Impedance
 *
 */

class ImpedanceBC : public Block {
 public:
  /**
   * @brief Construct a new ImpedanceBC object
   *
   * @param id Global ID of the block
   * @param model The model to which the block belongs
   */
  ImpedanceBC(int id, Model *model)
       : Block(id, model, BlockType::impedance_bc, BlockClass::boundary_condition,
              {{"t", InputParameter(false, true)},
               {"Z", InputParameter(false, true)},
               {"tree", InputParameter()},
               {"Pd", InputParameter()}}) {}
  
  bool return_array = true; // 

  /**
   * @brief Set up the degrees of freedom (DOF) of the block
   *
   * Set \ref global_var_ids and \ref global_eqn_ids of the element based on the
   * number of equations and the number of internal variables of the
   * element.
   *
   * @param dofhandler Degree-of-freedom handler to register variables and
   * equations at
   */
  void setup_dofs(DOFHandler &dofhandler);

  /**
   * @brief Update the constant contributions of the element in a sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_constant(SparseSystem &system, std::vector<double> &parameters, std::map<int, std::vector<double>> &parameter_arrays);

  /**
   * @brief Update the time-dependent contributions of the element in a sparse
   * system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   */
  void update_time(SparseSystem &system, std::vector<double> &parameters, std::map<int, std::vector<double>> &parameter_arrays);

  /**
   * @brief Update the solution of the element in a sparse system
   *
   * @param system System to update contributions at
   * @param parameters Parameters of the model
   * @param y Solution vector
   * @param dy Solution derivative vector
   */
  void update_solution(
      SparseSystem &system, std::vector<double> &parameters,
      std::map<int, std::vector<double>> &parameter_arrays,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
      const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy,
      bool &converged);

  /**
   * @brief Number of triplets of element
   *
   * Number of triplets that the element contributes to the global system
   * (relevant for sparse memory reservation)
   */
  TripletsContributions num_triplets{1, 0, 0};


 private:
   double zq_conv; // z and q convolved over the previous period
   double z_0; // impedance at the beginning of the period
   int num_timesteps; // number of timesteps in the period
   std::vector<double> z; // impedance at each timestep
   std::vector<double> q; // flow at each timestep
   std::vector<double> times_1per; // vector of times for debugging
   // std::vector<double> times_2per; // vector of times for debugging
   std::vector<double> z_interp;
   std::vector<double> times; // vector of times for debugging
   std::vector<double> zq_conv_vec;



   /**
   * @brief get the convolution of z and q
   *
   * @param parameters Parameters of the model
   * @param parameter_arrays Arrays of time-dependent model parameters
   */
   void convolve_zq(std::vector<double> &parameters, std::map<int, std::vector<double>> &parameter_arrays);

   /**
   * @brief Print the first five elements of a vector
   * 
   * @param vec Vector to print
   */
   void printvec(std::vector<double> &vec, std::string name, int n);

   /**
    * @brief Interpolate the impedance based on the actual simulation timesteps
    * 
    * @param parameter_arrays Parameters of the model
    */
   void interpolate_zq(double time, std::map<int, std::vector<double>> &parameter_arrays);
   
   /**
    *     * @brief Write a value to a file
    *         * 
    *             * @param value Value to write
    *                 * @param filename Name of the file to write to
    *                     */
   void writevalue(double &value, std::string filename);

    /**
      * @brief Write an array to a file
      * 
      * @param vec Vector to write
      * @param filename Name of the file to write to
      */
   void writearray(std::vector<double> &vec, std::string filename);

};

#endif  // SVZERODSOLVER_MODEL_IMPEDANCEBC_HPP_

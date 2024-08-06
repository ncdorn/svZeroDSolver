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

#include "ImpedanceBC.h"
#include "Model.h"

void ImpedanceBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 1, {});
}

void ImpedanceBC::update_constant(SparseSystem &system,
                                          std::vector<double> &parameters) {
  // eqn 0: P_in - zq_conv = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  // system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = z_0;
}

void ImpedanceBC::update_time(SparseSystem &system,
                                      std::vector<double> &parameters) {

  // eqn 0: P_in - zq_conv = 0
  
  // printf("zq_conv: %f\n", zq_conv);

  // convolve_zq(parameters);
  
  system.C(global_eqn_ids[0]) = -zq_conv;

}

void ImpedanceBC::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy,
    bool &converged) {
  
  q.push_back(y[global_var_ids[1]]); // add the new value to the q array

  convolve_zq(parameters, converged);
  
  system.C(global_eqn_ids[0]) = -zq_conv;

  if (!converged) {
    q.pop_back(); // remove the most recent value in the q array
  }

  if (model->time >= model->cardiac_cycle_period) {
    if (converged) {
      // printf("updating solution, converged \n");
      q.erase(q.begin()); // remove the first element in the q array (oldest value)
      printf("solution converged, q.size(): %d, z.size(): %d\n", q.size(), z.size());
    } else {
      // q.pop_back(); // remove the most recent value in the q array
      printf("solution NOT converged, q.size(): %d, z.size(): %d\n", q.size(), z.size());
    };
  
  }
  // printf("size of q: %d, size of z: %d \n", q.size(), z.size());
    
  }


void ImpedanceBC::setup_model_dependent_params() {
  z_0 = model->get_parameter_value(global_param_ids[0]); // I do not think we need this!
  // num_timesteps = model->sim
}

void ImpedanceBC::convolve_zq(std::vector<double> &parameters, bool &converged) {
  
  auto T_cardiac = model->cardiac_cycle_period;
  auto t = model->time;

  // if t < T_per, zq_conv = 0
  zq_conv = 0.0;
  // if t >= T_per, zq_conv = \sum_k=1^N-1(q[(n-k)%N] * z[k])

  if (t < T_cardiac) {
    z.push_back(parameters[global_param_ids[0]]);
    if (!converged) {
      z.pop_back();
    }
    } else {
    std::reverse(q.begin(), q.end()); // reverse the vector
    // int N =  // number of time steps in the period

    for (int k = 0; k < z.size(); ++k) {
      zq_conv += q[k] * z[k];
      // if (k % 100 == 0) {
      //   printf("k: %d, q[k]: %f, z[k]: %f, zq_conv: %f\n", k, q[k], z[k], zq_conv);
      // }
    };
    std::reverse(q.begin(), q.end()); // reverse back
  };
  // where N = number of time steps in the period
  //       n = current time step

}

void ImpedanceBC::printfive(std::vector<double> &vec) {
  // printf("q_rev[0:5] = [");
  // for (int i = 0; i < 5; ++i) {
  //   printf("%f,", i, vec[i]);
  // }
  // printf("]\n");

  // count zeros in the vector
  int count = 0;
  for (int i = 0; i < vec.size(); ++i) {
    if (abs(vec[i]) <= 0.00001) {
      count++;
    }
  }


  // check the integral of flow, to make sure it is staying the same?
  double sum = 0;
  for (int i = 0; i < vec.size(); ++i) {
    sum += vec[i];
  }
  printf("sum of q: %f\n", sum);
}

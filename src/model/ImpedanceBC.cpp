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
  // eqn 0: P_in - zq_conv - P_d = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  // system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = z_0;
}

void ImpedanceBC::update_time(SparseSystem &system,
                                      std::vector<double> &parameters,
                                      std::map<int, std::vector<double>> &parameter_arrays) {
  convolve_zq(parameters, parameter_arrays);
  if (model->time <= model->cardiac_cycle_period) {
    times_1per.push_back(model->time);
  }
  // else if (model->time < model->cardiac_cycle_period * 2) {
  //   times_2per.push_back(fmod(model->time, model->cardiac_cycle_period));

  //   printf("times_2per[0:10]: ");
  //   for (int i = 0; i < 10; i++) {
  //     printf("%f, ", times_2per[i]);
  //   }
  //   printf("\n");
  //   printf("times_2per.size(): %d\n", times_2per.size());
  // }

  // printf("times_1per[0:10]: ");
  // for (int i = 0; i < 10; i++) {
  //   printf("%f, ", times_1per[i]);
  // }
  // printf("\n");

  if (model->time > model->cardiac_cycle_period) {
    interpolate_zq(parameter_arrays);
  }

}

void ImpedanceBC::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy,
    bool &converged) {
  
  q.push_back(y[global_var_ids[1]]); // add the new value to the q array

  system.C(global_eqn_ids[0]) = -zq_conv - parameters[global_param_ids[1]];

  if (!converged) {
    // printf("solution not converged for q = %f\n", y[global_var_ids[1]]);
    q.pop_back(); // remove the most recent value in the q array
  }
  else {
    // printf("solution converged for q = %f\n", y[global_var_ids[1]]);
  }

  if (model->time > model->cardiac_cycle_period) {
    if (converged) {
      q.erase(q.begin()); // remove the first element in the q array (oldest value)
      // printfive(q);
    } else {
      // q.pop_back(); // remove the most recent value in the q array
    };
  
  }
    
  }

void ImpedanceBC::convolve_zq(std::vector<double> &parameters, std::map<int, std::vector<double>> &parameter_arrays) {
  
  auto T_cardiac = model->cardiac_cycle_period;
  auto t = model->time;
  
  std::vector<double> z = parameter_arrays[global_param_ids[0]];

  // printf("size of z: %d | size of q: %d\n", z.size(), q.size());

  // TODO: PRINT THINGE FROM HERE TO FIGURE OUT HOW TO BEST ACCESS THE PARAMETER ARRAY
  // AND THEN FIGURE OUT HOW TO SAVE THE FLOW RESULT FROM THE 3D SIMULATION TO CONVOLVE STUFF

  if (t < T_cardiac) {
    zq_conv = parameters[global_param_ids[0]];
  } else {
  
  std::reverse(q.begin(), q.end()); // reverse the vector
  // int N =  // number of time steps in the period

  // printf("q, z with NEW CODE\n");
  // printf("q_rev[0]: %f, q_rev[1]: %f, q_rev[2]: %f, q_rev[3]: %f, q_rev[4]: %f \n", q[0], q[1], q[2], q[3], q[4]);
  // printf("z[0]: %f, z[1]: %f, z[2]: %f, z[3]: %f, z[4]: %f \n", z[0], z[1], z[2], z[3], z[4]);
  

  for (int k = 0; k < z.size(); ++k) {
    zq_conv += q[k] * z[k]; // NEED TO GET TIMESTEP AND MULTIPLY BY THIS
    // if (k % 100 == 0) {
    //   printf("k: %d, q[k]: %f, z[k]: %f, zq_conv: %f\n", k, q[k], z[k], zq_conv);
    // }
  };
  std::reverse(q.begin(), q.end()); // reverse back

  float per = t / T_cardiac;

  printf("zq_conv for per = %f: %f \n", per, zq_conv);

  // where N = number of time steps in the period
  //       n = current time step
  }

}

void ImpedanceBC::interpolate_zq(std::map<int, std::vector<double>> &parameter_arrays) {

  double T_cardiac = model->cardiac_cycle_period;
  double t = model->time;

  std::vector<double> z = parameter_arrays[global_param_ids[0]];

  double tstep_size = T_cardiac / z.size();


  for (int i = 0; i < z.size(); ++i) {
    double t = i * tstep_size;

    z_times.push_back(t);
  }

  for (int i = 0; i < z.size(); ++i) {
    z_interp.push_back(z[i] + (z[i+1] - z[i]) * (times_1per[i] - z_times[i]) / (z_times[i+1] - z_times[i]));
  }
  
}

void ImpedanceBC::printfive(std::vector<double> &vec) {
  printf("z_interp[0:5] = [");
  for (int i = 0; i < 5; ++i) {
    printf("%f,", vec[i]);
  }
  printf("]\n");

  // printf("for z of len %d, z[:10] = [", z.size());
  // int len = z.size();
  // for (int i = 0; i < 10; ++i) {
  //   printf("%f,", vec[i]);
  // };
  // printf("]\n");

  // for (int i = 0; i < z.size(); ++i) {
  //   if (trunc(z[i]) == 1160) {
  //   // printf("z(t = 0) is in this array! at index %d \n", i);
  //   }
  // };
}

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

// for debugging
#include <iostream>
#include <fstream>
void ImpedanceBC::setup_dofs(DOFHandler &dofhandler) {
  Block::setup_dofs_(dofhandler, 1, {});
}

void ImpedanceBC::update_constant(SparseSystem &system,
                                  std::vector<double> &parameters,
                                  std::map<int, std::vector<double>> &parameter_arrays) {
  // eqn 0: P_in - zq_conv - P_d = 0
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[0]) = 1.0;
  // system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = z_0;

}

void ImpedanceBC::update_time(SparseSystem &system,
                                      std::vector<double> &parameters,
                                      std::map<int, std::vector<double>> &parameter_arrays) {
  
  if (model->time <= model->cardiac_cycle_period) {
    times_1per.push_back(model->time);
  }

  /// TESTING OUT ACTUAL DISCRETIZATION? ///
  system.F.coeffRef(global_eqn_ids[0], global_var_ids[1]) = -parameter_arrays[global_param_ids[0]][0];

  system.C(global_eqn_ids[0]) = -zq_conv - parameters[global_param_ids[1]];


}

void ImpedanceBC::update_solution(
    SparseSystem &system, std::vector<double> &parameters,
    std::map<int, std::vector<double>> &parameter_arrays,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &y,
    const Eigen::Matrix<double, Eigen::Dynamic, 1> &dy,
    bool &converged) {

  // we need to make sure the simulation is unique for each timestep

  q.push_back(y[global_var_ids[1]]); // add the new value to the q array

  if (times.size() == 0) {
    times.push_back(0.0);
  } else {
    if (model->time > times.back()) {
      // std::cout << "we are on a new timestep!" << std::endl;
      times.push_back(model->time);
    }
  }

  if (q.size() > times.size()) {
    // each q needs to be unique for each t
    if (q.size() >= 2) {  // Ensure the vector has at least two elements
        q.erase(q.end() - 2);  // Erase the second-to-last element
    }
    // q.pop_back();
  }

  convolve_zq(parameters, parameter_arrays);
  // zq_conv_vec.push_back(zq_conv);

  // THIS NEXT SECTION IS FOR DEBUGGING ZQ CONV VEC
  // if (zq_conv_vec.size() > times.size()) {
  //   // each q needs to be unique for each t
  //   if (zq_conv_vec.size() >= 2) {  // Ensure the vector has at least two elements
  //       zq_conv_vec.erase(zq_conv_vec.end() - 2);
  //   }
  //   // q.pop_back();
  // }


  if (!converged) {
    // printf("solution not converged for q = %f\n", y[global_var_ids[1]]);
    q.pop_back(); // remove the most recent value in the q array
    // zq_conv_vec.pop_back();
  }
  else {
    double time_diff = model->time - (2 * model->cardiac_cycle_period) - 0.5;

    if (time_diff > 0 && time_diff < 0.001) {
      // std::cout << time_diff << std::endl;
      // std::cout << "time: " << model->time << std::endl;
      std::string qfile = "debug/q.txt";
      writearray(q, qfile);

      std::string zfile = "debug/z_interp.txt";
      writearray(z_interp, zfile);

      std::string zqfile = "debug/zq_conv.txt";
      writearray(zq_conv_vec, zqfile);

    }

  }

  if (model->time > model->cardiac_cycle_period) {
    if (converged) {
      // THIS MIGHT BE THE BUG...
      // q.erase(q.begin()); // remove the first element in the q array (oldest value)
      // printfive(q);
    } else {
      // q.pop_back(); // remove the most recent value in the q array
    };
  }
    
}

void ImpedanceBC::convolve_zq(std::vector<double> &parameters, std::map<int, std::vector<double>> &parameter_arrays) {
  
  auto T_cardiac = model->cardiac_cycle_period;
  double t = model->time;
  double tstep = times.rbegin()[0] - times.rbegin()[1];

  // std::cout << "convolving z and q for t = " << t << std::endl;
  
  std::vector<double> z = parameter_arrays[global_param_ids[0]];

  if (t <= model->cardiac_cycle_period) {
    if (z_interp.size() < times.size()) {
      if (z_interp.size() == 0) {
        z_interp.push_back(parameter_arrays[global_param_ids[0]][0]);

      } else {
        interpolate_zq(t, parameter_arrays);
      }
    }
  }
  // print out z and q
  
//  printvec(times, t_name, times.size());

  // printf("size of z: %d | size of q: %d\n", z.size(), q.size());

  // TODO: PRINT THINGE FROM HERE TO FIGURE OUT HOW TO BEST ACCESS THE PARAMETER ARRAY
  // AND THEN FIGURE OUT HOW TO SAVE THE FLOW RESULT FROM THE 3D SIMULATION TO CONVOLVE STUFF

  // if (t < T_cardiac) {

  //   zq_conv = 0.0;

  //   double system_c = -parameters[global_param_ids[1]];

  //   // std::string debugcfile = "debug/system_c.txt";
  //   // writevalue(system_c, debugcfile);

  //   // std::cout << "global_param_ids[0]" << parameters[global_param_ids[0]] << std::endl;
  //   // std::cout << "global_param_ids[1]" << parameters[global_param_ids[1]] << std::endl;
  // } else {


  // std::cout << "convolving z and q for tstep = " << tstep << std::endl;

  std::string imp = "impedance interpolated";
  // printvec(z_interp, imp, z_interp.size());
  std::string q_name = "Q";
  // printvec(q, q_name, q.size());

  std::reverse(q.begin(), q.end()); // reverse the vector
  // int N =  // number of time steps in the period

  float per = t / T_cardiac;

  // std::cout << "z_interp size: " << z_interp.size() << std::endl;
  // std::cout << "q size: " << q.size() << std::endl;
  zq_conv = 0;

  if (z_interp.size() < q.size()) {
    for (int k = 1; k < z_interp.size(); ++k) {
      zq_conv += q[k] * z_interp[k]; // NEED TO GET TIMESTEP AND MULTIPLY BY THIS
      // if (k % 100 == 0) {
      //   printf("k: %d, q[k]: %f, z[k]: %f, zq_conv: %f\n", k, q[k], z[k], zq_conv);
      // }
    };
  } else {
    for (int k = 1; k < q.size(); ++k) {
      zq_conv += q[k] * z_interp[k]; // NEED TO GET TIMESTEP AND MULTIPLY BY THIS
      // if (k % 100 == 0) {
      //   printf("k: %d, q[k]: %f, z[k]: %f, zq_conv: %f\n", k, q[k], z[k], zq_conv);
      // }
    };
  }

  // std::cout << "q_guess: " << q[0] << std::endl;

  std::reverse(q.begin(), q.end()); // reverse back


  // }

}

void ImpedanceBC::interpolate_zq(double time, std::map<int, std::vector<double>> &parameter_arrays) {

  double T_cardiac = model->cardiac_cycle_period;

  std::vector<double> z = parameter_arrays[global_param_ids[0]];

  int tstep = floor(time / T_cardiac * (z.size() - 1));


  // get the next and previous timesteps for interpolation
  // double tstep_size = 1 / z.size() * T_cardiac
  double t_last = tstep / (double)z.size() * T_cardiac;
  double t_next = (tstep + 1) / (double)z.size() * T_cardiac;

  double z_interped = z[tstep] + (z[tstep+1] - z[tstep]) * (time - t_last) / (t_next - t_last);

  /////// DEBUGGING ///////
  if (z_interped != z_interped) {
    std::cout << "z_interped is nan!" << std::endl;
    std::cout << "tstep: " << tstep << std::endl;
    std::cout << "t_last: " << t_last << std::endl;
    std::cout << "t_next: " << t_next << std::endl;
    std::cout << "time: " << time << std::endl;
    std::cout << "z[tstep]: " << z[tstep] << std::endl;
    std::cout << "z[tstep+1]: " << z[tstep+1] << std::endl;
    std::cout << "z.size(): " << z.size() << std::endl;
    std::cout << "T_cardiac: " << T_cardiac << std::endl;
  }

  z_interp.push_back(z[tstep] + (z[tstep+1] - z[tstep]) * (time - t_last) / (t_next - t_last));
  
}

void ImpedanceBC::printvec(std::vector<double> &vec, std::string name, int n) {

  if (n > vec.size()) {
    n = vec.size();
  }
  std::cout << name << "[0:" << n << "] = [";
  for (int i = vec.size() - n; i < vec.size(); ++i) {
    printf("%f ", vec[i]);
  }
  std::cout << "]" << std::endl;

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

void ImpedanceBC::writevalue(double &value, std::string filename) {
  std::ofstream file;
  file.open(filename, std::ios::app);
  file << value << std::endl;
  file.close();
}

void ImpedanceBC::writearray(std::vector<double> &vec, std::string filename) {
  std::ofstream file;
  file.open(filename, std::ios::app);
  for (int i = 0; i < vec.size(); ++i) {
    file << vec[i] << std::endl;
  }
  file.close();
}

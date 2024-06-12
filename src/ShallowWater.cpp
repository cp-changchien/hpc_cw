#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../include/ShallowWater.h"
#include "/opt/homebrew/opt/openblas/include/cblas.h"


// Progress bar function
void Progress_Bar(int& iteration, const int& timestep){
    double progress = static_cast<double>(iteration) / timestep;
    const int barWidth = 55;
    int pos = static_cast<int>(barWidth * progress);

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "#";
        else if (i == pos)
            std::cout << "#";
        else
            std::cout << " ";
    }
    
    std::cout << "] " << std::fixed << std::setprecision(3)
                    << progress * 100.0 << "%\r";
    std::cout.flush();
}


/**
 * @brief Takes the command line input as parsed from the boost_option function
 * and store their values in the class. Allocates memory to solutions fields.
 *
 * @param arg_dt    The time step for integration
 * @param arg_T     The total time of integration
 * @param arg_Nx    The number of grid points for x
 * @param arg_Ny    The number of grid points for y
 * @param arg_ic    Options/Index for initial conditions
 */
void ShallowWater::SetParameters(const double &arg_dt, const int &arg_T,
                                 const int &arg_Nx, const int &arg_Ny,
                                 const int &arg_ic, const int &arg_choice) {

  dt = arg_dt;
  T = arg_T;
  Nx = arg_Nx;
  Ny = arg_Ny;
  ic = arg_ic;
  choice = arg_choice;

  u = new double[Nx * Ny];
  v = new double[Nx * Ny];
  h = new double[Nx * Ny];
  dudx = new double[Nx * Ny];
  dudy = new double[Nx * Ny];
  dvdx = new double[Nx * Ny];
  dvdy = new double[Nx * Ny];
  dhdx = new double[Nx * Ny];
  dhdy = new double[Nx * Ny];
  u_next = new double[Nx * Ny];
  v_next = new double[Nx * Ny];
  h_next = new double[Nx * Ny];

  time_step = T / dt;
  g = 9.81;
  dx = 1.0;
  dy = 1.0;

  // Printing values of all the parameters to the terminal.
  std::cout << std::endl;
  std::cout << "List of Commands Line Input for PDE:" << std::endl;
  std::cout << "    - dt       :  " << dt << std::endl;
  std::cout << "    - T        :  " << T << std::endl;
  std::cout << "    - Nx       :  " << Nx << std::endl;
  std::cout << "    - Ny       :  " << Ny << std::endl;
  std::cout << "    - option   :  " << ic << std::endl;
  std::cout << "    - Method   :  " << choice << std::endl;
  std::cout << std::endl;
};

/**
 * @brief Set the initial condition and boundary conditions. where u(x,y,0) =
 * v(x,y,0) = 0 and the h(x,y,0) = g(x,y). g(x,y) is the function prescribing
 * the initial surface height where the user can select options from the command
 * line input (option 1-4)
 *
 * Periodic boundary conditions are also used on all boundaries, so that waves
 * which propagate out of one of the boundaries re-enter on the opposite
 * boundary.
 *
 */
void ShallowWater::SetInitialConditions() {
  // 'i' is used for the column index, 'j' for the row index,
  // 'Nx' is the number of colume and 'Ny' the number of rows.
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      double x = i * dx;
      double y = j * dy;

      u[i * Ny + j] = 0.0;
      v[i * Ny + j] = 0.0;

      switch (ic) {
      case 1:
        h[i * Ny + j] = 10.0 + exp(-(x - 50.0) * (x - 50.0) / 25.0);
        break;
      case 2:
        h[i * Ny + j] = 10.0 + exp(-(y - 50.0) * (y - 50.0) / 25.0);
        break;
      case 3:
        h[i * Ny + j] =
            10.0 +
            exp(-((x - 50.0) * (x - 50.0) + (y - 50.0) * (y - 50.0)) / 25.0);
        break;
      case 4:
        h[i * Ny + j] = 10.0 +
                        exp(-(pow(x - 25.0, 2.0) + pow(y - 25.0, 2.0)) / 25.0) +
                        exp(-(pow(x - 75.0, 2.0) + pow(y - 75.0, 2.0)) / 25.0);
        break;
      }
    }
  }
};

/**
 * @brief Calculate the partial derivative terms using the 6th-order central
 * difference stencil. Applying periodic boundary conditions to equalise
 * opposite sides
 *
 * @param u
 * @param dudx
 * @param dudy
 * @param choice
 */
void ShallowWater::Partial_Derivative(double *u, double *dudx, double *dudy,
                                      int choice) {


    #pragma omp parallel
    {
    // Calculate the derivative wrt to x
    // Pre-defining the boundary for x (right and left three columns)
        #pragma omp for nowait
        for (int j = 0; j < Ny; ++j) {

            dudx[0 * Ny + j] = (1.0 / dx) * ((-1.0 / 60.0) * u[(Nx - 3) * Ny + j] + (3.0 / 20.0) * u[(Nx - 2) * Ny + j] - 
                                (3.0 / 4.0) * u[(Nx - 1) * Ny + j] + (3.0 / 4.0) * u[1 * Ny + j] - (3.0 / 20.0) * u[2 * Ny + j] + (1.0 / 60.0) * u[3 * Ny + j]);
            dudx[1 * Ny + j] = (1.0 / dx) * ((-1.0 / 60.0) * u[(Nx - 2) * Ny + j] + (3.0 / 20.0) * u[(Nx - 1) * Ny + j] - 
                                (3.0 / 4.0) * u[0 * Ny + j] + (3.0 / 4.0) * u[2 * Ny + j] - (3.0 / 20.0) * u[3 * Ny + j] + (1.0 / 60.0) * u[4 * Ny + j]);
            dudx[2 * Ny + j] = (1.0 / dx) * ((-1.0 / 60.0) * u[(Nx - 1) * Ny + j] + (3.0 / 20.0) * u[0 * Ny + j] - (3.0 / 4.0) * u[1 * Ny + j] + 
                                (3.0 / 4.0) * u[3 * Ny + j] - (3.0 / 20.0) * u[4 * Ny + j] + (1.0 / 60.0) * u[5 * Ny + j]);

            dudx[(Nx - 3) * Ny + j] = (1.0 / dx) * ((-1.0 / 60.0) * u[(Nx - 6) * Ny + j] + (3.0 / 20.0) * u[(Nx - 5) * Ny + j] - 
                                        (3.0 / 4.0) * u[(Nx - 4) * Ny + j] + (3.0 / 4.0) * u[(Nx - 2) * Ny + j] - (3.0 / 20.0) * u[(Nx - 1) * Ny + j] + (1.0 / 60.0) * u[0 * Ny + j]);
            dudx[(Nx - 2) * Ny + j] = (1.0 / dx) * ((-1.0 / 60.0) * u[(Nx - 5) * Ny + j] + (3.0 / 20.0) * u[(Nx - 4) * Ny + j] - 
                                        (3.0 / 4.0) * u[(Nx - 3) * Ny + j] + (3.0 / 4.0) * u[(Nx - 1) * Ny + j] - (3.0 / 20.0) * u[0 * Ny + j] + (1.0 / 60.0) * u[1 * Ny + j]);
            dudx[(Nx - 1) * Ny + j] = (1.0 / dx) * ((-1.0 / 60.0) * u[(Nx - 4) * Ny + j] + (3.0 / 20.0) * u[(Nx - 3) * Ny + j] - 
                                        (3.0 / 4.0) * u[(Nx - 2) * Ny + j] + (3.0 / 4.0) * u[0 * Ny + j] - (3.0 / 20.0) * u[1 * Ny + j] + (1.0 / 60.0) * u[2 * Ny + j]);
        }

        // calculate the middle regions
        #pragma omp for nowait collapse(2)
        for (int i = 3; i < (Nx - 3); ++i) {
            for (int j = 0; j < Ny; ++j) {
                // calculate du/dx
                dudx[i * Ny + j] = (1.0 / dx) * ((-1.0 / 60.0) * u[(i - 3) * Ny + j] + (3.0 / 20.0) * u[(i - 2) * Ny + j] - (3.0 / 4.0) * u[(i - 1) * Ny + j] + 
                                    (3.0 / 4.0) * u[(i + 1) * Ny + j] - (3.0 / 20.0) * u[(i + 2) * Ny + j] + (1.0 / 60.0) * u[(i + 3) * Ny + j]);
            }
        }

        // Calculate the derivative wrt to y
        // Pre-defining the boundary for y (top and bottom three rows)
        #pragma omp for nowait
        for (int i = 0; i < Nx; ++i) {

            dudy[i * Ny + 0] = (1.0 / dy) * ((-1.0 / 60.0) * u[i * Ny + (Ny - 3)] + (3.0 / 20.0) * u[i * Ny + (Ny - 2)] -
                                (3.0 / 4.0) * u[i * Ny + (Ny - 1)] + (3.0 / 4.0) * u[i * Ny + 1] - (3.0 / 20.0) * u[i * Ny + 2] + (1.0 / 60.0) * u[i * Ny + 3]);
            dudy[i * Ny + 1] = (1.0 / dy) * ((-1.0 / 60.0) * u[i * Ny + (Ny - 2)] + (3.0 / 20.0) * u[i * Ny + (Ny - 1)] - (3.0 / 4.0) * u[i * Ny + 0] + 
                                (3.0 / 4.0) * u[i * Ny + 2] - (3.0 / 20.0) * u[i * Ny + 3] + (1.0 / 60.0) * u[i * Ny + 4]);
            dudy[i * Ny + 2] = (1.0 / dy) * ((-1.0 / 60.0) * u[i * Ny + (Ny - 1)] + (3.0 / 20.0) * u[i * Ny + 0] - (3.0 / 4.0) * u[i * Ny + 1] +
                                (3.0 / 4.0) * u[i * Ny + 3] - (3.0 / 20.0) * u[i * Ny + 4] + (1.0 / 60.0) * u[i * Ny + 5]);

            dudy[i * Ny + (Ny - 3)] = (1.0 / dy) * ((-1.0 / 60.0) * u[i * Ny + (Ny - 6)] + (3.0 / 20.0) * u[i * Ny + (Ny - 5)] - (3.0 / 4.0) * u[i * Ny + (Ny - 4)] + 
                                        (3.0 / 4.0) * u[i * Ny + (Ny - 2)] - (3.0 / 20.0) * u[i * Ny + (Ny - 1)] + (1.0 / 60.0) * u[i * Ny + 0]);
            dudy[i * Ny + (Ny - 2)] = (1.0 / dy) * ((-1.0 / 60.0) * u[i * Ny + (Ny - 5)] + (3.0 / 20.0) * u[i * Ny + (Ny - 4)] - 
                                        (3.0 / 4.0) * u[i * Ny + (Ny - 3)] + (3.0 / 4.0) * u[i * Ny + (Ny - 1)] - (3.0 / 20.0) * u[i * Ny + 0] + (1.0 / 60.0) * u[i * Ny + 1]);
            dudy[i * Ny + (Ny - 1)] = (1.0 / dy) * ((-1.0 / 60.0) * u[i * Ny + (Ny - 4)] + (3.0 / 20.0) * u[i * Ny + (Ny - 3)] - 
                                        (3.0 / 4.0) * u[i * Ny + (Ny - 2)] + (3.0 / 4.0) * u[i * Ny + 0] - (3.0 / 20.0) * u[i * Ny + 1] + (1.0 / 60.0) * u[i * Ny + 2]);
        }

        // calculate the middle regions
        #pragma omp for nowait collapse(2)
        for (int i = 0; i < Nx; ++i) {
            for (int j = 3; j < (Ny - 3); ++j) {

            // calculate du/dy
            dudy[i * Ny + j] = (1.0 / dy) * ((-1.0 / 60.0) * u[i * Ny + j - 3] + (3.0 / 20.0) * u[i * Ny + j - 2] - (3.0 / 4.0) * u[i * Ny + j - 1] + 
                                    (3.0 / 4.0) * u[i * Ny + j + 1] - (3.0 / 20.0) * u[i * Ny + j + 2] + (1.0 / 60.0) * u[i * Ny + j + 3]);
            }
        }
    }
};





/**
 * @brief For the time-integration using the 4th-order Runge-Kutta explicit
 * scheme.
 *
 */
void ShallowWater::TimeIntegrate() {

  double *k1u = new double[Nx * Ny];
  double *k2u = new double[Nx * Ny];
  double *k3u = new double[Nx * Ny];
  double *k4u = new double[Nx * Ny];
  double *k1v = new double[Nx * Ny];
  double *k2v = new double[Nx * Ny];
  double *k3v = new double[Nx * Ny];
  double *k4v = new double[Nx * Ny];
  double *k1h = new double[Nx * Ny];
  double *k2h = new double[Nx * Ny];
  double *k3h = new double[Nx * Ny];
  double *k4h = new double[Nx * Ny];

  double *temp_u = new double[Nx * Ny];
  double *temp_v = new double[Nx * Ny];
  double *temp_h = new double[Nx * Ny];

  std::cout << "Start of Numerical Solving PDE.\n" << std::endl;
  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();

  for (int t = 0; t < time_step; t++) {

    // Start of pragma parallelisation region
    // #pragma omp parallel
    // {

    // k1 = f(yn)
    Partial_Derivative(u, dudx, dudy, choice);

    Partial_Derivative(v, dvdx, dvdy, choice);

    Partial_Derivative(h, dhdx, dhdy, choice);

    // k1u = fu(u, v, dudx, dudy, dhdx);
    // k1v = fv(u, v, dvdx, dvdy, dhdy);
    // k1h = fh(u, v, h, dudx, dvdy, dhdx, dhdy);

#pragma omp parallel for collapse(2)
    for (int i = 0; i < Nx; ++i) {
      for (int j = 0; j < Ny; ++j) {
        int index = i * Ny + j;
        k1u[index] = -u[index] * dudx[index] - v[index] * dudy[index] - g * dhdx[index];
        k1v[index] = -u[index] * dvdx[index] - v[index] * dvdy[index] - g * dhdy[index];
        k1h[index] = -u[index] * dhdx[index] - h[index] * dudx[index] - v[index] * dhdy[index] - h[index] * dvdy[index];

        temp_u[index] = u[index] + (dt * k1u[index]) / 2.0;
        temp_v[index] = v[index] + (dt * k1v[index]) / 2.0;
        temp_h[index] = h[index] + (dt * k1h[index]) / 2.0;
      }
    }

    // k2 = f(yn + k1*dt/2)
    Partial_Derivative(temp_u, dudx, dudy, choice);

    Partial_Derivative(temp_v, dvdx, dvdy, choice);

    Partial_Derivative(temp_h, dhdx, dhdy, choice);

    // k2u = fu(temp_u, temp_v, dudx, dudy, dhdx);
    // k2v = fv(temp_u, temp_v, dvdx, dvdy, dhdy);
    // k2h = fh(temp_u, temp_v, temp_h, dudx, dvdy, dhdx, dhdy);

#pragma omp parallel for collapse(2)
    for (int i = 0; i < Nx; ++i) {
      for (int j = 0; j < Ny; ++j) {
        int index = i * Ny + j;
        k2u[index] = -u[index] * dudx[index] - v[index] * dudy[index] - g * dhdx[index];
        k2v[index] = -u[index] * dvdx[index] - v[index] * dvdy[index] - g * dhdy[index];
        k2h[index] = -u[index] * dhdx[index] - h[index] * dudx[index] - v[index] * dhdy[index] - h[index] * dvdy[index];

        temp_u[index] = u[index] + (dt * k2u[index]) / 2.0;
        temp_v[index] = v[index] + (dt * k2v[index]) / 2.0;
        temp_h[index] = h[index] + (dt * k2h[index]) / 2.0;
      }
    }

    // k3 = f(yn + k2*dt/2)
    Partial_Derivative(temp_u, dudx, dudy, choice);

    Partial_Derivative(temp_v, dvdx, dvdy, choice);

    Partial_Derivative(temp_h, dhdx, dhdy, choice);

    // k3u = fu(temp_u, temp_v, dudx, dudy, dhdx);
    // k3v = fv(temp_u, temp_v, dvdx, dvdy, dhdy);
    // k3h = fh(temp_u, temp_v, temp_h, dudx, dvdy, dhdx, dhdy);

#pragma omp parallel for collapse(2)
    for (int i = 0; i < Nx; ++i) {
      for (int j = 0; j < Ny; ++j) {
        int index = i * Ny + j;
        k3u[index] = -u[index] * dudx[index] - v[index] * dudy[index] - g * dhdx[index];
        k3v[index] = -u[index] * dvdx[index] - v[index] * dvdy[index] - g * dhdy[index];
        k3h[index] = -u[index] * dhdx[index] - h[index] * dudx[index] - v[index] * dhdy[index] - h[index] * dvdy[index];

        temp_u[index] = u[index] + (dt * k3u[index]);
        temp_v[index] = v[index] + (dt * k3v[index]);
        temp_h[index] = h[index] + (dt * k3h[index]);
      }
    }

    // k4 = f(yn + k3*dt)
    Partial_Derivative(temp_u, dudx, dudy, choice);

    Partial_Derivative(temp_v, dvdx, dvdy, choice);

    Partial_Derivative(temp_h, dhdx, dhdy, choice);

    // k4u = fu(temp_u, temp_v, dudx, dudy, dhdx);
    // k4v = fv(temp_u, temp_v, dvdx, dvdy, dhdy);
    // k4h = fh(temp_u, temp_v, temp_h, dudx, dvdy, dhdx, dhdy);

#pragma omp parallel for collapse(2)
    for (int i = 0; i < Nx; ++i) {
      for (int j = 0; j < Ny; ++j) {
        int index = i * Ny + j;
        k4u[index] = -u[index] * dudx[index] - v[index] * dudy[index] - g * dhdx[index];
        k4v[index] = -u[index] * dvdx[index] - v[index] * dvdy[index] - g * dhdy[index];
        k4h[index] = -u[index] * dhdx[index] - h[index] * dudx[index] - v[index] * dhdy[index] - h[index] * dvdy[index];

        u_next[index] = u[index] + (1.0 / 6.0) * (k1u[index] + 2.0 * k2u[index] + 2.0 * k3u[index] + k4u[index]) * dt;
        v_next[index] = v[index] + (1.0 / 6.0) * (k1v[index] + 2.0 * k2v[index] + 2.0 * k3v[index] + k4v[index]) * dt;
        h_next[index] = h[index] + (1.0 / 6.0) * (k1h[index] + 2.0 * k2h[index] + 2.0 * k3h[index] + k4h[index]) * dt;
      }
    }

    // end of pragma parallisation region
    // implicit parallel barrier for all the process to reach to this point and
    // swap to new conditions for next step
    // }

    // swap u_next and u to update to next time step
    cblas_dcopy(Nx * Ny, u_next, 1, u, 1);
    cblas_dcopy(Nx * Ny, v_next, 1, v, 1);
    cblas_dcopy(Nx * Ny, h_next, 1, h, 1);

    // Providing feedback to the user of how fast the program is running during
    // execution.
    Progress_Bar(t, time_step);

    // if (t % 80 == 0) {
    //   std::cout << "Time-step: " << std::setw(6) << t << "/" << time_step
    //             << " (" << std::setw(2) << 100 * t / time_step << "%)"
    //             << std::endl;
    // }
  }
  // release memory for defined arrays

  delete[] k1u;
  delete[] k2u;
  delete[] k3u;
  delete[] k4u;
  delete[] k1v;
  delete[] k2v;
  delete[] k3v;
  delete[] k4v;
  delete[] k1h;
  delete[] k2h;
  delete[] k3h;
  delete[] k4h;

  delete[] temp_u;
  delete[] temp_v;
  delete[] temp_h;

  std::cout << "\n\nFinished solving PDE (from t_i = 0 to t_f = T).\n"
            << std::endl;
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time Spent (sec) = "
            << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) / 1000000.0
            << std::endl;
};



void ShallowWater::SaveToFile() {

  std::cout << "\nWriting output of simulation to file 'output.txt'."
            << std::endl;

  std::ofstream vOut("output.txt", std::ios::out | std::ios::trunc);

  // Checking that file opened successfully.
  if (vOut.is_open()) {
    // Writing solution row-by-row (x y u v)
    // of storing matrices column-wise.
    for (int j = 0; j < Ny; ++j) {
      for (int i = 0; i < Nx; ++i) {
        vOut << i << " " << j << " " << u[i * Ny + j] << " " << v[i * Ny + j]
             << " " << h[i * Ny + j] << std::endl;
      }
      vOut << std::endl; // Empty line after each row of points
    }
  } else {
    std::cout << "Did not open vOut successfully!" << std::endl;
  }

  vOut.close(); // Closing file.
  std::cout << "Finished writing to file.\n" << std::endl;
};

/**
 * @brief Performs clean up duties.
 */
ShallowWater::~ShallowWater() {

  // De-allocating memory.
  delete[] u;
  delete[] v;
  delete[] h;
  delete[] u_next;
  delete[] v_next;
  delete[] h_next;

  delete[] dudx;
  delete[] dudy;
  delete[] dvdx;
  delete[] dvdy;
  delete[] dhdx;
  delete[] dhdy;
}





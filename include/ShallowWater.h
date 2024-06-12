// #ifndef SHA_WATER
// #define SHA_WATER

/**
 * @class encapsulates the solution fields u, v, and h. With built in 
 *        functions SetInitialConditions and TimeIntegrates to solve the 2D 
 *        shallow-water equations
 * 
 */
class ShallowWater {
    
    private: 

        /// define input variables and grid size
        double dt;
        int T;
        int Nx;
        int Ny;
        int ic;
        int choice;

        double g;
        int dx;
        int dy;
        int time_step;


        /// define solution fields
        double* u;
        double* v;
        double* h;

        // pre-define arrays for integration timesteps
        double* u_next;
        double* v_next;
        double* h_next;
        double* k1;
        double* k2;
        double* k3;
        double* k4;
        
        // pre-define arrays for partial derivatives
        double* dudx;
        double* dudy;
        double* dvdx;
        double* dvdy;
        double* dhdx;
        double* dhdy;


    public: 

        /// Takes parsed arguments, stores them and calculates various variables for performance improvements.
        void SetParameters(const double& arg_dt, const int& arg_T,
                        const int& arg_Nx, const int& arg_Ny,
                        const int& arg_ic, const int& arg_sol);
                            
        /// Set the Initial conditions defined from the question
        void SetInitialConditions();

        // Calculate the partial derivative terms with respect to x
        void Partial_Derivative(double* u, double* dudx, double* dudy, int choice);

        // Runge Kutta Integration
        // auto f(double* u, double* v, double* h);
        double* fu(double* u, double* v, double* dudx, double* dudy, double* dhdx);
        double* fv(double* u, double* v, double* dvdx, double* dvdy, double* dhdy);
        double* fh(double* u, double* v, double* h, double* dudx, double* dvdy, double* dhdx, double* dhdy);

        /// solve the equation by integrating through time till T
        void TimeIntegrate();

        /// Saves the result of the simulation to a .txt file.
        void SaveToFile();

        /// Deconstructor de-allocates all dynamically allocated memory.
        ~ShallowWater();   

}; 


// #endif

#include "utils.h"
#include "quantum_state.h"


int main()
{

    quantum_state state{-1, 0, 10}; // kappa (relativistic angular q. num.), nr (radial q. num.), Z (charge num.)
    // pot1_ver1_calc(state, 1e-9Q, 1e-8Q);
    // pot1_ver2_calc(state, 1e-9Q, 1e-8Q);
    pot1_calc(state, 1e-7Q, 1e-6Q);
    return 0;
}

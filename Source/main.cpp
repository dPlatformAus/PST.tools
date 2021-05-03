#include <algorithm>
#include <vector>
#include <deque>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <cstdlib>
#include <iostream> //gotta go through and clean this all up .. I think there may be some reduncancy here and also try to use consistent syntax!!
#include <fstream>
#include <sstream>
#if (defined(_WIN32) || defined(_WIN64) || defined(WIN32) || defined(WIN64)) && !defined(UNIX)
    #include <io.h>
#endif
#include <dirent.h>
#define NAMLEN(dirent) strlen((dirent)->d_name)

#include <boost/config.hpp>
#include <boost/timer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>


// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand base_generator_type;

//uncomment the following line if you want to compile for MPI
//#define USING_MPI 1
#if (defined(USING_MPI))
    #include <boost/filesystem.hpp>
    #include <boost/mpi.hpp>
    namespace mpi = boost::mpi;
    #include "MPI_tag_id_constants.h"
#endif

#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include "multi_array_serialize.hpp"

#if (defined(_WIN32) || defined(_WIN64) || defined(WIN32) || defined(WIN64)) && !defined(UNIX)
    #include <tnt.h>
    #include <jama_eig.h>
#else
    #include <tnt/tnt.h>
    #include <jama/jama_eig.h>
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#include "EasyBMP_1.06/EasyBMP.h"

#define DO_CLEAR_SCREEN 1
#define RENAME_FILES_ON_COMPLETION 0

const std::string input_suffix="in";
const std::string done_suffix="done";
const std::string input_directory="input/";
const std::string output_directory="output/";
const std::string scratch_directory="/scratch/mpiuser/PST/"; //this must be set up correctly if you compile with mpi enabled and try to use calculations that write to local scratch (like 3F calc version that writes all states to file)

std::ostringstream persistent_output, semipersistent_output;
unsigned job_number;

#include "misc.h"
#include "colour_scale.h"

#include "factorials.h"
#include "range.h"
#include "histogram.h"
#include "histogram_2d.h"
#include "histogram_3d.h"
#include "progress_bar.h"
Progress_Bar the_progress_bar; //we define the progress bar as global so that the highest level progress bar can take precedence ... not sure how to do this with MPI progress bar yet ... but I don't think it will be multi level, so I have not defined it here yet
#if (defined(USING_MPI))
    #include "MPI_progress_bar.h"
#endif
#include "atom.h"
#include "geom_molecule.h"
#include "fit_point.h"
#include "phasespace_result_data.h"
#include "phasespace_probe.h"
#include "parse.h"
#include "base_calculation.h"
#include "rebin_calculation.h"
#include "phasespace_molecule.h"
#include "phasespace_result.h"
#include "phasespace_input.h"
#include "jvE_histogram.h"
#include "phasespace_core.h"
#include "phasespace_calculation.h"
#include "phasespace_get_state_list_calculations.h"
#include "impulsive_molecule.h"
#include "impulsive_input.h"
#include "ho_wavefunction.h"
#include "impulsive_core.h"
#include "impulsive_calculation.h"
#include "triple_fragmentation_input.h"
#include "triple_fragmentation_core.h"
#include "triple_fragmentation_calculation.h"
#include "roaming_result_data.h"
#include "roaming_probe.h"
#include "roaming_molecule.h"
#include "roaming_result.h"
#include "roaming_input.h"
#include "roaming_core.h"
#include "roaming_calculation.h"
#include "roaming_fraction_energy_dependance_calculation.h"
#include "roaming_fit_delta_E_roam_and_P_roam_calculation.h"
#include "roaming_fit_trot_and_tvib_and_p_roam_calculation.h"
/*#include "tunnel_tests.h"*/
#include "standard_execution.h"
#if (defined(USING_MPI))
    #include "slave_execution.h"
#endif

int main()
{
    initialise_log_factorial();
#if (defined(USING_MPI))
    mpi::environment env;
    mpi::communicator world;
    if (world.rank()==0){
#endif
        do_standard_execution();
#if (defined(USING_MPI))
    }else{
        do_slave_execution();
    }
#endif
    return 0;
}

#define UNITS_GAUSSIAN 1
#define UNITS_SI 0

struct Wavefunction_Coord {
    double x, deg, rad;
};
struct Wavefunction_Value {
    double psi, prob;
};

class HO_Wavefunctions //this class makes probability (psi^2) functions for the bond angle of the breaking bond used for the impulse calculation
{
    public:
        unsigned n_max;
        double frequency, reduced_mass, arc_radius, equilibrium_angle, force_constant, alpha, xi_0, x_0, dx;
        std::vector<unsigned> trim_index;
        std::vector<Wavefunction_Coord> coord;
        std::vector< std::vector<Wavefunction_Value> > wavefunctions;
        Output_File out_file;

        ~HO_Wavefunctions();
        HO_Wavefunctions(const std::string& out_filename, const std::string &job_path, const Impulsive_General_Input *the_general_input, const unsigned units, const unsigned in_n_max, const Impulsive_Molecule *frag);
        void make_units_SI();
        double Hermite(const unsigned n, const double xi);
        void normalise_wavefunctions();
        void write_wavefunctions(FILE *the_file);
        unsigned start(const unsigned n);
        unsigned end(const unsigned n);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HO_Wavefunctions::~HO_Wavefunctions() {
	out_file.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HO_Wavefunctions::HO_Wavefunctions(const std::string& out_filename, const std::string &job_path, const Impulsive_General_Input *the_general_input, const unsigned units, const unsigned in_n_max, const Impulsive_Molecule *frag) {
    bool trim_index_found;
    unsigned n, i;
    double i_deg, i_rad, xi;
    Wavefunction_Coord this_coord;
    Wavefunction_Value this_val;
    std::vector<Wavefunction_Value> this_wavefunction;
    n_max=in_n_max;
    frequency=frag->TS_bend_frequency;
    reduced_mass=frag->TS_bend_reduced_mass;
    arc_radius=frag->TS_bend_arc_radius;
    equilibrium_angle=frag->TS_frag1_equilibrium_angle;
    if (units==UNITS_GAUSSIAN) make_units_SI(); //otherwise assume inputs already in SI units
    alpha=sqrt(reduced_mass*frequency/HBAR);
    for (i_deg=-180; i_deg<180; i_deg+=the_general_input->wavefunction_step_size) {
        i_rad=i_deg*RAD_IN_DEGREE;
        this_coord.x=i_rad*arc_radius;
        this_coord.rad=equilibrium_angle+i_rad;
        if (this_coord.rad<0) this_coord.rad+=2*PI;
        this_coord.deg=this_coord.rad/RAD_IN_DEGREE;
        coord.push_back(this_coord);
    }
    for (n=0; n<=n_max; n++){
        xi_0=sqrt(2*n+1);
        x_0=sqrt((2*n+1)*HBAR/(reduced_mass*frequency)); //the classical turning point, where E=V
        this_wavefunction.clear();
        trim_index_found=0;
        for (i=0; i<coord.size(); i++) {
            xi=alpha*coord[i].x;
            this_val.psi = sqrt( 1 / (sqrt(PI)*pow((double)2,(double)n)*factorial(n)) ) * Hermite(n,xi) * exp(-pow(xi,2)/2);
            this_val.prob=pow(this_val.psi,2);
            if (this_val.prob>=the_general_input->wavefunction_zero_threshold && !trim_index_found) {
                trim_index.push_back(i);
                trim_index_found=1;
            }
            this_wavefunction.push_back(this_val);
        }
        wavefunctions.push_back(this_wavefunction);
    }
    normalise_wavefunctions();
    if (the_general_input->HO_Wavefunctions_write_to_file) {
        out_file.set_output_directory(job_path+"HO_Wavefunctions/");
        out_file.open(out_filename);
        write_wavefunctions(out_file.fp);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double HO_Wavefunctions::Hermite(const unsigned n, const double xi) {
    std::vector<double> h;
    unsigned i;
    h.push_back(1);
    h.push_back(2*xi);
    for (i=1; i<n; i++) h.push_back( (2*xi*h[i]) - (2*i*h[i-1]) );
    return h[n];
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned HO_Wavefunctions::start(const unsigned n){
    return trim_index[n];
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned HO_Wavefunctions::end(const unsigned n){
    return coord.size()-trim_index[n];
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void HO_Wavefunctions::normalise_wavefunctions(){
    unsigned i, n;
    double sum;
    for (n=0; n<wavefunctions.size(); n++) {
        sum=0;
        for (i=start(n); i<end(n); i++) sum+=wavefunctions[n][i].prob;
        for (i=start(n); i<end(n); i++) wavefunctions[n][i].prob/=sum;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void HO_Wavefunctions::make_units_SI(){
    frequency*=RAD_PER_SEC_IN_WAVENUMBER;
    reduced_mass*=KG_IN_AMU;
    arc_radius*=METERS_IN_ANGSTROM;
    //equilibrium_angle*=RAD_IN_DEGREE;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void HO_Wavefunctions::write_wavefunctions(FILE *the_file){
    unsigned i, n;
    for (n=0; n<wavefunctions.size(); n++) {
        fprintf(the_file, "\ndisp,rad,deg,%d psi,%d prob\n",n,n);
        for (i=start(n); i<end(n); i++) fprintf(the_file, "%e,%e,%e,%e,%e\n",coord[i].x, coord[i].rad, coord[i].deg, wavefunctions[n][i].psi, wavefunctions[n][i].prob);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


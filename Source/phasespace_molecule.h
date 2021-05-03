#define COMBINATION_BANDS 0
#define NORMAL_MODES 1
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Spin_Orbit_State
{
    public:
        double energy;
        float j; ///not integer, electron spin is +- 0.5
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & energy & j;}
};
bool spin_orbit_state_sort_predicate(const Spin_Orbit_State &lhs, const Spin_Orbit_State &rhs) {return lhs.energy < rhs.energy;}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Vibrational_Level
{
    public:
        double energy;
        unsigned degeneracy;
        std::vector <unsigned> mode_tally; ///how many of each mode this combo band has (if it's a combo band)
        std::string name; /// like 2^1 3^4 etc or just the mode number if a normal mode
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & energy & degeneracy & mode_tally & name;}
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool vibrational_level_sort_predicate(const Vibrational_Level &lhs, const Vibrational_Level &rhs) {return lhs.energy < rhs.energy;}
//***********************************************************************************************************************************************************************************
class Combination_Band_List
{
    public:
        std::vector<Vibrational_Level> modes;    ///contains normal modes
        std::vector<Vibrational_Level> bands;    ///contains combination bands (which may or may not have been generated from norml modes using the harmonic oscillator approx)

        Vibrational_Level& operator[] (unsigned i) {return bands[i];}
        const Vibrational_Level& operator[] (unsigned i) const {return bands[i];}
        unsigned size() const {return bands.size();}
        void sort_bands() {std::sort(bands.begin(), bands.end(), vibrational_level_sort_predicate); }
        void sort_modes() {std::sort(modes.begin(), modes.end(), vibrational_level_sort_predicate); }
        void analyse_degeneracy();
        void initialise_blank(const bool type_energy_levels, const unsigned num_vib);
        void get_combination_bands(const double E_available);
        void enumerate(const double E_available, std::vector<unsigned> &combo_band, unsigned length_remaining, unsigned in_i);
        void make_name(std::vector<unsigned> &combo_band, Vibrational_Level &combo_result);
        void add_zero_band();
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & modes & bands;}
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Combination_Band_List::add_zero_band() {
    unsigned i;
    Vibrational_Level temp_level;
    temp_level.degeneracy=1;
    temp_level.energy=0;
    temp_level.name="0";
    for (i=0; i<modes.size(); i++) temp_level.mode_tally.push_back(0);
    bands.push_back(temp_level);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Combination_Band_List::analyse_degeneracy() { ///this compresses the list of modes by counting the number of modes with the same energy and storing each energy once with the appropriate degeneracy instead ... will speed up PST count significantly
    unsigned i, degen_i;
    double last_band_energy=-1;
    std::vector<Vibrational_Level> degenerate_bands;
    Vibrational_Level temp_level;
    temp_level.degeneracy=1;
    temp_level.energy=0;
    degen_i=0;
    for (i=0; i < bands.size(); i++) { ///this assumes that the bands list has already been sorted!
        if (bands[i].energy == last_band_energy) {
            degenerate_bands[degen_i].degeneracy+=bands[i].degeneracy; ///at the moment we have forced bands[i].degeneracy to be 1 if generated from modes, however if we have input bands directly a degeneracy other than 1 may have been defined
            //cout<<"incrementing";
        } else {
            temp_level.degeneracy=1;
            temp_level.energy=last_band_energy=bands[i].energy;
            temp_level.name=bands[i].name;
            temp_level.mode_tally=bands[i].mode_tally;
            degenerate_bands.push_back(temp_level);
            degen_i=degenerate_bands.size()-1;
        }
    }
    bands=degenerate_bands;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Combination_Band_List::initialise_blank(const bool type_energy_levels, const unsigned num_vib){
    unsigned i;
    Vibrational_Level temp_level;
    temp_level.degeneracy=1;
    temp_level.energy=0;
    modes.clear();
    bands.clear();
    if (type_energy_levels==NORMAL_MODES) for (i=0;i<num_vib;i++) modes.push_back(temp_level);
    else for (i=0;i<num_vib;i++) bands.push_back(temp_level);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Combination_Band_List::get_combination_bands(const double E_available) {
    unsigned longest_combo=0;
    std::vector <unsigned>  combo_band; ///a temporary Combination_Band that will be modified inside and then stored
    if (modes.size()){ //only do this if there are some modes in here
        sort_modes(); //ensure input sorted ascending so that first element is the smallest (see next statement)
        if (modes[0].energy) longest_combo=(unsigned)floor(E_available/modes[0].energy); //now many times does smallest vib fit into energy available (assumes first is smallest)
        enumerate(E_available, combo_band, longest_combo, 0);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Combination_Band_List::enumerate(const double E_available, std::vector<unsigned> &combo_band, unsigned length_remaining, unsigned in_i) {
    unsigned i, j;
    Vibrational_Level combo_result;
    length_remaining-=1; //length remaining is one less after adding this character
    for (i=in_i; i < modes.size(); i++) { //try adding all the characters in the selection
        combo_band.push_back(i); // add this energy level to the combination (*** store the index ***)
        combo_result.degeneracy=1;
        combo_result.energy=0;
        for (j=0; j<combo_band.size(); j++){ //for each energy level in this combination, calculate its contribution to the combination band
            combo_result.degeneracy*=modes[ combo_band[j] ].degeneracy;
            combo_result.energy+=modes[ combo_band[j] ].energy;
        }
        // can check this combination is within E_available - kill recursion when E_available is exceeded preventing too many unnessesary calculations
        if (combo_result.energy<=E_available) { //this combo still fits into energy
            make_name(combo_band, combo_result);
            bands.push_back(combo_result); // store it
            if (length_remaining) {//length_remaining > 0 so do next length too - kill recursion once max length met
                enumerate(E_available, combo_band, length_remaining, i); //enumerate all starting with this combo, start from this normal mode so we get all combos starting with it repeating
            }
        }
        combo_band.pop_back(); // remove this energy level from the combination
    }//move to next energy level
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Combination_Band_List::make_name(std::vector<unsigned> &combo_band, Vibrational_Level &combo_result){
    unsigned i;
    std::stringstream result;
    bool write_separator=0;
    combo_result.mode_tally.clear();
    for (i=0; i<modes.size(); i++) combo_result.mode_tally.push_back(0);
    for (i=0; i<combo_band.size(); i++) combo_result.mode_tally[ combo_band[i] ]++;
    for (i=0; i<combo_result.mode_tally.size(); i++) {
        if (combo_result.mode_tally[i]) {
            if (write_separator) result<<":";
            result<<i<<"^"<<combo_result.mode_tally[i];
            write_separator=1;
        }
    }
    combo_result.name=result.str();
}
//***********************************************************************************************************************************************************************************
class Phasespace_Molecule/**
This is the base molecule class. There will be several derivatives of this class that build on it's functionality for specific subtypes of Phasespase Theory (PST) calculations.
This class stores various data relating to a molecule or fragment.
It also has a number of functions that calculate values related to the fragment, which are used by the core PST state count class (and elsewhere, such as in the run functions of some calculation class derivatives).
This class knows how to parse an input file and populate all of the data in the files relating to itself.
This class also has some functions that write output specific to itself to a file.
*/{
    public:
        unsigned number;
        std::string name;
        double E_dissociation; //cm-1 //if this is non zero a triple frag calculation will occur on this fragment
        double mass;  // of this fragment (AMU)
        double A, B, C;  // rotational constants set_Bbar_and_prolate()
        double Bbar; //set in set_Bbar_and_prolate() before use externally
        double C_6;
        double b_max; // A (maximum impact parameter)
        bool single_atom, linear, open_shell, prolate, rot_from_geom, vib_levels_initialised; // 1=prolate 0=oblate, 1=open shell 0= closed shell, 1= linear 0= nonlinear etc
        unsigned initial_J; //the angular momentum of this molecule, used for parent/reactant molecules (this is what used to be J_parent)
        double vel_mass_factor;
        Combination_Band_List vib_levels;
        std::vector<Spin_Orbit_State> spin_orbit_states;    //contains spin orbit state energies, angular momenta
        std::vector<Probe_State> probe_states;
        Histogram_Use *is_on, *probe_is_on;
 //       Histogram vel_histogram_rad, E_trans_histogram_rad, E_int_histogram_rad, E_rot_histogram_rad, E_vib_histogram_rad, v_histogram_rad;
 //       Histogram_2d E_int_trans_histogram_rad, nk_histogram_rad;

        Phasespace_Molecule(){};
        virtual Phasespace_Molecule* clone(){return new Phasespace_Molecule(*this);};
        virtual ~Phasespace_Molecule(){delete is_on; delete probe_is_on;};
        virtual void construct();
        virtual void create_histogram_use();
        virtual void parse_input(const std::string &file_name, const std::string &section_name);
        void add_mode(const double energy, const unsigned degeneracy);
        void add_band(const double energy, const unsigned degeneracy);
        void add_spin_orbit_state(const double energy, const float j);
        void set_Bbar_and_prolate();
        unsigned get_n_max(const double E_available) const;
        double get_E_nk(const unsigned N, const unsigned K) const;
        double get_vel_mass_factor(const double other_frag_mass);
        void sort_spin_orbit_states() {std::sort(spin_orbit_states.begin(), spin_orbit_states.end(), spin_orbit_state_sort_predicate); }
        void initialise(const double other_frag_mass);
        void initialise_vib_levels(const double E_available);
        void write_fragment_state(FILE *the_file);

        unsigned get_v_max(const double E_max) const;
        unsigned get_so_max(const double E_max) const;
        unsigned get_degeneracy(const unsigned K, const unsigned V, const unsigned SO) const;
        double get_velocity(const double E_translation) const;
        double get_E_trans(const double velocity) const;
    //private:
        Phasespace_Molecule(const Phasespace_Molecule &rhs);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & number & E_dissociation & mass & A & B & C & Bbar & C_6 & b_max & single_atom & linear & open_shell & prolate & rot_from_geom & vib_levels_initialised & initial_J & vel_mass_factor;
            ar & vib_levels & spin_orbit_states & probe_states & is_on & probe_is_on;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Molecule::Phasespace_Molecule(const Phasespace_Molecule &rhs){
    number=rhs.number;
    E_dissociation=rhs.E_dissociation;
    mass=rhs.mass;
    A=rhs.A;
    B=rhs.B;
    C=rhs.C;
    Bbar=rhs.Bbar;
    C_6=rhs.C_6;
    b_max=rhs.b_max;
    single_atom=rhs.single_atom;
    linear=rhs.linear;
    open_shell=rhs.open_shell;
    prolate=rhs.prolate;
    rot_from_geom=rhs.rot_from_geom;
    vib_levels_initialised=rhs.vib_levels_initialised;
    initial_J=rhs.initial_J;
    vel_mass_factor=rhs.vel_mass_factor;
    vib_levels=rhs.vib_levels;
    spin_orbit_states=rhs.spin_orbit_states;
    probe_states=rhs.probe_states;
    is_on=rhs.is_on->clone();
    probe_is_on=rhs.probe_is_on->clone();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Molecule::construct(){
    name="untitled";
    E_dissociation=0;
    mass=0;
    initial_J=0;
    A=B=C=Bbar=0;
    C_6=b_max=0;
    single_atom=linear=0; //default to nonlinear polyatomic fragments for backwards compatibility
    open_shell=1; //default to radical fragments for backwards compatibility
    prolate=0;
    rot_from_geom=0;///legacy flag for legacy rrkm geom type calculation ... new flags stored in general_input
    vel_mass_factor=0;
    vib_levels_initialised=0;
    create_histogram_use();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Molecule::create_histogram_use(){
    is_on=new Histogram_Use();
    probe_is_on=new Histogram_Use();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Molecule::parse_input(const std::string &file_name, const std::string &section_name){
    const unsigned NO_SECTION=0;
    const unsigned GENERAL=1;
    const unsigned FREQUENCIES=2;
    const unsigned SPIN_ORBIT_ENERGIES=3;
    const unsigned ATOMS=4;
    const unsigned PROBE_STATES=5;
    const unsigned HISTOGRAM_USE=6;
    const unsigned PROBE_HISTOGRAM_USE=7;
	unsigned read_state=NO_SECTION;
	std::string file_name_and_path=input_directory+file_name;
    std::ifstream infile (file_name_and_path.c_str());
	std::string line, key, value, symbol;
    unsigned degeneracy;
    double frequency, energy, ax, ay, az;
    float j;
    Probe_State temp_probe_state;
    unsigned calculate_simple_properties(0); //flag to indicate that a geom_molecule object will be used to obtain information about this molecule from its atomic coordinates etc
    Geom_Molecule *mol_geom; //will be created and used to derive some data from imput coordinates and then deleted
    semipersistent_output << "Phasespace_Molecule "<< number <<" parsing" << std::endl;
    print_semipersistent_output();
	if (infile.is_open()) {
        if (calculate_simple_properties) mol_geom=new Geom_Molecule(section_name);
        while(getline(infile,line)) {
            if (line.substr(0,2) != "//"){ //if this is not a comment line (we should ignore comment lines)
                if (line.substr(0,8) == "--------"){ //this is the start of a new section
                    read_state=NO_SECTION; //set to none first in case the section is not recognised by this parse routine
                    if (line.find( std::string("-")+section_name+std::string("-") ) != std::string::npos) read_state=GENERAL;
                    else if (line.find( std::string("-")+section_name+std::string(" FREQUENCIES-") ) != std::string::npos) read_state=FREQUENCIES;
                    else if (line.find( std::string("-")+section_name+std::string(" SPIN ORBIT ENERGIES-") ) != std::string::npos) read_state=SPIN_ORBIT_ENERGIES;
                    else if (line.find( std::string("-")+section_name+std::string(" ATOMS-") ) != std::string::npos) read_state=ATOMS;
                    else if (line.find( std::string("-")+section_name+std::string(" PROBE STATES-") ) != std::string::npos) read_state=PROBE_STATES;
                    else if (line.find( std::string("-")+section_name+std::string(" HISTOGRAM USE-") ) != std::string::npos) read_state=HISTOGRAM_USE;
                    else if (line.find( std::string("-")+section_name+std::string(" PROBE HISTOGRAM USE-") ) != std::string::npos) read_state=PROBE_HISTOGRAM_USE;
                }
                if (read_state==GENERAL){
                    read_key_value(line, key, value);
                    if (key != "") {
                        if (key == "fragment_name") name=value;
                        else if (key == "dissociation_energy") E_dissociation=atof(value.c_str());
                        else if (key == "A") A=atof(value.c_str());
                        else if (key == "B") B=atof(value.c_str());
                        else if (key == "C") C=atof(value.c_str());
                        else if (key == "C_6") C_6=atof(value.c_str());
                        else if (key == "b_max") b_max=atof(value.c_str());
                        else if (key == "mass") mass=atof(value.c_str());
                        else if (key == "single_atom") single_atom=atoi(value.c_str());
                        else if (key == "linear") linear=atoi(value.c_str());
                        else if (key == "open_shell") open_shell=atoi(value.c_str());
                        else if (key == "initial_J") initial_J=atoi(value.c_str());
                        else if (key == "calculate_simple_properties"){
                            calculate_simple_properties=atoi(value.c_str());
                        }
                    }
                }else if (read_state==HISTOGRAM_USE){
                    read_key_value(line, key, value);
                    if (key != "") is_on->read_histogram_use(key, value);
                }else if (read_state==PROBE_HISTOGRAM_USE){
                    read_key_value(line, key, value);
                    if (key != "") probe_is_on->read_histogram_use(key, value);
                }else if (read_state==FREQUENCIES){
                    if (read_mode_value(line, frequency, degeneracy)) add_mode(frequency, degeneracy);
                }else if (read_state==SPIN_ORBIT_ENERGIES){ ///note: if values are entered in input, the 0 energy level must be in the input file (otherwise it will not be considered)
                    if (read_spin_orbit_value(line, energy, j)) add_spin_orbit_state(energy, j);
                }else if (read_state==ATOMS){ //in the input file the atoms section must come after the general section to ensure the Geom_Molecule has been created!
                    if (calculate_simple_properties && read_atom_value(line, symbol, ax, ay, az)) mol_geom->add_atom(symbol, ax, ay, az);
                }else if (read_state==PROBE_STATES){
                    if (temp_probe_state.readState(line)) probe_states.push_back(temp_probe_state);
                }
            }
        }
        infile.close();
        if (calculate_simple_properties){
            mol_geom->get_simple_properties(mass, A, B, C); //calculates fragment mass and rotational constants and sets the values passed by reference into the function
            delete mol_geom;
        }
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Molecule::initialise(const double other_frag_mass){
    set_Bbar_and_prolate();
    Spin_Orbit_State temp_state; //this is not copy constructor stuff ... this needs to be in a seperate function, initialise_spin_orbit_states() or something!
    if (!spin_orbit_states.size()) { //no spin orbit states defined in input assign default values ... this code assumes 0 energy is defined explicitly in input if any spin orbit energies HAVE been defined in the input file
         if (open_shell) { //this is the default, a radical with a single unpaired electron (set in Phasespace_Calc_Input constructor, can be overridden by explicitly setting it in input file) *** probably better to do this the other way round and call it closed_shell!
            temp_state.energy=0;
            temp_state.j=0.5;
            spin_orbit_states.push_back(temp_state); //add the spin +1/2 state (the -1/2 state comes out in the degeneracy calculation 2j+1)
         } else { //A closed shell molecule is a singlet so just store a single spin orbit state with no angular momentum and no degeneracy
            temp_state.energy=0;
            temp_state.j=0;
            spin_orbit_states.push_back(temp_state); //add the 0 energy 0 j singly degenerate spin orbit state
         }
    }
    sort_spin_orbit_states(); //sort spin_orbit_states ascending

    vel_mass_factor=get_vel_mass_factor(other_frag_mass);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Molecule::initialise_vib_levels(const double E_available) {
    if (!vib_levels_initialised){ //don't do this twice! (and make sure, if you do this to a molecule before it's passed to PST core/s, that you do it for the maximum possible E_available, otherwise you'll miss levels above that energy)
        vib_levels_initialised=1;
        if (!vib_levels.size()) vib_levels.get_combination_bands(E_available); //if there are no combination bands, assume normal modes are populated and get the combination bands
        vib_levels.add_zero_band();//add the 0 vib energy entry (makes the code for getting available energy neater)
        vib_levels.sort_bands(); //sort energy_levels ascending
        vib_levels.analyse_degeneracy(); //compress list by finding degenerate bands and storing them once with the appropriate degeneracy (dramatically reduces the number of loops in count_states();)
    }else persistent_output<<"ERROR: Phasespace_Molecule::initialise_vib_levels -- DOUBLE CALL!"<<std::endl;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Molecule::add_mode(const double energy, const unsigned degeneracy) {
    unsigned i;
    Vibrational_Level temp_level;
    temp_level.degeneracy=1;
    for (i=1; i<=degeneracy; i++) { ///THIS IS IMPORTANT! a doubly degenerate mode has to be treated as 2 singly degenerate modes, 2 quanta of a doubly degenerate mode has a degeneracy of 3 not 4!
        temp_level.energy=energy;
        vib_levels.modes.push_back(temp_level);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Molecule::add_band(const double energy, const unsigned degeneracy) {
    Vibrational_Level temp_level;
    temp_level.degeneracy=degeneracy;
    temp_level.energy=energy;
    vib_levels.bands.push_back(temp_level);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Molecule::add_spin_orbit_state(const double energy, const float j) { ///note: if values are entered in input, the 0 energy level must be in the input file (otherwise it will not be considered)
    Spin_Orbit_State temp_state;
    temp_state.j=j;
    temp_state.energy=energy;
    spin_orbit_states.push_back(temp_state);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Molecule::set_Bbar_and_prolate() {
    if (linear) {
        Bbar=B;//there is only one rotational constant for a linear molecule
        prolate=1;//nothing is more cigar shaped than a linear molecule ... we should expicitly check the linear flag, so this should be redundant, but better to be robust
    } else {
        Bbar=(B+C)/2;
        prolate=(A>Bbar);// prolate is cigar shape, so Principal I is smallest => Principal B is biggest -- oblate is pancake shape, so Principal I is biggest => Principal B is smallest
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Phasespace_Molecule::get_n_max(const double E_available) const { //caution: don't call this unless set_Bbar_and_prolate() has been called
    //   Calculate max N, occurs when: K=0 if prolate (like HCO), K=N if oblate (like CH3), 0 if single atom
    if (single_atom) return 0;
    else { // ceil is important here, don't want to exclude anything ... if the resultant nk state is too big it is rejected anyway, but if we don't test and it fits - we are wrong
        if (prolate || linear) return (unsigned) ceil( (-Bbar+sqrt(pow(Bbar,2)+4*Bbar*E_available))/(2*Bbar) );   //Max is BN(N+1)  ie K=0 ... if linear==1 then prolate==1, so || is redundant, but it reads better
        else return (unsigned) ceil( (-Bbar+sqrt(pow(Bbar,2)+4*A*E_available))/(2*A) ); //Max is BN(N+1)+(A-B)K^2   where N=K
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Phasespace_Molecule::get_E_nk(const unsigned N, const unsigned K) const { ///caution: don't call this unless set_Bbar_and_prolate() has been called
    if (single_atom) return 0;
    else return Bbar*N*(N+1)+(A-Bbar)*pow((double)K,2); // Calculate the energy of the (n,k) state
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Phasespace_Molecule::get_vel_mass_factor(const double other_frag_mass) { ///should probably errorcheck divide by 0 in case someone puts in 0 mass fragment or something
     return ( 2 / (mass*KG_IN_AMU) ) * ( other_frag_mass / (mass+other_frag_mass) ) * JOULES_IN_WAVENUMBER;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Molecule::write_fragment_state(FILE *the_file) {
    unsigned i;
    if (single_atom) fprintf(the_file, "\nFragment %d is a single atom\n", number);
    else if (linear) fprintf(the_file, "\nFragment %d is linear (prolate)\n", number);
    else if (prolate) fprintf(the_file, "\nFragment %d is prolate\n", number);
    else fprintf(the_file, "\nFragment %d is oblate\n", number);
    if (!single_atom) {
        fprintf(the_file, "\nFragment %d rotational constants\n", number);
        if (!linear) fprintf(the_file, "\nA = %lf\n", A);
        fprintf(the_file, "\nB = %lf\n", B);
        if (!linear) fprintf(the_file, "\nC = %lf\n", C);
    }
    fprintf(the_file, "\n%d Sorted Available Vibrational Energy Levels\nenergy    degeneracy    name\n", vib_levels.size()-1);
    if (vib_levels.size() < 250) for (i=1;i<vib_levels.size(); i++) fprintf(the_file, "%lf   %d   %s\n", vib_levels[i].energy, vib_levels[i].degeneracy, vib_levels[i].name.c_str());
    fprintf(the_file, "\n%d Sorted Spin Orbit State Energy Levels\nenergy    j\n", (int)spin_orbit_states.size());
    for (i=0;i<spin_orbit_states.size(); i++) fprintf(the_file, "%lf   %lf\n", spin_orbit_states[i].energy, spin_orbit_states[i].j);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Phasespace_Molecule::get_v_max(const double E_max) const { ///assumes list sorted ascending
    unsigned i=vib_levels.size()-1; ///important to start from back so that all within the range are included in the count (which loops from front)
    unsigned result=0;
    bool found_v_max=0;
    if (!vib_levels_initialised){
        persistent_output<<"vib_levels_initialised FALSE!"<<std::endl;
        print_persistent_output();
    }
    if (E_max>0){
        while (i>=0 && !found_v_max){ //set the max vib level that fits inside the remaining energy
            if (E_max>vib_levels[i].energy){
                result=i;
                found_v_max=1;
            }
            i-=1;
        }
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Phasespace_Molecule::get_so_max(const double E_max) const { ///assumes list sorted ascending
    int i=spin_orbit_states.size()-1; ///important to start from back so that all within the range are included in the count (which loops from front)
    unsigned result=0;
    bool found_so_max=0;
    while (i>=0 && !found_so_max){ //set the max spin orbit state that fits inside the remaining energy
        if (E_max>spin_orbit_states[i].energy){
            result=(unsigned)i;
            found_so_max=1;
        }
        i-=1;
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Phasespace_Molecule::get_degeneracy(const unsigned K, const unsigned V, const unsigned SO) const{
    unsigned degeneracy=1;
    degeneracy*=int(2*spin_orbit_states[SO].j+1); //account for degeneracy of the spin orbit state of fragment in Hund's case (a) the degeneracy of all j states is 2j+1
    if (K) degeneracy*=2; //double the degeneracy for non zero k, as k can be plus or minus
    degeneracy*=vib_levels[V].degeneracy; //account for degeneracy of the vibrational level of fragment
    return degeneracy;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Phasespace_Molecule::get_velocity(const double E_translation) const {
    return sqrt(vel_mass_factor*E_translation);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Phasespace_Molecule::get_E_trans(const double velocity) const {
    return 0.5*mass*KG_IN_AMU*pow(velocity,2)/JOULES_IN_WAVENUMBER;
}

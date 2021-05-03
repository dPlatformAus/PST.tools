//***********************************************************************************************************************************************************************************
class General_Input
{
    public:
        Number_Range<double> excitation_energy_range; // cm-1
        double excitation_energy_increment; // cm-1
        Number_Range<double> equilibrium_temerature_range, rotational_temerature_range, vibrational_temerature_range; // K
        double equilibrium_temerature_increment, rotational_temerature_increment, vibrational_temerature_increment; //K default is 5
        double equilibrium_temerature, rotational_temerature, vibrational_temerature; //K
        bool jvE_histogram_write_to_file, histogram_write_to_file, objects_write_to_file;
        unsigned Boltzmann_bin_size; // default is 5
        double Boltzmann_max_energy_factor; //default is 5
        double Boltzmann_probability_included; //% -- default is 100
        double get_min_energy_within_top_percent; //% -- default is 0
        double get_percent_within_min_energy; //cm-1 -- default is 0
        bool allow_Tvib_less_than_Trot; // default is 0
        bool Boltzmann_j_averaging;

        General_Input();
        virtual ~General_Input(){};
        virtual void parse_input(const std::string &file_name);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
General_Input::General_Input(){
    excitation_energy_increment=10;
    //histogram_size=100;
    jvE_histogram_write_to_file=histogram_write_to_file=objects_write_to_file=0;
    equilibrium_temerature=rotational_temerature=vibrational_temerature=0;
    equilibrium_temerature_increment=rotational_temerature_increment=vibrational_temerature_increment=5;
    Boltzmann_bin_size=5;
    Boltzmann_max_energy_factor=5;
    Boltzmann_probability_included=100;
    get_min_energy_within_top_percent=0;
    get_percent_within_min_energy=0;
    allow_Tvib_less_than_Trot=0;
    Boltzmann_j_averaging=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void General_Input::parse_input(const std::string &file_name){
    const unsigned NO_SECTION=0;
    const unsigned GENERAL=1;
	unsigned read_state=NO_SECTION;
	std::string file_name_and_path=input_directory+file_name;
    std::ifstream infile (file_name_and_path.c_str());
	std::string line, key, value;
    semipersistent_output << "General_Input parsing" << std::endl;
    print_semipersistent_output();
	if (infile.is_open()) {
        while(getline(infile,line)) {
            if (line.substr(0,2) != "//"){ //if this is not a comment line (we should ignore comment lines)
                if (line.substr(0,8) == "--------"){ //this is the start of a new section
                    read_state=NO_SECTION; //set to none first in case the section is not recognised by this parse routine
                    if (line.find("-GENERAL-") != std::string::npos) read_state=GENERAL;
                }
                if (read_state==GENERAL){
                    read_key_value(line, key, value);
                    if (key != "") {
                        if (key == "excitation_energy_range") excitation_energy_range.setRange(value);
                        else if (key == "excitation_energy_increment") excitation_energy_increment=atof(value.c_str());
                        else if (key == "histogram_size") HISTOGRAM_SIZE=atoi(value.c_str());
                        else if (key == "histogram_2d_size") HISTOGRAM_2D_SIZE=atoi(value.c_str());
                        else if (key == "jvE_histogram_write_to_file") jvE_histogram_write_to_file=atoi(value.c_str());
                        else if (key == "histogram_write_to_file") histogram_write_to_file=atoi(value.c_str());
                        else if (key == "objects_write_to_file") objects_write_to_file=atoi(value.c_str());
                        else if (key == "equilibrium_temerature") equilibrium_temerature=atof(value.c_str());
                        else if (key == "rotational_temerature") rotational_temerature=atof(value.c_str());
                        else if (key == "vibrational_temerature") vibrational_temerature=atof(value.c_str());
                        else if (key == "equilibrium_temerature_range") equilibrium_temerature_range.setRange(value);
                        else if (key == "equilibrium_temerature_increment") equilibrium_temerature_increment=atof(value.c_str());
                        else if (key == "rotational_temerature_range") rotational_temerature_range.setRange(value);
                        else if (key == "rotational_temerature_increment") rotational_temerature_increment=atof(value.c_str());
                        else if (key == "vibrational_temerature_range") vibrational_temerature_range.setRange(value);
                        else if (key == "vibrational_temerature_increment") vibrational_temerature_increment=atof(value.c_str());
                        else if (key == "Boltzmann_bin_size") Boltzmann_bin_size=atoi(value.c_str());
                        else if (key == "Boltzmann_max_energy_factor") Boltzmann_max_energy_factor=atof(value.c_str());
                        else if (key == "Boltzmann_probability_included") Boltzmann_probability_included=atof(value.c_str());
                        else if (key == "get_min_energy_within_top_percent") get_min_energy_within_top_percent=atof(value.c_str());
                        else if (key == "get_percent_within_min_energy") get_percent_within_min_energy=atof(value.c_str());
                        else if (key == "allow_Tvib_less_than_Trot") allow_Tvib_less_than_Trot=atoi(value.c_str());
                        else if (key == "Boltzmann_j_averaging") Boltzmann_j_averaging=atoi(value.c_str());
                    }
                }
            }
        }
        infile.close();
	}
}
//***********************************************************************************************************************************************************************************
class Phasespace_Core_Input
{
    public:
        std::string job_name, job_path;
        bool output_to_file, output_histograms, use_mpi;
        double E_total, E_available;

        Phasespace_Core_Input();
        virtual Phasespace_Core_Input* clone(){return new Phasespace_Core_Input(*this);};
        virtual ~Phasespace_Core_Input(){};
    //private:
        Phasespace_Core_Input(const Phasespace_Core_Input &rhs);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & job_name & job_path & output_to_file & output_histograms & use_mpi & E_total & E_available;}
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Core_Input::Phasespace_Core_Input(){
    output_to_file=output_histograms=use_mpi=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Core_Input::Phasespace_Core_Input(const Phasespace_Core_Input &rhs){
    job_name=rhs.job_name;
    job_path=rhs.job_path;
    output_to_file=rhs.output_to_file;
    output_histograms=rhs.output_histograms;
    use_mpi=rhs.use_mpi;
    E_total=rhs.E_total;
    E_available=rhs.E_available;
}
//***********************************************************************************************************************************************************************************



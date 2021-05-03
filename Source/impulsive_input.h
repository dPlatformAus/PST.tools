//***********************************************************************************************************************************************************************************
class Impulsive_General_Input : public General_Input
{
    public:
        unsigned wavefunction_units; //default UNITS_GAUSSIAN 1
        double wavefunction_step_size; //degrees
        double wavefunction_zero_threshold;
        bool HO_Wavefunctions_write_to_file;
        unsigned tunneling_model, impulsive_J_mode;

        Impulsive_General_Input();
        virtual ~Impulsive_General_Input(){};
        virtual void parse_input(const std::string &file_name);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Impulsive_General_Input::Impulsive_General_Input():General_Input(){
    tunneling_model=1; //default is Eckart tunneling model
    impulsive_J_mode=1; //default is adding in Quadrature
    wavefunction_units=1;
    wavefunction_step_size=2;
    wavefunction_zero_threshold=1E-5;
    HO_Wavefunctions_write_to_file=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_General_Input::parse_input(const std::string &file_name){
    General_Input::parse_input(file_name); //first the parent class has to parse the file
    const unsigned NO_SECTION=0;
    const unsigned GENERAL=1;
	unsigned read_state=NO_SECTION;
	std::string file_name_and_path=input_directory+file_name;
    std::ifstream infile (file_name_and_path.c_str());
	std::string line, key, value;
    semipersistent_output << "Impulsive_General_Input parsing" << std::endl;
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
                        if (key == "wavefunction_units") wavefunction_units=atoi(value.c_str());
                        else if (key == "wavefunction_step_size") wavefunction_step_size=atof(value.c_str());
                        else if (key == "wavefunction_zero_threshold") wavefunction_zero_threshold=atof(value.c_str());
                        else if (key == "HO_Wavefunctions_write_to_file") HO_Wavefunctions_write_to_file=atoi(value.c_str());
                        else if (key == "tunneling_model") tunneling_model=atoi(value.c_str());
                        else if (key == "impulsive_J_mode") impulsive_J_mode=atoi(value.c_str());
                    }
                }
            }
        }
        infile.close();
	}
}
//***********************************************************************************************************************************************************************************
class Impulsive_Velocity_Angle_Distribution_Element
{
    public:
        double radians, probability;
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & radians & probability;}
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Impulsive_Bond_Angle_Distribution_Element
{
    public:
        double f1_radians, f2_radians, f1_mo_inertia, f2_mo_inertia, probability, f1_E_nk_imp, f2_E_nk_imp;
        unsigned f1_n_imp, f2_n_imp, f1_k_imp, f2_k_imp;
        Impulsive_Bond_Angle_Distribution_Element():f1_radians(0),f2_radians(0),f1_mo_inertia(0),f2_mo_inertia(0),probability(0),f1_E_nk_imp(0),f2_E_nk_imp(0),f1_n_imp(0),f2_n_imp(0),f1_k_imp(0),f2_k_imp(0){}
        bool operator==(const Impulsive_Bond_Angle_Distribution_Element &rhs);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & f1_radians & f2_radians & f1_mo_inertia & f2_mo_inertia & probability & f1_E_nk_imp & f2_E_nk_imp;
            ar & f1_n_imp & f2_n_imp & f1_k_imp & f2_k_imp;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Impulsive_Bond_Angle_Distribution_Element::operator==(const Impulsive_Bond_Angle_Distribution_Element &rhs){
    if (f1_n_imp==rhs.f1_n_imp && f1_k_imp==rhs.f1_k_imp && f2_n_imp==rhs.f2_n_imp && f2_k_imp==rhs.f2_k_imp) return true;
    else return false;
}
//***********************************************************************************************************************************************************************************
class Impulsive_Core_Input : public Phasespace_Core_Input
{
    public:
        unsigned tunneling_model; //0=none, 1=Eckart
        unsigned impulsive_J_mode; //0= end to end. 1=quadrature (default)
        std::vector<Impulsive_Bond_Angle_Distribution_Element> bond_angle_distribution;

        Impulsive_Core_Input();
        virtual Impulsive_Core_Input* clone(){return new Impulsive_Core_Input(*this);};
        virtual ~Impulsive_Core_Input(){};
    //private:
        Impulsive_Core_Input(const Impulsive_Core_Input &rhs);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Phasespace_Core_Input>(*this);
            ar & tunneling_model & impulsive_J_mode & bond_angle_distribution;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Impulsive_Core_Input::Impulsive_Core_Input():Phasespace_Core_Input(){
    tunneling_model=1; //default is Eckart tunneling model
    impulsive_J_mode=1; //default is adding in Quadrature
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Impulsive_Core_Input::Impulsive_Core_Input(const Impulsive_Core_Input &rhs):Phasespace_Core_Input(rhs){
    tunneling_model=rhs.tunneling_model;
    impulsive_J_mode=rhs.impulsive_J_mode;
    bond_angle_distribution=rhs.bond_angle_distribution;
}
//***********************************************************************************************************************************************************************************

#define SPHERICAL_ANGLE_DISTRIBUTION 1
#define PLANAR_ANGLE_DISTRIBUTION 2

//***********************************************************************************************************************************************************************************
class Triple_Fragmentation_General_Input : public Impulsive_General_Input
{
    public:
        unsigned velocity_angle_distribution_size, velocity_angle_distribution_type;
        bool use_secondary_states_histograms;

        Triple_Fragmentation_General_Input();
        virtual ~Triple_Fragmentation_General_Input(){};
        virtual void parse_input(const std::string &file_name);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Triple_Fragmentation_General_Input::Triple_Fragmentation_General_Input():Impulsive_General_Input(){
    velocity_angle_distribution_size=720;
    use_secondary_states_histograms=0;
    velocity_angle_distribution_type=SPHERICAL_ANGLE_DISTRIBUTION; //this is only the default at the moment because it was the first method we used ... for backwards compatibility ... whether it should be the default long term remains to be seen
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Fragmentation_General_Input::parse_input(const std::string &file_name){
    Impulsive_General_Input::parse_input(file_name); //first the base class has to read the file
    const unsigned NO_SECTION=0;
    const unsigned GENERAL=1;
	unsigned read_state=NO_SECTION;
	std::string file_name_and_path=input_directory+file_name;
    std::ifstream infile (file_name_and_path.c_str());
	std::string line, key, value;
    semipersistent_output << "Triple_Fragmentation_General_Input parsing" << std::endl;
    print_semipersistent_output();
	if (infile.is_open()) {
        while(getline(infile,line)) {
            if (line.substr(0,2) != "//"){ //if this is not a comment line (we should ignore comment lines)
                if (line.substr(0,8) == "--------"){ //this is the start of a new section
                    read_state=NO_SECTION; //set to none first in case the section is not recognised by this parse routine
                    if (line.find("GENERAL") != std::string::npos) read_state=GENERAL;
                }
                if (read_state==GENERAL){
                    read_key_value(line, key, value);
                    if (key != "") {
                        if (key == "velocity_angle_distribution_size") velocity_angle_distribution_size=atoi(value.c_str());//error checking will need to confirm this is > 0
                        if (key == "velocity_angle_distribution_type") velocity_angle_distribution_type=atoi(value.c_str());//error checking will need to confirm this is > 0
                        else if (key == "use_secondary_states_histograms") use_secondary_states_histograms=atoi(value.c_str());
                    }
                }
            }
        }
        infile.close();
	}
}
//***********************************************************************************************************************************************************************************
class Triple_Fragmentation_Core_Input : public Impulsive_Core_Input
{
    public:
        std::vector<Impulsive_Velocity_Angle_Distribution_Element> velocity_angle_distribution;

        Triple_Fragmentation_Core_Input();
        virtual Triple_Fragmentation_Core_Input* clone(){return new Triple_Fragmentation_Core_Input(*this);};
        virtual ~Triple_Fragmentation_Core_Input(){};
        void Initialise_Default_Velocity_Angle_Distribution(Triple_Fragmentation_General_Input *t_general_input);
    //private:
        Triple_Fragmentation_Core_Input(const Triple_Fragmentation_Core_Input &rhs);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Impulsive_Core_Input>(*this);
            ar & velocity_angle_distribution;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Triple_Fragmentation_Core_Input::Triple_Fragmentation_Core_Input() : Impulsive_Core_Input(){
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Triple_Fragmentation_Core_Input::Triple_Fragmentation_Core_Input(const Triple_Fragmentation_Core_Input &rhs) : Impulsive_Core_Input(rhs){
    velocity_angle_distribution=rhs.velocity_angle_distribution;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Fragmentation_Core_Input::Initialise_Default_Velocity_Angle_Distribution(Triple_Fragmentation_General_Input *t_general_input){
    double theta, d_theta, sum;
    unsigned i;
    Impulsive_Velocity_Angle_Distribution_Element temp_angle_distribution_element;
    sum=0;
    d_theta=PI/(1+t_general_input->velocity_angle_distribution_size);
    for (theta=d_theta; theta<PI; theta+=d_theta) {
        temp_angle_distribution_element.radians=theta;
        if (t_general_input->velocity_angle_distribution_type==SPHERICAL_ANGLE_DISTRIBUTION) temp_angle_distribution_element.probability=sin(theta); //many more ways for vectors to be at 90 degrees with a spherical distro
        else if (t_general_input->velocity_angle_distribution_type==PLANAR_ANGLE_DISTRIBUTION) temp_angle_distribution_element.probability=1; //equal probability for all angles in a planar distro
        else persistent_output<<"ERROR: velocity_angle_distribution_type is undefined!!"<<std::endl;
        sum+=temp_angle_distribution_element.probability;
        velocity_angle_distribution.push_back(temp_angle_distribution_element);
    }
    for (i=0; i<velocity_angle_distribution.size(); i++) velocity_angle_distribution[i].probability /= sum;
}
//***********************************************************************************************************************************************************************************

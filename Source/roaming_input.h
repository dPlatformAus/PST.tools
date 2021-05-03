//***********************************************************************************************************************************************************************************
class Roaming_General_Input : public General_Input
{
    public:
        Number_Range<double> delta_E_roam_range; // cm-1 //I can't decide whether this belongs here in General (becuase it's a parameter of the model) or whether this is a property of the reactant and should go into the molecule class.
        double delta_E_roam_increment; // cm-1
        Number_Range<double> P_roam_range; // dimensionless between 0 and 1 -- default is 1
        double P_roam_increment; // dimensionless between 0 and 1 -- default is 0.1
        bool fit_quantum_yield;
        Fit_Point_List fit_points;

        Roaming_General_Input();
        virtual ~Roaming_General_Input(){};
        virtual void parse_input(const std::string &file_name);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_General_Input::Roaming_General_Input():General_Input(){
    delta_E_roam_increment=1;
    P_roam_range.setRange(1,1); //default to 1, not the default for a Double_Range, which would be 0
    P_roam_increment=0.1;
    fit_quantum_yield=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_General_Input::parse_input(const std::string &file_name){
    General_Input::parse_input(file_name); //first the parent class has to parse the file
    const unsigned NO_SECTION=0;
    const unsigned GENERAL=1;
    const unsigned ROAMING_FIT_POINTS=2;
    const unsigned RADICAL_FIT_POINTS=3;
	unsigned read_state=NO_SECTION;
	std::string file_name_and_path=input_directory+file_name;
    std::ifstream infile (file_name_and_path.c_str());
	std::string line, key, value;
	double energy, target;
    semipersistent_output << "Roaming_General_Input parsing" << std::endl;
	if (infile.is_open()) {
        while(getline(infile,line)) {
            if (line.substr(0,2) != "//"){ //if this is not a comment line (we should ignore comment lines)
                if (line.substr(0,8) == "--------"){ //this is the start of a new section
                    read_state=NO_SECTION; //set to none first in case the section is not recognised by this parse routine
                    if (line.find("-GENERAL-") != std::string::npos) read_state=GENERAL;
                    else if (line.find("-ROAMING FIT POINTS-") != std::string::npos) read_state=ROAMING_FIT_POINTS;
                    else if (line.find("-RADICAL FIT POINTS-") != std::string::npos) read_state=RADICAL_FIT_POINTS;
                }
                if (read_state==GENERAL){
                    read_key_value(line, key, value);
                    if (key != "") {
                        if (key == "delta_E_roam_range") delta_E_roam_range.setRange(value);
                        else if (key == "delta_E_roam_increment") delta_E_roam_increment=atof(value.c_str());
                        else if (key == "P_roam_range") P_roam_range.setRange(value);
                        else if (key == "P_roam_increment") P_roam_increment=atof(value.c_str());
                        else if (key == "fit_quantum_yield") fit_quantum_yield=atoi(value.c_str());
                    }
                }else if (read_state==ROAMING_FIT_POINTS){
                    if (read_fit_point_value(line, energy, target)) fit_points.add_roaming(energy, target);
                }else if (read_state==RADICAL_FIT_POINTS){
                    if (read_fit_point_value(line, energy, target)) fit_points.add_radical(energy, target);
                }
            }
        }
        infile.close();
	}
}
//***********************************************************************************************************************************************************************************
class Roaming_Core_Input : public Phasespace_Core_Input /**
Currently this derivation of Phasespace_Core_Input is unnessesary
However, we may write new roaming type calculations in the future that may need additional inputs for the roaming core
Obviously, this would require a modification of Roaming_Phasespace_Core also
If we never do, this should be removed and Roaming_Phasespace_Core should simply use Phasespace_Core_Input **/
{
    public:

        Roaming_Core_Input();
        virtual Roaming_Core_Input* clone(){return new Roaming_Core_Input(*this);};
        virtual ~Roaming_Core_Input(){};

    //private:
        Roaming_Core_Input(const Roaming_Core_Input &rhs);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Core_Input::Roaming_Core_Input() : Phasespace_Core_Input(){
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Core_Input::Roaming_Core_Input(const Roaming_Core_Input &rhs) : Phasespace_Core_Input(rhs){
}
//***********************************************************************************************************************************************************************************


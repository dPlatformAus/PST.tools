//***********************************************************************************************************************************************************************************
class Roaming_Molecule : public Phasespace_Molecule
{
    public:
        double delta_E_roam; //cm-1
        double P_roam;
        Roaming_Histogram_Use *r_is_on, *r_probe_is_on;

        Roaming_Molecule(){};
        virtual Roaming_Molecule* clone(){return new Roaming_Molecule(*this);};
        virtual ~Roaming_Molecule(){}; //delete is done by parent
        virtual void construct();
        virtual void create_histogram_use();
        virtual void parse_input(const std::string &file_name, const std::string &section_name);

    private:
        Roaming_Molecule(const Roaming_Molecule &rhs);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Phasespace_Molecule>(*this);
            ar & delta_E_roam & P_roam;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Molecule::construct(){
    Phasespace_Molecule::construct();
    delta_E_roam=0;
    P_roam=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Molecule::Roaming_Molecule(const Roaming_Molecule &rhs) : Phasespace_Molecule(rhs){
    delta_E_roam=rhs.delta_E_roam;
    P_roam=rhs.P_roam;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Molecule::create_histogram_use(){
    is_on = new Roaming_Histogram_Use();
    probe_is_on = new Roaming_Histogram_Use();
    r_is_on = dynamic_cast<Roaming_Histogram_Use*>(is_on);
    r_probe_is_on = dynamic_cast<Roaming_Histogram_Use*>(probe_is_on);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Molecule::parse_input(const std::string &file_name, const std::string &section_name){
    Phasespace_Molecule::parse_input(file_name, section_name); //first parse the parent class's members
    const unsigned NO_SECTION=0;
    const unsigned GENERAL=1;
	unsigned read_state=NO_SECTION;
	std::string file_name_and_path=input_directory+file_name;
    std::ifstream infile (file_name_and_path.c_str());
	std::string line, key, value;
    semipersistent_output << "Roaming_Molecule "<< number <<" parsing" << std::endl;
	if (infile.is_open()) {
        while(getline(infile,line)) {
            if (line.substr(0,2) != "//"){ //if this is not a comment line (we should ignore comment lines)
                if (line.substr(0,8) == "--------"){ //this is the start of a new section
                    read_state=NO_SECTION; //set to none first in case the section is not recognised by this parse routine
                    if (line.find( std::string("-")+section_name+std::string("-") ) != std::string::npos) read_state=GENERAL;
                }
                if (read_state==GENERAL){
                    read_key_value(line, key, value);
                    if (key != "") {
                        if (key == "delta_E_roam") delta_E_roam=atof(value.c_str());
                        else if (key == "P_roam") P_roam=atof(value.c_str());
                    }
                }
            }
        }
        infile.close();
	}
}
//***********************************************************************************************************************************************************************************


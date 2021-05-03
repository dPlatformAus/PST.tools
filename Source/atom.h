class Element
{
    public:
        int atomic_number;
        std::string symbol;
        std::string name;
		double mass;	// atomic mass of atom
		double radius;  //van der Waal's radius of atom

};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Element_List
{
    public:
        Element *elements;
        int num_elements;
        Element_List();
        ~Element_List();
        double get_radius(const std::string atom_symbol);
        Element get_element(const std::string atom_symbol);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Element_List::Element_List() {
    Input_File in_elements;
    int i;
    char buffer [100], *temp;
    std::vector<std::string> results;
    bool read_header=0;

    num_elements=0;
    in_elements.set_input_directory("config/");
    in_elements.open("elements.csv");   //count the number of elements in the file
    if (in_elements.fp == NULL) fprintf(stderr, "Couldn't open the elements.csv file, errno=%d\n", errno);
    else{
        while ( ! feof (in_elements.fp) ){
            temp=fgets (buffer , 100 , in_elements.fp);
            if (read_header) num_elements++;
            read_header=1;
        }
    }
    in_elements.close();
    elements = new Element[num_elements];
    in_elements.open("elements.csv");   //populate the element list
    if (in_elements.fp == NULL) fprintf(stderr, "Couldn't open the elements.csv file, errno=%d\n", errno);
    else{
        i=0;
        while ( ! feof (in_elements.fp) ){
            if (read_header) temp=fgets (buffer , 100 , in_elements.fp); //will be 1 because set that way in previous loop
            read_header=0;
            temp=fgets (buffer , 100 , in_elements.fp);
            boost::split(results, buffer, boost::is_any_of(","));
            elements[i].mass=strtod(results[0].c_str(),NULL);
            elements[i].name=results[1];
            elements[i].symbol=results[2];
            elements[i].atomic_number=atoi(results[3].c_str());
            elements[i].radius=strtod(results[4].c_str(),NULL);
            i++;
        }
    }
    in_elements.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Element_List::~Element_List() {
    delete [] elements;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Element_List::get_radius(const std::string atom_symbol) {
    int i;
    double temp_result=0;
    for (i=0;i<num_elements;i++) if (elements[i].symbol==atom_symbol) temp_result=elements[i].radius;
    return temp_result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Element Element_List::get_element(const std::string atom_symbol) {
    int i;
    Element temp_result;
    for (i=0;i<num_elements;i++) {
        if (elements[i].symbol==atom_symbol) {
            temp_result=elements[i];
        }
    }
    return temp_result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Element_List periodic_table;
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Atom : public Element
{
    public:
		double coord[3]; // the position coordinates, 0=x, 1=y, 2=z

		void assign_default_radius(); //from old file format parser
		void set_element(const std::string atom_symbol);
		void set_coord(const double in_x, const double in_y, const double in_z);
		void set_atom(const std::string atom_symbol, const double in_x, const double in_y, const double in_z);

};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Atom::set_element(const std::string atom_symbol) {
    Element this_element;
    this_element=periodic_table.get_element(atom_symbol);
    atomic_number=this_element.atomic_number;
	symbol=this_element.symbol;
	name=this_element.name;
	mass=this_element.mass;
	radius=this_element.radius;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Atom::set_coord(const double in_x, const double in_y, const double in_z) {
    coord[0]=in_x;
    coord[1]=in_y;
    coord[2]=in_z;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Atom::set_atom(const std::string atom_symbol, const double in_x, const double in_y, const double in_z){
    set_element(atom_symbol);
    set_coord(in_x,in_y,in_z);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Atom::assign_default_radius() {
    //assigns a default radius to the atom based on it's mass.
    //if user radii are read from the file the user raduis will overwrite the radius assigned here.
    radius=0; //when code has to use the radii it will check them and throw an error if one is zero
    if (mass < 3) radius=periodic_table.get_radius("H");
    if ((int)mass == 12) radius=periodic_table.get_radius("C");
    if ((int)mass == 14) radius=periodic_table.get_radius("N");
    if ((int)mass == 16) radius=periodic_table.get_radius("O");
    if ((int)mass == 19) radius=periodic_table.get_radius("F");
    if ((int)mass == 28) radius=periodic_table.get_radius("Si");
    if ((int)mass == 31) radius=periodic_table.get_radius("P");
    if ((int)mass == 32) radius=periodic_table.get_radius("S");
    if (mass > 35 && mass < 38) radius=periodic_table.get_radius("Cl");
    if ((int)mass == 73) radius=periodic_table.get_radius("Ge");
    if ((int)mass == 75) radius=periodic_table.get_radius("As");
    if ((int)mass == 79) radius=periodic_table.get_radius("Se");
    if ((int)mass == 80) radius=periodic_table.get_radius("Br");
    if ((int)mass == 122) radius=periodic_table.get_radius("Sb");
    if ((int)mass == 127) radius=periodic_table.get_radius("I");
}

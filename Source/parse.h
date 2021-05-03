//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void load_input_files(std::vector<std::string> &input_files){
	DIR *dp;
	struct dirent *ep;
	int temp;
	dp = opendir (input_directory.c_str());
	if (dp != NULL) {
		while ((ep = readdir (dp))){
			if ( !( (ep->d_name[0] == '.' && NAMLEN(ep) == 1) || (ep->d_name[1] == '.' && NAMLEN(ep) == 2) ) ) {
				if (get_file_extension(ep->d_name)==input_suffix) input_files.push_back(ep->d_name); //we will run all ".in" files
			}
		}
		closedir (dp);
	} else {
        perror ("Couldn't open the input directory\n");
		temp=system("PAUSE");
		exit(1);
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void read_key_value(const std::string line, std::string &key, std::string &value){
    size_t equals_position;
    equals_position = line.find("=");
    if (equals_position != std::string::npos){ // then this is a key value pair
        key = line.substr(0, equals_position);
        value = line.substr(equals_position + 1);
        boost::trim(key);
        boost::trim(value);
    }else{//not key value pair
        key=value="";
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool read_atom_value(const std::string line, std::string &symbol, double &ax, double &ay, double &az){
    bool result = true;
    std::stringstream ss;
    std::string value;
    ss << line;
    if (getline(ss,value,',')) {
        symbol=value;
    }else{
        result=false;
    }
    if (getline(ss,value,',')) {
        ax=atof(value.c_str());
    }else{
        result=false;
    }
    if (getline(ss,value,',')) {
        ay=atof(value.c_str());
    }else{
        result=false;
    }
    if (getline(ss,value,',')) {
        az=atof(value.c_str());
    }else{
        result=false;
    }
    /*if(result){

        persistent_output<<"symbol="<<symbol.c_str()<<" ax="<<ax<<" ay="<<ay<<" az="<<az<<std::endl;
        print_persistent_output();
    }*/
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool read_mode_value(const std::string line, double &frequency, unsigned &degeneracy) {
    bool result = true;
    std::stringstream ss;
    std::string value;
    ss << line;
    if (getline(ss,value,',')) {
        frequency=atof(value.c_str());
    }else{
        result=false;
    }
    if (getline(ss,value,',')) {
        degeneracy=atoi(value.c_str());
    }else{
        result=false;
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool read_spin_orbit_value(const std::string line, double &energy, float &j) {
    bool result = true;
    std::stringstream ss;
    std::string value;
    ss << line;
    if (getline(ss,value,',')) {
        energy=atof(value.c_str());
    }else{
        result=false;
    }
    if (getline(ss,value,',')) {
        j=atof(value.c_str());
    }else{
        result=false;
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool read_fit_point_value(const std::string line, double &energy, double &target) {
    bool result = true;
    std::stringstream ss;
    std::string value;
    ss << line;
    if (getline(ss,value,',')) {
        energy=atof(value.c_str());
    }else{
        result=false;
    }
    if (getline(ss,value,',')) {
        target=atof(value.c_str());
    }else{
        result=false;
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void read_index_list(const std::string line, std::vector<unsigned> &the_list) {
    std::stringstream ss;
    std::string value;
    ss << line;
    while (getline(ss,value,',')) the_list.push_back(atoi(value.c_str()));
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

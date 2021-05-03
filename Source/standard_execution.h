//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void standard_execution(const std::string in_file_name){
#if (defined(USING_MPI))
    mpi::communicator world;
#endif
    Base_Calculation *the_calculation;
    bool calc_type_recognised(1); //will be set false if cal_type isn't recognised
    std::string file_name_and_path=input_directory+in_file_name;
    std::string line, key, value, calc_type("");
	std::ifstream infile(file_name_and_path.c_str());
	if (infile.is_open()){
        persistent_output << "****** parsing: " << in_file_name<<std::endl;
        print_persistent_output();
        while(getline(infile,line)){
            read_key_value(line, key, value);
            if (key == "calc_type") calc_type=value;
        }
        infile.close();
        if (calc_type!=""){//make the appropriate calculation object
#if (defined(USING_MPI))
            broadcast(world, calc_type, 0);//tell slave processes which function to run
#endif
            if (calc_type == "phasespace") the_calculation = new Phasespace_Calculation();
            else if (calc_type == "impulsive") the_calculation = new Impulsive_Calculation();
            else if (calc_type == "triple fragmentation") the_calculation = new Triple_Fragmentation_Calculation();
            //else if (calc_type == "roaming") the_calculation = new Roaming_Calculation();
            else if (calc_type == "roaming fraction energy dependence") the_calculation = new Roaming_Fraction_Energy_Dependance_Calculation();
            else if (calc_type == "roaming fit delta_E_roam and P_roam") the_calculation = new Roaming_Fit_Delta_E_Roam_And_P_Roam_Calculation();
            else if (calc_type == "roaming fit Trot and Tvib and P_roam") the_calculation = new Roaming_Fit_Trot_And_Tvib_And_P_Roam_Calculation();
            else if (calc_type == "rebin") the_calculation = new Rebin_Calculation();
            /*else if (calc_type == "get vibrational states") the_calculation = new Phasespace_Get_Vibrational_States_Calculation();
            else if (calc_type == "get rotational states") the_calculation = new Phasespace_Get_Rotational_States_Calculation();*/
            //else if (calc_type == "Eckart plot") the_calculation = new Eckart_Plot_Calculation();
            else calc_type_recognised=0; //the calc_type in the file doesn't correspond to any of the objects this program knows about, so we didn't make a calculation object
            if (calc_type_recognised) {
                the_calculation->go(in_file_name);
                delete the_calculation; //if we made a calculation object delete it
            }
            else persistent_output << " - ERROR: invalid calc_type \""<<calc_type<<"\" in file \""<<file_name_and_path<<"\""<<std::endl;
        }else persistent_output << " - ERROR: calc_type missing from file \""<<file_name_and_path<<"\""<<std::endl;
    }else persistent_output << " - ERROR: unable to open file \""<<file_name_and_path<<"\""<<std::endl;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void do_standard_execution(){
#if (defined(USING_MPI))
    mpi::communicator world;
#endif
    unsigned i, num_files;
    std::vector<std::string> input_files;
    load_input_files(input_files);
    num_files=input_files.size();
#if (defined(USING_MPI))
    broadcast(world, num_files, 0);
#endif
    if (num_files){
        persistent_output<<"Processing "<<num_files<<" jobs:"<<std::endl;
        print_persistent_output();
        for (i=0; i<num_files; i++){
            job_number=i+1;
#if (defined(USING_MPI))
            broadcast(world, job_number, 0);
#endif
            persistent_output<<std::endl<<"job "<<(i+1)<<" ";
            standard_execution(input_files[i]);
        }
        print_persistent_output();
    }else std::cout << "no input files detected" << std::endl;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


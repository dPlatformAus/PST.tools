//***********************************************************************************************************************************************************************************
class Secondary_States_Line
{
    public:
        double primary_E_internal, primary_degeneracy, primary_vel;
        unsigned primary_J;
        Secondary_States_Line():primary_E_internal(0),primary_degeneracy(0),primary_vel(0),primary_J(0){}
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & primary_E_internal & primary_degeneracy & primary_vel & primary_J;}
};
//***********************************************************************************************************************************************************************************
class Triple_Fragmentation_Calculation : public Impulsive_Calculation{
    public:
        Triple_Fragmentation_General_Input *t_general_input;
        Impulsive_Molecule *frag1_S1, *frag1_S2, *frag2_S1, *frag2_S2;

        virtual ~Triple_Fragmentation_Calculation();
        virtual void create_general_input();
        virtual void create_general_input_pointer();
        virtual void create_fragments();
        virtual void initialise(const std::string &file_name);
        virtual bool run();
        virtual void secondary_fragmentation(Primary_Dissociation_Result *primary_PST_output, Impulsive_Molecule *dissociating_parent,
                                                               Impulsive_Molecule *dissociating_frag1, Impulsive_Molecule *dissociating_frag2, Histogram_3d *sec_states);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Triple_Fragmentation_Calculation::~Triple_Fragmentation_Calculation(){
    delete frag1_S1;
    delete frag1_S2;
    delete frag2_S1;
    delete frag2_S2;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Fragmentation_Calculation::create_general_input(){
    the_general_input = new Triple_Fragmentation_General_Input();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Fragmentation_Calculation::create_general_input_pointer(){
    Impulsive_Calculation::create_general_input_pointer();
    t_general_input = dynamic_cast<Triple_Fragmentation_General_Input*>(the_general_input);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Fragmentation_Calculation::create_fragments(){
    Impulsive_Calculation::create_fragments();
    frag1_S1 = new Impulsive_Molecule();
    frag1_S1->construct();
    frag1_S2 = new Impulsive_Molecule();
    frag1_S2->construct();
    frag2_S1 = new Impulsive_Molecule();
    frag2_S1->construct();
    frag2_S2 = new Impulsive_Molecule();
    frag2_S2->construct();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Fragmentation_Calculation::initialise(const std::string &file_name){
    Impulsive_Calculation::initialise(file_name);
    frag1_S1->number=3;
    frag1_S2->number=4;
    frag2_S1->number=5;
    frag2_S2->number=6;
    frag1_S1->parse_input(file_name, "FRAGMENT 1 SECONDARY FRAGMENT 1");
    frag1_S2->parse_input(file_name, "FRAGMENT 1 SECONDARY FRAGMENT 2");
    frag2_S1->parse_input(file_name, "FRAGMENT 2 SECONDARY FRAGMENT 1");
    frag2_S2->parse_input(file_name, "FRAGMENT 2 SECONDARY FRAGMENT 2");
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Fragmentation_Calculation::secondary_fragmentation(Primary_Dissociation_Result *primary_PST_output, Impulsive_Molecule *dissociating_parent,
                                                               Impulsive_Molecule *dissociating_frag1, Impulsive_Molecule *dissociating_frag2, Histogram_3d *sec_states){
#if (defined(USING_MPI))
    mpi::communicator world;
    unsigned dest_thread(world.size()-1), portion, last_portion;
    dMPI_Master_Progress_Bar the_master_progress_bar;
    bool sending_secondary_results(0);
#endif

    Secondary_States_Line line_values, last_line_values;
    std::vector<Secondary_States_Line> portion_lines;
    double secondary_threshold;
    unsigned line_number(0), i,j,k,w, progress_bar_length(1), num_lines;
    unsigned thread_rank(0), num_threads(1);
    double process_energy_tally(0), process_energy_increment(0);
    double final_primary_degeneracy(0), primary_vel_max(0), secondary_vel_max;
    bool not_first_line(0), temp_bool;
    Output_File compressed_secondary_states;
    std::ifstream secondary_state_list, compressed_secondary_state_list;
    std::string line, last_line(""), value, secondary_state_list_filename, dest_path_and_name, wf_file_name, frag_num_string;
    std::stringstream ss;
    std::string bitmap_path_and_name_start;
//    HO_Wavefunctions *bend_wf;
    Secondary_Dissociation_Result *temp_secondary_result, *secondary_result, *grand_secondary_result;
    Impulsive_Bond_Angle_Distribution_Element temp_angle_distribution_element;
    Triple_Fragmentation_Core_Input *secondary_PST_input;
    Histogram_Limits secondary_limits;
    Triple_Frag_Secondary *the_secondary_PST;
    secondary_PST_input = new Triple_Fragmentation_Core_Input();
    secondary_result = new Secondary_Dissociation_Result();
    grand_secondary_result = new Secondary_Dissociation_Result();
    temp_secondary_result = new Secondary_Dissociation_Result();
    bitmap_path_and_name_start=job_path+"3F_Secondary_";
    //these three can probably be in a function set_standard_PST_input_values() or something
    secondary_PST_input->job_path=job_path;
    secondary_PST_input->output_to_file=0;// no point overwriting constantly ... just takes much longer!  general_input->objects_write_to_file;
    secondary_PST_input->output_histograms=t_general_input->histogram_write_to_file;
    secondary_PST_input->job_name="secondary_fragmentation";
    secondary_PST_input->Initialise_Default_Velocity_Angle_Distribution(t_general_input);
    secondary_PST_input->tunneling_model=i_general_input->tunneling_model; //later we may allow different models for each fragment and then we would use the general value if the fragment(parent) didn't have this value set (and we would add this variable to the impulsive molecule class too)
    secondary_PST_input->impulsive_J_mode=t_general_input->impulsive_J_mode;
    secondary_threshold=i_reactant->E_dissociation+dissociating_parent->E_dissociation; //the threshold energy to get these secondary products from the bottom of the primary well
    if(t_general_input->excitation_energy_range.min > secondary_threshold){ //secondary dissociation is asseccible ... get results for this
        dissociating_frag1->initialise(dissociating_frag2->mass);
        dissociating_frag2->initialise(dissociating_frag1->mass);
        frag_num_string=boost::lexical_cast<std::string>(dissociating_parent->number);
        secondary_limits.E_max = t_general_input->excitation_energy_range.min - secondary_threshold;
        if (dissociating_parent->number == 1) primary_vel_max=primary_PST_output->f1_vel_max;
        else if (dissociating_parent->number == 2) primary_vel_max=primary_PST_output->f2_vel_max;
        else persistent_output<<"ERROR: incorrect dissociating_parent argument in call to Triple_Fragmentation_Calculation::secondary_fragmentation."<<std::endl; ///this needs to return 0 and make the result a fail! (or something)
        if (dissociating_frag1->mass==i_frag1->mass) secondary_limits.vel_max1=primary_PST_output->data->frag1_vel.x_max; //use the same histogram parameters as were used for this fragment type (sec1=prim1) in the primary dissociation
        else if (dissociating_frag1->mass==i_frag2->mass) secondary_limits.vel_max1=primary_PST_output->data->frag2_vel.x_max; //use the same histogram parameters as were used for this fragment type (sec1=prim2) in the primary dissociation
        else {
            secondary_vel_max=sqrt(secondary_limits.E_max*dissociating_frag1->get_vel_mass_factor(dissociating_frag2->mass)); //calculate based on secondary fragment masses and available energy
            secondary_limits.vel_max1=secondary_vel_max+primary_vel_max; //max will be vector sum of primary and secondary velocities in the same direction
        }
        if (dissociating_frag2->mass==i_frag1->mass) secondary_limits.vel_max2=primary_PST_output->data->frag1_vel.x_max; //use the same histogram parameters as were used for this fragment type (sec2=prim1) in the primary dissociation
        if (dissociating_frag2->mass==i_frag2->mass) secondary_limits.vel_max2=primary_PST_output->data->frag2_vel.x_max; //use the same histogram parameters as were used for this fragment type (sec2=prim2) in the primary dissociation
        else {
            secondary_vel_max=sqrt(secondary_limits.E_max*dissociating_frag2->get_vel_mass_factor(dissociating_frag1->mass)); //calculate based on secondary fragment masses and available energy
            secondary_limits.vel_max2=secondary_vel_max+primary_vel_max; //max will be vector sum of primary and secondary velocities in the same direction
        }
        secondary_limits.E_max1=dissociating_frag1->get_E_trans(secondary_limits.vel_max1);
        secondary_limits.E_max2=dissociating_frag2->get_E_trans(secondary_limits.vel_max2);
 /*//**************************************************************************        angle distro
        if (dissociating_parent->use_angle_distribution) {
            if (dissociating_parent->use_angle_distribution==BOND_ANGLE_DISTRIBUTION_GROUND) primary_bend_n=0;
            else if (dissociating_parent->use_angle_distribution==BOND_ANGLE_DISTRIBUTION_PARENT) {
                if (dissociating_parent->number == 1) primary_bend_n=primary_PST_output->f1_max_bend_n;
                else if (dissociating_parent->number == 2) primary_bend_n=primary_PST_output->f2_max_bend_n;
            }
            wf_file_name="frag"+frag_num_string+"_bend_wf.csv";
            bend_wf = new HO_Wavefunctions(wf_file_name, job_path, i_general_input, UNITS_GAUSSIAN, primary_bend_n, dissociating_parent); //get the wavefunctions for the maximum n value for the bend
        }
 //**************************************************************************        angle distro */
        secondary_limits.J_max1=dissociating_frag1->get_n_max(secondary_limits.E_max);
        secondary_limits.J_max2=dissociating_frag2->get_n_max(secondary_limits.E_max);
        dissociating_frag1->initialise_vib_levels(secondary_limits.E_max-dissociating_parent->E_barrier);
        dissociating_frag2->initialise_vib_levels(secondary_limits.E_max-dissociating_parent->E_barrier);
        secondary_limits.v_max1=dissociating_frag1->get_v_max(secondary_limits.E_max-dissociating_parent->E_barrier);
        secondary_limits.v_max2=dissociating_frag2->get_v_max(secondary_limits.E_max-dissociating_parent->E_barrier);
        for (i=0; i<=secondary_limits.v_max1; i++) secondary_limits.frag1_vib_names.push_back(dissociating_frag1->vib_levels.bands[i].name);
        for (i=0; i<=secondary_limits.v_max2; i++) secondary_limits.frag2_vib_names.push_back(dissociating_frag2->vib_levels.bands[i].name);
        secondary_result->initialise(secondary_limits, dissociating_parent);
        temp_secondary_result->initialise(secondary_limits, dissociating_parent);
        grand_secondary_result->initialise(secondary_limits, dissociating_parent);
#if (defined(USING_MPI))
        sending_secondary_results=true;
        broadcast(world, sending_secondary_results, 0);
        broadcast(world, secondary_PST_input, 0);
        broadcast(world, dissociating_parent, 0);
        broadcast(world, dissociating_frag1, 0);
        broadcast(world, dissociating_frag2, 0);
        broadcast(world, secondary_result, 0);
#endif
        if (t_general_input->use_secondary_states_histograms){ //loop over histogram and calculate
#if (defined(USING_MPI))
            thread_rank=world.rank();
            num_threads=world.size();
            the_progress_bar.block(); // block the global progress bar so that long state counts can't start it (slave processes cannot print to screen ... may need to have a silent variable or something that kills all output)
            the_master_progress_bar.initialise(100,1,sec_states->bin_length1, CALCULATION_PROGRESS_BAR,1); // initialise the progress bar (this has no effect on the calculation)
#else
            the_progress_bar.initialise(100,1,sec_states->bin_length1, CALCULATION_PROGRESS_BAR); // initialise the progress bar (this has no effect on the calculation)
#endif
            semipersistent_output << "Assigning Secondary Dissociation States To Processes" << std::endl;
            print_semipersistent_output();
            for (i=thread_rank; i<sec_states->bin_length1; i+=num_threads){
#if (defined(USING_MPI))
                the_master_progress_bar.go(i, CALCULATION_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
#else
                the_progress_bar.go(i, CALCULATION_PROGRESS_BAR);//update progress bar (this has no affect on the calculation result)
#endif
                for (j=0; j<sec_states->bin_length2; j++){
                    for (k=0; k<sec_states->bin_length3; k++){
                        if (sec_states->val(i,j,k)){
                            line_values.primary_degeneracy=sec_states->val(i,j,k);
                            line_values.primary_E_internal=sec_states->resolveIndex1(i)+dissociating_parent->E_dissociation; //add parent E_dissoc because we store EXCESS energy in the 3d histogram not internal energy like in the file
                            line_values.primary_vel=sec_states->resolveIndex2(j);
                            line_values.primary_J=sec_states->resolveIndex3(k); //not sure I like this ... perhaps we need a clone for the input, we can't make the frag const if we do this!
                            secondary_result->total_parent_states_considered+=line_values.primary_degeneracy;
                            portion_lines.push_back(line_values);
                        }
                    }
                }
            }
        }else{ //loop over secondary state file and calculate (for MPI the file will already be stored locally for each thread, so there is no need to broadcast anything else)
            //sort out progress bar stuff later (if it isn't automatically by the result output from primary_PST_output - I think we'll only be passing the result in by the time we're finnished here)
            if (dissociating_parent->number == 1) progress_bar_length=primary_PST_output->f1_num_secondary_master_thread; //mpi causes an issue here as this number is the reduced value but the file is shorter (un-reduced value)!
            else if (dissociating_parent->number == 2) progress_bar_length=primary_PST_output->f2_num_secondary_master_thread; //so we store the master thread count seperately so we can access it here after the reduction
#if (defined(USING_MPI))
            the_progress_bar.block(); // block the global progress bar so that long state counts can't start it (slave processes cannot print to screen ... may need to have a silent variable or something that kills all output)
            the_master_progress_bar.initialise(100,1,progress_bar_length, CALCULATION_PROGRESS_BAR,1); // initialise the progress bar (this has no effect on the calculation)
#else
            the_progress_bar.initialise(100,1,progress_bar_length, CALCULATION_PROGRESS_BAR); // initialise the progress bar (this has no effect on the calculation)
#endif
#if (defined(USING_MPI)) ///need to manage job directories too before we can implement checkpoint files well (although, storing lots of checkpoint files on each node is undesireable, so not a big priority perhaps)
            secondary_state_list_filename=scratch_directory+SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+"_t"+boost::lexical_cast<std::string>(world.rank())+".csv";
#else
            secondary_state_list_filename=job_path+SECONDARY_STATES_OUTPUT_DIRECTORY+SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+".csv";
#endif
            secondary_state_list.open(secondary_state_list_filename.c_str());
#if (defined(USING_MPI)) ///consider using binary files for this too - faster for large files etc
            compressed_secondary_states.set_output_directory(scratch_directory);
            compressed_secondary_states.open(SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+"_cmp_t"+boost::lexical_cast<std::string>(world.rank())+".csv");
#else
            compressed_secondary_states.set_output_directory(job_path+SECONDARY_STATES_OUTPUT_DIRECTORY);
            compressed_secondary_states.open(SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+"_cmp.csv");
#endif
            semipersistent_output << "Compressing Secondary States File" << std::endl; //probably better to just do this when writing the first file ... can't be that hard ... last_state = state etc etc
            print_semipersistent_output();
            getline(secondary_state_list,line); //read the header line in the file
            fprintf(compressed_secondary_states.fp, "%s\n", line.c_str()); //write the header
            num_lines=line_number=0;
            while (getline(secondary_state_list,line)){
#if (defined(USING_MPI))
                the_master_progress_bar.go(++line_number, CALCULATION_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
#else
                the_progress_bar.go(++line_number, CALCULATION_PROGRESS_BAR);//update progress bar (this has no affect on the calculation result)
#endif
                ss << line;
                getline(ss,value,',');
                line_values.primary_E_internal=atof(value.c_str());
                getline(ss,value,',');
                line_values.primary_degeneracy=atof(value.c_str());
                getline(ss,value,',');
                line_values.primary_J=atoi(value.c_str());
                getline(ss,value,',');
                line_values.primary_vel=atof(value.c_str());
                ss.str("");
                ss.clear();
                if (line==last_line) final_primary_degeneracy+=line_values.primary_degeneracy;
                else{
                    if (not_first_line) fprintf(compressed_secondary_states.fp, "%lf,%lf,%d,%lf\n", last_line_values.primary_E_internal, final_primary_degeneracy, last_line_values.primary_J, last_line_values.primary_vel);
                    not_first_line=true;
                    num_lines++;
                    last_line=line;
                    last_line_values=line_values;
                    final_primary_degeneracy=line_values.primary_degeneracy;
                }
            }
            if (not_first_line) {
                fprintf(compressed_secondary_states.fp, "%lf,%lf,%d,%lf\n", last_line_values.primary_E_internal, final_primary_degeneracy, last_line_values.primary_J, last_line_values.primary_vel);
                num_lines++;
            }
            compressed_secondary_states.close();
            remove(secondary_state_list_filename.c_str()); //delete the uncompressed file
#if (defined(USING_MPI))
            the_master_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
            the_master_progress_bar.initialise(100,1,num_lines, CALCULATION_PROGRESS_BAR,2); // initialise the progress bar (this has no effect on the calculation)
#else
            the_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
            the_progress_bar.initialise(100,1,num_lines, CALCULATION_PROGRESS_BAR); // initialise the progress bar (this has no effect on the calculation)
#endif
#if (defined(USING_MPI)) ///need to manage job directories too before we can implement checkpoint files well (although, storing lots of checkpoint files on each node is undesireable, so not a big priority perhaps)
            secondary_state_list_filename=scratch_directory+SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+"_cmp_t"+boost::lexical_cast<std::string>(world.rank())+".csv";
#else
            secondary_state_list_filename=job_path+SECONDARY_STATES_OUTPUT_DIRECTORY+SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+"_cmp.csv";
#endif
            compressed_secondary_state_list.open(secondary_state_list_filename.c_str());
            getline(compressed_secondary_state_list,line); //read the header line in the file (and ignore it)
            line_number=0;
            while (getline(compressed_secondary_state_list,line)){
#if (defined(USING_MPI))
                the_master_progress_bar.go(++line_number, CALCULATION_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
#else
                the_progress_bar.go(++line_number, CALCULATION_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
#endif
                ss << line;
                getline(ss,value,',');
                line_values.primary_E_internal=atof(value.c_str());
                getline(ss,value,',');
                line_values.primary_degeneracy=atof(value.c_str());
                getline(ss,value,',');
                line_values.primary_J=atoi(value.c_str());
                getline(ss,value,',');
                line_values.primary_vel=atof(value.c_str());
                ss.str("");
                ss.clear();
                secondary_result->total_parent_states_considered+=line_values.primary_degeneracy;
                portion_lines.push_back(line_values);
            }
        } //end if secondary states histogram or file
#if (defined(USING_MPI))
        the_master_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
        the_master_progress_bar.initialise(100,1,portion_lines.size(), CALCULATION_PROGRESS_BAR,3); // initialise the progress bar (this has no effect on the calculation)
#else
        the_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
        the_progress_bar.initialise(100,1,portion_lines.size(), CALCULATION_PROGRESS_BAR); // initialise the progress bar (this has no effect on the calculation)
#endif
        semipersistent_output << "Calculating Secondary Dissociation States Portions" << std::endl;
        print_semipersistent_output();
        for (i=0;i<portion_lines.size();i++){
#if (defined(USING_MPI))
            the_master_progress_bar.go(i, CALCULATION_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
#else
            the_progress_bar.go(i, CALCULATION_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
#endif
            dissociating_parent->initial_J=portion_lines[i].primary_J; //not sure I like this ... perhaps we need a clone for the input, we can't make the frag const if we do this!
            dissociating_parent->initial_velocity=portion_lines[i].primary_vel;
            secondary_PST_input->E_total=portion_lines[i].primary_E_internal;
            secondary_PST_input->E_available = secondary_PST_input->E_total-dissociating_parent->E_dissociation-dissociating_parent->E_barrier;
/*     ///**************************************************************************        angle distro  //// this bit can be moved outside because we are dumping the non n=0 wavefunction bit ... at least for now
            if (dissociating_parent->use_angle_distribution){
                secondary_PST_input->bond_angle_distribution.clear();
                for (w=bend_wf->start(primary_bend_n); w<bend_wf->end(primary_bend_n); w++){
                    temp_angle_distribution_element.radians=bend_wf->coord[w].rad;
                    temp_angle_distribution_element.probability=bend_wf->wavefunctions[primary_bend_n][w].prob;
                    secondary_PST_input->bond_angle_distribution.push_back(temp_angle_distribution_element);
                }
            }
  ///**************************************************************************        angle distro  */
            temp_secondary_result->reset();
            the_secondary_PST = new Triple_Frag_Secondary(secondary_PST_input,dissociating_parent,dissociating_frag1,dissociating_frag2,temp_secondary_result); //*********************** expensive to keep recreating objects inside a loop ... consider a reset function! ******************
            secondary_result->increment(portion_lines[i].primary_degeneracy, temp_secondary_result); //increment all the cumulative histograms
            delete the_secondary_PST; //************************************************************************* expensive to keep recreating objects inside a loop ... consider a reset function! **********************************
        }
#if (defined(USING_MPI))
        the_progress_bar.un_block();//release the global progress bar too because we blocked it when we initialised this MPI progress bar.
        the_master_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar (
        reduce(world, secondary_result, grand_secondary_result, Reduce_Secondary_Dissociation_Results(),0);
        persistent_output<<"Process "<<world.rank()<<" : grand_secondary_result->data->grand_total = "<<grand_secondary_result->data->grand_total<<std::endl;
        print_semipersistent_output();
#else
        (*grand_secondary_result)=(*secondary_result);
        the_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
#endif
        fprintf(out_file.fp, "\n\n\n\nSecondary Dissociation Histograms For Primary Fragment %s\n\n", frag_num_string.c_str());
        grand_secondary_result->write_histograms(out_file.fp);
        grand_secondary_result->write_bitmaps(bitmap_path_and_name_start);
        fprintf(out_file.fp, "\n\n\n\nTotal Primary Dissociation, %lf\n", primary_PST_output->data->grand_total);
        fprintf(out_file.fp, "Total F%s States Considered, %lf\n", frag_num_string.c_str(), grand_secondary_result->total_parent_states_considered); ///this needs reducing first I think!
        fprintf(out_file.fp, "Total F%s States Dissociated, %lf\n", frag_num_string.c_str(), grand_secondary_result->data->grand_total);
        fprintf(out_file.fp, "Percent F%s Dissociated, %lf\n", frag_num_string.c_str(), 100*grand_secondary_result->data->grand_total/primary_PST_output->data->grand_total);
        delete secondary_PST_input;
        delete secondary_result;
        delete grand_secondary_result;
        delete temp_secondary_result;
 //       if (dissociating_parent->use_angle_distribution) delete bend_wf;

#if (defined(USING_MPI))
        broadcast(world, job_path, 0);
        dest_path_and_name=job_path+SECONDARY_STATES_OUTPUT_DIRECTORY+SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+"_cmp_t"+boost::lexical_cast<std::string>(world.rank())+".csv";
        temp_bool=boost::filesystem::create_directories(job_path+SECONDARY_STATES_OUTPUT_DIRECTORY);
        try { //have to copy (instead of simple rename) ... because the file is (potentially) physically moving to a different HD
            boost::filesystem::copy_file( secondary_state_list_filename, dest_path_and_name );
        } catch (const boost::filesystem::filesystem_error& e) {
            std::cout<<"Error: " << e.what() << std::endl;
        }


    //if they selected to delete the secondary state files after the calculation has completed  ... do this
    //(default behaviour, big files that humans are unlikely to read - only really worth keeping as checkpoint files or for debugging)
        remove(secondary_state_list_filename.c_str());
#endif

    } else {
#if (defined(USING_MPI))
        sending_secondary_results=false;
        broadcast(world, sending_secondary_results, 0);
#endif
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Triple_Fragmentation_Calculation::run(){
#if (defined(USING_MPI)) //--------------------------------
    mpi::communicator world;
#endif //--------------------------------------------------
    Impulsive_Core_Input *primary_PST_input;
    Primary_Dissociation_Result *primary_PST_output, *temp_PST_output;
    Triple_Frag_Primary *the_primary_PST;
    Histogram_Limits sec_limits;
    Histogram_3d *sec_states1, *sec_states2;
    std::string bitmap_path_and_name_start;
    bool sending_primary_results(0);

    sec_states1 = new Histogram_3d();
    sec_states2 = new Histogram_3d();
    clear_semipersistent_output();
    persistent_output << " - Calculate Triple Fragmentation:" << std::endl;
    primary_PST_input = new Impulsive_Core_Input();
    temp_PST_output = new Primary_Dissociation_Result();
    primary_PST_output = new Primary_Dissociation_Result();
    bitmap_path_and_name_start=job_path+"3F_Primary_";
    //these three can probably be in a function set_standard_PST_input_values() or something
    primary_PST_input->job_path = job_path;
    primary_PST_input->output_to_file = t_general_input->objects_write_to_file;
    primary_PST_input->output_histograms = t_general_input->histogram_write_to_file;
    primary_PST_input->use_mpi=true; //without setting this expicitly in the input object, a core will default to a single thread count (doesn't start at threadnum and skip ahead by numthreads)
    //these will be manipulated differently by different calculations
    primary_PST_input->job_name="primary_fragmentation";
    primary_PST_input->E_total = E_total; //not sure we really need this any more ... looking more like the calculation object might do all IO ... consider removing as not used in core calculations only IO
    primary_PST_input->E_available = E_available_statistical;
    primary_PST_input->tunneling_model=t_general_input->tunneling_model; //later we may allow different models for each fragment and then we would use the general value if the fragment(parent) didn't have this value set (and we would add this variable to the impulsive molecule class too)
    primary_PST_input->impulsive_J_mode=t_general_input->impulsive_J_mode; //later we may allow different models for each fragment and then we would use the general value if the fragment(parent) didn't have this value set (and we would add this variable to the impulsive molecule class too)
    if (E_available_products>0){
        primary_PST_output->initialise(limits, i_reactant);
        temp_PST_output->initialise(limits, i_reactant);
#if (defined(USING_MPI)) //--------------------------------
        sending_primary_results=true;
        broadcast(world, sending_primary_results, 0);
        broadcast(world, primary_PST_input, 0);
        if (primary_PST_input->use_mpi) {
            broadcast(world, i_reactant, 0);
            broadcast(world, i_frag1, 0);
            broadcast(world, i_frag2, 0);
            broadcast(world, temp_PST_output, 0);
            broadcast(world, t_general_input->use_secondary_states_histograms, 0);
        }
#endif //--------------------------------------------------
        semipersistent_output << "Counting Primary Dissociation States" << std::endl;
        print_semipersistent_output();
        if (t_general_input->use_secondary_states_histograms) {
            sec_limits.E_max1 = E_total-(i_reactant->E_dissociation+i_frag1->E_dissociation);
            sec_limits.E_max2 = E_total-(i_reactant->E_dissociation+i_frag2->E_dissociation);
            sec_limits.vel_max1 = sqrt(sec_limits.E_max1*i_frag1->get_vel_mass_factor(i_frag2->mass));
            sec_limits.vel_max2 = sqrt(sec_limits.E_max2*i_frag2->get_vel_mass_factor(i_frag1->mass));
            sec_limits.J_max1 = i_frag1->get_n_max(E_available_products);
            sec_limits.J_max2 = i_frag2->get_n_max(E_available_products);
            //we don't need v_max1 and 2 here as sec_limits is just for the sec_states1 and 2 3d histograms, which are for storing the secondary PST calculation values, besides, we couldn't call get_v_max here anyway as the frag->initialise_vib_levels routine hasn't been called for these fragments
            //don't initialise the secondary states histogram unless we are going to use it, otherwise we waste potentially quite a lot of ram
            if (i_frag1->E_dissociation) {
                sec_states1->initialiseN("Secondary States Fragment 1","Count","Internal Energy","Parent Velocity","Parent J",HISTOGRAM_SIZE,2*HISTOGRAM_SIZE,(sec_limits.J_max1+1),sec_limits.E_max1,sec_limits.vel_max1,sec_limits.J_max1);
#if (defined(USING_MPI))
                broadcast(world, sec_states1, 0);
#endif
            }else if (i_frag2->E_dissociation) {
                sec_states2->initialiseN("Secondary States Fragment 2","Count","Internal Energy","Parent Velocity","Parent J",HISTOGRAM_SIZE,2*HISTOGRAM_SIZE,(sec_limits.J_max2+1),sec_limits.E_max2,sec_limits.vel_max2,sec_limits.J_max2);
#if (defined(USING_MPI))
                broadcast(world, sec_states2, 0);
#endif
            }
            the_primary_PST= new Triple_Frag_Primary(primary_PST_input,i_reactant,i_frag1,i_frag2,temp_PST_output,sec_states1,sec_states2);
        }else{ //Triple_Frag_Primary constructor without sec_states histogram arguments knows to write each secondary state to the file instead
            the_primary_PST= new Triple_Frag_Primary(primary_PST_input,i_reactant,i_frag1,i_frag2,temp_PST_output);
            temp_PST_output->store_master_thread_counts();
        }
        delete the_primary_PST;
#if (defined(USING_MPI))
        if (primary_PST_input->use_mpi){
            if (t_general_input->use_secondary_states_histograms) {
                if (i_frag1->E_dissociation) all_reduce(world, sec_states1, Reduce_Histogram_3d());
                else if (i_frag2->E_dissociation) all_reduce(world, sec_states2, Reduce_Histogram_3d());
            }
            reduce(world, temp_PST_output, primary_PST_output, Reduce_Primary_Results(),0);
        }else if (!world.rank()) (*primary_PST_output)=(*temp_PST_output);
#else
        (*primary_PST_output)=(*temp_PST_output);
#endif
        std::cout<<"first count complete"<<std::endl;
        fprintf(out_file.fp, "\n\n\n\nPrimary Dissociation Histograms\n\n");

        primary_PST_output->write_histograms(out_file.fp);
        primary_PST_output->write_bitmaps(bitmap_path_and_name_start);
        std::cout<<"primary histograms written"<<std::endl;
        //dissociating fragment can be defined in the input file either way round (and eventually both frags could possibly fragment [replace else if with if], but will have to be careful with degeneracy! degen/2?)
    /// IF WE FOUND ANY STATES TO CHECK -- we should know this here, which means we'd avoid the check inside???
            if (i_frag1->E_dissociation) secondary_fragmentation(primary_PST_output,i_frag1,frag1_S1,frag1_S2,sec_states1); //a secondary dissociation has been defined for fragment 1, calculate it
            else if (i_frag2->E_dissociation) secondary_fragmentation(primary_PST_output,i_frag2,frag2_S1,frag2_S2,sec_states2); //a secondary dissociation has been defined for fragment 2, calculate it
        ///for now only one fragment can be flagged for secondary dissociation in the file (by setting a dissociation energy) ... I doubt the output would be correct otherwise!!! (degeneracy)
    }else{
        persistent_output<<"Error: Insufficient energy, check input file!"<<std::endl;
        print_persistent_output();
#if (defined(USING_MPI)) //--------------------------------
        sending_primary_results=false;
        broadcast(world, sending_primary_results, 0);
#endif //--------------------------------------------------
    }
    delete sec_states1;
    delete sec_states2;
    delete primary_PST_input;
    delete primary_PST_output;
    delete temp_PST_output;
    return 1;
}
//***********************************************************************************************************************************************************************************
#if (defined(USING_MPI))
void Triple_Fragmentation_Calculation_Slave_Process(){
    mpi::communicator world;
    std::vector<Secondary_States_Line> portion_lines;
    Impulsive_Molecule *i_reactant, *i_frag1, *i_frag2;
    Impulsive_Phasespace_Core *the_primary_PST;
    Impulsive_Core_Input *primary_PST_input;
    Primary_Dissociation_Result *primary_PST_output;
    Triple_Fragmentation_Core_Input *secondary_PST_input;
    Impulsive_Molecule *dissociating_parent, *dissociating_frag1, *dissociating_frag2;
    Secondary_Dissociation_Result *secondary_result, *temp_secondary_result;
    Triple_Frag_Secondary *the_secondary_PST;
    Histogram_3d *sec_states1, *sec_states2, *sec_states;
    Output_File compressed_secondary_states;
    std::ifstream secondary_state_list, compressed_secondary_state_list;
    std::string line, last_line(""), value, secondary_state_list_filename, dest_path_and_name, job_path;
    Secondary_States_Line line_values, last_line_values;
    std::stringstream ss;
    bool not_first_line(0), temp_bool;
    double final_primary_degeneracy(0);
    bool sending_primary_results(0), sending_secondary_results(0), sending_sec_states(0);
    unsigned line_number(0), i,j,k, progress_bar_length(1), num_lines;
    dMPI_Slave_Progress_Bar the_slave_progress_bar;

    sec_states1 = new Histogram_3d();
    sec_states2 = new Histogram_3d();
    broadcast(world, sending_primary_results, 0);
    if (sending_primary_results){
        broadcast(world, primary_PST_input, 0);
        if (primary_PST_input->use_mpi) {
            broadcast(world, i_reactant, 0);
            broadcast(world, i_frag1, 0);
            broadcast(world, i_frag2, 0);
            broadcast(world, primary_PST_output, 0);
            broadcast(world, sending_sec_states, 0);
            if(sending_sec_states){
                if (i_frag1->E_dissociation) broadcast(world, sec_states1, 0);
                else if (i_frag2->E_dissociation) broadcast(world, sec_states2, 0);
                the_primary_PST = new Triple_Frag_Primary(primary_PST_input,i_reactant,i_frag1,i_frag2,primary_PST_output,sec_states1,sec_states2);
                if (i_frag1->E_dissociation) all_reduce(world, sec_states1, Reduce_Histogram_3d());
                else if (i_frag2->E_dissociation) all_reduce(world, sec_states2, Reduce_Histogram_3d());
                if (i_frag1->E_dissociation) sec_states=sec_states1;
                else if (i_frag2->E_dissociation) sec_states=sec_states2; //this feels a little hacky ... we'll need the secondary part to be a function like it is for the master thread I think ...
            }else{
                the_primary_PST = new Triple_Frag_Primary(primary_PST_input,i_reactant,i_frag1,i_frag2,primary_PST_output);
            }
            delete the_primary_PST;
            reduce(world, primary_PST_output, Reduce_Primary_Results(),0);
            //std::cout<<"Process "<<world.rank()<<" : primary reduced"<<std::endl;
            delete i_reactant;
            delete i_frag1;
            delete i_frag2;
        }//may want else here to broadcast data (primary_PST_output) required for the following steps?
        delete primary_PST_input;

        broadcast(world, sending_secondary_results, 0);
        if (sending_secondary_results){
            broadcast(world, secondary_PST_input, 0);
            broadcast(world, dissociating_parent, 0);
            broadcast(world, dissociating_frag1, 0);
            broadcast(world, dissociating_frag2, 0);
            broadcast(world, secondary_result, 0);
            if(sending_sec_states){
                the_slave_progress_bar.initialise(100,1,sec_states->bin_length1, CALCULATION_PROGRESS_BAR, 1);
                for (i=world.rank(); i<sec_states->bin_length1; i+=world.size()){
                    the_slave_progress_bar.go(i, CALCULATION_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
                    for (j=0; j<sec_states->bin_length2; j++){
                        for (k=0; k<sec_states->bin_length3; k++){
                            if (sec_states->val(i,j,k)){
                                line_values.primary_degeneracy=sec_states->val(i,j,k);
                                line_values.primary_E_internal=sec_states->resolveIndex1(i)+dissociating_parent->E_dissociation; //add parent E_dissoc because we store EXCESS energy in the 3d histogram not internal energy like in the file
                                line_values.primary_vel=sec_states->resolveIndex2(j);
                                line_values.primary_J=sec_states->resolveIndex3(k); //not sure I like this ... perhaps we need a clone for the input, we can't make the frag const if we do this!
                                secondary_result->total_parent_states_considered+=line_values.primary_degeneracy;
                                portion_lines.push_back(line_values);
                            }
                        }
                    }
                }
            }else{ //states to calculate are stored in a file locally //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                secondary_state_list_filename=scratch_directory+SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+"_t"+boost::lexical_cast<std::string>(world.rank())+".csv";
                secondary_state_list.open(secondary_state_list_filename.c_str());///consider using binary files for this too - faster for large files etc
                compressed_secondary_states.set_output_directory(scratch_directory);
                compressed_secondary_states.open(SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+"_cmp_t"+boost::lexical_cast<std::string>(world.rank())+".csv");
                getline(secondary_state_list,line); //read the header line in the file
                fprintf(compressed_secondary_states.fp, "%s\n", line.c_str()); //write the header
                num_lines=line_number=0;
                if (dissociating_parent->number == 1) progress_bar_length=primary_PST_output->f1_num_secondary;
                else if (dissociating_parent->number == 2) progress_bar_length=primary_PST_output->f2_num_secondary;
                the_slave_progress_bar.initialise(100,1,progress_bar_length, CALCULATION_PROGRESS_BAR, 1); // initialise the progress bar (this has no effect on the calculation)
                while (getline(secondary_state_list,line)){
                    the_slave_progress_bar.go(++line_number, CALCULATION_PROGRESS_BAR);
                    ss << line;
                    getline(ss,value,',');
                    line_values.primary_E_internal=atof(value.c_str());
                    getline(ss,value,',');
                    line_values.primary_degeneracy=atof(value.c_str());
                    getline(ss,value,',');
                    line_values.primary_J=atoi(value.c_str());
                    getline(ss,value,',');
                    line_values.primary_vel=atof(value.c_str());
                    ss.str("");
                    ss.clear();
                    if (line==last_line) final_primary_degeneracy+=line_values.primary_degeneracy;
                    else{
                        if (not_first_line) fprintf(compressed_secondary_states.fp, "%lf,%lf,%d,%lf\n", last_line_values.primary_E_internal, final_primary_degeneracy, last_line_values.primary_J, last_line_values.primary_vel);
                        not_first_line=true;
                        num_lines++;
                        last_line=line;
                        last_line_values=line_values;
                        final_primary_degeneracy=line_values.primary_degeneracy;
                    }
                }
                if (not_first_line) {
                    fprintf(compressed_secondary_states.fp, "%lf,%lf,%d,%lf\n", last_line_values.primary_E_internal, final_primary_degeneracy, last_line_values.primary_J, last_line_values.primary_vel);
                    num_lines++;
                }
                compressed_secondary_states.close();
                remove(secondary_state_list_filename.c_str()); //delete the uncompressed file
                the_slave_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
                the_slave_progress_bar.initialise(100,1,num_lines, CALCULATION_PROGRESS_BAR, 2);
                secondary_state_list_filename=scratch_directory+SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+"_cmp_t"+boost::lexical_cast<std::string>(world.rank())+".csv";
                compressed_secondary_state_list.open(secondary_state_list_filename.c_str());
                getline(compressed_secondary_state_list,line); //read the header line in the file (and ignore it)
                line_number=0;
                while (getline(compressed_secondary_state_list,line)){
                    the_slave_progress_bar.go(++line_number, CALCULATION_PROGRESS_BAR);
                    ss << line;
                    getline(ss,value,',');
                    line_values.primary_E_internal=atof(value.c_str());
                    getline(ss,value,',');
                    line_values.primary_degeneracy=atof(value.c_str());
                    getline(ss,value,',');
                    line_values.primary_J=atoi(value.c_str());
                    getline(ss,value,',');
                    line_values.primary_vel=atof(value.c_str());
                    ss.str("");
                    ss.clear();
                    secondary_result->total_parent_states_considered+=line_values.primary_degeneracy;
                    portion_lines.push_back(line_values);
                }
            }//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            temp_secondary_result = new Secondary_Dissociation_Result();
            (*temp_secondary_result)=(*secondary_result);
            the_slave_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
            the_slave_progress_bar.initialise(100,1,portion_lines.size(), CALCULATION_PROGRESS_BAR, 3);
            for (i=0;i<portion_lines.size();i++){
                the_slave_progress_bar.go(i, CALCULATION_PROGRESS_BAR);
                //std::cout<<"Process "<<world.rank()<<" : the_slave_progress_bar.go("<<i<<", CALCULATION_PROGRESS_BAR);"<<std::endl;
                dissociating_parent->initial_J=portion_lines[i].primary_J; //not sure I like this ... perhaps we need a clone for the input, we can't make the frag const if we do this!
                dissociating_parent->initial_velocity=portion_lines[i].primary_vel;
                secondary_PST_input->E_total=portion_lines[i].primary_E_internal;
                secondary_PST_input->E_available=secondary_PST_input->E_total-dissociating_parent->E_dissociation-dissociating_parent->E_barrier;
                temp_secondary_result->reset();
                the_secondary_PST = new Triple_Frag_Secondary(secondary_PST_input,dissociating_parent,dissociating_frag1,dissociating_frag2,temp_secondary_result); //*********************** expensive to keep recreating objects inside a loop ... consider a reset function! ******************
                secondary_result->increment(portion_lines[i].primary_degeneracy, temp_secondary_result); //increment all the cumulative histograms
                delete the_secondary_PST; //************************************************************************* expensive to keep recreating objects inside a loop ... consider a reset function! **********************************
            }
            //std::cout<<"Process "<<world.rank()<<" : secondary slave counted "<<secondary_result->data->grand_total<<std::endl;
            the_slave_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
            reduce(world, secondary_result, Reduce_Secondary_Dissociation_Results(),0);
            //std::cout<<"Process "<<world.rank()<<" : secondary reduced"<<std::endl;

            broadcast(world, job_path, 0);
            dest_path_and_name=job_path+SECONDARY_STATES_OUTPUT_DIRECTORY+SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(dissociating_parent->number)+"_cmp_t"+boost::lexical_cast<std::string>(world.rank())+".csv";
            temp_bool=boost::filesystem::create_directories(job_path+SECONDARY_STATES_OUTPUT_DIRECTORY);
            try { //have to copy (instead of simple rename) ... because the file is (potentially) physically moving to a different HD
                boost::filesystem::copy_file( secondary_state_list_filename, dest_path_and_name );
            } catch (const boost::filesystem::filesystem_error& e) {
                std::cout<<"Error: " << e.what() << std::endl;
            }


        //if they selected to delete the secondary state files after the calculation has completed  ... do this
        //(default behaviour, big files that humans are unlikely to read - only really worth keeping as checkpoint files or for debugging)
            remove(secondary_state_list_filename.c_str());

            delete secondary_PST_input;
            delete dissociating_parent;
            delete dissociating_frag1;
            delete dissociating_frag2;
            delete secondary_result;
            delete temp_secondary_result;
        }
        delete sec_states1;
        delete sec_states2;
        delete primary_PST_output;
    }
}
//***********************************************************************************************************************************************************************************
#endif

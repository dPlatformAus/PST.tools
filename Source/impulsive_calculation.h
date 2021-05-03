class Impulsive_Calculation : public Phasespace_Calculation/**
    There is a lot more work required before general impulsive calculations will work:
        two polyatomic fragments is not yet supported (both centres of mass, moments of inertia and division of impulse rotation etc)!
        the leaving fragment is specific to H2CO, as defined in our input
        obviously the run routine itself!
        and VERY IMPORTANTLY, the process state and bin state stuff in the core class needs to be updated to at least part of what the triple version does!
*/{
    public:
        Impulsive_General_Input *i_general_input;
        Impulsive_Molecule *i_reactant, *i_frag1, *i_frag2;
        double E_available_products, E_max_products;
        Impulsive_Calculation(){}; //C++ doesn't do virtual constructors, so we can't do much with this
        virtual ~Impulsive_Calculation(){}; //no new objects ... they should be deleted by ~Phasespace_Calculation()
        virtual void create_general_input();
        virtual void create_general_input_pointer();
        virtual void create_fragments();
        virtual void create_fragment_pointers();
        virtual void set_energies();
        virtual void set_limits();
        virtual void initialise(const std::string &file_name);
        virtual bool run();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Calculation::create_general_input(){
    the_general_input = new Impulsive_General_Input();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Calculation::create_general_input_pointer(){
    i_general_input = dynamic_cast<Impulsive_General_Input*>(the_general_input);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Calculation::create_fragments(){
    reactant = new Impulsive_Molecule();
    reactant->construct();
    frag1 = new Impulsive_Molecule();
    frag1->construct();
    frag2 = new Impulsive_Molecule();
    frag2->construct();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Calculation::create_fragment_pointers(){
    i_reactant = dynamic_cast<Impulsive_Molecule*>(reactant);
    i_frag1 = dynamic_cast<Impulsive_Molecule*>(frag1);
    i_frag2 = dynamic_cast<Impulsive_Molecule*>(frag2);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Calculation::set_energies(){
    E_total=the_general_input->excitation_energy_range.max;
    E_available_products=E_total-i_reactant->E_dissociation;
    E_available_statistical=E_available_products-i_reactant->E_barrier;
    E_max_statistical=E_available_statistical+parent_Erovib_max; //add the maximum parent internal energy so that adding that energy doesn' put the answers outside histogram limits
    E_max_products=E_available_products+parent_Erovib_max; //add the maximum parent internal energy so that adding that energy doesn' put the answers outside histogram limits
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Calculation::set_limits(){
    unsigned i;
    limits.E_max = E_max_products;
    limits.vel_max1 = sqrt(E_max_products*i_frag1->get_vel_mass_factor(i_frag2->mass));
    limits.vel_max2 = sqrt(E_max_products*i_frag2->get_vel_mass_factor(i_frag1->mass));
    limits.E_max1 = i_frag1->get_E_trans(limits.vel_max1);
    limits.E_max2 = i_frag2->get_E_trans(limits.vel_max2);
    limits.J_max1 = i_frag1->get_n_max(E_max_products); //max is from impulse and statistical, so use product max avail energy
    limits.J_max2 = i_frag2->get_n_max(E_max_products);
    limits.v_max1 = i_frag1->get_v_max(E_max_statistical); //no vibrational excitation from impulse yet in our model, so max is statistical max
    limits.v_max2 = i_frag2->get_v_max(E_max_statistical);
    for (i=0; i<=limits.v_max1; i++) limits.frag1_vib_names.push_back(i_frag1->vib_levels.bands[i].name);
    for (i=0; i<=limits.v_max2; i++) limits.frag2_vib_names.push_back(i_frag2->vib_levels.bands[i].name);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Calculation::initialise(const std::string &file_name){
    Phasespace_Calculation::initialise(file_name);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Impulsive_Calculation::run(){
#if (defined(USING_MPI))
    mpi::communicator world;
#endif
    Impulsive_Phasespace_Core *the_PST;
    Impulsive_Core_Input *the_PST_input;
    Phasespace_Result *temp_PST_output, *the_PST_output;
    Impulsive_Bond_Angle_Distribution_Element temp_angle_distribution_element;
    unsigned i, primary_bend_n;
    bool sending_results;

    std::string bitmap_path_and_name_start;
    //std::string wf_file_name;
    //HO_Wavefunctions *bend_wf;

    persistent_output << " - Impulsive Calculation:" << std::endl;
    the_PST_input = new Impulsive_Core_Input();
    the_PST_output = new Phasespace_Result();
    temp_PST_output = new Phasespace_Result();
    bitmap_path_and_name_start=job_path+"2F_";
    //these three can probably be in a function set_standard_PST_input_values() or something
    the_PST_input->job_path = job_path;
    the_PST_input->output_to_file = the_general_input->objects_write_to_file;
    the_PST_input->output_histograms = the_general_input->histogram_write_to_file;
    the_PST_input->use_mpi=true; //without setting this expicitly in the input file, a core will default to a single thread count (doesn't start at threadnum and skip ahead by numthreads)
    //these will be manipulated differently by different calculations
    the_PST_input->tunneling_model=i_general_input->tunneling_model; //later we may allow different models for each fragment and then we would use the general value if the fragment(parent) didn't have this value set (and we would add this variable to the impulsive molecule class too)
    the_PST_input->impulsive_J_mode=i_general_input->impulsive_J_mode; //later we may allow different models for each fragment and then we would use the general value if the fragment(parent) didn't have this value set (and we would add this variable to the impulsive molecule class too)
    the_PST_input->E_total = E_total; //not sure we really need this any more ... looking more like the calculation object might do all IO ... consider removing as not used in core calculations only IO
    the_PST_input->E_available = E_available_statistical;


    if (E_available_products>0){
        the_PST_input->job_name = "Impulsive_PST_"+boost::lexical_cast<std::string>(the_PST_input->E_total);
     //**************************************************************************        angle distro
        if (i_reactant->use_angle_distribution){
    /*        primary_bend_n=0; //if not 0 this would have to be input from file ... which has not been implemented and is unlikely to be either
            wf_file_name="bend_wf.csv";
            bend_wf = new HO_Wavefunctions(wf_file_name, job_path, i_general_input, UNITS_GAUSSIAN, primary_bend_n, i_reactant); //get the wavefunctions for the maximum n value for the bend
            the_PST_input->bond_angle_distribution.clear();
            for (i=bend_wf->start(primary_bend_n); i<bend_wf->end(primary_bend_n); i++){
                temp_angle_distribution_element.radians=bend_wf->coord[i].rad;
                temp_angle_distribution_element.probability=bend_wf->wavefunctions[primary_bend_n][i].prob;
                the_PST_input->bond_angle_distribution.push_back(temp_angle_distribution_element);
            }*/
        }
     //**************************************************************************        angle distro
        the_PST_output->initialise(limits, i_reactant);
        temp_PST_output->initialise(limits, i_reactant);

#if (defined(USING_MPI))
        sending_results=true;
        broadcast(world, sending_results, 0);
        broadcast(world, the_PST_input, 0);
        if (the_PST_input->use_mpi) {
            broadcast(world, i_reactant, 0);
            broadcast(world, i_frag1, 0);
            broadcast(world, i_frag2, 0);
            broadcast(world, temp_PST_output, 0);
        }
#endif
        the_PST = new Impulsive_Phasespace_Core(the_PST_input, i_reactant, i_frag1, i_frag2, temp_PST_output);

#if (defined(USING_MPI))
        //want to see if we can't make this binary function just work with pointers? Think we can!
        //polymorphism, stack size ... etc!
        if (the_PST_input->use_mpi) reduce(world, temp_PST_output, the_PST_output, Reduce_Phasespace_Results(),0);
        else if (!world.rank()) (*the_PST_output)=(*temp_PST_output);
#else
        (*the_PST_output)=(*temp_PST_output);
#endif
        the_PST_output->write_histograms(out_file.fp);
        the_PST_output->write_bitmaps(bitmap_path_and_name_start);
        delete the_PST;
    }else{
        persistent_output<<"Error: Insufficient energy, check input file!"<<std::endl;
        print_persistent_output();
#if (defined(USING_MPI))
        mpi::communicator world;
        sending_results=false;
        broadcast(world, sending_results, 0);
#endif
    }
    delete the_PST_input;
    delete the_PST_output;
    delete temp_PST_output;
    return 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
#if (defined(USING_MPI))
void Impulsive_Calculation_Slave_Process(){
    mpi::communicator world;
    Impulsive_Molecule *i_reactant, *i_frag1, *i_frag2;
    Impulsive_Phasespace_Core *the_PST;
    Impulsive_Core_Input *the_PST_input;
    Phasespace_Result *the_PST_output;
    bool sending_results;
    unsigned i;
    broadcast(world, sending_results, 0);
    if (sending_results){
        the_PST_input = new Impulsive_Core_Input();
        broadcast(world, the_PST_input, 0);
        if (the_PST_input->use_mpi) {
            broadcast(world, i_reactant, 0);
            broadcast(world, i_frag1, 0);
            broadcast(world, i_frag2, 0);
            broadcast(world, the_PST_output, 0);
            the_PST = new Impulsive_Phasespace_Core(the_PST_input,i_reactant,i_frag1,i_frag2,the_PST_output); //*********************** expensive to keep recreating objects inside a loop ... consider a reset function! ******************
            delete the_PST;//*********************** expensive to keep recreating objects inside a loop ... consider a reset function! (before, not after:)******************
            reduce(world, the_PST_output, Reduce_Phasespace_Results(),0);
            std::cout<<"Process "<<world.rank()<<" : reduced"<<std::endl;
            delete i_reactant;
            delete i_frag1;
            delete i_frag2;
            delete the_PST_output;
        }
        delete the_PST_input;
    }
}
//***********************************************************************************************************************************************************************************
#endif
//***********************************************************************************************************************************************************************************



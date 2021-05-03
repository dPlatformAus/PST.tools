//***********************************************************************************************************************************************************************************
class Phasespace_Calculation : public Base_Calculation/**
The broad strokes of this program are that it reads input files from a directory and executes a certain type of Phasespace Theory (PST) calculation, depending on the contents of the file.
When each file is read, there is a line indicating what type of calculation is to be run.
This program then creates a calculation object, corresponding to the calculation type requested in the file, and passes the file to the calculation object.
The calculation object then takes care of all the remaining processing, from parsing all other input from the file, doing the actual calculations and writing the output files.
In the future, a factory pattern may be used to create the appropriate calculation object.
But currently, this is simply achieved with a big switch statement (well I should change it to a switch, rather than else ifs, at the VERY least).

This is the base calculation class. All other calculations will be derived from this or one of it's deriatives.
It has members that are General_Input and Phasespace_Molecule type objects, which, for the purposes of this class, parse the required input from the input file and store that data.
This class facilitates running a calculation using those objects and a Phasespace_Core object that this class will create in it's run() function and pass parameters to (including the molecule objects).
The main output of the calculation will be written to an Output_File that is also a member of the calculation.

There will be many sub class versions of this class (and it's derivatives), which will override one or more functions;
At the very least the run() function will be overloaded, as this is where Phasespace_Core objects (or derivatives of the Phasespace_Core class) are created
to perform PST counts and the results of these counts are retrieved for use in further calculations and/or output.
initialise(const string &file_name) may also be overridden to allow additional initialisation, specific to the task of the derived calculation.
However, the overridden version of initialise(const string &file_name) must always call the base implementation first!
The base implementation of initialise(const string &file_name) is where the functions that create the General_Input and Phasespace_Molecule type objects are called, after which their parse routines are initiated.
create_objects() and create_general_input() are overridden in derived versions of this class when a more derived General_Input and/or Phasespace_Molecule type object is required by the calculation.

create_objects() and create_general_input() may, in the future, be changed to use an "abstract factory pattern" to create the General_Input and Phasespace_Molecule type objects.
However, I dont think that will change the special requirements outlined in the following paragraph because this base class will still be defined with pointers to the base General_Input and Phasespace_Molecule type objects.

Because the the_general_input, reactant, frag1 and frag2 members are pointers to their base classes, if they end up pointing to objects of more derived versions of their respective class types,
the overridden versions of initialise(const string &file_name) and run() will need to dynamically cast them down to their derived class type to be able to access the functionality added to the derived versions.
This dynamic_cast should always be safe because those member objects will have been created by the derived version of this class and so will always be of the correct type expected by that same derived version of this class.
pointers to fulfill this requirement have been added as members of all the derived versions of this and the core PST classes.

Derivatives of the Phasespace_Core class, which also inherit members that are pointers to the base Phasespace_Molecule class,
will also have to dynamically cast those pointers when their methods use features specific to the Phasespace_Molecule derivative.
This should also always be safe because such Phasespace_Core derivatives will only be created by Phasespace_Calculation derivatives that create that same Phasespace_Molecule derivative.
The Phasespace_Molecule derivatives created in the Phasespace_Calculation derivative are passed to the Phasespace_Core derivative.
Therefore the dynamic cast will always be on a pointer to the correct Phasespace_Molecule derivative.

This design has been chosen becuase we do not want the interfaces of all the (potentially many) derived classes, designed for specific niche calculations, percolating up to the base classes.
We want to be able to present a subset of our code that performs one niche type of PST calculation, without revealing all of the interfaces required for other niche calculations that may be completely unrelated.
At the same time, we do not want to duplicate code that will be common to all, or a subset of the calculation types etc. We recognise that dynamic_cast is potentially an indication of poor object oriented design.
However, our concern is that the complexity required to achieve the above goals without using dynamic_cast (carefully) will make this code less accessible to people who are not familiar with object oriented design patterns.
*/{
    public:
        General_Input *the_general_input;
        Phasespace_Molecule *reactant, *frag1, *frag2;
        double E_total, E_available_statistical, E_max_statistical;
        Histogram_Limits limits;
        //these are all for the parent temps
        double parent_Erovib_max;
        unsigned jvE_i_max;
        jvE_histogram *Erovib_histogram;
        //-----

        Phasespace_Calculation(){}; //C++ doesn't do virtual constructors, so we can't do much with this
        virtual ~Phasespace_Calculation();
        virtual void create_objects();
        virtual void create_general_input();
        virtual void create_general_input_pointer(){};
        virtual void create_fragments();
        virtual void create_fragment_pointers(){};
        virtual void set_energies();
        virtual void set_limits();
        virtual void jvE_initialise();
        virtual void initialise(const std::string &file_name);
        virtual bool run(); ///creates a Phasespace_Core object passes it some parameters, calls it's run() method and then retrieves some results and writes them to the Output_File
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Calculation::~Phasespace_Calculation(){
    delete the_general_input;
    delete reactant;
    delete frag1;
    delete frag2;
    delete Erovib_histogram;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Calculation::create_general_input(){
    the_general_input = new General_Input();
    persistent_output<<"create General_Input "<<std::endl;
    print_persistent_output();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Calculation::create_fragments(){
    reactant = new Phasespace_Molecule();
    reactant->construct();
    frag1 = new Phasespace_Molecule();
    frag1->construct();
    frag2 = new Phasespace_Molecule();
    frag2->construct();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Calculation::create_objects(){
    create_general_input();
    create_fragments();
    create_general_input_pointer();
    create_fragment_pointers();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Calculation::set_energies(){
    E_total=the_general_input->excitation_energy_range.max; //use max because either only one value specified, in which case min and max are same, else we want to initialise limits to the largest energy encountered in the calculation
    E_available_statistical=E_total-reactant->E_dissociation;
    E_max_statistical=E_available_statistical+parent_Erovib_max; //add the maximum parent internal energy so that adding that energy doesn' put the answers outside histogram limits
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Calculation::set_limits(){
    unsigned i;
    limits.E_max = E_max_statistical;
    limits.vel_max1 = sqrt(E_max_statistical*frag1->get_vel_mass_factor(frag2->mass));
    limits.vel_max2 = sqrt(E_max_statistical*frag2->get_vel_mass_factor(frag1->mass));
    limits.E_max1 = frag1->get_E_trans(limits.vel_max1);
    limits.E_max2 = frag2->get_E_trans(limits.vel_max2);
    limits.J_max1 = frag1->get_n_max(E_max_statistical);
    limits.J_max2 = frag2->get_n_max(E_max_statistical);
    limits.v_max1 = frag1->get_v_max(E_max_statistical);
    limits.v_max2 = frag2->get_v_max(E_max_statistical);
    for (i=0; i<=limits.v_max1; i++) limits.frag1_vib_names.push_back(frag1->vib_levels.bands[i].name);
    for (i=0; i<=limits.v_max2; i++) limits.frag2_vib_names.push_back(frag2->vib_levels.bands[i].name);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Calculation::jvE_initialise(){
    /*
    This is at the base PST calc level because eventually we want the parent internal energy options to be available for all calculations.
    However, this will require some reasonably detailed changes and I'm in write up mode now,
    documentation for the code has already started and these changes would require significant changes to it too, you have to draw a line somewhere.
    At the moment, the parent temp options are only utilised in the Roaming_Calculation derivatives.
    Roaming_Calculation itself is incomplete and is where I intended to implement it first ... although this class, Phasespace_Calculation, might be a more sensible place to start.
    There are a few more notes on my initial thoughts about implementing this in the Roaming_Calculation::run function ...
    Roaming_Calculation is intended for PST calcs at a single energy for the purpose of generating histogram output for product state distributions.
    There is even more work required (sudden approximation for moment of inertia change etc for roational distros etc) before the roaming histograms are really worthwhile.
    */
    std::string jvE_file_name;
    double current_prob_included;
    parent_Erovib_max=0;
    jvE_file_name="jvE_histogram.csv";
    Erovib_histogram=new jvE_histogram(jvE_file_name, job_path, the_general_input, reactant);
    current_prob_included=jvE_i_max=0;//calculate how far down the jvE_histogram to go to get to the_general_input->Boltzmann_probability_included, doing this first means we can define the progress bar better and also means the inner loop for each fit point is slightly more efficient
    while (jvE_i_max<Erovib_histogram->size() && (current_prob_included*100<=the_general_input->Boltzmann_probability_included || the_general_input->Boltzmann_probability_included==100)) {
        current_prob_included+=Erovib_histogram->jvE_hist[jvE_i_max].jvE_prob;
        if (parent_Erovib_max<Erovib_histogram->jvE_hist[jvE_i_max].energy) parent_Erovib_max=Erovib_histogram->jvE_hist[jvE_i_max].energy;
        jvE_i_max++;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Calculation::initialise(const std::string &file_name){
    Base_Calculation::initialise(file_name);
    create_objects();
    reactant->number=0;
    frag1->number=1;
    frag2->number=2;
    //set histogram defaults before parse sets any user overridden settings
    reactant->is_on->E_trans=true;
    reactant->is_on->frag1_vel=true;
    reactant->is_on->frag2_vel=true;
    reactant->probe_is_on=reactant->is_on->clone(); //same defaults for probes
    frag1->is_on=reactant->is_on->clone(); //same defaults for secondary histograms
    frag1->probe_is_on=reactant->is_on->clone(); //same defaults for secondary histograms
    frag2->is_on=reactant->is_on->clone(); //same defaults for secondary histograms
    frag2->probe_is_on=reactant->is_on->clone(); //same defaults for secondary histograms
    the_general_input->parse_input(file_name);
    reactant->parse_input(file_name, "REACTANT");
    frag1->parse_input(file_name, "FRAGMENT 1");
    frag2->parse_input(file_name, "FRAGMENT 2");
    frag1->initialise(frag2->mass);
    frag2->initialise(frag1->mass); //reactant->initialise(); is unnessesary, pst cores dont use bbar prolate or spin orbit coupling state energies of the parent
    jvE_initialise();
    set_energies();
    frag1->initialise_vib_levels(E_max_statistical);
    frag2->initialise_vib_levels(E_max_statistical);
    set_limits();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Phasespace_Calculation::run(){
#if (defined(USING_MPI))
    mpi::communicator world;
#endif
    Phasespace_Core *the_PST;
    Phasespace_Core_Input *the_PST_input;
    Phasespace_Result *temp_PST_output, *the_PST_output;
    std::string bitmap_path_and_name_start;
    bool sending_results;

    persistent_output << " - Phasespace_Calculation:" << std::endl;
    the_PST_input = new Phasespace_Core_Input();
    temp_PST_output = new Phasespace_Result();
    the_PST_output = new Phasespace_Result();

    bitmap_path_and_name_start=job_path+"2F_";
    //these three can probably be in a function set_standard_PST_input_values() or something
    the_PST_input->job_path = job_path;
    the_PST_input->output_to_file = the_general_input->objects_write_to_file;
    the_PST_input->output_histograms = the_general_input->histogram_write_to_file;
    the_PST_input->use_mpi=true; //without setting this expicitly in the input file, a core will default to a single thread count (doesn't start at threadnum and skip ahead by numthreads)

    if (E_available_statistical>0){
        the_PST_input->E_total = E_total; //not sure we really need this any more ... looking more like the calculation object might do all IO ... consider removing as not used in core calculations only IO
        the_PST_input->E_available = E_available_statistical;
        the_PST_input->job_name = "PST_"+boost::lexical_cast<std::string>(E_total);

        the_PST_output->initialise(limits, reactant);
        temp_PST_output->initialise(limits, reactant);

#if (defined(USING_MPI))
        sending_results=true;
        broadcast(world, sending_results, 0);
        broadcast(world, the_PST_input, 0);
        if (the_PST_input->use_mpi) {
            broadcast(world, reactant, 0);
            broadcast(world, frag1, 0);
            broadcast(world, frag2, 0);
            broadcast(world, temp_PST_output, 0);
        }
#endif
        the_PST = new Phasespace_Core(the_PST_input, reactant, frag1, frag2, temp_PST_output);

#if (defined(USING_MPI))
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
    delete the_PST_output;
    delete temp_PST_output;
    delete the_PST_input;
    return 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
#if (defined(USING_MPI))
void Phasespace_Calculation_Slave_Process(){
    mpi::communicator world;
    Phasespace_Molecule *reactant, *frag1, *frag2;
    Phasespace_Core *the_PST;
    Phasespace_Core_Input *the_PST_input;
    Phasespace_Result *the_PST_output;
    bool sending_results;
    unsigned i;
    broadcast(world, sending_results, 0);
    if (sending_results){
        the_PST_input = new Phasespace_Core_Input();
        broadcast(world, the_PST_input, 0);
        if (the_PST_input->use_mpi) {
            broadcast(world, reactant, 0);
            broadcast(world, frag1, 0);
            broadcast(world, frag2, 0);
            broadcast(world, the_PST_output, 0);
            the_PST = new Phasespace_Core(the_PST_input,reactant,frag1,frag2,the_PST_output); //*********************** expensive to keep recreating objects inside a loop ... consider a reset function! ******************
            delete the_PST;
            reduce(world, the_PST_output, Reduce_Phasespace_Results(),0);
            std::cout<<"Process "<<world.rank()<<" : reduced"<<std::endl;
            delete reactant;
            delete frag1;
            delete frag2;
            delete the_PST_output;
        }
        delete the_PST_input;
    }
}
//***********************************************************************************************************************************************************************************
#endif


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/*bool Phasespace_Calculation::run(){
    Phasespace_Core *the_PST;
    Phasespace_Core_Input *the_PST_input;
    Phasespace_Result *the_PST_output;
    std::string bitmap_path_and_name_start;

    persistent_output << " - Phasespace_Calculation:" << std::endl;
    the_PST_input = new Phasespace_Core_Input();
    the_PST_output = new Phasespace_Result();

    bitmap_path_and_name_start=job_path+"2F_";
    //these three can probably be in a function set_standard_PST_input_values() or something
    the_PST_input->job_path = job_path;
    the_PST_input->output_to_file = the_general_input->objects_write_to_file;
    the_PST_input->output_histograms = the_general_input->histogram_write_to_file;

    if (E_available_statistical>0){
        the_PST_input->E_total = E_total; //not sure we really need this any more ... looking more like the calculation object might do all IO ... consider removing as not used in core calculations only IO
        the_PST_input->E_available = E_available_statistical;
        the_PST_input->job_name = "PST_"+boost::lexical_cast<std::string>(E_total);

        the_PST_output->initialise(limits, reactant->probe_states);

        the_PST = new Phasespace_Core(the_PST_input, reactant, frag1, frag2, the_PST_output);

        the_PST_output->write_histograms(out_file.fp);
        the_PST_output->write_bitmaps(bitmap_path_and_name_start);
        delete the_PST;
    }else{
        persistent_output<<"Error: Insufficient energy, check input file!"<<std::endl;
        print_persistent_output();
    }
    delete the_PST_input;
    delete the_PST_output;
    return 1;
}*/
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


std::string PST_OUTPUT_DIRECTORY("PST_output/");
//***********************************************************************************************************************************************************************************
class Phasespace_Core/**
This is the base class containing the core logic required to to a single Phasespace Theory (PST) state count for a two fragment system, at a single available energy etc.
The count_states() function is where the actual counting of states occurs.
Inside the many nested loops of count_states(), quantum states that can exist (because they conserve both energy and angular momentum) are processed by process_state(const Double_Quantum_State &state) and then recorded in histograms in bin_state(const Double_Quantum_State &state).
Both process_state(const Double_Quantum_State &state) and bin_state(const Double_Quantum_State &state) are likely to be overridden in derivatives of this class which will require different treatments of allowed states that are specific to certain types of PST calculations.
Additionally, the run() function of this class, which is the function called to start the state count by the code that created an object of this class, may also be overridden.
*/{
    public:
        Output_File out_file;
        const Phasespace_Core_Input *input_data; ///will hold the passed Phasespace_Core_Input ... which will be cloned
        const Phasespace_Molecule *parent, *frag1, *frag2;
        Phasespace_Result *output_data;
        double reduced_mass; ///reduced mass of the fragments
        double grand_total_L_unconstrained;

        virtual ~Phasespace_Core();
        Phasespace_Core(){}; //C++ doesn't do virtual constructors, so we can't do much with this
        Phasespace_Core(const Phasespace_Core_Input *in_input_data, const Phasespace_Molecule *in_parent, const Phasespace_Molecule *in_frag1, const Phasespace_Molecule *in_frag2, Phasespace_Result *in_output_data);
        void copy_pointers(const Phasespace_Core_Input *in_input_data, const Phasespace_Molecule *in_parent, const Phasespace_Molecule *in_frag1, const Phasespace_Molecule *in_frag2, Phasespace_Result *in_output_data);
        virtual void create_input_pointer(){};///in derived classes this function dynamically casts pointers to the derived version of input_data, allowing access to the extra functionality
        virtual void create_output_pointer(){};///in derived classes this function dynamically casts pointers to the derived version of input_data, allowing access to the extra functionality
        virtual void create_fragment_pointers(){};///in derived classes this function dynamically casts pointers to the derived version of parent, frag1 and frag2, allowing access to the extra functionality
        virtual void initialise();
        virtual void run();
        void write_energy_state(FILE *the_file);
        void write_result(FILE *the_file);
        double get_relative_velocity(const double E_translation);
        bool pass_centrifugal_barrier(const unsigned L, const double E_translation);
        virtual void count_states();
        virtual void process_state(const Double_Quantum_State &state);
        virtual void increment_histograms(Phasespace_Result_Data *data, const Double_Quantum_State &state, const Processed_State_Data &state_data);

};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Core::~Phasespace_Core() {
    out_file.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Core::Phasespace_Core(const Phasespace_Core_Input *in_input_data, const Phasespace_Molecule *in_parent, const Phasespace_Molecule *in_frag1, const Phasespace_Molecule *in_frag2, Phasespace_Result *in_output_data){
    copy_pointers(in_input_data, in_parent, in_frag1, in_frag2, in_output_data);
    initialise();
    run();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Core::copy_pointers(const Phasespace_Core_Input *in_input_data, const Phasespace_Molecule *in_parent, const Phasespace_Molecule *in_frag1, const Phasespace_Molecule *in_frag2, Phasespace_Result *in_output_data){
    input_data=in_input_data;
    parent=in_parent;
    frag1=in_frag1;
    frag2=in_frag2;
    output_data=in_output_data;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Core::initialise(){
    std::string str_output;
    grand_total_L_unconstrained=0;
    str_output=input_data->job_name;
    if (input_data->output_to_file) {
        out_file.set_output_directory(input_data->job_path+PST_OUTPUT_DIRECTORY);
        str_output.append(".out");
        out_file.open(str_output);
    }
    reduced_mass=get_reduced_mass(frag1->mass, frag2->mass);
    if (input_data->output_to_file) fprintf(out_file.fp, "\nReduced Mass Of Fragments: %lf\n\n", reduced_mass);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------persistent_output
void Phasespace_Core::run(){
    if (input_data->E_available>0){
        count_states();
        if (input_data->output_to_file) write_result(out_file.fp);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Core::write_energy_state(FILE *the_file){
    fprintf(the_file, "\nParent J: %d\n", parent->initial_J);
    fprintf(the_file, "\nTotal Energy: %lf\n", input_data->E_total);
    fprintf(the_file, "\nAvailable Energy: %lf\n", input_data->E_available);
    fprintf(the_file, "\nDissociation Energy: %lf\n", parent->E_dissociation);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Core::write_result(FILE *the_file){
    fprintf(the_file, "\ntotal: %lf\n", output_data->data->grand_total);
    fprintf(the_file, "\ntotal (L unconstrained): %lf\n", grand_total_L_unconstrained);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Phasespace_Core::get_relative_velocity(const double E_translation){
    return frag1->get_velocity(E_translation)+frag2->get_velocity(E_translation); //The weighted average is the correct way to divide the translational energy!
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Phasespace_Core::pass_centrifugal_barrier(const unsigned L, const double E_translation){
    bool result=true; /// if no b_max or C_6 value has been assigned then do not apply centrifugal barrier
    if (parent->b_max>0) {/// apply centrifugal barrier
        result= (parent->b_max*METERS_IN_ANGSTROM) >= (   (L*HBAR)/(reduced_mass*KG_IN_AMU*get_relative_velocity(E_translation))  ); // b >= L / (Mu*V)
    }else if (parent->C_6>0) {
        result= (L*(L+1)*HBARSD)  <=  (3*reduced_mass*KG_IN_AMU*pow( (2*parent->C_6*WAVENUMBER_ANGSTROM6_TO_JOULE_METER6*pow(E_translation*JOULES_IN_WAVENUMBER,2)) ,1.0/3.0)); // (L(L+1)hbar^2) <= ( 3Mu(3C_6*E_t^2)^(1/3) )
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Core::count_states(){

    unsigned thread_rank(0), num_threads(1);
    #if (defined(USING_MPI))
        mpi::communicator world;
        dMPI_Slave_Progress_Bar the_slave_progress_bar;
        dMPI_Master_Progress_Bar the_master_progress_bar;
        if (input_data->use_mpi) {
            thread_rank=world.rank();
            num_threads=world.size();
        }
    #endif


    Double_Quantum_State state;
    unsigned to_k_f1, to_k_f2, v_max_f1, v_max_f2, n_max_f1, n_max_f2, so_max_f1, so_max_f2;
    unsigned Ja, JJ_min, JJ_max, JJ, Jb, L_max, L_min;
    double E_remaining, E_remaining_ex_vib_f2, E_remaining_ex_vib_f1, E_remaining_ex_rovib_f2, E_remaining_ex_rovib_f1;
    double num_total, degeneracy_total_L_unconstrained, num_total_L_unconstrained, num_f1, num_f2;
    if (input_data->output_to_file) fprintf(out_file.fp, "\nCounting States\n");
    v_max_f1=frag1->get_v_max(input_data->E_available); //set the max vib level that fits inside the available energy
#if (defined(USING_MPI))
    if (input_data->use_mpi) {
        the_progress_bar.block(); // block the global progress bar so that long state counts can't start it (slave processes cannot print to screen ... may need to have a silent variable or something that kills all output)
        if (world.rank()) the_slave_progress_bar.initialise(100,1,v_max_f1+1, CORE_PROGRESS_BAR, MPIPB_CORE_BAR_NUMBER);
        else the_master_progress_bar.initialise(100,1,v_max_f1+1, CORE_PROGRESS_BAR, MPIPB_CORE_BAR_NUMBER); // initialise the progress bar (this has no effect on the calculation)
    } else if (!world.rank()) the_progress_bar.initialise(100,1,v_max_f1+1, CORE_PROGRESS_BAR); // initialise the progress bar (this has no effect on the calculation)
#else
    the_progress_bar.initialise(100,1,v_max_f1+1, CORE_PROGRESS_BAR); // initialise the progress bar (this has no effect on the calculation)
#endif
    for (state.f1.v=thread_rank; state.f1.v<=v_max_f1; state.f1.v+=num_threads){ //Loop over v_f1 (the vibrational levels for fragment 1)
        //std::cout<<"v_max_f1="<<v_max_f1<<"state.f1.v="<<state.f1.v<<endl;
#if (defined(USING_MPI))
        if (input_data->use_mpi) {
            if (world.rank()) the_slave_progress_bar.go(state.f1.v, CORE_PROGRESS_BAR);
            else the_master_progress_bar.go(state.f1.v, CORE_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
        } else if (!world.rank()) the_progress_bar.go(state.f1.v, CORE_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
#else
        the_progress_bar.go(state.f1.v, CORE_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
#endif
        E_remaining_ex_vib_f1=input_data->E_available-frag1->vib_levels[state.f1.v].energy; //calculate energy remaining after vibrational energy is assigned to fragment 1
        n_max_f1=frag1->get_n_max(E_remaining_ex_vib_f1); // Calculate max N, occurs when: K=0 if prolate (like HCO), K=N if oblate (like CH3)
        for (state.f1.n=0;state.f1.n<=n_max_f1;state.f1.n++){
            //cout<<"n_max_f1="<<n_max_f1<<"state.f1.n="<<state.f1.n<<endl;
            Ja=std::min(parent->initial_J,state.f1.n); //  Establish the total degeneracy (with J conserved), from the relation: N(PST) proportional to. sum over JJ of (2*Ja+1)*(2*Jb+1)
            JJ_min=abs((int)parent->initial_J-(int)state.f1.n);  // parent->initial_J is the Parent molecule Angular Momentum (input parameter)
            JJ_max=parent->initial_J+state.f1.n; //********** this and 2 lines above are used in the most nested loop ***********
            if (frag1->linear) to_k_f1=0; //linear fragments have only 1 rotational constant, they cannot have a non-zero k
            else to_k_f1=state.f1.n; //for non-linear fragments k can be any value between 0 and n
            for (state.f1.k=0;state.f1.k<=to_k_f1;state.f1.k++){
                state.f1.E_nk=frag1->get_E_nk(state.f1.n, state.f1.k);
                if (state.f1.E_nk<E_remaining_ex_vib_f1){ //if the energy of the (n,k) state is within the available energy, it might exist
                    E_remaining_ex_rovib_f1=E_remaining_ex_vib_f1-state.f1.E_nk;
                    so_max_f1=frag1->get_so_max(E_remaining_ex_rovib_f1); //set the max spin orbit state energy level that fits inside the available energy
                    for (state.f1.so=0; state.f1.so<=so_max_f1; state.f1.so++){ //Loop over state.f1.so (the spin orbital states for fragment 1)
                        num_total_L_unconstrained=num_total=0;; // reset  total number of fragment 1 + fragment 2 states  counter  to zero
                        num_f1=frag1->get_degeneracy(state.f1.k, state.f1.v, state.f1.so); //Increment the number of fragment 1 states accessible at Eavail
                        E_remaining=E_remaining_ex_rovib_f1-frag1->spin_orbit_states[state.f1.so].energy;//Now calculate the energy available to fragment 2 and translation (the remaining energy that can be distributed to fragment 2)
                        v_max_f2=frag2->get_v_max(E_remaining); //set the max vib level that fits inside the remaining energy
                        for (state.f2.v=0; state.f2.v<=v_max_f2; state.f2.v++){ //Loop over the vibrational levels for fragment 2
                            E_remaining_ex_vib_f2=E_remaining-frag2->vib_levels[state.f2.v].energy; //calculate energy remaining after vibrational energy is assigned
                            n_max_f2=frag2->get_n_max(E_remaining_ex_vib_f2); // Calculate max N, occurs when: K=0 if prolate (like HCO), K=N if oblate (like CH3)
                            for (state.f2.n=0; state.f2.n<=n_max_f2; state.f2.n++){
                                if (frag2->linear) to_k_f2=0; //linear fragments have only 1 rotational constant, they cannot have a non-zero k
                                else to_k_f2=state.f2.n; //for non-linear fragments k can be any value between 0 and n
                                for (state.f2.k=0; state.f2.k<=to_k_f2; state.f2.k++){
                                    state.f2.E_nk=frag2->get_E_nk(state.f2.n, state.f2.k);
                                    if (state.f2.E_nk<=E_remaining_ex_vib_f2){ //if the energy of this (N,K) state fits inside the remaining energy
                                        E_remaining_ex_rovib_f2=E_remaining_ex_vib_f2-state.f2.E_nk;
                                        so_max_f2=frag2->get_so_max(E_remaining_ex_rovib_f2); //set the max spin orbit state energy level that fits inside the available energy
                                        for (state.f2.so=0; state.f2.so<=so_max_f2; state.f2.so++){ //Loop over state.f1.so (the spin orbital states for fragment 1)
                                            num_f2=frag2->get_degeneracy(state.f2.k, state.f2.v, state.f2.so); //  Count the number of F2 states (for each F1 state)
                                            state.E_translation=E_remaining_ex_rovib_f2-frag2->spin_orbit_states[state.f2.so].energy;	//Calculate the energy distributed into translation
                                            for (JJ=JJ_min; JJ<=JJ_max; JJ++){  //so this is for jj = abs(parent->initial_J - N) to (parent->initial_J + N)  This is how parent->initial_J is conserved (triangle inequality) ... ALSO IMPLEMENTS THE 2J+1 DEGENERACY OF ALL J STATES
                                                Jb=std::min(JJ,state.f2.n);
                                                degeneracy_total_L_unconstrained=(2*Ja+1)*(2*Jb+1)*num_f2*num_f1; // ***** THIS IS NOT CURRENTLY WORKING PROPERLY! *****    //Incorporate the f1 and f2 degeneracy into the total deg count:
                                                num_total_L_unconstrained+=degeneracy_total_L_unconstrained; 	//And thus the total number of f1 and f2 states without any L constraints
                                                L_max=JJ+state.f2.n; //This is the triangle inequality to conserve total JJ when getting all possible L
                                                L_min=abs((int)JJ-(int)state.f2.n);
                                                for (state.L=L_min; state.L<=L_max; state.L++){ //we have to do it this way for centrifugal barriers etc ... 2J+1 degeneracy accounted for in these loops
                                                    if (pass_centrifugal_barrier(state.L, state.E_translation)) {
                                                        state.degeneracy=num_f2*num_f1;
                                                        num_total+=state.degeneracy; //count num that can form radicals and the num that cant
                                                        process_state(state);
                                                    }
                                                }
                                            } //end for (JJ=JJ_min; JJ<=JJ_maxl JJ++)
                                        } //end for (state.f2.so=0; state.f2.so<=so_max_f2; state.f2.so++)
                                    } //end if (state.f2.E_nk<=E_remaining_ex_vib_f2)
                                } //end for (state.f2.k=0; state.f2.k<=state.f2.n; state.f2.k++)
                            } //end for (state.f2.n=0; state.f2.n<=n_max_f2; state.f2.n++)
                        } //end for (state.f2.v=0; state.f2.v<=v_max_f2; state.f2.v++)
                        grand_total_L_unconstrained+=num_total_L_unconstrained; //pretty sure we'll be getting rid of this one
                    } //end for (state.f1.so=0; state.f1.so<=so_max_f1; state.f1.so++)
                } //end if (state.f1.E_nk<input_data->E_available)
            } //end for (state.f1.k=0;state.f1.k<=state.f1.n;state.f1.k++)
        } //end for (state.f1.n=0;state.f1.n<=n_max_f1;state.f1.n++)
    }
#if (defined(USING_MPI))
    if (input_data->use_mpi) {
        the_progress_bar.un_block();//release the global progress bar too because we blocked it when we initialised this MPI progress bar.
        if (world.rank()) the_slave_progress_bar.un_initialise(CORE_PROGRESS_BAR);//release the progress bar
        else the_master_progress_bar.un_initialise(CORE_PROGRESS_BAR);//release the progress bar
    } else if (!world.rank()) the_progress_bar.un_initialise(CORE_PROGRESS_BAR);//release the progress bar
#else
    the_progress_bar.un_initialise(CORE_PROGRESS_BAR);//release the progress bar
#endif
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Core::process_state(const Double_Quantum_State &state){
    unsigned p;
    Processed_State_Data state_data;
    std::vector<unsigned> probe_ids; //store the ids of all relevant probes

    for (p=0; p<output_data->probes.size(); p++){
        if (output_data->probes[p]->state.applies(state)) {
            probe_ids.push_back(p);
            output_data->probes[p]->data->grand_total+=state.degeneracy;
        }
    }
    output_data->data->grand_total+=state.degeneracy;

    if (input_data->output_histograms){
        state_data.E_trans=state.E_translation;
        state_data.vel_f1=frag1->get_velocity(state.E_translation);
        state_data.vel_f2=frag2->get_velocity(state.E_translation);
        state_data.E_trans_f1=frag1->get_E_trans(state_data.vel_f1);
        state_data.E_trans_f2=frag2->get_E_trans(state_data.vel_f2);
        state_data.E_vib_f1=frag1->vib_levels[state.f1.v].energy;
        state_data.E_vib_f2=frag2->vib_levels[state.f2.v].energy;
        state_data.E_int_f1=state_data.E_vib_f1+state.f1.E_nk+frag1->spin_orbit_states[state.f1.so].energy; //E_total=E_vib+E_rot+E_soc
        state_data.E_int_f2=state_data.E_vib_f2+state.f2.E_nk+frag2->spin_orbit_states[state.f2.so].energy; //E_total=E_vib+E_rot+E_soc

        increment_histograms(output_data->data, state, state_data);
        for (p=0; p<probe_ids.size(); p++) increment_histograms(output_data->probes[probe_ids[p]]->data, state, state_data);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Core::increment_histograms(Phasespace_Result_Data *data, const Double_Quantum_State &state, const Processed_State_Data &state_data){
    data->E_trans.increment_value(state_data.E_trans, state.degeneracy);
    data->frag1_vel.increment_value(state_data.vel_f1, state.degeneracy);
    data->frag1_E_trans.increment_value(state_data.E_trans_f1, state.degeneracy);
    data->frag2_vel.increment_value(state_data.vel_f2, state.degeneracy);
    data->frag2_E_trans.increment_value(state_data.E_trans_f2, state.degeneracy);
    if (!frag1->single_atom){
        data->frag1_v.increment_value(state.f1.v, state.degeneracy);
        data->frag1_E_vib.increment_value(state_data.E_vib_f1, state.degeneracy);
        data->frag1_E_int.increment_value(state_data.E_int_f1, state.degeneracy);
        data->frag1_E_int_trans.increment_value(state_data.E_int_f1, state_data.E_trans_f1, state.degeneracy);
        data->frag1_nk.increment_value(state.f1.n, state.f1.k, state.degeneracy);
        data->frag1_E_rot.increment_value(state.f1.E_nk, state.degeneracy);
    }
    if (!frag2->single_atom){
        data->frag2_v.increment_value(state.f2.v, state.degeneracy);
        data->frag2_E_vib.increment_value(state_data.E_vib_f2, state.degeneracy);
        data->frag2_E_int.increment_value(state_data.E_int_f2, state.degeneracy);
        data->frag2_E_int_trans.increment_value(state_data.E_int_f2, state_data.E_trans_f2, state.degeneracy);
        data->frag2_nk.increment_value(state.f2.n, state.f2.k, state.degeneracy);
        data->frag2_E_rot.increment_value(state.f2.E_nk, state.degeneracy);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//***********************************************************************************************************************************************************************************

#if (defined(USING_MPI))
//this is the binary function object used for the MPI reduce routine at the end of a 3F calculation
struct Reduce_Phasespace_Results : public std::binary_function<Phasespace_Result*, Phasespace_Result*, Phasespace_Result*>
{
    Phasespace_Result* operator() (Phasespace_Result *a, Phasespace_Result *b){
        (*b)+=(*a);
        return b;
    }
};
namespace boost { namespace mpi {//this allows boost to use a more efficient gather algorithm, by telling boost that Reduce_Phasespace_Results commutes
    template<>
    struct is_commutative<Reduce_Phasespace_Results, Phasespace_Result*> : mpl::true_{};
}}
#endif
//***********************************************************************************************************************************************************************************


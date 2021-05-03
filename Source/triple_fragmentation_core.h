std::string SECONDARY_STATES_OUTPUT_DIRECTORY("secondary_states/");
std::string SECONDARY_STATES_FILE_NAME("frag_");
//***********************************************************************************************************************************************************************************
class Primary_Dissociation_Result : public Phasespace_Result
{
    public:
        double f1_vel_max, f2_vel_max; //the maximum velocity enumerated for states that will undergo secondary dissociation, for each fragment
        unsigned f1_num_secondary, f2_num_secondary, f1_num_secondary_master_thread, f2_num_secondary_master_thread, f1_max_bend_n, f2_max_bend_n;

        virtual ~Primary_Dissociation_Result(){};
        Primary_Dissociation_Result(){};
        Primary_Dissociation_Result(const Primary_Dissociation_Result &rhs);
        Primary_Dissociation_Result& operator=(const Primary_Dissociation_Result &rhs);
        Primary_Dissociation_Result& operator+=(const Primary_Dissociation_Result &rhs);
        virtual void initialise(const Histogram_Limits &limits, const Phasespace_Molecule *reactant);
        void store_master_thread_counts();
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Phasespace_Result>(*this);
            ar & f1_vel_max & f2_vel_max & f1_num_secondary & f2_num_secondary & f1_num_secondary_master_thread & f2_num_secondary_master_thread & f1_max_bend_n & f2_max_bend_n;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Primary_Dissociation_Result::Primary_Dissociation_Result(const Primary_Dissociation_Result &rhs):Phasespace_Result(rhs){
    f1_num_secondary=rhs.f1_num_secondary;
    f2_num_secondary=rhs.f2_num_secondary;
    f1_num_secondary_master_thread=rhs.f1_num_secondary_master_thread;
    f2_num_secondary_master_thread=rhs.f2_num_secondary_master_thread;
    f1_max_bend_n=rhs.f1_max_bend_n;
    f2_max_bend_n=rhs.f2_max_bend_n;
    f1_vel_max=rhs.f1_vel_max;
    f2_vel_max=rhs.f2_vel_max;
#if (defined(USING_MPI))
    mpi::communicator world;
//    std::cout<<"| Process "<<world.rank()<<" : IN THE Primary_Dissociation_Result COPY CONSTRUCTOR!"<<std::endl;
#endif
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Primary_Dissociation_Result& Primary_Dissociation_Result::operator=(const Primary_Dissociation_Result &rhs){
    Phasespace_Result::operator=(rhs);
    f1_num_secondary=rhs.f1_num_secondary;
    f2_num_secondary=rhs.f2_num_secondary;
    f1_num_secondary_master_thread=rhs.f1_num_secondary_master_thread;
    f2_num_secondary_master_thread=rhs.f2_num_secondary_master_thread;
    f1_max_bend_n=rhs.f1_max_bend_n;
    f2_max_bend_n=rhs.f2_max_bend_n;
    f1_vel_max=rhs.f1_vel_max;
    f2_vel_max=rhs.f2_vel_max;
#if (defined(USING_MPI))
    mpi::communicator world;
 //   std::cout<<"| Process "<<world.rank()<<" : IN THE Primary_Dissociation_Result ASSIGNMENT OPERATOR!"<<std::endl;
#endif
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Primary_Dissociation_Result& Primary_Dissociation_Result::operator+=(const Primary_Dissociation_Result &rhs){
    Phasespace_Result::operator+=(rhs);
    f1_num_secondary+=rhs.f1_num_secondary;
    f2_num_secondary+=rhs.f2_num_secondary;
    f1_num_secondary_master_thread+=rhs.f1_num_secondary_master_thread; //these 2 are only non zero in the master thread, so they will still hold the correct value after reduction
    f2_num_secondary_master_thread+=rhs.f2_num_secondary_master_thread;
    //the vel_max and max_bend_n store the highest value encountered for all states, so compare the two value and keep the highest
    if (f1_max_bend_n<rhs.f1_max_bend_n) f1_max_bend_n=rhs.f1_max_bend_n;
    if (f2_max_bend_n<rhs.f2_max_bend_n) f2_max_bend_n=rhs.f2_max_bend_n;
    if (f1_vel_max<rhs.f1_vel_max) f1_vel_max=rhs.f1_vel_max;
    if (f2_vel_max<rhs.f2_vel_max) f2_vel_max=rhs.f2_vel_max;
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Primary_Dissociation_Result::initialise(const Histogram_Limits &limits, const Phasespace_Molecule *reactant) {
    Phasespace_Result::initialise(limits, reactant);
    f1_vel_max=f2_vel_max=0;
    f1_num_secondary=f2_num_secondary=f1_num_secondary_master_thread=f2_num_secondary_master_thread=f1_max_bend_n=f2_max_bend_n=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Primary_Dissociation_Result::store_master_thread_counts() { //we store the master thread counts before we reduce so that the compress sec states file progress bar can be correct for the master thread
    f1_num_secondary_master_thread=f1_num_secondary; //if we are running in a single thread, then these will be the same anyway (no reduction will occur)
    f2_num_secondary_master_thread=f2_num_secondary;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
class Triple_Frag_Primary : public Impulsive_Phasespace_Core
{
    public:
        Primary_Dissociation_Result *t_output_data;
        Impulsive_Core_Input *t_input_data;
        Histogram_3d *sec_states1, *sec_states2;
        Output_File secondary_states_frag1, secondary_states_frag2;
        unsigned use_secondary_states_histograms;

        virtual ~Triple_Frag_Primary();
        Triple_Frag_Primary(){};
        Triple_Frag_Primary(const Impulsive_Core_Input *in_input_data, const Impulsive_Molecule *in_parent, const Impulsive_Molecule *in_frag1, const Impulsive_Molecule *in_frag2, Phasespace_Result *in_output_data);
        Triple_Frag_Primary(const Impulsive_Core_Input *in_input_data, const Impulsive_Molecule *in_parent, const Impulsive_Molecule *in_frag1, const Impulsive_Molecule *in_frag2, Phasespace_Result *in_output_data,
                            Histogram_3d *in_sec_states1, Histogram_3d *in_sec_states2);
        virtual void create_output_pointer();
        virtual void initialise();
        virtual void process_state(const Double_Quantum_State &state);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Triple_Frag_Primary::~Triple_Frag_Primary() {
	if (frag1->E_dissociation) secondary_states_frag1.close();
	if (frag2->E_dissociation) secondary_states_frag2.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Triple_Frag_Primary::Triple_Frag_Primary(const Impulsive_Core_Input *in_input_data, const Impulsive_Molecule *in_parent, const Impulsive_Molecule *in_frag1, const Impulsive_Molecule *in_frag2, Phasespace_Result *in_output_data):
    Impulsive_Phasespace_Core(){
    use_secondary_states_histograms=0;
    copy_pointers(in_input_data, in_parent, in_frag1, in_frag2, in_output_data);
    create_input_pointer();
    create_output_pointer();
    create_fragment_pointers();
    initialise();
    run();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Triple_Frag_Primary::Triple_Frag_Primary(const Impulsive_Core_Input *in_input_data, const Impulsive_Molecule *in_parent, const Impulsive_Molecule *in_frag1, const Impulsive_Molecule *in_frag2, Phasespace_Result *in_output_data, Histogram_3d *in_sec_states1, Histogram_3d *in_sec_states2):
    Impulsive_Phasespace_Core(){
    use_secondary_states_histograms=1;
    sec_states1=in_sec_states1;
    sec_states2=in_sec_states2;
    copy_pointers(in_input_data, in_parent, in_frag1, in_frag2, in_output_data);
    create_input_pointer();
    create_output_pointer();
    create_fragment_pointers();
    initialise();
    run();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Frag_Primary::create_output_pointer(){
    t_output_data = dynamic_cast<Primary_Dissociation_Result*>(output_data); ///create a pointer to the output_data object that allows access to the derived fuctionality
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Frag_Primary::initialise(){
    Impulsive_Phasespace_Core::initialise();
#if (defined(USING_MPI))
    mpi::communicator world;
#endif
    if (frag1->E_dissociation) { //if a dissociation energy has been defined, putput states that can undergo secondary dissociation for this fragment
#if (defined(USING_MPI)) ///consider using binary files for this too - faster for large files etc
        secondary_states_frag1.set_output_directory(scratch_directory);
        secondary_states_frag1.open(SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(frag1->number)+"_t"+boost::lexical_cast<std::string>(world.rank())+".csv");
#else
        secondary_states_frag1.set_output_directory(input_data->job_path+SECONDARY_STATES_OUTPUT_DIRECTORY);
        secondary_states_frag1.open(SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(frag1->number)+".csv");
#endif
        fprintf(secondary_states_frag1.fp, "E_internal,degen,J,E_trans,v_f1,bend_n\n");
    }///I SMELL A FUNCTION!!
    if (frag2->E_dissociation) { //if a dissociation energy has been defined, putput states that can undergo secondary dissociation for this fragment
#if (defined(USING_MPI)) ///consider using binary files for this too - faster for large files etc
        secondary_states_frag2.set_output_directory(scratch_directory);
        secondary_states_frag2.open(SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(frag2->number)+"_t"+boost::lexical_cast<std::string>(world.rank())+".csv");
#else
        secondary_states_frag2.set_output_directory(input_data->job_path+SECONDARY_STATES_OUTPUT_DIRECTORY);
        secondary_states_frag2.open(SECONDARY_STATES_FILE_NAME+boost::lexical_cast<std::string>(frag2->number)+".csv");
#endif
        fprintf(secondary_states_frag2.fp, "E_internal,degen,J,E_trans,v_f2,bend_n\n");
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Frag_Primary::process_state(const Double_Quantum_State &state){ //if there is enough internal energy for second dissociation write state to file
    Double_Quantum_State state_with_impulse;
    Processed_State_Data state_data;
    double E_wrongBy;
    unsigned i, p;
    std::vector<unsigned> probe_ids; //store the ids of all relevant probes
///the different bit (for writing the file) - start
    double E_excess; ///check whether we agree that spin orbit coupling state energy is ergodic (ie, is it really a "degree of freedom")?
    unsigned bend_n_f1, bend_n_f2;
///the different bit (for writing the file) - end
///SAME IN ALL 3 (impulsive_core, triple_primary, triple_secondary) - start

    for (i=0; i<bond_angle_distribution.size(); i++){ //if not using angle distribution, this will contain a single element with eq_angle and prob=1

        apply_impulse(state, state_with_impulse, E_wrongBy, i);

        if (E_wrongBy <= state.E_translation || state.tunneled){ //Energy can be conserved, this state can exist (energy correction must come from the statistical reservoir)

            state_with_impulse.degeneracy=state.degeneracy*bond_angle_distribution[i].probability*tunneling_probability;

            output_data->data->grand_total+=state_with_impulse.degeneracy;
            probe_ids.clear();
            for (p=0; p<output_data->probes.size(); p++){
                if (output_data->probes[p]->state.applies(state_with_impulse)){
                    probe_ids.push_back(p);
                    output_data->probes[p]->data->grand_total+=state_with_impulse.degeneracy;
                }
            }
            state_data.E_trans=state_with_impulse.E_translation;
            state_data.vel_f1=i_frag1->get_velocity(state_with_impulse.E_translation);
            state_data.vel_f2=i_frag2->get_velocity(state_with_impulse.E_translation);
            state_data.E_trans_f1=i_frag1->get_E_trans(state_data.vel_f1);
            state_data.E_trans_f2=i_frag2->get_E_trans(state_data.vel_f2);
            state_data.E_vib_f1=i_frag1->vib_levels[state_with_impulse.f1.v].energy;
            state_data.E_vib_f2=i_frag2->vib_levels[state_with_impulse.f2.v].energy;
            state_data.E_int_f1=state_data.E_vib_f1+state_with_impulse.f1.E_nk+i_frag1->spin_orbit_states[state_with_impulse.f1.so].energy; //E_total=E_vib+E_rot+E_soc
            state_data.E_int_f2=state_data.E_vib_f2+state_with_impulse.f2.E_nk+i_frag2->spin_orbit_states[state_with_impulse.f2.so].energy; //E_total=E_vib+E_rot+E_soc
            increment_histograms(output_data->data, state_with_impulse, state_data);
            for (p=0; p<probe_ids.size(); p++) increment_histograms(output_data->probes[probe_ids[p]]->data, state_with_impulse, state_data);

            ///the different bit (for writing the file) - start
            if (i_frag1->E_dissociation) {
                E_excess=state_data.E_int_f1-i_frag1->E_dissociation;
                if (E_excess>0) {
                    bend_n_f1=i_frag1->vib_levels[state_with_impulse.f1.v].mode_tally[i_frag1->bend_mode];
                    if(state_data.vel_f1 > t_output_data->f1_vel_max) t_output_data->f1_vel_max=state_data.vel_f1;
                    if(bend_n_f1>t_output_data->f1_max_bend_n) t_output_data->f1_max_bend_n=bend_n_f1;
                    t_output_data->f1_num_secondary++;
                    if (use_secondary_states_histograms) sec_states1->increment_value(E_excess, state_data.vel_f1, state_with_impulse.f1.n, state_with_impulse.degeneracy);
                    else fprintf(secondary_states_frag1.fp, "%lf,%lf,%d,%lf\n", state_data.E_int_f1, state_with_impulse.degeneracy, state_with_impulse.f1.n, state_data.vel_f1);
                }
            }///I SMELL A FUNCTION!!
            if (i_frag2->E_dissociation) {
                E_excess=state_data.E_int_f2-i_frag2->E_dissociation;
                if (E_excess>0) {
                    bend_n_f2=i_frag2->vib_levels[state_with_impulse.f2.v].mode_tally[i_frag2->bend_mode];
                    if(state_data.vel_f2>t_output_data->f2_vel_max) t_output_data->f2_vel_max=state_data.vel_f2;
                    if(bend_n_f2>t_output_data->f2_max_bend_n) t_output_data->f2_max_bend_n=bend_n_f2;
                    t_output_data->f2_num_secondary++;
                    if (use_secondary_states_histograms) sec_states2->increment_value(E_excess, state_data.vel_f2, state_with_impulse.f2.n, state_with_impulse.degeneracy);
                    else fprintf(secondary_states_frag2.fp, "%lf,%lf,%d,%lf\n", state_data.E_int_f2, state_with_impulse.degeneracy, state_with_impulse.f2.n, state_data.vel_f2);
                }
            }
            ///the different bit (for writing the file) - end
        }
    }
    ///below here use state_with_impulse ... also ... don't write states that don't exist (E_wrongBy <= state.E_translation)??

}

//***********************************************************************************************************************************************************************************
#if (defined(USING_MPI))
//this is the binary function object used for the MPI reduce routine at the end of a 3F calculation
struct Reduce_Primary_Results : public std::binary_function<Primary_Dissociation_Result*, Primary_Dissociation_Result*, Primary_Dissociation_Result*>
{
    Primary_Dissociation_Result* operator() (Primary_Dissociation_Result *a, Primary_Dissociation_Result *b){
        (*b)+=(*a);
        return b;
    }
};
namespace boost { namespace mpi {//this allows boost to use a more efficient gather algorithm, by telling boost that Reduce_Phasespace_Results commutes
    template<>
    struct is_commutative<Reduce_Primary_Results, Primary_Dissociation_Result*> : mpl::true_{};
}}
#endif
//***********************************************************************************************************************************************************************************
//***********************************************************************************************************************************************************************************
class Secondary_Dissociation_Result : public Phasespace_Result
{
    public:
        double grand_total_no_tunneling_probability, total_parent_states_considered; //total_parent_states_considered will be different to f1_num_secondary or f2_num_secondary because this will include the degeneracy of the states

        virtual ~Secondary_Dissociation_Result(){}
        Secondary_Dissociation_Result(){}
        Secondary_Dissociation_Result(const Secondary_Dissociation_Result &rhs);
        Secondary_Dissociation_Result& operator=(const Secondary_Dissociation_Result &rhs);
        Secondary_Dissociation_Result& operator+=(const Secondary_Dissociation_Result &rhs);
        virtual void initialise(const Histogram_Limits &limits, const Phasespace_Molecule *reactant);
        virtual void reset();
        void increment(const double primary_degeneracy, Secondary_Dissociation_Result* secondary_PST_result);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Phasespace_Result>(*this);
            ar & grand_total_no_tunneling_probability & total_parent_states_considered;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Secondary_Dissociation_Result::Secondary_Dissociation_Result(const Secondary_Dissociation_Result &rhs):Phasespace_Result(rhs){
    grand_total_no_tunneling_probability=rhs.grand_total_no_tunneling_probability;
    total_parent_states_considered=rhs.total_parent_states_considered;
#if (defined(USING_MPI))
    mpi::communicator world;
//    std::cout<<"| Process "<<world.rank()<<" : IN THE Secondary_Dissociation_Result COPY CONSTRUCTOR!"<<std::endl;
#endif
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Secondary_Dissociation_Result& Secondary_Dissociation_Result::operator=(const Secondary_Dissociation_Result &rhs){
    Phasespace_Result::operator=(rhs);
    grand_total_no_tunneling_probability=rhs.grand_total_no_tunneling_probability;
    total_parent_states_considered=rhs.total_parent_states_considered;
#if (defined(USING_MPI))
    mpi::communicator world;
//    std::cout<<"| Process "<<world.rank()<<" : IN THE Secondary_Dissociation_Result ASSIGNMENT OPERATOR!"<<std::endl;
#endif
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Secondary_Dissociation_Result& Secondary_Dissociation_Result::operator+=(const Secondary_Dissociation_Result &rhs){
    Phasespace_Result::operator+=(rhs);
    grand_total_no_tunneling_probability+=rhs.grand_total_no_tunneling_probability;
    total_parent_states_considered+=rhs.total_parent_states_considered;
#if (defined(USING_MPI))
    mpi::communicator world;
//    std::cout<<"| Process "<<world.rank()<<" : IN THE Secondary_Dissociation_Result INCREMENT OPERATOR!"<<std::endl;
#endif
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Secondary_Dissociation_Result::initialise(const Histogram_Limits &limits, const Phasespace_Molecule *reactant) {
    Phasespace_Result::initialise(limits, reactant);
    grand_total_no_tunneling_probability=0;
    total_parent_states_considered=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Secondary_Dissociation_Result::reset(){
    Phasespace_Result::reset();
    grand_total_no_tunneling_probability=0;
    total_parent_states_considered=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Secondary_Dissociation_Result::increment(const double primary_degeneracy, Secondary_Dissociation_Result* secondary_PST_result) {
    unsigned p; //The reason we divide here is that this way scaling factor will be very large and we get a result. Doing it the other way makes it VERY small, so it is morelikely to be too small to be distinguishable from 0 at double precision
    double scaling_factor=secondary_PST_result->grand_total_no_tunneling_probability/primary_degeneracy; //we divide by this factor, so it's like multiplying by primary degen and dividing by grand total, that makes sense!
    if (scaling_factor){ //if grand_total=0 then no need to update histograms
        data->grand_total+=secondary_PST_result->data->grand_total/scaling_factor; //grand_total_no_tunneling_probability is different from grand_total only when there is tunneling involved ...  grand_total is the total with tunneling probability applied, grand_total_no_tunneling_probability is without tunneling prob applied ... thus with no tunneling, this gives back primary_degeneracy (ie tunneling alone gives primary_degeneracy*tunneling_degeneracy)
        data->frag1_v.inner_range_increment(secondary_PST_result->data->frag1_v, scaling_factor);
        data->frag1_nk.inner_range_increment(secondary_PST_result->data->frag1_nk, scaling_factor);
        data->frag2_v.inner_range_increment(secondary_PST_result->data->frag2_v, scaling_factor);
        data->frag2_nk.inner_range_increment(secondary_PST_result->data->frag2_nk, scaling_factor);
        data->E_trans.range_increment(secondary_PST_result->data->E_trans, scaling_factor); //this is function I  have added, BUT
        data->frag1_vel.range_increment(secondary_PST_result->data->frag1_vel, scaling_factor);
        data->frag1_E_trans.range_increment(secondary_PST_result->data->frag1_E_trans, scaling_factor);
        data->frag1_E_int.range_increment(secondary_PST_result->data->frag1_E_int, scaling_factor);
        data->frag1_E_vib.range_increment(secondary_PST_result->data->frag1_E_vib, scaling_factor);
        data->frag1_E_rot.range_increment(secondary_PST_result->data->frag1_E_rot, scaling_factor);
        data->frag1_E_int_trans.inner_range_increment(secondary_PST_result->data->frag1_E_int_trans, scaling_factor);
        data->frag2_vel.range_increment(secondary_PST_result->data->frag2_vel, scaling_factor);
        data->frag2_E_trans.range_increment(secondary_PST_result->data->frag2_E_trans, scaling_factor);
        data->frag2_E_int.range_increment(secondary_PST_result->data->frag2_E_int, scaling_factor);
        data->frag2_E_vib.range_increment(secondary_PST_result->data->frag2_E_vib, scaling_factor);
        data->frag2_E_rot.range_increment(secondary_PST_result->data->frag2_E_rot, scaling_factor);
        data->frag2_E_int_trans.inner_range_increment(secondary_PST_result->data->frag2_E_int_trans, scaling_factor);
        for (p=0; p<probes.size(); p++) {
            probes[p]->data->grand_total+=secondary_PST_result->probes[p]->data->grand_total/scaling_factor;
            probes[p]->data->frag1_v.inner_range_increment(secondary_PST_result->probes[p]->data->frag1_v, scaling_factor);
            probes[p]->data->frag1_nk.inner_range_increment(secondary_PST_result->probes[p]->data->frag1_nk, scaling_factor);
            probes[p]->data->frag2_v.inner_range_increment(secondary_PST_result->probes[p]->data->frag2_v, scaling_factor);
            probes[p]->data->frag2_nk.inner_range_increment(secondary_PST_result->probes[p]->data->frag2_nk, scaling_factor);
            probes[p]->data->E_trans.range_increment(secondary_PST_result->probes[p]->data->E_trans, scaling_factor);
            probes[p]->data->frag1_vel.range_increment(secondary_PST_result->probes[p]->data->frag1_vel, scaling_factor);
            probes[p]->data->frag1_E_trans.range_increment(secondary_PST_result->probes[p]->data->frag1_E_trans, scaling_factor);
            probes[p]->data->frag1_E_int.range_increment(secondary_PST_result->probes[p]->data->frag1_E_int, scaling_factor);
            probes[p]->data->frag1_E_vib.range_increment(secondary_PST_result->probes[p]->data->frag1_E_vib, scaling_factor);
            probes[p]->data->frag1_E_rot.range_increment(secondary_PST_result->probes[p]->data->frag1_E_rot, scaling_factor);
            probes[p]->data->frag1_E_int_trans.inner_range_increment(secondary_PST_result->probes[p]->data->frag1_E_int_trans, scaling_factor);
            probes[p]->data->frag2_vel.range_increment(secondary_PST_result->probes[p]->data->frag2_vel, scaling_factor);
            probes[p]->data->frag2_E_trans.range_increment(secondary_PST_result->probes[p]->data->frag2_E_trans, scaling_factor);
            probes[p]->data->frag2_E_int.range_increment(secondary_PST_result->probes[p]->data->frag2_E_int, scaling_factor);
            probes[p]->data->frag2_E_vib.range_increment(secondary_PST_result->probes[p]->data->frag2_E_vib, scaling_factor);
            probes[p]->data->frag2_E_rot.range_increment(secondary_PST_result->probes[p]->data->frag2_E_rot, scaling_factor);
            probes[p]->data->frag2_E_int_trans.inner_range_increment(secondary_PST_result->probes[p]->data->frag2_E_int_trans, scaling_factor);
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
class Triple_Frag_Secondary : public Impulsive_Phasespace_Core
{
    public:
        Secondary_Dissociation_Result *t_output_data;
        const Triple_Fragmentation_Core_Input *t_input_data;

        virtual ~Triple_Frag_Secondary();
        Triple_Frag_Secondary(){};
        Triple_Frag_Secondary(const Triple_Fragmentation_Core_Input *in_input_data, const Impulsive_Molecule *in_parent,
                                             const Impulsive_Molecule *in_frag1, const Impulsive_Molecule *in_frag2, Secondary_Dissociation_Result *in_output_data);
        virtual void create_input_pointer();
        virtual void create_output_pointer();
        virtual void initialise();
        virtual void process_state(const Double_Quantum_State &state);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Triple_Frag_Secondary::~Triple_Frag_Secondary() {
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Triple_Frag_Secondary::Triple_Frag_Secondary(const Triple_Fragmentation_Core_Input *in_input_data, const Impulsive_Molecule *in_parent,
                                             const Impulsive_Molecule *in_frag1, const Impulsive_Molecule *in_frag2, Secondary_Dissociation_Result *in_output_data){
    copy_pointers(in_input_data, in_parent, in_frag1, in_frag2, in_output_data);
    create_input_pointer();
    create_output_pointer();
    create_fragment_pointers();
    initialise();
    run();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Frag_Secondary::create_input_pointer(){
    Impulsive_Phasespace_Core::create_input_pointer();
    t_input_data = dynamic_cast<const Triple_Fragmentation_Core_Input*>(input_data);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Frag_Secondary::create_output_pointer(){
    t_output_data = dynamic_cast<Secondary_Dissociation_Result*>(output_data); ///create a pointer to the output_data object that allows access to the derived fuctionality
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Frag_Secondary::initialise(){
    Impulsive_Phasespace_Core::initialise();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Triple_Frag_Secondary::process_state(const Double_Quantum_State &state){
    double E_wrongBy;
    double intermediate_frame_vel_f1, intermediate_frame_vel_f2, intermediate_frame_E_trans_f1, intermediate_frame_E_trans_f2;
    unsigned i, j, p;
    std::vector<unsigned> probe_ids; //store the ids of all relevant probes

    Double_Quantum_State state_with_impulse;
    Processed_State_Data state_data;
    for (i=0; i<bond_angle_distribution.size(); i++){ //if not using BOND angle distribution, this will contain a single element with eq_angle and prob=1

        apply_impulse(state, state_with_impulse, E_wrongBy, i);

        if (E_wrongBy <= state.E_translation || state.tunneled){ //Energy can be conserved, this state can exist (energy correction must come from the statistical reservoir)

            intermediate_frame_vel_f1=i_frag1->get_velocity(state_with_impulse.E_translation);
            intermediate_frame_vel_f2=i_frag2->get_velocity(state_with_impulse.E_translation);

            probe_ids.clear();
            for (p=0; p<output_data->probes.size(); p++) if (output_data->probes[p]->state.applies(state_with_impulse)) probe_ids.push_back(p);

            for (j=0; j<t_input_data->velocity_angle_distribution.size(); j++){

                state_with_impulse.degeneracy=state.degeneracy*bond_angle_distribution[i].probability*t_input_data->velocity_angle_distribution[j].probability*tunneling_probability;

                state_data.E_trans=state_with_impulse.E_translation;
                state_data.vel_f1=sqrt( pow(i_parent->initial_velocity,2)+pow(intermediate_frame_vel_f1,2)-2*i_parent->initial_velocity*intermediate_frame_vel_f1*cos(t_input_data->velocity_angle_distribution[j].radians) ); //magnitude of the vector sum of velocity primary and secondary, at angle j in velocity_angle_distribution (using cosine rule)
                state_data.vel_f2=sqrt( pow(i_parent->initial_velocity,2)+pow(intermediate_frame_vel_f2,2)-2*i_parent->initial_velocity*intermediate_frame_vel_f2*cos(t_input_data->velocity_angle_distribution[j].radians) ); //velocity_angle_distribution ( sin(theta) ) has a period of 180 degrees so this is the same as including the 180 degree difference in direction between frag 1 and 2
                state_data.E_trans_f1=i_frag1->get_E_trans(state_data.vel_f1);
                state_data.E_trans_f2=i_frag2->get_E_trans(state_data.vel_f2);
                state_data.E_vib_f1=i_frag1->vib_levels[state_with_impulse.f1.v].energy;
                state_data.E_vib_f2=i_frag2->vib_levels[state_with_impulse.f2.v].energy;
                state_data.E_int_f1=state_data.E_vib_f1+state_with_impulse.f1.E_nk+i_frag1->spin_orbit_states[state_with_impulse.f1.so].energy; //E_total=E_vib+E_rot+E_soc
                state_data.E_int_f2=state_data.E_vib_f2+state_with_impulse.f2.E_nk+i_frag2->spin_orbit_states[state_with_impulse.f2.so].energy; //E_total=E_vib+E_rot+E_soc

                t_output_data->grand_total_no_tunneling_probability+=(state.degeneracy*bond_angle_distribution[i].probability*t_input_data->velocity_angle_distribution[j].probability);

                output_data->data->grand_total+=state_with_impulse.degeneracy;
                for (p=0; p<probe_ids.size(); p++) output_data->probes[p]->data->grand_total+=state_with_impulse.degeneracy;


                increment_histograms(t_output_data->data, state_with_impulse, state_data);
                for (p=0; p<probe_ids.size(); p++) increment_histograms(output_data->probes[probe_ids[p]]->data, state_with_impulse, state_data);
            }
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
#if (defined(USING_MPI))
//this is the binary function object used for the MPI reduce routine at the end of a 3F calculation
struct Reduce_Secondary_Dissociation_Results : public std::binary_function<Secondary_Dissociation_Result*, Secondary_Dissociation_Result*, Secondary_Dissociation_Result*>
{
    Secondary_Dissociation_Result* operator() (Secondary_Dissociation_Result *a, Secondary_Dissociation_Result *b){
        (*b)+=(*a);
        return b;
    }
};
namespace boost { namespace mpi {//this allows boost to use a more efficient gather algorithm, by telling boost that Reduce_Secondary_Dissociation_Results commutes
    template<>
    struct is_commutative<Reduce_Secondary_Dissociation_Results, Secondary_Dissociation_Result*> : mpl::true_{};
}}
#endif
//***********************************************************************************************************************************************************************************

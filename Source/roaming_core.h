//***********************************************************************************************************************************************************************************
class Roaming_Phasespace_Core : public Phasespace_Core
{
    public:
        Roaming_Result *r_output_data;
        const Roaming_Core_Input *r_input_data; //these are additional pointers to the Core_Input and Fragment objects pointed to by the pointers in the base class
        const Roaming_Molecule *r_parent, *r_frag1, *r_frag2; //however, these pointers are dynamically cast to be the derived class pointers so the extra functionality is accessible



        virtual ~Roaming_Phasespace_Core();
        Roaming_Phasespace_Core(){};
        Roaming_Phasespace_Core(const Roaming_Core_Input *in_input_data, const Roaming_Molecule *in_parent, const Roaming_Molecule *in_frag1, const Roaming_Molecule *in_frag2, Roaming_Result *in_output_data);
        void create_objects(Roaming_Core_Input *in_data, Roaming_Molecule *in_parent, Roaming_Molecule *in_frag1, Roaming_Molecule *in_frag2);
        virtual void create_input_pointer();
        virtual void create_output_pointer();
        virtual void create_fragment_pointers();
        virtual void initialise();
        virtual void process_state(const Double_Quantum_State &state);
        //virtual void process_state(const Double_Quantum_State &state);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Phasespace_Core::~Roaming_Phasespace_Core() {
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Phasespace_Core::Roaming_Phasespace_Core(const Roaming_Core_Input *in_input_data, const Roaming_Molecule *in_parent, const Roaming_Molecule *in_frag1, const Roaming_Molecule *in_frag2, Roaming_Result *in_output_data){
    copy_pointers(in_input_data, in_parent, in_frag1, in_frag2, in_output_data);
    create_input_pointer();
    create_output_pointer();
    create_fragment_pointers();
    initialise();
    run();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Phasespace_Core::create_output_pointer(){
    r_output_data = dynamic_cast<Roaming_Result*>(output_data); ///create a pointer to the output_data object that allows access to the derived fuctionality
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Phasespace_Core::create_input_pointer(){
    r_input_data = dynamic_cast<const Roaming_Core_Input*>(input_data); ///create a pointer to the input_data object that allows access to the derived fuctionality
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Phasespace_Core::create_fragment_pointers(){
    r_parent = dynamic_cast<const Roaming_Molecule*>(parent); ///create pointers to the fragment objects that allow access to the derived fuctionality
    r_frag1 = dynamic_cast<const Roaming_Molecule*>(frag1);
    r_frag2 = dynamic_cast<const Roaming_Molecule*>(frag2);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Phasespace_Core::initialise(){
    Phasespace_Core::initialise();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Phasespace_Core::process_state(const Double_Quantum_State &state){///this is called inside many nested loops ... efficiency here counts!
    unsigned p;
    std::vector<unsigned> probe_ids; //store the ids of all relevant probes
    Processed_State_Data state_data;

    for (p=0; p<r_output_data->r_probes.size(); p++) {
        if (r_output_data->r_probes[p]->state.applies(state)) {
            probe_ids.push_back(p);
            r_output_data->r_probes[p]->r_data->grand_total+=state.degeneracy;
        }
    }
    r_output_data->r_data->grand_total+=state.degeneracy;
    if (state.E_translation>=r_parent->delta_E_roam) { // dissociation will occur
        r_output_data->r_data->num_radical+=state.degeneracy;
        for (p=0; p<probe_ids.size(); p++) r_output_data->r_probes[probe_ids[p]]->r_data->num_radical+=state.degeneracy;

        if (r_input_data->output_histograms){
            state_data.E_trans=state.E_translation-r_parent->delta_E_roam; //at the potential of dissociation products (deltaEroam higher in energ) there will be deltaEroam less translational energy (E_avail is deltaEroam less for a traditional PST count)
            state_data.vel_f1=r_frag1->get_velocity(state_data.E_trans);
            state_data.vel_f2=r_frag2->get_velocity(state_data.E_trans);
            state_data.E_trans_f1=r_frag1->get_E_trans(state_data.vel_f1);
            state_data.E_trans_f2=r_frag2->get_E_trans(state_data.vel_f2);
            state_data.E_vib_f1=r_frag1->vib_levels[state.f1.v].energy;
            state_data.E_vib_f2=r_frag2->vib_levels[state.f2.v].energy;
            state_data.E_int_f1=state_data.E_vib_f1+state.f1.E_nk+r_frag1->spin_orbit_states[state.f1.so].energy; //E_total=E_vib+E_rot+E_soc
            state_data.E_int_f2=state_data.E_vib_f2+state.f2.E_nk+r_frag2->spin_orbit_states[state.f2.so].energy; //E_total=E_vib+E_rot+E_soc

            increment_histograms(r_output_data->r_data, state, state_data);
            for (p=0; p<probe_ids.size(); p++) increment_histograms(r_output_data->r_probes[probe_ids[p]]->r_data, state, state_data);
        }

    }else{ // may roam
        r_output_data->r_data->num_may_roam+=state.degeneracy;
        for (p=0; p<probe_ids.size(); p++) r_output_data->r_probes[probe_ids[p]]->r_data->num_may_roam+=state.degeneracy;
        if (r_input_data->output_histograms){
            if (!r_frag1->single_atom){//doubt we need the if here as the histgram won't increment if it isn't initialised, which the user shouldn't want and the code shouldn't allow ... revisit this! efficiency here counts!
                r_output_data->r_data->frag1_nk.increment_value(state.f1.n, state.f1.k, state.degeneracy);
            }
            if (!r_frag2->single_atom){
                r_output_data->r_data->frag2_nk.increment_value(state.f2.n, state.f2.k, state.degeneracy);
            }
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#define BOND_ANGLE_DISTRIBUTION_NONE 0
#define BOND_ANGLE_DISTRIBUTION_GROUND 1
#define BOND_ANGLE_DISTRIBUTION_PARENT 2
//***********************************************************************************************************************************************************************************
class Impulsive_Phasespace_Core : public Phasespace_Core
{
    public:
        const Impulsive_Core_Input *i_input_data; //these are additional pointers to the Core_Input and Fragment objects pointed to by the pointers in the base class
        const Impulsive_Molecule *i_parent, *i_frag1, *i_frag2; //however, these pointers are dynamically cast to be the derived class pointers so the extra functionality is accessible
        double E_recoil, tunneling_probability;
        std::vector<Impulsive_Bond_Angle_Distribution_Element> bond_angle_distribution;

        virtual ~Impulsive_Phasespace_Core();
        Impulsive_Phasespace_Core(){};
        Impulsive_Phasespace_Core(const Impulsive_Core_Input *in_input_data, const Impulsive_Molecule *in_parent, const Impulsive_Molecule *in_frag1, const Impulsive_Molecule *in_frag2, Phasespace_Result *in_output_data);
        virtual void create_input_pointer();
        virtual void create_fragment_pointers();
        virtual void initialise();
        double Eckart();
        virtual void run();
        virtual void process_state(const Double_Quantum_State &state);
        void apply_impulse(const Double_Quantum_State &state, Double_Quantum_State &state_with_impulse, double &E_wrongBy, const unsigned &i);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Impulsive_Phasespace_Core::~Impulsive_Phasespace_Core() {
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Impulsive_Phasespace_Core::Impulsive_Phasespace_Core(const Impulsive_Core_Input *in_input_data, const Impulsive_Molecule *in_parent, const Impulsive_Molecule *in_frag1, const Impulsive_Molecule *in_frag2, Phasespace_Result *in_output_data){
    copy_pointers(in_input_data, in_parent, in_frag1, in_frag2, in_output_data);
    create_input_pointer();
    create_fragment_pointers();
    initialise();
    run();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Phasespace_Core::create_input_pointer(){
    i_input_data = dynamic_cast<const Impulsive_Core_Input*>(input_data); ///create a pointer to the input_data object that allows access to the derived fuctionality
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Phasespace_Core::create_fragment_pointers(){
    i_parent = dynamic_cast<const Impulsive_Molecule*>(parent); ///create pointers to the fragment objects that allow access to the derived fuctionality
    i_frag1 = dynamic_cast<const Impulsive_Molecule*>(frag1);
    i_frag2 = dynamic_cast<const Impulsive_Molecule*>(frag2);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Phasespace_Core::initialise(){
    Phasespace_Core::initialise();
    Impulsive_Bond_Angle_Distribution_Element temp_angle_distribution_element;
    std::vector<Impulsive_Bond_Angle_Distribution_Element> compressed_bond_angle_distribution;
    unsigned i, j;
    bool element_found;
    double reduced_mass_kg, f1_i_param_m, f2_i_param_m, f1_m_inertia_kgsqm, f2_m_inertia_kgsqm, f1_x, f2_x, E_iTrans;
    if (i_parent->use_angle_distribution) bond_angle_distribution=i_input_data->bond_angle_distribution; //first, copy the bond angle distribution (or create single element version if not passing a distribution)
    else{//just use equilibrium angle with 100% probability for single element option
        temp_angle_distribution_element.probability=1;
        temp_angle_distribution_element.f1_radians=i_parent->TS_frag1_equilibrium_angle;//*RAD_IN_DEGREE; //default is 0 (if no barrier, impulsive_rotation_factor will not matter because multiplied by 0 barrier)
        temp_angle_distribution_element.f2_radians=i_parent->TS_frag2_equilibrium_angle;
        temp_angle_distribution_element.f1_mo_inertia=i_parent->TS_frag1_moment_of_inertia;
        temp_angle_distribution_element.f2_mo_inertia=i_parent->TS_frag2_moment_of_inertia;
        bond_angle_distribution.push_back(temp_angle_distribution_element);
    }// then, work out tunneling probability and recoil energy
    if (i_input_data->E_available <= 0) { //if not enough energy to get over barrier (we may want to be more robust here and ensure E_total>=E_dissociation, although should be eliminated in caller hopefully)
        if (i_input_data->tunneling_model==0) tunneling_probability=0; //no tunneling
        else if (i_input_data->tunneling_model==1) tunneling_probability=Eckart(); //use the Eckart tunneling model
        E_recoil=i_input_data->E_total-i_parent->E_dissociation; //recoil energy is energy over dissociation energy
    }else {//there is enough energy to get over barrier, no tunneling required
        tunneling_probability=1;
        E_recoil=i_parent->E_barrier;
    }// next, calculate E_rot for f1 and f2 at each bond angle (the bond angle distro is really only relevent for triatomic parents, and is therefore undocumanted. Future versions will account for all normal modes of TS and use metropolis sampling to describe how those modes change the bond angle and moment of inertia)
    if (i_parent->E_barrier){ //this if added because otherwise if no barrier, no TS, no i_parent->TS_moment_of_inertia, so impulsive_rotation_factor will become nan if this calculation occurs (impulsive_rotation_factor should just stay default of 0 if no barrier)
        reduced_mass_kg=get_reduced_mass(i_frag1->mass, i_frag2->mass)*KG_IN_AMU;
        //persistent_output<<"*** full distro ***"<<endl;
        if (i_input_data->output_to_file){
            if (!i_frag1->single_atom) fprintf(out_file.fp, "\n*** Impulsive Parameters Fragment 1 ***\nTS_frag1_equilibrium_angle=%lf\nTS_frag1_pivot_to_com_length=%lf\nTS_frag1_moment_of_inertia=%lf\nTS_frag1_kn_ratio=%lf\n",
                    i_parent->TS_frag1_equilibrium_angle, i_parent->TS_frag1_pivot_to_com_length, i_parent->TS_frag1_moment_of_inertia, i_parent->TS_frag1_kn_ratio);
            if (!i_frag2->single_atom) fprintf(out_file.fp, "\n*** Impulsive Parameters Fragment 2 ***\nTS_frag2_equilibrium_angle=%lf\nTS_frag2_pivot_to_com_length=%lf\nTS_frag2_moment_of_inertia=%lf\nTS_frag2_kn_ratio=%lf\n",
                    i_parent->TS_frag2_equilibrium_angle, i_parent->TS_frag2_pivot_to_com_length, i_parent->TS_frag2_moment_of_inertia, i_parent->TS_frag2_kn_ratio);
            fprintf(out_file.fp, "\n*** Impulsive Angle Distribution ***\n radians, probability, fragment, i_n, i_k, i_E_nk\n");
        }
        for (i=0; i<bond_angle_distribution.size(); i++){
            f1_x=f2_x=0;
            if (!i_frag1->single_atom){
                f1_m_inertia_kgsqm=bond_angle_distribution[i].f1_mo_inertia*KG_SQM_IN_AMU_SQA;
                f1_i_param_m=i_parent->TS_frag1_pivot_to_com_length*sin(bond_angle_distribution[i].f1_radians)*METERS_IN_ANGSTROM;
                f1_x=pow(f1_i_param_m,2)*reduced_mass_kg/f1_m_inertia_kgsqm;
            }
            if (!frag2->single_atom){
                f2_m_inertia_kgsqm=bond_angle_distribution[i].f2_mo_inertia*KG_SQM_IN_AMU_SQA;
                f2_i_param_m=i_parent->TS_frag2_pivot_to_com_length*sin(bond_angle_distribution[i].f2_radians)*METERS_IN_ANGSTROM;
                f2_x=pow(f2_i_param_m,2)*reduced_mass_kg/f2_m_inertia_kgsqm;
            }
            E_iTrans=E_recoil/(1+f1_x+f2_x);
            if (!i_frag1->single_atom){
                i_frag1->get_impulse_nk(E_iTrans*f1_x, i_parent->TS_frag1_kn_ratio, bond_angle_distribution[i].f1_n_imp, bond_angle_distribution[i].f1_k_imp, bond_angle_distribution[i].f1_E_nk_imp);
                if (i_input_data->output_to_file) fprintf(out_file.fp, "%lf, %lf, 1, %u, %u, %lf\n",
                                                          bond_angle_distribution[i].f1_radians, bond_angle_distribution[i].probability,
                                                          bond_angle_distribution[i].f1_n_imp, bond_angle_distribution[i].f1_k_imp, bond_angle_distribution[i].f1_E_nk_imp);
            }
            if (!i_frag2->single_atom){
                i_frag2->get_impulse_nk(E_iTrans*f2_x, i_parent->TS_frag2_kn_ratio, bond_angle_distribution[i].f2_n_imp, bond_angle_distribution[i].f2_k_imp, bond_angle_distribution[i].f2_E_nk_imp);
                if (i_input_data->output_to_file) fprintf(out_file.fp, "%lf, %lf, 2, %u, %u, %lf\n",
                                                          bond_angle_distribution[i].f2_radians, bond_angle_distribution[i].probability,
                                                          bond_angle_distribution[i].f2_n_imp, bond_angle_distribution[i].f2_k_imp, bond_angle_distribution[i].f2_E_nk_imp);
            }
            //if (i_input_data->E_total < E_potential) persistent_output<<i_input_data->E_total<<", "<<E_recoil<<", "<<E_iTrans<<", "<<tunneling_probability<<", "<<i_parent->initial_velocity<<endl;
            //persistent_output<<" f1_n_imp="<<bond_angle_distribution[i].f1_n_imp<<" f1_k_imp="<<bond_angle_distribution[i].f1_k_imp<<" f1_E_nk_imp="<<bond_angle_distribution[i].f1_E_nk_imp;
            //persistent_output<<" f2_n_imp="<<bond_angle_distribution[i].f2_n_imp<<" f2_k_imp="<<bond_angle_distribution[i].f2_k_imp<<" f2_E_nk_imp="<<bond_angle_distribution[i].f2_E_nk_imp<<" PROB="<<bond_angle_distribution[i].probability<<endl;
        }//then compress to only unique f1_nk_imp/f2_nk_imp combos to reduce calculation time
        for (i=0; i<bond_angle_distribution.size(); i++){
            element_found=0;
            for (j=0; j<compressed_bond_angle_distribution.size(); j++){
                if (bond_angle_distribution[i]==compressed_bond_angle_distribution[j]){
                    element_found=1;
                    compressed_bond_angle_distribution[j].probability+=bond_angle_distribution[i].probability;
                }
            }
            if (!element_found) compressed_bond_angle_distribution.push_back(bond_angle_distribution[i]);
        }
        bond_angle_distribution=compressed_bond_angle_distribution;
        /*persistent_output<<"*** compressed distro ***"<<endl;
        for (i=0; i<bond_angle_distribution.size(); i++){
            persistent_output<<" f1_n_imp="<<bond_angle_distribution[i].f1_n_imp<<" f1_k_imp="<<bond_angle_distribution[i].f1_k_imp<<" f1_E_nk_imp="<<bond_angle_distribution[i].f1_E_nk_imp;
            persistent_output<<" f2_n_imp="<<bond_angle_distribution[i].f2_n_imp<<" f2_k_imp="<<bond_angle_distribution[i].f2_k_imp<<" f2_E_nk_imp="<<bond_angle_distribution[i].f2_E_nk_imp<<" PROB="<<bond_angle_distribution[i].probability<<endl;
        }*/
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Impulsive_Phasespace_Core::Eckart(){
    /*
    The following function is based on the script Eckart.py, which is part of the CanTherm program.
    It was found at this URL: https://github.com/GreenGroup/CanTherm_old/blob/master/source/Eckart.py
    The licence for that code is reproduced below:

    Copyright (c) 2002-2009 William H. Green and the CanTherm Team

    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
    the Software, and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    */
    double a1, a2, twoPIa, twoPIb, twoPId, xi, nuStar, kappa_E, numerator, denominator;
    double FStar = i_parent->TS_reaction_coord_force_constant*NEWTONS_PER_METER_IN_MDYNE_PER_ANGSTROM;//this is the force constant as calculated by gaussian for the TS (in mdyne/A) converted to SI units (N/m)
    nuStar=(1/(2*PI))*sqrt(FStar/(i_parent->TS_reaction_coord_reduced_mass*KG_IN_AMU));
    a1=2*PI*(i_parent->E_dissociation+i_parent->E_barrier)*JOULES_IN_WAVENUMBER/(PLANCKS_CONSTANT*nuStar); //dimensionless energy difference between TS and reactants
    a2=2*PI*i_parent->E_barrier*JOULES_IN_WAVENUMBER/(PLANCKS_CONSTANT*nuStar); //dimensionless energy difference between TS and products
    xi=i_input_data->E_total/(i_parent->E_dissociation+i_parent->E_barrier);
    twoPIa=2*sqrt(a1*xi)/(1/sqrt(a1)+1/sqrt(a2));
    twoPIb=2*sqrt(abs((xi-1)*a1+a2))/(1/sqrt(a1)+1/sqrt(a2)); //diff from one ref ... this is how code did ... corrected in another ref
    twoPId=2*sqrt(abs(a1*a2-4*pow(PI,2)/16));
    if (twoPIa<200 && twoPIb<200 && twoPId<200){ //all cosh arguments are small enough to use the Eckart equation for tunneling probability
        kappa_E = 1-(cosh(twoPIa-twoPIb)+cosh(twoPId)) / (cosh(twoPIa+twoPIb)+cosh(twoPId));
    }else{//if at least one of the following expressions is greater than 5, we can eliminate most of the exponential terms after writing out the definition of cosh() and dividing all terms by exp(twopid)
        if ((twoPIa-twoPIb-twoPId > 10) || (twoPIb-twoPIa-twoPId > 10) || (twoPIa+twoPIb-twoPId > 10)){
            kappa_E = 1 - exp(-2*twoPIa) - exp(-2*twoPIb) - exp(-twoPIa-twoPIb+twoPId) - exp(-twoPIa-twoPIb-twoPId);
        }else{ //If all of the arguments are less than 5, then evaluate the kappa_E expression normally, except use the expanded definition - expanding the cosh argument and dividing all terms by exp(twopid)
            numerator = exp(twoPIa-twoPIb-twoPId)+exp(-twoPIa+twoPIb-twoPId)+1+exp(-2*twoPId);
            denominator = exp(twoPIa+twoPIb-twoPId)+exp(-twoPIa-twoPIb-twoPId)+1+exp(-2*twoPId);
            kappa_E = 1 - numerator/denominator;
        }
    }
    return kappa_E;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Phasespace_Core::run(){
    Double_Quantum_State tunneling_state;
    tunneling_state.tunneled=1;
    tunneling_state.f1.n=tunneling_state.f2.n=0;
    tunneling_state.f1.k=tunneling_state.f2.k=0;
    tunneling_state.f1.E_nk=tunneling_state.f2.E_nk=0;
    if (input_data->E_available>0){ //if there is energy available for statistical distribution
        count_states();
        if (i_input_data->output_to_file) write_result(out_file.fp);
    }else if (tunneling_probability){ //only tunneling states, which are fully impulsive
        tunneling_state.degeneracy=i_frag1->get_degeneracy(0,0,0)+i_frag2->get_degeneracy(0,0,0);
        process_state(tunneling_state);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Phasespace_Core::apply_impulse(const Double_Quantum_State &state, Double_Quantum_State &state_with_impulse, double &E_wrongBy, const unsigned &i){
        state_with_impulse=state;
        E_wrongBy=0;
        if (i_input_data->impulsive_J_mode==0){ // the old add vectors end to end version ... I think this will be removed eventually as it's probably the wrongest way
            if (!i_frag1->single_atom){
                state_with_impulse.f1.n+=bond_angle_distribution[i].f1_n_imp;//add the impulsive angular momentum
                state_with_impulse.f1.k+=bond_angle_distribution[i].f1_k_imp;//add the impulsive angular momentum component on principal axis
                state_with_impulse.f1.E_nk=i_frag1->get_E_nk(state_with_impulse.f1.n, state_with_impulse.f1.k);
                E_wrongBy+=state_with_impulse.f1.E_nk-state.f1.E_nk-bond_angle_distribution[i].f1_E_nk_imp; //calculate the energy difference between the impulsive and PST angular momentum seperately and together (for correction below to conserve energy
            }
            if (!i_frag2->single_atom){
                state_with_impulse.f2.n+=bond_angle_distribution[i].f2_n_imp;//add the impulsive angular momentum
                state_with_impulse.f2.k+=bond_angle_distribution[i].f2_k_imp;//add the impulsive angular momentum component on principal axis
                state_with_impulse.f2.E_nk=i_frag2->get_E_nk(state_with_impulse.f2.n, state_with_impulse.f2.k);
                E_wrongBy+=state_with_impulse.f2.E_nk-state.f2.E_nk-bond_angle_distribution[i].f2_E_nk_imp;
            }
            state_with_impulse.E_translation+=(E_recoil-bond_angle_distribution[i].f1_E_nk_imp-bond_angle_distribution[i].f2_E_nk_imp-E_wrongBy);
        }else if (i_input_data->impulsive_J_mode==1){ // add in quadrature (default) ... other option, distribution of angles over 180 degrees, will require a loop so the calling function will have to implement the if statement
            if (!i_frag1->single_atom){ //another option would be to add the impulsive energy to each state here (IE find the correct NK state here then add, I think that's probably to go!)
                state_with_impulse.f1.n=sqrt(pow(state.f1.n, 2) + pow(bond_angle_distribution[i].f1_n_imp, 2));//add the impulsive angular momentum
                state_with_impulse.f1.k=sqrt(pow(state.f1.k, 2) + pow(bond_angle_distribution[i].f1_k_imp, 2));//add the impulsive angular momentum component on principal axis
                state_with_impulse.f1.E_nk=i_frag1->get_E_nk(state_with_impulse.f1.n, state_with_impulse.f1.k);
                E_wrongBy+=state_with_impulse.f1.E_nk-state.f1.E_nk-bond_angle_distribution[i].f1_E_nk_imp; //calculate the energy difference between the impulsive and PST angular momentum seperately and together (for correction below to conserve energy
            }
            if (!i_frag2->single_atom){
                state_with_impulse.f2.n=(unsigned)sqrt(pow(state.f2.n, 2) + pow(bond_angle_distribution[i].f2_n_imp, 2));//add the impulsive angular momentum
                state_with_impulse.f2.k=(unsigned)sqrt(pow(state.f2.k, 2) + pow(bond_angle_distribution[i].f2_k_imp, 2));//add the impulsive angular momentum component on principal axis
                state_with_impulse.f2.E_nk=i_frag2->get_E_nk(state_with_impulse.f2.n, state_with_impulse.f2.k);
                E_wrongBy+=state_with_impulse.f2.E_nk-state.f2.E_nk-bond_angle_distribution[i].f2_E_nk_imp;
            }
            state_with_impulse.E_translation+=(E_recoil-bond_angle_distribution[i].f1_E_nk_imp-bond_angle_distribution[i].f2_E_nk_imp-E_wrongBy);
        }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Phasespace_Core::process_state(const Double_Quantum_State &state){
    Double_Quantum_State state_with_impulse;
    Processed_State_Data state_data;
    double E_wrongBy;
    unsigned i, p;
    std::vector<unsigned> probe_ids; //store the ids of all relevant probes
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

            if (i_input_data->output_histograms){

                //are these generic?
                state_data.E_trans=state.E_translation;
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
            }
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************

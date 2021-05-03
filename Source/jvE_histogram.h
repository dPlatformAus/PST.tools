struct jvE_parent_state
{
    unsigned j;
    double energy, jE, vE, vBoltzmann_prob, jBoltzmann_prob, j_prob, jvE_prob; ///don't really need to store jE, vE, vBoltzmann_prob, jBoltzmann_prob or j_prob, these are only for testing! get rid of them!
};
bool jvE_parent_state_sort_predicate(const jvE_parent_state &lhs, const jvE_parent_state &rhs) {return lhs.jvE_prob > rhs.jvE_prob;}
bool jvE_parent_state_sort_energy_desc(const jvE_parent_state &lhs, const jvE_parent_state &rhs) {return lhs.energy > rhs.energy;}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class jvE_histogram
{
    public:
        Histogram vE_hist, jE_hist;
        std::vector<Histogram> j_hist;
        std::vector<jvE_parent_state> jvE_hist;
        Output_File out_file;
        double e_energy_within, e_energy_outside, e_prob_within, e_prob_outside, p_energy_within, p_energy_outside, p_prob_within, p_prob_outside, get_min_energy_within_top_percent, get_percent_within_min_energy, final_prob_included;

        ~jvE_histogram();
        jvE_histogram(const std::string& out_filename, const std::string &job_path, const General_Input *the_general_input, Phasespace_Molecule *in_parent);
        unsigned size() {return jvE_hist.size();}
        void incrementJ(const double in_x, const double in_y, const unsigned in_j);
        void set_0K_jE_hist(const General_Input *the_general_input);
        void set_0K_vE_hist(const General_Input *the_general_input);
        void initialiseJ(const General_Input *the_general_input, Phasespace_Molecule *in_parent);
        void initialiseV(const General_Input *the_general_input, Phasespace_Molecule *in_parent);
        void initialiseJV(const General_Input *the_general_input, Phasespace_Molecule *in_parent);
        void write_histograms(FILE *the_file);
        void get_min_energy_within_specified_probability(const General_Input *the_general_input);
        void get_probability_within_specified_min_energy(const General_Input *the_general_input);
        void sort_jvE_hist() {std::sort(jvE_hist.begin(), jvE_hist.end(), jvE_parent_state_sort_predicate);}

};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
jvE_histogram::~jvE_histogram() {
	out_file.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
jvE_histogram::jvE_histogram(const std::string& out_filename, const std::string &job_path, const General_Input *the_general_input, Phasespace_Molecule *in_parent) {
    unsigned i, k, l, index; //no j counter because its too hard to read with all the .j in there too =)
    jvE_parent_state temp_jvE_parent_state;
    double jvE_prob;
    if (the_general_input->equilibrium_temerature == 0) {
        if (the_general_input->rotational_temerature == 0) set_0K_jE_hist(the_general_input); //then we just need one jE state 0,0
        else initialiseJ(the_general_input, in_parent);
        if (the_general_input->vibrational_temerature == 0) set_0K_vE_hist(the_general_input); //put the v=0 state in there
        else initialiseV(the_general_input, in_parent);
    } else initialiseJV(the_general_input, in_parent);
    index=0;
    for (l=0; l<vE_hist.bin_length; l++){
        for (i=0; i<jE_hist.bin_length; i++){
            j_hist[i].get_average(); //will print the average in jvE ouput file regardless of the_general_input->Boltzmann_j_averaging
            if (the_general_input->Boltzmann_j_averaging) { //only add a single entry for this jv bin combination, using the average Jparent value for the j bin
                if (j_hist[i].data_average) {//only store j state if it has non zero probability
                    jvE_prob=vE_hist.get_Boltzmann_prob(l)*jE_hist.get_Boltzmann_prob(i);
                    if (jvE_prob>0) { ///only store if nonzero probability
                        jvE_hist.push_back(temp_jvE_parent_state);
                        jvE_hist[index].jvE_prob=jvE_prob;
                        jvE_hist[index].energy=jE_hist.resolveIndex(i)+vE_hist.resolveIndex(l);
                        jvE_hist[index].j=boost::math::iround(j_hist[i].data_average);

                        ///the 4 values below are only for testing! get rid of them!
                        jvE_hist[index].jE=jE_hist.resolveIndex(i);
                        jvE_hist[index].vE=vE_hist.resolveIndex(l);
                        jvE_hist[index].vBoltzmann_prob=vE_hist.get_Boltzmann_prob(l);
                        jvE_hist[index].jBoltzmann_prob=jE_hist.get_Boltzmann_prob(i);

                        index++;
                    }
                }
            }else{ //not parent j averaging, store an entry for each Jparent value for this jv bin combination
                for (k=0; k<j_hist[i].bin_length; k++){
                    if (j_hist[i].data[k]) { //only store j state if it has non zero probability
                        jvE_prob=vE_hist.get_Boltzmann_prob(l)*jE_hist.get_Boltzmann_prob(i)*j_hist[i].get_normalised_area_prob(k);
                        if (jvE_prob>0) { ///only store if nonzero probability
                            jvE_hist.push_back(temp_jvE_parent_state);
                            jvE_hist[index].jvE_prob=jvE_prob;
                            jvE_hist[index].energy=jE_hist.resolveIndex(i)+vE_hist.resolveIndex(l);
                            jvE_hist[index].j=(int)j_hist[i].resolveIndex(k);

                            ///the 5 values below are only for testing! get rid of them!
                            jvE_hist[index].jE=jE_hist.resolveIndex(i);
                            jvE_hist[index].vE=vE_hist.resolveIndex(l);
                            jvE_hist[index].vBoltzmann_prob=vE_hist.get_Boltzmann_prob(l);
                            jvE_hist[index].jBoltzmann_prob=jE_hist.get_Boltzmann_prob(i);
                            jvE_hist[index].j_prob=j_hist[i].get_normalised_area_prob(k);

                            index++;
                        }
                    }
                }
            }
        }
    }
    get_percent_within_min_energy=get_min_energy_within_top_percent=p_energy_within=p_energy_outside=p_prob_within=p_prob_outside=e_energy_within=e_energy_outside=e_prob_within=e_prob_outside=0;
    if (the_general_input->get_percent_within_min_energy) get_probability_within_specified_min_energy(the_general_input);
    if (the_general_input->get_min_energy_within_top_percent) get_min_energy_within_specified_probability(the_general_input); // before we sort by prob for the output because it requires a sort by energy decending!
    sort_jvE_hist();
    if (the_general_input->jvE_histogram_write_to_file) {
        out_file.set_output_directory(job_path+"jvE_histogram/");
        out_file.open(out_filename);
        write_histograms(out_file.fp);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void jvE_histogram::set_0K_jE_hist(const General_Input *the_general_input) {//populate jE_hist(one bin E=0 prob=1) and j_hist(one bin j=0 prob=1)
    Histogram temp_j_hist;
    jE_hist.initialiseS("Binned Rotational Energies", "Energy", "Degeneracy", the_general_input->Boltzmann_bin_size, 0);
    jE_hist.increment_value(0,1);
    temp_j_hist.initialiseS("Binned j for energy bin: 0", "j", "Count", 1, 0);
    temp_j_hist.increment_value(0,1);
    j_hist.push_back(temp_j_hist);
    jE_hist.Z_component.push_back(1);///Can't call set_Boltzmann_temp(0) ... will divide by 0! ... so artificially populate Z_component and Z
    jE_hist.Z=1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void jvE_histogram::set_0K_vE_hist(const General_Input *the_general_input) {//populate vE_hist(one bin E=0 prob=1)
    vE_hist.initialiseS("Binned Vibrational Energies", "Energy", "Degeneracy", the_general_input->Boltzmann_bin_size, 0);
    vE_hist.increment_value(0,1);
    vE_hist.Z_component.push_back(1);///Can't call set_Boltzmann_temp(0) ... will divide by 0! ... so artificially populate Z_component and Z
    vE_hist.Z=1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void jvE_histogram::incrementJ(const double in_x, const double in_y, const unsigned in_j){
    jE_hist.increment_value(in_x, in_y);
    j_hist[jE_hist.getIndex(in_x)].increment_value(in_j, in_y);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void jvE_histogram::initialiseJ(const General_Input *the_general_input, Phasespace_Molecule *in_parent) {
    double E_nk, max_rotational_energy;
    unsigned i, n, n_max, k, to_k, degeneracy;
    max_rotational_energy=the_general_input->Boltzmann_max_energy_factor*the_general_input->rotational_temerature*WAVENUMBER_IN_KELVIN;
    Phasespace_Molecule *the_parent;
    the_parent=in_parent->clone();
    the_parent->initialise(0);
    the_parent->initialise_vib_levels(0);
    n_max=the_parent->get_n_max(max_rotational_energy);
    Histogram temp_j_hist;
    temp_j_hist.initialiseS("Binned j for energy bin: ", "j", "Count", 1, n_max);
    jE_hist.initialiseS("Binned Rotational Energies", "Energy", "Degeneracy", the_general_input->Boltzmann_bin_size, max_rotational_energy);
    for (i=0;i<jE_hist.bin_length;i++) {
        j_hist.push_back(temp_j_hist);
        j_hist[i].title.append(boost::lexical_cast<std::string>(jE_hist.resolveIndex(i)));
    }
    for (n=0;n<=n_max;n++){
        if (the_parent->linear) to_k=0;
        else to_k=n;
        for (k=0;k<=to_k;k++){
            E_nk=the_parent->get_E_nk(n, k);
            if (E_nk<=max_rotational_energy){
                degeneracy=2*n+1; // all J states have a degeneracy of 2J+1  *** this doesn't allow triangle inequality or centrifugal barrier as we don't interrogate L ... perhaps we need to do a full phasespace calc considering both fragments from N2O5 pyrolysis (but what if we have hot molecules that are not dissociation products?)
                if (k) degeneracy*=2; //double the degeneracy for non zero k, as k can be plus or minus
                incrementJ(E_nk, degeneracy, n);
            }
        }
    }
    delete the_parent;
    jE_hist.set_Boltzmann_temp(the_general_input->rotational_temerature);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void jvE_histogram::initialiseJV(const General_Input *the_general_input, Phasespace_Molecule *in_parent) {
    double max_energy, vib_energy, E_nk;
    unsigned i, n, n_max, k, to_k, v_degeneracy, j_degeneracy;
    max_energy=the_general_input->Boltzmann_max_energy_factor*the_general_input->equilibrium_temerature*WAVENUMBER_IN_KELVIN;
    Phasespace_Molecule *the_parent;
    the_parent=in_parent->clone();
    the_parent->initialise(0);
    the_parent->initialise_vib_levels(max_energy);
    Histogram temp_j_hist;
    n_max=the_parent->get_n_max(max_energy);
    temp_j_hist.initialiseS("Binned j for energy bin: ", "j", "Count", 1, n_max);
    jE_hist.initialiseS("Binned Rovibrational Energies", "Energy", "Degeneracy", the_general_input->Boltzmann_bin_size, max_energy);
    for (i=0;i<jE_hist.bin_length;i++) {
        j_hist.push_back(temp_j_hist);
        j_hist[i].title.append(boost::lexical_cast<std::string>(jE_hist.resolveIndex(i)));
    }
    for (i=0;i<the_parent->vib_levels.size(); i++) {
        vib_energy=the_parent->vib_levels[i].energy;
        v_degeneracy=the_parent->vib_levels[i].degeneracy;
        n_max=the_parent->get_n_max(max_energy-vib_energy);
        for (n=0;n<=n_max;n++){
            if (the_parent->linear) to_k=0;
            else to_k=n;
            for (k=0;k<=to_k;k++){
                E_nk=the_parent->get_E_nk(n, k);
                if ((vib_energy+E_nk)<=max_energy){
                    j_degeneracy=2*n+1; // all J states have a degeneracy of 2J+1
                    if (k) j_degeneracy*=2; //double the degeneracy for non zero k, as k can be plus or minus
                    incrementJ(vib_energy+E_nk, v_degeneracy*j_degeneracy, n);
                }
            }
        }
    }
    delete the_parent;
    jE_hist.set_Boltzmann_temp(the_general_input->equilibrium_temerature);
    set_0K_vE_hist(the_general_input); //all goes into jE_hist but doing it for full rovibrational state count
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void jvE_histogram::initialiseV(const General_Input *the_general_input, Phasespace_Molecule *in_parent) {
    double max_vibrational_energy;
    unsigned i;
    max_vibrational_energy=the_general_input->Boltzmann_max_energy_factor*the_general_input->vibrational_temerature*WAVENUMBER_IN_KELVIN;
    Phasespace_Molecule *the_parent;
    the_parent=in_parent->clone();
    the_parent->initialise(0);
    the_parent->initialise_vib_levels(max_vibrational_energy);
    vE_hist.initialiseS("Binned Vibrational Energies", "Energy", "Degeneracy", the_general_input->Boltzmann_bin_size, max_vibrational_energy);
    for (i=0;i<the_parent->vib_levels.size(); i++) {
        vE_hist.increment_value(the_parent->vib_levels[i].energy, the_parent->vib_levels[i].degeneracy);
    }
    delete the_parent;
    vE_hist.set_Boltzmann_temp(the_general_input->vibrational_temerature);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void jvE_histogram::get_min_energy_within_specified_probability(const General_Input *the_general_input){
    unsigned i=0;
    double current_prob_included=0;
    get_min_energy_within_top_percent=the_general_input->get_min_energy_within_top_percent;
    std::sort(jvE_hist.begin(), jvE_hist.end(), jvE_parent_state_sort_energy_desc); //sort be energy descending
    while (i<jvE_hist.size() && (current_prob_included*100.0<=get_min_energy_within_top_percent || get_min_energy_within_top_percent==100)) {
        p_energy_within=jvE_hist[i].energy;
        p_prob_within=current_prob_included;
        current_prob_included+=jvE_hist[i].jvE_prob;
        i++;
    }
    p_prob_outside=current_prob_included;
    final_prob_included=current_prob_included;
    if (i<jvE_hist.size()) p_energy_outside=jvE_hist[i].energy;
    else p_energy_outside=p_energy_within; //we are at the end, they asked for 100% (pointless)
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void jvE_histogram::get_probability_within_specified_min_energy(const General_Input *the_general_input){
    unsigned i=0;
    double current_prob_included=0;
    get_percent_within_min_energy=the_general_input->get_percent_within_min_energy;
    std::sort(jvE_hist.begin(), jvE_hist.end(), jvE_parent_state_sort_energy_desc); //sort be energy descending
    while (i<jvE_hist.size() && jvE_hist[i].energy>get_percent_within_min_energy) {
        e_energy_within=jvE_hist[i].energy;
        e_prob_within=current_prob_included;
        current_prob_included+=jvE_hist[i].jvE_prob;
        i++;
    }
    e_prob_outside=current_prob_included;
    if (i<jvE_hist.size()) e_energy_outside=jvE_hist[i].energy;
    else e_energy_outside=e_energy_within; //we are at the end, they asked for 100% (pointless)
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void jvE_histogram::write_histograms(FILE *the_file){
    unsigned i;
    if (get_min_energy_within_top_percent) fprintf(the_file, "\nLowest energy within top %lf%% = %lf\nHighest energy outside top %lf%% = %lf\nFirst probability above target=%lf%%\n", get_min_energy_within_top_percent, p_energy_within, get_min_energy_within_top_percent, p_energy_outside, final_prob_included);
    fprintf(the_file, "\nThe jvE histogram\nenergy, j, jE, vE, vBoltzmann_prob, jBoltzmann_prob, j_prob, jvE_prob\n");
    for (i=0; i<jvE_hist.size(); i++) {
        fprintf(the_file, "%lf, %d, %lf, %lf, %lf, %lf, %lf, %.17lf \n", jvE_hist[i].energy, jvE_hist[i].j, jvE_hist[i].jE, jvE_hist[i].vE, jvE_hist[i].vBoltzmann_prob, jvE_hist[i].jBoltzmann_prob, jvE_hist[i].j_prob, jvE_hist[i].jvE_prob);
    }
    vE_hist.write_Boltzmann(the_file);
    jE_hist.write_Boltzmann(the_file);
    for (i=0; i<jE_hist.bin_length; i++) {
        j_hist[i].write_histogram(the_file);
        j_hist[i].write_average(the_file);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

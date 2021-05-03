//***********************************************************************************************************************************************************************************
class Phasespace_Get_Rotational_States_Calculation : public Phasespace_Calculation
{
    public:
        virtual ~Phasespace_Get_Rotational_States_Calculation(){};
        virtual bool run();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Phasespace_Get_Rotational_States_Calculation::run(){ ///gets a rigid rotor rotational state list, up to an input energy, for input rotational constants
    double total_energy, E_nk;
    unsigned n, n_max, k, to_k, degeneracy, num_states;
    persistent_output << " - Get Rotational States:" << std::endl;
    reactant->initialise(0);
    for (total_energy=the_general_input->excitation_energy_range.min; total_energy<=the_general_input->excitation_energy_range.max; total_energy+=the_general_input->excitation_energy_increment){
        clear_semipersistent_output();
        semipersistent_output << "total_energy (" << the_general_input->excitation_energy_range.min << "-" << the_general_input->excitation_energy_range.max << ") = " << total_energy << std::endl;
        print_semipersistent_output();
        fprintf(out_file.fp, "\nat total_energy=%lf there are the following Rotational Energy Levels\nN,K,energy (cm-1),degeneracy\n", total_energy);
        num_states=0;
        n_max=reactant->get_n_max(total_energy);
        for (n=0;n<=n_max;n++){
            if (reactant->linear) to_k=0;
            else to_k=n;
            for (k=0;k<=to_k;k++){
                E_nk=reactant->get_E_nk(n, k);
                if (E_nk<=total_energy){
                    num_states++;
                    degeneracy=2*n+1; // all J states have a degeneracy of 2J+1  *** this doesn't allow triangle inequality or centrifugal barrier as we don't interrogate L ... perhaps we need to do a full phasespace calc considering both fragments from N2O5 pyrolysis (but what if we have hot molecules that are not dissociation products?)
                    if (k) degeneracy*=2; //double the degeneracy for non zero k, as k can be plus or minus
                    fprintf(out_file.fp, "%d,%d,%lf,%d\n",n,k,E_nk,degeneracy);
                }
            }
        }
        fprintf(out_file.fp, "\nat total_energy=%lf the total number of Rotational Energy Levels was %d\n", total_energy, num_states);
    }
    return true;
}
//***********************************************************************************************************************************************************************************
class Phasespace_Get_Vibrational_States_Calculation : public Phasespace_Calculation
{
    public:
        Phasespace_Get_Vibrational_States_Calculation(const std::string &file_name){go(file_name);};
        virtual ~Phasespace_Get_Vibrational_States_Calculation(){};
        virtual bool run();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Phasespace_Get_Vibrational_States_Calculation::run(){ ///gets a rigid rotor rotational state list, up to an input energy, for input rotational constants
    double total_energy;
    unsigned i, num_states;
    persistent_output << " - Get Vibrational States:" << std::endl;
    reactant->initialise_vib_levels(the_general_input->excitation_energy_range.max);
    for (total_energy=the_general_input->excitation_energy_range.min; total_energy<=the_general_input->excitation_energy_range.max; total_energy+=the_general_input->excitation_energy_increment){
        clear_semipersistent_output();
        semipersistent_output << "total_energy (" << the_general_input->excitation_energy_range.min << "-" << the_general_input->excitation_energy_range.max << ") = " << total_energy << std::endl;
        print_semipersistent_output();
        fprintf(out_file.fp, "\nat total_energy=%lf there are the following Vibrational Energy Levels\ncomponents,energy (cm-1),degeneracy\n", total_energy);
        num_states=0;
        for (i=0;i<reactant->vib_levels.size(); i++){
            if (reactant->vib_levels[i].energy<=total_energy) {
                fprintf(out_file.fp, "%s,%lf,%d\n",reactant->vib_levels[i].name.c_str(),reactant->vib_levels[i].energy, reactant->vib_levels[i].degeneracy);
                num_states++;
            }
        }
        fprintf(out_file.fp, "\nat total_energy=%lf the total number of Vibrational Energy Levels was %d\n", total_energy, num_states);
    }
    return true;
}
//***********************************************************************************************************************************************************************************




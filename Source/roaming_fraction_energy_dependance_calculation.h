//***********************************************************************************************************************************************************************************
class Roaming_Fraction_Energy_Dependance_Calculation : public Roaming_Calculation
{
    public:
        virtual ~Roaming_Fraction_Energy_Dependance_Calculation(){};
        virtual bool run();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Roaming_Fraction_Energy_Dependance_Calculation::run(){
    double total_energy;
    unsigned i, p;

    Roaming_Core_Input *the_pst_input;
    Roaming_Phasespace_Core *the_pst;
    Roaming_Result *the_PST_output; //stores the output from a single roaming PST count
    Roaming_Fraction_Result the_roaming_fraction; //stores the overall fractions, for a single photolysis energy, for however may PST counts are required to acquire them (which could be many with parent internal energy considered)
    std::vector<Roaming_Fraction_Result> probe_results(r_reactant->probe_states.size()); //eventually the overall result will be done with one of these too

    r_general_input->histogram_write_to_file=0; ///We force no histogram output for this calculation type, for roaming histograms, run roaming_calculation

    the_pst_input = new Roaming_Core_Input();
    the_PST_output = new Roaming_Result();
    //these three can probably be in a function set_standard_PST_input_values() or something
    the_pst_input->job_path = job_path;
    the_pst_input->output_to_file = r_general_input->objects_write_to_file;
    the_pst_input->output_histograms = r_general_input->histogram_write_to_file;
    //these will be manipulated differently by different calculations
    the_pst_input->job_name="Roaming_Fraction_Energy_Dependance_Calculation";
    the_PST_output->initialise(limits, r_reactant); //no setting limits here ... for histograms, use a different calculation I think ... this one outputs only totals
    the_roaming_fraction.initialise(out_file); //initialise the overall Roaming_Fraction_Result so that it writes to the main calculation output file.
    for(p=0; p<r_reactant->probe_states.size(); p++) probe_results[p].initialise(job_path, r_reactant->probe_states[p].name);
    persistent_output << " - Calculate Roaming Fraction Energy Dependence:" << std::endl;
    the_roaming_fraction.write_header();
    for(p=0; p<probe_results.size(); p++) probe_results[p].write_header();
    for (total_energy=r_general_input->excitation_energy_range.min; total_energy<=r_general_input->excitation_energy_range.max; total_energy+=r_general_input->excitation_energy_increment){
        clear_semipersistent_output();
        semipersistent_output << "excitation_energy (" << r_general_input->excitation_energy_range.min << "-" << r_general_input->excitation_energy_range.max << ") = " << total_energy << std::endl;
        print_semipersistent_output();
        the_progress_bar.initialise(100,1,jvE_i_max, CALCULATION_PROGRESS_BAR); // initialise the progress bar (this has no effect on the calculation)
        the_roaming_fraction.reset();
        for(p=0; p<probe_results.size(); p++) probe_results[p].reset();
        i=0;
        while (i<jvE_i_max) {
            the_progress_bar.go(i, CALCULATION_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
            r_reactant->initial_J=Erovib_histogram->jvE_hist[i].j;
            the_pst_input->E_total=total_energy+Erovib_histogram->jvE_hist[i].energy;
            the_pst_input->E_available=the_pst_input->E_total-r_reactant->E_dissociation+r_reactant->delta_E_roam;
            the_pst_input->job_name="PST_"+boost::lexical_cast<std::string>(r_reactant->delta_E_roam)+"_"+boost::lexical_cast<std::string>(total_energy);
            the_PST_output->reset();
            the_pst=new Roaming_Phasespace_Core(the_pst_input,r_reactant,r_frag1,r_frag2,the_PST_output);
            delete the_pst;
            the_roaming_fraction.i_reset();
            for(p=0; p<probe_results.size(); p++) probe_results[p].i_reset();
            if (the_PST_output->r_data->grand_total) {
                the_roaming_fraction.i_increment(the_PST_output->r_data->num_may_roam, the_PST_output->r_data->num_radical, the_PST_output->r_data->grand_total, r_reactant->P_roam, Erovib_histogram->jvE_hist[i].jvE_prob);
                for(p=0; p<probe_results.size(); p++) probe_results[p].i_increment(the_PST_output->r_probes[p]->r_data->num_may_roam, the_PST_output->r_probes[p]->r_data->num_radical, the_PST_output->r_probes[p]->r_data->grand_total, r_reactant->P_roam, Erovib_histogram->jvE_hist[i].jvE_prob);
            }
            i++;
        }
        the_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
        the_roaming_fraction.write_result(total_energy);
        for(p=0; p<probe_results.size(); p++) probe_results[p].write_result(total_energy);
    }

    delete the_pst_input;
    delete the_PST_output;
    return true;
}
//***********************************************************************************************************************************************************************************



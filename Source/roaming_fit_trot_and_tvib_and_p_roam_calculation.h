//***********************************************************************************************************************************************************************************
class Roaming_Fit_Trot_And_Tvib_And_P_Roam_Calculation : public Roaming_Calculation
{
    public:
        virtual ~Roaming_Fit_Trot_And_Tvib_And_P_Roam_Calculation(){};
        virtual void jvE_initialise();
        virtual bool run();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Fit_Trot_And_Tvib_And_P_Roam_Calculation::jvE_initialise(){
    std::string jvE_file_name;
    double current_prob_included;
    parent_Erovib_max=0;
    jvE_file_name="jvE_histogram.csv";
    r_general_input->rotational_temerature=r_general_input->rotational_temerature_range.max; //we want jvE_i_max to be the maximum so that fragmet vib levels and limits (if we implement histograms in ths calc) are initialised correctly
    r_general_input->vibrational_temerature=r_general_input->vibrational_temerature_range.max;
    Erovib_histogram=new jvE_histogram(jvE_file_name, job_path, r_general_input, reactant);
    current_prob_included=jvE_i_max=0;//calculate how far down the jvE_histogram to go to get to r_general_input->Boltzmann_probability_included, doing this first means we can define the progress bar better and also means the inner loop for each fit point is slightly more efficient
    while (jvE_i_max<Erovib_histogram->size() && (current_prob_included*100<=r_general_input->Boltzmann_probability_included || r_general_input->Boltzmann_probability_included==100)) {
        current_prob_included+=Erovib_histogram->jvE_hist[jvE_i_max].jvE_prob;
        if (parent_Erovib_max<Erovib_histogram->jvE_hist[jvE_i_max].energy) parent_Erovib_max=Erovib_histogram->jvE_hist[jvE_i_max].energy;
        jvE_i_max++;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Roaming_Fit_Trot_And_Tvib_And_P_Roam_Calculation::run(){
    double Trot, Tvib, num_roam, num_rad, num_total, f_roam, f_rad, if_roam, if_rad, bwf_roam, bwf_rad, current_prob_included;
    unsigned i, k, i_table, table_size;
    Roaming_Core_Input *the_pst_input;
    Roaming_Phasespace_Core *the_pst;
    Roaming_Result *the_PST_output;
    std::string file_name, jvE_file_name;
    std::vector<Roaming_Fraction> best_fits(r_general_input->fit_points.size());
    std::vector< std::vector<Roaming_BoltzmannFraction> > fit_result;
    std::vector<Roaming_BoltzmannFraction> temp_fit_point_vector;
    Roaming_BoltzmannFraction temp_PST_Boltzmann_fraction;
    the_pst_input = new Roaming_Core_Input();
    the_PST_output = new Roaming_Result();
    //these three can probably be in a function set_standard_PST_input_values() or something
    the_pst_input->job_path = job_path;
    the_pst_input->output_to_file = r_general_input->objects_write_to_file;
    the_pst_input->output_histograms = r_general_input->histogram_write_to_file;
    persistent_output << " - fit Trot and Tvib and P_roam:" << std::endl;
    i_table=0;
    table_size = (int)((r_general_input->rotational_temerature_range.max-r_general_input->rotational_temerature_range.min)/r_general_input->rotational_temerature_increment+1);
    table_size *= (int)((r_general_input->vibrational_temerature_range.max-r_general_input->vibrational_temerature_range.min)/r_general_input->vibrational_temerature_increment+1);
    boost::multi_array<double,2> out_table(boost::extents[table_size][13]);
    for (Trot=r_general_input->rotational_temerature_range.min; Trot<=r_general_input->rotational_temerature_range.max; Trot+=r_general_input->rotational_temerature_increment){
        r_general_input->rotational_temerature=Trot;
        for (Tvib=r_general_input->vibrational_temerature_range.min; Tvib<=r_general_input->vibrational_temerature_range.max; Tvib+=r_general_input->vibrational_temerature_increment){
            if(r_general_input->allow_Tvib_less_than_Trot || Tvib >= Trot) {
                r_general_input->vibrational_temerature=Tvib;
                jvE_file_name="jvE_histogram_Tr"+boost::lexical_cast<std::string>(Trot)+"_Tv"+boost::lexical_cast<std::string>(Tvib)+".csv";
                delete Erovib_histogram; //the Phasespace_Calculation created an object at Erovib_histogram, so it's ok to call delete in the first loop before the new on the next line ... the final delete will be handled by Phasespace_Calculation too
                Erovib_histogram=new jvE_histogram(jvE_file_name, job_path, r_general_input, reactant);
                current_prob_included=jvE_i_max=0;//calculate how far down the jvE_histogram to go to get to r_general_input->Boltzmann_probability_included, doing this first means we can define the progress bar better and also means the inner loop for each fit point is slightly more efficient
                while (jvE_i_max<Erovib_histogram->size() && (current_prob_included*100<=r_general_input->Boltzmann_probability_included || r_general_input->Boltzmann_probability_included==100)) {
                    current_prob_included+=Erovib_histogram->jvE_hist[jvE_i_max].jvE_prob;
                    jvE_i_max++;
                }
                clear_semipersistent_output();
                semipersistent_output << "Trot (" << r_general_input->rotational_temerature_range.min << "-" << r_general_input->rotational_temerature_range.max << ") = " << Trot << std::endl;
                semipersistent_output << "Tvib (" << r_general_input->vibrational_temerature_range.min << "-" << r_general_input->vibrational_temerature_range.max << ") = " << Tvib << std::endl;
                fit_result.clear();
                for (i=0;i<r_general_input->fit_points.size();i++){
                    semipersistent_output <<"  excitation_energy="<<r_general_input->fit_points[i].energy<<std::endl;
                    print_semipersistent_output();
                    the_progress_bar.initialise(100,1,jvE_i_max, CALCULATION_PROGRESS_BAR); // initialise the progress bar (this has no effect on the calculation)
                    fit_result.push_back(temp_fit_point_vector);
                    k=num_roam=num_rad=num_total=bwf_roam=bwf_rad=f_roam=f_rad=0;
                    while (k<jvE_i_max) {
                        the_progress_bar.go(k, CALCULATION_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
                        r_reactant->initial_J=Erovib_histogram->jvE_hist[k].j;
                        the_pst_input->E_total=r_general_input->fit_points[i].energy+Erovib_histogram->jvE_hist[k].energy;
                        the_pst_input->E_available=the_pst_input->E_total-r_reactant->E_dissociation+r_reactant->delta_E_roam;
                        the_pst_input->job_name="PST_"+boost::lexical_cast<std::string>(r_general_input->fit_points[i].energy);
                        the_PST_output->reset();
                        the_pst=new Roaming_Phasespace_Core(the_pst_input,r_reactant,r_frag1,r_frag2,the_PST_output);
                        delete the_pst;
                        if (the_PST_output->r_data->grand_total) {
                            if_roam=the_PST_output->r_data->num_may_roam/the_PST_output->r_data->grand_total;
                            if_rad=the_PST_output->r_data->num_radical/the_PST_output->r_data->grand_total;
                            temp_PST_Boltzmann_fraction.f_roam=if_roam;
                            temp_PST_Boltzmann_fraction.f_rad=if_rad;
                            temp_PST_Boltzmann_fraction.jvE_prob=Erovib_histogram->jvE_hist[k].jvE_prob;
                            fit_result[i].push_back(temp_PST_Boltzmann_fraction);
                        }
                        k++;
                    }
                    the_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
                }
                out_table[i_table][5]=Erovib_histogram->p_energy_within;
                out_table[i_table][6]=Erovib_histogram->p_prob_within;
                out_table[i_table][7]=Erovib_histogram->p_energy_outside;
                out_table[i_table][8]=Erovib_histogram->p_prob_outside;
                out_table[i_table][9]=Erovib_histogram->e_energy_within;
                out_table[i_table][10]=Erovib_histogram->e_prob_within;
                out_table[i_table][11]=Erovib_histogram->e_energy_outside;
                out_table[i_table][12]=Erovib_histogram->e_prob_outside;
                out_table[i_table][4]=Trot;
                out_table[i_table][0]=Tvib;
                selectBestProam(out_table, i_table, fit_result, best_fits);
                i_table++;
            }
        }
    }
    if (table_size<i_table) persistent_output << " SHOULD HAVE USED CEIL!" << std::endl;
    else if (table_size>i_table) table_size=i_table; // because if we exclude temp combos where Tvib < Trot we won't use the whole out table, giving us lots of 0, 0, 0 rows in the output!
    fprintf(out_file.fp, "Trot,Tvib,p_roam,rmse,rmse_ex_p_roam,p_energy_within,p_prob_within,p_energy_outside,p_prob_outside,e_energy_within,e_prob_within,e_energy_outside,e_prob_outside\n");
    for (i_table=0;i_table<table_size;i_table++){
        fprintf(out_file.fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", out_table[i_table][4], out_table[i_table][0], out_table[i_table][1], out_table[i_table][2], out_table[i_table][3], out_table[i_table][5], out_table[i_table][6], out_table[i_table][7], out_table[i_table][8], out_table[i_table][9], out_table[i_table][10], out_table[i_table][11], out_table[i_table][12]); //index 4 doesn't exist for other fit function ... 1,2,3 set in getBestProam
    }
    delete the_pst_input;
    return true;
}
//***********************************************************************************************************************************************************************************


//***********************************************************************************************************************************************************************************
class Roaming_Fit_Delta_E_Roam_And_P_Roam_Calculation : public Roaming_Calculation
{
    public:
        virtual ~Roaming_Fit_Delta_E_Roam_And_P_Roam_Calculation(){};
        virtual void set_energies();
        virtual bool run();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Fit_Delta_E_Roam_And_P_Roam_Calculation::set_energies(){
    //what is important for this calc is that E_max_statistical is set to the maximum possible value ... so vib levels are initialised properly ... so highest energy fit point and largest dEroam shoudl give us that
    unsigned i;
    double fit_point_E_max(0);
    for (i=0;i<r_general_input->fit_points.size();i++) if (fit_point_E_max<r_general_input->fit_points[i].energy) fit_point_E_max=r_general_input->fit_points[i].energy;
    E_total=fit_point_E_max;
    E_available_statistical=E_total-r_reactant->E_dissociation+r_general_input->delta_E_roam_range.max;
    E_max_statistical=E_available_statistical+parent_Erovib_max; //add the maximum parent internal energy so that adding that energy doesn' put the answers outside histogram limits
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Roaming_Fit_Delta_E_Roam_And_P_Roam_Calculation::run(){
    double delta_E_roam, num_roam, num_rad, num_total, f_roam, f_rad, if_roam, if_rad, bwf_roam, bwf_rad, current_prob_included;
    unsigned i, k, i_table, table_size;
    bool result=false;
    Roaming_Core_Input *the_pst_input;
    Roaming_Phasespace_Core *the_pst;
    Roaming_Result *the_PST_output;
    std::string file_name;
    std::vector<Roaming_Fraction> best_fits(r_general_input->fit_points.size());
    std::vector< std::vector<Roaming_BoltzmannFraction> > fit_result;
    std::vector<Roaming_BoltzmannFraction> temp_fit_point_vector;
    Roaming_BoltzmannFraction temp_PST_Boltzmann_fraction;
    Output_File best_fits_file;
    the_pst_input = new Roaming_Core_Input();
    the_PST_output = new Roaming_Result();
    //these three can probably be in a function set_standard_PST_input_values() or something
    the_pst_input->job_path = job_path;
    the_pst_input->output_to_file = r_general_input->objects_write_to_file;
    the_pst_input->output_histograms = r_general_input->histogram_write_to_file;
    best_fits_file.set_output_directory(job_path+"best_fits/");
    best_fits_file.open("best.csv");
    fprintf(best_fits_file.fp, "delta_E_roam,product");
    for (i=0;i<r_general_input->fit_points.size();i++) fprintf(best_fits_file.fp, ",%lf", r_general_input->fit_points[i].energy);
    fprintf(best_fits_file.fp, "\n");
    persistent_output << " - Fit delta_E_roam and P_roam:" << std::endl;
    i_table=0;

    table_size = (int)((r_general_input->delta_E_roam_range.max-r_general_input->delta_E_roam_range.min)/r_general_input->delta_E_roam_increment+1);
    boost::multi_array<double,2> out_table(boost::extents[table_size][4]); // [index][deroam, proam, mean squared error, proam ex mean squared error]
    for (delta_E_roam=r_general_input->delta_E_roam_range.min; delta_E_roam<=r_general_input->delta_E_roam_range.max; delta_E_roam+=r_general_input->delta_E_roam_increment){
        clear_semipersistent_output();
        semipersistent_output << "delta_E_roam (" << r_general_input->delta_E_roam_range.min << "-" << r_general_input->delta_E_roam_range.max << ") = " << delta_E_roam << std::endl;
        print_semipersistent_output();
        r_reactant->delta_E_roam = delta_E_roam; //the core will access this value to distinguish between roaming and radical states ... so while (at a glance) it might look like this is redundant it isn't
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
                the_pst_input->job_name="PST_"+boost::lexical_cast<std::string>(delta_E_roam)+"_"+boost::lexical_cast<std::string>(r_general_input->fit_points[i].energy);
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
        out_table[i_table][0]=delta_E_roam;
        selectBestProam(out_table, i_table, fit_result, best_fits);
        i_table++;
        fprintf(best_fits_file.fp, "%lf,roam", delta_E_roam);
        for (i=0;i<r_general_input->fit_points.size();i++) fprintf(best_fits_file.fp, ",%lf", best_fits[i].f_roam);
        fprintf(best_fits_file.fp, "\n");
        fprintf(best_fits_file.fp, "%lf,rad", delta_E_roam);
        for (i=0;i<r_general_input->fit_points.size();i++) fprintf(best_fits_file.fp, ",%lf", best_fits[i].f_rad);
        fprintf(best_fits_file.fp, "\n");
    }
    if (table_size<i_table) persistent_output << " SHOULD HAVE USED CEIL!" << std::endl;
    fprintf(out_file.fp, "delta_E_roam,p_roam,rmse,rmse_ex_p_roam\n");
    for (i_table=0;i_table<table_size;i_table++){
        fprintf(out_file.fp, "%lf,%lf,%lf,%lf\n", out_table[i_table][0], out_table[i_table][1], out_table[i_table][2], out_table[i_table][3]);
    }
    result=true;
    best_fits_file.close();
    delete the_pst_input;
    delete the_PST_output;
    return result;
}
//***********************************************************************************************************************************************************************************



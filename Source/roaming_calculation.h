//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Roaming_Fraction
{
        double f_roam, f_rad;
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Roaming_BoltzmannFraction
{
        double f_roam, f_rad, jvE_prob;
};
//***********************************************************************************************************************************************************************************
class Roaming_Fraction_Result
{
    public:
        Output_File out_file;
        bool opened_file;
        double num_roam, num_rad, num_total, f_roam, f_rad, if_roam, if_rad, bwf_roam, bwf_rad;
        Roaming_Fraction_Result();
        void close() {out_file.close();}
        ~Roaming_Fraction_Result(){if(opened_file) close();}
        void initialise(Output_File &commandeered_out_file); //this initialise will take over another output file that must have been opened (and must be closed) outside of this object ... I'm not sure I like this really, feels kinda hacky :-/
        void initialise(const std::string &job_path, const std::string &probe_name);
        void write_header();
        void reset() {num_roam=num_rad=num_total=bwf_roam=bwf_rad=f_roam=f_rad=0;}
        void i_reset() {if_roam=if_rad=0;}
        void i_increment(const double &num_may_roam, const double &num_radical, const double &grand_total, const double &P_roam, const double &jvE_prob);
        void write_result(const double &total_energy);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Fraction_Result::Roaming_Fraction_Result(){
    num_roam=num_rad=num_total=f_roam=f_rad=if_roam=if_rad=bwf_roam=bwf_rad=0;
    opened_file=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Fraction_Result::write_header(){
    fprintf(out_file.fp, "total_energy,Total Num States,Num Can Form Radicals, Num Can't Form Radicals, f_rad, f_roam, q_rad, q_roam\n");
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Fraction_Result::initialise(Output_File &commandeered_out_file){ //simply sets this object's Output_File file pointer to the address of the pointer of file being commandeered ... playing games with pointers is not awesome ... this may need a rethink
    out_file.fp = commandeered_out_file.fp; //copy pointer addresss!
    opened_file=0; //make sure this object doesn't try to close the file when it destructs because the commandeered_out_file will do that when it destructs and doing it twice would crash the program!
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Fraction_Result::initialise(const std::string &job_path, const std::string &probe_name){
    out_file.set_output_directory(job_path);
    out_file.open(probe_name+".csv");
    opened_file=1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Fraction_Result::i_increment(const double &num_may_roam, const double &num_radical, const double &grand_total, const double &P_roam, const double &jvE_prob){
    if (grand_total){
        if_roam=num_may_roam/grand_total;
        if_rad=num_radical/grand_total;
        if_roam*=P_roam; //apply the P_roam value -- master eqn stuff will have to happen here :-/
        normaliseFractions(if_roam, if_rad);
        bwf_roam+=if_roam*jvE_prob; ///must weight fraction, not state count! ... otherwise the marge state count values are like having more molecules in higher energy states!
        bwf_rad+=if_rad*jvE_prob;
        num_total+=grand_total*jvE_prob;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Fraction_Result::write_result(const double &total_energy){
    if (num_total) {
        f_roam=bwf_roam; ///bwf_roam & rad are quantum yield (ie they are the fraction of roaming and radical and molecules that do nothing, that is molecules that do nothing other than fluoresce etc) *** ONLY FOR OUR CURRENT SIMPLE P_ROAM SCALING FACTOR!!! RETHINK THIS WHEN I DO MASTER EQN!
        f_rad=bwf_rad;
        normaliseFractions(f_roam, f_rad); ///normalise quantum yields to get branching fractions
        num_roam=bwf_roam*num_total; //these numbers are meaningless! they are NOT state counts after Boltzmann weighting
        num_rad=bwf_rad*num_total;
    }else{
        num_roam=num_rad=f_roam=f_rad=-1;
    }
    fprintf(out_file.fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", total_energy, num_total, num_rad, num_roam, f_rad, f_roam, bwf_rad, bwf_roam);
}
//***********************************************************************************************************************************************************************************
class Roaming_Calculation : public Phasespace_Calculation
{
    public:
        Roaming_General_Input *r_general_input;
        Roaming_Molecule *r_reactant, *r_frag1, *r_frag2;

        Roaming_Calculation(){}; //C++ doesn't do virtual constructors, so we can't do much with this
        Roaming_Calculation(const std::string &file_name){go(file_name);};
        virtual ~Roaming_Calculation(){}; //no new objects ... they should be deleted by ~Phasespace_Calculation()
        virtual void create_general_input();
        virtual void create_general_input_pointer();
        virtual void create_fragments();
        virtual void create_fragment_pointers();
        virtual void set_energies();
        virtual void initialise(const std::string &file_name);
        void selectBestProam(boost::multi_array<double,2> &out_table, unsigned i_table, const std::vector< std::vector<Roaming_BoltzmannFraction> > &fit_result, std::vector<Roaming_Fraction> &best_fits);
        virtual bool run();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Calculation::create_general_input(){
    the_general_input = new Roaming_General_Input();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Calculation::create_general_input_pointer(){
    r_general_input = dynamic_cast<Roaming_General_Input*>(the_general_input);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Calculation::create_fragments(){
    reactant = new Roaming_Molecule();
    reactant->construct();
    frag1 = new Roaming_Molecule();
    frag1->construct();
    frag2 = new Roaming_Molecule();
    frag2->construct();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Calculation::create_fragment_pointers(){
    r_reactant = dynamic_cast<Roaming_Molecule*>(reactant);
    r_frag1 = dynamic_cast<Roaming_Molecule*>(frag1);
    r_frag2 = dynamic_cast<Roaming_Molecule*>(frag2);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Calculation::set_energies(){
    E_total=r_general_input->excitation_energy_range.max;
    E_available_statistical=E_total-r_reactant->E_dissociation+r_reactant->delta_E_roam; //note that this is E_avail to PST count ... E_avail at the potential corresponding to deltaEroam ... roaming PST core will divide these states in to roaming and dissoc ...
    E_max_statistical=E_available_statistical+parent_Erovib_max; //add the maximum parent internal energy so that adding that energy doesn' put the answers outside histogram limits
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Calculation::initialise(const std::string &file_name){
    Phasespace_Calculation::initialise(file_name);
    //set histogram defaults before parse sets any user overridden settings -- these override the default settings set in Phasespace_Calculation::initialise()
    reactant->is_on->E_trans=false;
    reactant->is_on->frag1_vel=false;
    reactant->is_on->frag2_vel=false;
    reactant->probe_is_on=reactant->is_on->clone(); //same defaults for probes
    frag1->is_on=reactant->is_on->clone(); //same defaults for secondary histograms
    frag1->probe_is_on=reactant->is_on->clone(); //same defaults for secondary histograms
    frag2->is_on=reactant->is_on->clone(); //same defaults for secondary histograms
    frag2->probe_is_on=reactant->is_on->clone(); //same defaults for secondary histograms
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Calculation::selectBestProam(boost::multi_array<double,2> &out_table, unsigned i_table, const std::vector< std::vector<Roaming_BoltzmannFraction> > &fit_result, std::vector<Roaming_Fraction> &best_fits) {
    unsigned i, j, num_results;
    bool first=1; //the first time we calculate the rmse and rmse_exProam we should automatically store it
    double jf_roam, jf_rad, xf_roam, xf_rad, P_roam, rmse, rmse_exProam, se, se_exProam;
    std::vector<double> f_rad(r_general_input->fit_points.size());
    std::vector<double> f_roam(r_general_input->fit_points.size());
    for (P_roam=r_general_input->P_roam_range.min; P_roam<=r_general_input->P_roam_range.max; P_roam+=r_general_input->P_roam_increment){
        num_results=0;
        se=se_exProam=0;
        for (i=0;i<r_general_input->fit_points.size();i++){
            xf_roam=xf_rad=f_roam[i]=f_rad[i]=0;
            for (j=0;j<fit_result[i].size();j++) {
                if (first) {
                    xf_roam+=fit_result[i][j].f_roam*fit_result[i][j].jvE_prob; //increment ex proam fraction for this boltzmann weighted energy
                    xf_rad+=fit_result[i][j].f_rad*fit_result[i][j].jvE_prob;
                }
                jf_roam=fit_result[i][j].f_roam*P_roam; //apply the P_roam value we are testing
                jf_rad=fit_result[i][j].f_rad;
                normaliseFractions(jf_roam, jf_rad);
                f_roam[i]+=jf_roam*fit_result[i][j].jvE_prob;//increment with proam fraction for this boltzmann weighted energy
                f_rad[i]+=jf_rad*fit_result[i][j].jvE_prob;
            }
            if (!r_general_input->fit_points.is_quantum) normaliseFractions(f_roam[i], f_rad[i]); //normalise to convert quantum yield to branching fraction
            if (r_general_input->fit_points[i].is_roaming) {
                for (j=0; j<r_general_input->fit_points[i].roaming.size(); j++) {
                    if (first) se_exProam+=pow((r_general_input->fit_points[i].roaming[j]-xf_roam),2);
                    se+=pow((r_general_input->fit_points[i].roaming[j]-f_roam[i]),2);
                    num_results++;
                }
            }
            if (r_general_input->fit_points[i].is_radical) {
                for (j=0; j<r_general_input->fit_points[i].radical.size(); j++) {
                    if (first) se_exProam+=pow((r_general_input->fit_points[i].radical[j]-xf_rad),2);
                    se+=pow((r_general_input->fit_points[i].radical[j]-f_rad[i]),2);
                    num_results++;
                }
            }
        }
        rmse=sqrt(se/num_results);
        if (out_table[i_table][2]>rmse || first) {
            out_table[i_table][1]=P_roam;
            out_table[i_table][2]=rmse;
            for (i=0;i<r_general_input->fit_points.size();i++){
                best_fits[i].f_rad=f_rad[i];
                best_fits[i].f_roam=f_roam[i];
            }
        }
        if (first) {
            rmse_exProam=sqrt(se_exProam/num_results);
            out_table[i_table][3]=rmse_exProam;
            first=0;
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Roaming_Calculation::run(){ //this function is not complete ... see block comment inside loop for details ...
    unsigned i, p;

    Roaming_Core_Input *the_pst_input;
    Roaming_Phasespace_Core *the_pst;
    Roaming_Result *the_PST_output; //stores the output from a single roaming PST count
    Roaming_Fraction_Result the_roaming_fraction; //stores the overall fractions, for a single photolysis energy, for however may PST counts are required to acquire them (which could be many with parent internal energy considered)
    std::vector<Roaming_Fraction_Result> probe_results(r_reactant->probe_states.size()); //eventually the overall result will be done with one of these too


    the_pst_input = new Roaming_Core_Input();
    the_PST_output = new Roaming_Result();
    //these three can probably be in a function set_standard_PST_input_values() or something
    the_pst_input->job_path = job_path;
    the_pst_input->output_to_file = r_general_input->objects_write_to_file;
    the_pst_input->output_histograms = r_general_input->histogram_write_to_file;
    //these will be manipulated differently by different calculations
    the_pst_input->job_name="Roaming_Calculation";
    the_PST_output->initialise(limits, r_reactant);
    the_roaming_fraction.initialise(out_file); //initialise the overall Roaming_Fraction_Result so that it writes to the main calculation output file.
    for(p=0; p<r_reactant->probe_states.size(); p++) probe_results[p].initialise(job_path, r_reactant->probe_states[p].name);
    persistent_output << " - Roaming Calculation:" << std::endl;
    print_persistent_output();
    the_roaming_fraction.write_header();
    for(p=0; p<probe_results.size(); p++) probe_results[p].write_header();
    i=0;
    while (i<jvE_i_max) {
        the_progress_bar.go(i, CALCULATION_PROGRESS_BAR);//update progress bar (this has no effect on the calculation)
        r_reactant->initial_J=Erovib_histogram->jvE_hist[i].j;
        the_pst_input->E_total=E_total+Erovib_histogram->jvE_hist[i].energy;
        the_pst_input->E_available=the_pst_input->E_total-r_reactant->E_dissociation+r_reactant->delta_E_roam;
        the_pst_input->job_name="PST_"+boost::lexical_cast<std::string>(r_reactant->delta_E_roam)+"_"+boost::lexical_cast<std::string>(E_total);
        the_PST_output->reset();
        the_pst=new Roaming_Phasespace_Core(the_pst_input,r_reactant,r_frag1,r_frag2,the_PST_output);
        delete the_pst;
        the_roaming_fraction.i_reset();
        for(p=0; p<probe_results.size(); p++) probe_results[p].i_reset();
        if (the_PST_output->r_data->grand_total) {

            /*
            for histograms to work with the jvE stuff, we'll need to increment the histograms etc for each jvE result ... this will require some more thought!
            We do something like this with the Secondary_Dissociation_Result::increment() function, except that function calculates the actual scaling factor inside

            Once we write a scaled_increment function for Phasespace_Result and it's derivatives, Secondary_Dissociation_Result::increment() should be renamed to 3F_increment and
            should calculate the scaling factor and then call the general scaled_increment passing the scaling factor (keeping the divide by large number rather than multiply by small number makes sence)
            Once we have the function, we can have an overall result object that gets incremented for each jvE energy ...

            We also need a way to choose between MPI for each PST count or MPI for the list of jvE energyies ... this could be automatic (ie: if jvE_i_max > world.size() then MPI jvE else MPI PST counts)
            or it could be a user input (perhaps even both, if not specified then decide automatically as above otherwise allow the user to force one option or the other).
            */

            the_roaming_fraction.i_increment(the_PST_output->r_data->num_may_roam, the_PST_output->r_data->num_radical, the_PST_output->r_data->grand_total, r_reactant->P_roam, Erovib_histogram->jvE_hist[i].jvE_prob);
            for(p=0; p<probe_results.size(); p++) probe_results[p].i_increment(the_PST_output->r_probes[p]->r_data->num_may_roam, the_PST_output->r_probes[p]->r_data->num_radical, the_PST_output->r_probes[p]->r_data->grand_total, r_reactant->P_roam, Erovib_histogram->jvE_hist[i].jvE_prob);
        }
        i++;
    }
    the_progress_bar.un_initialise(CALCULATION_PROGRESS_BAR);//release the progress bar
    the_roaming_fraction.write_result(E_total);
    for(p=0; p<probe_results.size(); p++) probe_results[p].write_result(E_total);


    delete the_pst_input;
    delete the_PST_output;
    return true;
}
//***********************************************************************************************************************************************************************************



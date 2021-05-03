//***********************************************************************************************************************************************************************************
class Impulsive_Molecule : public Phasespace_Molecule
{
    public:
        double E_barrier; //cm-1 //if this is non zero, an impulsive treatment will occur
        bool use_angle_distribution, write_TS_inertia_xyz;
        unsigned bend_mode; //only used to identify how many quanta of the bend mode exist in the parent state (not sure this is nessesary as that number of quanta is before getting to the TS config, which will use energy so may not still have vibs still)
        double TS_reaction_coord_force_constant, TS_reaction_coord_reduced_mass; //for Eckart tunneling (these apply to the parent, where the TS coordinates are defined)
        double TS_bend_frequency, TS_bend_reduced_mass, TS_bend_arc_radius; //for HO wavefunction (FOR NOW ... these apply to the parent, where the TS coordinates are defined) I THINK WE NEED TO RE THINK THIS FOR BIGGER MOLECULES!!
        double TS_frag1_equilibrium_angle, TS_frag2_equilibrium_angle;//the TS angle between the breaking bond and the centre of mass of this fragment (this value applies to children, calculated)
        double TS_frag1_pivot_to_com_length, TS_frag2_pivot_to_com_length;//the TS distance from the pivot atom to the centre of mass, for each fragment (TS defined once for parent)
        double TS_frag1_moment_of_inertia, TS_frag2_moment_of_inertia;//the moment of inertia in the TS configuration, for each fragment (TS defined once for parent)
        double TS_frag1_kn_ratio, TS_frag2_kn_ratio;//the value of cos(theta) where theta is the angle between the rotation and principal axes of rotation of the fragment due to the impulse
        double initial_velocity;

        Impulsive_Molecule(){};
        virtual void construct();
        virtual Impulsive_Molecule* clone(){return new Impulsive_Molecule(*this);};
        virtual ~Impulsive_Molecule(){};//persistent_output << endl << "Killing an Impulsive_Molecule!" << endl;};
        virtual void parse_input(const std::string &file_name, const std::string &section_name);
        void get_impulse_values(FILE *the_file, const unsigned &num_axis_atoms, Geom_Molecule *TS_geom, const std::vector<unsigned> &TS_frag_atoms, const double &TS_this_frag_pivot_atom, const double &TS_other_frag_pivot_atom, double &TS_equilibrium_angle, double &TS_pivot_to_com_length, double &TS_moment_of_inertia, double &TS_kn_ratio);
        void get_impulse_nk(const double &E_iRot, const double &kn_ratio, unsigned &i_n, unsigned &i_k, double &i_E_nk) const; //works out which
    private:
        Impulsive_Molecule(const Impulsive_Molecule &rhs);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Phasespace_Molecule>(*this);
            ar & E_barrier & use_angle_distribution & write_TS_inertia_xyz & bend_mode & TS_reaction_coord_force_constant & TS_reaction_coord_reduced_mass & TS_bend_frequency & TS_bend_reduced_mass & TS_bend_arc_radius;
            ar & TS_frag1_equilibrium_angle & TS_frag2_equilibrium_angle & TS_frag1_pivot_to_com_length & TS_frag2_pivot_to_com_length & TS_frag1_moment_of_inertia & TS_frag2_moment_of_inertia & TS_frag1_kn_ratio & TS_frag2_kn_ratio & initial_velocity;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Molecule::construct(){
    Phasespace_Molecule::construct();
    E_barrier=0;
    use_angle_distribution=write_TS_inertia_xyz=0;
    bend_mode=0;
    TS_reaction_coord_force_constant=TS_reaction_coord_reduced_mass=0;
    TS_bend_frequency=TS_bend_reduced_mass=TS_bend_arc_radius=0;
    TS_frag1_equilibrium_angle=TS_frag2_equilibrium_angle=TS_frag1_pivot_to_com_length=TS_frag2_pivot_to_com_length=TS_frag1_moment_of_inertia=TS_frag2_moment_of_inertia=TS_frag1_kn_ratio=TS_frag2_kn_ratio=0;
    initial_velocity=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Impulsive_Molecule::Impulsive_Molecule(const Impulsive_Molecule &rhs) : Phasespace_Molecule(rhs){
    E_barrier=rhs.E_barrier;
    use_angle_distribution=rhs.use_angle_distribution;
    write_TS_inertia_xyz=rhs.write_TS_inertia_xyz;
    bend_mode=rhs.bend_mode;
    TS_bend_frequency=rhs.TS_bend_frequency;
    TS_bend_reduced_mass=rhs.TS_bend_reduced_mass;
    TS_bend_arc_radius=rhs.TS_bend_arc_radius;
    TS_frag1_equilibrium_angle=rhs.TS_frag1_equilibrium_angle;
    TS_frag2_equilibrium_angle=rhs.TS_frag2_equilibrium_angle;
    TS_frag1_pivot_to_com_length=rhs.TS_frag1_pivot_to_com_length;
    TS_frag2_pivot_to_com_length=rhs.TS_frag2_pivot_to_com_length;
    TS_frag1_moment_of_inertia=rhs.TS_frag1_moment_of_inertia;
    TS_frag2_moment_of_inertia=rhs.TS_frag2_moment_of_inertia;
    TS_frag1_kn_ratio=rhs.TS_frag1_kn_ratio;
    TS_frag2_kn_ratio=rhs.TS_frag2_kn_ratio;
    TS_reaction_coord_force_constant=rhs.TS_reaction_coord_force_constant;
    TS_reaction_coord_reduced_mass=rhs.TS_reaction_coord_reduced_mass;
    initial_velocity=rhs.initial_velocity;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Molecule::parse_input(const std::string &file_name, const std::string &section_name){
    Phasespace_Molecule::parse_input(file_name, section_name); //first parse the parent class's members
    const unsigned NO_SECTION=0;
    const unsigned GENERAL=1;
    const unsigned TS_ATOMS=2;
	unsigned read_state(NO_SECTION);
	std::string file_name_and_path(input_directory+file_name);
    std::ifstream infile(file_name_and_path.c_str());
	std::string line, key, value, symbol;
	double ax, ay, az;
	std::vector<unsigned> TS_frag1_atoms, TS_frag2_atoms;//the indicies for this fragment, as defined in the TS
    unsigned TS_frag1_pivot_atom(0), TS_frag2_pivot_atom(0); //the index of the atom of this fragment in the breaking bond, as defined in the TS atoms input
    Geom_Molecule *TS_geom; //will be created and used to derive some data from imput coordinates and then deleted
    TS_geom=new Geom_Molecule(section_name); //if there is a barrier, there is a TS, so we need a Geom_Molecule
    Output_File xyz_file;
    std::string job_name;
    unsigned i, xyz_count, num_axis_atoms(6);
    semipersistent_output << "Impulsive_Molecule "<< number <<" parsing" << std::endl;
    print_semipersistent_output();
	if (infile.is_open()) {
        while(getline(infile,line)) {
            if (line.substr(0,2) != "//"){ //if this is not a comment line (we should ignore comment lines)
                if (line.substr(0,8) == "--------"){ //this is the start of a new section
                    read_state=NO_SECTION; //set to none first in case the section is not recognised by this parse routine
                    if (line.find( std::string("-")+section_name+std::string("-") ) != std::string::npos) read_state=GENERAL;
                    else if (line.find( std::string("-")+section_name+std::string(" TS ATOMS-") ) != std::string::npos) read_state=TS_ATOMS;
                }
                if (read_state==GENERAL){
                    read_key_value(line, key, value);
                    if (key != "") {
                        if (key == "barrier_energy"){
                            E_barrier=atof(value.c_str());
                        }
                        else if (key == "bend_mode") bend_mode=atoi(value.c_str());
                        else if (key == "use_angle_distribution") use_angle_distribution=atoi(value.c_str());
                        else if (key == "write_TS_inertia_xyz") write_TS_inertia_xyz=atoi(value.c_str());
                        else if (key == "TS_bend_frequency") TS_bend_frequency=atof(value.c_str()); //do we need more than one now???
                        else if (key == "TS_bend_reduced_mass") TS_bend_reduced_mass=atof(value.c_str()); //do we need more than one now???
                        else if (key == "TS_bend_arc_radius") TS_bend_arc_radius=atof(value.c_str()); //do we need more than one now???
                        else if (key == "TS_reaction_coord_force_constant") TS_reaction_coord_force_constant=atof(value.c_str());
                        else if (key == "TS_reaction_coord_reduced_mass") TS_reaction_coord_reduced_mass=atof(value.c_str());
                        else if (key == "TS_frag1_pivot_atom") TS_frag1_pivot_atom=atoi(value.c_str());
                        else if (key == "TS_frag2_pivot_atom") TS_frag2_pivot_atom=atoi(value.c_str());
                        else if (key == "TS_frag1_atoms") read_index_list(value, TS_frag1_atoms);
                        else if (key == "TS_frag2_atoms") read_index_list(value, TS_frag2_atoms);
                    }
                }else if (read_state==TS_ATOMS){ //if a barrier energy is defined, then TS Atoms must be defined in the file! The barrier energy must be defined before the atoms in the file!
                    if (E_barrier && read_atom_value(line, symbol, ax, ay, az)) TS_geom->add_atom(symbol, ax, ay, az);
                }
            }
        }
        infile.close();
        if (E_barrier){ //if there is a barrier defined then calculate the values required for the impulse calculation
            if (TS_frag1_atoms.size() && TS_frag2_atoms.size() && TS_geom->atoms.size()){ //check that the atoms in each fragment have been assigned (there are lots more error checking things that we could/should eventually do here)
                if (write_TS_inertia_xyz){
                    job_name=file_name;
                    job_name.erase(job_name.find_last_of("."), std::string::npos);
                    xyz_file.set_output_directory(output_directory+job_name+"/xyz/");
                    xyz_file.open(section_name+"_TS_inertia_axes.xyz");
                    xyz_count=TS_geom->atoms.size();
                    if (TS_frag1_atoms.size()>1) xyz_count+=1+(4*num_axis_atoms);//for each fragment, +1 for COM and +(4*num_axis_atoms) 2 points either side for each axis, principal and rotation, repeated num_axis_atoms times along axis
                    if (TS_frag2_atoms.size()>1) xyz_count+=1+(4*num_axis_atoms);
                    fprintf(xyz_file.fp,"%d\ncomment: TS geometry with centre of mass and principal and rotation axes illustrated with additional atoms\n",xyz_count);
                    for (i=0; i<TS_geom->atoms.size(); i++) fprintf(xyz_file.fp,"%s   %e   %e   %e\n", TS_geom->atoms[i].symbol.c_str(), TS_geom->atoms[i].coord[0], TS_geom->atoms[i].coord[1], TS_geom->atoms[i].coord[2]);
                }
                if (TS_frag1_atoms.size()>1) get_impulse_values(xyz_file.fp, num_axis_atoms, TS_geom, TS_frag1_atoms, TS_frag1_pivot_atom, TS_frag2_pivot_atom, TS_frag1_equilibrium_angle, TS_frag1_pivot_to_com_length, TS_frag1_moment_of_inertia, TS_frag1_kn_ratio);
                if (TS_frag2_atoms.size()>1) get_impulse_values(xyz_file.fp, num_axis_atoms, TS_geom, TS_frag2_atoms, TS_frag2_pivot_atom, TS_frag1_pivot_atom, TS_frag2_equilibrium_angle, TS_frag2_pivot_to_com_length, TS_frag2_moment_of_inertia, TS_frag2_kn_ratio);
                if (write_TS_inertia_xyz) xyz_file.close();
            }else{
                if (!TS_frag1_atoms.size()) persistent_output<<"ERROR: "<<section_name<<" TS_frag1_atoms have not been defined in input!"<<std::endl;
                if (!TS_frag2_atoms.size()) persistent_output<<"ERROR: "<<section_name<<" TS_frag2_atoms have not been defined in input!"<<std::endl;
                if (!TS_geom->atoms.size()) persistent_output<<"ERROR: --------"<<section_name<<" TS ATOMS -------- have not been defined in input!"<<std::endl;
                print_persistent_output();
            }
        }
	}
	delete TS_geom;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Molecule::get_impulse_values(FILE *the_file, const unsigned &num_axis_atoms, Geom_Molecule *TS_geom, const std::vector<unsigned> &TS_frag_atoms,
                                            const double &TS_this_frag_pivot_atom, const double &TS_other_frag_pivot_atom, double &TS_equilibrium_angle, double &TS_pivot_to_com_length,
                                            double &TS_moment_of_inertia, double &TS_kn_ratio){
    double rot_princ_angle;
    std::vector<Geom_Eigen_Data> principal_inertia_data;
	std::vector<double> frag_com, bond_vector, com_vector, rot_vector, rot_hat, princ_hat;
	unsigned i;
    frag_com=TS_geom->get_sub_centre_of_mass_coordinates(TS_frag_atoms);
    TS_equilibrium_angle=TS_geom->get_angle(TS_other_frag_pivot_atom, TS_this_frag_pivot_atom, frag_com);
    bond_vector=TS_geom->get_vector(TS_this_frag_pivot_atom, TS_other_frag_pivot_atom);
    com_vector=TS_geom->get_vector(TS_this_frag_pivot_atom, frag_com);
    TS_pivot_to_com_length=vector_magnitude(com_vector);
    rot_vector=cross_product(bond_vector,com_vector,3);
    TS_moment_of_inertia=TS_geom->get_sub_inertia_axis_vector(TS_frag_atoms, frag_com, rot_vector);
    principal_inertia_data=TS_geom->get_sub_inertia_tensor(TS_frag_atoms, frag_com);
    rot_princ_angle=vector_angle(rot_vector, principal_inertia_data[2].vect); //get the angle between the rotation axis and the principal axis
    TS_kn_ratio=cos(rot_princ_angle);
    if (write_TS_inertia_xyz){
        rot_hat=vector_normalise(rot_vector);
        princ_hat=vector_normalise(principal_inertia_data[2].vect);
        fprintf(the_file,"Si   %e   %e   %e\n", frag_com[0], frag_com[1], frag_com[2]);
        for (i=0;i<num_axis_atoms;i++){
            fprintf(the_file,"F   %e   %e   %e\n", frag_com[0]+(princ_hat[0]*i), frag_com[1]+(princ_hat[1]*i), frag_com[2]+(princ_hat[2]*i));
            fprintf(the_file,"F   %e   %e   %e\n", frag_com[0]+(-1*princ_hat[0]*i), frag_com[1]+(-1*princ_hat[1]*i), frag_com[2]+(-1*princ_hat[2]*i));
            fprintf(the_file,"Cl   %e   %e   %e\n", frag_com[0]+(rot_hat[0]*i), frag_com[1]+(rot_hat[1]*i), frag_com[2]+(rot_hat[2]*i));
            fprintf(the_file,"Cl   %e   %e   %e\n", frag_com[0]+(-1*rot_hat[0]*i), frag_com[1]+(-1*rot_hat[1]*i), frag_com[2]+(-1*rot_hat[2]*i));
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Impulsive_Molecule::get_impulse_nk(const double &E_iRot, const double &kn_ratio, unsigned &i_n, unsigned &i_k, double &i_E_nk) const {/// This function sets the 3 parameters i_n, i_k & i_E_nk to the n and k values that best fit E_iRot and kn_ratio (min RMSE)
    unsigned n_max, to_k, n, k;
    double E_nk, E_percent_error, nk_percent_error, this_total_percent_error, total_percent_error;
    total_percent_error=E_nk=i_E_nk=i_n=i_k=0; //start with 0 and if there is a better fit these will be reset in the loop
    if (E_iRot>0) total_percent_error=200;//can't calculate deviation of kn_ratio from k/n as n=0 so the RMSE is just the energy difference for n=k=0, which will be 100% for non zero E_iRot ... we need 200 because two percent errors are added
    n_max=get_n_max(E_iRot);
    for (n=1;n<=n_max;n++){// from 1 because we don't want to divide by 0!
        if (linear) to_k=0;
        else to_k=n;
        for (k=0;k<=to_k;k++){
            E_nk=get_E_nk(n, k);
            //persistent_output<<endl<<"n="<<n<<" k="<<k<<" ";
            if (E_nk<=E_iRot){//this state is within the available energy (this check is nesesary because get_E_nk() uses ceil to ensure nothing is omitted)
                E_percent_error=100*(E_iRot-E_nk)/E_iRot;
                nk_percent_error=100*abs(kn_ratio-((double)k/(double)n))/kn_ratio;
                this_total_percent_error=E_percent_error+nk_percent_error;
                //persistent_output<<E_nk<<"<="<<E_iRot<<" this_total_percent_error="<<this_total_percent_error<<" E_percent_error="<<E_percent_error<<" nk_percent_error="<<nk_percent_error<<" ";
                if (this_total_percent_error<=total_percent_error){//smallest RMSE fit wins!
                    total_percent_error=this_total_percent_error;
                    i_n=n;
                    i_k=k;
                    i_E_nk=E_nk;
                    //persistent_output<<"current_min_error ";
                }
            }
        }
    }
    //persistent_output<<endl<<"chosen_n="<<i_n<<" chosen_k="<<i_k<<" E="<<i_E_nk<<endl<<endl;
}
//***********************************************************************************************************************************************************************************


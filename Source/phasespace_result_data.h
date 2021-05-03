
//***********************************************************************************************************************************************************************************
class Histogram_Limits
{
    public:
        double E_max, vel_max1, vel_max2, E_max1, E_max2;
        unsigned J_max1, J_max2, v_max1, v_max2;
        std::vector<std::string> frag1_vib_names, frag2_vib_names;

        Histogram_Limits();
        virtual Histogram_Limits* clone(){return new Histogram_Limits(*this);};
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & E_max & vel_max1 & vel_max2 & E_max1 & E_max2 & J_max1 & J_max2 & v_max1 & v_max2 & frag1_vib_names & frag2_vib_names;}
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Histogram_Limits::Histogram_Limits(){
    E_max=vel_max1=vel_max2=E_max1=E_max2=0;
    J_max1=J_max2=v_max1=v_max2=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
class Histogram_Use
{
    public:
        bool E_trans, frag1_vel, frag2_vel, frag1_E_trans, frag2_E_trans, frag1_E_int, frag2_E_int, frag1_E_vib, frag2_E_vib, frag1_E_rot, frag2_E_rot, frag1_v, frag2_v;
        bool frag1_nk, frag2_nk, frag1_E_int_trans, frag2_E_int_trans;

        Histogram_Use();
        virtual ~Histogram_Use() {}
        virtual Histogram_Use* clone(){return new Histogram_Use(*this);};
        virtual void read_histogram_use(const std::string &key, const std::string &value);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){
            ar & E_trans & frag1_vel & frag2_vel & frag1_E_trans & frag2_E_trans & frag1_E_int & frag2_E_int & frag1_E_vib & frag2_E_vib & frag1_E_rot & frag2_E_rot & frag1_v & frag2_v;
            ar & frag1_nk & frag2_nk & frag1_E_int_trans & frag2_E_int_trans;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Histogram_Use::Histogram_Use(){
    E_trans=frag1_vel=frag2_vel=frag1_E_trans=frag2_E_trans=frag1_E_int=frag2_E_int=frag1_E_vib=frag2_E_vib=frag1_E_rot=frag2_E_rot=frag1_v=frag2_v=0;
    frag1_nk=frag2_nk=frag1_E_int_trans=frag2_E_int_trans=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_Use::read_histogram_use(const std::string &key, const std::string &value) {
    if (key == "E_trans") E_trans=atoi(value.c_str());
    else if (key == "frag1_vel") frag1_vel=atoi(value.c_str());
    else if (key == "frag2_vel") frag2_vel=atoi(value.c_str());
    else if (key == "frag1_E_trans") frag1_E_trans=atoi(value.c_str());
    else if (key == "frag2_E_trans") frag2_E_trans=atoi(value.c_str());
    else if (key == "frag1_E_int") frag1_E_int=atoi(value.c_str());
    else if (key == "frag2_E_int") frag2_E_int=atoi(value.c_str());
    else if (key == "frag1_E_vib") frag1_E_vib=atoi(value.c_str());
    else if (key == "frag2_E_vib") frag2_E_vib=atoi(value.c_str());
    else if (key == "frag1_E_rot") frag1_E_rot=atoi(value.c_str());
    else if (key == "frag2_E_rot") frag2_E_rot=atoi(value.c_str());
    else if (key == "frag1_v") frag1_v=atoi(value.c_str());
    else if (key == "frag2_v") frag2_v=atoi(value.c_str());
    else if (key == "frag1_nk") frag1_nk=atoi(value.c_str());
    else if (key == "frag2_nk") frag2_nk=atoi(value.c_str());
    else if (key == "frag1_E_int_trans") frag1_E_int_trans=atoi(value.c_str());
    else if (key == "frag2_E_int_trans") frag2_E_int_trans=atoi(value.c_str());
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
class Phasespace_Result_Data
{
    public:
        Histogram E_trans, frag1_vel, frag2_vel, frag1_E_trans, frag2_E_trans, frag1_E_int, frag2_E_int, frag1_E_vib, frag2_E_vib, frag1_E_rot, frag2_E_rot, frag1_v, frag2_v;
        Histogram_2d frag1_nk, frag2_nk, frag1_E_int_trans, frag2_E_int_trans;
        double grand_total;

        virtual ~Phasespace_Result_Data() {}
        Phasespace_Result_Data& operator=(const Phasespace_Result_Data &rhs);
        Phasespace_Result_Data& operator+=(const Phasespace_Result_Data &rhs);
        virtual Phasespace_Result_Data* clone(){return new Phasespace_Result_Data(*this);};
        virtual void initialise(const Histogram_Limits &limits, const Histogram_Use *is_on);
        virtual void reset();
        virtual void write_totals(FILE *the_file);
        virtual void write_histograms(FILE *the_file);
        virtual void write_bitmaps(const std::string job_path_and_name_start);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){
            ar & E_trans & frag1_vel & frag2_vel & frag1_E_trans & frag2_E_trans & frag1_E_int & frag2_E_int & frag1_E_vib & frag2_E_vib;
            ar & frag1_E_rot & frag2_E_rot & frag1_v & frag2_v & frag1_nk & frag2_nk & frag1_E_int_trans & frag2_E_int_trans & grand_total;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Result_Data& Phasespace_Result_Data::operator=(const Phasespace_Result_Data &rhs){
    grand_total=rhs.grand_total;
    E_trans=rhs.E_trans;
    frag1_vel=rhs.frag1_vel;
    frag2_vel=rhs.frag2_vel;
    frag1_E_trans=rhs.frag1_E_trans;
    frag2_E_trans=rhs.frag2_E_trans;
    frag1_E_int=rhs.frag1_E_int;
    frag2_E_int=rhs.frag2_E_int;
    frag1_E_vib=rhs.frag1_E_vib;
    frag2_E_vib=rhs.frag2_E_vib;
    frag1_E_rot=rhs.frag1_E_rot;
    frag2_E_rot=rhs.frag2_E_rot;
    frag1_v=rhs.frag1_v;
    frag2_v=rhs.frag2_v;
    frag1_nk=rhs.frag1_nk;
    frag2_nk=rhs.frag2_nk;
    frag1_E_int_trans=rhs.frag1_E_int_trans;
    frag2_E_int_trans=rhs.frag2_E_int_trans;
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Result_Data& Phasespace_Result_Data::operator+=(const Phasespace_Result_Data &rhs){
    grand_total+=rhs.grand_total;
    E_trans+=rhs.E_trans;
    frag1_vel+=rhs.frag1_vel;
    frag2_vel+=rhs.frag2_vel;
    frag1_E_trans+=rhs.frag1_E_trans;
    frag2_E_trans+=rhs.frag2_E_trans;
    frag1_E_int+=rhs.frag1_E_int;
    frag2_E_int+=rhs.frag2_E_int;
    frag1_E_vib+=rhs.frag1_E_vib;
    frag2_E_vib+=rhs.frag2_E_vib;
    frag1_E_rot+=rhs.frag1_E_rot;
    frag2_E_rot+=rhs.frag2_E_rot;
    frag1_v+=rhs.frag1_v;
    frag2_v+=rhs.frag2_v;
    frag1_nk+=rhs.frag1_nk;
    frag2_nk+=rhs.frag2_nk;
    frag1_E_int_trans+=rhs.frag1_E_int_trans;
    frag2_E_int_trans+=rhs.frag2_E_int_trans;
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result_Data::initialise(const Histogram_Limits &limits, const Histogram_Use *is_on) {
    //probably should do something about the names!
    grand_total=0;
    if (is_on->E_trans) E_trans.initialiseN("Relative Translational Energies", "Translational Energy (cm-1)", "Count", HISTOGRAM_SIZE, limits.E_max);
    if (is_on->frag1_vel) frag1_vel.initialiseN("Velocities For Fragment 1", "Velocity (m/s)", "Count", HISTOGRAM_SIZE, limits.vel_max1);
    if (is_on->frag2_vel) frag2_vel.initialiseN("Velocities For Fragment 2", "Velocity (m/s)", "Count", HISTOGRAM_SIZE, limits.vel_max2);
    if (is_on->frag1_E_trans) frag1_E_trans.initialiseN("Translational Energies For Fragment 1", "Translational Energy (cm-1)", "Count", HISTOGRAM_SIZE, limits.E_max1);
    if (is_on->frag2_E_trans) frag2_E_trans.initialiseN("Translational Energies For Fragment 2", "Translational Energy (cm-1)", "Count", HISTOGRAM_SIZE, limits.E_max2);
    if (is_on->frag1_E_int) frag1_E_int.initialiseN("Internal Energies For Fragment 1", "Internal Energy (cm-1)", "Count", HISTOGRAM_SIZE, limits.E_max);
    if (is_on->frag2_E_int) frag2_E_int.initialiseN("Internal Energies For Fragment 2", "Internal Energy (cm-1)", "Count", HISTOGRAM_SIZE, limits.E_max);
    if (is_on->frag1_E_vib) frag1_E_vib.initialiseN("Vibrational Energies For Fragment 1", "Vibrational Energy (cm-1)", "Count", HISTOGRAM_SIZE, limits.E_max);
    if (is_on->frag2_E_vib) frag2_E_vib.initialiseN("Vibrational Energies For Fragment 2", "Vibrational Energy (cm-1)", "Count", HISTOGRAM_SIZE, limits.E_max);
    if (is_on->frag1_E_rot) frag1_E_rot.initialiseN("Rotational Energies For Fragment 1", "Rotational Energy (cm-1)", "Count", HISTOGRAM_SIZE, limits.E_max);
    if (is_on->frag2_E_rot) frag2_E_rot.initialiseN("Rotational Energies For Fragment 2", "Rotational Energy (cm-1)", "Count", HISTOGRAM_SIZE, limits.E_max);
    if (is_on->frag1_v) frag1_v.initialiseS("Vibrational Levels For Fragment 1", "Index", "Count", 1, limits.v_max1);
    if (is_on->frag2_v) frag2_v.initialiseS("Vibrational Levels For Fragment 2", "Index", "Count", 1, limits.v_max2);
    if (is_on->frag1_v) frag1_v.set_names("Combination Band", limits.frag1_vib_names);
    if (is_on->frag2_v) frag2_v.set_names("Combination Band", limits.frag2_vib_names);
    if (is_on->frag1_nk) frag1_nk.initialiseS("J And K For Fragment 1", "Count", "J", "k", 1, 1, limits.J_max1, limits.J_max1);
    if (is_on->frag2_nk) frag2_nk.initialiseS("J And K For Fragment 2", "Count", "J", "k", 1, 1, limits.J_max2, limits.J_max2);
    if (is_on->frag1_E_int_trans) frag1_E_int_trans.initialiseN("E_int And E_trans For Fragment 1", "Count", "Internal Energy (cm-1)", "Translational Energy (cm-1)", HISTOGRAM_2D_SIZE, HISTOGRAM_2D_SIZE, limits.E_max, limits.E_max1);
    if (is_on->frag2_E_int_trans) frag2_E_int_trans.initialiseN("E_int And E_trans For Fragment 2", "Count", "Internal Energy (cm-1)", "Translational Energy (cm-1)", HISTOGRAM_2D_SIZE, HISTOGRAM_2D_SIZE, limits.E_max, limits.E_max2);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result_Data::reset() {
    grand_total=0;
    E_trans.clear();
    frag1_vel.clear();
    frag2_vel.clear();
    frag1_E_trans.clear();
    frag2_E_trans.clear();
    frag1_E_int.clear();
    frag2_E_int.clear();
    frag1_E_vib.clear();
    frag2_E_vib.clear();
    frag1_E_rot.clear();
    frag2_E_rot.clear();
    frag1_v.clear();
    frag2_v.clear();
    frag1_nk.clear();
    frag2_nk.clear();
    frag1_E_int_trans.clear();
    frag2_E_int_trans.clear();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result_Data::write_totals(FILE *the_file) {
    fprintf(the_file, "Total States Counted %lf\n", grand_total);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result_Data::write_histograms(FILE *the_file) {
    write_totals(the_file);
    if (E_trans.on) E_trans.write_histogram(the_file);
    if (frag1_vel.on) frag1_vel.write_histogram(the_file);
    if (frag1_E_trans.on) frag1_E_trans.write_histogram(the_file);
    if (frag1_E_int.on) frag1_E_int.write_histogram(the_file);
    if (frag1_E_vib.on) frag1_E_vib.write_histogram(the_file);
    if (frag1_v.on) frag1_v.write_histogram(the_file);
    if (frag1_E_rot.on) frag1_E_rot.write_histogram(the_file);
    if (frag1_nk.on) frag1_nk.write_histogram(the_file);
    if (frag2_vel.on) frag2_vel.write_histogram(the_file);
    if (frag2_E_trans.on) frag2_E_trans.write_histogram(the_file);
    if (frag2_E_int.on) frag2_E_int.write_histogram(the_file);
    if (frag2_E_vib.on) frag2_E_vib.write_histogram(the_file);
    if (frag2_v.on) frag2_v.write_histogram(the_file);
    if (frag2_E_rot.on) frag2_E_rot.write_histogram(the_file);
    if (frag2_nk.on) frag2_nk.write_histogram(the_file);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result_Data::write_bitmaps(const std::string job_path_and_name_start) {
    if (frag1_E_int_trans.on) frag1_E_int_trans.write_bitmap(job_path_and_name_start);
    if (frag2_E_int_trans.on) frag2_E_int_trans.write_bitmap(job_path_and_name_start);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




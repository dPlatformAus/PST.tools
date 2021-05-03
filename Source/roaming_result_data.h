
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
class Roaming_Histogram_Use : public Histogram_Use
{
    public:
        bool r_frag1_nk, r_frag2_nk;

        Roaming_Histogram_Use();
        virtual ~Roaming_Histogram_Use() {}
        virtual Roaming_Histogram_Use* clone(){return new Roaming_Histogram_Use(*this);};
        virtual void read_histogram_use(const std::string &key, const std::string &value);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){
            ar & boost::serialization::base_object<Histogram_Use>(*this);
            ar & r_frag1_nk & r_frag2_nk;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Histogram_Use::Roaming_Histogram_Use() : Histogram_Use(){
    r_frag1_nk=r_frag2_nk=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Histogram_Use::read_histogram_use(const std::string &key, const std::string &value) {
    Histogram_Use::read_histogram_use(key, value); //first check if the key belongs to the parent
    if (key == "r_frag1_nk") frag1_nk=atoi(value.c_str());
    else if (key == "r_frag2_nk") frag2_nk=atoi(value.c_str());
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
class Roaming_Result_Data : public Phasespace_Result_Data
{
    public:
        Histogram_2d r_frag1_nk, r_frag2_nk;
        double num_radical, num_may_roam;

        virtual ~Roaming_Result_Data() {}
        Roaming_Result_Data& operator=(const Roaming_Result_Data &rhs);
        Roaming_Result_Data& operator+=(const Roaming_Result_Data &rhs);
        virtual void initialise(const Histogram_Limits &limits, const Roaming_Histogram_Use *is_on);
        virtual void reset();
        virtual void write_totals(FILE *the_file);
        virtual void write_histograms(FILE *the_file);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){
            ar & boost::serialization::base_object<Phasespace_Result_Data>(*this);
            ar & r_frag1_nk & r_frag2_nk & num_radical & num_may_roam;
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Result_Data& Roaming_Result_Data::operator=(const Roaming_Result_Data &rhs){
    Phasespace_Result_Data::operator=(rhs);
    num_radical=rhs.num_radical;
    num_may_roam=rhs.num_may_roam;
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Result_Data& Roaming_Result_Data::operator+=(const Roaming_Result_Data &rhs){
    Phasespace_Result_Data::operator+=(rhs);
    num_radical+=rhs.num_radical;
    num_may_roam+=rhs.num_may_roam;
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Result_Data::initialise(const Histogram_Limits &limits, const Roaming_Histogram_Use *is_on) {
    Phasespace_Result_Data::initialise(limits, is_on);
    //probably should do something about the names! and the .on check should almost certainly be handled internally by Histogram!!!
    num_radical=num_may_roam=0;
    if (is_on->r_frag1_nk) r_frag1_nk.initialiseS("Roaming J And K For Fragment 1", "Count", "J", "k", 1, 1, limits.J_max1, limits.J_max1);
    if (is_on->r_frag2_nk) r_frag2_nk.initialiseS("Roaming J And K For Fragment 2", "Count", "J", "k", 1, 1, limits.J_max2, limits.J_max2);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Result_Data::reset() {
    Phasespace_Result_Data::reset();
    num_radical=num_may_roam=0;
    r_frag1_nk.clear();
    r_frag2_nk.clear();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Result_Data::write_totals(FILE *the_file) {
    Phasespace_Result_Data::write_totals(the_file);
    fprintf(the_file, "Total Radical States Counted %lf\n", num_radical);
    fprintf(the_file, "Total Potential Roaming States Counted %lf\n", num_may_roam);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Result_Data::write_histograms(FILE *the_file) { //write totals is called first, by Phasespace_Result_Data::write_histograms();, then the Phasespace_Result_Data histograms write .. then these do
    Phasespace_Result_Data::write_histograms(the_file);
    if (r_frag1_nk.on) r_frag1_nk.write_histogram(the_file);
    if (r_frag2_nk.on) r_frag2_nk.write_histogram(the_file);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

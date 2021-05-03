//***********************************************************************************************************************************************************************************
class Phasespace_Result
{
    public:
        Phasespace_Result_Data *data; ///THIS WILL NEED TO BE A POINTER SO THAT DERIVED VERSIONS CAN BE USED IN Phasespace_Result DERIVATIVES (WHICH WE NEED FOR THE ADDITIONAL ROAMING HISTOGRAMS?)
        std::vector<Phasespace_Probe*> probes; ///< list of the states to be probed

        virtual ~Phasespace_Result();
        Phasespace_Result(){create_data();}
        Phasespace_Result(const Phasespace_Result &rhs);
        Phasespace_Result& operator=(const Phasespace_Result &rhs);
        Phasespace_Result& operator+=(const Phasespace_Result &rhs);
        virtual void initialise(const Histogram_Limits &limits, const Phasespace_Molecule *reactant);
        virtual void create_data();
        virtual void create_probes(const Histogram_Limits &limits, const Phasespace_Molecule *reactant);
        virtual void reset();
        virtual void write_histograms(FILE *the_file);
        virtual void write_bitmaps(const std::string job_path_and_name_start);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){
            ar & data; /// IF this fails because the derived members aren't serialized (due to the inherited base class pointer for data, analogous to the reason we dynamically cast the molecule pointers etc).
            ar & probes;/// Then consider NOT calling the base_class serialize and instead repeating the serialization of probes and data in the derivatives serialize function but with dynamically cast pointers ...?
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Result::~Phasespace_Result(){
    unsigned i;
    for (i=0;i<probes.size();i++) delete probes[i];
    delete data;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Result::Phasespace_Result(const Phasespace_Result &rhs):probes(rhs.probes.size()){
    unsigned i;
    for (i=0;i<rhs.probes.size();i++) probes[i]=rhs.probes[i]->clone();
    data=rhs.data->clone();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Result& Phasespace_Result::operator=(const Phasespace_Result &rhs){
    unsigned i;
    for (i=0;i<probes.size();i++) delete probes[i];
    probes.resize(rhs.probes.size());
    for (i=0;i<rhs.probes.size();i++) probes[i]=rhs.probes[i]->clone();
    data=rhs.data->clone();
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Result& Phasespace_Result::operator+=(const Phasespace_Result &rhs){
    unsigned i;
    for (i=0;i<probes.size();i++) {
        (*probes[i])+=(*rhs.probes[i]);
    }
    (*data)+=(*rhs.data);
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result::initialise(const Histogram_Limits &limits, const Phasespace_Molecule *reactant) {
    unsigned i;
    for (i=0;i<probes.size();i++) delete probes[i];
    probes.clear();//clear any old data, in case this is a re-initialisation
    data->initialise(limits, reactant->is_on); //this may also need to be a function for the inheritance to work (call from dynamically cast pointer)...
    create_probes(limits, reactant);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result::create_data() {
    data = new Phasespace_Result_Data();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result::create_probes(const Histogram_Limits &limits, const Phasespace_Molecule *reactant) {
    unsigned i;
    for(i=0; i<reactant->probe_states.size(); i++) probes.push_back(new Phasespace_Probe(reactant->probe_states[i], limits, reactant->probe_is_on));
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result::reset() {
    unsigned i;
    for (i=0;i<probes.size();i++) probes[i]->reset();
    data->reset();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result::write_histograms(FILE *the_file) {
    unsigned p;
    data->write_histograms(the_file);
    for (p=0; p<probes.size(); p++) probes[p]->write_histograms(the_file);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Result::write_bitmaps(const std::string job_path_and_name_start) {
    unsigned p;
    std::string new_job_path_and_name_start;
    new_job_path_and_name_start=job_path_and_name_start+"All_";
    data->write_bitmaps(new_job_path_and_name_start);
    for (p=0; p<probes.size(); p++) {
        new_job_path_and_name_start=job_path_and_name_start+probes[p]->state.name+"_";
        probes[p]->data->write_bitmaps(new_job_path_and_name_start);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



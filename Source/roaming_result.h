//***********************************************************************************************************************************************************************************
class Roaming_Result : public Phasespace_Result
{
    public:
        Roaming_Result_Data *r_data;
        std::vector<Roaming_Probe*> r_probes;

        virtual ~Roaming_Result(){}//no new objects ... data and probes should be deleted by ~Phasespace_Result()
        Roaming_Result(){create_data();} //only difference is create the derived data type
        virtual void create_data();
        virtual void create_data(const Histogram_Limits &limits, const Roaming_Molecule *reactant);
        virtual void create_probes(const Histogram_Limits &limits, const Roaming_Molecule *reactant);

        /*Roaming_Result(const Roaming_Result &rhs);
        Roaming_Result& operator=(const Roaming_Result &rhs);
        Roaming_Result& operator+=(const Roaming_Result &rhs);
        virtual void initialise(const Histogram_Limits &limits, const Roaming_Molecule *reactant);
        virtual void reset();*/
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & r_data; /// Just doing it (due to the inherited base class pointer for data, analogous to the reason we dynamically cast the molecule pointers etc).
            ar & r_probes;/// NOT calling the base_class serialize and instead repeating the serialization of probes and data in this derivative serialize function but with dynamically cast pointers
        }
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Result::create_data(){
    data = new Roaming_Result_Data();
    r_data = dynamic_cast<Roaming_Result_Data*>(data);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Result::create_data(const Histogram_Limits &limits, const Roaming_Molecule *reactant) {
    create_data();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Roaming_Result::create_probes(const Histogram_Limits &limits, const Roaming_Molecule *reactant) {
    unsigned i;
    for(i=0; i<reactant->probe_states.size(); i++) {
        probes.push_back(new Roaming_Probe(reactant->probe_states[i], limits, reactant->r_probe_is_on));
        r_probes.push_back(dynamic_cast<Roaming_Probe*>(probes[i]));
    }
}/*
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Result::Roaming_Result(const Roaming_Result &rhs):Phasespace_Result(rhs){
    grand_total_no_tunneling_probability=rhs.grand_total_no_tunneling_probability;
    total_parent_states_considered=rhs.total_parent_states_considered;
#if (defined(USING_MPI))
    mpi::communicator world;
//    std::cout<<"| Process "<<world.rank()<<" : IN THE Secondary_Dissociation_Result COPY CONSTRUCTOR!"<<std::endl;
#endif
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Secondary_Dissociation_Result& Secondary_Dissociation_Result::operator=(const Secondary_Dissociation_Result &rhs){
    Phasespace_Result::operator=(rhs);
    grand_total_no_tunneling_probability=rhs.grand_total_no_tunneling_probability;
    total_parent_states_considered=rhs.total_parent_states_considered;
#if (defined(USING_MPI))
    mpi::communicator world;
//    std::cout<<"| Process "<<world.rank()<<" : IN THE Secondary_Dissociation_Result ASSIGNMENT OPERATOR!"<<std::endl;
#endif
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Secondary_Dissociation_Result& Secondary_Dissociation_Result::operator+=(const Secondary_Dissociation_Result &rhs){
    Phasespace_Result::operator+=(rhs);
    grand_total_no_tunneling_probability+=rhs.grand_total_no_tunneling_probability;
    total_parent_states_considered+=rhs.total_parent_states_considered;
#if (defined(USING_MPI))
    mpi::communicator world;
//    std::cout<<"| Process "<<world.rank()<<" : IN THE Secondary_Dissociation_Result INCREMENT OPERATOR!"<<std::endl;
#endif
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Secondary_Dissociation_Result::initialise(const Histogram_Limits &limits, const Phasespace_Molecule *reactant) {
    Phasespace_Result::initialise(limits, reactant);
    grand_total_no_tunneling_probability=0;
    total_parent_states_considered=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Secondary_Dissociation_Result::reset(){
    Phasespace_Result::reset();
    grand_total_no_tunneling_probability=0;
    total_parent_states_considered=0;
}*/
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

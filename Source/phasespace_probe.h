struct Fragment_Quantum_State {
    unsigned v, n, k, so;
    double E_nk;
    Fragment_Quantum_State():v(0),n(0),k(0),so(0),E_nk(0) {} //all defaults are 0
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Double_Quantum_State{
    bool tunneled;
    Fragment_Quantum_State f1, f2;
    unsigned L;
    double E_translation, degeneracy;
    Double_Quantum_State():tunneled(0),f1(),f2(),L(0),E_translation(0),degeneracy(0) {} //all defaults are 0
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
struct Processed_State_Data{
    double E_trans, vel_f1, vel_f2, E_trans_f1, E_trans_f2, E_vib_f1, E_vib_f2, E_int_f1, E_int_f2;
    Processed_State_Data():vel_f1(0),vel_f2(0),E_trans_f1(0),E_trans_f2(0),E_vib_f1(0),E_vib_f2(0),E_int_f1(0),E_int_f2(0) {} //all defaults are 0
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
class Probe_State {
    public:
        int fragment; // identify which fragment state we are probing 1 or 2 // later, I think we'll get rid of this and store the probe states under the fragment being probed
        Number_Range<int> vibRange, nRange, kRange;
        Number_Range<double> soRange; //if unspecified, default 0,0, YET TO WORK THIS ONE OUT, perhaps interpret spin range as [ (nRange.min + soRange.min) -> (nRange.max + soRange.max) ]
        std::string name;
        Probe_State();
        void reset();
        bool readState(const std::string line);
        bool satisfies_criteria(const Fragment_Quantum_State &state);
        bool applies(const Double_Quantum_State &state);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & fragment & vibRange & nRange & kRange & soRange & name;}
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Probe_State::Probe_State(){
    reset();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Probe_State::reset(){
    fragment=0;
    vibRange.setRange(-1); //we use -1 to indicate any value will be accepted for this quantum number
    nRange.setRange(-1);
    kRange.setRange(-1);
    soRange.setRange(-1);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Probe_State::readState(const std::string line){
    bool result = true;
    std::stringstream ss;
    std::string value;
    reset();
    ss << line;
    if (getline(ss,value,',')) {
        fragment=atoi(value.c_str()); //set fragment id
    }else{
        result=false;
    }
    if (getline(ss,value,',')) {
        vibRange.setRange(value); //set vibrational state range
    }else{
        result=false;
    }
    if (getline(ss,value,',')) {
        nRange.setRange(value); //set rotational n state range
    }else{
        result=false;
    }
    if (getline(ss,value,',')) {
        kRange.setRange(value); //set rotational k state range
    }else{
        result=false;
    }
    if (getline(ss,value,',')) {
        soRange.setRange(value); //set spin orbit coupling state spin range
    }//spin orbit range not mandatory don't return false if not entered
    name="F"+boost::lexical_cast<std::string>(fragment);
    name+="_v"+boost::lexical_cast<std::string>(vibRange.min)+";"+boost::lexical_cast<std::string>(vibRange.max);
    name+="_n"+boost::lexical_cast<std::string>(nRange.min)+";"+boost::lexical_cast<std::string>(nRange.max);
    name+="_k"+boost::lexical_cast<std::string>(kRange.min)+";"+boost::lexical_cast<std::string>(kRange.max);
    name+="_so"+boost::lexical_cast<std::string>(soRange.min)+";"+boost::lexical_cast<std::string>(soRange.max);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Probe_State::satisfies_criteria(const Fragment_Quantum_State &state){ ///CURRENTLY NOT IMPLEMENTING SPIN ORBITAL PART ... I SUSPECT WE NEED TO IMPLEMENT HUND'S CASE B FIRST
    bool result(true);
    if (vibRange.min >=0  &&  vibRange.max >=0  &&  (state.v < (unsigned)vibRange.min || state.v > (unsigned)vibRange.max) ) result=false;
    if (nRange.min >=0  &&  nRange.max >=0  &&  (state.n < (unsigned)nRange.min || state.n > (unsigned)nRange.max) ) result=false;
    if (kRange.min >=0  &&  kRange.max >=0  &&  (state.k < (unsigned)kRange.min || state.k > (unsigned)kRange.max) ) result=false;
    //if (state.n == nRange.min && state.so <...  YET TO SORT OUT HOW TO DO THIS SINCE SO DEGEN INCORPORATES ALL 2J+1 SPIN VALUES, RATHER THEN LOOPING THROUGH ALL SPINS ... needs some more thought, somehow need to subtract out the excess degen
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Probe_State::applies(const Double_Quantum_State &state){
    bool result;
    if (fragment == 1) result=satisfies_criteria(state.f1);
    else result=satisfies_criteria(state.f2); // Double_Quantum_State only has 2 fragments so if it's not 1 it must be 2
    return result;
}
//***********************************************************************************************************************************************************************************
class Phasespace_Probe {
    public:
        Probe_State state;
        Phasespace_Result_Data *data;

        virtual ~Phasespace_Probe() {delete data;}
        Phasespace_Probe(){data = new Phasespace_Result_Data();}
        Phasespace_Probe(const Probe_State &in_state, const Histogram_Limits &limits, const Histogram_Use *probe_is_on);
        virtual Phasespace_Probe* clone(){return new Phasespace_Probe(*this);};
        Phasespace_Probe(const Phasespace_Probe &rhs);
        Phasespace_Probe& operator+=(const Phasespace_Probe &rhs);
        void reset();
        void write_histograms(FILE *the_file);
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & state & data;}
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Probe::Phasespace_Probe(const Probe_State &in_state, const Histogram_Limits &limits, const Histogram_Use *probe_is_on){
    state=in_state;
    data = new Phasespace_Result_Data();
    data->initialise(limits, probe_is_on);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Probe::Phasespace_Probe(const Phasespace_Probe &rhs){
    state=rhs.state;
    data=rhs.data->clone();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Phasespace_Probe& Phasespace_Probe::operator+=(const Phasespace_Probe &rhs){

#if (defined(USING_MPI))
    mpi::communicator world;
    std::cout<<"| Process "<<world.rank()<<" : IN THE Phasespace_Probe Increment Operator!"<<std::endl;
#endif
    (*data)+=(*rhs.data);
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Probe::reset() {
    data->reset();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Phasespace_Probe::write_histograms(FILE *the_file) {
    fprintf(the_file, "\n\n\n%s\n\n", state.name.c_str());
    data->write_histograms(the_file);
}
//***********************************************************************************************************************************************************************************

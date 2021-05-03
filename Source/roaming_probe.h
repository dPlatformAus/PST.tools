//***********************************************************************************************************************************************************************************
class Roaming_Probe : public Phasespace_Probe
{
    public:
        Roaming_Result_Data *r_data;
        //double num_radical, num_may_roam; ///BUT INT HE DATA NOW!  specific to roaming (in the dtata now?)

        Roaming_Probe(){data = new Roaming_Result_Data();}
        Roaming_Probe(const Probe_State &in_state, const Histogram_Limits &limits, const Roaming_Histogram_Use *probe_is_on);
        virtual Roaming_Probe* clone(){return new Roaming_Probe(*this);};
        //Roaming_Probe(const Roaming_Probe &rhs);
        virtual ~Roaming_Probe(){};
/*    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Probe>(*this);
            ar & num_radical & num_may_roam;
        }*/
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Roaming_Probe::Roaming_Probe(const Probe_State &in_state, const Histogram_Limits &limits, const Roaming_Histogram_Use *probe_is_on){
    state=in_state;
    data = new Roaming_Result_Data(); //do we need a dynamic cast here? ... if so, we'll need a new pointer for it and probably need to override (uncomment/edit) serialize stuff too
    r_data = dynamic_cast<Roaming_Result_Data*>(data);
    r_data->initialise(limits, probe_is_on);
}

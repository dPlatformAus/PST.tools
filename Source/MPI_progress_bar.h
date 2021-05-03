
unsigned dMPI_bar_number(0);

inline unsigned dMPI_get_progress_bar_tag(const int process){
    return (process*MPIPB_PROCESS_NUMBER_FACTOR)+(job_number*MPIPB_JOB_NUMBER_FACTOR)+dMPI_bar_number;
    ///we should do an experiment or two to see if this needs to be so stringent in the single isend message mode (which I think was the whole problem)
    //return dMPI_bar_number*MPIPB_SEND_PROGRESS+process; //this is required because boost::mpi isend and irecv get mixed up if you send multiple requests with the same tag (even though this is allowed with MPI)
    ///however, it concerns me that there may be something the mpi progress bars aren't cleaning up that made this nessesary ... and that this might cause a memory leak because all the requests are being left unfulfilled?
    ///I seem to remember doing some reading about this and finding out that it was to do with the 2 posts boost must do to complete a single post (serial data length, then serialised data) ... better review and take some notes!
}
//***********************************************************************************************************************************************************************************//***********************************************************************************************************************************************************************************
class dMPI_Slave_Progress_Bar_State{
    public:
        unsigned job_length, progress, position, frame;
        bool finnished;

        dMPI_Slave_Progress_Bar_State();
        void reset();
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & job_length & progress & position & frame & finnished;}
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dMPI_Slave_Progress_Bar_State::dMPI_Slave_Progress_Bar_State(){
    reset();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void dMPI_Slave_Progress_Bar_State::reset(){
    job_length=progress=position=frame=0;
    finnished=0;
}
//***********************************************************************************************************************************************************************************//***********************************************************************************************************************************************************************************
//This line tells boost::mpi that the datatype dMPI_Slave_Progress_Bar_State is fixed size  and can be mapped to an mpi data type
//this means boost only sends a single isend message with the data, preventing the segmentation faults due to message overtaking on the same tag id with the 2 message skeleton mode
BOOST_IS_MPI_DATATYPE(dMPI_Slave_Progress_Bar_State);
//***********************************************************************************************************************************************************************************//***********************************************************************************************************************************************************************************
class dMPI_Master_Progress_Bar_Slave_Node{
    public:
        dMPI_Slave_Progress_Bar_State state;
        mpi::request irecv_request;
};
//***********************************************************************************************************************************************************************************//***********************************************************************************************************************************************************************************
class dMPI_Master_Progress_Bar : public Progress_Bar
{
    public:
        mpi::communicator world;
        std::vector<dMPI_Master_Progress_Bar_Slave_Node> slaves;

        void initialise(const unsigned in_bar_size, const int in_update_interval, const unsigned in_job_length, const unsigned caller, const unsigned bar_number);
        void un_initialise(const unsigned caller);
        std::string get_slave_bars_string();
        virtual void print();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void dMPI_Master_Progress_Bar::initialise(const unsigned in_bar_size, const int in_update_interval, const unsigned in_job_length, const unsigned caller, const unsigned bar_number){
    Progress_Bar::initialise(in_bar_size, in_update_interval, in_job_length, caller);
    //the_progress_bar.block(); // block the global progress bar so that long state counts can't start it (that would mess up this MPI progress bar's output)
    slaves.clear();
    slaves.resize(world.size()-1);
    dMPI_bar_number=bar_number;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void dMPI_Master_Progress_Bar::un_initialise(const unsigned caller){
    unsigned i;
    bool all_processes_not_finnished;
    std::string slave_bars_string;
    do{
        slave_bars_string=get_slave_bars_string();
        print_semipersistent_output();
        std::cout<<"Process 0:|JOB COMPLETE\n"<<slave_bars_string.c_str();
        all_processes_not_finnished=false;
        for (i=0; i<slaves.size(); i++) if (slaves[i].state.finnished==false) all_processes_not_finnished=true;
        if (all_processes_not_finnished) sleep(update_interval);
    }while (all_processes_not_finnished);
    Progress_Bar::un_initialise(caller);
    for (i=0; i<slaves.size(); i++){
        if (!slaves[i].irecv_request.test()) slaves[i].irecv_request.cancel();
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::string dMPI_Master_Progress_Bar::get_slave_bars_string() {
    unsigned i;
    std::ostringstream bar_stream;
    std::string slave_bars_string;
    for (i=0; i<slaves.size(); i++){
        if (slaves[i].state.finnished==false){
            if (slaves[i].irecv_request.test()) {
                slaves[i].irecv_request=world.irecv(i+1, dMPI_get_progress_bar_tag(i+1), slaves[i].state);
            }

        }
        slave_bars_string="";
        if (slaves[i].state.finnished==false) slave_bars_string=get_bar_string(slaves[i].state.progress, slaves[i].state.job_length, slaves[i].state.position, slaves[i].state.frame);
        else slave_bars_string="|JOB COMPLETE\n";
        bar_stream<<"Process "<<i+1<<":"<<slave_bars_string.c_str();
    }
    return bar_stream.str();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void dMPI_Master_Progress_Bar::print(){
    std::string bar_string=get_bar_string(initialiser_progress, job_length, position, frame);
    std::string slave_bars_string=get_slave_bars_string();
    print_semipersistent_output();
    std::cout<<"Process 0:"<<bar_string.c_str()<<slave_bars_string.c_str();
}
//***********************************************************************************************************************************************************************************
class dMPI_Slave_Progress_Bar : public Progress_Bar
{
    public:
        mpi::communicator world;
        mpi::request isend_request;
        dMPI_Slave_Progress_Bar_State state;

        void initialise(const unsigned in_bar_size, const int in_update_interval, const unsigned in_job_length, const unsigned caller, const unsigned bar_number);
        void un_initialise(const unsigned caller);
        virtual void print();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void dMPI_Slave_Progress_Bar::initialise(const unsigned in_bar_size, const int in_update_interval, const unsigned in_job_length, const unsigned caller, const unsigned bar_number){
    Progress_Bar::initialise(in_bar_size, in_update_interval, in_job_length, caller);
    dMPI_bar_number=bar_number;
    state.reset();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void dMPI_Slave_Progress_Bar::un_initialise(const unsigned caller){
    mpi::status wait_result;
    Progress_Bar::un_initialise(caller);
    state.finnished=true;
    if (!isend_request.test()){
        isend_request.cancel();
    }
    isend_request=world.isend(0, dMPI_get_progress_bar_tag(world.rank()), state);
    wait_result=isend_request.wait();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void dMPI_Slave_Progress_Bar::print(){
    state.job_length=job_length;
    state.progress=initialiser_progress;
    state.position=position;
    state.frame=frame;
    if (!isend_request.test()){ //if the previous isend_request has not been received yet, cancel it before sending the new one
        isend_request.cancel();
    }
    isend_request=world.isend(0, dMPI_get_progress_bar_tag(world.rank()), state);
}
//***********************************************************************************************************************************************************************************


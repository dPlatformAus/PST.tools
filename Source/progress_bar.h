#define CORE_PROGRESS_BAR 1
#define CALCULATION_PROGRESS_BAR 2
class Progress_Bar
{
    public:
        bool initialised, blocked;
        unsigned job_length, bar_size, frame, position, initialiser, initialiser_progress;
        int last_time, update_interval;
        std::string frames;

        Progress_Bar();
        void initialise(const unsigned in_bar_size, const int in_update_interval, const unsigned in_job_length, const unsigned caller);
        void un_initialise(const unsigned caller);
        void block(){blocked=1;};
        void un_block(){blocked=0;};
        void go(const unsigned progress, const unsigned caller);
        std::string get_bar_string(const unsigned progress, const unsigned job_length, const unsigned bar_position, const unsigned bar_frame);
        virtual void print();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Progress_Bar::Progress_Bar(){
    blocked=initialised=0;
    initialiser_progress=job_length=bar_size=frame=position=initialiser=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Progress_Bar::initialise(const unsigned in_bar_size, const int in_update_interval, const unsigned in_job_length, const unsigned caller){
    if (!initialised){
        initialised=1;
        initialiser=caller;
        bar_size=in_bar_size;
        update_interval=in_update_interval;
        job_length=in_job_length;
        frame=position=0;
        frames="|/-\\";
        last_time=time(NULL);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Progress_Bar::un_initialise(const unsigned caller){
    if (initialiser==caller) initialised=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Progress_Bar::go(const unsigned progress, const unsigned caller) {
    int current_time;
    current_time=time(NULL);
    if (initialiser==caller && job_length) {//only change the position if the initialising process is calling go with the correct progress value, always update position so if printed by non-caller the position still progresses
        position=int(bar_size*(1.0*progress/job_length)); // leave the 1.0 in there so the rounding isn't integer rounding (divide first to avoid looping over back past 0)
        initialiser_progress=progress; // only print out the number passed by the initialiser
    }
    if (!blocked && update_interval<current_time-last_time) {
        last_time=current_time;
        if (frame==3) frame=0;
        else frame++;
        print();
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::string Progress_Bar::get_bar_string(const unsigned progress, const unsigned job_length, const unsigned bar_position, const unsigned bar_frame) {
    unsigned i;
    std::ostringstream bar_stream;
    if (bar_position<=bar_size){
        bar_stream<<"|";
        for (i=0; i<bar_position; i++) bar_stream << "*";
        bar_stream << frames[bar_frame];
        for (i=bar_position; i<bar_size; i++) bar_stream << "_";
        bar_stream<<"|";
        bar_stream<<"("<<progress<<"/"<<job_length<<")";
        bar_stream<<std::endl;
    }else bar_stream<<"Error: Progress_Bar::get_bar_string - bar_position > bar_size !!"<<std::endl;
    return bar_stream.str();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Progress_Bar::print() {
    std::string bar_string;
    if (job_length) bar_string=get_bar_string(initialiser_progress, job_length, position, frame);
    print_semipersistent_output();
    std::cout<<bar_string.c_str();
}

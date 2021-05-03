//***********************************************************************************************************************************************************************************
//This class encapsulates all the tasks that are common to ALL calculations!
class Base_Calculation{
    public:
        std::string input_file;
        Output_File out_file;
        std::string job_path;

        Base_Calculation(){}; //C++ doesn't do virtual constructors, so we can't do much with this
        virtual ~Base_Calculation();
        virtual void create_output_directory_and_file();
        void rename_input_file(const std::string &file_name);
        virtual void initialise(const std::string &file_name);
        virtual bool run()=0;
        void go(const std::string &file_name);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Base_Calculation::~Base_Calculation(){
    out_file.close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Base_Calculation::create_output_directory_and_file(){
    std::string csv_name, job_name;
    job_name = input_file;
    job_name.erase(job_name.find_last_of("."), std::string::npos);
    job_path=output_directory+job_name+"/";
    out_file.set_output_directory(job_path);
    csv_name=job_name;
    csv_name.append(".csv");
    out_file.open(csv_name);
    print_persistent_output();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Base_Calculation::initialise(const std::string &file_name){
    input_file=file_name;
    create_output_directory_and_file();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Base_Calculation::rename_input_file(const std::string &file_name){ //only called if the job calculation ran sucessfully
    std::string old_file_name, new_file_name;
    size_t found;
    int status;
    if (RENAME_FILES_ON_COMPLETION) { //if we want files renamed if the job calculation ran sucessfully, rename the input file from .in to .done
        old_file_name=new_file_name=input_directory+file_name;
        found=new_file_name.rfind(input_suffix);
        new_file_name.replace(found,input_suffix.length(),done_suffix);
        status=rename(old_file_name.c_str(), new_file_name.c_str());
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Base_Calculation::go(const std::string &file_name){
    //(later, create an errors list that is passed by reference to calculations that will contain the errors that occured, or something)
    boost::timer the_timer;//start timing execution time
    initialise(file_name); //parse the input file, etc
    if (run()){
        persistent_output<<" - Executed in "<<pretty_time(the_timer.elapsed()).c_str()<<std::endl; //print execution time to screen
        rename_input_file(file_name);
    }
    //else handle_error(); //there is lots of work to do, providing proper error checking etc!*/
}


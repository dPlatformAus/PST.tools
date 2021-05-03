void do_slave_execution(){
    mpi::communicator world;
    unsigned i, num_files;
    the_progress_bar.block(); // block the global progress bar so that long state counts can't start it (slave processes cannot print to screen ... may need to have a silent variable or something that kills all output)
    std::string calc_type("");
    broadcast(world, num_files, 0);
    for (i=0; i<num_files; i++){
        broadcast(world, job_number, 0);
        broadcast(world, calc_type, 0);
        if (calc_type == "phasespace") Phasespace_Calculation_Slave_Process();
        else if (calc_type == "impulsive") Impulsive_Calculation_Slave_Process();
        else if (calc_type == "triple fragmentation") Triple_Fragmentation_Calculation_Slave_Process();
    }
}


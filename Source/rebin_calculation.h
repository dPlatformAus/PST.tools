//***********************************************************************************************************************************************************************************
struct XY_Data_Point{
    double x, y;
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool rebin_source_sort_predicate(const XY_Data_Point &lhs, const XY_Data_Point &rhs) {return lhs.x < rhs.x;}
//***********************************************************************************************************************************************************************************
class Rebin_Calculation : public Base_Calculation{
    public:
        unsigned method, specify_shape, bin_length, convolving_function;
        const gsl_interp_type *gsl_interpolation_type;
        double bin_size, new_x_min, new_x_max, convolving_function_width;
        std::vector<XY_Data_Point> data;
        std::string histogram_title, x_title, y_title;

        Rebin_Calculation(){}; //C++ doesn't do virtual constructors, so we can't do much with this
        Rebin_Calculation(const std::string &file_name){go(file_name);};
        virtual ~Rebin_Calculation(){};
        virtual void initialise(const std::string &file_name);
        virtual bool run();
        bool simple_rebin();
        bool interpolating_rebin();
        bool convolving_rebin();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Rebin_Calculation::initialise(const std::string &file_name){
    Base_Calculation::initialise(file_name);

    const unsigned NO_SECTION=0;
    const unsigned GENERAL=1;
    const unsigned DATA=2;
	unsigned read_state=NO_SECTION;
	size_t comma_position;
	std::string file_name_and_path=input_directory+file_name;
    std::ifstream infile (file_name_and_path.c_str());
	std::string line, key, value;
    std::string x_string, y_string;
    XY_Data_Point temp_data_point;
    unsigned interpolation_type(0); //default interpolation method is linear
    convolving_function=new_x_min=new_x_max=0;
	if (infile.is_open()) {
        while(getline(infile,line)) {
            if (line.substr(0,2) != "//"){ //if this is not a comment line (we should ignore comment lines)
                if (line.substr(0,8) == "--------"){ //this is the start of a new section
                    read_state=NO_SECTION; //set to none first in case the section is not recognised by this parse routine
                    if (line.find("-GENERAL-") != std::string::npos) read_state=GENERAL;
                    if (line.find("-DATA-") != std::string::npos) read_state=DATA;
                }
                if (read_state==GENERAL){
                    read_key_value(line, key, value);
                    if (key != "") {
                        if (key == "method") method=atoi(value.c_str());
                        else if (key == "specify_shape") specify_shape=atoi(value.c_str());
                        else if (key == "bin_length") bin_length=atoi(value.c_str());
                        else if (key == "bin_size") bin_size=atof(value.c_str());
                        else if (key == "new_x_min") new_x_min=atof(value.c_str());
                        else if (key == "new_x_max") new_x_max=atof(value.c_str());
                        else if (key == "histogram_title") histogram_title=value;
                        else if (key == "x_title") x_title=value;
                        else if (key == "y_title") y_title=value;
                        else if (key == "interpolation_type") interpolation_type=atoi(value.c_str());
                        else if (key == "convolving_function") convolving_function=atoi(value.c_str());
                        else if (key == "convolving_function_width") convolving_function_width=atof(value.c_str());
                    }
                }else if (read_state==DATA){
                    comma_position = line.find(",");
                    if (comma_position != std::string::npos){
                        x_string=line.substr(0, comma_position);
                        y_string=line.substr(comma_position + 1);
                        temp_data_point.x=atof(x_string.c_str());
                        temp_data_point.y=atof(y_string.c_str());
                        data.push_back(temp_data_point);
                    }
                }
            }
        }
        infile.close();
	}
	if (method==1){
        if (interpolation_type==0) gsl_interpolation_type=gsl_interp_linear;
        if (interpolation_type==1) gsl_interpolation_type=gsl_interp_cspline;
        if (interpolation_type==2) gsl_interpolation_type=gsl_interp_cspline_periodic;
        if (interpolation_type==3) gsl_interpolation_type=gsl_interp_akima;
        if (interpolation_type==4) gsl_interpolation_type=gsl_interp_akima_periodic;
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Rebin_Calculation::run(){
    persistent_output << " - ReBin_Calculation:";
    if (method==0) return simple_rebin();
    if (method==1) return interpolating_rebin();
    if (method==2) return convolving_rebin();
    return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Rebin_Calculation::simple_rebin(){
    unsigned i,j;
    double x_min, x_max, xi, l_bound_i, u_bound_i, l_bound_j, u_bound_j, half_w, overlap_weight;
    Histogram result;
    persistent_output << " Simple"<< std::endl;
    std::sort(data.begin(), data.end(), rebin_source_sort_predicate); //SORT THE INPUT DATA BY x ASCENDING!
    x_min=x_max=data[0].x;
    for (i=0;i<data.size();i++){ //find the extents of the source data (which will be used as the extents of the output unless new, valid, extents have been specified)
        if (data[i].x<x_min) x_min=data[i].x;
        if (data[i].x>x_max) x_max=data[i].x;
    }
    if (new_x_min<new_x_max && x_min<=new_x_min && new_x_max<=x_max){ //if new extents have been specified, use them (provided they are within the extents of the source data)
        x_min=new_x_min;
        x_max=new_x_max;
    }
    if (specify_shape==1) result.initialiseN(histogram_title, x_title, y_title, bin_length, x_max, x_min);
    if (specify_shape==2) result.initialiseS(histogram_title, x_title, y_title, bin_size, x_max, x_min); //NOTE: this method requires the min and max to be an integer multiple of bin_size (I think it will crash otherwise)
    half_w=0.5*result.bin_size;
    for (i=0;i<result.data.size();i++){
        xi=result.resolveIndex(i);
        l_bound_i=xi-half_w;
        u_bound_i=xi+half_w;
        l_bound_j=data[0].x;
        for (j=0;j<data.size();j++){
            if (j) l_bound_j=data[j-1].x+(data[j].x-data[j-1].x);
            if (j<data.size()-1) u_bound_j=data[j].x+(data[j+1].x-data[j].x);
            else u_bound_j=data[j].x;
            if (u_bound_j >= x_min && l_bound_j <= x_max) { //only count if within the extents of the new histogram
                if (l_bound_j >= l_bound_i && u_bound_j <= u_bound_i) { // source bin is entirely within this destination bin
                    result.increment_value(data[j].x, data[j].y); //simply increment
                }else if (l_bound_j > l_bound_i && l_bound_j < u_bound_i){ //only the lower bound of the source bin is within the bounds of this bin
                    overlap_weight=(u_bound_i-l_bound_j)/(u_bound_j-l_bound_j); //find the fraction of the source bin that overlaps this destination bin
                    result.increment_value(xi, overlap_weight*data[j].y); //increment this destination bin by the source bin value weighted by the overlap
                }else if (u_bound_j > l_bound_i && u_bound_j < u_bound_i){ //only the upper bound of the source bin is within the bounds of this bin
                    overlap_weight=(u_bound_j-l_bound_i)/(u_bound_j-l_bound_j); //find the fraction of the source bin that overlaps this destination bin
                    result.increment_value(xi, overlap_weight*data[j].y); //increment this destination bin by the source bin value weighted by the overlap
                }
            }
        }
    }
    result.write_histogram(out_file.fp);
    return 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Rebin_Calculation::interpolating_rebin(){
    unsigned i;
    double x_min, x_max, xi;
    double *x = new double[data.size()];
    double *y = new double[data.size()];
    Histogram result;
    persistent_output << " Interpolating"<< std::endl;
    std::sort(data.begin(), data.end(), rebin_source_sort_predicate); //SORT THE INPUT DATA BY x ASCENDING!
    x_min=x_max=data[0].x;
    for (i=0;i<data.size();i++){ //copy the data into the format required for gsl and find the extents of the source data (which will be used as the extents of the output unless new, valid, extents have been specified)
        if (data[i].x<x_min) x_min=data[i].x;
        if (data[i].x>x_max) x_max=data[i].x;
        x[i]=data[i].x;
        y[i]=data[i].y;
    }
    if (new_x_min<new_x_max && x_min<=new_x_min && new_x_max<=x_max){ //if new extents have been specified, use them (provided they are within the extents of the source data)
        x_min=new_x_min;
        x_max=new_x_max;
    }
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc (gsl_interpolation_type, data.size());
    gsl_spline_init (spline, x, y, data.size());
    if (specify_shape==1) result.initialiseN(histogram_title, x_title, y_title, bin_length, x_max, x_min);
    if (specify_shape==2) result.initialiseS(histogram_title, x_title, y_title, bin_size, x_max, x_min);
    for (i=0;i<result.data.size();i++){
        xi=result.resolveIndex(i);
        result.set_value(xi, gsl_spline_eval(spline, xi, acc));
    }
    result.write_histogram(out_file.fp);
    delete x;
    delete y;
    return 1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Rebin_Calculation::convolving_rebin(){
    unsigned i,j;
    double x_min, x_max, xi, norm_factor(0), other_factor(0), dx, dy(0);
    Histogram result;
    persistent_output << " Convolving"<< std::endl;
    x_min=x_max=data[0].x;
    for (i=0;i<data.size();i++){ //find the extents of the source data (which will be used as the extents of the output unless new, valid, extents have been specified)
        if (data[i].x<x_min) x_min=data[i].x;
        if (data[i].x>x_max) x_max=data[i].x;
    }
    if (new_x_min<new_x_max && x_min<=new_x_min && new_x_max<=x_max){ //if new extents have been specified, use them (provided they are within the extents of the source data)
        x_min=new_x_min;
        x_max=new_x_max;
    }
    if (specify_shape==1) result.initialiseN(histogram_title, x_title, y_title, bin_length, x_max, x_min);
    if (specify_shape==2) result.initialiseS(histogram_title, x_title, y_title, bin_size, x_max, x_min);
    if (convolving_function==0){ //Gaussian
        other_factor = -4*log(2)/pow(convolving_function_width,2);
        norm_factor = 2*sqrt(log(2))/(convolving_function_width*sqrt(PI));
    }else if(convolving_function==1){ //Lorentzian
        other_factor = pow(convolving_function_width,2)/4.0;
        norm_factor = convolving_function_width/(2*PI);
    }
    for (i=0;i<result.data.size();i++){
        xi=result.resolveIndex(i);
        for (j=0;j<data.size();j++){
            dx=xi-data[j].x;
            if (convolving_function==0) dy=norm_factor*data[j].y*exp(other_factor*pow(dx,2)); //Gaussian
            else if(convolving_function==1) dy=norm_factor*data[j].y/(pow(dx,2)+other_factor); //Lorentzian
            result.increment_value(xi, dy);
        }
    }
    result.write_histogram(out_file.fp);
    return 1;
}

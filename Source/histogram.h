#define GAUSSIAN_CONVOLVE 1
#define LORENTZIAN_CONVOLVE 2
unsigned HISTOGRAM_SIZE=100;
//***********************************************************************************************************************************************************************************
class Histogram
{
    public:
        unsigned bin_length;
        bool convolved, has_names, on;
        std::vector<std::string> bin_names;
        std::vector<double> data, Z_component, convolved_data;
        double x_max, x_min, bin_size, data_count, data_total, data_average, data_max, Z, T;
        std::string title, x_title, y_title, bin_names_title;

        Histogram();
        double& operator[] (unsigned i) {return data[i];}
        Histogram& operator+=(const Histogram &rhs);
        void clear();
        void initialiseN(const std::string& in_title, const std::string& in_x_title, const std::string& in_y_title, const unsigned in_length, const double in_x_max, const double in_x_min=0);
        void initialiseS(const std::string& in_title, const std::string& in_x_title, const std::string& in_y_title, const double in_size, const double in_x_max, const double in_x_min=0);
        void set_names(const std::string& in_bin_names_title, const std::vector<std::string> &in_names);
        void set_value(const double in_x, const double in_y);
        void increment_value(const double in_x, const double in_y);
        void range_increment(Histogram &in_hist, double scaling_factor);
        void inner_range_increment(Histogram &in_hist, double scaling_factor);
        unsigned getIndex(const double in_x);
        double resolveIndex(const unsigned in_index);
        void get_count();
        void get_max();
        void get_average();
        double get_normalised_area_prob(const unsigned i); //have to call get_count first!
        double get_normalised_height_prob(const unsigned i); //have to call get_max first!
        void write_average(FILE *the_file);
        void write_histogram(FILE *the_file);
        void write_Boltzmann(FILE *the_file); //only relevent if histogram has binned energies, which are assumed to be in cm-1 *** set_Boltzmann_temp must be called first!
        void set_Boltzmann_temp(const double in_T);
        double get_Boltzmann_prob(const unsigned i); //only relevent if histogram has binned energies, which are assumed to be in cm-1 *** set_Boltzmann_temp must be called first!
        void convolve(const double &width, const unsigned &method);
        void Gaussian_convolve(const double &width);
        void Lorentzian_convolve(const double &width);

    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){
            ar & bin_length & convolved & has_names & on & bin_names & data & Z_component & convolved_data & x_max & x_min;
            ar & bin_size & data_count & data_total & data_average & data_max & Z & T & title & x_title & y_title & bin_names_title;
        }
        void initialise(const std::string& in_title, const std::string& in_x_title, const std::string& in_y_title, const double in_x_max, const double in_x_min);
        double get_bin_size(const unsigned &in_length, const double &in_x_max, const double &in_x_min);
        unsigned get_bin_length(const double &in_size, const double &in_x_max, const double &in_x_min);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Histogram::Histogram(){
    bin_length=0; //gotta initialise these otherwise the = operator will go crazy when it tries to resize the array!!
    x_max=x_min=bin_size=data_count=data_total=data_average=data_max=Z=T=0;
    convolved=has_names=on=false;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Histogram& Histogram::operator+=(const Histogram &rhs){ //this function assumes that both histograms are defined as the same (ie, titles etc), it only checks that their dimensions are the same
    ///NOTE: no Boltzmann or convolution data is considered here currently!
    unsigned i; //currently this is only used by the Dissociation_Result += operator, so simply incrementing the data and recalculating the total count and average is considered sufficient
    if (bin_length==rhs.bin_length && bin_size==rhs.bin_size){
        for (i=0;i<bin_length;i++) data[i]+=rhs.data[i];
        get_average();
        get_max();
    }else{
        std::cout<<"Histogram::operator+=("<<title<<", "<<rhs.title<<") HISTOGRAMS INCOMPATIBLE!"<<std::endl;
        //persistent_output<<"Histogram::operator+=("<<title<<", "<<rhs.title<<") HISTOGRAMS INCOMPATIBLE!"<<std::endl;
        //print_persistent_output();
    }
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::initialise(const std::string& in_title, const std::string& in_x_title, const std::string& in_y_title, const double in_x_max, const double in_x_min){
    unsigned i;
    convolved=false;
    has_names=false;
    on=true;
    title=in_title;
    x_title=in_x_title;
    y_title=in_y_title;
    x_min=in_x_min;
    x_max=in_x_max;
    data.clear(); //in case histogram is reinitialised
    for (i=0;i<bin_length;i++) data.push_back(0);
    Z_component.clear();
    convolved_data.clear();
    bin_names.clear();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::initialiseN(const std::string& in_title, const std::string& in_x_title, const std::string& in_y_title, const unsigned in_length, const double in_x_max, const double in_x_min){
    bin_length=in_length;
    bin_size=get_bin_size(bin_length, in_x_max, in_x_min);
    initialise(in_title, in_x_title, in_y_title, in_x_max, in_x_min);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::initialiseS(const std::string& in_title, const std::string& in_x_title, const std::string& in_y_title, const double in_size, const double in_x_max, const double in_x_min){
    bin_size=in_size;
    bin_length=get_bin_length(bin_size, in_x_max, in_x_min);
    initialise(in_title, in_x_title, in_y_title, in_x_max, in_x_min);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::set_names(const std::string& in_bin_names_title, const std::vector<std::string> &in_names){
    unsigned i;
    if (bin_length <= in_names.size()){ //in_names must have at least as many names as the histogram has bins
        bin_names.clear();
        for (i=0; i<bin_length; i++) bin_names.push_back(in_names[i]);
        has_names=true;
        bin_names_title=in_bin_names_title;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram::get_bin_size(const unsigned &in_length, const double &in_x_max, const double &in_x_min){
    double result(1);
    if ((in_x_max-in_x_min)>0 && in_length>1) result=(in_x_max-in_x_min)/(in_length-1);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Histogram::get_bin_length(const double &in_size, const double &in_x_max, const double &in_x_min){
    unsigned result(1);
    if ((in_x_max-in_x_min)>0) result=(unsigned)ceil(1+(in_x_max-in_x_min)/in_size);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::clear(){
    unsigned i;
    for (i=0;i<bin_length;i++) data[i]=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Histogram::getIndex(const double in_x){
    unsigned result(0);
    //initialise shouldn't let bin_size=0, but checking is more robust
    if (bin_size) result=int(0.5 + (in_x-x_min)/bin_size);
    if (result>=bin_length) {
        persistent_output<<"Histogram::getIndex("<<in_x<<">="<<resolveIndex(bin_length-1)<<");"<<title<<" index out of range! - ("<<result<<">="<<bin_length<<")"<<std::endl;
        print_persistent_output();
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram::resolveIndex(const unsigned in_index){
    return x_min+in_index*bin_size;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::write_average(FILE *the_file){
    get_average();
    fprintf(the_file, "\ndata_total=%lf  , data_count=%lf  , data_average=%lf\n", data_total, data_count, data_average);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::write_histogram(FILE *the_file){
    unsigned i;
    std::string y_title_norm_area, y_title_norm_height;
    y_title_norm_area=y_title_norm_height=y_title;
    y_title_norm_area.append(" normalised (area)");
    y_title_norm_height.append(" normalised (height)");
    get_max();
    get_count();
    fprintf(the_file, "\n%s\n", title.c_str());
    if (has_names) fprintf(the_file, "%s  ,", bin_names_title.c_str());
    fprintf(the_file, "%s  ,  %s  ,  %s,  %s\n", x_title.c_str(), y_title.c_str(), y_title_norm_area.c_str(), y_title_norm_height.c_str());
    for (i=0; i<bin_length; i++){
        if (has_names) fprintf(the_file, "%s  ,", bin_names[i].c_str());
        fprintf(the_file, "%lf  , %lf  , %lf, %lf\n", resolveIndex(i), data[i], get_normalised_area_prob(i), get_normalised_height_prob(i));
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::write_Boltzmann(FILE *the_file) {
    unsigned i;
    std::string y_title_norm;
    y_title_norm=y_title;
    y_title_norm.append(" normalised");
    get_count();
    fprintf(the_file, "\n%s\n%s  ,  %s  ,  %s  ,  %s\n", title.c_str(), x_title.c_str(), y_title.c_str(), y_title_norm.c_str(), "Ni/N");
    for (i=0; i<bin_length; i++) fprintf(the_file, "%lf  , %lf  , %lf  , %lf\n", resolveIndex(i), data[i], data[i]/data_count, Z_component[i]/Z);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::set_value(const double in_x, const double in_y){
    data[getIndex(in_x)]=in_y;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::get_max(){
    unsigned i;
    data_max=0;
    for (i=0; i<bin_length; i++) if (data_max<data[i]) data_max=data[i];
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::get_count(){
    unsigned i;
    data_count=0;
    for (i=0; i<bin_length; i++) data_count+=data[i];
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::get_average(){ //incorporates get_count, so this can be called as an alternative to get_count if you also want the average
    unsigned i;
    data_average=data_total=data_count=0;
    for (i=0; i<bin_length; i++){
        data_count+=data[i];
        data_total+=(data[i]*resolveIndex(i));
    }
    if (data_count) data_average=data_total/data_count;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::set_Boltzmann_temp(const double in_T){
    unsigned i;
    double zi;
    T=in_T;
    Z=0;
    Z_component.clear(); //incase get_Z has been called before
    for (i=0;i<bin_length;i++){
        zi=0;//***** need to check if T==0 ... divide by 0!!
        if (T) zi=data[i]*exp(-1*resolveIndex(i)*JOULES_IN_WAVENUMBER/(BOLTZMANN_CONSTANT*T)); //zi=gi*e^(-E/KT) assuming the energy bins are in cm-1 and T is in Kelvin
        Z_component.push_back(zi); //probably not very efficient to clear and repopulate every time rather than just reseting the values ... might want to come back to this later if we end up recalculating the boltzmann distribution multiple times!
        Z+=zi;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::convolve(const double &width, const unsigned &method){
    switch(method){
        case GAUSSIAN_CONVOLVE: Gaussian_convolve(width); break;
        case LORENTZIAN_CONVOLVE: Lorentzian_convolve(width); break;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::Gaussian_convolve(const double &width){/*
    Gaussian: f(x) = a exp(-(x-b)^2/(2c^2))
    where a is height, b is centre position and c controls width
    the full width at half maximum = 2 sqrt(2 ln(2)) c ... so c = width / (2 sqrt(2 ln(2)))
    and the integral is 1 only if a = 1/(c sqrt(2 PI))
    so norm_factor = 2*sqrt(log(2))/(width*sqrt(PI)) normalises the gaussian to unit area, which then gets multiplied by the stick height
    exponent_factor = -4*log(2)/(width*width) comes from substituting the FWHM eqn for c into the Gaussian eqn */
    double exponent_factor, norm_factor, freq, df;
    unsigned i,j;
    convolved_data.clear();
    for (i=0;i<bin_length;i++) convolved_data.push_back(0);
    convolved=true;
    exponent_factor = -4*log(2)/pow(width,2);
    norm_factor = 2*sqrt(log(2))/(width*sqrt(PI));
    for (i=0;i<bin_length;i++){
        freq=resolveIndex(i);
        for (j=0;j<bin_length;j++){
             df=freq-resolveIndex(j);
             convolved_data[i]+=norm_factor*data[j]*exp(exponent_factor*pow(df,2));
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::Lorentzian_convolve(const double &width){ ///better check through this ... have not looked carefully at this yet because I have not actually used it yet
    double denominator_factor, norm_factor, freq, df;
    unsigned i,j;
    convolved_data.clear();
    for (i=0;i<bin_length;i++) convolved_data.push_back(0);
    convolved=true;
    denominator_factor = pow(width,2)/4.0;
    norm_factor = width/(2*PI);
    for (i=0;i<bin_length;i++){
        freq=resolveIndex(i);
        for (j=0;j<bin_length;j++){
             df=freq-resolveIndex(j);
             convolved_data[i]+=norm_factor*data[j]/(pow(df,2)+denominator_factor);
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram::get_normalised_area_prob(const unsigned i) {
    if (data_count) return data[i]/data_count;
    else return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram::get_normalised_height_prob(const unsigned i) {
    if (data_max) return data[i]/data_max;
    else return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram::get_Boltzmann_prob(const unsigned i) {
    if (Z) return Z_component[i]/Z;
    else return 0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::increment_value(const double in_x, const double in_y){
    if (on) data[getIndex(in_x)]+=in_y;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::range_increment(Histogram &in_hist, double scaling_factor){
    unsigned i;
    if (scaling_factor && on) { //possibly add (... && on) here ... but it's nice for debugging (I think histograms are not being initialised properly or something)
        if (abs(x_max-in_hist.x_max)>std::numeric_limits<double>::epsilon()) persistent_output << "Histogram::increment, Histograms misaligned ("<<x_max<<"!="<<in_hist.x_max<<") difference"<<(x_max-in_hist.x_max)<< std::endl;
        else { //this check needs to be more thorough ... should also check length and x_min
            for (i=0; i<bin_length; i++) {
                data[i]+=(in_hist.data[i]/scaling_factor);//The reason we divide here is that this way scaling factor will be very large and we get a result. Doing it the other way makes it VERY small, so it is morelikely to be too small to be distinguishable from 0 at double precision
            }
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram::inner_range_increment(Histogram &in_hist, double scaling_factor){ //this function assumes that the indexes of the data arrays in both histograms start with the same vaule and are separated by the same quantity!
    unsigned i, min_length;
    if (scaling_factor && on) {
        min_length=std::min(bin_length, in_hist.bin_length);
        for (i=0; i<min_length; i++) {
            data[i]+=(in_hist.data[i]/scaling_factor);
        }
    }
}

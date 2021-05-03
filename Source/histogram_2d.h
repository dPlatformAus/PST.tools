unsigned HISTOGRAM_2D_SIZE=1000;
//***********************************************************************************************************************************************************************************
class Histogram_2d
{
    public:
        bool on;
        unsigned bin_length1, bin_length2;
        boost::multi_array<double, 2> data;
        double x1_max, x2_max, x1_min, x2_min, bin_size1, bin_size2;
        std::string title, x1_title, x2_title, y_title;

        Histogram_2d();
        double& operator() (unsigned i, unsigned j) {return data[i][j];}
        Histogram_2d& operator=(const Histogram_2d &rhs);
        Histogram_2d& operator+=(const Histogram_2d &rhs);
        double& val(unsigned i, unsigned j) {return data[i][j];}
        void clear();
        void initialiseN(const std::string& in_title, const std::string& in_y_title,
                         const std::string& in_x1_title, const std::string& in_x2_title,
                         const unsigned &in_length1, const unsigned &in_length2,
                         const double &in_x1_max, const double &in_x2_max,
                         const double in_x1_min=0, const double in_x2_min=0);
        void initialiseS(const std::string& in_title, const std::string& in_y_title,
                         const std::string& in_x1_title, const std::string& in_x2_title,
                         const double &in_size1, const double &in_size2,
                         const double &in_x1_max, const double &in_x2_max,
                         const double in_x1_min=0, const double in_x2_min=0);
        unsigned getIndex1(const double &in_x);
        unsigned getIndex2(const double &in_x);
        double resolveIndex1(const unsigned &in_index);
        double resolveIndex2(const unsigned &in_index);
        void set_value(const double &in_x1, const double &in_x2, const double &in_y);
        void increment_value(const double &in_x1, const double &in_x2, const double &in_y);
        double get_y_max();
        void inner_range_increment(Histogram_2d &in_hist, double scaling_factor);
        void write_histogram(FILE *the_file);
        void write_bitmap(const std::string path_and_name_start);

    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & on & bin_length1 & bin_length2 & data & x1_max & x2_max & x1_min & x2_min & bin_size1 & bin_size2 & title & x1_title & x2_title & y_title;}
        void initialise(const std::string& in_title, const std::string& in_y_title,
                        const std::string& in_x1_title, const std::string& in_x2_title,
                        const double &in_x1_max, const double &in_x2_max,
                        const double &in_x1_min, const double &in_x2_min);
        double get_bin_size(const unsigned &in_length, const double &in_x_max, const double &in_x_min);
        unsigned get_bin_length(const double &in_size, const double &in_x_max, const double &in_x_min);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Histogram_2d::Histogram_2d(){
    bin_length1=bin_length2=0; //gotta initialise these otherwise the = operator will go crazy when it tries to resize the array!!
    x1_max=x2_max=x1_min=x2_min=bin_size1=bin_size2=0;
    on=false;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Histogram_2d& Histogram_2d::operator=(const Histogram_2d &rhs){
    on=rhs.on;
    title=rhs.title;
    x1_title=rhs.x1_title;
    x2_title=rhs.x2_title;
    y_title=rhs.y_title;
    x1_max=rhs.x1_max;
    x2_max=rhs.x2_max;
    x1_min=rhs.x1_min;
    x2_min=rhs.x2_min;
    bin_size1=rhs.bin_size1;
    bin_size2=rhs.bin_size2;
    bin_length1=rhs.bin_length1;
    bin_length2=rhs.bin_length2;
    data.resize(boost::extents[bin_length1][bin_length2]); //multi array is a little annoying, as the assignment operator insists that the dimensions must be the same ... cant initialise one with another ... you must resize it first
    data=rhs.data;
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Histogram_2d& Histogram_2d::operator+=(const Histogram_2d &rhs){ //this function assumes that both histograms are defined as the same (ie, titles etc), it only checks that their dimensions are the same
    unsigned i,j; //currently this is only used by the Dissociation_Result += operator, so simply incrementing the data and recalculating the total count and average is considered sufficient
    if (bin_length1==rhs.bin_length1 && bin_length2==rhs.bin_length2 && bin_size1==rhs.bin_size1 && bin_size2==rhs.bin_size2){
        for (i=0;i<bin_length1;i++) for (j=0;j<bin_length2;j++) data[i][j]+=rhs.data[i][j];
    }else{
        std::cout<<"Histogram_2d::operator+=("<<title<<", "<<rhs.title<<") HISTOGRAMS INCOMPATIBLE!"<<std::endl;
    }
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_2d::clear(){
    unsigned i,j;
    for (i=0;i<bin_length1;i++) for (j=0;j<bin_length2;j++) data[i][j]=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_2d::initialise(const std::string& in_title, const std::string& in_y_title,
                              const std::string& in_x1_title, const std::string& in_x2_title,
                              const double &in_x1_max, const double &in_x2_max,
                              const double &in_x1_min, const double &in_x2_min){
    unsigned i,j;
    on=true;
    title=in_title;
    y_title=in_y_title;
    x1_title=in_x1_title;
    x2_title=in_x2_title;
    x1_min=in_x1_min;
    x2_min=in_x2_min;
    x1_max=in_x1_max;
    x2_max=in_x2_max;
    data.resize(boost::extents[bin_length1][bin_length2]);
    for (i=0;i<bin_length1;i++) for (j=0;j<bin_length2;j++) data[i][j]=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_2d::initialiseN(const std::string& in_title, const std::string& in_y_title,
                               const std::string& in_x1_title, const std::string& in_x2_title,
                               const unsigned &in_length1, const unsigned &in_length2,
                               const double &in_x1_max, const double &in_x2_max,
                               const double in_x1_min, const double in_x2_min){
    bin_length1=in_length1;
    bin_length2=in_length2;
    bin_size1=get_bin_size(bin_length1, in_x1_max, in_x1_min);
    bin_size2=get_bin_size(bin_length2, in_x2_max, in_x2_min);
    initialise(in_title, in_y_title, in_x1_title, in_x2_title, in_x1_max, in_x2_max, in_x1_min, in_x2_min);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_2d::initialiseS(const std::string& in_title, const std::string& in_y_title,
                               const std::string& in_x1_title, const std::string& in_x2_title,
                               const double &in_size1, const double &in_size2,
                               const double &in_x1_max, const double &in_x2_max,
                               const double in_x1_min, const double in_x2_min){
    bin_size1=in_size1;
    bin_size2=in_size2;
    bin_length1=get_bin_length(bin_size1, in_x1_max, in_x1_min);
    bin_length2=get_bin_length(bin_size2, in_x2_max, in_x2_min);
    initialise(in_title, in_y_title, in_x1_title, in_x2_title, in_x1_max, in_x2_max, in_x1_min, in_x2_min);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_2d::get_bin_size(const unsigned &in_length, const double &in_x_max, const double &in_x_min){
    double result(1);
    if ((in_x_max-in_x_min)>0 && in_length>1) result=(in_x_max-in_x_min)/(in_length-1);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Histogram_2d::get_bin_length(const double &in_size, const double &in_x_max, const double &in_x_min){
    unsigned result(1);
    if ((in_x_max-in_x_min)>0) result=(unsigned)ceil(1+(in_x_max-in_x_min)/in_size);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Histogram_2d::getIndex1(const double &in_x){
    unsigned result(0);
    if (bin_size1) result=int(0.5 + (in_x-x1_min)/bin_size1);
    if (result>=bin_length1) {
        persistent_output<<"Histogram_2d::getIndex1("<<in_x<<")="<<result<<" "<<title<<"-"<<x1_title<<" index out of range!"<<std::endl;
        print_persistent_output();
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Histogram_2d::getIndex2(const double &in_x){
    unsigned result(0);
    if (bin_size2) result=int(0.5 + (in_x-x2_min)/bin_size2);
    if (result>=bin_length2) {
        persistent_output<<"Histogram_2d::getIndex2("<<in_x<<")="<<result<<" "<<title<<"-"<<x2_title<<" index out of range!"<<std::endl;
        print_persistent_output();
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_2d::resolveIndex1(const unsigned &in_index){
    return x1_min+in_index*bin_size1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_2d::resolveIndex2(const unsigned &in_index){
    return x2_min+in_index*bin_size2;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_2d::set_value(const double &in_x1, const double &in_x2, const double &in_y){
    data[getIndex1(in_x1)][getIndex2(in_x2)]=in_y;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_2d::increment_value(const double &in_x1, const double &in_x2, const double &in_y){
    if (on) data[getIndex1(in_x1)][getIndex2(in_x2)]+=in_y;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_2d::get_y_max(){
    unsigned i, j;
    double result(0);
    for(i=0;i<bin_length1;i++) for(j=0;j<bin_length2;j++) if(result<data[i][j]) result=data[i][j];
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_2d::inner_range_increment(Histogram_2d &in_hist, double scaling_factor){ //this function assumes that the indexes of the data arrays in both histograms start with the same vaule and are separated by the same quantity!
    unsigned i, j, min_length1, min_length2;
    if (scaling_factor && on) {
        min_length1=std::min(bin_length1, in_hist.bin_length1);
        min_length2=std::min(bin_length2, in_hist.bin_length2);
        for (i=0; i<min_length1; i++) {
            for (j=0; j<min_length2; j++) {
                data[i][j]+=(in_hist.data[i][j]/scaling_factor);
            }
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_2d::write_histogram(FILE *the_file){
    unsigned i,j;
    double row_total;
    std::vector<double>  column_totals;

    fprintf(the_file, "\n%s\n\n,%s,All %s\n%s,%s,total\nAll %s,total\n\n", title.c_str(), x2_title.c_str(), x2_title.c_str(), x1_title.c_str(), y_title.c_str(), x1_title.c_str());
    for (j=0; j<bin_length2; j++){
        fprintf(the_file, ", %lf", resolveIndex2(j));
        column_totals.push_back(0);
    }
    fprintf(the_file, ", Row Total\n");
    for (i=0; i<bin_length1; i++){
        fprintf(the_file, "%lf", resolveIndex1(i));
        row_total=0;
        for (j=0; j<bin_length2; j++){
            fprintf(the_file, ", %lf", data[i][j]);
            row_total+=data[i][j];
            column_totals[j]+=data[i][j];
        }
        fprintf(the_file, ", %lf\n", row_total);
    }
    fprintf(the_file, "Column Total");
    for (j=0; j<bin_length2; j++) fprintf(the_file, ", %lf", column_totals[j]);
    fprintf(the_file, "\n");
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_2d::write_bitmap(const std::string path_and_name_start){
    std::stringstream ss;
    std::string path_and_name;
    ss<<path_and_name_start<<title<<".bmp";
    path_and_name=ss.str();
    BMP AnImage;
    Pixel_Colour pixel;
    double y_max;
    unsigned i,j;
    y_max=get_y_max();
    /*AnImage.SetSize(bin_length1, bin_length2);
    AnImage.SetBitDepth(32);
    for (i=0;i<bin_length1;i++){
        for (j=0; j<bin_length2;j++){
            pixel=get_colour(data[i][j], y_max);
            set_pixel(AnImage(i,bin_length2-j-1), pixel);
        }
    } that was the old way with x1 horizontal and x2 vertical */
    AnImage.SetSize(bin_length2, bin_length1); //(width, height) so x1 is vertical and x2 is horizontal
    AnImage.SetBitDepth(32);
    for (i=0;i<bin_length1;i++){
        for (j=0; j<bin_length2;j++){
            pixel=get_colour(data[i][j], y_max);
            set_pixel(AnImage(j, bin_length1-i-1), pixel);
        }
    }
    AnImage.WriteToFile(path_and_name.c_str());
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

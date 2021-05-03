//***********************************************************************************************************************************************************************************
class Histogram_3d
{
    public:
        bool on;
        unsigned bin_length1, bin_length2, bin_length3;
        boost::multi_array<double, 3> data;
        double x1_max, x2_max, x3_max, x1_min, x2_min, x3_min, bin_size1, bin_size2, bin_size3;
        std::string title, x1_title, x2_title, x3_title, y_title;

        Histogram_3d();
        double& operator() (unsigned i, unsigned j, unsigned k) {return data[i][j][k];}
        Histogram_3d& operator=(const Histogram_3d &rhs);
        Histogram_3d& operator+=(const Histogram_3d &rhs);
        double& val(unsigned i, unsigned j, unsigned k) {return data[i][j][k];}
        void initialiseN(const std::string& in_title, const std::string& in_y_title,
                         const std::string& in_x1_title, const std::string& in_x2_title, const std::string& in_x3_title,
                         const unsigned &in_length1, const unsigned &in_length2, const unsigned &in_length3,
                         const double &in_x1_max, const double &in_x2_max, const double &in_x3_max,
                         const double in_x1_min=0, const double in_x2_min=0, const double in_x3_min=0);
        void initialiseS(const std::string& in_title, const std::string& in_y_title,
                         const std::string& in_x1_title, const std::string& in_x2_title, const std::string& in_x3_title,
                         const double &in_size1, const double &in_size2, const double &in_size3,
                         const double &in_x1_max, const double &in_x2_max, const double &in_x3_max,
                         const double in_x1_min=0, const double in_x2_min=0, const double in_x3_min=0);
        unsigned getIndex1(const double &in_x);
        unsigned getIndex2(const double &in_x);
        unsigned getIndex3(const double &in_x);
        double resolveIndex1(const unsigned &in_index);
        double resolveIndex2(const unsigned &in_index);
        double resolveIndex3(const unsigned &in_index);
        void set_value(const double &in_x1, const double &in_x2, const double &in_x3, const double &in_y);
        void increment_value(const double &in_x1, const double &in_x2, const double &in_x3, const double &in_y);
        unsigned getNonZeroCount();
        double getNonZeroIndex1Sum();
        double getNonZeroIndex2Sum();
        double getNonZeroIndex3Sum();

    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & on & bin_length1 & bin_length2 & bin_length3 & data & x1_max & x2_max & x3_max & x1_min & x2_min & x3_min & bin_size1 & bin_size2 & bin_size3 & title & x1_title & x2_title & x3_title & y_title;}
        void initialise(const std::string& in_title, const std::string& in_y_title,
                        const std::string& in_x1_title, const std::string& in_x2_title, const std::string& in_x3_title,
                        const double &in_x1_max, const double &in_x2_max, const double &in_x3_max,
                        const double &in_x1_min, const double &in_x2_min, const double &in_x3_min);
        double get_bin_size(const unsigned &in_length, const double &in_x_max, const double &in_x_min);
        unsigned get_bin_length(const double &in_size, const double &in_x_max, const double &in_x_min);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Histogram_3d::Histogram_3d(){
    bin_length1=bin_length2=bin_length3=0;
    x1_max=x2_max=x3_max=x1_min=x2_min=x3_min=bin_size1=bin_size2=bin_size3=0;
    on=false;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Histogram_3d& Histogram_3d::operator=(const Histogram_3d &rhs){
    on=rhs.on;
    title=rhs.title;
    x1_title=rhs.x1_title;
    x2_title=rhs.x2_title;
    x3_title=rhs.x3_title;
    y_title=rhs.y_title;
    x1_max=rhs.x1_max;
    x2_max=rhs.x2_max;
    x3_max=rhs.x3_max;
    x1_min=rhs.x1_min;
    x2_min=rhs.x2_min;
    x3_min=rhs.x3_min;
    bin_size1=rhs.bin_size1;
    bin_size2=rhs.bin_size2;
    bin_size3=rhs.bin_size3;
    bin_length1=rhs.bin_length1;
    bin_length2=rhs.bin_length2;
    bin_length3=rhs.bin_length3;
    data.resize(boost::extents[bin_length1][bin_length2][bin_length3]); //multi array is a little annoying, as the assignment operator insists that the dimensions must be the same ... cant initialise one with another ... you must resize it first
    data=rhs.data;
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Histogram_3d& Histogram_3d::operator+=(const Histogram_3d &rhs){ //this function assumes that both histograms are defined as the same (ie, titles etc), it only checks that their dimensions are the same
    unsigned i,j,k; //currently this is only used by the Dissociation_Result += operator, so simply incrementing the data and recalculating the total count and average is considered sufficient
    if (bin_length1==rhs.bin_length1 && bin_length2==rhs.bin_length2 && bin_length3==rhs.bin_length3 && bin_size1==rhs.bin_size1 && bin_size2==rhs.bin_size2 && bin_size3==rhs.bin_size3){
        for (i=0;i<bin_length1;i++) for (j=0;j<bin_length2;j++) for (k=0;k<bin_length3;k++) data[i][j][k]+=rhs.data[i][j][k];
    }else{
        std::cout<<"Histogram_3d::operator+=("<<title<<", "<<rhs.title<<") HISTOGRAMS INCOMPATIBLE!"<<std::endl;
    }
    return *this;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_3d::initialise(const std::string& in_title, const std::string& in_y_title,
                              const std::string& in_x1_title, const std::string& in_x2_title, const std::string& in_x3_title,
                              const double &in_x1_max, const double &in_x2_max, const double &in_x3_max,
                              const double &in_x1_min, const double &in_x2_min, const double &in_x3_min){
    unsigned i,j,k;
    on=true;
    title=in_title;
    y_title=in_y_title;
    x1_title=in_x1_title;
    x2_title=in_x2_title;
    x3_title=in_x3_title;
    x1_min=in_x1_min;
    x2_min=in_x2_min;
    x3_min=in_x3_min;
    x1_max=in_x1_max;
    x2_max=in_x2_max;
    x3_max=in_x3_max;
    //persistent_output<<in_title<<" resize "<<bin_length1<<", "<<bin_length2<<", "<<bin_length3<<endl;
    //print_semipersistent_output();
    data.resize(boost::extents[bin_length1][bin_length2][bin_length3]);
    for (i=0;i<bin_length1;i++) for (j=0;j<bin_length2;j++) for (k=0;k<bin_length3;k++) data[i][j][k]=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_3d::initialiseN(const std::string& in_title, const std::string& in_y_title,
                               const std::string& in_x1_title, const std::string& in_x2_title, const std::string& in_x3_title,
                               const unsigned &in_length1, const unsigned &in_length2, const unsigned &in_length3,
                               const double &in_x1_max, const double &in_x2_max, const double &in_x3_max,
                               const double in_x1_min, const double in_x2_min, const double in_x3_min){ //need to add initialiseS too later, see histogram class
    bin_length1=in_length1;
    bin_length2=in_length2;
    bin_length3=in_length3;
    bin_size1=get_bin_size(bin_length1, in_x1_max, in_x1_min);
    bin_size2=get_bin_size(bin_length2, in_x2_max, in_x2_min);
    bin_size3=get_bin_size(bin_length3, in_x3_max, in_x3_min);
    initialise(in_title, in_y_title, in_x1_title, in_x2_title, in_x3_title, in_x1_max, in_x2_max, in_x3_max, in_x1_min, in_x2_min, in_x3_min);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_3d::initialiseS(const std::string& in_title, const std::string& in_y_title,
                               const std::string& in_x1_title, const std::string& in_x2_title, const std::string& in_x3_title,
                               const double &in_size1, const double &in_size2, const double &in_size3,
                               const double &in_x1_max, const double &in_x2_max, const double &in_x3_max,
                               const double in_x1_min, const double in_x2_min, const double in_x3_min){
    bin_size1=in_size1;
    bin_size2=in_size2;
    bin_size3=in_size3;
    bin_length1=get_bin_length(bin_size1, in_x1_max, in_x1_min);
    bin_length2=get_bin_length(bin_size2, in_x2_max, in_x2_min);
    bin_length3=get_bin_length(bin_size3, in_x3_max, in_x3_min);
    initialise(in_title, in_y_title, in_x1_title, in_x2_title, in_x3_title, in_x1_max, in_x2_max, in_x3_max, in_x1_min, in_x2_min, in_x3_min);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_3d::get_bin_size(const unsigned &in_length, const double &in_x_max, const double &in_x_min){
    double result(1);
    if ((in_x_max-in_x_min)>0 && in_length>1) result=(in_x_max-in_x_min)/(in_length-1);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Histogram_3d::get_bin_length(const double &in_size, const double &in_x_max, const double &in_x_min){
    unsigned result(1);
    if ((in_x_max-in_x_min)>0) result=(unsigned)ceil(1+(in_x_max-in_x_min)/in_size);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Histogram_3d::getIndex1(const double &in_x){
    unsigned result(0);
    if (bin_size1) result=int(0.5 + (in_x-x1_min)/bin_size1);
    if (result>=bin_length1) {
        persistent_output<<"Histogram_3d::getIndex1("<<in_x<<")="<<result<<" "<<title<<"-"<<x1_title<<" index out of range!"<<std::endl;
        print_persistent_output();
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Histogram_3d::getIndex2(const double &in_x){
    unsigned result(0);
    if (bin_size2) result=int(0.5 + (in_x-x2_min)/bin_size2);
    if (result>=bin_length2) {
        persistent_output<<"Histogram_3d::getIndex2("<<in_x<<")="<<result<<" "<<title<<"-"<<x2_title<<" index out of range!"<<std::endl;
        print_persistent_output();
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Histogram_3d::getIndex3(const double &in_x){
    unsigned result(0);
    if (bin_size3) result=int(0.5 + (in_x-x3_min)/bin_size3);
    if (result>=bin_length3) {
        persistent_output<<"Histogram_3d::getIndex3("<<in_x<<")="<<result<<" >= "<<bin_length3<<title<<" -- "<<x3_title<<" index out of range!"<<std::endl;
        print_persistent_output();
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_3d::resolveIndex1(const unsigned &in_index){
    return x1_min+in_index*bin_size1;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_3d::resolveIndex2(const unsigned &in_index){
    return x2_min+in_index*bin_size2;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_3d::resolveIndex3(const unsigned &in_index){
    return x3_min+in_index*bin_size3;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_3d::set_value(const double &in_x1, const double &in_x2, const double &in_x3, const double &in_y){
    data[getIndex1(in_x1)][getIndex2(in_x2)][getIndex3(in_x3)]=in_y;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Histogram_3d::increment_value(const double &in_x1, const double &in_x2, const double &in_x3, const double &in_y){
    if (on) data[getIndex1(in_x1)][getIndex2(in_x2)][getIndex3(in_x3)]+=in_y;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned Histogram_3d::getNonZeroCount(){ //returns the number of bins with non zero count
    unsigned result(0), i, j, k;
    for (i=0;i<bin_length1;i++) for (j=0;j<bin_length2;j++) for (k=0;k<bin_length3;k++) if(data[i][j][k]) result++;
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_3d::getNonZeroIndex1Sum(){ //returns the sum of resolved index values for index 1 for all non zero count bins
    unsigned i, j, k;
    double result(0);
    for (i=0;i<bin_length1;i++) for (j=0;j<bin_length2;j++) for (k=0;k<bin_length3;k++) if(data[i][j][k]) result+=resolveIndex1(i);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_3d::getNonZeroIndex2Sum(){ //returns the sum of resolved index values for index 2 for all non zero count bins
    unsigned i, j, k;
    double result(0);
    for (i=0;i<bin_length1;i++) for (j=0;j<bin_length2;j++) for (k=0;k<bin_length3;k++) if(data[i][j][k]) result+=resolveIndex2(j);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Histogram_3d::getNonZeroIndex3Sum(){ //returns the sum of resolved index values for index 3 for all non zero count bins
    unsigned i, j, k;
    double result(0);
    for (i=0;i<bin_length1;i++) for (j=0;j<bin_length2;j++) for (k=0;k<bin_length3;k++) if(data[i][j][k]) result+=resolveIndex3(k);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
#if (defined(USING_MPI))
//this is the binary function object used for the MPI gather routine
struct Reduce_Histogram_3d : public std::binary_function<Histogram_3d*, Histogram_3d*, Histogram_3d*>
{
    Histogram_3d* operator() (Histogram_3d *a, Histogram_3d *b){
        (*b)+=(*a);
        return b;
    }
};
namespace boost { namespace mpi {//this allows boost to use a more efficient gather algorithm, by telling boost that Reduce_Histogram_3d commutes
    template<>
    struct is_commutative<Reduce_Histogram_3d, Histogram_3d*> : mpl::true_{};
}}
#endif
//***********************************************************************************************************************************************************************************

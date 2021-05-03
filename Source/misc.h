#define HOURS_IN_DAY 24
#define MINUTES_IN_HOUR 60
#define SECONDS_IN_MINUTE 60

#define PI 3.14159265358979
#define PLANCKS_CONSTANT 6.62606876E-34
#define SPEED_OF_LIGHT 2.99792458E8
#define AVOGADROS_NUMBER 6.02214199E23
#define METERS_IN_ANGSTROM 1E-10
#define WAVENUMBER_IN_KJ_PER_MOL 83.59346093
#define KJ_IN_KCAL 4.184
#define BOLTZMANN_CONSTANT 1.3806488E-23
#define WAVENUMBER_IN_KELVIN 0.695
#define NEWTONS_PER_METER_IN_MDYNE_PER_ANGSTROM 100

const double HBAR=PLANCKS_CONSTANT/(2*PI);
const double HBARSD=pow(HBAR,2);
const double KG_IN_AMU=1/(AVOGADROS_NUMBER*1000);
const double JOULES_IN_WAVENUMBER=PLANCKS_CONSTANT*SPEED_OF_LIGHT*100;
const double JOULES_PER_MOLE_IN_WAVENUMBER=JOULES_IN_WAVENUMBER*AVOGADROS_NUMBER;
const double KJOULES_PER_MOLE_IN_WAVENUMBER=JOULES_PER_MOLE_IN_WAVENUMBER/1000;
const double KG_SQM_IN_AMU_SQA=KG_IN_AMU*pow(METERS_IN_ANGSTROM,2);
const double I_AMU_SQANGSTROM_TO_B_WAVENUMBER = pow(PLANCKS_CONSTANT,2)/(8*pow(PI,2)*KG_IN_AMU*pow(METERS_IN_ANGSTROM,2)*JOULES_IN_WAVENUMBER);
const double I_AMU_SQANGSTROM_TO_B_HZ = PLANCKS_CONSTANT/(8*pow(PI,2)*KG_IN_AMU*pow(METERS_IN_ANGSTROM,2));
const double WAVENUMBER_ANGSTROM6_TO_JOULE_METER6 = JOULES_IN_WAVENUMBER/pow(1E10,6);
const double HZ_IN_WAVENUMBER=100*SPEED_OF_LIGHT;
const double RAD_IN_DEGREE=PI/180;
const double RAD_PER_SEC_IN_WAVENUMBER=HZ_IN_WAVENUMBER*2*PI;

std::string pretty_time(const double in_seconds){
    std::string result;
    int minutes(0), hours(0), days(0);
    float f_seconds(0), f_minutes(0), f_hours(0);
    f_seconds=fmod(in_seconds,SECONDS_IN_MINUTE);
    f_minutes=in_seconds/SECONDS_IN_MINUTE;
    minutes=fmod(f_minutes,MINUTES_IN_HOUR);
    f_hours=f_minutes/MINUTES_IN_HOUR;
    hours=fmod(f_hours,HOURS_IN_DAY);
    days=f_hours/HOURS_IN_DAY;
    if (days) result += (boost::lexical_cast<std::string>(days)+" days ");
    if (hours) result += (boost::lexical_cast<std::string>(hours)+" hours ");
    if (minutes) result += (boost::lexical_cast<std::string>(minutes)+" minutes ");
    if (f_seconds) result += (boost::lexical_cast<std::string>(f_seconds)+" seconds ");
    if (result=="") result = "much less than a second";
    return result;
}
//***********************************************************************************************************************************************************************************
class System_Info {
    public:
        std::string platform;
        std::string boost_platform;
        char str_time_stamp[16];
        int start_time, finish_time;

        System_Info();
        void clear_screen();
        void make_time_stamp_string();
        void check_execution_time();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
System_Info::System_Info() {
    start_time = time(NULL);
    boost_platform=BOOST_PLATFORM;
    if (boost_platform == "Win32" || boost_platform == "Win64" || boost_platform == "Cygwin" ) {
        platform="Windows";
    }else{
        platform="POSIX";
    }
    make_time_stamp_string();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void System_Info::clear_screen(){
	int temp(0);
    if ( platform == "Windows" ){
        temp = std::system ( "CLS" );
    }else if ( platform  == "POSIX" ){
        temp = std::system ( "clear" );
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void System_Info::check_execution_time(){
	finish_time = time(NULL);
	fprintf(stderr, "Execution time: %d\n", finish_time-start_time);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void System_Info::make_time_stamp_string(){
	time_t current_time;
	time(&current_time);
	struct tm *ptm = localtime(&current_time);
	int the_year, the_month, the_day, the_hour, the_minute, the_second;
	char strmonth[3], strday[3], strhr[3], strmin[3], strsec[3];
	the_year = 1900 + ptm->tm_year;
	the_month = 1 + ptm->tm_mon;

	if (the_month>9){
	  sprintf(strmonth, "%d", the_month);
	}else{
	  sprintf(strmonth, "%s%d", "0", the_month);
	}

	the_day = ptm->tm_mday;
	if (the_day>9){
	  sprintf(strday, "%d", the_day);
	}else{
	  sprintf(strday, "%s%d", "0", the_day);
	}

	the_hour = ptm->tm_hour;
	if (the_hour>9){
	  sprintf(strhr, "%d", the_hour);
	}else{
	  sprintf(strhr, "%s%d", "0", the_hour);
	}

	the_minute = ptm->tm_min;
	if (the_minute>9){
	  sprintf(strmin, "%d", the_minute);
	}else{
	  sprintf(strmin, "%s%d", "0", the_minute);
	}

	the_second = ptm->tm_sec;
	if (the_second>9){
	  sprintf(strsec, "%d", the_second);
	}else{
	  sprintf(strsec, "%s%d", "0", the_second);
	}

	sprintf(str_time_stamp, "%d%s%s%s%s%s_", the_year, strmonth, strday, strhr, strmin, strsec);
}
//***********************************************************************************************************************************************************************************
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
System_Info sys_info;
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int make_directory(const char* dir_name) {
#if (defined(_WIN32) || defined(_WIN64) || defined(WIN32) || defined(WIN64)) && !defined(UNIX)
    return mkdir(dir_name);
#else
    return mkdir(dir_name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector <std::string> read_directory( const std::string& path = std::string() )
{
  std::vector <std::string> result;
  dirent* de;
  DIR* dp;
  //err_no= 0;
  dp = opendir( path.empty() ? "." : path.c_str() );
  if (dp)
    {
    while (true)
      {
      //err_no= 0;
      de = readdir( dp );
      if (de == NULL) break;
      if ( !( (de->d_name[0] == '.' && NAMLEN(de) == 1) || (de->d_name[1] == '.' && NAMLEN(de) == 2) ) ) result.push_back( std::string( de->d_name ) );
      }
    closedir( dp );
    std::sort( result.begin(), result.end() );
    }
  return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool directory_exists(std::string dir_path) {
    bool result=true;
    DIR *dp;
    dp = opendir(dir_path.c_str());
    if (dp == NULL) result=false;
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::string get_file_extension(const std::string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void clear_semipersistent_output(){
    semipersistent_output.str("");
    semipersistent_output.clear();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void print_semipersistent_output(){
    if (DO_CLEAR_SCREEN) sys_info.clear_screen();
    std::cout << persistent_output.str() << semipersistent_output.str();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void print_persistent_output(){
    clear_semipersistent_output();
    print_semipersistent_output();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double old_factorial(int n) { //this is the slow way! use the datastructure log_factorial defined in factorials.h
    double answer=1.0;
    while (n>1) {
        answer*=(double)n;
        n--;
    }
    return answer;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double get_reduced_mass(const double &a, const double &b){
    return a*b/(a+b);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<double> cross_product(const std::vector<double> &a, const std::vector<double> &b, const unsigned L){
    std::vector<double> result;
    unsigned i;
    for (i=0; i<L; i++) result.push_back( a[(i+1)%L] * b[(i+2)%L] - ( a[(i+L-1)%L] * b[(i+L-2)%L] ) );
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double vector_magnitude(const std::vector<double> &a){
    double sq_sum=0;
    unsigned i;
    for (i=0; i<a.size(); i++) sq_sum+=pow(a[i],2);
    return sqrt(sq_sum);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<double> vector_normalise(const std::vector<double> &a){
    unsigned i;
    std::vector<double> result;
    double a_mag=vector_magnitude(a);
    for (i=0; i<a.size(); i++) result.push_back(a[i]/a_mag);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double vector_angle(const std::vector<double> &a, const std::vector<double> &b){ //assumes both vectors are 3 dimensional (or at least equidimensional)
    unsigned i;
    double a_mag, b_mag, dot_product(0);
    a_mag=vector_magnitude(a);
    b_mag=vector_magnitude(b);
    for (i=0; i<a.size(); i++) dot_product+=(a[i]/a_mag)*(b[i]/b_mag);
    return acos(dot_product);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void normaliseFractions(double &f_1, double &f_2) {
    double factor;
    if (f_1+f_2) {
        factor=1/(f_1+f_2);
        f_1*=factor;
        f_2*=factor;
    }else{
        f_1=f_2=0; // I suppose mathematically this is the most correct answer, HOWEVER if there is nothing there is nothing so I think PHYSICALLY 0 is the correct answer :-/
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//***********************************************************************************************************************************************************************************
class Data_File
{
	public:
		FILE *fp;
		std::string file_name;
		int file_open;
		bool use_timestamps;

		Data_File();
		~Data_File();
		void close();
		void warning_close();
};
class Input_File : public Data_File
{
    public:
        std::string input_directory;

		void set_input_directory(const std::string in_dir);
		void open(const std::string suffix);

};
class Output_File : public Data_File
{
    public:
        std::string output_directory;

		void set_output_directory(const std::string out_dir);
		void open(const std::string suffix);

};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Data_File::Data_File() {
	file_open=0;
	use_timestamps=0; //if true timestamp added to all filenames for openeing and closing (to use, may need to differentiate between opening with closing with timestamp added)
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Data_File::~Data_File() {
	warning_close();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Data_File::close() {
	if (file_open){
		fclose(fp);
		file_open=0;
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Data_File::warning_close() {
	if (file_open) {
		close();
		fprintf(stderr, "\nwarning: file \"%s\" closed implicitly\n", file_name.c_str());
	}
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Input_File::set_input_directory(const std::string in_dir) {
	input_directory=in_dir;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Input_File::open(const std::string suffix) {
	warning_close();
	if (use_timestamps) {
        file_name=input_directory+sys_info.str_time_stamp+suffix;
	}else{
        file_name=input_directory+suffix;
	}
	if ( (fp=fopen(file_name.c_str(), "r")) != NULL) {
		file_open=1;
	} else {
		fprintf(stderr, "Error: failed to load file: %s\n", file_name.c_str());
		exit(1);
   }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Output_File::set_output_directory(const std::string out_dir) {
    int status(0);
    std::stringstream ss;
    std::string path_segment, folder;
    path_segment="";
	output_directory=out_dir;
    ss << output_directory;
    while(getline(ss,folder,'/')) {
        path_segment += (folder+"/");
        if (!directory_exists(path_segment)) status=make_directory(path_segment.c_str());
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Output_File::open(const std::string suffix) {
	warning_close();
	if (use_timestamps) {
        file_name=output_directory+sys_info.str_time_stamp+suffix;
	}else{
        file_name=output_directory+suffix;
	}
	fp = fopen(file_name.c_str(), "w");
	file_open=1;
}
//***********************************************************************************************************************************************************************************



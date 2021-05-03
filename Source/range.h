template <class T>
class Number_Range {
    public:
        T min, max;
        void setRange(const T in_min, const T in_max);
        void setRange(const T in_val){setRange(in_val,in_val);};
        void setRange(const std::string values);
        void clear();
        Number_Range();
    protected:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version){ar & min & max;}
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <class T>
Number_Range<T>::Number_Range() {
    clear();
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <class T>
void Number_Range<T>::clear(){
    min=max=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <class T>
void Number_Range<T>::setRange(const T in_min, const T in_max){
    min=in_min;
    max=in_max;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template <class T>
void Number_Range<T>::setRange(const std::string values){
    size_t colon_position, comma_position, delimiter_position;
    std::string min_string, max_string;
    colon_position = values.find(":");
    comma_position = values.find(",");
    if (colon_position != std::string::npos) delimiter_position=colon_position; //I originally coded this to use colon as delimiter, because the syntax for probes was to delimit different ranges with comma and the min and max within that with colon
    else if (comma_position != std::string::npos) delimiter_position=comma_position; //however, other ranges are probably/maybe neater in the input file if they are delimited by comma ... so if no colon found, look for a comma.
    else delimiter_position=std::string::npos;
    if (delimiter_position != std::string::npos){ // then two values entered
        min_string=values.substr(0, delimiter_position);
        max_string=values.substr(delimiter_position + 1);
        boost::trim(min_string); //my first version of this was not templated and I had separate double, float and int versions that used either atoi of atof appropriately
        boost::trim(max_string); //atoi and atof don't get confused by white space (actually they appear to find numeric content inside all sorts of strings), however, the lexical cast baulks at white space, so trim it!
        try {
            min=boost::lexical_cast<T>(min_string.c_str());
        } catch (const boost::bad_lexical_cast &e) {
            persistent_output<<"Number_Range<T> Caught bad lexical cast with error "<<e.what()<<" for "<<min_string.c_str()<<std::endl;
        }
        try {
            max=boost::lexical_cast<T>(max_string.c_str());
        } catch (const boost::bad_lexical_cast &e) {
            persistent_output<<"Number_Range<T> Caught bad lexical cast with error "<<e.what()<<" for "<<max_string.c_str()<<std::endl;
        }
    }else{//single value specified, so set both min and max to that value and a loop from min to max is going to produce the simgle value result expected
        min_string=values;
        boost::trim(min_string);
        try {
            min=max=boost::lexical_cast<T>(min_string.c_str());
        } catch (const boost::bad_lexical_cast &e) {
            persistent_output<<"Number_Range<T> Caught bad lexical cast with error "<<e.what()<<" for "<<min_string.c_str()<<std::endl;
        }
    }
}

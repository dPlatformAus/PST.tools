class Fit_Point
{
    public:
        double energy;
        std::vector<double> radical, roaming;
        bool is_radical, is_roaming;

        Fit_Point();
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Fit_Point::Fit_Point() {
    is_radical=is_roaming=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class Fit_Point_List
{
    public:
        bool is_quantum;
        std::vector<Fit_Point> fit_energies;

        Fit_Point_List();
        const Fit_Point& operator[] (const unsigned i) const {return fit_energies[i];}
        unsigned size() const {return fit_energies.size();}
        bool energy_index(const double target_energy, unsigned &index);
        void add_radical(const double target_energy, const double target_value);
        void add_roaming(const double target_energy, const double target_value);
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Fit_Point_List::Fit_Point_List(){
    is_quantum=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bool Fit_Point_List::energy_index(const double target_energy, unsigned &index) {
    unsigned i;
    bool result=0;
    for (i=0; i<fit_energies.size(); i++) {
        if (fit_energies[i].energy == target_energy) {
            index=i;
            result=true;
        }
    }
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Fit_Point_List::add_radical(const double target_energy, const double target_value) {
    unsigned index;
    Fit_Point temp_fit_point;
    if (energy_index(target_energy, index)) {
        fit_energies[index].radical.push_back(target_value);
        fit_energies[index].is_radical=1;
    } else {
        temp_fit_point.is_radical=1;
        temp_fit_point.energy=target_energy;
        temp_fit_point.radical.push_back(target_value);
        fit_energies.push_back(temp_fit_point);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Fit_Point_List::add_roaming(const double target_energy, const double target_value) {
    unsigned index;
    Fit_Point temp_fit_point;
    if (energy_index(target_energy, index)) {
        fit_energies[index].roaming.push_back(target_value);
        fit_energies[index].is_roaming=1;
    } else {
        temp_fit_point.is_roaming=1;
        temp_fit_point.energy=target_energy;
        temp_fit_point.roaming.push_back(target_value);
        fit_energies.push_back(temp_fit_point);
    }
}
